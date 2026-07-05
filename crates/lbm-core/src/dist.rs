//! MPI distribution (feature `mpi`, docs/MPI_GUIDE.md): the [`MpiExchange`]
//! halo implementation and the [`MpiSolver`] driver for one-part-per-rank
//! runs (docs/ARCHITECTURE_V2.md §2.3, docs/HPC_SCALING.md step 3).
//!
//! ## Exchange plan
//!
//! [`MpiExchange`] runs the *identical* x → y → z phase plan and pack/unpack
//! layer maths as the in-process exchanges (`crate::halo`): per axis, both
//! side layers are posted as `MPI_Isend`/`MPI_Irecv` pairs and completed
//! before the next axis packs, so corner/edge halo data is forwarded through
//! face neighbours exactly like `InProcess` hands buffers across. A message
//! for the receiver's face `F` carries tag `base + F.index()`; a periodic
//! axis with one part wraps locally without touching MPI.
//!
//! Buffers travel as raw bytes of the scalar type (`f32`/`f64` are POD), so
//! distributed fields are **bit-identical** to a single-rank run; only the
//! `f64` diagnostics reassociate (rank partial sums → `MPI_Allreduce`).
//!
//! ## Rank topology
//!
//! Rank `r` owns part `r` of [`crate::solver::partition`]'s Cartesian
//! decomposition (part id = `(pz·dy + py)·dx + px`); neighbour ids stored in
//! the local [`Subdomain`] are therefore already rank numbers, and no
//! `MPI_Cart` communicator is needed.
//!
//! ## Collectives contract
//!
//! Everything that exchanges halos or reduces diagnostics is collective over
//! the communicator: `step`, `init_with`, `update_shan_chen_force`,
//! `total_mass` / `total_momentum` / `nonfinite_count`, the `gather_*`
//! family, and mask edits (`set_solid` must be *called on every rank* with
//! the same global coordinates — remote owners just mark their halo dirty).

use mpi::collective::SystemOperation;
use mpi::topology::SimpleCommunicator;
use mpi::traits::{Communicator, CommunicatorCollectives, Destination, Root, Source};
use mpi::{Rank, Tag};
use std::cell::RefCell;

use crate::backend::Backend;
use crate::fields::SoaFields;
use crate::halo::{
    layer_cell_count, layer_indices, pack_f_layer, pack_scalar_layer, unpack_f_layer,
    unpack_scalar_layer, ExchangeScope, HaloExchange,
};
use crate::lattice::{Face, Lattice};
use crate::params::{CollisionKind, FaceBC};
use crate::real::Real;
use crate::solver::{partition, GlobalSpec, Solver};
use crate::subdomain::Subdomain;

const TAG_F: Tag = 100;
const TAG_SCALAR: Tag = 200;
const TAG_MASK_META: Tag = 300;
const TAG_MASK_WALL: Tag = 400;
const TAG_GATHER: Tag = 500;
const TAG_MASS_ROWS: Tag = 600;

#[derive(Clone)]
struct AxisBuffers<E> {
    send: [Vec<E>; 2],
    recv: [Vec<E>; 2],
}

impl<E> Default for AxisBuffers<E> {
    fn default() -> Self {
        Self {
            send: [Vec::new(), Vec::new()],
            recv: [Vec::new(), Vec::new()],
        }
    }
}

struct MpiBuffers<T: Real> {
    f: [AxisBuffers<T>; 3],
    scalar: [AxisBuffers<T>; 3],
    mask_meta: [AxisBuffers<u8>; 3],
    mask_wall: [AxisBuffers<T>; 3],
}

impl<T: Real> Default for MpiBuffers<T> {
    fn default() -> Self {
        Self {
            f: std::array::from_fn(|_| AxisBuffers::default()),
            scalar: std::array::from_fn(|_| AxisBuffers::default()),
            mask_meta: std::array::from_fn(|_| AxisBuffers::default()),
            mask_wall: std::array::from_fn(|_| AxisBuffers::default()),
        }
    }
}

/// Raw-byte view of a POD slice (`f32` / `f64` / `u8`: no padding, every bit
/// pattern valid). Private to this module; only instantiated with those.
fn as_bytes<E: Copy>(v: &[E]) -> &[u8] {
    // SAFETY: E is POD (module invariant), len·size fits (it came from a slice).
    unsafe { std::slice::from_raw_parts(v.as_ptr().cast::<u8>(), std::mem::size_of_val(v)) }
}

fn as_bytes_mut<E: Copy>(v: &mut [E]) -> &mut [u8] {
    // SAFETY: as above; every byte pattern is a valid E for the types used.
    unsafe { std::slice::from_raw_parts_mut(v.as_mut_ptr().cast::<u8>(), std::mem::size_of_val(v)) }
}

/// Choose a Cartesian rank decomposition whose part count is `ranks` and
/// whose approximate block surface area is minimal for `dims`.
pub fn choose_decomp(d: usize, dims: [usize; 3], ranks: usize) -> [usize; 3] {
    assert!(ranks > 0, "rank count must be positive");
    assert!(d == 2 || d == 3, "dimension must be 2 or 3");
    let mut best = [ranks, 1, 1];
    let mut best_score = f64::INFINITY;
    for dx in 1..=ranks {
        if ranks % dx != 0 {
            continue;
        }
        let rem = ranks / dx;
        for dy in 1..=rem {
            if rem % dy != 0 {
                continue;
            }
            let dz = rem / dy;
            if d == 2 && dz != 1 {
                continue;
            }
            let lx = dims[0] as f64 / dx as f64;
            let ly = dims[1] as f64 / dy as f64;
            let lz = if d == 2 {
                1.0
            } else {
                dims[2] as f64 / dz as f64
            };
            let surface = if d == 2 {
                2.0 * (lx + ly)
            } else {
                2.0 * (lx * ly + lx * lz + ly * lz)
            };
            let balance = lx.max(ly).max(lz) / lx.min(ly).min(lz);
            let score = surface + balance * 1e-9;
            let decomp = [dx, dy, dz];
            if score < best_score
                || (score == best_score && decomp.iter().rev().cmp(best.iter().rev()).is_lt())
            {
                best_score = score;
                best = decomp;
            }
        }
    }
    best
}

fn fnv1a64(bytes: &[u8]) -> u64 {
    let mut h = 0xcbf29ce484222325u64;
    for &b in bytes {
        h ^= b as u64;
        h = h.wrapping_mul(0x100000001b3);
    }
    h
}

fn hash_bytes(h: &mut u64, bytes: &[u8]) {
    for &b in bytes {
        *h ^= b as u64;
        *h = h.wrapping_mul(0x100000001b3);
    }
}

fn hash_f64(h: &mut u64, v: f64) {
    hash_bytes(h, &v.to_bits().to_le_bytes());
}

fn spec_item_hashes<T: Real>(
    spec: &GlobalSpec<T>,
    solid: &[bool],
    wall_u: &[[T; 3]],
) -> [(&'static str, u64); 8] {
    let mut dims = 0xcbf29ce484222325u64;
    for v in spec.dims {
        hash_bytes(&mut dims, &(v as u64).to_le_bytes());
    }
    let mut periodic = 0xcbf29ce484222325u64;
    for v in spec.periodic {
        hash_bytes(&mut periodic, &[v as u8]);
    }
    let mut collision = 0xcbf29ce484222325u64;
    match spec.collision {
        CollisionKind::Bgk => hash_bytes(&mut collision, &[0]),
        CollisionKind::Trt { magic } => {
            hash_bytes(&mut collision, &[1]);
            hash_f64(&mut collision, magic);
        }
    }
    let mut faces = 0xcbf29ce484222325u64;
    for bc in spec.faces {
        match bc {
            FaceBC::Closed => hash_bytes(&mut faces, &[0]),
            FaceBC::Velocity { u } => {
                hash_bytes(&mut faces, &[1]);
                for v in u {
                    hash_f64(&mut faces, v.as_f64());
                }
            }
            FaceBC::Pressure { rho } => {
                hash_bytes(&mut faces, &[2]);
                hash_f64(&mut faces, rho.as_f64());
            }
            FaceBC::Outflow => hash_bytes(&mut faces, &[3]),
            FaceBC::Convective { u_conv } => {
                hash_bytes(&mut faces, &[4]);
                hash_f64(&mut faces, u_conv.as_f64());
            }
        }
    }
    let mut force = 0xcbf29ce484222325u64;
    for v in spec.force {
        hash_f64(&mut force, v.as_f64());
    }
    let mut solid_mask = 0xcbf29ce484222325u64;
    hash_bytes(&mut solid_mask, &(solid.len() as u64).to_le_bytes());
    for &v in solid {
        hash_bytes(&mut solid_mask, &[v as u8]);
    }
    let mut wall_mask = 0xcbf29ce484222325u64;
    hash_bytes(&mut wall_mask, &(wall_u.len() as u64).to_le_bytes());
    for u in wall_u {
        for &v in u {
            hash_f64(&mut wall_mask, v.as_f64());
        }
    }
    [
        ("dims", dims),
        ("nu", fnv1a64(&spec.nu.to_bits().to_le_bytes())),
        ("periodic", periodic),
        ("collision", collision),
        ("faces", faces),
        ("force", force),
        ("solid mask", solid_mask),
        ("wall velocity mask", wall_mask),
    ]
}

fn assert_rank_specs_match<T: Real>(
    world: &SimpleCommunicator,
    spec: &GlobalSpec<T>,
    solid: &[bool],
    wall_u: &[[T; 3]],
) {
    for (name, hash) in spec_item_hashes(spec, solid, wall_u) {
        let local = [hash];
        let mut min = [0u64];
        let mut max = [0u64];
        world.all_reduce_into(&local, &mut min, SystemOperation::min());
        world.all_reduce_into(&local, &mut max, SystemOperation::max());
        assert_eq!(
            min[0], max[0],
            "MPI rank specification mismatch: {name} differs across ranks"
        );
    }
}

/// MPI implementation of [`HaloExchange`]: serves exactly the local part of a
/// [`Solver::new_local_part`] decomposition, interpreting subdomain neighbour
/// ids as ranks of its (duplicated) communicator.
pub struct MpiExchange<T: Real> {
    comm: SimpleCommunicator,
    rank: usize,
    buffers: RefCell<MpiBuffers<T>>,
}

impl<T: Real> MpiExchange<T> {
    /// Duplicate `world` (collective) and bind the exchange to it, isolating
    /// halo traffic from the caller's communicator.
    pub fn new(world: &SimpleCommunicator) -> Self {
        let comm = world.duplicate();
        let rank = comm.rank() as usize;
        Self {
            comm,
            rank,
            buffers: RefCell::new(MpiBuffers::default()),
        }
    }

    /// This process's rank (= part id).
    pub fn rank(&self) -> usize {
        self.rank
    }

    /// Communicator size (= part count).
    pub fn size(&self) -> usize {
        self.comm.size() as usize
    }

    /// One axis phase: hand the two side layers to the face neighbours and
    /// return what arrived for (low, high). `payloads[s]` is the message the
    /// neighbour behind face `Face::ALL[2·axis+s]` of *the receiver* unpacks;
    /// it therefore goes to this part's *opposite*-side neighbour. A periodic
    /// self-wrap short-circuits to a local hand-off (no MPI).
    fn transfer_axis<E: Copy + Default>(
        comm: &SimpleCommunicator,
        me: usize,
        sub: &Subdomain,
        axis: usize,
        tag_base: Tag,
        buffers: &mut AxisBuffers<E>,
        recv_counts: [usize; 2],
        mut unpack: impl FnMut(usize, &[E]),
    ) {
        let faces = [Face::ALL[2 * axis], Face::ALL[2 * axis + 1]];
        let nb = [
            sub.neighbors[faces[0].index()],
            sub.neighbors[faces[1].index()],
        ];
        if nb[0] == Some(me) || nb[1] == Some(me) {
            // decomp == 1 on a periodic axis: both faces wrap onto this rank.
            assert!(
                nb[0] == Some(me) && nb[1] == Some(me),
                "self-wrap must be symmetric on axis {axis}"
            );
            debug_assert_eq!(buffers.send[0].len(), recv_counts[0]);
            debug_assert_eq!(buffers.send[1].len(), recv_counts[1]);
            unpack(0, &buffers.send[0]);
            unpack(1, &buffers.send[1]);
            return;
        }
        for s in 0..2 {
            if nb[s].is_some() {
                buffers.recv[s].resize(recv_counts[s], E::default());
            }
        }
        let (recv_lo, recv_hi) = buffers.recv.split_at_mut(1);
        let recv_lo = &mut recv_lo[0];
        let recv_hi = &mut recv_hi[0];
        let send = &buffers.send;
        mpi::request::scope(|sc| {
            let mut reqs = Vec::with_capacity(4);
            // Receives first, then sends; all are posted before any wait, so
            // the phase cannot deadlock regardless of the neighbour graph.
            if let Some(r) = nb[0] {
                reqs.push(
                    comm.process_at_rank(r as Rank)
                        .immediate_receive_into_with_tag(
                            sc,
                            as_bytes_mut(recv_lo.as_mut_slice()),
                            tag_base + faces[0].index() as Tag,
                        ),
                );
            }
            if let Some(r) = nb[1] {
                reqs.push(
                    comm.process_at_rank(r as Rank)
                        .immediate_receive_into_with_tag(
                            sc,
                            as_bytes_mut(recv_hi.as_mut_slice()),
                            tag_base + faces[1].index() as Tag,
                        ),
                );
            }
            if let Some(r) = nb[1] {
                // The +side neighbour unpacks this at its low face.
                reqs.push(comm.process_at_rank(r as Rank).immediate_send_with_tag(
                    sc,
                    as_bytes(send[0].as_slice()),
                    tag_base + faces[0].index() as Tag,
                ));
            }
            if let Some(r) = nb[0] {
                reqs.push(comm.process_at_rank(r as Rank).immediate_send_with_tag(
                    sc,
                    as_bytes(send[1].as_slice()),
                    tag_base + faces[1].index() as Tag,
                ));
            }
            for req in reqs {
                req.wait();
            }
        });
        for s in 0..2 {
            if nb[s].is_some() {
                unpack(s, &buffers.recv[s]);
            }
        }
    }
}

impl<T: Real> HaloExchange<T> for MpiExchange<T> {
    const SCOPE: ExchangeScope = ExchangeScope::Remote;

    fn exchange_f<L: Lattice>(&self, subs: &[Subdomain], parts: &mut [SoaFields<T>]) {
        assert_eq!(parts.len(), 1, "MpiExchange serves exactly the local part");
        let sub = &subs[0];
        let mut buffers = self.buffers.borrow_mut();
        for axis in 0..sub.geom.d {
            let faces = [Face::ALL[2 * axis], Face::ALL[2 * axis + 1]];
            let mut counts = [0usize; 2];
            for s in 0..2 {
                pack_f_layer::<L, T>(&parts[0], faces[s], &mut buffers.f[axis].send[s]);
                counts[s] =
                    layer_cell_count(&parts[0].geom, faces[s]) * L::unknowns(faces[s]).len();
            }
            Self::transfer_axis(
                &self.comm,
                self.rank,
                sub,
                axis,
                TAG_F,
                &mut buffers.f[axis],
                counts,
                |s, buf| unpack_f_layer::<L, T>(&mut parts[0], faces[s], buf),
            );
        }
    }

    fn exchange_masks(&self, subs: &[Subdomain], parts: &mut [SoaFields<T>]) {
        assert_eq!(parts.len(), 1, "MpiExchange serves exactly the local part");
        let sub = &subs[0];
        let geom = parts[0].geom;
        let mut buffers = self.buffers.borrow_mut();
        for axis in 0..sub.geom.d {
            let faces = [Face::ALL[2 * axis], Face::ALL[2 * axis + 1]];
            // Round A (u8): [solid n][probe-flag 1][probe n] — fixed size so
            // the receiver can post the buffer without a size handshake.
            let mut meta_counts = [0usize; 2];
            for s in 0..2 {
                let idx = layer_indices(&geom, faces[s], axis, false);
                let n = idx.len();
                let buf = &mut buffers.mask_meta[axis].send[s];
                buf.clear();
                buf.reserve(2 * n + 1);
                for &c in &idx {
                    buf.push(parts[0].solid[c] as u8);
                }
                buf.push(parts[0].probe.is_some() as u8);
                match &parts[0].probe {
                    Some(m) => buf.extend(idx.iter().map(|&c| m[c] as u8)),
                    None => buf.extend(std::iter::repeat(0u8).take(n)),
                }
                meta_counts[s] = 2 * layer_cell_count(&geom, faces[s]) + 1;
            }
            Self::transfer_axis(
                &self.comm,
                self.rank,
                sub,
                axis,
                TAG_MASK_META,
                &mut buffers.mask_meta[axis],
                meta_counts,
                |_, _| {},
            );
            // Round B (T): wall_u, cell-major with the 3 components inner.
            let mut wall_counts = [0usize; 2];
            for s in 0..2 {
                let idx = layer_indices(&geom, faces[s], axis, false);
                let buf = &mut buffers.mask_wall[axis].send[s];
                buf.clear();
                buf.reserve(3 * idx.len());
                for &c in &idx {
                    buf.extend_from_slice(&parts[0].wall_u[c]);
                }
                wall_counts[s] = 3 * layer_cell_count(&geom, faces[s]);
            }
            Self::transfer_axis(
                &self.comm,
                self.rank,
                sub,
                axis,
                TAG_MASK_WALL,
                &mut buffers.mask_wall[axis],
                wall_counts,
                |_, _| {},
            );
            for s in 0..2 {
                if sub.neighbors[faces[s].index()].is_none() {
                    continue;
                }
                let mbuf = if sub.neighbors[faces[s].index()] == Some(self.rank) {
                    &buffers.mask_meta[axis].send[s]
                } else {
                    &buffers.mask_meta[axis].recv[s]
                };
                let wbuf = if sub.neighbors[faces[s].index()] == Some(self.rank) {
                    &buffers.mask_wall[axis].send[s]
                } else {
                    &buffers.mask_wall[axis].recv[s]
                };
                let idx = layer_indices(&geom, faces[s], axis, true);
                let n = idx.len();
                for (k, &c) in idx.iter().enumerate() {
                    parts[0].solid[c] = mbuf[k] != 0;
                    parts[0].wall_u[c] = [wbuf[3 * k], wbuf[3 * k + 1], wbuf[3 * k + 2]];
                }
                if mbuf[n] != 0 {
                    let probe = parts[0]
                        .probe
                        .as_mut()
                        .expect("probe mask must be materialised on every rank");
                    for (k, &c) in idx.iter().enumerate() {
                        probe[c] = mbuf[n + 1 + k] != 0;
                    }
                }
            }
        }
    }

    fn exchange_scalar(&self, subs: &[Subdomain], planes: &mut [&mut [T]]) {
        assert_eq!(planes.len(), 1, "MpiExchange serves exactly the local part");
        let sub = &subs[0];
        let geom = &sub.geom;
        let mut buffers = self.buffers.borrow_mut();
        for axis in 0..sub.geom.d {
            let faces = [Face::ALL[2 * axis], Face::ALL[2 * axis + 1]];
            let mut counts = [0usize; 2];
            for s in 0..2 {
                pack_scalar_layer(geom, planes[0], faces[s], &mut buffers.scalar[axis].send[s]);
                counts[s] = layer_cell_count(geom, faces[s]);
            }
            Self::transfer_axis(
                &self.comm,
                self.rank,
                sub,
                axis,
                TAG_SCALAR,
                &mut buffers.scalar[axis],
                counts,
                |s, buf| unpack_scalar_layer(geom, planes[0], faces[s], buf),
            );
        }
    }
}

/// Distributed driver: rank `r` of the communicator runs part `r` of the
/// Cartesian decomposition through a [`Solver`] wired to [`MpiExchange`],
/// with the V1 diagnostics reduced across ranks and rank-0 field gathers.
///
/// See the module docs for the collectives contract.
pub struct MpiSolver<L, T, B>
where
    L: Lattice,
    T: Real,
    B: Backend<L, T, Fields = SoaFields<T>>,
{
    inner: Solver<L, T, B, MpiExchange<T>>,
    comm: SimpleCommunicator,
    rank: usize,
    size: usize,
    /// Full decomposition metadata (origins/extents of every rank's part).
    subs_meta: Vec<Subdomain>,
    probe_active: bool,
    probed_force: [T; 3],
}

impl<L, T, B> MpiSolver<L, T, B>
where
    L: Lattice,
    T: Real,
    B: Backend<L, T, Fields = SoaFields<T>>,
{
    /// Build the local part of `decomp` on this rank (collective).
    /// `decomp[0]·decomp[1]·decomp[2]` must equal the communicator size;
    /// `solid` / `wall_u` are the global compact arrays (see [`Solver::new`]).
    pub fn new(
        world: &SimpleCommunicator,
        spec: &GlobalSpec<T>,
        solid: &[bool],
        wall_u: &[[T; 3]],
        decomp: [usize; 3],
        backend: B,
    ) -> Self {
        let size = world.size() as usize;
        let rank = world.rank() as usize;
        assert_eq!(
            decomp[0] * decomp[1] * decomp[2],
            size,
            "decomp {decomp:?} must cover exactly the communicator size {size}"
        );
        assert_rank_specs_match(world, spec, solid, wall_u);
        let exchange = MpiExchange::<T>::new(world);
        let comm = world.duplicate();
        let subs_meta = partition(L::D, spec.dims, spec.periodic, decomp);
        let inner = Solver::new_local_part(spec, solid, wall_u, decomp, rank, backend, exchange);
        Self {
            inner,
            comm,
            rank,
            size,
            subs_meta,
            probe_active: false,
            probed_force: [T::zero(); 3],
        }
    }

    /// This process's rank (= part id).
    pub fn rank(&self) -> usize {
        self.rank
    }

    /// Number of ranks (= parts).
    pub fn size(&self) -> usize {
        self.size
    }

    /// The local [`Solver`] (one part). Cell accessors on it address global
    /// coordinates and are only valid for cells this rank owns.
    pub fn local(&self) -> &Solver<L, T, B, MpiExchange<T>> {
        &self.inner
    }

    /// Mutable access to the local solver (see [`MpiSolver::local`] caveats;
    /// mask edits through this must be followed by [`Solver::mark_masks_dirty`]
    /// *on every rank*).
    pub fn local_mut(&mut self) -> &mut Solver<L, T, B, MpiExchange<T>> {
        &mut self.inner
    }

    /// Whether this rank owns global cell `(x, y, z)`.
    pub fn owns(&self, x: usize, y: usize, z: usize) -> bool {
        let s = &self.subs_meta[self.rank];
        (0..3).all(|a| {
            let c = [x, y, z][a];
            c >= s.origin[a] && c < s.origin[a] + s.geom.core[a]
        })
    }

    /// Advance one step (collective). Refreshes the global probed force when
    /// a probe is active.
    pub fn step(&mut self) {
        self.inner.step();
        if self.probe_active {
            let lf = self.inner.probed_force();
            let local = [lf[0].as_f64(), lf[1].as_f64(), lf[2].as_f64()];
            let mut global = [0.0f64; 3];
            self.comm
                .all_reduce_into(&local[..], &mut global[..], SystemOperation::sum());
            self.probed_force = [T::r(global[0]), T::r(global[1]), T::r(global[2])];
        }
    }

    /// Advance `steps` steps (collective).
    pub fn run(&mut self, steps: usize) {
        for _ in 0..steps {
            self.step();
        }
    }

    /// Advance `steps` steps with a periodic non-finite watchdog (A-9):
    /// the distributed counterpart of `Solver::run_guarded`. Collective —
    /// every rank must call it with the same arguments; each check is one
    /// 2-double `MPI_SUM` Allreduce of the mass partials, so every rank sees
    /// the same global sum and takes the same `Ok`/`Err` branch (no
    /// divergence in control flow). `check_every == 0` is treated as 1.
    pub fn run_guarded(
        &mut self,
        steps: usize,
        check_every: usize,
    ) -> Result<(), crate::solver::Diverged> {
        let check_every = check_every.max(1);
        let mut since_check = 0usize;
        for _ in 0..steps {
            self.step();
            since_check += 1;
            if since_check == check_every {
                since_check = 0;
                self.check_mass_finite()?;
            }
        }
        if since_check > 0 {
            self.check_mass_finite()?;
        }
        Ok(())
    }

    fn check_mass_finite(&self) -> Result<(), crate::solver::Diverged> {
        let (fluid, m) = self.inner.local_mass_partials();
        let mut out = [0.0f64; 2];
        self.allreduce_sum(&[fluid, m], &mut out);
        if (out[0] + out[1]).is_finite() {
            Ok(())
        } else {
            Err(crate::solver::Diverged {
                step: self.inner.time(),
            })
        }
    }

    /// Initialise from `(rho, u) = init(x, y, z)` (collective; global
    /// coordinates, evaluated only on owned cells — pass the same closure on
    /// every rank).
    pub fn init_with(&mut self, init: impl Fn(usize, usize, usize) -> (T, [T; 3])) {
        self.inner.init_with(init);
    }

    /// Mark a global cell solid. **Call on every rank** with the same
    /// coordinates: the owner stores it, everyone schedules the collective
    /// mask-halo refresh.
    pub fn set_solid(&mut self, x: usize, y: usize, z: usize) {
        if self.owns(x, y, z) {
            self.inner.set_solid(x, y, z);
        }
        self.inner.mark_masks_dirty();
    }

    /// Select probed solids (collective; same predicate on every rank).
    pub fn set_force_probe(&mut self, pred: impl Fn(usize, usize, usize) -> bool) {
        self.inner.set_force_probe(pred);
        self.probe_active = true;
    }

    /// Per-node inlet profile on a `Velocity` face (collective; every rank
    /// passes the full global profile and slices out its own face cells).
    pub fn set_inlet_profile(&mut self, face: Face, values: &[[T; 3]]) {
        self.inner.set_inlet_profile(face, values);
    }

    /// Closure form of [`MpiSolver::set_inlet_profile`].
    pub fn set_inlet_profile_with(&mut self, face: Face, profile: impl Fn(usize, usize) -> [T; 3]) {
        self.inner.set_inlet_profile_with(face, profile);
    }

    /// Single-component Shan–Chen force refresh (collective): ψ halos travel
    /// through [`MpiExchange::exchange_scalar`].
    pub fn update_shan_chen_force(&mut self, g: T, psi: impl Fn(T) -> T) {
        self.inner.update_shan_chen_force(g, psi);
    }

    /// Wall-adhesion Shan–Chen refresh (collective); see
    /// [`Solver::update_shan_chen_force_with_walls`].
    pub fn update_shan_chen_force_with_walls(
        &mut self,
        g: T,
        g_wall: T,
        psi_wall: T,
        psi: impl Fn(T) -> T,
    ) {
        self.inner
            .update_shan_chen_force_with_walls(g, g_wall, psi_wall, psi);
    }

    /// Toggle two-pass streaming on the local part.
    pub fn set_two_pass(&mut self, on: bool) {
        self.inner.set_two_pass(on);
    }

    /// Completed steps.
    pub fn time(&self) -> u64 {
        self.inner.time()
    }

    /// Global grid extents.
    pub fn dims(&self) -> [usize; 3] {
        self.inner.dims()
    }

    /// Barrier over the solver's communicator (bench timing fences).
    pub fn barrier(&self) {
        self.comm.barrier();
    }

    // ------------------------------------------------------------------
    // Distributed diagnostics (Allreduce over rank partial sums)
    // ------------------------------------------------------------------

    fn allreduce_sum(&self, local: &[f64], global: &mut [f64]) {
        self.comm
            .all_reduce_into(local, global, SystemOperation::sum());
    }

    /// Global total mass (collective): rank partials of V1's `(fluid_cells,
    /// mass_deviation)` f64 sums, `MPI_SUM`-combined, then `fluid + m`.
    /// Differs from a single-rank run by f64 reassociation only.
    pub fn total_mass(&self) -> T {
        let (fluid, m) = self.inner.local_mass_partials();
        let mut out = [0.0f64; 2];
        self.allreduce_sum(&[fluid, m], &mut out);
        T::r(out[0] + out[1])
    }

    /// Global total mass with fixed-order composition. Rank blocks are first
    /// reduced to global `(z, y)` row partials, rank 0 folds rows in
    /// lexicographic order, and the scalar result is broadcast to all ranks.
    pub fn total_mass_deterministic(&self) -> T {
        let dims = self.inner.dims();
        let rows = dims[1] * dims[2];
        let mut local = vec![0.0f64; rows * 2];
        let sub = &self.subs_meta[self.rank];
        let fields = self.inner.fields(0);
        let g = fields.geom;
        let np = g.n_padded();
        for z in 0..g.core[2] {
            for y in 0..g.core[1] {
                let row = (sub.origin[2] + z) * dims[1] + sub.origin[1] + y;
                for x in 0..g.core[0] {
                    let pi = g.pidx(x, y, z);
                    if fields.solid[pi] {
                        continue;
                    }
                    local[2 * row] += 1.0;
                    for q in 0..L::Q {
                        local[2 * row + 1] += fields.f[q * np + pi].as_f64();
                    }
                }
            }
        }

        let mut total = 0.0f64;
        if self.rank == 0 {
            let mut rows_by_rank = vec![0.0f64; rows * 2];
            for row in 0..rows {
                rows_by_rank[2 * row] = local[2 * row];
                rows_by_rank[2 * row + 1] = local[2 * row + 1];
            }
            let mut staging = vec![0.0f64; rows * 2];
            let mut all = vec![vec![0.0f64; rows * 2]; self.size];
            all[0].copy_from_slice(&rows_by_rank);
            for r in 1..self.size {
                staging.fill(0.0);
                self.comm
                    .process_at_rank(r as Rank)
                    .receive_into_with_tag(as_bytes_mut(staging.as_mut_slice()), TAG_MASS_ROWS);
                all[r].copy_from_slice(&staging);
            }
            for row in 0..rows {
                let mut fluid = 0.0f64;
                let mut mass_dev = 0.0f64;
                for rank_rows in &all {
                    fluid += rank_rows[2 * row];
                    mass_dev += rank_rows[2 * row + 1];
                }
                total += fluid + mass_dev;
            }
        } else {
            self.comm
                .process_at_rank(0)
                .send_with_tag(as_bytes(local.as_slice()), TAG_MASS_ROWS);
        }
        self.comm.process_at_rank(0).broadcast_into(&mut total);
        T::r(total)
    }

    /// Global total momentum (collective; same contract as `total_mass`).
    pub fn total_momentum(&self) -> [T; 3] {
        let local = self.inner.local_momentum_partials();
        let mut out = [0.0f64; 3];
        self.allreduce_sum(&local, &mut out);
        [T::r(out[0]), T::r(out[1]), T::r(out[2])]
    }

    /// Momentum-exchange force on probed solids during the most recent step,
    /// summed over all ranks (valid after a `step` with an active probe).
    pub fn probed_force(&self) -> [T; 3] {
        self.probed_force
    }

    /// Global count of non-finite values in populations/moments (collective).
    /// `0` on a healthy run — the distributed NaN check.
    pub fn nonfinite_count(&self) -> u64 {
        let local = [self.inner.local_nonfinite_count() as f64];
        let mut out = [0.0f64];
        self.allreduce_sum(&local, &mut out);
        out[0] as u64
    }

    /// Global fluid-cell count (collective).
    pub fn fluid_cell_count(&self) -> usize {
        let local = [self.inner.local_mass_partials().0];
        let mut out = [0.0f64];
        self.allreduce_sum(&local, &mut out);
        out[0] as usize
    }

    // ------------------------------------------------------------------
    // Rank-0 gathers (output / verification)
    // ------------------------------------------------------------------

    /// Assemble a global compact field on rank 0 from per-rank compact-core
    /// blocks (collective). `Some(global)` on rank 0, `None` elsewhere.
    fn gather_compact(&self, local: &[T]) -> Option<Vec<T>> {
        let nc = self.subs_meta[self.rank].geom.n_core();
        assert_eq!(local.len(), nc, "local block must be the compact core");
        if self.rank != 0 {
            self.comm
                .process_at_rank(0)
                .send_with_tag(as_bytes(local), TAG_GATHER);
            return None;
        }
        let dims = self.inner.dims();
        let mut out = vec![T::zero(); dims[0] * dims[1] * dims[2]];
        let mut staging: Vec<T> = Vec::new();
        for r in 0..self.size {
            let sub = &self.subs_meta[r];
            let g = sub.geom;
            let block: &[T] = if r == 0 {
                local
            } else {
                staging.clear();
                staging.resize(g.n_core(), T::zero());
                self.comm
                    .process_at_rank(r as Rank)
                    .receive_into_with_tag(as_bytes_mut(staging.as_mut_slice()), TAG_GATHER);
                &staging
            };
            for z in 0..g.core[2] {
                for y in 0..g.core[1] {
                    for x in 0..g.core[0] {
                        let gi = ((sub.origin[2] + z) * dims[1] + (sub.origin[1] + y)) * dims[0]
                            + (sub.origin[0] + x);
                        out[gi] = block[g.cidx(x, y, z)];
                    }
                }
            }
        }
        Some(out)
    }

    /// Assemble a global compact vector with `ncomp` scalar components per
    /// cell on rank 0. `local` is cell-major: `cell * ncomp + component`.
    fn gather_compact_components(&self, local: &[T], ncomp: usize) -> Option<Vec<T>> {
        let nc = self.subs_meta[self.rank].geom.n_core();
        assert_eq!(
            local.len(),
            nc * ncomp,
            "local block must be compact core times component count"
        );
        if self.rank != 0 {
            self.comm
                .process_at_rank(0)
                .send_with_tag(as_bytes(local), TAG_GATHER);
            return None;
        }
        let dims = self.inner.dims();
        let mut out = vec![T::zero(); dims[0] * dims[1] * dims[2] * ncomp];
        let mut staging: Vec<T> = Vec::new();
        for r in 0..self.size {
            let sub = &self.subs_meta[r];
            let g = sub.geom;
            let block: &[T] = if r == 0 {
                local
            } else {
                staging.clear();
                staging.resize(g.n_core() * ncomp, T::zero());
                self.comm
                    .process_at_rank(r as Rank)
                    .receive_into_with_tag(as_bytes_mut(staging.as_mut_slice()), TAG_GATHER);
                &staging
            };
            for z in 0..g.core[2] {
                for y in 0..g.core[1] {
                    for x in 0..g.core[0] {
                        let gi = ((sub.origin[2] + z) * dims[1] + (sub.origin[1] + y)) * dims[0]
                            + (sub.origin[0] + x);
                        let lc = g.cidx(x, y, z);
                        let dst = gi * ncomp;
                        let src = lc * ncomp;
                        out[dst..dst + ncomp].copy_from_slice(&block[src..src + ncomp]);
                    }
                }
            }
        }
        Some(out)
    }

    /// Global density field on rank 0 (collective).
    pub fn gather_rho(&self) -> Option<Vec<T>> {
        self.gather_compact(&self.inner.fields(0).rho)
    }
    /// Global x-velocity field on rank 0 (collective).
    pub fn gather_ux(&self) -> Option<Vec<T>> {
        self.gather_compact(&self.inner.fields(0).ux)
    }
    /// Global y-velocity field on rank 0 (collective).
    pub fn gather_uy(&self) -> Option<Vec<T>> {
        self.gather_compact(&self.inner.fields(0).uy)
    }
    /// Global z-velocity field on rank 0 (collective).
    pub fn gather_uz(&self) -> Option<Vec<T>> {
        self.gather_compact(&self.inner.fields(0).uz)
    }

    /// Global strain-rate tensor on rank 0 (collective).
    ///
    /// Components are `[S_xx, S_yy, S_zz, S_xy, S_xz, S_yz]`; solid cells are
    /// zeros. See [`Solver::gather_strain_rate`] for the stage convention,
    /// force correction and TRT `omega_plus` relaxation note.
    pub fn gather_strain_rate(&self) -> Option<Vec<[T; 6]>> {
        let full = self.inner.gather_strain_rate();
        let sub = &self.subs_meta[self.rank];
        let g = sub.geom;
        let dims = self.inner.dims();
        let mut local = vec![T::zero(); g.n_core() * 6];
        for z in 0..g.core[2] {
            for y in 0..g.core[1] {
                for x in 0..g.core[0] {
                    let gi = ((sub.origin[2] + z) * dims[1] + (sub.origin[1] + y)) * dims[0]
                        + (sub.origin[0] + x);
                    let dst = g.cidx(x, y, z) * 6;
                    local[dst..dst + 6].copy_from_slice(&full[gi]);
                }
            }
        }
        self.gather_compact_components(&local, 6).map(|flat| {
            flat.chunks_exact(6)
                .map(|s| [s[0], s[1], s[2], s[3], s[4], s[5]])
                .collect()
        })
    }

    /// Global shear-rate invariant on rank 0 (collective).
    pub fn gather_shear_rate(&self) -> Option<Vec<T>> {
        let full = self.inner.gather_shear_rate();
        let sub = &self.subs_meta[self.rank];
        let g = sub.geom;
        let dims = self.inner.dims();
        let mut local = vec![T::zero(); g.n_core()];
        for z in 0..g.core[2] {
            for y in 0..g.core[1] {
                for x in 0..g.core[0] {
                    let gi = ((sub.origin[2] + z) * dims[1] + (sub.origin[1] + y)) * dims[0]
                        + (sub.origin[0] + x);
                    local[g.cidx(x, y, z)] = full[gi];
                }
            }
        }
        self.gather_compact(&local)
    }

    /// Global deviation-population plane `q` on rank 0 (collective) — the
    /// strongest equality statement for verification.
    pub fn gather_f(&self, q: usize) -> Option<Vec<T>> {
        let fields = self.inner.fields(0);
        let g = fields.geom;
        let np = g.n_padded();
        let mut local = vec![T::zero(); g.n_core()];
        for z in 0..g.core[2] {
            for y in 0..g.core[1] {
                for x in 0..g.core[0] {
                    local[g.cidx(x, y, z)] = fields.f[q * np + g.pidx(x, y, z)];
                }
            }
        }
        self.gather_compact(&local)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::backend::CpuScalar;
    use crate::halo::InProcess;
    use crate::lattice::{D2Q9, D3Q19};
    use crate::solver::WallSpec;

    fn patterned_parts<L: Lattice>(dims: [usize; 3], decomp: [usize; 3]) -> Vec<SoaFields<f64>> {
        let subs = partition(L::D, dims, [true, true, true], decomp);
        let mut parts: Vec<SoaFields<f64>> =
            subs.iter().map(|s| SoaFields::new(L::Q, s.geom)).collect();
        for (p, fields) in parts.iter_mut().enumerate() {
            let np = fields.plane_len();
            for q in 0..L::Q {
                for cell in 0..np {
                    fields.f[q * np + cell] = (p * 10_000 + q * 1_000 + cell) as f64;
                }
            }
        }
        parts
    }

    #[test]
    fn choose_decomp_minimizes_surface_for_common_rank_counts() {
        assert_eq!(choose_decomp(2, [96, 64, 1], 3), [3, 1, 1]);
        assert_eq!(choose_decomp(2, [96, 64, 1], 5), [5, 1, 1]);
        assert_eq!(choose_decomp(2, [96, 64, 1], 6), [3, 2, 1]);
        assert_eq!(choose_decomp(3, [30, 24, 18], 6), [3, 2, 1]);
    }

    #[test]
    fn mpi_phase_payloads_match_inprocess_exchange_buffers() {
        let dims = [8usize, 9, 1];
        let decomp = [2usize, 3, 1];
        let subs = partition(D2Q9::D, dims, [true, true, false], decomp);
        let before = patterned_parts::<D2Q9>(dims, decomp);
        let mut inprocess = before.clone();
        InProcess.exchange_f::<D2Q9>(&subs, &mut inprocess);

        let mut rebuilt = before.clone();
        for axis in 0..D2Q9::D {
            let phase_start = rebuilt.clone();
            for side in 0..2 {
                let recv_face = Face::ALL[2 * axis + side];
                for dst in 0..rebuilt.len() {
                    let Some(src) = subs[dst].neighbors[recv_face.index()] else {
                        continue;
                    };
                    let mut buf = Vec::new();
                    pack_f_layer::<D2Q9, f64>(&phase_start[src], recv_face, &mut buf);
                    let want_len = layer_cell_count(&before[dst].geom, recv_face)
                        * D2Q9::unknowns(recv_face).len();
                    assert_eq!(buf.len(), want_len);
                    unpack_f_layer::<D2Q9, f64>(&mut rebuilt[dst], recv_face, &buf);
                }
            }
        }
        for (a, b) in rebuilt.iter().zip(&inprocess) {
            assert_eq!(a.f, b.f);
        }
    }

    #[test]
    fn scalar_phase_payloads_match_inprocess_exchange_buffers() {
        let dims = [6usize, 4, 6];
        let decomp = [2usize, 1, 2];
        let subs = partition(D3Q19::D, dims, [true, true, true], decomp);
        let mut planes: Vec<Vec<f64>> = subs
            .iter()
            .enumerate()
            .map(|(p, s)| {
                (0..s.geom.n_padded())
                    .map(|cell| (p * 10_000 + cell) as f64)
                    .collect()
            })
            .collect();
        let before = planes.clone();
        let mut refs: Vec<&mut [f64]> = planes.iter_mut().map(|p| p.as_mut_slice()).collect();
        InProcess.exchange_scalar(&subs, &mut refs);

        let mut rebuilt = before.clone();
        for axis in 0..D3Q19::D {
            let phase_start = rebuilt.clone();
            for side in 0..2 {
                let recv_face = Face::ALL[2 * axis + side];
                for dst in 0..rebuilt.len() {
                    let Some(src) = subs[dst].neighbors[recv_face.index()] else {
                        continue;
                    };
                    let mut buf = Vec::new();
                    pack_scalar_layer(&subs[src].geom, &phase_start[src], recv_face, &mut buf);
                    assert_eq!(buf.len(), layer_cell_count(&subs[dst].geom, recv_face));
                    unpack_scalar_layer(&subs[dst].geom, &mut rebuilt[dst], recv_face, &buf);
                }
            }
        }
        assert_eq!(rebuilt, planes);
    }

    #[test]
    fn rank_spec_hash_names_changed_item() {
        let spec = GlobalSpec::<f64>::default();
        let mut changed = spec.clone();
        changed.nu = 0.03;
        let base = spec_item_hashes(&spec, &[], &[]);
        let other = spec_item_hashes(&changed, &[], &[]);
        let changed_names: Vec<&str> = base
            .iter()
            .zip(other.iter())
            .filter_map(|((name, a), (_, b))| (a != b).then_some(*name))
            .collect();
        assert_eq!(changed_names, vec!["nu"]);
    }

    #[test]
    fn one_rank_self_exchange_smoke_if_mpi_is_available() {
        let Some(universe) = mpi::initialize() else {
            return;
        };
        let world = universe.world();
        if world.size() != 1 {
            return;
        }
        let spec = GlobalSpec::<f64> {
            dims: [8, 6, 1],
            nu: 0.02,
            periodic: [true, true, false],
            ..Default::default()
        };
        let (solid, wall_u) =
            crate::solver::build_wall_rims(D2Q9::D, spec.dims, &WallSpec::default());
        let mut solver: MpiSolver<D2Q9, f64, CpuScalar> = MpiSolver::new(
            &world,
            &spec,
            &solid,
            &wall_u,
            [1, 1, 1],
            CpuScalar::default(),
        );
        solver.step();
        assert_eq!(solver.nonfinite_count(), 0);
        drop(solver);
        drop(world);
        drop(universe);
    }
}
