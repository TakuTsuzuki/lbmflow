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
use mpi::traits::{Communicator, CommunicatorCollectives, Destination, Source};
use mpi::{Rank, Tag};

use crate::backend::Backend;
use crate::fields::SoaFields;
use crate::halo::{
    layer_cell_count, layer_indices, pack_f_layer, pack_scalar_layer, unpack_f_layer,
    unpack_scalar_layer, ExchangeScope, HaloExchange,
};
use crate::lattice::{Face, Lattice};
use crate::real::Real;
use crate::solver::{partition, GlobalSpec, Solver};
use crate::subdomain::Subdomain;

const TAG_F: Tag = 100;
const TAG_SCALAR: Tag = 200;
const TAG_MASK_META: Tag = 300;
const TAG_MASK_WALL: Tag = 400;
const TAG_GATHER: Tag = 500;

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

/// MPI implementation of [`HaloExchange`]: serves exactly the local part of a
/// [`Solver::new_local_part`] decomposition, interpreting subdomain neighbour
/// ids as ranks of its (duplicated) communicator.
pub struct MpiExchange {
    comm: SimpleCommunicator,
    rank: usize,
}

impl MpiExchange {
    /// Duplicate `world` (collective) and bind the exchange to it, isolating
    /// halo traffic from the caller's communicator.
    pub fn new(world: &SimpleCommunicator) -> Self {
        let comm = world.duplicate();
        let rank = comm.rank() as usize;
        Self { comm, rank }
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
        &self,
        sub: &Subdomain,
        axis: usize,
        tag_base: Tag,
        payloads: [Vec<E>; 2],
        recv_counts: [usize; 2],
    ) -> [Option<Vec<E>>; 2] {
        let faces = [Face::ALL[2 * axis], Face::ALL[2 * axis + 1]];
        let nb = [
            sub.neighbors[faces[0].index()],
            sub.neighbors[faces[1].index()],
        ];
        let me = self.rank;
        if nb[0] == Some(me) || nb[1] == Some(me) {
            // decomp == 1 on a periodic axis: both faces wrap onto this rank.
            assert!(
                nb[0] == Some(me) && nb[1] == Some(me),
                "self-wrap must be symmetric on axis {axis}"
            );
            let [lo, hi] = payloads;
            debug_assert_eq!(lo.len(), recv_counts[0]);
            debug_assert_eq!(hi.len(), recv_counts[1]);
            return [Some(lo), Some(hi)];
        }
        let mut rb_lo: Option<Vec<E>> = nb[0].map(|_| vec![E::default(); recv_counts[0]]);
        let mut rb_hi: Option<Vec<E>> = nb[1].map(|_| vec![E::default(); recv_counts[1]]);
        let (pay_lo, pay_hi) = (&payloads[0], &payloads[1]);
        mpi::request::scope(|sc| {
            let mut reqs = Vec::with_capacity(4);
            // Receives first, then sends; all are posted before any wait, so
            // the phase cannot deadlock regardless of the neighbour graph.
            if let (Some(r), Some(buf)) = (nb[0], rb_lo.as_mut()) {
                reqs.push(
                    self.comm
                        .process_at_rank(r as Rank)
                        .immediate_receive_into_with_tag(
                            sc,
                            as_bytes_mut(buf.as_mut_slice()),
                            tag_base + faces[0].index() as Tag,
                        ),
                );
            }
            if let (Some(r), Some(buf)) = (nb[1], rb_hi.as_mut()) {
                reqs.push(
                    self.comm
                        .process_at_rank(r as Rank)
                        .immediate_receive_into_with_tag(
                            sc,
                            as_bytes_mut(buf.as_mut_slice()),
                            tag_base + faces[1].index() as Tag,
                        ),
                );
            }
            if let Some(r) = nb[1] {
                // The +side neighbour unpacks this at its low face.
                reqs.push(
                    self.comm
                        .process_at_rank(r as Rank)
                        .immediate_send_with_tag(
                            sc,
                            as_bytes(pay_lo.as_slice()),
                            tag_base + faces[0].index() as Tag,
                        ),
                );
            }
            if let Some(r) = nb[0] {
                reqs.push(
                    self.comm
                        .process_at_rank(r as Rank)
                        .immediate_send_with_tag(
                            sc,
                            as_bytes(pay_hi.as_slice()),
                            tag_base + faces[1].index() as Tag,
                        ),
                );
            }
            for req in reqs {
                req.wait();
            }
        });
        [rb_lo, rb_hi]
    }
}

impl<T: Real> HaloExchange<T> for MpiExchange {
    const SCOPE: ExchangeScope = ExchangeScope::Remote;

    fn exchange_f<L: Lattice>(&self, subs: &[Subdomain], parts: &mut [SoaFields<T>]) {
        assert_eq!(parts.len(), 1, "MpiExchange serves exactly the local part");
        let sub = &subs[0];
        for axis in 0..sub.geom.d {
            let faces = [Face::ALL[2 * axis], Face::ALL[2 * axis + 1]];
            let mut payloads = [Vec::new(), Vec::new()];
            let mut counts = [0usize; 2];
            for s in 0..2 {
                pack_f_layer::<L, T>(&parts[0], faces[s], &mut payloads[s]);
                counts[s] =
                    layer_cell_count(&parts[0].geom, faces[s]) * L::unknowns(faces[s]).len();
            }
            let recvd = self.transfer_axis(sub, axis, TAG_F, payloads, counts);
            for s in 0..2 {
                if let Some(buf) = &recvd[s] {
                    unpack_f_layer::<L, T>(&mut parts[0], faces[s], buf);
                }
            }
        }
    }

    fn exchange_masks(&self, subs: &[Subdomain], parts: &mut [SoaFields<T>]) {
        assert_eq!(parts.len(), 1, "MpiExchange serves exactly the local part");
        let sub = &subs[0];
        let geom = parts[0].geom;
        for axis in 0..sub.geom.d {
            let faces = [Face::ALL[2 * axis], Face::ALL[2 * axis + 1]];
            // Round A (u8): [solid n][probe-flag 1][probe n] — fixed size so
            // the receiver can post the buffer without a size handshake.
            let mut meta = [Vec::new(), Vec::new()];
            let mut meta_counts = [0usize; 2];
            for s in 0..2 {
                let idx = layer_indices(&geom, faces[s], axis, false);
                let n = idx.len();
                let mut buf = Vec::with_capacity(2 * n + 1);
                for &c in &idx {
                    buf.push(parts[0].solid[c] as u8);
                }
                buf.push(parts[0].probe.is_some() as u8);
                match &parts[0].probe {
                    Some(m) => buf.extend(idx.iter().map(|&c| m[c] as u8)),
                    None => buf.extend(std::iter::repeat(0u8).take(n)),
                }
                meta[s] = buf;
                meta_counts[s] = 2 * layer_cell_count(&geom, faces[s]) + 1;
            }
            let meta_recv = self.transfer_axis(sub, axis, TAG_MASK_META, meta, meta_counts);
            // Round B (T): wall_u, cell-major with the 3 components inner.
            let mut wall = [Vec::new(), Vec::new()];
            let mut wall_counts = [0usize; 2];
            for s in 0..2 {
                let idx = layer_indices(&geom, faces[s], axis, false);
                let mut buf: Vec<T> = Vec::with_capacity(3 * idx.len());
                for &c in &idx {
                    buf.extend_from_slice(&parts[0].wall_u[c]);
                }
                wall[s] = buf;
                wall_counts[s] = 3 * layer_cell_count(&geom, faces[s]);
            }
            let wall_recv = self.transfer_axis(sub, axis, TAG_MASK_WALL, wall, wall_counts);
            for s in 0..2 {
                let (Some(mbuf), Some(wbuf)) = (&meta_recv[s], &wall_recv[s]) else {
                    continue;
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
        for axis in 0..sub.geom.d {
            let faces = [Face::ALL[2 * axis], Face::ALL[2 * axis + 1]];
            let mut payloads = [Vec::new(), Vec::new()];
            let mut counts = [0usize; 2];
            for s in 0..2 {
                pack_scalar_layer(geom, planes[0], faces[s], &mut payloads[s]);
                counts[s] = layer_cell_count(geom, faces[s]);
            }
            let recvd = self.transfer_axis(sub, axis, TAG_SCALAR, payloads, counts);
            for s in 0..2 {
                if let Some(buf) = &recvd[s] {
                    unpack_scalar_layer(geom, planes[0], faces[s], buf);
                }
            }
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
    inner: Solver<L, T, B, MpiExchange>,
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
        let exchange = MpiExchange::new(world);
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
    pub fn local(&self) -> &Solver<L, T, B, MpiExchange> {
        &self.inner
    }

    /// Mutable access to the local solver (see [`MpiSolver::local`] caveats;
    /// mask edits through this must be followed by [`Solver::mark_masks_dirty`]
    /// *on every rank*).
    pub fn local_mut(&mut self) -> &mut Solver<L, T, B, MpiExchange> {
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
