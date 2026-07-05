//! `CpuSimd`: the fused, band-parallel CPU backend — the V2 port of V1's
//! `step_band()` kernel (docs/ARCHITECTURE_V2.md §2.4, docs/PERFORMANCE.md).
//!
//! ## Phase mapping (same observable state evolution as `CpuScalar`)
//!
//! The `Backend` phase split exists because the halo exchange runs *between*
//! collision and streaming. The fused backend keeps that contract with a
//! minimal collide phase:
//!
//! - [`Backend::collide`] collides **only the one-cell core boundary shell of
//!   faces that have a halo neighbour** in place in `f` — exactly the layers
//!   `pack_f_layer` reads, so the exchange ships post-collide populations
//!   just as it does under `CpuScalar`. Everything else stays pre-collide.
//! - [`Backend::stream`] is the fused pass: destination slabs (2D: rows,
//!   3D: z-planes) are processed in independent bands; each band collides its
//!   source slabs just-in-time into a 3-slab cache-resident ring (skipping
//!   the shell cells already collided above), pull-streams from the ring
//!   (solid-free spans as plain shifted copies, solid runs as half-way
//!   bounce-back with moving-wall momentum injection), and writes the new
//!   moments into double buffers while in-flight collides keep reading the
//!   previous step's moments (V1 `rho2` mechanics). Band-edge source slabs
//!   are collided redundantly by the neighbouring band (same inputs, same
//!   results, no synchronisation).
//! - [`Backend::swap`] swaps the population ping-pong *and* the moment
//!   double buffers.
//! - [`Backend::apply_open_faces`] first restores V1's stale-slot convention
//!   (see below), then delegates to the shared BC pass, byte-for-byte the
//!   `CpuScalar` code.
//! - [`Backend::update_moments`] right after a fused pass only recomputes
//!   the open-face boundary cells (whose populations the BC pass just
//!   rewrote) — V1 `fix_boundary_moments`; any other call (build / init)
//!   falls back to the shared full recompute.
//!
//! ## Stale-slot convention (ConvectiveOutflow and friends)
//!
//! Under `CpuScalar`, collision is in-place, so the out-buffer slots that
//! streaming skips at open faces still hold the **previous step's
//! post-collide** populations after the swap — the memory term the
//! convective BC reads, and the values Zou–He would see at exotic
//! open∩open corners. The fused pass never materialises post-collide state
//! in `f`, so it captures those populations per open face during `stream`
//! (re-colliding the face cells from the untouched `f` + old moments — the
//! V1 `capture_conv_stale` generalisation) and `apply_open_faces` writes the
//! *previous* step's capture into the unknown slots before the BC pass.
//! Zero-initialised: bitwise-matches `CpuScalar`'s all-zero first-step
//! ping-pong buffer.
//!
//! ## Equivalence to `CpuScalar`
//!
//! Every per-cell floating-point expression replicates the corresponding
//! `kernels.rs` expression operand-for-operand (the lattice constants are
//! compile-time promoted, which folds the ±1/0 multiplies without changing
//! IEEE results), so fields agree with `CpuScalar` **bitwise up to the sign
//! of exact zeros** — with one documented exception: when forcing is off,
//! this backend skips `kernels::collide_row`'s `+ cp*0 + cm*0` no-op terms,
//! which can flip a `-0.0` to `+0.0`. No downstream arithmetic observes
//! either difference (asserted at |Δ| = 0 by `tests/backend_simd_equiv.rs`).
//! The probed-force diagnostic sums link contributions in a different order
//! (direction-major per row, V1 fused order) and may differ by f64/f32
//! reassociation only; it never feeds back into the fields.

use crate::backend::{
    apply_open_faces_impl, read_moments_impl, reduce_impl, update_moments_impl, Backend,
    CellRange, HostMoments, PARALLEL_MIN_CELLS,
};
use crate::fields::{FusedScratch, LocalGeom, SoaFields};
use crate::kernels::{for_face_cells, RawSlice};
use crate::lattice::{Face, Lattice, Q_MAX};
use crate::params::{KParams, Reduction, StepParams};
use crate::real::Real;
use crate::subdomain::Subdomain;

#[cfg(feature = "parallel")]
use rayon::prelude::*;

/// Fused band-parallel CPU backend (drop-in for `CpuScalar`; same
/// `SoaFields` storage, so it composes with every `HaloExchange` and with
/// `MpiSolver`).
#[derive(Clone, Copy, Debug)]
pub struct CpuSimd {
    /// Band-parallel dispatch threshold in core cells (V1
    /// `PARALLEL_MIN_CELLS` semantics).
    pub parallel_min_cells: usize,
}

impl Default for CpuSimd {
    fn default() -> Self {
        Self {
            parallel_min_cells: PARALLEL_MIN_CELLS,
        }
    }
}

impl CpuSimd {
    fn use_parallel(&self, sub: &Subdomain) -> bool {
        cfg!(feature = "parallel") && sub.geom.n_core() >= self.parallel_min_cells
    }
}

// ---------------------------------------------------------------------------
// Small helpers
// ---------------------------------------------------------------------------

/// `c_q · v` seeded with the `d = 0` term — exactly `kernels::collide_row`'s
/// association. The velocity components come from the compile-time lattice
/// table, so after unrolling the ±1/0 multiplies fold to moves/negations
/// with IEEE-identical results.
#[inline(always)]
fn dotc<L: Lattice, T: Real>(q: usize, v: [T; 3]) -> T {
    let mut acc = T::r(L::C[q][0] as f64) * v[0];
    for d in 1..L::D {
        acc = acc + T::r(L::C[q][d] as f64) * v[d];
    }
    acc
}

/// Maximal runs of consecutive solid cells over one padded row.
fn solid_runs_row(solid_row: &[bool], out: &mut Vec<(u32, u32)>) {
    out.clear();
    let n = solid_row.len();
    let mut x = 0;
    while x < n {
        if solid_row[x] {
            let start = x;
            while x < n && solid_row[x] {
                x += 1;
            }
            out.push((start as u32, x as u32));
        } else {
            x += 1;
        }
    }
}

/// Walk the fluid (run-complement) sub-spans of `[w0, w1)` in padded row
/// coordinates. `runs` must be sorted and disjoint (as produced by
/// [`solid_runs_row`]).
#[inline]
fn for_fluid_spans(
    runs: &[(u32, u32)],
    w0: usize,
    w1: usize,
    mut body: impl FnMut(usize, usize),
) {
    let mut cursor = w0;
    for &(a, b) in runs {
        let (a, b) = (a as usize, b as usize);
        if b <= w0 {
            continue;
        }
        if a >= w1 {
            break;
        }
        if cursor < a.min(w1) {
            body(cursor, a.min(w1));
        }
        cursor = cursor.max(b.min(w1));
    }
    if cursor < w1 {
        body(cursor, w1);
    }
}

// ---------------------------------------------------------------------------
// Span kernels (branch-free inner loops over solid-free spans)
// ---------------------------------------------------------------------------

/// TRT collision with Guo forcing over core cells `x0..x1` of one row stored
/// as `L::Q` planes at stride `q_stride` inside `planes` (`base` = index of
/// the row's core `x = 0` cell in plane 0).
///
/// Per-cell arithmetic replicates `kernels::collide_row`
/// expression-for-expression; the caller decomposes rows into solid-free
/// spans so this loop is branch-free and auto-vectorizes (V1 `collide_span`
/// shape). With `FORCE = false` the source terms are compiled out (skipping
/// `collide_row`'s exact-zero adds — a sign-of-zero-only deviation).
///
/// # Safety
/// The caller must be the only concurrent accessor of the addressed cells,
/// and `q * q_stride + base + x` must be in bounds for all `q < L::Q`,
/// `x0 <= x < x1`.
#[allow(clippy::too_many_arguments)]
unsafe fn collide_span_fused<L: Lattice, T: Real, const FORCE: bool, const FF: bool>(
    planes: RawSlice<T>,
    q_stride: usize,
    base: usize,
    x0: usize,
    x1: usize,
    rho: &[T],
    ux: &[T],
    uy: &[T],
    uz: &[T],
    field: &[[T; 3]],
    kp: &KParams<T>,
) {
    if x0 >= x1 {
        return;
    }
    let three = T::r(3.0);
    let f45 = T::r(4.5);
    let f15 = T::r(1.5);
    let nine = T::r(9.0);
    let half = T::r(0.5);
    let (op, om, cp, cm) = (kp.omega_p, kp.omega_m, kp.cp, kp.cm);
    for x in x0..x1 {
        let r = rho[x];
        let u = [ux[x], uy[x], uz[x]];
        let fv = if FORCE {
            if FF {
                [
                    kp.force[0] + field[x][0],
                    kp.force[1] + field[x][1],
                    kp.force[2] + field[x][2],
                ]
            } else {
                kp.force
            }
        } else {
            [T::zero(); 3]
        };
        let mut usq = u[0] * u[0];
        for d in 1..L::D {
            usq = usq + u[d] * u[d];
        }
        let uf = if FORCE {
            let mut uf = u[0] * fv[0];
            for d in 1..L::D {
                uf = uf + u[d] * fv[d];
            }
            uf
        } else {
            T::zero()
        };
        let drho = r - T::one();
        let i = base + x;
        // Rest direction: f0 - op*(f0 - feq0) + cp*src0.
        {
            let q = L::REST;
            let cu = dotc::<L, T>(q, u);
            let feq = T::r(L::W[q]) * (drho + r * (three * cu + f45 * cu * cu - f15 * usq));
            let i0 = q * q_stride + i;
            // SAFETY: caller contract (disjoint cells, in bounds).
            let f0 = unsafe { planes.get(i0) };
            if FORCE {
                let cf = dotc::<L, T>(q, fv);
                let src = T::r(L::W[q]) * (three * (cf - uf) + nine * cu * cf);
                unsafe { planes.set(i0, f0 - op * (f0 - feq) + cp * src) };
            } else {
                unsafe { planes.set(i0, f0 - op * (f0 - feq)) };
            }
        }
        for &(a, b) in L::PAIRS {
            let (wa, wb) = (T::r(L::W[a]), T::r(L::W[b]));
            let cu_a = dotc::<L, T>(a, u);
            let cu_b = dotc::<L, T>(b, u);
            let feq_a = wa * (drho + r * (three * cu_a + f45 * cu_a * cu_a - f15 * usq));
            let feq_b = wb * (drho + r * (three * cu_b + f45 * cu_b * cu_b - f15 * usq));
            let (ia, ib) = (a * q_stride + i, b * q_stride + i);
            // SAFETY: caller contract.
            let fa = unsafe { planes.get(ia) };
            let fb = unsafe { planes.get(ib) };
            let fp = half * (fa + fb);
            let fm = half * (fa - fb);
            let ep = half * (feq_a + feq_b);
            let em = half * (feq_a - feq_b);
            let rp = op * (fp - ep);
            let rm = om * (fm - em);
            if FORCE {
                let cf_a = dotc::<L, T>(a, fv);
                let cf_b = dotc::<L, T>(b, fv);
                let src_a = wa * (three * (cf_a - uf) + nine * cu_a * cf_a);
                let src_b = wb * (three * (cf_b - uf) + nine * cu_b * cf_b);
                let sp = half * (src_a + src_b);
                let sm = half * (src_a - src_b);
                unsafe {
                    planes.set(ia, fa - rp - rm + cp * sp + cm * sm);
                    planes.set(ib, fb - rp + rm + cp * sp - cm * sm);
                }
            } else {
                unsafe {
                    planes.set(ia, fa - rp - rm);
                    planes.set(ib, fb - rp + rm);
                }
            }
        }
    }
}

/// Force-flavour dispatch for [`collide_span_fused`] (V1 pattern: the
/// `Option` is resolved outside the hot loop, keeping every flavour
/// branch-free).
///
/// # Safety
/// See [`collide_span_fused`]. `field` must cover the row's core cells when
/// present.
#[allow(clippy::too_many_arguments)]
unsafe fn collide_span_dispatch<L: Lattice, T: Real>(
    force_on: bool,
    field: Option<&[[T; 3]]>,
    planes: RawSlice<T>,
    q_stride: usize,
    base: usize,
    x0: usize,
    x1: usize,
    rho: &[T],
    ux: &[T],
    uy: &[T],
    uz: &[T],
    kp: &KParams<T>,
) {
    // SAFETY: forwarded caller contract.
    unsafe {
        match (force_on, field) {
            (true, Some(fr)) => collide_span_fused::<L, T, true, true>(
                planes, q_stride, base, x0, x1, rho, ux, uy, uz, fr, kp,
            ),
            (true, None) => collide_span_fused::<L, T, true, false>(
                planes, q_stride, base, x0, x1, rho, ux, uy, uz, &[], kp,
            ),
            (false, _) => collide_span_fused::<L, T, false, false>(
                planes, q_stride, base, x0, x1, rho, ux, uy, uz, &[], kp,
            ),
        }
    }
}

/// Macroscopic moments over core cells `x0..x1` of one just-streamed row
/// (planes at stride `q_stride`, `base` = core `x = 0`), written to the
/// spare moment buffers at compact offset `c0`. Per-cell arithmetic
/// replicates `kernels::moments_row` exactly.
///
/// # Safety
/// Caller must be the only concurrent writer of the `[c0+x0, c0+x1)` ranges
/// of the moment buffers; plane reads must be in bounds.
#[allow(clippy::too_many_arguments)]
unsafe fn moments_span_fused<L: Lattice, T: Real, const FF: bool>(
    planes: RawSlice<T>,
    q_stride: usize,
    base: usize,
    x0: usize,
    x1: usize,
    c0: usize,
    rho2: RawSlice<T>,
    ux2: RawSlice<T>,
    uy2: RawSlice<T>,
    uz2: RawSlice<T>,
    field: &[[T; 3]],
    kp: &KParams<T>,
) {
    let half = T::r(0.5);
    for x in x0..x1 {
        let i = base + x;
        let mut dr = T::zero();
        let mut m = [T::zero(); 3];
        for q in 0..L::Q {
            // SAFETY: caller contract (row written by this thread only).
            let fq = unsafe { planes.get(q * q_stride + i) };
            dr = dr + fq;
            m[0] = m[0] + T::r(L::C[q][0] as f64) * fq;
            m[1] = m[1] + T::r(L::C[q][1] as f64) * fq;
            if L::D == 3 {
                m[2] = m[2] + T::r(L::C[q][2] as f64) * fq;
            }
        }
        let fv = if FF {
            [
                kp.force[0] + field[x][0],
                kp.force[1] + field[x][1],
                kp.force[2] + field[x][2],
            ]
        } else {
            kp.force
        };
        let r = T::one() + dr;
        let inv = T::one() / r;
        // SAFETY: caller contract (disjoint moment rows).
        unsafe {
            rho2.set(c0 + x, r);
            ux2.set(c0 + x, (m[0] + half * fv[0]) * inv);
            uy2.set(c0 + x, (m[1] + half * fv[1]) * inv);
            if L::D == 3 {
                uz2.set(c0 + x, (m[2] + half * fv[2]) * inv);
            }
        }
    }
}

// ---------------------------------------------------------------------------
// The fused band kernel
// ---------------------------------------------------------------------------

/// Read-only step context shared by all bands (raw out/moment slices carry
/// the row-disjoint write contract).
struct FusedCtx<'a, L: Lattice, T: Real> {
    g: LocalGeom,
    /// Padded plane length (`q` stride of `f` / `out`).
    np: usize,
    /// Padded x extent.
    pnx: usize,
    halo: [bool; 6],
    f: &'a [T],
    out: RawSlice<T>,
    rho_old: &'a [T],
    ux_old: &'a [T],
    uy_old: &'a [T],
    uz_old: &'a [T],
    rho2: RawSlice<T>,
    ux2: RawSlice<T>,
    uy2: RawSlice<T>,
    uz2: RawSlice<T>,
    solid: &'a [bool],
    wall_u: &'a [[T; 3]],
    probe: Option<&'a [bool]>,
    ff: Option<&'a [[T; 3]]>,
    force_on: bool,
    kp: &'a KParams<T>,
    range: CellRange,
    /// Solid runs of the two halo slabs (`[0]` = low, `[1]` = high), one
    /// `Vec` per padded slab row; empty when the face has no halo.
    halo_runs: [Vec<Vec<(u32, u32)>>; 2],
    _l: std::marker::PhantomData<L>,
}

impl<'a, L: Lattice, T: Real> FusedCtx<'a, L, T> {
    /// Slab axis: rows for 2D lattices, z-planes for 3D.
    #[inline(always)]
    fn slab_axis() -> usize {
        L::D - 1
    }

    /// Number of padded rows per slab (1 for 2D, padded-y extent for 3D).
    #[inline(always)]
    fn rows_per_slab(&self) -> usize {
        if L::D == 3 {
            self.g.padded()[1]
        } else {
            1
        }
    }

    /// One slab's plane length (per direction).
    #[inline(always)]
    fn slab_plane(&self) -> usize {
        self.pnx * self.rows_per_slab()
    }

    /// Padded `f`-offset of slab `s` (may be -1 or `n_slabs`: the halo
    /// slabs), i.e. the start of its plane-0 region.
    #[inline(always)]
    fn slab_base_f(&self, s: isize) -> usize {
        ((s + 1) as usize) * self.slab_plane()
    }

    /// Whether the cell at core position `pos` was already collided in
    /// place by the shell pass (its value in `f` is post-collide).
    #[inline(always)]
    fn precollided(&self, pos: [usize; 3]) -> bool {
        let c = self.g.core;
        (pos[0] == 0 && self.halo[0])
            || (pos[0] == c[0] - 1 && self.halo[1])
            || (pos[1] == 0 && self.halo[2])
            || (pos[1] == c[1] - 1 && self.halo[3])
            || (L::D == 3 && ((pos[2] == 0 && self.halo[4]) || (pos[2] == c[2] - 1 && self.halo[5])))
    }
}

/// Ring of just-in-time collided source slabs (3 slots), with per-row solid
/// runs computed at insertion.
struct Ring<T> {
    /// Slot-major slab copies: `data[slot * Q * slab_plane + q * slab_plane + row * pnx + px]`.
    data: Vec<T>,
    /// Slot-major per-row solid runs (padded row coordinates).
    runs: Vec<Vec<(u32, u32)>>,
    tags: [usize; 3],
}

impl<T: Real> Ring<T> {
    fn new(q: usize, slab_plane: usize, rows_per_slab: usize) -> Self {
        Self {
            data: vec![T::zero(); 3 * q * slab_plane],
            runs: vec![Vec::new(); 3 * rows_per_slab],
            tags: [usize::MAX; 3],
        }
    }

    #[inline(always)]
    fn slot_of(&self, s: usize) -> usize {
        self.tags.iter().position(|&t| t == s).expect("resident slab")
    }
}

/// Copy slab `s` into a ring slot, rebuild its solid runs, and collide the
/// core cells that are not already post-collide (shell cells / halo rows are
/// plain copies). No-op when already resident.
fn ensure_slab<L: Lattice, T: Real>(
    ctx: &FusedCtx<'_, L, T>,
    ring: &mut Ring<T>,
    s: usize,
    needed: &[usize; 3],
) {
    if ring.tags.contains(&s) {
        return;
    }
    let slot = (0..3)
        .find(|&k| ring.tags[k] == usize::MAX || !needed.contains(&ring.tags[k]))
        .expect("three ring slots cover at most two other needed slabs");
    ring.tags[slot] = s;
    let sp = ctx.slab_plane();
    let rows = ctx.rows_per_slab();
    let pnx = ctx.pnx;
    let base_f = ctx.slab_base_f(s as isize);
    let region = &mut ring.data[slot * L::Q * sp..(slot + 1) * L::Q * sp];
    for q in 0..L::Q {
        region[q * sp..(q + 1) * sp].copy_from_slice(&ctx.f[q * ctx.np + base_f..][..sp]);
    }
    // Solid runs per padded row.
    for rr in 0..rows {
        solid_runs_row(
            &ctx.solid[base_f + rr * pnx..][..pnx],
            &mut ring.runs[slot * rows + rr],
        );
    }
    // Collide the not-yet-collided core cells (with the previous step's
    // moments — the primaries; the fused pass writes only the spares).
    let c = ctx.g.core;
    let (nx, ny) = (c[0], c[1]);
    let slab_precollided = if L::D == 3 {
        (s == 0 && ctx.halo[4]) || (s == c[2] - 1 && ctx.halo[5])
    } else {
        (s == 0 && ctx.halo[2]) || (s == ny - 1 && ctx.halo[3])
    };
    if slab_precollided {
        return;
    }
    let x_lo = usize::from(ctx.halo[0]);
    let x_hi = nx - usize::from(ctx.halo[1]);
    let region = RawSlice::new(&mut ring.data[slot * L::Q * sp..(slot + 1) * L::Q * sp]);
    let row_iter = if L::D == 3 { 0..ny } else { 0..1 };
    for y in row_iter {
        if L::D == 3 && ((y == 0 && ctx.halo[2]) || (y == ny - 1 && ctx.halo[3])) {
            continue;
        }
        // (y, z) of this row in core coordinates.
        let (cy, cz) = if L::D == 3 { (y, s) } else { (s, 0) };
        let rr = if L::D == 3 { y + 1 } else { 0 };
        let c0 = ctx.g.cidx(0, cy, cz);
        let rho = &ctx.rho_old[c0..c0 + nx];
        let ux = &ctx.ux_old[c0..c0 + nx];
        let uy = &ctx.uy_old[c0..c0 + nx];
        let uz = &ctx.uz_old[c0..c0 + nx];
        let ffrow = ctx.ff.map(|v| &v[c0..c0 + nx]);
        let row_base = rr * pnx + 1; // core x = 0 inside the slab copy
        for_fluid_spans(
            &ring.runs[slot * rows + rr],
            1 + x_lo,
            1 + x_hi,
            |a, b| {
                // SAFETY: this slot region is owned by this band; spans are
                // in bounds of the slab copy.
                unsafe {
                    collide_span_dispatch::<L, T>(
                        ctx.force_on,
                        ffrow,
                        region,
                        sp,
                        row_base,
                        a - 1,
                        b - 1,
                        rho,
                        ux,
                        uy,
                        uz,
                        ctx.kp,
                    )
                }
            },
        );
    }
}

/// One fused collide+stream+moments band over destination slabs
/// `[s0, s1)` ∩ the step range. Returns the momentum-exchange force over
/// probed solid links as **per-destination-row partials** in row order; each
/// partial replays its link contributions in `CpuScalar`'s `stream_row`
/// order (destination x ascending, direction ascending), so the flat fold
/// over all rows reproduces the scalar backend's probe sum bitwise.
fn fused_band<L: Lattice, T: Real>(
    ctx: &FusedCtx<'_, L, T>,
    s0: usize,
    s1: usize,
) -> Vec<[T; 3]> {
    let mut row_partials = Vec::new();
    if s0 >= s1 {
        return row_partials;
    }
    let g = ctx.g;
    let c = g.core;
    let (nx, ny) = (c[0], c[1]);
    let n_slabs = c[FusedCtx::<L, T>::slab_axis()];
    let pnx = ctx.pnx;
    let sp = ctx.slab_plane();
    let rows = ctx.rows_per_slab();
    let six = T::r(6.0);
    let two = T::r(2.0);
    let mut ring = Ring::<T>::new(L::Q, sp, rows);
    // Probed bounce-back links of the current row, gathered in the fused
    // pass's direction-major order and replayed in CpuScalar's cell-major
    // order: (destination x, direction q, ftot).
    let mut links: Vec<(u32, u32, T)> = Vec::new();
    let (xlo, xhi) = (ctx.range.lo[0], ctx.range.hi[0]);
    // Destination rows within a slab (3D: the range's y span; 2D: the slab).
    let (ylo, yhi) = if L::D == 3 {
        (ctx.range.lo[1], ctx.range.hi[1])
    } else {
        (0, 1)
    };
    for s in s0..s1 {
        // Source slabs feeding destination slab `s` (halo slabs excluded:
        // they are read straight from `f`).
        let needed = std::array::from_fn(|k| {
            let ss = s as isize + k as isize - 1;
            if ss < 0 || ss >= n_slabs as isize {
                usize::MAX
            } else {
                ss as usize
            }
        });
        for &ss in &needed {
            if ss != usize::MAX {
                ensure_slab(ctx, &mut ring, ss, &needed);
            }
        }
        let dest_slot = ring.slot_of(s);
        for yy in ylo..yhi {
            // Core coordinates of this destination row.
            let (y, z) = if L::D == 3 { (yy, s) } else { (s, 0) };
            let mut pf_row = [T::zero(); 3];
            links.clear();
            let pb = g.pidx(0, y, z);
            let c0 = g.cidx(0, y, z);
            let rho_row = &ctx.rho_old[c0..c0 + nx];
            let rr_dest = if L::D == 3 { y + 1 } else { 0 };
            for q in 0..L::Q {
                let cq = L::C[q];
                let (cx, cy, cz) = (cq[0] as isize, cq[1] as isize, cq[2] as isize);
                // Resolve the source slab and the row within it.
                let ss = if L::D == 3 {
                    z as isize - cz
                } else {
                    y as isize - cy
                };
                let sy = y as isize - cy;
                if L::D == 3 {
                    if sy < 0 && !ctx.halo[2] {
                        continue;
                    }
                    if sy >= ny as isize && !ctx.halo[3] {
                        continue;
                    }
                }
                let (lo_face, hi_face) = if L::D == 3 { (4, 5) } else { (2, 3) };
                // (src row data, its solid runs, padded f-offset of the row)
                let (src_row, runs_src, src_row_f): (&[T], &[(u32, u32)], usize) = if ss < 0
                    || ss >= n_slabs as isize
                {
                    let hf = if ss < 0 { lo_face } else { hi_face };
                    if !ctx.halo[hf] {
                        continue;
                    }
                    let hidx = usize::from(ss >= 0);
                    let rr = if L::D == 3 { (sy + 1) as usize } else { 0 };
                    let base = ctx.slab_base_f(ss) + rr * pnx;
                    (
                        &ctx.f[q * ctx.np + base..][..pnx],
                        &ctx.halo_runs[hidx][rr],
                        base,
                    )
                } else {
                    let slot = ring.slot_of(ss as usize);
                    let rr = if L::D == 3 { (sy + 1) as usize } else { 0 };
                    let base = ctx.slab_base_f(ss) + rr * pnx;
                    (
                        &ring.data[slot * L::Q * sp + q * sp + rr * pnx..][..pnx],
                        &ring.runs[slot * rows + rr],
                        base,
                    )
                };
                // Destination clamps: the step range, tightened at the x
                // boundary when the crossed face has no halo (open-face
                // unknown slots keep their prior out-buffer contents).
                let dlo = xlo.max(if cx == 1 && !ctx.halo[0] { 1 } else { 0 });
                let dhi = xhi.min(if cx == -1 && !ctx.halo[1] { nx - 1 } else { nx });
                if dlo >= dhi {
                    continue;
                }
                let out_base = q * ctx.np + pb;
                let mut cursor = 0usize;
                let sentinel = (pnx as u32, pnx as u32);
                for &(a, b) in runs_src.iter().chain(std::iter::once(&sentinel)) {
                    let (a, b) = (a as usize, b as usize);
                    // Fluid span [cursor, a): shifted copy ring/halo -> out.
                    let d0 = (cursor as isize - 1 + cx).max(dlo as isize);
                    let d1 = (a as isize - 1 + cx).min(dhi as isize);
                    if d0 < d1 {
                        let (d0, d1) = (d0 as usize, d1 as usize);
                        let s0x = (d0 as isize + 1 - cx) as usize;
                        // SAFETY: this band is the only writer of its
                        // destination rows; ranges are in bounds.
                        unsafe {
                            ctx.out.copy_from(out_base + d0, &src_row[s0x..s0x + (d1 - d0)]);
                        }
                    }
                    // Solid run [a, b): half-way bounce-back into the
                    // destination cells (the reflected populations come from
                    // the destination cell's own collided slab, resident in
                    // the ring).
                    for spx in a..b {
                        let dx = spx as isize - 1 + cx;
                        if dx < dlo as isize || dx >= dhi as isize {
                            continue;
                        }
                        let x = dx as usize;
                        if ctx.solid[pb + x] {
                            continue;
                        }
                        let fout = ring.data
                            [dest_slot * L::Q * sp + L::OPP[q] * sp + rr_dest * pnx + x + 1];
                        let gpi = src_row_f + spx;
                        let wu = ctx.wall_u[gpi];
                        let cu = dotc::<L, T>(q, wu);
                        let fin = fout + six * T::r(L::W[q]) * rho_row[x] * cu;
                        // SAFETY: row-disjoint dispatch.
                        unsafe { ctx.out.set(out_base + x, fin) };
                        if let Some(mask) = ctx.probe {
                            if mask[gpi] {
                                // Momentum given to the wall through this
                                // link, physical populations (V1 probe
                                // convention). Replayed after the direction
                                // loop in CpuScalar's summation order.
                                let ftot = fout + fin + two * T::r(L::W[q]);
                                links.push((x as u32, q as u32, ftot));
                            }
                        }
                    }
                    cursor = b;
                }
            }
            // Replay the probed links in CpuScalar's `stream_row` order:
            // destination cell ascending, direction ascending within the
            // cell (the sort is stable and the fused pass pushes ascending
            // q per destination cell). Bitwise-identical running sum.
            if !links.is_empty() {
                links.sort_by_key(|&(x, _, _)| x);
                for &(_, q, ftot) in &links {
                    let q = q as usize;
                    pf_row[0] = pf_row[0] - T::r(L::C[q][0] as f64) * ftot;
                    pf_row[1] = pf_row[1] - T::r(L::C[q][1] as f64) * ftot;
                    pf_row[2] = pf_row[2] - T::r(L::C[q][2] as f64) * ftot;
                }
            }
            // Fused moments: the destination row is still cache-resident.
            // Written to the spare buffers so in-flight collides of other
            // slabs keep reading the previous step's moments; solid cells
            // are refreshed from the primaries (V1 double-buffer sync).
            let dest_runs = &ring.runs[dest_slot * rows + rr_dest];
            let ffrow = ctx.ff.map(|v| &v[c0..c0 + nx]);
            for_fluid_spans(dest_runs, 1 + xlo, 1 + xhi, |a, b| {
                let (x0, x1) = (a - 1, b - 1);
                // SAFETY: this band owns the destination row (out and the
                // spare moment rows alike).
                unsafe {
                    match ffrow {
                        Some(fr) => moments_span_fused::<L, T, true>(
                            ctx.out, ctx.np, pb, x0, x1, c0, ctx.rho2, ctx.ux2, ctx.uy2,
                            ctx.uz2, fr, ctx.kp,
                        ),
                        None => moments_span_fused::<L, T, false>(
                            ctx.out, ctx.np, pb, x0, x1, c0, ctx.rho2, ctx.ux2, ctx.uy2,
                            ctx.uz2, &[], ctx.kp,
                        ),
                    }
                }
            });
            // Solid cells: keep the double buffers coherent (never computed
            // by the pass; multiphase wall densities live here).
            for &(a, b) in dest_runs.iter() {
                let (a, b) = (a as usize, b as usize);
                for px in a.max(1 + xlo)..b.min(1 + xhi) {
                    let x = px - 1;
                    // SAFETY: row-disjoint dispatch.
                    unsafe {
                        ctx.rho2.set(c0 + x, ctx.rho_old[c0 + x]);
                        ctx.ux2.set(c0 + x, ctx.ux_old[c0 + x]);
                        ctx.uy2.set(c0 + x, ctx.uy_old[c0 + x]);
                        if L::D == 3 {
                            ctx.uz2.set(c0 + x, ctx.uz_old[c0 + x]);
                        }
                    }
                }
            }
            row_partials.push(pf_row);
        }
    }
    row_partials
}

// ---------------------------------------------------------------------------
// Backend implementation
// ---------------------------------------------------------------------------

impl<L: Lattice, T: Real> Backend<L, T> for CpuSimd {
    type Fields = SoaFields<T>;

    fn alloc(&self, sub: &Subdomain) -> SoaFields<T> {
        SoaFields::new(L::Q, sub.geom)
    }

    /// Collide, in place in `f`, exactly the one-cell core boundary shell of
    /// faces with a halo neighbour — the layers the exchange packs. All
    /// other cells stay pre-collide (the fused pass collides them
    /// just-in-time).
    fn collide(&mut self, sub: &Subdomain, fields: &mut SoaFields<T>, p: &StepParams<T>) {
        let g = fields.geom;
        assert_eq!(g.halo, 1, "CpuSimd assumes one-cell halos");
        let halo = sub.halo_flags();
        if !halo.iter().any(|&h| h) {
            return;
        }
        let kp = KParams::new::<L>(p);
        let np = g.n_padded();
        let [nx, ny, nz] = g.core;
        let pnx = g.padded()[0];
        let force_on = p.force[0] != T::zero()
            || p.force[1] != T::zero()
            || p.force[2] != T::zero()
            || fields.force_field.is_some();
        let f = RawSlice::new(&mut fields.f);
        let (rho, ux, uy, uz) = (&fields.rho, &fields.ux, &fields.uy, &fields.uz);
        let solid = &fields.solid;
        let ff = fields.force_field.as_deref();
        let rows = ny * nz;
        let body = |r: usize| {
            let y = r % ny;
            let z = r / ny;
            let full = (z == 0 && halo[4])
                || (z == nz - 1 && halo[5])
                || (y == 0 && halo[2])
                || (y == ny - 1 && halo[3]);
            if !full && !halo[0] && !halo[1] {
                return;
            }
            let pb = g.pidx(0, y, z);
            let c0 = g.cidx(0, y, z);
            let rho_row = &rho[c0..c0 + nx];
            let ux_row = &ux[c0..c0 + nx];
            let uy_row = &uy[c0..c0 + nx];
            let uz_row = &uz[c0..c0 + nx];
            let ffrow = ff.map(|v| &v[c0..c0 + nx]);
            if full {
                let mut runs = Vec::new();
                solid_runs_row(&solid[pb - 1..pb - 1 + pnx], &mut runs);
                for_fluid_spans(&runs, 1, 1 + nx, |a, b| {
                    // SAFETY: rows are dispatched disjointly; spans in bounds.
                    unsafe {
                        collide_span_dispatch::<L, T>(
                            force_on,
                            ffrow,
                            f,
                            np,
                            pb,
                            a - 1,
                            b - 1,
                            rho_row,
                            ux_row,
                            uy_row,
                            uz_row,
                            &kp,
                        )
                    }
                });
            } else {
                let mut cols: [Option<usize>; 2] = [None, None];
                if halo[0] {
                    cols[0] = Some(0);
                }
                if halo[1] && (nx > 1 || !halo[0]) {
                    cols[1] = Some(nx - 1);
                }
                for xc in cols.into_iter().flatten() {
                    if solid[pb + xc] {
                        continue;
                    }
                    // SAFETY: single-cell span, row-disjoint dispatch.
                    unsafe {
                        collide_span_dispatch::<L, T>(
                            force_on,
                            ffrow,
                            f,
                            np,
                            pb,
                            xc,
                            xc + 1,
                            rho_row,
                            ux_row,
                            uy_row,
                            uz_row,
                            &kp,
                        )
                    }
                }
            }
        };
        #[cfg(feature = "parallel")]
        if self.use_parallel(sub) {
            (0..rows).into_par_iter().for_each(body);
            return;
        }
        (0..rows).for_each(body);
    }

    fn stream(
        &mut self,
        sub: &Subdomain,
        fields: &mut SoaFields<T>,
        p: &StepParams<T>,
        range: CellRange,
    ) -> [T; 3] {
        if range.is_empty() {
            return [T::zero(); 3];
        }
        let g = fields.geom;
        assert_eq!(g.halo, 1, "CpuSimd assumes one-cell halos");
        let kp = KParams::new::<L>(p);
        let halo = sub.halo_flags();
        let np = g.n_padded();
        let pnx = g.padded()[0];
        let force_on = p.force[0] != T::zero()
            || p.force[1] != T::zero()
            || p.force[2] != T::zero()
            || fields.force_field.is_some();
        let mut scratch = fields
            .fused
            .take()
            .unwrap_or_else(|| Box::new(FusedScratch::new(g.n_core())));
        let sa = FusedCtx::<L, T>::slab_axis();
        let n_slabs = g.core[sa];
        // Solid runs of the two halo slabs (streamed straight from `f`).
        let rows_per_slab = if L::D == 3 { g.padded()[1] } else { 1 };
        let slab_plane = pnx * rows_per_slab;
        let mk_halo_runs = |present: bool, base: usize| -> Vec<Vec<(u32, u32)>> {
            let mut v = vec![Vec::new(); if present { rows_per_slab } else { 0 }];
            for (rr, runs) in v.iter_mut().enumerate() {
                solid_runs_row(&fields.solid[base + rr * pnx..][..pnx], runs);
            }
            v
        };
        let (lo_face, hi_face) = if L::D == 3 { (4, 5) } else { (2, 3) };
        let halo_runs = [
            mk_halo_runs(halo[lo_face], 0),
            mk_halo_runs(halo[hi_face], (n_slabs + 1) * slab_plane),
        ];
        let ctx = FusedCtx::<L, T> {
            g,
            np,
            pnx,
            halo,
            f: &fields.f,
            out: RawSlice::new(&mut fields.ftmp),
            rho_old: &fields.rho,
            ux_old: &fields.ux,
            uy_old: &fields.uy,
            uz_old: &fields.uz,
            rho2: RawSlice::new(&mut scratch.rho2),
            ux2: RawSlice::new(&mut scratch.ux2),
            uy2: RawSlice::new(&mut scratch.uy2),
            uz2: RawSlice::new(&mut scratch.uz2),
            solid: &fields.solid,
            wall_u: &fields.wall_u,
            probe: fields.probe.as_deref(),
            ff: fields.force_field.as_deref(),
            force_on,
            kp: &kp,
            range,
            halo_runs,
            _l: std::marker::PhantomData,
        };
        let (slo, shi) = (range.lo[sa], range.hi[sa]);
        let n_range = shi - slo;
        #[cfg(feature = "parallel")]
        let nbands = if self.use_parallel(sub) {
            let max = if L::D == 3 {
                (n_range / 4).max(1)
            } else {
                (n_range / 16).max(1)
            };
            rayon::current_num_threads().clamp(1, max)
        } else {
            1
        };
        #[cfg(not(feature = "parallel"))]
        let nbands = 1;
        let band_size = n_range.div_ceil(nbands);
        let body = |band: usize| -> Vec<[T; 3]> {
            let s0 = slo + band * band_size;
            let s1 = (slo + (band + 1) * band_size).min(shi);
            fused_band(&ctx, s0, s1)
        };
        // Flat fold of the per-row partials in global row order — the exact
        // shape of CpuScalar's deterministic probe fold, so the diagnostic
        // is bitwise band-count-independent and backend-independent.
        let fold = |bands: Vec<Vec<[T; 3]>>| -> [T; 3] {
            bands
                .into_iter()
                .flatten()
                .fold([T::zero(); 3], |a, b| [a[0] + b[0], a[1] + b[1], a[2] + b[2]])
        };
        #[cfg(feature = "parallel")]
        let pf = if nbands > 1 {
            fold((0..nbands).into_par_iter().map(body).collect())
        } else {
            fold((0..nbands).map(body).collect())
        };
        #[cfg(not(feature = "parallel"))]
        let pf = fold((0..nbands).map(body).collect());
        // Capture the stale-slot memory for the *next* step's BC pass: the
        // post-collide unknown populations of every open-face cell in range
        // (V1 capture_conv_stale, generalised to every open face).
        capture_stale::<L, T>(sub, &ctx, p, &mut scratch);
        scratch.fresh = true;
        drop(ctx);
        fields.fused = Some(scratch);
        pf
    }

    /// Swap the population ping-pong pair and the moment double buffers.
    fn swap(&mut self, fields: &mut SoaFields<T>) {
        fields.swap_f();
        if let Some(s) = fields.fused.as_deref_mut() {
            std::mem::swap(&mut fields.rho, &mut s.rho2);
            std::mem::swap(&mut fields.ux, &mut s.ux2);
            std::mem::swap(&mut fields.uy, &mut s.uy2);
            std::mem::swap(&mut fields.uz, &mut s.uz2);
        }
    }

    /// Restore the previous step's post-collide populations into the open
    /// faces' unknown slots (the values a `CpuScalar` swap would have left
    /// there), run the shared BC pass, then rotate the stash pair.
    fn apply_open_faces(&mut self, sub: &Subdomain, fields: &mut SoaFields<T>, p: &StepParams<T>) {
        if let Some(mut scratch) = fields.fused.take() {
            let g = fields.geom;
            let np = g.n_padded();
            for face in Face::ALL {
                if face.axis() >= L::D
                    || !sub.touches_global_face(face)
                    || !p.faces[face.index()].is_open()
                {
                    continue;
                }
                let Some(stash) = scratch.stale[0][face.index()].as_ref() else {
                    continue;
                };
                let unknowns = L::unknowns(face);
                for_face_cells(&g, face, |coord, pos| {
                    let pi = g.pidx(pos[0], pos[1], pos[2]);
                    if fields.solid[pi] {
                        return;
                    }
                    for (k, &q) in unknowns.iter().enumerate() {
                        fields.f[q * np + pi] = stash[coord * unknowns.len() + k];
                    }
                });
            }
            apply_open_faces_impl::<L, T>(sub, fields, p);
            let [a, b] = &mut scratch.stale;
            std::mem::swap(a, b);
            fields.fused = Some(scratch);
        } else {
            apply_open_faces_impl::<L, T>(sub, fields, p);
        }
    }

    /// Right after a fused pass, only the open-face boundary cells need a
    /// moment recompute (their populations were just rewritten by the BC
    /// pass; every interior moment was produced in-pass from the identical
    /// post-swap populations). Otherwise (build / init), full recompute.
    fn update_moments(&mut self, sub: &Subdomain, fields: &mut SoaFields<T>, p: &StepParams<T>) {
        let fresh = fields.fused.as_deref().is_some_and(|s| s.fresh);
        if !fresh {
            update_moments_impl::<L, T>(fields, p, self.use_parallel(sub));
            return;
        }
        if let Some(s) = fields.fused.as_deref_mut() {
            s.fresh = false;
        }
        fix_open_face_moments::<L, T>(sub, fields, p);
    }

    fn reduce(
        &self,
        sub: &Subdomain,
        fields: &SoaFields<T>,
        p: &StepParams<T>,
        kind: Reduction,
    ) -> f64 {
        reduce_impl::<L, T>(sub, fields, p, kind)
    }

    fn read_moments(&self, fields: &SoaFields<T>, out: &mut HostMoments<T>) {
        read_moments_impl(fields, out);
    }
}

/// Capture this step's post-collide unknown populations at every open-face
/// cell within `range` into `scratch.stale[1]` (next step's memory term).
/// Shell-precollided cells read their value straight from `f`; every other
/// cell is re-collided from the untouched `f` + old moments — bitwise the
/// value `CpuScalar`'s in-place collide leaves in the ping-pong buffer.
fn capture_stale<L: Lattice, T: Real>(
    sub: &Subdomain,
    ctx: &FusedCtx<'_, L, T>,
    p: &StepParams<T>,
    scratch: &mut FusedScratch<T>,
) {
    let g = ctx.g;
    let np = ctx.np;
    for face in Face::ALL {
        if face.axis() >= L::D
            || !sub.touches_global_face(face)
            || !p.faces[face.index()].is_open()
        {
            continue;
        }
        let unknowns = L::unknowns(face);
        let (t1, t2) = face.tangents();
        let cells = g.core[t1] * g.core[t2];
        // Both stash generations exist from the first fused pass on; the
        // zero-initialised `stale[0]` reproduces CpuScalar's all-zero
        // first-step ping-pong buffer.
        if scratch.stale[0][face.index()].is_none() {
            scratch.stale[0][face.index()] = Some(vec![T::zero(); cells * unknowns.len()]);
        }
        let stash = scratch.stale[1][face.index()]
            .get_or_insert_with(|| vec![T::zero(); cells * unknowns.len()]);
        for_face_cells(&g, face, |coord, pos| {
            if (0..3).any(|a| pos[a] < ctx.range.lo[a] || pos[a] >= ctx.range.hi[a]) {
                return;
            }
            let pi = g.pidx(pos[0], pos[1], pos[2]);
            if ctx.solid[pi] {
                return;
            }
            if ctx.precollided(pos) {
                for (k, &q) in unknowns.iter().enumerate() {
                    stash[coord * unknowns.len() + k] = ctx.f[q * np + pi];
                }
                return;
            }
            let mut cell = [T::zero(); Q_MAX];
            for (q, v) in cell.iter_mut().enumerate().take(L::Q) {
                *v = ctx.f[q * np + pi];
            }
            let cidx = g.cidx(pos[0], pos[1], pos[2]);
            let ffcell = ctx.ff.map(|v| &v[cidx..cidx + 1]);
            // SAFETY: `cell` is thread-local; single-cell span.
            unsafe {
                collide_span_dispatch::<L, T>(
                    ctx.force_on,
                    ffcell,
                    RawSlice::new(&mut cell[..L::Q]),
                    1,
                    0,
                    0,
                    1,
                    &ctx.rho_old[cidx..cidx + 1],
                    &ctx.ux_old[cidx..cidx + 1],
                    &ctx.uy_old[cidx..cidx + 1],
                    &ctx.uz_old[cidx..cidx + 1],
                    ctx.kp,
                );
            }
            for (k, &q) in unknowns.iter().enumerate() {
                stash[coord * unknowns.len() + k] = cell[q];
            }
        });
    }
}

/// Recompute the moments of every open-face cell from the post-BC
/// populations (V1 `fix_boundary_moments`, restricted to the faces the BC
/// pass touches; idempotent per cell). Per-cell arithmetic replicates
/// `kernels::moments_row`.
fn fix_open_face_moments<L: Lattice, T: Real>(
    sub: &Subdomain,
    fields: &mut SoaFields<T>,
    p: &StepParams<T>,
) {
    let kp = KParams::new::<L>(p);
    let g = fields.geom;
    let np = g.n_padded();
    let half = T::r(0.5);
    for face in Face::ALL {
        if face.axis() >= L::D
            || !sub.touches_global_face(face)
            || !p.faces[face.index()].is_open()
        {
            continue;
        }
        let ff = fields.force_field.as_deref();
        let f = &fields.f;
        let solid = &fields.solid;
        let (rho, ux, uy, uz) = (
            RawSlice::new(&mut fields.rho),
            RawSlice::new(&mut fields.ux),
            RawSlice::new(&mut fields.uy),
            RawSlice::new(&mut fields.uz),
        );
        for_face_cells(&g, face, |_, pos| {
            let pi = g.pidx(pos[0], pos[1], pos[2]);
            if solid[pi] {
                return;
            }
            let c = g.cidx(pos[0], pos[1], pos[2]);
            let mut dr = T::zero();
            let mut m = [T::zero(); 3];
            for q in 0..L::Q {
                let fq = f[q * np + pi];
                dr = dr + fq;
                m[0] = m[0] + T::r(L::C[q][0] as f64) * fq;
                m[1] = m[1] + T::r(L::C[q][1] as f64) * fq;
                if L::D == 3 {
                    m[2] = m[2] + T::r(L::C[q][2] as f64) * fq;
                }
            }
            let fv = match ff {
                Some(field) => [
                    kp.force[0] + field[c][0],
                    kp.force[1] + field[c][1],
                    kp.force[2] + field[c][2],
                ],
                None => kp.force,
            };
            let r = T::one() + dr;
            let inv = T::one() / r;
            // SAFETY: sequential pass, exclusive access.
            unsafe {
                rho.set(c, r);
                ux.set(c, (m[0] + half * fv[0]) * inv);
                uy.set(c, (m[1] + half * fv[1]) * inv);
                if L::D == 3 {
                    uz.set(c, (m[2] + half * fv[2]) * inv);
                }
            }
        });
    }
}
