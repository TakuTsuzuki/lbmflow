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
//! Streaming (copies), bounce-back, moments, the boundary fix, the BC pass
//! and the diagnostics replicate the corresponding `kernels.rs` expressions
//! operand-for-operand (the lattice constants are compile-time promoted,
//! which folds the ±1/0 multiplies without changing IEEE results), and the
//! probed-force diagnostic replays its link contributions in `stream_row`'s
//! cell-major order with `CpuScalar`'s flat per-row fold.
//!
//! The collision uses **V1 `collide_span`'s pair-shared equilibrium form**
//! (`base`/`r3`/`r45` hoisted per cell — on D2Q9 this is V1's fused kernel
//! operand-for-operand), which is an exact algebraic regrouping of
//! `kernels::collide_row`'s TRT update but rounds differently at the last
//! ulp. The trade is deliberate and measured: the literal `collide_row` DAG
//! costs ~16% single-thread throughput (evaluated 2026-07: 182 vs 219 MLUPS,
//! f32 512² 1T) for zero physical benefit, while the resulting
//! `CpuScalar`↔`CpuSimd` drift is the same reassociation noise the
//! `v1_match` suite absorbs: ≤ ~1e-13 (f64) / ~1e-7 (f32) over hundreds of
//! steps, asserted at 1e-11 / 1e-6 by `tests/backend_simd_equiv.rs`.

use crate::backend::{
    apply_open_faces_impl, read_moments_impl, reduce_impl, update_moments_impl, Backend, CellRange,
    CpuScalar, HostMoments, PARALLEL_MIN_CELLS,
};
use crate::fields::{FusedScratch, LocalGeom, SoaFields};
use crate::halo::HaloExchange;
use crate::kernels::{central_basis, central_phi, for_face_cells, solve_moment_system, RawSlice};
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
fn for_fluid_spans(runs: &[(u32, u32)], w0: usize, w1: usize, mut body: impl FnMut(usize, usize)) {
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

/// Cells per kernel block for the high-`Q` (3D) span kernels: they stage
/// shared per-cell pieces in `[T; BLOCK]` stack rows and sweep the block
/// once per direction pair / plane, so every inner loop is a unit-stride
/// sweep with loop-invariant lattice constants — vectorization no longer
/// depends on `L::Q`. The flat per-cell form (used for D2Q9, where LLVM
/// fully unrolls the 4 pairs and vectorizes the cell loop) left D3Q19
/// essentially scalar: measured 18 vs 285 vector/scalar instructions and
/// 22 vs 43 MLUPS at 128³ f32 1T. Conversely, block staging costs D2Q9
/// ~1.2x (extra scratch traffic against only 4 pair sweeps), so each
/// lattice gets the form that measures faster.
const BLOCK: usize = 64;

/// Whether lattice `L` uses the blocked span kernels (see [`BLOCK`]).
#[inline(always)]
fn use_blocked<L: Lattice>() -> bool {
    L::Q > 9
}

/// TRT collision with Guo forcing over core cells `x0..x1` of one row stored
/// as `L::Q` planes at stride `q_stride` inside `planes` (`base` = index of
/// the row's core `x = 0` cell in plane 0).
///
/// Per-cell arithmetic is **V1 `collide_span` operand-for-operand** (the
/// pair decomposition works directly on the shared `base`/`r3`/`r45`
/// equilibrium pieces — on D2Q9, `dotc`'s folded ±1/0 constants reproduce
/// V1's hand-written `vx`/`vy`/`vx+vy`/`vy-vx` bitwise), generalised to any
/// lattice via `L::PAIRS` and staged through per-block scratch rows (which
/// does not change any rounding). It is an exact algebraic regrouping of
/// `kernels::collide_row`, differing from it by last-ulp rounding only (see
/// the module docs for the measured cost of the literal DAG). The caller
/// decomposes rows into solid-free spans, so the loops are branch-free.
///
/// # Safety
/// The caller must be the only concurrent accessor of the addressed cells,
/// and `q * q_stride + base + x` must be in bounds for all `q < L::Q`,
/// `x0 <= x < x1`.
#[allow(clippy::too_many_arguments)]
unsafe fn collide_span_fused<L: Lattice, T: Real, const FORCE: bool, const FF: bool>(
    src: PlaneView<T>,
    dst: PlaneView<T>,
    x0: usize,
    x1: usize,
    rho: &[T],
    ux: &[T],
    uy: &[T],
    uz: &[T],
    field: &[[T; 3]],
    omega: Option<&[T]>,
    kp: &KParams<T>,
) {
    if x0 >= x1 {
        return;
    }
    // SAFETY: forwarded caller contract.
    unsafe {
        if use_blocked::<L>() {
            // The blocked form is in-place only (its callers never fuse the
            // ring copy — see `ensure_slab`); a single view keeps the pair
            // sweeps free of aliasing checks, which is what lets them
            // vectorize (a src/dst pair measured ~30% slower end to end).
            debug_assert!(std::ptr::eq(src.planes.as_ptr(), dst.planes.as_ptr()));
            collide_span_blocked::<L, T, FORCE, FF>(dst, x0, x1, rho, ux, uy, uz, field, omega, kp);
        } else {
            collide_span_flat::<L, T, FORCE, FF>(
                src, dst, x0, x1, rho, ux, uy, uz, field, omega, kp,
            );
        }
    }
}

/// A `(planes, q_stride, base)` triple addressing one row of a q-major
/// plane set: direction `q`'s cell `x` lives at `q * stride + base + x`.
/// The collide span kernels read populations from a source view and write
/// results to a destination view — distinct views fuse the ring copy into
/// the collision (source `f`, destination ring); identical views collide in
/// place (each cell is fully read before it is written).
#[derive(Clone, Copy)]
struct PlaneView<T> {
    planes: RawSlice<T>,
    stride: usize,
    base: usize,
}

impl<T: Real> PlaneView<T> {
    #[inline(always)]
    fn idx(&self, q: usize, x: usize) -> usize {
        q * self.stride + self.base + x
    }
}

/// Flat per-cell form of [`collide_span_fused`] (D2Q9: LLVM unrolls the
/// pair loop into the cell body and vectorizes the cell loop — V1
/// `collide_span`'s exact shape and arithmetic).
///
/// # Safety
/// See [`collide_span_fused`].
#[allow(clippy::too_many_arguments)]
unsafe fn collide_span_flat<L: Lattice, T: Real, const FORCE: bool, const FF: bool>(
    src: PlaneView<T>,
    dst: PlaneView<T>,
    x0: usize,
    x1: usize,
    rho: &[T],
    ux: &[T],
    uy: &[T],
    uz: &[T],
    field: &[[T; 3]],
    omega: Option<&[T]>,
    kp: &KParams<T>,
) {
    let three = T::r(3.0);
    let f45 = T::r(4.5);
    let f15 = T::r(1.5);
    let nine = T::r(9.0);
    let half = T::r(0.5);
    let one = T::one();
    let (op0, om, cp0, cm) = (kp.omega_p, kp.omega_m, kp.cp, kp.cm);
    for x in x0..x1 {
        let op = omega.map_or(op0, |v| v[x]);
        let cp = if omega.is_some() {
            one - op / T::r(2.0)
        } else {
            cp0
        };
        let r = rho[x];
        let u = [ux[x], uy[x], uz[x]];
        let fv = if FORCE {
            kp.force_at(if FF { Some(field) } else { None }, x, r)
        } else {
            [T::zero(); 3]
        };
        let mut usq = u[0] * u[0];
        for d in 1..L::D {
            usq = usq + u[d] * u[d];
        }
        let drho = r - one;
        // Shared equilibrium pieces (V1: base / r3 / r45 / uf3).
        let eq_base = drho - f15 * r * usq;
        let r3 = three * r;
        let r45 = f45 * r;
        let uf3 = if FORCE {
            let mut uf = u[0] * fv[0];
            for d in 1..L::D {
                uf = uf + u[d] * fv[d];
            }
            three * uf
        } else {
            T::zero()
        };
        // Rest population: feq0 = w0 * base, src0 = -w0 * uf3.
        {
            let w0 = T::r(L::W[L::REST]);
            // SAFETY: caller contract (disjoint cells, in bounds).
            let f0 = unsafe { src.planes.get(src.idx(L::REST, x)) };
            let v = if FORCE {
                f0 - op * (f0 - w0 * eq_base) + cp * (-w0 * uf3)
            } else {
                f0 - op * (f0 - w0 * eq_base)
            };
            unsafe { dst.planes.set(dst.idx(L::REST, x), v) };
        }
        for &(a, b) in L::PAIRS {
            let wa = T::r(L::W[a]);
            let cu = dotc::<L, T>(a, u);
            let ep = wa * (eq_base + r45 * cu * cu);
            let em = wa * (r3 * cu);
            // SAFETY: caller contract.
            let fa = unsafe { src.planes.get(src.idx(a, x)) };
            let fb = unsafe { src.planes.get(src.idx(b, x)) };
            let fp = half * (fa + fb);
            let fm = half * (fa - fb);
            let rp = op * (fp - ep);
            let rm = om * (fm - em);
            if FORCE {
                let cf = dotc::<L, T>(a, fv);
                let sp = wa * (nine * cu * cf - uf3);
                let sm = wa * (three * cf);
                unsafe {
                    dst.planes
                        .set(dst.idx(a, x), fa - rp - rm + cp * sp + cm * sm);
                    dst.planes
                        .set(dst.idx(b, x), fb - rp + rm + cp * sp - cm * sm);
                }
            } else {
                unsafe {
                    dst.planes.set(dst.idx(a, x), fa - rp - rm);
                    dst.planes.set(dst.idx(b, x), fb - rp + rm);
                }
            }
        }
    }
}

/// Blocked form of [`collide_span_fused`] (D3Q19; see [`BLOCK`]).
///
/// # Safety
/// See [`collide_span_fused`].
#[allow(clippy::too_many_arguments)]
unsafe fn collide_span_blocked<L: Lattice, T: Real, const FORCE: bool, const FF: bool>(
    planes: PlaneView<T>,
    x0: usize,
    x1: usize,
    rho: &[T],
    ux: &[T],
    uy: &[T],
    uz: &[T],
    field: &[[T; 3]],
    omega: Option<&[T]>,
    kp: &KParams<T>,
) {
    let three = T::r(3.0);
    let f45 = T::r(4.5);
    let f15 = T::r(1.5);
    let nine = T::r(9.0);
    let half = T::r(0.5);
    let one = T::one();
    let (op0, om, cp0, cm) = (kp.omega_p, kp.omega_m, kp.cp, kp.cm);
    let mut xb = x0;
    while xb < x1 {
        let blen = (x1 - xb).min(BLOCK);
        // Pass A: shared per-cell equilibrium/forcing pieces
        // (V1: base / r3 / r45 / uf3).
        let mut eqb = [T::zero(); BLOCK];
        let mut r3v = [T::zero(); BLOCK];
        let mut r45v = [T::zero(); BLOCK];
        let mut uf3v = [T::zero(); BLOCK];
        let mut fvr = [[T::zero(); 3]; BLOCK];
        for j in 0..blen {
            let x = xb + j;
            let r = rho[x];
            let u = [ux[x], uy[x], uz[x]];
            let mut usq = u[0] * u[0];
            for d in 1..L::D {
                usq = usq + u[d] * u[d];
            }
            let drho = r - one;
            eqb[j] = drho - f15 * r * usq;
            r3v[j] = three * r;
            r45v[j] = f45 * r;
            if FORCE {
                let fv = kp.force_at(if FF { Some(field) } else { None }, x, r);
                fvr[j] = fv;
                let mut uf = u[0] * fv[0];
                for d in 1..L::D {
                    uf = uf + u[d] * fv[d];
                }
                uf3v[j] = three * uf;
            }
        }
        // Rest population: feq0 = w0 * base, src0 = -w0 * uf3.
        {
            let w0 = T::r(L::W[L::REST]);
            let i0 = planes.idx(L::REST, xb);
            for j in 0..blen {
                // SAFETY: caller contract (disjoint cells, in bounds).
                let x = xb + j;
                let op = omega.map_or(op0, |v| v[x]);
                let cp = if omega.is_some() {
                    one - op / T::r(2.0)
                } else {
                    cp0
                };
                let f0 = unsafe { planes.planes.get(i0 + j) };
                let v = if FORCE {
                    f0 - op * (f0 - w0 * eqb[j]) + cp * (-w0 * uf3v[j])
                } else {
                    f0 - op * (f0 - w0 * eqb[j])
                };
                unsafe { planes.planes.set(i0 + j, v) };
            }
        }
        // One sweep per TRT pair; c_a is loop-invariant, so `dotc` becomes
        // broadcast multiply-adds.
        for &(a, b) in L::PAIRS {
            let wa = T::r(L::W[a]);
            let (ia0, ib0) = (planes.idx(a, xb), planes.idx(b, xb));
            for j in 0..blen {
                let x = xb + j;
                let cu = dotc::<L, T>(a, [ux[x], uy[x], uz[x]]);
                let ep = wa * (eqb[j] + r45v[j] * cu * cu);
                let em = wa * (r3v[j] * cu);
                // SAFETY: caller contract.
                let fa = unsafe { planes.planes.get(ia0 + j) };
                let fb = unsafe { planes.planes.get(ib0 + j) };
                let fp = half * (fa + fb);
                let fm = half * (fa - fb);
                let op = omega.map_or(op0, |v| v[x]);
                let cp = if omega.is_some() {
                    one - op / T::r(2.0)
                } else {
                    cp0
                };
                let rp = op * (fp - ep);
                let rm = om * (fm - em);
                if FORCE {
                    let cf = dotc::<L, T>(a, fvr[j]);
                    let sp = wa * (nine * cu * cf - uf3v[j]);
                    let sm = wa * (three * cf);
                    unsafe {
                        planes.planes.set(ia0 + j, fa - rp - rm + cp * sp + cm * sm);
                        planes.planes.set(ib0 + j, fb - rp + rm + cp * sp - cm * sm);
                    }
                } else {
                    unsafe {
                        planes.planes.set(ia0 + j, fa - rp - rm);
                        planes.planes.set(ib0 + j, fb - rp + rm);
                    }
                }
            }
        }
        xb += blen;
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
    omega: Option<&[T]>,
    src: PlaneView<T>,
    dst: PlaneView<T>,
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
        if kp.cumulant {
            collide_span_central_moment::<L, T>(
                field, omega, src, dst, x0, x1, rho, ux, uy, uz, kp,
            );
            return;
        }
        match (force_on, field) {
            (true, Some(fr)) => collide_span_fused::<L, T, true, true>(
                src, dst, x0, x1, rho, ux, uy, uz, fr, omega, kp,
            ),
            (true, None) => collide_span_fused::<L, T, true, false>(
                src,
                dst,
                x0,
                x1,
                rho,
                ux,
                uy,
                uz,
                &[],
                omega,
                kp,
            ),
            (false, _) => collide_span_fused::<L, T, false, false>(
                src,
                dst,
                x0,
                x1,
                rho,
                ux,
                uy,
                uz,
                &[],
                omega,
                kp,
            ),
        }
    }
}

#[allow(clippy::too_many_arguments)]
unsafe fn collide_span_central_moment<L: Lattice, T: Real>(
    field: Option<&[[T; 3]]>,
    omega: Option<&[T]>,
    src: PlaneView<T>,
    dst: PlaneView<T>,
    x0: usize,
    x1: usize,
    rho: &[T],
    ux: &[T],
    uy: &[T],
    uz: &[T],
    kp: &KParams<T>,
) {
    let basis = central_basis::<L>();
    for x in x0..x1 {
        let r_t = rho[x];
        let r = r_t.as_f64();
        let u = [ux[x].as_f64(), uy[x].as_f64(), uz[x].as_f64()];
        let fv_t = kp.force_at(field, x, r_t);
        let fv = [fv_t[0].as_f64(), fv_t[1].as_f64(), fv_t[2].as_f64()];
        let force_active = fv[0] != 0.0 || fv[1] != 0.0 || fv[2] != 0.0;
        let mut phys = [0.0f64; Q_MAX];
        let mut src_pop = [0.0f64; Q_MAX];
        let feq_rest = u.iter().take(L::D).all(|&v| v == 0.0) && !force_active;
        let mut exact_rest_equilibrium = feq_rest;
        for q in 0..L::Q {
            let fq = unsafe { src.planes.get(src.idx(q, x)) };
            phys[q] = fq.as_f64() + L::W[q];
            if exact_rest_equilibrium {
                let expected = T::r(L::W[q] * (r - 1.0));
                exact_rest_equilibrium &= fq.as_f64().to_bits() == expected.as_f64().to_bits();
            }
        }
        if exact_rest_equilibrium {
            if !std::ptr::eq(src.planes.as_ptr(), dst.planes.as_ptr())
                || src.stride != dst.stride
                || src.base != dst.base
            {
                for q in 0..L::Q {
                    let fq = unsafe { src.planes.get(src.idx(q, x)) };
                    unsafe { dst.planes.set(dst.idx(q, x), fq) };
                }
            }
            continue;
        }
        if force_active {
            let mut uf = u[0] * fv[0];
            for d in 1..L::D {
                uf += u[d] * fv[d];
            }
            for q in 0..L::Q {
                let mut cu = L::C[q][0] as f64 * u[0];
                let mut cf = L::C[q][0] as f64 * fv[0];
                for d in 1..L::D {
                    cu += L::C[q][d] as f64 * u[d];
                    cf += L::C[q][d] as f64 * fv[d];
                }
                src_pop[q] = L::W[q] * (3.0 * (cf - uf) + 9.0 * cu * cf);
            }
        }
        let mut usq = u[0] * u[0];
        for d in 1..L::D {
            usq += u[d] * u[d];
        }
        let mut feq_phys = [0.0f64; Q_MAX];
        for q in 0..L::Q {
            let mut cu = L::C[q][0] as f64 * u[0];
            for d in 1..L::D {
                cu += L::C[q][d] as f64 * u[d];
            }
            feq_phys[q] = L::W[q] * r * (1.0 + 3.0 * cu + 4.5 * cu * cu - 1.5 * usq);
        }
        let mut mom = [0.0f64; Q_MAX];
        let mut src_mom = [0.0f64; Q_MAX];
        let mut eq = [0.0f64; Q_MAX];
        for m in 0..L::Q {
            for q in 0..L::Q {
                let phi = central_phi::<L>(q, basis[m], u);
                mom[m] += phi * phys[q];
                eq[m] += phi * feq_phys[q];
                src_mom[m] += phi * src_pop[q];
            }
        }
        let os_base = omega.map_or(kp.omega_shear.as_f64(), |v| v[x].as_f64());
        let usq_full = u[0] * u[0] + u[1] * u[1] + u[2] * u[2];
        let d3q19_lattice_viscosity_offset = if L::D == 3 && L::Q == 19 { 0.0025 } else { 0.0 };
        let os = (os_base * (1.0 + d3q19_lattice_viscosity_offset - 0.16 * usq_full)).min(2.0);
        let mut post = [0.0f64; Q_MAX];
        let mut diag = [usize::MAX; 3];
        for m in 0..L::Q {
            let e = basis[m];
            let order = e[0] as usize + e[1] as usize + e[2] as usize;
            if order == 2 {
                for a in 0..L::D {
                    let mut de = [0u8; 3];
                    de[a] = 2;
                    if e == de {
                        diag[a] = m;
                    }
                }
            }
            let rate = match order {
                0 | 1 => 0.0,
                2 => os,
                _ => 1.0,
            };
            post[m] = mom[m] - rate * (mom[m] - eq[m]) + (1.0 - 0.5 * rate) * src_mom[m];
        }
        let inv_d = 1.0 / L::D as f64;
        let mut trace_neq = 0.0;
        let mut trace_src = 0.0;
        for &idx in diag.iter().take(L::D) {
            debug_assert_ne!(idx, usize::MAX);
            trace_neq += mom[idx] - eq[idx];
            trace_src += src_mom[idx];
        }
        let bulk_neq = trace_neq * inv_d;
        let bulk_src = trace_src * inv_d;
        for &idx in diag.iter().take(L::D) {
            let dev_neq = mom[idx] - eq[idx] - bulk_neq;
            let dev_src = src_mom[idx] - bulk_src;
            post[idx] =
                eq[idx] + (1.0 - os) * dev_neq + 0.5 * bulk_src + (1.0 - 0.5 * os) * dev_src;
        }
        let out_phys = solve_moment_system::<L>(&basis, u, &post);
        for q in 0..L::Q {
            unsafe { dst.planes.set(dst.idx(q, x), T::r(out_phys[q] - L::W[q])) };
        }
    }
}

/// Macroscopic moments over core cells `x0..x1` of one just-streamed row
/// (planes at stride `q_stride`, `base` = core `x = 0`), written to the
/// spare moment buffers at compact offset `c0`. Per-cell arithmetic
/// replicates `kernels::moments_row` exactly (block staging accumulates the
/// same q-ascending sums per cell), one vectorizable sweep per direction.
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
    // SAFETY: forwarded caller contract.
    unsafe {
        if use_blocked::<L>() {
            moments_span_blocked::<L, T, FF>(
                planes, q_stride, base, x0, x1, c0, rho2, ux2, uy2, uz2, field, kp,
            );
        } else {
            moments_span_flat::<L, T, FF>(
                planes, q_stride, base, x0, x1, c0, rho2, ux2, uy2, uz2, field, kp,
            );
        }
    }
}

/// Flat per-cell form of [`moments_span_fused`] (D2Q9).
///
/// # Safety
/// See [`moments_span_fused`].
#[allow(clippy::too_many_arguments)]
unsafe fn moments_span_flat<L: Lattice, T: Real, const FF: bool>(
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
    let one = T::one();
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
        let r = one + dr;
        let fv = kp.force_at(if FF { Some(field) } else { None }, x, r);
        let inv = one / r;
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

/// Blocked form of [`moments_span_fused`] (D3Q19; see [`BLOCK`]).
///
/// # Safety
/// See [`moments_span_fused`].
#[allow(clippy::too_many_arguments)]
unsafe fn moments_span_blocked<L: Lattice, T: Real, const FF: bool>(
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
    let one = T::one();
    let mut xb = x0;
    while xb < x1 {
        let blen = (x1 - xb).min(BLOCK);
        let mut dr = [T::zero(); BLOCK];
        let mut m0 = [T::zero(); BLOCK];
        let mut m1 = [T::zero(); BLOCK];
        let mut m2 = [T::zero(); BLOCK];
        for q in 0..L::Q {
            let cq = [
                T::r(L::C[q][0] as f64),
                T::r(L::C[q][1] as f64),
                T::r(L::C[q][2] as f64),
            ];
            let iq = q * q_stride + base + xb;
            for j in 0..blen {
                // SAFETY: caller contract (row written by this thread only).
                let fq = unsafe { planes.get(iq + j) };
                dr[j] = dr[j] + fq;
                m0[j] = m0[j] + cq[0] * fq;
                m1[j] = m1[j] + cq[1] * fq;
                if L::D == 3 {
                    m2[j] = m2[j] + cq[2] * fq;
                }
            }
        }
        for j in 0..blen {
            let x = xb + j;
            let r = one + dr[j];
            let fv = kp.force_at(if FF { Some(field) } else { None }, x, r);
            let inv = one / r;
            // SAFETY: caller contract (disjoint moment rows).
            unsafe {
                rho2.set(c0 + x, r);
                ux2.set(c0 + x, (m0[j] + half * fv[0]) * inv);
                uy2.set(c0 + x, (m1[j] + half * fv[1]) * inv);
                if L::D == 3 {
                    uz2.set(c0 + x, (m2[j] + half * fv[2]) * inv);
                }
            }
        }
        xb += blen;
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
    /// Padded y extent (1 for 2D — matches `LocalGeom` padding).
    pny: usize,
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
    omega: Option<&'a [T]>,
    force_on: bool,
    kp: &'a KParams<T>,
    range: CellRange,
    _l: std::marker::PhantomData<L>,
}

impl<'a, L: Lattice, T: Real> FusedCtx<'a, L, T> {
    /// Number of slabs along the sweep axis (2D: y rows; 3D: z planes).
    #[inline(always)]
    fn n_slabs(&self) -> usize {
        self.g.core[L::D - 1]
    }

    /// Padded `f`-offset of row `sy` of slab `slab` (either index may be -1
    /// or the extent: the halo layers). For 2D lattices the slab *is* the
    /// row and `sy` is ignored.
    #[inline(always)]
    fn row_base_f(&self, slab: isize, sy: isize) -> usize {
        if L::D == 3 {
            (((slab + 1) as usize) * self.pny + (sy + 1) as usize) * self.pnx
        } else {
            ((slab + 1) as usize) * self.pnx
        }
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
            || (L::D == 3
                && ((pos[2] == 0 && self.halo[4]) || (pos[2] == c[2] - 1 && self.halo[5])))
    }
}

/// The window of slab rows one ring pass works on, in padded row indices.
///
/// 2D slabs are single rows (`row0 = 0`, `rows = 1`). A 3D z-band processes
/// its destination planes in **y-strips**: destination rows `[y0, y1)` need
/// source rows `[y0-1, y1]` of the three z-slabs around each destination
/// plane — the strip `[y0, y1+2)` in padded indices, swept over the band's
/// z range so each (strip × slab) is collided once per strip pass. Ringing
/// strips instead of whole planes keeps the ring cache-resident for any
/// grid (a 128³ f64 D3Q19 full-plane ring is 7.7 MB per band — 92 MB across
/// 12 bands, past even the SLC; a 32-row strip ring is ~2 MB) at the cost
/// of the two strip-edge rows being collided by both adjacent strip passes
/// of the *same* band. Bands still partition whole z-planes, so no two
/// threads ever write neighbouring rows of `out` (measured: y-partitioned
/// bands lose ~20% at 12 threads to exactly that boundary sharing).
#[derive(Clone, Copy, Debug, Default)]
pub(crate) struct Window {
    /// First padded row of the strip within a slab.
    row0: usize,
    /// Padded rows per strip.
    rows: usize,
}

/// Destination rows per 3D y-strip pass (see [`Window`]). `usize::MAX`
/// disables striping: one full-plane window per band. Measured on M5 Max
/// (128³, 12 threads, same-window A/B): 32-row strips lose ~20% against
/// full-plane rings in both precisions — the system-level cache absorbs the
/// per-band plane rings (3.9/7.7 MB × 12), so the strips' extra pass
/// overhead and +6% edge-row recollides never pay off. The machinery stays
/// for cache-poorer targets; re-tune there before enabling.
const STRIP_ROWS: usize = usize::MAX;

/// Ring of just-in-time collided source slab strips (3 slots), with
/// per-row solid runs computed at insertion. Owned per (part, band) by
/// [`FusedScratch`] and reused across steps and strip passes — the buffers
/// are only (re)zeroed when they grow, so the steady state pays no
/// allocation or memset traffic (a 128³ f64 D3Q19 ring is 7.7 MB per band;
/// zeroing 12 of them every step measurably costs ~10%).
#[derive(Clone, Debug, Default)]
pub(crate) struct Ring<T> {
    /// Slot-major strip copies:
    /// `data[(slot * Q + q) * cap_len + row * pnx + px]`.
    data: Vec<T>,
    /// Slot-major per-row solid runs (padded row coordinates), stride
    /// `cap_rows` per slot.
    runs: Vec<Vec<(u32, u32)>>,
    tags: [usize; 3],
    /// Capacity in padded rows per strip (allocation stride).
    cap_rows: usize,
    /// One strip slot's plane capacity (`cap_rows * pnx`).
    cap_len: usize,
    /// Active window (`win.rows <= cap_rows`).
    win: Window,
}

impl<T: Real> Ring<T> {
    /// Grow (never shrink) the buffers for the given geometry; contents are
    /// don't-care (every slot is fully rewritten before any read).
    fn prepare(&mut self, q: usize, pnx: usize, cap_rows: usize) {
        let len = 3 * q * cap_rows * pnx;
        if self.data.len() < len {
            self.data.resize(len, T::zero());
        }
        if self.runs.len() < 3 * cap_rows {
            self.runs.resize(3 * cap_rows, Vec::new());
        }
        self.cap_rows = cap_rows;
        self.cap_len = cap_rows * pnx;
        self.tags = [usize::MAX; 3];
    }

    /// Begin a strip pass: invalidate all slots and set the active window.
    fn set_window(&mut self, win: Window) {
        debug_assert!(win.rows <= self.cap_rows);
        self.tags = [usize::MAX; 3];
        self.win = win;
    }

    #[inline(always)]
    fn slot_of(&self, s: usize) -> usize {
        self.tags
            .iter()
            .position(|&t| t == s)
            .expect("resident slab")
    }
}

/// Bring the window of slab `s` into a ring slot: rebuild its solid runs,
/// collide the not-yet-collided core cells, and copy everything else
/// verbatim (halo rows and columns, shell-precollided cells, solid cells).
/// No-op when already resident.
///
/// For the flat-form lattices (D2Q9) the copy is **fused into the
/// collision** — the kernel reads `f` and writes the ring, saving a full
/// read+write pass per cell (measured +15% at 512² f32 1T, ahead of V1's
/// copy-then-collide). The blocked lattices (D3Q19) keep copy-then-collide:
/// their per-pair sweeps would otherwise walk 19 far-apart source streams
/// per block and lose more to prefetch thrash than the copy costs
/// (measured −10% at 128³ 12T).
fn ensure_slab<L: Lattice, T: Real>(
    ctx: &FusedCtx<'_, L, T>,
    ring: &mut Ring<T>,
    s: usize,
    needed: &[usize; 3],
) {
    if ring.tags.contains(&s) {
        return;
    }
    let win = ring.win;
    let slot = (0..3)
        .find(|&k| ring.tags[k] == usize::MAX || !needed.contains(&ring.tags[k]))
        .expect("three ring slots cover at most two other needed slabs");
    ring.tags[slot] = s;
    let (cap_len, cap_rows) = (ring.cap_len, ring.cap_rows);
    let pnx = ctx.pnx;
    let base_f = ctx.row_base_f(s as isize, win.row0 as isize - 1);
    // Solid runs per padded strip row (also the copy/collide span plan).
    for rr in 0..win.rows {
        solid_runs_row(
            &ctx.solid[base_f + rr * pnx..][..pnx],
            &mut ring.runs[slot * cap_rows + rr],
        );
    }
    let c = ctx.g.core;
    let (nx, ny) = (c[0], c[1]);
    let slab_precollided = if L::D == 3 {
        (s == 0 && ctx.halo[4]) || (s == c[2] - 1 && ctx.halo[5])
    } else {
        (s == 0 && ctx.halo[2]) || (s == ny - 1 && ctx.halo[3])
    };
    let x_lo = usize::from(ctx.halo[0]);
    let x_hi = nx - usize::from(ctx.halo[1]);
    let fuse_copy = !use_blocked::<L>();
    let dst_raw = RawSlice::new(&mut ring.data);
    if !fuse_copy {
        // Copy-then-collide (blocked lattices): one bulk copy per direction
        // plane brings the whole strip in, then the ring is collided in
        // place. Fusing the copy here would make every per-pair sweep walk
        // `L::Q` far-apart source streams per block — measured ~35% slower
        // at 128³ (prefetch thrash beats the saved pass).
        let sl = win.rows * pnx;
        for q in 0..L::Q {
            // SAFETY: this band owns the slot; ranges in bounds.
            unsafe {
                dst_raw.copy_from(
                    slot * L::Q * cap_len + q * cap_len,
                    &ctx.f[q * ctx.np + base_f..][..sl],
                );
            }
        }
    }
    if slab_precollided && !fuse_copy {
        // Blocked path: the bulk copy above already brought the
        // shell-collided slab in; nothing to collide.
        return;
    }
    let mut spans: Vec<(usize, usize)> = Vec::with_capacity(8);
    for rr in 0..win.rows {
        let row_f = base_f + rr * pnx;
        let dst_row = slot * L::Q * cap_len + rr * pnx;
        // Core row behind this strip row, if any (3D strips carry the two
        // y-halo rows; 2D slabs are their own single core row).
        let core_y = if L::D == 3 {
            let yy = (win.row0 + rr) as isize - 1;
            if yy < 0 || yy >= ny as isize {
                None
            } else {
                Some(yy as usize)
            }
        } else {
            Some(s)
        };
        let collide_row_here = !slab_precollided
            && core_y.is_some_and(|y| {
                !(L::D == 3 && ((y == 0 && ctx.halo[2]) || (y == ny - 1 && ctx.halo[3])))
            });
        if fuse_copy {
            spans.clear();
            if collide_row_here {
                for_fluid_spans(
                    &ring.runs[slot * cap_rows + rr],
                    1 + x_lo,
                    1 + x_hi,
                    |a, b| spans.push((a, b)),
                );
            }
            // Copy the complement of the collide spans (halo columns, shell
            // columns, solid runs — or, for copy-only rows, everything).
            let mut cursor = 0usize;
            for &(a, b) in spans.iter().chain(std::iter::once(&(pnx, pnx))) {
                if cursor < a.min(pnx) {
                    let (c0p, c1p) = (cursor, a.min(pnx));
                    for q in 0..L::Q {
                        // SAFETY: this band owns the slot; ranges in bounds.
                        unsafe {
                            dst_raw.copy_from(
                                q * cap_len + dst_row + c0p,
                                &ctx.f[q * ctx.np + row_f + c0p..][..c1p - c0p],
                            );
                        }
                    }
                }
                cursor = b;
            }
        }
        if !collide_row_here {
            continue;
        }
        let y = core_y.expect("collide rows lie in the core");
        let (cy, cz) = if L::D == 3 { (y, s) } else { (s, 0) };
        let c0 = ctx.g.cidx(0, cy, cz);
        let rho = &ctx.rho_old[c0..c0 + nx];
        let ux = &ctx.ux_old[c0..c0 + nx];
        let uy = &ctx.uy_old[c0..c0 + nx];
        let uz = &ctx.uz_old[c0..c0 + nx];
        let ffrow = ctx.ff.map(|v| &v[c0..c0 + nx]);
        let omega_row = ctx.omega.map(|v| &v[c0..c0 + nx]);
        let dst_view = PlaneView {
            planes: dst_raw,
            stride: cap_len,
            base: dst_row + 1,
        };
        if fuse_copy {
            // Collide straight from `f` into the ring over the planned spans.
            let src_view = PlaneView {
                planes: RawSlice::new_ref(ctx.f),
                stride: ctx.np,
                base: row_f + 1,
            };
            for &(a, b) in &spans {
                // SAFETY: this slot region is owned by this band; the source
                // view is read-only; spans are in bounds.
                unsafe {
                    collide_span_dispatch::<L, T>(
                        ctx.force_on,
                        ffrow,
                        omega_row,
                        src_view,
                        dst_view,
                        a - 1,
                        b - 1,
                        rho,
                        ux,
                        uy,
                        uz,
                        ctx.kp,
                    );
                }
            }
        } else {
            // Collide the ring in place over its fluid core spans.
            for_fluid_spans(
                &ring.runs[slot * cap_rows + rr],
                1 + x_lo,
                1 + x_hi,
                |a, b| {
                    // SAFETY: this slot region is owned by this band.
                    unsafe {
                        collide_span_dispatch::<L, T>(
                            ctx.force_on,
                            ffrow,
                            omega_row,
                            dst_view,
                            dst_view,
                            a - 1,
                            b - 1,
                            rho,
                            ux,
                            uy,
                            uz,
                            ctx.kp,
                        );
                    }
                },
            );
        }
    }
}

/// One fused collide+stream+moments band over destination slabs
/// `[s0, s1)` (2D: y rows; 3D: z planes — already intersected with the step
/// range by the dispatcher). A 3D band processes its planes in y-strip
/// passes (see [`Window`]), sweeping its z range once per strip so each
/// (strip × slab) is collided exactly once per pass.
///
/// Returns the momentum-exchange force over probed solid links as
/// per-destination-row partials tagged `(z, y)`, only for rows that touched
/// probed links; each partial replays its link contributions in
/// `CpuScalar`'s `stream_row` order (destination x ascending, direction
/// ascending), and the dispatcher folds the tagged rows in `CpuScalar`'s
/// global row order — the probe diagnostic is therefore bitwise band- and
/// strip-independent.
fn fused_band<L: Lattice, T: Real>(
    ctx: &FusedCtx<'_, L, T>,
    ring: &mut Ring<T>,
    s0: usize,
    s1: usize,
) -> Vec<(u32, u32, [T; 3])> {
    let mut partials = Vec::new();
    if s0 >= s1 {
        return partials;
    }
    if L::D == 3 {
        let (ylo, yhi) = (ctx.range.lo[1], ctx.range.hi[1]);
        let cap = (yhi - ylo).min(STRIP_ROWS) + 2;
        ring.prepare(L::Q, ctx.pnx, cap);
        let mut y0 = ylo;
        while y0 < yhi {
            let y1 = (y0 + STRIP_ROWS.min(yhi - y0)).min(yhi);
            ring.set_window(Window {
                row0: y0,
                rows: y1 - y0 + 2,
            });
            sweep_window(ctx, ring, s0, s1, y0, y1, &mut partials);
            y0 = y1;
        }
    } else {
        ring.prepare(L::Q, ctx.pnx, 1);
        ring.set_window(Window { row0: 0, rows: 1 });
        sweep_window(ctx, ring, s0, s1, 0, 1, &mut partials);
    }
    partials
}

/// Sweep destination slabs `[slo, shi)` for the ring's active window
/// (destination rows `[ylo, yhi)` for 3D lattices; the slab itself for 2D).
#[allow(clippy::too_many_arguments)]
fn sweep_window<L: Lattice, T: Real>(
    ctx: &FusedCtx<'_, L, T>,
    ring: &mut Ring<T>,
    slo: usize,
    shi: usize,
    ylo: usize,
    yhi: usize,
    partials: &mut Vec<(u32, u32, [T; 3])>,
) {
    let g = ctx.g;
    let c = g.core;
    let (nx, ny) = (c[0], c[1]);
    let n_slabs = ctx.n_slabs();
    let pnx = ctx.pnx;
    let win = ring.win;
    let (cap_len, cap_rows) = (ring.cap_len, ring.cap_rows);
    let six = T::r(6.0);
    let two = T::r(2.0);
    // Solid runs of the two sweep-halo layers' window rows
    // (`[0]` = slab -1, `[1]` = slab `n_slabs`), streamed straight from
    // `f`; empty when the face has no halo.
    let (lo_face, hi_face) = if L::D == 3 { (4, 5) } else { (2, 3) };
    let mk_halo_runs = |present: bool, slab: isize| -> Vec<Vec<(u32, u32)>> {
        let mut v = vec![Vec::new(); if present { win.rows } else { 0 }];
        let base = ctx.row_base_f(slab, win.row0 as isize - 1);
        for (rr, runs) in v.iter_mut().enumerate() {
            solid_runs_row(&ctx.solid[base + rr * pnx..][..pnx], runs);
        }
        v
    };
    let halo_runs = [
        mk_halo_runs(ctx.halo[lo_face], -1),
        mk_halo_runs(ctx.halo[hi_face], n_slabs as isize),
    ];
    // Probed bounce-back links of the current row, gathered in the fused
    // pass's direction-major order and replayed in CpuScalar's cell-major
    // order: (destination x, direction q, ftot).
    let mut links: Vec<(u32, u32, T)> = Vec::new();
    let (xlo, xhi) = (ctx.range.lo[0], ctx.range.hi[0]);
    for s in slo..shi {
        // Source slabs feeding destination slab `s` (halo slabs excluded:
        // they are read straight from `f`).
        let needed: [usize; 3] = std::array::from_fn(|k| {
            let ss = s as isize + k as isize - 1;
            if ss < 0 || ss >= n_slabs as isize {
                usize::MAX
            } else {
                ss as usize
            }
        });
        for &ss in &needed {
            if ss != usize::MAX {
                ensure_slab(ctx, ring, ss, &needed);
            }
        }
        let dest_slot = ring.slot_of(s);
        // Destination rows (3D: the strip's y rows; 2D: the slab itself).
        for yy in ylo..yhi {
            // Core coordinates of this destination row.
            let (y, z) = if L::D == 3 { (yy, s) } else { (s, 0) };
            let mut pf_row = [T::zero(); 3];
            links.clear();
            let pb = g.pidx(0, y, z);
            let c0 = g.cidx(0, y, z);
            let rho_row = &ctx.rho_old[c0..c0 + nx];
            let rr_dest = if L::D == 3 { y + 1 - win.row0 } else { 0 };
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
                // Strip-local row of the source (2D: the single row).
                let rr = if L::D == 3 {
                    (sy + 1) as usize - win.row0
                } else {
                    0
                };
                // (src row data, its solid runs, padded f-offset of the row)
                let (src_row, runs_src, src_row_f): (&[T], &[(u32, u32)], usize) =
                    if ss < 0 || ss >= n_slabs as isize {
                        let hf = if ss < 0 { lo_face } else { hi_face };
                        if !ctx.halo[hf] {
                            continue;
                        }
                        let hidx = usize::from(ss >= 0);
                        let base = ctx.row_base_f(ss, sy);
                        (
                            &ctx.f[q * ctx.np + base..][..pnx],
                            &halo_runs[hidx][rr],
                            base,
                        )
                    } else {
                        let slot = ring.slot_of(ss as usize);
                        (
                            &ring.data[(slot * L::Q + q) * cap_len + rr * pnx..][..pnx],
                            &ring.runs[slot * cap_rows + rr],
                            ctx.row_base_f(ss, sy),
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
                            ctx.out
                                .copy_from(out_base + d0, &src_row[s0x..s0x + (d1 - d0)]);
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
                            [(dest_slot * L::Q + L::OPP[q]) * cap_len + rr_dest * pnx + x + 1];
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
            let dest_runs = &ring.runs[dest_slot * cap_rows + rr_dest];
            let ffrow = ctx.ff.map(|v| &v[c0..c0 + nx]);
            for_fluid_spans(dest_runs, 1 + xlo, 1 + xhi, |a, b| {
                let (x0, x1) = (a - 1, b - 1);
                // SAFETY: this band owns the destination row (out and the
                // spare moment rows alike).
                unsafe {
                    match ffrow {
                        Some(fr) => moments_span_fused::<L, T, true>(
                            ctx.out, ctx.np, pb, x0, x1, c0, ctx.rho2, ctx.ux2, ctx.uy2, ctx.uz2,
                            fr, ctx.kp,
                        ),
                        None => moments_span_fused::<L, T, false>(
                            ctx.out,
                            ctx.np,
                            pb,
                            x0,
                            x1,
                            c0,
                            ctx.rho2,
                            ctx.ux2,
                            ctx.uy2,
                            ctx.uz2,
                            &[],
                            ctx.kp,
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
            if !links.is_empty() {
                partials.push((z as u32, y as u32, pf_row));
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Backend implementation
// ---------------------------------------------------------------------------

impl<L: Lattice, T: Real> Backend<L, T> for CpuSimd {
    type Fields = SoaFields<T>;

    fn alloc(&self, sub: &Subdomain) -> SoaFields<T> {
        SoaFields::new(L::Q, sub.geom)
    }

    fn stage_in(&self, _sub: &Subdomain, fields: &mut SoaFields<T>, host: &SoaFields<T>) {
        *fields = host.clone();
    }

    fn stage_out(&self, _sub: &Subdomain, fields: &SoaFields<T>, host: &mut SoaFields<T>) {
        *host = fields.clone();
    }

    fn supports_gravity_body_force(&self) -> bool {
        true
    }

    fn exchange_f<H: HaloExchange<T>>(
        &mut self,
        exchange: &H,
        subs: &[Subdomain],
        fields: &mut [SoaFields<T>],
    ) {
        exchange.exchange_f::<L>(subs, fields);
    }

    /// Collide, in place in `f`, exactly the one-cell core boundary shell of
    /// faces with a halo neighbour — the layers the exchange packs. All
    /// other cells stay pre-collide (the fused pass collides them
    /// just-in-time).
    fn collide(&mut self, sub: &Subdomain, fields: &mut SoaFields<T>, p: &StepParams<T>) {
        if fields.bouzidi.is_some() {
            fields.fused = None;
            let mut scalar = CpuScalar::default();
            return <CpuScalar as Backend<L, T>>::collide(&mut scalar, sub, fields, p);
        }
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
        let force_on = kp.force_on(fields.force_field.is_some());
        let f = RawSlice::new(&mut fields.f);
        let (rho, ux, uy, uz) = (&fields.rho, &fields.ux, &fields.uy, &fields.uz);
        let solid = &fields.solid;
        let ff = fields.force_field.as_deref();
        let omega = fields.omega_field.as_deref();
        let full_row = |y: usize, z: usize| {
            (z == 0 && halo[4])
                || (z == nz - 1 && halo[5])
                || (y == 0 && halo[2])
                || (y == ny - 1 && halo[3])
        };
        // Full boundary rows (y/z-face layers), span-collided per row.
        let full_rows: Vec<(usize, usize)> = (0..ny * nz)
            .map(|r| (r % ny, r / ny))
            .filter(|&(y, z)| full_row(y, z))
            .collect();
        let row_body = |&(y, z): &(usize, usize)| {
            let pb = g.pidx(0, y, z);
            let c0 = g.cidx(0, y, z);
            let ffrow = ff.map(|v| &v[c0..c0 + nx]);
            let omega_row = omega.map(|v| &v[c0..c0 + nx]);
            let mut runs = Vec::new();
            solid_runs_row(&solid[pb - 1..pb - 1 + pnx], &mut runs);
            let view = PlaneView {
                planes: f,
                stride: np,
                base: pb,
            };
            for_fluid_spans(&runs, 1, 1 + nx, |a, b| {
                // SAFETY: rows are dispatched disjointly; spans in bounds.
                unsafe {
                    collide_span_dispatch::<L, T>(
                        force_on,
                        ffrow,
                        omega_row,
                        view,
                        view,
                        a - 1,
                        b - 1,
                        &rho[c0..c0 + nx],
                        &ux[c0..c0 + nx],
                        &uy[c0..c0 + nx],
                        &uz[c0..c0 + nx],
                        &kp,
                    )
                }
            });
        };
        // The row layers are tiny for 2D grids; only fork for 3D-scale work.
        #[cfg(feature = "parallel")]
        if self.use_parallel(sub) && full_rows.len() * nx >= PARALLEL_MIN_CELLS {
            full_rows.par_iter().for_each(row_body);
        } else {
            full_rows.iter().for_each(row_body);
        }
        #[cfg(not(feature = "parallel"))]
        full_rows.iter().for_each(row_body);
        // x-face layers: the remaining rows contribute one cell per column.
        // Strided single cells collide poorly, so gather blocks of them into
        // a contiguous scratch pseudo-row, span-collide it once, and scatter
        // back (fluid cells only; the amortized cost is the two copies).
        let mut cols: [Option<usize>; 2] = [None, None];
        if halo[0] {
            cols[0] = Some(0);
        }
        if halo[1] && (nx > 1 || !halo[0]) {
            cols[1] = Some(nx - 1);
        }
        let col_body = |z: usize| {
            let mut cell_pi = [0usize; BLOCK];
            let mut fg = vec![T::zero(); L::Q * BLOCK];
            let mut rho_g = [T::zero(); BLOCK];
            let mut ux_g = [T::zero(); BLOCK];
            let mut uy_g = [T::zero(); BLOCK];
            let mut uz_g = [T::zero(); BLOCK];
            let mut ff_g = [[T::zero(); 3]; BLOCK];
            let mut omega_g = [T::zero(); BLOCK];
            let mut cnt = 0usize;
            let flush = |cnt: &mut usize,
                         cell_pi: &[usize; BLOCK],
                         fg: &mut Vec<T>,
                         rho_g: &[T; BLOCK],
                         ux_g: &[T; BLOCK],
                         uy_g: &[T; BLOCK],
                         uz_g: &[T; BLOCK],
                         ff_g: &[[T; 3]; BLOCK],
                         omega_g: &[T; BLOCK]| {
                if *cnt == 0 {
                    return;
                }
                // SAFETY: the scratch is exclusive; gathered cells are
                // pairwise distinct padded indices, and z-tasks touch
                // disjoint cell sets.
                unsafe {
                    let view = PlaneView {
                        planes: RawSlice::new(fg),
                        stride: BLOCK,
                        base: 0,
                    };
                    collide_span_dispatch::<L, T>(
                        force_on,
                        ff.map(|_| &ff_g[..*cnt]),
                        omega.map(|_| &omega_g[..*cnt]),
                        view,
                        view,
                        0,
                        *cnt,
                        &rho_g[..*cnt],
                        &ux_g[..*cnt],
                        &uy_g[..*cnt],
                        &uz_g[..*cnt],
                        &kp,
                    );
                    for k in 0..*cnt {
                        for q in 0..L::Q {
                            f.set(q * np + cell_pi[k], fg[q * BLOCK + k]);
                        }
                    }
                }
                *cnt = 0;
            };
            for xc in cols.into_iter().flatten() {
                for y in 0..ny {
                    if full_row(y, z) {
                        continue;
                    }
                    let pi = g.pidx(xc, y, z);
                    if solid[pi] {
                        continue;
                    }
                    let c = g.cidx(xc, y, z);
                    cell_pi[cnt] = pi;
                    for q in 0..L::Q {
                        // SAFETY: this z-task's gather; in bounds.
                        fg[q * BLOCK + cnt] = unsafe { f.get(q * np + pi) };
                    }
                    rho_g[cnt] = rho[c];
                    ux_g[cnt] = ux[c];
                    uy_g[cnt] = uy[c];
                    uz_g[cnt] = uz[c];
                    if let Some(v) = ff {
                        ff_g[cnt] = v[c];
                    }
                    if let Some(v) = omega {
                        omega_g[cnt] = v[c];
                    }
                    cnt += 1;
                    if cnt == BLOCK {
                        flush(
                            &mut cnt, &cell_pi, &mut fg, &rho_g, &ux_g, &uy_g, &uz_g, &ff_g,
                            &omega_g,
                        );
                    }
                }
                flush(
                    &mut cnt, &cell_pi, &mut fg, &rho_g, &ux_g, &uy_g, &uz_g, &ff_g, &omega_g,
                );
            }
        };
        let col_cells = 2 * ny * nz;
        #[cfg(feature = "parallel")]
        if self.use_parallel(sub) && col_cells >= PARALLEL_MIN_CELLS / 2 {
            (0..nz).into_par_iter().for_each(col_body);
        } else {
            (0..nz).for_each(col_body);
        }
        #[cfg(not(feature = "parallel"))]
        {
            let _ = col_cells;
            (0..nz).for_each(col_body);
        }
    }

    fn stream(
        &mut self,
        sub: &Subdomain,
        fields: &mut SoaFields<T>,
        p: &StepParams<T>,
        range: CellRange,
    ) -> [T; 3] {
        if fields.bouzidi.is_some() {
            fields.fused = None;
            let mut scalar = CpuScalar::default();
            return <CpuScalar as Backend<L, T>>::stream(&mut scalar, sub, fields, p, range);
        }
        if range.is_empty() {
            return [T::zero(); 3];
        }
        let g = fields.geom;
        assert_eq!(g.halo, 1, "CpuSimd assumes one-cell halos");
        let kp = KParams::new::<L>(p);
        let halo = sub.halo_flags();
        let np = g.n_padded();
        let pnx = g.padded()[0];
        let force_on = kp.force_on(fields.force_field.is_some());
        let mut scratch = fields
            .fused
            .take()
            .unwrap_or_else(|| Box::new(FusedScratch::new(g.n_core())));
        let ctx = FusedCtx::<L, T> {
            g,
            np,
            pnx,
            pny: g.padded()[1],
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
            omega: fields.omega_field.as_deref(),
            force_on,
            kp: &kp,
            range,
            _l: std::marker::PhantomData,
        };
        // Bands partition the sweep slabs of the range (2D: y rows, 3D: z
        // planes — bands own whole `out` planes, so threads never write
        // neighbouring rows of any plane; 3D bands run in y-strip passes
        // internally, see `Window`).
        let sa = L::D - 1;
        let (blo, bhi) = (range.lo[sa], range.hi[sa]);
        let n_range = bhi - blo;
        #[cfg(feature = "parallel")]
        let threads = if self.use_parallel(sub) {
            rayon::current_num_threads()
        } else {
            1
        };
        #[cfg(not(feature = "parallel"))]
        let threads = 1;
        // One band per worker, capped on asymmetric Apple Silicon: the last
        // E-core workers extend the 3D tail more than their bands repay, and
        // every extra band also recollides two shared edge slabs. On the
        // M5 Max 128^3 f32 bench, 16 bands beat the 18-band split in the
        // same overloaded window (146.7 vs 121.7 MLUPS; see TESTING_NOTES).
        let nbands = if L::D == 3 {
            threads.min(16).clamp(1, (n_range / 4).max(1))
        } else {
            threads.clamp(1, (n_range / 16).max(1))
        };
        let band_size = n_range.div_ceil(nbands);
        // Per-band rings persist across steps in the scratch (see `Ring`).
        if scratch.rings.len() < nbands {
            scratch.rings.resize_with(nbands, Ring::default);
        }
        let rings = &mut scratch.rings[..nbands];
        let body = |band: usize, ring: &mut Ring<T>| -> Vec<(u32, u32, [T; 3])> {
            let s0 = blo + band * band_size;
            let s1 = (blo + (band + 1) * band_size).min(bhi);
            fused_band(&ctx, ring, s0, s1)
        };
        // Fold the tagged row partials in CpuScalar's global row order
        // ((z, y) lexicographic over the range). Rows without probed links
        // contribute exact zeros to CpuScalar's running sum, so folding
        // only the touched rows in that order reproduces its probe total
        // bitwise (up to signs of exact zeros) for any band/strip shape.
        let fold = |bands: Vec<Vec<(u32, u32, [T; 3])>>| -> [T; 3] {
            let mut rows: Vec<(u32, u32, [T; 3])> = bands.into_iter().flatten().collect();
            rows.sort_by_key(|&(z, y, _)| (z, y));
            rows.into_iter().fold([T::zero(); 3], |a, (_, _, b)| {
                [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
            })
        };
        #[cfg(feature = "parallel")]
        let pf = if nbands > 1 {
            fold(
                rings
                    .par_iter_mut()
                    .enumerate()
                    .map(|(b, r)| body(b, r))
                    .collect(),
            )
        } else {
            fold(
                rings
                    .iter_mut()
                    .enumerate()
                    .map(|(b, r)| body(b, r))
                    .collect(),
            )
        };
        #[cfg(not(feature = "parallel"))]
        let pf = fold(
            rings
                .iter_mut()
                .enumerate()
                .map(|(b, r)| body(b, r))
                .collect(),
        );
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

    fn apply_bouzidi(
        &mut self,
        _sub: &Subdomain,
        fields: &mut SoaFields<T>,
        p: &StepParams<T>,
    ) -> [T; 3] {
        crate::bouzidi::apply_bouzidi_impl::<L, T>(fields, p)
    }

    fn apply_volume_sources(
        &mut self,
        sub: &Subdomain,
        fields: &mut SoaFields<T>,
        p: &StepParams<T>,
    ) {
        crate::backend::apply_volume_sources_impl::<L, T>(sub, fields, p);
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
                    || !face_has_open_bc(p, face)
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
        if !p.sources.is_empty() {
            update_moments_impl::<L, T>(fields, p, self.use_parallel(sub));
            if let Some(s) = fields.fused.as_deref_mut() {
                s.fresh = false;
            }
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
/// cell is re-collided from the untouched `f` + old moments with the fused
/// pass's own collide arithmetic — the value `CpuScalar`'s in-place collide
/// leaves in the ping-pong buffer, up to TRT-pair reassociation.
fn capture_stale<L: Lattice, T: Real>(
    sub: &Subdomain,
    ctx: &FusedCtx<'_, L, T>,
    p: &StepParams<T>,
    scratch: &mut FusedScratch<T>,
) {
    let g = ctx.g;
    let np = ctx.np;
    for face in Face::ALL {
        if face.axis() >= L::D || !sub.touches_global_face(face) || !face_has_open_bc(p, face) {
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
            let omega_cell = ctx.omega.map(|v| &v[cidx..cidx + 1]);
            // SAFETY: `cell` is thread-local; single-cell span.
            unsafe {
                let view = PlaneView {
                    planes: RawSlice::new(&mut cell[..L::Q]),
                    stride: 1,
                    base: 0,
                };
                collide_span_dispatch::<L, T>(
                    ctx.force_on,
                    ffcell,
                    omega_cell,
                    view,
                    view,
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
        if face.axis() >= L::D || !sub.touches_global_face(face) || !face_has_open_bc(p, face) {
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
            let r = T::one() + dr;
            let fv = kp.force_at(ff, c, r);
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

fn face_has_open_bc<T: Real>(p: &StepParams<T>, face: Face) -> bool {
    p.faces[face.index()].is_open()
        || p.face_patches
            .iter()
            .any(|patch| patch.face == face.index() && patch.bc.is_open())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::halo::LocalPeriodic;
    use crate::lattice::{D3Q19, D3Q27};
    use crate::params::CollisionKind;
    use crate::solver::{GlobalSpec, Solver};
    use std::f64::consts::PI;

    fn omega_from_nu(nu: f64) -> f64 {
        1.0 / (3.0 * nu + 0.5)
    }

    fn cumulant_tgv<L: Lattice, B: Backend<L, f64>>(
        n: usize,
        backend: B,
    ) -> Solver<L, f64, B, LocalPeriodic> {
        let nu = 0.02;
        let spec = GlobalSpec {
            dims: [n, n, n],
            nu,
            collision: CollisionKind::Cumulant {
                omega_shear: omega_from_nu(nu),
            },
            periodic: [true, true, true],
            ..Default::default()
        };
        let mut s = Solver::new(&spec, &[], &[], [1, 1, 1], backend, LocalPeriodic);
        let u0 = 1.28e-4 / n as f64;
        s.init_with(move |x, y, z| {
            let k = 2.0 * PI / n as f64;
            let (xf, yf, zf) = (k * x as f64, k * y as f64, k * z as f64);
            let p =
                u0 * u0 / 16.0 * (((2.0 * xf).cos() + (2.0 * yf).cos()) * ((2.0 * zf).cos() + 2.0));
            (
                1.0 + 3.0 * p,
                [
                    u0 * xf.sin() * yf.cos() * zf.cos(),
                    -u0 * xf.cos() * yf.sin() * zf.cos(),
                    0.0,
                ],
            )
        });
        s
    }

    fn max_delta(a: &[f64], b: &[f64]) -> f64 {
        a.iter()
            .zip(b)
            .map(|(x, y)| (x - y).abs())
            .fold(0.0, f64::max)
    }

    fn cumulant_simd_scalar_delta<L: Lattice>() -> f64 {
        let n = 16;
        let mut scalar = cumulant_tgv::<L, _>(n, CpuScalar::default());
        let mut simd = cumulant_tgv::<L, _>(n, CpuSimd::default());
        scalar.run(200);
        simd.run(200);
        max_delta(&scalar.gather_rho(), &simd.gather_rho())
            .max(max_delta(&scalar.gather_ux(), &simd.gather_ux()))
            .max(max_delta(&scalar.gather_uy(), &simd.gather_uy()))
            .max(max_delta(&scalar.gather_uz(), &simd.gather_uz()))
    }

    #[test]
    fn cumulant_simd_matches_scalar_measured_tgv3d_tolerance() {
        let d3q19 = cumulant_simd_scalar_delta::<D3Q19>();
        let d3q27 = cumulant_simd_scalar_delta::<D3Q27>();
        eprintln!("cumulant SIMD vs scalar 200-step TGV3D: D3Q19={d3q19:e} D3Q27={d3q27:e}");
        // Measured 2026-07-06 on the Stage-3 span kernel:
        // D3Q19 = 0.0, D3Q27 = 0.0 over a 200-step f64 TGV3D.
        // The protocol's measured*10 headroom is therefore zero; keep a
        // 1e-15 floor so harmless toolchain codegen noise does not create a
        // false failure while still staying far below the 1e-8 bug line.
        assert!(d3q19 <= 1.0e-15, "D3Q19 cumulant SIMD delta {d3q19:e}");
        assert!(d3q27 <= 1.0e-15, "D3Q27 cumulant SIMD delta {d3q27:e}");
    }
}
