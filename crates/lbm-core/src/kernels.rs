//! Physics kernels, generic over [`Lattice`] and [`Real`].
//!
//! **Faithful port of V1 (`lbm-core/src/sim.rs`, retired 2026-07-05)**:
//! every floating-point expression keeps V1's exact operand order and
//! grouping. At the port's departure point the 2D (D2Q9) `f64`
//! specialisation reproduced the *pre-fusion* V1 trajectories bit-for-bit
//! (`v1_match_*` tests; the final frozen measurements live in that test's
//! header in branch history, deleted with V1). The live constraint today is
//! the backend-equivalence gate `tests/backend_simd_equiv.rs` against the
//! fused V1 `step_band` port: `f64` max |Δ| ≤ 1e-11 per cell, observed
//! ~1.6e-14 over 50 steps (TRT-pair last-ulp reassociation only).
//! Dimension loops are seeded with the first term
//! (`acc = c[0]*u[0]; for d in 1.. { acc = acc + ... }`) precisely so the
//! D = 2 case compiles to V1's `cx*vx + cy*vy` association.
//!
//! Kernels address q-major SoA planes: value of direction `q` at padded cell
//! `i` lives at `f[q * np + i]` (`np` = padded plane length).
//!
//! Row-parallel dispatch hands rows to threads; the mutable population
//! buffer is passed as a [`RawSlice`], whose safety contract is that
//! concurrent callers touch disjoint cell columns (each row is written by
//! exactly one call).

use crate::fields::LocalGeom;
use crate::lattice::{Face, Lattice, Q_MAX};
use crate::params::KParams;
use crate::real::Real;

/// Unsafely shareable mutable view for row-parallel kernels.
///
/// Invariant (enforced by the dispatchers in `backend.rs`): concurrent users
/// write disjoint index sets. All accesses are bounds-checked in debug
/// builds.
#[derive(Clone, Copy)]
pub(crate) struct RawSlice<T> {
    ptr: *mut T,
    len: usize,
}

unsafe impl<T: Send> Send for RawSlice<T> {}
unsafe impl<T: Send> Sync for RawSlice<T> {}

impl<T: Copy> RawSlice<T> {
    pub(crate) fn new(s: &mut [T]) -> Self {
        Self {
            ptr: s.as_mut_ptr(),
            len: s.len(),
        }
    }

    /// Read-only view over a shared slice.
    ///
    /// # Safety contract (caller)
    /// The returned handle must never be written through (`set`/`copy_from`);
    /// it exists so read-side kernel parameters can share the `RawSlice`
    /// plumbing.
    pub(crate) fn new_ref(s: &[T]) -> Self {
        Self {
            ptr: s.as_ptr() as *mut T,
            len: s.len(),
        }
    }

    /// Base pointer (identity comparisons only).
    #[inline(always)]
    pub(crate) fn as_ptr(self) -> *const T {
        self.ptr
    }

    /// # Safety
    /// `i < len` and no concurrent writer of index `i`.
    #[inline(always)]
    pub(crate) unsafe fn get(self, i: usize) -> T {
        debug_assert!(i < self.len);
        unsafe { *self.ptr.add(i) }
    }

    /// # Safety
    /// `i < len` and this call is the only concurrent accessor of index `i`.
    #[inline(always)]
    pub(crate) unsafe fn set(self, i: usize, v: T) {
        debug_assert!(i < self.len);
        unsafe { *self.ptr.add(i) = v }
    }

    /// Bulk copy `src` into the buffer starting at `dst_start`.
    ///
    /// # Safety
    /// `dst_start + src.len() <= len`, `src` does not overlap the target
    /// range, and this call is the only concurrent accessor of that range.
    #[inline(always)]
    pub(crate) unsafe fn copy_from(self, dst_start: usize, src: &[T]) {
        debug_assert!(dst_start + src.len() <= self.len);
        unsafe { std::ptr::copy_nonoverlapping(src.as_ptr(), self.ptr.add(dst_start), src.len()) }
    }
}

/// Equilibrium distribution in deviation form (`feq_q - w_q`), written in
/// terms of `drho = rho - 1` so no large-magnitude cancellation occurs
/// (essential for f32 deviation storage; V1 `equilibrium`).
pub(crate) fn equilibrium<L: Lattice, T: Real>(p: &KParams<T>, r: T, u: [T; 3]) -> [T; Q_MAX] {
    let three = T::r(3.0);
    let f45 = T::r(4.5);
    let f15 = T::r(1.5);
    let mut usq = u[0] * u[0];
    for d in 1..L::D {
        usq = usq + u[d] * u[d];
    }
    let drho = r - T::one();
    let mut feq = [T::zero(); Q_MAX];
    for q in 0..L::Q {
        let mut cu = p.cr[q][0] * u[0];
        for d in 1..L::D {
            cu = cu + p.cr[q][d] * u[d];
        }
        feq[q] = p.wr[q] * (drho + r * (three * cu + f45 * cu * cu - f15 * usq));
    }
    feq
}

/// TRT collision (BGK when `omega_m == omega_p`) with Guo forcing, one row of
/// core cells (V1 `collide_row`).
///
/// # Safety
/// This call must be the only concurrent writer of this row's cells in `f`.
#[allow(clippy::too_many_arguments)]
pub(crate) unsafe fn collide_row<L: Lattice, T: Real>(
    f: RawSlice<T>,
    np: usize,
    pb: usize,
    rho: &[T],
    ux: &[T],
    uy: &[T],
    uz: &[T],
    solid: &[bool],
    ff: Option<&[[T; 3]]>,
    omega: Option<&[T]>,
    p: &KParams<T>,
) {
    let three = T::r(3.0);
    let f45 = T::r(4.5);
    let f15 = T::r(1.5);
    let nine = T::r(9.0);
    let half = T::r(0.5);
    let force_on = p.force[0] != T::zero()
        || p.force[1] != T::zero()
        || p.force[2] != T::zero()
        || ff.is_some();
    for x in 0..rho.len() {
        if solid[x] {
            continue;
        }
        let r = rho[x];
        let u = [ux[x], uy[x], uz[x]];
        let fv = match ff {
            Some(field) => [
                p.force[0] + field[x][0],
                p.force[1] + field[x][1],
                p.force[2] + field[x][2],
            ],
            None => p.force,
        };
        let mut usq = u[0] * u[0];
        for d in 1..L::D {
            usq = usq + u[d] * u[d];
        }
        let mut uf = u[0] * fv[0];
        for d in 1..L::D {
            uf = uf + u[d] * fv[d];
        }
        let drho = r - T::one();
        let mut feq = [T::zero(); Q_MAX];
        let mut src = [T::zero(); Q_MAX];
        for q in 0..L::Q {
            let mut cu = p.cr[q][0] * u[0];
            for d in 1..L::D {
                cu = cu + p.cr[q][d] * u[d];
            }
            feq[q] = p.wr[q] * (drho + r * (three * cu + f45 * cu * cu - f15 * usq));
            if force_on {
                let mut cf = p.cr[q][0] * fv[0];
                for d in 1..L::D {
                    cf = cf + p.cr[q][d] * fv[d];
                }
                src[q] = p.wr[q] * (three * (cf - uf) + nine * cu * cf);
            }
        }
        let i = pb + x;
        // SAFETY: row-disjoint dispatch (see RawSlice contract).
        unsafe {
            let op = omega.map_or(p.omega_p, |v| v[x]);
            let cp = if omega.is_some() {
                T::one() - op / T::r(2.0)
            } else {
                p.cp
            };
            let i0 = L::REST * np + i;
            let f0 = f.get(i0);
            f.set(i0, f0 - op * (f0 - feq[L::REST]) + cp * src[L::REST]);
            for &(a, b) in L::PAIRS {
                let (ia, ib) = (a * np + i, b * np + i);
                let (fa, fb) = (f.get(ia), f.get(ib));
                let fp = half * (fa + fb);
                let fm = half * (fa - fb);
                let ep = half * (feq[a] + feq[b]);
                let em = half * (feq[a] - feq[b]);
                let sp = half * (src[a] + src[b]);
                let sm = half * (src[a] - src[b]);
                let rp = op * (fp - ep);
                let rm = p.omega_m * (fm - em);
                f.set(ia, fa - rp - rm + cp * sp + p.cm * sm);
                f.set(ib, fb - rp + rm + cp * sp - p.cm * sm);
            }
        }
    }
}

fn central_basis<L: Lattice>() -> [[u8; 3]; Q_MAX] {
    let mut basis = [[0u8; 3]; Q_MAX];
    let mut n = 0usize;
    for ax in 0..=2 {
        for ay in 0..=2 {
            for az in 0..=2 {
                if L::D == 2 && az != 0 {
                    continue;
                }
                if L::D == 3 && L::Q == 19 && ax > 0 && ay > 0 && az > 0 {
                    continue;
                }
                if n < L::Q {
                    basis[n] = [ax, ay, az];
                    n += 1;
                }
            }
        }
    }
    debug_assert_eq!(n, L::Q);
    basis
}

#[inline]
fn pow_upto2(x: f64, e: u8) -> f64 {
    match e {
        0 => 1.0,
        1 => x,
        2 => x * x,
        _ => unreachable!("central basis only uses powers 0..=2"),
    }
}

fn central_phi<L: Lattice>(q: usize, exp: [u8; 3], u: [f64; 3]) -> f64 {
    let c = L::C[q];
    let mut v = pow_upto2(c[0] as f64 - u[0], exp[0]);
    v *= pow_upto2(c[1] as f64 - u[1], exp[1]);
    if L::D == 3 {
        v *= pow_upto2(c[2] as f64 - u[2], exp[2]);
    }
    v
}

fn central_equilibrium(exp: [u8; 3], rho: f64, d: usize) -> f64 {
    let mut v = rho;
    for &e in exp.iter().take(d) {
        v *= match e {
            0 => 1.0,
            1 => 0.0,
            2 => 1.0 / 3.0,
            _ => unreachable!("central basis only uses powers 0..=2"),
        };
    }
    v
}

fn solve_moment_system<L: Lattice>(
    basis: &[[u8; 3]; Q_MAX],
    u: [f64; 3],
    rhs: &[f64; Q_MAX],
) -> [f64; Q_MAX] {
    let n = L::Q;
    let mut a = [[0.0f64; Q_MAX + 1]; Q_MAX];
    for m in 0..n {
        for q in 0..n {
            a[m][q] = central_phi::<L>(q, basis[m], u);
        }
        a[m][n] = rhs[m];
    }
    for col in 0..n {
        let mut pivot = col;
        for row in col + 1..n {
            if a[row][col].abs() > a[pivot][col].abs() {
                pivot = row;
            }
        }
        assert!(
            a[pivot][col].abs() > 1.0e-14,
            "singular central-moment basis for lattice Q{}",
            L::Q
        );
        if pivot != col {
            a.swap(pivot, col);
        }
        let inv = 1.0 / a[col][col];
        for j in col..=n {
            a[col][j] *= inv;
        }
        for row in 0..n {
            if row == col {
                continue;
            }
            let factor = a[row][col];
            if factor == 0.0 {
                continue;
            }
            for j in col..=n {
                a[row][j] -= factor * a[col][j];
            }
        }
    }
    let mut out = [0.0f64; Q_MAX];
    for i in 0..n {
        out[i] = a[i][n];
    }
    out
}

/// Cascaded central-moment collision with Guo forcing, one row of core cells.
///
/// This is the stage-2 CPU scalar reference for `CollisionKind::Cumulant`.
/// It is not a logarithmic cumulant implementation: populations are
/// transformed to central moments, second-order deviatoric moments relax with
/// the per-cell shear rate, the second-order trace and all higher-order
/// moments relax at 1.0, and the Guo source is transformed through the same
/// central-moment basis with `(I - S/2)` prefactors.
///
/// # Safety
/// This call must be the only concurrent writer of this row's cells in `f`.
#[allow(clippy::too_many_arguments)]
pub(crate) unsafe fn collide_row_central_moment<L: Lattice, T: Real>(
    f: RawSlice<T>,
    np: usize,
    pb: usize,
    rho: &[T],
    ux: &[T],
    uy: &[T],
    uz: &[T],
    solid: &[bool],
    ff: Option<&[[T; 3]]>,
    omega: Option<&[T]>,
    p: &KParams<T>,
) {
    let basis = central_basis::<L>();
    for x in 0..rho.len() {
        if solid[x] {
            continue;
        }
        let r = rho[x].as_f64();
        let u = [ux[x].as_f64(), uy[x].as_f64(), uz[x].as_f64()];
        let fv_t = match ff {
            Some(field) => [
                p.force[0] + field[x][0],
                p.force[1] + field[x][1],
                p.force[2] + field[x][2],
            ],
            None => p.force,
        };
        let fv = [fv_t[0].as_f64(), fv_t[1].as_f64(), fv_t[2].as_f64()];
        let force_on = fv[0] != 0.0 || fv[1] != 0.0 || fv[2] != 0.0;
        let i = pb + x;
        let mut phys = [0.0f64; Q_MAX];
        let mut src = [0.0f64; Q_MAX];
        let feq_rest = u.iter().take(L::D).all(|&v| v == 0.0) && !force_on;
        let mut exact_rest_equilibrium = feq_rest;
        for q in 0..L::Q {
            // SAFETY: row-disjoint dispatch (see RawSlice contract).
            let fq = unsafe { f.get(q * np + i) };
            phys[q] = fq.as_f64() + L::W[q];
            if exact_rest_equilibrium {
                let expected = T::r(L::W[q] * (r - 1.0));
                exact_rest_equilibrium &= fq.as_f64().to_bits() == expected.as_f64().to_bits();
            }
        }
        if exact_rest_equilibrium {
            continue;
        }
        if force_on {
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
                src[q] = L::W[q] * (3.0 * (cf - uf) + 9.0 * cu * cf);
            }
        }

        let mut mom = [0.0f64; Q_MAX];
        let mut src_mom = [0.0f64; Q_MAX];
        let mut eq = [0.0f64; Q_MAX];
        for m in 0..L::Q {
            eq[m] = central_equilibrium(basis[m], r, L::D);
            for q in 0..L::Q {
                let phi = central_phi::<L>(q, basis[m], u);
                mom[m] += phi * phys[q];
                src_mom[m] += phi * src[q];
            }
        }

        let os = omega.map_or(p.omega_shear.as_f64(), |v| v[x].as_f64());
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
            // SAFETY: row-disjoint dispatch (see RawSlice contract).
            unsafe { f.set(q * np + i, T::r(out_phys[q] - L::W[q])) };
        }
    }
}

/// Pull-scheme streaming for one destination row segment `x0..x1` (V1
/// `stream_row`). Returns the momentum-exchange force accumulated over
/// probed solid links.
///
/// A source cell beyond the core is read from the halo when the crossed face
/// has a neighbour; otherwise the direction is skipped and the out-buffer
/// slot keeps its prior content (the unknown-population mechanics the
/// open-face BCs rely on).
///
/// # Safety
/// This call must be the only concurrent writer of this row segment in `out`.
#[allow(clippy::too_many_arguments)]
pub(crate) unsafe fn stream_row<L: Lattice, T: Real>(
    out: RawSlice<T>,
    f: &[T],
    np: usize,
    geom: &LocalGeom,
    halo: [bool; 6],
    y: usize,
    z: usize,
    x0: usize,
    x1: usize,
    solid: &[bool],
    wall_u: &[[T; 3]],
    rho_row: &[T],
    probe: Option<&[bool]>,
    p: &KParams<T>,
) -> [T; 3] {
    let six = T::r(6.0);
    let two = T::r(2.0);
    let mut pf = [T::zero(); 3];
    let pb = geom.pidx(0, y, z);
    for x in x0..x1 {
        let i = pb + x;
        if solid[i] {
            continue;
        }
        'dirs: for q in 0..L::Q {
            let c = L::C[q];
            let s = [
                x as isize - c[0] as isize,
                y as isize - c[1] as isize,
                z as isize - c[2] as isize,
            ];
            // V1's wrap-or-skip per axis, generalised to halo-or-skip per
            // face (the halo holds the wrapped/neighbour values).
            for a in 0..L::D {
                if s[a] < 0 {
                    if !halo[2 * a] {
                        continue 'dirs;
                    }
                } else if s[a] >= geom.core[a] as isize && !halo[2 * a + 1] {
                    continue 'dirs;
                }
            }
            let si = geom.pidx_i(s[0], s[1], s[2]);
            if solid[si] {
                // Half-way bounce-back off the wall between cells si and i,
                // with momentum injection for moving walls. In deviation
                // storage the formula is unchanged (w_q == w_opp(q)).
                let fout = f[L::OPP[q] * np + i];
                let wu = wall_u[si];
                let mut cu = p.cr[q][0] * wu[0];
                for d in 1..L::D {
                    cu = cu + p.cr[q][d] * wu[d];
                }
                let fin = fout + six * p.wr[q] * rho_row[x] * cu;
                // SAFETY: row-disjoint dispatch.
                unsafe { out.set(q * np + i, fin) };
                if let Some(mask) = probe {
                    if mask[si] {
                        // Momentum given to the wall through this link, using
                        // physical populations (deviation + weight); V1
                        // probe convention.
                        let ftot = fout + fin + two * p.wr[q];
                        pf[0] = pf[0] - p.cr[q][0] * ftot;
                        pf[1] = pf[1] - p.cr[q][1] * ftot;
                        pf[2] = pf[2] - p.cr[q][2] * ftot;
                    }
                }
            } else {
                // SAFETY: row-disjoint dispatch.
                unsafe { out.set(q * np + i, f[q * np + si]) };
            }
        }
    }
    pf
}

/// Recompute macroscopic fields from the populations for one row (V1
/// `moments_row`). Writes the compact (halo-free) moment rows.
#[allow(clippy::too_many_arguments)]
pub(crate) fn moments_row<L: Lattice, T: Real>(
    f: &[T],
    np: usize,
    pb: usize,
    rho: &mut [T],
    ux: &mut [T],
    uy: &mut [T],
    uz: &mut [T],
    solid: &[bool],
    ff: Option<&[[T; 3]]>,
    p: &KParams<T>,
) {
    let half = T::r(0.5);
    for x in 0..rho.len() {
        if solid[x] {
            continue;
        }
        let i = pb + x;
        // Deviation storage: rho = 1 + sum(f_dev); sum(w c) = 0 so the
        // momentum needs no correction.
        let mut dr = T::zero();
        let mut m = [T::zero(); 3];
        for q in 0..L::Q {
            let fq = f[q * np + i];
            dr = dr + fq;
            m[0] = m[0] + p.cr[q][0] * fq;
            m[1] = m[1] + p.cr[q][1] * fq;
            if L::D == 3 {
                m[2] = m[2] + p.cr[q][2] * fq;
            }
        }
        let fv = match ff {
            Some(field) => [
                p.force[0] + field[x][0],
                p.force[1] + field[x][1],
                p.force[2] + field[x][2],
            ],
            None => p.force,
        };
        let r = T::one() + dr;
        rho[x] = r;
        let inv = T::one() / r;
        ux[x] = (m[0] + half * fv[0]) * inv;
        uy[x] = (m[1] + half * fv[1]) * inv;
        if L::D == 3 {
            uz[x] = (m[2] + half * fv[2]) * inv;
        }
    }
}

// ---------------------------------------------------------------------------
// Open-face boundary conditions (run on the post-swap buffer, core cells of
// the face plane only; V1 apply_open_edges family)
// ---------------------------------------------------------------------------

/// Zou–He kind: prescribed velocity or prescribed density.
pub(crate) enum ZhKind<T: Real> {
    Velocity([T; 3]),
    Pressure(T),
}

/// Iterate the core cells of a face plane in canonical order: the two
/// tangent axes ascending, lower axis innermost. For 2D faces this is V1's
/// `side_cells` order (`y` ascending on X faces, `x` ascending on Y faces).
/// The along-face index passed to the callback enumerates cells in this
/// order (2D: the single tangent coordinate — V1's profile coordinate).
pub(crate) fn for_face_cells(
    geom: &LocalGeom,
    face: Face,
    mut body: impl FnMut(usize, [usize; 3]),
) {
    for_face_cells_selected(geom, face, FaceCellSelection::All, |coord, pos| {
        body(coord, pos)
    });
}

#[derive(Clone, Copy)]
pub(crate) enum FaceCellSelection<'a> {
    All,
    Rect {
        lo: [usize; 2],
        hi: [usize; 2],
    },
    Excluding {
        rects: &'a [([usize; 2], [usize; 2])],
    },
}

fn selected(sel: FaceCellSelection<'_>, c1: usize, c2: usize) -> bool {
    match sel {
        FaceCellSelection::All => true,
        FaceCellSelection::Rect { lo, hi } => {
            lo[0] <= c1 && c1 <= hi[0] && lo[1] <= c2 && c2 <= hi[1]
        }
        FaceCellSelection::Excluding { rects } => !rects
            .iter()
            .any(|(lo, hi)| lo[0] <= c1 && c1 <= hi[0] && lo[1] <= c2 && c2 <= hi[1]),
    }
}

pub(crate) fn for_face_cells_selected(
    geom: &LocalGeom,
    face: Face,
    selection: FaceCellSelection<'_>,
    mut body: impl FnMut(usize, [usize; 3]),
) {
    let a = face.axis();
    let fixed = if face.is_neg() { 0 } else { geom.core[a] - 1 };
    let (t1, t2) = face.tangents();
    let mut coord = 0;
    for c2 in 0..geom.core[t2] {
        for c1 in 0..geom.core[t1] {
            if !selected(selection, c1, c2) {
                coord += 1;
                continue;
            }
            let mut pos = [0usize; 3];
            pos[a] = fixed;
            pos[t1] = c1;
            pos[t2] = c2;
            body(coord, pos);
            coord += 1;
        }
    }
}

/// Zou–He boundary parameterised by the face normal (V1 `zou_he`,
/// generalised to D3Q19 faces).
///
/// ## Derivation (both lattices, any axis-aligned face)
///
/// After streaming, a face cell is missing exactly the populations entering
/// the domain, `{q : c_q·n = +1}` with `n` the inward normal (3 for D2Q9,
/// 5 for D3Q19: `n` itself plus `n ± t` per tangent axis `t`). Split the
/// known populations into `S0 = Σ f (c·n = 0)` and `S⁻ = Σ f (c·n = −1)`.
/// Imposing mass `Σ f = rho` and normal momentum `Σ (c·n) f = rho u·n` and
/// eliminating the unknown sums gives the closure (independent of `D`):
///
/// ```text
/// rho (1 - u·n) = S0 + 2 S⁻      → rho   (velocity BC)
///                                → u·n   (pressure BC, u.t = 0)
/// ```
///
/// Each unknown is then reconstructed as bounce-back of its opposite's
/// non-equilibrium part (Zou & He 1997; D3Q19: Hecht & Harting 2010,
/// J. Stat. Mech. P01018):
///
/// ```text
/// f_q = f_q̄ + (feq_q - feq_q̄) + corr_q,   feq_q - feq_q̄ = 6 w_q rho (c_q·u)
/// ```
///
/// Pure NEBB (`corr = 0`) already satisfies mass and normal momentum via the
/// closure identity (`Σ_unknown 6 w_q rho c_q·u = rho u·n` exactly, by the
/// 2nd-moment weight identities). The tangential momentum constraint per
/// tangent axis `t`, `Σ (c·t) f = rho u_t`, fixes an antisymmetric
/// correction `± N_t` on the diagonal pair `n ± t`. Substituting NEBB into
/// the constraint, the known negative-diagonal terms `f_{-(n±t)}` cancel
/// against the bounced ones and what remains is
///
/// ```text
/// N_t = (1/3) rho u_t - (1/2) Q_t,   Q_t = Σ_{c·n = 0} (c·t) f
/// ```
///
/// (`Q_t` is the transverse momentum carried by the face-parallel
/// populations: in 2D the pair `±t`, in 3D `±t` and the four in-plane
/// diagonals `±t ± s`.) The reconstruction, with `w_axis`, `w_diag` the
/// axis/diagonal weights:
///
/// ```text
/// f_n     = f_-n     + 6 w_axis rho (u·n)             (2D: 2/3, 3D: 1/3)
/// f_{n±t} = f_{-n∓t} + 6 w_diag rho (u·n ± u_t) ± N_t (w_diag = 1/36 both)
/// ```
///
/// For D2Q9 this regroups exactly into V1's form
/// `f_{n±t} = f_{-n∓t} + (1/6) rho (u·n) ± [(1/2) rho u_t - (1/2) Q_t]`;
/// the 2D branch below keeps V1's operand order bit-for-bit.
///
/// Deviation storage: the physical `S0 + 2 S⁻` equals the deviation sums
/// plus `+1` (the per-face closure constant, tested in `lattice.rs`); the
/// reconstruction and `Q_t` are weight-neutral (`w_q = w_q̄`, and the
/// face-parallel weights cancel pairwise in `Q_t`), so the formulas apply
/// to deviations unchanged.
///
/// Under a z-invariant field with `u_z = 0`, the D3Q19 reconstruction
/// projects exactly onto the D2Q9 one (sum populations over `c_z`): the
/// z-corrections vanish by z-reflection symmetry and the x/y formulas add up
/// to V1's — the T15 degeneracy tests pin this down numerically.
pub(crate) fn zou_he_face_selected<L: Lattice, T: Real>(
    f: &mut [T],
    np: usize,
    geom: &LocalGeom,
    solid: &[bool],
    face: Face,
    kind: &ZhKind<T>,
    profile: Option<&[[T; 3]]>,
    selection: FaceCellSelection<'_>,
) {
    if L::D == 3 {
        return zou_he_face_3d::<L, T>(f, np, geom, solid, face, kind, profile, selection);
    }
    let n = face.n_in();
    let (nxi, nyi) = (n[0] as i32, n[1] as i32);
    let (tx, ty) = (-nyi, nxi);
    let q_n = L::dir_index([nxi as i8, nyi as i8, 0]);
    let q_d1 = L::dir_index([(nxi + tx) as i8, (nyi + ty) as i8, 0]);
    let q_d2 = L::dir_index([(nxi - tx) as i8, (nyi - ty) as i8, 0]);
    let q_t = L::dir_index([tx as i8, ty as i8, 0]);
    let q_mt = L::dir_index([-tx as i8, -ty as i8, 0]);
    let (half, c23, c16, two) = (T::r(0.5), T::r(2.0 / 3.0), T::r(1.0 / 6.0), T::r(2.0));
    let (nxr, nyr) = (T::r(nxi as f64), T::r(nyi as f64));
    let (txr, tyr) = (T::r(tx as f64), T::r(ty as f64));
    for_face_cells_selected(geom, face, selection, |coord, pos| {
        let i = geom.pidx(pos[0], pos[1], pos[2]);
        if solid[i] {
            return;
        }
        let s0 = f[L::REST * np + i] + f[q_t * np + i] + f[q_mt * np + i];
        let sneg = f[L::OPP[q_n] * np + i] + f[L::OPP[q_d1] * np + i] + f[L::OPP[q_d2] * np + i];
        let closure = s0 + two * sneg + T::one();
        let (r, un, ut) = match *kind {
            ZhKind::Velocity(u) => {
                let u = profile.map_or(u, |p| p[coord]);
                let un = u[0] * nxr + u[1] * nyr;
                let ut = u[0] * txr + u[1] * tyr;
                (closure / (T::one() - un), un, ut)
            }
            ZhKind::Pressure(rho_bc) => {
                // From the closure rho (1 - u.n) = S0 + 2 S-.
                let un = T::one() - closure / rho_bc;
                (rho_bc, un, T::zero())
            }
        };
        let tcorr = half * (r * ut - (f[q_t * np + i] - f[q_mt * np + i]));
        f[q_n * np + i] = f[L::OPP[q_n] * np + i] + c23 * r * un;
        f[q_d1 * np + i] = f[L::OPP[q_d1] * np + i] + c16 * r * un + tcorr;
        f[q_d2 * np + i] = f[L::OPP[q_d2] * np + i] + c16 * r * un - tcorr;
    });
}

/// D3Q19 branch of [`zou_he_face`] (see its docs for the shared derivation).
///
/// Unknowns for inward normal `n` and tangent axes `t1 < t2`:
/// `n`, `n±t1`, `n±t2`. Reconstruction:
///
/// ```text
/// f_n      = f_-n      + (1/3) rho (u·n)
/// f_{n±tk} = f_{-n∓tk} + (1/6) rho (u·n ± u_tk) ± N_k
/// N_k      = (1/3) rho u_tk - (1/2) Q_tk
/// Q_t1     = f_{t1} - f_{-t1} + f_{t1+t2} + f_{t1-t2} - f_{-t1+t2} - f_{-t1-t2}
/// Q_t2     = f_{t2} - f_{-t2} + f_{t1+t2} - f_{t1-t2} + f_{-t1+t2} - f_{-t1-t2}
/// ```
#[allow(clippy::too_many_arguments)]
fn zou_he_face_3d<L: Lattice, T: Real>(
    f: &mut [T],
    np: usize,
    geom: &LocalGeom,
    solid: &[bool],
    face: Face,
    kind: &ZhKind<T>,
    profile: Option<&[[T; 3]]>,
    selection: FaceCellSelection<'_>,
) {
    // This branch reconstructs exactly 5 unknowns (the D3Q19 in-face set).
    // A lattice with a different unknown count (e.g. a future D3Q27, planned
    // with Q_MAX = 27) would silently leave slots unreconstructed here — and
    // `dir_index` would still resolve the 5 named directions — so guard it
    // rather than trust the caller (A-8). `unknowns(face)` is a const-derived
    // table, so this is effectively a compile-time expectation checked once
    // per face.
    assert_eq!(
        L::unknowns(face).len(),
        5,
        "zou_he_face_3d hardcodes 5 unknowns; lattice {:?} face has {}",
        std::any::type_name::<L>(),
        L::unknowns(face).len()
    );
    let a = face.axis();
    let (t1, t2) = face.tangents();
    let n = face.n_in();
    let unit = |axis: usize, s: i8| -> [i8; 3] {
        let mut v = [0i8; 3];
        v[axis] = s;
        v
    };
    let add = |p: [i8; 3], q: [i8; 3]| [p[0] + q[0], p[1] + q[1], p[2] + q[2]];
    // The 5 unknowns (c·n = +1).
    let q_n = L::dir_index(n);
    let q_p1 = L::dir_index(add(n, unit(t1, 1)));
    let q_m1 = L::dir_index(add(n, unit(t1, -1)));
    let q_p2 = L::dir_index(add(n, unit(t2, 1)));
    let q_m2 = L::dir_index(add(n, unit(t2, -1)));
    // The 8 non-rest face-parallel directions (c·n = 0): S0 and Q_t terms.
    let q_t1 = L::dir_index(unit(t1, 1));
    let q_mt1 = L::dir_index(unit(t1, -1));
    let q_t2 = L::dir_index(unit(t2, 1));
    let q_mt2 = L::dir_index(unit(t2, -1));
    let q_pp = L::dir_index(add(unit(t1, 1), unit(t2, 1)));
    let q_pm = L::dir_index(add(unit(t1, 1), unit(t2, -1)));
    let q_mp = L::dir_index(add(unit(t1, -1), unit(t2, 1)));
    let q_mm = L::dir_index(add(unit(t1, -1), unit(t2, -1)));
    let (half, c13, c16, two) = (T::r(0.5), T::r(1.0 / 3.0), T::r(1.0 / 6.0), T::r(2.0));
    // u·n_in for the single non-zero normal component (±1).
    let nsign = T::r(n[a] as f64);
    for_face_cells_selected(geom, face, selection, |coord, pos| {
        let i = geom.pidx(pos[0], pos[1], pos[2]);
        if solid[i] {
            return;
        }
        let s0 = f[L::REST * np + i]
            + f[q_t1 * np + i]
            + f[q_mt1 * np + i]
            + f[q_t2 * np + i]
            + f[q_mt2 * np + i]
            + f[q_pp * np + i]
            + f[q_pm * np + i]
            + f[q_mp * np + i]
            + f[q_mm * np + i];
        let sneg = f[L::OPP[q_n] * np + i]
            + f[L::OPP[q_p1] * np + i]
            + f[L::OPP[q_m1] * np + i]
            + f[L::OPP[q_p2] * np + i]
            + f[L::OPP[q_m2] * np + i];
        let closure = s0 + two * sneg + T::one();
        let (r, un, ut1, ut2) = match *kind {
            ZhKind::Velocity(u) => {
                let u = profile.map_or(u, |p| p[coord]);
                let un = u[a] * nsign;
                (closure / (T::one() - un), un, u[t1], u[t2])
            }
            ZhKind::Pressure(rho_bc) => {
                // From the closure rho (1 - u·n) = S0 + 2 S⁻.
                let un = T::one() - closure / rho_bc;
                (rho_bc, un, T::zero(), T::zero())
            }
        };
        // Transverse fluxes of the face-parallel populations.
        let qt1 = f[q_t1 * np + i] - f[q_mt1 * np + i] + f[q_pp * np + i] + f[q_pm * np + i]
            - f[q_mp * np + i]
            - f[q_mm * np + i];
        let qt2 = f[q_t2 * np + i] - f[q_mt2 * np + i] + f[q_pp * np + i] - f[q_pm * np + i]
            + f[q_mp * np + i]
            - f[q_mm * np + i];
        let n1 = c13 * r * ut1 - half * qt1;
        let n2 = c13 * r * ut2 - half * qt2;
        f[q_n * np + i] = f[L::OPP[q_n] * np + i] + c13 * r * un;
        f[q_p1 * np + i] = f[L::OPP[q_p1] * np + i] + c16 * r * (un + ut1) + n1;
        f[q_m1 * np + i] = f[L::OPP[q_m1] * np + i] + c16 * r * (un - ut1) - n1;
        f[q_p2 * np + i] = f[L::OPP[q_p2] * np + i] + c16 * r * (un + ut2) + n2;
        f[q_m2 * np + i] = f[L::OPP[q_m2] * np + i] + c16 * r * (un - ut2) - n2;
    });
}

/// Zero-gradient outflow: copy the unknown populations from the cell one
/// step inward along the face normal (V1 `outflow`). Generic over D.
pub(crate) fn outflow_face_selected<L: Lattice, T: Real>(
    f: &mut [T],
    np: usize,
    geom: &LocalGeom,
    solid: &[bool],
    face: Face,
    selection: FaceCellSelection<'_>,
) {
    let n = face.n_in();
    let unknowns = L::unknowns(face);
    for_face_cells_selected(geom, face, selection, |_, pos| {
        let i = geom.pidx(pos[0], pos[1], pos[2]);
        let j = geom.pidx_i(
            pos[0] as isize + n[0] as isize,
            pos[1] as isize + n[1] as isize,
            pos[2] as isize + n[2] as isize,
        );
        if solid[i] || solid[j] {
            return;
        }
        for &q in unknowns {
            f[q * np + i] = f[q * np + j];
        }
    });
}

/// Convective (radiation) outflow with mass pinning (V1 `convective_outflow`).
///
/// In the pull scheme the unknown slots at the face still hold the previous
/// step's post-collide populations after streaming, so
/// `f(edge,t+1) = (f(edge,t) + Uc f(interior,t+1)) / (1 + Uc)` needs no extra
/// storage. The mass correction pins `rho(edge)` to `rho(neighbour)` by
/// distributing the deficit over the unknowns by weight.
pub(crate) fn convective_face_selected<L: Lattice, T: Real>(
    f: &mut [T],
    np: usize,
    geom: &LocalGeom,
    solid: &[bool],
    face: Face,
    u_conv: T,
    selection: FaceCellSelection<'_>,
) {
    let n = face.n_in();
    let unknowns = L::unknowns(face);
    let lam = u_conv;
    let inv = T::one() / (T::one() + lam);
    // Weight share for the mass correction over the unknown links (f64 sum
    // first, then one conversion — V1 order).
    let mut ws = 0.0f64;
    for &q in unknowns {
        ws += L::W[q];
    }
    let wsum = T::r(ws);
    for_face_cells_selected(geom, face, selection, |_, pos| {
        let i = geom.pidx(pos[0], pos[1], pos[2]);
        let j = geom.pidx_i(
            pos[0] as isize + n[0] as isize,
            pos[1] as isize + n[1] as isize,
            pos[2] as isize + n[2] as isize,
        );
        if solid[i] || solid[j] {
            return;
        }
        for &q in unknowns {
            let prev = f[q * np + i];
            f[q * np + i] = (prev + lam * f[q * np + j]) * inv;
        }
        let mut di = T::zero();
        let mut dj = T::zero();
        for q in 0..L::Q {
            di = di + f[q * np + i];
            dj = dj + f[q * np + j];
        }
        let corr = dj - di;
        for &q in unknowns {
            f[q * np + i] = f[q * np + i] + corr * T::r(L::W[q]) / wsum;
        }
    });
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::lattice::{D2Q9, D3Q19, D3Q27};
    use crate::params::{FaceBC, StepParams};
    use crate::real::Real;

    /// A-10f: `equilibrium()` and the feq computed inline by `collide_row`
    /// must be bit-identical. Property: a state seeded exactly with
    /// `equilibrium()`'s output is a **bit-exact fixed point** of forceless
    /// collision — `f - ω (f − feq_inline)` returns `f` only when
    /// `feq_inline` matches the stored `equilibrium()` values bitwise (the
    /// pair half-sums then cancel exactly, so the relaxation terms are 0.0,
    /// not merely small). Any drift between the two formula copies (operand
    /// order, factored constants) breaks the equality for some sample.
    fn collide_fixed_point_on_equilibrium<L: Lattice, T: Real>() {
        let ncells = 64usize;
        let params = StepParams::<T> {
            omega_p: 1.25,
            omega_m: 0.7,
            force: [T::zero(); 3],
            faces: [FaceBC::Closed; 6],
            sources: Vec::new(),
            face_patches: Vec::new(),
        };
        let kp = KParams::new::<L>(&params);

        // Deterministic varied samples across the low-Mach envelope.
        let sample = |i: usize| -> (T, [T; 3]) {
            let h = |k: f64| ((i as f64 + 1.0) * k).sin();
            let r = T::r(1.0 + 0.5 * h(12.9898));
            let mut u = [T::r(0.12 * h(78.233)), T::r(0.12 * h(39.425)), T::zero()];
            if L::D == 3 {
                u[2] = T::r(0.12 * h(93.989));
            }
            (r, u)
        };

        let np = ncells; // single unpadded row: pb = 0
        let mut f = vec![T::zero(); L::Q * np];
        let mut rho = vec![T::zero(); ncells];
        let mut ux = vec![T::zero(); ncells];
        let mut uy = vec![T::zero(); ncells];
        let mut uz = vec![T::zero(); ncells];
        for i in 0..ncells {
            let (r, u) = sample(i);
            rho[i] = r;
            ux[i] = u[0];
            uy[i] = u[1];
            uz[i] = u[2];
            let feq = equilibrium::<L, T>(&kp, r, u);
            for q in 0..L::Q {
                f[q * np + i] = feq[q];
            }
        }
        let before = f.clone();
        let solid = vec![false; ncells];
        let raw = RawSlice::new(&mut f);
        // SAFETY: single-threaded, sole writer of the row.
        unsafe {
            collide_row::<L, T>(raw, np, 0, &rho, &ux, &uy, &uz, &solid, None, None, &kp);
        }
        for (k, (a, b)) in before.iter().zip(&f).enumerate() {
            // Bit equality via the exact f64 promotion (exact for f32 too,
            // and sign-of-zero preserving).
            assert!(
                a.as_f64().to_bits() == b.as_f64().to_bits(),
                "equilibrium fixed point broken at slot {k} (q = {}, cell = {}): \
                 collide's inline feq differs from equilibrium()",
                k / np,
                k % np
            );
        }
    }

    #[test]
    fn collide_feq_matches_equilibrium_bitwise_d2q9() {
        collide_fixed_point_on_equilibrium::<D2Q9, f64>();
        collide_fixed_point_on_equilibrium::<D2Q9, f32>();
    }

    #[test]
    fn collide_feq_matches_equilibrium_bitwise_d3q19() {
        collide_fixed_point_on_equilibrium::<D3Q19, f64>();
        collide_fixed_point_on_equilibrium::<D3Q19, f32>();
    }

    fn cumulant_params<T: Real>() -> KParams<T> {
        let params = StepParams::<T> {
            omega_p: 1.25,
            omega_m: -1.25,
            force: [T::zero(); 3],
            faces: [FaceBC::Closed; 6],
            sources: Vec::new(),
            face_patches: Vec::new(),
        };
        KParams::new::<D3Q19>(&params)
    }

    fn collide_cumulant_rest_fixed_point<L: Lattice>() {
        let ncells = 8usize;
        let params = StepParams::<f64> {
            omega_p: 1.25,
            omega_m: -1.25,
            force: [0.0; 3],
            faces: [FaceBC::Closed; 6],
            sources: Vec::new(),
            face_patches: Vec::new(),
        };
        let kp = KParams::new::<L>(&params);
        let mut f = vec![0.0; L::Q * ncells];
        let rho = vec![1.0; ncells];
        let ux = vec![0.0; ncells];
        let uy = vec![0.0; ncells];
        let uz = vec![0.0; ncells];
        let solid = vec![false; ncells];
        let before = f.clone();
        let raw = RawSlice::new(&mut f);
        // SAFETY: single-threaded, sole writer of the row.
        unsafe {
            collide_row_central_moment::<L, f64>(
                raw, ncells, 0, &rho, &ux, &uy, &uz, &solid, None, None, &kp,
            );
        }
        assert_eq!(before, f);
    }

    #[test]
    fn cumulant_rest_equilibrium_is_exact_fixed_point() {
        collide_cumulant_rest_fixed_point::<D3Q19>();
        collide_cumulant_rest_fixed_point::<D3Q27>();
    }

    #[test]
    fn cumulant_uniform_velocity_stays_uniform_after_collide() {
        let ncells = 16usize;
        let params = StepParams::<f64> {
            omega_p: 1.1,
            omega_m: -1.1,
            force: [0.0; 3],
            faces: [FaceBC::Closed; 6],
            sources: Vec::new(),
            face_patches: Vec::new(),
        };
        let kp = KParams::new::<D3Q19>(&params);
        let rho = vec![1.03; ncells];
        let ux = vec![0.04; ncells];
        let uy = vec![-0.03; ncells];
        let uz = vec![0.02; ncells];
        let solid = vec![false; ncells];
        let feq = equilibrium::<D3Q19, f64>(&kp, rho[0], [ux[0], uy[0], uz[0]]);
        let mut f = vec![0.0; D3Q19::Q * ncells];
        for x in 0..ncells {
            for q in 0..D3Q19::Q {
                f[q * ncells + x] = feq[q];
            }
        }
        let raw = RawSlice::new(&mut f);
        // SAFETY: single-threaded, sole writer of the row.
        unsafe {
            collide_row_central_moment::<D3Q19, f64>(
                raw, ncells, 0, &rho, &ux, &uy, &uz, &solid, None, None, &kp,
            );
        }
        for q in 0..D3Q19::Q {
            let v0 = f[q * ncells];
            for x in 1..ncells {
                assert_eq!(v0.to_bits(), f[q * ncells + x].to_bits());
            }
        }
    }

    #[test]
    fn cumulant_d3q19_tgv3d_short_run_finite_and_decays() {
        let n = 8usize;
        let np = n * n * n;
        let kp = cumulant_params::<f64>();
        let mut f = vec![0.0; D3Q19::Q * np];
        let mut rho = vec![1.0; np];
        let mut ux = vec![0.0; np];
        let mut uy = vec![0.0; np];
        let mut uz = vec![0.0; np];
        let solid = vec![false; np];
        for z in 0..n {
            for y in 0..n {
                for x in 0..n {
                    let i = (z * n + y) * n + x;
                    let sx = (2.0 * std::f64::consts::PI * x as f64 / n as f64).sin();
                    let cx = (2.0 * std::f64::consts::PI * x as f64 / n as f64).cos();
                    let sy = (2.0 * std::f64::consts::PI * y as f64 / n as f64).sin();
                    let cy = (2.0 * std::f64::consts::PI * y as f64 / n as f64).cos();
                    ux[i] = 0.02 * sx * cy;
                    uy[i] = -0.02 * cx * sy;
                    let feq = equilibrium::<D3Q19, f64>(&kp, 1.0, [ux[i], uy[i], 0.0]);
                    for q in 0..D3Q19::Q {
                        f[q * np + i] = feq[q];
                    }
                }
            }
        }
        let energy = |fx: &[f64]| -> f64 {
            let mut rho_m = vec![0.0; np];
            let mut ux_m = vec![0.0; np];
            let mut uy_m = vec![0.0; np];
            let mut uz_m = vec![0.0; np];
            moments_row::<D3Q19, f64>(
                fx, np, 0, &mut rho_m, &mut ux_m, &mut uy_m, &mut uz_m, &solid, None, &kp,
            );
            ux_m.iter()
                .zip(&uy_m)
                .zip(&uz_m)
                .map(|((&a, &b), &c)| a * a + b * b + c * c)
                .sum::<f64>()
        };
        let e0 = energy(&f);
        for _ in 0..8 {
            for z in 0..n {
                for y in 0..n {
                    let pb = (z * n + y) * n;
                    let raw = RawSlice::new(&mut f);
                    // SAFETY: serial row update.
                    unsafe {
                        collide_row_central_moment::<D3Q19, f64>(
                            raw,
                            np,
                            pb,
                            &rho[pb..pb + n],
                            &ux[pb..pb + n],
                            &uy[pb..pb + n],
                            &uz[pb..pb + n],
                            &solid[pb..pb + n],
                            None,
                            None,
                            &kp,
                        );
                    }
                }
            }
            let mut out = vec![0.0; D3Q19::Q * np];
            for q in 0..D3Q19::Q {
                let c = D3Q19::C[q];
                for z in 0..n {
                    let sz = (z + n - (c[2] as isize).rem_euclid(n as isize) as usize) % n;
                    for y in 0..n {
                        let sy = (y + n - (c[1] as isize).rem_euclid(n as isize) as usize) % n;
                        for x in 0..n {
                            let sx = (x + n - (c[0] as isize).rem_euclid(n as isize) as usize) % n;
                            let dst = (z * n + y) * n + x;
                            let src = (sz * n + sy) * n + sx;
                            out[q * np + dst] = f[q * np + src];
                        }
                    }
                }
            }
            f = out;
            let mut uz_next = vec![0.0; np];
            moments_row::<D3Q19, f64>(
                &f,
                np,
                0,
                &mut rho,
                &mut ux,
                &mut uy,
                &mut uz_next,
                &solid,
                None,
                &kp,
            );
            uz = uz_next;
        }
        let e1 = energy(&f);
        assert!(e1.is_finite());
        assert!(e1 < e0, "TGV kinetic energy did not decay: {e0} -> {e1}");
    }
}
