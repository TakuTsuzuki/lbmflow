//! Physics kernels, generic over [`Lattice`] and [`Real`].
//!
//! **Faithful port of V1 (`lbm-core/src/sim.rs`)**: every floating-point
//! expression keeps V1's exact operand order and grouping, so the 2D (D2Q9)
//! specialisation reproduces V1 trajectories bit-for-bit in `f64` (verified
//! by the `v1_match_*` tests). Dimension loops are seeded with the first term
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
            let i0 = L::REST * np + i;
            let f0 = f.get(i0);
            f.set(
                i0,
                f0 - p.omega_p * (f0 - feq[L::REST]) + p.cp * src[L::REST],
            );
            for &(a, b) in L::PAIRS {
                let (ia, ib) = (a * np + i, b * np + i);
                let (fa, fb) = (f.get(ia), f.get(ib));
                let fp = half * (fa + fb);
                let fm = half * (fa - fb);
                let ep = half * (feq[a] + feq[b]);
                let em = half * (feq[a] - feq[b]);
                let sp = half * (src[a] + src[b]);
                let sm = half * (src[a] - src[b]);
                let rp = p.omega_p * (fp - ep);
                let rm = p.omega_m * (fm - em);
                f.set(ia, fa - rp - rm + p.cp * sp + p.cm * sm);
                f.set(ib, fb - rp + rm + p.cp * sp - p.cm * sm);
            }
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
    let a = face.axis();
    let fixed = if face.is_neg() { 0 } else { geom.core[a] - 1 };
    let (t1, t2) = face.tangents();
    let mut coord = 0;
    for c2 in 0..geom.core[t2] {
        for c1 in 0..geom.core[t1] {
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
pub(crate) fn zou_he_face<L: Lattice, T: Real>(
    f: &mut [T],
    np: usize,
    geom: &LocalGeom,
    solid: &[bool],
    face: Face,
    kind: &ZhKind<T>,
    profile: Option<&[[T; 3]]>,
) {
    if L::D == 3 {
        return zou_he_face_3d::<L, T>(f, np, geom, solid, face, kind, profile);
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
    for_face_cells(geom, face, |coord, pos| {
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
    for_face_cells(geom, face, |coord, pos| {
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
pub(crate) fn outflow_face<L: Lattice, T: Real>(
    f: &mut [T],
    np: usize,
    geom: &LocalGeom,
    solid: &[bool],
    face: Face,
) {
    let n = face.n_in();
    let unknowns = L::unknowns(face);
    for_face_cells(geom, face, |_, pos| {
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
pub(crate) fn convective_face<L: Lattice, T: Real>(
    f: &mut [T],
    np: usize,
    geom: &LocalGeom,
    solid: &[bool],
    face: Face,
    u_conv: T,
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
    for_face_cells(geom, face, |_, pos| {
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
