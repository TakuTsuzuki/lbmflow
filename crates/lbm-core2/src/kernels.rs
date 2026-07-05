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
            f.set(i0, f0 - p.omega_p * (f0 - feq[L::REST]) + p.cp * src[L::REST]);
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
fn for_face_cells(geom: &LocalGeom, face: Face, mut body: impl FnMut(usize, [usize; 3])) {
    let a = face.axis();
    let fixed = if face.is_neg() { 0 } else { geom.core[a] - 1 };
    let (t1, t2) = match a {
        0 => (1, 2),
        1 => (0, 2),
        _ => (0, 1),
    };
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

/// Zou–He boundary parameterised by the face normal (V1 `zou_he`).
///
/// With inward normal `n` and tangent `t = rot90(n)`, the three unknown
/// populations after streaming are `n`, `n+t`, `n-t`:
///
/// ```text
/// rho     = (S0 + 2 S-) / (1 - u.n)          (velocity BC)
/// u.n     = 1 - (S0 + 2 S-) / rho            (pressure BC, u.t = 0)
/// f_n     = f_-n     + (2/3) rho (u.n)
/// f_{n±t} = f_{-n∓t} + (1/6) rho (u.n) ± [ (1/2) rho (u.t) - (1/2) T ]
/// ```
///
/// Deviation storage: the physical `S0 + 2 S-` equals the deviation sums
/// plus `+1` (the per-face closure constant, tested in `lattice.rs`).
///
/// The D3Q19 face has 5 unknowns and needs transverse-momentum corrections;
/// it lands with M-C (T15). Until then this asserts `D == 2`.
pub(crate) fn zou_he_face<L: Lattice, T: Real>(
    f: &mut [T],
    np: usize,
    geom: &LocalGeom,
    solid: &[bool],
    face: Face,
    kind: &ZhKind<T>,
    profile: Option<&[[T; 3]]>,
) {
    assert!(
        L::D == 2,
        "Zou-He closure for 3D lattices lands with M-C (T15)"
    );
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
        let sneg =
            f[L::OPP[q_n] * np + i] + f[L::OPP[q_d1] * np + i] + f[L::OPP[q_d2] * np + i];
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
