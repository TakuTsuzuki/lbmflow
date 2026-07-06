//! Lattice definitions (D2Q9, D3Q19, D3Q27) behind a compile-time `Lattice` trait.
//!
//! The direction ordering of each lattice is the single source of truth for
//! the whole project (V1 principle, now per-lattice). Derived tables — the
//! opposite-direction map, TRT pairs and per-face unknown sets — are computed
//! from `C` by `const fn` at compile time, so they can never drift from the
//! velocity table.
//!
//! Design deviation from docs/ARCHITECTURE_V2.md §2.1 (documented per the
//! "improvements allowed" clause): the sketch `const C: [[i8; 3]; Self::Q]`
//! requires the unstable `generic_const_exprs` feature, so the tables are
//! exposed as `&'static` slices instead, each promoted from a `const fn`
//! computed array. Lengths are checked at compile time (const asserts) and by
//! unit tests.

/// One of the six axis-aligned faces of a rectangular domain.
///
/// `index()` ordering: XNeg=0, XPos=1, YNeg=2, YPos=3, ZNeg=4, ZPos=5.
/// In V1's 2D terms: Left=XNeg, Right=XPos, Bottom=YNeg, Top=YPos.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Face {
    /// `x = 0` plane (V1 "Left").
    XNeg,
    /// `x = nx - 1` plane (V1 "Right").
    XPos,
    /// `y = 0` plane (V1 "Bottom").
    YNeg,
    /// `y = ny - 1` plane (V1 "Top").
    YPos,
    /// `z = 0` plane.
    ZNeg,
    /// `z = nz - 1` plane.
    ZPos,
}

impl Face {
    /// All six faces in index order.
    pub const ALL: [Face; 6] = [
        Face::XNeg,
        Face::XPos,
        Face::YNeg,
        Face::YPos,
        Face::ZNeg,
        Face::ZPos,
    ];

    /// Table index (0..6).
    #[inline]
    pub const fn index(self) -> usize {
        self as usize
    }

    /// Axis this face is perpendicular to (0 = x, 1 = y, 2 = z).
    #[inline]
    pub const fn axis(self) -> usize {
        self.index() / 2
    }

    /// Whether this is the low-coordinate face of its axis.
    #[inline]
    pub const fn is_neg(self) -> bool {
        self.index() % 2 == 0
    }

    /// Inward-pointing unit normal.
    #[inline]
    pub const fn n_in(self) -> [i8; 3] {
        let mut n = [0i8; 3];
        n[self.axis()] = if self.is_neg() { 1 } else { -1 };
        n
    }

    /// The two tangent axes of this face, ascending (the canonical
    /// along-face ordering used by face-cell iteration, inlet profiles and
    /// halo layers: `t1` varies fastest).
    #[inline]
    pub const fn tangents(self) -> (usize, usize) {
        match self.axis() {
            0 => (1, 2),
            1 => (0, 2),
            _ => (0, 1),
        }
    }

    /// The opposing face of the same axis.
    #[inline]
    pub const fn opposite(self) -> Face {
        match self {
            Face::XNeg => Face::XPos,
            Face::XPos => Face::XNeg,
            Face::YNeg => Face::YPos,
            Face::YPos => Face::YNeg,
            Face::ZNeg => Face::ZPos,
            Face::ZPos => Face::ZNeg,
        }
    }
}

/// A discrete velocity set (DdQq lattice), fully known at compile time.
///
/// All tables are `&'static` slices of length `Q` (see module docs for why
/// they are not fixed-size associated arrays). Invariants, enforced by
/// `const fn` construction and unit tests:
///
/// - `C[REST] == [0, 0, 0]` and `W[REST]` is the largest weight.
/// - `C[OPP[q]] == -C[q]`, `OPP[OPP[q]] == q`, `W[OPP[q]] == W[q]`.
/// - `PAIRS` lists every non-rest direction exactly once as `(q, OPP[q])`
///   with `q < OPP[q]`, in ascending `q` order (the TRT sweep order).
/// - `unknowns(face)` is `{ q : c_q · n_in(face) > 0 }` in ascending order —
///   the populations entering the domain through that face, i.e. the ones a
///   boundary condition must reconstruct after streaming (3 per face for
///   D2Q9, 5 for D3Q19, 9 for D3Q27) and the ones a halo exchange must
///   transfer.
/// - For every face: `Σ_{c·n=0} w + 2 Σ_{c·n<0} w == 1` exactly — the
///   closure constant that lets deviation-form Zou–He add `+1` instead of
///   summing weights (V1 convention, generalises to D3Q19).
pub trait Lattice: Copy + Send + Sync + 'static {
    /// Spatial dimension (2 or 3).
    const D: usize;
    /// Number of discrete velocities.
    const Q: usize;
    /// Discrete velocities; `C[q][2] == 0` for 2D lattices. Length `Q`.
    const C: &'static [[i8; 3]];
    /// Quadrature weights. Length `Q`.
    const W: &'static [f64];
    /// Opposite-direction indices. Length `Q`.
    const OPP: &'static [usize];
    /// Index of the rest velocity.
    const REST: usize = 0;
    /// TRT direction pairs `(q, OPP[q])`, `q < OPP[q]`, ascending. Length `(Q-1)/2`.
    const PAIRS: &'static [(usize, usize)];
    /// Squared lattice speed of sound, `cs^2 = 1/3`.
    const CS2: f64 = 1.0 / 3.0;
    /// Per-face unknown sets, indexed by [`Face::index`]. Faces along axes
    /// `>= D` have empty sets.
    const FACE_UNKNOWNS: [&'static [usize]; 6];

    /// Populations entering through `face` (unknown after streaming there).
    #[inline]
    fn unknowns(face: Face) -> &'static [usize] {
        Self::FACE_UNKNOWNS[face.index()]
    }

    /// Find the direction index with the given velocity.
    ///
    /// Panics if `c` is not a velocity of this lattice.
    fn dir_index(c: [i8; 3]) -> usize {
        (0..Self::Q)
            .find(|&q| Self::C[q] == c)
            .unwrap_or_else(|| panic!("({},{},{}) is not a lattice direction", c[0], c[1], c[2]))
    }
}

// ---------------------------------------------------------------------------
// const fn derivations
// ---------------------------------------------------------------------------

const fn opp_table<const Q: usize>(c: &[[i8; 3]; Q]) -> [usize; Q] {
    let mut opp = [usize::MAX; Q];
    let mut q = 0;
    while q < Q {
        let mut r = 0;
        while r < Q {
            if c[r][0] == -c[q][0] && c[r][1] == -c[q][1] && c[r][2] == -c[q][2] {
                opp[q] = r;
            }
            r += 1;
        }
        assert!(opp[q] != usize::MAX, "direction has no opposite");
        q += 1;
    }
    opp
}

const fn pairs_table<const Q: usize, const P: usize>(
    opp: &[usize; Q],
    rest: usize,
) -> [(usize, usize); P] {
    let mut pairs = [(0usize, 0usize); P];
    let mut n = 0;
    let mut q = 0;
    while q < Q {
        if q != rest && q < opp[q] {
            pairs[n] = (q, opp[q]);
            n += 1;
        }
        q += 1;
    }
    assert!(n == P, "TRT pair count mismatch");
    pairs
}

const fn face_unknowns<const Q: usize, const N: usize>(
    c: &[[i8; 3]; Q],
    n_in: [i8; 3],
) -> [usize; N] {
    let mut out = [0usize; N];
    let mut n = 0;
    let mut q = 0;
    while q < Q {
        let dot = c[q][0] as i32 * n_in[0] as i32
            + c[q][1] as i32 * n_in[1] as i32
            + c[q][2] as i32 * n_in[2] as i32;
        if dot > 0 {
            out[n] = q;
            n += 1;
        }
        q += 1;
    }
    assert!(n == N, "face unknown count mismatch");
    out
}

// ---------------------------------------------------------------------------
// D2Q9
// ---------------------------------------------------------------------------

/// D2Q9 lattice. Direction ordering (identical to V1 `lbm_core::lattice`):
///
/// ```text
///   6  2  5
///   3  0  1
///   7  4  8
/// ```
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct D2Q9;

const D2Q9_C: [[i8; 3]; 9] = [
    [0, 0, 0],
    [1, 0, 0],
    [0, 1, 0],
    [-1, 0, 0],
    [0, -1, 0],
    [1, 1, 0],
    [-1, 1, 0],
    [-1, -1, 0],
    [1, -1, 0],
];

const D2Q9_W: [f64; 9] = [
    4.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
];

const D2Q9_OPP: [usize; 9] = opp_table(&D2Q9_C);
const D2Q9_PAIRS: [(usize, usize); 4] = pairs_table(&D2Q9_OPP, 0);
const D2Q9_UNK_XN: [usize; 3] = face_unknowns(&D2Q9_C, [1, 0, 0]);
const D2Q9_UNK_XP: [usize; 3] = face_unknowns(&D2Q9_C, [-1, 0, 0]);
const D2Q9_UNK_YN: [usize; 3] = face_unknowns(&D2Q9_C, [0, 1, 0]);
const D2Q9_UNK_YP: [usize; 3] = face_unknowns(&D2Q9_C, [0, -1, 0]);

impl Lattice for D2Q9 {
    const D: usize = 2;
    const Q: usize = 9;
    const C: &'static [[i8; 3]] = &D2Q9_C;
    const W: &'static [f64] = &D2Q9_W;
    const OPP: &'static [usize] = &D2Q9_OPP;
    const PAIRS: &'static [(usize, usize)] = &D2Q9_PAIRS;
    const FACE_UNKNOWNS: [&'static [usize]; 6] = [
        &D2Q9_UNK_XN,
        &D2Q9_UNK_XP,
        &D2Q9_UNK_YN,
        &D2Q9_UNK_YP,
        &[],
        &[],
    ];
}

// ---------------------------------------------------------------------------
// D3Q19
// ---------------------------------------------------------------------------

/// D3Q19 lattice, standard ordering (Krüger et al., *The Lattice Boltzmann
/// Method*, 2017): rest, the 6 axis directions (+x, −x, +y, −y, +z, −z),
/// then the 12 edge diagonals grouped in opposite pairs:
///
/// ```text
/// q:  0    1    2    3    4    5    6
/// c: 000  +00  -00  0+0  0-0  00+  00-
/// q:  7    8    9   10   11   12   13   14   15   16   17   18
/// c: ++0  --0  +0+  -0-  0++  0--  +-0  -+0  +0-  -0+  0+-  0-+
/// ```
///
/// `OPP[q]` is `q±1` throughout, which the const derivation reproduces.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct D3Q19;

const D3Q19_C: [[i8; 3]; 19] = [
    [0, 0, 0],
    [1, 0, 0],
    [-1, 0, 0],
    [0, 1, 0],
    [0, -1, 0],
    [0, 0, 1],
    [0, 0, -1],
    [1, 1, 0],
    [-1, -1, 0],
    [1, 0, 1],
    [-1, 0, -1],
    [0, 1, 1],
    [0, -1, -1],
    [1, -1, 0],
    [-1, 1, 0],
    [1, 0, -1],
    [-1, 0, 1],
    [0, 1, -1],
    [0, -1, 1],
];

const D3Q19_W: [f64; 19] = [
    1.0 / 3.0,
    1.0 / 18.0,
    1.0 / 18.0,
    1.0 / 18.0,
    1.0 / 18.0,
    1.0 / 18.0,
    1.0 / 18.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
];

const D3Q19_OPP: [usize; 19] = opp_table(&D3Q19_C);
const D3Q19_PAIRS: [(usize, usize); 9] = pairs_table(&D3Q19_OPP, 0);
const D3Q19_UNK_XN: [usize; 5] = face_unknowns(&D3Q19_C, [1, 0, 0]);
const D3Q19_UNK_XP: [usize; 5] = face_unknowns(&D3Q19_C, [-1, 0, 0]);
const D3Q19_UNK_YN: [usize; 5] = face_unknowns(&D3Q19_C, [0, 1, 0]);
const D3Q19_UNK_YP: [usize; 5] = face_unknowns(&D3Q19_C, [0, -1, 0]);
const D3Q19_UNK_ZN: [usize; 5] = face_unknowns(&D3Q19_C, [0, 0, 1]);
const D3Q19_UNK_ZP: [usize; 5] = face_unknowns(&D3Q19_C, [0, 0, -1]);

impl Lattice for D3Q19 {
    const D: usize = 3;
    const Q: usize = 19;
    const C: &'static [[i8; 3]] = &D3Q19_C;
    const W: &'static [f64] = &D3Q19_W;
    const OPP: &'static [usize] = &D3Q19_OPP;
    const PAIRS: &'static [(usize, usize)] = &D3Q19_PAIRS;
    const FACE_UNKNOWNS: [&'static [usize]; 6] = [
        &D3Q19_UNK_XN,
        &D3Q19_UNK_XP,
        &D3Q19_UNK_YN,
        &D3Q19_UNK_YP,
        &D3Q19_UNK_ZN,
        &D3Q19_UNK_ZP,
    ];
}

// ---------------------------------------------------------------------------
// D3Q27
// ---------------------------------------------------------------------------

/// D3Q27 lattice. The first 19 directions are exactly [`D3Q19`]; the eight
/// body-diagonal corner directions are appended as adjacent opposite pairs:
///
/// ```text
/// q: 19   20   21   22   23   24   25   26
/// c: +++  ---  ++-  --+  +-+  -+-  +--  -++
/// ```
///
/// This preserves the project-wide `OPP[q] == q±1` convention for every
/// moving direction while adding the corner links used by the full tensor
/// product quadrature. Halo exchange forwards those corner populations through
/// the existing extended face layers.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct D3Q27;

const D3Q27_C: [[i8; 3]; 27] = [
    [0, 0, 0],
    [1, 0, 0],
    [-1, 0, 0],
    [0, 1, 0],
    [0, -1, 0],
    [0, 0, 1],
    [0, 0, -1],
    [1, 1, 0],
    [-1, -1, 0],
    [1, 0, 1],
    [-1, 0, -1],
    [0, 1, 1],
    [0, -1, -1],
    [1, -1, 0],
    [-1, 1, 0],
    [1, 0, -1],
    [-1, 0, 1],
    [0, 1, -1],
    [0, -1, 1],
    [1, 1, 1],
    [-1, -1, -1],
    [1, 1, -1],
    [-1, -1, 1],
    [1, -1, 1],
    [-1, 1, -1],
    [1, -1, -1],
    [-1, 1, 1],
];

const D3Q27_W: [f64; 27] = [
    8.0 / 27.0,
    2.0 / 27.0,
    2.0 / 27.0,
    2.0 / 27.0,
    2.0 / 27.0,
    2.0 / 27.0,
    2.0 / 27.0,
    1.0 / 54.0,
    1.0 / 54.0,
    1.0 / 54.0,
    1.0 / 54.0,
    1.0 / 54.0,
    1.0 / 54.0,
    1.0 / 54.0,
    1.0 / 54.0,
    1.0 / 54.0,
    1.0 / 54.0,
    1.0 / 54.0,
    1.0 / 54.0,
    1.0 / 216.0,
    1.0 / 216.0,
    1.0 / 216.0,
    1.0 / 216.0,
    1.0 / 216.0,
    1.0 / 216.0,
    1.0 / 216.0,
    1.0 / 216.0,
];

const D3Q27_OPP: [usize; 27] = opp_table(&D3Q27_C);
const D3Q27_PAIRS: [(usize, usize); 13] = pairs_table(&D3Q27_OPP, 0);
const D3Q27_UNK_XN: [usize; 9] = face_unknowns(&D3Q27_C, [1, 0, 0]);
const D3Q27_UNK_XP: [usize; 9] = face_unknowns(&D3Q27_C, [-1, 0, 0]);
const D3Q27_UNK_YN: [usize; 9] = face_unknowns(&D3Q27_C, [0, 1, 0]);
const D3Q27_UNK_YP: [usize; 9] = face_unknowns(&D3Q27_C, [0, -1, 0]);
const D3Q27_UNK_ZN: [usize; 9] = face_unknowns(&D3Q27_C, [0, 0, 1]);
const D3Q27_UNK_ZP: [usize; 9] = face_unknowns(&D3Q27_C, [0, 0, -1]);

impl Lattice for D3Q27 {
    const D: usize = 3;
    const Q: usize = 27;
    const C: &'static [[i8; 3]] = &D3Q27_C;
    const W: &'static [f64] = &D3Q27_W;
    const OPP: &'static [usize] = &D3Q27_OPP;
    const PAIRS: &'static [(usize, usize)] = &D3Q27_PAIRS;
    const FACE_UNKNOWNS: [&'static [usize]; 6] = [
        &D3Q27_UNK_XN,
        &D3Q27_UNK_XP,
        &D3Q27_UNK_YN,
        &D3Q27_UNK_YP,
        &D3Q27_UNK_ZN,
        &D3Q27_UNK_ZP,
    ];
}

/// Maximum `Q` over the supported lattices; kernel-local scratch arrays are
/// sized with this so they never allocate.
pub const Q_MAX: usize = 27;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::backend::CpuScalar;
    use crate::backend_simd::CpuSimd;
    use crate::halo::LocalPeriodic;
    use crate::params::CollisionKind;
    use crate::solver::{GlobalSpec, Solver};
    use std::f64::consts::PI;

    fn check_basic_invariants<L: Lattice>() {
        assert_eq!(L::C.len(), L::Q);
        assert_eq!(L::W.len(), L::Q);
        assert_eq!(L::OPP.len(), L::Q);
        assert_eq!(L::PAIRS.len(), (L::Q - 1) / 2);
        assert_eq!(L::C[L::REST], [0, 0, 0]);
        for q in 0..L::Q {
            for a in 0..3 {
                assert_eq!(L::C[L::OPP[q]][a], -L::C[q][a]);
                if a >= L::D {
                    assert_eq!(L::C[q][a], 0, "component {a} of dir {q} must be 0");
                }
            }
            assert_eq!(L::OPP[L::OPP[q]], q);
            assert_eq!(L::W[L::OPP[q]], L::W[q]);
            assert!(L::W[q] > 0.0);
        }
        // PAIRS cover every non-rest direction exactly once, ascending.
        let mut seen = vec![false; L::Q];
        seen[L::REST] = true;
        let mut last = 0;
        for &(a, b) in L::PAIRS {
            assert!(a < b, "pair ({a},{b}) not ordered");
            assert_eq!(b, L::OPP[a]);
            assert!(a >= last, "pairs not in ascending order");
            last = a;
            assert!(!seen[a] && !seen[b], "direction listed twice");
            seen[a] = true;
            seen[b] = true;
        }
        assert!(seen.iter().all(|&s| s), "pairs must cover all directions");
    }

    fn check_moments<L: Lattice>() {
        // 0th: sum w = 1.
        let sw: f64 = L::W.iter().sum();
        assert!((sw - 1.0).abs() < 1e-15, "sum w = {sw}");
        // 1st and 3rd (odd) moments vanish; 2nd = cs^2 I.
        for a in 0..3 {
            let m1: f64 = (0..L::Q).map(|q| L::W[q] * L::C[q][a] as f64).sum();
            assert!(m1.abs() < 1e-15, "first moment [{a}] = {m1}");
            for b in 0..3 {
                let m2: f64 = (0..L::Q)
                    .map(|q| L::W[q] * (L::C[q][a] as i32 * L::C[q][b] as i32) as f64)
                    .sum();
                let expect = if a == b && a < L::D { L::CS2 } else { 0.0 };
                assert!(
                    (m2 - expect).abs() < 1e-15,
                    "second moment [{a}][{b}] = {m2}, expect {expect}"
                );
                for c in 0..3 {
                    let m3: f64 = (0..L::Q)
                        .map(|q| {
                            L::W[q]
                                * (L::C[q][a] as i32 * L::C[q][b] as i32 * L::C[q][c] as i32) as f64
                        })
                        .sum();
                    assert!(m3.abs() < 1e-15, "third moment [{a}][{b}][{c}] = {m3}");
                }
            }
        }
    }

    /// 4th-order isotropy: sum w c_a c_b c_c c_d = cs^4 (d_ab d_cd + d_ac d_bd
    /// + d_ad d_bc) over the active axes — the condition for a correct
    /// Navier–Stokes viscous stress.
    fn check_isotropy<L: Lattice>() {
        let d = |i: usize, j: usize| if i == j { 1.0 } else { 0.0 };
        for a in 0..L::D {
            for b in 0..L::D {
                for c in 0..L::D {
                    for e in 0..L::D {
                        let m4: f64 = (0..L::Q)
                            .map(|q| {
                                let cv = L::C[q];
                                L::W[q]
                                    * (cv[a] as i32 * cv[b] as i32 * cv[c] as i32 * cv[e] as i32)
                                        as f64
                            })
                            .sum();
                        let expect = L::CS2
                            * L::CS2
                            * (d(a, b) * d(c, e) + d(a, c) * d(b, e) + d(a, e) * d(b, c));
                        assert!(
                            (m4 - expect).abs() < 1e-15,
                            "4th moment [{a}{b}{c}{e}] = {m4}, expect {expect}"
                        );
                    }
                }
            }
        }
    }

    fn check_face_unknowns<L: Lattice>() {
        for face in Face::ALL {
            let unk = L::unknowns(face);
            let n = face.n_in();
            if face.axis() >= L::D {
                assert!(unk.is_empty());
                continue;
            }
            let expected_count = (0..L::Q)
                .filter(|&q| (0..3).map(|a| L::C[q][a] as i32 * n[a] as i32).sum::<i32>() > 0)
                .count();
            assert_eq!(unk.len(), expected_count, "{face:?}");
            let mut last = 0;
            for &q in unk {
                let dot: i32 = (0..3).map(|a| L::C[q][a] as i32 * n[a] as i32).sum();
                assert!(dot > 0, "{face:?}: dir {q} does not enter the domain");
                assert!(q >= last, "{face:?}: unknowns not ascending");
                last = q;
            }
            // Deviation-form Zou–He closure constant: for any straight face,
            // sum(w | c.n = 0) + 2 sum(w | c.n < 0) == 1 exactly.
            let mut s = 0.0;
            for q in 0..L::Q {
                let dot: i32 = (0..3).map(|a| L::C[q][a] as i32 * n[a] as i32).sum();
                if dot == 0 {
                    s += L::W[q];
                } else if dot < 0 {
                    s += 2.0 * L::W[q];
                }
            }
            assert!((s - 1.0).abs() < 1e-15, "{face:?}: closure constant {s}");
        }
    }

    fn moment2<L: Lattice>(a: usize, b: usize, c: usize) -> f64 {
        (0..L::Q)
            .map(|q| {
                L::W[q]
                    * (L::C[q][a] as i32
                        * L::C[q][a] as i32
                        * L::C[q][b] as i32
                        * L::C[q][b] as i32
                        * L::C[q][c] as i32
                        * L::C[q][c] as i32) as f64
            })
            .sum()
    }

    fn moment_pow<L: Lattice>(powers: [u32; 3]) -> f64 {
        (0..L::Q)
            .map(|q| {
                let mut p = 1.0;
                for (a, pow) in powers.into_iter().enumerate() {
                    p *= (L::C[q][a] as f64).powi(pow as i32);
                }
                L::W[q] * p
            })
            .sum()
    }

    #[test]
    fn d2q9_invariants() {
        check_basic_invariants::<D2Q9>();
        check_moments::<D2Q9>();
        check_isotropy::<D2Q9>();
        check_face_unknowns::<D2Q9>();
    }

    #[test]
    fn d3q19_invariants() {
        check_basic_invariants::<D3Q19>();
        check_moments::<D3Q19>();
        check_isotropy::<D3Q19>();
        check_face_unknowns::<D3Q19>();
    }

    #[test]
    fn d3q27_invariants() {
        check_basic_invariants::<D3Q27>();
        check_moments::<D3Q27>();
        check_isotropy::<D3Q27>();
        check_face_unknowns::<D3Q27>();
    }

    #[test]
    fn d3q27_opposites_are_adjacent_and_d3q19_is_prefix() {
        for q in 0..D3Q19::Q {
            assert_eq!(D3Q27::C[q], D3Q19::C[q]);
        }
        for q in (1..27).step_by(2) {
            assert_eq!(D3Q27::OPP[q], q + 1);
            assert_eq!(D3Q27::OPP[q + 1], q);
        }
    }

    #[test]
    fn d3q27_corner_moments_distinguish_from_d3q19() {
        let cs6 = D3Q27::CS2 * D3Q27::CS2 * D3Q27::CS2;

        // D3Q27's tensor-product corners restore the diagonal xyz identity
        // Σ w cx² cy² cz² = cs^6. D3Q19 has no body diagonals, so this moment
        // is exactly zero there.
        assert!((moment2::<D3Q27>(0, 1, 2) - cs6).abs() < 1e-15);
        assert_eq!(moment2::<D3Q19>(0, 1, 2), 0.0);

        // Sixth-order isotropy at the diagonal-discrete level has the tensor
        // ratio <x^6>:<x^4 y^2>:<x^2 y^2 z^2> = 9:3:1. D3Q19 matches the
        // first two levels but fails the third because the corners are absent.
        assert!((moment_pow::<D3Q27>([6, 0, 0]) - 9.0 * cs6).abs() < 1e-15);
        assert!((moment_pow::<D3Q27>([4, 2, 0]) - 3.0 * cs6).abs() < 1e-15);
        assert!((moment_pow::<D3Q27>([2, 2, 2]) - cs6).abs() < 1e-15);
        assert_eq!(moment_pow::<D3Q19>([2, 2, 2]), 0.0);
    }

    fn run_d3q27_periodic_tgv<B>(backend: B)
    where
        B: crate::backend::Backend<D3Q27, f64>,
    {
        let dims = [12, 12, 12];
        let spec = GlobalSpec::<f64> {
            dims,
            nu: 1.0 / 30.0,
            collision: CollisionKind::Trt {
                magic: CollisionKind::MAGIC_STD,
            },
            periodic: [true, true, true],
            ..Default::default()
        };
        let mut solver: Solver<D3Q27, f64, B, LocalPeriodic> =
            Solver::new(&spec, &[], &[], [1, 1, 1], backend, LocalPeriodic);
        solver.init_with(|x, y, z| {
            let sx = (2.0 * PI * x as f64 / dims[0] as f64).sin();
            let cx = (2.0 * PI * x as f64 / dims[0] as f64).cos();
            let sy = (2.0 * PI * y as f64 / dims[1] as f64).sin();
            let cy = (2.0 * PI * y as f64 / dims[1] as f64).cos();
            let cz = (2.0 * PI * z as f64 / dims[2] as f64).cos();
            let amp = 0.02;
            (1.0, [amp * sx * cy * cz, -amp * cx * sy * cz, 0.0])
        });
        solver.run_guarded(240, 40).unwrap();
        assert!(solver.total_mass().is_finite());
    }

    #[test]
    fn d3q27_periodic_tgv_runs_cpu_scalar_and_simd() {
        run_d3q27_periodic_tgv(CpuScalar::default());
        run_d3q27_periodic_tgv(CpuSimd::default());
    }

    /// Lock the D2Q9 tables to the retired V1 engine's (`lbm_core::lattice`,
    /// deleted 2026-07-05; constants embedded verbatim from its source) —
    /// the direction ordering is project-wide load-bearing and must never
    /// drift.
    #[test]
    fn d2q9_matches_v1_tables() {
        const V1_CX: [i32; 9] = [0, 1, 0, -1, 0, 1, -1, -1, 1];
        const V1_CY: [i32; 9] = [0, 0, 1, 0, -1, 1, 1, -1, -1];
        #[rustfmt::skip]
        const V1_W: [f64; 9] = [
            4.0 / 9.0,
            1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0,
            1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0,
        ];
        const V1_OPP: [usize; 9] = [0, 3, 4, 1, 2, 7, 8, 5, 6];
        const V1_CS2: f64 = 1.0 / 3.0;
        const V1_PAIRS: [(usize, usize); 4] = [(1, 3), (2, 4), (5, 7), (6, 8)];
        assert_eq!(D2Q9::Q, 9);
        assert_eq!(D2Q9::CS2, V1_CS2);
        for q in 0..9 {
            assert_eq!(D2Q9::C[q][0] as i32, V1_CX[q]);
            assert_eq!(D2Q9::C[q][1] as i32, V1_CY[q]);
            assert_eq!(D2Q9::C[q][2], 0);
            assert_eq!(D2Q9::W[q], V1_W[q]);
            assert_eq!(D2Q9::OPP[q], V1_OPP[q]);
        }
        assert_eq!(D2Q9::PAIRS, &V1_PAIRS);
        assert_eq!(D2Q9::REST, 0);
    }

    #[test]
    fn d3q19_opposites_are_adjacent() {
        // The chosen standard ordering pairs opposites as (odd, odd+1).
        for q in (1..19).step_by(2) {
            assert_eq!(D3Q19::OPP[q], q + 1);
            assert_eq!(D3Q19::OPP[q + 1], q);
        }
    }

    #[test]
    fn dir_index_roundtrips() {
        for q in 0..D2Q9::Q {
            assert_eq!(D2Q9::dir_index(D2Q9::C[q]), q);
        }
        for q in 0..D3Q19::Q {
            assert_eq!(D3Q19::dir_index(D3Q19::C[q]), q);
        }
        for q in 0..D3Q27::Q {
            assert_eq!(D3Q27::dir_index(D3Q27::C[q]), q);
        }
    }

    #[test]
    #[should_panic(expected = "not a lattice direction")]
    fn dir_index_rejects_non_directions() {
        D2Q9::dir_index([2, 0, 0]);
    }

    #[test]
    fn face_geometry() {
        assert_eq!(Face::XNeg.n_in(), [1, 0, 0]);
        assert_eq!(Face::YPos.n_in(), [0, -1, 0]);
        assert_eq!(Face::ZPos.n_in(), [0, 0, -1]);
        for f in Face::ALL {
            assert_eq!(f.opposite().opposite(), f);
            assert_eq!(f.opposite().axis(), f.axis());
            assert_ne!(f.opposite().is_neg(), f.is_neg());
            assert_eq!(Face::ALL[f.index()], f);
        }
    }
}
