//! D2Q9 lattice constants in V1 (`lbm_core::lattice`) shape.
//!
//! Direction ordering (the single source of truth for the whole project):
//!
//! ```text
//!   6  2  5
//!   3  0  1
//!   7  4  8
//! ```
//!
//! Every table is derived from [`crate::lattice::D2Q9`] by `const fn`, so
//! the facade can never drift from the V2 core tables (which are themselves
//! locked to V1 by test).

use crate::lattice::{D2Q9, Lattice as _};

/// Number of discrete velocities.
pub const Q: usize = 9;

const fn c_component(a: usize) -> [i32; Q] {
    let mut t = [0i32; Q];
    let mut q = 0;
    while q < Q {
        t[q] = D2Q9::C[q][a] as i32;
        q += 1;
    }
    t
}

const fn w_table() -> [f64; Q] {
    let mut t = [0.0f64; Q];
    let mut q = 0;
    while q < Q {
        t[q] = D2Q9::W[q];
        q += 1;
    }
    t
}

const fn opp_table() -> [usize; Q] {
    let mut t = [0usize; Q];
    let mut q = 0;
    while q < Q {
        t[q] = D2Q9::OPP[q];
        q += 1;
    }
    t
}

const fn pairs_table() -> [(usize, usize); 4] {
    let mut t = [(0usize, 0usize); 4];
    let mut i = 0;
    while i < 4 {
        t[i] = D2Q9::PAIRS[i];
        i += 1;
    }
    t
}

/// x-components of the discrete velocities.
pub const CX: [i32; Q] = c_component(0);
/// y-components of the discrete velocities.
pub const CY: [i32; Q] = c_component(1);

/// Quadrature weights.
pub const W: [f64; Q] = w_table();

/// Index of the opposite direction: `C[OPP[q]] == -C[q]`.
pub const OPP: [usize; Q] = opp_table();

/// Squared lattice speed of sound, `cs^2 = 1/3`.
pub const CS2: f64 = D2Q9::CS2;

/// TRT direction pairs `(q, OPP[q])` excluding the rest population.
pub const PAIRS: [(usize, usize); 4] = pairs_table();

/// Find the direction index with the given velocity components.
///
/// Panics if `(cx, cy)` is not a D2Q9 velocity.
pub fn dir_index(cx: i32, cy: i32) -> usize {
    (0..Q)
        .find(|&q| CX[q] == cx && CY[q] == cy)
        .unwrap_or_else(|| panic!("({cx},{cy}) is not a D2Q9 direction"))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Retired-V1 (`lbm_core::lattice`, deleted 2026-07-05) constants,
    /// embedded verbatim from its source.
    #[test]
    fn facade_tables_match_v1() {
        #[rustfmt::skip]
        const V1_W: [f64; 9] = [
            4.0 / 9.0,
            1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0,
            1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0,
        ];
        assert_eq!(Q, 9);
        assert_eq!(CX, [0, 1, 0, -1, 0, 1, -1, -1, 1]);
        assert_eq!(CY, [0, 0, 1, 0, -1, 1, 1, -1, -1]);
        assert_eq!(W, V1_W);
        assert_eq!(OPP, [0, 3, 4, 1, 2, 7, 8, 5, 6]);
        assert_eq!(CS2, 1.0 / 3.0);
        assert_eq!(PAIRS, [(1, 3), (2, 4), (5, 7), (6, 8)]);
        for q in 0..Q {
            assert_eq!(dir_index(CX[q], CY[q]), q);
        }
    }
}
