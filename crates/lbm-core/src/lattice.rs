//! D2Q9 lattice constants.
//!
//! Direction ordering (the single source of truth for the whole project):
//!
//! ```text
//!   6  2  5
//!   3  0  1
//!   7  4  8
//! ```

/// Number of discrete velocities.
pub const Q: usize = 9;

/// x-components of the discrete velocities.
pub const CX: [i32; Q] = [0, 1, 0, -1, 0, 1, -1, -1, 1];
/// y-components of the discrete velocities.
pub const CY: [i32; Q] = [0, 0, 1, 0, -1, 1, 1, -1, -1];

/// Quadrature weights.
pub const W: [f64; Q] = [
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

/// Index of the opposite direction: `C[OPP[q]] == -C[q]`.
pub const OPP: [usize; Q] = [0, 3, 4, 1, 2, 7, 8, 5, 6];

/// Squared lattice speed of sound, `cs^2 = 1/3`.
pub const CS2: f64 = 1.0 / 3.0;

/// TRT direction pairs `(q, OPP[q])` excluding the rest population.
pub const PAIRS: [(usize, usize); 4] = [(1, 3), (2, 4), (5, 7), (6, 8)];

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

    #[test]
    fn opposites_are_consistent() {
        for q in 0..Q {
            assert_eq!(CX[OPP[q]], -CX[q]);
            assert_eq!(CY[OPP[q]], -CY[q]);
            assert_eq!(OPP[OPP[q]], q);
            assert_eq!(W[OPP[q]], W[q]);
        }
    }

    #[test]
    fn weights_reproduce_lattice_moments() {
        // sum w = 1, sum w c = 0, sum w c c = cs^2 I
        let sw: f64 = W.iter().sum();
        assert!((sw - 1.0).abs() < 1e-15);
        let (mut mx, mut my) = (0.0, 0.0);
        let (mut mxx, mut mxy, mut myy) = (0.0, 0.0, 0.0);
        for q in 0..Q {
            mx += W[q] * CX[q] as f64;
            my += W[q] * CY[q] as f64;
            mxx += W[q] * (CX[q] * CX[q]) as f64;
            mxy += W[q] * (CX[q] * CY[q]) as f64;
            myy += W[q] * (CY[q] * CY[q]) as f64;
        }
        assert!(mx.abs() < 1e-15 && my.abs() < 1e-15);
        assert!((mxx - CS2).abs() < 1e-15);
        assert!(mxy.abs() < 1e-15);
        assert!((myy - CS2).abs() < 1e-15);
    }

    #[test]
    fn dir_index_roundtrips() {
        for q in 0..Q {
            assert_eq!(dir_index(CX[q], CY[q]), q);
        }
    }
}
