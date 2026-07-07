//! Passive scalar advection-diffusion LBE helpers.

use crate::lattice::Lattice;
use crate::real::Real;

pub fn scalar_equilibrium<L: Lattice, T: Real>(concentration: T, u: [T; 3]) -> Vec<T> {
    let mut out = vec![T::zero(); L::Q];
    let three = T::r(3.0);
    for q in 0..L::Q {
        let c = L::C[q];
        let cu = T::r(c[0] as f64) * u[0] + T::r(c[1] as f64) * u[1] + T::r(c[2] as f64) * u[2];
        out[q] = T::r(L::W[q]) * concentration * (T::one() + three * cu);
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::lattice::{D2Q9, D3Q19};

    #[test]
    fn scalar_equilibrium_conserves_concentration() {
        let h = scalar_equilibrium::<D3Q19, f64>(2.5, [0.02, -0.01, 0.0]);
        let sum: f64 = h.iter().sum();
        assert!((sum - 2.5).abs() < 1.0e-12);
    }

    #[test]
    fn scalar_equilibrium_first_moment_matches_advection_velocity() {
        let u = [0.03, -0.02, 0.0];
        let h = scalar_equilibrium::<D2Q9, f64>(4.0, u);
        let mut mom = [0.0; 3];
        for q in 0..D2Q9::Q {
            mom[0] += D2Q9::C[q][0] as f64 * h[q];
            mom[1] += D2Q9::C[q][1] as f64 * h[q];
        }
        assert!((mom[0] - 4.0 * u[0]).abs() < 1.0e-12);
        assert!((mom[1] - 4.0 * u[1]).abs() < 1.0e-12);
    }
}
