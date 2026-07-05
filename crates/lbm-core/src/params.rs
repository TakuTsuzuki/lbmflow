//! Step parameters shared by every backend.

use crate::lattice::{Lattice, Q_MAX};
use crate::real::Real;

/// Maximum prescribed speed (lattice units) a configuration may request
/// before validation rejects it — the low-Mach limit shared by the V2
/// `GlobalSpec::validate` and the compat facade (`compat::domain::MAX_SPEED`
/// re-exports this so there is a single source of truth).
pub const MAX_SPEED: f64 = 0.3;

/// Collision operator selection (identical semantics to V1
/// `lbm_core::domain::Collision`).
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum CollisionKind {
    /// Single-relaxation-time BGK.
    Bgk,
    /// Two-relaxation-time with the given magic parameter Λ.
    Trt {
        /// Magic parameter Λ = (1/ω+ − 1/2)(1/ω− − 1/2).
        magic: f64,
    },
}

impl CollisionKind {
    /// The standard "magic" value 3/16 (exact half-way walls for Poiseuille).
    pub const MAGIC_STD: f64 = 3.0 / 16.0;

    /// Derive `(omega_p, omega_m)` from the viscosity, exactly as V1 does.
    pub fn omegas(self, nu: f64) -> (f64, f64) {
        let tau = 3.0 * nu + 0.5;
        let omega_p = 1.0 / tau;
        let omega_m = match self {
            CollisionKind::Bgk => omega_p,
            CollisionKind::Trt { magic } => {
                let lam_p = tau - 0.5;
                1.0 / (magic / lam_p + 0.5)
            }
        };
        (omega_p, omega_m)
    }
}

/// Open boundary condition on one global face. Walls and periodic wraps are
/// *not* face BCs here: walls are one-cell solid rims (data in the fields)
/// and periodicity is a halo-exchange concern, so both map to `Closed`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum FaceBC<T: Real> {
    /// No open-BC pass on this face (periodic axis or wall rim).
    Closed,
    /// Zou–He velocity boundary (prescribed `u`).
    Velocity {
        /// Prescribed velocity (z component ignored for 2D lattices).
        u: [T; 3],
    },
    /// Zou–He pressure boundary (prescribed density, `p = cs^2 rho`).
    Pressure {
        /// Prescribed density.
        rho: T,
    },
    /// Zero-gradient outflow (copies unknowns from one cell inward).
    Outflow,
    /// Convective (radiation) outflow with mass pinning (V1 semantics).
    Convective {
        /// Advection speed of the outgoing characteristics.
        u_conv: T,
    },
}

impl<T: Real> FaceBC<T> {
    /// Whether this face runs an open-BC pass after streaming.
    #[inline]
    pub fn is_open(&self) -> bool {
        !matches!(self, FaceBC::Closed)
    }
}

/// Per-step scalar parameters (uniform over the grid). Relaxation rates are
/// kept in `f64` and converted to `T` when kernel constants are built, so a
/// step sees exactly the values V1 computes.
#[derive(Clone, Copy, Debug)]
pub struct StepParams<T: Real> {
    /// Symmetric relaxation rate `1/tau`.
    pub omega_p: f64,
    /// Antisymmetric relaxation rate (TRT); equals `omega_p` for BGK.
    pub omega_m: f64,
    /// Uniform body force (Guo forcing).
    pub force: [T; 3],
    /// Open BC per global face, `Face::index()` order.
    pub faces: [FaceBC<T>; 6],
}

/// Kernel constants in working precision, rebuilt from [`StepParams`] each
/// step exactly like V1 `Simulation::params()`.
#[derive(Clone, Copy)]
pub struct KParams<T: Real> {
    /// `T::r(omega_p)`.
    pub omega_p: T,
    /// `T::r(omega_m)`.
    pub omega_m: T,
    /// Guo prefactor `1 - omega_p/2` (computed in f64, then converted).
    pub cp: T,
    /// Guo prefactor `1 - omega_m/2`.
    pub cm: T,
    /// Uniform body force.
    pub force: [T; 3],
    /// Discrete velocities as `T` (first `Q` entries valid).
    pub cr: [[T; 3]; Q_MAX],
    /// Weights as `T` (first `Q` entries valid).
    pub wr: [T; Q_MAX],
}

impl<T: Real> KParams<T> {
    /// Build kernel constants for lattice `L` (V1 `params()` equivalent).
    pub fn new<L: Lattice>(p: &StepParams<T>) -> Self {
        let mut cr = [[T::zero(); 3]; Q_MAX];
        let mut wr = [T::zero(); Q_MAX];
        for q in 0..L::Q {
            for a in 0..3 {
                cr[q][a] = T::r(L::C[q][a] as f64);
            }
            wr[q] = T::r(L::W[q]);
        }
        Self {
            omega_p: T::r(p.omega_p),
            omega_m: T::r(p.omega_m),
            cp: T::r(1.0 - p.omega_p / 2.0),
            cm: T::r(1.0 - p.omega_m / 2.0),
            force: p.force,
            cr,
            wr,
        }
    }
}

/// Backend-side reduction kinds (accumulated in `f64` regardless of `T`,
/// V1 diagnostic convention).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Reduction {
    /// Number of fluid (non-solid) core cells.
    FluidCells,
    /// `sum f_dev` over fluid cells (physical mass = fluid_cells + this).
    MassDeviation,
    /// `sum (m_a + f_a/2)` over fluid cells: physical momentum component,
    /// including the Guo half-force correction. Axis 0/1/2.
    Momentum(usize),
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn omegas_match_v1_derivation() {
        let nu = 0.02;
        let tau = 3.0 * nu + 0.5;
        let (op, om) = CollisionKind::Trt { magic: 3.0 / 16.0 }.omegas(nu);
        assert_eq!(op, 1.0 / tau);
        assert_eq!(om, 1.0 / ((3.0 / 16.0) / (tau - 0.5) + 0.5));
        let (bp, bm) = CollisionKind::Bgk.omegas(nu);
        assert_eq!(bp, bm);
    }
}
