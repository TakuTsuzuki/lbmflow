//! Step parameters shared by every backend.

use crate::lattice::{Lattice, Q_MAX};
use crate::real::Real;

/// Maximum prescribed speed (lattice units) a configuration may request
/// before validation rejects it — the low-Mach limit shared by the V2
/// `GlobalSpec::validate` and the compat facade (`compat::domain::MAX_SPEED`
/// re-exports this so there is a single source of truth).
pub const MAX_SPEED: f64 = 0.3;

/// Compile-time ablation switch for the unresolved central-moment
/// velocity-dependent shear-rate modifier. Default is off, meaning normal
/// builds keep the pending `-0.16 |u|^2` term; set true only for the
/// ANOM-P4-008 E1 ablation rerun.
pub const CENTRAL_MOMENT_DISABLE_VELOCITY_CORRECTION_FOR_ABLATION: bool = false;

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
    /// Cascaded central-moment collision. The shear-rate field controls the
    /// second-order deviatoric central moments, the second-order trace and
    /// all higher-order central moments relax to equilibrium at rate 1.0.
    CentralMoment {
        /// Relaxation rate for second-order shear central moments.
        omega_shear: f64,
    },
}

impl Default for CollisionKind {
    fn default() -> Self {
        CollisionKind::Trt {
            magic: CollisionKind::MAGIC_STD,
        }
    }
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
            CollisionKind::CentralMoment { .. } => omega_p,
        };
        (omega_p, omega_m)
    }

    /// Uniform shear relaxation rate used by the central-moment branch.
    pub fn omega_shear(self, nu: f64) -> f64 {
        match self {
            CollisionKind::CentralMoment { omega_shear } => omega_shear,
            CollisionKind::Bgk | CollisionKind::Trt { .. } => self.omegas(nu).0,
        }
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

/// Inclusive global-cell box for a localized volume source/sink.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct SourceRegion {
    pub lo: [usize; 3],
    pub hi: [usize; 3],
}

/// Per-step mass source over a [`SourceRegion`].
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum SourceKind<T: Real> {
    /// Total mass per lattice step, distributed uniformly over the region.
    /// Negative values are sinks.
    MassFlow { q_lu: T },
    /// Total mass per lattice step carrying prescribed velocity `u`.
    Jet { q_lu: T, u: [T; 3] },
}

/// Localized interior volume source/sink.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct VolumeSource<T: Real> {
    pub region: SourceRegion,
    pub kind: SourceKind<T>,
}

/// Rectangular boundary-condition override on one global face.
///
/// `lo`/`hi` are inclusive in-face coordinates on the two remaining axes in
/// ascending axis order, matching [`crate::lattice::Face::tangents`].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FacePatch<T: Real> {
    pub face: usize,
    pub lo: [usize; 2],
    pub hi: [usize; 2],
    pub bc: FaceBC<T>,
}

/// Per-step scalar parameters (uniform over the grid). Relaxation rates are
/// kept in `f64` and converted to `T` when kernel constants are built, so a
/// step sees exactly the values V1 computes.
#[derive(Clone, Debug)]
pub struct StepParams<T: Real> {
    /// Collision operator selected by the validated global spec.
    pub collision: CollisionKind,
    /// Symmetric relaxation rate `1/tau`.
    pub omega_p: f64,
    /// Antisymmetric relaxation rate (TRT); equals `omega_p` for BGK.
    pub omega_m: f64,
    /// Uniform body force (Guo forcing).
    pub force: [T; 3],
    /// Open BC per global face, `Face::index()` order.
    pub faces: [FaceBC<T>; 6],
    /// Localized interior volume sources/sinks.
    pub sources: Vec<VolumeSource<T>>,
    /// Per-cell open-BC patches on global faces.
    pub face_patches: Vec<FacePatch<T>>,
}

/// Kernel constants in working precision, rebuilt from [`StepParams`] each
/// step exactly like V1 `Simulation::params()`.
#[derive(Clone, Copy)]
pub struct KParams<T: Real> {
    /// Whether this step uses the central-moment branch.
    pub central_moment: bool,
    /// `T::r(omega_p)`.
    pub omega_p: T,
    /// `T::r(omega_m)`.
    pub omega_m: T,
    /// Central-moment shear relaxation rate.
    pub omega_shear: T,
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
            central_moment: matches!(p.collision, CollisionKind::CentralMoment { .. }),
            omega_p: T::r(p.omega_p),
            omega_m: T::r(p.omega_m),
            omega_shear: T::r(match p.collision {
                CollisionKind::CentralMoment { omega_shear } => omega_shear,
                CollisionKind::Bgk | CollisionKind::Trt { .. } => p.omega_p,
            }),
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
