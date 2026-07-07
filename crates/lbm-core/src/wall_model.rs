//! Wall-treatment diagnostics shared by the W1 observable and later closures.
//!
//! W1 is read-only instrumentation: these routines compute wall metrics from
//! the current velocity/mask state and do not write populations or relaxation
//! fields.

use crate::real::Real;

/// Von Karman constant used by the equilibrium log law.
///
/// Provenance: standard smooth-wall law of the wall, as frozen in
/// `docs/proposals/LES_WALL_TREATMENT_SPEC.md` from the Malaspinas-Sagaut
/// LBM wall-model class.
pub const LOG_LAW_KAPPA: f64 = 0.41;

/// Additive smooth-wall log-law constant paired with [`LOG_LAW_KAPPA`].
///
/// Provenance: standard smooth-wall law of the wall,
/// `u+ = ln(y+) / kappa + B`, frozen in the wall-treatment W1/W2 spec.
pub const LOG_LAW_B: f64 = 5.2;

/// Viscous/log-law branch switch y+.
///
/// The value is the matched-layer intersection used by the frozen W1 spec:
/// below this point the diagnostic reports the viscous branch
/// `u_tau = sqrt(nu * u_parallel / y_w)`.
pub const WALL_LAW_SWITCH_Y_PLUS: f64 = 11.6;

/// Source geometry used for a wall-cell metric.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum WallMetricSource {
    /// Half-way bounce-back against the one-cell solid rim.
    HalfwayRim,
    /// Bouzidi interpolated wall link; distance comes from `qd * |c_q|`.
    Bouzidi,
}

/// Read-only wall diagnostic for one wall-adjacent fluid cell.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WallCellMetric<T: Real> {
    /// Global compact index: `(z * ny + y) * nx + x`.
    pub cell_index: usize,
    /// Wall-normal distance from the fluid node to the physical wall.
    pub y_w: T,
    /// Fluid tangential speed relative to the local wall velocity.
    pub u_parallel: T,
    /// Friction velocity from the W1 wall-law diagnostic.
    pub u_tau: T,
    /// Wall coordinate `y+ = y_w * u_tau / nu`.
    pub y_plus: T,
    /// Wall shear stress per unit density: `tau_w / rho = u_tau^2`.
    pub tau_w: T,
    /// Geometry source used for `y_w`.
    pub source: WallMetricSource,
}

/// Friction velocity from the W1 diagnostic formulation.
///
/// The turbulent branch solves
/// `u_parallel / u_tau = ln(y_w * u_tau / nu) / kappa + B` by Newton iteration
/// in `s = ln(y+)`, which enforces positive `y+` without a numerical clamp.
/// If the resulting `y+` is below the matched-layer switch, the viscous branch
/// `u_tau = sqrt(nu * u_parallel / y_w)` is returned.
pub fn friction_velocity<T: Real>(u_parallel: T, y_w: T, nu: T) -> T {
    assert!(y_w > T::zero(), "wall distance must be positive");
    assert!(nu > T::zero(), "viscosity must be positive");
    assert!(
        u_parallel >= T::zero(),
        "tangential speed must be non-negative"
    );
    if u_parallel == T::zero() {
        return T::zero();
    }

    let kappa = T::r(LOG_LAW_KAPPA);
    let b = T::r(LOG_LAW_B);
    let switch = T::r(WALL_LAW_SWITCH_Y_PLUS);
    let u_y_over_nu = u_parallel * y_w / nu;
    let mut s = switch.ln();
    let tol = T::epsilon().sqrt();
    // 32 Newton iterations is a numerical solve allowance, not a model
    // constant; it gives wide headroom for f32/f64 convergence from y+=11.6.
    for _ in 0..32 {
        let y_plus = s.exp();
        let law = s / kappa + b;
        let f = y_plus * law - u_y_over_nu;
        if f.abs() <= tol * (u_y_over_nu.abs() + T::one()) {
            return friction_velocity_from_y_plus(u_parallel, y_w, nu, y_plus, switch);
        }
        let df = y_plus * (law + T::one() / kappa);
        assert!(df != T::zero(), "log-law Newton derivative vanished");
        s = s - f / df;
    }

    let y_plus = s.exp();
    let residual = (y_plus * (s / kappa + b) - u_y_over_nu).abs();
    assert!(
        residual <= tol * (u_y_over_nu.abs() + T::one()),
        "log-law Newton solve failed to converge: y_plus={y_plus}, residual={residual}"
    );
    friction_velocity_from_y_plus(u_parallel, y_w, nu, y_plus, switch)
}

fn friction_velocity_from_y_plus<T: Real>(u_parallel: T, y_w: T, nu: T, y_plus: T, switch: T) -> T {
    if y_plus < switch {
        (nu * u_parallel / y_w).sqrt()
    } else {
        y_plus * nu / y_w
    }
}
