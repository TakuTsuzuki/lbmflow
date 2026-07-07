//! Conservative Allen-Cahn phase-field transport for the W-VOF O1 stage.
//!
//! This module implements only the prescribed-velocity phase-field LBE:
//! no density feedback, surface-tension force, gravity edit, wetting boundary,
//! or hydrodynamic momentum coupling is active here.

use crate::fields::LocalGeom;
use crate::lattice::{Lattice, D3Q19};
use crate::real::Real;

/// Conservative Allen-Cahn phase-field parameters in lattice units.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PhaseFieldParams<T: Real> {
    /// Diffuse-interface width `W` in lattice cells.
    pub interface_width: T,
    /// Mobility `M`; `tau_phi = 3M + 0.5`.
    pub mobility: T,
}

/// Phase-field construction/validation error.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PhaseFieldError {
    pub message: String,
}

impl PhaseFieldError {
    fn new(message: impl Into<String>) -> Self {
        Self {
            message: message.into(),
        }
    }
}

impl std::fmt::Display for PhaseFieldError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.message)
    }
}

impl std::error::Error for PhaseFieldError {}

/// Lightweight behavior diagnostics from a completed phase-field step.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct PhaseFieldDiagnostics {
    /// Sum of `phi` over fluid core cells.
    pub total_phi: f64,
    /// Minimum `phi` over fluid core cells.
    pub min_phi: f64,
    /// Maximum `phi` over fluid core cells.
    pub max_phi: f64,
}

impl<T: Real> PhaseFieldParams<T> {
    /// Validate the frozen O1 validity domain:
    /// `W in [4,5]`, `M in (0,1/6]`, `tau_phi = 3M + 0.5`.
    pub fn validate(self) -> Result<Self, PhaseFieldError> {
        if self.interface_width < T::r(4.0) || self.interface_width > T::r(5.0) {
            return Err(PhaseFieldError::new(format!(
                "phase-field interface width W={} is outside the W-VOF O1 validity domain [4,5]",
                self.interface_width
            )));
        }
        if self.mobility <= T::zero() || self.mobility > T::r(1.0 / 6.0) {
            return Err(PhaseFieldError::new(format!(
                "phase-field mobility M={} is outside the W-VOF O1 validity domain (0,1/6]",
                self.mobility
            )));
        }
        Ok(self)
    }

    #[inline]
    pub(crate) fn omega(self) -> T {
        T::one() / (T::r(3.0) * self.mobility + T::r(0.5))
    }
}

#[inline]
pub(crate) fn equilibrium<T: Real>(phi: T, u: [T; 3]) -> [T; 19] {
    let three = T::r(3.0);
    let f45 = T::r(4.5);
    let f15 = T::r(1.5);
    let usq = u[0] * u[0] + u[1] * u[1] + u[2] * u[2];
    let mut geq = [T::zero(); 19];
    for q in 0..D3Q19::Q {
        let c = D3Q19::C[q];
        let cu = T::r(c[0] as f64) * u[0] + T::r(c[1] as f64) * u[1] + T::r(c[2] as f64) * u[2];
        geq[q] = phi * T::r(D3Q19::W[q]) * (T::one() + three * cu + f45 * cu * cu - f15 * usq);
    }
    geq
}

#[inline]
pub(crate) fn grad_lap<T: Real>(
    geom: LocalGeom,
    phi_plane: &[T],
    x: usize,
    y: usize,
    z: usize,
) -> ([T; 3], T) {
    let p0 = geom.pidx(x, y, z);
    let phi0 = phi_plane[p0];
    let inv_cs2 = T::r(3.0);
    let two_inv_cs2 = T::r(6.0);
    let mut grad = [T::zero(); 3];
    let mut lap = T::zero();
    for q in 1..D3Q19::Q {
        let c = D3Q19::C[q];
        let pj = geom.pidx_i(
            x as isize + c[0] as isize,
            y as isize + c[1] as isize,
            z as isize + c[2] as isize,
        );
        let w = T::r(D3Q19::W[q]);
        let phij = phi_plane[pj];
        grad[0] = grad[0] + inv_cs2 * w * T::r(c[0] as f64) * phij;
        grad[1] = grad[1] + inv_cs2 * w * T::r(c[1] as f64) * phij;
        grad[2] = grad[2] + inv_cs2 * w * T::r(c[2] as f64) * phij;
        lap = lap + two_inv_cs2 * w * (phij - phi0);
    }
    (grad, lap)
}

#[inline]
pub(crate) fn normal<T: Real>(grad: [T; 3]) -> [T; 3] {
    let mag = (grad[0] * grad[0] + grad[1] * grad[1] + grad[2] * grad[2]).sqrt();
    let denom = mag + T::epsilon();
    [grad[0] / denom, grad[1] / denom, grad[2] / denom]
}

/// Single discrete interface flux source path:
/// `J_phi = -M [grad(phi) - (4/W) phi(1-phi) n_hat]`.
#[inline]
fn sharpening_vector<T: Real>(params: PhaseFieldParams<T>, phi: T, grad: [T; 3]) -> [T; 3] {
    let n = normal(grad);
    let sharpen = T::r(4.0) / params.interface_width * phi * (T::one() - phi);
    [sharpen * n[0], sharpen * n[1], sharpen * n[2]]
}

#[inline]
pub(crate) fn phase_flux_jphi<T: Real>(
    params: PhaseFieldParams<T>,
    phi: T,
    grad: [T; 3],
) -> [T; 3] {
    let sharpen = sharpening_vector(params, phi, grad);
    [
        -params.mobility * (grad[0] - sharpen[0]),
        -params.mobility * (grad[1] - sharpen[1]),
        -params.mobility * (grad[2] - sharpen[2]),
    ]
}

#[inline]
pub(crate) fn collide_source<T: Real>(
    params: PhaseFieldParams<T>,
    phi: T,
    grad: [T; 3],
    q: usize,
) -> T {
    let sharpen = sharpening_vector(params, phi, grad);
    let c = D3Q19::C[q];
    let dot = T::r(c[0] as f64) * sharpen[0]
        + T::r(c[1] as f64) * sharpen[1]
        + T::r(c[2] as f64) * sharpen[2];
    T::r(D3Q19::W[q]) * T::r(3.0) * dot
}
