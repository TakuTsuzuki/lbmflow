//! Conservative Allen-Cahn phase-field helpers.

use crate::fields::LocalGeom;
use crate::lattice::{Lattice, D3Q19};
use crate::real::Real;
use serde::{Deserialize, Serialize};

/// Explicit boundedness policy for transported phase fraction.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum ClippingPolicy<T: Real> {
    Off,
    ClipToBounds { min: T, max: T },
}

impl<T: Real> Default for ClippingPolicy<T> {
    fn default() -> Self {
        Self::Off
    }
}

/// Conservative Allen-Cahn parameters in lattice units.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PhaseFieldParams<T: Real> {
    pub interface_width: T,
    pub mobility: T,
    pub clipping_policy: ClippingPolicy<T>,
}

impl<T: Real> PhaseFieldParams<T> {
    pub fn new(interface_width: T, mobility: T) -> Self {
        Self {
            interface_width,
            mobility,
            clipping_policy: ClippingPolicy::Off,
        }
    }

    pub fn validate(self) -> Result<Self, PhaseFieldError> {
        let w = self.interface_width.as_f64();
        let m = self.mobility.as_f64();
        if !w.is_finite() || w <= 0.0 {
            return Err(PhaseFieldError {
                message: format!("interface_width must be finite and > 0 (got {w})"),
            });
        }
        if w < 3.0 {
            return Err(PhaseFieldError {
                message: format!(
                    "interface_width/dx must be >= 3 for a resolved diffuse interface (got {w})"
                ),
            });
        }
        if !m.is_finite() || m <= 0.0 {
            return Err(PhaseFieldError {
                message: format!("mobility must be finite and > 0 (got {m})"),
            });
        }
        if let ClippingPolicy::ClipToBounds { min, max } = self.clipping_policy {
            let (lo, hi) = (min.as_f64(), max.as_f64());
            if !lo.is_finite() || !hi.is_finite() || lo >= hi {
                return Err(PhaseFieldError {
                    message: format!("clipping bounds must be finite and ordered (got {lo}..{hi})"),
                });
            }
        }
        Ok(self)
    }

    pub fn omega(self) -> T {
        let tau = T::r(0.5) + self.mobility / (T::r(D3Q19::CS2));
        T::one() / tau
    }
}

impl<T: Real> Default for PhaseFieldParams<T> {
    fn default() -> Self {
        Self::new(T::r(4.0), T::r(0.1))
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct PhaseFieldError {
    pub message: String,
}

impl std::fmt::Display for PhaseFieldError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.message)
    }
}

impl std::error::Error for PhaseFieldError {}

#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct PhaseFieldDiagnostics {
    pub total_phi: f64,
    pub min_phi: f64,
    pub max_phi: f64,
    pub clipped_fraction: f64,
    pub interface_cells: usize,
    pub total_cells: usize,
}

impl PhaseFieldDiagnostics {
    pub fn empty() -> Self {
        Self {
            min_phi: 0.0,
            max_phi: 0.0,
            ..Self::default()
        }
    }
}

#[inline]
pub fn sum_populations<T: Real>(populations: [T; 19]) -> T {
    populations.into_iter().sum()
}

#[inline]
pub fn equilibrium<T: Real>(phi: T, u: [T; 3]) -> [T; 19] {
    let mut out = [T::zero(); 19];
    for q in 0..D3Q19::Q {
        let c = D3Q19::C[q];
        let cu = T::r(c[0] as f64) * u[0] + T::r(c[1] as f64) * u[1] + T::r(c[2] as f64) * u[2];
        out[q] = T::r(D3Q19::W[q]) * phi * (T::one() + cu / T::r(D3Q19::CS2));
    }
    out
}

pub fn grad_lap<T: Real>(
    geom: LocalGeom,
    plane: &[T],
    x: usize,
    y: usize,
    z: usize,
) -> ([T; 3], T) {
    let center = plane[geom.pidx(x, y, z)];
    let half = T::r(0.5);
    let mut grad = [T::zero(); 3];
    let mut lap = T::zero();
    for a in 0..geom.d {
        let mut plus = [0isize; 3];
        let mut minus = [0isize; 3];
        plus[a] = 1;
        minus[a] = -1;
        let p = plane[geom.pidx_i(
            x as isize + plus[0],
            y as isize + plus[1],
            z as isize + plus[2],
        )];
        let m = plane[geom.pidx_i(
            x as isize + minus[0],
            y as isize + minus[1],
            z as isize + minus[2],
        )];
        grad[a] = (p - m) * half;
        lap = lap + p - T::r(2.0) * center + m;
    }
    (grad, lap)
}

pub fn phase_flux_jphi<T: Real>(params: PhaseFieldParams<T>, phi: T, grad: [T; 3]) -> [T; 3] {
    let grad_sq = grad[0] * grad[0] + grad[1] * grad[1] + grad[2] * grad[2];
    let mag = grad_sq.sqrt();
    let eps = T::r(1.0e-12);
    let inv = T::one() / (mag + eps);
    let sharpen = T::r(4.0) / params.interface_width * phi * (T::one() - phi);
    [
        -params.mobility * (grad[0] - sharpen * grad[0] * inv),
        -params.mobility * (grad[1] - sharpen * grad[1] * inv),
        -params.mobility * (grad[2] - sharpen * grad[2] * inv),
    ]
}

pub fn collide_source<T: Real>(params: PhaseFieldParams<T>, phi: T, grad: [T; 3], q: usize) -> T {
    let j = phase_flux_jphi(params, phi, grad);
    let c = D3Q19::C[q];
    T::r(D3Q19::W[q] / D3Q19::CS2)
        * (T::r(c[0] as f64) * j[0] + T::r(c[1] as f64) * j[1] + T::r(c[2] as f64) * j[2])
}

pub fn apply_clipping<T: Real>(policy: ClippingPolicy<T>, value: T) -> (T, bool) {
    match policy {
        ClippingPolicy::Off => (value, false),
        ClippingPolicy::ClipToBounds { min, max } => {
            if value < min {
                (min, true)
            } else if value > max {
                (max, true)
            } else {
                (value, false)
            }
        }
    }
}
