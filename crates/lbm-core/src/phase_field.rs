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
        if m > 1.0 / 6.0 {
            return Err(PhaseFieldError {
                message: format!(
                    "mobility must be <= 1/6 for the explicit D3Q19 phase update (got {m})"
                ),
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
    let mut values = populations;
    for i in 1..values.len() {
        let key = values[i];
        let mut j = i;
        while j > 0 && values[j - 1] > key {
            values[j] = values[j - 1];
            j -= 1;
        }
        values[j] = key;
    }
    let mut sum = T::zero();
    for v in values {
        sum = sum + v;
    }
    sum
}

#[inline]
pub fn equilibrium<T: Real>(phi: T, u: [T; 3]) -> [T; 19] {
    let three = T::r(3.0);
    let f45 = T::r(4.5);
    let f15 = T::r(1.5);
    let usq = u[0] * u[0] + u[1] * u[1] + u[2] * u[2];
    let mut out = [T::zero(); 19];
    for q in 0..D3Q19::Q {
        let c = D3Q19::C[q];
        let cu = T::r(c[0] as f64) * u[0] + T::r(c[1] as f64) * u[1] + T::r(c[2] as f64) * u[2];
        out[q] = phi * T::r(D3Q19::W[q]) * (T::one() + three * cu + f45 * cu * cu - f15 * usq);
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
    let at = |dx: isize, dy: isize, dz: isize| -> T {
        plane[geom.pidx_i(x as isize + dx, y as isize + dy, z as isize + dz)]
    };
    let phi0 = at(0, 0, 0);
    let sixth = T::one() / T::r(6.0);
    let twelfth = T::one() / T::r(12.0);
    let third = T::one() / T::r(3.0);

    let xp = at(1, 0, 0);
    let xm = at(-1, 0, 0);
    let yp = at(0, 1, 0);
    let ym = at(0, -1, 0);
    let zp = at(0, 0, 1);
    let zm = at(0, 0, -1);
    let xpyp = at(1, 1, 0);
    let xmym = at(-1, -1, 0);
    let xpzp = at(1, 0, 1);
    let xmzm = at(-1, 0, -1);
    let ypzp = at(0, 1, 1);
    let ymzm = at(0, -1, -1);
    let xpym = at(1, -1, 0);
    let xmyp = at(-1, 1, 0);
    let xpzm = at(1, 0, -1);
    let xmzp = at(-1, 0, 1);
    let ypzm = at(0, 1, -1);
    let ymzp = at(0, -1, 1);

    let grad = [
        sixth * (xp - xm)
            + twelfth
                * (sum4_sorted([xpyp, xpym, xpzp, xpzm]) - sum4_sorted([xmyp, xmym, xmzp, xmzm])),
        sixth * (yp - ym)
            + twelfth
                * (sum4_sorted([xpyp, xmyp, ypzp, ypzm]) - sum4_sorted([xpym, xmym, ymzp, ymzm])),
        sixth * (zp - zm)
            + twelfth
                * (sum4_sorted([xpzp, xmzp, ypzp, ymzp]) - sum4_sorted([xpzm, xmzm, ypzm, ymzm])),
    ];
    let lap = third
        * ((xp - phi0) + (xm - phi0) + (yp - phi0) + (ym - phi0) + (zp - phi0) + (zm - phi0))
        + sixth
            * ((xpyp - phi0)
                + (xmym - phi0)
                + (xpzp - phi0)
                + (xmzm - phi0)
                + (ypzp - phi0)
                + (ymzm - phi0)
                + (xpym - phi0)
                + (xmyp - phi0)
                + (xpzm - phi0)
                + (xmzp - phi0)
                + (ypzm - phi0)
                + (ymzp - phi0));
    (grad, lap)
}

#[inline]
fn sum4_sorted<T: Real>(mut values: [T; 4]) -> T {
    for i in 1..values.len() {
        let key = values[i];
        let mut j = i;
        while j > 0 && values[j - 1] > key {
            values[j] = values[j - 1];
            j -= 1;
        }
        values[j] = key;
    }
    (values[0] + values[1]) + (values[2] + values[3])
}

#[inline]
fn normal<T: Real>(grad: [T; 3]) -> [T; 3] {
    let mag = (grad[0] * grad[0] + grad[1] * grad[1] + grad[2] * grad[2]).sqrt();
    let denom = mag + T::epsilon();
    [grad[0] / denom, grad[1] / denom, grad[2] / denom]
}

#[inline]
fn sharpening_vector<T: Real>(params: PhaseFieldParams<T>, phi: T, grad: [T; 3]) -> [T; 3] {
    let n = normal(grad);
    let sharpen = T::r(4.0) / params.interface_width * phi * (T::one() - phi);
    [sharpen * n[0], sharpen * n[1], sharpen * n[2]]
}

pub fn phase_flux_jphi<T: Real>(params: PhaseFieldParams<T>, phi: T, grad: [T; 3]) -> [T; 3] {
    let sharpen = sharpening_vector(params, phi, grad);
    [
        -params.mobility * (grad[0] - sharpen[0]),
        -params.mobility * (grad[1] - sharpen[1]),
        -params.mobility * (grad[2] - sharpen[2]),
    ]
}

pub fn collide_source<T: Real>(params: PhaseFieldParams<T>, phi: T, grad: [T; 3], q: usize) -> T {
    let sharpen = sharpening_vector(params, phi, grad);
    let sixth = T::one() / T::r(6.0);
    let twelfth = T::one() / T::r(12.0);
    match q {
        0 => T::zero(),
        1 => sixth * sharpen[0],
        2 => -sixth * sharpen[0],
        3 => sixth * sharpen[1],
        4 => -sixth * sharpen[1],
        5 => sixth * sharpen[2],
        6 => -sixth * sharpen[2],
        7 => twelfth * (sharpen[0] + sharpen[1]),
        8 => -twelfth * (sharpen[0] + sharpen[1]),
        9 => twelfth * (sharpen[0] + sharpen[2]),
        10 => -twelfth * (sharpen[0] + sharpen[2]),
        11 => twelfth * (sharpen[1] + sharpen[2]),
        12 => -twelfth * (sharpen[1] + sharpen[2]),
        13 => twelfth * (sharpen[0] - sharpen[1]),
        14 => -twelfth * (sharpen[0] - sharpen[1]),
        15 => twelfth * (sharpen[0] - sharpen[2]),
        16 => -twelfth * (sharpen[0] - sharpen[2]),
        17 => twelfth * (sharpen[1] - sharpen[2]),
        18 => -twelfth * (sharpen[1] - sharpen[2]),
        _ => unreachable!("D3Q19 source direction index out of range"),
    }
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn uniform_phi_remains_uniform_no_advection() {
        let params = PhaseFieldParams::new(4.0, 0.1);
        let j = phase_flux_jphi(params, 0.42, [0.0, 0.0, 0.0]);
        assert_eq!(j, [-0.0, -0.0, -0.0]);
    }

    #[test]
    fn rejects_invalid_interface_width() {
        assert!(PhaseFieldParams::new(-1.0, 0.1).validate().is_err());
        assert!(PhaseFieldParams::new(2.99, 0.1).validate().is_err());
    }

    #[test]
    fn rejects_invalid_mobility() {
        assert!(PhaseFieldParams::new(4.0, 0.0).validate().is_err());
        assert!(PhaseFieldParams::new(4.0, -0.1).validate().is_err());
    }

    #[test]
    fn clipping_diagnostics_reported_after_forcing_out_of_bounds() {
        let policy = ClippingPolicy::ClipToBounds { min: 0.0, max: 1.0 };
        assert_eq!(apply_clipping(policy, -0.25), (0.0, true));
        assert_eq!(apply_clipping(policy, 1.25), (1.0, true));
        assert_eq!(apply_clipping(policy, 0.5), (0.5, false));
    }
}
