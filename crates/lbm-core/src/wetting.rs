//! Static contact-angle boundary helpers for phase-field runs.

use crate::real::Real;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ContactAngleParams<T: Real> {
    pub theta_deg: T,
}

impl<T: Real> ContactAngleParams<T> {
    pub fn validate(self) -> Result<Self, WettingError> {
        let theta = self.theta_deg.as_f64();
        if !theta.is_finite() || theta <= 0.0 || theta >= 180.0 {
            return Err(WettingError {
                message: format!("static contact angle must be in (0, 180) degrees (got {theta})"),
            });
        }
        Ok(self)
    }

    pub fn cot_theta(self) -> T {
        T::one() / T::r(self.theta_deg.as_f64().to_radians().tan())
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct WettingError {
    pub message: String,
}

impl std::fmt::Display for WettingError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.message)
    }
}

impl std::error::Error for WettingError {}

pub fn wall_normal_gradient<T: Real>(grad_parallel_abs: T, params: ContactAngleParams<T>) -> T {
    -grad_parallel_abs * params.cot_theta()
}

pub fn spherical_cap_phi<T: Real>(
    x: f64,
    y: f64,
    z: f64,
    center: [f64; 3],
    radius: f64,
    interface_width: f64,
) -> T {
    let dx = x - center[0];
    let dy = y - center[1];
    let dz = z - center[2];
    let r = (dx * dx + dy * dy + dz * dz).sqrt();
    T::r(0.5 * (1.0 - (2.0 * (r - radius) / interface_width).tanh()))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn neutral_90_degrees_case_stays_symmetric() {
        let params = ContactAngleParams { theta_deg: 90.0f64 }
            .validate()
            .unwrap();
        assert!(wall_normal_gradient(0.25, params).abs() < 1.0e-12);
    }

    #[test]
    fn hydrophobic_and_hydrophilic_boundary_signs_differ() {
        let hydrophilic = ContactAngleParams { theta_deg: 60.0f64 }
            .validate()
            .unwrap();
        let hydrophobic = ContactAngleParams {
            theta_deg: 120.0f64,
        }
        .validate()
        .unwrap();
        assert!(wall_normal_gradient(0.25, hydrophilic) < 0.0);
        assert!(wall_normal_gradient(0.25, hydrophobic) > 0.0);
    }

    #[test]
    fn rejects_invalid_static_contact_angle() {
        assert!(ContactAngleParams { theta_deg: 0.0f64 }.validate().is_err());
        assert!(ContactAngleParams {
            theta_deg: 180.0f64
        }
        .validate()
        .is_err());
    }
}
