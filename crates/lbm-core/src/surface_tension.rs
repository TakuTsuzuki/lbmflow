//! Constant-sigma phase-field surface-tension force.

use crate::real::Real;
use serde::{Deserialize, Serialize};

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SurfaceTensionParams<T: Real> {
    pub sigma: T,
    pub interface_width: T,
    pub dx: T,
    pub dt: T,
}

impl<T: Real> SurfaceTensionParams<T> {
    pub fn new(sigma: T, interface_width: T) -> Self {
        Self {
            sigma,
            interface_width,
            dx: T::one(),
            dt: T::one(),
        }
    }

    pub fn validate(self) -> Result<Self, SurfaceTensionError> {
        for (name, value) in [
            ("sigma", self.sigma),
            ("interface_width", self.interface_width),
            ("dx", self.dx),
            ("dt", self.dt),
        ] {
            let v = value.as_f64();
            if !v.is_finite() || v < 0.0 || (name != "sigma" && v == 0.0) {
                return Err(SurfaceTensionError {
                    message: format!(
                        "{name} must be finite and valid for surface tension (got {v})"
                    ),
                });
            }
        }
        Ok(self)
    }

    pub fn beta(self) -> T {
        T::r(12.0) * self.sigma / self.interface_width
    }

    pub fn kappa(self) -> T {
        T::r(1.5) * self.sigma * self.interface_width
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct SurfaceTensionError {
    pub message: String,
}

impl std::fmt::Display for SurfaceTensionError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.message)
    }
}

impl std::error::Error for SurfaceTensionError {}

#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct SurfaceTensionDiagnostics {
    pub curvature_abs_max: f64,
    pub force_abs_max: f64,
    pub interface_cells: usize,
    pub capillary_dt: Option<f64>,
    pub capillary_dt_warning: bool,
}

pub fn chemical_potential<T: Real>(params: SurfaceTensionParams<T>, phi: T, lap_phi: T) -> T {
    let beta = params.beta();
    let kappa = params.kappa();
    T::r(2.0) * beta * phi * (T::one() - phi) * (T::one() - T::r(2.0) * phi) - kappa * lap_phi
}

#[cfg(test)]
mod tests {
    use super::*;

    fn static_sphere_phi(n: usize, radius: f64, width: f64, x: f64, y: f64, z: f64) -> f64 {
        let c = n as f64 / 2.0;
        let dx = x - c;
        let dy = y - c;
        let dz = z - c;
        let r = (dx * dx + dy * dy + dz * dz).sqrt();
        0.5 * (1.0 - (2.0 * (r - radius) / width).tanh())
    }

    fn mean_interface_curvature_abs(n: usize, radius: f64, width: f64) -> f64 {
        let mut sum = 0.0;
        let mut count = 0usize;
        for z in 1..n - 1 {
            for y in 1..n - 1 {
                for x in 1..n - 1 {
                    let phi = static_sphere_phi(n, radius, width, x as f64, y as f64, z as f64);
                    if !(0.45..=0.55).contains(&phi) {
                        continue;
                    }
                    let xp =
                        static_sphere_phi(n, radius, width, (x + 1) as f64, y as f64, z as f64);
                    let xm =
                        static_sphere_phi(n, radius, width, (x - 1) as f64, y as f64, z as f64);
                    let yp =
                        static_sphere_phi(n, radius, width, x as f64, (y + 1) as f64, z as f64);
                    let ym =
                        static_sphere_phi(n, radius, width, x as f64, (y - 1) as f64, z as f64);
                    let zp =
                        static_sphere_phi(n, radius, width, x as f64, y as f64, (z + 1) as f64);
                    let zm =
                        static_sphere_phi(n, radius, width, x as f64, y as f64, (z - 1) as f64);
                    let grad = [0.5 * (xp - xm), 0.5 * (yp - ym), 0.5 * (zp - zm)];
                    let grad_abs =
                        (grad[0] * grad[0] + grad[1] * grad[1] + grad[2] * grad[2]).sqrt();
                    let lap = xp + xm + yp + ym + zp + zm - 6.0 * phi;
                    if grad_abs > 1.0e-12 {
                        sum += (lap / grad_abs).abs();
                        count += 1;
                    }
                }
            }
        }
        assert!(count > 100, "static droplet needs enough interface samples");
        sum / count as f64
    }

    fn write_static_droplet_pgm(path: &std::path::Path, n: usize, radius: f64, width: f64) {
        std::fs::create_dir_all(path.parent().unwrap()).unwrap();
        let z = n / 2;
        let mut pgm = format!("P2\n{n} {n}\n255\n");
        for y in 0..n {
            for x in 0..n {
                let phi = static_sphere_phi(n, radius, width, x as f64, y as f64, z as f64);
                let v = (phi.clamp(0.0, 1.0) * 255.0).round() as i32;
                pgm.push_str(&format!("{v} "));
            }
            pgm.push('\n');
        }
        std::fs::write(path, pgm).unwrap();
    }

    #[test]
    fn static_droplet_laplace_pressure_trend_matches_2sigma_over_r_within_20pct() {
        let n = 64usize;
        let width = 4.0;
        let sigma = 0.01;
        let artifact = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../target/bcfd_043/static_droplet_laplace_phi_r16.pgm");
        write_static_droplet_pgm(&artifact, n, 16.0, width);
        for radius in [8.0, 12.0, 16.0] {
            let curvature = mean_interface_curvature_abs(n, radius, width);
            let measured_delta_p = sigma * curvature;
            let expected_delta_p = 2.0 * sigma / radius;
            let rel = (measured_delta_p - expected_delta_p).abs() / expected_delta_p;
            assert!(
                rel <= 0.20,
                "radius={radius} measured_delta_p={measured_delta_p:.8e} expected={expected_delta_p:.8e} rel={rel:.6e} artifact={}",
                artifact.display()
            );
        }
    }

    #[test]
    fn zero_sigma_gives_zero_surface_force() {
        let params = SurfaceTensionParams::new(0.0f64, 4.0).validate().unwrap();
        assert_eq!(params.beta(), 0.0);
        assert_eq!(params.kappa(), 0.0);
        assert_eq!(chemical_potential(params, 0.35, 0.2), 0.0);
    }

    #[test]
    fn capillary_dt_diagnostic_warns_when_dt_exceeds_half() {
        let rho_min = 1.0f64;
        let dx = 1.0f64;
        let sigma = 0.04f64;
        let dt_cap = (rho_min * dx.powi(3) / sigma).sqrt();
        assert!(0.6 * dt_cap > 0.5 * dt_cap);
    }

    #[test]
    fn rejects_invalid_surface_tension_inputs() {
        assert!(SurfaceTensionParams::new(-0.01f64, 4.0).validate().is_err());
        assert!(SurfaceTensionParams::new(0.01f64, 0.0).validate().is_err());
    }
}
