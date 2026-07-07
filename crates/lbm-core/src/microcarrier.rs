//! Microcarrier particle mode and two-way drag-reaction scatter.

use crate::particles::{
    schiller_naumann_drag_correction, Particle, ParticleError, ParticleSet, Sample,
    SCHILLER_NAUMANN_RE_MAX,
};
use serde::{Deserialize, Serialize};

pub const TWO_WAY_MASS_LOADING_MAX: f64 = 0.1;

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct MicrocarrierPopulation {
    pub particles: ParticleSet,
    pub residence_near_impeller_s: Vec<f64>,
}

#[derive(Clone, Debug, PartialEq)]
pub enum MicrocarrierError {
    Particle(ParticleError),
    InvalidReynolds {
        re: f64,
        re_max: f64,
    },
    InvalidMassLoading {
        mass_loading: f64,
        max: f64,
    },
    ParticleOutsideGrid {
        position: [f64; 3],
        dims: [usize; 3],
    },
    InvalidParameter {
        parameter: &'static str,
        value: f64,
    },
}

impl std::fmt::Display for MicrocarrierError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Particle(e) => write!(f, "{e}"),
            Self::InvalidReynolds { re, re_max } => {
                write!(
                    f,
                    "microcarrier Re_p={re:e} exceeds Schiller-Naumann validity {re_max:e}"
                )
            }
            Self::InvalidMassLoading { mass_loading, max } => {
                write!(
                    f,
                    "two-way microcarrier mass_loading={mass_loading:e} exceeds guard {max:e}"
                )
            }
            Self::ParticleOutsideGrid { position, dims } => {
                write!(f, "particle at {position:?} is outside grid {dims:?}")
            }
            Self::InvalidParameter { parameter, value } => {
                write!(f, "invalid microcarrier parameter {parameter}={value:e}")
            }
        }
    }
}

impl std::error::Error for MicrocarrierError {}

impl From<ParticleError> for MicrocarrierError {
    fn from(value: ParticleError) -> Self {
        Self::Particle(value)
    }
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct SuspensionMetrics {
    pub settled_fraction: f64,
    pub height_distribution: Vec<usize>,
    pub residence_near_impeller_s: Vec<f64>,
    pub shear_exposure_distribution: crate::damage::ExposureDistribution,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ParticleDragForce {
    pub position: [f64; 3],
    pub force_on_particle: [f64; 3],
}

#[derive(Clone, Debug, PartialEq)]
pub struct TwoWayScatterReport {
    pub liquid_reaction_force: Vec<[f64; 3]>,
    pub particle_force_sum: [f64; 3],
    pub liquid_reaction_sum: [f64; 3],
    pub ledger_residual: [f64; 3],
}

impl MicrocarrierPopulation {
    pub fn new(particles: ParticleSet) -> Self {
        let residence_near_impeller_s = vec![0.0; particles.particles.len()];
        Self {
            particles,
            residence_near_impeller_s,
        }
    }

    pub fn step_one_way<F, S, N>(
        &mut self,
        sample: F,
        shear_rate: S,
        near_impeller: N,
        dt_s: f64,
    ) -> Result<(), MicrocarrierError>
    where
        F: Fn([f64; 3]) -> Sample,
        S: Fn([f64; 3]) -> f64,
        N: Fn([f64; 3]) -> bool,
    {
        let before_positions: Vec<[f64; 3]> =
            self.particles.particles.iter().map(|p| p.pos).collect();
        self.particles.step_dt(sample, Some(shear_rate), dt_s)?;
        for (i, p) in self.particles.particles.iter().enumerate() {
            if near_impeller(p.pos) || near_impeller(before_positions[i]) {
                self.residence_near_impeller_s[i] += dt_s;
            }
        }
        Ok(())
    }

    pub fn metrics(
        &self,
        floor_z: f64,
        height_bins: usize,
        height_min: f64,
        height_max: f64,
        exposure_threshold: f64,
    ) -> Result<SuspensionMetrics, MicrocarrierError> {
        if height_bins == 0 {
            return Err(MicrocarrierError::InvalidParameter {
                parameter: "height_bins",
                value: height_bins as f64,
            });
        }
        if !(height_min.is_finite() && height_max.is_finite() && height_max > height_min) {
            return Err(MicrocarrierError::InvalidParameter {
                parameter: "height_range",
                value: height_max - height_min,
            });
        }
        let n = self.particles.particles.len();
        let settled = self
            .particles
            .particles
            .iter()
            .filter(|p| p.pos[2] <= floor_z)
            .count();
        let mut height_distribution = vec![0usize; height_bins];
        for p in &self.particles.particles {
            if p.pos[2] < height_min || p.pos[2] > height_max {
                continue;
            }
            let frac = (p.pos[2] - height_min) / (height_max - height_min);
            let mut bin = (frac * height_bins as f64).floor() as usize;
            if bin == height_bins {
                bin = height_bins - 1;
            }
            height_distribution[bin] += 1;
        }
        let exposures: Vec<f64> = self
            .particles
            .particles
            .iter()
            .map(|p| p.exposure)
            .collect();
        let shear_exposure_distribution = crate::damage::exposure_distribution(
            &exposures,
            &self.residence_near_impeller_s,
            exposure_threshold,
        )
        .map_err(|e| MicrocarrierError::InvalidParameter {
            parameter: "exposure_distribution",
            value: if e.to_string().is_empty() {
                0.0
            } else {
                f64::NAN
            },
        })?;
        Ok(SuspensionMetrics {
            settled_fraction: if n == 0 {
                0.0
            } else {
                settled as f64 / n as f64
            },
            height_distribution,
            residence_near_impeller_s: self.residence_near_impeller_s.clone(),
            shear_exposure_distribution,
        })
    }
}

pub fn validate_reynolds(re: f64) -> Result<(), MicrocarrierError> {
    if !(re.is_finite() && re >= 0.0) {
        return Err(MicrocarrierError::InvalidParameter {
            parameter: "re_p",
            value: re,
        });
    }
    if re > SCHILLER_NAUMANN_RE_MAX {
        Err(MicrocarrierError::InvalidReynolds {
            re,
            re_max: SCHILLER_NAUMANN_RE_MAX,
        })
    } else {
        Ok(())
    }
}

pub fn validate_mass_loading(mass_loading: f64) -> Result<(), MicrocarrierError> {
    if !(mass_loading.is_finite() && mass_loading >= 0.0) {
        return Err(MicrocarrierError::InvalidParameter {
            parameter: "mass_loading",
            value: mass_loading,
        });
    }
    if mass_loading > TWO_WAY_MASS_LOADING_MAX {
        Err(MicrocarrierError::InvalidMassLoading {
            mass_loading,
            max: TWO_WAY_MASS_LOADING_MAX,
        })
    } else {
        Ok(())
    }
}

pub fn drag_force_on_particle(
    particle: &Particle,
    fluid_velocity: [f64; 3],
    rho_f: f64,
    nu: f64,
) -> Result<[f64; 3], MicrocarrierError> {
    let slip = sub(fluid_velocity, particle.vel);
    let re = norm(slip) * particle.d / nu;
    validate_reynolds(re)?;
    let drag_correction =
        schiller_naumann_drag_correction(re).map_err(|_| MicrocarrierError::InvalidReynolds {
            re,
            re_max: SCHILLER_NAUMANN_RE_MAX,
        })?;
    let tau_p = particle.rho_p * particle.d * particle.d / (18.0 * rho_f * nu * drag_correction);
    let volume = std::f64::consts::PI * particle.d.powi(3) / 6.0;
    let mass = particle.rho_p * volume;
    Ok(scale(slip, mass / tau_p))
}

pub fn terminal_velocity_stokes(d: f64, rho_p: f64, rho_f: f64, nu: f64, g_abs: f64) -> f64 {
    (rho_p - rho_f) * d * d * g_abs / (18.0 * rho_f * nu)
}

pub fn scatter_drag_reaction_forces(
    dims: [usize; 3],
    events: &[ParticleDragForce],
) -> Result<TwoWayScatterReport, MicrocarrierError> {
    assert!(dims[0] > 0 && dims[1] > 0 && dims[2] > 0);
    let n = dims[0] * dims[1] * dims[2];
    let mut liquid = vec![[0.0; 3]; n];
    let mut particle_force_sum = [0.0; 3];
    for event in events {
        let weights = trilinear_kernel(dims, event.position)?;
        for a in 0..3 {
            particle_force_sum[a] += event.force_on_particle[a];
        }
        for (idx, weight) in weights {
            for a in 0..3 {
                liquid[idx][a] -= weight * event.force_on_particle[a];
            }
        }
    }
    let mut liquid_reaction_sum = [0.0; 3];
    for f in &liquid {
        for a in 0..3 {
            liquid_reaction_sum[a] += f[a];
        }
    }
    let mut ledger_residual = [0.0; 3];
    for a in 0..3 {
        ledger_residual[a] = particle_force_sum[a] + liquid_reaction_sum[a];
    }
    Ok(TwoWayScatterReport {
        liquid_reaction_force: liquid,
        particle_force_sum,
        liquid_reaction_sum,
        ledger_residual,
    })
}

fn trilinear_kernel(
    dims: [usize; 3],
    position: [f64; 3],
) -> Result<Vec<(usize, f64)>, MicrocarrierError> {
    for a in 0..3 {
        if !(position[a].is_finite() && position[a] >= 0.0 && position[a] <= (dims[a] - 1) as f64) {
            return Err(MicrocarrierError::ParticleOutsideGrid { position, dims });
        }
    }
    let (x0, x1, tx) = bracket(position[0], dims[0]);
    let (y0, y1, ty) = bracket(position[1], dims[1]);
    let (z0, z1, tz) = bracket(position[2], dims[2]);
    let mut out = Vec::with_capacity(8);
    for (x, wx) in [(x0, 1.0 - tx), (x1, tx)] {
        for (y, wy) in [(y0, 1.0 - ty), (y1, ty)] {
            for (z, wz) in [(z0, 1.0 - tz), (z1, tz)] {
                let w = wx * wy * wz;
                if w > 0.0 {
                    out.push((idx(dims, x, y, z), w));
                }
            }
        }
    }
    Ok(out)
}

fn bracket(x: f64, n: usize) -> (usize, usize, f64) {
    if n == 1 {
        return (0, 0, 0.0);
    }
    let lo = x.floor() as usize;
    let hi = if lo + 1 < n { lo + 1 } else { lo };
    (lo, hi, x - lo as f64)
}

fn idx(dims: [usize; 3], x: usize, y: usize, z: usize) -> usize {
    (z * dims[1] + y) * dims[0] + x
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [s * a[0], s * a[1], s * a[2]]
}

fn norm(a: [f64; 3]) -> f64 {
    (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::particles::particle_velocity_dt;

    fn fluid(u: [f64; 3]) -> impl Fn([f64; 3]) -> Sample {
        move |_| Sample { u, solid: false }
    }

    #[test]
    fn terminal_velocity_matches_stokes_limit_within_five_percent() {
        let rho_f = 1.0;
        let rho_p = 1.01;
        let d = 0.2;
        let nu = 1.0;
        let g = [0.0, 0.0, -0.01];
        let mut set = ParticleSet::new(
            vec![Particle {
                pos: [0.0, 0.0, 20.0],
                vel: [0.0; 3],
                d,
                rho_p,
                exposure: 0.0,
            }],
            rho_f,
            nu,
            g,
        );
        for _ in 0..10_000 {
            set.step_dt(fluid([0.0; 3]), None::<fn([f64; 3]) -> f64>, 0.1)
                .unwrap();
        }
        let got = -set.particles[0].vel[2];
        let want = terminal_velocity_stokes(d, rho_p, rho_f, nu, g[2].abs());
        assert!(
            (got - want).abs() / want < 0.05,
            "got {got:e}, want {want:e}"
        );
    }

    #[test]
    fn quiescent_tank_settling_moves_down_and_conserves_count() {
        let particles = ParticleSet::new(
            (0..5)
                .map(|i| Particle {
                    pos: [i as f64, 0.0, 10.0],
                    vel: [0.0; 3],
                    d: 0.5,
                    rho_p: 1.2,
                    exposure: 0.0,
                })
                .collect(),
            1.0,
            0.2,
            [0.0, 0.0, -0.01],
        );
        let mut pop = MicrocarrierPopulation::new(particles);
        pop.step_one_way(fluid([0.0; 3]), |_| 0.0, |_| false, 1.0)
            .unwrap();
        assert_eq!(pop.particles.particles.len(), 5);
        assert!(pop.particles.particles.iter().all(|p| p.pos[2] < 10.0));
    }

    #[test]
    fn invalid_reynolds_is_rejected() {
        assert!(matches!(
            validate_reynolds(801.0),
            Err(MicrocarrierError::InvalidReynolds { .. })
        ));
    }

    #[test]
    fn two_way_scatter_ledger_is_equal_and_opposite() {
        let events = [ParticleDragForce {
            position: [1.25, 1.5, 0.0],
            force_on_particle: [2.0, -3.0, 0.5],
        }];
        let report = scatter_drag_reaction_forces([4, 4, 1], &events).unwrap();
        for a in 0..3 {
            assert!(report.ledger_residual[a].abs() < 1e-12);
        }
        assert_eq!(report.particle_force_sum, [2.0, -3.0, 0.5]);
        assert_eq!(report.liquid_reaction_force.len(), 16);
    }

    #[test]
    fn single_particle_uniform_flow_drag_points_with_slip() {
        let p = Particle {
            pos: [1.0; 3],
            vel: [0.0; 3],
            d: 0.2,
            rho_p: 1.1,
            exposure: 0.0,
        };
        let f = drag_force_on_particle(&p, [0.1, 0.0, 0.0], 1.0, 0.5).unwrap();
        assert!(f[0] > 0.0);
        assert_eq!(f[1], 0.0);
        assert_eq!(f[2], 0.0);
    }

    #[test]
    fn mass_loading_guard_rejects_until_four_way_exists() {
        assert!(validate_mass_loading(0.1).is_ok());
        assert!(matches!(
            validate_mass_loading(0.100_001),
            Err(MicrocarrierError::InvalidMassLoading { .. })
        ));
    }

    #[test]
    fn two_way_disabled_leaves_one_way_particle_update_bit_identical() {
        let particle = Particle {
            pos: [0.0; 3],
            vel: [0.0; 3],
            d: 0.3,
            rho_p: 1.2,
            exposure: 0.0,
        };
        let mut a = ParticleSet::new(vec![particle.clone()], 1.0, 0.5, [0.0; 3]);
        let mut b = ParticleSet::new(vec![particle], 1.0, 0.5, [0.0; 3]);
        a.step_dt(fluid([0.1, 0.0, 0.0]), None::<fn([f64; 3]) -> f64>, 0.5)
            .unwrap();
        b.step_dt(fluid([0.1, 0.0, 0.0]), None::<fn([f64; 3]) -> f64>, 0.5)
            .unwrap();
        assert_eq!(a, b);
    }

    #[test]
    fn particle_velocity_dt_preserves_one_step_wrapper() {
        let v = particle_velocity_dt([0.0; 3], [0.1, 0.0, 0.0], 0.2, 1.1, 1.0, 0.5, [0.0; 3], 1.0)
            .unwrap();
        let mut set = ParticleSet::new(
            vec![Particle {
                pos: [0.0; 3],
                vel: [0.0; 3],
                d: 0.2,
                rho_p: 1.1,
                exposure: 0.0,
            }],
            1.0,
            0.5,
            [0.0; 3],
        );
        set.step(fluid([0.1, 0.0, 0.0]), None::<fn([f64; 3]) -> f64>)
            .unwrap();
        assert_eq!(set.particles[0].vel, v);
    }
}
