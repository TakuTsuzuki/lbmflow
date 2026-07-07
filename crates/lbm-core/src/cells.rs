//! Cell tracer population and checkpoint state for bioprocess exposure QOIs.

use crate::damage::{DamageModelError, ShearDamageModel};
use crate::particles::{Particle, ParticleError, ParticleSet, Sample};
use serde::{Deserialize, Serialize};

/// Checkpoint payload producer for cell-tracer-like particle state.
pub trait CellCheckpointSection {
    fn checkpoint_bytes(&self) -> Option<Vec<u8>>;
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct NoCellTracers;

impl CellCheckpointSection for NoCellTracers {
    fn checkpoint_bytes(&self) -> Option<Vec<u8>> {
        None
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct CellFieldSample {
    pub velocity: [f64; 3],
    pub solid: bool,
    pub gamma_dot_1_s: f64,
    pub viscous_stress_pa: f64,
    pub oxygen_cl: Option<f64>,
}

impl From<CellFieldSample> for Sample {
    fn from(value: CellFieldSample) -> Self {
        Self {
            u: value.velocity,
            solid: value.solid,
        }
    }
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct CellTracer {
    pub position: [f64; 3],
    pub velocity: [f64; 3],
    pub gamma_dot_1_s: f64,
    pub viscous_stress_pa: f64,
    pub oxygen_cl: Option<f64>,
    pub shear_exposure: f64,
    pub residence_time_above_threshold_s: f64,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct CellTracerPopulation {
    pub tracers: Vec<CellTracer>,
    pub seed: u64,
}

#[derive(Clone, Debug, PartialEq)]
pub enum CellTracerError {
    NoFluidCells,
    Particle(ParticleError),
    Damage(DamageModelError),
}

impl std::fmt::Display for CellTracerError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::NoFluidCells => write!(f, "cell tracer seeding requires at least one fluid cell"),
            Self::Particle(e) => write!(f, "{e}"),
            Self::Damage(e) => write!(f, "{e}"),
        }
    }
}

impl std::error::Error for CellTracerError {}

impl From<ParticleError> for CellTracerError {
    fn from(value: ParticleError) -> Self {
        Self::Particle(value)
    }
}

impl From<DamageModelError> for CellTracerError {
    fn from(value: DamageModelError) -> Self {
        Self::Damage(value)
    }
}

impl CellTracerPopulation {
    pub fn seed_deterministic(
        count: usize,
        seed: u64,
        dims: [usize; 3],
        solid: impl Fn(usize, usize, usize) -> bool,
    ) -> Result<Self, CellTracerError> {
        assert!(dims[0] > 0 && dims[1] > 0 && dims[2] > 0);
        let mut fluid = Vec::new();
        for z in 0..dims[2] {
            for y in 0..dims[1] {
                for x in 0..dims[0] {
                    if !solid(x, y, z) {
                        fluid.push([x as f64 + 0.5, y as f64 + 0.5, z as f64 + 0.5]);
                    }
                }
            }
        }
        if count > 0 && fluid.is_empty() {
            return Err(CellTracerError::NoFluidCells);
        }

        let mut rng = SplitMix64::new(seed);
        let mut tracers = Vec::with_capacity(count);
        for _ in 0..count {
            let pos = fluid[(rng.next_u64() as usize) % fluid.len()];
            tracers.push(CellTracer {
                position: pos,
                velocity: [0.0; 3],
                gamma_dot_1_s: 0.0,
                viscous_stress_pa: 0.0,
                oxygen_cl: None,
                shear_exposure: 0.0,
                residence_time_above_threshold_s: 0.0,
            });
        }
        Ok(Self { tracers, seed })
    }

    pub fn len(&self) -> usize {
        self.tracers.len()
    }

    pub fn is_empty(&self) -> bool {
        self.tracers.is_empty()
    }

    pub fn step<F>(
        &mut self,
        sample: F,
        dt_s: f64,
        damage: Option<&ShearDamageModel>,
    ) -> Result<(), CellTracerError>
    where
        F: Fn([f64; 3]) -> CellFieldSample,
    {
        assert!(
            dt_s.is_finite() && dt_s > 0.0,
            "dt_s must be finite and positive"
        );
        let particles: Vec<Particle> = self
            .tracers
            .iter()
            .map(|tracer| Particle {
                pos: tracer.position,
                vel: tracer.velocity,
                d: 1.0,
                rho_p: 1.0,
                exposure: tracer.shear_exposure,
            })
            .collect();
        let mut set = ParticleSet::new(particles, 1.0, 1.0, [0.0; 3]);
        set.step_massless(|pos| sample(pos).into(), None::<fn([f64; 3]) -> f64>, dt_s)?;

        for (tracer, particle) in self.tracers.iter_mut().zip(set.particles) {
            tracer.position = particle.pos;
            tracer.velocity = particle.vel;
            let fields = sample(tracer.position);
            tracer.gamma_dot_1_s = fields.gamma_dot_1_s;
            tracer.viscous_stress_pa = fields.viscous_stress_pa;
            tracer.oxygen_cl = fields.oxygen_cl;
            if let Some(model) = damage {
                let inc =
                    model.increment(fields.viscous_stress_pa, fields.gamma_dot_1_s, None, dt_s)?;
                tracer.shear_exposure += inc.exposure_increment;
                if inc.above_threshold {
                    tracer.residence_time_above_threshold_s += dt_s;
                }
            }
        }
        Ok(())
    }
}

impl CellCheckpointSection for CellTracerPopulation {
    fn checkpoint_bytes(&self) -> Option<Vec<u8>> {
        Some(serde_json::to_vec(self).expect("cell tracer checkpoint serialization is infallible"))
    }
}

#[derive(Clone, Copy, Debug)]
struct SplitMix64 {
    state: u64,
}

impl SplitMix64 {
    fn new(seed: u64) -> Self {
        Self { state: seed }
    }

    fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn open(_: usize, _: usize, _: usize) -> bool {
        false
    }

    #[test]
    fn uniform_flow_advects_tracers_analytically() {
        let mut pop = CellTracerPopulation::seed_deterministic(4, 7, [8, 8, 2], open).unwrap();
        let initial: Vec<[f64; 3]> = pop.tracers.iter().map(|t| t.position).collect();
        let u = [0.2, -0.1, 0.05];
        pop.step(
            |_| CellFieldSample {
                velocity: u,
                solid: false,
                gamma_dot_1_s: 3.0,
                viscous_stress_pa: 0.006,
                oxygen_cl: Some(0.8),
            },
            0.5,
            None,
        )
        .unwrap();
        for (tracer, p0) in pop.tracers.iter().zip(initial) {
            for a in 0..3 {
                let want = p0[a] + 0.5 * u[a];
                assert!((tracer.position[a] - want).abs() < 1e-12);
            }
            assert_eq!(tracer.velocity, u);
        }
    }

    #[test]
    fn deterministic_seeding_matches_seed_and_changes_with_seed() {
        let a = CellTracerPopulation::seed_deterministic(10, 123, [5, 5, 2], open).unwrap();
        let b = CellTracerPopulation::seed_deterministic(10, 123, [5, 5, 2], open).unwrap();
        let c = CellTracerPopulation::seed_deterministic(10, 124, [5, 5, 2], open).unwrap();
        assert_eq!(a.tracers, b.tracers);
        assert_ne!(a.tracers, c.tracers);
    }

    #[test]
    fn tracer_count_preserved_and_sampled_fields_recorded() {
        let mut pop = CellTracerPopulation::seed_deterministic(6, 9, [6, 6, 1], open).unwrap();
        pop.step(
            |p| CellFieldSample {
                velocity: [0.0; 3],
                solid: false,
                gamma_dot_1_s: p[0] + p[1],
                viscous_stress_pa: 2.0 * (p[0] + p[1]),
                oxygen_cl: Some(p[2]),
            },
            1.0,
            None,
        )
        .unwrap();
        assert_eq!(pop.len(), 6);
        for tracer in &pop.tracers {
            let g = tracer.position[0] + tracer.position[1];
            assert_eq!(tracer.gamma_dot_1_s, g);
            assert_eq!(tracer.viscous_stress_pa, 2.0 * g);
            assert_eq!(tracer.oxygen_cl, Some(tracer.position[2]));
        }
    }

    #[test]
    fn closed_vessel_solid_boundary_reflects_tracer() {
        let mut pop = CellTracerPopulation {
            tracers: vec![CellTracer {
                position: [0.75, 0.5, 0.5],
                velocity: [0.0; 3],
                gamma_dot_1_s: 0.0,
                viscous_stress_pa: 0.0,
                oxygen_cl: None,
                shear_exposure: 0.0,
                residence_time_above_threshold_s: 0.0,
            }],
            seed: 0,
        };
        pop.step(
            |p| CellFieldSample {
                velocity: [-1.0, 0.0, 0.0],
                solid: p[0] < 0.0,
                gamma_dot_1_s: 0.0,
                viscous_stress_pa: 0.0,
                oxygen_cl: None,
            },
            1.0,
            None,
        )
        .unwrap();
        assert!(pop.tracers[0].position[0] >= 0.0);
    }
}
