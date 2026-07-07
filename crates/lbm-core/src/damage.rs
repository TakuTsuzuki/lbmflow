//! Shear exposure and damage-risk reductions for cell tracers.

use crate::stress::{percentile_summary, PercentileSummary};
use serde::{Deserialize, Serialize};

#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum DamageThreshold {
    ViscousStressPa { threshold_pa: f64 },
    ShearRate { threshold_1_s: f64 },
    EnergyDissipationPlaceholder { threshold_w_kg: f64 },
}

#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct ShearDamageModel {
    pub threshold: DamageThreshold,
    pub exponent: f64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DamageIncrement {
    pub exposure_increment: f64,
    pub above_threshold: bool,
}

#[derive(Clone, Debug, PartialEq)]
pub enum DamageModelError {
    InvalidParameter { parameter: &'static str, value: f64 },
    EnergyDissipationUnavailable,
    EmptyDistribution,
}

impl std::fmt::Display for DamageModelError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InvalidParameter { parameter, value } => {
                write!(f, "invalid damage-model parameter {parameter}={value:e}")
            }
            Self::EnergyDissipationUnavailable => write!(
                f,
                "energy-dissipation damage threshold requires an epsilon field"
            ),
            Self::EmptyDistribution => write!(f, "exposure distribution has no finite values"),
        }
    }
}

impl std::error::Error for DamageModelError {}

#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct ExposureDistribution {
    pub p50: f64,
    pub p90: f64,
    pub p95: f64,
    pub p99: f64,
    pub max: f64,
    pub fraction_above_threshold: f64,
    pub residence_time_above_threshold_s: f64,
}

impl ShearDamageModel {
    pub fn stress_threshold(threshold_pa: f64, exponent: f64) -> Result<Self, DamageModelError> {
        let model = Self {
            threshold: DamageThreshold::ViscousStressPa { threshold_pa },
            exponent,
        };
        model.validate()?;
        Ok(model)
    }

    pub fn validate(&self) -> Result<(), DamageModelError> {
        if !(self.exponent.is_finite() && self.exponent > 0.0) {
            return Err(DamageModelError::InvalidParameter {
                parameter: "exponent",
                value: self.exponent,
            });
        }
        match self.threshold {
            DamageThreshold::ViscousStressPa { threshold_pa } => {
                validate_non_negative("threshold_pa", threshold_pa)
            }
            DamageThreshold::ShearRate { threshold_1_s } => {
                validate_non_negative("threshold_1_s", threshold_1_s)
            }
            DamageThreshold::EnergyDissipationPlaceholder { threshold_w_kg } => {
                validate_non_negative("threshold_w_kg", threshold_w_kg)
            }
        }
    }

    pub fn increment(
        &self,
        viscous_stress_pa: f64,
        gamma_dot_1_s: f64,
        epsilon_w_kg: Option<f64>,
        dt_s: f64,
    ) -> Result<DamageIncrement, DamageModelError> {
        self.validate()?;
        if !(dt_s.is_finite() && dt_s > 0.0) {
            return Err(DamageModelError::InvalidParameter {
                parameter: "dt_s",
                value: dt_s,
            });
        }
        let (value, threshold) = match self.threshold {
            DamageThreshold::ViscousStressPa { threshold_pa } => (viscous_stress_pa, threshold_pa),
            DamageThreshold::ShearRate { threshold_1_s } => (gamma_dot_1_s, threshold_1_s),
            DamageThreshold::EnergyDissipationPlaceholder { threshold_w_kg } => (
                epsilon_w_kg.ok_or(DamageModelError::EnergyDissipationUnavailable)?,
                threshold_w_kg,
            ),
        };
        if !value.is_finite() {
            return Err(DamageModelError::InvalidParameter {
                parameter: "damage_threshold_value",
                value,
            });
        }
        let excess = (value - threshold).max(0.0);
        Ok(DamageIncrement {
            exposure_increment: excess.powf(self.exponent) * dt_s,
            above_threshold: excess > 0.0,
        })
    }
}

pub fn exposure_distribution(
    exposures: &[f64],
    residence_times_s: &[f64],
    exposure_threshold: f64,
) -> Result<ExposureDistribution, DamageModelError> {
    if !(exposure_threshold.is_finite() && exposure_threshold >= 0.0) {
        return Err(DamageModelError::InvalidParameter {
            parameter: "exposure_threshold",
            value: exposure_threshold,
        });
    }
    if exposures.len() != residence_times_s.len() {
        return Err(DamageModelError::InvalidParameter {
            parameter: "residence_times_s.len",
            value: residence_times_s.len() as f64,
        });
    }
    let summary = percentile_summary(exposures, Some(exposure_threshold))
        .ok_or(DamageModelError::EmptyDistribution)?;
    Ok(distribution_from_summary(
        summary,
        residence_times_s.iter().copied().sum(),
    ))
}

fn distribution_from_summary(
    summary: PercentileSummary,
    residence_time_above_threshold_s: f64,
) -> ExposureDistribution {
    ExposureDistribution {
        p50: summary.p50,
        p90: summary.p90,
        p95: summary.p95,
        p99: summary.p99,
        max: summary.max,
        fraction_above_threshold: summary.fraction_above_threshold.unwrap_or(0.0),
        residence_time_above_threshold_s,
    }
}

fn validate_non_negative(parameter: &'static str, value: f64) -> Result<(), DamageModelError> {
    if value.is_finite() && value >= 0.0 {
        Ok(())
    } else {
        Err(DamageModelError::InvalidParameter { parameter, value })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn constant_shear_exposure_matches_analytic() {
        let model = ShearDamageModel::stress_threshold(0.2, 2.0).unwrap();
        let mut exposure = 0.0;
        for _ in 0..10 {
            exposure += model
                .increment(0.7, 0.0, None, 0.1)
                .unwrap()
                .exposure_increment;
        }
        let want = (0.7_f64 - 0.2).powi(2);
        assert!((exposure - want).abs() < 1e-12);
    }

    #[test]
    fn below_threshold_gives_zero_exactly() {
        let model = ShearDamageModel::stress_threshold(1.0, 1.7).unwrap();
        let inc = model.increment(0.25, 0.0, None, 3.0).unwrap();
        assert_eq!(inc.exposure_increment, 0.0);
        assert!(!inc.above_threshold);
    }

    #[test]
    fn percentile_reducer_reports_distribution_not_max_alone() {
        let d = exposure_distribution(&[0.0, 1.0, 2.0, 3.0, 4.0], &[0.0, 0.5, 1.0, 1.5, 2.0], 2.0)
            .unwrap();
        assert_eq!(d.p50, 2.0);
        assert_eq!(d.max, 4.0);
        assert_eq!(d.fraction_above_threshold, 2.0 / 5.0);
        assert_eq!(d.residence_time_above_threshold_s, 5.0);
        assert!(d.p90 > d.p50);
    }

    #[test]
    fn halved_dt_is_time_step_invariant_for_constant_shear() {
        let model = ShearDamageModel::stress_threshold(0.1, 1.3).unwrap();
        let integrate = |dt: f64| {
            let mut e = 0.0;
            let n = (2.0 / dt) as usize;
            for _ in 0..n {
                e += model
                    .increment(0.4, 0.0, None, dt)
                    .unwrap()
                    .exposure_increment;
            }
            e
        };
        let coarse = integrate(0.2);
        let fine = integrate(0.1);
        assert!((coarse - fine).abs() / fine < 0.05);
    }
}
