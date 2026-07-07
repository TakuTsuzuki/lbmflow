//! Oxygen uptake and reaction-source hooks.

use crate::oxygen::{clip_negative_concentrations, OxygenDiagnostics};
use crate::solver::UnsupportedReason;
use serde::{Deserialize, Serialize};

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum OurModel {
    Constant {
        our_kmol_m3_s: f64,
    },
    Monod {
        our_max: f64,
        ks: f64,
        c_ref: f64,
    },
    CellDensityScaled {
        specific_our: f64,
        cell_density_field: String,
    },
}

#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct OurLedger {
    pub cumulative_kmol_consumed: f64,
    pub last_step_kmol_consumed: f64,
}

#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct ReactionError {
    pub code: &'static str,
    pub message: String,
    pub reason: UnsupportedReason,
}

impl ReactionError {
    fn invalid(message: impl Into<String>, detail: impl Into<String>) -> Self {
        Self {
            code: "reaction_out_of_validity_range",
            message: message.into(),
            reason: UnsupportedReason::OutOfValidityRange {
                detail: detail.into(),
            },
        }
    }
}

impl std::fmt::Display for ReactionError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{} ({:?})", self.message, self.reason)
    }
}

impl std::error::Error for ReactionError {}

pub type ReactionResult<T> = Result<T, ReactionError>;

pub fn oxygen_uptake_rate(
    model: &OurModel,
    c_liquid: f64,
    cell_density: Option<f64>,
) -> ReactionResult<f64> {
    if !(c_liquid.is_finite() && c_liquid >= 0.0) {
        return Err(ReactionError::invalid(
            "oxygen concentration must be finite and non-negative",
            "C_L must be >= 0 before OUR evaluation",
        ));
    }
    match model {
        OurModel::Constant { our_kmol_m3_s } => {
            validate_nonnegative("our_kmol_m3_s", *our_kmol_m3_s)
        }
        OurModel::Monod { our_max, ks, c_ref } => {
            validate_nonnegative("our_max", *our_max)?;
            if !(ks.is_finite() && *ks > 0.0 && c_ref.is_finite() && *c_ref > 0.0) {
                return Err(ReactionError::invalid(
                    "Monod parameters must be finite and positive where dimensional",
                    "ks and c_ref must be > 0",
                ));
            }
            let c_scaled = c_liquid / *c_ref;
            Ok(*our_max * c_scaled / (*ks + c_scaled))
        }
        OurModel::CellDensityScaled { specific_our, .. } => {
            validate_nonnegative("specific_our", *specific_our)?;
            let density = cell_density.ok_or_else(|| {
                ReactionError::invalid(
                    "cell-density-scaled OUR requires a cell density field value",
                    "cell_density must be supplied for CellDensityScaled",
                )
            })?;
            validate_nonnegative("cell_density", density)?;
            Ok(*specific_our * density)
        }
    }
}

pub fn apply_oxygen_reaction_source(
    c_liquid: &mut [f64],
    cell_volumes_m3: &[f64],
    dt_s: f64,
    model: &OurModel,
    cell_density: Option<&[f64]>,
    ledger: &mut OurLedger,
) -> ReactionResult<OxygenDiagnostics> {
    if c_liquid.len() != cell_volumes_m3.len()
        || cell_density.is_some_and(|density| density.len() != c_liquid.len())
    {
        return Err(ReactionError::invalid(
            "reaction-source arrays must have equal length",
            "C_L, cell volumes, and optional cell density must match",
        ));
    }
    if !(dt_s.is_finite() && dt_s >= 0.0) {
        return Err(ReactionError::invalid(
            "reaction time step must be finite and non-negative",
            "dt_s must be >= 0",
        ));
    }
    let mut consumed = 0.0;
    for i in 0..c_liquid.len() {
        if !(cell_volumes_m3[i].is_finite() && cell_volumes_m3[i] > 0.0) {
            return Err(ReactionError::invalid(
                "cell volume must be finite and positive",
                "cell_volumes_m3 entries must be > 0",
            ));
        }
        let density = cell_density.map(|d| d[i]);
        let rate = oxygen_uptake_rate(model, c_liquid[i], density)?;
        c_liquid[i] -= rate * dt_s;
        consumed += rate * cell_volumes_m3[i] * dt_s;
    }
    let diag = clip_negative_concentrations(c_liquid);
    ledger.last_step_kmol_consumed = consumed;
    ledger.cumulative_kmol_consumed += consumed;
    Ok(diag)
}

fn validate_nonnegative(name: &'static str, value: f64) -> ReactionResult<f64> {
    if value.is_finite() && value >= 0.0 {
        Ok(value)
    } else {
        Err(ReactionError::invalid(
            format!("{name} must be finite and non-negative"),
            format!("{name} must be >= 0"),
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn constant_our_gives_linear_depletion_closed_uniform_field() {
        let mut c = vec![1.0; 4];
        let volume = vec![1.0; 4];
        let mut ledger = OurLedger::default();
        apply_oxygen_reaction_source(
            &mut c,
            &volume,
            2.0,
            &OurModel::Constant { our_kmol_m3_s: 0.1 },
            None,
            &mut ledger,
        )
        .unwrap();
        for value in c {
            assert!((value - 0.8).abs() < 1.0e-12);
        }
        assert!((ledger.last_step_kmol_consumed - 0.8).abs() < 1.0e-12);
    }

    #[test]
    fn zero_our_has_no_effect() {
        let mut c = vec![0.3, 0.4];
        let volume = vec![1.0, 1.0];
        let mut ledger = OurLedger::default();
        apply_oxygen_reaction_source(
            &mut c,
            &volume,
            1.0,
            &OurModel::Constant { our_kmol_m3_s: 0.0 },
            None,
            &mut ledger,
        )
        .unwrap();
        assert_eq!(c, vec![0.3, 0.4]);
        assert_eq!(ledger.last_step_kmol_consumed, 0.0);
    }

    #[test]
    fn negative_concentration_guard_fires() {
        let mut c = vec![0.1];
        let volume = vec![1.0];
        let mut ledger = OurLedger::default();
        let diag = apply_oxygen_reaction_source(
            &mut c,
            &volume,
            2.0,
            &OurModel::Constant { our_kmol_m3_s: 0.2 },
            None,
            &mut ledger,
        )
        .unwrap();
        assert_eq!(c, vec![0.0]);
        assert_eq!(diag.clipped_cells, 1);
    }

    #[test]
    fn source_ledger_matches_consumption_integral() {
        let mut c = vec![1.0, 1.0];
        let volume = vec![2.0, 3.0];
        let mut ledger = OurLedger::default();
        apply_oxygen_reaction_source(
            &mut c,
            &volume,
            0.5,
            &OurModel::Constant { our_kmol_m3_s: 0.2 },
            None,
            &mut ledger,
        )
        .unwrap();
        assert!((ledger.last_step_kmol_consumed - 0.5).abs() < 1.0e-12);
    }
}
