//! Runtime divergence guards for optional physics paths.

use crate::phase_field::PhaseFieldDiagnostics;

#[derive(Clone, Debug, PartialEq)]
pub enum DivergenceError {
    Nan {
        step: u64,
    },
    PhiOutOfBounds {
        step: u64,
        min_phi: f64,
        max_phi: f64,
    },
    MassDriftExcessive {
        step: u64,
        initial_total_phi: f64,
        current_total_phi: f64,
        relative_drift: f64,
    },
}

impl std::fmt::Display for DivergenceError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Nan { step } => write!(f, "phase field contains NaN at step {step}"),
            Self::PhiOutOfBounds {
                step,
                min_phi,
                max_phi,
            } => write!(
                f,
                "phase field out of bounds at step {step}: min_phi={min_phi}, max_phi={max_phi}"
            ),
            Self::MassDriftExcessive {
                step,
                initial_total_phi,
                current_total_phi,
                relative_drift,
            } => write!(
                f,
                "phase mass drift excessive at step {step}: initial={initial_total_phi}, current={current_total_phi}, rel={relative_drift}"
            ),
        }
    }
}

impl std::error::Error for DivergenceError {}

pub type PhaseDiag = PhaseFieldDiagnostics;
