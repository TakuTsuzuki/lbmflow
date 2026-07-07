//! Hybrid resolved-interface plus point-bubble gas bookkeeping.

use crate::bubbles::{BubbleError, BubbleSet};
use crate::kla::KlaProvenance;
use crate::solver::UnsupportedReason;
use serde::Serialize;

pub const RESOLVED_INTERFACE_PHI_CENTER: f64 = 0.5;
pub const RESOLVED_INTERFACE_HALF_WIDTH: f64 = 0.1;

#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct HybridGasReport {
    pub epsilon_g_resolved: f64,
    pub epsilon_g_bubble: f64,
    pub epsilon_g_total: f64,
    pub a_resolved_1_m: f64,
    pub a_bubble_1_m: f64,
    pub metadata: HybridGasMetadata,
}

#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct HybridGasMetadata {
    pub method: &'static str,
    pub double_count_policy: &'static str,
    pub validation_tier: &'static str,
}

pub fn hybrid_gas_bookkeeping(
    phi_liquid_fraction: &[f64],
    resolved_area_density_1_m: &[f64],
    bubbles: &BubbleSet,
    dims: [usize; 3],
    dx_m: f64,
) -> HybridGasResult<HybridGasReport> {
    let n = dims[0] * dims[1] * dims[2];
    if phi_liquid_fraction.len() != n || resolved_area_density_1_m.len() != n {
        return Err(HybridGasError::out_of_validity_range(
            "hybrid gas fields must match grid dimensions",
            "phi_liquid_fraction and resolved_area_density_1_m length must equal cell count",
        ));
    }
    if !(dx_m.is_finite() && dx_m > 0.0) {
        return Err(HybridGasError::out_of_validity_range(
            "grid spacing must be finite and positive",
            "dx_m must be > 0",
        ));
    }
    let cell_volume = dx_m * dx_m * dx_m;
    let total_volume = cell_volume * n as f64;
    let mut resolved_volume = 0.0;
    let mut bubble_volume = 0.0;
    let mut a_resolved = 0.0;
    let mut a_bubble = 0.0;

    for (&phi, &area) in phi_liquid_fraction
        .iter()
        .zip(resolved_area_density_1_m.iter())
    {
        if !(phi.is_finite() && phi >= 0.0 && phi <= 1.0 && area.is_finite() && area >= 0.0) {
            return Err(HybridGasError::out_of_validity_range(
                "resolved gas fields must be finite and bounded",
                "phi must be in [0,1] and area density must be >= 0",
            ));
        }
        resolved_volume += (1.0 - phi) * cell_volume;
        a_resolved += area * cell_volume;
    }

    for bubble in &bubbles.bubbles {
        if let Some(i) = crate::bubbles::cell_index_for_position(dims, dx_m, bubble.position) {
            if is_resolved_cell(phi_liquid_fraction[i]) {
                continue;
            }
            bubble_volume += bubble.gas_volume_m3;
            a_bubble += std::f64::consts::PI * bubble.diameter_m * bubble.diameter_m;
        }
    }

    Ok(HybridGasReport {
        epsilon_g_resolved: resolved_volume / total_volume,
        epsilon_g_bubble: bubble_volume / total_volume,
        epsilon_g_total: (resolved_volume + bubble_volume) / total_volume,
        a_resolved_1_m: a_resolved / total_volume,
        a_bubble_1_m: a_bubble / total_volume,
        metadata: HybridGasMetadata {
            method: "hybrid_resolved_point_bubble",
            double_count_policy: "ignore_point_bubbles_inside_resolved_interface_cells",
            validation_tier: "experimental",
        },
    })
}

pub fn is_resolved_cell(phi_liquid_fraction: f64) -> bool {
    let diff = (phi_liquid_fraction - RESOLVED_INTERFACE_PHI_CENTER).abs();
    diff < RESOLVED_INTERFACE_HALF_WIDTH
}

pub fn reject_hybrid_evidence_tier() -> HybridGasResult<()> {
    Err(HybridGasError {
        code: "hybrid_evidence_gate_rejected",
        message: "hybrid gas bookkeeping requires VB-05 and VB-06 Engineering GREEN".to_string(),
        reason: UnsupportedReason::EvidenceGateFailed {
            missing: vec![
                "VB-05 gas holdup Engineering GREEN".to_string(),
                "VB-06 kLa Engineering GREEN".to_string(),
            ],
        },
    })
}

pub fn hybrid_kla_provenance() -> KlaProvenance {
    KlaProvenance {
        method: "hybrid_resolved_point_bubble",
        k_l_source: "hybrid",
        d32_source: "pbm",
        units: "1/s",
    }
}

#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct HybridGasError {
    pub code: &'static str,
    pub message: String,
    pub reason: UnsupportedReason,
}

impl HybridGasError {
    pub fn out_of_validity_range(message: impl Into<String>, detail: impl Into<String>) -> Self {
        Self {
            code: "hybrid_gas_out_of_validity_range",
            message: message.into(),
            reason: UnsupportedReason::OutOfValidityRange {
                detail: detail.into(),
            },
        }
    }
}

impl From<BubbleError> for HybridGasError {
    fn from(value: BubbleError) -> Self {
        Self {
            code: "bubble_error",
            message: value.message,
            reason: value.reason,
        }
    }
}

impl std::fmt::Display for HybridGasError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{} ({:?})", self.message, self.reason)
    }
}

impl std::error::Error for HybridGasError {}

pub type HybridGasResult<T> = Result<T, HybridGasError>;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bubbles::{Bubble, BubbleSet};

    #[test]
    fn synthetic_phi_and_bubbles_report_expected_totals() {
        let mut bubbles = BubbleSet::new();
        bubbles
            .bubbles
            .push(Bubble::new([1.5, 0.5, 0.5], [0.0; 3], 0.1, 10.0, 1).unwrap());
        let phi = vec![1.0, 1.0];
        let area = vec![0.0, 0.0];
        let report = hybrid_gas_bookkeeping(&phi, &area, &bubbles, [2, 1, 1], 1.0).unwrap();
        let expected = bubbles.bubbles[0].gas_volume_m3 / 2.0;
        assert!((report.epsilon_g_total - expected).abs() < 1.0e-12);
        assert_eq!(report.epsilon_g_resolved, 0.0);
    }

    #[test]
    fn double_count_excludes_bubbles_in_resolved_interface_cell() {
        let mut bubbles = BubbleSet::new();
        bubbles
            .bubbles
            .push(Bubble::new([0.5, 0.5, 0.5], [0.0; 3], 0.1, 10.0, 1).unwrap());
        let phi = vec![0.5];
        let area = vec![2.0];
        let report = hybrid_gas_bookkeeping(&phi, &area, &bubbles, [1, 1, 1], 1.0).unwrap();
        assert_eq!(report.epsilon_g_bubble, 0.0);
        assert_eq!(report.epsilon_g_resolved, 0.5);
        assert_eq!(report.a_resolved_1_m, 2.0);
    }

    #[test]
    fn method_metadata_and_evidence_rejection_are_explicit() {
        let bubbles = BubbleSet::new();
        let report = hybrid_gas_bookkeeping(&[1.0], &[0.0], &bubbles, [1, 1, 1], 1.0).unwrap();
        assert_eq!(report.metadata.method, "hybrid_resolved_point_bubble");
        assert_eq!(
            report.metadata.double_count_policy,
            "ignore_point_bubbles_inside_resolved_interface_cells"
        );
        assert_eq!(
            reject_hybrid_evidence_tier().unwrap_err().code,
            "hybrid_evidence_gate_rejected"
        );
    }
}
