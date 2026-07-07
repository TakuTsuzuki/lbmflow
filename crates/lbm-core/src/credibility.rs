//! Calibration, holdout, and evidence-claim registry checks.

use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};

use crate::qoi::QoiInterval;

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub struct CalibrationDataset {
    pub id: String,
    pub qoi: String,
    pub source: String,
    pub date: String,
    pub scale: String,
    pub operating_condition: String,
    pub measurement_uncertainty: f64,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub struct HoldoutDataset {
    pub id: String,
    pub qoi: String,
    pub source: String,
    pub date: String,
    pub scale: String,
    pub operating_condition: String,
    pub measurement_uncertainty: f64,
}

#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub struct DatasetRegistry {
    pub calibration: HashMap<String, CalibrationDataset>,
    pub holdout: HashMap<String, HoldoutDataset>,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub struct SensitivityRecord {
    pub additional_cases: usize,
    pub variation_fraction: f64,
    pub summary: String,
}

impl SensitivityRecord {
    pub fn passes_evidence_band(&self) -> bool {
        self.additional_cases >= 2
            && self.variation_fraction.is_finite()
            && self.variation_fraction <= 0.05
    }
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub struct LimitationReport {
    pub manifest_path: String,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub struct EvidenceGateInput {
    pub qoi_id: String,
    pub vb_id: String,
    pub vb_engineering_green: bool,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub vb_reason: Option<String>,
    pub datasets: DatasetRegistry,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub uq_interval: Option<QoiInterval>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub mesh_sensitivity: Option<SensitivityRecord>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub time_step_sensitivity: Option<SensitivityRecord>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub limitation_report: Option<LimitationReport>,
}

pub struct EvidenceGate;

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(tag = "status", rename_all = "snake_case")]
pub enum EvidenceGateResult {
    Ready {
        qoi_id: String,
        calibration_ids: Vec<String>,
        holdout_ids: Vec<String>,
        uq_interval: QoiInterval,
        sensitivity_summary: String,
    },
    Blocked {
        qoi_id: String,
        missing: Vec<MissingArtefact>,
    },
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum MissingArtefact {
    ValidationMatrixFail { vb_id: String, reason: String },
    NoCalibrationDataset,
    NoHoldoutDataset,
    DatasetReuseConflict { id: String },
    NoUqInterval,
    NoMeshSensitivity,
    NoTimeStepSensitivity,
    NoLimitationReport,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum CredibilityError {
    EvidenceGateFailed { qoi: String, missing: Vec<String> },
}

impl std::fmt::Display for CredibilityError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::EvidenceGateFailed { qoi, missing } => {
                write!(f, "evidence gate failed for {qoi}: {}", missing.join(", "))
            }
        }
    }
}

impl std::error::Error for CredibilityError {}

impl DatasetRegistry {
    pub fn validate_no_reuse_for_qoi(&self, qoi: &str) -> Result<(), CredibilityError> {
        let calibration: HashSet<&str> = self
            .calibration
            .values()
            .filter(|dataset| dataset.qoi == qoi)
            .map(|dataset| dataset.id.as_str())
            .collect();
        let reused = self
            .holdout
            .values()
            .any(|dataset| dataset.qoi == qoi && calibration.contains(dataset.id.as_str()));
        if reused {
            return Err(CredibilityError::EvidenceGateFailed {
                qoi: qoi.to_string(),
                missing: vec!["dataset_reuse_conflict".to_string()],
            });
        }
        Ok(())
    }

    pub fn require_holdout_for_evidence(&self, qoi: &str) -> Result<(), CredibilityError> {
        self.validate_no_reuse_for_qoi(qoi)?;
        if self.holdout.values().any(|dataset| dataset.qoi == qoi) {
            Ok(())
        } else {
            Err(CredibilityError::EvidenceGateFailed {
                qoi: qoi.to_string(),
                missing: vec!["holdout_dataset".to_string()],
            })
        }
    }

    pub fn engineering_calibration_only_warnings(&self, qoi: &str) -> Vec<String> {
        let calibration_ids: Vec<_> = self
            .calibration
            .values()
            .filter(|dataset| dataset.qoi == qoi)
            .map(|dataset| dataset.id.clone())
            .collect();
        let has_holdout = self.holdout.values().any(|dataset| dataset.qoi == qoi);
        if calibration_ids.is_empty() || has_holdout {
            Vec::new()
        } else {
            vec![format!(
                "{qoi} calibrated to {}, not validated against holdout",
                calibration_ids.join(",")
            )]
        }
    }
}

impl EvidenceGate {
    pub fn evaluate(input: EvidenceGateInput) -> EvidenceGateResult {
        let qoi = input.qoi_id.clone();
        let mut missing = Vec::new();

        if !input.vb_engineering_green {
            missing.push(MissingArtefact::ValidationMatrixFail {
                vb_id: input.vb_id.clone(),
                reason: input.vb_reason.clone().unwrap_or_else(|| {
                    "validation matrix entry is not Engineering GREEN".to_string()
                }),
            });
        }

        let calibration_ids = dataset_ids_for_qoi(
            input
                .datasets
                .calibration
                .values()
                .map(DatasetView::Calibration),
            &qoi,
        );
        let holdout_ids = dataset_ids_for_qoi(
            input.datasets.holdout.values().map(DatasetView::Holdout),
            &qoi,
        );
        if calibration_ids.is_empty() {
            missing.push(MissingArtefact::NoCalibrationDataset);
        }
        if holdout_ids.is_empty() {
            missing.push(MissingArtefact::NoHoldoutDataset);
        }
        for id in calibration_ids
            .iter()
            .filter(|id| holdout_ids.iter().any(|holdout| holdout == *id))
        {
            missing.push(MissingArtefact::DatasetReuseConflict { id: id.clone() });
        }

        let uq_interval = match input.uq_interval.clone() {
            Some(interval) => Some(interval),
            None => {
                missing.push(MissingArtefact::NoUqInterval);
                None
            }
        };

        let mesh_summary = match &input.mesh_sensitivity {
            Some(record) if record.passes_evidence_band() => Some(record.summary.clone()),
            _ => {
                missing.push(MissingArtefact::NoMeshSensitivity);
                None
            }
        };
        let time_summary = match &input.time_step_sensitivity {
            Some(record) if record.passes_evidence_band() => Some(record.summary.clone()),
            _ => {
                missing.push(MissingArtefact::NoTimeStepSensitivity);
                None
            }
        };
        if input.limitation_report.is_none() {
            missing.push(MissingArtefact::NoLimitationReport);
        }

        if missing.is_empty() {
            EvidenceGateResult::Ready {
                qoi_id: qoi,
                calibration_ids,
                holdout_ids,
                uq_interval: uq_interval.expect("missing is empty only when UQ interval exists"),
                sensitivity_summary: format!(
                    "mesh: {}; time_step: {}",
                    mesh_summary.unwrap_or_default(),
                    time_summary.unwrap_or_default()
                ),
            }
        } else {
            EvidenceGateResult::Blocked {
                qoi_id: qoi,
                missing,
            }
        }
    }
}

enum DatasetView<'a> {
    Calibration(&'a CalibrationDataset),
    Holdout(&'a HoldoutDataset),
}

impl DatasetView<'_> {
    fn qoi(&self) -> &str {
        match self {
            Self::Calibration(dataset) => &dataset.qoi,
            Self::Holdout(dataset) => &dataset.qoi,
        }
    }

    fn id(&self) -> String {
        match self {
            Self::Calibration(dataset) => dataset.id.clone(),
            Self::Holdout(dataset) => dataset.id.clone(),
        }
    }
}

fn dataset_ids_for_qoi<'a>(
    datasets: impl Iterator<Item = DatasetView<'a>>,
    qoi: &str,
) -> Vec<String> {
    let mut ids: Vec<_> = datasets
        .filter(|dataset| dataset.qoi() == qoi)
        .map(|dataset| dataset.id())
        .collect();
    ids.sort();
    ids.dedup();
    ids
}

#[cfg(test)]
mod tests {
    use super::*;

    fn calibration(id: &str, qoi: &str) -> CalibrationDataset {
        CalibrationDataset {
            id: id.to_string(),
            qoi: qoi.to_string(),
            source: "fixture".to_string(),
            date: "2026-07-07".to_string(),
            scale: "bench".to_string(),
            operating_condition: "N=1/s".to_string(),
            measurement_uncertainty: 0.1,
        }
    }

    fn holdout(id: &str, qoi: &str) -> HoldoutDataset {
        HoldoutDataset {
            id: id.to_string(),
            qoi: qoi.to_string(),
            source: "fixture".to_string(),
            date: "2026-07-07".to_string(),
            scale: "pilot".to_string(),
            operating_condition: "N=1/s".to_string(),
            measurement_uncertainty: 0.1,
        }
    }

    #[test]
    fn same_id_rejected_for_same_qoi() {
        let mut registry = DatasetRegistry::default();
        registry
            .calibration
            .insert("d1".to_string(), calibration("d1", "kla"));
        registry
            .holdout
            .insert("d1".to_string(), holdout("d1", "kla"));
        let err = registry.validate_no_reuse_for_qoi("kla").unwrap_err();
        assert_eq!(
            err,
            CredibilityError::EvidenceGateFailed {
                qoi: "kla".to_string(),
                missing: vec!["dataset_reuse_conflict".to_string()]
            }
        );
    }

    #[test]
    fn missing_holdout_rejected_at_evidence_tier() {
        let mut registry = DatasetRegistry::default();
        registry
            .calibration
            .insert("d1".to_string(), calibration("d1", "kla"));
        let err = registry.require_holdout_for_evidence("kla").unwrap_err();
        assert_eq!(
            err,
            CredibilityError::EvidenceGateFailed {
                qoi: "kla".to_string(),
                missing: vec!["holdout_dataset".to_string()]
            }
        );
    }

    #[test]
    fn engineering_tier_calibration_only_allowed_with_warning() {
        let mut registry = DatasetRegistry::default();
        registry
            .calibration
            .insert("d1".to_string(), calibration("d1", "kla"));
        assert!(registry.validate_no_reuse_for_qoi("kla").is_ok());
        let warnings = registry.engineering_calibration_only_warnings("kla");
        assert_eq!(warnings.len(), 1);
        assert!(warnings[0].contains("not validated against holdout"));
    }

    fn sensitivity(summary: &str) -> SensitivityRecord {
        SensitivityRecord {
            additional_cases: 2,
            variation_fraction: 0.04,
            summary: summary.to_string(),
        }
    }

    fn interval() -> QoiInterval {
        QoiInterval {
            q_hat: 1.0,
            q_lo: 0.9,
            q_hi: 1.1,
            method: "fixture_uq".to_string(),
        }
    }

    fn gate_input() -> EvidenceGateInput {
        let mut datasets = DatasetRegistry::default();
        datasets
            .calibration
            .insert("cal-1".to_string(), calibration("cal-1", "kla"));
        datasets
            .holdout
            .insert("hold-1".to_string(), holdout("hold-1", "kla"));
        EvidenceGateInput {
            qoi_id: "kla".to_string(),
            vb_id: "VB-06".to_string(),
            vb_engineering_green: true,
            vb_reason: None,
            datasets,
            uq_interval: Some(interval()),
            mesh_sensitivity: Some(sensitivity("three grids, finest-pair variation 4%")),
            time_step_sensitivity: Some(sensitivity("three dt values, finest-pair variation 4%")),
            limitation_report: Some(LimitationReport {
                manifest_path: "limitations.md".to_string(),
            }),
        }
    }

    #[test]
    fn evidence_tier_without_holdout_is_blocked() {
        let mut input = gate_input();
        input.datasets.holdout.clear();
        let result = EvidenceGate::evaluate(input);
        match result {
            EvidenceGateResult::Blocked { missing, .. } => {
                assert!(missing.contains(&MissingArtefact::NoHoldoutDataset));
            }
            other => panic!("expected blocked without holdout, got {other:?}"),
        }
    }

    #[test]
    fn evidence_tier_without_mesh_sensitivity_is_blocked() {
        let mut input = gate_input();
        input.mesh_sensitivity = None;
        let result = EvidenceGate::evaluate(input);
        match result {
            EvidenceGateResult::Blocked { missing, .. } => {
                assert!(missing.contains(&MissingArtefact::NoMeshSensitivity));
            }
            other => panic!("expected blocked without mesh sensitivity, got {other:?}"),
        }
    }

    #[test]
    fn same_id_calibration_and_holdout_is_blocked() {
        let mut input = gate_input();
        input.datasets.holdout.clear();
        input
            .datasets
            .holdout
            .insert("cal-1".to_string(), holdout("cal-1", "kla"));
        let result = EvidenceGate::evaluate(input);
        match result {
            EvidenceGateResult::Blocked { missing, .. } => {
                assert!(missing.contains(&MissingArtefact::DatasetReuseConflict {
                    id: "cal-1".to_string()
                }));
            }
            other => panic!("expected blocked on dataset reuse, got {other:?}"),
        }
    }

    #[test]
    fn full_artefact_set_is_evidence_ready() {
        let result = EvidenceGate::evaluate(gate_input());
        match result {
            EvidenceGateResult::Ready {
                calibration_ids,
                holdout_ids,
                uq_interval,
                sensitivity_summary,
                ..
            } => {
                assert_eq!(calibration_ids, vec!["cal-1"]);
                assert_eq!(holdout_ids, vec!["hold-1"]);
                assert_eq!(uq_interval.q_hat, 1.0);
                assert!(sensitivity_summary.contains("mesh"));
            }
            other => panic!("expected ready with full artefact set, got {other:?}"),
        }
    }
}
