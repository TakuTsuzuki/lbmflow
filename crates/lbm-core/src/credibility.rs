//! Calibration, holdout, and evidence-claim registry checks.

use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};

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
}
