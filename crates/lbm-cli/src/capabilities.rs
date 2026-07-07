use anyhow::Result;
use serde::{Deserialize, Serialize};

pub use lbm_core::solver::UnsupportedReason;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum CapabilityStatus {
    Unsupported,
    Experimental,
    Engineering,
    EvidenceBlocked,
    EvidenceReady,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum CapabilityTier {
    Screening,
    Engineering,
    Evidence,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct CapabilityEntry {
    pub id: &'static str,
    pub label: &'static str,
    pub status: CapabilityStatus,
    pub tier_ceiling: CapabilityTier,
    pub reason: Option<UnsupportedReason>,
    pub docs: &'static str,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct CapabilityRegistry {
    capabilities: Vec<CapabilityEntry>,
}

impl CapabilityRegistry {
    pub fn new() -> Self {
        Self {
            capabilities: vec![
                experimental(
                    "single_phase_stirred_tank",
                    "Single-phase stirred tank",
                    "BCFD-030 runner path implemented; VB-01 not yet green",
                    "docs/PLAN.md#bcfd-030",
                ),
                experimental(
                    "rotating_ibm",
                    "Rotating IBM impeller",
                    "BCFD-030 integration implemented; stirred-tank validation pending",
                    "docs/PLAN.md#bcfd-021",
                ),
                experimental(
                    "passive_scalar",
                    "Passive scalar transport",
                    "BCFD-034 ADE path implemented; VB-02 not yet green",
                    "docs/PLAN.md#bcfd-034",
                ),
                unsupported(
                    "phase_field_vof",
                    "Phase-field VOF",
                    "Conservative Allen-Cahn path is not yet implemented",
                    "docs/PLAN.md#bcfd-040",
                ),
                unsupported(
                    "oxygen_kla",
                    "Oxygen transport and kLa",
                    "Oxygen scalar and kLa QOI are not yet implemented",
                    "docs/PLAN.md#bcfd-050",
                ),
                unsupported(
                    "point_bubbles",
                    "Point bubbles",
                    "Point-bubble entity store is not yet implemented",
                    "docs/PLAN.md#bcfd-070",
                ),
                unsupported(
                    "pbm",
                    "Population balance model",
                    "PBM bins and kernels are not yet implemented",
                    "docs/PLAN.md#bcfd-073",
                ),
                unsupported(
                    "cell_exposure",
                    "Cell and microcarrier exposure",
                    "Cell tracer and exposure QOIs are not yet implemented",
                    "docs/PLAN.md#bcfd-060",
                ),
                CapabilityEntry {
                    id: "evidence_tier_report",
                    label: "Evidence-tier report",
                    status: CapabilityStatus::Unsupported,
                    tier_ceiling: CapabilityTier::Evidence,
                    reason: Some(UnsupportedReason::EvidenceGateFailed {
                        missing: vec![
                            "validation matrix pass".to_string(),
                            "calibration/holdout separation".to_string(),
                            "mesh/time-step sensitivity".to_string(),
                            "QOI uncertainty interval".to_string(),
                            "limitation report".to_string(),
                        ],
                    }),
                    docs: "docs/PLAN.md#bcfd-091",
                },
            ],
        }
    }

    #[allow(dead_code)]
    pub fn get(&self, id: &str) -> Option<&CapabilityEntry> {
        self.capabilities.iter().find(|entry| entry.id == id)
    }

    pub fn iter(&self) -> impl Iterator<Item = &CapabilityEntry> {
        self.capabilities.iter()
    }

    pub fn to_json(&self) -> serde_json::Value {
        serde_json::to_value(self).expect("capability registry must serialize")
    }
}

impl Default for CapabilityRegistry {
    fn default() -> Self {
        Self::new()
    }
}

fn unsupported(
    id: &'static str,
    label: &'static str,
    detail: &'static str,
    docs: &'static str,
) -> CapabilityEntry {
    CapabilityEntry {
        id,
        label,
        status: CapabilityStatus::Unsupported,
        tier_ceiling: CapabilityTier::Screening,
        reason: Some(UnsupportedReason::OutOfValidityRange {
            detail: detail.to_string(),
        }),
        docs,
    }
}

fn experimental(
    id: &'static str,
    label: &'static str,
    detail: &'static str,
    docs: &'static str,
) -> CapabilityEntry {
    CapabilityEntry {
        id,
        label,
        status: CapabilityStatus::Experimental,
        tier_ceiling: CapabilityTier::Screening,
        reason: Some(UnsupportedReason::OutOfValidityRange {
            detail: detail.to_string(),
        }),
        docs,
    }
}

pub fn run(json: bool) -> Result<()> {
    let registry = CapabilityRegistry::new();
    if json {
        println!("{}", serde_json::to_string_pretty(&registry.to_json())?);
    } else {
        print_human(&registry);
    }
    Ok(())
}

fn print_human(registry: &CapabilityRegistry) {
    println!("LBMFlow bioprocess capability registry");
    println!(
        "{:<30} {:<14} {:<12} {:<24} Docs",
        "id", "status", "tier ceiling", "reason"
    );
    for entry in registry.iter() {
        let reason = match &entry.reason {
            Some(UnsupportedReason::NotImplemented) => "not implemented".to_string(),
            Some(UnsupportedReason::OutOfValidityRange { detail }) => detail.clone(),
            Some(UnsupportedReason::MissingDependency { depends_on }) => {
                format!("missing {depends_on}")
            }
            Some(UnsupportedReason::EvidenceGateFailed { missing }) => {
                format!("evidence gate missing {}", missing.join(", "))
            }
            Some(UnsupportedReason::DemoOnly { rationale }) => rationale.clone(),
            None => "-".to_string(),
        };
        println!(
            "{:<30} {:<14?} {:<12?} {:<24} {}",
            entry.id,
            entry.status,
            entry.tier_ceiling,
            truncate(&reason, 24),
            entry.docs
        );
    }
}

fn truncate(text: &str, width: usize) -> String {
    if text.len() <= width {
        text.to_string()
    } else {
        format!("{}...", &text[..width.saturating_sub(3)])
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn capability_registry_json_round_trips() {
        let registry = CapabilityRegistry::new();
        let value = registry.to_json();
        let encoded = serde_json::to_string(&value).expect("registry JSON should encode");
        let decoded: serde_json::Value =
            serde_json::from_str(&encoded).expect("registry JSON should decode");
        assert_eq!(decoded, value);
        assert!(registry.get("single_phase_stirred_tank").is_some());
    }
}
