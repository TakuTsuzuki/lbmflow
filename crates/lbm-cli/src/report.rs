use crate::manifest::MANIFEST_PATH;
use anyhow::Result;
use lbm_core::credibility::{EvidenceGate, EvidenceGateInput, EvidenceGateResult, MissingArtefact};
use lbm_core::qoi::{CapabilityStatus, QoiBundle, QoiScalar};
use serde::Serialize;
use serde_json::Value;
use std::fmt::Write as _;
use std::fs;
use std::path::{Path, PathBuf};

pub const REPORT_PATH: &str = "report.md";

#[derive(Clone, Debug, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub struct ReportError {
    pub code: &'static str,
    pub message: String,
    pub path: String,
}

impl ReportError {
    fn missing_qoi(path: &Path) -> Self {
        Self {
            code: "missing_qoi_json",
            message: "bioprocess report requires qoi.json".to_string(),
            path: path.display().to_string(),
        }
    }
}

impl std::fmt::Display for ReportError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let text = serde_json::to_string(self).map_err(|_| std::fmt::Error)?;
        write!(f, "{text}")
    }
}

impl std::error::Error for ReportError {}

pub fn generate_report(run_dir: &Path) -> Result<PathBuf> {
    let qoi_path = run_dir.join(lbm_scenario::QOI_BUNDLE_JSON);
    if !qoi_path.exists() {
        return Err(ReportError::missing_qoi(&qoi_path).into());
    }
    let qoi: QoiBundle = serde_json::from_slice(&fs::read(&qoi_path)?)?;
    let manifest_path = run_dir.join(MANIFEST_PATH);
    let manifest: Value = serde_json::from_slice(&fs::read(&manifest_path)?)?;
    let scenario_path = run_dir.join("scenario.json");
    let scenario: Option<Value> = fs::read(&scenario_path)
        .ok()
        .and_then(|bytes| serde_json::from_slice(&bytes).ok());

    let mut out = String::new();
    let gate_results = evidence_gate_results(&qoi, &manifest);
    if !evidence_gate_ready(&gate_results) {
        out.push_str("# NOT EVIDENCE-GRADE\n\n");
    }
    out.push_str("# Bioprocess CFD Report\n\n");
    out.push_str("## Intended use\n\n");
    out.push_str(
        "LBMFlow is a bioprocess-specific CFD core for stirred-tank cell-culture / bioreactor process design: hydrodynamics, mixing, gas-liquid oxygen transfer, shear exposure, cell or microcarrier exposure, and scale-up operating-window evaluation.\n\n",
    );
    out.push_str("## Forbidden use\n\n");
    out.push_str(
        "Do not use this report for GMP / CMC filing claims unless the associated QOIs pass the evidence-tier gate. Do not claim validated kLa without calibration and independent holdout data. Do not use Shan-Chen gas-liquid output for production decisions. Do not make decisions from QOIs missing provenance metadata or from max shear/exposure alone.\n\n",
    );
    out.push_str("## Scenario summary\n\n");
    write_scenario_summary(&mut out, scenario.as_ref())?;
    out.push_str("\n## Unit feasibility\n\n");
    write_json_or_placeholder(
        &mut out,
        manifest.get("unitReport"),
        "No UnitReport found in manifest.",
    )?;
    out.push_str("\n## Active models\n\n");
    write_list(&mut out, manifest.get("activeModels"));
    out.push_str("\n## QOI summary\n\n");
    write_qoi_summary(&mut out, &qoi)?;
    out.push_str("\n## Validation status\n\n");
    out.push_str("| QOI | Status | Tier |\n|---|---|---|\n");
    for status in &qoi.validation_status {
        writeln!(
            out,
            "| {} | {:?} | {:?} |",
            status.qoi, status.status, status.tier
        )?;
    }
    out.push_str("\n## Evidence gate\n\n");
    write_evidence_gate_summary(&mut out, &gate_results)?;
    out.push_str("\n## Limitations\n\n");
    out.push_str("See [docs/LIMITATIONS.md](../../../docs/LIMITATIONS.md).\n\n");
    out.push_str("## Provenance\n\n");
    let scenario_hash = manifest
        .get("scenarioHash")
        .and_then(Value::as_str)
        .unwrap_or("unknown");
    let git_sha = manifest
        .get("provenance")
        .and_then(|p| p.get("gitSha"))
        .and_then(Value::as_str)
        .unwrap_or("not recorded");
    writeln!(out, "- Manifest: `{}`", MANIFEST_PATH)?;
    writeln!(out, "- Scenario: `scenario.json`")?;
    writeln!(out, "- scenario_hash: `{scenario_hash}`")?;
    writeln!(out, "- git_sha: `{git_sha}`")?;

    let report_path = run_dir.join(REPORT_PATH);
    fs::write(&report_path, out)?;
    Ok(report_path)
}

pub fn evidence_check(run_dir: &Path) -> Result<Vec<EvidenceGateResult>> {
    let qoi_path = run_dir.join(lbm_scenario::QOI_BUNDLE_JSON);
    if !qoi_path.exists() {
        return Err(ReportError::missing_qoi(&qoi_path).into());
    }
    let qoi: QoiBundle = serde_json::from_slice(&fs::read(&qoi_path)?)?;
    let manifest_path = run_dir.join(MANIFEST_PATH);
    let manifest: Value = serde_json::from_slice(&fs::read(&manifest_path)?)?;
    Ok(evidence_gate_results(&qoi, &manifest))
}

fn evidence_gate_ready(results: &[EvidenceGateResult]) -> bool {
    !results.is_empty()
        && results
            .iter()
            .all(|result| matches!(result, EvidenceGateResult::Ready { .. }))
}

fn evidence_gate_results(qoi: &QoiBundle, manifest: &Value) -> Vec<EvidenceGateResult> {
    if let Some(inputs) = manifest_evidence_inputs(manifest) {
        return inputs.into_iter().map(EvidenceGate::evaluate).collect();
    }
    qoi.validation_status
        .iter()
        .map(|status| {
            let vb_id = vb_id_for_qoi(&status.qoi).to_string();
            let mut missing = vec![
                MissingArtefact::NoCalibrationDataset,
                MissingArtefact::NoHoldoutDataset,
                MissingArtefact::NoUqInterval,
                MissingArtefact::NoMeshSensitivity,
                MissingArtefact::NoTimeStepSensitivity,
                MissingArtefact::NoLimitationReport,
            ];
            if !matches!(
                status.status,
                CapabilityStatus::Engineering | CapabilityStatus::EvidenceReady
            ) {
                missing.insert(
                    0,
                    MissingArtefact::ValidationMatrixFail {
                        vb_id,
                        reason: format!("QOI status is {:?}", status.status),
                    },
                );
            }
            EvidenceGateResult::Blocked {
                qoi_id: status.qoi.clone(),
                missing,
            }
        })
        .collect()
}

fn manifest_evidence_inputs(manifest: &Value) -> Option<Vec<EvidenceGateInput>> {
    let value = manifest
        .get("evidenceGate")
        .or_else(|| manifest.get("evidence_gate"))?;
    if let Some(array) = value.as_array() {
        let inputs: Vec<EvidenceGateInput> = array
            .iter()
            .filter_map(|entry| serde_json::from_value(entry.clone()).ok())
            .collect();
        return Some(inputs);
    }
    serde_json::from_value(value.clone())
        .ok()
        .map(|input| vec![input])
}

fn vb_id_for_qoi(qoi: &str) -> &'static str {
    match qoi {
        "power" => "VB-01",
        "mixing" => "VB-02",
        "shear" | "shear_rate" => "VB-03",
        "gas" | "gas_holdup" => "VB-05",
        "kla" | "oxygen" => "VB-06",
        "cells" | "cell_exposure" => "VB-07",
        "scaleup" | "scale_up" => "VB-08",
        _ => "VB-UNKNOWN",
    }
}

fn write_evidence_gate_summary(
    out: &mut String,
    results: &[EvidenceGateResult],
) -> std::fmt::Result {
    if results.is_empty() {
        out.push_str("No evidence-gate artefacts were found.\n");
        return Ok(());
    }
    out.push_str("| QOI | Gate | Detail |\n|---|---|---|\n");
    for result in results {
        match result {
            EvidenceGateResult::Ready {
                qoi_id,
                calibration_ids,
                holdout_ids,
                ..
            } => writeln!(
                out,
                "| {qoi_id} | Ready | calibration={} holdout={} |",
                calibration_ids.join(","),
                holdout_ids.join(",")
            )?,
            EvidenceGateResult::Blocked { qoi_id, missing } => writeln!(
                out,
                "| {qoi_id} | Blocked | {} |",
                missing
                    .iter()
                    .map(|item| format!("{item:?}"))
                    .collect::<Vec<_>>()
                    .join("; ")
            )?,
        }
    }
    Ok(())
}

fn write_scenario_summary(out: &mut String, scenario: Option<&Value>) -> std::fmt::Result {
    let Some(scenario) = scenario else {
        out.push_str("Scenario JSON was not found in the run directory.\n");
        return Ok(());
    };
    writeln!(
        out,
        "- Name: `{}`",
        scenario
            .get("name")
            .and_then(Value::as_str)
            .unwrap_or("unknown")
    )?;
    writeln!(
        out,
        "- Version: `{}`",
        scenario
            .get("version")
            .and_then(Value::as_str)
            .unwrap_or("unknown")
    )?;
    writeln!(
        out,
        "- Credibility tier requested: `{}`",
        scenario
            .get("credibility_tier")
            .and_then(Value::as_str)
            .unwrap_or("unknown")
    )
}

fn write_json_or_placeholder(
    out: &mut String,
    value: Option<&Value>,
    placeholder: &str,
) -> std::fmt::Result {
    if let Some(value) = value {
        writeln!(out, "```json")?;
        writeln!(
            out,
            "{}",
            serde_json::to_string_pretty(value).unwrap_or_else(|_| "null".to_string())
        )?;
        writeln!(out, "```")
    } else {
        writeln!(out, "{placeholder}")
    }
}

fn write_list(out: &mut String, value: Option<&Value>) {
    if let Some(models) = value.and_then(Value::as_array) {
        for model in models {
            if let Some(model) = model.as_str() {
                let _ = writeln!(out, "- `{model}`");
            }
        }
    } else {
        out.push_str("- none recorded\n");
    }
}

fn write_qoi_summary(out: &mut String, qoi: &QoiBundle) -> std::fmt::Result {
    out.push_str("| Section | QOI | Value | Units | Method | Tier |\n|---|---|---:|---|---|---|\n");
    if let Some(power) = &qoi.power {
        write_scalar_row(out, "power", "torque_n_m", &power.torque_n_m)?;
        write_scalar_row(out, "power", "power_w", &power.power_w)?;
        write_scalar_row(out, "power", "np", &power.np)?;
        write_scalar_row(out, "power", "p_over_v_w_m3", &power.p_over_v_w_m3)?;
        write_scalar_row(out, "power", "nq", &power.nq)?;
    }
    if let Some(mixing) = &qoi.mixing {
        write_scalar_row(out, "mixing", "cv0", &mixing.cv0)?;
        write_scalar_row(out, "mixing", "t95_s", &mixing.t95_s)?;
        write_scalar_row(out, "mixing", "t99_s", &mixing.t99_s)?;
    }
    if let Some(shear) = &qoi.shear {
        writeln!(
            out,
            "| shear | gamma_dot_1_s.p95 | {} | 1/s | percentile distribution | Screening |",
            shear.gamma_dot_1_s.p95
        )?;
        writeln!(
            out,
            "| shear | viscous_stress_pa.p95 | {} | Pa | percentile distribution | Screening |",
            shear.viscous_stress_pa.p95
        )?;
    }
    Ok(())
}

fn write_scalar_row(
    out: &mut String,
    section: &str,
    qoi: &str,
    scalar: &QoiScalar,
) -> std::fmt::Result {
    let value = scalar
        .value
        .map(|v| v.to_string())
        .unwrap_or_else(|| "skipped".to_string());
    let units = scalar.provenance.units.as_deref().unwrap_or("unknown");
    let method = scalar.provenance.method.as_deref().unwrap_or("unknown");
    let tier = scalar
        .provenance
        .validation_tier
        .map(|t| format!("{t:?}"))
        .unwrap_or_else(|| "unknown".to_string());
    writeln!(
        out,
        "| {section} | {qoi} | {value} | {units} | {method} | {tier} |"
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fixture_qoi(status: CapabilityStatus) -> QoiBundle {
        let provenance = lbm_core::qoi::QoiProvenance::new(
            vec!["ibm_marker_force".to_string()],
            "last_half_of_run",
            "impeller_marker_set",
            "W",
            "P = omega*Tq",
            lbm_core::qoi::ValidationTier::Screening,
        );
        QoiBundle {
            power: Some(lbm_core::qoi::PowerQoiSection {
                torque_n_m: QoiScalar::measured(1.0, provenance.clone()),
                power_w: QoiScalar::measured(2.0, provenance.clone()),
                rotational_speed_hz: QoiScalar::measured(3.0, provenance.clone()),
                np: QoiScalar::measured(4.0, provenance.clone()),
                p_over_v_w_m3: QoiScalar::measured(5.0, provenance.clone()),
                nq: QoiScalar::skipped("nq", "no discharge surface", provenance),
            }),
            validation_status: vec![lbm_core::qoi::QoiValidationStatus {
                qoi: "power".to_string(),
                status,
                tier: lbm_core::qoi::ValidationTier::Screening,
            }],
            ..QoiBundle::default()
        }
    }

    fn write_manifest(dir: &Path) {
        write_manifest_with_evidence(dir, None);
    }

    fn write_manifest_with_evidence(dir: &Path, evidence_gate: Option<Value>) {
        let mut manifest = serde_json::json!({
            "scenario": "fixture",
            "scenarioHash": "sha256:abc",
            "activeModels": ["single_phase"],
            "unitReport": {"verdict": "screening"},
            "provenance": {}
        });
        if let Some(evidence_gate) = evidence_gate {
            manifest["evidenceGate"] = evidence_gate;
        }
        fs::write(
            dir.join(MANIFEST_PATH),
            serde_json::to_vec_pretty(&manifest).unwrap(),
        )
        .unwrap();
        fs::write(
            dir.join("scenario.json"),
            serde_json::to_vec_pretty(&serde_json::json!({
                "version": "bioprocess-1.0",
                "name": "fixture",
                "credibility_tier": "screening"
            }))
            .unwrap(),
        )
        .unwrap();
    }

    fn ready_evidence_gate_json() -> Value {
        serde_json::json!({
            "qoi_id": "power",
            "vb_id": "VB-01",
            "vb_engineering_green": true,
            "datasets": {
                "calibration": {
                    "cal-power": {
                        "id": "cal-power",
                        "qoi": "power",
                        "source": "fixture",
                        "date": "2026-07-07",
                        "scale": "bench",
                        "operating_condition": "N=1/s",
                        "measurement_uncertainty": 0.1
                    }
                },
                "holdout": {
                    "hold-power": {
                        "id": "hold-power",
                        "qoi": "power",
                        "source": "fixture",
                        "date": "2026-07-07",
                        "scale": "pilot",
                        "operating_condition": "N=2/s",
                        "measurement_uncertainty": 0.1
                    }
                }
            },
            "uq_interval": {
                "q_hat": 1.0,
                "q_lo": 0.9,
                "q_hi": 1.1,
                "method": "fixture"
            },
            "mesh_sensitivity": {
                "additional_cases": 2,
                "variation_fraction": 0.04,
                "summary": "fixture mesh sensitivity"
            },
            "time_step_sensitivity": {
                "additional_cases": 2,
                "variation_fraction": 0.04,
                "summary": "fixture time-step sensitivity"
            },
            "limitation_report": {
                "manifest_path": "limitations.md"
            }
        })
    }

    #[test]
    fn report_generated_from_fixture_qoi_json() {
        let dir = std::env::temp_dir().join(format!("lbm_report_fixture_{}", std::process::id()));
        let _ = fs::remove_dir_all(&dir);
        fs::create_dir_all(&dir).unwrap();
        write_manifest(&dir);
        fs::write(
            dir.join(lbm_scenario::QOI_BUNDLE_JSON),
            serde_json::to_vec_pretty(&fixture_qoi(CapabilityStatus::Experimental)).unwrap(),
        )
        .unwrap();
        let report = generate_report(&dir).unwrap();
        let text = fs::read_to_string(report).unwrap();
        assert!(text.contains("NOT EVIDENCE-GRADE"));
        assert!(text.contains("| power | power_w | 2"));
    }

    #[test]
    fn missing_qoi_json_returns_structured_error() {
        let dir = std::env::temp_dir().join(format!("lbm_report_missing_{}", std::process::id()));
        let _ = fs::remove_dir_all(&dir);
        fs::create_dir_all(&dir).unwrap();
        write_manifest(&dir);
        let err = generate_report(&dir).unwrap_err().to_string();
        assert!(err.contains("missing_qoi_json"));
    }

    #[test]
    fn evidence_ready_qoi_omits_not_evidence_banner() {
        let dir = std::env::temp_dir().join(format!("lbm_report_evidence_{}", std::process::id()));
        let _ = fs::remove_dir_all(&dir);
        fs::create_dir_all(&dir).unwrap();
        write_manifest_with_evidence(&dir, Some(ready_evidence_gate_json()));
        fs::write(
            dir.join(lbm_scenario::QOI_BUNDLE_JSON),
            serde_json::to_vec_pretty(&fixture_qoi(CapabilityStatus::EvidenceReady)).unwrap(),
        )
        .unwrap();
        let report = generate_report(&dir).unwrap();
        let text = fs::read_to_string(report).unwrap();
        assert!(!text.contains("NOT EVIDENCE-GRADE"));
    }

    #[test]
    fn evidence_ready_label_without_artefacts_still_gets_banner() {
        let dir = std::env::temp_dir().join(format!(
            "lbm_report_evidence_blocked_{}",
            std::process::id()
        ));
        let _ = fs::remove_dir_all(&dir);
        fs::create_dir_all(&dir).unwrap();
        write_manifest(&dir);
        fs::write(
            dir.join(lbm_scenario::QOI_BUNDLE_JSON),
            serde_json::to_vec_pretty(&fixture_qoi(CapabilityStatus::EvidenceReady)).unwrap(),
        )
        .unwrap();
        let report = generate_report(&dir).unwrap();
        let text = fs::read_to_string(report).unwrap();
        assert!(text.contains("NOT EVIDENCE-GRADE"));
        assert!(text.contains("NoCalibrationDataset"));
    }
}
