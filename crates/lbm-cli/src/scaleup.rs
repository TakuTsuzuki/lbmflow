use anyhow::{Context, Result};
use lbm_core::scaleup::{
    evaluate_operating_window, ConstraintSet, OperatingPoint, ScaleUpEvaluation, ScaleUpMode,
    ScaleUpQois,
};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::fs;
use std::path::{Path, PathBuf};

pub const SCALEUP_WINDOW_JSON: &str = "scaleup_window.json";

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub struct ScaleUpRequest {
    pub sweep_summary: PathBuf,
    pub constraints: ConstraintSet,
    pub mode: ScaleUpMode,
}

pub fn run(path: &Path, out_dir: Option<PathBuf>) -> Result<PathBuf> {
    let request: ScaleUpRequest = serde_json::from_slice(
        &fs::read(path).with_context(|| format!("cannot read: {}", path.display()))?,
    )?;
    let summary: crate::sweep::SweepSummary = serde_json::from_slice(
        &fs::read(&request.sweep_summary)
            .with_context(|| format!("cannot read: {}", request.sweep_summary.display()))?,
    )?;
    let points = operating_points_from_summary(&summary);
    let evaluation = evaluate_operating_window(&points, &request.constraints, request.mode);
    let root = out_dir.unwrap_or_else(|| {
        request
            .sweep_summary
            .parent()
            .map(Path::to_path_buf)
            .unwrap_or_else(|| PathBuf::from("."))
    });
    fs::create_dir_all(&root)?;
    let path = root.join(SCALEUP_WINDOW_JSON);
    fs::write(&path, serde_json::to_vec_pretty(&evaluation)?)?;
    Ok(path)
}

pub fn operating_points_from_summary(summary: &crate::sweep::SweepSummary) -> Vec<OperatingPoint> {
    summary
        .cases
        .iter()
        .filter_map(|case| {
            let qoi = case.qoi.as_ref()?;
            let mut parameters = BTreeMap::new();
            for parameter in &case.parameters {
                parameters.insert(parameter.parameter_path.clone(), parameter.value);
            }
            Some(OperatingPoint {
                id: case.case_id.clone(),
                parameters,
                qois: qois_from_bundle(qoi),
            })
        })
        .collect()
}

fn qois_from_bundle(qoi: &lbm_core::qoi::QoiBundle) -> ScaleUpQois {
    ScaleUpQois {
        kla_1_s: qoi
            .kla
            .as_ref()
            .and_then(|kla| kla.dynamic_gassing.as_ref().or(kla.pbm.as_ref()))
            .and_then(|q| q.value),
        p_over_v_w_m3: qoi
            .power
            .as_ref()
            .and_then(|power| power.p_over_v_w_m3.value),
        p95_shear_pa: qoi.shear.as_ref().map(|shear| shear.viscous_stress_pa.p95),
        mixing_time_s: qoi.mixing.as_ref().and_then(|mixing| mixing.t95_s.value),
        gas_holdup: qoi
            .gas
            .as_ref()
            .and_then(|gas| gas.gas_holdup.as_ref())
            .and_then(|q| q.value),
    }
}

#[allow(dead_code)]
pub fn evaluate_summary(
    summary: &crate::sweep::SweepSummary,
    constraints: &ConstraintSet,
    mode: ScaleUpMode,
) -> ScaleUpEvaluation {
    evaluate_operating_window(&operating_points_from_summary(summary), constraints, mode)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn provenance(units: &str) -> lbm_core::qoi::QoiProvenance {
        lbm_core::qoi::QoiProvenance::new(
            vec!["field".to_string()],
            "window",
            "region",
            units,
            "synthetic",
            lbm_core::qoi::ValidationTier::Screening,
        )
    }

    fn scalar(value: f64, units: &str) -> lbm_core::qoi::QoiScalar {
        lbm_core::qoi::QoiScalar::measured(value, provenance(units))
    }

    fn bundle() -> lbm_core::qoi::QoiBundle {
        lbm_core::qoi::QoiBundle {
            power: Some(lbm_core::qoi::PowerQoiSection {
                torque_n_m: scalar(1.0, "N*m"),
                power_w: scalar(2.0, "W"),
                rotational_speed_hz: scalar(1.0, "1/s"),
                np: scalar(1.0, "dimensionless"),
                p_over_v_w_m3: scalar(80.0, "W/m^3"),
                nq: lbm_core::qoi::QoiScalar::skipped("nq", "none", provenance("dimensionless")),
            }),
            shear: Some(lbm_core::qoi::ShearQoiSection {
                gamma_dot_1_s: lbm_core::qoi::QoiPercentiles {
                    p50: 1.0,
                    p90: 1.0,
                    p95: 1.0,
                    p99: 1.0,
                    max: 1.0,
                    fraction_above_threshold: 0.0,
                    provenance: provenance("1/s"),
                },
                viscous_stress_pa: lbm_core::qoi::QoiPercentiles {
                    p50: 1.0,
                    p90: 1.5,
                    p95: 1.8,
                    p99: 1.9,
                    max: 2.0,
                    fraction_above_threshold: 0.0,
                    provenance: provenance("Pa"),
                },
                exposure_pa_s: None,
            }),
            validation_status: Vec::new(),
            ..lbm_core::qoi::QoiBundle::default()
        }
    }

    #[test]
    fn summary_cases_convert_to_operating_points() {
        let summary = crate::sweep::SweepSummary {
            method: lbm_scenario::SweepMethod::Grid,
            method_status: "engineering".to_string(),
            cases: vec![crate::sweep::SweepCaseResult {
                case_id: "case_0000".to_string(),
                parameters: vec![lbm_scenario::SweepCaseParameter {
                    parameter_path: "operation.duration_s".to_string(),
                    value: 1.0,
                }],
                qoi: Some(bundle()),
                error: None,
            }],
            aggregate: crate::sweep::AggregatedQoiSummary::default(),
        };
        let points = operating_points_from_summary(&summary);
        assert_eq!(points.len(), 1);
        assert_eq!(points[0].qois.p_over_v_w_m3, Some(80.0));
        assert_eq!(points[0].qois.p95_shear_pa, Some(1.8));
    }
}
