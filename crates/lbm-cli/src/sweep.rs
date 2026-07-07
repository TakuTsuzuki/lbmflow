use anyhow::{Context, Result};
use lbm_core::qoi::QoiBundle;
use lbm_scenario::{expand_grid, SweepCase, SweepCaseParameter, SweepMethod, SweepScenario};
use serde::{Deserialize, Serialize};
use serde_json::Value;
use std::collections::BTreeMap;
use std::fs;
use std::path::{Path, PathBuf};

pub const SWEEP_SUMMARY_JSON: &str = "sweep_summary.json";

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub struct SweepSummary {
    pub method: SweepMethod,
    pub method_status: String,
    pub cases: Vec<SweepCaseResult>,
    pub aggregate: AggregatedQoiSummary,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub struct SweepCaseResult {
    pub case_id: String,
    pub parameters: Vec<SweepCaseParameter>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub qoi: Option<QoiBundle>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub error: Option<String>,
}

#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub struct AggregatedQoiSummary {
    pub qoi: BTreeMap<String, AggregatedScalar>,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub struct AggregatedScalar {
    pub count: usize,
    pub min: f64,
    pub mean: f64,
    pub max: f64,
}

pub fn run(path: &Path, out_dir: Option<PathBuf>) -> Result<PathBuf> {
    let text =
        fs::read_to_string(path).with_context(|| format!("cannot read: {}", path.display()))?;
    let sweep: SweepScenario = serde_json::from_str(&text)?;
    let mut base_value: Value = serde_json::to_value(&sweep.base)?;
    let root =
        out_dir.unwrap_or_else(|| PathBuf::from("out").join(format!("{}_sweep", sweep.base.name)));
    fs::create_dir_all(&root)?;
    let cases = match sweep.method {
        SweepMethod::Grid | SweepMethod::LatinHypercubePlaceholder => {
            expand_grid(&sweep.sweep_params)
        }
    };
    let mut results = Vec::new();
    for (index, case) in cases.iter().enumerate() {
        let case_id = format!("case_{index:04}");
        let case_dir = root.join(&case_id);
        fs::create_dir_all(&case_dir)?;
        let result = run_case(&mut base_value, case, &case_id, &case_dir);
        results.push(result);
    }
    let summary = SweepSummary {
        method: sweep.method,
        method_status: sweep.method.capability_status().to_string(),
        aggregate: aggregate_qoi_bundles(results.iter().filter_map(|r| r.qoi.as_ref())),
        cases: results,
    };
    let path = root.join(SWEEP_SUMMARY_JSON);
    fs::write(&path, serde_json::to_vec_pretty(&summary)?)?;
    Ok(path)
}

fn run_case(
    base_value: &mut Value,
    case: &SweepCase,
    case_id: &str,
    case_dir: &Path,
) -> SweepCaseResult {
    let mut value = base_value.clone();
    for parameter in &case.parameters {
        if let Err(err) = set_json_path(&mut value, &parameter.parameter_path, parameter.value) {
            return failed_case(case_id, case, err);
        }
    }
    match serde_json::from_value::<lbm_scenario::BioprocessScenario>(value) {
        Ok(scenario) => match crate::runner::run_bioprocess_single_phase(&scenario, case_dir) {
            Ok(_) => match read_qoi(case_dir) {
                Ok(qoi) => SweepCaseResult {
                    case_id: case_id.to_string(),
                    parameters: case.parameters.clone(),
                    qoi: Some(qoi),
                    error: None,
                },
                Err(err) => failed_case(case_id, case, err.to_string()),
            },
            Err(err) => failed_case(case_id, case, err.to_string()),
        },
        Err(err) => failed_case(case_id, case, err.to_string()),
    }
}

fn read_qoi(case_dir: &Path) -> Result<QoiBundle> {
    let path = case_dir.join(lbm_scenario::QOI_BUNDLE_JSON);
    Ok(serde_json::from_slice(&fs::read(path)?)?)
}

fn failed_case(case_id: &str, case: &SweepCase, error: impl Into<String>) -> SweepCaseResult {
    SweepCaseResult {
        case_id: case_id.to_string(),
        parameters: case.parameters.clone(),
        qoi: None,
        error: Some(error.into()),
    }
}

fn set_json_path(root: &mut Value, path: &str, value: f64) -> Result<(), String> {
    let trimmed = path.strip_prefix("$.").unwrap_or(path);
    let parts: Vec<_> = trimmed.split('.').filter(|part| !part.is_empty()).collect();
    if parts.is_empty() {
        return Err("parameter_path is empty".to_string());
    }
    let mut cursor = root;
    for part in &parts[..parts.len() - 1] {
        cursor = cursor
            .get_mut(*part)
            .ok_or_else(|| format!("parameter_path component not found: {part}"))?;
    }
    let leaf = parts[parts.len() - 1];
    let slot = cursor
        .get_mut(leaf)
        .ok_or_else(|| format!("parameter_path leaf not found: {leaf}"))?;
    *slot = Value::from(value);
    Ok(())
}

pub fn aggregate_qoi_bundles<'a>(
    bundles: impl Iterator<Item = &'a QoiBundle>,
) -> AggregatedQoiSummary {
    let mut values: BTreeMap<String, Vec<f64>> = BTreeMap::new();
    for bundle in bundles {
        for (key, value) in bundle.scalar_values() {
            values.entry(key).or_default().push(value);
        }
    }
    let qoi = values
        .into_iter()
        .map(|(key, values)| {
            let count = values.len();
            let min = values.iter().copied().fold(f64::INFINITY, f64::min);
            let max = values.iter().copied().fold(f64::NEG_INFINITY, f64::max);
            let mean = values.iter().sum::<f64>() / count as f64;
            (
                key,
                AggregatedScalar {
                    count,
                    min,
                    mean,
                    max,
                },
            )
        })
        .collect();
    AggregatedQoiSummary { qoi }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn provenance() -> lbm_core::qoi::QoiProvenance {
        lbm_core::qoi::QoiProvenance::new(
            vec!["field".to_string()],
            "window",
            "region",
            "W",
            "synthetic",
            lbm_core::qoi::ValidationTier::Screening,
        )
    }

    fn bundle(power: f64) -> QoiBundle {
        QoiBundle {
            power: Some(lbm_core::qoi::PowerQoiSection {
                torque_n_m: lbm_core::qoi::QoiScalar::measured(1.0, provenance()),
                power_w: lbm_core::qoi::QoiScalar::measured(power, provenance()),
                rotational_speed_hz: lbm_core::qoi::QoiScalar::measured(1.0, provenance()),
                np: lbm_core::qoi::QoiScalar::measured(power / 2.0, provenance()),
                p_over_v_w_m3: lbm_core::qoi::QoiScalar::measured(power / 3.0, provenance()),
                nq: lbm_core::qoi::QoiScalar::skipped("nq", "none", provenance()),
            }),
            validation_status: Vec::new(),
            ..QoiBundle::default()
        }
    }

    #[test]
    fn aggregation_roundtrip() {
        let summary = aggregate_qoi_bundles([bundle(2.0), bundle(4.0)].iter());
        let text = serde_json::to_string(&summary).unwrap();
        let parsed: AggregatedQoiSummary = serde_json::from_str(&text).unwrap();
        assert_eq!(parsed.qoi["power.power_w"].count, 2);
        assert_eq!(parsed.qoi["power.power_w"].mean, 3.0);
    }

    #[test]
    fn failed_case_recorded() {
        let case = SweepCase {
            parameters: vec![SweepCaseParameter {
                parameter_path: "missing.path".to_string(),
                value: 1.0,
            }],
        };
        let result = failed_case("case_0000", &case, "bad path");
        assert_eq!(result.case_id, "case_0000");
        assert!(result.error.unwrap().contains("bad path"));
    }

    #[test]
    fn json_path_setter_updates_nested_number() {
        let mut value = serde_json::json!({"run": {"dt_s": 0.1}});
        set_json_path(&mut value, "$.run.dt_s", 0.2).unwrap();
        assert_eq!(value["run"]["dt_s"], 0.2);
    }
}
