//! Performance-regression baseline comparison support for benchmark examples.

use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};

pub const DEFAULT_REGRESSION_THRESHOLD: f64 = 0.10;

#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct BenchBaseline {
    pub schema: u32,
    pub host_tag: String,
    pub generated_at: String,
    pub git_commit: String,
    pub note: String,
    pub cases: Vec<BenchCaseBaseline>,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct BenchCaseBaseline {
    pub case: String,
    pub mlups: f64,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub regression_threshold: Option<f64>,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct BenchMeasurement {
    pub case: String,
    pub mlups: f64,
}

#[derive(Clone, Debug, PartialEq, Serialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum ComparisonFailure {
    Regression {
        case: String,
        baseline_mlups: f64,
        measured_mlups: f64,
        threshold: f64,
    },
    MissingMeasurement {
        case: String,
    },
    MissingBaseline {
        case: String,
    },
}

pub fn compare_measurements(
    baseline: &[BenchCaseBaseline],
    measured: &[BenchMeasurement],
    default_threshold: f64,
) -> Vec<ComparisonFailure> {
    let baseline_by_case: BTreeMap<&str, &BenchCaseBaseline> = baseline
        .iter()
        .map(|case| (case.case.as_str(), case))
        .collect();
    let measured_by_case: BTreeMap<&str, &BenchMeasurement> = measured
        .iter()
        .map(|case| (case.case.as_str(), case))
        .collect();

    let mut failures = Vec::new();
    let cases: BTreeSet<&str> = baseline_by_case
        .keys()
        .chain(measured_by_case.keys())
        .copied()
        .collect();

    for case in cases {
        match (baseline_by_case.get(case), measured_by_case.get(case)) {
            (Some(baseline), Some(measured)) => {
                let threshold = baseline
                    .regression_threshold
                    .unwrap_or(default_threshold)
                    .max(0.0);
                let minimum = baseline.mlups * (1.0 - threshold);
                if !measured.mlups.is_finite() || measured.mlups < minimum {
                    failures.push(ComparisonFailure::Regression {
                        case: case.to_string(),
                        baseline_mlups: baseline.mlups,
                        measured_mlups: measured.mlups,
                        threshold,
                    });
                }
            }
            (Some(_), None) => failures.push(ComparisonFailure::MissingMeasurement {
                case: case.to_string(),
            }),
            (None, Some(_)) => failures.push(ComparisonFailure::MissingBaseline {
                case: case.to_string(),
            }),
            (None, None) => unreachable!("case set is built from both maps"),
        }
    }

    failures
}

#[cfg(test)]
mod tests {
    use super::*;

    fn baseline(case: &str, mlups: f64, threshold: Option<f64>) -> BenchCaseBaseline {
        BenchCaseBaseline {
            case: case.to_string(),
            mlups,
            regression_threshold: threshold,
        }
    }

    fn measured(case: &str, mlups: f64) -> BenchMeasurement {
        BenchMeasurement {
            case: case.to_string(),
            mlups,
        }
    }

    #[test]
    fn comparison_detects_regression() {
        let failures = compare_measurements(
            &[baseline("d2q9-simd-f32", 100.0, Some(0.05))],
            &[measured("d2q9-simd-f32", 94.9)],
            DEFAULT_REGRESSION_THRESHOLD,
        );

        assert_eq!(
            failures,
            vec![ComparisonFailure::Regression {
                case: "d2q9-simd-f32".to_string(),
                baseline_mlups: 100.0,
                measured_mlups: 94.9,
                threshold: 0.05,
            }]
        );
    }

    #[test]
    fn comparison_accepts_within_threshold() {
        let failures = compare_measurements(
            &[baseline("d2q9-simd-f32", 100.0, None)],
            &[measured("d2q9-simd-f32", 90.0)],
            DEFAULT_REGRESSION_THRESHOLD,
        );

        assert!(failures.is_empty());
    }

    #[test]
    fn comparison_treats_missing_cases_as_failures() {
        let failures = compare_measurements(
            &[
                baseline("d2q9-simd-f32", 100.0, None),
                baseline("d3q19-simd-f32", 50.0, None),
            ],
            &[
                measured("d2q9-simd-f32", 100.0),
                measured("extra-case", 10.0),
            ],
            DEFAULT_REGRESSION_THRESHOLD,
        );

        assert_eq!(
            failures,
            vec![
                ComparisonFailure::MissingMeasurement {
                    case: "d3q19-simd-f32".to_string(),
                },
                ComparisonFailure::MissingBaseline {
                    case: "extra-case".to_string(),
                },
            ]
        );
    }
}
