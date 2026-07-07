//! Bioprocess sweep scenario schema and deterministic grid expansion.

use crate::BioprocessScenario;
use serde::{Deserialize, Serialize};

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case", deny_unknown_fields)]
pub struct SweepScenario {
    pub base: BioprocessScenario,
    #[serde(default)]
    pub method: SweepMethod,
    pub sweep_params: Vec<SweepParam>,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum SweepMethod {
    #[default]
    Grid,
    LatinHypercubePlaceholder,
}

impl SweepMethod {
    pub fn capability_status(self) -> &'static str {
        match self {
            Self::Grid => "engineering",
            Self::LatinHypercubePlaceholder => "experimental",
        }
    }
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case", deny_unknown_fields)]
pub struct SweepParam {
    pub parameter_path: String,
    pub values: Vec<f64>,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub struct SweepCase {
    pub parameters: Vec<SweepCaseParameter>,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub struct SweepCaseParameter {
    pub parameter_path: String,
    pub value: f64,
}

pub fn expand_grid(params: &[SweepParam]) -> Vec<SweepCase> {
    let mut cases = vec![SweepCase {
        parameters: Vec::new(),
    }];
    for param in params {
        let mut next = Vec::new();
        for case in &cases {
            for &value in &param.values {
                let mut parameters = case.parameters.clone();
                parameters.push(SweepCaseParameter {
                    parameter_path: param.parameter_path.clone(),
                    value,
                });
                next.push(SweepCase { parameters });
            }
        }
        cases = next;
    }
    cases
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn grid_sweep_expands_expected_cases() {
        let cases = expand_grid(&[
            SweepParam {
                parameter_path: "operation.duration_s".to_string(),
                values: vec![1.0, 2.0],
            },
            SweepParam {
                parameter_path: "run.dt_s".to_string(),
                values: vec![0.1, 0.2, 0.3],
            },
        ]);
        assert_eq!(cases.len(), 6);
        assert_eq!(cases[0].parameters[0].value, 1.0);
        assert_eq!(cases[5].parameters[1].value, 0.3);
    }
}
