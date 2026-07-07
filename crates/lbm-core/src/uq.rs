//! Uncertainty interval and local sensitivity helpers for QOI sweeps.

use crate::qoi::QoiInterval;
use serde::{Deserialize, Serialize};

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum UqComponentKind {
    ModelForm,
    Parameter,
    Numerical,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub struct UqComponent {
    pub kind: UqComponentKind,
    pub q_lo: f64,
    pub q_hi: f64,
    pub method: String,
}

pub fn combine_interval(q_hat: f64, components: &[UqComponent]) -> Option<QoiInterval> {
    if !(q_hat.is_finite()
        && components
            .iter()
            .all(|c| c.q_lo.is_finite() && c.q_hi.is_finite()))
    {
        return None;
    }
    let mut lo = q_hat;
    let mut hi = q_hat;
    let mut methods = Vec::new();
    for component in components {
        lo = lo.min(component.q_lo);
        hi = hi.max(component.q_hi);
        methods.push(format!("{:?}:{}", component.kind, component.method));
    }
    Some(QoiInterval {
        q_hat,
        q_lo: lo,
        q_hi: hi,
        method: methods.join(";"),
    })
}

pub fn one_factor_local_sensitivity(
    base_parameter: f64,
    base_qoi: f64,
    perturbed_parameter: f64,
    perturbed_qoi: f64,
) -> Option<f64> {
    if !(base_parameter.is_finite()
        && base_qoi.is_finite()
        && perturbed_parameter.is_finite()
        && perturbed_qoi.is_finite())
    {
        return None;
    }
    let dp = perturbed_parameter - base_parameter;
    if dp == 0.0 {
        return None;
    }
    Some((perturbed_qoi - base_qoi) / dp)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn synthetic_sensitivity_matches_analytical() {
        let derivative = one_factor_local_sensitivity(2.0, 9.0, 2.5, 11.0).unwrap();
        assert!((derivative - 4.0).abs() < 1.0e-12);
    }

    #[test]
    fn interval_spans_model_parameter_and_numerical_components() {
        let interval = combine_interval(
            10.0,
            &[
                UqComponent {
                    kind: UqComponentKind::ModelForm,
                    q_lo: 8.0,
                    q_hi: 11.0,
                    method: "closure_set".to_string(),
                },
                UqComponent {
                    kind: UqComponentKind::Numerical,
                    q_lo: 9.5,
                    q_hi: 12.0,
                    method: "mesh_dt".to_string(),
                },
            ],
        )
        .unwrap();
        assert_eq!(interval.q_lo, 8.0);
        assert_eq!(interval.q_hi, 12.0);
        assert!(interval.method.contains("ModelForm"));
    }
}
