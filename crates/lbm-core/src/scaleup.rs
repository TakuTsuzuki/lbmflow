//! Scale-up operating-window feasibility evaluation from QOI sweep points.

use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, HashMap};

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub struct ConstraintSet {
    pub kla_min_1_s: Option<f64>,
    pub p_over_v_max_w_m3: Option<f64>,
    pub p95_shear_max_pa: Option<f64>,
    pub mixing_time_max_s: Option<f64>,
    pub gas_holdup_range: Option<[f64; 2]>,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum ScaleUpMode {
    ConstantPoverV,
    ConstantTipSpeed,
    ConstantKla,
    ConstantMixingTime,
    CustomWeighted { weights: HashMap<String, f64> },
}

#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub struct ScaleUpQois {
    pub kla_1_s: Option<f64>,
    pub p_over_v_w_m3: Option<f64>,
    pub p95_shear_pa: Option<f64>,
    pub mixing_time_s: Option<f64>,
    pub gas_holdup: Option<f64>,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub struct OperatingPoint {
    pub id: String,
    pub parameters: BTreeMap<String, f64>,
    pub qois: ScaleUpQois,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub struct ScaleUpEvaluation {
    pub mode: ScaleUpMode,
    pub feasible_operating_window: Vec<OperatingPoint>,
    pub conflict_table: Vec<ConstraintConflict>,
    pub constraint_ranking: Vec<ConstraintTightness>,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub struct ConstraintConflict {
    pub constraint: String,
    pub tightest_violation: f64,
    pub points_near_boundary: Vec<String>,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub struct ConstraintTightness {
    pub constraint: String,
    pub tightness: f64,
}

pub fn evaluate_operating_window(
    points: &[OperatingPoint],
    constraints: &ConstraintSet,
    mode: ScaleUpMode,
) -> ScaleUpEvaluation {
    let feasible_operating_window: Vec<_> = points
        .iter()
        .filter(|point| point_feasible(point, constraints))
        .cloned()
        .collect();
    let constraint_ranking = rank_constraints(points, constraints, &mode);
    let conflict_table = if feasible_operating_window.is_empty() {
        conflict_table(points, constraints, &mode)
    } else {
        Vec::new()
    };
    ScaleUpEvaluation {
        mode,
        feasible_operating_window,
        conflict_table,
        constraint_ranking,
    }
}

fn point_feasible(point: &OperatingPoint, constraints: &ConstraintSet) -> bool {
    constraint_violations(point, constraints)
        .into_iter()
        .all(|(_, violation)| violation <= 0.0)
}

fn constraint_violations(
    point: &OperatingPoint,
    constraints: &ConstraintSet,
) -> Vec<(String, f64)> {
    let mut out = Vec::new();
    if let Some(limit) = constraints.kla_min_1_s {
        out.push((
            "kla_min_1_s".to_string(),
            min_constraint_violation(point.qois.kla_1_s, limit),
        ));
    }
    if let Some(limit) = constraints.p_over_v_max_w_m3 {
        out.push((
            "p_over_v_max_w_m3".to_string(),
            max_constraint_violation(point.qois.p_over_v_w_m3, limit),
        ));
    }
    if let Some(limit) = constraints.p95_shear_max_pa {
        out.push((
            "p95_shear_max_pa".to_string(),
            max_constraint_violation(point.qois.p95_shear_pa, limit),
        ));
    }
    if let Some(limit) = constraints.mixing_time_max_s {
        out.push((
            "mixing_time_max_s".to_string(),
            max_constraint_violation(point.qois.mixing_time_s, limit),
        ));
    }
    if let Some([lo, hi]) = constraints.gas_holdup_range {
        let violation = match point.qois.gas_holdup {
            Some(v) if v < lo => (lo - v) / lo.max(1.0e-12),
            Some(v) if v > hi => (v - hi) / hi.max(1.0e-12),
            Some(_) => 0.0,
            None => f64::INFINITY,
        };
        out.push(("gas_holdup_range".to_string(), violation));
    }
    out
}

fn min_constraint_violation(value: Option<f64>, target: f64) -> f64 {
    match value {
        Some(value) if value >= target => 0.0,
        Some(value) => (target - value) / target.max(1.0e-12),
        None => f64::INFINITY,
    }
}

fn max_constraint_violation(value: Option<f64>, limit: f64) -> f64 {
    match value {
        Some(value) if value <= limit => 0.0,
        Some(value) => (value - limit) / limit.max(1.0e-12),
        None => f64::INFINITY,
    }
}

fn rank_constraints(
    points: &[OperatingPoint],
    constraints: &ConstraintSet,
    mode: &ScaleUpMode,
) -> Vec<ConstraintTightness> {
    let mut maxima: HashMap<String, f64> = HashMap::new();
    for point in points {
        for (constraint, violation) in constraint_violations(point, constraints) {
            let tightness = if violation <= 0.0 { 0.0 } else { violation };
            maxima
                .entry(constraint)
                .and_modify(|v| *v = v.max(tightness))
                .or_insert(tightness);
        }
    }
    let mut ranking: Vec<_> = maxima
        .into_iter()
        .map(|(constraint, tightness)| ConstraintTightness {
            constraint,
            tightness,
        })
        .collect();
    ranking.sort_by(|a, b| {
        constraint_priority(mode, &a.constraint)
            .cmp(&constraint_priority(mode, &b.constraint))
            .then_with(|| b.tightness.total_cmp(&a.tightness))
    });
    ranking
}

fn conflict_table(
    points: &[OperatingPoint],
    constraints: &ConstraintSet,
    mode: &ScaleUpMode,
) -> Vec<ConstraintConflict> {
    let ranking = rank_constraints(points, constraints, mode);
    ranking
        .into_iter()
        .map(|rank| {
            let mut near = Vec::new();
            for point in points {
                let violation = constraint_violations(point, constraints)
                    .into_iter()
                    .find(|(name, _)| name == &rank.constraint)
                    .map(|(_, violation)| violation)
                    .unwrap_or(0.0);
                if violation.is_finite() && violation <= rank.tightness.max(0.05) {
                    near.push(point.id.clone());
                }
            }
            ConstraintConflict {
                constraint: rank.constraint,
                tightest_violation: rank.tightness,
                points_near_boundary: near,
            }
        })
        .collect()
}

fn constraint_priority(mode: &ScaleUpMode, constraint: &str) -> usize {
    if let ScaleUpMode::CustomWeighted { weights } = mode {
        if let Some(weight) = weights.get(constraint) {
            let scaled = (*weight * 1_000_000.0).round();
            if scaled.is_finite() && scaled > 0.0 {
                return usize::MAX.saturating_sub(scaled as usize);
            }
        }
    }
    match constraint {
        "kla_min_1_s" => 0,
        "p_over_v_max_w_m3" => 1,
        "p95_shear_max_pa" => 2,
        "mixing_time_max_s" => 3,
        "gas_holdup_range" => 4,
        _ => 100,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn point(id: &str, kla: f64, pv: f64, shear: f64, mixing: f64, gas: f64) -> OperatingPoint {
        OperatingPoint {
            id: id.to_string(),
            parameters: BTreeMap::new(),
            qois: ScaleUpQois {
                kla_1_s: Some(kla),
                p_over_v_w_m3: Some(pv),
                p95_shear_pa: Some(shear),
                mixing_time_s: Some(mixing),
                gas_holdup: Some(gas),
            },
        }
    }

    fn constraints() -> ConstraintSet {
        ConstraintSet {
            kla_min_1_s: Some(0.01),
            p_over_v_max_w_m3: Some(100.0),
            p95_shear_max_pa: Some(2.0),
            mixing_time_max_s: Some(60.0),
            gas_holdup_range: Some([0.02, 0.10]),
        }
    }

    #[test]
    fn synthetic_feasible_window_recovered() {
        let points = vec![
            point("bad", 0.005, 80.0, 1.0, 50.0, 0.04),
            point("good", 0.02, 90.0, 1.5, 40.0, 0.05),
        ];
        let eval = evaluate_operating_window(&points, &constraints(), ScaleUpMode::ConstantPoverV);
        assert_eq!(eval.feasible_operating_window.len(), 1);
        assert_eq!(eval.feasible_operating_window[0].id, "good");
        assert!(eval.conflict_table.is_empty());
    }

    #[test]
    fn synthetic_no_feasible_set_gives_conflict_table() {
        let points = vec![
            point("low_kla", 0.005, 80.0, 1.0, 50.0, 0.04),
            point("high_shear", 0.02, 90.0, 3.0, 40.0, 0.05),
        ];
        let eval = evaluate_operating_window(&points, &constraints(), ScaleUpMode::ConstantKla);
        assert!(eval.feasible_operating_window.is_empty());
        assert!(!eval.conflict_table.is_empty());
    }

    #[test]
    fn constraint_ranking_correct() {
        let points = vec![point("bad", 0.005, 300.0, 1.0, 50.0, 0.04)];
        let eval =
            evaluate_operating_window(&points, &constraints(), ScaleUpMode::ConstantTipSpeed);
        let constraints: Vec<_> = eval
            .constraint_ranking
            .iter()
            .map(|rank| rank.constraint.as_str())
            .collect();
        assert_eq!(
            constraints,
            vec![
                "kla_min_1_s",
                "p_over_v_max_w_m3",
                "p95_shear_max_pa",
                "mixing_time_max_s",
                "gas_holdup_range"
            ]
        );
    }
}
