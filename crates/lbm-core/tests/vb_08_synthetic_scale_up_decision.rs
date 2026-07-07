// VB-08 adversarial validation skeleton.
// Source of truth: docs/VALIDATION_BIOPROCESS.md#vb-08--synthetic-scale-up-decision

const FEASIBLE_SET_RELATIVE_TOLERANCE: f64 = 0.02;
const KLA_TARGET_1_PER_S: f64 = 0.025;
const P_OVER_V_LIMIT: f64 = 400.0;
const P95_SHEAR_LIMIT: f64 = 0.20;
const MIXING_TIME_LIMIT_S: f64 = 120.0;
const GAS_FLOW_QG: f64 = 1.0e-5;
const N_MIN_HZ: f64 = 0.5;
const N_MAX_HZ: f64 = 3.0;
const N_STEP_HZ: f64 = 0.1;
const LARGE_TANK_NP_INTERCEPT: f64 = 1.5;
const LARGE_TANK_NP_SLOPE: f64 = 0.1;
const LARGE_TANK_KLA_INTERCEPT: f64 = 0.010;
const LARGE_TANK_KLA_N_SLOPE: f64 = 0.009;
const LARGE_TANK_KLA_QG_SLOPE: f64 = 100.0;
const LARGE_TANK_P_OVER_V_COEFFICIENT: f64 = 45.0;
const LARGE_TANK_P95_SHEAR_INTERCEPT: f64 = 0.040;
const LARGE_TANK_P95_SHEAR_N_SLOPE: f64 = 0.045;
const LARGE_TANK_MIXING_TIME_COEFFICIENT: f64 = 180.0;

use lbm_core::scaleup::{
    evaluate_operating_window, ConstraintSet, OperatingPoint as CoreOperatingPoint, ScaleUpMode,
    ScaleUpQois,
};
use std::collections::BTreeMap;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum ConstraintId {
    Kla,
    PowerPerVolume,
    MixingTime,
    P95Shear,
    GasHoldup,
}

#[derive(Clone, Copy, Debug)]
struct OperatingPoint {
    n_hz: f64,
    qg: f64,
    np: f64,
    kla_1_per_s: f64,
    p_over_v: f64,
    p95_shear: f64,
    mixing_time_s: f64,
}

#[derive(Clone, Debug)]
struct ScaleUpDecision {
    feasible_large_tank_points: Vec<OperatingPoint>,
    conflict_table: Vec<ConstraintId>,
    tightest_constraint: Option<ConstraintId>,
    ranked_constraints: Vec<ConstraintId>,
}

#[test]
fn evaluator_recovers_analytic_large_tank_feasible_set_from_synthetic_maps() {
    let decision = pending_scale_up_decision_for_feasible_maps();

    assert_feasible_set_matches_analytic_window(&decision);
}

#[test]
fn infeasible_case_emits_explicit_conflict_table() {
    let decision = pending_scale_up_decision_for_infeasible_maps();

    assert_infeasible_case_has_explicit_conflict_table(&decision);
}

#[test]
fn tightest_constraint_and_default_priority_are_reported() {
    let decision = pending_scale_up_decision_for_infeasible_maps();

    assert_tightest_constraint_correctly_identified(&decision);
    assert_default_priority_order_is_documented_order(&decision);
}

fn pending_scale_up_decision_for_feasible_maps() -> ScaleUpDecision {
    let points: Vec<_> = analytic_large_tank_points()
        .into_iter()
        .map(core_operating_point)
        .collect();
    decision_from_core(evaluate_operating_window(
        &points,
        &core_constraints(P95_SHEAR_LIMIT),
        ScaleUpMode::ConstantKla,
    ))
}

fn pending_scale_up_decision_for_infeasible_maps() -> ScaleUpDecision {
    let points: Vec<_> = [1.0, 1.5, 2.0]
        .into_iter()
        .map(|n_hz| {
            let mut point = analytic_large_tank_map(n_hz, GAS_FLOW_QG);
            point.kla_1_per_s = KLA_TARGET_1_PER_S * 1.2;
            point.p_over_v = P_OVER_V_LIMIT * 0.5;
            point.mixing_time_s = MIXING_TIME_LIMIT_S * 0.5;
            point.p95_shear = P95_SHEAR_LIMIT * 1.5;
            core_operating_point(point)
        })
        .collect();
    decision_from_core(evaluate_operating_window(
        &points,
        &core_constraints(P95_SHEAR_LIMIT),
        ScaleUpMode::ConstantKla,
    ))
}

fn core_constraints(p95_shear_limit: f64) -> ConstraintSet {
    ConstraintSet {
        kla_min_1_s: Some(KLA_TARGET_1_PER_S),
        p_over_v_max_w_m3: Some(P_OVER_V_LIMIT),
        p95_shear_max_pa: Some(p95_shear_limit),
        mixing_time_max_s: Some(MIXING_TIME_LIMIT_S),
        gas_holdup_range: None,
    }
}

fn core_operating_point(point: OperatingPoint) -> CoreOperatingPoint {
    let mut parameters = BTreeMap::new();
    parameters.insert("N_hz".to_string(), point.n_hz);
    parameters.insert("Qg".to_string(), point.qg);
    CoreOperatingPoint {
        id: format!("N={:.1}", point.n_hz),
        parameters,
        qois: ScaleUpQois {
            kla_1_s: Some(point.kla_1_per_s),
            p_over_v_w_m3: Some(point.p_over_v),
            p95_shear_pa: Some(point.p95_shear),
            mixing_time_s: Some(point.mixing_time_s),
            gas_holdup: None,
        },
    }
}

fn decision_from_core(evaluation: lbm_core::scaleup::ScaleUpEvaluation) -> ScaleUpDecision {
    let feasible_large_tank_points = evaluation
        .feasible_operating_window
        .iter()
        .map(local_operating_point)
        .collect();
    let conflict_table = evaluation
        .conflict_table
        .iter()
        .filter_map(|conflict| constraint_id(&conflict.constraint))
        .collect();
    let tightest_constraint = evaluation
        .conflict_table
        .iter()
        .max_by(|a, b| a.tightest_violation.total_cmp(&b.tightest_violation))
        .and_then(|conflict| constraint_id(&conflict.constraint));
    let ranked_constraints = evaluation
        .constraint_ranking
        .iter()
        .filter_map(|rank| constraint_id(&rank.constraint))
        .collect();
    ScaleUpDecision {
        feasible_large_tank_points,
        conflict_table,
        tightest_constraint,
        ranked_constraints,
    }
}

fn local_operating_point(point: &CoreOperatingPoint) -> OperatingPoint {
    OperatingPoint {
        n_hz: point.parameters["N_hz"],
        qg: point.parameters["Qg"],
        np: np_from_n(point.parameters["N_hz"]),
        kla_1_per_s: point.qois.kla_1_s.unwrap(),
        p_over_v: point.qois.p_over_v_w_m3.unwrap(),
        p95_shear: point.qois.p95_shear_pa.unwrap(),
        mixing_time_s: point.qois.mixing_time_s.unwrap(),
    }
}

fn constraint_id(name: &str) -> Option<ConstraintId> {
    match name {
        "kla_min_1_s" => Some(ConstraintId::Kla),
        "p_over_v_max_w_m3" => Some(ConstraintId::PowerPerVolume),
        "p95_shear_max_pa" => Some(ConstraintId::P95Shear),
        "mixing_time_max_s" => Some(ConstraintId::MixingTime),
        "gas_holdup_range" => Some(ConstraintId::GasHoldup),
        _ => None,
    }
}

fn assert_feasible_set_matches_analytic_window(decision: &ScaleUpDecision) {
    let expected = analytic_feasible_large_tank_points();
    assert!(
        !decision.feasible_large_tank_points.is_empty(),
        "VB-08 feasible synthetic map must produce a non-empty large-tank operating window"
    );
    assert_eq!(
        decision.feasible_large_tank_points.len(),
        expected.len(),
        "VB-08 feasible-set cardinality mismatch: measured={}, expected={}",
        decision.feasible_large_tank_points.len(),
        expected.len()
    );
    for (measured, expected) in decision.feasible_large_tank_points.iter().zip(expected) {
        assert_operating_point_within_tolerance(measured, &expected);
    }
}

fn assert_infeasible_case_has_explicit_conflict_table(decision: &ScaleUpDecision) {
    assert!(
        decision.feasible_large_tank_points.is_empty(),
        "VB-08 infeasible synthetic map must return an empty feasible set"
    );
    assert!(
        !decision.conflict_table.is_empty(),
        "VB-08 infeasible scale-up case must emit an explicit conflict table"
    );
}

fn assert_tightest_constraint_correctly_identified(decision: &ScaleUpDecision) {
    assert_eq!(
        decision.tightest_constraint,
        Some(ConstraintId::P95Shear),
        "VB-08 tightest constraint should be P95 shear for the synthetic infeasible map"
    );
}

fn assert_default_priority_order_is_documented_order(decision: &ScaleUpDecision) {
    let expected_order = [
        ConstraintId::Kla,
        ConstraintId::PowerPerVolume,
        ConstraintId::P95Shear,
        ConstraintId::MixingTime,
    ];
    assert!(
        decision.ranked_constraints.starts_with(&expected_order),
        "VB-08 default constraint priority must be constant kLa -> P/V -> tip speed -> mixing time; \
         measured={:?}",
        decision.ranked_constraints
    );
}

fn assert_operating_point_within_tolerance(measured: &OperatingPoint, expected: &OperatingPoint) {
    assert_relative_agreement(measured.n_hz, expected.n_hz, "N");
    assert_relative_agreement(measured.qg, expected.qg, "Qg");
    assert_relative_agreement(measured.np, expected.np, "Np");
    assert_relative_agreement(measured.kla_1_per_s, expected.kla_1_per_s, "kLa");
    assert_relative_agreement(measured.p_over_v, expected.p_over_v, "P/V");
    assert_relative_agreement(measured.p95_shear, expected.p95_shear, "P95 shear");
    assert_relative_agreement(
        measured.mixing_time_s,
        expected.mixing_time_s,
        "mixing time",
    );
}

fn assert_relative_agreement(measured: f64, expected: f64, label: &str) {
    let relative_error = (measured - expected).abs() / expected.abs();
    assert!(
        relative_error <= FEASIBLE_SET_RELATIVE_TOLERANCE,
        "VB-08 feasible operating point {label}: measured={measured}, expected={expected}, \
         relative_error={relative_error}, tolerance={FEASIBLE_SET_RELATIVE_TOLERANCE}; \
         denominator is analytic synthetic-map value"
    );
}

fn analytic_feasible_large_tank_points() -> Vec<OperatingPoint> {
    analytic_large_tank_points()
        .into_iter()
        .filter(|point| {
            point.kla_1_per_s >= KLA_TARGET_1_PER_S
                && point.p_over_v <= P_OVER_V_LIMIT
                && point.p95_shear <= P95_SHEAR_LIMIT
                && point.mixing_time_s <= MIXING_TIME_LIMIT_S
        })
        .collect()
}

fn analytic_large_tank_points() -> Vec<OperatingPoint> {
    let mut points = Vec::new();
    let mut n_hz = N_MIN_HZ;
    while n_hz <= N_MAX_HZ {
        points.push(analytic_large_tank_map(n_hz, GAS_FLOW_QG));
        n_hz += N_STEP_HZ;
    }
    points
}

fn np_from_n(n_hz: f64) -> f64 {
    LARGE_TANK_NP_INTERCEPT + LARGE_TANK_NP_SLOPE * n_hz
}

fn analytic_large_tank_map(n_hz: f64, qg: f64) -> OperatingPoint {
    let np = np_from_n(n_hz);
    let kla_1_per_s =
        LARGE_TANK_KLA_INTERCEPT + LARGE_TANK_KLA_N_SLOPE * n_hz + LARGE_TANK_KLA_QG_SLOPE * qg;
    let p_over_v = LARGE_TANK_P_OVER_V_COEFFICIENT * np * n_hz.powi(3);
    let p95_shear = LARGE_TANK_P95_SHEAR_INTERCEPT + LARGE_TANK_P95_SHEAR_N_SLOPE * n_hz;
    let mixing_time_s = LARGE_TANK_MIXING_TIME_COEFFICIENT / n_hz;
    OperatingPoint {
        n_hz,
        qg,
        np,
        kla_1_per_s,
        p_over_v,
        p95_shear,
        mixing_time_s,
    }
}
