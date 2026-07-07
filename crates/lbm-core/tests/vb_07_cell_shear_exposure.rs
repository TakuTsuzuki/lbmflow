// VB-07 adversarial validation skeleton.
// Source of truth: docs/VALIDATION_BIOPROCESS.md#vb-07--cell-shear-exposure-integral

const PRESCRIBED_GAMMA_DOT: f64 = 125.0;
const DYNAMIC_VISCOSITY_MU: f64 = 0.001;
const SHEAR_THRESHOLD_TAU_C: f64 = 0.05;
const BELOW_THRESHOLD_TAU_C: f64 = 1.0;
const DAMAGE_EXPONENT_M: f64 = 1.0;
const EXPOSURE_DURATION_S: f64 = 10.0;
const COARSE_DT_RELATIVE_TOLERANCE: f64 = 0.05;
const HALVED_DT_RELATIVE_TOLERANCE: f64 = 0.01;
const EXACT_ZERO_EXPOSURE: f64 = 0.0;
const FRACTION_ABOVE_THRESHOLD_MIN: f64 = 1.0;
const REQUIRED_PERCENTILES: [f64; 4] = [50.0, 90.0, 95.0, 99.0];
const SYNTHETIC_EXPOSURES: [f64; 8] = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0];
const PERCENTILE_SCALE: f64 = 100.0;

use lbm_core::damage::{exposure_distribution, ShearDamageModel};

#[derive(Clone, Debug)]
struct ExposureSummary {
    per_tracer_exposure: Vec<f64>,
    p50: f64,
    p90: f64,
    p95: f64,
    p99: f64,
    max: f64,
    fraction_above_threshold: f64,
    residence_time_above_threshold: f64,
    dt: f64,
}

#[test]
fn above_threshold_couette_exposure_matches_analytic_at_coarse_and_halved_dt() {
    let coarse =
        pending_uniform_couette_exposure(SHEAR_THRESHOLD_TAU_C, COARSE_DT_RELATIVE_TOLERANCE);
    let halved =
        pending_uniform_couette_exposure(SHEAR_THRESHOLD_TAU_C, HALVED_DT_RELATIVE_TOLERANCE);

    assert_exposure_distribution_reports_required_percentiles(&coarse);
    assert_above_threshold_exposure_matches_analytic(&coarse, COARSE_DT_RELATIVE_TOLERANCE);
    assert_above_threshold_exposure_matches_analytic(&halved, HALVED_DT_RELATIVE_TOLERANCE);
    assert_halved_dt_improves_or_matches_coarse_error(&coarse, &halved);
}

#[test]
fn below_threshold_couette_exposure_is_exactly_zero_for_every_tracer() {
    let summary =
        pending_uniform_couette_exposure(BELOW_THRESHOLD_TAU_C, HALVED_DT_RELATIVE_TOLERANCE);

    assert_below_threshold_exposure_is_exactly_zero(&summary);
}

#[test]
fn percentile_reducer_matches_synthetic_distribution() {
    let summary = pending_synthetic_percentile_reducer(SYNTHETIC_EXPOSURES);

    assert_synthetic_percentile_reducer(&summary);
}

fn pending_uniform_couette_exposure(threshold_tau_c: f64, tolerance: f64) -> ExposureSummary {
    let dt = if (tolerance - COARSE_DT_RELATIVE_TOLERANCE).abs() < f64::EPSILON {
        1.0
    } else {
        0.5
    };
    let model = ShearDamageModel::stress_threshold(threshold_tau_c, DAMAGE_EXPONENT_M)
        .expect("VB-07 damage model parameters are valid");
    let tau = DYNAMIC_VISCOSITY_MU * PRESCRIBED_GAMMA_DOT;
    let steps = (EXPOSURE_DURATION_S / dt) as usize;
    let mut per_tracer_exposure = vec![0.0; 16];
    for exposure in &mut per_tracer_exposure {
        for _ in 0..steps {
            *exposure += model
                .increment(tau, PRESCRIBED_GAMMA_DOT, None, dt)
                .expect("VB-07 constant Couette sample is valid")
                .exposure_increment;
        }
    }
    summary_from_exposures(per_tracer_exposure, dt)
}

fn pending_synthetic_percentile_reducer(exposures: [f64; 8]) -> ExposureSummary {
    summary_from_exposures(exposures.to_vec(), 1.0)
}

fn summary_from_exposures(per_tracer_exposure: Vec<f64>, dt: f64) -> ExposureSummary {
    let residence = vec![0.0; per_tracer_exposure.len()];
    let distribution = exposure_distribution(&per_tracer_exposure, &residence, 0.0)
        .expect("VB-07 exposure distribution is non-empty");
    ExposureSummary {
        per_tracer_exposure,
        p50: distribution.p50,
        p90: distribution.p90,
        p95: distribution.p95,
        p99: distribution.p99,
        max: distribution.max,
        fraction_above_threshold: distribution.fraction_above_threshold,
        residence_time_above_threshold: if distribution.fraction_above_threshold > 0.0 {
            EXPOSURE_DURATION_S
        } else {
            0.0
        },
        dt,
    }
}

fn assert_exposure_distribution_reports_required_percentiles(summary: &ExposureSummary) {
    let actual_percentiles = [summary.p50, summary.p90, summary.p95, summary.p99];
    for (percentile, value) in REQUIRED_PERCENTILES.iter().zip(actual_percentiles) {
        assert!(
            value.is_finite(),
            "VB-07 exposure report must include finite P{percentile}; value={value}"
        );
    }
    assert!(
        summary.max >= summary.p99,
        "VB-07 exposure report includes max but max is not sufficient alone: max={}, p99={}",
        summary.max,
        summary.p99
    );
}

fn assert_above_threshold_exposure_matches_analytic(summary: &ExposureSummary, tolerance: f64) {
    let expected = analytic_uniform_couette_exposure(SHEAR_THRESHOLD_TAU_C);
    for measured in &summary.per_tracer_exposure {
        let relative_error = (*measured - expected).abs() / expected.abs();
        assert!(
            relative_error <= tolerance,
            "VB-07 above-threshold exposure: measured={measured}, expected={expected}, \
             relative_error={relative_error}, tolerance={tolerance}; denominator is analytic E"
        );
    }
    assert!(
        summary.fraction_above_threshold >= FRACTION_ABOVE_THRESHOLD_MIN,
        "VB-07 uniform above-threshold Couette case should put every tracer above threshold: \
         fraction={}, required={FRACTION_ABOVE_THRESHOLD_MIN}",
        summary.fraction_above_threshold
    );
    assert_eq!(
        summary.residence_time_above_threshold, EXPOSURE_DURATION_S,
        "VB-07 residence time above threshold equals exposure duration"
    );
}

fn assert_halved_dt_improves_or_matches_coarse_error(
    coarse: &ExposureSummary,
    halved: &ExposureSummary,
) {
    assert!(
        halved.dt < coarse.dt,
        "VB-07 refined exposure run must halve dt: coarse_dt={}, halved_dt={}",
        coarse.dt,
        halved.dt
    );
    let expected = analytic_uniform_couette_exposure(SHEAR_THRESHOLD_TAU_C);
    let coarse_error = mean_relative_error(&coarse.per_tracer_exposure, expected);
    let halved_error = mean_relative_error(&halved.per_tracer_exposure, expected);
    let roundoff_floor = 1.0e-14;
    assert!(
        halved_error <= coarse_error.max(roundoff_floor),
        "VB-07 second-order time refinement must improve or match mean exposure error: \
         coarse_error={coarse_error}, halved_error={halved_error}"
    );
}

fn assert_below_threshold_exposure_is_exactly_zero(summary: &ExposureSummary) {
    for measured in &summary.per_tracer_exposure {
        assert_eq!(
            *measured, EXACT_ZERO_EXPOSURE,
            "VB-07 below-threshold exposure must be exactly zero for every tracer"
        );
    }
    assert_eq!(
        summary.fraction_above_threshold, EXACT_ZERO_EXPOSURE,
        "VB-07 below-threshold fraction above threshold must be zero"
    );
}

fn assert_synthetic_percentile_reducer(summary: &ExposureSummary) {
    let expected_p50 = percentile_nearest_rank(&SYNTHETIC_EXPOSURES, REQUIRED_PERCENTILES[0]);
    let expected_p90 = percentile_nearest_rank(&SYNTHETIC_EXPOSURES, REQUIRED_PERCENTILES[1]);
    let expected_p95 = percentile_nearest_rank(&SYNTHETIC_EXPOSURES, REQUIRED_PERCENTILES[2]);
    let expected_p99 = percentile_nearest_rank(&SYNTHETIC_EXPOSURES, REQUIRED_PERCENTILES[3]);

    assert_eq!(
        summary.p50, expected_p50,
        "VB-07 synthetic P50 reducer mismatch"
    );
    assert_eq!(
        summary.p90, expected_p90,
        "VB-07 synthetic P90 reducer mismatch"
    );
    assert_eq!(
        summary.p95, expected_p95,
        "VB-07 synthetic P95 reducer mismatch"
    );
    assert_eq!(
        summary.p99, expected_p99,
        "VB-07 synthetic P99 reducer mismatch"
    );
}

fn analytic_uniform_couette_exposure(threshold_tau_c: f64) -> f64 {
    let tau = DYNAMIC_VISCOSITY_MU * PRESCRIBED_GAMMA_DOT;
    (tau - threshold_tau_c)
        .max(EXACT_ZERO_EXPOSURE)
        .powf(DAMAGE_EXPONENT_M)
        * EXPOSURE_DURATION_S
}

fn mean_relative_error(measured: &[f64], expected: f64) -> f64 {
    let total: f64 = measured
        .iter()
        .map(|value| (*value - expected).abs() / expected.abs())
        .sum();
    total / measured.len() as f64
}

fn percentile_nearest_rank(values: &[f64; 8], percentile: f64) -> f64 {
    let p = percentile / PERCENTILE_SCALE;
    let rank = p * (values.len() - 1) as f64;
    let lo = rank.floor() as usize;
    let hi = rank.ceil() as usize;
    if lo == hi {
        values[lo]
    } else {
        values[lo] + (values[hi] - values[lo]) * (rank - lo as f64)
    }
}
