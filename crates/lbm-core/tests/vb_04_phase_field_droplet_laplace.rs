// VB-04 adversarial validation skeleton.
// Source of truth: docs/VALIDATION_BIOPROCESS.md#vb-04--phase-field-droplet-and-laplace-law

const VB04_IGNORE_REASON: &str =
    "VB-04: waits on Bundle P BCFD-048 validation suite; BCFD-046 sparger dependency also not landed";

const DROPLET_RADII_LU: [f64; 3] = [8.0, 16.0, 32.0];
const SURFACE_TENSION_SIGMA: f64 = 0.01;
const LAPLACE_RELATIVE_TOLERANCE: f64 = 0.10;
const SLOPE_CONSTANT_OFFSET_TOLERANCE: f64 = 0.05;
const TOTAL_PHI_DRIFT_TOLERANCE: f64 = 0.001;
const PHI_DRIFT_WINDOW_STEPS: usize = 10_000;
const TWO_DIMENSIONAL_LAPLACE_FACTOR: f64 = 1.0;
const THREE_DIMENSIONAL_LAPLACE_FACTOR: f64 = 2.0;

#[derive(Clone, Copy, Debug)]
struct LaplaceMeasurement {
    radius_lu: f64,
    pressure_jump: f64,
}

#[derive(Clone, Copy, Debug)]
struct PhiLedger {
    initial_total_phi: f64,
    final_total_phi: f64,
    steps: usize,
}

#[ignore = "VB-04: waits on Bundle P BCFD-048 validation suite; BCFD-046 sparger dependency also not landed"]
#[test]
fn static_2d_droplets_follow_laplace_pressure_jump() {
    let measurements = pending_laplace_measurements(TWO_DIMENSIONAL_LAPLACE_FACTOR);

    assert_laplace_pressure_jump_for_each_radius(&measurements, TWO_DIMENSIONAL_LAPLACE_FACTOR);
    assert_laplace_slope_through_origin(&measurements, TWO_DIMENSIONAL_LAPLACE_FACTOR);
}

#[ignore = "VB-04: waits on Bundle P BCFD-048 validation suite; BCFD-046 sparger dependency also not landed"]
#[test]
fn static_3d_droplets_follow_laplace_pressure_jump() {
    let measurements = pending_laplace_measurements(THREE_DIMENSIONAL_LAPLACE_FACTOR);

    assert_laplace_pressure_jump_for_each_radius(&measurements, THREE_DIMENSIONAL_LAPLACE_FACTOR);
    assert_laplace_slope_through_origin(&measurements, THREE_DIMENSIONAL_LAPLACE_FACTOR);
}

#[ignore = "VB-04: waits on Bundle P BCFD-048 validation suite; BCFD-046 sparger dependency also not landed"]
#[test]
fn total_phi_drift_stays_below_point_one_percent_for_10000_steps() {
    let ledgers = pending_phi_ledgers(PHI_DRIFT_WINDOW_STEPS);

    for ledger in &ledgers {
        assert_total_phi_drift_within_band(*ledger);
    }
}

fn pending_laplace_measurements(laplace_factor: f64) -> [LaplaceMeasurement; 3] {
    panic!(
        "{VB04_IGNORE_REASON}: run real static droplets for radii {DROPLET_RADII_LU:?}, \
         laplace_factor={laplace_factor}; no mocked phase-field data"
    )
}

fn pending_phi_ledgers(steps: usize) -> [PhiLedger; 3] {
    panic!("{VB04_IGNORE_REASON}: run real conservative Allen-Cahn droplets for steps={steps}")
}

fn assert_laplace_pressure_jump_for_each_radius(
    measurements: &[LaplaceMeasurement; 3],
    laplace_factor: f64,
) {
    for (measurement, expected_radius) in measurements.iter().zip(DROPLET_RADII_LU) {
        assert_eq!(
            measurement.radius_lu, expected_radius,
            "VB-04 requires droplet radii {DROPLET_RADII_LU:?}"
        );
        let expected_delta_p = laplace_factor * SURFACE_TENSION_SIGMA / measurement.radius_lu;
        let relative_error =
            (measurement.pressure_jump - expected_delta_p).abs() / expected_delta_p.abs();
        assert!(
            relative_error <= LAPLACE_RELATIVE_TOLERANCE,
            "VB-04 Laplace law: radius={}, measured_delta_p={}, expected_delta_p={}, \
             relative_error={relative_error}, tolerance={LAPLACE_RELATIVE_TOLERANCE}; \
             denominator is analytic sigma/R or 2sigma/R",
            measurement.radius_lu,
            measurement.pressure_jump,
            expected_delta_p
        );
    }
}

fn assert_laplace_slope_through_origin(
    measurements: &[LaplaceMeasurement; 3],
    laplace_factor: f64,
) {
    let expected_slope = laplace_factor * SURFACE_TENSION_SIGMA;
    let fitted_slope = slope_through_origin_for_delta_p_vs_inverse_radius(measurements);
    let implied_offsets = measurements.map(|m| m.pressure_jump - fitted_slope / m.radius_lu);
    let max_offset = implied_offsets
        .into_iter()
        .map(f64::abs)
        .fold(f64::NEG_INFINITY, f64::max);
    let relative_slope_error = (fitted_slope - expected_slope).abs() / expected_slope.abs();
    let offset_limit = SLOPE_CONSTANT_OFFSET_TOLERANCE * expected_slope / DROPLET_RADII_LU[0];

    assert!(
        relative_slope_error <= LAPLACE_RELATIVE_TOLERANCE,
        "VB-04 Laplace slope: fitted_slope={fitted_slope}, expected_slope={expected_slope}, \
         relative_error={relative_slope_error}, tolerance={LAPLACE_RELATIVE_TOLERANCE}"
    );
    assert!(
        max_offset <= offset_limit,
        "VB-04 Laplace line must pass through origin: max_constant_offset={max_offset}, \
         offset_limit={offset_limit}, tolerance_fraction={SLOPE_CONSTANT_OFFSET_TOLERANCE}"
    );
}

fn assert_total_phi_drift_within_band(ledger: PhiLedger) {
    assert_eq!(
        ledger.steps, PHI_DRIFT_WINDOW_STEPS,
        "VB-04 total-phi drift window must be {PHI_DRIFT_WINDOW_STEPS} steps"
    );
    let relative_drift =
        (ledger.final_total_phi - ledger.initial_total_phi).abs() / ledger.initial_total_phi.abs();
    assert!(
        relative_drift <= TOTAL_PHI_DRIFT_TOLERANCE,
        "VB-04 total-phi drift: initial={}, final={}, relative_drift={relative_drift}, \
         tolerance={TOTAL_PHI_DRIFT_TOLERANCE}; denominator is initial total phi",
        ledger.initial_total_phi,
        ledger.final_total_phi
    );
}

fn slope_through_origin_for_delta_p_vs_inverse_radius(
    measurements: &[LaplaceMeasurement; 3],
) -> f64 {
    let numerator: f64 = measurements
        .iter()
        .map(|m| (1.0 / m.radius_lu) * m.pressure_jump)
        .sum();
    let denominator: f64 = measurements
        .iter()
        .map(|m| (1.0 / m.radius_lu).powi(2))
        .sum();
    numerator / denominator
}
