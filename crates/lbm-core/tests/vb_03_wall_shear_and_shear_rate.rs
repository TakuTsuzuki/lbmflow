// VB-03 adversarial validation skeleton.
// Source of truth: docs/VALIDATION_BIOPROCESS.md#vb-03--wall-shear-and-shear-rate-fields

const VB03_IGNORE_REASON: &str = "VB-03: waits on BCFD-032/033";

const COUETTE_WALL_SPEED_U: f64 = 0.02;
const CHANNEL_HEIGHT_H: f64 = 1.0;
const POISEUILLE_CENTERLINE_UMAX: f64 = 0.02;
const L2_ERROR_MAX_AT_N64: f64 = 1.0e-3;
const SECOND_ORDER_MIN_OBSERVED_ORDER: f64 = 1.8;
const SECOND_ORDER_MAX_OBSERVED_ORDER: f64 = 2.2;
const CONVERGENCE_RESOLUTIONS: [usize; 3] = [32, 64, 128];
const POISEUILLE_GRADIENT_COEFFICIENT: f64 = -6.0;

#[derive(Clone, Debug)]
struct ShearFieldRun {
    resolution: usize,
    samples: Vec<ShearSample>,
    l2_error: f64,
}

#[derive(Clone, Copy, Debug)]
struct ShearSample {
    y_from_center: f64,
    measured_du_dy: f64,
}

#[ignore = "VB-03: waits on BCFD-032/033"]
#[test]
fn couette_gamma_dot_matches_u_over_h_at_n64() {
    let run = pending_couette_shear_run(CONVERGENCE_RESOLUTIONS[1]);

    assert_l2_error_within_n64_band(&run);
    assert_couette_profile_matches_analytic(&run);
}

#[ignore = "VB-03: waits on BCFD-032/033"]
#[test]
fn poiseuille_gradient_matches_signed_analytic_profile_at_n64() {
    let run = pending_poiseuille_shear_run(CONVERGENCE_RESOLUTIONS[1]);

    assert_l2_error_within_n64_band(&run);
    assert_poiseuille_profile_matches_analytic(&run);
}

#[ignore = "VB-03: waits on BCFD-032/033"]
#[test]
fn shear_fields_converge_second_order_on_32_64_128() {
    let errors = pending_shear_convergence_errors(CONVERGENCE_RESOLUTIONS);

    assert_second_order_convergence(&errors);
}

fn pending_couette_shear_run(resolution: usize) -> ShearFieldRun {
    panic!(
        "{VB03_IGNORE_REASON}: run real Couette shear field at resolution={resolution}; \
         no mocked fluid solver data"
    )
}

fn pending_poiseuille_shear_run(resolution: usize) -> ShearFieldRun {
    panic!(
        "{VB03_IGNORE_REASON}: run real Poiseuille shear field at resolution={resolution}; \
         no mocked fluid solver data"
    )
}

fn pending_shear_convergence_errors(resolutions: [usize; 3]) -> [(usize, f64); 3] {
    panic!(
        "{VB03_IGNORE_REASON}: run real Couette/Poiseuille shear convergence for {resolutions:?}"
    )
}

fn assert_l2_error_within_n64_band(run: &ShearFieldRun) {
    assert_eq!(
        run.resolution, CONVERGENCE_RESOLUTIONS[1],
        "VB-03 L2 band is specified at N=64"
    );
    assert!(
        run.l2_error <= L2_ERROR_MAX_AT_N64,
        "VB-03 shear L2 error at N=64: measured={}, tolerance={L2_ERROR_MAX_AT_N64}",
        run.l2_error
    );
}

fn assert_couette_profile_matches_analytic(run: &ShearFieldRun) {
    let expected_du_dy = COUETTE_WALL_SPEED_U / CHANNEL_HEIGHT_H;
    for sample in &run.samples {
        assert_relative_agreement(
            sample.measured_du_dy,
            expected_du_dy,
            L2_ERROR_MAX_AT_N64,
            "VB-03 Couette du/dy=U/H",
        );
    }
}

fn assert_poiseuille_profile_matches_analytic(run: &ShearFieldRun) {
    for sample in &run.samples {
        let expected_du_dy =
            POISEUILLE_GRADIENT_COEFFICIENT * POISEUILLE_CENTERLINE_UMAX * sample.y_from_center
                / CHANNEL_HEIGHT_H.powi(2);
        assert_relative_agreement(
            sample.measured_du_dy,
            expected_du_dy,
            L2_ERROR_MAX_AT_N64,
            "VB-03 Poiseuille du/dy=-6*Umax*y/H^2",
        );
    }
}

fn assert_second_order_convergence(errors: &[(usize, f64); 3]) {
    let resolutions = [errors[0].0, errors[1].0, errors[2].0];
    assert_eq!(
        resolutions, CONVERGENCE_RESOLUTIONS,
        "VB-03 convergence must use {CONVERGENCE_RESOLUTIONS:?}"
    );
    let observed_order = (errors[0].1 / errors[2].1).log2()
        / (CONVERGENCE_RESOLUTIONS[2] as f64 / CONVERGENCE_RESOLUTIONS[0] as f64).log2();
    assert!(
        (SECOND_ORDER_MIN_OBSERVED_ORDER..=SECOND_ORDER_MAX_OBSERVED_ORDER)
            .contains(&observed_order),
        "VB-03 second-order convergence: observed_order={observed_order}, \
         expected_range=[{SECOND_ORDER_MIN_OBSERVED_ORDER}, {SECOND_ORDER_MAX_OBSERVED_ORDER}], \
         errors={errors:?}"
    );
}

fn assert_relative_agreement(measured: f64, expected: f64, tolerance: f64, label: &str) {
    let denominator = expected.abs().max(f64::EPSILON);
    let relative_error = (measured - expected).abs() / denominator;
    assert!(
        relative_error <= tolerance,
        "{label}: measured={measured}, expected={expected}, relative_error={relative_error}, \
         tolerance={tolerance}; denominator is analytic shear magnitude"
    );
}
