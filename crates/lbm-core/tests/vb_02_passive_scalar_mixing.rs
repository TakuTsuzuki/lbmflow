// VB-02 adversarial validation skeleton.
// Source of truth: docs/VALIDATION_BIOPROCESS.md#vb-02--passive-scalar-mixing

const VB02_IGNORE_REASON: &str = "VB-02: waits on BCFD-034/035";

const CV95_FRACTION_OF_INITIAL: f64 = 0.05;
const CV99_FRACTION_OF_INITIAL: f64 = 0.01;
const TIMESTEP_INVARIANCE_TOLERANCE: f64 = 0.05;
const MIXING_TRANSIENT_SAMPLE_COUNT: usize = 8;
const PUBLISHED_NTHETA_MIN: f64 = 3.0;
const PUBLISHED_NTHETA_MAX: f64 = 7.0;
const BASE_TIMESTEP_DT: f64 = 1.0;
const HALVED_TIMESTEP_DT: f64 = 0.5;

#[derive(Clone, Debug)]
struct MixingRun {
    rotational_speed_hz: f64,
    dt: f64,
    initial_cv: f64,
    cv_history: Vec<(f64, f64)>,
    t95: f64,
    t99: f64,
}

#[ignore = "VB-02: waits on BCFD-034/035"]
#[test]
fn point_pulse_scalar_cv_decays_and_reports_t95_t99_n_theta() {
    let run = pending_point_pulse_mixing_run(BASE_TIMESTEP_DT);

    assert_cv_decay_monotonic_after_transient(&run);
    assert_t95_and_t99_cross_declared_cv_thresholds(&run);
    assert_dimensionless_n_theta_within_published_band(&run);
}

#[ignore = "VB-02: waits on BCFD-034/035"]
#[test]
fn mixing_time_is_invariant_to_halved_timestep() {
    let base = pending_point_pulse_mixing_run(BASE_TIMESTEP_DT);
    let halved = pending_point_pulse_mixing_run(HALVED_TIMESTEP_DT);

    assert_halved_timestep_agrees_with_base(&base, &halved);
}

fn pending_point_pulse_mixing_run(dt: f64) -> MixingRun {
    panic!(
        "{VB02_IGNORE_REASON}: run real passive-scalar ADE point-pulse mixing at dt={dt}; \
         no mocked fluid solver data"
    )
}

fn assert_cv_decay_monotonic_after_transient(run: &MixingRun) {
    for window in run.cv_history[MIXING_TRANSIENT_SAMPLE_COUNT..].windows(2) {
        let (time_a, cv_a) = window[0];
        let (time_b, cv_b) = window[1];
        assert!(
            cv_b <= cv_a,
            "VB-02 CV must decay monotonically after transient: time_a={time_a}, \
             cv_a={cv_a}, time_b={time_b}, cv_b={cv_b}"
        );
    }
}

fn assert_t95_and_t99_cross_declared_cv_thresholds(run: &MixingRun) {
    let cv95 = cv_at_or_after(&run.cv_history, run.t95);
    let cv99 = cv_at_or_after(&run.cv_history, run.t99);
    let t95_threshold = CV95_FRACTION_OF_INITIAL * run.initial_cv;
    let t99_threshold = CV99_FRACTION_OF_INITIAL * run.initial_cv;

    assert!(
        cv95 <= t95_threshold,
        "VB-02 t95 threshold: cv_at_t95={cv95}, threshold={t95_threshold}, \
         initial_cv={}, fraction={CV95_FRACTION_OF_INITIAL}",
        run.initial_cv
    );
    assert!(
        cv99 <= t99_threshold,
        "VB-02 t99 threshold: cv_at_t99={cv99}, threshold={t99_threshold}, \
         initial_cv={}, fraction={CV99_FRACTION_OF_INITIAL}",
        run.initial_cv
    );
    assert!(
        run.t95 <= run.t99,
        "VB-02 t95 must not exceed t99: t95={}, t99={}",
        run.t95,
        run.t99
    );
}

fn assert_dimensionless_n_theta_within_published_band(run: &MixingRun) {
    let n_theta = run.rotational_speed_hz * run.t95;
    assert!(
        (PUBLISHED_NTHETA_MIN..=PUBLISHED_NTHETA_MAX).contains(&n_theta),
        "VB-02 Ntheta=N*t95 outside published band: Ntheta={n_theta}, \
         band=[{PUBLISHED_NTHETA_MIN}, {PUBLISHED_NTHETA_MAX}], \
         N={}, t95={}",
        run.rotational_speed_hz,
        run.t95
    );
}

fn assert_halved_timestep_agrees_with_base(base: &MixingRun, halved: &MixingRun) {
    assert_eq!(base.dt, BASE_TIMESTEP_DT, "VB-02 base run uses declared dt");
    assert_eq!(halved.dt, HALVED_TIMESTEP_DT, "VB-02 refined run halves dt");

    assert_relative_agreement(
        base.t95,
        halved.t95,
        TIMESTEP_INVARIANCE_TOLERANCE,
        "VB-02 t95 timestep invariance",
    );
    assert_relative_agreement(
        base.t99,
        halved.t99,
        TIMESTEP_INVARIANCE_TOLERANCE,
        "VB-02 t99 timestep invariance",
    );
}

fn assert_relative_agreement(base: f64, refined: f64, tolerance: f64, label: &str) {
    let relative_difference = (base - refined).abs() / refined.abs();
    assert!(
        relative_difference <= tolerance,
        "{label}: base={base}, refined={refined}, relative_difference={relative_difference}, \
         tolerance={tolerance}; denominator is refined value"
    );
}

fn cv_at_or_after(history: &[(f64, f64)], target_time: f64) -> f64 {
    history
        .iter()
        .find(|(time, _)| *time >= target_time)
        .map(|(_, cv)| *cv)
        .expect("VB-02 CV history must include samples at or after the reported mixing time")
}
