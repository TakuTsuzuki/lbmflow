// VB-06 adversarial validation skeleton.
// Source of truth: docs/VALIDATION_BIOPROCESS.md#vb-06--oxygen-kla-synthetic

const VB06_IGNORE_REASON: &str = "VB-06: waits on BCFD-050/051/052";

const INPUT_KLA_1_PER_S: f64 = 0.04;
const KLA_RELATIVE_TOLERANCE: f64 = 0.05;
const FIT_R2_MIN: f64 = 0.99;
const EQUILIBRIUM_KLA_ABSOLUTE_TOLERANCE: f64 = 1.0e-9;
const INITIAL_OXYGEN_C0: f64 = 0.0;
const SATURATION_OXYGEN_C_STAR: f64 = 1.0;
const FIXED_INTERFACIAL_AREA_A: f64 = 100.0;
const FIXED_KL: f64 = INPUT_KLA_1_PER_S / FIXED_INTERFACIAL_AREA_A;

#[derive(Clone, Copy, Debug)]
struct KlaFit {
    input_kla_1_per_s: f64,
    recovered_kla_1_per_s: f64,
    fit_r2: f64,
    c_initial: f64,
    c_star: f64,
    interfacial_area_a: f64,
    k_l: f64,
}

#[ignore = "VB-06: waits on BCFD-050/051/052"]
#[test]
fn dynamic_gassing_fit_recovers_input_kla_with_high_r2() {
    let fit = pending_dynamic_gassing_fit();

    assert_dynamic_gassing_setup_matches_synthetic_case(fit);
    assert_recovered_kla_within_5_percent_of_input(fit);
    assert_fit_r2_at_least_099(fit);
}

#[ignore = "VB-06: waits on BCFD-050/051/052"]
#[test]
fn equilibrium_case_fits_zero_kla() {
    let fit = pending_equilibrium_gassing_fit();

    assert_equilibrium_kla_is_zero_within_tolerance(fit);
}

fn pending_dynamic_gassing_fit() -> KlaFit {
    panic!(
        "{VB06_IGNORE_REASON}: run real oxygen scalar dynamic gassing with fixed a={FIXED_INTERFACIAL_AREA_A} \
         and kL={FIXED_KL}; no mocked solver data"
    )
}

fn pending_equilibrium_gassing_fit() -> KlaFit {
    panic!("{VB06_IGNORE_REASON}: run real equilibrium C=C* oxygen scalar case")
}

fn assert_dynamic_gassing_setup_matches_synthetic_case(fit: KlaFit) {
    assert_eq!(
        fit.c_initial, INITIAL_OXYGEN_C0,
        "VB-06 dynamic gassing starts from C0={INITIAL_OXYGEN_C0}"
    );
    assert_eq!(
        fit.c_star, SATURATION_OXYGEN_C_STAR,
        "VB-06 dynamic gassing uses C*={SATURATION_OXYGEN_C_STAR}"
    );
    assert_eq!(
        fit.interfacial_area_a, FIXED_INTERFACIAL_AREA_A,
        "VB-06 synthetic case uses fixed interfacial area a={FIXED_INTERFACIAL_AREA_A}"
    );
    assert_eq!(
        fit.k_l, FIXED_KL,
        "VB-06 synthetic case uses fixed kL={FIXED_KL}"
    );
}

fn assert_recovered_kla_within_5_percent_of_input(fit: KlaFit) {
    let relative_error =
        (fit.recovered_kla_1_per_s - fit.input_kla_1_per_s).abs() / fit.input_kla_1_per_s.abs();
    assert!(
        relative_error <= KLA_RELATIVE_TOLERANCE,
        "VB-06 kLa fit: recovered={}, input={}, relative_error={relative_error}, \
         tolerance={KLA_RELATIVE_TOLERANCE}; denominator is input kLa",
        fit.recovered_kla_1_per_s,
        fit.input_kla_1_per_s
    );
}

fn assert_fit_r2_at_least_099(fit: KlaFit) {
    assert!(
        fit.fit_r2 >= FIT_R2_MIN,
        "VB-06 kLa transient fit R^2: measured={}, minimum={FIT_R2_MIN}",
        fit.fit_r2
    );
}

fn assert_equilibrium_kla_is_zero_within_tolerance(fit: KlaFit) {
    assert!(
        fit.recovered_kla_1_per_s.abs() <= EQUILIBRIUM_KLA_ABSOLUTE_TOLERANCE,
        "VB-06 equilibrium case must fit kLa≈0: recovered={}, absolute_tolerance={EQUILIBRIUM_KLA_ABSOLUTE_TOLERANCE}",
        fit.recovered_kla_1_per_s
    );
}
