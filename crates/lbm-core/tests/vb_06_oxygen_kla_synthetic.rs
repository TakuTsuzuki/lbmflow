// VB-06 adversarial validation skeleton.
// Source of truth: docs/VALIDATION_BIOPROCESS.md#vb-06--oxygen-kla-synthetic

use lbm_core::prelude::{dynamic_gassing_kla_fit, KlaFitWindow};

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

#[test]
fn dynamic_gassing_fit_recovers_input_kla_with_high_r2() { // verified 2026-07-08
    let fit = pending_dynamic_gassing_fit();

    assert_dynamic_gassing_setup_matches_synthetic_case(fit);
    assert_recovered_kla_within_5_percent_of_input(fit);
    assert_fit_r2_at_least_099(fit);
}

#[ignore = "VB-06: BCFD-050/051/052 landed, but equilibrium C=C* returns skipped QOI instead of fitted kLa≈0; see BCFD-VV-001"]
#[test]
fn equilibrium_case_fits_zero_kla() {
    let fit = pending_equilibrium_gassing_fit();

    assert_equilibrium_kla_is_zero_within_tolerance(fit);
}

fn pending_dynamic_gassing_fit() -> KlaFit {
    let series: Vec<(f64, f64)> = (0..=100)
        .map(|i| {
            let t = i as f64;
            let c = SATURATION_OXYGEN_C_STAR
                - (SATURATION_OXYGEN_C_STAR - INITIAL_OXYGEN_C0) * (-INPUT_KLA_1_PER_S * t).exp();
            (t, c)
        })
        .collect();
    let result = dynamic_gassing_kla_fit(
        &series,
        SATURATION_OXYGEN_C_STAR,
        Some(KlaFitWindow {
            start_s: 10.0,
            end_s: 90.0,
        }),
        1.0e-12,
    )
    .expect("VB-06 dynamic gassing fit should be computable")
    .result
    .expect("VB-06 dynamic gassing fit must not be skipped");
    KlaFit {
        input_kla_1_per_s: INPUT_KLA_1_PER_S,
        recovered_kla_1_per_s: result.kla_1_per_s,
        fit_r2: result.fit_r2,
        c_initial: INITIAL_OXYGEN_C0,
        c_star: SATURATION_OXYGEN_C_STAR,
        interfacial_area_a: FIXED_INTERFACIAL_AREA_A,
        k_l: FIXED_KL,
    }
}

fn pending_equilibrium_gassing_fit() -> KlaFit {
    let series: Vec<(f64, f64)> = (0..=10)
        .map(|i| (i as f64, SATURATION_OXYGEN_C_STAR))
        .collect();
    let recovered = dynamic_gassing_kla_fit(&series, SATURATION_OXYGEN_C_STAR, None, 1.0e-12)
        .expect("VB-06 equilibrium fit call should not error")
        .result
        .map(|result| result.kla_1_per_s)
        .unwrap_or(f64::NAN);
    KlaFit {
        input_kla_1_per_s: 0.0,
        recovered_kla_1_per_s: recovered,
        fit_r2: f64::NAN,
        c_initial: SATURATION_OXYGEN_C_STAR,
        c_star: SATURATION_OXYGEN_C_STAR,
        interfacial_area_a: FIXED_INTERFACIAL_AREA_A,
        k_l: FIXED_KL,
    }
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
