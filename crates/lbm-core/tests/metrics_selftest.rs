//! Self-tests for the shared agreement-metrics library (tests/common/metrics.rs).
//!
//! Fixtures are analytic, so every expected value below is exact up to
//! floating-point rounding — no simulation is run.

mod common;

use common::metrics::*;
use std::f64::consts::PI;

#[test]
fn l2_rel_and_linf_rel_on_known_perturbation() {
    let reference = [1.0, 2.0, 3.0, 4.0];
    // Perturb one entry by +0.3: L2 = 0.3/sqrt(30), Linf = 0.3/4.
    let actual = [1.0, 2.3, 3.0, 4.0];
    let l2 = l2_rel(&actual, &reference);
    assert!((l2 - 0.3 / 30.0f64.sqrt()).abs() < 1e-15, "l2={l2}");
    let li = linf_rel(&actual, &reference, 0.0);
    assert!((li - 0.3 / 4.0).abs() < 1e-15, "linf={li}");
    // Floor dominates a small reference.
    let li_floored = linf_rel(&[0.1], &[0.0], 1.0);
    assert!((li_floored - 0.1).abs() < 1e-15);
}

#[test]
fn order_fit_recovers_exact_second_order() {
    // err = 3 h^2 exactly -> slope 2, intercept ln 3, r2 = 1.
    let h = [0.1, 0.05, 0.025, 0.0125];
    let err: Vec<f64> = h.iter().map(|&x| 3.0 * x * x).collect();
    let fit = order_fit(&h, &err);
    assert!((fit.slope - 2.0).abs() < 1e-12, "slope={}", fit.slope);
    assert!((fit.intercept - 3.0f64.ln()).abs() < 1e-12);
    assert!(fit.r2 > 1.0 - 1e-12);
}

#[test]
fn order_fit_r2_flags_non_powerlaw_data() {
    // Error saturating at a plateau is NOT a power law; r2 must drop so a
    // test asserting r2 >= 0.98 rejects the fit instead of trusting the slope.
    let h = [0.1, 0.05, 0.025, 0.0125];
    let err = [3e-2, 8e-3, 5e-3, 4.8e-3];
    let fit = order_fit(&h, &err);
    assert!(fit.r2 < 0.98, "r2={} should flag the plateau", fit.r2);
}

#[test]
fn envelope_fit_recovers_exponential_decay() {
    // amp = 0.7 exp(-0.35 y) -> intercept ln 0.7, slope -0.35.
    let y: [f64; 5] = [0.5, 1.0, 2.0, 3.5, 5.0];
    let amp: Vec<f64> = y.iter().map(|&v| 0.7 * (-0.35 * v).exp()).collect();
    let fit = envelope_fit(&y, &amp);
    assert!((fit.slope + 0.35).abs() < 1e-12);
    assert!((fit.intercept.exp() - 0.7).abs() < 1e-12);
    assert!(fit.r2 > 1.0 - 1e-12);
}

#[test]
fn phase_fit_recovers_amplitude_and_phase() {
    // s = 0.02 sin(omega t + 0.6), sampled over exactly 4 periods.
    let omega = 2.0 * PI / 100.0;
    let t: Vec<f64> = (0..400).map(|i| i as f64).collect();
    let s: Vec<f64> = t.iter().map(|&ti| 0.02 * (omega * ti + 0.6).sin()).collect();
    let (amp, phase) = phase_fit(&t, &s, omega);
    assert!((amp - 0.02).abs() < 1e-6, "amp={amp}");
    assert!((phase - 0.6).abs() < 1e-6, "phase={phase}");
}

#[test]
fn monotonicity_counts_decreasing_pairs() {
    assert_eq!(monotonicity(&[4.0, 3.0, 2.0, 1.0]), 1.0);
    assert_eq!(monotonicity(&[4.0, 3.0, 3.5, 1.0]), 2.0 / 3.0);
    assert_eq!(monotonicity(&[1.0, 2.0]), 0.0);
}

#[test]
fn curve_agreement_localizes_the_worst_deviation() {
    // theory: y = x^2; samples exact except x=3 which is 10% off.
    let samples = [(1.0, 1.0), (2.0, 4.0), (3.0, 9.9), (4.0, 16.0)];
    let agree = curve_agreement(|x| x * x, &samples, 0.05, 0.0);
    assert!((agree.max_rel_dev - 0.1).abs() < 1e-12);
    assert_eq!(agree.worst_x, 3.0);
    assert!((agree.frac_in_band - 0.75).abs() < 1e-12);
}
