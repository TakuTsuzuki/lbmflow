//! ACC-AUDIT P5: checkerboard / odd-even density-mode decay.
//!
//! This closes radar row 15 and pitfall-checklist #5. The test is deliberately
//! small but adversarial: the staggered `(pi, pi)` mode lives at the Brillouin
//! corner, so a weakly damped ghost mode can persist while ordinary smooth-flow
//! validation bands remain green.

mod common;

use common::metrics::{envelope_fit, monotonicity};
use lbm_core::compat::prelude::*;

const N: usize = 64;
const EPS: f64 = 1.0e-4;
const STEPS: usize = 500;
const SAMPLE_EVERY: usize = 50;
const LEAKAGE_LIMIT: f64 = 1.0e-12;
const ENVELOPE_WINDOW: usize = 3;
const ENVELOPE_MONOTONICITY_LIMIT: f64 = 0.9;

#[derive(Clone, Copy)]
struct Case {
    label: &'static str,
    tau: f64,
    collision: Collision,
    final_decay_factor_limit: f64,
}

#[derive(Debug)]
struct Sample {
    step: usize,
    amp_pi_pi_signed: f64,
    amp_pi_pi_abs: f64,
    leak_pi_0: f64,
    leak_0_pi: f64,
}

fn nu_from_tau(tau: f64) -> f64 {
    (tau - 0.5) / 3.0
}

fn checker_sign(x: usize, y: usize) -> f64 {
    if (x + y) & 1 == 0 {
        1.0
    } else {
        -1.0
    }
}

fn x_stagger_sign(x: usize) -> f64 {
    if x & 1 == 0 {
        1.0
    } else {
        -1.0
    }
}

fn y_stagger_sign(y: usize) -> f64 {
    if y & 1 == 0 {
        1.0
    } else {
        -1.0
    }
}

fn mode_projection(sim: &Simulation<f64>, sign: impl Fn(usize, usize) -> f64) -> f64 {
    let nx = sim.nx();
    let ny = sim.ny();
    let norm = (nx * ny) as f64;
    let mut sum = 0.0;
    for y in 0..ny {
        for x in 0..nx {
            sum += sim.rho(x, y) * sign(x, y);
        }
    }
    sum / norm
}

fn sample_modes(sim: &Simulation<f64>, step: usize) -> Sample {
    let signed = mode_projection(sim, checker_sign);
    Sample {
        step,
        amp_pi_pi_signed: signed,
        amp_pi_pi_abs: signed.abs(),
        leak_pi_0: mode_projection(sim, |x, _| x_stagger_sign(x)).abs(),
        leak_0_pi: mode_projection(sim, |_, y| y_stagger_sign(y)).abs(),
    }
}

fn make_sim(case: Case) -> Simulation<f64> {
    let mut sim = SimConfig {
        nx: N,
        ny: N,
        nu: nu_from_tau(case.tau),
        collision: case.collision,
        edges: Edges {
            left: EdgeBC::Periodic,
            right: EdgeBC::Periodic,
            bottom: EdgeBC::Periodic,
            top: EdgeBC::Periodic,
        },
        ..Default::default()
    }
    .build()
    .unwrap();
    sim.init_with(|x, y| (1.0 + EPS * checker_sign(x, y), 0.0, 0.0));
    sim
}

fn run_case(case: Case) -> Vec<Sample> {
    let mut sim = make_sim(case);
    let mut samples = vec![sample_modes(&sim, 0)];
    for step in (SAMPLE_EVERY..=STEPS).step_by(SAMPLE_EVERY) {
        sim.run(SAMPLE_EVERY);
        samples.push(sample_modes(&sim, step));
    }
    samples
}

fn rolling_max_envelope(xs: &[f64], window: usize) -> Vec<f64> {
    assert!(window > 0);
    (0..xs.len())
        .map(|i| {
            let start = (i + 1).saturating_sub(window);
            xs[start..=i].iter().copied().fold(0.0, f64::max)
        })
        .collect()
}

fn non_increasing_fraction(xs: &[f64]) -> f64 {
    assert!(xs.len() >= 2);
    let non_inc = xs.windows(2).filter(|w| w[1] <= w[0]).count();
    non_inc as f64 / (xs.len() - 1) as f64
}

fn assert_case(case: Case) {
    let samples = run_case(case);
    let abs_pi_pi: Vec<f64> = samples.iter().map(|s| s.amp_pi_pi_abs).collect();
    let envelope = rolling_max_envelope(&abs_pi_pi, ENVELOPE_WINDOW);
    let signed_pi_pi: Vec<(usize, f64)> = samples
        .iter()
        .map(|s| (s.step, s.amp_pi_pi_signed))
        .collect();
    let times: Vec<f64> = samples.iter().map(|s| s.step as f64).collect();
    let raw_monotonicity = monotonicity(&abs_pi_pi);
    let envelope_monotonicity = non_increasing_fraction(&envelope);
    let initial = abs_pi_pi[0];
    let final_amp = *abs_pi_pi.last().unwrap();
    let decay_factor = final_amp / initial;
    let fit = envelope_fit(&times, &abs_pi_pi);
    let lambda = -fit.slope;
    let max_leak = samples
        .iter()
        .map(|s| s.leak_pi_0 + s.leak_0_pi)
        .fold(0.0, f64::max);
    let max_growth = abs_pi_pi.iter().copied().fold(0.0, f64::max) / initial;

    println!(
        "P5 checkerboard {label}: tau={tau:.5}, samples={samples:?}, \
         signed_pi_pi={signed_pi_pi:?}, \
         raw_monotonicity={raw_monotonicity:.3}, envelope_window={ENVELOPE_WINDOW}, \
         envelope={envelope:?}, envelope_monotonicity={envelope_monotonicity:.3}, \
         max_growth={max_growth:.6e}, \
         decay_factor={decay_factor:.6e}, lambda_fit={lambda:.6e}, r2={r2:.6}, \
         max_leak_sum={max_leak:.6e}",
        label = case.label,
        tau = case.tau,
        r2 = fit.r2
    );

    // Linear reference derivation for the P5 checkerboard mode:
    //
    // Seed rho(x,y) = 1 + eps*(-1)^(x+y), u = 0 on an even periodic box. The
    // perturbation is the Fourier mode k = (pi, pi), because
    // exp(i*pi*x + i*pi*y) = (-1)^(x+y). For every D2Q9 velocity c_q, streaming
    // multiplies this Fourier component by exp(-i c_q.k). Since c_q has integer
    // components, exp(-i*pi*(c_qx+c_qy)) is exactly +/-1. Opposite pairs satisfy
    // c_q.k + c_opp.k = 0, so the pair phases are reciprocal and the BGK/TRT
    // Fourier update reduces to an exact collision-spectrum eigenproblem at this
    // k. Hydrodynamic modes are damped by viscosity; non-hydrodynamic "ghost"
    // eigenvectors are damped by their relaxation rates. In BGK that ghost rate
    // is tied to omega = 1/tau. In TRT the odd-parity ghost damping is set by
    // omega_minus, with Lambda = (1/omega_plus - 1/2)(1/omega_minus - 1/2). The
    // observable consequence is exponential decay of the checkerboard amplitude,
    // amp(t) ~= amp(0)*exp(-lambda_ghost*t), with no transfer to the orthogonal
    // staggered modes (pi,0) or (0,pi) on a symmetric periodic grid. Close to
    // round-off and near tau=0.5 the signed ghost-mode coefficient may oscillate;
    // the physical anchor is decay of its short-window absolute-amplitude envelope.
    assert!(
        (initial - EPS).abs() <= 2.0e-15,
        "{label} initial checkerboard projection changed: initial={initial:.12e}, \
         expected={EPS:.12e}, abs_tol=2e-15",
        label = case.label
    );
    assert!(
        envelope_monotonicity >= ENVELOPE_MONOTONICITY_LIMIT,
        "{label} checkerboard envelope is not sufficiently non-increasing: \
         envelope_monotonicity={envelope_monotonicity:.6}, \
         limit={ENVELOPE_MONOTONICITY_LIMIT:.6}, envelope_window={ENVELOPE_WINDOW}, \
         envelope={envelope:?}, raw_abs_samples={abs_pi_pi:?}, samples={samples:?}",
        label = case.label
    );
    assert!(
        max_growth < 2.0,
        "{label} checkerboard mode grew beyond the no-growth ceiling: \
         max_growth={max_growth:.12e}, ceiling=2.0, samples={samples:?}",
        label = case.label
    );
    assert!(
        max_leak < LEAKAGE_LIMIT,
        "{label} leaked into orthogonal staggered modes: max |amp(pi,0)|+|amp(0,pi)| \
         = {max_leak:.12e}, limit={LEAKAGE_LIMIT:.12e}, samples={samples:?}",
        label = case.label
    );
    assert!(
        decay_factor < case.final_decay_factor_limit,
        "{label} checkerboard floor too high after {STEPS} steps: final={final_amp:.12e}, \
         initial={initial:.12e}, decay_factor={decay_factor:.12e}, \
         limit={limit:.12e}, lambda_fit={lambda:.12e}, \
         samples={samples:?}",
        label = case.label,
        limit = case.final_decay_factor_limit
    );
}

#[test]
fn bgk_tau_0_6_checkerboard_mode_decays_without_leakage() {
    assert_case(Case {
        label: "BGK tau=0.6",
        tau: 0.6,
        collision: Collision::Bgk,
        final_decay_factor_limit: 1.0e-8,
    });
}

#[test]
fn trt_tau_0_51_checkerboard_mode_decays_without_leakage() {
    assert_case(Case {
        label: "TRT tau=0.51 magic=3/16",
        tau: 0.51,
        collision: Collision::Trt {
            magic: Collision::MAGIC_STD,
        },
        final_decay_factor_limit: 1.0e-3,
    });
}
