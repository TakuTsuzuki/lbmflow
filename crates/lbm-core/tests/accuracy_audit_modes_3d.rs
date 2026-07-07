//! ACC-AUDIT P5-3D: Zanetti checkerboard / odd-even density-mode decay.
//!
//! This is the 3D extension of `accuracy_audit_modes.rs`: the seeded
//! Brillouin-corner mode is `k = (pi, pi, pi)`. D3Q19 and D3Q27 have different
//! ghost-mode spectra at this corner because D3Q27 carries the body diagonals
//! that D3Q19 omits, so both lattices are pinned explicitly.

mod common;

use common::metrics::{envelope_fit, monotonicity};
use lbm_core::prelude::*;

const N: usize = 32;
const EPS: f64 = 1.0e-4;
const STEPS: usize = 500;
const SAMPLE_EVERY: usize = 50;
const LEAKAGE_LIMIT: f64 = 1.0e-12;
const ENVELOPE_WINDOW: usize = 3;

type CpuPeriodic<L> = Solver<L, f64, CpuScalar, LocalPeriodic>;

#[derive(Clone, Copy)]
struct Case {
    label: &'static str,
    tau: f64,
    collision: CollisionKind,
    final_decay_factor_limit: f64,
}

#[derive(Clone, Copy, Debug)]
struct Sample {
    step: usize,
    amp_pi_pi_pi_signed: f64,
    amp_pi_pi_pi_abs: f64,
    leak_pi_pi_0: f64,
    leak_pi_0_pi: f64,
    leak_0_pi_pi: f64,
}

fn nu_from_tau(tau: f64) -> f64 {
    (tau - 0.5) / 3.0
}

fn stagger_sign(sum: usize) -> f64 {
    if sum & 1 == 0 {
        1.0
    } else {
        -1.0
    }
}

fn checker_sign(x: usize, y: usize, z: usize) -> f64 {
    stagger_sign(x + y + z)
}

fn mode_projection(rho: &[f64], sign: impl Fn(usize, usize, usize) -> f64) -> f64 {
    let norm = (N * N * N) as f64;
    let mut sum = 0.0;
    for z in 0..N {
        for y in 0..N {
            for x in 0..N {
                let i = (z * N + y) * N + x;
                sum += rho[i] * sign(x, y, z);
            }
        }
    }
    sum / norm
}

fn sample_modes<L: Lattice>(sim: &CpuPeriodic<L>, step: usize) -> Sample {
    let rho = sim.gather_rho();
    let signed = mode_projection(&rho, checker_sign);
    Sample {
        step,
        amp_pi_pi_pi_signed: signed,
        amp_pi_pi_pi_abs: signed.abs(),
        leak_pi_pi_0: mode_projection(&rho, |x, y, _| stagger_sign(x + y)).abs(),
        leak_pi_0_pi: mode_projection(&rho, |x, _, z| stagger_sign(x + z)).abs(),
        leak_0_pi_pi: mode_projection(&rho, |_, y, z| stagger_sign(y + z)).abs(),
    }
}

fn make_sim<L: Lattice>(case: Case) -> CpuPeriodic<L> {
    let spec = GlobalSpec::<f64> {
        dims: [N, N, N],
        nu: nu_from_tau(case.tau),
        periodic: [true, true, true],
        collision: case.collision,
        ..Default::default()
    };
    let mut sim: CpuPeriodic<L> = Solver::new(
        &spec,
        &[],
        &[],
        [1, 1, 1],
        CpuScalar::default(),
        LocalPeriodic,
    );
    sim.init_with(|x, y, z| (1.0 + EPS * checker_sign(x, y, z), [0.0; 3]));
    sim
}

fn run_case<L: Lattice>(case: Case) -> Vec<Sample> {
    let mut sim = make_sim::<L>(case);
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

fn assert_case<L: Lattice>(case: Case) {
    let samples = run_case::<L>(case);
    let abs_corner: Vec<f64> = samples.iter().map(|s| s.amp_pi_pi_pi_abs).collect();
    let envelope = rolling_max_envelope(&abs_corner, ENVELOPE_WINDOW);
    let signed_corner: Vec<(usize, f64)> = samples
        .iter()
        .map(|s| (s.step, s.amp_pi_pi_pi_signed))
        .collect();
    let times: Vec<f64> = samples.iter().map(|s| s.step as f64).collect();
    let raw_monotonicity = monotonicity(&abs_corner);
    let envelope_monotonicity = non_increasing_fraction(&envelope);
    let initial = abs_corner[0];
    let final_amp = *abs_corner.last().unwrap();
    let decay_factor = final_amp / initial;
    let positive_fit_points: Vec<(f64, f64)> = times
        .iter()
        .copied()
        .zip(abs_corner.iter().copied())
        .filter(|(_, amp)| *amp > 0.0)
        .collect();
    let positive_times: Vec<f64> = positive_fit_points.iter().map(|(t, _)| *t).collect();
    let positive_amp: Vec<f64> = positive_fit_points.iter().map(|(_, amp)| *amp).collect();
    let fit = envelope_fit(&positive_times, &positive_amp);
    let lambda = -fit.slope;
    let max_leak = samples
        .iter()
        .map(|s| s.leak_pi_pi_0.max(s.leak_pi_0_pi).max(s.leak_0_pi_pi))
        .fold(0.0, f64::max);
    let max_growth = abs_corner.iter().copied().fold(0.0, f64::max) / initial;

    println!(
        "P5-3D checkerboard {label}: lattice={lattice}, tau={tau:.5}, \
         samples={samples:?}, signed_pi_pi_pi={signed_corner:?}, \
         raw_monotonicity={raw_monotonicity:.3}, envelope_window={ENVELOPE_WINDOW}, \
         envelope={envelope:?}, envelope_monotonicity={envelope_monotonicity:.3}, \
         max_growth={max_growth:.6e}, decay_factor={decay_factor:.6e}, \
         lambda_fit={lambda:.6e}, r2={r2:.6}, max_planar_leak={max_leak:.6e}",
        label = case.label,
        lattice = std::any::type_name::<L>(),
        tau = case.tau,
        r2 = fit.r2
    );

    // Linear reference derivation for the P5-3D checkerboard mode:
    //
    // Seed rho(x,y,z) = 1 + eps*(-1)^(x+y+z), u = 0 on an even periodic box.
    // The perturbation is the Fourier mode k = (pi, pi, pi). Streaming
    // multiplies each population by exp(-i c_q.k), exactly +/-1 for integer
    // D3Q19/D3Q27 velocities. The periodic, zero-velocity setup preserves the
    // Brillouin-corner parity class, so a correct collision/streaming sequence
    // damps the corner ghost mode without transferring energy into the planar
    // corner modes (pi,pi,0), (pi,0,pi), or (0,pi,pi). D3Q19 and D3Q27 can
    // damp at different rates because the D3Q27 body diagonals add a distinct
    // phase subset at this k. Near roundoff, signed samples may oscillate; the
    // pinned behavior is non-increase of the short-window absolute envelope.
    assert!(
        (initial - EPS).abs() <= 5.0e-15,
        "{label} {lattice} initial checkerboard projection changed: \
         initial={initial:.12e}, expected={EPS:.12e}, abs_tol=5e-15",
        label = case.label,
        lattice = std::any::type_name::<L>()
    );
    assert!(
        envelope_monotonicity == 1.0,
        "{label} {lattice} checkerboard envelope is not non-increasing: \
         envelope_monotonicity={envelope_monotonicity:.6}, required=1.0, \
         envelope_window={ENVELOPE_WINDOW}, envelope={envelope:?}, \
         raw_abs_samples={abs_corner:?}, samples={samples:?}",
        label = case.label,
        lattice = std::any::type_name::<L>()
    );
    assert!(
        max_growth <= 2.0,
        "{label} {lattice} checkerboard mode grew beyond the no-growth ceiling: \
         max_growth={max_growth:.12e}, ceiling=2.0, samples={samples:?}",
        label = case.label,
        lattice = std::any::type_name::<L>()
    );
    assert!(
        max_leak <= LEAKAGE_LIMIT,
        "{label} {lattice} leaked into planar-corner modes: \
         max_planar_leak={max_leak:.12e}, limit={LEAKAGE_LIMIT:.12e}, \
         samples={samples:?}",
        label = case.label,
        lattice = std::any::type_name::<L>()
    );
    assert!(
        decay_factor < case.final_decay_factor_limit,
        "{label} {lattice} checkerboard floor too high after {STEPS} steps: \
         final={final_amp:.12e}, initial={initial:.12e}, \
         decay_factor={decay_factor:.12e}, limit={limit:.12e}, \
         lambda_fit={lambda:.12e}, samples={samples:?}",
        label = case.label,
        lattice = std::any::type_name::<L>(),
        limit = case.final_decay_factor_limit
    );
}

fn bgk_case() -> Case {
    Case {
        label: "BGK tau=0.6",
        tau: 0.6,
        collision: CollisionKind::Bgk,
        final_decay_factor_limit: 1.0e-8,
    }
}

fn trt_case() -> Case {
    Case {
        label: "TRT tau=0.51 magic=3/16",
        tau: 0.51,
        collision: CollisionKind::Trt {
            magic: CollisionKind::MAGIC_STD,
        },
        final_decay_factor_limit: 1.0e-3,
    }
}

#[test]
fn d3q19_bgk_tau_0_6_checkerboard_mode_decays_without_leakage() {
    assert_case::<D3Q19>(bgk_case());
}

#[test]
fn d3q19_trt_tau_0_51_checkerboard_mode_decays_without_leakage() {
    assert_case::<D3Q19>(trt_case());
}

#[test]
fn d3q27_bgk_tau_0_6_checkerboard_mode_decays_without_leakage() {
    assert_case::<D3Q27>(bgk_case());
}

#[test]
fn d3q27_trt_tau_0_51_checkerboard_mode_decays_without_leakage() {
    assert_case::<D3Q27>(trt_case());
}
