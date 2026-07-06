//! Accuracy audit for the cumulant/central-moment viscosity-defect correction.
//!
//! These tests audit the closure property, not the implementation constant:
//! a valid correction must remove the finite-frame lattice viscosity defect
//! across velocity amplitude, resolution, and orientation. Failing assertions
//! are triage findings, not invitations to loosen bands.

mod common;

use common::metrics::{linear_fit, LinFit};
use common::tgv_analysis::{energy_decay_rate, ke3d, tgv_nu_eff};
use lbm_core::prelude::*;
use std::f64::consts::PI;

type S3<L> = Solver<L, f64, CpuScalar, LocalPeriodic>;

const NU: f64 = 0.02;
const TRT_MAGIC: f64 = 3.0 / 16.0;

#[derive(Clone, Copy)]
enum TgvPlane {
    Xy,
    Yz,
    Zx,
}

fn omega_from_nu(nu: f64) -> f64 {
    1.0 / (3.0 * nu + 0.5)
}

fn periodic_spec<L: Lattice>(n: usize, nu: f64, collision: CollisionKind) -> GlobalSpec<f64> {
    let _ = L::D;
    GlobalSpec {
        dims: [n, n, n],
        nu,
        collision,
        periodic: [true, true, true],
        ..Default::default()
    }
}

fn pressure_consistent_rho(n: usize, u0: f64, x: usize, y: usize, z: usize) -> f64 {
    let k = 2.0 * PI / n as f64;
    let (xf, yf, zf) = (k * x as f64, k * y as f64, k * z as f64);
    // For the classic 3D TGV, the incompressible pressure satisfying
    // grad p = -u·grad u at t=0 is
    // p = u0^2/16 * (cos 2x + cos 2y) * (cos 2z + 2).
    // LBM uses p = cs^2 (rho - 1), cs^2 = 1/3, so rho = 1 + 3p.
    let p = u0 * u0 / 16.0 * (((2.0 * xf).cos() + (2.0 * yf).cos()) * ((2.0 * zf).cos() + 2.0));
    1.0 + 3.0 * p
}

fn tgv_velocity(n: usize, u0: f64, plane: TgvPlane, x: usize, y: usize, z: usize) -> [f64; 3] {
    let k = 2.0 * PI / n as f64;
    let (xf, yf, zf) = (k * x as f64, k * y as f64, k * z as f64);
    match plane {
        TgvPlane::Xy => [
            u0 * xf.sin() * yf.cos() * zf.cos(),
            -u0 * xf.cos() * yf.sin() * zf.cos(),
            0.0,
        ],
        TgvPlane::Yz => [
            0.0,
            u0 * yf.sin() * zf.cos() * xf.cos(),
            -u0 * yf.cos() * zf.sin() * xf.cos(),
        ],
        TgvPlane::Zx => [
            -u0 * zf.cos() * xf.sin() * yf.cos(),
            0.0,
            u0 * zf.sin() * xf.cos() * yf.cos(),
        ],
    }
}

fn make_tgv<L: Lattice>(
    n: usize,
    nu: f64,
    u0: f64,
    collision: CollisionKind,
    plane: TgvPlane,
) -> S3<L> {
    let spec = periodic_spec::<L>(n, nu, collision);
    let mut s: S3<L> = Solver::new(
        &spec,
        &[],
        &[],
        [1, 1, 1],
        CpuScalar::default(),
        LocalPeriodic,
    );
    s.init_with(move |x, y, z| {
        (
            pressure_consistent_rho(n, u0, x, y, z),
            tgv_velocity(n, u0, plane, x, y, z),
        )
    });
    s
}

fn kinetic_energy<L: Lattice>(s: &S3<L>) -> f64 {
    ke3d(&s.gather_ux(), &s.gather_uy(), &s.gather_uz())
}

fn sample_steps(n: usize, nu: f64) -> (usize, usize) {
    let k = 2.0 * PI / n as f64;
    let t1 = (0.02 / (nu * k * k)).round() as usize;
    let t2 = (0.06 / (nu * k * k)).round() as usize;
    assert!(t2 > t1);
    (t1, t2)
}

fn measure_nu_eff<L: Lattice>(
    n: usize,
    nu: f64,
    u0: f64,
    collision: CollisionKind,
    plane: TgvPlane,
) -> f64 {
    let (t1, t2) = sample_steps(n, nu);
    let mut s = make_tgv::<L>(n, nu, u0, collision, plane);
    s.run(t1);
    let ke1 = kinetic_energy(&s);
    s.run(t2 - t1);
    let ke2 = kinetic_energy(&s);
    assert!(
        ke1.is_finite() && ke1 > 0.0 && ke2.is_finite() && ke2 > 0.0,
        "TGV energy must be finite and positive: ke1={ke1:e}, ke2={ke2:e}"
    );
    let k = 2.0 * PI / n as f64;
    let k2_sum = 3.0 * k * k;
    let rate = energy_decay_rate(ke1, ke2, (t2 - t1) as f64);
    let nu_eff = tgv_nu_eff(ke1, ke2, k2_sum, (t2 - t1) as f64);
    assert!(
        rate.is_finite() && nu_eff.is_finite(),
        "TGV decay observables must be finite: rate={rate:e}, nu_eff={nu_eff:e}"
    );
    nu_eff
}

fn cumulant_for_nu(nu: f64) -> CollisionKind {
    CollisionKind::Cumulant {
        omega_shear: omega_from_nu(nu),
    }
}

fn cumulant() -> CollisionKind {
    cumulant_for_nu(NU)
}

fn trt() -> CollisionKind {
    CollisionKind::Trt { magic: TRT_MAGIC }
}

fn fit_defect_vs_u2<L: Lattice>(
    n: usize,
    nu: f64,
    u0s: &[f64],
    collision: CollisionKind,
) -> (LinFit, Vec<f64>) {
    let x: Vec<f64> = u0s.iter().map(|u0| u0 * u0).collect();
    let y: Vec<f64> = u0s
        .iter()
        .map(|&u0| measure_nu_eff::<L>(n, nu, u0, collision, TgvPlane::Xy) / nu - 1.0)
        .collect();
    let fit = linear_fit(&x, &y);
    (fit, y)
}

fn correction_tau_fingerprint(nu1: f64, nu2: f64) -> f64 {
    let omega1 = omega_from_nu(nu1);
    let omega2 = omega_from_nu(nu2);
    0.16 * (2.0 / (2.0 - omega1) - 2.0 / (2.0 - omega2))
}

fn assert_e1_tau_fingerprint(n: usize, u0s: &[f64], label: &str) {
    // Tau-fingerprint derivation. The cumulant shear-rate correction is
    // applied in omega space. With the lattice relation
    //
    //     nu = (1/omega - 1/2) / 3,
    //
    // a small omega perturbation gives
    //
    //     dnu = -(1 / (3 omega^2)) domega
    //     dnu/nu = [-domega/(3 omega^2)] / [(1/omega - 1/2)/3]
    //            = -(domega/omega) * 2/(2 - omega).
    //
    // Therefore a relative omega perturbation
    //
    //     domega/omega = -0.16 u^2
    //
    // contributes
    //
    //     dnu/nu = +0.16 u^2 * 2/(2 - omega).
    //
    // The measured slope c(tau) in
    //
    //     nu_eff/nu - 1 = c(tau) u0^2 + b
    //
    // must carry this known tau fingerprint if the correction reaches the
    // shear viscosity. Compressibility error and the intrinsic cubic-frame
    // viscosity defect are tau-independent at leading order, so a two-tau
    // Cumulant-minus-TRT double difference isolates the omega-space correction
    // without treating the raw u0^2 slope itself as proof.
    let nu1 = 0.02;
    let nu2 = 0.10;
    let omega1 = omega_from_nu(nu1);
    let omega2 = omega_from_nu(nu2);

    let (cum1_fit, cum1_y) = fit_defect_vs_u2::<D3Q19>(n, nu1, u0s, cumulant_for_nu(nu1));
    let (cum2_fit, cum2_y) = fit_defect_vs_u2::<D3Q19>(n, nu2, u0s, cumulant_for_nu(nu2));
    let (trt1_fit, trt1_y) = fit_defect_vs_u2::<D3Q19>(n, nu1, u0s, trt());
    let (trt2_fit, trt2_y) = fit_defect_vs_u2::<D3Q19>(n, nu2, u0s, trt());

    for (name, fit) in [
        ("cum_nu0.02", cum1_fit),
        ("cum_nu0.10", cum2_fit),
        ("trt_nu0.02", trt1_fit),
        ("trt_nu0.10", trt2_fit),
    ] {
        assert!(
            fit.r2 >= 0.99,
            "ACC CUM E1 {label}: {name} u0^2 slope fit r2={:.6e} below 0.99; slope={:.6e}, intercept={:.6e}, N={n}, u0={u0s:?}",
            fit.r2,
            fit.slope,
            fit.intercept
        );
    }

    let cum_tau_delta = cum1_fit.slope - cum2_fit.slope;
    let trt_tau_delta = trt1_fit.slope - trt2_fit.slope;
    let measured = cum_tau_delta - trt_tau_delta;
    let predicted = correction_tau_fingerprint(nu1, nu2);
    let allowed = 0.30 * predicted.abs();
    println!(
        "ACC CUM E1 {label}: N={n} u0={u0s:?} omega=({omega1:e},{omega2:e}) fingerprint_factors=({:e},{:e}) Cumulant rel nu0.02={cum1_y:?} c={:e} b={:e} r2={:e}; Cumulant rel nu0.10={cum2_y:?} c={:e} b={:e} r2={:e}; TRT rel nu0.02={trt1_y:?} c={:e} b={:e} r2={:e}; TRT rel nu0.10={trt2_y:?} c={:e} b={:e} r2={:e}; tau_delta_control_trt={trt_tau_delta:e}; measured_double_diff={measured:e}; predicted_correction={predicted:e}; band_abs={allowed:e}",
        2.0 / (2.0 - omega1),
        2.0 / (2.0 - omega2),
        cum1_fit.slope,
        cum1_fit.intercept,
        cum1_fit.r2,
        cum2_fit.slope,
        cum2_fit.intercept,
        cum2_fit.r2,
        trt1_fit.slope,
        trt1_fit.intercept,
        trt1_fit.r2,
        trt2_fit.slope,
        trt2_fit.intercept,
        trt2_fit.r2,
    );
    assert!(
        (measured - predicted).abs() <= allowed,
        "ACC CUM E1 {label}: tau-fingerprint double difference measured={measured:.6e} differs from correction prediction {predicted:.6e} by more than +/-30% ({allowed:.6e}); Cumulant delta={cum_tau_delta:.6e}, TRT tau-dependence control={trt_tau_delta:.6e}, omega=({omega1:.6e},{omega2:.6e}), normalization=d(nu_eff/nu-1)/d(u0^2)"
    );
}

// SPEC-GAP (triage 2026-07-06, ANOM-P4-008): the tau-fingerprint double
// difference cannot isolate the -0.16 u^2 omega modulation from the public
// API, for two measured reasons (N=48 data, kept for the record):
//   (1) the intrinsic cubic-defect slope c(tau) is strongly tau-dependent
//       (TRT control: c = 18.90 at nu=0.02 vs 0.207 at nu=0.10, delta 18.69),
//       so the "tau-independent baseline" assumption of the double
//       difference is wrong at the size of the effect being sought;
//   (2) the correction acts on the LOCAL u^2, so its nu_eff footprint is
//       weighted by the spatial average <u^2> = u0^2/4 on the classic TGV —
//       the prediction 0.16*(2/(2-w1) - 2/(2-w2)) = 1.067 overstates the
//       observable by ~4x (correct ~0.267); measured cum-TRT single
//       differences were +0.057 (nu=0.02) and +0.081 (nu=0.10), i.e. the
//       u^2 modulation shows NO clear nu_eff footprint at either predicted
//       size, but operator-intrinsic c differences confound the residual.
// Deciding the u^2 term needs a core-side ablation toggle (not exposed).
// The offset (+0.0025) verdict does NOT depend on this row - see E2/E3.
#[test]
#[ignore = "SPEC-GAP: isolating the -0.16 u^2 omega modulation requires a core-side ablation toggle; measured N=48 data recorded in the comment above"]
fn e1_cumulant_omega_space_correction_tau_fingerprint_light() {
    assert_e1_tau_fingerprint(32, &[0.02, 0.04, 0.08], "light");
}

#[test]
#[ignore = "SPEC-GAP: same as the light row; heavy N=48 variant kept for the post-ablation rerun"]
fn e1_cumulant_omega_space_correction_tau_fingerprint_heavy() {
    assert_e1_tau_fingerprint(48, &[0.02, 0.04, 0.08], "heavy");
}

fn h2_extrapolated_offset<L: Lattice>(
    ns: &[usize],
    collision_for_nu: impl Fn(f64) -> CollisionKind,
) -> (f64, f64, Vec<(usize, f64)>) {
    // Resolution-separation derivation. The TGV decay-rate fit inherits the
    // O(h^2) spatial discretization error of the second-order lattice
    // discretization, so the raw finite-N defect is
    //
    //     d(N) = nu_eff(N)/nu - 1 = a + b h^2 + higher order
    //
    // and h is proportional to 1/N for the fixed periodic domain. Fitting
    // d(N) = a + b/N^2 separates the resolution-independent closure offset
    // residual a from the spatial-error floor. A finite-N band alone would
    // confound these two terms, especially at N=24.
    let nu = NU;
    let samples: Vec<(f64, f64)> = ns
        .iter()
        .map(|&n| {
            let u0 = 1.28e-4 / n as f64;
            let nu_eff = measure_nu_eff::<L>(n, nu, u0, collision_for_nu(nu), TgvPlane::Xy);
            (1.0 / ((n * n) as f64), nu_eff / nu - 1.0)
        })
        .collect();
    let x: Vec<f64> = samples.iter().map(|&(x, _)| x).collect();
    let y: Vec<f64> = samples.iter().map(|&(_, d)| d).collect();
    let fit = linear_fit(&x, &y);
    let defects_by_n: Vec<(usize, f64)> = ns.iter().copied().zip(y).collect();
    (fit.intercept, fit.slope, defects_by_n)
}

fn assert_e2_h2_extrapolated_d3q19_offset(ns: &[usize], band: f64, label: &str) -> f64 {
    let (a, b, defects_by_n) = h2_extrapolated_offset::<D3Q19>(ns, |nu| cumulant_for_nu(nu));
    println!(
        "ACC CUM E2 {label}: N={ns:?} defects_nu_eff_over_nu_minus_1={defects_by_n:?}; h2_intercept_a={a:e}; h2_slope_b={b:e}; band_abs={band:e}"
    );
    assert!(
        a.abs() <= band,
        "ACC CUM E2 {label}: h^2-extrapolated D3Q19 offset residual |a|={:.6e} exceeds band {band:.6e}; a={a:.6e}, b={b:.6e}, defects_by_N={defects_by_n:?}, normalization=resolution-independent intercept of nu_eff/nu - 1",
        a.abs()
    );
    a
}

#[test]
fn e2_d3q19_h2_extrapolated_offset_residual_canary_light() {
    assert_e2_h2_extrapolated_d3q19_offset(&[24, 32], 4.0e-3, "light-canary");
}

#[test]
#[ignore = "heavy h^2-extrapolated offset audit: D3Q19 cumulant at N=24,32,48"]
fn e2_d3q19_h2_extrapolated_offset_residual_heavy() {
    assert_e2_h2_extrapolated_d3q19_offset(&[24, 32, 48], 2.0e-3, "heavy");
}

#[test]
fn e3_d3q27_h2_extrapolated_intrinsic_bias_canary_light() {
    assert_e3_h2_extrapolated_d3q27_bias(&[24, 32], 4.0e-3, "light-canary");
}

#[test]
#[ignore = "heavy h^2-extrapolated intrinsic-bias audit: D3Q27 cumulant at N=24,32,48"]
fn e3_d3q27_h2_extrapolated_intrinsic_bias_heavy() {
    assert_e3_h2_extrapolated_d3q27_bias(&[24, 32, 48], 2.0e-3, "heavy");
}

fn assert_e3_h2_extrapolated_d3q27_bias(ns: &[usize], band: f64, label: &str) -> f64 {
    // Cross-lattice control after removing the spatial-error confound. D3Q27
    // has the full tensor-product velocity set and the cumulant offset is zero
    // by construction, so the "D3Q19-only offset" explanation predicts a
    // near-zero h^2 intercept a27. A raw finite-N defect, such as the N=32
    // value, cannot classify operator bias until the O(h^2) spatial floor is
    // extrapolated away.
    let (a27, b27, defects27_by_n) = h2_extrapolated_offset::<D3Q27>(ns, |nu| cumulant_for_nu(nu));
    let (a19, b19, defects19_by_n) = h2_extrapolated_offset::<D3Q19>(ns, |nu| cumulant_for_nu(nu));
    println!(
        "ACC CUM E3 {label}: N={ns:?} D3Q27 defects_nu_eff_over_nu_minus_1={defects27_by_n:?}; a27={a27:e}; b27={b27:e}; D3Q19 control defects={defects19_by_n:?}; a19={a19:e}; b19={b19:e}; band_abs={band:e}"
    );
    assert!(
        a27.abs() <= band,
        "ACC CUM E3 {label}: h^2-extrapolated D3Q27 intrinsic offset |a27|={:.6e} exceeds band {band:.6e}; a27={a27:.6e}, b27={b27:.6e}, D3Q27 defects_by_N={defects27_by_n:?}, D3Q19 control a19={a19:.6e}, b19={b19:.6e}, D3Q19 defects_by_N={defects19_by_n:?}, normalization=resolution-independent intercept of nu_eff/nu - 1",
        a27.abs()
    );
    a27
}

fn orientation_spread(n: usize, planes: &[TgvPlane], band: f64, label: &str) {
    // Orientation argument. A scalar relaxation correction can only cancel the
    // isotropic part of the viscosity defect. If the residual depends on the
    // vortex plane, the same nominal nu produces different decay rates for
    // rotated copies of the same analytic flow. Therefore max pairwise
    // |nu_eff_i - nu_eff_j|/nu directly measures the non-scalar anisotropic
    // residue that a scalar offset cannot represent.
    let u0 = 1.28e-4 / n as f64;
    let values: Vec<f64> = planes
        .iter()
        .map(|&plane| measure_nu_eff::<D3Q19>(n, NU, u0, cumulant(), plane))
        .collect();
    let mut spread = 0.0f64;
    for i in 0..values.len() {
        for j in i + 1..values.len() {
            spread = spread.max((values[i] - values[j]).abs() / NU);
        }
    }
    println!(
        "ACC CUM E4 {label}: N={n} u0={u0:e} nu_eff_by_orientation={values:?} spread_over_nu={spread:e}"
    );
    assert!(
        spread <= band,
        "ACC CUM E4 {label}: orientation spread max |nu_eff_i - nu_eff_j|/nu={spread:.6e} exceeds band {band:.6e}; values={values:?}, normalization=nominal nu"
    );
}

#[test]
fn e4_d3q19_scalar_correction_orientation_canary_light() {
    orientation_spread(32, &[TgvPlane::Xy, TgvPlane::Yz], 2.0e-3, "light-canary");
}

#[test]
#[ignore = "heavy orientation audit: all three D3Q19 vortex planes at N=48"]
fn e4_d3q19_scalar_correction_orientation_invariance_heavy() {
    orientation_spread(
        48,
        &[TgvPlane::Xy, TgvPlane::Yz, TgvPlane::Zx],
        1.0e-3,
        "heavy",
    );
}
