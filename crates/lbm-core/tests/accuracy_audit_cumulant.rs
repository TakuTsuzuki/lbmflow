//! Accuracy audit for the cumulant/central-moment viscosity-defect correction.
//!
//! These tests audit the closure property, not the implementation constant:
//! a valid correction must remove the finite-frame lattice viscosity defect
//! across velocity amplitude, resolution, and orientation. Failing assertions
//! are triage findings, not invitations to loosen bands.

mod common;

use common::metrics::{curve_agreement, linear_fit, monotonicity};
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

fn cumulant() -> CollisionKind {
    CollisionKind::Cumulant {
        omega_shear: omega_from_nu(NU),
    }
}

fn trt() -> CollisionKind {
    CollisionKind::Trt { magic: TRT_MAGIC }
}

fn fit_defect_vs_u2<L: Lattice>(
    n: usize,
    u0s: &[f64],
    collision: CollisionKind,
) -> (f64, f64, f64, Vec<f64>) {
    let x: Vec<f64> = u0s.iter().map(|u0| u0 * u0).collect();
    let y: Vec<f64> = u0s
        .iter()
        .map(|&u0| measure_nu_eff::<L>(n, NU, u0, collision, TgvPlane::Xy) / NU - 1.0)
        .collect();
    let fit = linear_fit(&x, &y);
    // With y = c*u0^2 + b, the residual scatter around the least-squares line
    // estimates the noise floor of this probe. The TRT guard requires
    // |c_trt| >= 5*sigma_fit so the audit can actually see the cubic-frame
    // viscosity defect before crediting the cumulant correction.
    let dof = (x.len() as f64 - 2.0).max(1.0);
    let sigma_fit = x
        .iter()
        .zip(&y)
        .map(|(&xi, &yi)| {
            let r = yi - (fit.slope * xi + fit.intercept);
            r * r
        })
        .sum::<f64>()
        .sqrt()
        / dof.sqrt();
    (fit.slope, fit.intercept, sigma_fit, y)
}

fn assert_e1_amplitude_independence(n: usize, u0s: &[f64], label: &str) {
    // Defect relation. A finite-frame cubic lattice error enters the measured
    // decay as nu_eff = nu * (1 + c_lat*u0^2 + O(u0^4)) on a TGV because the
    // leading non-Galilean stress term is cubic in velocity and therefore
    // changes the viscous decay rate in proportion to u0^2. A correct closure
    // drives the residual coefficient toward zero across amplitudes, not at
    // one tuned point.
    let (c_trt, b_trt, sigma_fit, trt_y) = fit_defect_vs_u2::<D3Q19>(n, u0s, trt());
    let (c_cum, b_cum, _, cum_y) = fit_defect_vs_u2::<D3Q19>(n, u0s, cumulant());
    println!(
        "ACC CUM E1 {label}: N={n} u0={u0s:?} TRT rel={trt_y:?} c_trt={c_trt:e} b_trt={b_trt:e} sigma_fit={sigma_fit:e}; Cumulant rel={cum_y:?} c_res={c_cum:e} b_res={b_cum:e}"
    );
    assert!(
        c_trt.abs() >= 5.0 * sigma_fit.abs(),
        "ACC CUM E1 {label}: TRT sensitivity guard failed: |c_trt|={:.6e}, band >= 5*sigma_fit={:.6e}, sigma_fit={:.6e}, normalization=(nu_eff/nu-1)/u0^2",
        c_trt.abs(),
        5.0 * sigma_fit.abs(),
        sigma_fit
    );
    assert!(
        c_cum.abs() <= 0.2 * c_trt.abs(),
        "ACC CUM E1 {label}: cumulant residual defect |c_res|={:.6e} exceeds 20% of TRT control band {:.6e}; c_trt={:.6e}, normalization=(nu_eff/nu-1)/u0^2",
        c_cum.abs(),
        0.2 * c_trt.abs(),
        c_trt
    );
}

#[test]
fn e1_cumulant_nu_eff_independent_of_velocity_amplitude_light() {
    assert_e1_amplitude_independence(32, &[0.02, 0.04, 0.08], "light");
}

#[test]
#[ignore = "heavy asymptotic audit: D3Q19 cumulant amplitude sweep at N=48"]
fn e1_cumulant_nu_eff_independent_of_velocity_amplitude_heavy() {
    assert_e1_amplitude_independence(48, &[0.015, 0.03, 0.06], "heavy");
}

fn resolution_defects(ns: &[usize], band: f64, label: &str) {
    let samples: Vec<(f64, f64)> = ns
        .iter()
        .map(|&n| {
            let u0 = 1.28e-4 / n as f64;
            let nu_eff = measure_nu_eff::<D3Q19>(n, NU, u0, cumulant(), TgvPlane::Xy);
            (n as f64, nu_eff / NU - 1.0)
        })
        .collect();
    let defects: Vec<f64> = samples.iter().map(|&(_, d)| d).collect();
    let abs_defects: Vec<f64> = defects.iter().map(|d| d.abs()).collect();
    let max_abs = abs_defects.iter().copied().fold(0.0, f64::max);
    let neg_abs_defects: Vec<f64> = abs_defects.iter().map(|d| -d).collect();
    let growth = monotonicity(&neg_abs_defects);
    let zero_curve = curve_agreement(|_| 0.0, &samples, band, 1.0);
    println!(
        "ACC CUM E2 {label}: defects nu_eff/nu-1 by N = {samples:?}; max_abs={max_abs:e}; increasing_abs_fraction={growth:e}; zero_curve_max_dev={:e}",
        zero_curve.max_rel_dev
    );
    assert!(
        max_abs <= band,
        "ACC CUM E2 {label}: max |nu_eff/nu - 1|={max_abs:.6e} exceeds band {band:.6e}; values by N={samples:?}, normalization=relative viscosity defect"
    );
    assert!(
        growth < 1.0,
        "ACC CUM E2 {label}: |defect| grows monotonically with N, indicating a resolution trend in the resolution-independent offset; increasing_abs_fraction={growth:.6e}, values by N={samples:?}"
    );
}

#[test]
fn e2_d3q19_offset_resolution_canary_light() {
    resolution_defects(&[24, 32], 4.0e-3, "light-canary");
}

#[test]
#[ignore = "heavy resolution audit: D3Q19 cumulant at N=24,32,48"]
fn e2_d3q19_offset_resolution_validity_heavy() {
    resolution_defects(&[24, 32, 48], 2.0e-3, "heavy");
}

#[test]
fn e3_d3q27_needs_no_offset_light() {
    // Cross-lattice control. D3Q27 has the full tensor-product velocity set,
    // so the D3Q19 reduced-set bias story predicts no comparable small-u0
    // viscosity offset when the scalar correction offset is absent.
    let n = 32usize;
    let u0 = 1.28e-4 / n as f64;
    let nu_eff = measure_nu_eff::<D3Q27>(n, NU, u0, cumulant(), TgvPlane::Xy);
    let defect = nu_eff / NU - 1.0;
    println!("ACC CUM E3: D3Q27 N={n} u0={u0:e} nu_eff={nu_eff:e} defect={defect:e}");
    assert!(
        defect.abs() <= 2.0e-3,
        "ACC CUM E3: D3Q27 |nu_eff/nu - 1|={:.6e} exceeds band 2.000000e-3; nu_eff={nu_eff:.6e}, nu={NU:.6e}, normalization=relative viscosity defect",
        defect.abs()
    );
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
