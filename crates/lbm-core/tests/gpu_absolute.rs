#![cfg(feature = "gpu")]
#![allow(deprecated)]
//! GPU-direct absolute physics gates.
//!
//! These tests intentionally do not compare GPU fields to CPU fields. They
//! compare the f32 Wgpu backend directly against analytic solutions or the
//! Ghia et al. cavity reference, closing the CPU-equivalence blind spot called
//! out by VALIDATION.md T14 / SOLVER_IMPROVEMENT_SPEC D-3.

use lbm_core::gpu::GpuInitError;
use lbm_core::prelude::*;
use std::f64::consts::PI;
use std::sync::{Arc, OnceLock};

type Gpu2 = GpuSolver<D2Q9>;
type Gpu3 = GpuSolver<D3Q19>;

const TGV2D_N64_L2_BAND: f64 = 1.5e-3;
const TGV2D_ORDER_BAND: f64 = 1.7;
const TGV2D_NU_EFF_REL_BAND: f64 = 0.02;
const CAVITY_RE100_RMS_BAND: f64 = 0.02 * CAVITY_U;
const TGV3D_RATE_REL_BAND: f64 = 0.02;
const POISEUILLE_F32_LINF_REL_BAND: f64 = 1.0e-5;

const CAVITY_N: usize = 129;
const CAVITY_L: f64 = (CAVITY_N - 2) as f64;
const CAVITY_U: f64 = 0.1;

// Ghia, U., Ghia, K. N., & Shin, C. T. (1982), JCP 48(3), 387-411.
const GHIA_Y: [f64; 17] = [
    1.0000, 0.9766, 0.9688, 0.9609, 0.9531, 0.8516, 0.7344, 0.6172, 0.5000, 0.4531, 0.2813, 0.1719,
    0.1016, 0.0703, 0.0625, 0.0547, 0.0000,
];
const GHIA_X: [f64; 17] = [
    1.0000, 0.9688, 0.9609, 0.9531, 0.9453, 0.9063, 0.8594, 0.8047, 0.5000, 0.2344, 0.2266, 0.1563,
    0.0938, 0.0781, 0.0703, 0.0625, 0.0000,
];
const U_RE100: [f64; 17] = [
    1.00000, 0.84123, 0.78871, 0.73722, 0.68717, 0.23151, 0.00332, -0.13641, -0.20581, -0.21090,
    -0.15662, -0.10150, -0.06434, -0.04775, -0.04192, -0.03717, 0.00000,
];
const V_RE100: [f64; 17] = [
    0.00000, -0.05906, -0.07391, -0.08864, -0.10313, -0.16914, -0.22445, -0.24533, 0.05454,
    0.17527, 0.17507, 0.16077, 0.12317, 0.10890, 0.10091, 0.09233, 0.00000,
];

fn gpu_ctx_or_skip() -> Option<Arc<GpuContext>> {
    static CTX: OnceLock<Result<Arc<GpuContext>, GpuInitError>> = OnceLock::new();
    match CTX.get_or_init(GpuContext::new) {
        Ok(ctx) => Some(ctx.clone()),
        Err(e @ GpuInitError::NoAdapter) => {
            if std::env::var_os("LBM_REQUIRE_GPU").is_some() {
                panic!("GPU absolute physics tests require an adapter: {e}");
            }
            eprintln!("PENDING-NATIVE-RUN: skipping GPU absolute physics test; no adapter ({e})");
            None
        }
        Err(e) => panic!("GPU absolute physics tests require a usable adapter: {e}"),
    }
}

fn l2_rel(actual: &[f64], reference: &[f64]) -> f64 {
    assert_eq!(actual.len(), reference.len());
    let mut num = 0.0;
    let mut den = 0.0;
    for (&a, &r) in actual.iter().zip(reference) {
        let d = a - r;
        num += d * d;
        den += r * r;
    }
    (num / den.max(1.0e-30)).sqrt()
}

fn assert_strictly_decreasing(values: &[f64], what: &str) {
    assert!(values.len() >= 2);
    for (i, w) in values.windows(2).enumerate() {
        assert!(
            w[1] < w[0],
            "{what} must decrease monotonically, sample {i}: {} -> {}",
            w[0],
            w[1]
        );
    }
}

fn ke2d(s: &mut Gpu2) -> f64 {
    let ux = s.gather_ux();
    let uy = s.gather_uy();
    ux.iter()
        .zip(&uy)
        .map(|(&a, &b)| {
            let (a, b) = (a as f64, b as f64);
            a * a + b * b
        })
        .sum()
}

fn ke3d(s: &mut Gpu3) -> f64 {
    let ux = s.gather_ux();
    let uy = s.gather_uy();
    let uz = s.gather_uz();
    ux.iter()
        .zip(&uy)
        .zip(&uz)
        .map(|((&a, &b), &c)| {
            let (a, b, c) = (a as f64, b as f64, c as f64);
            a * a + b * b + c * c
        })
        .sum()
}

fn run_tgv2d_case(ctx: Arc<GpuContext>, n: usize) -> (f64, f64, Vec<f64>) {
    let nu = 0.02f64;
    let u0 = 1.28 / n as f64;
    let k = 2.0 * PI / n as f64;
    let steps = (1.0 / (2.0 * nu * k * k)).round() as usize;
    let spec = GlobalSpec::<f32> {
        dims: [n, n, 1],
        nu,
        periodic: [true, true, false],
        ..Default::default()
    };
    let mut s = Gpu2::new(&spec, &[], &[], ctx);
    s.init_with(move |x, y, _| {
        let (xf, yf) = (k * x as f64, k * y as f64);
        let rho = 1.0 - 3.0 * u0 * u0 / 4.0 * ((2.0 * xf).cos() + (2.0 * yf).cos());
        (
            rho as f32,
            [
                (-u0 * xf.cos() * yf.sin()) as f32,
                (u0 * xf.sin() * yf.cos()) as f32,
                0.0,
            ],
        )
    });

    let mut energy = vec![ke2d(&mut s)];
    let mut done = 0usize;
    for target in [steps / 4, steps / 2, 3 * steps / 4, steps] {
        s.run(target - done);
        done = target;
        energy.push(ke2d(&mut s));
    }
    assert_strictly_decreasing(&energy, "GPU D2Q9 TGV kinetic energy");

    let ux = s.gather_ux();
    let uy = s.gather_uy();
    let decay = (-2.0 * nu * k * k * steps as f64).exp();
    let mut actual = Vec::with_capacity(2 * n * n);
    let mut reference = Vec::with_capacity(2 * n * n);
    for y in 0..n {
        for x in 0..n {
            let (xf, yf) = (k * x as f64, k * y as f64);
            let i = y * n + x;
            actual.push(ux[i] as f64);
            actual.push(uy[i] as f64);
            reference.push(-u0 * xf.cos() * yf.sin() * decay);
            reference.push(u0 * xf.sin() * yf.cos() * decay);
        }
    }
    let err = l2_rel(&actual, &reference);
    let nu_eff = (energy[0] / *energy.last().unwrap()).ln() / (4.0 * k * k * steps as f64);
    (err, nu_eff, energy)
}

#[test]
fn gpu_d2q9_tgv_decay_matches_analytic_t1_f32() {
    let Some(ctx) = gpu_ctx_or_skip() else {
        return;
    };
    let (e32, _, _) = run_tgv2d_case(ctx.clone(), 32);
    let (e64, nu_eff, _) = run_tgv2d_case(ctx, 64);
    let order = (e32 / e64).log2();
    let nu_eff_rel = (nu_eff / 0.02 - 1.0).abs();

    // Bands cite validation_tgv.rs:
    // - t1_tgv_trt_accuracy_and_second_order_convergence: N=64 L2rel <= 1.5e-3, order >= 1.7.
    // - t1_tgv_effective_viscosity_within_two_percent: nu_eff within 2%.
    // The GPU path is f32-only; these are the same T1 band class, evaluated
    // directly against the analytic TGV field and decay law rather than CPU.
    assert!(
        e64 <= TGV2D_N64_L2_BAND,
        "GPU D2Q9 TGV N=64 L2rel {e64:.6e} > T1 band {TGV2D_N64_L2_BAND:.6e}"
    );
    assert!(
        order >= TGV2D_ORDER_BAND,
        "GPU D2Q9 TGV order {order:.6e} < T1 band {TGV2D_ORDER_BAND:.6e} (e32={e32:.6e}, e64={e64:.6e})"
    );
    assert!(
        nu_eff_rel <= TGV2D_NU_EFF_REL_BAND,
        "GPU D2Q9 TGV nu_eff relative error {nu_eff_rel:.6e} > T1 band {TGV2D_NU_EFF_REL_BAND:.6e}"
    );
}

fn cavity_re100(ctx: Arc<GpuContext>) -> Gpu2 {
    let re = 100.0f64;
    let nu = CAVITY_U * CAVITY_L / re;
    let mut walls = WallSpec::<f32>::default();
    for face in [Face::XNeg, Face::XPos, Face::YNeg, Face::YPos] {
        walls.is_wall[face.index()] = true;
    }
    walls.u[Face::YPos.index()] = [CAVITY_U as f32, 0.0, 0.0];
    let spec = GlobalSpec::<f32> {
        dims: [CAVITY_N, CAVITY_N, 1],
        nu,
        periodic: [false, false, false],
        ..Default::default()
    };
    let (solid, wall_u) = build_wall_rims(2, spec.dims, &walls);
    Gpu2::new(&spec, &solid, &wall_u, ctx)
}

fn run_to_steady_gpu2(s: &mut Gpu2, check_every: usize, tol: f64, max_steps: usize) -> bool {
    let mut prev = Vec::<f32>::new();
    let mut elapsed = 0usize;
    while elapsed < max_steps {
        s.run(check_every);
        elapsed += check_every;
        let ux = s.gather_ux();
        let uy = s.gather_uy();
        let mut dmax = 0.0f64;
        let mut umax = 0.0f64;
        if !prev.is_empty() {
            for (i, &x) in ux.iter().enumerate() {
                let y = uy[i];
                dmax = dmax
                    .max((x as f64 - prev[2 * i] as f64).abs())
                    .max((y as f64 - prev[2 * i + 1] as f64).abs());
                umax = umax.max((x as f64).abs()).max((y as f64).abs());
            }
            if umax > 0.0 && dmax <= tol * umax {
                return true;
            }
        }
        prev.clear();
        prev.reserve(2 * ux.len());
        for i in 0..ux.len() {
            prev.push(ux[i]);
            prev.push(uy[i]);
        }
    }
    false
}

fn sample_centerline_u(ux: &[f32], y_frac: f64) -> f64 {
    let x = CAVITY_N / 2;
    let pos = 0.5 + y_frac * CAVITY_L;
    let y0 = pos.floor().clamp(1.0, (CAVITY_N - 2) as f64) as usize;
    let y1 = (y0 + 1).min(CAVITY_N - 2);
    let t = pos - y0 as f64;
    (1.0 - t) * ux[y0 * CAVITY_N + x] as f64 + t * ux[y1 * CAVITY_N + x] as f64
}

fn sample_centerline_v(uy: &[f32], x_frac: f64) -> f64 {
    let y = CAVITY_N / 2;
    let pos = 0.5 + x_frac * CAVITY_L;
    let x0 = pos.floor().clamp(1.0, (CAVITY_N - 2) as f64) as usize;
    let x1 = (x0 + 1).min(CAVITY_N - 2);
    let t = pos - x0 as f64;
    (1.0 - t) * uy[y * CAVITY_N + x0] as f64 + t * uy[y * CAVITY_N + x1] as f64
}

fn cavity_rms_centerline_error(ux: &[f32], uy: &[f32]) -> f64 {
    let mut sum = 0.0;
    let mut n = 0usize;
    for i in 0..17 {
        let du = sample_centerline_u(ux, GHIA_Y[i]) - CAVITY_U * U_RE100[i];
        let dv = sample_centerline_v(uy, GHIA_X[i]) - CAVITY_U * V_RE100[i];
        sum += du * du + dv * dv;
        n += 2;
    }
    (sum / n as f64).sqrt()
}

fn primary_vortex_center_and_psi(ux: &[f32]) -> (f64, f64, f64) {
    let mut best = (0usize, 0usize, 0.0f64);
    for x in 1..=(CAVITY_N - 2) {
        let mut psi = 0.0f64;
        for y in 1..=(CAVITY_N - 2) {
            psi += ux[y * CAVITY_N + x] as f64;
            if psi.abs() > best.2.abs() {
                best = (x, y, psi);
            }
        }
    }
    (
        (best.0 as f64 - 0.5) / CAVITY_L,
        (best.1 as f64 - 0.5) / CAVITY_L,
        best.2,
    )
}

#[test]
#[ignore = "heavy: requires a native GPU adapter and runs 129^2 Re=100 cavity to steady/Ghia"]
fn gpu_d2q9_cavity_re100_matches_ghia_t7_f32() {
    let Some(ctx) = gpu_ctx_or_skip() else {
        return;
    };
    let mut s = cavity_re100(ctx);
    let steady = run_to_steady_gpu2(&mut s, 1_000, 1.0e-8, 300_000);
    let ux = s.gather_ux();
    let uy = s.gather_uy();
    let rms = cavity_rms_centerline_error(&ux, &uy);
    let (cx, cy, psi) = primary_vortex_center_and_psi(&ux);

    // Band cites validation_cavity.rs::t7_lid_driven_cavity_re100_matches_ghia:
    // RMS <= 0.02*U and primary vortex center within +/-0.02L. VALIDATION.md
    // T14/D-3 explicitly names the same f32 absolute-physics line for GPU.
    assert!(
        rms <= CAVITY_RE100_RMS_BAND,
        "GPU D2Q9 cavity Re=100 Ghia RMS {rms:.6e} > T7 band {CAVITY_RE100_RMS_BAND:.6e}, steady={steady}, time={}",
        s.time()
    );
    assert!(
        (cx - 0.6172).abs() <= 0.02 && (cy - 0.7344).abs() <= 0.02,
        "GPU D2Q9 cavity Re=100 vortex center ({cx:.6e}, {cy:.6e}) outside T7 +/-0.02L band, psi={psi:.6e}, steady={steady}, time={}",
        s.time()
    );
    assert!(
        psi < 0.0 && sample_centerline_u(&ux, 0.5) < 0.0 && sample_centerline_v(&uy, 0.5) > 0.0,
        "GPU D2Q9 cavity Re=100 primary vortex sign anchor failed: psi={psi:.6e}, u_mid={:.6e}, v_mid={:.6e}",
        sample_centerline_u(&ux, 0.5),
        sample_centerline_v(&uy, 0.5)
    );
}

fn run_tgv3d_case(ctx: Arc<GpuContext>) -> (f64, Vec<f64>) {
    let n = 64usize;
    let nu = 0.02f32;
    let u0_coef = 1.28e-4f32;
    let u0 = u0_coef / n as f32;
    let k64 = 2.0 * PI / n as f64;
    let k = k64 as f32;
    let steps = (0.1 / (nu as f64 * k64 * k64)).round() as usize;
    let spec = GlobalSpec::<f32> {
        dims: [n, n, n],
        nu: nu as f64,
        periodic: [true, true, true],
        ..Default::default()
    };
    let mut s = Gpu3::new(&spec, &[], &[], ctx);
    s.init_with(move |x, y, z| {
        let (xf, yf, zf) = (k * x as f32, k * y as f32, k * z as f32);
        let u = [
            u0 * xf.sin() * yf.cos() * zf.cos(),
            -u0 * xf.cos() * yf.sin() * zf.cos(),
            0.0,
        ];
        let p = u0 * u0 / 16.0 * (((2.0 * xf).cos() + (2.0 * yf).cos()) * ((2.0 * zf).cos() + 2.0));
        (1.0 + 3.0 * p, u)
    });

    let mut energy = vec![ke3d(&mut s)];
    let mut done = 0usize;
    for target in [steps / 4, steps / 2, 3 * steps / 4, steps] {
        s.run(target - done);
        done = target;
        energy.push(ke3d(&mut s));
    }
    assert_strictly_decreasing(&energy, "GPU D3Q19 TGV3D kinetic energy");
    let rate = -(energy.last().unwrap() / energy[0]).ln() / steps as f64;
    (rate, energy)
}

#[test]
fn gpu_d3q19_tgv3d_decay_matches_analytic_t15_4_f32() {
    let Some(ctx) = gpu_ctx_or_skip() else {
        return;
    };
    let n = 64usize;
    let nu = 0.02f64;
    let k = 2.0 * PI / n as f64;
    let (rate, _) = run_tgv3d_case(ctx);
    let rate_ref = 6.0 * nu * k * k;
    let rate_rel = (rate - rate_ref).abs() / rate_ref;

    // Band cites t15_3d_f32.rs::t15_4_f32_tgv3d_decay_rate_matches_diffusion_limit:
    // N=64, nu=0.02, u0_coef=1.28e-4, energy decay-rate relative error <= 2%.
    assert!(
        rate_rel <= TGV3D_RATE_REL_BAND,
        "GPU D3Q19 TGV3D decay-rate rel {rate_rel:.6e} > T15.4 f32 band {TGV3D_RATE_REL_BAND:.6e} (rate={rate:.6e}, ref={rate_ref:.6e})"
    );
}

fn run_poiseuille_case(ctx: Arc<GpuContext>) -> (Vec<f64>, Vec<f64>) {
    let ny = 10usize;
    let nu = 0.1f64;
    let g = 1.0e-6f64;
    let mut walls = WallSpec::<f32>::default();
    walls.is_wall[Face::YNeg.index()] = true;
    walls.is_wall[Face::YPos.index()] = true;
    let spec = GlobalSpec::<f32> {
        dims: [4, ny, 1],
        nu,
        periodic: [true, false, false],
        force: [g as f32, 0.0, 0.0],
        ..Default::default()
    };
    let (solid, wall_u) = build_wall_rims(2, spec.dims, &walls);
    let mut s = Gpu2::new(&spec, &solid, &wall_u, ctx);
    let steady = run_to_steady_gpu2(&mut s, 500, 1.0e-7, 200_000);
    assert!(
        steady,
        "GPU D2Q9 Poiseuille did not reach f32 steady state, time={}",
        s.time()
    );
    let ux = s.gather_ux();
    let h = (ny - 2) as f64;
    let actual = (1..=(ny - 2)).map(|y| ux[y * 4] as f64).collect::<Vec<_>>();
    let reference = (1..=(ny - 2))
        .map(|y| {
            let yw = y as f64 - 0.5;
            g / (2.0 * nu) * yw * (h - yw)
        })
        .collect::<Vec<_>>();
    (actual, reference)
}

#[test]
fn gpu_d2q9_poiseuille_body_force_matches_analytic_t2_f32() {
    let Some(ctx) = gpu_ctx_or_skip() else {
        return;
    };
    let (actual, reference) = run_poiseuille_case(ctx);
    let den = reference.iter().copied().map(f64::abs).fold(0.0, f64::max);
    let err = actual
        .iter()
        .zip(&reference)
        .map(|(a, r)| (a - r).abs())
        .fold(0.0, f64::max)
        / den;

    // Band cites validation_channel.rs::t2_trt_magic_poiseuille_is_exact_and_symmetric:
    // the f64 TRT half-way wall gate is L_inf_rel <= 1e-10 plus top/bottom
    // symmetry. Wgpu is f32-only, so this uses the f32 field band class already
    // used for validation-grade f32 gates (T14 field and t15_3d_f32 conventions):
    // 1e-5 relative, without changing the analytic profile or adding a CPU oracle.
    assert!(
        err <= POISEUILLE_F32_LINF_REL_BAND,
        "GPU D2Q9 Poiseuille L_inf_rel {err:.6e} > f32 absolute band {POISEUILLE_F32_LINF_REL_BAND:.6e}, actual={actual:?}, reference={reference:?}"
    );

    let n = actual.len();
    for j in 0..n / 2 {
        let d = (actual[j] - actual[n - 1 - j]).abs() / den;
        assert!(
            d <= POISEUILLE_F32_LINF_REL_BAND,
            "GPU D2Q9 Poiseuille symmetry rel diff {d:.6e} > {POISEUILLE_F32_LINF_REL_BAND:.6e}, row={j}, actual={actual:?}"
        );
    }
    let center = n / 2;
    assert!(
        actual[center] > actual[0] && actual[center - 1] > actual[0],
        "GPU D2Q9 Poiseuille behavior anchor failed: center rows must exceed wall-adjacent rows, actual={actual:?}"
    );
}
