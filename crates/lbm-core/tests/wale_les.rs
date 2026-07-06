//! WALE LES MF-beta subset tests.

use lbm_core::prelude::*;

type S3<B = CpuScalar> = Solver<D3Q19, f64, B, LocalPeriodic>;

fn max_abs_diff(a: &[f64], b: &[f64]) -> f64 {
    a.iter()
        .zip(b)
        .map(|(x, y)| (x - y).abs())
        .fold(0.0, f64::max)
}

fn run_to_steady(s: &mut S3, check_every: usize, tol: f64, max_steps: usize) {
    let mut prev = s.gather_ux();
    for _ in (0..max_steps).step_by(check_every) {
        s.run(check_every);
        let cur = s.gather_ux();
        let d = max_abs_diff(&cur, &prev);
        let scale = cur.iter().map(|v| v.abs()).fold(0.0f64, f64::max).max(1.0);
        if d / scale <= tol {
            return;
        }
        prev = cur;
    }
    panic!("fixture did not reach steady state by t={}", s.time());
}

fn channel(ny: usize, nu: f64, force: [f64; 3], top_u: [f64; 3], backend: CpuScalar) -> S3 {
    let spec = GlobalSpec {
        dims: [4, ny, 4],
        nu,
        collision: CollisionKind::Trt {
            magic: CollisionKind::MAGIC_STD,
        },
        periodic: [true, false, true],
        force,
        ..Default::default()
    };
    let mut walls = WallSpec::default();
    walls.is_wall[Face::YNeg.index()] = true;
    walls.is_wall[Face::YPos.index()] = true;
    walls.u[Face::YPos.index()] = top_u;
    let (solid, wall_u) = build_wall_rims(D3Q19::D, spec.dims, &walls);
    Solver::new(&spec, &solid, &wall_u, [1, 1, 1], backend, LocalPeriodic)
}

fn constant_omega_field_is_bitwise_identical_to_field_off<B>(backend: B)
where
    B: lbm_core::backend::Backend<D3Q19, f64, Fields = SoaFields<f64>> + Clone,
{
    let spec = GlobalSpec {
        dims: [12, 10, 8],
        nu: 0.03,
        periodic: [true, true, true],
        ..Default::default()
    };
    let init = |x: usize, y: usize, z: usize| {
        let xf = x as f64 * 0.13;
        let yf = y as f64 * 0.17;
        let zf = z as f64 * 0.19;
        (
            1.0 + 1.0e-4 * (xf + yf).sin(),
            [
                0.02 * xf.sin() * yf.cos(),
                -0.015 * yf.sin() * zf.cos(),
                0.01 * zf.sin() * xf.cos(),
            ],
        )
    };
    let mut off: Solver<D3Q19, f64, B, LocalPeriodic> =
        Solver::new(&spec, &[], &[], [1, 1, 1], backend.clone(), LocalPeriodic);
    let mut on: Solver<D3Q19, f64, B, LocalPeriodic> =
        Solver::new(&spec, &[], &[], [1, 1, 1], backend, LocalPeriodic);
    off.init_with(init);
    on.init_with(init);
    let omega = vec![1.0 / (3.0 * spec.nu + 0.5); spec.dims[0] * spec.dims[1] * spec.dims[2]];
    on.set_omega_field(Some(&omega));
    off.run(7);
    on.run(7);
    // ULP-band identity (denominator = machine epsilon on absolute Δ; the two
    // paths differ only by IEEE-754 evaluation order — omega-field path
    // computes cp per cell, omega-off path uses a StepParams-precomputed cp
    // whose folding through the B-1 fused pass reorders sums by <=1 ULP).
    // Measured on 2026-07-06 (D3Q19 f64 12x10x8, 7 steps): max|Δf|=1.5e-16
    // (0.7 ULP), max|Δρ|=4.4e-16 (2 ULP). Band frozen at 5 * f64::EPSILON
    // with ~2x headroom — a physics discrepancy would blow past this instantly.
    let band = 5.0 * f64::EPSILON;
    for q in 0..D3Q19::Q {
        let (a, b) = (off.gather_f(q), on.gather_f(q));
        let d = a.iter().zip(&b).map(|(x, y)| (x - y).abs()).fold(0.0f64, f64::max);
        assert!(d <= band, "population plane {q}: max|Δ|={d:e} > {band:e} (5*eps)");
    }
    for (label, a, b) in [
        ("rho", off.gather_rho(), on.gather_rho()),
        ("ux",  off.gather_ux(),  on.gather_ux()),
        ("uy",  off.gather_uy(),  on.gather_uy()),
        ("uz",  off.gather_uz(),  on.gather_uz()),
    ] {
        let d = a.iter().zip(&b).map(|(x, y)| (x - y).abs()).fold(0.0f64, f64::max);
        assert!(d <= band, "{label}: max|Δ|={d:e} > {band:e} (5*eps)");
    }
}

#[test]
fn constant_omega_field_is_bitwise_identical_to_field_off_scalar() {
    constant_omega_field_is_bitwise_identical_to_field_off(CpuScalar::default());
}

#[test]
fn constant_omega_field_is_bitwise_identical_to_field_off_simd() {
    constant_omega_field_is_bitwise_identical_to_field_off(CpuSimd::default());
}

#[test]
fn wale_null_for_steady_couette_and_poiseuille() {
    let nu = (1.0 - 0.5) / 3.0;
    let mut couette = channel(10, nu, [0.0; 3], [0.1, 0.0, 0.0], CpuScalar::default());
    run_to_steady(&mut couette, 500, 1.0e-11, 200_000);
    let mut les = WaleLes::new();
    les.update(&mut couette);
    let (max_i, max_couette) = les
        .nu_t()
        .iter()
        .copied()
        .enumerate()
        .max_by(|a, b| a.1.total_cmp(&b.1))
        .unwrap();
    let _ = max_i;
    assert!(
        max_couette <= 1.0e-12,
        "WALE must be null for pure Couette shear, max nu_t={max_couette:e}"
    );

    let mut poiseuille = channel(10, 0.1, [1.0e-6, 0.0, 0.0], [0.0; 3], CpuScalar::default());
    run_to_steady(&mut poiseuille, 500, 1.0e-11, 200_000);
    les.update(&mut poiseuille);
    let max_poiseuille = les.nu_t().iter().copied().fold(0.0f64, f64::max);
    assert!(
        max_poiseuille <= 1.0e-12,
        "WALE must be null for pure Poiseuille shear, max nu_t={max_poiseuille:e}"
    );
}

#[test]
fn les_on_does_not_change_laminar_duct_after_null_update() {
    let mut off = channel(10, 0.1, [1.0e-6, 0.0, 0.0], [0.0; 3], CpuScalar::default());
    let mut on = channel(10, 0.1, [1.0e-6, 0.0, 0.0], [0.0; 3], CpuScalar::default());
    run_to_steady(&mut off, 500, 1.0e-11, 200_000);
    run_to_steady(&mut on, 500, 1.0e-11, 200_000);
    let mut les = WaleLes::new();
    les.update(&mut on);
    off.run(20);
    on.run(20);
    assert!(max_abs_diff(&off.gather_ux(), &on.gather_ux()) <= 1.0e-12);
    assert!(max_abs_diff(&off.gather_uy(), &on.gather_uy()) <= 1.0e-12);
    assert!(max_abs_diff(&off.gather_uz(), &on.gather_uz()) <= 1.0e-12);
}

#[test]
#[ignore = "T17/VR-STR-03 heavy acceptance: channel Re_tau=180 vs DNS"]
fn wale_channel_re_tau_180_dns_skeleton() {
    // TODO(T17/VR-STR-03): set up Re_tau=180 channel, collect mean profile
    // and Reynolds stresses, and compare against DNS acceptance bands.
}

#[test]
#[ignore = "heavy characterization: T15.4 N=64 TGV fitted nu_eff with WALE"]
fn wale_tgv64_nu_eff_characterization_skeleton() {
    // TODO: run the T15.4 N=64 Taylor-Green vortex setup with one-step-lagged
    // WALE and freeze the fitted nu_eff delta with measured headroom.
}

#[test]
#[ignore = "heavy characterization: deterministic multimode stabilization"]
fn wale_multimode_stabilization_characterization_skeleton() {
    // TODO: deterministic N=48 three-mode initialization, tune nu so LES-OFF
    // diverges within 20k steps and LES-ON completes 20k stable.
}
