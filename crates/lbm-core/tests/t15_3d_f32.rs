//! T15 f32 x 3D validation for the product D3Q19 path.
//!
//! `lbm-scenario::Sim3Handle::F32` wraps this concrete solver type; the core
//! crate cannot import the downstream scenario crate without changing the
//! dependency graph, so these tests pin the wrapped engine directly.

use lbm_core::prelude::*;
use std::f64::consts::PI;

type S2 = Solver<D2Q9, f32, CpuScalar, LocalPeriodic>;
type S3 = Solver<D3Q19, f32, CpuScalar, LocalPeriodic>;

fn zslice(f: &[f32], nx: usize, ny: usize, z: usize) -> &[f32] {
    &f[z * nx * ny..(z + 1) * nx * ny]
}

fn linf_rel(a: &[f32], b: &[f32], floor: f64) -> f64 {
    assert_eq!(a.len(), b.len());
    let mut d = 0.0f64;
    let mut m = 0.0f64;
    for (&x, &y) in a.iter().zip(b) {
        d = d.max((x as f64 - y as f64).abs());
        m = m.max(y.abs() as f64);
    }
    d / m.max(floor)
}

fn linf_rel_scale(a: &[f32], b: &[f32], scale: f64) -> f64 {
    assert_eq!(a.len(), b.len());
    let d = a
        .iter()
        .zip(b)
        .map(|(&x, &y)| (x as f64 - y as f64).abs())
        .fold(0.0f64, f64::max);
    d / scale
}

fn mass_f64(s: &S3) -> f64 {
    let (fluid, deviation) = s.local_mass_partials();
    fluid + deviation
}

fn init_tgv3d(s: &mut S3, n: usize, nu: f32, u0_coef: f32) {
    let u0 = u0_coef / n as f32;
    let k = (2.0 * PI / n as f64) as f32;
    assert_eq!(s.dims(), [n, n, n]);
    assert!(nu > 0.0);
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
}

fn tgv3d_spec(n: usize, nu: f32) -> GlobalSpec<f32> {
    GlobalSpec::<f32> {
        dims: [n, n, n],
        nu: nu as f64,
        periodic: [true, true, true],
        ..Default::default()
    }
}

fn tgv3d_ke(s: &S3) -> f64 {
    let (ux, uy, uz) = (s.gather_ux(), s.gather_uy(), s.gather_uz());
    ux.iter()
        .zip(&uy)
        .zip(&uz)
        .map(|((&a, &b), &c)| {
            let (a, b, c) = (a as f64, b as f64, c as f64);
            a * a + b * b + c * c
        })
        .sum()
}

#[test]
fn t15_1_f32_tgv_z_invariant_degenerates_to_d2q9() {
    let n = 32;
    let nz = 4;
    let nu = 0.02f32;
    let u0 = 1.28 / n as f32;
    let k64 = 2.0 * PI / n as f64;
    let k = k64 as f32;
    let init = move |x: usize, y: usize| {
        let (xf, yf) = (k * x as f32, k * y as f32);
        let rho = 1.0 - 3.0 * u0 * u0 / 4.0 * ((2.0 * xf).cos() + (2.0 * yf).cos());
        (
            rho,
            [-u0 * xf.cos() * yf.sin(), u0 * xf.sin() * yf.cos(), 0.0],
        )
    };
    let steps = (1.0 / (2.0 * nu as f64 * k64 * k64)).round() as usize;

    let spec2 = GlobalSpec::<f32> {
        dims: [n, n, 1],
        nu: nu as f64,
        periodic: [true, true, false],
        ..Default::default()
    };
    let mut s2: S2 = Solver::new(
        &spec2,
        &[],
        &[],
        [1, 1, 1],
        CpuScalar::default(),
        LocalPeriodic,
    );
    s2.init_with(move |x, y, _| init(x, y));
    s2.run(steps);

    let spec3 = GlobalSpec::<f32> {
        dims: [n, n, nz],
        nu: nu as f64,
        periodic: [true, true, true],
        ..Default::default()
    };
    let mut s3: S3 = Solver::new(
        &spec3,
        &[],
        &[],
        [1, 1, 1],
        CpuScalar::default(),
        LocalPeriodic,
    );
    s3.init_with(move |x, y, _| init(x, y));
    s3.run(steps);

    let (r2, x2, y2) = (s2.gather_rho(), s2.gather_ux(), s2.gather_uy());
    let (r3, x3, y3, z3) = (
        s3.gather_rho(),
        s3.gather_ux(),
        s3.gather_uy(),
        s3.gather_uz(),
    );
    let mut rrho = 0.0f64;
    let mut rux = 0.0f64;
    let mut ruy = 0.0f64;
    let u_scale = u0 as f64;
    for z in 0..nz {
        rrho = rrho.max(linf_rel(zslice(&r3, n, n, z), &r2, 1.0));
        rux = rux.max(linf_rel_scale(zslice(&x3, n, n, z), &x2, u_scale));
        ruy = ruy.max(linf_rel_scale(zslice(&y3, n, n, z), &y2, u_scale));
    }
    let ruz = z3.iter().map(|&v| (v as f64).abs()).fold(0.0f64, f64::max) / u_scale;
    let rel = rrho.max(rux).max(ruy).max(ruz);
    println!(
        "T15-1 f32 z-invariant TGV over {steps} steps: rho rel {rrho:.3e}, ux rel {rux:.3e}, uy rel {ruy:.3e}, |uz|/u0 {ruz:.3e}"
    );
    assert!(
        rel <= 1.0e-5,
        "T15-1 f32 z-invariant D3Q19 vs D2Q9 rel {rel:.3e} > 1.0e-5 \
         (rho {rrho:.3e}, ux {rux:.3e}, uy {ruy:.3e}, |uz|/u0 {ruz:.3e})"
    );
}

#[test]
fn t15_4_f32_tgv3d_decay_rate_matches_diffusion_limit() {
    let n = 64;
    let nu = 0.02f32;
    let u0_coef = 1.28e-4f32;
    let k64 = 2.0 * PI / n as f64;
    let steps = (0.1 / (nu as f64 * k64 * k64)).round() as usize;
    let spec = tgv3d_spec(n, nu);
    let mut s: S3 = Solver::new(
        &spec,
        &[],
        &[],
        [1, 1, 1],
        CpuScalar::default(),
        LocalPeriodic,
    );
    init_tgv3d(&mut s, n, nu, u0_coef);
    let e0 = tgv3d_ke(&s);
    s.run(steps);
    let e1 = tgv3d_ke(&s);
    let rate = -(e1 / e0).ln() / steps as f64;
    let rate_ref = 6.0 * nu as f64 * k64 * k64;
    let rate_rel = (rate - rate_ref).abs() / rate_ref;
    println!(
        "T15-4 f32 TGV3D N={n}: steps={steps}, rate={rate:.6e}, ref={rate_ref:.6e}, rel={rate_rel:.3e}, u0={:.3e}",
        u0_coef / n as f32
    );
    assert!(
        rate_rel <= 0.02,
        "T15-4 f32 N={n} decay-rate relative error {rate_rel:.3e} > 2% \
         (rate {rate:.6e}, ref {rate_ref:.6e}, steps {steps})"
    );
}

#[test]
fn t15_4_f32_tgv3d_mass_drift_stays_below_1e_minus_5_per_1000_steps() {
    let n = 64;
    let nu = 0.02f32;
    let u0_coef = 1.28e-4f32;
    let steps = 1000;
    let spec = tgv3d_spec(n, nu);
    let mut s: S3 = Solver::new(
        &spec,
        &[],
        &[],
        [1, 1, 1],
        CpuScalar::default(),
        LocalPeriodic,
    );
    init_tgv3d(&mut s, n, nu, u0_coef);
    let m0 = mass_f64(&s);
    s.run(steps);
    let m1 = mass_f64(&s);
    let rel_per_1000 = ((m1 - m0).abs() / m0.abs()) * (1000.0 / steps as f64);
    println!(
        "T15-4 f32 TGV3D mass drift N={n}: m0={m0:.9e}, m1={m1:.9e}, rel/1000={rel_per_1000:.3e}"
    );
    assert!(
        rel_per_1000 <= 1.0e-5,
        "T15-4 f32 TGV3D mass drift {rel_per_1000:.3e}/1000 steps > 1.0e-5 \
         (m0 {m0:.9e}, m1 {m1:.9e}, steps {steps})"
    );
}
