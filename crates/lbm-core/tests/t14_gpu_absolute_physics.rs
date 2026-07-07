#![cfg(feature = "gpu")]
//! T14 absolute GPU physics sentinels.
//!
//! These tests do not compare the GPU to the CPU backend. They compare GPU
//! fields directly with analytic low-Mach references so a shared CPU/GPU
//! implementation error cannot pass as "equivalent".

use lbm_core::prelude::*;
use std::f64::consts::PI;
use std::sync::{Arc, OnceLock};

fn gpu_ctx_or_skip() -> Option<Arc<GpuContext>> {
    static CTX: OnceLock<Result<Arc<GpuContext>, String>> = OnceLock::new();
    match CTX.get_or_init(|| GpuContext::new().map_err(|e| e.to_string())) {
        Ok(ctx) => Some(ctx.clone()),
        Err(e) => {
            if std::env::var_os("LBM_REQUIRE_GPU").is_some() {
                panic!("T14 absolute GPU tests require an adapter: {e}");
            }
            eprintln!("skipping T14 absolute GPU test: no adapter ({e})");
            None
        }
    }
}

fn l2_rel(actual: &[f64], reference: &[f64]) -> f64 {
    assert_eq!(actual.len(), reference.len());
    let mut num = 0.0;
    let mut den = 0.0;
    for (&a, &r) in actual.iter().zip(reference) {
        num += (a - r) * (a - r);
        den += r * r;
    }
    (num / den).sqrt()
}

fn tgv_gpu_l2(n: usize) -> Option<f64> {
    let ctx = gpu_ctx_or_skip()?;
    let nu = 0.02f64;
    let u0 = 1.28f64 / n as f64;
    let k = 2.0 * PI / n as f64;
    let spec = GlobalSpec::<f32> {
        dims: [n, n, 1],
        nu,
        periodic: [true, true, false],
        collision: CollisionKind::Trt { magic: 3.0 / 16.0 },
        ..Default::default()
    };
    let mut gpu = GpuSolver::<D2Q9>::new(&spec, &[], &[], ctx);
    gpu.init_with(|x, y, _| {
        let xf = k * x as f64;
        let yf = k * y as f64;
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
    let steps = (1.0 / (2.0 * nu * k * k)).round() as usize;
    gpu.run(steps);
    let decay = (-2.0 * nu * k * k * steps as f64).exp();
    let ux = gpu.gather_ux();
    let uy = gpu.gather_uy();
    let mut actual = Vec::with_capacity(2 * n * n);
    let mut reference = Vec::with_capacity(2 * n * n);
    for y in 0..n {
        for x in 0..n {
            let xf = k * x as f64;
            let yf = k * y as f64;
            let i = y * n + x;
            actual.push(ux[i] as f64);
            actual.push(uy[i] as f64);
            reference.push(-u0 * xf.cos() * yf.sin() * decay);
            reference.push(u0 * xf.sin() * yf.cos() * decay);
        }
    }
    Some(l2_rel(&actual, &reference))
}

#[test]
fn t14_gpu_tgv_matches_analytic_decay_and_converges() {
    let Some(e32) = tgv_gpu_l2(32) else {
        return;
    };
    let e64 = tgv_gpu_l2(64).expect("same GPU context should be reusable");
    let order = (e32 / e64).log2();
    println!("T14 GPU absolute TGV: e32={e32:.6e}, e64={e64:.6e}, order={order:.3}");
    assert!(
        e64 <= 2.0e-3,
        "T14 GPU absolute TGV N=64 L2rel={e64:e} > 2e-3 analytic band"
    );
    assert!(
        order >= 1.6,
        "T14 GPU absolute TGV convergence order={order:e} < 1.6, e32={e32:e}, e64={e64:e}"
    );
}

#[test]
fn t14_gpu_pressure_channel_matches_poiseuille_bulk_shape() {
    let Some(ctx) = gpu_ctx_or_skip() else {
        return;
    };
    let (nx, ny) = (96usize, 34usize);
    let h = (ny - 2) as f64;
    let nu = 0.04f64;
    let u_peak = 0.02f64;
    let cs2 = 1.0 / 3.0;
    let delta_rho = 8.0 * nu * u_peak * (nx - 1) as f64 / (cs2 * h * h);
    let rho_mid = 1.0f64;
    let mut walls = WallSpec::<f32>::default();
    walls.is_wall[Face::YNeg.index()] = true;
    walls.is_wall[Face::YPos.index()] = true;
    let mut faces = [FaceBC::Closed; 6];
    faces[Face::XNeg.index()] = FaceBC::Pressure {
        rho: (rho_mid + 0.5 * delta_rho) as f32,
    };
    faces[Face::XPos.index()] = FaceBC::Pressure {
        rho: (rho_mid - 0.5 * delta_rho) as f32,
    };
    let spec = GlobalSpec::<f32> {
        dims: [nx, ny, 1],
        nu,
        periodic: [false, false, false],
        faces,
        collision: CollisionKind::Trt { magic: 3.0 / 16.0 },
        ..Default::default()
    };
    let (solid, wall_u) = build_wall_rims(2, spec.dims, &walls);
    let mut gpu = GpuSolver::<D2Q9>::new(&spec, &solid, &wall_u, ctx);
    gpu.init_with(|x, y, _| {
        if y == 0 || y == ny - 1 {
            return (1.0, [0.0, 0.0, 0.0]);
        }
        let rho = rho_mid + 0.5 * delta_rho - delta_rho * x as f64 / (nx - 1) as f64;
        let yw = y as f64 - 0.5;
        let ux = 4.0 * u_peak * yw * (h - yw) / (h * h);
        (rho as f32, [ux as f32, 0.0, 0.0])
    });
    gpu.run(2_000);

    let ux = gpu.gather_ux();
    let x = nx / 2;
    let mut actual = Vec::with_capacity(ny - 2);
    let mut reference = Vec::with_capacity(ny - 2);
    for y in 1..ny - 1 {
        let yw = y as f64 - 0.5;
        actual.push(ux[y * nx + x] as f64);
        reference.push(4.0 * u_peak * yw * (h - yw) / (h * h));
    }
    let profile_l2 = l2_rel(&actual, &reference);
    let center_ux = actual[actual.len() / 2];
    let near_wall_ux = actual[0];
    println!(
        "T14 GPU absolute pressure channel: delta_rho={delta_rho:.6e}, \
         center profile L2rel={profile_l2:.6e}, center_ux={center_ux:.6e}, \
         near_wall_ux={near_wall_ux:.6e}"
    );
    assert!(
        profile_l2 <= 3.0e-2,
        "T14 GPU pressure-channel center profile L2rel={profile_l2:e} > 3e-2 analytic band"
    );
    assert!(
        center_ux > 0.0 && center_ux > near_wall_ux && near_wall_ux > 0.0,
        "T14 GPU pressure-channel behavior anchor failed: center_ux={center_ux:e}, near_wall_ux={near_wall_ux:e}"
    );
}
