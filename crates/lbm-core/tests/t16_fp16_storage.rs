#![cfg(feature = "gpu")]
//! T16: FP16 distribution-storage characterization.
//!
//! Compute remains f32. These ignored tests compare the same GPU scenario with
//! f32 vs f16 distribution buffers and print the measured degradation values
//! PM will freeze into final bands.

use lbm_core::gpu::{GpuInitError, GpuStorage, KernelCfg};
use lbm_core::prelude::*;
use std::f64::consts::PI;
use std::sync::Arc;

type Gpu2 = Solver<D2Q9, f32, WgpuBackend<D2Q9>, LocalPeriodic>;

const T16_BAND_PLACEHOLDER: f64 = 1.0e-2; // BAND-FREEZE-PENDING(PM)

fn shader_f16_ctx_or_skip() -> Option<Arc<GpuContext>> {
    match GpuContext::new_with_shader_f16(true) {
        Ok(ctx) => Some(ctx),
        Err(GpuInitError::NoAdapter) => {
            eprintln!("BENCH-PENDING: skipping T16 FP16 characterization; no usable GPU adapter");
            None
        }
        Err(e) => panic!("T16 FP16 characterization requires SHADER_F16: {e}"),
    }
}

fn gpu_solver(
    ctx: &Arc<GpuContext>,
    spec: &GlobalSpec<f32>,
    solid: &[bool],
    wall_u: &[[f32; 3]],
    storage: GpuStorage,
) -> Gpu2 {
    let backend = WgpuBackend::<D2Q9>::with_config(ctx.clone(), KernelCfg { storage });
    Solver::new(spec, solid, wall_u, [1, 1, 1], backend, LocalPeriodic)
}

fn l2_rel(a: &[f64], b: &[f64]) -> f64 {
    assert_eq!(a.len(), b.len());
    let mut num = 0.0f64;
    let mut den = 0.0f64;
    for (&x, &y) in a.iter().zip(b) {
        let d = x - y;
        num += d * d;
        den += y * y;
    }
    num.sqrt() / den.sqrt().max(1.0e-30)
}

fn gather_u2(s: &Gpu2) -> Vec<f64> {
    let ux = s.gather_ux();
    let uy = s.gather_uy();
    let mut out = Vec::with_capacity(2 * ux.len());
    for i in 0..ux.len() {
        out.push(ux[i] as f64);
        out.push(uy[i] as f64);
    }
    out
}

#[test]
#[ignore = "heavy: requires SHADER_F16 GPU adapter and runs 256^2 TGV to one decay time"]
fn t16_tgv2d_f16_storage_degradation_vs_f32_gpu() {
    let Some(ctx) = shader_f16_ctx_or_skip() else {
        return;
    };
    let n = 256usize;
    let nu = 0.02f64;
    let u0 = 0.05f64;
    let k = 2.0 * PI / n as f64;
    let steps = (1.0 / (2.0 * nu * k * k)).round() as usize;
    let spec = GlobalSpec::<f32> {
        dims: [n, n, 1],
        nu,
        periodic: [true, true, false],
        ..Default::default()
    };
    let init = |x: usize, y: usize, _z: usize| {
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
    };

    let mut f32s = gpu_solver(&ctx, &spec, &[], &[], GpuStorage::F32);
    let mut f16s = gpu_solver(&ctx, &spec, &[], &[], GpuStorage::F16);
    f32s.init_with(init);
    f16s.init_with(init);
    f32s.run(steps);
    f16s.run(steps);

    let f32_u = gather_u2(&f32s);
    let f16_u = gather_u2(&f16s);
    let decay = (-2.0 * nu * k * k * steps as f64).exp();
    let mut analytic = Vec::with_capacity(2 * n * n);
    for y in 0..n {
        for x in 0..n {
            let (xf, yf) = (k * x as f64, k * y as f64);
            analytic.push(-u0 * xf.cos() * yf.sin() * decay);
            analytic.push(u0 * xf.sin() * yf.cos() * decay);
        }
    }
    let f16_vs_f32 = l2_rel(&f16_u, &f32_u);
    let f16_vs_analytic = l2_rel(&f16_u, &analytic);
    println!(
        "T16 TGV2D f16 storage: steps={steps}, f16_vs_f32_u_l2rel={f16_vs_f32:.9e}, \
         f16_vs_analytic_u_l2rel={f16_vs_analytic:.9e}"
    );
    assert!(
        f16_vs_f32 <= T16_BAND_PLACEHOLDER,
        "T16 TGV2D f16-vs-f32 u-field L2rel measured {f16_vs_f32:.9e} > placeholder band {T16_BAND_PLACEHOLDER:.9e}"
    );
    assert!(
        f16_vs_analytic <= T16_BAND_PLACEHOLDER,
        "T16 TGV2D f16-vs-analytic u-field L2rel measured {f16_vs_analytic:.9e} > placeholder band {T16_BAND_PLACEHOLDER:.9e}"
    );
}

fn cavity_centerlines(s: &Gpu2, n: usize, u_lid: f64) -> Vec<f64> {
    let ux = s.gather_ux();
    let uy = s.gather_uy();
    let c = n / 2;
    let mut out = Vec::with_capacity(2 * (n - 2));
    for y in 1..n - 1 {
        out.push(ux[y * n + c] as f64 / u_lid);
    }
    for x in 1..n - 1 {
        out.push(uy[c * n + x] as f64 / u_lid);
    }
    out
}

#[test]
#[ignore = "heavy: requires SHADER_F16 GPU adapter and runs 128^2 lid cavity for 40k steps"]
fn t16_cavity2d_f16_storage_degradation_vs_f32_gpu() {
    let Some(ctx) = shader_f16_ctx_or_skip() else {
        return;
    };
    let n = 128usize;
    let re = 100.0f64;
    let u_lid = 0.1f64;
    let nu = u_lid * (n - 2) as f64 / re;
    let steps = 40_000usize;
    let mut walls = WallSpec::default();
    for face in [Face::XNeg, Face::XPos, Face::YNeg, Face::YPos] {
        walls.is_wall[face.index()] = true;
    }
    walls.u[Face::YPos.index()] = [u_lid as f32, 0.0, 0.0];
    let spec = GlobalSpec::<f32> {
        dims: [n, n, 1],
        nu,
        periodic: [false, false, false],
        ..Default::default()
    };
    let (solid, wall_u) = build_wall_rims(2, spec.dims, &walls);
    let mut f32s = gpu_solver(&ctx, &spec, &solid, &wall_u, GpuStorage::F32);
    let mut f16s = gpu_solver(&ctx, &spec, &solid, &wall_u, GpuStorage::F16);
    f32s.run(steps);
    f16s.run(steps);

    let f32_profiles = cavity_centerlines(&f32s, n, u_lid);
    let f16_profiles = cavity_centerlines(&f16s, n, u_lid);
    let f16_vs_f32 = l2_rel(&f16_profiles, &f32_profiles);
    println!(
        "T16 cavity2D f16 storage: steps={steps}, Re={re:.1}, f16_vs_f32_centerline_l2rel={f16_vs_f32:.9e}"
    );
    assert!(
        f16_vs_f32 <= T16_BAND_PLACEHOLDER,
        "T16 cavity2D f16-vs-f32 centerline L2rel measured {f16_vs_f32:.9e} > placeholder band {T16_BAND_PLACEHOLDER:.9e}"
    );
}
