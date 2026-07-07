#![cfg(feature = "gpu")]
#![allow(deprecated)]
//! T14 D3Q19 absolute GPU physics sentinels.
//!
//! These tests compare the WGPU D3Q19 path directly with analytic references,
//! not with the CPU backend, so a shared implementation error cannot pass as
//! backend equivalence.

use lbm_core::prelude::*;
use std::f64::consts::PI;
use std::sync::{Arc, OnceLock};

fn gpu_ctx_or_skip() -> Option<Arc<GpuContext>> {
    static CTX: OnceLock<Result<Arc<GpuContext>, String>> = OnceLock::new();
    match CTX.get_or_init(|| GpuContext::new().map_err(|e| e.to_string())) {
        Ok(ctx) => Some(ctx.clone()),
        Err(e) => {
            if std::env::var_os("LBM_REQUIRE_GPU").is_some() {
                panic!("T14 D3Q19 absolute GPU tests require an adapter: {e}");
            }
            eprintln!(
                "skipping T14 D3Q19 absolute GPU test: no adapter ({e}); #[ignore = \"gpu-required\"]"
            );
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

fn idx(x: usize, y: usize, z: usize, nx: usize, ny: usize) -> usize {
    (z * ny + y) * nx + x
}

struct TgvMetrics {
    l2: f64,
    rate_rel: f64,
    uz_rel: f64,
    symmetry_rel: f64,
    steps: usize,
}

fn tgv3d_gpu_metrics(n: usize) -> Option<TgvMetrics> {
    let ctx = gpu_ctx_or_skip()?;
    let nu = 0.02f64;
    let u0_coef = 1.28e-4f64;
    let u0 = u0_coef / n as f64;
    let k = 2.0 * PI / n as f64;
    let spec = GlobalSpec::<f32> {
        dims: [n, n, n],
        nu,
        periodic: [true, true, true],
        collision: CollisionKind::Trt { magic: 3.0 / 16.0 },
        ..Default::default()
    };
    let mut gpu = GpuSolver::<D3Q19>::new(&spec, &[], &[], ctx);
    gpu.set_submit_chunk(25);
    let vel = move |x: usize, y: usize, z: usize| -> [f64; 3] {
        let (xf, yf, zf) = (k * x as f64, k * y as f64, k * z as f64);
        [
            u0 * xf.sin() * yf.cos() * zf.cos(),
            -u0 * xf.cos() * yf.sin() * zf.cos(),
            0.0,
        ]
    };
    gpu.init_with(move |x, y, z| {
        let (xf, yf, zf) = (k * x as f64, k * y as f64, k * z as f64);
        let p = u0 * u0 / 16.0 * (((2.0 * xf).cos() + (2.0 * yf).cos()) * ((2.0 * zf).cos() + 2.0));
        let u = vel(x, y, z);
        (1.0 + 3.0 * p as f32, [u[0] as f32, u[1] as f32, 0.0])
    });

    let e0 = {
        let ux = gpu.gather_ux();
        let uy = gpu.gather_uy();
        let uz = gpu.gather_uz();
        ux.iter()
            .zip(&uy)
            .zip(&uz)
            .map(|((&a, &b), &c)| {
                let (a, b, c) = (a as f64, b as f64, c as f64);
                a * a + b * b + c * c
            })
            .sum::<f64>()
    };
    let steps = (0.1 / (nu * k * k)).round() as usize;
    gpu.run(steps);
    let (ux, uy, uz) = (gpu.gather_ux(), gpu.gather_uy(), gpu.gather_uz());
    let e1 = ux
        .iter()
        .zip(&uy)
        .zip(&uz)
        .map(|((&a, &b), &c)| {
            let (a, b, c) = (a as f64, b as f64, c as f64);
            a * a + b * b + c * c
        })
        .sum::<f64>();
    let rate = -(e1 / e0).ln() / steps as f64;
    let rate_ref = 6.0 * nu * k * k;
    let rate_rel = (rate - rate_ref).abs() / rate_ref;

    let decay = (-3.0 * nu * k * k * steps as f64).exp();
    let mut actual = Vec::with_capacity(3 * n * n * n);
    let mut reference = Vec::with_capacity(3 * n * n * n);
    let mut uz_max = 0.0f64;
    let mut sym = 0.0f64;
    let mut sym_scale = 0.0f64;
    for z in 0..n {
        let zm = (n - z) % n;
        for y in 0..n {
            let ym = (n - y) % n;
            for x in 0..n {
                let i = idx(x, y, z, n, n);
                let v = vel(x, y, z);
                let ref_u = [v[0] * decay, v[1] * decay, 0.0];
                actual.push(ux[i] as f64);
                actual.push(uy[i] as f64);
                actual.push(uz[i] as f64);
                reference.extend(ref_u);
                uz_max = uz_max.max((uz[i] as f64).abs());

                let im = idx(x, ym, zm, n, n);
                sym = sym.max((ux[i] as f64 - ux[im] as f64).abs());
                sym = sym.max((uy[i] as f64 + uy[im] as f64).abs());
                sym_scale = sym_scale
                    .max((ux[i] as f64).abs())
                    .max((uy[i] as f64).abs());
            }
        }
    }
    Some(TgvMetrics {
        l2: l2_rel(&actual, &reference),
        rate_rel,
        uz_rel: uz_max / (u0 * decay).max(1e-30),
        symmetry_rel: sym / sym_scale.max(1e-30),
        steps,
    })
}

#[test]
fn t14_gpu_d3q19_tgv3d_matches_analytic_convergence() {
    let Some(m32) = tgv3d_gpu_metrics(32) else {
        return;
    };
    let m64 = tgv3d_gpu_metrics(64).expect("same GPU context should be reusable");
    let order = (m32.l2 / m64.l2).log2();
    println!(
        "T14 GPU D3Q19 TGV3D: e32={:.6e} (steps {}), e64={:.6e} (steps {}), \
         order={order:.3}, rate64_rel={:.6e}, uz64/u0={:.6e}, symmetry64={:.6e}",
        m32.l2, m32.steps, m64.l2, m64.steps, m64.rate_rel, m64.uz_rel, m64.symmetry_rel
    );
    assert!(
        m32.l2 <= 5.0e-3 && m64.l2 <= 5.0e-3,
        "T14 GPU D3Q19 TGV3D L2rel out of band: e32={:.6e}, e64={:.6e}, band=5e-3",
        m32.l2,
        m64.l2
    );
    assert!(
        order >= 1.7,
        "T14 GPU D3Q19 TGV3D convergence order={order:.3} < 1.7, e32={:.6e}, e64={:.6e}",
        m32.l2,
        m64.l2
    );
    assert!(
        m64.rate_rel <= 2.0e-2,
        "T14 GPU D3Q19 TGV3D N=64 decay-rate relative error {:.6e} > 2%",
        m64.rate_rel
    );
    assert!(
        m64.uz_rel <= 1.0e-2,
        "T14 GPU D3Q19 TGV3D behavior anchor failed: |uz|/(u0 decay)={:.6e}",
        m64.uz_rel
    );
    assert!(
        m64.symmetry_rel <= 2.0e-2,
        "T14 GPU D3Q19 TGV3D mirror symmetry drift {:.6e} > 2%",
        m64.symmetry_rel
    );
}

/// Exact series for axial velocity in a rectangular duct with half-widths
/// `a`, `b` and uniform body force `g` along x. This is the T15.2 reference.
fn duct_series(y: f64, zt: f64, a: f64, b: f64, g: f64, nu: f64, nmax: usize) -> f64 {
    let pref = 16.0 * a * a * g / (nu * PI * PI * PI);
    let mut sum = 0.0;
    let mut n = 1;
    while n <= nmax {
        let nf = n as f64;
        let kn = nf * PI / (2.0 * a);
        let ratio =
            ((kn * zt.abs()).exp() + (-kn * zt.abs()).exp()) / ((kn * b).exp() + (-kn * b).exp());
        sum += (1.0 - ratio) * (kn * y).sin() / (nf * nf * nf);
        n += 2;
    }
    pref * sum
}

fn run_duct_until_steady(
    gpu: &mut GpuSolver<D3Q19>,
    check_every: usize,
    tol: f64,
    max_steps: usize,
) -> bool {
    let mut prev: Option<Vec<f32>> = None;
    let mut elapsed = 0;
    while elapsed < max_steps {
        gpu.run(check_every);
        elapsed += check_every;
        let mut cur = gpu.gather_ux();
        cur.extend(gpu.gather_uy());
        cur.extend(gpu.gather_uz());
        if let Some(p) = &prev {
            let mut dmax = 0.0f64;
            let mut umax = 0.0f64;
            for (&a, &b) in cur.iter().zip(p) {
                dmax = dmax.max((a as f64 - b as f64).abs());
                umax = umax.max((a as f64).abs());
            }
            if umax > 0.0 && dmax <= tol * umax {
                return true;
            }
        }
        prev = Some(cur);
    }
    false
}

#[test]
fn t14_gpu_d3q19_duct_poiseuille_matches_series() {
    let Some(ctx) = gpu_ctx_or_skip() else {
        return;
    };
    let (nx, ny, nz) = (4usize, 34usize, 34usize);
    let (hy, hz) = ((ny - 2) as f64, (nz - 2) as f64);
    let (a, b) = (hy / 2.0, hz / 2.0);
    let tau = 0.8f64;
    let nu = (tau - 0.5) / 3.0;
    let g = 5.0e-6f64;
    let mut walls = WallSpec::<f32>::default();
    for face in [Face::YNeg, Face::YPos, Face::ZNeg, Face::ZPos] {
        walls.is_wall[face.index()] = true;
    }
    let spec = GlobalSpec::<f32> {
        dims: [nx, ny, nz],
        nu,
        periodic: [true, false, false],
        force: [g as f32, 0.0, 0.0],
        collision: CollisionKind::Trt { magic: 3.0 / 16.0 },
        ..Default::default()
    };
    let (solid, wall_u) = build_wall_rims(3, spec.dims, &walls);
    let mut gpu = GpuSolver::<D3Q19>::new(&spec, &solid, &wall_u, ctx);
    gpu.set_submit_chunk(100);
    let steady = run_duct_until_steady(&mut gpu, 1_000, 2.0e-6, 80_000);
    assert!(steady, "T14 GPU D3Q19 duct did not reach f32 steady state");

    let ux = gpu.gather_ux();
    let uy = gpu.gather_uy();
    let uz = gpu.gather_uz();
    let mut actual = Vec::with_capacity((ny - 2) * (nz - 2));
    let mut reference = Vec::with_capacity((ny - 2) * (nz - 2));
    let mut center = 0.0f64;
    let mut near_wall = f64::INFINITY;
    let mut cross = 0.0f64;
    let mut mirror = 0.0f64;
    for z in 1..=(nz - 2) {
        for y in 1..=(ny - 2) {
            let yy = y as f64 - 0.5;
            let zt = z as f64 - 0.5 - b;
            let ana = duct_series(yy, zt, a, b, g, nu, 99);
            let mut u = 0.0;
            for x in 0..nx {
                u += ux[idx(x, y, z, nx, ny)] as f64;
            }
            u /= nx as f64;
            actual.push(u);
            reference.push(ana);
            center = center.max(u);
            near_wall = near_wall.min(u);
            cross = cross
                .max((uy[idx(0, y, z, nx, ny)] as f64).abs())
                .max((uz[idx(0, y, z, nx, ny)] as f64).abs());
            let ym = ny - 1 - y;
            let zm = nz - 1 - z;
            mirror = mirror.max((u - ux[idx(0, ym, z, nx, ny)] as f64).abs());
            mirror = mirror.max((u - ux[idx(0, y, zm, nx, ny)] as f64).abs());
        }
    }
    let l2 = l2_rel(&actual, &reference);
    let umax_ref = duct_series(a, 0.0, a, b, g, nu, 99);
    let cross_rel = cross / umax_ref;
    let mirror_rel = mirror / umax_ref;
    println!(
        "T14 GPU D3Q19 duct Poiseuille: L2rel={l2:.6e}, tau={tau:.3}, nu={nu:.6e}, \
         g={g:.6e}, umax_ref={umax_ref:.6e}, center={center:.6e}, \
         near_wall_min={near_wall:.6e}, cross/umax={cross_rel:.6e}, mirror/umax={mirror_rel:.6e}"
    );
    assert!(
        l2 <= 1.0e-2,
        "T14 GPU D3Q19 duct Poiseuille L2rel={l2:.6e} > 1e-2"
    );
    assert!(
        center > 0.5 * umax_ref && near_wall >= 0.0 && center > near_wall,
        "T14 GPU D3Q19 duct behavior anchor failed: center={center:.6e}, near_wall={near_wall:.6e}, umax_ref={umax_ref:.6e}"
    );
    assert!(
        cross_rel <= 1.0e-3,
        "T14 GPU D3Q19 duct cross-flow ratio {cross_rel:.6e} > 1e-3"
    );
    assert!(
        mirror_rel <= 1.0e-2,
        "T14 GPU D3Q19 duct mirror symmetry drift {mirror_rel:.6e} > 1e-2"
    );
}
