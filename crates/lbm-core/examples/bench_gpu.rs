//! GPU backend benchmark: MLUPS at 512²/1024²/2048² vs the CpuScalar
//! backend and vs the lbm-gpu-proto reference numbers (abstraction-overhead
//! check; docs/GPU_EVALUATION.md §1).
//!
//! Run: `cargo run -p lbm-core --release --features gpu --example bench_gpu`
//! (`--gpu-only` skips the CPU baselines.)
//!
//! Timing convention matches the proto: wall time of encode + submit +
//! wait-until-idle after a warmup, so nothing hides in queue latency.
//!
//! **Measurement hygiene** (GPU_EVALUATION.md's own caveat): this box runs
//! other agents' suites; on unified memory a saturated CPU eats the DRAM
//! bandwidth the (bandwidth-bound) GPU kernel needs, so compare against a
//! *same-window* proto run (`cd crates/lbm-gpu-proto && cargo run --release
//! -- --gpu-only`), not only the frozen table. Measurements 2026-07-05:
//! quieter window (load ~13-16): 11,365 / 6,808 / 5,857 MLUPS = −6.5% /
//! −10.2% / −16.0% vs the *frozen* proto table (−8.9% / −7.9% vs proto run
//! back-to-back, which itself hit 12,478 / 7,395). Heavily loaded window
//! (load ~38): 10,581 vs proto-same-window 11,307 (−6.4%) / 5,435 vs 6,086
//! (−10.7%) / 4,623 vs 5,349 (−13.6%). Both inside ±20%. The residual gap
//! is the push-form fused kernel's scatter writes (vs the proto's gather
//! reads) — the deliberate trade that preserves the CPU's S∘C operator
//! order and makes the boundary-condition set + direct T14 equivalence
//! possible. Measured non-causes: the mask reads (skipping all 8 neighbour
//! mask loads changes nothing — cache-served) and the workgroup shape
//! (256×1 beat 128×1 / 64×2 / 32×4 for the push kernel).

use lbm_core::lattice::D2Q9;
use lbm_core::prelude::*;
use std::time::Instant;

const NU: f64 = 0.02;
const U0: f64 = 0.05;

/// lbm-gpu-proto reference MLUPS on this machine (GPU_EVALUATION.md §1,
/// M5 Max / Metal, wgpu 26, TRT f32, submit→wait inclusive).
const PROTO_REF: [(usize, f64); 3] = [(512, 12_152.0), (1024, 7_584.0), (2048, 6_975.0)];

fn tgv_ic(x: usize, y: usize, n: usize) -> (f32, [f32; 3]) {
    let k = 2.0 * std::f64::consts::PI / n as f64;
    let (xf, yf) = (k * x as f64, k * y as f64);
    let rho = 1.0 - 3.0 * U0 * U0 / 4.0 * ((2.0 * xf).cos() + (2.0 * yf).cos());
    (
        rho as f32,
        [
            (-U0 * xf.cos() * yf.sin()) as f32,
            (U0 * xf.sin() * yf.cos()) as f32,
            0.0,
        ],
    )
}

fn spec(n: usize) -> GlobalSpec<f32> {
    GlobalSpec {
        dims: [n, n, 1],
        nu: NU,
        periodic: [true, true, false],
        ..Default::default()
    }
}

fn gpu_solver(ctx: &std::sync::Arc<GpuContext>, n: usize) -> GpuSolver<D2Q9> {
    let mut s = GpuSolver::<D2Q9>::new(&spec(n), &[], &[], ctx.clone());
    s.init_with(|x, y, _| tgv_ic(x, y, n));
    s
}

/// Warm up, calibrate the step count to ~`target_s` seconds, measure.
fn bench_gpu(ctx: &std::sync::Arc<GpuContext>, n: usize, target_s: f64) -> f64 {
    let mut s = gpu_solver(ctx, n);
    s.run(50);
    s.wait_idle(); // absorb pipeline compile + first touch
    let cells = (n * n) as f64;
    let t = Instant::now();
    s.run(100);
    s.wait_idle();
    let rate = 100.0 * cells / t.elapsed().as_secs_f64().max(1e-9);
    let steps = ((target_s * rate / cells) as usize).clamp(100, 60_000);
    let t = Instant::now();
    s.run(steps);
    s.wait_idle();
    let mlups = cells * steps as f64 / t.elapsed().as_secs_f64() / 1e6;
    // Sanity: the field must still be finite (catches a silently broken
    // kernel masquerading as a fast one).
    let m = s.total_mass();
    assert!(m.is_finite(), "diverged during benchmark");
    mlups
}

/// CpuScalar baseline (rayon all-cores when the `parallel` feature is on).
fn bench_cpu(n: usize, target_s: f64) -> f64 {
    let mut s: Solver<D2Q9, f32, CpuScalar, LocalPeriodic> = Solver::new(
        &spec(n),
        &[],
        &[],
        [1, 1, 1],
        CpuScalar::default(),
        LocalPeriodic,
    );
    s.init_with(|x, y, _| tgv_ic(x, y, n));
    s.run(10);
    let cells = (n * n) as f64;
    let t = Instant::now();
    s.run(20);
    let rate = 20.0 * cells / t.elapsed().as_secs_f64().max(1e-9);
    let steps = ((target_s * rate / cells) as usize).clamp(20, 10_000);
    let mut best = 0.0f64;
    for _ in 0..2 {
        let t = Instant::now();
        s.run(steps);
        best = best.max(cells * steps as f64 / t.elapsed().as_secs_f64() / 1e6);
    }
    best
}

fn main() {
    let gpu_only = std::env::args().any(|a| a == "--gpu-only");
    let ctx = match GpuContext::new() {
        Ok(ctx) => ctx,
        Err(e) => {
            eprintln!("bench_gpu requires a usable GPU adapter: {e}");
            std::process::exit(2);
        }
    };
    println!("# lbm-core Wgpu backend benchmark\n");
    println!(
        "- adapter: {} / backend {:?}",
        ctx.adapter_info.name, ctx.adapter_info.backend
    );
    println!(
        "- physics: D2Q9 TRT (magic 3/16), f32 deviation storage, periodic TGV \
         (nu={NU}, u0={U0}); timing = encode+submit+wait after warmup"
    );
    println!(
        "- proto reference: crates/lbm-gpu-proto on the same machine \
         (GPU_EVALUATION.md §1); acceptance line ±20%\n"
    );
    if gpu_only {
        println!("| grid | GPU MLUPS | proto MLUPS | vs proto |");
        println!("|---|---|---|---|");
    } else {
        println!("| grid | GPU MLUPS | proto MLUPS | vs proto | CPU MLUPS | GPU/CPU |");
        println!("|---|---|---|---|---|---|");
    }
    for (n, proto) in PROTO_REF {
        let g = bench_gpu(&ctx, n, 1.2);
        let vs = 100.0 * (g / proto - 1.0);
        if gpu_only {
            println!("| {n}x{n} | {g:.0} | {proto:.0} | {vs:+.1}% |");
        } else {
            let c = bench_cpu(n, 1.0);
            println!(
                "| {n}x{n} | {g:.0} | {proto:.0} | {vs:+.1}% | {c:.0} | {:.1}x |",
                g / c
            );
        }
    }
}
