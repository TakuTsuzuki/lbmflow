//! Phase 9c evaluation harness (docs/REVIEW_2026-07-05.md §3).
//!
//! One-shot run that (a) verifies the fused wgpu D2Q9 kernel against the
//! lbm-core CPU engine on a periodic Taylor–Green vortex — same f32 initial
//! state, 2000 steps, L∞ of the velocity difference — and (b) benchmarks
//! MLUPS at 512² / 1024² / 2048² including a workgroup-size sweep, a
//! submit-granularity sweep, a fresh CPU baseline on the same machine, and
//! the velocity readback cost. Prints markdown tables consumed by
//! docs/GPU_EVALUATION.md.
//!
//! Run: `cd crates/lbm-gpu-proto && cargo run --release`
//!
//! GPU timing convention: wall time of encode + submit + wait-until-idle
//! ("effective" MLUPS — nothing hidden in queue latency), measured after a
//! warmup that absorbs pipeline compilation and first-touch costs.

mod gpu;
mod hostinit;

use gpu::{GpuContext, GpuLbm};
use lbm_core::compat::prelude::{SimConfig, Simulation};
use std::f64::consts::PI;
use std::time::Instant;

/// Kinematic viscosity (tau = 0.56), same as the CPU TGV suite.
const NU: f64 = 0.02;
/// TGV amplitude. Fixed (not diffusive-scaled): the check compares two f32
/// *implementations*, so a healthy signal-to-roundoff ratio matters more
/// than compressibility error, which cancels between the two anyway.
const U0: f64 = 0.05;
const VALIDATION_STEPS: usize = 2000;
/// L∞(Δu)/max‖u_cpu‖ acceptance threshold (f32 vs f32).
const LINF_TOL: f64 = 1e-4;

/// Build a GPU sim seeded with the pre-collided CPU-identical TGV state
/// (see hostinit.rs for the operator-ordering argument).
fn gpu_from_host<'a>(ctx: &'a GpuContext, n: usize, wg: (u32, u32)) -> GpuLbm<'a> {
    let mut h = hostinit::init_tgv_f32(n, U0, NU);
    let (op, om) = hostinit::omegas_f32(NU);
    hostinit::collide_trt_f32(&mut h.f_aos, &h.rho, &h.ux, &h.uy, op, om);
    let soa = hostinit::aos_to_soa(&h.f_aos, n * n);
    let g = GpuLbm::new(ctx, n as u32, n as u32, NU, wg);
    g.upload(&soa);
    g
}

fn cpu_sim(n: usize) -> Simulation<f32> {
    let mut sim: Simulation<f32> = SimConfig {
        nx: n,
        ny: n,
        nu: NU,
        ..Default::default()
    }
    .build()
    .expect("valid config");
    sim.init_with(|x, y| hostinit::tgv_ic(x, y, n, U0));
    sim
}

struct Validation {
    n: usize,
    linf_rel: f64,
    l2_diff_rel: f64,
    gpu_vs_ana: f64,
    cpu_vs_ana: f64,
    pass: bool,
}

fn validate(ctx: &GpuContext, n: usize, steps: usize) -> Validation {
    let mut sim = cpu_sim(n);
    sim.run(steps);

    let mut g = gpu_from_host(ctx, n, (16, 8));
    g.run(steps, 250);
    let vel = g.velocity();

    let k = 2.0 * PI / n as f64;
    let decay = (-2.0 * NU * k * k * steps as f64).exp();
    let (ux_c, uy_c) = (sim.ux_field(), sim.uy_field());
    let mut linf = 0.0f64;
    let mut umax = 0.0f64;
    let (mut nd, mut dd) = (0.0f64, 0.0f64);
    let (mut ng, mut nc, mut da) = (0.0f64, 0.0f64, 0.0f64);
    for y in 0..n {
        for x in 0..n {
            let i = y * n + x;
            let (cux, cuy) = (ux_c[i] as f64, uy_c[i] as f64);
            let (gux, guy) = (vel[2 * i] as f64, vel[2 * i + 1] as f64);
            linf = linf.max((gux - cux).abs()).max((guy - cuy).abs());
            umax = umax.max((cux * cux + cuy * cuy).sqrt());
            nd += (gux - cux).powi(2) + (guy - cuy).powi(2);
            dd += cux * cux + cuy * cuy;
            let (xf, yf) = (k * x as f64, k * y as f64);
            let uxa = -U0 * xf.cos() * yf.sin() * decay;
            let uya = U0 * xf.sin() * yf.cos() * decay;
            ng += (gux - uxa).powi(2) + (guy - uya).powi(2);
            nc += (cux - uxa).powi(2) + (cuy - uya).powi(2);
            da += uxa * uxa + uya * uya;
        }
    }
    let linf_rel = linf / umax;
    Validation {
        n,
        linf_rel,
        l2_diff_rel: (nd / dd).sqrt(),
        gpu_vs_ana: (ng / da).sqrt(),
        cpu_vs_ana: (nc / da).sqrt(),
        pass: linf_rel < LINF_TOL,
    }
}

/// Warm up, auto-calibrate the step count to roughly `target_s` seconds,
/// then measure. Timing includes encode + submit + final wait.
fn bench_gpu(
    ctx: &GpuContext,
    n: usize,
    wg: (u32, u32),
    chunk: usize,
    target_s: f64,
    wait_each_submit: bool,
) -> f64 {
    let mut g = gpu_from_host(ctx, n, wg);
    g.run(50, chunk); // warmup: pipeline compile, first touch
    let cells = (n * n) as f64;
    let t = Instant::now();
    g.run_opts(100, chunk, wait_each_submit);
    let rate = 100.0 * cells / t.elapsed().as_secs_f64().max(1e-9);
    let steps = ((target_s * rate / cells) as usize).clamp(100, 60_000);
    let t = Instant::now();
    g.run_opts(steps, chunk, wait_each_submit);
    cells * steps as f64 / t.elapsed().as_secs_f64() / 1e6
}

/// CPU baseline on the same machine: lbm-core f32/TRT, rayon default threads.
/// Best of two reps — the max is the least contaminated by other processes
/// (this box also runs other agents' test suites; see docs/GPU_EVALUATION.md).
fn bench_cpu(n: usize, target_s: f64) -> f64 {
    let mut sim = cpu_sim(n);
    sim.run(10); // warmup
    let cells = (n * n) as f64;
    let t = Instant::now();
    sim.run(20);
    let rate = 20.0 * cells / t.elapsed().as_secs_f64().max(1e-9);
    let steps = ((target_s * rate / cells) as usize).clamp(20, 10_000);
    let mut best = 0.0f64;
    for _ in 0..2 {
        let t = Instant::now();
        sim.run(steps);
        best = best.max(cells * steps as f64 / t.elapsed().as_secs_f64() / 1e6);
    }
    best
}

/// Average blocking cost of one velocity-field readback
/// (moments dispatch + buffer copy + map + memcpy to host Vec).
fn bench_readback(ctx: &GpuContext, n: usize, wg: (u32, u32)) -> f64 {
    let mut g = gpu_from_host(ctx, n, wg);
    g.run(10, 10);
    std::hint::black_box(g.velocity()); // warm
    let reps = 10;
    let t = Instant::now();
    for _ in 0..reps {
        std::hint::black_box(g.velocity());
    }
    t.elapsed().as_secs_f64() / reps as f64 * 1e3
}

fn main() {
    // `--cpu-only` measures just the lbm-core baseline table (useful to grab
    // a clean sample on a busy box); `--gpu-only` skips the CPU baselines.
    let args: Vec<String> = std::env::args().skip(1).collect();
    let cpu_only = args.iter().any(|a| a == "--cpu-only");
    let gpu_only = args.iter().any(|a| a == "--gpu-only");

    if cpu_only {
        println!("## CPU baseline (lbm-core f32/TRT, rayon default threads)\n");
        println!("| grid | CPU MLUPS |");
        println!("|---|---|");
        for n in [512usize, 1024, 2048] {
            println!("| {n}x{n} | {:.0} |", bench_cpu(n, 1.0));
        }
        return;
    }

    println!("# lbm-gpu-proto — Phase 9c wgpu D2Q9 evaluation\n");
    let Some(ctx) = GpuContext::new() else {
        println!(
            "RESULT: no usable GPU adapter (see stderr). This is itself the \
             evaluation outcome for this machine: the wgpu backend cannot run here."
        );
        std::process::exit(2);
    };
    let info = &ctx.adapter_info;
    println!("- adapter: {} / backend {:?}", info.name, info.backend);
    println!(
        "- physics: D2Q9 TRT (magic 3/16), f32 deviation storage, periodic TGV \
         (nu={NU}, u0={U0})"
    );
    println!(
        "- GPU timing: encode+submit+wait inclusive, after warmup; \
         CPU reference: lbm-core f32/TRT, rayon default threads"
    );

    // ---------------------------------------------------------------- 1
    println!("\n## Verification vs lbm-core ({VALIDATION_STEPS} steps, f32 vs f32)\n");
    println!("| grid | Linf(du)/max|u| | L2 rel diff | GPU vs analytic L2 | CPU vs analytic L2 | verdict |");
    println!("|---|---|---|---|---|---|");
    let mut all_pass = true;
    for n in [256usize, 512] {
        let v = validate(&ctx, n, VALIDATION_STEPS);
        all_pass &= v.pass;
        println!(
            "| {}x{} | {:.2e} | {:.2e} | {:.3e} | {:.3e} | {} |",
            v.n,
            v.n,
            v.linf_rel,
            v.l2_diff_rel,
            v.gpu_vs_ana,
            v.cpu_vs_ana,
            if v.pass { "PASS" } else { "FAIL" }
        );
    }
    println!("\n(threshold: Linf(du)/max|u| < {LINF_TOL:.0e})");

    // ---------------------------------------------------------------- 2
    println!("\n## Workgroup-size sweep (1024x1024, chunk=100)\n");
    println!("| workgroup | MLUPS |");
    println!("|---|---|");
    let mut best = (16u32, 8u32);
    let mut best_m = 0.0f64;
    for wg in [
        (8u32, 8u32),
        (16, 4),
        (32, 2),
        (64, 1),
        (16, 8),
        (32, 4),
        (128, 1),
        (16, 16),
        (32, 8),
        (256, 1),
    ] {
        let m = bench_gpu(&ctx, 1024, wg, 100, 0.4, false);
        if m > best_m {
            best_m = m;
            best = wg;
        }
        println!("| {}x{} | {m:.0} |", wg.0, wg.1);
    }
    println!("\nbest: {}x{} ({best_m:.0} MLUPS)", best.0, best.1);

    // ---------------------------------------------------------------- 3
    println!("\n## Submit granularity (1024x1024, wg {}x{})\n", best.0, best.1);
    println!("| dispatches per submit | wait per submit | MLUPS |");
    println!("|---|---|---|");
    for (chunk, wait_each) in [
        (1usize, true),
        (1, false),
        (10, false),
        (100, false),
        (500, false),
    ] {
        let m = bench_gpu(&ctx, 1024, best, chunk, 0.4, wait_each);
        println!("| {chunk} | {} | {m:.0} |", if wait_each { "yes" } else { "no" });
    }

    // ---------------------------------------------------------------- 4
    println!(
        "\n## MLUPS by grid (TRT f32, wg {}x{}, chunk=100, submit->wait inclusive)\n",
        best.0, best.1
    );
    if gpu_only {
        println!("| grid | GPU MLUPS |");
        println!("|---|---|");
        for n in [512usize, 1024, 2048] {
            let gm = bench_gpu(&ctx, n, best, 100, 1.2, false);
            println!("| {n}x{n} | {gm:.0} |");
        }
    } else {
        println!("| grid | GPU MLUPS | CPU MLUPS (same machine) | speedup |");
        println!("|---|---|---|---|");
        for n in [512usize, 1024, 2048] {
            let gm = bench_gpu(&ctx, n, best, 100, 1.2, false);
            let cm = bench_cpu(n, 1.0);
            println!("| {n}x{n} | {gm:.0} | {cm:.0} | {:.1}x |", gm / cm);
        }
    }

    // ---------------------------------------------------------------- 5
    println!("\n## Velocity-field readback (moments dispatch + copy + map, blocking)\n");
    println!("| grid | payload | time |");
    println!("|---|---|---|");
    for n in [512usize, 1024, 2048] {
        let ms = bench_readback(&ctx, n, best);
        println!(
            "| {n}x{n} | {:.1} MB | {ms:.2} ms |",
            (n * n * 8) as f64 / 1e6
        );
    }

    if all_pass {
        println!("\nAll verification checks passed.");
    } else {
        println!("\nVERIFICATION FAILED — see table above.");
        std::process::exit(1);
    }
}
