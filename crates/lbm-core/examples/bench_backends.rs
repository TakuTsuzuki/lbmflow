//! MLUPS benchmark: `CpuScalar` vs `CpuSimd` on the same TGV-style periodic
//! scenario. Prints a markdown table for docs/PERFORMANCE.md.
//!
//! Run: `cargo run --release -p lbm-core --example bench_backends`
//!
//! Single-config mode (for A/B comparisons under varying machine load —
//! alternate the runs in the same time window, best-of-N):
//! `bench_backends <scalar|simd> <f32|f64> <n> <threads> <steps> [nz]`
//! prints one MLUPS value. `nz` > 1 selects the 3D (D3Q19) `n x n x nz`
//! grid.
//!
//! V1 comparison: the live V1 column was retired with `crates/lbm-core`
//! (2026-07-05). The frozen same-window measurements against the V1 fused
//! kernel are documented in docs/PERFORMANCE.md ("V2 CpuSimd バックエンド"
//! table: e.g. 2D 512²/1T f32 V1 232 vs CpuSimd 273 MLUPS, 1024²/12T f32
//! V1 1084 vs 1183; target "V1 − 10% 以内" met on all configurations).

use lbm_core::lattice::{D2Q9, D3Q19};
use lbm_core::prelude::*;
use std::f64::consts::PI;
use std::time::Instant;

fn spec2d<T: Real>(n: usize) -> GlobalSpec<T> {
    GlobalSpec {
        dims: [n, n, 1],
        nu: 0.02,
        periodic: [true, true, false],
        ..Default::default()
    }
}

fn spec3d<T: Real>(n: usize, nz: usize) -> GlobalSpec<T> {
    GlobalSpec {
        dims: [n, n, nz],
        nu: 0.02,
        periodic: [true, true, true],
        ..Default::default()
    }
}

/// V1 bench_mlups initial condition (2D) / its z-modulated 3D analogue.
fn init<T: Real>(n: usize) -> impl Fn(usize, usize, usize) -> (T, [T; 3]) + Copy {
    move |x, y, z| {
        let k = 2.0 * PI / n as f64;
        let (xf, yf, zf) = (k * x as f64, k * y as f64, k * z as f64);
        (
            T::one(),
            [
                T::r(0.03 * yf.sin()),
                T::r(0.03 * (2.0 * xf).sin()),
                T::r(if z == 0 { 0.0 } else { 0.02 * zf.sin() }),
            ],
        )
    }
}

fn bench_v2<L: Lattice, T: Real, B: Backend<L, T, Fields = SoaFields<T>>>(
    spec: &GlobalSpec<T>,
    backend: B,
    steps: usize,
) -> f64 {
    let mut s: Solver<L, T, B, LocalPeriodic> =
        Solver::new(spec, &[], &[], [1, 1, 1], backend, LocalPeriodic);
    s.init_with(init::<T>(spec.dims[0]));
    s.run(10); // warmup
    let cells = spec.dims[0] * spec.dims[1] * spec.dims[2];
    let t0 = Instant::now();
    s.run(steps);
    (cells * steps) as f64 / t0.elapsed().as_secs_f64() / 1e6
}

fn run_one(engine: &str, prec: &str, n: usize, steps: usize, nz: usize) -> f64 {
    match (engine, prec, nz > 1) {
        ("scalar", "f32", false) => {
            bench_v2::<D2Q9, f32, _>(&spec2d(n), CpuScalar::default(), steps)
        }
        ("scalar", "f64", false) => {
            bench_v2::<D2Q9, f64, _>(&spec2d(n), CpuScalar::default(), steps)
        }
        ("simd", "f32", false) => bench_v2::<D2Q9, f32, _>(&spec2d(n), CpuSimd::default(), steps),
        ("simd", "f64", false) => bench_v2::<D2Q9, f64, _>(&spec2d(n), CpuSimd::default(), steps),
        ("scalar", "f32", true) => {
            bench_v2::<D3Q19, f32, _>(&spec3d(n, nz), CpuScalar::default(), steps)
        }
        ("scalar", "f64", true) => {
            bench_v2::<D3Q19, f64, _>(&spec3d(n, nz), CpuScalar::default(), steps)
        }
        ("simd", "f32", true) => {
            bench_v2::<D3Q19, f32, _>(&spec3d(n, nz), CpuSimd::default(), steps)
        }
        ("simd", "f64", true) => {
            bench_v2::<D3Q19, f64, _>(&spec3d(n, nz), CpuSimd::default(), steps)
        }
        other => panic!("unsupported combination {other:?}"),
    }
}

fn pool(threads: usize) -> rayon::ThreadPool {
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .unwrap()
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    if args.len() >= 5 {
        let engine = args[0].as_str();
        let prec = args[1].as_str();
        let n: usize = args[2].parse().expect("grid size");
        let threads: usize = args[3].parse().expect("thread count");
        let steps: usize = args[4].parse().expect("step count");
        let nz: usize = args.get(5).map_or(1, |s| s.parse().expect("nz"));
        let mlups = pool(threads).install(|| run_one(engine, prec, n, steps, nz));
        println!("{mlups:.1}");
        return;
    }

    // Full table. Interleave the engines per configuration so shared-machine
    // load shifts hit all engines alike (PERFORMANCE.md measurement note).
    // V1-fused reference values are frozen in docs/PERFORMANCE.md (the live
    // column left with crates/lbm-core, 2026-07-05).
    println!("## 2D D2Q9 (TGV-style periodic, TRT) — MLUPS, best of 3\n");
    println!("| grid | threads | prec | CpuScalar | CpuSimd | Simd/Scalar |");
    println!("|---|---|---|---|---|---|");
    for &n in &[512usize, 1024] {
        for &threads in &[1usize, 12] {
            for prec in ["f32", "f64"] {
                let steps = ((100_000_000 / (n * n)) * threads.min(4)).max(30);
                let p = pool(threads);
                let mut best = [0.0f64; 2];
                for _ in 0..3 {
                    for (i, engine) in ["scalar", "simd"].iter().enumerate() {
                        let m = p.install(|| run_one(engine, prec, n, steps, 1));
                        if m > best[i] {
                            best[i] = m;
                        }
                    }
                }
                println!(
                    "| {n}² | {threads} | {prec} | {:.0} | {:.0} | {:.2} |",
                    best[0],
                    best[1],
                    best[1] / best[0],
                );
            }
        }
    }
    println!("\n## 3D D3Q19 (128³ periodic TGV-style, TRT) — MLUPS, best of 3\n");
    println!("| grid | threads | prec | CpuScalar | CpuSimd | Simd/Scalar |");
    println!("|---|---|---|---|---|---|");
    for &threads in &[1usize, 12] {
        for prec in ["f32", "f64"] {
            let n = 128usize;
            let steps = (threads.min(4) * 20).max(10);
            let p = pool(threads);
            let mut best = [0.0f64; 2];
            for _ in 0..3 {
                for (i, engine) in ["scalar", "simd"].iter().enumerate() {
                    let m = p.install(|| run_one(engine, prec, n, steps, n));
                    if m > best[i] {
                        best[i] = m;
                    }
                }
            }
            println!(
                "| {n}³ | {threads} | {prec} | {:.0} | {:.0} | {:.2} |",
                best[0],
                best[1],
                best[1] / best[0],
            );
        }
    }
}
