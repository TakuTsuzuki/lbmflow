//! MLUPS benchmark across precision, collision operator, grid size and
//! thread count. Prints a markdown table for docs/PERFORMANCE.md.
//!
//! Run: `cargo run --release --example bench_mlups`

use lbm_core::prelude::*;
use std::f64::consts::PI;
use std::time::Instant;

fn build<T: Real>(n: usize, collision: Collision) -> Simulation<T> {
    let mut sim: Simulation<T> = SimConfig {
        nx: n,
        ny: n,
        nu: 0.02,
        collision,
        ..Default::default()
    }
    .build()
    .unwrap();
    let k = 2.0 * PI / n as f64;
    sim.init_with(|x, y| {
        let (xf, yf) = (k * x as f64, k * y as f64);
        (
            T::one(),
            T::r(0.03 * yf.sin()),
            T::r(0.03 * (2.0 * xf).sin()),
        )
    });
    sim
}

fn bench<T: Real>(n: usize, collision: Collision, steps: usize) -> f64 {
    let mut sim = build::<T>(n, collision);
    sim.run(20); // warmup
    let t0 = Instant::now();
    sim.run(steps);
    let dt = t0.elapsed().as_secs_f64();
    (n * n * steps) as f64 / dt / 1e6
}

fn main() {
    let trt = Collision::default();
    let bgk = Collision::Bgk;

    println!("## スレッドスケーリング (512^2, TRT)\n");
    println!("| threads | f32 MLUPS | f64 MLUPS |");
    println!("|---|---|---|");
    for threads in [1usize, 2, 4, 8, 12, 18] {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap();
        let (a, b) = pool.install(|| (bench::<f32>(512, trt, 300), bench::<f64>(512, trt, 300)));
        println!("| {threads} | {a:.0} | {b:.0} |");
    }

    println!("\n## 格子サイズ・演算子・精度 (全スレッド)\n");
    println!("| grid | collision | f32 MLUPS | f64 MLUPS |");
    println!("|---|---|---|---|");
    for n in [256usize, 512, 1024] {
        let steps = (400 * 512 * 512 / (n * n)).max(50);
        for (c, name) in [(bgk, "BGK"), (trt, "TRT")] {
            let a = bench::<f32>(n, c, steps);
            let b = bench::<f64>(n, c, steps);
            println!("| {n}x{n} | {name} | {a:.0} | {b:.0} |");
        }
    }

    println!("\n## シングルスレッド比較用 (256^2, TRT, serial path)\n");
    // below PARALLEL_MIN_CELLS the engine runs serially; emulate by 1-thread pool
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(1)
        .build()
        .unwrap();
    let (a, b) = pool.install(|| (bench::<f32>(256, trt, 400), bench::<f64>(256, trt, 400)));
    println!("1-thread 256^2: f32 {a:.0} MLUPS / f64 {b:.0} MLUPS");
}
