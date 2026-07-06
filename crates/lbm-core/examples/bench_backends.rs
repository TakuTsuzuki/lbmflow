//! MLUPS benchmark: `CpuScalar` vs `CpuSimd` on the same TGV-style periodic
//! scenario. Prints a markdown table for docs/PERFORMANCE.md.
//!
//! Run: `cargo run --release -p lbm-core --example bench_backends`
//! Check a committed host baseline:
//! `cargo run --release -p lbm-core --example bench_backends -- --check --host-tag m5-max-sandbox --timestamp 2026-07-07T00:00:00Z`
//! Seed/update a baseline explicitly:
//! `cargo run --release -p lbm-core --example bench_backends -- --update-baseline --host-tag m5-max-sandbox --timestamp 2026-07-07T00:00:00Z`
//!
//! Single-config mode (for A/B comparisons under varying machine load —
//! alternate the runs in the same time window, best-of-N):
//! `bench_backends <scalar|simd> <f32|f64> <n> <threads> <steps> [nz]`
//! prints one MLUPS value. `nz` > 1 selects the 3D (D3Q19) `n x n x nz`
//! grid.
//!
//! V1 comparison: the live V1 column was retired with `crates/lbm-core`
//! (2026-07-05). The frozen same-window measurements against the V1 fused
//! kernel are documented in docs/PERFORMANCE.md ("V2 CpuSimd backend"
//! table: e.g. 2D 512²/1T f32 V1 232 vs CpuSimd 273 MLUPS, 1024²/12T f32
//! V1 1084 vs 1183; target "within V1 − 10%" met on all configurations).

use lbm_core::bench_regression::{
    compare_measurements, BenchBaseline, BenchCaseBaseline, BenchMeasurement, ComparisonFailure,
    DEFAULT_REGRESSION_THRESHOLD,
};
use lbm_core::lattice::{D2Q9, D3Q19};
use lbm_core::prelude::*;
use std::f64::consts::PI;
use std::fs::{self, OpenOptions};
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::Command;
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

#[derive(Clone, Copy, Debug)]
struct BenchCase {
    engine: &'static str,
    precision: &'static str,
    n: usize,
    threads: usize,
    steps: usize,
    nz: usize,
    repeats: usize,
}

impl BenchCase {
    fn id(self) -> String {
        let dim = if self.nz > 1 {
            format!("{}x{}x{}", self.n, self.n, self.nz)
        } else {
            format!("{}x{}", self.n, self.n)
        };
        let lattice = if self.nz > 1 { "d3q19" } else { "d2q9" };
        format!(
            "{lattice}-{}-{}-{dim}-{}t",
            self.engine, self.precision, self.threads
        )
    }
}

fn standard_cases() -> Vec<BenchCase> {
    let mut cases = Vec::new();
    for &n in &[512usize, 1024] {
        for &threads in &[1usize, 12] {
            for precision in ["f32", "f64"] {
                let steps = ((100_000_000 / (n * n)) * threads.min(4)).max(30);
                for engine in ["scalar", "simd"] {
                    cases.push(BenchCase {
                        engine,
                        precision,
                        n,
                        threads,
                        steps,
                        nz: 1,
                        repeats: 3,
                    });
                }
            }
        }
    }
    for &threads in &[1usize, 12] {
        for precision in ["f32", "f64"] {
            for engine in ["scalar", "simd"] {
                cases.push(BenchCase {
                    engine,
                    precision,
                    n: 128,
                    threads,
                    steps: (threads.min(4) * 20).max(10),
                    nz: 128,
                    repeats: 3,
                });
            }
        }
    }
    cases
}

fn run_standard_cases(verbose: bool) -> Vec<BenchMeasurement> {
    let mut measured = Vec::new();
    for cases in standard_cases().chunks_exact(2) {
        let p = pool(cases[0].threads);
        let mut best = [0.0f64; 2];
        for _ in 0..cases[0].repeats {
            for (i, case) in cases.iter().enumerate() {
                let mlups =
                    p.install(|| run_one(case.engine, case.precision, case.n, case.steps, case.nz));
                if mlups > best[i] {
                    best[i] = mlups;
                }
            }
        }
        for (case, mlups) in cases.iter().zip(best) {
            if verbose {
                println!("{} {:.1} MLUPS", case.id(), mlups);
            }
            measured.push(BenchMeasurement {
                case: case.id(),
                mlups,
            });
        }
    }
    measured
}

#[derive(Debug)]
struct BenchOptions {
    mode: BenchMode,
    host_tag: String,
    timestamp: String,
    baseline_path: PathBuf,
    history_path: PathBuf,
    results_dir: PathBuf,
    note: String,
}

#[derive(Clone, Copy, Debug)]
enum BenchMode {
    Check,
    UpdateBaseline,
}

fn usage() -> &'static str {
    "usage:\n  bench_backends\n  bench_backends <scalar|simd> <f32|f64> <n> <threads> <steps> [nz]\n  bench_backends --check --host-tag <tag> --timestamp <ts> [--baseline <path>] [--history <path>]\n  bench_backends --update-baseline --host-tag <tag> --timestamp <ts> [--baseline <path>] [--history <path>] [--note <text>]\n\nTimestamp may also be provided with LBM_BENCH_TIMESTAMP."
}

fn parse_flag_value(args: &[String], index: &mut usize, flag: &str) -> String {
    *index += 1;
    args.get(*index)
        .unwrap_or_else(|| panic!("{flag} requires a value"))
        .clone()
}

fn parse_options(args: &[String]) -> BenchOptions {
    let mut check = false;
    let mut update_baseline = false;
    let mut host_tag: Option<String> = None;
    let mut timestamp: Option<String> = None;
    let mut baseline_path: Option<PathBuf> = None;
    let mut history_path = PathBuf::from("bench/history.csv");
    let mut results_dir = PathBuf::from("bench/results");
    let mut note = "Generated by bench_backends. Thresholds default to 10% unless case entries override them; the PM should tighten per-host baselines from real variance data.".to_string();

    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--check" => check = true,
            "--update-baseline" => update_baseline = true,
            "--host-tag" => host_tag = Some(parse_flag_value(args, &mut i, "--host-tag")),
            "--timestamp" => timestamp = Some(parse_flag_value(args, &mut i, "--timestamp")),
            "--baseline" => {
                baseline_path = Some(PathBuf::from(parse_flag_value(args, &mut i, "--baseline")))
            }
            "--history" => {
                history_path = PathBuf::from(parse_flag_value(args, &mut i, "--history"));
            }
            "--results-dir" => {
                results_dir = PathBuf::from(parse_flag_value(args, &mut i, "--results-dir"));
            }
            "--note" => note = parse_flag_value(args, &mut i, "--note"),
            "--help" | "-h" => {
                println!("{}", usage());
                std::process::exit(0);
            }
            other => panic!("unknown benchmark flag {other:?}\n{}", usage()),
        }
        i += 1;
    }

    let mode = match (check, update_baseline) {
        (true, false) => BenchMode::Check,
        (false, true) => BenchMode::UpdateBaseline,
        (true, true) => panic!("choose only one of --check or --update-baseline"),
        (false, false) => panic!("missing --check or --update-baseline"),
    };
    let host_tag = host_tag.expect("--host-tag is required for check/update modes");
    validate_host_tag(&host_tag);
    let timestamp = timestamp
        .or_else(|| std::env::var("LBM_BENCH_TIMESTAMP").ok())
        .expect("--timestamp or LBM_BENCH_TIMESTAMP is required");
    let baseline_path =
        baseline_path.unwrap_or_else(|| PathBuf::from(format!("bench/baselines/{host_tag}.json")));

    BenchOptions {
        mode,
        host_tag,
        timestamp,
        baseline_path,
        history_path,
        results_dir,
        note,
    }
}

fn validate_host_tag(host_tag: &str) {
    assert!(!host_tag.is_empty(), "--host-tag must not be empty");
    assert!(
        host_tag
            .bytes()
            .all(|b| b.is_ascii_alphanumeric() || matches!(b, b'.' | b'_' | b'-')),
        "--host-tag may only contain ASCII letters, digits, '.', '_', and '-'"
    );
}

fn current_git_commit() -> String {
    let output = Command::new("git")
        .args(["rev-parse", "HEAD"])
        .output()
        .ok()
        .filter(|output| output.status.success())
        .and_then(|output| String::from_utf8(output.stdout).ok())
        .map(|hash| hash.trim().to_string())
        .filter(|hash| !hash.is_empty());
    output
        .or_else(|| std::env::var("GITHUB_SHA").ok())
        .unwrap_or_else(|| "unknown".to_string())
}

fn append_history(
    path: &Path,
    timestamp: &str,
    host_tag: &str,
    git_commit: &str,
    measured: &[BenchMeasurement],
) {
    if let Some(parent) = path.parent() {
        fs::create_dir_all(parent).expect("create benchmark history directory");
    }
    let write_header = !path.exists();
    let mut file = OpenOptions::new()
        .create(true)
        .append(true)
        .open(path)
        .expect("open benchmark history");
    if write_header {
        writeln!(file, "timestamp,host_tag,git_commit,case,mlups")
            .expect("write benchmark history header");
    }
    for measurement in measured {
        writeln!(
            file,
            "{},{},{},{},{:.3}",
            csv_escape(timestamp),
            csv_escape(host_tag),
            csv_escape(git_commit),
            csv_escape(&measurement.case),
            measurement.mlups
        )
        .expect("append benchmark history row");
    }
}

fn csv_escape(value: &str) -> String {
    if value
        .bytes()
        .any(|byte| matches!(byte, b',' | b'"' | b'\n' | b'\r'))
    {
        format!("\"{}\"", value.replace('"', "\"\""))
    } else {
        value.to_string()
    }
}

#[derive(serde::Serialize)]
struct BenchRunResult<'a> {
    schema: u32,
    mode: &'a str,
    host_tag: &'a str,
    timestamp: &'a str,
    git_commit: &'a str,
    cases: &'a [BenchMeasurement],
    failures: &'a [ComparisonFailure],
}

fn write_results(
    options: &BenchOptions,
    mode: &str,
    git_commit: &str,
    measured: &[BenchMeasurement],
    failures: &[ComparisonFailure],
) {
    fs::create_dir_all(&options.results_dir).expect("create benchmark results directory");
    let path = options
        .results_dir
        .join(format!("{}-latest.json", options.host_tag));
    let result = BenchRunResult {
        schema: 1,
        mode,
        host_tag: &options.host_tag,
        timestamp: &options.timestamp,
        git_commit,
        cases: measured,
        failures,
    };
    let text = serde_json::to_string_pretty(&result).expect("serialize benchmark result");
    fs::write(path, format!("{text}\n")).expect("write benchmark result");
}

fn update_baseline(options: &BenchOptions, git_commit: &str, measured: &[BenchMeasurement]) {
    if let Some(parent) = options.baseline_path.parent() {
        fs::create_dir_all(parent).expect("create benchmark baseline directory");
    }
    let baseline = BenchBaseline {
        schema: 1,
        host_tag: options.host_tag.clone(),
        generated_at: options.timestamp.clone(),
        git_commit: git_commit.to_string(),
        note: options.note.clone(),
        cases: measured
            .iter()
            .map(|measurement| BenchCaseBaseline {
                case: measurement.case.clone(),
                mlups: measurement.mlups,
                regression_threshold: Some(DEFAULT_REGRESSION_THRESHOLD),
            })
            .collect(),
    };
    let text = serde_json::to_string_pretty(&baseline).expect("serialize benchmark baseline");
    fs::write(&options.baseline_path, format!("{text}\n")).expect("write benchmark baseline");
}

fn run_check_or_update(options: BenchOptions) {
    let git_commit = current_git_commit();
    if matches!(options.mode, BenchMode::Check) && !options.baseline_path.exists() {
        eprintln!(
            "benchmark baseline does not exist: {}",
            options.baseline_path.display()
        );
        std::process::exit(2);
    }

    let measured = run_standard_cases(true);
    append_history(
        &options.history_path,
        &options.timestamp,
        &options.host_tag,
        &git_commit,
        &measured,
    );

    match options.mode {
        BenchMode::UpdateBaseline => {
            update_baseline(&options, &git_commit, &measured);
            write_results(&options, "update_baseline", &git_commit, &measured, &[]);
            println!(
                "updated benchmark baseline {} with {} cases",
                options.baseline_path.display(),
                measured.len()
            );
        }
        BenchMode::Check => {
            let text = fs::read_to_string(&options.baseline_path).expect("read benchmark baseline");
            let baseline: BenchBaseline =
                serde_json::from_str(&text).expect("parse benchmark baseline");
            let failures =
                compare_measurements(&baseline.cases, &measured, DEFAULT_REGRESSION_THRESHOLD);
            write_results(&options, "check", &git_commit, &measured, &failures);
            if failures.is_empty() {
                println!(
                    "benchmark check passed against {}",
                    options.baseline_path.display()
                );
            } else {
                eprintln!(
                    "benchmark check failed against {}:",
                    options.baseline_path.display()
                );
                for failure in &failures {
                    match failure {
                        ComparisonFailure::Regression {
                            case,
                            baseline_mlups,
                            measured_mlups,
                            threshold,
                        } => eprintln!(
                            "  {case}: {:.1} MLUPS < {:.1} MLUPS baseline with {:.1}% threshold",
                            measured_mlups,
                            baseline_mlups,
                            threshold * 100.0
                        ),
                        ComparisonFailure::MissingMeasurement { case } => {
                            eprintln!("  {case}: baseline case was not measured")
                        }
                        ComparisonFailure::MissingBaseline { case } => {
                            eprintln!("  {case}: measured case is missing from baseline")
                        }
                    }
                }
                std::process::exit(1);
            }
        }
    }
}

fn print_standard_table() {
    let measured = run_standard_cases(false);

    println!("\n## 2D D2Q9 (TGV-style periodic, TRT) - MLUPS, best of 3\n");
    println!("| grid | threads | prec | CpuScalar | CpuSimd | Simd/Scalar |");
    println!("|---|---|---|---|---|---|");
    for &n in &[512usize, 1024] {
        for &threads in &[1usize, 12] {
            for precision in ["f32", "f64"] {
                let scalar = measurement(
                    &measured,
                    &format!("d2q9-scalar-{precision}-{n}x{n}-{threads}t"),
                );
                let simd = measurement(
                    &measured,
                    &format!("d2q9-simd-{precision}-{n}x{n}-{threads}t"),
                );
                println!(
                    "| {n}^2 | {threads} | {precision} | {:.0} | {:.0} | {:.2} |",
                    scalar,
                    simd,
                    simd / scalar,
                );
            }
        }
    }
    println!("\n## 3D D3Q19 (128^3 periodic TGV-style, TRT) - MLUPS, best of 3\n");
    println!("| grid | threads | prec | CpuScalar | CpuSimd | Simd/Scalar |");
    println!("|---|---|---|---|---|---|");
    for &threads in &[1usize, 12] {
        for precision in ["f32", "f64"] {
            let scalar = measurement(
                &measured,
                &format!("d3q19-scalar-{precision}-128x128x128-{threads}t"),
            );
            let simd = measurement(
                &measured,
                &format!("d3q19-simd-{precision}-128x128x128-{threads}t"),
            );
            println!(
                "| 128^3 | {threads} | {precision} | {:.0} | {:.0} | {:.2} |",
                scalar,
                simd,
                simd / scalar,
            );
        }
    }
}

fn measurement(measured: &[BenchMeasurement], case: &str) -> f64 {
    measured
        .iter()
        .find(|measurement| measurement.case == case)
        .unwrap_or_else(|| panic!("missing benchmark measurement {case}"))
        .mlups
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    if args.first().is_some_and(|arg| arg.starts_with("--")) {
        run_check_or_update(parse_options(&args));
        return;
    }

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
    if !args.is_empty() {
        panic!("{}", usage());
    }

    // Full table. Interleave the engines per configuration so shared-machine
    // load shifts hit all engines alike (PERFORMANCE.md measurement note).
    // V1-fused reference values are frozen in docs/PERFORMANCE.md (the live
    // column left with crates/lbm-core, 2026-07-05).
    print_standard_table();
}
