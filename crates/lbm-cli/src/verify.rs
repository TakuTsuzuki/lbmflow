use anyhow::Result;
use clap::ValueEnum;
use lbm_core::compat::prelude::*;
use serde::Serialize;
use std::f64::consts::PI;

#[derive(Clone, Copy, Debug, ValueEnum)]
pub enum VerifyTier {
    Quick,
    Gpu,
    Mpi,
    Full,
}

#[derive(Clone, Copy, Debug, Serialize, PartialEq, Eq)]
#[serde(rename_all = "lowercase")]
enum CheckStatus {
    Pass,
    Fail,
    Skip,
}

#[derive(Clone, Debug, Serialize)]
#[serde(rename_all = "camelCase")]
struct CheckOutcome {
    name: &'static str,
    tier: &'static str,
    status: CheckStatus,
    metric: &'static str,
    value: Option<f64>,
    band: &'static str,
    source: &'static str,
    message: Option<String>,
}

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
struct VerifySummary {
    tier: String,
    ran: usize,
    passed: usize,
    failed: usize,
    skipped: usize,
    checks: Vec<CheckOutcome>,
}

pub fn run(tier: VerifyTier) -> Result<i32> {
    let mut checks = Vec::new();
    match tier {
        VerifyTier::Quick => checks.extend(run_quick_cpu()),
        VerifyTier::Gpu => checks.extend(run_gpu_tier()),
        VerifyTier::Mpi => checks.extend(run_mpi_tier()),
        VerifyTier::Full => {
            checks.extend(run_quick_cpu());
            checks.extend(run_gpu_tier());
        }
    }
    let summary = summarize(format!("{tier:?}").to_lowercase(), checks);
    print_summary(&summary);
    Ok(exit_code(&summary))
}

fn summarize(tier: String, checks: Vec<CheckOutcome>) -> VerifySummary {
    let ran = checks
        .iter()
        .filter(|c| c.status != CheckStatus::Skip)
        .count();
    let passed = checks
        .iter()
        .filter(|c| c.status == CheckStatus::Pass)
        .count();
    let failed = checks
        .iter()
        .filter(|c| c.status == CheckStatus::Fail)
        .count();
    let skipped = checks
        .iter()
        .filter(|c| c.status == CheckStatus::Skip)
        .count();
    VerifySummary {
        tier,
        ran,
        passed,
        failed,
        skipped,
        checks,
    }
}

fn exit_code(summary: &VerifySummary) -> i32 {
    if summary.failed > 0 {
        1
    } else if summary.ran == 0 {
        2
    } else {
        0
    }
}

fn print_summary(summary: &VerifySummary) {
    println!(
        "lbm verify tier={} ran={} passed={} failed={} skipped={}",
        summary.tier, summary.ran, summary.passed, summary.failed, summary.skipped
    );
    for check in &summary.checks {
        let status = match check.status {
            CheckStatus::Pass => "PASS",
            CheckStatus::Fail => "FAIL",
            CheckStatus::Skip => "SKIP",
        };
        let value = check
            .value
            .map(|v| format!("{v:.6e}"))
            .unwrap_or_else(|| "-".to_string());
        println!(
            "{status}: {} metric={} value={} band={} source={}",
            check.name, check.metric, value, check.band, check.source
        );
        if let Some(message) = &check.message {
            println!("  {message}");
        }
    }
}

fn run_quick_cpu() -> Vec<CheckOutcome> {
    let mut checks = Vec::new();
    checks.push(tgv_cpu_check());
    checks.extend(poiseuille_cpu_checks());
    checks
}

fn tgv_cpu_check() -> CheckOutcome {
    let n = 64;
    let nu = 0.02;
    // Formula, initial density, N=64 setup, and band copied from
    // docs/VALIDATION.md T1 and crates/lbm-core/tests/validation_tgv.rs
    // t1_tgv_trt_accuracy_and_second_order_convergence.
    const BAND: f64 = 1.5e-3;
    let err = tgv_l2_cpu(n, nu);
    pass_fail(
        "T1 Taylor-Green vortex decay (CPU)",
        "quick",
        "L2rel",
        err,
        BAND,
        "<= 1.5e-3",
        "docs/VALIDATION.md T1; crates/lbm-core/tests/validation_tgv.rs",
    )
}

fn tgv_l2_cpu(n: usize, nu: f64) -> f64 {
    let u0 = 1.28 / n as f64;
    let k = 2.0 * PI / n as f64;
    let mut sim: Simulation<f64> = SimConfig {
        nx: n,
        ny: n,
        nu,
        collision: Collision::default(),
        ..Default::default()
    }
    .build()
    .expect("T1 quick check config must build");
    sim.init_with(|x, y| {
        let (xf, yf) = (k * x as f64, k * y as f64);
        let rho = 1.0 - 3.0 * u0 * u0 / 4.0 * ((2.0 * xf).cos() + (2.0 * yf).cos());
        (rho, -u0 * xf.cos() * yf.sin(), u0 * xf.sin() * yf.cos())
    });
    let t_star = (1.0 / (2.0 * nu * k * k)).round() as usize;
    sim.run(t_star);
    let decay = (-2.0 * nu * k * k * t_star as f64).exp();
    let mut actual = Vec::with_capacity(2 * n * n);
    let mut reference = Vec::with_capacity(2 * n * n);
    for y in 0..n {
        for x in 0..n {
            let (xf, yf) = (k * x as f64, k * y as f64);
            actual.push(sim.ux(x, y));
            actual.push(sim.uy(x, y));
            reference.push(-u0 * xf.cos() * yf.sin() * decay);
            reference.push(u0 * xf.sin() * yf.cos() * decay);
        }
    }
    l2_rel(&actual, &reference)
}

fn poiseuille_cpu_checks() -> Vec<CheckOutcome> {
    let ny = 10;
    let nu = 0.1;
    let g = 1.0e-6;
    let mut sim: Simulation<f64> = SimConfig {
        nx: 4,
        ny,
        nu,
        collision: Collision::default(),
        edges: Edges {
            left: EdgeBC::Periodic,
            right: EdgeBC::Periodic,
            bottom: EdgeBC::BounceBack,
            top: EdgeBC::BounceBack,
        },
        force: [g, 0.0],
        ..Default::default()
    }
    .build()
    .expect("T2 quick check config must build");
    // Steady criterion copied from docs/VALIDATION.md notation and
    // crates/lbm-core/tests/validation_channel.rs::poiseuille_horizontal.
    if !run_to_steady(&mut sim, 500, 1.0e-11, 200_000) {
        return vec![failed_without_value(
            "T2 body-force Poiseuille exactness (CPU)",
            "quick",
            "Linf_rel",
            "<= 1.0e-10",
            "docs/VALIDATION.md T2; crates/lbm-core/tests/validation_channel.rs",
            "steady-state criterion was not reached",
        )];
    }
    let h = (ny - 2) as f64;
    let exact: Vec<f64> = (1..=(ny - 2))
        .map(|j| {
            let yw = j as f64 - 0.5;
            g / (2.0 * nu) * yw * (h - yw)
        })
        .collect();
    let got: Vec<f64> = (1..=(ny - 2)).map(|y| sim.ux(0, y)).collect();
    // Bands copied from docs/VALIDATION.md T2 and
    // crates/lbm-core/tests/validation_channel.rs
    // t2_trt_magic_poiseuille_is_exact_and_symmetric.
    let err = linf_rel(&got, &exact);
    let sym = symmetry_abs(&got);
    vec![
        pass_fail(
            "T2 body-force Poiseuille exactness (CPU)",
            "quick",
            "Linf_rel",
            err,
            1.0e-10,
            "<= 1.0e-10",
            "docs/VALIDATION.md T2; crates/lbm-core/tests/validation_channel.rs",
        ),
        pass_fail(
            "T2 body-force Poiseuille top/bottom symmetry (CPU)",
            "quick",
            "max_abs",
            sym,
            1.0e-13,
            "<= 1.0e-13",
            "docs/VALIDATION.md T2; crates/lbm-core/tests/validation_channel.rs",
        ),
    ]
}

fn run_to_steady(
    sim: &mut Simulation<f64>,
    check_every: usize,
    tol: f64,
    max_steps: usize,
) -> bool {
    let mut prev: Vec<f64> = Vec::new();
    let mut elapsed = 0;
    while elapsed < max_steps {
        sim.run(check_every);
        elapsed += check_every;
        let cur: Vec<f64> = sim
            .ux_field()
            .iter()
            .chain(sim.uy_field())
            .copied()
            .collect();
        if !prev.is_empty() {
            let mut dmax = 0.0f64;
            let mut umax = 0.0f64;
            for (c, p) in cur.iter().zip(&prev) {
                dmax = dmax.max((c - p).abs());
                umax = umax.max(c.abs());
            }
            if umax > 0.0 && dmax <= tol * umax {
                return true;
            }
        }
        prev = cur;
    }
    false
}

fn l2_rel(actual: &[f64], reference: &[f64]) -> f64 {
    let mut num = 0.0;
    let mut den = 0.0;
    for (a, r) in actual.iter().zip(reference) {
        num += (a - r) * (a - r);
        den += r * r;
    }
    (num / den).sqrt()
}

fn linf_rel(actual: &[f64], reference: &[f64]) -> f64 {
    let den = reference.iter().map(|v| v.abs()).fold(0.0f64, f64::max);
    actual
        .iter()
        .zip(reference)
        .map(|(a, r)| (a - r).abs())
        .fold(0.0f64, f64::max)
        / den
}

fn symmetry_abs(values: &[f64]) -> f64 {
    let h = values.len();
    let mut max_abs = 0.0f64;
    for j in 0..h / 2 {
        max_abs = max_abs.max((values[j] - values[h - 1 - j]).abs());
    }
    max_abs
}

#[cfg(feature = "gpu")]
fn symmetry_rel(values: &[f64]) -> f64 {
    let h = values.len();
    let mut max_rel = 0.0f64;
    for j in 0..h / 2 {
        let denom = values[j]
            .abs()
            .max(values[h - 1 - j].abs())
            .max(f64::MIN_POSITIVE);
        max_rel = max_rel.max((values[j] - values[h - 1 - j]).abs() / denom);
    }
    max_rel
}

#[cfg(feature = "gpu")]
#[allow(deprecated)]
fn run_gpu_tier() -> Vec<CheckOutcome> {
    use lbm_core::prelude::{
        build_wall_rims, CollisionKind, FaceBC, GlobalSpec, GpuContext, GpuSolver, WallSpec, D2Q9,
    };

    let ctx = match GpuContext::new() {
        Ok(ctx) => ctx,
        Err(e) => {
            return vec![skipped(
                "GPU verification",
                "gpu",
                "adapter",
                "available",
                "docs/VALIDATION.md T14",
                format!("SKIPPED: no usable GPU adapter was found ({e})"),
            )];
        }
    };

    let mut checks = Vec::new();
    checks.push(tgv_gpu_check(&ctx));

    let ny = 10usize;
    let nu = 0.1f64;
    let g = 1.0e-6f32;
    let dims = [4usize, ny, 1usize];
    let mut walls = WallSpec::<f32>::default();
    walls.is_wall[2] = true;
    walls.is_wall[3] = true;
    let (solid, wall_u) = build_wall_rims::<f32>(2, dims, &walls);
    let spec = GlobalSpec::<f32> {
        dims,
        nu,
        collision: CollisionKind::default(),
        periodic: [true, false, true],
        faces: [FaceBC::Closed; 6],
        force: [g, 0.0, 0.0],
        ..Default::default()
    };
    let mut sim =
        match GpuSolver::<D2Q9>::try_new(&spec, &solid, &wall_u, std::sync::Arc::clone(&ctx)) {
            Ok(sim) => sim,
            Err(e) => {
                checks.push(failed_without_value(
                    "T2 body-force Poiseuille exactness (GPU)",
                    "gpu",
                    "Linf_rel",
                    "<= 1.0e-10",
                    "docs/VALIDATION.md T2; crates/lbm-core/tests/validation_channel.rs",
                    format!("GPU solver initialization failed: {e}"),
                ));
                return checks;
            }
        };
    // f32 arithmetic cannot reach the f64 T2 exactness thresholds; use the
    // f32 band class frozen in crates/lbm-core/tests/gpu_absolute.rs
    // (steadiness 1e-7, L_inf_rel <= 1e-5 = T14 field-agreement class).
    if !run_gpu_to_steady(&mut sim, 500, 1.0e-7, 200_000) {
        checks.push(failed_without_value(
            "T2 body-force Poiseuille exactness (GPU, f32 band)",
            "gpu",
            "Linf_rel",
            "<= 1.0e-5",
            "docs/VALIDATION.md T2/T14; crates/lbm-core/tests/gpu_absolute.rs",
            "steady-state criterion was not reached",
        ));
        return checks;
    }
    let h = (ny - 2) as f64;
    let exact: Vec<f64> = (1..=(ny - 2))
        .map(|j| {
            let yw = j as f64 - 0.5;
            g as f64 / (2.0 * nu) * yw * (h - yw)
        })
        .collect();
    let ux = sim.gather_ux();
    let got: Vec<f64> = (1..=(ny - 2)).map(|y| ux[y * dims[0]] as f64).collect();
    checks.push(pass_fail(
        "T2 body-force Poiseuille exactness (GPU, f32 band)",
        "gpu",
        "Linf_rel",
        linf_rel(&got, &exact),
        1.0e-5,
        "<= 1.0e-5",
        "docs/VALIDATION.md T2/T14; crates/lbm-core/tests/gpu_absolute.rs",
    ));
    checks.push(pass_fail(
        "T2 body-force Poiseuille top/bottom symmetry (GPU, f32 band)",
        "gpu",
        "max_rel",
        symmetry_rel(&got),
        1.0e-5,
        "<= 1.0e-5",
        "docs/VALIDATION.md T2/T14; crates/lbm-core/tests/gpu_absolute.rs",
    ));
    checks
}

#[cfg(not(feature = "gpu"))]
fn run_gpu_tier() -> Vec<CheckOutcome> {
    vec![skipped(
        "GPU verification",
        "gpu",
        "feature",
        "compiled",
        "docs/VALIDATION.md T14",
        "SKIPPED: built without --features gpu",
    )]
}

#[cfg(feature = "gpu")]
#[allow(deprecated)]
fn tgv_gpu_check(ctx: &std::sync::Arc<lbm_core::prelude::GpuContext>) -> CheckOutcome {
    use lbm_core::prelude::{CollisionKind, GlobalSpec, GpuSolver, D2Q9};

    let n = 64usize;
    let nu = 0.02f64;
    let u0 = 1.28f32 / n as f32;
    let k = (2.0 * PI / n as f64) as f32;
    let spec = GlobalSpec::<f32> {
        dims: [n, n, 1],
        nu,
        collision: CollisionKind::default(),
        periodic: [true, true, true],
        ..Default::default()
    };
    let solid = vec![false; n * n];
    let wall_u = vec![[0.0f32; 3]; n * n];
    let mut sim =
        match GpuSolver::<D2Q9>::try_new(&spec, &solid, &wall_u, std::sync::Arc::clone(ctx)) {
            Ok(sim) => sim,
            Err(e) => {
                return failed_without_value(
                    "T1 Taylor-Green vortex decay (GPU)",
                    "gpu",
                    "L2rel",
                    "<= 1.5e-3",
                    "docs/VALIDATION.md T1; crates/lbm-core/tests/validation_tgv.rs",
                    format!("GPU solver initialization failed: {e}"),
                )
            }
        };
    sim.init_with(|x, y, _| {
        let (xf, yf) = (k * x as f32, k * y as f32);
        let rho = 1.0 - 3.0 * u0 * u0 / 4.0 * ((2.0 * xf).cos() + (2.0 * yf).cos());
        (
            rho,
            [-u0 * xf.cos() * yf.sin(), u0 * xf.sin() * yf.cos(), 0.0],
        )
    });
    let k64 = 2.0 * PI / n as f64;
    let t_star = (1.0 / (2.0 * nu * k64 * k64)).round() as usize;
    sim.run(t_star);
    sim.sync();
    let ux = sim.gather_ux();
    let uy = sim.gather_uy();
    let decay = (-2.0 * nu * k64 * k64 * t_star as f64).exp();
    let mut actual = Vec::with_capacity(2 * n * n);
    let mut reference = Vec::with_capacity(2 * n * n);
    for y in 0..n {
        for x in 0..n {
            let (xf, yf) = (k64 * x as f64, k64 * y as f64);
            actual.push(ux[y * n + x] as f64);
            actual.push(uy[y * n + x] as f64);
            reference.push(-(u0 as f64) * xf.cos() * yf.sin() * decay);
            reference.push((u0 as f64) * xf.sin() * yf.cos() * decay);
        }
    }
    pass_fail(
        "T1 Taylor-Green vortex decay (GPU)",
        "gpu",
        "L2rel",
        l2_rel(&actual, &reference),
        1.5e-3,
        "<= 1.5e-3",
        "docs/VALIDATION.md T1; crates/lbm-core/tests/validation_tgv.rs",
    )
}

#[cfg(feature = "gpu")]
#[allow(deprecated)]
fn run_gpu_to_steady(
    sim: &mut lbm_core::prelude::GpuSolver<lbm_core::prelude::D2Q9>,
    check_every: usize,
    tol: f64,
    max_steps: usize,
) -> bool {
    let mut prev: Vec<f64> = Vec::new();
    let mut elapsed = 0;
    while elapsed < max_steps {
        sim.run(check_every);
        sim.sync();
        elapsed += check_every;
        let ux = sim.gather_ux();
        let uy = sim.gather_uy();
        let cur: Vec<f64> = ux.iter().chain(&uy).map(|v| *v as f64).collect();
        if !prev.is_empty() {
            let mut dmax = 0.0f64;
            let mut umax = 0.0f64;
            for (c, p) in cur.iter().zip(&prev) {
                dmax = dmax.max((c - p).abs());
                umax = umax.max(c.abs());
            }
            if umax > 0.0 && dmax <= tol * umax {
                return true;
            }
        }
        prev = cur;
    }
    false
}

fn run_mpi_tier() -> Vec<CheckOutcome> {
    if cfg!(feature = "mpi") {
        vec![skipped(
            "MPI verification",
            "mpi",
            "launcher",
            "scripts/test_mpi.sh",
            "docs/VALIDATION.md T13-MPI",
            "SKIPPED: run ./scripts/test_mpi.sh from an MPI-enabled shell",
        )]
    } else {
        vec![skipped(
            "MPI verification",
            "mpi",
            "feature",
            "compiled",
            "docs/VALIDATION.md T13-MPI",
            "UNAVAILABLE: built without --features mpi",
        )]
    }
}

fn pass_fail(
    name: &'static str,
    tier: &'static str,
    metric: &'static str,
    value: f64,
    band_value: f64,
    band: &'static str,
    source: &'static str,
) -> CheckOutcome {
    CheckOutcome {
        name,
        tier,
        status: if value <= band_value {
            CheckStatus::Pass
        } else {
            CheckStatus::Fail
        },
        metric,
        value: Some(value),
        band,
        source,
        message: None,
    }
}

fn failed_without_value(
    name: &'static str,
    tier: &'static str,
    metric: &'static str,
    band: &'static str,
    source: &'static str,
    message: impl Into<String>,
) -> CheckOutcome {
    CheckOutcome {
        name,
        tier,
        status: CheckStatus::Fail,
        metric,
        value: None,
        band,
        source,
        message: Some(message.into()),
    }
}

fn skipped(
    name: &'static str,
    tier: &'static str,
    metric: &'static str,
    band: &'static str,
    source: &'static str,
    message: impl Into<String>,
) -> CheckOutcome {
    CheckOutcome {
        name,
        tier,
        status: CheckStatus::Skip,
        metric,
        value: None,
        band,
        source,
        message: Some(message.into()),
    }
}
