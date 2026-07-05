//! A-9: `run_guarded` — the standard run-time non-finite watchdog.
//!
//! The physics kernels are deliberately guard-free (V1 equivalence), so
//! divergence detection lives in the driver: `run_guarded(steps,
//! check_every)` periodically inspects the f64 mass aggregation
//! (`local_mass_partials`), into which any NaN/±Inf population propagates.
//!
//! Acceptance (SOLVER_IMPROVEMENT_SPEC A-9): a 1-cell NaN injection is
//! detected within one check interval, with the step number; a healthy run
//! is bit-identical to `run`; overhead < 1% at 512².

use lbm_core::lattice::D2Q9;
use lbm_core::prelude::*;

type Cpu = Solver<D2Q9, f64, CpuScalar, LocalPeriodic>;

fn tgv(n: usize) -> Cpu {
    let spec = GlobalSpec::<f64> {
        dims: [n, n, 1],
        nu: 0.02,
        periodic: [true, true, false],
        ..Default::default()
    };
    let mut s: Cpu = Solver::new(
        &spec,
        &[],
        &[],
        [1, 1, 1],
        CpuScalar::default(),
        LocalPeriodic,
    );
    let k = 2.0 * std::f64::consts::PI / n as f64;
    s.init_with(move |x, y, _| {
        let (xx, yy) = (k * x as f64, k * y as f64);
        (
            1.0,
            [0.03 * xx.sin() * yy.cos(), -0.03 * xx.cos() * yy.sin(), 0.0],
        )
    });
    s
}

/// A single NaN population injected mid-run is caught at the next check,
/// with the completed-step count in the error.
#[test]
fn nan_injection_detected_within_one_interval() {
    let mut s = tgv(32);
    s.run(5); // healthy prefix
    assert_eq!(s.time(), 5);

    // Inject NaN into one population of one interior fluid cell.
    {
        let fields = s.fields_mut(0);
        let g = fields.geom;
        let np = g.n_padded();
        let pi = g.pidx(7, 9, 0);
        fields.f[3 * np + pi] = f64::NAN;
    }

    let check_every = 4;
    let err = s
        .run_guarded(50, check_every)
        .expect_err("NaN must be detected");
    // First check runs at step 5 + 4 = 9; detection must not be later.
    assert_eq!(err.step, 9, "detected at the first check after injection");
    assert_eq!(s.time(), 9, "run_guarded stops at the failing check");
    // The error carries a readable message with the step.
    let msg = format!("{err}");
    assert!(msg.contains("step 9"), "{msg}");
}

/// An initial-state NaN is caught even when `steps` is smaller than
/// `check_every` (the tail check).
#[test]
fn tail_check_catches_short_runs() {
    let mut s = tgv(16);
    {
        let fields = s.fields_mut(0);
        let g = fields.geom;
        let np = g.n_padded();
        let pi = g.pidx(3, 3, 0);
        fields.f[1 * np + pi] = f64::INFINITY; // ±Inf must be caught too
    }
    let err = s.run_guarded(3, 100).expect_err("Inf must be detected");
    assert_eq!(err.step, 3);
}

/// A healthy run returns Ok and produces a trajectory bit-identical to
/// `run` — the watchdog only reads, never writes (the legal-configuration
/// bit-invariance DoD).
#[test]
fn healthy_run_is_ok_and_bit_identical() {
    let mut a = tgv(32);
    let mut b = tgv(32);
    a.run(50);
    b.run_guarded(50, 7).expect("healthy run must be Ok");
    assert_eq!(a.time(), b.time());
    let (fa, fb) = (a.fields(0), b.fields(0));
    assert_eq!(fa.f, fb.f, "populations must be bit-identical");
    assert_eq!(fa.rho, fb.rho);
    assert_eq!(fa.ux, fb.ux);
    assert_eq!(fa.uy, fb.uy);
}

/// Overhead of the watchdog at 512² with check_every = 100 (the A-9
/// acceptance line: < 1%). Component-timed for robustness — the per-check
/// reduction cost is measured directly against the per-step cost, then an
/// end-to-end comparison is printed for the record. Heavy: run with
/// `--include-ignored`.
#[test]
#[ignore]
fn overhead_under_one_percent_at_512sq() {
    use std::time::Instant;
    let n = 512;
    let check_every = 100;

    // Warm-up + per-step cost.
    let mut s = tgv(n);
    s.run(20);
    let steps = 200;
    let t0 = Instant::now();
    s.run(steps);
    let per_step = t0.elapsed().as_secs_f64() / steps as f64;

    // Per-check cost (the mass reduction the watchdog adds).
    let reps = 50;
    let t1 = Instant::now();
    let mut acc = 0.0f64;
    for _ in 0..reps {
        let (fluid, m) = s.local_mass_partials();
        acc += fluid + m;
    }
    let per_check = t1.elapsed().as_secs_f64() / reps as f64;
    assert!(acc.is_finite());

    let overhead = per_check / (check_every as f64 * per_step);
    println!(
        "run_guarded overhead @512^2: per_step {:.3e}s, per_check {:.3e}s, \
         check_every {check_every} -> overhead {:.4}%",
        per_step,
        per_check,
        overhead * 100.0
    );

    // End-to-end record (informational; the assert is on the component ratio,
    // which is robust against scheduler noise).
    let mut plain = tgv(n);
    let mut guarded = tgv(n);
    plain.run(20);
    guarded.run(20);
    let t2 = Instant::now();
    plain.run(steps);
    let t_plain = t2.elapsed().as_secs_f64();
    let t3 = Instant::now();
    guarded.run_guarded(steps, check_every).unwrap();
    let t_guarded = t3.elapsed().as_secs_f64();
    println!(
        "run_guarded end-to-end @512^2 x{steps}: run {t_plain:.3}s, \
         run_guarded {t_guarded:.3}s ({:+.3}%)",
        (t_guarded / t_plain - 1.0) * 100.0
    );

    assert!(
        overhead < 0.01,
        "watchdog overhead {:.4}% exceeds the 1% acceptance line",
        overhead * 100.0
    );
}

// ---------------------------------------------------------------------------
// GPU counterpart (feature `gpu`): same API through the readback path.
// ---------------------------------------------------------------------------
#[cfg(feature = "gpu")]
mod gpu {
    use super::*;
    use std::sync::{Arc, OnceLock};

    fn ctx() -> Arc<GpuContext> {
        static CTX: OnceLock<Arc<GpuContext>> = OnceLock::new();
        CTX.get_or_init(|| GpuContext::new().expect("run_guarded GPU test requires a GPU adapter"))
            .clone()
    }

    /// A NaN seeded into the initial state is detected by the readback check
    /// within one interval, with the step number; a healthy GPU run is Ok.
    #[test]
    fn gpu_run_guarded_detects_nan_and_passes_healthy() {
        let n = 32usize;
        let spec = GlobalSpec::<f32> {
            dims: [n, n, 1],
            nu: 0.02,
            periodic: [true, true, false],
            ..Default::default()
        };

        // Healthy: Ok and time advances.
        let mut healthy: GpuSolver<D2Q9> = GpuSolver::new(&spec, &[], &[], ctx());
        let k = 2.0 * std::f32::consts::PI / n as f32;
        healthy.init_with(move |x, y, _| {
            let (xx, yy) = (k * x as f32, k * y as f32);
            (
                1.0,
                [0.03 * xx.sin() * yy.cos(), -0.03 * xx.cos() * yy.sin(), 0.0],
            )
        });
        healthy
            .run_guarded(60, 25)
            .expect("healthy GPU run must be Ok");
        assert_eq!(healthy.time(), 60);

        // Poisoned: one cell's density is NaN from the start.
        let mut bad: GpuSolver<D2Q9> = GpuSolver::new(&spec, &[], &[], ctx());
        bad.init_with(|x, y, _| {
            if (x, y) == (5, 5) {
                (f32::NAN, [0.0; 3])
            } else {
                (1.0, [0.0; 3])
            }
        });
        let err = bad.run_guarded(60, 25).expect_err("NaN must be detected");
        assert_eq!(err.step, 25, "caught at the first readback check");
    }
}
