//! Public per-cell body-force API (`Solver::set_body_force_field`).
//!
//! Three guarantees:
//!   1. momentum accounting — a uniform force in a fully periodic box injects
//!      exactly `F * N_cells` of momentum per step (tested as a rate, so the
//!      constant Guo half-force offset cancels);
//!   2. `clear_body_force_field` actually stops the injection;
//!   3. decomposition invariance — a spatially varying force yields a
//!      bit-identical global field under any `decomp` (the closure is
//!      evaluated at global coordinates and stored in each part's compact
//!      layout), matching the T13 partition-invariance guarantee.

use lbm_core::prelude::*;

fn spec(nx: usize, ny: usize) -> GlobalSpec<f64> {
    GlobalSpec::<f64> {
        dims: [nx, ny, 1],
        nu: 1.0 / 6.0,
        periodic: [true, true, false],
        collision: CollisionKind::Trt { magic: CollisionKind::MAGIC_STD },
        ..Default::default()
    }
}

fn local(nx: usize, ny: usize) -> Solver<D2Q9, f64, CpuScalar, LocalPeriodic> {
    let n = nx * ny;
    Solver::new(&spec(nx, ny), &vec![false; n], &vec![[0.0; 3]; n], [1, 1, 1],
        CpuScalar::default(), LocalPeriodic)
}

fn in_proc(nx: usize, ny: usize, decomp: [usize; 3]) -> Solver<D2Q9, f64, CpuScalar, InProcess> {
    let n = nx * ny;
    Solver::new(&spec(nx, ny), &vec![false; n], &vec![[0.0; 3]; n], decomp,
        CpuScalar::default(), InProcess)
}

#[test]
fn uniform_force_injects_exact_momentum_rate() {
    let (nx, ny) = (16, 16);
    let mut s = local(nx, ny);
    let fx = 1.0e-5;
    s.set_body_force_field(|_, _, _| [fx, 0.0, 0.0]);

    for _ in 0..20 {
        s.step();
    }
    let p1 = s.total_momentum()[0];
    for _ in 0..20 {
        s.step();
    }
    let p2 = s.total_momentum()[0];

    // Per-step injection == F * N_cells (Guo forcing, fully periodic: no loss).
    // Differencing cancels the constant half-force offset and the init state.
    let rate = (p2 - p1) / 20.0;
    let expect = fx * (nx * ny) as f64;
    assert!((rate - expect).abs() / expect < 1e-9, "rate={rate} expected={expect}");
}

#[test]
fn clear_stops_injection() {
    let mut s = local(16, 16);
    s.set_body_force_field(|_, _, _| [1.0e-4, 0.0, 0.0]);
    for _ in 0..10 {
        s.step();
    }
    s.clear_body_force_field();
    for _ in 0..5 {
        s.step();
    }
    let a = s.total_momentum()[0];
    for _ in 0..10 {
        s.step();
    }
    let b = s.total_momentum()[0];
    // No force + periodic + no walls => bulk momentum is frozen step to step.
    assert!((b - a).abs() / a.abs() < 1e-12, "injection continued after clear: {a} -> {b}");
}

#[test]
fn spatially_varying_force_is_decomposition_invariant() {
    let (nx, ny) = (32, 32);
    // A smooth, sign-varying force so an indexing/origin slip would show up.
    let force = |x: usize, y: usize, _z: usize| -> [f64; 3] {
        [1.0e-4 * ((x as f64) * 0.19).sin(), 1.0e-4 * ((y as f64) * 0.23).cos(), 0.0]
    };

    let mut mono = in_proc(nx, ny, [1, 1, 1]);
    let mut part = in_proc(nx, ny, [2, 2, 1]);
    for s in [&mut mono, &mut part] {
        for _ in 0..30 {
            s.set_body_force_field(force);
            s.step();
        }
    }

    let (ax, ay) = (mono.gather_ux(), mono.gather_uy());
    let (bx, by) = (part.gather_ux(), part.gather_uy());
    let max_diff = ax.iter().zip(&bx).chain(ay.iter().zip(&by))
        .map(|(a, b)| (a - b).abs())
        .fold(0.0f64, f64::max);
    // T13-style: decomposition must be bit-exact, not merely close.
    assert_eq!(max_diff, 0.0, "decomposition variance: {max_diff:e}");
}
