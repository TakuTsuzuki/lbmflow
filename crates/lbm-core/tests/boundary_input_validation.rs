use lbm_core::params::MAX_SPEED;
use lbm_core::prelude::*;
use std::panic::{catch_unwind, AssertUnwindSafe};

type Sim = Solver<D2Q9, f64, CpuScalar, LocalPeriodic>;

fn closed_box_spec(dims: [usize; 3]) -> GlobalSpec<f64> {
    GlobalSpec {
        dims,
        nu: 0.08,
        collision: CollisionKind::Bgk,
        periodic: [false; 3],
        faces: [FaceBC::Closed; 6],
        force: [0.0; 3],
        sources: Vec::new(),
        face_patches: Vec::new(),
    }
}

fn periodic_spec(dims: [usize; 3]) -> GlobalSpec<f64> {
    GlobalSpec {
        dims,
        nu: 0.08,
        collision: CollisionKind::Bgk,
        periodic: [true, true, false],
        faces: [FaceBC::Closed; 6],
        force: [0.0; 3],
        sources: Vec::new(),
        face_patches: Vec::new(),
    }
}

fn wall_box_with_top_velocity(u: [f64; 3]) -> (GlobalSpec<f64>, Vec<bool>, Vec<[f64; 3]>) {
    let dims = [6, 6, 1];
    let mut walls = WallSpec::<f64>::default();
    walls.is_wall = [true; 6];
    walls.u[Face::YPos.index()] = u;
    let (solid, wall_u) = build_wall_rims(2, dims, &walls);
    (closed_box_spec(dims), solid, wall_u)
}

fn try_build(
    spec: &GlobalSpec<f64>,
    solid: &[bool],
    wall_u: &[[f64; 3]],
) -> Result<Sim, SpecError> {
    Solver::try_new(
        spec,
        solid,
        wall_u,
        [1, 1, 1],
        CpuScalar::default(),
        LocalPeriodic,
    )
}

fn periodic_solver() -> Sim {
    let spec = periodic_spec([5, 5, 1]);
    Solver::new(
        &spec,
        &[],
        &[],
        [1, 1, 1],
        CpuScalar::default(),
        LocalPeriodic,
    )
}

fn panic_message(err: Box<dyn std::any::Any + Send>) -> String {
    if let Some(msg) = err.downcast_ref::<String>() {
        msg.clone()
    } else if let Some(msg) = err.downcast_ref::<&'static str>() {
        (*msg).to_string()
    } else {
        String::from("<non-string panic>")
    }
}

fn assert_init_panics(init: impl Fn(usize, usize, usize) -> (f64, [f64; 3]), needle: &str) {
    let result = catch_unwind(AssertUnwindSafe(|| {
        let mut sim = periodic_solver();
        sim.init_with(init);
    }));
    let err = result.expect_err("init_with should reject the invalid seed");
    let msg = panic_message(err);
    assert!(
        msg.contains(needle),
        "panic message {msg:?} did not contain {needle:?}"
    );
}

#[test]
fn wall_u_rejects_speeds_above_low_mach_limit() {
    let (spec, solid, wall_u) = wall_box_with_top_velocity([MAX_SPEED, 0.0, 0.0]);
    assert!(
        try_build(&spec, &solid, &wall_u).is_ok(),
        "wall speed exactly at MAX_SPEED must remain legal"
    );

    let (spec, solid, wall_u) = wall_box_with_top_velocity([MAX_SPEED + 1.0e-12, 0.0, 0.0]);
    assert!(
        matches!(
            try_build(&spec, &solid, &wall_u),
            Err(SpecError::VelocityTooHigh { speed }) if speed > MAX_SPEED
        ),
        "wall speed just above MAX_SPEED must return VelocityTooHigh"
    );
}

#[test]
fn wall_u_rejects_non_finite_components() {
    let (spec, solid, mut wall_u) = wall_box_with_top_velocity([0.02, 0.0, 0.0]);
    wall_u[0][0] = f64::NAN;

    assert!(
        matches!(
            try_build(&spec, &solid, &wall_u),
            Err(SpecError::NonFiniteParameter { what: "velocity" })
        ),
        "non-finite wall_u components must be rejected"
    );
}

#[test]
fn init_with_rejects_non_physical_density() {
    assert_init_panics(|_, _, _| (0.0, [0.0; 3]), "density at (0,0,0)");
    assert_init_panics(|_, _, _| (f64::INFINITY, [0.0; 3]), "density at (0,0,0)");
}

#[test]
fn init_with_rejects_non_finite_or_too_fast_velocity() {
    assert_init_panics(
        |_, _, _| (1.0, [f64::NAN, 0.0, 0.0]),
        "velocity component 0 at (0,0,0)",
    );
    assert_init_panics(
        |_, _, _| (1.0, [MAX_SPEED + 1.0e-12, 0.0, 0.0]),
        "exceeds the low-Mach limit",
    );
}

#[test]
fn init_with_accepts_physical_seed_and_limit_boundary() {
    let mut physical = periodic_solver();
    physical.init_with(|_, _, _| (1.0, [0.02, -0.01, 0.0]));

    let mut at_limit = periodic_solver();
    at_limit.init_with(|_, _, _| (1.0, [MAX_SPEED, 0.0, 0.0]));
}
