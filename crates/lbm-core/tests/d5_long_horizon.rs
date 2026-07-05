//! D-5 validation-horizon extension: long-horizon compat/native equivalence
//! plus native Solver absolute TGV convergence.

use lbm_core::backend::CpuScalar;
use lbm_core::compat::prelude as compat;
use lbm_core::halo::LocalPeriodic;
use lbm_core::lattice::{Face, D2Q9};
use lbm_core::prelude::{build_wall_rims, CollisionKind, GlobalSpec, Solver, WallSpec};
use std::f64::consts::PI;

type NativeSolver = Solver<D2Q9, f64, CpuScalar, LocalPeriodic>;

const CAVITY_N: usize = 129;
const CAVITY_U: f64 = 0.1;
const CAVITY_RE: f64 = 100.0;
const CAVITY_STEPS: usize = 20_000;
const MAGIC: f64 = 3.0 / 16.0;

fn compat_cavity() -> compat::Simulation<f64> {
    let l = (CAVITY_N - 2) as f64;
    compat::SimConfig {
        nx: CAVITY_N,
        ny: CAVITY_N,
        nu: CAVITY_U * l / CAVITY_RE,
        collision: compat::Collision::Trt { magic: MAGIC },
        edges: compat::Edges {
            left: compat::EdgeBC::BounceBack,
            right: compat::EdgeBC::BounceBack,
            bottom: compat::EdgeBC::BounceBack,
            top: compat::EdgeBC::MovingWall { u: [CAVITY_U, 0.0] },
        },
        ..Default::default()
    }
    .build()
    .unwrap()
}

fn native_cavity() -> NativeSolver {
    let l = (CAVITY_N - 2) as f64;
    let mut walls = WallSpec::<f64>::default();
    for face in [Face::XNeg, Face::XPos, Face::YNeg, Face::YPos] {
        walls.is_wall[face.index()] = true;
    }
    walls.u[Face::YPos.index()] = [CAVITY_U, 0.0, 0.0];

    let spec = GlobalSpec {
        dims: [CAVITY_N, CAVITY_N, 1],
        nu: CAVITY_U * l / CAVITY_RE,
        collision: CollisionKind::Trt { magic: MAGIC },
        periodic: [false, false, false],
        ..Default::default()
    };
    let (solid, wall_u) = build_wall_rims(2, spec.dims, &walls);
    Solver::new(
        &spec,
        &solid,
        &wall_u,
        [1, 1, 1],
        CpuScalar::default(),
        LocalPeriodic,
    )
}

fn max_mixed_relative(a: &[f64], b: &[f64]) -> f64 {
    assert_eq!(a.len(), b.len());
    a.iter()
        .zip(b)
        .map(|(x, y)| {
            let scale = x.abs().max(y.abs()).max(1.0);
            (x - y).abs() / scale
        })
        .fold(0.0, f64::max)
}

#[test]
#[ignore]
fn d5_cavity_re100_compat_matches_native_after_20k_steps() {
    let mut facade = compat_cavity();
    let mut native = native_cavity();

    facade.run(CAVITY_STEPS);
    native.run(CAVITY_STEPS);

    let rho = max_mixed_relative(facade.rho_field(), &native.gather_rho());
    let ux = max_mixed_relative(facade.ux_field(), &native.gather_ux());
    let uy = max_mixed_relative(facade.uy_field(), &native.gather_uy());
    let worst = rho.max(ux).max(uy);

    println!(
        "D-5 cavity Re={CAVITY_RE:.0} {CAVITY_N}x{CAVITY_N} after {CAVITY_STEPS} steps: \
         rho={rho:.6e}, ux={ux:.6e}, uy={uy:.6e}, worst={worst:.6e}"
    );
    assert!(
        worst <= 1.0e-12,
        "D-5 long-horizon compat/native max mixed relative field deviation = {worst:e} \
         (rho={rho:e}, ux={ux:e}, uy={uy:e})"
    );
}

fn native_tgv_l2(n: usize) -> f64 {
    let nu = 0.02;
    let u0 = 1.28 / n as f64;
    let k = 2.0 * PI / n as f64;
    let spec = GlobalSpec {
        dims: [n, n, 1],
        nu,
        collision: CollisionKind::Trt { magic: MAGIC },
        periodic: [true, true, false],
        ..Default::default()
    };
    let mut solver: NativeSolver = Solver::new(
        &spec,
        &[],
        &[],
        [1, 1, 1],
        CpuScalar::default(),
        LocalPeriodic,
    );
    solver.init_with(|x, y, _| {
        let (xf, yf) = (k * x as f64, k * y as f64);
        let rho = 1.0 - 3.0 * u0 * u0 / 4.0 * ((2.0 * xf).cos() + (2.0 * yf).cos());
        (
            rho,
            [-u0 * xf.cos() * yf.sin(), u0 * xf.sin() * yf.cos(), 0.0],
        )
    });

    let steps = (1.0 / (2.0 * nu * k * k)).round() as usize;
    solver.run(steps);
    let ux = solver.gather_ux();
    let uy = solver.gather_uy();
    let decay = (-2.0 * nu * k * k * steps as f64).exp();

    let mut num = 0.0;
    let mut den = 0.0;
    for y in 0..n {
        for x in 0..n {
            let i = y * n + x;
            let (xf, yf) = (k * x as f64, k * y as f64);
            let ux_ref = -u0 * xf.cos() * yf.sin() * decay;
            let uy_ref = u0 * xf.sin() * yf.cos() * decay;
            num += (ux[i] - ux_ref).powi(2) + (uy[i] - uy_ref).powi(2);
            den += ux_ref.powi(2) + uy_ref.powi(2);
        }
    }
    (num / den).sqrt()
}

#[test]
fn d5_native_solver_tgv_converges_second_order() {
    let e32 = native_tgv_l2(32);
    let e64 = native_tgv_l2(64);
    let order = (e32 / e64).log2();
    println!("D-5 native TGV convergence: e32={e32:.6e}, e64={e64:.6e}, order={order:.3}");
    assert!(
        order >= 1.7,
        "D-5 native Solver TGV order = {order:e}, e32 = {e32:e}, e64 = {e64:e}"
    );
}
