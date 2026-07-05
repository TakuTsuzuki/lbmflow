//! Stage 2 gate: the V2 D2Q9 CpuScalar specialisation must reproduce V1
//! (`lbm-core`) trajectories at bit-exact level (`f64` max |Δ| ≤ 1e-14; in
//! practice the port is operand-order faithful, so the observed difference
//! is exactly 0.0 — asserted).
//!
//! Every configuration steps V1 and V2 side by side and compares the full
//! rho/ux/uy fields plus the f64 diagnostics (total mass / momentum /
//! probed force) at every step for the first ten steps and periodically
//! afterwards.

use lbm_core::prelude::{
    Collision as V1Collision, EdgeBC as V1EdgeBC, Edges as V1Edges, SimConfig as V1SimConfig,
    Simulation as V1Simulation,
};
use lbm_core2::lattice::D2Q9;
use lbm_core2::prelude::*;

type V2Solver<T> = Solver<D2Q9, T, CpuScalar, LocalPeriodic>;

/// Map a V1 config onto the V2 core spec + wall rims (what the compat facade
/// does; duplicated here so this test does not depend on Stage 3).
fn v2_from_v1_edges<T: Real + lbm_core::real::Real>(
    nx: usize,
    ny: usize,
    nu: f64,
    collision: V1Collision,
    edges: &V1Edges<T>,
    force: [T; 2],
) -> V2Solver<T> {
    let faces4 = [
        (Face::XNeg, edges.left),
        (Face::XPos, edges.right),
        (Face::YNeg, edges.bottom),
        (Face::YPos, edges.top),
    ];
    let mut walls = WallSpec::<T>::default();
    let mut faces = [FaceBC::Closed; 6];
    let mut periodic = [false, false, false];
    for (face, bc) in faces4 {
        match bc {
            V1EdgeBC::Periodic => periodic[face.axis()] = true,
            V1EdgeBC::BounceBack => walls.is_wall[face.index()] = true,
            V1EdgeBC::MovingWall { u } => {
                walls.is_wall[face.index()] = true;
                walls.u[face.index()] = [u[0], u[1], T::zero()];
            }
            V1EdgeBC::VelocityInlet { u } => {
                faces[face.index()] = FaceBC::Velocity {
                    u: [u[0], u[1], T::zero()],
                }
            }
            V1EdgeBC::PressureOutlet { rho } => faces[face.index()] = FaceBC::Pressure { rho },
            V1EdgeBC::Outflow => faces[face.index()] = FaceBC::Outflow,
            V1EdgeBC::ConvectiveOutflow { u_conv } => {
                faces[face.index()] = FaceBC::Convective { u_conv }
            }
        }
    }
    let spec = GlobalSpec {
        dims: [nx, ny, 1],
        nu,
        collision: match collision {
            V1Collision::Bgk => CollisionKind::Bgk,
            V1Collision::Trt { magic } => CollisionKind::Trt { magic },
        },
        periodic,
        faces,
        force: [force[0], force[1], T::zero()],
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

fn v2_from_v1_config(cfg: &V1SimConfig<f64>) -> V2Solver<f64> {
    v2_from_v1_edges(cfg.nx, cfg.ny, cfg.nu, cfg.collision, &cfg.edges, cfg.force)
}

fn max_abs_diff(a: &[f64], b: &[f64]) -> f64 {
    assert_eq!(a.len(), b.len());
    a.iter()
        .zip(b)
        .map(|(x, y)| (x - y).abs())
        .fold(0.0, f64::max)
}

/// Compare fields + diagnostics; returns the max difference seen.
fn compare(v1: &V1Simulation<f64>, v2: &V2Solver<f64>, tol: f64, what: &str) -> f64 {
    let mut worst = 0.0f64;
    let pairs: [(&str, &[f64], Vec<f64>); 3] = [
        ("rho", v1.rho_field(), v2.gather_rho()),
        ("ux", v1.ux_field(), v2.gather_ux()),
        ("uy", v1.uy_field(), v2.gather_uy()),
    ];
    for (name, a, b) in pairs {
        let d = max_abs_diff(a, &b);
        assert!(d <= tol, "{what}: {name} max|Δ| = {d:e} > {tol:e}");
        worst = worst.max(d);
    }
    let dm = (v1.total_mass() - v2.total_mass()).abs();
    assert!(dm <= tol, "{what}: total_mass Δ = {dm:e}");
    let p1 = v1.total_momentum();
    let p2 = v2.total_momentum();
    for a in 0..2 {
        let dp = (p1[a] - p2[a]).abs();
        assert!(dp <= tol, "{what}: total_momentum[{a}] Δ = {dp:e}");
        worst = worst.max(dp);
    }
    let f1 = v1.probed_force();
    let f2 = v2.probed_force();
    for a in 0..2 {
        let df = (f1[a] - f2[a]).abs();
        assert!(df <= tol, "{what}: probed_force[{a}] Δ = {df:e}");
        worst = worst.max(df);
    }
    worst
}

/// Step both engines `steps` times, comparing every step for the first 10
/// and every 25 afterwards. Returns the worst difference observed.
fn run_and_compare(
    v1: &mut V1Simulation<f64>,
    v2: &mut V2Solver<f64>,
    steps: usize,
    tol: f64,
    what: &str,
) -> f64 {
    let mut worst = compare(v1, v2, tol, &format!("{what} t=0"));
    for s in 1..=steps {
        v1.step();
        v2.step();
        if s <= 10 || s % 25 == 0 || s == steps {
            worst = worst.max(compare(v1, v2, tol, &format!("{what} t={s}")));
        }
    }
    println!("{what}: worst |Δ| over {steps} steps = {worst:e}");
    worst
}

const TOL: f64 = 1e-14;

// ---------------------------------------------------------------------------

#[test]
fn tgv_periodic_trt_matches_v1_bitwise() {
    let cfg = V1SimConfig::<f64> {
        nx: 64,
        ny: 48,
        nu: 0.02,
        ..Default::default()
    };
    let mut v2 = v2_from_v1_config(&cfg);
    let mut v1 = cfg.build().unwrap();
    let (nx, ny) = (64usize, 48usize);
    let init = move |x: usize, y: usize| {
        let kx = 2.0 * std::f64::consts::PI / nx as f64;
        let ky = 2.0 * std::f64::consts::PI / ny as f64;
        let u0 = 0.04;
        let ux = -u0 * (kx * x as f64).cos() * (ky * y as f64).sin();
        let uy = u0 * (kx * x as f64).sin() * (ky * y as f64).cos();
        let rho = 1.0 + 0.01 * (kx * x as f64).cos() * (ky * y as f64).cos();
        (rho, ux, uy)
    };
    v1.init_with(init);
    v2.init_with(move |x, y, _| {
        let (r, ux, uy) = init(x, y);
        (r, [ux, uy, 0.0])
    });
    let worst = run_and_compare(&mut v1, &mut v2, 500, TOL, "tgv f64");
    // The port is operand-order exact: expect literal zero, not just <=1e-14.
    assert_eq!(worst, 0.0, "expected bit-exact TGV trajectory");
}

#[test]
fn tgv_periodic_bgk_f32_matches_v1_bitwise() {
    // Identical f32 operations round identically: bit-exact in f32 too.
    let cfg = V1SimConfig::<f32> {
        nx: 48,
        ny: 32,
        nu: 0.05,
        collision: V1Collision::Bgk,
        ..Default::default()
    };
    let mut v2 = v2_from_v1_edges::<f32>(48, 32, 0.05, cfg.collision, &cfg.edges, [0.0, 0.0]);
    let mut v1 = cfg.build().unwrap();
    let init = |x: usize, y: usize| {
        let kx = 2.0 * std::f64::consts::PI / 48.0;
        let ky = 2.0 * std::f64::consts::PI / 32.0;
        let ux = -0.03 * (kx * x as f64).cos() * (ky * y as f64).sin();
        let uy = 0.03 * (kx * x as f64).sin() * (ky * y as f64).cos();
        (1.0f32, ux as f32, uy as f32)
    };
    v1.init_with(init);
    v2.init_with(move |x, y, _| {
        let (r, ux, uy) = init(x, y);
        (r, [ux, uy, 0.0])
    });
    for s in 1..=300 {
        v1.step();
        v2.step();
        if s <= 5 || s % 50 == 0 {
            let ux2 = v2.gather_ux();
            let uy2 = v2.gather_uy();
            let rho2 = v2.gather_rho();
            for i in 0..48 * 32 {
                assert_eq!(v1.ux_field()[i], ux2[i], "f32 ux t={s} i={i}");
                assert_eq!(v1.uy_field()[i], uy2[i], "f32 uy t={s} i={i}");
                assert_eq!(v1.rho_field()[i], rho2[i], "f32 rho t={s} i={i}");
            }
        }
    }
    println!("tgv f32/bgk: bit-exact over 300 steps");
}

#[test]
fn lid_driven_cavity_matches_v1_bitwise() {
    let cfg = V1SimConfig::<f64> {
        nx: 48,
        ny: 48,
        nu: 0.02,
        edges: V1Edges {
            left: V1EdgeBC::BounceBack,
            right: V1EdgeBC::BounceBack,
            bottom: V1EdgeBC::BounceBack,
            top: V1EdgeBC::MovingWall { u: [0.1, 0.0] },
        },
        ..Default::default()
    };
    let mut v2 = v2_from_v1_config(&cfg);
    let mut v1 = cfg.build().unwrap();
    // Sanity: identical rims.
    for y in 0..48 {
        for x in 0..48 {
            assert_eq!(v1.is_solid(x, y), v2.is_solid(x, y, 0), "solid ({x},{y})");
        }
    }
    let worst = run_and_compare(&mut v1, &mut v2, 400, TOL, "cavity f64");
    assert_eq!(worst, 0.0, "expected bit-exact cavity trajectory");
}

#[test]
fn poiseuille_body_force_matches_v1_bitwise() {
    // Periodic x, walls y, uniform Guo force: covers the forcing path.
    let cfg = V1SimConfig::<f64> {
        nx: 32,
        ny: 24,
        nu: 0.1,
        edges: V1Edges {
            left: V1EdgeBC::Periodic,
            right: V1EdgeBC::Periodic,
            bottom: V1EdgeBC::BounceBack,
            top: V1EdgeBC::BounceBack,
        },
        force: [1e-5, 0.0],
        ..Default::default()
    };
    let mut v2 = v2_from_v1_config(&cfg);
    let mut v1 = cfg.build().unwrap();
    let worst = run_and_compare(&mut v1, &mut v2, 400, TOL, "poiseuille f64");
    assert_eq!(worst, 0.0, "expected bit-exact poiseuille trajectory");
}

#[test]
fn channel_inlet_profile_outflow_matches_v1_bitwise() {
    // Zou-He velocity inlet with a parabolic per-node profile + zero-gradient
    // outflow + wall rims.
    let cfg = V1SimConfig::<f64> {
        nx: 80,
        ny: 40,
        nu: 0.05,
        edges: V1Edges {
            left: V1EdgeBC::VelocityInlet { u: [0.05, 0.0] },
            right: V1EdgeBC::Outflow,
            bottom: V1EdgeBC::BounceBack,
            top: V1EdgeBC::BounceBack,
        },
        ..Default::default()
    };
    let mut v2 = v2_from_v1_config(&cfg);
    let mut v1 = cfg.build().unwrap();
    let ny = 40usize;
    let prof = |y: usize| {
        let yy = y as f64 / (ny - 1) as f64;
        [0.08 * 4.0 * yy * (1.0 - yy), 0.0]
    };
    v1.set_inlet_profile(lbm_core::prelude::Edge::Left, prof);
    let values: Vec<[f64; 3]> = (0..ny).map(|y| [prof(y)[0], prof(y)[1], 0.0]).collect();
    v2.set_inlet_profile(Face::XNeg, &values);
    let worst = run_and_compare(&mut v1, &mut v2, 400, TOL, "channel profile f64");
    assert_eq!(worst, 0.0, "expected bit-exact channel trajectory");
}

#[test]
fn pressure_outlet_matches_v1_bitwise() {
    // Uniform velocity inlet + Zou-He pressure outlet.
    let cfg = V1SimConfig::<f64> {
        nx: 64,
        ny: 32,
        nu: 0.05,
        edges: V1Edges {
            left: V1EdgeBC::VelocityInlet { u: [0.04, 0.0] },
            right: V1EdgeBC::PressureOutlet { rho: 1.0 },
            bottom: V1EdgeBC::BounceBack,
            top: V1EdgeBC::BounceBack,
        },
        ..Default::default()
    };
    let mut v2 = v2_from_v1_config(&cfg);
    let mut v1 = cfg.build().unwrap();
    let worst = run_and_compare(&mut v1, &mut v2, 400, TOL, "pressure outlet f64");
    assert_eq!(worst, 0.0, "expected bit-exact pressure-outlet trajectory");
}

#[test]
fn convective_outflow_matches_v1_bitwise() {
    // Convective outflow depends on the ping-pong stale-slot mechanics; this
    // is the regression that pins that implicit V1 contract.
    let cfg = V1SimConfig::<f64> {
        nx: 80,
        ny: 40,
        nu: 0.03,
        edges: V1Edges {
            left: V1EdgeBC::VelocityInlet { u: [0.06, 0.0] },
            right: V1EdgeBC::ConvectiveOutflow { u_conv: 0.06 },
            bottom: V1EdgeBC::BounceBack,
            top: V1EdgeBC::BounceBack,
        },
        ..Default::default()
    };
    let mut v2 = v2_from_v1_config(&cfg);
    let mut v1 = cfg.build().unwrap();
    let worst = run_and_compare(&mut v1, &mut v2, 400, TOL, "convective f64");
    assert_eq!(worst, 0.0, "expected bit-exact convective-outflow trajectory");
}

#[test]
fn cylinder_probe_matches_v1_bitwise() {
    // Momentum-exchange force probe on an obstacle (drag/lift path).
    let (nx, ny) = (96usize, 48usize);
    let cfg = V1SimConfig::<f64> {
        nx,
        ny,
        nu: 0.02,
        edges: V1Edges {
            left: V1EdgeBC::VelocityInlet { u: [0.05, 0.0] },
            right: V1EdgeBC::Outflow,
            bottom: V1EdgeBC::BounceBack,
            top: V1EdgeBC::BounceBack,
        },
        ..Default::default()
    };
    let mut v2 = v2_from_v1_config(&cfg);
    let mut v1 = cfg.build().unwrap();
    let (cx, cy, r) = (24.0, 23.7, 6.3);
    let inside = move |x: usize, y: usize| {
        let (dx, dy) = (x as f64 - cx, y as f64 - cy);
        dx * dx + dy * dy < r * r
    };
    v1.set_solid_region(inside);
    for y in 0..ny {
        for x in 0..nx {
            if inside(x, y) {
                v2.set_solid(x, y, 0);
            }
        }
    }
    v1.set_force_probe(inside);
    v2.set_force_probe(move |x, y, _| inside(x, y));
    let mut worst = 0.0f64;
    for s in 1..=300 {
        v1.step();
        v2.step();
        // Probe force compared every step (it is the sensitive diagnostic).
        let (f1, f2) = (v1.probed_force(), v2.probed_force());
        for a in 0..2 {
            let d = (f1[a] - f2[a]).abs();
            assert!(d <= TOL, "t={s}: probed_force[{a}] Δ = {d:e}");
            worst = worst.max(d);
        }
        if s % 50 == 0 || s == 300 {
            worst = worst.max(compare(&v1, &v2, TOL, &format!("cylinder t={s}")));
        }
    }
    println!("cylinder: worst |Δ| = {worst:e}");
    assert_eq!(worst, 0.0, "expected bit-exact cylinder trajectory");
}

#[test]
fn per_cell_force_field_matches_v1_bitwise() {
    // The per-cell force path is the mechanism multiphase models drive;
    // exercise it with a synthetic rotating force pattern rewritten each
    // step, on top of a uniform force.
    let cfg = V1SimConfig::<f64> {
        nx: 48,
        ny: 32,
        nu: 0.04,
        force: [2e-6, -1e-6],
        ..Default::default()
    };
    let mut v2 = v2_from_v1_config(&cfg);
    let mut v1 = cfg.build().unwrap();
    let (nx, ny) = (48usize, 32usize);
    let pat = |x: usize, y: usize, t: usize| {
        let kx = 2.0 * std::f64::consts::PI / nx as f64;
        let ky = 2.0 * std::f64::consts::PI / ny as f64;
        let ph = t as f64 * 0.01;
        [
            1e-5 * ((kx * x as f64) + ph).sin() * (ky * y as f64).cos(),
            1e-5 * (kx * x as f64).cos() * ((ky * y as f64) - ph).sin(),
        ]
    };
    let mut worst = 0.0f64;
    for s in 0..250 {
        {
            let ff1 = v1.force_field_mut();
            for y in 0..ny {
                for x in 0..nx {
                    ff1[y * nx + x] = pat(x, y, s);
                }
            }
        }
        {
            let f2 = v2.fields_mut(0);
            let ff2 = f2
                .force_field
                .get_or_insert_with(|| vec![[0.0; 3]; nx * ny]);
            for y in 0..ny {
                for x in 0..nx {
                    let p = pat(x, y, s);
                    ff2[y * nx + x] = [p[0], p[1], 0.0];
                }
            }
        }
        v1.step();
        v2.step();
        if s % 25 == 0 || s == 249 {
            worst = worst.max(compare(&v1, &v2, TOL, &format!("force field t={s}")));
        }
    }
    println!("force field: worst |Δ| = {worst:e}");
    assert_eq!(worst, 0.0, "expected bit-exact force-field trajectory");
}

#[test]
fn build_time_state_matches_v1() {
    // V1 implicit spec: from_config ends with update_moments, so u(t=0)
    // includes the half-force correction even before any step.
    let cfg = V1SimConfig::<f64> {
        nx: 16,
        ny: 12,
        nu: 0.1,
        force: [4e-4, -2e-4],
        edges: V1Edges {
            left: V1EdgeBC::Periodic,
            right: V1EdgeBC::Periodic,
            bottom: V1EdgeBC::BounceBack,
            top: V1EdgeBC::BounceBack,
        },
        ..Default::default()
    };
    let v2 = v2_from_v1_config(&cfg);
    let v1 = cfg.build().unwrap();
    assert_eq!(compare(&v1, &v2, TOL, "build-time"), 0.0);
    // Fluid cells carry u = F/2 at t = 0 (documented V1 behaviour).
    assert_eq!(v1.ux(4, 4), 2e-4);
    assert_eq!(v2.u(4, 4, 0)[0], 2e-4);
}
