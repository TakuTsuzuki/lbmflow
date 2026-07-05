use lbm_core::compat::multiphase::ShanChen;
use lbm_core::compat::prelude::*;
use lbm_core::prelude::{
    CollisionKind, CpuScalar, GlobalSpec, LocalPeriodic, Solver, WallSpec, D3Q19,
};

fn periodic_sim(nx: usize, ny: usize) -> Simulation<f64> {
    SimConfig {
        nx,
        ny,
        nu: 1.0 / 6.0,
        collision: Collision::Trt {
            magic: CollisionKind::MAGIC_STD,
        },
        ..Default::default()
    }
    .build()
    .unwrap()
}

fn closed_sim(nx: usize, ny: usize) -> Simulation<f64> {
    SimConfig {
        nx,
        ny,
        nu: 1.0 / 6.0,
        collision: Collision::Trt {
            magic: CollisionKind::MAGIC_STD,
        },
        edges: Edges {
            left: EdgeBC::BounceBack,
            right: EdgeBC::BounceBack,
            bottom: EdgeBC::BounceBack,
            top: EdgeBC::BounceBack,
        },
        ..Default::default()
    }
    .build()
    .unwrap()
}

#[test]
fn periodic_uniform_density_gravity_injects_mass_weighted_momentum() {
    let (nx, ny) = (32, 32);
    let g = [1.0e-6, 0.0];
    let mut sim = periodic_sim(nx, ny);
    assert_eq!(sim.gravity(), None);
    sim.set_gravity(g);
    assert_eq!(sim.gravity(), Some(g));

    sim.run(20);
    let p0 = sim.total_momentum()[0];
    let steps = 200usize;
    sim.run(steps);
    let gained = sim.total_momentum()[0] - p0;
    let expect = steps as f64 * sim.fluid_cell_count() as f64 * g[0];
    let rel = ((gained - expect) / expect).abs();
    assert!(
        rel <= 1.0e-10,
        "gravity momentum rel={rel:e}, gained={gained:e}, expected={expect:e}"
    );
}

#[test]
fn native_3d_gravity_is_mass_weighted_and_additive() {
    let dims = [12, 10, 8];
    let n = dims[0] * dims[1] * dims[2];
    let spec = GlobalSpec::<f64> {
        dims,
        nu: 1.0 / 6.0,
        periodic: [true, true, true],
        collision: CollisionKind::Trt {
            magic: CollisionKind::MAGIC_STD,
        },
        ..Default::default()
    };
    let mut solver = Solver::<D3Q19, f64, CpuScalar, LocalPeriodic>::new(
        &spec,
        &vec![false; n],
        &vec![[0.0; 3]; n],
        [1, 1, 1],
        CpuScalar::default(),
        LocalPeriodic,
    );
    solver.set_body_force_field(|_, _, _| [2.0e-7, 0.0, 0.0]);
    solver.set_gravity([0.0, -1.0e-7, 3.0e-7]);

    for _ in 0..20 {
        solver.step();
    }
    let p0 = solver.total_momentum();
    let steps = 200usize;
    for _ in 0..steps {
        solver.step();
    }
    let p1 = solver.total_momentum();
    let mass = n as f64;
    let expect = [
        steps as f64 * mass * 2.0e-7,
        steps as f64 * mass * -1.0e-7,
        steps as f64 * mass * 3.0e-7,
    ];
    for a in 0..3 {
        let gained = p1[a] - p0[a];
        let rel = ((gained - expect[a]) / expect[a]).abs();
        assert!(
            rel <= 1.0e-10,
            "axis={a}, rel={rel:e}, gained={gained:e}, expected={:e}",
            expect[a]
        );
    }
}

#[test]
fn closed_box_gravity_forms_stable_hydrostatic_stratification() {
    let (nx, ny) = (48, 48);
    let mut sim = closed_sim(nx, ny);
    sim.set_gravity([0.0, -1.0e-6]);
    sim.run(20_000);

    let mut max_u = 0.0f64;
    for y in 0..ny {
        for x in 0..nx {
            if sim.is_solid(x, y) {
                continue;
            }
            let u = sim.ux(x, y).hypot(sim.uy(x, y));
            assert!(u.is_finite(), "non-finite velocity at ({x},{y})");
            max_u = max_u.max(u);
        }
    }
    assert!(max_u <= 6.0e-14, "hydrostatic residual max|u|={max_u:e}");

    let q = ny / 4;
    let mut bottom = (0.0, 0usize);
    let mut top = (0.0, 0usize);
    for y in 0..ny {
        for x in 0..nx {
            if sim.is_solid(x, y) {
                continue;
            }
            if y < q {
                bottom.0 += sim.rho(x, y);
                bottom.1 += 1;
            } else if y >= ny - q {
                top.0 += sim.rho(x, y);
                top.1 += 1;
            }
        }
    }
    let rho_bottom = bottom.0 / bottom.1 as f64;
    let rho_top = top.0 / top.1 as f64;
    assert!(
        rho_bottom > rho_top,
        "expected bottom density > top density, got {rho_bottom:e} <= {rho_top:e}"
    );
}

fn phase_center_of_mass_y(sim: &Simulation<f64>, light: bool) -> f64 {
    let mut m = 0.0;
    let mut my = 0.0;
    let rho_cut = 0.5 * (0.15 + 2.0);
    for y in 0..sim.ny() {
        for x in 0..sim.nx() {
            if sim.is_solid(x, y) {
                continue;
            }
            let rho = sim.rho(x, y);
            let w = if light {
                (rho_cut - rho).max(0.0)
            } else {
                (rho - rho_cut).max(0.0)
            };
            if w > 0.0 {
                m += w;
                my += w * y as f64;
            }
        }
    }
    my / m
}

#[test]
fn shan_chen_gravity_composes_with_force_overwrite_and_creates_buoyancy() {
    let (nx, ny) = (80, 80);
    let (cx, cy, r) = (nx / 2, ny / 2, 16.0f64);
    let inside = |x: usize, y: usize| {
        let dx = x as f64 - cx as f64;
        let dy = y as f64 - cy as f64;
        dx * dx + dy * dy <= r * r
    };

    let run_case = |rho_in: f64, rho_out: f64| {
        let mut sim = SimConfig {
            nx,
            ny,
            nu: 1.0 / 6.0,
            collision: Collision::Trt {
                magic: CollisionKind::MAGIC_STD,
            },
            edges: Edges {
                left: EdgeBC::Periodic,
                right: EdgeBC::Periodic,
                bottom: EdgeBC::BounceBack,
                top: EdgeBC::BounceBack,
            },
            ..Default::default()
        }
        .build()
        .unwrap();
        sim.init_with(|x, y| {
            let rho = if inside(x, y) { rho_in } else { rho_out };
            (rho, 0.0, 0.0)
        });
        sim.set_gravity([0.0, -5.0e-6]);
        let sc = ShanChen::new(-5.0);
        let y0 = phase_center_of_mass_y(&sim, rho_in < rho_out);
        for _ in 0..20_000 {
            sc.update_force(&mut sim);
            sim.step();
        }
        let y1 = phase_center_of_mass_y(&sim, rho_in < rho_out);
        (y0, y1)
    };

    let (light0, light1) = run_case(0.15, 2.0);
    assert!(
        light1 - light0 >= 2.0,
        "light blob did not rise enough: {light0:e} -> {light1:e}"
    );
    let (heavy0, heavy1) = run_case(2.0, 0.15);
    assert!(
        heavy0 - heavy1 >= 2.0,
        "heavy blob did not sink enough: {heavy0:e} -> {heavy1:e}"
    );
}

#[test]
fn gravity_skips_solid_cells_in_momentum_accounting() {
    let (nx, ny) = (32, 32);
    let mut sim = periodic_sim(nx, ny);
    sim.set_solid_region(|x, y| (12..20).contains(&x) && (12..20).contains(&y));
    sim.set_gravity([1.0e-6, 0.0]);
    let gained = sim.total_momentum()[0];
    let mass = sim.total_mass();
    let expect = 0.5 * mass * 1.0e-6;
    let rel = ((gained - expect) / expect).abs();
    assert!(
        rel <= 1.0e-8,
        "solid gravity rel={rel:e}, gained={gained:e}, expected={expect:e}"
    );
}

#[test]
fn native_gravity_skips_solid_cells() {
    let dims = [10, 9, 8];
    let mut walls = WallSpec::<f64>::default();
    walls.is_wall = [true; 6];
    let (solid, wall_u) = lbm_core::prelude::build_wall_rims(3, dims, &walls);
    let spec = GlobalSpec::<f64> {
        dims,
        nu: 1.0 / 6.0,
        periodic: [false, false, false],
        collision: CollisionKind::Trt {
            magic: CollisionKind::MAGIC_STD,
        },
        ..Default::default()
    };
    let mut solver = Solver::<D3Q19, f64, CpuScalar, LocalPeriodic>::new(
        &spec,
        &solid,
        &wall_u,
        [1, 1, 1],
        CpuScalar::default(),
        LocalPeriodic,
    );
    solver.set_gravity([0.0, 1.0e-6, 0.0]);
    let gained = solver.total_momentum()[1];
    let expect = solver.total_mass() * 1.0e-6;
    let expect = 0.5 * expect;
    let rel = ((gained - expect) / expect).abs();
    assert!(
        rel <= 1.0e-8,
        "native solid gravity rel={rel:e}, gained={gained:e}, expected={expect:e}"
    );
}
