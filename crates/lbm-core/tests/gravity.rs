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

fn closed_box_quiescence_2d<T: lbm_core::prelude::Real>(name: &str, bound: f64) -> f64 {
    let (nx, ny) = (48, 48);
    let mut sim: Simulation<T> = SimConfig {
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
    .unwrap();
    sim.set_gravity([T::zero(), T::r(-1.0e-6)]);
    sim.run(5_000);

    let mut max_u = 0.0f64;
    for y in 0..ny {
        for x in 0..nx {
            if sim.is_solid(x, y) {
                continue;
            }
            let ux = sim.ux(x, y).as_f64();
            let uy = sim.uy(x, y).as_f64();
            assert!(
                ux.is_finite() && uy.is_finite(),
                "{name}: non-finite velocity at ({x},{y})"
            );
            max_u = max_u.max(ux.hypot(uy));
        }
    }
    assert!(
        max_u <= bound,
        "{name}: VR-STR-06 residual max|u|={max_u:e}, bound={bound:e}"
    );
    max_u
}

fn closed_box_quiescence_3d<T: lbm_core::prelude::Real>(name: &str, bound: f64) -> f64 {
    let dims = [24, 24, 16];
    let mut walls = WallSpec::<T>::default();
    walls.is_wall = [true; 6];
    let (solid, wall_u) = lbm_core::prelude::build_wall_rims(3, dims, &walls);
    let spec = GlobalSpec::<T> {
        dims,
        nu: 1.0 / 6.0,
        collision: CollisionKind::Trt {
            magic: CollisionKind::MAGIC_STD,
        },
        periodic: [false, false, false],
        ..Default::default()
    };
    let mut solver = Solver::<D3Q19, T, CpuScalar, LocalPeriodic>::new(
        &spec,
        &solid,
        &wall_u,
        [1, 1, 1],
        CpuScalar::default(),
        LocalPeriodic,
    );
    solver.set_gravity([T::r(3.0e-7), T::r(-8.0e-7), T::r(2.0e-7)]);
    solver.run(5_000);

    let ux = solver.gather_ux();
    let uy = solver.gather_uy();
    let uz = solver.gather_uz();
    let mut max_u = 0.0f64;
    for z in 0..dims[2] {
        for y in 0..dims[1] {
            for x in 0..dims[0] {
                if solver.is_solid(x, y, z) {
                    continue;
                }
                let i = (z * dims[1] + y) * dims[0] + x;
                let u = ux[i].as_f64().hypot(uy[i].as_f64()).hypot(uz[i].as_f64());
                assert!(
                    u.is_finite(),
                    "{name}: non-finite velocity at ({x},{y},{z})"
                );
                max_u = max_u.max(u);
            }
        }
    }
    assert!(
        max_u <= bound,
        "{name}: VR-STR-06 residual max|u|={max_u:e}, bound={bound:e}"
    );
    max_u
}

#[test]
fn vr_str_06_static_stratification_quiescent_all_lattices_and_precisions() {
    let d2_f64 = closed_box_quiescence_2d::<f64>("D2Q9/f64", 2.0e-9);
    let d2_f32 = closed_box_quiescence_2d::<f32>("D2Q9/f32", 5.0e-7);
    let d3_f64 = closed_box_quiescence_3d::<f64>("D3Q19/f64", 1.0e-13);
    let d3_f32 = closed_box_quiescence_3d::<f32>("D3Q19/f32", 5.0e-7);
    eprintln!(
        "measured VR-STR-06 residuals: D2/f64={d2_f64:e}, D2/f32={d2_f32:e}, D3/f64={d3_f64:e}, D3/f32={d3_f32:e}"
    );
    assert!(
        d2_f64 < 2.0e-9 && d3_f64 < 1.0e-13 && d2_f32 < 5.0e-7 && d3_f32 < 5.0e-7,
        "measured VR-STR-06 residuals: D2/f64={d2_f64:e}, D2/f32={d2_f32:e}, D3/f64={d3_f64:e}, D3/f32={d3_f32:e}"
    );
}

#[test]
fn gravity_channel_is_bit_identical_to_raw_rho_g_force_field() {
    let (nx, ny) = (48, 18);
    let g = [8.0e-7, 0.0];
    let build = || -> Simulation<f64> {
        SimConfig {
            nx,
            ny,
            nu: 0.1,
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
        .unwrap()
    };
    let mut grav = build();
    grav.set_gravity(g);
    let mut raw = build();
    for _ in 0..600 {
        let mut force = vec![[0.0; 2]; nx * ny];
        for y in 0..ny {
            for x in 0..nx {
                if !raw.is_solid(x, y) {
                    force[y * nx + x] = [raw.rho(x, y) * g[0], raw.rho(x, y) * g[1]];
                }
            }
        }
        raw.force_field_mut().copy_from_slice(&force);
        grav.step();
        raw.step();
    }
    for y in 0..ny {
        for x in 0..nx {
            assert_eq!(
                grav.rho(x, y).to_bits(),
                raw.rho(x, y).to_bits(),
                "rho differs at ({x},{y})"
            );
            assert_eq!(
                grav.ux(x, y).to_bits(),
                raw.ux(x, y).to_bits(),
                "ux differs at ({x},{y})"
            );
            assert_eq!(
                grav.uy(x, y).to_bits(),
                raw.uy(x, y).to_bits(),
                "uy differs at ({x},{y})"
            );
        }
    }
}

#[test]
fn tilted_gravity_channel_matches_poiseuille_profile() {
    let (nx, ny) = (64, 34);
    let nu = 0.1;
    let gx = 1.0e-6;
    let mut sim = SimConfig {
        nx,
        ny,
        nu,
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
    sim.set_gravity([gx, -2.0e-7]);
    sim.run(25_000);

    let h = (ny - 2) as f64;
    let mut linf_rel = 0.0f64;
    for y in 1..ny - 1 {
        let d = y as f64 - 0.5;
        let theory = gx * d * (h - d) / (2.0 * nu);
        let measured = sim.ux(nx / 2, y);
        let rel = ((measured - theory) / theory.max(1.0e-30)).abs();
        linf_rel = linf_rel.max(rel);
    }
    assert!(
        linf_rel <= 2.0e-2,
        "tilted gravity Poiseuille L_inf_rel={linf_rel:e}, bound=2e-2"
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
fn shan_chen_gravity_composes_with_additive_force_field_and_creates_buoyancy() {
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
            sim.force_field_mut().fill([0.0; 2]);
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
