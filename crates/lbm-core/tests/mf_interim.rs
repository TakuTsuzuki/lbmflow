//! Adversarial acceptance tests for the MF-interim Wave-1 contracts.
//!
//! The implementation-facing tests are gated behind `mf-interim` because the
//! implementation branches are developed concurrently. The always-on tests
//! freeze the independent algebra/reference helpers used by those gated tests.

const EPS: f64 = 1.0e-12;

#[derive(Clone, Copy, Debug)]
struct RotorGeometry {
    center: [f64; 2],
    n_blades: usize,
    r_hub: f64,
    r_blade: f64,
    blade_thickness: f64,
    theta: f64,
}

impl RotorGeometry {
    fn contains(&self, p: [f64; 2]) -> bool {
        let dx = p[0] - self.center[0];
        let dy = p[1] - self.center[1];
        let r = dx.hypot(dy);
        if !(self.r_hub..=self.r_blade).contains(&r) {
            return false;
        }
        let blade_period = std::f64::consts::TAU / self.n_blades as f64;
        let phi = (dy.atan2(dx) - self.theta).rem_euclid(blade_period);
        let signed_phi = if phi > 0.5 * blade_period {
            phi - blade_period
        } else {
            phi
        };
        (r * signed_phi.sin()).abs() <= 0.5 * self.blade_thickness
    }
}

fn penalization_force(rho: f64, chi: f64, u_target: [f64; 2], u_star: [f64; 2]) -> [f64; 2] {
    [
        2.0 * rho * chi * (u_target[0] - u_star[0]),
        2.0 * rho * chi * (u_target[1] - u_star[1]),
    ]
}

fn physical_velocity_after_penalization(
    rho: f64,
    chi: f64,
    u_target: [f64; 2],
    u_star: [f64; 2],
) -> [f64; 2] {
    let f = penalization_force(rho, chi, u_target, u_star);
    [u_star[0] + 0.5 * f[0] / rho, u_star[1] + 0.5 * f[1] / rho]
}

fn schiller_naumann_factor(re: f64) -> f64 {
    assert!(
        re <= 800.0,
        "Schiller-Naumann reference factor is outside its validity domain: Re_p={re:e}"
    );
    1.0 + 0.15 * re.powf(0.687)
}

fn sn_tau_p(speed: f64, d: f64, rho_p: f64, rho_f: f64, nu: f64) -> f64 {
    let re = speed * d / nu;
    rho_p * d * d / (18.0 * rho_f * nu * schiller_naumann_factor(re))
}

fn terminal_speed_schiller_naumann(d: f64, rho_p: f64, rho_f: f64, nu: f64, g: f64) -> f64 {
    let g_eff = (1.0 - rho_f / rho_p) * g;
    let mut v = g_eff.abs() * sn_tau_p(0.0, d, rho_p, rho_f, nu);
    for _ in 0..80 {
        v = g_eff.abs() * sn_tau_p(v, d, rho_p, rho_f, nu);
    }
    v
}

fn particle_free_update(v: [f64; 3], u_f: [f64; 3], tau_p: f64, g_eff: [f64; 3]) -> [f64; 3] {
    [
        (v[0] + u_f[0] / tau_p + g_eff[0]) / (1.0 + 1.0 / tau_p),
        (v[1] + u_f[1] / tau_p + g_eff[1]) / (1.0 + 1.0 / tau_p),
        (v[2] + u_f[2] / tau_p + g_eff[2]) / (1.0 + 1.0 / tau_p),
    ]
}

#[test]
fn sn_terminal_velocity_solver_satisfies_fixed_point() {
    let (d, rho_p, rho_f, nu, g) = (0.37, 2500.0, 1000.0, 0.012, 9.81e-4);
    let vt = terminal_speed_schiller_naumann(d, rho_p, rho_f, nu, g);
    let rhs = (1.0 - rho_f / rho_p) * g * sn_tau_p(vt, d, rho_p, rho_f, nu);
    assert!((vt - rhs).abs() / vt < 1.0e-13, "vt={vt} rhs={rhs}");

    let stokes = (1.0 - rho_f / rho_p) * g * sn_tau_p(0.0, d, rho_p, rho_f, nu);
    assert!(
        vt < stokes,
        "Schiller-Naumann correction must reduce the Stokes terminal speed"
    );
}

#[test]
fn rotor_blade_indicator_matches_axis_aligned_analytic_point_set() {
    let geom = RotorGeometry {
        center: [10.0, 10.0],
        n_blades: 4,
        r_hub: 1.5,
        r_blade: 4.5,
        blade_thickness: 0.75,
        theta: 0.0,
    };

    let expected_inside = [
        [12.0, 10.0],
        [14.0, 10.0],
        [8.0, 10.0],
        [6.0, 10.0],
        [10.0, 12.0],
        [10.0, 14.0],
        [10.0, 8.0],
        [10.0, 6.0],
    ];
    for p in expected_inside {
        assert!(geom.contains(p), "axis blade point {p:?} was missed");
    }

    let expected_outside = [
        [11.0, 10.0], // hub exclusion.
        [15.0, 10.0], // beyond r_blade.
        [13.0, 11.0], // outside thickness.
        [11.0, 13.0], // outside thickness on the vertical blade.
        [13.0, 13.0], // diagonal gap between blades.
    ];
    for p in expected_outside {
        assert!(!geom.contains(p), "non-blade point {p:?} was included");
    }
}

#[test]
fn semi_implicit_particle_update_matches_hand_derived_two_step_sequence() {
    let tau_p = 3.0;
    let u = [1.2, -0.6, 0.3];
    let g_eff = [0.0, -0.2, 0.05];
    let v0 = [0.3, 0.9, -0.1];

    let a = 1.0 / (1.0 + 1.0 / tau_p);
    let b = [
        (u[0] / tau_p + g_eff[0]) * a,
        (u[1] / tau_p + g_eff[1]) * a,
        (u[2] / tau_p + g_eff[2]) * a,
    ];
    let manual_v1 = [a * v0[0] + b[0], a * v0[1] + b[1], a * v0[2] + b[2]];
    let manual_v2 = [
        a * a * v0[0] + (1.0 + a) * b[0],
        a * a * v0[1] + (1.0 + a) * b[1],
        a * a * v0[2] + (1.0 + a) * b[2],
    ];

    let v1 = particle_free_update(v0, u, tau_p, g_eff);
    let v2 = particle_free_update(v1, u, tau_p, g_eff);
    for i in 0..3 {
        assert!((v1[i] - manual_v1[i]).abs() < EPS);
        assert!((v2[i] - manual_v2[i]).abs() < EPS);
    }
}

#[test]
fn penalization_force_formula_produces_exact_physical_velocity_contract() {
    let rho = 1.7;
    let u_star = [-0.04, 0.11];
    let u_target = [0.23, -0.02];
    for chi in [0.0, 0.125, 0.5, 1.0] {
        let u = physical_velocity_after_penalization(rho, chi, u_target, u_star);
        let expected = [
            u_star[0] + chi * (u_target[0] - u_star[0]),
            u_star[1] + chi * (u_target[1] - u_star[1]),
        ];
        assert!((u[0] - expected[0]).abs() < EPS, "chi={chi} u={u:?}");
        assert!((u[1] - expected[1]).abs() < EPS, "chi={chi} u={u:?}");
    }
}

#[cfg(feature = "mf-interim")]
mod mf_contract_tests {
    use super::*;
    use lbm_core::compat::multiphase::ShanChen;
    use lbm_core::compat::prelude::*;
    use lbm_core::compat::rotor::Rotor;
    use lbm_core::particles::{sample_grid, Particle, ParticleSet, Sample};
    use lbm_core::prelude::{
        build_wall_rims, CollisionKind, CpuScalar, FaceBC, GlobalSpec, LocalPeriodic, Solver,
        WallSpec, D2Q9, D3Q19,
    };

    fn periodic_compat(nx: usize, ny: usize) -> Simulation<f64> {
        SimConfig::<f64> {
            nx,
            ny,
            nu: 1.0 / 6.0,
            collision: Collision::Trt { magic: 0.1875 },
            edges: Edges {
                left: EdgeBC::Periodic,
                right: EdgeBC::Periodic,
                bottom: EdgeBC::Periodic,
                top: EdgeBC::Periodic,
            },
            ..Default::default()
        }
        .build()
        .unwrap()
    }

    fn closed_compat(nx: usize, ny: usize) -> Simulation<f64> {
        SimConfig::<f64> {
            nx,
            ny,
            nu: 1.0 / 6.0,
            collision: Collision::Trt { magic: 0.1875 },
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
    fn gravity_survives_shan_chen_force_overwrite_each_step() {
        let mut sim = periodic_compat(64, 96);
        sim.init_with(|_, y| {
            let light_blob = (y as isize - 32).abs() <= 10;
            (if light_blob { 0.15 } else { 2.0 }, 0.0, 0.0)
        });
        sim.set_gravity([0.0, 2.0e-6]);
        let sc = ShanChen::new(-5.0);

        let mut y0 = 0.0;
        let mut m0 = 0.0;
        for y in 0..sim.ny() {
            for x in 0..sim.nx() {
                let deficit = (2.0 - sim.rho(x, y)).max(0.0);
                y0 += deficit * y as f64;
                m0 += deficit;
            }
        }
        for _ in 0..250 {
            sc.update_force(&mut sim);
            sim.step();
        }
        let mut y1 = 0.0;
        let mut m1 = 0.0;
        for y in 0..sim.ny() {
            for x in 0..sim.nx() {
                let deficit = (2.0 - sim.rho(x, y)).max(0.0);
                y1 += deficit * y as f64;
                m1 += deficit;
            }
        }
        assert!(y1 / m1 > y0 / m0 + 0.25, "light blob did not rise");
    }

    #[test]
    fn gravity_momentum_growth_excludes_solid_obstacles() {
        // Same-path twin: set_gravity(g) must equal a hand-maintained
        // per-cell force field rho(x)*g written on FLUID cells only, every
        // step, through the same force-field kernel. This pins solid
        // exclusion, rho weighting and composition bit-tightly and is immune
        // to ANOM-P2-001 (the uniform-force vs force-field transient
        // inconsistency, docs/qa/anomaly-log.md), which blocks any tight
        // CROSS-path twin.
        let obstacle = |x: usize, y: usize| (10..15).contains(&x) && (6..12).contains(&y);
        let g = [3.0e-7, -4.0e-7];
        let mut sim_g = periodic_compat(32, 24);
        sim_g.set_solid_region(obstacle);
        sim_g.set_gravity(g);
        let mut sim_f = periodic_compat(32, 24);
        sim_f.set_solid_region(obstacle);

        for _ in 0..50 {
            let (nx, ny) = (sim_f.nx(), sim_f.ny());
            let mut forces = vec![[0.0f64; 2]; nx * ny];
            for y in 0..ny {
                for x in 0..nx {
                    if !sim_f.is_solid(x, y) {
                        let rho = sim_f.rho(x, y);
                        forces[y * nx + x] = [rho * g[0], rho * g[1]];
                    }
                }
            }
            sim_f.force_field_mut().copy_from_slice(&forces);
            sim_g.step();
            sim_f.step();
        }
        let mut max_du = 0.0f64;
        let mut max_u = 0.0f64;
        for y in 0..sim_g.ny() {
            for x in 0..sim_g.nx() {
                if sim_g.is_solid(x, y) {
                    continue;
                }
                max_du = max_du
                    .max((sim_g.ux(x, y) - sim_f.ux(x, y)).abs())
                    .max((sim_g.uy(x, y) - sim_f.uy(x, y)).abs())
                    .max((sim_g.rho(x, y) - sim_f.rho(x, y)).abs());
                max_u = max_u.max(sim_g.ux(x, y).abs());
            }
        }
        assert!(max_u > 1.0e-6, "flow did not develop: {max_u}");
        assert!(
            max_du < 1.0e-13,
            "gravity != rho*g force-field twin: {max_du}"
        );
    }

    #[test]
    fn gravity_2d_matches_z_thin_3d_periodic_degeneracy() {
        let (nx, ny, nz) = (16, 12, 4);
        let mut sim2 = periodic_compat(nx, ny);
        sim2.set_gravity([1.0e-7, -2.0e-7]);
        let spec = GlobalSpec::<f64> {
            dims: [nx, ny, nz],
            nu: 1.0 / 6.0,
            collision: CollisionKind::Trt { magic: 0.1875 },
            periodic: [true, true, true],
            faces: [FaceBC::Closed; 6],
            ..Default::default()
        };
        let n = nx * ny * nz;
        let mut sim3: Solver<D3Q19, f64, CpuScalar, LocalPeriodic> = Solver::new(
            &spec,
            &vec![false; n],
            &vec![[0.0; 3]; n],
            [1, 1, 1],
            CpuScalar::default(),
            LocalPeriodic,
        );
        sim3.set_gravity([1.0e-7, -2.0e-7, 0.0]);
        for _ in 0..40 {
            sim2.step();
            sim3.step();
        }
        let rho3 = sim3.gather_rho();
        let ux3 = sim3.gather_ux();
        let uy3 = sim3.gather_uy();
        for y in 0..ny {
            for x in 0..nx {
                let i2 = y * nx + x;
                let i3 = y * nx + x;
                assert!((sim2.rho_field()[i2] - rho3[i3]).abs() < 1.0e-12);
                assert!((sim2.ux_field()[i2] - ux3[i3]).abs() < 1.0e-12);
                assert!((sim2.uy_field()[i2] - uy3[i3]).abs() < 1.0e-12);
            }
        }
    }

    #[test]
    fn rotor_chi_one_edge_cells_never_overshoot_target_velocity() {
        // The chi = 1 pinning identity u_phys = u_target holds at the forcing
        // stage; streaming then mixes in populations from NON-penalized
        // neighbors, so post-step blade-EDGE cells legitimately deviate.
        // Assert the frozen tracking band on cells whose full neighborhood is
        // penalized (catches u_star bookkeeping bugs) and only boundedness at
        // the edge.
        let mut sim = closed_compat(48, 48);
        let mut rotor = Rotor::new(24.0, 24.0)
            .n_blades(3)
            .r_hub(3.0)
            .r_blade(17.0)
            .blade_thickness(1.2)
            .omega(0.012)
            .chi(1.0)
            .omega_ramp_steps(50)
            .theta0(0.0);
        let u_tip = 0.012 * 17.0;
        for step in 0..2_000 {
            sim.clear_force_field();
            rotor.update_force(&mut sim);
            sim.step();
            if step < 200 {
                continue; // motor ramp + startup transient
            }
            for y in 1..sim.ny() - 1 {
                for x in 1..sim.nx() - 1 {
                    let p = [x as f64, y as f64];
                    if !rotor.contains(p) {
                        continue;
                    }
                    let interior = (-1..=1).all(|dy: isize| {
                        (-1..=1)
                            .all(|dx: isize| rotor.contains([p[0] + dx as f64, p[1] + dy as f64]))
                    });
                    if !interior {
                        // Freshly swept-in edge cells legitimately carry
                        // ambient flow (inter-blade jets can locally exceed
                        // u_tip); their only guarantee is the global 0.3
                        // ceiling, asserted elsewhere.
                        continue;
                    }
                    let u_t = rotor.target_velocity(p);
                    let u = [sim.ux(x, y), sim.uy(x, y)];
                    let err = (u[0] - u_t[0]).hypot(u[1] - u_t[1]);
                    let band = 0.02 * u_tip;
                    assert!(
                        err <= band,
                        "tracking error {err} > {band} at ({x},{y}) \
                         u={u:?} target={u_t:?}"
                    );
                }
            }
        }
    }

    #[test]
    fn rotor_torque_integral_matches_fluid_angular_momentum_gain_order_and_sign() {
        let mut sim = closed_compat(64, 64);
        let c = [32.0, 32.0];
        let mut rotor = Rotor::new(c[0], c[1])
            .n_blades(4)
            .r_hub(4.0)
            .r_blade(20.0)
            .blade_thickness(1.5)
            .omega(0.01)
            .chi(0.6)
            .omega_ramp_steps(100)
            .theta0(0.0);
        let l0 = fluid_angular_momentum_2d(&sim, c);
        let mut torque_integral: f64 = 0.0;
        for _ in 0..300 {
            sim.clear_force_field();
            rotor.update_force(&mut sim);
            torque_integral += rotor.torque();
            sim.step();
        }
        let dl = fluid_angular_momentum_2d(&sim, c) - l0;
        // torque() is the REACTION on the rotor (module doc); the torque on
        // the FLUID is its negation.
        let fluid_torque_integral = -torque_integral;
        assert!(
            dl.signum() == fluid_torque_integral.signum(),
            "dl={dl} fluid_torque={fluid_torque_integral}"
        );
        let ratio = (dl / fluid_torque_integral).abs();
        assert!(
            (0.5..=2.0).contains(&ratio),
            "angular momentum gain and torque integral mismatch: dl={dl} torque={torque_integral}"
        );
    }

    fn fluid_angular_momentum_2d(sim: &Simulation<f64>, c: [f64; 2]) -> f64 {
        let mut l = 0.0;
        for y in 0..sim.ny() {
            for x in 0..sim.nx() {
                if !sim.is_solid(x, y) {
                    let r = [x as f64 - c[0], y as f64 - c[1]];
                    l += sim.rho(x, y) * (r[0] * sim.uy(x, y) - r[1] * sim.ux(x, y));
                }
            }
        }
        l
    }

    #[test]
    fn rotor_and_gravity_force_writes_are_additive_and_tracking_is_unchanged() {
        let mut no_g = closed_compat(40, 40);
        let mut with_g = closed_compat(40, 40);
        with_g.set_gravity([2.0e-7, -3.0e-7]);
        let mut rotor0 = Rotor::new(20.0, 20.0)
            .n_blades(2)
            .r_hub(3.0)
            .r_blade(13.0)
            .blade_thickness(1.0)
            .omega(0.01)
            .chi(1.0)
            .omega_ramp_steps(1)
            .theta0(0.0);
        let mut rotor1 = Rotor::new(20.0, 20.0)
            .n_blades(2)
            .r_hub(3.0)
            .r_blade(13.0)
            .blade_thickness(1.0)
            .omega(0.01)
            .chi(1.0)
            .omega_ramp_steps(1)
            .theta0(0.0);
        for _ in 0..80 {
            no_g.clear_force_field();
            with_g.clear_force_field();
            rotor0.update_force(&mut no_g);
            rotor1.update_force(&mut with_g);
            no_g.step();
            with_g.step();
        }
        let e0 = rotor_tracking_error(&no_g, &rotor0);
        let e1 = rotor_tracking_error(&with_g, &rotor1);
        // Gravity stratifies rho hydrostatically (delta rho ~ g*H/cs^2), and
        // the penalization force is rho-weighted, so tracking WILL shift at
        // that order; a clobbered force write would shift it at O(u_tip).
        assert!(
            (e0 - e1).abs() < 0.1 * e0.max(1.0e-30),
            "no_g={e0} with_g={e1}"
        );
    }

    fn rotor_tracking_error(sim: &Simulation<f64>, rotor: &Rotor<f64>) -> f64 {
        let mut max = 0.0f64;
        for y in 1..sim.ny() - 1 {
            for x in 1..sim.nx() - 1 {
                let p = [x as f64, y as f64];
                if rotor.contains(p) {
                    let u = [sim.ux(x, y), sim.uy(x, y)];
                    let t = rotor.target_velocity(p);
                    max = max.max((u[0] - t[0]).hypot(u[1] - t[1]));
                }
            }
        }
        max
    }

    #[test]
    fn particle_on_lattice_node_and_solid_corner_is_finite_and_deterministic() {
        let make = || ParticleSet {
            particles: vec![Particle {
                pos: [1.0, 1.0, 0.0],
                vel: [0.2, -0.3, 0.0],
                d: 0.2,
                rho_p: 2.0,
                exposure: 0.0,
            }],
            rho_f: 1.0,
            nu: 0.1,
            g: [0.0, -0.01, 0.0],
            restitution: 0.3,
        };
        let sample = |p: [f64; 3]| Sample {
            u: [0.0, 0.0, 0.0],
            solid: p[0] <= 1.0 && p[1] <= 1.0,
        };
        let mut a = make();
        let mut b = make();
        a.step(sample, None::<fn([f64; 3]) -> f64>).unwrap();
        b.step(sample, None::<fn([f64; 3]) -> f64>).unwrap();
        assert_eq!(a.particles[0].pos, b.particles[0].pos);
        assert_eq!(a.particles[0].vel, b.particles[0].vel);
        assert!(a.particles[0].pos.iter().all(|v| v.is_finite()));
        assert!(a.particles[0].vel.iter().all(|v| v.is_finite()));
    }

    #[test]
    fn particle_tracer_in_solid_body_rotation_preserves_radius_heavy_particle_drifts_outward() {
        let omega = 0.02;
        let center = [32.0, 32.0, 0.0];
        let sample = |p: [f64; 3]| {
            let r = [p[0] - center[0], p[1] - center[1], 0.0];
            Sample {
                u: [-omega * r[1], omega * r[0], 0.0],
                solid: false,
            }
        };
        let mut tracer = ParticleSet {
            particles: vec![Particle {
                pos: [42.0, 32.0, 0.0],
                vel: [0.0, omega * 10.0, 0.0],
                d: 1.0e-4,
                rho_p: 1.0001,
                exposure: 0.0,
            }],
            rho_f: 1.0,
            nu: 0.1,
            g: [0.0, 0.0, 0.0],
            restitution: 0.0,
        };
        let mut heavy = ParticleSet {
            particles: vec![Particle {
                pos: [42.0, 32.0, 0.0],
                vel: [0.0, omega * 10.0, 0.0],
                d: 1.0,
                rho_p: 1000.0,
                exposure: 0.0,
            }],
            rho_f: 1.0,
            nu: 0.1,
            g: [0.0, 0.0, 0.0],
            restitution: 0.0,
        };
        for _ in 0..((std::f64::consts::TAU / omega).round() as usize) {
            tracer.step(sample, None::<fn([f64; 3]) -> f64>).unwrap();
            heavy.step(sample, None::<fn([f64; 3]) -> f64>).unwrap();
        }
        let r_tracer = radius(tracer.particles[0].pos, center);
        let r_heavy = radius(heavy.particles[0].pos, center);
        // The contract pins forward-Euler position updates (dt = 1), whose
        // exact solid-body-rotation radius growth is (1 + omega^2)^(n/2):
        // a tracer must reproduce it, not beat it. The systematic drift is
        // filed as an S2 improvement (2nd-order particle advection) in
        // docs/qa/anomaly-log.md.
        let n_steps = (std::f64::consts::TAU / omega).round();
        let r_euler = 10.0 * (1.0 + omega * omega).powf(n_steps / 2.0);
        assert!(
            (r_tracer - r_euler).abs() < 0.01 * r_euler,
            "tracer radius={r_tracer} euler={r_euler}"
        );
        assert!(
            r_heavy > r_tracer + 0.1,
            "heavy={r_heavy} tracer={r_tracer}"
        );
    }

    fn radius(p: [f64; 3], c: [f64; 3]) -> f64 {
        ((p[0] - c[0]).powi(2) + (p[1] - c[1]).powi(2)).sqrt()
    }

    #[test]
    fn settled_particle_resuspends_when_imposed_shear_exceeds_force_balance() {
        let rho_f = 1.0;
        let rho_p = 2.0;
        let d = 0.5;
        let nu = 0.1;
        let g = [0.0, -1.0e-3, 0.0];
        let g_eff_y = (1.0 - rho_f / rho_p) * g[1];
        let tau = sn_tau_p(0.0, d, rho_p, rho_f, nu);
        let above_threshold_uy = -2.0 * tau * g_eff_y;
        // The upward flow must reach the rested particle (it sits just above
        // the floor at y ~ 0), and the lift speed is ~tau*|g_eff| per step,
        // so crossing y = 0.5 needs O(0.5 / (tau*|g_eff|)) ~ 3.6k steps.
        let sample = move |p: [f64; 3]| Sample {
            u: [0.0, if p[1] > 0.0 { above_threshold_uy } else { 0.0 }, 0.0],
            solid: p[1] <= 0.0,
        };
        let mut ps = ParticleSet {
            particles: vec![Particle {
                pos: [4.0, 0.05, 0.0],
                vel: [0.0, 0.0, 0.0],
                d,
                rho_p,
                exposure: 0.0,
            }],
            rho_f,
            nu,
            g,
            restitution: 0.0,
        };
        for _ in 0..4_000 {
            ps.step(sample, None::<fn([f64; 3]) -> f64>).unwrap();
        }
        assert!(
            ps.particles[0].pos[1] > 0.5,
            "particle remained stuck: {:?}",
            ps.particles[0]
        );
    }

    #[test]
    fn rested_particle_on_moving_floor_acquires_horizontal_velocity_from_drag() {
        let sample = |p: [f64; 3]| Sample {
            u: if p[1] <= 0.0 {
                [0.0, 0.0, 0.0]
            } else {
                [0.1, 0.0, 0.0]
            },
            solid: p[1] <= 0.0,
        };
        let mut ps = ParticleSet {
            particles: vec![Particle {
                pos: [3.0, 0.02, 0.0],
                vel: [0.0, 0.0, 0.0],
                d: 0.2,
                rho_p: 10.0,
                exposure: 0.0,
            }],
            rho_f: 1.0,
            nu: 0.1,
            g: [0.0, -1.0e-3, 0.0],
            restitution: 0.0,
        };
        for _ in 0..8 {
            ps.step(sample, None::<fn([f64; 3]) -> f64>).unwrap();
        }
        assert!(
            ps.particles[0].vel[0] > 0.0,
            "drag ignored on rested particle"
        );
    }

    #[test]
    fn cfl_guard_prevents_tunneling_through_one_cell_wall() {
        let sample = |p: [f64; 3]| Sample {
            u: [0.0, 0.0, 0.0],
            solid: p[0] >= 5.0 && p[0] < 6.0,
        };
        let mut ps = ParticleSet {
            particles: vec![Particle {
                pos: [3.25, 2.0, 0.0],
                vel: [2.5, 0.0, 0.0],
                d: 0.2,
                rho_p: 1000.0,
                exposure: 0.0,
            }],
            rho_f: 1.0,
            nu: 0.1,
            g: [0.0, 0.0, 0.0],
            restitution: 0.0,
        };
        ps.step(sample, None::<fn([f64; 3]) -> f64>).unwrap();
        assert!(
            ps.particles[0].pos[0] < 5.0,
            "particle tunneled through wall: {:?}",
            ps.particles[0]
        );
    }

    #[test]
    fn sample_grid_uses_zero_velocity_for_solid_nodes() {
        let u = vec![[1.0, 0.0, 0.0]; 8];
        let solid = vec![false, true, false, false, false, false, false, false];
        let s = sample_grid([0.75, 0.25, 0.0], [2, 2, 2], |x, y, z| {
            let i = (z * 2 + y) * 2 + x;
            (u[i], solid[i])
        });
        assert!(
            s.u[0] < 1.0,
            "solid node velocity was included instead of zeroed: {:?}",
            s
        );
    }

    #[test]
    #[ignore = "SPEC-GAP: contract must define particle-particle collisions or explicitly forbid overlap"]
    fn spec_gap_two_particles_occupying_one_cell() {
        let mut ps = ParticleSet {
            particles: vec![
                Particle {
                    pos: [2.0, 2.0, 0.0],
                    vel: [1.0, 0.0, 0.0],
                    d: 0.5,
                    rho_p: 10.0,
                    exposure: 0.0,
                },
                Particle {
                    pos: [2.0, 2.0, 0.0],
                    vel: [-1.0, 0.0, 0.0],
                    d: 0.5,
                    rho_p: 10.0,
                    exposure: 0.0,
                },
            ],
            rho_f: 1.0,
            nu: 0.1,
            g: [0.0, 0.0, 0.0],
            restitution: 0.0,
        };
        ps.step(
            |_| Sample {
                u: [0.0; 3],
                solid: false,
            },
            None::<fn([f64; 3]) -> f64>,
        )
        .unwrap();
        assert!(ps.particles[0].pos != ps.particles[1].pos);
    }

    #[test]
    #[ignore = "SPEC-GAP: contract must define whether a particle starting inside solid is projected, reflected, or rejected"]
    fn spec_gap_particle_starting_inside_solid() {
        let mut ps = ParticleSet {
            particles: vec![Particle {
                pos: [0.0, 0.0, 0.0],
                vel: [0.0, 0.0, 0.0],
                d: 0.5,
                rho_p: 10.0,
                exposure: 0.0,
            }],
            rho_f: 1.0,
            nu: 0.1,
            g: [0.0, 0.0, 0.0],
            restitution: 0.0,
        };
        ps.step(
            |_| Sample {
                u: [0.0; 3],
                solid: true,
            },
            None::<fn([f64; 3]) -> f64>,
        )
        .unwrap();
        assert!(!ps.particles[0].pos.iter().any(|v| v.is_nan()));
    }

    #[test]
    #[ignore = "SPEC-GAP: rotor contract says chi in (0,1], but must pin whether chi=0 is rejected or is a no-op"]
    fn spec_gap_rotor_chi_zero() {
        let mut sim = closed_compat(16, 16);
        let mut rotor = Rotor::new(8.0, 8.0)
            .n_blades(2)
            .r_hub(1.0)
            .r_blade(5.0)
            .blade_thickness(1.0)
            .omega(0.1)
            .chi(0.0)
            .omega_ramp_steps(1)
            .theta0(0.0);
        rotor.update_force(&mut sim);
        sim.step();
    }

    #[test]
    #[ignore = "SPEC-GAP: native Solver::set_gravity must define how it composes with existing force fields across subdomains"]
    fn spec_gap_native_gravity_after_body_force_field_rewrite() {
        let spec = GlobalSpec::<f64> {
            dims: [8, 8, 1],
            periodic: [true, true, false],
            ..Default::default()
        };
        let mut s: Solver<D2Q9, f64, CpuScalar, LocalPeriodic> = Solver::new(
            &spec,
            &vec![false; 64],
            &vec![[0.0; 3]; 64],
            [1, 1, 1],
            CpuScalar::default(),
            LocalPeriodic,
        );
        s.set_gravity([0.0, 1.0e-6, 0.0]);
        s.set_body_force_field(|_, _, _| [1.0e-6, 0.0, 0.0]);
        s.step();
    }
}
