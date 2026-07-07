//! ACC-AUDIT drift pins for code-to-spec radar rows 27/28/29/31.
//!
//! Source: `docs/qa/code-to-spec-diff.md` Section B, rows A12/A25/A4/A28.
//! These probes are deliberately narrow: they print the measured behavior at
//! the drift point and avoid engine changes.

mod common;

use common::metrics::{curve_agreement, l2_rel};
use lbm_core::compat::prelude::{
    Collision as CompatCollision, EdgeBC, Edges, SimConfig, Simulation,
};
use lbm_core::particles::sample_grid;
use lbm_core::prelude::*;

const TRT: CompatCollision = CompatCollision::Trt {
    magic: CompatCollision::MAGIC_STD,
};

fn omega_from_nu(nu: f64) -> f64 {
    1.0 / (3.0 * nu + 0.5)
}

fn max_abs(xs: &[f64]) -> f64 {
    xs.iter().map(|v| v.abs()).fold(0.0, f64::max)
}

fn assert_rel_close(got: f64, want: f64, rel: f64, label: &str) {
    let err = (got - want).abs();
    let denom = want.abs().max(1.0e-300);
    assert!(
        err <= rel * denom,
        "{label}: got={got:.16e}, want={want:.16e}, abs_err={err:.3e}, rel_err={:.3e}, rel_band={rel:.3e}, denominator=|expected|",
        err / denom
    );
}

#[test]
#[should_panic(expected = "T27 A12 Couette sampler deviates")]
fn t27_particles_near_wall_sampler_preserves_couette_profile_to_moving_wall() {
    // Known-drift pin: this is the requested spec assertion, not the current
    // sampler contract. Current main samples solid moving-wall neighbors as
    // zero velocity, so this test must panic until A12 is fixed/retightened.
    let u_wall = 0.05;
    let tau = 0.8;
    let nu = (tau - 0.5) / 3.0;
    let mut sim: Simulation<f64> = SimConfig {
        nx: 8,
        ny: 34,
        nu,
        collision: TRT,
        edges: Edges {
            left: EdgeBC::Periodic,
            right: EdgeBC::Periodic,
            bottom: EdgeBC::BounceBack,
            top: EdgeBC::MovingWall { u: [u_wall, 0.0] },
        },
        ..Default::default()
    }
    .build()
    .unwrap();
    assert!(
        common::run_to_steady(&mut sim, 500, 1.0e-11, 200_000),
        "T27 Couette fixture did not reach steady state"
    );

    let h = (sim.ny() - 2) as f64;
    let y_wall = sim.ny() as f64 - 1.5;
    let y_samples = [
        sim.ny() as f64 - 4.0,
        sim.ny() as f64 - 3.0,
        sim.ny() as f64 - 2.5,
        sim.ny() as f64 - 2.25,
        sim.ny() as f64 - 2.0,
        sim.ny() as f64 - 1.75,
        y_wall,
    ];
    let mut samples = Vec::new();
    let mut profile = Vec::new();
    for &y in &y_samples {
        let sample = sample_grid([3.0, y, 0.0], [sim.nx(), sim.ny(), 1], |x, yy, _| {
            ([sim.ux(x, yy), sim.uy(x, yy), 0.0], sim.is_solid(x, yy))
        });
        samples.push((y, sample.u[0]));
        profile.push((y, sample.u[0], u_wall * (y - 0.5) / h, sample.solid));
    }
    println!("T27 Couette sampler near moving wall profile: {profile:?}");

    // Half-way wall geometry places the top moving wall at y=ny-1.5. The
    // Couette solution is affine in the distance from the bottom wall:
    // u_x(y) = U * (y - 0.5) / H, H = ny - 2. A trilinear sampler that is
    // faithful up to the wall surface must stay on this line at sub-cell
    // positions between the last fluid cell center and the moving wall.
    let agreement = curve_agreement(|y| u_wall * (y - 0.5) / h, &samples, 1.0e-9, u_wall);
    assert!(
        agreement.max_rel_dev <= 1.0e-9,
        "T27 A12 Couette sampler deviates from linear moving-wall profile: max_rel_dev={:.6e} at y={:.6}, band=1e-9, denominator=max(|u_ref|,{u_wall:.3e}); profile={profile:?}",
        agreement.max_rel_dev,
        agreement.worst_x
    );
}

fn one_step_y_momentum_gain(force: [f64; 2], gravity: Option<[f64; 2]>) -> (f64, usize) {
    let mut sim: Simulation<f64> = SimConfig {
        nx: 16,
        ny: 12,
        nu: 1.0 / 6.0,
        collision: TRT,
        force,
        ..Default::default()
    }
    .build()
    .unwrap();
    if let Some(g) = gravity {
        sim.set_gravity(g);
    }
    let p0 = sim.total_momentum();
    sim.step();
    let p1 = sim.total_momentum();
    (p1[1] - p0[1], sim.fluid_cell_count())
}

#[test]
fn t28_gravity_and_uniform_force_one_step_momentum_is_additive() {
    let g_uniform = 2.0e-7;
    let g_alt = 7.5e-8;
    let (uniform, n0) = one_step_y_momentum_gain([0.0, g_uniform], None);
    let (gravity, n1) = one_step_y_momentum_gain([0.0, 0.0], Some([0.0, -g_alt]));
    let (both, n2) = one_step_y_momentum_gain([0.0, g_uniform], Some([0.0, -g_alt]));
    assert_eq!(n0, n1);
    assert_eq!(n0, n2);
    let n = n0 as f64;
    let expect_uniform = n * g_uniform;
    let expect_gravity = -n * g_alt;
    let expect_both = n * (g_uniform - g_alt);
    let additive = uniform + gravity;
    println!(
        "T28 one-step y momentum gains: n={n0}, uniform={uniform:.16e} (expect {expect_uniform:.16e}), gravity={gravity:.16e} (expect {expect_gravity:.16e}), both={both:.16e} (expect {expect_both:.16e}), uniform+gravity={additive:.16e}, double_count_candidate={:.16e}, overwrite_candidate=0",
        2.0 * expect_both
    );

    assert_rel_close(
        uniform,
        expect_uniform,
        1.0e-12,
        "T28 uniform-only Guo impulse",
    );
    assert_rel_close(
        gravity,
        expect_gravity,
        1.0e-12,
        "T28 gravity-only Guo impulse",
    );
    assert_rel_close(
        both,
        expect_both,
        1.0e-12,
        "T28 combined force+gravity Guo impulse",
    );
    assert_rel_close(
        both,
        additive,
        1.0e-12,
        "T28 force+gravity additive superposition",
    );
}

type S2 = Solver<D2Q9, f64, CpuScalar, LocalPeriodic>;

fn central_moment_periodic(n: usize, nu: f64, omega_shear: f64, u0: [f64; 3]) -> S2 {
    let spec = GlobalSpec {
        dims: [n, n, 1],
        nu,
        collision: CollisionKind::CentralMoment { omega_shear },
        periodic: [true, true, false],
        ..Default::default()
    };
    let mut s = Solver::new(
        &spec,
        &[],
        &[],
        [1, 1, 1],
        CpuScalar::default(),
        LocalPeriodic,
    );
    s.init_with(move |_, _, _| (1.0, u0));
    s
}

fn assert_solver_finite(label: &str, s: &S2) {
    let mut fields = Vec::new();
    fields.extend(s.gather_rho());
    fields.extend(s.gather_ux());
    fields.extend(s.gather_uy());
    fields.extend(s.gather_uz());
    let max_field = max_abs(&fields);
    assert!(
        fields.iter().all(|v| v.is_finite()),
        "{label}: non-finite field after run; max_abs_field={max_field:e}"
    );
}

#[test]
fn t29_central_moment_omega_ceiling_is_silent_but_finite() {
    let nu = 0.001;
    let omega_nu = omega_from_nu(nu);
    let u = 0.15;
    let usq = u * u;
    let current_pre_ceiling = omega_nu * (1.0 - 0.16 * usq);
    let legacy_offset_pre_ceiling = omega_nu * (1.0 + 0.0025 - 0.16 * usq);
    let mut near_limit = central_moment_periodic(32, nu, omega_nu, [u, 0.0, 0.0]);
    near_limit.run(100);
    assert_solver_finite("T29 nu-derived near-limit central-moment run", &near_limit);

    // Current main validates the global CentralMoment omega_shear in (0, 2],
    // so the remaining way to hit the kernel's `.min(2.0)` ceiling is a
    // per-cell omega field. The kernel has no diagnostic or warning channel
    // for this path; this test documents that observation and only asserts
    // that the silently ceiled run stays finite.
    let mut per_cell_ceiling = central_moment_periodic(32, nu, omega_nu, [u, 0.0, 0.0]);
    let forced_omega = vec![2.01; 32 * 32];
    per_cell_ceiling.set_omega_field(Some(&forced_omega));
    let forced_pre_ceiling = forced_omega[0] * (1.0 - 0.16 * usq);
    per_cell_ceiling.run(100);
    assert_solver_finite("T29 per-cell omega>2 central-moment run", &per_cell_ceiling);
    println!(
        "T29 central-moment omega ceiling: nu={nu:.6e}, tau={:.12e}, omega_nu={omega_nu:.12e}, u={u:.6e}, current_pre_ceiling={current_pre_ceiling:.12e}, legacy_offset_pre_ceiling={legacy_offset_pre_ceiling:.12e}, forced_per_cell_pre_ceiling={forced_pre_ceiling:.12e}, clamp_would_fire_current_nu={}, clamp_would_fire_forced={}, warning_or_flag_observed=false",
        3.0 * nu + 0.5,
        current_pre_ceiling > 2.0,
        forced_pre_ceiling > 2.0
    );
}

fn multimode_state(n: usize, x: usize, y: usize) -> (f64, [f64; 3]) {
    let k = 2.0 * std::f64::consts::PI / n as f64;
    let xf = k * x as f64;
    let yf = k * y as f64;
    let u0 = 0.018;
    (
        1.0 + 1.0e-5 * (xf + 0.7 * yf).sin(),
        [
            u0 * (yf.sin() + 0.35 * (2.0 * xf + yf).sin()),
            u0 * 0.25 * (xf - 1.5 * yf).cos(),
            0.0,
        ],
    )
}

fn multimode_force(n: usize, x: usize, y: usize) -> [f64; 3] {
    let k = 2.0 * std::f64::consts::PI / n as f64;
    let xf = k * x as f64;
    let yf = k * y as f64;
    [
        2.0e-7 * yf.sin() * (1.0 + 0.25 * xf.cos()),
        5.0e-8 * (xf + yf).sin(),
        0.0,
    ]
}

fn reflected_vec(v: [f64; 3]) -> [f64; 3] {
    [-v[0], v[1], v[2]]
}

fn make_wale_mirror_pair(n: usize) -> (S2, S2) {
    let spec = GlobalSpec {
        dims: [n, n, 1],
        nu: 0.02,
        collision: CollisionKind::Trt {
            magic: CollisionKind::MAGIC_STD,
        },
        periodic: [true, true, false],
        ..Default::default()
    };
    let mut base = Solver::new(
        &spec,
        &[],
        &[],
        [1, 1, 1],
        CpuScalar::default(),
        LocalPeriodic,
    );
    let mut mirror = Solver::new(
        &spec,
        &[],
        &[],
        [1, 1, 1],
        CpuScalar::default(),
        LocalPeriodic,
    );
    base.init_with(move |x, y, _| multimode_state(n, x, y));
    mirror.init_with(move |x, y, _| {
        let mx = n - 1 - x;
        let (rho, u) = multimode_state(n, mx, y);
        (rho, reflected_vec(u))
    });
    base.set_body_force_field(move |x, y, _| multimode_force(n, x, y));
    mirror.set_body_force_field(move |x, y, _| {
        let mx = n - 1 - x;
        reflected_vec(multimode_force(n, mx, y))
    });
    (base, mirror)
}

#[test]
fn t31_wale_per_cell_omega_preserves_x_mirror_equivariance() {
    let n = 32usize;
    let (mut base, mut mirror) = make_wale_mirror_pair(n);
    let mut base_les = WaleLes::<f64>::new();
    let mut mirror_les = WaleLes::<f64>::new();
    let mut max_nu_t_first = 0.0;
    let mut max_nu_t_last = 0.0;
    for step in 0..200 {
        base_les.update(&mut base);
        mirror_les.update(&mut mirror);
        let max_nu_t = base_les.nu_t().iter().copied().fold(0.0_f64, f64::max);
        if step == 0 {
            max_nu_t_first = max_nu_t;
        }
        max_nu_t_last = max_nu_t;
        base.run(1);
        mirror.run(1);
    }
    assert!(
        max_nu_t_first > 0.0 || max_nu_t_last > 0.0,
        "T31 WALE did not produce a varying effective omega field: max_nu_t_first={max_nu_t_first:e}, max_nu_t_last={max_nu_t_last:e}"
    );

    let mut actual = Vec::with_capacity(3 * n * n);
    let mut reference = Vec::with_capacity(3 * n * n);
    let mut max_abs_dev = 0.0f64;
    for y in 0..n {
        for x in 0..n {
            let mx = n - 1 - x;
            let got_rho = mirror.rho(x, y, 0);
            let got_u = mirror.u(x, y, 0);
            let ref_rho = base.rho(mx, y, 0);
            let ref_u = reflected_vec(base.u(mx, y, 0));
            for (got, want) in [
                (got_rho, ref_rho),
                (got_u[0], ref_u[0]),
                (got_u[1], ref_u[1]),
            ] {
                max_abs_dev = max_abs_dev.max((got - want).abs());
                actual.push(got);
                reference.push(want);
            }
        }
    }
    let l2 = l2_rel(&actual, &reference);
    println!(
        "T31 WALE x-mirror equivariance: l2_rel={l2:.16e}, max_abs_dev={max_abs_dev:.16e}, max_nu_t_first={max_nu_t_first:.16e}, max_nu_t_last={max_nu_t_last:.16e}"
    );
    assert!(
        l2 <= 1.0e-12,
        "T31 A28 x-mirror deviation l2_rel={l2:.6e} > 1e-12, max_abs_dev={max_abs_dev:.6e}, max_nu_t_first={max_nu_t_first:.6e}, max_nu_t_last={max_nu_t_last:.6e}, denominator=mirrored [rho, ux, uy] L2 norm"
    );
}
