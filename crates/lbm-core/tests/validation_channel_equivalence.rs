//! Radar row 21: force-driven and pressure-driven Poiseuille equivalence.

mod common;

use common::{l2_rel, run_to_steady};
use lbm_core::compat::prelude::*;

const NX: usize = 52;
const NY: usize = 34;
const H: f64 = (NY - 2) as f64;
const NU: f64 = 0.02;
const G: f64 = 1.0e-6;
const CS2: f64 = 1.0 / 3.0;
const UMAX: f64 = G * H * H / (8.0 * NU);

fn poiseuille_u(y: usize) -> f64 {
    let yw = y as f64 - 0.5;
    G / (2.0 * NU) * yw * (H - yw)
}

fn analytic_profile() -> Vec<f64> {
    (1..=(NY - 2)).map(poiseuille_u).collect()
}

fn profile_at(sim: &Simulation<f64>, x: usize) -> Vec<f64> {
    (1..=(NY - 2)).map(|y| sim.ux(x, y)).collect()
}

fn run_profile_to_steady(
    sim: &mut Simulation<f64>,
    x: usize,
    check_every: usize,
    tol: f64,
    min_steps: usize,
    max_steps: usize,
) -> (bool, f64) {
    let mut prev: Vec<f64> = Vec::new();
    let mut elapsed = 0;
    let mut last_rel = f64::INFINITY;
    while elapsed < max_steps {
        sim.run(check_every);
        elapsed += check_every;
        let cur = profile_at(sim, x);
        if !prev.is_empty() {
            let mut dmax = 0.0f64;
            let mut umax = 0.0f64;
            for (c, p) in cur.iter().zip(&prev) {
                dmax = dmax.max((c - p).abs());
                umax = umax.max(c.abs());
            }
            last_rel = if umax > 0.0 {
                dmax / umax
            } else {
                f64::INFINITY
            };
            if elapsed >= min_steps && umax > 0.0 && dmax <= tol * umax {
                return (true, last_rel);
            }
        }
        prev = cur;
    }
    (false, last_rel)
}

fn init_linear_pressure_poiseuille(sim: &mut Simulation<f64>, rho_out: f64) {
    let length = (NX - 1) as f64;
    sim.init_with(|x, y| {
        let rho = 1.0 + (rho_out - 1.0) * (x as f64 / length);
        let ux = if y == 0 || y == NY - 1 {
            0.0
        } else {
            poiseuille_u(y)
        };
        (rho, ux, 0.0)
    });
}

fn run_force_driven() -> Simulation<f64> {
    let mut sim: Simulation<f64> = SimConfig {
        nx: NX,
        ny: NY,
        nu: NU,
        collision: Collision::default(),
        edges: Edges {
            left: EdgeBC::Periodic,
            right: EdgeBC::Periodic,
            bottom: EdgeBC::BounceBack,
            top: EdgeBC::BounceBack,
        },
        force: [G, 0.0],
        ..Default::default()
    }
    .build()
    .unwrap();
    init_linear_pressure_poiseuille(&mut sim, 1.0);
    assert!(
        run_to_steady(&mut sim, 500, 1.0e-11, 300_000),
        "force-driven Poiseuille did not reach steady state, time={}",
        sim.time()
    );
    sim
}

fn run_pressure_driven() -> Simulation<f64> {
    // Continuous steady x-momentum for the body-force channel is
    // 0 = nu d2u/dy2 + g.  The pressure-driven channel satisfies
    // 0 = nu d2u/dy2 - (1/rho) dp/dx.  Matching profiles therefore requires
    // -(1/rho) dp/dx = g.  With the isothermal LBM equation of state
    // p = cs^2 rho and boundary nodes separated by L = nx - 1 (T5), the
    // matching density drop is Delta rho = g L rho / cs^2.  Taking rho_in=1
    // gives rho_out = 1 - g (nx - 1) / cs^2, and both drives have
    // u_max = g H^2 / (8 nu) at steady state.
    let rho_out = 1.0 - G * (NX - 1) as f64 / CS2;
    let mut sim: Simulation<f64> = SimConfig {
        nx: NX,
        ny: NY,
        nu: NU,
        collision: Collision::default(),
        edges: Edges {
            left: EdgeBC::VelocityInlet { u: [UMAX, 0.0] },
            right: EdgeBC::PressureOutlet { rho: rho_out },
            bottom: EdgeBC::BounceBack,
            top: EdgeBC::BounceBack,
        },
        ..Default::default()
    }
    .build()
    .unwrap();
    sim.set_inlet_profile(Edge::Left, |y| [poiseuille_u(y), 0.0]);
    let (steady, rel) = run_profile_to_steady(&mut sim, NX / 2, 500, 1.0e-7, 300_000, 300_000);
    assert!(
        steady,
        "pressure-driven Poiseuille mid-channel profile did not reach steady state, rel={rel:.9e}, time={}",
        sim.time()
    );
    sim
}

fn assert_centerline_shape(label: &str, profile: &[f64]) {
    let lower_wall = profile[0];
    let upper_wall = profile[profile.len() - 1];
    let c0 = profile[profile.len() / 2 - 1];
    let c1 = profile[profile.len() / 2];
    assert!(
        c0 > lower_wall && c1 > upper_wall,
        "{label}: center rows must exceed wall-adjacent rows, profile={profile:?}"
    );
    for j in 0..profile.len() / 2 {
        let diff = (profile[j] - profile[profile.len() - 1 - j]).abs();
        assert!(
            diff / UMAX <= 3.0e-3,
            "{label}: profile symmetry diff/umax={:.6e} at row pair {j}, profile={profile:?}",
            diff / UMAX
        );
    }
}

#[test]
fn force_driven_and_pressure_driven_poiseuille_profiles_are_equivalent() {
    let force = run_force_driven();
    let pressure = run_pressure_driven();
    let x_mid = NX / 2;
    let analytic = analytic_profile();
    let force_profile = profile_at(&force, x_mid);
    let pressure_profile = profile_at(&pressure, x_mid);

    let force_l2 = l2_rel(&force_profile, &analytic);
    let pressure_l2 = l2_rel(&pressure_profile, &analytic);
    let profile_l2 = l2_rel(&pressure_profile, &force_profile);
    let force_umax = force_profile
        .iter()
        .copied()
        .fold(f64::NEG_INFINITY, f64::max);
    let pressure_umax = pressure_profile
        .iter()
        .copied()
        .fold(f64::NEG_INFINITY, f64::max);
    let umax_rel = (force_umax - pressure_umax).abs() / UMAX;

    println!(
        "row21 Poiseuille equivalence: force_l2={force_l2:.9e}, \
         pressure_l2={pressure_l2:.9e}, umax_ref={UMAX:.9e}, \
         force_umax={force_umax:.9e}, pressure_umax={pressure_umax:.9e}, \
         umax_rel={umax_rel:.9e}, profile_l2={profile_l2:.9e}, \
         force_steps={}, pressure_steps={}",
        force.time(),
        pressure.time()
    );

    assert!(
        force_l2 <= 1.0e-4,
        "force-driven analytic profile L2rel={force_l2:.9e} > 1e-4"
    );
    assert!(
        pressure_l2 <= 1.0e-4,
        "pressure-driven analytic profile L2rel={pressure_l2:.9e} > 1e-4"
    );
    assert!(
        umax_rel <= 2.0e-2,
        "force/pressure centerline umax rel diff={umax_rel:.9e} > 2e-2, \
         force_umax={force_umax:.9e}, pressure_umax={pressure_umax:.9e}, umax_ref={UMAX:.9e}"
    );
    assert!(
        profile_l2 <= 3.0e-3,
        "force/pressure full-profile L2rel={profile_l2:.9e} > 3e-3"
    );
    assert_centerline_shape("force-driven", &force_profile);
    assert_centerline_shape("pressure-driven", &pressure_profile);
}
