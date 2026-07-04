//! Smoke test: Couette flow with a moving top wall. The steady profile is
//! linear and exact for half-way bounce-back at any tau; also sanity-checks
//! the momentum-exchange force probe against the analytical wall shear.

mod common;
use common::run_to_steady;

use lbm_core::prelude::*;

fn couette(tau: f64, collision: Collision) -> Simulation<f64> {
    let nu = (tau - 0.5) / 3.0;
    let mut sim: Simulation<f64> = SimConfig {
        nx: 4,
        ny: 10,
        nu,
        collision,
        edges: Edges {
            left: EdgeBC::Periodic,
            right: EdgeBC::Periodic,
            bottom: EdgeBC::BounceBack,
            top: EdgeBC::MovingWall { u: [0.1, 0.0] },
        },
        ..Default::default()
    }
    .build()
    .unwrap();
    assert!(
        run_to_steady(&mut sim, 500, 1e-11, 200_000),
        "did not reach steady state (tau = {tau})"
    );
    sim
}

#[test]
fn linear_profile_all_taus() {
    let u_wall = 0.1;
    for tau in [0.6, 1.0, 1.4] {
        for collision in [Collision::Bgk, Collision::default()] {
            let sim = couette(tau, collision);
            let h = (sim.ny() - 2) as f64;
            for j in 1..=(sim.ny() - 2) {
                let yw = j as f64 - 0.5;
                let expect = u_wall * yw / h;
                let got = sim.ux(0, j);
                assert!(
                    (got - expect).abs() < 1e-10 * u_wall,
                    "tau={tau} {collision:?} row {j}: {got} vs {expect}"
                );
            }
        }
    }
}

#[test]
fn wall_drag_matches_shear_stress() {
    let tau = 1.0;
    let nu = (tau - 0.5) / 3.0;
    let u_wall = 0.1;
    let mut sim = couette(tau, Collision::default());
    let ny = sim.ny();
    sim.set_force_probe(move |_, y| y == ny - 1);
    sim.run(10);
    let f = sim.probed_force();
    let h = (ny - 2) as f64;
    let expect = -nu * u_wall / h * sim.nx() as f64; // fluid resists lid motion
    let rel = (f[0] - expect).abs() / expect.abs();
    assert!(
        rel < 0.15,
        "lid drag = {}, expected ~{expect} (rel = {rel})",
        f[0]
    );
    assert!(f[0] < 0.0, "drag must oppose lid motion, got {}", f[0]);
}
