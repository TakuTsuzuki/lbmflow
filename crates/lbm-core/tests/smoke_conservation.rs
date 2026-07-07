// Inherited verbatim from the retired V1 suite at its retirement (2026-07-05,
// scripts/sync-tests.sh mechanical retarget); now the canonical facade tests.
//! Smoke test: exact conservation properties.

mod common;

use lbm_core::compat::prelude::*;
use std::f64::consts::PI;

fn smooth_init(n: usize) -> impl Fn(usize, usize) -> (f64, f64, f64) {
    let k = 2.0 * PI / n as f64;
    move |x, y| {
        let (xf, yf) = (k * x as f64, k * y as f64);
        (
            1.0 + 0.01 * (xf + 2.0 * yf).cos(),
            0.03 * yf.sin(),
            0.03 * (2.0 * xf).sin(),
        )
    }
}

#[test]
fn mass_is_conserved_periodic() {
    let n = 64;
    let mut sim: Simulation<f64> = SimConfig {
        nx: n,
        ny: n,
        nu: 0.02,
        ..Default::default()
    }
    .build()
    .unwrap();
    sim.init_with(smooth_init(n));
    let m0 = sim.total_mass();
    sim.run(1000);
    let drift = ((sim.total_mass() - m0) / m0).abs();
    assert!(drift < 1e-12, "mass drift = {drift:e}");
}

#[test]
fn mass_is_conserved_bounce_back_box() {
    let n = 48;
    let mut sim: Simulation<f64> = SimConfig {
        nx: n,
        ny: n,
        nu: 0.05,
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
    sim.init_with(smooth_init(n));
    let m0 = sim.total_mass();
    sim.run(1000);
    let drift = ((sim.total_mass() - m0) / m0).abs();
    assert!(drift < 1e-12, "mass drift = {drift:e}");
}

#[test]
fn mass_is_conserved_moving_wall() {
    let mut sim: Simulation<f64> = SimConfig {
        nx: 32,
        ny: 16,
        nu: 0.05,
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
    let m0 = sim.total_mass();
    sim.run(10_000);
    let drift = ((sim.total_mass() - m0) / m0).abs();
    assert!(drift < 1e-12, "mass drift = {drift:e}");
}

#[test]
fn uniform_force_adds_exact_momentum() {
    let n = 32;
    let force = [1e-6, 2e-6];
    let mut sim: Simulation<f64> = SimConfig {
        nx: n,
        ny: n,
        nu: 0.05,
        force,
        ..Default::default()
    }
    .build()
    .unwrap();
    let p0 = sim.total_momentum();
    let steps = 1000usize;
    sim.run(steps);
    let p1 = sim.total_momentum();
    let nf = sim.fluid_cell_count() as f64;
    for axis in 0..2 {
        let gained = p1[axis] - p0[axis];
        let expect = steps as f64 * nf * force[axis];
        let rel = ((gained - expect) / expect).abs();
        println!(
            "smoke uniform-force momentum axis {axis}: rel={rel:.3e}, gained={gained:.12e}, expected={expect:.12e}"
        );
        assert!(
            rel <= 7.0e-13,
            "axis {axis}: gained {gained}, expected {expect} (rel = {rel:e}, band = 7.0e-13)"
        );
    }
}
