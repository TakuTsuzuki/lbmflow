//! Smoke test: Taylor–Green vortex decay (periodic box).
//!
//! Diffusive scaling (u0 ∝ 1/N) so the compressibility error shrinks together
//! with the spatial error and the measured convergence order is ~2.

mod common;

use lbm_core::prelude::*;
use std::f64::consts::PI;

fn tgv_l2(n: usize, collision: Collision) -> f64 {
    let nu = 0.02;
    let u0 = 1.28 / n as f64;
    let k = 2.0 * PI / n as f64;
    let mut sim: Simulation<f64> = SimConfig {
        nx: n,
        ny: n,
        nu,
        collision,
        ..Default::default()
    }
    .build()
    .unwrap();
    sim.init_with(|x, y| {
        let (xf, yf) = (k * x as f64, k * y as f64);
        // Consistent initial state includes the analytic pressure field;
        // a flat density would radiate slowly-decaying acoustic waves and
        // pollute the solution at O(u0) (see docs/PHYSICS.md).
        let rho = 1.0 - 3.0 * u0 * u0 / 4.0 * ((2.0 * xf).cos() + (2.0 * yf).cos());
        (rho, -u0 * xf.cos() * yf.sin(), u0 * xf.sin() * yf.cos())
    });
    let t_star = (1.0 / (2.0 * nu * k * k)).round() as usize;
    sim.run(t_star);
    let decay = (-2.0 * nu * k * k * t_star as f64).exp();
    let mut num = 0.0;
    let mut den = 0.0;
    for y in 0..n {
        for x in 0..n {
            let (xf, yf) = (k * x as f64, k * y as f64);
            let uxa = -u0 * xf.cos() * yf.sin() * decay;
            let uya = u0 * xf.sin() * yf.cos() * decay;
            num += (sim.ux(x, y) - uxa).powi(2) + (sim.uy(x, y) - uya).powi(2);
            den += uxa * uxa + uya * uya;
        }
    }
    (num / den).sqrt()
}

#[test]
fn tgv_accuracy_and_convergence_trt() {
    let e32 = tgv_l2(32, Collision::default());
    let e64 = tgv_l2(64, Collision::default());
    assert!(e64 < 1.5e-3, "L2rel(64) = {e64}");
    let order = (e32 / e64).log2();
    assert!(order > 1.7, "order = {order} (e32 = {e32}, e64 = {e64})");
}

#[test]
fn tgv_accuracy_bgk() {
    let e64 = tgv_l2(64, Collision::Bgk);
    assert!(e64 < 2e-3, "L2rel(64) = {e64}");
}

#[test]
fn tgv_effective_viscosity() {
    let n = 64;
    let nu = 0.02;
    let u0 = 0.02;
    let k = 2.0 * PI / n as f64;
    let mut sim: Simulation<f64> = SimConfig {
        nx: n,
        ny: n,
        nu,
        ..Default::default()
    }
    .build()
    .unwrap();
    sim.init_with(|x, y| {
        let (xf, yf) = (k * x as f64, k * y as f64);
        (1.0, -u0 * xf.cos() * yf.sin(), u0 * xf.sin() * yf.cos())
    });
    let energy = |s: &Simulation<f64>| -> f64 {
        s.ux_field()
            .iter()
            .zip(s.uy_field())
            .map(|(ux, uy)| ux * ux + uy * uy)
            .sum()
    };
    let warmup = 200; // let the feq-only init settle
    let dt = 1000;
    sim.run(warmup);
    let e1 = energy(&sim);
    sim.run(dt);
    let e2 = energy(&sim);
    let nu_eff = (e1 / e2).ln() / (4.0 * k * k * dt as f64);
    let rel = (nu_eff / nu - 1.0).abs();
    assert!(rel < 0.02, "nu_eff = {nu_eff}, nominal = {nu}, rel = {rel}");
}
