//! Validation T1: Taylor-Green vortex decay, convergence, viscosity, and
//! 90-degree rotational symmetry.

mod common;

use common::l2_rel;
use lbm_core::prelude::*;
use std::f64::consts::PI;

#[derive(Clone, Copy)]
enum TgvMode {
    Base,
    Rot90,
}

fn tgv_init(n: usize, nu: f64, mode: TgvMode) -> (Simulation<f64>, usize, f64, f64) {
    let u0 = 1.28 / n as f64;
    let k = 2.0 * PI / n as f64;
    let mut sim: Simulation<f64> = SimConfig {
        nx: n,
        ny: n,
        nu,
        collision: Collision::default(),
        ..Default::default()
    }
    .build()
    .unwrap();
    sim.init_with(|x, y| {
        let (sx, sy, rotate_vec) = match mode {
            TgvMode::Base => (x, y, false),
            TgvMode::Rot90 => (y, (n - x) % n, true),
        };
        let (xf, yf) = (k * sx as f64, k * sy as f64);
        let rho = 1.0 - 3.0 * u0 * u0 / 4.0 * ((2.0 * xf).cos() + (2.0 * yf).cos());
        let ux = -u0 * xf.cos() * yf.sin();
        let uy = u0 * xf.sin() * yf.cos();
        if rotate_vec {
            (rho, -uy, ux)
        } else {
            (rho, ux, uy)
        }
    });
    let t_star = (1.0 / (2.0 * nu * k * k)).round() as usize;
    (sim, t_star, u0, k)
}

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
        let rho = 1.0 - 3.0 * u0 * u0 / 4.0 * ((2.0 * xf).cos() + (2.0 * yf).cos());
        (rho, -u0 * xf.cos() * yf.sin(), u0 * xf.sin() * yf.cos())
    });
    let t_star = (1.0 / (2.0 * nu * k * k)).round() as usize;
    sim.run(t_star);
    let decay = (-2.0 * nu * k * k * t_star as f64).exp();
    let mut actual = Vec::with_capacity(2 * n * n);
    let mut reference = Vec::with_capacity(2 * n * n);
    for y in 0..n {
        for x in 0..n {
            let (xf, yf) = (k * x as f64, k * y as f64);
            actual.push(sim.ux(x, y));
            actual.push(sim.uy(x, y));
            reference.push(-u0 * xf.cos() * yf.sin() * decay);
            reference.push(u0 * xf.sin() * yf.cos() * decay);
        }
    }
    l2_rel(&actual, &reference)
}

#[test]
fn t1_tgv_trt_accuracy_and_second_order_convergence() {
    let e32 = tgv_l2(32, Collision::default());
    let e64 = tgv_l2(64, Collision::default());
    let order = (e32 / e64).log2();
    assert!(e64 <= 1.5e-3, "T1 TRT N=64 L2rel = {e64:e}");
    assert!(
        order >= 1.7,
        "T1 TRT order = {order:e}, e32 = {e32:e}, e64 = {e64:e}"
    );
}

#[test]
fn t1_tgv_bgk_accuracy_is_comparable() {
    let e64 = tgv_l2(64, Collision::Bgk);
    assert!(e64 <= 2.0e-3, "T1 BGK N=64 L2rel = {e64:e}");
}

#[test]
fn t1_tgv_effective_viscosity_within_two_percent() {
    let n = 64;
    let nu = 0.02;
    let (mut sim, _, _, k) = tgv_init(n, nu, TgvMode::Base);
    let energy = |s: &Simulation<f64>| -> f64 {
        s.ux_field()
            .iter()
            .zip(s.uy_field())
            .map(|(ux, uy)| ux * ux + uy * uy)
            .sum()
    };
    sim.run(200);
    let e1 = energy(&sim);
    sim.run(1000);
    let e2 = energy(&sim);
    let nu_eff = (e1 / e2).ln() / (4.0 * k * k * 1000.0);
    let rel = (nu_eff / nu - 1.0).abs();
    assert!(
        rel <= 0.02,
        "T1 nu_eff rel = {rel:e}, nu_eff = {nu_eff:e}, nominal = {nu:e}"
    );
}

#[test]
fn t1_tgv_rotated_initial_field_stays_rotationally_symmetric() {
    let n = 64;
    let nu = 0.02;
    let (mut base, steps, _, _) = tgv_init(n, nu, TgvMode::Base);
    let (mut rot, _, _, _) = tgv_init(n, nu, TgvMode::Rot90);
    base.run(steps);
    rot.run(steps);
    let mut linf = 0.0f64;
    for y in 0..n {
        for x in 0..n {
            let sx = y;
            let sy = (n - x) % n;
            linf = linf.max((rot.ux(x, y) + base.uy(sx, sy)).abs());
            linf = linf.max((rot.uy(x, y) - base.ux(sx, sy)).abs());
        }
    }
    assert!(linf <= 1.0e-12, "T1 90deg rotation L_inf = {linf:e}");
}
