// Inherited verbatim from the retired V1 suite at its retirement (2026-07-05,
// scripts/sync-tests.sh mechanical retarget); now the canonical facade tests.
//! Validation T4/T5: open boundary channels.
//!
//! T4 uses a prescribed parabolic inlet profile and checks bulk mass-flux
//! constancy away from the known Zou-He pressure-outlet boundary layer.

mod common;

use common::{l2_rel, run_to_steady};
use lbm_core::compat::lattice::CS2;
use lbm_core::compat::prelude::*;

#[derive(Clone, Copy, Debug)]
enum Orientation {
    LeftToRight,
    RightToLeft,
    BottomToTop,
    TopToBottom,
}

fn channel_edges(orientation: Orientation, u: f64, rho_out: f64) -> Edges<f64> {
    match orientation {
        Orientation::LeftToRight => Edges {
            left: EdgeBC::VelocityInlet { u: [u, 0.0] },
            right: EdgeBC::PressureOutlet { rho: rho_out },
            bottom: EdgeBC::BounceBack,
            top: EdgeBC::BounceBack,
        },
        Orientation::RightToLeft => Edges {
            left: EdgeBC::PressureOutlet { rho: rho_out },
            right: EdgeBC::VelocityInlet { u: [-u, 0.0] },
            bottom: EdgeBC::BounceBack,
            top: EdgeBC::BounceBack,
        },
        Orientation::BottomToTop => Edges {
            left: EdgeBC::BounceBack,
            right: EdgeBC::BounceBack,
            bottom: EdgeBC::VelocityInlet { u: [0.0, u] },
            top: EdgeBC::PressureOutlet { rho: rho_out },
        },
        Orientation::TopToBottom => Edges {
            left: EdgeBC::BounceBack,
            right: EdgeBC::BounceBack,
            bottom: EdgeBC::PressureOutlet { rho: rho_out },
            top: EdgeBC::VelocityInlet { u: [0.0, -u] },
        },
    }
}

fn inlet_edge(orientation: Orientation) -> Edge {
    match orientation {
        Orientation::LeftToRight => Edge::Left,
        Orientation::RightToLeft => Edge::Right,
        Orientation::BottomToTop => Edge::Bottom,
        Orientation::TopToBottom => Edge::Top,
    }
}

fn dims(orientation: Orientation) -> (usize, usize) {
    match orientation {
        Orientation::LeftToRight | Orientation::RightToLeft => (96, 34),
        Orientation::BottomToTop | Orientation::TopToBottom => (34, 96),
    }
}

fn channel_width(sim: &Simulation<f64>, orientation: Orientation) -> usize {
    match orientation {
        Orientation::LeftToRight | Orientation::RightToLeft => sim.ny() - 2,
        Orientation::BottomToTop | Orientation::TopToBottom => sim.nx() - 2,
    }
}

fn parabolic_speed(coord: usize, h: usize, umax: f64) -> f64 {
    if coord == 0 || coord == h + 1 {
        0.0
    } else {
        let h_f = h as f64;
        let yw = coord as f64 - 0.5;
        4.0 * umax * yw * (h_f - yw) / (h_f * h_f)
    }
}

fn inlet_velocity(orientation: Orientation, coord: usize, h: usize, umax: f64) -> [f64; 2] {
    let u = parabolic_speed(coord, h, umax);
    match orientation {
        Orientation::LeftToRight => [u, 0.0],
        Orientation::RightToLeft => [-u, 0.0],
        Orientation::BottomToTop => [0.0, u],
        Orientation::TopToBottom => [0.0, -u],
    }
}

fn normal_velocity(sim: &Simulation<f64>, orientation: Orientation, x: usize, y: usize) -> f64 {
    match orientation {
        Orientation::LeftToRight => sim.ux(x, y),
        Orientation::RightToLeft => -sim.ux(x, y),
        Orientation::BottomToTop => sim.uy(x, y),
        Orientation::TopToBottom => -sim.uy(x, y),
    }
}

fn section_flow(sim: &Simulation<f64>, orientation: Orientation, s: usize) -> f64 {
    match orientation {
        Orientation::LeftToRight | Orientation::RightToLeft => (1..=(sim.ny() - 2))
            .map(|y| sim.rho(s, y) * normal_velocity(sim, orientation, s, y))
            .sum(),
        Orientation::BottomToTop | Orientation::TopToBottom => (1..=(sim.nx() - 2))
            .map(|x| sim.rho(x, s) * normal_velocity(sim, orientation, x, s))
            .sum(),
    }
}

fn central_profile(sim: &Simulation<f64>, orientation: Orientation) -> Vec<f64> {
    match orientation {
        Orientation::LeftToRight | Orientation::RightToLeft => {
            let x = sim.nx() / 2;
            (1..=(sim.ny() - 2))
                .map(|y| normal_velocity(sim, orientation, x, y))
                .collect()
        }
        Orientation::BottomToTop | Orientation::TopToBottom => {
            let y = sim.ny() / 2;
            (1..=(sim.nx() - 2))
                .map(|x| normal_velocity(sim, orientation, x, y))
                .collect()
        }
    }
}

fn parabolic_profile(h: usize, umax: f64) -> Vec<f64> {
    (1..=h).map(|j| parabolic_speed(j, h, umax)).collect()
}

fn bulk_sections(sim: &Simulation<f64>, orientation: Orientation) -> Vec<usize> {
    match orientation {
        Orientation::LeftToRight => (1..=(sim.nx() - 25)).collect(),
        Orientation::RightToLeft => (24..=(sim.nx() - 2)).collect(),
        Orientation::BottomToTop => (1..=(sim.ny() - 25)).collect(),
        Orientation::TopToBottom => (24..=(sim.ny() - 2)).collect(),
    }
}

#[test]
fn t4_velocity_inlet_pressure_outlet_channel_all_four_orientations() {
    let umax = 0.05;
    for orientation in [
        Orientation::LeftToRight,
        Orientation::RightToLeft,
        Orientation::BottomToTop,
        Orientation::TopToBottom,
    ] {
        let (nx, ny) = dims(orientation);
        let mut sim: Simulation<f64> = SimConfig {
            nx,
            ny,
            nu: 0.02,
            collision: Collision::default(),
            edges: channel_edges(orientation, 0.0, 1.0),
            ..Default::default()
        }
        .build()
        .unwrap();
        let h = channel_width(&sim, orientation);
        sim.set_inlet_profile(inlet_edge(orientation), |c| {
            inlet_velocity(orientation, c, h, umax)
        });
        assert!(
            run_to_steady(&mut sim, 500, 1.0e-11, 160_000),
            "T4 steady = false, orientation = {orientation:?}, time = {}",
            sim.time()
        );
        let sections = bulk_sections(&sim, orientation);
        let flows: Vec<f64> = sections
            .iter()
            .copied()
            .map(|s| section_flow(&sim, orientation, s))
            .collect();
        let q_bar = flows.iter().sum::<f64>() / flows.len() as f64;
        let q_span = flows
            .iter()
            .map(|q| (q - q_bar).abs())
            .fold(0.0f64, f64::max)
            / q_bar.abs();
        assert!(
            q_span <= 1.0e-4,
            "T4 bulk mass-flux const rel = {q_span:e}, orientation = {orientation:?}, q_bar = {q_bar:e}, sections = {sections:?}, flows = {flows:?}"
        );
        let profile = central_profile(&sim, orientation);
        let reference = parabolic_profile(profile.len(), umax);
        let err = l2_rel(&profile, &reference);
        assert!(
            err <= 2.0e-3,
            "T4 central profile L2rel = {err:e}, orientation = {orientation:?}, profile = {profile:?}, reference = {reference:?}"
        );
        let m0 = sim.total_mass();
        sim.run(10_000);
        let drift = ((sim.total_mass() - m0) / m0).abs();
        assert!(
            drift <= 1.0e-11,
            "T4 steady mass drift = {drift:e}, orientation = {orientation:?}, time = {}",
            sim.time()
        );
    }
}

fn pressure_channel(delta_rho: f64) -> Simulation<f64> {
    let rho_mid = 1.0;
    let mut sim: Simulation<f64> = SimConfig {
        nx: 96,
        ny: 34,
        nu: 0.04,
        collision: Collision::default(),
        edges: Edges {
            left: EdgeBC::PressureOutlet {
                rho: rho_mid + 0.5 * delta_rho,
            },
            right: EdgeBC::PressureOutlet {
                rho: rho_mid - 0.5 * delta_rho,
            },
            bottom: EdgeBC::BounceBack,
            top: EdgeBC::BounceBack,
        },
        ..Default::default()
    }
    .build()
    .unwrap();
    assert!(
        run_to_steady(&mut sim, 500, 1.0e-11, 180_000),
        "T5 steady = false, delta_rho = {delta_rho:e}, time = {}",
        sim.time()
    );
    sim
}

fn pressure_channel_flow(sim: &Simulation<f64>) -> f64 {
    let x = sim.nx() / 2;
    (1..=(sim.ny() - 2)).map(|y| sim.ux(x, y)).sum()
}

fn pressure_linearity_r2(sim: &Simulation<f64>) -> f64 {
    let ys = sim.ny() / 2;
    let xs = 8..(sim.nx() - 8);
    let n = xs.len() as f64;
    let x_mean = xs.clone().map(|x| x as f64).sum::<f64>() / n;
    let p_mean = xs.clone().map(|x| CS2 * sim.rho(x, ys)).sum::<f64>() / n;
    let mut sxx = 0.0;
    let mut sxy = 0.0;
    let mut syy = 0.0;
    for x in xs.clone() {
        let dx = x as f64 - x_mean;
        let dp = CS2 * sim.rho(x, ys) - p_mean;
        sxx += dx * dx;
        sxy += dx * dp;
        syy += dp * dp;
    }
    let ss_res = xs
        .map(|x| {
            let fit = p_mean + sxy / sxx * (x as f64 - x_mean);
            let res = CS2 * sim.rho(x, ys) - fit;
            res * res
        })
        .sum::<f64>();
    1.0 - ss_res / syy
}

#[test]
fn t5_pressure_pressure_channel_flow_rate_and_linear_pressure() {
    let delta = 2.0e-3;
    let sim = pressure_channel(delta);
    let q = pressure_channel_flow(&sim);
    let h = sim.ny() - 2;
    let g = CS2 * delta / (sim.nx() - 1) as f64;
    let q_ref: f64 = (1..=h)
        .map(|j| {
            let yw = j as f64 - 0.5;
            g / (2.0 * sim.nu()) * yw * (h as f64 - yw)
        })
        .sum();
    let rel = ((q - q_ref) / q_ref).abs();
    assert!(
        rel <= 0.02,
        "T5 pressure-driven Q rel = {rel:e}, q = {q:e}, q_ref = {q_ref:e}"
    );
    let r2 = pressure_linearity_r2(&sim);
    assert!(r2 >= 0.999, "T5 pressure R2 = {r2:e}");
}

#[test]
fn t5_pressure_reversal_with_x_mirror_is_exact() {
    let fwd = pressure_channel(2.0e-3);
    let rev = pressure_channel(-2.0e-3);
    let mut linf = 0.0f64;
    for y in 1..=(fwd.ny() - 2) {
        for x in 0..fwd.nx() {
            let mx = fwd.nx() - 1 - x;
            linf = linf.max((rev.rho(x, y) - fwd.rho(mx, y)).abs());
            linf = linf.max((rev.ux(x, y) + fwd.ux(mx, y)).abs());
            linf = linf.max((rev.uy(x, y) - fwd.uy(mx, y)).abs());
        }
    }
    assert!(
        linf <= 1.0e-12,
        "T5 mirrored pressure reversal L_inf = {linf:e}"
    );
}

#[test]
fn t5_plain_pressure_reversal_is_approximately_antisymmetric() {
    let fwd = pressure_channel(2.0e-3);
    let rev = pressure_channel(-2.0e-3);
    let mut num = 0.0f64;
    let mut den = 0.0f64;
    for y in 1..=(fwd.ny() - 2) {
        for x in 0..fwd.nx() {
            num = num.max((fwd.ux(x, y) + rev.ux(x, y)).abs());
            num = num.max((fwd.uy(x, y) + rev.uy(x, y)).abs());
            den = den.max(fwd.ux(x, y).abs());
            den = den.max(fwd.uy(x, y).abs());
        }
    }
    let rel = num / den;
    assert!(
        rel <= 5.0e-3,
        "T5 plain pressure reversal relative L_inf = {rel:e}, abs = {num:e}, scale = {den:e}"
    );
}
