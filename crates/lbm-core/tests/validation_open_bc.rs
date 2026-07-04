//! Validation T4/T5: open boundary channels.
//!
//! The current public API prescribes one velocity vector per inlet edge, not a
//! spatially varying parabolic profile. T4 therefore uses the available
//! uniform velocity inlet and checks the spec-observable consequences:
//! orientation coverage, flow-rate constancy, outlet flux balance, and
//! developed central profile against the Poiseuille profile with matching
//! flow rate.

mod common;

use common::{l2_rel, run_to_steady};
use lbm_core::lattice::CS2;
use lbm_core::prelude::*;

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

fn dims(orientation: Orientation) -> (usize, usize) {
    match orientation {
        Orientation::LeftToRight | Orientation::RightToLeft => (96, 34),
        Orientation::BottomToTop | Orientation::TopToBottom => (34, 96),
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
            .map(|y| normal_velocity(sim, orientation, s, y))
            .sum(),
        Orientation::BottomToTop | Orientation::TopToBottom => (1..=(sim.nx() - 2))
            .map(|x| normal_velocity(sim, orientation, x, s))
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

fn poiseuille_from_flux(h: usize, q: f64) -> Vec<f64> {
    let h_f = h as f64;
    let shape: Vec<f64> = (1..=h)
        .map(|j| {
            let yw = j as f64 - 0.5;
            yw * (h_f - yw)
        })
        .collect();
    let scale = q / shape.iter().sum::<f64>();
    shape.into_iter().map(|v| scale * v).collect()
}

#[test]
fn t4_velocity_inlet_pressure_outlet_channel_all_four_orientations() {
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
            edges: channel_edges(orientation, 0.05, 1.0),
            ..Default::default()
        }
        .build()
        .unwrap();
        assert!(
            run_to_steady(&mut sim, 500, 1.0e-11, 160_000),
            "T4 steady = false, orientation = {orientation:?}, time = {}",
            sim.time()
        );
        let sections = match orientation {
            Orientation::LeftToRight | Orientation::RightToLeft => 1..(sim.nx() - 1),
            Orientation::BottomToTop | Orientation::TopToBottom => 1..(sim.ny() - 1),
        };
        let flows: Vec<f64> = sections
            .map(|s| section_flow(&sim, orientation, s))
            .collect();
        let q_bar = flows.iter().sum::<f64>() / flows.len() as f64;
        let q_span = flows
            .iter()
            .map(|q| (q - q_bar).abs())
            .fold(0.0f64, f64::max)
            / q_bar.abs();
        assert!(
            q_span <= 1.0e-6,
            "T4 flow const rel = {q_span:e}, orientation = {orientation:?}, q_bar = {q_bar:e}, flows = {flows:?}"
        );
        let q_in = flows[0];
        let q_out = flows[flows.len() - 1];
        let flux_rel = ((q_in - q_out) / q_bar).abs();
        assert!(
            flux_rel <= 1.0e-6,
            "T4 inlet/outlet flux rel = {flux_rel:e}, orientation = {orientation:?}, q_in = {q_in:e}, q_out = {q_out:e}, q_bar = {q_bar:e}"
        );
        let profile = central_profile(&sim, orientation);
        let reference = poiseuille_from_flux(profile.len(), q_bar);
        let err = l2_rel(&profile, &reference);
        assert!(
            err <= 2.0e-3,
            "T4 central profile L2rel = {err:e}, orientation = {orientation:?}, q_bar = {q_bar:e}, profile = {profile:?}, reference = {reference:?}"
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
    let n = sim.nx() as f64;
    let x_mean = (sim.nx() - 1) as f64 / 2.0;
    let p_mean = (0..sim.nx()).map(|x| CS2 * sim.rho(x, ys)).sum::<f64>() / n;
    let mut sxx = 0.0;
    let mut sxy = 0.0;
    let mut syy = 0.0;
    for x in 0..sim.nx() {
        let dx = x as f64 - x_mean;
        let dp = CS2 * sim.rho(x, ys) - p_mean;
        sxx += dx * dx;
        sxy += dx * dp;
        syy += dp * dp;
    }
    let ss_res = (0..sim.nx())
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
fn t5_pressure_sign_reversal_is_antisymmetric() {
    let fwd = pressure_channel(2.0e-3);
    let rev = pressure_channel(-2.0e-3);
    let mut linf = 0.0f64;
    for y in 1..=(fwd.ny() - 2) {
        for x in 0..fwd.nx() {
            linf = linf.max((fwd.ux(x, y) + rev.ux(x, y)).abs());
            linf = linf.max((fwd.uy(x, y) + rev.uy(x, y)).abs());
        }
    }
    assert!(linf <= 1.0e-12, "T5 sign reversal L_inf = {linf:e}");
}
