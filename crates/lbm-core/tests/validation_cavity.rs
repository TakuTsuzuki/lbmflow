//! Validation T7: lid-driven cavity against Ghia et al. centreline data.

mod common;

use common::run_to_steady;
use lbm_core::prelude::*;

const N: usize = 129;
const L: f64 = (N - 2) as f64;
const U: f64 = 0.1;
const MAGIC: f64 = 3.0 / 16.0;

// Ghia, U., Ghia, K. N., & Shin, C. T. (1982).
// High-Re solutions for incompressible flow using the Navier-Stokes equations
// and a multigrid method. Journal of Computational Physics, 48(3), 387-411.
const GHIA_Y: [f64; 17] = [
    1.0000, 0.9766, 0.9688, 0.9609, 0.9531, 0.8516, 0.7344, 0.6172, 0.5000, 0.4531, 0.2813, 0.1719,
    0.1016, 0.0703, 0.0625, 0.0547, 0.0000,
];
const GHIA_X: [f64; 17] = [
    1.0000, 0.9688, 0.9609, 0.9531, 0.9453, 0.9063, 0.8594, 0.8047, 0.5000, 0.2344, 0.2266, 0.1563,
    0.0938, 0.0781, 0.0703, 0.0625, 0.0000,
];
const U_RE100: [f64; 17] = [
    1.00000, 0.84123, 0.78871, 0.73722, 0.68717, 0.23151, 0.00332, -0.13641, -0.20581, -0.21090,
    -0.15662, -0.10150, -0.06434, -0.04775, -0.04192, -0.03717, 0.00000,
];
const V_RE100: [f64; 17] = [
    0.00000, -0.05906, -0.07391, -0.08864, -0.10313, -0.16914, -0.22445, -0.24533, 0.05454,
    0.17527, 0.17507, 0.16077, 0.12317, 0.10890, 0.10091, 0.09233, 0.00000,
];
const U_RE400: [f64; 17] = [
    1.00000, 0.75837, 0.68439, 0.61756, 0.55892, 0.29093, 0.16256, 0.02135, -0.11477, -0.17119,
    -0.32726, -0.24299, -0.14612, -0.10338, -0.09266, -0.08186, 0.00000,
];
const V_RE400: [f64; 17] = [
    0.00000, -0.12146, -0.15663, -0.19254, -0.22847, -0.23827, -0.44993, -0.38598, 0.05186,
    0.30174, 0.30203, 0.28124, 0.22965, 0.20920, 0.19713, 0.18360, 0.00000,
];
const U_RE1000: [f64; 17] = [
    1.00000, 0.65928, 0.57492, 0.51117, 0.46604, 0.33304, 0.18719, 0.05702, -0.06080, -0.10648,
    -0.27805, -0.38289, -0.29730, -0.22220, -0.20196, -0.18109, 0.00000,
];
const V_RE1000: [f64; 17] = [
    0.00000, -0.21388, -0.27669, -0.33714, -0.39188, -0.51550, -0.42665, -0.31966, 0.02526,
    0.32235, 0.33075, 0.37095, 0.32627, 0.30353, 0.29012, 0.27485, 0.00000,
];

#[derive(Clone, Copy, Debug)]
enum Lid {
    Top,
    Left,
    Bottom,
    Right,
}

fn cavity(re: f64, lid: Lid) -> Simulation<f64> {
    let nu = U * L / re;
    let edges = match lid {
        Lid::Top => Edges {
            left: EdgeBC::BounceBack,
            right: EdgeBC::BounceBack,
            bottom: EdgeBC::BounceBack,
            top: EdgeBC::MovingWall { u: [U, 0.0] },
        },
        Lid::Left => Edges {
            left: EdgeBC::MovingWall { u: [0.0, -U] },
            right: EdgeBC::BounceBack,
            bottom: EdgeBC::BounceBack,
            top: EdgeBC::BounceBack,
        },
        Lid::Bottom => Edges {
            left: EdgeBC::BounceBack,
            right: EdgeBC::BounceBack,
            bottom: EdgeBC::MovingWall { u: [-U, 0.0] },
            top: EdgeBC::BounceBack,
        },
        Lid::Right => Edges {
            left: EdgeBC::BounceBack,
            right: EdgeBC::MovingWall { u: [0.0, U] },
            bottom: EdgeBC::BounceBack,
            top: EdgeBC::BounceBack,
        },
    };
    SimConfig {
        nx: N,
        ny: N,
        nu,
        collision: Collision::Trt { magic: MAGIC },
        edges,
        ..Default::default()
    }
    .build()
    .unwrap()
}

fn run_cavity_to_spec_limit(sim: &mut Simulation<f64>) -> bool {
    run_to_steady(sim, 1_000, 1.0e-8, 200_000)
}

fn sample_centerline_u(sim: &Simulation<f64>, y_frac: f64) -> f64 {
    let x = N / 2;
    let pos = 0.5 + y_frac * L;
    let y0 = pos.floor().clamp(1.0, (N - 2) as f64) as usize;
    let y1 = (y0 + 1).min(N - 2);
    let t = pos - y0 as f64;
    (1.0 - t) * sim.ux(x, y0) + t * sim.ux(x, y1)
}

fn sample_centerline_v(sim: &Simulation<f64>, x_frac: f64) -> f64 {
    let y = N / 2;
    let pos = 0.5 + x_frac * L;
    let x0 = pos.floor().clamp(1.0, (N - 2) as f64) as usize;
    let x1 = (x0 + 1).min(N - 2);
    let t = pos - x0 as f64;
    (1.0 - t) * sim.uy(x0, y) + t * sim.uy(x1, y)
}

fn rms_centerline_error(sim: &Simulation<f64>, u_ref: &[f64; 17], v_ref: &[f64; 17]) -> f64 {
    let mut sum = 0.0;
    let mut n = 0usize;
    for i in 0..17 {
        let du = sample_centerline_u(sim, GHIA_Y[i]) - U * u_ref[i];
        let dv = sample_centerline_v(sim, GHIA_X[i]) - U * v_ref[i];
        sum += du * du + dv * dv;
        n += 2;
    }
    (sum / n as f64).sqrt()
}

fn primary_vortex_center(sim: &Simulation<f64>) -> (f64, f64, f64) {
    let mut best = (0usize, 0usize, 0.0f64);
    for x in 1..=(N - 2) {
        let mut psi = 0.0;
        for y in 1..=(N - 2) {
            psi += sim.ux(x, y);
            if psi.abs() > best.2.abs() {
                best = (x, y, psi);
            }
        }
    }
    ((best.0 as f64 - 0.5) / L, (best.1 as f64 - 0.5) / L, best.2)
}

fn assert_ghia_case(
    re: f64,
    u_ref: &[f64; 17],
    v_ref: &[f64; 17],
    rms_limit: f64,
    center: (f64, f64),
) {
    let mut sim = cavity(re, Lid::Top);
    let steady = run_cavity_to_spec_limit(&mut sim);
    let rms = rms_centerline_error(&sim, u_ref, v_ref);
    assert!(
        rms <= rms_limit,
        "T7 Re={re} Ghia RMS = {rms:e}, limit = {rms_limit:e}, steady = {steady}, time = {}",
        sim.time()
    );
    let (cx, cy, psi) = primary_vortex_center(&sim);
    assert!(
        (cx - center.0).abs() <= 0.02 && (cy - center.1).abs() <= 0.02,
        "T7 Re={re} vortex center = ({cx:e}, {cy:e}), expected = ({:e}, {:e}), psi = {psi:e}, steady = {steady}, time = {}",
        center.0,
        center.1,
        sim.time()
    );
}

#[test]
fn t7_lid_driven_cavity_re100_matches_ghia() {
    assert_ghia_case(100.0, &U_RE100, &V_RE100, 0.02 * U, (0.6172, 0.7344));
}

#[test]
fn t7_lid_driven_cavity_re400_matches_ghia() {
    assert_ghia_case(400.0, &U_RE400, &V_RE400, 0.02 * U, (0.5547, 0.6055));
}

#[test]
#[ignore]
fn t7_lid_driven_cavity_re1000_matches_ghia() {
    assert_ghia_case(1000.0, &U_RE1000, &V_RE1000, 0.03 * U, (0.5313, 0.5625));
}

fn mapped_velocity(sim: &Simulation<f64>, lid: Lid, x: usize, y: usize) -> (f64, f64) {
    match lid {
        Lid::Top => (sim.ux(x, y), sim.uy(x, y)),
        Lid::Left => {
            let xr = N - 1 - y;
            let yr = x;
            (-sim.uy(xr, yr), sim.ux(xr, yr))
        }
        Lid::Bottom => {
            let xr = N - 1 - x;
            let yr = N - 1 - y;
            (-sim.ux(xr, yr), -sim.uy(xr, yr))
        }
        Lid::Right => {
            let xr = y;
            let yr = N - 1 - x;
            (sim.uy(xr, yr), -sim.ux(xr, yr))
        }
    }
}

#[test]
fn t7_re100_cavity_is_exact_under_four_lid_orientations() {
    let steps = 2_000;
    let mut sims = [
        (Lid::Top, cavity(100.0, Lid::Top)),
        (Lid::Left, cavity(100.0, Lid::Left)),
        (Lid::Bottom, cavity(100.0, Lid::Bottom)),
        (Lid::Right, cavity(100.0, Lid::Right)),
    ];
    for (_, sim) in &mut sims {
        sim.run(steps);
    }
    let base = &sims[0].1;
    for (lid, sim) in &sims[1..] {
        let mut linf = 0.0f64;
        for y in 1..=(N - 2) {
            for x in 1..=(N - 2) {
                let (ux, uy) = mapped_velocity(sim, *lid, x, y);
                linf = linf.max((ux - base.ux(x, y)).abs());
                linf = linf.max((uy - base.uy(x, y)).abs());
            }
        }
        assert!(
            linf <= 1.0e-10,
            "T7 Re=100 orientation {lid:?} mapped L_inf = {linf:e}, steps = {steps}"
        );
    }
}
