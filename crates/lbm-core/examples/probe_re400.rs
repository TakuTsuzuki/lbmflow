//! Dev probe: Ghia Re=400 RMS vs lid speed (Mach dependence) and vs
//! convergence depth, to set the T7 threshold on evidence.

use lbm_core::prelude::*;

const N: usize = 129;
// Ghia, Ghia & Shin (1982), Re=400, u on vertical centreline (17 y-stations)
// and v on horizontal centreline (17 x-stations).
const Y_FRAC: [f64; 17] = [
    0.0000, 0.0547, 0.0625, 0.0703, 0.1016, 0.1719, 0.2813, 0.4531, 0.5000, 0.6172, 0.7344,
    0.8516, 0.9531, 0.9609, 0.9688, 0.9766, 1.0000,
];
const U_RE400: [f64; 17] = [
    0.00000, -0.08186, -0.09266, -0.10338, -0.14612, -0.24299, -0.32726, -0.17119, -0.11477,
    0.02135, 0.16256, 0.29093, 0.55892, 0.61756, 0.68439, 0.75837, 1.00000,
];
const X_FRAC: [f64; 17] = [
    0.0000, 0.0625, 0.0703, 0.0781, 0.0938, 0.1563, 0.2266, 0.2344, 0.5000, 0.8047, 0.8594,
    0.9063, 0.9453, 0.9531, 0.9609, 0.9688, 1.0000,
];
const V_RE400: [f64; 17] = [
    0.00000, 0.18360, 0.19713, 0.20920, 0.22965, 0.28124, 0.30203, 0.30174, 0.05186, -0.38598,
    -0.44993, -0.23827, -0.22847, -0.19254, -0.15663, -0.12146, 0.00000,
];

/// Sample a centreline velocity at wall-based fraction `frac`, treating the
/// walls as known sample points (u = 0 at stationary walls, u = U at the lid)
/// instead of clamping into the fluid.
fn sample_u(sim: &Simulation<f64>, u_lid: f64, frac: f64) -> f64 {
    let l = (N - 2) as f64;
    let pos = 0.5 + frac * l; // cell-index coordinate; walls at 0.5 and N-1.5
    let x = N / 2;
    // sample points: wall(0.5)=0, cells 1..=N-2, wall(N-1.5)=u_lid
    let grab = |j: f64| -> f64 {
        if j <= 0.5 {
            0.0
        } else if j >= (N - 1) as f64 - 0.5 {
            u_lid
        } else {
            sim.ux(x, j as usize)
        }
    };
    let j0 = pos.floor();
    let t = pos - j0;
    // piecewise-linear between the nearest sample points, including walls
    let a = if j0 < 1.0 { 0.5 } else { j0 };
    let b = if j0 + 1.0 > (N - 2) as f64 {
        (N - 1) as f64 - 0.5
    } else {
        j0 + 1.0
    };
    let (va, vb) = (grab(a), grab(b));
    let tt = if (b - a) > 0.0 { (pos - a) / (b - a) } else { 0.0 };
    let _ = t;
    va + (vb - va) * tt.clamp(0.0, 1.0)
}

fn sample_v(sim: &Simulation<f64>, frac: f64) -> f64 {
    let l = (N - 2) as f64;
    let pos = 0.5 + frac * l;
    let y = N / 2;
    let grab = |j: f64| -> f64 {
        if j <= 0.5 || j >= (N - 1) as f64 - 0.5 {
            0.0
        } else {
            sim.uy(j as usize, y)
        }
    };
    let j0 = pos.floor();
    let a = if j0 < 1.0 { 0.5 } else { j0 };
    let b = if j0 + 1.0 > (N - 2) as f64 {
        (N - 1) as f64 - 0.5
    } else {
        j0 + 1.0
    };
    let (va, vb) = (grab(a), grab(b));
    let tt = if (b - a) > 0.0 { (pos - a) / (b - a) } else { 0.0 };
    va + (vb - va) * tt.clamp(0.0, 1.0)
}

fn rms(sim: &Simulation<f64>, u: f64) -> f64 {
    let mut se = 0.0;
    let mut count = 0.0;
    for (yf, uref) in Y_FRAC.iter().zip(&U_RE400) {
        let val = sample_u(sim, u, *yf);
        se += (val / u - uref).powi(2);
        count += 1.0;
    }
    for (xf, vref) in X_FRAC.iter().zip(&V_RE400) {
        let val = sample_v(sim, *xf);
        se += (val / u - vref).powi(2);
        count += 1.0;
    }
    (se / count).sqrt()
}

fn run_case(u: f64, eps: f64, max_steps: usize) {
    let l = (N - 2) as f64;
    let mut sim: Simulation<f64> = SimConfig {
        nx: N,
        ny: N,
        nu: u * l / 400.0,
        edges: Edges {
            left: EdgeBC::BounceBack,
            right: EdgeBC::BounceBack,
            bottom: EdgeBC::BounceBack,
            top: EdgeBC::MovingWall { u: [u, 0.0] },
        },
        ..Default::default()
    }
    .build()
    .unwrap();
    let mut prev: Vec<f64> = Vec::new();
    let mut steps = 0;
    loop {
        sim.run(1000);
        steps += 1000;
        let cur: Vec<f64> = sim.ux_field().to_vec();
        if !prev.is_empty() {
            let dmax = cur
                .iter()
                .zip(&prev)
                .map(|(c, p)| (c - p).abs())
                .fold(0.0f64, f64::max);
            let umax = cur.iter().fold(0.0f64, |a, v| a.max(v.abs()));
            if dmax <= eps * umax || steps >= max_steps {
                break;
            }
        }
        prev = cur;
    }
    println!(
        "U={u} eps={eps:.0e}: steps={steps} RMS(in units of U)={:.4e}",
        rms(&sim, u)
    );
}

fn main() {
    run_case_verbose(0.1, 1e-8, 300_000);
}

fn run_case_verbose(u: f64, eps: f64, max_steps: usize) {
    let l = (N - 2) as f64;
    let mut sim: Simulation<f64> = SimConfig {
        nx: N,
        ny: N,
        nu: u * l / 400.0,
        edges: Edges {
            left: EdgeBC::BounceBack,
            right: EdgeBC::BounceBack,
            bottom: EdgeBC::BounceBack,
            top: EdgeBC::MovingWall { u: [u, 0.0] },
        },
        ..Default::default()
    }
    .build()
    .unwrap();
    let mut prev: Vec<f64> = Vec::new();
    let mut steps = 0;
    loop {
        sim.run(1000);
        steps += 1000;
        let cur: Vec<f64> = sim.ux_field().to_vec();
        if !prev.is_empty() {
            let dmax = cur
                .iter()
                .zip(&prev)
                .map(|(c, p)| (c - p).abs())
                .fold(0.0f64, f64::max);
            let umax = cur.iter().fold(0.0f64, |a, v| a.max(v.abs()));
            if dmax <= eps * umax || steps >= max_steps {
                break;
            }
        }
        prev = cur;
    }
    println!("steps = {steps}");
    println!("u-centreline (y_frac, sim/U, ghia, diff):");
    for (yf, uref) in Y_FRAC.iter().zip(&U_RE400) {
        let val = sample_u(&sim, u, *yf) / u;
        println!("  {yf:.4}  {val:+.5}  {uref:+.5}  {:+.4}", val - uref);
    }
    println!("v-centreline (x_frac, sim/U, ghia, diff):");
    for (xf, vref) in X_FRAC.iter().zip(&V_RE400) {
        let val = sample_v(&sim, *xf) / u;
        println!("  {xf:.4}  {val:+.5}  {vref:+.5}  {:+.4}", val - vref);
    }
}
