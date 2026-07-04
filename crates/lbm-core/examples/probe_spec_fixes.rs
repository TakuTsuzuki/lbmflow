//! Dev probes to pin down spec numbers for the 4 TESTING_NOTES items.

use lbm_core::prelude::*;

fn probe_f32_momentum() {
    let n = 32;
    let force = [1e-5f32, 0.0];
    let mut sim: Simulation<f32> = SimConfig {
        nx: n,
        ny: n,
        nu: 0.05,
        force,
        ..Default::default()
    }
    .build()
    .unwrap();
    let p0 = sim.total_momentum();
    let steps = 100;
    sim.run(steps);
    let p1 = sim.total_momentum();
    let nf = sim.fluid_cell_count() as f64;
    let gained = (p1[0] - p0[0]) as f64;
    let expect = steps as f64 * nf * force[0] as f64;
    println!(
        "[f32 momentum] gained={gained:.7} expect={expect:.7} rel={:.3e}",
        ((gained - expect) / expect).abs()
    );
}

fn probe_mass_flux_constancy() {
    let (nx, ny) = (96, 34); // H = 32
    let h = (ny - 2) as f64;
    let umax = 0.05;
    let nu = 0.05;
    let mut sim: Simulation<f64> = SimConfig {
        nx,
        ny,
        nu,
        edges: Edges {
            left: EdgeBC::VelocityInlet { u: [0.0, 0.0] },
            right: EdgeBC::PressureOutlet { rho: 1.0 },
            bottom: EdgeBC::BounceBack,
            top: EdgeBC::BounceBack,
        },
        ..Default::default()
    }
    .build()
    .unwrap();
    sim.set_inlet_profile(Edge::Left, |y| {
        if y == 0 || y == ny - 1 {
            return [0.0, 0.0];
        }
        let yw = y as f64 - 0.5;
        [4.0 * umax * yw * (h - yw) / (h * h), 0.0]
    });
    // run to steady with the same criterion as the tests
    let mut prev: Vec<f64> = Vec::new();
    let mut steps = 0usize;
    loop {
        sim.run(2000);
        steps += 2000;
        let cur: Vec<f64> = sim.ux_field().to_vec();
        if !prev.is_empty() {
            let dmax = cur
                .iter()
                .zip(&prev)
                .map(|(c, p)| (c - p).abs())
                .fold(0.0f64, f64::max);
            let umax = cur.iter().fold(0.0f64, |a, v| a.max(v.abs()));
            if dmax <= 1e-11 * umax || steps >= 400_000 {
                println!("[mass flux] steady after {steps} steps (dmax/umax={:.1e})", dmax / umax);
                break;
            }
        }
        prev = cur;
    }
    // mass flux per column
    let mut qmin = f64::INFINITY;
    let mut qmax = f64::NEG_INFINITY;
    for x in 0..nx {
        let mut q = 0.0;
        for y in 1..ny - 1 {
            q += sim.rho(x, y) * sim.ux(x, y);
        }
        qmin = qmin.min(q);
        qmax = qmax.max(q);
    }
    let qbar = 0.5 * (qmin + qmax);
    println!(
        "[mass flux] qmin={qmin:.8} qmax={qmax:.8} spread_rel={:.3e}",
        (qmax - qmin) / qbar
    );
    // also volume flux spread for comparison
    let mut vmin = f64::INFINITY;
    let mut vmax = f64::NEG_INFINITY;
    for x in 0..nx {
        let q: f64 = (1..ny - 1).map(|y| sim.ux(x, y)).sum();
        vmin = vmin.min(q);
        vmax = vmax.max(q);
    }
    println!(
        "[volume flux] spread_rel={:.3e}",
        (vmax - vmin) / (0.5 * (vmin + vmax))
    );
}

fn probe_cavity_stability() {
    for (u, magic, label) in [
        (0.1, 3.0 / 16.0, "U=0.10 L=3/16"),
        (0.1, 0.25, "U=0.10 L=1/4 "),
        (0.05, 3.0 / 16.0, "U=0.05 L=3/16"),
        (0.05, 0.25, "U=0.05 L=1/4 "),
    ] {
        let n = 128;
        let nu = (0.51 - 0.5) / 3.0;
        let mut sim: Simulation<f64> = SimConfig {
            nx: n,
            ny: n,
            nu,
            collision: Collision::Trt { magic },
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
        let mut nan_at = None;
        for k in 0..20 {
            sim.run(500);
            let bad = sim
                .ux_field()
                .iter()
                .chain(sim.uy_field())
                .any(|v| !v.is_finite());
            if bad {
                nan_at = Some((k + 1) * 500);
                break;
            }
        }
        let re = u * (n as f64 - 2.0) / nu;
        match nan_at {
            Some(t) => println!("[cavity tau=0.51 {label}] Re={re:.0} NaN at step {t}"),
            None => {
                let umax = sim
                    .ux_field()
                    .iter()
                    .chain(sim.uy_field())
                    .fold(0.0f64, |a, v| a.max(v.abs()));
                println!("[cavity tau=0.51 {label}] Re={re:.0} OK (10k steps, max|u|={umax:.3})");
            }
        }
    }
}

fn main() {
    probe_f32_momentum();
    probe_mass_flux_constancy();
    probe_cavity_stability();
}
