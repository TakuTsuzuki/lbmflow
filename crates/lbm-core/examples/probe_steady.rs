//! Dev probe: watch steady-state convergence + step timing for the tiny
//! Poiseuille grid (investigating BGK steady detection + rayon overhead).

use lbm_core::prelude::*;
use std::time::Instant;

fn main() {
    let mut sim: Simulation<f64> = SimConfig {
        nx: 4,
        ny: 10,
        nu: 0.1,
        collision: Collision::Bgk,
        edges: Edges {
            left: EdgeBC::Periodic,
            right: EdgeBC::Periodic,
            bottom: EdgeBC::BounceBack,
            top: EdgeBC::BounceBack,
        },
        force: [1e-6, 0.0],
        ..Default::default()
    }
    .build()
    .unwrap();

    let t0 = Instant::now();
    let mut prev: Vec<f64> = Vec::new();
    for it in 0..40 {
        sim.run(500);
        let cur: Vec<f64> = sim
            .ux_field()
            .iter()
            .chain(sim.uy_field())
            .copied()
            .collect();
        if !prev.is_empty() {
            let mut dmax = 0.0f64;
            let mut umax = 0.0f64;
            for (c, p) in cur.iter().zip(&prev) {
                dmax = dmax.max((c - p).abs());
                umax = umax.max(c.abs());
            }
            println!(
                "steps {:>6}  dmax {:.3e}  umax {:.3e}  dmax/umax {:.3e}",
                (it + 1) * 500,
                dmax,
                umax,
                dmax / umax
            );
        }
        prev = cur;
    }
    let dt = t0.elapsed();
    println!(
        "20000 steps took {:?} ({:.1} us/step on 4x10 grid)",
        dt,
        dt.as_micros() as f64 / 20_000.0
    );
}
