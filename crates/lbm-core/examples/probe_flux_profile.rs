//! Dev probe: where does the mass-flux spread live? Print Q(x) along the
//! channel + check global mass balance at steady state.

use lbm_core::prelude::*;

fn main() {
    for (collision, label) in [
        (Collision::Trt { magic: 3.0 / 16.0 }, "TRT 3/16"),
        (Collision::Trt { magic: 0.25 }, "TRT 1/4"),
        (Collision::Bgk, "BGK"),
    ] {
        println!("=== {label} ===");
        run_case(collision);
    }
}

fn run_case(collision: Collision) {
    let (nx, ny) = (96, 34);
    let h = (ny - 2) as f64;
    let umax = 0.05;
    let mut sim: Simulation<f64> = SimConfig {
        nx,
        ny,
        nu: 0.05,
        collision,
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
    sim.run(150_000);

    let q = |x: usize| -> f64 { (1..ny - 1).map(|y| sim.rho(x, y) * sim.ux(x, y)).sum() };
    // exact discrete face flux between columns x and x+1 is not directly
    // accessible via public API; use moments per column instead.
    print!("Q(x): ");
    for x in [0usize, 1, 2, 3, 8, 24, 48, 72, 92, 93, 94, 95] {
        print!("x{x}={:.6} ", q(x));
    }
    println!();

    // bulk constancy (exclude 6 columns next to the pressure outlet)
    let bulk: Vec<f64> = (0..nx - 6).map(q).collect();
    let (bmin, bmax) = bulk
        .iter()
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(a, b), &v| (a.min(v), b.max(v)));
    println!(
        "bulk (x<{}) spread_rel={:.3e}",
        nx - 6,
        (bmax - bmin) / (0.5 * (bmax + bmin))
    );
}
