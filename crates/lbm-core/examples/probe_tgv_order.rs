//! Dev probe: print TGV convergence numbers for the PHYSICS.md record.

use lbm_core::prelude::*;
use std::f64::consts::PI;

fn tgv_l2(n: usize, collision: Collision, with_pressure: bool) -> f64 {
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
        let rho = if with_pressure {
            1.0 - 3.0 * u0 * u0 / 4.0 * ((2.0 * xf).cos() + (2.0 * yf).cos())
        } else {
            1.0
        };
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

fn main() {
    for with_p in [false, true] {
        let e32 = tgv_l2(32, Collision::default(), with_p);
        let e64 = tgv_l2(64, Collision::default(), with_p);
        let e128 = tgv_l2(128, Collision::default(), with_p);
        println!(
            "pressure_init={with_p}: e32={e32:.3e} e64={e64:.3e} e128={e128:.3e} order(32->64)={:.2} order(64->128)={:.2}",
            (e32 / e64).log2(),
            (e64 / e128).log2()
        );
    }
}
