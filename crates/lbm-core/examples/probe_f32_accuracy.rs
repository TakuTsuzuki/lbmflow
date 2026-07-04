//! Dev probe: f32 accuracy with deviation storage — TGV L2 error vs f64.

use lbm_core::prelude::*;
use std::f64::consts::PI;

fn tgv_l2<T: Real>(n: usize) -> f64 {
    let nu = 0.02;
    let u0 = 1.28 / n as f64;
    let k = 2.0 * PI / n as f64;
    let mut sim: Simulation<T> = SimConfig {
        nx: n,
        ny: n,
        nu,
        ..Default::default()
    }
    .build()
    .unwrap();
    sim.init_with(|x, y| {
        let (xf, yf) = (k * x as f64, k * y as f64);
        let rho = 1.0 - 3.0 * u0 * u0 / 4.0 * ((2.0 * xf).cos() + (2.0 * yf).cos());
        (
            T::r(rho),
            T::r(-u0 * xf.cos() * yf.sin()),
            T::r(u0 * xf.sin() * yf.cos()),
        )
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
            num += (sim.ux(x, y).as_f64() - uxa).powi(2) + (sim.uy(x, y).as_f64() - uya).powi(2);
            den += uxa * uxa + uya * uya;
        }
    }
    (num / den).sqrt()
}

fn main() {
    println!("TGV L2rel  f64: N=64 {:.3e}", tgv_l2::<f64>(64));
    println!("TGV L2rel  f32: N=64 {:.3e}", tgv_l2::<f32>(64));
    println!("TGV L2rel  f64: N=128 {:.3e}", tgv_l2::<f64>(128));
    println!("TGV L2rel  f32: N=128 {:.3e}", tgv_l2::<f32>(128));
}
