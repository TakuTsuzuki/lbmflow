//! Dev probe: Shan-Chen first physics — flat-interface coexistence,
//! pressure equality (SC EOS), spurious currents, and the Laplace law.

use lbm_core::multiphase::ShanChen;
use lbm_core::prelude::*;

fn flat_interface() {
    let (nx, ny) = (64, 128);
    let mut sim: Simulation<f64> = SimConfig {
        nx,
        ny,
        nu: 1.0 / 6.0, // tau = 1
        ..Default::default()
    }
    .build()
    .unwrap();
    sim.init_with(|_, y| {
        let liquid = y >= ny / 4 && y < 3 * ny / 4;
        (if liquid { 2.0 } else { 0.15 }, 0.0, 0.0)
    });
    let sc = ShanChen::new(-5.0);
    for _ in 0..30_000 {
        sc.update_force(&mut sim);
        sim.step();
    }
    let rho_l = sim.rho(nx / 2, ny / 2);
    let rho_v = sim.rho(nx / 2, 0);
    let p_l = sc.pressure(rho_l);
    let p_v = sc.pressure(rho_v);
    let umax = sim
        .ux_field()
        .iter()
        .chain(sim.uy_field())
        .fold(0.0f64, |a, v| a.max(v.abs()));
    // interface width: cells with 1.1*rho_v < rho < 0.9*rho_l along a column
    let mut width = 0;
    for y in 0..ny {
        let r = sim.rho(nx / 2, y);
        if r > rho_v * 1.5 && r < rho_l * 0.9 {
            width += 1;
        }
    }
    println!("[flat] rho_l={rho_l:.5} rho_v={rho_v:.5} ratio={:.1}", rho_l / rho_v);
    println!("[flat] p_l={p_l:.6} p_v={p_v:.6} |dp|/p={:.3e}", ((p_l - p_v) / p_l).abs());
    println!("[flat] max|u| (spurious) = {umax:.3e}, interface width ~ {} cells (2 interfaces)", width);
    println!("[flat] mass finite = {}", sim.total_mass().is_finite());
}

fn laplace() {
    println!("[laplace] r_init -> r_fit, dp, sigma=dp*r");
    let mut sigmas = Vec::new();
    for r0 in [12.0f64, 16.0, 20.0, 24.0] {
        let n = 128;
        let mut sim: Simulation<f64> = SimConfig {
            nx: n,
            ny: n,
            nu: 1.0 / 6.0,
            ..Default::default()
        }
        .build()
        .unwrap();
        let c = n as f64 / 2.0;
        sim.init_with(|x, y| {
            let d = ((x as f64 - c).powi(2) + (y as f64 - c).powi(2)).sqrt();
            (if d < r0 { 2.0 } else { 0.15 }, 0.0, 0.0)
        });
        let sc = ShanChen::new(-5.0);
        for _ in 0..40_000 {
            sc.update_force(&mut sim);
            sim.step();
        }
        let rho_in = sim.rho(n / 2, n / 2);
        let rho_out = sim.rho(2, 2);
        let dp = sc.pressure(rho_in) - sc.pressure(rho_out);
        // measure droplet radius: count cells above mean density -> area
        let rho_mid = 0.5 * (rho_in + rho_out);
        let area = sim
            .rho_field()
            .iter()
            .filter(|&&r| r > rho_mid)
            .count() as f64;
        let r_fit = (area / std::f64::consts::PI).sqrt();
        let sigma = dp * r_fit;
        println!(
            "  r0={r0:>4.0}  r_fit={r_fit:6.2}  dp={dp:.6e}  sigma={sigma:.6e}"
        );
        sigmas.push((1.0 / r_fit, dp));
    }
    // linear fit dp = sigma * (1/r): R^2
    let n = sigmas.len() as f64;
    let sx: f64 = sigmas.iter().map(|(x, _)| x).sum::<f64>() / n;
    let sy: f64 = sigmas.iter().map(|(_, y)| y).sum::<f64>() / n;
    let mut sxx = 0.0;
    let mut sxy = 0.0;
    let mut syy = 0.0;
    for (x, y) in &sigmas {
        sxx += (x - sx) * (x - sx);
        sxy += (x - sx) * (y - sy);
        syy += (y - sy) * (y - sy);
    }
    let slope = sxy / sxx;
    let r2 = (sxy * sxy) / (sxx * syy);
    println!("[laplace] sigma(slope) = {slope:.5e}, R^2 = {r2:.5}");
}

fn main() {
    flat_interface();
    laplace();
}
