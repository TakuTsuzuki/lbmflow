//! Dev probe: (1) convective outflow pressure reflection vs zero-gradient,
//! (2) contact angle sweep via virtual wall density.

use lbm_core::multiphase::ShanChen;
use lbm_core::prelude::*;

/// Cylinder wake, measure near-outlet pressure RMS relative to mid-domain
/// (the T9 metric), comparing Outflow vs ConvectiveOutflow.
fn outlet_reflection(convective: bool) -> f64 {
    let (nx, ny) = (440, 164);
    let d = 40.0;
    let umax = 0.15;
    let right = if convective {
        EdgeBC::ConvectiveOutflow { u_conv: 0.1 }
    } else {
        EdgeBC::Outflow
    };
    let mut sim: Simulation<f64> = SimConfig {
        nx,
        ny,
        nu: 0.04,
        edges: Edges {
            left: EdgeBC::VelocityInlet { u: [0.0, 0.0] },
            right,
            bottom: EdgeBC::BounceBack,
            top: EdgeBC::BounceBack,
        },
        ..Default::default()
    }
    .build()
    .unwrap();
    let h = (ny - 2) as f64;
    sim.set_inlet_profile(Edge::Left, move |c| {
        if c == 0 || c as f64 >= h + 1.0 {
            return [0.0, 0.0];
        }
        let yw = c as f64 - 0.5;
        [4.0 * umax * yw * (h - yw) / (h * h), 0.0]
    });
    let (cx, cy, r) = (80.0, 80.0, d / 2.0);
    sim.set_solid_region(|x, y| {
        let dx = x as f64 - cx;
        let dy = y as f64 - cy;
        dx * dx + dy * dy <= r * r
    });
    sim.run(30_000); // develop shedding
                     // collect pressure time series at near-outlet and mid points
    let mut near = Vec::new();
    let mut mid = Vec::new();
    for _ in 0..4000 {
        sim.step();
        near.push(sim.rho(nx - 8, ny / 2));
        mid.push(sim.rho(nx / 2, ny / 4));
    }
    let rms = |v: &[f64]| {
        let m = v.iter().sum::<f64>() / v.len() as f64;
        (v.iter().map(|x| (x - m) * (x - m)).sum::<f64>() / v.len() as f64).sqrt()
    };
    let ratio = rms(&near) / rms(&mid).max(1e-30);
    println!(
        "[outlet {}] near_rms/mid_rms = {ratio:.2}",
        if convective {
            "convective"
        } else {
            "zero-grad "
        }
    );
    ratio
}

/// Contact angle via spherical cap (theta = 2 atan(2h/w)) for a droplet on
/// the bottom wall, sweeping the virtual wall density.
fn contact_angle(wall_rho: f64) -> f64 {
    let (nx, ny) = (160, 100);
    let mut sim: Simulation<f64> = SimConfig {
        nx,
        ny,
        nu: 1.0 / 6.0,
        edges: Edges {
            left: EdgeBC::Periodic,
            right: EdgeBC::Periodic,
            bottom: EdgeBC::BounceBack,
            top: EdgeBC::BounceBack,
        },
        ..Default::default()
    }
    .build()
    .unwrap();
    let (cx, cy, r0) = (nx as f64 / 2.0, 1.0, 22.0);
    sim.init_with(|x, y| {
        let d2 = (x as f64 - cx).powi(2) + (y as f64 - cy).powi(2);
        (if d2 < r0 * r0 { 2.0 } else { 0.15 }, 0.0, 0.0)
    });
    let sc = ShanChen::new(-5.0).with_wall_rho(wall_rho);
    for _ in 0..30_000 {
        sc.update_force(&mut sim);
        sim.step();
    }
    // half-density contour: contact width on the first fluid row + cap height
    let rho_mid = 1.0;
    let mut w = 0.0;
    for x in 0..nx {
        if sim.rho(x, 1) > rho_mid {
            w += 1.0;
        }
    }
    let mut hgt = 0.0;
    for y in 1..ny - 1 {
        if sim.rho(nx / 2, y) > rho_mid {
            hgt = y as f64 - 0.5; // wall at y = 0.5
        }
    }
    let theta = if w > 0.0 {
        2.0 * (2.0 * hgt / w).atan() * 180.0 / std::f64::consts::PI
    } else {
        180.0
    };
    println!("[angle wall_rho={wall_rho}] w={w:.0} h={hgt:.1} theta={theta:.0} deg");
    theta
}

fn main() {
    for wr in [1.6, 1.0, 0.6, 0.3] {
        contact_angle(wr);
    }
    outlet_reflection(false);
    outlet_reflection(true);
}
