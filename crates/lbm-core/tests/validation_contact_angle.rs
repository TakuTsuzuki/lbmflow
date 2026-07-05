//! Validation T11b: Shan-Chen wall adhesion contact-angle regression.

use lbm_core::multiphase::ShanChen;
use lbm_core::prelude::*;
use std::f64::consts::PI;

const G: f64 = -5.0;
const NU_TAU_1: f64 = 1.0 / 6.0;
const NX: usize = 96;
const NY: usize = 72;
const CONTACT_STEPS: usize = 20_000;

#[derive(Clone, Copy, Debug)]
struct ContactStats {
    g_wall: f64,
    theta_deg: f64,
    width: f64,
    height: f64,
}

fn run_wall_droplet(g_wall: f64, steps: usize) -> ContactStats {
    // The requested "bottom BounceBack, others Periodic" is not a valid
    // configuration in this engine because periodic boundaries must be paired
    // by axis. Use the valid bottom-wall setup: left/right periodic and a top
    // BounceBack rim far from the droplet.
    let mut sim: Simulation<f64> = SimConfig {
        nx: NX,
        ny: NY,
        nu: NU_TAU_1,
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
    let (cx, cy, r0) = (NX as f64 / 2.0, 15.0, 18.0);
    sim.init_with(|x, y| {
        if y == 0 || y + 1 == NY {
            return (0.15, 0.0, 0.0);
        }
        let d = ((x as f64 - cx).powi(2) + (y as f64 - cy).powi(2)).sqrt();
        (if d < r0 { 2.0 } else { 0.15 }, 0.0, 0.0)
    });
    let sc = ShanChen::new(G).with_wall(g_wall);
    for _ in 0..steps {
        sc.update_force(&mut sim);
        sim.step();
    }
    measure_contact_angle(&sim, g_wall)
}

fn measure_contact_angle(sim: &Simulation<f64>, g_wall: f64) -> ContactStats {
    let rho_l = sim.rho(NX / 2, 2).max(sim.rho(NX / 2, NY / 4));
    let rho_v = sim.rho(2, NY - 2);
    let rho_half = 0.5 * (rho_l + rho_v);

    // Spherical-cap fit from the half-density contour:
    // the contact width w is the span of columns whose half-density contour
    // reaches the first fluid row above the wall; the cap height h is the
    // distance from the physical wall (halfway between y=0 and y=1) to the
    // highest half-density cell. The fitted angle is theta = 2 atan(2h / w).
    let wall_y = 0.5;
    let mut min_y = NY;
    let mut max_y = 1usize;
    for x in 0..NX {
        for y in 1..NY - 1 {
            if sim.rho(x, y) > rho_half {
                min_y = min_y.min(y);
                max_y = max_y.max(y);
            }
        }
    }
    assert!(
        min_y < NY,
        "T11b no half-density liquid contour detected for G_w = {g_wall}, rho_half = {rho_half:e}, rho_l = {rho_l:e}, rho_v = {rho_v:e}"
    );

    let mut min_x = NX;
    let mut max_x = 0usize;
    for x in 0..NX {
        let touches_lowest_contour =
            (min_y..=(min_y + 1).min(NY - 2)).any(|y| sim.rho(x, y) > rho_half);
        if touches_lowest_contour {
            min_x = min_x.min(x);
            max_x = max_x.max(x);
        }
    }
    assert!(
        min_x <= max_x,
        "T11b no wall contact detected for G_w = {g_wall}, rho_half = {rho_half:e}, rho_l = {rho_l:e}, rho_v = {rho_v:e}"
    );
    let width = (max_x - min_x + 1) as f64;
    let height = max_y as f64 - wall_y;
    let theta_deg = (2.0 * (2.0 * height / width).atan()) * 180.0 / PI;
    ContactStats {
        g_wall,
        theta_deg,
        width,
        height,
    }
}

fn assert_angle(stats: ContactStats, expected: f64) {
    let err = (stats.theta_deg - expected).abs();
    assert!(
        err <= 8.0,
        "T11b G_w = {:+.1}, theta = {:.3} deg, frozen = {:.3} deg, err = {:.3} deg, width = {:.3}, height = {:.3}",
        stats.g_wall,
        stats.theta_deg,
        expected,
        err,
        stats.width,
        stats.height
    );
}

#[test]
fn t11b_wall_adhesion_contact_angles_are_monotone_and_frozen() {
    let wet = run_wall_droplet(-1.5, CONTACT_STEPS);
    let neutral = run_wall_droplet(0.0, CONTACT_STEPS);
    let dry = run_wall_droplet(1.5, CONTACT_STEPS);
    eprintln!(
        "T11b measured contact angles: G_w=-1.5 theta={:.3} width={:.3} height={:.3}; G_w=0 theta={:.3} width={:.3} height={:.3}; G_w=+1.5 theta={:.3} width={:.3} height={:.3}",
        wet.theta_deg,
        wet.width,
        wet.height,
        neutral.theta_deg,
        neutral.width,
        neutral.height,
        dry.theta_deg,
        dry.width,
        dry.height
    );

    assert!(
        wet.theta_deg < neutral.theta_deg && neutral.theta_deg < dry.theta_deg,
        "T11b contact angle monotonicity failed: G_w=-1.5 theta={:.3}, G_w=0 theta={:.3}, G_w=+1.5 theta={:.3}",
        wet.theta_deg,
        neutral.theta_deg,
        dry.theta_deg
    );

    // regression values measured 2026-07-05
    assert_angle(wet, 133.191);
    assert_angle(neutral, 160.435);
    assert_angle(dry, 163.740);
}
