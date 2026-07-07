//! ACC INIT-RINGING: coverage-gap radar #20, pitfall #3.
//!
//! The T1 Taylor-Green vortex must be initialized with the matching O(u0^2)
//! pressure field, not flat density. This test quantifies the acoustic residual
//! that a flat-density TGV launches and pins the pressure-consistent reduction.

use lbm_core::compat::prelude::*;
use std::f64::consts::PI;

const N: usize = 64;
const U0: f64 = 0.05;
const NU: f64 = 0.02;
const STEPS: usize = 200;
const TRT_MAGIC: f64 = 3.0 / 16.0;

#[derive(Clone, Copy)]
enum InitDensity {
    Flat,
    PressureConsistent,
}

fn make_tgv(density: InitDensity) -> Simulation<f64> {
    let k = 2.0 * PI / N as f64;
    let mut sim: Simulation<f64> = SimConfig {
        nx: N,
        ny: N,
        nu: NU,
        collision: Collision::Trt { magic: TRT_MAGIC },
        ..Default::default()
    }
    .build()
    .unwrap();
    sim.init_with(|x, y| {
        let xf = k * x as f64;
        let yf = k * y as f64;
        let rho = match density {
            InitDensity::Flat => 1.0,
            InitDensity::PressureConsistent => {
                1.0 - 3.0 * U0 * U0 / 4.0 * ((2.0 * xf).cos() + (2.0 * yf).cos())
            }
        };
        (rho, -U0 * xf.cos() * yf.sin(), U0 * xf.sin() * yf.cos())
    });
    sim
}

fn density_mode_coefficient(sim: &Simulation<f64>) -> f64 {
    let k2 = 4.0 * PI / N as f64;
    let mut num = 0.0;
    let mut den = 0.0;
    for y in 0..N {
        let cy = (k2 * y as f64).cos();
        for x in 0..N {
            let phi = (k2 * x as f64).cos() + cy;
            num += (sim.rho(x, y) - 1.0) * phi;
            den += phi * phi;
        }
    }
    num / den
}

fn pressure_consistent_mode_at_step(step: usize) -> f64 {
    let k = 2.0 * PI / N as f64;
    // TGV is divergence-free at O(u0), so the incompressible pressure follows
    // from the Poisson equation
    //
    //     ∇²p = -ρ0 ∂i uj ∂j ui.
    //
    // For ux = -u0 cos(kx) sin(ky), uy = u0 sin(kx) cos(ky),
    // ∂i uj ∂j ui = 2 u0² k² sin²(kx) sin²(ky)
    //                 - 2 u0² k² cos²(kx) cos²(ky)
    //               = -u0² k² * (cos 2kx + cos 2ky).
    // Since ∇²[cos 2kx + cos 2ky] = -4 k²[cos 2kx + cos 2ky],
    // p = -ρ0 u0² / 4 * [cos 2kx + cos 2ky]. In the weakly-compressible LBM
    // equation of state p = cs² rho' with cs² = 1/3, so
    //
    //     rho' = 3p = -(3u0²/4) [cos 2kx + cos 2ky].
    //
    // The mode decays as u(t)^2 under viscous TGV decay, so the pressure
    // coefficient decays as exp(-4 nu k² t).
    -3.0 * U0 * U0 / 4.0 * (-4.0 * NU * k * k * step as f64).exp()
}

fn max_acoustic_residual(mut sim: Simulation<f64>) -> (f64, f64) {
    let mut max_residual = 0.0f64;
    let mut max_raw = 0.0f64;
    for step in 0..=STEPS {
        let raw = density_mode_coefficient(&sim);
        let residual = raw - pressure_consistent_mode_at_step(step);
        max_raw = max_raw.max(raw.abs());
        max_residual = max_residual.max(residual.abs());
        if step < STEPS {
            sim.step();
        }
    }
    (max_residual, max_raw)
}

#[test]
fn pressure_consistent_tgv_init_reduces_acoustic_ringing() {
    let (flat, flat_raw) = max_acoustic_residual(make_tgv(InitDensity::Flat));
    let (consistent, consistent_raw) =
        max_acoustic_residual(make_tgv(InitDensity::PressureConsistent));
    let ratio = flat / consistent;
    let u0_sq = U0 * U0;

    println!(
        "ACC INIT-RINGING radar20 pitfall3: N={N} u0={U0:.6e} nu={NU:.6e} steps={STEPS} \
         flat_residual_max={flat:.6e} consistent_residual_max={consistent:.6e} \
         ratio={ratio:.6e} flat_raw_mode_max={flat_raw:.6e} \
         consistent_raw_mode_max={consistent_raw:.6e} u0_sq={u0_sq:.6e}"
    );
    assert!(
        flat > 0.5 * u0_sq,
        "flat-density acoustic residual max={flat:.6e} must be > 0.5*u0^2={:.6e}; \
         denominator is u0^2={u0_sq:.6e}",
        0.5 * u0_sq
    );
    assert!(
        consistent < 0.02 * u0_sq,
        "pressure-consistent acoustic residual max={consistent:.6e} must be < \
         0.02*u0^2={:.6e}; denominator is u0^2={u0_sq:.6e}",
        0.02 * u0_sq
    );
    assert!(
        ratio > 20.0,
        "flat/pressure-consistent acoustic residual ratio={ratio:.6e} must be > 20; \
         flat={flat:.6e}, consistent={consistent:.6e}, denominator is consistent residual"
    );
}
