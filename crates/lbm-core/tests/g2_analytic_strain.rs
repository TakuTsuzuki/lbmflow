//! G2 analytic-strain gate: resolved dissipation from native shear-rate output.

mod common;

use common::run_to_steady;
use lbm_core::compat::prelude::*;

#[derive(Debug)]
struct EpsilonMetrics {
    linf_rel: f64,
    mean_rel: f64,
    mean_measured: f64,
    mean_reference: f64,
    max_abs_error: f64,
}

fn channel_sim(ny: usize, nu: f64, force: [f64; 2], top: EdgeBC<f64>) -> Simulation<f64> {
    SimConfig {
        nx: 4,
        ny,
        nu,
        collision: Collision::default(),
        edges: Edges {
            left: EdgeBC::Periodic,
            right: EdgeBC::Periodic,
            bottom: EdgeBC::BounceBack,
            top,
        },
        force,
        ..Default::default()
    }
    .build()
    .unwrap()
}

fn interior_epsilon(sim: &Simulation<f64>, y_range: std::ops::Range<usize>) -> Vec<(usize, f64)> {
    let shear = sim.shear_rate_field();
    let mut out = Vec::new();
    for y in y_range {
        for x in 0..sim.nx() {
            let eps = sim.nu() * shear[y * sim.nx() + x].powi(2);
            assert!(
                eps.is_finite(),
                "G2 epsilon must be finite at solid-adjacent/interior cell ({x},{y}): {eps:e}"
            );
            out.push((y, eps));
        }
    }
    out
}

fn compare_epsilon(
    sim: &Simulation<f64>,
    y_range: std::ops::Range<usize>,
    reference_at_y: impl Fn(usize) -> f64,
    mean_reference: f64,
) -> EpsilonMetrics {
    let values = interior_epsilon(sim, y_range);
    let max_ref = values
        .iter()
        .map(|(y, _)| reference_at_y(*y).abs())
        .fold(0.0f64, f64::max);
    let mut max_abs_error = 0.0f64;
    let mut sum = 0.0f64;
    for (y, eps) in &values {
        let reference = reference_at_y(*y);
        max_abs_error = max_abs_error.max((eps - reference).abs());
        sum += eps;
    }
    let mean_measured = sum / values.len() as f64;
    EpsilonMetrics {
        linf_rel: max_abs_error / max_ref,
        mean_rel: ((mean_measured - mean_reference) / mean_reference).abs(),
        mean_measured,
        mean_reference,
        max_abs_error,
    }
}

#[test]
fn g2_couette_dissipation_is_uniform_and_matches_analytic_gradient() {
    let u = 0.1;
    let ny = 10;
    let nu = (1.0 - 0.5) / 3.0;
    let h = (ny - 2) as f64;
    let mut sim = channel_sim(ny, nu, [0.0, 0.0], EdgeBC::MovingWall { u: [u, 0.0] });
    assert!(
        run_to_steady(&mut sim, 500, 1.0e-11, 200_000),
        "G2 Couette fixture did not reach steady state, time={}",
        sim.time()
    );

    let expected = nu * (u / h).powi(2);
    let metrics = compare_epsilon(&sim, 2..sim.ny() - 2, |_| expected, expected);
    eprintln!(
        "G2 Couette epsilon: L_inf_rel={:.3e}, mean_rel={:.3e}, max_abs={:.3e}, mean={:.12e}, ref={:.12e}",
        metrics.linf_rel,
        metrics.mean_rel,
        metrics.max_abs_error,
        metrics.mean_measured,
        metrics.mean_reference
    );
    assert!(
        metrics.linf_rel <= 2.0e-4,
        "G2 Couette epsilon L_inf_rel={:.3e} > 2.0e-4, mean_rel={:.3e}, max_abs={:.3e}",
        metrics.linf_rel,
        metrics.mean_rel,
        metrics.max_abs_error
    );
    assert!(
        metrics.mean_rel <= 1.0e-4,
        "G2 Couette epsilon mean_rel={:.3e} > 1.0e-4, L_inf_rel={:.3e}, measured={:.12e}, ref={:.12e}",
        metrics.mean_rel,
        metrics.linf_rel,
        metrics.mean_measured,
        metrics.mean_reference
    );
}

#[test]
fn g2_forced_poiseuille_dissipation_profile_and_volume_mean_match_analytic() {
    let ny = 34;
    let nu = 0.1;
    let g = 1.0e-6;
    let h = (ny - 2) as f64;
    let mut sim = channel_sim(ny, nu, [g, 0.0], EdgeBC::BounceBack);
    assert!(
        run_to_steady(&mut sim, 500, 1.0e-11, 200_000),
        "G2 Poiseuille fixture did not reach steady state, time={}",
        sim.time()
    );

    let eps_at_y = |y: usize| {
        let y_w = y as f64 - 0.5;
        nu * (g / (2.0 * nu) * (h - 2.0 * y_w)).powi(2)
    };
    let continuous_mean = g * g * h * h / (12.0 * nu);
    let metrics = compare_epsilon(&sim, 1..sim.ny() - 1, eps_at_y, continuous_mean);
    eprintln!(
        "G2 Poiseuille epsilon: L_inf_rel={:.3e}, mean_rel={:.3e}, max_abs={:.3e}, mean={:.12e}, integral_ref={:.12e}",
        metrics.linf_rel,
        metrics.mean_rel,
        metrics.max_abs_error,
        metrics.mean_measured,
        metrics.mean_reference
    );
    assert!(
        metrics.linf_rel <= 1.0e-6,
        "G2 Poiseuille epsilon profile L_inf_rel={:.3e} > 1.0e-6, mean_rel={:.3e}, max_abs={:.3e}",
        metrics.linf_rel,
        metrics.mean_rel,
        metrics.max_abs_error
    );
    assert!(
        metrics.mean_rel <= 1.1e-3,
        "G2 Poiseuille epsilon volume mean_rel={:.3e} > 1.1e-3, L_inf_rel={:.3e}, measured={:.12e}, integral_ref={:.12e}",
        metrics.mean_rel,
        metrics.linf_rel,
        metrics.mean_measured,
        metrics.mean_reference
    );
}

#[test]
fn g2_solid_adjacent_cells_report_finite_shear() {
    let mut sim = channel_sim(10, 0.1, [1.0e-6, 0.0], EdgeBC::BounceBack);
    assert!(
        run_to_steady(&mut sim, 500, 1.0e-11, 200_000),
        "G2 finite-adjacent fixture did not reach steady state, time={}",
        sim.time()
    );
    let shear = sim.shear_rate_field();
    for y in [0, 1, sim.ny() - 2, sim.ny() - 1] {
        for x in 0..sim.nx() {
            let gamma = shear[y * sim.nx() + x];
            let eps = sim.nu() * gamma * gamma;
            assert!(
                gamma.is_finite() && eps.is_finite(),
                "G2 shear/epsilon must be finite at x={x}, y={y}: gamma={gamma:e}, eps={eps:e}"
            );
        }
    }
}
