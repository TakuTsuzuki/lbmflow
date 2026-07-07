//! G2 analytic-strain gate: resolved dissipation from native shear-rate output.

mod common;

use common::run_to_steady;
use lbm_core::compat::prelude::*;
use lbm_core::prelude::{
    build_wall_rims, CollisionKind, CpuScalar, Face, GlobalSpec, LocalPeriodic, Solver, WaleLes,
    WallSpec, D2Q9, WALE_CW,
};

type CoreCouette = Solver<D2Q9, f64, CpuScalar, LocalPeriodic>;

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

fn core_couette(nx: usize, ny: usize, nu: f64, top_u: f64) -> CoreCouette {
    let spec = GlobalSpec {
        dims: [nx, ny, 1],
        nu,
        collision: CollisionKind::Trt {
            magic: CollisionKind::MAGIC_STD,
        },
        periodic: [true, false, false],
        ..Default::default()
    };
    let mut walls = WallSpec::default();
    walls.is_wall[Face::YNeg.index()] = true;
    walls.is_wall[Face::YPos.index()] = true;
    walls.u[Face::YPos.index()] = [top_u, 0.0, 0.0];
    let (solid, wall_u) = build_wall_rims(2, spec.dims, &walls);
    let mut solver = Solver::new(
        &spec,
        &solid,
        &wall_u,
        [1, 1, 1],
        CpuScalar::default(),
        LocalPeriodic,
    );

    let h = (ny - 2) as f64;
    solver.init_with(|_, y, _| {
        let ux = if y == 0 {
            0.0
        } else if y == ny - 1 {
            top_u
        } else {
            top_u * (y as f64 - 0.5) / h
        };
        (1.0, [ux, 0.0, 0.0])
    });
    solver
}

fn max_abs_pair(a: &[f64], b: &[f64]) -> f64 {
    a.iter()
        .zip(b)
        .map(|(x, y)| (x - y).abs())
        .fold(0.0, f64::max)
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
    // The remaining Couette error is the steady-state residual of the native
    // strain reconstruction at the 1e-11 fixture tolerance; keep only modest
    // headroom because measured*20 would loosen the old band.
    assert!(
        metrics.linf_rel <= 1.5e-4,
        "G2 Couette epsilon L_inf_rel={:.3e} > 1.5e-4, mean_rel={:.3e}, max_abs={:.3e}",
        metrics.linf_rel,
        metrics.mean_rel,
        metrics.max_abs_error
    );
    assert!(
        metrics.mean_rel <= 7.0e-5,
        "G2 Couette epsilon mean_rel={:.3e} > 7.0e-5, L_inf_rel={:.3e}, measured={:.12e}, ref={:.12e}",
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
        metrics.linf_rel <= 3.7e-7,
        "G2 Poiseuille epsilon profile L_inf_rel={:.3e} > 3.7e-7, mean_rel={:.3e}, max_abs={:.3e}",
        metrics.linf_rel,
        metrics.mean_rel,
        metrics.max_abs_error
    );
    // The reference is the continuous volume mean while the test samples cell
    // midpoints; the exact midpoint-quadrature offset is 1/H^2 = 9.766e-4.
    assert!(
        metrics.mean_rel <= 1.0e-3,
        "G2 Poiseuille epsilon volume mean_rel={:.3e} > 1.0e-3, L_inf_rel={:.3e}, measured={:.12e}, integral_ref={:.12e}",
        metrics.mean_rel,
        metrics.linf_rel,
        metrics.mean_measured,
        metrics.mean_reference
    );
}

#[test]
fn wale_pure_shear_couette_has_zero_eddy_viscosity_and_identical_velocity() {
    let nx = 32;
    let ny = 34;
    let nu = 0.02;
    let top_u = 0.02;
    let steps = 5_000;
    assert_eq!(
        WALE_CW.to_bits(),
        0.325f64.to_bits(),
        "P21 fixture requires Nicoud-Ducros WALE Cw=0.325"
    );

    let mut off = core_couette(nx, ny, nu, top_u);
    let mut on = core_couette(nx, ny, nu, top_u);
    let mut les = WaleLes::<f64>::new();
    for _ in 0..steps {
        les.update(&mut on);
        on.run(1);
        off.run(1);
    }
    les.update(&mut on);

    let max_nu_t_interior = les
        .nu_t()
        .iter()
        .enumerate()
        .filter_map(|(i, nu_t)| {
            let y = (i / nx) % ny;
            (2..ny - 2).contains(&y).then_some(nu_t.abs())
        })
        .fold(0.0_f64, f64::max);
    let max_du = max_abs_pair(&off.gather_ux(), &on.gather_ux());
    let max_dv = max_abs_pair(&off.gather_uy(), &on.gather_uy());
    let velocity_band = 2.0 * f64::EPSILON * top_u;
    let diagnostics = les.diagnostics();
    eprintln!(
        "P21 WALE pure-shear Couette: Cw={:.3}, steps={}, max|nu_t|_interior={:.12e}, \
         max_raw_nu_t={:.12e}, clipped_cells={}, max|du|={:.12e}, max|dv|={:.12e}, \
         velocity_band={:.12e}",
        WALE_CW,
        steps,
        max_nu_t_interior,
        diagnostics.max_nu_t_before_clipping,
        diagnostics.clipped_cells,
        max_du,
        max_dv,
        velocity_band
    );

    // WALE uses
    //   nu_t = (Cw Delta)^2 (S^d:S^d)^(3/2)
    //          / ((S:S)^(5/2) + (S^d:S^d)^(5/4)),
    // where S is the symmetric velocity-gradient tensor and
    //   S^d_ij = 0.5 * (g_ik g_kj + g_jk g_ki) - delta_ij tr(g^2)/3.
    // For pure Couette shear, u=(a y,0,0), so the only non-zero gradient is
    // g_xy=a. This nilpotent gradient has g^2=0, hence tr(g^2)=0 and every
    // S^d_ij is exactly zero: no eigenvalue collision or limiting argument is
    // involved. Therefore S^d:S^d=0 and WALE's defining property is
    // nu_t=0 in pure shear. This is the reason to prefer WALE over
    // Smagorinsky for wall-bounded shear (Nicoud-Ducros 1999); if eddy
    // viscosity leaks into this field, the WALE closure is broken.
    assert!(
        max_nu_t_interior <= 1.0e-14,
        "WALE pure-shear property failed: max|nu_t|_interior={max_nu_t_interior:e} > 1e-14"
    );
    assert!(
        max_du <= velocity_band && max_dv <= velocity_band,
        "WALE-on pure-shear velocity changed vs WALE-off: max|du|={max_du:e}, \
         max|dv|={max_dv:e}, band={velocity_band:e}"
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
