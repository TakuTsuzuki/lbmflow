// Inherited verbatim from the retired V1 suite at its retirement (2026-07-05,
// scripts/sync-tests.sh mechanical retarget); now the canonical facade tests.
//! Validation T11: Shan-Chen single-component multiphase behaviour.

use lbm_core::compat::multiphase::ShanChen;
use lbm_core::compat::prelude::*;
use std::f64::consts::PI;

const G: f64 = -5.0;
const NU_TAU_1: f64 = 1.0 / 6.0;
const RHO_L_REF: f64 = 1.888;
const RHO_V_REF: f64 = 0.1194;
const SIGMA_REF: f64 = 3.32e-2;

#[derive(Clone, Copy, Debug)]
struct FlatStats {
    rho_l: f64,
    rho_v: f64,
    p_rel: f64,
    umax: f64,
    mass_drift: f64,
}

#[derive(Clone, Copy, Debug)]
struct DropStats {
    r0: f64,
    r_fit: f64,
    dp: f64,
    sigma_local: f64,
}

fn rel_err(actual: f64, expected: f64) -> f64 {
    ((actual - expected) / expected).abs()
}

fn run_flat_f64(steps: usize) -> FlatStats {
    let (nx, ny) = (64, 128);
    let mut sim: Simulation<f64> = SimConfig {
        nx,
        ny,
        nu: NU_TAU_1,
        ..Default::default()
    }
    .build()
    .unwrap();
    sim.init_with(|_, y| {
        let liquid = y >= ny / 4 && y < 3 * ny / 4;
        (if liquid { 2.0 } else { 0.15 }, 0.0, 0.0)
    });
    let m0 = sim.total_mass();
    let sc = ShanChen::new(G);
    for _ in 0..steps {
        sim.force_field_mut().fill([0.0; 2]);
        sc.update_force(&mut sim);
        sim.step();
    }

    let rho_l = sim.rho(nx / 2, ny / 2);
    let rho_v = sim.rho(nx / 2, 0);
    let p_l = sc.pressure(rho_l);
    let p_v = sc.pressure(rho_v);
    let p_rel = ((p_l - p_v) / p_l).abs();
    let umax = sim
        .ux_field()
        .iter()
        .chain(sim.uy_field())
        .fold(0.0f64, |a, v| a.max(v.abs()));
    let mass_drift = ((sim.total_mass() - m0) / m0).abs();

    FlatStats {
        rho_l,
        rho_v,
        p_rel,
        umax,
        mass_drift,
    }
}

fn run_flat_f32(steps: usize) -> (f64, f64, bool) {
    let (nx, ny) = (64, 128);
    let mut sim: Simulation<f32> = SimConfig {
        nx,
        ny,
        nu: NU_TAU_1,
        ..Default::default()
    }
    .build()
    .unwrap();
    sim.init_with(|_, y| {
        let liquid = y >= ny / 4 && y < 3 * ny / 4;
        (if liquid { 2.0 } else { 0.15 }, 0.0, 0.0)
    });
    let sc = ShanChen::new(G);
    for _ in 0..steps {
        sim.force_field_mut().fill([0.0; 2]);
        sc.update_force(&mut sim);
        sim.step();
    }
    let finite = sim.rho_field().iter().all(|v| v.is_finite())
        && sim.ux_field().iter().all(|v| v.is_finite())
        && sim.uy_field().iter().all(|v| v.is_finite());
    (
        sim.rho(nx / 2, ny / 2) as f64,
        sim.rho(nx / 2, 0) as f64,
        finite,
    )
}

fn run_droplet(r0: f64, steps: usize) -> DropStats {
    let n = 128;
    let mut sim: Simulation<f64> = SimConfig {
        nx: n,
        ny: n,
        nu: NU_TAU_1,
        ..Default::default()
    }
    .build()
    .unwrap();
    let c = n as f64 / 2.0;
    sim.init_with(|x, y| {
        let d = ((x as f64 - c).powi(2) + (y as f64 - c).powi(2)).sqrt();
        (if d < r0 { 2.0 } else { 0.15 }, 0.0, 0.0)
    });
    let sc = ShanChen::new(G);
    for _ in 0..steps {
        sim.force_field_mut().fill([0.0; 2]);
        sc.update_force(&mut sim);
        sim.step();
    }

    let rho_in = sim.rho(n / 2, n / 2);
    let rho_out = sim.rho(2, 2);
    let dp = sc.pressure(rho_in) - sc.pressure(rho_out);
    let rho_mid = 0.5 * (rho_in + rho_out);
    let area = sim.rho_field().iter().filter(|&&r| r > rho_mid).count() as f64;
    let r_fit = (area / PI).sqrt();
    DropStats {
        r0,
        r_fit,
        dp,
        sigma_local: dp * r_fit,
    }
}

#[test]
fn t11_shan_chen_adds_to_existing_force_field_anom_p4_022() {
    let (nx, ny) = (32, 32);
    let build = || {
        let mut sim: Simulation<f64> = SimConfig {
            nx,
            ny,
            nu: NU_TAU_1,
            ..Default::default()
        }
        .build()
        .unwrap();
        sim.init_with(|x, y| {
            let kx = 2.0 * PI * x as f64 / nx as f64;
            let ky = 2.0 * PI * y as f64 / ny as f64;
            (1.0 + 0.10 * kx.cos() + 0.05 * ky.sin(), 0.0, 0.0)
        });
        sim
    };
    let gravity = [1.0e-6, -2.0e-6];
    let gravity_force = |sim: &Simulation<f64>| {
        sim.rho_field()
            .iter()
            .map(|&rho| [rho * gravity[0], rho * gravity[1]])
            .collect::<Vec<_>>()
    };

    let sc = ShanChen::new(G);
    let mut sc_only = build();
    sc_only.force_field_mut().fill([0.0; 2]);
    sc.update_force(&mut sc_only);
    let sc_force = sc_only.force_field_mut().to_vec();

    let mut gravity_only = build();
    let g_force = gravity_force(&gravity_only);
    gravity_only.force_field_mut().copy_from_slice(&g_force);

    let mut composed = build();
    composed.force_field_mut().copy_from_slice(&g_force);
    sc.update_force(&mut composed);
    let composed_force = composed.force_field_mut().to_vec();

    let mut max_sum_err = 0.0f64;
    let mut max_sc = 0.0f64;
    let mut max_gravity = 0.0f64;
    let mut max_delta_from_sc_only = 0.0f64;
    let mut max_delta_from_gravity_only = 0.0f64;
    for ((got, sc), g) in composed_force.iter().zip(&sc_force).zip(&g_force) {
        for c in 0..2 {
            max_sum_err = max_sum_err.max((got[c] - sc[c] - g[c]).abs());
            max_sc = max_sc.max(sc[c].abs());
            max_gravity = max_gravity.max(g[c].abs());
            max_delta_from_sc_only = max_delta_from_sc_only.max((got[c] - sc[c]).abs());
            max_delta_from_gravity_only = max_delta_from_gravity_only.max((got[c] - g[c]).abs());
        }
    }
    assert!(max_sc > 1.0e-6, "SC contribution was too small: {max_sc:e}");
    assert!(
        max_gravity > 1.0e-6,
        "gravity contribution was too small: {max_gravity:e}"
    );
    assert!(
        max_sum_err < 1.0e-14,
        "SC did not add into existing force field: max_sum_err={max_sum_err:e}"
    );
    assert!(
        max_delta_from_sc_only > 0.5 * max_gravity,
        "composed field collapsed to SC alone: delta={max_delta_from_sc_only:e}"
    );
    assert!(
        max_delta_from_gravity_only > 0.5 * max_sc,
        "composed field collapsed to gravity alone: delta={max_delta_from_gravity_only:e}"
    );
}

fn fit_slope_r2(drops: &[DropStats]) -> (f64, f64) {
    let n = drops.len() as f64;
    let xbar = drops.iter().map(|d| 1.0 / d.r_fit).sum::<f64>() / n;
    let ybar = drops.iter().map(|d| d.dp).sum::<f64>() / n;
    let mut sxx = 0.0;
    let mut sxy = 0.0;
    let mut syy = 0.0;
    for d in drops {
        let x = 1.0 / d.r_fit;
        sxx += (x - xbar) * (x - xbar);
        sxy += (x - xbar) * (d.dp - ybar);
        syy += (d.dp - ybar) * (d.dp - ybar);
    }
    (sxy / sxx, (sxy * sxy) / (sxx * syy))
}

#[test]
fn t11_flat_interface_coexistence_pressure_currents_and_mass() {
    let s = run_flat_f64(30_000);
    assert!(
        rel_err(s.rho_l, RHO_L_REF) <= 0.02,
        "T11 flat rho_l = {:.8}, ref = {RHO_L_REF:.8}, rel = {:e}",
        s.rho_l,
        rel_err(s.rho_l, RHO_L_REF)
    );
    assert!(
        rel_err(s.rho_v, RHO_V_REF) <= 0.03,
        "T11 flat rho_v = {:.8}, ref = {RHO_V_REF:.8}, rel = {:e}",
        s.rho_v,
        rel_err(s.rho_v, RHO_V_REF)
    );
    assert!(
        s.p_rel <= 1.0e-4,
        "T11 flat SC-EOS pressure rel = {:e}",
        s.p_rel
    );
    assert!(s.umax <= 5.0e-3, "T11 flat spurious max|u| = {:e}", s.umax);
    assert!(
        s.mass_drift <= 1.0e-10,
        "T11 flat relative mass drift = {:e}",
        s.mass_drift
    );
}

#[test]
fn t11_laplace_single_radius_smoke() {
    let d = run_droplet(16.0, 40_000);
    assert!(
        rel_err(d.sigma_local, SIGMA_REF) <= 0.15,
        "T11 Laplace smoke r0 = {}, r_fit = {:.6}, dp = {:.8e}, sigma = {:.8e}, ref = {:.8e}, rel = {:e}",
        d.r0,
        d.r_fit,
        d.dp,
        d.sigma_local,
        SIGMA_REF,
        rel_err(d.sigma_local, SIGMA_REF)
    );
}

#[test]
#[ignore = "full four-radius 128^2 x 40k Laplace sweep is intentionally outside the default runtime budget"]
fn t11_laplace_four_radius_sweep_is_linear() {
    let drops: Vec<_> = [12.0, 16.0, 20.0, 24.0]
        .into_iter()
        .map(|r0| run_droplet(r0, 40_000))
        .collect();
    let (slope, r2) = fit_slope_r2(&drops);
    assert!(
        r2 >= 0.999,
        "T11 Laplace R^2 = {:.8}, slope = {:.8e}, drops = {:?}",
        r2,
        slope,
        drops
    );
    assert!(
        rel_err(slope, SIGMA_REF) <= 0.10,
        "T11 Laplace slope sigma = {:.8e}, ref = {:.8e}, rel = {:e}, drops = {:?}",
        slope,
        SIGMA_REF,
        rel_err(slope, SIGMA_REF),
        drops
    );
    for d in drops {
        assert!(
            rel_err(d.sigma_local, slope) <= 0.05,
            "T11 Laplace local sigma r0 = {}, r_fit = {:.6}, dp = {:.8e}, sigma = {:.8e}, slope = {:.8e}, rel = {:e}",
            d.r0,
            d.r_fit,
            d.dp,
            d.sigma_local,
            slope,
            rel_err(d.sigma_local, slope)
        );
    }
}

#[test]
fn t11_f32_flat_interface_stays_finite_and_close_to_f64() {
    let (rho_l, rho_v, finite) = run_flat_f32(30_000);
    assert!(finite, "T11 f32 flat interface produced NaN/Inf");
    assert!(
        rel_err(rho_l, RHO_L_REF) <= 0.05,
        "T11 f32 flat rho_l = {:.8}, f64 ref = {RHO_L_REF:.8}, rel = {:e}",
        rho_l,
        rel_err(rho_l, RHO_L_REF)
    );
    assert!(
        rel_err(rho_v, RHO_V_REF) <= 0.05,
        "T11 f32 flat rho_v = {:.8}, f64 ref = {RHO_V_REF:.8}, rel = {:e}",
        rho_v,
        rel_err(rho_v, RHO_V_REF)
    );
}
