//! ACC-AUDIT SC-R: Shan-Chen pressure-tensor sigma(R) referee.
//!
//! ANOM-P4-023 measured the direct pressure-tensor Young-Laplace tension as
//! radius-dependent: sigma_YL was 13.5% low at R=12, 9.1% low at R=16, and
//! 6.1% low at R=20 relative to the T11 Laplace calibration
//! sigma_Laplace = 3.32e-2. ANOM-P4-014 measured Jurin capillary rise in a
//! slot as 1.54 * sigma_Laplace = 5.11e-2 using the differential form
//! h proportional to sigma cos(theta) * (1/w_slot - 1/w_out). A slot meniscus
//! has finite radius r_m = w_slot / (2 cos(theta)), around 12-30 cells in the
//! gap sweep, which overlaps exactly the P4-023 finite-R deviations.
//!
//! This file tests the Tolman-length-like hypothesis
//!
//!     sigma_YL(R) = sigma_inf + C / R.
//!
//! It does not tune a production model. It characterizes which static
//! pressure-tensor sigma the existing Shan-Chen discretization delivers to
//! finite-curvature menisci, then prints the extrapolated sigma at r_m=12 for
//! comparison with Jurin's 1.54x inferred sigma.

mod common;

use common::metrics::linear_fit;
use lbm_core::compat::lattice::{CS2, CX, CY, Q, W};
use lbm_core::compat::multiphase::ShanChen;
use lbm_core::compat::prelude::*;
use std::f64::consts::PI;

const G: f64 = -5.0;
const NU_TAU_1: f64 = 1.0 / 6.0;
const RHO_INIT_L: f64 = 2.0;
const RHO_INIT_V: f64 = 0.15;
const SIGMA_LAPLACE_T11: f64 = 3.32e-2;
const JURIN_SIGMA_EFFECTIVE: f64 = 1.54 * SIGMA_LAPLACE_T11;
const JURIN_GAP_24_MENISCUS_R: f64 = 12.0;

#[derive(Clone, Copy, Debug)]
struct Tensor {
    xx: f64,
    xy: f64,
    yy: f64,
}

#[derive(Clone, Copy, Debug)]
struct SigmaPoint {
    r0: f64,
    r_fit: f64,
    sigma_yl: f64,
    rel_to_laplace: f64,
}

#[derive(Debug)]
struct DropletCase {
    r_fit: f64,
    sim: Simulation<f64>,
}

fn psi(rho: f64) -> f64 {
    1.0 - (-rho).exp()
}

fn rel_err(actual: f64, expected: f64) -> f64 {
    (actual - expected).abs() / expected.abs()
}

fn periodic_index(i: isize, n: usize) -> usize {
    i.rem_euclid(n as isize) as usize
}

fn run_t11_droplet(r0: f64, steps: usize) -> DropletCase {
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
        (if d < r0 { RHO_INIT_L } else { RHO_INIT_V }, 0.0, 0.0)
    });
    let sc = ShanChen::new(G);
    for _ in 0..steps {
        sim.force_field_mut().fill([0.0; 2]);
        sc.update_force(&mut sim);
        sim.step();
    }

    let rho_mid = 0.5 * (sim.rho(n / 2, n / 2) + sim.rho(2, 2));
    let area = sim.rho_field().iter().filter(|&&r| r > rho_mid).count() as f64;
    DropletCase {
        r_fit: (area / PI).sqrt(),
        sim,
    }
}

fn pressure_tensor_field(sim: &Simulation<f64>) -> Vec<Tensor> {
    let (nx, ny) = (sim.nx(), sim.ny());
    let rho = sim.rho_field();
    let psi_field: Vec<f64> = rho.iter().map(|&r| psi(r)).collect();
    let mut out = vec![
        Tensor {
            xx: 0.0,
            xy: 0.0,
            yy: 0.0,
        };
        nx * ny
    ];

    for y in 0..ny {
        for x in 0..nx {
            let i = y * nx + x;
            let psi_i = psi_field[i];
            let eos = CS2 * rho[i] + 0.5 * G * CS2 * psi_i * psi_i;
            let mut sxx = 0.0;
            let mut sxy = 0.0;
            let mut syy = 0.0;
            for q in 1..Q {
                let xp = periodic_index(x as isize + CX[q] as isize, nx);
                let yp = periodic_index(y as isize + CY[q] as isize, ny);
                let psi_j = psi_field[yp * nx + xp];
                let cx = CX[q] as f64;
                let cy = CY[q] as f64;
                sxx += W[q] * cx * cx * psi_j;
                sxy += W[q] * cx * cy * psi_j;
                syy += W[q] * cy * cy * psi_j;
            }
            let pref = -0.5 * G * psi_i;
            out[i] = Tensor {
                xx: eos + pref * sxx,
                xy: pref * sxy,
                yy: eos + pref * syy,
            };
        }
    }
    out
}

fn sample_tensor_periodic(field: &[Tensor], nx: usize, ny: usize, x: f64, y: f64) -> Tensor {
    let x0 = x.floor();
    let y0 = y.floor();
    let tx = x - x0;
    let ty = y - y0;
    let x0i = periodic_index(x0 as isize, nx);
    let y0i = periodic_index(y0 as isize, ny);
    let x1i = periodic_index(x0 as isize + 1, nx);
    let y1i = periodic_index(y0 as isize + 1, ny);
    let t00 = field[y0i * nx + x0i];
    let t10 = field[y0i * nx + x1i];
    let t01 = field[y1i * nx + x0i];
    let t11 = field[y1i * nx + x1i];
    let blend = |a00: f64, a10: f64, a01: f64, a11: f64| {
        (1.0 - tx) * (1.0 - ty) * a00
            + tx * (1.0 - ty) * a10
            + (1.0 - tx) * ty * a01
            + tx * ty * a11
    };
    Tensor {
        xx: blend(t00.xx, t10.xx, t01.xx, t11.xx),
        xy: blend(t00.xy, t10.xy, t01.xy, t11.xy),
        yy: blend(t00.yy, t10.yy, t01.yy, t11.yy),
    }
}

fn radial_components(t: Tensor, theta: f64) -> (f64, f64) {
    let er = [theta.cos(), theta.sin()];
    let et = [-theta.sin(), theta.cos()];
    let prr = er[0] * er[0] * t.xx + 2.0 * er[0] * er[1] * t.xy + er[1] * er[1] * t.yy;
    let ptt = et[0] * et[0] * t.xx + 2.0 * et[0] * et[1] * t.xy + et[1] * et[1] * t.yy;
    (prr, ptt)
}

fn integrate_profile(samples: &[(f64, f64)]) -> f64 {
    samples
        .windows(2)
        .map(|w| 0.5 * (w[0].1 + w[1].1) * (w[1].0 - w[0].0))
        .sum()
}

fn sigma_yl_from_pressure_tensor(case: &DropletCase, dirs: usize) -> f64 {
    let tensor = pressure_tensor_field(&case.sim);
    let (nx, ny) = (case.sim.nx(), case.sim.ny());
    let c = [nx as f64 / 2.0, ny as f64 / 2.0];
    let r_min = (case.r_fit - 12.0).max(1.5);
    let r_max = case.r_fit + 12.0;
    let dr = 0.25;
    let n = ((r_max - r_min) / dr).round() as usize;
    let mut samples = Vec::with_capacity(n + 1);
    for k in 0..=n {
        let r = r_min + k as f64 * dr;
        let mut prr = 0.0;
        let mut ptt = 0.0;
        for d in 0..dirs {
            let theta = 2.0 * PI * d as f64 / dirs as f64;
            let t = sample_tensor_periodic(
                &tensor,
                nx,
                ny,
                c[0] + r * theta.cos(),
                c[1] + r * theta.sin(),
            );
            let (rcomp, tcomp) = radial_components(t, theta);
            prr += rcomp;
            ptt += tcomp;
        }
        samples.push((r, ptt / dirs as f64 - prr / dirs as f64));
    }
    integrate_profile(&samples)
}

fn measure_sigma_points(radii: &[f64], steps: usize) -> Vec<SigmaPoint> {
    radii
        .iter()
        .map(|&r0| {
            let case = run_t11_droplet(r0, steps);
            let sigma_yl = sigma_yl_from_pressure_tensor(&case, 64);
            let rel_to_laplace = rel_err(sigma_yl, SIGMA_LAPLACE_T11);
            let point = SigmaPoint {
                r0,
                r_fit: case.r_fit,
                sigma_yl,
                rel_to_laplace,
            };
            println!(
                "VAL SC-R point: r0={:.1}, r_fit={:.8}, inv_r={:.8e}, sigma_YL={:.8e}, sigma_Laplace={SIGMA_LAPLACE_T11:.8e}, rel_to_Laplace={:.6e}",
                point.r0,
                point.r_fit,
                1.0 / point.r_fit,
                point.sigma_yl,
                point.rel_to_laplace
            );
            point
        })
        .collect()
}

// Derivation of the measured observable:
// The Shan-Chen interaction is a nearest/diagonal-link pair force
// F_a(x) = -G psi(x) sum_q w_q psi(x+c_q) c_qa. The mechanical pressure
// tensor is the link second moment whose lattice divergence balances that
// force, plus the ideal and bulk-interaction isotropic EOS:
//
// P_ab = (rho cs^2 + G cs^2 psi^2 / 2) delta_ab
//        - G psi(x) / 2 sum_q w_q c_qa c_qb psi(x+c_q).
//
// Across a locally circular interface, the normal component is P_rr and the
// tangential component is P_theta_theta. The pressure jump satisfies
// Delta p = sigma / R at leading order, so the mechanical Young-Laplace
// tension measured directly from the tensor is
//
// sigma_YL(R) = integral_{interface} (P_theta_theta - P_rr) dr.
//
// A finite-curvature correction with one length scale has the first-order
// Tolman form sigma_YL(R) = sigma_inf + C/R. This audit fits sigma_YL against
// 1/R_fit and gates the functional form by r2, not by closeness to the T11
// Laplace calibration. The sign/trend is printed so the static meniscus
// referee is visible without turning this characterization into a tuned model.
fn run_sigma_r_sweep(label: &str, radii: &[f64], steps: usize) {
    let points = measure_sigma_points(radii, steps);
    let inv_r: Vec<_> = points.iter().map(|p| 1.0 / p.r_fit).collect();
    let sigma: Vec<_> = points.iter().map(|p| p.sigma_yl).collect();
    let fit = linear_fit(&inv_r, &sigma);
    let sigma_inf = fit.intercept;
    let c = fit.slope;
    let sigma_at_jurin_gap_24 = sigma_inf + c / JURIN_GAP_24_MENISCUS_R;
    let rel_to_jurin = rel_err(sigma_at_jurin_gap_24, JURIN_SIGMA_EFFECTIVE);
    let rel_to_laplace = rel_err(sigma_at_jurin_gap_24, SIGMA_LAPLACE_T11);
    let monotone_non_decreasing = points.windows(2).all(|w| w[1].sigma_yl >= w[0].sigma_yl);

    println!(
        "VAL SC-R {label}: sigma_fit = sigma_inf + C/R, sigma_inf={sigma_inf:.8e}, C={c:.8e}, r2={:.8}, monotone_non_decreasing={monotone_non_decreasing}, points={points:?}",
        fit.r2
    );
    println!(
        "VAL SC-R {label}: sigma_YL(r_m=12 for Jurin gap-24 meniscus)={sigma_at_jurin_gap_24:.8e}, Jurin_1p54x_sigma_Laplace={JURIN_SIGMA_EFFECTIVE:.8e}, rel_to_Jurin={rel_to_jurin:.6e}, rel_to_Laplace={rel_to_laplace:.6e}, sigma_Laplace={SIGMA_LAPLACE_T11:.8e}"
    );

    assert!(
        fit.r2 >= 0.90,
        "VAL SC-R {label}: inverse-R fit r2={:.8} < band=0.90000000, sigma_inf={sigma_inf:.8e}, C={c:.8e}, points={points:?}",
        fit.r2
    );
}

#[test]
fn sc_sigma_r_light_three_point_inverse_r_characterization() {
    run_sigma_r_sweep("LIGHT", &[12.0, 16.0, 20.0], 40_000);
}

#[test]
#[ignore = "full nine-radius SC pressure-tensor R sweep is intentionally outside the default runtime budget"]
fn sc_sigma_r_full_nine_radius_inverse_r_referee() {
    run_sigma_r_sweep(
        "HEAVY",
        &[6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 20.0, 24.0, 32.0],
        40_000,
    );
}
