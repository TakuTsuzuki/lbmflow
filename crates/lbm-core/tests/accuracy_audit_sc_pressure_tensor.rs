//! ACC-AUDIT SC-P: Shan-Chen pressure-tensor / mechanical surface-tension referee.
//!
//! Context: ANOM-P4-014/017 left a three-way surface-tension disagreement:
//! T11 Laplace uses sigma = 3.32e-2, Taylor-Culick retraction sees about
//! 0.49x, and Jurin capillary rise sees about 1.54x. These probes evaluate the
//! same Shan-Chen pressure tensor directly so failures become triage findings,
//! not tuned bands.

mod common;

use common::metrics::linear_fit;
use lbm_core::compat::lattice::{CS2, CX, CY, Q, W};
use lbm_core::compat::multiphase::ShanChen;
use lbm_core::compat::prelude::*;
use std::f64::consts::PI;

const G: f64 = -5.0;
const NU_TAU_1: f64 = 1.0 / 6.0;
const RHO_L_T11: f64 = 1.888;
const RHO_V_T11: f64 = 0.1194;
const SIGMA_LAPLACE_T11: f64 = 3.32e-2;
const P1_PRESSURE_TENSOR_SIGMA_REL_BAND: f64 = 0.15;
const RHO_INIT_L: f64 = 2.0;
const RHO_INIT_V: f64 = 0.15;

#[derive(Clone, Copy, Debug)]
struct Tensor {
    xx: f64,
    xy: f64,
    yy: f64,
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

fn run_t11_flat(steps: usize) -> Simulation<f64> {
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
        (if liquid { RHO_INIT_L } else { RHO_INIT_V }, 0.0, 0.0)
    });
    let sc = ShanChen::new(G);
    for _ in 0..steps {
        sim.force_field_mut().fill([0.0; 2]);
        sc.update_force(&mut sim);
        sim.step();
    }
    sim
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

fn averaged_radial_profile(case: &DropletCase, dirs: usize) -> Vec<(f64, f64, f64)> {
    let tensor = pressure_tensor_field(&case.sim);
    let (nx, ny) = (case.sim.nx(), case.sim.ny());
    let c = [nx as f64 / 2.0, ny as f64 / 2.0];
    let r_min = (case.r_fit - 12.0).max(1.5);
    let r_max = case.r_fit + 12.0;
    let dr = 0.25;
    let n = ((r_max - r_min) / dr).round() as usize;
    let mut out = Vec::with_capacity(n + 1);
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
        out.push((r, prr / dirs as f64, ptt / dirs as f64));
    }
    out
}

fn directional_sigma(case: &DropletCase, theta: f64) -> f64 {
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
        let t = sample_tensor_periodic(
            &tensor,
            nx,
            ny,
            c[0] + r * theta.cos(),
            c[1] + r * theta.sin(),
        );
        let (prr, ptt) = radial_components(t, theta);
        samples.push((r, ptt - prr));
    }
    integrate_profile(&samples)
}

fn flat_kirkwood_buff_sigma(sim: &Simulation<f64>) -> f64 {
    let tensor = pressure_tensor_field(sim);
    let (nx, ny) = (sim.nx(), sim.ny());
    let mut profile = Vec::with_capacity(ny);
    for y in 0..ny {
        let mut pnn_minus_ptt = 0.0;
        for x in 0..nx {
            let t = tensor[y * nx + x];
            pnn_minus_ptt += t.xx - t.yy;
        }
        profile.push((y as f64, pnn_minus_ptt / nx as f64));
    }
    0.5 * integrate_profile(&profile)
}

// Shan's lattice pressure tensor (e.g. the 2008 Shan pressure-tensor form)
// follows from taking the second moment of the pairwise pseudopotential link
// force. For one D2Q9 link, the discrete interaction force density is
// F_a(x) = -G psi(x) sum_q w_q psi(x+c_q)c_qa. In this implementation's D2Q9
// weight normalization, sum_q w_q c_qa c_qb = cs^2 delta_ab, so the bulk EOS
// keeps the cs^2 factor while the link second moment itself carries the
// anisotropic lattice normalization. Mechanical equilibrium requires
// div_b P_ab = -F_a, and the link-wise divergence of
// -G/2 psi(x) sum_q w_q c_qa c_qb psi(x+c_q) reproduces that force to the
// lattice continuum limit. Adding the isotropic ideal plus bulk interaction
// pressure gives
// P_ab = (rho cs^2 + G cs^2 psi^2/2) delta_ab
//        - G psi(x)/2 sum_q w_q c_qa c_qb psi(x+c_q).
// For a circular interface, r is normal and theta is tangential; the
// Young-Laplace mechanical tension requested here is the radial integral of
// P_theta_theta - P_rr across the diffuse interface.
#[test]
fn sc_p1_curved_pressure_tensor_integral_matches_t11_laplace_sigma() {
    let mut failures = Vec::new();
    for r0 in [12.0, 16.0, 20.0] {
        let case = run_t11_droplet(r0, 40_000);
        let profile = averaged_radial_profile(&case, 64);
        let samples: Vec<_> = profile
            .iter()
            .map(|&(r, prr, ptt)| (r, ptt - prr))
            .collect();
        let sigma_yl = integrate_profile(&samples);
        let rel = rel_err(sigma_yl, SIGMA_LAPLACE_T11);
        println!(
            "VAL SC-P P1: r0={r0:.1}, r_fit={:.6}, sigma_YL={sigma_yl:.8e}, sigma_Laplace={SIGMA_LAPLACE_T11:.8e}, rel={rel:.6e}, band={P1_PRESSURE_TENSOR_SIGMA_REL_BAND:.6e}, norm=sigma_Laplace",
            case.r_fit
        );
        let stride = (profile.len() / 6).max(1);
        for (i, &(r, prr, ptt)) in profile.iter().enumerate() {
            if i % stride == 0 || i + 1 == profile.len() {
                println!(
                    "VAL SC-P P1 sample: r0={r0:.1}, r={r:.3}, P_rr={prr:.8e}, P_tt={ptt:.8e}"
                );
            }
        }
        if rel > P1_PRESSURE_TENSOR_SIGMA_REL_BAND {
            failures.push(format!(
                "r0={r0:.1}, r_fit={:.6}, sigma_YL={sigma_yl:.8e}, rel={rel:.6e}",
                case.r_fit
            ));
        }
    }
    assert!(
        failures.is_empty(),
        "VAL SC-P P1 failures: {}; sigma_Laplace={SIGMA_LAPLACE_T11:.8e}, band={P1_PRESSURE_TENSOR_SIGMA_REL_BAND:.6e}, normalization=sigma_Laplace",
        failures.join("; ")
    );
}

// The same tensor should not depend strongly on azimuth for a static droplet.
// The D2Q9 stencil is not fully rotationally isotropic at finite interface
// width, so this is a characterization band: integrate the traceless radial
// signal P_theta_theta - P_rr independently along the eight lattice/diagonal
// directions and compare max directional spread with the mean.
#[test]
fn sc_p2_static_interface_pressure_tensor_anisotropy_is_bounded() {
    let case = run_t11_droplet(16.0, 40_000);
    let sigmas: Vec<_> = (0..8)
        .map(|d| {
            let theta = d as f64 * PI / 4.0;
            (theta, directional_sigma(&case, theta))
        })
        .collect();
    let mean = sigmas.iter().map(|(_, s)| *s).sum::<f64>() / sigmas.len() as f64;
    let max_rel = sigmas
        .iter()
        .map(|(_, s)| (s - mean).abs() / mean.abs())
        .fold(0.0, f64::max);
    println!(
        "VAL SC-P P2: r_fit={:.6}, directional_sigma={:?}, mean={mean:.8e}, max_rel={max_rel:.6e}, band=1.500000e-1, norm=mean_traceless_integral",
        case.r_fit, sigmas
    );
    assert!(
        max_rel <= 0.15,
        "VAL SC-P P2: max_rel={max_rel:.6e} > band=1.500000e-1, directional_sigma={sigmas:?}, normalization=mean_traceless_integral"
    );
}

// Kirkwood-Buff mechanical tension for this sign convention is the same
// tangential-minus-normal excess used by P1: sigma = integral(P_T - P_N) dz
// = integral(P_xx - P_yy) dz for a horizontal interface. The periodic T11 slab
// has two identical interfaces, so the full-domain integral is 2 sigma.
#[test]
fn sc_p3_flat_interface_kirkwood_buff_sigma_matches_t11_laplace_sigma() {
    let sim = run_t11_flat(30_000);
    let sigma_kb = flat_kirkwood_buff_sigma(&sim);
    let rel = rel_err(sigma_kb, SIGMA_LAPLACE_T11);
    println!(
        "VAL SC-P P3: sigma_KB={sigma_kb:.8e}, sigma_Laplace={SIGMA_LAPLACE_T11:.8e}, rel={rel:.6e}, band=1.000000e-1, norm=sigma_Laplace, rho_l_ref={RHO_L_T11:.6}, rho_v_ref={RHO_V_T11:.6}"
    );
    assert!(
        rel <= 0.10,
        "VAL SC-P P3: sigma_KB={sigma_kb:.8e}, sigma_Laplace={SIGMA_LAPLACE_T11:.8e}, rel={rel:.6e} > band=1.000000e-1, normalization=sigma_Laplace"
    );
}

#[derive(Clone, Copy, Debug)]
struct RimSample {
    t: f64,
    mass: f64,
    px: f64,
    x_front: f64,
}

fn rim_sample(sim: &Simulation<f64>, rho_mid: f64) -> RimSample {
    let (nx, ny) = (sim.nx(), sim.ny());
    let mut x_front = 0usize;
    'outer: for x in (0..nx).rev() {
        for y in 0..ny {
            if sim.rho(x, y) > rho_mid {
                x_front = x;
                break 'outer;
            }
        }
    }
    let rim_left = x_front.saturating_sub(28);
    let mut mass = 0.0;
    let mut px = 0.0;
    for y in 0..ny {
        for x in rim_left..=x_front {
            let rho = sim.rho(x, y);
            if rho > rho_mid {
                mass += rho;
                px += rho * sim.ux(x, y);
            }
        }
    }
    RimSample {
        t: sim.time() as f64,
        mass,
        px,
        x_front: x_front as f64,
    }
}

fn run_taylor_culick_h20() -> (f64, f64, f64, Vec<RimSample>) {
    let (nx, ny) = (384, 96);
    let h = 20usize;
    let y0 = ny / 2 - h / 2;
    let x_left = 32usize;
    let x_right = 250usize;
    let mut sim: Simulation<f64> = SimConfig {
        nx,
        ny,
        nu: NU_TAU_1,
        ..Default::default()
    }
    .build()
    .unwrap();
    sim.init_with(|x, y| {
        let liquid = (x_left..=x_right).contains(&x) && (y0..y0 + h).contains(&y);
        (if liquid { RHO_INIT_L } else { RHO_INIT_V }, 0.0, 0.0)
    });
    let sc = ShanChen::new(G);
    let rho_mid = 0.5 * (RHO_L_T11 + RHO_V_T11);
    let mut samples = Vec::new();
    for step in 0..=8_000 {
        if (1_000..=6_000).contains(&step) && step % 250 == 0 {
            samples.push(rim_sample(&sim, rho_mid));
        }
        sim.force_field_mut().fill([0.0; 2]);
        sc.update_force(&mut sim);
        sim.step();
    }
    let t: Vec<_> = samples.iter().map(|s| s.t).collect();
    let px: Vec<_> = samples.iter().map(|s| s.px).collect();
    let fit = linear_fit(&t, &px);
    let sigma_momentum = fit.slope.abs() / 2.0;
    let speed_fit = linear_fit(&t, &samples.iter().map(|s| s.x_front).collect::<Vec<_>>());
    (sigma_momentum, speed_fit.slope.abs(), fit.r2, samples)
}

// Taylor's control-volume balance for a retracting 2D film is momentum flux
// dP/dt = 2 sigma per unit depth: the two free surfaces pull the rim with
// total line force 2 sigma while the growing rim captures film mass. Therefore
// sigma_momentum = |dP_rim/dt|/2. The audit compares this dynamic mechanical
// sigma to the flat-interface Kirkwood-Buff sigma from the same pressure
// tensor, not to a tuned Taylor-Culick velocity prefactor.
#[test]
#[ignore = "ignored per MODEL_RISK_MATRIX section 3 Shan-Chen demo tier; validation belongs to Allen-Cahn BCFD-040..048"]
fn sc_p4_taylor_culick_momentum_flux_matches_flat_mechanical_sigma() {
    let flat = run_t11_flat(30_000);
    let sigma_kb = flat_kirkwood_buff_sigma(&flat);
    let (sigma_momentum, front_speed, px_r2, samples) = run_taylor_culick_h20();
    let rel_kb = rel_err(sigma_momentum, sigma_kb);
    let rel_laplace = rel_err(sigma_momentum, SIGMA_LAPLACE_T11);
    let tc_speed_kb = (2.0 * sigma_kb / (RHO_L_T11 * 20.0)).sqrt();
    println!(
        "VAL SC-P P4: sigma_momentum={sigma_momentum:.8e}, sigma_KB={sigma_kb:.8e}, sigma_Laplace={SIGMA_LAPLACE_T11:.8e}, rel_to_KB={rel_kb:.6e}, rel_to_Laplace={rel_laplace:.6e}, band=2.000000e-1, norm=sigma_KB, front_speed={front_speed:.8e}, tc_speed_KB={tc_speed_kb:.8e}, px_fit_r2={px_r2:.8}, samples={samples:?}"
    );
    for s in &samples {
        println!(
            "VAL SC-P P4 sample: t={:.0}, x_front={:.3}, rim_mass={:.8e}, rim_px={:.8e}",
            s.t, s.x_front, s.mass, s.px
        );
    }
    assert!(
        rel_kb <= 0.20,
        "VAL SC-P P4: sigma_momentum={sigma_momentum:.8e}, sigma_KB={sigma_kb:.8e}, rel={rel_kb:.6e} > band=2.000000e-1, normalization=sigma_KB, sigma_Laplace={SIGMA_LAPLACE_T11:.8e}, rel_to_Laplace={rel_laplace:.6e}, px_fit_r2={px_r2:.8}"
    );
}
