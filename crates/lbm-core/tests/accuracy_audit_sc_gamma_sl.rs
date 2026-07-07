//! ACC-AUDIT SC-GSL: direct solid-liquid interfacial-tension measurement.
//!
//! Radar residual for ANOM-P4-023 / ANOM-P4-014: the flat liquid-vapor
//! Kirkwood-Buff pressure-tensor tension is about gamma_lv = 3.6e-2, while
//! Jurin capillary rise in a wetting slot inferred a larger effective tension.
//! This test measures the wall-side solid-liquid contribution directly, using
//! the same lattice pressure-tensor construction as the SC pressure-tensor
//! audits, but with solid neighbors represented by the virtual wall density
//! used by `ShanChen::with_wall_rho`.

use lbm_core::compat::lattice::{CS2, CX, CY, Q, W};
use lbm_core::compat::multiphase::ShanChen;
use lbm_core::compat::prelude::*;
use std::f64::consts::PI;
use std::fs::{create_dir_all, File};
use std::io::{Result as IoResult, Write};
use std::path::PathBuf;

const G: f64 = -5.0;
const NU_TAU_1: f64 = 1.0 / 6.0;
const NX: usize = 128;
const NY: usize = 128;
const SLAB_THICKNESS: usize = 16;
const STEPS: usize = 30_000;

const RHO_INIT_L: f64 = 2.0;
const RHO_L_T11: f64 = 1.888;
const WALL_RHO_WET: f64 = 1.0;
const WALL_RHO_NEUTRAL: f64 = RHO_L_T11;
const THETA_T11C_DEG: f64 = 63.0;

const GAMMA_LV_KB_FLAT_P4_023: f64 = 3.6e-2;
const SIGMA_LAPLACE_T11: f64 = 3.32e-2;
const JURIN_SIGMA_EFF_P4_014: f64 = 1.54 * SIGMA_LAPLACE_T11;

#[derive(Clone, Copy, Debug)]
struct Tensor {
    xx: f64,
    yy: f64,
}

#[derive(Debug)]
struct GammaSlStats {
    label: &'static str,
    wall_rho: f64,
    gamma_sl: f64,
    bulk_anisotropy: f64,
    rho_first: f64,
    rho_bulk: f64,
    min_rho: f64,
    max_rho: f64,
    max_u: f64,
    profile: Vec<(f64, f64, f64)>,
    artifact_path: String,
}

fn psi(rho: f64) -> f64 {
    1.0 - (-rho).exp()
}

fn rel_err(actual: f64, expected: f64) -> f64 {
    (actual - expected).abs() / expected.abs().max(1.0e-30)
}

fn wrapped_neighbor(sim: &Simulation<f64>, x: usize, y: usize, dx: i32, dy: i32) -> Option<usize> {
    let nx = sim.nx();
    let ny = sim.ny();
    let mut xp = x as isize + dx as isize;
    let mut yp = y as isize + dy as isize;
    if xp < 0 || xp >= nx as isize {
        if sim.is_periodic_x() {
            xp = xp.rem_euclid(nx as isize);
        } else {
            return None;
        }
    }
    if yp < 0 || yp >= ny as isize {
        if sim.is_periodic_y() {
            yp = yp.rem_euclid(ny as isize);
        } else {
            return None;
        }
    }
    Some(yp as usize * nx + xp as usize)
}

fn pressure_tensor_field_with_virtual_wall(sim: &Simulation<f64>, wall_rho: f64) -> Vec<Tensor> {
    let (nx, ny) = (sim.nx(), sim.ny());
    let rho = sim.rho_field();
    let solid = sim.solid_field();
    let psi_wall = psi(wall_rho);
    let psi_field: Vec<f64> = rho
        .iter()
        .zip(solid)
        .map(|(&r, &s)| if s { 0.0 } else { psi(r) })
        .collect();
    let mut out = vec![Tensor { xx: 0.0, yy: 0.0 }; nx * ny];

    for y in 0..ny {
        for x in 0..nx {
            let i = y * nx + x;
            if solid[i] {
                continue;
            }
            let psi_i = psi_field[i];
            let eos = CS2 * rho[i] + 0.5 * G * CS2 * psi_i * psi_i;
            let mut sxx = 0.0;
            let mut syy = 0.0;
            for q in 1..Q {
                let Some(j) = wrapped_neighbor(sim, x, y, CX[q], CY[q]) else {
                    continue;
                };
                let psi_j = if solid[j] { psi_wall } else { psi_field[j] };
                let cx = CX[q] as f64;
                let cy = CY[q] as f64;
                sxx += W[q] * cx * cx * psi_j;
                syy += W[q] * cy * cy * psi_j;
            }
            let pref = -0.5 * G * psi_i;
            out[i] = Tensor {
                xx: eos + pref * sxx,
                yy: eos + pref * syy,
            };
        }
    }
    out
}

fn integrate_profile(samples: &[(f64, f64)]) -> f64 {
    samples
        .windows(2)
        .map(|w| 0.5 * (w[0].1 + w[1].1) * (w[1].0 - w[0].0))
        .sum()
}

fn run_slab_case(wall_rho: f64) -> Simulation<f64> {
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
    sim.set_solid_region(|_, y| y < SLAB_THICKNESS);
    sim.init_with(|_, y| {
        let rho = if y < SLAB_THICKNESS {
            wall_rho
        } else {
            RHO_INIT_L
        };
        (rho, 0.0, 0.0)
    });
    let sc = ShanChen::new(G).with_wall_rho(wall_rho);
    for _ in 0..STEPS {
        sim.force_field_mut().fill([0.0; 2]);
        sc.update_force(&mut sim);
        sim.step();
    }
    sim
}

fn row_mean_density(sim: &Simulation<f64>, y: usize) -> f64 {
    (0..sim.nx()).map(|x| sim.rho(x, y)).sum::<f64>() / sim.nx() as f64
}

fn row_mean_anisotropy(sim: &Simulation<f64>, tensor: &[Tensor], y: usize) -> f64 {
    let nx = sim.nx();
    let mut sum = 0.0;
    let mut count = 0usize;
    for x in 0..nx {
        let i = y * nx + x;
        if !sim.solid_field()[i] {
            let t = tensor[i];
            sum += t.yy - t.xx;
            count += 1;
        }
    }
    sum / count as f64
}

fn dump_density_pgm(sim: &Simulation<f64>, label: &str) -> IoResult<String> {
    let dir = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("target/accuracy_audit");
    create_dir_all(&dir)?;
    let path = dir.join(format!("sc_gamma_sl_{label}.pgm"));
    let mut file = File::create(&path)?;
    writeln!(file, "P2")?;
    writeln!(file, "{} {}", sim.nx(), sim.ny())?;
    writeln!(file, "255")?;
    let (min_rho, max_rho) = sim
        .rho_field()
        .iter()
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), &r| {
            (lo.min(r), hi.max(r))
        });
    for y in (0..sim.ny()).rev() {
        for x in 0..sim.nx() {
            let i = y * sim.nx() + x;
            let px = if sim.solid_field()[i] {
                0
            } else {
                let scaled =
                    ((sim.rho(x, y) - min_rho) / (max_rho - min_rho).max(1.0e-30)).clamp(0.0, 1.0);
                (32.0 + 223.0 * scaled).round() as u8
            };
            write!(file, "{px} ")?;
        }
        writeln!(file)?;
    }
    Ok(path.display().to_string())
}

fn measure_gamma_sl(label: &'static str, wall_rho: f64) -> GammaSlStats {
    let sim = run_slab_case(wall_rho);
    let tensor = pressure_tensor_field_with_virtual_wall(&sim, wall_rho);

    let bulk_start = SLAB_THICKNESS + 36;
    let bulk_end = SLAB_THICKNESS + 48;
    let bulk_anisotropy = (bulk_start..bulk_end)
        .map(|y| row_mean_anisotropy(&sim, &tensor, y))
        .sum::<f64>()
        / (bulk_end - bulk_start) as f64;

    let band_end = SLAB_THICKNESS + 24;
    let profile: Vec<_> = (SLAB_THICKNESS..=band_end)
        .map(|y| {
            let raw = row_mean_anisotropy(&sim, &tensor, y);
            (y as f64, raw, raw - bulk_anisotropy)
        })
        .collect();
    let integral_samples: Vec<_> = profile.iter().map(|&(y, _, excess)| (y, excess)).collect();
    let gamma_sl = integrate_profile(&integral_samples);

    let mut min_rho = f64::INFINITY;
    let mut max_rho = f64::NEG_INFINITY;
    let mut max_u = 0.0f64;
    for y in 0..sim.ny() {
        for x in 0..sim.nx() {
            if sim.is_solid(x, y) {
                continue;
            }
            let rho = sim.rho(x, y);
            min_rho = min_rho.min(rho);
            max_rho = max_rho.max(rho);
            max_u = max_u.max((sim.ux(x, y).powi(2) + sim.uy(x, y).powi(2)).sqrt());
        }
    }
    let artifact_path =
        dump_density_pgm(&sim, label).unwrap_or_else(|e| format!("PGM dump failed: {e}"));

    GammaSlStats {
        label,
        wall_rho,
        gamma_sl,
        bulk_anisotropy,
        rho_first: row_mean_density(&sim, SLAB_THICKNESS),
        rho_bulk: row_mean_density(&sim, SLAB_THICKNESS + 48),
        min_rho,
        max_rho,
        max_u,
        profile,
        artifact_path,
    }
}

fn assert_physical_measurement(stats: &GammaSlStats) {
    assert!(
        stats.gamma_sl.is_finite()
            && stats.bulk_anisotropy.is_finite()
            && stats.rho_first.is_finite()
            && stats.rho_bulk.is_finite()
            && stats.min_rho.is_finite()
            && stats.max_rho.is_finite()
            && stats.max_u.is_finite(),
        "SC gamma_sl {} produced a non-finite measurement: {stats:?}",
        stats.label
    );
    assert!(
        stats.min_rho > 0.0 && stats.max_rho < 3.0,
        "SC gamma_sl {} density bounds are nonphysical: min_rho={:.8e}, max_rho={:.8e}, stats={stats:?}",
        stats.label,
        stats.min_rho,
        stats.max_rho
    );
    assert!(
        stats.gamma_sl.abs() < 10.0 * GAMMA_LV_KB_FLAT_P4_023,
        "SC gamma_sl {} magnitude is outside the characterization guard: gamma_sl={:.8e}, guard=10*gamma_lv={:.8e}, stats={stats:?}",
        stats.label,
        stats.gamma_sl,
        10.0 * GAMMA_LV_KB_FLAT_P4_023
    );
}

// Young's law at a solid-liquid-vapor contact line is
//
//     gamma_sv - gamma_sl = gamma_lv cos(theta).
//
// The T11c virtual-wall-density measurement gives theta(wall_rho=1.0) = 63 deg.
// The P4-023 flat Kirkwood-Buff pressure-tensor measurement gives
// gamma_lv = 3.6e-2, so the Young-predicted solid-side tension difference is
//
//     gamma_sv - gamma_sl = 3.6e-2 * cos(63 deg)
//                         = 0.4540 * gamma_lv
//                         = 1.634e-2.
//
// Jurin in a narrow slot uses this same liquid-vapor tension:
//
//     h = 2 gamma_lv cos(theta) / (Delta rho g w)
//
// for two wetting walls. The textbook formula therefore does not contain an
// independent gamma_sl term. This audit checks whether the discrete SC
// solid-fluid interaction creates a wall-side excess of the same scale. Because
// gamma_sv is not measured in this file, the printed Young comparison uses the
// neutral-wall solid-liquid measurement as the local unchanged-baseline proxy:
//
//     gamma_sl,Young = gamma_sl,neutral - gamma_lv cos(theta).
//
// The sign is meaningful in Young's convention; the absolute delta is also
// printed because the KB solid-side integral sign depends on the chosen
// normal direction. This is characterization, not a production validation gate.
#[test]
fn sc_gamma_sl_direct_solid_liquid_kb_measurement_radar_residual() {
    let wet = measure_gamma_sl("wall_rho_1p0", WALL_RHO_WET);
    let neutral = measure_gamma_sl("wall_rho_1p888", WALL_RHO_NEUTRAL);

    assert_physical_measurement(&wet);
    assert_physical_measurement(&neutral);

    let young_delta = GAMMA_LV_KB_FLAT_P4_023 * (THETA_T11C_DEG * PI / 180.0).cos();
    let young_gamma_sl_signed = neutral.gamma_sl - young_delta;
    let measured_delta = wet.gamma_sl - neutral.gamma_sl;
    let measured_delta_abs = measured_delta.abs();
    let young_match_signed_rel = rel_err(wet.gamma_sl, young_gamma_sl_signed);
    let young_match_abs_rel = rel_err(measured_delta_abs, young_delta);
    let gamma_sl_extra = measured_delta_abs - young_delta;
    let measurable_delta = measured_delta_abs > 1.0e-5;

    println!(
        "VAL SC-GSL wet: wall_rho={:.6}, gamma_sl={:.8e}, gamma_sl/gamma_lv={:.8}, bulk_anisotropy={:.8e}, rho_first={:.8}, rho_bulk={:.8}, min_rho={:.8}, max_rho={:.8}, max_u={:.8e}, artifact={}",
        wet.wall_rho,
        wet.gamma_sl,
        wet.gamma_sl / GAMMA_LV_KB_FLAT_P4_023,
        wet.bulk_anisotropy,
        wet.rho_first,
        wet.rho_bulk,
        wet.min_rho,
        wet.max_rho,
        wet.max_u,
        wet.artifact_path
    );
    println!(
        "VAL SC-GSL neutral: wall_rho={:.6}, gamma_sl={:.8e}, gamma_sl/gamma_lv={:.8}, bulk_anisotropy={:.8e}, rho_first={:.8}, rho_bulk={:.8}, min_rho={:.8}, max_rho={:.8}, max_u={:.8e}, artifact={}",
        neutral.wall_rho,
        neutral.gamma_sl,
        neutral.gamma_sl / GAMMA_LV_KB_FLAT_P4_023,
        neutral.bulk_anisotropy,
        neutral.rho_first,
        neutral.rho_bulk,
        neutral.min_rho,
        neutral.max_rho,
        neutral.max_u,
        neutral.artifact_path
    );
    println!(
        "VAL SC-GSL Young/Jurin: gamma_lv_KB={GAMMA_LV_KB_FLAT_P4_023:.8e}, theta_T11c={THETA_T11C_DEG:.3}, cos_theta={:.8}, Young_delta_gamma=gamma_lv*cos(theta)={young_delta:.8e}, Young_delta/gamma_lv={:.8}, gamma_sl_Young_signed_baseline_proxy={young_gamma_sl_signed:.8e}, measured_delta=wet-neutral={measured_delta:.8e}, abs_delta={measured_delta_abs:.8e}, abs_delta/gamma_lv={:.8}, rel_signed_to_Young={young_match_signed_rel:.8e}, rel_abs_delta_to_Young={young_match_abs_rel:.8e}, gamma_sl_extra_abs_minus_Young={gamma_sl_extra:.8e}, Jurin_sigma_eff_P4_014={JURIN_SIGMA_EFF_P4_014:.8e}, Jurin_sigma_eff/gamma_lv_KB={:.8}, measurable_delta_gt_1e-5={measurable_delta}",
        (THETA_T11C_DEG * PI / 180.0).cos(),
        young_delta / GAMMA_LV_KB_FLAT_P4_023,
        measured_delta_abs / GAMMA_LV_KB_FLAT_P4_023,
        JURIN_SIGMA_EFF_P4_014 / GAMMA_LV_KB_FLAT_P4_023
    );

    for stats in [&wet, &neutral] {
        let stride = (stats.profile.len() / 8).max(1);
        for (i, &(y, raw, excess)) in stats.profile.iter().enumerate() {
            if i % stride == 0 || i + 1 == stats.profile.len() {
                println!(
                    "VAL SC-GSL profile {}: y={y:.3}, Pyy_minus_Pxx_raw={raw:.8e}, bulk_subtracted={excess:.8e}",
                    stats.label
                );
            }
        }
    }
}
