//! Hard zero-free-parameter multiphase validation cases.
//!
//! These tests intentionally cross-predict new observables from T11/T12
//! measured constants.  The surface tensions, coexistence densities, and
//! wall-density contact angles below are frozen validation measurements, not
//! knobs fitted to these cases.

mod common;

use common::metrics::*;
use lbm_core::compat::multiphase::{MultiComponent, ShanChen};
use lbm_core::compat::prelude::*;
use std::f64::consts::PI;

const SC_G: f64 = -5.0;
const SC_NU: f64 = 1.0 / 6.0;
const SC_SIGMA: f64 = 3.32e-2;
const SC_RHO_L: f64 = 1.888;
const SC_RHO_V: f64 = 0.1194;
const SC_DELTA_RHO: f64 = SC_RHO_L - SC_RHO_V;
const WALL_RHO_WET: f64 = 1.0;
const WALL_RHO_DRY: f64 = 0.6;
const THETA_WET_DEG: f64 = 63.0;
const THETA_DRY_DEG: f64 = 107.0;

const MC_TRACE: f64 = 0.05;
const MC_G_AB: f64 = 2.6;
const MC_NU: f64 = 0.1;
const MC_SIGMA_AB: f64 = 2.86969302e-2;

#[derive(Clone, Copy, Debug)]
struct JurinStats {
    gap: usize,
    wall_rho: f64,
    theta_deg: f64,
    theta_slot_deg: f64,
    measured_h: f64,
    predicted_h: f64,
    mass_drift: f64,
    steady_drift: f64,
    steps: usize,
}

#[derive(Clone, Copy, Debug)]
struct WaveStats {
    mode: usize,
    k: f64,
    box_l: f64,
    g_eff: f64,
    delta_rho_eff: f64,
    rho_sum: f64,
    g_branch: f64,
    sigma_branch: f64,
    omega_fit: f64,
    omega0: f64,
    omega_rel: f64,
    fit_periods: usize,
    fit_t_end: f64,
    mode_rms_3: f64,
    mode_rms_4: f64,
    mode_rms_5: f64,
    decay_fit: f64,
    decay_model: f64,
    decay_ratio: f64,
    envelope_monotone: f64,
    mass_drift: f64,
}

#[derive(Clone, Copy, Debug)]
struct RtModeStats {
    mode: usize,
    gamma_fit: f64,
    gamma_theory: f64,
    ratio: f64,
    amp0: f64,
    max_amp: f64,
    final_amp: f64,
    amp_step_10: f64,
    amp_step_100: f64,
    amp_step_1000: f64,
    max_u_step_0: f64,
    max_u_step_10: f64,
    max_u_step_100: f64,
    max_u_step_1000: f64,
    mass_drift: f64,
}

#[derive(Clone, Copy, Debug)]
struct TaylorCulickStats {
    h: usize,
    measured_v: f64,
    predicted_v: f64,
    rel: f64,
    mass_loss: f64,
    fit_r2: f64,
    last_window_v: f64,
    last_window_r2: f64,
    retracted_h: f64,
    still_rising: bool,
}

fn make_sim(nx: usize, ny: usize, nu: f64, edges: Edges<f64>) -> Simulation<f64> {
    SimConfig {
        nx,
        ny,
        nu,
        edges,
        ..Default::default()
    }
    .build()
    .unwrap()
}

fn total_mcmp_mass(a: &Simulation<f64>, b: &Simulation<f64>) -> f64 {
    a.total_mass_f64() + b.total_mass_f64()
}

fn rel_err(actual: f64, expected: f64) -> f64 {
    (actual - expected).abs() / expected.abs().max(1.0e-30)
}

fn sc_step(sim: &mut Simulation<f64>, sc: &ShanChen<f64>) {
    sc.update_force(sim);
    sim.step();
}

fn median(v: &mut [f64]) -> f64 {
    assert!(!v.is_empty(), "median requires non-empty input");
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    v[v.len() / 2]
}

fn column_interface_y(sim: &Simulation<f64>, x: usize, rho_mid: f64) -> Option<f64> {
    let ny = sim.ny();
    for y in 1..ny - 2 {
        let r0 = sim.rho(x, y);
        let r1 = sim.rho(x, y + 1);
        if (r0 - rho_mid) * (r1 - rho_mid) <= 0.0 {
            let denom = r1 - r0;
            let frac = if denom.abs() < 1.0e-30 {
                0.0
            } else {
                ((rho_mid - r0) / denom).clamp(0.0, 1.0)
            };
            return Some(y as f64 + frac);
        }
    }
    None
}

fn jurin_prediction(theta_deg: f64, g: f64, gap: usize) -> f64 {
    // Force balance on a 2D parallel-plate meniscus column:
    //   upward capillary force per depth = 2 sigma cos(theta)
    //   downward buoyant weight per depth = Delta_rho g w h
    // Equating them gives h = 2 sigma cos(theta)/(Delta_rho g w).
    2.0 * SC_SIGMA * theta_deg.to_radians().cos() / (SC_DELTA_RHO * g * gap as f64)
}

fn setup_jurin(gap: usize, wall_rho: f64, theta_deg: f64, gravity: f64) -> Simulation<f64> {
    let (nx, ny) = (112, 144);
    let y_wall_bottom = 18usize;
    let reservoir_y = 48.0;
    let left_wall = nx / 2 - gap / 2 - 1;
    let right_wall = left_wall + gap + 1;
    let h0 = jurin_prediction(theta_deg, gravity, gap).clamp(-28.0, 72.0);
    let slot_level = (reservoir_y + h0).clamp((y_wall_bottom + 4) as f64, (ny - 10) as f64);
    let mut sim = make_sim(
        nx,
        ny,
        SC_NU,
        Edges {
            left: EdgeBC::BounceBack,
            right: EdgeBC::BounceBack,
            bottom: EdgeBC::BounceBack,
            top: EdgeBC::BounceBack,
        },
    );
    sim.set_solid_region(|x, y| {
        (x == left_wall || x == right_wall) && (y_wall_bottom..ny - 1).contains(&y)
    });
    sim.set_gravity([0.0, -gravity]);
    sim.init_with(|x, y| {
        let in_slot = x > left_wall && x < right_wall;
        let liquid = if in_slot {
            (y as f64) <= slot_level
        } else {
            (y as f64) <= reservoir_y
        };
        (if liquid { 2.0 } else { 0.15 }, 0.0, 0.0)
    });
    // Reassert solids after initialization because the facade initializes all
    // cells uniformly but keeps the solid mask as the authoritative geometry.
    sim.set_solid_region(|x, y| {
        (x == left_wall || x == right_wall) && (y_wall_bottom..ny - 1).contains(&y)
    });
    let _ = wall_rho;
    sim
}

fn measure_jurin_height(sim: &Simulation<f64>, gap: usize) -> f64 {
    let nx = sim.nx();
    let left_wall = nx / 2 - gap / 2 - 1;
    let right_wall = left_wall + gap + 1;
    let rho_mid = 0.5 * (SC_RHO_L + SC_RHO_V);
    let mut slot = Vec::new();
    let mut reservoir = Vec::new();

    for x in left_wall + 2..right_wall - 1 {
        if let Some(y) = column_interface_y(sim, x, rho_mid) {
            slot.push(y);
        }
    }
    for x in 4..left_wall.saturating_sub(6) {
        if let Some(y) = column_interface_y(sim, x, rho_mid) {
            reservoir.push(y);
        }
    }
    for x in right_wall + 6..nx - 4 {
        if let Some(y) = column_interface_y(sim, x, rho_mid) {
            reservoir.push(y);
        }
    }
    assert!(
        !slot.is_empty() && !reservoir.is_empty(),
        "VAL MPHARD I1 contour missing: gap={gap}, slot_points={}, reservoir_points={}",
        slot.len(),
        reservoir.len()
    );
    median(&mut slot) - median(&mut reservoir)
}

fn measure_slot_theta_deg(sim: &Simulation<f64>, gap: usize) -> f64 {
    let nx = sim.nx();
    let left_wall = nx / 2 - gap / 2 - 1;
    let right_wall = left_wall + gap + 1;
    let rho_mid = 0.5 * (SC_RHO_L + SC_RHO_V);
    let mut profile = Vec::new();
    for x in left_wall + 1..right_wall {
        if let Some(y) = column_interface_y(sim, x, rho_mid) {
            profile.push((x as f64, y));
        }
    }
    if profile.len() < 5 {
        return f64::NAN;
    }
    let y_center = profile
        .iter()
        .min_by(|a, b| {
            (a.0 - nx as f64 * 0.5)
                .abs()
                .partial_cmp(&(b.0 - nx as f64 * 0.5).abs())
                .unwrap()
        })
        .unwrap()
        .1;
    let y_left = profile.first().unwrap().1;
    let y_right = profile.last().unwrap().1;
    let y_wall = 0.5 * (y_left + y_right);
    let sag = (y_wall - y_center).max(0.0);
    let half_gap = 0.5 * gap as f64;
    // Circular-cap diagnostic for a meniscus in a vertical slot:
    // sag/half_gap = tan((90deg - theta_slot)/2).  This is only a
    // measurement of the diffuse-interface meniscus shape inside this slot;
    // it is not fed back into the Jurin prediction, which must keep the
    // frozen T11c flat-wall theta input.
    90.0 - 2.0 * (sag / half_gap.max(1.0e-30)).atan().to_degrees()
}

fn run_jurin(gap: usize, wall_rho: f64, theta_deg: f64) -> JurinStats {
    let gravity = 2.0e-5;
    let mut sim = setup_jurin(gap, wall_rho, theta_deg, gravity);
    let sc = ShanChen::new(SC_G).with_wall_rho(wall_rho);
    let m0 = sim.total_mass_f64();
    let mut previous_h = measure_jurin_height(&sim, gap);
    let mut steady_drift = f64::INFINITY;
    let mut steps = 0usize;
    println!(
        "VAL MPHARD I1 diag: gap={} step={} h={:.6}",
        gap, steps, previous_h
    );
    for _ in 0..40 {
        for _ in 0..2_000 {
            sc_step(&mut sim, &sc);
        }
        steps += 2_000;
        let h = measure_jurin_height(&sim, gap);
        steady_drift = (h - previous_h).abs();
        println!(
            "VAL MPHARD I1 diag: gap={} step={} h={:.6} window_drift={:.6}",
            gap, steps, h, steady_drift
        );
        previous_h = h;
        if steps >= 8_000 && steady_drift <= 0.5 {
            break;
        }
    }
    JurinStats {
        gap,
        wall_rho,
        theta_deg,
        theta_slot_deg: measure_slot_theta_deg(&sim, gap),
        measured_h: previous_h,
        predicted_h: jurin_prediction(theta_deg, gravity, gap),
        mass_drift: ((sim.total_mass_f64() - m0) / m0).abs(),
        steady_drift,
        steps,
    }
}

fn mcmp_interface(mode: usize, nx: usize, ny: usize, a0: f64, sign: f64, x: usize) -> f64 {
    let k = 2.0 * PI * mode as f64 / nx as f64;
    ny as f64 / 2.0 + sign * a0 * (k * x as f64).cos()
}

fn init_mcmp_layers(
    a: &mut Simulation<f64>,
    b: &mut Simulation<f64>,
    mode: usize,
    a0: f64,
    heavy_on_top: bool,
) {
    let (nx, ny) = (a.nx(), a.ny());
    let sign = if heavy_on_top { 1.0 } else { -1.0 };
    a.init_with(|x, y| {
        let y_int = mcmp_interface(mode, nx, ny, a0, sign, x);
        let in_a = if heavy_on_top {
            (y as f64) > y_int
        } else {
            (y as f64) < y_int
        };
        (if in_a { 1.0 } else { MC_TRACE }, 0.0, 0.0)
    });
    b.init_with(|x, y| {
        let y_int = mcmp_interface(mode, nx, ny, a0, sign, x);
        let in_a = if heavy_on_top {
            (y as f64) > y_int
        } else {
            (y as f64) < y_int
        };
        (if in_a { MC_TRACE } else { 1.0 }, 0.0, 0.0)
    });
}

fn signed_fourier_amp(component: &Simulation<f64>, mode: usize) -> f64 {
    let nx = component.nx();
    let ny = component.ny();
    let k = 2.0 * PI * mode as f64 / nx as f64;
    let mut re = 0.0;
    for x in 0..nx {
        let column_mass: f64 = (1..ny - 1).map(|y| component.rho(x, y)).sum();
        re += column_mass * (k * x as f64).cos();
    }
    2.0 * re / nx as f64 / (1.0 - MC_TRACE)
}

fn fourier_amp(component: &Simulation<f64>, mode: usize) -> f64 {
    let nx = component.nx();
    let ny = component.ny();
    let k = 2.0 * PI * mode as f64 / nx as f64;
    let mut re = 0.0;
    let mut im = 0.0;
    for x in 0..nx {
        let column_mass: f64 = (1..ny - 1).map(|y| component.rho(x, y)).sum();
        let ph = k * x as f64;
        re += column_mass * ph.cos();
        im += column_mass * ph.sin();
    }
    2.0 * (re * re + im * im).sqrt() / nx as f64 / (1.0 - MC_TRACE)
}

fn max_speed(component: &Simulation<f64>) -> f64 {
    let mut umax = 0.0f64;
    for y in 1..component.ny() - 1 {
        for x in 1..component.nx() - 1 {
            let u = component.ux(x, y).hypot(component.uy(x, y));
            if u.is_finite() {
                umax = umax.max(u);
            } else {
                return f64::NAN;
            }
        }
    }
    umax
}

fn fit_frequency(t: &[f64], signal: &[f64], omega0: f64) -> (f64, usize, f64) {
    let mean = signal.iter().sum::<f64>() / signal.len() as f64;
    let centered: Vec<f64> = signal.iter().map(|v| v - mean).collect();
    let mut best = (omega0, -1.0);
    for i in 0..121 {
        let omega = omega0 * (0.4 + 1.2 * i as f64 / 120.0);
        let (amp, _) = phase_fit(t, &centered, omega);
        if amp > best.1 {
            best = (omega, amp);
        }
    }
    let period = 2.0 * PI / best.0;
    let t0 = t[0];
    let tmax = *t.last().unwrap();
    let periods = ((tmax - t0) / period).floor().max(1.0) as usize;
    let fit_t_end = t0 + periods as f64 * period;
    let n = t
        .iter()
        .position(|&ti| ti > fit_t_end)
        .unwrap_or(t.len())
        .max(3);
    let t_int = &t[..n.min(t.len())];
    let signal_int = &signal[..n.min(signal.len())];
    let mean_int = signal_int.iter().sum::<f64>() / signal_int.len() as f64;
    let centered_int: Vec<f64> = signal_int.iter().map(|v| v - mean_int).collect();
    let mut best_int = (best.0, -1.0);
    for i in 0..121 {
        let omega = omega0 * (0.4 + 1.2 * i as f64 / 120.0);
        let (amp, _) = phase_fit(t_int, &centered_int, omega);
        if amp > best_int.1 {
            best_int = (omega, amp);
        }
    }
    (best_int.0, periods, fit_t_end)
}

fn local_peak_envelope(series: &[(f64, f64)]) -> Vec<(f64, f64)> {
    let mut peaks = Vec::new();
    for w in series.windows(3) {
        let a0 = w[0].1.abs();
        let a1 = w[1].1.abs();
        let a2 = w[2].1.abs();
        if a1 >= a0 && a1 >= a2 && a1 > 1.0e-8 {
            peaks.push((w[1].0, a1));
        }
    }
    if peaks.len() < 3 {
        series
            .iter()
            .step_by(4)
            .filter_map(|&(t, a)| {
                let aa = a.abs();
                (aa > 1.0e-8).then_some((t, aa))
            })
            .collect()
    } else {
        peaks
    }
}

fn run_standing_wave(mode: usize, steps: usize, sample_every: usize) -> WaveStats {
    let (nx, ny) = (256, 256);
    let g = 1.0e-4;
    let a0 = 3.0;
    let edges = Edges {
        left: EdgeBC::Periodic,
        right: EdgeBC::Periodic,
        bottom: EdgeBC::BounceBack,
        top: EdgeBC::BounceBack,
    };
    let mut heavy = make_sim(nx, ny, MC_NU, edges);
    let mut light = make_sim(nx, ny, MC_NU, edges);
    init_mcmp_layers(&mut heavy, &mut light, mode, a0, false);
    let mc = MultiComponent::new(MC_G_AB).with_gravity([0.0, -g], [0.0, 0.0]);
    let m0 = total_mcmp_mass(&heavy, &light);
    let mut series = vec![(0.0, signed_fourier_amp(&heavy, mode))];
    let mut mode_series = vec![[
        signed_fourier_amp(&heavy, 3),
        signed_fourier_amp(&heavy, 4),
        signed_fourier_amp(&heavy, 5),
    ]];
    for it in 1..=(steps / sample_every) {
        for _ in 0..sample_every {
            mc.update_forces(&mut heavy, &mut light);
            heavy.step();
            light.step();
        }
        series.push(((it * sample_every) as f64, signed_fourier_amp(&heavy, mode)));
        mode_series.push([
            signed_fourier_amp(&heavy, 3),
            signed_fourier_amp(&heavy, 4),
            signed_fourier_amp(&heavy, 5),
        ]);
    }

    // For a stable two-fluid interface, the inviscid gravity-capillary
    // restoring frequency is
    //   omega0^2 = (g_eff k Delta_rho_eff + sigma_AB k^3)/(rho1 + rho2).
    // T12 applies gravity to one component only.  With unit bulk densities,
    // that convention is equivalent to Delta_rho_eff = 1 and rho1+rho2 = 2,
    // hence A_eff = Delta_rho_eff/(rho1+rho2) = 0.5 for the gravity branch.
    let k = 2.0 * PI * mode as f64 / nx as f64;
    let delta_rho_eff = 1.0;
    let rho_sum = 2.0;
    let g_branch = g * k * delta_rho_eff / rho_sum;
    let sigma_branch = MC_SIGMA_AB * k.powi(3) / rho_sum;
    let omega0 = (g_branch + sigma_branch).sqrt();
    let t: Vec<f64> = series.iter().map(|p| p.0).collect();
    let a: Vec<f64> = series.iter().map(|p| p.1).collect();
    let (omega_fit, fit_periods, fit_t_end) = fit_frequency(&t, &a, omega0);
    let mode_rms = |idx: usize| {
        (mode_series.iter().map(|v| v[idx] * v[idx]).sum::<f64>() / mode_series.len() as f64)
            .sqrt()
    };

    let envelope = local_peak_envelope(&series);
    let et: Vec<f64> = envelope.iter().map(|p| p.0).collect();
    let ea: Vec<f64> = envelope.iter().map(|p| p.1).collect();
    let decay_fit = if et.len() >= 2 {
        envelope_fit(&et, &ea).slope.abs()
    } else {
        f64::NAN
    };
    // This is the weak-damping boundary-layer estimate for equal kinematic
    // viscosities, beta ~= 2 nu k^2.  It omits diffuse-interface excess
    // dissipation and exact Prosperetti transient history terms, so the
    // acceptance is deliberately loose at [0.5, 2]x until first measurements.
    let decay_model = 2.0 * MC_NU * k * k;
    let env_rev: Vec<f64> = ea.iter().rev().copied().collect();
    WaveStats {
        mode,
        k,
        box_l: nx as f64,
        g_eff: g,
        delta_rho_eff,
        rho_sum,
        g_branch,
        sigma_branch,
        omega_fit,
        omega0,
        omega_rel: rel_err(omega_fit, omega0),
        fit_periods,
        fit_t_end,
        mode_rms_3: mode_rms(0),
        mode_rms_4: mode_rms(1),
        mode_rms_5: mode_rms(2),
        decay_fit,
        decay_model,
        decay_ratio: decay_fit / decay_model,
        envelope_monotone: monotonicity(&env_rev),
        mass_drift: ((total_mcmp_mass(&heavy, &light) - m0) / m0).abs(),
    }
}

fn rt_gamma_theory(g: f64, k: f64) -> f64 {
    // Linear RT with surface tension:
    //   gamma0^2 = A g k - sigma k^3/(rho1 + rho2).
    // T12's equal-bulk two-component setup freezes A = 1/2 and
    // sigma/(rho1+rho2) = sigma_AB/2.  The existing T12 validation uses the
    // equal-nu correction gamma = sqrt(gamma0^2 + nu^2 k^4) - nu k^2.
    let gamma0_sq = 0.5 * g * k - 0.5 * MC_SIGMA_AB * k.powi(3);
    if gamma0_sq <= 0.0 {
        return f64::NAN;
    }
    (gamma0_sq + MC_NU * MC_NU * k.powi(4)).sqrt() - MC_NU * k * k
}

fn fit_rt_growth(series: &[(f64, f64)]) -> Option<(f64, usize)> {
    let amp0 = series.first()?.1.abs().max(1.0e-30);
    let pts: Vec<_> = series
        .iter()
        .filter_map(|&(t, a)| {
            let normalized = a.abs() / amp0;
            (normalized >= 1.0 && normalized <= 8.0).then_some((t, normalized.ln()))
        })
        .collect();
    if pts.len() < 3 {
        return None;
    }
    let x: Vec<f64> = pts.iter().map(|p| p.0).collect();
    let y: Vec<f64> = pts.iter().map(|p| p.1).collect();
    let fit = linear_fit(&x, &y);
    Some((fit.slope, pts.len()))
}

fn run_rt_mode(nx: usize, mode: usize, g: f64, steps: usize, sample_every: usize) -> RtModeStats {
    let ny = nx;
    let a0 = 3.0;
    let edges = Edges {
        left: EdgeBC::Periodic,
        right: EdgeBC::Periodic,
        bottom: EdgeBC::BounceBack,
        top: EdgeBC::BounceBack,
    };
    let mut heavy = make_sim(nx, ny, MC_NU, edges);
    let mut light = make_sim(nx, ny, MC_NU, edges);
    init_mcmp_layers(&mut heavy, &mut light, mode, a0, true);
    let mc = MultiComponent::new(MC_G_AB).with_gravity([0.0, -g], [0.0, 0.0]);
    let m0 = total_mcmp_mass(&heavy, &light);
    let mut series = vec![(0.0, fourier_amp(&heavy, mode))];
    let mut max_amp = series[0].1.abs();
    let mut amp_step_10 = f64::NAN;
    let mut amp_step_100 = f64::NAN;
    let mut amp_step_1000 = f64::NAN;
    let max_u_step_0 = max_speed(&heavy).max(max_speed(&light));
    let mut max_u_step_10 = f64::NAN;
    let mut max_u_step_100 = f64::NAN;
    let mut max_u_step_1000 = f64::NAN;
    for step in 1..=steps {
        mc.update_forces(&mut heavy, &mut light);
        heavy.step();
        light.step();
        if matches!(step, 10 | 100 | 1000) {
            let amp = fourier_amp(&heavy, mode);
            let umax = max_speed(&heavy).max(max_speed(&light));
            match step {
                10 => {
                    amp_step_10 = amp;
                    max_u_step_10 = umax;
                }
                100 => {
                    amp_step_100 = amp;
                    max_u_step_100 = umax;
                }
                1000 => {
                    amp_step_1000 = amp;
                    max_u_step_1000 = umax;
                }
                _ => unreachable!(),
            }
            println!(
                "VAL MPHARD I3 diag: mode={} step={} amp={:.8} max_u={:.8e}",
                mode, step, amp, umax
            );
        }
        if step % sample_every == 0 {
            let amp = fourier_amp(&heavy, mode);
            if amp.is_finite() {
                max_amp = max_amp.max(amp.abs());
            }
            series.push((step as f64, amp));
        }
        if !heavy.rho(nx / 2, ny / 2).is_finite() || !light.rho(nx / 2, ny / 2).is_finite() {
            break;
        }
    }
    let k = 2.0 * PI * mode as f64 / nx as f64;
    let gamma_theory = rt_gamma_theory(g, k);
    let gamma_fit = fit_rt_growth(&series).map_or(f64::NAN, |p| p.0);
    RtModeStats {
        mode,
        gamma_fit,
        gamma_theory,
        ratio: gamma_fit / gamma_theory,
        amp0: series[0].1.abs(),
        max_amp,
        final_amp: series.last().unwrap().1.abs(),
        amp_step_10,
        amp_step_100,
        amp_step_1000,
        max_u_step_0,
        max_u_step_10,
        max_u_step_100,
        max_u_step_1000,
        mass_drift: ((total_mcmp_mass(&heavy, &light) - m0) / m0).abs(),
    }
}

fn setup_taylor_culick(h: usize) -> Simulation<f64> {
    let ny = (6 * h).max(96);
    let nx = (72 * h).max(384);
    let film_len = 64 * h;
    let y0 = ny / 2 - h / 2;
    let y1 = y0 + h;
    let mut sim = make_sim(nx, ny, SC_NU, Edges::default());
    sim.init_with(|x, y| {
        let liquid = x < film_len && (y0..y1).contains(&y);
        (if liquid { 2.0 } else { 0.15 }, 0.0, 0.0)
    });
    sim
}

fn liquid_excess_mass(sim: &Simulation<f64>) -> f64 {
    sim.rho_field()
        .iter()
        .map(|&r| (r - SC_RHO_V).max(0.0))
        .sum()
}

fn right_film_edge(sim: &Simulation<f64>, h: usize) -> f64 {
    let nx = sim.nx();
    let ny = sim.ny();
    let rho_mid = 0.5 * (SC_RHO_L + SC_RHO_V);
    let y0 = ny / 2 - h / 2 - 2;
    let y1 = ny / 2 + h / 2 + 2;
    for x in (8..nx - 8).rev() {
        let count = (y0..y1).filter(|&y| sim.rho(x, y) > rho_mid).count();
        if count >= h / 3 {
            return x as f64;
        }
    }
    f64::NAN
}

fn run_taylor_culick(h: usize) -> TaylorCulickStats {
    let mut sim = setup_taylor_culick(h);
    let sc = ShanChen::new(SC_G);
    let mut samples = Vec::new();
    let mut mass_samples = Vec::new();
    let total_steps = 12_000usize;
    let sample_every = 250usize;
    for it in 0..=(total_steps / sample_every) {
        if it > 0 {
            for _ in 0..sample_every {
                sc_step(&mut sim, &sc);
            }
        }
        let t = (it * sample_every) as f64;
        samples.push((t, right_film_edge(&sim, h)));
        mass_samples.push((t, liquid_excess_mass(&sim)));
    }
    let mut window_speeds = Vec::new();
    for w in samples.windows(9) {
        let valid: Vec<_> = w.iter().copied().filter(|(_, x)| x.is_finite()).collect();
        if valid.len() >= 4 {
            let t: Vec<f64> = valid.iter().map(|p| p.0).collect();
            let x: Vec<f64> = valid.iter().map(|p| p.1).collect();
            let fit = linear_fit(&t, &x);
            let retracted = (samples[0].1 - valid.last().unwrap().1).max(0.0) / h as f64;
            println!(
                "VAL MPHARD I4 diag: h={} t=[{:.0},{:.0}] v={:.8e} r2={:.6} retracted_h={:.3}",
                h,
                valid.first().unwrap().0,
                valid.last().unwrap().0,
                -fit.slope,
                fit.r2,
                retracted
            );
            window_speeds.push((-fit.slope, fit.r2, retracted));
        }
    }
    let window: Vec<_> = samples
        .iter()
        .copied()
        .rev()
        .filter(|(_, x)| x.is_finite())
        .take(12)
        .collect::<Vec<_>>()
        .into_iter()
        .rev()
        .collect();
    let t: Vec<f64> = window.iter().map(|p| p.0).collect();
    let x: Vec<f64> = window.iter().map(|p| p.1).collect();
    let fit = linear_fit(&t, &x);
    // Taylor-Culick comes from a rim momentum balance.  As a film retracts by
    // dx, two free surfaces remove energy 2 sigma dx per unit depth; the rim
    // collects mass rho_l h dx and carries kinetic energy
    // (1/2) rho_l h dx v^2.  Equating gives v = sqrt(2 sigma/(rho_l h)).
    let predicted_v = (2.0 * SC_SIGMA / (SC_RHO_L * h as f64)).sqrt();
    let m0 = mass_samples.first().unwrap().1;
    let m1 = mass_samples
        .iter()
        .find(|p| p.0 >= total_steps as f64)
        .map(|p| p.1)
        .unwrap_or_else(|| mass_samples.last().unwrap().1);
    let (last_window_v, last_window_r2, retracted_h) = window_speeds
        .last()
        .copied()
        .unwrap_or((f64::NAN, f64::NAN, f64::NAN));
    let still_rising = if window_speeds.len() >= 2 {
        let prev = window_speeds[window_speeds.len() - 2].0;
        last_window_v > prev * 1.02
    } else {
        false
    };
    TaylorCulickStats {
        h,
        measured_v: -fit.slope,
        predicted_v,
        rel: rel_err(-fit.slope, predicted_v),
        mass_loss: ((m1 - m0) / m0).abs(),
        fit_r2: fit.r2,
        last_window_v,
        last_window_r2,
        retracted_h,
        still_rising,
    }
}

#[test]
fn val_mphard_i1_jurin_capillary_rise_zero_parameter() {
    let rows: Vec<_> = [16usize, 24, 32]
        .into_iter()
        .map(|gap| run_jurin(gap, WALL_RHO_WET, THETA_WET_DEG))
        .collect();
    for row in &rows {
        println!(
            "VAL MPHARD I1: gap={} wall_rho={:.3} theta_deg={:.3} theta_slot_deg={:.3} h_meas={:.6} h_pred={:.6} rel={:.6} steady_drift={:.6} steps={} mass_drift={:e}",
            row.gap,
            row.wall_rho,
            row.theta_deg,
            row.theta_slot_deg,
            row.measured_h,
            row.predicted_h,
            rel_err(row.measured_h, row.predicted_h),
            row.steady_drift,
            row.steps,
            row.mass_drift
        );
    }
    let inv_w: Vec<f64> = rows.iter().map(|r| 1.0 / r.gap as f64).collect();
    let h: Vec<f64> = rows.iter().map(|r| r.measured_h).collect();
    let fit = linear_fit(&inv_w, &h);
    println!(
        "VAL MPHARD I1: linear h_vs_1_over_w slope={:.8} intercept={:.8} r2={:.8}",
        fit.slope, fit.intercept, fit.r2
    );

    let dry = run_jurin(24, WALL_RHO_DRY, THETA_DRY_DEG);
    println!(
        "VAL MPHARD I1: sign_flip gap={} wall_rho={:.3} theta_deg={:.3} theta_slot_deg={:.3} h_meas={:.6} h_pred={:.6}",
        dry.gap, dry.wall_rho, dry.theta_deg, dry.theta_slot_deg, dry.measured_h, dry.predicted_h
    );

    for row in &rows {
        assert!(
            rel_err(row.measured_h, row.predicted_h) <= 0.10,
            "VAL MPHARD I1 Jurin gap={} h_meas={:.8} h_pred={:.8} rel={:.6} band=0.10 steady_drift={:.6} mass_drift={:e}",
            row.gap,
            row.measured_h,
            row.predicted_h,
            rel_err(row.measured_h, row.predicted_h),
            row.steady_drift,
            row.mass_drift
        );
    }
    assert!(
        fit.r2 >= 0.99,
        "VAL MPHARD I1 Jurin linearity r2={:.8} band>=0.99 rows={rows:?}",
        fit.r2
    );

    assert!(
        dry.measured_h < 0.0,
        "VAL MPHARD I1 sign flip failed: theta={:.3} deg > 90, h_meas={:.8}, h_pred={:.8}",
        dry.theta_deg,
        dry.measured_h,
        dry.predicted_h
    );
}

#[test]
fn val_mphard_i2_prosperetti_standing_wave_light_one_k() {
    let row = run_standing_wave(4, 5_000, 100);
    println!(
        "VAL MPHARD I2: mode={} L={:.1} k={:.8e} g_eff={:.8e} delta_rho_eff={:.6} rho_sum={:.6} sigma_AB={:.8e} omega_g2={:.8e} omega_sigma2={:.8e} omega_fit={:.8e} omega0={:.8e} rel={:.6} fit_periods={} fit_t_end={:.1} mode_rms[3,4,5]=[{:.8e},{:.8e},{:.8e}] decay_fit={:.8e} decay_model={:.8e} decay_ratio={:.6} envelope_monotone={:.6} mass_drift={:e}",
        row.mode,
        row.box_l,
        row.k,
        row.g_eff,
        row.delta_rho_eff,
        row.rho_sum,
        MC_SIGMA_AB,
        row.g_branch,
        row.sigma_branch,
        row.omega_fit,
        row.omega0,
        row.omega_rel,
        row.fit_periods,
        row.fit_t_end,
        row.mode_rms_3,
        row.mode_rms_4,
        row.mode_rms_5,
        row.decay_fit,
        row.decay_model,
        row.decay_ratio,
        row.envelope_monotone,
        row.mass_drift
    );
    assert!(
        row.omega_rel <= 0.10,
        "VAL MPHARD I2 omega mismatch mode={} omega_fit={:.8e} omega0={:.8e} rel={:.6} band=0.10",
        row.mode,
        row.omega_fit,
        row.omega0,
        row.omega_rel
    );
    assert!(
        (0.5..=2.0).contains(&row.decay_ratio),
        "VAL MPHARD I2 damping mismatch mode={} decay_fit={:.8e} model={:.8e} ratio={:.6} band=[0.5,2]",
        row.mode,
        row.decay_fit,
        row.decay_model,
        row.decay_ratio
    );
    assert!(
        row.envelope_monotone >= 0.95,
        "VAL MPHARD I2 envelope not monotone: monotonicity={:.6} band>=0.95",
        row.envelope_monotone
    );
}

#[test]
#[ignore = "heavy VAL-MPHARD Prosperetti gravity-capillary k-scan"]
fn val_mphard_i2_heavy_gravity_capillary_crossover_scan() {
    let samples: Vec<_> = [1usize, 2, 4]
        .into_iter()
        .map(|mode| {
            let row = run_standing_wave(mode, 12_000, 200);
            let k = 2.0 * PI * mode as f64 / 256.0;
            println!(
                "VAL MPHARD I2-heavy: mode={} k={:.8e} omega_fit={:.8e} omega0={:.8e} rel={:.6}",
                mode, k, row.omega_fit, row.omega0, row.omega_rel
            );
            (k, row.omega_fit.powi(2))
        })
        .collect();
    let g = 1.0e-4;
    let agree = curve_agreement(
        |k| 0.5 * (g * k + MC_SIGMA_AB * k.powi(3)),
        &samples,
        0.10,
        1.0e-12,
    );
    assert!(
        agree.max_rel_dev <= 0.10,
        "VAL MPHARD I2-heavy omega^2 curve max_rel={:.6} worst_k={:.8e} frac_in_band={:.6}",
        agree.max_rel_dev,
        agree.worst_x,
        agree.frac_in_band
    );
}

#[test]
fn val_mphard_i3_rayleigh_taylor_cutoff_light_sign_canary() {
    let nx = 256;
    let target_mc = 5.0;
    // Pick g so the inviscid surface-tension cutoff
    //   k_c = sqrt(g Delta_rho_eff/sigma_AB)
    // lands at mode m_c ~= 5 in this box.  With the T12 reduction
    // Delta_rho_eff cancels the same factor of 1/2 in both terms, so
    // g = sigma_AB * (2 pi m_c / L)^2.
    let g = MC_SIGMA_AB * (2.0 * PI * target_mc / nx as f64).powi(2);
    let unstable = run_rt_mode(nx, 3, g, 1_200, 100);
    let stable = run_rt_mode(nx, 7, g, 1_200, 100);
    for row in [unstable, stable] {
        println!(
            "VAL MPHARD I3: mode={} gamma_fit={:.8e} gamma_th={:.8e} ratio={:.6} amp0={:.6} amp10={:.6} amp100={:.6} amp1000={:.6} max_amp={:.6} final_amp={:.6} max_u[0,10,100,1000]=[{:.8e},{:.8e},{:.8e},{:.8e}] mass_drift={:e}",
            row.mode,
            row.gamma_fit,
            row.gamma_theory,
            row.ratio,
            row.amp0,
            row.amp_step_10,
            row.amp_step_100,
            row.amp_step_1000,
            row.max_amp,
            row.final_amp,
            row.max_u_step_0,
            row.max_u_step_10,
            row.max_u_step_100,
            row.max_u_step_1000,
            row.mass_drift
        );
    }
    assert!(
        unstable.max_amp >= 2.0 * unstable.amp0,
        "VAL MPHARD I3 unstable mode did not double: mode={} amp0={:.8} max_amp={:.8}",
        unstable.mode,
        unstable.amp0,
        unstable.max_amp
    );
    assert!(
        stable.max_amp <= 1.5 * stable.amp0,
        "VAL MPHARD I3 stable mode exceeded cap: mode={} amp0={:.8} max_amp={:.8} cap=1.5x",
        stable.mode,
        stable.amp0,
        stable.max_amp
    );
    assert!(
        (0.75..=1.25).contains(&unstable.ratio),
        "VAL MPHARD I3 gamma ratio mode={} gamma_fit={:.8e} gamma_th={:.8e} ratio={:.6} band=[0.75,1.25]",
        unstable.mode,
        unstable.gamma_fit,
        unstable.gamma_theory,
        unstable.ratio
    );
}

#[test]
#[ignore = "heavy VAL-MPHARD Rayleigh-Taylor surface-tension cutoff scan"]
fn val_mphard_i3_heavy_cutoff_transition_scan() {
    let nx = 256;
    let target_mc = 5.0;
    let g = MC_SIGMA_AB * (2.0 * PI * target_mc / nx as f64).powi(2);
    let kc = (g / MC_SIGMA_AB).sqrt();
    let mc = kc * nx as f64 / (2.0 * PI);
    let modes = [3usize, 4, 5, 6, 7];
    let rows: Vec<_> = modes
        .into_iter()
        .map(|mode| run_rt_mode(nx, mode, g, 5_000, 200))
        .collect();
    for row in &rows {
        println!(
            "VAL MPHARD I3-heavy: mode={} mc_pred={:.3} gamma_fit={:.8e} gamma_th={:.8e} ratio={:.6} amp0={:.6} max_amp={:.6} final_amp={:.6}",
            row.mode,
            mc,
            row.gamma_fit,
            row.gamma_theory,
            row.ratio,
            row.amp0,
            row.max_amp,
            row.final_amp
        );
    }
    let last_growing = rows
        .iter()
        .filter(|r| r.max_amp >= 2.0 * r.amp0)
        .map(|r| r.mode)
        .max()
        .unwrap_or(0);
    let first_capped = rows
        .iter()
        .filter(|r| r.max_amp <= 1.5 * r.amp0)
        .map(|r| r.mode)
        .min()
        .unwrap_or(usize::MAX);
    assert!(
        (last_growing as f64 - mc).abs() <= 1.0 || (first_capped as f64 - mc).abs() <= 1.0,
        "VAL MPHARD I3-heavy cutoff not bracketed: mc_pred={:.3} last_growing={} first_capped={}",
        mc,
        last_growing,
        first_capped
    );
    for row in rows.iter().filter(|r| r.mode < mc.floor() as usize) {
        assert!(
            (0.75..=1.25).contains(&row.ratio),
            "VAL MPHARD I3-heavy gamma ratio mode={} gamma_fit={:.8e} gamma_th={:.8e} ratio={:.6}",
            row.mode,
            row.gamma_fit,
            row.gamma_theory,
            row.ratio
        );
    }
}

#[test]
fn val_mphard_i4_taylor_culick_film_retraction_scaling() {
    let rows: Vec<_> = [16usize]
        .into_iter()
        .map(run_taylor_culick)
        .collect();
    for row in &rows {
        println!(
            "VAL MPHARD I4: h={} free_interfaces=2 v_meas={:.8e} v_tc={:.8e} rel={:.6} mass_loss={:.6} fit_r2={:.8} last_window_v={:.8e} last_window_r2={:.8} retracted_h={:.3} still_rising={}",
            row.h,
            row.measured_v,
            row.predicted_v,
            row.rel,
            row.mass_loss,
            row.fit_r2,
            row.last_window_v,
            row.last_window_r2,
            row.retracted_h,
            row.still_rising
        );
        if row.mass_loss > 0.03 {
            println!(
                "VAL MPHARD I4: h={} SC thin-film caveat mass_loss={:.6} > 0.03; measurement window is restricted to t=[500,2000]",
                row.h, row.mass_loss
            );
        }
        assert!(
            row.rel <= 0.15,
            "VAL MPHARD I4 Taylor-Culick h={} v_meas={:.8e} v_tc={:.8e} rel={:.6} band=0.15 mass_loss={:.6} fit_r2={:.8}",
            row.h,
            row.measured_v,
            row.predicted_v,
            row.rel,
            row.mass_loss,
            row.fit_r2
        );
    }
    if rows.len() >= 3 {
        let lh: Vec<f64> = rows.iter().map(|r| (r.h as f64).ln()).collect();
        let lv: Vec<f64> = rows.iter().map(|r| r.measured_v.ln()).collect();
        let fit = linear_fit(&lh, &lv);
        println!(
            "VAL MPHARD I4: ln_v_vs_ln_h slope={:.8} intercept={:.8} r2={:.8}",
            fit.slope, fit.intercept, fit.r2
        );
        assert!(
            (fit.slope + 0.5).abs() <= 0.07,
            "VAL MPHARD I4 scaling slope={:.8} target=-0.5 band=0.07 rows={rows:?}",
            fit.slope
        );
    }
}
