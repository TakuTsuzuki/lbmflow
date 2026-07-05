//! Validation T12: two-component MCMP separation and Rayleigh-Taylor growth.

use lbm_core::multiphase::MultiComponent;
use lbm_core::prelude::*;
use std::f64::consts::PI;

const TRACE: f64 = 0.05;
const G_AB_QUANT: f64 = 2.6;
const SIGMA_AB_FROZEN_2026_07_05: f64 = 2.86969302e-2;
const SIGMA_MIN: f64 = 2.0e-2;
const SIGMA_MAX: f64 = 3.5e-2;

#[derive(Clone, Copy, Debug)]
struct SeparationStats {
    g_ab: f64,
    contrast: f64,
    mass_drift: f64,
}

#[derive(Clone, Copy, Debug)]
struct RtStats {
    sigma: f64,
    gamma_fit: f64,
    gamma_theory: f64,
    ratio: f64,
    max_amp: f64,
    mass_drift: f64,
    fit_points: usize,
}

fn make(nx: usize, ny: usize, walls_y: bool, nu: f64) -> Simulation<f64> {
    let edges = if walls_y {
        Edges {
            left: EdgeBC::Periodic,
            right: EdgeBC::Periodic,
            bottom: EdgeBC::BounceBack,
            top: EdgeBC::BounceBack,
        }
    } else {
        Edges::default()
    };
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

fn total_mass(a: &Simulation<f64>, b: &Simulation<f64>) -> f64 {
    a.total_mass() + b.total_mass()
}

fn run_separation(g_ab: f64) -> SeparationStats {
    let (nx, ny) = (96, 96);
    let mut a = make(nx, ny, false, 1.0 / 6.0);
    let mut b = make(nx, ny, false, 1.0 / 6.0);
    a.init_with(|x, _| (if x < nx / 2 { 1.0 } else { TRACE }, 0.0, 0.0));
    b.init_with(|x, _| (if x < nx / 2 { TRACE } else { 1.0 }, 0.0, 0.0));
    let mc = MultiComponent::new(g_ab);
    let m0 = total_mass(&a, &b);
    for _ in 0..5_000 {
        mc.update_forces(&mut a, &mut b);
        a.step();
        b.step();
    }
    let left = a.rho(nx / 4, ny / 2);
    let right = a.rho(3 * nx / 4, ny / 2);
    let ratio = left / right.max(1.0e-30);
    SeparationStats {
        g_ab,
        contrast: ratio.max(1.0 / ratio.max(1.0e-30)),
        mass_drift: ((total_mass(&a, &b) - m0) / m0).abs(),
    }
}

fn measure_sigma_ab() -> f64 {
    let n = 128;
    let r0 = 24.0;
    let c = n as f64 / 2.0;
    let mut a = make(n, n, false, 0.1);
    let mut b = make(n, n, false, 0.1);
    a.init_with(|x, y| {
        let d2 = (x as f64 - c).powi(2) + (y as f64 - c).powi(2);
        (if d2 < r0 * r0 { 1.0 } else { TRACE }, 0.0, 0.0)
    });
    b.init_with(|x, y| {
        let d2 = (x as f64 - c).powi(2) + (y as f64 - c).powi(2);
        (if d2 < r0 * r0 { TRACE } else { 1.0 }, 0.0, 0.0)
    });
    let mc = MultiComponent::new(G_AB_QUANT);
    for _ in 0..20_000 {
        mc.update_forces(&mut a, &mut b);
        a.step();
        b.step();
    }

    let cs2 = 1.0 / 3.0;
    let pressure = |x: usize, y: usize| {
        let (ra, rb) = (a.rho(x, y), b.rho(x, y));
        cs2 * (ra + rb + G_AB_QUANT * ra * rb)
    };
    let dp = pressure(n / 2, n / 2) - pressure(4, 4);
    let rho_mid = 0.5 * (a.rho(n / 2, n / 2) + a.rho(4, 4));
    let area = a.rho_field().iter().filter(|&&r| r > rho_mid).count() as f64;
    dp * (area / PI).sqrt()
}

fn fourier_amp(heavy: &Simulation<f64>, tr: f64) -> f64 {
    let nx = heavy.nx();
    let ny = heavy.ny();
    let k = 2.0 * PI / nx as f64;
    let mut re = 0.0;
    let mut im = 0.0;
    for x in 0..nx {
        let m: f64 = (1..ny - 1).map(|y| heavy.rho(x, y)).sum();
        let ph = k * x as f64;
        re += m * ph.cos();
        im += m * ph.sin();
    }
    let mode = (re * re + im * im).sqrt() * 2.0 / nx as f64;
    mode / (1.0 - tr)
}

fn fit_gamma(series: &[(f64, f64)]) -> Option<(f64, usize)> {
    let mut best = (0usize, 0usize);
    let mut start = 0usize;
    for i in 1..=series.len() {
        let ok = i < series.len()
            && series[i].1 > series[i - 1].1
            && (1.0..=10.0).contains(&series[i].1);
        if !ok {
            if i - start > best.1 - best.0 {
                best = (start, i);
            }
            start = i;
        }
    }
    let pts: Vec<_> = series[best.0..best.1]
        .iter()
        .filter(|&&(_, a)| a > 0.0)
        .map(|&(t, a)| (t, a.ln()))
        .collect();
    if pts.len() < 3 {
        return None;
    }
    let n = pts.len() as f64;
    let sx = pts.iter().map(|p| p.0).sum::<f64>() / n;
    let sy = pts.iter().map(|p| p.1).sum::<f64>() / n;
    let sxy = pts.iter().map(|p| (p.0 - sx) * (p.1 - sy)).sum::<f64>();
    let sxx = pts.iter().map(|p| (p.0 - sx) * (p.0 - sx)).sum::<f64>();
    Some((sxy / sxx, pts.len()))
}

fn run_rt(nx: usize, ny: usize, steps: usize, sample_every: usize, sigma: f64) -> RtStats {
    let nu = 0.1;
    let g = 1.0e-4;
    let y0 = ny as f64 / 2.0;
    let a0 = 6.0;
    let k = 2.0 * PI / nx as f64;
    let interface = |x: usize| y0 + a0 * (k * x as f64).cos();
    let mut heavy = make(nx, ny, true, nu);
    let mut light = make(nx, ny, true, nu);
    heavy.init_with(|x, y| {
        let up = (y as f64) > interface(x);
        (if up { 1.0 } else { TRACE }, 0.0, 0.0)
    });
    light.init_with(|x, y| {
        let up = (y as f64) > interface(x);
        (if up { TRACE } else { 1.0 }, 0.0, 0.0)
    });
    let mc = MultiComponent::new(G_AB_QUANT).with_gravity([0.0, -g], [0.0, 0.0]);
    let m0 = total_mass(&heavy, &light);
    let mut series = Vec::new();
    let mut max_amp = fourier_amp(&heavy, TRACE);
    for it in 1..=(steps / sample_every) {
        for _ in 0..sample_every {
            mc.update_forces(&mut heavy, &mut light);
            heavy.step();
            light.step();
        }
        let amp = fourier_amp(&heavy, TRACE);
        assert!(
            amp.is_finite(),
            "T12 RT non-finite amp at step {}, nx = {nx}, ny = {ny}",
            it * sample_every
        );
        max_amp = max_amp.max(amp);
        series.push(((it * sample_every) as f64, amp));
    }
    let gamma0_sq = 0.5 * g * k - 0.5 * sigma * k.powi(3);
    let gamma_theory = (gamma0_sq + nu * nu * k.powi(4)).sqrt() - nu * k * k;
    let (gamma_fit, fit_points) = fit_gamma(&series).unwrap_or((f64::NAN, 0));
    RtStats {
        sigma,
        gamma_fit,
        gamma_theory,
        ratio: gamma_fit / gamma_theory,
        max_amp,
        mass_drift: ((total_mass(&heavy, &light) - m0) / m0).abs(),
        fit_points,
    }
}

#[test]
fn t12_mcmp_separation_threshold_smoke() {
    let separated = run_separation(2.2);
    let mixed = run_separation(1.8);
    assert!(
        separated.contrast >= 3.0,
        "T12 G_ab = {}, contrast = {:.6}, mass_drift = {:e}",
        separated.g_ab,
        separated.contrast,
        separated.mass_drift
    );
    assert!(
        mixed.contrast <= 1.5,
        "T12 G_ab = {}, contrast = {:.6}, mass_drift = {:e}",
        mixed.g_ab,
        mixed.contrast,
        mixed.mass_drift
    );
    for s in [separated, mixed] {
        assert!(
            s.mass_drift <= 1.0e-10,
            "T12 G_ab = {}, relative mass drift = {:e}, contrast = {:.6}",
            s.g_ab,
            s.mass_drift,
            s.contrast
        );
    }
}

#[test]
fn t12_mcmp_sigma_ab_laplace_regression() {
    let sigma = measure_sigma_ab();
    eprintln!("T12 measured sigma_AB = {sigma:.8e}");
    assert!(
        (SIGMA_MIN..=SIGMA_MAX).contains(&sigma),
        "T12 sigma_AB = {:.8e}, expected [{:.1e}, {:.1e}]",
        sigma,
        SIGMA_MIN,
        SIGMA_MAX
    );
}

#[test]
fn t12_rt_default_growth_smoke_reaches_amp_8() {
    let sigma = SIGMA_AB_FROZEN_2026_07_05;
    let stats = run_rt(128, 128, 8_000, 500, sigma);
    eprintln!(
        "T12 smoke sigma={:.8e} max_amp={:.6} mass_drift={:e}",
        stats.sigma, stats.max_amp, stats.mass_drift
    );
    assert!(
        stats.max_amp >= 8.0,
        "T12 smoke max_amp = {:.6}, sigma = {:.8e}, mass_drift = {:e}",
        stats.max_amp,
        stats.sigma,
        stats.mass_drift
    );
    assert!(
        stats.mass_drift <= 1.0e-10,
        "T12 smoke mass_drift = {:e}, max_amp = {:.6}, sigma = {:.8e}",
        stats.mass_drift,
        stats.max_amp,
        stats.sigma
    );
}

#[test]
#[ignore = "256^2 x 12k RT rate fit is intentionally outside the default runtime budget"]
fn t12_rt_growth_rate_matches_corrected_dispersion() {
    let sigma = measure_sigma_ab();
    let stats = run_rt(256, 256, 12_000, 500, sigma);
    eprintln!(
        "T12 full sigma={:.8e} gamma_fit={:.8e} gamma_th={:.8e} ratio={:.6} max_amp={:.6} mass_drift={:e} fit_points={}",
        stats.sigma,
        stats.gamma_fit,
        stats.gamma_theory,
        stats.ratio,
        stats.max_amp,
        stats.mass_drift,
        stats.fit_points
    );
    assert!(
        stats.max_amp > 10.0,
        "T12 RT max_amp = {:.6}, sigma = {:.8e}, mass_drift = {:e}",
        stats.max_amp,
        stats.sigma,
        stats.mass_drift
    );
    assert!(
        stats.mass_drift <= 1.0e-10,
        "T12 RT mass_drift = {:e}, max_amp = {:.6}, sigma = {:.8e}",
        stats.mass_drift,
        stats.max_amp,
        stats.sigma
    );
    assert!(
        (0.75..=1.25).contains(&stats.ratio),
        "T12 RT gamma_fit = {:.8e}, gamma_th = {:.8e}, ratio = {:.6}, sigma = {:.8e}, max_amp = {:.6}, fit_points = {}",
        stats.gamma_fit,
        stats.gamma_theory,
        stats.ratio,
        stats.sigma,
        stats.max_amp,
        stats.fit_points
    );
}
