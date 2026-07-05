//! Dev probe: two-component Shan-Chen — separation sanity + Rayleigh-Taylor
//! growth rate vs linear theory gamma = sqrt(A g k).

use lbm_core::multiphase::MultiComponent;
use lbm_core::prelude::*;

fn make(nx: usize, ny: usize, walls_y: bool, nu: f64) -> Simulation<f64> {
    let e = if walls_y {
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
        edges: e,
        ..Default::default()
    }
    .build()
    .unwrap()
}

fn separation_smoke(g_ab: f64, trace: f64) {
    let (nx, ny) = (96, 96);
    let mut a = make(nx, ny, false, 1.0 / 6.0);
    let mut b = make(nx, ny, false, 1.0 / 6.0);
    a.init_with(|x, _| (if x < nx / 2 { 1.0 } else { trace }, 0.0, 0.0));
    b.init_with(|x, _| (if x < nx / 2 { trace } else { 1.0 }, 0.0, 0.0));
    let mc = MultiComponent::new(g_ab);
    let m0 = a.total_mass() + b.total_mass();
    for _ in 0..5000 {
        mc.update_forces(&mut a, &mut b);
        a.step();
        b.step();
    }
    let m1 = a.total_mass() + b.total_mass();
    let a_left = a.rho(nx / 4, ny / 2);
    let a_right = a.rho(3 * nx / 4, ny / 2);
    let umax = a
        .ux_field()
        .iter()
        .chain(a.uy_field())
        .fold(0.0f64, |acc, v| acc.max(v.abs()));
    println!(
        "[sep G={g_ab} tr={trace}] rhoA left={a_left:.4} right={a_right:.4} contrast={:.1} mass_drift={:.1e} max|u|={umax:.1e} finite={}",
        a_left / a_right.max(1e-9),
        ((m1 - m0) / m0).abs(),
        m1.is_finite()
    );
}

/// Measure sigma_AB via an MCMP droplet: p = cs^2 (rhoA + rhoB + G rhoA rhoB),
/// dp = sigma / R.
fn mcmp_sigma(g_ab: f64, tr: f64) -> f64 {
    let n = 128;
    let mut a = make(n, n, false, 0.1);
    let mut b = make(n, n, false, 0.1);
    let c = n as f64 / 2.0;
    let r0 = 24.0;
    a.init_with(|x, y| {
        let d2 = (x as f64 - c).powi(2) + (y as f64 - c).powi(2);
        (if d2 < r0 * r0 { 1.0 } else { tr }, 0.0, 0.0)
    });
    b.init_with(|x, y| {
        let d2 = (x as f64 - c).powi(2) + (y as f64 - c).powi(2);
        (if d2 < r0 * r0 { tr } else { 1.0 }, 0.0, 0.0)
    });
    let mc = MultiComponent::new(g_ab);
    for _ in 0..20_000 {
        mc.update_forces(&mut a, &mut b);
        a.step();
        b.step();
    }
    let cs2 = 1.0 / 3.0;
    let p = |x: usize, y: usize| {
        let (ra, rb) = (a.rho(x, y), b.rho(x, y));
        cs2 * (ra + rb + g_ab * ra * rb)
    };
    let dp = p(n / 2, n / 2) - p(4, 4);
    // effective radius from heavy-phase area
    let rho_mid = 0.5 * (a.rho(n / 2, n / 2) + a.rho(4, 4));
    let area = a.rho_field().iter().filter(|&&r| r > rho_mid).count() as f64;
    let r_fit = (area / std::f64::consts::PI).sqrt();
    let sigma = dp * r_fit;
    println!("[mcmp sigma] G={g_ab}: dp={dp:.4e} r={r_fit:.1} sigma={sigma:.4e}");
    sigma
}

/// Matched-density RT: both components rho=1 in their bulk (trace `tr`),
/// gravity only on the heavy component. Reference growth rate from the
/// dispersion relation with tension + viscosity corrections:
/// gamma = sqrt(g k / 2 - sigma k^3 / 2 + nu^2 k^4) - nu k^2.
fn rt_growth(g: f64, g_ab: f64, tr: f64, sigma: f64, label: &str) {
    let (nx, ny) = (256, 256);
    let nu = 0.1;
    let mut heavy = make(nx, ny, true, nu);
    let mut light = make(nx, ny, true, nu);
    let y0 = ny as f64 / 2.0;
    let a0 = 6.0;
    let k = 2.0 * std::f64::consts::PI / nx as f64;
    let interface = |x: usize| y0 + a0 * (k * x as f64).cos();
    heavy.init_with(|x, y| {
        let up = (y as f64) > interface(x);
        (if up { 1.0 } else { tr }, 0.0, 0.0)
    });
    light.init_with(|x, y| {
        let up = (y as f64) > interface(x);
        (if up { tr } else { 1.0 }, 0.0, 0.0)
    });
    let mc = MultiComponent::new(g_ab).with_gravity([0.0, -g], [0.0, 0.0]);

    // amplitude via Fourier projection of the per-column heavy mass:
    // M(x) = sum_y rho_A ~ linear in interface height, so |k-mode of M|
    // recovers the perturbation amplitude robustly (no contour glitches).
    let amp = |h: &Simulation<f64>| -> f64 {
        let mut re = 0.0;
        let mut im = 0.0;
        for x in 0..nx {
            let m: f64 = (1..ny - 1).map(|y| h.rho(x, y)).sum();
            let ph = k * x as f64;
            re += m * ph.cos();
            im += m * ph.sin();
        }
        let mode = (re * re + im * im).sqrt() * 2.0 / nx as f64;
        // interface height amplitude = mode / (rho_heavy - rho_trace)
        mode / (1.0 - tr)
    };

    let gamma0_sq = 0.5 * g * k - 0.5 * sigma * k * k * k;
    let gamma_theory = if gamma0_sq > 0.0 {
        (gamma0_sq + nu * nu * k.powi(4)).sqrt() - nu * k * k
    } else {
        f64::NAN
    };
    print!("[rt {label}] amp every 500 steps: ");
    let mut series: Vec<(f64, f64)> = Vec::new();
    for it in 0..24 {
        for _ in 0..500 {
            mc.update_forces(&mut heavy, &mut light);
            heavy.step();
            light.step();
        }
        let a = amp(&heavy);
        print!("{a:.1} ");
        if !a.is_finite() {
            break;
        }
        series.push((((it + 1) as f64) * 500.0, a));
    }
    println!();
    // exponential fit on the longest monotonically-increasing run with
    // amplitudes in [1.0, 10.0] (the linear regime before k-mode saturation)
    let mut best: (usize, usize) = (0, 0);
    let mut start = 0;
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
    let pts: Vec<(f64, f64)> = series[best.0..best.1]
        .iter()
        .filter(|&&(_, a)| a > 0.0)
        .map(|&(t, a)| (t, a.ln()))
        .collect();
    if pts.len() >= 3 {
        let n = pts.len() as f64;
        let sx: f64 = pts.iter().map(|p| p.0).sum::<f64>() / n;
        let sy: f64 = pts.iter().map(|p| p.1).sum::<f64>() / n;
        let sxy: f64 = pts.iter().map(|p| (p.0 - sx) * (p.1 - sy)).sum();
        let sxx: f64 = pts.iter().map(|p| (p.0 - sx) * (p.0 - sx)).sum();
        let gamma_fit = sxy / sxx;
        println!(
            "[rt {label}] gamma_fit={gamma_fit:.4e} theory={gamma_theory:.4e} ratio={:.3}",
            gamma_fit / gamma_theory
        );
    } else {
        println!("[rt {label}] not enough points in the linear window; theory={gamma_theory:.4e}");
    }
}

fn main() {
    let sigma = mcmp_sigma(2.6, 0.05);
    rt_growth(1e-4, 2.6, 0.05, sigma, "g=1e-4 G=2.6 L=256");
}
