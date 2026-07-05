//! ACC-AUDIT: adversarial accuracy probes for physics-bending approximations.
//!
//! These tests intentionally look beyond the steady-state validation bands in
//! `docs/VALIDATION.md`: transient wall layers, dispersion, Galilean defects,
//! compressibility scaling, staircase sensitivity, slip laws, and force-path
//! consistency. Each reference formula is derived locally in comments so the
//! test is reviewable without consulting implementation internals.

use lbm_core::compat::prelude::*;
use lbm_core::compat::rotor::Rotor;
use std::f64::consts::{PI, SQRT_2};

const TRT: Collision = Collision::Trt {
    magic: Collision::MAGIC_STD,
};

fn l2_rel(actual: &[f64], reference: &[f64]) -> f64 {
    let (mut num, mut den) = (0.0, 0.0);
    for (&a, &r) in actual.iter().zip(reference) {
        num += (a - r) * (a - r);
        den += r * r;
    }
    (num / den).sqrt()
}

fn assert_close(value: f64, expected: f64, tol: f64, label: &str) {
    assert!(
        (value - expected).abs() <= tol,
        "{label}: value={value:.12e}, expected={expected:.12e}, tol={tol:.3e}"
    );
}

// Abramowitz-Stegun 7.1.26. Maximum absolute error is ~1.5e-7, far below the
// audit's measured LBM errors. The Stokes reference uses erfc(z) for z >= 0.
fn erf_approx(x: f64) -> f64 {
    let sign = if x < 0.0 { -1.0 } else { 1.0 };
    let x = x.abs();
    let t = 1.0 / (1.0 + 0.3275911 * x);
    let poly = (((((1.061405429 * t - 1.453152027) * t) + 1.421413741) * t - 0.284496736) * t
        + 0.254829592)
        * t;
    sign * (1.0 - poly * (-x * x).exp())
}

fn erfc_approx(x: f64) -> f64 {
    1.0 - erf_approx(x)
}

fn channel_periodic_x(
    ny: usize,
    nu: f64,
    bottom: EdgeBC<f64>,
    top: EdgeBC<f64>,
) -> Simulation<f64> {
    SimConfig {
        nx: 8,
        ny,
        nu,
        collision: TRT,
        edges: Edges {
            left: EdgeBC::Periodic,
            right: EdgeBC::Periodic,
            bottom,
            top,
        },
        ..Default::default()
    }
    .build()
    .unwrap()
}

fn stokes_first_error(ny: usize, steps: usize) -> f64 {
    let u_wall = 0.02;
    let nu = 0.1;
    let mut sim = channel_periodic_x(
        ny,
        nu,
        EdgeBC::MovingWall { u: [u_wall, 0.0] },
        EdgeBC::BounceBack,
    );
    sim.init_with(|_, _| (1.0, 0.0, 0.0));
    sim.run(steps);

    let layer = 4.0 * (nu * steps as f64).sqrt();
    let mut actual = Vec::new();
    let mut reference = Vec::new();
    for y in 1..ny - 1 {
        let y_w = y as f64 - 0.5;
        if y_w > layer {
            continue;
        }
        // Stokes first problem is diffusion on a half-space:
        // du/dt = nu d2u/dy2, u(0,t)=U, u(y,0)=0. Similarity variable
        // eta = y/(2 sqrt(nu t)) gives u/U = erfc(eta).
        reference.push(u_wall * erfc_approx(y_w / (2.0 * (nu * steps as f64).sqrt())));
        actual.push(sim.ux(4, y));
    }
    l2_rel(&actual, &reference)
}

#[test]
fn stokes_first_impulsive_plate_light_transient_profile() {
    let e100 = stokes_first_error(96, 100);
    let e400 = stokes_first_error(96, 400);
    let ratio = e100 / e400;
    println!("ACC Stokes-I light: e100={e100:.6e}, e400={e400:.6e}, ratio={ratio:.3}");
    assert!(e100 <= 8.0e-2, "Stokes-I t=100 L2rel={e100:e}");
    assert!(e400 <= 6.0e-2, "Stokes-I t=400 L2rel={e400:e}");
    assert!(
        (2.5..=5.5).contains(&ratio),
        "Stokes-I error ratio e100/e400={ratio:e}"
    );
}

#[test]
#[ignore = "heavy ACC-AUDIT full Stokes-I ny=192, t={200,800,3200}"]
fn stokes_first_impulsive_plate_full_spec() {
    let e200 = stokes_first_error(192, 200);
    let e800 = stokes_first_error(192, 800);
    let e3200 = stokes_first_error(192, 3200);
    println!("ACC Stokes-I full: e200={e200:.6e}, e800={e800:.6e}, e3200={e3200:.6e}");
    assert!(e200 <= 8.0e-2);
    assert!(e800 <= 6.0e-2);
    assert!(e3200 <= 6.0e-2);
}

#[test]
#[ignore = "SPEC-GAP: compat API has no runtime setter for MovingWall velocity; oscillating-wall probe must land with that contract"]
fn stokes_second_oscillating_wall_spec_gap() {
    // Analytic derivation for the future implementation:
    // For wall velocity U sin(omega t), the linear Stokes solution is
    // u(y,t) = U exp(-k y) sin(omega t - k y), where k = sqrt(omega/(2 nu)).
    // With omega = 2*pi/T, k = sqrt(pi/(nu T)). The fitted amplitude envelope
    // must be U exp(-k y_w), and phase lag must be k y_w.
}

fn density_mode_coeff(sim: &Simulation<f64>) -> (f64, f64) {
    let nx = sim.nx();
    let k = 2.0 * PI / nx as f64;
    let mut s = 0.0;
    let mut c = 0.0;
    for x in 0..nx {
        let rho = sim.rho(x, 2) - 1.0;
        s += rho * (k * x as f64).sin();
        c += rho * (k * x as f64).cos();
    }
    (2.0 * s / nx as f64, 2.0 * c / nx as f64)
}

#[test]
fn acoustic_sound_speed_and_damping_periodic_density_wave() {
    let (nx, ny) = (512, 4);
    let nu = 0.02;
    let amp = 1.0e-4;
    let k = 2.0 * PI / nx as f64;
    let cs = 1.0 / 3.0_f64.sqrt();
    let period = nx as f64 / cs;
    let steps = (4.0 * period).round() as usize;
    let mut sim = SimConfig {
        nx,
        ny,
        nu,
        collision: TRT,
        ..Default::default()
    }
    .build()
    .unwrap();
    sim.init_with(|x, _| (1.0 + amp * (k * x as f64).sin(), 0.0, 0.0));

    let mut crossings = Vec::new();
    let mut amps = Vec::new();
    let mut prev = density_mode_coeff(&sim).0;
    for t in 1..=steps {
        sim.step();
        let (a_sin, a_cos) = density_mode_coeff(&sim);
        let a = (a_sin * a_sin + a_cos * a_cos).sqrt();
        amps.push((t as f64, a));
        if t > 2 {
            let cur = a_sin;
            if prev > 0.0 && cur <= 0.0 {
                crossings.push(t as f64);
            }
            prev = cur;
        }
    }
    let measured_period =
        (crossings.last().unwrap() - crossings.first().unwrap()) / (crossings.len() - 1) as f64;
    let c_measured = nx as f64 / measured_period;
    let c_rel = (c_measured / cs - 1.0).abs();

    let mut maxima = Vec::new();
    for w in amps.windows(3) {
        if w[1].1 >= w[0].1 && w[1].1 >= w[2].1 {
            maxima.push(w[1]);
        }
    }
    let a0 = maxima.first().unwrap().1;
    let a1 = maxima.last().unwrap().1;
    let dt = maxima.last().unwrap().0 - maxima.first().unwrap().0;
    let gamma_measured = (a0 / a1).ln() / dt;
    // Linearized isothermal BGK/TRT acoustics has longitudinal attenuation
    // Gamma = nu + zeta/2 in 2D. The D2Q9 bulk viscosity is zeta = 2 nu / 3,
    // hence modal amplitude decays as exp[-(4/3) nu k^2 t].
    let gamma_ref = (4.0 / 3.0) * nu * k * k;
    let gamma_rel = (gamma_measured / gamma_ref - 1.0).abs();
    println!(
        "ACC acoustics: c={c_measured:.9e}, rel={c_rel:.3e}, gamma={gamma_measured:.9e}, ref={gamma_ref:.9e}, rel={gamma_rel:.3e}"
    );
    assert!(c_rel <= 1.0e-3, "sound speed rel={c_rel:e}");
    assert!(gamma_rel <= 3.0e-1, "acoustic damping rel={gamma_rel:e}");
}

fn tgv_sim(n: usize, nu: f64, u0: f64, u_adv: f64, collision: Collision) -> Simulation<f64> {
    let k = 2.0 * PI / n as f64;
    let mut sim = SimConfig {
        nx: n,
        ny: n,
        nu,
        collision,
        ..Default::default()
    }
    .build()
    .unwrap();
    sim.init_with(|x, y| {
        let xf = k * x as f64;
        let yf = k * y as f64;
        let rho = 1.0 - 3.0 * u0 * u0 / 4.0 * ((2.0 * xf).cos() + (2.0 * yf).cos());
        (
            rho,
            u_adv - u0 * xf.cos() * yf.sin(),
            u0 * xf.sin() * yf.cos(),
        )
    });
    sim
}

fn galilean_defect(collision: Collision) -> f64 {
    let n = 64;
    let u0 = 0.02;
    let u_adv = 0.05;
    let steps = 1280usize;
    let shift = (u_adv * steps as f64).round() as usize % n;
    assert_eq!(
        shift, 0,
        "chosen Galilean shift must be an integer cell wrap"
    );
    let nu = n as f64 * n as f64 / (8.0 * PI * PI * steps as f64);
    let mut static_frame = tgv_sim(n, nu, u0, 0.0, collision);
    let mut moving_frame = tgv_sim(n, nu, u0, u_adv, collision);
    static_frame.run(steps);
    moving_frame.run(steps);

    let mut actual = Vec::with_capacity(2 * n * n);
    let mut reference = Vec::with_capacity(2 * n * n);
    for y in 0..n {
        for x in 0..n {
            let sx = (x + n - shift) % n;
            actual.push(moving_frame.ux(x, y) - u_adv);
            actual.push(moving_frame.uy(x, y));
            reference.push(static_frame.ux(sx, y));
            reference.push(static_frame.uy(sx, y));
        }
    }
    l2_rel(&actual, &reference)
}

#[test]
fn galilean_invariance_tgv_defect_bgk_and_trt() {
    let bgk = galilean_defect(Collision::Bgk);
    let trt = galilean_defect(TRT);
    println!("ACC Galilean TGV defect: BGK={bgk:.6e}, TRT={trt:.6e}");
    assert!(bgk <= 1.5e-1, "BGK Galilean defect={bgk:e}");
    assert!(trt <= 1.5e-1, "TRT Galilean defect={trt:e}");
}

fn cavity_profile(n: usize, u_lid: f64, nu: f64, steps: usize) -> Vec<f64> {
    let mut sim = SimConfig {
        nx: n,
        ny: n,
        nu,
        collision: TRT,
        edges: Edges {
            left: EdgeBC::BounceBack,
            right: EdgeBC::BounceBack,
            bottom: EdgeBC::BounceBack,
            top: EdgeBC::MovingWall { u: [u_lid, 0.0] },
        },
        ..Default::default()
    }
    .build()
    .unwrap();
    sim.run(steps);
    let cx = n / 2;
    (1..n - 1).map(|y| sim.ux(cx, y) / u_lid).collect()
}

#[test]
fn cavity_same_re_half_mach_scaling_light() {
    let n = 65;
    let l = (n - 2) as f64;
    let re = 100.0;
    let p_fast = cavity_profile(n, 0.1, 0.1 * l / re, 8_000);
    let p_slow = cavity_profile(n, 0.05, 0.05 * l / re, 16_000);
    let diff = l2_rel(&p_fast, &p_slow);
    // Same Re with Ma halved should leave the incompressible solution fixed;
    // the leading compressibility defect is O(Ma^2). cs = 1/sqrt(3), so the
    // fast-vs-slow expected scale is (0.1/cs)^2 - (0.05/cs)^2 = 0.0225.
    println!("ACC cavity Re=100 same-Re half-Ma profile diff={diff:.6e}");
    assert!(diff <= 1.8e-1, "same-Re half-Ma cavity diff={diff:e}");
}

#[test]
#[ignore = "heavy ACC-AUDIT full cavity N=129, nu={0.127,0.0635}"]
fn cavity_same_re_half_mach_scaling_full_spec() {
    let p_fast = cavity_profile(129, 0.1, 0.127, 80_000);
    let p_slow = cavity_profile(129, 0.05, 0.0635, 160_000);
    let diff = l2_rel(&p_fast, &p_slow);
    println!("ACC cavity full same-Re half-Ma profile diff={diff:.6e}");
    assert!(diff <= 8.0e-2, "full same-Re half-Ma cavity diff={diff:e}");
}

fn rotor_mean_torque(center: [f64; 2], steps: usize) -> f64 {
    let mut sim = SimConfig {
        nx: 129,
        ny: 129,
        nu: 0.02,
        collision: TRT,
        edges: Edges {
            left: EdgeBC::BounceBack,
            right: EdgeBC::BounceBack,
            bottom: EdgeBC::BounceBack,
            top: EdgeBC::BounceBack,
        },
        ..Default::default()
    }
    .build()
    .unwrap();
    let omega = 0.05 / 28.0;
    let mut rotor = Rotor::new(center[0], center[1])
        .n_blades(4)
        .r_hub(4.0)
        .r_blade(28.0)
        .blade_thickness(3.0)
        .omega(omega)
        .chi(1.0)
        .omega_ramp_steps(200);
    let mut sum = 0.0;
    let mut count = 0usize;
    for t in 0..steps {
        sim.force_field_mut().fill([0.0; 2]);
        rotor.update_force(&mut sim);
        sim.step();
        if t >= steps / 2 {
            sum += rotor.torque();
            count += 1;
        }
    }
    sum / count as f64
}

#[test]
fn rotor_half_cell_translation_staircase_sensitivity_light() {
    let a = rotor_mean_torque([64.0, 64.0], 2_000);
    let b = rotor_mean_torque([64.5, 64.5], 2_000);
    let rel = ((a - b) / (0.5 * (a.abs() + b.abs()))).abs();
    println!("ACC rotor light: torque_center={a:.6e}, torque_half={b:.6e}, rel={rel:.3e}");
    assert!(
        rel <= 5.0e-1,
        "rotor half-cell staircase sensitivity={rel:e}"
    );
}

#[test]
#[ignore = "heavy ACC-AUDIT full rotor 20k steps"]
fn rotor_half_cell_translation_staircase_sensitivity_full_spec() {
    let a = rotor_mean_torque([64.0, 64.0], 20_000);
    let b = rotor_mean_torque([64.5, 64.5], 20_000);
    let rel = ((a - b) / (0.5 * (a.abs() + b.abs()))).abs();
    println!("ACC rotor full: torque_center={a:.6e}, torque_half={b:.6e}, rel={rel:.3e}");
    assert!(rel <= 3.0e-1, "full rotor half-cell sensitivity={rel:e}");
}

fn poiseuille_bgk_slip_offset(tau: f64) -> f64 {
    let h = 16.0;
    let ny = h as usize + 2;
    let nu = (tau - 0.5) / 3.0;
    let g = 1.0e-6;
    let mut sim = SimConfig {
        nx: 32,
        ny,
        nu,
        collision: Collision::Bgk,
        force: [g, 0.0],
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
    sim.run(20_000);
    let a = g / (2.0 * nu);
    let mut offsets = Vec::new();
    for y in 1..ny - 1 {
        let yw = y as f64 - 0.5;
        let parabolic = a * yw * (h - yw);
        offsets.push(sim.ux(16, y) - parabolic);
    }
    offsets.iter().sum::<f64>() / offsets.len() as f64
}

#[test]
#[ignore = "SPEC-GAP: freeze BGK half-way bounce-back slip law after choosing the exact Ginzburg/d'Humieres convention for this force/discrete-time path"]
fn bgk_poiseuille_slip_matches_theory_spec_gap() {
    // Probe values are available today through `poiseuille_bgk_slip_offset`.
    // Theory to pin: half-way bounce-back with BGK has a tau-dependent
    // apparent-wall displacement / slip term proportional to the body-force
    // curvature. The common TRT parametrization cancels the term at Lambda =
    // 3/16. This ignored test is deliberately a SPEC-GAP until the project
    // records the exact convention (force path, channel-height convention, and
    // fitted observable) to avoid freezing the wrong closed-form coefficient.
    let _ = [
        poiseuille_bgk_slip_offset(0.7),
        poiseuille_bgk_slip_offset(1.0),
        poiseuille_bgk_slip_offset(1.5),
    ];
}

fn one_step_momentum_gain_force_path(path: &str) -> f64 {
    let (nx, ny) = (16, 16);
    let f = [1.0e-7, -2.0e-7];
    let mut sim = SimConfig {
        nx,
        ny,
        nu: 1.0 / 6.0,
        collision: TRT,
        force: if path == "uniform" { f } else { [0.0, 0.0] },
        ..Default::default()
    }
    .build()
    .unwrap();
    if path == "field" {
        sim.force_field_mut().fill(f);
    } else if path == "gravity" {
        sim.set_gravity(f);
    }
    let p0 = sim.total_momentum();
    sim.step();
    let p1 = sim.total_momentum();
    let gain = [p1[0] - p0[0], p1[1] - p0[1]];
    gain[0] / (nx * ny) as f64
}

#[test]
fn forcing_path_gravity_and_force_field_one_step_match() {
    let field = one_step_momentum_gain_force_path("field");
    let gravity = one_step_momentum_gain_force_path("gravity");
    println!("ACC force path: force_field_gain={field:.12e}, gravity_gain={gravity:.12e}");
    assert_close(
        gravity,
        field,
        1.0e-14,
        "gravity vs force_field one-step gain",
    );
}

#[test]
#[ignore = "expected failure until R2-C fixes ANOM-P2-001; then this current-wrong-value pin must fail loudly and be retightened"]
fn uniform_force_impulse_current_wrong_value_pin_anom_p2_001() {
    let uniform = one_step_momentum_gain_force_path("uniform");
    let field = one_step_momentum_gain_force_path("field");
    // ANOM-P2-001 calibration: at tau=1 TRT, the uniform-force path currently
    // injects a different step-1 impulse than the per-cell force field /
    // gravity path, even though steady slopes match. In this observable
    // (post-step momentum gain after subtracting each path's own half-force
    // initial momentum), the current wrong uniform/field ratio is 7/3. This
    // assertion pins that wrong value so a correct R2-C implementation breaks
    // the test and forces the band to be tightened to equality.
    let ratio = uniform / field;
    println!("ACC force path ANOM-P2-001: uniform={uniform:.12e}, field={field:.12e}, ratio={ratio:.12e}");
    assert_close(
        ratio,
        7.0 / 3.0,
        2.0e-2,
        "current wrong uniform/field impulse ratio",
    );
}

#[test]
#[ignore = "native 3D Solver probe reserved for ACC-AUDIT extension once a non-steady 3D analytic target is selected"]
fn native_3d_accuracy_audit_placeholder_spec_gap() {
    // SPEC-GAP: The order requests native 3D Solver use "where noted", but no
    // specific 3D item in this audit list identifies a native-only analytic
    // target. Keep this ignored marker so the gap is explicit instead of
    // silently claiming 3D coverage from unrelated T15 steady tests.
    let _ = SQRT_2;
}
