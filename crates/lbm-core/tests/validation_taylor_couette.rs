//! VAL-TAYC / V&V Axis 9.8: Taylor-Couette wavy-vortex onset.
//!
//! Physics anchor: concentric cylinders with stationary outer wall and rotating
//! inner wall. The laminar circular-Couette solution is
//! `u_theta(r) = a r + b/r`. For the narrow-gap onset scan we use
//! `Ta = Omega_i^2 R_i (R_o - R_i)^3 / nu^2` and the textbook critical value
//! `Ta_c ~= 3390` as an onset bracket, not an exact wavelength/growth-rate
//! benchmark. Exact transition curves need larger aspect-ratio domains.
//!
//! Caveat: ANOM-P4-010 documents divergence for filled solid-disc compat
//! rotor penalization. This test deliberately uses one-cell THIN CYLINDRICAL
//! SHELLS, not filled discs, so it stays in the marginally stable rotating
//! boundary regime until the P4-001/P4-010 core fixes land.

use lbm_core::prelude::*;
use std::f64::consts::TAU;

type TaycSolver = Solver<D3Q19, f64, CpuSimd, LocalPeriodic>;

const NX: usize = 128;
const NY: usize = 64;
const NZ: usize = 64;
const CENTER_X: f64 = 64.0;
const CENTER_Y: f64 = 32.0;
const R_I: f64 = 12.0;
const R_O: f64 = 28.0;
const R_MID: f64 = 20.0;
const NU: f64 = 1.0 / 60.0;
const TA_CRIT: f64 = 3390.0;
const SHELL_HALF_WIDTH: f64 = 0.55;
const DIVERGED_SPEED: f64 = 0.5;

#[derive(Clone, Debug)]
struct VelocityFields {
    ux: Vec<f64>,
    uy: Vec<f64>,
    uz: Vec<f64>,
}

#[derive(Clone, Debug)]
struct Spectrum {
    ta: f64,
    omega: f64,
    samples: usize,
    axial_energy: Vec<f64>,
    axisymmetric_energy: f64,
    max_speed_trajectory: Vec<(usize, f64)>,
}

fn idx(x: usize, y: usize, z: usize) -> usize {
    (z * NY + y) * NX + x
}

fn radius_xy(x: usize, y: usize) -> (f64, f64, f64) {
    let dx = x as f64 - CENTER_X;
    let dy = y as f64 - CENTER_Y;
    (dx, dy, (dx * dx + dy * dy).sqrt())
}

fn omega_for_ta(ta: f64) -> f64 {
    let gap = R_O - R_I;
    (ta * NU * NU / (R_I * gap * gap * gap)).sqrt()
}

fn circular_couette_coeffs(omega: f64) -> (f64, f64) {
    let denom = R_O * R_O - R_I * R_I;
    let a = -omega * R_I * R_I / denom;
    let b = omega * R_I * R_I * R_O * R_O / denom;
    (a, b)
}

fn shell_geometry(omega: f64) -> (Vec<bool>, Vec<[f64; 3]>) {
    let mut solid = vec![false; NX * NY * NZ];
    let mut wall_u = vec![[0.0; 3]; NX * NY * NZ];
    for z in 0..NZ {
        for y in 0..NY {
            for x in 0..NX {
                let (dx, dy, r) = radius_xy(x, y);
                let in_inner_shell = (r - R_I).abs() <= SHELL_HALF_WIDTH;
                let in_outer_shell = (r - R_O).abs() <= SHELL_HALF_WIDTH;
                if in_inner_shell || in_outer_shell {
                    let i = idx(x, y, z);
                    solid[i] = true;
                    if in_inner_shell {
                        wall_u[i] = [-omega * dy, omega * dx, 0.0];
                    }
                }
            }
        }
    }
    (solid, wall_u)
}

fn build_solver(omega: f64, perturb: bool) -> TaycSolver {
    let (solid, wall_u) = shell_geometry(omega);
    let spec = GlobalSpec::<f64> {
        dims: [NX, NY, NZ],
        nu: NU,
        periodic: [true, true, true],
        collision: CollisionKind::Trt {
            magic: CollisionKind::MAGIC_STD,
        },
        ..Default::default()
    };
    let mut sim: TaycSolver = Solver::new(
        &spec,
        &solid,
        &wall_u,
        [1, 1, 1],
        CpuSimd::default(),
        LocalPeriodic,
    );
    if perturb {
        seed_radial_impulse(&mut sim);
    }
    sim
}

fn seed_radial_impulse(sim: &mut TaycSolver) {
    // A deterministic infinitesimal radial impulse exposes the axial onset
    // mode. The subcritical case must damp it; supercritical cases should
    // amplify it. This is a one-step initial perturbation, not a persistent
    // closure term, and is far below the circular-Couette speed.
    sim.set_body_force_field(|x, y, z| {
        let (dx, dy, r) = radius_xy(x, y);
        if !(R_I + 2.0..R_O - 2.0).contains(&r) {
            return [0.0; 3];
        }
        let radial = 1.0e-6 * (TAU * z as f64 / NZ as f64).cos();
        [radial * dx / r, radial * dy / r, 0.0]
    });
    sim.step();
    sim.clear_body_force_field();
}

fn gather_velocity(sim: &TaycSolver) -> VelocityFields {
    VelocityFields {
        ux: sim.gather_ux(),
        uy: sim.gather_uy(),
        uz: sim.gather_uz(),
    }
}

fn max_speed(fields: &VelocityFields) -> f64 {
    fields
        .ux
        .iter()
        .zip(&fields.uy)
        .zip(&fields.uz)
        .map(|((ux, uy), uz)| (ux * ux + uy * uy + uz * uz).sqrt())
        .fold(0.0, f64::max)
}

fn assert_finite_or_stop(label: &str, step: usize, fields: &VelocityFields, traj: &[(usize, f64)]) {
    let max_u = max_speed(fields);
    let finite = fields
        .ux
        .iter()
        .chain(&fields.uy)
        .chain(&fields.uz)
        .all(|v| v.is_finite());
    if !finite || !max_u.is_finite() || max_u > DIVERGED_SPEED {
        panic!(
            "STOP-RULE: VAL-TAYC {label} thin-shell rotating-boundary run diverged at step {step}; \
             ANOM-P4-010 caveat still active. max|u|={max_u:.6e}, trajectory={traj:?}"
        );
    }
}

fn run_checked(
    sim: &mut TaycSolver,
    steps: usize,
    check_every: usize,
    label: &str,
) -> Vec<(usize, f64)> {
    let mut traj = Vec::new();
    let mut done = 0usize;
    while done < steps {
        let chunk = (steps - done).min(check_every);
        sim.run(chunk);
        done += chunk;
        let fields = gather_velocity(sim);
        let max_u = max_speed(&fields);
        traj.push((done, max_u));
        assert_finite_or_stop(label, done, &fields, &traj);
    }
    traj
}

fn laminar_profile_error(fields: &VelocityFields, omega: f64) -> (usize, f64, f64) {
    let (a, b) = circular_couette_coeffs(omega);
    let mut l2 = 0.0;
    let mut ref2 = 0.0;
    let mut linf = 0.0f64;
    let mut count = 0usize;
    for z in 0..NZ {
        for y in 0..NY {
            for x in 0..NX {
                let (dx, dy, r) = radius_xy(x, y);
                if !(R_I + 4.0..=R_O - 4.0).contains(&r) {
                    continue;
                }
                let i = idx(x, y, z);
                let utheta = (-fields.ux[i] * dy + fields.uy[i] * dx) / r;
                let reference = a * r + b / r;
                let d = utheta - reference;
                l2 += d * d;
                ref2 += reference * reference;
                linf = linf.max(d.abs());
                count += 1;
            }
        }
    }
    (count, (l2 / ref2).sqrt(), linf / (omega * R_I))
}

fn radial_line_signal(fields: &VelocityFields) -> (Vec<f64>, f64) {
    let x = (CENTER_X + R_MID).round() as usize;
    let y = CENTER_Y.round() as usize;
    let mut ur = Vec::with_capacity(NZ);
    let mut utheta_sum = 0.0;
    for z in 0..NZ {
        let i = idx(x, y, z);
        ur.push(fields.ux[i]);
        utheta_sum += fields.uy[i];
    }
    let mean_utheta = utheta_sum / NZ as f64;
    (ur, mean_utheta * mean_utheta)
}

fn axial_energy(signal: &[f64]) -> Vec<f64> {
    let n = signal.len();
    let mut out = vec![0.0; n / 2 + 1];
    for (k, slot) in out.iter_mut().enumerate() {
        let mut re = 0.0;
        let mut im = 0.0;
        for (j, &s) in signal.iter().enumerate() {
            let phase = TAU * k as f64 * j as f64 / n as f64;
            re += s * phase.cos();
            im -= s * phase.sin();
        }
        *slot = (re * re + im * im) / (n * n) as f64;
    }
    out
}

fn measure_spectrum(ta: f64, final_window: usize, sample_every: usize) -> Spectrum {
    let omega = omega_for_ta(ta);
    let mut sim = build_solver(omega, true);
    let settle_steps = 5000;
    let mut traj = run_checked(&mut sim, settle_steps, 500, &format!("Ta={ta:.3} settle"));
    let mut spectrum = vec![0.0; NZ / 2 + 1];
    let mut axisymmetric_energy = 0.0;
    let mut samples = 0usize;
    let mut elapsed = 0usize;
    while elapsed < final_window {
        let chunk = (final_window - elapsed).min(sample_every);
        sim.run(chunk);
        elapsed += chunk;
        if elapsed % sample_every == 0 || elapsed == final_window {
            let fields = gather_velocity(&sim);
            let global_step = settle_steps + elapsed;
            let max_u = max_speed(&fields);
            traj.push((global_step, max_u));
            assert_finite_or_stop(&format!("Ta={ta:.3} final"), global_step, &fields, &traj);
            let (signal, axis_e) = radial_line_signal(&fields);
            for (dst, src) in spectrum.iter_mut().zip(axial_energy(&signal)) {
                *dst += src;
            }
            axisymmetric_energy += axis_e;
            samples += 1;
        }
    }
    assert!(
        samples > 0,
        "spectrum measurement needs at least one sample"
    );
    for e in &mut spectrum {
        *e /= samples as f64;
    }
    axisymmetric_energy /= samples as f64;
    Spectrum {
        ta,
        omega,
        samples,
        axial_energy: spectrum,
        axisymmetric_energy,
        max_speed_trajectory: traj,
    }
}

fn nonzero_mode_ratio(spectrum: &Spectrum) -> f64 {
    let nonzero = spectrum.axial_energy[1..].iter().sum::<f64>();
    nonzero / spectrum.axisymmetric_energy.max(1.0e-30)
}

fn high_mode_ratio(spectrum: &Spectrum) -> f64 {
    let high = spectrum.axial_energy[2..].iter().sum::<f64>();
    high / spectrum.axisymmetric_energy.max(1.0e-30)
}

fn print_spectrum(label: &str, spectrum: &Spectrum) {
    let ratios: Vec<f64> = spectrum
        .axial_energy
        .iter()
        .map(|e| e / spectrum.axisymmetric_energy.max(1.0e-30))
        .collect();
    println!(
        "VAL-TAYC {label}: Ta={:.6e}, Omega={:.6e}, samples={}, axisymmetric_energy={:.6e}, \
         nonzero_ratio={:.6e}, high_mode_ratio={:.6e}, spectrum_abs={:?}, spectrum_ratio={:?}, max|u|_trajectory={:?}",
        spectrum.ta,
        spectrum.omega,
        spectrum.samples,
        spectrum.axisymmetric_energy,
        nonzero_mode_ratio(spectrum),
        high_mode_ratio(spectrum),
        spectrum.axial_energy,
        ratios,
        spectrum.max_speed_trajectory
    );
}

#[test]
fn val_tayc_laminar_circular_couette_profile_light() {
    let ta = 0.10 * TA_CRIT;
    let omega = omega_for_ta(ta);
    let mut sim = build_solver(omega, false);
    let traj = run_checked(&mut sim, 5000, 500, "laminar profile");
    let fields = gather_velocity(&sim);
    let (count, l2_rel, linf_rel_ui) = laminar_profile_error(&fields, omega);
    println!(
        "VAL-TAYC laminar: Ta={ta:.6e}, Omega={omega:.6e}, samples={count}, \
         L2_rel={l2_rel:.6e}, Linf/U_i={linf_rel_ui:.6e}, max|u|_trajectory={traj:?}"
    );
    assert!(count > 0, "bulk profile anchor had no samples");
    assert!(
        l2_rel <= 0.10,
        "laminar circular-Couette L2_rel={l2_rel:.6e} exceeds 10% band"
    );
}

#[test]
#[ignore = "heavy VAL-TAYC Axis 9.8: 128x64x64 D3Q19 x three Taylor numbers"]
fn val_tayc_wavy_vortex_onset_heavy() {
    let subcritical = measure_spectrum(0.5 * TA_CRIT, 2000, 20);
    print_spectrum("Ta=0.5Ta_c", &subcritical);
    let onset = measure_spectrum(1.5 * TA_CRIT, 2000, 20);
    print_spectrum("Ta=1.5Ta_c", &onset);
    let high = measure_spectrum(3.0 * TA_CRIT, 2000, 20);
    print_spectrum("Ta=3.0Ta_c", &high);

    let sub_ratio = nonzero_mode_ratio(&subcritical);
    let onset_ratio = nonzero_mode_ratio(&onset);
    let onset_high = high_mode_ratio(&onset);
    let high_high = high_mode_ratio(&high);

    assert!(
        sub_ratio < 0.01,
        "subcritical Ta=0.5Ta_c grew axial modes: nonzero/axisymmetric={sub_ratio:.6e}"
    );
    assert!(
        onset_ratio > 0.05,
        "supercritical Ta=1.5Ta_c did not show visible axial-mode growth: nonzero/axisymmetric={onset_ratio:.6e}"
    );
    assert!(
        high_high > onset_high,
        "Ta=3Ta_c did not move more energy into higher modes: high={high_high:.6e}, onset={onset_high:.6e}"
    );
}
