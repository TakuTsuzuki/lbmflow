//! VAL-TAYC-IBM / V&V Axis 9.8 rev 3: Taylor-Couette onset with a
//! marker-based rotating inner cylinder.
//!
//! ANOM-P4-024 established that the rev-2 thin cylindrical shells driven by
//! volume penalization filter the axial wavy-vortex mode. This revision keeps
//! the stationary outer cylinder as an Eulerian solid mask, but applies the
//! rotating inner cylinder through the direct-forcing IBM path validated by
//! ANOM-P4-001.
//!
//! Physics anchor: inner-cylinder angular speed is selected from the
//! narrow-gap Rayleigh Taylor number
//! `Ta = Omega_i^2 R_i (R_o - R_i)^3 / nu^2`. The textbook critical value
//! `Ta_c ~= 3390` is used as an onset discriminator, not as an exact
//! growth-rate benchmark for this finite-aspect-ratio lattice. Linear
//! Taylor-Couette theory gives the critical axial wavelength
//! `lambda_c ~= 2 (R_o - R_i) = 32`; with `L_z = 64`, this is Fourier mode
//! `m = 2`, so the imposed perturbation is
//! `u_r' = 1e-4 cos(4 pi z / L_z)`.

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
const R_O_EFF: f64 = R_O + 1.0;
const R_MID: f64 = 20.0;
const NU: f64 = 1.0 / 60.0;
const TA_CRIT: f64 = 3390.0;
const AXIAL_SEED_EPS: f64 = 1.0e-4;
const LIGHT_PROFILE_STEPS: usize = 200;
const WARMUP_STEPS: usize = 3000;
const MEASURE_STEPS: usize = 2000;
const SAMPLE_EVERY: usize = 20;
const DIVERGED_SPEED: f64 = 0.5;
const MARKERS_PER_PLANE: usize = 96;

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
    ratio_trajectory: Vec<RatioSample>,
    max_speed_trajectory: Vec<(usize, f64)>,
    ibm_slip_trajectory: Vec<(usize, f64, f64, f64)>,
}

#[derive(Clone, Copy, Debug)]
struct RatioSample {
    step: usize,
    nonzero_ratio: f64,
    high_mode_ratio: f64,
    axisymmetric_energy: f64,
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

fn ibm_cfg() -> DirectForcingConfig {
    DirectForcingConfig {
        max_iterations: 3,
        slip_tolerance: 1.0e-3,
        kernel_radius: 1,
        relaxation: 1.0,
    }
}

fn circular_couette_coeffs(omega: f64) -> (f64, f64) {
    let denom = R_O_EFF * R_O_EFF - R_I * R_I;
    let a = -omega * R_I * R_I / denom;
    let b = omega * R_I * R_I * R_O_EFF * R_O_EFF / denom;
    (a, b)
}

fn outer_cylinder_solid() -> Vec<bool> {
    let mut solid = vec![false; NX * NY * NZ];
    for z in 0..NZ {
        for y in 0..NY {
            for x in 0..NX {
                let (_, _, r) = radius_xy(x, y);
                solid[idx(x, y, z)] = r > R_O + 0.5;
            }
        }
    }
    solid
}

fn extruded_inner_cylinder(omega: f64) -> RotatingBody {
    let ds = TAU * R_I / MARKERS_PER_PLANE as f64;
    let mut markers = Vec::with_capacity(NZ * MARKERS_PER_PLANE);
    for z in 0..NZ {
        for i in 0..MARKERS_PER_PLANE {
            let th = TAU * i as f64 / MARKERS_PER_PLANE as f64;
            markers.push(IbmMarker {
                position: [
                    CENTER_X + R_I * th.cos(),
                    CENTER_Y + R_I * th.sin(),
                    z as f64,
                ],
                weight: ds,
            });
        }
    }
    RotatingBody::from_markers([CENTER_X, CENTER_Y, 0.0], [0.0, 0.0, omega], markers)
}

fn build_solver() -> TaycSolver {
    let spec = GlobalSpec::<f64> {
        dims: [NX, NY, NZ],
        nu: NU,
        periodic: [true, true, true],
        collision: CollisionKind::Trt {
            magic: CollisionKind::MAGIC_STD,
        },
        ..Default::default()
    };
    Solver::new(
        &spec,
        &outer_cylinder_solid(),
        &vec![[0.0; 3]; NX * NY * NZ],
        [1, 1, 1],
        CpuSimd::default(),
        LocalPeriodic,
    )
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
            "STOP-RULE: VAL-TAYC-IBM {label} rotating-IBM run diverged at step {step}; \
             max|u|={max_u:.6e}, trajectory={traj:?}"
        );
    }
}

fn ibm_step(sim: &mut TaycSolver, body: &RotatingBody, cfg: DirectForcingConfig) -> IbmDiagnostics {
    sim.clear_body_force_field();
    let diag = sim.apply_rotating_ibm(body, cfg);
    sim.step();
    diag
}

fn run_checked(
    sim: &mut TaycSolver,
    body: &RotatingBody,
    steps: usize,
    check_every: usize,
    label: &str,
) -> (Vec<(usize, f64)>, Vec<(usize, f64, f64, f64)>) {
    let mut max_u_traj = Vec::new();
    let mut slip_traj = Vec::new();
    let cfg = ibm_cfg();
    for step in 1..=steps {
        let diag = ibm_step(sim, body, cfg);
        if step % check_every == 0 || step == steps {
            let fields = gather_velocity(sim);
            let max_u = max_speed(&fields);
            max_u_traj.push((step, max_u));
            slip_traj.push((step, diag.slip_max_rel, diag.slip_rms_rel, diag.torque[2]));
            assert_finite_or_stop(label, step, &fields, &max_u_traj);
        }
    }
    (max_u_traj, slip_traj)
}

fn axial_seed_delta_ur(x: usize, y: usize, z: usize) -> f64 {
    let (_, _, r) = radius_xy(x, y);
    if !(R_I + 2.0..R_O - 2.0).contains(&r) {
        return 0.0;
    }
    let sigma = (R_O - R_I) / 4.0;
    let radial_envelope = (-((r - R_MID) * (r - R_MID)) / (sigma * sigma)).exp();
    let kz = 4.0 * std::f64::consts::PI / NZ as f64;
    AXIAL_SEED_EPS * radial_envelope * (kz * z as f64).cos()
}

fn inject_axial_seed(sim: &mut TaycSolver) {
    let rho = sim.gather_rho();
    let ux = sim.gather_ux();
    let uy = sim.gather_uy();
    let uz = sim.gather_uz();
    sim.init_with(|x, y, z| {
        let i = idx(x, y, z);
        let (dx, dy, r) = radius_xy(x, y);
        let delta_ur = axial_seed_delta_ur(x, y, z);
        let radial = if r > 0.0 {
            [delta_ur * dx / r, delta_ur * dy / r, 0.0]
        } else {
            [0.0; 3]
        };
        (
            rho[i],
            [ux[i] + radial[0], uy[i] + radial[1], uz[i] + radial[2]],
        )
    });
}

fn init_laminar_couette(sim: &mut TaycSolver, omega: f64) {
    let (a, b) = circular_couette_coeffs(omega);
    sim.init_with(|x, y, _z| {
        let (dx, dy, r) = radius_xy(x, y);
        let utheta = if r < 1.0e-12 {
            0.0
        } else if r < R_I {
            omega * r
        } else if r <= R_O_EFF {
            a * r + b / r
        } else {
            0.0
        };
        let u = if r > 0.0 {
            [-utheta * dy / r, utheta * dx / r, 0.0]
        } else {
            [0.0; 3]
        };
        (1.0, u)
    });
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
                if !(R_I + 4.0..=R_O - 5.0).contains(&r) {
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

fn ratio_sample(step: usize, spectrum: &[f64], axisymmetric_energy: f64) -> RatioSample {
    RatioSample {
        step,
        nonzero_ratio: spectrum[1..].iter().sum::<f64>() / axisymmetric_energy.max(1.0e-30),
        high_mode_ratio: spectrum[2..].iter().sum::<f64>() / axisymmetric_energy.max(1.0e-30),
        axisymmetric_energy,
    }
}

fn sample_spectrum(sim: &TaycSolver, step: usize) -> (Vec<f64>, f64, RatioSample, VelocityFields) {
    let fields = gather_velocity(sim);
    let (signal, axis_e) = radial_line_signal(&fields);
    let spectrum = axial_energy(&signal);
    let ratio = ratio_sample(step, &spectrum, axis_e);
    (spectrum, axis_e, ratio, fields)
}

fn measure_spectrum(ta: f64) -> Spectrum {
    let omega = omega_for_ta(ta);
    let body = extruded_inner_cylinder(omega);
    let mut sim = build_solver();
    let (mut max_u_traj, mut slip_traj) = run_checked(
        &mut sim,
        &body,
        WARMUP_STEPS,
        500,
        &format!("Ta={ta:.3} warmup"),
    );
    inject_axial_seed(&mut sim);

    let mut spectrum = vec![0.0; NZ / 2 + 1];
    let mut axisymmetric_energy = 0.0;
    let mut ratio_trajectory = Vec::new();
    let mut samples = 0usize;

    let (seed_spectrum, seed_axis_e, seed_ratio, seed_fields) = sample_spectrum(&sim, WARMUP_STEPS);
    max_u_traj.push((WARMUP_STEPS, max_speed(&seed_fields)));
    assert_finite_or_stop(
        &format!("Ta={ta:.3} seeded"),
        WARMUP_STEPS,
        &seed_fields,
        &max_u_traj,
    );
    for (dst, src) in spectrum.iter_mut().zip(seed_spectrum) {
        *dst += src;
    }
    axisymmetric_energy += seed_axis_e;
    ratio_trajectory.push(seed_ratio);
    samples += 1;

    let cfg = ibm_cfg();
    for elapsed in 1..=MEASURE_STEPS {
        let diag = ibm_step(&mut sim, &body, cfg);
        if elapsed % SAMPLE_EVERY == 0 || elapsed == MEASURE_STEPS {
            let global_step = WARMUP_STEPS + elapsed;
            let (sample_spectrum, sample_axis_e, ratio, fields) =
                sample_spectrum(&sim, global_step);
            max_u_traj.push((global_step, max_speed(&fields)));
            slip_traj.push((
                global_step,
                diag.slip_max_rel,
                diag.slip_rms_rel,
                diag.torque[2],
            ));
            assert_finite_or_stop(
                &format!("Ta={ta:.3} measure"),
                global_step,
                &fields,
                &max_u_traj,
            );
            for (dst, src) in spectrum.iter_mut().zip(sample_spectrum) {
                *dst += src;
            }
            axisymmetric_energy += sample_axis_e;
            ratio_trajectory.push(ratio);
            samples += 1;
        }
    }

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
        ratio_trajectory,
        max_speed_trajectory: max_u_traj,
        ibm_slip_trajectory: slip_traj,
    }
}

fn nonzero_mode_ratio(spectrum: &Spectrum) -> f64 {
    spectrum.axial_energy[1..].iter().sum::<f64>() / spectrum.axisymmetric_energy.max(1.0e-30)
}

fn high_mode_ratio(spectrum: &Spectrum) -> f64 {
    spectrum.axial_energy[2..].iter().sum::<f64>() / spectrum.axisymmetric_energy.max(1.0e-30)
}

fn ratio_trajectory_values(spectrum: &Spectrum) -> Vec<(usize, f64, f64, f64)> {
    spectrum
        .ratio_trajectory
        .iter()
        .map(|s| {
            (
                s.step,
                s.nonzero_ratio,
                s.high_mode_ratio,
                s.axisymmetric_energy,
            )
        })
        .collect()
}

fn print_spectrum(label: &str, spectrum: &Spectrum) {
    let ratios: Vec<f64> = spectrum
        .axial_energy
        .iter()
        .map(|e| e / spectrum.axisymmetric_energy.max(1.0e-30))
        .collect();
    println!(
        "VAL-TAYC-IBM {label}: Ta={:.6e}, Omega={:.6e}, samples={}, axisymmetric_energy={:.6e}, \
         nonzero_ratio={:.6e}, high_mode_ratio={:.6e}, ratio_trajectory={:?}, spectrum_abs={:?}, \
         spectrum_ratio={:?}, max|u|_trajectory={:?}, ibm_slip_trajectory={:?}",
        spectrum.ta,
        spectrum.omega,
        spectrum.samples,
        spectrum.axisymmetric_energy,
        nonzero_mode_ratio(spectrum),
        high_mode_ratio(spectrum),
        ratio_trajectory_values(spectrum),
        spectrum.axial_energy,
        ratios,
        spectrum.max_speed_trajectory,
        spectrum.ibm_slip_trajectory
    );
}

#[test]
fn val_tayc_ibm_laminar_profile_light() {
    let ta = 0.5 * TA_CRIT;
    let omega = omega_for_ta(ta);
    let body = extruded_inner_cylinder(omega);
    let mut sim = build_solver();
    init_laminar_couette(&mut sim, omega);
    let (max_u_traj, slip_traj) = run_checked(
        &mut sim,
        &body,
        LIGHT_PROFILE_STEPS,
        100,
        "light laminar profile",
    );
    let fields = gather_velocity(&sim);
    let (count, l2_rel, linf_rel_ui) = laminar_profile_error(&fields, omega);
    println!(
        "VAL-TAYC-IBM laminar: Ta={ta:.6e}, Omega={omega:.6e}, samples={count}, \
         L2_rel={l2_rel:.6e}, Linf/U_i={linf_rel_ui:.6e}, max|u|_trajectory={max_u_traj:?}, \
         ibm_slip_trajectory={slip_traj:?}"
    );
    assert!(count > 0, "bulk profile anchor had no samples");
    assert!(
        l2_rel <= 0.15,
        "VAL-TAYC-IBM laminar circular-Couette L2_rel={l2_rel:.6e} exceeds 15% band"
    );
    assert!(
        linf_rel_ui <= 0.35,
        "VAL-TAYC-IBM laminar circular-Couette Linf/U_i={linf_rel_ui:.6e} exceeds behavior band 35%"
    );
}

#[test]
#[ignore = "heavy VAL-TAYC-IBM Axis 9.8 / ANOM-P4-024 route: 128x64x64 D3Q19 x three Taylor numbers"]
fn val_tayc_ibm_wavy_vortex_onset_heavy() {
    let subcritical = measure_spectrum(0.5 * TA_CRIT);
    print_spectrum("Ta=0.5Ta_c", &subcritical);
    let onset = measure_spectrum(1.5 * TA_CRIT);
    print_spectrum("Ta=1.5Ta_c", &onset);
    let high = measure_spectrum(3.0 * TA_CRIT);
    print_spectrum("Ta=3.0Ta_c", &high);

    let sub_start = subcritical.ratio_trajectory.first().unwrap();
    let sub_end = subcritical.ratio_trajectory.last().unwrap();
    let onset_start = onset.ratio_trajectory.first().unwrap();
    let onset_end = onset.ratio_trajectory.last().unwrap();
    let high_end = high.ratio_trajectory.last().unwrap();

    println!(
        "VAL-TAYC-IBM damping/growth: seed_mode=2, eps={AXIAL_SEED_EPS:.6e}, \
         sub_start={:.6e}, sub_end={:.6e}, onset_start={:.6e}, onset_end={:.6e}, \
         high_end={:.6e}, onset_high_end={:.6e}, high_high_end={:.6e}",
        sub_start.nonzero_ratio,
        sub_end.nonzero_ratio,
        onset_start.nonzero_ratio,
        onset_end.nonzero_ratio,
        high_end.nonzero_ratio,
        onset_end.high_mode_ratio,
        high_end.high_mode_ratio
    );

    assert!(
        sub_end.nonzero_ratio < sub_start.nonzero_ratio,
        "VAL-TAYC-IBM subcritical Ta=0.5Ta_c did not damp the axial seed: start={:.6e}, end={:.6e}, trajectory={:?}",
        sub_start.nonzero_ratio,
        sub_end.nonzero_ratio,
        ratio_trajectory_values(&subcritical)
    );
    assert!(
        onset_end.nonzero_ratio > onset_start.nonzero_ratio && onset_end.nonzero_ratio > 0.05,
        "VAL-TAYC-IBM residual finding: Ta=1.5Ta_c rotating IBM still filters/damps the wavy-vortex seed. \
         start={:.6e}, end={:.6e}, trajectory={:?}",
        onset_start.nonzero_ratio,
        onset_end.nonzero_ratio,
        ratio_trajectory_values(&onset)
    );
    assert!(
        high_end.high_mode_ratio > onset_end.high_mode_ratio,
        "VAL-TAYC-IBM Ta=3Ta_c did not move more energy into higher modes than Ta=1.5Ta_c: \
         high={:.6e}, onset={:.6e}, high_trajectory={:?}, onset_trajectory={:?}",
        high_end.high_mode_ratio,
        onset_end.high_mode_ratio,
        ratio_trajectory_values(&high),
        ratio_trajectory_values(&onset)
    );
}
