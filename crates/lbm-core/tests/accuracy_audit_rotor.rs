//! ACC ROT: adversarial accuracy audit for the compat penalized-rotor path.
//!
//! This file is gated like `tests/mf_interim.rs` because the compat rotor
//! belongs to the MF-interim surface. Run with:
//! `cargo test -p lbm-core --release --features mf-interim --test accuracy_audit_rotor -- --nocapture`
//!
//! Heavy/blocked probes are visible with:
//! `cargo test -p lbm-core --release --features mf-interim --test accuracy_audit_rotor -- --include-ignored --list`

#![cfg(feature = "mf-interim")]

mod common;

use common::metrics::*;
use lbm_core::compat::prelude::*;
use lbm_core::compat::rotor::Rotor;
use std::f64::consts::PI;

const NX: usize = 80;
const NY: usize = 80;
const CENTER: [f64; 2] = [40.0, 40.0];
const NU: f64 = 1.0 / 6.0;
const R_DISC: f64 = 10.0;
const R_OUTER_SOLID_CUT: f64 = 30.5;
const R_OUTER_EFF: f64 = 31.0;
const OMEGA_MID: f64 = 1.5e-4;
const STEADY_WINDOW: usize = 200;
const STEADY_REL: f64 = 1.0e-4;
const MAX_STEPS: usize = 6_000;

fn periodic_tank(center: [f64; 2]) -> Simulation<f64> {
    let mut sim = SimConfig {
        nx: NX,
        ny: NY,
        nu: NU,
        collision: Collision::Trt {
            magic: Collision::MAGIC_STD,
        },
        edges: Edges {
            left: EdgeBC::Periodic,
            right: EdgeBC::Periodic,
            bottom: EdgeBC::Periodic,
            top: EdgeBC::Periodic,
        },
        ..Default::default()
    }
    .build()
    .unwrap();
    sim.set_solid_region(|x, y| {
        let dx = x as f64 - center[0];
        let dy = y as f64 - center[1];
        dx.hypot(dy) > R_OUTER_SOLID_CUT
    });
    sim
}

fn disc_rotor(center: [f64; 2], omega: f64, ramp_steps: u64) -> Rotor<f64> {
    Rotor::new(center[0], center[1])
        .n_blades(1)
        .r_hub(R_DISC)
        .r_blade(R_DISC)
        .blade_thickness(2.0 * R_DISC)
        .omega(omega)
        .chi(1.0)
        .omega_ramp_steps(ramp_steps)
        .theta0(0.0)
}

fn annular_couette_torque_abs(omega: f64) -> f64 {
    annular_couette_torque_abs_with_radii(omega, R_DISC + 0.5, R_OUTER_EFF)
}

fn annular_couette_torque_abs_with_radii(omega: f64, r_i_eff: f64, r_o_eff: f64) -> f64 {
    let mu = NU; // rho ~= 1 in the Stokes, weakly-compressible limit.
    4.0 * PI * mu * omega * r_i_eff * r_i_eff * r_o_eff * r_o_eff
        / (r_o_eff * r_o_eff - r_i_eff * r_i_eff)
}

#[derive(Debug)]
struct SteadyTorque {
    steps: usize,
    mean: f64,
    previous_mean: f64,
    rel_change: f64,
}

fn run_disc_to_steady(center: [f64; 2], omega: f64) -> (Simulation<f64>, Rotor<f64>, SteadyTorque) {
    let mut sim = periodic_tank(center);
    let mut rotor = disc_rotor(center, omega, 0);
    let mut samples = Vec::with_capacity(MAX_STEPS);

    for step in 1..=MAX_STEPS {
        rotor.update_force(&mut sim);
        samples.push(rotor.torque());
        sim.step();
        if step >= 2 * STEADY_WINDOW && step % STEADY_WINDOW == 0 {
            let a = mean(&samples[step - STEADY_WINDOW..step]);
            let b = mean(&samples[step - 2 * STEADY_WINDOW..step - STEADY_WINDOW]);
            let rel_change = (a - b).abs() / a.abs().max(1.0e-30);
            if rel_change < STEADY_REL {
                let steady = SteadyTorque {
                    steps: step,
                    mean: a,
                    previous_mean: b,
                    rel_change,
                };
                return (sim, rotor, steady);
            }
        }
    }

    let n = samples.len();
    let a = mean(&samples[n - STEADY_WINDOW..n]);
    let b = mean(&samples[n - 2 * STEADY_WINDOW..n - STEADY_WINDOW]);
    let rel_change = (a - b).abs() / a.abs().max(1.0e-30);
    (
        sim,
        rotor,
        SteadyTorque {
            steps: n,
            mean: a,
            previous_mean: b,
            rel_change,
        },
    )
}

fn mean(xs: &[f64]) -> f64 {
    xs.iter().sum::<f64>() / xs.len() as f64
}

#[test]
fn f1_penalized_rotating_disc_torque_matches_annular_couette_canary() {
    // Annular Couette derivation, independent of implementation details:
    // Steady, axisymmetric, incompressible Stokes flow between concentric
    // cylinders has u_r = 0 and u_theta(r). The theta momentum equation reduces
    // to 0 = nu * (d2u/dr2 + (1/r) du/dr - u/r^2), whose general solution is
    // u_theta = a r + b/r. No slip gives u(r_i)=Omega*r_i and u(r_o)=0, hence
    // a = -Omega*r_i^2/(r_o^2-r_i^2) and
    // b = Omega*r_i^2*r_o^2/(r_o^2-r_i^2). The shear traction on the inner
    // cylinder is tau_rtheta = mu * (du/dr - u/r)|r_i = -2*mu*b/r_i^2. The
    // torque magnitude per unit depth is |T| = |tau| * (2*pi*r_i) * r_i =
    // 4*pi*mu*Omega*r_i^2*r_o^2/(r_o^2-r_i^2).
    //
    // Penalization acts on a staircase set of cell centers. For a disc of input
    // radius R=10 the effective no-slip radius carries a +-0.5 cell ambiguity;
    // this canary uses r_i_eff = R+0.5 because cell centers up to R are driven
    // and the no-slip surface is enforced at the rim. Around this geometry,
    // d ln(T)/d ln(r_i) = 2*r_o^2/(r_o^2-r_i^2) ~= 2.25, so a +-5% radius
    // ambiguity induces roughly +-11% torque ambiguity before staircase and
    // diffuse-interface errors are counted.
    let (_, _, steady) = run_disc_to_steady(CENTER, OMEGA_MID);
    let analytic = annular_couette_torque_abs(OMEGA_MID);
    let measured_abs = steady.mean.abs();
    let rel = (measured_abs - analytic).abs() / analytic;
    let ratio = steady.mean / analytic;
    println!(
        "ACC ROT F1: steps={} torque_mean={:.12e} prev_mean={:.12e} rel_change={:.3e} analytic_abs={:.12e} signed_ratio={:.6e} rel_err={:.6e}",
        steady.steps, steady.mean, steady.previous_mean, steady.rel_change, analytic, ratio, rel
    );
    assert!(
        steady.rel_change < STEADY_REL,
        "ACC ROT F1 steady rel_change={:.6e}, band={:.6e}, denominator=last-window mean",
        steady.rel_change,
        STEADY_REL
    );
    // Sign convention pinned by this measurement and the independent angular
    // momentum balance: rotor.torque() reports reaction torque on the rotor.
    // For positive omega, the rotor applies positive torque to the fluid and
    // the reported reaction torque must be negative.
    assert!(
        steady.mean < 0.0,
        "ACC ROT F1 sign: measured torque={:.12e}; expected negative reaction torque for omega={:.12e}",
        steady.mean,
        OMEGA_MID
    );
    assert!(
        rel <= 0.25,
        "ACC ROT F1 torque rel_err={:.6e}, band=2.500000e-1, denominator analytic={:.12e}, measured_abs={:.12e}",
        rel,
        analytic,
        measured_abs
    );
}

#[test]
#[ignore = "heavy ACC ROT F1 omega sweep {0.75e-4,1.5e-4,3e-4}"]
fn f1_penalized_rotating_disc_torque_annular_couette_heavy_sweep() {
    let omegas = [0.75e-4, 1.5e-4, 3.0e-4];
    let mut samples = Vec::new();
    let mut measured_abs = Vec::new();
    for omega in omegas {
        let (_, _, steady) = run_disc_to_steady(CENTER, omega);
        let t = steady.mean.abs();
        println!(
            "ACC ROT F1-heavy point: omega={:.12e} steps={} torque_mean={:.12e} rel_change={:.3e} analytic_abs={:.12e}",
            omega,
            steady.steps,
            steady.mean,
            steady.rel_change,
            annular_couette_torque_abs(omega)
        );
        samples.push((omega, t));
        measured_abs.push(t);
    }
    let fit = linear_fit(&omegas, &measured_abs);
    let agreement = curve_agreement(
        annular_couette_torque_abs,
        &samples,
        0.25,
        annular_couette_torque_abs(OMEGA_MID) * 1.0e-12,
    );
    let mid = annular_couette_torque_abs(OMEGA_MID);
    println!(
        "ACC ROT F1-heavy: slope={:.12e} intercept={:.12e} r2={:.12e} max_rel_dev={:.6e} worst_omega={:.12e} frac_in_band={:.3}",
        fit.slope,
        fit.intercept,
        fit.r2,
        agreement.max_rel_dev,
        agreement.worst_x,
        agreement.frac_in_band
    );
    assert!(
        fit.r2 >= 0.999,
        "ACC ROT F1-heavy linear_fit r2={:.12e}, band=9.990000e-1",
        fit.r2
    );
    assert!(
        fit.intercept.abs() <= 0.02 * mid,
        "ACC ROT F1-heavy intercept_abs={:.12e}, band={:.12e}, denominator T_mid={:.12e}",
        fit.intercept.abs(),
        0.02 * mid,
        mid
    );
    assert!(
        agreement.max_rel_dev <= 0.25,
        "ACC ROT F1-heavy curve max_rel_dev={:.6e}, band=2.500000e-1",
        agreement.max_rel_dev
    );
}

#[test]
fn f2_chi_one_disc_interior_tracks_rigid_body_after_f1_steady_state() {
    let (sim, rotor, steady) = run_disc_to_steady(CENTER, OMEGA_MID);
    assert!(
        steady.rel_change.is_finite() && steady.rel_change < STEADY_REL,
        "ACC ROT F2 prerequisite steady rel_change={:.6e}, band={:.6e}, denominator=last-window mean",
        steady.rel_change,
        STEADY_REL
    );
    let band = 1.0e-3 * OMEGA_MID * R_DISC;
    let mut max_dev = 0.0f64;
    let mut max_radius = 0.0f64;
    let mut max_point = [0usize; 2];
    let mut count = 0usize;
    for y in 0..sim.ny() {
        for x in 0..sim.nx() {
            if sim.is_solid(x, y) {
                continue;
            }
            let p = [x as f64, y as f64];
            let dx = p[0] - CENTER[0];
            let dy = p[1] - CENTER[1];
            let r = dx.hypot(dy);
            if r > R_DISC - 2.0 {
                continue;
            }
            let target = rotor.target_velocity(p);
            assert!(
                sim.ux(x, y).is_finite() && sim.uy(x, y).is_finite(),
                "ACC ROT F2 non-finite velocity at ({x},{y}): ux={:.12e} uy={:.12e}",
                sim.ux(x, y),
                sim.uy(x, y)
            );
            let dux = (sim.ux(x, y) - target[0]).abs();
            let duy = (sim.uy(x, y) - target[1]).abs();
            let dev = dux.max(duy);
            if dev > max_dev {
                max_dev = dev;
                max_radius = r;
                max_point = [x, y];
            }
            count += 1;
        }
    }
    println!(
        "ACC ROT F2: steps={} samples={} max_component_dev={:.12e} radius={:.6e} point=({},{}) band={:.12e}",
        steady.steps, count, max_dev, max_radius, max_point[0], max_point[1], band
    );
    assert!(count > 0, "ACC ROT F2 sampled no interior fluid cells");
    assert!(
        max_dev <= band,
        "ACC ROT F2 max_component_dev={:.12e}, band={:.12e}, denominator Omega*R={:.12e}",
        max_dev,
        band,
        OMEGA_MID * R_DISC
    );
}

#[test]
fn f3_disc_torque_subcell_translation_sensitivity_is_bounded() {
    let centers = [[40.0, 40.0], [40.3, 40.17], [40.5, 40.5]];
    let mut torques = Vec::new();
    for center in centers {
        let (_, _, steady) = run_disc_to_steady(center, OMEGA_MID);
        let t = steady.mean.abs();
        println!(
            "ACC ROT F3 point: center=({:.3},{:.3}) steps={} torque_mean={:.12e} abs={:.12e} rel_change={:.3e}",
            center[0], center[1], steady.steps, steady.mean, t, steady.rel_change
        );
        torques.push(t);
    }
    let min_t = torques.iter().copied().fold(f64::INFINITY, f64::min);
    let max_t = torques.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let mean_t = mean(&torques);
    let spread = (max_t - min_t) / mean_t;
    println!(
        "ACC ROT F3: torques={:?} spread={:.6e} band=1.500000e-1",
        torques, spread
    );
    assert!(
        spread <= 0.15,
        "ACC ROT F3 spread={:.6e}, band=1.500000e-1, denominator mean_torque={:.12e}",
        spread,
        mean_t
    );
}

#[test]
#[ignore = "blocked on ANOM-P4-001: cross-path torque referee runs after the core IBM fix"]
fn f4_penalized_disc_and_rotating_ibm_torques_referee_same_annular_couette_geometry() {
    use lbm_core::prelude::{
        CollisionKind, CpuScalar, DirectForcingConfig, GlobalSpec, IbmDiagnostics, InProcess,
        RotatingBody, Solver, D2Q9,
    };

    let (_, _, penalized) = run_disc_to_steady(CENTER, OMEGA_MID);

    let spec = GlobalSpec {
        dims: [NX, NY, 1],
        nu: NU,
        periodic: [true, true, false],
        collision: CollisionKind::Trt {
            magic: CollisionKind::MAGIC_STD,
        },
        ..Default::default()
    };
    let mut solid = vec![false; NX * NY];
    for y in 0..NY {
        for x in 0..NX {
            let dx = x as f64 - CENTER[0];
            let dy = y as f64 - CENTER[1];
            solid[y * NX + x] = dx.hypot(dy) > R_OUTER_SOLID_CUT;
        }
    }
    let mut ibm: Solver<D2Q9, f64, CpuScalar, InProcess> = Solver::new(
        &spec,
        &solid,
        &vec![[0.0; 3]; NX * NY],
        [1, 1, 1],
        CpuScalar::default(),
        InProcess,
    );
    let body = RotatingBody::circle_2d(CENTER, R_DISC + 0.5, OMEGA_MID, 160);
    let cfg = DirectForcingConfig {
        max_iterations: 1,
        slip_tolerance: 1.0,
        kernel_radius: 1,
        relaxation: 0.05,
    };
    let mut last = IbmDiagnostics::default();
    let mut ibm_torques = Vec::with_capacity(MAX_STEPS);
    for step in 1..=MAX_STEPS {
        ibm.clear_body_force_field();
        last = ibm.apply_rotating_ibm(&body, cfg);
        ibm_torques.push(last.torque[2]);
        ibm.step();
        if step >= 2 * STEADY_WINDOW && step % STEADY_WINDOW == 0 {
            let a = mean(&ibm_torques[step - STEADY_WINDOW..step]);
            let b = mean(&ibm_torques[step - 2 * STEADY_WINDOW..step - STEADY_WINDOW]);
            if (a - b).abs() / a.abs().max(1.0e-30) < STEADY_REL {
                break;
            }
        }
    }
    let ibm_mean = mean(&ibm_torques[ibm_torques.len() - STEADY_WINDOW..]);
    let analytic = annular_couette_torque_abs(OMEGA_MID);
    let penalized_abs = penalized.mean.abs();
    let ibm_abs = ibm_mean.abs();
    let cross_rel = (penalized_abs - ibm_abs).abs() / ((penalized_abs + ibm_abs) * 0.5);
    let pen_rel = (penalized_abs - analytic).abs() / analytic;
    let ibm_rel = (ibm_abs - analytic).abs() / analytic;
    println!(
        "ACC ROT F4: penalized={:.12e} ibm={:.12e} analytic={:.12e} cross_rel={:.6e} pen_rel={:.6e} ibm_rel={:.6e} ibm_slip_max_rel={:.6e}",
        penalized.mean, ibm_mean, analytic, cross_rel, pen_rel, ibm_rel, last.slip_max_rel
    );
    assert!(
        cross_rel <= 0.15,
        "ACC ROT F4 cross_rel={:.6e}, band=1.500000e-1, denominator mean_abs_pair",
        cross_rel
    );
    assert!(
        pen_rel <= 0.25,
        "ACC ROT F4 penalized analytic rel={:.6e}, band=2.500000e-1, denominator analytic={:.12e}",
        pen_rel,
        analytic
    );
    assert!(
        ibm_rel <= 0.25,
        "ACC ROT F4 IBM analytic rel={:.6e}, band=2.500000e-1, denominator analytic={:.12e}",
        ibm_rel,
        analytic
    );
}

#[test]
fn f5_ramp_torque_integral_matches_discrete_per_step_sum() {
    let mut sim = periodic_tank(CENTER);
    let mut rotor = Rotor::new(CENTER[0], CENTER[1])
        .n_blades(4)
        .r_hub(4.0)
        .r_blade(18.0)
        .blade_thickness(1.5)
        .omega(0.002)
        .chi(0.7)
        .omega_ramp_steps(200)
        .theta0(0.0);
    let mut sum = 0.0;
    for _ in 0..400 {
        rotor.update_force(&mut sim);
        sum += rotor.torque();
        sim.step();
    }
    let integral = rotor.torque_integral();
    let denom = integral.abs();
    let rel = (integral - sum).abs() / denom;
    println!(
        "ACC ROT F5: torque_integral={:.12e} accumulated_sum={:.12e} abs_diff={:.12e} rel={:.6e}",
        integral,
        sum,
        (integral - sum).abs(),
        rel
    );
    assert!(
        rel <= 1.0e-12,
        "ACC ROT F5 rel={:.6e}, band=1.000000e-12, denominator |torque_integral|={:.12e}, integral={:.12e}, sum={:.12e}",
        rel,
        denom,
        integral,
        sum
    );
}
