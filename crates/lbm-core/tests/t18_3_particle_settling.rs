//! T18.3 runnable acceptance tests for the existing one-way particle layer.
//!
//! These tests deliberately cover only the public `ParticleSet::step` API that
//! exists before CR-3. Deposition-specific checks live in
//! `t18_3_particle_deposition.rs` and are expected not to compile until the
//! frozen CR-3 API lands.

use approx::relative_eq;
use lbm_core::particles::{Particle, ParticleSet, Sample};

const RHO_F: f64 = 1.0;
const RHO_P: f64 = 2.0;
const NU: f64 = 0.1;
const TERMINAL_TOL: f64 = 0.02;
const REFERENCE_CONV_TOL: f64 = 1.0e-12;

#[derive(Clone, Copy)]
struct SettlingCase {
    name: &'static str,
    d: f64,
    g: f64,
    require_stokes_limit: bool,
}

fn fluid(u: [f64; 3]) -> impl Fn([f64; 3]) -> Sample {
    move |_| Sample { u, solid: false }
}

fn particle(d: f64) -> Particle {
    Particle {
        pos: [7.0, 11.0, 1_000.0],
        vel: [0.0; 3],
        d,
        rho_p: RHO_P,
        exposure: 0.0,
    }
}

fn stokes_terminal_velocity(d: f64, g: f64) -> f64 {
    (RHO_P / RHO_F - 1.0) * g * d * d / (18.0 * NU)
}

/// Independent Schiller-Naumann terminal settling solve.
///
/// At terminal velocity in still fluid,
/// `v = tau_p(Re(v)) * g * (1 - rho_f/rho_p)`, with
/// `tau_p = rho_p d^2 / (18 rho_f nu (1 + 0.15 Re^0.687))`.
/// The fixed-point iteration below starts from the Stokes velocity and is run
/// until the absolute update is <= 1e-12 lattice units; the test asserts that
/// convergence happened instead of silently accepting a loose reference.
fn schiller_naumann_terminal_velocity(d: f64, g: f64) -> (f64, f64, usize, f64) {
    let g_eff = g * (1.0 - RHO_F / RHO_P);
    let tau_stokes = RHO_P * d * d / (18.0 * RHO_F * NU);
    let mut v = tau_stokes * g_eff;
    let mut last_delta = f64::INFINITY;

    for iter in 1..=10_000 {
        let re = v.abs() * d / NU;
        let drag_correction = 1.0 + 0.15 * re.min(800.0).powf(0.687);
        let next = tau_stokes * g_eff / drag_correction;
        last_delta = (next - v).abs();
        v = next;
        if last_delta <= REFERENCE_CONV_TOL {
            return (v, v.abs() * d / NU, iter, last_delta);
        }
    }

    panic!(
        "SN fixed-point reference did not converge for d={d:e}, g={g:e}; last_delta={last_delta:e}"
    );
}

fn run_until_terminal(d: f64, g: f64) -> f64 {
    let mut set = ParticleSet::new(vec![particle(d)], RHO_F, NU, [0.0, 0.0, -g]);
    for _ in 0..20_000 {
        set.step(fluid([0.0; 3]), None::<fn([f64; 3]) -> f64>);
    }
    -set.particles[0].vel[2]
}

#[test]
fn t18_3_single_particle_terminal_velocity_matches_schiller_naumann() {
    let cases = [
        // Stokes response time tau_p spans 1e-4..1e1, i.e. five decades.
        SettlingCase {
            name: "low-Re tau=1e-4",
            d: 9.486_832_980_505_138e-3,
            g: 1.0e-6,
            require_stokes_limit: true,
        },
        SettlingCase {
            name: "low-Re tau=1e-3",
            d: 3.0e-2,
            g: 1.0e-6,
            require_stokes_limit: true,
        },
        SettlingCase {
            name: "low-Re tau=1e-1",
            d: 3.0e-1,
            g: 1.0e-6,
            require_stokes_limit: true,
        },
        SettlingCase {
            name: "low-Re tau=1e1",
            d: 3.0,
            g: 1.0e-6,
            require_stokes_limit: true,
        },
        SettlingCase {
            name: "moderate-Re SN",
            d: 4.0,
            g: 3.0e-2,
            require_stokes_limit: false,
        },
    ];

    for case in cases {
        let (want, re, iters, residual) = schiller_naumann_terminal_velocity(case.d, case.g);
        let got = run_until_terminal(case.d, case.g);
        let rel = (got - want).abs() / want.abs();
        assert!(
            relative_eq!(got, want, max_relative = TERMINAL_TOL),
            "{}: terminal velocity rel_err={rel:e}, got={got:e}, want={want:e}, Re_p={re:e}, \
             reference_iters={iters}, reference_residual={residual:e}",
            case.name
        );

        if case.require_stokes_limit {
            let stokes = stokes_terminal_velocity(case.d, case.g);
            let stokes_rel = (want - stokes).abs() / stokes.abs();
            assert!(
                re < 0.1 && stokes_rel < 1.0e-3,
                "{}: low-Re Stokes check failed; Re_p={re:e}, SN_vs_Stokes_rel={stokes_rel:e}, \
                 SN={want:e}, Stokes={stokes:e}",
                case.name
            );
        } else {
            assert!(
                (1.0..=10.0).contains(&re),
                "{}: expected moderate Re_p in [1, 10], got Re_p={re:e}, v={want:e}",
                case.name
            );
        }

        println!(
            "{}: got={got:.12e}, SN={want:.12e}, rel_err={rel:.6e}, Re_p={re:.6e}, \
             reference_iters={iters}",
            case.name
        );
    }
}

#[test]
fn t18_3_particle_step_is_bit_deterministic() {
    fn build() -> ParticleSet {
        ParticleSet::new(
            vec![
                Particle {
                    pos: [0.125, 4.0, 9.0],
                    vel: [0.02, -0.01, 0.0],
                    d: 0.07,
                    rho_p: 1.5,
                    exposure: 0.0,
                },
                Particle {
                    pos: [2.5, 1.25, 3.75],
                    vel: [-0.03, 0.015, -0.004],
                    d: 0.4,
                    rho_p: 2.1,
                    exposure: 0.0,
                },
            ],
            1.0,
            0.08,
            [0.0, -2.0e-5, -1.0e-5],
        )
        .with_restitution(0.25)
    }

    let mut a = build();
    let mut b = build();
    let sample = |p: [f64; 3]| Sample {
        u: [
            0.01 + 1.0e-4 * p[1],
            -0.02 + 2.0e-4 * p[2],
            0.005 - 1.0e-4 * p[0],
        ],
        solid: p[1] < -0.25,
    };
    let exposure = |p: [f64; 3]| 0.125 + 0.01 * p[0] - 0.02 * p[2];

    for step in 0..256 {
        a.step(sample, Some(exposure));
        b.step(sample, Some(exposure));
        assert_eq!(
            a.particles, b.particles,
            "bit-identical replay diverged at step {step}: left={:?}, right={:?}",
            a.particles, b.particles
        );
    }
}
