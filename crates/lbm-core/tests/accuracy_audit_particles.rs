//! CR-3 adversarial accuracy-audit probes for the Lagrangian particle layer.
//!
//! These tests intentionally use only the public particle API plus analytic
//! references. They are not regression tests for implementation details.

mod common;

use common::metrics::*;
use lbm_core::particles::{sample_grid, DepositEvent, Particle, ParticleSet, Sample};

const RHO_F: f64 = 1.0;
const RHO_P: f64 = 2.0;
const NU: f64 = 0.1;
// Triage 2026-07-06 (ANOM-P4-004): v0 = 1e-4 violated the Stokes-regime
// assumption because hitting tau_p = 2 at rho_p = 2, nu = 0.1 needs
// d = sqrt(18*0.1*2/2) ~= 1.34 lattice units, so Re_p = v0*d/nu = 1.34e-3
// and the SN factor 1 + 0.15*Re^0.687 shifts lambda by 1.3e-3 relative —
// exactly the first-pass failure. With v0 = 1e-10, Re_p = 1.34e-9 and the
// SN residual on lambda is 0.15*Re^0.687 = 1.2e-7, an ~8x margin under the
// 1e-6 identification band. The integrator identity itself is v0-invariant
// (linear ODE), so shrinking v0 only purifies the regime.
const V0: f64 = 1.0e-10;
const REL_BAND: f64 = 1.0e-6;
const ABS_BAND: f64 = 1.0e-12;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Scheme {
    ExactExponential,
    ExplicitEuler,
    BackwardEuler,
}

impl Scheme {
    fn name(self) -> &'static str {
        match self {
            Scheme::ExactExponential => "exact exponential",
            Scheme::ExplicitEuler => "explicit Euler",
            Scheme::BackwardEuler => "backward Euler",
        }
    }
}

fn still_fluid(_: [f64; 3]) -> Sample {
    Sample {
        u: [0.0; 3],
        solid: false,
    }
}

fn particle_for_tau(tau_p: f64, vel: [f64; 3], rho_p: f64) -> Particle {
    // Stokes relaxation time is tau_p = rho_p d^2 / (18 rho_f nu), so the
    // requested tau_p fixes d = sqrt(tau_p 18 rho_f nu / rho_p).
    let d = (tau_p * 18.0 * RHO_F * NU / rho_p).sqrt();
    Particle {
        pos: [3.0, 4.0, 5.0],
        vel,
        d,
        rho_p,
        exposure: 0.0,
    }
}

fn one_step_lambda(tau_p: f64) -> f64 {
    let mut set = ParticleSet::new(
        vec![particle_for_tau(tau_p, [V0, 0.0, 0.0], RHO_P)],
        RHO_F,
        NU,
        [0.0; 3],
    );
    set.step(still_fluid, None::<fn([f64; 3]) -> f64>);
    -(set.particles[0].vel[0] / V0).ln()
}

fn lambda_curve(scheme: Scheme, tau_p: f64) -> f64 {
    match scheme {
        // The exact ODE in still fluid and the Stokes limit is dv/dt = -v/tau_p.
        // Integrating over dt=1 gives v_{n+1} = v_n exp(-1/tau_p), therefore
        // lambda = -ln(v_{n+1}/v_n) = 1/tau_p.
        Scheme::ExactExponential => 1.0 / tau_p,
        // Forward Euler gives v_{n+1} = v_n + dt(-v_n/tau_p)
        // = v_n(1 - 1/tau_p), therefore lambda = -ln(1 - 1/tau_p).
        // For large tau_p this is 1/tau_p + 1/(2 tau_p^2) + O(tau_p^-3).
        Scheme::ExplicitEuler => -(1.0 - 1.0 / tau_p).ln(),
        // Backward Euler evaluates drag at the new velocity:
        // v_{n+1} = v_n - v_{n+1}/tau_p, so
        // v_{n+1}/v_n = 1/(1 + 1/tau_p). Hence
        // lambda = ln(1 + 1/tau_p)
        // = 1/tau_p - 1/(2 tau_p^2) + O(tau_p^-3).
        Scheme::BackwardEuler => (1.0 + 1.0 / tau_p).ln(),
    }
}

fn identify_integrator_signature() -> (Scheme, [(f64, f64); 4], [(Scheme, f64); 3]) {
    let samples = [2.0, 5.0, 10.0, 40.0].map(|tau| (tau, one_step_lambda(tau)));
    let agreements = [
        (
            Scheme::ExactExponential,
            curve_agreement(
                |tau| lambda_curve(Scheme::ExactExponential, tau),
                &samples,
                REL_BAND,
                0.0,
            )
            .max_rel_dev,
        ),
        (
            Scheme::ExplicitEuler,
            curve_agreement(
                |tau| lambda_curve(Scheme::ExplicitEuler, tau),
                &samples,
                REL_BAND,
                0.0,
            )
            .max_rel_dev,
        ),
        (
            Scheme::BackwardEuler,
            curve_agreement(
                |tau| lambda_curve(Scheme::BackwardEuler, tau),
                &samples,
                REL_BAND,
                0.0,
            )
            .max_rel_dev,
        ),
    ];
    let winner = agreements
        .iter()
        .min_by(|a, b| a.1.total_cmp(&b.1))
        .map(|(scheme, _)| *scheme)
        .unwrap();
    (winner, samples, agreements)
}

fn affine_u(pos: [f64; 3]) -> [f64; 3] {
    [
        0.01 + 2.0e-3 * pos[0] + 1.0e-3 * pos[1] - 5.0e-4 * pos[2],
        -0.02 + 1.0e-3 * pos[0],
        0.03 - 2.0e-3 * pos[2],
    ]
}

fn affine_sample(pos: [f64; 3]) -> Sample {
    sample_grid(pos, [12, 10, 8], |x, y, z| {
        (affine_u([x as f64, y as f64, z as f64]), false)
    })
}

#[test]
fn c1_drag_time_integrator_signature_curve() {
    let (winner, samples, agreements) = identify_integrator_signature();
    println!(
        "ACC PART C1: identified_scheme={}, samples={samples:?}, agreements={agreements:?}",
        winner.name()
    );

    for &(tau, lambda) in &samples {
        let want = lambda_curve(winner, tau);
        let rel = (lambda - want).abs() / lambda.abs();
        assert!(
            rel <= REL_BAND,
            "C1 winning curve mismatch at tau_p={tau:e}: measured_lambda={lambda:e}, \
             expected_lambda={want:e}, rel={rel:e}, band={REL_BAND:e}, denominator=measured lambda"
        );
    }

    let tau_sharp = 2.0;
    let measured_sharp = samples
        .iter()
        .find(|(tau, _)| *tau == tau_sharp)
        .map(|(_, lambda)| *lambda)
        .unwrap();
    let winner_dev =
        (measured_sharp - lambda_curve(winner, tau_sharp)).abs() / measured_sharp.abs();
    for scheme in [
        Scheme::ExactExponential,
        Scheme::ExplicitEuler,
        Scheme::BackwardEuler,
    ] {
        if scheme == winner {
            continue;
        }
        let loser_dev =
            (measured_sharp - lambda_curve(scheme, tau_sharp)).abs() / measured_sharp.abs();
        assert!(
            loser_dev >= 10.0 * winner_dev,
            "C1 sharpness failed at tau_p=2: winner={}, loser={}, measured_lambda={measured_sharp:e}, \
             winner_rel_dev={winner_dev:e}, loser_rel_dev={loser_dev:e}, required loser >= 10x winner",
            winner.name(),
            scheme.name()
        );
    }
}

#[test]
fn c2_sub_relaxation_time_behavior() {
    let (scheme, _, _) = identify_integrator_signature();
    let tau_p = 0.4;
    let mut set = ParticleSet::new(
        vec![particle_for_tau(tau_p, [V0, 0.0, 0.0], RHO_P)],
        RHO_F,
        NU,
        [0.0; 3],
    );

    let mut velocities = vec![set.particles[0].vel[0]];
    for _ in 0..20 {
        set.step(still_fluid, None::<fn([f64; 3]) -> f64>);
        velocities.push(set.particles[0].vel[0]);
    }
    let abs_velocities: Vec<f64> = velocities.iter().map(|v| v.abs()).collect();
    let mono = monotonicity(&abs_velocities);
    let sign_flips = velocities.windows(2).filter(|w| w[0] * w[1] < 0.0).count();
    let measured_ratio = velocities[1] / velocities[0];
    let expected_ratio = match scheme {
        // Exact integration of dv/dt = -v/tau_p over dt=1 gives
        // v_{n+1}/v_n = exp(-1/tau_p).
        Scheme::ExactExponential => (-1.0 / tau_p).exp(),
        // Forward Euler gives v_{n+1}/v_n = 1 - 1/tau_p, which is negative
        // for tau_p < 1 and would flip sign here.
        Scheme::ExplicitEuler => 1.0 - 1.0 / tau_p,
        // Backward Euler gives v_{n+1}/v_n = 1/(1 + 1/tau_p), monotone for
        // positive tau_p.
        Scheme::BackwardEuler => 1.0 / (1.0 + 1.0 / tau_p),
    };
    let rel = (measured_ratio - expected_ratio).abs() / expected_ratio.abs();
    println!(
        "ACC PART C2: identified_scheme={}, velocities={velocities:?}, abs_monotonicity={mono:e}, \
         sign_flips={sign_flips}, measured_ratio={measured_ratio:e}, expected_ratio={expected_ratio:e}",
        scheme.name()
    );

    assert_eq!(
        sign_flips, 0,
        "C2 sign-flip stability failed: sign_flips={sign_flips}, velocities={velocities:?}"
    );
    assert!(
        mono == 1.0,
        "C2 monotone |v| decay failed: monotonicity={mono:e}, band exact 1.0, values={abs_velocities:?}"
    );
    assert!(
        rel <= REL_BAND,
        "C2 ratio mismatch for {} at tau_p=0.4: measured_ratio={measured_ratio:e}, \
         expected_ratio={expected_ratio:e}, rel={rel:e}, band={REL_BAND:e}, denominator=expected ratio",
        scheme.name()
    );
}

#[test]
fn c3_grid_sampler_exactness_on_affine_fields() {
    // Trilinear interpolation is a tensor product of linear interpolation along
    // x, y, z. Any affine field u(x) = A x + b is therefore reproduced exactly
    // up to round-off because each component is a sum of constants plus terms
    // linear in one coordinate.
    let positions = [
        [1.25, 1.50, 1.75],
        [2.00, 2.25, 2.50],
        [3.75, 3.00, 2.25],
        [4.50, 4.75, 3.00],
        [5.25, 5.50, 4.75],
        [6.75, 6.25, 5.50],
        [7.00, 7.75, 6.25],
        [8.25, 2.75, 3.50],
        [9.75, 8.50, 5.25],
    ];
    let mut max_abs = 0.0f64;

    for pos in positions {
        let sample = affine_sample(pos);
        let want = affine_u(pos);
        for axis in 0..3 {
            let err = (sample.u[axis] - want[axis]).abs();
            max_abs = max_abs.max(err);
            assert!(
                err <= ABS_BAND,
                "C3 affine sampler mismatch at pos={pos:?}, axis={axis}: sampled={:e}, \
                 analytic={:e}, abs_err={err:e}, band={ABS_BAND:e} absolute",
                sample.u[axis],
                want[axis]
            );
        }
    }
    println!("ACC PART C3: positions={positions:?}, max_abs_err={max_abs:e}, band={ABS_BAND:e}");
}

#[test]
fn c4_sampler_contract_at_solid_nodes() {
    // For a stationary no-slip wall, letting solid nodes contribute zero
    // velocity makes the interpolated velocity approach zero at the wall
    // surface, which is defensible. It is not a complete moving-wall particle
    // coupling contract because the accessor has no wall-velocity channel:
    // SPEC-GAP for moving-wall particle coupling, recorded here without failing
    // this stationary-wall probe.
    let xs = [4.5, 5.0, 5.25, 5.5, 5.75];
    let mut measured = Vec::new();
    for x in xs {
        let sample = sample_grid([x, 2.0, 2.0], [12, 4, 4], |ix, _, _| {
            if ix >= 6 {
                ([0.0, 0.0, 0.0], true)
            } else {
                ([0.02, 0.0, 0.0], false)
            }
        });
        measured.push((x, sample.u[0], sample.solid));
    }
    println!("ACC PART C4: measured_wall_probe={measured:?}");

    for &(x, ux, solid) in &measured {
        let expected_weight = if x < 5.0 { 1.0 } else { 6.0 - x };
        let expected_ux = 0.02 * expected_weight;
        let expected_solid = false;
        let err = (ux - expected_ux).abs();
        assert!(
            err <= ABS_BAND,
            "C4 solid-node velocity contract drift at x={x:e}: measured_ux={ux:e}, \
             expected_ux={expected_ux:e}, fluid_weight={expected_weight:e}, abs_err={err:e}, \
             band={ABS_BAND:e} absolute"
        );
        assert_eq!(
            solid, expected_solid,
            "C4 solid flag contract drift at x={x:e}: measured_solid={solid}, \
             expected_solid={expected_solid}; sampled values={:?}",
            measured
        );
    }
}

#[test]
fn c5_buoyancy_sign_antisymmetry() {
    // In the Stokes limit, terminal velocity is
    // v_t = tau_p g (1 - rho_f/rho_p), with
    // tau_p = rho_p d^2 / (18 rho_f nu). Multiplying out gives
    // v_t = d^2 g (rho_p - rho_f) / (18 rho_f nu), so choosing rho_p=1.5
    // and rho_p=0.5 makes the two terminal velocities exactly opposite.
    // Schiller-Naumann drag depends on |v_t| through Re_p, so the identity is
    // exact only as Re_p -> 0; here Re_p is tiny enough that the correction is
    // below the stated Stokes-limit tolerance.
    let d = 1.0;
    let mut down = ParticleSet::new(
        vec![Particle {
            pos: [4.0, 4.0, 4.0],
            vel: [0.0; 3],
            d,
            rho_p: 1.5,
            exposure: 0.0,
        }],
        RHO_F,
        NU,
        [0.0, 0.0, -1.0e-6],
    );
    let mut up = ParticleSet::new(
        vec![Particle {
            pos: [4.0, 4.0, 4.0],
            vel: [0.0; 3],
            d,
            rho_p: 0.5,
            exposure: 0.0,
        }],
        RHO_F,
        NU,
        [0.0, 0.0, -1.0e-6],
    );
    for _ in 0..2_000 {
        down.step(still_fluid, None::<fn([f64; 3]) -> f64>);
        up.step(still_fluid, None::<fn([f64; 3]) -> f64>);
    }

    let v_down = down.particles[0].vel[2];
    let v_up = up.particles[0].vel[2];
    let antisym = (v_up + v_down).abs();
    let denom = v_down.abs();
    let rel = antisym / denom;
    println!(
        "ACC PART C5: v_down={v_down:e}, v_up={v_up:e}, antisym_abs={antisym:e}, \
         rel={rel:e}, denominator=|v_down|"
    );
    assert!(
        rel <= 1.0e-6,
        "C5 buoyancy antisymmetry failed: v_down={v_down:e}, v_up={v_up:e}, \
         |sum|={antisym:e}, rel={rel:e}, band=1e-6, denominator=|v_down|={denom:e}"
    );
}

#[test]
fn c6_step_vs_step_depositing_degeneracy() {
    let particles = vec![
        Particle {
            pos: [1.25, 2.50, 3.75],
            vel: [0.015, -0.002, 0.004],
            d: 0.4,
            rho_p: 1.4,
            exposure: 0.0,
        },
        Particle {
            pos: [4.50, 5.25, 2.25],
            vel: [-0.008, 0.006, -0.003],
            d: 0.7,
            rho_p: 2.2,
            exposure: 0.0,
        },
        Particle {
            pos: [8.00, 3.75, 5.50],
            vel: [0.002, -0.005, 0.001],
            d: 1.1,
            rho_p: 0.8,
            exposure: 0.0,
        },
    ];
    let mut stepping =
        ParticleSet::new(particles.clone(), RHO_F, NU, [0.0, 0.0, -2.0e-6]).with_restitution(0.25);
    let mut depositing =
        ParticleSet::new(particles, RHO_F, NU, [0.0, 0.0, -2.0e-6]).with_restitution(0.25);
    let exposure = |pos: [f64; 3]| 0.01 + 1.0e-4 * pos[0] - 2.0e-4 * pos[1] + 3.0e-4 * pos[2];
    let mut deposits = Vec::<DepositEvent>::new();

    for step in 0..200 {
        stepping.step(affine_sample, Some(exposure));
        depositing.step_depositing(affine_sample, Some(exposure), -1.0e9, &mut deposits);
        assert_eq!(
            stepping.particles, depositing.particles,
            "C6 step/step_depositing trajectory drift at step {step}: step={:?}, depositing={:?}",
            stepping.particles, depositing.particles
        );
        assert!(
            deposits.is_empty(),
            "C6 unreachable floor emitted deposits at step {step}: deposits={deposits:?}"
        );
    }
    println!(
        "ACC PART C6: steps=200, final_particles={:?}, deposits={:?}",
        stepping.particles, deposits
    );
}
