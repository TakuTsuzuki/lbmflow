// Inherited verbatim from the retired V1 suite at its retirement (2026-07-05,
// scripts/sync-tests.sh mechanical retarget); now the canonical facade tests.
//! Validation T9: zero-gradient outflow robustness for a cylinder wake.

use lbm_core::compat::lattice::CS2;
use lbm_core::compat::prelude::*;

#[derive(Clone, Copy, Debug)]
struct OutflowCase {
    nx: usize,
    ny: usize,
    d: f64,
    cx: f64,
    cy: f64,
    u: f64,
    re: f64,
    steps: usize,
}

fn build_case(case: OutflowCase) -> Simulation<f64> {
    build_case_with_right(case, EdgeBC::Outflow)
}

fn build_case_with_right(case: OutflowCase, right: EdgeBC<f64>) -> Simulation<f64> {
    let nu = case.u * case.d / case.re;
    let mut sim: Simulation<f64> = SimConfig {
        nx: case.nx,
        ny: case.ny,
        nu,
        collision: Collision::Trt { magic: 3.0 / 16.0 },
        edges: Edges {
            left: EdgeBC::VelocityInlet { u: [case.u, 0.0] },
            right,
            bottom: EdgeBC::BounceBack,
            top: EdgeBC::BounceBack,
        },
        ..Default::default()
    }
    .build()
    .unwrap();
    sim.set_inlet_profile(Edge::Left, |y| {
        if y == 0 || y == case.ny - 1 {
            [0.0, 0.0]
        } else {
            [case.u, 0.0]
        }
    });
    let r = 0.5 * case.d;
    let is_cylinder = |x: usize, y: usize| {
        let dx = x as f64 - case.cx;
        let dy = y as f64 - case.cy;
        dx * dx + dy * dy <= r * r
    };
    sim.set_solid_region(is_cylinder);
    sim.init_with(|x, y| {
        if is_cylinder(x, y) || x == 0 || x == case.nx - 1 || y == 0 || y == case.ny - 1 {
            (1.0, 0.0, 0.0)
        } else {
            (1.0, case.u, 0.0)
        }
    });
    sim
}

fn assert_finite(sim: &Simulation<f64>, label: &str) {
    for y in 0..sim.ny() {
        for x in 0..sim.nx() {
            let rho = sim.rho(x, y);
            let ux = sim.ux(x, y);
            let uy = sim.uy(x, y);
            assert!(
                rho.is_finite() && ux.is_finite() && uy.is_finite(),
                "T9 {label}: non-finite at ({x},{y}), rho = {rho:e}, ux = {ux:e}, uy = {uy:e}, time = {}",
                sim.time()
            );
        }
    }
}

fn backflow_fraction(sim: &Simulation<f64>) -> (f64, f64, f64) {
    let x_in = 1;
    let x_out = sim.nx() - 2;
    let mut inflow = 0.0;
    let mut backflow = 0.0;
    for y in 1..=(sim.ny() - 2) {
        inflow += (sim.rho(x_in, y) * sim.ux(x_in, y)).max(0.0);
        backflow += (-sim.rho(x_out, y) * sim.ux(x_out, y)).max(0.0);
    }
    (backflow / inflow, backflow, inflow)
}

fn pressure_rms_ratio(sim: &Simulation<f64>) -> (f64, f64, f64) {
    let near_start = sim.nx() * 9 / 10;
    let mid_start = sim.nx() * 45 / 100;
    let mid_end = sim.nx() * 55 / 100;
    let near = pressure_rms(sim, near_start, sim.nx() - 2);
    let mid = pressure_rms(sim, mid_start, mid_end);
    (near / mid, near, mid)
}

fn pressure_rms(sim: &Simulation<f64>, x0: usize, x1: usize) -> f64 {
    let mut vals = Vec::new();
    for x in x0..=x1 {
        for y in 1..=(sim.ny() - 2) {
            if !sim.is_solid(x, y) {
                vals.push(CS2 * sim.rho(x, y));
            }
        }
    }
    let mean = vals.iter().sum::<f64>() / vals.len() as f64;
    (vals.iter().map(|v| (v - mean) * (v - mean)).sum::<f64>() / vals.len() as f64).sqrt()
}

fn run_outflow_case(case: OutflowCase, right: EdgeBC<f64>, label: &str) -> (f64, f64, f64, f64) {
    let mut sim = build_case_with_right(case, right);
    sim.run(case.steps);
    assert_finite(&sim, label);
    let (frac, backflow, inflow) = backflow_fraction(&sim);
    let (ratio, near, mid) = pressure_rms_ratio(&sim);
    eprintln!(
        "{label}: backflow_frac={frac:.8e}, pressure_rms_ratio={ratio:.8e}, near={near:.8e}, mid={mid:.8e}"
    );
    (frac, backflow, inflow, ratio)
}

#[test]
fn t9_outflow_cylinder_wake_smoke_stays_finite_with_limited_backflow() {
    let case = OutflowCase {
        nx: 260,
        ny: 88,
        d: 16.0,
        cx: 64.0,
        cy: 44.0,
        u: 0.05,
        re: 100.0,
        steps: 20_000,
    };
    let mut sim = build_case(case);
    sim.run(case.steps);
    assert_finite(&sim, "smoke");
    let (frac, backflow, inflow) = backflow_fraction(&sim);
    let (ratio, near, mid) = pressure_rms_ratio(&sim);
    eprintln!(
        "T9 outflow smoke: backflow_frac={frac:.8e}, pressure_rms_ratio={ratio:.8e}, near={near:.8e}, mid={mid:.8e}"
    );
    assert!(
        frac <= 0.05,
        "T9 smoke backflow fraction = {frac:e}, backflow = {backflow:e}, inflow = {inflow:e}, steps = {}",
        case.steps
    );
}

#[test]
#[ignore]
fn t9_outflow_cylinder_wake_long_run_stays_sane() {
    let case = OutflowCase {
        nx: 440,
        ny: 160,
        d: 20.0,
        cx: 110.0,
        cy: 81.0,
        u: 0.05,
        re: 100.0,
        steps: 100_000,
    };
    let mut sim = build_case(case);
    sim.run(case.steps);
    assert_finite(&sim, "long");
    let (frac, backflow, inflow) = backflow_fraction(&sim);
    assert!(
        frac <= 0.05,
        "T9 long backflow fraction = {frac:e}, backflow = {backflow:e}, inflow = {inflow:e}, steps = {}",
        case.steps
    );
    let (ratio, near, mid) = pressure_rms_ratio(&sim);
    // VALIDATION.md T9 allows up to 15x: zero-gradient outflow partially
    // reflects pressure waves, with measured ratio around 11.3 on this case.
    assert!(
        ratio <= 15.0,
        "T9 pressure RMS ratio = {ratio:e}, near = {near:e}, mid = {mid:e}, steps = {}",
        case.steps
    );
}

#[test]
fn t9b_convective_outflow_cylinder_wake_stays_finite_with_limited_backflow() {
    let case = OutflowCase {
        nx: 260,
        ny: 88,
        d: 16.0,
        cx: 64.0,
        cy: 44.0,
        u: 0.05,
        re: 100.0,
        steps: 20_000,
    };
    let (frac, backflow, inflow, ratio) = run_outflow_case(
        case,
        EdgeBC::ConvectiveOutflow { u_conv: case.u },
        "T9b convective smoke",
    );
    assert!(
        frac <= 0.05,
        "T9b convective backflow fraction = {frac:e}, backflow = {backflow:e}, inflow = {inflow:e}, steps = {}",
        case.steps
    );
    // Measured 2026-07-05 on the T9 smoke geometry: ConvectiveOutflow
    // ratio = 7.96079611, while the matching Outflow ratio = 5.81418778.
    // This freezes observed behaviour without requiring convective to win.
    assert!(
        (ratio - 7.96079611).abs() <= 0.25,
        "T9b convective pressure RMS ratio = {ratio:e}, frozen = 7.96079611, steps = {}",
        case.steps
    );
}

#[test]
#[ignore]
fn t9b_convective_outflow_long_run_pressure_ratio_is_frozen() {
    let case = OutflowCase {
        nx: 440,
        ny: 160,
        d: 20.0,
        cx: 110.0,
        cy: 81.0,
        u: 0.05,
        re: 100.0,
        steps: 100_000,
    };
    let mut outflow = build_case(case);
    outflow.run(case.steps);
    assert_finite(&outflow, "T9b outflow comparison long");
    let (out_ratio, out_near, out_mid) = pressure_rms_ratio(&outflow);
    eprintln!(
        "T9b outflow comparison long: pressure_rms_ratio={out_ratio:.8e}, near={out_near:.8e}, mid={out_mid:.8e}"
    );
    let (frac, backflow, inflow, conv_ratio) = run_outflow_case(
        case,
        EdgeBC::ConvectiveOutflow { u_conv: case.u },
        "T9b convective long",
    );
    assert!(
        frac <= 0.05,
        "T9b long convective backflow fraction = {frac:e}, backflow = {backflow:e}, inflow = {inflow:e}, steps = {}",
        case.steps
    );
    // Measured 2026-07-05 on this T9 long geometry: ConvectiveOutflow
    // ratio = 0.714127152, while the matching Outflow ratio = 11.3253818.
    // This records reality for this geometry without making a general claim
    // that convective is always better.
    assert!(
        (conv_ratio - 0.714127152).abs() <= 0.05,
        "T9b long convective pressure RMS ratio = {conv_ratio:e}, frozen = 0.714127152; Outflow ratio = {out_ratio:e}, near = {out_near:e}, mid = {out_mid:e}, steps = {}",
        case.steps
    );
}

#[derive(Clone, Copy, Debug)]
struct ReflectionMeasurement {
    raw_r: f64,
    emitted_amp: f64,
    transmitted_amp: f64,
    reflected_amp: f64,
}

#[derive(Clone, Copy, Debug)]
struct ReflectionCurvePoint {
    u_conv: f64,
    r: f64,
    raw_r: f64,
    transmitted_amp: f64,
}

const REFLECTION_NX: usize = 128;
const REFLECTION_NY: usize = 32;
const REFLECTION_TAU: f64 = 0.8;
const REFLECTION_NU: f64 = (REFLECTION_TAU - 0.5) / 3.0;
const REFLECTION_PULSE_X: f64 = 40.0;
const REFLECTION_PULSE_WIDTH: f64 = 8.0;
const REFLECTION_ETA: f64 = 0.01;
const REFLECTION_WARMUP_STEPS: usize = 500;
const REFLECTION_SAMPLE_STEPS: usize = 420;

fn build_reflection_channel(right: EdgeBC<f64>) -> Simulation<f64> {
    let mut sim: Simulation<f64> = SimConfig {
        nx: REFLECTION_NX,
        ny: REFLECTION_NY,
        nu: REFLECTION_NU,
        collision: Collision::Trt { magic: 3.0 / 16.0 },
        edges: Edges {
            left: EdgeBC::VelocityInlet { u: [0.0, 0.0] },
            right,
            bottom: EdgeBC::BounceBack,
            top: EdgeBC::BounceBack,
        },
        ..Default::default()
    }
    .build()
    .unwrap();
    sim.set_inlet_profile(Edge::Left, |_| [0.0, 0.0]);
    sim.init_with(|_, _| (1.0, 0.0, 0.0));
    sim
}

fn reflection_pulse_mode(x: usize, center: f64) -> f64 {
    let dx = (x as f64 - center) / REFLECTION_PULSE_WIDTH;
    (-dx * dx).exp()
}

fn inject_right_going_density_pulse(sim: &mut Simulation<f64>) {
    // Linear acoustic eigenmode for a right-going small-amplitude packet:
    // p' = cs^2 rho', and the 1-D characteristic relation is u' = p'/(rho0 cs)
    // = cs rho' for rho0=1. This avoids launching an equal left-going packet
    // that would reflect from the inlet before the outlet-reflection window.
    let cs = CS2.sqrt();
    sim.init_with(|x, y| {
        if y == 0 || y == REFLECTION_NY - 1 {
            (1.0, 0.0, 0.0)
        } else {
            let delta_rho = REFLECTION_ETA * reflection_pulse_mode(x, REFLECTION_PULSE_X);
            (1.0 + delta_rho, cs * delta_rho, 0.0)
        }
    });
}

fn projected_density_amplitude(sim: &Simulation<f64>, center: f64, x0: usize, x1: usize) -> f64 {
    let mut count = 0.0;
    let mut rho_sum = 0.0;
    let mut mode_sum = 0.0;
    for x in x0..=x1 {
        let mode = reflection_pulse_mode(x, center);
        for y in 1..=(sim.ny() - 2) {
            rho_sum += sim.rho(x, y) - 1.0;
            mode_sum += mode;
            count += 1.0;
        }
    }
    let rho_mean = rho_sum / count;
    let mode_mean = mode_sum / count;

    let mut num = 0.0;
    let mut den = 0.0;
    for x in x0..=x1 {
        let mode = reflection_pulse_mode(x, center) - mode_mean;
        for y in 1..=(sim.ny() - 2) {
            let rho = sim.rho(x, y) - 1.0 - rho_mean;
            num += rho * mode;
            den += mode * mode;
        }
    }
    num / den
}

fn measure_reflection(right: EdgeBC<f64>, label: &str) -> ReflectionMeasurement {
    let mut sim = build_reflection_channel(right);
    sim.run(REFLECTION_WARMUP_STEPS);
    inject_right_going_density_pulse(&mut sim);

    let emitted_amp = projected_density_amplitude(&sim, REFLECTION_PULSE_X, 12, 68).abs();
    let cs = CS2.sqrt();
    let outlet_probe_x = 112.0;
    let t_transmit = ((outlet_probe_x - REFLECTION_PULSE_X) / cs).round() as usize;
    let t_return = (2.0 * ((REFLECTION_NX - 1) as f64 - REFLECTION_PULSE_X) / cs).round() as usize;
    let mut transmitted_amp = 0.0_f64;
    let mut reflected_amp = 0.0_f64;

    for step in 1..=REFLECTION_SAMPLE_STEPS {
        sim.step();
        assert_finite(&sim, label);
        let transmitted = projected_density_amplitude(&sim, outlet_probe_x, 84, 126).abs();
        if step.abs_diff(t_transmit) <= 60 {
            transmitted_amp = transmitted_amp.max(transmitted);
        }
        let reflected = projected_density_amplitude(&sim, REFLECTION_PULSE_X, 12, 68).abs();
        if step.abs_diff(t_return) <= 70 {
            reflected_amp = reflected_amp.max(reflected);
        }
    }

    ReflectionMeasurement {
        raw_r: reflected_amp / emitted_amp,
        emitted_amp,
        transmitted_amp,
        reflected_amp,
    }
}

#[test]
#[ignore = "lane 1.7 adversarial probe: current ConvectiveOutflow reflects near hard-wall in this rest-channel pulse setup"]
fn t9b_convective_outflow_reflection_coefficient_curve_has_better_than_outflow_regime() {
    // Boundary reflection coefficient derivation:
    // a right-going linear acoustic pulse has delta u = cs * delta rho. The
    // outlet reflection is the returning packet amplitude divided by the
    // emitted packet amplitude. Because this finite-viscosity packet damps
    // during the round trip, we divide each raw ratio by the same-channel hard
    // wall raw ratio; a hard wall has physical R=1, so this calibration removes
    // propagation damping without using implementation internals.
    let hard_wall = measure_reflection(EdgeBC::BounceBack, "T9b reflection hard wall");
    let outflow = measure_reflection(EdgeBC::Outflow, "T9b reflection Outflow");
    let baseline_r = outflow.raw_r / hard_wall.raw_r;
    let u_convs = [0.05, 0.1, 0.2, 0.3, 0.5, 1.0];
    let mut curve = Vec::new();

    for u_conv in u_convs {
        let measured = measure_reflection(
            EdgeBC::ConvectiveOutflow { u_conv },
            "T9b reflection ConvectiveOutflow",
        );
        curve.push(ReflectionCurvePoint {
            u_conv,
            r: measured.raw_r / hard_wall.raw_r,
            raw_r: measured.raw_r,
            transmitted_amp: measured.transmitted_amp,
        });
    }

    let min_point = curve
        .iter()
        .min_by(|a, b| a.r.total_cmp(&b.r))
        .copied()
        .unwrap();
    println!(
        "T9b reflection curve: hard_wall_raw_R={:.8e}, emitted_amp={:.8e}, transmitted_amp={:.8e}, reflected_amp={:.8e}",
        hard_wall.raw_r, hard_wall.emitted_amp, hard_wall.transmitted_amp, hard_wall.reflected_amp
    );
    println!(
        "T9b reflection baseline Outflow: baseline_R={:.8e}, raw_R={:.8e}, transmitted_amp={:.8e}",
        baseline_r, outflow.raw_r, outflow.transmitted_amp
    );
    for point in &curve {
        println!(
            "T9b reflection ConvectiveOutflow: u_conv={:.2}, R={:.8e}, raw_R={:.8e}, transmitted_amp={:.8e}",
            point.u_conv, point.r, point.raw_r, point.transmitted_amp
        );
    }
    println!(
        "T9b reflection summary: baseline_R={:.8e}, min_R={:.8e}, u_conv_at_min={:.2}",
        baseline_r, min_point.r, min_point.u_conv
    );

    assert!(
        min_point.r < baseline_r,
        "T9b reflection curve has no better-than-zero-gradient regime: min_R={:.8e} at u_conv={:.2}, baseline_R={:.8e}",
        min_point.r,
        min_point.u_conv,
        baseline_r
    );
    assert!(
        curve[0].r >= curve[1].r && curve[1].r >= curve[2].r,
        "T9b reflection low-u branch does not approach hard-wall R=1 monotonically as u_conv -> 0: R(0.05)={:.8e}, R(0.10)={:.8e}, R(0.20)={:.8e}",
        curve[0].r,
        curve[1].r,
        curve[2].r
    );
    assert!(
        curve[0].r > 0.5,
        "T9b reflection smallest u_conv is not in the high-reflection approach-to-wall regime: R(0.05)={:.8e}",
        curve[0].r
    );
}
