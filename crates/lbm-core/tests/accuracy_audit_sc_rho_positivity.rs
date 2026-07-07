//! ACC-AUDIT radar #36 / code-to-spec-diff A22: Shan-Chen density positivity.
//!
//! A22 flags that `compat/multiphase.rs` may evaluate the Shan-Chen
//! pseudopotential at non-positive density. That is outside the model domain:
//! `Psi::Exponential` divides by rho, and a non-positive rho can silently turn
//! the force into zero, infinity, or NaN. These probes characterize both the
//! validated T11 flat-interface envelope and an intentionally aggressive
//! coexistence-boundary setup, then check that the source carries an explicit
//! positivity guard.

use lbm_core::compat::multiphase::ShanChen;
use lbm_core::compat::prelude::*;

const NX: usize = 64;
const NY: usize = 128;
const G: f64 = -5.0;
const NU_TAU_1: f64 = 1.0 / 6.0;
const RHO_L_INIT: f64 = 2.0;
const RHO_V_T11: f64 = 0.1194;
const RHO_V_STABLE_INIT: f64 = 0.15;
const RHO_V_AGGRESSIVE_INIT: f64 = 0.05;
const G1_STEPS: usize = 30_000;
const G2_STEPS: usize = 5_000;
const G1_MIN_RHO_THRESHOLD: f64 = 0.5 * RHO_V_T11;
const G2_POSITIVITY_FLOOR: f64 = 0.01;

#[derive(Clone, Copy, Debug)]
struct RhoMinimum {
    value: f64,
    step: usize,
    x: usize,
    y: usize,
}

impl RhoMinimum {
    fn new() -> Self {
        Self {
            value: f64::INFINITY,
            step: 0,
            x: 0,
            y: 0,
        }
    }
}

#[derive(Debug)]
struct RhoTrace {
    min: RhoMinimum,
    first_non_finite: Option<RhoMinimum>,
    first_non_positive: Option<RhoMinimum>,
    trajectory: Vec<RhoMinimum>,
}

fn build_flat_interface(vapor_rho: f64) -> Simulation<f64> {
    let mut sim: Simulation<f64> = SimConfig {
        nx: NX,
        ny: NY,
        nu: NU_TAU_1,
        ..Default::default()
    }
    .build()
    .unwrap();
    sim.init_with(|_, y| {
        let liquid = y >= NY / 4 && y < 3 * NY / 4;
        (if liquid { RHO_L_INIT } else { vapor_rho }, 0.0, 0.0)
    });
    sim
}

fn field_min(sim: &Simulation<f64>, step: usize) -> RhoMinimum {
    let mut out = RhoMinimum::new();
    for (i, &rho) in sim.rho_field().iter().enumerate() {
        if rho < out.value || (!rho.is_finite() && out.value.is_finite()) {
            out = RhoMinimum {
                value: rho,
                step,
                x: i % NX,
                y: i / NX,
            };
        }
    }
    out
}

fn push_minimum_sample(samples: &mut Vec<RhoMinimum>, sample: RhoMinimum) {
    samples.push(sample);
    const KEEP: usize = 16;
    if samples.len() > KEEP {
        samples.remove(0);
    }
}

fn run_and_trace_minimum(mut sim: Simulation<f64>, steps: usize) -> RhoTrace {
    let sc = ShanChen::new(G);
    let mut min = field_min(&sim, 0);
    let mut first_non_finite = (!min.value.is_finite()).then_some(min);
    let mut first_non_positive = (min.value <= 0.0).then_some(min);
    let mut trajectory = vec![min];

    for step in 1..=steps {
        sim.force_field_mut().fill([0.0; 2]);
        sc.update_force(&mut sim);
        sim.step();

        let candidate = field_min(&sim, step);
        if !candidate.value.is_finite() && first_non_finite.is_none() {
            first_non_finite = Some(candidate);
        }
        if candidate.value <= 0.0 && first_non_positive.is_none() {
            first_non_positive = Some(candidate);
        }
        if candidate.value < min.value || !candidate.value.is_finite() {
            min = candidate;
            push_minimum_sample(&mut trajectory, candidate);
        }
    }

    RhoTrace {
        min,
        first_non_finite,
        first_non_positive,
        trajectory,
    }
}

#[test]
fn g1_t11_flat_interface_keeps_density_above_half_vapor_reference() {
    let trace = run_and_trace_minimum(build_flat_interface(RHO_V_STABLE_INIT), G1_STEPS);
    println!(
        "G1 T11 flat interface rho positivity: steps={G1_STEPS}, rho_v_ref={RHO_V_T11:.8}, threshold={G1_MIN_RHO_THRESHOLD:.8}, min_rho={:.12e} at step={} loc=({},{}), first_non_positive={:?}, first_non_finite={:?}, trajectory={:?}",
        trace.min.value,
        trace.min.step,
        trace.min.x,
        trace.min.y,
        trace.first_non_positive,
        trace.first_non_finite,
        trace.trajectory
    );
    assert!(
        trace.first_non_finite.is_none(),
        "G1 T11 flat interface produced non-finite rho: {:?}; trajectory={:?}",
        trace.first_non_finite,
        trace.trajectory
    );
    assert!(
        trace.min.value >= G1_MIN_RHO_THRESHOLD,
        "G1 T11 flat interface rho positivity failed: min_rho={:.12e} at step={} loc=({},{}), threshold=0.5*rho_v_ref={G1_MIN_RHO_THRESHOLD:.12e}, rho_v_ref={RHO_V_T11:.8}, trajectory={:?}",
        trace.min.value,
        trace.min.step,
        trace.min.x,
        trace.min.y,
        trace.trajectory
    );
}

#[test]
fn g2_aggressive_flat_interface_does_not_cross_hard_positivity_floor() {
    let trace = run_and_trace_minimum(build_flat_interface(RHO_V_AGGRESSIVE_INIT), G2_STEPS);
    println!(
        "G2 aggressive flat interface rho positivity: steps={G2_STEPS}, initial_liquid={RHO_L_INIT:.8}, initial_vapor={RHO_V_AGGRESSIVE_INIT:.8}, floor={G2_POSITIVITY_FLOOR:.8}, min_rho={:.12e} at step={} loc=({},{}), first_non_positive={:?}, first_non_finite={:?}, trajectory={:?}",
        trace.min.value,
        trace.min.step,
        trace.min.x,
        trace.min.y,
        trace.first_non_positive,
        trace.first_non_finite,
        trace.trajectory
    );
    assert!(
        trace.first_non_finite.is_none(),
        "G2 aggressive flat interface produced non-finite rho before/after crossing the SC psi domain: {:?}; trajectory={:?}",
        trace.first_non_finite,
        trace.trajectory
    );
    assert!(
        trace.first_non_positive.is_none(),
        "G2 aggressive flat interface crossed non-positive rho, which is a Shan-Chen positivity-guard failure: {:?}; trajectory={:?}",
        trace.first_non_positive,
        trace.trajectory
    );
    assert!(
        trace.min.value > G2_POSITIVITY_FLOOR,
        "G2 aggressive flat interface rho positivity failed: min_rho={:.12e} at step={} loc=({},{}), floor={G2_POSITIVITY_FLOOR:.12e}, trajectory={:?}",
        trace.min.value,
        trace.min.step,
        trace.min.x,
        trace.min.y,
        trace.trajectory
    );
}

fn parse_positive_numeric_literal(line: &str) -> Option<f64> {
    let bytes = line.as_bytes();
    let mut i = 0;
    while i < bytes.len() {
        let c = bytes[i] as char;
        let starts_number = c.is_ascii_digit()
            || (c == '.' && i + 1 < bytes.len() && (bytes[i + 1] as char).is_ascii_digit());
        if !starts_number {
            i += 1;
            continue;
        }
        let start = i;
        i += 1;
        while i < bytes.len() {
            let c = bytes[i] as char;
            if c.is_ascii_digit()
                || c == '.'
                || c == '_'
                || c == 'e'
                || c == 'E'
                || c == '+'
                || c == '-'
            {
                i += 1;
            } else {
                break;
            }
        }
        if let Ok(value) = line[start..i].replace('_', "").parse::<f64>() {
            if value > 0.0 {
                return Some(value);
            }
        }
    }
    None
}

#[test]
fn g3_multiphase_source_has_explicit_density_floor_or_guard() {
    let source = include_str!("../src/compat/multiphase.rs");
    let guard_lines: Vec<_> = source
        .lines()
        .enumerate()
        .filter_map(|(line_idx, line)| {
            let line_has_rho = line.contains("rho") || line.contains('ρ');
            let line_has_floor = line.contains(".max(") || line.contains(".clamp(");
            (line_has_rho && line_has_floor).then_some((line_idx + 1, line.trim()))
        })
        .collect();
    let debug_assert_lines: Vec<_> = source
        .lines()
        .enumerate()
        .filter_map(|(line_idx, line)| {
            let line_has_rho = line.contains("rho") || line.contains('ρ');
            let line_asserts_positive =
                line.contains("debug_assert") && line_has_rho && line.contains("> 0");
            line_asserts_positive.then_some((line_idx + 1, line.trim()))
        })
        .collect();

    println!("G3 SC rho guard candidates (.max/.clamp): {guard_lines:?}");
    println!("G3 SC rho debug-assert candidates: {debug_assert_lines:?}");

    let positive_floors: Vec<_> = guard_lines
        .iter()
        .filter_map(|&(line, text)| {
            parse_positive_numeric_literal(text).map(|floor| (line, text, floor))
        })
        .collect();
    let physically_meaningful_floor = positive_floors
        .iter()
        .any(|(_, _, floor)| *floor > 0.0 && *floor <= G2_POSITIVITY_FLOOR);

    assert!(
        !debug_assert_lines.is_empty() || physically_meaningful_floor,
        "G3 A22 SC rho guard missing or not at a positive floor: guard_lines={guard_lines:?}, parsed_positive_floors={positive_floors:?}, debug_assert_lines={debug_assert_lines:?}; expected a rho .max/.clamp floor in (0,{G2_POSITIVITY_FLOOR:e}] or a debug_assert!(rho > 0) before Psi::eval"
    );
}
