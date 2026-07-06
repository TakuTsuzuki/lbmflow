// Inherited verbatim from the retired V1 suite at its retirement (2026-07-05,
// scripts/sync-tests.sh mechanical retarget); now the canonical facade tests.
//! Shared helpers for validation tests.
#![allow(dead_code)]

use lbm_core::compat::prelude::*;

pub mod metrics;
pub mod tgv_analysis;
#[allow(unused_imports)]
pub use metrics::l2_rel;

/// Run until the velocity field stops changing: returns `true` once
/// `max |u - u_prev| <= tol * max |u|` between checks `check_every` steps
/// apart, `false` if `max_steps` elapsed first.
pub fn run_to_steady(
    sim: &mut Simulation<f64>,
    check_every: usize,
    tol: f64,
    max_steps: usize,
) -> bool {
    let mut prev: Vec<f64> = Vec::new();
    let mut elapsed = 0;
    while elapsed < max_steps {
        sim.run(check_every);
        elapsed += check_every;
        let cur: Vec<f64> = sim
            .ux_field()
            .iter()
            .chain(sim.uy_field())
            .copied()
            .collect();
        if !prev.is_empty() {
            let mut dmax = 0.0f64;
            let mut umax = 0.0f64;
            for (c, p) in cur.iter().zip(&prev) {
                dmax = dmax.max((c - p).abs());
                umax = umax.max(c.abs());
            }
            if umax > 0.0 && dmax <= tol * umax {
                return true;
            }
        }
        prev = cur;
    }
    false
}
