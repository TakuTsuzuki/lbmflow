//! TGV (Taylor-Green vortex) analysis observables — one-signal-in /
//! one-scalar-out extractors used by validation and characterization tests.
//!
//! Charter distinction (see .claude/skills/lbmflow-accuracy-audit/references/
//! metrics-api.md): these are OBSERVABLES, not agreement metrics — they take
//! a field/signal and return a scalar. The agreement metrics in
//! `common::metrics` compare two quantities to each other. That's why these
//! live under `common::tgv_analysis` and not `common::metrics`, and why
//! there is no Python mirror / drift-guard for them.
//!
//! Second-caller check (metrics-api.md promotion rule adapted): the
//! consolidated helpers below already have two independent callers in the
//! current suite (2D `Simulation<f64>` in validation_tgv.rs / smoke_tgv.rs;
//! 3D `Solver<D3Q19, ...>` in d3q19_smoke.rs), which is the trigger to lift
//! them out of inline copies — not a speculative abstraction.

#![allow(dead_code)]

/// Kinetic-energy-proxy of a velocity field: `Σ |u|²` over the field
/// samples. Units: whatever the caller's velocity is squared, times the
/// number of cells — so this is a PROXY suitable for RATIOS
/// (`ke(t2) / ke(t1)`), NOT for absolute-energy comparisons across grids.
///
/// Cross-caller expectations:
///
/// - 2D compat `Simulation<f64>`: pass `sim.ux_field()` and `sim.uy_field()`.
/// - 3D generic `Solver<D3Q19, f64, ...>`: pass `s.gather_ux()`,
///   `s.gather_uy()`, `s.gather_uz()`.
///
/// The variadic-in-dim signature is materialized as two functions to keep
/// callers explicit about which case they are in — LBM tests are always in
/// exactly 2D or 3D, so an over-general dispatch would only obscure the site.
pub fn ke2d(ux: &[f64], uy: &[f64]) -> f64 {
    assert_eq!(ux.len(), uy.len());
    ux.iter()
        .zip(uy)
        .map(|(a, b)| a * a + b * b)
        .sum()
}

pub fn ke3d(ux: &[f64], uy: &[f64], uz: &[f64]) -> f64 {
    assert_eq!(ux.len(), uy.len());
    assert_eq!(ux.len(), uz.len());
    ux.iter()
        .zip(uy)
        .zip(uz)
        .map(|((a, b), c)| a * a + b * b + c * c)
        .sum()
}

/// Max speed of a velocity field: `max sqrt(ux² + uy² [+ uz²])`. Used as a
/// stability observable (max_speed staying finite = no divergence) and as a
/// low-Ma sanity check (max_speed < 0.15 · cs per the tune-stability
/// thresholds).
pub fn max_speed2d(ux: &[f64], uy: &[f64]) -> f64 {
    assert_eq!(ux.len(), uy.len());
    ux.iter()
        .zip(uy)
        .map(|(a, b)| (a * a + b * b).sqrt())
        .fold(0.0f64, f64::max)
}

pub fn max_speed3d(ux: &[f64], uy: &[f64], uz: &[f64]) -> f64 {
    assert_eq!(ux.len(), uy.len());
    assert_eq!(ux.len(), uz.len());
    ux.iter()
        .zip(uy)
        .zip(uz)
        .map(|((a, b), c)| (a * a + b * b + c * c).sqrt())
        .fold(0.0f64, f64::max)
}

/// Fitted effective kinematic viscosity from the exponential decay of the
/// single-mode Taylor-Green vortex.
///
/// Derivation. For an incompressible flow initialized with a single Fourier
/// mode of wavevector K, the linear-viscous prediction is `u(x,t) =
/// u₀(x) · exp(-ν |K|² t)`, so the kinetic-energy proxy `ke = Σ |u|²`
/// decays at rate `-2 ν |K|² t`:
///
///   ke(t₁) / ke(t₂) = exp(2 ν |K|² (t₂ - t₁))
///   ⇒  ν = ln(ke(t₁) / ke(t₂)) / (2 |K|² Δt)
///
/// For the standard TGV mode taken with the single wavenumber `k` on each
/// axis in D dimensions, `|K|² = D · k²`:
///
///   D = 2:  ν = ln(ke₁/ke₂) / (4 k² Δt)     (matches validation_tgv.rs)
///   D = 3:  ν = ln(ke₁/ke₂) / (6 k² Δt)     (the QA-sweep W-LES freeze form)
///
/// Both cases share the same code path — pass the correct `k2_sum = D · k²`
/// (or the full |K|² of your mode set for multi-mode inits).
///
/// `dt` is the number of steps between the two energy samples in lattice
/// units; energies must come from `ke2d`/`ke3d` on the SAME field size (the
/// ratio cancels the cells-count factor, but only if the field size is
/// identical between samples).
pub fn tgv_nu_eff(ke1: f64, ke2: f64, k2_sum: f64, dt: f64) -> f64 {
    (ke1 / ke2).ln() / (2.0 * k2_sum * dt)
}

/// Exponential decay rate from two energy samples: `-ln(ke₂/ke₁) / Δt`.
/// This is the raw observable that `tgv_nu_eff` normalizes by `|K|²`.
/// Present here so callers that want the rate for other purposes
/// (LES-vs-DNS characterization tables, multi-run decay comparisons) do not
/// have to invert `tgv_nu_eff` back to the rate.
pub fn energy_decay_rate(ke1: f64, ke2: f64, dt: f64) -> f64 {
    -(ke2 / ke1).ln() / dt
}
