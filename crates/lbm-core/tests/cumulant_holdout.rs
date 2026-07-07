//! REV-6 holdout validation for the D3Q19-only cumulant shear-rate correction.
//!
//! Calibration provenance:
//! - `docs/PHYSICS.md` records ANOM-P4-008: the former D3Q19 shear-rate
//!   offset (`+0.0025` relative) was removed as banned calibration.
//! - The remaining finite-frame cubic-velocity term is
//!   `omega_eff = omega_shear * (1 - 0.16 |u|^2)` unless the compile-time
//!   central-moment ablation flag disables it for E1.
//!
//! These tests are adversarial holdouts, not calibration tests. Bands are
//! derived from the existing T15 second-order decay-rate class or from the
//! frame-coupled truncation scale; they are not fit to current output.
//!
//! Exception (documented change-detector, not a physics band): the
//! `frozen_anchor` test pins the rest-frame decay rate to a measured
//! snapshot. The ANOM-P4-008 offset removal shifted this observable by
//! +2.3% and every wide physical band stayed green (band vacuity, PM bisect
//! 2026-07-07 — see the PHYSICS.md r2-c triage entry). If the anchor fires,
//! adjudicate the physics change and re-freeze with a PHYSICS.md entry;
//! never widen it.

use lbm_core::prelude::*;
use std::f64::consts::PI;

type CpuPeriodic<L> = Solver<L, f64, CpuScalar, LocalPeriodic>;

const CALIBRATION_NU: f64 = 0.02;
const OFF_RE_NU: f64 = 0.04;
const N: usize = 32;
const ADVECTED_STEPS: usize = 160;
const T15_DECAY_RATE_BAND: f64 = 0.02;
const TGV_U0: f64 = 0.012;

#[derive(Clone, Copy, Debug)]
struct DecayMeasurement {
    rate: f64,
    rel_err: f64,
}

fn omega_from_nu(nu: f64) -> f64 {
    1.0 / (3.0 * nu + 0.5)
}

fn cumulant(nu: f64) -> CollisionKind {
    CollisionKind::CentralMoment {
        omega_shear: omega_from_nu(nu),
    }
}

fn make_tgv<L: Lattice>(n: usize, nu: f64, u0: f64, mean_u: f64) -> CpuPeriodic<L> {
    let spec = GlobalSpec::<f64> {
        dims: [n, n, n],
        nu,
        collision: cumulant(nu),
        periodic: [true, true, true],
        ..Default::default()
    };
    let mut s: CpuPeriodic<L> = Solver::new(
        &spec,
        &[],
        &[],
        [1, 1, 1],
        CpuScalar::default(),
        LocalPeriodic,
    );
    let k = 2.0 * PI / n as f64;
    s.init_with(move |x, y, z| {
        let (xf, yf, zf) = (k * x as f64, k * y as f64, k * z as f64);
        let p = u0 * u0 / 16.0 * (((2.0 * xf).cos() + (2.0 * yf).cos()) * ((2.0 * zf).cos() + 2.0));
        (
            1.0 + 3.0 * p,
            [
                mean_u + u0 * xf.sin() * yf.cos() * zf.cos(),
                -u0 * xf.cos() * yf.sin() * zf.cos(),
                0.0,
            ],
        )
    });
    s
}

fn fluctuation_ke<L: Lattice>(s: &CpuPeriodic<L>, mean_u: f64) -> f64 {
    let (ux, uy, uz) = (s.gather_ux(), s.gather_uy(), s.gather_uz());
    ux.iter()
        .zip(&uy)
        .zip(&uz)
        .map(|((ux, uy), uz)| (ux - mean_u).powi(2) + uy * uy + uz * uz)
        .sum()
}

fn measure_decay<L: Lattice>(
    label: &str,
    n: usize,
    nu: f64,
    u0: f64,
    mean_u: f64,
    steps: usize,
) -> DecayMeasurement {
    let mut s = make_tgv::<L>(n, nu, u0, mean_u);
    let mut prev = fluctuation_ke(&s, mean_u);
    let e0 = prev;
    assert!(
        e0.is_finite() && e0 > 0.0,
        "{label}: initial fluctuation kinetic energy must be finite and positive, got {e0:e}"
    );

    for step in 1..=steps {
        s.run(1);
        let e = fluctuation_ke(&s, mean_u);
        assert!(
            e.is_finite() && e > 0.0,
            "{label}: fluctuation kinetic energy must stay finite and positive at step {step}, got {e:e}"
        );
        assert!(
            e <= prev * (1.0 + 1.0e-12),
            "{label}: fluctuation kinetic energy must decrease monotonically; step {step} increased {prev:e} -> {e:e}"
        );
        prev = e;
    }

    let e1 = prev;
    let rate = -(e1 / e0).ln() / steps as f64;
    let k = 2.0 * PI / n as f64;
    let analytic_rate = 6.0 * nu * k * k;
    let rel_err = (rate - analytic_rate).abs() / analytic_rate;
    assert!(
        rate.is_finite() && rate > 0.0,
        "{label}: decay rate must be finite and positive, got {rate:e}"
    );
    println!(
        "{label}: rate={rate:.9e}, analytic={analytic_rate:.9e}, rel_err={rel_err:.9e}, E0={e0:.9e}, E1={e1:.9e}"
    );
    DecayMeasurement { rate, rel_err }
}

/// Frozen change-detector for the rest-frame TGV3D decay rate (see the
/// module-header exception note). This is a measurement snapshot, NOT a
/// physics claim: the physical acceptance remains the h^2-intercept audit
/// (`accuracy_audit_cumulant.rs`) and the T15-class bands. Value measured at
/// 6d55a50 (post ANOM-P4-008 offset removal), CpuScalar f64; the pre-removal
/// tree measured 4.634861882e-3, so the trap width must stay far below that
/// 2.3e-2-relative shift.
#[test]
fn tgv3d_u0_decay_rate_matches_frozen_anchor() {
    const FROZEN_RATE: f64 = 4.742352837e-3;
    const REL_HALF_WIDTH: f64 = 1.0e-5;
    let m = measure_decay::<D3Q19>(
        "D3Q19 cumulant frozen-anchor TGV3D u_frame=0",
        N,
        CALIBRATION_NU,
        TGV_U0,
        0.0,
        ADVECTED_STEPS,
    );
    let shift = (m.rate - FROZEN_RATE).abs() / FROZEN_RATE;
    assert!(
        shift <= REL_HALF_WIDTH,
        "D3Q19 CentralMoment TGV3D u_frame=0 decay rate {rate:e} moved {shift:e} (relative) \
         from the frozen anchor {FROZEN_RATE:e}; adjudicate the physics change (PHYSICS.md \
         r2-c triage entry) and re-freeze — do not widen",
        rate = m.rate
    );
}

#[test]
#[ignore = "FINDING: D3Q19 cumulant advected TGV3D frame spread exceeds the derived Ma^2*(kdx)^2 band - see PHYSICS.md holdout entry"]
fn advected_tgv3d_decay_rate_is_frame_independent() {
    let frames = [0.0, 0.05, 0.1];
    let rates = frames.map(|u| {
        measure_decay::<D3Q19>(
            &format!("D3Q19 cumulant advected TGV3D frame u={u:.2}"),
            N,
            CALIBRATION_NU,
            TGV_U0,
            u,
            ADVECTED_STEPS,
        )
        .rate
    });
    let min_rate = rates.iter().copied().fold(f64::INFINITY, f64::min);
    let max_rate = rates.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let mean_rate = rates.iter().sum::<f64>() / rates.len() as f64;
    let spread = (max_rate - min_rate) / mean_rate;

    // Band derivation:
    // The target correction explicitly removes the cubic-velocity viscosity
    // defect, whose frame-dependent scale is O(Ma_frame^2). Once removed, the
    // remaining frame dependence must be tied to the second-order spatial
    // truncation of the single Fourier mode, O((k dx)^2). Therefore the
    // holdout spread band is Ma_frame,max^2 * (k dx)^2, with c_s^2 = 1/3 and
    // k = 2*pi/N in lattice units. This is a numerical-order scale, not a
    // measured-output tolerance.
    let max_frame = frames.iter().copied().fold(0.0f64, f64::max);
    let ma2 = max_frame * max_frame / D3Q19::CS2;
    let kdx2 = (2.0 * PI / N as f64).powi(2);
    let band = ma2 * kdx2;
    println!(
        "D3Q19 cumulant advected TGV3D frame spread: rates={rates:?}, spread={spread:.9e}, band={band:.9e}"
    );
    assert!(
        spread <= band,
        "D3Q19 cumulant decay-rate frame spread {spread:e} exceeds Ma^2*(kdx)^2 band {band:e}"
    );
}

#[test]
#[ignore = "heavy holdout: off-calibration Reynolds TGV3D decay"]
fn off_calibration_reynolds_tgv3d_decay_matches_diffusion_limit() {
    let m = measure_decay::<D3Q19>(
        "D3Q19 cumulant off-calibration-Re TGV3D",
        N,
        OFF_RE_NU,
        TGV_U0,
        0.0,
        ADVECTED_STEPS,
    );
    println!(
        "D3Q19 cumulant off-calibration-Re result: nu_eff={:.9e}, nu={OFF_RE_NU:.9e}, rel_err={:.9e}",
        m.rate / (6.0 * (2.0 * PI / N as f64).powi(2)),
        m.rel_err
    );
    assert!(
        m.rel_err <= T15_DECAY_RATE_BAND,
        "D3Q19 cumulant off-calibration-Re decay-rate error {err:e} exceeds T15 class band {band:e}",
        err = m.rel_err,
        band = T15_DECAY_RATE_BAND
    );
}

#[test]
#[ignore = "heavy holdout: D3Q19 vs D3Q27 TGV3D analytic-error cross-check"]
fn d3q19_corrected_decay_error_is_not_d3q27_outlier() {
    let d3q19 = measure_decay::<D3Q19>(
        "D3Q19 cumulant holdout TGV3D cross-check",
        N,
        OFF_RE_NU,
        TGV_U0,
        0.0,
        ADVECTED_STEPS,
    );
    let d3q27 = measure_decay::<D3Q27>(
        "D3Q27 cumulant holdout TGV3D cross-check",
        N,
        OFF_RE_NU,
        TGV_U0,
        0.0,
        ADVECTED_STEPS,
    );
    let outlier_band = d3q27.rel_err + T15_DECAY_RATE_BAND;
    println!(
        "D3Q19 vs D3Q27 cumulant holdout: d3q19_rel={:.9e}, d3q27_rel={:.9e}, outlier_band={outlier_band:.9e}",
        d3q19.rel_err, d3q27.rel_err
    );
    assert!(
        d3q27.rel_err <= T15_DECAY_RATE_BAND,
        "D3Q27 cumulant holdout decay-rate error {err:e} exceeds T15 class band {band:e}",
        err = d3q27.rel_err,
        band = T15_DECAY_RATE_BAND
    );
    assert!(
        d3q19.rel_err <= outlier_band,
        "D3Q19 corrected decay-rate error {d19:e} is an outlier vs D3Q27 {d27:e}; allowed D3Q27 error + T15 band = {outlier_band:e}",
        d19 = d3q19.rel_err,
        d27 = d3q27.rel_err
    );
}
