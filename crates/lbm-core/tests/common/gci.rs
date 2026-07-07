//! Roache Grid Convergence Index helpers for three-grid solution verification.
//!
//! Notation follows the usual Roache convention: `f_fine = f_h`,
//! `f_medium = f_2h`, `f_coarse = f_4h`. The refinement ratios are
//! `r_21 = h_medium / h_fine` and `r_32 = h_coarse / h_medium`.

/// Roache's recommended safety factor for a three-grid family in the
/// asymptotic range.
pub const ASYMPTOTIC_THREE_GRID_SAFETY_FACTOR: f64 = 1.25;

/// Complete three-grid GCI summary for one quantity of interest.
#[derive(Clone, Copy, Debug)]
pub struct GciResult {
    /// Observed order of accuracy `p`.
    pub observed_order: f64,
    /// Fine-grid Richardson estimate of the zero-grid-spacing limit.
    pub richardson_limit: f64,
    /// Fine-grid GCI, `GCI_21`, reported as a percent.
    pub gci_21_pct: f64,
    /// Medium-grid GCI, `GCI_32`, reported as a percent.
    pub gci_32_pct: f64,
    /// Roache asymptotic-range consistency ratio; values near 1 indicate
    /// that the three grids are in the asymptotic range.
    pub asymptotic_ratio: f64,
}

fn assert_valid_ratio(name: &str, r: f64) {
    assert!(
        r.is_finite() && r > 1.0,
        "{name} must be finite and > 1, got {r:e}"
    );
}

fn assert_finite(name: &str, value: f64) {
    assert!(value.is_finite(), "{name} must be finite, got {value:e}");
}

/// Observed order of accuracy from three grids.
///
/// With equal refinement ratio `r`, the assumed expansion is
///
/// ```text
/// f_h = f_0 + C h^p
/// ```
///
/// so
///
/// ```text
/// f_4h - f_2h = C h^p (r_32^p r_21^p - r_21^p)
/// f_2h - f_h  = C h^p (r_21^p - 1)
/// ```
///
/// and for `r_21 == r_32 == r` this reduces to Roache's direct form
/// `p = ln((f_4h - f_2h) / (f_2h - f_h)) / ln(r)`. For unequal ratios this
/// uses the standard fixed-point form
///
/// ```text
/// p = [ln(|e32 / e21|) + ln((r_21^p - s) / (r_32^p - s))] / ln(r_21)
/// ```
///
/// where `e32 = f_4h - f_2h`, `e21 = f_2h - f_h`, and
/// `s = sign(e32 / e21)`.
pub fn observed_order(f_coarse: f64, f_medium: f64, f_fine: f64, r_21: f64, r_32: f64) -> f64 {
    assert_valid_ratio("r_21", r_21);
    assert_valid_ratio("r_32", r_32);
    let e32 = f_coarse - f_medium;
    let e21 = f_medium - f_fine;
    assert_finite("e32", e32);
    assert_finite("e21", e21);
    assert!(
        e32 != 0.0 && e21 != 0.0,
        "observed_order needs nonzero grid differences, got e32={e32:e}, e21={e21:e}"
    );

    let ratio = e32 / e21;
    if (r_21 - r_32).abs() <= 1.0e-12 * r_21.max(r_32) {
        assert!(
            ratio > 0.0,
            "equal-ratio observed_order needs monotone convergence, got e32/e21={ratio:e}"
        );
        return ratio.ln() / r_21.ln();
    }

    let s = ratio.signum();
    let mut p = ratio.abs().ln() / r_21.ln();
    assert!(
        p.is_finite() && p > 0.0,
        "initial observed_order iterate must be finite and positive, got {p:e}"
    );
    for _ in 0..50 {
        let numerator = r_21.powf(p) - s;
        let denominator = r_32.powf(p) - s;
        assert!(
            numerator > 0.0 && denominator > 0.0,
            "observed_order fixed-point term invalid: p={p:e}, numerator={numerator:e}, denominator={denominator:e}"
        );
        let next = (ratio.abs().ln() + (numerator / denominator).ln()) / r_21.ln();
        assert_finite("observed_order iterate", next);
        if (next - p).abs() <= 1.0e-12 {
            return next;
        }
        p = next;
    }
    p
}

/// Richardson extrapolation from the fine and medium grids.
///
/// From `f_h = f_0 + C h^p` and `f_2h = f_0 + C (r_21 h)^p`,
/// eliminating `C h^p` gives
///
/// ```text
/// f_0 ~= f_h + (f_h - f_2h) / (r_21^p - 1)
/// ```
pub fn richardson_extrapolate(f_medium: f64, f_fine: f64, r_21: f64, p: f64) -> f64 {
    assert_valid_ratio("r_21", r_21);
    assert_finite("p", p);
    let denom = r_21.powf(p) - 1.0;
    assert!(
        denom > 0.0,
        "Richardson denominator must be positive, got {denom:e} for r_21={r_21:e}, p={p:e}"
    );
    f_fine + (f_fine - f_medium) / denom
}

/// Roache Grid Convergence Index for one adjacent grid pair.
///
/// The fractional GCI is
///
/// ```text
/// GCI_21 = Fs * |f_h - f_2h| / (|f_h| * (r_21^p - 1))
/// ```
///
/// with safety factor `Fs = 1.25` for a three-grid family in the asymptotic
/// range. Use `Fs = 3.0` for a two-grid estimate or when the asymptotic range
/// has not been demonstrated.
pub fn gci(f_medium: f64, f_fine: f64, r_21: f64, p: f64, safety_factor: f64) -> f64 {
    assert_valid_ratio("r_21", r_21);
    assert_finite("p", p);
    assert_finite("safety_factor", safety_factor);
    assert!(
        safety_factor > 0.0,
        "safety_factor must be positive, got {safety_factor:e}"
    );
    assert!(
        f_fine != 0.0,
        "GCI fractional normalization needs nonzero fine-grid value"
    );
    let denom = r_21.powf(p) - 1.0;
    assert!(
        denom > 0.0,
        "GCI denominator must be positive, got {denom:e} for r_21={r_21:e}, p={p:e}"
    );
    safety_factor * (f_fine - f_medium).abs() / (f_fine.abs() * denom)
}

/// Roache asymptotic-range consistency check.
///
/// In the asymptotic range, adjacent-grid GCI values should shrink by
/// approximately `r_21^p`, so
///
/// ```text
/// GCI_32 / (r_21^p * GCI_21) ~= 1
/// ```
pub fn asymptotic_range_check(gci_21: f64, gci_32: f64, r_21: f64, p: f64) -> f64 {
    assert_valid_ratio("r_21", r_21);
    assert_finite("p", p);
    assert!(
        gci_21 > 0.0 && gci_32 > 0.0,
        "asymptotic_range_check needs positive GCI values, got GCI_21={gci_21:e}, GCI_32={gci_32:e}"
    );
    gci_32 / (r_21.powf(p) * gci_21)
}

/// Build a complete GCI summary from coarse, medium, and fine values.
pub fn gci_result(
    f_coarse: f64,
    f_medium: f64,
    f_fine: f64,
    r_21: f64,
    r_32: f64,
    safety_factor: f64,
) -> GciResult {
    let p = observed_order(f_coarse, f_medium, f_fine, r_21, r_32);
    let richardson_limit = richardson_extrapolate(f_medium, f_fine, r_21, p);
    let gci_21 = gci(f_medium, f_fine, r_21, p, safety_factor);
    let gci_32 = gci(f_coarse, f_medium, r_32, p, safety_factor);
    let asymptotic_ratio = asymptotic_range_check(gci_21, gci_32, r_21, p);
    GciResult {
        observed_order: p,
        richardson_limit,
        gci_21_pct: 100.0 * gci_21,
        gci_32_pct: 100.0 * gci_32,
        asymptotic_ratio,
    }
}

/// Convenience wrapper for an equal-ratio series ordered coarse, medium, fine.
///
/// Example: for grids `D = 20, 40, 80`, call
/// `gci_from_series([qoi_d20, qoi_d40, qoi_d80], 2.0)`.
pub fn gci_from_series(coarse_medium_fine: [f64; 3], r: f64) -> GciResult {
    gci_result(
        coarse_medium_fine[0],
        coarse_medium_fine[1],
        coarse_medium_fine[2],
        r,
        r,
        ASYMPTOTIC_THREE_GRID_SAFETY_FACTOR,
    )
}
