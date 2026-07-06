//! Agreement metrics for accuracy-audit and validation tests.
//!
//! Pure functions only: no simulation types, no I/O. Every accuracy assertion
//! in an adversarial test should go through one of these so the semantics are
//! reviewable in one place. Mirrored 1:1 (names + semantics) by
//! `scripts/qa/metrics.py` for post-hoc analysis of CLI output.

/// Relative L2 norm of the error between `actual` and `reference`:
/// `||a - r||_2 / ||r||_2`.
pub fn l2_rel(actual: &[f64], reference: &[f64]) -> f64 {
    assert_eq!(actual.len(), reference.len());
    let mut num = 0.0;
    let mut den = 0.0;
    for (a, r) in actual.iter().zip(reference) {
        num += (a - r) * (a - r);
        den += r * r;
    }
    (num / den).sqrt()
}

/// Relative L-infinity error: `max |a - r| / max(max |r|, floor)`.
/// `floor` guards near-zero reference fields (pass 0.0 to disable).
pub fn linf_rel(actual: &[f64], reference: &[f64], floor: f64) -> f64 {
    assert_eq!(actual.len(), reference.len());
    let mut dmax = 0.0f64;
    let mut rmax = floor;
    for (a, r) in actual.iter().zip(reference) {
        dmax = dmax.max((a - r).abs());
        rmax = rmax.max(r.abs());
    }
    dmax / rmax
}

/// Least-squares straight-line fit `y ≈ slope·x + intercept`, with the
/// coefficient of determination `r2` so callers can reject sloppy fits
/// (an order estimate from a bad fit is meaningless — assert `r2` too).
#[derive(Clone, Copy, Debug)]
pub struct LinFit {
    pub slope: f64,
    pub intercept: f64,
    pub r2: f64,
}

pub fn linear_fit(x: &[f64], y: &[f64]) -> LinFit {
    assert_eq!(x.len(), y.len());
    assert!(x.len() >= 2, "linear_fit needs >= 2 points");
    let n = x.len() as f64;
    let mx = x.iter().sum::<f64>() / n;
    let my = y.iter().sum::<f64>() / n;
    let mut sxx = 0.0;
    let mut sxy = 0.0;
    let mut syy = 0.0;
    for (&xi, &yi) in x.iter().zip(y) {
        sxx += (xi - mx) * (xi - mx);
        sxy += (xi - mx) * (yi - my);
        syy += (yi - my) * (yi - my);
    }
    let slope = sxy / sxx;
    let intercept = my - slope * mx;
    let r2 = if syy == 0.0 { 1.0 } else { (sxy * sxy) / (sxx * syy) };
    LinFit {
        slope,
        intercept,
        r2,
    }
}

/// Observed convergence order from resolutions `h` and errors `err`
/// (log-log fit; slope = order). Assert BOTH `fit.slope` against the
/// theoretical order band AND `fit.r2` (>= ~0.98) — points off a straight
/// log-log line mean the asymptotic regime was not reached.
pub fn order_fit(h: &[f64], err: &[f64]) -> LinFit {
    let lx: Vec<f64> = h.iter().map(|v| v.ln()).collect();
    let ly: Vec<f64> = err.iter().map(|v| v.ln()).collect();
    linear_fit(&lx, &ly)
}

/// Exponential-envelope fit `amp ≈ A·exp(-k·y)` (semilog fit of ln(amp) on y).
/// Returns the fit on ln(amp): amplitude `A = intercept.exp()`, decay
/// `k = -slope`. Used for e.g. Stokes-II wall-layer decay and acoustic
/// damping envelopes.
pub fn envelope_fit(y: &[f64], amp: &[f64]) -> LinFit {
    let ly: Vec<f64> = amp.iter().map(|v| v.ln()).collect();
    linear_fit(y, &ly)
}

/// Single-frequency quadrature projection of `signal(t)` on `sin/cos(omega t)`:
/// returns `(amplitude, phase)` with `signal ≈ amplitude·sin(omega·t + phase)`.
/// Sampling should cover an integer number of periods for clean orthogonality.
pub fn phase_fit(t: &[f64], signal: &[f64], omega: f64) -> (f64, f64) {
    assert_eq!(t.len(), signal.len());
    assert!(!t.is_empty());
    let n = t.len() as f64;
    let mut s = 0.0;
    let mut c = 0.0;
    for (&ti, &si) in t.iter().zip(signal) {
        s += si * (omega * ti).sin();
        c += si * (omega * ti).cos();
    }
    let a_sin = 2.0 * s / n; // coefficient of sin(omega t)
    let a_cos = 2.0 * c / n; // coefficient of cos(omega t)
    ((a_sin * a_sin + a_cos * a_cos).sqrt(), a_cos.atan2(a_sin))
}

/// Fraction of adjacent pairs that are strictly decreasing (1.0 = the whole
/// sequence decreases monotonically). Use for error-vs-resolution and
/// error-vs-time sequences where theory demands monotone decay.
pub fn monotonicity(xs: &[f64]) -> f64 {
    assert!(xs.len() >= 2, "monotonicity needs >= 2 points");
    let dec = xs.windows(2).filter(|w| w[1] < w[0]).count();
    dec as f64 / (xs.len() - 1) as f64
}

/// Result of comparing measured samples against a theoretical curve.
#[derive(Clone, Copy, Debug)]
pub struct CurveAgreement {
    /// Largest relative deviation |measured - theory| / max(|theory|, floor).
    pub max_rel_dev: f64,
    /// x at which the largest deviation occurs.
    pub worst_x: f64,
    /// Fraction of samples within `rel_band` of the curve.
    pub frac_in_band: f64,
}

/// A3-axis primitive: the measured error/observable must lie ON the known
/// theoretical curve, point by point — "small" is not a pass. `floor` guards
/// zero crossings of the theory curve (pass 0.0 to disable).
pub fn curve_agreement(
    theory: impl Fn(f64) -> f64,
    samples: &[(f64, f64)],
    rel_band: f64,
    floor: f64,
) -> CurveAgreement {
    assert!(!samples.is_empty());
    let mut max_rel_dev = 0.0f64;
    let mut worst_x = samples[0].0;
    let mut in_band = 0usize;
    for &(x, measured) in samples {
        let th = theory(x);
        let dev = (measured - th).abs() / th.abs().max(floor);
        if dev > max_rel_dev {
            max_rel_dev = dev;
            worst_x = x;
        }
        if dev <= rel_band {
            in_band += 1;
        }
    }
    CurveAgreement {
        max_rel_dev,
        worst_x,
        frac_in_band: in_band as f64 / samples.len() as f64,
    }
}
