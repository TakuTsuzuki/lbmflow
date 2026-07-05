# Shared agreement-metrics library — API

Two mirrored implementations, kept in sync in the same commit:

- **Rust (source of truth, in-suite assertions)**:
  `crates/lbm-core/tests/common/metrics.rs` — import in a test file with
  `mod common;` then `use common::metrics::*;`
- **Python (post-hoc analysis of CLI/scenario output, stdlib-only)**:
  `scripts/qa/metrics.py` — composes with the output parsers in
  `scripts/qa/qa_checks.py` (branch qa/anomaly-sweep). Self-test:
  `python3 scripts/qa/metrics.py`.

All functions are pure (no simulation types, no I/O). Self-test fixtures are
identical on both sides (`crates/lbm-core/tests/metrics_selftest.rs` /
`_selftest()`).

## Error norms

```rust
l2_rel(actual: &[f64], reference: &[f64]) -> f64
```
`||a − r||₂ / ||r||₂`. The default profile-agreement metric.

```rust
linf_rel(actual: &[f64], reference: &[f64], floor: f64) -> f64
```
`max|a − r| / max(max|r|, floor)`. Use for "no point escapes the band";
`floor` guards near-zero reference fields (0.0 disables).

## Fits (all return `LinFit { slope, intercept, r2 }`)

```rust
linear_fit(x: &[f64], y: &[f64]) -> LinFit
```
Least-squares line. **Always assert `r2` alongside the slope** — a slope
through scattered points is noise, and `order_fit` inherits this rule.

```rust
order_fit(h: &[f64], err: &[f64]) -> LinFit      // A1
```
Log-log fit; `slope` = observed convergence order. Assert the slope band
(e.g. `(1.8..=2.3).contains(&fit.slope)`) AND `fit.r2 >= 0.98`. Plateaued
(non-asymptotic) data fails the r2 gate by construction.

```rust
envelope_fit(y: &[f64], amp: &[f64]) -> LinFit   // A4
```
Semilog fit `amp ≈ A·exp(−k·y)`: `A = fit.intercept.exp()`, `k = −fit.slope`.
For Stokes-II wall layers, acoustic damping envelopes.

## Oscillations

```rust
phase_fit(t: &[f64], signal: &[f64], omega: f64) -> (f64, f64)  // A4
```
Quadrature projection at a known frequency → `(amplitude, phase)` with
`signal ≈ amplitude·sin(omega·t + phase)`. Sample an integer number of
periods (orthogonality). Frequency checks (e.g. sound speed cs = 1/√3) work
by projecting at the *theoretical* omega: a detuned actual frequency leaks
energy out of the projection, so asserting the fitted amplitude recovers the
expected signal amplitude doubles as the frequency assertion.

## Sequences

```rust
monotonicity(xs: &[f64]) -> f64
```
Fraction of adjacent pairs strictly decreasing; 1.0 = monotone decay. For
error-vs-resolution and error-vs-time sequences where theory demands decay.

## Curves (A3 primitive)

```rust
curve_agreement(theory: impl Fn(f64) -> f64,
                samples: &[(f64, f64)],
                rel_band: f64, floor: f64) -> CurveAgreement
// CurveAgreement { max_rel_dev, worst_x, frac_in_band }
```
Point-by-point relative deviation of measured samples from a theoretical
curve. Assert `max_rel_dev <= band` (strict) or `frac_in_band >= 0.9`
(tolerant tails), and print `worst_x` — it localizes the finding. This is
the "error must lie ON the curve; small is not a pass" instrument.

## Extension rule

A new metric goes into BOTH files with identical fixtures in both
self-tests, in the same commit. A metric needed by only one test does not
belong here until a second caller appears — inline it in the test with its
derivation until then.
