# S-Fingerprint spec (pre-authoring, PM-fixed gates — SCALEUP v1.1 + review corrections)

Status: SPEC ONLY. SU-fingerprint-author spawns after the ε/shear output channel
lands (reactor-demo order C item 2). Scope = GREEN/YELLOW per s0-capability-map.md:
{U_tip, Re_impeller, Fr, V, ε-field stats (mean/max/percentiles), η} with
Np / P / power-based ε̄ / θ95 / lifelines / t_c routed to
unsupported_or_spec_only_items. Input contract: N, D, T, H, g, ρ, μ are USER
inputs (the schema has no impeller); the ONLY sim-derived quantities are the
resolved ε field + η.

## Hard verification gates (replace the three original gates, which all
## depended on RED primitives — torque ×2, tracer ADE ×1)

- **G1 (formula recompute, analytic exact)**: U_tip = πND, Re = ND²/ν, Fr = N²D/g,
  V, η = (ν³/ε)^¼ recomputed from the declared inputs must match the emitted
  fingerprint exactly (dimensional/formula-correctness gate).
- **G2 (analytic-strain gate — LOAD-BEARING)**: validate ε-field extraction against
  exact solutions, not power:
  - Plane Couette (T3 setup): du/dy = U/H uniform ⇒ ε = 2νS:S = ν(U/H)² uniform.
    Pointwise L∞rel ≤ ~1e-2 (interior) + volume mean within a frozen band
    (characterize→freeze; strain gathers are machine-precision so expect tight).
  - Body-force Poiseuille (T2): ε(y) = ν·[(g/2ν)(H−2y_w)]²; volume mean matches
    the analytic integral (band frozen at characterization).
  Rationale: with N, D user-declared, G2 is the ONLY sim-derived verification in
  this pilot scope — the Skill must not ship without it. It also pre-validates the
  ε machinery for the day the M-F impeller lands.

## Estimator disambiguation (Minor 1)

Power-based ε̄ = P/(ρV): spec-only (torque RED). Volume-mean of the RESOLVED field
<ε>_vol: GREEN once the channel lands — emit it, labeled distinctly in the JSON
(e.g. eps_resolved_volmean vs eps_power_speconly). The old "resolved vs P/(ρV)
±15%" cross-check gate is deferred until torque lands.

## S-Match Gate-A precondition (Minor 2 — carried into match-specs)

The P/V and ε̄ matchers (N_low = N_high·S^(2/3)) assume Np ≈ const, valid only
fully turbulent (Re_impeller ≳ 1e4). Gate-A must assert Re_impeller ≳ 1e4 for
those two criteria; below it, flag into residuals_to_report instead of returning
a confident N_low (laminar: Np·Re ≈ const breaks the scaling silently).
Tip-speed / Re / Fr matchers are unconditional.
