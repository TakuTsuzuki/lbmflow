# Structural sweep — lanes 7.2 / 7.4 / 2.3 (2026-07-06, PM-triaged)

Read-only subagent sweep + PM verification. Raw agent output triaged below.

## 7.2 Numeric-literal provenance (physics paths)

- VERY_HIGH flags (kernels.rs 0.0025 / 0.16): already adjudicated —
  ANOM-P4-008 verdict (C), core removal queued; no separate action.
- Actionable P2 comment queue (fold into ONE doc-hygiene codex order,
  together with the ANOM-P4-009 rotor contract comment and the
  ANOM-P4-011 q-convention note):
  - kernels.rs ~282: 1e-14 Gaussian-elimination pivot floor — mark as
    numerical guard with rationale.
  - les.rs 104-106: WALE exponents 2.5/1.5/1.25 — cite Nicoud-Ducros
    operator definition inline.
  - rotating_ibm.rs 218-220: quadratic B-spline kernel coefficients —
    cite the kernel family.
  - particles.rs 205: SN 0.15/Re^0.687 — inline citation (validity note
    already present).
  - solver.rs / compat: rho.max(1e-30) underflow guards — one-line guard
    comments.
- No NEW undocumented physics constants found beyond the known ANOM set.

## 7.4 VALIDATION.md ↔ test drift

- 14/14 verifiable doc bands match the asserted bands EXACTLY (T1, T2, T3,
  T5, T6, T7×3, T11 spurious currents, TGV symmetry, and others). Zero
  numeric drift.
- The agent's "T4/T5 MISSING" rows were AGENT ERRORS (looked in
  validation_channel.rs; the tests live in validation_open_bc.rs:
  t4_velocity_inlet_pressure_outlet_channel_all_four_orientations,
  t5_pressure_pressure_channel_flow_rate_and_linear_pressure,
  t5_pressure_reversal_with_x_mirror_is_exact — PM-verified). No missing
  tests confirmed by this sweep.
- T7 vortex-center band: assertion exists (validation_cavity.rs) — the
  agent marked it "implicit"; verified present during the lane-2.2 scan.

## 2.3 Behavior-anchor coverage (two-layer rule)

- 17/23 T-sections carry explicit pattern/sign/monotonicity/symmetry
  anchors alongside scalar bands. Distribution: symmetry 7, pattern 8,
  monotonicity 5, sign 3 — healthy mix.
- Band-only: T4 (justified — the anchored artifact IS the documented
  outlet oscillation spec). Partial: T13 diagnostics, T16/T17
  (post-implementation freeze pending, consistent with the traceability
  matrix SPEC-ONLY rows). NO retrofit order needed now; T17 anchors are
  mandated by REQ rev.4 negative tests when implemented.

## Verdict

Structure is sound: zero spec/test numeric drift, anchor coverage healthy,
literal hygiene reduces to one P2 comment order. Lanes 7.2/7.4/2.3 DONE.
