# Credibility Policy — Bioprocess CFD

Lifecycle: living (owning doc for calibration / holdout / UQ / evidence-tier
policy).

Companion to [SPEC_BIOPROCESS_CORE.md](SPEC_BIOPROCESS_CORE.md) tiers and
[MODEL_RISK_MATRIX.md](MODEL_RISK_MATRIX.md). This document defines the
process gates a QOI must pass before it can carry an Engineering or
Evidence label in a `lbm report` output.

## 1. Credibility tiers (recap)

- **Tier 0 – Screening.** Qualitative comparison across configurations. Runs
  may be under-resolved. Never used as evidence artefact.
- **Tier 1 – Engineering.** Physics model has its VB validation entry
  green. Mesh / time-step resolved to documented convergence band. QOI
  emitted with full provenance. No calibration / holdout pinned. Suitable
  for design-of-experiments and internal decisions.
- **Tier 2 – Evidence.** Engineering + calibration dataset + independent
  holdout dataset + UQ interval on the QOI + mesh / time-step sensitivity
  records; enforced by BCFD-091 EvidenceGate. Suitable for external
  submissions where the tool's role must be defensible.

## 2. Calibration and holdout separation (BCFD-082)

A QOI is *calibrated* when at least one model parameter (kL model, drag
coefficient, breakup / coalescence kernel constants, wall function
constants, ...) is fitted to measurement data. A QOI is *validated* when
the model — including any calibrated parameters — reproduces measurements
that were **not seen during calibration**.

Rules:

- Every dataset used has an `id`, the QOI it reports, source (paper /
  internal / vendor), date, scale, operating condition, and measurement
  uncertainty (BCFD-082).
- **The same dataset id cannot appear in both the calibration and holdout
  registry.** BCFD-082 rejects the scenario at parse time.
- Evidence tier requires at least one holdout dataset per QOI being
  claimed.
- Engineering tier allows calibration-only, but the report label surfaces
  "calibrated to <ids>, not validated against holdout".
- Screening tier ignores the registry.
- When a holdout dataset is subsequently used to *retune* a parameter, it
  becomes a calibration dataset and can no longer serve as a holdout for
  that QOI — even if you register it with a new id.

## 3. UQ requirements (BCFD-083)

Evidence tier requires an uncertainty interval on each reported QOI. In
priority order:

1. **Model-form uncertainty.** Repeat the QOI with the alternative closure
   set enumerated in [MODEL_RISK_MATRIX.md](MODEL_RISK_MATRIX.md); the
   interval spans the plausible outcomes. This is the mandatory minimum.
2. **Parameter uncertainty.** Sweep BCFD-083 across the calibrated
   parameter's confidence interval; propagate through to the QOI.
3. **Numerical uncertainty.** Refined mesh + halved Δt as the reference;
   convergence exponent recorded (three-grid GCI when meshes admit it).

The report cites the interval as `Q = q̂ [q_lo, q_hi] (method: ...)`,
never as a bare number.

## 4. Sensitivity records (BCFD-083)

Evidence tier also requires:

- **Mesh sensitivity.** At least two additional grid resolutions; the QOI
  variation between the two finest is ≤5% (or the value used as the QOI
  uncertainty floor).
- **Time-step sensitivity.** At least two Δt values; QOI variation ≤5%
  between the finer two.
- **Windowing sensitivity.** For time-averaged QOIs, the averaging window
  is documented and the QOI is stable to ≥50% window changes.

Failure of any of these blocks evidence tier but not engineering tier.

## 5. Evidence gate (BCFD-091)

`EvidenceGate` returns one of:

- `EvidenceReady { qoi_id, calibration_ids, holdout_ids, uq_interval,
  sensitivity_summary }` — report may carry the "evidence-grade" label.
- `EvidenceBlocked { qoi_id, missing: [<reasons>] }` — report is still
  generated but marked "not evidence-grade" with the reasons enumerated.

The gate is mechanical, not editorial: humans cannot override to `Ready`
without adding the missing artefacts.

## 6. Retraction and drift

- When a model parameter is retuned (calibration data updated), all
  downstream evidence claims are automatically demoted to Engineering
  until BCFD-091 re-runs.
- When a code path referenced by a QOI's method string changes materially,
  the same demotion applies.
- Retracted evidence claims are recorded in
  `docs/archive/<date>-retraction-<slug>.md` with the reason; a pointer
  is added to the affected VB entry in
  [VALIDATION_BIOPROCESS.md](VALIDATION_BIOPROCESS.md).

## 7. Forbidden shortcuts

- Do NOT combine calibration and holdout data to produce "average"
  parameters; that erases the holdout and blocks evidence tier.
- Do NOT declare evidence tier on the basis of a "well-known correlation"
  without a dataset registry entry.
- Do NOT hide UQ intervals in an appendix; the report body shows the
  interval next to the point estimate.
- Do NOT calibrate on the reported operating point and validate on it —
  the holdout must be independent along the axis being decided (scale,
  regime, phase geometry).
