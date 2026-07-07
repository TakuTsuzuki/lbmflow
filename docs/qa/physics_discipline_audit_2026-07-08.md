# Physics-discipline audit — 2026-07-08

## Executive summary

- Bundles audited: X, Y, Z, W, V, U, T, Q, S
- PHYSICS.md entries checked: 69 dated/current entries
- Ban-list violations found: 0 direct case-keyed/calibrated/clamping/decorative physics violations; 5 provenance, registry, or coverage anomalies filed
- Provenance gate compliance: fail

Audit method: read `.claude/skills/lbmflow-physics-discipline/SKILL.md`,
`CLAUDE.md`, `docs/BIOPROCESS_PIVOT.md`, `docs/SPEC_BIOPROCESS_CORE.md`,
`docs/MODEL_RISK_MATRIX.md`, all `docs/PHYSICS.md`, plus the repo-required
`docs/PLAN.md` and `docs/VALIDATION_BIOPROCESS.md`; then grepped implementation
paths for calibrated constants, case identity, clamps/caps, fallback, GPU
unsupported handling, and QOI provenance.

## Per-bundle findings

### Bundle X (BCFD-010..012)

- PHYSICS.md entries: current physics stack distribution storage; 2026-07-08
  passive scalar ADE; 2026-07-08 single-phase stirred-tank IBM runner unit
  conversion.
- Grep results: GPU/backend distribution support returns structured
  `UnsupportedReason::NotImplemented` for phase/scalar on GPU; core
  `QoiProvenance` serializes six mandatory fields and has negative tests.
- Findings: core QOI provenance passes, but CLI manifest
  `QoiMethodDescriptor` uses a separate `QoiProvenance` with only
  `source_fields`, `averaging_window`, `units`, and `validation_tier`.
  It lacks mandatory `method` and `averaging_region`.
- Disposition: MAJOR-VIOLATIONS

### Bundle Y (BCFD-020..023)

- PHYSICS.md entries: 2026-07-08 single-phase stirred-tank IBM runner unit
  conversion; 2026-07-07 direct-forcing IBM full-step impulse and overlap
  mobility; 2026-07-08 wall-shear proxy diagnostics.
- Grep results: no `scenario.name`/case-keyed geometry branches found; STL
  import and unsupported geometry paths return structured errors; geometry
  `.min`/`.max` uses are grid/raster bounds, not transport clamps.
- Findings: no ban-list violation found in geometry/impeller paths.
- Disposition: CLEAN

### Bundle Z (BCFD-030..035)

- PHYSICS.md entries: 2026-07-08 single-phase stirred-tank IBM runner unit
  conversion; 2026-07-08 passive scalar ADE distribution; 2026-07-08 central
  finite-difference stress and shear QOI; 2026-07-08 wall-shear proxy
  diagnostics; 2026-07-08 behavior review — BCFD-030 stirred-tank smoke.
- Grep results: no calibrated comments or case-keyed branches; shear reports
  P50/P90/P95/P99/max plus fraction-above-threshold; mixing skips with reason.
- Findings: no direct ban-list violation. Validation tier remains screening,
  correctly below VB-01/VB-02/VB-03 Engineering.
- Disposition: CLEAN

### Bundle W (BCFD-040..045)

- PHYSICS.md entries: 2026-07-08 conservative Allen-Cahn phase-field
  transport; phase-field guards; density/viscosity and `J_rho`; constant-sigma
  surface tension; static contact-angle metadata; free-surface/degassing
  placeholder; BCFD-040..045 behavior review.
- Grep results: phase-field clipping is explicit policy and reports
  `clipped_fraction`; boundedness failures return structured diagnostics;
  degassing is evidence-blocked.
- Findings: no silent clamp or fallback found. Capability registry and
  LIMITATIONS still say phase-field is not implemented; logged under
  cross-cutting ANOM-PHY-2.
- Disposition: MINOR-FINDINGS

### Bundle V (BCFD-050..053)

- PHYSICS.md entries: 2026-07-08 oxygen scalar and Henry equilibrium;
  resolved-interface oxygen kL/a source; dynamic-gassing kLa fit; OUR
  reaction-source hooks.
- Grep results: non-negative oxygen clipping records `OxygenDiagnostics`;
  uncalibrated kL is rejected for Evidence; Henry constants and partial
  pressure are explicit inputs/defaulted only at model policy boundary.
- Findings: no silent transport clamp or evidence fallback found.
- Disposition: CLEAN

### Bundle U (BCFD-060..063)

- PHYSICS.md entries: 2026-07-08 cell tracer massless advection; shear damage
  exposure integral; microcarrier Schiller-Naumann mode; microcarrier drag
  reaction scatter; BCFD-U behavior review.
- Grep results: Schiller-Naumann `Re_p > 800` rejects; two-way mass loading
  rejects above 0.1; outside-grid microcarrier scatter errors instead of
  clamping; exposure QOI includes distributions.
- Findings: no ban-list violation found.
- Disposition: CLEAN

### Bundle T (BCFD-070..075)

- PHYSICS.md entries: point-bubble entity store; buoyancy; Schiller-Naumann
  drag; added mass; lift/wall/turbulent-dispersion placeholders; RK4 substep;
  bubble-to-liquid scatter; PBM disabled/constant kernels; PBM interfacial
  area and kLa; hybrid gas double-count policy; terminal-velocity behavior
  review.
- Grep results: placeholder coefficients are declared in PHYSICS.md and
  Evidence-blocked; high holdup and Re limits reject; bubble reaction scatter
  returns ledger error if force cannot balance.
- Findings: module-local `KlaProvenance` and `HybridGasMetadata` are not full
  BCFD-080 `QoiProvenance` metadata. If these reports are exposed as QOIs,
  they need wrapping in the core QOI schema before serialization.
- Disposition: MINOR-FINDINGS

### Bundle Q (BCFD-100..102)

- PHYSICS.md entries: current MPI model risk row; checkpoint and MPI behavior
  covered by LIMITATIONS rather than specific PHYSICS entries.
- Grep results: `MpiSolver::new_local` avoids global mask replication; phase
  and scalar slab output paths exist; checkpoint serializes scalar fields,
  phase field, and QOI accumulators.
- Findings: `MpiSolver::new_local` accepts `_initial_scalar_fn` but does not
  use it; checkpoint manifest still marks `rng=false` and `particles=false`,
  with no integrated cell-tracer or bubble checkpoint payload despite
  BCFD-060 and BCFD-070 now being present.
- Disposition: MAJOR-VIOLATIONS

### Bundle S (BCFD-080..084)

- PHYSICS.md entries: QOI schema behavior is mostly in `qoi.rs`; no dedicated
  PHYSICS.md entries for UQ interval combination or scale-up operating-window
  ranking.
- Grep results: core `QoiBundle`/`QoiScalar`/`QoiPercentiles` require full
  provenance; CLI manifest provenance does not; capability registry and
  LIMITATIONS are stale for several landed capabilities; scale-up evaluator
  ignores the selected mode and ranks constraints by largest relative
  violation, not by the documented priority.
- Findings: provenance gate fails for the manifest QOI-method surface; PHYSICS
  coverage is missing for UQ/scale-up decision rules; capability registry is
  not the machine-readable truth after the landed bundles.
- Disposition: MAJOR-VIOLATIONS

## Cross-cutting findings

- Duplicate PHYSICS entries: overlapping but not contradictory phase-field
  entries exist for W-VOF O1 and the 2026-07-08 conservative Allen-Cahn path.
  Keep as historical context unless future edits fork the formulas.
- Missing validity domains: scale-up operating-window ranking
  (`crates/lbm-core/src/scaleup.rs`) and UQ interval combination
  (`crates/lbm-core/src/uq.rs`) lack PHYSICS.md entries with validation
  anchors.
- Un-cited literature: no uncited active physical force closure found in the
  audited bundles. The scale-up ranking/UQ logic is process-decision logic,
  not a force closure, but still needs provenance because it affects reported
  operating windows.

## Anomalies for follow-up

- ANOM-PHY-1: CLI manifest QOI provenance omits mandatory `method` and
  `averaging_region`, `crates/lbm-cli/src/manifest.rs:93`; proposed ticket:
  BCFD-080-FIX-MANIFEST-PROVENANCE.
- ANOM-PHY-2: capability registry and LIMITATIONS still mark landed
  phase-field, point-bubble/PBM, and cell-exposure capabilities as
  unsupported/not implemented, `crates/lbm-cli/src/capabilities.rs:63`;
  proposed ticket: BCFD-002-FIX-CAPABILITY-REGISTRY-SYNC.
- ANOM-PHY-3: checkpoint manifest still reports `rng=false` and
  `particles=false` after cell tracers and bubbles landed,
  `crates/lbm-core/src/solver.rs:4587`; proposed ticket:
  BCFD-102-FIX-PARTICLE-RNG-CHECKPOINT.
- ANOM-PHY-4: `MpiSolver::new_local` accepts but ignores initial scalar
  callback, `crates/lbm-core/src/dist.rs:886`; proposed ticket:
  BCFD-100-FIX-LOCAL-SCALAR-INITIALIZATION.
- ANOM-PHY-5: scale-up evaluator ignores `ScaleUpMode` in ranking and has no
  PHYSICS.md provenance entry for the decision rule,
  `crates/lbm-core/src/scaleup.rs:68`; proposed ticket:
  BCFD-084-FIX-SCALEUP-PROVENANCE-AND-PRIORITY.
