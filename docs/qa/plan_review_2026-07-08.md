# M4 Plan Review — 2026-07-08

## Executive judgement

- Overall verdict on M4 as stated: **Unsound**. The plan contains necessary
  anomaly-closure work, but it conflates screening hardening with new M3
  product-path bring-up, weakens at least one VB acceptance target, and does
  not make the capability registry mechanically truthful before adding more
  runnable claims.
- Confidence: **High**.

## Ticket-by-ticket critique

### BCFD-200 — ANOM-DEMO-1 root-cause + fix

- Priority sanity check: **P0 is correct**, but only for the finite-runner
  root cause, finite-output gate, and realistic screening demo rerun. The
  VB-01 correlation target belongs in BCFD-204, not this ticket.
- Scope sanity check: **Too broad as written.** Root-causing a step-32/33 NaN,
  rerunning the end-to-end demo, and making VB-01 comparable to Rushton are
  three different gates. The DoD also says Np within **±20%**, contradicting
  VALIDATION_BIOPROCESS.md VB-01 **±15% plus ±5% mesh convergence**. That is
  an acceptance-band weakening and must be removed.
- Hidden dependencies not stated in the plan:
  - A minimal non-finite reproduction harness with step-local dumps for hydro,
    IBM force, scalar ADE, solid mask, and stress inputs.
  - A finite-QOI serialization rule; the current demo wrote NaN power values
    as JSON `null` before failing later.
  - A decision on screening feasibility warnings that are currently allowed
    through despite `TAU_NEAR_HALF`, `GRID_REYNOLDS_HIGH`, and
    `SCALAR_DIFFUSION_UNSTABLE`.
  - A behavior-validity review after the corrected demo, not just metric
    completion.
  - No-physics-hack stop rule if the fix would require force caps, scalar
    clipping, case-keyed branches, or a tuned resolution floor.
- Bundling recommendation: **Split.** BCFD-200 should close NaN root cause and
  finite E2E only. VB-01 correlation belongs to BCFD-204. Finite output
  serialization deserves a separate cross-cutting P0 ticket.
- Alternative approach the PM might have missed: first add a deterministic
  failing regression that stops on the first non-finite field with a compact
  state dump. Then classify the failure as unit-feasibility reject, IBM force
  staging bug, scalar stability bug, or geometry-mask bug before changing any
  model behavior.

### BCFD-201 — VB-VV-001 equilibrium kLa fit

- Priority sanity check: **P1, not P0.** It blocks VB-06 Engineering status,
  but it does not block single-phase screening recovery.
- Scope sanity check: Achievable, but the DoD must specify how fit metadata is
  represented when there is no fit variance. Reporting `kLa=0` with fake
  `R2=1` would be a physics-discipline smell.
- Hidden dependencies not stated in the plan:
  - A `steady_equilibrium` or equivalent method/status so downstream reports
    do not mistake a zero-transfer steady case for a successful dynamic fit.
  - MODEL_RISK_MATRIX update for the equilibrium branch and validity domain.
  - QOI provenance that records the branch, time window, and reason R2 is not
    applicable or is defined by policy.
- Bundling recommendation: Can be bundled with other `qoi.rs` reducer work
  only if the same order owns all tests. Do not bundle with runner work.
- Alternative approach the PM might have missed: return a distinct
  `KlaDynamicFitOutcome::steady_zero` semantics, preserving zero kLa while
  making the non-regression test assert metadata rather than a meaningless
  regression statistic.

### BCFD-202 — VB-VV-002 percentile method freeze

- Priority sanity check: **P1 is appropriate.** It is necessary before VB-07
  Engineering, but not a root E2E blocker.
- Scope sanity check: Achievable, but underspecified. The live reducer is in
  `crates/lbm-core/src/stress.rs`, not only `qoi.rs`, and the choice affects
  stress, damage, cell exposure, and any report surface that says P50/P90/P95/P99.
- Hidden dependencies not stated in the plan:
  - Update all provenance strings to include the percentile method.
  - Decide whether historical interpolated percentiles remain available for
    numerical diagnostics.
  - Audit `damage.rs` and report generation so the method is consistent across
    shear and exposure QOIs.
- Bundling recommendation: Keep with stress/damage reducer changes. Do not
  combine with BCFD-201 unless same-file conflicts are explicitly managed.
- Alternative approach the PM might have missed: expose both methods, but make
  `NearestRank` the bioprocess-report default and require the selected method
  in `QoiProvenance.method`.

### BCFD-203 — VB-VV-003 constraint ranking + tip-speed field

- Priority sanity check: **P1 is appropriate.** Scale-up is not usable until
  this is coherent, but it is downstream of finite single-phase QOIs.
- Scope sanity check: Partly achievable, but the plan contradicts
  VALIDATION_BIOPROCESS.md. M4 says priority is `constant P/V -> tip speed ->
  kLa -> mixing time`; VB-08 says `constant kLa -> P/V -> tip speed -> mixing
  time` unless weights override. The code currently has hard-coded
  `kla -> P/V -> shear -> mixing -> gas holdup` and no `tip_speed_max_m_per_s`.
- Hidden dependencies not stated in the plan:
  - Spec reconciliation before implementation. The priority order must have
    one home.
  - A source for tip speed in `OperatingPoint.parameters` or a first-class QOI
    field.
  - Preservation of the existing P95 shear constraint; adding tip speed must
    not silently replace shear.
  - Report provenance for the ranking rule.
- Bundling recommendation: **Combine with BCFD-213** because they touch the
  same decision rule and PHYSICS.md provenance.
- Alternative approach the PM might have missed: separate “constraint priority”
  from “violation magnitude”. Reports should show both: a deterministic
  process-priority ordering and the quantitative tightness/violation.

### BCFD-204 — VB-VV-004 VB-01 Np validation harness

- Priority sanity check: **P1 after BCFD-200.** It becomes P0 only if M4 exit
  insists on Engineering GREEN rather than honest screening.
- Scope sanity check: Too compressed. A public validation API, steady-state
  detection, two published operating points, and three-grid convergence are
  not one small harness task.
- Hidden dependencies not stated in the plan:
  - Stable single-phase runner from BCFD-200.
  - Defined steady/statistical window and torque averaging convention.
  - Reference geometry fixtures for Rushton and PBT with documented T/D,
    baffles, clearance, and Re.
  - Runtime budget and ignored-heavy-test policy.
  - Mesh sequence and acceptance against both ±15% correlation and ±5%
    finest-grid convergence.
- Bundling recommendation: Split into harness API/fixture generation and
  heavy validation unignore. Do not bundle with the NaN fix.
- Alternative approach the PM might have missed: make the validation surface a
  CLI-driven fixture that runs the same product command and reads `qoi.json`,
  instead of inventing a test-only API that may drift from the user path.

### BCFD-205 — VB-VV-005 VB-02 mixing-time validation harness

- Priority sanity check: **P1 after BCFD-200.**
- Scope sanity check: Achievable only after scalar stability is known. The
  current demo has scalar CV becoming non-finite at step 32, so VB-02 harness
  work before scalar root cause risks building tests around a broken path.
- Hidden dependencies not stated in the plan:
  - BCFD-200 scalar stability classification.
  - A published Ntheta correlation with geometry family, not a generic
    “published band”.
  - Pulse injection semantics, CV compartment masks, and uniform-scalar skip
    behavior.
  - Halved-time-step invariance runner.
- Bundling recommendation: Split from BCFD-204 unless a single validation
  runner framework is created first.
- Alternative approach the PM might have missed: first freeze the scalar CV
  reducer and monotonicity checks on saved fields, then add full stirred-tank
  transient validation once the product run is finite.

### BCFD-210 — ANOM-PHY-1 QOI provenance completion

- Priority sanity check: **P0.** SPEC §7 makes provenance non-negotiable, and
  the CLI manifest currently has a second, weaker `QoiProvenance` type that
  omits `method` and `averaging_region`.
- Scope sanity check: Achievable if scoped to one canonical provenance schema.
- Hidden dependencies not stated in the plan:
  - Reconcile `crates/lbm-cli/src/manifest.rs::QoiProvenance` with
    `crates/lbm-core/src/qoi.rs::QoiProvenance`.
  - Cover skipped QOIs, not just measured QOIs.
  - Verify manifest, `qoi.json`, CSV sidecars, and report generation agree.
  - Add drift guard against future “nearly provenance” structs.
- Bundling recommendation: Standalone P0. It should land before any registry
  promotion or user-facing report claim.
- Alternative approach the PM might have missed: delete the CLI-local
  provenance struct or make it a thin wrapper around the core provenance type.
  Two provenance schemas invite exactly this drift.

### BCFD-211 — ANOM-PHY-3 checkpoint manifest particle/RNG flags

- Priority sanity check: **P2 / defer from M4.** It matters for restart
  correctness, but it is not necessary for a local screening demo unless M4
  promises restartable particle/bubble runs.
- Scope sanity check: The wording is dangerous. “Flip flags when corresponding
  data is actually serialised” is right; simply flipping flags because modules
  landed would make the manifest less truthful.
- Hidden dependencies not stated in the plan:
  - Actual checkpoint payloads for cell tracers, point bubbles, PBM state, and
    RNG streams.
  - Restart equivalence tests for particle and bubble trajectories.
  - Versioning/backward compatibility for checkpoint v3 or successor format.
- Bundling recommendation: Do not bundle with M4 runner stabilization. Defer
  until BCFD-220/222 actually serialize particle-like state.
- Alternative approach the PM might have missed: keep flags false and add an
  explicit `unsupported_sections` list until payloads exist. Truth beats an
  optimistic flag.

### BCFD-212 — ANOM-PHY-4 MpiSolver::new_local scalar init

- Priority sanity check: **P2 / out of M4.** PLAN hard cut lines say not to
  start MPI parallel work before M0/M1 are green. M4 is trying to make M0/M1
  screening finite.
- Scope sanity check: Achievable, but it requires MPI feature/toolchain
  coverage and is orthogonal to the product runner that failed.
- Hidden dependencies not stated in the plan:
  - Native MPI toolchain availability.
  - Per-rank scalar allocation path and halo behavior.
  - Feature-gated tests that actually run in CI or a documented manual gate.
- Bundling recommendation: Defer to the MPI track. Do not spend M4 critical
  path on it.
- Alternative approach the PM might have missed: file it as an MPI correctness
  bug with a small feature-gated regression, but keep capability status and
  M4 exit criteria independent of it.

### BCFD-213 — ANOM-PHY-5 scale-up mode ranking + provenance

- Priority sanity check: **P1**, and it should be inseparable from BCFD-203.
- Scope sanity check: Achievable after spec reconciliation. As of this branch,
  `evaluate_operating_window` uses a hard-coded priority order except for
  custom weights, so the real issue is not only “ranks by violation magnitude”
  in the PM text; it is that `ScaleUpMode` semantics and the documented
  priority are inconsistent/incomplete.
- Hidden dependencies not stated in the plan:
  - One authoritative priority order in VALIDATION_BIOPROCESS.md, PLAN.md, and
    PHYSICS.md.
  - Clear meaning for `ConstantPoverV`, `ConstantTipSpeed`,
    `ConstantKla`, and `ConstantMixingTime`.
  - UQ provenance if scale-up windows are later used in reports.
- Bundling recommendation: Combine with BCFD-203.
- Alternative approach the PM might have missed: make `ScaleUpMode` choose the
  primary invariant used to generate candidate operating points, while
  constraint ranking remains a separate documented reporting rule.

### BCFD-220 — point_bubbles capability wire-up

- Priority sanity check: **P2 / remove from M4.** This is new product-path
  bring-up, not screening-tier hardening.
- Scope sanity check: Too broad. A point-bubble entity store plus closures in
  core do not equal a runnable bioprocess path. The runner currently rejects
  anything without `SinglePhase`.
- Hidden dependencies not stated in the plan:
  - Product runner dispatch for `PhysicsModel::PointBubble`.
  - Sparger-to-bubble injection mapping and gas-volume ledger.
  - Bubble force validity guards surfaced through scenario errors.
  - Gas holdup and d32 QOI provenance through the canonical QOI schema.
  - VB coverage and behavior review for the end-to-end path.
- Bundling recommendation: Defer. If kept, split into a discovery/design
  ticket and a minimal runnable-path ticket, with registry promotion last.
- Alternative approach the PM might have missed: keep `point_bubbles`
  Unsupported and add a non-user-facing reducer/demo test first. Do not move
  registry status until a product scenario runs end-to-end.

### BCFD-221 — PBM capability wire-up

- Priority sanity check: **P2 / remove from M4.** It depends on BCFD-220 and
  validated gas bookkeeping.
- Scope sanity check: Too broad. PBM kernels include declared placeholders;
  wiring them into product d32/kLa outputs risks overstating model maturity.
- Hidden dependencies not stated in the plan:
  - Point-bubble or gas holdup source of alpha_g.
  - Kernel validity declarations in user-facing scenario/report output.
  - Bin conservation diagnostics in `qoi.json`.
  - kLa provenance if PBM interfacial area feeds mass transfer.
- Bundling recommendation: Must follow BCFD-220, not run in parallel with it.
- Alternative approach the PM might have missed: expose PBM d32 as an
  offline reducer over a synthetic bin state first, then promote to runner
  integration after gas-source ledger validation.

### BCFD-222 — cell_exposure capability wire-up

- Priority sanity check: **P1 if M4 promises cell-exposure screening; P2 if
  M4 is narrowed to single-phase hydrodynamic screening.** It should not be
  a P0 before finite hydrodynamics and percentile semantics are fixed.
- Scope sanity check: Too broad as written. Binding `CellTracer` and
  `ShearDamageModel` into the runner covers only shear exposure; oxygen
  exposure depends on oxygen fields that the current product runner does not
  advance.
- Hidden dependencies not stated in the plan:
  - BCFD-202 percentile method.
  - Tracer seeding policy, RNG provenance, and checkpoint/restart decision.
  - Single-phase velocity field stability from BCFD-200.
  - Damage-model calibration warning in reports.
  - Oxygen field availability for oxygen exposure, or explicit skipped QOI.
- Bundling recommendation: Split into shear-only tracer exposure for the
  single-phase runner and oxygen-exposure integration later. Do not bundle
  with point bubbles/PBM.
- Alternative approach the PM might have missed: for a one-week usable subset,
  keep `cell_exposure` Unsupported but emit Eulerian shear percentiles with
  correct provenance. That is less ambitious but more truthful.

### BCFD-223 — phase_field_vof capability wire-up

- Priority sanity check: **P1 for aerated screening, not P0 for single-phase
  recovery.** It is necessary only if M4 exit keeps the aerated demo as a
  required output.
- Scope sanity check: Too broad unless staged. Binding Allen-Cahn,
  sparger/free-surface behavior, gas holdup, oxygen, and kLa into the product
  runner is not a small wire-up.
- Hidden dependencies not stated in the plan:
  - Product runner dispatch for `PhysicsModel::ResolvedPhaseField`.
  - Phase-field feasibility from BCFD-004: Cn, Pe_phi, density ratio, and
    interface width.
  - VB-04/VB-05 heavy validations and the gas-holdup QOI schema.
  - Oxygen/kLa coupling if aerated report requires kLa.
  - Registry promotion rule that distinguishes module implementation from
    runnable product path.
- Bundling recommendation: Split into phase-field-only stirred-tank runner,
  sparger gas ledger, and oxygen/kLa report integration. Do not bundle with
  point bubbles/PBM.
- Alternative approach the PM might have missed: first add an end-to-end
  phase-field validation runner for droplet/Laplace and sparger ledger. Only
  then wire the aerated stirred-tank scenario.

## Cross-cutting findings

### Missing tickets (things the plan does not cover but should)

- **Finite QOI serialization gate.** Current artifacts can contain non-finite
  values that serialize as JSON `null`, which is indistinguishable from an
  intentional skipped QOI. Priority: **P0**. Dependencies: BCFD-200 and
  BCFD-210.
- **Capability registry truth gate.** M4 promotes capabilities by wire-up
  tickets but lacks a mechanical rule that a registry status can only move to
  Experimental after an end-to-end product scenario, QOI provenance, and
  limitation text agree. Priority: **P0**. Dependencies: BCFD-210 and any
  specific runner wire-up.
- **Screening feasibility policy.** The failed demo validated with severe
  stability warnings. M4 needs a ticket deciding which warnings stay warnings
  at Screening and which become hard rejects for product examples. Priority:
  **P0**. Dependencies: BCFD-200 root-cause classification.
- **Validation status ledger for ignored heavy tests.** M4 says all VB groups
  are Engineering GREEN or ignored-with-calibration-reason, but several are
  ignored because no integrated harness exists, not calibration. Priority:
  **P1**. Dependencies: BCFD-204/205/223.
- **Runner dispatch matrix.** The CLI/MCP/sweep surfaces call
  `run_bioprocess_single_phase`; the plan lacks a design ticket for how
  multiple `PhysicsModel` combinations map to runner backends and structured
  Unsupported errors. Priority: **P1**. Dependencies: BCFD-200, BCFD-223,
  BCFD-222 if kept.
- **Behavior-validity review gate for M4 demos.** M4 requires reports but not
  the post-run behavior review mandated by BIOPROCESS_PIVOT.md. Priority:
  **P1**. Dependencies: corrected E2E demos.
- **Example scenario stability fixtures.** The example JSONs are now part of
  the trust boundary; they need pinned smoke expectations and explicit
  warnings/limitations. Priority: **P1**. Dependencies: BCFD-200 and
  screening feasibility policy.

### Mis-prioritised tickets (P0 that should be lower / lower that should be P0)

- BCFD-210 should be **P0**, because provenance is a SPEC non-negotiable.
- BCFD-211 should be **P2** unless M4 includes checkpoint/restart as a user
  promise.
- BCFD-212 should be **P2/out of M4** because MPI is behind a hard cut line.
- BCFD-220 and BCFD-221 should be **P2/out of M4**; they are M3 bring-up, not
  hardening.
- BCFD-222 is **P1 or P2**, depending on whether M4 is narrowed to
  single-phase screening. It is not a prerequisite for fixing the NaN demo.
- BCFD-223 is **P1** only if the aerated demo remains in M4 exit criteria.
- A new finite-QOI gate and registry-truth gate should be **P0**.

### Ordering / dependency graph issues

- BCFD-200 must precede BCFD-204 and BCFD-205. Validating Np or mixing time
  against correlations before the runner is finite is wasted work.
- BCFD-210 should precede all registry promotions and all user-facing M4
  reports.
- BCFD-202 should precede BCFD-222 so cell-exposure percentiles do not change
  after runner integration.
- BCFD-203 and BCFD-213 should be one dependency unit, after the VB-08 priority
  contradiction is resolved.
- BCFD-223 must precede any aerated stirred-tank report that claims gas holdup
  or kLa from resolved gas-liquid fields.
- BCFD-220 must precede BCFD-221; PBM product d32/kLa needs a gas-source
  ledger.
- BCFD-211 should follow BCFD-220/222 only if those tickets actually serialize
  particle/RNG state.
- BCFD-212 should be detached from M4 and routed to the MPI track.

### Physics-discipline gaps in the plan itself

- M4 does not explicitly carry the stop-rule into each ticket. BCFD-200 is
  especially risky because the likely failure modes invite force caps, scalar
  clipping, resolution-floor tuning, or scenario-specific branches.
- BCFD-200 weakens VB-01 from ±15% to ±20%. The plan must not lower any
  acceptance band under the “screening-tier hardening” label.
- BCFD-201 must not invent fake fit quality for the equilibrium branch. A
  zero kLa steady case needs explicit method/provenance semantics.
- BCFD-203 contradicts VB-08's documented priority order. The spec must be
  fixed before code is changed.
- BCFD-220..223 risk treating module-level PHYSICS.md provenance as product
  evidence. PHYSICS entries for point bubbles, PBM, cells, and phase field do
  not by themselves justify moving registry status or emitting user-facing QOIs.
- M4 exit criteria require “zero MAJOR-VIOLATION” in a physics-audit rerun but
  do not say that MINOR findings affecting exposed QOIs must be either fixed
  or explicitly blocked before registry promotion.

## Answer to "when should this be shipped to users"

- **Is M4-GREEN enough for a screening-tier user to make a decision?** Only if
  M4 is corrected. As stated, no: it can become green by rushing broad
  wire-ups while still lacking a truthful registry, finite QOI gate, and
  validation harnesses. A corrected M4-GREEN would support qualitative
  screening comparisons only, not Engineering decisions. A real decision-maker
  needs VB-01/VB-02/VB-04..VB-08 Engineering GREEN as applicable, mesh and
  time-step sensitivity, QOI uncertainty intervals, calibration/holdout for
  kLa and damage models, and reports that separate skipped, unsupported, and
  measured QOIs.
- **Minimum useful subset in 1 week:** ship a single-phase stirred-tank
  screening path only. Required: BCFD-200 root cause fixed; finite-QOI gate;
  BCFD-210 provenance; screening feasibility policy; capability registry truth
  gate; BCFD-202 percentile method; BCFD-203/213 scale-up rule reconciliation
  if scale-up is shown; BCFD-204/205 at least as runnable ignored-heavy
  validation harnesses with honest status. Keep phase field, point bubbles,
  PBM, and cell exposure Unsupported unless their product scenarios actually
  run.
- **Minimum useful subset in 2 days:** do not ship a user decision tool. Produce
  an internal screening demo candidate only: fix or hard-reject the NaN
  scenario, block non-finite QOI serialization, complete manifest provenance,
  keep the registry conservative, and generate one single-phase report with a
  behavior-validity review. No aerated report, no cell exposure, no point
  bubbles/PBM promotion.

## Concrete follow-up ticket suggestions

### BCFD-230-finite-qoi-gate — Non-finite QOI rejection and skipped semantics

Depends: BCFD-200, BCFD-210.
Targets: `crates/lbm-core/src/qoi.rs`, `crates/lbm-cli/src/runner.rs`,
`crates/lbm-cli/src/report.rs`, tests under `crates/lbm-cli/tests`.
Impl: reject or explicitly skip every non-finite QOI before JSON/CSV/report
serialization; distinguish `value: null` because skipped from serialization
failure because non-finite measured input.
Tests: runner fixture with NaN torque/scalar/stress inputs fails with a
structured error; skipped Nq still serializes with reason and provenance.
DoD: no M4 report artifact can contain NaN, infinity, or ambiguous JSON null
for a measured QOI.

### BCFD-231-screening-feasibility-policy — Warning-to-reject policy for product examples

Depends: BCFD-200 root-cause classification.
Targets: `crates/lbm-scenario/src/units.rs`,
`crates/lbm-scenario/src/bioprocess.rs`, example scenario JSONs,
`docs/LIMITATIONS.md`.
Impl: define which feasibility diagnostics are allowed at Screening and which
must reject product examples; surface remediation strings.
Tests: examples that violate hard screening limits fail validation; the
curated screening example validates without hard rejects.
DoD: `stirred_tank_screening.json` no longer proceeds into a known unstable
unit mapping, or the mapping is documented as a deliberate internal-only stress
case.

### BCFD-232-registry-runnable-truth — Mechanical capability promotion gate

Depends: BCFD-210 and each capability-specific runner ticket.
Targets: `crates/lbm-cli/src/capabilities.rs`,
`crates/lbm-cli/tests/capabilities_drift_guard.rs`, `docs/LIMITATIONS.md`,
`README.md`.
Impl: require each Experimental capability to cite a runnable product scenario,
QOI provenance coverage, validation status, and limitation reason; Unsupported
entries must not say “not implemented” if the module exists but runner path is
missing.
Tests: drift guard fails if registry, LIMITATIONS, and runnable scenario
fixtures disagree.
DoD: `lbm capabilities --json` is the machine-readable truth for every M4
status.

### BCFD-233-vb-heavy-status-ledger — Honest VB ignored-test accounting

Depends: BCFD-204, BCFD-205, BCFD-223 where applicable.
Targets: `docs/VALIDATION_BIOPROCESS.md`, validation tests, `lbm verify`
output if present.
Impl: distinguish ignored because calibration/evidence is missing from ignored
because product harness is missing or implementation anomaly remains.
Tests: machine-readable verification output lists VB groups as GREEN, blocked
by calibration, blocked by harness, or impl-anomaly.
DoD: M4 cannot call a VB group “ignored-with-calibration-reason” when the real
blocker is absent runner/harness code.

### BCFD-234-e2e-behavior-review-gate — Mandatory behavior review for M4 demos

Depends: BCFD-200 and any demo runner ticket.
Targets: `docs/PHYSICS.md`, `docs/qa/e2e_demo_2026-07-08.md` successor,
report-generation docs.
Impl: add a required behavior-review checklist for every M4 end-to-end demo:
dominant mechanism, resolved vs closure terms, artifacts, and verdict.
Tests: documentation/check script or review template check ensures demo reports
include the behavior review before status promotion.
DoD: every M4 demo report has a recorded behavior-validity review before being
used to justify capability status.

### BCFD-235-runner-dispatch-matrix — Bioprocess physics combination dispatch

Depends: BCFD-200, BCFD-232.
Targets: `crates/lbm-cli/src/runner.rs`, `crates/lbm-cli/src/main.rs`,
`crates/lbm-cli/src/mcp.rs`, `crates/lbm-cli/src/sweep.rs`,
`crates/lbm-scenario/src/bioprocess.rs`.
Impl: define a dispatch matrix for `SinglePhase`, `PassiveScalar`,
`Oxygen`, `ResolvedPhaseField`, `PointBubble`, `Hybrid`, and `CellTracer`
combinations; unsupported combinations return structured errors with
capability ids.
Tests: CLI, MCP, and sweep all route the same scenario to the same runner or
the same structured unsupported error.
DoD: adding BCFD-220..223 cannot silently bypass the product path or only work
through one surface.

### BCFD-236-example-scenario-stability-fixtures — Curated screening examples as tests

Depends: BCFD-200, BCFD-231.
Targets: `examples/bioprocess/*.json`, `crates/lbm-cli/tests`, `docs/README`
or example docs.
Impl: add short release-mode smoke fixtures for curated screening examples
with finite QOI assertions and explicit expected warnings.
Tests: `stirred_tank_screening` completes with finite Np/P/V and shear
percentiles; aerated example remains skipped/unsupported until BCFD-223.
DoD: examples cannot drift back into non-finite outputs while tests stay green.

### BCFD-237-scaleup-spec-reconciliation — Single source for scale-up priority semantics

Depends: none; must precede BCFD-203/213 implementation.
Targets: `docs/SPEC_BIOPROCESS_CORE.md`, `docs/VALIDATION_BIOPROCESS.md`,
`docs/PHYSICS.md`, `crates/lbm-core/src/scaleup.rs` tests.
Impl: choose and document one default ranking order; define what each
`ScaleUpMode` changes; preserve quantitative violation magnitude as separate
report data.
Tests: VB-08 adversarial test asserts the documented order and custom-weight
override.
DoD: PLAN, VALIDATION, code tests, and PHYSICS.md no longer disagree about
scale-up priority.
