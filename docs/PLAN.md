# LBMFlow Bioprocess CFD — Implementation Plan (BCFD-000..110)

Lifecycle: living (owning doc for the ticket ledger, milestone plan,
development protocol, merge-queue rules, and known traps).

**Product mission.** Turn LBMFlow into a bioprocess-specific CFD core (see
[BIOPROCESS_PIVOT.md](BIOPROCESS_PIVOT.md) and
[SPEC_BIOPROCESS_CORE.md](SPEC_BIOPROCESS_CORE.md)). Do not add generic
CFD features; do not prioritise GPU, FP16, or WASM GUI before the QOI
pipeline works.

**Ordering law.** Single-phase stirred tank → Np / mixing / shear QOI →
resolved gas-liquid → oxygen / kLa → cell exposure → UQ / evidence gate →
point bubble / PBM → MPI / parallel output. Breaking this ordering to
chase GPU, FP16, GUI, or generic geometry has been the failure mode of
prior planning iterations.

---

## 0. Milestones

Each milestone's exit criterion is that its tickets have landed AND their
VB validation groups have their acceptance bands green at Engineering tier
(no evidence-tier claim unless BCFD-091 passes).

### M0 — Bioprocess pivot and single-phase stirred tank

Tickets: BCFD-000..005, BCFD-010..012, BCFD-020..021, BCFD-030..035,
BCFD-080..081, BCFD-090, BCFD-110.

Exit criteria:

- 3D stirred-tank single-phase scenario runs end-to-end.
- Impeller creates flow; torque and Np are extracted with provenance.
- P/V is emitted with units and averaging window.
- Passive-scalar mixing-time QOI (t95 / t99) is emitted.
- Shear-rate field and P50/P90/P95/P99 exposure percentiles are emitted.
- `lbm report` produces a bioprocess report scaffold.
- Capability registry states this is screening / engineering only, not
  evidence-grade.

### M1 — Resolved gas-liquid and oxygen / kLa

Tickets: BCFD-022, BCFD-040..048, BCFD-050..053, BCFD-080..081, BCFD-090,
BCFD-110.

Exit criteria:

- Phase-field gas / liquid run possible (conservative Allen-Cahn, coupled
  to solver velocity, with mass-conservation and boundedness guards).
- Sparger gas injection with conservation ledger.
- Gas-holdup QOI with threshold / averaging metadata.
- Oxygen scalar transport with Henry interfacial flux.
- Synthetic kLa fit from dynamic gassing works and reports fit metadata.
- Report marks kLa as model-dependent and non-evidence unless a
  calibration + holdout pair exists.

### M2 — Cell exposure and scale-up decision layer

Tickets: BCFD-060..062, BCFD-082..084, BCFD-091..092, BCFD-110.

Exit criteria:

- Cell tracers record shear and oxygen exposure with percentile
  distributions.
- Microcarrier one-way mode with Schiller-Naumann validity enforced.
- Sweep runner aggregates QOI bundles with per-case provenance.
- Scale-up operating window can be computed; infeasible sets produce an
  explicit conflict table.
- Evidence gate blocks unsupported claims.
- CLI + MCP bioprocess surface complete.

### M3 — Engineering aeration with point bubbles / PBM

Tickets: BCFD-070..075, BCFD-090, BCFD-110.

Exit criteria:

- Point bubbles injected from sparger with volume ledger.
- Bubble force closures active with declared validity ranges.
- PBM computes d32 with bin-conservation diagnostics.
- kLa from PBM interfacial area emitted.
- Hybrid resolved-interface + point-bubble bookkeeping avoids
  double-counting.

### Hard cut lines

- Do NOT start BCFD-100..102 (MPI parallel work) before M0 and M1 are
  green.
- Do NOT start GPU bioprocess support before CPU correctness + QOI
  validation are green.
- Do NOT start GUI before BCFD-081 (report generator) is useful.
- Do NOT implement generic CAD meshing beyond BCFD-023 until the QOI
  pipeline proves its value on a stirred tank.
- Do NOT claim evidence-grade until BCFD-091 passes.
- Do NOT add model complexity without a QOI and a validation entry.

### M4 — Screening-tier hardening (2026-07-08, revised per plan review)

**Rationale.** BCFD-000..110 landed as scaffolding (100% of the ledger),
but the honest capability status after landing is:

- Registry: 4 Experimental (screening only), 5 Unsupported.
- VB validation: only VB-03 (shear reducer, synthetic Couette/Poiseuille)
  is Engineering GREEN. VB-01/02 lack integrated runners; VB-04/05 are
  Landed (quick smoke) but heavy validation gated `--include-ignored`;
  VB-06/07/08 have impl-anomalies discovered by adversarial VB-verify.
- End-to-end demo (2026-07-08) hit ANOM-DEMO-1: the single-phase
  stirred-tank runner produces NaN at step 32-33 on a realistic
  screening scenario (see `docs/qa/e2e_demo_2026-07-08.md`).
- Physics-audit filed 5 ban-list findings (`docs/qa/physics_discipline_audit_2026-07-08.md`).
- **Adversarial plan review filed 8 additional gates + corrected priorities**
  (`docs/qa/plan_review_2026-07-08.md`, 2026-07-08).

M4 closes the gap between "code landed" and "screening-tier decisions
possible with confidence" — **for the single-phase stirred-tank path
only**. Multi-physics wire-up (point bubbles / PBM / cell exposure /
phase-field-as-product) is out of M4 scope and moved to M5.

**M4 scope guard (per plan review):**
- No acceptance-band weakening. VB-01 remains ±15% + ±5% mesh
  convergence; do not restate as anything weaker under any ticket DoD.
- No silent scope creep. Anything that promotes a capability from
  Unsupported → Experimental other than `single_phase_stirred_tank`
  and its dependencies belongs to M5, not M4.
- Every M4 ticket carries the physics-discipline stop-rule + ban clause
  by reference to `.claude/skills/lbmflow-physics-discipline`.

Tickets (revised priorities per review):

**P0 mechanical gates (nothing downstream can be trusted until these):**

- **BCFD-200 — ANOM-DEMO-1 root-cause + finite E2E (narrow).**
  Single-phase runner NaN at step 32-33 for the realistic stirred-tank
  scenario (250 rpm, 2 L, 48³ D3Q19). Symptoms: `scalar_cv.csv:32` NaN,
  `torque_force.csv:33` NaN, stress QOI aborts at
  `crates/lbm-cli/src/runner.rs:1659`. First step: add a deterministic
  regression that stops on the first non-finite field with a compact
  state dump (hydro, IBM force, scalar ADE, solid mask, stress inputs).
  Then classify: unit-feasibility reject / IBM force staging bug /
  scalar stability bug / geometry-mask bug. Fix at the classified layer.
  **Do NOT change acceptance bands, add force caps, add scalar clipping,
  add case-keyed branches, or tune the resolution floor without a
  stop-rule report.** DoD: E2E demo (`stirred_tank_screening.json`)
  reruns to completion producing finite QOI JSON; PHYSICS.md entry
  records the root cause and the fix. VB-01 correlation validation
  belongs to BCFD-204, not this ticket.

- **BCFD-230 — Finite QOI serialisation gate.** Non-finite values must
  not serialise as ambiguous JSON `null`. Reject or explicitly skip
  every non-finite QOI before JSON/CSV/report; distinguish
  `value: null` because skipped vs serialisation failure because
  non-finite measured. DoD: no M4 report artefact can contain NaN,
  infinity, or ambiguous null for a measured QOI.

- **BCFD-231 — Screening feasibility policy.** The failed demo
  validated with `TAU_NEAR_HALF`, `GRID_REYNOLDS_HIGH`, and
  `SCALAR_DIFFUSION_UNSTABLE` warnings still passing. Define which
  feasibility diagnostics are warnings-only at Screening and which
  become hard rejects for product example scenarios. Update
  `stirred_tank_screening.json` so it either validates cleanly or is
  documented as a deliberate internal-only stress case with a marker.

- **BCFD-232 — Registry runnable-truth gate.** Registry status can only
  move to Experimental after (a) end-to-end product scenario runs,
  (b) QOI provenance is complete, (c) `LIMITATIONS.md` text agrees.
  Drift-guard test fails if registry ↔ LIMITATIONS ↔ runnable fixtures
  disagree. Applied to every future promotion in M5.

- **BCFD-210 (ANOM-PHY-1) — QOI provenance completion.** *Elevated to
  P0 per review* (SPEC §7 non-negotiable). CLI manifest QOI provenance
  omits mandatory `method` and `averaging_region` for several QOIs.
  Enforce serialisation failure for every QOI missing metadata; retrofit
  missing fields; extend the BCFD-110 drift guard.

- **BCFD-237 — Scale-up spec reconciliation.** *Must precede
  BCFD-203/213.* PLAN, VALIDATION_BIOPROCESS §VB-08, and code disagree
  on constraint priority order. Choose and document one default order
  (single home). Define what each `ScaleUpMode` variant changes. Keep
  quantitative violation magnitude as separate report data, not
  conflated with priority.

**P1 VB impl-anomaly closures (after P0 gates; VB-VV series):**

- **BCFD-201 (VB-VV-001) — VB-06 equilibrium fit semantics.**
  `dynamic_gassing_kla_fit` returns `SkippedQoi` on C=C*; adversarial
  VB-06 expects `kLa ≈ 0` within tolerance. Add a distinct
  `KlaDynamicFitOutcome::SteadyZero` variant that preserves zero kLa
  with explicit method+provenance; **do NOT fake `R²=1`** for a case
  with no fit variance. Update MODEL_RISK_MATRIX §4.

- **BCFD-202 (VB-VV-002) — Percentile method freeze.** Live reducer is
  in `crates/lbm-core/src/stress.rs` (not only `qoi.rs`) and affects
  stress, damage, cell exposure, and every P50/P90/P95/P99 surface.
  Expose `PercentileMethod::{ Interpolated, NearestRank }` on
  `PercentileSummary`; default to `NearestRank` for the bioprocess
  report path; require the selected method in `QoiProvenance.method`.
  Audit `damage.rs` and report generation for consistency.

- **BCFD-203 + BCFD-213 (bundled) — Scale-up constraint ranking + mode +
  tip-speed field + provenance.** *After BCFD-237.* Add
  `tip_speed_max_m_per_s` to `ConstraintSet`. Implement the reconciled
  priority order from BCFD-237. `ScaleUpMode` actually changes ranking.
  Reports show both process-priority ordering AND quantitative
  tightness/violation. PHYSICS.md entry cites the source of the
  priority order.

- **BCFD-204 (VB-VV-004) — VB-01 Np validation harness.**
  *Prerequisite: BCFD-200 landed.* Public integrated
  `bioprocess_np_validation_run(scenario_path)` API driving a
  stirred-tank scenario to steady state and extracting Np. **VB-01 band
  is ±15% correlation + ±5% mesh convergence (NOT ±20%).** Reference
  Rushton Np ≈ 5.0 (Re > 10⁴, T/D=3) and PBT-45 Np ≈ 1.3. Includes
  steady-state detection, torque averaging convention, reference
  geometry fixtures, mesh sequence, runtime-budget policy.

- **BCFD-205 (VB-VV-005) — VB-02 mixing-time validation harness.**
  Same shape as BCFD-204 for a point-pulse scalar. Reference published
  Nθ correlation. Includes Sc-appropriate diffusivity, injection
  volume/site convention, CV averaging region.

**P1 M4 quality gates (per plan review):**

- **BCFD-233 — VB heavy-status ledger.** Distinguish
  `ignored-with-calibration-reason` from `ignored-with-no-harness` from
  `ignored-with-impl-anomaly`. Update `docs/VALIDATION_BIOPROCESS.md`
  and `lbm verify` output to expose the honest reason per VB group.

- **BCFD-234 — Mandatory behavior-review gate for M4 demos.** Every
  M4 end-to-end demo (`stirred_tank_screening.json`, and any future
  M4 example) requires a recorded behavior-validity review (dominant
  mechanism named, resolved vs closure separated, artefacts cited,
  verdict) before that demo can justify a registry status.

- **BCFD-235 — Runner dispatch matrix.** Define the mapping from
  `PhysicsModel::{SinglePhase, PassiveScalar, Oxygen, ResolvedPhaseField,
  PointBubble, Hybrid, CellTracer}` to runner backends AND to
  structured-Unsupported errors. CLI + MCP + sweep must all route the
  same scenario to the same runner or the same error. Prevents M5
  wire-ups from silently bypassing the product path or working through
  only one surface.

- **BCFD-236 — Example scenario stability fixtures.** `examples/bioprocess/*.json`
  are now part of the trust boundary. Add release-mode smoke fixtures
  with finite-QOI assertions and explicit expected warnings/limitations.
  Prevents regression back into non-finite outputs.

**P2 (deferred out of M4 unless a user commitment attaches):**

- **BCFD-211 (ANOM-PHY-3) — Checkpoint manifest particle/RNG flags.**
  Only enters M4 if checkpoint/restart is a user promise this release.
  Otherwise: fix in M5 when the particle/bubble state actually needs
  cross-run persistence.

**Out of M4 (moved to M5 or deferred to their proper track):**

- **BCFD-212 (ANOM-PHY-4) — MpiSolver::new_local scalar init.** MPI is
  behind a hard cut line. Routed to the MPI track; not M4.
- **BCFD-220 (`point_bubbles` wire-up)** — M3 bring-up, not screening
  hardening. → M5.
- **BCFD-221 (`pbm` wire-up)** — Same. → M5.
- **BCFD-222 (`cell_exposure` wire-up)** — → M5.
- **BCFD-223 (`phase_field_vof` wire-up)** — → M5.

**Milestone exit criteria (revised).** M4 GREEN means:

- ANOM-DEMO-1 closed. `stirred_tank_screening.json` runs to completion,
  producing finite Np / P/V / t95 / shear percentiles (no NaN, no
  ambiguous null).
- BCFD-230 finite-QOI gate active. BCFD-231 feasibility policy applied
  to every example scenario. BCFD-232 registry-truth drift-guard
  active.
- BCFD-210 QOI provenance completion + drift guard extended.
- BCFD-237 scale-up spec reconciliation landed BEFORE any 203/213
  implementation.
- VB-01, VB-02, VB-06, VB-07, VB-08 either Engineering GREEN
  or ledger-tagged as `ignored-with-calibration-reason` /
  `ignored-with-no-harness` / `ignored-with-impl-anomaly` per BCFD-233
  (no false "GREEN" claims).
- Physics-audit re-run passes with zero MAJOR-VIOLATION and MINOR
  findings on exposed QOIs are either fixed or explicitly ledger-tagged.
- Corrected E2E demo `stirred_tank_screening` has a recorded
  behavior-validity review per BCFD-234.
- Registry: `single_phase_stirred_tank`, `rotating_ibm`, `passive_scalar`
  Experimental with `stirred_tank_screening.json` as their runnable
  fixture. `oxygen_kla` stays Experimental only if a synthetic oxygen
  demo also has a corrected E2E run + behavior review; otherwise it
  demotes to Unsupported until M5.

### M5 — Multi-physics wire-up (post-M4)

**Rationale.** After M4 hardens the single-phase screening path,
M5 promotes the additional physics stacks to Experimental in the
registry. Each ticket must pass the BCFD-232 registry-truth gate.

- **BCFD-223 — `phase_field_vof` capability wire-up (M5-P0).**
  Product runner dispatch for `PhysicsModel::ResolvedPhaseField`.
  Feasibility from BCFD-004 (Cn, Pe_φ, density ratio, W). VB-04 heavy
  validation runnable. Split from any point-bubble or aerated-oxygen
  work — one bring-up at a time.

- **BCFD-222 — `cell_exposure` wire-up (M5-P1).** Bind `CellTracer` +
  `ShearDamageModel` into the runner. **Blocked on BCFD-202 landing**
  (percentile method) so exposure percentiles don't change under
  users.

- **BCFD-220 — `point_bubbles` wire-up (M5-P1).** Product path for
  `PhysicsModel::PointBubble`. Precedes BCFD-221 (PBM d32 needs a
  gas-source ledger from bubbles).

- **BCFD-221 — `pbm` wire-up (M5-P1, after BCFD-220).**

- **BCFD-224 — Aerated stirred-tank product path (M5-P2, after
  BCFD-223).** End-to-end path for `aerated_stirred_tank_screening.json`
  including gas holdup + oxygen + kLa. VB-05 heavy validation.

Hard cut line: no capability moves to Engineering (Tier 1) in M5.
That requires calibration/holdout + UQ + sensitivity per
CREDIBILITY_BIOPROCESS.md — deferred to a future evidence campaign,
not scheduled here.

---

## 1. Ticket ledger (BCFD-000..110)

Every ticket carries: **Depends** · **Targets** (files) · **Impl** ·
**Tests** · **DoD** (definition of done). Multi-file tickets bundle
same-file items into one codex order. Concurrency and worktree isolation
are the codex-dispatch skill's job, not the ticket's.

### Foundation (M0-blocking)

**BCFD-000 — Product pivot and repository guardrails.**
Depends: none.
Targets: README.md, docs/PLAN.md, docs/LIMITATIONS.md, docs/BIOPROCESS_PIVOT.md,
crates/lbm-cli/src/main.rs.
Impl: pivot docs, unsupported-as-product-grade labels, legacy-preset
warning string.
DoD: docs no longer imply general-purpose LBM production-readiness;
LIMITATIONS distinguishes demo / engineering / evidence; existing tests
pass.

**BCFD-001 — Bioprocess core specification docs.**
Depends: BCFD-000.
Targets: docs/SPEC_BIOPROCESS_CORE.md, docs/VALIDATION_BIOPROCESS.md,
docs/CREDIBILITY_BIOPROCESS.md, docs/MODEL_RISK_MATRIX.md.
Impl: define intended / forbidden use, tiers, VB-01..VB-08 scaffolding,
credibility policy, model-risk table.
DoD: four docs exist; each has current-status and not-yet-validated
sections.

**BCFD-002 — Capability registry and unsupported-combination errors.**
Depends: BCFD-000.
Targets: crates/lbm-scenario/src/lib.rs, crates/lbm-cli/src/{main,capabilities}.rs,
crates/lbm-core/src/solver.rs.
Impl: `CapabilityRegistry`, `CapabilityStatus`, `CapabilityTier`,
`UnsupportedReason`; enumerate capabilities including single-phase
stirred, rotating IBM, passive scalar, phase-field VOF, oxygen kLa, point
bubbles, PBM, cell exposure, evidence-tier report; structured errors
replace string-only errors on unsupported paths; add
`lbm capabilities --json`.
DoD: users can query support; no silent fallback for bioprocess modes.

**BCFD-003 — Bioprocess scenario schema v1.**
Depends: BCFD-001, BCFD-002.
Targets: crates/lbm-scenario/src/{bioprocess,lib}.rs,
crates/lbm-cli/src/schema.rs.
Impl: `BioprocessScenario` alongside legacy `Scenario`; sections `reactor`,
`fluids`, `operation`, `physics`, `cells`, `qoi`, `run`, `outputs`;
`ReactorSpec::StirredTank`, impeller / baffle / sparger lists;
`PhysicsSpec` discriminated union; unknown-field rejection;
impossible-combination rejection.
DoD: `lbm schema --bioprocess` emits schema; old scenarios still parse;
strict validation on bioprocess.

**BCFD-004 — Unit conversion and dimensionless feasibility layer.**
Depends: BCFD-003.
Targets: crates/lbm-scenario/src/{units,bioprocess}.rs,
crates/lbm-cli/src/validate.rs.
Impl: SI ↔ lattice for D, T, H, V, N, rpm, gas flow, vvm, ρ, μ, σ;
compute Re, Fr, We, Eo, Mo, Sc, Pe, St, Ma_lattice, Cn, Pe_φ; matching
priority (Re → density/viscosity ratio + We/Eo → Fr → Sc/Pe/Da → St);
feasibility diagnostics (Ma_lattice warn > 0.1 hard-reject > configured;
τ→0.5; interface too thin; bubble under-resolved; scalar diffusion
unstable); emit `UnitReport`.
DoD: every bioprocess scenario emits a unit feasibility report; invalid
lattice mappings fail before simulation.

**BCFD-005 — Verification command and machine-readable capability output.**
Depends: BCFD-002.
Targets: crates/lbm-cli/src/{main,verify,capabilities}.rs.
Impl: `lbm verify --tier {quick|bioprocess|full}`;
`lbm capabilities --json`; output includes `tests_run`, `tests_skipped`,
`unsupported_capabilities`, `validation_tier`, `git_sha`, `build_features`.
DoD: agents can query validation state without reading docs.

### Backend generalisation (M0-blocking)

**BCFD-010 — Backend fields generalisation for multiple distributions.**
Depends: BCFD-002.
Targets: crates/lbm-core/src/{fields,backend,solver}.rs,
crates/lbm-core/src/gpu/backend.rs.
Impl: support hydrodynamic `f`, phase-field `g`, scalar `h[k]`;
`DistributionKind::{Hydro, Phase, Scalar}`; allocation API
`enable_phase_distribution()`, `enable_scalar_distribution(name, D)`;
hydro-only path bit-identical; CPU supports phase / scalar first; GPU
rejects phase / scalar explicitly (structured error).
DoD: multiple distributions coexist in CPU solver; existing tests pass
unchanged.

**BCFD-011 — Material-property fields: ρ, μ, σ, phase fractions.**
Depends: BCFD-010.
Targets: crates/lbm-core/src/{fields,solver,materials,params}.rs.
Impl: `MaterialFields` (ρ_phys, μ_phys, ν_phys, σ, α_liquid, α_gas);
`MaterialModel::{SinglePhase, PhaseFieldMixture, ActiveScalarFeedback}`;
interpolation `ρ(φ) = ρ_g + φ(ρ_l - ρ_g)`; harmonic μ default with
explicit-flag linear μ fallback; material update after phase/scalar
updates; no update in hydro-only path.
DoD: physical material fields exist independent of lattice ρ.

**BCFD-012 — Solver-state manifest and QOI provenance.**
Depends: BCFD-003, BCFD-010.
Targets: crates/lbm-core/src/solver.rs, crates/lbm-cli/src/{runner,manifest}.rs.
Impl: extend manifest with `scenario_hash`, `bioprocess_schema_version`,
`backend`, `lattice`, `precision`, `active_models`, `qoi_methods`,
`unit_report`, `capability_report`; add `QoiProvenance` (source field,
averaging window, units, validation tier); every output file lists
manifest path.
DoD: every bioprocess run is auditable.

### Geometry (M0)

**BCFD-020 — Stirred-tank geometry templates.**
Depends: BCFD-003, BCFD-004.
Targets: crates/lbm-scenario/src/bioprocess.rs, crates/lbm-core/src/geometry.rs,
crates/lbm-core/src/solver.rs, crates/lbm-cli/src/presets.rs.
Impl: cylindrical tank voxel mask (flat bottom; dished bottom rejected
for M0); baffle template (count, width, thickness, wall-attachment);
coordinate conventions (x,y horizontal, z vertical, tank centre at
domain centre); geometry validation (impeller inside liquid, baffles
inside tank, sparger inside tank, enough grid resolution); build solid
mask and wall velocity fields.
DoD: bioprocess scenario can generate a 3D stirred-tank solid mask.

**BCFD-021 — Impeller and baffle geometry realisation.**
Depends: BCFD-020.
Targets: crates/lbm-core/src/{geometry,rotating_ibm}.rs,
crates/lbm-scenario/src/bioprocess.rs.
Impl: `ImpellerSpec::{Rushton, PitchedBlade, Marine, CustomMarkerSet}`
(custom is placeholder); realise Rushton and pitched blade as IBM
marker sets with blade thickness / angle; rotating marker generation;
baffles remain solid walls.
DoD: single-phase stirred-tank scenario supports rotating impeller +
baffles.

**BCFD-022 — Sparger geometry schema and masks.**
Depends: BCFD-020.
Targets: crates/lbm-scenario/src/bioprocess.rs, crates/lbm-core/src/geometry.rs.
Impl: `SpargerSpec::{Ring, Pipe, PointOrifices}`; gas inputs
(volumetric flow, vvm, `inlet_phase = gas`); orifice mask; no gas
injection yet; validation (sparger below liquid surface, orifice
resolved, positive gas flow, no raw φ inlets, no gas model → reject).
DoD: sparger geometry + metadata exist as a mask.

**BCFD-023 — Optional STL voxel import MVP.**
Depends: BCFD-020.
Targets: crates/lbm-core/src/voxel_import.rs,
crates/lbm-scenario/src/bioprocess.rs, crates/lbm-cli/src/runner.rs.
Impl: feature `geometry-import`; binary STL only; voxelise → solid
mask; patch labels (wall, impeller, baffle, sparger, unknown); unknown
allowed for screening tier only; evidence tier requires explicit
labels.
DoD: non-template geometry can be ingested experimentally.

### Single-phase runner and QOI (M0)

**BCFD-030 — Single-phase stirred-tank runner path.**
Depends: BCFD-021.
Targets: crates/lbm-cli/src/runner.rs, crates/lbm-scenario/src/bioprocess.rs,
crates/lbm-core/src/solver.rs.
Impl: run path for `BioprocessScenario` single_phase; 3D D3Q19 or
D3Q27; install tank / baffle solids, rotating impeller IBM; record
torque, force, velocity snapshots; reject 2D scenarios.
DoD: minimal stirred-tank scenario runs; rest tank quiescent; impeller
non-zero velocity; old scenario tests pass.

**BCFD-031 — Torque, Np, P/V, Nq QOI extraction.**
Depends: BCFD-030.
Targets: crates/lbm-core/src/qoi.rs, crates/lbm-cli/src/{runner,output}.rs.
Impl: compute Tq, `P = ωTq`, `N = ω/(2π)`, `Np = P/(ρN³D⁵)`,
`P_over_V = P/V_working`, `Nq = Q/(ND³)` when discharge surface
defined (skip with reason otherwise); output qoi_power.{csv,json}.
DoD: power QOIs emitted with units and averaging window.

**BCFD-032 — Stress tensor and shear-rate field.**
Depends: BCFD-011, BCFD-030.
Targets: crates/lbm-core/src/{stress,solver}.rs, crates/lbm-cli/src/output.rs.
Impl: velocity-gradient `S`; `gamma_dot = √(2 S:S)`; `viscous_stress =
μ_eff · gamma_dot`; second invariant, von Mises proxy; percentile
reducer (P50/P90/P95/P99/max) — never max alone; non-equilibrium
stress optional later.
DoD: shear fields + percentiles exist for single-phase runs.

**BCFD-033 — Wall shear and y+ diagnostics.**
Depends: BCFD-032.
Targets: crates/lbm-core/src/{wall_model,stress}.rs, crates/lbm-cli/src/output.rs.
Impl: wall-adjacent detector; wall-distance estimate; wall-shear proxy
from tangential velocity gradient; y+ when LES on; wall_shear.csv;
label as "proxy" unless a validated mode enabled.
DoD: wall-shear diagnostics available and clearly labelled.

**BCFD-034 — Passive scalar ADE distribution.**
Depends: BCFD-010, BCFD-030.
Targets: crates/lbm-core/src/{scalar,solver}.rs,
crates/lbm-scenario/src/bioprocess.rs.
Impl: scalar `h`; diffusivity, initial field, point/region pulse,
inlet source, no-flux wall; scalar evolves after fluid step; closed
domain mass conservation.
DoD: passive scalar simulated in stirred tank.

**BCFD-035 — Mixing-time QOI.**
Depends: BCFD-034.
Targets: crates/lbm-core/src/qoi.rs, crates/lbm-cli/src/output.rs.
Impl: scalar CV(t); t95 (CV ≤ 0.05 CV₀), t99 (≤ 0.01 CV₀);
compartment-wise CV; skip when no pulse configured; outputs
mixing_time.json, scalar_cv.csv.
DoD: mixing-time QOI reproducible.

### Phase-field gas-liquid (M1)

**BCFD-040 — Conservative Allen-Cahn phase-field evolution.**
Depends: BCFD-010, BCFD-011.
Targets: crates/lbm-core/src/{phase_field,solver}.rs,
crates/lbm-scenario/src/bioprocess.rs.
Impl: solver-coupled phase update; `∂φ/∂t + ∇·(φu + J_φ) = 0`,
`J_φ = -M[∇φ - (4/W)φ(1-φ) n]`; φ ∈ [0,1] with diagnostics;
`PhaseFieldParams { interface_width, mobility, clipping }`; uses
current fluid velocity; no surface tension in this ticket.
DoD: phase field dynamically coupled to solver velocity.

**BCFD-041 — Phase-field mass conservation and boundedness guards.**
Depends: BCFD-040.
Targets: crates/lbm-core/src/{phase_field,solver,divergence}.rs.
Impl: diagnostics `total_φ`, `min_φ`, `max_φ`, `clipped_fraction`,
`interface_cells`; hard divergence guard (NaN, out-of-bounds after
clip, excessive drift); manifest fields; `run_guarded` support.
DoD: phase-field runs cannot silently diverge.

**BCFD-042 — High-density-ratio ρ/μ interpolation and J_ρ consistency.**
Depends: BCFD-040, BCFD-041.
Targets: crates/lbm-core/src/{materials,phase_field,solver}.rs.
Impl: compute ρ(φ), μ(φ); `J_ρ = (ρ_l - ρ_g) · J_φ`; same J_ρ used in
continuity diagnostics AND momentum advection correction (consistency
check); moderate density-ratio guard; high ratio > threshold remains
Experimental.
DoD: material property fields update from φ consistently.

**BCFD-043 — Surface tension force.**
Depends: BCFD-040, BCFD-042.
Targets: crates/lbm-core/src/{surface_tension,solver,params}.rs.
Impl: chemical-potential force for constant σ; curvature / normal
diagnostics; applied through existing Guo force composition; disabled
in single-phase; capillary Δt diagnostic.
DoD: resolved interface generates capillary force.

**BCFD-044 — Contact-angle and wettability boundary.**
Depends: BCFD-043.
Targets: crates/lbm-core/src/{wetting,phase_field,geometry}.rs.
Impl: wall contact-angle parameter; phase-field wall flux condition;
per-wall metadata; droplet-on-wall init helper; reject contact-angle
without phase field.
DoD: contact angle represented in phase-field path.

**BCFD-045 — Free-surface top boundary and degassing placeholder.**
Depends: BCFD-043.
Targets: crates/lbm-core/src/{free_surface,solver}.rs,
crates/lbm-scenario/src/bioprocess.rs.
Impl: top modes `closed_lid`, `free_surface` (engineering
experimental), `degassing_outlet` (placeholder: gas out, liquid in;
mass ledger records gas outflow); evidence tier rejects degassing.
DoD: top boundary is explicit and cannot be silently wrong.

**BCFD-046 — Sparger gas-injection boundary.**
Depends: BCFD-022, BCFD-040, BCFD-045.
Targets: crates/lbm-core/src/{sparger,phase_field}.rs,
crates/lbm-scenario/src/bioprocess.rs.
Impl: gas-only injection at sparger; `inlet_phase: gas` maps to φ = 0;
enforced conservation (gas volumetric flow, gas volume ledger,
pressure diagnostic); reject under-resolved orifice for resolved
injection; point-bubble mode is BCFD-070+.
DoD: resolved gas injection with conservation ledger.

**BCFD-047 — Resolved gas-volume and gas-holdup QOI.**
Depends: BCFD-041, BCFD-046.
Targets: crates/lbm-core/src/qoi.rs, crates/lbm-cli/src/output.rs.
Impl: `α_g = 1 - φ`; global / compartment / thresholded / raw holdup;
metadata (threshold, averaging volume, time window,
`method = resolved_phase_field`); outputs gas_holdup.{json,csv}.
DoD: resolved gas holdup is recomputable and metadata-complete.

**BCFD-048 — Phase-field validation suite.**
Depends: BCFD-040, BCFD-041, BCFD-043, BCFD-046.
Targets: crates/lbm-core/tests/bioprocess_phase_field.rs,
docs/VALIDATION_BIOPROCESS.md (VB-04 + VB-05).
Impl: static planar interface, advected droplet, Laplace law, gas
injection volume ledger, contact-angle smoke; validation bands added
to docs; tests written from docs, not implementation; heavy tests
`--include-ignored`.
DoD: phase-field has explicit validation status.

### Oxygen and kLa (M1)

**BCFD-050 — Oxygen scalar transport.**
Depends: BCFD-034, BCFD-040.
Targets: crates/lbm-core/src/{oxygen,scalar}.rs,
crates/lbm-scenario/src/bioprocess.rs.
Impl: oxygen as named scalar over the ADE infrastructure; C_L, C_star,
D_O2, boundary sources; may be liquid-phase-only initially.
DoD: oxygen scalar simulated as a named scalar.

**BCFD-051 — Henry equilibrium and interfacial oxygen flux.**
Depends: BCFD-043, BCFD-050.
Targets: crates/lbm-core/src/{oxygen,surface_tension,qoi}.rs.
Impl: Henry equilibrium; interfacial area density `a_local` from
phase-field δ-approximation; flux `S = kL·a_local·(C* - C)`; kL model
enum (constant | correlation_placeholder | calibrated); ledger;
evidence claims reject uncalibrated kL.
DoD: resolved-interface oxygen transfer path exists.

**BCFD-052 — kLa QOI from scalar transient.**
Depends: BCFD-050, BCFD-051.
Targets: crates/lbm-core/src/qoi.rs, crates/lbm-cli/src/output.rs.
Impl: fit `dC/dt = kLa(C* - C)`; output kla_1_per_s, kla_1_per_hr,
fit_r2, fitting_window, method, CI when available; skip with reason
on invalid fit.
DoD: kLa estimable from oxygen scalar transient.

**BCFD-053 — OUR and reaction-source hooks.**
Depends: BCFD-050.
Targets: crates/lbm-core/src/{reaction,oxygen}.rs,
crates/lbm-scenario/src/bioprocess.rs.
Impl: OUR options constant / Monod placeholder / cell-density-scaled;
source term to oxygen; non-negative concentration guard; source
ledger.
DoD: oxygen consumption representable.

### Cell / microcarrier exposure (M2)

**BCFD-060 — Cell tracer population.**
Depends: BCFD-030, BCFD-032.
Targets: crates/lbm-core/src/{cells,particles}.rs,
crates/lbm-scenario/src/bioprocess.rs.
Impl: massless tracers sampling flow; record position, velocity,
gamma_dot, viscous_stress, ε proxy if available, oxygen C_L if
enabled; deterministic seeding; serialisable state.
DoD: cell trajectory / exposure history collectible.

**BCFD-061 — Shear exposure integral and damage model.**
Depends: BCFD-032, BCFD-060.
Targets: crates/lbm-core/src/{cells,damage}.rs, crates/lbm-cli/src/output.rs.
Impl: `E = ∫ max(0, τ - τ_c)^m dt`; alternatives (γ̇ threshold, ε
threshold placeholder); output P50/P90/P95/P99/max exposure,
fraction_above_threshold, residence_time_above_threshold — never max
alone.
DoD: shear damage risk computable from tracer histories.

**BCFD-062 — Microcarrier particle mode.**
Depends: BCFD-060.
Targets: crates/lbm-core/src/{microcarrier,particles}.rs,
crates/lbm-scenario/src/bioprocess.rs.
Impl: finite-size particles (d, ρ, restitution, drag, buoyancy);
reuse Schiller-Naumann within Re_p ≤ 800; suspension metrics (settled
fraction, height distribution, residence-near-impeller, shear
exposure); no collision yet.
DoD: microcarrier-like particles trackable one-way.

**BCFD-063 — Two-way particle reaction-force scatter.**
Depends: BCFD-062.
Targets: crates/lbm-core/src/{microcarrier,particles,solver}.rs.
Impl: two-way coupling; regularised kernel scatter; momentum ledger;
mass-loading guard until four-way exists.
DoD: two-way particle coupling with momentum conservation
diagnostics.

### Point-bubble and PBM (M3)

**BCFD-070 — Point-bubble entity store.**
Depends: BCFD-011, BCFD-022.
Targets: crates/lbm-core/src/bubbles.rs, crates/lbm-scenario/src/bioprocess.rs.
Impl: `Bubble { position, velocity, diameter, gas_volume, age }`;
`BubbleSet`; sparger seeding; deterministic injection schedule; no
forces yet.
DoD: point bubbles trackable as entities.

**BCFD-071 — Point-bubble force closures.**
Depends: BCFD-070.
Targets: crates/lbm-core/src/{bubbles,bubble_forces}.rs.
Impl: buoyancy, drag, added mass placeholder, lift placeholder, wall
lubrication placeholder, turbulent dispersion placeholder; each
declares validity range; substep supported.
DoD: point bubbles move under basic forces.

**BCFD-072 — Bubble-to-liquid momentum coupling.**
Depends: BCFD-071.
Targets: crates/lbm-core/src/{bubbles,solver}.rs.
Impl: scatter bubble force reaction to liquid with regularised
kernel; momentum ledger; gas holdup from point bubbles; high-holdup
guard until continuum mode exists.
DoD: point bubbles couple momentum with ledger.

**BCFD-073 — PBM bins and breakup/coalescence.**
Depends: BCFD-070.
Targets: crates/lbm-core/src/{pbm,bubbles}.rs,
crates/lbm-scenario/src/bioprocess.rs.
Impl: binned distribution; breakup / coalescence kernel traits;
placeholder kernels (disabled / constant); hooks for Luo-Svendsen and
Prince-Blanch; d32 and interfacial area updated.
DoD: PBM infrastructure computes d32 / a.

**BCFD-074 — kLa from point-bubble interfacial area.**
Depends: BCFD-050, BCFD-073.
Targets: crates/lbm-core/src/{kla,pbm,oxygen}.rs.
Impl: `a = 6 α_g / d32`; kL enum (constant / penetration placeholder /
calibrated); oxygen source `kL·a·(C* - C)`; metadata includes method,
kL source, d32 source; evidence tier requires calibrated kL.
DoD: kLa computable from PBM interfacial area.

**BCFD-075 — Hybrid resolved-interface + point-bubble bookkeeping.**
Depends: BCFD-047, BCFD-074.
Targets: crates/lbm-core/src/{hybrid_gas,qoi}.rs.
Impl: track gas in resolved phase field + point bubbles; prevent
double-counting; ε_g_{resolved, bubble, total}, a_{resolved, bubble};
evidence tier rejects hybrid until validation exists.
DoD: hybrid gas bookkeeping explicit.

### Reporting, credibility, sweep, evidence (M2)

**BCFD-080 — QOI output schema.**
Depends: BCFD-031, BCFD-035, BCFD-047, BCFD-052, BCFD-061.
Targets: crates/lbm-core/src/qoi.rs, crates/lbm-cli/src/output.rs,
crates/lbm-scenario/src/qoi_schema.rs.
Impl: `QoiBundle` (power, mixing, gas, oxygen, kla, shear, cells,
microcarriers, validation_status); every QOI has units, method,
time_window, averaging_region, source_fields, validation_tier;
qoi.json + CSVs for time series.
DoD: bioprocess QOI output stable and machine-readable.

**BCFD-081 — Bioprocess report generator.**
Depends: BCFD-080.
Targets: crates/lbm-cli/src/report.rs, docs/report_templates/bioprocess.md.
Impl: `lbm report <run-dir>`; sections (intended use, forbidden use,
scenario summary, unit feasibility, active models, QOI summary,
validation status, limitations, provenance); marks "not
evidence-grade" unless BCFD-091 passes; no PDF/docx.
DoD: human-readable bioprocess report generated.

**BCFD-082 — Calibration and holdout registry.**
Depends: BCFD-080.
Targets: crates/lbm-core/src/credibility.rs,
crates/lbm-scenario/src/bioprocess.rs, crates/lbm-cli/src/validate.rs.
Impl: `CalibrationDataset`, `HoldoutDataset` (id, QOI, source, date,
scale, operating condition, measurement uncertainty); same id cannot
be used for both; evidence tier requires holdout.
DoD: calibration/holdout separation enforced.

**BCFD-083 — UQ and sweep runner.**
Depends: BCFD-080, BCFD-082.
Targets: crates/lbm-cli/src/sweep.rs, crates/lbm-core/src/uq.rs,
crates/lbm-scenario/src/sweep.rs.
Impl: sweep scenario schema; deterministic grid; Latin hypercube
placeholder; aggregate QoiBundles; one-factor local sensitivity;
sweep_summary.json; failed cases recorded, not lost.
DoD: parameter sweeps produce aggregated QOI summaries.

**BCFD-084 — Scale-up operating-window evaluator.**
Depends: BCFD-031, BCFD-052, BCFD-061, BCFD-083.
Targets: crates/lbm-core/src/scaleup.rs, crates/lbm-cli/src/scaleup.rs.
Impl: constraints (kLa ≥ target, P/V ≤ limit, P95_shear ≤ limit,
mixing_time ≤ limit, gas_holdup range); modes (constant P/V, tip
speed, kLa, mixing time, custom weighted); explicit conflict table
when no feasible point.
DoD: scale-up decision summaries produced from QOI sweeps.

### Validation, gate, CLI/MCP surface (M2)

**BCFD-090 — Bioprocess validation matrix.**
Depends: BCFD-031, BCFD-035, BCFD-048, BCFD-052, BCFD-061.
Targets: docs/VALIDATION_BIOPROCESS.md, crates/lbm-core/tests/bioprocess_*.rs.
Impl: expand VB-01..VB-08 with setup, QOI, acceptance, tier, current
status; quick tests default, heavy `--include-ignored`.
DoD: validation coverage maps to every product QOI.

**BCFD-091 — Evidence-tier gatekeeper.**
Depends: BCFD-082, BCFD-083, BCFD-090.
Targets: crates/lbm-core/src/credibility.rs, crates/lbm-cli/src/{report,validate}.rs.
Impl: `EvidenceGate`; evidence tier requires validation-matrix pass +
calibration/holdout separated + mesh/time-step sensitivity + QOI
uncertainty interval + limitation report; failure → report marked
"not evidence-grade"; `EvidenceGateResult`.
DoD: evidence claims mechanically gated.

**BCFD-092 — CLI/MCP bioprocess tools.**
Depends: BCFD-080, BCFD-081, BCFD-091.
Targets: crates/lbm-cli/src/{main,mcp}.rs.
Impl: `lbm bioprocess {validate, run, qoi, report, sweep,
evidence-check}`; MCP tools `validate_bioprocess_scenario`,
`run_bioprocess_scenario`, `get_bioprocess_qoi`,
`generate_bioprocess_report`, `check_evidence_gate`; structured JSON;
human errors include remediation.
DoD: agent operates the bioprocess workflow end-to-end.

### MPI and scale (post-M1/M2 only)

**BCFD-100 — MPI memory localisation for bioprocess-scale runs.**
Depends: BCFD-010, BCFD-012.
Targets: crates/lbm-core/src/{dist,solver}.rs, docs/LIMITATIONS.md.
Impl: local geometry construction path for MPI; avoid global compact
arrays; closures for solid mask, wall velocity, material fields,
initial scalar; keep legacy path for small tests; per-rank memory
estimate in manifest.
DoD: bioprocess large-grid path not blocked by global-array
replication.

**BCFD-101 — Parallel field output.**
Depends: BCFD-080, BCFD-100.
Targets: crates/lbm-core/src/dist.rs, crates/lbm-cli/src/output.rs.
Impl: per-rank output files; rank-0 writes manifest; fields (velocity,
φ, oxygen, γ̇, gas_holdup); simple reader metadata; rank-0 gather kept
for validation-small cases.
DoD: parallel output no longer requires rank-0 full field gather.

**BCFD-102 — Checkpoint includes particles, scalars, statistics.**
Depends: BCFD-012, BCFD-034, BCFD-060, BCFD-070.
Targets: crates/lbm-core/src/{solver,cells,bubbles,qoi}.rs.
Impl: extend checkpoint format to include scalar/phase fields, cell
tracers, point bubbles, QOI accumulators, RNG/deterministic-injection
state; explicit mismatch errors; manifest reserved flags become true.
DoD: long bioprocess runs restart without losing statistics.

### Release sync

**BCFD-110 — Product README and limitation sync.**
Depends: all release-target tickets.
Targets: README.md, docs/LIMITATIONS.md, docs/VALIDATION_BIOPROCESS.md,
docs/SPEC_BIOPROCESS_CORE.md, docs/CREDIBILITY_BIOPROCESS.md.
Impl: README describes actual supported bioprocess tiers; capability
matrix generated from registry where possible; LIMITATIONS updated;
no GMP/CMC claims; evidence-tier only if BCFD-091 passes.
DoD: README, LIMITATIONS, and capability registry agree.

---

## 2. Validation-driven development protocol (preserved)

1. Acceptance criteria are stated numerically in
   [VALIDATION_BIOPROCESS.md](VALIDATION_BIOPROCESS.md), one entry per
   VB-XX group.
2. codex (or Opus / Sonnet) writes the tests **from the spec, without
   looking at the implementation**. Test order and implementation order
   never share a worktree.
3. If a gate fails:
   - **Engine bug** → fix the engine. The test is the source of truth.
   - **Spec is physically wrong** → run an experiment, revise the spec,
     record rationale in [PHYSICS.md](PHYSICS.md). Never fake the
     physics to pass a gate.
   - **Test bug** → send back to the test author with a repro. The
     implementation is not changed.
4. Bands frozen after a validation lands = **measured value + declared
   headroom**, printed in the assert message. The denominator /
   normalisation of every relative tolerance must be stated in the
   assert (per-component vs scale-relative flipped pass/fail on
   identical physics; this bit the pre-pivot session).
5. Anomaly pins carry an ANOM id in the assert message; the fix flips
   the pin **in the same commit** as the fix. Pin retighten is never
   deferred.

## 3. Merge-queue rules (hard requirements)

Every BCFD landing follows these; violations block the merge:

(a) Landing gate = `cargo test --workspace --release --no-fail-fast`
    with **UNPIPED exit code** (`; echo EXIT:$?`). Never pipe a gate
    through `tail` / `grep` — pipe eats the exit code.

(b) GPU landing evidence is **PM-run outside the codex sandbox**.
    In-sandbox Metal adapter access is intermittent. Codex GPU orders
    say "build + CPU gates only; report BENCH/GPU-PENDING".

(c) Audit / adversarial test orders **never share a worktree** with
    implementation orders.

(d) Band freezes = measured value + stated headroom; measured values
    printed in assert messages. Loosening needs PHYSICS.md rationale
    and owner sign-off.

(e) Current-wrong-value pins carry their ANOM id in the assert
    message; the fix flips the pin in the same commit as the fix.

(f) [Anomaly log entry](archive/2026-07-07-pivot/qa/anomaly-log.md
    format) required before merge for every P3-confirmed finding,
    including test-side dispositions. Bioprocess-era finding log will
    be re-established under `docs/BIOPROCESS_FINDINGS.md` when the
    first V&V loop closes.

## 4. Known traps (learned pre-pivot, all still active)

- **Backticks in inline codex order strings** die in zsh command
  substitution. Pass orders via file:
  `codex exec ... "$(cat <order-file>)" < /dev/null > /tmp/codex-<tag>.log 2>&1 &`.
- **Sandbox `git commit` fails intermittently** with `index.lock`
  EPERM on the shared `.git` in worktrees. Order text must include the
  "committed-ready fallback" clause — PM commits on codex's behalf at
  merge time.
- **Metal GPU adapter denial** in-sandbox is intermittent. GPU tests
  / bench are PM-run.
- **Loaded-window MLUPS false-negative trap**: on unified memory,
  background cargo/codex load halves-to-thirds GPU MLUPS. NEVER flip
  a perf gate RED (or dispatch kernel-opt orders) from a loaded-window
  number; re-measure quiet with A/B/A interleave.
- **`cargo test` fail-fast** masks regressions after the first
  failing binary. Landing gates MUST use `--no-fail-fast`.
- **Piping a gate through `| tail`** eats the exit code. Always
  separate: run gate raw, then `echo EXIT:$?`.
- **Keep-both merge resolution** via naive regex on `<<<<<<< HEAD` /
  `=======` / `>>>>>>>` markers can drop the closing `}` of the HEAD
  block when both sides are two full functions concatenated. Always
  `cargo build --workspace --release` before committing the merge.
- **Cargo test binary ordering** (alphabetical): earlier tests run
  first, hiding later regressions when fail-fast is on. Landing gate
  is `--no-fail-fast`.

## 5. Master agent prompt (paste into every codex BCFD order)

```
You are working on TakuTsuzuki/lbmflow, the bioprocess-specific CFD core
(2026-07-07 pivot). Read docs/BIOPROCESS_PIVOT.md, docs/SPEC_BIOPROCESS_CORE.md,
docs/VALIDATION_BIOPROCESS.md, docs/LIMITATIONS.md, and CLAUDE.md before
touching code.

Your work must implement only the assigned BCFD ticket in docs/PLAN.md.
Do not broaden scope. Do not add unsupported product claims. Do not
prioritize GPU, FP16, WASM GUI, or generic CFD capabilities unless the
ticket asks.

Rules:
  1. Existing tests must pass. Landing gate = cargo test --workspace
     --release --no-fail-fast with UNPIPED exit code.
  2. New physics must have a validation test or be marked Experimental
     / Unsupported.
  3. Unsupported combinations must fail loudly with structured errors,
     never silently fall back.
  4. Every QOI must include units, method, time window, averaging
     region, and validation tier.
  5. Evidence-tier claims require calibration + holdout + UQ +
     sensitivity records (BCFD-091 gate).
  6. Do not use Shan-Chen as the production gas-liquid path.
  7. Do not report only max shear; percentiles and exposure
     distributions are required.
  8. Do not silently fallback to CPU/GPU/precision changes.
  9. Keep old scenarios backward-compatible unless the ticket
     explicitly deprecates them.
 10. Update docs/LIMITATIONS.md for every new experimental or
     unsupported capability.

After coding:
  - Run cargo test --workspace --release --no-fail-fast ; echo EXIT:$?
  - Add or update validation tests (VB-XX group).
  - Update capability registry (BCFD-002 surface).
  - Update manifest/provenance if output behaviour changes.

Report unverified as unverified, skipped as skipped, failures with
their output. Fabricated progress is the worst possible failure.
```
