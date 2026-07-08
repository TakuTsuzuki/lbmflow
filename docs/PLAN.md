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

### M4 — Screening-tier hardening (2026-07-08, post-BCFD-000..110 landing)

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

M4 closes the gap between "code landed" and "screening-tier decisions
possible with confidence". No new ambitious scope; only stabilise what
BCFD-000..110 already promised.

Tickets:

**P0 blocker (nothing else can be Engineering-tier until this lands):**

- **BCFD-200 — ANOM-DEMO-1 root-cause + fix.** Single-phase runner NaN at
  step 32-33 for the realistic stirred-tank scenario (250 rpm, 2 L, 48³
  D3Q19). Symptoms: `scalar_cv.csv:32` NaN, `torque_force.csv:33` NaN,
  stress QOI aborts at `crates/lbm-cli/src/runner.rs:1659`. Suspects:
  (a) Bundle Z's scalar step interacts with Bundle X's material fields
  in a divergent way, (b) impeller IBM produces spurious velocity
  outside the CFL / low-Ma limit, (c) geometry validation resolution
  floor `STIRRED_TANK_SCREENING_MIN_CELLS = 12` is too coarse for the
  scenario's actual physics.  DoD: E2E demo reruns to completion, VB-01
  Np comparable to Rushton correlation ± 20%, PHYSICS.md entry records
  the root cause and the fix.

**VB impl-anomaly closures (BCFD-VV series, from `docs/qa/vb_verification_2026-07-08.md`):**

- **BCFD-201 (VB-VV-001) — VB-06 equilibrium fit.** `dynamic_gassing_kla_fit`
  returns `SkippedQoi` on the equilibrium case (C=C*); adversarial VB-06
  expects `kLa ≈ 0` within tolerance. Decision: change the fit to accept
  the equilibrium branch and emit `kLa=0` when the residual variance is
  below the steady_epsilon floor; document the decision in
  MODEL_RISK_MATRIX §4.
- **BCFD-202 (VB-VV-002) — VB-07 percentile method freeze.**
  `percentile_summary` uses linear interpolation between rank samples;
  adversarial VB-07 asserts nearest-rank percentiles. Decision: expose
  the method as a `PercentileMethod::{ Interpolated, NearestRank }` enum
  on `PercentileSummary`, default to `NearestRank` for the shear/exposure
  path (matches BCFD-032 spec §7 and canonical bioprocess reporting),
  and update VB-07 to freeze against nearest-rank. PHYSICS.md entry
  records the choice.
- **BCFD-203 (VB-VV-003) — VB-08 constraint ranking + tip-speed field.**
  `evaluate_operating_window` ranks conflicting constraints by violation
  magnitude; adversarial VB-08 expects the documented priority order
  (constant P/V → tip speed → kLa → mixing time). Also, the
  `ConstraintSet` type lacks a `tip_speed_max_m_per_s` field entirely.
  Decision: add the field, restore the documented priority order, and
  put the "custom weighted" mode behind an explicit request.
- **BCFD-204 (VB-VV-004) — VB-01 Np validation harness.** Add a public
  integrated `bioprocess_np_validation_run(scenario_path)` API that
  drives a stirred-tank scenario to steady state, extracts Np, and
  returns a Result usable by VB-01. Reference Rushton correlation
  Np ≈ 5.0 at Re > 10⁴ with T/D=3.
- **BCFD-205 (VB-VV-005) — VB-02 mixing-time validation harness.** Same
  shape as BCFD-204 for a point-pulse scalar. Reference published Nθ
  correlation for the geometry family.

**Physics-audit follow-ups (ANOM-PHY series, from `docs/qa/physics_discipline_audit_2026-07-08.md`):**

- **BCFD-210 (ANOM-PHY-1) — QOI provenance completion.** CLI manifest
  QOI provenance omits mandatory `method` and `averaging_region` for
  several QOIs. Enforce serialisation failure per BCFD-012 §7 for every
  QOI, retrofit missing fields, add a drift-guard test similar to
  BCFD-110's capability drift guard.
- **BCFD-211 (ANOM-PHY-3) — Checkpoint manifest particle/RNG flags.**
  `Manifest.reserved.{rng, particles}` still `false` after BCFD-060
  cell tracers and BCFD-070 point bubbles landed. Flip flags when
  the corresponding data is actually serialised; add roundtrip test.
- **BCFD-212 (ANOM-PHY-4) — MpiSolver::new_local scalar init.** The
  closure argument is accepted but ignored, silently defaulting to
  zero. Wire the closure through to the per-rank scalar allocation
  path. Feature-gated test.
- **BCFD-213 (ANOM-PHY-5) — Scale-up mode ranking + provenance.**
  `ScaleUpEvaluation` currently ignores `ScaleUpMode` when ranking; also
  lacks a PHYSICS.md entry for the decision rule. Bundled fix with
  BCFD-203 (same file surface). Add PHYSICS.md entry citing the source
  of the priority order.

**M3 wire-up to product path (currently Unsupported in the registry
despite Bundle T + U + partial S landed):**

- **BCFD-220 — `point_bubbles` capability wire-up.** Extend
  `BioprocessScenario.physics::PointBubble` to a runnable path
  (currently only the entity store exists). DoD: registry moves
  `point_bubbles` from Unsupported to Experimental; a screening scenario
  with sparger + point bubbles runs end-to-end producing gas holdup and
  d32 QOIs.
- **BCFD-221 — `pbm` capability wire-up.** Bind the PBM bins to the
  bioprocess runner's step; QOI outputs d32 with provenance.
- **BCFD-222 — `cell_exposure` capability wire-up.** Bind `CellTracer`
  and `ShearDamageModel` into the runner's step; emit exposure
  distributions per BCFD-061.
- **BCFD-223 — `phase_field_vof` capability wire-up.** Bind the landed
  Allen-Cahn phase field (BCFD-040..048) as a runnable
  `physics::ResolvedPhaseField` scenario path.

**Milestone exit criteria.** M4 GREEN means:

- ANOM-DEMO-1 closed. `lbm bioprocess run examples/bioprocess/stirred_tank_screening.json`
  completes without NaN. QOI outputs land.
- All VB-01..08 groups either Engineering GREEN or explicitly
  Ignored-with-calibration-reason (Evidence gate blocked, not Impl-anomaly).
- Registry: at least `single_phase_stirred_tank`, `rotating_ibm`,
  `passive_scalar`, `oxygen_kla`, `phase_field_vof`, `cell_exposure`
  all Experimental with an end-to-end runnable scenario.
- Physics-audit re-run passes with zero MAJOR-VIOLATION.
- End-to-end demo report.md generated for both `stirred_tank_screening`
  and `aerated_stirred_tank_screening`.

No M5/M6 declared in this document. Evidence-tier claims (Tier 2)
remain calibration-dataset-blocked and are not scheduled here.

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
