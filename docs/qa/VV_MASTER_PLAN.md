# V&V Master Plan — every angle, in parallel (2026-07-06)

Active routing table. Findings flow to `docs/qa/anomaly-log.md`; core/demo/
spec routing per the CLAUDE.md V&V rule. Machine budget: ≤5-8 compile-bound
codex orders at once — lanes are wave-scheduled, "parallel" means every
INDEPENDENT lane proceeds without waiting unless a dependency is stated.

Principles (frozen): physical rigor prime directive; adversarial and
implementation-separated tests; derive before blaming the engine; every band
has a behavior anchor; every experiment leaves a visual artifact; bands
frozen on (C) terms are not authoritative.

Status: **RUN** running · **DONE** delivered (regression-pinned) ·
**GATE** blocked on named dependency · **W2/W3** wave.

---

## Axis 1 — Adversarial analytic accuracy audits (Skill: lbmflow-accuracy-audit)

- 1.1 Momentum-exchange probe + walls — **DONE** (cx/audit-probe, 5 probes).
- 1.2 Rotating IBM — **GATE (ANOM-P4-001 core fix)** — 8 probes on
  cx/audit-ibm are the acceptance suite.
- 1.3 Particles CR-3 — **DONE** (cx/audit-particles, 6 probes).
- 1.4 Sources + patches CR-1/2 — **DONE** (cx/audit-sources, 5+1 probes).
- 1.5 Cumulant viscosity correction — **GATE (ANOM-P4-008 core fix)** —
  cx/audit-cumulant e2 canary (|a| ≤ 2e-3).
- 1.6 Penalized rotor — **GATE (ANOM-P4-010 core fix)** — cx/audit-rotor F1-F5.
- 1.7 **W2** Remaining subsystems: Shan-Chen force accuracy (pressure-tensor
  form, spurious-current tau dependence beyond T11); WALE strain
  reconstruction vs analytic (extend g2); convective outflow reflection
  coefficient curve; FP16/deviation storage error model; Bouzidi MOVING
  walls; Bouzidi mixed-qd probed force (queued kill-case from ANOM-P4-011
  triage). One codex order per subsystem, A1-A5 loop.
- 1.8 Tracer-limit particle audit (near-neutral buoyancy) — **GATE
  (Phase B spec)** — v→u relaxation in nonuniform flow; BBO validity floor.

## Axis 2 — Meta-V&V: does the net catch fish?

- 2.1 **Mutation testing of the validation suite** — **W1**: inject
  deliberate physics bugs one at a time (sign flip in Guo source, factor-2
  in MW term, w_q swap, off-by-half wall position, dropped corner link,
  tau→omega typo, feq u² coefficient); every MISSED mutant = a coverage hole.
  Deliverable: `docs/qa/mutation-coverage.md`; target 0 silent mutants in
  physics kernels.
- 2.2 Band-vacuity scan — **DONE** (commit 2e121c8). Retighten queue
  survives in `docs/qa/band-vacuity-scan.md`.
- 2.3 Behavior-anchor coverage — **W2**: for every VALIDATION.md band,
  check a pattern/sign/monotonicity anchor exists (two-layer rule); list
  band-only gates → codex retrofit order.

## Axis 3 — Cross-intelligence review

- 3.1 codex independent physics review — **DONE** (commit 714da6a).
  Converged on P4-001 / P4-008 / P4-010; F19/F20 (P4-011) refuted by
  derivation. Cumulant naming routed to D-track.
- 3.2 **W2** Code→spec back-translation — subagent derives "what equations
  does this code actually solve" per module; diff against PHYSICS.md.
- 3.3 **W2** Claims-ledger cross-check — every PHYSICS.md / whitepaper
  claim mapped to the test that proves it; red rows.
- 3.4 Review panel on fix designs — **GATE (P4-001/010 core landing)** —
  2-3 agent judge panel on the fix design before gate-running.

## Axis 4 — Literature & community bad-practice sweep (web)

- 4.1 **LBM pitfall compendium** — **W1**: forcing-scheme inconsistencies,
  non-eq initialization, corner/edge BC, checkerboard modes, tau→0.5
  oscillations, viscosity-dependent BB slip, momentum-exchange Galilean
  (Wen et al.), interpolated-BB mass leakage, Shan-Chen thermodynamic
  inconsistency, LES wall, cumulant Galilean claims, Ma² contamination,
  staircase artifacts. Map every item → {covered / gap / N-A}.
  Deliverable: `docs/qa/pitfall-checklist.md`; gaps become audit rows.
- 4.2 Benchmark-matrix expansion — **DONE** (commit 4eac49e); outstanding
  queue in `docs/qa/benchmark-backlog.md`. Kovasznay / Sangani-Acrivos /
  Womersley landed in 2c78d85.
- 4.3 **W2** V&V methodology audit — ASME V&V 20 / Oberkampf-Roy vocabulary;
  solution verification = grid-convergence reporting discipline (GCI).

## Axis 5 — Systematic experiment matrices + visualization

- 5.1 **Feature-interaction conservation matrix** — **W1**: pairwise (later
  triple) combinations of {gravity, SC multiphase, rotor, particles,
  sources, patches, LES, f32, 3D} on small grids; per cell: mass/momentum
  ledger closure, NaN watchdog, symmetry spot-check. Composition bugs
  (force-field overwrite/add rules) live here.
- 5.2 **W2** Randomized property sweep — proptest-style legal configs;
  invariants = conservation, boundedness, determinism, mirror equivariance.
- 5.3 **W2** Visual anomaly trawl — scripted matrix × params → PNG/VTK
  → sim-anomaly-scan + qa-viewer; PM behavior-review flagged cases.
- 5.4 **W2** Long-horizon soak — 10⁶-step runs of cavity / channel+cylinder /
  SC droplet / post-fix penalized blade rotor; drift ledgers.
- 5.5 **W3** Stability-envelope cartography — automated (tau, Ma, grid-Re)
  boundary mapping per collision operator; measured vs tune-stability
  thresholds.
- 5.6 **W2** Degenerate-geometry gauntlet — 1-cell gaps, diagonal
  staircases, obstacles touching rims/patches/seams, zero-thickness,
  moving-wall corners.

## Axis 6 — External referees

- 6.1 **W2** Cross-solver comparison — OpenLB (~/projects/cfd-bench) for
  cavity/cylinder/TGV.
- 6.2 **W3** Spectral mini-referee — 64-line Python spectral NS for 2D TGV
  / Kolmogorov.
- 6.3 GPU absolute physics — **GATE (R2-D3 lands)** — TGV order + cavity
  Ghia directly on GPU; PM-run.

## Axis 7 — Static & structural verification

- 7.1 **W1** Ban-list grep sweep — case-identity branches, transported-
  quantity clamps, bare calibrated literals; recurring CI-able check.
- 7.2 **W2** Numeric-literal provenance scan — every float in src/ physics
  paths needs a provenance comment or named constant.
- 7.3 **W3** Unsafe/precision audit — SIMD unsafe UB; f32↔f64 cast audit;
  Kahan-summation needs.
- 7.4 **W2** Docs/spec consistency — VALIDATION.md ↔ test files drift.

## Axis 8 — Interface/tooling V&V (product surface)

- 8.1 **W2** Scenario schema round-trip + fuzz.
- 8.2 **W3** WASM/GUI parity vs native.
- 8.3 **W3** CLI/MCP contract tests (presets/gallery/VTK/manifest schema).
- 8.4 **W3** Units V&V extension — SI↔lattice across all constructors.

## Dependencies & sequencing

- **W1 (dispatchable now)**: 2.1 mutation pilot, 4.1 pitfall compendium,
  5.1 interaction matrix, 7.1 ban-list script.
- **W2** dispatches as W1 lanes return; **W3** after.
- GATE lanes unblock on: core fixes (P4-001 / 008 / 010), R-Phase 2 D-3,
  Phase B spec, MF-γ.
- Failures route per the standing rule (core / demo / spec) and log to
  anomaly-log with ANOM ids. Every experiment lane ships visual artifacts.

## Reporting

Dashboard artifact (pass-4 board, extended per wave):
https://claude.ai/code/artifact/ced42aaf-5832-4dbd-8954-be9900140519.
Per-lane deliverables under `docs/qa/`. Cross-session routing: core-engine
(defects), D-track PM (demo/spec), user (dashboard + summaries).
