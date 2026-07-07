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

- 1.9 **Dynamic multiphase exact suite** — **RUN** (codex cx/mp-dynamics in
  flight): capillary-wave dispersion (omega^2 = sigma k^3/(rho_l+rho_v),
  k^3-law fit), Lamb n=2 droplet oscillation (R^-3/2 scaling), two-layer
  stratified Couette exact piecewise profile under shear, Lucas-Washburn
  x^2 ~ t imbibition. Research agent enumerating the next tier (Tomotika,
  Taylor deformation D(Ca), spinodal exponents, KH threshold, Jurin) →
  docs/qa/multiphase-hard-cases.md.

## Axis 2 — Meta-V&V: does the net catch fish?

- 2.1 **Mutation testing** — **RUN (split)**: cx/vv-mutation (other V&V
  worker) landed the runner + 3 mutants, all KILLED (moving-wall sign,
  Zou-He pressure normal sign, outflow stale slot); extension order queued
  (task #8) for the remaining mutants: inject
  deliberate physics bugs one at a time (sign flip in Guo source, factor-2
  in MW term, w_q swap, off-by-half wall position, dropped corner link,
  tau→omega typo, feq u² coefficient); every MISSED mutant = a coverage hole.
  Deliverable: `docs/qa/mutation-coverage.md`; target 0 silent mutants in
  physics kernels.
- 2.2 Band-vacuity scan — **DONE** (commit 2e121c8). Retighten queue
  survives in `docs/qa/band-vacuity-scan.md`.
- 2.3 Behavior-anchor coverage — **DONE** (structural-sweep-2026-07-06.md):
  17/23 sections dual-layer, T4 band-only justified; no retrofit needed.

## Axis 3 — Cross-intelligence review

- 3.1 codex independent physics review — **DONE** (commit 714da6a).
  Converged on P4-001 / P4-008 / P4-010; F19/F20 (P4-011) refuted by
  derivation. Cumulant naming routed to D-track.
- 3.2 **W2** Code→spec back-translation — subagent derives "what equations
  does this code actually solve" per module; diff against PHYSICS.md.
- 3.3 Claims-ledger cross-check — **DONE (by the parallel V&V campaign)**:
  cx/vv-trace docs/qa/VV_TRACEABILITY.md — 27 VALIDATED / 4 VERIFIED-ONLY /
  3 BENCH-PENDING / 8 SPEC-ONLY / 1 UNSAFE-CLAIM (T17-03 stress/LES/IBM,
  consistent with our ANOM-P4-001 routing). Merge via task #12.
- 3.4 Review panel on fix designs — **GATE (P4-001/010 core landing)** —
  2-3 agent judge panel on the fix design before gate-running.

## Axis 4 — Literature & community bad-practice sweep (web)

- 4.1 LBM pitfall compendium — **DONE** (commit f1d9eb2,
  `docs/qa/pitfall-checklist.md`): 27 pitfalls mapped, 16 COVERED /
  8 PARTIAL / 1 GAP (checkerboard modes) / 1 N/A; 9 kill-cases queued in
  priority order → fold into lane 1.7 W2 orders (top 4: checkerboard-mode
  decay, co-moving-frame zero wall force (Wen-2014 decider), WALE
  pure-shear nu_t field check, single-delta streaming + OPP[] table test).
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
- 5.3 **W2** Visual anomaly trawl — **HARNESS LANDED 2026-07-07**:
  `scripts/qa/visual_trawl.py` scans existing PNG/VTK outputs and gallery
  `vtk_field` files; guide: `docs/qa/visual-trawl-guide.md`. Operator matrix
  run + qa-viewer behavior review remains the W2 campaign task.
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
- 6.3 GPU absolute physics — **PARTIAL**: cx/vv-backend (parallel campaign)
  landed D2Q9 GPU-direct TGV (order 1.883) + pressure-channel sentinels;
  remaining: D3Q19 direct, GPU Ghia cavity, GPU conservation — GATE (R2-D3)
  for the full matrix.

## Axis 7 — Static & structural verification

- 7.1 **W1** Ban-list grep sweep — case-identity branches, transported-
  quantity clamps, bare calibrated literals; recurring CI-able check.
- 7.2 Numeric-literal provenance scan — **DONE** (structural-sweep):
  no new undocumented physics constants; P2 comment order queued.
- 7.3 **W3** Unsafe/precision audit — SIMD unsafe UB; f32↔f64 cast audit;
  Kahan-summation needs.
- 7.4 Docs/spec consistency — **DONE** (structural-sweep): 14/14 verified
  bands match exactly, zero drift; agent's missing-test rows PM-refuted.

## Axis 8 — Interface/tooling V&V (product surface)

- 8.1 **W2** Scenario schema round-trip + fuzz.
- 8.2 **W3** WASM/GUI parity vs native — harness landed 2026-07-07:
  native exporter `wasm_native_parity_export`, comparator
  `scripts/qa/wasm_parity_check.py`, workflow guide
  `docs/qa/wasm-native-parity-guide.md`; operator WASM build/snapshot run
  remains pending.
- 8.3 **W3** CLI/MCP contract tests (presets/gallery/VTK/manifest schema).
- 8.4 **W3** Units V&V extension — SI↔lattice across all constructors.

## Axis 9 — Real-world scenario gauntlet (experiment-driven, visualized)

User directive 2026-07-06: V&V must also exercise REAL multiphysics
scenarios, not only analytic gates — run them, visualize them, behavior-
review them, and be honest about capability limits. Every case: at least one
quantitative anchor + behavior anchors (pattern/sign/trend), mandatory
visual artifacts (PNG series / viewer), behavior-validity review by the PM,
findings to anomaly-log. Capability tags: RUN-NOW / AFTER-FIX (named) /
GATED (named track).

- 9.1 **Sparger bubble injection** — **GATED (ANOM-P4-020)**: SCMP mass
  sources have no phase identity — no bubbles form at any rate (visual
  evidence on cx/vv-sparger). Needs MCMP component sources or MF-γ gas
  inflow (VR-STR-02). Original plan text: at the
  Shan-Chen achievable density ratio (~15; REAL air-water 1000:1 is GATED
  on MF-γ — say so in every report). 2D column, SC liquid pool + gas
  injected through a bottom masked patch / volume source; observables:
  bubble detachment period vs injection rate (monotone anchor), rise
  velocity vs the SC-consistent buoyancy-drag balance, Laplace-consistent
  bubble pressure, coalescence behavior; spurious-current context printed.
- 9.2 **Half-filled rotating drum** — RUN-NOW with caveats: SC liquid/vapor
  pool (diffuse interface, ratio ~15) + per-mass gravity + rotating THIN
  SHELL via penalization (thin features sit in the measured stable margin;
  ANOM-P4-010 caveat printed; rerun after cx/fix-p4-010 lands). Observables:
  static pool level, surface inclination vs rotation rate (solid-body limit
  at low Ro), onset of recirculation/cascading analog, vortical structures.
  HONESTY: moderate-Re unsteady flow within the grid-Re envelope — label it
  "unsteady vortical flow", never "turbulence"; free surface is a diffuse
  low-ratio interface, not a VOF surface (MF-γ).
- 9.3 Fully-filled drum spin-up vs analytic spin-up timescale — AFTER-FIX
  (P4-010 or P4-001; the rotating boundary needs a healthy path).
- 9.4 Sedimentation basin: CR-3 particles + gravity + weak crossflow,
  deposition-map visualization vs settling-length estimate — RUN-NOW.
- 9.5 Stirred-tank free-surface vortex — GATED (MF-γ free surface).
- 9.6 Single-bubble Grace-diagram + Hysing benchmark — GATED (MF-γ;
  already the MF-γ headline acceptance).
- 9.7 Lateral-oscillation sloshing — RUN-NOW as SC low-ratio analog
  (resonance-frequency anchor vs shallow-water estimate); real after MF-γ.
- 9.8 Taylor-Couette wavy-vortex transition (3D, LES) — READY-heavy
  (schedule against compute budget; onset Ta_c anchor).
- 9.9 Bubble-plume in a tank (SC vapor source composition) — RUN-NOW after
  9.1 validates the source+SC composition (compose rules are exactly the
  interaction-matrix lane 5.1 territory).

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
