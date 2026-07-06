# V&V Master Plan — every angle, in parallel (2026-07-06)

Owner: FSI V&V session (PM). Execution is delegated to sub-agents: codex CLI
orders (worktree-isolated, adversarial where applicable), Claude subagents
(research / review / synthesis), and background scripts. This plan is the
routing table; per-lane findings flow into `docs/qa/anomaly-log.md` and
core/demo/spec routing per the CLAUDE.md V&V rule. Machine budget: ≤5-8
compile-bound orders at once — lanes are scheduled in waves, not all at once;
"parallel" means every INDEPENDENT lane proceeds without waiting on another
lane's result unless a dependency is stated.

Principles (frozen): physical rigor prime directive; tests adversarial and
implementation-separated; derive before blaming the engine; every band has a
behavior anchor; every experiment leaves a visual artifact; bands frozen on
(C) terms are not authoritative.

## Status legend

RUN = running now · DONE = delivered (regression-pinned) · W1/W2/W3 = wave ·
GATE = blocked on a named dependency

---

## Axis 1 — Adversarial analytic accuracy audits (skill: lbmflow-accuracy-audit)

| Lane | Method | Executor | Deliverable / gate | Status |
|---|---|---|---|---|
| 1.1 Momentum-exchange probe + walls | A1-A5 audit, exact ledgers | codex (cx/audit-probe) | 5 probes, merged | DONE |
| 1.2 Rotating IBM | torque/analytic, stability sweep | codex (cx/audit-ibm) | 8 probes = ANOM-P4-001 gate | GATE (core fix) |
| 1.3 Particles CR-3 | integrator pinning, sampler contracts | codex (cx/audit-particles) | 6 probes, merged | DONE |
| 1.4 Sources + patches CR-1/2 | bookkeeping identities, dipole, node exactness | codex (cx/audit-sources) | 5+1 probes, merged | DONE |
| 1.5 Cumulant viscosity correction | h² extrapolation, tau fingerprint | codex (cx/audit-cumulant) | e2 = ANOM-P4-008 gate | GATE (core fix) |
| 1.6 Penalized rotor | annular-Couette torque, contracts | codex (cx/audit-rotor) | F1-F5 = ANOM-P4-010 gate | GATE (core fix) |
| 1.7 Remaining subsystems: Shan-Chen force accuracy (beyond T11 bands: pressure-tensor form, spurious-current tau-dependence), WALE strain reconstruction vs analytic strain (extend g2), convective outflow reflection coefficient curve, FP16/deviation storage error model, Bouzidi MOVING walls (if supported) | same A1-A5 loop | codex, one order per subsystem | new audit files + triage | W2 |
| 1.8 Tracer-limit particle audit (near-neutral buoyancy, Phase B) | v→u relaxation in nonuniform flow; BBO-term validity floor | codex | audit rows + validity-domain doc | GATE (Phase B spec) |

## Axis 2 — Meta-V&V: does the net catch fish? (NEW)

| Lane | Method | Executor | Deliverable / gate | Status |
|---|---|---|---|---|
| 2.1 **Mutation testing of the validation suite** | Inject deliberate physics bugs (sign flip in Guo source, factor-2 in MW term, w_q swap, off-by-half wall position, dropped corner link, tau→omega typo, feq u² coefficient) one at a time into a scratch worktree; run the default suite; record CAUGHT/MISSED per mutant. Every MISSED mutant = a coverage hole → new audit row. | codex order per mutant batch (scripted harness), PM triages misses | `docs/qa/mutation-coverage.md` matrix; target: 0 silent mutants in physics kernels | **W1** |
| 2.2 Band-vacuity scan | Grep all asserts for band/measured ratios ≥ 100x (the rotating_ibm L2<0.95 class); flag near-vacuous gates for retightening | Claude subagent (read-only) + PM ruling | list + retighten orders | **W1** |
| 2.3 Behavior-anchor coverage | For every VALIDATION.md band, check a pattern/sign/monotonicity anchor exists (two-layer rule); list band-only gates | Claude subagent | gap list → codex retrofit order | W2 |

## Axis 3 — Cross-intelligence review (NEW)

| Lane | Method | Executor | Deliverable | Status |
|---|---|---|---|---|
| 3.1 **codex independent physics review** | codex reads kernels.rs/boundary/forcing/multiphase COLD (no audit list, no anomaly log) and re-derives the physics from Chapman-Enskog + cited papers; reports every place the code disagrees with its own derivation. Different intelligence, different blind spots. | codex (read-only order, own worktree) | findings memo → PM triage | **W1** |
| 3.2 Code→spec back-translation | Subagent derives "what equations does this code actually solve" per module and diffs against PHYSICS.md claims; mismatches are doc bugs or physics bugs | Claude general agent | diff table | W2 |
| 3.3 Claims-ledger cross-check | Every PHYSICS.md / whitepaper claim mapped to the test that proves it; unproven claims flagged (ties into the sales-paper claims ledger) | Claude subagent | claim→test map, red rows | W2 |
| 3.4 Review panel on fix designs | When core lands the P4-001/010 family fix, run a 2-3 agent judge panel on the design before gate-running | Claude agents | panel verdict | GATE (core fix) |

## Axis 4 — Literature & community bad-practice sweep (NEW, web)

| Lane | Method | Executor | Deliverable | Status |
|---|---|---|---|---|
| 4.1 **LBM pitfall compendium** | Web-research the known failure catalog: forcing-scheme inconsistencies, non-equilibrium initialization, corner/edge BC handling, checkerboard modes, tau→0.5 oscillations, viscosity-dependent BB slip, momentum-exchange Galilean correction (Wen et al.), interpolated-BB mass leakage, Shan-Chen thermodynamic inconsistency, LES wall behavior, cumulant Galilean claims vs measurements, Ma² contamination of benchmarks, staircase artifacts. Map EVERY item → {covered by test X / gap / N-A with reason}. | Claude research agent (WebSearch/WebFetch) | `docs/qa/pitfall-checklist.md` with coverage column; gaps become audit rows | **W1** |
| 4.2 Benchmark-matrix expansion | Survey standard benchmarks we do NOT run: ten Cate sphere sedimentation (experimental data), DKT two-particle, Segré-Silberberg (two-way gate), Turek-Hron FSI (elastic — capability gap), Hysing bubble benchmark (MF-γ gate), Womersley pulsatile, oscillating-cylinder KC. Rank by (capability exists now / needs planned track / out of scope) | Claude research agent | ranked benchmark backlog with reference data sources | **W1** |
| 4.3 V&V methodology audit | Check our practice against ASME V&V 20 / Oberkampf-Roy vocabulary (verification vs validation vs calibration separation; solution verification = grid-convergence reporting discipline) | Claude research agent | gap memo (e.g. GCI reporting) | W2 |

## Axis 5 — Systematic experiment matrices + visualization (NEW)

| Lane | Method | Executor | Deliverable | Status |
|---|---|---|---|---|
| 5.1 **Feature-interaction conservation matrix** | Pairwise (later triple) combinations of {gravity, SC multiphase, rotor, particles, sources, patches, LES, f32, 3D} on small grids; per cell of the matrix: mass/momentum ledger closure, NaN watchdog, symmetry spot-check. Composition bugs (force-field overwrite/add rules) live exactly here. | codex order (test harness `tests/interaction_matrix.rs` or script) | matrix report; every red cell → anomaly entry | **W1** |
| 5.2 Randomized property sweep | proptest-style randomized legal configs (BC combos × obstacles × params in envelope): invariants = conservation, boundedness, determinism, mirror equivariance. Seeded, reproducible. | codex order | property test file (default-suite light + heavy fuzz tier) | W2 |
| 5.3 Visual anomaly trawl | Scripted scenario matrix (per preset × parameter grid incl. envelope corners) → PNG/VTK export → sim-anomaly-scan + qa-viewer dashboards; PM does behavior-validity review on flagged cases (patterns, walls, seams) | scripts/qa + Claude agent orchestration; PM reviews visuals | flagged-case gallery + reviews in findings ledger | W2 |
| 5.4 Long-horizon soak | 10⁶-step runs of 3-4 canonical cases (cavity, channel+cylinder, SC droplet, penalized blade rotor post-fix) watching drift ledgers (mass, momentum, energy-like monitors) | background scripts (nice'd) | drift report; extends d5_long_horizon | W2 |
| 5.5 Stability-envelope cartography | Automated (tau, Ma, grid-Re) boundary mapping per collision operator; publish measured envelope vs documented guidance (tune-stability thresholds) | codex order + script | envelope charts (artifact) + doc reconciliation | W3 |
| 5.6 Degenerate-geometry gauntlet | 1-cell gaps, diagonal staircases, obstacles touching rims/patches/seams, zero-thickness features, moving-wall corners | codex adversarial order | gauntlet test file | W2 |

## Axis 6 — External referees

| Lane | Method | Executor | Deliverable | Status |
|---|---|---|---|---|
| 6.1 Cross-solver comparison | Same scenarios in OpenLB (already built at ~/projects/cfd-bench) for cavity/cylinder/TGV; independent-implementation referee for absolute values (not just self-consistency) | codex order in cfd-bench + comparison script | cross-solver report with L2 diffs | W2 |
| 6.2 Spectral mini-referee | 64-line Python spectral NS solver for 2D TGV/Kolmogorov flow as an in-repo independent reference (no external dep) | codex order (scripts/qa) | reference generator + comparison | W3 |
| 6.3 GPU absolute physics (D-3 gap) | TGV order + cavity Ghia directly on GPU (not only CPU-relative) — planned R-Phase 2 item; PM-run evidence | PM (outside sandbox) after R2-D3 lands | T14 extension | GATE (R2) |

## Axis 7 — Static & structural verification

| Lane | Method | Executor | Deliverable | Status |
|---|---|---|---|---|
| 7.1 Ban-list grep sweep (physics-discipline) | Automated: case-identity branches, clamps on transported quantities, bare calibrated literals, silent physical defaults | script + Claude triage | recurring CI-able check | W1 (script), recurring |
| 7.2 Numeric-literal provenance scan | Every float literal in src/ physics paths must have a provenance comment or named constant; list violators | Claude subagent | violator list → doc order | W2 |
| 7.3 Unsafe/precision audit | unsafe blocks (SIMD) UB review; f32↔f64 cast audit; summation-order sensitivity spot checks (Kahan needed anywhere?) | codex read-only order | memo | W3 |
| 7.4 Docs/spec consistency | VALIDATION.md ↔ test files drift check (bands quoted in doc vs asserted in code) | Claude subagent | drift list | W2 |

## Axis 8 — Interface/tooling V&V (product surface)

| Lane | Method | Executor | Deliverable | Status |
|---|---|---|---|---|
| 8.1 Scenario schema round-trip + fuzz | JSON→build→export→re-import identity; malformed-input error quality | codex order | tests in lbm-scenario | W2 |
| 8.2 WASM/GUI parity | Same scenario native vs WASM: field parity bands; mock-engine badge (D-track order covers the badge) | codex order | parity test + CI note | W3 |
| 8.3 CLI/MCP contract tests | presets/gallery/VTK validity (parseable), manifest schema, MCP tool round-trips | codex order | contract test file | W3 |
| 8.4 Units V&V extension | SI→lattice conversions vs hand-computed anchors across all constructors incl. gravity/pressure edge ranges | codex order | extended unit tests | W3 |

## Dependencies & sequencing

- **W1 (dispatch now, parallel)**: 2.1 mutation pilot, 2.2 band-vacuity scan,
  3.1 codex cold review, 4.1 pitfall compendium, 4.2 benchmark survey,
  5.1 interaction matrix, 7.1 ban-list script. (2 compile-bound codex + 3
  read-light agents + 1 script — fits beside the running full tier.)
- **W2** dispatches as W1 lanes return and free slots; W3 after.
- GATE lanes unblock on: core fixes (P4-001/008/010), R-Phase 2 D-3, Phase B
  spec, MF-γ.
- Every lane's failures route per the standing rule (core / demo / spec) and
  are logged in anomaly-log with ANOM ids. Every experiment lane ships visual
  artifacts; PM does the looking.

## Reporting

Single dashboard artifact (pass-4 board, extended per wave):
https://claude.ai/code/artifact/ced42aaf-5832-4dbd-8954-be9900140519
plus per-lane deliverable files under docs/qa/. Cross-session routing:
core-engine session (defects), D-track PM (demo/spec), user (dashboard +
summaries).
