# Ad-hoc physics inventory — 2026-07-06 (V&V session, ORDER Task 1)

Deliverable for `ORDER_vv_adhoc_extermination.md` Task 1. Read-only sweep of
`crates/` (core, scenario, cli incl. examples) and `web/`; three parallel
sub-sweeps (core / cli+examples / scenario+wasm+web). Classification:
**(A)** resolved physics · **(B)** literature-backed documented closure ·
**(C)** ad-hoc, to kill. Kill plans are proposals — PM triage (Task 2 gate)
decides before any code change is ordered.

Priority scale: P0 = distorts current headline results silently; P1 = wrong
physics in a supported path, bounded blast radius; P2 = hygiene/documentation.

---

## 1. `crates/lbm-cli/examples/dispersed_seeding/` — the (C) cluster

The example's Lagrangian step loop and `sample_tray` layer analytic closures
and case-identity branches ON TOP of the resolved LBM tray field. These set
the magnitude and spatial pattern of the headline deposition maps (T18.4
anchors ran through this path). All confirmed at exact locations.

| # | Item | Location | Class | Evidence | Kill plan | Priority |
|---|---|---|---|---|---|---|
| 1.1 | `harshness` composite switch (0.01, 0.15, 20.0) gating 3 regimes | particles.rs:172-174 | **C** | Case-identity branch; constants undocumented anywhere | Delete the branch; regime differences must emerge from the resolved jet/agitation flow | **P0** |
| 1.2 | `wall_jet_len` 0.42 vs 0.10 of tray width by harshness | particles.rs:254-257 | **C** | Case-identity constants, undocumented | Superseded by resolving the wall jet with core CR-1 `Jet` + CR-2 patches (already landed); delete the analytic wall-jet entirely | **P0** |
| 1.3 | Analytic jet-gaussian + wall-jet closure layered on LBM field (0.25, 0.75, 0.35 coefficients; exp(−r²/2σ²), exp(−r/L)) | particles.rs:265-271 | **C** | Undocumented closure duplicating what the resolved field should provide | Same as 1.2: sample the RESOLVED tray velocity (core interpolation) instead of the superposed closure. If tray resolution cannot support the near-field jet, that is a CORE CAPABILITY GAP → routing package | **P0** |
| 1.4 | Lateral dispersion constant 2.5e-5 m²/s (`gentle_k`) | particles.rs:175 | **C** | SPEC_FINDINGS admits calibrated to the acceptance gate | Replace with a cited turbulent/laminar dispersion closure with validity domain, or resolve dispersion (LES at tray scale = capability gap → route); until then the constant must NOT drive accepted results | **P0** |
| 1.5 | Side-wall position clamps trapping particles at the rim | particles.rs:218-222 | **C** | The behavior-validity origin finding: edge-ring deposition pattern is MADE by this clamp; a clamp is not a wall BC | Use the core CR-3 solid-contact model (substep crossing + restitution/deposit) for tray walls instead of clamping; the edge ring must disappear or be physically explained | **P0** |
| 1.6 | Deterministic `sin()` pseudo-agitation kicks (a·ω²·sin(ωt); harsh-mode 0.15·jet_sigma·sin(11.0·t + d·1e6)) | particles.rs:180, 210-213 | **C** | Decorative term; frequency keyed to particle diameter is pure invention | Agitation = oscillating body force / moving boundary resolved by the core (translational pattern is frozen in DISPERSED_DEPOSITION §3.2); delete the kicks | **P0** |
| 1.7 | Ejection jitter 0.35·jet_sigma; jet_sigma floor 1.25·dx | particles.rs:86, 93-96 | **C** | Calibrated constants, no source | Nozzle exit distribution is a spec input (document in DISPERSED_DEPOSITION) or derive from nozzle geometry; floor only as numerical guard with comment | P1 |
| 1.8 | Reservoir extraction heuristics: band (0.08, 0.10, 2000.0, 1.5), `settled_bonus` size-gate, score weight 0.18 | reservoir.rs:23-36 | **C** | Ad-hoc scoring, undocumented | Replace with resolved settling column (core particles + Stokes/SN settling already validated in T18.3) sampled at the withdrawal port | P1 |
| 1.9 | Mystery reservoir force [0,0,−1e-7] | main.rs:173 | **C** | Purpose unstated; 8 orders below g | Delete or justify with a comment + doc entry | P2 |
| 1.10 | Fixed spin-up step counts 40 / 90 | main.rs:40, 48 | **C** | No justification | Replace with the documented steady/quasi-steady criterion or record rationale | P2 |
| 1.11 | dt pickers 0.16·dx/u_max, 0.012·dx²/ν | protocol.rs:142-148 | **B-gap** | Conservative margins consistent with the frozen envelope (Ma soft 0.15, τ floor) but underivation not written | Document derivation from the frozen thresholds inline + DISPERSED_DEPOSITION note | P2 |

## 2. `crates/lbm-cli/examples/stirred_tank_3d.rs`

| # | Item | Location | Class | Evidence | Kill plan | Priority |
|---|---|---|---|---|---|---|
| 2.1 | Penalization α=0.32 + load-bearing per-cell force cap f_cap=0.25·u_tip | stirred_tank_3d.rs:240-242, 349-351 | **C** | PHYSICS.md (rotor penalization entry) already marks this scheme superseded: "stability comes from the implicit-style Guo force balance, not from clipping"; the example still runs the old scheme | Migrate the example to the landed core penalization; delete α/f_cap | **P1** |
| 2.2 | Rushton geometry ratios (0.66, 0.22, 0.12, 0.30, 1.2, 0.2, clearance 0.35) | stirred_tank_3d.rs:144-153 | **B-gap** | Geometry choices ≈ standard Rushton proportions but source not cited; two mystery lattice constants (1.2, shaft floor 1.5) | Cite the standard tank configuration; state lattice-resolution floors as numerical guards | P2 |
| 2.3 | Ma_tip ≤ 0.3 abort, spin-up ramp 1500 | stirred_tank_3d.rs:241, runner.rs:286 | **B** | Frozen in PHYSICS.md stability envelope | none | — |

## 3. `crates/lbm-core/` — clean, one flagged closure

Core kernels (equilibrium, TRT+Guo, Zou-He with in-code 60-line derivation,
outflow/convective, Bouzidi with half-way degeneracy test, WALE Cw=0.325 with
null-property gates, IBM Uhlmann/Wang with characterization freeze, SN drag
with Re<800 validity guard, per-mass gravity, CR-1 equilibrium-shaped
injection, CR-3 contact/tunneling guard, probe momentum-exchange): **(A)/(B)
compliant** — formulas literature-standard, PHYSICS.md entries + validation
tests exist. Numerical-safety floors (ρ.max(1e-30), interpolation bracket
clamps) documented as guards.

| # | Item | Location | Class | Evidence | Action | Priority |
|---|---|---|---|---|---|---|
| 3.1 | Cumulant D3Q19 viscosity-offset correction: `omega_eff = omega_shear·(1 + 0.0025 − 0.16·u²)`, `.min(2.0)`, D3Q19-only branch | kernels.rs:415-416 | **C — CONFIRMED by measurement (Order E verdict, ANOM-P4-008, commit 5eae598)** | 3-point h² extrapolation: uncorrected operator has ~zero continuum bias (D3Q27 control intercept −6.9e-5); the offset's footprint −0.0025·2/(2−ω) = −2.333e-2 matches the D3Q19 N=32 intercept −2.3275e-2 to 0.2% — it is a resolution-point calibration that cancels O(h²) spatial error at exactly N=32 and corrupts every other resolution. PHYSICS.md documents it (Geier et al. refs) BUT the 0.0025 was fitted to the TGV3D decay — a constant calibrated against the acceptance observable, and the O(u²) defect coefficient 0.16 should be derivable analytically for the lattice rather than fitted | (a) inline comment citing PHYSICS.md + refs; (b) **adversarial audit row queued (accuracy-audit Order E)**: ν_eff(u, N) functional-form sweep — the correction must cancel the defect across u AND resolutions, not only at the calibration point; VERDICT: resolution-dependent → (C) — CLOSED 2026-07-06: offset removed by core; the −0.16·u² companion measured (B) by ablation (final inventory 3.1 verdict, ANOM-P4-008 closed) | **P0 CLOSED** |

## 4. `crates/lbm-scenario`, `crates/lbm-wasm`, `web/`

| # | Item | Location | Class | Evidence | Action | Priority |
|---|---|---|---|---|---|---|
| 4.1 | Unit conversions (gravity/pressure/force/τ↔ν/resolution constructors) | lbm-scenario/src/units.rs | **A** | Exact dimensional analysis; 3-constructor equivalence + anchor tests | none | — |
| 4.2 | Inlet parabola, obstacle rasterization, maskToRects | lbm-scenario/src/lib.rs, web/src/scenario.ts | **A** | Exact forms/lossless | none | — |
| 4.3 | Mock engine synthetic flow (~10 look-tuned constants, BGK-noise "flourish") | web/src/engine/mock.ts:120-400 | **C-quarantined** | Candidly self-documented as non-physics UI fallback; never used for real runs | Runtime `console.warn` + visible UI badge when mock is active; forbid presenting mock output as physics | P2 |
| 4.4 | Vorticity display scale 0.7 + EMA range smoothing | web/src/render.ts:120-148 | display-only | Colorbar normalization only; vorticity field itself is exact central-difference | inline comment marking as UI scaling | P2 |

## 5. Summary and routing

- **(C) to kill: 11 items**, ALL in `dispersed_seeding` + the deprecated
  stirred-tank penalization pathway + quarantined mock. None in lbm-core.
- **P0 cluster (1.1-1.6)**: one coherent kill = rebuild the dispersed_seeding
  step loop on the landed core capabilities (CR-1 Jet, CR-2 patches, CR-3
  particles/contact/deposition) and delete the closure layer. Anything the
  resolved field cannot supply at current resolution (near-field jet,
  turbulent dispersion) becomes a CORE CAPABILITY GAP routing package, not a
  constant.
- **Verification queued**: cumulant offset functional-form audit (3.1) joins
  the accuracy-audit fan-out as Order E.
- T18.4 trend anchors and the P1.1 bands were measured through the (C) layer;
  after the kill they must be re-measured and re-frozen (expected to move —
  that is the point).

Adversarial-accuracy audit running in parallel (same session): codex orders
A (momentum-exchange probe + walls) and B (rotating IBM) dispatched; C
(particles) and D (sources/patches) queued. Findings will be triaged into
docs/qa/anomaly-log.md per protocol.

---

## 6. PM triage (D-track PM, 2026-07-06)

- **1.1–1.11 APPROVED as ONE kill order** (`cx/kill-deposition-closures`):
  rebuild dispersed_seeding on co-evolved resolved physics — solver stepped
  WITH the particles (point-sampled live field), closures/clamps/kicks/jitter
  deleted, agitation = oscillating body force on the fluid + matched
  density-weighted pseudo-force on particles (same frozen §3.2 pattern),
  reservoir extraction from the VALIDATED core settling column, spin-up by
  criterion. 1.4 decision: the dispersion constant is DELETED, not replaced —
  if the resolved field under-disperses, that is a finding (literature
  closure with the four artifacts, or capability-gap routing), not a knob.
  All P1.1/T18.4 bands measured through the (C) layer are INVALIDATED; the
  order reports new values + visual artifacts, PM re-freezes.
- **2.1–2.2 APPROVED as a second order** (`cx/stirred-penalization-migrate`):
  migrate the example to the landed core penalization, delete α/f_cap, cite
  the standard Rushton proportions.
- **3.1 ENDORSED**: accuracy-audit Order E (functional-form sweep across u
  AND N) decides (B) vs (C); core session notified. No comment-only fix
  before the audit verdict.
- **4.3 APPROVED as a light order** (`cx/mock-engine-warning`): runtime
  console.warn + visible UI badge when the mock engine is active.
- **4.4, 1.11, 2.2 doc items**: folded into the respective orders above.
- Sequencing: kill order runs at current sample parameters (100 µm);
  the 20 µm real-problem reparameterization (Phase B: exact ballistic-settle
  switch, budget redesign) follows as a separate order once the bead
  material/density is confirmed by the user.
