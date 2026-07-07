# Physics Anomaly Sweep — log

Living register of physics anomalies. Rows drop when superseded by a landed
fix + regression pin (or refuted). Test-cited stubs stay: they document
current-wrong-value pins that must flip in the same commit as the fix.

Taxonomy: **S0** silently-wrong · **S1** divergence-leak · **S2**
below-expected accuracy · **S3** minor. Historical triage narrative is in
git history (2e121c8 band-vacuity scan, 714da6a cold-review triage, pass 1–4
context in commits 82946d3 / cd3999f / 20d0e10).

---

## RESOLVED — fixed and regression-pinned

### ANOM-P4-001 — time-stepped direct-forcing IBM diverges in default config
- Severity: **S1**.
- Status: **RESOLVED 2026-07-07**.
- Fix: marker force sizing now targets the realized full Guo force-field
  impulse `F/rho` instead of the half-force diagnostic increment, and dense
  marker overlap uses the row-sum mobility of the interpolation-spreading Gram
  operator instead of the self-mobility `M_jj` alone.
- Gate: `cargo test -p lbm-core --release --test accuracy_audit_ibm -- --nocapture`
  passed with default `relaxation=1.0`; `rotating_ibm.rs` quantitative bands
  were tightened or marked superseded by the audit. Measured highlights:
  B1 torque ratio `1.06598`, B2 spread `1.641e-3`, B5 kernel/relaxation
  torques mutually within 5%, B6 torque antisymmetry `1.041e-16`, B8
  Taylor-Couette `L2_rel=4.522e-2`.
- Follow-up boundary: ANOM-P4-010 remains open. Volume penalization has no
  marker interpolation-spreading Gram operator and still needs its own derived
  implicit/relaxation treatment.

## OPEN — core-engine fix pending

### ANOM-P4-008 — cumulant D3Q19 "+0.0025 viscosity offset" is a resolution-point calibration (verdict C)
- Severity: **S2** leaning S0 (silent tau-dependent bias at every N except
  calibration point).
- Gate: `cx/audit-cumulant` e2 canary (|a| ≤ 2e-3 after removal).
- Measured (N ∈ {24,32,48}, diffusive u0): D3Q19 defects fit
  d = a + b/N² exactly with a = −2.3275e-2, b = +23.22. The intercept matches
  the offset's own nu-space footprint −0.0025·2/(2−ω) at ω=1.7857 to 99.8%.
  D3Q27 control (offset=0): a27 ≈ 0. Confirms the uncorrected cumulant has
  ~zero resolution-independent bias — the offset cancels ordinary O(h²) at
  one resolution only.
- Recommended: remove offset; re-freeze TGV3D acceptance with h²-intercept
  criterion; ablate the −0.16·u² term separately.

### ANOM-P4-010 — compat volume penalization diverges for solid disc at Re 0.09
- Severity: **S1**. Same family as P4-001 (target-the-half-force sizing).
- Gate: `cx/audit-rotor` F1/F2/F3 (F4 cross-path after both fixes).
- Config: 80×80 periodic compat, TRT Λ=3/16, nu=1/6, solid disc
  n_blades=2/r_hub=0/thickness=2R, chi=1, Ω=1.5e-4.
- Observed: step-1 torque correct; then sign-flipping growth
  (u* → 2u_t − u* per step); density catastrophe by step ~120, e131 by
  step 400. Thin blades are marginally damped by streaming exchange with
  neighbors — the old stirred example's f_cap clamp was load-bearing here.
- Disposition: one core fix should cover both family members
  (full-step-consistent force sizing).

---

## OPEN — audit-side revision

### ANOM-P4-014 — Jurin SC coefficient is linear but high after wetted-wall audit
- Severity: **S2 characterization**.
- Gate: `validation_multiphase_hard.rs::val_mphard_i1_jurin_capillary_rise_zero_parameter`
  rev 5 freezes linearity plus measured slope while the coefficient question
  remains open.
- Measured (2026-07-07, rev 5): `h` vs two-wetted-channel inverse-width
  contrast has slope 1312.71105871, intercept -1.38619080, R² 0.99998340.
  Flat-wall Jurin theory from T11/T11c constants gives 852.22687953, so the
  measured coefficient is 1.54033050× theory.
- Wetted-wall audit: `ShanChen::with_wall_rho` applies to all solid
  neighbours, including the domain rim and the inserted slot walls, because
  `compat::multiphase::ShanChen::update_force` branches on
  `sim.solid_field()[j]` and adds `psi_wall`.
- Gap-24 measured meniscus angles: slot 66.416148°, outside 65.505811°
  (left 65.504972°, right 65.506650°), so rim and inserted walls do not show
  a distinct wall-density class in this setup.
- Behavior review: the slot rise/depression pattern is monotone with gap and
  flips sign for the dry wall-density case. The dominant mechanism is
  capillary rise in finite outside channels connected through the reservoir.
  The unresolved part is coefficient magnitude, plausibly SC diffuse-meniscus
  curvature vs flat-wall contact-angle calibration. Artifacts checked:
  generated density maps at `target/vv_jurin/jurin_gap16_wallrho1.000.pgm`,
  `target/vv_jurin/jurin_gap24_wallrho1.000.pgm`,
  `target/vv_jurin/jurin_gap32_wallrho1.000.pgm`, and
  `target/vv_jurin/jurin_gap24_wallrho0.600.pgm`; no disconnected liquid
  column was observed by the connectivity diagnostic.
- Disposition: keep as characterization until a curvature-aware SC capillary
  coefficient or a matched meniscus/flat-wall calibration gate is derived.

### ANOM-P4-007 — cumulant viscosity-offset audit design confounded
- Severity: audit design.
- Standing: orientation consistency PASSED (spread 2.2e-10); calibration
  residual −5.9e-4. E1 u0-sweep at fixed N conflates O(Ma²) with cubic
  defect; E2 N=24 band ignored O(h²) spatial-error floor.
- Revised probes queued: tau-sweep, N-sweep at fixed Ma for D3Q27,
  spatial-error-modeled bands. Superseded in verdict by ANOM-P4-008 (offset
  is a calibration hack); this row stays until the revised audit lands.

---

## Test-cited stubs (must flip in the fix commit)

- **ANOM-P4-004** — `accuracy_audit_particles.rs:14`. Test-side v0=1e-10
  fix (Stokes limit). No engine work.
- **ANOM-P4-005** — `accuracy_audit_sources.rs:239`. Semantics pin: q_lu is
  REGION TOTAL (not per-cell). Doc gap queued in DISPERSED_DEPOSITION.md §5.
- **ANOM-P4-006** — `accuracy_audit_sources.rs:431`. Patch BC nodes = face
  layer (exact); adjacent interior carries developed flow.
- **ANOM-P4-009** — `crates/lbm-core/src/compat/rotor.rs:163,208`. Two
  contracts: hub region r<r_hub is a HOLE (not solid); `update_force` ADDS
  into the field — caller must `clear_force_field()`.

---

## Runtime/tooling proposals (open, non-core-engine)

- **ANOM-P1-001** — S2 tooling: mass/momentum drift is not first-class in
  the manifest. Proposal: periodic diagnostics series
  (step, totalMass, totalMomentum[, maxSpeed]) at `checkEvery`.
- **ANOM-P1-003** — S3 monitoring: runtime `maxSpeed` can exceed the
  compressibility advisory (0.15 / 0.3) with no signal. Proposal: promote
  validate thresholds to run-end runtime check on `maxSpeed`.
- **A1 (pass 1)** — S0 monitoring: out-of-envelope stirred runs (Ma>0.3 or
  grid-Re≫15) stay bounded and report STABLE with no runtime signal.
  Proposal: echo `max_Ma` and `grid_Re` into the manifest + runtime warn.

---

## Dropped (resolved / superseded / refuted / test-side landed)

Retained pointers only; details in git history:

- **ANOM-P2-002** rotor blade mirror-arms (odd blade counts) — FIXED
  in `qa/mf-integration` (along-blade sign check).
- **ANOM-P4-002 / P4-003** — test-side, fixed in-worktree (probe A2 sign,
  A1/A5 band floors). Commit 2024a52.
- **ANOM-P4-011** — cold-review F19/F20 (Bouzidi probe-force sign) REFUTED
  by derivation (two modules use opposite q conventions consistently).
  Residues: doc comment on convention; W2 Bouzidi mixed-qd kill-case queued.
- **ANOM-DRY-001** — S3 test-side (Bouzidi audit dry run): convergence-fit
  x-axis reversal (`width` vs `1/width`). Fixed in-file with inline note.
- **ANOM-P1-002** — S3 spec: T4 profile band 2e-3 is calibrated to the
  frozen ν=0.02; doc footnote proposal only. Pin ν=0.02 in matrix.py in the
  meantime.
- Cold-review S3 doc routing: `CollisionKind::Cumulant` implements a
  central-moment operator + velocity-dependent relaxation, not full Geier
  cumulant; rename or implement true cumulants (D-track owns the routing).
- Particle SN validity clamp at Re_p=800 is silent; add
  debug_assert→warn/documented saturation at API surface.
- **ANOM-P2-001** — FIXED in R2-C v2 for equivalent uniform per-cell forcing.
  Exactly uniform force-field installs/clears, plus gravity on all-fluid
  domains, refresh moments through the same Guo half-force definition before
  collision as uniform force. Regression:
  `crates/lbm-core/tests/accuracy_audit.rs::uniform_force_impulse_matches_force_field_anom_p2_001`.

---

## Pass 5 — 2026-07-07, hard-multiphase first measurements
(cx/mp-hard rev 2; the diagnostic protocol REFUTED all four test-side
hypotheses — these are real finding candidates)

- **ANOM-P4-014 — Jurin capillary rise: gap-dependent collapse steeper than
  1/w** — S2 candidate, UNDER INVESTIGATION (visuals required). gap16:
  +3.0% (good); gap24: −22.5%; gap32: −67% vs h = 2σcosθ/(Δρgw) with
  T11/T11c inputs; θ_slot ≈ 66-68° stable; h·w falls 878 → 660 → 280
  instead of constant. Known mechanisms predict the OPPOSITE trend. Field
  dumps ordered (rev 3).
- **ANOM-P4-015 — MCMP interfacial-wave frequency +29% at kW ≈ 0.5** — S3
  test-design: bookkeeping fully verified; sharp-interface dispersion
  applied at kW ≈ 0.5 where O((kW)²) ~ 25%. Rev 3: k-sweep (kW ≤ 0.25)
  must recover theory as kW → 0; σ(k) trend = diffuse-interface
  characterization.
- **ANOM-P4-016 — MCMP RT-cutoff seeding NaN** — OPEN: max|u| = 0.3 by
  step 10, NaN < step 1000, at BOTH stable and unstable modes, while T12
  runs the same class healthily. Rev 3 diffs the initialization against
  validation_rt.rs before any engine claim.
- **ANOM-P4-017 — Taylor-Culick: scaling law holds, prefactor 0.49** — S2
  characterization candidate: v ∝ h^(-1/2) confirmed (slope −0.5, r² 0.96,
  h ∈ {16,24,32}) but v/v_TC plateaus at 0.49-0.54 (20h retraction, two
  interfaces confirmed, vapor drag ~6%). Candidate mechanism: Laplace-σ
  (static, T11) vs MECHANICAL σ (momentum flux) discrepancy of the SC
  pressure tensor — the lane-1.7 "SC pressure-tensor form" audit row and
  this measurement now referee each other. If confirmed: documented SC
  validity limit (PHYSICS.md), core FYI.

### Pass 5 addendum (cx/mp-dynamics rev 2) — the SC dynamic-limitations cluster

- **ANOM-P4-015 UPDATED**: H1 (SCMP capillary wave, damped-oscillator fit)
  shows omega +38% / gamma 2.1x at Q ~ 0.8 — same direction as I2's MCMP
  +29% at kW = 0.5. CROSS-MODEL consistency: diffuse-interface sigma(k)
  stiffening at kW ≳ 0.25 is a characterization of both SC variants, not a
  bookkeeping error. kW-sweep (I rev 3) remains the decider.
- **ANOM-P4-018 — SC near-wall shear artifact (H3)** — S2 candidate: in the
  stratified two-layer Couette the LIQUID side carries l2 = 0.56 (vapor
  0.13) with a visible lower-wall density/velocity artifact — the SC wall
  interaction (psi = 0 solid exclusion) generates a near-wall layer that
  distorts the sheared profile. Connects to the T11b known wall-scheme
  limitation; lane 1.7 SC audit takes it.
- **ANOM-P4-019 — SC contact-line immobility (H4 + I1 unified hypothesis)**
  — S2 model-validity candidate: H4 (verified contact + curved meniscus at
  t = 500) STALLS at 5 cells / 6000 steps vs predicted ~38 (Washburn).
  Unified reading of I1: its columns were initialized FILLED and partially
  drain — gap16's +3% is then start-near-answer luck, and gap32's collapse
  is drainage past a pinned contact line. Decisive test (queue into the
  next mp-hard rev): initialize the gap-16 column at two heights (above and
  below prediction) — convergence to DIFFERENT plateaus confirms
  pinning/hysteresis. If confirmed: SC statics (T11/T11c) remain valid;
  any WICKING/SPREADING dynamics claim is out of the SC validity domain
  (strengthens the MF-gamma case). PHYSICS.md entry drafted on
  cx/mp-dynamics (merge with care).

### ANOM-P4-013 VERDICT (Kovasznay, cx/exact-benchmarks 0257e65)
Wet-node Zou-He velocity faces with tangentially-varying profiles carry an
O(h) boundary error dominating the O(h²) interior: three-point ladder
eu = 3.202e-3 / 1.671e-3 / 8.618e-4 (N = 48/96/192), order 0.94-0.96,
r² ≈ 1. Interior order 2 separately proven (T1/T2). Characterization
FROZEN (heavy pin [0.85, 1.15] — a future 2nd-order open-BC must fail it
upward); light band 5e-3 (measured 2.611e-3). Pressure-field validation
passes (slope → 0.991, r² 0.99999). Benchmark-design note: Zou-He
PRESSURE outlets force v = 0 and are NOT valid closures for exact
solutions with tangential outlet velocity (first pass measured the
resulting O(1)-fraction elliptic contamination). S2 characterization;
core FYI: 2nd-order open BC (NEBB / non-equilibrium extrapolation) is the
improvement candidate.

### ANOM-P4-014 VERDICT: TEST-SIDE (reference-datum error, found by the
mandatory visual review). The PGM dumps show healthy connected columns and
menisci at all gaps; but the domain's OUTSIDE channels are themselves
finite-width wetting slots, so the "flat reservoir level" datum is
capillary-elevated by 2σcosθ/(Δρ g w_out). Differential prediction
h_meas = (2σcosθ/Δρg)(1/w_slot − 1/w_out) reconciles the pattern
(gap32: w_out ≈ w_slot ⇒ h_meas → 0, measured 8.8; gap16: w_out ≫ w_slot
⇒ near-ideal, measured +3%). Rev 4: widen the true reservoir (or assert
the differential formula). COROLLARY for ANOM-P4-019: gap-16 liquid ROSE
~55 cells to equilibrium — SC contact lines DO move in a 16-slot; the
"immobility" hypothesis needs re-examination against H4's specific setup
(20-slit + side reservoir) before any model-limit claim.

### ANOM-P4-016 ESCALATED: with the T12-IDENTICAL protocol, the STABLE
orientation (heavy below) diverges by step 400 (max|u| = 4e3,
rho_min = −6.0) while T12's unstable orientation passes CI. Engine-finding
candidate: MCMP + per-component gravity, stable stratification divergence.
Routing package = the rev-3 test + printed trajectory (cx/mp-hard).

2026-07-07 codex characterization after ANOM-P4-022 fix: Item 1 overwrite is
not the P4-016 mechanism. `MultiComponent::update_forces` builds
cross-repulsion and per-component gravity in the same local force array before
adding it to each component's caller-owned field, so there is no prior gravity
field for SC to clobber in the i3 path. Rerun command:
`cargo test --release -p lbm-core --test validation_multiphase_hard val_mphard_i3_rayleigh_taylor_cutoff_light_sign_canary -- --nocapture`.
Result remains red. Mode 3: `p_total_y=-1.44e2` by step 10,
`-1.11e3` by step 100; max-speed locus moves wall-adjacent by step 350
(`heavy:(85,9)`, `max|u|=4.51e-1`), then fails at step 400
(`max|u|=4.013e3` at `heavy:(10,12)`, `rho_min=-6.009` at
`heavy:(160,10)`). Mode 7 stays finite through step 400 but shows the same
large negative bulk momentum and lower-wall high-speed locus by steps
375-400 (`heavy:(1,9..10)`). Artifact density maps:
`target/vv_rt_i3/rt_mode3_step400_heavy.pgm`,
`target/vv_rt_i3/rt_mode3_step400_light.pgm`,
`target/vv_rt_i3/rt_mode7_step400_heavy.pgm`,
`target/vv_rt_i3/rt_mode7_step400_light.pgm`.
Mechanism hypothesis: the hard i3 setup applies gravity to the heavy
component only in a closed box, producing a large nonzero bulk body-force
impulse. The resulting wall-mediated return flow, not the mid-height RT
interface mode, reaches low-Mach-violating velocities and drives the
wall-adjacent density failure. Fixing this requires a derived MCMP buoyancy
forcing model or a spec change to the i3 body-force protocol; no ad-hoc
mean-force subtraction was applied.

### ANOM-P4-008 RESOLVED (core merge 15adfdd; V&V gate verified 2026-07-07)
Offset removed (CPU/SIMD/WGSL); Cumulant→CentralMoment rename; finite-N
band re-frozen to uncorrected measurement; e2 h²-intercept canary flipped
GREEN in the same commit and INDEPENDENTLY re-verified by this session on
main (3 passed / 5 ignored, 68 s). Remaining: the −0.16|u|² term sits
behind CENTRAL_MOMENT_DISABLE_VELOCITY_CORRECTION_FOR_ABLATION (default
active); E1 ablation A/B measurement dispatched (cx/e1-ablation) — verdict
rule: ON−OFF slope difference vs 0.16·W·2/(2−ω) within 30% at two
viscosities ⇒ (B), else (C) candidate. Core also CONFIRMED the P4-001/010
family diagnosis (per-cell gain fix delays 120→1450 steps but a collective
mode remains; ANOM-P2-001 is the deeper defect) — R2-C re-dispatch in
flight on the core side; our gates hold.

### ANOM-P4-008 FULLY CLOSED — the −0.16|u|² term verdict: (B)
Ablation A/B (cx/e1-ablation, docs/qa/e1-ablation-report.md + PNG): ON−OFF
slope difference matches the derived footprint δc = 0.16·(1/4)·2/(2−ω)
within 26.5% (ν=0.02) and 24.9% (ν=0.10); the τ-fingerprint RATIO across
the two viscosities is 3.42 vs 3.50 predicted — the term is a real
Galilean-class correction, retained with measured provenance (inventory
3.1: offset was (C), removed; u²-term is (B), validated). The audit file's
E1 SPEC-GAP stays as documentation (the toggle is a compile-time const, so
the A/B cannot live in one test binary); the report is its evidence.

### ANOM-P4-020 — Axis 9.1 sparger: SCMP interior mass source cannot
express GAS injection — CORE CAPABILITY GAP (by design, not a bug)
Visual evidence (out/vv_sparger_2d/*/density_030000.png, PM-reviewed): no
bubbles at any rate; the pool just densifies (rho_max 2.385 super-liquid
at the source). Mechanism: SC phase identity is local density — a MassFlow
source inside liquid adds liquid; there is no phase channel to inject.
All-rate-identical observables (events=3 detector noise, rise 0.1882,
negative Laplace) are the no-bubble consequence, not an aggregation bug
(mass ledgers differ per rate and pass at 1e-8). Disposition: Axis 9.1
re-tagged GATED — needs MCMP per-component sources (capability gap) or
MF-gamma phase-field gas inflow (VR-STR-02 sparger unit test is the
planned home). Routing package = branch cx/vv-sparger (example + report +
PNGs). The STOP-rule/honesty machinery worked as designed.

### ANOM-P4-014 CLOSED as characterization (cx/mp-hard rev 5)
Wetted-wall bookkeeping settled (wall_rho applies to ALL solids incl. rim —
compat/multiphase.rs:356, in-code comment added); with the correct
two-wetted-channel contrast the intercept collapses −20 → −1.4 and the law
is exact (r² = 0.99998). Remaining coefficient anomaly FROZEN: measured
slope 1312.7 = 1.54× the flat-wall theory 852.2 (in-situ θ ≈ 66° both
channels, so it is not an angle error). THREE-WAY REFEREE now stands on
the SC interface tension: Laplace σ (1.0 by construction, T11), Taylor-
Culick mechanical σ (0.49×, P4-017), Jurin meniscus σcosθ (1.54×, this) —
the lane-1.7 SC pressure-tensor audit is promoted to the highest-value W2
item; its job is to derive which σ the SC pressure tensor actually
delivers on curved menisci vs flat interfaces vs retracting rims.

### ANOM-P4-021 — body force × Zou-He face patch: secular mass leak — FIXED
in `cx/fix-p4-021`
Original measurement: steady-state discriminator confirmed persistent leak
after hydrostatic equilibration: uniform-force × patch +2.47e-5 mass/step,
gravity × patch −7.42e-5/step (rel 2.2e-9 / 6.7e-9 vs band 1e-9), scale
~F·A_patch. Root cause: the Zou-He patch reconstruction of unknown
populations ignored the Guo half-force contribution, so it imposed raw
momentum `rho*u_prescribed`; the subsequent `moments_row` half-force shift
reported physical velocity `u_prescribed + F/(2 rho)` and created secular
mass drift.

Fix: D2Q9, D3Q19, D3Q27, and generated WGSL `bc` closures now reconstruct on
raw Guo momentum `rho*u_phys - F/2` (implemented as the equivalent
pre-force velocity `v = u_phys - F/(2 rho)`, with analytic handling for
`F = F0 + rho*g`). Whole-face Zou-He and T18.2 face patches share this
closure. Zero-force branch remains bit-identical to the legacy D2Q9 formula.

Evidence: `zou_he_force.rs` passes for D2Q9/D3Q19/D3Q27 with uniform force
plus gravity; `kernels::tests::zou_he_d2q9_zero_force_matches_legacy_formula_bitwise`
passes; `feature_interaction_conservation_matrix` passes, flipping the two
documented red cells green while preserving the other cells' PASS/SKIP
status. GPU kernel text changed only in the `bc` reconstruction; F32
collision byte-identity scope was not widened or narrowed. Native GPU suites
remain PM-run.

### Lane 2.1 mutation-testing extension CLOSED (merge 8c26f14)
10-mutant matrix all CAUGHT: Guo-source sign, moving-wall factor 2, D2Q9
weight non-opposite swap, half-way source-cell off-by-one, probe-corner
link drop, TRT ω⁻ tau/λ typo, feq u² coefficient 4.5→4.0 — plus the
baseline 3 (MW sign, Zou-He pressure sign, outflow stale slot). Zero
survivors ⇒ zero silent physics-kernel mutants in the current validation
net. Complements ANOM-P4-021 (interaction-matrix): the pair-only defect
class is NOT reachable by single-mutation coverage — the two lanes are
complementary and both are needed.

### ANOM-P4-022 — Shan-Chen force-field OVERWRITE breaks additive composition — RESOLVED 2026-07-07 (S2, found by code-to-spec back-translation, lane 3.2)
`ShanChen::update_force` and `MultiComponent::update_forces`
`copy_from_slice` into `force_field` (compat/multiphase.rs:387), silently
discarding any prior rotor/gravity/user contribution — contradicts the
W-GRAV additive composition-point invariant. Rotor + SC coexistence needs
the caller to call SC FIRST then add rotor; the reverse order silently
zeros the rotor force. This is invisible in the current lane-5.1 matrix
(SC × rotor SKIPs by API-incompatibility). Fix landed in this worktree:
both SCMP and MCMP now add into the caller-owned per-cell force field, with
the same "caller clears once per step" contract as rotor. Regression:
`validation_multiphase.rs::t11_shan_chen_adds_to_existing_force_field_anom_p4_022`
asserts cell-by-cell `gravity + SC` composition and rejects either
contribution alone. Call-site audit: CLI scenario runner, WASM stepping,
interaction matrix helper, T11/T11b/T12/T13/pressure-tensor/hard-multiphase
tests now reset/zero-fill before composing transient SC/MCMP sources; MCMP's own
per-component gravity was already built in the same local force array and was
not a prior field overwritten by SC.

### ANOM-P4-021 DERIVATION CONFIRMED (from lane 3.2 code-to-spec) — CLOSED
Zou-He closures at kernels.rs:754-773 (D2Q9), 941-954 (D3Q19), 828-867
(D3Q27) reconstructed unknowns from raw populations and used caller-passed u
verbatim. Applied macroscopic velocity was therefore
u_prescribed + F/(2ρ), and mass leaked as F·A_patch per step — matching the
interaction-matrix measurement. Closed by the raw-Guo-momentum closure above.

### Documentation drift found (lane 3.2 residuals, S3)
- PHYSICS.md T11 σ entry does not acknowledge P4-014/017 mechanical-σ
  three-way referee — one-paragraph update recommended.
- Bouzidi "2nd-order" claim: three of four branches (qd≥1/2 and halo-edge
  cases without a second fluid node) silently degrade to 1st order —
  PHYSICS.md missing the validity clause.
- CLAUDE.md pass-order invariant + ARCHITECTURE_V2.md §3.4 elide
  apply_bouzidi / apply_volume_sources / swap; the real sequence is
  collide → halo → stream(interior/shell) → bouzidi → swap → open-BCs →
  volume-source → moments → end_step (backend.rs:258-320).

### V&V methodology audit findings (lane 4.3, docs/qa/vv-methodology-audit.md)
10 conforms / 10 gaps against ASME V&V 20 + Oberkampf-Roy. Top gaps (all
S3 process, none S0 physics): (1) no headline result carries a stated
U_num (Roache GCI never computed even though T8 is one step away —
common/gci.rs helper is the cheapest fix); (2) no vocabulary crosswalk
(U_num/U_input/U_D/U_val/δ_model vs our SPEC-GAP/STOP-RULE/frozen-band);
(3) T11/T11b/T11c frozen bands are regression pins masquerading as
validation gates (label VALIDATION_GATE vs REGRESSION_PIN + add one
external Maxwell-rule holdout for T11). No hidden calibration hack
survived the +0.0025 cleanup — process risk only.

## Pass 7 — 2026-07-07, SC pressure-tensor + band-retighten

### ANOM-P4-023 — SC MECHANICAL σ is 6x below the Laplace σ (lane 1.7 flagship result)
Direct pressure-tensor Kirkwood-Buff and Young-Laplace integrations
FINALLY QUANTIFY the three-way referee (P4-014 Jurin, P4-017 TC, P4-018/019):
- **P1 Young-Laplace σ**: 2.87e-2 at r=12 (rel error 13.5% vs T11), r=16
  9.1%, r=20 6.1% — **converges to sigma_Laplace as R→∞** — the discrete
  pressure-tensor IS the Laplace-consistent one. (Fails the ±10% band at
  the smallest R by 3.5 pp — freeze as R-dependence characterization.)
- **P2 SC pressure-tensor anisotropy**: 1.07% (well below the 15% band).
- **P3 KB flat-interface**: sigma_KB = 3.64e-2, rel 9.5% vs Laplace (in band).
- **P4 momentum-flux σ from Taylor-Culick rim**: **6.09e-3, rel 83% vs KB**.
The SC pressure-tensor integrates to a Laplace-consistent σ on flat and
curved interfaces (P1/P3 consistent within band), but the MOMENTUM
delivered by a moving interface is 6x smaller (P4). Verdict candidate:
**it is not σ that varies — it is the momentum coupling of the moving
interface**. The Taylor-Culick rim deficit (P4-017 v/v_TC=0.49) is
consistent with this — sqrt(0.49) is not 1/6, so the mapping is not
trivial, but the direction matches: the moving interface transmits
sub-nominal momentum. Related to P4-018 (near-wall shear artifact) and
P4-019 (contact-line immobility) — same SC "static-vs-dynamic" limitation
family. Rename Jurin 1.54× (P4-014) into this framework: static σ from a
STATIC meniscus is 1.54× flat-wall theory, but the KB flat σ is
Laplace-consistent, so P4-014's 1.54× may be the CURVATURE-dependent
static σ, not a dynamic issue. Hypothesis for lane 1.7 rev 2: P4-023's
r-dependent P1 (5.1e-3 shift over R∈{12,20}) IS the R-dependence that
projects onto Jurin's meniscus geometry. Investigation: measure σ_YL(R)
across a wider R range and cross-plot Jurin's h·w vs the theory using the
R-dependent σ.
- Documented-red gate: cx/sc-pressure-tensor.

### Lane 2.2 band-retighten COMPLETED (cx/band-retighten 8374f9d, merged 62cdc4a)
All 9 rows retightened using measured×20 (or physical-model caps with
stated headroom). Full workspace gate green (EXIT:0). Prints every
measured value in every retightened assert. Two-layer rule retrofit
applied to D3Q19/D3Q27 lid smoke gates (added upper bounds + magnitude
ratios). Zero tightened bands failed — the gates were purely vacuous.

### Merge-queue completion (task #12)
All 9 cx/vv-* branches merged into main (26137ae, 8d05b5a, 5b18a53,
fcc1eea, 08cdaf1, 47b4e83, ecb247c, 2050e50, 62cdc4a). CI cron / traceability
matrix / mutation runner / FSI safe-downgrade / multiphase anchors / GPU
absolute physics D2Q9 / scenario hardening / visual plots / evidence
templates all landed.

### Axis 9.4 sedimentation experiment observation (main 23a876d)
Harness LANDED and RUNS but its physics anchor was too loose: with the
posted parameters (d=1.5, rho_p=2, nu=1/6, g=5e-5) v_stokes=3.75e-5 gives
settling length 3.28e4 cells >> 128-cell basin, so 499/500 particles stay
suspended after 10k steps and the "factor-of-3" anchor is trivially met.
The visual (out/vv_sedim_2d/deposition_map.png) shows particles piled at
the inlet column only — a valid honesty artifact, not a physics defect.
Rev-2 needed: d=6+ (v_stokes scales as d²), longer run, or gravity ×10
(state Reynolds-limit tradeoffs). Mass conservation exact (n_deposited +
n_suspended = 500). Not filing an ANOM — this is an experiment-design
issue on my side to be fixed in the next revision.

### T9b convective outflow rev 1: rest-frame setup mismatch (self-triaged)
`main adcadcb`: convective-outflow reflection-coefficient sweep in a REST
CHANNEL measured baseline Outflow R = 0.26 vs ConvectiveOutflow R ≈ 0.998
at every u_conv ∈ {0.05..1.0}. Root cause (physics-honest self-triage,
not a defect): ConvectiveOutflow advects populations at u_conv; the
incoming disturbance is an ACOUSTIC wave traveling at cs = 1/√3 ≈ 0.577,
so any u_conv well below cs under-advects and reflects near-hard-wall.
Convective is designed for MEAN-FLOW use (u_conv ≈ local flow speed).
Rev 2 dispatches with the correct setup (uniform-inlet channel + u_conv
centered on u_in). Rev-1 test lands #[ignore]'d with the rationale
recorded in-code.

### Axis 9.4 rev 2: STOP-RULE self-triage (main 69d4ef2 addendum)
Rev-2 raised d=6 and g=1e-4: v_stokes=1.2e-3, deposition_fraction=1.0 —
BUT mean_deposition_x = 500 (out-of-domain) because particles advect past
the pressure outlet with the crossflow and keep sampling the clamped
boundary state until they hit the floor. Verdict per behavior-validity
review: ARTIFACT of an open-outlet + particle-domain-clamp interaction,
not a physics defect. Rev-3 dispatches with PM option (c): CLOSED BASIN
(all-BB, no crossflow, top-line seed, quiescent settling — the canonical
Stokes settling geometry). Cross-tests: physics is unchanged; this is
the third experiment-design pass on a single Axis-9.4 lane and each
revision is producing a stronger physics anchor. Filed in PHYSICS.md
(rev 2 dirty tree) — will fold into the rev-3 commit.

## Pass 8 — 2026-07-07, core fix landings (P4-001 / P4-021 / P4-022 CLOSED)

### ANOM-P4-001 CLOSED (merge 8a0c546)
IBM force sizing targets realized full-step F/rho with row-sum overlap
mobility + simultaneous Richardson sweeps (the collective-gain treatment
verified by our diagnosis). Default relaxation=1.0 stable at ds/h ∈
[0.39, 1.0]. Gate cx/audit-ibm: 8/8 GREEN in our environment.
Independent re-verify: 7 passed, 3 ignored, 1.65s.

### ANOM-P4-021 CLOSED (merge 434091f)
Zou-He closes on raw Guo momentum ρ·u_phys − F/2 (uniform force + gravity
analytic); whole-face + T18.2 patches share one corrected reconstruction;
GPU BC mirrored. Independent re-verify: interaction_matrix 1 passed
(previously 2 documented-red cells).

### ANOM-P4-022 CLOSED (merge 084dee6)
SC/MCMP force ADDS into the field with caller zero-fill; gravity+SC
composition regression added. Additive composition convention now
documented.

### ANOM-P4-016 STOP-RULE'D BY CORE (documented-red pin)
Core verified: i3 spec applies gravity to heavy component ONLY in a closed
box, creating net bulk downward momentum that drives wall-adjacent
failure before the RT cutoff mode is measurable. Additive-force fix did
NOT cure it. Options to us: (a) revise i3 forcing spec, (b) pressure-
balanced/zero-net-force buoyancy protocol, (c) literature MCMP RT
closure. Route: keep our documented-red pin, revise spec in next
mp-hard rev (option b — Boussinesq-analog with per-component gravity
adjusted so net = 0 at t=0 — most physical fix).

### ANOM-P4-010 PM RULING ACCEPTED (V&V concurs)
"Volume penalization's validity domain = thin/porous structures.
Coherent solid interiors → route to rotating IBM (now validated) or
Bouzidi." Physics-principled: penalization approximates a distributed
Darcy drag; coherent solid regions are outside its derivation. F1-F3
re-scope order dispatched (thin-blade IBM cross-referee stays as F4;
new F6 = paired forbidden-disc/valid-IBM domain-boundary witness).

**Open ledger:** P4-023 (SC σ referee, characterization), P4-016 (i3, our
spec revision), P4-018/019 (SC dynamic wetting family, characterization).
Core routings CLEARED: **P4-001, P4-021, P4-022 all closed**; V&V loop
has now driven 4 total core fixes to landing (P4-008/001/021/022).

### Axis 9.8 Taylor-Couette wavy-vortex — heavy STOP-RULE (rev 1)
Light laminar-Couette PASSES (bulk profile L2_rel 8.5e-2 within 10% band).
Heavy wavy-vortex onset: z-invariant seed + z-invariant forcing kept
axial-mode energy at 1e-14..1e-18 for all Ta (0.5, 1.5, 3.0 × Ta_c). This
is a well-known DNS-practice matter: axial-mode instabilities need a
finite axial seed to trigger the linear regime. Rev 2 dispatches with an
explicit small cos(k_z z) seed on u_r (eps=1e-4, k_z=2π/nz or 4π/nz,
chosen to match λ_c=2(R_o-R_i)). If the seed still fails to grow at
1.5*Ta_c, that IS a finding — rotating-boundary path filters the
instability, and the light Couette pass doesn't cover this regime.

### ANOM-P4-024 — thin-shell penalization filters axial modes (Taylor-Couette rev 2 real finding)
Explicit axial seed (eps=1e-4, cos(k_z z) at k_z = 4π/nz matching λ_c=32
= 2(R_o-R_i)) DAMPS through 5000 spectrum-window steps at every Ta:
- Ta=0.5*Ta_c: 1.24e-5 → 2.10e-9 (5000x damping)
- Ta=1.5*Ta_c: 4.15e-6 → 2.04e-8 (200x damping — should GROW here)
- Ta=3.0*Ta_c: 2.08e-6 → 2.32e-7 (10x damping — should grow more)
The rotating-boundary path (thin cylindrical shells via volume
penalization) FILTERS axial instabilities at every Ta up to 3*Ta_c. S2
characterization candidate (not a defect — physics-motivated: the
penalization sizing damps small-amplitude perturbations near the moving
boundary faster than Rayleigh instability can grow). Fits the P4-010
disposition family (penalization ≈ Darcy drag): the "thin structure" is
still coherent-solid at the axial-perturbation scale.
Route: rev 3 = replace thin cylindrical shells with rotating IBM
cylinders (now validated post-P4-001) as the physically appropriate
rotating boundary for the wavy-vortex study; if IBM path grows the
axial mode at 1.5*Ta_c, the finding confirms the penalization filter
class limit + provides the corrected acceptance route. Documented-red
until then; gate = cx/vv-tayc heavy.

### ANOM-P4-025 (formerly L1_7-001) — Bouzidi moving-wall qd<0.5 imposes σ·U_wall instead of U_wall
Native Bouzidi supports moving walls (wall_u[wall_ref] path); qd=0.5
degenerates to half-way MW bitwise; qd≥0.5 Couette profile matches
analytic within 2e-3. BUT qd<0.5 branch scales the wall speed by
sigma = 2*qd ≈ 0.5 at qd=0.25 (measured ~0.5·U_wall vs expected 1.0·U_wall).
This is the missing (sigma_i · 2·w_q · rho · c·u_wall / cs²) correction on
the qd<0.5 second-point interpolation branch — see Bouzidi 2001 §4;
symmetric qd>0.5 gets it right, qd<0.5 was skipped. Current-wrong-value
pin: qd=0.25 imposes σ·U_wall; #[ignore]'d all-qd exact test flips green
when the fix lands. S3 (bounded impact - most Bouzidi records fall on
qd>=0.5 side for well-resolved obstacles). Route: core-engine, small
one-file fix in bouzidi.rs qd<0.5 branch. Gate = cx/vv-bmw
qd_sweep_moving_wall_couette_should_match_offgrid_linear_profile_all_qd.

### Wen-2014 GALILEAN-INVARIANT PROBE — COVERED (radar #16 closed, pitfall #10 upgraded)
Decisive one-step test on cx/vv-wenmxg (main): co-moving frame with u_0 =
0.05 everywhere and BOTH walls MovingWall at u_0 measured F_probe_x =
4.44e-16, vs Ladd conventional prediction ρ·u_0·nx = 1.6, vs Wen-invariant
prediction 0. Signal/floor = 3.6e15. The current LBMFlow momentum-exchange
implementation is Galilean-invariant (Wen 2014-class) — pitfall #10
upgrades from PARTIAL to COVERED. Rotor/rotating-boundary force diagnostics
in moving frames are therefore trustworthy from a Galilean-invariance
standpoint; this closes a longstanding open concern in the pitfall list.

### ANOM-P4-024 STATUS UPDATE — rev 3 rotating IBM light PASSES, heavy still no onset (Ta_c estimate error)
Rotating IBM cylinders WORK for the LAMINAR case: annular Couette profile
L2_rel 6.87e-3, Linf/U_i 1.15e-2, IBM slip 2.18e-4 — well within the 15%
band. This CONFIRMS rotating IBM (P4-001 fix) as a valid coherent-solid
route. However, heavy wavy-vortex test still fails to grow the seed at
Ta=1.5*Ta_c with damping 1.66e-5 -> 4.56e-8 (200x). Diagnosis on my side:
Ta_c ~ 3390 is the NARROW-GAP Rayleigh estimate (valid for R_i/R_o -> 1);
our geometry is R_i/R_o = 12/28 = 0.43 (WIDE gap). Chandrasekhar 1961
Table X shows Ta_c climbs sharply for wide gaps — the true Ta_c(0.43) is
much larger than 3390 (~4-5x higher), so our "1.5*Ta_c" was still BELOW
critical. This is a TEST-DESIGN error, not a rotating-boundary filter.
Rev 4 queued: use R_i=20, R_o=24 (narrow gap 0.83) with proper Ta_c ~
1750 and the SAME 1.5*Ta_c and 3.0*Ta_c multipliers. ANOM-P4-024
stays OPEN as characterization; the underlying rotating-boundary path is
NOT filtering the mode — the test just used a wrong critical value.
Also: rev 3 CONFIRMED rotating IBM cylinders work in Taylor-Couette
setups (light PASS), which strengthens the P4-010 disposition endorsement
in a NEW use case.

### ANOM-P4-023 DEEPENED (sigma(R) sweep — radar #1, main 79d0cb0)
9-point R-sweep + 2-point light: sigma_YL(R) fits EXACTLY the inverse-R
Tolman form sigma_YL = sigma_inf + C/R with r² = 0.99977 (heavy 9-point:
r² = 0.99971). Fitted constants (heavy):
  sigma_inf = 3.610e-2, C = -1.085e-1
- At R→∞: sigma_YL → 3.61e-2, ~9% ABOVE sigma_Laplace(T11) = 3.32e-2.
- At R=12: sigma_YL = 2.71e-2, matches P4-023 P1 measurement (2.87e-2)
  within 6% (consistency check).
- At Jurin's meniscus r_m ~ 12 (gap 24): sigma_YL(12) = 2.71e-2
- Jurin-inferred sigma_eff (P4-014 slope 1.54× Laplace) = 5.11e-2
- **rel error 47%**: the R-dependence of sigma_YL DOES NOT explain
  Jurin's 1.54× enhancement.

Three-way referee status (Laplace T11: 1.00×, Jurin P4-014: 1.54×,
Taylor-Culick P4-017: 0.49×, sigma_KB flat P4-023 P3: 1.10×,
sigma_YL(R→∞) new: 1.09×) — flat-interface pressure-tensor sigma is
~10% above the T11 Laplace calibration, static curved-drop sigma matches
that in the R→∞ limit, both flat and asymptotic-curved cases converge to
~3.6e-2. But the WETTING slot (Jurin, sigma_eff 5.1e-2) is 40% above
even the asymptotic sigma_∞. Remaining hypothesis: WALL AFFINITY at the
solid rim adds a solid-fluid interface tension γ_sl that Young's law
naturally couples into cosθ_slot, effectively enhancing sigma·cosθ
without changing the bulk sigma. This is a CLOSURE hypothesis, not a
free-parameter fit — Young's γ_sl - γ_sv = γ_lv cosθ. Next investigation
(radar rev): measure γ_sl/γ_lv via the contact-line curvature in T11c
setup and see if 1.54× emerges from the correct Young-Laplace formulation
in a bounded slot. LEFT AS OPEN for a targeted follow-up; the σ(R) test
lands as a definitive characterization of the SC bulk sigma R-dependence.
Model-domain status: SC static sigma DOES obey a Tolman-length-type
correction, well-fit r²>0.9997 across R∈[6,32]. This is a physics finding
worth PHYSICS.md documentation.

### ANOM-P4-023 SC σ REFEREE — CLOSED WITH RESIDUAL (γ_sl direct measurement)
γ_sl KB-integral measurement (main dfc7146 → f018f8c):
- Wet wall (wall_rho=1.0, θ_T11c=63°): γ_sl = −2.517e-2, γ_sl/γ_lv = −0.70
- Neutral wall (wall_rho=1.888): γ_sl = −1.84e-3 (~0, as expected)
- **Wet−neutral shift**: 2.33e-2 (absolute); direct measurement of the
  SC wall-interaction contribution to the interfacial band.
- Young's law prediction γ_lv·cos(63°) = 1.63e-2
- **Wet−neutral shift is 43% ABOVE Young's prediction** — the SC solid
  interaction adds an EXTRA 7.0e-3 of interfacial-tension-band energy
  beyond what Young's γ_sv−γ_sl = γ_lv·cos(θ) alone would predict.
- Jurin-inferred σ_eff (P4-014) = 5.11e-2 = 1.42 × γ_lv (KB flat)

**Physical picture (final)**: SC static physics at solid-liquid contacts
is NOT purely Young's-law-consistent. There is a documented extra
contact-band tension of order 40% γ_lv that becomes visible in wetting
problems. This IS the physical origin of Jurin's 1.54× enhancement: the
SC contact line does more work per unit length than γ_lv·cos(θ) alone.
Verdict: (B) documented closure — this behavior is a known SC-model
characteristic (wall interaction Ψ_wall term acts on top of the bulk
cohesion), not a defect; the T11c contact-angle measurements themselves
build in the wall interaction, so downstream tests using T11c θ are
self-consistent. What this measurement OPENS is the option to (a) accept
the SC contact-line as a validated closure with 1.42× enhancement factor,
(b) route future high-precision wetting work to MF-γ phase-field.

**Three-way referee CLOSED** with the residual explained:
- σ_bulk (KB flat, Laplace R→∞) = 3.6e-2 (bulk sigma, Tolman-consistent)
- σ_wetting_effective (Jurin, T11c contact line) = 5.1e-2 = 1.42× σ_bulk
  → EXPLAINED by the extra wall-interaction contact-band tension
- σ_dynamic (Taylor-Culick rim) = 1.6e-2 = 0.49× σ_bulk → still open,
  moving-interface momentum coupling limit (P4-018/019 family)

Radar #1 CLOSED, ANOM-P4-023 CLOSED as characterization. PHYSICS.md entry
draft prepared in worktree; will land with next mp-hard rev.
