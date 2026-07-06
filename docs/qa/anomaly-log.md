# Physics Anomaly Sweep — log

Living register of physics anomalies. Rows drop when superseded by a landed
fix + regression pin (or refuted). Test-cited stubs stay: they document
current-wrong-value pins that must flip in the same commit as the fix.

Taxonomy: **S0** silently-wrong · **S1** divergence-leak · **S2**
below-expected accuracy · **S3** minor. Historical triage narrative is in
git history (2e121c8 band-vacuity scan, 714da6a cold-review triage, pass 1–4
context in commits 82946d3 / cd3999f / 20d0e10).

---

## OPEN — core-engine fix pending

### ANOM-P2-001 — uniform-force vs per-cell force-field transient impulse mismatch
- Severity: **S2** (steady-state invisible; transients wrong).
- Cited in: `crates/lbm-core/tests/mf_interim.rs:265`,
  `crates/lbm-core/tests/accuracy_audit.rs:471` (`#[ignore]` current-wrong-value
  pin — flips at R2-C fix).
- Measured (32×24, tau=1, TRT Λ=3/16, F=3e-7): uniform u(1)=1.5F (exact Guo);
  force-field u(1)=0.9286F — 1/(2 tau_minus)·F = 4/7 F impulse deficit.
- Disposition: fold into R2-C mechanical TRT port (order text in
  `scratchpad/order-r2c.txt`).

### ANOM-P4-001 — time-stepped direct-forcing IBM diverges in default config
- Severity: **S1**.
- Gate: `cx/audit-ibm` B1–B8 (currently RED / NaN).
- Config: 80×80 periodic D2Q9 TRT Λ=3/16 nu=1/6, IBM circle r_i=10 Ω=1.5e-4
  (Re_r=0.09), relaxation=1.0 (module DEFAULT). NaN at n_markers∈{63,160}.
- Root cause (derived): marker force targets the Guo half-force velocity
  (Σ W·F/(2ρ) = slip), but the full step realizes F/ρ — 2× overshoot per apply;
  overlapping kernels (ds<h) amplify collectively. Same family as P4-010.
- Positive control: at stable point (relax=0.25, n=160) T ratio 1.075 —
  spatial discretization is right; the coupling loop is broken.
- Disposition: core-engine routing.

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
