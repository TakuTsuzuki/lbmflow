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
