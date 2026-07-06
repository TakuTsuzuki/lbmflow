# Physics Anomaly Sweep — log

Append-only. Filed by the QA/viewer session (autonomous pass; `send_message` to
PM is pending the auto-approve hook, so S0/S1 are recorded here for PM pickup
instead of messaged). Sink + taxonomy per the agreed protocol
(S0 silently-wrong physics / S1 divergence-leak in a supported config /
S2 below-expected-accuracy / S3 minor). Anomaly references are ONLY the
documented `lbmflow-user-tune-stability` thresholds and analytic constraints —
no invented thresholds (band governance).

---

## Pass 1 — 2026-07-05 — stability-envelope sweep (3D stirred-tank example)

**Config**: `crates/lbm-cli/examples/stirred_tank_3d.rs`, D3Q19 TRT, n=64, 2000
steps, penalized Rushton impeller. Sweep of (tip speed `u_tip`, viscosity `nu`)
to probe the documented stability envelope. Detector: draft
`sim-anomaly-scan` (universal checks + the tune-stability thresholds
tau≥0.55, |u|≤0.15 Ma / hard 0.3, grid-Re U/ν≤15).

**Raw sweep** (example built-in verdict `||` scanner verdict):

| case    | u_tip | nu     | tau   | Ma_tip | grid-Re U/ν | example | scanner |
|---------|-------|--------|-------|--------|-------------|---------|---------|
| ctrl    | 0.08  | 0.02   | 0.560 | 0.14   | ~4          | STABLE  | CLEAN |
| ma_hi   | 0.20  | 0.02   | 0.560 | 0.35   | ~10         | STABLE  | CRITICAL (mach) |
| ma_xhi  | 0.30  | 0.02   | 0.560 | 0.52   | ~15         | STABLE  | CRITICAL (mach) |
| tau_lo  | 0.08  | 0.004  | 0.512 | 0.14   | ~22         | STABLE  | FLAGGED (tau, grid-Re) |
| tau_xlo | 0.08  | 0.0015 | 0.504 | 0.14   | ~65         | STABLE  | CRITICAL (tau, grid-Re) |
| gridre  | 0.12  | 0.002  | 0.506 | 0.21   | ~1160→NaN   | DIVERGED| CRITICAL (mach, tau, grid-Re) |

### A1 [S0] — out-of-envelope runs stay bounded and report STABLE with no runtime signal
`ma_hi` (Ma_tip 0.35) and `ma_xhi` (Ma 0.52) exceed the low-Mach limit
(|u|≤0.15 Ma, hard 0.3) yet stay bounded (the penalization cap holds |u|<0.3)
and the run reports **STABLE**. Above Ma 0.3 the method no longer approximates
the incompressible NSE — the field is silently compressible/wrong. Same class:
`tau_xlo` (tau 0.504 at the ν→0 floor, grid-Re 65) runs bounded for 2000 steps
and passes naive finiteness/divergence, but is grossly under-resolved →
physically meaningless.
**Not** an "expected-limitation": the limits ARE documented, but the run path
gives **no runtime signal** — no max-Ma / grid-Re echo, no warn — so a user who
skips pre-run `lbm validate` gets bounded, plausible-looking, wrong output.
**Ask (core)**: echo `max_Ma` and `grid_Re` into the run manifest and emit a
runtime warn (or opt-in abort) when they cross the documented thresholds, so the
silently-wrong regime is not silent. Ref: tune-stability |u|≤0.15 Ma, U/ν≤15,
tau≥0.55.

### A2 [tooling, mine] — example STABLE/DIVERGED gate too permissive → FIXED
The gate was `final_max|u| < 0.5` (caught only full divergence, missed A1). Now a
three-state verdict: `DIVERGED` (non-finite or |u|≥0.5) / `OUT-OF-ENVELOPE`
(bounded but Ma_field>0.3 or grid-Re U/ν>15) / `STABLE`; the SUMMARY line also
prints `Ma_field` and `grid_Re`. Verified: ma_hi and tau_xlo now report
OUT-OF-ENVELOPE (were STABLE). Commit-side change in the example.

### B [tooling fix, mine] — draft anomaly-scan missed the under-resolved case → FIXED
First run returned CLEAN for `tau_xlo` (no tau/grid-Re check). Added the
documented-threshold checks (tau≥0.55, U/ν≤15); it now flags `tau_xlo` CRITICAL.
Also added `nu` to the example's `volume.json` so the scan auto-applies them with
no `--nu` flag. Verified. **Hand-off to the core-owned `sim-anomaly-scan`**: adopt
the tau-floor + grid-Re checks from the frozen thresholds; draft script at
`scratchpad/draft_anomaly_scan_for_worker.py`.

### Divergence boundary observed
Stable at (0.08, 0.02); the first hard NaN divergence in this sweep is
(0.12, 0.002) — tau 0.506 with grid-Re ~1160. Consistent with the documented
envelope (tau near the floor AND grid-Re ≫ 15).

### Coverage NOT run this pass (deferred to the worker + CLI collection surface)
Single-phase 2D analytic cases (Poiseuille/Couette profile L2, Ghia cavity
centrelines), cylinder T8, 2D Shan-Chen spurious-current vs T11 bands. These
need the CLI preset/scenario + VTK/manifest surface the PM described; the worker
owns `sim-run`/`sim-anomaly-scan`/`sim-qa-report`. This pass used the
self-contained 3D example to calibrate the detector against the documented
stability thresholds.

**Pass 1 verdict**: detector calibrated and catching S0 correctly. Open for PM:
**A1 (S0)** — add a runtime Ma / grid-Re guard so out-of-envelope runs are not
silent (core). Closed by me: A2 (example 3-state gate) + B (scanner tau/grid-Re
+ nu-in-export), both verified. STOP for PM go per protocol (no unattended loop).

Next pass (needs the CLI collection surface + worker): 2D analytic cases
(Poiseuille/Couette L2, Ghia cavity centrelines), cylinder T8, 2D Shan-Chen
spurious-current vs T11 — the worker owns run/scan/report; it should request the
`lbmflow-qa-viewer` skill for any spatially-flagged case rather than rebuild.

---

## Pass 2 — 2026-07-06, MF-interim Wave 1 (gravity / rotor / particles) + resuspension observation

Context: user directive to take over the stirred-tank resuspension capability.
Wave-1 codex orders landed per-mass gravity (`cx/mf-grav`), rotor volume
penalization (`cx/mf-rotor`), one-way Lagrangian particles (`cx/mf-particles`)
and an adversarial suite (`cx/mf-tests`); integrated + scenario/runner wiring
on `qa/mf-integration`. Adversarial suite: **16 pass / 0 fail / 4 SPEC-GAP**
(`cargo test -p lbm-core --release --features mf-interim --test mf_interim`).

### Anomalies

**ANOM-P2-001 — uniform-force vs per-cell force-field transient impulse
mismatch** — S2 (correctness of transients; steady-state invisible),
disposition-proposal: collision-kernel owners (B-1/R2-C in flight) unify the
source-term weighting.
- Scenario+config: any TRT run driving the same F through `SimConfig::force`
  vs the per-cell force field (gravity / Shan-Chen / rotor path).
- Expected: identical dynamics (Guo forcing, single definition — REQ rev.4
  "forcing second-moment single-definition" is the same family of issue).
- Observed (probe, 32x24 periodic + obstacle, tau=1, TRT Lambda=3/16,
  F=3e-7): uniform path u(1) = 1.5 F (exact Guo); force-field path
  u(1) = 0.9286 F — a one-time impulse deficit of 1/(2 tau_minus) * F = 4/7 F.
  Growth is F/step on both paths afterwards, so T2/T6/T11 steady gates cannot
  see it; the offset then seeds slowly diverging trajectories near obstacles
  (measured 2.1e-7 growth divergence over 50 steps).
- Impact: any transient force-driven measurement (SC droplet oscillation
  phase, rotor spin-up torque transient, gravity startup).
- Workaround in tests: same-path twins only (see mf_interim.rs).

**ANOM-P2-002 — rotor blade indicator produced mirror arms for odd blade
counts** — S2 (wrong geometry, silently plausible fields), **fixed** in
`qa/mf-integration` (along-blade sign check; 3 blades were 6 half-thickness
arms). Found by cross-reading the adversarial suite's independent geometry
against the implementation. The frozen stability envelope used 4 blades
(even) and is unaffected.

### SPEC-GAPs raised by the adversarial suite (S3, to pin in the contracts)

1. Native `Solver::set_gravity` composition with `set_body_force_field`
   rewrites across subdomains. 2. Particle starting inside solid: project /
   reflect / reject. 3. Rotor `chi = 0`: rejected vs no-op. 4. Particle
   overlap: no collision model (document one-way explicitly at the scenario
   surface).

### Behavior observation — stirred-tank flake resuspension (the target case)

2D reduction of the user case (128^2 closed tank, per-mass gravity
2e-4, 300 one-way flakes, 4-blade penalized rotor; configs in
`scripts/qa/resuspension/`, stats via `scripts/qa/observe_resuspension.py`):

- **High-clearance impeller (C/T = 0.5), tip 0.10**: flakes settle to the
  floor and are NEVER resuspended (laminar near-floor vertical velocity
  1.9e-3 total speed at row 1, mostly horizontal, vs settling velocity
  1.8e-3). Permanent deposition — consistent with mixing practice (high
  clearance is bad for solids suspension).
- **Low-clearance impeller (C/T = 0.25, the textbook solids-suspension
  geometry), d=6, rho_p=1.03**:
  - tip 0.01: 100% settled at 100k steps (one transient single-particle
    pickup to y=60 that re-settled — threshold intermittency).
  - tip 0.10: **sustained partial suspension — 34% of flakes above y=16,
    excursions to y=117/128, 66% remaining as a bed** (the classic below-N_js
    bed + cloud equilibrium). Flow field: `out/.../speed_100000.png`.
- Verdict: settling, threshold behavior and rotation-driven resuspension all
  emerge from the implemented force balance (no stochastic kicks). Honest
  limits: laminar only — realistic resuspension at industrial Re needs W-LES
  (MF-beta, this session's next order) + lift forces (FR-PART roadmap); free
  surface (half-filled vessel) still rigid-lid until MF-gamma.

### Process errata (pass 2)

- Stale-binary trap: `cargo test -p lbm-core` does not rebuild lbm-cli; the
  first observation ran a pre-rotor binary (maxSpeed 0.0000 exposed it
  immediately — the observation gate works).
- codex sandbox cannot commit in shared-.git worktrees (`index.lock` EPERM):
  2 of 4 orders needed PM-side commits. Fold into the dispatch Skill notes.

---

## Pass 3 — 2026-07-06, W-LES heavy characterization freeze

Ran on `origin/main` 99bb32a (post B-1 + cx/acc + cx/wles + W-GRAV + W-ROT
landings). No new anomalies; two frozen turbulence-tractability data points.

### Frozen values (both #[ignore], measured on this pass)

- **TGV64 nu_eff shift under WALE** (T15.4 setup, N=64, nu=0.02, u0=1.28e-4/N,
  tstar=832 steps, ~35 s wall):
  - nu_eff_off = 1.9977e-2, nu_eff_on = 1.9977e-2 → **dnu_rel = 6.60e-8**
  - max nu_t (on) = 1.39e-8 (essentially null under diffusive scaling — the
    intended WALE behavior for a small-strain resolved flow)
  - Band frozen at 1e-6 (~15x headroom over measured value); the original
    order allowed 1% which is far too loose for what WALE actually did here.
- **Multimode stabilization existence proof** (deterministic 3-mode init,
  N=48, nu=0.003, u0=0.10, U/nu=33, 20k steps, ~7 min wall total):
  - LES-OFF diverged at step 200 (max|u| > 0.3 or non-finite)
  - LES-ON completed 20000 steps: max|u| = 5.15e-4, max nu_t = 5.08e-6
    (~0.17% of nu_0 — a real, non-trivial modeling contribution)
  - Horizon extension: **100x** (200 → 20000 steps).

### Honest scoping (what these do NOT prove)

- **Turbulence ACCURACY is still open** — the Re_tau=180 channel vs DNS test
  (T17/VR-STR-03) remains a skeleton. The paper's turbulence-accuracy claim
  stays red on the claims ledger until that lands.
- The multimode case proves *stabilization exists*, not that the stabilized
  solution is quantitatively correct at that Re — it is a tractability seed,
  not a validation.

### Physics-honesty check on the WALE choice (bent-physics avoidance)

The steady-Couette/Poiseuille null gate (measured max nu_t <= 1e-12) already
proved WALE cannot leak into resolved pure-shear physics. TGV64 (small S^d
regime, dnu_rel = 6.6e-8) confirms it does not touch the diffusive-limit
scaling either. Multimode (large S^d) shows it activates as expected. Three
regimes, three consistent behaviors — the WALE-over-Smagorinsky ruling holds
by measurement, not just by cited theory.
## Dry run — 2026-07-06, lbmflow-accuracy-audit Skill on Bouzidi curved BC

**Skill**: `.claude/skills/lbmflow-accuracy-audit` (v1). **Test file**:
`crates/lbm-core/tests/accuracy_audit_bouzidi.rs`. **Branch**:
`qa/skill-accuracy-audit`. **Result**: 2 default light tests pass (G1 A1
convergence order slope 1.993 r² 1.0000; G2 A2 sub-cell translation spread
5.69% peak-relative), 3 ignored (G4/G6 SPEC-GAPs carrying derivations, G5
heavy tau sweep). **Zero engine bugs surfaced.** One P3 test-side finding:

**ANOM-DRY-001 — convergence-order fit x-axis reversal** — S3 (minor, test
tooling; taxonomy: fit-parameterization). **Disposition: test-fix (self,
in-worktree).**
- Scenario+config: G1 in accuracy_audit_bouzidi.rs, off-grid Poiseuille,
  ny ∈ {22, 30, 42}, diffusive scaling F ∝ 1/width², nu=0.04.
- Expected: `order_fit` slope ≈ +2 (Bouzidi second-order at fractional wall).
- Observed (first pass): slope = −1.993, r² = 1.0000. Raw errors
  `[2.50e-3, 1.29e-3, 6.41e-4]` at width `[20.4, 28.4, 40.4]`.
- Root cause: passed `width` (grows with resolution) to `order_fit` instead
  of `1/width` (mesh spacing → 0 as we refine); fit landed at
  `err ∝ width^{-p}` with p ≈ 2.
- Impact (bounded, this file only): none in engine. Trap documented in-line
  in the test's `hs.push(...)` comment so the next audit avoids it.
- **Calibration**: this is one more test-side derivation error on top of the
  6/6 from the accuracy_audit.rs pass; the Skill's P3 discipline (derive
  before blaming) caught it in one pass and saved a wasted engine-fix order.

---

## Pass 4 — 2026-07-06, FSI coupling-surface accuracy audit (V&V session)

Context: goal "FSI simulator V&V from every angle". Audit list:
session scratchpad `audit-list-fsi.md`; adversarial suites authored by codex
orders on `cx/audit-probe` (momentum-exchange probe + walls, 5 probes) and
`cx/audit-ibm` (rotating IBM, 8 probes). Parallel ad-hoc-physics inventory
delivered in `docs/proposals/adhoc-inventory-2026-07-06.md` (ORDER Task 1).

### Verified-exact results (probe suite 5/5 green after triage)

Momentum-exchange force measurement on the compat 2D path is VERIFIED at
round-off precision: (a) per-step global momentum ledger
`dp = N_fluid*F − F_probe(all solid)` closes to ~1e-14 abs over 10
consecutive steps with an interior obstacle + mixed shear field (probe
completeness AND pre/post-collision timing are both exact); (b) steady
Poiseuille wall-friction balance F_top+F_bottom = +g*N_fluid to 2.7e-15
with exact per-wall split; (c) static-pressure normal push = nx*rho*cs^2
per wall to 12 digits (O(u^2) equilibrium corrections cancel exactly —
derivation in test); (d) x-mirror equivariance of probed force at 1e-16;
(e) Couette probed shear = rho*nu*U/H*nx to ~1e-10 at tau∈{0.6,1.0};
(f) moving-wall momentum term uses LOCAL density (rho0=1.05 Couette exact,
linf_rel 1.4e-11 — a rho=1 hardcoding would read 0.952).

### ANOM-P4-001 — time-stepped direct-forcing IBM diverges in legal (incl.
DEFAULT) configurations — S1 (divergence leak), disposition: **core-engine
routing** (per the 2026-07-06 V&V routing rule; adversarial suite
cx/audit-ibm B1–B8 is the acceptance gate)

- Scenario+config: 80×80 periodic D2Q9, TRT Λ=3/16, nu=1/6; IBM circle
  r_i=10, Ω=1.5e-4 (Re_r=0.09, maximally benign), outer solid rim at r>30.5;
  per-step pattern `clear_body_force_field(); apply_rotating_ibm(body,cfg);
  step()` — the pattern used by rotating_ibm.rs itself.
- Expected: Uhlmann/Wang direct forcing is stable in time-stepped use at
  relaxation 1.0 (that is the literature's direct-forcing value, and the
  module DEFAULT: max_iterations=3, relaxation=1.0).
- Observed (scratch sweep, 2000 steps, torque vs analytic annular-Couette
  T=4πμΩr_i²r_o²/(r_o²−r_i²)=3.506e-2):
  - relax=1.0: NaN at n_markers∈{63,160} (ds/h∈{1.0,0.39}), iters 1 or 4
  - relax=0.5: NaN at n=160; stable at n=63 with T ratio 1.27
  - relax=0.25, n=160: stable, T ratio 1.075 (+7.5% — consistent with the
    ±11% diffuse-interface radius ambiguity at r_i=10), slip_max_rel 3.2e-5
  - empirical stability threshold ≈ relax ≲ 0.5·(ds/h); n=31 (ds=2.0) at
    relax=1.0 stays finite but T is garbage (ratio 76, boundary leakage)
  - before NaN the fields are silently wrong: at step 200 torque is 19×
    steady scale and Ω→−Ω antisymmetry is broken at O(1) (audit B6)
- Root-cause hypothesis (derived, two compounding terms): (i) the marker
  force is sized so the interpolated Guo HALF-force velocity increment
  Σ W·F/(2ρ) equals the sweep slip (comment at solver.rs:2110-2116), but the
  realized full-step Guo momentum change is F/ρ — a 2× overshoot per apply;
  (ii) overlapping kernels of neighboring markers (ds<h) amplify the
  collective gain further (Uhlmann's known spreading-interpolation
  eigenvalue growth). (i) alone predicts neutral oscillation at relax=1,
  (i)+(ii) predicts the observed ds-dependent divergence.
- Impact: any time-stepped IBM use (VR-STR-01 Np/torque, MF-δ rotating
  boundary) with default or literature-standard config; existing
  rotating_ibm.rs tests mask it by running relaxation=0.05 with 1 sweep and
  near-vacuous bands (L2_rel<0.95, Linf/U_i<5.8, torque never asserted).
- Positive control: at the stable point the steady physics is RIGHT
  (torque +7.5% of analytic, inside the diffuse-radius ambiguity) — the
  defect is the temporal coupling loop, not the spatial discretization.
- Repro: `cx/audit-ibm` branch, tests b1/b2/b5/b6/b8 (currently failing,
  NaN or asymmetry); scratch sweep preserved in the routing package.

### ANOM-P4-002 — A2 wall-balance sign expectation was wrong in the audit
row — S3, test-side (P1 derivation error by the PM), **fixed in-worktree**
- Expected (wrongly): ΣF_probe = −g·N_fluid. Correct: probed_force is the
  force ON the probed solid, so steady walls absorb +g·N_fluid (Newton pair
  of the body-force injection). Measured +gN/2 per wall exactly.
- Also mis-specified: "no normal wall force" — a resting wall in
  near-equilibrium fluid receives the static-pressure push nx·ρ·cs² (the
  O(u²) equilibrium corrections cancel in the c_y>0 link sum). Measured
  exact to 12 digits; the fixed test asserts the derived value.

### ANOM-P4-003 — A1/A5 tolerance floors under-modeled — S3, test-side,
**fixed in-worktree**
- A1: band scaled by max(N|F|,|F_probe|) under-floors the cancellation error
  of differencing O(N)-term sums of magnitude |p| (measured 1.6e-14 abs vs
  3.8e-15 band). Fixed: denominator |p_t|+N|F|+|F_probe|, factor 1e-11.
- A5: the API exposes one probe set at a time so the two walls are read on
  consecutive steps; near the 1e-11 steady criterion the residual drift
  bounds cancellation at ~3e-10 relative. Fixed: band 1e-8 with the
  measurement-protocol comment.

### Observation (no anomaly) — IBM force spreading is conservative even
under domain-edge kernel truncation: a body tangent to a wall rim measured
momentum_error_rel = 6.2e-15 (audit B3 "conservative surprise" branch) —
the mobility normalization renormalizes truncated kernels. Recorded as a
contract observation; the SPEC-GAP candidate is closed as not-a-gap.

### Pass 4 continued — particle (cx/audit-particles) and sources
(cx/audit-sources) suites

**Verified-exact results (both suites green after triage)**: trilinear
sampler exact on affine fields (6.9e-18); buoyancy sign antisymmetry exact
(rel 0); step vs step_depositing bit-identical over 200 steps; the drag
integrator is PINNED as semi-implicit (backward) Euler with SN drag —
agreement 3.3e-7 across tau_p ∈ {2,5,10,40}, unconditionally stable at
tau_p = 0.4 (monotone, no sign flips) — consistent with T18.3 passing at
tau_p = 1e-4. Near-wall sampler contract pinned: solid-node zeros are
blended linearly (u decays toward a stationary wall; SPEC-GAP noted for
moving-wall particle coupling — the Sample has no wall-velocity channel).
Jet per-step momentum ledger EXACT (abs err ~1e-24 vs band 6e-22, all axes,
8 steps) once q_lu semantics are read as region-total; source-dipole far
field r^-3 verified (slope 3.04, r² 0.9995); source and patch mirror
equivariance at 1e-16..1e-25; masked-patch BC nodes exact at machine
precision ON THE FACE LAYER (vel 3.5e-18, lid 4.1e-20, rho 1.1e-16,
patch-lid seam 4.1e-20).

**ANOM-P4-004** — C1/C2 first-pass failures (rel 1.3e-3) were TEST-SIDE:
hitting tau_p = 2 needs d = 1.34 lattice units, so v0 = 1e-4 gave
Re_p = 1.3e-3 and the SN factor polluted the Stokes-limit identity by
exactly the observed deviation. Fixed v0 = 1e-10 (SN residual 1.2e-7 <
band 1e-6). S3, test-fix (in-worktree), derivation in the test comment.

**ANOM-P4-005** — D1 first-pass failure was TEST-SIDE + a DOC GAP: the
audit order read `SourceKind` q_lu as per-cell; the implementation (and the
T18.1 mass-ledger wording Σ q_lu) means REGION TOTAL. First measurement
over-predicted by exactly N_region = 64. Fixed; the identity then holds to
round-off. S3 doc action: make the per-region-total semantics explicit in
DISPERSED_DEPOSITION.md §5 (queued with the D-track PM).

**ANOM-P4-006** — D4 first-pass failure was TEST-SIDE (sampling layer): the
patch BC nodes are the face layer z = nz−1 (exact at machine precision);
the adjacent interior layer carries developed flow (7.4e-3), not BC error.
Convention pinned in the test.

**ANOM-P4-007** — cumulant viscosity-offset audit (cx/audit-cumulant,
Order E) first pass: verdict **OPEN, order being revised**. What stands:
orientation consistency of the scalar correction PASSED (spread 2.2e-10);
at the calibration point (N=32, diffusive u0) residual −5.9e-4. What is
confounded (audit-design side, not yet engine evidence): E1's u0-sweep at
fixed N conflates O(Ma²) compressibility with the cubic defect (TRT control
c=7.98 vs cumulant c=8.22 — both dominated by compressibility); E2's N=24
band ignored the O(h²) spatial-error floor of the nu_eff fit (~1.7e-2 at
N=24). GENUINE flag to resolve: D3Q27 with offset=0 shows defect 9.1e-3 at
N=32 where D3Q19 shows −5.9e-4 — 15x, unexplained by spatial error at that
N; the "D3Q19-only bias" story needs the revised probes. Revision order
queued: tau-sweep to expose the omega-space correction's tau-dependent
nu-space footprint, N-sweep at fixed Ma for D3Q27, spatial-error-modeled
bands. 3.1 (B) vs (C) reclassification stays undecided until then.

**ANOM-P4-008 — cumulant D3Q19 "+0.0025 viscosity offset" is a
resolution-point calibration that corrupts the continuum limit — VERDICT
(C), S2 (leaning S0: silent systematic viscosity bias at every N except the
calibration point). Disposition: core-engine routing; acceptance gate =
cx/audit-cumulant e2 canary (must pass with |a| <= 2e-3 after removal).**
- Measured (Order E rev 2, heavy 3-point, N ∈ {24,32,48}, diffusive u0):
  D3Q19 defects d(N) = nu_eff/nu − 1 = [+1.7035e-2, −5.943e-4, −1.3200e-2]
  fit d = a + b/N² EXACTLY (all three points on the line): a = −2.3275e-2,
  b = +23.22. D3Q27 control (offset = 0): d(N) = [+1.6224e-2, +9.080e-3,
  +4.010e-3] → a27 = −6.9e-5 ≈ 0, b27 = 9.38 (pure spatial error).
- The smoking gun: the D3Q19 intercept equals the offset's own nu-space
  footprint −0.0025·2/(2−ω) = −2.333e-2 at ω = 1.7857 (nu = 0.02) — a
  99.8% match. Therefore the UNCORRECTED D3Q19 cumulant has ~zero
  resolution-independent bias; what the offset "fixed" at N=32 was the
  ordinary O(h²) spatial error (b/N² = +2.27e-2 there), i.e. a constant
  calibrated to cancel discretization error at one resolution. It
  under-corrects at N=24 (+1.7e-2), over-corrects at N=48 (−1.3e-2), and
  injects a tau-dependent −0.0025·2/(2−ω) viscosity error in the refined
  limit. This is exactly the prime directive's banned class.
- The −0.16·u² term: no clear nu_eff footprint at the predicted size
  (cum−TRT slope differences +0.057/+0.081 vs ~0.37/0.11 predicted with the
  <u²> = u0²/4 TGV weighting), but operator-intrinsic differences confound
  the residual — needs a core-side ablation toggle; E1 is SPEC-GAP'd with
  the N=48 dataset recorded. Its provenance inherits the same calibration
  concern; recommend re-deriving or removing together with the offset.
- Recommended core action: remove the +0.0025 offset; re-freeze TGV3D
  acceptance with a resolution-aware criterion (the h²-intercept |a| ≤ 2e-3
  in the audit file IS that criterion); decide the u² term by ablation.
- E4 note: the correction family is orientation-consistent (spread
  2.2e-10) — the defect is provenance/magnitude, not anisotropy.
