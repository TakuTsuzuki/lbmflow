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
