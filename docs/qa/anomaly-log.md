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
