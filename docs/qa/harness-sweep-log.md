# Physics Anomaly Sweep — harness pass 1 summary

Automated physics-QA over runnable (scenario × config) pairs. Reference
values come ONLY from `docs/VALIDATION.md` frozen bands or analytic
solutions. Harness: `scripts/qa/` (matrix.py, qa_checks.py, run_sweep.py).

## Pass 1 — 2026-07-06, `qa/anomaly-sweep`

Matrix: 18 configs across 8 tracks (conservation / channel-analytic /
open-boundary / cavity / cylinder / 3D / multiphase / robustness), all via
the `lbm run` scenario path. f64: 2D D2Q9 CpuSimd (compat) + 3D D3Q19
CpuScalar. Evidence:
[pass1/results-main.json](pass1/results-main.json) (18 configs) and
[pass1/results-rerun-fixed.json](pass1/results-rerun-fixed.json)
(4-config rerun after harness fixes). Reproduce via
`python3 scripts/qa/run_sweep.py --bin target/release/lbm --out out/qa-pass1`.

**Headline: no engine physics anomaly (no S0/S1).** All 18 configs meet
their frozen-band or analytic checks. All four initial failures were
harness / matrix-fidelity defects (`xy_mirror_symmetry` axis, snapshot
cadence vs steady-stop, T4 ν parameter drift, staircase parity) — fixed
in-harness and re-verified in the rerun.

Findings landed to `docs/qa/anomaly-log.md`:
- ANOM-P1-001 (manifest lacks first-class conservation diagnostics — S2
  tooling proposal).
- ANOM-P1-002 (T4 profile band calibrated to frozen ν=0.02 — S3 spec
  footnote; matrix pins ν=0.02).
- ANOM-P1-003 (runtime max|u| can exceed compressibility advisory silently
  — S3 monitoring proposal).

Coverage gaps recorded in the pass 1 log (git history): T1 Taylor-Green
not expressible via scenario JSON (no analytic-init surface); buoyant
droplet gap (`physics.force` is constant density; no per-mass gravity
schema); no Uz FieldKind for 3D output; T9 horizon capped at 60k vs
spec 1e5; f32 configs are pass-2. These are matrix/schema surface
proposals, not engine findings.

For the raw per-config metrics table, harness errata, and
expected-limitation notes see commit history around this file.
