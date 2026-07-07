# Worked example — Bouzidi curved BC (dry run, 2026-07-06)

The full excavate→encode→triage loop applied to the Bouzidi interpolated
bounce-back subsystem (`crates/lbm-core/src/bouzidi.rs`), landed 2026-07-06,
as the Skill's own calibration and worked example. Test file:
`crates/lbm-core/tests/accuracy_audit_bouzidi.rs`.

## Scope & lane safety

- Subsystem: Bouzidi interpolated bounce-back — the qd < 1/2 and qd > 1/2
  interpolation stencils, the wall-velocity term, and the ray-intersection
  qd builders. Approximation to bound: how far Bouzidi departs from the
  "half-way BB is second-order at the wall midpoint" baseline as the wall
  moves off the midpoint.
- Lane safety (PM ruling, 2026-07-06): `r2-bouzidi` is freshly merged to
  main; sub-cell translation invariance is precisely the axis that would
  have caught the cylinder-centre-convention bug fixed in 52eaf85 a priori
  (+8% Cd bias). Constraint: **CPU-only** — GPU Bouzidi does not exist
  until ME-1 lands, and P2 template convention 5 keeps audits CPU-only in
  general so T14 gates every backend transitively.

## P1 — EXCAVATE (audit list)

| # | Approximation | Analytic ref & derivation sketch | Order/band | Axis | Cost |
|---|---|---|---|---|---|
| G1 | Interpolation-stencil convergence at fractional straight wall | Off-grid Poiseuille `u(y) = F·yy·(h−yy) / (2ν)` (integrate `ν d²u/dy²=−F` with u=0 at fractional walls). Refine cells-per-width; diffusive scaling `F ∝ 1/width²` fixes `u_peak_ref`. Expected p≈2. | slope ≥ 1.6, r² ≥ 0.9 | A1 | light |
| G2 | Sub-cell translation invariance | Same Poiseuille reference across wall_lo ∈ {0.20, 0.50, 0.80}. `u_peak_ref` invariant under translation; measured peak spread must be smooth (small, no jumps). Denominator = peak (peak-relative). | spread ≤ 6% of peak | A2 | light canary |
| G3 | qd=1/2 degeneracy with half-way BB | Already covered by `tests/bouzidi.rs::qd_half_records_are_bitwise_half_way_bounce_back`. Cross-reference only; not duplicated. | bitwise | A5 | (existing) |
| G4 | Off-grid Couette linear-profile recovery | `u(y) = U·yy/h` — linear, so bulk scheme contributes zero; residual isolates boundary error. Needs runtime moving-wall velocity setter for fractional Bouzidi walls. | L2rel ≤ 1e-3 | A1 | SPEC-GAP |
| G5 | Effective wall position vs τ | Half-way BB has known τ-dependent slip; Bouzidi's design promise is τ-independence. Fit quadratic zeros across τ ∈ {0.55…2.0}; drift bounded by channel width. | drift ≤ 2% of width | A3 | heavy |
| G6 | Rotational anisotropy | Diagonal-oriented channel vs axis-aligned at matched h, Ma. Needs a slanted-domain builder. | peak-relative ≤ 5% | A2 | SPEC-GAP |

## P2 — ENCODE

Written in-place (rather than dispatched, since P1 audit was authored in this
session — dispatch is used when P1 is done in a frontier session and P2 fans
out mechanically). Test file follows every P2 hard convention:

- References derived analytically in comments in the test file itself.
- Bands set from theory with 10x float-noise headroom (G1 asserts slope ≥ 1.6
  vs a theoretical 2; G2 asserts spread ≤ 6% vs a staircase-BB baseline of
  O(20–30%)).
- Two SPEC-GAP `#[ignore]`d stubs (G4, G6) carrying the analytic derivations.
- One heavy `#[ignore]`d test (G5) with the tau-sweep derivation embedded.
- CPU scalar backend only; no per-backend variants.
- Metrics: `common::metrics::l2_rel` and `common::metrics::order_fit` (no
  inline reimplementations).
- No PIN rows in this audit — the initial pass surfaced zero calibrated
  engine bugs. If a triage confirms one later it lands in this file with its
  ANOM id in the assert message.

## P3 — TRIAGE (calibration finding — recorded per the "log even test-side
disposition" rule)

**Initial G1 failure**: slope = −1.993, r² = 1.0000 — assertion failed with
message "convergence slope −1.993 < 1.6".

**Derive-before-blaming**: raw errors were `2.5e-3 → 1.3e-3 → 6.4e-4` as
`ny = 22 → 30 → 42`. Error decreases as resolution grows, so the physics IS
converging. r² = 1.0000 means the fit is on a clean line. The slope sign is
negative — a red flag pointing at the fit's x-axis, not the engine.

**Root cause (test-side)**: `order_fit` assumes `err ∝ h^p` with `h → 0`
as we refine. I passed `h = width` (channel width in lattice units), which
**grows** with resolution — so the fit landed at `err ∝ width^{-p}` with
`p ≈ 2`. The engine is second-order; the fit's x-axis was reversed.

**Disposition**: **test-fix** (S3, taxonomy: convergence-fit x-axis reversal).
Logged in `docs/qa/anomaly-log.md` as `ANOM-DRY-001`. Fix: pass `1/width` as
the mesh spacing so `h → 0` as ny grows. Added an inline comment in the
test naming the trap so the next audit's P1 does not repeat it. No engine
change.

**After fix**: slope = +1.993, r² = 1.0000, both light tests pass; three
ignored (2 SPEC-GAPs + 1 heavy).

**Calibration match**: this is one more test-side derivation trap on top of
the 6/6 from the accuracy_audit.rs pass on 2026-07-06. The dry run *found*
the same class of failure the Skill was built to catch and disposed of it
via P3 (test-fix) rather than a wasted engine-fix order. That IS the
justification for the P3 discipline.

## P4 — FIX

Not exercised — no engine bug in this pass. If G5 (τ-sweep) or G4/G6
(SPEC-GAP resolutions) surface an engine bug on a follow-up run, the
`references/order-templates.md` P4 template drives that dispatch verbatim.

## Verification

```
cargo test -p lbm-core --release --test accuracy_audit_bouzidi \
  -- --nocapture
# 2 passed; 0 failed; 3 ignored — the passing tests print measured slope
# 1.993 / r² 1.0000 and spread 5.69% (denominator: peak).
```

Full landing-gate check happens at commit time via
lbmflow-build-verify (see the commit for the qa/skill-accuracy-audit branch).
