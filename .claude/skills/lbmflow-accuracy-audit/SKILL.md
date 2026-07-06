---
name: lbmflow-accuracy-audit
description: >-
  Run the LBMFlow accuracy-audit loop: excavate every discretization/modeling
  approximation in a target subsystem, encode adversarial analytic tests for
  each, triage the failures (derive before blaming the engine), dispatch fixes,
  and pin regressions so they can never silently return. Use whenever the plan
  is to "audit accuracy", "hunt bent physics", "find approximation errors",
  "write adversarial accuracy tests for <subsystem>", "check the convergence
  order", "verify against the analytic solution", or a freshly landed physics
  feature needs a systematic accuracy shakedown beyond the VALIDATION.md bands.
  This Skill owns the five attack axes (A1-A5), the phase/model-tier split, the
  P2/P4 codex order templates, the triage decision table, and the shared
  agreement-metrics library. Do NOT use it to run the build gates themselves
  (lbmflow-build-verify), to dispatch/monitor codex mechanically
  (lbmflow-codex-dispatch owns the invocation), or for routine VALIDATION.md
  acceptance tests — this Skill looks for errors those bands cannot see.
---

# LBMFlow accuracy audit — excavate → encode → triage → fix → pin

Steady-state validation bands pass while the physics is bent. Field-proven on
2026-07-06: a force-path transient impulse error of 4/7·F (ANOM-P2-001) was
invisible to every steady gate (T2/T6/T11); a rotor geometry bug produced
mirror arms for odd blade counts while all fields looked plausible; an FD
shear reconstruction under-reported peaks by up to −35%. This Skill packages
the loop that found them. Its premise: **every approximation has an analytic
case where its error is first-order visible — find that case and encode it.**

## Phase overview and model tiers

The loop is deliberately split so only ideation/judgment needs a frontier
model; the mechanical phases are runnable verbatim by codex/Sonnet-class
agents from the templates in `references/order-templates.md`.

| Phase | What | Model tier | Output |
|---|---|---|---|
| P1 EXCAVATE | Enumerate approximations, pick the analytic kill-case for each | Frontier (Opus/Fable) — or Sonnet scaffolded by `references/axis-taxonomy.md` | Audit list (schema below) |
| P2 ENCODE | Write adversarial tests from the audit list, one template per item | codex / Sonnet (mechanical) | Test file(s) + first failure report |
| P3 TRIAGE | Classify each failure: test bug / engine bug / spec gap | Frontier (judgment) | anomaly-log entries + dispositions |
| P4 FIX | Implement engine fixes as file-disjoint orders; retighten pins | codex (mechanical) | Landed fixes, pins retightened |

Never collapse P2 into P1 ("write tests while excavating") — the audit list is
the reviewable artifact that makes P2 mechanical, and skipping it produces
tests whose reference physics nobody derived. Never let P4 start before P3 has
a disposition — on 2026-07-06 the **first six** adversarial failures were all
test-side physics-design errors, zero engine bugs.

## P1 — EXCAVATE (frontier model, or Sonnet + taxonomy)

For the target subsystem, enumerate **every** discretization and modeling
approximation it makes (interpolation stencils, link truncations, forcing
weights, boundary placement, moment corrections, sub-grid closures…). For
EACH approximation, name the analytic case where its error is
**first-order visible** — not a case where it merely contributes.

Work through the five attack axes in `references/axis-taxonomy.md` (each has a
field-proven worked example):

- **A1 Convergence order** — measured order vs theoretical order, r2-gated.
- **A2 Invariances** — Galilean, rotation/reflection, sub-cell translation,
  diffusive-scaling (same Re, different Ma).
- **A3 Functional form** — the error must lie ON the known curve;
  "small" is not a pass (e.g. bounce-back slip law vs tau).
- **A4 Transient fidelity** — Stokes I/II, acoustic cs + damping, startup
  impulse bookkeeping.
- **A5 Cross-path consistency** — force paths, backends, BC variants,
  2D/3D degeneracy.

Output: an **audit list**, one row per approximation:

```
{approximation, analytic reference + derivation sketch, expected order/band,
 axis (A1-A5), cost (light <1 s | heavy #[ignore])}
```

Done-check for P1: every approximation row has a *derivable* reference (a
formula you can re-derive in comments, not "compare with a finer grid"), and
at least one row per applicable axis. If no analytic case exists for an
approximation, that row becomes a SPEC-GAP candidate, not a dropped row.

## P2 — ENCODE (codex/Sonnet, mechanical)

Dispatch as a codex order using the **P2 template** in
`references/order-templates.md` (via lbmflow-codex-dispatch: own worktree,
`< /dev/null`, never shared with an implementation order). Hard conventions
the template enforces — review the returned tests against ALL of them:

1. **References derived analytically in comments.** Every reference formula is
   re-derived in the test file so the test is reviewable without consulting
   implementation internals.
2. **Freeze measured values with ~10x headroom.** Bands are set from the
   theory, not tuned to pass; where a measured constant must be frozen, leave
   ~10x headroom above float noise so refactors don't flake it.
3. **Ambiguity ⇒ `#[ignore]` + SPEC-GAP comment.** If the API cannot express
   the probe (e.g. compat has no runtime MovingWall velocity setter, so
   Stokes-II is unwritable), land an ignored test carrying the full analytic
   derivation and a `SPEC-GAP:` comment — never silently skip coverage.
4. **Known-anomaly pins assert the CURRENT WRONG value.** For a confirmed
   engine bug awaiting a fix, pin the wrong value (e.g. ANOM-P2-001 pins the
   uniform/field impulse ratio at 7/3) so the fix FAILS the pin loudly and
   forces the band to be retightened to the correct value.
5. **CPU-reference-only.** Audit tests target the CPU scalar reference; T14
   backend-equivalence gates all other backends transitively. Do not write
   per-backend audit variants.
6. **Metrics are pure functions from the shared library** —
   `crates/lbm-core/tests/common/metrics.rs` (Rust, source of truth) and
   `scripts/qa/metrics.py` (Python mirror): `l2_rel`, `linf_rel`, `order_fit`,
   `envelope_fit`, `phase_fit`, `monotonicity`, `curve_agreement`. API doc:
   `references/metrics-api.md`. No inline reimplementations.

Done-check for P2: the test file compiles, every audit-list row maps to a
test (or an ignored SPEC-GAP), and `cargo test --release` has been run once to
collect the initial pass/fail report — failures at this stage are *expected*
and go to P3, not to a fix order.

## P3 — TRIAGE (frontier model)

**Derive before blaming the engine.** For each failing test, re-derive the
reference physics independently before proposing an engine fix — check the
wall placement convention (half-way: wall surface at rim center + 0.5), the
forcing half-step (velocity moments include Guo F/2), the initialization
transient, and the asymptotic-regime assumption. Calibration: on 2026-07-06
the first 6 adversarial failures were ALL test-side physics-design errors.

Use the decision table in `references/triage-decision-table.md`. Severity is
the house taxonomy (S0 silently-wrong physics / S1 divergence-leak in a
supported config / S2 below-expected-accuracy / S3 minor). Disposition is one
of: **test-fix** (fix the test yourself, in-worktree) / **engine-fix order**
(becomes a P4 order) / **SPEC-GAP pin** (ignored test + contract note).

Every finding — including test-side errors, so the next audit doesn't repeat
them — is logged in `docs/qa/anomaly-log.md` (append-only):

```
{id: ANOM-<pass>-<nnn>, config, expected (with source), observed (excerpt),
 severity, disposition}
```

Done-check for P3: zero failing tests without a logged disposition.

## P4 — FIX (codex, mechanical)

Dispatch engine-fix orders using the **P4 template** in
`references/order-templates.md`. Rules:

- **File-disjoint worktree orders** per the PM/QA ownership map — check with
  the PM that the target file is not already churning in another order (e.g.
  collision-kernel fixes queue behind in-flight backend orders).
- **Regression pins stay after fixes.** A fix order must flip the
  known-anomaly pin to the correct value (retighten), never delete it. The
  retightened pin IS the regression test.
- Verification is the lbmflow-build-verify gate tier for the touched files —
  a physics change triggers the full `--include-ignored` tier plus a
  PHYSICS.md rationale entry.

## Merge-queue hard requirements (PM ruling — every P4 order encodes these)

All six are mandatory in P4 fix orders and their landing gates. The
`references/order-templates.md` P4 template already embeds them; changes to
this list happen only by PM ruling. Each row also names the recent incident
that motivates it, so a future audit cannot loosen a rule without
understanding what it costs.

1. **Landing gate = `cargo test --workspace --release --no-fail-fast`
   with UNPIPED exit code** — invoke as
   `cargo test --workspace --release --no-fail-fast; echo EXIT:$?`.
   Never pipe the gate through `tail`/`grep`/`awk` — the pipe eats the exit
   code and reports red as green. (Bit us for ~4h on a `backend_equiv`
   regression that fail-fast masked; a piped tail reported the run green.)
2. **GPU landing evidence is PM-run outside the codex sandbox.** In-sandbox
   Metal adapter access is intermittent (2 of 5 orders had it on
   2026-07-05). Fix orders touching GPU say "build + CPU gates only; report
   `BENCH/GPU-PENDING` for PM" — they do NOT gate on GPU benches.
3. **Adversarial/audit test orders NEVER share a worktree with implementation
   or fix orders.** Enforced by the P2 template stating it and by the
   dispatcher assigning worktrees.
4. **Band freeze = measured value + stated headroom.** Measured values are
   printed in every assert message (`println!` before `assert!`, and quoted
   in the assert body). Loosening a band needs a `docs/PHYSICS.md` rationale
   entry AND PM sign-off. **The denominator / normalization of every
   relative tolerance is stated in the assert message** — per-component
   relative vs force-scale relative flipped pass/fail on identical physics
   on 2026-07-05.
5. **Current-wrong-value pins carry their ANOM id in the assert message**
   (e.g. `"ANOM-P2-001 current wrong uniform/field impulse ratio"`) — so a
   fix order flips the pin instead of silently "fixing" the test. The pin
   is retightened to the correct value **in the same commit as the fix**;
   pin retighten is never deferred to a follow-up.
6. **Anomaly-log entry required before merge for every P3-confirmed finding**
   — including test-side dispositions. On 2026-07-06 the first 6 adversarial
   failures were 6/6 test-side; recording those in the log is what stops the
   next audit's P1 from stepping on the same rakes.

## Machine etiquette (before any heavy run)

- Check `TESTING_NOTES.md` (main repo copy) for an open
  `MEASUREMENT WINDOW **OPEN**` notice without a matching `**CLOSED**` — if a
  quiet window is open, hold all heavy runs and codex dispatches.
- `cargo test` ALWAYS `--release` (LBM is ~50x slower in debug).
- ≤ 5-8 compile-bound codex orders machine-wide — count the PM's fleet, not
  just your own orders.

## Metrics-library promotion rule

The metrics library today lives as a tests/common module in `lbm-core` +
a `scripts/qa` Python mirror (PM placement ruling). **Promote it to a shared
dev-dependency crate ONLY when a second crate independently needs the same
functions.** A speculative crate is abstraction-for-hypothetical-futures,
which CLAUDE.md's minimal-scope discipline bans. A drift-guard cross-check
test pins the Rust and Python semantics to each other (see
`references/metrics-api.md`).

## Worked example

`references/worked-example-bouzidi.md` — the full loop applied to the Bouzidi
curved-BC subsystem (interpolation order + sub-cell translation invariance),
recorded from the Skill's own dry run.

## Top failure modes (and the fix)

- **Blaming the engine first.** A failing adversarial test is a *claim*, not a
  finding. Re-derive the reference (P3) — most failures are test-side.
- **"Error is small" accepted as a pass (A3 violation).** Small-but-off-curve
  means the functional form is wrong and will bite at other parameters.
  Assert `curve_agreement`/`order_fit` against the curve, not a magnitude.
- **Order fit trusted without r2.** A slope through non-asymptotic points is
  noise; assert `fit.r2 >= 0.98` alongside the slope band.
- **Fix lands, pin deleted.** The pin must be retightened to the correct
  value in the same change — deleting it erases the only regression witness.
- **Coverage silently dropped.** An unwritable probe must land as
  `#[ignore]` + SPEC-GAP with the derivation in comments, so the gap is
  visible in `cargo test -- --ignored` listings.
- **Audit tests written against a non-reference backend.** CPU scalar only;
  everything else is T14's job.
