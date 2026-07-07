# ORDER -- Documentation / QA Rationalization and Cleanup

Commissioned 2026-07-07 after implementation-vs-documentation audit. Goal:
make LBMFlow's documentation trustworthy for agents and users by assigning one
owner per fact, regenerating traceability from current implementation, and
cleaning old QA reports so resolved/stale findings cannot be mistaken for the
active queue.

This is a docs/tooling order. Do not change solver physics. If the cleanup
uncovers a new physics defect, record it in `docs/qa/anomaly-log.md` and stop
at a routing package; do not fix physics in this order.

## Required context

Read these first:

- `AGENTS.md` and `CLAUDE.md` documentation principles.
- `docs/qa/documentation-vv-rationalization-2026-07-07.md`.
- `docs/VALIDATION.md`.
- `docs/LIMITATIONS.md`.
- `docs/PHYSICS.md`.
- `docs/qa/anomaly-log.md`.
- `docs/qa/VV_TRACEABILITY.md`.
- `docs/qa/code-to-spec-diff.md`.

Respect the current owner hierarchy unless this order explicitly changes it:

| Fact class | Owner |
|---|---|
| Current release capability and unsupported combinations | `docs/LIMITATIONS.md` plus `lbm capabilities` |
| Validation acceptance bands | `docs/VALIDATION.md` |
| Governing equations, closures, validity domains | `docs/PHYSICS.md` |
| Roadmap and work queue | `docs/PLAN.md` |
| Known anomalies and dispositions | `docs/qa/anomaly-log.md` |
| Historical/generated audits | `docs/archive/` or dated `docs/qa/` snapshots |

## Task A -- Confirm P0 cleanup state

Verify whether commit `a41c3b4` or equivalent changes are present:

1. `AGENTS.md` and `CLAUDE.md` docs index says `VALIDATION.md` covers T1-T18.
2. `docs/VALIDATION.md` T18 heading distinguishes:
   - T18.1-T18.3 validated,
   - T18.4 example-level,
   - T18.5 stub.
3. `README.md` advertises `central_moment` as canonical and `cumulant` only as
   a deprecated alias.
4. `docs/qa/VV_TRACEABILITY.md` is clearly marked as a 2026-07-06 snapshot if
   it has not been regenerated.
5. `docs/qa/code-to-spec-diff.md` has a first-screen banner that ANOM-P4-021
   and ANOM-P4-022 are closed in `anomaly-log.md`, so its priority list is not
   the active queue.

If any item is missing, apply the smallest docs-only patch to restore it.

## Task B -- Clean old QA documents as resolved/stale

Sweep every file under `docs/qa/` and classify it in the first 10 lines as one
of:

- `Lifecycle: living-register` -- current active register; must not contain
  stale resolved work as active queue.
- `Lifecycle: generated-snapshot (YYYY-MM-DD)` -- historical output; can
  contain old findings, but must not claim to be current.
- `Lifecycle: method-doc` -- process guidance.
- `Lifecycle: superseded` -- no longer active; header must name the current
  owner or replacement.

Cleanup rule:

- Resolved findings must be removed from active queue sections.
- Historical evidence may remain only in generated snapshots or archive files.
- If a file is mostly stale generated output, either:
  - move it to `docs/archive/` with its date in the filename, or
  - keep it in place but add a prominent snapshot/superseded banner and a link
    to the current owner.
- Do not delete unique anomaly evidence unless the same evidence is preserved
  in `docs/qa/anomaly-log.md`, `docs/PHYSICS.md`, or an archived snapshot.
- Do not leave ANOM-P4-021 or ANOM-P4-022 in any active "open", "priority",
  "highest risk", or "to fix" list. They are closed by the current
  `anomaly-log.md`.

Expected direct cleanup targets:

- `docs/qa/code-to-spec-diff.md`: convert to snapshot or archive; extract only
  unresolved residuals into a new active register.
- `docs/qa/VV_TRACEABILITY.md`: regenerate or keep as explicit snapshot; remove
  any current-state language if not regenerated.
- Dated or generated sweep outputs: add lifecycle headers if absent.
- Any `docs/qa/ext-codex/` or imported external reports: keep as external
  snapshots unless independently verified and promoted into current owner docs.

## Task C -- Extract a live residual register

Create `docs/qa/code-to-spec-residuals.md` as the live follow-up register for
unresolved doc/code drift. Each row must include:

| ID | Source row/report | Current status | Owner doc | Owner code path | Required action | Gate |
|---|---|---|---|---|---|---|

Populate it by reconciling `docs/qa/code-to-spec-diff.md` against:

- `docs/qa/anomaly-log.md`.
- `docs/PHYSICS.md`.
- Current tests under `crates/lbm-core/tests/`.
- Current implementation paths cited by the source row.

Rules:

- Closed items become `closed` with a pointer to the closing anomaly-log row;
  do not put them in the active-action section.
- Still-open S2/S3 items become actionable residuals.
- If a row cannot be verified cheaply, mark `needs-recheck` and include the
  exact grep/file evidence needed for the next agent.
- The register should start with a short "Active residuals only" summary, so
  readers do not have to parse the archived full report.

Seed residuals to re-check from the prior report:

- Silent central-moment omega ceiling / diagnostic gap.
- Bouzidi second-order claim vs fallback branches.
- Particle sampler near-wall and out-of-grid semantics.
- SC surface-tension referee mismatch, unless fully closed elsewhere.
- Pass-order docs vs actual backend pass order.
- Backend trait docs vs actual trait surface.
- Gravity composition debug/type-state guard.
- Per-cell omega odd-relaxation latent asymmetry.

Do not assume these are still open. Verify each against current docs/code.

## Task D -- Make traceability reproducible

Add a small script or `cargo xtask` only if it is clearly cheaper than manual
regeneration. Minimal acceptable output is a checked-in generated markdown or
CSV with:

- Validation ID.
- Owning spec section.
- Test file(s).
- Command.
- Run class: default, ignored, gpu, mpi, cluster, manual.
- Backend/precision/lattice where relevant.
- Status: validated, verified-only, spec-only, bench-pending, unsafe-claim.

Required drift checks:

- Every T-section in `docs/VALIDATION.md` has a traceability row.
- Every traceability row points to an existing file or explicitly says
  "manual/external".
- README/LIMITATIONS advertised capabilities match `lbm capabilities` output
  or have a named reason for being docs-only.
- `AGENTS.md` and `CLAUDE.md` docs indexes stay synchronized.

If a full generator is too large for one order, deliver:

1. A manually refreshed `VV_TRACEABILITY.md`.
2. A script stub that validates file existence and T-section coverage.
3. A follow-up order with exact remaining work.

## Task E -- Update public trust boundary without duplicating facts

Review README, `docs/LIMITATIONS.md`, `docs/PLAN.md`, and
`docs/VALIDATION.md` for repeated current-state prose. Replace repeated status
with links where possible.

Rules:

- README may summarize capability for users, but `docs/LIMITATIONS.md` owns
  current unsupported combinations.
- PLAN may describe roadmap and queue, but must not be the release trust
  boundary.
- VALIDATION owns acceptance bands, not product marketing status.
- QA snapshots must not be cited as current evidence unless refreshed in this
  order.

## Verification commands

Run at least:

```bash
rg -n "ANOM-P4-021|ANOM-P4-022" docs/qa
rg -n "not yet implemented|P2 in progress|T1.*T15|cumulant" README.md AGENTS.md CLAUDE.md docs/VALIDATION.md docs/qa
rg -n "^Lifecycle:" docs/qa
git diff --check
```

Interpretation:

- ANOM-P4-021/022 may appear in snapshots or closed-history sections, but not
  in active open/priority lists.
- "cumulant" may appear only as deprecated alias, historical text, or cited
  snapshot evidence; user-facing current docs should prefer `central_moment`.
- Any QA file without a lifecycle header must be intentionally fixed or listed
  in the final report as a follow-up.

If scripts or Rust tooling are changed, also run the relevant targeted tests.
For docs-only edits, cargo validation is not required; explicitly report that
it was skipped because no code behavior changed.

## Deliverables

1. Cleaned QA docs with lifecycle headers and stale/resolved findings removed
   from active sections.
2. `docs/qa/code-to-spec-residuals.md` containing only active unresolved
   doc/code residuals.
3. Refreshed or clearly-snapshotted `docs/qa/VV_TRACEABILITY.md`.
4. Any drift-check script or documented follow-up if automation is too large.
5. Final report listing:
   - files changed,
   - old QA docs archived/superseded,
   - resolved items removed from active queues,
   - active residual count by severity,
   - verification commands and results.

## Stop conditions

Stop and ask the PM only if:

- A QA document contains unique evidence for an unresolved S0/S1 anomaly and
  there is no clear owner document to preserve it.
- Current implementation contradicts both `docs/PHYSICS.md` and
  `docs/qa/anomaly-log.md`.
- Cleanup would require deleting large external artifacts whose provenance is
  unclear.

Otherwise proceed with conservative cleanup and leave precise follow-up rows.
