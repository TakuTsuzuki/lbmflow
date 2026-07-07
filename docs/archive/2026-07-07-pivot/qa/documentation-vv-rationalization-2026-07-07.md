# Documentation and V&V Rationalization Proposal

Lifecycle: snapshot (2026-07-07) -- read-only implementation/documentation cross-check.

This report compares the current implementation surface, validation tests, QA
registers, and release-facing docs. It is not a test run. It records
reproducible documentation and V&V-quality findings, then proposes a smaller
source-of-truth structure that preserves the useful audit evidence without
letting stale snapshots drive current decisions.

## Method

Commands used:

```bash
rg --files docs crates web
rg -n "not yet implemented|P2 in progress|SPEC-ONLY|UNSAFE-CLAIM|ANOM-P4-021|ANOM-P4-022" docs crates
nl -ba docs/VALIDATION.md
nl -ba docs/qa/VV_TRACEABILITY.md
nl -ba docs/qa/code-to-spec-diff.md
nl -ba docs/qa/anomaly-log.md
nl -ba crates/lbm-scenario/src/lib.rs
```

Scope inspected:

- `docs/PLAN.md`, `docs/VALIDATION.md`, `docs/LIMITATIONS.md`,
  `docs/PHYSICS.md`, `README.md`, `AGENTS.md`, `CLAUDE.md`.
- QA docs under `docs/qa/`, especially `VV_TRACEABILITY.md`,
  `code-to-spec-diff.md`, and `anomaly-log.md`.
- Implementation entry points in `crates/lbm-scenario/src/lib.rs`,
  `crates/lbm-cli/src/{runner,capabilities,verify}.rs`, and selected
  regression tests in `crates/lbm-core/tests/`.

No cargo validation gate was run for this report.

## Reproducible Findings

### DOC-VV-001 -- Stale QA snapshots contradict newer source-of-truth docs

Evidence:

- `docs/qa/VV_TRACEABILITY.md:42-44` says T16 is internally inconsistent
  because `VALIDATION.md` still says FP16 is "not yet implemented".
- Current `docs/VALIDATION.md:289-307` instead says T16 is implemented,
  gated, and records frozen FP16 bands.
- `docs/qa/code-to-spec-diff.md:71-75` lists ANOM-P4-021 and ANOM-P4-022 as
  open priority fixes.
- Current `docs/qa/anomaly-log.md:323-346`, `docs/qa/anomaly-log.md:358-374`,
  and `docs/qa/anomaly-log.md:496-505` mark ANOM-P4-021 and ANOM-P4-022
  fixed/closed with regression evidence.
- `docs/PHYSICS.md:486-499` and `docs/PHYSICS.md:827-838` also record the
  resolved Guo/Zou-He and force-field timing contracts.

Impact: a reader following the QA docs first can route already-fixed items as
active P0 defects and miss the real remaining work.

Proposal: mark generated QA cross-check files with one of
`current`, `superseded-by:<path>`, or `snapshot-only` in the header. Move
superseded generated reports to `docs/archive/` or add a first-screen banner
linking to the active owner. Do not keep stale "priority action items" in
living QA files after the anomaly log has closed them.

### DOC-VV-002 -- The docs index underspecifies the actual validation scope

Evidence:

- `AGENTS.md:41-42` says `VALIDATION.md` covers "T1...T15.x".
- Current `docs/VALIDATION.md:309-397` contains T17 and T18 status and
  acceptance sections.
- `AGENTS.md:63-65` summarizes the whole QA tree as generic `VV_*` docs and
  anomaly logs, but the tree currently contains 27 QA files, including
  dated audits, generated coverage matrices, hard-case inventories, and
  machine-readable pass data.

Impact: the top-level onboarding path does not tell an agent where the live
coverage boundary is. This encourages duplicate audits and stale file reads.

Proposal: update `AGENTS.md` and `CLAUDE.md` docs index to say
`VALIDATION.md` owns T1-T18 acceptance, while `docs/qa/` is divided into
active registers, generated snapshots, and method/process docs. Keep the
index high-level, but name the active registers explicitly:
`anomaly-log.md`, `VV_TRACEABILITY.md` or its replacement, and
`VV_MASTER_PLAN.md`.

### DOC-VV-003 -- T18 status is split across plan, validation, traceability, and examples

Evidence:

- `docs/PLAN.md` states D-track P0-P2 are done and T18.1/.2/.3 are green.
- `docs/VALIDATION.md:354` still headlines T18 as "spec wired, P2 in progress".
- `docs/VALIDATION.md:390-397` correctly keeps T18.4/T18.5 as example/stub
  future work.
- `docs/qa/VV_TRACEABILITY.md` records T18.4 as no committed CLI example path,
  but `crates/lbm-cli/examples/dispersed_seeding/main.rs` exists in the tree.

Impact: the same feature reads as done, in progress, and missing depending on
the document. The useful distinction is T18.1-T18.3 validated, T18.4
example-level, T18.5 stub.

Proposal: change the T18 header in `VALIDATION.md` to "T18.1-T18.3 validated;
T18.4 example-level; T18.5 stub". Regenerate or retire the traceability row
that says no example path exists.

### DOC-VV-004 -- User-facing collision naming is still mixed

Evidence:

- `README.md:146` says the scenario path exposes `bgk`, `trt`, and
  `cumulant`.
- `crates/lbm-scenario/src/lib.rs:292-334` shows the canonical scenario value
  is `central_moment`; `cumulant` is accepted as `DeprecatedCumulantAlias`.
- `crates/lbm-scenario/src/lib.rs:343-360` maps both names to the central
  moment path but panics for 2D compat because central moment is a 3D/native
  route.

Impact: new users are pointed at a deprecated alias and may infer true Geier
cumulants where the implementation is a central-moment operator.

Proposal: make `central_moment` the only advertised name in README and
capability docs. Mention `cumulant` only as a backward-compatible deprecated
alias in schema reference material.

### DOC-VV-005 -- Valuable code-to-spec findings are not reconciled after fixes

Evidence:

- `docs/qa/code-to-spec-diff.md:39-63` records useful arithmetic drifts,
  including silent omega ceiling, Bouzidi order fallback, particle sampling
  semantics, pass-order drift, and backend trait drift.
- Some high-severity rows in the same file are now fixed, as shown by
  `docs/qa/anomaly-log.md:496-505`.

Impact: the report is too valuable to delete, but too stale to treat as live.

Proposal: split code-to-spec output into two artifacts:

- A frozen dated snapshot in `docs/archive/`.
- A live residual register containing only unresolved doc/code drifts, each
  with `status`, `owner doc`, `owner code path`, and `next gate`.

### DOC-VV-006 -- V&V coverage is broad, but the public trust boundary is not generated

Evidence:

- The implementation has many focused gates: `accuracy_audit_*.rs`, T13/T14
  equivalence tests, T15/T16/T18 tests, `lbm verify`, and `lbm capabilities`.
- `docs/LIMITATIONS.md` is the intended release-facing trust boundary, but
  capability facts are manually repeated in README, PLAN, VALIDATION,
  traceability, and proposals.

Impact: product claims drift faster than tests. The repo already has the data
needed for a generated trust boundary, but still depends on manual edits.

Proposal: make `lbm capabilities --json` and a small QA manifest the canonical
machine-readable source for product-path capabilities. Generate the README
capability matrix and the `LIMITATIONS.md` capability rows from that manifest,
or at minimum add a CI drift test that compares them.

## Improvement Plan

### P0 -- Stop stale documents from looking current

1. Add a header status to every `docs/qa/*.md`: `living-register`,
   `generated-snapshot`, `method-doc`, or `superseded`.
2. Move or banner stale generated reports:
   - `code-to-spec-diff.md`: preserve as snapshot; extract unresolved residuals.
   - `VV_TRACEABILITY.md`: regenerate from current tests or mark as
     2026-07-06 snapshot.
3. Fix direct contradictions:
   - T18 header in `VALIDATION.md`.
   - README collision naming.
   - `AGENTS.md`/`CLAUDE.md` validation index T1-T18.

### P1 -- Define one owner per fact

Recommended ownership:

| Fact class | Owner | Consumers |
|---|---|---|
| Current release capability and unsupported combinations | `docs/LIMITATIONS.md` plus `lbm capabilities` | README, GUI, MCP help |
| Validation acceptance bands | `docs/VALIDATION.md` | tests, traceability, claims ledger |
| Governing equations, closures, validity domains | `docs/PHYSICS.md` | validation docs, code comments |
| Roadmap and work queue | `docs/PLAN.md` | PM handoffs, proposals |
| Known anomalies and dispositions | `docs/qa/anomaly-log.md` | PLAN, PHYSICS, tests |
| Historical/generated audits | `docs/archive/` or dated `docs/qa/` snapshots | forensic reference only |

Rule: no document outside the owner should restate status in prose unless it
also names the owner and a checked date.

### P2 -- Make traceability reproducible

Build a small script or `cargo xtask` that emits:

- Validation ID.
- Owning spec file/section.
- Test file(s).
- Command.
- Run class: default, ignored, gpu, mpi, cluster, manual.
- Status: validated, verified-only, spec-only, bench-pending, unsafe-claim.

The script should fail CI when:

- A `VALIDATION.md` T-row has no traceability row.
- A traceability row points at a missing file.
- README or LIMITATIONS advertises a capability absent from `lbm capabilities`.
- `AGENTS.md` and `CLAUDE.md` docs indexes diverge.

### P3 -- Prioritize remaining V&V/user-quality work by risk

1. T17 full stirred-reactor claims: keep all full-coupled stirred, aeration,
   scalar/reaction, and relaxation-mode claims blocked until executable
   VR-STR gates exist.
2. Dynamic multiphase surface-tension behavior: keep ANOM-P4-014/P4-017
   characterized until the pressure-tensor audit resolves which sigma is
   valid for statics, menisci, and retracting rims.
3. Rotating boundary/penalization: keep ANOM-P4-010 and the Taylor-Couette
   heavy gate visible; do not promote impeller claims beyond current IBM
   subsystem evidence.
4. Bouzidi moving-wall audit: preserve the current wrong-value pin
   `ANOM-L1_7-001` and route a focused fix/order only after the wall-velocity
   derivation is reviewed.
5. Product surface quality: promote runtime `max_Ma`, grid-Re, mass/momentum
   drift, and selected backend/lattice/collision/storage into scenario
   manifests so users see when a run is outside the trust boundary.

## Proposed Acceptance Criteria

- A new reader can answer "what is supported today?" from README ->
  LIMITATIONS without reading PLAN.
- A validation owner can answer "what tests prove T18.2?" from
  VALIDATION -> traceability without reading archived handoffs.
- Every open anomaly in `anomaly-log.md` has one of: regression pin,
  ignored SPEC-GAP test, or explicit owner doc section.
- No current QA doc lists ANOM-P4-021 or ANOM-P4-022 as open.
- `cargo test --workspace --release` is not required to update docs, but the
  docs drift checks above run in the default CI tier.
