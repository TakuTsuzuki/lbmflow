# LBMFlow — Bioprocess-Specific CFD Core

**Product mission (2026-07-07 pivot; supersedes prior "commercial-grade
general-purpose LBM simulator" framing)**: a bioprocess-specific CFD core
whose QOI pipeline drives stirred-tank cell-culture / bioreactor process
decisions (Np, P/V, mixing time, gas holdup, d32, kLa, shear exposure,
oxygen exposure, cell/microcarrier damage risk, scale-up operating window).
Generic CFD parity is NOT a goal.

**Required reading before working on this repo**:

- [docs/BIOPROCESS_PIVOT.md](docs/BIOPROCESS_PIVOT.md) — the pivot and what
  it retracts / preserves.
- [docs/SPEC_BIOPROCESS_CORE.md](docs/SPEC_BIOPROCESS_CORE.md) — intended
  use, forbidden use, credibility tiers, QOI catalog, scenario schema,
  non-negotiables.
- [docs/PLAN.md](docs/PLAN.md) — BCFD-000..110 tickets, M0..M3 milestones,
  validation-driven development protocol, merge-queue rules, known traps.
- [docs/VALIDATION_BIOPROCESS.md](docs/VALIDATION_BIOPROCESS.md) — VB-01..VB-08
  validation groups (bioprocess acceptance criteria).

Keep this file in sync with [AGENTS.md](AGENTS.md) — the English mirror read
by codex and other non-Claude agents.

## Build & test

LBM is ~50× slower in debug — every gate uses `--release`.

```bash
cargo build --workspace --release
cargo test --workspace --release --no-fail-fast          # default gate
cargo test --release -- --include-ignored                # heavy validation (~5 min)
# WASM for the web GUI (lbm-wasm is outside the workspace):
wasm-pack build crates/lbm-wasm --target web --release --out-dir ../../web/src/engine/pkg
#   then delete pkg/.gitignore and commit pkg
cd web && npm run build                                  # GUI (tsc strict + vite)
./target/release/lbm presets run cavity                  # CLI smoke test (legacy demo)
```

**Landing gates use `--no-fail-fast` and UNPIPED exit code** (`; echo EXIT:$?`).
Piping through `tail`/`grep` eats the exit code — never do it. `cargo test`
without `--no-fail-fast` masks regressions after the first failing binary.

## Repository map

- `crates/lbm-core` — the single core (V2): D2Q9/D3Q19/D3Q27 lattices; CPU
  scalar/SIMD backends; wgpu GPU backend (feature `gpu`, off by default —
  the workspace suite does NOT cover it, CI runs `--features gpu` on GPU
  hosts); MPI halo exchange (feature `mpi`, off by default, needs a native
  MPI toolchain — scripts/test_mpi.sh); legacy 2D facade in `compat/`.
- `crates/lbm-scenario` — JSON scenario schema + runner. The legacy
  `Scenario` type continues to parse; `BioprocessScenario` is added
  alongside it (BCFD-003).
- `crates/lbm-cli` — `lbm` binary: presets (legacy demos, emit
  "not bioprocess decision-grade" warning), gallery, schema, scenario run,
  MCP server (7 tools incl. async start_run / run_status / list_runs; new
  bioprocess-tool surface added per BCFD-092).
- `crates/lbm-wasm` — WASM bindings for the web GUI (outside the workspace).
- `crates/lbm-gpu-proto` — wgpu evaluation prototype (outside the workspace;
  measurement record, superseded by the in-core `gpu` module).
- `web/` — TypeScript GUI (Vite; engine WASM committed under
  `web/src/engine/pkg`). Not on the bioprocess critical path until BCFD-081
  report generator is useful.

Version: crates are `0.2.0-bioprocess.0`. The pre-pivot final state is
`git tag v1-lbm-general-final`.

## Docs index (read on demand)

Bioprocess product docs (living):

- [BIOPROCESS_PIVOT.md](docs/BIOPROCESS_PIVOT.md) — pivot announcement,
  retracted vs preserved claims.
- [SPEC_BIOPROCESS_CORE.md](docs/SPEC_BIOPROCESS_CORE.md) — intended use,
  tiers, QOI catalog, scenario schema, non-negotiables.
- [VALIDATION_BIOPROCESS.md](docs/VALIDATION_BIOPROCESS.md) — VB-01..VB-08.
- [CREDIBILITY_BIOPROCESS.md](docs/CREDIBILITY_BIOPROCESS.md) — calibration/
  holdout separation, evidence gate policy.
- [MODEL_RISK_MATRIX.md](docs/MODEL_RISK_MATRIX.md) — per-model risk table.
- [PLAN.md](docs/PLAN.md) — BCFD-000..110 tickets, M0..M3 milestones,
  development protocol, merge-queue rules, known traps.
- [LIMITATIONS.md](docs/LIMITATIONS.md) — machine-readable capability status
  and unsupported combinations.

Preserved engineering references (living):

- [PHYSICS.md](docs/PHYSICS.md) — physics decisions and experiment log
  (update whenever you change physics; carried forward across the pivot).
- [ARCHITECTURE_V2.md](docs/ARCHITECTURE_V2.md) — dimension × lattice ×
  precision × backend × partition design (still describes the code).
- [KERNEL_EXTENSION_POINTS.md](docs/KERNEL_EXTENSION_POINTS.md) — kernel
  extension contracts.
- [MPI_GUIDE.md](docs/MPI_GUIDE.md) · [CLUSTER_OPTIONS.md](docs/CLUSTER_OPTIONS.md)
  · [CLUSTER_RUNBOOK.md](docs/CLUSTER_RUNBOOK.md) — distributed runs.
- [REQ_STIRRED_REACTOR.md](docs/REQ_STIRRED_REACTOR.md) — pre-pivot
  requirements text, still useful as bioprocess reference for BCFD-020.
- [T15_5_CAVITY3D_REFERENCE.md](docs/T15_5_CAVITY3D_REFERENCE.md) —
  3D cavity reference data (single-phase validation reference).

Archive:

- [archive/2026-07-07-pivot/](docs/archive/2026-07-07-pivot/) — pre-pivot
  PLAN, VALIDATION, LIMITATIONS; T1..T18 matrix; M-A..M-F track; R-Phase
  spec; V&V ledger; whitepaper; claims-ledger; competitor analysis. Read
  only when you need the pre-pivot history.

## Documentation principles (preserved)

- **One home per fact.** Each topic has one owning doc; others link to it.
  When facts change, update the owning doc in place — never fork a new file
  or append a contradicting section.
- **Index or archive.** Every file directly under `docs/` has a one-line
  entry in the Docs index above; update docs and index (this file AND
  AGENTS.md) in the same commit. Superseded or one-shot docs (handoffs,
  dated audits, drafts) move to `docs/archive/<date>-<slug>/` with a pointer.
- **Sections over files.** Only a new subsystem or track justifies a new
  top-level doc; findings, reports, and experiment records go into the
  owning doc.
- **Lifecycle in the header**: living (updated in place), frozen spec
  (change = recorded revision), or snapshot (dated, never edited).
- **Docs move with the change.** A behavior-altering change updates the
  affected docs (physics + rationale + experiment results → PHYSICS.md;
  claims → LIMITATIONS.md; plan → PLAN.md).
- **English, absolute dates, repo-relative links.** Write 2026-07-07, never
  "today".

## Team & conventions (preserved)

- Fable is PM. **Implementation workhorse = codex CLI, fanned out in parallel
  (user directive 2026-07-05).** One order = one focused BCFD-item bundle =
  one dedicated git worktree + background `codex exec -C <worktree>`; the PM
  merges landed branches in dependency order. Bundle same-file items into
  ONE order — parallelize across files, never within one. Concurrency is
  machine-bounded: ~5-8 compile-bound orders at once on this M5 Max.
- **Validation tests are written adversarially** (by codex or Opus/Sonnet,
  from the spec in VALIDATION_BIOPROCESS.md); a test order and an
  implementation order never share a worktree.
- codex invocation:
  `codex exec --sandbox workspace-write --skip-git-repo-check "<task>" < /dev/null`
  (model gpt-5.5). **`< /dev/null` is mandatory** — with a pipe as stdin,
  codex waits for EOF and hangs forever. Monitor progress via
  `~/.codex/sessions/<date>/rollout-*.jsonl`.
- **ALL artifacts in English (user directive 2026-07-05)**: code,
  identifiers, commits, docs, and user-facing strings (GUI / CLI / errors).
  Conversation with the user may remain Japanese — the product is English.
- Commit at each phase completion; never commit red tests (explicit WIP,
  flagged in the commit message, is the only exception).

## Working discipline (every agent on this repo)

- **Physical rigor is the prime directive; ad-hoc physics is BANNED (user
  directive 2026-07-06).** Every physical behavior anywhere in the stack —
  core, scenario, CLI, examples, GUI — is either resolved from the governing
  equations or a literature-backed closure with a recorded derivation,
  validity domain, and its own validation test (PHYSICS.md entry mandatory).
  Banned outright: constants calibrated to pass an acceptance band; branches
  keyed to sample/case identity; clamps or caps that silently absorb
  transport; decorative physics terms. If a gate cannot be met without such
  a hack, STOP and report — the spec gets revised, not the physics faked.
  The executable procedure (provenance decision table, ban-list greps,
  two-layer gates, stop-rule template, behavior-review record) is
  `.claude/skills/lbmflow-physics-discipline`.
- **Behavior-validity review (user directive 2026-07-06).** After every
  experiment/demo run, before reporting: judge whether the OBSERVED behavior
  — spatial patterns, trends, signs, not just gated metrics — is physically
  plausible; name the dominant mechanism; separate resolved physics from
  closures and boundary artifacts; record the review in PHYSICS.md or the
  track's findings file. A metric passing its band does NOT validate a
  pattern no band covers.
- **Finding routing (V&V loop).** Core-engine defect → phenomenon report +
  data package (scenario JSON, exported fields, metrics, repro command) to
  the core-engine session. Demo/example defect → PM dispatches a codex order.
  Spec defect → spec revision, rationale in PHYSICS.md.
- **Evidence-based progress.** Match every "done" claim against a tool
  result from THIS session (test output, file diff, run log); "done" means
  the relevant Build & test gate passed HERE — a codex order finishing is
  not evidence its branch is green. Report unverified as unverified, skipped
  as skipped, failures with their output. Fabricated progress is the worst
  possible failure.
- **Minimal scope.** Implement exactly what the BCFD ticket asks: no
  drive-by refactors, no helpers for one-off operations, no abstractions
  for hypothetical futures, no defensive code for impossible states
  (validate at system boundaries only: scenario JSON, CLI args, user input).
- **Finish, don't announce.** Never end a turn on a plan or "I'll now do X"
  — do X, then end. Stop only when the task is complete or blocked on input
  only the user can provide. When unsure between options, pick one and
  proceed; don't re-litigate settled decisions.

## Bioprocess-specific discipline (added by pivot)

- **QOI provenance is mandatory.** Every QOI ships `units`, `method`,
  `time_window`, `averaging_region`, `source_fields`, `validation_tier`
  (BCFD-080). Missing metadata → serialisation fails. `max` alone is never
  the report — percentile distributions (P50/P90/P95/P99) plus
  fraction-above-threshold are required for shear, exposure, and
  stress-derived QOIs.
- **Capability registry says what's supported (BCFD-002).** Every
  bioprocess capability carries a status
  (`Unsupported`/`Experimental`/`Engineering`/`EvidenceBlocked`/`EvidenceReady`).
  Unsupported combinations fail with a structured error, never silently
  fall back. `lbm capabilities --json` is the machine-readable truth.
- **Evidence-tier claims are mechanically gated (BCFD-091).** Reports
  cannot carry the evidence label unless `EvidenceGate` returns ready:
  validation matrix status pass + calibration/holdout separated + mesh /
  time-step sensitivity + QOI uncertainty interval + limitation report.
- **No production Shan-Chen.** Shan-Chen SCMP/MCMP remains in
  `crates/lbm-core` as an engineering / demo capability. The production
  gas-liquid path is conservative Allen-Cahn phase field (BCFD-040..048).

## Core design invariants (breaking these kills the validation suite)

Preserved unchanged across the pivot — these are properties of the code,
not the product framing:

- The single core is `crates/lbm-core`. The legacy V1 API lives on as
  `lbm_core::compat`, used by the scenario / CLI / wasm 2D paths.
- **D2Q9 direction ordering** (single source of truth: the `Lattice` impl
  in `lattice.rs`): 0:(0,0), 1:(1,0), 2:(0,1), 3:(-1,0), 4:(0,-1), 5:(1,1),
  6:(-1,1), 7:(-1,-1), 8:(1,-1).
- **f layout** is q-major SoA with halo padding (`fields.rs`):
  `f[q*plane + cell]`, `cell = z·(nx·ny) + y·nx + x` — identical to the GPU
  coalescing assumption; never exposed in the public API.
- **Step order**: collide → halo exchange → streaming → open-boundary BCs →
  boundary moments correction (CpuSimd fuses collide+stream+moments in
  `step_band`). Any change to pass structure or storage order must pass the
  backend-equivalence and partition-invariance gates before it lands.
- **Walls**: edges are a 1-cell solid rim; wall surfaces are half-way
  (midpoint between rim center and fluid center).
- **Guo forcing F/2 correction** included in velocity moments — `sim.ux()`
  etc. return physical velocity.
- **`tau = 3·nu + 0.5`** (cs² = 1/3).
- **Deviation storage** (`f - w`) is the f32-precision scheme: bounce-back
  is linear so it's invariant in deviation space; Zou-He needs the constant
  term folded in; `ρ = 1 + Σ dev`.
- **Backend contract**: fused CpuSimd is equivalent to CpuScalar to
  ≤1e-5 rel on TGV / cavity; the GPU backend to the same tolerance on the
  subset of scenarios it supports. Bit-reproducibility across partition
  counts is a gate.
