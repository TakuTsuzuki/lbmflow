# LBMFlow — Lattice Boltzmann Method Fluid Simulator

Commercial-grade LBM simulator: Rust core + TypeScript GUI + Agent mode.
**Required reading**: [docs/PLAN.md](docs/PLAN.md) (phase plan / team structure) ·
[docs/VALIDATION.md](docs/VALIDATION.md) (validation specs = acceptance criteria).
Keep this file in sync with AGENTS.md (English mirror read by codex and other
non-Claude agents).

## Build & test

```bash
cargo build --workspace --release
cargo test --workspace --release          # normal suite (LBM is ~50x slower in debug — ALWAYS --release)
cargo test --release -- --include-ignored # full validation incl. heavy benches (~5 min)
# WASM for the web GUI (lbm-wasm is outside the workspace):
wasm-pack build crates/lbm-wasm --target web --release --out-dir ../../web/src/engine/pkg
#   then delete pkg/.gitignore and commit pkg
cd web && npm run build                   # GUI (tsc strict + vite)
./target/release/lbm presets run cavity   # CLI smoke test
```

## Repository map

- `crates/lbm-core` — the single core (V2): D2Q9/D3Q19 lattices; CPU
  scalar/SIMD backends; wgpu GPU backend (feature `gpu`, off by default —
  the workspace suite does NOT cover it, CI runs `--features gpu` on GPU
  hosts); MPI halo exchange (feature `mpi`, off by default, needs a native
  MPI toolchain — scripts/test_mpi.sh); legacy 2D facade in `compat/`
- `crates/lbm-scenario` — JSON scenario schema + runner (2D compat path)
- `crates/lbm-cli` — `lbm` binary: presets, gallery, schema, scenario run,
  MCP server (7 tools incl. async start_run / run_status / list_runs)
- `crates/lbm-wasm` — WASM bindings for the web GUI (outside the workspace)
- `crates/lbm-gpu-proto` — wgpu evaluation prototype (outside the workspace;
  measurement record, superseded by the in-core `gpu` module)
- `web/` — TypeScript GUI (Vite; engine WASM committed under `web/src/engine/pkg`)

## Docs index (read on demand)

- [PLAN.md](docs/PLAN.md) — milestones M-A…M-F, current queue ·
  [VALIDATION.md](docs/VALIDATION.md) — acceptance criteria (T1…T18)
- [LIMITATIONS.md](docs/LIMITATIONS.md) — release-facing trust boundary
  (what is NOT supported/validated today; keep in sync with claims)
- [PHYSICS.md](docs/PHYSICS.md) — physics decisions + experiment log
  (update whenever you change physics)
- [ARCHITECTURE_V2.md](docs/ARCHITECTURE_V2.md) — dimension × lattice ×
  precision × backend × partition design
- [SOLVER_IMPROVEMENT_SPEC.md](docs/SOLVER_IMPROVEMENT_SPEC.md) — R-Phase spec ·
  [KERNEL_EXTENSION_POINTS.md](docs/KERNEL_EXTENSION_POINTS.md) — B-8 kernel extension contracts
- [REQ_STIRRED_REACTOR.md](docs/REQ_STIRRED_REACTOR.md) — M-F requirements ·
  [T15_5_CAVITY3D_REFERENCE.md](docs/T15_5_CAVITY3D_REFERENCE.md) — 3D cavity reference data
- [DISPERSED_DEPOSITION.md](docs/DISPERSED_DEPOSITION.md) — D-track dispersed-phase
  deposition tool (frozen spec, P0–P4, CR-1/2/3, acceptance = T18)
- [PERFORMANCE.md](docs/PERFORMANCE.md) / [GPU_EVALUATION.md](docs/GPU_EVALUATION.md) /
  [BENCH_COMPARISON_DRAFT.md](docs/BENCH_COMPARISON_DRAFT.md) — perf measurements
- [MPI_GUIDE.md](docs/MPI_GUIDE.md) / [CLUSTER_OPTIONS.md](docs/CLUSTER_OPTIONS.md) /
  [CLUSTER_RUNBOOK.md](docs/CLUSTER_RUNBOOK.md) — distributed runs
- [MULTIPHASE_DESIGN.md](docs/MULTIPHASE_DESIGN.md) /
  [WASM_BRIDGE_DESIGN.md](docs/WASM_BRIDGE_DESIGN.md) /
  [AGENT_MODE_DESIGN.md](docs/AGENT_MODE_DESIGN.md) /
  [COMPETITIVE_SPEC.md](docs/COMPETITIVE_SPEC.md) — subsystem designs
- `docs/qa/` — V&V / QA track: [VV_MASTER_PLAN.md](docs/qa/VV_MASTER_PLAN.md),
  VV_* process docs (claim-status vocabulary, traceability matrix, mutation
  plan, report template, visual-anomaly guide), dated audit snapshots, anomaly log
- `docs/archive/` — superseded snapshots (PM handoffs, dated one-shot docs)

## Documentation principles

- **One home per fact.** Each topic has one owning doc; others link to it.
  When facts change, update the owning doc in place — never fork a new file
  or append a contradicting section.
- **Index or archive.** Every file directly under `docs/` has a one-line
  entry in the Docs index above; update docs and index (CLAUDE.md AND
  AGENTS.md) in the same commit. Superseded or one-shot docs (handoffs,
  dated audits, drafts) move to `docs/archive/` with a pointer.
- **Sections over files.** Only a new subsystem or track justifies a new
  top-level doc; findings, reports, and experiment records go into the
  owning doc or the track's subdirectory (`docs/qa/`, `docs/proposals/`, …).
- **Lifecycle in the header**: living (updated in place), frozen spec
  (change = recorded revision), or snapshot (dated, never edited).
- **Docs move with the change.** A behavior-altering change updates the
  affected docs (physics + rationale + experiment results → PHYSICS.md;
  claims → LIMITATIONS.md; plan → PLAN.md). Delegating the doc work is fine
  if the report says so; neither written nor delegated = not done.
- **English, absolute dates, repo-relative links.** Write 2026-07-07, never
  "today".

## Team & conventions

- Fable is PM. **Implementation workhorse = codex CLI, fanned out in parallel
  (user directive 2026-07-05: "use codex fully, even ~100 parallel, finish
  fastest").** One order = one focused item bundle = one dedicated git
  worktree + background `codex exec -C <worktree>`; the PM merges landed
  branches in dependency order. Bundle same-file items into ONE order —
  parallelize across files, never within one. Concurrency is machine-bounded:
  ~5-8 compile-bound orders at once on this M5 Max (cargo is CPU/RAM heavy);
  more only for read/write-light work. Claude subagents handle design-heavy
  or cross-cutting work.
- **Validation tests are written adversarially** (by codex or Opus/Sonnet,
  from the spec in VALIDATION.md); a test order and an implementation order
  never share a worktree.
- codex invocation:
  `codex exec --sandbox workspace-write --skip-git-repo-check "<task>" < /dev/null`
  (model gpt-5.5). **`< /dev/null` is mandatory** — with a pipe as stdin,
  codex waits for EOF and hangs forever. Monitor progress via
  `~/.codex/sessions/<date>/rollout-*.jsonl`.
- **ALL artifacts in English (user directive 2026-07-05)**: code, identifiers,
  commits, docs, and user-facing strings (GUI / CLI / errors). Legacy Japanese
  content is being translated by a dedicated session; write new content in
  English only. Conversation with the user may remain Japanese — the product
  is English.
- Commit at each phase completion; never commit red tests (explicit WIP,
  flagged in the commit message, is the only exception).

## Working discipline (every agent on this repo)

- **Physical rigor is the prime directive; ad-hoc physics is BANNED (user
  directive 2026-07-06).** Every physical behavior anywhere in the stack —
  core, scenario, CLI, examples, demos, GUI — is either resolved from the
  governing equations or a literature-backed closure with a recorded
  derivation, validity domain, and its own validation test (PHYSICS.md entry
  mandatory). Banned outright: constants calibrated to pass an acceptance
  band; branches keyed to sample/case identity; clamps or caps that silently
  absorb transport; decorative physics terms. If a gate cannot be met without
  such a hack, STOP and report — the spec gets revised, not the physics
  faked. Existing violations are inventoried and hedged by the V&V track,
  not grandfathered. The executable procedure (provenance decision table,
  ban-list greps, two-layer gates, stop-rule template, behavior-review
  record) is `.claude/skills/lbmflow-physics-discipline`; every
  physics-affecting codex order embeds its clauses (lbmflow-codex-dispatch
  Step 1.5).
- **Behavior-validity review (user directive 2026-07-06).** After every
  experiment/demo run, before reporting: judge whether the OBSERVED behavior
  — spatial patterns, trends, signs, not just gated metrics — is physically
  plausible; name the dominant mechanism; separate resolved physics from
  closures and boundary artifacts; record the review in PHYSICS.md or the
  track's findings file. A metric passing its band does NOT validate a
  pattern no band covers (origin: D-track T18 edge-ring finding, 2026-07-06).
- **Finding routing (V&V loop).** Core-engine defect → phenomenon report +
  data package (scenario JSON, exported fields, metrics, repro command) to
  the core-engine session. Demo/example defect → PM dispatches a codex order.
  Spec defect → spec revision, rationale in PHYSICS.md.
- **Evidence-based progress.** Match every "done" claim against a tool result
  from THIS session (test output, file diff, run log); "done" means the
  relevant Build & test gate passed HERE — a codex order finishing is not
  evidence its branch is green. Report unverified as unverified, skipped as
  skipped, failures with their output. Fabricated progress is the worst
  possible failure.
- **Minimal scope.** Implement exactly what the order asks: no drive-by
  refactors, no helpers for one-off operations, no abstractions for
  hypothetical futures, no defensive code for impossible states (validate at
  system boundaries only: scenario JSON, CLI args, user input).
- **Finish, don't announce.** Never end a turn on a plan or "I'll now do X" —
  do X, then end. Stop only when the task is complete or blocked on input
  only the user can provide. When unsure between options, pick one and
  proceed; don't re-litigate settled decisions.

## Core design invariants (breaking these kills the validation suite)

- The single core is `crates/lbm-core` (V2; V1 retired 2026-07-05, freeze
  values in the `tests/v1_match.rs` header in branch history). The legacy V1
  API lives on as `lbm_core::compat`, used by the scenario / CLI / wasm 2D
  paths.
- D2Q9 direction ordering (single source of truth: the `Lattice` impl in
  lattice.rs): 0:(0,0), 1:(1,0), 2:(0,1), 3:(-1,0), 4:(0,-1), 5:(1,1),
  6:(-1,1), 7:(-1,-1), 8:(1,-1).
- f layout is q-major SoA with halo padding (fields.rs): `f[q*plane + cell]`,
  cell = z·(nx·ny) + y·nx + x — identical to the GPU coalescing assumption;
  never exposed in the public API.
- One step = collide → halo exchange → streaming → open-boundary BCs →
  boundary moments correction (CpuSimd fuses collide+stream+moments in
  step_band). Any change to pass structure or storage order must pass the
  gates in `tests/backend_simd_equiv.rs` and T13 (partition invariance)
  before it lands.
- Wall edges are a 1-cell solid rim; wall surfaces are half-way (midpoint
  between rim center and fluid center).
- Velocity moments include the Guo forcing F/2 correction — `sim.ux()` etc.
  return physical velocity.
- tau = 3·nu + 0.5 (cs² = 1/3).
