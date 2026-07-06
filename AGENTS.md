# LBMFlow — Lattice Boltzmann Method Fluid Simulator

Agent instructions for this repository (read by codex and any non-Claude coding
agent; Claude agents receive the same invariants via CLAUDE.md — keep the two
files in sync when editing either).

Commercial-grade LBM simulator. Rust core + TypeScript GUI + Agent mode.
**Required reading**: [docs/PLAN.md](docs/PLAN.md) (phase plan / team structure),
[docs/VALIDATION.md](docs/VALIDATION.md) (validation test specs = acceptance criteria).

## Build & test

```bash
cargo build --workspace --release
cargo test --workspace --release          # normal suite (LBM is ~50x slower in debug — ALWAYS use --release)
cargo test --release -- --include-ignored # full validation incl. heavy benches (~5 min)
# WASM (for the web GUI; lbm-wasm is outside the workspace):
wasm-pack build crates/lbm-wasm --target web --release --out-dir ../../web/src/engine/pkg
#   (after generating, delete pkg/.gitignore and commit pkg)
cd web && npm run build                   # GUI (tsc strict + vite)
./target/release/lbm presets run cavity   # CLI smoke test
```

## Repository map

- `crates/lbm-core` — the single core (V2): D2Q9/D3Q19 lattices, CPU
  scalar/SIMD backends, wgpu GPU backend (feature `gpu`, off by default —
  `cargo test --workspace` does NOT cover it; CI runs `--features gpu` on GPU
  hosts), MPI halo exchange (feature `mpi`, off by default, needs a native MPI
  toolchain — see scripts/test_mpi.sh), legacy 2D facade in `compat/`
- `crates/lbm-scenario` — JSON scenario schema + runner (2D compat path)
- `crates/lbm-cli` — `lbm` binary: presets (list/show/run), gallery, schema,
  scenario run, MCP server (7 tools incl. async start_run/run_status/list_runs)
- `crates/lbm-wasm` — WASM bindings for the web GUI (excluded from the workspace)
- `crates/lbm-gpu-proto` — standalone wgpu evaluation prototype (excluded from
  the workspace; measurement record, superseded by the in-core `gpu` module)
- `web/` — TypeScript GUI (Vite; engine WASM committed under `web/src/engine/pkg`)

## Docs index (read on demand)

- [PLAN.md](docs/PLAN.md) — milestones M-A…M-F, current queue ·
  [VALIDATION.md](docs/VALIDATION.md) — acceptance criteria (T1…T15.x)
- [PHYSICS.md](docs/PHYSICS.md) — physics decisions + experiment log
  (update whenever you change physics)
- [ARCHITECTURE_V2.md](docs/ARCHITECTURE_V2.md) — dimension × lattice ×
  precision × backend × partition design
- [SOLVER_IMPROVEMENT_SPEC.md](docs/SOLVER_IMPROVEMENT_SPEC.md) — R-Phase spec
- [REQ_STIRRED_REACTOR.md](docs/REQ_STIRRED_REACTOR.md) — M-F requirements ·
  [T15_5_CAVITY3D_REFERENCE.md](docs/T15_5_CAVITY3D_REFERENCE.md) — 3D cavity reference data
- [DISPERSED_DEPOSITION.md](docs/DISPERSED_DEPOSITION.md) — D-track dispersed-phase
  deposition tool (frozen spec, P0–P4, CR-1/2/3, acceptance = T18)
- [PERFORMANCE.md](docs/PERFORMANCE.md) / [GPU_EVALUATION.md](docs/GPU_EVALUATION.md) /
  [BENCH_COMPARISON_DRAFT.md](docs/BENCH_COMPARISON_DRAFT.md) — perf measurements
- [MPI_GUIDE.md](docs/MPI_GUIDE.md) / [CLUSTER_OPTIONS.md](docs/CLUSTER_OPTIONS.md) — distributed runs
- [MULTIPHASE_DESIGN.md](docs/MULTIPHASE_DESIGN.md) /
  [WASM_BRIDGE_DESIGN.md](docs/WASM_BRIDGE_DESIGN.md) /
  [AGENT_MODE_DESIGN.md](docs/AGENT_MODE_DESIGN.md) /
  [COMPETITIVE_SPEC.md](docs/COMPETITIVE_SPEC.md) — subsystem designs

## Team & conventions

- Fable is PM. **Implementation workhorse = codex CLI, fanned out aggressively
  (user directive 2026-07-05: "use codex fully, even ~100 parallel, finish
  fastest")**: one order = one focused item bundle = one dedicated git worktree
  + background `codex exec -C <worktree>`; PM merges landed branches in
  dependency order. Bundle same-file items into ONE order (parallelize across
  files, not within a file). Effective concurrency is machine-bounded — cargo
  is CPU/RAM heavy, ~5-8 compile-bound orders at once on this M5 Max; add more
  only for read/write-light work. Claude subagents remain for design-heavy or
  cross-cutting work. **Validation tests are written adversarially by codex or
  Opus/Sonnet from the spec (VALIDATION.md)**, kept separate from the
  implementation (a test order and an implementation order never share a
  worktree).
- **Language policy (user directive 2026-07-05): ALL artifacts in English** —
  code, identifiers, commit messages, documentation, and user-facing strings
  (docs / GUI / CLI / error messages). Legacy Japanese content is being
  translated by a dedicated session; write new content in English only.
- When you change the physics spec, record the rationale and experimental
  results in docs/PHYSICS.md.
- Commit at each phase completion. Never commit with red tests (WIP is the
  exception; say so explicitly in the commit message).

## Working discipline (applies to every agent on this repo)

- **Physical rigor is the prime directive; ad-hoc physics is BANNED (user
  directive 2026-07-06).** Every physical behavior anywhere in the stack —
  core, scenario, CLI, examples, demos, GUI — must be either resolved from
  the governing equations or a literature-backed closure with a recorded
  derivation, validity domain, and its own validation test (PHYSICS.md entry
  mandatory). Prohibited outright: constants calibrated to pass an acceptance
  band; branches keyed to sample/case identity (e.g. "harshness" switches);
  position clamps or caps that silently absorb transport; decorative physics
  terms. If a gate cannot be met without such a hack, STOP and report — the
  spec gets revised, not the physics faked. Existing violations are being
  inventoried and risk-hedged with minimal effort (V&V sweep commissioned
  2026-07-06), not grandfathered. **The executable procedure (provenance
  decision table, ban-list greps, two-layer gates, stop-rule template,
  behavior-review record) is `.claude/skills/lbmflow-physics-discipline` —
  every developer agent follows it mechanically; every physics-affecting
  codex order embeds its clauses (see lbmflow-codex-dispatch Step 1.5).**
- **V&V loop and finding routing.** V&V continuously runs experiment
  matrices, visualizes the results (qa-viewer), and behavior-reviews them;
  every anomaly becomes a finding. Routing: core-engine defect → send the
  phenomenon report + data package (scenario JSON, exported fields, metrics,
  repro command) to the core-engine session and request the fix;
  demo/example defect → the PM dispatches a codex order; spec defect → spec
  revision with the rationale recorded in PHYSICS.md.
- **Behavior-validity review (user directive 2026-07-06).** After every
  experiment/demo run, before reporting results: review whether the OBSERVED
  behavior — spatial patterns, trends, signs, not just the gated metrics — is
  physically plausible. Identify the dominant mechanism, separate resolved
  physics from model closures / boundary artifacts (clamps, ad-hoc branches,
  calibrated constants), and record the review in PHYSICS.md or the track's
  findings file. A metric passing its band does NOT validate a pattern no
  band covers. (Origin: the dispersed-deposition gentle case deposited in an
  edge ring despite center dispensing — mechanism plausible, but the
  magnitude was set by an uncalibrated closure branch + a side-wall position
  clamp, and no gate looked at the spatial pattern.)
- **Evidence-based progress.** Before reporting anything as done, match each
  claim against a tool result from THIS session (test output, file diff, run
  log). Report unverified work as unverified; report skipped steps as skipped;
  report failing tests with their output. "Done" means the relevant gate in
  Build & test actually passed here — a codex order finishing is not evidence
  its branch is green. Fabricated progress is the worst possible failure.
- **Minimal scope.** Implement exactly what the order/task asks: no drive-by
  refactors, no helpers for one-off operations, no abstractions for
  hypothetical future requirements, no defensive code for impossible states
  (validate at system boundaries only — scenario JSON, CLI args, user input).
  A bug fix does not need surrounding cleanup.
- **Finish, don't announce.** Never end a turn on a plan, a checklist, or
  "I'll now do X" — do X first, then end. Stopping is allowed only when the
  task is complete or blocked on input only the user can provide. When unsure
  between options, pick one recommendation and proceed; don't enumerate
  alternatives you won't take or re-litigate decisions already made.

## Core design invariants (breaking these kills the whole validation suite)

- The single core is `crates/lbm-core` (formerly lbm-core2 = V2 architecture.
  V1 retired 2026-07-05; equivalence freeze values are in the `tests/v1_match.rs`
  header in branch history). The legacy V1 API is provided by `lbm_core::compat`
  (public facade), used by the scenario / CLI / wasm 2D paths.
- The D2Q9 direction ordering defined in lattice.rs (the `Lattice` trait impl)
  is the single source of truth: 0:(0,0), 1:(1,0), 2:(0,1), 3:(-1,0), 4:(0,-1),
  5:(1,1), 6:(-1,1), 7:(-1,-1), 8:(1,-1).
- f layout is q-major SoA (fields.rs, with halo padding): `f[q*plane + cell]`,
  cell = z·(nx·ny) + y·nx + x. Identical to the GPU coalescing assumption.
  Never exposed in the public API.
- One step = collide → halo exchange → streaming → open-boundary BCs → boundary
  moments correction (CpuSimd fuses collide+stream+moments in step_band).
  Any change to the pass structure or storage order must pass the bit/threshold
  gates in `tests/backend_simd_equiv.rs` and T13 (partition invariance) before
  it lands.
- Wall edges are a 1-cell solid rim. Wall surfaces are half-way (midpoint
  between rim center and fluid center).
- Velocity moments include the Guo forcing F/2 correction (`sim.ux()` etc.
  return physical velocity).
- tau = 3*nu + 0.5 (cs² = 1/3).
