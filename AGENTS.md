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
  scalar/SIMD backends, wgpu GPU backend (feature `gpu`), MPI halo exchange
  (feature `mpi`), legacy 2D facade in `compat/`
- `crates/lbm-scenario` — JSON scenario schema + runner (2D compat path)
- `crates/lbm-cli` — `lbm` binary: presets, gallery, schema, scenario run, MCP server
- `crates/lbm-wasm` — WASM bindings for the web GUI (outside the workspace)
- `crates/lbm-gpu-proto` — wgpu evaluation prototype (measurement record;
  superseded by the in-core `gpu` module)
- `web/` — TypeScript GUI (Vite; engine WASM committed under `web/src/engine/pkg`)

## Docs index (read on demand)

- [PLAN.md](docs/PLAN.md) — milestones M-A…M-F, current queue ·
  [VALIDATION.md](docs/VALIDATION.md) — acceptance criteria (T1…T15.x)
- [PHYSICS.md](docs/PHYSICS.md) — physics decisions + experiment log
  (update whenever you change physics)
- [ARCHITECTURE_V2.md](docs/ARCHITECTURE_V2.md) — dimension × lattice ×
  precision × backend × partition design
- [SOLVER_IMPROVEMENT_SPEC.md](docs/SOLVER_IMPROVEMENT_SPEC.md) — R-Phase spec ·
  [REVIEW_2026-07-05.md](docs/REVIEW_2026-07-05.md) (+ `_2`) — solver review findings
- [REQ_STIRRED_REACTOR.md](docs/REQ_STIRRED_REACTOR.md) — M-F requirements ·
  [T15_5_CAVITY3D_REFERENCE.md](docs/T15_5_CAVITY3D_REFERENCE.md) — 3D cavity reference data
- [PERFORMANCE.md](docs/PERFORMANCE.md) / [GPU_EVALUATION.md](docs/GPU_EVALUATION.md) /
  [BENCH_COMPARISON_DRAFT.md](docs/BENCH_COMPARISON_DRAFT.md) — perf measurements
- [MPI_GUIDE.md](docs/MPI_GUIDE.md) / [HPC_SCALING.md](docs/HPC_SCALING.md) /
  [CLUSTER_OPTIONS.md](docs/CLUSTER_OPTIONS.md) — distributed runs
- [MULTIPHASE_DESIGN.md](docs/MULTIPHASE_DESIGN.md) /
  [WASM_BRIDGE_DESIGN.md](docs/WASM_BRIDGE_DESIGN.md) /
  [AGENT_MODE_DESIGN.md](docs/AGENT_MODE_DESIGN.md) /
  [COMPETITIVE_SPEC.md](docs/COMPETITIVE_SPEC.md) — subsystem designs

## Team & conventions

- Fable is PM. Implementation is delegated to Opus/Sonnet subagents / codex.
  **Validation tests are written adversarially by codex or Opus/Sonnet from the
  spec (VALIDATION.md)**, kept separate from the implementation.
- **Language policy (user directive 2026-07-05): ALL artifacts in English** —
  code, identifiers, commit messages, documentation, and user-facing strings
  (docs / GUI / CLI / error messages). Legacy Japanese content is being
  translated by a dedicated session; write new content in English only.
- When you change the physics spec, record the rationale and experimental
  results in docs/PHYSICS.md.
- Commit at each phase completion. Never commit with red tests (WIP is the
  exception; say so explicitly in the commit message).

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
