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
