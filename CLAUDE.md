# LBMFlow — Lattice Boltzmann Method Fluid Simulator

An LBM simulator aiming for commercial grade. Rust core + TypeScript GUI + Agent mode.
**Required reading**: [docs/PLAN.md](docs/PLAN.md) (phase plan, team structure),
[docs/VALIDATION.md](docs/VALIDATION.md) (validation test spec = acceptance criteria).

## Build & test

```bash
cargo build --workspace --release
cargo test --workspace --release          # regular suite (LBM is ~50x slower in debug — always --release)
cargo test --release -- --include-ignored # full validation incl. heavy benchmarks (~5 min)
# WASM (for the web GUI; lbm-wasm is outside the workspace):
wasm-pack build crates/lbm-wasm --target web --release --out-dir ../../web/src/engine/pkg
#   (after generating, delete pkg/.gitignore and commit pkg — standing procedure)
cd web && npm run build                   # GUI (tsc strict + vite)
./target/release/lbm presets run cavity   # CLI smoke test
```

## Team structure & conventions

- Fable is PM. Implementation is delegated to Opus/Sonnet subagents / codex.
  **Validation tests are written adversarially by codex or Opus/Sonnet from the spec
  (VALIDATION.md)**, kept separate from the implementation.
- codex invocation example: `codex exec --sandbox workspace-write --skip-git-repo-check "<task>" < /dev/null`
  (model gpt-5.5. **`< /dev/null` is REQUIRED** — with a pipe as stdin it waits for
  EOF and hangs forever. Progress can be monitored via updates to
  `~/.codex/sessions/<date>/rollout-*.jsonl`)
- **Language policy (user directive 2026-07-05): ALL artifacts in English** — code,
  identifiers, commit messages, documentation, and user-facing strings (docs / GUI /
  CLI / error messages). Legacy Japanese content is being translated by a dedicated
  session; write new content in English only. (Conversation with the user may remain
  Japanese — the user is a Japanese speaker; the product is English.)
- Whenever a physics spec changes, record the reason and experimental results in
  docs/PHYSICS.md.
- git commit at the end of each phase. Never commit with tests red (WIP is the
  exception — state it explicitly in the message).

## Core design commitments (breaking these wipes out validation)

- The single core is `crates/lbm-core` (formerly lbm-core2 = the V2 architecture.
  V1 retired 2026-07-05; the frozen equivalence values live in the
  `tests/v1_match.rs` header in branch history). The old V1 API is provided by
  `lbm_core::compat` (public facade), used by the 2D paths of scenario / CLI / wasm.
- The D2Q9 direction ordering defined in lattice.rs (the `Lattice` trait impl) is
  the single source of truth. 0:(0,0), 1:(1,0), 2:(0,1), 3:(-1,0), 4:(0,-1),
  5:(1,1), 6:(-1,1), 7:(-1,-1), 8:(1,-1).
- f layout is q-major SoA (fields.rs, halo-padded): `f[q*plane + cell]`,
  cell = z·(nx·ny) + y·nx + x. Identical to the GPU-coalescing assumption. Never
  exposed in the public API.
- 1 step = collision → halo exchange → streaming → open-boundary BCs → boundary-line
  moments correction (CpuSimd fuses collide+stream+moments in step_band). Any change
  that alters the pass structure or storage order must pass the bit/threshold gates
  of `tests/backend_simd_equiv.rs` and T13 (partition invariance) before landing.
- Wall edges are a 1-cell solid rim. The wall surface is half-way (midpoint between
  rim center and fluid center).
- Velocity moments include the Guo forcing F/2 correction (`sim.ux()` etc. are
  physical velocities).
- tau = 3*nu + 0.5 (cs² = 1/3).
