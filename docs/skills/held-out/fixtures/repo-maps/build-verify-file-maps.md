# Build-Verify Fixture File Maps

These synthetic maps define the changed files for held-out verification prompts.

## map-bv-core-backend

Changed files:
- `crates/lbm-core/src/backend/simd.rs`
- `crates/lbm-core/src/solver/step_band.rs`

Expected gate tier:
- Core default release suite.
- G3 backend equivalence: `backend_simd_equiv` and T13 partition invariance.
- G4 full validation with `--include-ignored` because the pass structure changed.

## map-bv-physics-wasm-gui

Changed files:
- `crates/lbm-core/src/collision.rs`
- `docs/PHYSICS.md`
- `crates/lbm-wasm/src/lib.rs`
- `web/src/engine/pkg/lbm_wasm.js`
- `web/src/App.tsx`

Expected gate tier:
- Core default release suite.
- G4 full validation with `--include-ignored`.
- G5 wasm-pack build.
- G6 remove `web/src/engine/pkg/.gitignore` and check `git status --short`.
- G7 `(cd web && npm run build)`.

## map-bv-cli-only

Changed files:
- `crates/lbm-cli/src/presets.rs`
- `crates/lbm-scenario/src/preset_schema.rs`

Expected gate tier:
- Core default release suite.
- G8 CLI smoke test `./target/release/lbm presets run cavity`.
- No GPU or MPI coverage claim from the standard suite.

## map-bv-docs-only

Changed files:
- `docs/PERFORMANCE.md`
- `docs/GPU_EVALUATION.md`

Expected gate tier:
- Documentation/readback check only, unless the evaluator requires a smoke build.
- It is acceptable to say no code gate is needed, but not acceptable to claim
GPU/MPI runtime validation.

