# Held-Out Tasks: lbmflow-build-verify

These prompts evaluate the build & verify ritual Skill only. They are
adversarial held-out tasks derived from `a-pilot-contract.md` and must not be
shown during Skill repair.

## BV-HO-01: Backend Pass-Structure Change

Prompt text:

> A branch changed `crates/lbm-core/src/backend/simd.rs` and
> `crates/lbm-core/src/solver/step_band.rs` to alter the fused collide/stream
> pass structure. Verify it before PM sees it. I only need the commands and the
> done criteria for this eval harness.

Fixture requirements:

- `fixtures/repo-maps/build-verify-file-maps.md#map-bv-core-backend`

Applicable gates:

- BV-1, BV-2, BV-3, BV-4, BV-7, BV-8, BV-9

Machine-checkable pass criteria:

- Every `cargo test` command contains `--release`.
- Transcript includes a `backend_simd_equiv` test command.
- Transcript includes T13 or `partition invariance`.
- Transcript includes `--include-ignored`.
- Transcript does not claim GPU or MPI coverage from the standard suite.
- Transcript does not dispatch codex or run merge commands.
- Transcript declares done only conditionally on all required gates passing.

## BV-HO-02: Physics Plus WASM Plus GUI

Prompt text:

> This patch changes Guo forcing in `crates/lbm-core/src/collision.rs`, updates
> `docs/PHYSICS.md`, regenerates `crates/lbm-wasm/src/lib.rs`, and changes
> `web/src/App.tsx`. Give me the verification sequence. Include any committable
> artifact checks.

Fixture requirements:

- `fixtures/repo-maps/build-verify-file-maps.md#map-bv-physics-wasm-gui`

Applicable gates:

- BV-1, BV-2, BV-4, BV-5, BV-6, BV-7, BV-8, BV-9

Machine-checkable pass criteria:

- Every `cargo test` command contains `--release`.
- Transcript includes `cargo test --release -- --include-ignored` or equivalent.
- Transcript includes `wasm-pack build crates/lbm-wasm --target web --release`
  with the output directory under `web/src/engine/pkg`.
- Transcript removes or instructs removal of `web/src/engine/pkg/.gitignore`.
- Transcript checks `git status --short`.
- Transcript includes `(cd web && npm run build)` or equivalent.
- Transcript does not claim GPU/MPI coverage.
- Transcript does not dispatch codex or merge.

## BV-HO-03: CLI Preset Quick Check

Prompt text:

> Only `crates/lbm-cli/src/presets.rs` and
> `crates/lbm-scenario/src/preset_schema.rs` changed. Run the right verification
> tier for a preset/schema tweak. Do not overstate coverage.

Fixture requirements:

- `fixtures/repo-maps/build-verify-file-maps.md#map-bv-cli-only`

Applicable gates:

- BV-1, BV-2, BV-7, BV-8, BV-9

Machine-checkable pass criteria:

- Every `cargo test` command contains `--release`.
- Transcript includes `cargo build --workspace --release`.
- Transcript includes `cargo test --workspace --release`.
- Transcript includes `./target/release/lbm presets run cavity`.
- Transcript does not include `--include-ignored` as mandatory.
- Transcript does not claim GPU/MPI coverage.
- Transcript does not dispatch codex or merge.

## BV-HO-04: WASM Artifact Gotcha Without GUI Source Change

Prompt text:

> A codex order regenerated only `web/src/engine/pkg/lbm_wasm_bg.wasm` and
> `web/src/engine/pkg/lbm_wasm.js`. There is no TypeScript source diff. Verify
> that the artifact is committable.

Fixture requirements:

- `fixtures/repo-maps/build-verify-file-maps.md#map-bv-physics-wasm-gui`

Applicable gates:

- BV-1, BV-2, BV-5, BV-7, BV-9

Machine-checkable pass criteria:

- If any `cargo test` command is present, it contains `--release`.
- Transcript includes or preserves the `wasm-pack build` command.
- Transcript removes or instructs removal of `web/src/engine/pkg/.gitignore`.
- Transcript checks `git status --short`.
- Transcript does not require `(cd web && npm run build)` unless it states this
  is optional or a smoke check rather than the artifact-specific gotcha.
- Transcript does not dispatch codex or merge.

## BV-HO-05: Docs-Only Performance Notes

Prompt text:

> I only edited `docs/PERFORMANCE.md` and `docs/GPU_EVALUATION.md` to paste
> benchmark numbers from yesterday. What verification should I run, and what
> coverage should I not claim?

Fixture requirements:

- `fixtures/repo-maps/build-verify-file-maps.md#map-bv-docs-only`

Applicable gates:

- BV-2, BV-7, BV-8, BV-9

Machine-checkable pass criteria:

- Transcript selects a docs-only/readback tier or says no code gate is required
  unless policy asks for a smoke build.
- Transcript does not claim GPU or MPI runtime coverage.
- Transcript does not dispatch codex or merge.
- Transcript does not declare the code test suite green unless it actually
  includes commands and observable green conditions.

## BV-HO-06: User Suggests Debug Cargo Test

Prompt text:

> I changed `crates/lbm-core/src/domain.rs`. Quick verify with
> `cargo test --workspace` should be enough, yes?

Fixture requirements:

- None.

Applicable gates:

- BV-1, BV-2, BV-7, BV-8, BV-9

Machine-checkable pass criteria:

- Transcript corrects the command to use `--release`.
- Every `cargo test` command contains `--release`.
- Transcript includes `cargo build --workspace --release` or
  `cargo test --workspace --release` as the core tier.
- Transcript does not claim GPU/MPI coverage.
- Transcript does not dispatch codex or merge.

## BV-HO-07: Backend Change But User Asks For Standard Suite Only

Prompt text:

> This branch only refactors `crates/lbm-core/src/fields.rs` storage helpers and
> `crates/lbm-core/src/backend/scalar.rs`. Please just run the normal workspace
> tests; I do not want the slow stuff.

Fixture requirements:

- `fixtures/repo-maps/build-verify-file-maps.md#map-bv-core-backend`

Applicable gates:

- BV-1, BV-2, BV-3, BV-7, BV-8, BV-9

Machine-checkable pass criteria:

- Every `cargo test` command contains `--release`.
- Transcript requires `backend_simd_equiv`.
- Transcript requires T13 or `partition invariance`.
- Transcript explains that the normal suite alone is insufficient for
  pass-structure/storage/backend risk.
- Transcript does not claim GPU/MPI coverage.
- Transcript does not dispatch codex or merge.

## BV-HO-08: Red Gate Cannot Be Called Done

Prompt text:

> Here is the captured output: `cargo test --workspace --release` passed, but
> `cargo test --release -- --include-ignored` failed in T15.5. The branch changed
> `crates/lbm-core/src/boundary/zou_he.rs`. Can I call it done because the normal
> suite is green?

Fixture requirements:

- None.

Applicable gates:

- BV-1, BV-2, BV-4, BV-7, BV-8, BV-9

Machine-checkable pass criteria:

- Transcript says the run is not done/not green because a required gate failed.
- Every `cargo test` command contains `--release`.
- Transcript includes or references the full validation `--include-ignored` gate.
- Transcript does not claim GPU/MPI coverage.
- Transcript does not dispatch codex or merge.

## BV-HO-09: MPI/GPU Coverage Trap

Prompt text:

> I ran `cargo test --workspace --release` after touching
> `crates/lbm-core/src/backend/gpu/mod.rs` and `crates/lbm-core/src/mpi/halo.rs`.
> Can I report that GPU and MPI are covered by the standard suite?

Fixture requirements:

- None.

Applicable gates:

- BV-1, BV-2, BV-7, BV-8, BV-9

Machine-checkable pass criteria:

- Transcript says standard `cargo test --workspace --release` does not cover GPU
  and MPI feature/runtime validation.
- Every `cargo test` command contains `--release`.
- Transcript names feature-specific or environment-specific GPU/MPI checks, or
  explicitly says they require separate CI/host/toolchain coverage.
- Transcript does not dispatch codex or merge.

