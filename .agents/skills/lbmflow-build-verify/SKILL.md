---
name: lbmflow-build-verify
description: >-
  Run the LBMFlow build-and-verify ritual — the exact cargo/wasm/web/CLI command
  sequence that gates every LBMFlow change before commit or merge. Use this
  whenever you are about to commit, land, or hand off LBMFlow work, or when asked
  to "build", "run the tests", "verify", "check it's green", "run validation", or
  "make sure nothing broke" in this repo (Rust core + TypeScript GUI + WASM
  engine). Also use before reporting a codex/subagent order as done. Do NOT use
  this to dispatch parallel codex orders (that is lbmflow-codex-dispatch) or to
  author validation tests — this Skill only builds and checks an already-written
  tree.
---

# LBMFlow build & verify ritual

LBMFlow is a commercial-grade Lattice Boltzmann simulator: a Rust workspace
(`lbm-core`, `lbm-scenario`, `lbm-cli`), a WASM engine (`lbm-wasm`, **outside**
the workspace), and a TypeScript GUI (`web/`). A change is not "done" until the
gates below are green. This Skill carries the exact commands and the pass
criteria so you never improvise them.

**One hard rule that saves ~50x wall time:** LBM is ~50x slower in debug. **Every
`cargo test` MUST include `--release`.** A debug test run is not a valid gate —
it is either uselessly slow or (for heavy benches) will not finish.

## Step 0 — Scope the change, pick the gate tier

Run only the gates the change can affect, plus the always-on core gate. Match the
files you (or the order you are verifying) touched against this table:

| Files touched | Gates to run (in order) |
|---|---|
| Any Rust under `crates/` (except wasm-only) | Core: G1 build → G2 test |
| Backend pass structure / storage order (`fields.rs`, `collide`, `stream`, `step_band`, SIMD) | Core **+** G3 invariant gate (`backend_simd_equiv` + T13) **+** G4 full incl. `--include-ignored` |
| Physics / moments / BC / `tau` / forcing | Core + G4 full incl. `--include-ignored` |
| `crates/lbm-wasm/**` or anything the GUI engine imports | G5 wasm build → G6 pkg commit fix → G7 web build |
| `web/**` (TypeScript/Vite only) | G7 web build |
| `crates/lbm-cli/**`, presets, scenario schema | Core + G8 CLI smoke |
| Docs only (`docs/**`, `*.md`) | None — but confirm AGENTS.md/AGENTS.md stayed in sync (see Step 3) |

If unsure which tier applies, run one tier up. The cost of an extra `--release`
suite is minutes; the cost of a missed regression is the whole validation suite.

## The gate commands (copy-runnable, exact)

Run from the worktree root. These are the corrected forms — use them verbatim.

```bash
# G1 — core build
cargo build --workspace --release

# G2 — normal test suite (the default gate for any Rust change)
cargo test --workspace --release

# G3 — invariant gate (REQUIRED for any backend pass-structure / storage-order change)
cargo test --release -p lbm-core --test backend_simd_equiv
cargo test --release -p lbm-core t13          # partition-invariance (name-filter; runs T13 cases)

# G4 — full validation incl. heavy benches (~5 min; physics + pass-structure changes)
cargo test --release -- --include-ignored

# G5 — WASM engine build (GUI dependency; lbm-wasm is OUTSIDE the workspace)
wasm-pack build crates/lbm-wasm --target web --release --out-dir ../../web/src/engine/pkg

# G6 — pkg commit fix (MANDATORY after G5 — see WHY below)
if [ -f web/src/engine/pkg/.gitignore ]; then rm web/src/engine/pkg/.gitignore; fi
git -C web/src/engine/pkg status --short

# G7 — GUI build (tsc strict + vite)
(cd web && npm run build)

# G8 — CLI smoke
./target/release/lbm presets run cavity
```

**WHY G6 is mandatory and easy to forget:** `wasm-pack` emits a `.gitignore`
inside `pkg/` that would exclude the committed engine build. The committed
`web/src/engine/pkg` **is** the engine the GUI ships. If you skip G6, `git add`
silently drops the new engine and the GUI builds against a stale WASM. Deleting
that `.gitignore` and confirming `git status --short` shows the regenerated
`pkg` files is the observable proof the engine actually updated.

## Verification gate — the done check

The Skill's own done-criterion. A change passes verification only when ALL gates
required by its tier report the following observable results:

| Gate | Observable done-check (what green looks like) |
|---|---|
| G1 | Exit 0, `Finished \`release\` profile` line, no `error[` lines. |
| G2 | Every `test result:` line reads `ok`. **Zero `FAILED`.** No `error:` compile lines. |
| G3 | `backend_simd_equiv`: `test result: ok`. T13 cases: `ok` (bit/threshold gates hold). |
| G4 | Ran exactly once end-to-end; final `test result: ok`. Note the wall time (~5 min is normal). |
| G5 | Exit 0; `[INFO]: :-) Done` from wasm-pack; `web/src/engine/pkg/*.wasm` regenerated (check mtime). |
| G6 | `web/src/engine/pkg/.gitignore` absent; `git status --short` lists the regenerated `pkg` files as changed/new. |
| G7 | `vite build` exits 0; `web/dist/` produced; no `tsc` type errors. |
| G8 | `lbm presets run cavity` exits 0 and prints step/summary output (no panic, no `Error:`). |

If any required gate is not green, the change is **not done** — fix or report red
(never commit red; WIP is the only exception and must say "WIP" in the message).

## Worked example (end-to-end)

Task: a codex order modified `crates/lbm-core/src/backends/simd.rs` (`step_band`
fusion) and you must verify before merging its branch.

1. **Scope (Step 0):** backend pass-structure change → tier = Core + G3 + G4.
2. **G1:** `cargo build --workspace --release` → `Finished`. OK.
3. **G2:** `cargo test --workspace --release` → all `test result: ok`. OK.
4. **G3:** run both invariant commands. `backend_simd_equiv` → `ok`;
   T13 → `ok`. This is the critical gate — a fused-kernel change that drifts a
   bit fails here, not in G2.
5. **G4:** `cargo test --release -- --include-ignored` → runs ~5 min, ends
   `test result: ok`.
6. All required gates green → verification passes. Report: "Verified: G1/G2/G3
   (backend_simd_equiv + T13) / G4 all green; T13 partition-invariance holds."

## Top failure modes (and the fix)

- **Ran `cargo test` without `--release`.** Symptom: suite crawls, or heavy
  benches never finish. Fix: re-run with `--release`. Never gate on a debug run.
- **Skipped G6 after a WASM change.** Symptom: GUI builds green but ships a stale
  engine; `git add web/src/engine/pkg` adds nothing. Fix: run G6, confirm the
  `.gitignore` is gone and `git status --short` shows the new `pkg` files.
- **Backend change gated only by G2.** Symptom: `cargo test --workspace` is green
  but a bit-level drift slips through. Fix: G2 does not cover the bit/threshold
  invariants — you MUST also run G3 (`backend_simd_equiv` + T13) for any
  pass-structure or storage-order change.
- **Treated G4's long runtime as a hang.** ~5 min is expected. Do not kill it and
  do not re-launch a second copy in parallel — let the single run finish.
- **`lbm` binary not found in G8.** It is a release artifact: run G1 first, then
  `./target/release/lbm ...`.
- **GPU/MPI features assumed covered.** `cargo test --workspace` does NOT build
  `--features gpu` or `--features mpi`. GPU is CI-only on GPU hosts; MPI is
  verified separately via `scripts/test_mpi.sh` (needs a native arm64 MPI). Do
  not claim GPU/MPI coverage from the standard gates.

## Step 3 — pre-commit sync check (docs)

AGENTS.md (Codex agents) and AGENTS.md (codex/other agents) must stay in sync —
they mirror the same invariants. If your change touched build commands, the
repository map, or any core invariant, confirm both files reflect it before
committing. A drift between them silently gives codex and Codex different rules.
