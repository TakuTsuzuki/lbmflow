# V&V Evidence Report Template

Lifecycle: living — report template, updated in place.

Use this template for adversarial validation campaign reports and per-branch
evidence packs. Do not mark a claim as done unless the evidence is from the
current session or from a linked artifact whose provenance is recorded here.

## Report Metadata

| Field | Value |
|---|---|
| Report ID |  |
| Date/time |  |
| Coordinator |  |
| Repository |  |
| Branch |  |
| Commit SHA |  |
| Base commit SHA |  |
| Worktree path |  |
| OS |  |
| Hardware |  |
| Rust toolchain |  |
| Node/npm toolchain |  |
| MPI toolchain |  |
| GPU adapter/driver |  |

## Scope

| Item | Value |
|---|---|
| Suborder / campaign area |  |
| Files intentionally touched |  |
| Files intentionally not touched |  |
| Physics-affecting change? | yes / no |
| New constants, closures, branches, clamps, or fallbacks? | yes / no |
| Required validation tier | docs-only / core / invariant / full / wasm / web / CLI / GPU / MPI |

## Command Log

Record every command used as evidence. If a command was skipped, put it in
`Tests Skipped` or `BENCH-PENDING Items`, not here.

| # | Command | Working directory | Exit | Evidence / notes |
|---:|---|---|---:|---|
| 1 |  |  |  |  |

## Tests Run

| Gate | Command | Exit | Result summary | Log/artifact path |
|---|---|---:|---|---|
| G1 build | `cargo build --workspace --release` |  |  |  |
| G2 tests | `cargo test --workspace --release` |  |  |  |
| G3 SIMD/backend invariants | `cargo test --release -p lbm-core --test backend_simd_equiv` |  |  |  |
| G3 T13 partition invariance | `cargo test --release -p lbm-core t13` |  |  |  |
| G4 full validation | `cargo test --release -- --include-ignored` |  |  |  |
| G5 WASM build | `wasm-pack build crates/lbm-wasm --target web --release --out-dir ../../web/src/engine/pkg` |  |  |  |
| G7 web build | `(cd web && npm run build)` |  |  |  |
| G8 CLI smoke | `./target/release/lbm presets run cavity` |  |  |  |
| GPU build/runtime |  |  |  |  |
| MPI build/runtime |  |  |  |  |

## Tests Skipped

Every skipped gate needs a reason. Use `BENCH-PENDING` for hardware- or
environment-dependent validation that cannot be run in this worktree.

| Gate / test | Reason skipped | Risk accepted | Follow-up owner |
|---|---|---|---|
|  |  |  |  |

## Artifacts Generated

Every experiment that produces fields or spatial behavior must include a visual
artifact path. Scalar-only physical validation is not reportable.

| Artifact | Type | Producer command | Purpose |
|---|---|---|---|
|  | manifest / CSV / PNG / VTK / log / report |  |  |

## Scalar Metrics

| Metric | Scenario / test | Measured value | Acceptance band | Status | Evidence |
|---|---|---:|---|---|---|
|  |  |  |  |  |  |

## Behavior-Review Summary

Complete this for every experiment/demo run before reporting the result. If no
experiment/demo was run, write `not applicable - docs-only` and cite the command
or diff evidence instead.

| Run ID | Visual artifact path | Dominant mechanism | Resolved physics vs closure | Boundary/seam/outlet sweep | Verdict | Routing |
|---|---|---|---|---|---|---|
|  |  |  |  |  | PHYSICAL / CLOSURE-DRIVEN / ARTIFACT / UNKNOWN | none / core / codex / spec |

## Claim Status

Allowed statuses: `VALIDATED`, `VERIFIED-ONLY`, `SPEC-ONLY`, `MISSING`,
`BENCH-PENDING`, `UNSAFE-CLAIM`, `STOP-RULE`.

| Claim area | Status | Evidence | Missing evidence | Downgrade recommendation |
|---|---|---|---|---|
| D2Q9 fluid core |  |  |  |  |
| D3Q19 fluid core |  |  |  |  |
| D3Q27 |  |  |  |  |
| GPU/wgpu |  |  |  |  |
| MPI/partition |  |  |  |  |
| f32/f16 precision |  |  |  |  |
| multiphase |  |  |  |  |
| particle/deposition |  |  |  |  |
| FSI/IBM/rotating boundary |  |  |  |  |
| scenario/CLI/MCP |  |  |  |  |

## Known Limitations

| Limitation | Scope affected | User-visible impact | Mitigation / next gate |
|---|---|---|---|
|  |  |  |  |

## Claims Allowed

List only claims supported by evidence in this report.

- 

## Claims Forbidden

List claims that must not be made from this evidence pack. Include claims that
would require GPU, MPI, long-horizon, visual, FSI, or cluster evidence not
present here.

- 

## BENCH-PENDING Items

| Item | Why pending | Required environment | Required command / artifact |
|---|---|---|---|
| GPU runtime |  |  |  |
| MPI runtime |  |  |  |
| Cluster scale |  |  |  |
| Full heavy validation |  |  |  |
| Long FSI benchmarks |  |  |  |

## STOP-RULE Reports

Use this exact shape when a gate cannot be met without unphysical terms.

```text
STOP-RULE: gate <T-id / band> is unreachable without unphysical terms.
Attempted: <what physical approaches were tried>
Blocking mechanism: <one-sentence physics of why it cannot pass>
Options for the PM: (a) spec/band revision, (b) resolved-physics capability
needed in core: <what>, (c) validated closure exists in literature: <ref>.
```

| Gate | Report | Routing |
|---|---|---|
|  |  |  |

## Findings

| Severity | Finding | Evidence | Routing |
|---|---|---|---|
| S0 correctness / false assurance |  |  |  |
| S1 high-risk physics |  |  |  |
| S2 validation gap |  |  |  |
| S3 docs/process |  |  |  |

## Merge Recommendation

| Branch / worktree | Recommendation | Conditions before merge |
|---|---|---|
|  | MERGE / DO NOT MERGE / BENCH-PENDING / COMMIT-READY-BUT-NOT-COMMITTED |  |

