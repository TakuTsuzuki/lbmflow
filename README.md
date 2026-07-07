# LBMFlow

**A bioprocess-specific CFD core.** Rust engine, scenario CLI, and a native
MCP server so AI agents can drive it. Its purpose is to compute the QOIs
that drive stirred-tank cell-culture and bioreactor process decisions:
power number, mixing time, gas holdup, kLa, shear exposure, oxygen
exposure, cell / microcarrier damage risk, and scale-up operating windows.

Generic CFD parity is *not* a goal. See [docs/BIOPROCESS_PIVOT.md](docs/BIOPROCESS_PIVOT.md)
for the pivot rationale and what has been retracted / preserved.

Licensed under MIT OR Apache-2.0.

**Version:** `0.2.0-bioprocess.0`. The pre-pivot general-purpose LBM
snapshot is `git tag v1-lbm-general-final`; its docs live under
`docs/archive/2026-07-07-pivot/`.

## Product mission

LBMFlow is the CFD side of a bioprocess design workflow, not a
general-purpose simulator. Concretely:

- **Impeller and vessel hydrodynamics** — Np, P/V, discharge Nq for
  Rushton and pitched-blade in stirred tanks.
- **Shear-rate and viscous-stress fields** with percentile distributions
  (P50 / P90 / P95 / P99); wall-shear diagnostics.
- **Passive-scalar mixing** — CV-based t95 and t99 mixing time.
- **Resolved gas-liquid** via conservative Allen-Cahn phase field, with
  surface tension, contact angle, and a sparger conservation ledger.
- **Oxygen transport** with Henry-equilibrium interfacial flux; kLa QOI
  from dynamic gassing.
- **Cell / microcarrier trajectories** with shear-exposure integrals and
  Schiller-Naumann validity enforcement.
- **Point-bubble + PBM engineering mode** for aeration with d32 and
  interfacial-area-based kLa.
- **UQ / sweep runner** and **scale-up operating-window evaluator**.
- **Evidence-tier gate** requiring calibration + holdout + UQ +
  sensitivity records before any evidence claim can be labelled.

The plan lives in [docs/PLAN.md](docs/PLAN.md) as tickets BCFD-000..110
across milestones M0..M3.

## Design principles

- **Explicit accuracy–speed control.** Every trade-off — collision
  operator, precision, backend, resolution — is a first-class knob, not a
  hidden default. Pick your point on the curve; the trade-off is measured,
  not asserted.
- **Physically rigorous.** Every model term is derived from the governing
  equations or a literature-backed closure with a recorded derivation,
  validity domain, and its own validation test
  ([docs/PHYSICS.md](docs/PHYSICS.md)). Constants calibrated to pass a
  band, case-keyed branches, silent clamps that absorb transport, and
  decorative physics are prohibited — if a gate cannot be met without a
  hack, the spec is revised, not the physics.
- **Validated adversarially.** The bioprocess validation suite
  ([docs/VALIDATION_BIOPROCESS.md](docs/VALIDATION_BIOPROCESS.md),
  VB-01..VB-08) is authored independently of the engine from the public
  spec. The engine is fixed until the tests pass — never the other way
  around.
- **Nothing silent.** Unsupported combinations fail with structured
  errors carrying an `UnsupportedReason` (BCFD-002 capability registry).
  QOIs without units / method / averaging metadata fail serialisation.
  `max` alone is never a report — distributions are required.
- **Evidence claims are mechanical.** BCFD-091 gates evidence-tier labels
  on calibration/holdout separation, UQ intervals, and mesh / time-step
  sensitivity records.

## Getting started

### CLI (bioprocess workflow — target surface, ticket-gated)

```bash
cargo build --workspace --release
./target/release/lbm capabilities --json                # what is supported
./target/release/lbm schema --bioprocess                # bioprocess scenario JSON schema
./target/release/lbm bioprocess validate my-tank.json   # unit feasibility + capability check
./target/release/lbm bioprocess run my-tank.json        # run + QOI extraction
./target/release/lbm bioprocess qoi out/my-tank/        # print qoi.json + provenance
./target/release/lbm bioprocess report out/my-tank/     # human-readable report
./target/release/lbm bioprocess sweep my-sweep.json     # UQ / parameter sweep
./target/release/lbm bioprocess evidence-check out/     # evidence-gate result
```

The `bioprocess` subcommands land per BCFD-092. Until then the M0 subset
(single-phase Np / P/V / mixing / shear) is what the CLI can actually
produce.

### CLI (legacy demos, still runnable)

```bash
./target/release/lbm presets list             # cavity, cylinder-karman, two-phase-droplet, droplet-on-wall
./target/release/lbm presets run cavity       # → out/cavity/ (emits legacy-preset warning)
```

Legacy presets are not bioprocess decision-grade and emit a warning to
stderr on run.

### MCP server (AI-agent integration)

```bash
claude mcp add lbmflow -- /path/to/target/release/lbm mcp
```

Legacy tools (`run_scenario`, `start_run` / `run_status` / `list_runs`,
`validate_scenario`, `list_presets`, `get_schema`) plus the bioprocess
tool surface added per BCFD-092 (`validate_bioprocess_scenario`,
`run_bioprocess_scenario`, `get_bioprocess_qoi`,
`generate_bioprocess_report`, `check_evidence_gate`).

### Library

The Rust API is documented in [docs/ARCHITECTURE_V2.md](docs/ARCHITECTURE_V2.md).
The v2 core lives in `lbm_core::prelude`; the legacy 2D facade lives in
`lbm_core::compat::prelude`.

## Building and testing

```bash
cargo build --workspace --release
cargo test  --workspace --release --no-fail-fast   # default gate — always --release --no-fail-fast
cargo test  --release -- --include-ignored         # + heavy bioprocess validation (~5 min)

# WebAssembly for the browser GUI (lbm-wasm is outside the workspace):
wasm-pack build crates/lbm-wasm --target web --release --out-dir ../../web/src/engine/pkg
cd web && npm run build

# Optional features:
cargo test  --workspace --release --features gpu   # wgpu backend (GPU hosts only)
cargo test  --workspace --release --features mpi   # requires a native MPI toolchain
```

LBM is roughly 50× slower in debug; `--release` is not optional. The
default gate is `cargo test --workspace --release --no-fail-fast`.
Piping the gate through `tail` / `grep` eats the exit code.

## Repository map

- `crates/lbm-core` — the engine: D2Q9 / D3Q19 / D3Q27 lattices, CPU
  scalar / SIMD backends, wgpu GPU backend (feature `gpu`), MPI halo
  exchange (feature `mpi`), WALE LES, rotating IBM, Bouzidi curved walls.
  Legacy 2D facade in `compat/`. Bioprocess physics modules (phase field,
  materials, sparger, oxygen, cells, bubbles, PBM, damage) land per BCFD
  tickets.
- `crates/lbm-scenario` — JSON scenario schema (legacy `Scenario` +
  `BioprocessScenario`) and runner.
- `crates/lbm-cli` — `lbm` binary: presets (legacy demos), gallery,
  schema, scenario run, MCP server, and the bioprocess CLI surface (per
  BCFD-092).
- `crates/lbm-wasm` — WASM bindings for the web GUI (outside the
  workspace). Not on the bioprocess critical path until BCFD-081.
- `crates/lbm-gpu-proto` — wgpu evaluation prototype (measurement
  record, superseded by the in-core `gpu` module).
- `web/` — TypeScript / Vite GUI.

## Documentation

Bioprocess product docs (living):

- [docs/BIOPROCESS_PIVOT.md](docs/BIOPROCESS_PIVOT.md) — pivot, retracted
  vs preserved claims.
- [docs/SPEC_BIOPROCESS_CORE.md](docs/SPEC_BIOPROCESS_CORE.md) — intended
  use, tiers, QOI catalog, scenario schema.
- [docs/VALIDATION_BIOPROCESS.md](docs/VALIDATION_BIOPROCESS.md) —
  VB-01..VB-08.
- [docs/CREDIBILITY_BIOPROCESS.md](docs/CREDIBILITY_BIOPROCESS.md) —
  calibration / holdout / UQ policy.
- [docs/MODEL_RISK_MATRIX.md](docs/MODEL_RISK_MATRIX.md) — per-model risk.
- [docs/PLAN.md](docs/PLAN.md) — BCFD tickets, milestones, dev protocol,
  merge-queue rules, known traps.
- [docs/LIMITATIONS.md](docs/LIMITATIONS.md) — machine-readable capability
  status.

Preserved engineering references:

- [docs/PHYSICS.md](docs/PHYSICS.md) — physics decisions + experiment
  log.
- [docs/ARCHITECTURE_V2.md](docs/ARCHITECTURE_V2.md) — code architecture.
- [docs/KERNEL_EXTENSION_POINTS.md](docs/KERNEL_EXTENSION_POINTS.md).
- [docs/MPI_GUIDE.md](docs/MPI_GUIDE.md) ·
  [docs/CLUSTER_OPTIONS.md](docs/CLUSTER_OPTIONS.md) ·
  [docs/CLUSTER_RUNBOOK.md](docs/CLUSTER_RUNBOOK.md).
- [docs/REQ_STIRRED_REACTOR.md](docs/REQ_STIRRED_REACTOR.md) — pre-pivot
  stirred-reactor requirements text, useful as bioprocess reference.
- [docs/T15_5_CAVITY3D_REFERENCE.md](docs/T15_5_CAVITY3D_REFERENCE.md).

Archive:

- [docs/archive/2026-07-07-pivot/](docs/archive/2026-07-07-pivot/) —
  pre-pivot PLAN, VALIDATION, LIMITATIONS; T1..T18 matrix; M-A..M-F track;
  R-Phase spec; V&V ledger; whitepaper; competitor analysis. Read only for
  the pre-pivot history.

## License

Dual-licensed under MIT OR Apache-2.0. Contributions are accepted under the
same terms.
