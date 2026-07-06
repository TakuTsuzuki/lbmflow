# LBMFlow

A commercial-grade lattice Boltzmann fluid simulator. Rust core, WebAssembly
GUI, a scenario CLI, and a native MCP server so AI agents can drive it.

D2Q9 / D3Q19 / D3Q27 lattices · BGK / TRT / cumulant collisions · scalar,
SIMD-fused, and wgpu GPU backends · FP16 storage · MPI halo exchange · Shan-Chen
multiphase · WALE LES · rotating immersed boundary · Bouzidi curved walls ·
Guo forcing.

Licensed under MIT OR Apache-2.0.

## Design principles

- **Explicit accuracy–speed control.** Every trade-off — collision operator,
  precision, backend, resolution — is a first-class knob, not a hidden default.
  BGK is fast; TRT with magic parameter Λ = 3/16 reproduces plane Poiseuille
  exactly on the half-way bounce-back grid; cumulant restores Galilean
  invariance at high Re. Pick your point on the curve; the trade-off is
  measured, not asserted.
- **Physically rigorous.** Every model term is derived from the governing
  equations or a literature-backed closure with a recorded derivation, validity
  domain, and its own validation test (`docs/PHYSICS.md`). Constants calibrated
  to pass a band, case-keyed branches, silent clamps that absorb transport, and
  decorative physics are prohibited by policy — if a gate cannot be met without
  a hack, the spec is revised, not the physics.
- **Validated adversarially.** The validation suite (`docs/VALIDATION.md`,
  T1–T18.x) is authored independently of the engine from a public spec. The
  engine is fixed until the tests pass — not the other way around.
- **Three front-ends over one core.** Browser GUI (Rust → WASM), scenario CLI,
  and an MCP server. The physics kernels are written once, generic over
  lattice and precision, and specialised at compile time.

## Highlights

- **Second-order Taylor–Green convergence** (T1); exact Poiseuille with
  TRT(Λ=3/16) (T2); the Ghia et al. (1982) lid-driven cavity benchmark to
  Re = 1000 (T7); Schäfer–Turek cylinder drag & Strouhal (T8).
- **3D GPU (D3Q19)**: quiet-window measurement 2 791–2 813 MLUPS at 192³,
  2 778–2 880 MLUPS at 128³ on an Apple M5 Max (A/B/A interleaved), against a
  1 500 MLUPS acceptance gate. T14 backend-equivalence (CPU ↔ GPU) verified
  to ≤ 1 × 10⁻⁵ on 32³ TGV and 24³ cavity.
- **FP16 storage** doubles the addressable grid and delivers ≈ 2× MLUPS at
  2048² with validation bands frozen to measured error (TGV transient
  1.401 × 10⁻¹, cavity steady 2.579 × 10⁻³); D3Q19 f16 exceeds 5 GLUPS.
- **Multiphase**: Shan-Chen single-component (droplets, Laplace law, full
  contact-angle range via virtual wall density) and two-component MCMP
  (Rayleigh–Taylor growth rate) — measurement-calibrated, T11/T11b/T11c/T12.
- **Turbulence**: WALE subgrid-scale model with near-wall damping recovered
  by construction; MKM 1999 channel-flow reference profiles landed for
  Re_τ = 180 characterisation.
- **Rotating machinery**: rotating immersed-boundary method for impellers and
  stirred-reactor geometries; dispersed-phase deposition tracking (D-track)
  with adhesion-capture and resuspension closures.
- **Curved walls**: Bouzidi second-order interpolated bounce-back.
- **Rich boundary catalogue**: periodic, half-way bounce-back (static /
  moving), Zou-He velocity inlet (uniform or `set_inlet_profile`), Zou-He
  pressure, zero-gradient and convective outflow, arbitrary internal
  obstacles, momentum-exchange force probes.
- **Bit-reproducible across backends and partitions**: T13 (partition
  invariance) and T14 (backend equivalence) are gate-tested every commit.

## Getting started

### Browser (no build required for users, engine ships as committed WASM)

```bash
cd web && npm install && npm run dev   # → http://localhost:5173
```

Pick a preset (lid-driven cavity / Kármán vortex street / two-phase droplet /
droplet-on-wall / free canvas) and press run. Obstacles can be drawn with
the mouse. This is the same Rust LBM engine, compiled to WebAssembly.

### CLI

```bash
cargo build --release -p lbm-cli
./target/release/lbm presets list             # cavity, cylinder-karman, two-phase-droplet, droplet-on-wall
./target/release/lbm presets show cavity      # print the preset's scenario JSON
./target/release/lbm presets run cylinder-karman   # → out/cylinder-karman/ (PNG + CSV + VTK + manifest.json)
./target/release/lbm gallery                  # run all presets, emit an HTML report
./target/release/lbm schema                   # scenario JSON schema
./target/release/lbm run my-scenario.json     # your own scenario
```

### MCP server (AI-agent integration)

```bash
claude mcp add lbmflow -- /path/to/target/release/lbm mcp
```

Seven tools: `run_scenario` (synchronous), `start_run` / `run_status` /
`list_runs` (async jobs for long runs and parallel sweeps), plus
`validate_scenario`, `list_presets`, `get_schema`. Results are structured
(manifest + PNG / CSV / VTK).

### Library

```rust
use lbm_core::compat::prelude::*;   // stable 2D facade

let mut sim: Simulation<f64> = SimConfig {
    nx: 128, ny: 128,
    nu: 0.02,
    collision: Collision::Trt { magic: 0.1875 },
    edges: Edges {
        left:   EdgeBC::BounceBack,
        right:  EdgeBC::BounceBack,
        bottom: EdgeBC::BounceBack,
        top:    EdgeBC::MovingWall { u: [0.1, 0.0] },
    },
    ..Default::default()
}.build()?;

sim.run(10_000);
println!("centre velocity = {}", sim.ux(64, 64));
```

The native V2 core API (`lbm_core::prelude` — `Solver`, `GlobalSpec`,
D2Q9 / D3Q19 / D3Q27, backend selection) is documented in
[docs/ARCHITECTURE_V2.md](docs/ARCHITECTURE_V2.md).

## Trade-off axes

| Axis        | Choices                                                    | Notes |
|-------------|------------------------------------------------------------|-------|
| Dimension   | 2D (D2Q9), 3D (D3Q19, D3Q27)                               | Compile-time specialisation over `Lattice`. |
| Collision   | BGK, TRT (magic Λ=3/16 default), cumulant (central-moment) | TRT for accuracy on Poiseuille and BCs; cumulant for high-Re Galilean invariance. |
| Precision   | `f32` (deviation storage, validation-grade), `f64`, `f16`  | `f16` storage doubles capacity at ≈ 2× MLUPS with frozen bands. |
| Backend     | `CpuScalar`, `CpuSimd` (fused collide+stream+moments), `Wgpu` | `--features gpu` for wgpu; T14 verifies CPU ↔ GPU equivalence. |
| Parallelism | rayon threads (auto-serial on small grids), MPI ranks      | `--features mpi` for domain-decomposed halo exchange. |

## Capability matrix

The core API is broader than the product path exposed through scenario JSON,
the CLI, and MCP. Use this matrix to choose combinations that are usable today.

| Capability | Core API | Scenario JSON + MCP path | Notes |
|------------|----------|--------------------------|-------|
| Lattices | D2Q9 and D3Q19 support periodic, closed-wall, and open-face runs; D3Q27 is available for periodic / closed-wall runs only. | `grid.nz <= 1` uses D2Q9; `grid.nz > 1` uses D3Q19. D3Q27 is not exposed in scenario JSON. | D3Q27 open faces are rejected as `UnsupportedOpenFaceLattice`. |
| Collision | `CollisionKind` exposes BGK, TRT, and cumulant / central-moment collision. | `CollisionSpec` exposes only `bgk` and `trt`. | Cumulant is core-only today. |
| Precision / storage | CPU solvers run with `f32` or `f64`; the GPU backend runs `f32` compute and can store distributions as `f32` or `f16`. | `Precision` exposes only `f32` and `f64`; GPU scenario dispatch rejects `f64`. | `f16` is a core GPU storage mode, not a scenario JSON precision. |
| Backends | `CpuScalar`, `CpuSimd`, `WgpuBackend` under feature `gpu`, and `MpiSolver` / `MpiExchange` under feature `mpi`. | CPU is the default scenario path. `compute.backend:"gpu"` is wired for 2D f32 scenarios when built with feature `gpu`; the CLI/MCP runner does not run 3D GPU scenarios yet. No JSON/MCP selector exists for MPI. | 3D GPU scenario capability rejects f64, multiphase, rotor, particles, non-rest initialization, and force probes; the core GPU backend also rejects localized volume sources and masked face patches. |
| Checkpoint / restart | `Solver::save` / `load` / `restore` are single-local-part only. Serialized state covers populations, stale buffer, moments, solid mask, and optional force field. | CLI run options can save and restore checkpoints for CPU-backed 2D and 3D scenario runs; MCP tools do not expose checkpoint options. | Multi-part / multi-rank checkpoint returns `CKPT_UNSUPPORTED`; manifest reserves `rng`, `particles`, and `stats` as not serialized. |
| Particles | `ParticleSet` is deterministic one-way Lagrangian coupling only. | Scenario particles are optional in 2D CPU runs; 3D and GPU scenario paths reject particles. | Particles feel sampled flow, gravity / buoyancy, and drag; they do not apply reaction force back to the fluid. |
| Multi-node scaling | MPI halo exchange and `MpiSolver` exist under feature `mpi`. | Not exposed through scenario JSON or MCP. | 64-rank weak-scaling acceptance remains RED pending cluster measurement. |

<!-- Capability matrix verification:
- Lattices/core: crates/lbm-core/src/lattice.rs defines D2Q9, D3Q19, D3Q27; crates/lbm-core/src/solver.rs validate_lattice rejects open faces unless the lattice has 3 or 5 unknowns, surfacing UnsupportedOpenFaceLattice for D3Q27; crates/lbm-core/src/solver.rs tests cover D3Q27 periodic / walled boxes and open-face rejection.
- Lattices/scenario: crates/lbm-scenario/src/lib.rs Grid::is_3d, Solver3<T> = Solver<D3Q19, ...>, build3d, and build_t/Simulation<D2Q9 compat path.
- Collision/core: crates/lbm-core/src/params.rs CollisionKind::{Bgk, Trt, Cumulant}; Collision/scenario: crates/lbm-scenario/src/lib.rs CollisionSpec::{Bgk, Trt}.
- Precision/core: crates/lbm-core/src/real.rs implements Real for f32/f64; crates/lbm-core/src/gpu/backend.rs GpuStorage::{F32, F16}, KernelCfg, WgpuBackend::with_config. Precision/scenario: crates/lbm-scenario/src/lib.rs Precision::{F32, F64} and gpu_capability_error rejects f64 GPU dispatch.
- Backends/core: crates/lbm-core/src/backend.rs CpuScalar Backend impl; crates/lbm-core/src/backend_simd.rs CpuSimd Backend impl; crates/lbm-core/src/gpu/backend.rs WgpuBackend Backend<L, f32> impl; crates/lbm-core/src/dist.rs MpiExchange/MpiSolver; crates/lbm-core/src/solver.rs rejects unsupported localized GPU features.
- Backends/scenario: crates/lbm-scenario/src/lib.rs BackendSpec and gpu_capability_error; crates/lbm-cli/src/runner.rs dispatches GPU only through run_gpu2d and bails for 3D GPU.
- Checkpoint/core: crates/lbm-core/src/solver.rs Solver::save/load/restore, CKPT_UNSUPPORTED single-part guard, serialized sections, and reserved rng/particles/stats false. Checkpoint/scenario: crates/lbm-cli/src/runner.rs RunOptions restore/save paths in run_t and run3d_t.
- Particles: crates/lbm-core/src/particles.rs module docs and ParticleSet; crates/lbm-scenario/src/lib.rs ParticlesSpec comment plus build3d/build_gpu2d rejections; crates/lbm-cli/src/runner.rs one-way particle stepping.
- Multi-node: crates/lbm-core/src/dist.rs feature-gated MPI module; docs/paper/claims-ledger.md keeps Multi-node scaling RED; README measured status retains the RED row.
-->

## Measured status

Working snapshot from `docs/paper/claims-ledger.md`:

| Capability                                     | Gate                                          | Status |
|------------------------------------------------|-----------------------------------------------|:------:|
| 3D GPU D3Q19 (T14-3D + ≥ 1 500 MLUPS)          | 32³ TGV3D u ≤ 1 × 10⁻⁵ · MLUPS quiet-window   | GREEN  |
| Explicit `backend:"gpu"` in 2D scenarios       | End-to-end honoured                           | GREEN  |
| FP16 storage, × 2 grid capacity                | T16 bands frozen · ≥ 1.5× MLUPS @ 2048²       | GREEN  |
| 2D GPU GLUPS · CPU MLUPS · T13 bit-exact       | Landed and measured                           | GREEN  |
| WASM bit-identity · agent-native MCP + Skills  | Landed                                        | GREEN  |
| Multi-node weak scaling ≥ 80 % @ 64 rank       | 64-rank weak measurement                      |  RED   |
| Full-physics stirred workload                  | Degradation ratio vs single-phase             |  RED   |

RED rows track external inputs (cluster access, M-F integration
completion), not implementation velocity.

## Building and testing

```bash
cargo build --workspace --release
cargo test  --workspace --release              # the default gate — always --release
cargo test  --release -- --include-ignored     # + heavy validation and benches (~5 min)

# WebAssembly for the browser GUI (lbm-wasm is outside the workspace):
wasm-pack build crates/lbm-wasm --target web --release --out-dir ../../web/src/engine/pkg
cd web && npm run build

# Optional features:
cargo test  --workspace --release --features gpu   # wgpu backend (GPU hosts only)
cargo test  --workspace --release --features mpi   # requires a native MPI toolchain
```

LBM is roughly 50× slower in debug builds; `--release` is not optional. The
default gate is `cargo test --workspace --release --no-fail-fast`.

## Repository map

- `crates/lbm-core` — the engine: D2Q9 / D3Q19 / D3Q27 lattices, CPU scalar /
  SIMD backends, wgpu GPU backend (feature `gpu`), MPI halo exchange
  (feature `mpi`), WALE LES, rotating IBM, Bouzidi. The legacy 2D public
  facade lives in `compat/`.
- `crates/lbm-scenario` — JSON scenario schema, runner, and unit system.
- `crates/lbm-cli` — the `lbm` binary: presets, gallery, schema, scenario
  run, and the MCP server (7 tools including async `start_run` /
  `run_status` / `list_runs`).
- `crates/lbm-wasm` — WASM bindings for the web GUI (excluded from the
  workspace; the built `pkg/` is committed under `web/src/engine/pkg`).
- `crates/lbm-gpu-proto` — standalone wgpu evaluation prototype (measurement
  record, superseded by the in-core `gpu` module).
- `web/` — the TypeScript / Vite GUI.

## Documentation

Physics and validation

- [docs/PHYSICS.md](docs/PHYSICS.md) — physics decisions and experiment log.
- [docs/VALIDATION.md](docs/VALIDATION.md) — validation-test specification (T1–T18.x).
- [docs/LIMITATIONS.md](docs/LIMITATIONS.md) — explicit capability limits and unsupported combinations.
- [docs/T15_5_CAVITY3D_REFERENCE.md](docs/T15_5_CAVITY3D_REFERENCE.md) — 3D cavity reference data.

Architecture and design

- [docs/ARCHITECTURE_V2.md](docs/ARCHITECTURE_V2.md) — dimension × lattice × precision × backend × partition.
- [docs/SOLVER_IMPROVEMENT_SPEC.md](docs/SOLVER_IMPROVEMENT_SPEC.md) — R-Phase solver spec.
- [docs/KERNEL_EXTENSION_POINTS.md](docs/KERNEL_EXTENSION_POINTS.md) — extending the kernel.
- [docs/MULTIPHASE_DESIGN.md](docs/MULTIPHASE_DESIGN.md) · [docs/WASM_BRIDGE_DESIGN.md](docs/WASM_BRIDGE_DESIGN.md) · [docs/AGENT_MODE_DESIGN.md](docs/AGENT_MODE_DESIGN.md).
- [docs/DISPERSED_DEPOSITION.md](docs/DISPERSED_DEPOSITION.md) — dispersed-phase deposition track.
- [docs/REQ_STIRRED_REACTOR.md](docs/REQ_STIRRED_REACTOR.md) — stirred-reactor requirements.

Performance and scale

- [docs/PERFORMANCE.md](docs/PERFORMANCE.md) — measured MLUPS, thread scaling, mode-selection guide.
- [docs/GPU_EVALUATION.md](docs/GPU_EVALUATION.md) — wgpu evaluation and kernel notes.
- [docs/BENCH_COMPARISON_DRAFT.md](docs/BENCH_COMPARISON_DRAFT.md) — external-comparison working draft.
- [docs/MPI_GUIDE.md](docs/MPI_GUIDE.md) · [docs/CLUSTER_OPTIONS.md](docs/CLUSTER_OPTIONS.md).

Programme

- [docs/PLAN.md](docs/PLAN.md) — milestones M-A … M-F and the active queue.
- [docs/paper/LBMFlow-whitepaper.md](docs/paper/LBMFlow-whitepaper.md) — living technical paper.
- [docs/paper/claims-ledger.md](docs/paper/claims-ledger.md) — measurement-status snapshot.

## License

Dual-licensed under MIT OR Apache-2.0. Contributions are accepted under the
same terms.
