# LBMFlow — Lattice Boltzmann Method Fluid Simulator

A fast Rust LBM engine (D2Q9 / D3Q19) with a browser GUI, a scenario CLI, and
AI-agent integration (MCP).

## Features

- **Validated physics**: the validation suite is authored adversarially,
  separate from the implementation — second-order Taylor-Green convergence,
  exact Poiseuille reproduction with TRT(Λ=3/16), the Ghia et al. (1982)
  lid-driven cavity benchmark, Shan-Chen multiphase (Laplace law, contact
  angles), and more (spec: [docs/VALIDATION.md](docs/VALIDATION.md)).
- **Explicit accuracy/speed trade-offs**:
  - Collision operator: BGK (fast) ⇔ TRT (accurate and stable, recommended)
  - Precision: `f32` (deviation storage, validation-grade) ⇔ `f64`
  - Parallelism: rayon multithreading (small grids fall back to serial),
    SIMD-fused CPU backend
  - GPU: wgpu compute backend (feature `gpu`, D2Q9/f32) — measured figures
    in [docs/PERFORMANCE.md](docs/PERFORMANCE.md) and
    [docs/GPU_EVALUATION.md](docs/GPU_EVALUATION.md)
- **2D and 3D**: compile-time lattice selection (D2Q9 / D3Q19); 3D lid-driven
  cavity validation is in progress
  ([docs/VALIDATION.md](docs/VALIDATION.md), T15.x).
- **Scale-out**: MPI domain decomposition behind feature `mpi`
  ([docs/MPI_GUIDE.md](docs/MPI_GUIDE.md)).
- **Rich boundary conditions**: periodic, half-way bounce-back (static and
  moving walls), Zou-He velocity inlet (uniform or arbitrary profile via
  `set_inlet_profile`), Zou-He pressure, zero-gradient and convective outflow,
  arbitrary internal obstacles, momentum-exchange drag probes.
- **Multiphase**: Shan-Chen single-component two-phase (droplets, contact
  angles — validated).
- **Three front-ends**: browser GUI (WASM) / CLI (JSON scenarios) /
  MCP server (driven by AI agents).

## Usage 1: browser GUI (beginner-friendly)

```bash
cd web && npm install && npm run dev   # → http://localhost:5173
```

Pick a preset (lid-driven cavity / Kármán vortex street / Poiseuille /
two-phase droplet / free canvas) and press run. Obstacles can be drawn with
the mouse. This is the real LBM engine (Rust→WASM) running in your browser.

## Usage 2: CLI (scenario runner)

```bash
cargo build --release -p lbm-cli
./target/release/lbm presets list                 # built-in presets
./target/release/lbm presets show cavity          # print a preset's scenario JSON
./target/release/lbm presets run cylinder-karman  # run → PNG/CSV/VTK/manifest.json in out/
./target/release/lbm gallery                      # run all presets → HTML report
./target/release/lbm schema                       # scenario JSON format
./target/release/lbm run my-scenario.json --json  # your own scenario
```

## Usage 3: MCP server (AI-agent integration)

```bash
claude mcp add lbmflow -- /path/to/target/release/lbm mcp
```

Agents drive simulations through seven tools: `run_scenario` (synchronous),
`start_run` / `run_status` / `list_runs` (async jobs for long runs and
parallel sweeps), plus `validate_scenario` / `list_presets` / `get_schema`.
Results are structured (manifest + PNG/CSV).

## Quickstart (library)

```rust
use lbm_core::compat::prelude::*;   // stable 2D facade

// Lid-driven cavity
let mut sim: Simulation<f64> = SimConfig {
    nx: 128, ny: 128,
    nu: 0.02,
    edges: Edges {
        left: EdgeBC::BounceBack,
        right: EdgeBC::BounceBack,
        bottom: EdgeBC::BounceBack,
        top: EdgeBC::MovingWall { u: [0.1, 0.0] },
    },
    ..Default::default()
}.build()?;

sim.run(10_000);
println!("centre velocity = {}", sim.ux(64, 64));
```

The native V2 core API (`lbm_core::prelude` — `Solver`, `GlobalSpec`,
D2Q9/D3Q19, backend selection) is documented in
[docs/ARCHITECTURE_V2.md](docs/ARCHITECTURE_V2.md).

## Development

```bash
cargo test --workspace --release           # validation suite (ALWAYS --release)
cargo test --release -- --include-ignored  # full validation incl. heavy benches
```

- Plan & team structure: [docs/PLAN.md](docs/PLAN.md)
- Validation spec (acceptance criteria): [docs/VALIDATION.md](docs/VALIDATION.md)
- Physics models & experiment log: [docs/PHYSICS.md](docs/PHYSICS.md)
- Core architecture (V2): [docs/ARCHITECTURE_V2.md](docs/ARCHITECTURE_V2.md)
- Performance measurements: [docs/PERFORMANCE.md](docs/PERFORMANCE.md)
- MPI / distributed runs: [docs/MPI_GUIDE.md](docs/MPI_GUIDE.md)
- Agent mode design: [docs/AGENT_MODE_DESIGN.md](docs/AGENT_MODE_DESIGN.md)
- Multiphase design: [docs/MULTIPHASE_DESIGN.md](docs/MULTIPHASE_DESIGN.md)

## License

MIT OR Apache-2.0
