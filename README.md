# LBMFlow — Lattice Boltzmann Method Fluid Simulator

A fast LBM (Lattice Boltzmann Method, D2Q9) fluid simulation engine written in Rust,
with an integrated browser GUI and agent integration (CLI / MCP).

## Features

- **Validated physics**: passes an adversarially written validation suite including
  2nd-order convergence of the Taylor-Green vortex, exact Poiseuille flow with
  TRT (Λ=3/16), and the Ghia+1982 lid-driven cavity benchmark
  (spec: [docs/VALIDATION.md](docs/VALIDATION.md))
- **Explicit control over the accuracy/speed trade-off**:
  - Collision operator: BGK (fast) ⇔ TRT (accurate, more stable, recommended)
  - Numeric precision: `f32` (validation-grade via deviation storage) ⇔ `f64`
  - Parallelism: rayon multithreading (small grids automatically run serial)
  - GPU: wgpu/Metal measured at **6,975–12,152 MLUPS (16–42x over CPU)**;
    full integration in progress ([docs/GPU_EVALUATION.md](docs/GPU_EVALUATION.md))
- **Rich boundary conditions**: periodic, half-way bounce-back (stationary and
  moving walls), Zou-He velocity inlet (uniform or arbitrary profile), Zou-He
  pressure, zero-gradient outflow, arbitrarily shaped internal obstacles,
  momentum-exchange drag measurement
- **Three ways to use it**: browser GUI (WASM) / CLI (JSON scenarios) /
  MCP server (driven by AI agents)
- **Multiphase support**: Shan-Chen single-component multiphase (droplets,
  contact angles — validated)

## Usage 1: Browser GUI (beginner-friendly)

```bash
cd web && npm install && npm run dev   # → http://localhost:5173
```

Pick a preset (lid-driven cavity / Kármán vortex street / Poiseuille /
two-phase droplet / free canvas) and press "▶ Run". Obstacles can be painted
with the mouse. A real LBM engine (Rust→WASM) runs in the browser at roughly
600k lattice-site updates per second.

## Usage 2: CLI (scenario runner)

```bash
cargo build --release -p lbm-cli
./target/release/lbm presets list                 # built-in presets (4)
./target/release/lbm presets run cylinder-karman  # run → PNG/CSV/VTK/manifest.json in out/
./target/release/lbm gallery                      # run all presets → HTML report
./target/release/lbm schema                       # scenario JSON format
./target/release/lbm run my-scenario.json --json  # your own scenario
```

## Usage 3: MCP server (AI agent integration)

```bash
claude mcp add lbmflow -- /path/to/target/release/lbm mcp
```

Agents can run simulations through the 4 tools `run_scenario` /
`validate_scenario` / `list_presets` / `get_schema` and receive structured
results (manifest + PNG/CSV).

## Quick start (library)

```rust
use lbm_core::prelude::*;

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

## Development

```bash
cargo test --release                       # validation suite (always --release)
cargo test --release -- --include-ignored  # full validation incl. heavy benchmarks
```

- Plan & team structure: [docs/PLAN.md](docs/PLAN.md)
- Validation spec (acceptance criteria): [docs/VALIDATION.md](docs/VALIDATION.md)
- Physics models & experiment records: [docs/PHYSICS.md](docs/PHYSICS.md)
- Agent-mode design: [docs/AGENT_MODE_DESIGN.md](docs/AGENT_MODE_DESIGN.md)
- Multiphase design: [docs/MULTIPHASE_DESIGN.md](docs/MULTIPHASE_DESIGN.md)

## License

MIT OR Apache-2.0
