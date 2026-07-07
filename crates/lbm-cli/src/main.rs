//! `lbm` — the LBMFlow scenario-runner CLI (entry point of Agent mode).
//!
//! Design principles for agents: self-description (`lbm schema` / `lbm presets`),
//! structured errors (JSON), determinism.

mod capabilities;
mod gallery;
mod manifest;
mod mcp;
mod output;
mod render;
mod report;
mod runner;
mod schema;
mod sweep;
mod validate;
mod verify;

use anyhow::{Context, Result};
use clap::{Parser, Subcommand};
use lbm_scenario::{BioprocessScenario, Scenario};
use std::path::PathBuf;

#[derive(Parser)]
#[command(
    name = "lbm",
    about = "LBMFlow lattice Boltzmann method simulator CLI",
    version
)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// Print the supported lattices, collisions, precision/storage modes, and backends
    Capabilities {
        /// Emit stable machine-readable JSON
        #[arg(long)]
        json: bool,
    },
    /// Run built-in verification checks
    Verify {
        /// Verification tier to run
        #[arg(long, value_enum, default_value_t = verify::VerifyTier::Quick)]
        tier: verify::VerifyTier,
        /// Emit stable machine-readable JSON
        #[arg(long)]
        json: bool,
    },
    /// Run a scenario JSON and write results to the output directory
    Run {
        /// Scenario JSON file (`-` for stdin)
        scenario: String,
        /// Output directory (default: out/<scenario name>)
        #[arg(long)]
        out: Option<PathBuf>,
        /// Save checkpoints every N steps (0 disables checkpoint writes)
        #[arg(long)]
        save_every: Option<usize>,
        /// Directory where ckpt_<step>/ checkpoint folders are written
        #[arg(long)]
        checkpoint_dir: Option<PathBuf>,
        /// Restore from an existing checkpoint directory before running
        #[arg(long)]
        restore: Option<PathBuf>,
        /// Output strategy: gather writes legacy whole-field files; per-rank is MPI-only
        #[arg(long, value_enum)]
        output_mode: Option<output::OutputMode>,
        /// Print the result manifest to stdout as JSON
        #[arg(long)]
        json: bool,
    },
    /// Validate a scenario without running it (report errors/warnings as JSON)
    Validate {
        /// Scenario JSON file (`-` for stdin)
        scenario: String,
        /// Emit stable machine-readable JSON
        #[arg(long)]
        json: bool,
    },
    /// List, show, or run the built-in presets
    Presets {
        #[command(subcommand)]
        action: PresetAction,
    },
    /// Run all presets in sequence and generate a self-contained HTML gallery (index.html)
    Gallery {
        /// Output directory (default: out/gallery)
        #[arg(long)]
        out: Option<PathBuf>,
    },
    /// Print the scenario JSON format reference (for agent self-discovery)
    Schema {
        /// Emit the BioprocessScenario schema instead of the legacy Scenario reference
        #[arg(long)]
        bioprocess: bool,
    },
    /// Serve as an MCP server on stdio (AI agent integration)
    Mcp,
    /// Bioprocess reporting and decision utilities
    Bioprocess {
        #[command(subcommand)]
        action: BioprocessAction,
    },
}

#[derive(Subcommand)]
enum PresetAction {
    /// List the presets
    List,
    /// Show a preset's scenario JSON
    Show { name: String },
    /// Run a preset
    Run {
        name: String,
        #[arg(long)]
        out: Option<PathBuf>,
    },
}

#[derive(Subcommand)]
enum BioprocessAction {
    /// Generate report.md from a bioprocess run directory
    Report {
        /// Run directory containing qoi.json and manifest.json
        run_dir: PathBuf,
    },
    /// Run a deterministic bioprocess parameter sweep
    Sweep {
        /// Sweep JSON file
        sweep: PathBuf,
        /// Output directory for case folders and sweep_summary.json
        #[arg(long)]
        out: Option<PathBuf>,
    },
}

enum LoadedRunScenario {
    Legacy(Scenario),
    Bioprocess(BioprocessScenario),
}

fn load_run_scenario(path: &str) -> Result<LoadedRunScenario> {
    let text = if path == "-" {
        std::io::read_to_string(std::io::stdin())?
    } else {
        std::fs::read_to_string(path).with_context(|| format!("cannot read: {path}"))?
    };
    let value: serde_json::Value = serde_json::from_str(&text).map_err(|e| {
        anyhow::anyhow!(serde_json::to_string_pretty(&serde_json::json!({
            "error": "invalid-scenario-json",
            "message": e.to_string(),
            "hint": "see `lbm schema` or `lbm schema --bioprocess`"
        }))
        .unwrap())
    })?;
    if value.get("version").and_then(|v| v.as_str()) == Some("bioprocess-1.0") {
        let sc = BioprocessScenario::from_json_str(&text).map_err(|e| {
            anyhow::anyhow!(serde_json::to_string_pretty(&serde_json::json!({
                "error": "invalid-bioprocess-scenario",
                "message": e.message,
                "reason": e.reason
            }))
            .unwrap())
        })?;
        Ok(LoadedRunScenario::Bioprocess(sc))
    } else {
        let sc: Scenario = serde_json::from_value(value)?;
        Ok(LoadedRunScenario::Legacy(sc))
    }
}

fn run_and_report(
    sc: &Scenario,
    out: Option<PathBuf>,
    json: bool,
    options: runner::RunOptions,
) -> Result<()> {
    let out_dir = out.unwrap_or_else(|| PathBuf::from("out").join(&sc.name));
    let manifest = runner::run_with_options(sc, &out_dir, &options)?;
    if json {
        println!("{}", serde_json::to_string_pretty(&manifest)?);
    } else {
        println!(
            "status={} steps={} wall={:.1}s mlups={:.0} out={}",
            manifest.status,
            manifest.steps_run,
            manifest.wall_seconds,
            manifest.mlups,
            out_dir.display()
        );
        for w in &manifest.warnings {
            eprintln!("warning[{}]: {}", w.field, w.message);
        }
        for f in &manifest.files {
            println!("  {f}");
        }
    }
    Ok(())
}

fn warn_legacy_scenario_dispatch() {
    eprintln!("legacy LBM demo preset; not bioprocess decision-grade");
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    match cli.command {
        Command::Capabilities { json } => {
            capabilities::run(json)?;
        }
        Command::Verify { tier, json } => {
            let code = verify::run(tier, json)?;
            std::process::exit(code);
        }
        Command::Run {
            scenario,
            out,
            save_every,
            checkpoint_dir,
            restore,
            output_mode,
            json,
        } => match load_run_scenario(&scenario)? {
            LoadedRunScenario::Legacy(sc) => {
                warn_legacy_scenario_dispatch();
                run_and_report(
                    &sc,
                    out,
                    json,
                    runner::RunOptions {
                        save_every,
                        checkpoint_dir,
                        restore,
                        output_mode,
                    },
                )?;
            }
            LoadedRunScenario::Bioprocess(sc) => {
                anyhow::ensure!(
                    save_every.is_none() && checkpoint_dir.is_none() && restore.is_none(),
                    "bioprocess runner checkpoint options are not implemented in BCFD-030"
                );
                let out_dir = out.unwrap_or_else(|| PathBuf::from("out").join(&sc.name));
                let manifest = runner::run_bioprocess_single_phase(&sc, &out_dir)?;
                if json {
                    println!("{}", serde_json::to_string_pretty(&manifest)?);
                } else {
                    println!(
                        "status={} steps={} wall={:.1}s mlups={:.0} out={}",
                        manifest.status,
                        manifest.steps_run,
                        manifest.wall_seconds,
                        manifest.mlups,
                        out_dir.display()
                    );
                    for f in &manifest.files {
                        println!("  {f}");
                    }
                }
            }
        },
        Command::Validate { scenario, json } => {
            let code = validate::run(&scenario, json)?;
            std::process::exit(code);
        }
        Command::Presets { action } => match action {
            PresetAction::List => {
                for (name, desc, _) in lbm_scenario::presets() {
                    println!("{name:<20} {desc}");
                }
            }
            PresetAction::Show { name } => {
                let all = lbm_scenario::presets();
                let found = all.iter().find(|(n, _, _)| *n == name).ok_or_else(|| {
                    anyhow::anyhow!(
                        "no such preset '{name}'. See `lbm presets list` for available presets"
                    )
                })?;
                println!("{}", serde_json::to_string_pretty(&found.2)?);
            }
            PresetAction::Run { name, out } => {
                let all = lbm_scenario::presets();
                let found = all
                    .iter()
                    .find(|(n, _, _)| *n == name)
                    .ok_or_else(|| anyhow::anyhow!("no such preset '{name}'"))?;
                warn_legacy_scenario_dispatch();
                run_and_report(&found.2, out, false, runner::RunOptions::default())?;
            }
        },
        Command::Gallery { out } => {
            let out_root = out.unwrap_or_else(|| PathBuf::from("out").join("gallery"));
            warn_legacy_scenario_dispatch();
            gallery::run(&out_root)?;
        }
        Command::Schema { bioprocess } => {
            if bioprocess {
                println!("{}", schema::bioprocess_schema_json());
            } else {
                println!("{}", SCHEMA_DOC);
            }
        }
        Command::Mcp => {
            mcp::serve()?;
        }
        Command::Bioprocess { action } => match action {
            BioprocessAction::Report { run_dir } => {
                let path = report::generate_report(&run_dir)?;
                println!("{}", path.display());
            }
            BioprocessAction::Sweep { sweep, out } => {
                let path = sweep::run(&sweep, out)?;
                println!("{}", path.display());
            }
        },
    }
    Ok(())
}

const SCHEMA_DOC: &str = r#"Scenario JSON (v0) — lbm run <file.json>

{
  "version": 0,
  "name": "my-sim",                       // output directory name
  "grid": { "nx": 128, "ny": 128 },        // add "nz": 64 to run in 3D
                                           //   optional 3D "lattice": "d3q19" | "d3q27"
                                           //   absent lattice = d3q19; omitted/1 nz = 2D D2Q9
  "physics": {
    "nu": 0.02,                            // kinematic viscosity (lattice units); tau = 3*nu + 0.5
    "collision": { "type": "trt" },        // "trt" (recommended) | "bgk" | "central_moment"
                                           //   central_moment is currently exposed only on 3D CPU
    "force": [0.0, 0.0],                   // uniform body force (e.g. gravity; z component 0 in 3D)
    "precision": "f64"                     // "f32" | "f64"
  },
  "units": {                               // optional SI boundary converter; omitted = raw lattice units
    "constructor": "FromResolutionAndLatticeVelocity",
    "characteristicLength": 0.1,           // m
    "characteristicVelocity": 1.0,         // m/s
    "kinematicViscosity": 1.0e-6,          // m^2/s
    "density": 998.2,                      // kg/m^3, required when units is present
    "resolution": 200,                     // N
    "latticeVelocity": 0.1                 // or use relaxationTime; third constructor derives N
  },
  "compute": {                             // optional
    "backend": "auto",                     // "auto" | "cpu" | "gpu"
    "storage": "f32"                       // "f32" | "f16"; f16 requires 2D GPU + SHADER_F16
  },
  "edges": {                               // boundary conditions on the 4 edges
    "left":   { "type": "velocityInlet", "u": [0.1, 0.0] },
    "right":  { "type": "pressureOutlet", "rho": 1.0 },
    "bottom": { "type": "bounceBack" },
    "top":    { "type": "bounceBack" }
    // In 3D, "front" (z=0) / "back" (z=nz-1) may be added (omitted = periodic).
    // Others: {"type":"periodic"} (must be paired on opposite edges), {"type":"movingWall","u":[ux,uy]},
    //     {"type":"outflow"},
    //     {"type":"convectiveOutflow","uConv":0.1}  // convective outflow; less pressure reflection
    //       than outflow. uConv = expected mean outflow velocity (0 < uConv <= 1; ~inlet velocity is a good guide)
    // Constraint: open boundaries (velocityInlet/pressureOutlet/outflow/convectiveOutflow) must not
    //       be orthogonally adjacent to each other (in 3D, open boundaries on one axis only)
  },
  "inletProfile": {                        // optional: parabolic inflow
    "edge": "left", "kind": "parabolic", "umax": 0.15
    // 3D: parabolic along each wall-bounded tangential axis, uniform along periodic axes
    //     (with 4 walls: duct-like u = umax·f(y)·f(z))
  },
  "obstacles": [                           // optional
    { "shape": "circle", "cx": 80, "cy": 80, "r": 20 },   // in 3D, extruded along z (a cylinder)
    { "shape": "rect", "x0": 10, "y0": 10, "x1": 20, "y1": 40 },
    { "shape": "sphere", "cx": 60, "cy": 32, "cz": 32, "r": 12 }  // 3D only
  ],
  "init": { "kind": "rest" },              // rest | droplet{cx,cy,r,rhoLiquid,rhoVapor}
                                           // | pool{heightFrac,rhoLiquid,rhoVapor}
  "multiphase": {                          // optional: Shan-Chen single-component multiphase
    "g": -5.0,                             // cohesion strength (negative; -5.0 is the validated default)
    "gWall": 0.0,                          // wall adhesion (negative = wetting); prefer wallRho instead
    "wallRho": 1.0                         // optional: contact-angle control via virtual wall density
                                           // (toward liquid density → wetting. 0.3:~180°, 0.6:~107°, 1.0:~63°)
  },
  "run": {
    "steps": 20000,
    "stopWhenSteady": { "epsilon": 1e-8, "checkEvery": 500 }  // optional
  },
  "probes": [                              // optional: time-series CSV
    { "type": "force", "every": 10 },      // force on obstacles (force.csv; fx,fy,fz in 3D)
    { "type": "point", "x": 220, "y": 80, "every": 100 }  // in 3D "z" may also be set (omitted = nz/2)
  ],
  "outputs": [                             // optional: field snapshots
    { "field": "speed", "format": "png", "every": 0 }   // every=0 = only at the end
    // field: speed | ux | uy | rho | vorticity
    // format: png | csv | vtk (VTK legacy structured points; opens in ParaView etc.)
    // 3D: png/csv are the z mid-plane slice, vtk is the full 3D volume (DIMENSIONS nx ny nz)
  ]
}

3D (nz > 1) restrictions: single-phase only (no multiphase), init must be rest.
grid.lattice may be "d3q19" (default) or "d3q27"; D3Q27 supports CPU periodic, closed-wall,
velocity-inlet, pressure-outlet, outflow, and convective faces. D3Q27 GPU open-face scenarios
are rejected explicitly.
compute.backend must be cpu/auto in the current CLI runner. Engine is the V2 core (D3Q19/D3Q27).
CentralMoment collision is exposed on this 3D CPU path.
compute.storage f16 is GPU-storage-only and is rejected for CPU, 3D, and no-gpu-feature builds.

Results: <out>/manifest.json (status/steps/mlups/diagnostics/provenance/warnings/units/file list)
status: completed | steady (early stop on steady-state detection) | diverged (NaN detected)
Units constructors:
- FromResolutionAndRelaxationTime: resolution + relaxationTime; derives latticeVelocity.
- FromResolutionAndLatticeVelocity: resolution + latticeVelocity; derives relaxationTime.
- FromRelaxationTimeAndLatticeVelocity: relaxationTime + latticeVelocity; derives resolution.
`lbm validate` and MCP validate_scenario echo units{lattice,conversionFactors,dimensionless,diagnostics}.
Unit diagnostics use verdict "ok" | "warn" | "error"; validate exits non-zero on error.
Examples: lbm presets show cavity | cylinder-karman | two-phase-droplet | droplet-on-wall
Batch run: lbm gallery --out DIR (all presets + self-contained HTML gallery)

MCP: lbm mcp serves an MCP server (stdio). run_scenario runs synchronously (blocks until done).
For long runs or sweeps use the async API:
  start_run { scenario, outDir? } -> { runId }   ... responds immediately, runs in background
  run_status { runId } -> { state: running|completed|failed, manifest?, error? }
  list_runs {} -> list of all runs
runId is "run-<seq>-<scenario name>" (deterministic). At most 4 concurrent runs
(excess requests are rejected immediately with "failed: too many concurrent runs").
Runs live inside the server process, so keep the MCP connection open until completion.
"#;
