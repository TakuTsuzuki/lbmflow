//! Scenario execution: probes, snapshots, divergence/steady detection,
//! machine-readable manifest.

use crate::render::write_png;
use anyhow::{Context, Result};
use lbm_core::multiphase::ShanChen;
use lbm_core::prelude::*;
use lbm_scenario::{
    FieldKind, OutputFormat, OutputSpec, ProbeSpec, Scenario, SimHandle,
};
use serde::Serialize;
use std::fs;
use std::io::Write as _;
use std::path::{Path, PathBuf};
use std::time::Instant;

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
pub struct Manifest {
    pub scenario: String,
    pub status: String, // completed | steady | diverged
    pub steps_run: u64,
    pub wall_seconds: f64,
    pub mlups: f64,
    pub diagnostics: Diagnostics,
    pub warnings: Vec<lbm_scenario::Warning>,
    pub files: Vec<String>,
}

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
pub struct Diagnostics {
    pub total_mass: f64,
    pub max_speed: f64,
    pub tau: f64,
}

pub fn run(sc: &Scenario, out_dir: &Path) -> Result<Manifest> {
    fs::create_dir_all(out_dir)
        .with_context(|| format!("出力ディレクトリを作成できません: {}", out_dir.display()))?;
    match lbm_scenario::build(sc)? {
        SimHandle::F64(sim, mp) => run_t(sc, sim, mp, out_dir),
        SimHandle::F32(sim, mp) => run_t(sc, sim, mp, out_dir),
    }
}

struct CsvProbe {
    every: usize,
    file: fs::File,
    kind: &'static str,
    point: (usize, usize),
}

fn run_t<T: Real>(
    sc: &Scenario,
    mut sim: Simulation<T>,
    mp: Option<ShanChen<T>>,
    out_dir: &Path,
) -> Result<Manifest> {
    let mut files: Vec<String> = Vec::new();
    let mut probes: Vec<CsvProbe> = Vec::new();
    for p in &sc.probes {
        match *p {
            ProbeSpec::Force { every } => {
                let path = out_dir.join("force.csv");
                let mut file = fs::File::create(&path)?;
                writeln!(file, "step,fx,fy")?;
                files.push("force.csv".into());
                probes.push(CsvProbe {
                    every: every.max(1),
                    file,
                    kind: "force",
                    point: (0, 0),
                });
            }
            ProbeSpec::Point { x, y, every } => {
                let name = format!("point_{x}_{y}.csv");
                let path = out_dir.join(&name);
                let mut file = fs::File::create(&path)?;
                writeln!(file, "step,ux,uy,rho")?;
                files.push(name);
                probes.push(CsvProbe {
                    every: every.max(1),
                    file,
                    kind: "point",
                    point: (x, y),
                });
            }
        }
    }

    let t0 = Instant::now();
    let mut status = "completed";
    let mut prev_steady: Option<Vec<f64>> = None;
    let total = sc.run.steps;
    let mut executed: usize = 0;

    'main: for step in 1..=total {
        if let Some(mp) = &mp {
            mp.update_force(&mut sim);
        }
        sim.step();
        executed = step;

        for p in &mut probes {
            if step % p.every == 0 {
                match p.kind {
                    "force" => {
                        let f = sim.probed_force();
                        writeln!(p.file, "{step},{},{}", f[0].as_f64(), f[1].as_f64())?;
                    }
                    _ => {
                        let (x, y) = p.point;
                        writeln!(
                            p.file,
                            "{step},{},{},{}",
                            sim.ux(x, y).as_f64(),
                            sim.uy(x, y).as_f64(),
                            sim.rho(x, y).as_f64()
                        )?;
                    }
                }
            }
        }

        for (i, o) in sc.outputs.iter().enumerate() {
            if o.every > 0 && step % o.every == 0 {
                files.push(write_output(&sim, o, i, step, out_dir)?);
            }
        }

        if step % 1000 == 0 {
            let bad = sim
                .rho_field()
                .iter()
                .any(|v| !v.as_f64().is_finite());
            if bad {
                status = "diverged";
                break 'main;
            }
        }

        if let Some(ss) = &sc.run.stop_when_steady {
            if step % ss.check_every == 0 {
                let cur: Vec<f64> = sim.ux_field().iter().map(|v| v.as_f64()).collect();
                if let Some(prev) = &prev_steady {
                    let mut dmax = 0.0f64;
                    let mut umax = 0.0f64;
                    for (c, p) in cur.iter().zip(prev) {
                        dmax = dmax.max((c - p).abs());
                        umax = umax.max(c.abs());
                    }
                    if umax > 0.0 && dmax <= ss.epsilon * umax {
                        status = "steady";
                        break 'main;
                    }
                }
                prev_steady = Some(cur);
            }
        }
    }

    // end-of-run outputs (every == 0)
    for (i, o) in sc.outputs.iter().enumerate() {
        if o.every == 0 {
            files.push(write_output(&sim, o, i, executed, out_dir)?);
        }
    }

    let wall = t0.elapsed().as_secs_f64();
    let cells = (sim.nx() * sim.ny()) as f64;
    let max_speed = sim
        .ux_field()
        .iter()
        .zip(sim.uy_field())
        .map(|(a, b)| (a.as_f64().powi(2) + b.as_f64().powi(2)).sqrt())
        .fold(0.0f64, f64::max);
    let manifest = Manifest {
        scenario: sc.name.clone(),
        status: status.to_string(),
        steps_run: executed as u64,
        wall_seconds: wall,
        mlups: cells * executed as f64 / wall.max(1e-9) / 1e6,
        diagnostics: Diagnostics {
            total_mass: sim.total_mass().as_f64(),
            max_speed,
            tau: sim.tau(),
        },
        warnings: lbm_scenario::validate(sc),
        files,
    };
    fs::write(
        out_dir.join("manifest.json"),
        serde_json::to_string_pretty(&manifest)?,
    )?;
    Ok(manifest)
}

fn field_values<T: Real>(sim: &Simulation<T>, kind: FieldKind) -> Vec<f64> {
    let (nx, ny) = (sim.nx(), sim.ny());
    let idx = |x: usize, y: usize| y * nx + x;
    match kind {
        FieldKind::Ux => sim.ux_field().iter().map(|v| v.as_f64()).collect(),
        FieldKind::Uy => sim.uy_field().iter().map(|v| v.as_f64()).collect(),
        FieldKind::Rho => sim.rho_field().iter().map(|v| v.as_f64()).collect(),
        FieldKind::Speed => sim
            .ux_field()
            .iter()
            .zip(sim.uy_field())
            .map(|(a, b)| (a.as_f64().powi(2) + b.as_f64().powi(2)).sqrt())
            .collect(),
        FieldKind::Vorticity => {
            let ux: Vec<f64> = sim.ux_field().iter().map(|v| v.as_f64()).collect();
            let uy: Vec<f64> = sim.uy_field().iter().map(|v| v.as_f64()).collect();
            let mut w = vec![0.0; nx * ny];
            for y in 1..ny - 1 {
                for x in 1..nx - 1 {
                    w[idx(x, y)] = 0.5
                        * ((uy[idx(x + 1, y)] - uy[idx(x - 1, y)])
                            - (ux[idx(x, y + 1)] - ux[idx(x, y - 1)]));
                }
            }
            w
        }
    }
}

fn write_output<T: Real>(
    sim: &Simulation<T>,
    o: &OutputSpec,
    index: usize,
    step: usize,
    out_dir: &Path,
) -> Result<String> {
    let kind_name = format!("{:?}", o.field).to_lowercase();
    let values = field_values(sim, o.field);
    let solid: Vec<bool> = sim.solid_field().to_vec();
    let (nx, ny) = (sim.nx(), sim.ny());
    match o.format {
        OutputFormat::Png => {
            let name = format!("{kind_name}_{step}.png");
            let path: PathBuf = out_dir.join(&name);
            let diverging = matches!(o.field, FieldKind::Vorticity | FieldKind::Ux | FieldKind::Uy);
            write_png(&path, &values, &solid, nx, ny, diverging)?;
            Ok(name)
        }
        OutputFormat::Csv => {
            let name = format!("{kind_name}_{step}.csv");
            let mut file = fs::File::create(out_dir.join(&name))?;
            writeln!(file, "# {kind_name}, nx={nx}, ny={ny}, row-major y*nx+x, step={step}")?;
            for y in 0..ny {
                let row: Vec<String> = (0..nx)
                    .map(|x| format!("{}", values[y * nx + x]))
                    .collect();
                writeln!(file, "{}", row.join(","))?;
            }
            let _ = index;
            Ok(name)
        }
        OutputFormat::Vtk => {
            // VTK legacy structured points (ASCII). Point order is x-fastest,
            // which matches the row-major field layout (y up, like the sim).
            let name = format!("{kind_name}_{step}.vtk");
            let mut file = std::io::BufWriter::new(fs::File::create(out_dir.join(&name))?);
            writeln!(file, "# vtk DataFile Version 3.0")?;
            writeln!(file, "LBMFlow {kind_name} step={step}")?;
            writeln!(file, "ASCII")?;
            writeln!(file, "DATASET STRUCTURED_POINTS")?;
            writeln!(file, "DIMENSIONS {nx} {ny} 1")?;
            writeln!(file, "ORIGIN 0 0 0")?;
            writeln!(file, "SPACING 1 1 1")?;
            writeln!(file, "POINT_DATA {}", nx * ny)?;
            writeln!(file, "SCALARS {kind_name} double 1")?;
            writeln!(file, "LOOKUP_TABLE default")?;
            for chunk in values.chunks(9) {
                let line: Vec<String> = chunk.iter().map(|v| format!("{v}")).collect();
                writeln!(file, "{}", line.join(" "))?;
            }
            Ok(name)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn vtk_output_is_legacy_structured_points() {
        let sc: Scenario = serde_json::from_str(
            r#"{
                "name": "vtk-smoke",
                "grid": { "nx": 16, "ny": 12 },
                "physics": { "nu": 0.05 },
                "edges": {
                    "left": { "type": "bounceBack" },
                    "right": { "type": "bounceBack" },
                    "bottom": { "type": "bounceBack" },
                    "top": { "type": "movingWall", "u": [0.05, 0.0] }
                },
                "run": { "steps": 5 },
                "outputs": [ { "field": "rho", "format": "vtk", "every": 0 } ]
            }"#,
        )
        .unwrap();
        let dir = std::env::temp_dir().join(format!("lbm_vtk_test_{}", std::process::id()));
        let manifest = run(&sc, &dir).unwrap();
        assert!(manifest.files.contains(&"rho_5.vtk".to_string()), "{:?}", manifest.files);
        let text = fs::read_to_string(dir.join("rho_5.vtk")).unwrap();
        let mut lines = text.lines();
        assert_eq!(lines.next(), Some("# vtk DataFile Version 3.0"));
        assert!(text.contains("\nASCII\n"), "missing ASCII marker");
        assert!(text.contains("\nDATASET STRUCTURED_POINTS\n"));
        assert!(text.contains("\nDIMENSIONS 16 12 1\n"));
        assert!(text.contains("\nPOINT_DATA 192\n"));
        assert!(text.contains("\nSCALARS rho double 1\n"));
        assert!(text.contains("\nLOOKUP_TABLE default\n"));
        // all 16*12 values present after the header, and all parse as f64
        let data: Vec<f64> = text
            .lines()
            .skip_while(|l| *l != "LOOKUP_TABLE default")
            .skip(1)
            .flat_map(|l| l.split_whitespace())
            .map(|v| v.parse().unwrap())
            .collect();
        assert_eq!(data.len(), 192);
        assert!(data.iter().all(|v| v.is_finite()));
        fs::remove_dir_all(&dir).ok();
    }
}
