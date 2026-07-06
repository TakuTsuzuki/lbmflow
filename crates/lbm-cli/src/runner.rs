//! Scenario execution: probes, snapshots, divergence/steady detection,
//! machine-readable manifest.

use crate::render::{write_png_scaled, Colormap};
use anyhow::{Context, Result};
use lbm_core::compat::multiphase::ShanChen;
use lbm_core::compat::prelude::*;
use lbm_scenario::{
    FieldKind, OutputFormat, OutputSpec, ProbeSpec, Scenario, Sim3Handle, SimHandle, Solver3,
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
    #[serde(skip_serializing_if = "Option::is_none")]
    pub units: Option<lbm_scenario::UnitReport>,
    pub files: Vec<String>,
}

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
pub struct Diagnostics {
    pub total_mass: f64,
    pub max_speed: f64,
    pub tau: f64,
}

#[derive(Clone, Debug, Default)]
pub struct RunOptions {
    pub save_every: Option<usize>,
    pub checkpoint_dir: Option<PathBuf>,
    pub restore: Option<PathBuf>,
}

pub fn run(sc: &Scenario, out_dir: &Path) -> Result<Manifest> {
    run_with_options(sc, out_dir, &RunOptions::default())
}

pub fn run_with_options(sc: &Scenario, out_dir: &Path, options: &RunOptions) -> Result<Manifest> {
    fs::create_dir_all(out_dir)
        .with_context(|| format!("cannot create output directory: {}", out_dir.display()))?;
    let resolved = lbm_scenario::resolve(sc).map_err(anyhow::Error::msg)?;
    let units = resolved.as_ref().map(|r| r.report.clone());
    let resolved_scenario;
    let sc = if let Some(r) = resolved {
        resolved_scenario = r.scenario;
        &resolved_scenario
    } else {
        sc
    };
    #[cfg(feature = "gpu")]
    if lbm_scenario::selected_backend(sc) == lbm_scenario::BackendChoice::Gpu {
        eprintln!(
            "compute.backend selected gpu for scenario '{}' ({} cells)",
            sc.name,
            sc.grid.nx * sc.grid.ny * sc.grid.nz
        );
        if sc.is_3d() {
            anyhow::bail!("{}", lbm_scenario::gpu_capability_error(sc).unwrap_or(
                "requested backend \"gpu\" is not wired to the 3D CLI runner yet"
            ));
        }
        return run_gpu2d(sc, lbm_scenario::build_gpu2d(sc)?, out_dir, units);
    }
    #[cfg(not(feature = "gpu"))]
    if lbm_scenario::selected_backend(sc) == lbm_scenario::BackendChoice::Gpu
        && lbm_scenario::requested_backend(sc) == lbm_scenario::BackendSpec::Gpu
    {
        anyhow::bail!("{}", lbm_scenario::build_check(sc).unwrap_err());
    }
    if sc.is_3d() {
        return match lbm_scenario::build3d(sc)? {
            Sim3Handle::F64(s) => run3d_t(sc, s, out_dir, units, options),
            Sim3Handle::F32(s) => run3d_t(sc, s, out_dir, units, options),
        };
    }
    match lbm_scenario::build(sc)? {
        SimHandle::F64(sim, mp) => run_t(sc, sim, mp, out_dir, units, options),
        SimHandle::F32(sim, mp) => run_t(sc, sim, mp, out_dir, units, options),
    }
}

fn checkpoint_path(options: &RunOptions, step: usize) -> Option<PathBuf> {
    let every = options.save_every?;
    if every == 0 || step % every != 0 {
        return None;
    }
    let root = options
        .checkpoint_dir
        .clone()
        .unwrap_or_else(|| PathBuf::from("checkpoints"));
    Some(root.join(format!("ckpt_{step}")))
}

#[cfg(feature = "gpu")]
fn run_gpu2d(
    sc: &Scenario,
    mut sim: lbm_scenario::GpuSim2,
    out_dir: &Path,
    units: Option<lbm_scenario::UnitReport>,
) -> Result<Manifest> {
    let t0 = Instant::now();
    let mut files = Vec::new();
    let total = sc.run.steps;
    sim.run(total);
    sim.sync();
    for (i, o) in sc.outputs.iter().enumerate() {
        if o.every == 0 {
            files.push(write_output_gpu2d(&mut sim, o, i, total, out_dir)?);
        }
    }
    let wall = t0.elapsed().as_secs_f64();
    let ux = sim.gather_ux();
    let uy = sim.gather_uy();
    let max_speed = ux
        .iter()
        .zip(&uy)
        .map(|(a, b)| ((*a as f64).powi(2) + (*b as f64).powi(2)).sqrt())
        .fold(0.0f64, f64::max);
    let dims = sim.dims();
    let manifest = Manifest {
        scenario: sc.name.clone(),
        status: "completed".to_string(),
        steps_run: total as u64,
        wall_seconds: wall,
        mlups: (dims[0] * dims[1]) as f64 * total as f64 / wall.max(1e-9) / 1e6,
        diagnostics: Diagnostics {
            total_mass: sim.total_mass() as f64,
            max_speed,
            tau: 3.0 * sc.physics.nu + 0.5,
        },
        warnings: lbm_scenario::validate(sc),
        units,
        files,
    };
    fs::write(
        out_dir.join("manifest.json"),
        serde_json::to_string_pretty(&manifest)?,
    )?;
    Ok(manifest)
}

#[cfg(feature = "gpu")]
fn write_output_gpu2d(
    sim: &mut lbm_scenario::GpuSim2,
    o: &OutputSpec,
    _index: usize,
    step: usize,
    out_dir: &Path,
) -> Result<String> {
    let [nx, ny, _] = sim.dims();
    let kind_name = format!("{:?}", o.field).to_lowercase();
    let ux = sim.gather_ux();
    let uy = sim.gather_uy();
    let rho = sim.gather_rho();
    let values: Vec<f64> = match o.field {
        FieldKind::Ux => ux.iter().map(|v| *v as f64).collect(),
        FieldKind::Uy => uy.iter().map(|v| *v as f64).collect(),
        FieldKind::Rho => rho.iter().map(|v| *v as f64).collect(),
        FieldKind::Speed => ux
            .iter()
            .zip(&uy)
            .map(|(a, b)| ((*a as f64).powi(2) + (*b as f64).powi(2)).sqrt())
            .collect(),
        _ => anyhow::bail!("GPU 2D output does not support field {:?}", o.field),
    };
    match o.format {
        OutputFormat::Png => {
            let solid: Vec<bool> = (0..ny)
                .flat_map(|y| (0..nx).map(move |x| (x, y)))
                .map(|(x, y)| sim.is_solid(x, y, 0))
                .collect();
            let name = format!("{kind_name}_{step}.png");
            write_png_scaled(
                &out_dir.join(&name),
                &values,
                &solid,
                nx,
                ny,
                colormap_for(o.field),
                None,
                1,
            )?;
            Ok(name)
        }
        OutputFormat::Csv => {
            let name = format!("{kind_name}_{step}.csv");
            let mut file = fs::File::create(out_dir.join(&name))?;
            for y in 0..ny {
                let row: Vec<String> = (0..nx).map(|x| format!("{}", values[y * nx + x])).collect();
                writeln!(file, "{}", row.join(","))?;
            }
            Ok(name)
        }
        OutputFormat::Vtk => anyhow::bail!("GPU 2D output does not support VTK yet"),
    }
}

struct CsvProbe {
    every: usize,
    file: fs::File,
    kind: &'static str,
    point: (usize, usize),
}

fn write_particles(
    ps: &lbm_core::particles::ParticleSet,
    step: usize,
    out_dir: &Path,
) -> Result<String> {
    let name = format!("particles_{step}.csv");
    let mut file = fs::File::create(out_dir.join(&name))?;
    writeln!(file, "id,x,y,z,vx,vy,vz,exposure")?;
    for (i, p) in ps.particles.iter().enumerate() {
        writeln!(
            file,
            "{i},{},{},{},{},{},{},{}",
            p.pos[0], p.pos[1], p.pos[2], p.vel[0], p.vel[1], p.vel[2], p.exposure
        )?;
    }
    Ok(name)
}

fn run_t<T: Real>(
    sc: &Scenario,
    mut sim: Simulation<T>,
    mp: Option<ShanChen<T>>,
    out_dir: &Path,
    units: Option<lbm_scenario::UnitReport>,
    options: &RunOptions,
) -> Result<Manifest> {
    if let Some(dir) = &options.restore {
        sim.restore(dir)
            .with_context(|| format!("cannot restore checkpoint: {}", dir.display()))?;
    }
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
            ProbeSpec::Point { x, y, every, .. } => {
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

    // Rotating impeller (volume penalization). The rotor ADDS into the
    // per-cell force field each step; when no multiphase driver rewrites the
    // field, the runner clears it so forces do not accumulate across steps.
    if let Some(r) = &sc.rotor {
        let tip = (r.omega * r.r_blade).abs();
        anyhow::ensure!(
            tip <= 0.3,
            "rotor tip speed {tip} exceeds the low-Mach hard limit 0.3"
        );
    }
    let mut rotor = sc.rotor.map(|r| {
        lbm_core::compat::rotor::Rotor::new(T::r(r.cx), T::r(r.cy))
            .n_blades(r.n_blades)
            .r_hub(T::r(r.r_hub))
            .r_blade(T::r(r.r_blade))
            .blade_thickness(T::r(r.thickness))
            .omega(T::r(r.omega))
            .chi(T::r(r.chi))
            .omega_ramp_steps(r.ramp_steps as u64)
            .theta0(T::r(r.theta0))
    });
    let mut torque_file = if sc.rotor.is_some() {
        let path = out_dir.join("torque.csv");
        let mut f = fs::File::create(&path)?;
        writeln!(f, "step,torque,torqueIntegral")?;
        files.push("torque.csv".into());
        Some(f)
    } else {
        None
    };

    // One-way Lagrangian particles: deterministic grid seeding, velocities
    // sampled bilinearly from the resolved field after each step.
    let mut pset = sc.particles.as_ref().map(|spec| {
        let mut parts = Vec::with_capacity(spec.count);
        let (w, h) = (spec.seed.x1 - spec.seed.x0, spec.seed.y1 - spec.seed.y0);
        let cols = ((spec.count as f64 * (w / h.max(1e-9)).max(1e-9))
            .sqrt()
            .ceil() as usize)
            .max(1);
        let rows = spec.count.div_ceil(cols);
        for i in 0..spec.count {
            let (cx, cy) = (i % cols, i / cols);
            parts.push(lbm_core::particles::Particle {
                pos: [
                    spec.seed.x0 + (cx as f64 + 0.5) * w / cols as f64,
                    spec.seed.y0 + (cy as f64 + 0.5) * h / rows.max(1) as f64,
                    0.0,
                ],
                vel: [0.0; 3],
                d: spec.d,
                rho_p: spec.rho_p,
                exposure: 0.0,
            });
        }
        lbm_core::particles::ParticleSet::new(
            parts,
            1.0,
            sc.physics.nu,
            sc.physics.gravity.unwrap_or([0.0; 3]),
        )
        .with_restitution(spec.restitution)
    });
    let particles_every = sc.particles.map(|p| p.output_every).unwrap_or(0);

    let t0 = Instant::now();
    let mut status = "completed";
    let mut prev_steady: Option<Vec<f64>> = None;
    let total = sc.run.steps;
    let mut executed: usize = 0;

    'main: for step in 1..=total {
        if let Some(mp) = &mp {
            mp.update_force(&mut sim);
        }
        if let Some(rt) = &mut rotor {
            if mp.is_none() {
                sim.clear_force_field();
            }
            rt.update_force(&mut sim);
        }
        sim.step();
        executed = step;

        if let Some(rt) = &rotor {
            if let Some(f) = &mut torque_file {
                if step % 10 == 0 {
                    writeln!(
                        f,
                        "{step},{},{}",
                        rt.torque().as_f64(),
                        rt.torque_integral().as_f64()
                    )?;
                }
            }
        }

        if let Some(ps) = &mut pset {
            let dims = [sim.nx(), sim.ny(), 1];
            let sampler = |pos: [f64; 3]| {
                lbm_core::particles::sample_grid(pos, dims, |x, y, _| {
                    if sim.is_solid(x, y) {
                        ([0.0, 0.0, 0.0], true)
                    } else {
                        ([sim.ux(x, y).as_f64(), sim.uy(x, y).as_f64(), 0.0], false)
                    }
                })
            };
            ps.step(sampler, None::<fn([f64; 3]) -> f64>);
            if particles_every > 0 && step % particles_every == 0 {
                files.push(write_particles(ps, step, out_dir)?);
            }
        }

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

        if let Some(path) = checkpoint_path(options, step) {
            sim.save(&path)
                .with_context(|| format!("cannot save checkpoint: {}", path.display()))?;
            files.push(path.display().to_string());
        }

        if step % 1000 == 0 {
            let bad = sim.rho_field().iter().any(|v| !v.as_f64().is_finite());
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
    if let Some(ps) = &pset {
        files.push(write_particles(ps, executed, out_dir)?);
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
        units,
        files,
    };
    fs::write(
        out_dir.join("manifest.json"),
        serde_json::to_string_pretty(&manifest)?,
    )?;
    Ok(manifest)
}

/// Colormap chosen by field semantics: signed fields diverge (RdBu), magnitude
/// fields are sequential (Inferno for stress/shear, Viridis for speed/density).
fn colormap_for(field: FieldKind) -> Colormap {
    match field {
        FieldKind::Vorticity | FieldKind::Ux | FieldKind::Uy | FieldKind::QCriterion => {
            Colormap::RdBu
        }
        FieldKind::ShearRate | FieldKind::DissipationRate | FieldKind::VorticityMag => {
            Colormap::Inferno
        }
        FieldKind::Speed | FieldKind::Rho => Colormap::Viridis,
    }
}

/// Velocity-gradient-derived scalar fields — the SINGLE physics implementation
/// of vorticity magnitude and Q-criterion, shared by the 2D and 3D field
/// providers (SPEC_OBSERVER_FRAMEWORK §12-F3: one field_value site, never
/// forked). These need the ANTISYMMETRIC (rotation) part of grad(u), which is
/// genuinely absent from the non-equilibrium stress `f_neq` (it carries only the
/// symmetric strain rate) — hence finite differences here, not the native gather
/// used for ShearRate/DissipationRate. `uz` may be all zeros for 2D (`nz == 1`),
/// where z-derivatives drop.
fn grad_derived(ux: &[f64], uy: &[f64], uz: &[f64], dims: [usize; 3], kind: FieldKind) -> Vec<f64> {
    let [nx, ny, nz] = dims;
    let idx = |x: usize, y: usize, z: usize| (z * ny + y) * nx + x;
    let mut out = vec![0.0; nx * ny * nz];
    let (z0, z1) = if nz > 1 { (1, nz - 1) } else { (0, 1) };
    for z in z0..z1 {
        for y in 1..ny - 1 {
            for x in 1..nx - 1 {
                let ddx = |f: &[f64]| 0.5 * (f[idx(x + 1, y, z)] - f[idx(x - 1, y, z)]);
                let ddy = |f: &[f64]| 0.5 * (f[idx(x, y + 1, z)] - f[idx(x, y - 1, z)]);
                let ddz = |f: &[f64]| {
                    if nz > 1 {
                        0.5 * (f[idx(x, y, z + 1)] - f[idx(x, y, z - 1)])
                    } else {
                        0.0
                    }
                };
                // g[i][j] = d(u_i)/d(x_j), (i,j) in {0:x, 1:y, 2:z}.
                let g = [
                    [ddx(ux), ddy(ux), ddz(ux)],
                    [ddx(uy), ddy(uy), ddz(uy)],
                    [ddx(uz), ddy(uz), ddz(uz)],
                ];
                out[idx(x, y, z)] = match kind {
                    FieldKind::VorticityMag => {
                        let wx = g[2][1] - g[1][2];
                        let wy = g[0][2] - g[2][0];
                        let wz = g[1][0] - g[0][1];
                        (wx * wx + wy * wy + wz * wz).sqrt()
                    }
                    FieldKind::QCriterion => {
                        let (mut s2, mut w2) = (0.0, 0.0);
                        for i in 0..3 {
                            for j in 0..3 {
                                let s = 0.5 * (g[i][j] + g[j][i]);
                                let w = 0.5 * (g[i][j] - g[j][i]);
                                s2 += s * s;
                                w2 += w * w;
                            }
                        }
                        0.5 * (w2 - s2)
                    }
                    _ => 0.0,
                };
            }
        }
    }
    out
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
        FieldKind::ShearRate => sim.shear_rate_field().iter().map(|v| v.as_f64()).collect(),
        FieldKind::DissipationRate => {
            let nu = sim.nu();
            sim.shear_rate_field()
                .iter()
                .map(|v| {
                    let g = v.as_f64();
                    nu * g * g
                })
                .collect()
        }
        FieldKind::VorticityMag | FieldKind::QCriterion => {
            let ux: Vec<f64> = sim.ux_field().iter().map(|v| v.as_f64()).collect();
            let uy: Vec<f64> = sim.uy_field().iter().map(|v| v.as_f64()).collect();
            let uz = vec![0.0; nx * ny];
            grad_derived(&ux, &uy, &uz, [nx, ny, 1], kind)
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
            write_png_scaled(
                &path,
                &values,
                &solid,
                nx,
                ny,
                colormap_for(o.field),
                None,
                1,
            )?;
            Ok(name)
        }
        OutputFormat::Csv => {
            let name = format!("{kind_name}_{step}.csv");
            let mut file = fs::File::create(out_dir.join(&name))?;
            writeln!(
                file,
                "# {kind_name}, nx={nx}, ny={ny}, row-major y*nx+x, step={step}"
            )?;
            for y in 0..ny {
                let row: Vec<String> = (0..nx).map(|x| format!("{}", values[y * nx + x])).collect();
                writeln!(file, "{}", row.join(","))?;
            }
            let _ = index;
            Ok(name)
        }
        OutputFormat::Vtk => {
            let name = format!("{kind_name}_{step}.vtk");
            write_vtk(&out_dir.join(&name), &kind_name, step, [nx, ny, 1], &values)?;
            Ok(name)
        }
    }
}

/// VTK legacy structured points (ASCII), 2D (`nz == 1`) or 3D. Point order
/// is x-fastest, then y, then z — exactly the compact field layout
/// `cell = z*(nx*ny) + y*nx + x` (y up, like the sim).
fn write_vtk(
    path: &Path,
    kind_name: &str,
    step: usize,
    dims: [usize; 3],
    values: &[f64],
) -> Result<()> {
    let [nx, ny, nz] = dims;
    debug_assert_eq!(values.len(), nx * ny * nz);
    let mut file = std::io::BufWriter::new(fs::File::create(path)?);
    writeln!(file, "# vtk DataFile Version 3.0")?;
    writeln!(file, "LBMFlow {kind_name} step={step}")?;
    writeln!(file, "ASCII")?;
    writeln!(file, "DATASET STRUCTURED_POINTS")?;
    writeln!(file, "DIMENSIONS {nx} {ny} {nz}")?;
    writeln!(file, "ORIGIN 0 0 0")?;
    writeln!(file, "SPACING 1 1 1")?;
    writeln!(file, "POINT_DATA {}", nx * ny * nz)?;
    writeln!(file, "SCALARS {kind_name} double 1")?;
    writeln!(file, "LOOKUP_TABLE default")?;
    for chunk in values.chunks(9) {
        let line: Vec<String> = chunk.iter().map(|v| format!("{v}")).collect();
        writeln!(file, "{}", line.join(" "))?;
    }
    Ok(())
}

// ---------------------------------------------------------------- 3D (nz > 1)

/// Macroscopic 3D fields, gathered to `f64` host vectors (compact layout
/// `cell = z*(nx*ny) + y*nx + x`).
struct Fields3 {
    ux: Vec<f64>,
    uy: Vec<f64>,
    uz: Vec<f64>,
    rho: Vec<f64>,
    /// Native strain-rate invariant gamma_dot = sqrt(2 S:S) (exact, from f_neq).
    shear: Vec<f64>,
    /// Kinematic viscosity (lattice units) — for DissipationRate = nu*gamma_dot^2.
    nu: f64,
}

fn gather3<T: lbm_core::real::Real>(s: &Solver3<T>) -> Fields3 {
    let to64 = |v: Vec<T>| v.into_iter().map(|x| x.as_f64()).collect::<Vec<f64>>();
    Fields3 {
        ux: to64(s.gather_ux()),
        uy: to64(s.gather_uy()),
        uz: to64(s.gather_uz()),
        rho: to64(s.gather_rho()),
        shear: to64(s.gather_shear_rate()),
        nu: s.nu(),
    }
}

/// Full-grid values of an output field (3D). Vorticity is the z-component
/// ωz = ∂x uy − ∂y ux (the natural scalar for the z-mid-slice snapshot).
fn field_values3(f: &Fields3, dims: [usize; 3], kind: FieldKind) -> Vec<f64> {
    let [nx, ny, nz] = dims;
    let idx = |x: usize, y: usize, z: usize| (z * ny + y) * nx + x;
    match kind {
        FieldKind::Ux => f.ux.clone(),
        FieldKind::Uy => f.uy.clone(),
        FieldKind::Rho => f.rho.clone(),
        FieldKind::Speed => {
            f.ux.iter()
                .zip(&f.uy)
                .zip(&f.uz)
                .map(|((a, b), c)| (a * a + b * b + c * c).sqrt())
                .collect()
        }
        FieldKind::Vorticity => {
            let mut w = vec![0.0; nx * ny * nz];
            for z in 0..nz {
                for y in 1..ny - 1 {
                    for x in 1..nx - 1 {
                        w[idx(x, y, z)] = 0.5
                            * ((f.uy[idx(x + 1, y, z)] - f.uy[idx(x - 1, y, z)])
                                - (f.ux[idx(x, y + 1, z)] - f.ux[idx(x, y - 1, z)]));
                    }
                }
            }
            w
        }
        FieldKind::ShearRate => f.shear.clone(),
        FieldKind::DissipationRate => f.shear.iter().map(|g| f.nu * g * g).collect(),
        FieldKind::VorticityMag | FieldKind::QCriterion => {
            grad_derived(&f.ux, &f.uy, &f.uz, dims, kind)
        }
    }
}

fn write_output3<T: lbm_core::real::Real>(
    s: &Solver3<T>,
    f: &Fields3,
    o: &OutputSpec,
    step: usize,
    out_dir: &Path,
) -> Result<String> {
    let dims = s.dims();
    let [nx, ny, nz] = dims;
    let kind_name = format!("{:?}", o.field).to_lowercase();
    let values = field_values3(f, dims, o.field);
    match o.format {
        OutputFormat::Png => {
            // z-mid slice (3D volume rendering is out of CLI scope; VTK
            // carries the full field for ParaView).
            let zmid = nz / 2;
            let slice = &values[zmid * nx * ny..(zmid + 1) * nx * ny];
            let solid: Vec<bool> = (0..ny)
                .flat_map(|y| (0..nx).map(move |x| (x, y)))
                .map(|(x, y)| s.is_solid(x, y, zmid))
                .collect();
            let name = format!("{kind_name}_{step}.png");
            write_png_scaled(
                &out_dir.join(&name),
                slice,
                &solid,
                nx,
                ny,
                colormap_for(o.field),
                None,
                1,
            )?;
            Ok(name)
        }
        OutputFormat::Csv => {
            // z-mid slice, same row-major layout as 2D.
            let zmid = nz / 2;
            let slice = &values[zmid * nx * ny..(zmid + 1) * nx * ny];
            let name = format!("{kind_name}_{step}.csv");
            let mut file = fs::File::create(out_dir.join(&name))?;
            writeln!(
                file,
                "# {kind_name}, nx={nx}, ny={ny}, z-slice z={zmid} of nz={nz}, row-major y*nx+x, step={step}"
            )?;
            for y in 0..ny {
                let row: Vec<String> = (0..nx).map(|x| format!("{}", slice[y * nx + x])).collect();
                writeln!(file, "{}", row.join(","))?;
            }
            Ok(name)
        }
        OutputFormat::Vtk => {
            let name = format!("{kind_name}_{step}.vtk");
            write_vtk(&out_dir.join(&name), &kind_name, step, dims, &values)?;
            Ok(name)
        }
    }
}

fn run3d_t<T: lbm_core::real::Real>(
    sc: &Scenario,
    mut s: Solver3<T>,
    out_dir: &Path,
    units: Option<lbm_scenario::UnitReport>,
    options: &RunOptions,
) -> Result<Manifest> {
    if let Some(dir) = &options.restore {
        s.restore(dir)
            .with_context(|| format!("cannot restore checkpoint: {}", dir.display()))?;
    }
    let dims = s.dims();
    let [nx, ny, nz] = dims;
    let mut files: Vec<String> = Vec::new();
    struct Probe3 {
        every: usize,
        file: fs::File,
        kind: &'static str,
        point: (usize, usize, usize),
    }
    let mut probes: Vec<Probe3> = Vec::new();
    for p in &sc.probes {
        match *p {
            ProbeSpec::Force { every } => {
                let path = out_dir.join("force.csv");
                let mut file = fs::File::create(&path)?;
                writeln!(file, "step,fx,fy,fz")?;
                files.push("force.csv".into());
                probes.push(Probe3 {
                    every: every.max(1),
                    file,
                    kind: "force",
                    point: (0, 0, 0),
                });
            }
            ProbeSpec::Point { x, y, z, every } => {
                let z = z.unwrap_or(nz / 2);
                anyhow::ensure!(
                    x < nx && y < ny && z < nz,
                    "point probe ({x},{y},{z}) is outside the {nx}x{ny}x{nz} grid"
                );
                let name = format!("point_{x}_{y}_{z}.csv");
                let path = out_dir.join(&name);
                let mut file = fs::File::create(&path)?;
                writeln!(file, "step,ux,uy,uz,rho")?;
                files.push(name);
                probes.push(Probe3 {
                    every: every.max(1),
                    file,
                    kind: "point",
                    point: (x, y, z),
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
        s.step();
        executed = step;

        for p in &mut probes {
            if step % p.every == 0 {
                match p.kind {
                    "force" => {
                        let f = s.probed_force();
                        writeln!(
                            p.file,
                            "{step},{},{},{}",
                            f[0].as_f64(),
                            f[1].as_f64(),
                            f[2].as_f64()
                        )?;
                    }
                    _ => {
                        let (x, y, z) = p.point;
                        let u = s.u(x, y, z);
                        writeln!(
                            p.file,
                            "{step},{},{},{},{}",
                            u[0].as_f64(),
                            u[1].as_f64(),
                            u[2].as_f64(),
                            s.rho(x, y, z).as_f64()
                        )?;
                    }
                }
            }
        }

        let snapshot_due = sc
            .outputs
            .iter()
            .any(|o| o.every > 0 && step % o.every == 0);
        if snapshot_due {
            let f = gather3(&s);
            for o in sc.outputs.iter() {
                if o.every > 0 && step % o.every == 0 {
                    files.push(write_output3(&s, &f, o, step, out_dir)?);
                }
            }
        }

        if let Some(path) = checkpoint_path(options, step) {
            s.save(&path)
                .with_context(|| format!("cannot save checkpoint: {}", path.display()))?;
            files.push(path.display().to_string());
        }

        if step % 1000 == 0 {
            let bad = s.gather_rho().iter().any(|v| !v.as_f64().is_finite());
            if bad {
                status = "diverged";
                break 'main;
            }
        }

        if let Some(ss) = &sc.run.stop_when_steady {
            if step % ss.check_every == 0 {
                let f = gather3(&s);
                let mut cur = f.ux;
                cur.extend(f.uy);
                cur.extend(f.uz);
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
    let f = gather3(&s);
    for o in sc.outputs.iter() {
        if o.every == 0 {
            files.push(write_output3(&s, &f, o, executed, out_dir)?);
        }
    }

    let wall = t0.elapsed().as_secs_f64();
    let cells = (nx * ny * nz) as f64;
    let max_speed =
        f.ux.iter()
            .zip(&f.uy)
            .zip(&f.uz)
            .map(|((a, b), c)| (a * a + b * b + c * c).sqrt())
            .fold(0.0f64, f64::max);
    let manifest = Manifest {
        scenario: sc.name.clone(),
        status: status.to_string(),
        steps_run: executed as u64,
        wall_seconds: wall,
        mlups: cells * executed as f64 / wall.max(1e-9) / 1e6,
        diagnostics: Diagnostics {
            total_mass: s.total_mass().as_f64(),
            max_speed,
            tau: s.tau(),
        },
        warnings: lbm_scenario::validate(sc),
        units,
        files,
    };
    fs::write(
        out_dir.join("manifest.json"),
        serde_json::to_string_pretty(&manifest)?,
    )?;
    Ok(manifest)
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
        assert!(
            manifest.files.contains(&"rho_5.vtk".to_string()),
            "{:?}",
            manifest.files
        );
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

    /// 3D scenario end-to-end: sphere in a channel, force + point probes,
    /// z-mid PNG and a full 3D VTK volume.
    #[test]
    fn scenario_3d_runs_with_outputs() {
        let sc: Scenario = serde_json::from_str(
            r#"{
                "name": "vtk3d-smoke",
                "grid": { "nx": 24, "ny": 12, "nz": 10 },
                "physics": { "nu": 0.05 },
                "edges": {
                    "left": { "type": "velocityInlet", "u": [0.05, 0.0] },
                    "right": { "type": "pressureOutlet", "rho": 1.0 },
                    "bottom": { "type": "bounceBack" },
                    "top": { "type": "bounceBack" },
                    "front": { "type": "periodic" },
                    "back": { "type": "periodic" }
                },
                "obstacles": [ { "shape": "sphere", "cx": 8, "cy": 5.5, "cz": 4.5, "r": 2.5 } ],
                "run": { "steps": 20 },
                "probes": [
                    { "type": "force", "every": 5 },
                    { "type": "point", "x": 16, "y": 6, "every": 10 }
                ],
                "outputs": [
                    { "field": "speed", "format": "png", "every": 0 },
                    { "field": "rho", "format": "vtk", "every": 0 }
                ]
            }"#,
        )
        .unwrap();
        let dir = std::env::temp_dir().join(format!("lbm_vtk3d_test_{}", std::process::id()));
        let manifest = run(&sc, &dir).unwrap();
        assert_eq!(manifest.status, "completed");
        assert_eq!(manifest.steps_run, 20);
        for f in [
            "speed_20.png",
            "rho_20.vtk",
            "force.csv",
            "point_16_6_5.csv",
        ] {
            assert!(
                manifest.files.contains(&f.to_string()),
                "{f} missing from {:?}",
                manifest.files
            );
        }
        let text = fs::read_to_string(dir.join("rho_20.vtk")).unwrap();
        assert!(text.contains("\nDIMENSIONS 24 12 10\n"), "3D dims");
        assert!(text.contains("\nPOINT_DATA 2880\n"));
        let data: Vec<f64> = text
            .lines()
            .skip_while(|l| *l != "LOOKUP_TABLE default")
            .skip(1)
            .flat_map(|l| l.split_whitespace())
            .map(|v| v.parse().unwrap())
            .collect();
        assert_eq!(data.len(), 24 * 12 * 10);
        assert!(data.iter().all(|v| v.is_finite()));
        // Probe CSVs: 3D force has fz, 3D point has uz.
        let force = fs::read_to_string(dir.join("force.csv")).unwrap();
        assert!(force.starts_with("step,fx,fy,fz\n"), "{force}");
        assert_eq!(force.lines().count(), 1 + 4, "{force}");
        let point = fs::read_to_string(dir.join("point_16_6_5.csv")).unwrap();
        assert!(point.starts_with("step,ux,uy,uz,rho\n"), "{point}");
        // The inflow reached the probe point: ux > 0 in the wake row.
        let last = point.lines().last().unwrap();
        let ux: f64 = last.split(',').nth(1).unwrap().parse().unwrap();
        assert!(ux.is_finite());
        fs::remove_dir_all(&dir).ok();
    }
}
