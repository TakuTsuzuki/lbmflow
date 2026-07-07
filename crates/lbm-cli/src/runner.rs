//! Scenario execution: probes, snapshots, divergence/steady detection,
//! machine-readable manifest.

pub use crate::manifest::Manifest;
use crate::manifest::{
    active_models_for_legacy, capability_report, lattice_id, precision_id, scenario_hash,
    ActiveModelTag, BackendId, CollisionProvenance, Diagnostics, LatticeId, Provenance,
    QoiMethodDescriptor, QoiProvenance, MANIFEST_PATH,
};
use crate::output::{select_output_mode, OutputMode};
use crate::render::{write_png_scaled, Colormap};
use anyhow::{Context, Result};
use lbm_core::backend::CpuScalar;
use lbm_core::compat::multiphase::ShanChen;
use lbm_core::compat::prelude::*;
use lbm_core::halo::LocalPeriodic;
use lbm_core::lattice::{D3Q19, D3Q27};
use lbm_core::prelude::{GlobalSpec, Lattice};
use lbm_core::solver::Solver as CoreSolver;
use lbm_scenario::bioprocess::{PhysicsModel, PulseSpec};
use lbm_scenario::{
    BioprocessScenario, FieldKind, OutputFormat, OutputSpec, ProbeSpec, Scenario, Sim3Handle,
    SimHandle, Solver3,
};
use serde::Serialize;
use std::fs;
use std::io::Write as _;
use std::path::{Path, PathBuf};
use std::time::Instant;

fn provenance(
    sc: &Scenario,
    backend: lbm_scenario::BackendChoice,
    lattice: &str,
    storage: lbm_scenario::StorageSpec,
) -> Provenance {
    let collision = match sc.physics.collision {
        lbm_scenario::CollisionSpec::Bgk => CollisionProvenance {
            kind: "bgk".to_string(),
            magic: None,
            omega_shear: None,
        },
        lbm_scenario::CollisionSpec::Trt => CollisionProvenance {
            kind: "trt".to_string(),
            magic: Some(lbm_core::params::CollisionKind::MAGIC_STD),
            omega_shear: None,
        },
        lbm_scenario::CollisionSpec::CentralMoment
        | lbm_scenario::CollisionSpec::DeprecatedCumulantAlias => CollisionProvenance {
            kind: "central_moment".to_string(),
            magic: None,
            omega_shear: Some(lbm_scenario::central_moment_omega_shear(sc.physics.nu)),
        },
    };
    Provenance {
        backend,
        lattice: lattice.to_string(),
        collision,
        precision: sc.physics.precision,
        storage,
    }
}

#[derive(Clone, Debug, Default)]
pub struct RunOptions {
    pub save_every: Option<usize>,
    pub checkpoint_dir: Option<PathBuf>,
    pub restore: Option<PathBuf>,
    pub output_mode: Option<OutputMode>,
}

pub fn run(sc: &Scenario, out_dir: &Path) -> Result<Manifest> {
    run_with_options(sc, out_dir, &RunOptions::default())
}

pub fn prepare_bioprocess_geometry(
    sc: &lbm_scenario::BioprocessScenario,
) -> Result<lbm_core::geometry::StirredTankGeometry> {
    sc.import_stl_geometry().map_err(|err| {
        anyhow::anyhow!(
            "{}",
            serde_json::to_string(&err).unwrap_or_else(|_| err.to_string())
        )
    })
}

pub fn run_bioprocess_single_phase(sc: &BioprocessScenario, out_dir: &Path) -> Result<Manifest> {
    fs::create_dir_all(out_dir)
        .with_context(|| format!("cannot create output directory: {}", out_dir.display()))?;
    if sc.run.grid_nz <= 1 || sc.run.lattice == Some(lbm_scenario::bioprocess::LatticeSpec::D2q9) {
        anyhow::bail!(
            "{}",
            serde_json::to_string_pretty(&lbm_scenario::BioprocessScenarioError {
                message: "single-phase stirred-tank runner requires a 3D lattice".to_string(),
                reason: lbm_scenario::UnsupportedReason::OutOfValidityRange {
                    detail: "run.grid_nz must be > 1 and lattice must be D3Q19 or D3Q27"
                        .to_string(),
                },
            })?
        );
    }
    if !sc
        .physics
        .models
        .iter()
        .any(|m| matches!(m, PhysicsModel::SinglePhase))
    {
        anyhow::bail!("bioprocess runner currently supports physics.kind single_phase");
    }
    let unit_report = sc.compute_unit_report().map_err(|err| {
        anyhow::anyhow!(
            "{}",
            serde_json::to_string_pretty(&err).unwrap_or_else(|_| err.to_string())
        )
    })?;
    let geometry = prepare_bioprocess_geometry(sc)?;
    match (
        sc.run
            .lattice
            .unwrap_or(lbm_scenario::bioprocess::LatticeSpec::D3q19),
        sc.run
            .precision
            .unwrap_or(lbm_scenario::bioprocess::Precision::F64),
    ) {
        (
            lbm_scenario::bioprocess::LatticeSpec::D3q19,
            lbm_scenario::bioprocess::Precision::F64,
        ) => {
            run_bioprocess_single_phase_t::<D3Q19, f64>(sc, unit_report, geometry, "D3Q19", out_dir)
        }
        (
            lbm_scenario::bioprocess::LatticeSpec::D3q19,
            lbm_scenario::bioprocess::Precision::F32,
        ) => {
            run_bioprocess_single_phase_t::<D3Q19, f32>(sc, unit_report, geometry, "D3Q19", out_dir)
        }
        (
            lbm_scenario::bioprocess::LatticeSpec::D3q27,
            lbm_scenario::bioprocess::Precision::F64,
        ) => {
            run_bioprocess_single_phase_t::<D3Q27, f64>(sc, unit_report, geometry, "D3Q27", out_dir)
        }
        (
            lbm_scenario::bioprocess::LatticeSpec::D3q27,
            lbm_scenario::bioprocess::Precision::F32,
        ) => {
            run_bioprocess_single_phase_t::<D3Q27, f32>(sc, unit_report, geometry, "D3Q27", out_dir)
        }
        (lbm_scenario::bioprocess::LatticeSpec::D2q9, _) => unreachable!("D2Q9 rejected above"),
    }
}

pub fn run_with_options(sc: &Scenario, out_dir: &Path, options: &RunOptions) -> Result<Manifest> {
    if select_output_mode(1, options.output_mode) == OutputMode::PerRank {
        anyhow::bail!("--output-mode per-rank requires an MPI run; the current legacy runner has world_size=1");
    }
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
            anyhow::bail!(
                "{}",
                lbm_scenario::gpu_capability_error(sc)
                    .map(|e| e.to_string())
                    .unwrap_or_else(|| {
                        "requested backend \"gpu\" is not wired to the 3D CLI runner yet"
                            .to_string()
                    })
            );
        }
        return run_gpu2d(sc, lbm_scenario::build_gpu2d(sc)?, out_dir, units);
    }
    #[cfg(not(feature = "gpu"))]
    if lbm_scenario::selected_backend(sc) == lbm_scenario::BackendChoice::Gpu {
        anyhow::bail!("{}", lbm_scenario::build_check(sc).unwrap_err());
    }
    if sc.is_3d() {
        return match lbm_scenario::build3d(sc)? {
            Sim3Handle::D3Q19F64(s) => run3d_t(sc, s, "D3Q19", out_dir, units, options),
            Sim3Handle::D3Q19F32(s) => run3d_t(sc, s, "D3Q19", out_dir, units, options),
            Sim3Handle::D3Q27F64(s) => run3d_t(sc, s, "D3Q27", out_dir, units, options),
            Sim3Handle::D3Q27F32(s) => run3d_t(sc, s, "D3Q27", out_dir, units, options),
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
    units: Option<lbm_scenario::LegacyUnitReport>,
) -> Result<Manifest> {
    let t0 = Instant::now();
    let mut files = Vec::new();
    let total = sc.run.steps;
    sim.run(total);
    sim.sync_host();
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
        scenario_hash: scenario_hash(sc)?,
        manifest_path: MANIFEST_PATH.to_string(),
        bioprocess_schema_version: None,
        backend: BackendId::Gpu,
        lattice: LatticeId::D2q9,
        precision: precision_id(sc.physics.precision),
        active_models: active_models_for_legacy(sc),
        qoi_methods: Vec::new(),
        unit_report: None,
        capability_report: capability_report(),
        status: "completed".to_string(),
        steps_run: total as u64,
        wall_seconds: wall,
        mlups: (dims[0] * dims[1]) as f64 * total as f64 / wall.max(1e-9) / 1e6,
        diagnostics: Diagnostics {
            total_mass: sim.total_mass() as f64,
            max_speed,
            tau: 3.0 * sc.physics.nu + 0.5,
            phase_field: None,
        },
        mpi_ranks: Vec::new(),
        provenance: provenance(
            sc,
            lbm_scenario::BackendChoice::Gpu,
            "D2Q9",
            lbm_scenario::requested_storage(sc),
        ),
        warnings: lbm_scenario::validate(sc),
        units,
        files,
    };
    fs::write(
        out_dir.join(MANIFEST_PATH),
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
                Some(MANIFEST_PATH),
            )?;
            Ok(name)
        }
        OutputFormat::Csv => {
            let name = format!("{kind_name}_{step}.csv");
            let mut file = fs::File::create(out_dir.join(&name))?;
            writeln!(file, "# manifest_path={MANIFEST_PATH}")?;
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
    writeln!(file, "# manifest_path={MANIFEST_PATH}")?;
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
    units: Option<lbm_scenario::LegacyUnitReport>,
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
                writeln!(file, "# manifest_path={MANIFEST_PATH}")?;
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
                writeln!(file, "# manifest_path={MANIFEST_PATH}")?;
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
        writeln!(f, "# manifest_path={MANIFEST_PATH}")?;
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
        if mp.is_some() || rotor.is_some() {
            sim.force_field_mut().fill([T::zero(); 2]);
        }
        if let Some(mp) = &mp {
            mp.update_force(&mut sim);
        }
        if let Some(rt) = &mut rotor {
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
            ps.step(sampler, None::<fn([f64; 3]) -> f64>)?;
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
        scenario_hash: scenario_hash(sc)?,
        manifest_path: MANIFEST_PATH.to_string(),
        bioprocess_schema_version: None,
        backend: BackendId::Cpu,
        lattice: LatticeId::D2q9,
        precision: precision_id(sc.physics.precision),
        active_models: active_models_for_legacy(sc),
        qoi_methods: Vec::new(),
        unit_report: None,
        capability_report: capability_report(),
        status: status.to_string(),
        steps_run: executed as u64,
        wall_seconds: wall,
        mlups: cells * executed as f64 / wall.max(1e-9) / 1e6,
        diagnostics: Diagnostics {
            total_mass: sim.total_mass_f64(),
            max_speed,
            tau: sim.tau(),
            phase_field: None,
        },
        mpi_ranks: Vec::new(),
        provenance: provenance(
            sc,
            lbm_scenario::BackendChoice::Cpu,
            "D2Q9",
            lbm_scenario::StorageSpec::F32,
        ),
        warnings: lbm_scenario::validate(sc),
        units,
        files,
    };
    fs::write(
        out_dir.join(MANIFEST_PATH),
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
                Some(MANIFEST_PATH),
            )?;
            Ok(name)
        }
        OutputFormat::Csv => {
            let name = format!("{kind_name}_{step}.csv");
            let mut file = fs::File::create(out_dir.join(&name))?;
            writeln!(
                file,
                "# {kind_name}, nx={nx}, ny={ny}, row-major y*nx+x, step={step}, manifest_path={MANIFEST_PATH}"
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
            write_vtk(
                &out_dir.join(&name),
                &kind_name,
                step,
                [nx, ny, 1],
                &values,
                MANIFEST_PATH,
            )?;
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
    manifest_path: &str,
) -> Result<()> {
    let [nx, ny, nz] = dims;
    debug_assert_eq!(values.len(), nx * ny * nz);
    let mut file = std::io::BufWriter::new(fs::File::create(path)?);
    writeln!(file, "# vtk DataFile Version 3.0")?;
    writeln!(
        file,
        "LBMFlow {kind_name} step={step} manifest_path={manifest_path}"
    )?;
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

#[derive(Clone, Copy, Debug)]
struct TorqueRecord {
    step: u64,
    torque_lu: f64,
    torque_n_m: f64,
    force_lu: [f64; 3],
    force_n: [f64; 3],
}

#[derive(Serialize)]
#[serde(rename_all = "snake_case")]
struct PowerQoiFile {
    provenance: QoiFileProvenance,
    torque_n_m: f64,
    power_w: f64,
    rotational_speed_hz: f64,
    np: f64,
    p_over_v_w_m3: f64,
    nq: Option<f64>,
    skipped: Vec<SkippedQoiFile>,
}

#[derive(Serialize)]
#[serde(rename_all = "snake_case")]
struct StressQoiFile {
    provenance: QoiFileProvenance,
    gamma_dot_1_s: stress_summary::SerdePercentiles,
    viscous_stress_pa: stress_summary::SerdePercentiles,
}

mod stress_summary {
    use lbm_core::stress::PercentileSummary;
    use serde::Serialize;

    #[derive(Serialize)]
    #[serde(rename_all = "snake_case")]
    pub struct SerdePercentiles {
        pub p50: f64,
        pub p90: f64,
        pub p95: f64,
        pub p99: f64,
        pub max: f64,
        pub fraction_above_threshold: Option<f64>,
    }

    impl From<PercentileSummary> for SerdePercentiles {
        fn from(value: PercentileSummary) -> Self {
            Self {
                p50: value.p50,
                p90: value.p90,
                p95: value.p95,
                p99: value.p99,
                max: value.max,
                fraction_above_threshold: value.fraction_above_threshold,
            }
        }
    }
}

#[derive(Serialize)]
#[serde(rename_all = "snake_case")]
struct MixingQoiFile {
    provenance: QoiFileProvenance,
    cv0: Option<f64>,
    t95_s: Option<f64>,
    t99_s: Option<f64>,
    skipped: Option<SkippedQoiFile>,
    compartments: Vec<CompartmentCvFile>,
}

#[derive(Serialize)]
#[serde(rename_all = "snake_case")]
struct CompartmentCvFile {
    name: String,
    cv: Option<f64>,
    cell_count: usize,
}

#[derive(Serialize)]
#[serde(rename_all = "snake_case")]
struct SkippedQoiFile {
    qoi: String,
    reason: String,
}

#[derive(Serialize)]
#[serde(rename_all = "snake_case")]
struct QoiFileProvenance {
    units: String,
    method: String,
    time_window: String,
    averaging_region: String,
    source_fields: Vec<String>,
    validation_tier: String,
}

fn run_bioprocess_single_phase_t<L, T>(
    sc: &BioprocessScenario,
    unit_report: lbm_scenario::UnitReport,
    mut geometry: lbm_core::geometry::StirredTankGeometry,
    lattice: &str,
    out_dir: &Path,
) -> Result<Manifest>
where
    L: Lattice,
    T: lbm_core::real::Real,
{
    let dims = [
        sc.run.grid_nx as usize,
        sc.run.grid_ny as usize,
        sc.run.grid_nz as usize,
    ];
    let mut wall_u_t = Vec::with_capacity(geometry.wall_velocity.len());
    for u in &geometry.wall_velocity {
        wall_u_t.push([
            T::r(u[0] * sc.run.dt_s),
            T::r(u[1] * sc.run.dt_s),
            T::r(u[2] * sc.run.dt_s),
        ]);
    }
    scale_impellers_to_lattice_time(&mut geometry, sc.run.dt_s);
    let solid = geometry.solid.clone();
    let spec = GlobalSpec::<T> {
        dims,
        nu: unit_report.lattice.nu_lu,
        collision: lbm_core::params::CollisionKind::Trt {
            magic: lbm_core::params::CollisionKind::MAGIC_STD,
        },
        periodic: [false, false, false],
        faces: [lbm_core::params::FaceBC::Closed; 6],
        force: [T::zero(); 3],
        sources: Vec::new(),
        face_patches: Vec::new(),
    };
    let solid_t = solid.clone();
    let mut solver: CoreSolver<L, T, CpuScalar, LocalPeriodic> = CoreSolver::try_new(
        &spec,
        &solid_t,
        &wall_u_t,
        [1, 1, 1],
        CpuScalar::default(),
        LocalPeriodic,
    )
    .map_err(|e| anyhow::anyhow!("invalid bioprocess solver spec: {e}"))?;

    let scalar = primary_scalar(sc);
    if let Some((name, diffusivity_lu, initial)) =
        scalar_initial_condition(sc, &unit_report, &solid)?
    {
        solver
            .enable_scalar_distribution(name, diffusivity_lu)
            .map_err(|e| anyhow::anyhow!("cannot enable scalar distribution: {e}"))?;
        let init_t: Vec<T> = initial.into_iter().map(T::r).collect();
        solver
            .set_scalar_concentration(name, &init_t)
            .map_err(|e| anyhow::anyhow!("cannot initialize scalar distribution: {e}"))?;
    }

    let probe_every = sc.outputs.probes_every_n_steps.unwrap_or(1).max(1);
    let field_every = sc
        .outputs
        .fields_every_n_steps
        .unwrap_or(sc.run.steps)
        .max(1);
    let mut files = Vec::new();
    let mut qoi_methods = Vec::new();
    let mut active_models = vec![ActiveModelTag::SinglePhase, ActiveModelTag::RotatingIbm];
    if sc
        .physics
        .models
        .iter()
        .any(|m| matches!(m, PhysicsModel::PassiveScalar { .. }))
    {
        active_models.push(ActiveModelTag::PassiveScalar);
    }
    if sc
        .physics
        .models
        .iter()
        .any(|m| matches!(m, PhysicsModel::Oxygen { .. }))
    {
        active_models.push(ActiveModelTag::Oxygen);
    }

    let mut torque_csv = fs::File::create(out_dir.join("torque_force.csv"))?;
    writeln!(
        torque_csv,
        "step,time_s,torque_lu,torque_n_m,fx_lu,fy_lu,fz_lu,fx_n,fy_n,fz_n"
    )?;
    writeln!(torque_csv, "# manifest_path={MANIFEST_PATH}")?;
    files.push("torque_force.csv".to_string());

    let mut scalar_cv_csv = if scalar.is_some() {
        let mut f = fs::File::create(out_dir.join("scalar_cv.csv"))?;
        writeln!(f, "step,time_s,cv")?;
        writeln!(f, "# manifest_path={MANIFEST_PATH}")?;
        files.push("scalar_cv.csv".to_string());
        Some(f)
    } else {
        None
    };

    let mut torque_records = Vec::new();
    let mut cv_series = Vec::new();
    let t0 = Instant::now();
    let mut status = "completed";
    let mut executed = 0u64;
    for step in 1..=sc.run.steps {
        solver.clear_body_force_field();
        let diagnostics = solver.apply_impeller_marker_sets(
            &geometry.impellers,
            lbm_core::rotating_ibm::DirectForcingConfig::default(),
        );
        solver.step();
        executed = step;
        if step % probe_every == 0 || step == sc.run.steps {
            let rec = torque_record_from_ibm(&diagnostics, &unit_report, sc);
            writeln!(
                torque_csv,
                "{},{},{},{},{},{},{},{},{},{}",
                step,
                step as f64 * sc.run.dt_s,
                rec.torque_lu,
                rec.torque_n_m,
                rec.force_lu[0],
                rec.force_lu[1],
                rec.force_lu[2],
                rec.force_n[0],
                rec.force_n[1],
                rec.force_n[2]
            )?;
            torque_records.push(TorqueRecord { step, ..rec });
        }
        if let (Some((name, _, _)), Some(csv)) = (scalar, scalar_cv_csv.as_mut()) {
            if step % probe_every == 0 || step == sc.run.steps {
                if let Some(c) = solver.gather_scalar_concentration(name) {
                    let values: Vec<f64> = c.into_iter().map(|v| v.as_f64()).collect();
                    let include: Vec<bool> = solid.iter().map(|&s| !s).collect();
                    if let Some(cv) = lbm_core::qoi::scalar_cv(&values, &include) {
                        let t = step as f64 * sc.run.dt_s;
                        writeln!(csv, "{step},{t},{cv}")?;
                        cv_series.push((t, cv));
                    }
                }
            }
        }
        if step % field_every == 0 || step == sc.run.steps {
            files.extend(write_bioprocess_velocity_artifacts(
                &solver,
                step as usize,
                out_dir,
            )?);
        }
        if step % 1000 == 0 && !solver.total_mass_f64().is_finite() {
            status = "diverged";
            break;
        }
    }
    if !files.iter().any(|f| f.ends_with(".png")) {
        files.extend(write_bioprocess_velocity_artifacts(
            &solver,
            executed as usize,
            out_dir,
        )?);
    }

    if sc.qoi.power.is_some() {
        files.extend(write_power_qois(
            sc,
            &unit_report,
            &torque_records,
            out_dir,
        )?);
        qoi_methods.push(QoiMethodDescriptor {
            qoi: "power".to_string(),
            method: "IBM marker drive torque: P = omega*Tq, Np = P/(rho*N^3*D^5)".to_string(),
            input_fields: vec!["ibm_marker_force".to_string()],
            provenance: QoiProvenance::new(
                vec!["ibm_marker_force".to_string()],
                "last_half_of_run",
                "SI",
                "screening",
            ),
        });
    }

    files.extend(write_stress_and_wall_outputs(
        &solver,
        sc,
        &unit_report,
        &solid,
        &geometry.wall_velocity,
        out_dir,
    )?);
    qoi_methods.push(QoiMethodDescriptor {
        qoi: "shear_rate".to_string(),
        method: "central finite differences on velocity moments".to_string(),
        input_fields: vec!["ux".to_string(), "uy".to_string(), "uz".to_string()],
        provenance: QoiProvenance::new(
            vec!["ux".to_string(), "uy".to_string(), "uz".to_string()],
            "final_step",
            "1/s and Pa",
            "screening",
        ),
    });

    if sc.qoi.mixing.is_some() {
        files.extend(write_mixing_qois(
            sc,
            &mut solver,
            &unit_report,
            scalar,
            &solid,
            &cv_series,
            out_dir,
        )?);
        qoi_methods.push(QoiMethodDescriptor {
            qoi: "mixing_time".to_string(),
            method: "scalar coefficient-of-variation threshold".to_string(),
            input_fields: vec!["C:tracer".to_string()],
            provenance: QoiProvenance::new(
                vec!["C:tracer".to_string()],
                "post_pulse_time_series",
                "s",
                "screening",
            ),
        });
    }

    let fields = gather_core3(&solver);
    let max_speed = fields
        .ux
        .iter()
        .zip(&fields.uy)
        .zip(&fields.uz)
        .map(|((a, b), c)| (a * a + b * b + c * c).sqrt())
        .fold(0.0, f64::max);
    let wall = t0.elapsed().as_secs_f64();
    let manifest = Manifest {
        scenario: sc.name.clone(),
        scenario_hash: scenario_hash(sc)?,
        manifest_path: MANIFEST_PATH.to_string(),
        bioprocess_schema_version: Some(sc.version.clone()),
        backend: BackendId::Cpu,
        lattice: lattice_id(lattice),
        precision: match std::mem::size_of::<T>() {
            4 => crate::manifest::PrecisionId::F32,
            _ => crate::manifest::PrecisionId::F64,
        },
        active_models,
        qoi_methods,
        unit_report: Some(unit_report.clone()),
        capability_report: capability_report(),
        status: status.to_string(),
        steps_run: executed,
        wall_seconds: wall,
        mlups: (dims[0] * dims[1] * dims[2]) as f64 * executed as f64 / wall.max(1.0e-9) / 1.0e6,
        mpi_ranks: Vec::new(),
        diagnostics: Diagnostics {
            total_mass: solver.total_mass_f64(),
            max_speed,
            tau: solver.tau(),
            phase_field: None,
        },
        provenance: Provenance {
            backend: lbm_scenario::BackendChoice::Cpu,
            lattice: lattice.to_string(),
            collision: CollisionProvenance {
                kind: "trt".to_string(),
                magic: Some(lbm_core::params::CollisionKind::MAGIC_STD),
                omega_shear: None,
            },
            precision: match std::mem::size_of::<T>() {
                4 => lbm_scenario::Precision::F32,
                _ => lbm_scenario::Precision::F64,
            },
            storage: lbm_scenario::StorageSpec::F32,
        },
        warnings: Vec::new(),
        units: None,
        files,
    };
    fs::write(
        out_dir.join(MANIFEST_PATH),
        serde_json::to_string_pretty(&manifest)?,
    )?;
    Ok(manifest)
}

fn scale_impellers_to_lattice_time(
    geometry: &mut lbm_core::geometry::StirredTankGeometry,
    dt_s: f64,
) {
    for impeller in &mut geometry.impellers {
        for a in 0..3 {
            impeller.omega[a] *= dt_s;
            for u in &mut impeller.wall_velocity {
                u[a] *= dt_s;
            }
        }
    }
}

fn primary_scalar(sc: &BioprocessScenario) -> Option<(&'static str, &PulseSpec, f64)> {
    sc.physics.models.iter().find_map(|model| {
        if let PhysicsModel::PassiveScalar {
            diffusivity_m2_per_s,
            initial_pulse: Some(pulse),
        } = model
        {
            Some(("tracer", pulse, *diffusivity_m2_per_s))
        } else {
            None
        }
    })
}

fn scalar_initial_condition(
    sc: &BioprocessScenario,
    unit_report: &lbm_scenario::UnitReport,
    solid: &[bool],
) -> Result<Option<(&'static str, f64, Vec<f64>)>> {
    let Some((name, pulse, diffusivity_m2_s)) = primary_scalar(sc) else {
        if sc
            .physics
            .models
            .iter()
            .any(|m| matches!(m, PhysicsModel::PassiveScalar { .. }))
        {
            let n = (sc.run.grid_nx * sc.run.grid_ny * sc.run.grid_nz) as usize;
            let mut values = vec![0.0; n];
            for (v, &s) in values.iter_mut().zip(solid) {
                if !s {
                    *v = 1.0;
                }
            }
            let d_lu = diffusivity_lu_from_model(sc, unit_report).ok_or_else(|| {
                anyhow::anyhow!("passive scalar model is missing diffusivity_m2_s")
            })?;
            return Ok(Some(("tracer", d_lu, values)));
        }
        return Ok(None);
    };
    let dims = [
        sc.run.grid_nx as usize,
        sc.run.grid_ny as usize,
        sc.run.grid_nz as usize,
    ];
    let mut values = vec![0.0; dims[0] * dims[1] * dims[2]];
    let dx = unit_report.lattice.dx_m;
    for z in 0..dims[2] {
        for y in 0..dims[1] {
            for x in 0..dims[0] {
                let i = (z * dims[1] + y) * dims[0] + x;
                if solid[i] {
                    continue;
                }
                let p = [
                    (x as f64 + 0.5) * dx,
                    (y as f64 + 0.5) * dx,
                    (z as f64 + 0.5) * dx,
                ];
                let r = ((p[0] - pulse.center_m[0]).powi(2)
                    + (p[1] - pulse.center_m[1]).powi(2)
                    + (p[2] - pulse.center_m[2]).powi(2))
                .sqrt();
                if r <= pulse.radius_m {
                    values[i] = pulse.concentration;
                }
            }
        }
    }
    Ok(Some((
        name,
        diffusivity_m2_s * unit_report.lattice.dt_s / unit_report.lattice.dx_m.powi(2),
        values,
    )))
}

fn diffusivity_lu_from_model(
    sc: &BioprocessScenario,
    unit_report: &lbm_scenario::UnitReport,
) -> Option<f64> {
    sc.physics.models.iter().find_map(|model| {
        if let PhysicsModel::PassiveScalar {
            diffusivity_m2_per_s,
            ..
        } = model
        {
            Some(diffusivity_m2_per_s * unit_report.lattice.dt_s / unit_report.lattice.dx_m.powi(2))
        } else {
            None
        }
    })
}

fn torque_record_from_ibm(
    diagnostics: &[lbm_core::rotating_ibm::IbmDiagnostics],
    unit_report: &lbm_scenario::UnitReport,
    sc: &BioprocessScenario,
) -> TorqueRecord {
    let mut torque_lu = 0.0;
    let mut force_lu = [0.0; 3];
    for diag in diagnostics {
        torque_lu += -diag.torque[2];
        for (a, dst) in force_lu.iter_mut().enumerate() {
            *dst += diag.fluid_force[a];
        }
    }
    let rho = sc.fluids.liquid_density_kg_m3;
    let dx = unit_report.lattice.dx_m;
    let dt = unit_report.lattice.dt_s;
    let force_scale = rho * dx.powi(4) / dt.powi(2);
    let torque_scale = rho * dx.powi(5) / dt.powi(2);
    TorqueRecord {
        step: 0,
        torque_lu,
        torque_n_m: torque_lu * torque_scale,
        force_lu,
        force_n: [
            force_lu[0] * force_scale,
            force_lu[1] * force_scale,
            force_lu[2] * force_scale,
        ],
    }
}

struct CoreFields3 {
    ux: Vec<f64>,
    uy: Vec<f64>,
    uz: Vec<f64>,
}

fn gather_core3<L, T>(solver: &CoreSolver<L, T, CpuScalar, LocalPeriodic>) -> CoreFields3
where
    L: Lattice,
    T: lbm_core::real::Real,
{
    CoreFields3 {
        ux: solver.gather_ux().into_iter().map(|v| v.as_f64()).collect(),
        uy: solver.gather_uy().into_iter().map(|v| v.as_f64()).collect(),
        uz: solver.gather_uz().into_iter().map(|v| v.as_f64()).collect(),
    }
}

fn write_bioprocess_velocity_artifacts<L, T>(
    solver: &CoreSolver<L, T, CpuScalar, LocalPeriodic>,
    step: usize,
    out_dir: &Path,
) -> Result<Vec<String>>
where
    L: Lattice,
    T: lbm_core::real::Real,
{
    let dims = solver.dims();
    let [nx, ny, nz] = dims;
    let fields = gather_core3(solver);
    let speed: Vec<f64> = fields
        .ux
        .iter()
        .zip(&fields.uy)
        .zip(&fields.uz)
        .map(|((a, b), c)| (a * a + b * b + c * c).sqrt())
        .collect();
    let vtk = format!("velocity_speed_{step}.vtk");
    write_vtk(
        &out_dir.join(&vtk),
        "velocity_speed",
        step,
        dims,
        &speed,
        MANIFEST_PATH,
    )?;
    let zmid = nz / 2;
    let png = format!("velocity_speed_{step}.png");
    let solid: Vec<bool> = (0..ny)
        .flat_map(|y| (0..nx).map(move |x| solver.is_solid(x, y, zmid)))
        .collect();
    write_png_scaled(
        &out_dir.join(&png),
        &speed[zmid * nx * ny..(zmid + 1) * nx * ny],
        &solid,
        nx,
        ny,
        Colormap::Viridis,
        None,
        1,
        Some(MANIFEST_PATH),
    )?;
    Ok(vec![vtk, png])
}

fn reactor_power_inputs(sc: &BioprocessScenario) -> Result<(f64, f64, f64)> {
    match &sc.reactor {
        lbm_scenario::bioprocess::ReactorSpec::StirredTank {
            working_volume_m3,
            impellers,
            ..
        } => {
            let impeller = impellers
                .first()
                .ok_or_else(|| anyhow::anyhow!("power QOI requires at least one impeller"))?;
            let diameter = impeller.diameter_m().ok_or_else(|| {
                anyhow::anyhow!("power QOI requires parametric impeller diameter")
            })?;
            Ok((
                *working_volume_m3,
                diameter,
                impeller.rotational_speed_rpm() * std::f64::consts::TAU / 60.0,
            ))
        }
    }
}

fn write_power_qois(
    sc: &BioprocessScenario,
    unit_report: &lbm_scenario::UnitReport,
    records: &[TorqueRecord],
    out_dir: &Path,
) -> Result<Vec<String>> {
    if records.is_empty() {
        return Ok(Vec::new());
    }
    let start = sc.run.steps / 2;
    let window: Vec<_> = records.iter().filter(|r| r.step >= start).collect();
    let used = if window.is_empty() {
        records.iter().collect::<Vec<_>>()
    } else {
        window
    };
    let torque = used.iter().map(|r| r.torque_n_m).sum::<f64>() / used.len() as f64;
    let (volume, diameter, omega) = reactor_power_inputs(sc)?;
    let power = lbm_core::qoi::power_qois(lbm_core::qoi::PowerQoiInput {
        torque_n_m: torque,
        omega_rad_s: omega,
        rho_kg_m3: sc.fluids.liquid_density_kg_m3,
        impeller_diameter_m: diameter,
        working_volume_m3: volume,
        discharge_flow_m3_s: None,
    })
    .map_err(|e| anyhow::anyhow!(e))?;
    let provenance = QoiFileProvenance {
        units: "SI".to_string(),
        method: "IBM marker shaft torque, P = omega*Tq".to_string(),
        time_window: "last_half_of_run".to_string(),
        averaging_region: "impeller_marker_set".to_string(),
        source_fields: vec!["ibm_marker_force".to_string()],
        validation_tier: "screening".to_string(),
    };
    let skipped = vec![SkippedQoiFile {
        qoi: "nq".to_string(),
        reason: "discharge surface is not defined in the M0 scenario schema".to_string(),
    }];
    let json = PowerQoiFile {
        provenance,
        torque_n_m: power.torque_n_m,
        power_w: power.power_w,
        rotational_speed_hz: power.rotational_speed_hz,
        np: power.np,
        p_over_v_w_m3: power.p_over_v_w_m3,
        nq: power.nq,
        skipped,
    };
    fs::write(
        out_dir.join("qoi_power.json"),
        serde_json::to_string_pretty(&json)?,
    )?;
    let mut csv = fs::File::create(out_dir.join("qoi_power.csv"))?;
    writeln!(csv, "qoi,value,units,method,time_window")?;
    writeln!(csv, "# manifest_path={MANIFEST_PATH}")?;
    writeln!(
        csv,
        "torque,{},{},IBM marker shaft torque,last_half_of_run",
        power.torque_n_m, "N*m"
    )?;
    writeln!(
        csv,
        "power,{},{},P=omega*Tq,last_half_of_run",
        power.power_w, "W"
    )?;
    writeln!(
        csv,
        "np,{},dimensionless,Np=P/(rho*N^3*D^5),last_half_of_run",
        power.np
    )?;
    writeln!(
        csv,
        "p_over_v,{},{},P/V,last_half_of_run",
        power.p_over_v_w_m3, "W/m^3"
    )?;
    let _ = unit_report;
    Ok(vec![
        "qoi_power.json".to_string(),
        "qoi_power.csv".to_string(),
    ])
}

fn write_stress_and_wall_outputs<L, T>(
    solver: &CoreSolver<L, T, CpuScalar, LocalPeriodic>,
    sc: &BioprocessScenario,
    unit_report: &lbm_scenario::UnitReport,
    solid: &[bool],
    wall_velocity_si: &[[f64; 3]],
    out_dir: &Path,
) -> Result<Vec<String>>
where
    L: Lattice,
    T: lbm_core::real::Real,
{
    let dims = solver.dims();
    let n = dims[0] * dims[1] * dims[2];
    let fields = gather_core3(solver);
    let velocity_scale = unit_report.lattice.dx_m / unit_report.lattice.dt_s;
    let ux_si: Vec<f64> = fields.ux.iter().map(|v| v * velocity_scale).collect();
    let uy_si: Vec<f64> = fields.uy.iter().map(|v| v * velocity_scale).collect();
    let uz_si: Vec<f64> = fields.uz.iter().map(|v| v * velocity_scale).collect();
    let mu = vec![sc.fluids.liquid_viscosity_pa_s; n];
    let stress = lbm_core::stress::compute_stress_field(
        dims,
        &ux_si,
        &uy_si,
        &uz_si,
        solid,
        unit_report.lattice.dx_m,
        &mu,
    );
    let fluid_gamma: Vec<f64> = stress
        .iter()
        .zip(solid)
        .filter_map(|(s, &is_solid)| (!is_solid).then_some(s.gamma_dot))
        .collect();
    let fluid_tau: Vec<f64> = stress
        .iter()
        .zip(solid)
        .filter_map(|(s, &is_solid)| (!is_solid).then_some(s.viscous_stress_pa))
        .collect();
    let gamma = lbm_core::stress::percentile_summary(&fluid_gamma, Some(0.0))
        .ok_or_else(|| anyhow::anyhow!("stress QOI has no fluid cells"))?;
    let tau = lbm_core::stress::percentile_summary(&fluid_tau, Some(0.0))
        .ok_or_else(|| anyhow::anyhow!("stress QOI has no fluid cells"))?;
    let stress_json = StressQoiFile {
        provenance: QoiFileProvenance {
            units: "1/s for gamma_dot, Pa for viscous_stress".to_string(),
            method: "central finite differences; fraction_above_threshold uses threshold 0 pending BCFD-080 threshold schema".to_string(),
            time_window: "final_step".to_string(),
            averaging_region: "tank_fluid_cells".to_string(),
            source_fields: vec!["ux".to_string(), "uy".to_string(), "uz".to_string()],
            validation_tier: "screening".to_string(),
        },
        gamma_dot_1_s: gamma.into(),
        viscous_stress_pa: tau.into(),
    };
    fs::write(
        out_dir.join("stress_summary.json"),
        serde_json::to_string_pretty(&stress_json)?,
    )?;

    let wall = lbm_core::stress::wall_shear_proxy(
        dims,
        &ux_si,
        &uy_si,
        &uz_si,
        solid,
        wall_velocity_si,
        unit_report.lattice.dx_m,
        sc.fluids.liquid_density_kg_m3,
        sc.fluids.liquid_viscosity_pa_s,
    );
    let mut csv = fs::File::create(out_dir.join("wall_shear.csv"))?;
    writeln!(
        csv,
        "cell_index,tau_w_pa_proxy,y_plus,u_parallel_m_s,y_m,nx,ny,nz,label"
    )?;
    writeln!(csv, "# manifest_path={MANIFEST_PATH}")?;
    for w in wall {
        writeln!(
            csv,
            "{},{},{},{},{},{},{},{},proxy",
            w.cell_index,
            w.tau_w_pa,
            w.y_plus.unwrap_or(f64::NAN),
            w.u_parallel_m_s,
            w.y_m,
            w.normal[0],
            w.normal[1],
            w.normal[2]
        )?;
    }
    Ok(vec![
        "stress_summary.json".to_string(),
        "wall_shear.csv".to_string(),
    ])
}

fn write_mixing_qois<L, T>(
    sc: &BioprocessScenario,
    solver: &mut CoreSolver<L, T, CpuScalar, LocalPeriodic>,
    unit_report: &lbm_scenario::UnitReport,
    scalar: Option<(&'static str, &PulseSpec, f64)>,
    solid: &[bool],
    cv_series: &[(f64, f64)],
    out_dir: &Path,
) -> Result<Vec<String>>
where
    L: Lattice,
    T: lbm_core::real::Real,
{
    let dims = [
        sc.run.grid_nx as usize,
        sc.run.grid_ny as usize,
        sc.run.grid_nz as usize,
    ];
    let (result, skipped, compartments) = if let Some((name, _, _)) = scalar {
        let final_scalar = solver
            .gather_scalar_concentration(name)
            .ok_or_else(|| anyhow::anyhow!("scalar distribution {name} missing"))?;
        let values: Vec<f64> = final_scalar.into_iter().map(|v| v.as_f64()).collect();
        let impeller_center_z = first_impeller_center_z_lu(sc, unit_report);
        let compartments = lbm_core::qoi::compartment_cv(dims, &values, solid, impeller_center_z)
            .into_iter()
            .map(|c| CompartmentCvFile {
                name: c.name,
                cv: c.cv,
                cell_count: c.cell_count,
            })
            .collect();
        (
            lbm_core::qoi::mixing_time_from_cv(cv_series),
            None,
            compartments,
        )
    } else {
        (
            None,
            Some(SkippedQoiFile {
                qoi: "mixing_time".to_string(),
                reason: "scalar pulse is not configured".to_string(),
            }),
            Vec::new(),
        )
    };
    let json = MixingQoiFile {
        provenance: QoiFileProvenance {
            units: "s".to_string(),
            method: "CV(t)=std(C)/mean(C); t95/t99 use 0.05*CV0 and 0.01*CV0".to_string(),
            time_window: "post_pulse_time_series".to_string(),
            averaging_region: "tank_fluid_cells".to_string(),
            source_fields: vec!["C:tracer".to_string()],
            validation_tier: "screening".to_string(),
        },
        cv0: result.as_ref().map(|r| r.cv0),
        t95_s: result.as_ref().and_then(|r| r.t95_s),
        t99_s: result.as_ref().and_then(|r| r.t99_s),
        skipped,
        compartments,
    };
    fs::write(
        out_dir.join("mixing_time.json"),
        serde_json::to_string_pretty(&json)?,
    )?;
    Ok(vec!["mixing_time.json".to_string()])
}

fn first_impeller_center_z_lu(
    sc: &BioprocessScenario,
    unit_report: &lbm_scenario::UnitReport,
) -> Option<f64> {
    match &sc.reactor {
        lbm_scenario::bioprocess::ReactorSpec::StirredTank { impellers, .. } => {
            impellers.first().and_then(|impeller| match impeller {
                lbm_scenario::bioprocess::ImpellerSpec::Rushton {
                    clearance_from_bottom_m,
                    ..
                }
                | lbm_scenario::bioprocess::ImpellerSpec::PitchedBlade {
                    clearance_from_bottom_m,
                    ..
                }
                | lbm_scenario::bioprocess::ImpellerSpec::Marine {
                    clearance_from_bottom_m,
                    ..
                } => Some(*clearance_from_bottom_m / unit_report.lattice.dx_m),
                lbm_scenario::bioprocess::ImpellerSpec::CustomMarkerSet { .. } => None,
            })
        }
    }
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

fn gather3<L, T>(s: &Solver3<L, T>) -> Fields3
where
    L: Lattice,
    T: lbm_core::real::Real,
{
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

fn write_output3<L, T>(
    s: &Solver3<L, T>,
    f: &Fields3,
    o: &OutputSpec,
    step: usize,
    out_dir: &Path,
) -> Result<String>
where
    L: Lattice,
    T: lbm_core::real::Real,
{
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
                Some(MANIFEST_PATH),
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
                "# {kind_name}, nx={nx}, ny={ny}, z-slice z={zmid} of nz={nz}, row-major y*nx+x, step={step}, manifest_path={MANIFEST_PATH}"
            )?;
            for y in 0..ny {
                let row: Vec<String> = (0..nx).map(|x| format!("{}", slice[y * nx + x])).collect();
                writeln!(file, "{}", row.join(","))?;
            }
            Ok(name)
        }
        OutputFormat::Vtk => {
            let name = format!("{kind_name}_{step}.vtk");
            write_vtk(
                &out_dir.join(&name),
                &kind_name,
                step,
                dims,
                &values,
                MANIFEST_PATH,
            )?;
            Ok(name)
        }
    }
}

fn run3d_t<L, T>(
    sc: &Scenario,
    mut s: Solver3<L, T>,
    lattice: &str,
    out_dir: &Path,
    units: Option<lbm_scenario::LegacyUnitReport>,
    options: &RunOptions,
) -> Result<Manifest>
where
    L: Lattice,
    T: lbm_core::real::Real,
{
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
                writeln!(file, "# manifest_path={MANIFEST_PATH}")?;
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
                writeln!(file, "# manifest_path={MANIFEST_PATH}")?;
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
        scenario_hash: scenario_hash(sc)?,
        manifest_path: MANIFEST_PATH.to_string(),
        bioprocess_schema_version: None,
        backend: BackendId::Cpu,
        lattice: lattice_id(lattice),
        precision: precision_id(sc.physics.precision),
        active_models: active_models_for_legacy(sc),
        qoi_methods: Vec::new(),
        unit_report: None,
        capability_report: capability_report(),
        status: status.to_string(),
        steps_run: executed as u64,
        wall_seconds: wall,
        mlups: cells * executed as f64 / wall.max(1e-9) / 1e6,
        diagnostics: Diagnostics {
            total_mass: s.total_mass_f64(),
            max_speed,
            tau: s.tau(),
            phase_field: None,
        },
        mpi_ranks: Vec::new(),
        provenance: provenance(
            sc,
            lbm_scenario::BackendChoice::Cpu,
            lattice,
            lbm_scenario::StorageSpec::F32,
        ),
        warnings: lbm_scenario::validate(sc),
        units,
        files,
    };
    fs::write(
        out_dir.join(MANIFEST_PATH),
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
        assert!(force.contains("manifest_path=manifest.json"), "{force}");
        assert_eq!(force.lines().count(), 1 + 1 + 4, "{force}");
        let point = fs::read_to_string(dir.join("point_16_6_5.csv")).unwrap();
        assert!(point.starts_with("step,ux,uy,uz,rho\n"), "{point}");
        assert!(point.contains("manifest_path=manifest.json"), "{point}");
        // The inflow reached the probe point: ux > 0 in the wake row.
        let last = point.lines().last().unwrap();
        let ux: f64 = last.split(',').nth(1).unwrap().parse().unwrap();
        assert!(ux.is_finite());
        fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn manifest_snapshot_for_legacy_scenario_matches_baseline() {
        let sc: Scenario = serde_json::from_str(
            r#"{
                "name": "legacy-manifest-snapshot",
                "grid": { "nx": 10, "ny": 8 },
                "physics": { "nu": 0.05, "precision": "f64" },
                "edges": {
                    "left": { "type": "periodic" }, "right": { "type": "periodic" },
                    "bottom": { "type": "bounceBack" }, "top": { "type": "bounceBack" }
                },
                "run": { "steps": 1 }
            }"#,
        )
        .unwrap();
        let dir =
            std::env::temp_dir().join(format!("lbm_legacy_manifest_test_{}", std::process::id()));
        let manifest = run(&sc, &dir).unwrap();
        let value = serde_json::to_value(&manifest).unwrap();
        assert_eq!(value["scenario"], "legacy-manifest-snapshot");
        assert_eq!(value["manifestPath"], "manifest.json");
        assert!(value["scenarioHash"]
            .as_str()
            .unwrap()
            .starts_with("sha256:"));
        assert_eq!(value["bioprocessSchemaVersion"], serde_json::Value::Null);
        assert_eq!(value["backend"], "cpu");
        assert_eq!(value["lattice"], "d2q9");
        assert_eq!(value["precision"], "f64");
        assert_eq!(value["activeModels"], serde_json::json!(["single_phase"]));
        assert_eq!(value["qoiMethods"], serde_json::json!([]));
        assert!(value.get("unitReport").is_none());
        assert!(value["capabilityReport"].as_array().unwrap().len() >= 8);
        fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn output_files_reference_manifest_path() {
        let sc: Scenario = serde_json::from_str(
            r#"{
                "name": "manifest-path-outputs",
                "grid": { "nx": 12, "ny": 8 },
                "physics": { "nu": 0.05 },
                "edges": {
                    "left": { "type": "periodic" }, "right": { "type": "periodic" },
                    "bottom": { "type": "bounceBack" }, "top": { "type": "bounceBack" }
                },
                "run": { "steps": 1 },
                "outputs": [
                    { "field": "rho", "format": "csv", "every": 0 },
                    { "field": "speed", "format": "vtk", "every": 0 },
                    { "field": "rho", "format": "png", "every": 0 }
                ]
            }"#,
        )
        .unwrap();
        let dir =
            std::env::temp_dir().join(format!("lbm_manifest_path_test_{}", std::process::id()));
        let manifest = run(&sc, &dir).unwrap();
        assert_eq!(manifest.manifest_path, "manifest.json");
        let csv = fs::read_to_string(dir.join("rho_1.csv")).unwrap();
        assert!(csv.contains("manifest_path=manifest.json"), "{csv}");
        let vtk = fs::read_to_string(dir.join("speed_1.vtk")).unwrap();
        assert!(vtk.contains("manifest_path=manifest.json"), "{vtk}");

        let decoder = png::Decoder::new(fs::File::open(dir.join("rho_1.png")).unwrap());
        let reader = decoder.read_info().unwrap();
        let info = reader.info();
        let found = info
            .uncompressed_latin1_text
            .iter()
            .any(|chunk| chunk.keyword == "manifest_path" && chunk.text == "manifest.json");
        assert!(found, "PNG manifest_path text chunk missing");
        fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn manifest_records_actual_provenance() {
        let sc: Scenario = serde_json::from_str(
            r#"{
                "name": "manifest-provenance",
                "grid": { "nx": 10, "ny": 8, "nz": 6 },
                "physics": {
                    "nu": 0.02,
                    "collision": { "type": "central_moment" },
                    "precision": "f64"
                },
                "compute": { "backend": "cpu" },
                "edges": {
                    "left": { "type": "periodic" }, "right": { "type": "periodic" },
                    "bottom": { "type": "bounceBack" }, "top": { "type": "bounceBack" },
                    "front": { "type": "bounceBack" }, "back": { "type": "bounceBack" }
                },
                "run": { "steps": 2 }
            }"#,
        )
        .unwrap();
        let dir = std::env::temp_dir().join(format!(
            "lbm_manifest_provenance_test_{}",
            std::process::id()
        ));
        let manifest = run(&sc, &dir).unwrap();
        assert_eq!(
            manifest.provenance.backend,
            lbm_scenario::BackendChoice::Cpu
        );
        assert_eq!(manifest.provenance.lattice, "D3Q19");
        assert_eq!(manifest.provenance.collision.kind, "central_moment");
        assert_eq!(
            manifest.provenance.collision.omega_shear,
            Some(lbm_scenario::central_moment_omega_shear(sc.physics.nu))
        );
        assert_eq!(manifest.provenance.precision, lbm_scenario::Precision::F64);
        assert_eq!(manifest.provenance.storage, lbm_scenario::StorageSpec::F32);

        let text = fs::read_to_string(dir.join("manifest.json")).unwrap();
        assert!(text.contains("\"provenance\""), "{text}");
        assert!(text.contains("\"lattice\": \"D3Q19\""), "{text}");
        assert!(text.contains("\"kind\": \"central_moment\""), "{text}");
        fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn manifest_records_d3q27_lattice_provenance() {
        let sc: Scenario = serde_json::from_str(
            r#"{
                "name": "manifest-d3q27",
                "grid": { "nx": 10, "ny": 8, "nz": 6, "lattice": "d3q27" },
                "physics": { "nu": 0.05, "precision": "f64" },
                "compute": { "backend": "cpu" },
                "edges": {
                    "left": { "type": "velocityInlet", "u": [0.02, 0.0] },
                    "right": { "type": "pressureOutlet", "rho": 1.0 },
                    "bottom": { "type": "bounceBack" }, "top": { "type": "bounceBack" },
                    "front": { "type": "bounceBack" }, "back": { "type": "bounceBack" }
                },
                "run": { "steps": 3 },
                "outputs": [ { "field": "speed", "format": "vtk", "every": 0 } ]
            }"#,
        )
        .unwrap();
        let dir =
            std::env::temp_dir().join(format!("lbm_manifest_d3q27_test_{}", std::process::id()));
        let manifest = run(&sc, &dir).unwrap();
        assert_eq!(manifest.status, "completed");
        assert_eq!(manifest.provenance.lattice, "D3Q27");
        assert_eq!(
            manifest.provenance.backend,
            lbm_scenario::BackendChoice::Cpu
        );
        assert!(manifest.files.contains(&"speed_3.vtk".to_string()));

        let text = fs::read_to_string(dir.join("manifest.json")).unwrap();
        assert!(text.contains("\"lattice\": \"D3Q27\""), "{text}");
        fs::remove_dir_all(&dir).ok();
    }
}
