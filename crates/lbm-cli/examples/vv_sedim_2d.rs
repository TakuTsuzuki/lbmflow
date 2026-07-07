//! Axis 9.4 sedimentation-basin visual experiment, rev 3.
//!
//! RUN-NOW capability exercise: closed-box quiescent settling with one-way
//! CR-3 Schiller-Naumann particles. The CR-3 floor coordinate is the particle
//! `z` coordinate, so this example maps the physical basin coordinates as
//! particle `(x, dummy, z=y)` while sampling the 2D fluid at `(x, y=z)`.

use anyhow::{Context, Result};
use lbm_core::compat::prelude::*;
use lbm_core::particles::{sample_grid, DepositEvent, Particle, ParticleSet, Sample};
use std::fs;
use std::fs::File;
use std::io::BufWriter;
use std::path::{Path, PathBuf};

const NX: usize = 128;
const NY: usize = 64;
const STEPS: usize = 60_000;
const SNAPSHOTS: [usize; 3] = [0, 30_000, 60_000];
const N_PARTICLES: usize = 500;
const NU: f64 = 1.0 / 6.0;
const RHO_F: f64 = 1.0;
const RHO_P: f64 = 2.0;
// Rev 3 keeps the in-domain Stokes particle parameters from the rev 2 attempt
// and removes the crossflow/outlet bookkeeping artifact. With no fluid
// crossflow, a seed at z=60 reaches the z=1 deposition plane in about
// (60 - 1) / 1.2e-3 = 4.92e4 terminal-speed time steps; 60k steps includes
// the finite response-time transient while keeping Re_p safely below 1.
// Re_p = v_s d / nu = 1.2e-3 * 6 / (1/6) = 0.0432, safely in the Stokes regime.
const D_PARTICLE: f64 = 6.0;
const G_Y: f64 = -1.0e-4;
// `step_depositing` records floor crossings before solid-contact handling.
// With the 2D compat wall rim, `sample_grid` marks the containing lower node as
// solid for contact tests, so particles entering z < 1 contact the bottom rim
// before they can cross the half-way wall coordinate z=0.5. Use the first
// fluid-row coordinate as the deposition counting plane for this visual run.
const FLOOR_Z: f64 = 1.0;
const SEED_Z: f64 = 60.0;
const SEED_X_MIN: f64 = 10.0;
const SEED_X_MAX: f64 = 118.0;
const SEED_X_MEAN: f64 = 0.5 * (SEED_X_MIN + SEED_X_MAX);
const OUT_DIR: &str = "out/vv_sedim_2d";
const HIST_BINS: usize = 16;
const MAX_MEAN_DEPOSITION_X_ERROR: f64 = 5.0;
const MAX_LATERAL_SCATTER_STD: f64 = 5.0;

const VIRIDIS: [[u8; 3]; 16] = [
    [68, 1, 84],
    [72, 26, 108],
    [71, 47, 125],
    [65, 68, 135],
    [57, 86, 140],
    [49, 104, 142],
    [42, 120, 142],
    [35, 136, 142],
    [31, 152, 139],
    [34, 168, 132],
    [53, 183, 121],
    [84, 197, 104],
    [122, 209, 81],
    [165, 219, 54],
    [210, 226, 27],
    [253, 231, 37],
];

const INFERNO: [[u8; 3]; 16] = [
    [0, 0, 4],
    [12, 8, 38],
    [36, 12, 79],
    [66, 10, 104],
    [93, 18, 110],
    [120, 28, 109],
    [147, 38, 103],
    [174, 48, 92],
    [199, 62, 76],
    [220, 81, 57],
    [237, 105, 37],
    [246, 133, 17],
    [251, 163, 12],
    [249, 195, 41],
    [240, 226, 96],
    [252, 255, 164],
];

fn main() -> Result<()> {
    let out_dir = PathBuf::from(OUT_DIR);
    prepare_out_dir(&out_dir)?;

    let mut sim = build_sim()?;
    let mut particles = build_particles();
    let mut deposits = Vec::<DepositEvent>::new();
    let mut artifacts = Vec::<PathBuf>::new();

    write_density_png(&sim, &out_dir, 0, &mut artifacts)?;
    for step in 1..=STEPS {
        sim.step();
        step_particles(&mut particles, &sim, &mut deposits)
            .with_context(|| format!("advance particles at step {step}"))?;
        if SNAPSHOTS.contains(&step) {
            write_density_png(&sim, &out_dir, step, &mut artifacts)?;
        }
    }

    let hist = deposition_histogram(&deposits);
    let map_path = out_dir.join("deposition_map.png");
    write_deposition_png(&map_path, &hist)?;
    artifacts.push(map_path);

    let v_stokes = stokes_settling_speed();
    let re_particle = v_stokes * D_PARTICLE / NU;
    let settling_height = SEED_Z - FLOOR_Z;
    let terminal_settling_steps = settling_height / v_stokes;
    let mean_deposit_x = mean_deposition_x(&deposits);
    let deposit_x_std = std_deposition_x(&deposits, mean_deposit_x);
    let lateral_scatter_std = std_lateral_deposition_error(&deposits);
    let deposited = deposits.len();
    let suspended = particles.particles.len();
    let deposition_fraction = deposited as f64 / N_PARTICLES as f64;
    let x_mean = mean_deposit_x.context("no deposited particles; mean deposit x is undefined")?;
    let x_std = deposit_x_std.context("no deposited particles; deposit x std is undefined")?;
    let lateral_std =
        lateral_scatter_std.context("no deposited particles; lateral scatter std is undefined")?;

    println!("VV_SEDIM_2D Axis 9.4 closed-basin quiescent settling rev 3");
    println!(
        "SETUP nx={NX} ny={NY} steps={STEPS} edges=all_bounce_back fluid_force=[0,0] particle_gravity=[0,0,{G_Y:.6e}] particles={N_PARTICLES} d={D_PARTICLE:.6e} rho_p={RHO_P:.6e} rho_f={RHO_F:.6e} seed_line_x=[{SEED_X_MIN:.6e},{SEED_X_MAX:.6e}] seed_z={SEED_Z:.6e} floor_z={FLOOR_Z:.6e} nu={NU:.6e}"
    );
    println!(
        "ANCHOR settling_height={settling_height:.6e} v_stokes={v_stokes:.6e} re_particle={re_particle:.6e} terminal_settling_steps={terminal_settling_steps:.6e} seed_x_mean={SEED_X_MEAN:.6e}"
    );
    println!(
        "MEASURED deposition_fraction={deposition_fraction:.6e} mean_deposition_x={x_mean:.6e} deposit_x_std={x_std:.6e} lateral_scatter_std={lateral_std:.6e} n_deposited={deposited} n_suspended={suspended}"
    );
    println!(
        "MASS_CHECK deposited_plus_suspended={} expected={N_PARTICLES}",
        deposited + suspended
    );
    println!("HIST bins={HIST_BINS} counts={hist:?}");
    println!("ARTIFACTS");
    for path in &artifacts {
        println!("  {}", path.display());
    }

    assert_eq!(
        deposited + suspended,
        N_PARTICLES,
        "mass conservation failed: deposited {deposited} + suspended {suspended} != {N_PARTICLES}"
    );
    assert!(
        deposition_fraction == 1.0,
        "deposition fraction is not exact: {deposition_fraction:.6e} != 1.0"
    );
    let mean_x_error = (x_mean - SEED_X_MEAN).abs();
    assert!(
        mean_x_error <= MAX_MEAN_DEPOSITION_X_ERROR,
        "mean deposition x drift too large: mean={x_mean:.6e} seed_mean={SEED_X_MEAN:.6e} abs_err={mean_x_error:.6e} > {MAX_MEAN_DEPOSITION_X_ERROR:.6e}"
    );
    assert!(
        lateral_std < MAX_LATERAL_SCATTER_STD,
        "lateral deposition scatter too large: {lateral_std:.6e} >= {MAX_LATERAL_SCATTER_STD:.6e}"
    );
    Ok(())
}

fn prepare_out_dir(out_dir: &Path) -> Result<()> {
    fs::create_dir_all(out_dir).with_context(|| format!("create {}", out_dir.display()))?;
    for entry in fs::read_dir(out_dir).with_context(|| format!("read {}", out_dir.display()))? {
        let entry = entry.with_context(|| format!("read entry in {}", out_dir.display()))?;
        let name = entry.file_name();
        let name = name.to_string_lossy();
        if name == "deposition_map.png" || (name.starts_with("density_") && name.ends_with(".png"))
        {
            fs::remove_file(entry.path())
                .with_context(|| format!("remove stale {}", entry.path().display()))?;
        }
    }
    Ok(())
}

fn build_sim() -> Result<Simulation<f64>> {
    let sim = SimConfig {
        nx: NX,
        ny: NY,
        nu: NU,
        edges: Edges {
            left: EdgeBC::BounceBack,
            right: EdgeBC::BounceBack,
            bottom: EdgeBC::BounceBack,
            top: EdgeBC::BounceBack,
        },
        collision: Collision::Trt {
            magic: Collision::MAGIC_STD,
        },
        ..Default::default()
    }
    .build()
    .context("build 2D compat closed sedimentation basin")?;
    Ok(sim)
}

fn build_particles() -> ParticleSet {
    let mut particles = Vec::with_capacity(N_PARTICLES);
    for i in 0..N_PARTICLES {
        let frac = i as f64 / (N_PARTICLES - 1) as f64;
        let x = SEED_X_MIN + frac * (SEED_X_MAX - SEED_X_MIN);
        particles.push(Particle {
            pos: [x, 0.0, SEED_Z],
            vel: [0.0; 3],
            d: D_PARTICLE,
            rho_p: RHO_P,
            exposure: x,
        });
    }
    ParticleSet::new(particles, RHO_F, NU, [0.0, 0.0, G_Y])
}

fn step_particles(
    particles: &mut ParticleSet,
    sim: &Simulation<f64>,
    deposits: &mut Vec<DepositEvent>,
) -> Result<()> {
    let sample = |p: [f64; 3]| sample_basin(sim, p);
    particles
        .step_depositing(sample, None::<fn([f64; 3]) -> f64>, FLOOR_Z, deposits)
        .map_err(anyhow::Error::new)
}

fn sample_basin(sim: &Simulation<f64>, p: [f64; 3]) -> Sample {
    sample_grid(p, [NX, 1, NY], |x, _, z| {
        ([sim.ux(x, z), 0.0, sim.uy(x, z)], sim.is_solid(x, z))
    })
}

fn stokes_settling_speed() -> f64 {
    ((RHO_P / RHO_F) - 1.0) * G_Y.abs() * D_PARTICLE * D_PARTICLE / (18.0 * NU)
}

fn mean_deposition_x(deposits: &[DepositEvent]) -> Option<f64> {
    if deposits.is_empty() {
        return None;
    }
    Some(deposits.iter().map(|e| e.pos[0]).sum::<f64>() / deposits.len() as f64)
}

fn std_deposition_x(deposits: &[DepositEvent], mean: Option<f64>) -> Option<f64> {
    let mean = mean?;
    Some(
        (deposits
            .iter()
            .map(|e| (e.pos[0] - mean).powi(2))
            .sum::<f64>()
            / deposits.len() as f64)
            .sqrt(),
    )
}

fn std_lateral_deposition_error(deposits: &[DepositEvent]) -> Option<f64> {
    if deposits.is_empty() {
        return None;
    }
    let mean = deposits
        .iter()
        .map(|e| e.pos[0] - e.particle.exposure)
        .sum::<f64>()
        / deposits.len() as f64;
    Some(
        (deposits
            .iter()
            .map(|e| {
                let dx = e.pos[0] - e.particle.exposure;
                (dx - mean).powi(2)
            })
            .sum::<f64>()
            / deposits.len() as f64)
            .sqrt(),
    )
}

fn deposition_histogram(deposits: &[DepositEvent]) -> Vec<usize> {
    let mut hist = vec![0usize; HIST_BINS];
    for event in deposits {
        let bin = ((event.pos[0] / NX as f64) * HIST_BINS as f64).floor() as isize;
        let bin = bin.clamp(0, HIST_BINS as isize - 1) as usize;
        hist[bin] += 1;
    }
    hist
}

fn density_field(sim: &Simulation<f64>) -> Vec<f64> {
    let mut field = vec![0.0; NX * NY];
    for y in 0..NY {
        for x in 0..NX {
            field[y * NX + x] = sim.rho(x, y) - 1.0;
        }
    }
    field
}

fn solid_field(sim: &Simulation<f64>) -> Vec<bool> {
    let mut solid = vec![false; NX * NY];
    for y in 0..NY {
        for x in 0..NX {
            solid[y * NX + x] = sim.is_solid(x, y);
        }
    }
    solid
}

fn write_density_png(
    sim: &Simulation<f64>,
    out_dir: &Path,
    step: usize,
    artifacts: &mut Vec<PathBuf>,
) -> Result<()> {
    let path = out_dir.join(format!("density_{step:05}.png"));
    write_png(
        &path,
        &density_field(sim),
        &solid_field(sim),
        NX,
        NY,
        1.0e-2,
        &VIRIDIS,
        4,
    )?;
    artifacts.push(path);
    Ok(())
}

fn write_deposition_png(path: &Path, hist: &[usize]) -> Result<()> {
    let w = hist.len();
    let h = 24usize;
    let max_count = hist.iter().copied().max().unwrap_or(0).max(1) as f64;
    let mut field = vec![0.0; w * h];
    let solid = vec![false; w * h];
    for (x, &count) in hist.iter().enumerate() {
        let filled = ((count as f64 / max_count) * h as f64).round() as usize;
        for y in 0..filled.min(h) {
            field[y * w + x] = count as f64;
        }
    }
    write_png(path, &field, &solid, w, h, max_count, &INFERNO, 12)
}

fn write_png(
    path: &Path,
    field: &[f64],
    solid: &[bool],
    w: usize,
    h: usize,
    vmax: f64,
    anchors: &[[u8; 3]],
    scale: usize,
) -> Result<()> {
    let (ow, oh) = (w * scale, h * scale);
    let mut buf = vec![0u8; ow * oh * 3];
    for oy in 0..oh {
        let y = oy / scale;
        for ox in 0..ow {
            let x = ox / scale;
            let i = y * w + x;
            let rgb = if solid[i] {
                [90u8, 94, 100]
            } else {
                lut(anchors, field[i] / vmax.max(1.0e-30))
            };
            let px = ((oh - 1 - oy) * ow + ox) * 3;
            buf[px..px + 3].copy_from_slice(&rgb);
        }
    }
    let file = File::create(path).with_context(|| format!("create {}", path.display()))?;
    let mut enc = png::Encoder::new(BufWriter::new(file), ow as u32, oh as u32);
    enc.set_color(png::ColorType::Rgb);
    enc.set_depth(png::BitDepth::Eight);
    let mut writer = enc.write_header().context("write png header")?;
    writer.write_image_data(&buf).context("write png data")?;
    Ok(())
}

fn lut(anchors: &[[u8; 3]], t: f64) -> [u8; 3] {
    let t = t.clamp(0.0, 1.0) * (anchors.len() - 1) as f64;
    let i = (t as usize).min(anchors.len() - 2);
    let f = t - i as f64;
    let (a, b) = (anchors[i], anchors[i + 1]);
    [
        (a[0] as f64 + (b[0] as f64 - a[0] as f64) * f) as u8,
        (a[1] as f64 + (b[1] as f64 - a[1] as f64) * f) as u8,
        (a[2] as f64 + (b[2] as f64 - a[2] as f64) * f) as u8,
    ]
}
