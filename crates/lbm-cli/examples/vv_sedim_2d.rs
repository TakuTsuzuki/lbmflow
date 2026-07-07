//! Axis 9.4 sedimentation-basin visual experiment.
//!
//! RUN-NOW capability exercise: 2D compat channel flow with gravity and
//! one-way CR-3 Schiller-Naumann particles. The CR-3 floor coordinate is the
//! particle `z` coordinate, so this example maps the physical basin coordinates
//! as particle `(x, dummy, z=y)` while sampling the 2D fluid at `(x, y=z)`.

use anyhow::{bail, Context, Result};
use lbm_core::compat::prelude::*;
use lbm_core::particles::{sample_grid, DepositEvent, Particle, ParticleSet, Sample};
use std::fs;
use std::fs::File;
use std::io::BufWriter;
use std::path::{Path, PathBuf};

const NX: usize = 128;
const NY: usize = 64;
const STEPS: usize = 10_000;
const SNAPSHOTS: [usize; 3] = [0, 5_000, 10_000];
const N_PARTICLES: usize = 500;
const U_IN: f64 = 0.02;
const NU: f64 = 1.0 / 6.0;
const RHO_F: f64 = 1.0;
const RHO_P: f64 = 2.0;
const D_PARTICLE: f64 = 1.5;
const G_Y: f64 = -5.0e-5;
const FLOOR_Z: f64 = 0.5;
const TOP_Z: f64 = (NY - 2) as f64;
const X_INJECT: f64 = 1.0;
const OUT_DIR: &str = "out/vv_sedim_2d";
const HIST_BINS: usize = 16;

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
    fs::create_dir_all(&out_dir).with_context(|| format!("create {}", out_dir.display()))?;

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
    let full_height = TOP_Z - FLOOR_Z;
    let full_height_xs = U_IN * full_height / v_stokes;
    let mean_deposit_x = mean_deposition_x(&deposits);
    let deposited = deposits.len();
    let suspended = particles.particles.len();

    assert_eq!(
        deposited + suspended,
        N_PARTICLES,
        "mass conservation failed: deposited {deposited} + suspended {suspended} != {N_PARTICLES}"
    );
    assert_monotone_nonincreasing(&hist)?;
    assert_censored_settling_anchor(&deposits, v_stokes)?;

    println!("VV_SEDIM_2D Axis 9.4 sedimentation basin RUN-NOW");
    println!(
        "SETUP nx={NX} ny={NY} steps={STEPS} u_in={U_IN:.6e} nu={NU:.6e} gravity=[0,{G_Y:.6e}] particles={N_PARTICLES} d={D_PARTICLE:.6e} rho_p={RHO_P:.6e} rho_f={RHO_F:.6e} floor_z={FLOOR_Z:.6e}"
    );
    println!(
        "ANCHOR full_height_H={full_height:.6e} v_stokes={v_stokes:.6e} x_s_full_height={full_height_xs:.6e}"
    );
    match mean_deposit_x {
        Some(x) => println!(
            "MEASURED mean_deposition_x={x:.6e} n_deposited={deposited} n_suspended={suspended}"
        ),
        None => println!("MEASURED mean_deposition_x=none n_deposited=0 n_suspended={suspended}"),
    }
    println!("HIST bins={HIST_BINS} counts={hist:?}");
    println!(
        "CENSORED_NOTE full-height settling length exceeds the 128-cell basin; the executable assertion uses each deposited particle's inlet height because only near-floor particles can deposit in {STEPS} steps with the requested g."
    );
    println!("ARTIFACTS");
    for path in artifacts {
        println!("  {}", path.display());
    }
    Ok(())
}

fn build_sim() -> Result<Simulation<f64>> {
    let mut sim = SimConfig {
        nx: NX,
        ny: NY,
        nu: NU,
        edges: Edges {
            left: EdgeBC::VelocityInlet { u: [U_IN, 0.0] },
            right: EdgeBC::PressureOutlet { rho: 1.0 },
            bottom: EdgeBC::BounceBack,
            top: EdgeBC::BounceBack,
        },
        collision: Collision::Trt {
            magic: Collision::MAGIC_STD,
        },
        ..Default::default()
    }
    .build()
    .context("build 2D compat sedimentation channel")?;
    sim.set_gravity([0.0, G_Y]);
    Ok(sim)
}

fn build_particles() -> ParticleSet {
    let span = TOP_Z - FLOOR_Z;
    let mut particles = Vec::with_capacity(N_PARTICLES);
    for i in 0..N_PARTICLES {
        let frac = i as f64 / (N_PARTICLES - 1) as f64;
        let z = FLOOR_Z + 1.0e-6 + frac * (span - 1.0e-6);
        particles.push(Particle {
            pos: [X_INJECT, 0.0, z],
            vel: [U_IN, 0.0, 0.0],
            d: D_PARTICLE,
            rho_p: RHO_P,
            exposure: z,
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

fn deposition_histogram(deposits: &[DepositEvent]) -> Vec<usize> {
    let mut hist = vec![0usize; HIST_BINS];
    for event in deposits {
        let bin = ((event.pos[0] / NX as f64) * HIST_BINS as f64).floor() as isize;
        let bin = bin.clamp(0, HIST_BINS as isize - 1) as usize;
        hist[bin] += 1;
    }
    hist
}

fn assert_monotone_nonincreasing(hist: &[usize]) -> Result<()> {
    for (i, pair) in hist.windows(2).enumerate() {
        if pair[1] > pair[0] {
            bail!(
                "deposition histogram is not monotone decreasing at bins {}->{}: {} -> {}",
                i,
                i + 1,
                pair[0],
                pair[1]
            );
        }
    }
    Ok(())
}

fn assert_censored_settling_anchor(deposits: &[DepositEvent], v_stokes: f64) -> Result<()> {
    if deposits.is_empty() {
        bail!("no deposited particles; deposition-map behavior is unobservable");
    }

    let mut rel_sum = 0.0;
    for event in deposits {
        let inlet_height = event.particle.exposure - FLOOR_Z;
        let expected_x = X_INJECT + U_IN * inlet_height / v_stokes;
        let ratio = event.pos[0] / expected_x.max(1.0e-30);
        rel_sum += ratio;
        if !(1.0 / 3.0..=3.0).contains(&ratio) {
            bail!(
                "deposition x outside factor-3 censored Stokes anchor: measured={:.6e} expected={:.6e} ratio={:.6e} initial_height={:.6e}",
                event.pos[0],
                expected_x,
                ratio,
                inlet_height
            );
        }
    }
    println!(
        "CENSORED_ANCHOR mean_x_over_stokes_expected={:.6e}",
        rel_sum / deposits.len() as f64
    );
    Ok(())
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
