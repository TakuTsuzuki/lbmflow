//! Axis 9.7 lateral-oscillation sloshing V&V run.
//!
//! This is deliberately labeled a Shan-Chen low-density-ratio analog. It is not
//! an air-water free-surface validation: the interface is diffuse, the initial
//! density ratio is only 2.0 / 0.15 = 13.3, and the T11 SC coexistence densities
//! are rho_l ~= 1.888 and rho_v ~= 0.1194. The real free-surface case belongs to
//! the MF-gamma phase-field track.

use anyhow::{bail, ensure, Context, Result};
use lbm_core::compat::multiphase::ShanChen;
use lbm_core::compat::prelude::*;
use std::f64::consts::PI;
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

const NX: usize = 128;
const NY: usize = 96;
const H: f64 = 48.0;
const L: f64 = 128.0;
const STEPS: usize = 20_000;
const SNAPSHOT_STEPS: [usize; 3] = [0, 10_000, 20_000];
const SAMPLE_EVERY: usize = 20;
const INTERFACE_X_MARGIN: usize = 8;

const G_SC: f64 = -5.0;
const NU_TAU_1: f64 = 1.0 / 6.0;
const RHO_INIT_L: f64 = 2.0;
const RHO_INIT_V: f64 = 0.15;
const RHO_T11_L: f64 = 1.888;
const RHO_T11_V: f64 = 0.1194;

const G_Y: f64 = -2.0e-5;
const DRIVE_A: f64 = 1.0e-6;
const MASS_TOL_REL: f64 = 1.0e-6;

#[derive(Clone, Debug)]
struct Sample {
    step: usize,
    modal_eta: f64,
}

#[derive(Clone, Debug)]
struct RunResult {
    label: &'static str,
    omega: f64,
    amplitude: f64,
    mass_rel_drift: f64,
    artifacts: Vec<PathBuf>,
}

fn main() -> Result<()> {
    let out_root = PathBuf::from("out/vv_slosh_2d");
    fs::create_dir_all(&out_root).with_context(|| format!("create {}", out_root.display()))?;

    // Hydrostatic estimate used only to choose the drive frequencies:
    // shallow-water first mode in a closed tank has
    // omega_0^2 = g k tanh(k h), k = pi / L, so
    // omega_0 = sqrt(pi * g_eff * tanh(pi*h/L) / L).
    // The vertical gravity is applied as a per-mass acceleration g_y through
    // `set_gravity`, therefore g_eff = |g_y|.
    let g_eff = -G_Y;
    let omega0 = (PI * g_eff * (PI * H / L).tanh() / L).sqrt();

    println!("Axis 9.7 lateral-oscillation sloshing");
    println!("SCMP: G={G_SC}, tau=1, psi=1-exp(-rho), rho_l/rho_v init = 13.3");
    println!(
        "T11 coexistence anchors: rho_l={RHO_T11_L}, rho_v={RHO_T11_V}, ratio={:.2}",
        RHO_T11_L / RHO_T11_V
    );
    println!("hydrostatic g_y={G_Y:.3e}, g_eff={g_eff:.3e}, omega0={omega0:.6e}");
    println!(
        "HONESTY: density ratio ~13-16 is NOT air-water 1000; this is an unsteady 2D low-ratio diffuse-interface analog, not a sharp free surface. MF-gamma phase-field is required for the real free surface."
    );

    let cases = [
        ("omega_0p5", 0.5 * omega0),
        ("omega_1p0", omega0),
        ("omega_2p0", 2.0 * omega0),
    ];
    let mut results = Vec::new();
    for (label, omega) in cases {
        results.push(run_case(&out_root, label, omega)?);
    }

    let amp_half = results[0].amplitude;
    let amp_res = results[1].amplitude;
    let amp_double = results[2].amplitude;
    ensure!(
        amp_res > amp_half && amp_res > amp_double,
        "resonance anchor failed: amp(omega0)={amp_res:.6e}, amp(0.5omega0)={amp_half:.6e}, amp(2omega0)={amp_double:.6e}"
    );

    println!();
    println!("Summary");
    println!("label, omega, amplitude, mass_rel_drift");
    for r in &results {
        println!(
            "{}, {:.9e}, {:.9e}, {:.9e}",
            r.label, r.omega, r.amplitude, r.mass_rel_drift
        );
    }
    println!("RESONANCE anchor PASS: amp(omega0) > amp(0.5omega0) and amp(omega0) > amp(2omega0)");
    println!("mass conservation PASS: each relative drift <= {MASS_TOL_REL:.1e}");

    println!();
    println!("Artifacts");
    for r in &results {
        for path in &r.artifacts {
            println!("{}", path.display());
        }
    }

    println!();
    println!("Behavior-validity review");
    println!(
        "Pattern: the median-density interface responds most strongly near the shallow-water first-mode estimate."
    );
    println!(
        "Mechanism: lateral oscillation drives the k=pi/L antisymmetric interface mode against hydrostatic restoring gravity."
    );
    println!(
        "Resolved vs closure: Guo forcing, closed bounce-back walls, and SC cohesion are active; the free surface itself is a diffuse Shan-Chen closure, not a validated sharp-interface model."
    );
    println!(
        "Artifacts checked: density PNGs and interface modal time-series were emitted for every frequency; wall/contact-line artifacts require visual review of the listed PNGs."
    );
    println!("Verdict: CLOSURE-DRIVEN low-ratio analog; routing none if anchors pass.");

    Ok(())
}

fn run_case(out_root: &Path, label: &'static str, omega: f64) -> Result<RunResult> {
    let case_dir = out_root.join(label);
    fs::create_dir_all(&case_dir).with_context(|| format!("create {}", case_dir.display()))?;

    let mut sim: Simulation<f64> = SimConfig {
        nx: NX,
        ny: NY,
        nu: NU_TAU_1,
        edges: Edges {
            left: EdgeBC::BounceBack,
            right: EdgeBC::BounceBack,
            bottom: EdgeBC::BounceBack,
            top: EdgeBC::BounceBack,
        },
        ..Default::default()
    }
    .build()
    .context("build closed-box SC sloshing simulation")?;
    sim.init_with(|_, y| {
        if y < NY / 2 {
            (RHO_INIT_L, 0.0, 0.0)
        } else {
            (RHO_INIT_V, 0.0, 0.0)
        }
    });
    sim.set_gravity([0.0, G_Y]);

    let sc = ShanChen::new(G_SC);
    let m0 = sim.total_mass_f64();
    let mut artifacts = Vec::new();
    let mut samples = Vec::new();

    write_density_snapshot(&case_dir, label, 0, &sim, &mut artifacts)?;
    samples.push(Sample {
        step: 0,
        modal_eta: interface_mode(&sim)?,
    });

    for step in 0..STEPS {
        sc.update_force(&mut sim);

        // Non-inertial-frame derivation: if the container displacement is
        // X_b(t), the tank frame adds acceleration -X_b''(t). A sinusoidal
        // acceleration a_x(t)=A sin(omega t) appears in the momentum equation
        // as a uniform horizontal body-force density. SC also owns this force
        // field, so we add the drive after `update_force`; a literal
        // `force_field_mut().fill([A*sin(...), 0])` here would erase cohesion.
        let fx = DRIVE_A * (omega * step as f64).sin();
        for f in sim.force_field_mut() {
            f[0] += fx;
        }

        sim.step();
        let done_step = step + 1;

        if done_step % SAMPLE_EVERY == 0 {
            samples.push(Sample {
                step: done_step,
                modal_eta: interface_mode(&sim)?,
            });
        }
        if SNAPSHOT_STEPS.contains(&done_step) {
            write_density_snapshot(&case_dir, label, done_step, &sim, &mut artifacts)?;
        }
    }

    let mass_rel_drift = relative_abs(sim.total_mass_f64() - m0, m0);
    ensure!(
        mass_rel_drift <= MASS_TOL_REL,
        "{label} mass drift {mass_rel_drift:.6e} exceeds {MASS_TOL_REL:.1e}"
    );
    ensure_finite(&sim, label)?;

    let amplitude = modal_amplitude(&samples)?;
    let ascii_path = case_dir.join(format!("{label}_interface_mode.txt"));
    write_ascii_plot(&ascii_path, label, omega, amplitude, &samples)?;
    artifacts.push(ascii_path);

    println!(
        "{label}: omega={omega:.9e}, amplitude={amplitude:.9e}, mass_rel_drift={mass_rel_drift:.9e}"
    );

    Ok(RunResult {
        label,
        omega,
        amplitude,
        mass_rel_drift,
        artifacts,
    })
}

fn interface_mode(sim: &Simulation<f64>) -> Result<f64> {
    let threshold = 0.5 * (RHO_T11_L + RHO_T11_V);
    let mut numerator = 0.0;
    let mut denominator = 0.0;
    let mut count = 0usize;

    for x in INTERFACE_X_MARGIN..NX - INTERFACE_X_MARGIN {
        let Some(y_int) = interface_y(sim, x, threshold) else {
            continue;
        };
        let xi = x as f64 + 0.5;
        let phi = (PI * xi / L).cos();
        numerator += y_int * phi;
        denominator += phi * phi;
        count += 1;
    }

    let span = NX - 2 * INTERFACE_X_MARGIN;
    if count < span * 4 / 5 || denominator <= 0.0 {
        bail!("interface projection has only {count} usable contour columns");
    }
    Ok(numerator / denominator)
}

fn interface_y(sim: &Simulation<f64>, x: usize, threshold: f64) -> Option<f64> {
    let mut prev_y = 1usize;
    let mut prev = sim.rho(x, prev_y) - threshold;
    for y in 2..NY - 1 {
        let cur = sim.rho(x, y) - threshold;
        if (prev >= 0.0 && cur <= 0.0) || (prev <= 0.0 && cur >= 0.0) {
            let denom = prev - cur;
            if denom.abs() < 1.0e-30 {
                return Some((prev_y as f64 + y as f64) * 0.5);
            }
            let frac = prev / denom;
            return Some(prev_y as f64 + frac * (y - prev_y) as f64);
        }
        prev_y = y;
        prev = cur;
    }
    None
}

fn modal_amplitude(samples: &[Sample]) -> Result<f64> {
    let start = STEPS / 2;
    let mut lo = f64::INFINITY;
    let mut hi = f64::NEG_INFINITY;
    let mut n = 0usize;
    for s in samples {
        if s.step >= start {
            if s.modal_eta < lo {
                lo = s.modal_eta;
            }
            if s.modal_eta > hi {
                hi = s.modal_eta;
            }
            n += 1;
        }
    }
    if n < 2 || !lo.is_finite() || !hi.is_finite() {
        bail!("not enough finite second-half modal samples");
    }
    Ok(0.5 * (hi - lo).abs())
}

fn write_density_snapshot(
    case_dir: &Path,
    label: &str,
    step: usize,
    sim: &Simulation<f64>,
    artifacts: &mut Vec<PathBuf>,
) -> Result<()> {
    let path = case_dir.join(format!("{label}_density_{step:05}.png"));
    let field: Vec<f64> = sim.rho_field().to_vec();
    write_density_png(&path, &field, sim.solid_field(), NX, NY, 5)?;
    artifacts.push(path);
    Ok(())
}

fn write_density_png(
    path: &Path,
    field: &[f64],
    solid: &[bool],
    nx: usize,
    ny: usize,
    scale: usize,
) -> Result<()> {
    let mut lo = f64::INFINITY;
    let mut hi = f64::NEG_INFINITY;
    for (v, s) in field.iter().zip(solid) {
        if !*s && v.is_finite() {
            if *v < lo {
                lo = *v;
            }
            if *v > hi {
                hi = *v;
            }
        }
    }
    if !lo.is_finite() || !hi.is_finite() || hi <= lo {
        lo = 0.0;
        hi = 1.0;
    }

    let sc = if scale == 0 { 1 } else { scale };
    let ow = nx * sc;
    let oh = ny * sc;
    let mut buf = vec![0u8; ow * oh * 3];
    for oy in 0..oh {
        let y = oy / sc;
        for ox in 0..ow {
            let x = ox / sc;
            let i = y * nx + x;
            let p = ((oh - 1 - oy) * ow + ox) * 3;
            let rgb = if solid[i] {
                [90, 94, 100]
            } else {
                density_rgb((field[i] - lo) / (hi - lo))
            };
            buf[p..p + 3].copy_from_slice(&rgb);
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

fn density_rgb(t: f64) -> [u8; 3] {
    let u = if t.is_finite() {
        if t < 0.0 {
            0.0
        } else if t > 1.0 {
            1.0
        } else {
            t
        }
    } else {
        0.0
    };
    let r = (30.0 + 220.0 * u) as u8;
    let g = (55.0 + 175.0 * (1.0 - (2.0 * u - 1.0).abs())) as u8;
    let b = (180.0 - 150.0 * u) as u8;
    [r, g, b]
}

fn write_ascii_plot(
    path: &Path,
    label: &str,
    omega: f64,
    amplitude: f64,
    samples: &[Sample],
) -> Result<()> {
    let width = 80usize;
    let height = 21usize;
    let mut lo = f64::INFINITY;
    let mut hi = f64::NEG_INFINITY;
    for s in samples {
        if s.modal_eta.is_finite() {
            if s.modal_eta < lo {
                lo = s.modal_eta;
            }
            if s.modal_eta > hi {
                hi = s.modal_eta;
            }
        }
    }
    if !lo.is_finite() || !hi.is_finite() || hi <= lo {
        lo = -1.0;
        hi = 1.0;
    }

    let mut grid = vec![vec![b' '; width]; height];
    for col in 0..width {
        let idx = col * samples.len() / width;
        let s = &samples[idx];
        let y_norm = (s.modal_eta - lo) / (hi - lo);
        let mut row = ((1.0 - y_norm) * (height - 1) as f64).round() as isize;
        if row < 0 {
            row = 0;
        }
        if row >= height as isize {
            row = height as isize - 1;
        }
        grid[row as usize][col] = b'*';
    }

    let mut file =
        BufWriter::new(File::create(path).with_context(|| format!("create {}", path.display()))?);
    writeln!(file, "label={label}")?;
    writeln!(file, "omega={omega:.12e}")?;
    writeln!(
        file,
        "amplitude_second_half_half_peak_to_peak={amplitude:.12e}"
    )?;
    writeln!(file, "step,modal_eta")?;
    for s in samples {
        writeln!(file, "{}, {:.12e}", s.step, s.modal_eta)?;
    }
    writeln!(file)?;
    writeln!(file, "ASCII plot of modal_eta over the full run:")?;
    for row in grid {
        writeln!(file, "{}", String::from_utf8(row).expect("ascii plot row"))?;
    }
    Ok(())
}

fn ensure_finite(sim: &Simulation<f64>, label: &str) -> Result<()> {
    for (i, rho) in sim.rho_field().iter().enumerate() {
        ensure!(rho.is_finite(), "{label}: non-finite rho at cell {i}");
    }
    for (i, ux) in sim.ux_field().iter().enumerate() {
        ensure!(ux.is_finite(), "{label}: non-finite ux at cell {i}");
    }
    for (i, uy) in sim.uy_field().iter().enumerate() {
        ensure!(uy.is_finite(), "{label}: non-finite uy at cell {i}");
    }
    Ok(())
}

fn relative_abs(delta: f64, reference: f64) -> f64 {
    if reference == 0.0 {
        delta.abs()
    } else {
        (delta / reference).abs()
    }
}
