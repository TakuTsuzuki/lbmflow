mod particles;
mod protocol;
mod readout;
mod reservoir;

use lbm_core::prelude::*;
use particles::{deposit_batch, make_reservoir_particles};
use protocol::ProtocolInput;
use reservoir::extract_by_depth;
use std::path::PathBuf;

type Sim3 = Solver<D3Q19, f64, CpuSimd, LocalPeriodic>;

fn main() -> anyhow::Result<()> {
    let path = std::env::args()
        .nth(1)
        .map(PathBuf::from)
        .ok_or_else(|| anyhow::anyhow!("usage: dispersed_seeding <protocol.json>"))?;
    let input = ProtocolInput::from_path(&path)?;
    validate_protocol(&input)?;
    let regime = input.regime()?;
    println!(
        "REGIME dx={:.3e} m dt={:.3e} s nu*={:.5} u_jet*={:.4} Ma={:.3} Re_jet={:.1} St={:.3e} Fr={:.3} tau={:.4} nozzle_d={:.4e} m tau_p={:.4e} s settling={:.4e} m/s",
        regime.dx,
        regime.dt,
        regime.nu_lattice,
        regime.u_jet_lattice,
        regime.ma,
        regime.re_jet,
        regime.st,
        regime.fr,
        regime.tau,
        regime.nozzle_d_m,
        regime.particle_tau_s,
        regime.settling_m_s
    );

    let mut reservoir_sim = build_reservoir_sim(&input, regime.nu_lattice);
    run_until_quasi_steady(&mut reservoir_sim, 0.0, 1)?;
    let reservoir_velocity = gather_velocity(&reservoir_sim);

    let reservoir_particles = make_reservoir_particles(&input);
    let extraction = extract_by_depth(&input, &reservoir_particles);
    let batch = extraction.batch;
    let n_extracted = batch.len();

    let mut tray_sim = build_tray_sim(&input, &regime, true);
    run_until_quasi_steady(&mut tray_sim, regime.u_jet_lattice, 512)?;
    let (suspended, deposits) = deposit_batch(&input, &regime, batch, &mut tray_sim)?;
    let tray_velocity = gather_velocity_si(&tray_sim, &regime);
    if deposits.len() + suspended.particles.len() != n_extracted {
        anyhow::bail!(
            "particle count ledger failed: deposited {} + suspended {} != extracted {}",
            deposits.len(),
            suspended.particles.len(),
            n_extracted
        );
    }
    let metrics = readout::write_outputs(
        &input,
        &regime,
        &deposits,
        suspended.particles.len(),
        n_extracted,
        &reservoir_velocity,
        &tray_velocity,
    )?;
    println!(
        "RESULT CV={:.6} max_over_mean={:.6} empty_bin_fraction={:.6} n_deposited={} n_suspended={} n_extracted={} out={}",
        metrics.cv,
        metrics.max_over_mean,
        metrics.empty_bin_fraction,
        metrics.n_deposited,
        metrics.n_suspended,
        metrics.n_extracted,
        input.output.dir
    );
    if !extraction.histogram.is_empty() {
        let summary = extraction
            .histogram
            .iter()
            .map(|(d, n)| format!("{:.0}um:{n}", d * 1.0e6))
            .collect::<Vec<_>>()
            .join(", ");
        println!("EXTRACTED_DIAMETER_HIST {summary}");
    }
    Ok(())
}

fn validate_protocol(input: &ProtocolInput) -> anyhow::Result<()> {
    validate_grid_counts(input)?;
    if !(0.0..=1.0).contains(&input.reservoir.initial_conc) {
        anyhow::bail!("reservoir.initial_conc must be in [0, 1]");
    }
    input
        .op("withdraw")
        .ok_or_else(|| anyhow::anyhow!("protocol requires a withdraw operation"))?;
    input
        .op("eject")
        .ok_or_else(|| anyhow::anyhow!("protocol requires an eject operation"))?;
    for op in &input.protocol {
        if op.op == "settle" && op.duration_s.is_none() {
            anyhow::bail!("settle.duration_s is required");
        }
        if op.op == "withdraw" {
            let depth = op
                .depth_frac
                .ok_or_else(|| anyhow::anyhow!("withdraw.depth_frac is required"))?;
            if !(0.0..=1.0).contains(&depth) {
                anyhow::bail!("withdraw.depth_frac must be in [0, 1]");
            }
            let volume = op
                .volume_frac
                .ok_or_else(|| anyhow::anyhow!("withdraw.volume_frac is required"))?;
            if !(0.0..=1.0).contains(&volume) {
                anyhow::bail!("withdraw.volume_frac must be in [0, 1]");
            }
            if op.rate_ul_s.is_none() {
                anyhow::bail!("withdraw.rate_uLs is required");
            }
        }
        if op.op == "agitate" && op.pattern.as_deref() != Some("translational") {
            anyhow::bail!(
                "unsupported agitate pattern {:?}; only translational is implemented",
                op.pattern
            );
        }
        if op.op == "eject" {
            let nozzle_d = op
                .nozzle_diameter_m
                .ok_or_else(|| anyhow::anyhow!("eject.nozzle_diameter_m is required"))?;
            if nozzle_d <= 0.0 {
                anyhow::bail!("eject.nozzle_diameter_m must be positive");
            }
            if op.rate_ul_s.is_none() {
                anyhow::bail!("eject.rate_uLs is required");
            }
            let h = op
                .height_m
                .ok_or_else(|| anyhow::anyhow!("eject.height_m is required"))?;
            if !(0.0..=input.target.height_m).contains(&h) {
                anyhow::bail!("eject.height_m must be inside the target height");
            }
            let points = op
                .points_xy_frac
                .as_ref()
                .ok_or_else(|| anyhow::anyhow!("eject.points_xy_frac is required"))?;
            let radius = 0.5 * nozzle_d;
            for &pt in points {
                let cx = pt[0] * input.target.width_m;
                let cy = pt[1] * input.target.depth_m;
                if cx - radius < 0.0
                    || cx + radius > input.target.width_m
                    || cy - radius < 0.0
                    || cy + radius > input.target.depth_m
                {
                    anyhow::bail!(
                        "eject nozzle disk centered at ({:.3}, {:.3}) must lie inside the target",
                        pt[0],
                        pt[1]
                    );
                }
            }
        }
    }
    Ok(())
}

fn validate_grid_counts(input: &ProtocolInput) -> anyhow::Result<()> {
    let res_dx = input.grid.dx_m;
    let tray_dx = input.grid.tray_dx_m.unwrap_or(res_dx);
    if res_dx <= 0.0 {
        anyhow::bail!("grid.dx_m must be positive");
    }
    if tray_dx <= 0.0 {
        anyhow::bail!("grid.tray_dx_m must be positive");
    }
    let expected = [
        (
            "grid.res_nx",
            input.grid.res_nx,
            input.reservoir.width_m / res_dx,
        ),
        (
            "grid.res_ny",
            input.grid.res_ny,
            input.reservoir.width_m / res_dx,
        ),
        (
            "grid.res_nz",
            input.grid.res_nz,
            input.reservoir.height_m / res_dx,
        ),
        (
            "grid.tray_nx",
            input.grid.tray_nx,
            input.target.width_m / tray_dx,
        ),
        (
            "grid.tray_ny",
            input.grid.tray_ny,
            input.target.depth_m / tray_dx,
        ),
        (
            "grid.tray_nz",
            input.grid.tray_nz,
            input.target.height_m / tray_dx,
        ),
    ];
    for (name, actual, expected_f) in expected {
        let expected_n = expected_f.round();
        if (actual as f64 - expected_n).abs() > 1.0 {
            anyhow::bail!(
                "{name}={actual} disagrees with SI/dx={expected_f:.3} by more than 1 cell"
            );
        }
    }
    Ok(())
}

fn build_reservoir_sim(input: &ProtocolInput, nu: f64) -> Sim3 {
    let dims = [input.grid.res_nx, input.grid.res_ny, input.grid.res_nz];
    let mut walls = WallSpec::<f64>::default();
    for face in Face::ALL {
        walls.is_wall[face.index()] = true;
    }
    let spec = GlobalSpec::<f64> {
        dims,
        nu,
        periodic: [false, false, false],
        collision: CollisionKind::Trt {
            magic: CollisionKind::MAGIC_STD,
        },
        ..Default::default()
    };
    let (solid, wall_u) = build_wall_rims(3, dims, &walls);
    Sim3::new(
        &spec,
        &solid,
        &wall_u,
        [1, 1, 1],
        CpuSimd::default(),
        LocalPeriodic,
    )
}

fn build_tray_sim(input: &ProtocolInput, regime: &protocol::Regime, jets_on: bool) -> Sim3 {
    let dims = [input.grid.tray_nx, input.grid.tray_ny, input.grid.tray_nz];
    let mut walls = WallSpec::<f64>::default();
    for face in [Face::XNeg, Face::XPos, Face::YNeg, Face::YPos, Face::ZNeg] {
        walls.is_wall[face.index()] = true;
    }
    let faces = [FaceBC::Closed; 6];
    let face_patches = if jets_on {
        nozzle_face_patches(input, regime)
    } else {
        Vec::new()
    };
    let spec = GlobalSpec::<f64> {
        dims,
        nu: regime.nu_lattice,
        periodic: [false, false, false],
        faces,
        face_patches,
        collision: CollisionKind::Trt {
            magic: CollisionKind::MAGIC_STD,
        },
        ..Default::default()
    };
    let (solid, wall_u) = build_wall_rims(3, dims, &walls);
    Sim3::new(
        &spec,
        &solid,
        &wall_u,
        [1, 1, 1],
        CpuSimd::default(),
        LocalPeriodic,
    )
}

fn nozzle_face_patches(input: &ProtocolInput, regime: &protocol::Regime) -> Vec<FacePatch<f64>> {
    let points = input
        .op("eject")
        .and_then(|op| op.points_xy_frac.clone())
        .unwrap_or_else(|| vec![[0.5, 0.5]]);
    let nozzle_radius_m = 0.5 * regime.nozzle_d_m;
    let nx = input.grid.tray_nx;
    let ny = input.grid.tray_ny;
    let mut patches = Vec::new();
    for y in 0..ny {
        let mut x = 0usize;
        while x < nx {
            while x < nx && !inside_any_nozzle(input, regime, &points, nozzle_radius_m, x, y) {
                x += 1;
            }
            if x == nx {
                break;
            }
            let lo_x = x;
            while x + 1 < nx && inside_any_nozzle(input, regime, &points, nozzle_radius_m, x + 1, y)
            {
                x += 1;
            }
            patches.push(FacePatch {
                face: Face::ZPos.index(),
                lo: [lo_x, y],
                hi: [x, y],
                bc: FaceBC::Velocity {
                    u: [0.0, 0.0, -regime.u_jet_lattice.min(0.17)],
                },
            });
            x += 1;
        }
    }
    patches
}

fn inside_any_nozzle(
    input: &ProtocolInput,
    regime: &protocol::Regime,
    points: &[[f64; 2]],
    radius_m: f64,
    x: usize,
    y: usize,
) -> bool {
    let xp = x as f64 * regime.dx;
    let yp = y as f64 * regime.dx;
    points.iter().any(|pt| {
        let cx = pt[0] * input.target.width_m;
        let cy = pt[1] * input.target.depth_m;
        let dx = xp - cx;
        let dy = yp - cy;
        dx * dx + dy * dy <= radius_m * radius_m
    })
}

fn run_until_quasi_steady(
    sim: &mut Sim3,
    speed_scale: f64,
    max_steps: usize,
) -> anyhow::Result<()> {
    if max_steps == 0 || speed_scale == 0.0 {
        return Ok(());
    }
    // Numerical spin-up criterion: once the max node-wise velocity change per
    // step is below 0.1% of the imposed jet speed, the resolved jet field is
    // quasi-steady for particle insertion. This is a convergence residual, not
    // a physical closure term.
    let tolerance = 1.0e-3 * speed_scale;
    let mut prev = gather_velocity(sim);
    for step in 1..=max_steps {
        sim.step();
        let current = gather_velocity(sim);
        let max_delta = current
            .iter()
            .zip(prev.iter())
            .map(|(a, b)| {
                let dx = a[0] - b[0];
                let dy = a[1] - b[1];
                let dz = a[2] - b[2];
                (dx * dx + dy * dy + dz * dz).sqrt()
            })
            .fold(0.0, f64::max);
        if max_delta <= tolerance {
            return Ok(());
        }
        prev = current;
        if step == max_steps {
            anyhow::bail!(
                "tray flow did not reach the quasi-steady residual criterion within {max_steps} steps; max_delta={max_delta:.3e}, tolerance={tolerance:.3e}"
            );
        }
    }
    Ok(())
}

fn gather_velocity(sim: &Sim3) -> Vec<[f64; 3]> {
    let ux = sim.gather_ux();
    let uy = sim.gather_uy();
    let uz = sim.gather_uz();
    (0..ux.len()).map(|i| [ux[i], uy[i], uz[i]]).collect()
}

fn gather_velocity_si(sim: &Sim3, regime: &protocol::Regime) -> Vec<[f64; 3]> {
    gather_velocity(sim)
        .into_iter()
        .map(|u| {
            [
                u[0] * regime.dx / regime.dt,
                u[1] * regime.dx / regime.dt,
                u[2] * regime.dx / regime.dt,
            ]
        })
        .collect()
}
