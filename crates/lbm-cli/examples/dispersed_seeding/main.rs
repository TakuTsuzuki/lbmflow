mod particles;
mod protocol;
mod readout;
mod reservoir;

use lbm_core::prelude::*;
use particles::{deposit_batch, make_reservoir_particles};
use protocol::ProtocolInput;
use reservoir::extract_by_depth;
use std::path::PathBuf;

type Sim3 = Solver<D3Q19, f64, CpuScalar, LocalPeriodic>;

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
    run_lbm_steps(&mut reservoir_sim, 40);
    let reservoir_velocity = gather_velocity(&reservoir_sim);

    let reservoir_particles = make_reservoir_particles(&input);
    let extraction = extract_by_depth(&input, &reservoir_particles);
    let mut batch = extraction.batch;

    let mut tray_sim = build_tray_sim(&input, &regime);
    run_lbm_steps(&mut tray_sim, 90);
    let tray_velocity = gather_velocity(&tray_sim)
        .into_iter()
        .map(|u| {
            [
                u[0] * regime.dx / regime.dt,
                u[1] * regime.dx / regime.dt,
                u[2] * regime.dx / regime.dt,
            ]
        })
        .collect::<Vec<_>>();

    deposit_batch(&input, &regime, &mut batch, &tray_velocity)?;
    let metrics = readout::write_outputs(
        &input,
        &regime,
        &batch,
        batch.len(),
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
    for op in &input.protocol {
        if op.op == "agitate" && op.pattern.as_deref() != Some("translational") {
            anyhow::bail!(
                "unsupported agitate pattern {:?}; only translational is implemented",
                op.pattern
            );
        }
        if op.op == "eject" && op.nozzle_diameter_m.is_none() {
            anyhow::bail!("eject.nozzle_diameter_m is required");
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
        force: [0.0, 0.0, -1.0e-7],
        ..Default::default()
    };
    let (solid, wall_u) = build_wall_rims(3, dims, &walls);
    Sim3::new(
        &spec,
        &solid,
        &wall_u,
        [1, 1, 1],
        CpuScalar::default(),
        LocalPeriodic,
    )
}

fn build_tray_sim(input: &ProtocolInput, regime: &protocol::Regime) -> Sim3 {
    let dims = [input.grid.tray_nx, input.grid.tray_ny, input.grid.tray_nz];
    let mut walls = WallSpec::<f64>::default();
    for face in [Face::XNeg, Face::XPos, Face::YNeg, Face::YPos, Face::ZNeg] {
        walls.is_wall[face.index()] = true;
    }
    let mut faces = [FaceBC::Closed; 6];
    faces[Face::ZPos.index()] = FaceBC::Velocity {
        u: [0.0, 0.0, -regime.u_jet_lattice.min(0.08)],
    };
    let spec = GlobalSpec::<f64> {
        dims,
        nu: regime.nu_lattice,
        periodic: [false, false, false],
        faces,
        collision: CollisionKind::Trt {
            magic: CollisionKind::MAGIC_STD,
        },
        ..Default::default()
    };
    let (solid, wall_u) = build_wall_rims(3, dims, &walls);
    let mut sim = Sim3::new(
        &spec,
        &solid,
        &wall_u,
        [1, 1, 1],
        CpuScalar::default(),
        LocalPeriodic,
    );
    let points = input
        .op("eject")
        .and_then(|op| op.points_xy_frac.clone())
        .unwrap_or_else(|| vec![[0.5, 0.5]]);
    let nozzle_radius_m = 0.5 * regime.nozzle_d_m;
    sim.set_inlet_profile_with(Face::ZPos, |x, y| {
        let xp = x as f64 * regime.dx;
        let yp = y as f64 * regime.dx;
        let mut inside_patch = false;
        for pt in &points {
            let cx = pt[0] * input.target.width_m;
            let cy = pt[1] * input.target.depth_m;
            let dx = xp - cx;
            let dy = yp - cy;
            inside_patch |= dx * dx + dy * dy <= nozzle_radius_m * nozzle_radius_m;
        }
        if inside_patch {
            [0.0, 0.0, -regime.u_jet_lattice.min(0.17)]
        } else {
            [0.0, 0.0, 0.0]
        }
    });
    sim
}

fn run_lbm_steps(sim: &mut Sim3, steps: usize) {
    for _ in 0..steps {
        sim.step();
    }
}

fn gather_velocity(sim: &Sim3) -> Vec<[f64; 3]> {
    let ux = sim.gather_ux();
    let uy = sim.gather_uy();
    let uz = sim.gather_uz();
    (0..ux.len()).map(|i| [ux[i], uy[i], uz[i]]).collect()
}
