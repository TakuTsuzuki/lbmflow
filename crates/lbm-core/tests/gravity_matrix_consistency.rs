//! Radar-remainder follow-up to V&V lane 5.1: gravity/force composition matrix.
//!
//! This pins the documented Guo composition point:
//! `F_total(x,t) = F_uniform + F_cell(x,t) + rho(x,t) * g`.

use lbm_core::compat::prelude::{Collision, EdgeBC, Edges, SimConfig, Simulation};

const NX: usize = 32;
const NY: usize = 32;
const STEPS: usize = 500;
const MAGNITUDES: [f64; 4] = [0.0, 1.0e-6, 1.0e-5, 1.0e-4];
const MOMENTUM_REL_BAND: f64 = 1.0e-11;
const MOMENTUM_ABS_FLOOR: f64 = 5.0e-14;
const MASS_REL_BAND: f64 = 1.0e-12;
const OMEGA: f64 = 1.0e-2;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum ForcePath {
    UniformConfigForce,
    SetGravity,
    PerCellForceField,
    TimeVaryingForceField,
}

impl ForcePath {
    const ALL: [Self; 4] = [
        Self::UniformConfigForce,
        Self::SetGravity,
        Self::PerCellForceField,
        Self::TimeVaryingForceField,
    ];

    fn name(self) -> &'static str {
        match self {
            Self::UniformConfigForce => "uniform SimConfig::force",
            Self::SetGravity => "uniform via sim.set_gravity",
            Self::PerCellForceField => "per-cell force_field",
            Self::TimeVaryingForceField => "time-varying force_field",
        }
    }
}

#[derive(Debug)]
enum Cell {
    Pass(Snapshot),
    Fail(String),
    Skip(String),
}

impl Cell {
    fn label(&self) -> &'static str {
        match self {
            Self::Pass(_) => "PASS",
            Self::Fail(_) => "FAIL",
            Self::Skip(_) => "SKIP",
        }
    }

    fn reason(&self) -> &str {
        match self {
            Self::Pass(_) => "",
            Self::Fail(reason) | Self::Skip(reason) => reason,
        }
    }

    fn snapshot(&self) -> Option<&Snapshot> {
        match self {
            Self::Pass(snapshot) => Some(snapshot),
            Self::Fail(_) | Self::Skip(_) => None,
        }
    }

    fn is_fail(&self) -> bool {
        matches!(self, Self::Fail(_))
    }
}

#[derive(Clone, Debug)]
struct Snapshot {
    rho: Vec<u64>,
    ux: Vec<u64>,
    uy: Vec<u64>,
    mass: u64,
    momentum: [u64; 2],
}

#[test]
fn gravity_force_matrix_conserves_mass_and_momentum() {
    let mut failures = Vec::new();

    for magnitude in MAGNITUDES {
        let mut row = Vec::new();

        for path in ForcePath::ALL {
            let cell = run_cell(magnitude, path)
                .map(Cell::Pass)
                .unwrap_or_else(Cell::Fail);
            println!(
                "GRAVITY_MATRIX | magnitude={magnitude:.1e} | path={} | {} | {}",
                path.name(),
                cell.label(),
                cell.reason()
            );
            if cell.is_fail() {
                failures.push(format!(
                    "magnitude={magnitude:.1e}, path={}: {}",
                    path.name(),
                    cell.reason()
                ));
            }
            row.push((path, cell));
        }

        if magnitude != 0.0 {
            let uniform = row
                .iter()
                .find(|(path, _)| *path == ForcePath::UniformConfigForce)
                .and_then(|(_, cell)| cell.snapshot());
            let gravity = row
                .iter()
                .find(|(path, _)| *path == ForcePath::SetGravity)
                .and_then(|(_, cell)| cell.snapshot());
            let bit_cell = match (uniform, gravity) {
                (Some(a), Some(b)) => compare_snapshots(a, b),
                _ => Cell::Skip(
                    "uniform-force or set_gravity cell failed; bit-identical comparison skipped"
                        .to_string(),
                ),
            };
            println!(
                "GRAVITY_MATRIX_BITWISE | magnitude={magnitude:.1e} | uniform SimConfig::force vs sim.set_gravity | {} | {}",
                bit_cell.label(),
                bit_cell.reason()
            );
            if bit_cell.is_fail() {
                failures.push(format!(
                    "magnitude={magnitude:.1e}, uniform force vs set_gravity bit identity: {}",
                    bit_cell.reason()
                ));
            }
        }
    }

    assert!(
        failures.is_empty(),
        "gravity/force consistency matrix failures:\n{}",
        failures.join("\n")
    );
}

fn run_cell(magnitude: f64, path: ForcePath) -> Result<Snapshot, String> {
    let mut sim = build_sim(magnitude, path)?;
    let n_fluid = sim.fluid_cell_count() as f64;
    let mass0 = sim.total_mass_f64();
    let mut measured_impulse = [0.0, 0.0];
    let mut expected_impulse = [0.0, 0.0];

    for step in 0..STEPS {
        let force = force_for_step(magnitude, path, step + 1);
        if path == ForcePath::TimeVaryingForceField {
            set_uniform_force_field(&mut sim, force);
        }

        assert_finite(&sim, magnitude, path, step, "before step")?;
        let mass_before = sim.total_mass_f64();
        let p0 = sim.total_momentum();
        sim.step();
        assert_finite(&sim, magnitude, path, step, "after step")?;
        let p1 = sim.total_momentum();

        let expected = match path {
            ForcePath::SetGravity => [0.0, -mass_before * magnitude],
            _ => [n_fluid * force[0], n_fluid * force[1]],
        };
        measured_impulse[0] += p1[0] - p0[0];
        measured_impulse[1] += p1[1] - p0[1];
        expected_impulse[0] += expected[0];
        expected_impulse[1] += expected[1];
        assert_momentum_cumulative(magnitude, path, step, measured_impulse, expected_impulse)?;
    }

    let mass1 = sim.total_mass_f64();
    let mass_rel = ((mass1 - mass0) / mass0).abs();
    if mass_rel > MASS_REL_BAND {
        return Err(format!(
            "mass drift rel={mass_rel:.3e} > {MASS_REL_BAND:.1e}, m0={mass0:.12e}, m1={mass1:.12e}"
        ));
    }

    Ok(snapshot(&sim))
}

fn build_sim(magnitude: f64, path: ForcePath) -> Result<Simulation<f64>, String> {
    let mut sim = SimConfig {
        nx: NX,
        ny: NY,
        nu: 1.0 / 6.0,
        collision: Collision::Trt {
            magic: Collision::MAGIC_STD,
        },
        edges: Edges {
            left: EdgeBC::Periodic,
            right: EdgeBC::Periodic,
            bottom: EdgeBC::Periodic,
            top: EdgeBC::Periodic,
        },
        force: match path {
            ForcePath::UniformConfigForce => [0.0, -magnitude],
            ForcePath::SetGravity
            | ForcePath::PerCellForceField
            | ForcePath::TimeVaryingForceField => [0.0, 0.0],
        },
    }
    .build()
    .map_err(|err| format!("build failed: {err}"))?;

    match path {
        ForcePath::UniformConfigForce => {}
        ForcePath::SetGravity => sim.set_gravity([0.0, -magnitude]),
        ForcePath::PerCellForceField => set_uniform_force_field(&mut sim, [0.0, -magnitude]),
        ForcePath::TimeVaryingForceField => set_uniform_force_field(&mut sim, [0.0, 0.0]),
    }

    Ok(sim)
}

fn force_for_step(magnitude: f64, path: ForcePath, step: usize) -> [f64; 2] {
    match path {
        ForcePath::UniformConfigForce | ForcePath::SetGravity | ForcePath::PerCellForceField => {
            [0.0, -magnitude]
        }
        ForcePath::TimeVaryingForceField => {
            let amp = -magnitude * (OMEGA * step as f64).sin();
            [0.0, amp]
        }
    }
}

fn set_uniform_force_field(sim: &mut Simulation<f64>, force: [f64; 2]) {
    sim.force_field_mut().fill(force);
}

fn assert_finite(
    sim: &Simulation<f64>,
    magnitude: f64,
    path: ForcePath,
    step: usize,
    when: &str,
) -> Result<(), String> {
    if !sim.total_mass_f64().is_finite() {
        return Err(format!(
            "non-finite total mass {when}, step={step}, magnitude={magnitude:.1e}, path={}",
            path.name()
        ));
    }
    let p = sim.total_momentum();
    if !p[0].is_finite() || !p[1].is_finite() {
        return Err(format!(
            "non-finite total momentum {when}, step={step}, magnitude={magnitude:.1e}, path={}, p=({:.12e},{:.12e})",
            path.name(),
            p[0],
            p[1]
        ));
    }
    for (i, &rho) in sim.rho_field().iter().enumerate() {
        if !rho.is_finite() {
            return Err(format!(
                "non-finite rho[{i}] {when}, step={step}, magnitude={magnitude:.1e}, path={}",
                path.name()
            ));
        }
    }
    for (i, &ux) in sim.ux_field().iter().enumerate() {
        if !ux.is_finite() {
            return Err(format!(
                "non-finite ux[{i}] {when}, step={step}, magnitude={magnitude:.1e}, path={}",
                path.name()
            ));
        }
    }
    for (i, &uy) in sim.uy_field().iter().enumerate() {
        if !uy.is_finite() {
            return Err(format!(
                "non-finite uy[{i}] {when}, step={step}, magnitude={magnitude:.1e}, path={}",
                path.name()
            ));
        }
    }
    Ok(())
}

fn assert_momentum_cumulative(
    magnitude: f64,
    path: ForcePath,
    step: usize,
    measured: [f64; 2],
    expected: [f64; 2],
) -> Result<(), String> {
    for axis in 0..2 {
        let err = (measured[axis] - expected[axis]).abs();
        let band = (MOMENTUM_REL_BAND * expected[axis].abs()).max(MOMENTUM_ABS_FLOOR);
        if err > band {
            let rel = err / expected[axis].abs().max(MOMENTUM_ABS_FLOOR);
            return Err(format!(
                "cumulative momentum ledger failed through step={step}, axis={axis}, rel={rel:.3e}, abs={err:.12e}, band={band:.12e}, measured={:.12e}, expected={:.12e}, magnitude={magnitude:.1e}, path={}",
                measured[axis],
                expected[axis],
                path.name()
            ));
        }
    }
    Ok(())
}

fn snapshot(sim: &Simulation<f64>) -> Snapshot {
    Snapshot {
        rho: sim.rho_field().iter().map(|v| v.to_bits()).collect(),
        ux: sim.ux_field().iter().map(|v| v.to_bits()).collect(),
        uy: sim.uy_field().iter().map(|v| v.to_bits()).collect(),
        mass: sim.total_mass_f64().to_bits(),
        momentum: sim.total_momentum().map(f64::to_bits),
    }
}

fn compare_snapshots(a: &Snapshot, b: &Snapshot) -> Cell {
    if a.mass != b.mass {
        return Cell::Fail(format!(
            "mass bits differ: uniform=0x{:016x}, gravity=0x{:016x}",
            a.mass, b.mass
        ));
    }
    if a.momentum != b.momentum {
        return Cell::Fail(format!(
            "momentum bits differ: uniform=[0x{:016x},0x{:016x}], gravity=[0x{:016x},0x{:016x}]",
            a.momentum[0], a.momentum[1], b.momentum[0], b.momentum[1]
        ));
    }
    compare_bits("rho", &a.rho, &b.rho)
        .or_else(|| compare_bits("ux", &a.ux, &b.ux))
        .or_else(|| compare_bits("uy", &a.uy, &b.uy))
        .map(Cell::Fail)
        .unwrap_or_else(|| Cell::Pass(a.clone()))
}

fn compare_bits(name: &str, a: &[u64], b: &[u64]) -> Option<String> {
    a.iter()
        .zip(b.iter())
        .enumerate()
        .find(|(_, (av, bv))| av != bv)
        .map(|(i, (av, bv))| {
            format!("{name}[{i}] bits differ: uniform=0x{av:016x}, gravity=0x{bv:016x}")
        })
}
