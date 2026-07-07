use anyhow::Result;
use lbm_core::prelude::{Backend, CpuScalar, CpuSimd, D2Q9};
use lbm_scenario::CollisionSpec;
use serde::Serialize;

// Static product facts: keep in sync with docs/LIMITATIONS.md.
pub const STATIC_FACTS: StaticFacts = StaticFacts {
    d3q27_open_face_restriction:
        "D3Q27 supports CPU full boundary coverage including periodic, closed-wall, velocity-inlet, pressure-outlet, outflow, and convective open faces; GPU open faces are rejected explicitly; scenario JSON exposure landed 2026-07-07",
    checkpoint_scope: "multi-part and per-rank",
    particle_coupling: "one-way",
};

pub struct StaticFacts {
    pub d3q27_open_face_restriction: &'static str,
    pub checkpoint_scope: &'static str,
    pub particle_coupling: &'static str,
}

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
pub struct CapabilityMatrix {
    pub lattices: Vec<LatticeCapability>,
    pub collisions: CollisionCapabilities,
    pub precisions: PrecisionCapabilities,
    pub backends: Vec<BackendCapability>,
    pub backend_gravity_fallback: &'static str,
    pub checkpoint: CheckpointCapability,
    pub particle_coupling: &'static str,
}

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
pub struct LatticeCapability {
    pub name: &'static str,
    pub dimension: u8,
    pub status: &'static str,
    pub restrictions: Vec<&'static str>,
}

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
pub struct CollisionCapabilities {
    pub core: Vec<&'static str>,
    pub scenario_path: Vec<&'static str>,
}

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
pub struct PrecisionCapabilities {
    pub compute: Vec<&'static str>,
    pub storage: Vec<&'static str>,
    pub notes: Vec<&'static str>,
}

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
pub struct BackendCapability {
    pub name: &'static str,
    pub compiled: bool,
    pub gravity_body_force: bool,
    pub notes: Vec<&'static str>,
}

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
pub struct CheckpointCapability {
    pub scope: &'static str,
}

pub fn run(json: bool) -> Result<()> {
    let matrix = matrix();
    if json {
        println!("{}", serde_json::to_string_pretty(&matrix)?);
    } else {
        print_human(&matrix);
    }
    Ok(())
}

pub fn matrix() -> CapabilityMatrix {
    CapabilityMatrix {
        lattices: vec![
            LatticeCapability {
                name: "d2q9",
                dimension: 2,
                status: "supported",
                restrictions: Vec::new(),
            },
            LatticeCapability {
                name: "d3q19",
                dimension: 3,
                status: "supported",
                restrictions: Vec::new(),
            },
            LatticeCapability {
                name: "d3q27",
                dimension: 3,
                status: "restricted",
                restrictions: vec![STATIC_FACTS.d3q27_open_face_restriction],
            },
        ],
        collisions: CollisionCapabilities {
            core: vec!["bgk", "trt", "central_moment"],
            scenario_path: scenario_collision_names(),
        },
        precisions: PrecisionCapabilities {
            compute: vec!["f32", "f64"],
            storage: vec!["f32", "f64", "f16"],
            notes: vec!["f16 is GPU storage mode; arithmetic remains f32"],
        },
        backends: vec![
            BackendCapability {
                name: "cpu-scalar",
                compiled: true,
                gravity_body_force: cpu_scalar_gravity_body_force(),
                notes: vec!["CpuScalar is always available in the core"],
            },
            BackendCapability {
                name: "cpu-simd",
                compiled: true,
                gravity_body_force: cpu_simd_gravity_body_force(),
                notes: vec!["CpuSimd is always available in the core"],
            },
            BackendCapability {
                name: "gpu",
                compiled: cfg!(feature = "gpu"),
                gravity_body_force: true,
                notes: vec![
                    "WgpuBackend is compiled only with --features gpu",
                    "WgpuBackend supports backend-side rho*g composition when compiled",
                ],
            },
            BackendCapability {
                name: "mpi",
                compiled: cfg!(feature = "mpi"),
                gravity_body_force: false,
                notes: vec![
                    "MPI support requires a native MPI toolchain and --features mpi",
                    "MPI is a distribution layer; gravity composition follows the selected compute backend",
                ],
            },
        ],
        backend_gravity_fallback:
            "Backends that do not support backend-side gravity use the host-staged force-field fallback",
        checkpoint: CheckpointCapability {
            scope: STATIC_FACTS.checkpoint_scope,
        },
        particle_coupling: STATIC_FACTS.particle_coupling,
    }
}

fn scenario_collision_names() -> Vec<&'static str> {
    [
        CollisionSpec::Bgk,
        CollisionSpec::Trt,
        CollisionSpec::CentralMoment,
        CollisionSpec::DeprecatedCumulantAlias,
    ]
    .into_iter()
    .map(|collision| match collision {
        CollisionSpec::Bgk => "bgk",
        CollisionSpec::Trt => "trt",
        // Honored on the 3D CPU scenario path only; other paths reject it
        // explicitly. `cumulant` is still accepted as a deprecated schema
        // alias, so the capability output lists both accepted spellings.
        CollisionSpec::CentralMoment => "central_moment",
        CollisionSpec::DeprecatedCumulantAlias => "cumulant",
    })
    .collect()
}

fn cpu_scalar_gravity_body_force() -> bool {
    <CpuScalar as Backend<D2Q9, f64>>::supports_gravity_body_force(&CpuScalar::default())
}

fn cpu_simd_gravity_body_force() -> bool {
    <CpuSimd as Backend<D2Q9, f64>>::supports_gravity_body_force(&CpuSimd::default())
}

fn print_human(matrix: &CapabilityMatrix) {
    println!("LBMFlow capability matrix");
    println!();
    println!("Lattices");
    println!(
        "{:<10} {:<10} {:<12} Restrictions",
        "name", "dimension", "status"
    );
    for l in &matrix.lattices {
        let restrictions = if l.restrictions.is_empty() {
            "-".to_string()
        } else {
            l.restrictions.join("; ")
        };
        println!(
            "{:<10} {:<10} {:<12} {}",
            l.name, l.dimension, l.status, restrictions
        );
    }
    println!();
    println!("Collisions");
    println!("  core:          {}", matrix.collisions.core.join(", "));
    println!(
        "  scenario path: {}",
        matrix.collisions.scenario_path.join(", ")
    );
    println!();
    println!("Precision/storage");
    println!("  compute: {}", matrix.precisions.compute.join(", "));
    println!("  storage: {}", matrix.precisions.storage.join(", "));
    for note in &matrix.precisions.notes {
        println!("  note: {note}");
    }
    println!();
    println!("Backends");
    println!(
        "{:<12} {:<8} {:<15} Notes",
        "name", "compiled", "gravity force"
    );
    for b in &matrix.backends {
        println!(
            "{:<12} {:<8} {:<15} {}",
            b.name,
            b.compiled,
            b.gravity_body_force,
            b.notes.join("; ")
        );
    }
    println!("  fallback: {}", matrix.backend_gravity_fallback);
    println!();
    println!("Checkpoint scope: {}", matrix.checkpoint.scope);
    println!("Particle coupling: {}", matrix.particle_coupling);
}
