use anyhow::Result;
use lbm_scenario::CollisionSpec;
use serde::Serialize;

// Static product facts: keep in sync with docs/LIMITATIONS.md.
const STATIC_FACTS: StaticFacts = StaticFacts {
    d3q27_open_face_restriction:
        "D3Q27 supports periodic and closed-wall cases only; open faces are unsupported",
    checkpoint_scope: "single-rank",
    particle_coupling: "one-way",
};

struct StaticFacts {
    d3q27_open_face_restriction: &'static str,
    checkpoint_scope: &'static str,
    particle_coupling: &'static str,
}

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
pub struct CapabilityMatrix {
    lattices: Vec<LatticeCapability>,
    collisions: CollisionCapabilities,
    precisions: PrecisionCapabilities,
    backends: Vec<BackendCapability>,
    checkpoint: CheckpointCapability,
    particle_coupling: &'static str,
}

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
struct LatticeCapability {
    name: &'static str,
    dimension: u8,
    status: &'static str,
    restrictions: Vec<&'static str>,
}

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
struct CollisionCapabilities {
    core: Vec<&'static str>,
    scenario_path: Vec<&'static str>,
}

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
struct PrecisionCapabilities {
    compute: Vec<&'static str>,
    storage: Vec<&'static str>,
    notes: Vec<&'static str>,
}

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
struct BackendCapability {
    name: &'static str,
    compiled: bool,
    notes: Vec<&'static str>,
}

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
struct CheckpointCapability {
    scope: &'static str,
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

fn matrix() -> CapabilityMatrix {
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
                name: "cpu",
                compiled: true,
                notes: vec!["CpuScalar and CpuSimd are available in the core"],
            },
            BackendCapability {
                name: "gpu",
                compiled: cfg!(feature = "gpu"),
                notes: vec!["wgpu backend is compiled only with --features gpu"],
            },
            BackendCapability {
                name: "mpi",
                compiled: cfg!(feature = "mpi"),
                notes: vec!["MPI support requires a native MPI toolchain and --features mpi"],
            },
        ],
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
    ]
    .into_iter()
    .map(|collision| match collision {
        CollisionSpec::Bgk => "bgk",
        CollisionSpec::Trt => "trt",
        // Honored on the 3D D3Q19 CPU scenario path only; other paths
        // reject it explicitly.
        CollisionSpec::CentralMoment | CollisionSpec::DeprecatedCumulantAlias => "central_moment",
    })
    .collect()
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
    println!("{:<10} {:<8} Notes", "name", "compiled");
    for b in &matrix.backends {
        println!("{:<10} {:<8} {}", b.name, b.compiled, b.notes.join("; "));
    }
    println!();
    println!("Checkpoint scope: {}", matrix.checkpoint.scope);
    println!("Particle coupling: {}", matrix.particle_coupling);
}
