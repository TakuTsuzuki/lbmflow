use lbm_core::prelude::*;

fn closed_box(dims: [usize; 3]) -> (Vec<bool>, Vec<[f64; 3]>) {
    let mut walls = WallSpec::<f64>::default();
    walls.is_wall = [true; 6];
    build_wall_rims(3, dims, &walls)
}

fn base_spec(dims: [usize; 3]) -> GlobalSpec<f64> {
    GlobalSpec {
        dims,
        nu: 0.08,
        collision: CollisionKind::Bgk,
        periodic: [false; 3],
        faces: [FaceBC::Closed; 6],
        force: [0.0; 3],
        sources: Vec::new(),
        face_patches: Vec::new(),
    }
}

#[test]
fn volume_source_validation_errors() {
    let dims = [8, 8, 8];
    let (solid, _) = closed_box(dims);

    let mut touching = base_spec(dims);
    touching.sources.push(VolumeSource {
        region: SourceRegion {
            lo: [0, 2, 2],
            hi: [1, 2, 2],
        },
        kind: SourceKind::MassFlow { q_lu: 1e-4 },
    });
    assert!(matches!(
        touching.validate(3, &solid),
        Err(SpecError::SourceRegionNotInterior { .. })
    ));

    let mut overlap = base_spec(dims);
    overlap.sources.push(VolumeSource {
        region: SourceRegion {
            lo: [2, 2, 2],
            hi: [3, 3, 3],
        },
        kind: SourceKind::MassFlow { q_lu: 1e-4 },
    });
    overlap.sources.push(VolumeSource {
        region: SourceRegion {
            lo: [3, 3, 3],
            hi: [4, 4, 4],
        },
        kind: SourceKind::MassFlow { q_lu: 1e-4 },
    });
    assert!(matches!(
        overlap.validate(3, &solid),
        Err(SpecError::SourceOverlap { .. })
    ));

    let mut fast_jet = base_spec(dims);
    fast_jet.sources.push(VolumeSource {
        region: SourceRegion {
            lo: [2, 2, 2],
            hi: [2, 2, 2],
        },
        kind: SourceKind::Jet {
            q_lu: 1e-4,
            u: [0.31, 0.0, 0.0],
        },
    });
    assert!(matches!(
        fast_jet.validate(3, &solid),
        Err(SpecError::VelocityTooHigh { .. })
    ));

    let mut strong_sink = base_spec(dims);
    strong_sink.sources.push(VolumeSource {
        region: SourceRegion {
            lo: [2, 2, 2],
            hi: [2, 2, 2],
        },
        kind: SourceKind::MassFlow { q_lu: -1.0 },
    });
    assert!(matches!(
        strong_sink.validate(3, &solid),
        Err(SpecError::SourceSinkTooStrong { .. })
    ));
}

#[test]
fn mass_ledger_matches_sum_q_lu() {
    let dims = [8, 8, 8];
    let (solid, wall_u) = closed_box(dims);
    let mut spec = base_spec(dims);
    spec.sources.push(VolumeSource {
        region: SourceRegion {
            lo: [3, 3, 3],
            hi: [3, 3, 3],
        },
        kind: SourceKind::MassFlow { q_lu: 2.5e-5 },
    });

    let mut sim: Solver<D3Q19, f64, CpuScalar, LocalPeriodic> = Solver::new(
        &spec,
        &solid,
        &wall_u,
        [1, 1, 1],
        CpuScalar::default(),
        LocalPeriodic,
    );
    let m0 = sim.total_mass() as f64;
    let steps = 8;
    sim.run(steps);
    let dm = sim.total_mass() as f64 - m0;
    let expected = steps as f64 * 2.5e-5;
    assert!(
        (dm - expected).abs() <= 1e-12,
        "mass ledger drift: got {dm:e}, expected {expected:e}"
    );
}
