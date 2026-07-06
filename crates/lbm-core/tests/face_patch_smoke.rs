use lbm_core::prelude::*;

fn top_patched_box(dims: [usize; 3], patches: &[FacePatch<f64>]) -> (Vec<bool>, Vec<[f64; 3]>) {
    let n = dims[0] * dims[1] * dims[2];
    let mut solid = vec![false; n];
    let wall_u = vec![[0.0; 3]; n];
    for y in 0..dims[1] {
        for x in 0..dims[0] {
            let on_wall = x == 0 || x + 1 == dims[0] || y == 0 || y + 1 == dims[1];
            if !on_wall {
                continue;
            }
            let in_top_patch = y + 1 == dims[1]
                && patches.iter().any(|p| {
                    p.face == Face::YPos.index()
                        && p.lo[0] <= x
                        && x <= p.hi[0]
                        && p.lo[1] == 0
                        && p.hi[1] == 0
                });
            if !in_top_patch {
                solid[y * dims[0] + x] = true;
            }
        }
    }
    (solid, wall_u)
}

fn base_spec(dims: [usize; 3]) -> GlobalSpec<f64> {
    GlobalSpec {
        dims,
        nu: 0.08,
        collision: CollisionKind::Bgk,
        periodic: [false, false, false],
        faces: [FaceBC::Closed; 6],
        force: [0.0; 3],
        sources: Vec::new(),
        face_patches: Vec::new(),
    }
}

#[test]
fn face_patch_validation_errors() {
    let dims = [12, 10, 1];

    let mut out_of_bounds = base_spec(dims);
    out_of_bounds.face_patches.push(FacePatch {
        face: Face::YPos.index(),
        lo: [10, 0],
        hi: [12, 0],
        bc: FaceBC::Outflow,
    });
    let mut walls = WallSpec::<f64>::default();
    walls.is_wall = [true; 6];
    let (solid, _) = build_wall_rims(2, dims, &walls);
    let err = out_of_bounds.validate(2, &solid);
    assert!(
        matches!(err, Err(SpecError::FacePatchOutOfBounds { .. })),
        "{err:?}"
    );

    let mut overlap = base_spec(dims);
    overlap.face_patches.push(FacePatch {
        face: Face::YPos.index(),
        lo: [4, 0],
        hi: [6, 0],
        bc: FaceBC::Velocity {
            u: [0.0, -0.04, 0.0],
        },
    });
    overlap.face_patches.push(FacePatch {
        face: Face::YPos.index(),
        lo: [6, 0],
        hi: [8, 0],
        bc: FaceBC::Outflow,
    });
    let (solid, _) = top_patched_box(dims, &overlap.face_patches);
    assert!(matches!(
        overlap.validate(2, &solid),
        Err(SpecError::FacePatchOverlap { .. })
    ));

    let mut two_axes = base_spec(dims);
    two_axes.face_patches.push(FacePatch {
        face: Face::YPos.index(),
        lo: [4, 0],
        hi: [6, 0],
        bc: FaceBC::Outflow,
    });
    two_axes.face_patches.push(FacePatch {
        face: Face::XPos.index(),
        lo: [3, 0],
        hi: [5, 0],
        bc: FaceBC::Outflow,
    });
    let (solid, _) = top_patched_box(dims, &two_axes.face_patches);
    assert!(matches!(
        two_axes.validate(2, &solid),
        Err(SpecError::OpenFacesOnMultipleAxes)
    ));
}

#[test]
fn top_velocity_patch_drives_local_flow() {
    let dims = [16, 12, 1];
    let patches = vec![
        FacePatch {
            face: Face::YPos.index(),
            lo: [6, 0],
            hi: [8, 0],
            bc: FaceBC::Velocity {
                u: [0.0, -0.04, 0.0],
            },
        },
        FacePatch {
            face: Face::YPos.index(),
            lo: [10, 0],
            hi: [11, 0],
            bc: FaceBC::Outflow,
        },
    ];
    let (solid, wall_u) = top_patched_box(dims, &patches);
    let mut spec = base_spec(dims);
    spec.face_patches = patches;
    spec.validate(2, &solid).unwrap();

    let mut sim: Solver<D2Q9, f64, CpuScalar, LocalPeriodic> = Solver::new(
        &spec,
        &solid,
        &wall_u,
        [1, 1, 1],
        CpuScalar::default(),
        LocalPeriodic,
    );
    sim.run(20);
    let u = sim.u(7, dims[1] - 1, 0);
    assert!(
        u[1] < -1e-3,
        "patched inlet did not drive downward flow: {u:?}"
    );
}
