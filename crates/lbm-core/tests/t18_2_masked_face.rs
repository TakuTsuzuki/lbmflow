//! T18.2 adversarial acceptance tests for CR-2 per-cell masked face BCs.
//!
//! These tests intentionally target the frozen public API from
//! docs/DISPERSED_DEPOSITION.md: `FacePatch { face, lo, hi, bc }` with
//! in-face coordinates ordered by the remaining axes ascending, plus
//! `GlobalSpec::face_patches`. The implementation is not expected to exist in
//! this branch.

use lbm_core::lattice::D3Q19;
use lbm_core::params::{FacePatch, MAX_SPEED};
use lbm_core::prelude::*;

type Cpu3<H> = Solver<D3Q19, f64, CpuScalar, H>;

const IMPINGING_DIMS: [usize; 3] = [32, 32, 24];
const IMPINGING_STEPS: usize = 1_200;
const IMPINGING_MASS_DRIFT_REL_BAND: f64 = 2.0e-10;
const EQUIV_FIELD_ABS_BAND: f64 = 1.0e-14;

fn closed_faces() -> [FaceBC<f64>; 6] {
    [FaceBC::Closed; 6]
}

fn walls_for_faces(faces: &[Face]) -> WallSpec<f64> {
    let mut walls = WallSpec::default();
    for &face in faces {
        walls.is_wall[face.index()] = true;
    }
    walls
}

fn build<H: HaloExchange<f64>>(
    spec: &GlobalSpec<f64>,
    walls: &WallSpec<f64>,
    decomp: [usize; 3],
    exchange: H,
) -> Cpu3<H> {
    let (solid, wall_u) = build_wall_rims(3, spec.dims, walls);
    Solver::new(
        spec,
        &solid,
        &wall_u,
        decomp,
        CpuScalar::default(),
        exchange,
    )
}

fn z_top_patch(lo: [usize; 2], hi: [usize; 2], bc: FaceBC<f64>) -> FacePatch<f64> {
    // Z face in-face coordinates are [x, y].
    FacePatch {
        face: Face::ZPos.index(),
        lo,
        hi,
        bc,
    }
}

fn x_neg_patch(lo: [usize; 2], hi: [usize; 2], bc: FaceBC<f64>) -> FacePatch<f64> {
    // X face in-face coordinates are [y, z].
    FacePatch {
        face: Face::XNeg.index(),
        lo,
        hi,
        bc,
    }
}

fn impinging_jet_spec() -> (GlobalSpec<f64>, WallSpec<f64>) {
    let mut faces = closed_faces();
    faces[Face::ZPos.index()] = FaceBC::Closed;

    let jet = FaceBC::Velocity {
        u: [0.0, 0.0, -0.05],
    };
    let outlet = FaceBC::Pressure { rho: 1.0 };

    // A rectangular four-piece annulus around the central inlet. The CR-2 API
    // only exposes rectangular patches, so the annulus is represented as
    // non-overlapping rectangles on the same top face.
    let face_patches = vec![
        z_top_patch([14, 14], [17, 17], jet),
        z_top_patch([11, 11], [20, 13], outlet),
        z_top_patch([11, 18], [20, 20], outlet),
        z_top_patch([11, 14], [13, 17], outlet),
        z_top_patch([18, 14], [20, 17], outlet),
    ];

    let spec = GlobalSpec {
        dims: IMPINGING_DIMS,
        nu: 0.04,
        periodic: [false, false, false],
        faces,
        face_patches,
        ..Default::default()
    };
    let walls = walls_for_faces(&[Face::XNeg, Face::XPos, Face::YNeg, Face::YPos, Face::ZNeg]);
    (spec, walls)
}

fn radial_wall_profile(s: &Cpu3<LocalPeriodic>, radii: &[usize]) -> Vec<f64> {
    let [nx, ny, _nz] = IMPINGING_DIMS;
    let cx = (nx - 1) as f64 * 0.5;
    let cy = (ny - 1) as f64 * 0.5;
    let z = 1usize;
    radii
        .iter()
        .map(|&r| {
            let mut values = Vec::new();
            for y in 1..(ny - 1) {
                for x in 1..(nx - 1) {
                    let dx = x as f64 - cx;
                    let dy = y as f64 - cy;
                    let rr = (dx * dx + dy * dy).sqrt();
                    if (rr - r as f64).abs() <= 0.5 {
                        let u = s.u(x, y, z);
                        let ur = (u[0] * dx + u[1] * dy) / rr.max(1.0e-12);
                        values.push(ur);
                    }
                }
            }
            assert!(
                !values.is_empty(),
                "radial bin r={r} unexpectedly empty for dims={IMPINGING_DIMS:?}"
            );
            values.iter().sum::<f64>() / values.len() as f64
        })
        .collect()
}

fn max_field_delta<HA: HaloExchange<f64>, HB: HaloExchange<f64>>(
    a: &Cpu3<HA>,
    b: &Cpu3<HB>,
) -> f64 {
    let mut max = 0.0f64;
    for (va, vb) in [
        (a.gather_rho(), b.gather_rho()),
        (a.gather_ux(), b.gather_ux()),
        (a.gather_uy(), b.gather_uy()),
        (a.gather_uz(), b.gather_uz()),
    ] {
        max = max.max(
            va.iter()
                .zip(&vb)
                .map(|(x, y)| (x - y).abs())
                .fold(0.0, f64::max),
        );
    }
    for q in 0..D3Q19::Q {
        max = max.max(
            a.gather_f(q)
                .iter()
                .zip(&b.gather_f(q))
                .map(|(x, y)| (x - y).abs())
                .fold(0.0, f64::max),
        );
    }
    max
}

#[test]
fn t18_2_impinging_jet_same_top_face_conserves_mass_and_has_radial_wall_jet() {
    let (spec, walls) = impinging_jet_spec();
    assert!(
        spec.validate(3, &[]).is_ok(),
        "impinging jet masked-face spec must validate: {:?}",
        spec.validate(3, &[])
    );

    let mut s = build(&spec, &walls, [1, 1, 1], LocalPeriodic);
    s.run(IMPINGING_STEPS);
    let m0 = s.total_mass();
    s.run(200);
    let m1 = s.total_mass();
    let drift = (m1 - m0).abs() / m0.abs();
    assert!(
        drift <= IMPINGING_MASS_DRIFT_REL_BAND,
        "masked impinging jet steady mass drift rel={drift:e}, band={IMPINGING_MASS_DRIFT_REL_BAND:e}, m0={m0:.12e}, m1={m1:.12e}"
    );

    let radii = [0usize, 2, 4, 6, 8, 10, 12];
    let profile = radial_wall_profile(&s, &radii);
    let axis = profile[0].abs();
    let peak_idx = profile
        .iter()
        .enumerate()
        .skip(1)
        .max_by(|(_, a), (_, b)| a.total_cmp(b))
        .map(|(i, _)| i)
        .unwrap();
    let peak = profile[peak_idx];
    assert!(
        axis <= 0.25 * peak.abs(),
        "floor wall-jet axis is not a stagnation minimum: radii={radii:?}, profile={profile:?}, axis_abs={axis:e}, peak={peak:e}"
    );
    assert!(
        peak_idx > 0 && peak > 1.0e-4,
        "floor wall-jet lacks a positive off-axis peak: radii={radii:?}, profile={profile:?}, peak_idx={peak_idx}, peak={peak:e}"
    );
    for i in (peak_idx + 1)..profile.len() {
        assert!(
            profile[i] <= profile[i - 1] + 5.0e-5,
            "floor wall-jet does not decay monotonically after peak: radii={radii:?}, profile={profile:?}, violation_pair=({}, {})",
            radii[i - 1],
            radii[i]
        );
    }
}

#[test]
fn t18_2_rejects_patch_out_of_face_bounds() {
    let spec = GlobalSpec {
        dims: [12, 10, 8],
        nu: 0.04,
        periodic: [false, false, false],
        faces: closed_faces(),
        face_patches: vec![z_top_patch([0, 0], [12, 9], FaceBC::Outflow)],
        ..Default::default()
    };
    let err = spec.validate(3, &[]);
    assert!(
        err.is_err(),
        "out-of-bounds top-face patch must be rejected, got {err:?}"
    );
}

#[test]
fn t18_2_rejects_overlapping_patches_on_same_face() {
    let spec = GlobalSpec {
        dims: [16, 16, 8],
        nu: 0.04,
        periodic: [false, false, false],
        faces: closed_faces(),
        face_patches: vec![
            z_top_patch([4, 4], [8, 8], FaceBC::Outflow),
            z_top_patch([8, 6], [12, 10], FaceBC::Pressure { rho: 1.0 }),
        ],
        ..Default::default()
    };
    let err = spec.validate(3, &[]);
    assert!(
        err.is_err(),
        "overlapping top-face patches must be rejected, got {err:?}"
    );
}

#[test]
fn t18_2_rejects_base_patch_union_open_on_more_than_one_axis() {
    let mut faces = closed_faces();
    faces[Face::XPos.index()] = FaceBC::Outflow;
    let spec = GlobalSpec {
        dims: [16, 16, 8],
        nu: 0.04,
        periodic: [false, false, false],
        faces,
        face_patches: vec![z_top_patch(
            [6, 6],
            [9, 9],
            FaceBC::Velocity {
                u: [0.0, 0.0, -0.03],
            },
        )],
        ..Default::default()
    };
    let err = spec.validate(3, &[]);
    assert!(
        err.is_err(),
        "base+patch union with open BCs on x and z axes must be rejected, got {err:?}"
    );
}

#[test]
fn t18_2_rejects_velocity_patch_above_max_speed() {
    let spec = GlobalSpec {
        dims: [16, 16, 8],
        nu: 0.04,
        periodic: [false, false, false],
        faces: closed_faces(),
        face_patches: vec![z_top_patch(
            [6, 6],
            [9, 9],
            FaceBC::Velocity {
                u: [0.0, 0.0, -(MAX_SPEED + 1.0e-12)],
            },
        )],
        ..Default::default()
    };
    let err = spec.validate(3, &[]);
    assert!(
        matches!(err, Err(SpecError::VelocityTooHigh { speed }) if speed > MAX_SPEED),
        "velocity patch above MAX_SPEED must return VelocityTooHigh, got {err:?}"
    );
}

#[test]
fn t18_2_full_face_patch_matches_base_face_bc_within_roundoff() {
    let dims = [18usize, 12usize, 10usize];
    let mut base_faces = closed_faces();
    base_faces[Face::XNeg.index()] = FaceBC::Velocity {
        u: [0.035, 0.0, 0.0],
    };
    base_faces[Face::XPos.index()] = FaceBC::Pressure { rho: 1.0 };
    let base = GlobalSpec {
        dims,
        nu: 0.04,
        periodic: [false, false, false],
        faces: base_faces,
        ..Default::default()
    };

    let mut patch_faces = closed_faces();
    patch_faces[Face::XPos.index()] = FaceBC::Pressure { rho: 1.0 };
    let patched = GlobalSpec {
        dims,
        nu: 0.04,
        periodic: [false, false, false],
        faces: patch_faces,
        face_patches: vec![x_neg_patch(
            [0, 0],
            [dims[1] - 1, dims[2] - 1],
            FaceBC::Velocity {
                u: [0.035, 0.0, 0.0],
            },
        )],
        ..Default::default()
    };
    let walls = walls_for_faces(&[Face::YNeg, Face::YPos, Face::ZNeg, Face::ZPos]);

    let mut a = build(&base, &walls, [1, 1, 1], LocalPeriodic);
    let mut b = build(&patched, &walls, [1, 1, 1], LocalPeriodic);
    a.run(80);
    b.run(80);
    let max_delta = max_field_delta(&a, &b);
    assert!(
        max_delta <= EQUIV_FIELD_ABS_BAND,
        "full-face patch/base equivalence max field delta={max_delta:e}, band={EQUIV_FIELD_ABS_BAND:e}"
    );
}

#[test]
fn t18_2_patch_straddling_inprocess_seam_is_bit_exact() {
    let dims = [20usize, 20usize, 12usize];
    let mut faces = closed_faces();
    faces[Face::ZPos.index()] = FaceBC::Closed;
    let spec = GlobalSpec {
        dims,
        nu: 0.04,
        periodic: [false, false, false],
        faces,
        face_patches: vec![
            z_top_patch(
                [8, 8],
                [12, 12],
                FaceBC::Velocity {
                    u: [0.0, 0.0, -0.035],
                },
            ),
            z_top_patch([13, 8], [16, 12], FaceBC::Pressure { rho: 1.0 }),
        ],
        ..Default::default()
    };
    let walls = walls_for_faces(&[Face::XNeg, Face::XPos, Face::YNeg, Face::YPos, Face::ZNeg]);

    let mut mono = build(&spec, &walls, [1, 1, 1], LocalPeriodic);
    let mut split = build(&spec, &walls, [2, 2, 2], InProcess);
    for t in 1..=120 {
        mono.step();
        split.step();
        if t <= 3 || t % 30 == 0 {
            let max_delta = max_field_delta(&mono, &split);
            assert_eq!(
                max_delta, 0.0,
                "masked patch straddling 2x2x2 seam differs at t={t}: max_delta={max_delta:e}"
            );
        }
    }
}

#[cfg(feature = "gpu")]
#[test]
fn t18_2_gpu_rejects_nonempty_face_patches() {
    use std::sync::{Arc, OnceLock};

    fn ctx() -> Arc<GpuContext> {
        static CTX: OnceLock<Arc<GpuContext>> = OnceLock::new();
        CTX.get_or_init(|| GpuContext::new().expect("T18.2 GPU rejection requires an adapter"))
            .clone()
    }

    let spec = GlobalSpec::<f32> {
        dims: [12, 12, 8],
        nu: 0.04,
        periodic: [false, false, false],
        faces: [FaceBC::Closed; 6],
        face_patches: vec![FacePatch {
            face: Face::ZPos.index(),
            lo: [4, 4],
            hi: [7, 7],
            bc: FaceBC::Velocity {
                u: [0.0, 0.0, -0.03],
            },
        }],
        ..Default::default()
    };
    let err = GpuSolver::<D3Q19>::try_new(&spec, &[], &[], ctx())
        .expect_err("GPU solver construction must reject non-empty face_patches");
    assert!(
        err.to_string().contains("SpecError") || err.to_string().contains("face_patches"),
        "GPU face_patches rejection must be surfaced as a SpecError, got {err:?}"
    );
}
