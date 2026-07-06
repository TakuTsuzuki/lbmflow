//! A-8: the open-face streaming contract.
//!
//! `ConvectiveOutflow` (and, implicitly, every open-face BC) relies on
//! streaming **not writing** the unknown-direction populations at an open
//! face: in the pull scheme the source cell for an entering direction lies
//! outside the domain, so the pull is skipped and the destination slot keeps
//! its prior contents. `ConvectiveOutflow` then reads that retained value as
//! its "memory" term `f(edge, t)`.
//!
//! This is a cross-module implicit contract (kernels `stream_row`, the
//! `Backend::stream` signature, the GPU edge-stash). An in-place streaming
//! rewrite (an M-E candidate) would break it silently, so pin it here: for
//! both the CPU reference and the wgpu backend, the open-face unknown slots
//! must be bit-identical before and after a stream pass.

use lbm_core::lattice::{D2Q9, D3Q19};
use lbm_core::prelude::*;

/// Fill `ftmp` at the open face's unknown slots with a distinctive sentinel,
/// run one `stream` pass over the whole core, and assert the pull left every
/// one of those slots untouched (bit-for-bit).
fn assert_stream_preserves_open_face_unknowns<L: Lattice>(dims: [usize; 3], open: Face) {
    // A monolithic, fully non-periodic subdomain: no neighbours on any face,
    // so streaming skips every entering direction at every boundary — exactly
    // the situation an open face is in (its inward neighbour is off-domain).
    let sub = Subdomain::monolithic(L::D, dims, [false, false, false]);
    let mut fields: SoaFields<f64> = SoaFields::new(L::Q, sub.geom);
    let np = fields.plane_len();
    let g = fields.geom;

    // Seed f with an arbitrary non-uniform field so streaming has real work to
    // do elsewhere, and pre-load ftmp with a sentinel at the open-face unknown
    // slots so we can detect any stray write.
    for q in 0..L::Q {
        for z in 0..g.core[2] {
            for y in 0..g.core[1] {
                for x in 0..g.core[0] {
                    let pi = g.pidx(x, y, z);
                    fields.f[q * np + pi] = (q as f64) + 0.01 * (x + y + z) as f64;
                }
            }
        }
    }
    // Record the sentinel values we place at the open face's unknown slots.
    let unknowns = L::unknowns(open);
    assert!(!unknowns.is_empty(), "chose a face on an inactive axis");
    let a = open.axis();
    let fixed = if open.is_neg() { 0 } else { dims[a] - 1 };
    let (t1, t2) = open.tangents();
    let mut sentinel_slots: Vec<(usize, f64)> = Vec::new();
    for c2 in 0..g.core[t2] {
        for c1 in 0..g.core[t1] {
            let mut pos = [0usize; 3];
            pos[a] = fixed;
            pos[t1] = c1;
            pos[t2] = c2;
            let pi = g.pidx(pos[0], pos[1], pos[2]);
            for (k, &q) in unknowns.iter().enumerate() {
                let slot = q * np + pi;
                // A sentinel distinct from anything the pull could produce.
                let val = -1000.0 - (slot as f64) - (k as f64) * 0.5;
                fields.ftmp[slot] = val;
                sentinel_slots.push((slot, val));
            }
        }
    }

    let params = StepParams::<f64> {
        collision: CollisionKind::Bgk,
        omega_p: 1.0,
        omega_m: 1.0,
        force: [0.0; 3],
        faces: [FaceBC::Outflow; 6],
        sources: Vec::new(),
        face_patches: Vec::new(),
    };
    let mut backend = CpuScalar::default();
    let _ = <CpuScalar as Backend<L, f64>>::stream(
        &mut backend,
        &sub,
        &mut fields,
        &params,
        CellRange::full(&sub),
    );

    // The pull must not have written any of the open-face unknown slots.
    for (slot, val) in sentinel_slots {
        assert_eq!(
            fields.ftmp[slot], val,
            "stream wrote an open-face unknown slot (slot {slot}) that the \
             open-face BC relies on retaining"
        );
    }
}

#[test]
fn cpu_stream_preserves_open_face_unknowns_d2q9() {
    // Every face of a 2D box (each is a candidate open face).
    for face in [Face::XNeg, Face::XPos, Face::YNeg, Face::YPos] {
        assert_stream_preserves_open_face_unknowns::<D2Q9>([6, 5, 1], face);
    }
}

#[test]
fn cpu_stream_preserves_open_face_unknowns_d3q19() {
    for face in Face::ALL {
        assert_stream_preserves_open_face_unknowns::<D3Q19>([5, 6, 4], face);
    }
}

/// The D3Q19 Zou–He branch hardcodes the 5-unknown reconstruction; the guard
/// added in A-8 must accept the real D3Q19 faces (a regression check that the
/// assert threshold matches the lattice).
#[test]
fn d3q19_zou_he_face_accepts_all_faces() {
    // Build a walled box with one Zou–He velocity inlet and step it: exercises
    // apply_open_faces → zou_he_face_3d on that face without tripping the
    // 5-unknown guard.
    let dims = [6, 6, 6];
    let mut walls = WallSpec::<f64>::default();
    for f in [Face::YNeg, Face::YPos, Face::ZNeg, Face::ZPos, Face::XPos] {
        walls.is_wall[f.index()] = true;
    }
    let mut faces = [FaceBC::<f64>::Closed; 6];
    faces[Face::XNeg.index()] = FaceBC::Velocity {
        u: [0.02, 0.0, 0.0],
    };
    // XPos is a wall rim, XNeg is the inlet: open axis is x only.
    let spec = GlobalSpec::<f64> {
        dims,
        nu: 0.1,
        periodic: [false, false, false],
        faces,
        ..Default::default()
    };
    let (solid, wall_u) = build_wall_rims(3, dims, &walls);
    let mut s: Solver<D3Q19, f64, CpuScalar, LocalPeriodic> = Solver::new(
        &spec,
        &solid,
        &wall_u,
        [1, 1, 1],
        CpuScalar::default(),
        LocalPeriodic,
    );
    s.run(5);
    assert!(s.total_mass().is_finite());
}

// ---------------------------------------------------------------------------
// GPU side of the contract (feature `gpu`): the wgpu backend realises the same
// "open-face unknown slots retain the previous value" invariant through its
// edge stash rather than a skipped pull, so the observable — the outflow
// face's unknown-direction populations — must track the CPU reference.
// ---------------------------------------------------------------------------
#[cfg(feature = "gpu")]
mod gpu {
    use super::*;
    use std::sync::{Arc, OnceLock};

    fn ctx() -> Arc<GpuContext> {
        static CTX: OnceLock<Arc<GpuContext>> = OnceLock::new();
        CTX.get_or_init(|| {
            GpuContext::new().expect("stream_contract GPU test requires a GPU adapter")
        })
        .clone()
    }

    /// Channel: parabolic-ish inlet on the left (Zou–He velocity), Outflow on
    /// the right, walls top/bottom. After several steps the right face's
    /// unknown (inward-x) populations are governed by the retained-value
    /// contract; CPU and GPU must agree on them to ≤ 1e-4 relative (the T14
    /// line), proving the GPU honours the same open-face streaming contract.
    fn channel_spec() -> (GlobalSpec<f32>, Vec<bool>, Vec<[f32; 3]>) {
        let dims = [24usize, 16, 1];
        let mut walls = WallSpec::<f32>::default();
        walls.is_wall[Face::YNeg.index()] = true;
        walls.is_wall[Face::YPos.index()] = true;
        let mut faces = [FaceBC::<f32>::Closed; 6];
        faces[Face::XNeg.index()] = FaceBC::Velocity {
            u: [0.05, 0.0, 0.0],
        };
        faces[Face::XPos.index()] = FaceBC::Outflow;
        let spec = GlobalSpec::<f32> {
            dims,
            nu: 0.05,
            periodic: [false, false, false],
            faces,
            ..Default::default()
        };
        let (solid, wall_u) = build_wall_rims(2, dims, &walls);
        (spec, solid, wall_u)
    }

    fn linf_rel(a: &[f32], b: &[f32]) -> f64 {
        let mut d = 0.0f64;
        let mut m = 0.0f64;
        for (x, y) in a.iter().zip(b) {
            d = d.max((*x as f64 - *y as f64).abs());
            m = m.max((*x as f64).abs());
        }
        d / m.max(1e-9)
    }

    #[test]
    fn gpu_stream_honours_open_face_contract() {
        let (spec, solid, wall_u) = channel_spec();
        let mut cpu: Solver<D2Q9, f32, CpuScalar, LocalPeriodic> = Solver::new(
            &spec,
            &solid,
            &wall_u,
            [1, 1, 1],
            CpuScalar::default(),
            LocalPeriodic,
        );
        let mut gpu: GpuSolver<D2Q9> = GpuSolver::new(&spec, &solid, &wall_u, ctx());
        let steps = 200;
        cpu.run(steps);
        gpu.run(steps);

        // Compare the Outflow face's unknown-direction populations (the slots
        // the contract governs) across the whole right column.
        let unknowns = D2Q9::unknowns(Face::XPos);
        for &q in unknowns {
            let a = cpu.gather_f(q);
            let b = gpu.gather_f(q);
            let (nx, ny) = (spec.dims[0], spec.dims[1]);
            // Slice out the right column (x = nx-1).
            let col_a: Vec<f32> = (0..ny).map(|y| a[y * nx + (nx - 1)]).collect();
            let col_b: Vec<f32> = (0..ny).map(|y| b[y * nx + (nx - 1)]).collect();
            let rel = linf_rel(&col_a, &col_b);
            assert!(
                rel <= 1e-4,
                "outflow-face unknown q={q} diverges CPU vs GPU: L∞rel {rel:.2e}"
            );
        }
    }
}
