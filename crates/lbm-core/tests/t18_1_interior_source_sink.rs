//! T18.1 adversarial acceptance tests for CR-1 localized interior
//! volume sources/sinks.
//!
//! Authored from docs/VALIDATION.md T18.1 and docs/DISPERSED_DEPOSITION.md
//! CR-1 only. The CR-1 public API is expected to be absent in this worktree;
//! this file should compile only after the implementation branch lands.

use lbm_core::lattice::D3Q19;
use lbm_core::params::{SourceKind, SourceRegion, VolumeSource, MAX_SPEED};
use lbm_core::prelude::*;

type Sol<H> = Solver<D3Q19, f64, CpuScalar, H>;

fn closed_box_spec(dims: [usize; 3], nu: f64, sources: Vec<VolumeSource<f64>>) -> GlobalSpec<f64> {
    GlobalSpec {
        dims,
        nu,
        periodic: [false, false, false],
        faces: [FaceBC::Closed; 6],
        sources,
        ..Default::default()
    }
}

fn all_walls() -> WallSpec<f64> {
    let mut walls = WallSpec::default();
    for face in Face::ALL {
        walls.is_wall[face.index()] = true;
    }
    walls
}

fn build<H: HaloExchange<f64>>(spec: &GlobalSpec<f64>, decomp: [usize; 3], exchange: H) -> Sol<H> {
    let (solid, wall_u) = build_wall_rims(3, spec.dims, &all_walls());
    Solver::new(
        spec,
        &solid,
        &wall_u,
        decomp,
        CpuScalar::default(),
        exchange,
    )
}

fn source(region: SourceRegion, kind: SourceKind<f64>) -> VolumeSource<f64> {
    VolumeSource { region, kind }
}

fn mass_flow(lo: [usize; 3], hi: [usize; 3], q_lu: f64) -> VolumeSource<f64> {
    source(SourceRegion { lo, hi }, SourceKind::MassFlow { q_lu })
}

fn jet(lo: [usize; 3], hi: [usize; 3], q_lu: f64, u: [f64; 3]) -> VolumeSource<f64> {
    source(SourceRegion { lo, hi }, SourceKind::Jet { q_lu, u })
}

fn assert_spec_err(spec: &GlobalSpec<f64>, solid: &[bool], label: &str) {
    let result = spec.validate(3, solid);
    assert!(result.is_err(), "{label}: expected Err(SpecError), got Ok");
}

fn cell_index(dims: [usize; 3], x: usize, y: usize, z: usize) -> usize {
    z * dims[0] * dims[1] + y * dims[0] + x
}

fn radial_profile_error(
    s: &Sol<LocalPeriodic>,
    dims: [usize; 3],
    q_abs: f64,
    center: [f64; 3],
    r_min: f64,
    r_max: f64,
) -> f64 {
    let ux = s.gather_ux();
    let uy = s.gather_uy();
    let uz = s.gather_uz();
    let mut sum_inward_ur = 0.0;
    let mut sum_expected = 0.0;
    let mut n = 0usize;
    for z in 1..dims[2] - 1 {
        for y in 1..dims[1] - 1 {
            for x in 1..dims[0] - 1 {
                let dx = x as f64 - center[0];
                let dy = y as f64 - center[1];
                let dz = z as f64 - center[2];
                let r = dx.hypot(dy).hypot(dz);
                if r < r_min || r > r_max {
                    continue;
                }
                let i = cell_index(dims, x, y, z);
                let outward_ur = (ux[i] * dx + uy[i] * dy + uz[i] * dz) / r;
                sum_inward_ur += -outward_ur;
                sum_expected += q_abs / (4.0 * std::f64::consts::PI * r * r);
                n += 1;
            }
        }
    }
    assert!(n > 128, "radial shell sample count too small: n={n}");
    let measured = sum_inward_ur / n as f64;
    let expected = sum_expected / n as f64;
    let rel = (measured - expected).abs() / expected.abs();
    assert!(
        rel <= 0.10,
        "sink radial far-field rel={rel:.3e}, measured={measured:.6e}, expected={expected:.6e}, samples={n}"
    );
    rel
}

fn assert_fields_bit_equal<HA: HaloExchange<f64>, HB: HaloExchange<f64>>(
    a: &Sol<HA>,
    b: &Sol<HB>,
    label: &str,
) {
    for (name, va, vb) in [
        ("rho", a.gather_rho(), b.gather_rho()),
        ("ux", a.gather_ux(), b.gather_ux()),
        ("uy", a.gather_uy(), b.gather_uy()),
        ("uz", a.gather_uz(), b.gather_uz()),
    ] {
        let d = va
            .iter()
            .zip(&vb)
            .map(|(x, y)| (x - y).abs())
            .fold(0.0f64, f64::max);
        assert_eq!(d, 0.0, "{label}: field {name} max|delta|={d:e}");
    }
    for q in 0..D3Q19::Q {
        let fa = a.gather_f(q);
        let fb = b.gather_f(q);
        let d = fa
            .iter()
            .zip(&fb)
            .map(|(x, y)| (x - y).abs())
            .fold(0.0f64, f64::max);
        assert_eq!(d, 0.0, "{label}: f[{q}] max|delta|={d:e}");
    }
}

#[test]
fn t18_1_sink_far_field_default_shell_profile() {
    let q_abs = 2.0e-5;
    let dims = [32, 32, 32];
    let region = SourceRegion {
        lo: [15, 15, 15],
        hi: [16, 16, 16],
    };
    let mut s = build(
        &closed_box_spec(
            dims,
            0.04,
            vec![source(region, SourceKind::MassFlow { q_lu: -q_abs })],
        ),
        [1, 1, 1],
        LocalPeriodic,
    );
    s.run(360);
    // Shell band: r in [5.5, 10.5], excluding the 2^3 source near-field and
    // leaving at least 5 cells to the wall rim. This band is frozen at first
    // measurement and checks the incompressible sink law u_r=q/(4*pi*r^2).
    let rel = radial_profile_error(&s, dims, q_abs, [15.5, 15.5, 15.5], 5.5, 10.5);
    assert!(
        rel <= 0.10,
        "default sink far-field rel={rel:.3e} after 360 steps"
    );
}

#[test]
#[ignore = "64^3 quasi-steady T18.1 far-field acceptance case"]
fn t18_1_sink_far_field_64_class_ignored() {
    let q_abs = 8.0e-5;
    let dims = [64, 64, 64];
    let region = SourceRegion {
        lo: [31, 31, 31],
        hi: [32, 32, 32],
    };
    let mut s = build(
        &closed_box_spec(
            dims,
            0.04,
            vec![source(region, SourceKind::MassFlow { q_lu: -q_abs })],
        ),
        [1, 1, 1],
        LocalPeriodic,
    );
    s.run(1200);
    // Shell band: r in [9.5, 22.5], well outside the 2^3 source stencil and
    // still at least 9 cells from every wall. Frozen at first measurement.
    let rel = radial_profile_error(&s, dims, q_abs, [31.5, 31.5, 31.5], 9.5, 22.5);
    assert!(
        rel <= 0.10,
        "64^3 sink far-field rel={rel:.3e} after 1200 steps"
    );
}

#[test]
fn t18_1_mass_ledger_matches_sum_q_per_step() {
    let q_sum = 3.0e-6;
    let spec = closed_box_spec(
        [18, 18, 18],
        0.05,
        vec![mass_flow([8, 8, 8], [9, 9, 9], q_sum)],
    );
    let mut s = build(&spec, [1, 1, 1], LocalPeriodic);
    let mut m0 = s.total_mass();
    for step in 1..=32 {
        s.step();
        let m1 = s.total_mass();
        let dm = m1 - m0;
        let rel = (dm - q_sum).abs() / q_sum.abs();
        assert!(
            rel <= 1.0e-12,
            "mass ledger step={step}: rel={rel:.3e}, dm={dm:.12e}, sum_q={q_sum:.12e}"
        );
        m0 = m1;
    }
}

#[test]
fn t18_1_jet_delivers_prescribed_momentum_flux() {
    let q_lu = 2.0e-5;
    let u = [0.04, 0.012, 0.0];
    let spec = closed_box_spec(
        [24, 24, 24],
        0.04,
        vec![jet([11, 11, 11], [12, 12, 12], q_lu, u)],
    );
    let mut s = build(&spec, [1, 1, 1], LocalPeriodic);
    let p0 = s.total_momentum();
    s.run(24);
    let p1 = s.total_momentum();
    for axis in 0..2 {
        let measured = (p1[axis] - p0[axis]) / 24.0;
        let expected = q_lu * u[axis];
        let rel = (measured - expected).abs() / expected.abs();
        assert!(
            rel <= 0.02,
            "jet momentum axis={axis}: rel={rel:.3e}, measured_flux={measured:.12e}, expected_flux={expected:.12e}"
        );
    }
    let z_flux = (p1[2] - p0[2]) / 24.0;
    assert!(
        z_flux.abs() <= 1.0e-14,
        "jet momentum z leakage measured_flux={z_flux:.12e}"
    );
}

#[test]
fn t18_1_rejects_region_touching_face() {
    let spec = closed_box_spec(
        [12, 12, 12],
        0.05,
        vec![mass_flow([0, 5, 5], [1, 6, 6], 1.0e-6)],
    );
    assert_spec_err(&spec, &[], "source region touching x-neg face");
}

#[test]
fn t18_1_rejects_overlapping_sources() {
    let spec = closed_box_spec(
        [14, 14, 14],
        0.05,
        vec![
            mass_flow([5, 5, 5], [7, 7, 7], 1.0e-6),
            mass_flow([7, 7, 7], [8, 8, 8], -1.0e-6),
        ],
    );
    assert_spec_err(&spec, &[], "inclusive source-region overlap at [7,7,7]");
}

#[test]
fn t18_1_rejects_source_overlapping_solid() {
    let dims = [14, 14, 14];
    let spec = closed_box_spec(dims, 0.05, vec![mass_flow([6, 6, 6], [7, 7, 7], 1.0e-6)]);
    let mut solid = vec![false; dims[0] * dims[1] * dims[2]];
    solid[cell_index(dims, 6, 6, 6)] = true;
    assert_spec_err(
        &spec,
        &solid,
        "source region overlapping build-time solid mask",
    );
}

#[test]
fn t18_1_rejects_jet_speed_above_max_speed() {
    let spec = closed_box_spec(
        [14, 14, 14],
        0.05,
        vec![jet(
            [6, 6, 6],
            [7, 7, 7],
            1.0e-6,
            [MAX_SPEED + 1.0e-12, 0.0, 0.0],
        )],
    );
    assert_spec_err(&spec, &[], "jet speed above MAX_SPEED");
}

#[test]
fn t18_1_rejects_sink_that_can_drive_local_density_nonpositive() {
    let spec = closed_box_spec(
        [14, 14, 14],
        0.05,
        vec![mass_flow([6, 6, 6], [7, 7, 7], -9.0)],
    );
    assert_spec_err(
        &spec,
        &[],
        "sink draining more than one rho unit per source cell per step",
    );
}

#[test]
fn t18_1_partition_invariance_source_straddles_2x2x2_seams() {
    let spec = closed_box_spec(
        [16, 16, 16],
        0.04,
        vec![mass_flow([7, 7, 7], [8, 8, 8], 1.0e-6)],
    );
    let mut base = build(&spec, [1, 1, 1], LocalPeriodic);
    let mut split = build(&spec, [2, 2, 2], InProcess);
    assert_fields_bit_equal(&base, &split, "t=0");
    for t in 1..=80 {
        base.step();
        split.step();
        if t <= 3 || t % 20 == 0 {
            assert_fields_bit_equal(&base, &split, &format!("t={t}"));
        }
    }
}

#[cfg(feature = "gpu")]
#[test]
fn t18_1_gpu_rejects_nonempty_sources_with_spec_error() {
    use std::sync::{Arc, OnceLock};

    fn ctx() -> Arc<GpuContext> {
        static CTX: OnceLock<Arc<GpuContext>> = OnceLock::new();
        CTX.get_or_init(|| {
            GpuContext::new().expect("T18.1 GPU rejection test requires a GPU adapter")
        })
        .clone()
    }

    let spec = closed_box_spec(
        [16, 16, 16],
        0.05,
        vec![mass_flow([7, 7, 7], [8, 8, 8], 1.0e-6)],
    );
    let (solid, wall_u) = build_wall_rims(3, spec.dims, &all_walls());
    let result = GpuSolver::<D3Q19>::build(&spec, &solid, &wall_u, ctx());
    assert!(
        matches!(result, Err(_)),
        "GPU solver must reject non-empty sources with Err(SpecError), got Ok"
    );
}
