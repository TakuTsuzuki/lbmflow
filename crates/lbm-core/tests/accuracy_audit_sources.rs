//! ACC-AUDIT: adversarial source and masked-face probes for CR-1/CR-2.
//!
//! These tests are authored from the public V2 API plus continuum/discrete
//! conservation laws only. They intentionally sharpen the broader T18 bands:
//! exact step ledgers, symmetry maps, far-field scaling, and per-node BC
//! exactness.

mod common;

use common::metrics::{linf_rel, monotonicity, order_fit};
use lbm_core::lattice::D3Q19;
use lbm_core::params::{FacePatch, SourceKind, SourceRegion, VolumeSource, MAX_SPEED};
use lbm_core::prelude::*;

type Cpu3 = Solver<D3Q19, f64, CpuScalar, LocalPeriodic>;

fn closed_faces() -> [FaceBC<f64>; 6] {
    [FaceBC::Closed; 6]
}

fn all_walls() -> WallSpec<f64> {
    let mut walls = WallSpec::default();
    for face in Face::ALL {
        walls.is_wall[face.index()] = true;
    }
    walls
}

fn walls_for_faces(faces: &[Face]) -> WallSpec<f64> {
    let mut walls = WallSpec::default();
    for &face in faces {
        walls.is_wall[face.index()] = true;
    }
    walls
}

fn build(spec: &GlobalSpec<f64>, walls: &WallSpec<f64>) -> Cpu3 {
    let (solid, wall_u) = build_wall_rims(3, spec.dims, walls);
    Solver::new(
        spec,
        &solid,
        &wall_u,
        [1, 1, 1],
        CpuScalar::default(),
        LocalPeriodic,
    )
}

fn idx(dims: [usize; 3], x: usize, y: usize, z: usize) -> usize {
    z * dims[0] * dims[1] + y * dims[0] + x
}

fn volume_source(lo: [usize; 3], hi: [usize; 3], kind: SourceKind<f64>) -> VolumeSource<f64> {
    VolumeSource {
        region: SourceRegion { lo, hi },
        kind,
    }
}

fn region_cell_count(lo: [usize; 3], hi: [usize; 3]) -> usize {
    (hi[0] - lo[0] + 1) * (hi[1] - lo[1] + 1) * (hi[2] - lo[2] + 1)
}

fn momentum_from_gather(s: &Cpu3, dims: [usize; 3]) -> [f64; 3] {
    let rho = s.gather_rho();
    let ux = s.gather_ux();
    let uy = s.gather_uy();
    let uz = s.gather_uz();
    let mut p = [0.0; 3];
    for z in 0..dims[2] {
        for y in 0..dims[1] {
            for x in 0..dims[0] {
                if s.is_solid(x, y, z) {
                    continue;
                }
                let i = idx(dims, x, y, z);
                p[0] += rho[i] * ux[i];
                p[1] += rho[i] * uy[i];
                p[2] += rho[i] * uz[i];
            }
        }
    }
    p
}

fn z_top_patch(lo: [usize; 2], hi: [usize; 2], bc: FaceBC<f64>) -> FacePatch<f64> {
    FacePatch {
        face: Face::ZPos.index(),
        lo,
        hi,
        bc,
    }
}

fn top_annulus_patches() -> Vec<FacePatch<f64>> {
    vec![
        z_top_patch(
            [10, 10],
            [13, 13],
            FaceBC::Velocity {
                u: [0.0, 0.0, -0.02],
            },
        ),
        z_top_patch([7, 7], [16, 9], FaceBC::Pressure { rho: 1.0 }),
        z_top_patch([7, 14], [16, 16], FaceBC::Pressure { rho: 1.0 }),
        z_top_patch([7, 10], [9, 13], FaceBC::Pressure { rho: 1.0 }),
        z_top_patch([14, 10], [16, 13], FaceBC::Pressure { rho: 1.0 }),
    ]
}

fn top_patch_spec(dims: [usize; 3], patches: Vec<FacePatch<f64>>) -> (GlobalSpec<f64>, WallSpec<f64>) {
    let mut faces = closed_faces();
    faces[Face::ZPos.index()] = FaceBC::Closed;
    let spec = GlobalSpec {
        dims,
        nu: 0.05,
        periodic: [false, false, false],
        faces,
        face_patches: patches,
        ..Default::default()
    };
    let walls = walls_for_faces(&[Face::XNeg, Face::XPos, Face::YNeg, Face::YPos, Face::ZNeg]);
    (spec, walls)
}

fn in_rect(p: [usize; 2], lo: [usize; 2], hi: [usize; 2]) -> bool {
    lo[0] <= p[0] && p[0] <= hi[0] && lo[1] <= p[1] && p[1] <= hi[1]
}

fn adjacent_to_rect(p: [usize; 2], lo: [usize; 2], hi: [usize; 2]) -> bool {
    let x_near = p[0] + 1 >= lo[0] && p[0] <= hi[0] + 1;
    let y_near = p[1] + 1 >= lo[1] && p[1] <= hi[1] + 1;
    x_near && y_near && !in_rect(p, lo, hi)
}

fn max_mirror_delta(a: &Cpu3, b: &Cpu3, dims: [usize; 3]) -> [f64; 4] {
    let ar = a.gather_rho();
    let ax = a.gather_ux();
    let ay = a.gather_uy();
    let az = a.gather_uz();
    let br = b.gather_rho();
    let bx = b.gather_ux();
    let by = b.gather_uy();
    let bz = b.gather_uz();
    let mut d = [0.0f64; 4];
    for z in 0..dims[2] {
        for y in 0..dims[1] {
            for x in 0..dims[0] {
                let ia = idx(dims, x, y, z);
                let ib = idx(dims, dims[0] - 1 - x, y, z);
                d[0] = d[0].max((ar[ia] - br[ib]).abs());
                d[1] = d[1].max((ax[ia] + bx[ib]).abs());
                d[2] = d[2].max((ay[ia] - by[ib]).abs());
                d[3] = d[3].max((az[ia] - bz[ib]).abs());
            }
        }
    }
    d
}

fn dipole_projection(dims: [usize; 3], s: &Cpu3, radii: &[usize]) -> Vec<f64> {
    let ux = s.gather_ux();
    let uy = s.gather_uy();
    let uz = s.gather_uz();
    let c = [
        (dims[0] - 1) as f64 * 0.5,
        (dims[1] - 1) as f64 * 0.5,
        (dims[2] - 1) as f64 * 0.5,
    ];
    radii
        .iter()
        .map(|&r0| {
            let mut num = 0.0;
            let mut den = 0.0;
            for z in 0..dims[2] {
                for y in 0..dims[1] {
                    for x in 0..dims[0] {
                        let dx = x as f64 - c[0];
                        let dy = y as f64 - c[1];
                        let dz = z as f64 - c[2];
                        let r = dx.hypot(dy).hypot(dz);
                        if (r - r0 as f64).abs() > 0.5 || r == 0.0 {
                            continue;
                        }
                        let i = idx(dims, x, y, z);
                        let ur = (ux[i] * dx + uy[i] * dy + uz[i] * dz) / r;
                        let cos_theta = dz / r;
                        num += ur * cos_theta;
                        den += cos_theta * cos_theta;
                    }
                }
            }
            assert!(den > 0.0, "empty dipole shell r={r0}, dims={dims:?}");
            (num / den).abs()
        })
        .collect()
}

fn dipole_spec(dims: [usize; 3], q: f64) -> GlobalSpec<f64> {
    let c = [dims[0] / 2, dims[1] / 2, dims[2] / 2];
    GlobalSpec {
        dims,
        nu: 0.05,
        periodic: [true, true, true],
        faces: closed_faces(),
        sources: vec![
            volume_source(
                [c[0] - 1, c[1] - 1, c[2] + 3],
                [c[0], c[1], c[2] + 4],
                SourceKind::MassFlow { q_lu: q },
            ),
            volume_source(
                [c[0] - 1, c[1] - 1, c[2] - 4],
                [c[0], c[1], c[2] - 3],
                SourceKind::MassFlow { q_lu: -q },
            ),
        ],
        ..Default::default()
    }
}

#[test]
fn d1_jet_per_step_momentum_ledger_identity() {
    let dims = [24, 24, 24];
    let lo = [10, 10, 10];
    let hi = [13, 13, 13];
    let q_lu = 2.0e-7;
    let u = [0.03, 0.01, -0.02];
    assert!(u.iter().map(|v| v * v).sum::<f64>().sqrt() < MAX_SPEED);
    let spec = GlobalSpec {
        dims,
        nu: 0.05,
        periodic: [false, false, false],
        faces: closed_faces(),
        sources: vec![volume_source(lo, hi, SourceKind::Jet { q_lu, u })],
        ..Default::default()
    };
    let mut s = build(&spec, &all_walls());
    // Semantics pin (triage 2026-07-06, ANOM-P4-005): q_lu is the REGION
    // TOTAL mass flow, not a per-cell rate — consistent with the T18.1 mass
    // ledger d(total_mass)/step = SUM over SOURCES of q_lu. The injected
    // momentum per step is therefore q_lu * u (first measurement: the
    // per-cell reading over-predicted by exactly N_region = 64, and the
    // region-total identity then holds to 12 printed digits). The region
    // cell count must NOT appear in J; keep it computed only to document
    // the distinction.
    let n = region_cell_count(lo, hi) as f64;
    let _ = n;
    let j = [q_lu * u[0], q_lu * u[1], q_lu * u[2]];
    let mut p0 = momentum_from_gather(&s, dims);
    // A Jet source adds equilibrium-shaped populations with zeroth moment q
    // and first moment q*u at each source cell. Summing the discrete first
    // moment over all fluid nodes before wall contact gives
    // P(t+1)-P(t)=sum_cells q_cell*u_jet. The public velocity convention adds
    // the Guo half-force shift for body forces, but these sources are not a
    // force field, so gather-summed rho*u_phys should equal the raw first
    // moment and show the clean recurrence, not a half-step offset.
    for step in 1..=8 {
        s.step();
        let p1 = momentum_from_gather(&s, dims);
        for a in 0..3 {
            let delta = p1[a] - p0[a];
            let den = p0[a].abs() + j[a].abs();
            let err = (delta - j[a]).abs();
            let half_err = (delta - 0.5 * j[a]).abs();
            println!(
                "ACC SRC D1: step={step} axis={a} delta={delta:.12e} expected={:.12e} abs_err={err:.3e} band={:.3e} den={den:.12e} half_step_err={half_err:.3e}",
                j[a],
                1.0e-13 * den
            );
            assert!(
                err <= 1.0e-13 * den,
                "ACC SRC D1 axis={a} step={step}: measured_delta={delta:.12e}, expected_J={:.12e}, abs_err={err:.3e}, band={:.3e}, denominator=|P(t)|+|J|={den:.12e}, half_step_err={half_err:.3e}",
                j[a],
                1.0e-13 * den
            );
        }
        p0 = p1;
    }
}

#[test]
fn d2_source_dipole_far_field_light_canary() {
    let dims = [32, 32, 32];
    let radii = [6usize, 9, 12];
    let spec = dipole_spec(dims, 1.0e-7);
    let mut s = build(&spec, &WallSpec::default());
    s.run(900);
    // For a source/sink pair separated along z, the leading incompressible
    // far field is a dipole: u_r(r,theta)=C*cos(theta)/r^3. The least-squares
    // shell projection m(r)=sum(u_r*cos(theta))/sum(cos^2(theta)) therefore
    // returns C/r^3 for a pure dipole; the constant and discrete shell area
    // cancel in the log-log slope.
    let m = dipole_projection(dims, &s, &radii);
    let h: Vec<f64> = radii.iter().map(|&r| 1.0 / r as f64).collect();
    let fit = order_fit(&h, &m);
    let mono = monotonicity(&m);
    println!(
        "ACC SRC D2 light: radii={radii:?} projected={m:?} slope={:.6} r2={:.6} monotonicity={mono:.3}",
        fit.slope, fit.r2
    );
    assert!(
        mono == 1.0,
        "ACC SRC D2 light monotonicity={mono:.3}, band=1.0, projected={m:?}"
    );
    assert!(
        (2.0..=4.0).contains(&fit.slope) && fit.r2 >= 0.90,
        "ACC SRC D2 light slope={:.6}, band=[2.0,4.0], r2={:.6}, r2_band>=0.90, projected={m:?}, radii={radii:?}",
        fit.slope,
        fit.r2
    );
}

#[test]
#[ignore = "heavy ACC-AUDIT D2 dipole far-field 64^3"]
fn d2_source_dipole_far_field_full_64() {
    let dims = [64, 64, 64];
    let radii = [8usize, 11, 14, 17, 20];
    let spec = dipole_spec(dims, 1.0e-7);
    let mut s = build(&spec, &WallSpec::default());
    s.run(3000);
    // Same projection derivation as the light canary: any shellwise linear
    // projection of C*cos(theta)/r^3 onto cos(theta) inherits r^-3, so the
    // constant, source strength, and shell population do not enter the slope.
    let m = dipole_projection(dims, &s, &radii);
    let h: Vec<f64> = radii.iter().map(|&r| 1.0 / r as f64).collect();
    let fit = order_fit(&h, &m);
    println!(
        "ACC SRC D2 heavy: radii={radii:?} projected={m:?} slope={:.6} r2={:.6}",
        fit.slope, fit.r2
    );
    assert!(
        (fit.slope - 3.0).abs() <= 0.5 && fit.r2 >= 0.98,
        "ACC SRC D2 heavy slope={:.6}, band=3.0+/-0.5, r2={:.6}, r2_band>=0.98, projected={m:?}, radii={radii:?}",
        fit.slope,
        fit.r2
    );
}

#[test]
fn d3_sink_mirror_equivariance() {
    let dims = [20, 20, 20];
    let q_lu = -5.0e-8;
    let plus_lo = [14, 9, 9];
    let plus_hi = [15, 10, 10];
    // Mirroring an inclusive cell interval [lo,hi] on an N-cell axis maps
    // every x to N-1-x, so the mirrored interval is [N-1-hi, N-1-lo].
    let minus_lo = [dims[0] - 1 - plus_hi[0], plus_lo[1], plus_lo[2]];
    let minus_hi = [dims[0] - 1 - plus_lo[0], plus_hi[1], plus_hi[2]];
    let spec_plus = GlobalSpec {
        dims,
        nu: 0.05,
        periodic: [false, false, false],
        faces: closed_faces(),
        sources: vec![volume_source(
            plus_lo,
            plus_hi,
            SourceKind::MassFlow { q_lu },
        )],
        ..Default::default()
    };
    let spec_minus = GlobalSpec {
        sources: vec![volume_source(
            minus_lo,
            minus_hi,
            SourceKind::MassFlow { q_lu },
        )],
        ..spec_plus.clone()
    };
    let mut a = build(&spec_plus, &all_walls());
    let mut b = build(&spec_minus, &all_walls());
    a.run(200);
    b.run(200);
    let d = max_mirror_delta(&a, &b, dims);
    println!(
        "ACC SRC D3: max_delta rho={:.3e} ux_mirror={:.3e} uy={:.3e} uz={:.3e} band=1.000e-12",
        d[0], d[1], d[2], d[3]
    );
    assert!(
        d.iter().all(|v| *v <= 1.0e-12),
        "ACC SRC D3 mirror max deltas={d:?}, band=1e-12"
    );
}

#[test]
fn d4_per_cell_masked_patch_bc_node_exactness() {
    let dims = [24, 24, 24];
    let (spec, walls) = top_patch_spec(dims, top_annulus_patches());
    let mut s = build(&spec, &walls);
    s.run(50);
    let rho = s.gather_rho();
    let ux = s.gather_ux();
    let uy = s.gather_uy();
    let uz = s.gather_uz();
    let layer_error = |z: usize| {
        let mut vel_max = 0.0f64;
        let mut lid_max = 0.0f64;
        let mut pressure_max = 0.0f64;
        let mut seam_lid_max = 0.0f64;
        for y in 1..dims[1] - 1 {
            for x in 1..dims[0] - 1 {
                let p = [x, y];
                let i = idx(dims, x, y, z);
                let in_vel = in_rect(p, [10, 10], [13, 13]);
                let in_pressure = in_rect(p, [7, 7], [16, 9])
                    || in_rect(p, [7, 14], [16, 16])
                    || in_rect(p, [7, 10], [9, 13])
                    || in_rect(p, [14, 10], [16, 13]);
                if in_vel {
                    vel_max = vel_max
                        .max(ux[i].abs())
                        .max(uy[i].abs())
                        .max((uz[i] + 0.02).abs());
                } else if in_pressure {
                    pressure_max = pressure_max.max((rho[i] - 1.0).abs());
                } else {
                    let u_abs = ux[i].abs().max(uy[i].abs()).max(uz[i].abs());
                    lid_max = lid_max.max(u_abs);
                    if adjacent_to_rect(p, [7, 7], [16, 9])
                        || adjacent_to_rect(p, [7, 14], [16, 16])
                        || adjacent_to_rect(p, [7, 10], [9, 13])
                        || adjacent_to_rect(p, [14, 10], [16, 13])
                    {
                        seam_lid_max = seam_lid_max.max(u_abs);
                    }
                }
            }
        }
        (vel_max, lid_max, pressure_max, seam_lid_max)
    };
    // Convention pin (triage 2026-07-06, ANOM-P4-006): the patch BC nodes
    // are the FACE LAYER itself (z = nz-1) — first measurement there:
    // vel 3.5e-18, lid 4.1e-20, pressure 1.1e-16, seam 4.1e-20, i.e. exact
    // at machine precision, matching the whole-face Zou-He node exactness
    // class (T15 measured 6.9e-18). The adjacent interior layer z = nz-2
    // shows O(1e-3..1e-2) deviations that are developed FLOW, not BC error;
    // it stays printed as an informational probe only.
    let rim_z = dims[2] - 1;
    let z = dims[2] - 2;
    let interior = layer_error(z);
    let (vel_max, lid_max, pressure_max, seam_lid_max) = layer_error(rim_z);
    println!(
        "ACC SRC D4: bc_node_z={rim_z} vel_max={vel_max:.3e} lid_max={lid_max:.3e} pressure_max={pressure_max:.3e} seam_lid_max={seam_lid_max:.3e} interior_z={z} interior_flow_errors={interior:?} band=1.000e-13"
    );
    assert!(
        vel_max <= 1.0e-13,
        "ACC SRC D4 velocity patch max_abs={vel_max:.3e}, band=1e-13, denominator=componentwise absolute prescribed velocity"
    );
    assert!(
        lid_max <= 1.0e-13,
        "ACC SRC D4 zero-velocity lid max_abs={lid_max:.3e}, band=1e-13, denominator=absolute velocity"
    );
    assert!(
        seam_lid_max <= 1.0e-13,
        "ACC SRC D4 patch-lid seam max_abs={seam_lid_max:.3e}, band=1e-13, denominator=absolute velocity"
    );
    assert!(
        pressure_max <= 1.0e-13,
        "ACC SRC D4 pressure patch max_abs_rho_error={pressure_max:.3e}, band=1e-13, denominator=absolute rho"
    );
}

#[test]
fn d5_patch_mirror_equivariance() {
    let dims = [24, 24, 24];
    let u = [0.0, 0.0, -0.02];
    let left = z_top_patch([4, 10], [7, 13], FaceBC::Velocity { u });
    // In-face ZPos coordinates are [x,y]. Mirroring x in an N=24 face maps
    // inclusive [4,7] to [24-1-7, 24-1-4] = [16,19].
    let right = z_top_patch([16, 10], [19, 13], FaceBC::Velocity { u });
    let (spec_left, walls) = top_patch_spec(dims, vec![left]);
    let (spec_right, _) = top_patch_spec(dims, vec![right]);
    let mut a = build(&spec_left, &walls);
    let mut b = build(&spec_right, &walls);
    a.run(100);
    b.run(100);
    let d = max_mirror_delta(&a, &b, dims);
    let rel = linf_rel(&[d[0], d[1], d[2], d[3]], &[0.0, 0.0, 0.0, 0.0], 1.0);
    println!(
        "ACC SRC D5: max_delta rho={:.3e} ux_mirror={:.3e} uy={:.3e} uz={:.3e} linf_rel={rel:.3e} band=1.000e-12",
        d[0], d[1], d[2], d[3]
    );
    assert!(
        d.iter().all(|v| *v <= 1.0e-12),
        "ACC SRC D5 mirror max deltas={d:?}, band=1e-12, linf_rel_denominator=1 -> {rel:.3e}"
    );
}
