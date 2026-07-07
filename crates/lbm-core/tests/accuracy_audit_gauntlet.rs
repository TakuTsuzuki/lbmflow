//! ACC-AUDIT lane 5.6: degenerate-geometry gauntlet.
//!
//! These tests are intentionally adversarial geometry sentinels rather than
//! smooth analytic benchmarks. Each case names its configuration, invariant,
//! and current expected disposition in comments; every executed case runs at
//! least 100 steps, checks NaN-freeness, prints the measured quantities, and
//! asserts a concrete invariant.

use lbm_core::compat::prelude::{
    Collision as CCollision, EdgeBC as CEdgeBC, Edges as CEdges, SimConfig as CSimConfig,
    Simulation as CSimulation,
};
use lbm_core::prelude::*;

const CTRT: CCollision = CCollision::Trt {
    magic: CCollision::MAGIC_STD,
};
const NATIVE_TRT: CollisionKind = CollisionKind::Trt {
    magic: CollisionKind::MAGIC_STD,
};

fn cperiodic(nx: usize, ny: usize, force: [f64; 2]) -> CSimulation<f64> {
    CSimConfig {
        nx,
        ny,
        nu: 0.05,
        collision: CTRT,
        force,
        ..Default::default()
    }
    .build()
    .unwrap()
}

fn assert_compat_finite(sim: &CSimulation<f64>, label: &str) -> f64 {
    let mut max_speed = 0.0f64;
    for y in 0..sim.ny() {
        for x in 0..sim.nx() {
            let rho = sim.rho(x, y);
            let ux = sim.ux(x, y);
            let uy = sim.uy(x, y);
            assert!(
                rho.is_finite() && ux.is_finite() && uy.is_finite(),
                "{label}: non-finite at ({x},{y}): rho={rho:e}, ux={ux:e}, uy={uy:e}"
            );
            max_speed = max_speed.max((ux * ux + uy * uy).sqrt());
        }
    }
    max_speed
}

fn assert_native_finite<L: Lattice, H: HaloExchange<f64>>(
    s: &Solver<L, f64, CpuScalar, H>,
    label: &str,
) -> f64 {
    let rho = s.gather_rho();
    let ux = s.gather_ux();
    let uy = s.gather_uy();
    let uz = if L::D == 3 {
        s.gather_uz()
    } else {
        vec![0.0; rho.len()]
    };
    let mut max_speed = 0.0f64;
    for i in 0..rho.len() {
        assert!(
            rho[i].is_finite() && ux[i].is_finite() && uy[i].is_finite() && uz[i].is_finite(),
            "{label}: non-finite at compact index {i}: rho={:.6e}, ux={:.6e}, uy={:.6e}, uz={:.6e}",
            rho[i],
            ux[i],
            uy[i],
            uz[i]
        );
        max_speed = max_speed.max((ux[i] * ux[i] + uy[i] * uy[i] + uz[i] * uz[i]).sqrt());
    }
    max_speed
}

fn compat_mass_drift_rel(sim: &CSimulation<f64>, m0: f64) -> f64 {
    (sim.total_mass_f64() - m0).abs() / m0.abs()
}

fn idx3(dims: [usize; 3], x: usize, y: usize, z: usize) -> usize {
    (z * dims[1] + y) * dims[0] + x
}

fn solid_with_box(dims: [usize; 3], lo: [usize; 3], hi_exclusive: [usize; 3]) -> Vec<bool> {
    let mut solid = vec![false; dims[0] * dims[1] * dims[2]];
    for z in lo[2]..hi_exclusive[2] {
        for y in lo[1]..hi_exclusive[1] {
            for x in lo[0]..hi_exclusive[0] {
                solid[idx3(dims, x, y, z)] = true;
            }
        }
    }
    solid
}

fn max_d2q9_delta<HA: HaloExchange<f64>, HB: HaloExchange<f64>>(
    a: &Solver<D2Q9, f64, CpuScalar, HA>,
    b: &Solver<D2Q9, f64, CpuScalar, HB>,
) -> f64 {
    let mut max = 0.0f64;
    for (va, vb) in [
        (a.gather_rho(), b.gather_rho()),
        (a.gather_ux(), b.gather_ux()),
        (a.gather_uy(), b.gather_uy()),
    ] {
        max = max.max(
            va.iter()
                .zip(&vb)
                .map(|(x, y)| (x - y).abs())
                .fold(0.0, f64::max),
        );
    }
    for q in 0..D2Q9::Q {
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

fn kinetic_centroid_distance_to_origin(sim: &CSimulation<f64>) -> f64 {
    let mut wsum = 0.0;
    let mut xsum = 0.0;
    let mut ysum = 0.0;
    for y in 0..sim.ny() {
        for x in 0..sim.nx() {
            if sim.is_solid(x, y) {
                continue;
            }
            let ux = sim.ux(x, y);
            let uy = sim.uy(x, y);
            let w = ux * ux + uy * uy;
            wsum += w;
            xsum += w * x as f64;
            ysum += w * y as f64;
        }
    }
    if wsum == 0.0 {
        return f64::INFINITY;
    }
    let cx = xsum / wsum;
    let cy = ysum / wsum;
    (cx * cx + cy * cy).sqrt()
}

fn compat_max_vorticity(sim: &CSimulation<f64>) -> f64 {
    let mut max_vort = 0.0f64;
    for y in 1..sim.ny() - 1 {
        for x in 1..sim.nx() - 1 {
            if sim.is_solid(x, y) {
                continue;
            }
            let dvdx = 0.5 * (sim.uy(x + 1, y) - sim.uy(x - 1, y));
            let dudy = 0.5 * (sim.ux(x, y + 1) - sim.ux(x, y - 1));
            max_vort = max_vort.max((dvdx - dudy).abs());
        }
    }
    max_vort
}

// CASE G1: {config: 32x32, two 4x4 blocks at x=10..14 and 18..22,
// y=14..18, 4-cell fluid slot, vertical body force; invariant:
// left-right mirror equivariance; expected: PASS}.
#[test]
fn g1_one_cell_gap_left_right_mirror_equivariant() {
    let mut sim = cperiodic(32, 32, [0.0, 1.0e-6]);
    for y in 14..18 {
        for x in 10..14 {
            sim.set_solid(x, y);
        }
        for x in 18..22 {
            sim.set_solid(x, y);
        }
    }
    sim.run(500);
    let max_speed = assert_compat_finite(&sim, "G1");

    let mut max_mirror = 0.0f64;
    for y in 0..sim.ny() {
        for x in 0..sim.nx() {
            let mx = sim.nx() - 1 - x;
            assert_eq!(
                sim.is_solid(x, y),
                sim.is_solid(mx, y),
                "G1 solid mirror at ({x},{y})"
            );
            max_mirror = max_mirror.max((sim.rho(x, y) - sim.rho(mx, y)).abs());
            max_mirror = max_mirror.max((sim.ux(x, y) + sim.ux(mx, y)).abs());
            max_mirror = max_mirror.max((sim.uy(x, y) - sim.uy(mx, y)).abs());
        }
    }
    println!(
        "ACC GAUNTLET G1: one-cell-gap mirror max={max_mirror:.6e}, band=1.0e-12, max_speed={max_speed:.6e}, expected=PASS"
    );
    assert!(
        max_mirror <= 1.0e-12,
        "G1 mirror equivariance max={max_mirror:.6e}, band=1.0e-12"
    );
}

// CASE G2: {config: diagonal staircase solid y>x, x body force; invariant:
// solid-cell velocity is zero and kinetic-energy centroid does not drift
// monotonically into the acute corner; expected: PASS}.
#[test]
fn g2_diagonal_staircase_solid_enforcement_no_corner_drift() {
    let mut sim = cperiodic(32, 32, [0.0, 1.0e-6]);
    sim.set_solid_region(|x, y| y > x);

    let mut centroid_dist = Vec::new();
    for _ in 0..5 {
        sim.run(100);
        centroid_dist.push(kinetic_centroid_distance_to_origin(&sim));
    }
    let max_speed = assert_compat_finite(&sim, "G2");

    let mut max_solid_ux = 0.0f64;
    for y in 0..sim.ny() {
        for x in 0..sim.nx() {
            if y > x {
                assert!(sim.is_solid(x, y), "G2 y>x cell ({x},{y}) is not solid");
                max_solid_ux = max_solid_ux.max(sim.ux(x, y).abs());
            }
        }
    }
    let monotone_toward_corner = centroid_dist.windows(2).all(|w| w[1] < w[0] - 1.0e-12);
    println!(
        "ACC GAUNTLET G2: diagonal-staircase max_solid_ux={max_solid_ux:.6e}, band=1.0e-14, centroid_dist={centroid_dist:?}, monotone_toward_corner={monotone_toward_corner}, max_speed={max_speed:.6e}, expected=PASS"
    );
    assert!(
        max_solid_ux <= 1.0e-14,
        "G2 solid enforcement max_solid_ux={max_solid_ux:.6e}, band=1.0e-14"
    );
    assert!(
        !monotone_toward_corner,
        "G2 kinetic centroid drifts monotonically toward corner: {centroid_dist:?}"
    );
}

// CASE G3: {config: 4x4 obstacle at x=0..4, y=10..14 touching the left wall
// rim; invariant: build succeeds, finite closed-box evolution, mass drift;
// expected: PASS}.
#[test]
fn g3_obstacle_touching_rim_builds_and_conserves_mass() {
    let mut sim = CSimConfig {
        nx: 32,
        ny: 32,
        nu: 0.05,
        collision: CTRT,
        edges: CEdges {
            left: CEdgeBC::BounceBack,
            right: CEdgeBC::BounceBack,
            bottom: CEdgeBC::BounceBack,
            top: CEdgeBC::BounceBack,
        },
        force: [1.0e-6, 0.0],
    }
    .build()
    .unwrap();
    for y in 10..14 {
        for x in 0..4 {
            sim.set_solid(x, y);
        }
    }
    let m0 = sim.total_mass_f64();
    sim.run(200);
    let max_speed = assert_compat_finite(&sim, "G3");
    let drift = compat_mass_drift_rel(&sim, m0);
    println!(
        "ACC GAUNTLET G3: rim-touching obstacle mass_drift_rel={drift:.6e}, band=1.0e-11, max_speed={max_speed:.6e}, expected=PASS"
    );
    assert!(
        drift <= 1.0e-11,
        "G3 mass drift rel={drift:.6e}, band=1.0e-11"
    );
}

// CASE G4: {config: native D2Q9 32x32, 4x4 obstacle straddling the 2x2
// partition seam; invariant: split is bit-identical to unsplit; expected:
// PASS}.
#[test]
fn g4_obstacle_straddling_partition_seam_bit_identical() {
    let dims = [32, 32, 1];
    let solid = solid_with_box(dims, [14, 14, 0], [18, 18, 1]);
    let wall_u = vec![[0.0; 3]; solid.len()];
    let spec = GlobalSpec {
        dims,
        nu: 0.05,
        collision: NATIVE_TRT,
        force: [1.0e-6, 0.0, 0.0],
        ..Default::default()
    };
    let mut base: Solver<D2Q9, f64, CpuScalar, LocalPeriodic> = Solver::new(
        &spec,
        &solid,
        &wall_u,
        [1, 1, 1],
        CpuScalar::default(),
        LocalPeriodic,
    );
    let mut split: Solver<D2Q9, f64, CpuScalar, InProcess> = Solver::new(
        &spec,
        &solid,
        &wall_u,
        [2, 2, 1],
        CpuScalar::default(),
        InProcess,
    );
    for _ in 0..150 {
        base.step();
        split.step();
    }
    let max_speed = assert_native_finite(&split, "G4 split");
    let max_delta = max_d2q9_delta(&base, &split);
    println!(
        "ACC GAUNTLET G4: seam-obstacle split max_field_delta={max_delta:.6e}, band=0.0, max_speed={max_speed:.6e}, expected=PASS"
    );
    assert_eq!(
        max_delta, 0.0,
        "G4 T13-style split field delta={max_delta:.6e}, band=0.0"
    );
}

fn patch_open_box_solid_with_block(dims: [usize; 3]) -> (Vec<bool>, Vec<[f64; 3]>) {
    let mut solid = vec![false; dims[0] * dims[1] * dims[2]];
    for z in 0..dims[2] {
        for y in 0..dims[1] {
            for x in 0..dims[0] {
                let rim = x == 0
                    || y == 0
                    || z == 0
                    || x + 1 == dims[0]
                    || y + 1 == dims[1]
                    || z + 1 == dims[2];
                let top_patch =
                    z + 1 == dims[2] && (10..=14).contains(&x) && (10..=14).contains(&y);
                let block = (11..14).contains(&x) && (11..14).contains(&y) && (20..23).contains(&z);
                solid[idx3(dims, x, y, z)] = (rim && !top_patch) || block;
            }
        }
    }
    let wall_u = vec![[0.0; 3]; solid.len()];
    (solid, wall_u)
}

// CASE G5: {config: native D3Q19 24^3, top velocity patch x/y=10..14 with a
// 3x3x3 obstacle immediately below its center; invariant: finite run and
// prescribed patch-node velocity; expected: PASS unless the open-patch/solid
// contact contract is reclassified as SPEC-GAP}.
#[test]
fn g5_obstacle_touching_face_patch_preserves_patch_velocity() {
    let dims = [24, 24, 24];
    let patch_u = [0.0, 0.0, -0.02];
    let spec = GlobalSpec {
        dims,
        nu: 0.05,
        collision: NATIVE_TRT,
        periodic: [false, false, false],
        face_patches: vec![FacePatch {
            face: Face::ZPos.index(),
            lo: [10, 10],
            hi: [14, 14],
            bc: FaceBC::Velocity { u: patch_u },
        }],
        ..Default::default()
    };
    let (solid, wall_u) = patch_open_box_solid_with_block(dims);
    spec.validate_lattice::<D3Q19>(&solid)
        .expect("G5 expected PASS: face patch touching an interior obstacle should validate");
    let mut s: Solver<D3Q19, f64, CpuScalar, LocalPeriodic> = Solver::new(
        &spec,
        &solid,
        &wall_u,
        [1, 1, 1],
        CpuScalar::default(),
        LocalPeriodic,
    );
    s.run(120);
    let max_speed = assert_native_finite(&s, "G5");
    let mut max_patch_err = 0.0f64;
    for y in 10..=14 {
        for x in 10..=14 {
            let u = s.u(x, y, 23);
            for a in 0..3 {
                max_patch_err = max_patch_err.max((u[a] - patch_u[a]).abs());
            }
        }
    }
    println!(
        "ACC GAUNTLET G5: face-patch contact max_patch_velocity_err={max_patch_err:.6e}, band=1.0e-12, max_speed={max_speed:.6e}, expected=PASS"
    );
    assert!(
        max_patch_err <= 1.0e-12,
        "G5 patch prescribed velocity max_err={max_patch_err:.6e}, band=1.0e-12"
    );
}

// CASE G6: {config: single-cell-wide diagonal fin from (5,5) to (26,26);
// invariant: build succeeds, finite run, bounded velocity/vorticity; expected:
// PASS. If a future geometry validator rejects zero-thickness fins, replace
// this with an exact ConfigError/SpecError assertion}.
#[test]
fn g6_zero_thickness_diagonal_fin_bounded_vorticity() {
    let mut sim = cperiodic(32, 32, [1.0e-6, 0.0]);
    for i in 5..=26 {
        sim.set_solid(i, i);
    }

    sim.run(100);
    let vort100 = compat_max_vorticity(&sim);
    sim.run(200);
    let vort300 = compat_max_vorticity(&sim);
    sim.run(200);
    let vort500 = compat_max_vorticity(&sim);
    let max_speed = assert_compat_finite(&sim, "G6");
    println!(
        "ACC GAUNTLET G6: zero-thickness fin vort100={vort100:.6e}, vort300={vort300:.6e}, vort500={vort500:.6e}, growth_band=10x+1e-12, max_speed={max_speed:.6e}, speed_band=5.0e-2, expected=PASS"
    );
    assert!(
        vort500 <= 10.0 * vort100.max(1.0e-12),
        "G6 no-diverging-vortex vorticity: vort500={vort500:.6e}, vort100={vort100:.6e}, band=10x"
    );
    assert!(
        max_speed <= 5.0e-2,
        "G6 max_speed={max_speed:.6e}, band=5.0e-2"
    );
}

// CASE G7: {config: moving top wall U=0.05, stationary bottom/side walls,
// 3x3 solid block at x/y=1..4 touching the bottom-left rim corner; invariant:
// finite Couette-like bulk shear and mass drift; expected: PASS}.
#[test]
fn g7_moving_wall_corner_with_adjacent_obstacle_forms_bulk_shear() {
    let mut sim = CSimConfig {
        nx: 96,
        ny: 32,
        nu: 0.1,
        collision: CTRT,
        edges: CEdges {
            left: CEdgeBC::BounceBack,
            right: CEdgeBC::BounceBack,
            bottom: CEdgeBC::BounceBack,
            top: CEdgeBC::MovingWall { u: [0.05, 0.0] },
        },
        ..Default::default()
    }
    .build()
    .unwrap();
    for y in 1..4 {
        for x in 1..4 {
            sim.set_solid(x, y);
        }
    }
    let m0 = sim.total_mass_f64();
    sim.run(4_000);
    let max_speed = assert_compat_finite(&sim, "G7");
    let drift = compat_mass_drift_rel(&sim, m0);

    let h = (sim.ny() - 2) as f64;
    let ref_slope = 0.05 / h;
    let mut n = 0.0;
    let mut sum_y = 0.0;
    let mut sum_u = 0.0;
    let mut sum_yy = 0.0;
    let mut sum_yu = 0.0;
    let mut max_uy = 0.0f64;
    for y in 1..sim.ny() - 1 {
        let yw = y as f64 - 0.5;
        for x in 40..56 {
            if sim.is_solid(x, y) {
                continue;
            }
            let ux = sim.ux(x, y);
            n += 1.0;
            sum_y += yw;
            sum_u += ux;
            sum_yy += yw * yw;
            sum_yu += yw * ux;
            max_uy = max_uy.max(sim.uy(x, y).abs());
        }
    }
    let slope = (n * sum_yu - sum_y * sum_u) / (n * sum_yy - sum_y * sum_y);
    let slope_ratio = slope / ref_slope;
    println!(
        "ACC GAUNTLET G7: moving-wall corner mass_drift_rel={drift:.6e}, mass_band=1.0e-11, bulk_slope={slope:.6e}, ref_slope={ref_slope:.6e}, slope_ratio={slope_ratio:.3}, ratio_band=[0.30,1.40], max_uy={max_uy:.6e}, uy_band=1.0e-2, max_speed={max_speed:.6e}, expected=PASS"
    );
    assert!(
        drift <= 1.0e-11,
        "G7 mass drift rel={drift:.6e}, band=1.0e-11"
    );
    assert!(
        (0.30..=1.40).contains(&slope_ratio),
        "G7 Couette-like bulk slope ratio={slope_ratio:.6e}, band=[0.30,1.40], slope={slope:.6e}, ref={ref_slope:.6e}"
    );
    assert!(
        max_uy <= 1.0e-2,
        "G7 bulk transverse velocity max_uy={max_uy:.6e}, band=1.0e-2"
    );
}

// CASE G8: {config: single-cell obstacle at box center with force probe;
// invariant: finite run and no accumulating force asymmetry, i.e. force norm
// stays below 100x its initial post-step norm; expected: PASS}.
#[test]
fn g8_single_cell_obstacle_force_probe_bounded() {
    let mut sim = cperiodic(32, 32, [1.0e-6, 0.0]);
    let center = (16usize, 16usize);
    sim.set_solid(center.0, center.1);
    sim.set_force_probe(move |x, y| (x, y) == center);

    sim.step();
    let f0 = sim.probed_force();
    let f0_norm = (f0[0] * f0[0] + f0[1] * f0[1]).sqrt().max(1.0e-18);
    let mut max_norm = f0_norm;
    for _ in 1..500 {
        sim.step();
        let f = sim.probed_force();
        max_norm = max_norm.max((f[0] * f[0] + f[1] * f[1]).sqrt());
    }
    let max_speed = assert_compat_finite(&sim, "G8");
    let ratio = max_norm / f0_norm;
    println!(
        "ACC GAUNTLET G8: single-cell obstacle force0_norm={f0_norm:.6e}, max_force_norm={max_norm:.6e}, ratio={ratio:.3}, band=100, max_speed={max_speed:.6e}, expected=PASS"
    );
    assert!(
        ratio <= 100.0,
        "G8 force asymmetry ratio={ratio:.6e}, band=100, f0_norm={f0_norm:.6e}, max_norm={max_norm:.6e}"
    );
}
