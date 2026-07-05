//! Bouzidi curved-wall regression tests.

use lbm_core::lattice::D2Q9;
use lbm_core::prelude::*;

fn channel_spec(nx: usize, ny: usize) -> (GlobalSpec<f64>, WallSpec<f64>) {
    let mut walls = WallSpec::default();
    walls.is_wall[Face::YNeg.index()] = true;
    walls.is_wall[Face::YPos.index()] = true;
    let mut faces = [FaceBC::Closed; 6];
    faces[Face::XNeg.index()] = FaceBC::Velocity {
        u: [0.04, 0.0, 0.0],
    };
    faces[Face::XPos.index()] = FaceBC::Pressure { rho: 1.0 };
    (
        GlobalSpec {
            dims: [nx, ny, 1],
            nu: 0.04,
            periodic: [false, false, false],
            faces,
            ..Default::default()
        },
        walls,
    )
}

fn build_scalar() -> Solver<D2Q9, f64, CpuScalar, LocalPeriodic> {
    let (spec, walls) = channel_spec(80, 42);
    let (solid, wall_u) = build_wall_rims(2, spec.dims, &walls);
    Solver::new(
        &spec,
        &solid,
        &wall_u,
        [1, 1, 1],
        CpuScalar::default(),
        LocalPeriodic,
    )
}

fn build_simd() -> Solver<D2Q9, f64, CpuSimd, LocalPeriodic> {
    let (spec, walls) = channel_spec(80, 42);
    let (solid, wall_u) = build_wall_rims(2, spec.dims, &walls);
    Solver::new(
        &spec,
        &solid,
        &wall_u,
        [1, 1, 1],
        CpuSimd::default(),
        LocalPeriodic,
    )
}

fn add_cylinder<B, H>(s: &mut Solver<D2Q9, f64, B, H>, bouzidi: bool)
where
    B: Backend<D2Q9, f64, Fields = SoaFields<f64>>,
    H: HaloExchange<f64>,
{
    let (cx, cy, r) = (31.5, 20.0, 6.5);
    let inside = |x: usize, y: usize| {
        let dx = x as f64 - cx;
        let dy = y as f64 - cy;
        dx * dx + dy * dy <= r * r
    };
    for y in 0..42 {
        for x in 0..80 {
            if inside(x, y) {
                s.set_solid(x, y, 0);
            }
        }
    }
    s.set_force_probe(move |x, y, _| inside(x, y));
    if bouzidi {
        s.set_bouzidi_circle(cx, cy, r);
    }
    s.init_with(|x, y, _| {
        if inside(x, y) || x == 0 || x == 79 || y == 0 || y == 41 {
            (1.0, [0.0, 0.0, 0.0])
        } else {
            (1.0, [0.04, 1.0e-5 * (y as f64 - cy), 0.0])
        }
    });
}

fn assert_bitwise_same(
    a: &Solver<D2Q9, f64, CpuScalar, LocalPeriodic>,
    b: &Solver<D2Q9, f64, CpuScalar, LocalPeriodic>,
) {
    assert_eq!(a.gather_rho(), b.gather_rho(), "rho differs");
    assert_eq!(a.gather_ux(), b.gather_ux(), "ux differs");
    assert_eq!(a.gather_uy(), b.gather_uy(), "uy differs");
    assert_eq!(a.probed_force(), b.probed_force(), "probe differs");
    for q in 0..D2Q9::Q {
        assert_eq!(a.gather_f(q), b.gather_f(q), "f[{q}] differs");
    }
}

#[test]
fn qd_half_records_are_bitwise_half_way_bounce_back() {
    let mut half = build_scalar();
    let mut bz = build_scalar();
    add_cylinder(&mut half, false);
    add_cylinder(&mut bz, false);
    bz.set_bouzidi_half_way_links();
    for _ in 0..20 {
        half.step();
        bz.step();
        assert_bitwise_same(&half, &bz);
    }
}

#[test]
fn analytic_circle_records_are_sorted_and_nonempty() {
    let mut s = build_scalar();
    add_cylinder(&mut s, true);
    let records = &s.fields(0).bouzidi.as_ref().unwrap().records;
    assert!(!records.is_empty());
    assert!(records
        .windows(2)
        .all(|w| (w[0].cell, w[0].q) <= (w[1].cell, w[1].q)));
    assert!(records.iter().all(|r| r.qd > 0.0 && r.qd < 1.0));
    assert!(records.iter().any(|r| r.qd != 0.5));
}

#[test]
fn cpu_scalar_and_simd_match_with_bouzidi_circle() {
    let mut scalar = build_scalar();
    let mut simd = build_simd();
    add_cylinder(&mut scalar, true);
    add_cylinder(&mut simd, true);
    for _ in 0..120 {
        scalar.step();
        simd.step();
    }
    for (name, a, b) in [
        ("rho", scalar.gather_rho(), simd.gather_rho()),
        ("ux", scalar.gather_ux(), simd.gather_ux()),
        ("uy", scalar.gather_uy(), simd.gather_uy()),
    ] {
        let d = a
            .iter()
            .zip(&b)
            .map(|(x, y)| (x - y).abs())
            .fold(0.0, f64::max);
        assert!(d <= 1e-11, "{name} max delta {d:e}");
    }
    let fa = scalar.probed_force();
    let fb = simd.probed_force();
    for c in 0..2 {
        let d = (fa[c] - fb[c]).abs();
        assert!(d <= 1e-10, "probe[{c}] delta {d:e}");
    }
}
