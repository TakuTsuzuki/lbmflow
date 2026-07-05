//! Adversarial T13 tests for the V2 decomposition layer.
//!
//! These cases intentionally put solids, open boundaries, moving walls,
//! per-cell forces, and D3Q19 table assumptions on decomposition seams.

use std::panic::{self, AssertUnwindSafe};

use lbm_core::lattice::{D2Q9, D3Q19};
use lbm_core::prelude::*;

type Sol<L, H> = Solver<L, f64, CpuScalar, H>;

fn build<L: Lattice, H: HaloExchange<f64>>(
    spec: &GlobalSpec<f64>,
    walls: &WallSpec<f64>,
    decomp: [usize; 3],
    exchange: H,
    two_pass: bool,
) -> Sol<L, H> {
    let (solid, wall_u) = build_wall_rims(L::D, spec.dims, walls);
    let mut s = Solver::new(spec, &solid, &wall_u, decomp, CpuScalar::default(), exchange);
    s.set_two_pass(two_pass);
    s
}

fn max_abs(a: &[f64], b: &[f64]) -> f64 {
    assert_eq!(a.len(), b.len());
    a.iter()
        .zip(b)
        .map(|(x, y)| (x - y).abs())
        .fold(0.0, f64::max)
}

fn assert_close<L: Lattice, HA: HaloExchange<f64>, HB: HaloExchange<f64>>(
    a: &Sol<L, HA>,
    b: &Sol<L, HB>,
    field_tol: f64,
    probe_tol: f64,
    what: &str,
) {
    for (name, va, vb) in [
        ("rho", a.gather_rho(), b.gather_rho()),
        ("ux", a.gather_ux(), b.gather_ux()),
        ("uy", a.gather_uy(), b.gather_uy()),
        ("uz", a.gather_uz(), b.gather_uz()),
    ] {
        let d = max_abs(&va, &vb);
        assert!(
            d <= field_tol,
            "{what}: {name} max|delta| = {d:e} > {field_tol:e}"
        );
    }
    for q in 0..L::Q {
        let d = max_abs(&a.gather_f(q), &b.gather_f(q));
        assert!(
            d <= field_tol,
            "{what}: f[{q}] max|delta| = {d:e} > {field_tol:e}"
        );
    }
    let (fa, fb) = (a.probed_force(), b.probed_force());
    for c in 0..L::D {
        let d = (fa[c] - fb[c]).abs();
        assert!(
            d <= probe_tol,
            "{what}: probed_force[{c}] delta = {d:e} > {probe_tol:e} (mono {:e}, split {:e})",
            fa[c],
            fb[c]
        );
    }
}

fn channel_spec(nx: usize, ny: usize, outlet: FaceBC<f64>) -> (GlobalSpec<f64>, WallSpec<f64>) {
    let mut walls = WallSpec::default();
    walls.is_wall[Face::YNeg.index()] = true;
    walls.is_wall[Face::YPos.index()] = true;
    let mut faces = [FaceBC::Closed; 6];
    faces[Face::XNeg.index()] = FaceBC::Velocity {
        u: [0.045, 0.0, 0.0],
    };
    faces[Face::XPos.index()] = outlet;
    (
        GlobalSpec {
            dims: [nx, ny, 1],
            nu: 0.035,
            periodic: [false, false, false],
            faces,
            ..Default::default()
        },
        walls,
    )
}

fn parabolic_profile(ny: usize, umax: f64) -> Vec<[f64; 3]> {
    (0..ny)
        .map(|y| {
            let yy = y as f64 / (ny - 1) as f64;
            [umax * 4.0 * yy * (1.0 - yy), 0.0, 0.0]
        })
        .collect()
}

fn set_obstacle_pair(
    base: &mut Sol<D2Q9, LocalPeriodic>,
    split: &mut Sol<D2Q9, InProcess>,
    nx: usize,
    ny: usize,
    inside: impl Fn(usize, usize) -> bool + Copy,
) {
    for y in 0..ny {
        for x in 0..nx {
            if inside(x, y) {
                base.set_solid(x, y, 0);
                split.set_solid(x, y, 0);
            }
        }
    }
    base.set_force_probe(move |x, y, _| inside(x, y));
    split.set_force_probe(move |x, y, _| inside(x, y));
}

#[test]
fn t13_cylinder_centered_exactly_on_split_line_probe_matches() {
    let (nx, ny) = (64usize, 40usize);
    let (spec, walls) = channel_spec(nx, ny, FaceBC::Pressure { rho: 1.0 });
    let profile = parabolic_profile(ny, 0.07);
    let cx = nx as f64 / 2.0;
    let cy = ny as f64 / 2.0;
    let r2 = 5.5f64 * 5.5;
    let inside = move |x: usize, y: usize| {
        let dx = x as f64 - cx;
        let dy = y as f64 - cy;
        dx * dx + dy * dy <= r2
    };

    for decomp in [[2, 1, 1], [2, 2, 1]] {
        let mut base = build::<D2Q9, _>(&spec, &walls, [1, 1, 1], LocalPeriodic, false);
        let mut split = build::<D2Q9, _>(&spec, &walls, decomp, InProcess, true);
        base.set_inlet_profile(Face::XNeg, &profile);
        split.set_inlet_profile(Face::XNeg, &profile);
        set_obstacle_pair(&mut base, &mut split, nx, ny, inside);
        for t in 1..=180 {
            base.step();
            split.step();
            if t <= 3 || t % 45 == 0 || t == 180 {
                assert_close(
                    &base,
                    &split,
                    1e-12,
                    1e-11,
                    &format!("exact split-line cylinder {decomp:?} t={t}"),
                );
            }
        }
    }
}

#[test]
fn t13_l_shaped_obstacle_spanning_three_subdomains_matches() {
    let (nx, ny) = (64usize, 48usize);
    let (spec, walls) = channel_spec(nx, ny, FaceBC::Outflow);
    let inside = |x: usize, y: usize| {
        let vertical = (28..32).contains(&x) && (8..36).contains(&y);
        let horizontal = (28..54).contains(&x) && (24..28).contains(&y);
        vertical || horizontal
    };
    let mut base = build::<D2Q9, _>(&spec, &walls, [1, 1, 1], LocalPeriodic, false);
    let mut split = build::<D2Q9, _>(&spec, &walls, [2, 2, 1], InProcess, true);
    set_obstacle_pair(&mut base, &mut split, nx, ny, inside);
    for t in 1..=160 {
        base.step();
        split.step();
        if t <= 3 || t % 40 == 0 || t == 160 {
            assert_close(
                &base,
                &split,
                1e-12,
                1e-11,
                &format!("L obstacle 2x2 t={t}"),
            );
        }
    }
}

#[test]
fn t13_open_and_moving_boundary_faces_split_across_ranks_match() {
    let (spec_p, walls_p) = channel_spec(72, 34, FaceBC::Pressure { rho: 1.0 });
    let profile = parabolic_profile(34, 0.065);
    for decomp in [[1, 2, 1], [2, 2, 1]] {
        let mut base = build::<D2Q9, _>(&spec_p, &walls_p, [1, 1, 1], LocalPeriodic, false);
        let mut split = build::<D2Q9, _>(&spec_p, &walls_p, decomp, InProcess, true);
        base.set_inlet_profile(Face::XNeg, &profile);
        split.set_inlet_profile(Face::XNeg, &profile);
        for t in 1..=160 {
            base.step();
            split.step();
            if t <= 3 || t % 40 == 0 || t == 160 {
                assert_close(
                    &base,
                    &split,
                    1e-12,
                    1e-11,
                    &format!("velocity-inlet/pressure-outlet {decomp:?} t={t}"),
                );
            }
        }
    }

    let mut walls = WallSpec::default();
    for f in [Face::XNeg, Face::XPos, Face::YNeg, Face::YPos] {
        walls.is_wall[f.index()] = true;
    }
    walls.u[Face::YPos.index()] = [0.085, 0.0, 0.0];
    let spec = GlobalSpec {
        dims: [50, 42, 1],
        nu: 0.03,
        periodic: [false, false, false],
        ..Default::default()
    };
    for decomp in [[2, 1, 1], [2, 2, 1]] {
        let mut base = build::<D2Q9, _>(&spec, &walls, [1, 1, 1], LocalPeriodic, false);
        let mut split = build::<D2Q9, _>(&spec, &walls, decomp, InProcess, true);
        base.set_force_probe(|_, y, _| y == 41);
        split.set_force_probe(|_, y, _| y == 41);
        for t in 1..=180 {
            base.step();
            split.step();
            if t <= 3 || t % 45 == 0 || t == 180 {
                assert_close(
                    &base,
                    &split,
                    1e-12,
                    1e-11,
                    &format!("moving lid spans ranks {decomp:?} t={t}"),
                );
            }
        }
    }
}

#[test]
fn t13_uneven_311_split_and_min_width_guard() {
    let spec = GlobalSpec {
        dims: [53, 31, 1],
        nu: 0.04,
        periodic: [true, true, false],
        ..Default::default()
    };
    let walls = WallSpec::default();
    let init = |x: usize, y: usize, _z: usize| {
        let kx = 2.0 * std::f64::consts::PI / 53.0;
        let ky = 2.0 * std::f64::consts::PI / 31.0;
        (
            1.0 + 0.003 * (kx * x as f64).cos(),
            [
                0.025 * (ky * y as f64).sin(),
                -0.02 * (kx * x as f64).sin(),
                0.0,
            ],
        )
    };
    let mut base = build::<D2Q9, _>(&spec, &walls, [1, 1, 1], LocalPeriodic, false);
    let mut split = build::<D2Q9, _>(&spec, &walls, [3, 1, 1], InProcess, true);
    base.init_with(init);
    split.init_with(init);
    for t in 1..=140 {
        base.step();
        split.step();
        if t <= 3 || t % 35 == 0 || t == 140 {
            assert_close(&base, &split, 1e-12, 1e-11, &format!("uneven [3,1,1] t={t}"));
        }
    }

    let bad = panic::catch_unwind(AssertUnwindSafe(|| {
        let spec = GlobalSpec {
            dims: [5, 8, 1],
            periodic: [true, true, false],
            ..Default::default()
        };
        let _: Sol<D2Q9, InProcess> =
            build::<D2Q9, _>(&spec, &WallSpec::default(), [3, 1, 1], InProcess, false);
    }));
    assert!(bad.is_err(), "expected min-width guard to reject 5 cells split over 3 parts");
}

fn set_force_field<L: Lattice, H: HaloExchange<f64>>(
    s: &mut Sol<L, H>,
    force: impl Fn(usize, usize, usize) -> [f64; 3],
) {
    for pi in 0..s.part_count() {
        let sub = s.sub(pi).clone();
        let g = sub.geom;
        let fields = s.fields_mut(pi);
        let ff = fields
            .force_field
            .get_or_insert_with(|| vec![[0.0; 3]; g.n_core()]);
        for z in 0..g.core[2] {
            for y in 0..g.core[1] {
                for x in 0..g.core[0] {
                    ff[g.cidx(x, y, z)] =
                        force(sub.origin[0] + x, sub.origin[1] + y, sub.origin[2] + z);
                }
            }
        }
    }
}

#[test]
fn t13_hand_rolled_per_cell_force_droplet_on_four_rank_corner_matches() {
    let (nx, ny) = (48usize, 48usize);
    let spec = GlobalSpec {
        dims: [nx, ny, 1],
        nu: 0.045,
        periodic: [true, true, false],
        force: [1e-6, -2e-6, 0.0],
        ..Default::default()
    };
    let walls = WallSpec::default();
    let init = move |x: usize, y: usize, _z: usize| {
        let dx = x as f64 + 0.5 - nx as f64 / 2.0;
        let dy = y as f64 + 0.5 - ny as f64 / 2.0;
        let r2 = dx * dx + dy * dy;
        (1.0 + 0.08 * (-r2 / 60.0).exp(), [0.0, 0.0, 0.0])
    };
    let force = move |x: usize, y: usize, _z: usize| {
        let dx = x as f64 + 0.5 - nx as f64 / 2.0;
        let dy = y as f64 + 0.5 - ny as f64 / 2.0;
        let amp = 2.5e-5 * (-(dx * dx + dy * dy) / 80.0).exp();
        [-amp * dx / nx as f64, -amp * dy / ny as f64, 0.0]
    };
    let mut base = build::<D2Q9, _>(&spec, &walls, [1, 1, 1], LocalPeriodic, false);
    let mut split = build::<D2Q9, _>(&spec, &walls, [2, 2, 1], InProcess, true);
    base.init_with(init);
    split.init_with(init);
    for t in 0..180 {
        set_force_field(&mut base, force);
        set_force_field(&mut split, force);
        base.step();
        split.step();
        if t < 3 || t % 45 == 44 || t == 179 {
            assert_close(
                &base,
                &split,
                1e-12,
                1e-11,
                &format!("per-cell force droplet corner t={}", t + 1),
            );
        }
    }
}

#[test]
fn d3q19_lattice_properties_from_all_angles() {
    for q in 0..D3Q19::Q {
        let opp = D3Q19::OPP[q];
        assert_eq!(D3Q19::OPP[opp], q, "OPP involution failed at q={q}");
        assert_eq!(D3Q19::W[opp], D3Q19::W[q], "W/OPP mismatch at q={q}");
        for a in 0..3 {
            assert_eq!(
                D3Q19::C[opp][a],
                -D3Q19::C[q][a],
                "C/OPP mismatch q={q} axis={a}"
            );
        }
    }

    let delta = |a: usize, b: usize| if a == b { 1.0 } else { 0.0 };
    for a in 0..3 {
        for b in 0..3 {
            for c in 0..3 {
                for d in 0..3 {
                    let got: f64 = (0..D3Q19::Q)
                        .map(|q| {
                            D3Q19::W[q]
                                * (D3Q19::C[q][a] as i32
                                    * D3Q19::C[q][b] as i32
                                    * D3Q19::C[q][c] as i32
                                    * D3Q19::C[q][d] as i32) as f64
                        })
                        .sum();
                    let want = D3Q19::CS2
                        * D3Q19::CS2
                        * (delta(a, b) * delta(c, d)
                            + delta(a, c) * delta(b, d)
                            + delta(a, d) * delta(b, c));
                    assert!(
                        (got - want).abs() <= 1e-15,
                        "D3Q19 4th moment [{a}{b}{c}{d}] = {got:e}, want {want:e}"
                    );
                }
            }
        }
    }

    for face in Face::ALL {
        let n = face.n_in();
        let unknowns = D3Q19::unknowns(face);
        assert_eq!(unknowns.len(), 5, "{face:?} unknown count");
        for &q in unknowns {
            let dot: i32 = (0..3).map(|a| D3Q19::C[q][a] as i32 * n[a] as i32).sum();
            assert!(dot > 0, "{face:?} unknown q={q} has dot={dot}");
        }
        let mut closure = 0.0;
        for q in 0..D3Q19::Q {
            let dot: i32 = (0..3).map(|a| D3Q19::C[q][a] as i32 * n[a] as i32).sum();
            if dot == 0 {
                closure += D3Q19::W[q];
            } else if dot < 0 {
                closure += 2.0 * D3Q19::W[q];
            }
        }
        // Analytically the constant is exactly 1 for every face; the engine
        // hardcodes the analytic value (kernels.rs zou_he: `+ T::one()`).
        // A runtime f64 sum of the weights depends on summation order and can
        // land 1 ulp off (XNeg measured 1.0000000000000002), so assert
        // consistency to 4 ulp rather than exactness of this test's own sum.
        assert!(
            (closure - 1.0f64).abs() <= 4.0 * f64::EPSILON,
            "{face:?} closure constant {closure} deviates from the analytic 1"
        );
    }

    let rotations: [[usize; 3]; 6] = [
        [0, 1, 2],
        [1, 2, 0],
        [2, 0, 1],
        [1, 0, 2],
        [0, 2, 1],
        [2, 1, 0],
    ];
    for perm in rotations {
        for q in 0..D3Q19::Q {
            let c = D3Q19::C[q];
            let mapped = [c[perm[0]], c[perm[1]], c[perm[2]]];
            let r = D3Q19::dir_index(mapped);
            assert_eq!(
                D3Q19::W[r],
                D3Q19::W[q],
                "rotation {perm:?}: weight changed q={q} -> {r}"
            );
            assert_eq!(
                D3Q19::dir_index([
                    -mapped[0],
                    -mapped[1],
                    -mapped[2],
                ]),
                D3Q19::OPP[r],
                "rotation {perm:?}: opposite closure failed q={q}"
            );
        }
    }
}

#[test]
#[ignore = "20k-step adversarial stability run; execute manually when validating T13 depth"]
fn t13_long_run_cavity_2x2_matches_after_20k_steps() {
    let mut walls = WallSpec::default();
    for f in [Face::XNeg, Face::XPos, Face::YNeg, Face::YPos] {
        walls.is_wall[f.index()] = true;
    }
    walls.u[Face::YPos.index()] = [0.05, 0.0, 0.0];
    let spec = GlobalSpec {
        dims: [40, 40, 1],
        nu: 0.04,
        periodic: [false, false, false],
        ..Default::default()
    };
    let mut base = build::<D2Q9, _>(&spec, &walls, [1, 1, 1], LocalPeriodic, false);
    let mut split = build::<D2Q9, _>(&spec, &walls, [2, 2, 1], InProcess, true);
    for _ in 0..20_000 {
        base.step();
        split.step();
    }
    for v in split
        .gather_rho()
        .into_iter()
        .chain(split.gather_ux())
        .chain(split.gather_uy())
    {
        assert!(v.is_finite(), "long-run split field contains non-finite value {v}");
    }
    assert_close(&base, &split, 1e-10, 1e-10, "20k cavity 2x2");
}
