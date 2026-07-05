//! T13 seed: decomposition invariance (COMPETITIVE_SPEC.md §4).
//!
//! The same scenario run as one monolithic domain and as 2×1 / 1×2 / 2×2
//! in-process subdomain decompositions must produce identical physics:
//! per-cell arithmetic is decomposition-independent and the halo carries
//! exact copies, so the *fields* are required to match bit-for-bit (== 0.0,
//! stronger than the ≤1e-12 acceptance line). Global f64 diagnostics
//! (mass / momentum / probe force) are summed per part and then combined,
//! so they may differ by reassociation only: ≤1e-12 relative.
//!
//! Each split case also runs with two-pass streaming (interior first, then
//! the boundary shell — the communication-overlap seam) and must again match
//! bit-for-bit.

use lbm_core::lattice::D2Q9;
use lbm_core::prelude::*;

type S<H> = Solver<D2Q9, f64, CpuScalar, H>;

struct Case {
    spec: GlobalSpec<f64>,
    walls: WallSpec<f64>,
}

fn build<H: HaloExchange<f64>>(case: &Case, decomp: [usize; 3], ex: H, two_pass: bool) -> S<H> {
    let (solid, wall_u) = build_wall_rims(2, case.spec.dims, &case.walls);
    let mut s = Solver::new(
        &case.spec,
        &solid,
        &wall_u,
        decomp,
        CpuScalar::default(),
        ex,
    );
    s.set_two_pass(two_pass);
    s
}

fn assert_fields_equal<HA: HaloExchange<f64>, HB: HaloExchange<f64>>(
    a: &S<HA>,
    b: &S<HB>,
    what: &str,
) {
    let pairs = [
        ("rho", a.gather_rho(), b.gather_rho()),
        ("ux", a.gather_ux(), b.gather_ux()),
        ("uy", a.gather_uy(), b.gather_uy()),
    ];
    for (name, va, vb) in pairs {
        let d = va
            .iter()
            .zip(&vb)
            .map(|(x, y)| (x - y).abs())
            .fold(0.0f64, f64::max);
        assert_eq!(d, 0.0, "{what}: field {name} differs by {d:e}");
    }
    // Populations too (the strongest statement).
    for q in 0..9 {
        let (fa, fb) = (a.gather_f(q), b.gather_f(q));
        let d = fa
            .iter()
            .zip(&fb)
            .map(|(x, y)| (x - y).abs())
            .fold(0.0f64, f64::max);
        assert_eq!(d, 0.0, "{what}: f[{q}] differs by {d:e}");
    }
}

/// Global diagnostics are per-part f64 partial sums combined in part order;
/// versus the monolithic single sum they may differ by reassociation only:
/// `|Δ| ≤ atol + rtol·|ref|` with both at 1e-12 (observed ~1e-15).
fn assert_diagnostics_close<HA: HaloExchange<f64>, HB: HaloExchange<f64>>(
    a: &S<HA>,
    b: &S<HB>,
    tol: f64,
    what: &str,
) {
    let close = |x: f64, y: f64, name: &str| {
        let d = (x - y).abs();
        assert!(
            d <= tol + tol * x.abs(),
            "{what}: {name} Δ = {d:e} (ref {x:e})"
        );
    };
    close(a.total_mass(), b.total_mass(), "total_mass");
    let (pa, pb) = (a.total_momentum(), b.total_momentum());
    for c in 0..2 {
        close(pa[c], pb[c], &format!("momentum[{c}]"));
    }
    let (fa, fb) = (a.probed_force(), b.probed_force());
    for c in 0..2 {
        close(fa[c], fb[c], &format!("probed_force[{c}]"));
    }
}

const DECOMPS: [[usize; 3]; 3] = [[2, 1, 1], [1, 2, 1], [2, 2, 1]];

/// Run `steps`, comparing the split runs against the monolithic baseline at
/// several checkpoints (bit-exact fields at every checkpoint).
fn check_case(case: &Case, steps: usize, init: Option<&dyn Fn(usize, usize) -> (f64, [f64; 3])>, what: &str) {
    let mut base = build(case, [1, 1, 1], LocalPeriodic, false);
    if let Some(f) = init {
        base.init_with(|x, y, _| f(x, y));
    }
    let mut splits: Vec<(String, S<InProcess>)> = Vec::new();
    for d in DECOMPS {
        for two_pass in [false, true] {
            let mut s = build(case, d, InProcess, two_pass);
            if let Some(f) = init {
                s.init_with(|x, y, _| f(x, y));
            }
            splits.push((format!("{what} {d:?} two_pass={two_pass}"), s));
        }
    }
    // t = 0 (init path must already be split-invariant).
    for (name, s) in &splits {
        assert_fields_equal(&base, s, &format!("{name} t=0"));
    }
    let checkpoints = [1, 2, 3, 5, steps / 2, steps];
    let mut t = 0;
    for &cp in checkpoints.iter().filter(|&&c| c > 0) {
        while t < cp {
            base.step();
            for (_, s) in splits.iter_mut() {
                s.step();
            }
            t += 1;
        }
        for (name, s) in &splits {
            assert_fields_equal(&base, s, &format!("{name} t={t}"));
            assert_diagnostics_close(&base, s, 1e-12, &format!("{name} t={t}"));
        }
    }
    println!("{what}: fields bit-exact across {:?} (+two-pass) over {steps} steps", DECOMPS);
}

// ---------------------------------------------------------------------------

#[test]
fn t13_tgv_periodic_split_invariant() {
    let case = Case {
        spec: GlobalSpec {
            dims: [48, 36, 1],
            nu: 0.02,
            periodic: [true, true, false],
            ..Default::default()
        },
        walls: WallSpec::default(),
    };
    let init = |x: usize, y: usize| {
        let kx = 2.0 * std::f64::consts::PI / 48.0;
        let ky = 2.0 * std::f64::consts::PI / 36.0;
        let u0 = 0.04;
        (
            1.0 + 0.01 * (kx * x as f64).cos() * (ky * y as f64).cos(),
            [
                -u0 * (kx * x as f64).cos() * (ky * y as f64).sin(),
                u0 * (kx * x as f64).sin() * (ky * y as f64).cos(),
                0.0,
            ],
        )
    };
    check_case(&case, 200, Some(&init), "TGV");
}

#[test]
fn t13_cavity_split_invariant() {
    let mut walls = WallSpec::default();
    for f in [Face::XNeg, Face::XPos, Face::YNeg, Face::YPos] {
        walls.is_wall[f.index()] = true;
    }
    walls.u[Face::YPos.index()] = [0.1, 0.0, 0.0];
    let case = Case {
        spec: GlobalSpec {
            dims: [40, 40, 1],
            nu: 0.02,
            periodic: [false, false, false],
            ..Default::default()
        },
        walls,
    };
    check_case(&case, 200, None, "cavity");
}

#[test]
fn t13_channel_profile_outflow_split_invariant() {
    // Open faces + per-node inlet profile distributed across parts: the
    // profile slicing and BC-ownership logic must be split-invariant too.
    let mut walls = WallSpec::default();
    walls.is_wall[Face::YNeg.index()] = true;
    walls.is_wall[Face::YPos.index()] = true;
    let mut faces = [FaceBC::Closed; 6];
    faces[Face::XNeg.index()] = FaceBC::Velocity {
        u: [0.05, 0.0, 0.0],
    };
    faces[Face::XPos.index()] = FaceBC::Convective { u_conv: 0.05 };
    let case = Case {
        spec: GlobalSpec {
            dims: [64, 32, 1],
            nu: 0.05,
            periodic: [false, false, false],
            faces,
            ..Default::default()
        },
        walls,
    };
    let ny = 32usize;
    let profile: Vec<[f64; 3]> = (0..ny)
        .map(|y| {
            let yy = y as f64 / (ny - 1) as f64;
            [0.08 * 4.0 * yy * (1.0 - yy), 0.0, 0.0]
        })
        .collect();
    let mut base = build(&case, [1, 1, 1], LocalPeriodic, false);
    base.set_inlet_profile(Face::XNeg, &profile);
    let mut splits: Vec<(String, S<InProcess>)> = Vec::new();
    for d in DECOMPS {
        let mut s = build(&case, d, InProcess, false);
        s.set_inlet_profile(Face::XNeg, &profile);
        splits.push((format!("channel {d:?}"), s));
    }
    for t in 1..=200 {
        base.step();
        for (_, s) in splits.iter_mut() {
            s.step();
        }
        if t <= 3 || t % 50 == 0 {
            for (name, s) in &splits {
                assert_fields_equal(&base, s, &format!("{name} t={t}"));
                assert_diagnostics_close(&base, s, 1e-12, &format!("{name} t={t}"));
            }
        }
    }
    println!("channel: fields bit-exact across {DECOMPS:?} over 200 steps");
}

#[test]
fn t13_cylinder_probe_split_invariant() {
    // Obstacle + probe crossing part boundaries in a 2x2 split: bounce-back
    // and momentum-exchange bookkeeping must not care where the seams are.
    let (nx, ny) = (64usize, 40usize);
    let mut walls = WallSpec::default();
    walls.is_wall[Face::YNeg.index()] = true;
    walls.is_wall[Face::YPos.index()] = true;
    let mut faces = [FaceBC::Closed; 6];
    faces[Face::XNeg.index()] = FaceBC::Velocity {
        u: [0.05, 0.0, 0.0],
    };
    faces[Face::XPos.index()] = FaceBC::Outflow;
    let case = Case {
        spec: GlobalSpec {
            dims: [nx, ny, 1],
            nu: 0.02,
            periodic: [false, false, false],
            faces,
            ..Default::default()
        },
        walls,
    };
    // Cylinder centred on the 2x2 seam so its cells span all four parts.
    let (cx, cy, r) = (nx as f64 / 2.0, ny as f64 / 2.0 - 0.3, 5.4);
    let inside = move |x: usize, y: usize| {
        let (dx, dy) = (x as f64 - cx, y as f64 - cy);
        dx * dx + dy * dy < r * r
    };
    let mut base = build(&case, [1, 1, 1], LocalPeriodic, false);
    let mut split = build(&case, [2, 2, 1], InProcess, true);
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
    for t in 1..=250 {
        base.step();
        split.step();
        let (fa, fb) = (base.probed_force(), split.probed_force());
        for c in 0..2 {
            let d = (fa[c] - fb[c]).abs();
            assert!(
                d <= 1e-12 + 1e-12 * fa[c].abs(),
                "t={t}: probed_force[{c}] Δ = {d:e} (ref {:e})",
                fa[c]
            );
        }
        if t % 50 == 0 || t == 250 {
            assert_fields_equal(&base, &split, &format!("cylinder t={t}"));
        }
    }
    println!("cylinder probe: split-invariant over 250 steps (2x2, two-pass)");
}

#[test]
fn t13_tgv3d_2x2x2_split_invariant() {
    // 3D decomposition (COMPETITIVE_SPEC T13: "+3D: 2x2x2"): exercises the
    // z exchange phase and 3D edge/corner halo forwarding, which no 2D case
    // touches. D3Q19, triply periodic TGV.
    use lbm_core::lattice::D3Q19;
    type S3<H> = Solver<D3Q19, f64, CpuScalar, H>;
    let n = 16usize;
    let spec = GlobalSpec::<f64> {
        dims: [n, n, n],
        nu: 0.02,
        periodic: [true, true, true],
        ..Default::default()
    };
    let k = 2.0 * std::f64::consts::PI / n as f64;
    let init = move |x: usize, y: usize, z: usize| {
        let (xx, yy, zz) = (k * x as f64, k * y as f64, k * z as f64);
        (
            1.0,
            [
                0.03 * xx.sin() * yy.cos() * zz.cos(),
                -0.03 * xx.cos() * yy.sin() * zz.cos(),
                0.0,
            ],
        )
    };
    for (decomp, two_pass) in [([2, 2, 2], false), ([2, 2, 2], true), ([2, 1, 2], false)] {
        let mut s: S3<InProcess> =
            Solver::new(&spec, &[], &[], decomp, CpuScalar::default(), InProcess);
        s.set_two_pass(two_pass);
        s.init_with(init);
        let mut b: S3<LocalPeriodic> = Solver::new(
            &spec,
            &[],
            &[],
            [1, 1, 1],
            CpuScalar::default(),
            LocalPeriodic,
        );
        b.init_with(init);
        for _ in 0..100 {
            s.step();
            b.step();
        }
        for (name, va, vb) in [
            ("rho", b.gather_rho(), s.gather_rho()),
            ("ux", b.gather_ux(), s.gather_ux()),
            ("uy", b.gather_uy(), s.gather_uy()),
            ("uz", b.gather_uz(), s.gather_uz()),
        ] {
            let d = va
                .iter()
                .zip(&vb)
                .map(|(x, y)| (x - y).abs())
                .fold(0.0f64, f64::max);
            assert_eq!(
                d, 0.0,
                "3D {decomp:?} two_pass={two_pass}: {name} differs by {d:e}"
            );
        }
        for q in 0..19 {
            let (fa, fb) = (b.gather_f(q), s.gather_f(q));
            let d = fa
                .iter()
                .zip(&fb)
                .map(|(x, y)| (x - y).abs())
                .fold(0.0f64, f64::max);
            assert_eq!(
                d, 0.0,
                "3D {decomp:?} two_pass={two_pass}: f[{q}] differs by {d:e}"
            );
        }
        let dm = (b.total_mass() - s.total_mass()).abs();
        assert!(dm <= 1e-12 * b.total_mass(), "3D mass Δ = {dm:e}");
    }
    println!("TGV3D 16^3: bit-exact across 2x2x2 (+two-pass) and 2x1x2 over 100 steps");
}

#[test]
fn t13_shan_chen_droplet_native_split_invariant() {
    // Single-component Shan–Chen droplet sitting on the 2x2 corner, driven
    // by the native Solver::update_shan_chen_force (ψ halo exchange via
    // exchange_scalar) — closes the TESTING_NOTES gap where T13 could only
    // hand-roll per-cell forces. Fields must match bit-for-bit.
    let n = 48usize;
    let case = Case {
        spec: GlobalSpec {
            dims: [n, n, 1],
            nu: 1.0 / 6.0,
            periodic: [true, true, false],
            ..Default::default()
        },
        walls: WallSpec::default(),
    };
    // Smooth droplet: liquid ~1.9 inside r0, vapour ~0.16 outside (classic
    // ψ = 1 - exp(-rho), G = -5 two-phase state).
    let init = move |x: usize, y: usize| {
        let (dx, dy) = (x as f64 + 0.5 - n as f64 / 2.0, y as f64 + 0.5 - n as f64 / 2.0);
        let r = (dx * dx + dy * dy).sqrt();
        let rho = 0.16 + (1.90 - 0.16) * 0.5 * (1.0 - ((r - 10.0) / 2.0).tanh());
        (rho, [0.0, 0.0, 0.0])
    };
    let g = -5.0f64;
    let psi = |rho: f64| 1.0 - (-rho).exp();
    let mut base = build(&case, [1, 1, 1], LocalPeriodic, false);
    base.init_with(|x, y, _| init(x, y));
    for decomp in DECOMPS {
        let mut s = build(&case, decomp, InProcess, false);
        s.init_with(|x, y, _| init(x, y));
        let mut b = build(&case, [1, 1, 1], LocalPeriodic, false);
        b.init_with(|x, y, _| init(x, y));
        for t in 0..150 {
            b.update_shan_chen_force(g, psi);
            s.update_shan_chen_force(g, psi);
            b.step();
            s.step();
            if t < 3 || t % 50 == 49 {
                assert_fields_equal(&b, &s, &format!("shan-chen droplet {decomp:?} t={}", t + 1));
                assert_diagnostics_close(&b, &s, 1e-12, &format!("shan-chen droplet {decomp:?} t={}", t + 1));
            }
        }
    }
    println!("Shan-Chen droplet (native update_shan_chen_force): bit-exact across {DECOMPS:?}");
}

#[test]
fn t13_shan_chen_wall_adhesion_native_matches_compat_and_split() {
    // Wall-adhesion Shan–Chen (g_wall and virtual wall density wall_rho):
    // the native Solver::update_shan_chen_force_with_walls must
    //   (a) reproduce the facade compat::ShanChen::update_force bit-exactly
    //       on a monolithic domain (the facade carries the V1 wetting
    //       numerics that T11b/T11c freeze), and
    //   (b) stay bit-exact under decomposition (solid rims cross the seams,
    //       so halo-solid adhesion is exercised).
    use lbm_core::compat::multiphase::ShanChen as CompatShanChen;
    use lbm_core::compat::prelude::{
        EdgeBC as CEdgeBC, Edges as CEdges, SimConfig as CSimConfig, Simulation as CompatSim,
    };

    let (nx, ny) = (64usize, 40usize);
    let g = -5.0f64;
    let psi = |rho: f64| 1.0 - (-rho).exp();
    // Droplet resting on the bottom wall rim (T11b geometry, shrunk).
    let (cx, cy, r0) = (nx as f64 / 2.0, 9.0, 11.0);
    let init2 = move |x: usize, y: usize| {
        if y == 0 || y + 1 == ny {
            return (0.15, 0.0, 0.0);
        }
        let d = ((x as f64 - cx).powi(2) + (y as f64 - cy).powi(2)).sqrt();
        (if d < r0 { 2.0 } else { 0.15 }, 0.0, 0.0)
    };
    let steps = 150usize;

    // (g_wall, wall_rho): wetting, de-wetting, virtual-wall-density control.
    for (g_wall, wall_rho) in [(-1.5f64, None), (0.9, None), (0.0, Some(1.2f64))] {
        let psi_wall = wall_rho.map_or(0.0, psi);
        let what = format!("g_wall={g_wall} wall_rho={wall_rho:?}");

        // Facade reference (V1 numerics, proven by the synced T11 suite).
        let mut sim: CompatSim<f64> = CSimConfig {
            nx,
            ny,
            nu: 1.0 / 6.0,
            edges: CEdges {
                left: CEdgeBC::Periodic,
                right: CEdgeBC::Periodic,
                bottom: CEdgeBC::BounceBack,
                top: CEdgeBC::BounceBack,
            },
            ..Default::default()
        }
        .build()
        .unwrap();
        sim.init_with(init2);
        let mut sc = CompatShanChen::new(g).with_wall(g_wall);
        if let Some(r) = wall_rho {
            sc = sc.with_wall_rho(r);
        }

        // Native monolithic + split runs of the same scenario.
        let case = Case {
            spec: GlobalSpec {
                dims: [nx, ny, 1],
                nu: 1.0 / 6.0,
                periodic: [true, false, false],
                ..Default::default()
            },
            walls: WallSpec {
                is_wall: {
                    let mut w = [false; 6];
                    w[Face::YNeg.index()] = true;
                    w[Face::YPos.index()] = true;
                    w
                },
                ..Default::default()
            },
        };
        let mut mono = build(&case, [1, 1, 1], LocalPeriodic, false);
        mono.init_with(|x, y, _| {
            let (r, ux, uy) = init2(x, y);
            (r, [ux, uy, 0.0])
        });
        let mut splits: Vec<S<InProcess>> = DECOMPS
            .iter()
            .map(|&d| {
                let mut s = build(&case, d, InProcess, false);
                s.init_with(|x, y, _| {
                    let (r, ux, uy) = init2(x, y);
                    (r, [ux, uy, 0.0])
                });
                s
            })
            .collect();

        for t in 0..steps {
            sc.update_force(&mut sim);
            sim.step();
            mono.update_shan_chen_force_with_walls(g, g_wall, psi_wall, psi);
            mono.step();
            for s in splits.iter_mut() {
                s.update_shan_chen_force_with_walls(g, g_wall, psi_wall, psi);
                s.step();
            }
            if t < 3 || t % 50 == 49 || t + 1 == steps {
                // (a) native == facade, bit-exact.
                let pairs = [
                    ("rho", sim.rho_field().to_vec(), mono.gather_rho()),
                    ("ux", sim.ux_field().to_vec(), mono.gather_ux()),
                    ("uy", sim.uy_field().to_vec(), mono.gather_uy()),
                ];
                for (name, va, vb) in pairs {
                    let d = va
                        .iter()
                        .zip(&vb)
                        .map(|(x, y)| (x - y).abs())
                        .fold(0.0f64, f64::max);
                    assert_eq!(d, 0.0, "{what} t={}: native vs compat {name}", t + 1);
                }
                // (b) split-invariant.
                for (s, d) in splits.iter().zip(DECOMPS.iter()) {
                    assert_fields_equal(&mono, s, &format!("{what} {d:?} t={}", t + 1));
                    assert_diagnostics_close(&mono, s, 1e-12, &format!("{what} {d:?} t={}", t + 1));
                }
            }
        }
    }
    println!("Shan-Chen wall adhesion (g_wall/wall_rho): native == compat bit-exact, split-invariant");
}

#[test]
fn t13_uneven_split_and_deeper_decomp() {
    // Remainder handling: 50 cells over 3 parts (17/17/16) and a 4x1 strip
    // decomposition, against the monolithic baseline.
    let case = Case {
        spec: GlobalSpec {
            dims: [50, 30, 1],
            nu: 0.03,
            periodic: [true, true, false],
            ..Default::default()
        },
        walls: WallSpec::default(),
    };
    let init = |x: usize, y: usize| {
        let kx = 2.0 * std::f64::consts::PI / 50.0;
        let ky = 2.0 * std::f64::consts::PI / 30.0;
        (
            1.0,
            [
                0.03 * (ky * y as f64).sin(),
                0.03 * (kx * x as f64).sin(),
                0.0,
            ],
        )
    };
    let mut base = build(&case, [1, 1, 1], LocalPeriodic, false);
    base.init_with(|x, y, _| init(x, y));
    for decomp in [[3, 1, 1], [4, 1, 1], [3, 2, 1]] {
        let mut s = build(&case, decomp, InProcess, false);
        s.init_with(|x, y, _| init(x, y));
        let mut b = build(&case, [1, 1, 1], LocalPeriodic, false);
        b.init_with(|x, y, _| init(x, y));
        for _ in 0..120 {
            s.step();
            b.step();
        }
        assert_fields_equal(&b, &s, &format!("uneven {decomp:?}"));
        assert_diagnostics_close(&b, &s, 1e-12, &format!("uneven {decomp:?}"));
    }
    println!("uneven splits (3x1, 4x1, 3x2 over 50x30): bit-exact");
}
