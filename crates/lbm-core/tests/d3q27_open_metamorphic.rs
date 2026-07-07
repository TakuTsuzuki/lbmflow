//! Adversarial metamorphic gates for D3Q27 open-face NEBB.
//!
//! These tests do not validate a new flow model. They pin coordinate
//! equivariance, mirror equivariance, sign reversal, open-kind acceptance
//! coverage, and D3Q27 GPU open-face CPU-vs-GPU equivalence.

use lbm_core::prelude::*;
use std::f64::consts::PI;

type S27 = Solver<D3Q27, f64, CpuScalar, LocalPeriodic>;

const U_PEAK: f64 = 0.032;
const STEPS: usize = 10_000;
const EQUIVARIANCE_TOL: f64 = 1.0e-10;

fn duct_shape(c1: usize, c2: usize, n1: usize, n2: usize) -> f64 {
    if c1 == 0 || c2 == 0 || c1 + 1 == n1 || c2 + 1 == n2 {
        return 0.0;
    }
    let h = (n1 - 2) as f64;
    let w = (n2 - 2) as f64;
    let (a, b) = (h / 2.0, w / 2.0);
    let yy = c1 as f64 - 0.5;
    let zt = c2 as f64 - 0.5 - b;
    let pref = 16.0 * a * a / (PI * PI * PI);
    let mut sum = 0.0;
    let mut n = 1;
    while n <= 99 {
        let nf = n as f64;
        let kn = nf * PI / (2.0 * a);
        let ratio =
            ((kn * zt.abs()).exp() + (-kn * zt.abs()).exp()) / ((kn * b).exp() + (-kn * b).exp());
        sum += (1.0 - ratio) * (kn * yy).sin() / (nf * nf * nf);
        n += 2;
    }
    pref * sum
}

fn profile_speed(c1: usize, c2: usize, n1: usize, n2: usize, u_peak: f64, skew: f64) -> f64 {
    let center = duct_shape(n1 / 2, n2 / 2, n1, n2);
    let base = u_peak * duct_shape(c1, c2, n1, n2) / center;
    let skew_factor = 1.0 + skew * (2.0 * c1 as f64 / (n1 - 1) as f64 - 1.0);
    base * skew_factor
}

fn velocity(axis: usize, speed: f64) -> [f64; 3] {
    let mut u = [0.0; 3];
    u[axis] = speed;
    u
}

fn duct_solver(dims: [usize; 3], axis: usize, inlet: Face, u_peak: f64, skew: f64) -> S27 {
    assert_eq!(inlet.axis(), axis);
    let outlet = inlet.opposite();
    let sign = if inlet.is_neg() { 1.0 } else { -1.0 };
    let mut walls = WallSpec::<f64>::default();
    for face in Face::ALL {
        if face.axis() != axis {
            walls.is_wall[face.index()] = true;
        }
    }
    let mut faces = [FaceBC::Closed; 6];
    faces[inlet.index()] = FaceBC::Velocity {
        u: velocity(axis, sign * u_peak),
    };
    faces[outlet.index()] = FaceBC::Pressure { rho: 1.0 };
    let spec = GlobalSpec::<f64> {
        dims,
        nu: 0.05,
        periodic: [false, false, false],
        faces,
        ..Default::default()
    };
    let (solid, wall_u) = build_wall_rims(3, dims, &walls);
    let mut s = S27::new(
        &spec,
        &solid,
        &wall_u,
        [1, 1, 1],
        CpuScalar::default(),
        LocalPeriodic,
    );
    let (t1, t2) = inlet.tangents();
    s.set_inlet_profile_with(inlet, |c1, c2| {
        velocity(
            axis,
            sign * profile_speed(c1, c2, dims[t1], dims[t2], u_peak, skew),
        )
    });
    s.init_with(|x, y, z| {
        let p = [x, y, z];
        (
            1.0,
            velocity(
                axis,
                sign * profile_speed(p[t1], p[t2], dims[t1], dims[t2], u_peak, skew),
            ),
        )
    });
    s
}

fn mean_rho_plane(s: &S27, axis: usize, plane: usize) -> f64 {
    let dims = s.dims();
    let mut sum = 0.0;
    let mut n = 0.0;
    for z in 1..dims[2] - 1 {
        for y in 1..dims[1] - 1 {
            for x in 1..dims[0] - 1 {
                let p = [x, y, z];
                if p[axis] == plane {
                    sum += s.rho(x, y, z);
                    n += 1.0;
                }
            }
        }
    }
    sum / n
}

fn assert_duct_behavior(s: &S27, axis: usize, sign: f64, label: &str) {
    let dims = s.dims();
    let lo = 1;
    let mid = dims[axis] / 2;
    let hi = dims[axis] - 2;
    let rho_lo = mean_rho_plane(s, axis, lo);
    let rho_mid = mean_rho_plane(s, axis, mid);
    let rho_hi = mean_rho_plane(s, axis, hi);
    if sign > 0.0 {
        assert!(
            rho_lo > rho_mid && rho_mid > rho_hi,
            "{label}: pressure must drop monotonically along +axis, rho_lo={rho_lo:e}, rho_mid={rho_mid:e}, rho_hi={rho_hi:e}"
        );
    } else {
        assert!(
            rho_hi > rho_mid && rho_mid > rho_lo,
            "{label}: pressure must drop monotonically along -axis, rho_hi={rho_hi:e}, rho_mid={rho_mid:e}, rho_lo={rho_lo:e}"
        );
    }

    let ux = s.gather_ux();
    let uy = s.gather_uy();
    let uz = s.gather_uz();
    let fields = [&ux, &uy, &uz];
    let mut primary_min = f64::INFINITY;
    let mut primary_max = 0.0f64;
    let mut transverse_max = 0.0f64;
    for z in 1..dims[2] - 1 {
        for y in 1..dims[1] - 1 {
            for x in 1..dims[0] - 1 {
                let i = (z * dims[1] + y) * dims[0] + x;
                let up = sign * fields[axis][i];
                primary_min = primary_min.min(up);
                primary_max = primary_max.max(up.abs());
                for (a, field) in fields.iter().enumerate() {
                    if a != axis {
                        transverse_max = transverse_max.max(field[i].abs());
                    }
                }
            }
        }
    }
    assert!(
        primary_min >= -1.0e-12,
        "{label}: primary velocity changed sign, sign-adjusted min={primary_min:e}"
    );
    let transverse_rel = transverse_max / primary_max.max(1.0e-12);
    // Same unidirectional-flow anchor scale as
    // d3q27_open_bc.rs::d3q27_open_duct_matches_series_shape_and_d3q19
    // (`cross_rel <= 1e-3`) for the same Fourier duct inlet profile.
    assert!(
        transverse_rel <= 1.0e-3,
        "{label}: transverse velocity ratio {transverse_rel:e} too large for unidirectional duct, transverse_max={transverse_max:e}, primary_max={primary_max:e}"
    );
}

fn max_rotation_delta(base: &S27, y_axis: &S27, z_axis: &S27) -> f64 {
    let [nx, ny, nz] = base.dims();
    let mut max_delta = 0.0f64;
    for z in 0..nz {
        for y in 0..ny {
            for x in 0..nx {
                let u = base.u(x, y, z);
                let uy = y_axis.u(y, x, z);
                let uz = z_axis.u(y, z, x);
                let mapped_y = [uy[1], uy[0], uy[2]];
                let mapped_z = [uz[2], uz[0], uz[1]];
                for a in 0..3 {
                    max_delta = max_delta.max((u[a] - mapped_y[a]).abs());
                    max_delta = max_delta.max((u[a] - mapped_z[a]).abs());
                }
            }
        }
    }
    max_delta
}

#[test]
fn d3q27_open_duct_velocity_field_rotates_between_axes() {
    let mut x_axis = duct_solver([18, 10, 8], 0, Face::XNeg, U_PEAK, 0.0);
    let mut y_axis = duct_solver([10, 18, 8], 1, Face::YNeg, U_PEAK, 0.0);
    let mut z_axis = duct_solver([10, 8, 18], 2, Face::ZNeg, U_PEAK, 0.0);
    x_axis.run(STEPS);
    y_axis.run(STEPS);
    z_axis.run(STEPS);

    assert_duct_behavior(&x_axis, 0, 1.0, "x-axis duct");
    assert_duct_behavior(&y_axis, 1, 1.0, "y-axis duct");
    assert_duct_behavior(&z_axis, 2, 1.0, "z-axis duct");

    let max_delta = max_rotation_delta(&x_axis, &y_axis, &z_axis);
    // Tolerance basis: validation_cavity.rs::t7_re100_cavity_is_exact_under_four_lid_orientations
    // gates orientation equivariance with L_inf <= 1e-10 after 2000 steps.
    assert!(
        max_delta <= EQUIVARIANCE_TOL,
        "D3Q27 open duct rotation equivariance max |delta|={max_delta:e} > {EQUIVARIANCE_TOL:e}"
    );
}

fn max_y_mirror_delta(a: &S27, b: &S27) -> f64 {
    let [nx, ny, nz] = a.dims();
    let mut max_delta = 0.0f64;
    for z in 0..nz {
        for y in 0..ny {
            for x in 0..nx {
                let u = a.u(x, y, z);
                let m = b.u(x, ny - 1 - y, z);
                let mapped = [m[0], -m[1], m[2]];
                for c in 0..3 {
                    max_delta = max_delta.max((u[c] - mapped[c]).abs());
                }
            }
        }
    }
    max_delta
}

#[test]
fn d3q27_open_duct_mirrors_across_cross_axis_midplane() {
    let dims = [18, 10, 8];
    let mut original = duct_solver(dims, 0, Face::XNeg, U_PEAK, 0.0);
    let mut mirrored = duct_solver(dims, 0, Face::XNeg, U_PEAK, 0.0);
    original.run(STEPS);
    mirrored.run(STEPS);

    assert_duct_behavior(&original, 0, 1.0, "original mirrored duct");
    assert_duct_behavior(&mirrored, 0, 1.0, "reflected mirrored duct");

    let max_delta = max_y_mirror_delta(&original, &mirrored);
    // Tolerance basis: same orientation-equivariance L_inf <= 1e-10 gate as
    // validation_cavity.rs::t7_re100_cavity_is_exact_under_four_lid_orientations.
    assert!(
        max_delta <= EQUIVARIANCE_TOL,
        "D3Q27 open duct y-mirror equivariance max |delta|={max_delta:e} > {EQUIVARIANCE_TOL:e}"
    );
}

fn max_x_reversal_delta(forward: &S27, reverse: &S27) -> f64 {
    let [nx, ny, nz] = forward.dims();
    let mut max_delta = 0.0f64;
    for z in 0..nz {
        for y in 0..ny {
            for x in 0..nx {
                let u = forward.u(x, y, z);
                let r = reverse.u(nx - 1 - x, y, z);
                let mapped = [-r[0], r[1], r[2]];
                for c in 0..3 {
                    max_delta = max_delta.max((u[c] - mapped[c]).abs());
                }
            }
        }
    }
    max_delta
}

#[test]
fn d3q27_open_duct_inlet_profile_sign_anchor_reverses_flow() {
    let dims = [18, 10, 8];
    let mut forward = duct_solver(dims, 0, Face::XNeg, U_PEAK, 0.0);
    let mut reverse = duct_solver(dims, 0, Face::XPos, U_PEAK, 0.0);
    forward.run(STEPS);
    reverse.run(STEPS);

    assert_duct_behavior(&forward, 0, 1.0, "forward duct");
    assert_duct_behavior(&reverse, 0, -1.0, "reverse duct");

    let max_delta = max_x_reversal_delta(&forward, &reverse);
    // Swapping the velocity and pressure faces while reversing the prescribed
    // inlet normal velocity is the x-reflection of this symmetric duct. Vector
    // components therefore map as ux -> -ux while transverse components keep
    // their sign; this is the exact reversed-flow anchor for the admissible
    // velocity-inlet/pressure-outlet D3Q27 configuration.
    assert!(
        max_delta <= EQUIVARIANCE_TOL,
        "D3Q27 open duct sign-reversal max |delta|={max_delta:e} > {EQUIVARIANCE_TOL:e}"
    );
}

#[test]
fn d3q27_open_face_kinds_are_accepted_on_every_face() {
    for face in Face::ALL {
        for bc in [
            FaceBC::Velocity {
                u: velocity(face.axis(), if face.is_neg() { 0.02 } else { -0.02 }),
            },
            FaceBC::Pressure { rho: 1.0 },
            FaceBC::Outflow,
            FaceBC::Convective { u_conv: 0.04 },
        ] {
            let mut faces = [FaceBC::Closed; 6];
            faces[face.index()] = bc;
            faces[face.opposite().index()] = FaceBC::Velocity {
                u: velocity(
                    face.axis(),
                    if face.opposite().is_neg() {
                        0.02
                    } else {
                        -0.02
                    },
                ),
            };
            let mut periodic = [true, true, true];
            periodic[face.axis()] = false;
            let spec = GlobalSpec::<f64> {
                dims: [8, 7, 6],
                periodic,
                faces,
                ..Default::default()
            };
            assert!(
                spec.validate_lattice::<D3Q27>(&[]).is_ok(),
                "D3Q27 should accept {face:?} {bc:?}"
            );
        }
    }
}

#[cfg(feature = "gpu")]
fn gpu_ctx_or_skip() -> Option<std::sync::Arc<GpuContext>> {
    use std::sync::OnceLock;

    static CTX: OnceLock<Result<std::sync::Arc<GpuContext>, String>> = OnceLock::new();
    match CTX.get_or_init(|| GpuContext::new().map_err(|e| e.to_string())) {
        Ok(ctx) => Some(ctx.clone()),
        Err(e) => {
            if std::env::var_os("LBM_REQUIRE_GPU").is_some() {
                panic!("D3Q27 GPU open-face equivalence test requires an adapter: {e}");
            }
            eprintln!("skipping D3Q27 GPU open-face equivalence test: no adapter ({e})");
            None
        }
    }
}

#[cfg(feature = "gpu")]
#[test]
fn d3q27_gpu_open_faces_match_cpu_for_all_supported_kinds() {
    let Some(ctx) = gpu_ctx_or_skip() else {
        return;
    };

    for (outlet, label) in [
        (FaceBC::Pressure { rho: 1.0 }, "pressure"),
        (FaceBC::Outflow, "outflow"),
        (
            FaceBC::Convective {
                u_conv: U_PEAK as f32,
            },
            "convective",
        ),
    ] {
        let dims = [18, 10, 8];
        let mut walls = WallSpec::<f32>::default();
        for face in [Face::YNeg, Face::YPos, Face::ZNeg, Face::ZPos] {
            walls.is_wall[face.index()] = true;
        }
        let mut faces = [FaceBC::Closed; 6];
        faces[Face::XNeg.index()] = FaceBC::Velocity {
            u: [U_PEAK as f32, 0.0, 0.0],
        };
        faces[Face::XPos.index()] = outlet;
        let spec = GlobalSpec::<f32> {
            dims,
            nu: 0.05,
            periodic: [false, false, false],
            faces,
            ..Default::default()
        };
        let (solid, wall_u) = build_wall_rims(3, dims, &walls);
        let mut cpu = Solver::<D3Q27, f32, CpuScalar, LocalPeriodic>::new(
            &spec,
            &solid,
            &wall_u,
            [1, 1, 1],
            CpuScalar::default(),
            LocalPeriodic,
        );
        let mut gpu = Solver::<D3Q27, f32, WgpuBackend<D3Q27>, LocalPeriodic>::new(
            &spec,
            &solid,
            &wall_u,
            [1, 1, 1],
            WgpuBackend::<D3Q27>::new(ctx.clone()),
            LocalPeriodic,
        );
        let [_nx, ny, nz] = dims;
        let profile = (0..nz)
            .flat_map(|z| {
                (0..ny).map(move |y| [profile_speed(y, z, ny, nz, U_PEAK, 0.0) as f32, 0.0, 0.0])
            })
            .collect::<Vec<_>>();
        cpu.set_inlet_profile(Face::XNeg, &profile);
        gpu.set_inlet_profile(Face::XNeg, &profile);
        cpu.init_with(|_, y, z| {
            (
                1.0,
                [profile_speed(y, z, ny, nz, U_PEAK, 0.0) as f32, 0.0, 0.0],
            )
        });
        gpu.init_with(|_, y, z| {
            (
                1.0,
                [profile_speed(y, z, ny, nz, U_PEAK, 0.0) as f32, 0.0, 0.0],
            )
        });
        for _ in 0..3 {
            cpu.run(25);
            gpu.run(25);
        }
        let mut dr = 0.0f64;
        let mut du = 0.0f64;
        let cr = cpu.gather_rho();
        let gr = gpu.gather_rho();
        let cu = [cpu.gather_ux(), cpu.gather_uy(), cpu.gather_uz()];
        let gu = [gpu.gather_ux(), gpu.gather_uy(), gpu.gather_uz()];
        for i in 0..cr.len() {
            dr = dr.max((cr[i] as f64 - gr[i] as f64).abs());
            for a in 0..3 {
                du = du.max((cu[a][i] as f64 - gu[a][i] as f64).abs());
            }
        }
        eprintln!("D3Q27 GPU {label}: rho_abs={dr:.3e} u_abs={du:.3e}");
        assert!(dr <= 1.0e-4, "{label}: D3Q27 GPU rho abs {dr:e}");
        assert!(du <= 1.0e-5, "{label}: D3Q27 GPU ux abs {du:e}");
    }
}
