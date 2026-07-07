use lbm_core::prelude::*;

fn periodic_spec(nx: usize, ny: usize, nu: f64) -> GlobalSpec<f64> {
    GlobalSpec {
        dims: [nx, ny, 1],
        nu,
        periodic: [true, true, false],
        collision: CollisionKind::Trt {
            magic: CollisionKind::MAGIC_STD,
        },
        ..Default::default()
    }
}

fn periodic_solver(
    nx: usize,
    ny: usize,
    decomp: [usize; 3],
) -> Solver<D2Q9, f64, CpuScalar, InProcess> {
    let n = nx * ny;
    Solver::new(
        &periodic_spec(nx, ny, 1.0 / 6.0),
        &vec![false; n],
        &vec![[0.0; 3]; n],
        decomp,
        CpuScalar::default(),
        InProcess,
    )
}

fn walled_channel(nx: usize, ny: usize, top_u: f64) -> Solver<D2Q9, f64, CpuScalar, LocalPeriodic> {
    let mut walls = WallSpec::default();
    walls.is_wall[Face::YNeg.index()] = true;
    walls.is_wall[Face::YPos.index()] = true;
    walls.u[Face::YPos.index()] = [top_u, 0.0, 0.0];
    let spec = GlobalSpec {
        dims: [nx, ny, 1],
        nu: 1.0 / 6.0,
        periodic: [true, false, false],
        collision: CollisionKind::Trt {
            magic: CollisionKind::MAGIC_STD,
        },
        ..Default::default()
    };
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

fn ibm_couette_body(nx: usize, y: f64, speed: f64) -> RotatingBody {
    let markers = (0..nx)
        .map(|x| IbmMarker {
            position: [x as f64, y, 0.0],
            weight: 1.0,
        })
        .collect();
    // A far-away center with omega_z chosen so Omega x r gives nearly uniform
    // +x velocity on the marker line; the tiny y-variation is below test bands.
    let omega = 1.0e-6;
    let cy = y + speed / omega;
    RotatingBody::from_markers([0.0, cy, 0.0], [0.0, 0.0, omega], markers)
}

#[test]
fn direct_forcing_iterations_reduce_marker_slip_and_conserve_momentum() {
    let mut one = periodic_solver(48, 48, [1, 1, 1]);
    let body = RotatingBody::circle_2d([24.0, 24.0], 8.0, 0.003, 96);
    let single = one.apply_rotating_ibm(
        &body,
        DirectForcingConfig {
            max_iterations: 1,
            slip_tolerance: 0.0,
            kernel_radius: 2,
            relaxation: 1.0,
        },
    );

    let mut multi = periodic_solver(48, 48, [1, 1, 1]);
    let iterated = multi.apply_rotating_ibm(
        &body,
        DirectForcingConfig {
            max_iterations: 4,
            slip_tolerance: 1.0e-3,
            kernel_radius: 2,
            relaxation: 1.0,
        },
    );
    println!(
        "IBM rotating cylinder force update: single slip_max_rel={:.6e}, multi slip_max_rel={:.6e}, slip_rms_rel={:.6e}, momentum_error_rel={:.6e}, torque_z={:.6e}",
        single.slip_max_rel,
        iterated.slip_max_rel,
        iterated.slip_rms_rel,
        iterated.momentum_error_rel,
        iterated.torque[2]
    );
    assert!(iterated.iterations > 1);
    assert!(iterated.slip_max_rel < single.slip_max_rel * 0.13);
    assert!(iterated.slip_max_rel < 3.0e-3);
    assert!(iterated.slip_rms_rel < 2.0e-3);
    assert!(iterated.momentum_error_rel < 1.0e-12);
}

#[test]
fn rotating_cylinder_torque_flips_with_omega_sign() {
    let body_pos = RotatingBody::circle_2d([24.0, 24.0], 8.0, 0.003, 96);
    let body_neg = RotatingBody::circle_2d([24.0, 24.0], 8.0, -0.003, 96);
    let cfg = DirectForcingConfig {
        max_iterations: 4,
        slip_tolerance: 1.0e-3,
        kernel_radius: 2,
        relaxation: 1.0,
    };

    let mut pos = periodic_solver(48, 48, [1, 1, 1]);
    let mut neg = periodic_solver(48, 48, [1, 1, 1]);
    let d_pos = pos.apply_rotating_ibm(&body_pos, cfg);
    let d_neg = neg.apply_rotating_ibm(&body_neg, cfg);

    println!(
        "IBM torque sign symmetry: torque_pos={:.6e}, torque_neg={:.6e}, slip_pos={:.6e}, slip_neg={:.6e}",
        d_pos.torque[2], d_neg.torque[2], d_pos.slip_max_rel, d_neg.slip_max_rel
    );
    assert!(
        d_pos.torque[2] * d_neg.torque[2] < 0.0,
        "opposite rotations must produce opposite torque signs"
    );
    let rel = (d_pos.torque[2] + d_neg.torque[2]).abs() / d_pos.torque[2].abs();
    assert!(rel < 1.0e-12, "torque sign symmetry rel={rel:e}");
    assert!(
        (d_pos.slip_max_rel - d_neg.slip_max_rel).abs() < 1.0e-15,
        "opposite rotations should have identical slip magnitudes"
    );
}

#[test]
fn marker_straddling_partition_seam_is_bit_identical() {
    let body = RotatingBody::circle_2d([32.0, 32.0], 9.0, 0.0025, 120);
    let cfg = DirectForcingConfig {
        max_iterations: 3,
        slip_tolerance: 1.0e-4,
        kernel_radius: 1,
        relaxation: 1.0,
    };
    let mut mono = periodic_solver(64, 64, [1, 1, 1]);
    let mut split = periodic_solver(64, 64, [2, 2, 1]);
    for _ in 0..20 {
        mono.clear_body_force_field();
        split.clear_body_force_field();
        let da = mono.apply_rotating_ibm(&body, cfg);
        let db = split.apply_rotating_ibm(&body, cfg);
        assert_eq!(da.iterations, db.iterations);
        mono.step();
        split.step();
    }
    let ax = mono.gather_ux();
    let bx = split.gather_ux();
    let ay = mono.gather_uy();
    let by = split.gather_uy();
    let max_diff = ax
        .iter()
        .zip(&bx)
        .chain(ay.iter().zip(&by))
        .map(|(a, b)| (a - b).abs())
        .fold(0.0f64, f64::max);
    println!("IBM T13 straddle max velocity diff = {max_diff:.6e}");
    assert_eq!(max_diff, 0.0);
}

#[test]
fn ibm_moving_wall_couette_matches_native_moving_wall_profile() {
    let (nx, ny, top_u) = (48, 34, 0.002);
    let mut native = walled_channel(nx, ny, top_u);
    native.run(3000);

    let mut ibm = walled_channel(nx, ny, 0.0);
    let body = ibm_couette_body(nx, (ny - 3) as f64, top_u);
    let cfg = DirectForcingConfig {
        max_iterations: 1,
        slip_tolerance: 1.0,
        kernel_radius: 1,
        relaxation: 0.05,
    };
    let mut last = IbmDiagnostics::default();
    for _ in 0..800 {
        ibm.clear_body_force_field();
        last = ibm.apply_rotating_ibm(&body, cfg);
        ibm.step();
    }
    let a = native.gather_ux();
    let b = ibm.gather_ux();
    let mut l2 = 0.0;
    let mut ref2 = 0.0;
    let mut linf = 0.0f64;
    for y in 2..ny - 3 {
        let ia = y * nx + nx / 2;
        let d = b[ia] - a[ia];
        l2 += d * d;
        ref2 += a[ia] * a[ia];
        linf = linf.max(d.abs());
    }
    let l2_rel = (l2 / ref2).sqrt();
    let linf_rel = linf / top_u;
    println!(
        "IBM Couette vs native moving wall: slip_max_rel={:.6e}, slip_rms_rel={:.6e}, L2_rel={:.6e}, Linf/U={:.6e}",
        last.slip_max_rel, last.slip_rms_rel, l2_rel, linf_rel
    );
    assert!(last.slip_max_rel < 1.0e-2);
    assert!(l2_rel < 0.65);
    assert!(linf_rel < 0.55);
}

#[test]
fn taylor_couette_marker_profile_matches_analytic_band() {
    let (nx, ny) = (80, 80);
    let center = [40.0, 40.0];
    let r_i = 10.0;
    let r_o = 30.0;
    let omega = 0.00015;
    let body = RotatingBody::circle_2d(center, r_i, omega, 160);
    let spec = periodic_spec(nx, ny, 1.0 / 6.0);
    let mut solid = vec![false; nx * ny];
    for y in 0..ny {
        for x in 0..nx {
            let dx = x as f64 - center[0];
            let dy = y as f64 - center[1];
            let r = (dx * dx + dy * dy).sqrt();
            solid[y * nx + x] = r < r_i - 1.5 || r > r_o + 0.5;
        }
    }
    let mut sim: Solver<D2Q9, f64, CpuScalar, InProcess> = Solver::new(
        &spec,
        &solid,
        &vec![[0.0; 3]; nx * ny],
        [1, 1, 1],
        CpuScalar::default(),
        InProcess,
    );
    let cfg = DirectForcingConfig {
        max_iterations: 1,
        slip_tolerance: 1.0,
        kernel_radius: 1,
        relaxation: 0.05,
    };
    let mut last = IbmDiagnostics::default();
    for _ in 0..800 {
        sim.clear_body_force_field();
        last = sim.apply_rotating_ibm(&body, cfg);
        sim.step();
    }
    let ux = sim.gather_ux();
    let uy = sim.gather_uy();
    let ui = omega * r_i;
    let a = -ui * r_i * r_i / (r_o * r_o - r_i * r_i);
    let b = ui * r_i * r_i * r_o * r_o / (r_o * r_o - r_i * r_i);
    let mut l2 = 0.0;
    let mut ref2 = 0.0;
    let mut linf = 0.0f64;
    let mut count = 0usize;
    for y in 0..ny {
        for x in 0..nx {
            let dx = x as f64 - center[0];
            let dy = y as f64 - center[1];
            let r = (dx * dx + dy * dy).sqrt();
            if !(r_i + 4.0..=r_o - 5.0).contains(&r) {
                continue;
            }
            let utheta = (-ux[y * nx + x] * dy + uy[y * nx + x] * dx) / r;
            let exact = a * r + b / r;
            let d = utheta - exact;
            l2 += d * d;
            ref2 += exact * exact;
            linf = linf.max(d.abs());
            count += 1;
        }
    }
    let l2_rel = (l2 / ref2).sqrt();
    let linf_rel = linf / ui;
    println!(
        "IBM Taylor-Couette: samples={count}, slip_max_rel={:.6e}, slip_rms_rel={:.6e}, torque_z={:.6e}, momentum_error_rel={:.6e}, L2_rel={:.6e}, Linf/U_i={:.6e}",
        last.slip_max_rel,
        last.slip_rms_rel,
        last.torque[2],
        last.momentum_error_rel,
        l2_rel,
        linf_rel
    );
    assert!(last.slip_max_rel < 1.0e-2);
    assert!(last.momentum_error_rel < 1.0e-10);
    assert!(l2_rel < 0.95);
    assert!(linf_rel < 5.8);
}
