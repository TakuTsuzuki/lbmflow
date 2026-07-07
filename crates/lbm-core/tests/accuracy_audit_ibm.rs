mod common;

use common::metrics::*;
use lbm_core::prelude::*;

const NU: f64 = 1.0 / 6.0;
const NX_B1: usize = 80;
const NY_B1: usize = 80;
const R_I_B1: f64 = 10.0;
const R_O_B1: f64 = 30.0;
const R_O_EFF_B1: f64 = R_O_B1 + 1.0;
const N_MARKERS_B1: usize = 160;
const OMEGA_MID: f64 = 1.5e-4;

fn converged_cfg() -> DirectForcingConfig {
    DirectForcingConfig {
        max_iterations: 4,
        slip_tolerance: 1.0e-3,
        kernel_radius: 2,
        relaxation: 1.0,
    }
}

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

fn solver_from_solid(
    nx: usize,
    ny: usize,
    nu: f64,
    solid: Vec<bool>,
) -> Solver<D2Q9, f64, CpuScalar, InProcess> {
    Solver::new(
        &periodic_spec(nx, ny, nu),
        &solid,
        &vec![[0.0; 3]; nx * ny],
        [1, 1, 1],
        CpuScalar::default(),
        InProcess,
    )
}

fn b1_outer_wall_solid(nx: usize, ny: usize, center: [f64; 2], r_o: f64) -> Vec<bool> {
    let mut solid = vec![false; nx * ny];
    for y in 0..ny {
        for x in 0..nx {
            let dx = x as f64 - center[0];
            let dy = y as f64 - center[1];
            solid[y * nx + x] = (dx * dx + dy * dy).sqrt() > r_o + 0.5;
        }
    }
    solid
}

fn run_torque_to_steady(
    nx: usize,
    ny: usize,
    center: [f64; 2],
    r_i: f64,
    omega: f64,
    n_markers: usize,
    solid: Vec<bool>,
    cfg: DirectForcingConfig,
    cap: usize,
) -> (f64, IbmDiagnostics, usize) {
    let mut solver = solver_from_solid(nx, ny, NU, solid);
    let body = RotatingBody::circle_2d(center, r_i, omega, n_markers);
    let mut last = IbmDiagnostics::default();
    let mut window = Vec::with_capacity(200);
    let mut previous_avg: Option<f64> = None;
    for step in 1..=cap {
        solver.clear_body_force_field();
        last = solver.apply_rotating_ibm(&body, cfg);
        solver.step();
        window.push(last.torque[2]);
        if window.len() > 200 {
            window.remove(0);
        }
        if step % 200 == 0 && window.len() == 200 {
            let avg = window.iter().sum::<f64>() / window.len() as f64;
            if let Some(prev) = previous_avg {
                let rel = (avg - prev).abs() / avg.abs().max(1.0e-30);
                if rel < 1.0e-4 {
                    return (avg, last, step);
                }
            }
            previous_avg = Some(avg);
        }
    }
    let avg = window.iter().sum::<f64>() / window.len().max(1) as f64;
    (avg, last, cap)
}

fn analytic_couette_torque(omega: f64, r_i: f64, r_o_eff: f64) -> f64 {
    // Steady Stokes annular Couette has u_theta(r) = a r + b/r. Enforcing
    // u_theta(r_i) = Omega r_i and u_theta(r_o_eff) = 0 gives
    // a = -Omega r_i^2/(r_o_eff^2 - r_i^2) and
    // b = Omega r_i^2 r_o_eff^2/(r_o_eff^2 - r_i^2). The shear is
    // tau_rtheta = mu r d(u_theta/r)/dr = -2 mu b/r^2, so the torque per
    // unit depth on the inner cylinder is
    // |T| = 2 pi r_i^2 |tau_rtheta(r_i)|
    //     = 4 pi mu Omega r_i^2 r_o_eff^2/(r_o_eff^2 - r_i^2), with rho=1
    // and mu=rho*nu. For r_o_eff in the staircase half-way ambiguity band
    // [r_o+0.5, r_o+1.5], T changes by 0.41% at r_i=10,r_o=30, below 1%.
    let mu = NU;
    4.0 * std::f64::consts::PI * mu * omega * r_i * r_i * r_o_eff * r_o_eff
        / (r_o_eff * r_o_eff - r_i * r_i)
}

#[test]
fn b1_steady_rotating_cylinder_torque_matches_annular_couette_canary() {
    let center = [40.0, 40.0];
    let solid = b1_outer_wall_solid(NX_B1, NY_B1, center, R_O_B1);
    let (measured, diag, steps) = run_torque_to_steady(
        NX_B1,
        NY_B1,
        center,
        R_I_B1,
        OMEGA_MID,
        N_MARKERS_B1,
        solid,
        converged_cfg(),
        6000,
    );
    let measured_abs = measured.abs();
    let reference = analytic_couette_torque(OMEGA_MID, R_I_B1, R_O_EFF_B1);
    let rel = (measured_abs - reference).abs() / reference;
    println!(
        "ACC IBM B1: omega={OMEGA_MID:.6e}, steps={steps}, measured_T={measured:.6e}, |T|={measured_abs:.6e}, T_analytic={reference:.6e}, ratio={:.6e}, rel_err={rel:.6e}, slip_max_rel={:.6e}",
        measured_abs / reference,
        diag.slip_max_rel
    );
    assert!(
        rel <= 0.25,
        "ACC IBM B1 measured rel_err={rel:.6e} exceeds band 2.5e-1 using |T_measured-T_analytic|/T_analytic; measured_T={measured:.6e}, T_analytic={reference:.6e}"
    );
}

#[test]
#[ignore = "heavy ACC-AUDIT B1 torque/Omega linear sweep"]
fn b1_steady_rotating_cylinder_torque_linear_heavy_sweep() {
    let center = [40.0, 40.0];
    let omegas = [0.75e-4, 1.5e-4, 3.0e-4];
    let mut samples = Vec::new();
    let mut measured = Vec::new();
    for omega in omegas {
        let solid = b1_outer_wall_solid(NX_B1, NY_B1, center, R_O_B1);
        let (torque, _diag, steps) = run_torque_to_steady(
            NX_B1,
            NY_B1,
            center,
            R_I_B1,
            omega,
            N_MARKERS_B1,
            solid,
            converged_cfg(),
            6000,
        );
        let t_abs = torque.abs();
        let reference = analytic_couette_torque(omega, R_I_B1, R_O_EFF_B1);
        println!(
            "ACC IBM B1 heavy: omega={omega:.6e}, steps={steps}, measured_T={torque:.6e}, |T|={t_abs:.6e}, T_analytic={reference:.6e}, ratio={:.6e}",
            t_abs / reference
        );
        samples.push((omega, t_abs));
        measured.push(t_abs);
    }
    let fit = linear_fit(&omegas, &measured);
    let mid_t = measured[1].abs();
    let agreement = curve_agreement(
        |omega| analytic_couette_torque(omega, R_I_B1, R_O_EFF_B1),
        &samples,
        0.25,
        0.0,
    );
    println!(
        "ACC IBM B1 heavy: slope={:.6e}, intercept={:.6e}, r2={:.6e}, max_rel_dev={:.6e}, worst_omega={:.6e}, frac_in_band={:.6e}",
        fit.slope, fit.intercept, fit.r2, agreement.max_rel_dev, agreement.worst_x, agreement.frac_in_band
    );
    assert!(
        fit.r2 >= 0.999,
        "ACC IBM B1 heavy measured r2={:.6e} below band 9.99e-1 for linear_fit(T,Omega)",
        fit.r2
    );
    assert!(
        fit.intercept.abs() <= 0.02 * mid_t,
        "ACC IBM B1 heavy measured |intercept|={:.6e} exceeds band 2e-2*T(mid)={:.6e}",
        fit.intercept.abs(),
        0.02 * mid_t
    );
    assert!(
        agreement.max_rel_dev <= 0.25,
        "ACC IBM B1 heavy measured max_rel_dev={:.6e} exceeds band 2.5e-1 using |T-T_analytic|/T_analytic at omega={:.6e}",
        agreement.max_rel_dev,
        agreement.worst_x
    );
}

#[test]
fn b2_sub_cell_translation_invariance_of_steady_torque() {
    let centers = [[40.0, 40.0], [40.3, 40.17], [40.5, 40.5]];
    let mut torques = Vec::new();
    for center in centers {
        let solid = b1_outer_wall_solid(NX_B1, NY_B1, center, R_O_B1);
        let (torque, diag, steps) = run_torque_to_steady(
            NX_B1,
            NY_B1,
            center,
            R_I_B1,
            OMEGA_MID,
            N_MARKERS_B1,
            solid,
            converged_cfg(),
            6000,
        );
        println!(
            "ACC IBM B2: center=({:.3},{:.3}), steps={steps}, measured_T={torque:.6e}, |T|={:.6e}, slip_max_rel={:.6e}",
            center[0],
            center[1],
            torque.abs(),
            diag.slip_max_rel
        );
        torques.push(torque.abs());
    }
    let min_t = torques.iter().copied().fold(f64::INFINITY, f64::min);
    let max_t = torques.iter().copied().fold(0.0f64, f64::max);
    let mean = torques.iter().sum::<f64>() / torques.len() as f64;
    let spread = (max_t - min_t) / mean;
    println!(
        "ACC IBM B2: min_T={min_t:.6e}, max_T={max_t:.6e}, mean_T={mean:.6e}, spread_over_mean={spread:.6e}"
    );
    assert!(
        spread <= 0.05,
        "ACC IBM B2 measured spread={spread:.6e} exceeds band 5e-2 using (max-min)/mean steady |torque|"
    );
}

#[test]
fn b3_force_spreading_conservation_domain_edge_kernel_truncation_pin() {
    // Spreading and interpolation share the same kernel, so sum_cells F_cell =
    // sum_markers F_marker is a quadrature identity only when each marker's
    // full kernel support lies inside the domain. If a marker sits within the
    // kernel radius of a non-periodic edge, out-of-domain stencil points are
    // dropped and the physical, pre-truncation marker force need not equal the
    // spread fluid force. SPEC-GAP candidate: near-wall marker support
    // truncation behavior is not documented by the public API.
    let (nx, ny) = (48, 48);
    let mut walls = WallSpec::default();
    walls.is_wall[Face::XNeg.index()] = true;
    walls.is_wall[Face::XPos.index()] = true;
    let spec = GlobalSpec {
        dims: [nx, ny, 1],
        nu: NU,
        periodic: [false, true, false],
        collision: CollisionKind::Trt {
            magic: CollisionKind::MAGIC_STD,
        },
        ..Default::default()
    };
    let (solid, wall_u) = build_wall_rims(2, spec.dims, &walls);
    let mut solver = Solver::<D2Q9, f64, CpuScalar, LocalPeriodic>::new(
        &spec,
        &solid,
        &wall_u,
        [1, 1, 1],
        CpuScalar::default(),
        LocalPeriodic,
    );
    let body = RotatingBody::circle_2d([10.0, 24.0], 8.0, OMEGA_MID, 128);
    solver.clear_body_force_field();
    let diag = solver.apply_rotating_ibm(&body, converged_cfg());
    println!(
        "ACC IBM B3: momentum_error_rel={:.6e}, fluid_force=({:.6e},{:.6e}), marker_force=({:.6e},{:.6e}), slip_max_rel={:.6e}",
        diag.momentum_error_rel,
        diag.fluid_force[0],
        diag.fluid_force[1],
        diag.marker_force[0],
        diag.marker_force[1],
        diag.slip_max_rel
    );
    if diag.momentum_error_rel <= 1.0e-12 {
        assert!(
            diag.momentum_error_rel <= 1.0e-12,
            "ACC IBM B3 conservative surprise branch failed: measured momentum_error_rel={:.6e} should be <= 1e-12",
            diag.momentum_error_rel
        );
    } else {
        assert!(
            diag.momentum_error_rel > 1.0e-10 && diag.momentum_error_rel < 0.5,
            "ACC IBM B3 measured momentum_error_rel={:.6e} outside sanity band (1e-10,5e-1) for near-edge kernel truncation",
            diag.momentum_error_rel
        );
    }
}

fn slip_after_steady_resolution(n: usize) -> (f64, IbmDiagnostics, usize) {
    let center = [0.5 * n as f64, 0.5 * n as f64];
    let r_i = n as f64 / 8.0;
    let r_o = 3.0 * n as f64 / 8.0;
    let u_tip = 0.02 * (48.0 / n as f64);
    let omega = u_tip / r_i;
    let solid = b1_outer_wall_solid(n, n, center, r_o);
    let (_torque, diag, steps) = run_torque_to_steady(
        n,
        n,
        center,
        r_i,
        omega,
        (16.0 * r_i).round() as usize,
        solid,
        converged_cfg(),
        4000 * (n / 48),
    );
    (diag.slip_max_rel, diag, steps)
}

#[test]
fn b4_marker_slip_decreases_under_resolution_refinement_canary() {
    let mut errs = Vec::new();
    for n in [48usize, 96] {
        let (err, diag, steps) = slip_after_steady_resolution(n);
        println!(
            "ACC IBM B4: N={n}, steps={steps}, slip_max_rel={err:.6e}, slip_rms_rel={:.6e}, iterations={}",
            diag.slip_rms_rel, diag.iterations
        );
        errs.push(err);
    }
    let mono = monotonicity(&errs);
    println!(
        "ACC IBM B4: monotonicity={mono:.6e}, err48={:.6e}, err96={:.6e}",
        errs[0], errs[1]
    );
    assert!(
        mono == 1.0,
        "ACC IBM B4 measured monotonicity={mono:.6e} below band 1.0 for slip_max_rel under N=48->96 refinement"
    );
}

#[test]
#[ignore = "heavy ACC-AUDIT B4 three-resolution slip convergence order"]
fn b4_marker_slip_convergence_order_heavy() {
    let ns = [48usize, 96, 192];
    let mut h = Vec::new();
    let mut err = Vec::new();
    for n in ns {
        let (e, diag, steps) = slip_after_steady_resolution(n);
        println!(
            "ACC IBM B4 heavy: N={n}, h={:.6e}, steps={steps}, slip_max_rel={e:.6e}, slip_rms_rel={:.6e}",
            1.0 / n as f64,
            diag.slip_rms_rel
        );
        h.push(1.0 / n as f64);
        err.push(e);
    }
    let fit = order_fit(&h, &err);
    println!(
        "ACC IBM B4 heavy: order_slope={:.6e}, intercept={:.6e}, r2={:.6e}",
        fit.slope, fit.intercept, fit.r2
    );
    assert!(
        (0.8..=2.5).contains(&fit.slope),
        "ACC IBM B4 heavy measured order={:.6e} outside band [8e-1,2.5] from order_fit(h=1/N, slip_max_rel)",
        fit.slope
    );
    assert!(
        fit.r2 >= 0.98,
        "ACC IBM B4 heavy measured r2={:.6e} below band 9.8e-1 for order_fit(h=1/N, slip_max_rel)",
        fit.r2
    );
}

#[test]
fn b5_kernel_relaxation_cross_path_consistency_of_steady_torque() {
    let center = [40.0, 40.0];
    let configs = [
        DirectForcingConfig {
            max_iterations: 4,
            slip_tolerance: 1.0e-3,
            kernel_radius: 1,
            relaxation: 1.0,
        },
        DirectForcingConfig {
            max_iterations: 4,
            slip_tolerance: 1.0e-3,
            kernel_radius: 2,
            relaxation: 1.0,
        },
        DirectForcingConfig {
            max_iterations: 6,
            slip_tolerance: 1.0e-3,
            kernel_radius: 2,
            relaxation: 0.7,
        },
    ];
    let mut torques = Vec::new();
    for cfg in configs {
        let solid = b1_outer_wall_solid(NX_B1, NY_B1, center, R_O_B1);
        let (torque, diag, steps) = run_torque_to_steady(
            NX_B1,
            NY_B1,
            center,
            R_I_B1,
            OMEGA_MID,
            N_MARKERS_B1,
            solid,
            cfg,
            6000,
        );
        let t_abs = torque.abs();
        println!(
            "ACC IBM B5: kernel_radius={}, relaxation={:.3}, max_iterations={}, steps={steps}, measured_T={torque:.6e}, |T|={t_abs:.6e}, slip_max_rel={:.6e}",
            cfg.kernel_radius, cfg.relaxation, cfg.max_iterations, diag.slip_max_rel
        );
        torques.push(t_abs);
    }
    for i in 0..torques.len() {
        for j in i + 1..torques.len() {
            let denom = 0.5 * (torques[i].abs() + torques[j].abs());
            let rel = (torques[i] - torques[j]).abs() / denom;
            assert!(
                rel <= 0.05,
                "ACC IBM B5 measured pairwise rel_diff={rel:.6e} exceeds band 5e-2 using |Ti-Tj|/mean(|Ti|,|Tj|), i={i}, j={j}, Ti={:.6e}, Tj={:.6e}",
                torques[i],
                torques[j]
            );
        }
    }
}

fn run_b6_transient(omega: f64) -> (Vec<f64>, Vec<f64>, f64) {
    let center = [40.0, 40.0];
    let solid = b1_outer_wall_solid(NX_B1, NY_B1, center, R_O_B1);
    let mut solver = solver_from_solid(NX_B1, NY_B1, NU, solid);
    let body = RotatingBody::circle_2d(center, R_I_B1, omega, N_MARKERS_B1);
    let mut torque = 0.0;
    for _ in 0..200 {
        solver.clear_body_force_field();
        let diag = solver.apply_rotating_ibm(&body, converged_cfg());
        torque = diag.torque[2];
        solver.step();
    }
    (solver.gather_ux(), solver.gather_uy(), torque)
}

#[test]
fn b6_rotation_antisymmetry_transient_torque_and_field() {
    let (ux_p, uy_p, tq_p) = run_b6_transient(OMEGA_MID);
    let (ux_m, uy_m, tq_m) = run_b6_transient(-OMEGA_MID);
    let torque_diff = (tq_p + tq_m).abs();
    let torque_scale = tq_p.abs().max(tq_m.abs()).max(1.0e-30);

    // Reflection x -> 2 xc - x maps a -Omega solution to the +Omega solution.
    // The reflected vector basis changes the normal component sign and keeps
    // the tangential y component: (ux, uy)_plus(x,y) =
    // (-ux_minus(mirror_x,y), uy_minus(mirror_x,y)).
    let center = [40.0, 40.0];
    let solid = b1_outer_wall_solid(NX_B1, NY_B1, center, R_O_B1);
    let mirror_twice = (2.0 * center[0]).round() as isize;
    let mut max_field_diff = 0.0f64;
    for y in 0..NY_B1 {
        for x in 0..NX_B1 {
            let i = y * NX_B1 + x;
            if solid[i] {
                continue;
            }
            let xm = (mirror_twice - x as isize).rem_euclid(NX_B1 as isize) as usize;
            let im = y * NX_B1 + xm;
            max_field_diff = max_field_diff.max((ux_p[i] + ux_m[im]).abs());
            max_field_diff = max_field_diff.max((uy_p[i] - uy_m[im]).abs());
        }
    }
    println!(
        "ACC IBM B6: torque_plus={tq_p:.6e}, torque_minus={tq_m:.6e}, |sum|={torque_diff:.6e}, torque_scale={torque_scale:.6e}, max_mapped_field_diff={max_field_diff:.6e}"
    );
    assert!(
        torque_diff <= 1.0e-12 * torque_scale,
        "ACC IBM B6 measured torque antisym diff={torque_diff:.6e} exceeds band 1e-12*|T|={:.6e} using |T(+Omega)+T(-Omega)|",
        1.0e-12 * torque_scale
    );
    assert!(
        max_field_diff <= 1.0e-12,
        "ACC IBM B6 measured mapped velocity max_abs_diff={max_field_diff:.6e} exceeds band 1e-12 over fluid cells"
    );
}

#[test]
#[ignore = "heavy ACC-AUDIT B7 spin-up transient Stokes-layer torque"]
fn b7_spin_up_transient_torque_matches_stokes_layer_asymptote() {
    // Impulsive cylinder rotation has an early-time viscous layer
    // delta=sqrt(nu t) thin compared with r_i. Locally the surface is the
    // Stokes-I plane-wall problem, tau_w(t)=mu Omega r_i/sqrt(pi nu t).
    // Multiplying by the circumference moment arm gives
    // T(t) ~= 2 pi r_i^2 tau_w
    //      = 2 pi mu Omega r_i^3 / sqrt(pi nu t), so T(t)*sqrt(t) is constant.
    // The window begins after delta>=2 cells: t>=4/nu=24. It ends at
    // delta<=r_i/3: t<=r_i^2/(9 nu)=170 for r_i=16, nu=1/6. ANOM-P2-001
    // (per-cell force-field one-step impulse deficit, owned by
    // accuracy_audit.rs) pollutes only the first few steps; this window avoids
    // re-pinning it here.
    let (nx, ny) = (128, 128);
    let center = [64.0, 64.0];
    let r_i = 16.0;
    let omega = 0.01 / r_i;
    let mut solver = solver_from_solid(nx, ny, NU, vec![false; nx * ny]);
    let body = RotatingBody::circle_2d(center, r_i, omega, 256);
    let mut samples = Vec::new();
    let mut torques = Vec::new();
    for t in 1..=170 {
        solver.clear_body_force_field();
        let diag = solver.apply_rotating_ibm(&body, converged_cfg());
        solver.step();
        if t >= 24 {
            samples.push((t as f64, diag.torque[2].abs()));
            torques.push(diag.torque[2].abs());
        }
    }
    let mu = NU;
    let theory = |t: f64| {
        2.0 * std::f64::consts::PI * mu * omega * r_i * r_i * r_i
            / (std::f64::consts::PI * NU * t).sqrt()
    };
    let agreement = curve_agreement(theory, &samples, 0.25, 0.0);
    let mono = monotonicity(&torques);
    println!(
        "ACC IBM B7: samples={}, first_T={:.6e}, last_T={:.6e}, max_rel_dev={:.6e}, worst_t={:.6e}, monotonicity={mono:.6e}",
        samples.len(),
        torques[0],
        torques[torques.len() - 1],
        agreement.max_rel_dev,
        agreement.worst_x
    );
    assert!(
        agreement.max_rel_dev <= 0.25,
        "ACC IBM B7 measured max_rel_dev={:.6e} exceeds band 2.5e-1 using |T-T_stokes|/T_stokes at t={:.6e}",
        agreement.max_rel_dev,
        agreement.worst_x
    );
    assert!(
        mono == 1.0,
        "ACC IBM B7 measured monotonicity={mono:.6e} below band 1.0 for monotone decay of |T| in t=24..170"
    );
}

#[test]
fn b8_converged_taylor_couette_profile_accuracy() {
    let (nx, ny) = (80, 80);
    let center = [40.0, 40.0];
    let r_i = 10.0;
    let r_o = 30.0;
    let omega = OMEGA_MID;
    let body = RotatingBody::circle_2d(center, r_i, omega, 160);
    // IBM supplies the rotating inner boundary. Adding a stationary Eulerian
    // solid core inside the marker radius creates a separate half-way
    // bounce-back wall inside the IBM kernel support and contaminates this
    // profile probe with a narrow stationary-wall gap artifact.
    let solid = b1_outer_wall_solid(nx, ny, center, r_o);
    let mut solver = solver_from_solid(nx, ny, NU, solid);
    let mut last = IbmDiagnostics::default();
    for _ in 0..1500 {
        solver.clear_body_force_field();
        last = solver.apply_rotating_ibm(&body, converged_cfg());
        solver.step();
    }
    let ux = solver.gather_ux();
    let uy = solver.gather_uy();
    // The same Stokes annular-Couette derivation as B1 gives
    // u_theta(r)=a r + b/r, a=-Omega r_i^2/(r_o_eff^2-r_i^2),
    // b=Omega r_i^2 r_o_eff^2/(r_o_eff^2-r_i^2), with stationary outer wall
    // at r_o_eff=r_o+1.0 and inner no-slip speed U_i=Omega r_i.
    let r_o_eff = r_o + 1.0;
    let a = -omega * r_i * r_i / (r_o_eff * r_o_eff - r_i * r_i);
    let b = omega * r_i * r_i * r_o_eff * r_o_eff / (r_o_eff * r_o_eff - r_i * r_i);
    let mut actual = Vec::new();
    let mut reference = Vec::new();
    for y in 0..ny {
        for x in 0..nx {
            let dx = x as f64 - center[0];
            let dy = y as f64 - center[1];
            let r = (dx * dx + dy * dy).sqrt();
            if !(r_i + 4.0..=r_o - 5.0).contains(&r) {
                continue;
            }
            actual.push((-ux[y * nx + x] * dy + uy[y * nx + x] * dx) / r);
            reference.push(a * r + b / r);
        }
    }
    let l2 = l2_rel(&actual, &reference);
    let linf = linf_rel(&actual, &reference, omega * r_i);
    println!(
        "ACC IBM B8: samples={}, slip_max_rel={:.6e}, slip_rms_rel={:.6e}, torque_z={:.6e}, l2_rel={l2:.6e}, linf_rel_floor_Ui={linf:.6e}",
        actual.len(),
        last.slip_max_rel,
        last.slip_rms_rel,
        last.torque[2]
    );
    assert!(
        l2 <= 0.10,
        "ACC IBM B8 measured l2_rel={l2:.6e} exceeds band 1e-1 using ||u_theta-u_ref||2/||u_ref||2"
    );
    assert!(
        linf <= 0.25,
        "ACC IBM B8 measured linf_rel={linf:.6e} exceeds band 2.5e-1 using max|u_theta-u_ref|/max(max|u_ref|,U_i)"
    );
}
