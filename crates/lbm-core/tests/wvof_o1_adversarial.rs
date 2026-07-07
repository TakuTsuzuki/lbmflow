//! Independent adversarial W-VOF O1 phase-field transport tests.
//!
//! These tests are written from `docs/proposals/WVOF_IMPL_SPEC.md` and
//! `docs/REQ_STIRRED_REACTOR.md`, using only the public O1 phase-field API.
//! Bands are derived here rather than fitted to implementation output.

use lbm_core::prelude::*;
use std::f64::consts::PI;
use std::path::Path;

type Sim = Solver<D3Q19, f64, CpuScalar, LocalPeriodic>;

const W_TARGET: f64 = 4.0;
const M_DEFAULT: f64 = 0.05;
const PHI_DRIFT_BAND: f64 = 1.0e-3;
const WIDTH_LO: f64 = 4.0;
const WIDTH_HI: f64 = 5.0;
const ROTATION_LINF_BAND: f64 = 1.0e-10;

fn params(mobility: f64) -> PhaseFieldParams<f64> {
    PhaseFieldParams::new(W_TARGET, mobility)
}

fn periodic_solver(dims: [usize; 3]) -> Sim {
    let spec = GlobalSpec::<f64> {
        dims,
        nu: 0.04,
        periodic: [true, true, true],
        ..Default::default()
    };
    Solver::new(
        &spec,
        &[],
        &[],
        [1, 1, 1],
        CpuScalar::default(),
        LocalPeriodic,
    )
}

fn idx(dims: [usize; 3], x: usize, y: usize, z: usize) -> usize {
    (z * dims[1] + y) * dims[0] + x
}

fn periodic_delta(a: f64, b: f64, n: usize) -> f64 {
    let mut d = a - b;
    let nf = n as f64;
    if d > 0.5 * nf {
        d -= nf;
    } else if d < -0.5 * nf {
        d += nf;
    }
    d
}

fn tanh_sphere_phi(dims: [usize; 3], center: [f64; 3], radius: f64, width: f64) -> Vec<f64> {
    let mut phi = vec![0.0; dims[0] * dims[1] * dims[2]];
    for z in 0..dims[2] {
        for y in 0..dims[1] {
            for x in 0..dims[0] {
                let dx = periodic_delta(x as f64, center[0], dims[0]);
                let dy = periodic_delta(y as f64, center[1], dims[1]);
                let dz = periodic_delta(z as f64, center[2], dims[2]);
                let r = (dx * dx + dy * dy + dz * dz).sqrt();
                phi[idx(dims, x, y, z)] = 0.5 * (1.0 - ((2.0 * (r - radius)) / width).tanh());
            }
        }
    }
    phi
}

fn two_droplet_phi(dims: [usize; 3], width: f64) -> Vec<f64> {
    let a = tanh_sphere_phi(dims, [20.0, 16.0, 16.0], 8.0, width);
    let b = tanh_sphere_phi(dims, [44.0, 16.0, 16.0], 8.0, width);
    let mut out = vec![0.0; a.len()];
    for i in 0..out.len() {
        out[i] = if a[i] > b[i] { a[i] } else { b[i] };
    }
    out
}

fn enable(mut sim: Sim, phase: &[f64], velocity: [f64; 3], p: PhaseFieldParams<f64>) -> Sim {
    sim.enable_phase_field_prescribed_velocity(p, phase, move |_, _, _| velocity)
        .expect("phase field should enable");
    sim
}

fn run_phase(
    mut sim: Sim,
    steps: usize,
    velocity: [f64; 3],
    p: PhaseFieldParams<f64>,
) -> (Sim, PhaseFieldDiagnostics) {
    let mut diag = PhaseFieldDiagnostics::default();
    for _ in 0..steps {
        diag = sim
            .phase_field_step_prescribed_velocity(p, move |_, _, _| velocity)
            .expect("phase field step should succeed");
    }
    (sim, diag)
}

fn total_phi(phi: &[f64]) -> f64 {
    phi.iter().sum()
}

fn circular_center(phi: &[f64], dims: [usize; 3]) -> [f64; 3] {
    let mut center = [0.0; 3];
    for axis in 0..3 {
        let n = dims[axis] as f64;
        let mut s = 0.0;
        let mut c = 0.0;
        for z in 0..dims[2] {
            for y in 0..dims[1] {
                for x in 0..dims[0] {
                    let p = [x, y, z];
                    let angle = 2.0 * PI * p[axis] as f64 / n;
                    let w = phi[idx(dims, x, y, z)];
                    s += w * angle.sin();
                    c += w * angle.cos();
                }
            }
        }
        let mut angle = s.atan2(c);
        if angle < 0.0 {
            angle += 2.0 * PI;
        }
        center[axis] = angle * n / (2.0 * PI);
    }
    center
}

fn fit_tanh_width(phi: &[f64], dims: [usize; 3]) -> f64 {
    let center = circular_center(phi, dims);
    let mut n = 0.0;
    let mut sr = 0.0;
    let mut sy = 0.0;
    let mut srr = 0.0;
    let mut sry = 0.0;
    for z in 0..dims[2] {
        for y in 0..dims[1] {
            for x in 0..dims[0] {
                let p = phi[idx(dims, x, y, z)];
                if p <= 0.1 || p >= 0.9 {
                    continue;
                }
                let dx = periodic_delta(x as f64, center[0], dims[0]);
                let dy = periodic_delta(y as f64, center[1], dims[1]);
                let dz = periodic_delta(z as f64, center[2], dims[2]);
                let r = (dx * dx + dy * dy + dz * dz).sqrt();
                let yfit = (1.0 - 2.0 * p).atanh();
                n += 1.0;
                sr += r;
                sy += yfit;
                srr += r * r;
                sry += r * yfit;
            }
        }
    }
    let denom = n * srr - sr * sr;
    assert!(
        n > 20.0 && denom.abs() > 1.0e-12,
        "not enough interface samples for tanh fit"
    );
    let slope = (n * sry - sr * sy) / denom;
    2.0 / slope
}

fn l2_rel(a: &[f64], b: &[f64]) -> f64 {
    let mut num = 0.0;
    let mut den = 0.0;
    for i in 0..a.len() {
        let d = a[i] - b[i];
        num += d * d;
        den += b[i] * b[i];
    }
    (num / den).sqrt()
}

fn max_abs_delta(a: &[f64], b: &[f64]) -> f64 {
    let mut m = 0.0;
    for i in 0..a.len() {
        let d = (a[i] - b[i]).abs();
        if d > m {
            m = d;
        }
    }
    m
}

fn write_midplane_pgm(path: impl AsRef<Path>, phi: &[f64], dims: [usize; 3]) {
    let path = path.as_ref();
    std::fs::create_dir_all(path.parent().expect("artifact path should have parent"))
        .expect("artifact directory should be creatable");
    let z = dims[2] / 2;
    let mut bytes = format!("P2\n{} {}\n255\n", dims[0], dims[1]);
    for y in 0..dims[1] {
        for x in 0..dims[0] {
            let p = phi[idx(dims, x, y, z)];
            let scaled = if p <= 0.0 {
                0
            } else if p >= 1.0 {
                255
            } else {
                (p * 255.0).round() as i32
            };
            bytes.push_str(&format!("{scaled} "));
        }
        bytes.push('\n');
    }
    std::fs::write(path, bytes).expect("artifact should be writable");
}

fn write_width_csv(path: impl AsRef<Path>, rows: &[(usize, f64, f64)]) {
    let path = path.as_ref();
    std::fs::create_dir_all(path.parent().expect("artifact path should have parent"))
        .expect("artifact directory should be creatable");
    let mut csv = String::from("step,width,target_error\n");
    for (step, width, target_error) in rows {
        csv.push_str(&format!("{step},{width:.12e},{target_error:.12e}\n"));
    }
    std::fs::write(path, csv).expect("width CSV artifact should be writable");
}

fn assert_width_band(label: &str, width: f64) {
    assert!(
        (WIDTH_LO..=WIDTH_HI).contains(&width),
        "{label}: fitted tanh width {width:.6e} outside spec validity band [{WIDTH_LO}, {WIDTH_HI}]"
    );
}

#[test]
fn interface_width_survives_multiple_box_transits_at_low_and_high_mach() {
    let dims = [32, 24, 24];
    let initial = tanh_sphere_phi(dims, [16.0, 12.0, 12.0], 8.0, W_TARGET);
    let p = params(M_DEFAULT);
    let mut failures = Vec::new();
    for (speed, steps) in [(0.02, 4_800usize), (0.08, 1_200usize)] {
        let initial_width = fit_tanh_width(&initial, dims);
        let sim = enable(periodic_solver(dims), &initial, [speed, 0.0, 0.0], p);
        let (mut sim, diag) = run_phase(sim, steps, [speed, 0.0, 0.0], p);
        let final_phi = sim.gather_phi().expect("phi should gather");
        let final_width = fit_tanh_width(&final_phi, dims);
        let drift = ((total_phi(&final_phi) - total_phi(&initial)) / total_phi(&initial)).abs();
        write_midplane_pgm(
            format!("/tmp/lbmflow_wvof_o1_adversarial/width_u_{speed:.2}.pgm"),
            &final_phi,
            dims,
        );
        println!(
            "WVOF-O1 width transit speed={speed:.2} steps={steps} initial_width={initial_width:.9e} final_width={final_width:.9e} phi_drift={drift:.9e} min_phi={:.9e} max_phi={:.9e}",
            diag.min_phi, diag.max_phi
        );
        assert_width_band("initial droplet", initial_width);
        if !(WIDTH_LO..=WIDTH_HI).contains(&final_width) {
            failures.push(format!(
                "speed={speed:.2}: fitted tanh width {final_width:.6e} outside [{WIDTH_LO}, {WIDTH_HI}]"
            ));
        }
        if drift > PHI_DRIFT_BAND {
            failures.push(format!(
                "speed={speed:.2}: phi drift {drift:.6e} exceeds {PHI_DRIFT_BAND:.6e}"
            ));
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("; "));
}

#[test]
fn mobility_domain_edges_conserve_phi_and_validate_rejections() {
    let dims = [32, 24, 24];
    let initial = tanh_sphere_phi(dims, [16.0, 12.0, 12.0], 8.0, W_TARGET);
    for m in [1.0e-4, 1.0 / 6.0] {
        let p = params(m).validate().expect("mobility edge should be valid");
        let sim = enable(periodic_solver(dims), &initial, [0.04, 0.0, 0.0], p);
        let (mut sim, diag) = run_phase(sim, 800, [0.04, 0.0, 0.0], p);
        let final_phi = sim.gather_phi().expect("phi should gather");
        let width = fit_tanh_width(&final_phi, dims);
        let drift = ((total_phi(&final_phi) - total_phi(&initial)) / total_phi(&initial)).abs();
        write_midplane_pgm(
            format!("/tmp/lbmflow_wvof_o1_adversarial/mobility_{m:.6}.pgm"),
            &final_phi,
            dims,
        );
        println!(
            "WVOF-O1 mobility edge M={m:.9e} steps=800 width={width:.9e} phi_drift={drift:.9e} min_phi={:.9e} max_phi={:.9e}",
            diag.min_phi, diag.max_phi
        );
        assert_width_band("mobility edge droplet", width);
        assert!(
            drift <= PHI_DRIFT_BAND,
            "M={m:.9e}: phi drift {drift:.6e} exceeds {PHI_DRIFT_BAND:.6e}"
        );
    }

    for m in [-1.0e-6, 0.0, (1.0 / 6.0) + 1.0e-6] {
        assert!(
            params(m).validate().is_err(),
            "out-of-domain mobility M={m:.9e} must be rejected by spec domain (0, 1/6]"
        );
    }
}

#[test]
fn two_droplet_zalesak_style_one_period_conserves_and_returns_shape() {
    let dims = [64, 32, 32];
    let speed = 0.08;
    let steps = 800usize;
    let initial = two_droplet_phi(dims, W_TARGET);
    let p = params(M_DEFAULT);
    let sim = enable(periodic_solver(dims), &initial, [speed, 0.0, 0.0], p);
    let (mut sim, diag) = run_phase(sim, steps, [speed, 0.0, 0.0], p);
    let final_phi = sim.gather_phi().expect("phi should gather");
    let drift = ((total_phi(&final_phi) - total_phi(&initial)) / total_phi(&initial)).abs();
    let shape = l2_rel(&final_phi, &initial);
    let gap_phi = final_phi[idx(dims, 32, 16, 16)];
    let left_core = final_phi[idx(dims, 20, 16, 16)];
    let right_core = final_phi[idx(dims, 44, 16, 16)];
    write_midplane_pgm(
        "/tmp/lbmflow_wvof_o1_adversarial/zalesak_two_droplet_period.pgm",
        &final_phi,
        dims,
    );
    println!(
        "WVOF-O1 two-droplet one-period steps={steps} l2_rel={shape:.9e} phi_drift={drift:.9e} gap_phi={gap_phi:.9e} left_core={left_core:.9e} right_core={right_core:.9e} min_phi={:.9e} max_phi={:.9e}",
        diag.min_phi, diag.max_phi
    );
    assert!(
        drift <= PHI_DRIFT_BAND,
        "two-droplet phi drift {drift:.6e} exceeds {PHI_DRIFT_BAND:.6e}"
    );
    // Conservative Allen-Cahn O1 uses a second-order lattice stencil. For
    // W=4, the expected normalized shape error scale after a reversible full
    // period is O((dx/W)^2)=1/16. A factor of two covers finite interface
    // curvature and the two disconnected interfaces without fitting output.
    let band = 2.0 / (W_TARGET * W_TARGET);
    assert!(
        shape <= band,
        "two-droplet one-period L2rel {shape:.6e} exceeds second-order band {band:.6e}"
    );
    assert!(
        gap_phi < 0.25 && left_core > 0.75 && right_core > 0.75,
        "two-droplet behavior anchor failed: gap_phi={gap_phi:.6e}, left_core={left_core:.6e}, right_core={right_core:.6e}"
    );
}

#[test]
fn rotation_metamorphic_phi_fields_are_index_permuted() {
    let dims = [32, 32, 32];
    let speed = 0.04;
    let steps = 240usize;
    let initial = tanh_sphere_phi(dims, [16.0, 16.0, 16.0], 8.0, W_TARGET);
    let p = params(M_DEFAULT);
    let x_sim = enable(periodic_solver(dims), &initial, [speed, 0.0, 0.0], p);
    let y_sim = enable(periodic_solver(dims), &initial, [0.0, speed, 0.0], p);
    let z_sim = enable(periodic_solver(dims), &initial, [0.0, 0.0, speed], p);
    let (mut x_sim, _) = run_phase(x_sim, steps, [speed, 0.0, 0.0], p);
    let (mut y_sim, _) = run_phase(y_sim, steps, [0.0, speed, 0.0], p);
    let (mut z_sim, _) = run_phase(z_sim, steps, [0.0, 0.0, speed], p);
    let px = x_sim.gather_phi().expect("x phi should gather");
    let py = y_sim.gather_phi().expect("y phi should gather");
    let pz = z_sim.gather_phi().expect("z phi should gather");

    let mut max_delta = 0.0;
    for z in 0..dims[2] {
        for y in 0..dims[1] {
            for x in 0..dims[0] {
                let base = px[idx(dims, x, y, z)];
                let y_map = py[idx(dims, y, x, z)];
                let z_map = pz[idx(dims, z, y, x)];
                let dy = (base - y_map).abs();
                let dz = (base - z_map).abs();
                if dy > max_delta {
                    max_delta = dy;
                }
                if dz > max_delta {
                    max_delta = dz;
                }
            }
        }
    }
    write_midplane_pgm(
        "/tmp/lbmflow_wvof_o1_adversarial/rotation_x_reference.pgm",
        &px,
        dims,
    );
    println!(
        "WVOF-O1 rotation metamorphic steps={steps} max_linf_delta={max_delta:.9e} band={ROTATION_LINF_BAND:.9e}"
    );
    // Tolerance basis follows the existing orientation convention:
    // validation_cavity.rs and d3q27_open_metamorphic.rs use L_inf <= 1e-10
    // for exact coordinate equivariance after finite time.
    assert!(
        max_delta <= ROTATION_LINF_BAND,
        "phase-field rotation equivariance max |delta|={max_delta:.6e} exceeds {ROTATION_LINF_BAND:.6e}"
    );
}

#[test]
fn anti_diffusion_counter_term_sharpens_a_slightly_diffused_interface() {
    let dims = [40, 32, 32];
    let diffused_width = 5.0;
    let initial = tanh_sphere_phi(dims, [20.0, 16.0, 16.0], 10.0, diffused_width);
    let p = params(M_DEFAULT);
    let initial_width = fit_tanh_width(&initial, dims);
    let initial_target_error = (initial_width - W_TARGET).abs();
    let mut sim = enable(periodic_solver(dims), &initial, [0.0, 0.0, 0.0], p);
    let mut rows = vec![(0usize, initial_width, initial_target_error)];
    let mut diag = PhaseFieldDiagnostics::default();
    for step in 1..=500 {
        diag = sim
            .phase_field_step_prescribed_velocity(p, |_, _, _| [0.0, 0.0, 0.0])
            .expect("phase field step should succeed");
        if step % 50 == 0 || step == 500 {
            let phi = sim.gather_phi().expect("phi should gather");
            let width = fit_tanh_width(&phi, dims);
            rows.push((step, width, (width - W_TARGET).abs()));
        }
    }
    let final_phi = sim.gather_phi().expect("phi should gather");
    let final_width = fit_tanh_width(&final_phi, dims);
    let final_target_error = (final_width - W_TARGET).abs();
    write_width_csv(
        Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../target/wvof_o1/counter_term_width_vs_time.csv"),
        &rows,
    );
    write_midplane_pgm(
        "/tmp/lbmflow_wvof_o1_adversarial/counter_term_sharpening.pgm",
        &final_phi,
        dims,
    );
    println!(
        "WVOF-O1 counter-term sign initial_width={initial_width:.9e} final_width={final_width:.9e} initial_target_error={initial_target_error:.9e} final_target_error={final_target_error:.9e} min_phi={:.9e} max_phi={:.9e}",
        diag.min_phi, diag.max_phi
    );
    assert!(
        final_width < initial_width && final_target_error < initial_target_error,
        "counter-term sign anchor failed: diffused interface should sharpen toward W={W_TARGET}, initial_width={initial_width:.6e}, final_width={final_width:.6e}"
    );
}

#[test]
fn g_none_canonical_hydrodynamic_run_is_bit_identical() {
    fn solver() -> Sim {
        let dims = [18, 16, 14];
        let spec = GlobalSpec::<f64> {
            dims,
            nu: 0.04,
            periodic: [true, true, true],
            force: [2.0e-7, -1.0e-7, 3.0e-7],
            ..Default::default()
        };
        let mut s = Solver::new(
            &spec,
            &[],
            &[],
            [1, 1, 1],
            CpuScalar::default(),
            LocalPeriodic,
        );
        s.init_with(|x, y, z| {
            let kx = 2.0 * PI * x as f64 / dims[0] as f64;
            let ky = 2.0 * PI * y as f64 / dims[1] as f64;
            let kz = 2.0 * PI * z as f64 / dims[2] as f64;
            (
                1.0 + 1.0e-4 * (kx.cos() + ky.sin() * kz.cos()),
                [0.01 * ky.sin(), 0.008 * kz.cos(), -0.006 * kx.sin()],
            )
        });
        s
    }

    let mut reference = solver();
    let mut candidate = solver();
    assert!(
        candidate.gather_phi().is_err(),
        "g=None regression precondition failed: phi must not be enabled"
    );
    reference.run(120);
    candidate.run(120);

    let rho_ref = reference.gather_rho();
    let rho_candidate = candidate.gather_rho();
    let ux_ref = reference.gather_ux();
    let ux_candidate = candidate.gather_ux();
    let uy_ref = reference.gather_uy();
    let uy_candidate = candidate.gather_uy();
    let uz_ref = reference.gather_uz();
    let uz_candidate = candidate.gather_uz();
    assert_eq!(
        rho_candidate, rho_ref,
        "rho field must be bit-identical with g=None"
    );
    assert_eq!(
        ux_candidate, ux_ref,
        "ux field must be bit-identical with g=None"
    );
    assert_eq!(
        uy_candidate, uy_ref,
        "uy field must be bit-identical with g=None"
    );
    assert_eq!(
        uz_candidate, uz_ref,
        "uz field must be bit-identical with g=None"
    );
    let d_rho = max_abs_delta(&rho_candidate, &rho_ref);
    let d_ux = max_abs_delta(&ux_candidate, &ux_ref);
    let d_uy = max_abs_delta(&uy_candidate, &uy_ref);
    let d_uz = max_abs_delta(&uz_candidate, &uz_ref);
    println!(
        "WVOF-O1 g=None hydrodynamic bit identity steps=120 max_delta rho={d_rho:.1e} ux={d_ux:.1e} uy={d_uy:.1e} uz={d_uz:.1e}"
    );
}
