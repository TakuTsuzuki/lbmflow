//! W-VOF O1: conservative Allen-Cahn phase-field transport under prescribed
//! velocity, with no hydrodynamic coupling.

use lbm_core::lattice::D3Q19;
use lbm_core::prelude::*;
use std::io::Write;
use std::path::PathBuf;

type Sim = Solver<D3Q19, f64, CpuScalar, LocalPeriodic>;

fn periodic_spec(n: usize) -> GlobalSpec<f64> {
    GlobalSpec {
        dims: [n, n, n],
        periodic: [true, true, true],
        ..Default::default()
    }
}

fn build(n: usize) -> Sim {
    Solver::new(
        &periodic_spec(n),
        &[],
        &[],
        [1, 1, 1],
        CpuScalar::default(),
        LocalPeriodic,
    )
}

fn params() -> PhaseFieldParams<f64> {
    PhaseFieldParams::new(4.0, 0.04)
}

fn slab_phi(nx: usize, w: f64, x: usize) -> f64 {
    let xc = x as f64 + 0.5;
    let x1 = nx as f64 * 0.25;
    let x2 = nx as f64 * 0.75;
    0.5 * ((2.0 * (xc - x1) / w).tanh() - (2.0 * (xc - x2) / w).tanh())
}

fn periodic_abs_delta(a: usize, b: usize, n: usize) -> f64 {
    let d = if a >= b { a - b } else { b - a };
    let wrapped = n - d;
    if d < wrapped {
        d as f64
    } else {
        wrapped as f64
    }
}

fn droplet_phi(n: usize, w: f64, radius: f64, x: usize, y: usize, z: usize) -> f64 {
    let c = n / 2;
    let dx = periodic_abs_delta(x, c, n);
    let dy = periodic_abs_delta(y, c, n);
    let dz = periodic_abs_delta(z, c, n);
    let r = (dx * dx + dy * dy + dz * dz).sqrt();
    0.5 * (1.0 - (2.0 * (r - radius) / w).tanh())
}

fn l2_rel(a: &[f64], b: &[f64]) -> f64 {
    let mut num = 0.0;
    let mut den = 0.0;
    for (&x, &y) in a.iter().zip(b) {
        num += (x - y) * (x - y);
        den += y * y;
    }
    (num / den).sqrt()
}

fn x_line(field: &[f64], n: usize, y: usize, z: usize) -> Vec<f64> {
    (0..n).map(|x| field[(z * n + y) * n + x]).collect()
}

fn interface_width_right(line: &[f64], center: usize) -> f64 {
    let crossing = |target: f64| -> f64 {
        let mut prev = line[center];
        for x in center + 1..line.len() {
            let cur = line[x];
            if prev >= target && cur <= target {
                let denom = prev - cur;
                let frac = if denom > 0.0 {
                    (prev - target) / denom
                } else {
                    0.0
                };
                return (x - 1) as f64 + frac;
            }
            prev = cur;
        }
        line.len() as f64
    };
    crossing(0.1) - crossing(0.9)
}

fn artifact_path() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../target/wvof_o1/droplet_profile_before_after.csv")
}

#[test]
fn phase_field_params_enforce_frozen_validity_domain() {
    assert!(params().validate().is_ok());
    assert!(PhaseFieldParams {
        interface_width: 3.9,
        mobility: 0.04,
        clipping_policy: ClippingPolicy::Off,
    }
    .validate()
    .is_err());
    assert!(PhaseFieldParams {
        interface_width: 4.0,
        mobility: 0.0,
        clipping_policy: ClippingPolicy::Off,
    }
    .validate()
    .is_err());
}

#[test]
fn g_none_keeps_existing_hydrodynamic_path_bit_identical() {
    let spec = periodic_spec(12);
    let mut a: Sim = Solver::new(
        &spec,
        &[],
        &[],
        [1, 1, 1],
        CpuScalar::default(),
        LocalPeriodic,
    );
    let mut b: Sim = Solver::new(
        &spec,
        &[],
        &[],
        [1, 1, 1],
        CpuScalar::default(),
        LocalPeriodic,
    );
    a.init_with(|x, y, z| {
        let s = ((x + 2 * y + 3 * z) as f64).sin();
        (1.0 + 1.0e-4 * s, [1.0e-3 * s, -7.0e-4 * s, 5.0e-4 * s])
    });
    b.init_with(|x, y, z| {
        let s = ((x + 2 * y + 3 * z) as f64).sin();
        (1.0 + 1.0e-4 * s, [1.0e-3 * s, -7.0e-4 * s, 5.0e-4 * s])
    });
    a.run(25);
    b.run(25);
    assert_eq!(a.gather_rho(), b.gather_rho());
    assert_eq!(a.gather_ux(), b.gather_ux());
    assert_eq!(a.gather_uy(), b.gather_uy());
    assert_eq!(a.gather_uz(), b.gather_uz());
    for q in 0..D3Q19::Q {
        assert_eq!(a.gather_f(q), b.gather_f(q), "f[{q}] changed with g=None");
    }
}

#[test]
fn flat_interface_at_rest_holds_tanh_profile() {
    let n = 40;
    let p = params();
    let mut sim = build(n);
    let phi0: Vec<f64> = (0..n)
        .flat_map(|_z| (0..n).flat_map(move |_y| (0..n).map(move |x| slab_phi(n, 4.0, x))))
        .collect();
    sim.enable_phase_field_prescribed_velocity(p, &phi0, |_, _, _| [0.0; 3])
        .unwrap();
    let before = x_line(&sim.gather_phi().unwrap(), n, n / 2, n / 2);
    for _ in 0..120 {
        sim.phase_field_step_prescribed_velocity(p, |_, _, _| [0.0; 3])
            .unwrap();
    }
    let after_field = sim.gather_phi().unwrap();
    let after = x_line(&after_field, n, n / 2, n / 2);
    // The profile is sampled on a D3Q19 second-order stencil with W=4, so
    // O((dx/W)^2) gives a 6.25% scale. The 8% band leaves room for the two
    // periodic interfaces interacting weakly through their exponential tails.
    let err = l2_rel(&after, &before);
    assert!(err < 0.08, "flat-interface profile L2_rel={err:.4e}");
    let diag = sim
        .phase_field_step_prescribed_velocity(p, |_, _, _| [0.0; 3])
        .unwrap();
    assert!(
        diag.min_phi > -0.03 && diag.max_phi < 1.03,
        "phi overshoot [{}, {}]",
        diag.min_phi,
        diag.max_phi
    );
}

#[test]
fn diagonal_periodic_droplet_advection_conserves_mass_and_profile() {
    let n = 32;
    let p = params();
    let velocity = [0.032, 0.032, 0.032];
    let mut sim = build(n);
    let phi0: Vec<f64> = (0..n)
        .flat_map(|z| {
            (0..n).flat_map(move |y| (0..n).map(move |x| droplet_phi(n, 4.0, 8.0, x, y, z)))
        })
        .collect();
    sim.enable_phase_field_prescribed_velocity(p, &phi0, |_, _, _| velocity)
        .unwrap();
    let initial = sim.gather_phi().unwrap();
    let mass0: f64 = initial.iter().sum();
    let before_line = x_line(&initial, n, n / 2, n / 2);

    let mut sign_flips = 0usize;
    let mut last_sign = 0i8;
    for _ in 0..1000 {
        let diag = sim
            .phase_field_step_prescribed_velocity(p, |_, _, _| velocity)
            .unwrap();
        let drift = diag.total_phi - mass0;
        let sign = if drift > 0.0 {
            1
        } else if drift < 0.0 {
            -1
        } else {
            0
        };
        let active = drift.abs() > 1.0e-10 * mass0;
        if active && sign != 0 && last_sign != 0 && sign != last_sign {
            sign_flips += 1;
        }
        if active && sign != 0 {
            last_sign = sign;
        }
        assert!(
            diag.min_phi > -0.04 && diag.max_phi < 1.04,
            "phi overshoot [{}, {}]",
            diag.min_phi,
            diag.max_phi
        );
    }

    let final_phi = sim.gather_phi().unwrap();
    let mass1: f64 = final_phi.iter().sum();
    let mass_drift = (mass1 - mass0).abs() / mass0;
    assert!(
        mass_drift < 1.0e-3,
        "droplet phi mass drift {mass_drift:.4e} exceeds 0.1%/1000 steps"
    );
    assert!(
        sign_flips < 25,
        "mass drift sign flipped {sign_flips} times, indicating non-monotone growing drift"
    );

    let after_line = x_line(&final_phi, n, n / 2, n / 2);
    // Return-to-start under U=(0.032,0.032,0.032) over 1000 steps translates
    // by one 32-cell period. The same second-order W=4 profile allowance as
    // the flat test is widened to 14% for diagonal lattice anisotropy and
    // finite-curvature correction on R/W=2.
    let profile_err = l2_rel(&after_line, &before_line);
    assert!(
        profile_err < 0.14,
        "returned droplet profile L2_rel={profile_err:.4e}"
    );
    let width = interface_width_right(&after_line, n / 2);
    assert!(
        width >= 0.5 * p.interface_width && width <= 2.0 * p.interface_width,
        "interface width {width:.3} outside [0.5W,2W]"
    );

    let path = artifact_path();
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent).unwrap();
    }
    let mut file = std::fs::File::create(&path).unwrap();
    writeln!(file, "x,phi_before,phi_after").unwrap();
    for x in 0..n {
        writeln!(file, "{x},{},{}", before_line[x], after_line[x]).unwrap();
    }
    println!(
        "W-VOF O1 droplet: mass_drift={mass_drift:.4e}, profile_l2={profile_err:.4e}, width={width:.3}, artifact={}",
        path.display()
    );
}
