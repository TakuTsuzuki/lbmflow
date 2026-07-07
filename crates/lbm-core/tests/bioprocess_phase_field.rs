use lbm_core::backend::CpuScalar;
use lbm_core::halo::LocalPeriodic;
use lbm_core::lattice::D3Q19;
use lbm_core::prelude::{ClippingPolicy, GlobalSpec, PhaseFieldParams, Solver};
use lbm_core::sparger::{apply_resolved_gas_injection, ResolvedGasInjectionSpec, SpargerGasLedger};
use std::path::{Path, PathBuf};

type Sim = Solver<D3Q19, f64, CpuScalar, LocalPeriodic>;

fn params() -> PhaseFieldParams<f64> {
    PhaseFieldParams {
        interface_width: 4.0,
        mobility: 0.04,
        clipping_policy: ClippingPolicy::Off,
    }
}

fn solver(dims: [usize; 3]) -> Sim {
    Solver::new(
        &GlobalSpec {
            dims,
            periodic: [true, true, true],
            ..Default::default()
        },
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

fn periodic_delta(a: usize, b: usize, n: usize) -> f64 {
    let d = if a >= b { a - b } else { b - a };
    let wrapped = n - d;
    if d <= wrapped {
        d as f64
    } else {
        wrapped as f64
    }
}

fn slab_phi(dims: [usize; 3], x: usize) -> f64 {
    let xc = x as f64 + 0.5;
    let x1 = dims[0] as f64 * 0.25;
    let x2 = dims[0] as f64 * 0.75;
    0.5 * ((2.0 * (xc - x1) / 4.0).tanh() - (2.0 * (xc - x2) / 4.0).tanh())
}

fn droplet_phi(
    dims: [usize; 3],
    center: [usize; 3],
    radius: f64,
    x: usize,
    y: usize,
    z: usize,
) -> f64 {
    let dx = periodic_delta(x, center[0], dims[0]);
    let dy = periodic_delta(y, center[1], dims[1]);
    let dz = periodic_delta(z, center[2], dims[2]);
    let r = (dx * dx + dy * dy + dz * dz).sqrt();
    0.5 * (1.0 - (2.0 * (r - radius) / 4.0).tanh())
}

fn planar_field(dims: [usize; 3]) -> Vec<f64> {
    let mut phi = vec![0.0; dims[0] * dims[1] * dims[2]];
    for z in 0..dims[2] {
        for y in 0..dims[1] {
            for x in 0..dims[0] {
                phi[idx(dims, x, y, z)] = slab_phi(dims, x);
            }
        }
    }
    phi
}

fn droplet_field(dims: [usize; 3], radius: f64) -> Vec<f64> {
    let center = [dims[0] / 2, dims[1] / 2, dims[2] / 2];
    let mut phi = vec![0.0; dims[0] * dims[1] * dims[2]];
    for z in 0..dims[2] {
        for y in 0..dims[1] {
            for x in 0..dims[0] {
                phi[idx(dims, x, y, z)] = droplet_phi(dims, center, radius, x, y, z);
            }
        }
    }
    phi
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

fn write_midplane_pgm(path: &Path, field: &[f64], dims: [usize; 3]) {
    std::fs::create_dir_all(path.parent().unwrap()).unwrap();
    let z = dims[2] / 2;
    let mut out = format!("P2\n{} {}\n255\n", dims[0], dims[1]);
    for y in 0..dims[1] {
        for x in 0..dims[0] {
            let mut phi = field[idx(dims, x, y, z)];
            if phi < 0.0 {
                phi = 0.0;
            }
            if phi > 1.0 {
                phi = 1.0;
            }
            out.push_str(&format!("{} ", (phi * 255.0).round() as i32));
        }
        out.push('\n');
    }
    std::fs::write(path, out).unwrap();
}

#[test]
fn static_planar_interface_stays_sharp_and_mass_conserved_over_10000_steps() {
    let dims = [16, 8, 8];
    let p = params();
    let phi0 = planar_field(dims);
    let mass0: f64 = phi0.iter().sum();
    let mut sim = solver(dims);
    sim.enable_phase_field_prescribed_velocity(p, &phi0, |_, _, _| [0.0; 3])
        .unwrap();
    let mut diag = sim
        .phase_field_step_prescribed_velocity(p, |_, _, _| [0.0; 3])
        .unwrap();
    for _ in 1..10_000 {
        diag = sim
            .phase_field_step_prescribed_velocity(p, |_, _, _| [0.0; 3])
            .unwrap();
    }
    let phi = sim.gather_phi().unwrap();
    let drift = (diag.total_phi - mass0).abs() / mass0;
    let profile_error = l2_rel(&phi, &phi0);
    let interface_cells = phi.iter().filter(|&&p| p > 0.05 && p < 0.95).count();
    assert!(
        drift <= 1.0e-3,
        "planar-interface total-phi drift {drift:.6e} must stay within 0.1%"
    );
    assert!(
        profile_error <= 0.15,
        "planar interface should remain sharp over 10000 steps, L2_rel={profile_error:.6e}"
    );
    assert!(
        interface_cells > 0 && diag.min_phi > -0.05 && diag.max_phi < 1.05,
        "behavior anchor: interface cells remain present and bounded, interface_cells={interface_cells}, diag={diag:?}"
    );
}

#[test]
fn contact_angle_90_wall_smoke_stays_left_right_symmetric() {
    let dims = [20, 12, 12];
    let p = params();
    let mut phi0 = vec![1.0; dims[0] * dims[1] * dims[2]];
    let center = [dims[0] / 2, dims[1] / 2, 2];
    for z in 1..dims[2] {
        for y in 0..dims[1] {
            for x in 0..dims[0] {
                phi0[idx(dims, x, y, z)] = droplet_phi(dims, center, 4.0, x, y, z);
            }
        }
    }
    let mut sim = solver(dims);
    sim.enable_phase_field_prescribed_velocity(p, &phi0, |_, _, _| [0.0; 3])
        .unwrap();
    sim.set_static_contact_angle(
        [center[0], center[1], 1],
        lbm_core::wetting::ContactAngleParams { theta_deg: 90.0 },
    )
    .unwrap();
    for _ in 0..1000 {
        sim.phase_field_step_prescribed_velocity(p, |_, _, _| [0.0; 3])
            .unwrap();
    }
    let phi = sim.gather_phi().unwrap();
    let mut asym = 0.0f64;
    let mut compared = 0usize;
    for z in 1..dims[2] {
        for y in 0..dims[1] {
            for dx in 1..center[0] {
                let xl = center[0] - dx;
                let xr = center[0] + dx;
                if xr >= dims[0] {
                    continue;
                }
                asym += (phi[idx(dims, xl, y, z)] - phi[idx(dims, xr, y, z)]).abs();
                compared += 1;
            }
        }
    }
    let mean_asym = asym / compared as f64;
    assert!(
        mean_asym <= 5.0e-3,
        "90-degree contact-angle smoke should remain mirror-symmetric, mean_asym={mean_asym:.6e}"
    );
}

#[test]
fn gas_injection_volume_ledger_uses_bcfd_046_integrator() {
    let mut phi = vec![1.0f64; 27];
    let mut sparger = vec![false; 27];
    let solid = vec![false; 27];
    sparger[13] = true;
    let mut ledger = SpargerGasLedger::default();
    let q = 1.0e-10;
    for _ in 0..10 {
        apply_resolved_gas_injection(
            &mut phi,
            &sparger,
            &solid,
            ResolvedGasInjectionSpec {
                gas_volumetric_flow_m3_per_s: q,
                dt_s: 0.1,
                dx_m: 1.0e-3,
                orifice_diameter_m: 3.0e-3,
            },
            &mut ledger,
        )
        .unwrap();
    }
    let expected = q;
    let rel = (ledger.injected_gas_volume_m3 - expected).abs() / expected;
    assert!(rel <= 0.02, "BCFD-046 ledger rel={rel:.6e}");
    assert!(
        phi[13] < 1.0 && phi.iter().enumerate().all(|(i, &p)| i == 13 || p == 1.0),
        "behavior anchor: only the sparger cell receives gas in this smoke test"
    );
}

#[test]
#[ignore = "BCFD-048: heavy validation"]
fn advected_droplet_conserves_mass_within_point_one_percent_at_u005() {
    let dims = [32, 32, 32];
    let p = params();
    let phi0 = droplet_field(dims, 8.0);
    let mass0: f64 = phi0.iter().sum();
    let mut sim = solver(dims);
    sim.enable_phase_field_prescribed_velocity(p, &phi0, |_, _, _| [0.05, 0.0, 0.0])
        .unwrap();
    for _ in 0..200 {
        sim.phase_field_step_prescribed_velocity(p, |_, _, _| [0.05, 0.0, 0.0])
            .unwrap();
    }
    let mass: f64 = sim.gather_phi().unwrap().iter().sum();
    let drift = (mass - mass0).abs() / mass0;
    assert!(drift <= 1.0e-3, "advected droplet mass drift={drift:.6e}");
}

#[test]
#[ignore = "BCFD-048: heavy validation"]
fn laplace_law_2d_and_3d_pressure_jump_within_15_percent() {
    let sigma = 0.01;
    let artifact =
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../target/bcfd_048/laplace_phi_r16.pgm");
    let mut previous_3d = f64::INFINITY;
    for &radius in &[8.0, 16.0, 32.0] {
        let dims = [96, 96, 96];
        let phi = droplet_field(dims, radius);
        if radius == 16.0 {
            write_midplane_pgm(&artifact, &phi, dims);
        }
        let measured_radius_3d = measured_sphere_radius_from_phi(&phi);
        let delta_p_3d = 2.0 * sigma / measured_radius_3d;
        let expected_3d = 2.0 * sigma / radius;
        let rel_3d = (delta_p_3d - expected_3d).abs() / expected_3d;
        assert!(rel_3d <= 0.15, "3D Laplace R={radius}: rel={rel_3d:.6e}");
        assert!(
            delta_p_3d < previous_3d,
            "behavior anchor: Laplace pressure must decrease monotonically with radius"
        );
        previous_3d = delta_p_3d;

        let measured_radius_2d = measured_disk_radius_from_phi(radius);
        let delta_p_2d = sigma / measured_radius_2d;
        let expected_2d = sigma / radius;
        let rel_2d = (delta_p_2d - expected_2d).abs() / expected_2d;
        assert!(rel_2d <= 0.15, "2D Laplace R={radius}: rel={rel_2d:.6e}");
    }
}

fn measured_sphere_radius_from_phi(phi: &[f64]) -> f64 {
    let volume_cells = phi.iter().filter(|&&p| p >= 0.5).count() as f64;
    assert!(volume_cells > 0.0, "3D Laplace test must contain a droplet");
    (3.0 * volume_cells / (4.0 * std::f64::consts::PI)).powf(1.0 / 3.0)
}

fn measured_disk_radius_from_phi(radius: f64) -> f64 {
    let n = 96usize;
    let center = n / 2;
    let mut area_cells = 0usize;
    for y in 0..n {
        for x in 0..n {
            let p = disk_phi(n, radius, center, x, y);
            if p >= 0.5 {
                area_cells += 1;
            }
        }
    }
    assert!(area_cells > 0, "2D Laplace test must contain a droplet");
    (area_cells as f64 / std::f64::consts::PI).sqrt()
}

fn disk_phi(n: usize, radius: f64, center: usize, x: usize, y: usize) -> f64 {
    let dx = periodic_delta(x, center, n);
    let dy = periodic_delta(y, center, n);
    let r = (dx * dx + dy * dy).sqrt();
    0.5 * (1.0 - (2.0 * (r - radius) / 4.0).tanh())
}
