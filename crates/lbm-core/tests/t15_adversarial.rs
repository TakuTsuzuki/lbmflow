//! T15 adversarial D3Q19/f64 physics attacks plus R-Phase guard boundaries.
//!
//! Public API only: no population-layout assumptions, no kernel mirroring.

use lbm_core::params::MAX_SPEED;
use lbm_core::prelude::*;
use std::f64::consts::PI;

type S3 = Solver<D3Q19, f64, CpuScalar, LocalPeriodic>;

fn run_to_steady3(s: &mut S3, check_every: usize, tol: f64, max_steps: usize) -> bool {
    let mut prev: Option<Vec<f64>> = None;
    let mut elapsed = 0;
    while elapsed < max_steps {
        s.run(check_every);
        elapsed += check_every;
        let mut cur = s.gather_ux();
        cur.extend(s.gather_uy());
        cur.extend(s.gather_uz());
        if let Some(p) = &prev {
            let d = cur
                .iter()
                .zip(p)
                .map(|(a, b)| (a - b).abs())
                .fold(0.0, f64::max);
            let m = cur.iter().map(|v| v.abs()).fold(0.0, f64::max);
            if m > 0.0 && d <= tol * m {
                return true;
            }
        }
        prev = Some(cur);
    }
    false
}

fn duct_series(y: f64, zt: f64, a: f64, b: f64, g: f64, nu: f64, nmax: usize) -> f64 {
    let pref = 16.0 * a * a * g / (nu * PI.powi(3));
    let mut sum = 0.0;
    let mut n = 1;
    while n <= nmax {
        let nf = n as f64;
        let kn = nf * PI / (2.0 * a);
        let ratio =
            ((kn * zt.abs()).exp() + (-kn * zt.abs()).exp()) / ((kn * b).exp() + (-kn * b).exp());
        sum += (1.0 - ratio) * (kn * y).sin() / nf.powi(3);
        n += 2;
    }
    pref * sum
}

fn rectangular_duct_case(
    dims: [usize; 3],
    axial: usize,
    walls: &[Face],
    force: [f64; 3],
    what: &str,
) -> (f64, f64) {
    let nu = 0.12;
    let g = force[axial].abs();
    let mut wall_spec = WallSpec::<f64>::default();
    for &face in walls {
        wall_spec.is_wall[face.index()] = true;
    }
    let mut periodic = [false; 3];
    periodic[axial] = true;
    let spec = GlobalSpec {
        dims,
        nu,
        periodic,
        force,
        ..Default::default()
    };
    let (solid, wall_u) = build_wall_rims(3, dims, &wall_spec);
    let mut s = S3::new(
        &spec,
        &solid,
        &wall_u,
        [1, 1, 1],
        CpuScalar::default(),
        LocalPeriodic,
    );
    assert!(
        run_to_steady3(&mut s, 400, 3e-10, 40_000),
        "{what}: duct did not reach steady state"
    );

    let cross_axes: Vec<usize> = (0..3).filter(|&a| a != axial).collect();
    let (a0, a1) = (cross_axes[0], cross_axes[1]);
    let (n0, n1) = (dims[a0], dims[a1]);
    let (half0, half1) = ((n0 - 2) as f64 / 2.0, (n1 - 2) as f64 / 2.0);
    let umax_ref = duct_series(half0, 0.0, half0, half1, g, nu, 99);
    let ux = s.gather_ux();
    let uy = s.gather_uy();
    let uz = s.gather_uz();
    let vel = [&ux, &uy, &uz][axial];
    let idx = |p: [usize; 3]| (p[2] * dims[1] + p[1]) * dims[0] + p[0];
    let mut err = 0.0f64;
    for i1 in 1..=(n1 - 2) {
        for i0 in 1..=(n0 - 2) {
            let y = i0 as f64 - 0.5;
            let zt = i1 as f64 - 0.5 - half1;
            let ana = duct_series(y, zt, half0, half1, g, nu, 99);
            let mut avg = 0.0;
            for ia in 0..dims[axial] {
                let mut p = [0usize; 3];
                p[axial] = ia;
                p[a0] = i0;
                p[a1] = i1;
                avg += vel[idx(p)];
            }
            avg /= dims[axial] as f64;
            err = err.max((avg - ana).abs());
        }
    }
    let rel = err / umax_ref;
    println!("{what}: L_inf_rel={rel:.3e}, umax_ref={umax_ref:.6e}, dims={dims:?}");
    (rel, umax_ref)
}

#[test]
fn t15_z_degeneracy_breaker_keeps_z_dependent_mode() {
    let (nx, ny, nz) = (32usize, 32usize, 8usize);
    let nu = 0.03;
    let eps = 1.0e-7;
    let kx = 2.0 * PI / nx as f64;
    let kz = 2.0 * PI / nz as f64;
    let spec = GlobalSpec::<f64> {
        dims: [nx, ny, nz],
        nu,
        periodic: [true, true, true],
        ..Default::default()
    };
    let mut s = S3::new(
        &spec,
        &[],
        &[],
        [1, 1, 1],
        CpuScalar::default(),
        LocalPeriodic,
    );
    s.init_with(move |x, y, z| {
        let xf = kx * x as f64;
        let yf = kx * y as f64;
        let zf = kz * z as f64;
        (
            1.0,
            [
                0.01 * xf.sin() * yf.cos() + eps * zf.sin(),
                -0.01 * xf.cos() * yf.sin(),
                eps * (xf + zf).cos(),
            ],
        )
    });
    s.run(120);
    let ux = s.gather_ux();
    let uz = s.gather_uz();
    let idx = |x: usize, y: usize, z: usize| (z * ny + y) * nx + x;
    let mut z_spread = 0.0f64;
    for y in 0..ny {
        for x in 0..nx {
            z_spread = z_spread.max((ux[idx(x, y, 2)] - ux[idx(x, y, 6)]).abs());
        }
    }
    let uzmax = uz.iter().map(|v| v.abs()).fold(0.0, f64::max);
    println!("z-breaker: z_spread={z_spread:.3e}, max|uz|={uzmax:.3e}, eps={eps:.1e}");
    assert!(
        z_spread > 1.0e-9 || uzmax > 1.0e-9,
        "z-dependent perturbation was silently projected away: z_spread={z_spread:.3e}, max|uz|={uzmax:.3e}"
    );
}

#[test]
fn t15_extreme_aspect_ratio_ducts_match_series_light() {
    let (rel_x, _) = rectangular_duct_case(
        [64, 8, 8],
        0,
        &[Face::YNeg, Face::YPos, Face::ZNeg, Face::ZPos],
        [2.0e-6, 0.0, 0.0],
        "64x8x8 x-duct",
    );
    let (rel_z, _) = rectangular_duct_case(
        [8, 8, 64],
        2,
        &[Face::XNeg, Face::XPos, Face::YNeg, Face::YPos],
        [0.0, 0.0, 2.0e-6],
        "8x8x64 z-duct",
    );
    assert!(
        rel_x <= 1.5e-2 && rel_z <= 1.5e-2,
        "extreme-aspect light ducts exceeded frozen 1.5% band: x={rel_x:.3e}, z={rel_z:.3e}"
    );
}

#[test]
#[ignore = "heavier D-8 aspect ratio attack"]
fn t15_extreme_aspect_ratio_ducts_match_series_spec_size() {
    let (rel_x, _) = rectangular_duct_case(
        [128, 8, 8],
        0,
        &[Face::YNeg, Face::YPos, Face::ZNeg, Face::ZPos],
        [2.0e-6, 0.0, 0.0],
        "128x8x8 x-duct",
    );
    let (rel_z, _) = rectangular_duct_case(
        [8, 8, 128],
        2,
        &[Face::XNeg, Face::XPos, Face::YNeg, Face::YPos],
        [0.0, 0.0, 2.0e-6],
        "8x8x128 z-duct",
    );
    assert!(
        rel_x <= 1.5e-2 && rel_z <= 1.5e-2,
        "extreme-aspect spec ducts exceeded frozen 1.5% band: x={rel_x:.3e}, z={rel_z:.3e}"
    );
}

fn schiller_naumann(re: f64) -> f64 {
    24.0 / re * (1.0 + 0.15 * re.powf(0.687))
}

fn offcenter_sphere_drag(d: usize, dims: [usize; 3], re: f64, u_in: f64) -> f64 {
    let [nx, ny, nz] = dims;
    let r = d as f64 / 2.0;
    let nu = u_in * d as f64 / re;
    let mut faces = [FaceBC::Closed; 6];
    faces[Face::XNeg.index()] = FaceBC::Velocity {
        u: [u_in, 0.0, 0.0],
    };
    faces[Face::XPos.index()] = FaceBC::Pressure { rho: 1.0 };
    let spec = GlobalSpec {
        dims,
        nu,
        periodic: [false, true, true],
        faces,
        ..Default::default()
    };
    let mut s = S3::new(
        &spec,
        &[],
        &[],
        [1, 1, 1],
        CpuScalar::default(),
        LocalPeriodic,
    );
    let (cx, cy, cz) = (2.5 * d as f64, ny as f64 * 0.45, nz as f64 * 0.56);
    let r2 = r * r;
    let inside = move |x: usize, y: usize, z: usize| {
        let (dx, dy, dz) = (x as f64 - cx, y as f64 - cy, z as f64 - cz);
        dx * dx + dy * dy + dz * dz <= r2
    };
    for z in 0..nz {
        for y in 0..ny {
            for x in 0..nx {
                if inside(x, y, z) {
                    s.set_solid(x, y, z);
                }
            }
        }
    }
    s.set_force_probe(inside);
    s.init_with(move |x, y, z| (1.0, [if inside(x, y, z) { 0.0 } else { u_in }, 0.0, 0.0]));
    let r_h = r + 0.5;
    let cd_of = |fx: f64| fx / (0.5 * u_in * u_in * PI * r_h * r_h);
    let mut cd_prev = f64::NAN;
    let mut cd = f64::NAN;
    for chunk in 0..80 {
        let mut fx = 0.0;
        for _ in 0..300 {
            s.step();
            fx += s.probed_force()[0];
        }
        cd = cd_of(fx / 300.0);
        assert!(
            cd.is_finite(),
            "off-center sphere diverged at chunk {chunk}"
        );
        if (cd - cd_prev).abs() <= 1.0e-3 * cd.abs() {
            break;
        }
        cd_prev = cd;
    }
    println!(
        "off-center sphere D={d} Re={re} dims={dims:?}: Cd={cd:.4}, steps={}",
        s.time()
    );
    cd
}

#[test]
fn t15_offcenter_sphere_drag_light() {
    let d = 10;
    let re = 20.0;
    let cd = offcenter_sphere_drag(d, [80, 56, 56], re, 0.06);
    let sn_h = schiller_naumann(re * (d as f64 + 1.0) / d as f64);
    let rel = (cd - sn_h).abs() / sn_h;
    assert!(
        rel <= 0.15,
        "off-center light sphere Cd={cd:.4} vs SN_h={sn_h:.4}, rel={rel:.3e} > 15%"
    );
}

#[test]
#[ignore = "heavy off-center sphere D=24"]
fn t15_offcenter_sphere_drag_spec_size() {
    let d = 24;
    let re = 20.0;
    let cd = offcenter_sphere_drag(d, [192, 128, 128], re, 0.05);
    let sn_h = schiller_naumann(re * (d as f64 + 1.0) / d as f64);
    let rel = (cd - sn_h).abs() / sn_h;
    assert!(
        rel <= 0.10,
        "off-center spec sphere Cd={cd:.4} vs SN_h={sn_h:.4}, rel={rel:.3e} > 10%"
    );
}

#[test]
fn t15_mass_conservation_with_all_six_domain_faces_walled() {
    let dims = [18, 14, 12];
    let mut walls = WallSpec::<f64>::default();
    for face in Face::ALL {
        walls.is_wall[face.index()] = true;
    }
    let spec = GlobalSpec {
        dims,
        nu: 0.05,
        periodic: [false, false, false],
        ..Default::default()
    };
    let (solid, wall_u) = build_wall_rims(3, dims, &walls);
    let mut s = S3::new(
        &spec,
        &solid,
        &wall_u,
        [1, 1, 1],
        CpuScalar::default(),
        LocalPeriodic,
    );
    let (m0_fluid, m0_dev) = s.local_mass_partials();
    let m0 = m0_fluid + m0_dev;
    s.run(600);
    let (m1_fluid, m1_dev) = s.local_mass_partials();
    let m1 = m1_fluid + m1_dev;
    let rel = (m1 - m0).abs() / m0.abs();
    println!("six-face walled mass: m0={m0:.12e}, m1={m1:.12e}, rel={rel:.3e}");
    assert!(
        rel <= 1.0e-12,
        "six-face closed cavity mass drift rel={rel:.3e} (m0={m0:.12e}, m1={m1:.12e})"
    );
}

fn x_open_spec(bc: FaceBC<f64>) -> GlobalSpec<f64> {
    let mut faces = [FaceBC::Closed; 6];
    faces[Face::XNeg.index()] = bc;
    faces[Face::XPos.index()] = FaceBC::Outflow;
    GlobalSpec {
        dims: [8, 6, 6],
        nu: 0.05,
        periodic: [false, true, true],
        faces,
        ..Default::default()
    }
}

#[test]
fn t15_guard_boundary_velocity_density_and_convective_pairs() {
    let just_inside_v = x_open_spec(FaceBC::Velocity {
        u: [MAX_SPEED, 0.0, 0.0],
    });
    assert!(
        just_inside_v.validate(3, &[]).is_ok(),
        "velocity exactly at MAX_SPEED must remain legal"
    );
    let just_outside_v = x_open_spec(FaceBC::Velocity {
        u: [MAX_SPEED + 1.0e-12, 0.0, 0.0],
    });
    assert!(
        matches!(
            just_outside_v.validate(3, &[]),
            Err(SpecError::VelocityTooHigh { speed }) if speed > MAX_SPEED
        ),
        "velocity just above MAX_SPEED must Err with VelocityTooHigh"
    );

    let just_inside_rho = x_open_spec(FaceBC::Pressure {
        rho: f64::MIN_POSITIVE,
    });
    assert!(
        just_inside_rho.validate(3, &[]).is_ok(),
        "positive pressure density boundary must remain legal"
    );
    let just_outside_rho = x_open_spec(FaceBC::Pressure { rho: 0.0 });
    assert!(
        matches!(
            just_outside_rho.validate(3, &[]),
            Err(SpecError::NonPositiveDensity { rho }) if rho == 0.0
        ),
        "zero pressure density must Err with NonPositiveDensity"
    );

    let just_inside_conv = x_open_spec(FaceBC::Convective { u_conv: 1.0 });
    assert!(
        just_inside_conv.validate(3, &[]).is_ok(),
        "convective u=1 must remain legal"
    );
    let just_outside_conv = x_open_spec(FaceBC::Convective {
        u_conv: 1.0 + 1.0e-12,
    });
    assert!(
        matches!(
            just_outside_conv.validate(3, &[]),
            Err(SpecError::InvalidConvectiveSpeed { u_conv }) if u_conv > 1.0
        ),
        "convective u just above 1 must Err with InvalidConvectiveSpeed"
    );
}
