//! 3D Poiseuille exact-solution sentinels (radar remainder + T15.2 extension).
//!
//! Public API only: native D3Q19, f64, CPU.

use lbm_core::prelude::*;
use std::f64::consts::PI;

type S3 = Solver<D3Q19, f64, CpuScalar, LocalPeriodic>;

fn idx(x: usize, y: usize, z: usize, nx: usize, ny: usize) -> usize {
    (z * ny + y) * nx + x
}

fn l2_rel(actual: &[f64], reference: &[f64]) -> f64 {
    assert_eq!(actual.len(), reference.len());
    let (num, den) = actual
        .iter()
        .zip(reference)
        .fold((0.0, 0.0), |(num, den), (&a, &r)| {
            let d = a - r;
            (num + d * d, den + r * r)
        });
    (num / den).sqrt()
}

fn max_abs_diff(a: &[f64], b: &[f64]) -> f64 {
    assert_eq!(a.len(), b.len());
    a.iter()
        .zip(b)
        .map(|(x, y)| (x - y).abs())
        .fold(0.0, f64::max)
}

fn run_to_steady3(s: &mut S3, check_every: usize, tol: f64, max_steps: usize, what: &str) {
    let mut prev: Option<Vec<f64>> = None;
    let mut elapsed = 0;
    while elapsed < max_steps {
        s.run(check_every);
        elapsed += check_every;
        let mut cur = s.gather_ux();
        cur.extend(s.gather_uy());
        cur.extend(s.gather_uz());
        if let Some(p) = &prev {
            let dmax = max_abs_diff(&cur, p);
            let umax = cur.iter().map(|v| v.abs()).fold(0.0, f64::max);
            if umax > 0.0 && dmax <= tol * umax {
                println!(
                    "{what}: steady after {} steps, rel_delta={:.3e}",
                    s.time(),
                    dmax / umax
                );
                return;
            }
        }
        prev = Some(cur);
    }
    panic!("{what}: did not reach steady state after {max_steps} steps");
}

/// Exact series for the axial velocity in a rectangular duct
/// (-a <= y~ <= a, -b <= z~ <= b) driven by uniform body force `g`.
fn duct_series(y: f64, zt: f64, a: f64, b: f64, g: f64, nu: f64, nmax: usize) -> f64 {
    let pref = 16.0 * a * a * g / (nu * PI * PI * PI);
    let mut sum = 0.0;
    let mut n = 1;
    while n <= nmax {
        let nf = n as f64;
        let kn = nf * PI / (2.0 * a);
        let ratio =
            ((kn * zt.abs()).exp() + (-kn * zt.abs()).exp()) / ((kn * b).exp() + (-kn * b).exp());
        sum += (1.0 - ratio) * (kn * y).sin() / (nf * nf * nf);
        n += 2;
    }
    pref * sum
}

#[test]
fn g1_circular_pipe_poiseuille_matches_bulk_parabola() {
    let (nx, n) = (8usize, 48usize);
    let dims = [nx, n, n];
    let tau = 0.8f64;
    let nu = (tau - 0.5) / 3.0;
    let g = 3.0e-6f64;
    let r_pipe = 23.0f64;
    let center = (n as f64 - 1.0) / 2.0;
    let r2_pipe = r_pipe * r_pipe;
    let radius2 = |y: usize, z: usize| {
        let dy = y as f64 - center;
        let dz = z as f64 - center;
        dy * dy + dz * dz
    };
    let analytic = |r2: f64| g * (r2_pipe - r2) / (4.0 * nu);

    let mut solid = vec![false; dims[0] * dims[1] * dims[2]];
    let wall_u = vec![[0.0; 3]; solid.len()];
    for z in 0..n {
        for y in 0..n {
            let outside = radius2(y, z) > r2_pipe;
            for x in 0..nx {
                solid[idx(x, y, z, nx, n)] = outside;
            }
        }
    }

    let spec = GlobalSpec::<f64> {
        dims,
        nu,
        periodic: [true, false, false],
        force: [g, 0.0, 0.0],
        collision: CollisionKind::Trt {
            magic: CollisionKind::MAGIC_STD,
        },
        ..Default::default()
    };
    let mut s = S3::new(
        &spec,
        &solid,
        &wall_u,
        [1, 1, 1],
        CpuScalar::default(),
        LocalPeriodic,
    );
    s.init_with(|_, y, z| {
        let r2 = radius2(y, z);
        let u = if r2 <= r2_pipe { analytic(r2) } else { 0.0 };
        (1.0, [u, 0.0, 0.0])
    });
    run_to_steady3(&mut s, 500, 2.0e-10, 40_000, "circular pipe");

    let ux = s.gather_ux();
    let uy = s.gather_uy();
    let uz = s.gather_uz();
    let mut actual = Vec::new();
    let mut reference = Vec::new();
    let mut cross = 0.0f64;
    let mut center_u = 0.0f64;
    let mut mid_u = 0.0f64;
    let mut near_wall_u = 0.0f64;
    let mut mirror = 0.0f64;

    for z in 0..n {
        for y in 0..n {
            let r2 = radius2(y, z);
            if r2 <= r2_pipe {
                let mut u = 0.0;
                for x in 0..nx {
                    let i = idx(x, y, z, nx, n);
                    u += ux[i];
                    cross = cross.max(uy[i].abs()).max(uz[i].abs());
                }
                u /= nx as f64;
                if r2 < (r_pipe - 2.0).powi(2) {
                    actual.push(u);
                    reference.push(analytic(r2));
                }
                if r2 < 1.0 {
                    center_u = center_u.max(u);
                } else if (r2.sqrt() - 0.5 * r_pipe).abs() < 0.75 {
                    mid_u = mid_u.max(u);
                } else if (r2.sqrt() - (r_pipe - 1.5)).abs() < 0.75 {
                    near_wall_u = near_wall_u.max(u);
                }
                let ym = n - 1 - y;
                let zm = n - 1 - z;
                if radius2(ym, z) <= r2_pipe {
                    mirror = mirror.max((u - ux[idx(0, ym, z, nx, n)]).abs());
                }
                if radius2(y, zm) <= r2_pipe {
                    mirror = mirror.max((u - ux[idx(0, y, zm, nx, n)]).abs());
                }
            }
        }
    }

    let l2 = l2_rel(&actual, &reference);
    let umax_ref = analytic(0.0);
    let cross_rel = cross / umax_ref;
    let mirror_rel = mirror / umax_ref;
    println!(
        "G1 circular pipe N={n} L={nx} R={r_pipe:.1}: L2rel={l2:.6e}, samples={}, \
         umax_ref={umax_ref:.6e}, cross_rel={cross_rel:.3e}, mirror_rel={mirror_rel:.3e}",
        actual.len()
    );
    assert!(l2 <= 3.0e-2, "circular pipe bulk L2rel={l2:.6e} > 3%");
    assert!(
        center_u > mid_u && mid_u > near_wall_u,
        "pipe behavior anchor failed: center={center_u:.6e}, mid={mid_u:.6e}, near_wall={near_wall_u:.6e}"
    );
    assert!(
        cross_rel <= 1.0e-4,
        "pipe cross-flow too large: cross/umax={cross_rel:.3e}"
    );
    assert!(
        mirror_rel <= 1.0e-3,
        "pipe mirror symmetry drift too large: {mirror_rel:.3e}"
    );
}

#[test]
fn g2_square_duct_t15_2_canary_matches_series() {
    let (nx, n) = (4usize, 24usize);
    let (ny, nz) = (n + 2, n + 2);
    let (a, b) = (n as f64 / 2.0, n as f64 / 2.0);
    let tau = 0.8f64;
    let nu = (tau - 0.5) / 3.0;
    let g = 5.0e-6f64;
    let mut walls = WallSpec::<f64>::default();
    for face in [Face::YNeg, Face::YPos, Face::ZNeg, Face::ZPos] {
        walls.is_wall[face.index()] = true;
    }
    let spec = GlobalSpec::<f64> {
        dims: [nx, ny, nz],
        nu,
        periodic: [true, false, false],
        force: [g, 0.0, 0.0],
        collision: CollisionKind::Trt {
            magic: CollisionKind::MAGIC_STD,
        },
        ..Default::default()
    };
    let (solid, wall_u) = build_wall_rims(3, spec.dims, &walls);
    let mut s = S3::new(
        &spec,
        &solid,
        &wall_u,
        [1, 1, 1],
        CpuScalar::default(),
        LocalPeriodic,
    );
    s.init_with(|_, y, z| {
        if y == 0 || y == ny - 1 || z == 0 || z == nz - 1 {
            return (1.0, [0.0; 3]);
        }
        let yy = y as f64 - 0.5;
        let zt = z as f64 - 0.5 - b;
        (1.0, [duct_series(yy, zt, a, b, g, nu, 99), 0.0, 0.0])
    });
    run_to_steady3(&mut s, 500, 1.0e-10, 40_000, "square duct canary");

    let ux = s.gather_ux();
    let uy = s.gather_uy();
    let uz = s.gather_uz();
    let mut actual = Vec::with_capacity(n * n);
    let mut reference = Vec::with_capacity(n * n);
    let mut cross = 0.0f64;
    let mut center = 0.0f64;
    let mut near_wall = f64::INFINITY;
    let mut mirror = 0.0f64;

    for z in 1..=n {
        for y in 1..=n {
            let yy = y as f64 - 0.5;
            let zt = z as f64 - 0.5 - b;
            let ana = duct_series(yy, zt, a, b, g, nu, 99);
            let mut u = 0.0;
            for x in 0..nx {
                let i = idx(x, y, z, nx, ny);
                u += ux[i];
                cross = cross.max(uy[i].abs()).max(uz[i].abs());
            }
            u /= nx as f64;
            actual.push(u);
            reference.push(ana);
            if y == n / 2 && z == n / 2 {
                center = u;
            }
            if y == 1 || y == n || z == 1 || z == n {
                near_wall = near_wall.min(u);
            }
            mirror = mirror.max((u - ux[idx(0, ny - 1 - y, z, nx, ny)]).abs());
            mirror = mirror.max((u - ux[idx(0, y, nz - 1 - z, nx, ny)]).abs());
        }
    }

    let l2 = l2_rel(&actual, &reference);
    let umax_ref = duct_series(a, 0.0, a, b, g, nu, 99);
    let cross_rel = cross / umax_ref;
    let mirror_rel = mirror / umax_ref;
    println!(
        "G2 square duct N={n} tau={tau:.1}: L2rel={l2:.6e}, umax_ref={umax_ref:.6e}, \
         cross_rel={cross_rel:.3e}, mirror_rel={mirror_rel:.3e}"
    );
    assert!(l2 <= 5.0e-3, "square duct canary L2rel={l2:.6e} > 5e-3");
    assert!(
        center > 10.0 * near_wall,
        "duct behavior anchor failed: center={center:.6e}, near_wall={near_wall:.6e}"
    );
    assert!(
        cross_rel <= 1.0e-4,
        "duct cross-flow too large: cross/umax={cross_rel:.3e}"
    );
    assert!(
        mirror_rel <= 1.0e-12,
        "duct mirror symmetry drift too large: {mirror_rel:.3e}"
    );
}
