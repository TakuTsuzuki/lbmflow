//! Native strain-rate observable checks for FR-STRESS-01.

use lbm_core::prelude::*;

type S2<H> = Solver<D2Q9, f64, CpuScalar, H>;

fn max_abs(a: &[f64], b: &[f64]) -> f64 {
    a.iter()
        .zip(b)
        .map(|(x, y)| (x - y).abs())
        .fold(0.0, f64::max)
}

fn run_to_steady(s: &mut S2<LocalPeriodic>, check_every: usize, tol: f64, max_steps: usize) {
    let mut prev = s.gather_ux();
    for _ in (0..max_steps).step_by(check_every) {
        s.run(check_every);
        let cur = s.gather_ux();
        let scale = cur.iter().map(|v| v.abs()).fold(0.0f64, f64::max).max(1.0);
        if max_abs(&cur, &prev) / scale < tol {
            return;
        }
        prev = cur;
    }
    panic!("strain-rate fixture did not reach steady state");
}

fn channel_solver(
    nx: usize,
    ny: usize,
    nu: f64,
    force: [f64; 3],
    top_u: [f64; 3],
) -> S2<LocalPeriodic> {
    let spec = GlobalSpec {
        dims: [nx, ny, 1],
        nu,
        collision: CollisionKind::Trt {
            magic: CollisionKind::MAGIC_STD,
        },
        periodic: [true, false, false],
        force,
        ..Default::default()
    };
    let mut walls = WallSpec::default();
    walls.is_wall[Face::YNeg.index()] = true;
    walls.is_wall[Face::YPos.index()] = true;
    walls.u[Face::YPos.index()] = top_u;
    let (solid, wall_u) = build_wall_rims(D2Q9::D, spec.dims, &walls);
    Solver::new(
        &spec,
        &solid,
        &wall_u,
        [1, 1, 1],
        CpuScalar::default(),
        LocalPeriodic,
    )
}

#[test]
fn couette_strain_rate_matches_half_way_wall_gradient() {
    let u_wall = 0.1;
    let ny = 10;
    let h = (ny - 2) as f64;
    let mut s = channel_solver(4, ny, (1.0 - 0.5) / 3.0, [0.0; 3], [u_wall, 0.0, 0.0]);
    run_to_steady(&mut s, 500, 1.0e-11, 200_000);

    let strain = s.gather_strain_rate();
    let shear = s.gather_shear_rate();
    let expect_sxy = 0.5 * u_wall / h;
    let mut err_sxy = 0.0f64;
    let mut err_gamma_consistency = 0.0f64;
    for y in 1..ny - 1 {
        let i = y * 4;
        err_sxy = err_sxy.max((strain[i][3] - expect_sxy).abs());
        let s2 = strain[i][0] * strain[i][0]
            + strain[i][1] * strain[i][1]
            + strain[i][2] * strain[i][2]
            + 2.0
                * (strain[i][3] * strain[i][3]
                    + strain[i][4] * strain[i][4]
                    + strain[i][5] * strain[i][5]);
        err_gamma_consistency = err_gamma_consistency.max((shear[i] - (2.0 * s2).sqrt()).abs());
        assert_eq!(strain[i][2], 0.0);
        assert_eq!(strain[i][4], 0.0);
        assert_eq!(strain[i][5], 0.0);
    }
    eprintln!(
        "couette strain errors: Sxy={err_sxy:e}, gamma consistency={err_gamma_consistency:e}"
    );
    assert!(err_sxy <= 3.0e-15, "Couette Sxy error = {err_sxy:e}");
    assert!(
        err_gamma_consistency <= 1.0e-15,
        "Couette gamma_dot consistency error = {err_gamma_consistency:e}"
    );
    assert_eq!(
        strain[0], [0.0; 6],
        "solid rim cell must report zero strain"
    );
    assert_eq!(shear[0], 0.0, "solid rim cell must report zero shear rate");
}

#[test]
fn forced_poiseuille_shear_rate_uses_rev4_force_sign() {
    let ny = 10;
    let nu = 0.1;
    let g = 1.0e-6;
    let h = (ny - 2) as f64;
    let mut s = channel_solver(4, ny, nu, [g, 0.0, 0.0], [0.0; 3]);
    run_to_steady(&mut s, 500, 1.0e-11, 200_000);

    let shear = s.gather_shear_rate();
    let mut err = 0.0f64;
    for y in 1..ny - 1 {
        let yw = y as f64 - 0.5;
        let expect = (g / (2.0 * nu) * (h - 2.0 * yw)).abs();
        err = err.max((shear[y * 4] - expect).abs());
    }
    eprintln!("poiseuille shear-rate error with Pi_force = -0.5(uF+Fu): {err:e}");
    assert!(err <= 3.1e-13, "Poiseuille gamma_dot error = {err:e}");
}

#[test]
fn strain_rate_gather_is_bit_identical_across_inprocess_split() {
    let spec = GlobalSpec {
        dims: [12, 10, 1],
        nu: 0.08,
        collision: CollisionKind::Bgk,
        periodic: [true, true, false],
        force: [1.0e-6, -2.0e-6, 0.0],
        ..Default::default()
    };
    let mut mono: S2<LocalPeriodic> = Solver::new(
        &spec,
        &[],
        &[],
        [1, 1, 1],
        CpuScalar::default(),
        LocalPeriodic,
    );
    let mut split: S2<InProcess> =
        Solver::new(&spec, &[], &[], [2, 2, 1], CpuScalar::default(), InProcess);
    mono.set_body_force_field(|x, y, _| [1.0e-7 * x as f64, -2.0e-7 * y as f64, 0.0]);
    split.set_body_force_field(|x, y, _| [1.0e-7 * x as f64, -2.0e-7 * y as f64, 0.0]);
    for _ in 0..25 {
        mono.step();
        split.step();
    }

    let a = mono.gather_strain_rate();
    let b = split.gather_strain_rate();
    for (i, (x, y)) in a.iter().zip(&b).enumerate() {
        for c in 0..6 {
            assert_eq!(
                x[c].to_bits(),
                y[c].to_bits(),
                "strain[{i}][{c}] differs: {:e} vs {:e}",
                x[c],
                y[c]
            );
        }
    }
    let (ga, gb) = (mono.gather_shear_rate(), split.gather_shear_rate());
    for (i, (x, y)) in ga.iter().zip(&gb).enumerate() {
        assert_eq!(
            x.to_bits(),
            y.to_bits(),
            "gamma_dot[{i}] differs: {x:e} vs {y:e}"
        );
    }
}
