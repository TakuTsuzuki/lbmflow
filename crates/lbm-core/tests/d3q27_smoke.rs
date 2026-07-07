//! D3Q27 stage-1 smoke tests: closed-wall BGK/TRT kernels, wall rims,
//! moving-wall bounce-back, and open-face construction smoke coverage.

use lbm_core::prelude::*;

type Solver27 = Solver<D3Q27, f64, CpuScalar, LocalPeriodic>;

fn all_wall_box(n: usize, lid_u: Option<[f64; 3]>) -> (GlobalSpec<f64>, Vec<bool>, Vec<[f64; 3]>) {
    let spec = GlobalSpec::<f64> {
        dims: [n, n, n],
        nu: 0.05,
        periodic: [false, false, false],
        ..Default::default()
    };
    let mut walls = WallSpec::<f64>::default();
    for face in Face::ALL {
        walls.is_wall[face.index()] = true;
    }
    if let Some(u) = lid_u {
        walls.u[Face::YPos.index()] = u;
    }
    let (solid, wall_u) = build_wall_rims(3, spec.dims, &walls);
    (spec, solid, wall_u)
}

#[test]
fn d3q27_closed_box_stays_finite_and_conserves_mass() {
    let n = 14usize;
    let (spec, solid, wall_u) = all_wall_box(n, None);
    let mut s: Solver27 = Solver::new(
        &spec,
        &solid,
        &wall_u,
        [1, 1, 1],
        CpuScalar::default(),
        LocalPeriodic,
    );
    let k = 2.0 * std::f64::consts::PI / n as f64;
    s.init_with(move |x, y, z| {
        let (xf, yf, zf) = (k * x as f64, k * y as f64, k * z as f64);
        (
            1.0 + 1.0e-4 * xf.cos() * yf.cos() * zf.cos(),
            [
                0.01 * xf.sin() * yf.cos() * zf.cos(),
                -0.01 * xf.cos() * yf.sin() * zf.cos(),
                0.005 * xf.cos() * yf.cos() * zf.sin(),
            ],
        )
    });
    let m0 = s.total_mass();
    s.run(400);
    let m1 = s.total_mass();
    let rel = (m1 - m0).abs() / m0.abs();
    for (name, field) in [
        ("rho", s.gather_rho()),
        ("ux", s.gather_ux()),
        ("uy", s.gather_uy()),
        ("uz", s.gather_uz()),
    ] {
        let max_abs = field.iter().map(|v| v.abs()).fold(0.0f64, f64::max);
        assert!(
            max_abs.is_finite(),
            "D3Q27 closed box {name} field is not finite: max_abs={max_abs:e}"
        );
    }
    println!("D3Q27 closed box mass drift rel={rel:.3e} over 400 steps");
    assert!(
        rel <= 1.0e-13,
        "D3Q27 closed-box mass drift rel={rel:.3e} (m0={m0:.12e}, m1={m1:.12e})"
    );
}

#[test]
fn d3q27_outflow_open_face_constructs_and_runs() {
    let mut faces = [FaceBC::Closed; 6];
    faces[Face::XNeg.index()] = FaceBC::Velocity {
        u: [0.02, 0.0, 0.0],
    };
    faces[Face::XPos.index()] = FaceBC::Outflow;
    let spec = GlobalSpec::<f64> {
        dims: [8, 6, 6],
        periodic: [false, true, true],
        faces,
        ..Default::default()
    };
    let mut s = Solver::<D3Q27, f64, CpuScalar, LocalPeriodic>::new(
        &spec,
        &[],
        &[],
        [1, 1, 1],
        CpuScalar::default(),
        LocalPeriodic,
    );
    s.run(20);
    assert!(s.gather_rho().iter().all(|v| v.is_finite()));
}

#[test]
fn d3q27_moving_lid_box_produces_nonzero_circulation() {
    let n = 16usize;
    let lid_speed = 0.08;
    let (spec, solid, wall_u) = all_wall_box(n, Some([lid_speed, 0.0, 0.0]));
    let mut s: Solver27 = Solver::new(
        &spec,
        &solid,
        &wall_u,
        [1, 1, 1],
        CpuScalar::default(),
        LocalPeriodic,
    );
    let m0 = s.total_mass();
    s.run(350);
    let m1 = s.total_mass();
    let rel = (m1 - m0).abs() / m0.abs();
    let z = n / 2;
    let lo = 2usize;
    let hi = n - 3;
    let mut circulation = 0.0;
    for x in lo..hi {
        circulation += s.u(x, hi, z)[0] - s.u(x, lo, z)[0];
    }
    for y in lo..hi {
        circulation += s.u(hi, y, z)[1] - s.u(lo, y, z)[1];
    }
    let below_lid = s.u(n / 2, n - 2, n / 2);
    let ux_rel_lid = below_lid[0] / lid_speed;
    let circulation_rel_lid = circulation.abs() / lid_speed;
    println!(
        "D3Q27 lid box: circulation={circulation:.6e}, circulation/lid={circulation_rel_lid:.3e}, ux below lid={:.6e}, ux/lid={ux_rel_lid:.3e}, lid speed={lid_speed:.3e}, mass rel={rel:.3e}",
        below_lid[0]
    );
    assert!(
        rel <= 1.0e-13,
        "D3Q27 moving-lid mass drift rel={rel:.3e} (m0={m0:.12e}, m1={m1:.12e})"
    );
    assert!(
        below_lid.iter().all(|v| v.is_finite()) && circulation.is_finite(),
        "D3Q27 moving-lid fields must stay finite: u={below_lid:?}, circulation={circulation:e}"
    );
    assert!(
        circulation.abs() > 1.0e-4,
        "D3Q27 moving lid must produce nonzero circulation, got {circulation:e}"
    );
    assert!(
        (0.25..=1.0).contains(&circulation_rel_lid),
        "D3Q27 moving-lid circulation must be O(lid speed): circulation={circulation:e}, lid_speed={lid_speed:e}, circulation/lid={circulation_rel_lid:e}"
    );
    assert!(
        below_lid[0] > 1.0e-4,
        "D3Q27 moving lid must drag fluid in +x, got ux={}",
        below_lid[0]
    );
    assert!(
        (0.25..=1.0).contains(&ux_rel_lid),
        "D3Q27 moving-lid ux must be O(lid speed) and not exceed the lid: ux={}, lid_speed={}, ux/lid={}",
        below_lid[0],
        lid_speed,
        ux_rel_lid
    );
}
