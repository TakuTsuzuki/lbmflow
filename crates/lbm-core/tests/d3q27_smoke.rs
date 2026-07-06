//! D3Q27 stage-1 smoke tests: closed-wall BGK/TRT kernels, wall rims,
//! moving-wall bounce-back, and the explicit unimplemented-open-kind guard.

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
fn d3q27_outflow_still_returns_unsupported_kind_error() {
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
    let err = match Solver::<D3Q27, f64, CpuScalar, LocalPeriodic>::try_new(
        &spec,
        &[],
        &[],
        [1, 1, 1],
        CpuScalar::default(),
        LocalPeriodic,
    ) {
        Ok(_) => panic!("D3Q27 open faces must be rejected during construction"),
        Err(err) => err,
    };
    assert!(
        matches!(
            err,
            SpecError::UnsupportedOpenFaceKind {
                lattice: "D3Q27",
                face,
                ..
            }
            if face == Face::XPos.index()
        ),
        "expected UnsupportedOpenFaceKind(D3Q27, XPos Outflow), got {err:?}"
    );
}

#[test]
fn d3q27_moving_lid_box_produces_nonzero_circulation() {
    let n = 16usize;
    let (spec, solid, wall_u) = all_wall_box(n, Some([0.08, 0.0, 0.0]));
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
    println!(
        "D3Q27 lid box: circulation={circulation:.6e}, ux below lid={:.6e}, mass rel={rel:.3e}",
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
        below_lid[0] > 1.0e-4,
        "D3Q27 moving lid must drag fluid in +x, got ux={}",
        below_lid[0]
    );
}
