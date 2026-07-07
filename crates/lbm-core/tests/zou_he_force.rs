//! ANOM-P4-021: force-consistent Zou-He reconstruction.
//!
//! Velocity faces prescribe the physical Guo velocity. The boundary
//! populations must therefore close on raw momentum `rho*u - F/2`, otherwise
//! the subsequent moment refresh reports `u + F/(2*rho)`.

use lbm_core::prelude::*;

fn opposite(face: Face) -> Face {
    face.opposite()
}

fn face_positions(face: Face, dims: [usize; 3]) -> Vec<[usize; 3]> {
    let a = face.axis();
    let fixed = if face.is_neg() { 0 } else { dims[a] - 1 };
    let (t1, t2) = face.tangents();
    let mut out = Vec::new();
    for c2 in 0..dims[t2] {
        for c1 in 0..dims[t1] {
            let mut p = [0usize; 3];
            p[a] = fixed;
            p[t1] = c1;
            p[t2] = c2;
            out.push(p);
        }
    }
    out
}

fn velocity_for_face(face: Face, d: usize) -> [f64; 3] {
    let n = face.n_in();
    let (t1, t2) = face.tangents();
    let mut u = [0.0; 3];
    u[face.axis()] = 0.035 * n[face.axis()] as f64;
    u[t1] = -0.011;
    if d == 3 {
        u[t2] = 0.007;
    }
    u
}

fn assert_forced_velocity_face<L: Lattice>(dims: [usize; 3], force: [f64; 3], gravity: [f64; 3]) {
    let face = Face::XNeg;
    let out = opposite(face);
    let u_bc = velocity_for_face(face, L::D);
    let mut faces = [FaceBC::Closed; 6];
    faces[face.index()] = FaceBC::Velocity { u: u_bc };
    faces[out.index()] = FaceBC::Pressure { rho: 1.0 };
    let mut periodic = [true; 3];
    periodic[face.axis()] = false;
    if L::D == 2 {
        periodic[2] = false;
    }
    let spec = GlobalSpec::<f64> {
        dims,
        nu: 0.05,
        periodic,
        faces,
        force,
        ..Default::default()
    };
    let mut s: Solver<L, f64, CpuScalar, LocalPeriodic> = Solver::new(
        &spec,
        &[],
        &[],
        [1, 1, 1],
        CpuScalar::default(),
        LocalPeriodic,
    );
    s.set_gravity(gravity);
    s.run(3);

    let mut max_velocity_err = 0.0f64;
    for p in face_positions(face, dims) {
        let got = s.u(p[0], p[1], p[2]);
        for a in 0..L::D {
            max_velocity_err = max_velocity_err.max((got[a] - u_bc[a]).abs());
        }
    }
    println!(
        "ANOM-P4-021 {} forced velocity face: max |u-u_bc|={max_velocity_err:.3e}",
        std::any::type_name::<L>()
    );
    assert!(
        max_velocity_err <= 2.0e-14,
        "{} forced Zou-He velocity face applied physical u with error {max_velocity_err:.3e}",
        std::any::type_name::<L>()
    );
}

#[test]
fn d2q9_velocity_face_uses_guo_half_force_velocity() {
    assert_forced_velocity_face::<D2Q9>([10, 8, 1], [1.0e-6, -2.0e-6, 0.0], [-3.0e-7, 4.0e-7, 0.0]);
}

#[test]
fn d3q19_velocity_face_uses_guo_half_force_velocity() {
    assert_forced_velocity_face::<D3Q19>(
        [10, 8, 6],
        [1.0e-6, -2.0e-6, 3.0e-6],
        [-3.0e-7, 4.0e-7, -5.0e-7],
    );
}

#[test]
fn d3q27_velocity_face_uses_guo_half_force_velocity() {
    assert_forced_velocity_face::<D3Q27>(
        [10, 8, 6],
        [1.0e-6, -2.0e-6, 3.0e-6],
        [-3.0e-7, 4.0e-7, -5.0e-7],
    );
}
