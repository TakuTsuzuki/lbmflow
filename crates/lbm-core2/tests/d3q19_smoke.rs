//! D3Q19 smoke tests: the same generic kernels that reproduce V1 in 2D must
//! run a 3D lattice with sane physics (conservation, viscous decay,
//! bounce-back walls). Full 3D validation (T15: TGV3D convergence, sphere
//! drag, 3D cavity) is M-C scope.

use lbm_core2::lattice::D3Q19;
use lbm_core2::prelude::*;

type Solver3<T> = Solver<D3Q19, T, CpuScalar, LocalPeriodic>;

fn kinetic_energy(s: &Solver3<f64>) -> f64 {
    let (ux, uy, uz) = (s.gather_ux(), s.gather_uy(), s.gather_uz());
    ux.iter()
        .zip(&uy)
        .zip(&uz)
        .map(|((a, b), c)| a * a + b * b + c * c)
        .sum()
}

fn tgv3d(n: usize, nu: f64) -> Solver3<f64> {
    let spec = GlobalSpec::<f64> {
        dims: [n, n, n],
        nu,
        periodic: [true, true, true],
        ..Default::default()
    };
    let mut s: Solver3<f64> = Solver::new(
        &spec,
        &[],
        &[],
        [1, 1, 1],
        CpuScalar::default(),
        LocalPeriodic,
    );
    let k = 2.0 * std::f64::consts::PI / n as f64;
    let u0 = 0.03;
    s.init_with(move |x, y, z| {
        let (xx, yy, zz) = (k * x as f64, k * y as f64, k * z as f64);
        // Classic 3D Taylor-Green vortex initial field.
        let ux = u0 * xx.sin() * yy.cos() * zz.cos();
        let uy = -u0 * xx.cos() * yy.sin() * zz.cos();
        (1.0, [ux, uy, 0.0])
    });
    s
}

#[test]
fn d3q19_tgv_conserves_mass_and_momentum() {
    let mut s = tgv3d(24, 0.02);
    let m0 = s.total_mass();
    let p0 = s.total_momentum();
    let e0 = kinetic_energy(&s);
    s.run(200);
    let m1 = s.total_mass();
    let p1 = s.total_momentum();
    let e1 = kinetic_energy(&s);
    let n3 = (24.0f64).powi(3);
    assert!(
        (m1 - m0).abs() / n3 < 1e-14,
        "mass drift: {m0} -> {m1}"
    );
    for a in 0..3 {
        assert!(
            (p1[a] - p0[a]).abs() / n3 < 1e-14,
            "momentum[{a}] drift: {} -> {}",
            p0[a],
            p1[a]
        );
    }
    // Viscous decay: energy strictly drops, and everything stays finite.
    assert!(e1 < e0, "kinetic energy must decay: {e0} -> {e1}");
    assert!(e1 > 0.0 && e1.is_finite());
    println!("TGV3D: E {e0:.6e} -> {e1:.6e}, mass drift {:.2e}", (m1 - m0).abs());
}

#[test]
fn d3q19_tgv_decay_rate_tracks_viscosity() {
    // For the 3D TGV at low Re the early-time energy decay follows
    // E(t) ≈ E0 exp(-2 nu k_eff^2 t). Check the measured rate against the
    // analytic one within a coarse band — a physics sanity check, not a
    // validation-grade convergence study (that is T15/M-C).
    let n = 24;
    let nu = 0.04;
    let mut s = tgv3d(n, nu);
    let e0 = kinetic_energy(&s);
    let steps = 150;
    s.run(steps);
    let e1 = kinetic_energy(&s);
    let k = 2.0 * std::f64::consts::PI / n as f64;
    // The classic TGV field has wavevector components in x, y, z with
    // |k|^2 = 3 k^2 for the energy-carrying mode.
    let rate_expect = 2.0 * nu * 3.0 * k * k;
    let rate_measured = -(e1 / e0).ln() / steps as f64;
    let rel = (rate_measured - rate_expect).abs() / rate_expect;
    println!(
        "TGV3D decay: measured {rate_measured:.4e}, analytic {rate_expect:.4e}, rel {rel:.3}"
    );
    assert!(
        rel < 0.15,
        "decay rate off by {rel:.3} (measured {rate_measured:e}, expect {rate_expect:e})"
    );
}

#[test]
fn d3q19_moving_lid_box_is_stable_and_conserves_mass() {
    // 3D box, all faces walls, top lid (y = ny-1) moving in +x: exercises
    // 3D rims, bounce-back with wall velocity, and the momentum probe.
    let n = 16;
    let spec = GlobalSpec::<f64> {
        dims: [n, n, n],
        nu: 0.05,
        periodic: [false, false, false],
        ..Default::default()
    };
    let mut walls = WallSpec::<f64>::default();
    for f in Face::ALL {
        walls.is_wall[f.index()] = true;
    }
    walls.u[Face::YPos.index()] = [0.08, 0.0, 0.0];
    let (solid, wall_u) = build_wall_rims(3, spec.dims, &walls);
    let mut s: Solver3<f64> = Solver::new(
        &spec,
        &solid,
        &wall_u,
        [1, 1, 1],
        CpuScalar::default(),
        LocalPeriodic,
    );
    s.set_force_probe(move |_, y, _| y == n - 1);
    let m0 = s.total_mass();
    s.run(300);
    let m1 = s.total_mass();
    assert!(
        (m1 - m0).abs() / m0 < 1e-13,
        "mass drift in closed box: {m0} -> {m1}"
    );
    // The lid drags fluid: interior velocity below the lid is positive-x.
    let u = s.u(n / 2, n - 2, n / 2);
    assert!(u[0] > 1e-4, "lid must drag fluid, got ux = {}", u[0]);
    // The lid feels a reaction force opposing its motion.
    let f = s.probed_force();
    assert!(f[0] < 0.0, "drag on moving lid must be -x, got {:?}", f);
    assert!(f.iter().all(|c| c.is_finite()));
    println!("3D lid box: ux(below lid) = {:.3e}, lid force = {:?}", u[0], f);
}
