//! W1 wall-treatment diagnostics: read-only y+ / u_tau observable.

use lbm_core::prelude::*;

type S2 = Solver<D2Q9, f64, CpuScalar, LocalPeriodic>;

fn channel(ny: usize, force: f64, nu: f64) -> S2 {
    let spec = GlobalSpec {
        dims: [8, ny, 1],
        nu,
        collision: CollisionKind::Trt {
            magic: CollisionKind::MAGIC_STD,
        },
        periodic: [true, false, true],
        force: [force, 0.0, 0.0],
        ..Default::default()
    };
    let mut walls = WallSpec::default();
    walls.is_wall[Face::YNeg.index()] = true;
    walls.is_wall[Face::YPos.index()] = true;
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

fn assert_close(actual: f64, expected: f64, tol: f64, label: &str) {
    let err = (actual - expected).abs();
    assert!(
        err <= tol,
        "{label}: actual={actual:e}, expected={expected:e}, err={err:e}, tol={tol:e}"
    );
}

#[test]
fn poiseuille_laminar_wall_metric_matches_analytic_shear_band() {
    let ny = 66;
    let nu = 0.1;
    let force = 1.0e-6;
    let mut solver = channel(ny, force, nu);
    let h = (ny - 2) as f64;
    solver.init_with(|_, y, _| {
        if y == 0 || y + 1 == ny {
            (1.0, [0.0; 3])
        } else {
            let y_w = y as f64 - 0.5;
            (1.0, [force * y_w * (h - y_w) / (2.0 * nu), 0.0, 0.0])
        }
    });

    let metrics = solver.gather_wall_metrics();
    assert_eq!(metrics.len(), 2 * solver.dims()[0]);

    let analytic_u_tau = (force * h / 2.0).sqrt();
    let first_node_u_tau = (force * (h - 0.5) / 2.0).sqrt();
    // The observable's laminar branch estimates du/dy at the half-way wall
    // from the first fluid node: u(0.5)/0.5 = G(h-0.5)/(2nu), while the
    // continuum wall shear is G h/(2nu). The relative discretization deficit
    // is therefore 1 - sqrt((h-0.5)/h), used here as the band with roundoff
    // headroom. No fitted constant is introduced.
    let rel_band = (analytic_u_tau - first_node_u_tau) / analytic_u_tau + 1.0e-12;
    for m in metrics {
        assert_eq!(m.source, WallMetricSource::HalfwayRim);
        assert_close(m.y_w, 0.5, 1.0e-15, "half-way wall distance");
        let rel = (m.u_tau - analytic_u_tau).abs() / analytic_u_tau;
        assert!(
            rel <= rel_band,
            "Poiseuille u_tau rel={rel:e} > band={rel_band:e}, metric={m:?}"
        );
    }
}

#[test]
fn quiescent_wall_metrics_have_zero_friction_velocity() {
    let solver = channel(10, 0.0, 0.1);
    let metrics = solver.gather_wall_metrics();
    assert_eq!(metrics.len(), 2 * solver.dims()[0]);
    for m in metrics {
        assert_close(m.u_parallel, 0.0, 0.0, "quiescent u_parallel");
        assert_close(m.u_tau, 0.0, 0.0, "quiescent u_tau");
        assert_close(m.y_plus, 0.0, 0.0, "quiescent y_plus");
        assert_close(m.tau_w, 0.0, 0.0, "quiescent tau_w");
    }
}

#[test]
fn fully_periodic_box_reports_no_wall_metrics() {
    let spec = GlobalSpec {
        dims: [6, 6, 1],
        nu: 0.1,
        periodic: [true, true, true],
        ..Default::default()
    };
    let solver: S2 = Solver::new(
        &spec,
        &[],
        &[],
        [1, 1, 1],
        CpuScalar::default(),
        LocalPeriodic,
    );
    assert!(solver.gather_wall_metrics().is_empty());
}

#[test]
fn bouzidi_wall_distance_controls_y_plus_linearly_for_fixed_u_tau() {
    let nu = 0.1;
    let target_u_tau = 0.02;
    let spec = GlobalSpec {
        dims: [7, 5, 1],
        nu,
        periodic: [true, true, true],
        ..Default::default()
    };
    let mut solver: S2 = Solver::new(
        &spec,
        &[],
        &[],
        [1, 1, 1],
        CpuScalar::default(),
        LocalPeriodic,
    );
    solver.set_solid(2, 2, 0);
    solver.set_solid(5, 2, 0);
    let geom = solver.fields(0).geom;
    let cell_a = geom.pidx(1, 2, 0);
    let wall_a = geom.pidx(2, 2, 0);
    let cell_b = geom.pidx(4, 2, 0);
    let wall_b = geom.pidx(5, 2, 0);
    solver.set_bouzidi_links(
        0,
        Some(BouzidiLinks::new(vec![
            BouzidiLink {
                cell: cell_a as u32,
                q: 1,
                qd: 0.25,
                has_second: false,
                wall_ref: wall_a as u32,
            },
            BouzidiLink {
                cell: cell_b as u32,
                q: 1,
                qd: 0.5,
                has_second: false,
                wall_ref: wall_b as u32,
            },
        ])),
    );
    solver.init_with(|x, y, _| {
        let speed = if x == 1 && y == 2 {
            target_u_tau * target_u_tau * 0.25 / nu
        } else if x == 4 && y == 2 {
            target_u_tau * target_u_tau * 0.5 / nu
        } else {
            0.0
        };
        (1.0, [0.0, speed, 0.0])
    });

    let metrics = solver.gather_wall_metrics();
    let gi_a = 2 * spec.dims[0] + 1;
    let gi_b = 2 * spec.dims[0] + 4;
    let a = metrics.iter().find(|m| m.cell_index == gi_a).unwrap();
    let b = metrics.iter().find(|m| m.cell_index == gi_b).unwrap();
    assert_eq!(a.source, WallMetricSource::Bouzidi);
    assert_eq!(b.source, WallMetricSource::Bouzidi);
    assert_close(a.y_w, 0.25, 1.0e-15, "Bouzidi qd distance A");
    assert_close(b.y_w, 0.5, 1.0e-15, "Bouzidi qd distance B");
    assert_close(a.u_tau, target_u_tau, 1.0e-14, "Bouzidi u_tau A");
    assert_close(b.u_tau, target_u_tau, 1.0e-14, "Bouzidi u_tau B");
    assert_close(
        b.y_plus / a.y_plus,
        2.0,
        1.0e-13,
        "Bouzidi y+ distance ratio",
    );
}
