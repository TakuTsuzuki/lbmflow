//! Validation T10: configuration errors, near-limit TRT stability, and
//! unsupported solid placement on open edges.

use lbm_core::prelude::*;

fn expect_err(cfg: SimConfig<f64>, want: fn(&ConfigError) -> bool, label: &str) {
    match cfg.build() {
        Ok(_) => panic!("T10 {label}: expected ConfigError, got Ok"),
        Err(err) => assert!(want(&err), "T10 {label}: wrong error = {err:?}"),
    }
}

#[test]
fn t10_every_config_error_path_is_reported() {
    expect_err(
        SimConfig {
            nu: 0.0,
            ..Default::default()
        },
        |e| matches!(e, ConfigError::NonPositiveViscosity { nu } if *nu == 0.0),
        "tau<=0.5",
    );
    expect_err(
        SimConfig {
            nx: 2,
            ny: 64,
            ..Default::default()
        },
        |e| matches!(e, ConfigError::DomainTooSmall { nx, ny } if *nx == 2 && *ny == 64),
        "nx<3",
    );
    expect_err(
        SimConfig {
            edges: Edges {
                left: EdgeBC::Periodic,
                right: EdgeBC::BounceBack,
                bottom: EdgeBC::BounceBack,
                top: EdgeBC::BounceBack,
            },
            ..Default::default()
        },
        |e| matches!(e, ConfigError::UnpairedPeriodic { axis } if *axis == "x"),
        "unpaired periodic x",
    );
    expect_err(
        SimConfig {
            edges: Edges {
                left: EdgeBC::VelocityInlet { u: [0.02, 0.0] },
                right: EdgeBC::PressureOutlet { rho: 1.0 },
                bottom: EdgeBC::VelocityInlet { u: [0.0, 0.02] },
                top: EdgeBC::BounceBack,
            },
            ..Default::default()
        },
        |e| matches!(e, ConfigError::AdjacentOpenEdges),
        "Zou-He orthogonal open edge violation",
    );
    expect_err(
        SimConfig {
            edges: Edges {
                left: EdgeBC::BounceBack,
                right: EdgeBC::BounceBack,
                bottom: EdgeBC::MovingWall {
                    u: [MAX_SPEED + 0.01, 0.0],
                },
                top: EdgeBC::BounceBack,
            },
            ..Default::default()
        },
        |e| matches!(e, ConfigError::VelocityTooHigh { speed } if *speed > MAX_SPEED),
        "velocity too high",
    );
    expect_err(
        SimConfig {
            edges: Edges {
                left: EdgeBC::PressureOutlet { rho: 0.0 },
                right: EdgeBC::PressureOutlet { rho: 1.0 },
                bottom: EdgeBC::BounceBack,
                top: EdgeBC::BounceBack,
            },
            ..Default::default()
        },
        |e| matches!(e, ConfigError::NonPositiveDensity { rho } if *rho == 0.0),
        "non-positive density",
    );
}

#[test]
fn t10_tau_051_trt_cavity_runs_without_nan() {
    let mut sim: Simulation<f64> = SimConfig {
        nx: 34,
        ny: 34,
        nu: (0.51 - 0.5) / 3.0,
        collision: Collision::default(),
        edges: Edges {
            left: EdgeBC::BounceBack,
            right: EdgeBC::BounceBack,
            bottom: EdgeBC::BounceBack,
            top: EdgeBC::MovingWall { u: [0.1, 0.0] },
        },
        ..Default::default()
    }
    .build()
    .unwrap();
    sim.run(10_000);
    let mut max_abs = 0.0f64;
    for y in 0..sim.ny() {
        for x in 0..sim.nx() {
            let rho = sim.rho(x, y);
            let ux = sim.ux(x, y);
            let uy = sim.uy(x, y);
            max_abs = max_abs.max(rho.abs()).max(ux.abs()).max(uy.abs());
            assert!(
                rho.is_finite() && ux.is_finite() && uy.is_finite(),
                "T10 tau=0.51 non-finite at ({x},{y}), rho = {rho:e}, ux = {ux:e}, uy = {uy:e}"
            );
        }
    }
    assert!(max_abs < 10.0, "T10 tau=0.51 max_abs = {max_abs:e}");
}

#[test]
fn t10_set_solid_panics_on_all_open_edges() {
    let cases = [
        (
            "left",
            (0usize, 5usize),
            Edges {
                left: EdgeBC::VelocityInlet { u: [0.02, 0.0] },
                right: EdgeBC::PressureOutlet { rho: 1.0 },
                bottom: EdgeBC::BounceBack,
                top: EdgeBC::BounceBack,
            },
        ),
        (
            "right",
            (15usize, 5usize),
            Edges {
                left: EdgeBC::VelocityInlet { u: [0.02, 0.0] },
                right: EdgeBC::PressureOutlet { rho: 1.0 },
                bottom: EdgeBC::BounceBack,
                top: EdgeBC::BounceBack,
            },
        ),
        (
            "bottom",
            (5usize, 0usize),
            Edges {
                left: EdgeBC::BounceBack,
                right: EdgeBC::BounceBack,
                bottom: EdgeBC::VelocityInlet { u: [0.0, 0.02] },
                top: EdgeBC::PressureOutlet { rho: 1.0 },
            },
        ),
        (
            "top",
            (5usize, 15usize),
            Edges {
                left: EdgeBC::BounceBack,
                right: EdgeBC::BounceBack,
                bottom: EdgeBC::VelocityInlet { u: [0.0, 0.02] },
                top: EdgeBC::PressureOutlet { rho: 1.0 },
            },
        ),
    ];
    for (label, (x, y), edges) in cases {
        let result = std::panic::catch_unwind(move || {
            let mut sim: Simulation<f64> = SimConfig {
                nx: 16,
                ny: 16,
                nu: 0.05,
                edges,
                ..Default::default()
            }
            .build()
            .unwrap();
            sim.set_solid(x, y);
        });
        assert!(
            result.is_err(),
            "T10 set_solid open edge label = {label}, x = {x}, y = {y}, panicked = {}",
            result.is_err()
        );
    }
}
