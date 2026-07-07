//! Native-side snapshot exporter for the WASM/native parity lane.
//!
//! This test intentionally does not build or execute WASM. It produces the
//! canonical native f64 field at `test/tmp/native_scenario_field.json` so an
//! operator can compare it with a separately captured WASM f32 snapshot.

mod common;

use common::run_to_steady;
use lbm_core::compat::prelude::*;
use serde::{Deserialize, Serialize};
use std::fs;
use std::path::PathBuf;

const N: usize = 32;
const NU: f64 = 0.1;
const G: f64 = 1.0e-6;
const STEADY_CHECK_EVERY: usize = 500;
const STEADY_TOL: f64 = 1.0e-11;
const STEADY_MAX_STEPS: usize = 200_000;

#[derive(Debug, Serialize, Deserialize)]
struct FieldSnapshot {
    schema_version: u32,
    producer: String,
    scenario_id: String,
    validation_case: String,
    field_order: String,
    scalar: String,
    lattice: String,
    collision: String,
    storage_note: String,
    nx: usize,
    ny: usize,
    steps: u64,
    nu: f64,
    tau: f64,
    force: [f64; 2],
    rho: Vec<f64>,
    ux: Vec<f64>,
    uy: Vec<f64>,
    solid: Vec<bool>,
}

fn workspace_export_path() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .join("test/tmp/native_scenario_field.json")
}

fn canonical_poiseuille_snapshot() -> FieldSnapshot {
    let mut sim: Simulation<f64> = SimConfig {
        nx: N,
        ny: N,
        nu: NU,
        collision: Collision::default(),
        edges: Edges {
            left: EdgeBC::Periodic,
            right: EdgeBC::Periodic,
            bottom: EdgeBC::BounceBack,
            top: EdgeBC::BounceBack,
        },
        force: [G, 0.0],
    }
    .build()
    .unwrap();

    assert!(
        run_to_steady(&mut sim, STEADY_CHECK_EVERY, STEADY_TOL, STEADY_MAX_STEPS),
        "canonical T2 Poiseuille snapshot did not reach steady state; time = {}",
        sim.time()
    );

    FieldSnapshot {
        schema_version: 1,
        producer: "lbm-core native compat facade".to_string(),
        scenario_id: "t2_poiseuille_periodic_x_bounceback_y_n32".to_string(),
        validation_case: "VALIDATION.md T2 body-force Poiseuille".to_string(),
        field_order: "row-major y*nx+x, full domain including solid rim".to_string(),
        scalar: "f64".to_string(),
        lattice: "D2Q9".to_string(),
        collision: "TRT magic=3/16".to_string(),
        storage_note: "macroscopic fields reconstructed from native deviation-storage populations"
            .to_string(),
        nx: sim.nx(),
        ny: sim.ny(),
        steps: sim.time(),
        nu: NU,
        tau: 3.0 * NU + 0.5,
        force: [G, 0.0],
        rho: sim.rho_field().to_vec(),
        ux: sim.ux_field().to_vec(),
        uy: sim.uy_field().to_vec(),
        solid: sim.solid_field().to_vec(),
    }
}

fn assert_snapshot_shape(s: &FieldSnapshot) {
    assert_eq!(s.schema_version, 1);
    assert_eq!(s.nx, N);
    assert_eq!(s.ny, N);
    let cells = s.nx * s.ny;
    assert_eq!(s.rho.len(), cells);
    assert_eq!(s.ux.len(), cells);
    assert_eq!(s.uy.len(), cells);
    assert_eq!(s.solid.len(), cells);
    assert!(s.rho.iter().all(|v| v.is_finite()));
    assert!(s.ux.iter().all(|v| v.is_finite()));
    assert!(s.uy.iter().all(|v| v.is_finite()));
}

fn assert_float_roundtrip(label: &str, actual: &[f64], expected: &[f64]) {
    assert_eq!(actual.len(), expected.len(), "{label} length mismatch");
    let max_abs = actual
        .iter()
        .zip(expected)
        .map(|(a, e)| (a - e).abs())
        .fold(0.0f64, f64::max);
    assert!(
        max_abs <= 1.0e-18,
        "{label} JSON roundtrip max_abs={max_abs:e} exceeds 1e-18"
    );
}

#[test]
fn export_native_canonical_poiseuille_snapshot() {
    let snapshot = canonical_poiseuille_snapshot();
    assert_snapshot_shape(&snapshot);

    let path = workspace_export_path();
    fs::create_dir_all(path.parent().unwrap()).unwrap();
    fs::write(&path, serde_json::to_vec_pretty(&snapshot).unwrap()).unwrap();

    let roundtrip: FieldSnapshot =
        serde_json::from_slice(&fs::read(&path).expect("snapshot should be readable")).unwrap();
    assert_snapshot_shape(&roundtrip);
    assert_eq!(roundtrip.scenario_id, snapshot.scenario_id);
    assert_eq!(roundtrip.steps, snapshot.steps);
    assert_eq!(roundtrip.solid, snapshot.solid);
    assert_float_roundtrip("rho", &roundtrip.rho, &snapshot.rho);
    assert_float_roundtrip("ux", &roundtrip.ux, &snapshot.ux);
    assert_float_roundtrip("uy", &roundtrip.uy, &snapshot.uy);

    println!("exported {}", path.display());
}
