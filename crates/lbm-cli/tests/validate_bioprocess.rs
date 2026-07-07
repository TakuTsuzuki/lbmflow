use serde_json::{json, Value};
use std::fs;
use std::process::Command;
use std::time::{SystemTime, UNIX_EPOCH};

fn lbm_command() -> Command {
    Command::new(env!("CARGO_BIN_EXE_lbm"))
}

fn write_json(prefix: &str, value: &Value) -> String {
    let nonce = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let path = std::env::temp_dir().join(format!("{prefix}-{nonce}.json"));
    fs::write(&path, serde_json::to_vec_pretty(value).unwrap()).unwrap();
    path.to_string_lossy().into_owned()
}

fn base_bioprocess() -> Value {
    json!({
        "version": "bioprocess-1.0",
        "name": "validate-bcfd-004",
        "credibility_tier": "screening",
        "reactor": {
            "kind": "stirred_tank",
            "vessel_diameter_m": 0.2,
            "liquid_height_m": 0.2,
            "working_volume_m3": 0.006,
            "impellers": [{
                "kind": "rushton",
                "diameter_m": 0.05,
                "clearance_from_bottom_m": 0.06,
                "rotational_speed_rpm": 60.0,
                "blade_count": 6
            }],
            "baffles": [],
            "spargers": []
        },
        "fluids": {
            "liquid_density_kg_m3": 1000.0,
            "liquid_viscosity_pa_s": 0.1,
            "gas_density_kg_m3": null,
            "gas_viscosity_pa_s": null,
            "surface_tension_n_m": null,
            "oxygen_diffusivity_m2_per_s": null,
            "henry_constant": null
        },
        "operation": {
            "duration_s": 30.0,
            "gas_inlet_temp_c": null,
            "initial_condition": { "kind": "quiescent" }
        },
        "physics": [{ "kind": "single_phase" }],
        "qoi": {
            "power": {},
            "mixing": null,
            "gas_holdup": null,
            "bubble_size": null,
            "kla": null,
            "shear_exposure": null,
            "oxygen_exposure": null,
            "calibration_dataset_id": null,
            "holdout_dataset_id": null
        },
        "run": {
            "steps": 100,
            "dt_s": 0.0008,
            "grid_nx": 100,
            "grid_ny": 100,
            "grid_nz": 100,
            "backend": "cpu",
            "precision": "f64",
            "lattice": "d3q19"
        },
        "outputs": {
            "manifest_path": "out/validate-bcfd-004/manifest.json",
            "fields_every_n_steps": null,
            "probes_every_n_steps": null,
            "emit_qoi_json": true,
            "emit_qoi_csv": false
        }
    })
}

#[test]
fn validate_command_emits_unit_report_json_for_bioprocess_scenario() {
    let path = write_json("bioprocess-ok", &base_bioprocess());
    let output = lbm_command()
        .args(["validate", &path, "--json"])
        .output()
        .expect("lbm validate should run");
    assert!(
        output.status.success(),
        "validate failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let stdout = String::from_utf8(output.stdout).unwrap();
    let report: Value = serde_json::from_str(&stdout).unwrap();
    assert_eq!(report["lattice"]["grid_nx"], 100);
    let reynolds = report["groups"]["reynolds"].as_f64().unwrap();
    assert!(
        (reynolds - 25.0).abs() < 1.0e-12,
        "unexpected Reynolds number: {reynolds:e}"
    );
    assert!(report["feasibility"]["rejections"]
        .as_array()
        .unwrap()
        .is_empty());
}

#[test]
fn validate_command_rejects_bioprocess_with_ma_over_cap() {
    let mut value = base_bioprocess();
    value["run"]["dt_s"] = json!(0.007159143792496811);
    let path = write_json("bioprocess-ma-reject", &value);
    let output = lbm_command()
        .args(["validate", &path, "--json"])
        .output()
        .expect("lbm validate should run");
    assert!(!output.status.success());
    let stdout = String::from_utf8(output.stdout).unwrap();
    let report: Value = serde_json::from_str(&stdout).unwrap();
    let rejections = report["feasibility"]["rejections"].as_array().unwrap();
    assert!(rejections
        .iter()
        .any(|issue| issue["code"] == "MA_LATTICE_TOO_HIGH"));
}

#[test]
fn validate_command_still_runs_for_legacy_scenario() {
    let (_, _, scenario) = lbm_scenario::presets()
        .into_iter()
        .next()
        .expect("legacy preset should exist");
    let path = write_json("legacy", &serde_json::to_value(scenario).unwrap());
    let output = lbm_command()
        .args(["validate", &path, "--json"])
        .output()
        .expect("lbm validate should run");
    assert!(
        output.status.success(),
        "legacy validate failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let stdout = String::from_utf8(output.stdout).unwrap();
    let report: Value = serde_json::from_str(&stdout).unwrap();
    assert_eq!(report["unit_report"], Value::Null);
    assert_eq!(report["message"], "legacy scenario; UnitReport unavailable");
}
