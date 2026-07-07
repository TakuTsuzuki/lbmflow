use lbm_scenario::{BioprocessScenario, Scenario, UnsupportedReason};
use serde_json::{json, Value};

fn base_scenario() -> Value {
    json!({
        "version": "bioprocess-1.0",
        "name": "m0-single-phase",
        "credibility_tier": "screening",
        "reactor": {
            "kind": "stirred_tank",
            "vessel_diameter_m": 0.2,
            "liquid_height_m": 0.2,
            "working_volume_m3": 0.006,
            "impellers": [{
                "kind": "rushton",
                "diameter_m": 0.07,
                "clearance_from_bottom_m": 0.07,
                "rotational_speed_rpm": 120.0,
                "blade_count": 6
            }],
            "baffles": [{
                "count": 4,
                "width_m": 0.02,
                "thickness_m": 0.002,
                "wall_attached": true
            }],
            "spargers": []
        },
        "fluids": {
            "liquid_density_kg_m3": 1000.0,
            "liquid_viscosity_pa_s": 0.001,
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
            "dt_s": 0.01,
            "backend": "cpu",
            "precision": "f64",
            "lattice": "d3q19"
        },
        "outputs": {
            "manifest_path": "out/m0-single-phase/manifest.json",
            "fields_every_n_steps": null,
            "probes_every_n_steps": null,
            "emit_qoi_json": true,
            "emit_qoi_csv": false
        }
    })
}

fn parse_value(value: Value) -> Result<BioprocessScenario, lbm_scenario::BioprocessScenarioError> {
    BioprocessScenario::from_json_str(&serde_json::to_string(&value).unwrap())
}

#[test]
fn parses_valid_stirred_tank_single_phase() {
    let scenario = parse_value(base_scenario()).expect("valid single-phase scenario should parse");
    assert_eq!(scenario.name, "m0-single-phase");
    assert!(scenario
        .physics
        .models
        .iter()
        .any(|model| { matches!(model, lbm_scenario::bioprocess::PhysicsModel::SinglePhase) }));
}

#[test]
fn parses_valid_aerated_stirred_tank_with_resolved_gas() {
    let mut scenario = base_scenario();
    scenario["reactor"]["spargers"] = json!([{
        "kind": "ring",
        "center_z_m": 0.015,
        "outer_radius_m": 0.035,
        "orifice_count": 12,
        "orifice_diameter_m": 0.001,
        "gas_volumetric_flow_m3_per_s": 1.0e-6,
        "vvm": null,
        "inlet_phase": "gas"
    }]);
    scenario["fluids"]["gas_density_kg_m3"] = json!(1.2);
    scenario["fluids"]["gas_viscosity_pa_s"] = json!(1.8e-5);
    scenario["fluids"]["surface_tension_n_m"] = json!(0.072);
    scenario["physics"] = json!([
        {
            "kind": "resolved_phase_field",
            "interface_width_m": 0.002,
            "mobility_m2_per_s": 1.0e-8,
            "contact_angle_deg": null
        },
        {
            "kind": "oxygen",
            "henry_constant": 1.0,
            "interfacial_flux_model": "henry_equilibrium",
            "our_model": "none"
        }
    ]);
    scenario["qoi"]["gas_holdup"] = json!({});
    scenario["qoi"]["kla"] = json!({});

    parse_value(scenario).expect("valid resolved-gas scenario should parse");
}

#[test]
fn rejects_kla_without_oxygen() {
    let mut scenario = base_scenario();
    scenario["qoi"]["kla"] = json!({});

    let err = parse_value(scenario).expect_err("kLa without oxygen should reject");
    assert_eq!(
        err.reason,
        UnsupportedReason::MissingDependency {
            depends_on: "oxygen".to_string()
        }
    );
}

#[test]
fn rejects_sparger_without_gas_model() {
    let mut scenario = base_scenario();
    scenario["reactor"]["spargers"] = json!([{
        "kind": "pipe",
        "center_z_m": 0.01,
        "length_m": 0.08,
        "diameter_m": 0.004,
        "gas_volumetric_flow_m3_per_s": null,
        "vvm": 0.05,
        "inlet_phase": "gas"
    }]);

    let err = parse_value(scenario).expect_err("sparger without gas model should reject");
    assert_eq!(
        err.reason,
        UnsupportedReason::MissingDependency {
            depends_on: "gas_model_in_physics".to_string()
        }
    );
}

#[test]
fn rejects_evidence_tier_without_dataset_reference() {
    let mut scenario = base_scenario();
    scenario["credibility_tier"] = json!("evidence");

    let err = parse_value(scenario).expect_err("evidence tier without registry should reject");
    assert_eq!(
        err.reason,
        UnsupportedReason::EvidenceGateFailed {
            missing: vec!["calibration_and_holdout_dataset_registry".to_string()]
        }
    );
}

#[test]
fn rejects_unknown_top_level_field() {
    let mut scenario = base_scenario();
    scenario["unexpected"] = json!(true);

    let err = serde_json::from_value::<BioprocessScenario>(scenario)
        .expect_err("unknown top-level field should reject");
    assert!(
        err.to_string().contains("unknown field"),
        "unexpected error: {err}"
    );
}

#[test]
fn rejects_missing_required_field() {
    let mut scenario = base_scenario();
    scenario.as_object_mut().unwrap().remove("version");

    let err = serde_json::from_value::<BioprocessScenario>(scenario)
        .expect_err("missing version should reject");
    assert!(
        err.to_string().contains("missing field `version`"),
        "unexpected error: {err}"
    );
}

#[test]
fn rejects_non_gas_sparger_inlet_phase() {
    let mut scenario = base_scenario();
    scenario["reactor"]["spargers"] = json!([{
        "kind": "point_orifices",
        "center_z_m": 0.01,
        "positions": [[0.0, 0.0, 0.01]],
        "gas_volumetric_flow_m3_per_s": 1.0e-6,
        "vvm": null,
        "inlet_phase": "liquid"
    }]);
    scenario["physics"] = json!([{
        "kind": "resolved_phase_field",
        "interface_width_m": 0.002,
        "mobility_m2_per_s": 1.0e-8,
        "contact_angle_deg": null
    }]);

    let err = parse_value(scenario).expect_err("non-gas sparger inlet should reject");
    assert_eq!(
        err.reason,
        UnsupportedReason::OutOfValidityRange {
            detail: "sparger inlet_phase must be gas".to_string()
        }
    );
}

#[test]
fn rejects_2d_lattice_for_bioprocess_scenario() {
    let mut scenario = base_scenario();
    scenario["run"]["lattice"] = json!("d2q9");

    let err = parse_value(scenario).expect_err("2D bioprocess scenario should reject");
    assert_eq!(err.reason, UnsupportedReason::NotImplemented);
}

#[test]
fn legacy_scenario_still_parses() {
    let (_, _, preset) = lbm_scenario::presets()
        .into_iter()
        .next()
        .expect("at least one legacy preset should exist");
    let json = serde_json::to_string(&preset).expect("legacy preset should serialize");
    serde_json::from_str::<Scenario>(&json).expect("legacy Scenario should still parse");
}
