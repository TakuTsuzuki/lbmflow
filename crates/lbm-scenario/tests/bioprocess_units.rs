use lbm_scenario::units::{
    BIOPROCESS_GRID_REYNOLDS_WARN_THRESHOLD, MACH_LATTICE_REJECT_THRESHOLD,
    MACH_LATTICE_WARN_THRESHOLD, TAU_NEAR_HALF_WARN_THRESHOLD,
};
use lbm_scenario::{BioprocessScenario, UnsupportedReason};
use serde_json::{json, Value};

fn base_scenario() -> Value {
    json!({
        "version": "bioprocess-1.0",
        "name": "bcfd-004-units",
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
            "manifest_path": "out/bcfd-004-units/manifest.json",
            "fields_every_n_steps": null,
            "probes_every_n_steps": null,
            "emit_qoi_json": true,
            "emit_qoi_csv": false
        }
    })
}

fn parse(value: Value) -> BioprocessScenario {
    BioprocessScenario::from_json_str(&serde_json::to_string(&value).unwrap()).unwrap()
}

fn codes(report: &lbm_scenario::UnitReport, rejections: bool) -> Vec<&str> {
    let issues = if rejections {
        &report.feasibility.rejections
    } else {
        &report.feasibility.warnings
    };
    issues.iter().map(|issue| issue.code.as_str()).collect()
}

fn assert_close(actual: f64, expected: f64) {
    let scale = actual.abs().max(expected.abs()).max(1.0);
    assert!(
        (actual - expected).abs() <= 1.0e-12 * scale,
        "actual={actual:e} expected={expected:e}"
    );
}

#[test]
fn dimensionless_groups_reynolds_matches_analytical_impeller() {
    let scenario = parse(base_scenario());
    let report = scenario.compute_unit_report().unwrap();
    assert_close(report.groups.reynolds, 25.0);
}

#[test]
fn dimensionless_groups_optional_absent_when_no_surface_tension() {
    let mut value = base_scenario();
    value["fluids"]["gas_density_kg_m3"] = json!(1.2);
    value["physics"] = json!([{
        "kind": "resolved_phase_field",
        "interface_width_m": 0.01,
        "mobility_m2_per_s": 1.0e-5,
        "contact_angle_deg": null
    }]);
    let report = parse(value).unit_report_with_diagnostics().unwrap();
    assert_eq!(report.groups.weber, None);
    assert_eq!(report.groups.eotvos, None);
    assert_eq!(report.groups.morton, None);
    assert!(codes(&report, false).contains(&"WEBER_UNAVAILABLE"));
}

#[test]
fn mach_lattice_warn_at_0_11_reject_at_0_31() {
    let mut warn = base_scenario();
    warn["run"]["dt_s"] = json!(0.0025403411844343537);
    let warn_report = parse(warn).compute_unit_report().unwrap();
    assert!(warn_report.groups.mach_lattice > MACH_LATTICE_WARN_THRESHOLD);
    assert!(codes(&warn_report, false).contains(&"MA_LATTICE_HIGH"));

    let mut reject = base_scenario();
    reject["run"]["dt_s"] = json!(0.007159143792496811);
    let scenario = parse(reject);
    let raw = scenario.unit_report_with_diagnostics().unwrap();
    assert!(raw.groups.mach_lattice > MACH_LATTICE_REJECT_THRESHOLD);
    assert!(codes(&raw, true).contains(&"MA_LATTICE_TOO_HIGH"));
    let err = scenario.compute_unit_report().unwrap_err();
    assert!(matches!(
        err.reason,
        UnsupportedReason::OutOfValidityRange { .. }
    ));
}

#[test]
fn tau_near_half_warn_and_below_half_rejects() {
    let mut warn = base_scenario();
    warn["run"]["dt_s"] = json!(0.0002);
    let report = parse(warn).compute_unit_report().unwrap();
    assert!(report.lattice.tau_lu <= TAU_NEAR_HALF_WARN_THRESHOLD);
    assert!(codes(&report, false).contains(&"TAU_NEAR_HALF"));

    let mut bad = base_scenario();
    bad["fluids"]["liquid_viscosity_pa_s"] = json!(0.0);
    let err = parse(bad).compute_unit_report().unwrap_err();
    assert!(matches!(
        err.reason,
        UnsupportedReason::OutOfValidityRange { .. }
    ));
}

#[test]
fn cahn_thin_and_thick_reject_only_when_phase_field_on() {
    let mut single_phase = base_scenario();
    single_phase["physics"] = json!([{ "kind": "single_phase" }]);
    let report = parse(single_phase).compute_unit_report().unwrap();
    assert_eq!(report.groups.cahn, None);

    let mut thin = base_scenario();
    thin["physics"] = json!([{
        "kind": "resolved_phase_field",
        "interface_width_m": 0.001,
        "mobility_m2_per_s": 1.0e-5,
        "contact_angle_deg": null
    }]);
    assert!(
        codes(&parse(thin).unit_report_with_diagnostics().unwrap(), true)
            .contains(&"INTERFACE_WIDTH_TOO_THIN")
    );

    let mut thick = base_scenario();
    thick["physics"] = json!([{
        "kind": "resolved_phase_field",
        "interface_width_m": 0.12,
        "mobility_m2_per_s": 1.0e-5,
        "contact_angle_deg": null
    }]);
    assert!(
        codes(&parse(thick).unit_report_with_diagnostics().unwrap(), true)
            .contains(&"INTERFACE_WIDTH_TOO_LARGE")
    );
}

#[test]
fn bubble_under_resolved_rejected_for_resolved_gas_injection() {
    let mut value = base_scenario();
    value["reactor"]["spargers"] = json!([{
        "kind": "ring",
        "center_z_m": 0.02,
        "outer_radius_m": 0.03,
        "orifice_count": 8,
        "orifice_diameter_m": 0.001,
        "gas_volumetric_flow_m3_per_s": 1.0e-6,
        "vvm": null,
        "inlet_phase": "gas"
    }]);
    value["fluids"]["gas_density_kg_m3"] = json!(1.2);
    value["fluids"]["surface_tension_n_m"] = json!(0.072);
    value["physics"] = json!([{
        "kind": "resolved_phase_field",
        "interface_width_m": 0.01,
        "mobility_m2_per_s": 1.0e-5,
        "contact_angle_deg": null
    }]);
    let report = parse(value).unit_report_with_diagnostics().unwrap();
    assert!(codes(&report, true).contains(&"BUBBLE_UNDER_RESOLVED"));
}

#[test]
fn grid_reynolds_warning_fires_at_configured_threshold() {
    let mut value = base_scenario();
    value["fluids"]["liquid_viscosity_pa_s"] = json!(0.001);
    let report = parse(value).unit_report_with_diagnostics().unwrap();
    let grid_re = report.lattice.u_ref_lu / report.lattice.nu_lu;
    assert!(grid_re > BIOPROCESS_GRID_REYNOLDS_WARN_THRESHOLD);
    assert!(codes(&report, false).contains(&"GRID_REYNOLDS_HIGH"));
}

#[test]
fn matching_priority_declares_only_re_as_mandatory() {
    let report = parse(base_scenario()).compute_unit_report().unwrap();
    assert_eq!(report.matching_priority.achieved, vec!["Re"]);
    assert!(report
        .matching_priority
        .missing
        .iter()
        .any(|entry| entry.contains("Sc/Pe/Da")));
}

#[test]
fn unit_report_serialization_roundtrip() {
    let report = parse(base_scenario()).compute_unit_report().unwrap();
    let encoded = serde_json::to_string(&report).unwrap();
    assert!(!encoded.contains("NaN"));
    let decoded: lbm_scenario::UnitReport = serde_json::from_str(&encoded).unwrap();
    assert_eq!(decoded, report);
}
