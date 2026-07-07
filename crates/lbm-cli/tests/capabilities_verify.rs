use serde_json::Value;
use std::process::Command;

#[test]
fn capabilities_json_is_machine_readable() {
    let output = Command::new(env!("CARGO_BIN_EXE_lbm"))
        .args(["capabilities", "--json"])
        .output()
        .expect("run lbm capabilities --json");
    assert!(
        output.status.success(),
        "lbm capabilities --json failed: status={} stderr={}",
        output.status,
        String::from_utf8_lossy(&output.stderr)
    );
    let value: Value = serde_json::from_slice(&output.stdout)
        .unwrap_or_else(|e| panic!("capabilities output must be JSON ({e})"));
    let capabilities = value["capabilities"]
        .as_array()
        .unwrap_or_else(|| panic!("capabilities must be an array: {value}"));
    let expected = [
        "single_phase_stirred_tank",
        "rotating_ibm",
        "passive_scalar",
        "phase_field_vof",
        "oxygen_kla",
        "point_bubbles",
        "pbm",
        "cell_exposure",
        "evidence_tier_report",
    ];
    assert_eq!(capabilities.len(), expected.len(), "{value}");
    for id in expected {
        let entry = capabilities
            .iter()
            .find(|entry| entry["id"] == id)
            .unwrap_or_else(|| panic!("missing capability id {id}: {value}"));
        assert_eq!(entry["status"], "unsupported", "{entry}");
        assert!(entry["docs"]
            .as_str()
            .unwrap()
            .starts_with("docs/PLAN.md#bcfd-"));
    }
}

#[test]
fn verify_quick_json_has_machine_readable_shape() {
    let output = Command::new(env!("CARGO_BIN_EXE_lbm"))
        .args(["verify", "--tier", "quick", "--json"])
        .output()
        .expect("run lbm verify --tier quick --json");
    assert!(
        output.status.success(),
        "lbm verify --tier quick --json failed: status={} stdout={} stderr={}",
        output.status,
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
    let value: Value = serde_json::from_slice(&output.stdout)
        .unwrap_or_else(|e| panic!("verify output must be JSON ({e})"));
    assert_eq!(value["tier"], "quick", "{value}");
    assert!(value["tests_run"].as_u64().unwrap() >= 1, "{value}");
    assert!(value["tests_skipped"].is_u64(), "{value}");
    assert_eq!(value["validation_tier"], "screening", "{value}");
    assert!(value["git_sha"].is_string() || value["git_sha"].is_null(), "{value}");
    assert!(value["build_features"]
        .as_array()
        .unwrap()
        .iter()
        .any(|feature| feature == "default"));
    assert_eq!(
        value["unsupported_capabilities"].as_array().unwrap().len(),
        9,
        "{value}"
    );
    assert!(value["failure"].is_null(), "{value}");
}
