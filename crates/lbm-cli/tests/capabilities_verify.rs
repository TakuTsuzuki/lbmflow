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
fn verify_quick_passes_on_cpu() {
    let output = Command::new(env!("CARGO_BIN_EXE_lbm"))
        .args(["verify", "--tier", "quick"])
        .output()
        .expect("run lbm verify --tier quick");
    assert!(
        output.status.success(),
        "lbm verify --tier quick failed: status={} stdout={} stderr={}",
        output.status,
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(stdout.contains("PASS: T1 Taylor-Green vortex decay (CPU)"));
    assert!(stdout.contains("PASS: T2 body-force Poiseuille exactness (CPU)"));
}
