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
    for key in [
        "lattices",
        "collisions",
        "precisions",
        "backends",
        "backendGravityFallback",
        "checkpoint",
        "particleCoupling",
    ] {
        assert!(value.get(key).is_some(), "missing key {key}: {value}");
    }
    assert!(value["lattices"]
        .as_array()
        .unwrap()
        .iter()
        .any(|l| l["name"] == "d3q27"
            && l["restrictions"]
                .as_array()
                .unwrap()
                .iter()
                .any(|r| r.as_str().unwrap().contains("open faces"))));
    assert!(value["collisions"]["scenarioPath"]
        .as_array()
        .unwrap()
        .iter()
        .any(|c| c == "bgk"));
    assert!(value["backends"]
        .as_array()
        .unwrap()
        .iter()
        .any(|b| b["name"] == "cpu-scalar"
            && b["compiled"] == true
            && b["gravityBodyForce"] == true));
    assert!(value["backends"]
        .as_array()
        .unwrap()
        .iter()
        .any(|b| b["name"] == "cpu-simd"
            && b["compiled"] == true
            && b["gravityBodyForce"] == true));
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
