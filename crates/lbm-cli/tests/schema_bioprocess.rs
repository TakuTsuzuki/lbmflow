use std::process::Command;

fn lbm_command() -> Command {
    Command::new(env!("CARGO_BIN_EXE_lbm"))
}

#[test]
fn schema_bioprocess_flag_emits_bioprocess_schema() {
    let output = lbm_command()
        .args(["schema", "--bioprocess"])
        .output()
        .expect("lbm schema --bioprocess should run");

    assert!(
        output.status.success(),
        "schema command failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let stdout = String::from_utf8(output.stdout).expect("schema should be UTF-8");
    let schema: serde_json::Value =
        serde_json::from_str(&stdout).expect("bioprocess schema should be JSON");
    assert_eq!(schema["title"], "BioprocessScenario");
    assert_eq!(schema["properties"]["version"]["const"], "bioprocess-1.0");
    assert!(stdout.contains("credibility_tier"));
    assert!(stdout.contains("stirred_tank"));
}

#[test]
fn schema_default_still_emits_legacy_schema() {
    let output = lbm_command()
        .arg("schema")
        .output()
        .expect("lbm schema should run");

    assert!(
        output.status.success(),
        "schema command failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let stdout = String::from_utf8(output.stdout).expect("schema should be UTF-8");
    assert!(stdout.contains("Scenario JSON (v0)"));
    assert!(stdout.contains("\"grid\""));
    assert!(stdout.contains("\"physics\""));
    assert!(!stdout.contains("\"title\": \"BioprocessScenario\""));
}
