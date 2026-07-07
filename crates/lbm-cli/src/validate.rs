use anyhow::{Context, Result};
use lbm_scenario::{BioprocessScenario, Scenario};
use serde_json::json;

pub fn run(path: &str, json_output: bool) -> Result<i32> {
    let text = if path == "-" {
        std::io::read_to_string(std::io::stdin())?
    } else {
        std::fs::read_to_string(path).with_context(|| format!("cannot read: {path}"))?
    };
    let value: serde_json::Value = serde_json::from_str(&text).with_context(|| {
        serde_json::to_string_pretty(&json!({
            "error": "invalid-scenario-json",
            "hint": "see `lbm schema` or `lbm schema --bioprocess`"
        }))
        .expect("error JSON must serialize")
    })?;

    if value.get("version").and_then(|v| v.as_str()) == Some("bioprocess-1.0") {
        return validate_bioprocess(&text, json_output);
    }
    validate_legacy(&text, json_output)
}

fn validate_bioprocess(text: &str, json_output: bool) -> Result<i32> {
    let scenario = match BioprocessScenario::from_json_str(text) {
        Ok(scenario) => scenario,
        Err(err) => {
            if json_output {
                println!(
                    "{}",
                    serde_json::to_string_pretty(&json!({
                        "ok": false,
                        "error": err,
                        "unit_report": null
                    }))?
                );
            } else {
                eprintln!("bioprocess scenario invalid: {}", err.message);
            }
            return Ok(1);
        }
    };
    let mut report = match scenario.unit_report_with_diagnostics() {
        Ok(report) => report,
        Err(err) => {
            if json_output {
                println!(
                    "{}",
                    serde_json::to_string_pretty(&json!({
                        "ok": false,
                        "error": err,
                        "unit_report": null
                    }))?
                );
            } else {
                eprintln!("unit feasibility failed: {}", err.message);
            }
            return Ok(1);
        }
    };
    if scenario.credibility_tier == lbm_scenario::bioprocess::CredibilityTier::Engineering {
        for warning in scenario.qoi.engineering_calibration_only_warnings() {
            report
                .feasibility
                .warnings
                .push(lbm_scenario::FeasibilityIssue {
                    code: "engineering_calibration_only".to_string(),
                    message: warning,
                    value: None,
                    threshold: None,
                });
        }
    }
    if scenario.has_stl_import() {
        if let Err(err) = crate::runner::prepare_bioprocess_geometry(&scenario) {
            if json_output {
                println!(
                    "{}",
                    serde_json::to_string_pretty(&json!({
                        "ok": false,
                        "error": err.to_string(),
                        "unit_report": report
                    }))?
                );
            } else {
                eprintln!("geometry preparation failed: {err}");
            }
            return Ok(1);
        }
    }
    let ok = report.feasibility.rejections.is_empty();
    if json_output {
        println!("{}", serde_json::to_string_pretty(&report)?);
    } else {
        print_bioprocess_human(&scenario, &report);
    }
    Ok(if ok { 0 } else { 1 })
}

fn validate_legacy(text: &str, json_output: bool) -> Result<i32> {
    let scenario: Scenario = serde_json::from_str(text)?;
    let warnings = lbm_scenario::validate(&scenario);
    let build_result = lbm_scenario::build_check(&scenario);
    let ok = build_result.is_ok();
    if json_output {
        println!(
            "{}",
            serde_json::to_string_pretty(&json!({
                "ok": ok,
                "error": build_result.err(),
                "warnings": warnings,
                "unit_report": null,
                "message": "legacy scenario; UnitReport unavailable"
            }))?
        );
    } else {
        println!("legacy scenario; UnitReport unavailable");
        for warning in &warnings {
            println!("warning[{}]: {}", warning.field, warning.message);
        }
        if let Err(error) = build_result {
            eprintln!("error: {error}");
        }
    }
    Ok(if ok { 0 } else { 1 })
}

fn print_bioprocess_human(scenario: &BioprocessScenario, report: &lbm_scenario::UnitReport) {
    println!("Bioprocess unit feasibility: {}", scenario.name);
    println!("{:<24} {:>16.8e}", "reynolds", report.groups.reynolds);
    println!("{:<24} {:>16.8e}", "froude", report.groups.froude);
    println!(
        "{:<24} {:>16.8e}",
        "mach_lattice", report.groups.mach_lattice
    );
    println!("{:<24} {:>16.8e}", "tau_lu", report.lattice.tau_lu);
    println!("{:<24} {:>16.8e}", "dx_m", report.lattice.dx_m);
    for warning in &report.feasibility.warnings {
        println!("warning[{}]: {}", warning.code, warning.message);
    }
    for rejection in &report.feasibility.rejections {
        println!("rejection[{}]: {}", rejection.code, rejection.message);
    }
}
