use crate::capabilities::{CapabilityRegistry, CapabilityStatus};
use anyhow::Result;
use clap::ValueEnum;
use serde_json::{json, Value};
use std::process::Command;

#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
pub enum VerifyTier {
    Quick,
    Bioprocess,
    Full,
}

impl VerifyTier {
    fn as_str(self) -> &'static str {
        match self {
            VerifyTier::Quick => "quick",
            VerifyTier::Bioprocess => "bioprocess",
            VerifyTier::Full => "full",
        }
    }
}

pub fn run(tier: VerifyTier, json_output: bool) -> Result<i32> {
    let (report, code) = report(tier);
    if json_output {
        println!("{}", serde_json::to_string_pretty(&report)?);
    } else {
        print_human(&report);
    }
    Ok(code)
}

pub fn report(tier: VerifyTier) -> (Value, i32) {
    let unsupported = unsupported_capabilities();
    let mut tests_run = 0u32;
    let mut tests_skipped = 0u32;
    let mut failure = Value::Null;
    let mut exit_code = 0;

    match tier {
        VerifyTier::Quick => match run_quick_checks() {
            Ok(n) => tests_run = n,
            Err(detail) => {
                tests_run = 1;
                failure = detail;
                exit_code = 1;
            }
        },
        VerifyTier::Bioprocess => {}
        VerifyTier::Full => {
            let result = run_full_cargo_test();
            tests_run = result.tests_run;
            tests_skipped = result.tests_skipped;
            if !result.success {
                failure = result.failure;
                exit_code = 1;
            }
        }
    }

    let value = json!({
        "tier": tier.as_str(),
        "tests_run": tests_run,
        "tests_skipped": tests_skipped,
        "unsupported_capabilities": unsupported,
        "validation_tier": "screening",
        "git_sha": git_sha(),
        "build_features": build_features(),
        "failure": failure,
    });
    (value, exit_code)
}

fn run_quick_checks() -> Result<u32, Value> {
    let registry = CapabilityRegistry::new();
    let ids = unsupported_capabilities();
    if registry.iter().count() != 9 || ids.len() != 5 {
        return Err(json!({
            "check": "capability_registry",
            "message": "capability registry must contain 9 BCFD entries and 5 unsupported entries",
            "capabilities": ids,
        }));
    }

    let cavity = lbm_scenario::presets()
        .into_iter()
        .find(|(name, _, _)| *name == "cavity")
        .map(|(_, _, scenario)| scenario)
        .ok_or_else(|| {
            json!({
                "check": "legacy_cavity_preset",
                "message": "cavity preset is missing"
            })
        })?;
    if let Err(message) = lbm_scenario::build_check(&cavity) {
        return Err(json!({
            "check": "legacy_cavity_build_check",
            "message": message
        }));
    }

    Ok(2)
}

fn unsupported_capabilities() -> Vec<&'static str> {
    CapabilityRegistry::new()
        .iter()
        .filter(|entry| entry.status == CapabilityStatus::Unsupported)
        .map(|entry| entry.id)
        .collect()
}

fn build_features() -> Vec<&'static str> {
    let mut features = vec!["default"];
    if cfg!(feature = "gpu") {
        features.push("gpu");
    }
    if cfg!(feature = "mpi") {
        features.push("mpi");
    }
    features
}

fn git_sha() -> Option<String> {
    let output = Command::new("git")
        .args(["rev-parse", "HEAD"])
        .output()
        .ok()?;
    if !output.status.success() {
        return None;
    }
    let sha = String::from_utf8_lossy(&output.stdout).trim().to_string();
    if sha.is_empty() {
        None
    } else {
        Some(sha)
    }
}

struct FullResult {
    success: bool,
    tests_run: u32,
    tests_skipped: u32,
    failure: Value,
}

fn run_full_cargo_test() -> FullResult {
    let command = [
        "cargo",
        "test",
        "--workspace",
        "--release",
        "--no-fail-fast",
    ];
    let output = Command::new(command[0]).args(&command[1..]).output();
    match output {
        Ok(output) => {
            let stdout = String::from_utf8_lossy(&output.stdout);
            let stderr = String::from_utf8_lossy(&output.stderr);
            let (tests_run, tests_skipped) = parse_cargo_test_counts(&stdout, &stderr);
            FullResult {
                success: output.status.success(),
                tests_run,
                tests_skipped,
                failure: if output.status.success() {
                    Value::Null
                } else {
                    json!({
                        "command": command.join(" "),
                        "status": output.status.code(),
                        "stdout_tail": tail_lines(&stdout, 80),
                        "stderr_tail": tail_lines(&stderr, 80),
                    })
                },
            }
        }
        Err(e) => FullResult {
            success: false,
            tests_run: 0,
            tests_skipped: 0,
            failure: json!({
                "command": command.join(" "),
                "message": e.to_string(),
            }),
        },
    }
}

fn parse_cargo_test_counts(stdout: &str, stderr: &str) -> (u32, u32) {
    let mut tests_run = 0u32;
    let mut tests_skipped = 0u32;
    for line in stdout.lines().chain(stderr.lines()) {
        let Some(rest) = line.trim().strip_prefix("test result:") else {
            continue;
        };
        for part in rest.split(';') {
            let words: Vec<&str> = part
                .trim()
                .trim_end_matches('.')
                .split_whitespace()
                .collect();
            for pair in words.windows(2) {
                let Ok(n) = pair[0].parse::<u32>() else {
                    continue;
                };
                match pair[1] {
                    "passed" | "failed" => tests_run = tests_run.saturating_add(n),
                    "ignored" => tests_skipped = tests_skipped.saturating_add(n),
                    _ => {}
                }
            }
        }
    }
    (tests_run, tests_skipped)
}

fn tail_lines(text: &str, max_lines: usize) -> String {
    let lines: Vec<&str> = text.lines().collect();
    let start = lines.len().saturating_sub(max_lines);
    lines[start..].join("\n")
}

fn print_human(report: &Value) {
    println!(
        "lbm verify tier={} tests_run={} tests_skipped={} validation_tier={}",
        report["tier"].as_str().unwrap_or("unknown"),
        report["tests_run"].as_u64().unwrap_or(0),
        report["tests_skipped"].as_u64().unwrap_or(0),
        report["validation_tier"].as_str().unwrap_or("screening")
    );
    if let Some(ids) = report["unsupported_capabilities"].as_array() {
        let ids: Vec<&str> = ids.iter().filter_map(|id| id.as_str()).collect();
        println!("unsupported_capabilities={}", ids.join(","));
    }
    if !report["failure"].is_null() {
        println!("failure={}", report["failure"]);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bioprocess_report_shape_has_empty_test_subset() {
        let (report, code) = report(VerifyTier::Bioprocess);
        assert_eq!(code, 0);
        assert_eq!(report["tier"], "bioprocess");
        assert_eq!(report["tests_run"], 0);
        assert_eq!(
            report["unsupported_capabilities"]
                .as_array()
                .expect("unsupported_capabilities should be an array")
                .len(),
            5
        );
        assert_eq!(report["validation_tier"], "screening");
        assert!(report["build_features"]
            .as_array()
            .expect("build_features should be an array")
            .iter()
            .any(|feature| feature == "default"));
    }

    #[test]
    fn parses_cargo_test_result_counts() {
        let stdout = "test result: ok. 3 passed; 1 failed; 2 ignored; 0 measured; 0 filtered out;";
        assert_eq!(parse_cargo_test_counts(stdout, ""), (4, 2));
    }
}
