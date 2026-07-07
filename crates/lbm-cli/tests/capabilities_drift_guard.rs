use lbm_cli::capabilities::{CapabilityRegistry, CapabilityStatus};
use std::fs;
use std::path::PathBuf;

fn repo_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..")
}

fn read_repo_file(path: &str) -> String {
    fs::read_to_string(repo_root().join(path)).unwrap_or_else(|e| panic!("read {path}: {e}"))
}

fn lower(text: &str) -> String {
    text.to_ascii_lowercase()
}

fn limitation_alias(id: &str, label: &str) -> String {
    match id {
        "rotating_ibm" => "rotating ibm".to_string(),
        "oxygen_kla" => "oxygen transport / kla".to_string(),
        "pbm" => "pbm".to_string(),
        "cell_exposure" => "cell / microcarrier exposure".to_string(),
        _ => lower(label),
    }
}

#[test]
fn limitations_lists_every_registry_capability_label() {
    let registry = CapabilityRegistry::new();
    let limitations = lower(&read_repo_file("docs/LIMITATIONS.md"));

    for entry in registry.iter() {
        let expected_status = match entry.id {
            "single_phase_stirred_tank" | "rotating_ibm" | "passive_scalar" => {
                CapabilityStatus::Experimental
            }
            _ => CapabilityStatus::Unsupported,
        };
        assert_eq!(
            entry.status, expected_status,
            "BCFD-002 registry status drift for {}",
            entry.id
        );
        let alias = limitation_alias(entry.id, entry.label);
        assert!(
            limitations.contains(&alias),
            "docs/LIMITATIONS.md is missing capability alias `{}` from registry id `{}`",
            alias,
            entry.id
        );
    }
}

#[test]
fn registry_docs_point_to_plan_bcfd_anchors() {
    let registry = CapabilityRegistry::new();
    let plan = lower(&read_repo_file("docs/PLAN.md"));

    for entry in registry.iter() {
        let anchor = entry
            .docs
            .strip_prefix("docs/PLAN.md#")
            .unwrap_or_else(|| panic!("{} docs must point into docs/PLAN.md", entry.id));
        assert!(
            plan.contains(anchor),
            "docs/PLAN.md is missing anchor target `{anchor}` for registry id `{}`",
            entry.id
        );
    }
}
