use lbm_cli::capabilities::{CapabilityRegistry, CapabilityStatus, CapabilityTier};
use std::collections::HashMap;
use std::fs;
use std::path::PathBuf;

#[derive(Clone, Debug, PartialEq, Eq)]
struct DocCapabilityRow {
    status: String,
    tier_ceiling: String,
}

fn repo_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..")
}

fn read_repo_file(path: &str) -> String {
    fs::read_to_string(repo_root().join(path)).unwrap_or_else(|e| panic!("read {path}: {e}"))
}

fn status_name(status: CapabilityStatus) -> &'static str {
    match status {
        CapabilityStatus::Unsupported => "Unsupported",
        CapabilityStatus::Experimental => "Experimental",
        CapabilityStatus::Engineering => "Engineering",
        CapabilityStatus::EvidenceBlocked => "EvidenceBlocked",
        CapabilityStatus::EvidenceReady => "EvidenceReady",
    }
}

fn tier_name(tier: CapabilityTier) -> &'static str {
    match tier {
        CapabilityTier::Screening => "Screening",
        CapabilityTier::Engineering => "Engineering",
        CapabilityTier::Evidence => "Evidence",
    }
}

fn strip_backticks(text: &str) -> &str {
    text.trim().trim_matches('`')
}

fn extract_capability_table(markdown: &str) -> HashMap<String, DocCapabilityRow> {
    let mut rows = HashMap::new();
    for line in markdown.lines() {
        let trimmed = line.trim();
        if !trimmed.starts_with("| `") {
            continue;
        }
        let cells: Vec<_> = trimmed
            .trim_matches('|')
            .split('|')
            .map(str::trim)
            .collect();
        if cells.len() < 4 {
            continue;
        }
        let id = strip_backticks(cells[0]).to_string();
        rows.insert(
            id,
            DocCapabilityRow {
                status: cells[2].to_string(),
                tier_ceiling: cells[3].to_string(),
            },
        );
    }
    rows
}

#[test]
fn limitations_registry_statuses_match_live_registry() {
    let registry = CapabilityRegistry::new();
    let rows = extract_capability_table(&read_repo_file("docs/LIMITATIONS.md"));

    assert_eq!(
        rows.len(),
        registry.iter().count(),
        "docs/LIMITATIONS.md capability table must contain exactly the registry ids: {rows:?}"
    );

    for entry in registry.iter() {
        let row = rows
            .get(entry.id)
            .unwrap_or_else(|| panic!("docs/LIMITATIONS.md is missing registry id `{}`", entry.id));
        assert_eq!(
            row.status,
            status_name(entry.status),
            "docs/LIMITATIONS.md status drift for registry id `{}`",
            entry.id
        );
        assert_eq!(
            row.tier_ceiling,
            tier_name(entry.tier_ceiling),
            "docs/LIMITATIONS.md tier-ceiling drift for registry id `{}`",
            entry.id
        );
    }
}

#[test]
fn readme_capability_matrix_mentions_every_registry_id_with_status() {
    let registry = CapabilityRegistry::new();
    let rows = extract_capability_table(&read_repo_file("README.md"));

    for entry in registry.iter() {
        let row = rows
            .get(entry.id)
            .unwrap_or_else(|| panic!("README.md capability matrix is missing `{}`", entry.id));
        assert_eq!(
            row.status,
            status_name(entry.status),
            "README.md status drift for registry id `{}`",
            entry.id
        );
        assert_eq!(
            row.tier_ceiling,
            tier_name(entry.tier_ceiling),
            "README.md tier-ceiling drift for registry id `{}`",
            entry.id
        );
    }
}

#[test]
fn registry_docs_point_to_plan_bcfd_anchors() {
    let registry = CapabilityRegistry::new();
    let plan = read_repo_file("docs/PLAN.md").to_ascii_lowercase();

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
