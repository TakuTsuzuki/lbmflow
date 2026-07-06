use lbm_cli::capabilities::matrix;
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

fn mentions(text: &str, fact: &str) -> bool {
    let text = lower(text);
    let fact = lower(fact);
    text.contains(&fact) || text.contains(&fact.replace('_', "-"))
}

fn section_between<'a>(text: &'a str, start: &str, next_prefix: &str) -> &'a str {
    let start_idx = text
        .find(start)
        .unwrap_or_else(|| panic!("README.md is missing section marker {start:?}"));
    let after_start = &text[start_idx..];
    if let Some(next_idx) = after_start[start.len()..].find(next_prefix) {
        &after_start[..start.len() + next_idx]
    } else {
        after_start
    }
}

fn table_row_containing<'a>(section: &'a str, needle: &str) -> Option<&'a str> {
    let needle = lower(needle);
    section
        .lines()
        .find(|line| line.trim_start().starts_with('|') && lower(line).contains(&needle))
}

fn restriction_keywords(restriction: &str) -> Vec<&'static str> {
    let mut keywords = Vec::new();
    let lower = lower(restriction);
    for keyword in [
        "outflow",
        "convective",
        "rejected",
        "gpu",
        "open faces",
        "velocity-inlet",
        "pressure-outlet",
    ] {
        if mentions(&lower, keyword) {
            keywords.push(keyword);
        }
    }
    keywords
}

#[test]
fn readme_capability_matrix_matches_cli_capability_facts() {
    let matrix = matrix();
    let readme = read_repo_file("README.md");
    let capability_section = section_between(&readme, "## Capability matrix", "\n## ");

    for lattice in &matrix.lattices {
        let row = table_row_containing(capability_section, lattice.name).unwrap_or_else(|| {
            panic!(
                "README.md capability matrix diverged for lattice fact `{}`; fix is to update README.md or STATIC_FACTS/capability matrix together",
                lattice.name
            )
        });
        for restriction in &lattice.restrictions {
            for keyword in restriction_keywords(restriction) {
                assert!(
                    mentions(row, keyword),
                    "README.md capability matrix diverged for lattice `{}` restriction keyword `{}`; fix is to update README.md or STATIC_FACTS/capability matrix together",
                    lattice.name,
                    keyword
                );
            }
        }
    }

    let collision_row =
        table_row_containing(capability_section, "| Collision |").unwrap_or_else(|| {
            panic!(
                "README.md capability matrix diverged for collision row; fix is to update README.md or STATIC_FACTS/capability matrix together"
            )
        });
    for collision in &matrix.collisions.scenario_path {
        assert!(
            mentions(collision_row, collision),
            "README.md capability matrix diverged for scenario collision `{}`; fix is to update README.md or STATIC_FACTS/capability matrix together",
            collision
        );
    }
}

#[test]
fn limitations_matches_static_cli_capability_facts() {
    let matrix = matrix();
    let limitations = read_repo_file("docs/LIMITATIONS.md");

    assert!(
        mentions(&limitations, matrix.checkpoint.scope),
        "docs/LIMITATIONS.md diverged for checkpoint scope `{}`; fix is to update docs/LIMITATIONS.md or STATIC_FACTS together",
        matrix.checkpoint.scope
    );
    assert!(
        mentions(&limitations, matrix.particle_coupling),
        "docs/LIMITATIONS.md diverged for particle coupling `{}`; fix is to update docs/LIMITATIONS.md or STATIC_FACTS together",
        matrix.particle_coupling
    );
}
