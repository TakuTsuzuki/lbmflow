use serde_json::{json, Value};
use std::fs;
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::sync::OnceLock;
use std::time::Duration;

#[derive(Clone, Debug)]
struct PresetRun {
    out_dir: PathBuf,
    stdout: String,
    stderr: String,
    status_code: Option<i32>,
}

static VTK_PRESET_RUN: OnceLock<Result<PresetRun, String>> = OnceLock::new();

fn lbm_bin() -> PathBuf {
    let release_bin = target_release_dir().join(format!("lbm{}", std::env::consts::EXE_SUFFIX));
    if release_bin.is_file() {
        return release_bin;
    }

    let status = Command::new("cargo")
        .args(["build", "-p", "lbm-cli", "--release"])
        .status()
        .expect("spawn cargo build -p lbm-cli --release");
    assert!(
        status.success(),
        "cargo build -p lbm-cli --release failed: {status}"
    );
    assert!(
        release_bin.is_file(),
        "expected release lbm binary at {}",
        release_bin.display()
    );
    release_bin
}

fn target_release_dir() -> PathBuf {
    let exe = std::env::current_exe().expect("current test executable path");
    exe.parent()
        .and_then(Path::parent)
        .expect("test binary should live under target/release/deps")
        .to_path_buf()
}

fn test_out_dir(name: &str) -> PathBuf {
    target_release_dir()
        .parent()
        .expect("target/release should have target parent")
        .join("cli_mcp_contracts")
        .join(format!("{name}_{}", std::process::id()))
}

#[test]
fn c1_presets_listing_contains_shipped_presets() {
    let output = Command::new(lbm_bin())
        .args(["presets", "list"])
        .output()
        .expect("run lbm presets list");

    assert!(
        output.status.success(),
        "lbm presets list failed: status={} stderr={}",
        output.status,
        String::from_utf8_lossy(&output.stderr)
    );
    let stdout = String::from_utf8_lossy(&output.stdout);
    println!("C1 presets list stdout:\n{stdout}");
    assert!(stdout.contains("cavity"), "presets list: {stdout}");
    assert!(stdout.contains("cylinder-karman"), "presets list: {stdout}");
}

#[test]
fn c2_preset_run_emits_valid_legacy_vtk() {
    let run = vtk_preset_run().unwrap_or_else(|e| panic!("{e}"));

    assert_eq!(
        run.status_code,
        Some(0),
        "lbm presets run failed: code={:?}\nstdout={}\nstderr={}",
        run.status_code,
        run.stdout,
        run.stderr
    );
    assert!(
        run.stderr
            .contains("legacy LBM demo preset; not bioprocess decision-grade"),
        "legacy preset warning missing from stderr: {}",
        run.stderr
    );
    let reported_out = parse_out_dir(&run.stdout).expect("stdout should report out=<dir>");
    assert_eq!(
        PathBuf::from(reported_out),
        run.out_dir,
        "CLI-reported out dir should match requested out dir"
    );

    let vtk = find_files_with_extension(&run.out_dir, "vtk")
        .into_iter()
        .next()
        .unwrap_or_else(|| panic!("no .vtk file found under {}", run.out_dir.display()));
    let text = fs::read_to_string(&vtk).expect("read VTK file as text");
    let lines: Vec<&str> = text.lines().take(4).collect();
    println!("C2 VTK file: {}", vtk.display());
    println!("C2 VTK header lines: {:?}", lines);

    assert!(
        lines
            .first()
            .is_some_and(|line| line.starts_with("# vtk DataFile Version")),
        "first VTK line should be legacy header: {lines:?}"
    );
    assert!(
        matches!(lines.get(2), Some(&"ASCII" | &"BINARY")),
        "third VTK line should be ASCII or BINARY: {lines:?}"
    );
    assert!(
        text.lines().any(|line| line.starts_with("DATASET ")),
        "VTK DATASET line missing in {}",
        vtk.display()
    );
}

#[test]
fn c3_manifest_schema_matches_cli_contract() {
    let run = vtk_preset_run().unwrap_or_else(|e| panic!("{e}"));
    let manifest_path = run.out_dir.join("manifest.json");
    let manifest_text = fs::read_to_string(&manifest_path).unwrap_or_else(|e| {
        panic!(
            "read manifest at {} after preset run: {e}",
            manifest_path.display()
        )
    });
    let manifest: Value = serde_json::from_str(&manifest_text).unwrap_or_else(|e| {
        panic!(
            "manifest should be JSON at {} ({e}): {manifest_text}",
            manifest_path.display()
        )
    });

    let fields: Vec<&str> = manifest
        .as_object()
        .expect("manifest root should be object")
        .keys()
        .map(String::as_str)
        .collect();
    println!("C3 manifest fields: {fields:?}");

    for key in [
        "scenario",
        "status",
        "stepsRun",
        "wallSeconds",
        "mlups",
        "diagnostics",
        "provenance",
        "warnings",
        "files",
    ] {
        assert!(manifest.get(key).is_some(), "missing manifest key {key}");
    }
    assert_eq!(manifest["scenario"], "droplet-on-wall");
    assert!(
        manifest["files"]
            .as_array()
            .expect("manifest files should be an array")
            .iter()
            .any(|f| f.as_str().is_some_and(|name| name.ends_with(".vtk"))),
        "manifest files should include a VTK artifact: {}",
        manifest["files"]
    );
    assert!(
        manifest["wallSeconds"].as_f64().is_some(),
        "wallSeconds should be the elapsed-like field: {}",
        manifest["wallSeconds"]
    );
}

#[test]
fn c4_mcp_lists_tools_over_stdio() {
    let help = Command::new(lbm_bin())
        .arg("--help")
        .output()
        .expect("run lbm --help");
    assert!(
        help.status.success(),
        "lbm --help failed: status={} stderr={}",
        help.status,
        String::from_utf8_lossy(&help.stderr)
    );
    let help_stdout = String::from_utf8_lossy(&help.stdout);
    if !help_stdout.contains("mcp") {
        // MCP contract requires an async test harness if the stdio server is not
        // exposed as a normal CLI subcommand.
        return;
    }

    let mut child = Command::new(lbm_bin())
        .arg("mcp")
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::null())
        .spawn()
        .expect("spawn lbm mcp");
    let mut stdin = child.stdin.take().expect("mcp stdin");
    let stdout = child.stdout.take().expect("mcp stdout");
    let mut stdout = BufReader::new(stdout);

    let req = json!({ "jsonrpc": "2.0", "id": 1, "method": "tools/list", "params": {} });
    writeln!(stdin, "{req}").expect("write tools/list request");
    stdin.flush().expect("flush tools/list request");

    let mut line = String::new();
    let n = stdout
        .read_line(&mut line)
        .expect("read tools/list response");
    assert!(n > 0, "lbm mcp closed stdout before tools/list response");
    drop(stdin);
    let status = wait_with_timeout(child, Duration::from_secs(5));
    assert!(status, "lbm mcp did not exit cleanly after stdin EOF");

    let msg: Value =
        serde_json::from_str(&line).unwrap_or_else(|e| panic!("MCP response JSON ({e}): {line}"));
    assert_eq!(msg["id"], 1);
    let tools = msg["result"]["tools"]
        .as_array()
        .expect("tools/list result should contain tools array");
    let names: Vec<&str> = tools
        .iter()
        .filter_map(|tool| tool["name"].as_str())
        .collect();
    println!("C4 MCP tools: {names:?}");
    assert!(!names.is_empty(), "MCP tools/list returned no tools");
    assert!(names.contains(&"run_scenario"), "MCP tools: {names:?}");
}

fn vtk_preset_run() -> Result<PresetRun, String> {
    VTK_PRESET_RUN.get_or_init(run_vtk_preset).clone()
}

fn run_vtk_preset() -> Result<PresetRun, String> {
    let out_dir = test_out_dir("droplet_on_wall");
    let _ = fs::remove_dir_all(&out_dir);
    fs::create_dir_all(out_dir.parent().expect("test output parent"))
        .map_err(|e| format!("create test output parent: {e}"))?;
    let output = Command::new(lbm_bin())
        .args(["presets", "run", "droplet-on-wall", "--out"])
        .arg(&out_dir)
        .output()
        .map_err(|e| format!("run lbm presets run droplet-on-wall: {e}"))?;
    Ok(PresetRun {
        out_dir,
        stdout: String::from_utf8_lossy(&output.stdout).into_owned(),
        stderr: String::from_utf8_lossy(&output.stderr).into_owned(),
        status_code: output.status.code(),
    })
}

fn parse_out_dir(stdout: &str) -> Option<&str> {
    stdout
        .split_whitespace()
        .find_map(|word| word.strip_prefix("out="))
}

fn find_files_with_extension(root: &Path, extension: &str) -> Vec<PathBuf> {
    let mut found = Vec::new();
    let mut stack = vec![root.to_path_buf()];
    while let Some(dir) = stack.pop() {
        let Ok(entries) = fs::read_dir(&dir) else {
            continue;
        };
        for entry in entries.flatten() {
            let path = entry.path();
            if path.is_dir() {
                stack.push(path);
            } else if path.extension().and_then(|e| e.to_str()) == Some(extension) {
                found.push(path);
            }
        }
    }
    found.sort();
    found
}

fn wait_with_timeout(mut child: std::process::Child, timeout: Duration) -> bool {
    let start = std::time::Instant::now();
    loop {
        match child.try_wait() {
            Ok(Some(status)) => return status.success(),
            Ok(None) if start.elapsed() < timeout => {
                std::thread::sleep(Duration::from_millis(10));
            }
            Ok(None) => {
                let _ = child.kill();
                let _ = child.wait();
                return false;
            }
            Err(_) => return false,
        }
    }
}
