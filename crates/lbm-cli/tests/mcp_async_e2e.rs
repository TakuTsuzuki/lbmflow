//! E2E: validate the async job API of `lbm mcp` against a real process.
//!
//! Drives an actual JSON-RPC sequence (initialize → start_run → run_status polling →
//! completed confirmation → list_runs) over stdio, and also confirms compatibility
//! of the synchronous tools. In case the server hangs, a watchdog forcibly
//! kills the child process.

use serde_json::{json, Value};
use std::io::{BufRead, BufReader, Write};
use std::process::{Child, ChildStdin, ChildStdout, Command, Stdio};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::Arc;
use std::time::{Duration, Instant};

struct McpClient {
    child: Child,
    stdin: ChildStdin,
    stdout: BufReader<ChildStdout>,
    next_id: u64,
    test_done: Arc<AtomicBool>,
}

impl McpClient {
    fn spawn() -> Self {
        let mut child = Command::new(env!("CARGO_BIN_EXE_lbm"))
            .arg("mcp")
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .stderr(Stdio::null())
            .spawn()
            .expect("failed to spawn `lbm mcp`");
        let stdin = child.stdin.take().unwrap();
        let stdout = BufReader::new(child.stdout.take().unwrap());
        // Watchdog: kill the server if the test has not finished within the
        // deadline, so a hung server turns into EOF instead of a stuck CI job.
        let test_done = Arc::new(AtomicBool::new(false));
        let done = Arc::clone(&test_done);
        let pid = child.id();
        std::thread::spawn(move || {
            for _ in 0..120 {
                if done.load(Ordering::Relaxed) {
                    return;
                }
                std::thread::sleep(Duration::from_secs(1));
            }
            let _ = Command::new("kill").args(["-9", &pid.to_string()]).status();
        });
        Self {
            child,
            stdin,
            stdout,
            next_id: 0,
            test_done,
        }
    }

    /// Send one JSON-RPC request and read its response; returns `result`.
    fn request(&mut self, method: &str, params: Value) -> Value {
        self.next_id += 1;
        let id = self.next_id;
        let req = json!({ "jsonrpc": "2.0", "id": id, "method": method, "params": params });
        writeln!(self.stdin, "{req}").expect("write request");
        self.stdin.flush().unwrap();
        let mut line = String::new();
        let n = self.stdout.read_line(&mut line).expect("read response");
        assert!(
            n > 0,
            "server closed stdout while waiting for {method} response"
        );
        let msg: Value = serde_json::from_str(&line).expect("response should be JSON");
        assert_eq!(msg["id"], id, "response id should match request: {msg}");
        assert!(
            msg.get("error").is_none(),
            "unexpected JSON-RPC error for {method}: {msg}"
        );
        msg["result"].clone()
    }

    /// Fire-and-forget notification (no id, no response expected).
    fn notify(&mut self, method: &str) {
        let req = json!({ "jsonrpc": "2.0", "method": method });
        writeln!(self.stdin, "{req}").unwrap();
        self.stdin.flush().unwrap();
    }

    fn call_tool(&mut self, name: &str, arguments: Value) -> Value {
        self.request(
            "tools/call",
            json!({ "name": name, "arguments": arguments }),
        )
    }

    fn shutdown(mut self) {
        drop(self.stdin); // EOF -> server loop exits
        let status = self.child.wait().expect("wait for server exit");
        self.test_done.store(true, Ordering::Relaxed);
        assert!(status.success(), "server should exit cleanly: {status}");
    }
}

/// Text payload of a tool result.
fn tool_text(result: &Value) -> &str {
    result["content"][0]["text"]
        .as_str()
        .expect("tool result should carry text content")
}

/// Tool results that contain JSON text, parsed.
fn tool_json(result: &Value) -> Value {
    serde_json::from_str(tool_text(result))
        .unwrap_or_else(|e| panic!("tool text should be JSON ({e}): {result}"))
}

fn is_error(result: &Value) -> bool {
    result["isError"].as_bool().unwrap_or(false)
}

fn cavity_scenario(name: &str, nx: u64, steps: u64) -> Value {
    json!({
        "name": name,
        "grid": { "nx": nx, "ny": nx },
        "physics": { "nu": 0.05 },
        "edges": {
            "left":   { "type": "bounceBack" },
            "right":  { "type": "bounceBack" },
            "bottom": { "type": "bounceBack" },
            "top":    { "type": "movingWall", "u": [0.05, 0.0] }
        },
        "run": { "steps": steps },
        "outputs": [ { "field": "speed", "format": "png", "every": 0 } ]
    })
}

#[test]
fn mcp_async_job_lifecycle() {
    let tmp = std::env::temp_dir().join(format!("lbm_mcp_e2e_{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&tmp);
    let mut mcp = McpClient::spawn();

    // --- handshake ---------------------------------------------------------
    let init = mcp.request(
        "initialize",
        json!({ "protocolVersion": "2024-11-05", "capabilities": {} }),
    );
    assert_eq!(init["serverInfo"]["name"], "lbmflow");
    mcp.notify("notifications/initialized");

    // --- tools/list exposes the async trio alongside the legacy four -------
    let tools = mcp.request("tools/list", json!({}));
    let names: Vec<&str> = tools["tools"]
        .as_array()
        .unwrap()
        .iter()
        .map(|t| t["name"].as_str().unwrap())
        .collect();
    for expected in [
        "run_scenario",
        "start_run",
        "run_status",
        "list_runs",
        "validate_scenario",
        "list_presets",
        "get_schema",
    ] {
        assert!(
            names.contains(&expected),
            "missing tool {expected}: {names:?}"
        );
    }

    // --- start_run: small 64^2 / 500-step cavity, returns immediately ------
    let out_dir = tmp.join("cavity-out");
    let started = Instant::now();
    let res = mcp.call_tool(
        "start_run",
        json!({
            "scenario": cavity_scenario("mcp-e2e-cavity", 64, 500),
            "outDir": out_dir.display().to_string(),
        }),
    );
    assert!(!is_error(&res), "start_run should succeed: {res}");
    let ack = tool_json(&res);
    assert_eq!(ack["runId"], "run-1-mcp-e2e-cavity", "deterministic run id");
    assert_eq!(ack["outDir"], out_dir.display().to_string());
    assert!(
        started.elapsed() < Duration::from_secs(5),
        "start_run must not block until run completion"
    );

    // --- poll run_status until terminal ------------------------------------
    let deadline = Instant::now() + Duration::from_secs(60);
    let status = loop {
        let res = mcp.call_tool("run_status", json!({ "runId": "run-1-mcp-e2e-cavity" }));
        assert!(!is_error(&res), "run_status should succeed: {res}");
        let st = tool_json(&res);
        match st["state"].as_str().unwrap() {
            "running" => {
                assert!(
                    st["elapsedNote"].as_str().is_some(),
                    "running state should include elapsedNote: {st}"
                );
                assert!(Instant::now() < deadline, "run did not finish in time");
                std::thread::sleep(Duration::from_millis(50));
            }
            _ => break st,
        }
    };
    assert_eq!(
        status["state"], "completed",
        "run should complete: {status}"
    );
    let manifest = &status["manifest"];
    assert_eq!(manifest["scenario"], "mcp-e2e-cavity");
    assert_eq!(manifest["status"], "completed");
    assert_eq!(manifest["stepsRun"], 500);
    let files: Vec<&str> = manifest["files"]
        .as_array()
        .unwrap()
        .iter()
        .map(|f| f.as_str().unwrap())
        .collect();
    assert!(files.contains(&"speed_500.png"), "files: {files:?}");
    assert!(
        out_dir.join("manifest.json").is_file(),
        "manifest.json should be written to outDir"
    );

    // --- async failure path: broken scenario lands in `failed` -------------
    let res = mcp.call_tool(
        "start_run",
        json!({
            "scenario": {
                "name": "mcp-e2e-bad",
                "grid": { "nx": 16, "ny": 16 },
                "physics": { "nu": 0.05 },
                "edges": {
                    "left":   { "type": "periodic" },
                    "right":  { "type": "bounceBack" },
                    "bottom": { "type": "bounceBack" },
                    "top":    { "type": "bounceBack" }
                },
                "run": { "steps": 10 }
            },
            "outDir": tmp.join("bad-out").display().to_string(),
        }),
    );
    assert!(!is_error(&res));
    assert_eq!(tool_json(&res)["runId"], "run-2-mcp-e2e-bad");
    let deadline = Instant::now() + Duration::from_secs(30);
    let failed = loop {
        let st = tool_json(&mcp.call_tool("run_status", json!({ "runId": "run-2-mcp-e2e-bad" })));
        if st["state"] != "running" {
            break st;
        }
        assert!(
            Instant::now() < deadline,
            "failed run never became terminal"
        );
        std::thread::sleep(Duration::from_millis(20));
    };
    assert_eq!(failed["state"], "failed", "{failed}");
    assert!(
        !failed["error"].as_str().unwrap().is_empty(),
        "failed state should explain why: {failed}"
    );

    // --- run_status on an unknown id is a tool error ------------------------
    let res = mcp.call_tool("run_status", json!({ "runId": "run-99-ghost" }));
    assert!(is_error(&res), "unknown runId should be an error: {res}");

    // --- list_runs shows both runs in start order ---------------------------
    let list = tool_json(&mcp.call_tool("list_runs", json!({})));
    let runs = list.as_array().unwrap();
    assert_eq!(runs.len(), 2, "{list}");
    assert_eq!(runs[0]["runId"], "run-1-mcp-e2e-cavity");
    assert_eq!(runs[0]["state"], "completed");
    assert_eq!(runs[0]["scenarioName"], "mcp-e2e-cavity");
    assert_eq!(runs[1]["runId"], "run-2-mcp-e2e-bad");
    assert_eq!(runs[1]["state"], "failed");

    // --- legacy tools stay intact -------------------------------------------
    let res = mcp.call_tool(
        "run_scenario",
        json!({
            "scenario": cavity_scenario("mcp-e2e-sync", 16, 50),
            "outDir": tmp.join("sync-out").display().to_string(),
        }),
    );
    assert!(
        !is_error(&res),
        "sync run_scenario should still work: {res}"
    );
    assert_eq!(tool_json(&res)["manifest"]["stepsRun"], 50);

    let res = mcp.call_tool(
        "validate_scenario",
        json!({ "scenario": cavity_scenario("v", 16, 10) }),
    );
    assert_eq!(tool_json(&res)["ok"], true);

    let mut gpu_2d = cavity_scenario("gpu-2d", 16, 10);
    gpu_2d["compute"] = json!({ "backend": "gpu" });
    let res = mcp.call_tool("validate_scenario", json!({ "scenario": gpu_2d }));
    let report = tool_json(&res);
    assert_eq!(report["ok"], false, "{report}");
    let error = report["error"].as_str().unwrap();
    assert!(
        error.contains("requested backend \"gpu\" is unavailable")
            && error.contains("2D compat scenario path")
            && error.contains("--features gpu"),
        "{error}"
    );

    let res = mcp.call_tool("get_schema", json!({}));
    assert!(
        tool_text(&res).contains("start_run"),
        "schema doc should mention the async API"
    );

    let res = mcp.call_tool("list_presets", json!({}));
    assert!(!tool_json(&res).as_array().unwrap().is_empty());

    mcp.shutdown();
    let _ = std::fs::remove_dir_all(&tmp);
}
