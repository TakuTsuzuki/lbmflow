//! Minimal MCP (Model Context Protocol) server over stdio.
//!
//! JSON-RPC 2.0, one message per line. Implements the subset needed for
//! tool serving: initialize / notifications/initialized / ping /
//! tools/list / tools/call.
//!
//! Register with e.g. `claude mcp add lbmflow -- lbm mcp` (or the codex
//! equivalent).
//!
//! Two execution models:
//! - `run_scenario` is synchronous and blocks the connection until the run
//!   finishes — fine for small runs.
//! - `start_run` spawns the run on a background thread and returns a
//!   deterministic run id (`run-<seq>-<scenario name>`) immediately. Poll
//!   `run_status { runId }` until `completed` (result carries the manifest)
//!   or `failed`; `list_runs` enumerates all runs of this server. At most
//!   [`MAX_CONCURRENT_RUNS`] runs execute at once — further `start_run`
//!   calls fail fast ("failed: too many concurrent runs") instead of
//!   saturating the CPU. Background runs live inside this process: keep the
//!   MCP connection open until `run_status` reports a terminal state.

use crate::runner;
use anyhow::Result;
use lbm_scenario::Scenario;
use serde_json::{json, Value};
use std::collections::HashMap;
use std::io::{BufRead, Write};
use std::path::PathBuf;
use std::sync::{Arc, Mutex, MutexGuard};
use std::time::Instant;

/// Upper bound on simultaneously running background jobs (safety rail so an
/// agent loop cannot fork-bomb the host CPU).
pub const MAX_CONCURRENT_RUNS: usize = 4;

// ---------------------------------------------------------------------------
// Background-run registry
// ---------------------------------------------------------------------------

/// State of one background run.
enum RunState {
    Running,
    Completed { manifest: Value },
    Failed { error: String },
}

impl RunState {
    fn label(&self) -> &'static str {
        match self {
            RunState::Running => "running",
            RunState::Completed { .. } => "completed",
            RunState::Failed { .. } => "failed",
        }
    }
}

struct RunEntry {
    seq: u64,
    scenario_name: String,
    out_dir: String,
    started: Instant,
    state: RunState,
}

#[derive(Default)]
struct RegistryState {
    next_seq: u64,
    runs: HashMap<String, RunEntry>,
}

/// `start_run` rejection: the concurrency cap is reached.
#[derive(Debug)]
struct TooManyRuns {
    running: usize,
    max: usize,
}

/// Registry of background runs, shared between the request loop and worker
/// threads (`Arc<Mutex<HashMap<runId, RunEntry>>>` underneath).
pub struct RunRegistry {
    max_concurrent: usize,
    state: Mutex<RegistryState>,
}

impl RunRegistry {
    fn new(max_concurrent: usize) -> Arc<Self> {
        Arc::new(Self {
            max_concurrent,
            state: Mutex::new(RegistryState::default()),
        })
    }

    /// Lock helper; a poisoned mutex only means a worker panicked between
    /// two consistent updates, so keep serving.
    fn lock(&self) -> MutexGuard<'_, RegistryState> {
        self.state.lock().unwrap_or_else(|p| p.into_inner())
    }

    /// Register a run and execute `job` on a background thread. Returns the
    /// deterministic run id `run-<seq>-<scenario name>` (seq starts at 1; no
    /// clock or randomness involved). Fails fast when `max_concurrent` jobs
    /// are already running — a rejected start consumes no sequence number.
    fn start<F>(
        self: &Arc<Self>,
        scenario_name: &str,
        out_dir: &str,
        job: F,
    ) -> std::result::Result<String, TooManyRuns>
    where
        F: FnOnce() -> Result<Value> + Send + 'static,
    {
        let run_id = {
            let mut st = self.lock();
            let running = st
                .runs
                .values()
                .filter(|e| matches!(e.state, RunState::Running))
                .count();
            if running >= self.max_concurrent {
                return Err(TooManyRuns {
                    running,
                    max: self.max_concurrent,
                });
            }
            st.next_seq += 1;
            let seq = st.next_seq;
            let run_id = format!("run-{seq}-{scenario_name}");
            st.runs.insert(
                run_id.clone(),
                RunEntry {
                    seq,
                    scenario_name: scenario_name.to_string(),
                    out_dir: out_dir.to_string(),
                    started: Instant::now(),
                    state: RunState::Running,
                },
            );
            run_id
        };
        let registry = Arc::clone(self);
        let id = run_id.clone();
        std::thread::spawn(move || {
            // A panicking run (or a panic escaping rayon inside lbm-core)
            // must land in `failed`, never kill the server silently.
            let outcome = std::panic::catch_unwind(std::panic::AssertUnwindSafe(job));
            let state = match outcome {
                Ok(Ok(manifest)) => RunState::Completed { manifest },
                Ok(Err(e)) => RunState::Failed {
                    error: e.to_string(),
                },
                Err(panic) => RunState::Failed {
                    // `.as_ref()`: inspect the payload, not the Box (a bare
                    // `&panic` would unsize-coerce the Box itself to dyn Any).
                    error: format!("panic: {}", panic_message(panic.as_ref())),
                },
            };
            if let Some(entry) = registry.lock().runs.get_mut(&id) {
                entry.state = state;
            }
        });
        Ok(run_id)
    }

    /// `run_status` payload for one run id, or `None` if unknown.
    fn status_json(&self, run_id: &str) -> Option<Value> {
        let st = self.lock();
        st.runs.get(run_id).map(|e| match &e.state {
            RunState::Running => json!({
                "runId": run_id,
                "state": "running",
                "elapsedNote": format!(
                    "{:.1} s since start. Wait a while and poll run_status again until completion.",
                    e.started.elapsed().as_secs_f64()
                ),
            }),
            RunState::Completed { manifest } => json!({
                "runId": run_id,
                "state": "completed",
                "outDir": e.out_dir,
                "manifest": manifest,
            }),
            RunState::Failed { error } => json!({
                "runId": run_id,
                "state": "failed",
                "error": error,
            }),
        })
    }

    /// `list_runs` payload: all runs in start order (seq ascending).
    fn list_json(&self) -> Value {
        let st = self.lock();
        let mut entries: Vec<&RunEntry> = st.runs.values().collect();
        entries.sort_by_key(|e| e.seq);
        Value::Array(
            entries
                .iter()
                .map(|e| {
                    json!({
                        "runId": format!("run-{}-{}", e.seq, e.scenario_name),
                        "state": e.state.label(),
                        "scenarioName": e.scenario_name,
                        "outDir": e.out_dir,
                    })
                })
                .collect(),
        )
    }
}

fn panic_message(panic: &(dyn std::any::Any + Send)) -> String {
    if let Some(s) = panic.downcast_ref::<&str>() {
        (*s).to_string()
    } else if let Some(s) = panic.downcast_ref::<String>() {
        s.clone()
    } else {
        "unknown panic".to_string()
    }
}

// ---------------------------------------------------------------------------
// Server loop
// ---------------------------------------------------------------------------

pub fn serve() -> Result<()> {
    let registry = RunRegistry::new(MAX_CONCURRENT_RUNS);
    let stdin = std::io::stdin();
    let mut stdout = std::io::stdout().lock();
    for line in stdin.lock().lines() {
        let line = line?;
        if line.trim().is_empty() {
            continue;
        }
        let msg: Value = match serde_json::from_str(&line) {
            Ok(v) => v,
            Err(_) => continue,
        };
        let id = msg.get("id").cloned();
        let method = msg.get("method").and_then(|m| m.as_str()).unwrap_or("");
        // Notifications (no id) need no reply.
        let Some(id) = id else { continue };
        let params = msg.get("params").cloned().unwrap_or(Value::Null);
        let reply = match method {
            "initialize" => json!({
                "protocolVersion": params.get("protocolVersion").and_then(|v| v.as_str()).unwrap_or("2024-11-05"),
                "capabilities": { "tools": {} },
                "serverInfo": { "name": "lbmflow", "version": env!("CARGO_PKG_VERSION") }
            }),
            "ping" => json!({}),
            "tools/list" => tools_list(),
            "tools/call" => match tools_call(&params, &registry) {
                Ok(v) => v,
                Err(e) => json!({
                    "content": [{ "type": "text", "text": format!("error: {e}") }],
                    "isError": true
                }),
            },
            _ => {
                write_msg(
                    &mut stdout,
                    &json!({ "jsonrpc": "2.0", "id": id,
                        "error": { "code": -32601, "message": format!("method not found: {method}") } }),
                )?;
                continue;
            }
        };
        write_msg(
            &mut stdout,
            &json!({ "jsonrpc": "2.0", "id": id, "result": reply }),
        )?;
    }
    Ok(())
}

fn write_msg(out: &mut impl Write, v: &Value) -> Result<()> {
    let s = serde_json::to_string(v)?;
    writeln!(out, "{s}")?;
    out.flush()?;
    Ok(())
}

fn scenario_schema() -> Value {
    json!({
        "type": "object",
        "description": "LBMFlow scenario (v0). Use the get_schema tool for the full format reference.",
        "required": ["name", "grid", "physics", "edges", "run"]
    })
}

fn tools_list() -> Value {
    json!({ "tools": [
        {
            "name": "run_scenario",
            "description": "Run an LBM fluid-simulation scenario synchronously and return the manifest (diagnostics, output file list). Outputs are written as PNG/CSV to outDir (default out/<name>). Consider checking the format first with validate_scenario or get_schema. Blocks until completion, so for long runs (tens of thousands of steps, large grids) or parallel sweeps use start_run instead.",
            "inputSchema": {
                "type": "object",
                "properties": {
                    "scenario": scenario_schema(),
                    "outDir": { "type": "string", "description": "output directory (default out/<name>)" }
                },
                "required": ["scenario"]
            }
        },
        {
            "name": "start_run",
            "description": "Run a scenario asynchronously in the background and return a runId immediately (non-blocking). Use this for long runs, parameter sweeps, and optimization loops, then poll run_status { runId } until completion. At most 4 concurrent runs (excess requests are rejected immediately with 'failed: too many concurrent runs'). Runs live inside this server process, so keep the MCP connection open until completion is confirmed.",
            "inputSchema": {
                "type": "object",
                "properties": {
                    "scenario": scenario_schema(),
                    "outDir": { "type": "string", "description": "output directory (default out/<name>)" }
                },
                "required": ["scenario"]
            }
        },
        {
            "name": "run_status",
            "description": "Return the state of a run started with start_run. state is running / completed / failed. If completed, includes the manifest (diagnostics, output file list). If failed, error carries the cause. While running, re-poll every few seconds.",
            "inputSchema": {
                "type": "object",
                "properties": {
                    "runId": { "type": "string", "description": "runId returned by start_run (run-<seq>-<scenario name>)" }
                },
                "required": ["runId"]
            }
        },
        {
            "name": "list_runs",
            "description": "List all runs started via start_run on this server (runId, state, scenario name, output dir) in start order. Useful for tracking sweep progress.",
            "inputSchema": { "type": "object", "properties": {} }
        },
        {
            "name": "validate_scenario",
            "description": "Validate a scenario without running it; returns configuration errors and stability warnings (tau, Mach number, grid Reynolds number).",
            "inputSchema": {
                "type": "object",
                "properties": { "scenario": scenario_schema() },
                "required": ["scenario"]
            }
        },
        {
            "name": "list_presets",
            "description": "Return the names and full scenario JSON of the built-in presets (lid-driven cavity, Kármán vortex street, two-phase droplet, droplet on wall). Useful as working examples for writing scenarios.",
            "inputSchema": { "type": "object", "properties": {} }
        },
        {
            "name": "get_schema",
            "description": "Return the complete format reference for scenario JSON (v0).",
            "inputSchema": { "type": "object", "properties": {} }
        }
    ]})
}

fn text_result(text: String) -> Value {
    json!({ "content": [{ "type": "text", "text": text }] })
}

/// Extract `{ scenario, outDir? }` shared by run_scenario / start_run.
fn scenario_args(args: &Value) -> Result<(Scenario, PathBuf)> {
    let sc: Scenario = serde_json::from_value(
        args.get("scenario")
            .cloned()
            .ok_or_else(|| anyhow::anyhow!("missing scenario"))?,
    )?;
    let out_dir = args
        .get("outDir")
        .and_then(|v| v.as_str())
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from("out").join(&sc.name));
    Ok((sc, out_dir))
}

fn tools_call(params: &Value, registry: &Arc<RunRegistry>) -> Result<Value> {
    let name = params
        .get("name")
        .and_then(|v| v.as_str())
        .unwrap_or_default();
    let args = params.get("arguments").cloned().unwrap_or(json!({}));
    match name {
        "run_scenario" => {
            let (sc, out_dir) = scenario_args(&args)?;
            let manifest = runner::run(&sc, &out_dir)?;
            Ok(text_result(serde_json::to_string_pretty(&json!({
                "manifest": manifest,
                "outDir": out_dir.display().to_string(),
            }))?))
        }
        "start_run" => {
            let (sc, out_dir) = scenario_args(&args)?;
            let out_dir_str = out_dir.display().to_string();
            let scenario_name = sc.name.clone();
            let job = move || -> Result<Value> {
                let manifest = runner::run(&sc, &out_dir)?;
                Ok(serde_json::to_value(&manifest)?)
            };
            match registry.start(&scenario_name, &out_dir_str, job) {
                Ok(run_id) => Ok(text_result(serde_json::to_string_pretty(&json!({
                    "runId": run_id,
                    "outDir": out_dir_str,
                    "note": "Running in the background. Poll run_status { runId } until it reports completed/failed.",
                }))?)),
                Err(TooManyRuns { running, max }) => Ok(json!({
                    "content": [{ "type": "text", "text": format!(
                        "failed: too many concurrent runs ({running} running / limit {max}). Wait for existing runs to finish via run_status, or check list_runs, then retry."
                    ) }],
                    "isError": true
                })),
            }
        }
        "run_status" => {
            let run_id = args
                .get("runId")
                .and_then(|v| v.as_str())
                .ok_or_else(|| anyhow::anyhow!("missing runId"))?;
            let status = registry.status_json(run_id).ok_or_else(|| {
                anyhow::anyhow!(
                    "runId '{run_id}' not found. Check list_runs for the run list"
                )
            })?;
            Ok(text_result(serde_json::to_string_pretty(&status)?))
        }
        "list_runs" => Ok(text_result(serde_json::to_string_pretty(
            &registry.list_json(),
        )?)),
        "validate_scenario" => {
            let sc: Scenario = serde_json::from_value(
                args.get("scenario")
                    .cloned()
                    .ok_or_else(|| anyhow::anyhow!("missing scenario"))?,
            )?;
            let warnings = lbm_scenario::validate(&sc);
            let build = lbm_scenario::build_check(&sc);
            Ok(text_result(serde_json::to_string_pretty(&json!({
                "ok": build.is_ok(),
                "error": build.err(),
                "warnings": warnings,
            }))?))
        }
        "list_presets" => {
            let list: Vec<Value> = lbm_scenario::presets()
                .into_iter()
                .map(
                    |(name, desc, sc)| json!({ "name": name, "description": desc, "scenario": sc }),
                )
                .collect();
            Ok(text_result(serde_json::to_string_pretty(&Value::Array(
                list,
            ))?))
        }
        "get_schema" => Ok(text_result(crate::SCHEMA_DOC.to_string())),
        other => anyhow::bail!("unknown tool: {other}"),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::mpsc;
    use std::time::Duration;

    /// Poll until the run reaches `want` (registry updates are asynchronous).
    fn wait_for_state(reg: &Arc<RunRegistry>, run_id: &str, want: &str) -> Value {
        for _ in 0..1000 {
            let v = reg.status_json(run_id).expect("run should be registered");
            if v["state"] == want {
                return v;
            }
            std::thread::sleep(Duration::from_millis(5));
        }
        panic!("timeout: {run_id} never reached state {want}");
    }

    /// Job that blocks until the returned sender is dropped or signalled,
    /// then finishes with the given result. Lets tests hold a run in
    /// `running` deterministically.
    fn gated_job(
        result: Result<Value>,
    ) -> (
        mpsc::Sender<()>,
        impl FnOnce() -> Result<Value> + Send + 'static,
    ) {
        let (tx, rx) = mpsc::channel::<()>();
        (tx, move || {
            let _ = rx.recv(); // released on send() or sender drop
            result
        })
    }

    #[test]
    fn run_ids_are_deterministic_and_sequential() {
        let reg = RunRegistry::new(4);
        let id1 = reg.start("cavity", "out/cavity", || Ok(json!({}))).unwrap();
        let id2 = reg
            .start("cavity", "out/cavity2", || Ok(json!({})))
            .unwrap();
        let id3 = reg.start("karman", "out/karman", || Ok(json!({}))).unwrap();
        assert_eq!(id1, "run-1-cavity");
        assert_eq!(id2, "run-2-cavity");
        assert_eq!(id3, "run-3-karman");
        // A fresh registry restarts the sequence: ids depend only on call
        // order, never on clock or randomness.
        let reg2 = RunRegistry::new(4);
        let again = reg2
            .start("cavity", "out/cavity", || Ok(json!({})))
            .unwrap();
        assert_eq!(again, "run-1-cavity");
    }

    #[test]
    fn transitions_running_to_completed_with_manifest() {
        let reg = RunRegistry::new(4);
        let (gate, job) = gated_job(Ok(json!({ "status": "completed", "stepsRun": 42 })));
        let id = reg.start("demo", "out/demo", job).unwrap();
        let st = reg.status_json(&id).unwrap();
        assert_eq!(st["state"], "running");
        assert!(
            st["elapsedNote"].as_str().is_some(),
            "running status should carry elapsedNote"
        );
        assert!(st.get("manifest").is_none());
        gate.send(()).unwrap();
        let done = wait_for_state(&reg, &id, "completed");
        assert_eq!(done["manifest"]["stepsRun"], 42);
        assert_eq!(done["outDir"], "out/demo");
    }

    #[test]
    fn transitions_running_to_failed_on_error() {
        let reg = RunRegistry::new(4);
        let id = reg
            .start("bad", "out/bad", || anyhow::bail!("boom: invalid config"))
            .unwrap();
        let failed = wait_for_state(&reg, &id, "failed");
        assert!(
            failed["error"].as_str().unwrap().contains("boom"),
            "error should carry the job message: {failed}"
        );
        assert!(failed.get("manifest").is_none());
    }

    #[test]
    fn panic_is_caught_as_failed() {
        let reg = RunRegistry::new(4);
        let id = reg
            .start("crashy", "out/crashy", || panic!("kaboom in solver"))
            .unwrap();
        let failed = wait_for_state(&reg, &id, "failed");
        let msg = failed["error"].as_str().unwrap();
        assert!(
            msg.contains("panic") && msg.contains("kaboom"),
            "panic payload should be reported: {msg}"
        );
    }

    #[test]
    fn concurrency_cap_rejects_then_recovers() {
        let reg = RunRegistry::new(4);
        let mut gates = Vec::new();
        for i in 0..4 {
            let (gate, job) = gated_job(Ok(json!({})));
            reg.start(&format!("s{i}"), "out", job).unwrap();
            gates.push(gate);
        }
        // 5th start must be rejected immediately (fail-fast safety rail).
        let (_gate5, job5) = gated_job(Ok(json!({})));
        let err = reg.start("s4", "out", job5).unwrap_err();
        assert_eq!((err.running, err.max), (4, 4));
        // Release one slot; the next start succeeds and — because the
        // rejected attempt consumed no sequence number — gets seq 5.
        gates.remove(0).send(()).unwrap();
        wait_for_state(&reg, "run-1-s0", "completed");
        let (gate6, job6) = gated_job(Ok(json!({})));
        let id = reg.start("s5", "out", job6).unwrap();
        assert_eq!(id, "run-5-s5");
        drop(gate6);
        drop(gates);
    }

    #[test]
    fn list_runs_is_sorted_by_start_order() {
        let reg = RunRegistry::new(4);
        let (gate, job) = gated_job(Ok(json!({})));
        reg.start("alpha", "out/a", job).unwrap();
        let id2 = reg.start("beta", "out/b", || Ok(json!({}))).unwrap();
        wait_for_state(&reg, &id2, "completed");
        let list = reg.list_json();
        let arr = list.as_array().unwrap();
        assert_eq!(arr.len(), 2);
        assert_eq!(arr[0]["runId"], "run-1-alpha");
        assert_eq!(arr[0]["state"], "running");
        assert_eq!(arr[0]["scenarioName"], "alpha");
        assert_eq!(arr[0]["outDir"], "out/a");
        assert_eq!(arr[1]["runId"], "run-2-beta");
        assert_eq!(arr[1]["state"], "completed");
        gate.send(()).unwrap();
        wait_for_state(&reg, "run-1-alpha", "completed");
    }

    #[test]
    fn unknown_run_id_yields_none() {
        let reg = RunRegistry::new(4);
        assert!(reg.status_json("run-99-ghost").is_none());
    }
}
