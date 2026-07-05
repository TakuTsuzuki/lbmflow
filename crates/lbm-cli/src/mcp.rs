//! Minimal MCP (Model Context Protocol) server over stdio.
//!
//! JSON-RPC 2.0, one message per line. Implements the subset needed for
//! tool serving: initialize / notifications/initialized / ping /
//! tools/list / tools/call.
//!
//! Register with e.g. `claude mcp add lbmflow -- lbm mcp` (or the codex
//! equivalent). Long runs block the connection until finished — size
//! scenarios accordingly or launch several servers.

use crate::runner;
use anyhow::Result;
use lbm_scenario::Scenario;
use serde_json::{json, Value};
use std::io::{BufRead, Write};
use std::path::PathBuf;

pub fn serve() -> Result<()> {
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
            "tools/call" => match tools_call(&params) {
                Ok(v) => v,
                Err(e) => json!({
                    "content": [{ "type": "text", "text": format!("エラー: {e}") }],
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
        "description": "LBMFlow シナリオ (v0)。詳細な書式は get_schema ツールで取得。",
        "required": ["name", "grid", "physics", "edges", "run"]
    })
}

fn tools_list() -> Value {
    json!({ "tools": [
        {
            "name": "run_scenario",
            "description": "LBM 流体シミュレーションのシナリオを実行し、manifest（診断・出力ファイル一覧）を返す。出力は outDir（既定 out/<name>）に PNG/CSV で書き出される。まず validate_scenario か get_schema で書式を確認するとよい。",
            "inputSchema": {
                "type": "object",
                "properties": {
                    "scenario": scenario_schema(),
                    "outDir": { "type": "string", "description": "出力ディレクトリ（省略時 out/<name>）" }
                },
                "required": ["scenario"]
            }
        },
        {
            "name": "validate_scenario",
            "description": "シナリオを実行せずに検証し、構成エラーと安定性警告（tau・マッハ数・グリッドレイノルズ数）を返す。",
            "inputSchema": {
                "type": "object",
                "properties": { "scenario": scenario_schema() },
                "required": ["scenario"]
            }
        },
        {
            "name": "list_presets",
            "description": "組み込みプリセット（キャビティ流れ・カルマン渦列・二相液滴）の名前と完全なシナリオ JSON を返す。シナリオ作成の実例として使える。",
            "inputSchema": { "type": "object", "properties": {} }
        },
        {
            "name": "get_schema",
            "description": "シナリオ JSON (v0) の完全な書式リファレンスを返す。",
            "inputSchema": { "type": "object", "properties": {} }
        }
    ]})
}

fn text_result(text: String) -> Value {
    json!({ "content": [{ "type": "text", "text": text }] })
}

fn tools_call(params: &Value) -> Result<Value> {
    let name = params
        .get("name")
        .and_then(|v| v.as_str())
        .unwrap_or_default();
    let args = params.get("arguments").cloned().unwrap_or(json!({}));
    match name {
        "run_scenario" => {
            let sc: Scenario = serde_json::from_value(
                args.get("scenario")
                    .cloned()
                    .ok_or_else(|| anyhow::anyhow!("scenario がありません"))?,
            )?;
            let out_dir = args
                .get("outDir")
                .and_then(|v| v.as_str())
                .map(PathBuf::from)
                .unwrap_or_else(|| PathBuf::from("out").join(&sc.name));
            let manifest = runner::run(&sc, &out_dir)?;
            Ok(text_result(serde_json::to_string_pretty(&json!({
                "manifest": manifest,
                "outDir": out_dir.display().to_string(),
            }))?))
        }
        "validate_scenario" => {
            let sc: Scenario = serde_json::from_value(
                args.get("scenario")
                    .cloned()
                    .ok_or_else(|| anyhow::anyhow!("scenario がありません"))?,
            )?;
            let warnings = lbm_scenario::validate(&sc);
            let build = lbm_scenario::build(&sc);
            Ok(text_result(serde_json::to_string_pretty(&json!({
                "ok": build.is_ok(),
                "error": build.err().map(|e| e.to_string()),
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
