//! `lbm` — LBMFlow のシナリオ実行 CLI（Agent モードの入口）。
//!
//! エージェント向け設計原則: 自己記述（`lbm schema` / `lbm presets`）、
//! 構造化エラー（JSON）、決定論。

mod mcp;
mod render;
mod runner;

use anyhow::{Context, Result};
use clap::{Parser, Subcommand};
use lbm_scenario::Scenario;
use std::path::PathBuf;

#[derive(Parser)]
#[command(
    name = "lbm",
    about = "LBMFlow 格子ボルツマン法シミュレータ CLI",
    version
)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// シナリオ JSON を実行し、結果を出力ディレクトリに書き出す
    Run {
        /// シナリオ JSON ファイル（`-` で stdin）
        scenario: String,
        /// 出力ディレクトリ（既定: out/<scenario name>）
        #[arg(long)]
        out: Option<PathBuf>,
        /// 結果 manifest を stdout に JSON で出す
        #[arg(long)]
        json: bool,
    },
    /// シナリオを実行せずに検証する（エラー/警告を JSON で報告）
    Validate {
        /// シナリオ JSON ファイル（`-` で stdin）
        scenario: String,
    },
    /// 組み込みプリセットの一覧・表示・実行
    Presets {
        #[command(subcommand)]
        action: PresetAction,
    },
    /// シナリオ JSON の書式説明を出力する（エージェントの自己発見用）
    Schema,
    /// MCP サーバーとして stdio で待ち受ける（AI エージェント連携）
    Mcp,
}

#[derive(Subcommand)]
enum PresetAction {
    /// 一覧を表示
    List,
    /// プリセットのシナリオ JSON を表示
    Show { name: String },
    /// プリセットを実行
    Run {
        name: String,
        #[arg(long)]
        out: Option<PathBuf>,
    },
}

fn load_scenario(path: &str) -> Result<Scenario> {
    let text = if path == "-" {
        std::io::read_to_string(std::io::stdin())?
    } else {
        std::fs::read_to_string(path).with_context(|| format!("読めません: {path}"))?
    };
    let sc: Scenario = serde_json::from_str(&text).map_err(|e| {
        anyhow::anyhow!(serde_json::to_string_pretty(&serde_json::json!({
            "error": "invalid-scenario-json",
            "message": e.to_string(),
            "hint": "lbm schema で書式を、lbm presets show <name> で実例を確認できます"
        }))
        .unwrap())
    })?;
    Ok(sc)
}

fn run_and_report(sc: &Scenario, out: Option<PathBuf>, json: bool) -> Result<()> {
    let out_dir = out.unwrap_or_else(|| PathBuf::from("out").join(&sc.name));
    let manifest = runner::run(sc, &out_dir)?;
    if json {
        println!("{}", serde_json::to_string_pretty(&manifest)?);
    } else {
        println!(
            "status={} steps={} wall={:.1}s mlups={:.0} out={}",
            manifest.status,
            manifest.steps_run,
            manifest.wall_seconds,
            manifest.mlups,
            out_dir.display()
        );
        for w in &manifest.warnings {
            eprintln!("警告[{}]: {}", w.field, w.message);
        }
        for f in &manifest.files {
            println!("  {f}");
        }
    }
    Ok(())
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    match cli.command {
        Command::Run {
            scenario,
            out,
            json,
        } => {
            let sc = load_scenario(&scenario)?;
            run_and_report(&sc, out, json)?;
        }
        Command::Validate { scenario } => {
            let sc = load_scenario(&scenario)?;
            let warnings = lbm_scenario::validate(&sc);
            let build_result = lbm_scenario::build(&sc);
            let report = serde_json::json!({
                "ok": build_result.is_ok(),
                "error": build_result.err().map(|e| e.to_string()),
                "warnings": warnings,
            });
            println!("{}", serde_json::to_string_pretty(&report)?);
        }
        Command::Presets { action } => match action {
            PresetAction::List => {
                for (name, desc, _) in lbm_scenario::presets() {
                    println!("{name:<20} {desc}");
                }
            }
            PresetAction::Show { name } => {
                let all = lbm_scenario::presets();
                let found = all.iter().find(|(n, _, _)| *n == name).ok_or_else(|| {
                    anyhow::anyhow!(
                        "プリセット '{name}' はありません。lbm presets list で一覧を確認してください"
                    )
                })?;
                println!("{}", serde_json::to_string_pretty(&found.2)?);
            }
            PresetAction::Run { name, out } => {
                let all = lbm_scenario::presets();
                let found = all
                    .iter()
                    .find(|(n, _, _)| *n == name)
                    .ok_or_else(|| anyhow::anyhow!("プリセット '{name}' はありません"))?;
                run_and_report(&found.2, out, false)?;
            }
        },
        Command::Schema => {
            println!("{}", SCHEMA_DOC);
        }
        Command::Mcp => {
            mcp::serve()?;
        }
    }
    Ok(())
}

const SCHEMA_DOC: &str = r#"シナリオ JSON (v0) — lbm run <file.json>

{
  "version": 0,
  "name": "my-sim",                       // 出力ディレクトリ名
  "grid": { "nx": 128, "ny": 128 },
  "physics": {
    "nu": 0.02,                            // 動粘性係数（格子単位）。tau = 3*nu + 0.5
    "collision": { "type": "trt" },        // "trt"（推奨） | "bgk"
    "force": [0.0, 0.0],                   // 一様体積力（重力など）
    "precision": "f64"                     // "f32" | "f64"
  },
  "edges": {                               // 4辺の境界条件
    "left":   { "type": "velocityInlet", "u": [0.1, 0.0] },
    "right":  { "type": "pressureOutlet", "rho": 1.0 },
    "bottom": { "type": "bounceBack" },
    "top":    { "type": "bounceBack" }
    // 他: {"type":"periodic"}（対辺ペア必須）, {"type":"movingWall","u":[ux,uy]},
    //     {"type":"outflow"}
    // 制約: 開境界(velocityInlet/pressureOutlet/outflow)同士は直交して隣接不可
  },
  "inletProfile": {                        // 省略可: 放物線流入
    "edge": "left", "kind": "parabolic", "umax": 0.15
  },
  "obstacles": [                           // 省略可
    { "shape": "circle", "cx": 80, "cy": 80, "r": 20 },
    { "shape": "rect", "x0": 10, "y0": 10, "x1": 20, "y1": 40 }
  ],
  "init": { "kind": "rest" },              // rest | droplet{cx,cy,r,rhoLiquid,rhoVapor}
                                           // | pool{heightFrac,rhoLiquid,rhoVapor}
  "multiphase": { "g": -5.0, "gWall": 0.0 }, // 省略可: Shan-Chen 単成分多相
  "run": {
    "steps": 20000,
    "stopWhenSteady": { "epsilon": 1e-8, "checkEvery": 500 }  // 省略可
  },
  "probes": [                              // 省略可: 時系列 CSV
    { "type": "force", "every": 10 },      // 障害物への力（force.csv）
    { "type": "point", "x": 220, "y": 80, "every": 100 }
  ],
  "outputs": [                             // 省略可: 場のスナップショット
    { "field": "speed", "format": "png", "every": 0 }   // every=0 は終了時のみ
    // field: speed | ux | uy | rho | vorticity, format: png | csv
  ]
}

結果: <out>/manifest.json（status/steps/mlups/診断/警告/ファイル一覧）
status: completed | steady（定常判定で早期終了）| diverged（NaN 検出）
実例: lbm presets show cavity | cylinder-karman | two-phase-droplet
"#;
