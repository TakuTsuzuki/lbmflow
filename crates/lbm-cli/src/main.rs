//! `lbm` — LBMFlow のシナリオ実行 CLI（Agent モードの入口）。
//!
//! エージェント向け設計原則: 自己記述（`lbm schema` / `lbm presets`）、
//! 構造化エラー（JSON）、決定論。

mod gallery;
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
    /// 全プリセットを順に実行し、自己完結 HTML ギャラリー（index.html）を生成する
    Gallery {
        /// 出力ディレクトリ（既定: out/gallery）
        #[arg(long)]
        out: Option<PathBuf>,
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
            let build_result = lbm_scenario::build_check(&sc);
            let report = serde_json::json!({
                "ok": build_result.is_ok(),
                "error": build_result.err(),
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
        Command::Gallery { out } => {
            let out_root = out.unwrap_or_else(|| PathBuf::from("out").join("gallery"));
            gallery::run(&out_root)?;
        }
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
  "grid": { "nx": 128, "ny": 128 },        // "nz": 64 を足すと 3D (D3Q19) で実行
                                           //   （省略 or 1 = 2D。3D の制約は末尾参照）
  "physics": {
    "nu": 0.02,                            // 動粘性係数（格子単位）。tau = 3*nu + 0.5
    "collision": { "type": "trt" },        // "trt"（推奨） | "bgk"
    "force": [0.0, 0.0],                   // 一様体積力（重力など。3D では z 成分 0）
    "precision": "f64"                     // "f32" | "f64"
  },
  "compute": { "backend": "auto" },        // 省略可: "auto" | "cpu" | "gpu"（gpu は未提供）
  "edges": {                               // 4辺の境界条件
    "left":   { "type": "velocityInlet", "u": [0.1, 0.0] },
    "right":  { "type": "pressureOutlet", "rho": 1.0 },
    "bottom": { "type": "bounceBack" },
    "top":    { "type": "bounceBack" }
    // 3D では "front"（z=0）/"back"（z=nz-1）を追加可（省略 = periodic）。
    // 他: {"type":"periodic"}（対辺ペア必須）, {"type":"movingWall","u":[ux,uy]},
    //     {"type":"outflow"},
    //     {"type":"convectiveOutflow","uConv":0.1}  // 対流流出。outflow より圧力反射が
    //       小さい。uConv = 期待平均流出速度（0 < uConv <= 1、流入速度と同程度が目安）
    // 制約: 開境界(velocityInlet/pressureOutlet/outflow/convectiveOutflow)同士は
    //       直交して隣接不可（3D では開境界は 1 軸のみ）
  },
  "inletProfile": {                        // 省略可: 放物線流入
    "edge": "left", "kind": "parabolic", "umax": 0.15
    // 3D: 壁に挟まれた接線軸ごとに放物線、periodic 軸は一様
    //     （4壁ならダクト型 u = umax·f(y)·f(z)）
  },
  "obstacles": [                           // 省略可
    { "shape": "circle", "cx": 80, "cy": 80, "r": 20 },   // 3D では z 方向に押し出し（円柱）
    { "shape": "rect", "x0": 10, "y0": 10, "x1": 20, "y1": 40 },
    { "shape": "sphere", "cx": 60, "cy": 32, "cz": 32, "r": 12 }  // 3D 専用
  ],
  "init": { "kind": "rest" },              // rest | droplet{cx,cy,r,rhoLiquid,rhoVapor}
                                           // | pool{heightFrac,rhoLiquid,rhoVapor}
  "multiphase": {                          // 省略可: Shan-Chen 単成分多相
    "g": -5.0,                             // 凝集強度（負。-5.0 が検証済み既定）
    "gWall": 0.0,                          // 壁付着（負=濡れ性）。wallRho の方を推奨
    "wallRho": 1.0                         // 省略可: 仮想壁密度による接触角制御
                                           // （液密度側 → 濡れる。0.3:~180°, 0.6:~107°, 1.0:~63°）
  },
  "run": {
    "steps": 20000,
    "stopWhenSteady": { "epsilon": 1e-8, "checkEvery": 500 }  // 省略可
  },
  "probes": [                              // 省略可: 時系列 CSV
    { "type": "force", "every": 10 },      // 障害物への力（force.csv。3D は fx,fy,fz）
    { "type": "point", "x": 220, "y": 80, "every": 100 }  // 3D は "z" も指定可（省略 = nz/2）
  ],
  "outputs": [                             // 省略可: 場のスナップショット
    { "field": "speed", "format": "png", "every": 0 }   // every=0 は終了時のみ
    // field: speed | ux | uy | rho | vorticity
    // format: png | csv | vtk（VTK legacy structured points。ParaView 等で開ける）
    // 3D: png/csv は z 中央断面、vtk は 3D 全体（DIMENSIONS nx ny nz）
  ]
}

3D (nz > 1) の制約: 単相のみ（multiphase 不可）、init は rest のみ、
compute.backend は cpu/auto（gpu は未提供）。エンジンは V2 コア（D3Q19）。

結果: <out>/manifest.json（status/steps/mlups/診断/警告/ファイル一覧）
status: completed | steady（定常判定で早期終了）| diverged（NaN 検出）
実例: lbm presets show cavity | cylinder-karman | two-phase-droplet | droplet-on-wall
一括実行: lbm gallery --out DIR（全プリセット + 自己完結 HTML ギャラリー）

MCP: lbm mcp で MCP サーバー（stdio）。run_scenario は同期実行（完了までブロック）。
長時間ランやスイープは非同期 API を使う:
  start_run { scenario, outDir? } -> { runId }   … 即応答、バックグラウンド実行
  run_status { runId } -> { state: running|completed|failed, manifest?, error? }
  list_runs {} -> 全ランの一覧
runId は "run-<連番>-<シナリオ名>"（決定論）。同時実行は最大 4 ラン
（超過は "failed: too many concurrent runs" で即時拒否）。ランはサーバー
プロセス内で動くため、完了確認まで MCP 接続を維持すること。
"#;
