# Agent モード設計（Phase 6）

エージェント（Claude/Codex 等）と CLI ユーザーが同じ入口でシミュレーションを
実行できるようにする。**シナリオ JSON が唯一の実行契約**であり、GUI のプリセットも
内部的に同じスキーマを生成する（3 モード統一）。

## クレート構成

- `crates/lbm-cli` → バイナリ名 `lbm`
  - `lbm run <scenario.json> [--out DIR]` : 実行し、成果物 + manifest.json を出力
  - `lbm validate <scenario.json>` : スキーマ・物理妥当性チェックのみ（実行しない）
  - `lbm presets [list|show NAME|run NAME]` : 組み込みプリセット
  - `lbm schema` : シナリオ JSON Schema を stdout に出力（エージェントの自己発見用）
  - `lbm mcp` : stdio で MCP サーバーを起動
- 依存: serde/serde_json, clap, png（可視化出力）, （MCP は rmcp か手書き JSON-RPC）

## シナリオ JSON スキーマ（v0）

```jsonc
{
  "version": 0,
  "name": "cylinder-re100",                 // 出力ディレクトリ名にも使用
  "grid": { "nx": 440, "ny": 160 },
  "physics": {
    "nu": 0.01,                              // 格子単位
    "collision": { "type": "trt", "magic": 0.1875 },   // {"type":"bgk"} も可
    "force": [0.0, 0.0],
    "precision": "f64"                       // "f32" | "f64"（精度/速度トレードオフ）
  },
  "edges": {                                 // lbm-core の EdgeBC と 1:1
    "left":   { "type": "velocityInlet", "u": [0.05, 0.0] },
    "right":  { "type": "outflow" },
    "bottom": { "type": "bounceBack" },
    "top":    { "type": "bounceBack" }
  },
  "inletProfile": { "edge": "left", "kind": "parabolic", "umax": 0.05 },  // 省略可
  "obstacles": [                             // set_solid_region に展開
    { "shape": "circle", "cx": 110, "cy": 80, "r": 12 },
    { "shape": "rect", "x0": 0, "y0": 0, "x1": 10, "y1": 5 }
  ],
  "init": { "kind": "rest" },                // rest | taylorGreen | custom(将来)
  "run": {
    "steps": 100000,
    "stopWhenSteady": { "epsilon": 1e-11, "checkEvery": 500 }   // 省略可
  },
  "probes": [                                // 時系列の記録
    { "type": "force", "target": "obstacles", "every": 10 },
    { "type": "point", "x": 220, "y": 80, "fields": ["ux","uy","rho"], "every": 100 }
  ],
  "outputs": [
    { "field": "speed", "at": "end", "format": "png", "colormap": "viridis" },
    { "field": "ux",    "at": "end", "format": "csv" },
    { "snapshotEvery": 10000, "field": "vorticity", "format": "png" }
  ]
}
```

### 実行結果（out ディレクトリ）

- `manifest.json`: 機械可読サマリ
  ```jsonc
  {
    "scenario": "cylinder-re100",
    "status": "completed" | "diverged" | "error",
    "steps": 100000, "wallSeconds": 42.1, "mlups": 168.3,
    "diagnostics": { "totalMass": ..., "maxSpeed": ..., "reynolds": ... },
    "warnings": ["..."],
    "files": [ {"path": "speed_100000.png", "kind": "field", ...} ]
  }
  ```
- probes は CSV（`force.csv`: step,fx,fy）
- 途中で NaN 検出 → status="diverged" + 直前の診断を残して即終了（エージェントが
  パラメータを直して再試行できるよう、原因ヒント文字列を含める）

### 設計原則（エージェント UX）

1. **自己記述**: `lbm schema` と `lbm presets list` だけで使い方を発見できる
2. **失敗が構造化**: バリデーションエラーは JSON で「どのフィールドが・なぜ・
   どう直すか」を返す（例: "nu must be > 0; tau = 3*nu + 0.5 must exceed 0.5"）
3. **数値の安全柵**: validate 段階で安定性ヒューリスティック（|u|>0.15 警告、
   グリッドレイノルズ数 U/ν > 15 警告など）を warnings に出す
4. **決定論**: 同じシナリオ → 同じ結果（シード不要の設計を維持）

## MCP サーバー（`lbm mcp`）

ツール:
- `run_scenario(scenario: object, outDir?: string) -> manifest`（進捗 notification 対応）
- `validate_scenario(scenario: object) -> {ok, errors[], warnings[]}`
- `list_presets() -> [{name, description, scenario}]`
- `get_schema() -> JSON Schema`
- `read_field(runDir, field, format="csv") -> data`（実行済み結果の取得）

## GUI との関係

- GUI のプリセット = シナリオ JSON（web/src/presets.ts から生成 or 共有 JSON）
- 将来: GUI に「シナリオを書き出す」ボタン → Agent モードで再現可能
