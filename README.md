# LBMFlow — 格子ボルツマン法 流体シミュレータ

Rust 製の高速な LBM（Lattice Boltzmann Method, D2Q9）流体シミュレーションエンジンと、
ブラウザ GUI・エージェント連携（CLI / MCP）を備えた統合環境。

## 特徴

- **検証済みの物理**: Taylor-Green 渦の2次収束、TRT(Λ=3/16) による Poiseuille 流の
  厳密再現、Ghia+1982 キャビティベンチマーク等、敵対的に作成された検証スイートを
  全て通過（仕様: [docs/VALIDATION.md](docs/VALIDATION.md)）
- **精度と速度のトレードオフを明示的に制御**:
  - 衝突演算子: BGK（速い）⇔ TRT（高精度・高安定、推奨）
  - 数値精度: `f32`（速い・省メモリ）⇔ `f64`（検証グレード）
  - 並列度: rayon マルチスレッド（小格子は自動でシリアル実行）
- **豊富な境界条件**: 周期、half-way bounce-back（静止壁・移動壁）、
  Zou-He 速度流入（一様・任意プロファイル）、Zou-He 圧力、ゼロ勾配流出、
  任意形状の内部障害物、momentum-exchange 抗力測定
- **3 つの使い方**: ブラウザ GUI（WASM）/ CLI（JSON シナリオ）/
  MCP サーバー（AI エージェントから操作）
- **混相流対応**: Shan-Chen 単成分多相（液滴・接触角、検証済み）

## 使い方 1: ブラウザ GUI（初学者向け）

```bash
cd web && npm install && npm run dev   # → http://localhost:5173
```

プリセット（キャビティ流れ / カルマン渦列 / ポアズイユ / 二相液滴 / 自由キャンバス）を
選んで「▶実行」するだけ。障害物はマウスで描けます。本物の LBM（Rust→WASM）が
ブラウザ内で毎秒約 60 万格子点更新で動きます。

## 使い方 2: CLI（シナリオ実行）

```bash
cargo build --release -p lbm-cli
./target/release/lbm presets list                 # 組み込みプリセット
./target/release/lbm presets run cylinder-karman  # 実行 → out/ に PNG/CSV/manifest.json
./target/release/lbm schema                       # シナリオ JSON の書式
./target/release/lbm run my-scenario.json --json  # 自作シナリオ
```

## 使い方 3: MCP サーバー（AI エージェント連携）

```bash
claude mcp add lbmflow -- /path/to/target/release/lbm mcp
```

エージェントは `run_scenario` / `validate_scenario` / `list_presets` /
`get_schema` の 4 ツールでシミュレーションを実行し、構造化された結果
（manifest + PNG/CSV）を受け取れます。

## クイックスタート（ライブラリ）

```rust
use lbm_core::prelude::*;

// リッド駆動キャビティ
let mut sim: Simulation<f64> = SimConfig {
    nx: 128, ny: 128,
    nu: 0.02,
    edges: Edges {
        left: EdgeBC::BounceBack,
        right: EdgeBC::BounceBack,
        bottom: EdgeBC::BounceBack,
        top: EdgeBC::MovingWall { u: [0.1, 0.0] },
    },
    ..Default::default()
}.build()?;

sim.run(10_000);
println!("centre velocity = {}", sim.ux(64, 64));
```

## 開発

```bash
cargo test --release                       # 検証スイート（必ず --release で）
cargo test --release -- --include-ignored  # 重いベンチマーク込みフル検証
```

- 計画・体制: [docs/PLAN.md](docs/PLAN.md)
- 検証仕様（受入基準）: [docs/VALIDATION.md](docs/VALIDATION.md)
- 物理モデルと実験記録: [docs/PHYSICS.md](docs/PHYSICS.md)
- Agent モード設計: [docs/AGENT_MODE_DESIGN.md](docs/AGENT_MODE_DESIGN.md)
- 混相流設計: [docs/MULTIPHASE_DESIGN.md](docs/MULTIPHASE_DESIGN.md)

## ライセンス

MIT OR Apache-2.0
