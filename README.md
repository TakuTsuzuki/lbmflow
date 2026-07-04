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
- **3 つの使い方**（開発中）: ブラウザ GUI（WASM）/ CLI（JSON シナリオ）/
  MCP サーバー（AI エージェントから操作）
- 混相流（Shan-Chen）対応は Phase 4 で追加予定

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
