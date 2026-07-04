# LBMFlow 実装プラン

**目標**: 商用グレードの格子ボルツマン法（LBM）流体シミュレータ。
精度と速度のトレードオフを明示的に制御でき、混相流に対応し、初学者でも迷わず使える
GUI モードと、エージェントから操作できる Agent モードの両方を備える。

**技術スタック**: Rust（コアエンジン / CLI / MCP）+ TypeScript（GUI, Vite）+ WASM（ブラウザ実行）

**体制**: Fable = PM / アーキテクト / 統合検証。実装は Claude サブエージェント・Codex に
実行可能な粒度で委任。**境界条件テストスイートは Codex が敵対的に作成**し、
エンジン側はそれを全て通過するまで修正する（テスト作者と実装者を分離して仕様バグを検出）。

---

## アーキテクチャ

```
流体シミュレータ/
├── Cargo.toml              # workspace
├── crates/
│   ├── lbm-core/           # コアエンジン（純Rustライブラリ、no I/O）
│   │   ├── src/
│   │   │   ├── lattice.rs      # D2Q9 定数（速度・重み・反対方向）
│   │   │   ├── real.rs         # f32/f64 ジェネリック（Real トレイト）
│   │   │   ├── domain.rs       # 領域・エッジ境界条件・障害物マスク
│   │   │   ├── collision.rs    # BGK / TRT（/ 将来 MRT・cumulant）
│   │   │   ├── sim.rs          # Simulation: collide→stream→BC の1ステップ
│   │   │   ├── multiphase.rs   # Shan-Chen（Phase 4）
│   │   │   └── analysis.rs     # 誤差ノルム・保存量・力測定ヘルパ
│   │   └── tests/          # 検証テスト（codex 作成の敵対的スイート含む）
│   ├── lbm-cli/            # JSONシナリオ実行 CLI（Agent モードの土台）
│   └── lbm-wasm/           # wasm-bindgen バインディング
├── web/                    # TypeScript GUI (Vite)
├── mcp/                    # MCP サーバー（Agent モード）
└── docs/
    ├── PLAN.md             # 本ファイル
    ├── VALIDATION.md       # 検証テスト仕様マトリクス（codex への発注仕様）
    └── PHYSICS.md          # 採用した物理モデル・数式の根拠（随時追記）
```

### コア設計の要点

- **格子**: D2Q9（2D）から開始。データ配置は cell-major AoS `f[cell*9 + q]`
  （rayon 行並列と安全に両立、WASM でも同一コード）。
- **ストリーミング**: pull 方式（gather）。collide（in-place）→ stream（f→f_tmp）→ swap。
- **衝突演算子**: BGK（速い/低安定）と TRT（magic Λ=3/16 で壁位置厳密・推奨デフォルト）。
  精度と速度のトレードオフ軸①。将来 MRT / cumulant を追加。
- **精度**: `Simulation<T: Real>` で f32/f64 を切替（トレードオフ軸②）。
- **並列**: rayon（feature "parallel"、WASM では off）。トレードオフ軸③（スレッド数）。
- **体積力**: Guo forcing（2次精度、u は F/2 補正込み）。Shan-Chen でも使う。
- **壁**: half-way bounce-back（静止壁・移動壁）。エッジ指定の壁は 1 セルのソリッド
  リムとして実現 → Zou-He とのコーナー特殊処理が不要になる。
- **開境界**: Zou-He（速度流入・圧力流出）を面法線パラメタ化した単一実装で 4 辺対応。
  Outflow（ゼロ勾配コピー）も提供。
- **力測定**: momentum-exchange 法（円柱 Cd/St ベンチに必要）。

### 境界条件の組合せ規則（仕様）

- Periodic は対向エッジでペア必須。
- Zou-He / Outflow エッジの直交エッジは Wall（リム）か Periodic であること
  （エッジ同士が裸で接するコーナーは非サポート、構築時にエラー）。
- τ ≤ 0.5 は構築時エラー。|u| > 0.3（格子単位）は警告。

---

## フェーズ計画

| Phase | 内容 | 完了条件 |
|---|---|---|
| 0 | 基盤: git / workspace / PLAN / VALIDATION / CLAUDE.md | 文書コミット済み |
| 1 | lbm-core 縦切り（D2Q9, BGK/TRT, Guo力, BB/移動壁/周期/Zou-He/Outflow, 力測定）+ スモークテスト | TGV 収束次数≈2、Poiseuille(TRT)厳密、Couette 厳密、保存則テスト green |
| 2 | **codex に VALIDATION.md の敵対的テストスイートを実装させる** → 全通過までエンジン修正。キャビティ(Ghia)・円柱(Cd/St) ベンチ含む | `cargo test --release` 全 green |
| 3 | 精度×速度: f32 検証、MLUPS ベンチ、スレッドスケーリング、モード選択ガイド、**偏差格納方式**（f−w を保持して f32 の有効精度を引き上げる。BB は線形なので偏差空間で不変、Zou-He は定数項の折込みが必要、ρ = 1+Σdev） | トレードオフ実測表が docs に載る。f32 の T6 誤差が改善 |
| 4 | 混相流: Shan-Chen（単成分多相 + 二成分）+ codex 検証（Laplace 則、接触角、Rayleigh-Taylor） | 混相テスト全 green |
| 5 | GUI: wasm-pack + Vite + TS。プリセット駆動、障害物ペイント、リアルタイム可視化、日本語 UI | プリセット 5 種がワンクリックで動く |
| 6 | Agent モード: lbm-cli（JSON シナリオ→構造化出力）+ MCP サーバー | エージェントがシナリオ実行→結果取得できる |
| 7 | 総合レビュー → 次期プラン（3D/GPU/LES/cumulant 等）策定 → 継続改善 | レビュー文書 + 次期プラン |

### 検証駆動の開発プロトコル

1. 仕様（VALIDATION.md）に受入基準を数値で明記。
2. codex がテストを書く（実装を見ずに仕様から書くことを原則とする）。
3. エンジンが落ちたら: (a) エンジンのバグ → 修正、(b) 仕様の物理が誤り → 実験して仕様を修正
   （修正理由を PHYSICS.md に記録）、(c) テストのバグ → codex に差し戻し。
4. 全 green で次フェーズへ。フェーズ末に `git commit`。

## 進捗メモ

- 2026-07-04: プロジェクト開始。Phase 0 着手。
- 2026-07-04: Phase 1 完了。lbm-core（D2Q9, BGK/TRT, Guo力, half-way BB/移動壁,
  Zou-He, Outflow, 力測定, f32/f64, rayon+小格子シリアルフォールバック）実装。
  スモークテスト 21 件 green: TGV 2次収束（1.91/1.98）、Poiseuille TRT 厳密
  (<1e-10)、Couette 厳密、保存則 ~1e-13。実験知見 4 件を PHYSICS.md に記録。
- 2026-07-05: **Phase 4a 実装完了**（検証は codex #4 待ち）。セル別力場 API +
  Shan-Chen SCMP。実測: 密度比 15.8、圧力平衡 8.5e-6、疑似速度 1.3e-3、
  Laplace R²=0.9999。**スコープ再編**: MCMP+RT（Phase 4b）は初回レビュー後へ。
  GUI（Phase 5）/ Agent モード（Phase 6）を先行し 3 モード統一を先に完成させる。
- 2026-07-05: **Phase 3 完了**。偏差格納方式（f−w）導入で f32 が検証グレードに
  （運動量誤差 4800 倍改善・TGV で f64 同等）。MLUPS 実測: 峰 381（f32/1024²/18T）、
  シングル 35、TRT は BGK と同速。PERFORMANCE.md にトレードオフガイド。
  全 49 テスト green 維持。
- 2026-07-05: **Phase 2 完了**。codex 敵対テスト 3 巡（計 9 指摘を triage:
  エンジンバグ 2 = Zou-He 圧力符号・リムコーナー異方性 / 仕様バグ 5 / 参照データ
  誤植 1 / f32 特性 1）。デフォルト 49 テスト・フル 53 テスト全 green。
  Ghia キャビティ Re=100/400/1000、Schäfer-Turek 円柱 2D-1/2D-2、
  厳密等変性（機械精度 4e-16）、Zou-He 4 方向、保存則を通過。
  GUI シェル（Vite+TS+モック）も先行完成（Phase 5a）。
