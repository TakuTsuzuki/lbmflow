# LBMFlow 実装プラン

> **2026-07-05 改定**: ユーザー指示により 3D・スパコンスケール・GPU を必須要件へ昇格。
> 以降の計画は [COMPETITIVE_SPEC.md](COMPETITIVE_SPEC.md)（勝ち筋4本柱 + R1-R5）と
> [ARCHITECTURE_V2.md](ARCHITECTURE_V2.md)（次元×格子×精度×バックエンド×分割の直交設計）
> に従う。マイルストーンは M-A（CPU SIMD/wgpu 評価, 進行中）→ M-B（コアV2）→
> M-C（3D）→ M-D（MPI分散）→ M-E（FP16/マルチGPU/公開ベンチ）→ M-F（垂直機能）。

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

## 現行キュー: 改善フェーズ（R-Phase）と M-F（2026-07-05 制定）

原典 2 本: [SOLVER_IMPROVEMENT_SPEC.md](SOLVER_IMPROVEMENT_SPEC.md)（並走レビューの
41 項目・実験 E1〜E10 で全主張検証済み・main 取込済み）と
[REQ_STIRRED_REACTOR.md](REQ_STIRRED_REACTOR.md)（M-F 要求 rev.1b、codex 敵対レビュー
48 件反映済み、受入は VALIDATION.md **T17**）。

### R-Phase（改善仕様書 §4 の実施順序に従う）

- **R-Phase 1（実行中 2026-07-05〜）**: A-2〜A-10 = 入口ガード・正しさ。
  worktree `r-phase1`、Opus 委任。共通 DoD = 既存テスト無修正 green・合法設定ビット不変。
  D-6/D-7（文書整合）は PM 直轄で適用済み（COMPETITIVE_SPEC 改定履歴・VALIDATION T13/T14 節）。
- **R-Phase 2（R-1 着地後・~1.5 週）**: B-1（Fields 一般化 + GpuSolver 統合）→ B-2
  （同期点契約）→ 並行 B-3 / B-5〜B-8、C-9〜C-11、D-1〜D-5。
  **M-F からの追加要求**: B-1 の設計は複数分布セット（相場 g・スカラー h）・per-cell
  物性場（B-6 の一般化）・Lagrangian バッファ（IBM マーカー/粒子/点気泡）を収容できる
  こと — M-F の構造前提であり、M-E（FP16/マルチ GPU）と共通の土台。
- **R-Phase 3（M-E と並走可）**: C-1（MPI セットアップのローカル化 = 10⁹ 格子 OOM 解消）、
  C-2（通信オーバーラップ。前提: E8 の probe 二重計上シェル修正）、C-4〜C-8、
  C-12〜C-16、D-8〜D-10。
- **M-E は B-1/B-2/C-9/C-12/C-13/D-9 完了が前提**（仕様書 §4 の依存関係）。

### M-F: 回転境界・高密度比二相・LES 連成 3D（REQ-M-F-STR rev.1b）

確定済み設計判断（プロジェクトオーナー決定）: **スコープ一括**（サブシステムを段階分割
せず同時実装）/ **忠実度既定**（IBM-inertial・resolved-phasefield・active スカラー・
two-way 粒子・uniform 格子・界面近傍 f64+バルク f32）、低コスト近似（MRF・point-bubble・
one-way・AMR・積極 f32）は同一 trait 背後の後付け拡張 / 物理競合モードは構成検証で
実行時相互排他（A-4 の `GlobalSpec::validate` を拡張）。

実装トラック（R-Phase 2 着地後に並列発注。worktree 分離・実装 Opus/Sonnet・
**検証テストは codex が REQ/T17 から敵対的に作成**の従来体制）:

| トラック | 内容 | 主要 FR | 検証 | 依存 |
|---|---|---|---|---|
| MF-α | D3Q27 格子 + 中心モーメント/cumulant 衝突 | FR-CORE-01/02 | モーメント等方性・TGV3D 次数・ガリレイ不変帯 | R-2 |
| MF-β | LES（WALE 既定）+ 非 Newton μ(γ̇) + 応力場評価（規約 FR-STRESS-01 固定済み） | FR-LES-*, FR-STRESS-* | VR-STR-03、チャネル Re_τ=180 vs DNS | R-2（B-6） |
| MF-γ | 保存型 Allen-Cahn 高密度比二相（10³）+ well-balanced 重力 + スパージャ/脱気 BC | FR-VOF-*, FR-BC-* | VR-STR-02/06、Laplace・寄生流 Ca<10⁻³・単一気泡 Grace | R-2（B-1 複数分布） |
| MF-δ | IBM-inertial 回転境界 + トルク/Np 測定 | FR-ROT-* | VR-STR-01、Taylor-Couette トルク・IBM 球抗力 vs T15 基準 | R-2（B-1 Lagrangian） |
| MF-ε | スカラー ADE（active 帰還）+ Lagrangian 粒子（せん断曝露記録） | FR-LES-04, FR-PART-* | VR-STR-04、Taylor-Aris・沈降終端速度 vs SN | R-2（B-1） |
| MF-ζ | 連成統合（§5 データフロー・dt 制約）+ 構成排他 + I/O/統計/GUI 3D 表示 + 受入ラン | FR-COUP-*, FR-INIT, FR-IO-* | VR-STR-05/07 + 01/02 の連成系 | MF-α〜ε |

残仕様詰めの担当割当: **active スカラー帰還式** = リサーチ委任中（成果 →
docs/proposals/active-scalar-feedback.md → PM レビュー → REQ rev.2 へ反映）/
**f64/f32 界面帯幅** = MF-γ 実装時に実験で凍結（characterize→freeze）/
**trait 境界 API** = R-Phase 2 の B-1 設計と同時に確定（ARCHITECTURE_V2 に追記）/
**REQ 第 2 次 codex 検証** = rev.1b に対して発注（原則追加後の残存齟齬確認）。

工数感（並列エージェント前提）: R-2 ~1.5 週 → MF-α〜ε 並列 ~1-2 週 → MF-ζ 統合 ~1 週。
1e9 格子級はメモリ予算表（REQ §7）により**クラスタ専用**（単機開発線は ≤256³）—
実測はクラスタ計画（CLUSTER_OPTIONS.md、ユーザー判断待ち）に統合。

## 進捗メモ

- 2026-07-05 深夜: **並走レビューセッションの成果を main へ統合**。改善仕様書 v1 +
  実験クレート取込（E2/E7 が改名後 main で数値一致再現）。PM 判断 4 件確定:
  (a) 仕様書取込 = PM 実施済み (b) R-Phase 2 は R-1 直後・M-E/M-F 共通前提
  (c) D-6 = PM 直轄で適用済み（COMPETITIVE_SPEC 改定履歴参照）(d) codex D-8 発注は
  R-1 着地後。**R-Phase 1 発注**（worktree r-phase1・Opus・A-2〜A-10）。
  **M-F 要求 rev.1b 確定**（表題中立化・メモリ予算表・T17 配線）と実装トラック計画
  （上表）。codex #7（T15.5 3D キャビティ）実行中。
- 2026-07-05 晩: **R1/R2/R3 全達成を公式判定**（REVIEW_2026-07-05_2.md。
  ※受入帯は 2026-07-05 D-6 改定を正とする: 球抗力 ±10%・D_h 正規化、弱スケ
  85% 線は n≤4 局所、3D キャビティ T15.5 は codex #7 で追実装中）。
  M-D MPI 完了（T13-MPI 全ケースビット一致、弱スケ 97-99% n≤4、MPI_GUIDE 完備）。
  CpuSimd 融合バックエンド（2D 新記録 1,183 MLUPS、3D 2-3x、等価 ~6e-14）。
  ワークスペース 205 テスト緑。進行中: V1 引退（v1-retirement）+ リサーチ 3 本
  （公開ベンチ比較 / 3D キャビティ参照 / クラスタ選択肢）。CI workflow 準備済み。
- 2026-07-05 夕: **M-C 3D物理 完了（R1達成）**。D3Q19 Zou-He 面境界（未知5本+接線補正、
  2D退化 8.9e-16）、ダクト厳密級数 2.3e-4、球抗力は流体力学的ペア (r+½, Re(D+1)/D)
  正規化で +7.1%/+0.6%/+2.3%（重量級D=24含む全合格）、TGV3D 次数1.910。
  シナリオ/CLI nz 対応（VTK 3D/断面PNG）。**Wgpuバックエンド統合（R2）**: push型融合で
  CPU と演算子順序厳密一致、T14 6構成 ≤1e-5、5.9〜11.4 GLUPS。ワークスペース184緑。
  進行中: M-D MPI分散（arm64 Open MPI 5.0.9 を ~/.local にソースビルド済み）。
- 2026-07-05 午後: **M-B コアV2 完了・統合**（V1と物理等価、T13分割不変は敵対攻撃
  8種＋3D 2×2×2までビット一致級、D3Q19スモーク動作、V1暗黙仕様8件をテスト凍結）。
  **Phase 9 完了**: CPU融合カーネル 3.2〜7倍（f32峰1,124 MLUPS）+ GPU実測
  6,975〜12,152 MLUPS。MCP非同期ジョブAPI（R4）。codex #6 の敵対T13は全て耐えた。
  進行中: Wgpuバックエンド本実装（m-b-wgpu）と3D物理 M-C（m-c-3d）の2並列。
- 2026-07-05 昼: **Phase 8 完了**（T12 RT γ比1.118・T11c 接触角フルレンジ・T9b 対流流出は
  反射ケースで16倍改善、67テスト緑）。**GPU 実測完了**: M5 Max Metal で 6,975〜12,152
  MLUPS（CPU比16〜42倍、検証L∞ 7e-6、GPU_EVALUATION.md、R2目標を4-8倍超過）。
  シナリオ/CLI に convectiveOutflow・wallRho・VTK・gallery。GUI にシナリオ書き出し
  （→lbm run E2E確認）・発散ガード・MLUPS。SoA/SIMD は WIP（phase9-perf、クォータ
  回復後に再開）。次: M-B コアV2 の並列発注。

- 2026-07-05: **Phase 4a/5/6 完了（三モード統一達成）**。codex #4 の多相検証
  全緑（共存密度・EOS 圧力平衡・Laplace・接触角回帰 133/160/164°・f32 強化 1e-5）。
  GUI: WASM 実エンジンでキャビティ/カルマン渦/二相液滴が動作（~600 steps/s）。
  Agent モード: lbm CLI（run/validate/presets/schema, manifest+PNG/CSV）+
  MCP サーバー（run_scenario 等 4 ツール）。ワークスペース 56 テスト全緑。
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
