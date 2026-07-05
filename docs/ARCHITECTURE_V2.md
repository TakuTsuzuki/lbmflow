# アーキテクチャ V2 — 3D・GPU・分散を単一コアで（2026-07-05）

COMPETITIVE_SPEC.md の必須要件（R1-R5）に整合する設計。現行 V1（2D/CPU/単一領域）を
**検証スイートを緑に保ったまま**段階置換する。

## 0. 設計原則

1. **1つの物理カーネル定義から全ターゲットへ**: 次元・格子・精度・バックエンド・分割は
   直交する軸。物理（衝突・力・境界）を一箇所に書き、軸の組合せはコンパイル時に展開する
2. **等価性が唯一の真実**: 新しい軸を足すたびに「既存構成と一致する」テスト（T13/T14/T16）
   を先に定義する。速さの主張はその後
3. **エージェント契約は不変**: シナリオ JSON / MCP が全構成の共通入口（R4）

## 1. レイヤ構成

```
┌─────────────────────────────────────────────────────┐
│  scenario / CLI / MCP / GUI（契約層: JSON, 変更最小）  │
├─────────────────────────────────────────────────────┤
│  Solver Orchestrator（時間発展・診断・プローブ・出力）  │
├───────────────┬───────────────────┬─────────────────┤
│ Decomposition │  Physics Kernels   │  Diagnostics    │
│ Subdomain     │  collide+stream    │  reductions     │
│ HaloExchange  │  BCs / forces      │  (backend-side) │
├───────────────┴───────────────────┴─────────────────┤
│  Backend trait: CpuSimd | Wgpu | (Cuda: feature)     │
├─────────────────────────────────────────────────────┤
│  Lattice trait: D2Q9 | D3Q19 | (D3Q27)  × Real/f16   │
└─────────────────────────────────────────────────────┘
```

## 2. 各抽象の定義

### 2.1 Lattice（コンパイル時定数）

```rust
pub trait Lattice: Copy + 'static {
    const D: usize;            // 2 | 3
    const Q: usize;            // 9 | 19 | 27
    const C: [[i8; 3]; Self::Q];   // 速度（2Dはz=0）
    const W: [f64; Self::Q];
    const OPP: [usize; Self::Q];
    const CS2: f64;            // 1/3
    // TRT ペア、面ごとの unknown 集合などの派生テーブルも const fn で提供
}
pub struct D2Q9; pub struct D3Q19;
```
- 現行 lattice.rs の「方向順序が唯一の正」原則を trait 化して維持
- Zou-He の面法線パラメタ化（V1 で実証済み）は D3 でも同型（面 unknown = c·n>0 の 5 本）

### 2.2 ストレージ（GPU コアレッシング前提で SoA 固定）

- `f[q][cell]`（q-major SoA）。cell = z·(nx·ny) + y·nx + x（2D は z=0 のみ）
- **偏差格納（f−w）を維持**（V1 で f32 検証グレードを実証。FP16 の前提条件でもある）
- 精度は「演算精度 × 格納精度」の2軸: (f32,f32) / (f64,f64) / **(f32,f16 格納)**（R2）。
  f16 は pack/unpack をバックエンド内に閉じ、API 上は f32 で見せる
- moments キャッシュ（rho,u）は診断・可視化・多相用。バックエンド側メモリに常駐し、
  ホストへは明示 readback（暗黙同期を作らない）

### 2.3 Subdomain / HaloExchange（R3）

```rust
pub struct Subdomain { global_box: Box3, local: Box3, halo: usize /*=1*/,
                       neighbors: [Option<RankId>; 6/*faces*/ + edges…] }
pub trait HaloExchange {
    /// 面ごとに「外向き分布のみ」（D2Q9: 3本/面, D3Q19: 5本/面）+ 多相時は ψ を交換
    fn exchange(&mut self, field: &mut BackendField, plan: &HaloPlan);
}
impl: LocalPeriodic（単一領域=現行動作） / InProcess（スレッド間, T13用） / Mpi（rsmpi）
- (R-Phase 1) `HaloExchange::SCOPE: ExchangeScope { Local, Remote }` — building a
  single-part owner (`only=Some(part)`) with a Local exchange is a construction
  error (silent self-wrap prevention, spec A-5).
- (R-Phase 1) `Backend::stream` contract, pinned by tests/stream_contract.rs:
  streaming must NOT write open-face unknown slots (ConvectiveOutflow memory
  depends on it; GPU realizes it via the edge stash). Any in-place streaming
  (M-E candidate) must preserve this contract or replace the mechanism.
```
- ステップ構造を**内部→境界の2パス**に分け、halo 通信と内部計算をオーバーラップ
  （V1 の行分割ループはこの分割と親和的）
- 大域診断（total_mass 等）は backend 内 reduce → ランク間 Allreduce
- rims/障害物/プローブは Subdomain ローカルなデータに分配（V1 で「BC=データ」に
  してあるため機械的に分割可能）

### 2.4 Backend

```rust
pub trait Backend<L: Lattice> {
    type Field;   // デバイス常駐の f / moments / mask
    fn step(&mut self, dom: &Subdomain, fields: &mut Fields<Self>, params: &StepParams);
    fn reduce(&self, kind: Reduction) -> f64;
    fn read_moments(&self, out: &mut HostMoments);   // 明示 readback
}
```
- **CpuSimd**: 現行 rayon 実装を SoA 化+SIMD（phase9-perf ブランチの成果を吸収）
- **Wgpu**: WGSL カーネル（collide+stream 融合、ping-pong）。phase9-wgpu の
  評価結果で採否・チューニング指針を確定。shader-f16 で FP16 格納
- **Cuda**（feature、後日・NVIDIA スパコン向け）: 同一 trait の追加実装。
  MPI+CUDA の GPUDirect は M-D 以降
- 境界条件パス（Zou-He 等）は「エッジセルのみの小カーネル」としてバックエンド毎に実装。
  数式は Lattice trait の面テーブルから生成し、CPU/GPU で同一定義を共有

### 2.5 移行戦略（R5: 2Dスイートを人質に取らない）— **全手順完了（2026-07-05）**

1. `lbm-core2` を新設し V2 抽象を実装（V1 は凍結・参照実装化）✅
2. **V1 API のファサード**を lbm-core2 上に実装（`Simulation<T>` の全公開メソッド）。
   既存 56+ テスト・wasm・CLI がそのまま通ることを最初のマイルストーンにする ✅
3. T13（分割不変）/ T14（バックエンド等価）を codex に発注、V2 の敵対検証を確立 ✅
4. 3D（D3Q19）を Lattice 追加として実装 → T15（3D物理）✅
5. 安定後、lbm-core2 → lbm-core に改名し V1 を削除 ✅
   （2026-07-05 実施。現在 `crates/lbm-core` が本文書の V2 実装そのもの。
   V1 との等価性証明は削除直前の `tests/v1_match.rs` ヘッダに凍結値として
   記録済み — ブランチ履歴参照。`compat` モジュールは公開 API として存続し、
   scenario / CLI / wasm の 2D 経路はこれを使う）

## 3. シナリオ契約 v1 への拡張（後方互換）

```jsonc
{
  "grid": { "nx": 256, "ny": 256, "nz": 128 },        // nz 追加（省略=2D）
  "physics": { "precision": "f32", "storage": "f16" }, // storage 追加
  "compute": {                                          // 新設（全て省略可）
    "backend": "auto | cpu | gpu",
    "decompose": { "ranks": [2, 2, 1] }                 // MPI 実行時のヒント
  },
  "outputs": [ { "field": "q-criterion", "format": "vtk" } ]  // 3D 可視化系を追加
}
```
- 既存フィールドは不変（deny_unknown_fields のまま追加のみ）
- GUI は 2D 専用のまま（3D はまず CLI/エージェント経由。断面/等値面 GUI は M-F 以降に判断）

## 4. 検証マップ（何を書けば「done」か）

| 軸 | テスト | 内容 |
|---|---|---|
| 分割 | T13 | 1×1 vs 2×2 vs 4×1 vs 1×4（+3D: 2×2×2）一致（f64 ≤1e-12） |
| バックエンド | T14 | CPU vs Wgpu 同一シナリオ f32 相対 ≤1e-5、診断値一致 |
| 3D 物理 | T15 | TGV3D 収束次数・球抗力（Re=100 Cd≈1.09 帯）・3D キャビティ・（LES後）チャネル Re_τ=180 |
| 精度 | T16 | f16 格納の劣化定量（TGV/キャビティで凍結帯） |
| 回帰 | 既存 T1-T12 | V1 ファサード経由で無修正緑 |

## 5. リスクと正直な見積り

- **wgpu の f64 なし** → f64 が要る検証は CPU バックエンド担当と割り切る（T14 は f32 比較）
- **MPI 実測の上限はローカル**（複数ランク機能検証まで）。真の弱スケーリングは
  クラスタアクセス取得後（COMPETITIVE_SPEC §5）
- **Esoteric-Pull 等の in-place ストリーミング**（メモリ半減）は FluidX3D 級を狙う
  M-E での導入候補。まず素直な ping-pong で正しさを固める
- 工数感（自律エージェント並列前提）: M-B コア V2 ≈ 2-4 晩、M-C 3D ≈ 2-3 晩、
  M-D MPI ≈ 2-3 晩 + クラスタ実測は別途
