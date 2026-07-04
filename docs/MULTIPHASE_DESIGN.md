# 混相流モジュール設計（Phase 4）

## 方針

Shan-Chen 疑似ポテンシャル法を採用する。理由:
- 実装が衝突・ストリーミングに対して**直交**（力場を注入するだけ）で、
  検証済みの単相コアを壊さない
- 界面追跡が不要（拡散界面）で初学者向けデモに強い
- 単成分多相（液滴・気泡・凝縮）と二成分混相（RT 不安定性・二相流）の
  両方に自然に拡張できる

制約（既知の弱点、docs に明記すること）:
- 密度比はクラシック ψ で ~50、CS-EOS 導入で ~100-1000（Phase 4 では前者から）
- 疑似速度（spurious currents）が界面近傍に O(1e-2) 出る
- 表面張力と密度比が G で連動する（独立制御は multi-range ψ が必要 → 将来）

## Phase 4a: 単成分多相（SCMP）

### エンジン変更（lbm-core）

1. **セル別力場**: `Simulation` に `force_field: Option<Vec<[T; 2]>>` を追加。
   - collide の Guo 項 / update_moments の F/2 補正で `F_local = force + force_field[i]`
   - 公開 API: `sim.set_force_field(Some(vec))` / `sim.force_field_mut() -> &mut [...]`
     （毎ステップ書き換える呼び出し側のためにアロケーション再利用）
2. 既存の一様 `force` はそのまま（重力用）。

### multiphase モジュール

```rust
pub struct ShanChen<T> {
    pub g: T,                 // 流体間相互作用強度（負で引力）
    pub g_wall: T,            // 壁付着強度（接触角制御）
    pub psi: Psi,             // ポテンシャル関数の選択
}
pub enum Psi { Classic /* 1 - exp(-rho) */, Exponential { rho0: f64 } }

impl ShanChen<T> {
    /// 現在の rho 場から SC 力場を計算して sim にセットする。
    /// 使い方: loop { sc.update_force(&mut sim); sim.step(); }
    pub fn update_force(&self, sim: &mut Simulation<T>);
}
```

- F(x) = −G ψ(x) Σ_q w_q ψ(x+c_q) c_q（流体隣接）
- 壁: F_ads(x) = −G_w ψ(x) Σ_q w_q s(x+c_q) c_q（s=1 if solid）
- 周期境界はラップ、開境界セルはゼロ勾配外挿（Phase 4 では SC と開境界の
  併用は非サポートでもよい — 仕様化する）

### 検証（T11 具体化）

- 平坦界面の共存密度: G=−5.0, ψ=Classic, 128×64 周期、上下半分 ρ_l/ρ_v 初期化
  → 定常密度が理論 Maxwell 構成（数値積分で別途計算した参照値）±3%
- Laplace 則: R ∈ {12, 16, 20, 24} の液滴、Δp = σ/R 線形フィット R² ≥ 0.99
- 疑似速度: max|u| ≤ 0.05（G=−5, τ=1）
- 接触角: G_w ∈ 掃引で θ ∈ {~60°, 90°, ~120°} を ±10°（測定は
  液滴の高さ/半径から球冠フィット）

## Phase 4b: 二成分混相（MCMP）

- `MultiComponent<T>`: 2 つの分布関数セット（内部的に Simulation を 2 つ持つか、
  f を [2][N*Q] で持つ専用構造体か → 実装時に決定。共有 solid/幾何は必須）
- 相互作用: F_σ = −G_AB ψ_σ(x) Σ w_q ψ_σ̄(x+c_q) c_q（σ̄ は相手成分）
- 共通速度 u' = (Σ_σ m_σ ω_σ + ...)/... （Shan-Chen 標準の合成速度で衝突）
- 検証（T12）: RT 不安定性の線形成長率 vs 理論（Atwood 数 0.5、±20%）、
  液滴分離・合体の質量保存

## 実装順

1. 力場 API（+単体テスト: 一様 force_field ≡ uniform force の一致）
2. SCMP + 平坦界面/Laplace（codex 検証発注）
3. 接触角
4. MCMP + RT（codex 検証発注）
