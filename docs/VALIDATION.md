# VALIDATION.md — 境界条件・物理検証テスト仕様マトリクス

このファイルは検証テストスイートの**発注仕様**である。テスト作者（codex）は
この仕様と公開 API のみを根拠にテストを書くこと（エンジン内部実装を写さない）。
受入基準は f64・`--release` 実行を前提とする。

## 公開 API（テストから使うもの）

```rust
use lbm_core::prelude::*;

let mut sim: Simulation<f64> = SimConfig {
    nx: 64, ny: 64,
    nu: 0.02,                        // 動粘性係数（格子単位）。tau = 3*nu + 0.5
    collision: Collision::Trt { magic: 0.1875 },   // または Collision::Bgk
    edges: Edges {
        left:   EdgeBC::Periodic,
        right:  EdgeBC::Periodic,
        bottom: EdgeBC::BounceBack,               // 静止壁（1セルのソリッドリム）
        top:    EdgeBC::MovingWall { u: [0.1, 0.0] },
    },
    force: [1e-6, 0.0],              // 一様体積力（Guo）
    ..Default::default()
}.build().unwrap();                  // 不正構成は Err(ConfigError)

sim.set_solid(x, y);                       // 内部障害物（リム外の任意セル）
sim.set_solid_region(|x, y| bool);         // 述語で一括指定
sim.init_with(|x, y| (rho, ux, uy));       // f = feq(rho,u) で初期化（未呼出時 rho=1, u=0）
sim.step();                                // 1 タイムステップ
sim.run(n);                                // n ステップ

sim.nx(); sim.ny(); sim.time();            // 形状・経過ステップ
sim.rho(x, y); sim.ux(x, y); sim.uy(x, y); // 巨視量（force 半補正込みの物理速度）
sim.rho_field(); sim.ux_field(); sim.uy_field(); // &[T]（cell = y*nx + x）
sim.is_solid(x, y);
sim.total_mass(); sim.total_momentum();    // Σρ, [Σρux, Σρuy]（流体セルのみ）
sim.set_force_probe(|x, y| bool);          // momentum-exchange 力測定の対象ソリッド集合
sim.probed_force();                        // 直近 step の [Fx, Fy]
sim.fluid_cell_count();                    // 非ソリッドセル数（運動量テスト用）
```

- エッジ種: `Periodic` / `BounceBack` / `MovingWall{u}` / `VelocityInlet{u}`（Zou-He）/
  `PressureOutlet{rho}`（Zou-He）/ `Outflow`（ゼロ勾配）
- 壁エッジはソリッドリム（1セル）として実現される。**壁面はリムセル中心と隣接流体セル
  中心の中間**（half-way）に位置する。例: 上下壁・格子高さ `Ny` のとき、リムは y=0 と
  y=Ny-1、流体行は y=1..=Ny-2、壁面は y=0.5 と y=Ny-1.5。したがって
  **チャネル幅 H = Ny-2（= 流体セル行数）**、流体セル中心の壁面からの距離は
  y_w = j - 0.5（j = 1..H）。Poiseuille 最大速度は g·H²/(8ν)。
- 構築時エラー（`ConfigError`）: tau ≤ 0.5（nu ≤ 0）、Periodic の非ペア、
  Zou-He/Outflow エッジの直交エッジが壁/Periodic 以外、nx/ny < 3 など。

## 記法

- 誤差ノルム: `L2rel(u, u_ref) = sqrt(Σ|u-u_ref|²) / sqrt(Σ|u_ref|²)`（流体セルのみ）
- 収束次数: `order = log2(err(N) / err(2N))`
- 定常判定: `max|u^{t+Δ} - u^t| / max|u| < ε`（Δ=500 step）。**ε = 1e-11 を推奨**:
  BGK は丸め誤差プラトー ~1e-12 で恒久振動するため 1e-13 は到達不能
  （docs/PHYSICS.md の実験記録参照）。

---

## テストマトリクス

### T1. Taylor–Green 渦（周期境界・粘性減衰・収束次数）
- 設定: 全辺 Periodic、N×N（N=32, 64）、ν=0.02、**拡散スケーリング u0 = 1.28/N**、
  k=2π/N。解析解 `ux = -u0 cos(kx) sin(ky) e^{-2νk²t}`, `uy = +u0 sin(kx) cos(ky) e^{-2νk²t}`。
  初期化は `init_with` で**圧力整合密度 ρ = 1 − (3u0²/4)(cos 2kx + cos 2ky)** を渡すこと
  （一様 ρ=1 だと音波残留で O(u0) 汚染される。docs/PHYSICS.md 参照）。
- 合格基準:
  - t = 1/(2νk²) 経過時点で速度場 L2rel ≤ 1.5e-3（N=64, TRT。実測 7.0e-4）
  - 収束次数 order ≥ 1.7（N=32→64。実測 1.91）
  - 減衰率フィットから実効粘性 ν_eff が公称 ν の ±2% 以内（N=64）
- 角度: BGK と TRT で同等、90°回転した初期場で結果が回転対称（L∞ ≤ 1e-12）。

### T2. 体積力駆動 Poiseuille 流（half-way BB の厳密性）
- 設定: 上下 BounceBack、左右 Periodic、力 F=[g,0]（g=1e-6 など）、ny 任意（H=ny-2）。
  解析解: `ux(y) = g/(2ν) * y_w (H - y_w)`、y_w = (セル中心の壁面からの距離) = j-0.5。
- 合格基準:
  - TRT（Λ=3/16）: 定常（ε=1e-11）で L∞rel ≤ 1e-10（**厳密**。H=8 でも成立）
  - BGK: H=8→16 で収束次数 ≥ 1.7（BGK は τ 依存スリップ誤差があるため厳密は要求しない）
  - プロファイルの上下対称性: |ux(j) - ux(H+1-j)| ≤ 1e-13
- 角度: 同じ設定を 90°回転（左右壁・F=[0,g]）しても同一プロファイル。

### T3. Couette 流(移動壁)
- 設定: 上 MovingWall{u:[U,0]}（U=0.1）、下 BounceBack、左右 Periodic。
  解析解: `ux(y_w) = U * y_w / H`（壁位置 half-way 基準）。
- 合格基準: 定常で L∞rel ≤ 1e-10（BGK/TRT・τ∈{0.6, 1.0, 1.4} すべて）。
- 角度: 下壁を動かす/左右壁で縦 Couette にしても同等。移動壁の質量保存（総質量ドリフト
  ≤ 1e-12 相対 / 10⁴ step）。

### T4. Zou-He 速度流入 + 圧力流出チャネル
- 設定: left VelocityInlet{放物線 u(y)}, right PressureOutlet{rho:1}, 上下 BounceBack。
  流入プロファイルは T2 の解析形（最大速度 0.05 程度）。
- 合格基準:
  - 定常で全断面の流量 Q(x) が一定（max|Q(x)-Q̄|/Q̄ ≤ 1e-6）
  - 中央断面プロファイル L2rel ≤ 2e-3（TRT, H=32）
  - 質量保存: 流入流束 = 流出流束（相対差 ≤ 1e-6）
- 角度: 4 方向（左→右、右→左、下→上、上→下）すべてで同等の結果。

### T5. 圧力差駆動チャネル（Zou-He 圧力-圧力）
- 設定: left PressureOutlet{rho_in}, right PressureOutlet{rho_out}（Δρ 小、例 2e-3）、
  上下 BounceBack。解析: dp/dx = cs²Δρ/L による Poiseuille。
  **L = nx−1**（圧力指定ノードは境界列上にあり、その間隔が有効チャネル長）。
- 合格基準: 定常流量が解析値の ±2%（TRT, H=32）。圧力場 p(x)=cs²ρ(x) の線形性
  R² ≥ 0.999。
- 角度: Δρ の符号反転で流れが正確に反転（ux 場の符号反転一致 L∞ ≤ 1e-12）。

### T6. 保存則・整合性
- 周期箱 + 任意初期場: 総質量が 10⁴ step で相対 1e-11 以内に一定
  （丸め誤差の蓄積 ~1e-13/10³step を実測済みのため。物理的には厳密保存）。
- BB 箱（全辺壁）: 同上。
- 一様力 F の周期箱: 総運動量が 1 step あたり `N_fluid * F` ずつ増える（相対 1e-10）。
- feq の 0,1,2 次モーメント恒等式（単体テスト）: Σfeq=ρ, Σfeq c=ρu, Σfeq cc = ρ(cs²I+uu)
  （|u|≤0.1 の数点で 1e-14）。
- 角度: f32 でも成立（緩めた許容 1e-4 相対）。

### T7. リッド駆動キャビティ（Ghia et al. 1982 比較）
- 設定: 全辺壁・上辺 MovingWall{[U,0]}、Re = U*L/ν ∈ {100, 400, 1000}、N=129
  （L=N-2）、U=0.1、TRT。定常まで（ε=1e-8 か 200k step 上限）。
- 合格基準: 幾何中心線 u(y)・v(x) を Ghia 表の 17 点と比較し RMS 誤差 ≤ 0.02·U
  （Re=100/400）、≤ 0.03·U（Re=1000）。主渦中心位置が文献値 ±0.02L。
- 角度: 蓋の向きを 4 方向に回して同一解（対称変換後 L∞ ≤ 1e-10、Re=100 のみで可）。

### T8. 円柱周り流れ（力測定・非定常渦放出）
- 設定: チャネル内円柱。D=20〜24、領域 ~(22D)×(8D)+ 上下 Periodic
  または壁、left VelocityInlet{一様 U=0.05〜0.1}, right Outflow または PressureOutlet。
  Re=UD/ν。momentum-exchange で円柱の力を毎 step 測定。
- 合格基準:
  - Re=20（定常）: Cd = 2Fx/(ρU²D) ∈ [1.8, 2.4]（文献 ~2.0-2.2、ブロッケージ込み許容）
  - Re=100（非定常）: Strouhal St = f·D/U ∈ [0.15, 0.19]、Cd ∈ [1.2, 1.5]、
    Cl 振幅 ∈ [0.2, 0.45]
  - 渦放出の周期性: Cl(t) の FFT ピークが明瞭（ピーク/背景 ≥ 10）
- 備考: staircase 近似のため許容帯は文献値より広い。将来（曲面境界導入時）狭める。

### T9. Outflow（ゼロ勾配）の健全性
- 設定: T8 と同じチャネルで右辺 Outflow。渦が流出面を通過しても発散しない。
- 合格基準: 10⁵ step で NaN/Inf なし、逆流質量流束が総流入の 5% 以下、
  流出面近傍 (x>0.9L) の圧力振動 rms が中央部の 3 倍以内。

### T10. ロバスト性・エラーパス
- τ ≤ 0.5、Periodic 非ペア、Zou-He 直交エッジ違反、nx<3 → `ConfigError`。
- τ=0.51・Re 高めのキャビティ（TRT）で 10⁴ step NaN なし。
- `set_solid` を Zou-He エッジ上に置いた場合の挙動が仕様どおり（エラー or 無視を仕様化）。
- u 指定が音速制限超（|u|>0.3）で警告 or エラー（仕様化されていること）。

### T11.（Phase 4 で有効化）Shan-Chen 多相: Laplace 則
- 液滴半径 R を変えて Δp = σ/R の線形性（R² ≥ 0.99）、σ の再現性 ±5%。
- 平坦界面の共存密度が Maxwell 構成 ±3%（EOS: Carnahan-Starling 等を採用時）。
- 疑似速度（spurious currents）max|u| ≤ 0.05（採用モデルの標準水準を仕様化）。

### T12.（Phase 4）接触角 / Rayleigh–Taylor
- 壁付着力パラメータ掃引で接触角 30°/90°/150° を ±10° で実現。
- RT 不安定性: 初期擾乱波長 λ の成長率が線形理論と ±20% で一致（Atwood 数 0.5）。

---

## テスト実装の規約（codex 向け）

- 置き場所: `crates/lbm-core/tests/validation_*.rs`（1 テーマ 1 ファイル）。
- 共有ヘルパは `crates/lbm-core/tests/common/mod.rs`。
- 重い計算（T7 の Re=1000, T8, T9）は `#[ignore]` を付け、CI 相当は
  `cargo test --release`、フルは `cargo test --release -- --include-ignored`。
- 乱数不使用（決定論）。`assert!` には実測値を含むメッセージを付ける
  （例: `assert!(err < 5e-3, "L2rel = {err}")`）。
- 外部クレート追加は `approx` のみ可（それ以外は要相談）。
- Ghia 参照データはテストファイル内に定数表として埋め込む（出典コメント付き）。
