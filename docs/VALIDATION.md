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
sim.set_inlet_profile(Edge::Left, |c| [ux, uy]); // VelocityInlet エッジの節点別プロファイル
                                           // c はエッジ沿い座標（左右エッジ=y, 上下=x）
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
- 設定: left VelocityInlet + `set_inlet_profile` で放物線プロファイル
  （u_max=0.05、T2 の解析形、リム座標は [0,0]）、right PressureOutlet{rho:1}、
  上下 BounceBack。96×34（H=32）、TRT。
- 合格基準（定常後）:
  - **バルク領域**（流出境界から 24 列以上離れた断面）の質量流束
    Q(x)=Σ_y ρ·ux が一定: max|Q−Q̄|/Q̄ ≤ 1e-4（実測 2.4e-5）
  - 中央断面プロファイル L2rel ≤ 2e-3 vs 放物線
  - 定常後 10⁴ step の総質量ドリフト ≤ 1e-11（実測 2e-13）
- **既知アーティファクト（仕様）**: 圧力流出境界の直前 ~4 列に O(Ma²) の
  スタッガード振動（±2% 程度、衝突演算子非依存、減衰長 ~4 セル）が出る。
  これは Zou-He 圧力境界の固有特性であり不合格条件ではない（PHYSICS.md 参照）。
- 角度: 4 方向（左→右、右→左、下→上、上→下）すべてで同等の結果。

### T5. 圧力差駆動チャネル（Zou-He 圧力-圧力）
- 設定: left PressureOutlet{rho_in}, right PressureOutlet{rho_out}（Δρ 小、例 2e-3）、
  上下 BounceBack。解析: dp/dx = cs²Δρ/L による Poiseuille。
  **L = nx−1**（圧力指定ノードは境界列上にあり、その間隔が有効チャネル長）。
- 合格基準: 定常流量が解析値の ±2%（TRT, H=32。実測 0.26%）。
  圧力場 p(x)=cs²ρ(x) の線形性（両端 8 列を除くバルク）R² ≥ 0.999。
- 角度:
  - **厳密**: Δρ 符号反転 + x 鏡映で場が厳密に鏡映一致（L∞ ≤ 1e-12。
    離散系の x 反転対称性による）
  - **近似**: 鏡映なしの単純符号反転は慣性項・圧縮性が O(Ma²) で対称性を破るため
    相対 L∞ ≤ 5e-3 まで（実測 1.7e-3）。厳密 1e-12 を要求するのは物理的に誤り。

### T6. 保存則・整合性
- 周期箱 + 任意初期場: 総質量が 10⁴ step で相対 1e-11 以内に一定
  （丸め誤差の蓄積 ~1e-13/10³step を実測済みのため。物理的には厳密保存）。
- BB 箱（全辺壁）: 同上。
- 一様力 F の周期箱: 総運動量が 1 step あたり `N_fluid * F` ずつ増える（相対 1e-10）。
- feq の 0,1,2 次モーメント恒等式（単体テスト）: Σfeq=ρ, Σfeq c=ρu, Σfeq cc = ρ(cs²I+uu)
  （|u|≤0.1 の数点で 1e-14）。
- 角度（f32）: 質量ドリフト ≤ 1e-5（10³step）、力による運動量成長の相対誤差
  ≤ 1e-5（10²step）。**偏差格納方式（2026-07-05 導入）後の実測は 2.8e-7**
  （導入前は一様場のコヒーレント丸めバイアスで 1.3e-3 だった。PHYSICS.md 参照）。

### T7. リッド駆動キャビティ（Ghia et al. 1982 比較）
- 設定: 全辺壁・上辺 MovingWall{[U,0]}、Re = U*L/ν ∈ {100, 400, 1000}、N=129
  （L=N-2）、U=0.1、TRT。定常まで（ε=1e-8 か 300k step 上限）。
- 合格基準: 幾何中心線 u(y)・v(x) を Ghia 表の 17 点と比較し RMS 誤差 ≤ 0.02·U
  （Re=100/400）、≤ 0.03·U（Re=1000）。主渦中心位置が文献値 ±0.02L。
  **既知の誤植**: Re=400 の v(x=0.9063)=−0.23827 は流通データの既知の誤り
  （隣接点と不連続、PHYSICS.md 2026-07-05 参照）。この 1 点は RMS から除外する。
- 角度: 蓋の向きを 4 方向に回して同一解。**正しい対称写像**（PHYSICS.md 記載。
  左蓋 [0,−U] は反対角鏡映 p'=(N−1−y,N−1−x), v=(−uy',−ux') 等）を用いること。
  合格基準 L∞ ≤ 1e-10（実測は機械精度 ~4e-16）。Re=100・2000 step で可。

### T8. 円柱周り流れ — Schäfer–Turek ベンチマーク（力測定・渦放出）
確定参照値を持つ標準ベンチマーク（Schäfer & Turek 1996, "Benchmark computations
of laminar flow around a cylinder"）を採用する。幾何（比率厳守）:
チャネル 22D × 4.1D、円柱中心は流入から 2D・下壁から 2D（**わずかに非対称**、
これが渦放出のトリガー。中心 y/H = 0.4878）。上下 BounceBack、
left VelocityInlet + `set_inlet_profile` 放物線 u(y) = 4 u_max y_w(H−y_w)/H²、
right PressureOutlet{1.0}。U_mean = (2/3) u_max、Re = U_mean·D/ν。
Cd = 2Fx/(ρ U_mean² D)、Cl = 2Fy/(ρ U_mean² D)。

- **2D-1（Re=20, 定常）** 参照値: Cd = 5.5795, Cl = 0.0106, Δp* = Δp/(ρU_mean²) = 2.9375
  - D=20（格子 440×82、u_max=0.075, ν=0.05 → Re = 0.05·20/0.05 = 20）
    デフォルトスイート: Cd ∈ [5.2, 6.0]、Cl ∈ [−0.05, 0.08]（staircase 粗格子帯）
  - D=40（格子 880×164、u_max=0.075, ν=0.1）#[ignore]:
    Cd ∈ [5.35, 5.85]、収束傾向（|Cd(40)−5.5795| < |Cd(20)−5.5795|）
- **2D-2（Re=100, 非定常）** 参照値: Cd_max ≈ 3.22–3.24, Cl_max ≈ 0.99–1.01,
  St ≈ 0.295–0.305。D=40, u_max=0.15（U_mean=0.1, ν=0.04）#[ignore]:
  - St ∈ [0.28, 0.32]（Cl のゼロ交差から測定）
  - Cd_max ∈ [3.0, 3.5]、Cl_max ∈ [0.8, 1.2]
  - 渦放出の周期性: 連続する Cl 周期の長さのばらつき ≤ 2%
- 備考: staircase 近似のため帯は参照値より広い。曲面境界（Phase 7 候補）導入時に
  タイト化する。旧仕様（周期境界・非拘束帯との比較）は幾何不整合のため廃止
  （PHYSICS.md 2026-07-05）。

### T9. Outflow（ゼロ勾配）の健全性
- 設定: T8-2D-2 相当のチャネルで右辺を Outflow に置換。渦が流出面を通過しても
  発散しない。
- 合格基準: 10⁵ step で NaN/Inf なし、逆流質量流束が総流入の 5% 以下、
  流出面近傍 (x>0.9L) の圧力振動 rms が中央部の **15 倍以内**（実測 11.3。
  ゼロ勾配流出は圧力波を部分反射する固有特性。改善は convective outlet を
  Phase 7 バックログで検討、PHYSICS.md 2026-07-05）。

### T10. ロバスト性・エラーパス
- τ ≤ 0.5、Periodic 非ペア、Zou-He 直交エッジ違反、nx<3 → `ConfigError`。
- **安定限界ケース（パラメータ確定済み）**: τ=0.51, N=128, U=0.05（Re≈1890）,
  TRT Λ=3/16 のキャビティが 10⁴ step NaN/Inf なし（実測 max|u|=0.046 で安定）。
  U=0.1（Re≈3780）は Λ=3/16, 1/4 とも ~3.5-7k step で発散する（既知の限界、
  グリッドレイノルズ数 U/ν ≈ 30 は超過。ガイドライン: τ→0.5 では U/ν ≤ 15）。
- `set_solid` を開境界エッジ上に置くと panic（仕様）。
- 移動壁/流入速度が |u| > MAX_SPEED(=0.3) なら `ConfigError::VelocityTooHigh`。
  `set_inlet_profile` の速度超過は panic。

### T11. Shan-Chen 単成分多相（Phase 4a・実測校正済み 2026-07-05）
共通設定: `ShanChen::new(-5.0)`（classic ψ = 1−e^{−ρ}）、τ=1（nu=1/6）、
初期化 液 ρ=2.0 / 蒸気 ρ=0.15、毎 step `sc.update_force(&mut sim); sim.step()`。
圧力は**必ず SC EOS**（`sc.pressure(rho)` = cs²ρ + (G cs²/2)ψ²）で比較する。

- **平坦界面**（64×128 周期、30k step）:
  - 共存密度 ρ_l = 1.888 ± 2%、ρ_v = 0.1194 ± 3%（実測回帰値）
  - 相間圧力平衡: |p_l − p_v|/p ≤ 1e-4（実測 8.5e-6）
  - 疑似速度 max|u| ≤ 5e-3（実測 1.26e-3）
  - 総質量ドリフト ≤ 1e-10 相対（SC 力は質量 0 次モーメントを持たない）
- **Laplace 則**（128²、R₀ ∈ {12,16,20,24}、40k step）:
  - Δp vs 1/R_fit の線形性 R² ≥ 0.999（実測 0.99988）
  - 傾き σ = 3.32e-2 ± 10%（実測回帰値）、各液滴の σ=Δp·R が傾きと ±5%
  - 半径測定は密度中央値の等値線面積から（area/π の平方根）
- **f32 角度**: 平坦界面ケースが f32 でも安定（NaN なし・共存密度 ±5%）

### T11b. 接触角（G_w 特性の凍結）
- 壁付き液滴（下辺 BounceBack、他 Periodic）で G_w ∈ {−1.5, 0, +1.5} を測定。
  本実装は solid の ψ=0（cohesion から除外）+ 別項 −G_w ψ Σw s c のため
  **G_w=0 は 90° にならない**（非湿潤側に寄る）。テストは:
  - θ(G_w) が単調（G_w が負に大きいほど湿潤 = θ 小）
  - 3 点の実測角を回帰凍結（±8°）し、測定法（球冠フィット: 接触幅と高さから
    θ = 2·atan(2h/w)）をコメントで文書化
- 将来: 仮想壁密度方式（solid に ψ(ρ_w) を与える）への切替を検討（Phase 7）。

### T12.（Phase 4b・二成分 MCMP 導入後）Rayleigh–Taylor
- RT 不安定性: 初期擾乱波長 λ の成長率が線形理論と ±20% で一致（Atwood 数 0.5）。
- 前提となる二成分 Shan-Chen（MCMP）は初回総合レビュー後に実装予定
  （PLAN.md 2026-07-05 のスコープ再編を参照）。接触角は T11b に移動済み。

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
