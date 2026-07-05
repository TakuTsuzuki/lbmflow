# 要求定義書（完成版）：回転境界・高密度比二相・LES 連成 3D マルチフィジックス LBM ソルバ

**文書ID**: REQ-M-F-STR / **改訂**: rev.3（rev.1=codex 敵対的レビュー全48件反映／rev.1a=PM 決定「既定は忠実度最優先・緩和は後付け拡張点」を反映／rev.1b=PM 統合: 表題のドメイン中立化・コア改名追随・§7 メモリ予算表追加・VALIDATION T17 配線／rev.2=codex 第2次レビュー 11 件全採択: 緩和同等性検証 VR-STR-RELAX 新設・スコープ語の明確化・変数 σ 表面張力規約・F_b^scalar 追加・予算表算術修正ほか — docs/proposals/req-round2-findings.md 参照／**rev.3** = competitive-review triage diff merged (authored as "rev.1c" against rev.1b, layered here onto rev.2): P1 population balance, P2 §4.8 extension contracts, P3 FR-IO-05/06, P4 reference datasets, P5 product-layer scope note, §11 implementation dependency DAG. New content in English per the 2026-07-05 language directive.／**rev.4** = external numerical-physics review (REV-CFD-*, filed against rev.1a) triaged and merged: all 4 Critical fixes (sparger phase inversion, Allen–Cahn/continuity mass-flux consistency, non-equilibrium stress stage convention, forcing second-moment single-definition), Ca_spurious dimensional fix, Pe_N/Pe_tip split, active-scalar predictor–corrector dataflow, precision-profile enum, conservative scalar forms, four-way contact contract (extension-gated), viscosity-interpolation & σ(κ,β) freeze, ε_g processing definitions, and **provisional numeric acceptance bands** under a band-governance rule. Disposition table in TESTING_NOTES.md; MJ-008 was already fixed in rev.2. This file is PM-owned for translation — excluded from the bulk translation session.）
**位置づけ**: `docs/PLAN.md` の M-F（垂直機能）／ `ARCHITECTURE_V2.md` への上位要求。
検証受入は [VALIDATION.md](VALIDATION.md) **T17**（VR-STR-01〜07 を配線）。
**対象コア**: `lbm-core`（旧 lbm-core2。D3Q19/D3Q27, CpuScalar/CpuSimd/wgpu, MPI 分割）
**代表適用問題**: 撹拌槽反応器（機能要求はドメイン中立に定義し、§2 と §8 の検証ベンチが
この適用を具体化する。M-Star CFD の中核ユースケースに正面から重なる領域）

## 0. レビュー反映方針

- Critical/Major の式・符号・定義バグ（非平衡応力評価、`τ_eff`、MRF 見かけ力、表面張力、保存型 Allen-Cahn、Np/N_Q、無次元数）を全て修正・係数固定した。
- 欠落物理（上部脱気、静水圧 well-balanced、乱流 Schmidt 数 SGS スカラー、粒子 SGS 分散、気泡誘起乱流 BIT、接触角、初期化・助走、メモリ予算）を新規要求として追加した。
- **納品スコープ（rev.2 で明確化、rev.4 で列挙を確定）**: **忠実度既定のサブシステム群（§1 表の「既定」列）を一括同時実装**する。緩和拡張（MRF・point-bubble・one-way・block-AMR・積極的 f32）は初版では **trait 境界・設定スキーマ・検証項目（VR-STR-RELAX）を予約するのみ**とし、実装は後付け。codex の「モード分割」は納品フェーズではなく **実行時モードの相互排他制約** として §1 の構成マトリクスに落とす。同一計算内で物理的に競合するモード（MRF＋IBM 同一ゾーン、位相平均＋MRF 等）はコンフィグ検証で棄却する。
  **Initial delivery (release gate)**: IBM-inertial, resolved-phasefield, active scalar
  (predictor–corrector coupling, §5), two-way particles, uniform grid, fidelity
  precision profile (`mixed_safe`, §7). **Phase 2 (API-reserved now, implemented
  later, accepted via VR-STR-RELAX)**: MRF-frozen-rotor, point-bubble (+PBM),
  one-way particles, four-way particle contact (FR-PART-04..06), block-AMR,
  aggressive f32 (`mixed_fast`), hybrid interface, thermal axis.
  "All subsystems simultaneously" means all physics axes exist in the initial
  delivery — not all modes of each axis. (REV-CFD-MJ-008; same substance as the
  rev.2 fix, list made explicit.)
- 検証不能語（「保証」「自然に許容」「安定に積分」）を測定可能な誤差・保存則・適用範囲へ置換した。
- AMR は「初期版は一様格子基準、AMR は上位オプション（coarse-fine 保存補間・時間ステップ比・検証問題を定義した上で有効化）」へ降格（#29）。

---

## 1. 実行時モード構成マトリクス（相互排他）／設計原則

**設計原則（忠実度既定・緩和は後付け拡張点）**: 各モード軸を strategy/trait で抽象化する。**既定は各軸とも忠実度最優先の実装（＝基準解）** とし、低コスト近似（MRF・point-bubble・one-way・AMR・積極的 f32）は同一 trait 背後の**後付け拡張**として追加する。コア連成ループ（§5）を変更せずに差し替え可能な構造とし、緩和モードは対応する忠実度基準解に対する許容誤差（§8 VR に閾値定義）で検証する。忠実度既定は計算コスト最大（§7 メモリ予算は忠実度構成でサイジングする）。

全モードを実装するが、1 計算あたり各軸から排他的に1つを選択。コンフィグ検証層が不整合な組合せを棄却する。

| 軸 | 既定（忠実度最優先） | 緩和拡張（後日）／基準級 | 排他制約・備考 |
|---|---|---|---|
| 回転 | `IBM-inertial`（非定常・時間精度） | 緩和: `MRF-frozen-rotor`（定常近似）／基準級: `sliding-overset` | MRF は IBM 移動翼と併用不可（#6）。位相平均統計は IBM/overset のみ（#37）。 |
| 界面 | `resolved-phasefield`（保存型 Allen-Cahn。界面・物質移動忠実度優先） | 緩和: `point-bubble`（Euler-Lagrange）／`hybrid` | 切替判定は `d_b/Δx, d_b/W, Eo, Re_b, α_g, We_b`（#12）。hybrid は相間 質量・運動量・スカラー保存則を定義（§5）。 |
| スカラー | `active`（σ・粘性・密度・[温度] への帰還を有効） | 緩和: `passive`（帰還オプトアウト） | active の帰還対象と安定化を明示（#13）。 |
| 粒子結合 | `two-way`（高 `α_p` で `four-way`） | 緩和: `one-way` | `α_p`／mass-loading 閾値（#16）。反力散布カーネルと運動量保存検証を伴う。 |
| 精度 | 忠実度プロファイル: 界面近傍・保存量・トルク・界面曲率・縮約は `f64`、遠方バルクのみ `f32` | 緩和: 積極的 `f32`／基準級: 全 `f64` | #32, §7。 |
| 格子 | `uniform`（要求解像度で完全解像） | 緩和: `block-AMR`（coarse-fine 保存補間・時間ステップ比・検証を満たす場合） | #29。AMR は実装リスクゆえ後付け。 |

---

## 2. 代表適用問題・代表量・無次元数（撹拌槽反応器）

本節は §8 検証ベンチの**代表適用**を具体化するものであり、§4 の機能要求は
回転境界・高密度比二相・LES・スカラー/粒子連成一般に適用する（表題の中立化と整合）。

3D 円筒（または角柱）容器。連続相（Newton/非 Newton 液）、下部スパージャからの分散気相（`ρ_l/ρ_g ≈ 10³`, `μ_l/μ_g ≈ 10²`）、一定角速度 `Ω` 剛体回転翼、中立浮力近傍の懸濁粒子、界面物質移動と液相反応を伴う複数スカラー。目的観測量: 時間／位相平均 3D 速度場、ひずみ速度第二不変量に基づくせん断応力場、粒子ラグランジュ累積せん断曝露、ガスホールドアップ `ε_g` と溶存スカラー濃度場。

### 2.1 無次元数定義（代表量を固定, #26, #23, #24）

代表回転数 `N = Ω/(2π)` [rev/s]、インペラ径 `D`、槽径 `T`、液深 `H`、重力 `g`、気泡径 `d_b`、粒子径 `d_p`、分子拡散 `D_m`、表面張力 `σ`、`Δρ = ρ_l − ρ_g`。

```
Re   = ρ_l N D² / μ_l                 (撹拌レイノルズ, 代表速度 U_tip = πND)
Fr   = N² D / g
We   = ρ_l N² D³ / σ
Eo   = Δρ g d_b² / σ                  (=Bond, 気泡スケール)
Mo   = g μ_l⁴ Δρ / (ρ_l² σ³)
Ca   = μ_l U / σ
Sc   = ν_l / D_m
Pe_N   = N D² / D_m = Re·Sc          (impeller velocity scale ND — REV-CFD-MJ-006)
Pe_tip = U_tip D / D_m = π·Re·Sc     (tip speed U_tip = πND; each use site must
                                      state which Pe it means — no bare "Pe")
Da_n = k C_ref^{n-1} · (L/U)          (反応 n 次, k は速度定数; 次数ごとに別記)
St   = τ_p / τ_f,   τ_p = ρ_p d_p² / (18 μ_l)
Np   = P / (ρ_l N³ D⁵),  P = 2π N T_q = Ω T_q   (T_q=トルク; N は rev/s, ρ は液相基準)
N_Q  = Q / (N D³)                     (Q=翼吐出面での正味体積流量)
```

格子側制約: `Ma_lattice = U_tip/c_s ≤ 0.1`、Cahn 数 `Cn = W/L`、界面 Péclet `Pe_φ = U W / M`、緩和時間 `τ ∈ [τ_min, τ_max]`。

### 2.2 マッチング優先順位（同時一致不能時, #25）

物理→格子の変換自由度は有限。全無次元数を同時一致できない場合の優先順位を固定する:
**(1) Re → (2) 密度比・粘性比＋We/Eo（界面動力学）→ (3) Fr（自由液面・浮力支配時）→ (4) Sc/Pe・Da（スカラー・反応）→ (5) St（粒子）**。
単位変換層は feasibility check を必須実行し、`Ma>0.1` / `τ∉[τ_min,τ_max]` / `Cn` 過大 / 拡散数・CFL 違反時は妥協した無次元数と誤差を明示して警告する。

---

## 3. 支配方程式系（修正版）

```
連続相（低 Mach LBM で回収, 相別密度で well-balanced 重力。
rev.4 / REV-CFD-CR-002: mass-flux consistency with the phase-field diffusion —
with ρ=ρ(φ) and a diffusive phase flux J_φ, the naive ∂ρ/∂t+∇·(ρu)=0 cannot hold;
the density flux J_ρ = (ρ_l−ρ_g) J_φ must appear in BOTH the continuity identity
and the momentum advection (consistent/AGG-type formulation — mandatory at
ρ_l/ρ_g ≈ 10³)）:
  ∂ρ/∂t + ∇·(ρu + J_ρ) = 0,        J_ρ = (ρ_l−ρ_g) J_φ
  ∂(ρu)/∂t + ∇·[(ρu + J_ρ) u] = -∇p + ∇·[ (μ(γ̇)+μ_t)(∇u+∇uᵀ) ]
                        + F_s + ρ g + F_b^{scalar} + F_g^{disp} + F_p + F_rot
  ・The SAME discrete J_ρ is used in both equations (single code path — verified
    by code review and the advected-droplet conservation test, §8 VR-STR-03/05).
  ・If a quasi-incompressible / pressure-evolution formulation is adopted instead,
    its continuity statement, divergence condition, and conservation checks must
    replace the above explicitly — silence is not an option.
  ・重力は全相に ρg を課し、静水圧 ∇p_hydro = ρg を well-balanced に離散化（#34）。
  ・F_b^{scalar}（rev.2, active 密度帰還）: 溶質浮力 F_b = ρ_0 β_C (C−C_0) g の
    Boussinesq 摂動力。C≡C_0 で厳密に 0 とし、ρ(φ)g の well-balanced 静水圧相殺とは
    **混ぜない**（独立の力源として合成。詳細は docs/proposals/active-scalar-feedback.md）。
  ・F_rot は MRF モードのみ（§4.3）。μ(γ̇) と μ_t の合成は §4.7 の陰的整合。

二相界面（保存型 Allen–Cahn 相場, Fakhari 2017 系に固定, #8。
rev.4: written in explicit conservative-flux form so J_φ is a first-class object）:
  ∂φ/∂t + ∇·(φu + J_φ) = 0,   J_φ = −M [ ∇φ − (4/W) φ(1−φ) n̂ ]
  n̂ = ∇φ / (|∇φ| + ε),  φ∈[0,1] (φ=1: liquid, φ=0: gas),  M[length²/time]
  Density interpolation: ρ(φ) = ρ_g + φ(ρ_l−ρ_g).
  Viscosity interpolation (rev.4 / REV-CFD-MJ-013 — default frozen):
    **harmonic in μ**:  1/μ(φ) = φ/μ_l + (1−φ)/μ_g
  Alternatives (linear-in-μ, linear-in-ν) are explicit config options, logged in
  run metadata, and are NOT covered by the default validation bands.

表面張力（化学ポテンシャル形式を基準に、σ 可変時は規約分岐。rev.2, #7）:
  μ_φ = 4β φ(φ−1)(φ−1/2) − κ ∇²φ
  σ = √(2κβ)/6,   W = 4√(κ/(2β))    ← these ARE the definitions for the adopted
  double-well free energy (rev.4 / REV-CFD-MJ-013: the "coefficients are
  model-defined" hedge is removed; internal consistency of {σ, W, κ, β, μ_φ} was
  verified in the codex round-2 review; the (σ,W) ↔ (κ,β) inversion is unique)
  ・σ 一定時（基準形）: F_s = μ_φ ∇φ （CSF 等価 σκn̂δ_s は検証項へ）
  ・σ が C_k / 温度に依存する active 時: F_s = μ_φ∇φ は直接使用せず、Marangoni 接線力
    との二重計上を避ける well-balanced CSF/化学ポテンシャル併用形に一本化する
    （docs/proposals/active-scalar-feedback.md 規約 D1/D2）。係数は本節の (κ,β,W,σ)
    規約へ導出してから凍結（要導出 — 実装前必須）。∇σ=0 への退化で σ 一定基準形と
    一致する退化テストを §8 に置く。

分散気相（point-bubble モード, Euler-Lagrange, #12）:
  m_b dv_b/dt = F_buoy + F_drag(Tomiyama) + F_lift + F_addedmass + F_walllub + F_TD
  BIT（気泡誘起乱流）生成項を LES に加算（§4.2, #46）

分散粒子（Euler-Lagrange, #16）:
  m_p dv_p/dt = F_drag(Schiller-Naumann, Re_p 範囲明示) + F_buoy
               + [高精度時] F_Saffman + F_Basset + F_Faxen
  two/four-way 時は反力を正則化カーネルで散布し運動量保存を検証

スカラー・反応（成分 k, active/passive 明示, #13, #14。
rev.4 / REV-CFD-MJ-011: conservative forms are normative for two-phase and
active-density cases — the non-conservative single-phase form is a special case）:
  Single-phase passive (simplified form, valid only when ρ, α uniform):
    ∂C_k/∂t + u·∇C_k = ∇·[ (D_k + ν_t/Sc_t) ∇C_k ] + R_k(C) + Ṡ_k^{if}
  Two-phase, phase-wise conservative (normative for gas–liquid scalars,
  q ∈ {gas, liquid}, α_liq = φ, α_gas = 1−φ):
    ∂(α_q C_{k,q})/∂t + ∇·(α_q u C_{k,q})
      = ∇·[ α_q (D_{k,q} + ν_t/Sc_t) ∇C_{k,q} ] + α_q R_{k,q}(C) + S_{k,q}^{if}
  Density-based active scalar (when the scalar feeds back into ρ):
    ∂(ρY_k)/∂t + ∇·(ρ u Y_k + J_k) = R_k + S_k^{if}
  ・SGS スカラー流束は乱流 Schmidt 数 Sc_t で閉じる（既定 Sc_t=0.7, 可変）。
  ・S^{if}: 解像界面は法線ジャンプ＋分配係数（Henry partition — the sign convention
    is: S_{k,liq}^{if} = −S_{k,gas}^{if}, interfacial flux positive into liquid）、
    point-bubble は k_L a(C*−C)（#35）。
  ・Conservation statement: Σ_q ∫ α_q C_{k,q} dV changes only by boundary fluxes
    and reactions — this is the quantity tested in VR-STR-05 scalar drift.

渦粘性（SGS, Smagorinsky と WALE を分離, #4）:
  Smagorinsky: ν_t = (C_s Δ)² |S̄|,  |S̄|=√(2 S̄:S̄)
  WALE:        ν_t = (C_w Δ)² (S^d:S^d)^{3/2} / [ (S̄:S̄)^{5/2} + (S^d:S^d)^{5/4} ]
               S^d は速度勾配テンソル二乗の deviatoric symmetric part（局所勾配復元が必要）
```

---

## 4. 機能要求（数値手法, 係数固定版）

### 4.1 基盤 LBM コア
- **FR-CORE-01**: D3Q19/D3Q27 選択可。**D3Q27 を既定とする条件を「多相 or 強 forcing or cumulant 使用時」と限定**（#30）。各格子で保持する平衡分布の Hermite 次数と回収精度を分けて定義（D3Q19 は 3 次等方性に制限あり）。**M-F 忠実度既定シナリオは多相・強 forcing 条件に該当するため常に D3Q27**。単相・弱 forcing の派生シナリオ（例: VR-STR-01 単相撹拌）では D3Q19 を許可（rev.2, 所見11）。
- **FR-CORE-02**: 中心モーメント（cascaded）／cumulant を実装。安定性は「保証」ではなく **対象ベンチでの許容緩和率範囲・positivity・regularization/filtering/entropic limiter の有無** で規定（#31）。
- **FR-CORE-03**: Guo forcing。速度モーメントは `ρu = Σ c_i f_i + Δt F/2`。
  Stress evaluation uses the forcing second-moment correction **as defined by the
  single equation in FR-STRESS-01** — prose words like "subtract"/"add" are banned
  from this topic; the equation is the only definition (rev.4 / REV-CFD-CR-004).
- **FR-CORE-04**: `Ma_lattice ≤ 0.1`、圧縮性誤差 `O(Ma²)` 制御。音響スケーリングと非圧縮性の整合を単位変換 feasibility に含める（#25）。

### 4.2 乱流モデル（LES-LBM）
- **FR-LES-01**: Smagorinsky（動的 Germano 含む）と WALE を**別式で**実装。WALE 既定（壁近傍 `ν_t→0`）。**WALE は速度勾配全体を要するため「有限差分不使用」要求は撤回**し、局所勾配の復元法（モーメント or コンパクト差分）を明示（#4）。
- **FR-LES-02**: 渦粘性の緩和時間反映は **`τ_eff = 1/2 + (ν_0+ν_t)/(c_s²Δt)`**（一般式）。格子単位 `c_s²=1/3, Δt=1` の簡略形 `Δτ_t = 3ν_t` を別記（#3 修正）。
- **FR-LES-03**: 壁せん断支配域は `y⁺` 壁関数 or 壁適合内挿境界。`τ_eff` は下限 `>1/2` だけでなく **上限クリッピングと診断**を設ける（過拡散・境界精度劣化回避, #27）。
- **FR-LES-04**: SGS スカラー流束（乱流 Schmidt 数 `Sc_t`）と SGS 熱流束（乱流 Prandtl 数）を ADE-LBM 緩和時間へ反映（#14）。

### 4.3 回転インペラ（モード排他, #5, #6, #21, #22）
- **FR-ROT-01**（IBM-inertial）: direct-forcing IBM（Uhlmann 型）。目標剛体速度 `U=Ω×r`。**「ガリレイ不変を保証」は削除**し、Taylor-Couette・回転円柱・移動壁 Couette での **すべり速度・トルク誤差・運動量保存誤差に閾値**を設定（multi-direct-forcing / implicit IBM 採用条件を含む）。
- **FR-ROT-02**（MRF-frozen-rotor）: 回転ゾーン内で **相対速度 `u_rel = u_abs − Ω×r` を解き**、体積力に コリオリ `−2ρ Ω×u_rel`、遠心 `−ρ Ω×(Ω×r)` を課す。**静止槽壁・バッフルには MRF を適用しない**。回転ゾーン境界の速度整合条件を定義。IBM 移動翼と同時起動不可。
- **FR-ROT-03**: 静止壁・バッフル＝内挿バウンスバック（Bouzidi/Ginzburg）、**移動翼＝IBM または moving-wall interpolated BB** と明確分離（#22）。STL 距離場の更新頻度・回転時幾何誤差を定義。
- **FR-ROT-04**: `Np = P/(ρ_l N³ D⁵)`, `P = Ω T_q`, `N = Ω/(2π)` を固定（2π 二重計上禁止）。ガス通気時は **非通気 `Np_0` と通気 `Np_g`、通気動力低下比**を別出力。`N_Q = Q/(ND³)`, `Q` の積分面・速度成分・時間/位相平均・逆流の扱いを定義（#23, #24）。
- **FR-ROT-05**（sliding-overset, 上位）: 重合格子 halo 補間を MPI と両立。基準級検証用。

### 4.4 高密度比二相流
- **FR-VOF-01**: 保存型 Allen-Cahn（§3 に固定）。質量保存誤差はベンチ別に規定 — **閉じた静止液滴 / 上昇単一気泡 / スパージャ開境界** で時間・格子解像度・流出入量込みの許容誤差を設定（#9）。Shan-Chen は本用途非採用。
- **FR-VOF-02** (rev.4 / REV-CFD-MJ-005 — dimensional fix): spurious currents on a
  static droplet are bounded by the (dimensionless) capillary number
  `Ca_spurious = μ_l |u|_spurious / σ < 10⁻³` (target We→0, resolution stated).
  The old `|u|·L/(σ/μ)` form carried a stray length dimension and is void. A
  length-bearing indicator, if wanted, is `Re_spurious = |u|_spurious L/ν_l` —
  a separate metric, never called Ca. well-balanced 化学ポテンシャル形式（#7 の係数関係を実装）。
- **FR-VOF-03**（スパージャ, rev.4 / REV-CFD-CR-001 — **phase-inversion fix**）:
  the sparger injects GAS; under the §3 definition (φ=1: liquid, φ=0: gas) the
  injected phase value is **φ=0**. The rev.1 text banned "plain `φ=1` + velocity
  Dirichlet" — that read as a liquid-injection ban and inverted the phase BC.
  Corrected requirements:
  - Choose from 気相体積流量境界 / 確率的気泡注入 / 解像オリフィス. A plain
    **`φ=0` + velocity Dirichlet alone is banned** — the injection model must
    simultaneously satisfy gas volumetric-flow conservation, pressure consistency,
    contact angle, and the `d_b/W`, `d_b/Δx` lower bounds（#10）.
  - **The scenario schema/API never exposes raw φ for inlets**: config says
    `inlet_phase: gas | liquid` and the core maps it (gas→φ=0, liquid→φ=1) —
    enforced by config validation (A-4 style). Outputs report `φ_liquid` and
    `α_g = 1−φ` with explicit names.
  - Acceptance: gas-inlet setting injects φ=0 (unit test); a sparger-only case
    balances injected gas volume vs. domain gas-volume increase within tolerance
    (VR-STR-02c precursor); no schema field accepts a raw φ boundary value.
  分裂・合体は「数値的に許容」までに弱め、実薄膜排液は解かない旨明記（#11）。
- **FR-VOF-04**（point-bubble）: 切替条件に `d_b/W, Eo, Re_b, α_g, We_b, 物質移動一貫性` を含める。hybrid 混在時の相間 質量・運動量・スカラー保存則を定義（#12）。
  **(rev.3, P1)** Population balance modelling (PBM) of the bubble-size distribution is
  required on the point-bubble path (breakup/coalescence kernels, e.g. Luo–Svendsen /
  Prince–Blanch): a mono-disperse point-bubble model cannot support the `d_32`
  acceptance of VR-STR-02 (internal consistency). Per-bubble gas-phase composition
  bookkeeping (component inventory and interfacial transfer budgets) must reconcile
  with FR-VOF-05. *Scope alignment (rev.2/§0)*: point-bubble is a relaxation extension
  (API-reserved in v1); this PBM requirement binds when that extension is implemented —
  in the resolved-phasefield default, `d_32` is measured from the resolved interface.
- **FR-VOF-05**: 界面物質移動を **解像界面（法線フラックス・分配係数・相別拡散）と point-bubble（`k_L a(C*−C)`）で分離**（#35）。Henry 則・Sherwood 数の適用範囲を明示。

### 4.5 分散粒子
- **FR-PART-01**: `α_p`/mass-loading で one/two/four-way 切替（閾値明示）。Schiller-Naumann の `Re_p` 適用範囲、反力散布カーネル、運動量保存検証を要求（#16）。中立浮力微粒子では Saffman/Basset/Faxen の要否を `d_p/Δx`・`St` で判定。
- **FR-PART-02**: 解像粒子法（PSM/Noble-Torczynski, Ladd/Aidun-Lu）へ切替可。
- **FR-PART-03**: 軌跡沿い `∫γ̇dt`・`max γ̇` を記録。**LES 追跡時は SGS 乱流分散（stochastic dispersion）を有効化**、または resolved-only を明記（曝露 PDF/CDF の格子依存を回避, #17）。
- **FR-PART-04 (rev.4 / REV-CFD-MJ-012 — four-way contact contract; Phase-2
  extension, API-reserved in v1)**: four-way coupling requires a soft-sphere
  normal-collision model with explicit parameters: restitution `e_n`, collision
  time `T_col`, spring `k`, dashpot `η`, max overlap `δ_max`, particle substep
  `Δt_p` (with `Δt_p ≲ T_col/10`).
- **FR-PART-05 (rev.4)**: when `d_p/Δx` does not resolve the lubrication gap, a
  lubrication correction (or calibrated implicit lubrication) is required; the
  applicability condition is stated with the model.
- **FR-PART-06 (rev.4 — config guard, initial delivery)**: while four-way is
  unimplemented/unvalidated, runs exceeding the `α_p` / mass-loading threshold of
  the two-way regime are **rejected at config validation** (A-4 style), with the
  threshold and its source stated in the error message. Initial delivery ships
  two-way + this guard; contact benches (particle–particle, particle–wall,
  settling, sheared suspension, overlap ≤ δ_max) gate the Phase-2 extension.

### 4.6 応力場評価（規約固定, #1, #2, #18, #19, #20）
- **FR-STRESS-01** (rev.4 / REV-CFD-CR-003, CR-004 — stage convention and forcing
  correction fixed by equations, not prose): strain rate is evaluated locally from
  non-equilibrium distributions. **The default stage is pre-collision /
  post-streaming** (the distribution as it arrives, before collide — the stage the
  standard coefficient below is derived for):
  ```
  f_i^{neq,pre} = f_i^{pre} − f_i^{eq}(ρ, u)        (u includes the F/2 correction)
  Π_neq_raw     = Σ_i c_iα c_iβ f_i^{neq,pre}
  Π_force       = −(Δt/2)(u_α F_β + u_β F_α)         (Guo forcing second moment,
                                                      for THIS engine's u/f_eq defs)
  Π_neq_corr    = Π_neq_raw − Π_force  =  Π_neq_raw + (Δt/2)(uF + Fu)
  S_αβ          = − Π_neq_corr / (2 ρ c_s² τ_eff Δt)
  ```
  `Π_neq_corr` is the ONLY normative definition; natural-language sign words are
  non-normative. The exact sign of Π_force is derivation-frozen against this
  engine's Guo discretisation **before implementation** and locked by a negative
  test (body-force Poiseuille must FAIL with the sign flipped — §8 VR-STR-03).
  **If the post-collision / pre-streaming stage is used instead** (e.g. inside a
  fused kernel where it is cheaper), the stage transform is mandatory:
  BGK: `Π_neq,pre = Π_neq,post / (1 − 1/τ_eff)`; MRT/cumulant: apply the inverse
  shear-moment relaxation `R(τ_shear)⁻¹` — then proceed with the equations above.
  The stress-evaluation API takes a required `neq_stage` enum
  (`PreCollision | PostCollision`) — no default-by-silence, misuse is a compile-
  or construct-time error (same philosophy as A-4/A-5 guards).
  cumulant/MRT ではせん断モーメント緩和率で係数補正。Smagorinsky 閉包の循環依存は
  **代数閉形式**で解く（`|Q|` から `τ_eff` を陽に求める; Hou et al. 型二次式）。
- **FR-STRESS-02**: 出力応力を **`resolved viscous` / `SGS` / `capillary` / `particle`** に分離定義。`γ̇=√(2S:S)`、第二不変量 `II_S`、von Mises は算出元テンソルを限定（#19）。
- **FR-STRESS-03**: 壁せん断はモード別に定義（**接線速度勾配再構成 / IBM forcing 積分 / MEM**）。内挿境界近傍の非平衡量が壁勾配を表さない場合の扱いを明記。検証は曲面移動壁を含む（#20）。
- **FR-STRESS-04**: 非 Newton `μ(γ̇)`（Carreau-Yasuda/Casson/power-law）と `μ_t` の合成則・反復手順・収束基準・`τ_min/τ_max`・LES 適用範囲を明示（二重計上・発散回避, #18）。

### 4.7 境界・重力・初期化（新規, #33, #34, #45, #47）
- **FR-BC-01**（上部境界, 必須指定）: `closed` / `free-surface` / `degassing-outlet` から選択。スパージング時は気相排出 outlet を要求（閉槽＋気相流入のみは気体蓄積で非物理, #33）。ヘッドスペース圧・液面変形・液面接触角を定義。
- **FR-BC-02**（重力）: 全相 `ρg`、動圧/静水圧分解、well-balanced hydrostatic test（静止成層で `|u|<ε`）を要求（#34）。
- **FR-BC-03**（濡れ性）: 壁ごとに接触角境界条件、滑り/非滑り、相場フラックス条件を定義（#47）。
- **FR-BC-04**（スカラー壁）: 無流束/吸着/反応壁を選択（#35）。
- **FR-INIT-01**: 初期 速度/圧力/相場/スカラー/粒子配置、インペラ ramp-up、ガス流量 ramp、統計サンプリング開始時刻、準定常判定基準を要求（#45）。

### 4.8 Extension & closure contracts (rev.3, P2)

- **FR-EXT-01**: Define explicit contracts for the trait/strategy extension points of
  §1 and for user-supplied closures — reaction rates `R_k`, non-Newtonian viscosity
  `μ(γ̇)`, body-force sources, and the relaxation-mode implementations
  (MRF / point-bubble / one-way / AMR):
  - input/output signatures with explicit physical vs. lattice units;
  - determinism (identical inputs → bit-identical outputs);
  - GPU evaluability (state-free, portable to wgpu);
  - error handling (NaN/divergence detection at the contract boundary);
  - schema versioning and backward compatibility.
  The primary boundary is Rust traits; foreign-language ABI/SDK is deferred to a
  separate API specification (see §10 product-layer note). The fidelity-default
  implementation is the default of each trait; relaxation implementations swap in
  under the same contract and are accepted via VR-STR-RELAX.
  *Implementation note*: this contract work is co-designed with the R-Phase 2 / B-1
  trait-boundary design (SOLVER_IMPROVEMENT_SPEC WP-B) — one design, two consumers.

---

## 5. 連成・時間積分（#28）

- **FR-COUP-01** (rev.4 / REV-CFD-MJ-007 — the dataflow is split by scalar mode so
  "active" is not silently one step lagged):
  **passive scalar**: 相場更新 → ρ/μ 場更新 → 力源合成(`F_s+ρg+F_g+F_p+F_rot`) →
  融合 collide-stream-moments → 境界 → スカラー ADE → 反応(split) → 粒子積分。
  **active scalar (fidelity default — predictor–corrector)**:
  scalar/reaction predictor → property update `ρ(C), μ(C), σ(C)[, T]` →
  力源合成(incl. `F_b^{scalar}`, Marangoni) → flow step → scalar ADE corrector →
  reaction corrector → property re-evaluation (→ optional flow–scalar iteration
  for stiff coupling). Time-lagged explicit feedback is allowed only as the
  flagged relaxation `active_scalar_lagged=true`, with stated applicability
  (weak feedback, non-stiff), stability conditions, and a lag-error benchmark —
  accepted via VR-STR-RELAX. Mode (coupled/lagged) is logged in run metadata.
  **強連成・剛直反応・表面張力波では演算子分割誤差・サブサイクリング・反復強連成を要求**。
  capillary time step `Δt_σ ≤ √(ρ̄ Δx³/(2πσ))`、粒子 `Δt_p`、反応 ODE `Δt_r` の各制約を課す。
  Acceptance: on the active-scalar standard bench (Marangoni or
  concentration-dependent viscosity), the feedback error converges under
  time-step halving (§8 VR-STR-06+/RELAX).
- **FR-COUP-02**: 反応ソルバは陽/陰/Rosenbrock-BDF を剛直性判定で切替。**負濃度制限・元素保存誤差・split 誤差の受入基準**を定義（#15）。
- **FR-COUP-03**: 無次元マッチングは §2.2 の優先順位＋feasibility check（#25）。
- **FR-COUP-04**: `probe_state_hash` ビット等価は**単一バックエンドの実装回帰限定**。物理妥当性・保存則は別基準（§8, #28, #42）。
- **FR-COUP-05**: AMR は上位オプション。有効化時は coarse-fine 保存補間・時間ステップ比・専用検証を要求（#29）。

---

## 6. 入出力・可視化

- **FR-IO-01**: 3D フィールド出力は **一様格子=VTI、構造曲線=VTS、非構造/AMR=VTU/AMR**（#43）。`φ` は拡散界面指標であってボイド率ではない（#36）。
  **ε_g processing definitions (rev.4 / REV-CFD-MN-014)** — every ε_g output
  carries filter width, averaging volume, and time window as metadata:
  - resolved-phasefield: `ε_g_raw = ⟨1−φ⟩_V` and
    `ε_g_thresholded(φ_c) = volume(φ<φ_c)/V`, default `φ_c = 0.5` — both output.
  - point-bubble: `ε_g_bubble = Σ_b V_b W_kernel(x−x_b) / V_filter`
    (kernel-smoothed void fraction).
  - hybrid: `ε_g_total = ε_g_resolved + ε_g_bubble` with double-count exclusion
    over the resolved region.
  Any ε_g indicator must be recomputable from a snapshot; experiment comparisons
  state which definition was used.
- **FR-IO-02**: 時間平均／位相平均統計（平均場・RMS・レイノルズ応力）。**位相平均は IBM/overset 非定常モードのみ**。MRF は回転座標平均/疑似定常として別出力（#37）。
- **FR-IO-03**: Web GUI に 3D 表示（スライス・等値面・せん断ヒートマップ・時系列プローブ）。既存 2D canvas を WebGL/WebGPU 拡張。
- **FR-IO-04**: 粒子累積せん断曝露のヒストグラム/CDF（SGS 分散の有無を明記）。
- **FR-IO-05 (rev.3, P3 — mixing metrics)**: Derived outputs for **blend time**
  (time until the coefficient of variation CoV of a tracer falls below a stated
  threshold) and **RTD** (tracer response `E(t)`, mean residence time, variance).
  The homogenisation threshold and the tracer injection/detection surfaces must be
  explicitly defined per scenario.
- **FR-IO-06 (rev.3, P3 — large-scale I/O & resilience)**: Full-field dumps are
  impractical at target scales (§7 budget); require **parallel I/O**
  (HDF5/ADIOS2-class) + compression + in-situ statistics / downsampling.
  **Deterministic checkpoint/restart with crash recovery** (bit-reproducible resume
  including RNG state, particle state, and statistics accumulators) is mandatory.
  Formats are sized against the §7 budget table.
  *Convergence note*: builds on SOLVER_IMPROVEMENT_SPEC B-5 (snapshot API),
  C-3 (per-rank parallel I/O), C-8 (distributed checkpoint) — reuse, don't duplicate.

---

## 7. 非機能要求

- **NFR-01（スケール・メモリ予算, #44）**: `O(10⁸–10⁹)` 格子。**格子あたりバイト数・分布数（D3Q27×相場×スカラー）・粒子数・GPU メモリ・I/O 量・チェックポイント頻度の予算表**を必須。1e9 格子×多分布は忠実度既定で **0.6 TB 級**、全 f64・複数スカラー・チェックポイント同時保持・I/O バッファ込みで **TB〜数 TB 級**（rev.2 修正）となるため wgpu マルチ GPU＋MPI 分割の見積りを添付。

  **予算表（rev.1b、忠実度既定構成・偏差格納・ping-pong ×2。1 セルあたり）**:

  | 構成要素 | 格子/型 | bytes/セル |
  |---|---|---|
  | 流体分布 f | D3Q27 × 2 × f32 | 216 |
  | 相場分布 g（保存型 Allen-Cahn） | D3Q19 × 2 × f32 | 152 |
  | スカラー分布 h（成分あたり） | D3Q7 × 2 × f32 | 56 |
  | moments・物性場（ρ, u×3, φ, μ_φ, ∇φ×3, ν_t, γ̇, τ_eff） | 12 × f32 | 48 |
  | マスク・フラグ | u8×2 | 2 |
  | 統計アキュムレータ（平均 u×3・RMS×3・レイノルズ応力 6 ほか） | ~13 × f32〜f64 | 52–104 |
  | 界面帯 f64 昇格（帯幅 ~2W、全セルの 5–10%、f+g の +368 B/帯セルを償却） | 償却 | +18–37 |
  | 界面帯の曲率・縮約作業域を含む場合の上乗せ（rev.2 検算） | 償却 | 〜+40 まで |
  | **合計（スカラー 1 成分）** | | **≈ 540–620 B/セル** |

  換算: **1e8 格子 ≈ 56–62 GB**（本機 M5 Max 128 GB で単ノード可、上限 ~1.5e8）／
  **1e9 格子 ≈ 0.56–0.62 TB**（f32 バルク）、全 f64 基準級で **≈ 1.1–1.2 TB**。
  粒子 10⁷ 個 × ~100 B = 1 GB（無視可能）。チェックポイントは分布の raw 保存で
  1 回あたりフィールド実体と同量（1e9 で ~0.5 TB/回）→ 頻度は I/O 帯域から逆算し
  ジョブあたり 2–5 回を既定とする。GPU: 8–16 GB/枚 → f32 構成で 1.3–2.6e7 セル/枚、
  1e9 格子は **40–80 枚のマルチ GPU or CPU クラスタ MPI が必須**（単 GPU では不成立）。
  結論: 忠実度既定での 1e9 はクラスタ専用。開発・検証は ≤256³（1.7e7 セル ≈ 10 GB）を
  標準とし、スケール実測は R3 クラスタ計画（CLUSTER_OPTIONS.md）に統合する。
- **NFR-02（精度ポリシー, #32, rev.2 語彙整理, rev.4 / REV-CFD-MJ-009 —
  enumerated so array design / GPU kernels / memory budget can bind to it）**:
  `precision_profile ∈ { full_f64, mixed_safe (default), mixed_fast }`:
  - **full_f64** (reference tier): all distributions, phase field, scalars,
    particle statistics, reductions in f64. High-density-ratio reference
    validations also run here.
  - **mixed_safe** (fidelity default = §1 profile): bulk distributions f32;
    **f64 fixed for**: `φ, ∇φ, κ(curvature), μ_φ, F_s, ρ(φ), μ(φ)`,
    distributions inside the interface band, all global reductions, torque,
    `Np`, `N_Q`, mass/volume counters, particle cumulative exposure.
    **interface_band = max(3W, 6Δx)** — provisional default; the band width is
    re-frozen by the W-VOF characterization (§10) and recorded in PHYSICS.md.
  - **mixed_fast** (relaxation extension): single-phase / weak-coupling only;
    permitted only when density ratio ≤ stated limit AND the Ca_spurious and
    mass-drift validations pass; config validation rejects out-of-range use.
    Accepted via VR-STR-RELAX-f32.
  Each profile has an array-type table and memory-budget column (§7).
  `ρ_l/ρ_g≈10³`・`Ca_spurious<10⁻³`・質量保存要求と整合させる。
- **NFR-03（性能）**: 融合 `step_band` に相場・スカラー・forcing 統合、リング二重化・SoA plane-major の 3D 拡張維持。
- **NFR-04（決定性）**: 縮約は決定的順序。GPU/MPI は許容誤差ベース回帰（ビット等価は単一バックエンド限定, #42）。

---

## 8. 検証・受入基準（定量化, VALIDATION.md **T17** として配線済み。閾値付き, #38–#42）

検証テストは codex/Opus が本仕様から敵対的に作成し実装と分離。

**Band governance (rev.4 / REV-CFD-MJ-010 — reconciling "numbers now" with the
experiment-driven freeze protocol)**: every VR-STR item carries a **provisional
numeric band from day one** (table below — these are the MVP gate). Bands are
finalized by the established protocol (implement → characterize → record rationale
in PHYSICS.md → freeze in VALIDATION.md T17) under one asymmetric rule:
**tightening a band is always allowed; loosening a provisional band requires a
recorded physical rationale in PHYSICS.md** (reference uncertainty, method order,
resolution limit — as exercised for T15.5). This removes both failure modes:
un-testable placeholder specs AND post-hoc self-serving thresholds.
Each test is specified with: metric / target·reference / tolerance / resolution /
time window / backend / pass-fail rule (the T17 row format).

**Provisional bands (MVP gate; supersede the "±許容%" placeholders)**:
- Rushton `Np` vs experimental correlation: **±10%**
- PIV/LDA velocity profiles (VR-STR-01): **L2_rel < 15%, L∞_rel < 30%** per line
- static droplet mass drift: **< 0.1% / 1000 steps**; advected droplet (one period,
  periodic box): total mass drift **< 0.1%** (CR-002 acceptance)
- single-bubble terminal velocity vs Grace (02a): **±10%**
- `k_L a` vs correlation (02c): **±25%**
- well-balanced static stratification (VR-STR-06): **max|u| < 10⁻⁶ (lattice units)**
  at ρ ratio 10³ — provisional; retighten after discretisation freeze
- GPU/MPI cross-backend drift (VR-STR-05): mean quantities **< 2%**, higher-order
  statistics **< 5%** (bit-equality stays single-backend-only)
- `Ca_spurious < 10⁻³` (already fixed, dimensionally corrected FR-VOF-02)

**Mandatory negative/consistency tests added by rev.4**:
- forcing-moment sign negative test: body-force Poiseuille FAILS with Π_force sign
  flipped (CR-004);
- stress stage-convention cross-check: pre-collision evaluation vs post-collision+
  transform agree within tolerance on Couette/Poiseuille/Taylor–Couette (CR-003);
- J_ρ consistency code-path check + droplet advection conservation (CR-002);
- sparger phase unit test: gas inlet injects φ=0, gas-volume balance closes
  (CR-001);
- scalar total-mass conservation in the phase-wise form (MJ-011);
- active-scalar dt-halving convergence (MJ-007).

- **VR-STR-01（単相撹拌）**: 標準 baffled tank（`D/T`, `C/T`, blade geometry, バッフル数を固定）、指定 Re 範囲・非通気。Rushton `Np`＝実験相関±許容%、翼吐出速度プロファイルを PIV/LDA 基準測線で `L2/L∞rel` 閾値照合（#38）。
  **(rev.3, P4) Reference datasets**: Wu & Patterson (1989) LDA; Deen et al. (2002)
  PIV (standard Rushton, D/T=1/3, 4 baffles); standard `Np` correlations. Numeric
  bands are frozen via the T17 experiment-driven protocol — not hardcoded here.
- **VR-STR-02（気液, rev.2 で 02a/b/c に分割）**: **02a 単一気泡** = `U_t` を Grace 線図 Eo-Mo-Re と相対誤差照合。**02b 気泡群** = `ε_g` 空間分布・群上昇速度（hindered rise）・合体/分裂を許す場合の `d_32`・BIT 使用時の乱流強度（`ν_t` 応答）。**02c 撹拌槽通気** = `ε_g, d_32, k_L a` の実験相関比（#39）。
  **(rev.3, P4) References**: single bubble = Grace diagram (Eo-Mo-Re); aerated tank =
  published `ε_g`/`d_32`/`k_L a` data and correlations. In point-bubble / RELAX-PB
  evaluations, `d_32` presupposes the FR-VOF-04 population balance (P1); in the
  resolved-phasefield default it is measured by interface segmentation.
- **VR-STR-03（せん断・応力）**: 製造解（MMS）単相、曲面 Couette、回転円柱、非 Newton Poiseuille、多相静止液滴を分け、**格子収束次数**と `L2/L∞` を設定。壁近傍 `L∞` の発散的厳しさを考慮した測線設計（#40）。
- **VR-STR-04（スカラー/反応）**: Taylor-Aris 分散、既知 `Da` 反応拡散前線、`k_L a`（算出式＝界面積分か相関か明示）。各々の許容誤差・対象 `Pe/Da/Sc`・境界条件を指定（#41）。
- **VR-STR-05（連成回帰・保存）**: `probe_state_hash` は単一バックエンド回帰限定。**質量・運動量・スカラー総量・気相体積・粒子数・エネルギー様量のドリフト閾値を個別設定**。エネルギー様量（運動エネルギー・界面自由エネルギー・粒子運動エネルギー）は**厳密保存でなく非物理ドリフトの監視量**として扱う（rev.2）。GPU/MPI は許容誤差ベース（#42）。
- **VR-STR-06（well-balanced）**: 静止成層で `|u|<ε`（#34）。**06+（rev.2）**: active スカラー ON かつ `C≡C_0` で同一の静止性を満たす（`F_b^{scalar}` の厳密ゼロ退化）。σ 可変形の `∇σ=0` 退化（σ 一定基準形との一致）も本群に置く。
- **VR-STR-07（初期化非依存性）**: 助走・統計開始条件を変えて準定常統計が閾値内で一致（#45）。
- **VR-STR-RELAX（緩和モード同等性, rev.2 新設 — 所見1）**: 緩和拡張は**対応する忠実度基準解**に対する相対劣化で受入する。各軸の比較対象・測定量（許容差は characterize→freeze）:
  - **RELAX-MRF**: 同一幾何・同一 Re の IBM-inertial（または sliding-overset）基準に対し `Np`・翼吐出速度測線・平均速度場・トルクの許容差。適用は定常近似が成り立つ構成に限定。
  - **RELAX-PB（point-bubble）**: resolved-phasefield 基準に対し `ε_g`・`d_32`・`k_L a`・運動量/スカラー収支の許容差。適用範囲を `d_b/Δx, d_b/W, α_g` で限定（FR-VOF-04 の切替条件と同一）。
  - **RELAX-1W（one-way）**: two-way 基準に対し粒子統計の許容差と、反作用無視が許される mass-loading 上限。
  - **RELAX-AMR**: uniform 基準に対し保存量・界面位置・トルク・速度場ノルム・coarse-fine 境界通過時の収支誤差。
  - **RELAX-f32（積極的 f32）**: 忠実度プロファイル（または全 f64）基準に対し保存量ドリフト・`Ca_spurious`・`Np`・界面曲率・縮約量の許容差。適用は単相/弱連成に限定（NFR-02）。

---

## 9. 主要技術リスク

| # | 難所 | リスク | 緩和策 |
|---|---|---|---|
| 1 | 高密度比二相 (`10³`) | 界面不安定・寄生流・f32 丸め | well-balanced 相場＋D3Q27＋f64 縮約、point-bubble 代替（§1） |
| 2 | 高 Re 安定性 | 発散・超粘性・positivity | cumulant＋WALE、代数閉包 `τ_eff`、limiter（§4.1,4.6） |
| 3 | 回転境界保存性 | IBM すべり・トルク誤差 | multi-direct-forcing、overset 基準検証、閾値化（§4.3） |
| 4 | 連成剛直性 | 反応・界面・回転の時間スケール乖離 | 演算子分割誤差評価・サブサイクリング・capillary dt（§5） |
| 5 | 計算コスト・メモリ | 1e9 格子×多分布=TB 級 | メモリ予算表・マルチ GPU＋MPI、AMR 上位オプション（§7） |
| 6 | 物質移動・BIT | 未解像通気の乱流・移動不足 | Sc_t SGS、BIT 生成項、解像/point 分離（§4.2,4.4） |

---

## 10. 設計判断（確定事項と残実装詳細）

**確定（rev.1a, PM 決定）**: 全軸の既定を**忠実度最優先**に振り、低コスト近似は後付け拡張点として実装する（§1 設計原則）。従来の未解決4件は次のとおり解消:

- 界面既定 = `resolved-phasefield`（界面・物質移動忠実度優先）。`point-bubble` は緩和拡張。
- スカラー既定 = `active`（物性帰還あり）。`passive` は緩和拡張。
- 格子既定 = `uniform`（完全解像）。`block-AMR` は緩和拡張。
- 精度既定 = 忠実度プロファイル（界面近傍・保存量・縮約 f64、遠方バルク f32）。全 f64 は基準級、積極的 f32 は緩和拡張。

**残実装詳細（決定ではなく仕様詰め）** — status as of rev.3:

- `active` スカラーの帰還対象（σ・粘性・密度・[温度]）の具体式と安定化（Marangoni 含む）。
  → researched: docs/proposals/active-scalar-feedback.md. **One derivation is mandatory
  before implementation** (Marangoni coefficient consistency with the (κ,β) convention,
  §3). Thermal axis recommended as API-reserved extension.
- 緩和拡張ごとの忠実度基準解に対する許容誤差閾値（§8 VR へ追記）。
  → structure defined as VR-STR-RELAX (rev.2); numeric bands frozen at relaxation
  implementation time.
- 忠実度プロファイルの f64/f32 境界（界面近傍の帯幅・縮約範囲）。
  → frozen experimentally during W-VOF implementation (characterize→freeze).
- 各モード軸の trait 境界（strategy 差し替え点）の API 定義。
  → contract requirements fixed as FR-EXT-01 (§4.8); concrete Rust API co-designed
  with R-Phase 2 / B-1.

**Product-layer scope note (rev.3, P5)**: GUI/CAD & STL import, materials DB,
Python/CLI SDK, parameter sweep & optimizer, cloud/cluster/queue integration,
packaged validation assets, and competitive benchmark tables are **out of scope for
this solver specification**. They are version-managed in separate volumes
(Product Requirements / API Specification / Validation Pack / Performance Benchmark).

---

## 11. Implementation dependency graph (rev.3 — priority + dependency DAG, not stage gates)

Items with no dependency edge between them are implemented **concurrently**
(parallel-agent worktrees, per the standing parallelization directive).
Mapping to the PLAN.md M-F delegation tracks is noted per row
(MF-α…ζ are the delegation bundles; W-items are the fine-grained DAG nodes).

| Item | Hard deps (must precede) | Parallel | Notes / PLAN track |
|---|---|---|---|
| W0 core basis (D3Q19/27, cumulant, Guo forcing) | — (strengthens M-C 3D basis) | — | prerequisite for all; = MF-α |
| W-EXT trait contracts (FR-EXT-01) | W0 | yes | early definition = prerequisite of all relaxation modes; low cost, high leverage; co-designed with R-Phase 2 B-1 |
| W-UNIT unit/nondimensional feasibility (§2.2) | W0 | yes | independent, early |
| W-STRESS stress fields (FR-STRESS) | W0 | yes | top priority (primary output + prerequisite of LES & particle exposure); ⊂ MF-β |
| W-ROT rotating IBM (FR-ROT-01) | W0 | yes | prerequisite of Np/N_Q; MRF/overset live behind W-EXT as relaxation/reference tiers; = MF-δ |
| W-GRAV well-balanced gravity (FR-BC-02) | W0 | yes | prerequisite of the interface track; ⊂ MF-γ |
| W-SCAL passive scalar ADE (§3 scalar eq.; SGS flux part waits on W-LES) | W0 | yes | ⊂ MF-ε |
| W-LES turbulence SGS (FR-LES) | W-STRESS | conditional | \|S\| closure needs the stress evaluation; ⊂ MF-β |
| W-VOF resolved interface (FR-VOF-01/02) | W-GRAV | conditional | fidelity default; hardest item; **critical path**; ⊂ MF-γ |
| W-PART particles + cumulative exposure (FR-PART) | W-STRESS (SGS dispersion: W-LES) | conditional | exposure integral needs the γ̇ field; ⊂ MF-ε |
| W-REACT reaction / active feedback (§3, FR-COUP-02; active feedback needs W-VOF) | W-SCAL | conditional | ⊂ MF-ε |
| W-BUB point bubbles + PBM + interfacial transfer (FR-VOF-03/04/05) | W0, W-SCAL, W-EXT | conditional | relaxation extension (API-reserved in v1, per §0) |
| W-BCTOP top boundary / degassing / contact angle (FR-BC-01/03) | W-VOF | conditional | ⊂ MF-γ |
| W-COUP coupling loop (FR-COUP) | active subsystem set | incremental | grows as tracks land; ⊂ MF-ζ |
| W-IO I/O & analysis (FR-IO incl. -05/-06) | each producing subsystem | incremental | Np←ROT, blend/RTD←SCAL, exposure←PART; ⊂ MF-ζ |
| W-VAL validation T17 (VR-STR-01–07, RELAX) | each subsystem | yes | codex adversarial authorship, separated from implementation |

**Parallel waves** (sets that start together):
1. After W0, mutually independent: **W-EXT / W-UNIT / W-STRESS / W-ROT / W-GRAV / W-SCAL** (6-way parallel).
2. After their deps: W-LES (←STRESS) / W-VOF (←GRAV) / W-PART (←STRESS) / W-REACT (←SCAL).
3. Later: W-BCTOP (←VOF) / W-BUB (←SCAL,EXT) / active feedback & interfacial transfer (←VOF).
4. Cross-cutting throughout: W-COUP / W-IO / W-VAL.

**Critical paths** (staff first):
`W0 → W-GRAV → W-VOF → W-BCTOP/interfacial transfer` (interface chain — longest, hardest) and
`W0 → W-STRESS → W-LES → W-PART` (stress/exposure chain).

*Boundary decisions upheld (rev.3)*: throughput/scaling KPIs stay delegated to
CLUSTER_OPTIONS.md (R3) — not duplicated here; no hardcoded numeric thresholds
(P4 adds dataset names only — bands freeze via the T17 protocol); the product
ecosystem (P5 list) lives in separate volumes.
