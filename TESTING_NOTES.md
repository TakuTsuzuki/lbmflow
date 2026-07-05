# TESTING_NOTES

## MPI ops bundle (2026-07-05)

1. Persistent MPI exchange buffers were verified with
   `PATH=$HOME/.local/openmpi/bin:$PATH cargo test -p lbm-core --release --features mpi dist::tests -- --nocapture`.
   The tested steady exchange path reuses `MpiExchange`'s per-axis typed send/receive buffers, and the hot
   population/scalar pack/unpack helpers now iterate layer indices without allocating a temporary index vector.
   The one-rank cargo smoke completed successfully; Open MPI printed a local TCP bind warning in this sandbox,
   but the test process exited green.

テスト作者（codex）とエンジン作者（PM/Fable）の連絡帳。
新しい不一致は末尾に追記する。処理済み項目は Disposition を残して保持。

## 処理済み（2026-07-05 PM triage、詳細: docs/PHYSICS.md）

1. `t6_f32_...`: f32 一様場の運動量成長誤差 5.3e-4 (>1e-4)
   → **仕様変更**: コヒーレント丸めバイアスは f32 の本質特性。許容 5e-3 に改定
   （VALIDATION.md T6 更新済み）。診断集計は f64 化した（エンジン変更）。
   **テスト側の対応が必要**: 閾値を新仕様に合わせる。

2. `t4_...flow_rate_constancy`: 流量一定性 4.8e-3 (>1e-6)
   → **仕様/API バグ**: 放物線流入 API が無かった → `set_inlet_profile` 追加。
   流量は質量流束 Q=Σρux で、流出境界の直前は Zou-He 固有のスタッガード層
   （O(Ma²), 減衰長~4セル）があるため**バルク領域（流出から24列以上）で ≤1e-4**
   に改定（VALIDATION.md T4 更新済み）。
   **テスト側の対応が必要**: プロファイル API 使用 + バルク判定に書き換え。

3. `t5_pressure_sign_reversal`: 反対称 4.5e-5 (>1e-12)
   → **仕様バグ**: 慣性項は2次なので厳密反対称は成立しない。
   厳密角度は「Δρ反転 + x鏡映 = 厳密一致 ≤1e-12」に置換、単純反転は ≤5e-3 相対
   （VALIDATION.md T5 更新済み）。
   **テスト側の対応が必要**: 2 つの角度に分割。

4. `t10_tau_051_cavity`: NaN
   → **仕様確定**: U=0.05（Re≈1890）なら安定、U=0.1 は発散（実測）。
   T10 パラメータを τ=0.51, N=128, U=0.05, Λ=3/16 に確定（VALIDATION.md 更新済み）。
   **テスト側の対応が必要**: パラメータ変更。

## 処理済み（2026-07-05 PM triage #2、詳細: docs/PHYSICS.md 2026-07-05 各節）

order #2 の 5 件の dispositions:

1. `t7_re400`: **参照データの既知の誤植**。Ghia Re=400 の v(0.9063)=−0.23827 は
   流通データ自体が誤り（隣接点と不連続・出典 gist にも注記あり・我々の解は
   −0.37657 で滑らか）。→ この 1 点を RMS から除外（VALIDATION T7 更新済み）。
   **テスト側対応**: 除外処理 + 出典コメント。
2. `t7_orientation`: **両側バグ**。(a) エンジン: リムコーナーの wall_u が適用順
   依存 → 「速い壁が勝つ」規則に修正済み。(b) テスト: Left/Right の対称写像が
   誤り（[0,−U] 左蓋は回転でなく反対角鏡映）。正しい写像は VALIDATION T7 /
   PHYSICS.md に記載。エンジンは正しい写像で L∞ ~4e-16 を実証済み
   （examples/probe_equivariance.rs）。**テスト側対応**: 写像修正（Bottom は正しい）。
3. `t8_re20`: **仕様の幾何不整合**。周期境界+ブロッケージの Cd=2.55 は物理的に
   正しい値。→ T8 を Schäfer-Turek 2D-1/2D-2 に全面再定義（VALIDATION T8 更新済み）。
   **テスト側対応**: validation_cylinder.rs を新仕様で書き直し。
4. `t8_re100`: 同上（2D-2 へ）。
5. `t9_outflow`: **仕様改定**。ゼロ勾配流出の圧力反射は固有特性（実測 ratio 11.3、
   衝突演算子非依存の見込み）。→ 上限 15 に改定（VALIDATION T9 更新済み）。
   convective outlet は Phase 7 バックログ。**テスト側対応**: 閾値変更。

## 新規不一致（2026-07-05 codex adversarial test order #2）

1. `t7_lid_driven_cavity_re400_matches_ghia`: N=129, U=0.1, TRT Λ=3/16, Re=400。
   `run_to_steady(1000, 1e-8, 200000)` は 99000 step で定常判定に到達したが、
   Ghia et al. 1982 中心線 u/v の RMS 誤差が 2.6577415383317194e-3 で、
   仕様上限 2.0e-3 (= 0.02U) を超過。Re=100 は同じテストで合格、Re=1000 ignored
   は 2:10.57 で合格。

2. `t7_re100_cavity_is_exact_under_four_lid_orientations`: N=129, U=0.1, TRT Λ=3/16,
   Re=100。蓋向きを左へ回したケースを 2000 step 後に回転写像で比較すると
   L_inf = 2.843743051315205e-2 で、仕様上限 1e-10 を超過。テスト側では座標写像と
   左右壁の接線速度符号を修正済み。左壁 MovingWall 経路または壁更新の回転対称性を
   要確認。

3. `t8_re20_cylinder_steady_drag_is_in_reference_band`: D=20, domain 440x160,
   left VelocityInlet U=0.05, right PressureOutlet rho=1, top/bottom Periodic,
   cylinder center (110,80), Re=20, TRT Λ=3/16。7000..10000 step の平均で
   Cd = 2.5454767275786616、仕様帯 [1.8, 2.4] を超過。

4. `t8_re100_cylinder_vortex_shedding_has_expected_st_cd_cl`（ignored）:
   D=20, domain 440x160, right Outflow, off-centre cylinder y=81, Re=100。
   80000 step 実行で mean Cd = 1.6199794592087982 となり、仕様帯 [1.2, 1.5] を超過。
   St 判定は先に通過。

5. `t9_outflow_cylinder_wake_long_run_stays_sane`（ignored）:
   D=20, domain 440x160, right Outflow, off-centre cylinder y=81, Re=100。
   100000 step 実行で NaN/Inf と backflow 判定は先に通過したが、
   near-outlet pressure RMS ratio = 11.32538182631078
   （near = 2.0001796481913235e-3, mid = 1.766103499967275e-4）で仕様上限 3 を超過。

## 新規メモ（2026-07-05 codex adversarial test order #4）

1. T11b の文面「bottom BounceBack, others Periodic」は現 API では構築不可。
   `SimConfig::validate` が周期境界の軸ペアを必須にしているため、bottom BounceBack
   + top Periodic は `ConfigError::UnpairedPeriodic { axis: "y" }` になる。
   追加テストは left/right Periodic + bottom/top BounceBack（上壁は液滴から遠い）で
   G_w 特性を凍結した。

## 新規メモ（2026-07-05 codex adversarial test order #5）

1. `cargo test --release -p lbm-core -- --include-ignored` はテスト本体を最後まで通過したが、
   doctest 段で `crates/lbm-core/src/multiphase.rs` の `MultiComponent` ignored doc snippet
   （line 77）がコンパイル対象になり失敗した。スニペットは `MultiComponent` import と
   `a`/`b` simulation 定義を省略した疑似コードで、通常の default run では ignored doctest
   としてスキップされる。今回の作業範囲は `crates/lbm-core/tests/**` と
   `TESTING_NOTES.md` のため、`src/**` doctest は修正せず証拠だけ記録する。

## 解決状況（2026-07-05 codex adversarial test order #3）

1. `t7_re400`: fixed-by-spec / fixed-in-test — 既知の誤植 datum を RMS から除外。
2. `t7_orientation`: fixed-by-engine / fixed-in-test — リムコーナー修正済み、Left/Right 写像を仕様通りに修正。
3. `t8_re20`: fixed-by-spec / fixed-in-test — Schäfer-Turek 2D-1 に全面更新。
4. `t8_re100`: fixed-by-spec / fixed-in-test — Schäfer-Turek 2D-2 に全面更新。
5. `t9_outflow`: fixed-by-spec / fixed-in-test — 圧力 RMS ratio 閾値を T9 の 15 に更新。

## 新規メモ（2026-07-05 codex adversarial test order #6）

1. Core V2 `Solver` には分割対応の single-component Shan-Chen driver が直接は露出していない。
   Shan-Chen は現状 V1 互換 `Simulation` facade 経由で利用できるが、T13 の 2x2 seam 上で
   `Solver<D2Q9, ..., InProcess>` に対して密度から force field を再計算する公開 API は未整備。
   そのため `t13_adversarial.rs` では gap を残しつつ、同じ下層経路である per-cell
   `force_field` を各 subdomain の compact core に直接設定し、四分割 corner 上の droplet 型
   force field が一枚岩と一致することを検査する。

## 新規不一致（2026-07-05 codex adversarial test order #6）

1. `d3q19_lattice_properties_from_all_angles`: D3Q19 の face closure constant を
   `assert_eq!(closure, 1.0)` で検査すると、`XNeg` で `closure = 1.0000000000000002`
   となり失敗する。既存 unit test は `abs <= 1e-15` で許容しているが、今回の発注条件
   「closure constant exactly 1」に対しては red。原因は `1/3, 1/18, 1/36` の f64 加算順に
   よる丸めの可能性が高いが、テーブル/API が「exact」を名乗るなら rational/整数式での
   定数化、または仕様文言の明確化が必要。

## 処理済み（2026-07-05 PM triage #3: codex order #6）

- `d3q19_lattice_properties_from_all_angles` の閉包定数「厳密 == 1.0」失敗
  → **テストの過剰厳密**（カテゴリ: テストのバグ）。エンジンは解析値 T::one() を
  ハードコード済み（kernels.rs zou_he）で物理は正しい。テスト自身の f64 総和が
  加算順で 1 ulp ずれるだけ（XNeg: 1.0000000000000002）。判定を 4 ulp 許容に修正。
  分割不変性への攻撃 8 種（分割線上円柱+プローブ/L字3分割跨ぎ/蓋・流入の分割跨ぎ/
  不均等[3,1,1]/最小幅ガード/4分割コーナー液滴/20k長時間）は**全て耐えた**。
  Shan-Chen V2 ネイティブAPI未整備のギャップは記録どおり（M-C/M-D で配線）。
## 新規メモ（2026-07-05 M-B Wgpu backend / T14 実装）

1. **T14 圧力 BC の許容線**: Zou–He 圧力面は `un = 1 - closure/rho_bc` が
   O(1) スケールの closure の丸め差（Metal fast-math の逆数除算・再結合、
   ~ulp(1) ≈ 1.2e-7）を**そのまま面の法線速度に**写像する（速度 BC では同じ
   除算誤差が rho に落ち、f への寄与は u_n 倍で減衰する — 非対称）。
   実測: CPU↔GPU 差は圧力面に固定（argmax が面に張り付き）、t=2 で ~2.2e-7、
   t=100 で ~2.5e-6（u0=0.1 の速度相対 2.5e-5）。**CPU-vs-CPU で rho_bc を
   1 ulp だけ摂動した対照実験が同じ成長曲線を再現**（t=100 で ~1.5e-6）した
   ため、バックエンド欠陥ではなく BC の条件数と確定。
   → **Disposition**: T14 は 6 構成（TGV/キャビティ/プロファイル流入チャネル/
   円柱+プローブ/セル別力/Convective）を厳格線 1e-5 で凍結し受入を満たす。
   圧力チャネルは第 7 構成として文書化済み緩和線 1e-4 + 恒久対照テスト
   `t14_pressure_bc_ulp_sensitivity_control`（1 ulp 摂動ドリフトが 1e-6..1e-5
   の帯にあることを常時検証; 帯を外れたら許容線を見直す）で凍結。

2. **GPU ベンチの測定衛生**: ユニファイドメモリでは並走する CPU スイート
   （本日: 3D エージェントの t15、load ~38）が GPU の DRAM 帯域を食い、
   帯域律速カーネルは 1024²/2048² で 15-25% 落ちる（SLC に乗る 512² は鈍感）。
   proto 凍結値との比較は**同一時間窓で proto を併走**させること
   （examples/bench_gpu.rs のヘッダに手順と 2026-07-05 同窓実測を記録:
   -6.4% / -10.7% / -13.6%、合格線 ±20% 内）。
## 新規（2026-07-05 M-C 3D 実装）

1. **T15.3 の参照値注記の誤記**: VALIDATION.md T15.3 は受入基準を
   「Schiller-Naumann 相関 Cd = (24/Re)(1 + 0.15 Re^0.687) の ±10%」と定義し、
   括弧書きで「Re=20: ≈2.09、Re=100: ≈1.09」を添えていたが、式の値は
   **Re=20 → 2.6095**（2.09 は Re≈28 の値。Re=100 → 1.0917 は正しい）。
   テスト（crates/lbm-core2/tests/t15_3d.rs）は一次基準である**式**に対して
   ±10% を判定する。VALIDATION.md の括弧書きは 2.61 に訂正済み（許容幅は不変更）。

2. **T15.3 球抗力 Re=20 (D=24) が公称 D 基準で仕様帯 ±10% を超過**（PM triage 依頼）。
   測定（momentum-exchange、窓平均 Cd = 500 step 窓の平均が相対 5e-4 で収束するまで。
   窓平均にしたのは、速度流入↔圧力流出間の弱減衰音響定在波（減衰 ~1/(νk²) ≈ 6e4 step）
   の O(Ma) リップルが瞬時サンプルの収束判定を ~1e5 step 停滞させるため）:
   | 構成 | Cd 実測 | SN(Re) | 公称D 誤差 |
   |---|---|---|---|
   | D=24, Re=100, 192×128×128, u=0.10 | 1.1698 (11k step) | 1.0917 | **+7.2% 合格** |
   | D=24, Re=20, 192×128×128, u=0.05 | 2.9551 (39k step) | 2.6095 | **+13.2% 不合格** |
   | D=12, Re=20, 96×64×64, u=0.06（軽量版・帯±25%） | 2.9790 (3k step) | 2.6095 | +14.2% 合格 |
   **原因分析**: half-way BB の staircase 球は流体力学的半径が公称より約半リンク大きい
   （Ladd 較正の古典的事実）。r_h = r + 0.5 で Cd と Re を再正規化すると誤差は
   **+0.6% / +7.1% / +2.3%** に潰れ、物理（エンジン）は正しい。残差 +7.1%（Re=20, D=24）
   は低 Re で遮蔽されない周期側面イメージ（D/L_y = 0.19、Stokes 的 O(D/L) 補正）が主。
   つまり「公称 D・±10%・D=24・ブロッケージ≤3%」の組は Re=20 では物理的に両立しない
   （半リンクバイアス ~ +2/D = +8.5% だけで帯をほぼ使い切る）。
   **triage 候補**: (a) D を流体力学的直径（D_h = D+1）で定義（Ladd 流儀。全ケース余裕で
   合格、帯の締め付けも可能）、(b) 公称 D のまま D ≥ 48 に引き上げ + 側面 ≥ 8D、
   (c) Re=20 のみ帯を +15% に拡大。テストは仕様どおり公称 D ±10% のまま**弱めずに**
   コミット（#[ignore] 重量級のみ red、デフォルトスイートは緑）。

## 処理済み（2026-07-05 PM triage #4: M-C 球抗力）

- T15.3 球 Re=20/D=24 の +13.2% 帯超過 → **仕様の正規化定義バグ**。
  half-way BB の staircase 球は流体力学的半径 r_h = r+0.5（Ladd 較正）を持ち、
  (Cd_h, Re_h) ペアで再正規化すると 3 ケースが +0.6%/+7.1%/+2.3% に収束
  （エンジン正常）。VALIDATION T15.3 を D_h 定義に改定、テストは sn_hydro 化。
  併せて仕様の誤記（Re=20 の SN 値 2.09 → 正しくは 2.6095）も訂正済み。

## 新規（2026-07-05 M-E CpuSimd 融合バックエンド）

1. **等価性ゲート実測（tests/backend_simd_equiv.rs）**: CpuScalar vs CpuSimd、
   8 シナリオ（2D TGV/キャビティ/プロファイル流入チャネル/円柱+プローブ/
   セル別力 BGK/Convective、3D TGV/ダクト）× f64/f32、150〜400 step。
   場（rho/u/流体セル f 平面）の実測 worst |Δ| は f64 で ~6e-14
   （ゲート 1e-11）、f32 は全構成ゲート 1e-6 内。probed_force はリンク寄与を
   CpuScalar の (x,q) セル順に並べ替えて再生する設計でビット等価。
   InProcess 2×2（円柱+プローブが縫い目跨ぎ/周期 TGV/3D 2×2×1 ダクト）
   vs 単一領域 CpuScalar も ≤6.4e-12（f64 部分和再結合のみ）。
2. **f32 の外延診断は次元整合ゲートが必要**: total momentum は N セルの
   f64 総和なので、バックエンド最終 ulp ドリフト（~1e-9/セル）が N 倍に
   蓄積する（実測 96×64 で 1.1e-6、48×20×20 で 3.8e-6 — N 線形）。
   v1_match の f32 ケースが場のみ比較していたのはこのため。ゲートは
   場 1e-6（絶対）/ 質量・運動量 1e-6·N_fluid / probe 1e-5（実測 ≤1.3e-6）
   に整理し、テスト頭書に測定根拠を記載。
3. **カーネル形状の測定断面**（詳細 docs/PERFORMANCE.md「V2 CpuSimd」節 +
   コード内 doc）: kernels.rs の逐語 DAG は V1 ペア形式比 1T −16%;
   D3Q19 flat 展開はスカラー化（vec/scalar 命令 18/285）→ blocked 化で
   3.0x; blocked への src/dst 別ビューはエイリアス検査でベクトル化崩壊
   （−30%）; y ストリップリングは本機では SLC が平面リングを吸収するため
   −20%; バンド 2 倍過剰分割 −8%。全て実装→実測→棄却の記録付き。
4. **3D 12T の 2 倍未達を記録**: 128³ 12T で f32 1.9x / f64 1.4x
   （1T は 3.0x/2.0x で達成）。支配要因はバンド端二重衝突（+19%）と
   P/E 異種コアでの粗粒度バンド不均衡（scalar は行粒度 stealing で
   スケール 7.5x、融合 5.2x）。改善候補: バンド端衝突の共有（要同期）、
   または nz を跨ぐ動的バンドサイズ。

## GPU ops bundle (2026-07-05 cx-gpu-ops)

1. **GPU adapter availability in this sandbox**: both the pre-change baseline
   attempt and the post-change benchmark attempt failed before measurement
   because wgpu could not acquire an adapter:
   `bench_gpu requires a usable GPU adapter: no usable GPU adapter was found`.
   Command used:
   `cargo run -p lbm-core --release --features gpu --example bench_gpu -- --gpu-only`.
   Therefore the requested 1024² MLUPS regression, 2048² sync+diagnostics
   speedup, and no-force/no-solid live allocation measurement are not
   available from this sandbox run.

2. **Implemented non-GPU-verifiable checks**:
   `cargo test -p lbm-core --release --features gpu gpu:: --lib` passed
   10 GPU module unit tests, covering submit-chunk calibration, poll-error
   conversion, resource-limit rejection, naga parse+validate of generated WGSL,
   and `BcParams` field-order consistency.

3. **Release gates**:
   `cargo test --workspace --release` passed.
   `cargo test --workspace --release --features gpu` passed, including all
   8 T14 GPU backend-equivalence tests.

## 新規（2026-07-05 V1 引退作業）

1. **sync-tests.sh の置換が macOS では無効だった**: `sed -E 's/\blbm_core\b/…/'` の
   `\b` は BSD sed 非対応で無置換コピーになっており、「compat へ再標的化済み」の
   複製テスト 16 ファイルは実際には dev-dependency の V1 を直接テストしていた
   （M-A の「56+ テストが compat 経由で緑」は未検証状態だった）。perl 置換に修正して
   再同期し、compat 実経由で全複製スイート緑を実測確認（T11b/T11c 含む）。
   結果として compat ファサードの欠陥は見つからず — 事後的に主張は正しかった。
2. **compat 切替で 2D 実行経路は CpuScalar になり V1 比で遅くなる**（要 triage）:
   compat ファサードは `Solver<D2Q9, T, CpuScalar, LocalPeriodic>` 固定
   （V1 ビット一致の根拠）。CLI/GUI の 2D は V1 融合カーネル → CpuScalar への
   置換になり、実測 `lbm presets run cavity` は 140 → 52 MLUPS（2.7x 減）。
   wasm も同様（V1 シリアル融合 → シリアル CpuScalar、俯瞰値で ~5x 減の見込み）。
   対処はファサードの backend を CpuSimd に差し替える 1 行だが、これは軌道を
   ulp レベルで変える挙動変更（backend_simd_equiv のゲートは f64 1e-11 /
   f32 1e-6）なので、本引退作業では**行わず**現状維持。複製スイート緑のまま
   差し替え可能なことは backend_simd_equiv が示唆しており、別途サインオフで。
3. **Shan-Chen 壁吸着の V2 ネイティブ配線完了**（M-D 申し送りの解消）:
   `Solver::update_shan_chen_force_with_walls(g, g_wall, psi_wall, psi)` を追加
   （`MpiSolver` にも同名ラッパ）。solid 隣接は cohesion 和に ψ_wall を寄与し
   `g_wall` の吸着項を加算、非周期域外は無寄与 — V1 と演算子順まで同一。
   受入: `t13_shan_chen_wall_adhesion_native_matches_compat_and_split`
   （3 ケース g_wall=-1.5/+0.9/wall_rho=1.2 × 150 step、native vs compat
   ビット一致 + 2x1/1x2/2x2 分割不変ビット一致）。中立壁の既存呼び出しは
   歴史的式を保持（ビット同一）。

## 新規（2026-07-05 M-D MPI 分散実装）

1. **T13-MPI 全 PASS（場はビット一致）**: mpirun -n {1,2,4} × {2D TGV/キャビティ
   （蓋が縫い目跨ぎ）/縫い目上円柱+プローブ+放物線流入/Shan-Chen 液滴（2×2 コーナー、
   ψ を exchange_scalar 経由）} と -n 8 × {3D TGV 24³ 2×2×2} で、rank-0 gather 場
   （rho/u/全 f 平面）が単一ランク基準と **max|Δ| = 0.0**。診断（mass/momentum/
   probed_force/NaN 数）は rank 部分和 → Allreduce の f64 再結合差のみ
   （≤9.1e-13 abs、液滴 mass は ≤4.5e-11 abs = 相対 ~3e-14）。判定線は T13 流儀
   atol+rtol 各 1e-12（場）/1e-11（診断）。再現: `./scripts/test_mpi.sh`。
2. **Shan-Chen V2 ネイティブ API ギャップ解消**（codex order #6 記載分）:
   `Solver::update_shan_chen_force`（単成分、ψ ハローを exchange_scalar で配線）を
   追加。InProcess の 2×2 コーナー液滴 T13（`t13_shan_chen_droplet_native_split_
   invariant`）もビット一致で緑。壁吸着（g_wall/wall_rho）は未配線 — 必要になった
   時点で compat::ShanChen から移植する。
3. **rsmpi/Open MPI の罠**（詳細 docs/MPI_GUIDE.md）: (a) x86_64 Homebrew MPI が
   PATH 先頭だと rsmpi ビルド/実行が壊れる（arm64 版を先頭に）。(b) 複製
   コミュニケータを持つ MpiSolver を Universe drop（MPI_Finalize）後に drop すると
   MPI_Comm_free で abort（exit 14）— bench_mpi.rs で実際に踏んだ。(c) マスク編集は
   collective: set_solid を所有ランクだけで呼ぶと exchange_masks の呼び出し回数が
   ずれてデッドロック（MpiSolver は非所有ランクも dirty マークを立てて回避）。
4. **弱スケーリング（単一ノード・共有メモリ経由の参考値）**: 512²/rank 直列
   バックエンドで n=1: 40.2 / n=2: 79.9 (99.4%) / n=4: 155.9 (97.0%) /
   n=8: 235.5 MLUPS (73.2%)。n=8 の低下は M5 Max の異種コア（6 Super + 12
   Performance）+帯域競合が主因: 通信ゼロの対照実験（独立 1 ランク×8 並走）でも
   84% 相当が天井で、MPI 化の追加損は ~12%（ロックステップのジッタ結合）。
   n≤4（均質コア内）は R3 ローカル線 ≥85% を満たす。真の測定はクラスタ待ち
   （測定リスト: docs/MPI_GUIDE.md §クラスタ）。

## 新規（2026-07-05 T15.5 3D cavity Re=1000）

1. **T15.5 既定スイートは N=64 qualitative sentinel に固定**:
   `cargo test -p lbm-core --release --test t15_5_cavity3d` は 47.32s wall で
   green（2 passed / 2 ignored）。N=48 は Re/(N-2)=21.7 で 20k step 内に NaN
   発散し、docs/T15_5_CAVITY3D_REFERENCE.md の Re/(N-2) ≲ 15 安定性警告と整合。
   N=64 は同制約をわずかに超えるが、20k step で mass_rel=1.2e-16、
   symmetry-plane max|v|/U≈2e-15、定性的 extrema signs/locations は通るため、
   default では profile 数値帯を要求しない。
2. **T15.5 N=72 spec-profile は red のまま凍結**:
   `cargo test -p lbm-core --release --test t15_5_cavity3d \
   t15_5_cavity3d_re1000_profiles_n72 -- --ignored --nocapture` は 1477.27s wall。
   steady=true at 324500 step、mass_rel=2.546e-15、midplane max|v|/U=1.700e-15、
   anti-2D RMS/U=0.1031、profile RMS/U は u=0.0153（limit 0.030）、
   w=0.0255（limit 0.035）で通過。失敗点は extremum band:
   u_min=-0.25084 at z=0.12925 vs A&K -0.2803833 at z=0.12419、rel=0.105
   （limit 0.06）。w_min=-0.39537 at x=0.90383、w_max=0.22148 at x=0.11181
   も A&K より浅い傾向。従って N=72 の中心線形状は合うが、渦強度は
   A&K/Ben Beya band より数値拡散側で、ignored validation は red evidence として保持。
3. **Endpoint sampling correction**:
   A&K/Ghia 型 17 点表の端点は境界条件値そのものなので、T15.5 sampler は
   u(z=0)=0, u(z=1)=U, w(x=0)=w(x=1)=0 を直接返す。隣接流体セルを端点として
   使うと N=72 で u-line RMS/U が 0.0374 まで悪化し、half-way moving-wall
   境界層を参照端点と混同する。

## PM 回答（2026-07-05 深夜）— レビューセッション判断依頼 4 件 + M-F 統合

- **(a) 仕様書の main 取込**: PM 実施済み（コミット 5cf7a97）。SOLVER_IMPROVEMENT_SPEC.md
  冒頭に main 用パス翻訳注記を追加。scripts/spec-experiments はパス翻訳
  （lbm_core2→lbm_core、V1→compat）のうえ **E2/E7 が改名後 main で仕様書の数値と
  厳密一致再現**することを確認済み。R-Phase 1 セッション側での取込作業は不要。
- **(b) R-Phase 2 の発注時期**: R-Phase 1 着地直後に発注する。M-E 前提に加え、
  M-F（REQ-M-F-STR rev.1b）の構造前提 = **複数分布セット（相場 g・スカラー h）・
  per-cell 物性場・Lagrangian バッファ**を B-1 の設計要求に追加した（PLAN.md 現行キュー参照）。
- **(c) D-6**: PM 直轄で適用済み — COMPETITIVE_SPEC R1/R3 を改定履歴付きで更新
  （球 ±10%・D_h 正規化、弱スケ n≤4 局所線）、PLAN の「R1/R3 達成」表記に注記、
  VALIDATION T15 の ±25%/±15% 不整合も解消（= A-10(c) の文書側は処理済み。
  R-Phase 1 エージェントはコード側の t15_3d.rs コメントのみ対処すればよい）。
- **(d) codex D-8（T14/T15 敵対発注）**: R-Phase 1 着地後に発注（入口ガードで
  不正構成の挙動が Err に変わるため、ガード後の仕様で攻撃させるのが正しい）。
  T15.5（3D キャビティ A&K 2005）は別途 codex order #7 として実行中。
- **R-Phase 1 の起動**: チップ task_f890716a の押下は不要 — PM が worktree
  `/Users/taku/projects/lbmflow-wt-rphase1`（ブランチ r-phase1、main 5cf7a97 ベース）で
  Opus に発注済み。スコープ A-2〜A-10（D-6/D-7 除外、A-1 残作業は現地判断）。
- **CpuSimd 切替**: 引き続き保留（B-1/B-2 の同期点契約整理後に判断）— 貴見解と一致。
- **M-F 統合完了**: REQ rev.1b（表題中立化・コア改名追随・§7 メモリ予算表・T17 配線）、
  VALIDATION.md に T13/T14 節（D-7）+ T16 プレースホルダ + **T17（VR-STR-01〜07）**新設、
  PLAN.md に R-Phase キューと MF-α〜ζ 実装トラック表を制定。
  残仕様詰め: active スカラー帰還式はリサーチ委任中（→ docs/proposals/）、
  REQ 第 2 次 codex 検証は rev.1b に対し発注。

### 進行中プロセスの注意（2026-07-05 深夜時点）
- codex #7 が main ツリーで実行中。副作用として `cargo fmt` がソース 17 ファイルに
  整形のみの diff を生成している（lattice.rs / kernels.rs を目視確認 — 意味変更なし）。
  codex 完了時に PM が triage（テスト成果物以外の整形 diff は revert 予定）。

## codex REQ 第 2 次レビュー triage（2026-07-05 深夜、PM）

11 件（Critical 1 / Major 6 / Minor 4）**全採択** → REQ rev.2 として適用済み。
所見原本: docs/proposals/req-round2-findings.md。要点:
- C1: 緩和拡張の忠実度基準比検証が未配線 → **VR-STR-RELAX 群を新設**（REQ §8 + T17。
  初版は trait/スキーマ/検証項目の予約のみ、帯は緩和実装時に凍結）。
- M2: 「一括実装」と「後付け拡張」の混線 → 納品スコープを「忠実度既定 = 一括実装、
  緩和 = API 予約のみ」に明文化。
- M3: 変数 σ 時の表面張力規約 → §3 を条件分岐化（σ 一定 = μ_φ∇φ 基準形 / active =
  well-balanced 併用形。**係数導出が実装前必須** — リサーチ提案の要導出と同一項目）。
- M4: F_b^scalar（Boussinesq）を運動量式・FR-COUP-01 力源合成に追加、
  C≡C_0 厳密ゼロ + VR-STR-06+ 退化検証。
- M5/M6: T17 転記漏れ（エネルギー様量 = 監視量扱い）と 02a/b/c 分割を REQ/T17 両面で修正。
- M7: NFR-02 の旧「f32 既定」語彙を忠実度プロファイル既定に統一。
- m8/m9: メモリ予算表の界面帯償却 +18–37 B に修正（帯 5–10% と整合）、
  「数 TB 級」→「0.6 TB 級（既定）/ TB〜数 TB（全 f64・複数スカラー・CP 込み）」。
- m10/m11: §2 見出しの中立化整合、M-F 忠実度既定 = 常に D3Q27 の明示。

## PM integration record (2026-07-05, late night — English from here on per user directive)

- **Language policy change (user directive)**: ALL artifacts in English going forward
  (code, docs, commits, UI/CLI strings). CLAUDE.md rule updated. A dedicated spawned
  session translates all legacy Japanese content (docs/*.md, TESTING_NOTES, GUI/CLI/
  wasm strings). Until it lands, documents are transitional mixed-language.
- **REQ rev.3 applied**: competitive-review triage diff (authored as "rev.1c" against
  rev.1b by the requirements session) merged on top of rev.2. P1 population balance
  (scope-aligned to point-bubble relaxation), P2 §4.8 FR-EXT-01 extension contracts
  (co-designed with R-Phase 2 B-1), P3 FR-IO-05 (blend time/RTD) + FR-IO-06 (parallel
  I/O, deterministic checkpoint — converges with spec B-5/C-3/C-8), P4 reference
  datasets (names only; bands stay experiment-frozen), P5 product-layer out-of-scope
  note, §11 implementation dependency DAG (W-items, 6-way wave-1, two critical paths).
  Boundary decisions upheld: no KPI duplication (CLUSTER_OPTIONS owns R3), no
  hardcoded thresholds, product ecosystem in separate volumes.
- Note to the requirements session: your "rev.1c" landed as **rev.3** because rev.2
  (codex round-2, 11 findings, all adopted) had already been applied on main. No
  content of yours was dropped; scope-alignment notes were added where rev.2's
  "fidelity-default = initial delivery, relaxations = API-reserved" decision
  interacts with P1/W-BUB.

## External review (REV-CFD-*, filed vs rev.1a) — PM triage → REQ rev.4 (2026-07-05)

All 14 findings critically verified against the CURRENT document (rev.3, which the
reviewer had not seen). Dispositions:

| ID | Verdict | Disposition |
|---|---|---|
| CR-001 sparger phase inversion | **valid bug** (φ=1 ban read as liquid-injection ban) | ADOPTED: FR-VOF-03 rewritten — gas inlet = φ=0, `inlet_phase: gas\|liquid` in schema (raw φ never exposed), volume-balance acceptance |
| CR-002 AC/continuity mass-flux inconsistency | **valid** — with ρ=ρ(φ) and diffusive J_φ, naive continuity fails at ratio 10³ | ADOPTED: consistent/AGG-type formulation normative (J_ρ=(ρ_l−ρ_g)J_φ in continuity AND momentum advection, same discrete path), droplet advection test |
| CR-003 neq stress stage/coefficient mismatch | **valid residual** of codex round-1 #2 fix — post-collision stage stated with pre-collision coefficient | ADOPTED: default = pre-collision/post-streaming; explicit BGK (1−1/τ) / MRT R(τ)⁻¹ transforms; required `neq_stage` enum; stage cross-check test |
| CR-004 forcing 2nd-moment sign contradiction | **valid** (prose "subtract" vs formula "+") | ADOPTED: Π_neq_raw/Π_force/Π_neq_corr single-equation definition; prose sign words banned; sign derivation-frozen pre-implementation + negative test (body-force Poiseuille) |
| MJ-005 Ca_spurious dimensional | **valid** (stray L) | ADOPTED: Ca_spurious = μ_l\|u\|/σ; Re_spurious separate. VALIDATION T17 synced |
| MJ-006 Pe vs U_tip | **valid** (π ambiguity) | ADOPTED: Pe_N = Re·Sc / Pe_tip = π·Re·Sc split; bare "Pe" banned |
| MJ-007 active scalar 1-step lag | **valid** vs fidelity-default principle | ADOPTED: dataflow split passive/active; predictor–corrector default; `active_scalar_lagged` = flagged relaxation via VR-STR-RELAX; dt-halving acceptance |
| MJ-008 一括 vs 後日 | already fixed in rev.2 (codex round-2 #2) | STRENGTHENED: explicit Initial-delivery / Phase-2 lists added to §0 |
| MJ-009 f32/f64 boundary undefined | **valid** (needed for array/GPU design now) | ADOPTED: precision_profile enum {full_f64, mixed_safe(default), mixed_fast}; interface_band = max(3W,6Δx) provisional, re-frozen at W-VOF |
| MJ-010 no numeric thresholds | conflicts with characterize→freeze protocol; concern (post-hoc band-fitting) legitimate | **ADAPTED**: provisional numeric bands added NOW (Np ±10% etc.) + asymmetric governance — tighten freely, loosen only with PHYSICS.md rationale (T15.5 precedent). Reviewer's per-test metadata format adopted for T17 rows |
| MJ-011 scalar non-conservative form | **valid** for two-phase/active | ADOPTED: phase-wise conservative + ρY forms normative; Henry flux sign convention; total-mass conservation test |
| MJ-012 four-way contact undefined | **valid gap** | ADOPTED as Phase-2 contract: FR-PART-04 (soft-sphere params), -05 (lubrication), -06 (config rejection beyond two-way regime — ships in initial delivery) |
| MJ-013 viscosity interp / σ coefficient hedges | **valid** ("固定版" claim violated) | ADOPTED: harmonic-in-μ default frozen (alternatives = logged options outside default bands); "(係数はモデル定義)" hedge removed — σ=√(2κβ)/6, W=4√(κ/(2β)) are THE definitions (internal consistency was verified by codex round-2) |
| MN-014 ε_g processing units | **valid refinement** | ADOPTED: ε_g_raw / ε_g_thresholded(φ_c=0.5) / kernel-smoothed / hybrid-dedup definitions + mandatory metadata |

Net: 13 adopted (1 adapted), 1 already-fixed-and-strengthened. REQ is now rev.4.
Reviewer read rev.1a — overlaps with rev.2/rev.3 noted above to avoid double-fixing.

## D-5 validation horizon (2026-07-05)

- Added `crates/lbm-core/tests/d5_long_horizon.rs`.
- Native `Solver<D2Q9, f64, CpuScalar, LocalPeriodic>` TGV convergence, default suite:
  `cargo test -p lbm-core --release --test d5_long_horizon d5_native_solver_tgv_converges_second_order -- --nocapture`
  measured `e32=2.622406e-3`, `e64=6.982198e-4`, `order=1.909`.
- Ignored long-horizon Re=100 cavity compat facade vs native Solver:
  `cargo test -p lbm-core --release --test d5_long_horizon d5_cavity_re100_compat_matches_native_after_20k_steps -- --ignored --nocapture`
  measured `rho=0.000000e0`, `ux=0.000000e0`, `uy=0.000000e0`, `worst=0.000000e0`
  after 20,000 steps on a 129x129 f64 TRT cavity. The frozen assert is `worst <= 1e-12`,
  below the D-5 `1e-9` ceiling.
## D-4 f32 x 3D validation measurements (2026-07-05, branch cx-d4)

New default-suite test file: `crates/lbm-core/tests/t15_3d_f32.rs`.

- T15-1 f32 z-invariant TGV degeneracy, D3Q19 `32x32x4` vs D2Q9 `32x32`,
  `nu=0.02`, `u0=1.28/N`, 648 steps: max relative agreement on the characteristic
  velocity scale is `4.400e-6` (`rho=5.958e-7`, `ux=4.400e-6`, `uy=3.795e-6`,
  `|uz|/u0=1.164e-8`). Test gate: `<= 1.0e-5`.
- T15-4 f32 TGV3D decay rate, D3Q19 `64^3`, `nu=0.02`, frozen scaling
  `u0=1.28e-4/N=2.000e-6`, 519 steps: measured rate `1.155265e-3`,
  diffusive reference `1.156594e-3`, relative error `1.149e-3`.
  Test gate: `<= 2.0e-2`.
- T15-4 f32 TGV3D mass drift, same `64^3` setup, 1000 steps: `m0=2.621440000e5`,
  `m1=2.621440000e5`, relative drift per 1000 steps `3.109e-15`.
  Test gate: `<= 1.0e-5`.

Note: `lbm-scenario::Sim3Handle::F32` is a thin wrapper around
`Solver<D3Q19, f32, CpuScalar, LocalPeriodic>`. These tests live in `lbm-core`,
so they pin that product engine type directly without adding a reverse
dependency from core tests to the scenario crate.
## D-11 wasm smoke record (2026-07-05, branch cx-wasm-smoke)

- Added a wasm-bindgen-test smoke in `crates/lbm-wasm` using a test-only
  Taylor-Green JSON initializer on the existing `WasmSim::init` JSON path:
  32x32, nu=0.02, BGK, periodic edges, u0=1.28/32, 100 steps.
- Native f32 characterization for the same compat path:
  - rho-view mass sum before: 1023.999993563
  - rho-view mass sum after 100 steps: 1023.999934435
  - relative mass drift: 5.7741999989150555e-8
  - frozen probe at (7, 11) after 100 steps:
    rho bits 0x3f8025b6, ux bits 0xbbb62bd2, uy bits 0xbc98d05a.
- `wasm-pack test --node crates/lbm-wasm` result: PASS
  (`tests::wasm::wasm_tgv_smoke_matches_compat_f32` passed; wasm rho/ux/uy
  views matched the compat f32 run bit-for-bit, and velocity views had no NaN).
- `wasm-pack build crates/lbm-wasm --target web --release --out-dir ../../web/src/engine/pkg`
  initially failed after Rust wasm compilation at wasm-pack's external optimizer/helper install
  step with: `Operation not permitted (os error 1)` and wasm-pack's hint
  `To disable wasm-opt, add wasm-opt = false to your package metadata in your Cargo.toml`.
  The crate now sets `[package.metadata.wasm-pack.profile.release] wasm-opt = false`.
  With `XDG_CACHE_HOME=/private/tmp/lbmflow-wasm-pack-cache`, the same build command passed.
- Cargo registry/network note: this sandbox cannot resolve crates.io/static.crates.io. The current
  `wasm-bindgen-test` release has a target-gated `minicov` coverage dependency in its lock graph;
  a local `minicov` stub is patched in under `crates/lbm-wasm/test-support/` so metadata resolves
  offline. The stub is not compiled for the normal wasm smoke.
- Added a Rust-only `lbm-scenario` test for the GUI-exported scenario JSON shape; it parses, builds,
  reserializes, reparses, and serializes byte-stably without node/web tooling.

## PM record — B1 approval, order-A/B/C triage, dispatch lesson (2026-07-05 late)

- **B1 capability map APPROVED and merged** (docs/skills/b1-capability-map.md, one-file
  branch). Highlights: 7 MCP tools empirically confirmed (async path driven end-to-end);
  **BUG found: explicit 2D backend:"gpu" silently runs on CPU** (status:completed,
  validate ok:true, no warning) → fix order dispatched (branch cx-gpu-fallback-guard:
  honored-or-error for explicit backend requests; "auto" may fall back by design).
  Other reds: no unit->lattice conversion anywhere; no user-facing accuracy-compare
  command (validation is cargo-test only); 3D limited to single-phase + init:rest + CPU.
  B2 session launched on the approved map (branch skills/b2).
- **Bioreactor session's follow-up orders triaged**: §1 body-force API = already in
  trunk (guard suites green). Order A (strain-rate observable per FR-STRESS-01) =
  ACCEPTED, dispatched (branch cx-strain-rate; W-STRESS pulled forward — hard dep is
  W0 only per REQ §11). Order B (moving no-slip boundary) = DEFERRED to MF-δ; their
  adversarial test seeds recorded: translating flat-plate drag vs analytic,
  Taylor-Couette interior azimuthal profile, mass conservation across mask motion,
  partition invariance. Order C (raster lift) = queued behind A.
- **Dispatch lesson (feeds lbmflow-codex-dispatch Skill v2)**: an inline codex order
  containing backticks dies in zsh command substitution (parse error near ')').
  Robust invocation: write the order to a file and pass "$(cat <file>)" — the
  substituted string is NOT re-parsed. The Skill's invocation section should make
  file-passing the default for any order containing backticks/code spans.

## Parity harness smoke SM-1 (2026-07-05 late) — harness defect found and fixed

CD-HO-01 on Sonnet: evaluee REFUSED — flagged the external fixture-file trust hop as a
prompt-injection pattern (defensible) and noted the hypothetical task IDs don't exist
in the repo. Meanwhile it had read lbmflow-codex-dispatch and cited the CD-3 same-file
bundling rule correctly — the Skill content reached the model; the harness framing
failed. Protocol amended (runner preamble v2 on branch skills/a-pilot-eval-tasks):
fixtures inlined into the prompt, exercise declared self-contained/hypothetical,
refusal-handling rule added. Full 96-run parity batch deferred to a dedicated
orchestration session with the v2 preamble; smoke rerun first.

## PM record — B2 approved & merged (2026-07-05 late)

Five green user Skills merged (.claude/skills/lbmflow-user-{run-preset, author-scenario,
tune-stability, run-monitor-mcp, postprocess}) + docs/skills/b2-skill-specs.md.
PM answers to B2's open questions: (1) obstacle-composition FOLD approved, no Y1 order;
(2) no defensive 2D-gpu warning line — the honored-or-error fix (cx-gpu-fallback-guard)
lands first and gates the Skill's assumption; (3) unit conversion stays routing-only
until W-UNIT (REQ §11) delivers the feasibility layer — the user-facing converter is
spec'd together with it; (4) run-preset / run-monitor-mcp split accepted, no 6th Skill.
Parity evaluation for user Skills follows the A-pilot protocol once that pipeline is
validated (runner preamble v2, smoke rerun pending).
## New (2026-07-05 R-Phase 1: entry guards A-2..A-10, branch r-phase1)

Written in English per the 2026-07-05 language policy (new notes English-only).
Engine-side changes that alter *rejection* semantics — adversarial tests
(codex order #7+) should target these seams. No numerical path changed:
probe_state_hash-equivalent bit invariance of legal configurations is pinned
by `healthy_run_is_ok_and_bit_identical` (run vs run_guarded) and the
untouched T13/T14/backend_simd gates.

1. **A-2/A-6 (compat `SimConfig::build`)**: NaN/Inf in edge velocities, body
   force, or TRT magic now `Err(NonFiniteParameter)` / `InvalidParameter`
   (speed test reversed to NaN-safe `!(s <= MAX_SPEED)`). A MovingWall with a
   wall-normal velocity component is rejected (E7: silent -56% mass / 500
   steps). `MAX_SPEED` moved to `params::MAX_SPEED`; compat re-exports it.
2. **A-7 (compat `init_with`)**: panics with the offending coordinate on
   rho <= 0, non-finite rho, or speed > MAX_SPEED. Closure purity documented
   (re-evaluated up to 5x per cell by the FD stencil).
3. **A-3 (compat/wasm `set_solid`)**: placing a solid on the cell directly
   inward of an open edge (x==1 / x==nx-2 / y==1 / y==ny-2 for open
   left/right/bottom/top) now panics — that neighbour feeds the open-face BC;
   a solid there froze the unknown slots (E5b: permanent ux=-0.115, no NaN).
   Non-panicking pre-check: `Simulation::set_solid_allowed(x, y)`. The wasm
   paint tool refuses such strokes silently.
4. **A-4 (`GlobalSpec::validate(d, solid) -> Result<(), SpecError>`)**: the
   V2-native gate. Rejects: non-finite/non-positive nu; bad TRT magic;
   non-finite force and (2D) force[2] != 0; active axis < 3 cells; periodic x
   open on one axis; open faces on more than one axis; a non-periodic face
   that is neither open nor a full solid rim (E2); inlet speed > MAX_SPEED
   (NaN-safe); outlet rho <= 0; u_conv outside (0,1]; open-face axis < 3.
   `Solver::build` enforces it (panic, defense-in-depth); lbm-scenario
   `build3d` calls it and maps `SpecError` -> `Build3Error::Spec`. Scenario
   keeps only: periodic *pairing* (two EdgeSpecs -> one bool), the
   `AdjacentOpenEdges` kind its guard test pins, and MovingWall speed (wall
   velocities live in WallSpec, invisible to GlobalSpec).
5. **A-5 (`HaloExchange::SCOPE`)**: `Solver::new_local_part` (single-part
   owner) now requires a `Remote`-scope exchange at construction;
   LocalPeriodic/InProcess are `Local` and panic (E4: silent self-wrap,
   rho off by 7.7e-2). MpiExchange declares `Remote`.
6. **A-8**: `zou_he_face_3d` asserts `unknowns(face).len() == 5` (D3Q27
   would otherwise silently skip 4 slots). New `tests/stream_contract.rs`
   pins the ConvectiveOutflow memory-term contract: streaming must not write
   open-face unknown slots — CPU: sentinel bits unchanged across a full
   stream pass (D2Q9 4 faces, D3Q19 6 faces); GPU: 200-step channel agrees
   with CPU on the outflow-face unknowns <= 1e-4 (M5 Max/Metal, green).
7. **A-9 (`run_guarded(steps, check_every)`)**: standard watchdog on
   Solver / GpuSolver (readback) / MpiSolver (collective 2-double
   Allreduce). NaN/Inf caught at the next check with the step number.
   Overhead at 512^2, check_every=100: 0.45-0.49% (per-check 1.6-1.9 ms vs
   per-step 3.6-3.9 ms; ignored test asserts the <1% line on the component
   ratio — end-to-end timing is machine-noise-dominated on a shared box).
   CLI drivers still use their own rho-scan (behaviour pinned by runner
   tests); rewiring them onto run_guarded is a PM follow-up.
8. **A-10f**: `equilibrium()` vs collide's inline feq pinned to bit identity
   via the fixed-point property (equilibrium state must survive forceless
   collision bit-exactly), D2Q9/D3Q19 x f32/f64.

## Stirred-tank demo — measured behavior (MF-δ precursor, 2026-07-05)

3D baffled Rushton stirred tank, `crates/lbm-cli/examples/stirred_tank_3d.rs`
(kept UNTRACKED per PM until the raster/product framing is resolved). Ran on the
primary checkout `feat/body-force-field-api @ d7c4053` — NO branch switch. Backend
CpuScalar, D3Q19, TRT (MAGIC_STD). The impeller is volume-penalization, NOT a
resolved moving solid: a Guo body force (public `set_body_force_field`, b74298e)
drags turbine-footprint cells toward `v = omega x r`. This is the sanctioned interim
before IBM-inertial (REQ §4.3 FR-ROT-01 / W-ROT, MF-δ). Baffles + round wall are
true no-slip solids (half-way bounce-back). Shear here is an EXAMPLE-SIDE finite-
difference proxy `nu*sqrt(2 S:S)` (central diff) — the core exposes no strain-rate
field yet; replace with the non-equilibrium-moment field when order A / FR-STRESS-01
(branch cx-strain-rate) lands, then re-baseline shear_max.

Config (n^3 default 80): tip_r=12.33 (D=T/3), 6 blades, 4 baffles, spin-up ramp
1500 steps, penalization gain alpha=0.32, force cap 0.02. Added backward-compatible
CLI args `u_tip` (arg5) and `nu` (arg6) + a SUMMARY line + divergence early-break
for the sweeps below (defaults 0.08 / 0.02 unchanged).

Measured (n=48 fast sweep + n=80 reference/edge):
- **omega / Ma_tip is NOT the binding limit.** STABLE across the whole u_tip sweep
  0.04..0.20 (Ma_tip 0.069..0.346) at nu=0.02, and at Ma_tip=0.277 (u_tip=0.16,
  nu=0.01). `final_max|u| ~= 0.96*u_tip` (penalization reaches ~96% of rigid tip
  speed). Compressibility error ~O(Ma^2) is the real cap: recommend Ma_tip <= 0.1
  (u_tip <= ~0.058) for quantitative use; default u_tip=0.08 (Ma_tip=0.139) is a
  visualization compromise (~2% Ma^2 error).
- **tau / Re edge IS the binding limit** (80^3, no SGS model): STABLE down to
  tau~=0.507 (nu=0.0025, Re~789); **DIVERGES at tau~=0.504 (nu=0.00125, Re~1579)**,
  max|u| -> 2.4e6 at step ~2500. Practical envelope: tau >~ 0.51 (nu >~ 0.0025),
  Re <~ ~800 at 80^3 without a subgrid model. Above this needs W-LES/cumulant —
  concrete motivation for REQ risk #2 (§4.2) before any high-Re stirred run.
  NB the default 80^3 config is Re~99 (laminar); Re scales with n (tip_r).
- **Shear-field sanity**: monotonic with u_tip and nu in the stable regime;
  explodes to shear_max=1772 at divergence (clean blow-up signature). Spatially
  correct: six blade-tip shear lobes decaying into the bulk (textbook Rushton
  discharge), velocity mirrors it with a six-lobe radial jet, baffles break the
  swirl. Reference (80^3, u_tip=0.08, nu=0.02, 4000 steps): speed_max=0.0767,
  shear_max=6.3e-4; PNG slices + subsampled volume.bin/json emitted.

Feeds MF-δ: penalization gives the right qualitative discharge/shear topology and a
bounded, well-characterized stable envelope; the tau-floor divergence at Re~1579 is
the numeric evidence that W-LES must precede high-Re stirred validation. Next: swap
the FD shear proxy for FR-STRESS-01 once cx-strain-rate lands, then IBM-inertial
(W-ROT) supersedes penalization for torque/Np fidelity.
9. **A-1 residual**: not needed — no AUTO-GENERATED headers remain under
   crates/lbm-core/tests/ (sync-tests.sh deleted; suites are compat-native).
   A-10a/b: not applicable on main (V1 deleted; facade carries neither the
   unused_mut nor the misleading solid-rho comment).

## Triage record — "flickering particles" report (2026-07-05, viewer session)

Reported visual artifact in the 3D stirred-tank viewer investigated by that session:
**NOT a solver defect, no core change**. Field evidence (60^3 subsampled export from
the D3Q19 n=80 TRT baffled-tank run, volume-penalization impeller on
set_body_force_field): no NaN/Inf; |u| decays smoothly from the impeller plane
(mean 0.021 / max 0.075 at simZ~27) to a near-quiescent headspace (mean ~1e-4,
max ~4e-3 at simZ>55); the small top-layer residual is the shaft's swirl sheath
(penalization region runs full height) — physical. Root cause was viewer tracer
reseeding + hard cull threshold; fixed viewer-side. Optional modeling note recorded:
cap shaft forcing a few cells below the lid if a dead headspace is ever wanted.
Value for MF-δ: first end-to-end field-sanity pass of the penalization-impeller
pipeline on the new body-force API.
## W-STRESS strain-rate observable (2026-07-05, branch cx-strain-rate)

- Implemented native `Solver::gather_strain_rate()` and `Solver::gather_shear_rate()` for CPU-backed D2Q9/D3Q19 solvers, plus MPI rank-0 gather wrappers and GPU CPU-readback facade wrappers.
- Verified FR-STRESS-01 rev.4 force correction sign in the native body-force Poiseuille check: `Pi_force = -0.5 * (uF + Fu)`, so `Pi_neq_corr = Pi_neq_raw + 0.5 * (uF + Fu)` for this engine's physical velocity and deviation-form equilibrium. Measured interior `gamma_dot(y)` max absolute error: `1.528708800171974e-14`.
- Plane Couette half-way-wall check: analytic `S_xy = 0.5 * U / H` with `U=0.1`, `H=8`; measured max absolute `S_xy` error: `1.3877787807814457e-16`. `gather_shear_rate()` matched `sqrt(2 S:S)` from the returned tensor exactly in this fixture.
- InProcess decomposition check `[1,1,1]` vs `[2,2,1]`: `gather_strain_rate()` and `gather_shear_rate()` matched bit-for-bit after 25 forced BGK steps with a spatial body-force field.

## Attribution correction (2026-07-05, per the requirements session / Taku)

docs/REQ_STIRRED_REACTOR.md + the orders A/B/C derivation = the requirements session
(correct as recorded). Commit b74298e (body-force field API) + the primary-checkout
switch to feat/body-force-field-api = a DIFFERENT worker (earlier PM notes attributed
both to one session — corrected). The stirred-tank demo is now owned by the
requirements session: builds against trunk (R-Phase 1 guards active), demo example
stays untracked, volume-penalization interim, measured behavior lands here in
English for the MF-δ record; no primary-checkout branch switches.
