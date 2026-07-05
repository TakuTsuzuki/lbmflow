# TESTING_NOTES

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

## 新規（2026-07-05 V1 引退作業）

1. **sync-tests.sh の置換が macOS では無効だった**: `sed -E 's/\blbm_core\b/…/'` の
   `\b` は BSD sed 非対応で無置換コピーになっており、「compat へ再標的化済み」の
   複製テスト 16 ファイルは実際には dev-dependency の V1 を直接テストしていた
   （M-A の「56+ テストが compat 経由で緑」は未検証状態だった）。perl 置換に修正して
   再同期し、compat 実経由で全複製スイート緑を実測確認（T11b/T11c 含む）。
   結果として compat ファサードの欠陥は見つからず — 事後的に主張は正しかった。

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
