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
