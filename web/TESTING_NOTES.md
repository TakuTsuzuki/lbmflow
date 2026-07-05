
## セッション間調整メモ（2026-07-05 晩）

- **並走セッション注意**: ソルバー全体レビュー（別セッション、memory:
  lbmflow-solver-review-2026-07）が docs/SOLVER_IMPROVEMENT_SPEC.md v1 を
  ブランチ `claude/amazing-mirzakhani-4060d3` に作成済み（main 未マージ）。
  R-Phase 1 実装セッションが発注済み（A-2〜A-10 + D-6/D-7）。
- **R-Phase 1 セッションへ**: main は 2026-07-05 晩に V1 引退を完了し
  **crates/lbm-core2 → crates/lbm-core に改名済み**（V1 削除、compat 存続）。
  ブランチが改名前ベースの場合は rebase 時にパス読み替えが必要。
  sync-tests.sh の sed バグ修正（S0）は v1-retirement 側 622bbb2 で同内容が
  main に入っている（二重適用不要、sync-tests.sh 自体も削除済み）。
- **PM保留事項**: compat の既定バックエンド CpuScalar→CpuSimd 切替
  （2D CLI/wasm の 2.7x 回復）は R-Phase 1 と重複しうるため、その着地後に判断。
