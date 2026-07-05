
## Cross-session coordination notes (evening of 2026-07-05)

- **Parallel-session note**: the full solver review (a separate session, memory:
  lbmflow-solver-review-2026-07) has created docs/SOLVER_IMPROVEMENT_SPEC.md v1 on
  branch `claude/amazing-mirzakhani-4060d3` (not yet merged to main).
  The R-Phase 1 implementation session has been commissioned (A-2–A-10 + D-6/D-7).
- **To the R-Phase 1 session**: on the evening of 2026-07-05 main completed V1
  retirement and **renamed crates/lbm-core2 → crates/lbm-core** (V1 removed, compat
  retained). If your branch is based on the pre-rename layout, path remapping is
  needed at rebase time. The sed bug fix in sync-tests.sh (S0) has the same content
  landed on main via 622bbb2 on the v1-retirement side (no double application needed;
  sync-tests.sh itself has also been removed).
- **PM pending item**: switching compat's default backend CpuScalar→CpuSimd
  (2.7x recovery for 2D CLI/wasm) may overlap with R-Phase 1, so decide after that lands.
