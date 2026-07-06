# SPEC: B-1 GPU run() semantics — LANDED

Status: **LANDED** (2026-07-06). This spec is retained as a one-page contract
statement; the operational trap is duplicated in `docs/HANDOFF-PM-2026-07-07.md §9`.

## Contract

On a GPU backend, `Solver::run(n)` submits recorded chunks per C-9 auto-calibrated
chunking (target 100–250 ms/submit, cap 200 steps) and returns after submission.
Device work completion is guaranteed at the next explicit `sync()` /
`gather_*` call, not at `run()` return.

- Tests that assert "`run()` blocks until device idle" are **wrong** — assert the
  submissions counter and a sync fence instead.
- The host staging buffer may be stale after `run()` returns; that is intentional.
- MCP `run_status` progress is driven off submission count, not device wall clock.

## Reference implementation

`crates/lbm-core/src/solver.rs` (GpuSolver path), C-9 chunk calibration in the
GPU backend. Landing regression tests live under
`crates/lbm-core/tests/` (T14 3D + chunk-execution regression).

## History

The original dispatch order (three-defect fix, 2026-07-06) landed via the `r2-b1`
branch; the mixed-BC defect (t14_mixed_force_field_moving_wall_and_open_faces) is
green and un-ignored. Full defect notes in git history / TESTING_NOTES.md.
