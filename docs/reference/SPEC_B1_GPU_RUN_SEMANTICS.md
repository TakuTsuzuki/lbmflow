# SPEC: B-1 completion — restore GPU chunked-execution semantics + fix the T14 mixed-BC defect

Status: dispatch-ready (2026-07-06). Owner: PM. Implements in the existing `r2-b1` worktree
(branch `r2-b1`, commits `ac1c8db` + `55dbccb` already present).

## 1. Context and evidence

B-1 stages 1+2 (Backend `Fields` generalization, `stage_in`/`stage_out` staging contract,
GpuSolver routed through the unified `Solver` orchestrator) are implemented and committed.
Gates green: workspace suite, T14 8/8, run_guarded-GPU, t13_split_invariance and
backend_simd_equiv byte-unmodified.

Two defects block landing; a third pre-existing GPU defect is bundled here because it
lives in the same files.

### Defect 1 (introduced by B-1): the GPU run path no longer executes work inside `run()`

Evidence (PM-measured, outside the codex sandbox, this machine, branch `r2-b1`):

```
| grid      | GPU MLUPS   | sane value (pre-B-1) |
| 512x512   | 7,514,878   | ~11,500              |
| 1024x1024 | 31,703,649  | 6,846                |
| 2048x2048 | 118,214,121 | ~5,900               |
```

MLUPS values 3-5 orders of magnitude above the memory-bandwidth ceiling are physically
impossible: the bench's timed region (encode+submit+wait per chunk, C-9 contract) is
returning without the GPU actually executing the steps. T14 stays green because `gather_*`
forces a device sync, so *correctness* is intact — only *when* work executes broke.

Consequences beyond the bench: C-9's auto-calibrated submit chunking (1 submit ≈ 100–250 ms,
TDR safety) is effectively bypassed; unboundedly deferred work risks device-lost on slow
GPUs and breaks any wall-clock-based progress reporting (MCP `run_status`).

### Defect 2 (environment, document-only): adapter not found inside the codex sandbox

`bench_gpu --gpu-only` reports `no usable GPU adapter was found` under the codex sandbox
while the same binary run by the PM finds `Apple M5 Max / Metal`. Cause: sandbox denies
Metal device acquisition for this invocation path. This is NOT a code bug. Consequence for
process: GPU bench numbers can only be collected by the PM outside the sandbox; a codex
order must build the bench and leave a `BENCH-PENDING (sandbox adapter)` note instead of
inventing numbers or looping on retries.

### Defect 3 (pre-existing on trunk, found by D-8): mixed force field + moving wall +
convective open faces exceeds the T14 tolerance

Test: `t14_mixed_force_field_moving_wall_and_open_faces` in
`crates/lbm-core/tests/t14_adversarial.rs` (landed on trunk, `#[ignore]`d with a defect
note; obtain in the worktree via `git checkout 24ef513 -- crates/lbm-core/tests/t14_adversarial.rs`
if not yet visible after the trunk merge).
Measured divergence at t=75: |Δux| = 8.516e-5, |Δuy| = 8.418e-5 vs the 1e-5 gate
(|Δrho| = 4.235e-6 is within its 1e-4 pressure gate). The error grows with t (systematic
boundary-cell accumulation, not noise).

Prime suspects, in order (verify, don't assume):
1. The GPU open-face BC pass does not apply the Guo F/2 velocity correction identically
   to the CPU boundary-moments-correction pass when a per-cell force field is active on
   open-face cells (CPU pass order: collide → halo → stream → open-BC → boundary moments
   correction; the GPU fused push kernel must reproduce this exactly at open faces).
2. Moving-wall bounce-back contribution `2 w_i rho (c_i · u_w)/cs²` computed with a
   different rho (pre- vs post-collision, or with/without F/2-corrected u) than CPU.
3. Pass-order divergence when BOTH a moving wall and an open face touch the same cell
   (corner cells) — CPU applies open-BC before the moments correction; check the GPU
   kernel's boundary branch order.

## 2. Required changes

### Stage 1 — restore chunked execution in the unified orchestrator
- `Solver::run()` (and the deprecated `GpuSolver` wrapper) on a GPU backend must submit
  and wait in C-9 auto-calibrated chunks exactly as the pre-B-1 `GpuSolver::run` did:
  after `run(n)` returns, all n steps have completed on the device (the host staging may
  remain stale — that is fine — but the DEVICE work must be done).
- Do not reintroduce a per-step sync; the C-9 chunk calibration (first-chunk measurement,
  target 100–250 ms/submit, cap 200 steps) is the contract. Reuse the existing calibration
  code; do not rewrite it.
- Add a regression test that cannot lie about this: e.g. run(k) on a grid sized so that
  k steps at the roofline take ≥50 ms and assert wall-clock ≥ a floor derived from a
  measured single-step time, OR assert via a wgpu timestamp/on-submitted-work-done
  callback that the queue is idle when run() returns. Pick the least flaky construction;
  print measured values in the assert message.

### Stage 2 — fix Defect 3 (GPU mixed-BC equivalence)
- Reproduce first: un-ignore `t14_mixed_force_field_moving_wall_and_open_faces` locally
  and bisect WHICH pass first diverges (dump per-pass max|Δ| at t=1,2,3 between CPU and
  GPU — a temporary diagnostic, delete before commit).
- Fix on the GPU side only. CPU is the reference; its results must not change
  (backend_simd_equiv and t13 stay byte-unmodified).
- Remove the `#[ignore]` from the test once green. No tolerance loosening anywhere.

## 3. Acceptance

1. `cargo test --workspace --release` green.
2. `cargo test --workspace --release --features gpu` green, INCLUDING the un-ignored
   `t14_mixed_force_field_moving_wall_and_open_faces` and the new chunk-execution
   regression test. T14 8/8 green.
3. `t13_split_invariance`, `backend_simd_equiv` byte-unmodified and green.
4. Bench: build `bench_gpu`; if the sandbox denies the adapter, leave `BENCH-PENDING
   (sandbox adapter)` in the final summary. Landing evidence (PM-run): 1024² within 3%
   of 6,846 MLUPS and all values in a physically sane band (10²–10⁴ MLUPS).
5. TESTING_NOTES.md append (English): per-pass divergence findings, the root cause of
   Defect 3, measured bench values or BENCH-PENDING.

## 4. Rules

One commit per stage, English. Never commit red. Do not touch CPU backends, collide
internals, or docs other than TESTING_NOTES.md. Sandbox git-commit failure → leave
committed-ready and say so.
