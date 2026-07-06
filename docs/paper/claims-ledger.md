# Claims Ledger — paper ↔ implementation convergence (owner strategy 2026-07-05)

Owner ruling (Taku, overrides PM ruling B): the technical paper is written in
present tense as the ideal state (paper = product spec). RELEASE GATE: the paper
stays an INTERNAL DRAFT until every present-tense claim is true; residual deltas
are trued up immediately before release; no external distribution before that.
This ledger is the gate instrument — release is blocked until every row is
green or the claim is edited at true-up.

| Claim (paper, present tense) | Implementing item | Hard gate | Status |
|---|---|---|---|
| 3D GPU acceleration (D3Q19) | ME-1 (B-1 → WGSL 3D + C-13 explicit path) | T14-3D ≤1e-5 + ≥1,500 MLUPS single GPU | **AMBER** — T14-3D GREEN (32³ TGV3D u 2.8e-6, 24³ cavity3D u 1.7e-6); MLUPS RED (peak 1353 @192³, three-sample max; kernel-perf follow-up dispatched) |
| explicit backend:"gpu" runs on GPU | C-13 (bundled into ME-1) | scenario gpu request honored end-to-end, guard flips from error to run | **GREEN** — landed 2026-07-06 (commit 1a14d90) |
| FP16 storage, ×2 grid capacity | ME-2 (C-12) | T16 bands frozen + ≥1.5× MLUPS @2048² | **AMBER** — scaffolding landed (commit 892efdd); T16 matrix ignored, awaiting SHADER_F16 adapter for freeze |
| Multi-node scaling | ME-3 (cluster campaign) | 64-rank weak ≥80% measured | RED — needs cluster access go |
| Full-physics stirred workload | ME-4 (MF-ζ) | degradation ratio vs single-phase published | RED — longest lead; true-up scope risk, flag early |
| 2D GPU GLUPS / CPU MLUPS / T13 bit-exact / weak-scaling n≤4 / wasm bit-identity / agent-native MCP+Skills | landed | measured tonight or earlier | GREEN |

## Aggressive resequencing (~1-week horizon, dispatched via codex-max)
Day 0 (tonight): consolidated gate → push → R-Phase 2 B-1 dispatch (staged orders).
Day 1-2: B-1 stages (Backend Fields generalization, GpuSolver unification) + B-2.
Day 2-4: ME-1 — 3D WGSL kernels + BC passes + T14-3D + C-13 explicit-gpu path.
Day 4-6: ME-2 FP16 (shares kernel plumbing) + bench_gpu 3D numbers.
Parallel: ME-3 prep (bench_mpi 3D/weak modes = C-16b) + cluster spend confirm (user);
campaign executes within the week of access. ME-4 tracks M-F (cannot honestly close
in 1 week — release true-up must scope or wait; flagged to the paper session).
