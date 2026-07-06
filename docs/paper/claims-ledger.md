# Claims Ledger — capability status snapshot

**What this is (2026-07-07 clarification, overrides the 2026-07-05 paper-first
framing):** a working table that maps each notable product claim to its
implementing item and the measurement that verifies it. It is a status
snapshot, NOT a release gate — the paper describes what is measured today and
gets updated when the measurements do. The product goal is to become the best
conceivable LBM simulator; this ledger just tracks where we are on that.

Update rules:
- When a measurement changes, update the Status cell here in the same commit.
- When a row becomes GREEN, edit the paper to reflect it (or add coverage).
- Do NOT hold implementation to a paper claim ahead of the measurement —
  implementation converges to physics/spec, the paper follows.

| Claim (paper, present tense) | Implementing item | Hard gate | Status |
|---|---|---|---|
| 3D GPU acceleration (D3Q19) | ME-1 (B-1 → WGSL 3D + C-13 explicit path) | T14-3D ≤1e-5 + ≥1,500 MLUPS single GPU | **GREEN** — T14-3D GREEN (32³ TGV3D u 2.8e-6, 24³ cavity3D u 1.7e-6); MLUPS GREEN 2026-07-06 quiet-window A/B/A: 192³ 2791-2813, 128³ 2778-2880 MLUPS on unmodified main (the earlier 1353 was a loaded-window artifact; the step_periodic kernel follow-up measured SLOWER and was rejected — TESTING_NOTES 2026-07-06) |
| explicit backend:"gpu" runs on GPU | C-13 (bundled into ME-1) | scenario gpu request honored end-to-end, guard flips from error to run | **GREEN** — landed 2026-07-06 (commit 1a14d90) |
| FP16 storage, ×2 grid capacity | ME-2 (C-12) | T16 bands frozen + ≥1.5× MLUPS @2048² | **GREEN** — 2026-07-06: bands frozen (TGV transient 2e-1 / measured 1.401e-1; cavity steady 5e-3 / measured 2.579e-3; steady-vs-transient dichotomy in PHYSICS.md), MLUPS ~2.0× @2048² (interleaved A/B ×2), D3Q19 f16 >5 GLUPS; capacity ×2 inherent |
| Multi-node scaling | ME-3 (cluster campaign) | 64-rank weak ≥80% measured | RED — needs cluster access go |
| Full-physics stirred workload | ME-4 (MF-ζ) | degradation ratio vs single-phase published | RED — longest lead; true-up scope risk, flag early |
| 2D GPU GLUPS / CPU MLUPS / T13 bit-exact / weak-scaling n≤4 / wasm bit-identity / agent-native MCP+Skills | landed | measured tonight or earlier | GREEN |

## Historical note (2026-07-05 → 2026-07-06)
The original 2026-07-05 sequencing had B-1 → ME-1 → ME-2 all land within
a ~1-week window. That schedule was met (ME-1/C-13/ME-2 all GREEN by
2026-07-06). ME-3 (cluster) and ME-4 (full-physics) remain the two RED
rows; both need external inputs (cluster access, M-F completion) rather
than more implementation velocity here. Track them via docs/PLAN.md, not
via a paper-first schedule.
