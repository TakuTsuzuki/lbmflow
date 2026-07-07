# Claims Ledger — capability status snapshot

**What this is (2026-07-07 clarification, overrides the 2026-07-05 paper-first
framing):** a working table that maps each notable product claim to its
implementing item and the measurement that verifies it. It is a status
snapshot, NOT a release gate — the paper describes what is measured today and
gets updated when the measurements do. The product goal is to become the best
conceivable LBM simulator; this ledger just tracks where we are on that.

Authoritative release-facing limitations and unsupported combinations are in
[`docs/LIMITATIONS.md`](../LIMITATIONS.md); this ledger is only the measurement
status snapshot behind paper claims.

Update rules:
- When a measurement changes, update the Status cell here in the same commit.
- When a row becomes GREEN, edit the paper to reflect it (or add coverage).
- Do NOT hold implementation to a paper claim ahead of the measurement —
  implementation converges to physics/spec, the paper follows.

| Claim (paper, present tense) | Implementing item | Hard gate | Status |
|---|---|---|---|
| 3D GPU acceleration (D3Q19) | ME-1 (B-1 → WGSL 3D + C-13 explicit path) | T14-3D ≤1e-5 + ≥1,500 MLUPS single GPU | **GREEN** — T14-3D GREEN (32³ TGV3D u 2.8e-6, 24³ cavity3D u 1.7e-6); MLUPS GREEN 2026-07-06 quiet-window A/B/A: 192³ 2791-2813, 128³ 2778-2880 MLUPS on unmodified main (the earlier 1353 was a loaded-window artifact; the step_periodic kernel follow-up measured SLOWER and was rejected — TESTING_NOTES 2026-07-06) |
| explicit backend:"gpu" runs on GPU | C-13 (bundled into ME-1) | 2D f32 scenario gpu request honored end-to-end, guard flips from error to run | **GREEN** — landed 2026-07-06 (commit 1a14d90); scope is the scenario path exposed today: 2D f32 GPU when built with `gpu`; README capability matrix keeps 3D scenario GPU, f64 GPU, MPI selector, and unsupported feature combinations out of this claim |
| FP16 storage, ×2 grid capacity | ME-2 (C-12) | T16 capacity/throughput-mode bands frozen + ≥1.5× MLUPS @2048² | **GREEN** — 2026-07-06: T16 implemented and gated for capacity/throughput mode only; D2Q9 frozen bands pass (TGV transient 2e-1 / measured 1.401e-1; cavity steady 5e-3 / measured 2.579e-3; steady-vs-transient dichotomy in PHYSICS.md), MLUPS ~2.0× @2048² (interleaved A/B ×2), D3Q19 f16 >5 GLUPS; capacity ×2 inherent; not a validation-grade long-transient reference mode |
| Multi-node scaling | ME-3 (cluster campaign) | 64-rank weak ≥80% measured | RED — true multi-node weak scaling awaits cluster measurement; current MPI evidence is multi-rank single-node functional coverage / n≤4 weak scaling, and MPI is not exposed through scenario JSON or MCP |
| Full-physics stirred workload | ME-4 (MF-ζ) | degradation ratio vs single-phase published | RED — T17 is mixed: W0/W-ROT/W-GRAV/W-LES and current W-PART scope are landed/validated or characterized as listed in VALIDATION.md; W-VOF, W-BCTOP, W-SCAL, W-REACT, W-BUB, and full coupled W-COUP/W-IO remain pending, so no full-physics degradation ratio is measured today |
| 2D GPU GLUPS / CPU MLUPS / T13 bit-exact / weak-scaling n≤4 / wasm bit-identity / agent-native MCP+Skills | landed | measured tonight or earlier | GREEN |

## Historical note (2026-07-05 → 2026-07-06)
The original 2026-07-05 sequencing had B-1 → ME-1 → ME-2 all land within
a ~1-week window. That schedule was met (ME-1/C-13/ME-2 all GREEN by
2026-07-06). ME-3 (cluster) and ME-4 (full-physics) remain the two RED
rows. Treat them as unmeasured claims here; use VALIDATION.md and
LIMITATIONS.md for the current supported scope.
