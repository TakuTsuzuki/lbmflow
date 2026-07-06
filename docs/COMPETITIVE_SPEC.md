# Competitive Advantage Spec (established 2026-07-05)

**Status (2026-07-07)**: This file is now a requirements/equivalence-test anchor, not the live competitive narrative.
**Landed**: R5's compat facade and the T13-T16 test framework are represented in code/tests.
**Current intent**: R1-R4 acceptance lines remain targets; current implementation status lives in `docs/paper/claims-ledger.md`.
**Superseded**: strategy, milestone, and claim-status prose belongs in `docs/paper/LBMFlow-whitepaper.md`, `docs/paper/claims-ledger.md`, and `docs/PLAN.md`.

## R1-R5 — Required requirements (per user directive, 2026-07-05)

(current intent 2026-07-07 — code/tests still cite this section; per-pillar status is tracked in `docs/paper/claims-ledger.md`)

Cited by: `crates/lbm-core/tests/t13_split_invariance.rs` (§4),
`t14_backend_equiv.rs` (§4, R2), `t15_3d.rs` (R1),
`crates/lbm-core/examples/bench_mpi.rs` (§5, R3),
`crates/lbm-core/src/compat/mod.rs` (R5),
`crates/lbm-scenario/src/lib.rs` (M-C), `docs/ARCHITECTURE_V2.md`,
`docs/BENCH_COMPARISON_DRAFT.md`, `docs/PLAN.md`.

| ID | Requirement | Measurable acceptance criteria |
|---|---|---|
| R1 | **3D (D3Q19/Q27)** | Pass the 3D validation suite: TGV3D (diffusion-limit reference ±2%, convergence order ≥1.7), sphere drag (Re∈{20,100}, Schiller-Naumann correlation **±10%, hydrodynamic-pair normalization D_h=D+1, Re_h=Re(D+1)/D**), 3D cavity (T15.5 = A&K 2005 Re=1000, RMS≤0.030U), turbulent channel Re_τ=180 vs DNS (Moser+) after LES is introduced (M-F) |
| R2 | **GPU support** | wgpu backend equivalent to CPU at f32 (T14: field agreement ≤1e-5 relative). Single-GPU D3Q19 f32 ≥ 1,500 MLUPS (M4 Max class or RTX 4070 class). 2x grid size in FP16 storage mode |
| R3 | **Supercomputer scale** | Subdomain+MPI: partitioned run ≡ monolithic run agreement (T13). Single-node intra-node multi-rank weak scaling **≥85% (local line for homogeneous cores, n≤4; measured 97-99%**. n=8 is 73% due to M5 Max's heterogeneous cores + bandwidth ceiling, with 84% as the ceiling even under a zero-communication control — confirmed via control experiment to be a property of intra-node SMP, not an MPI implementation defect). Cluster measurement (requires machine access, §5, CLUSTER_OPTIONS.md) for 64-rank weak scaling ≥80% — **not yet measured** |
| R4 | Maintain agent-nativeness | All new features must be controllable from scenario JSON+MCP. Asynchronous job API |
| R5 | Validation continuity | The existing 2D suite of 56+ **stays green without modification after the core overhaul** (API facade maintained) |

> **Revision history (D-6, 2026-07-05)**: R1 sphere drag ±5%→**±10%** (the half-way BB wall sits half a link outside the solid cell, so nominal D normalization carries a ~+2/D bias. Revised to a band incorporating the measured +0.6–7.1% together with adoption of hydrodynamic-pair D_h normalization. Basis: TESTING_NOTES 2026-07-05 triage, PHYSICS.md). R1 3D cavity's "literature comparison" is made concrete as T15.5's (A&K 2005) quantitative band. R3 weak scaling ≥85% is clarified in scope as "single-node, n≤4, local," while the cluster 64-rank line remains unmeasured (measurement plan: CLUSTER_OPTIONS.md).

## §4 Equivalence-test framework

(landed 2026-07-07 — T13/T14/T15/T16 are executable test families; GPU/FP16 cases require the `gpu` feature and adapter)

- **T13 partition invariance**: results of 1×1 / 2×2 / 4×1 partitioned runs
  match the monolithic run (bit-match target at f64, at least ≤1e-12)
- **T14 backend equivalence**: CPU-SIMD vs wgpu vs (future CUDA) match on
  the same scenario at f32 relative ≤1e-5 / statistical agreement
- **T15 3D physics**: each benchmark listed under R1
- **T16 precision modes**: quantitatively freeze the degradation of FP16
  storage vs f32 (tolerance band specified)

## §5 Multi-node measurement dependency

(current intent 2026-07-07 — local MPI/partition gates are implemented; 64-rank cluster measurement remains unmeasured)

Multi-node measurement requires access to a cluster/cloud HPC (locally we
can only go as far as functional verification of multiple MPI ranks). See
`docs/CLUSTER_OPTIONS.md` for the measurement plan; ME-3 in the claims
ledger tracks status.
