# LBMFlow Implementation Plan

> **2026-07-05 revision**: per user directive, 3D, supercomputer-scale, and GPU are mandatory.
> Requirements framework: [COMPETITIVE_SPEC.md](COMPETITIVE_SPEC.md) (R1–R5 + T13–T16 equivalence framework);
> strategic positioning: [paper/LBMFlow-whitepaper.md](paper/LBMFlow-whitepaper.md).
> Architecture: [ARCHITECTURE_V2.md](ARCHITECTURE_V2.md) (dimension × lattice × precision × backend × partitioning).
> Milestones M-A–M-D landed; M-E mostly GREEN (ME-3 cluster + ME-4 stirred RED); M-F in flight.

**Goal**: A commercial-grade lattice Boltzmann method (LBM) fluid simulator.
It can explicitly control the accuracy-vs-speed tradeoff, supports multiphase flow, and provides both a
GUI mode that even beginners can use without getting lost and an Agent mode operable from agents.

**Technology stack**: Rust (core engine / CLI / MCP) + TypeScript (GUI, Vite) + WASM (browser execution)

**Team structure**: Fable = PM / architect / integration validation. Implementation is delegated to Claude subagents and Codex
at an executable granularity. **The boundary-condition test suite is created adversarially by Codex**, and
the engine side is fixed until it passes them all (separating test authors from implementers to detect specification bugs).

---

## Architecture

```
LBMFlow/
├── Cargo.toml              # workspace
├── crates/
│   ├── lbm-core/           # V2 core engine (pure Rust library, no I/O)
│   │   ├── src/
│   │   │   ├── lattice.rs      # D2Q9 / D3Q19 / D3Q27 velocity sets and derived tables
│   │   │   ├── real.rs         # f32/f64 generic (Real trait)
│   │   │   ├── fields.rs       # q-major SoA deviation storage with halo padding
│   │   │   ├── params.rs       # BGK / TRT / cumulant collision and BC parameters
│   │   │   ├── backend.rs      # backend trait + CPU scalar reference backend
│   │   │   ├── backend_simd.rs # fused CPU SIMD/band-parallel backend
│   │   │   ├── gpu/            # wgpu backend (feature "gpu")
│   │   │   ├── halo.rs         # local/in-process halo exchange
│   │   │   ├── dist.rs         # MPI exchange/solver (feature "mpi")
│   │   │   ├── solver.rs       # V2 solver orchestrator
│   │   │   └── compat/         # supported legacy 2D facade over V2
│   │   └── tests/          # validation tests (including the codex-authored adversarial suite)
│   ├── lbm-cli/            # JSON scenario execution CLI + MCP stdio server
│   └── lbm-wasm/           # wasm-bindgen bindings
├── web/                    # TypeScript GUI (Vite)
└── docs/
    ├── PLAN.md             # this file
    ├── VALIDATION.md       # validation-test specification matrix (commissioning spec for codex)
    └── PHYSICS.md          # rationale for adopted physics models and formulas (appended as needed)
```

### Core design highlights

- **Lattice**: V2 is generic over the `Lattice` trait with D2Q9, D3Q19, and D3Q27 implementations.
  Population storage is q-major SoA deviation form over halo-padded local boxes:
  `f[q * n_padded + cell]`, with `cell = z * (pnx * pny) + y * pnx + x`.
- **Streaming / step order**: the common solver orchestrates
  `collide` → halo exchange → `stream` → `swap` → open-boundary BCs → moments update.
  `CpuSimd` fuses collide/stream internally while preserving the backend contract.
- **Collision operator**: BGK, TRT (with magic Λ=3/16, wall position is exact, recommended default),
  and the landed central-moment/cumulant branch (`CollisionKind::Cumulant`).
- **Precision**: f32/f64 are selected via `Real`; GPU execution is f32-only today.
- **Backends / distribution**: CPU scalar and fused CPU SIMD backends are available in the core;
  the wgpu backend is behind feature `gpu`, and MPI halo exchange / distributed runs are behind feature `mpi`.
- **Body force**: Guo forcing (2nd-order accurate, u includes F/2 correction). Used in Shan-Chen too.
- **Walls**: half-way bounce-back (stationary walls, moving walls). Edge-specified walls are realized as a 1-cell solid
  rim → the corner special-casing with Zou-He becomes unnecessary.
- **Open boundaries**: Zou-He velocity/pressure, zero-gradient outflow, and convective outflow are face-normal
  parameterized. Open-face kernels currently support D2Q9 and D3Q19; D3Q27 open-face support remains explicit follow-on work.
- **Force measurement**: momentum-exchange method (required for the cylinder Cd/St benchmark).

### Boundary-condition combination rules (specification)

- Periodic is selected per axis and cannot be combined with an open BC on that axis.
- Open faces may lie on at most one axis. Perpendicular non-periodic faces must be fully covered by a wall rim,
  or covered by validated face patches.
- Non-positive viscosity (`tau = 3*nu + 0.5 <= 0.5`) is a construction-time error.
  Prescribed speed above 0.3 lattice units is a construction-time error.

---

## Phase plan

(2026-07-07 historical record — describes the state at the original phase milestones; see ARCHITECTURE_V2.md for the current design.)

| Phase | Content | Completion criteria |
|---|---|---|
| 0 | Foundation: git / workspace / PLAN / VALIDATION / CLAUDE.md | Documents committed |
| 1 | lbm-core vertical slice (D2Q9, BGK/TRT, Guo force, BB/moving wall/periodic/Zou-He/Outflow, force measurement) + smoke tests | TGV convergence order ≈2, Poiseuille(TRT) exact, Couette exact, conservation-law tests green |
| 2 | **Have codex implement the adversarial test suite from VALIDATION.md** → fix the engine until all pass. Includes the cavity (Ghia) and cylinder (Cd/St) benchmarks | `cargo test --release` all green |
| 3 | Accuracy × speed: f32 validation, MLUPS benchmarks, thread scaling, mode-selection guide, **deviation-storage scheme** (keep f−w to raise the effective precision of f32. BB is linear so it is invariant in deviation space, Zou-He needs the constant term folded in, ρ = 1+Σdev) | The measured tradeoff table appears in docs. The T6 error for f32 improves |
| 4 | Multiphase flow: Shan-Chen (single-component multiphase + two-component) + codex validation (Laplace law, contact angle, Rayleigh-Taylor) | Multiphase tests all green |
| 5 | GUI: wasm-pack + Vite + TS. Preset-driven, obstacle painting, real-time visualization, Japanese UI | 5 presets run with one click |
| 6 | Agent mode: lbm-cli (JSON scenario → structured output) + MCP server | An agent can run a scenario and retrieve results |
| 7 | Comprehensive review → formulate the next plan (3D/GPU/LES/cumulant, etc.) → continuous improvement | Review document + next plan |

### Validation-driven development protocol

1. State acceptance criteria numerically in the specification (VALIDATION.md).
2. codex writes the tests (in principle written from the specification without looking at the implementation).
3. If the engine fails: (a) engine bug → fix, (b) specification physics is wrong → experiment and fix the specification
   (record the reason for the fix in PHYSICS.md), (c) test bug → send back to codex.
4. Move to the next phase when all green. Run `git commit` at the end of each phase.

## Current queue: improvement phase (R-Phase) and M-F (established 2026-07-05)

2 source documents: [SOLVER_IMPROVEMENT_SPEC.md](SOLVER_IMPROVEMENT_SPEC.md) (41 items from the
concurrent review, all claims validated by experiments E1–E10, merged into main) and
[REQ_STIRRED_REACTOR.md](REQ_STIRRED_REACTOR.md) (M-F requirement rev.1b, 48 codex adversarial-review
items reflected, acceptance is VALIDATION.md **T17**).

### R-Phase (follows the execution order in §4 of the improvement spec)

- **R-Phase 1 (in progress 2026-07-05–)**: A-2–A-10 = entry guards, correctness.
  worktree `r-phase1`, delegated to Opus. Common DoD = existing tests green without modification, legal-configuration bits invariant.
  D-6/D-7 (document consistency) already applied under direct PM control (COMPETITIVE_SPEC revision history, VALIDATION T13/T14 sections).
- **R-Phase 2 (after R-1 lands, ~1.5 weeks)**: B-1 (Fields generalization + GpuSolver integration) → B-2
  (sync-point contract) → concurrent B-3 / B-5–B-8, C-9–C-11, D-1–D-5.
  **Additional requirement from M-F**: the B-1 design must accommodate multiple distribution sets (phase-field g, scalar h), per-cell
  material-property fields (generalization of B-6), and Lagrangian buffers (IBM markers/particles/point bubbles)
  — this is a structural premise of M-F and a foundation shared with M-E (FP16/multi-GPU).
- **R-Phase 3 (can run concurrently with M-E)**: C-1 (localizing the MPI setup = resolving 10⁹-lattice OOM),
  C-2 (communication overlap. Premise: the E8 probe double-counting shell fix), C-4–C-8,
  C-12–C-16, D-8–D-10.
- **M-E is premised on completion of B-1/B-2/C-9/C-12/C-13/D-9** (dependency relations in §4 of the spec).

### External review response (2026-07-07) — claim hygiene + product-path gaps

An external review (static cross-check of docs vs implementation) was triaged by the PM.
Verdict: core findings confirmed; the FSI-overclaim and "supercomputer-class claim" points
were rejected (the repo never claims FSI; multi-node is already RED in the claims table and
"supercomputer-scale" appears only as a roadmap mandate in this file).

**P0 — dispatched 2026-07-07 as 7 parallel codex orders (branches `cx/*`):**
ci.yml nightly schedule + dispatch + artifacts (`cx/ci-schedule`); VALIDATION T16/T17
status sync + D3Q27 open-face note (`cx/validation-sync`); REQ particle status tags
(`cx/req-sync`); README capability matrix (`cx/readme-cap`); docs/LIMITATIONS.md
(`cx/limitations`); SOLVER_IMPROVEMENT_SPEC resolved/open audit (`cx/sis-audit`);
particles.rs silent Re>800 clamp → explicit validity-domain error (`cx/re-clamp`).

**P1 — queued (new items; overlapping ones fold into existing tracks):**

| Item | Content | Folds into |
|---|---|---|
| REV-1 | Scenario schema exposure: `collision: cumulant`, `storage: f16`, `backend: auto\|cpu\|gpu`, lattice selection; unsupported-combination errors; manifest records actual backend/lattice/collision/storage | new (agent-native product path) |
| REV-2 | `lbm verify --tier quick/full/gpu/mpi` + `lbm capabilities` (machine-readable capability/limitation output) | SOLVER_IMPROVEMENT_SPEC open item |
| REV-3 | Gravity force composition device-side (kill per-step host staging in `run()`; prerequisite for GPU VOF/buoyancy performance claims) | B-2 sync-point contract territory |
| REV-4 | LES tau_eff upper clipping + diagnostics, y+ output, wall treatment | MF-β scope (REQ W-LES already requires it) |
| REV-5 | D3Q27 open-face kernels + dedicated gates (Zou-He moment exactness, duct in/out, CPU/GPU equivalence, MPI seam) | MF-α follow-on |
| REV-6 | Cumulant D3Q19 shear-rate offset: classify as empirical calibration; holdout validation set (advected TGV at multiple frame velocities, rotating shear, channel, Ma sweep, D3Q19 vs D3Q27) | physics-discipline follow-up, PHYSICS.md |
| REV-7 | Distributed checkpoint/restart (per-rank, RNG/particle/statistics state, restart match test) | R-Phase 3 / C-track |

**P2:** UQ/sensitivity harness (parameter sweeps with CI-band QOIs, calibration/holdout
separation); absolute-physics GPU gates beyond equivalence (VALIDATION already flags the
blind spot); geometry ingestion (STL/voxelization) — deferred, not currently claimed.

### Performance roadmap — ALL FOUR gaps are committed implementation items (user directive 2026-07-05: "implement all of these; at minimum put them on the roadmap")

The four gaps identified by the sales-paper analysis (vs FluidX3D / waLBerla / OpenLB /
M-Star) are hereby **committed**, in this order:

| # | Item | Acceptance line | Depends on | Wave |
|---|---|---|---|---|
| ME-1 | **3D GPU** (D3Q19 WGSL kernels + BC passes) | R2: single-GPU D3Q19 f32 **≥1,500 MLUPS** (expect multi-GLUPS on M5 Max); T14 extended to 3D | R-Phase 2 B-1 (orchestrator unification) | immediately after R-2 |
| ME-2 | **FP16 storage** (C-12, shader-f16; compute stays f32) | T16 degradation bands frozen (TGV/cavity); **≥1.5×** MLUPS vs f32 at 2048²; grid capacity ×2 | B-1, ME-1 (shares kernel plumbing) | with ME-1 wave |
| ME-3 | **Multi-node / cluster measurement** (R3) | 64-rank weak scaling **≥80%** on a real cluster; the 8-item measurement list in MPI_GUIDE §cluster | C-1 (MPI setup localization) + cluster access — **AWS hpc7g×8 (~¥13k, CLUSTER_OPTIONS.md recommended) or Fugaku trial**; committed on the roadmap, instance spend gets a one-line user confirm at execution | R-Phase 3 window |
| ME-4 | **Full-physics benchmark** (stirred-tank workload: two-phase + particles + scalar + LES) | performance-degradation ratio vs single-phase measured and published (the M-Star-comparable number); runs as part of MF-ζ acceptance | M-F tracks | MF-ζ |

The paper (docs/paper/LBMFlow-whitepaper.md) describes what is measured today;
edit the paper as the implementation lands, not the other way around. The
acceptance lines above are the technical targets the roadmap converges on —
they are not paper claims to true up at release. See docs/paper/claims-ledger.md
for the current measurement-status snapshot.

### M-F: rotating boundary, high-density-ratio two-phase, LES-coupled 3D (REQ-M-F-STR rev.1b)

Confirmed design decisions (decided by the project owner): **scope all-at-once** (implement subsystems simultaneously
without staged splitting) / **fidelity defaults** (IBM-inertial, resolved-phasefield, active scalar,
two-way particles, uniform lattice, f64 near the interface + f32 in the bulk); low-cost approximations (MRF, point-bubble,
one-way, AMR, aggressive f32) are add-on extensions behind the same trait / physics-conflicting modes are
mutually exclusive at runtime via configuration validation (extending A-4's `GlobalSpec::validate`).

Implementation tracks (commissioned in parallel after R-Phase 2 lands. Conventional team structure of worktree separation, implementation by Opus/Sonnet,
**validation tests created adversarially by codex from REQ/T17**):

| Track | Content | Primary FR | Validation | Dependency |
|---|---|---|---|---|
| MF-α | D3Q27 lattice + central-moment/cumulant collision | FR-CORE-01/02 | moment isotropy, TGV3D order, Galilean-invariance band | R-2 |
| MF-β | LES (WALE default) + non-Newtonian μ(γ̇) + stress-field evaluation (convention FR-STRESS-01 fixed) | FR-LES-*, FR-STRESS-* | VR-STR-03, channel Re_τ=180 vs DNS | R-2 (B-6) |
| MF-γ | conservative Allen-Cahn high-density-ratio two-phase (10³) + well-balanced gravity + sparger/degassing BC | FR-VOF-*, FR-BC-* | VR-STR-02/06, Laplace, parasitic current Ca<10⁻³, single-bubble Grace | R-2 (B-1 multiple distributions) |
| MF-δ | IBM-inertial rotating boundary + torque/Np measurement | FR-ROT-* | VR-STR-01, Taylor-Couette torque, IBM sphere drag vs T15 reference | R-2 (B-1 Lagrangian) |
| MF-ε | scalar ADE (active feedback) + Lagrangian particles (shear-exposure recording) | FR-LES-04, FR-PART-* | VR-STR-04, Taylor-Aris, settling terminal velocity vs SN | R-2 (B-1) |
| MF-ζ | coupling integration (§5 dataflow, dt constraints) + configuration exclusivity + I/O/statistics/GUI 3D display + acceptance run | FR-COUP-*, FR-INIT, FR-IO-* | VR-STR-05/07 + coupled system of 01/02 | MF-α–ε |

Status of remaining spec details (late night 2026-07-05): **active scalar feedback formula** = research complete
(docs/proposals/active-scalar-feedback.md. **1 derivation required before implementation**: consistency of the
Liu et al. Marangoni-form W²↔(κ,β) coefficients — already stated in REQ §3) /
**f64/f32 interface band width** = frozen by experiment at MF-γ implementation time (characterize→freeze) /
**trait boundary API** = confirmed together with the B-1 design of R-Phase 2 (appended to ARCHITECTURE_V2) /
**REQ 2nd-round codex validation** = complete, all 11 items adopted → **rev.2 applied**
(new VR-STR-RELAX, delivery-scope clarification, variable-σ convention, and others.
Findings original: docs/proposals/req-round2-findings.md).

Effort sense (assuming parallel agents): R-2 ~1.5 weeks → MF-α–ε parallel ~1-2 weeks → MF-ζ integration ~1 week.
The 1e9-lattice class is **cluster-only** per the memory budget table (REQ §7) (the single-machine development line is ≤256³) —
measurements are consolidated into the cluster plan (CLUSTER_OPTIONS.md, awaiting user decision).

### D-track: dispersed-phase deposition-design tool (established 2026-07-06)

AI-agent-native forward/inverse design of deposited number-density fields
n(x,y) from a withdraw/eject/agitate/settle protocol. Spec + phasing (P0–P4)
frozen in [DISPERSED_DEPOSITION.md](DISPERSED_DEPOSITION.md); acceptance =
VALIDATION.md **T18**. Status: **P0–P2 done** (2026-07-06). P2 landed CR-1
interior volume source/sink, CR-2 per-cell masked face patches (incl. the
frozen zero-velocity-lid semantics and global-coords seam handling), and CR-3
deposition-aware particles (`step_depositing`), each gated by a
codex-adversarial T18 suite (T18.1/.2/.3 green; 2 impl bugs caught and fixed
at first measurement — see PHYSICS.md "T18 first-measurement
reconciliation"); the example now runs on the real core capabilities and
reproduces the P1.1 bands (gentle CV=1.138 ∈ [1.05, 1.30], deposition parity
0.4%). Next: P3 (VOF-on-LBM free surface) is evidence-gated, not speculative;
P4 = inverse solver (discrete-recipe comparison first, then CMA-ES/BO +
surrogate) — T18.4 trend anchors and T18.5 recovery are the gates.

**Fine-grained scheduling (rev.3)**: the authoritative dependency DAG is
**REQ_STIRRED_REACTOR.md §11** (W-items; MF-α〜ζ above are the delegation bundles,
each row maps to its track). Execution shape: after W0 (=MF-α core basis), wave 1
runs **6-way parallel** {W-EXT, W-UNIT, W-STRESS, W-ROT, W-GRAV, W-SCAL}; then
{W-LES, W-VOF, W-PART, W-REACT}; then {W-BCTOP, W-BUB, active feedback};
W-COUP / W-IO / W-VAL run cross-cutting. **Critical paths to staff first**:
`W0→W-GRAV→W-VOF→W-BCTOP` (interface chain) and `W0→W-STRESS→W-LES→W-PART`
(stress/exposure chain). W-EXT is co-designed with R-Phase 2 B-1 (one trait design,
two consumers). Validation (W-VAL) stays codex-adversarial and implementation-separated.

**Language policy (2026-07-05 user directive)**: all artifacts English from now on
(code, docs, commit messages, UI/CLI strings). Legacy Japanese content is being
translated by a dedicated spawned session; this file will be fully translated there.

## Progress notes

(2026-07-07 historical record — dated notes preserve what was true or planned at each milestone; see ARCHITECTURE_V2.md for the current design.)

- 2026-07-05 late night: **Integrated the results of the concurrent review session into main**. Improvement spec v1 +
  experiment crate merged (E2/E7 reproduce numerically matching values on main after renaming). 4 PM decisions confirmed:
  (a) spec merge = done by PM (b) R-Phase 2 is right after R-1, a common premise of M-E/M-F
  (c) D-6 = applied under direct PM control (see COMPETITIVE_SPEC revision history) (d) codex D-8 commissioning is
  after R-1 lands. **R-Phase 1 commissioned** (worktree r-phase1, Opus, A-2–A-10).
  **M-F requirement rev.1b confirmed** (neutral title, memory budget table, T17 wiring) and the implementation-track plan
  (table above). codex #7 (T15.5 3D cavity) in progress.
- 2026-07-05 evening: **Officially judged R1/R2/R3 all achieved** (REVIEW_2026-07-05_2.md.
  Note: the acceptance band takes the 2026-07-05 D-6 revision as authoritative: sphere drag ±10%, D_h normalization, the weak-scaling
  85% line is n≤4 local, 3D cavity T15.5 is being additionally implemented in codex #7).
  M-D MPI complete (T13-MPI bit-identical in all cases, weak scaling 97-99% n≤4, MPI_GUIDE complete).
  CpuSimd fused backend (2D new record 1,183 MLUPS, 3D 2-3x, equivalence ~6e-14).
  Workspace 205 tests green. In progress: V1 retirement (v1-retirement) + 3 research items
  (public benchmark comparison / 3D cavity reference / cluster options). CI workflow prepared.
- 2026-07-05 evening: **M-C 3D physics complete (R1 achieved)**. D3Q19 Zou-He face boundary (5 unknowns + tangential correction,
  2D degeneration 8.9e-16), duct exact series 2.3e-4, sphere drag with hydrodynamic-pair (r+½, Re(D+1)/D)
  normalization gives +7.1%/+0.6%/+2.3% (all pass including the heavyweight D=24), TGV3D order 1.910.
  Scenario/CLI nz support (VTK 3D / cross-section PNG). **Wgpu backend integration (R2)**: with push-type fusion,
  strict operator-order match with CPU, T14 6 configurations ≤1e-5, 5.9–11.4 GLUPS. Workspace 184 green.
  In progress: M-D MPI distribution (arm64 Open MPI 5.0.9 already source-built into ~/.local).
- 2026-07-05 afternoon: **M-B core V2 complete and integrated** (physically equivalent to V1, T13 partition invariance is bit-identical-class up to 8 adversarial attacks + 3D 2×2×2, D3Q19 smoke works, 8 V1 implicit-spec items frozen as tests).
  **Phase 9 complete**: CPU fused kernel 3.2–7x (f32 peak 1,124 MLUPS) + GPU measured
  6,975–12,152 MLUPS. MCP asynchronous job API (R4). Survived all of codex #6's adversarial T13.
  In progress: 2-way parallel of the Wgpu backend real implementation (m-b-wgpu) and 3D physics M-C (m-c-3d).
- 2026-07-05 midday: **Phase 8 complete** (T12 RT γ ratio 1.118, T11c contact-angle full range, T9b convective outflow
  improved 16x in the reflection case, 67 tests green). **GPU measurement complete**: 6,975–12,152
  MLUPS on M5 Max Metal (16–42x vs CPU, validation L∞ 7e-6, GPU_EVALUATION.md, exceeds the R2 target by 4-8x).
  Scenario/CLI gains convectiveOutflow, wallRho, VTK, gallery. GUI gains scenario export
  (→lbm run E2E confirmed), divergence guard, MLUPS. SoA/SIMD is WIP (phase9-perf, resumes after quota
  recovery). Next: parallel commissioning of M-B core V2.

- 2026-07-05: **Phase 4a/5/6 complete (three-mode unification achieved)**. codex #4 multiphase validation
  all green (coexistence densities, EOS pressure equilibrium, Laplace, contact-angle regression 133/160/164°, f32 hardening 1e-5).
  GUI: cavity/Karman vortex/two-phase droplet run on the real WASM engine (~600 steps/s).
  Agent mode: lbm CLI (run/validate/presets/schema, manifest + PNG/CSV) +
  MCP server (4 tools including run_scenario). Workspace 56 tests all green.
- 2026-07-04: Project started. Phase 0 begun.
- 2026-07-04: Phase 1 complete. Implemented lbm-core (D2Q9, BGK/TRT, Guo force, half-way BB/moving wall,
  Zou-He, Outflow, force measurement, f32/f64, rayon + small-lattice serial fallback).
  21 smoke tests green: TGV 2nd-order convergence (1.91/1.98), Poiseuille TRT exact
  (<1e-10), Couette exact, conservation laws ~1e-13. Recorded 4 experimental findings in PHYSICS.md.
- 2026-07-05: **Phase 4a implementation complete** (validation awaiting codex #4). Per-cell force-field API +
  Shan-Chen SCMP. Measured: density ratio 15.8, pressure equilibrium 8.5e-6, spurious velocity 1.3e-3,
  Laplace R²=0.9999. **Scope reorganization**: MCMP+RT (Phase 4b) moved to after the first review.
  Do GUI (Phase 5) / Agent mode (Phase 6) first and complete the three-mode unification first.
- 2026-07-05: **Phase 3 complete**. Introducing the deviation-storage scheme (f−w) brought f32 to validation grade
  (momentum error improved 4800x, on par with f64 in TGV). MLUPS measured: peak 381 (f32/1024²/18T),
  single 35, TRT is the same speed as BGK. Tradeoff guide in PERFORMANCE.md.
  Maintained all 49 tests green.
- 2026-07-05: **Phase 2 complete**. 3 rounds of codex adversarial tests (triaged 9 findings total:
  2 engine bugs = Zou-He pressure sign, rim-corner anisotropy / 5 specification bugs / 1 reference-data
  typo / 1 f32 characteristic). Default 49 tests, full 53 tests all green.
  Passed Ghia cavity Re=100/400/1000, Schäfer-Turek cylinder 2D-1/2D-2,
  exact equivariance (machine precision 4e-16), Zou-He 4 directions, conservation laws.
  The GUI shell (Vite+TS+mock) was also completed ahead (Phase 5a).
