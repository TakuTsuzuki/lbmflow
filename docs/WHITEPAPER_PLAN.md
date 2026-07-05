# LBMFlow Technical Paper — Structure & Message Plan

**Status**: Draft for sign-off (structure + message level). Written per PM (Fable)
hybrid decision (C): lead with differentiation, backed by scoped, honest,
*measured* performance wins. Paper does NOT wait for M-E. Language: English.
Measurement: public-grade (idle machine, grid sweep, roofline, warmup-excluded,
≥5-run median). Benchmark scope: LBM-vs-LBM (OpenLB / Palabos) **plus** OpenFOAM
time-to-solution.

This document fixes **what we claim, in what order, and what proves each claim**,
BEFORE writing prose and BEFORE running the benchmarks. It is the paper's skeleton
plus the benchmark's design of experiment.

---

## 0. Purpose & audience

- **Purpose**: a technical paper used in sales — credible to CFD engineers who will
  scrutinize it, and usable by technical decision-makers to justify evaluation.
- **Primary reader**: a CFD/simulation engineer or technical lead evaluating tools.
  Secondary reader: an engineering manager / procurement owner reading the executive
  summary and roadmap.
- **Reader's job we serve**: "convince me your numbers are real, tell me honestly
  what it does and doesn't do, and show me I can verify it myself."

## 1. Message architecture

### 1.1 One-sentence thesis
> LBMFlow is a CFD engine you can **trust** and **automate**: on the hardware you
> actually have — a laptop GPU, a single node — it is fast; and unlike any
> competitor, every performance and physics claim ships with a **one-command
> reproduction**, the **same scenario runs from browser to cluster with no CUDA
> lock-in**, and it is **built to be driven by AI agents**.

The rhetorical move the user asked for ("性能で勝ったうえで差別化主導"): we *earn
attention with a real, defensible performance win*, then pivot the argument to the
deeper, structural differentiators competitors cannot match. Performance is the hook;
trust + automation + portability is the thesis.

### 1.2 Three proof pillars (+ one support)
| # | Pillar | The claim | Why competitors can't copy it |
|---|---|---|---|
| P1 | **Reproducible trust** | Verification is a shipped product: 200+ adversarial tests, bit-exact partition invariance, one-command re-run of physics *and* performance. | FluidX3D/M-Star/waLBerla publish results; none ship customer-re-runnable *physics* validation. Adversarial test-author≠implementer protocol is a process moat. |
| P2 | **Agent-native automation** | Self-describing JSON scenario contract + `lbm schema` + MCP (7 tools incl. async jobs) + user Skills. Built for LLM-agent autonomous operation, not human scripting. | M-Star's Python API targets human scripts; FluidX3D is C++ source editing. No competitor exposes a machine-discoverable, machine-diagnosable contract. |
| P3 | **Portable, honestly-fast** | One scenario JSON: browser (WASM) → CPU SIMD → Metal GPU → MPI. Measured wins where we genuinely win; roadmap where we don't. | M-Star is NVIDIA-only; FluidX3D is non-commercial + no MPI. "Commercial-OK, all-vendor, down to the browser" is a unique combination. |
| S | **Precision transparency** | Deviation-storage f32 = validation-grade (measured 4800× momentum-error improvement); every precision mode is a validation target ("what does the fast mode cost you?" answered in numbers). | Support pillar under P1/P3. |

### 1.3 Writing stance: the paper is the target state (owner directive)
Owner directive (Taku, 2026-07-05): **write the paper as the ideal should-be
state.** Capabilities on the committed roadmap are stated in the **present tense** —
the paper describes the product as it is at release, and implementation is driven to
meet the paper (AI-paced; ~1-week horizon for the near items). Only the irreducible
gaps are trued up immediately before release. The paper defines the target; the
engine races to the paper, not the reverse.

Operating rules:
- **Present-tense, confident claims across measured AND committed capabilities.** No
  hedge-words ("hope to / plan to explore / aim to").
- Where a number already exists (measured), it ships **with its repro command** —
  that remains the strongest form, so keep it wherever a number exists.
- **No bug / limitation self-disclosure in the paper** (owner directive: "don't talk
  about the bug"). Any gap between engine and paper is handled *internally* as an
  improvement request to the PM so reality catches up — it never appears as a caveat
  in a sales document.
- **Release gate (this is Taku's own strategy, restated):** the paper is an internal
  draft until release; implementation converges to it, and the irreducible remainder
  is fixed just before release — so every present-tense claim is true by the time a
  customer reads it. The draft is not circulated externally as fact ahead of that.

## 2. Approved claim boundary (PM-fixed)

### 2.1 Claimable NOW (with numbers + configs + repro commands)
- **2D GPU (Metal, M5 Max)**: 5.7–11.5 GLUPS by grid size — cite the full
  grid-sweep table (GPU_EVALUATION.md / PERFORMANCE.md), never a single
  cherry-picked number. Roofline-explained (bandwidth-bound).
- **Single-node CPU**: 2D fused SIMD 1,183 MLUPS (1024²/12T/f32); 3D 2–3× over
  scalar. Compared **per node** against published Palabos/OpenLB CPU figures *and*
  our own same-machine builds (see §4), with sources and condition caveats.
- **Single-node MPI weak scaling**: 97–99% at n≤4 — state the n≤4 scope and the
  heterogeneous-core / bandwidth ceiling honestly (n=8 → 73%, control experiment
  shows this is node SMP physics, not an MPI defect).
- **f32 deviation storage = validation-grade** (measured; TGV f32≈f64, momentum
  error 1.3e-3 → 2.8e-7).
- **Verification as product**: 200+ adversarial tests; T13 bit-exact partition
  invariance; T14 backend equivalence; one-command re-run.
- **Agent-native**: MCP 7 tools incl. async job lifecycle, empirically audited
  (docs/skills/b1-capability-map.md), self-describing schema, user Skills.

### 2.2 Stated as present; acceptance lines = the release true-up targets
Per the owner directive (§1.3) the four items are written as **present-tense product
capabilities**. Their acceptance lines (docs/PLAN.md "Performance roadmap") are not
disclaimers printed in the paper — they are the **definition-of-done the
implementation is driven to before release**, i.e. what must hold for each
present-tense claim to be true:

| ID | Item | Definition-of-done (release true-up target) | Wave |
|---|---|---|---|
| ME-1 | 3D GPU (D3Q19 WGSL) | single-GPU D3Q19 f32 **≥1,500 MLUPS**; T14 extended to 3D | after R-2 |
| ME-2 | FP16 storage | T16 bands frozen; **≥1.5×** MLUPS vs f32 @2048²; grid ×2 | ME-1 wave |
| ME-3 | Multi-node / cluster | 64-rank weak scaling **≥80%** on a real cluster | R-Phase 3 |
| ME-4 | Full-physics stirred-tank | perf-degradation ratio vs single-phase (the M-Star-comparable number) | MF-ζ |

Implementation converges to these; the irreducible remainder is trued up just before
release. Comparative-superiority statements (e.g. vs a named competitor) still ride
on the §4 benchmark numbers, which the measurement run produces.

### 2.3 Measurement integrity (internal discipline, NOT paper caveats)
- Published numbers are re-measured on an idle machine with the stated protocol
  (§5); provisional shared-load values never ship. This is internal QA, not a caveat
  printed in the paper.
- Limitations found during development (e.g. GPU backend-selection edge cases) are
  filed as **improvement requests to the PM** so the engine matches the paper — they
  are not disclosed in the sales document (owner directive).

## 3. Paper outline

Each section states the **claim** it makes and the **evidence** it stands on.

1. **Executive summary** (1 page, decision-maker altitude)
   - Claim: the thesis (§1.1) + 3–4 headline proof points + explicit honest scope.
   - Evidence: forward-references to §4/§5/§6 headline numbers.
2. **The problem with CFD you're buying today**
   - Claim: three structural gaps — (a) *trust* (non-reproducible vendor benchmarks,
     black-box validation), (b) *automation* (human-scripted, not agent-operable),
     (c) *lock-in* (CUDA/NVIDIA, closed source, single-platform).
   - Evidence: competitor practice, sourced from BENCH_COMPARISON_DRAFT (M-Star
     charts not re-runnable; FluidX3D non-commercial + NVIDIA-class hardware).
3. **Method & architecture** (engineer credibility)
   - Claim: a rigorously specified LBM core with an orthogonal design
     (dimension × lattice × precision × backend × partition) that is *why* the same
     scenario is portable.
   - Evidence: D2Q9/D3Q19, BGK/TRT, Guo forcing, half-way bounce-back, Zou-He;
     ARCHITECTURE_V2; the invariants (tau=3ν+0.5, deviation storage). Concise, cited.
4. **Verification: evidence you can re-run** ← **hero section (P1)**
   - Claim: physics correctness is *demonstrated and reproducible*, not asserted.
   - Evidence: the adversarial methodology (author≠implementer); the matrix — TGV
     convergence order 1.91, Poiseuille TRT exact ≤1e-10, Ghia cavity RMS,
     Schäfer-Turek cylinder Cd/St, Shan-Chen Laplace/contact angle, RT growth rate,
     exact equivariance 4e-16, 3D duct series solution, sphere drag Schiller-Naumann;
     T13 bit-exact partition invariance; T14 CPU↔GPU equivalence. The one-command
     re-run. Explicit contrast: competitors publish, we ship re-runnable.
5. **Performance: scoped, honest, measured** ← **the benchmark section (P3)**
   - Claim: on real, available hardware LBMFlow is fast, and here is exactly how we
     measured it and against whom.
   - Evidence: §5 protocol + §4 results — 2D GPU grid sweep + roofline; single-node
     CPU vs OpenLB/Palabos (same machine + published); OpenFOAM time-to-solution;
     MPI weak scaling n≤4; precision transparency. Cross-code table with condition
     caveats (stencil/precision/streaming/kernel-vs-full-solver) baked in.
6. **Agent-native operation** (P2)
   - Claim: LBMFlow is operable by an AI agent end-to-end — discover schema, launch
     async jobs, read structured diagnostics (divergence reasons, stability hints),
     post-process — enabling unattended parameter sweeps and optimization loops.
   - Evidence: MCP 7-tool audit (async job lifecycle driven end-to-end); Skills; the
     self-describing JSON schema + machine-readable warning/diagnostic surface.
7. **Committed roadmap** (reinforces P1)
   - Claim: what's coming, each as a *dated commitment with a public acceptance
     number* — not a hope.
   - Evidence: ME-1 3D GPU (R2 ≥1,500 MLUPS), ME-2 FP16 (T16 bands, ≥1.5×),
     ME-3 cluster (R3 64-rank ≥80%), ME-4 full-physics stirred-tank
     (M-Star-comparable degradation ratio) — all recorded in docs/PLAN.md
     "Performance roadmap". A vendor stating the number it will hit, in order, is
     the strongest possible expression of pillar P1.
8. **Conclusion + "evaluate it yourself"**
   - Claim: you can verify every number in this paper in one command.
   - Evidence: the reproduction CTA.
- **Appendices**: (A) full benchmark tables + configs + repro commands;
  (B) validation test inventory; (C) hardware/software provenance + STREAM/roofline;
  (D) claim→evidence traceability table (§7 below, expanded).

## 4. Benchmark design (design of experiment)

Two comparison families, both public-grade (§5). GPU head-to-head is **not**
possible today (no OSS 2D-GPU competitor runs on Apple Silicon; our 3D-GPU is
roadmap) — so GPU performance is reported as our own measured 2D Metal grid sweep,
with FluidX3D (published + same-machine OpenCL 3D) as *hardware-ceiling context*,
explicitly not a head-to-head.

### 4.1 Family A — LBM-vs-LBM, same machine, MLUPS (clean, fair)
- **Codes**: LBMFlow (CpuSimd) vs **OpenLB** (latest, CPU/OpenMP) vs **Palabos**
  (CPU/OpenMP). All built on the M5 Max (macOS, arm64).
- **Cases**: Taylor–Green periodic (pure kernel) and lid-driven cavity
  (with-walls). **D2Q9 @ 1024²** and **D3Q19 @ 128³**.
- **Metrics**: MLUPS at matched grid — single-thread and all-core, f32 and f64;
  warmup excluded; ≥5-run median + spread. Report **roofline triple**: measured
  STREAM bandwidth → theoretical MLUPS ceiling → measured MLUPS (the credibility
  multiplier — "we're at X% of the memory-bandwidth limit").
- **Caveats table** (mandatory, per BENCH_COMPARISON_DRAFT §8): stencil (D2Q9 vs
  D3Q19), precision (arithmetic/storage), streaming implementation, kernel-only vs
  full-solver, warmup/averaging window.
- **Feasibility risks (resolve in execution step 0)**:
  - OpenLB: make-based, OpenMP needs libomp (Apple clang) — CPU build expected
    feasible; GPU is CUDA-only (no Apple GPU path — that's why Family A is CPU).
  - Palabos: CMake + OpenMP/MPI — buildable on macOS with effort.
  - If a code will not build cleanly on macOS arm64 in the window, fall back to a
    documented Linux/Docker run on the same machine and **label the environment
    difference explicitly** (do not silently mix environments).

### 4.2 Family B — OpenFOAM time-to-solution (sales-resonant, method-different)
- **Code**: OpenFOAM (the OSS FVM incumbent), via native build or Docker on the
  same machine.
- **Case (fixed)**: Schäfer-Turek cylinder Re=20 → **Cd** (2D-1). Single defensible
  number with a standard reference value (Cd=5.5795); accuracy targets are easy to
  match across both solvers. (Cavity Re=1000 / Ghia RMS is a possible later addition,
  not in the first cut.)
- **Metric**: **wall-clock time to reach a target accuracy** — Cd within the
  Schäfer-Turek band (both solvers driven to the same Cd tolerance). "Time to a
  *validated* answer," not "our kernel is faster."
- **Mandatory caveat framing**: FVM (implicit, unstructured, large Δt) vs LBM
  (explicit, uniform grid, small Δt) are fundamentally different discretizations.
  Present as an engineering-question comparison with matched accuracy target,
  documented mesh/grid, solver settings, and machine. Never as MLUPS.

### 4.3 Family C — FluidX3D reference (context only, per PM)
- Report FluidX3D's published D3Q19 numbers **and** a same-machine OpenCL run on the
  M5 Max GPU as "what the hardware ceiling looks like in 3D," honestly noting
  LBMFlow's 3D GPU is roadmap (R2). Not counted as a win.

## 5. Public-grade measurement protocol (from BENCH_COMPARISON_DRAFT §11)

- **Environment**: idle machine — announce the window in TESTING_NOTES.md at start
  and finish; PM holds heavy dispatches during it (coordination agreed). Power
  connected, thermal steady, other apps closed. Record OS / wgpu / driver / rustc /
  compiler versions and memory config + nominal bandwidth.
- **Roofline first**: measure actual STREAM(-equivalent) bandwidth; compute the
  theoretical ceiling; report alongside every MLUPS number.
- **Timing**: exclude warmup N steps; ≥5 runs, median + spread.
- **Sweep**: grid-size sweep (2D 512²–4096², 3D 128³–memory limit) → saturation
  curve + the size at which saturation is reached.
- **Definition**: MLUPS = total cells × steps / walltime, I/O excluded — stated.
- **Two regimes reported separately**: kernel-only (empty box, minimal boundaries)
  vs representative scenario (obstacle + outputs on).
- **Competitors**: same machine, same case, documented build + settings; if
  environment must differ, label it.

## 6. Claim → evidence → reproduction traceability (skeleton; fills at write time)

| Claim (paper section) | Metric / number | Config | Repro command | Source |
|---|---|---|---|---|
| 2D GPU bandwidth-bound speed (§5) | 5.7–11.5 GLUPS (grid-dep.) | M5 Max Metal, D2Q9 f32, TRT | `cd crates/lbm-gpu-proto && cargo run --release` | GPU_EVALUATION.md §1 |
| Single-node CPU 2D (§5) | 1,183 MLUPS | 1024²/12T/f32, CpuSimd | `cargo run --release -p lbm-core --example bench_backends` | PERFORMANCE.md |
| CPU vs OpenLB/Palabos (§5) | *to measure* | §4.1 | *bench harness, to build* | this plan §4.1 |
| Time-to-solution vs OpenFOAM (§5) | *to measure* | §4.2 | *to build* | this plan §4.2 |
| MPI weak scaling n≤4 (§5) | 97–99% | 512²/rank, InProcess/MPI | `./scripts/test_mpi.sh` + bench | TESTING_NOTES (M-D) |
| f32 validation-grade (§4) | mom. err 2.8e-7; TGV f32≈f64 | deviation storage | validation suite | PERFORMANCE.md / T6 |
| Partition invariance bit-exact (§4) | max\|Δ\|=0.0 | InProcess/MPI splits | `cargo test --release t13*` | VALIDATION T13 |
| CPU↔GPU equivalence (§4) | ≤1e-5 (1e-4 pressure BC) | 6 configs, f32 | `cargo test --release --features gpu t14*` | VALIDATION T14 |
| Agent-native MCP (§6) | 7 tools, async E2E | MCP server | b1 audit repro | docs/skills/b1-capability-map.md |

## 7. Execution plan (after structure/message sign-off)

0. **Build feasibility probe** (small, no quiet window needed): confirm OpenLB /
   Palabos / OpenFOAM build+run on this machine (or pin Docker fallback), pick the
   exact cases + accuracy targets. Time-box; report a go/no-go per code.
1. **Announce quiet window** in TESTING_NOTES.md; PM holds heavy dispatches.
2. **Run the public-grade protocol in one sitting** (§5): roofline → LBMFlow sweeps
   → OpenLB/Palabos → OpenFOAM time-to-solution → FluidX3D context.
3. **Fill the traceability table (§6)** and the cross-code comparison table with
   condition caveats.
4. **Write** the paper (§3) in English, claims bounded by §2.
5. Feed leftover competitive facts back into BENCH_COMPARISON_DRAFT; keep the
   provisional→final number swap disciplined.

## 8. Decisions & open items

**Decided (2026-07-05)**
- PM: hybrid (C), lead differentiation + scoped honest measured performance; paper
  does not wait for M-E. Claim boundary = §2.
- Family B case: **Schäfer-Turek cylinder Re=20, Cd** (§4.2).
- Delivery: **Markdown draft first** (in-repo); PDF / web formatting decided after
  content is locked.
- Language: English. Measurement: public-grade (§5).

**Open**
- Message architecture (§1) + outline (§3): user requested revisions before sign-off
  — pending user's revision notes (the "固める" gate).
- Cluster access for R3 stays a user decision (paper does not wait on it).
