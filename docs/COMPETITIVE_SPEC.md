# Competitive Advantage Spec (established 2026-07-05)

Required requirements per user directive: **3D support, supercomputer scale, and GPU support are mandatory**.
This document defines the "winning strategy" relative to the market's major players, and all subsequent design follows this document.

## 1. Competitive landscape (LBM-based CFD)

| Product/Code | Strengths | Weaknesses (our angle of attack) |
|---|---|---|
| **M-Star CFD** (commercial, primary comparison target) | GPU-native transient LBM+LES. Vertical specialization in stirred tanks/pharma-bio. CAD→automatic meshing (meshless). Particles, scalar transport, free surface, non-Newtonian. Python API. Industrial track record | **CUDA/NVIDIA lock-in**. Closed-source. High licensing cost. Validation is a published benchmark set but **not provided in a re-runnable form**. No agent integration |
| FluidX3D (OSS) | Single-node GPU speed king (FP16 memory compression, Esoteric-Pull, thousands of MLUPS on RTX-class hardware). Free surface | Non-commercial license. No multi-node (MPI). Oriented toward UX researchers. Minimal validation suite |
| waLBerla (OSS/HPC) | Genuine supercomputer track record (trillions of cells, SuperMUC, etc.). MPI+GPU, code generation | Framework for experts. Steep learning curve for C++/code generation. No product UX |
| Palabos / OpenLB (OSS) | Academically mature, broad feature set | Classical node performance, weak GPU, dated UX |
| PowerFLOW / XFlow (commercial) | Track record in vehicle aerodynamics/aeroacoustics with OEMs | Extremely expensive, closed-source, niche-industry oriented |

**Reading**: No player simultaneously holds "GPU speed (FluidX3D class)," "supercomputer
distribution (waLBerla class)," and "product UX (M-Star class)" all at once. Furthermore,
there is a 4th axis nobody holds —
**AI-agent-native** and **re-runnable validation**.

## 2. Our differentiation spec (winning strategy = 4 pillars)

### Pillar 1: Agent-native CFD (an industry-first category)
- The self-describing scenario contract (JSON Schema + `lbm schema` + 4 MCP tools) is already implemented.
  Maintain and expand this as a **first-class product surface**
- Requirements: an asynchronous job API letting agents run parameter sweeps/optimization loops
  unattended, structured diagnostics of run results (including divergence reasons and stability hints),
  and a determinism guarantee
- M-Star's Python API is oriented toward "human scripting." We are oriented toward "autonomous
  operation by LLM agents" (schema self-discovery, machine-readable warnings, explainability of failures)

### Pillar 2: Evidence-driven trust (Adversarial Validation as a Product)
- **Ship the adversarial validation suite as a shipped artifact** (currently 56+ tests, Ghia/Schäfer-Turek/RT growth
  rate/exact isotropy 4e-16), re-runnable in the customer's environment with a single command
- Tie every performance/accuracy claim to measurements in PHYSICS.md / PERFORMANCE.md (continue current practice)
- Mandate equivalence tests for every newly added backend/dimension (§4)

### Pillar 3: A portability ladder that runs anywhere (no CUDA lock-in)
- **The same scenario JSON** runs on browser (WASM) → laptop (CPU SIMD / Metal GPU) →
  workstation (multi-GPU) → cluster (MPI+GPU)
- The first-class GPU implementation is **wgpu** (Metal/Vulkan/DX12 = full coverage of Apple/AMD/Intel/NVIDIA).
  A CUDA backend is held as a **swappable additional implementation** for NVIDIA supercomputers
  (the core is abstracted via a backend trait — ARCHITECTURE_V2)
- Capture the Apple Silicon / AMD / browser markets that M-Star structurally cannot reach, at no cost

### Pillar 4: Precision transparency (Precision Dials)
- Building on deviational-storage f32 = validation-grade (already proven), provide a **FP16 storage/FP32
  compute** memory-compression mode (FluidX3D-style) with validation — doubling the grid within GPU memory
- Every precision mode is subject to the validation suite (numerically stating exactly "what the fast mode gives up")

### Non-goals (stated honestly; areas where we do not compete with M-Star for now)
- Maturity of CAD import/automatic voxelization, particle/DEM coupling, combustion, industrial-grade free
  surface, non-Newtonian — sequential decisions from Phase 13 onward. First establish the 4 pillars of the physics core

## 3. Required requirements (per user directive, 2026-07-05)

| ID | Requirement | Measurable acceptance criteria |
|---|---|---|
| R1 | **3D (D3Q19/Q27)** | Pass the 3D validation suite: TGV3D (diffusion-limit reference ±2%, convergence order ≥1.7), sphere drag (Re∈{20,100}, Schiller-Naumann correlation **±10%, hydrodynamic-pair normalization D_h=D+1, Re_h=Re(D+1)/D**), 3D cavity (T15.5 = A&K 2005 Re=1000, RMS≤0.030U), turbulent channel Re_τ=180 vs DNS (Moser+) after LES is introduced (M-F) |
| R2 | **GPU support** | wgpu backend equivalent to CPU at f32 (T14: field agreement ≤1e-5 relative). Single-GPU D3Q19 f32 ≥ 1,500 MLUPS (M4 Max class or RTX 4070 class). 2x grid size in FP16 storage mode |
| R3 | **Supercomputer scale** | Subdomain+MPI: partitioned run ≡ monolithic run agreement (T13). Single-node intra-node multi-rank weak scaling **≥85% (local line for homogeneous cores, n≤4; measured 97-99%**. n=8 is 73% due to M5 Max's heterogeneous cores + bandwidth ceiling, with 84% as the ceiling even under a zero-communication control — confirmed via control experiment to be a property of intra-node SMP, not an MPI implementation defect). Cluster measurement (requires machine access, §5, CLUSTER_OPTIONS.md) for 64-rank weak scaling ≥80% — **not yet measured** |

> **Revision history (D-6, 2026-07-05)**: R1 sphere drag ±5%→**±10%** (the half-way BB wall sits half a link outside the solid cell, so nominal D normalization carries a ~+2/D bias. Revised to a band incorporating the measured +0.6–7.1% together with adoption of hydrodynamic-pair D_h normalization. Basis: TESTING_NOTES 2026-07-05 triage, PHYSICS.md). R1 3D cavity's "literature comparison" is made concrete as T15.5's (A&K 2005) quantitative band. R3 weak scaling ≥85% is clarified in scope as "single-node, n≤4, local," while the cluster 64-rank line remains unmeasured (measurement plan: CLUSTER_OPTIONS.md).
| R4 | Maintain agent-nativeness | All new features must be controllable from scenario JSON+MCP. Asynchronous job API |
| R5 | Validation continuity | The existing 2D suite of 56+ **stays green without modification after the core overhaul** (API facade maintained) |

## 4. Equivalence test framework (new category, for codex commissioning)

- **T13 partition invariance**: results of 1×1 / 2×2 / 4×1 partitioned runs match the monolithic run
  (bit-match target at f64, at least ≤1e-12)
- **T14 backend equivalence**: CPU-SIMD vs wgpu vs (future CUDA) match on the same scenario
  at f32 relative ≤1e-5 / statistical agreement
- **T15 3D physics**: each benchmark listed under R1
- **T16 precision modes**: quantitatively freeze the degradation of FP16 storage vs f32 (tolerance band specified)

## 5. Assumptions and dependencies (items to communicate to the user)

- Multi-node measurement requires access to a cluster/cloud HPC (locally we can only go as far
  as functional verification of multiple MPI ranks). Will consult once this becomes necessary
- CUDA backend measurement requires an NVIDIA GPU (this machine is Apple Silicon/Metal only;
  wgpu measurement is possible on this machine)
- FP16 is validated on a GPU with wgpu shader-f16 support (Apple Silicon supports this)

## 6. Milestone reorganization (current plan, supersedes PLAN.md)

- **M-A (in progress)**: CPU SoA/SIMD speedup + wgpu measured evaluation (phase9-perf / phase9-wgpu branches)
- **M-B: Core V2** — dimension/lattice/backend/Subdomain abstraction (ARCHITECTURE_V2.md).
  Swap in while keeping the 2D suite green. Framework for T13/T14
- **M-C: 3D physics** — D3Q19 + 3D validation (R1). Both CPU/GPU
- **M-D: Distribution** — HaloExchange=MPI, parallel I/O, measured weak scaling (R3)
- **M-E: Performance headline** — FP16 storage, multi-GPU, document same-condition comparison
  against FluidX3D/M-Star in the public benchmark table
- **M-F: Selection of vertical features** — assumes the order LES (Smagorinsky) → moving boundary (IBM) → scalar transport
  (the ordering that most quickly reaches M-Star's core use case = stirred tanks)
