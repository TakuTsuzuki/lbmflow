# Architecture V2 - 3D, GPU, and distribution from a single core

The V2 architecture aligned with COMPETITIVE_SPEC.md's required requirements (R1-R5).
V2 has replaced V1 (2D/CPU/single-domain) as the single core. This document records
the design principles the core is built on and the current orthogonal-axis set.

Status terms used below:

- **Implemented in core**: available through `crates/lbm-core` APIs.
- **Product path exposed**: available through scenario JSON / CLI / MCP / GUI wiring.
- **Design target**: intended architecture, not yet fully implemented or exposed; see
  the PLAN.md queue before treating it as a product capability.

## 0. Design principles

1. **From one physics-kernel definition to every target**: dimension, lattice, precision,
   backend, and partitioning are orthogonal axes. Physics (collision, forces, boundaries)
   is written in one place, and axis combinations are expanded at compile time.
2. **Equivalence is the only source of truth**: every time a new axis is added, define a
   test first that "matches the existing configuration" (T13/T14/T16). Speed claims come after.
3. **The agent contract is staged**: scenario JSON / MCP is the product entry
   point for exposed configurations (R4). Core-only axes are exposed through
   that contract as REV-1 / follow-on schema work lands.

## 1. Layer structure

```
┌─────────────────────────────────────────────────────┐
│  scenario / CLI / MCP / GUI (contract layer: JSON, staged exposure) │
├─────────────────────────────────────────────────────┤
│  Solver Orchestrator (time evolution, diagnostics, probes, output)  │
├───────────────┬───────────────────┬─────────────────┤
│ Decomposition │  Physics Kernels   │  Diagnostics    │
│ Subdomain     │  collide+stream    │  reductions     │
│ HaloExchange  │  BCs / forces      │  (backend-side) │
├───────────────┴───────────────────┴─────────────────┤
│  Backend trait: CpuScalar | CpuSimd | Wgpu           │
├─────────────────────────────────────────────────────┤
│  Lattice trait: D2Q9 | D3Q19 | D3Q27  × f32/f64      │
│  GPU storage mode: f32 | f16 (compute remains f32)   │
└─────────────────────────────────────────────────────┘
```

## 2. Orthogonal axes (current set)

### 2.1 Core implementation reality

| Axis | Implemented in `crates/lbm-core` |
|---|---|
| Dimension × lattice | D2Q9, D3Q19, D3Q27 |
| Collision | BGK, TRT, Cumulant / cascaded central-moment (`CollisionKind`) |
| Arithmetic precision | f32, f64 (`Real`) |
| Distribution storage | CPU: f32/f64 deviation storage matching arithmetic; GPU: f32 or f16 distribution buffers with f32 arithmetic (`GpuStorage`) |
| Backend | CpuScalar, CpuSimd, Wgpu (`feature = "gpu"`) |
| Partition / halo | LocalPeriodic, InProcess, MPI (`feature = "mpi"`) |
| Open faces | D2Q9, D3Q19, and (since 2026-07-07) D3Q27 velocity-inlet / pressure-outlet closures are supported on CPU; D3Q27 outflow/convective and GPU open faces are rejected explicitly |

WALE LES is implemented as a core solver-level relaxation-field driver. Scenario-level
LES controls, Smagorinsky selection, and the Re_tau DNS acceptance line remain design
targets / validation queue items unless a caller wires the core API directly.

### 2.2 Product-path exposure today

Scenario JSON currently exposes a narrower surface than the core:

| Axis | Product-path exposure today |
|---|---|
| Dimension × lattice | `grid.nz <= 1` -> 2D D2Q9; `grid.nz > 1` -> 3D D3Q19. D3Q27 is not selectable from scenario JSON today. |
| Collision | `physics.collision` accepts `bgk`, `trt`, or `cumulant` (landed 2026-07-07). Cumulant is honored on the 3D D3Q19 native CPU path only; 2D compat and GPU paths reject it with explicit errors. |
| Precision / storage | `physics.precision` accepts `f32` or `f64`. `compute.storage` accepts `f32` or `f16` (landed 2026-07-07); `f16` is honored only for 2D D2Q9 GPU scenarios with a SHADER_F16 adapter, otherwise rejected with explicit errors. |
| Backend | `compute.backend` accepts `auto`, `cpu`, or `gpu`; explicit GPU requests are honored or rejected. Current GPU scenario dispatch is constrained to f32 and rejects unsupported 3D combinations. |
| Partition | Scenario JSON does not expose MPI decomposition hints today. MPI is a core/CLI validation path, not a general scenario contract field. |

## 3. Definition of each abstraction

### 3.1 Lattice (compile-time constants)

```rust
pub trait Lattice: Copy + 'static {
    const D: usize;            // 2 | 3
    const Q: usize;            // 9 | 19 | 27
    const C: &'static [[i8; 3]];   // velocities (z=0 for 2D), length Q
    const W: &'static [f64];       // length Q
    const OPP: &'static [usize];   // length Q
    const CS2: f64;            // 1/3
    // Derived tables such as TRT pairs and per-face unknown sets are also
    // provided via const fn
}
pub struct D2Q9; pub struct D3Q19; pub struct D3Q27;
```
- Turns the lattice.rs principle "direction order is the single source of truth"
  into a trait and preserves it.
- Face unknowns are derived from the lattice table as the directions where
  `c dot n_in > 0`: 3 per D2Q9 face, 5 per D3Q19 face, 9 per D3Q27 face.
  Implemented open-face kernels currently accept only the 3-unknown and
  5-unknown cases plus the D3Q27 9-unknown NEBB closure for velocity inlet /
  pressure outlet (landed 2026-07-07); D3Q27 outflow/convective remain open.

### 3.2 Storage (SoA fixed, for GPU coalescing)

- `f[q][cell]` (q-major SoA). cell = z·(nx·ny) + y·nx + x (2D uses only z=0).
- **Deviation storage (f-w)** - f32 is verification-grade and is the precondition
  for f16.
- Precision is 2 axes of "compute precision × storage precision": (f32,f32) /
  (f64,f64) / GPU-only (f32 compute, f16 distribution storage). f16 pack/unpack
  is kept inside the wgpu backend; the core API surface presents f32 for that
  mode. Scenario JSON exposes this axis as `compute.storage: f32 | f16`
  (2D D2Q9 GPU scenarios only, landed 2026-07-07).
- The moments cache (rho,u) is for diagnostics, visualization, and multiphase.
  It resides in backend-side memory, with explicit readback (no implicit sync).

### 3.3 Subdomain / HaloExchange (R3)

```rust
pub struct Subdomain { global_box: Box3, local: Box3, halo: usize /*=1*/,
                       neighbors: [Option<RankId>; 6/*faces*/ + edges…] }
pub trait HaloExchange {
    /// Per face, outgoing distributions only (D2Q9: 3/face, D3Q19: 5/face,
    /// D3Q27: 9/face).
    fn exchange(&mut self, field: &mut BackendField, plan: &HaloPlan);
}
impl: LocalPeriodic (single domain) / InProcess (inter-thread, for T13) / Mpi (rsmpi)
- (R-Phase 1) `HaloExchange::SCOPE: ExchangeScope { Local, Remote }` — building a
  single-part owner (`only=Some(part)`) with a Local exchange is a construction
  error (silent self-wrap prevention, spec A-5).
- (R-Phase 1) `Backend::stream` contract, pinned by tests/stream_contract.rs:
  streaming must NOT write open-face unknown slots (ConvectiveOutflow memory
  depends on it; GPU realizes it via the edge stash). Any in-place streaming
  must preserve this contract or replace the mechanism.
```
- Step structure is 2-pass internal→boundary so halo communication and internal
  computation can overlap.
- Global diagnostics do a backend-internal reduce → inter-rank Allreduce.
- Rims/obstacles/probes are distributed into Subdomain-local data (BCs are data).
- Extra scalar/multiphase exchanges are design targets for later multi-distribution
  coupling, not a general implemented halo-exchange field today.

### 3.4 Backend

```rust
pub trait Backend<L: Lattice, T: Real> {
    type Fields;   // backend-owned f / moments / mask / side buffers
    fn alloc(&self, sub: &Subdomain) -> Self::Fields;
    fn stage_in(&self, sub: &Subdomain, fields: &mut Self::Fields, host: &SoaFields<T>);
    fn stage_out(&self, sub: &Subdomain, fields: &Self::Fields, host: &mut SoaFields<T>);
    fn exchange_f<H: HaloExchange<T>>(&mut self, exchange: &H, subs: &[Subdomain],
                                      fields: &mut [Self::Fields]);
    fn collide(&mut self, sub: &Subdomain, fields: &mut Self::Fields, p: &StepParams<T>);
    fn stream(&mut self, sub: &Subdomain, fields: &mut Self::Fields, p: &StepParams<T>,
              range: CellRange) -> [T; 3];
    fn swap(&mut self, fields: &mut Self::Fields);
    fn apply_open_faces(&mut self, sub: &Subdomain, fields: &mut Self::Fields,
                        p: &StepParams<T>);
    fn update_moments(&mut self, sub: &Subdomain, fields: &mut Self::Fields,
                      p: &StepParams<T>);
    fn reduce(&self, sub: &Subdomain, fields: &Self::Fields, p: &StepParams<T>,
              kind: Reduction) -> f64;
    fn read_moments(&self, fields: &Self::Fields, out: &mut HostMoments<T>);
}
```
- **CpuScalar / CpuSimd**: reference and SIMD-optimized CPU backends (fused
  collide+stream+moments in step_band on CpuSimd).
- **Wgpu**: WGSL backend for `T = f32`, with generated kernels over the lattice
  and selectable f32/f16 distribution storage. Current product-path GPU dispatch
  is narrower than the core backend and rejects unsupported scenario combinations.
- The boundary-condition pass (Zou-He, etc.) is per-backend as "a small kernel
  touching only edge cells" for supported open-face lattices. D3Q27 open-face
  formulas are not implemented yet.

### 3.5 Migration history — all steps complete (2026-07-05)

V1 (2D/CPU) retired; `crates/lbm-core` is now the V2 implementation described here.
Proof of equivalence with V1 is recorded as frozen values in `tests/v1_match.rs`'s
header right before its removal (branch history). The `compat` module remains as
a public API, and the scenario / CLI / wasm 2D paths use it.

## 4. Scenario contract status

```jsonc
{
  "grid": { "nx": 256, "ny": 256, "nz": 128 },        // nz added (omit = 2D)
  "physics": { "precision": "f32" },                   // f32 | f64 today
  "compute": {                                          // all optional
    "backend": "auto | cpu | gpu"
  },
  "outputs": [ { "field": "q-criterion", "format": "vtk" } ]
}
```
- Existing fields are unchanged (stays deny_unknown_fields, additions only).
- `grid.nz` is implemented, but scenario selection maps 3D to D3Q19 today.
- `physics.collision` accepts `bgk`/`trt`/`cumulant` in the scenario schema
  (cumulant honored on 3D D3Q19 CPU only; landed 2026-07-07). `compute.storage:
  "f16"` is landed for 2D D2Q9 GPU scenarios. The run manifest records the
  actually-used backend/lattice/collision/precision/storage (`provenance`).
- Scenario-level lattice selection and MPI `decompose` hints remain design
  targets (PLAN.md REV-1 residual), not implemented product-path fields today.
- GUI remains 2D-only (3D goes via CLI/agent first).

## 5. Verification map (what to write to call it "done")

| Axis | Test | Content |
|---|---|---|
| Partition | T13 | 1×1 vs 2×2 vs 4×1 vs 1×4 (+3D: 2×2×2) match (f64 ≤1e-12) |
| Backend | T14 | CPU vs Wgpu same scenario, f32 relative ≤1e-5, diagnostics match |
| 3D physics | T15 | TGV3D convergence order, sphere drag (Re=100 Cd≈1.09 band), 3D cavity; D3Q27 currently periodic / closed-wall only |
| Precision | T16 | Quantified degradation of GPU f16 storage (frozen band for D2Q9 wgpu f32-vs-f16 scenarios; broader f16 coverage is capacity/perf characterization) |
| Regression | existing T1-T12 | Green unmodified via the V1 facade |

## 6. Risks and honest constraints

- **wgpu has no f64** → verification that needs f64 is the CPU backend's
  responsibility (T14 compares in f32).
- **MPI measurements are local-only** here (multi-rank functional verification).
  True weak scaling awaits cluster access (COMPETITIVE_SPEC §5).
- **In-place streaming (Esoteric-Pull)** halves memory and is a candidate for
  M-E; must preserve the stream contract (§3.3).
- **Product-path lag**: scenario JSON intentionally trails core capabilities.
  Cumulant, f16 storage, D3Q27 selection, and decomposition hints must not be
  documented as product-path capabilities until the schema and manifest wiring
  land.
