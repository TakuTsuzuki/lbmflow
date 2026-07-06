# Architecture V2 — 3D, GPU, and distribution from a single core

The V2 architecture aligned with COMPETITIVE_SPEC.md's required requirements (R1-R5).
V2 has replaced V1 (2D/CPU/single-domain) as the single core. This document records
the design principles the core is built on and the current orthogonal-axis set.

## 0. Design principles

1. **From one physics-kernel definition to every target**: dimension, lattice, precision,
   backend, and partitioning are orthogonal axes. Physics (collision, forces, boundaries)
   is written in one place, and axis combinations are expanded at compile time.
2. **Equivalence is the only source of truth**: every time a new axis is added, define a
   test first that "matches the existing configuration" (T13/T14/T16). Speed claims come after.
3. **The agent contract is invariant**: scenario JSON / MCP is the common entry point
   across all configurations (R4).

## 1. Layer structure

```
┌─────────────────────────────────────────────────────┐
│  scenario / CLI / MCP / GUI (contract layer: JSON, minimal change)  │
├─────────────────────────────────────────────────────┤
│  Solver Orchestrator (time evolution, diagnostics, probes, output)  │
├───────────────┬───────────────────┬─────────────────┤
│ Decomposition │  Physics Kernels   │  Diagnostics    │
│ Subdomain     │  collide+stream    │  reductions     │
│ HaloExchange  │  BCs / forces      │  (backend-side) │
├───────────────┴───────────────────┴─────────────────┤
│  Backend trait: CpuScalar | CpuSimd | Wgpu           │
├─────────────────────────────────────────────────────┤
│  Lattice trait: D2Q9 | D3Q19 | D3Q27  × f32/f64/f16   │
└─────────────────────────────────────────────────────┘
```

## 2. Orthogonal axes (current set)

| Axis | Values in trunk |
|---|---|
| Dimension × lattice | D2Q9, D3Q19, D3Q27 |
| Collision | BGK, TRT, Cumulant (central-moment) |
| Precision (compute × storage) | f32, f64, f16 storage (compute-in-f32) |
| Backend | CpuScalar, CpuSimd, Wgpu |
| Partition | Local, InProcess, MPI (feature) |

Turbulence closure: WALE LES (default; Smagorinsky retained as reference).

## 3. Definition of each abstraction

### 3.1 Lattice (compile-time constants)

```rust
pub trait Lattice: Copy + 'static {
    const D: usize;            // 2 | 3
    const Q: usize;            // 9 | 19 | 27
    const C: [[i8; 3]; Self::Q];   // velocities (z=0 for 2D)
    const W: [f64; Self::Q];
    const OPP: [usize; Self::Q];
    const CS2: f64;            // 1/3
    // Derived tables such as TRT pairs and per-face unknown sets are also
    // provided via const fn
}
pub struct D2Q9; pub struct D3Q19; pub struct D3Q27;
```
- Turns the lattice.rs principle "direction order is the single source of truth"
  into a trait and preserves it.
- Zou-He's face-normal parameterization is isomorphic in D3 (face unknowns = the
  directions where c·n>0).

### 3.2 Storage (SoA fixed, for GPU coalescing)

- `f[q][cell]` (q-major SoA). cell = z·(nx·ny) + y·nx + x (2D uses only z=0).
- **Deviation storage (f−w)** — f32 is verification-grade and is the precondition
  for f16.
- Precision is 2 axes of "compute precision × storage precision": (f32,f32) /
  (f64,f64) / (f32, f16 storage). f16 pack/unpack is kept inside the backend;
  the API surface presents f32.
- The moments cache (rho,u) is for diagnostics, visualization, and multiphase.
  It resides in backend-side memory, with explicit readback (no implicit sync).

### 3.3 Subdomain / HaloExchange (R3)

```rust
pub struct Subdomain { global_box: Box3, local: Box3, halo: usize /*=1*/,
                       neighbors: [Option<RankId>; 6/*faces*/ + edges…] }
pub trait HaloExchange {
    /// Per face, "outgoing distributions only" (D2Q9: 3/face, D3Q19: 5/face,
    /// D3Q27: 9/face) + also exchange ψ for multiphase.
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

### 3.4 Backend

```rust
pub trait Backend<L: Lattice> {
    type Field;   // device-resident f / moments / mask
    fn step(&mut self, dom: &Subdomain, fields: &mut Fields<Self>, params: &StepParams);
    fn reduce(&self, kind: Reduction) -> f64;
    fn read_moments(&self, out: &mut HostMoments);   // explicit readback
}
```
- **CpuScalar / CpuSimd**: reference and SIMD-optimized CPU backends (fused
  collide+stream+moments in step_band on CpuSimd).
- **Wgpu**: WGSL kernel (fused collide+stream, ping-pong), FP16 storage via
  shader-f16. Quiet-window MLUPS 2791-2813 @ 192³ D3Q19 (2026-07-06).
- The boundary-condition pass (Zou-He, etc.) is per-backend as "a small kernel
  touching only edge cells." The formulas are generated from the Lattice trait's
  face table, sharing the same definition across CPU/GPU.

### 3.5 Migration history — all steps complete (2026-07-05)

V1 (2D/CPU) retired; `crates/lbm-core` is now the V2 implementation described here.
Proof of equivalence with V1 is recorded as frozen values in `tests/v1_match.rs`'s
header right before its removal (branch history). The `compat` module remains as
a public API, and the scenario / CLI / wasm 2D paths use it.

## 4. Extension of scenario contract v1 (backward-compatible)

```jsonc
{
  "grid": { "nx": 256, "ny": 256, "nz": 128 },        // nz added (omit = 2D)
  "physics": { "precision": "f32", "storage": "f16" }, // storage added
  "compute": {                                          // all optional
    "backend": "auto | cpu | gpu",
    "decompose": { "ranks": [2, 2, 1] }                 // hint for MPI execution
  },
  "outputs": [ { "field": "q-criterion", "format": "vtk" } ]
}
```
- Existing fields are unchanged (stays deny_unknown_fields, additions only).
- GUI remains 2D-only (3D goes via CLI/agent first).

## 5. Verification map (what to write to call it "done")

| Axis | Test | Content |
|---|---|---|
| Partition | T13 | 1×1 vs 2×2 vs 4×1 vs 1×4 (+3D: 2×2×2) match (f64 ≤1e-12) |
| Backend | T14 | CPU vs Wgpu same scenario, f32 relative ≤1e-5, diagnostics match |
| 3D physics | T15 | TGV3D convergence order, sphere drag (Re=100 Cd≈1.09 band), 3D cavity, channel Re_τ=180 (LES) |
| Precision | T16 | Quantified degradation of f16 storage (frozen band for TGV/cavity) |
| Regression | existing T1-T12 | Green unmodified via the V1 facade |

## 6. Risks and honest constraints

- **wgpu has no f64** → verification that needs f64 is the CPU backend's
  responsibility (T14 compares in f32).
- **MPI measurements are local-only** here (multi-rank functional verification).
  True weak scaling awaits cluster access (COMPETITIVE_SPEC §5).
- **In-place streaming (Esoteric-Pull)** halves memory and is a candidate for
  M-E; must preserve the stream contract (§3.3).
