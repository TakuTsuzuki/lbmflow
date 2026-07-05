# Architecture V2 — 3D, GPU, and distribution from a single core (2026-07-05)

A design aligned with COMPETITIVE_SPEC.md's required requirements (R1-R5). Replaces the
current V1 (2D/CPU/single-domain) in stages **while keeping the verification suite green**.

## 0. Design principles

1. **From one physics-kernel definition to every target**: dimension, lattice, precision,
   backend, and partitioning are orthogonal axes. Physics (collision, forces, boundaries)
   is written in one place, and axis combinations are expanded at compile time
2. **Equivalence is the only source of truth**: every time a new axis is added, define a
   test first that "matches the existing configuration" (T13/T14/T16). Speed claims come after
3. **The agent contract is invariant**: scenario JSON / MCP is the common entry point
   across all configurations (R4)

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
│  Backend trait: CpuSimd | Wgpu | (Cuda: feature)     │
├─────────────────────────────────────────────────────┤
│  Lattice trait: D2Q9 | D3Q19 | (D3Q27)  × Real/f16   │
└─────────────────────────────────────────────────────┘
```

## 2. Definition of each abstraction

### 2.1 Lattice (compile-time constants)

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
pub struct D2Q9; pub struct D3Q19;
```
- Turns the current lattice.rs principle "direction order is the single source of truth"
  into a trait and preserves it
- Zou-He's face-normal parameterization (proven in V1) is isomorphic in D3 as well
  (face unknowns = the 5 directions where c·n>0)

### 2.2 Storage (SoA fixed, for GPU coalescing)

- `f[q][cell]` (q-major SoA). cell = z·(nx·ny) + y·nx + x (2D uses only z=0)
- **Keep deviation storage (f−w)** (V1 proved f32 verification-grade; also a
  precondition for FP16)
- Precision is 2 axes of "compute precision × storage precision": (f32,f32) /
  (f64,f64) / **(f32, f16 storage)** (R2). f16 pack/unpack is kept inside the backend;
  the API surface presents f32
- The moments cache (rho,u) is for diagnostics, visualization, and multiphase. It
  resides in backend-side memory, with explicit readback to the host (never create
  implicit synchronization)

### 2.3 Subdomain / HaloExchange (R3)

```rust
pub struct Subdomain { global_box: Box3, local: Box3, halo: usize /*=1*/,
                       neighbors: [Option<RankId>; 6/*faces*/ + edges…] }
pub trait HaloExchange {
    /// Per face, "outgoing distributions only" (D2Q9: 3/face, D3Q19: 5/face)
    /// + also exchange ψ for multiphase
    fn exchange(&mut self, field: &mut BackendField, plan: &HaloPlan);
}
impl: LocalPeriodic (single domain = current behavior) / InProcess (inter-thread,
      for T13) / Mpi (rsmpi)
- (R-Phase 1) `HaloExchange::SCOPE: ExchangeScope { Local, Remote }` — building a
  single-part owner (`only=Some(part)`) with a Local exchange is a construction
  error (silent self-wrap prevention, spec A-5).
- (R-Phase 1) `Backend::stream` contract, pinned by tests/stream_contract.rs:
  streaming must NOT write open-face unknown slots (ConvectiveOutflow memory
  depends on it; GPU realizes it via the edge stash). Any in-place streaming
  (M-E candidate) must preserve this contract or replace the mechanism.
```
- Split the step structure into a **2-pass internal→boundary** structure, so halo
  communication and internal computation can be overlapped (V1's row-partitioned loop
  is compatible with this split)
- Global diagnostics (total_mass, etc.) do a backend-internal reduce → inter-rank
  Allreduce
- Rims/obstacles/probes are distributed into Subdomain-local data (mechanically
  splittable because V1 already makes "BC = data")

### 2.4 Backend

```rust
pub trait Backend<L: Lattice> {
    type Field;   // device-resident f / moments / mask
    fn step(&mut self, dom: &Subdomain, fields: &mut Fields<Self>, params: &StepParams);
    fn reduce(&self, kind: Reduction) -> f64;
    fn read_moments(&self, out: &mut HostMoments);   // explicit readback
}
```
- **CpuSimd**: the current rayon implementation, made SoA + SIMD (absorbing the
  results of the phase9-perf branch)
- **Wgpu**: WGSL kernel (fused collide+stream, ping-pong). Adoption and tuning
  guidance are settled by the phase9-wgpu evaluation results. FP16 storage via
  shader-f16
- **Cuda** (feature, for later / NVIDIA supercomputers): an additional implementation
  of the same trait. MPI+CUDA GPUDirect is M-D or later
- The boundary-condition pass (Zou-He, etc.) is implemented per backend as "a small
  kernel touching only edge cells." The formulas are generated from the Lattice
  trait's face table, sharing the same definition across CPU/GPU

### 2.5 Migration strategy (R5: don't hold the 2D suite hostage) — **all steps complete (2026-07-05)**

1. Set up `lbm-core2` and implement the V2 abstractions (V1 is frozen, becomes the
   reference implementation) ✅
2. Implement a **V1 API facade** on top of lbm-core2 (all public methods of
   `Simulation<T>`). Make it the first milestone that the existing 56+ tests,
   wasm, and CLI all pass unchanged ✅
3. Commission T13 (partition invariance) / T14 (backend equivalence) to codex,
   establishing V2's adversarial verification ✅
4. Implement 3D (D3Q19) as an added Lattice → T15 (3D physics) ✅
5. Once stable, rename lbm-core2 → lbm-core and remove V1 ✅
   (Done on 2026-07-05. `crates/lbm-core` is now the V2 implementation described in
   this document. The proof of equivalence with V1 is recorded as frozen values in
   the `tests/v1_match.rs` header right before its removal — see branch history.
   The `compat` module remains as a public API, and the scenario / CLI / wasm 2D
   paths use it)

## 3. Extension of scenario contract v1 (backward-compatible)

```jsonc
{
  "grid": { "nx": 256, "ny": 256, "nz": 128 },        // nz added (omit = 2D)
  "physics": { "precision": "f32", "storage": "f16" }, // storage added
  "compute": {                                          // new (all optional)
    "backend": "auto | cpu | gpu",
    "decompose": { "ranks": [2, 2, 1] }                 // hint for MPI execution
  },
  "outputs": [ { "field": "q-criterion", "format": "vtk" } ]  // add 3D visualization fields
}
```
- Existing fields are unchanged (stays deny_unknown_fields, additions only)
- GUI remains 2D-only (3D goes via CLI/agent first. A cross-section/isosurface GUI is
  judged at M-F or later)

## 4. Verification map (what to write to call it "done")

| Axis | Test | Content |
|---|---|---|
| Partition | T13 | 1×1 vs 2×2 vs 4×1 vs 1×4 (+3D: 2×2×2) match (f64 ≤1e-12) |
| Backend | T14 | CPU vs Wgpu same scenario, f32 relative ≤1e-5, diagnostics match |
| 3D physics | T15 | TGV3D convergence order, sphere drag (Re=100 Cd≈1.09 band), 3D cavity, (post-LES) channel Re_τ=180 |
| Precision | T16 | Quantified degradation of f16 storage (frozen band for TGV/cavity) |
| Regression | existing T1-T12 | Green unmodified via the V1 facade |

## 5. Risks and honest estimate

- **wgpu has no f64** → verification that needs f64 is accepted as the CPU backend's
  responsibility (T14 compares in f32)
- **The upper bound of MPI measurements is local** (up to multi-rank functional
  verification only. True weak scaling awaits cluster access
  (COMPETITIVE_SPEC §5))
- **In-place streaming such as Esoteric-Pull** (halving memory) is a candidate for
  introduction at M-E, aiming for FluidX3D-class. First nail down correctness with
  plain ping-pong
- Effort estimate (assuming autonomous-agent parallelism): M-B core V2 ≈ 2-4 nights,
  M-C 3D ≈ 2-3 nights, M-D MPI ≈ 2-3 nights + cluster measurement separately
