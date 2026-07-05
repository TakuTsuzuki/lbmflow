# GPU_EVALUATION.md — Phase 9c: wgpu compute backend measured evaluation

**Conclusion: recommend adoption (conditional yes).** For D2Q9 f32, **1024² 7,584 MLUPS / 2048² 6,975 MLUPS**
(**19–16x** vs full-core CPU, 42x at 512²). Verified with the same initial conditions and 2000 steps as
lbm-core: **L∞ relative 7.0e-6** (1/14 of the 1e-4 tolerance) — the GPU's physics is effectively the same
trajectory as the CPU's. The conditions are two: "f32 only" and "data GPU-resident design" (see risks below).

- Measurement date: 2026-07-05. Environment: Apple M5 Max (GPU, Metal) / CPU 18 cores / 128 GB,
  macOS, rustc 1.93.0, wgpu 26.0.1, `--release` (thin LTO, cu=1).
- Prototype: `crates/lbm-gpu-proto` (a throwaway evaluation crate **outside the workspace**;
  already registered in the root Cargo.toml's exclude — the production build graph does not
  depend on wgpu).
- Reproduce: `cd crates/lbm-gpu-proto && cargo run --release` (runs verification+benchmark
  together, outputs a markdown table. Use `--cpu-only` / `--gpu-only` to measure separately).
- Note on measurement conditions: since other agents' verification suites run intermittently
  on this machine, the CPU baseline was taken in "a window confirmed to have zero heavy-load
  processes before and after" (confirmed to match the already-recorded values of 285/376-381
  MLUPS in PERFORMANCE.md at 512²/1024²). GPU values are nearly insensitive to CPU load,
  within ±10% across 5 runs.

## 1. Measured: MLUPS (TRT, f32, effective value including submit→wait, after warmup)

| grid | GPU MLUPS | CPU MLUPS (same day/machine, lbm-core f32/TRT full-core) | ratio |
|---|---|---|---|
| 512² | **12,152** | 290 | **41.9x** |
| 1024² | **7,584** | 397 | **19.1x** |
| 2048² | **6,975** | 441 | **15.8x** |

- Also **~20x** against PERFORMANCE.md's representative value (1024² ≈ 380 MLUPS).
- Effective memory traffic is 72 B/cell/step (read+write of 9 directions of f32), so
  that's **~502 GB/s** at 2048² and ~546 GB/s at 1024². The jump to 12.1 GLUPS
  (equivalent to 875 GB/s) at 512² alone is because the ~19 MB working set fits
  on-chip cache (SLC).
  In other words this is **entirely bandwidth-bound**, and compute (TRT) is essentially
  free — the same conclusion as on CPU also holds on GPU.
- ~7,200 steps/s at 1024². At GUI 60 fps, that's **~120 steps per frame**
  (the current WASM CPU does ~600 steps/s at 256×128).
- Reference (extrapolation to D3Q19): at 152 B/cell/step, ~500 GB/s → **~3,300 MLUPS**,
  ~1,500 steps/s at 128³ ≈ 2.1M cells. 3D (Phase 10) is exactly where GPU pays off most.

## 2. Verification: comparison against lbm-core with identical initial conditions, 2000 steps (f32 vs f32)

TGV (periodic boundary, nu=0.02, u0=0.05, TRT magic 3/16). The GPU starts from the
**same f32 initial distribution** as the CPU (faithfully reproducing the second-order
consistent feq+f_neq initialization on the host side).

| grid | L∞(Δu)/max‖u‖ | L2 relative diff | GPU vs analytical L2 | CPU vs analytical L2 | verdict |
|---|---|---|---|---|---|
| 256² | 6.21e-6 | 3.60e-6 | 8.834e-4 | 8.826e-4 | PASS |
| 512² | 7.01e-6 | 3.87e-6 | 4.487e-4 | 4.483e-4 | PASS |

- **14x margin** against the pass criterion L∞ < 1e-4. The difference is only
  f32 rounding accumulation (from Metal compiler FMA/reassociation), and accuracy
  against the analytical solution matches the CPU to the 4th significant digit.
  **Carrying the deviation-storage scheme (keeping f−w) over to the GPU as-is is what
  makes this work** — since the stationary background is exactly zero, f32 rounding
  only acts on the fluctuation scale.
- Equivalence of operator ordering: for the fused kernel (pull→collide), one step is
  C∘S, whose composition order is reversed from the CPU's S∘C. Applying collision once
  on the host **to the initial state before upload** makes the GPU state after k steps
  equal to C(cpu_k); since collision is invariant for density/momentum, the velocity
  field can be compared 1:1 ((C∘S)^k∘C = C∘(S∘C)^k). This identity can also be used
  for regression tests at adoption time.

## 3. Implementation insights (design decisions settled by the prototype)

Key points based on measurement. Kernel is `crates/lbm-gpu-proto/src/shader.wgsl`.

1. **Memory layout is SoA (per-direction planes `f[q*n + i]`)**. Bringing over the CPU's
   cell-major AoS as-is would completely kill coalescing. Sticking to the 72 B/cell/step
   transfer floor is thanks to the combination of SoA + fused kernel.
2. **Single fused collide+stream kernel** (pull method, ping-pong double buffer).
   Intermediate write-out disappears, halving bandwidth. Shared memory/tiling turned out
   to be **unnecessary** (locality of the gather is sufficient).
3. **Wide workgroups are fastest, but the difference is small**: 256×1 = 7,550 / 128×1 =
   7,394 / 8×8 = 6,567 MLUPS (at 1024², within a ±8% band). A shape consistent with SoA
   row-direction contiguous access is good. 128×1–256×1, within the WebGPU default limit
   (256 invocations), is enough as the default.
4. **Submit granularity is the biggest pitfall** (1024², wg 256×1):

   | dispatch/submit | wait on every submit | MLUPS |
   |---|---|---|
   | 1 | **yes** | **821** (9x slower) |
   | 1 | no | 7,036 |
   | 10 | no | 7,297 |
   | 100 | no | 7,416 |

   Synchronizing completion on the CPU every step ruins everything. **"Encode N steps
   together → 1 submit → wait only when needed" should be enforced as the API shape**
   (`run(steps)` is exactly this shape).
5. **Velocity-field readback is 1.3–1.9 ms** (moments kernel + copy + map, blocking;
   1.9 ms even for the 33.6 MB of 2048²). Once per frame is within the 60 fps budget.
   However, reading every step collapses along with (4), so moments should only be
   computed/fetched on demand.
6. **2048² is a single binding of 151 MB** → exceeds the WebGPU default limit
   (max_storage_buffer_binding_size = 128 MiB). On native, requesting the adapter's
   limit is fine (Metal is GB-scale). For browser support, either request a higher
   limit or split the direction planes across bindings.
7. One uniform (nx, ny, ω+, ω−) + 2 storage faces is sufficient. Neither push
   constants nor timestamp queries were needed.

## 4. Proposed architecture for adoption

Policy: **swap the backend behind the Simulation API**. scenario/CLI/MCP/GUI stay unchanged.

```
lbm-scenario / lbm-cli / lbm-wasm / GUI
        │  (existing public API surface = extracted into a trait)
        ▼
trait LbmEngine {            // run(steps), rho/ux/uy access, set_solid, init_with, …
    // Field access is shaped to make "read on demand" explicit (either &mut self or
    // returning Cow), to allow lazy readback of GPU-resident data
}
        ├─ CpuSimulation<T=f32|f64>   … current lbm-core::Simulation (wrap unchanged)
        └─ GpuSimulation (f32)        … new crate lbm-gpu (wgpu; kept separate from core,
                                          with the option to exclude it from the default
                                          workspace build)
```

- **Staged introduction**: (1) trait extraction (only the surface used by CLI/GUI:
  run / moments / solid / init)
  → (2) implement in lbm-gpu "periodic + wall bounce-back + BGK/TRT single-phase"
  (extension of this prototype; the solid mask is a single u8 buffer, walls are just
  a reversed read at pull time)
  → (3) add Zou–He / outflow as a small kernel touching only edge cells, appended as
  one more dispatch in a later stage
  → (4) Guo forcing / Shan–Chen (force-field computation is yet another kernel;
  multiphase is phase 2).
- **Capability flags**: scenario's validate step explicitly rejects or CPU-falls-back
  "features the GPU backend doesn't yet support" (f64, some BCs, probes, etc.).
  No silent degrade.
- **Verification**: freeze this evaluation's "identical initial conditions, 2000-step L∞"
  method as a CPU↔GPU regression test, and additionally run the existing
  VALIDATION.md's main cases (TGV convergence order, Poiseuille, cavity) on the GPU
  implementation too (should come out at the same level as the measured f32 CPU results
  — consistent with the deviation-storage track record).
- **GUI/WASM**: wgpu runs as-is on WebGPU. In the browser, the real prize is **rendering
  the density/velocity buffer directly in the same device's render pass without
  readback** (zero readback, putting a 1024²-class real-time GUI within reach).
- **wgpu version**: the prototype is pinned to 26.0.1. At adoption time, move to the
  latest (30.x) — the API surface for compute usage is small, so migration cost is low.

## 5. Risks and limitations (honestly)

| Risk | Measurement/fact | Mitigation |
|---|---|---|
| **No f64** (WGSL/Metal have no f64) | With f32 deviation storage, TGV 2000 steps gives L∞ 7e-6; accuracy against the analytical solution is identical to CPU f32. PHYSICS.md's f32 track record (momentum error 2.8e-7) is also a tailwind | Keep a CPU f64 path for verification-grade computation and long-time integration (coexisting via the trait). If needed, double-single (2×f32) emulation is also an option (~2-3x cost) |
| **Readback cost** | 1.3–1.9 ms/call (~33.6 MB). A 9x slowdown (821 MLUPS) if synced every step | Enforce GPU-resident data + batched execution via the API. Visualization avoids readback entirely via in-GPU rendering |
| **Doesn't win at small grids** | At 256² and below, the CPU runs in tens of µs/step, and dispatch/sync overhead is relatively large (see submit granularity table) | Auto-select backend by grid size (decide the ~256² threshold with measurement) |
| **No bit-reproducibility across devices** | Reproducible on the same machine (this evaluation's verification values match run to run). But FMA/reassociation is compiler-dependent, so the lowest bits change on a different GPU | Regression tests are tolerance-based (same policy as existing). No bit-exact comparison against golden values |
| **Path to WASM/WebGPU** | The wgpu API is expected to work as-is, but (1) the 128 MiB binding limit (→2048² needs a limit request or splitting) (2) Safari's WebGPU maturity (3) in the browser, map_async is genuinely asynchronous | Default to up to 1024² with a limit request, and route fields directly to rendering to avoid readback. The existing CPU-WASM fallback remains available |
| **Unimplemented scope of the prototype** | Only periodic boundary + TGV (walls, open boundaries, obstacles, forces, multiphase, probes are unimplemented). BGK/TRT are implemented | The staged introduction plan in §4. Wall bounce-back is conceptually simple as a single pull-side branch (performance impact needs remeasurement; the branch should be cheap if warp-uniform) |
| **Pace of wgpu API change** | 4 major releases per year. Prototype is pinned to 26 | Upgrade to latest at adoption; thereafter following once or twice a year should suffice (the compute surface is stable) |

## 6. Verdict

**Recommend conditional adoption (yes)**. Rationale:

1. **Performance**: measured 16–42x (a straightforward result of being bandwidth-bound).
   This leaps past the CPU optimization via SoA+SIMD (Phase 9a/9b, expected 2-3x), and is
   **the only path to run 3D (Phase 10) at a realistic speed**.
2. **Accuracy**: deviation-storage f32 carries over to verification-grade as-is on GPU.
   The difference from CPU is only rounding accumulation on the order of 1e-5, and
   "same physics" has been shown by measurement.
3. **Integration cost**: the kernel body is ~150 lines of WGSL. All risks are of a kind
   that design can avoid (residency, capability flags, f64 coexists on CPU), and no
   unknown showstoppers were found.

Conditions: (a) explicitly flag f32-only via a capability flag and keep the CPU f64 path,
(b) design the API assuming batched execution and GPU residency (disallow per-step sync),
(c) once walls/open boundaries are added, always run the same CPU↔GPU regression as this
evaluation.
