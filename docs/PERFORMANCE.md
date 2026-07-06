# PERFORMANCE.md — Measured Performance and Precision/Speed Trade-off Guide

Machine: Apple M5 Max, 18 cores, 128 GB, macOS. `--release` (thin LTO,
codegen-units=1). MLUPS = million lattice-point updates per second;
GLUPS = billion. Definition = total cells × steps / wall time, warm-up
excluded via best-of-N. Measurement path: `examples/bench_backends.rs`
(CPU) / `examples/bench_gpu.rs` (GPU) / T14-3D quiet-window A/B/A for the
current 3D GPU headline. All figures below are on the same machine.

## Current measured headline (M-E, 2026-07-07)

| backend | grid | precision | MLUPS / GLUPS | notes |
|---|---|---|---:|---|
| GPU (Metal, wgpu) D3Q19 f32 | 192³ | f32/f32 | **2 791 – 2 813** MLUPS | quiet-window A/B/A, 2026-07-06 |
| GPU (Metal, wgpu) D3Q19 f32 | 128³ | f32/f32 | **2 778 – 2 880** MLUPS | same window |
| GPU (Metal, wgpu) D3Q19 f16 | 2048² equiv | f32/f16 storage | **> 5 GLUPS**, ~2.0× f32 | ME-2 GREEN, T16 bands frozen |
| GPU (Metal, wgpu) D2Q9 f32 | 512² / 1024² / 2048² | f32/f32 | 12 205 / 7 073 / 6 720 MLUPS | measurement window #1 |
| CPU (CpuSimd) D2Q9 f32 | 2048², 18T | f32 | **1 480** MLUPS | 2D CPU peak |
| CPU (CpuSimd) D2Q9 f32 | 1024², 18T | f32 | 1 208 MLUPS | |
| CPU (CpuSimd) D3Q19 f32 | 192³, 18T | f32 | 302 MLUPS | 3D CPU peak |
| CPU (CpuSimd) D3Q19 f32 | 128³, 18T | f32 | 267 MLUPS | |

Roofline context (arm64 STREAM triad, 18 threads): 344 GB/s (STREAM
convention) / 459–464 GB/s (write-allocate) / ~546 GB/s (spec).
Roofline MLUPS ceilings (write-alloc 459 GB/s, f32 pull traffic
72 B/cell·step for D2Q9, 152 B/cell·step for D3Q19): 2D ≈ 6 375,
3D ≈ 3 020. The GPU D3Q19 192³ headline sits at ~93% of the 3D roofline.

## Loaded-window trap (live gotcha)

Unified-memory GPU MLUPS collapses under concurrent CPU load. Same D3Q19
192³ scenario measured **1 353 MLUPS under load vs 2 791 – 2 813 quiet**
(2026-07-06 — cost several hours of ME-1 false-RED before the artifact
was diagnosed). Rules that fall out of this:

- Never flip a perf gate RED, and never dispatch a kernel-optimization
  order, from a loaded-window number.
- Re-measure with A/B/A interleave in a quiet window before acting.
- Documented in the whitepaper memory trap list; the earlier
  `step_periodic` follow-up kernel measured SLOWER under the correct
  quiet-window baseline and was rejected on that ground alone.

## Mode-selection guide

Collision:
- **TRT** is the default. Same speed as BGK in the paired form (BGK is
  the ω+ = ω− special case), superior in accuracy and stability.
- **BGK** for teaching, comparison, and reproducing legacy results.
- **Cumulant / central-moment (MF-α stage 3)** for high-Re turbulent
  regimes where TRT stability margins shrink; costlier per step.

Precision (storage / compute):
- **f32** default for 2D and typical 3D. With deviation storage (f−w),
  f32 is verification-grade: momentum-growth error under uniform force
  1.34e-3 → 2.8e-7 (~4800×); TGV L2 (N=64) f32 7.1e-4 vs f64 7.0e-4.
- **f64** for long-time integration, extremely small gradients, or
  reference computations. ~1.9× slower than f32 in bandwidth-bound
  regimes.
- **f16 storage / f32 compute** (ME-2) for large 3D on-device. ~2.0×
  MLUPS at 2048²; T16 bands frozen (TGV transient 1.401e-1 vs band
  2e-1; cavity steady 2.579e-3 vs band 5e-3). Grid capacity doubles
  since each DDF is halved on-device.

Backend:
- **CpuSimd** on the CPU path. Auto-vectorized NEON; wins single-thread
  vs OpenLB on Apple Silicon (see `docs/paper/benchmark-results.md`).
- **Wgpu (Metal / Vulkan / DX12)** for GPU-covered core workloads. Prefer
  it for large f32 runs when the requested lattice/features are supported:
  2 800 MLUPS D3Q19 f32 on M5 Max is ~10× the CPU 3D peak. The
  scenario-runner `compute.backend:"gpu"` product path is narrower than the
  in-core backend: it is feature-gated and currently dispatches through the
  2D D2Q9 f32 builder; 3D GPU scenario dispatch is not wired yet.
- **MPI** (feature `mpi`, off by default) for multi-node. Weak-scaling
  campaign RED pending cluster access.

## Kernel-shape notes (unchanged since M-E port)

- **D2Q9 = flat form**. LLVM fully unrolls the 4 pairs and vectorizes
  the x loop. Blocking loses ~1.2×.
- **D3Q19 = blocked form**. For flat unroll of 9 pairs LLVM gives up
  and drops to scalar (measured vec/scalar instruction ratio 18 vs
  285). Shared quantities go in a `[T; 64]` scratch buffer; each pair
  swept with unit stride.
- **In-place views**: passing separate src/dst views into the blocked
  kernel breaks vectorization via alias checking (−30%). Blocked
  kernel uses a single view; flat kernel uses copy-fused src=f/dst=ring.
- **Ring persistence**: per-band ring carried in `FusedScratch`
  (eliminates 92 MB/step@f64 per-step allocation + zero-fill).
- One step = collide → halo exchange → streaming → open-boundary BCs
  → boundary moments correction (CpuSimd fuses collide+stream+moments
  in `step_band`). Any change to pass structure or storage order must
  pass `tests/backend_simd_equiv.rs` and T13 (partition invariance).

## Known optimization opportunities

- **Explicit SIMD** (`std::simd` / `wide`) on the CpuSimd path: ~2×
  theoretical headroom remains at f32 1T. Currently auto-vectorization
  only.
- **In-place streaming (Esoteric-pull)**: halves memory, roughly
  doubles effective bandwidth by eliminating `ftmp`. Candidate for
  the GPU path.
- **3D 12T band scaling**: 12T f32 falls short of 2× CpuScalar
  (1.9×) at 128³ because of (a) double collision at band-edge slabs
  (+19% collision work with 12 bands) and (b) coarse-grained bands
  can't load-balance across heterogeneous P/E cores. Both mitigations
  (y-strip ringing, over-splitting bands) were measured and rejected
  on this machine.
- **CpuSimd facade switch** (2D compat path): still on hold pending
  the remaining synchronization-point/facade contract work.
