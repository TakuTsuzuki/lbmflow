# GPU_EVALUATION.md — kernel-shape lessons from the wgpu prototype

**Prototype banner: superseded by `crates/lbm-core` gpu module — historical record.**

**Status: historical evaluation, superseded.** The 2026-07-05 throwaway
prototype in `crates/lbm-gpu-proto` motivated the in-core wgpu backend
(`crates/lbm-core/src/gpu`, landed and GREEN by 2026-07-06 — see M-E in
docs/PLAN.md and the claims ledger). Current headline numbers live in
[PERFORMANCE.md](PERFORMANCE.md); do not re-quote 2D-only prototype
figures as product claims. This file is retained only for the kernel-shape
design decisions the prototype pinned down, plus the risk list they informed.

## Kernel-shape decisions (still current in the in-core backend)

1. **SoA per-direction planes `f[q*n + i]`** — cell-major AoS kills coalescing
   and blows past the 72 B/cell/step (D2Q9 f32) memory-traffic floor.
2. **Single fused collide+stream kernel** (pull, ping-pong double buffer).
   Intermediate write-out disappears; halves bandwidth. Shared-memory tiling
   turned out to be unnecessary (gather locality is sufficient).
3. **Wide 1D workgroups (128×1–256×1)** are consistent with SoA
   row-contiguous access and fastest; the difference to 8×8 was ~15%.
   256-invocation default sits within the WebGPU limit.
4. **Batched submit is mandatory.** Per-step CPU sync at 1024² collapsed
   throughput 9× (821 vs 7,036 MLUPS). The `run(steps)` API contract must
   encode N steps → 1 submit → wait only when moments are read.
5. **Velocity-field readback is 1.3–1.9 ms** (moments kernel + copy + map,
   blocking). Fine at 60 fps; catastrophic if driven every step. Moments
   are computed/fetched on demand.
6. **2048² D2Q9 f32 is a 151 MB single binding** — over the WebGPU
   default 128 MiB. Native adapters (Metal is GB-scale) request the higher
   limit; for browser paths, split direction planes across bindings.
7. **f32 with deviation storage (f−w) is verification-grade.** TGV
   2000 steps: L∞(GPU vs CPU) 7.0e-6, 14× under the 1e-4 gate;
   accuracy vs analytical matches CPU f32 to the 4th significant digit.
   The stationary background is exactly zero, so rounding only hits the
   fluctuation scale.
8. **Operator-ordering identity for CPU↔GPU regressions.** Fused pull kernel
   is `C∘S` per step vs CPU's `S∘C`. Applying one collision on the host to
   the initial state before upload gives `(C∘S)^k ∘ C = C ∘ (S∘C)^k`;
   density/momentum are collision-invariant so velocity fields compare 1:1.

## Probe-force CAS accumulation (2026-07-07)

The in-core WGSL path keeps probe force on-device. During bounce-back into a
solid cell marked as a probe, each streamed link contributes its
momentum-exchange force component to the three-element `probe_acc` buffer
(`Fx`, `Fy`, `Fz`). The value added per link is the physical incoming plus
reflected population in deviation-storage form: the two deviations plus the
two equilibrium weights, signed by the lattice direction component. The
`clear_probe` kernel zeroes this accumulator at the start of a probed step, so
the buffer represents the most recent step until an explicit readback.

WGSL/WebGPU does not provide portable floating-point atomics, so
`probe_acc` is stored as `array<atomic<u32>, 3>` and the shader implements
`atomic_add_f32` with a compare-and-swap loop: load the old bits, bitcast to
`f32`, add the contribution, bitcast the new value back to `u32`, and retry
with `atomicCompareExchangeWeak` until the exchange succeeds. This preserves
device residency and avoids the per-step CPU synchronization that the prototype
identified as a 9× throughput trap.

The result is a diagnostic, not a bit-exact invariant. Contributions from
different invocations arrive in nondeterministic order, so the accumulated
`f32` sum has normal order-dependent roundoff. CPU↔GPU force checks therefore
use tolerances, matching the broader GPU regression policy.

## Risks the prototype flagged (still live)

| Risk | Mitigation in the current backend |
|---|---|
| No f64 in WGSL/Metal | f32 with deviation storage is default; CPU f64 path retained for reference/long-time integration. |
| Readback cost | GPU-resident data + batched execution; visualization renders in-device. |
| Doesn't win at small grids | Auto-select backend by grid size; ~256² CPU is faster than GPU dispatch overhead. |
| No bit-reproducibility across devices | Regression tests are tolerance-based; no bit-exact goldens against GPU output. |
| WebGPU browser path | 128 MiB binding cap → limit request or plane split; Safari WebGPU maturity is a moving target. |
| wgpu API cadence (~4 major/yr) | Prototype was pinned to 26.0.1; core has since tracked forward. Compute surface is stable. |

For the current in-core backend scope (generic `WgpuBackend<L>` over the core
lattice tables, exercised by D3Q19 and D3Q27 GPU tests; f32 compute with f32
or f16 distribution storage; WALE and cumulant kernel entries; monolithic
decomposition) see `crates/lbm-core/src/gpu/` and the ME-1/ME-2/C-13 rows in
`docs/paper/claims-ledger.md`. The scenario-runner product path is narrower:
`compute.backend:"gpu"` is feature-gated and currently dispatches through the
2D D2Q9 f32 builder; 3D GPU scenario dispatch is not wired yet, and unsupported
GPU combinations are rejected instead of silently falling back.
