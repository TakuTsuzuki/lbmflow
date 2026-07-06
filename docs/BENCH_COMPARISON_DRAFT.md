# Published Benchmark Comparison Draft (M-E headline)

**Status**: working draft (2026-07-07 refresh; ME-1 GREEN invalidates the
2026-07-05 "no 3D GPU yet" premise). Do not list numbers without a source
URL. Where a figure could not be found, note "no published figure" rather
than infer. Our own numbers are Tier-1 measured; the OpenLB CPU comparison
is a finalized single-window head-to-head, other OSS comparisons are
staged at `~/projects/cfd-bench/` and marked "measurement pending" until
their formal window closes.

Targets from COMPETITIVE_SPEC.md §1: FluidX3D / M-Star CFD / waLBerla /
Palabos / OpenLB.

---

## 1. Summary

- **3D GPU (D3Q19) on M5 Max**: LBMFlow measures **2 791 – 2 813 MLUPS
  @192³, 2 778 – 2 880 MLUPS @128³ f32** (quiet-window A/B/A,
  2026-07-06 — claims-ledger ME-1 GREEN). FluidX3D on the same M-series
  Apple GPU family reports 4 641 (M2 Max, FP16S) — a different device
  and precision, but the closest published anchor.
- **FP16 storage** delivers ~2.0× MLUPS at 2048² and D3Q19 f16 >5 GLUPS
  on the same machine (T16 bands frozen, ME-2 GREEN).
- **CPU 3D**: 302 MLUPS on M5 Max 18C (D3Q19 f32) sits in the same
  order of magnitude as 64–128-core-class server CPUs (204–330 MLUPS,
  §5).
- **M-Star publishes MLUPS as charts only** (no text values, §3).
- **Multi-node** is the domain of waLBerla / OpenLB (trillion-cell /
  TLUPS-class, §7). ME-3 cluster campaign RED, blocked on cluster
  access.

---

## 2. FluidX3D published benchmarks (single GPU, D3Q19)

Source: [FluidX3D GitHub README](https://github.com/ProjectPhysX/FluidX3D)
([raw](https://raw.githubusercontent.com/ProjectPhysX/FluidX3D/master/README.md),
retrieved 2026-07-05).

Conditions (from README): D3Q19 SRT, no extensions, empty cubic box,
typically 256³. Compute always FP32; storage FP32/FP16S/FP16C. Arithmetic
intensities 2.37 / 5.27 / 16.56 FLOPs/Byte; "performance is limited by
memory bandwidth alone." Memory per cell 55 B with Esoteric-Pull + FP16.

| GPU | FP32/FP32 | FP32/FP16S | FP32/FP16C |
|---|---:|---:|---:|
| AMD MI300X | 22 867 | **41 327** | 31 670 |
| NVIDIA H100 NVL | 20 303 | **32 922** | 18 424 |
| NVIDIA H100 SXM5 | 17 602 | **29 561** | 20 227 |
| NVIDIA RTX 5090 | 9 522 | 18 459 | **19 141** |
| NVIDIA A100 PCIe 80GB | 9 657 | **17 896** | 10 817 |
| NVIDIA A100 PCIe 40GB | 8 526 | **16 035** | 11 088 |
| NVIDIA RTX 4090 | 5 624 | 11 091 | **11 496** |
| NVIDIA RTX 3090 | 5 418 | **10 732** | 10 215 |
| Apple M2 Ultra (76-CU) | 4 629 | **8 769** | 7 972 |
| Apple M3 Ultra (60-CU) | 4 438 | **8 174** | 8 086 |
| Apple M2 Max (38-CU) | 2 405 | **4 641** | 2 444 |
| Apple M1 Max (24-CU) | 2 369 | **4 496** | 2 777 |

Notes:
- **No M4 or M5 Max/Pro/Ultra entries in the FluidX3D table** as of
  retrieval — the closest anchors on Apple silicon are M2 Max / M2
  Ultra / M3 Ultra.
- No dedicated multi-GPU benchmark table in the README. Multi-GPU
  exists; multi-node (MPI) not supported. Non-commercial license.

## 3. M-Star CFD published-figure status

**No text MLUPS values could be found; official documentation publishes
scaling charts (SVG) with stated conditions.**

| Item | Content | Source |
|---|---|---|
| Scaling benchmark | v3.3.123 results. Two case types: stirred tank (Rushton + particles, 1M–512M points) and baffled pipe (2.6M–970M points). Simulation-average MLUPS as charts for 1/2/4/8 GPU. | [Scaling Performance](https://docs.mstarcfd.com/19_Scaling_Performance/txt-files/Scaling-performance-index.html) |
| Platform | AWS p3.8xlarge (8× V100) / GCE a2-highgpu-8g (8× A100 SXM4 40GB) | same |
| Chart axis range | 8×A100: stirred-tank axis tops 16 000 MLUPS, pipe axis 25 000 MLUPS (no numeric labels; visual reading required) | [agitated SVG](https://docs.mstarcfd.com/_images/gce-a2-highgpu-8g_agitated.svg) / [pipe SVG](https://docs.mstarcfd.com/_images/gce-a2-highgpu-8g_pipe.svg) |
| Re-run availability | Benchmark package from v3.3.140+ | same |
| Sizing rule | "30–60M grid points per GPU as a guideline"; "1 GB GPU RAM ≈ 2–4M grid + 1M particles" | same / [Hardware](https://docs.mstarcfd.com/2_Installation/txt-files/hardware.html) |
| Hardware | NVIDIA-only (GeForce 40/50 through RTX 6000 Ada / H100/B100 SXM + NVLINK/NVSWITCH). "Bandwidth is what to compare next after capacity and compute." | [Hardware](https://docs.mstarcfd.com/2_Installation/txt-files/hardware.html) |

M-Star's benchmarks are full-physics ("stirred tank + particles"); do
not directly compare against empty-box kernel figures (§8.5).

## 4. Cross-code single-GPU D3Q19/Q27 (NVIDIA A100 anchor)

Each code's published values on A100. Case, grid, precision, and
streaming differ — read the conditions column.

| Code | A100 variant | Grid/case | Precision | MLUPS | Type | Source |
|---|---|---|---|---:|---|---|
| FluidX3D | PCIe 40GB | D3Q19 SRT, empty box, ~256³ | FP32/FP32 | 8 526 | measured (published table) | [README](https://github.com/ProjectPhysX/FluidX3D) |
| FluidX3D | PCIe 40GB | same | FP32/FP16S | 16 035 | measured (published table) | same |
| Palabos (GPU, C++ stdpar) | SXM4 40GB | D3Q19 BGK TGV, L=590 | FP32 | 75–85% of peak 9 481 ≈ **7 100 – 8 060** (converted) | converted from paper % | [arXiv:2506.09242](https://arxiv.org/html/2506.09242v1) |
| Palabos (GPU) | SXM4 40GB | D3Q19 BGK, L=480 | FP64 | same efficiency band of peak 4 921 | as stated | same |
| waLBerla-wind | JUWELS Booster | **D3Q27** cumulant, full solver + turbine | FP32 | 1 677 (22.3% of roofline 7 513) | measured (paper) | [arXiv:2402.13171](https://arxiv.org/html/2402.13171v1) |
| OpenLB 1.5 | 4× A100 node | D3Q19 BGK, 1000³ cavity, Periodic Shift | FP32 | 24 800/node ≈ **6 200/GPU** (converted) | measured official + conversion | [OpenLB 1.5](https://www.openlb.net/news/openlb-release-1-5-available-for-download/) |
| **LBMFlow (us)** | — (no A100 hardware) | — | — | **not measured on A100**; M5 Max D3Q19 f32 = 2 791 – 2 813 @192³ | measured on Apple hardware | claims-ledger ME-1 GREEN |

Reference (other GPUs): STLBM (Palabos-family) reports D3Q19 **FP64**
cavity N=128 at GTX 1080 Ti ≈820 / RTX 2080 Ti ≈1 100 / V100 PCIe
≈2 300 MLUPS (AA-pattern) —
[PLOS ONE](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0250306).

## 5. CPU (single node) 3D comparison

| Code | CPU | Grid/case | Precision | MLUPS | Source |
|---|---|---|---|---:|---|
| **LBMFlow (us)** | Apple M5 Max 18C | D3Q19 192³ | f32 | **302** | measured 2026-07-05, `benchmark-results.md` |
| **LBMFlow (us)** | Apple M5 Max 18C | D3Q19 128³ | f32 | **267** | same |
| **LBMFlow (us)** | Apple M5 Max, 1T | D3Q19 128³ | f32 | **52.0** | same window (single-thread win) |
| OpenLB 1.9 | Apple M5 Max 18C (CPU_SISD; ARM has no AVX/NEON path OOTB) | D3Q19 128³ | f32 | 298.8 (18r) / 44.6 (1r) | measured same window, `benchmark-results.md` |
| STLBM | AMD EPYC 64-core | D3Q19 cavity N=128 | FP64 | ≈300 (AA, SoA) | [PLOS ONE](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0250306) |
| STLBM | Intel Xeon 48-core | same | FP64 | ≈330 (swap-AoS) | same |
| waLBerla-wind | AMD EPYC 7763 128-core (1 node) | **D3Q27** cumulant (no turbine) | FP32 | 204 (roofline 461) | [arXiv:2402.13171](https://arxiv.org/html/2402.13171v1) |
| OpenLB 1.3 (2018) | Magnus 32 784-core | D3Q19 | — | total 142 479 (≈4.3 MLUPS/core) | [openlb.net/performance](https://www.openlb.net/performance/) |

Head-to-head reading (same machine + window vs OpenLB 1.9 128³ f32):
LBMFlow wins single-thread **+17%** (52.0 vs 44.6, NEON auto-vec vs
OpenLB CPU_SISD on ARM); OpenLB wins all-18-core **+12%** (298.8 vs
266.6) via MPI domain decomposition. STLBM/waLBerla comparisons need
precision and stencil discounts before "equivalent" — competitive not
equivalent (§8).

## 6. 2D (D2Q9) — separate category

| Code | Device | Grid | Precision | MLUPS | Source |
|---|---|---|---|---:|---|
| **LBMFlow (us)** | M5 Max GPU (Metal, wgpu) | D2Q9 512² / 1024² / 2048² | f32 | **12 205 / 7 073 / 6 720** | measured, `benchmark-results.md` |
| **LBMFlow (us)** | M5 Max 18C CPU | D2Q9 2048² | f32 | **1 480** (peak) | same |
| **LBMFlow (us)** | M5 Max 18C CPU | D2Q9 1024² | f32 | 1 208 | same |

**No published 2D figures from competitors were found** — FluidX3D
3D-only, M-Star 3D product, waLBerla/OpenLB/Palabos representative
benchmarks 3D. D2Q9 moves ~half the data per update of D3Q19; placing
D2Q9 MLUPS in the same table as D3Q19 inflates by roughly 2×. Keep
separate externally.

## 7. Multi-node / scaling published track record

| Code | Machine / scale | Track record | Source |
|---|---|---|---|
| waLBerla | JUQUEEN 458 752 cores / 1.8M threads | >1 trillion cells, up to 1.93 trillion cell updates/s (abstract; visual confirm needed). Strong scaling on SuperMUC 32 768 cores | [SC13 DOI:10.1145/2503210.2503273](https://dl.acm.org/doi/10.1145/2503210.2503273) |
| waLBerla | JUQUEEN | ">1 trillion cells; good scalability to >400 000 cores" (direct body-text quote) | [arXiv:1511.07261](https://ar5iv.labs.arxiv.org/html/1511.07261) |
| waLBerla-wind | JUWELS Booster 30 nodes / 120 A100 | Per-GPU perf ≈constant under weak scaling (17.5M cells/GPU, ~74.46 steps/s ≈ ~1 300 MLUPS/GPU) | [arXiv:2402.13171](https://arxiv.org/html/2402.13171v1) |
| OpenLB 1.5 | HoreKa 128 nodes / 512 A100 | ≈1.33 TLUPS (D3Q19 FP32 cavity). 64→128 GPU strong-scaling eff 0.64–0.81 (575³–2300³). 92% of bench perf on LES turbulent-nozzle case (224 GPU) | [openlb.net/performance](https://www.openlb.net/performance/) |
| OpenLB 1.5 | HoreKa 2 nodes / 8 A100 | 1000³ FP32 cavity 42.2 GLUPS (1-node 4×A100 24.8 GLUPS, 2-node CPU AVX-512 2.7 GLUPS → GPU 15.6×) | [OpenLB 1.5](https://www.openlb.net/news/openlb-release-1-5-available-for-download/) |
| OpenLB 1.9 | Aurora ~1 000 nodes (~10% of system) | Peak 21 120 GLUPS, 4 trillion cells (D3Q19 FP32) | [openlb.net/performance](https://www.openlb.net/performance/) |
| Palabos GPU | DGX 4× A100 40GB | 80–90% weak-scaling ideal, 65–80% strong | [arXiv:2506.09242](https://arxiv.org/html/2506.09242v1) |
| M-Star | 8× V100 / 8× A100 (single node) | 1/2/4/8 GPU scaling charts (values in figures, §3) | [Scaling Performance](https://docs.mstarcfd.com/19_Scaling_Performance/txt-files/Scaling-performance-index.html) |
| **LBMFlow (us)** | — | **measurement pending** (ME-3). MPI 3D bench + weak modes ready in code; blocked on cluster spend. Acceptance = 64-rank weak ≥80% | claims-ledger ME-3 RED |

## 8. Reading rules for the tables

1. **MLUPS definition varies**. Common: "million lattice-point updates
   per second." What is included (boundary, output, communication)
   varies by code.
2. **Stencil discount**. D3Q27 moves more per update than D3Q19
   (27 vs 19 DDFs); under bandwidth-bound conditions MLUPS is
   structurally lower. Reading waLBerla-wind 1 677 (D3Q27) against
   FluidX3D 8 526 (D3Q19) as "waLBerla is slow" is a mistake. D2Q9
   is a separate category (§6).
3. **Precision discount**. Halving storage precision roughly doubles
   MLUPS in bandwidth-bound LBM (FluidX3D RTX 4090 FP32 5 624 → FP16S
   11 091 — [README](https://github.com/ProjectPhysX/FluidX3D);
   Palabos FP32 peak 9 481 vs FP64 4 921 GLUPS —
   [arXiv:2506.09242](https://arxiv.org/html/2506.09242v1)). Lehmann
   et al.: FP64 vs FP32 negligible in almost all cases, 16-bit
   sufficient in many
   ([arXiv:2112.08926](https://arxiv.org/abs/2112.08926)). Always
   label "compute / storage" together.
4. **Grid-size dependence**. Small grids underperform. FluidX3D uses
   "typically 256³"; Palabos GPU paper states perf increases with
   resolution; M-Star recommends 30–60M cells/GPU. Report saturation
   and the size at which it is reached.
5. **Kernel-alone vs full solver**. Same code, same GPU: empty box
   vs full physics differs 4–5× (waLBerla-wind roofline 7 513 →
   measured 1 677 with turbine = 22.3%; OpenLB counterexample: 92%
   of bench perf on real LES). FluidX3D's table = empty box, M-Star
   chart = stirred tank + particles. **Do not compare directly.**
6. **Loaded-window / warm-up**. On unified memory, background CPU
   load halved GPU MLUPS in this repo (D3Q19 192³: 1 353 loaded vs
   2 791 – 2 813 quiet, PERFORMANCE.md). Publish quiet-window,
   warm-up-excluded, best-of-N + variance. M-Star publishes
   "simulation-average" — a different convention worth naming.
7. **Implementation variance**. On identical hardware, data layout
   and streaming implementation change perf by 2–3×
   ([arXiv:1711.11468](https://arxiv.org/abs/1711.11468); concrete
   example: 211 → 550 MLUPS
   [SAGE](https://journals.sagepub.com/doi/full/10.1177/10943420211016525)).
   Argue "code vs code" separately from hardware and conditions.

## 9. Where LBMFlow is losing today

1. **FP16C** (compressed 16-bit) not implemented. FluidX3D delivers
   another 30% headroom with FP16C on some GPUs. Awaiting P4.
2. **NVIDIA/AMD flagship headroom is unreachable on Apple silicon.**
   H100 NVL 32 922, MI300X 41 327 (FP16S) require a CUDA/HIP backend
   and access to real hardware (SPEC §5).
3. **Multi-node**: zero track record. waLBerla trillions of cells /
   400k+ cores; OpenLB 1.33 TLUPS / 4 trillion cells on Aurora. Our
   R3 target (64-rank weak ≥80%) is an initial step orders of
   magnitude below these.
4. **Full-physics stirred workload not benchmarked.** M-Star's chart
   is a "sellable workload." We know the empty-box + basic-boundary
   figure; the LES/moving/scalar layering drop is not yet measured
   (ME-4, tracks M-F).

## 10. Where LBMFlow leads today (measured, not marketing)

- **3D GPU on Apple silicon**: 2 791 – 2 813 MLUPS D3Q19 f32 192³
  quiet-window (this repo, ME-1 GREEN). OpenLB/Palabos GPU backends
  are CUDA-only — no comparable path on this hardware family.
- **Single-thread CPU 3D**: LBMFlow +17% vs OpenLB 1.9 on same
  machine (128³ f32: 52.0 vs 44.6 MLUPS), via NEON auto-vec.
- **Bundled verification suite, one-command re-run**: 56+
  adversarial tests (Ghia / Schäfer-Turek / RT / equivariance 4e-16).
  Competitors' verification is mostly papers or published benchmark
  collections; re-runnable physical-accuracy verification as part
  of the product spec is our differentiator. (M-Star ships a
  scaling-benchmark package from v3.3.140+, so re-running
  *performance* is partially covered there.)
- **Agent-native**: JSON Schema self-description + MCP. M-Star's
  Python API is human-oriented; FluidX3D requires C++ setup edit.
- **Portability**: wgpu (Metal/Vulkan/DX12) + WASM.
  M-Star assumes NVIDIA. FluidX3D uses OpenCL (all-vendor) but is
  non-commercial and no MPI. Commercially usable + all-vendor +
  browser is ours alone.

## 11. Measurement pending (live workstreams)

- **OpenLB 2D (cavity2d)** for the D2Q9 head-to-head vs 1 480 MLUPS.
- **OpenLB CPU_SIMD on ARM** (if buildable) — SIMD-vs-SIMD line.
- **Palabos** — cmake install + build at `~/projects/cfd-bench/`.
- **OpenFOAM via Colima** — cylinder Re=20 Cd time-to-solution
  (Family B).
- **LBMFlow FluidX3D same-machine head-to-head** (FluidX3D OpenCL
  builds on Apple silicon; the FluidX3D published Apple values are
  M2 Max / M2 Ultra / M3 Ultra, not M5).
- **A100 head-to-head** — pending cluster access; the LBMFlow row
  in §4 stays "not measured on A100" until then.
- **ME-3 multi-node weak scaling** — bench_mpi ready; blocked on
  cluster spend.
- **ME-4 full-physics stirred workload** — waits on M-F integration.

## 12. Source list

| # | Source | Use |
|---|---|---|
| 1 | https://github.com/ProjectPhysX/FluidX3D | FluidX3D table, conditions, intensity, memory/cell |
| 2 | https://raw.githubusercontent.com/ProjectPhysX/FluidX3D/master/README.md | same (raw, 2026-07-05) |
| 3 | https://docs.mstarcfd.com/19_Scaling_Performance/txt-files/Scaling-performance-index.html | M-Star scaling conditions |
| 4 | https://docs.mstarcfd.com/_images/gce-a2-highgpu-8g_agitated.svg / …_pipe.svg | M-Star chart originals |
| 5 | https://docs.mstarcfd.com/2_Installation/txt-files/hardware.html | M-Star NVIDIA requirements |
| 6 | https://arxiv.org/html/2402.13171v1 | waLBerla-wind A100 D3Q27 measured / roofline, EPYC 7763, weak scaling |
| 7 | https://ar5iv.labs.arxiv.org/html/1511.07261 | waLBerla trillion cells / 400k+ cores, MLUP/s definition |
| 8 | https://dl.acm.org/doi/10.1145/2503210.2503273 | waLBerla SC13 (1.93 trillion updates/s from abstract) |
| 9 | https://www.openlb.net/performance/ | OpenLB 512×A100 1.33 TLUPS, Aurora 21 120 GLUPS, Magnus 142 479 MLUPS |
| 10 | https://www.openlb.net/news/openlb-release-1-5-available-for-download/ | OpenLB 1000³ cavity 42.2/24.8/2.7 GLUPS |
| 11 | https://arxiv.org/html/2506.21804v1 | OpenLB heterogeneous HPC, 18G cells, eff 0.66–0.91 |
| 12 | https://arxiv.org/html/2506.09242v1 | Palabos GPU: A100 peak 9 481/4 921 GLUPS, measured 75–85% |
| 13 | https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0250306 | STLBM: GPU/CPU D3Q19 FP64 MLUPS |
| 14 | https://arxiv.org/abs/2112.08926 | FP64/FP32/16-bit precision (Phys. Rev. E 106, 015308) |
| 15 | https://arxiv.org/abs/1711.11468 | LBM benchmark kernels (implementation impact) |
| 16 | https://journals.sagepub.com/doi/full/10.1177/10943420211016525 | Optimization range 211→550 MLUPS |
| 17 | `docs/paper/benchmark-results.md` + `docs/paper/claims-ledger.md` + `TESTING_NOTES.md` | our measured figures (Tier-1 quiet-window) |
