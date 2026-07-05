# Published Benchmark Comparison Draft (M-E: for the performance headline)

**Status**: draft (results of web-published information collection, retrieved 2026-07-05)
**Rule**: do not list numbers without a source URL. Where a quantitative figure could not be found,
explicitly note "no published figure / unconfirmed." Our own numbers are provisional values under shared
load (before idle re-measurement). **External publication only after the §8 checklist is complete**.

Targets are the competitors from COMPETITIVE_SPEC.md §1: FluidX3D / M-Star CFD / waLBerla / Palabos / OpenLB.

---

## 1. Summary (with reading caveats)

- **The published single-GPU 3D (D3Q19) champion is FluidX3D**. On an A100 PCIe 40GB: 8,526 MLUPS with
  FP32 storage, 16,035 MLUPS with FP16S storage. The high end (MI300X 41,327 / H100 NVL 32,922, both
  FP16S) is an order of magnitude beyond that.
- **We have no 3D GPU implementation yet**, so we can't compete on this turf. Our 2D D2Q9 GPU figures of
  5,857–11,365 MLUPS are a different event — data movement is roughly half that of D3Q19 — and **must not
  be mixed into the 3D table** (§6).
- **CPU 3D holds up well**: 260 MLUPS on an M5 Max 18C is in the same order of magnitude (with caveats on
  conditions, §4) as published measurements on 64–128-core-class server CPUs (204–330 MLUPS).
- **M-Star has no published MLUPS figure as text** (charts only, §3).
- **Multi-node is the domain of waLBerla / OpenLB** (trillion-cell / TLUPS-class, §5). We have zero track
  record here since M-D is not yet reached.

---

## 2. FluidX3D published benchmarks (single GPU, D3Q19)

Source: [FluidX3D GitHub README](https://github.com/ProjectPhysX/FluidX3D)
(raw file: [raw README.md](https://raw.githubusercontent.com/ProjectPhysX/FluidX3D/master/README.md),
retrieved 2026-07-05)

Measurement conditions (summary of the README's original wording): "D3Q19 SRT, no extensions (pure LBM
with only implicit mid-grid bounce-back boundaries), empty cubic box, sufficient size (typically 256³)."
Compute is always FP32, with storage precision switchable between FP32 / FP16S / FP16C. The README states
arithmetic intensities of 2.37 (FP32/FP32), 5.27 (FP32/FP16S), and 16.56 (FP32/FP16C) FLOPs/Byte, and
states that "performance is limited by memory bandwidth alone." Memory usage with Esoteric-Pull + FP16
compression is 55 Bytes/cell (versus roughly 344 Bytes/cell for conventional FP64 LBM).

| GPU | FP32/FP32 | FP32/FP16S | FP32/FP16C | Source |
|---|---:|---:|---:|---|
| AMD MI300X | 22,867 | **41,327** | 31,670 | [README](https://github.com/ProjectPhysX/FluidX3D) |
| NVIDIA H100 NVL | 20,303 | **32,922** | 18,424 | same |
| NVIDIA H100 SXM5 | 17,602 | **29,561** | 20,227 | same |
| NVIDIA RTX 5090 | 9,522 | 18,459 | **19,141** | same |
| NVIDIA A100 PCIe 80GB | 9,657 | **17,896** | 10,817 | same |
| NVIDIA A100 PCIe 40GB | 8,526 | **16,035** | 11,088 | same |
| NVIDIA RTX 4090 | 5,624 | 11,091 | **11,496** | same |
| NVIDIA RTX 3090 | 5,418 | **10,732** | 10,215 | same |
| Apple M2 Ultra (76-CU) | 4,629 | **8,769** | 7,972 | same |
| Apple M3 Ultra (60-CU) | 4,438 | **8,174** | 8,086 | same |
| Apple M1 Ultra (64-CU) | (fastest-mode value 8,418) | | | same |
| Apple M2 Max (38-CU) | 2,405 | **4,641** | 2,444 | same |
| Apple M1 Max (24-CU) | 2,369 | **4,496** | 2,777 | same |
| Apple M5 (10-CU) | (fastest-mode value 1,613) | | | same |

Notes (from the README, as of 2026-07-05):
- A100 fastest-mode values differ by variant: SXM4 80GB 18,448 / PCIe 80GB 17,896 / PCIe 40GB 16,035 /
  SXM4 40GB 16,013.
- **There are no entries for the Apple M4 generation or M5 Max/Pro/Ultra in the table** (i.e., there is
  currently no published comparison figure for a machine identical to our M5 Max; the closest published
  values are M2 Max / M2·M3 Ultra).
- No dedicated multi-GPU benchmark table appears to exist in the README (as of this retrieval).
  Multi-GPU support itself exists, but multi-node (MPI) is not supported. Licensing is free for
  non-commercial use.

## 3. M-Star CFD (primary comparison target): status of published figures

**Conclusion: no MLUPS value published as text could be found (treated as "no published figure").
However, the official documentation does have scaling charts (SVG) with stated conditions.**

Confirmed facts (all from official sources):

| Item | Content | Source |
|---|---|---|
| Scaling benchmark | Results from v3.3.123. Two case types: stirred tank (1 Rushton impeller + particles, grid 1M–512M points) / baffled pipe (static geometry, 2.6M–970M points). Publishes "simulation-average MLUPS" as charts for 1/2/4/8 GPU configurations | [Scaling Performance](https://docs.mstarcfd.com/19_Scaling_Performance/txt-files/Scaling-performance-index.html) |
| Measurement platform | AWS p3.8xlarge (8× V100 SXM2 16GB, CUDA 11.5, restricted peer access) / GCE a2-highgpu-8g (8× A100 SXM4 40GB, CUDA 11.2, full peer access) | same |
| Chart axis range | For the 8×A100 configuration, the chart's vertical-axis scale tops out at 16,000 MLUPS for the stirred-tank case and 25,000 MLUPS for the pipe case (data points have no numeric labels; exact values require visual reading from the figure) | [agitated chart SVG](https://docs.mstarcfd.com/_images/gce-a2-highgpu-8g_agitated.svg) / [pipe chart SVG](https://docs.mstarcfd.com/_images/gce-a2-highgpu-8g_pipe.svg) |
| Re-run availability | A benchmark package is available from v3.3.140+ (customers can re-run it in their own environment) | [Scaling Performance](https://docs.mstarcfd.com/19_Scaling_Performance/txt-files/Scaling-performance-index.html) |
| Sizing rule of thumb | "30–60M grid points or more per GPU as a guideline" / "1GB of GPU RAM ≈ 2–4M grid points + 1M particles" | same / [Hardware Guide](https://docs.mstarcfd.com/2_Installation/txt-files/hardware.html) |
| Hardware requirements | Assumes NVIDIA GPU (recommended: GeForce 40/50 series through RTX 6000 Ada / RTX PRO 6000 through H100/B100 SXM + NVLINK/NVSWITCH). "Bandwidth is the spec to compare next, after memory capacity and compute performance" | [Hardware Guide](https://docs.mstarcfd.com/2_Installation/txt-files/hardware.html) |
| Marketing claims | "Built to run on millions and billions of lattice grid points" / "Detailed and accurate process simulation in minutes" | [mstarcfd.com/software](https://mstarcfd.com/software/) |

Note: we also attempted automated data-point extraction from the SVGs this time, but the coordinate
readings were unreliable (e.g., readings exceeding the axis maximum appeared), so **this draft does not
list data-point values**. If exact values are needed, they should be determined separately via visual
reading of the figures per §8.
Comparison implication: M-Star's benchmarks are **full-physics cases** such as "stirred tank + particles,"
and must not be compared directly against FluidX3D's "empty-box kernel" figures (§6.5).

## 4. Cross-code comparison table: single-GPU D3Q19/Q27 on the same device (NVIDIA A100)

A table aligning each code's published values on the A100. **Even so, case, grid, precision, and
streaming implementation differ, so this is not "same conditions."** Be sure to read the conditions
column.

| Code | A100 variant | Grid/case | Precision (compute/storage) | MLUPS | Type | Source |
|---|---|---|---|---:|---|---|
| FluidX3D | PCIe 40GB | D3Q19 SRT, empty box, typical 256³ | FP32/FP32 | 8,526 | measured (published table) | [README](https://github.com/ProjectPhysX/FluidX3D) |
| FluidX3D | PCIe 40GB | same | FP32/FP16S | 16,035 | measured (published table) | same |
| Palabos (GPU port, C++ stdpar) | SXM4 40GB | D3Q19 BGK, Taylor-Green, L=590 | FP32 | 75–85% of theoretical peak 9,481 ≈ **7,100–8,060** (converted) | converted from paper-stated % | [arXiv:2506.09242](https://arxiv.org/html/2506.09242v1) |
| Palabos (GPU port) | SXM4 40GB | D3Q19 BGK, L=480 | FP64 | same efficiency band of theoretical peak 4,921 | as stated in paper | same |
| waLBerla (waLBerla-wind) | JUWELS Booster | **D3Q27** cumulant, full solver with turbine | FP32 | 1,677 (22.3% of roofline ceiling 7,513) | measured (paper) | [arXiv:2402.13171](https://arxiv.org/html/2402.13171v1) |
| OpenLB 1.5 | 4× A100 node | D3Q19 BGK, 1000³ cavity, Periodic Shift | FP32 | 24,800 per node → **≈6,200 per GPU** (converted) | measured (official) + conversion | [OpenLB 1.5 release](https://www.openlb.net/news/openlb-release-1-5-available-for-download/) |
| **LBMFlow (us)** | — (no A100 measurement environment) | — | — | **not measured** (no 3D GPU implementation + no physical NVIDIA hardware) | — | COMPETITIVE_SPEC.md §5 |

Reference (cross-code points on other GPUs): STLBM (a Palabos-family research code) reports D3Q19
**FP64** cavity N=128 at GTX 1080 Ti ≈820 / RTX 2080 Ti ≈1,100 / V100 PCIe ≈2,300 MLUPS (AA-pattern,
approximate values derived from the paper's figure)
— [PLOS ONE 10.1371/journal.pone.0250306](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0250306).
The Palabos GPU paper also states, as related work, that "the waLBerla CUDA backend achieves about 85% of
theoretical peak on an A100-SXM4 40GB" — [arXiv:2506.09242](https://arxiv.org/html/2506.09242v1).

## 5. CPU (single node) 3D comparison

| Code | CPU | Grid/case | Precision | MLUPS | Source |
|---|---|---|---|---:|---|
| **LBMFlow (us)** | Apple M5 Max 18C | D3Q19 | f32 | **260** (under shared load, provisional) | measured in this repository (PERFORMANCE.md series) |
| STLBM | AMD EPYC 64-core | D3Q19 cavity N=128 | FP64 | ≈300 (AA-pattern, SoA) | [PLOS ONE](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0250306) |
| STLBM | Intel Xeon 48-core | same | FP64 | ≈330 (swap-AoS) | same |
| waLBerla-wind | AMD EPYC 7763 128-core (1 node) | **D3Q27** cumulant (no turbine) | FP32 | 204 (roofline ceiling 461) | [arXiv:2402.13171](https://arxiv.org/html/2402.13171v1) |
| OpenLB 1.3 (2018) | Magnus 32,784-core | D3Q19 | see original text | total 142,479 (≈4.3 MLUPS/core) | [openlb.net/performance](https://www.openlb.net/performance/) |

Honest reading: being in the same order of magnitude as published values for 64–128-core-class server
CPUs, on an 18-core laptop SoC, is a strong data point. However, STLBM is FP64 (roughly double the
bandwidth consumption → equivalent to roughly half the MLUPS), and waLBerla-wind is D3Q27 (data movement
≈27/19×), so **once corrected for precision and stencil, this is not "equivalent" but rather "fairly
viewed as competitive, depending on conditions."** A firm claim will be made after formal measurement.

## 6. 2D (D2Q9): where to place our GPU figures

| Code | Device | Grid | Precision | MLUPS | Source |
|---|---|---|---|---:|---|
| **LBMFlow (us)** | M5 Max GPU (Metal, wgpu) | D2Q9 | f32 | **5,857–11,365** (under shared load, provisional, range depends on conditions) | measured in this repository |
| **LBMFlow (us)** | M5 Max 18C CPU | D2Q9 | f32 | **1,183** (same as above) | measured in this repository |

**No published 2D figures from competitors were found this time**: FluidX3D is 3D-only, M-Star is a 3D
product, and the representative published benchmarks of waLBerla/OpenLB/Palabos are also 3D. In other
words, the 2D table currently stands "alone," and is not the main battleground for external speed claims.
D2Q9 moves roughly half the data per update compared to D3Q19 (9 vs. 19 distribution functions), so
**placing D2Q9 MLUPS alongside D3Q19 in the same table inflates it by roughly 2×**. Always keep them
separate in external material.

## 7. Published multi-node / scaling track record

| Code | Machine / scale | Track record | Source |
|---|---|---|---|
| waLBerla | JUQUEEN (BG/Q) 458,752 cores / 1.8M threads | Over 1 trillion cells, up to 1.93 trillion cell updates/s (abstract-stated value; direct retrieval from the ACM page was not possible, so visual confirmation is needed). Strong scaling demonstrated on SuperMUC with 32,768 cores | [SC13 DOI:10.1145/2503210.2503273](https://dl.acm.org/doi/10.1145/2503210.2503273) |
| waLBerla | JUQUEEN | "Largest simulation exceeds 1 trillion cells," "good scalability to over 400,000 cores" (confirmed by direct quotation from the body text) | [arXiv:1511.07261 (ar5iv)](https://ar5iv.labs.arxiv.org/html/1511.07261) |
| waLBerla-wind | JUWELS Booster 30 nodes / 120 A100 | "Per-GPU performance is nearly constant" under weak scaling (17.5M cells/GPU, average 74.46 steps/s ≈ roughly 1,300 MLUPS per GPU converted) | [arXiv:2402.13171](https://arxiv.org/html/2402.13171v1) |
| OpenLB 1.5 | HoreKa 128 nodes / 512 A100 | Total ≈1.33 TLUPS (D3Q19 FP32 cavity). 64→128 GPU strong-scaling efficiency 0.64–0.81 (grid 575³–2300³). 92% of benchmark performance on a real turbulent-nozzle case with LES (224 GPU) | [openlb.net/performance](https://www.openlb.net/performance/) |
| OpenLB 1.5 | HoreKa 2 nodes / 8 A100 | 1000³ FP32 cavity 42.2 GLUPS (1-node 4×A100 is 24.8 GLUPS, 2-node CPU AVX-512 is 2.7 GLUPS → GPU speedup 15.6×) | [OpenLB 1.5 release](https://www.openlb.net/news/openlb-release-1-5-available-for-download/) |
| OpenLB 1.9 | Aurora 1,000 nodes (≈10% of the system) | Peak 21,120 GLUPS, 4 trillion cells (D3Q19 FP32) | [openlb.net/performance](https://www.openlb.net/performance/) |
| OpenLB (2026) | HoreKa heterogeneous (3 CPU/GPU partitions) | Up to 18G cells, strong-scaling efficiency 0.66–0.91 (by segment), single GPU node capable of ~1e9 cells | [arXiv:2506.21804](https://arxiv.org/html/2506.21804v1) |
| Palabos GPU | DGX 4× A100 40GB | 80–90% of weak-scaling ideal, 65–80% strong scaling | [arXiv:2506.09242](https://arxiv.org/html/2506.09242v1) |
| M-Star | 8× V100 / 8× A100 (single cloud node) | Publishes 1/2/4/8 GPU scaling charts (values are in figures, §3) | [Scaling Performance](https://docs.mstarcfd.com/19_Scaling_Performance/txt-files/Scaling-performance-index.html) |
| **LBMFlow (us)** | — | **no track record** (MPI planned for M-D. Acceptance criterion is single-node multi-rank weak scaling ≥85%) | COMPETITIVE_SPEC.md §3 R3 |

## 8. Notes on differences in conditions (required reading before the comparison tables)

1. **The definition of MLUPS is shared, but the measurement method differs.** The definition is "million
   lattice-point updates per second" (the waLBerla paper defines "MLUP/s per core = number of cells one
   core updates per second"
   — [ar5iv:1511.07261](https://ar5iv.labs.arxiv.org/html/1511.07261)). However, what is included in one
   "update" (boundary handling, output, communication) varies by code.
2. **Discount for stencil**: D3Q27 moves more data per update than D3Q19 (27 vs. 19 DDFs), so under
   bandwidth-bound conditions MLUPS is structurally lower. Reading waLBerla-wind's 1,677 MLUPS
   (D3Q27 cumulant) alongside FluidX3D's 8,526 (D3Q19 SRT) as "waLBerla is slow" would be a mistake.
   2D D2Q9 is in an entirely separate category (§6).
3. **Discount for precision**: in bandwidth-bound LBM, halving storage precision roughly doubles MLUPS
   (FluidX3D: on RTX 4090, FP32 5,624 → FP16S 11,091 — [README](https://github.com/ProjectPhysX/FluidX3D).
   Palabos: FP32 theoretical peak 9,481 vs. FP64 4,921 GLUPS
   — [arXiv:2506.09242](https://arxiv.org/html/2506.09242v1)). On the precision side, Lehmann et al.
   report that "the accuracy difference between FP64 and FP32 is negligible in almost all cases" and that
   "16-bit is sufficient in many cases"
   — [arXiv:2112.08926 / Phys. Rev. E 106, 015308](https://arxiv.org/abs/2112.08926). Comparison tables
   must always state "compute precision / storage precision" together.
4. **Grid-size dependence**: small grids underperform. FluidX3D measures with an empty box of "sufficient
   size (typically 256³)." The Palabos GPU paper explicitly states "performance is an increasing function
   of mesh resolution" ([arXiv:2506.09242](https://arxiv.org/html/2506.09242v1)). M-Star, too, uses "30–60M
   grid points or more per GPU" as its rule of thumb for scaling efficiency
   ([docs](https://docs.mstarcfd.com/19_Scaling_Performance/txt-files/Scaling-performance-index.html)).
   → Our formal measurements should also sweep grid size to report both the saturation value and the size
   at which it is reached.
5. **Kernel-alone vs. full solver**: even on the same code and same GPU, an empty-box kernel and the full
   physics can differ by 4–5× (waLBerla-wind: roofline 7,513 → measured 1,677 MLUPS with turbine = 22.3%
   — [arXiv:2402.13171](https://arxiv.org/html/2402.13171v1). OpenLB offers a good counterexample too: 92%
   of benchmark performance on a real turbulent-nozzle case
   — [openlb.net/performance](https://www.openlb.net/performance/)). FluidX3D's table is the former,
   M-Star's chart is the latter (stirred tank + particles). **Do not compare these two directly.**
6. **Averaging window / warm-up**: M-Star publishes the "average MLUPS over the whole simulation"
   ([docs](https://docs.mstarcfd.com/19_Scaling_Performance/txt-files/Scaling-performance-index.html)).
   Whether initialization, JIT, and cache effects are included moves the number. Our formal measurements
   should explicitly state warm-up exclusion.
7. **The range of implementation optimization itself is large**: even on identical hardware, performance
   varies greatly with data layout and streaming implementation (this is the very point of LBM
   benchmark-kernel research — [arXiv:1711.11468](https://arxiv.org/abs/1711.11468). Concrete example:
   naive 211 → optimized 550 MLUPS
   — [multiphase code-gen paper](https://journals.sagepub.com/doi/full/10.1177/10943420211016525)).
   "Code vs. code" differences should be argued separately from hardware differences and condition
   differences.
8. **Handling of vendor-published bandwidth and roofline**: FluidX3D publishes arithmetic intensity and
   explains it via roofline ([README](https://github.com/ProjectPhysX/FluidX3D)). It would be consistent
   with our verification culture to publish our own figures as a three-part set: "measured bandwidth →
   theoretical ceiling → measured MLUPS."

## 9. Where we are honestly losing (as of now)

1. **No 3D GPU exists.** We cannot yet produce a comparable figure on the competitors' main battleground
   (single-GPU D3Q19). We only get onto the table once we meet R2's acceptance criterion (single-GPU D3Q19
   f32 ≥1,500 MLUPS).
2. **FP16 storage mode is not implemented.** FluidX3D has demonstrated roughly 2× over FP32 with FP16S
   (e.g., 4090: 11,091 vs. 5,624) ([README](https://github.com/ProjectPhysX/FluidX3D)). Awaiting Pillar-4
   implementation.
3. **No way to reach high-end NVIDIA/AMD figures.** Numbers at the level of H100 NVL 32,922 / MI300X
   41,327 MLUPS (FP16S, [README](https://github.com/ProjectPhysX/FluidX3D)) would be physically
   unreachable on Apple Silicon bandwidth even with a flawless 3D GPU implementation. A CUDA/HIP backend
   and access to real hardware are required (SPEC §5).
4. **Zero multi-node track record.** waLBerla reaches trillions of cells across 400,000+ cores
   ([ar5iv:1511.07261](https://ar5iv.labs.arxiv.org/html/1511.07261)); OpenLB reaches 1.33 TLUPS on 512
   GPUs / 4 trillion cells on Aurora ([openlb.net](https://www.openlb.net/performance/)). Our R3 target
   (64-rank weak scaling ≥80%) is merely an initial goal many orders of magnitude below these.
5. **No full-physics benchmark.** M-Star's chart is measured on a "sellable workload" of stirred tank +
   particles. We only have figures for empty-box-style kernels, and do not yet know the performance drop
   when LES, moving boundaries, and scalar transport are layered on (in waLBerla-wind's case, -78%).
6. **Our own figures themselves are provisional.** They were measured under shared load, and publishing
   them would itself violate our verification culture (all claims tied to actual measurement). Do not
   release externally until §10 is complete.

## 10. Points that are ours alone (differentiation to argue in the body text, not as a table footnote)

- **Bundled verification suite, one-command re-run** (56+ adversarial verification tests, Ghia/Schäfer-
  Turek/RT/equivariance 4e-16). Competitors' verification is mostly in the form of "papers / published
  benchmark collections"; making re-runnability in the customer's own environment part of the product spec
  is unique to us (M-Star does provide a scaling benchmark package from v3.3.140+, so "performance
  re-running" is partially possible
  — [docs](https://docs.mstarcfd.com/19_Scaling_Performance/txt-files/Scaling-performance-index.html).
  Bundling **re-runnable physical-accuracy verification** is our differentiator).
- **Agent-native** (JSON Schema self-description + MCP. M-Star's Python API is human-oriented; FluidX3D
  requires editing a C++ setup).
- **Portability**: wgpu (Metal/Vulkan/DX12) + WASM. M-Star assumes NVIDIA
  ([Hardware Guide](https://docs.mstarcfd.com/2_Installation/txt-files/hardware.html)). FluidX3D supports
  all vendors via OpenCL but has a non-commercial license and no MPI. → The combination of
  "commercially usable, all-vendor, down to the browser" is ours alone.
- **Transparency of precision**: even when introducing FP16, we quantify "what is lost" via the
  verification suite (translating Lehmann et al.'s findings
  [arXiv:2112.08926](https://arxiv.org/abs/2112.08926) into the product spec).

## 11. Formal measurement checklist on an idle machine (required before publication)

Environment:
- [ ] Physical idle machine with no shared load (this M5 Max). Record power connection, thermal
      steady-state, and other apps closed
- [ ] Record OS / wgpu / driver / compiler versions and memory configuration / rated bandwidth
- [ ] Take measured memory bandwidth (STREAM-equivalent), compute the roofline ceiling first, and report
      it alongside

Measurement protocol:
- [ ] A measurement window excluding N warm-up steps, median of 5+ runs, with variance recorded
- [ ] Publish saturation curves via a grid-size sweep (2D: 512²–4096², 3D: 128³–memory limit)
- [ ] State the MLUPS definition explicitly (total cells × steps / wall time, noting exclusion of output
      I/O)
- [ ] Report "kernel alone (empty box, minimal boundary)" and "representative scenario (with obstacle +
      output)" separately

Aligning comparison conditions (when going head-to-head with FluidX3D):
- [ ] Build and run FluidX3D on the **same machine** (possible via OpenCL on Apple Silicon), and compare
      against a same-environment measurement rather than the README values (their measurement conditions:
      D3Q19 SRT, empty box, typical 256³, FP32 compute)
- [ ] Always label both stencil (D3Q19 vs. D2Q9) and storage precision (f32 vs. FP16S/C) together
- [ ] Footnote the difference between our deviation storage f32 and their FP32 (noting that it has been
      precision-verified)

Publication:
- [ ] Replace our own provisional figures in this draft with formal values, and remove the "provisional"
      annotation
- [ ] If finalizing the visually-read values from the M-Star charts, state the reading method and error
      margin
- [ ] Primary confirmation of the waLBerla SC13 abstract figures (1.93 trillion cell updates/s, 1.8M
      threads) (visual check of the ACM page or the paper PDF)
- [ ] After 3D GPU implementation: compare the M5 Max single-GPU D3Q19 f32 measurement against the Apple
      rows in Table 2 (M2 Max/M2 Ultra)
- [ ] After FP16 storage implementation: re-compare against Table 2 in FP16S-equivalent mode + publish the
      verification suite's degradation figures at the same time

## 12. Source list

| # | Source | Use |
|---|---|---|
| 1 | https://github.com/ProjectPhysX/FluidX3D | FluidX3D benchmark table, measurement conditions, arithmetic intensity, memory/cell |
| 2 | https://raw.githubusercontent.com/ProjectPhysX/FluidX3D/master/README.md | same (raw data, retrieved 2026-07-05) |
| 3 | https://docs.mstarcfd.com/19_Scaling_Performance/txt-files/Scaling-performance-index.html | M-Star scaling benchmark conditions, platform, rules of thumb |
| 4 | https://docs.mstarcfd.com/_images/gce-a2-highgpu-8g_agitated.svg / …_pipe.svg | M-Star chart originals (axis range confirmation) |
| 5 | https://docs.mstarcfd.com/2_Installation/txt-files/hardware.html | M-Star NVIDIA requirements, VRAM rule of thumb |
| 6 | https://mstarcfd.com/software/ | M-Star marketing claims (qualitative) |
| 7 | https://arxiv.org/html/2402.13171v1 | waLBerla-wind: A100 D3Q27 measured 1,677 / roofline 7,513, EPYC 7763 measured, weak scaling |
| 8 | https://ar5iv.labs.arxiv.org/html/1511.07261 | waLBerla direct quotation on trillion cells / 400,000+ cores, MLUP/s definition |
| 9 | https://dl.acm.org/doi/10.1145/2503210.2503273 | waLBerla SC13 (1.93 trillion updates/s etc. from the abstract, requires visual confirmation) |
| 10 | https://www.openlb.net/performance/ | OpenLB 512×A100 1.33 TLUPS, Aurora 21,120 GLUPS, Magnus 142,479 MLUPS, strong-scaling efficiency |
| 11 | https://www.openlb.net/news/openlb-release-1-5-available-for-download/ | OpenLB 1000³ cavity 42.2/24.8/2.7 GLUPS, GPU 15.6× |
| 12 | https://arxiv.org/html/2506.21804v1 | OpenLB heterogeneous HPC, 18G cells, efficiency 0.66–0.91 |
| 13 | https://arxiv.org/html/2506.09242v1 | Palabos GPU: A100 theoretical peak 9.481/4.921 GLUPS, measured 75–85%, scaling, mentions waLBerla 85% |
| 14 | https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0250306 | STLBM: GPU/CPU D3Q19 FP64 MLUPS |
| 15 | https://arxiv.org/abs/2112.08926 | FP64/FP32/16-bit precision impact (Phys. Rev. E 106, 015308) |
| 16 | https://arxiv.org/abs/1711.11468 | LBM benchmark kernels (impact of implementation differences, methodology) |
| 17 | https://journals.sagepub.com/doi/full/10.1177/10943420211016525 | Concrete example of optimization range (211→550 MLUPS) |
| 18 | measured in this repository (PERFORMANCE.md lineage) | our own figures (under shared load, provisional) |

**Figures withheld from this listing as unconfirmed**: OpenLB single-A100's "8.3 GLUPS (Periodic Shift,
D3Q19 FP32)" — this appears in a search snippet in a
[third-party paper abstract on ScienceDirect](https://www.sciencedirect.com/science/article/abs/pii/S0010465522003228),
but neither that page nor the primary paper
([Wiley cpe.7509](https://onlinelibrary.wiley.com/doi/full/10.1002/cpe.7509)) could be directly retrieved
(403/402), so it was not adopted. waLBerla's "889,602 MLUPS" was also not adopted because its source
machine attribution could not be confirmed (it is inconsistent with the 1.93 trillion updates/s-series
figure).
