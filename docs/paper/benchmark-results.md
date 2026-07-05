# Benchmark Results — public-grade measurement window #1 (2026-07-05 23:11–23:35 JST)

Machine: Apple M5 Max, 18 cores, 128 GB, macOS. LBMFlow @ `feat/body-force-field-api`
HEAD `b262447` (CpuSimd fused backend; GPU = `lbm-gpu-proto`, Metal/wgpu).
Idle window (user stopped non-LBM jobs; load1≈3.4). Warmup excluded via best-of-N.
All numbers are **measured Tier-1** (claims-ledger GREEN). Raw CSVs + scripts under
`~/projects/cfd-bench/` (outside repo): `run_lbmflow_cpu_sweep.sh`,
`run_openlb_sweep.sh`, `bw_triad.c`.

## 0. Roofline — memory bandwidth (native arm64 triad, 18 threads)
| metric | GB/s |
|---|---|
| STREAM convention (24 B/elem) | **344** |
| write-allocate (32 B/elem) | **459–464** |
| nominal (spec) | ~546 |
| single-thread reference | 121 / 161 |

LBM per-cell traffic (f32 pull): D2Q9 = 72 B/cell·step, D3Q19 = 152 B/cell·step.
Roofline MLUPS ceiling (write-alloc 459 GB/s): **2D ≈ 6,375**, **3D ≈ 3,020**.

## 1. LBMFlow CPU — CpuSimd, best-of-5, MLUPS
### 2D D2Q9
| grid | 1T f32 | 1T f64 | 18T f32 | 18T f64 |
|---|---|---|---|---|
| 512² | 280 | 144 | 809 | 577 |
| 1024² | 258 | 149 | **1208** | 711 |
| 2048² | 279 | 126 | **1480** | 734 |
| 4096² | 198 | 130 | 1324 | 687 |

Peak 2D = **1,480 MLUPS** (2048²/18T/f32) = 23% of the write-alloc roofline
(31% of STREAM). Beats the prior documented peak (1,183).

### 3D D3Q19
| grid | 1T f32 | 1T f64 | 18T f32 | 18T f64 |
|---|---|---|---|---|
| 128³ | 52 | 33 | 267 | 167 |
| 192³ | 54 | 35 | **302** | 189 |

## 2. LBMFlow GPU — Metal (wgpu), TRT f32, MLUPS (submit→wait inclusive)
| grid | GPU MLUPS | effective BW |
|---|---|---|
| 512² | **12,205** | (fits SLC cache) |
| 1024² | **7,073** | ~509 GB/s |
| 2048² | **6,720** | ~484 GB/s |

Bandwidth-bound and near the GPU memory ceiling. All CPU↔GPU verification checks
passed (physics identical to CPU). **On Apple Silicon this is a laptop-GPU path that
OpenLB/Palabos cannot take (their GPU backends are CUDA-only).**

## 3. Head-to-head — LBMFlow vs OpenLB 1.9, 3D D3Q19, f32, same machine & window
OpenLB built native arm64 (arm64 Open MPI, Apple clang 21), `cavity3dBenchmark`
(float D3Q19), **CPU_SISD** platform (OpenLB's CPU SIMD targets x86 AVX; no ARM NEON
path out of the box — a disclosed, fair "best-available on this hardware" comparison).
LBMFlow = CpuSimd (NEON autovectorized). Best-of runs.

| grid | LBMFlow | OpenLB | winner |
|---|---|---|---|
| 128³, 1 thread / 1 rank | **52.0** | 44.6 | LBMFlow **+17%** |
| 128³, 18 thread / 18 rank | 266.6 | **298.8** | OpenLB **+12%** |

Reading (honest): LBMFlow's NEON kernels win single-thread; OpenLB's MPI domain
decomposition scales better to all 18 cores (LBMFlow's rayon 3D scaling is a known
limit — band-edge double-collision + P/E heterogeneous cores, PERFORMANCE.md).
**Same order of magnitude = competitive on CPU 3D**, with LBMFlow additionally
offering the GPU path OpenLB lacks here.

OpenLB per-config best (MLUPS): 64³ 1r 43.4 / 18r 208.7 · 100³ 1r 43.9 / 18r 224.3 ·
128³ 1r 44.6 / 18r 298.8.

## 4. Not yet measured (follow-up window)
- OpenLB **2D** (cavity2d) for the D2Q9 head-to-head vs LBMFlow 1,480 MLUPS.
- OpenLB **CPU_SIMD** attempt on ARM (if supported) — for a SIMD-vs-SIMD line.
- **Palabos** (Family A #2) — needs cmake install then build.
- **OpenFOAM via Colima** (Family B) — cylinder Re=20 Cd time-to-solution.
- LBMFlow **scalar** baseline rows (sweep-script bug produced 0; re-run cleanly).
