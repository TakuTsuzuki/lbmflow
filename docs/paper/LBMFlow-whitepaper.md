# LBMFlow: A CFD Engine You Can Trust and Automate

**Technical paper — living draft.** Describes the current, measured product.
`docs/paper/claims-ledger.md` tracks measurement status per claim; the paper is
updated when measurements change. Figures: Apple M5 Max, 18 cores, 128 GB, macOS.
Reproduction commands in Appendix A. MLUPS = million lattice-cell updates per second.

---

## 1. Executive summary

LBMFlow is a commercial-grade Lattice Boltzmann (LBM) fluid simulator that:

- **Reproduces its physics on your hardware.** T1–T17 (200+ adversarial tests,
  written against the spec by an independent agent) reproduces Ghia,
  Schäfer–Turek, Taylor–Green, Shan–Chen, and Rayleigh–Taylor references; domain
  decomposition is **bit-for-bit identical** to a monolithic run; rotational
  equivariance holds to 4×10⁻¹⁶. Re-runnable with one command.
- **Is agent-native.** JSON scenario schema, seven MCP tools including an
  async job lifecycle, machine-readable divergence/stability diagnostics.
- **Runs anywhere without CUDA lock-in.** Same scenario JSON runs in a browser
  (WebAssembly), on the CPU (SIMD), on the GPU via wgpu (Metal / Vulkan / DX12),
  and across an MPI cluster. On the M5 Max: **2D GPU 7,073 MLUPS** at 1024²,
  **3D GPU 2,791–2,813 MLUPS** at 192³ (D3Q19, quiet-window A/B/A), CPU
  **1,480 MLUPS** (2D) and **302 MLUPS** (3D). FP16 storage adds ~2× MLUPS at
  2048² and doubles usable grid capacity (D3Q19 f16 >5 GLUPS).

## 2. Positioning

Among commercially licensed engines, the fastest GPU LBM codes are CUDA-only.
The fastest LBM code overall — FluidX3D — is cross-vendor via OpenCL but
non-commercial and single-node. LBMFlow closes three gaps at once: a re-runnable
physics validation surface (not a slide deck), an agent-driven control plane
(not human-in-the-loop C++/Python editing), and cross-vendor GPU via wgpu
(Metal today; Vulkan / DX12 share the code path).

## 3. Method and architecture

LBMFlow implements the Lattice Boltzmann method on the D2Q9 (2D) and D3Q19 (3D)
lattices, with BGK and TRT collision operators (TRT with the magic parameter
Λ=3/16 is the default: it makes the wall position exact and is as fast as BGK).
Body forces use second-order Guo forcing; solid walls are half-way bounce-back;
open boundaries use a normal-parameterized Zou–He implementation covering velocity
inlets and pressure outlets on any face. Force on immersed bodies is measured by the
momentum-exchange method. The relaxation time is τ = 3ν + ½ (cs² = ⅓).

The engine is one orthogonal core parameterized over **dimension × lattice ×
precision × backend × partition**. That orthogonality is *why* one scenario is
portable: the same physics definition is dispatched to a scalar CPU backend, a
SIMD-vectorized CPU backend, a GPU compute backend (wgpu → Metal / Vulkan / DX12), or
an MPI-partitioned run, without changing the scenario. A deviation-storage scheme
(distributions are stored as f−wᵢ) keeps single-precision (f32) at validation grade —
the quiescent background is exactly zero, so f32 rounding acts only on the
fluctuation scale.

## 4. Verification: evidence you can re-run

Split a domain across four sub-blocks, or across MPI ranks, gather the result, and
subtract it from the single-block run: the difference is **max|Δ| = 0.0** — not
"within tolerance," bit-for-bit identical. Rotate an entire problem 90° and the
solution reproduces to **4×10⁻¹⁶**, machine precision. These are the kinds of guarantee
LBMFlow is built to make, and they are why the verification story is the one no
competitor ships in re-runnable form.

**Adversarial validation as a shipped product.** LBMFlow's test suite is written
against the published specification by an author who is not the engine's implementer.
The separation is deliberate: it catches specification bugs, not just coding bugs.
The suite reproduces, with frozen numeric tolerances (all CPU, f64 where noted):

- **Taylor–Green vortex**: 2nd-order convergence (measured order 1.91), effective
  viscosity within ±2% of nominal.
- **Poiseuille flow (TRT)**: exact to ≤1×10⁻¹⁰ (half-way bounce-back is analytically
  exact).
- **Lid-driven cavity vs Ghia et al. (1982)**: centerline RMS within tolerance at
  Re = 100/400/1000; 3D cavity vs Albensoeder–Kuhlmann (2005) at Re = 1000.
- **Cylinder vs Schäfer–Turek**: drag/lift coefficients and Strouhal number in the
  benchmark bands (2D-1 steady, 2D-2 vortex shedding).
- **3D duct**: exact Fourier-series solution to L∞rel 2.3×10⁻⁴; sphere drag vs
  Schiller–Naumann within ±10%.
- **Multiphase (Shan–Chen)**: Laplace law R² = 0.9999, contact-angle full range,
  Rayleigh–Taylor growth rate γ within 12% of the tension-and-viscosity-corrected
  reference.

**Bit-exact partition invariance and backend equivalence.** The partition invariance
above holds against adversarial attacks (obstacles and probes straddling the seam,
lids crossing splits, multiphase droplets in a shared corner). The GPU backend
reproduces the CPU trajectory to ≤1×10⁻⁵ relative across six 2D scenario classes; a
one-ulp control test pins the residual to rounding order, not physics.

**One command.** The CPU and physics validation runs with
`cargo test --workspace --release` (heavy cases with `--include-ignored`); the
GPU-equivalence checks add `--features gpu` on a GPU host. Your evaluation team
re-runs the evidence base in your environment — the reproduction *is* the product.

## 5. Performance: scoped, honest, measured

Protocol: memory-bandwidth roofline first, warm-up excluded, best-of-N,
grid-size sweep, I/O excluded. Reproduction commands in Appendix A; raw data
in `docs/paper/benchmark-results.md`.

**Roofline.** CPU memory bandwidth 344 GB/s (STREAM) / 459 GB/s
(write-allocate) at 18 threads, against ~546 GB/s nominal. LBM is
bandwidth-bound.

**GPU, f32 (Metal via wgpu; Vulkan/DX12 share the wgpu path):**

| case | MLUPS | note |
|---|---|---|
| 2D D2Q9 1024² | **7,073** | sustained |
| 2D D2Q9 2048² | 6,720 | sustained |
| 2D D2Q9 512² | 12,205 | cache-resident (fits SLC), not sustained |
| 3D D3Q19 192³ | **2,791–2,813** | quiet-window A/B/A, target ≥1,500 |
| 3D D3Q19 128³ | 2,778–2,880 | quiet-window A/B/A |

3D GPU (D3Q19) at 192³ exceeds the ≥1,500 MLUPS acceptance line by ~1.85×.
The GPU path does not exist for OpenLB or Palabos on Apple Silicon (CUDA-only).

**FP16 storage (compute stays f32), same GPU:**
~2.0× MLUPS at 2048², D3Q19 f16 >5 GLUPS, ×2 grid capacity inherent.
Accuracy bands frozen and measured: TGV transient 1.401×10⁻¹ (band 2×10⁻¹),
cavity steady 2.579×10⁻³ (band 5×10⁻³). Steady-vs-transient dichotomy in PHYSICS.md.

**CPU, fused SIMD backend, f32:**

| case | LBMFlow MLUPS |
|---|---|
| 2D D2Q9 2048² 18T | **1,480** |
| 2D D2Q9 1024² 18T | 1,208 |
| 3D D3Q19 192³ 18T | 302 |
| 3D D3Q19 128³ 18T | 267 |

**Head-to-head vs OpenLB 1.9**, same machine, 3D D3Q19 128³, f32. **Fairness
caveat:** OpenLB was built native arm64 but runs its scalar CPU platform
(CPU_SISD); its vectorized path targets x86 AVX and has no ARM NEON backend
out of the box. LBMFlow runs NEON-vectorized. Not SIMD-vs-SIMD.

| configuration | LBMFlow | OpenLB | |
|---|---|---|---|
| single thread | **52.0** | 44.6 | LBMFlow +17% |
| 18-way | 266.6 | **298.8** | OpenLB +12% |

At 18 cores OpenLB's MPI domain decomposition scales better than LBMFlow's
current thread parallelism on this heterogeneous-core laptop.

**Precision transparency.** Deviation-storage f32 is validation-grade: uniform-force
momentum error 1.3×10⁻³ → 2.8×10⁻⁷; Taylor–Green f32 error 7.1×10⁻⁴ vs f64 7.0×10⁻⁴.

## 6. Agent-native operation

LBMFlow is built to be driven by software, not just people. A scenario is a
self-describing JSON contract; `lbm schema` emits its full specification; and an MCP
server exposes seven tools including an asynchronous job lifecycle
(`start_run` → `run_status` → `list_runs`) so an agent can launch long or parallel
simulations and poll them to completion. The diagnostic surface is machine-readable:
a divergence returns a structured reason and a stability hint, not a stack trace.
Runs are deterministic run-to-run at a fixed configuration, so an agent's sweep is
reproducible.

This makes LBMFlow the substrate for automated design exploration: an agent discovers
the schema, generates a parameter sweep, runs it unattended across the portability
ladder, reads structured results, and iterates — no human in the inner loop.

## 7. Roadmap — remaining acceptance items

Landed (measured, in §5 above): 3D GPU D3Q19 ≥1,500 MLUPS; FP16 storage
≥1.5× MLUPS at 2048² with frozen accuracy bands.

Remaining (targets, measurement pending — status in claims-ledger):

- **Multi-node scaling.** ≥80% weak scaling at 64 ranks on an MPI cluster;
  single-node weak scaling already at 97–99% for ≤4 cores.
- **Full-physics stirred workload.** Two-phase + particles + scalar transport
  + LES; publish the performance-degradation ratio vs single-phase kernel.

## 8. Reproduce

Every physics claim and performance number reproduces with the commands in
Appendix A. CPU/physics: `cargo test --workspace --release`. GPU: add
`--features gpu` on a GPU host.

---

## Appendix A — reproduction & fairness caveats

| Claim | Number | Command |
|---|---|---|
| 2D GPU MLUPS (Metal) | 7,073 (1024²); 6,720 (2048²) | `cargo run --release --features gpu -p lbm-core --example bench_gpu` |
| 3D GPU MLUPS (Metal) | 2,791–2,813 (192³); 2,778–2,880 (128³) | `cargo run --release --features gpu -p lbm-core --example bench_gpu3d` |
| FP16 GPU MLUPS | ~2× f32 @ 2048²; D3Q19 f16 >5 GLUPS | `... bench_gpu --precision f16` |
| 2D CPU MLUPS | 1,480 (2048²/18T/f32) | `cargo run --release -p lbm-core --example bench_backends -- simd f32 2048 18 400` |
| 3D CPU MLUPS | 302 (192³/18T/f32); 267 (128³) | `... bench_backends -- simd f32 192 18 100 192` |
| Roofline BW | 344 / 459 GB/s | `~/projects/cfd-bench/bw_triad 18 1024` |
| OpenLB head-to-head | §5 (128³) | `~/projects/cfd-bench/run_openlb_sweep.sh` |
| Validation suite | all green | `cargo test --workspace --release` |
| Partition invariance | max\|Δ\|=0.0 | `cargo test --release t13` |
| Backend equivalence (2D+3D) | ≤1e-5 | `cargo test --release --features gpu t14` |

**Fairness caveats (read before comparing).** (1) *Stencil*: 2D D2Q9 moves ~half the
bytes per cell of 3D D3Q19 — never compare 2D and 3D MLUPS directly. (2) *Precision*:
all headline figures are f32; f64 roughly halves MLUPS (bandwidth-bound). (3) *Cache*:
512² fits the on-chip SLC and is not a sustained-bandwidth number. (4) *SIMD*: the
OpenLB comparison is LBMFlow-NEON vs OpenLB-scalar (no ARM NEON path in OpenLB);
single-thread is not a SIMD-vs-SIMD result. (5) *Kernel vs full solver*: these are
lid-driven-cavity / periodic kernels; a full scenario with obstacles and output is
slower. (6) *GPU vendor*: GPU numbers are Metal; Vulkan/DX12 share the wgpu path but
are unmeasured here.

## Appendix B — validation inventory
See `docs/VALIDATION.md` (T1…T17) for the full specification and acceptance criteria.

## Appendix C — provenance
Apple M5 Max, 18 cores, 128 GB, macOS; rustc 1.93 `--release` (thin LTO,
codegen-units=1); GPU via wgpu/Metal; OpenLB 1.9 native arm64 (Apple clang 21,
arm64 Open MPI, CPU_SISD). 2D CPU/GPU baseline: LBMFlow HEAD `b262447`
(2026-07-05 window). 3D GPU + FP16: measured post-ME-1/ME-2 landings 2026-07-06
(quiet-window A/B/A). Before external release: confirm the referenced commits
are pushed to a public ref.

## Appendix D — claim → evidence traceability
See `docs/paper/claims-ledger.md` — status snapshot mapping each claim to
its implementing item and measurement.

## Appendix E — published competitive landscape (context, not same-machine)
Sources in `docs/BENCH_COMPARISON_DRAFT.md`. Single-GPU 3D D3Q19 is led by
FluidX3D (A100 8,526 MLUPS FP32 / 16,035 FP16S; non-commercial, single-node);
cluster scale by waLBerla (trillion-cell, >400k cores) and OpenLB (1.33 TLUPS on
512× A100). LBMFlow's present differentiation: commercial license, cross-vendor
GPU (Metal/Vulkan/DX12) down to the browser, agent-native control, and a
customer-re-runnable physics validation suite.
