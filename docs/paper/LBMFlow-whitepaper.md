# LBMFlow: A CFD Engine You Can Trust and Automate

**Technical paper — living draft.** Update the paper text as the implementation
lands or changes; describe the current, measured product, not an aspirational target.
`docs/paper/claims-ledger.md` is a status snapshot of what has been measured and
what has not — a working reference, not a release gate. Measured figures: Apple M5
Max, 18 cores, 128 GB, macOS; LBMFlow HEAD `b262447`. Every measured number
reproduces from that commit
(Appendix A). MLUPS = million lattice-cell updates per second, the standard LBM
throughput unit.

---

## 1. Executive summary

On an Apple M5 Max **laptop**, LBMFlow's GPU solver sustains **7,073 MLUPS** in 2D —
on a GPU its open-source and commercial competitors cannot use at all, because their
GPU backends are CUDA-only. That one fact is the shape of this product: fast on the
hardware you actually have, and honest about every number behind the claim.

Computational fluid dynamics is powerful and, today, hard to trust and hard to
automate. Vendor benchmarks are not reproducible on your hardware; validation is a
slide deck, not a command you can run; and the fastest *commercially licensed*
engines lock you to one GPU vendor and one operating model — a human editing C++ or
Python by hand.

LBMFlow is a commercial-grade Lattice Boltzmann (LBM) fluid simulator built to remove
all three frictions:

- **Trust by reproduction.** Physical correctness is not asserted, it is
  *demonstrated and re-runnable*. The T1–T17 validation suite (200+ adversarial
  tests) — authored against the specification by an independent agent, not the
  engine's author — reproduces the classical references (Ghia lid-driven cavity,
  Schäfer–Turek cylinder drag, Taylor–Green decay, Shan–Chen Laplace law,
  Rayleigh–Taylor growth) and proves exact rotational equivariance to 4×10⁻¹⁶.
  Domain decomposition is **bit-for-bit identical** to a monolithic run. You re-run
  it with one command.

- **Automate by design.** LBMFlow is operable end-to-end by an AI agent: a
  self-describing JSON scenario schema, seven MCP tools including an asynchronous job
  lifecycle, and a machine-readable diagnostic surface (divergence reasons, stability
  hints). Parameter sweeps and optimization loops run unattended.

- **Run anywhere, honestly fast.** The *same scenario JSON* runs in a browser
  (WebAssembly), on the CPU (SIMD), on the GPU via wgpu (Metal today; Vulkan / DX12
  by the same code path), and across an MPI cluster — **no CUDA, no vendor lock-in.**
  On the M5 Max laptop the GPU path sustains **7,073 MLUPS** at 1024², and the CPU
  path reaches **1,480 MLUPS** in 2D and **302 MLUPS** in 3D — the 3D figure within
  the same order as OpenLB on the same machine.

Where LBMFlow claims speed, it ships the number and the command to reproduce it —
and where a mode is faster, it publishes what that mode costs, in numbers.

## 2. The problem with the CFD you are buying today

Three structural gaps define the incumbent tools.

**The trust gap.** Commercial LBM vendors publish performance charts and a list of
validation cases, but you cannot re-run their *physics* validation in your own
environment — you take the marketing on faith. Open-source research codes publish
papers, but reproducing them is a research project in itself. When the answer feeds a
real engineering decision, "trust us" is not good enough.

**The automation gap.** Today's engines are built for a human at a keyboard —
editing C++ source, hand-writing Python scripts. None expose a machine-discoverable
contract that an AI agent can read, drive, and diagnose on its own. As simulation
moves into automated design-exploration and optimization loops, a human-in-the-loop
API is a bottleneck.

**The lock-in gap.** Among *commercially licensed* engines, the fastest GPU LBM
codes are NVIDIA/CUDA-only. (The fastest LBM code overall, the research tool
FluidX3D, is cross-vendor via OpenCL but non-commercial and single-node — not an
option for a commercial, cluster-bound workflow.) If your hardware is Apple Silicon,
AMD, an Intel GPU, or a web browser, the commercial options lock you out of GPU
acceleration entirely.

LBMFlow is designed from the core outward to close all three.

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

All figures are measured on one Apple M5 Max in a single idle window with a stated
protocol: memory-bandwidth roofline first, warm-up excluded, best-of-N, grid-size
sweep, I/O excluded. Full tables, reproduction commands, and the fairness caveats
table are in Appendix A; raw data in `docs/paper/benchmark-results.md`.

**Roofline.** Measured CPU memory bandwidth is 344 GB/s (STREAM convention) to
459 GB/s (write-allocate) at 18 threads, against ~546 GB/s nominal. LBM is
bandwidth-bound; these set the CPU ceiling every CPU number below is measured against.

**GPU (2D, Metal via wgpu), f32:**

| grid | MLUPS | note |
|---|---|---|
| 1024² | **7,073** | sustained, bandwidth-bound — the headline number |
| 2048² | 6,720 | sustained |
| 512² | 12,205 | *cache-resident* (working set fits the SLC); not a sustained figure |

The sustained 7,073 MLUPS at 1024² is ~6× LBMFlow's own all-core 2D CPU
(1,208 MLUPS at the same grid). Its significance is not the raw ratio — it is that
**this GPU path does not exist for OpenLB or Palabos on Apple Silicon (their GPU
backends are CUDA-only).** On the laptop and workstation hardware most engineers
have, LBMFlow's portability turns GPU acceleration into something the CUDA-locked
competition cannot run at all. Numbers are Metal-measured; Vulkan and DX12 use the
same wgpu code path (not separately benchmarked here).

**CPU (single node), fused SIMD backend, f32:**

| case | LBMFlow MLUPS |
|---|---|
| 2D D2Q9, 2048², 18 threads | **1,480** |
| 2D D2Q9, 1024², 18 threads | 1,208 |
| 3D D3Q19, 192³, 18 threads | 302 |
| 3D D3Q19, 128³, 18 threads | 267 |

**Head-to-head vs OpenLB 1.9**, same machine, same window, 3D D3Q19, 128³, f32.
**Fairness caveat (stated plainly):** OpenLB was built native arm64 but runs its
scalar CPU platform (CPU_SISD) — its vectorized path targets x86 AVX and has no ARM
NEON backend out of the box — while LBMFlow runs NEON-vectorized. This is a
"best-available on this hardware" comparison, not SIMD-vs-SIMD.

| configuration | LBMFlow | OpenLB | |
|---|---|---|---|
| single thread | **52.0** | 44.6 | LBMFlow +17% (NEON vs scalar) |
| 18-way | 266.6 | **298.8** | **OpenLB leads by 12%** |

We state it plainly: at all 18 cores OpenLB is 12% faster here — its MPI domain
decomposition scales better than LBMFlow's current thread parallelism on this
heterogeneous-core laptop. LBMFlow leads single-thread and is the same order of
magnitude all-core, and it carries the GPU, browser, and single-scenario portability
OpenLB does not. (OpenLB's own performance benchmark is 3D-only, so the LBMFlow 2D
CPU figures above have no same-machine OpenLB baseline yet; that measurement is a
scoped follow-up. Published cross-code and cluster figures — FluidX3D single-GPU,
waLBerla/OpenLB at tera-scale — are compiled with sources in the competitive
appendix.)

**Precision transparency.** Deviation-storage f32 is validation-grade, measured: the
uniform-force momentum error improves from 1.3×10⁻³ to 2.8×10⁻⁷, and Taylor–Green f32
error (7.1×10⁻⁴) is indistinguishable from f64 (7.0×10⁻⁴). Every precision mode is a
validation target — the paper tells you what the fast mode costs, in numbers.

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

## 7. Committed roadmap

Each item below is a dated commitment with a public acceptance criterion recorded in
`docs/PLAN.md` — the number we hold ourselves to, in order.

- **3D GPU acceleration (D3Q19).** The GPU backend is dimension-agnostic by
  construction; 3D runs on the GPU under the same CPU↔GPU bit-equivalence gate that
  governs the 2D path, against an acceptance line of **≥1,500 MLUPS** on a single GPU.
- **FP16 storage.** Memory-halving f16 storage (compute stays f32) doubles the grid
  that fits in GPU memory, with its accuracy cost *measured and frozen* — a ≥1.5×
  throughput gain at 2048² and a published degradation band.
- **Multi-node scaling.** The same scenario JSON runs across an MPI cluster;
  full-cluster weak scaling is certified at **≥80% at 64 ranks** (single-node weak
  scaling is already measured at 97–99% for ≤4 cores).
- **Full-physics workload.** The stirred-tank benchmark (two-phase + particles +
  scalar transport + LES) publishes its performance-degradation ratio against the
  single-phase kernel — the number that makes a full-physics comparison honest.

A vendor that states the number it will hit, in order, with the architecture already
in place, is the strongest possible expression of the trust this product is built on.

## 8. Conclusion — evaluate it yourself

Every physics claim and every performance number in this paper reproduces in your
environment (CPU/physics in one command; GPU numbers add `--features gpu` on a GPU
host). That is the whole thesis: a CFD engine whose correctness you can verify, whose
operation you can automate, and whose speed you can confirm on the hardware you
already own — from a browser tab to an MPI cluster, without changing a line of your
scenario. The evidence is the product.

---

## Appendix A — reproduction & fairness caveats

| Claim | Number | Command |
|---|---|---|
| 2D GPU MLUPS (Metal) | 7,073 (1024²); 6,720 (2048²) | `cargo run --release --features gpu -p lbm-core --example bench_gpu` |
| 2D CPU MLUPS | 1,480 (2048²/18T/f32) | `cargo run --release -p lbm-core --example bench_backends -- simd f32 2048 18 400` |
| 3D CPU MLUPS | 302 (192³/18T/f32); 267 (128³) | `... bench_backends -- simd f32 192 18 100 192` |
| Roofline BW | 344 / 459 GB/s | `~/projects/cfd-bench/bw_triad 18 1024` |
| OpenLB head-to-head | §5 (128³) | `~/projects/cfd-bench/run_openlb_sweep.sh` |
| Validation suite | all green | `cargo test --workspace --release` |
| Partition invariance | max\|Δ\|=0.0 | `cargo test --release t13` |
| Backend equivalence (2D) | ≤1e-5 | `cargo test --release --features gpu t14` |

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
codegen-units=1); GPU via wgpu/Metal; OpenLB 1.9 native arm64 (Apple clang 21, arm64
Open MPI, CPU_SISD). LBMFlow HEAD `b262447`. Measurement window 2026-07-05 (idle
machine). **Before external release**: confirm `b262447` is pushed to a public ref so
the reproduction commands resolve for a customer clone.

## Appendix D — claim → evidence traceability
See `docs/paper/claims-ledger.md` — a working measurement-status snapshot
mapping each paper claim to its implementing item and the measurement that
verifies it. It is a status ledger, not a release gate: the paper describes
what is measured today, and is updated when the measurements change.

## Appendix E — published competitive landscape (context, not same-machine)
Full sourced figures in `docs/BENCH_COMPARISON_DRAFT.md`. Summary of the landscape
these numbers sit in (all from cited public sources; different hardware, not
same-machine): single-GPU 3D D3Q19 is led by FluidX3D (A100 8,526 MLUPS FP32 /
16,035 FP16S; non-commercial license, single-node); cluster scale is led by
waLBerla (trillion-cell, >400k cores) and OpenLB (1.33 TLUPS on 512× A100). LBMFlow
does not claim to beat these on their hardware today — §7 states the committed path
(3D GPU ≥1,500 MLUPS single-GPU; 64-rank weak scaling ≥80%). LBMFlow's present,
measured differentiation is the axis none of them occupy at once: commercial license,
all-vendor GPU (Metal/Vulkan/DX12) down to the browser, agent-native control, and a
customer-re-runnable physics validation suite.
