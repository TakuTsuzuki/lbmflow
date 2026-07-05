# PERFORMANCE.md — Measured Performance and Precision/Speed Trade-off Guide

Measurement environment: Apple Silicon 18-core (P/E hybrid), macOS, rustc 1.93,
`--release` (thin LTO, codegen-units=1). Measured with: `examples/bench_backends.rs`
(the V1-era `bench_mlups.rs` was retired along with V1 — its measurement records
from that time are frozen in the Phase 9 section of this document). MLUPS = million
lattice-point updates per second. As of 2026-07-05.

## V2 CpuSimd backend (M-E, 2026-07-05)

The V1 fused kernel (`step_band`) ported into V2 as a `Backend`
(`crates/lbm-core/src/backend_simd.rs`). D2Q9 uses the same fully-fused form as V1
(collide+stream+moments in 1 pass, row-band parallel, 3-slab ring); D3Q19 uses the
same structure (z-plane slabs, blocked span kernel). Halo exchange remains at the
phase boundary, so it composes as-is with Subdomain / InProcess / MPI. Equivalence is
verified by `tests/backend_simd_equiv.rs` (CpuScalar vs. 8 scenarios × f64/f32; fields
show a measured difference of ~1e-13/f64, all gates green).

Measurement: `cargo run --release -p lbm-core --example bench_backends`
(runs V1 and V2 alternately from the same binary at the same time, best-of-3. Since
this is a shared machine, absolute rates vary by time of day. Measured in a window at
load≈7).

### 2D D2Q9 (512²/1024², TGV-type periodic, TRT)

| grid | threads | prec | V1 fused | CpuScalar | CpuSimd | Simd/V1 | Simd/Scalar |
|---|---|---|---|---|---|---|---|
| 512² | 1 | f32 | 232 | 40 | 273 | **1.18** | 6.9 |
| 512² | 1 | f64 | 105 | 40 | 139 | **1.33** | 3.5 |
| 512² | 12 | f32 | 828 | 288 | 798 | 0.96 | 2.8 |
| 512² | 12 | f64 | 557 | 272 | 542 | 0.97 | 2.0 |
| 1024² | 1 | f32 | 216 | 39 | 252 | **1.17** | 6.4 |
| 1024² | 1 | f64 | 119 | 40 | 148 | **1.24** | 3.7 |
| 1024² | 12 | f32 | 1084 | 314 | **1183** | 1.09 | 3.8 |
| 1024² | 12 | f64 | 653 | 296 | 716 | 1.10 | 2.4 |

The target of "within −10% of V1" was achieved on all configurations. Single-thread
exceeds V1 (+17–33%): ensure_slab fuses the copy into the ring into the collision
(reads f, writes directly into the ring), so it has one pass less traffic than V1's
copy-then-collide. Only the 12T/512² case is slightly negative (−3–4%, from one
barrier's worth of phase splitting + shell collision). The **1183 MLUPS** for
1024²/12T/f32 is a new project CPU record (exceeding the V1 documented value of 1124).

### 3D D3Q19 (128³ periodic, TRT)

| grid | threads | prec | CpuScalar | CpuSimd | Simd/Scalar |
|---|---|---|---|---|---|
| 128³ | 1 | f32 | 17 | 52 | **3.0** |
| 128³ | 1 | f64 | 17 | 33 | **2.0** |
| 128³ | 12 | f32 | 139 | 260 | 1.9 |
| 128³ | 12 | f64 | 123 | 168 | 1.4 |

The target of "2× or more over CpuScalar" was achieved at 1T (3.0x / 2.0x). At 12T it
falls short, at 1.9x for f32 / 1.4x for f64. The cause has been isolated through
measurement: (a) double collision at band-edge slabs (+19% collision work with 12
bands), (b) coarse-grained bands cannot load-balance across heterogeneous P/E cores
(CpuScalar scales 7.5x with row-granularity work stealing, while the fused version
reaches only 5.2x). Both counter-measures — y-strip ringing and over-splitting the
bands — were implemented, measured, and rejected (−20%/−8% respectively: on this
machine the SLC absorbs the plane ring, so the re-collision of strips and degraded
prefetching outweigh the benefit of a smaller ring. Measured values are recorded in
the in-code doc). The remaining room for improvement is either "sharing band-edge
collisions (requires synchronization)" or "a fallback threshold to phase splitting,"
but even in absolute terms 12T f32 260 MLUPS is +87% over CpuScalar.

### Key points on kernel shape (recorded with measurements in the in-code doc)

- **D2Q9 = flat form**: within a cell, LLVM fully unrolls the 4 pairs and vectorizes
  the x loop (matches V1's `collide_span` down to the operator ordering). Blocking is
  counter-productive for D2Q9 (~1.2x loss).
- **D3Q19 = blocked form**: for the flat unroll of 9 pairs, LLVM gives up and falls
  back to scalar code (measured vec/scalar instruction ratio of 18 vs. 285). Shared
  quantities are placed in a `[T; 64]` scratch buffer, and each pair is swept with
  unit stride (lattice constants are loop-invariant broadcasts), making it
  vectorizable independent of Q.
- **In-place views**: passing separate src/dst views into the blocked kernel breaks
  vectorization due to alias checking (−30%). The blocked kernel uses a single view;
  only the flat kernel uses the copy-fused src=f/dst=ring form.
- **Ring persistence**: the per-band ring is carried in a `FusedScratch`
  (eliminating 92 MB/step@f64 of per-step allocation + zero-fill).

## Results of the Phase 9 kernel rework (relative comparison via same-time alternating runs)

**Historical record** (rework record for the V1 engine. V1, `bench_mlups`, and
`probe_state_hash` were removed with the V1 retirement on 2026-07-05 — preserved here
as a record of the numbers and method. The fused kernel itself was ported to the V2
`CpuSimd`, and the equivalence gate was taken over by `tests/backend_simd_equiv.rs` +
T13).

**Measurement caveat**: on the measurement machine that day, verification agents etc.
were running concurrently, and absolute values varied by up to 3× depending on time of
day. The table below is a relative comparison from **alternating runs of the old and
new binaries at the same time** (best-of-5, single-configuration mode of
`bench_mlups <prec> <n> <threads> <steps>`), which is robust to load variation in this
form.

| Configuration | Old (AoS scalar) | New (SoA fused kernel) | Ratio |
|---|---|---|---|
| f32 512² 1T | 38 | 220 | 5.7x |
| f32 512² 12T | 319 | 860 | 2.6x |
| f32 1024² 1T | 37 | 210 | 5.7x |
| f32 1024² 12T | 341 | **1124** | 3.2x |
| f64 512² 1T | 37 | 112 | 2.9x |
| f64 512² 12T | 303 | 506 | 1.6x |
| f64 1024² 1T | 37 | 120 | 3.2x |
| f64 1024² 12T | 323 | 602 | 1.8x |

The Phase 9 target of "2× or more for f32" was achieved (2.6–5.7×). The peak is
**~1100 MLUPS** at f32/1024²/12T (~3× the old document's peak of 381. An idle machine
would likely push this higher still). Grids under 16384 cells continue to run serially
automatically, as before (`PARALLEL_MIN_CELLS`).

## What was done (contribution by stage)

1. **SoA layout + span decomposition + paired-form kernel** (~4x at 1T): changed to
   plane-major storage of `f[q*(nx*ny) + y*nx + x]`. Each row is decomposed via a
   solid-run table into "spans containing no solid," completely eliminating branches
   from the inner loop → enables NEON auto-vectorization. Collision is strength-
   reduced into a form that directly computes the TRT pair quantities (ep/em) — an
   equivalent transformation. Because the expression form is x/y mirror-symmetric,
   strict equivariance is also preserved. The fluid span in streaming is simply a
   shifted slice copy.
2. **stream+moments fusion** (+15–30% at 12T, rayon barriers 3→2/step).
3. **Full single-pass fusion of collide+stream+moments** (1T +14–21%, 12T +31–95%):
   `step_band()` collides source rows just-in-time into a 3-row ring
   (cache-resident), pull-streams from the ring, and also writes moments in the same
   pass. The post-collision grid is not materialized to DRAM (memory traffic reduced
   ~40%, 1 barrier/step). Moments use double buffering (rho2/ux2/uy2) to preserve the
   convention that "collision reads the previous step's moments." The memory term of
   ConvectiveOutflow (the post-collision distribution from the previous step) is
   supplied bit-identically via a separate channel by capture_conv_stale().
4. **const-ifying the force-term loop** (FORCE/FF flags): fixed the issue where the
   per-cell Option branch on force_field was forcing the FORCE=true collision loop
   into fully scalar code. The collision loops for Poiseuille/gravity/multiphase are
   now all vectorized as well (confirmed via .4s/.2d instruction counts in the
   disassembly).

**Guarantee of equivalence**: alongside all-green `cargo test --workspace --release`
at every stage (57 tests, the heavy `--include-ignored` series also green), a
**bit-identical match with the old implementation** was confirmed via
`examples/probe_state_hash` (bit-exact hash of rho/u/probed force across 10
scenarios). With a single thread, even the probed force matches exactly; with
multiple threads, the fields match bit-for-bit while only the probed force can shift
in the final ulp due to reduction order. Any rework that changes the pass structure
or storage order must pass this probe's bit-exact match.

## Takeaways (= recommended settings)

1. **TRT is the same speed as BGK** (in the paired form, BGK is the special case
   ω+=ω− of TRT). TRT, which is superior in accuracy and stability, can always be the
   default. BGK is for education/comparison.
2. **f32 is ~1.9× f64** (2× vector width + half the bandwidth. The old implementation
   showed only a 5% difference because it was scalar-compute-bound). With deviation
   storage (f − w, introduced 2026-07-05), f32 precision is now verification-grade:
   - Momentum growth error under uniform force: 1.34e-3 → **2.8e-7** (~4800× improvement)
   - TGV L2 error (N=64): f32 7.1e-4 vs. f64 7.0e-4 (effectively equivalent)
   → **f32 is recommended for typical 2D simulations.** Use f64 for long-time
   integration, extremely small gradients, or verification computations.
3. Since single-thread went from ~35 → ~210 MLUPS, **the WASM GUI's resolution
   ceiling can potentially be raised** (requires actual WASM measurement; a similar
   gain is expected with SIMD128 enabled. Even in the automatic-serial region below
   `PARALLEL_MIN_CELLS`, still ~5x over the old version).
4. The smaller gain for 12T f64 (1.6–1.8x) versus f32 reflects approaching
   bandwidth-bound territory (shared-machine load may also be a factor; re-measurement
   while idle is recommended).

## Known remaining optimization opportunities (Phase 9 continuation / notes for V2)

- **Explicit SIMD (std::simd / wide)**: currently relies on auto-vectorization alone.
  In theoretical FLOPs there is still ~2x of headroom at f32 1T, but if this is
  pursued it should be in the V2 CpuSimd backend (`step_band()` is the reference
  kernel targeted for porting. ARCHITECTURE_V2.md §2.4). → **Porting completed in
  M-E** (see the "V2 CpuSimd backend" section above). 2D exceeds V1, and 3D is also
  fused with the same structure. Explicit SIMD remains untouched (reached via
  auto-vectorization alone).
- **In-place streaming such as Esoteric-pull**: room to halve memory and roughly
  double effective bandwidth by eliminating ftmp. A candidate for introduction in V2
  M-E.
- **Band-split tuning**: with the number of bands fixed to the number of threads,
  imbalance between bands remains for geometries where obstacles are unevenly
  distributed (work stealing is not possible within a band, since it is sequential).
  Double collision at halo rows is ~2/band-row-count (~2% at 12T/1024²).
- **GPU (wgpu, Phase 9c measured 2026-07-05)**: on the same physical machine, M5 Max
  reaches 7,584 MLUPS (1024²) / 6,975 MLUPS (2048²), verification L∞ 7e-6. See
  [GPU_EVALUATION.md](GPU_EVALUATION.md) for details and the adoption decision.
  Planned for integration into V2's Backend trait.
