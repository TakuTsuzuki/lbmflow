# MPI Distributed Guide

The lbm-core `mpi` feature provides an MPI implementation of `HaloExchange`
(`dist::MpiExchange`) and a 1-rank = 1-subdomain driver (`dist::MpiSolver`).
Design: docs/ARCHITECTURE_V2.md §2.3.

## Status

- **Verified**: multi-rank within a single node (Open MPI 5.0.9 / arm64
  macOS via shared-memory BTL). T13-MPI: distributed execution is
  **bit-identical** to single-rank in all cases (fields Δ = 0.0; diagnostics
  differ only by f64 partial-sum recombination, atol+rtol 1e-11).
- **Weak scaling** (single node, shmem BTL, not a real interconnect):
  99.4% / 97.0% at 2 / 4 ranks — meets R3's ≥85% local pass line.
  n=8 (73.2%) drops due to M5 Max E/P core heterogeneity + lockstep jitter,
  not communication volume (~50 KB/step/rank). True weak scaling awaits
  cluster measurement.
- **Supported after BCFD-100..101**:
  - Large-grid local geometry construction via `MpiSolver::new_local`, where
    ranks evaluate solid, wall-velocity, and material callbacks only for owned
    cells. The older `MpiSolver::new` global-array constructor remains for
    small compatibility tests.
  - Per-rank binary field slabs plus a rank-0 manifest for velocity, `phi`,
    named scalar concentrations such as oxygen `C_L`, shear rate, and gas
    holdup. Use `read_parallel_field` to reconstruct a compact global vector
    for small validation comparisons.
- **Not yet supported**:
  - `mpi` + `gpu`: builds, but `MpiSolver` assumes a `CpuScalar`-family
    `SoaFields` backend. Device-resident halo + GPUDirect is M-E or later.
  - Communication/computation overlap: `exchange_f` is synchronous per
    axis phase. The two-pass streaming seam exists; switching to
    Isend → interior → wait → boundary is M-E.
  - Parallel VTK / HDF5. BCFD-101 writes raw per-rank slabs; richer container
    formats remain future work.
  - Multi-node measurements: **ME-3 = RED**, blocked on cluster spend
    confirm (see [CLUSTER_OPTIONS.md](CLUSTER_OPTIONS.md)).

## Build

rsmpi (`mpi` 0.8) probes for `mpicc` at build time. **A native arm64
MPI must be first in `PATH`** (an x86_64 MPI in `/usr/local` breaks the
probe or forces Rosetta at launch; verify with `file $(which mpirun)`).

```bash
export PATH=$HOME/.local/openmpi/bin:$PATH
file $(which mpirun)   # confirm Mach-O 64-bit arm64
cargo build -p lbm-core --release --features mpi
```

Open MPI source build (5.0.9 / arm64, Fortran disabled, ~15 min):

```bash
mkdir -p ~/.local/src && cd ~/.local/src
curl -sL https://download.open-mpi.org/release/open-mpi/v5.0/openmpi-5.0.9.tar.bz2 | tar xj
cd openmpi-5.0.9
./configure --prefix=$HOME/.local/openmpi CC=clang CXX=clang++ --disable-mpi-fortran
make -j$(sysctl -n hw.ncpu) && make install
```

The default build has zero rsmpi dependency; `cargo test --workspace`
continues to pass without an MPI environment.

## Run

```bash
# T13-MPI verification (-n 1,2,4 for 2D; -n 8 for 3D TGV 2×2×2)
./scripts/test_mpi.sh

# Weak-scaling sweep (512²/rank, ranks {1,2,4,8})
./scripts/bench_mpi.sh          # override: LOCAL=512 STEPS=200 RANKS="1 2 4 8"
```

Minimal API (1 rank = 1 part; Cartesian decomposition must match rank count):

```rust
use lbm_core::dist::MpiSolver;
use lbm_core::prelude::*;

let universe = mpi::initialize().unwrap();
let world = universe.world();
let spec = GlobalSpec::<f64> { dims: [1024, 512, 1], ..Default::default() };
let mut s: MpiSolver<D2Q9, f64, CpuScalar> =
    MpiSolver::new(&world, &spec, &[], &[], [world.size() as usize, 1, 1],
                   CpuScalar::default());
s.init_with(|x, y, _| (1.0, [0.0, 0.0, 0.0]));
s.run(1000);                       // step/diagnostics/gather are collective
let mass = s.total_mass();         // Allreduce (same on all ranks)
let rho  = s.gather_rho();         // Some(field) on rank 0 only
drop(s);                           // drop the duplicated comm before finalize
```

Large-grid API (BCFD-100):

```rust
use lbm_core::prelude::*;

let mut s: MpiSolver<D3Q19, f64, CpuScalar> = MpiSolver::new_local(
    &world,
    &spec,
    [2, 2, 2],
    CpuScalar::default(),
    |g| tank_mask(g.x, g.y, g.z),
    |g| tank_wall_velocity(g.x, g.y, g.z),
    |_g| None,
    |_g, _name| None,
);
let local_bytes = s.mem_estimate_bytes();
let _manifest = s.write_velocity_slabs("out/fields")?; // Some(...) on rank 0
drop(s);
```

**Collective contract.** `step` / `init_with` / `update_shan_chen_force` /
diagnostics (`total_mass`, `total_momentum`, `nonfinite_count`) / `gather_*`
and mask edits must be called by all ranks in the same order. `set_solid`
must be called by every rank with the same coordinate sequence (owner
stores; non-owners only reserve halo re-exchange). Drop `MpiSolver`
**before** the `Universe` (i.e. before `MPI_Finalize`).

## Exchange protocol

- Identical x → y → z phase order and pack/unpack as the in-process
  exchange; both call the shared helper in `halo.rs`. Each phase transfers
  the halo-inclusive layer of preceding axes (only the 6 face links carry
  data even with MPI).
- Per axis phase: Irecv → Isend → wait-all → unpack. Tag for a message to
  face `F` is `base + F.index()` (f: 100, ψ: 200, mask: 300/400, gather:
  500). Periodic axis with decomp=1 self-wraps locally, bypassing MPI.
- Transferred payload is raw scalar bytes, so the partitioned field is
  bit-identical to the single-rank field; only diagnostics accumulate a
  partial-sum → Allreduce f64 recombination difference.
- Rank placement uses `solver::partition`'s Cartesian layout
  (part id = `(pz·dy + py)·dx + px` = rank). `MPI_Cart` not used.

## T13-MPI evidence (2026-07-05, M5 Max / Open MPI 5.0.9)

| Case | -n | decomp | field max\|Δ\| | diagnostic max\|Δ\| |
|---|---|---|---|---|
| 2D TGV 96×64 | 1/2/4 | 1×1 / 2×1 / 2×2 | **0.0** | ≤ 3.3e-14 |
| Cavity 64×64 (lid crosses seam) | 1/2/4 | " | **0.0** | ≤ 2.3e-14 |
| Cylinder + force probe (on seam) + parabolic inflow | 1/2/4 | " | **0.0** | ≤ 9.1e-13 |
| Shan-Chen droplet (ψ via `exchange_scalar`) | 1/2/4 | " | **0.0** | ≤ 4.5e-11 * |
| 3D TGV 24³ (D3Q19) | 8 | 2×2×2 | **0.0** | ≤ 4.6e-15 |

\* Relative to total_mass ≈ 1.5e3 (~3e-14 relative); pass judged by
`atol + rtol·|ref|`, both 1e-11.

## Weak scaling (single-node shmem, 2026-07-05)

512²/rank, D2Q9 f64 TGV, serial backend per rank, 200-step measurement,
20-step warmup:

| ranks | time | total MLUPS | MLUPS/rank | efficiency |
|---|---|---|---|---|
| 1 | 1.304 s | 40.2 | 40.2 | 100% |
| 2 | 1.313 s | 79.9 | 40.0 | 99.4% |
| 4 | 1.346 s | 155.9 | 39.0 | 97.0% |
| 8 | 1.781 s | 235.5 | 29.4 | 73.2% |

n=8 drop is hardware (6 E + 12 P heterogeneous cores + bandwidth
contention: 8 independent single-rank jobs only reach 84% of ideal) ×
lockstep jitter, not MPI volume. **True weak scaling requires a real
interconnect** — see the cluster checklist below.

## Cluster measurement checklist (R3 completion)

1. **Weak scaling**: 3D 128³/rank (D3Q19) from 1→64 ranks. R3 pass line:
   ≥80% at 64 ranks. Also 2D 512² for continuity with the table above.
2. **Strong scaling**: fixed 1024³ from 8→512 ranks.
3. **Communication/computation ratio**: `MPI_T` profile or mpiP;
   if `exchange_f` > 10%, bring forward two-pass overlap (M-E).
4. **BTL/MTL**: UCX/OFI selection, eager/rendezvous thresholds (layer
   messages are 10–200 KB), tag-matching contention.
5. **Rank × thread hybrid sweep** per node (reuse
   `CpuScalar::parallel_min_cells`).
6. **Placement/binding**: `--map-by` / `--bind-to`, NUMA pack/unpack.
7. **Cross-implementation correctness**: `scripts/test_mpi.sh` on MPICH /
   Cray MPICH as well (rsmpi supports both; bit-identity should be
   implementation-independent).
8. **Diagnostics/output cost at scale**: Allreduce and rank-0 gather
   scaling; decision material for HDF5 / ADIOS2 parallel I/O migration.

## Known pitfalls

- **Coexistence with x86_64 MPI**: Homebrew (x86_64) mpicc in `/usr/local`
  ahead in `PATH` makes rsmpi's probe pick x86_64 flags — link fails or
  launch goes through Rosetta. Keep arm64 first.
- **`MPI_Comm_free` after `MPI_Finalize`**: `MpiSolver` / `MpiExchange`
  release their duplicated communicator on `Drop`. Drop the solver **before**
  the `mpi::initialize()` `Universe`.
- **Collectiveness of mask edits**: calling `set_solid` only on the
  owning rank drifts the collective call count and deadlocks.
  `MpiSolver::set_solid` fixes this; direct `Solver` users must call
  `mark_masks_dirty()` on every rank.
- **Probe Allreduce timing**: `probed_force` Allreduces every step. Omitted
  when no probe is configured (to keep the benchmark loop clean).

## `run_guarded` (R-Phase 1 A-9)

`MpiSolver::run_guarded(steps, check_every)` is collective: every rank
invokes with identical arguments. Each check is a 2-double Allreduce over
mass partials (NaN propagates through the sum); the divergence branch
(`Err(Diverged { step })`) is taken uniformly on all ranks, so there is
no divergence-induced deadlock. Overhead at 512²/rank with `check_every=100`
is <0.5%. Trajectory is bit-identical to plain `run`.
