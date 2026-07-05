# MPI Distributed Guide (M-D, 2026-07-05)

The lbm-core `mpi` feature provides an MPI implementation of HaloExchange
(`dist::MpiExchange`) and a 1-rank = 1-subdomain driver (`dist::MpiSolver`).
The design corresponds to docs/ARCHITECTURE_V2.md §2.3 / docs/HPC_SCALING.md
staged plan 3.

## Current scope (honest current state)

- **Verified**: Multi-rank within a single node (Open MPI 5.0.9 / arm64 macOS,
  via shared-memory BTL). T13-MPI confirmed distributed execution ≡
  single-rank execution (fields are bit-identical, diagnostics show only
  f64 recombination differences. See below).
- **Not yet supported**:
  - Combined use with the GPU backend (`--features mpi,gpu` builds, but
    `MpiSolver` assumes a `CpuScalar`-family `SoaFields` backend.
    Device-resident halo transfer and GPUDirect are M-E or later).
  - Overlap of communication and computation (`exchange_f` completes
    synchronously within each axis phase. The two-pass streaming seam
    already exists in the Solver, so switching to Isend issue → internal
    computation → wait → boundary computation is an M-E candidate).
  - Parallel I/O (only full-field reconstruction via rank-0 gather.
    Parallel VTK/HDF5 not yet started).
  - Multi-node measurements (awaiting cluster access. See §Cluster
    measurement checklist).

## Build

rsmpi (crate `mpi` 0.8) looks for `mpicc` at build time. **A native arm64
MPI must be first in PATH** (in environments where an x86_64 MPI lives in
/usr/local, PATH ordering causes trouble. Confirm arm64 with
`file $(which mpirun)`).

```bash
# Use a source-built Open MPI (e.g. $HOME/.local/openmpi)
export PATH=$HOME/.local/openmpi/bin:$PATH
file $(which mpirun)   # → confirm Mach-O 64-bit executable arm64

cargo build -p lbm-core --release --features mpi
```

Open MPI source build steps (reference: 5.0.9 / arm64, Fortran disabled,
~15 minutes):

```bash
mkdir -p ~/.local/src && cd ~/.local/src
curl -sL https://download.open-mpi.org/release/open-mpi/v5.0/openmpi-5.0.9.tar.bz2 | tar xj
cd openmpi-5.0.9
./configure --prefix=$HOME/.local/openmpi CC=clang CXX=clang++ --disable-mpi-fortran
make -j$(sysctl -n hw.ncpu) && make install
```

The default build (no feature) has zero dependency on rsmpi.
`cargo test --workspace` continues to pass without an MPI environment,
as before.

## Run

```bash
# T13-MPI verification (-n 1,2,4: 2D 4 cases / -n 8: 3D TGV 2x2x2. nonzero exit = failure)
./scripts/test_mpi.sh

# Weak scaling (512^2 per rank, ranks {1,2,4,8}, table output)
./scripts/bench_mpi.sh          # tune via LOCAL=512 STEPS=200 RANKS="1 2 4 8"
```

Minimal API example (1 rank = 1 part. Cartesian decomposition must match
the rank count):

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
s.run(1000);                       // step/diagnostics/gather are all collective
let mass = s.total_mass();         // Allreduce (same value on all ranks)
let rho = s.gather_rho();          // Some(full field) only on rank 0
drop(s);                           // release the duplicated communicator before finalize
```

**Collective contract**: `step` / `init_with` / `update_shan_chen_force` /
diagnostics (`total_mass` / `total_momentum` / `nonfinite_count`) /
`gather_*` / mask edits must be called by all ranks in the same order.
`set_solid` must be **called by all ranks with the same coordinate
sequence** (the owning rank stores it, other ranks only reserve the halo
re-exchange). Since `MpiSolver` holds a duplicated communicator, drop the
solver **before** dropping the `mpi::initialize()` Universe (i.e. before
MPI_Finalize).

## Exchange protocol (implementation notes)

- **Identical x → y → z phase order and identical pack/unpack** as
  InProcess (both implementations call the shared helper in `halo.rs`).
  Corners/edges are forwarded via face adjacency: each phase transfers the
  halo-inclusive layer of preceding axes (only the 6 face links are used
  even with MPI).
- For each axis phase, the layers on both of the two faces are issued as
  Irecv → Isend → wait for all completions → unpack. The tag for a
  message addressed to receiving face `F` is `base + F.index()`
  (f: 100, ψ: 200, mask: 300/400, gather: 500). For a periodic axis with
  decomp=1, self-wrap is a local copy that bypasses MPI.
- Since the transfer is raw bytes of a scalar type (reversible for both
  f64/f32), the partitioned-execution field is **bit-identical** to the
  single-rank execution. Only diagnostics allow a difference from
  rank partial-sum → Allreduce f64 recombination (T13 convention:
  atol + rtol, 1e-11).
- Rank placement follows `solver::partition`'s Cartesian decomposition
  as-is (part id = `(pz·dy+py)·dx+px` = rank). MPI_Cart is not used.

## T13-MPI measurements (2026-07-05, M5 Max / Open MPI 5.0.9)

| Case | -n | decomp | field max\|Δ\| | diagnostic max\|Δ\| |
|---|---|---|---|---|
| 2D TGV 96×64 | 1/2/4 | 1×1 / 2×1 / 2×2 | **0.0** (bit-identical) | ≤ 3.3e-14 |
| Cavity 64×64 (lid crosses the seam) | 1/2/4 | same as above | **0.0** | ≤ 2.3e-14 |
| Cylinder + force probe (on the seam) + parabolic inflow | 1/2/4 | same as above | **0.0** | ≤ 9.1e-13 |
| Shan-Chen droplet (ψ via exchange_scalar, 2×2 corner) | 1/2/4 | same as above | **0.0** | ≤ 4.5e-11* |
| 3D TGV 24³ (D3Q19) | 8 | 2×2×2 | **0.0** | ≤ 4.6e-15 |

\* The droplet's diagnostic difference is the recombination difference
relative to total_mass ≈ 1.5e3 (relative ~3e-14). Pass/fail is judged by
`atol + rtol·|ref|` (both 1e-11) and all cases PASS.

## Weak scaling (single-node measurement, 2026-07-05)

512² per rank (D2Q9 f64 TGV, **serial** backend within each rank,
decomp [n,1,1], 200-step measurement, 20-step warmup):

| ranks | time | total MLUPS | MLUPS/rank | efficiency |
|---|---|---|---|---|
| 1 | 1.304 s | 40.2 | 40.2 | 100% |
| 2 | 1.313 s | 79.9 | 40.0 | 99.4% |
| 4 | 1.346 s | 155.9 | 39.0 | 97.0% |
| 8 | 1.781 s | 235.5 | 29.4 | 73.2% |

**How to read this (important)**: this measurement goes through Open MPI's
**shared memory**, not a real interconnect measurement. Furthermore, the
measurement machine (M5 Max) has a heterogeneous core configuration of
6 Efficiency + 12 Performance cores, and even a control experiment
(8 independent single-rank jobs run concurrently with zero communication)
only reaches 33.7 MLUPS/rank (= equivalent to 84%). In other words, the
breakdown of n=8's 73.2% is "the ceiling from hardware (heterogeneous
cores + bandwidth contention) at 84%" × "the remaining ~87.5% from
MPI-ization (mainly jitter coupling from lockstep synchronization to the
slowest rank; the communication volume itself is ~50 KB/step/rank and
negligible)." n≤4, which stays within the homogeneous cores, meets R3's
local pass line of ≥85% (97-99%). **True weak scaling requires cluster
measurements** (see below).

## Checklist of measurements to do on a cluster (R3 completion criteria)

1. **Full weak-scaling measurement**: 3D 128³ per rank (D3Q19) from
   1→64 ranks (within-node → multi-node). R3 pass line: efficiency
   ≥80% at 64 ranks. Also do a 2D 512² version for comparison
   (to connect with this guide's single-node table).
2. **Strong scaling**: fixed 1024³ from 8→512 ranks (the point where
   surface-area/volume ratio degrades).
3. **Measured communication/computation ratio**: the share occupied by
   exchange_f via an MPI_T profile or mpiP. If >10%, bring forward the
   two-pass overlap (M-E).
4. **Inter-node BTL/MTL confirmation**: UCX/OFI selection, eager/rendezvous
   thresholds (layer messages are in the ~10-200 KB band), and presence
   of tag-matching contention.
5. **Optimal rank × thread hybrid point**: this guide's measurements use
   serial execution within each rank. Sweep a grid of
   "rank count × rayon thread count" per node (the
   `CpuScalar::parallel_min_cells` threshold can be reused as-is).
6. **Process placement/binding**: `--map-by`/`--bind-to` (not possible on
   macOS, required on Linux clusters). Also check layer pack/unpack
   bandwidth across NUMA nodes.
7. **Re-confirm correctness**: make scripts/test_mpi.sh pass fully on the
   cluster's MPI implementation too (MPICH/Cray MPICH in addition to
   Open MPI) (rsmpi supports both. The bit-identity requirement should be
   implementation-independent).
8. **Diagnostic cost at scale**: the scaling limits of Allreduce
   (diagnostics) and rank-0 gather (output). For output, gather material
   for a decision on migrating to parallel I/O (HDF5/ADIOS2-family).

## Known pitfalls (hit this time)

- **Coexistence with x86_64 MPI**: if the Homebrew (x86_64) mpicc in
  /usr/local is first in PATH, rsmpi's probe picks up x86_64 flags and
  either fails to link, or launching breaks via a Rosetta-mediated
  mpirun. Always keep the arm64 version first in PATH.
- **MPI_Comm_free after MPI_Finalize**: `MpiSolver` (and `MpiExchange`)
  release the duplicated communicator on Drop. Design the scope so it is
  dropped **before** the `mpi::initialize()` Universe (see examples).
- **Collectiveness of mask edits**: calling `set_solid` only on the
  owning rank causes the halo re-exchange (collective) call count to
  drift between ranks, causing a deadlock. `MpiSolver::set_solid` is
  designed so that non-owning ranks only set the dirty mark. When using
  `Solver` directly, call `mark_masks_dirty()` on all ranks.
- **Allreduce timing for probes**: probed_force is an Allreduce on every
  step. It is omitted when no probe is configured (so as not to add
  extraneous collectives to the benchmark).
