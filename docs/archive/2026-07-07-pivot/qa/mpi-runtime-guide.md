# MPI Runtime Validation Guide

Lifecycle: living (updated in place). Owner: V&V / ME-3 runtime lane.

`docs/qa/VV_TRACEABILITY.md` marks `T13-MPI` as `BENCH-PENDING` because
runtime execution needs a native MPI launcher and is not covered by
`cargo test --workspace --release`. The ready-to-run harness is:

```bash
scripts/qa/mpi_runtime_check.sh --dry-run
scripts/qa/mpi_runtime_check.sh
```

The dry run prints the detected MPI environment and the commands that would be
executed. It does not call `mpirun`.

## What T13-MPI Validates

T13-MPI is a distributed code-verification gate for partition invariance:

- The MPI run uses `crates/lbm-core/examples/mpi_t13.rs`.
- Each rank owns one subdomain through `dist::MpiSolver`.
- Rank 0 also runs the monolithic baseline for the same scenario.
- Gathered fields (`rho`, velocity, and every population plane) must match the
  monolithic reference with max field difference `0.0`.
- Diagnostics (`total_mass`, momentum, probed force, nonfinite count) are
  checked with the tolerance in `docs/VALIDATION.md` because f64 partial sums
  are recombined in a different order across ranks.

This is not an absolute physics validation. It proves that the MPI halo,
gather, scalar exchange, seam probes, and distributed run contract reproduce
the stored T13 scenario reference bit-identically for fields. The reference is
the same class of monolithic-vs-partition comparison used by
`crates/lbm-core/tests/t13_split_invariance.rs`, exercised through the MPI
runtime path.

## Harness Behavior

`scripts/qa/mpi_runtime_check.sh` performs these steps:

1. Reports MPI environment: `MPI_HOME`, `mpirun`, `mpicc`, `mpirun --version`,
   rank counts, and artifact directory.
2. If MPI is missing, prints setup commands for local Open MPI, Linux/AWS, and
   cluster module environments.
3. Verifies launcher startup with `mpirun -n 2 echo hi`.
4. Runs `scripts/test_mpi.sh` with `MPI_RANKS="1 2 4 8"` by default.
5. Checks the T13 output for `ALL PASS` at every requested rank count and for
   bit-identical field comparisons in every PASS line.
6. Runs small weak/strong `bench_mpi` scaling smoke cases when `RUN_SCALING=1`
   and writes CSV summaries with `scripts/qa/aggregate_mpi_results.sh`.

Artifacts are written under `target/qa/mpi-runtime/<timestamp>/` unless
`LOG_DIR` is set. Important files:

- `hello.log` - launcher preflight.
- `t13-mpi.log` - T13-MPI correctness output.
- `scaling.log` - raw weak/strong benchmark output.
- `scaling-summary.csv` - aggregate performance table.
- `scaling-ranks.csv` - rank placement / environment evidence.

Useful overrides:

```bash
RANK_COUNTS="1 2 4 8" scripts/qa/mpi_runtime_check.sh
RUN_SCALING=0 scripts/qa/mpi_runtime_check.sh
SCALING_RANKS="1 2 4 8 16 32 64" WEAK_LOCAL_EDGE=128 \
  STRONG_GLOBAL_EDGE=1024 SCALING_STEPS=200 SCALING_WARMUP=20 \
  scripts/qa/mpi_runtime_check.sh
```

## Local Open MPI

For an Apple silicon workstation, keep a native arm64 MPI first in `PATH`.
The default harness path is `$HOME/.local/openmpi/bin`.

```bash
mkdir -p ~/.local/src && cd ~/.local/src
curl -sL https://download.open-mpi.org/release/open-mpi/v5.0/openmpi-5.0.9.tar.bz2 | tar xj
cd openmpi-5.0.9
./configure --prefix=$HOME/.local/openmpi CC=clang CXX=clang++ --disable-mpi-fortran
make -j$(sysctl -n hw.ncpu) && make install
export PATH=$HOME/.local/openmpi/bin:$PATH
scripts/qa/mpi_runtime_check.sh
```

If Homebrew or another MPI appears earlier in `PATH`, verify architecture with:

```bash
file "$(command -v mpirun)"
mpirun --version
```

## AWS hpc7g

AWS hpc7g is the recommended fast ME-3 cluster path in
`docs/CLUSTER_OPTIONS.md`.

Build and run on the cluster host or inside the job allocation:

```bash
sudo apt-get update
sudo apt-get install -y openmpi-bin libopenmpi-dev
export MPIRUN=mpirun MPICC=mpicc
scripts/qa/mpi_runtime_check.sh --dry-run
scripts/qa/mpi_runtime_check.sh
```

For the full ME-3 measurement, use larger benchmark settings after the T13
gate passes:

```bash
SCALING_RANKS="1 2 4 8 16 32 64" \
WEAK_LOCAL_EDGE=128 STRONG_GLOBAL_EDGE=1024 \
SCALING_STEPS=200 SCALING_WARMUP=20 \
scripts/qa/mpi_runtime_check.sh
```

Capture placement and network evidence from the scheduler or launcher output
with the resulting `scaling-ranks.csv`.

## Fugaku

On Fugaku, run inside a normal job allocation and load the site MPI/compiler
environment before building. The exact module names are site-managed; use the
local wrappers exposed by the allocation.

```bash
module avail mpi
module load <site-mpi-module>
export MPIRUN=<site-mpi-launcher>
export MPICC=<site-mpicc-wrapper>
scripts/qa/mpi_runtime_check.sh --dry-run
scripts/qa/mpi_runtime_check.sh
```

If the site uses a launcher other than `mpirun`, set `MPIRUN` to that wrapper
and keep `MPICC` pointed at the wrapper compiler that `rsmpi` should probe.
Fugaku results are valuable even at small rank counts because they exercise a
non-local MPI stack and interconnect.

## Known Issue: Sandboxed macOS Sessions

Local sandboxed sessions on macOS can block MPI socket binding. The observed
failure is `Operation not permitted` from the Open MPI / PRTE launch path before
the Rust example starts. Treat that as an environment block, not as T13-MPI
evidence.

Production MPI runs happen off-sandbox by the PM or on a real cluster. Keep the
traceability status `BENCH-PENDING` until a real `t13-mpi.log` from the harness
shows all requested rank counts passing.

## MPI-Only ANOM Format

Log MPI-only findings in `docs/qa/anomaly-log.md` with enough detail for a
cluster rerun. Use this shape:

```markdown
### ANOM-P4-NNN - MPI <short finding title>

- Date: 2026-07-07.
- Host / allocation: <cluster, node type, node count, scheduler job id>.
- MPI stack: <Open MPI / MPICH / Cray MPICH / Fujitsu MPI version>.
- Command: `scripts/qa/mpi_runtime_check.sh ...`
- Artifact directory: `target/qa/mpi-runtime/<timestamp>/` or durable path.
- T13 status: PASS / FAIL / NOT-RUN, rank counts, first failing rank.
- Scaling status: weak/strong summary, first failing rank or efficiency drop.
- Observed behavior: field mismatch, diagnostic drift, deadlock, launcher
  failure, placement issue, or performance regression.
- Expected behavior: T13 field max difference `0.0`; diagnostics within
  `docs/VALIDATION.md` tolerance; scaling trend consistent with ME-3 target.
- Routing: core-engine / runtime-infra / cluster-config / spec.
- Repro package: `hello.log`, `t13-mpi.log`, `scaling.log`,
  `scaling-summary.csv`, `scaling-ranks.csv`.
```

Do not file a core-engine MPI ANOM from a sandbox socket-bind failure alone.
File it as an environment block only if it prevents scheduled off-sandbox
execution.
