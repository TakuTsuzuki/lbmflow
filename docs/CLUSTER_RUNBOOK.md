# ME-3 Cluster Campaign Runbook

This runbook prepares the ME-3 cluster campaign without launching any cloud
resources. It cross-references the campaign rationale in
[CLUSTER_OPTIONS.md](CLUSTER_OPTIONS.md), the MPI implementation notes in
[MPI_GUIDE.md](MPI_GUIDE.md), and the scaling target in
[paper/claims-ledger.md](paper/claims-ledger.md): 64-rank weak scaling >=80%.

## Preflight

On the target cluster, build with the MPI feature using the cluster MPI wrapper:

```bash
module purge
module load openmpi   # or the site MPI module
which mpicc mpirun
mpirun --version
cargo build -p lbm-core --release --features mpi --example bench_mpi --example mpi_t13
```

Confirm the MPI thread level before any performance run. Every `bench_mpi`
mode prints:

```text
MPI_THREAD requested=Funneled provided=Funneled
```

`provided` must be `Funneled`, `Serialized`, or `Multiple`. `Single` is a stop
condition if `--parallel`/rayon hybrid measurements are planned.

Verify NUMA and affinity placement before timing:

```bash
mpirun --report-bindings --map-by ppr:64:node:PE=1 --bind-to core \
  -n 64 ./target/release/examples/bench_mpi --mode placement
```

Expected shape:

```text
RANK_RESULT mode=placement ... rank=0 ... hostname=... affinity=...
RESULT mode=placement ranks=64 hostname=... ompi_pml=... fi_provider=...
```

On Slurm sites, use the site equivalent, for example:

```bash
srun --mpi=pmix --cpu-bind=cores -n 64 ./target/release/examples/bench_mpi --mode placement
```

## The 8 Measurements

Common variables:

```bash
BIN=./target/release/examples/bench_mpi
MPI="mpirun --report-bindings"
```

1. Full weak scaling, primary 3D D3Q19 target:

```bash
$MPI -n 1  $BIN --mode weak3d --local-edge 128 --steps 200 --warmup 20 --decomp 1x1x1
$MPI -n 2  $BIN --mode weak3d --local-edge 128 --steps 200 --warmup 20 --decomp 2x1x1
$MPI -n 4  $BIN --mode weak3d --local-edge 128 --steps 200 --warmup 20 --decomp 2x2x1
$MPI -n 8  $BIN --mode weak3d --local-edge 128 --steps 200 --warmup 20 --decomp 2x2x2
$MPI -n 16 $BIN --mode weak3d --local-edge 128 --steps 200 --warmup 20 --decomp 4x2x2
$MPI -n 32 $BIN --mode weak3d --local-edge 128 --steps 200 --warmup 20 --decomp 4x4x2
$MPI -n 64 $BIN --mode weak3d --local-edge 128 --steps 200 --warmup 20 --decomp 4x4x4
```

Expected `RESULT`: `mode=weak3d lattice=D3Q19 ranks=... decomp=... global=...
time_s=... mlups_total=... mlups_per_rank=... nonfinite=0`.

Also run the 2D comparison from `MPI_GUIDE.md`:

```bash
$MPI -n 64 $BIN --mode weak2d --local-edge 512 --steps 200 --warmup 20 --decomp 64x1x1
```

2. Strong scaling, fixed 1024^3 D3Q19:

```bash
$MPI -n 8   $BIN --mode strong3d --global-edge 1024 --steps 200 --warmup 20 --decomp 2x2x2
$MPI -n 16  $BIN --mode strong3d --global-edge 1024 --steps 200 --warmup 20 --decomp 4x2x2
$MPI -n 32  $BIN --mode strong3d --global-edge 1024 --steps 200 --warmup 20 --decomp 4x4x2
$MPI -n 64  $BIN --mode strong3d --global-edge 1024 --steps 200 --warmup 20 --decomp 4x4x4
$MPI -n 128 $BIN --mode strong3d --global-edge 1024 --steps 200 --warmup 20 --decomp 8x4x4
$MPI -n 256 $BIN --mode strong3d --global-edge 1024 --steps 200 --warmup 20 --decomp 8x8x4
$MPI -n 512 $BIN --mode strong3d --global-edge 1024 --steps 200 --warmup 20 --decomp 8x8x8
```

3. Communication/computation ratio:

Run the same 64-rank weak3d command under the site profiler. Example mpiP hook:

```bash
LD_PRELOAD=/path/to/libmpiP.so $MPI -n 64 $BIN \
  --mode weak3d --local-edge 128 --steps 200 --warmup 20 --decomp 4x4x4
```

Expected artifacts: normal `RESULT` line plus profiler report with `MPI_Isend`,
`MPI_Irecv`, `MPI_Waitall`, and collective percentages. If MPI time exceeds 10%,
bring forward the M-E overlap work.

4. Inter-node BTL/MTL confirmation:

```bash
FI_PROVIDER=efa $MPI -n 64 $BIN --mode placement
FI_PROVIDER=tcp $MPI -n 64 $BIN --mode placement
```

For UCX/IB clusters, use the site equivalent:

```bash
UCX_TLS=rc,sm,self $MPI -n 64 $BIN --mode placement
```

Expected output: `RANK_RESULT` rows showing host distribution and a `RESULT
mode=placement` row showing the selected environment. Record launcher verbose
logs for eager/rendezvous thresholds.

5. Optimal rank x thread hybrid point:

```bash
RAYON_NUM_THREADS=1  $MPI -n 64 $BIN --mode weak3d --parallel --local-edge 128 --steps 200 --warmup 20 --decomp 4x4x4
RAYON_NUM_THREADS=2  $MPI -n 32 $BIN --mode weak3d --parallel --local-edge 128 --steps 200 --warmup 20 --decomp 4x4x2
RAYON_NUM_THREADS=4  $MPI -n 16 $BIN --mode weak3d --parallel --local-edge 128 --steps 200 --warmup 20 --decomp 4x2x2
RAYON_NUM_THREADS=8  $MPI -n 8  $BIN --mode weak3d --parallel --local-edge 128 --steps 200 --warmup 20 --decomp 2x2x2
```

Expected output: `parallel=true rayon_threads=<N>`, with the best point selected
by `mlups_total` and efficiency.

6. Process placement/binding and NUMA:

```bash
$MPI --map-by ppr:64:node:PE=1 --bind-to core -n 64 $BIN --mode placement
numactl --hardware
lstopo-no-graphics || true
```

Expected output: every rank has an affinity token and host name. Save launcher
binding reports.

7. Correctness under the cluster MPI:

```bash
./scripts/test_mpi.sh
$MPI -n 64 $BIN --mode correctness --local-edge 64 --steps 20 --warmup 5 --decomp 4x4x4
```

Expected output: `test_mpi.sh: ALL PASS` and `RESULT mode=correctness ...
nonfinite=0`. If the cluster provides MPICH/Cray/Fujitsu MPI, rebuild after
switching modules and repeat.

8. Diagnostic Allreduce and rank-0 gather cost:

```bash
$MPI -n 64 $BIN --mode diagnostics --local-edge 128 --steps 200 --warmup 20 \
  --diagnostics-every 1 --gather-rho --decomp 4x4x4
```

Expected output: `diag_calls=200 diag_time_s=... gather_rho=true
gather_time_s=...`. Use the trend to decide whether parallel I/O must move
ahead of schedule.

## Weak Scaling Procedure

Use the 3D commands in item 1 at ranks `1, 2, 4, 8, 16, 32, 64`. The fixed
per-rank grid is 128^3 D3Q19 f64 cells. Efficiency is:

```text
efficiency(rank) = (MLUPS_per_rank(rank) / MLUPS_per_rank(1)) * 100
```

The claims-ledger target is `efficiency(64) >= 80%`.

## Strong Scaling Procedure

Use the item 2 commands at ranks `8, 16, 32, 64, 128, 256, 512`. The fixed total
grid is 1024^3 D3Q19 f64 cells. Start at 8 ranks only when each rank has enough
memory; otherwise record the memory stop condition and start from the smallest
rank count that fits.

## Data Collection

Capture raw logs per command:

```bash
mkdir -p target/cluster-me3/logs
$MPI -n 64 $BIN --mode weak3d --local-edge 128 --steps 200 --warmup 20 --decomp 4x4x4 \
  2>&1 | tee target/cluster-me3/logs/weak3d-r64.log
```

Every run prints `RANK_RESULT` rows for per-rank host/affinity/environment and
one rank-0 `RESULT` row for the aggregate measurement.

## Aggregation + Reporting

Summarize logs into a claims-ledger-ready table:

```bash
RANK_CSV=target/cluster-me3/mpi-ranks.csv \
  ./scripts/qa/aggregate_mpi_results.sh target/cluster-me3/logs/*.log \
  > target/cluster-me3/mpi-summary.csv
```

Use `mpi-summary.csv` for the ME-3 row in `docs/paper/claims-ledger.md`. Keep
`mpi-ranks.csv` with the run artifacts.

## Reproducibility Manifest

Save this manifest next to the logs:

```bash
{
  git rev-parse HEAD
  hostname
  date -Is
  uname -a
  mpirun --version
  mpicc --showme:command 2>/dev/null || true
  mpicc --showme:compile 2>/dev/null || true
  mpicc --showme:link 2>/dev/null || true
  env | sort | grep -E '^(OMPI_|PMI_|PMIX_|FI_|UCX_|SLURM_|OMP_|RAYON_)' || true
  numactl --hardware 2>/dev/null || true
} > target/cluster-me3/manifest.txt
```

Also preserve the scheduler script, node list, network fabric description, rank
mapping policy, and per-rank affinity CSV.

## AWS hpc7g x8 Spend Estimate

`CLUSTER_OPTIONS.md` uses hpc7g.16xlarge x8 for 6 hours as the recommended
AWS campaign unit. Current checked pricing remains consistent with that
estimate: hpc7g.16xlarge in us-east-1 is $1.6832/hour, with 64 vCPUs, 128 GiB
RAM, 200 Gbit/s networking, and no Spot price listed.

Bench wall-time upper bound for one campaign at the recommended grids:

| Item | Upper bound |
|---|---:|
| 1. 3D weak scaling + 2D comparison | 0.5 h |
| 2. 1024^3 strong scaling | 2.0 h |
| 3. mpiP/MPI_T repeat of 64-rank weak run | 0.5 h |
| 4. EFA/TCP or UCX placement checks | 0.25 h |
| 5. rank x thread hybrid sweep | 1.0 h |
| 6. binding/NUMA audit | 0.25 h |
| 7. correctness matrix | 0.5 h |
| 8. diagnostic/gather scaling | 0.5 h |
| Buffer for module rebuilds / scheduler jitter | 0.5 h |
| Total campaign reservation | 6.0 h |

Pricing math:

```text
8 nodes * 6 h * $1.6832/h = $80.7936
$80.7936 * 160 JPY/USD = 12,927 JPY
```

Tokyo region reference:

```text
8 nodes * 6 h * $2.1117/h = $101.3616
$101.3616 * 160 JPY/USD = 16,218 JPY
```

The original ~13k JPY us-east-1 estimate is within 20%, so no correction is
needed. Add a small head-node/EBS buffer separately, and re-confirm pricing
immediately before the owner authorizes spend.
