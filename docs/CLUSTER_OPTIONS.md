# Cluster/Cloud HPC Options — R3 Multi-Node Measurement

Status snapshot (as of 2026-07-07): **ME-3 = RED, blocked on owner spend confirm**
(HANDOFF-PM-2026-07-07 §3). `bench_mpi` 3D + weak/strong modes are ready in code;
`test_mpi.sh` is green locally on Open MPI. What is missing is a multi-node run
against real RDMA fabric to show ≥80% weak scaling at 64 ranks.

Prices retrieved 2026-07-05, 1 USD = 160 JPY. Re-confirm before ordering.

## 1. Resource envelope

Extrapolated from local M5 Max ≈ 40 MLUPS/rank on D2Q9:

- 3D 128³/rank D3Q19 f64 ≈ 0.7–0.8 GiB/rank.
- Weak scaling 1→64 ranks: needs ≥2 nodes (~50 GiB aggregate).
- Strong scaling 1024³ 8→512 ranks: ~370 GiB aggregate (~45 GiB/rank at 8 ranks).
- 220 steps ≈ 25–60 s/config; one full campaign ≈ **50–150 node-hours** compute
  + a few hours setup, i.e. 8–16 nodes × 4–8 h.

## 2. Recommended options

| Option | Campaign cost | Lead time | Coverage of the 8 R3 items | Notes |
|---|---|---|---|---|
| **AWS hpc7g ×8–16 + EFA (ParallelCluster)** | **13,000–35,000 JPY** (us-east-1) | 2–5 days (HPC quota bump) | 8/8 (item 4 = OFI/EFA side only) | **Recommended.** arm64 = same ISA as local. hpc7g.16xlarge $1.6832/h, no Spot. |
| **Fugaku trial** | **0 JPY** (first-touch 1,000 NH) | 2–4 weeks (review + HPCI setup) | 8/8 (item 4 = Tofu-specific) | Free but slow; requires A64FX port + Fujitsu MPI wrapper. Run in parallel with AWS. 8-rank 1024³ start impossible on 32 GiB/node → 16-rank start or 768³. |
| Azure HBv4 ×3–4 Spot | ~5,000 JPY Spot / 28,000 JPY OD | 3–10 days (HB quota review, rejection risk) | 8/8 (best for UCX/IB) | Only worth it as an optional +5-point IB-vs-EFA comparison. |

All other surveyed options — GCP H3/H4D, ABCI 3.0, TSUBAME4, U Tokyo Wisteria,
FOCUS, Hetzner bare metal, single large instance — are dominated by AWS+Fugaku
on cost, coverage, or lead time for this specific campaign. Kept out of scope
here; ABCI/TSUBAME re-enter the picture when M-E moves to GPU cluster work.

## 3. Recommended plan

**Plan B + Fugaku in parallel**:

1. AWS quota bump (512–1,024 vCPU for "Running On-Demand HPC instances",
   us-east-1). 1–3 business days.
2. Bring up `hpc7g.16xlarge ×8–16` cluster via ParallelCluster (EFA on,
   cluster placement group). Measure all 8 R3 items in one half-day session.
   Cost 13k–35k JPY + head-node change.
3. In parallel, submit Fugaku trial (first-touch, 1,000 NH, free). Approval
   yields a second data point on a non-Open-MPI stack (item 7) and a Tofu-D
   fabric (item 4).

Skip: Azure Spot IB comparison unless the quota lands quickly.

## 4. Repo-side readiness (already landed)

- `crates/lbm-core/examples/bench_mpi.rs` — 3D D3Q19, arbitrary Cartesian
  decomposition, RESULT-line format.
- `scripts/bench_mpi.sh` — sweeps ranks × nodes × threads, aggregates CSV.
- `test_mpi.sh` — green on Open MPI locally.

**Not yet drafted** (do before dispatching the AWS run):
- ParallelCluster config YAML (head c7g.large, compute hpc7g×16, EFA,
  placement group, shared /home, auto scale-down).
- Slurm sbatch template that runs the sweep and formats results into the
  MPI_GUIDE table.
- Fugaku porting checklist (Fujitsu MPI mpicc wrapper compatibility, pjsub
  template).

## 5. Owner actions blocking ME-3

1. Confirm spend ceiling (Plan B = 15k–40k JPY).
2. Submit AWS HPC-instance quota request.
3. Submit Fugaku trial application (affiliation + LBM weak/strong-scaling
   summary). Confirm sole-proprietor eligibility with HPCI help desk.

Once (1) is confirmed and (2)/(3) submitted, the repo-side ParallelCluster
YAML and sbatch template can be finalized and the run dispatched the same
week the AWS quota lands.

## 6. Source list (retrieved 2026-07-05)

- AWS: [hpc7g.16xlarge (Vantage)](https://instances.vantage.sh/aws/ec2/hpc7g.16xlarge) ·
  [hpc7g + ParallelCluster setup](https://swsmith.cc/posts/hpc7g-parallelcluster.html) ·
  [ParallelCluster docs](https://docs.aws.amazon.com/parallelcluster/latest/ug/slurm-workload-manager-v3.html) ·
  [EFA docs](https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/efa.html)
- Fugaku: [Trial project](https://www.hpci-office.jp/using_hpci/proposal_submission_current/fugaku_trial) ·
  [Paid usage fees](https://www.hpci-office.jp/using_hpci/proposal_submission_current/fugaku_price)
- Azure: [HBv4 overview](https://learn.microsoft.com/en-us/azure/virtual-machines/hbv4-series-overview) ·
  [HB176rs_v4 (Vantage)](https://instances.vantage.sh/azure/vm/hb176)
- Exchange rate: [Trading Economics USD/JPY](https://tradingeconomics.com/japan/currency) (2026-07-03: 161.27)
