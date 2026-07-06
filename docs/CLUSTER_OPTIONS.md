# Cluster/Cloud HPC Options Survey — Toward R3 Multi-Node Measurement (2026-07-05)

A survey of computing resource options and a decision memo for measuring the 8 items in
docs/MPI_GUIDE.md's "list of measurements to run on a cluster" (R3 completion criteria).
Target audience: individuals to small organizations based in Japan.

**Disclaimer**: Prices and specs are cited per section with source URL and retrieval date (all retrieved
2026-07-05; pre-tax/post-tax follows the source's own notation unless stated otherwise).
Cloud prices fluctuate, so re-confirm before ordering. JPY conversion uses **1 USD = 160 JPY**
(2026-07-03 market rate 161.27 JPY, [Trading Economics](https://tradingeconomics.com/japan/currency))
and **1 EUR = 185 JPY (rough assumption)** for the approximate figures below.
Places marked "to be confirmed" have not been verified against a primary source.

---

## 0. Required Resource Estimate

Requirements back-calculated from the 8 items in MPI_GUIDE (estimates extrapolated from the repo's measured
M5 Max 40 MLUPS/rank, D2Q9 figure):

| Measurement | Scale | Memory required (estimate) | Cores/node required |
|---|---|---|---|
| 1. Weak scaling 1→64 ranks (3D 128³/rank D3Q19 f64) | 64 ranks | f double buffer 0.59 GiB + moments etc. ≈ **0.7–0.8 GiB/rank** (total ~50 GiB) | 64+ physical cores, **must be distributed across 2+ nodes** |
| 2. Strong scaling, fixed 1024³, 8→512 ranks | 512 ranks | Total **~340–370 GiB** (~45 GiB/rank at 8 ranks) | **512 physical cores** + aggregate memory ~400 GiB |
| 3. Communication/compute ratio (mpiP / MPI_T) | Piggybacks on the above | — | Profiler build permissions |
| 4. BTL/MTL (UCX/OFI, eager/rendezvous) | 2+ nodes | — | **The more RDMA-class fabric (IB/EFA etc.), the more valuable** |
| 5. Rank × thread hybrid grid | 1–2 nodes | — | More cores per node lets you sweep a wider grid |
| 6. --map-by/--bind-to, NUMA | 1+ node | — | Linux required (macOS not usable), 2+ NUMA domains preferred |
| 7. All-PASS of test_mpi.sh under another MPI implementation (MPICH family) | Small scale | — | Feasible if MPICH can be built |
| 8. Scaling limits of diagnostic Allreduce / rank-0 gather | 64–512 ranks | gather brings the whole field to rank0 (~170 GiB/field at 1024³ is unrealistic → use a reduced problem to capture the trend) | More ranks is better |

- One step at 3D 128³/rank takes ~0.1–0.3 s/rank (assuming 8–20 MLUPS/rank). Measuring 220 steps ≈ 25–60 s/configuration.
  Weak-scaling curve + hybrid grid + strong-scaling curve + profile capture amounts to **1–2 hours of net compute**;
  including environment setup and test runs, treat **1 session = 8–16 nodes × 4–8 hours (≈ 50–150 node-hours)**
  as the standard unit of a campaign.
- Note: the 8-rank starting point for 1024³ needs 45 GiB/rank → with 128 GiB/node this works if distributed as
  2 ranks/node × 4 nodes. On Fugaku-class hardware (32 GiB/node) an 8-rank starting point is physically impossible
  (start instead from 16 ranks at 1 rank/node ≈ 23 GiB/rank, or shrink to 768³).

---

## 1. Overview Comparison Table

Costs are approximate for "1 campaign (≈50–150 node-hours)" as defined above. ◎○△× indicate how far each of
the 8 items can be measured (details in §5).

| Option | Campaign cost estimate | Lead time | Coverage | Effort | One-liner |
|---|---|---|---|---|---|
| **AWS hpc7g ×8–16 (ParallelCluster + EFA)** | **13,000–35,000 JPY** (us-east-1) | 2–5 days (quota request) | **8/8 ◎** | Medium | The frontrunner. arm64, same ISA as local machine |
| Azure HBv4 ×3–4 (InfiniBand NDR 400G) | OD 28,000 JPY / Spot ~5,000 JPY | 3–10 days (HB quota review) | **8/8 ◎** (best for item 4 via UCX/IB) | Medium–High | Best for a controlled IB comparison experiment |
| GCP H3 ×8 / H4D ×4 | H3 38,000 JPY / H4D 53,000 JPY | 2–5 days | H3 7/8 (no RDMA) / H4D 8/8 | Medium | Only H4D Singapore is near Japan |
| Fugaku trial project (free tier) | **0 JPY** | **2–4 weeks** (1–2 weeks review + HPCI procedures) | 8/8 (item 4 is Tofu-specific) | Medium–High | Free, rolling acceptance. Requires A64FX porting |
| ABCI 3.0 (rt_HF multi-node) | Actual usage ~170,000 JPY / **minimum purchase 220,000 JPY** | 3–4 weeks | 8/8 | Low–Medium | Paying for H200×8 nodes for CPU-only measurement is overkill |
| TSUBAME4.0 (industrial use, published results) | **Minimum 1 unit = 110,000 JPY** (actual usage ~15,000 JPY worth) | Several weeks (to be confirmed) | 8/8 | Low–Medium | 192 cores/node + IB. A stepping stone toward GPU deployment |
| U Tokyo Wisteria/Miyabi (trial) | A few thousand JPY/set~ (corporate quota to be confirmed) | Several weeks (to be confirmed) | 8/8 (Odyssey is same family as Fugaku) | Medium–High | Cheapest in the academic tier. Corporates go via the trial quota |
| FOCUS supercomputer (industry-only) | **10,000–30,000 JPY** (annual fee 10,000 + usage, first-year free tier worth 10,000 JPY) | 1–3 weeks (inquiry needed) | 6.5/8 (512 ranks not possible) | Low–Medium | Domestic, cheap, has InfiniBand. Small scale |
| Hetzner dedicated servers ×4 + manual MPI | **~150,000–180,000 JPY/month + ~60,000 JPY setup** (monthly contract) | 3–10 days | 4/8 (the ≥80% pass line is hard due to network bottleneck) | **High** | Makes sense only for a permanent mini-cluster use case |
| Single large instance (c8g/c7a.48xlarge, 192 vCPU) | **3,000–13,000 JPY** | **Same day** | 4.5/8 (no inter-node) | **Low** | A skirmish achievable today. Cannot achieve R3 alone |

---

## 2. Cloud HPC

### 2.1 AWS — hpc7g + EFA + ParallelCluster (the frontrunner)

- **hpc7g.16xlarge**: Graviton3E **64 physical vCPU / 128 GiB / EFA 200 Gbps**.
  On-demand **$1.6832/h (us-east-1), $2.1117/h (Tokyo ap-northeast-1, +25%)**. **No Spot support** (true of all HPC-class instances).
  Source: [Vantage](https://instances.vantage.sh/aws/ec2/hpc7g.16xlarge), [aws-pricing.com](https://aws-pricing.com/hpc7g.16xlarge.html) (2026-07-05).
  Available regions are US East (N. Virginia) / Tokyo / Ireland / GovCloud ([AWS news, 2023-09](https://aws.amazon.com/about-aws/whats-new/2023/09/amazon-ec2-hpc7g-instances-additional-regions/)).
  As of 2026-07-05 the Graviton4-generation hpc8g has not been confirmed to exist.
- **hpc7a.96xlarge** (AMD EPYC 9R14 **192 cores / 768 GiB / EFA 300 Gbps**): $9.0793/h, no Spot support
  ([Vantage](https://instances.vantage.sh/aws/ec2/hpc7a.96xlarge), 2026-07-05). An alternative if you want fewer nodes.
- **Example configurations and cost** (us-east-1, 1 USD = 160 JPY):
  - hpc7g ×8 × 6 h = 48 NH × $1.6832 ≈ **$81 ≈ 13,000 JPY** (512 physical cores: sufficient for both 64-rank weak scaling and 512-rank strong scaling)
  - hpc7g ×16 × 8 h = 128 NH ≈ **$215 ≈ 35,000 JPY** (a generous upper bound)
  - Head node (c7g.large etc.) + shared EBS costs on the order of a few hundred JPY. Tokyo region adds +25%.
- **MPI support**: EFA goes through libfabric (OFI). ParallelCluster automatically provisions the EFA driver + Open MPI + Slurm
  ([ParallelCluster docs](https://docs.aws.amazon.com/parallelcluster/latest/ug/slurm-workload-manager-v3.html);
  hpc7g setup example: [Sean Smith's guide](https://swsmith.cc/posts/hpc7g-parallelcluster.html)).
  rsmpi is expected to pass its mpicc probe cleanly (Open MPI family). Item 4 lets you measure the **OFI/EFA side**,
  including a comparison against `FI_PROVIDER=tcp`. The UCX/InfiniBand side cannot be measured this way (→ supplement with Azure or a domestic IB machine).
- **Setup effort (medium)**: (1) account + **request a service quota increase for "Running On-Demand HPC instances"**
  (new accounts sometimes have a small/zero default; 1–3 business days from request, requires describing the use case)
  (2) `pip install aws-parallelcluster` → generate a Slurm cluster from a single YAML file (3) specify a placement group (cluster).
  There is also managed Slurm via AWS PCS, but ParallelCluster is sufficient for a one-off campaign.
- **Note**: HPC instances are limited to specific AZs within a region. After running, be sure to tear the cluster down with
  `pcluster delete-cluster` (to stop billing).

### 2.2 Azure — HBv4 (InfiniBand NDR 400 Gbps)

- **HB176rs_v4**: AMD EPYC Genoa-X **176 physical cores (SMT disabled) / 768 GiB / NDR InfiniBand 400 Gb/s**
  ([Microsoft Learn](https://learn.microsoft.com/en-us/azure/virtual-machines/hbv4-series-overview)).
  On-demand **$7.20/h, Spot $1.331/h** (reference region, [Vantage](https://instances.vantage.sh/azure/vm/hb176), 2026-07-05).
  Its biggest draw is that **Spot is available despite being an HPC VM** (the preemption risk is easy to tolerate for a short benchmark run).
- Regions: East US family, West Europe, **Southeast Asia, Korea Central**, etc. **Not available in East Japan**
  ([Spare Cores](https://sparecores.com/server/azure/Standard_HB176rs_v4), 2026-07-05).
- Example configuration: ×4 nodes (704 cores) × 6 h = OD $173 ≈ **28,000 JPY** / Spot $32 ≈ **5,000 JPY**.
  1024³ strong scaling fits within 3 nodes (528 cores, 2.3 TiB).
- MPI: the Azure HPC image (AlmaLinux/Ubuntu-HPC) bundles Mellanox OFED + HPC-X (Open MPI/UCX) + Intel MPI.
  **This is the most HPC-authentic way to measure item 4's UCX/IB side (eager/rendezvous, tag matching)**.
- Effort (medium–high): the HB family requires a vCPU quota request, and **new/pay-as-you-go subscriptions are sometimes rejected**
  (via a support request, taking a few days or more). Cluster assembly is via CycleCloud or manual setup (same PPG + IB verification).
  This is a notch more tedious than AWS.
- The next-generation HBv5 (HBM-equipped EPYC, 800 Gbps IB) has been announced but pricing was not confirmed in this survey (to be confirmed).

### 2.3 GCP — H3 / H4D + compact placement

- **h3-standard-88**: Sapphire Rapids **88 vCPU (SMT disabled) / 352 GB / 200 Gbps**, supports compact placement.
  **$4.9236/h (us-central1), no Spot support**. Available only in the 3 regions us-central1 / europe-west4 / northamerica-northeast1
  ([gcloud-compute.com](https://gcloud-compute.com/h3-standard-88.html), 2026-07-05).
  **No RDMA between nodes (gVNIC/TCP)** → the drawback is that this reduces the measurement value for item 4.
  ×8 × 6 h = $236 ≈ **38,000 JPY**.
- **h4d-standard-192**: EPYC Turin **192 vCPU (SMT disabled) / 720 GB / Cloud RDMA (Falcon)**.
  **GA as of 2026-03**, us-central1-a / europe-west4-b / **asia-southeast1-a (Singapore)**
  ([Google Cloud blog](https://cloud.google.com/blog/products/compute/h4d-vms-now-ga)).
  Reference price **~$13.74/h equivalent** ($10,033/month, [CloudPrice](https://cloudprice.net/gcp/compute/instances/h4d-standard-192), 2026-07-05.
  The blog post mentions DWS Flex Start pricing "from 3 cents/core-hour" → could be less than half depending on conditions).
  ×4 × 6 h ≈ $330 ≈ **53,000 JPY**.
- Effort (medium): Slurm can be built with Cluster Toolkit. Quota and zone stock availability need checking.
  Configuring RDMA-capable MPI (supported Intel MPI / Open MPI versions) takes some extra work per GCP documentation.
- Verdict: neither price nor available information gives a clear edge over AWS/Azure. Recommended **only if you already have
  existing assets on GCP**.

### 2.4 Reference — OCI (not surveyed)

Oracle Cloud's BM.HPC line reportedly offers bare metal + RDMA cluster networking at low cost, but
price and availability were not surveyed this time (to be investigated; worth keeping as a candidate).

---

## 3. Domestic Academic/Public Options

### 3.1 Fugaku — the trial project (free) is the top candidate

- **Trial project (general/industrial): free, rolling acceptance, review takes 1–2 weeks, up to 6 months, cap of 100,000 NH**.
  A simpler "first-touch option" is also available: a fixed 1,000 NH, up to 3 months
  ([HPCI Fugaku trial project](https://www.hpci-office.jp/using_hpci/proposal_submission_current/fugaku_trial), 2026-07-05).
  Even 1,000 NH is more than enough for this campaign (50–150 NH).
- If moving to paid usage: pay-as-you-go **98.64 JPY/NH (unpublished results), 49.32 JPY/NH (published results)**
  ([HPCI pricing page](https://www.hpci-office.jp/using_hpci/proposal_submission_current/fugaku_price), 2026-07-05).
  Example: 128 nodes × 10 h = 1,280 NH ≈ 126,000 JPY (unpublished) — even paid usage is on par with cloud pricing.
- Node: A64FX 48 cores / **32 GiB HBM2** / Tofu-D (nominal specs: [R-CCS](https://www.r-ccs.riken.jp/fugaku/)).
  - Weak scaling at 64 ranks: achievable with, e.g., 8 nodes × 8 ranks. 512 ranks also fits comfortably from 11 nodes up.
  - **The 8-rank starting point for 1024³ is not possible due to the 32 GiB/node constraint** (see §0 note; switch to a
    16-rank starting point or shrink to 768³).
- Effort (medium–high): (1) adds the overhead of HPCI account and digital certificate procedures (2) **aarch64 is fine, but
  MPI is Fujitsu MPI (Open MPI–based, wrapper name mpifcc etc.)**. There may be some friction around rsmpi's mpicc probe
  compatibility and bindgen (needs verification). (3) Jobs go through pjsub (Fujitsu TCS). This has high value for extending
  item 7's "implementation-independent bit match" beyond the Open MPI family.
- Eligibility: the call for proposals does not exclude corporations/universities, but **whether an unaffiliated individual
  can apply needs to be confirmed** (the application requires stating affiliation and research plan).

### 3.2 ABCI 3.0 (AIST)

- Pricing: **1 point = 220 JPY (tax included), minimum purchase 1,000 pt = 220,000 JPY**.
  **rt_HF (full node)**, which is required for multi-node MPI, costs **16 pt/h = 3,520 JPY/node-hour** (Spot/on-demand tiers)
  ([ABCI pricing, FY2026](https://abci.ai/ja/how_to_use/tariffs.html), 2026-07-05).
  Example: 8 nodes × 6 h = 768 pt ≈ **169,000 JPY worth** (within the 220,000 JPY minimum purchase).
- Node (H): Xeon Platinum 8558 48c ×2 = **96 physical cores / 2 TB / InfiniBand NDR200 ×8** + **H200 ×8**.
  Multi-node is rt_HF only, **up to 128 nodes** ([ABCI 3.0 docs: system overview](https://docs.abci.ai/v3/ja/system-overview/),
  [job execution](https://docs.abci.ai/v3/ja/job-execution/), 2026-07-05).
- Lead time: usage application review ~10 business days + up to 10 days for point allocation → **effectively 3–4 weeks**
  ([how to use](https://abci.ai/ja/how_to_use/), [point application guide](https://abci.ai/news/2025/01/23/ja_news_Point_Application.html), 2026-07-05).
- Verdict: for CPU-only R3 measurement, **paying 3,520 JPY/h to leave an H200×8 node idle is overpriced**, and the
  220,000 JPY minimum purchase is also a heavy commitment. However, it's the strongest domestic option if you're
  looking ahead to M-E (GPU/CUDA deployment, GPUDirect), so **"handle R3 elsewhere and apply for ABCI when M-E comes up"**
  is the sensible approach. Eligibility: the terms are operated on the premise of corporate affiliation
  (practically difficult for unaffiliated individuals; sole proprietors should inquire directly).

### 3.3 TSUBAME4.0 (Institute of Science Tokyo)

- Pricing (off-campus, pay-as-you-go): **275 JPY/node-hour with published results, 1,100 JPY/node-hour unpublished**.
  **For off-campus users, the minimum purchase unit is 1 unit = 400 node-hours = 110,000 JPY (published) / 440,000 JPY (unpublished)**.
  Points expire at fiscal year end (3/31)
  ([TSUBAME4 pricing overview](https://www.t4.cii.isct.ac.jp/fare_overview), 2026-07-05).
- Node: EPYC 9654 ×2 = **192 physical cores / 768 GiB / InfiniBand NDR200 ×4 / H100 ×4** (for the nominal configuration,
  see the [TSUBAME4 site](https://www.t4.cii.isct.ac.jp/)). Actual campaign usage would be ~40 NH ≈ 11,000 JPY worth,
  but the minimum purchase is 1 unit = 110,000 JPY.
- Verdict: a node including H100s at 275 JPY/h (published-results condition) is an exceptional deal, and **cheaper than
  ABCI if it also serves as a stepping stone toward GPU deployment**. However, factor in the obligation to publish
  results (usage report), the fiscal-year expiration, and the application lead time (several weeks, to be confirmed).
  Whether a trial-use (free tier) option currently exists is not stated on the pricing page and needs confirmation.

### 3.4 University of Tokyo Information Technology Center — Wisteria/BDEC-01, Miyabi

- Regular use (trial), FY2026: **Wisteria/BDEC-01 1 set = 2,250 JPY (720 tokens), Miyabi 1 set = 7,500 JPY (720 tokens)**,
  up to 12 sets ([usage fee page](https://www.cc.u-tokyo.ac.jp/guide/application/charge_trial.php), 2026-07-05.
  The token-to-node-hour conversion factor is machine-specific, so check the same page for details).
  **A separate "corporate use trial" quota (Wisteria)** exists for companies
  ([guide](https://www.cc.u-tokyo.ac.jp/guide/trial/company.php)).
- Wisteria-Odyssey is A64FX + Tofu-D (same family as Fugaku) → serves as a preliminary step or alternative to the Fugaku trial.
  Miyabi is primarily GH200-based (on the GPU side of the roadmap).
- Verdict: pricing is among the cheapest, but **academic use is the primary target audience, and companies/individuals must
  confirm quota and review requirements**. Lead time is several weeks (to be confirmed).

### 3.5 FOCUS Supercomputer (Foundation for Computational Science, Kobe) — industry-only

- Pricing (pay-as-you-go): **F system (CPU) 300 JPY/node-hour** (tapering down to 150 JPY depending on node count used),
  **X system (A64FX) 80 JPY/node-hour**, S system (EPYC 9654 192c, billed per 8-core VM unit) 60 JPY/VM-hour, etc.
  Account fee **10,000 JPY/user/fiscal year**, **first year includes a free tier worth 10,000 JPY of pay-as-you-go usage**.
  A free **trial-use** program also exists (for porting/benchmarking purposes)
  ([pricing](https://www.j-focus.or.jp/focus/fee.html), [usage types](https://www.j-focus.or.jp/focus/form.html),
  [trial use](https://www.j-focus.or.jp/focus/free-trial.html), 2026-07-05).
- Systems ([user guide](https://www.j-focus.jp/user_guide/ug0001000000/), 2026-07-05):
  F = Xeon E5-2698v4 40c/128GB/**InfiniBand FDR 56G** (12-node scale, deployed 2016),
  X = **A64FX 48c/32GB/InfiniBand EDR 100G**, Z = Xeon 40c/**EDR 100G**, S = EPYC 9654 192c/**100GbE**.
  Node counts for each system require direct inquiry (unyo@j-focus.or.jp).
- Verdict: **domestic, cheap, real InfiniBand, industry-only usage (small/medium businesses welcome)**, and items
  1/3/4/5/6/7 can all be measured. However, the scale is small, and **512-rank strong scaling (item 2) is not possible**
  (F tops out at 480 cores). The X system (A64FX+EDR) at 80 JPY/h also doubles as a rehearsal for Fugaku porting.
  Sole-proprietor eligibility and lead time need direct inquiry (roughly 1–3 weeks).

---

## 4. The Quick-and-Cheap Tier

### 4.1 Earning rank count from a single large-core-count instance

- Candidates (192-vCPU class, retrieved 2026-07-05):
  - **c8g.48xlarge** (Graviton4, 192 vCPU/384 GiB): OD $7.657/h, **Spot $2.727/h** ([Vantage](https://instances.vantage.sh/aws/ec2/c8g.48xlarge)).
    Being arm64, its advantage is matching the ISA of the local M5 Max and hpc7g.
  - **c7a.48xlarge** (EPYC, 192 vCPU/384 GiB): OD $9.853/h, **Spot $3.312/h** ([Vantage](https://instances.vantage.sh/aws/ec2/c7a.48xlarge)).
  - GCP's c3d/c4d 360-vCPU class, and a lone Azure HBv4 (176c, Spot $1.33/h), can serve the same purpose.
- Cost: **on Spot, 3,000–4,000 JPY for 8 hours; even on-demand, around 10,000 JPY. Can start the same day** (just watch
  the vCPU quota for the large instance size).
- What this achieves: item 1's "single-node version of 64 ranks" (re-running MPI_GUIDE's existing table on homogeneous
  cores under Linux, extended up to n≤64), plus **item 5 (hybrid grid), item 6 (map-by/bind-to, 2 NUMA domains), item 7
  (re-verifying the shared-memory path with MPICH installed), and item 8 (trend up to ~192 ranks)**.
- **Limitations (stated clearly)**: since all ranks are connected via shared memory (vader/xpmem), **the inter-node
  network (item 4) and R3's central goal of "multi-node weak scaling ≥80%" cannot be measured in principle**.
  Also, 384 GiB cannot hold 1024³ (~370 GiB), so item 2 requires either shrinking to 768³ or using r8g/m8g etc.
  (768 GiB–1.5 TiB, pricing to be confirmed). Its role is "calibration/debugging before the real cluster run and
  establishing a single-node baseline."

### 4.2 Multiple Hetzner-style bare-metal servers + manual MPI

- Price example: **AX162-R** (EPYC 9454P 48c/96t, 256 GB DDR5), **€199/month + €79 setup fee**
  ([Hetzner press release](https://www.hetzner.com/pressroom/new-ax162/); there's also a report of a 2026-06 price
  increase to €238/month: [Northflank roundup](https://northflank.com/blog/hetzner-cloud-server-price-increases),
  both retrieved 2026-07-05).
  4 units: **€796–952/month ≈ 150,000–180,000 JPY + €316 ≈ 60,000 JPY setup fee** (monthly contract, not billed hourly).
  Hetzner Cloud's CCX63 (48 vCPU = 24 physical cores, €1.37/h, [Spare Cores](https://sparecores.com/server/hcloud/ccx63))
  is billed hourly, but has half the physical cores and no network guarantee, so it is not recommended for MPI use.
- Network: standard 1 GbE; vSwitch/10G uplink is optional (configuration and extra cost to be confirmed). **No RDMA
  (TCP only)**.
- Quantitative estimate: 3D 128³/rank face-halo ≈ 0.6–2.5 MB/face/step. At 16 ranks/node, tens of MB/step cross node
  boundaries, with a step time of 0.1–0.3 s → **1 GbE (~120 MB/s) is certain to saturate, and even 10 GbE would struggle
  to hit ≥80% efficiency once TCP latency is included**. This tends to fail R3's pass line "because of the network,
  not because of the code."
- Verdict: you can capture item 4's "behavior under the TCP BTL" and item 1's lower bound, but **this is not the main
  battleground for achieving R3**. Worth reconsidering once you want a permanently standing cluster (e.g., to run
  CI-style) at a fixed monthly cost.

---

## 5. Coverage Table for the 8 Items

◎ = directly measurable ○ = measurable (with caveats) △ = partial × = not possible

| # | Measurement (MPI_GUIDE §measurement list) | AWS hpc7g | Azure HBv4 | GCP H3/H4D | Fugaku (trial) | ABCI 3.0 | TSUBAME4 | FOCUS | Hetzner | Single large node |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | Weak scaling 1→64 (multi-node) | ◎ | ◎ | ○/◎ | ◎ | ◎ | ◎ | ○(F/X/Z) | △(low efficiency) | △(intra-node only) |
| 2 | Strong scaling 1024³, 8→512 | ◎(8 nodes) | ◎(3 nodes) | ○(H3 is TCP-bound) | ○(16-rank start) | ◎ | ◎ | ×(≤480 cores) | × | ×(up to 192, memory-limited) |
| 3 | Communication/compute ratio (mpiP/MPI_T) | ◎ | ◎ | ◎ | ○(Fujitsu MPI convention) | ◎ | ◎ | ◎ | ◎ | ○(shared-memory ratio) |
| 4 | BTL/MTL and thresholds (UCX/OFI) | ○(OFI/EFA side) | ◎(UCX/IB side) | △(H3 TCP)/○(H4D RDMA) | ○(Tofu-specific) | ◎(IB NDR) | ◎(IB NDR) | ◎(IB FDR/EDR) | △(TCP only) | ×
| 5 | Rank × thread optimum | ◎(64c) | ◎(176c) | ◎ | ○(48c/32GB) | ◎(96c) | ◎(192c) | ○ | ○ | ◎(192c) |
| 6 | map-by/bind-to, NUMA | ◎ | ◎ | ◎ | ○(pjsub convention) | ◎ | ◎ | ◎ | ◎ | ◎(2 NUMA) |
| 7 | test_mpi.sh under other MPI implementations | ◎(MPICH addable) | ◎(HPC-X/Intel MPI bundled) | ◎ | ◎(Fujitsu MPI = separate lineage) | ◎ | ◎ | ○ | ◎ | ○(shared-memory path only) |
| 8 | Diagnostic Allreduce/gather scaling limits | ◎(512 ranks) | ◎ | ◎ | ◎(and beyond) | ◎ | ◎ | △(up to ~480) | △ | △(up to ~192) |

---

## 6. Recommended Scenarios

### Plan A — Cheapest, gets 60 points (total: a few thousand to 15,000 JPY, starting this week)

1. **Today or tomorrow**: rent **1 AWS c8g.48xlarge (Spot)** for a few hours and work through
   items 5, 6, 7 (shared-memory range), 8 (up to ~192 ranks), and the "single-node 64-rank baseline."
   Being Linux/Graviton, the local arm64 build should carry over almost unchanged. Cost 3,000–10,000 JPY.
   (If touching 1024³, either shrink to 768³ or pick a larger-memory r8g-family instance.)
2. **In parallel**: **apply for the Fugaku trial project (first-touch, 1,000 NH, free)** (review 1–2 weeks + HPCI procedures).
   If approved, the remaining core goal (multi-node weak scaling, items 2 and 4) can be measured **for free**.
- Outcome: until Fugaku is approved, only R3's central "≥80% at 64 ranks" multi-node measurement remains outstanding
  (= 60 points). After Fugaku access arrives you can push to 100 points, but this carries **the technical risk of
  porting to A64FX and confirming connectivity with Fujitsu MPI, plus a wait of several weeks**.

### Plan B — Solid, gets 90 points (total: 20,000–50,000 JPY, lead time 2–5 days)

1. On an AWS account, **request a quota increase for HPC instances** (512–1,024 vCPU, 1–3 business days).
2. **Stand up a Slurm cluster of hpc7g.16xlarge ×8–16 (us-east-1, EFA, cluster placement group) via ParallelCluster**
   for just half a day. Measure all 8 items in a single session (item 4 covers the OFI/EFA side).
   Cost 13,000–35,000 JPY + a small amount for the head node.
3. (Optional, +5 points) Add **Azure HBv4 Spot ×3–4 for a few hours** (~5,000 JPY) to run a controlled comparison
   of the UCX/InfiniBand side (item 4) and "how efficiency differs between EFA and IB" → 95 points.
   Note: since the HB quota review takes days and carries rejection risk, treat this as "do it if you can get it."
- Outcome: R3's completion criteria can be satisfied on your own schedule. The measurements are reproducible
  (the YAML + scripts can be kept in the repo).

### Recommended in parallel (regardless of plan)

- **Since the Fugaku trial project is free, there's no downside to submitting it even alongside Plan B**
  (it's the strongest sample for item 7's "non–Open MPI" case).
- **Position ABCI 3.0 / TSUBAME4 not for R3 but as a stepping stone toward M-E (GPU/CUDA, GPUDirect)**.
  Paying a 110,000–220,000 JPY minimum purchase for a CPU-only R3 measurement is not worth it.
- If you need everything to stay domestic or must pay by invoice, **FOCUS (trial use → pay-as-you-go)** is the
  realistic answer (though you'd have to give up on item 2's 512-rank measurement).

---

## 7. Next Actions

### Things for the user (Mr. Tsuzuki) to do

1. **Decide on a direction**: Plan A / Plan B / running both in parallel (recommendation: "Plan B + parallel Fugaku
   application"), and confirm the budget ceiling.
2. Prepare an AWS account and **submit the quota request** (Plan B: 512–1,024 vCPU for "Running On-Demand HPC
   instances" in us-east-1; Plan A: 192 vCPU of Standard/Spot). The request text can be drafted on our side.
3. **Apply for the Fugaku trial project** (including obtaining an HPCI account): requires stating affiliation and a
   project summary (a few hundred characters describing LBM solver weak/strong scaling measurement).
   Confirm eligibility as a sole proprietor with the HPCI help desk beforehand.
4. (Optional) Create an Azure subscription and request an HB-family quota (if rejected, give up and use AWS only).
5. If using FOCUS: inquire at unyo@j-focus.or.jp about usage eligibility (small business), the node count for the
   X/Z systems, and availability of trial use.

### Things we (the repo side) can prepare in advance

1. **3D support for bench_mpi**: the current `crates/lbm-core2/examples/bench_mpi.rs` is fixed to 2D 512²/rank with
   band decomposition. Extend it to 128³/rank, D3Q19, arbitrary Cartesian decomposition, and a common RESULT-line
   format (ranks/nodes/threads/MLUPS/efficiency).
2. **Slurm-based driver**: templatize `scripts/bench_mpi.sh` into an sbatch script
   (running the ranks × nodes × threads grid, sweeping `--map-by`/`--bind-to` variants, aggregating results into
   CSV, and formatting into MPI_GUIDE's table format).
3. **Draft ParallelCluster config YAML** (head: c7g.large, compute: hpc7g×16, EFA on, placement group, shared /home,
   automatic scale-down after job completion) and an execution runbook.
4. **Profiling procedure for item 3**: notes on building mpiP, or an alternative measurement script using Open MPI's
   `MPI_T`/OSU micro-benchmarks.
5. **MPI implementation matrix for item 7**: procedure for running test_mpi.sh under Open MPI / MPICH / Intel MPI
   (including rebuilding rsmpi).
6. **Fugaku porting checklist**: aarch64 build (done: already routine on M5 Max), Fujitsu MPI's mpicc wrapper
   compatibility, a pjsub job template.

### ME-3 preparation update (2026-07-07)

`docs/CLUSTER_RUNBOOK.md` now provides the no-spend handoff runbook for the
recommended hpc7g x8 campaign: exact `bench_mpi` commands for all 8 MPI_GUIDE
measurements, local preflight, aggregation, and reproducibility manifest. The
spend math was rechecked against the same hpc7g.16xlarge us-east-1 price
($1.6832/h): `8 nodes * 6 h * $1.6832/h * 160 JPY/USD = 12,927 JPY`. This is
consistent with the original ~13,000 JPY estimate, so no correction is needed.

---

## 8. Source List (all retrieved 2026-07-05)

- AWS: [hpc7g product page](https://aws.amazon.com/ec2/instance-types/hpc7g/) /
  [Vantage hpc7g.16xlarge](https://instances.vantage.sh/aws/ec2/hpc7g.16xlarge) /
  [aws-pricing.com hpc7g.16xlarge (by region)](https://aws-pricing.com/hpc7g.16xlarge.html) /
  [hpc7g additional regions](https://aws.amazon.com/about-aws/whats-new/2023/09/amazon-ec2-hpc7g-instances-additional-regions/) /
  [Vantage hpc7a.96xlarge](https://instances.vantage.sh/aws/ec2/hpc7a.96xlarge) /
  [Vantage c8g.48xlarge](https://instances.vantage.sh/aws/ec2/c8g.48xlarge) /
  [Vantage c7a.48xlarge](https://instances.vantage.sh/aws/ec2/c7a.48xlarge) /
  [ParallelCluster docs](https://docs.aws.amazon.com/parallelcluster/latest/ug/slurm-workload-manager-v3.html) /
  [hpc7g + ParallelCluster setup example](https://swsmith.cc/posts/hpc7g-parallelcluster.html) /
  [EFA docs](https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/efa.html)
- GCP: [gcloud-compute.com h3-standard-88](https://gcloud-compute.com/h3-standard-88.html) /
  [H3 announcement blog](https://cloud.google.com/blog/products/compute/new-h3-vm-instances-are-optimized-for-hpc) /
  [H4D GA blog](https://cloud.google.com/blog/products/compute/h4d-vms-now-ga) /
  [CloudPrice h4d-standard-192](https://cloudprice.net/gcp/compute/instances/h4d-standard-192)
- Azure: [HBv4 series overview (Microsoft Learn)](https://learn.microsoft.com/en-us/azure/virtual-machines/hbv4-series-overview) /
  [Vantage HB176rs_v4](https://instances.vantage.sh/azure/vm/hb176) /
  [Spare Cores HB176rs_v4 (regions)](https://sparecores.com/server/azure/Standard_HB176rs_v4)
- ABCI: [Pricing (FY2026)](https://abci.ai/ja/how_to_use/tariffs.html) / [How to use](https://abci.ai/ja/how_to_use/) /
  [Point application](https://abci.ai/news/2025/01/23/ja_news_Point_Application.html) /
  [ABCI 3.0 system overview](https://docs.abci.ai/v3/ja/system-overview/) / [Job execution](https://docs.abci.ai/v3/ja/job-execution/)
- Fugaku: [Trial project (rolling acceptance)](https://www.hpci-office.jp/using_hpci/proposal_submission_current/fugaku_trial) /
  [Paid project usage fees](https://www.hpci-office.jp/using_hpci/proposal_submission_current/fugaku_price) /
  [Fugaku (R-CCS)](https://www.r-ccs.riken.jp/fugaku/)
- TSUBAME4.0: [Pricing overview](https://www.t4.cii.isct.ac.jp/fare_overview)
- U Tokyo: [Regular use (trial) fees](https://www.cc.u-tokyo.ac.jp/guide/application/charge_trial.php) /
  [Corporate use (trial)](https://www.cc.u-tokyo.ac.jp/guide/trial/company.php)
- FOCUS: [Pricing](https://www.j-focus.or.jp/focus/fee.html) / [Usage types](https://www.j-focus.or.jp/focus/form.html) /
  [Trial use (free)](https://www.j-focus.or.jp/focus/free-trial.html) /
  [System overview (user guide)](https://www.j-focus.jp/user_guide/ug0001000000/)
- Hetzner: [AX162 press release](https://www.hetzner.com/pressroom/new-ax162/) /
  [AX162-R product page](https://www.hetzner.com/dedicated-rootserver/ax162-r/) /
  [2026-06 price increase roundup (Northflank)](https://northflank.com/blog/hetzner-cloud-server-price-increases) /
  [Spare Cores CCX63](https://sparecores.com/server/hcloud/ccx63)
- Exchange rate: [Trading Economics USD/JPY](https://tradingeconomics.com/japan/currency) (2026-07-03: 161.27)
