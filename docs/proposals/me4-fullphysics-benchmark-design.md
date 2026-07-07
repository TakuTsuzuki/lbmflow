# ME-4 measurement design — full-physics stirred-tank degradation benchmark

Produced 2026-07-07 by a PM-commissioned design survey (Plan agent) of
crates/lbm-core. Coordination-ready measurement spec (no code); the PM can
dispatch Stage A now or hand it to the QA-sweep session that owns M-F.
Line refs current as of main a375ba2.

## 1. Landed vs pending — what the benchmark can toggle today

Confirmed by source:

| M-F subsystem | State | Evidence | Runnable now? |
|---|---|---|---|
| Central-moment/cumulant collision (MF-α) | LANDED all 3 backends | params.rs:32 CollisionKind::CentralMoment; SIMD backend_simd.rs:510-551; GPU wgsl.rs:407, backend.rs:1018 | Yes: scalar/SIMD/GPU |
| LES (WALE) (MF-β) | LANDED, externally-driven omega | les.rs:47 WaleLes, les.rs:120 update installs next-step omega; solver.rs:3095 set_omega_field; SIMD reads backend_simd.rs:1535,1751 | Yes CPU; GPU on-device omega is a separate narrower cell (driver is CPU-side gather/scatter) |
| Rotating IBM rotor (MF-δ) | LANDED apply_rotating_ibm | solver.rs:2153, rotating_ibm.rs:24/146 | Yes CPU (Solver method on gathered fields; not a device kernel) |
| Two-phase interface | Shan-Chen only; Allen-Cahn W-VOF is a prescribed-velocity O1 transport stub, D3Q19+local-CPU only, no momentum coupling | solver.rs:2492, 2616-2720; REQ:594 CRITICAL PATH | No two-way multiphase |
| Scalar ADE | NOT PRESENT | no Scalar/ADE/D3Q7 in src | No (blocks on W-SCAL/MF-ε) |
| Particles | one-way only, standalone ParticleSet | particles.rs:73/104/7; solver.rs:3258 reserves particles:false | Yes bolt-on (zero flow cost; two/four-way pending) |
| D3Q27 | collision equivalence validated; open-face CPU landed | backend_simd.rs:2166, gpu/backend.rs:2507 | Yes periodic/closed; GPU 3D open-face restricted |

Runnable set today: {central-moment collision, WALE-LES, rotating IBM,
one-way particles, D3Q27}. NOT runnable: {resolved Allen-Cahn two-phase,
scalar ADE, active feedback, two/four-way particles} — block on MF-γ (W-VOF)
and MF-ε (W-SCAL/REACT). Matches claims-ledger.md:26.

CAVEAT: the shipped stirred_tank_3d example (crates/lbm-cli/examples/
stirred_tank_3d.rs:11-15) drives the impeller with volume PENALIZATION, not
the landed apply_rotating_ibm, and enables neither WALE nor phase field.
Reuse its GEOMETRY (Rushton D=T/3, 4 baffles), not its physics path; the ME-4
harness swaps in apply_rotating_ibm to measure the landed rotor's cost.

## 2. Metric

ratio(subsystem_set) = MLUPS_full(G,B,P,set) / MLUPS_baseline(G,B,P), reported
as a slowdown "1/ratio × slower". baseline = single-phase D3Q19 TRT (Λ=3/16),
pure fluid step. MLUPS = cells·measured_steps/wall_s/1e6, warmup excluded
(bench_backends.rs:77-81). NEVER report a bare ratio — always tag the enabled
set, e.g. ratio[cm+wale+ibm+part].

Incremental ladder (marginal cost per subsystem):

| Row | Collision | LES | Rotor | Particles |
|---|---|---|---|---|
| B0 | TRT | off | off | off (baseline) |
| C1 | CentralMoment | off | off | off |
| C2 | CentralMoment | WALE | off | off |
| C3 | CentralMoment | WALE | IBM | off |
| C4 | CentralMoment | WALE | IBM | one-way (full landed-today) |

Cross with grid G∈{96³,128³,192³} (128/192 are the anchored GPU sizes),
backend B∈{CpuSimd, GPU} (CpuScalar = correctness anchor only), precision
P∈{f32,f64} on CPU, f32-only on GPU.

Backend honesty gates (MUST be in the published table):
- GPU: only collision toggles are device kernels; WALE-update and
  apply_rotating_ibm are CPU-side host ops in the landed code. Restrict the
  GPU headline to B0→C1 (fair device-vs-device); run the WALE/IBM ladder on
  CpuSimd where it is a real in-process cost. Do NOT publish a GPU C4 ratio
  as kernel-fused.
- D3Q27: run a SEPARATE D3Q27 B0→C1 collision ladder (periodic/closed);
  keep apples-to-apples (D3Q27-baseline vs D3Q27-full, never mixed with
  D3Q19-baseline).

## 3. M-Star comparability + caveats

Correction to common framing: M-Star CFD is a GPU-resident LATTICE-BOLTZMANN
transient solver (not FVM; the FVM incumbents are Fluent/STAR-CCM+). It quotes
throughput as (million) cell-updates/s and wall-clock-per-simulated-second for
a named stirred-tank case on a named datacenter GPU (A100/H100), for a full
coupled config (LES/DES + free surface or Eulerian multiphase + Lagrangian
particles + species + rotating impeller).

Comparability: our absolute MLUPS is directly comparable to M-Star cell-
updates/s (LBM↔LBM, legitimate); our degradation ratio is the fairer cross-
vendor claim (cancels hardware + kernel-tuning maturity). Mandatory caveats to
print: (1) coverage asymmetry — our "full-today" is a SUBSET (no resolved
multiphase/scalars/two-way), so the ratio is optimistic and must be labeled
"landed-physics subset" until MF-γ/ε land (claims-ledger stays RED); (2)
LBM-vs-LBM cell-update parity only — never vs FVM iteration counts; (3)
hardware disclosure (M5 Max Metal vs A100/H100 — the ratio transfers, absolute
MLUPS does not); (4) rotor fidelity — ours is direct-forcing IBM, M-Star uses
sliding-mesh/immersed; MRF/overset pending.

## 4. Staged plan

Stage A — "ME-4a landed-physics degradation ratio (subset)", RUNNABLE NOW:
the B0→C4 ladder on CpuSimd (f32+f64) + B0→C1 on GPU (f32). New example
crates/lbm-core/examples/bench_stirred.rs (peer of bench_backends.rs, NOT an
edit to it — keep the pure-fluid regression clean): reuse stirred_tank_3d
geometry (parameterized by n), build Solver<D3Q19,T,CpuSimd,LocalPeriodic>
with tank walls+baffles as static solids, per row wire the landed APIs
(GlobalSpec.collision toggle; WaleLes::new().update(&mut solver) each step
before step(); apply_rotating_ibm(&RotatingBody, &DirectForcingConfig) per
step; ParticleSet stepped after step()). Warmup-excluded MLUPS. Publishes the
metric + protocol + a real subset number; de-RED-ifies the METHODOLOGY while
claims-ledger.md:26 correctly stays RED.

Stage B — full ME-4 (blocked): C5 += resolved Allen-Cahn two-phase (MF-γ /
W-VOF, critical path), C6 += scalar ADE active (MF-ε / W-SCAL+REACT), C7 +=
two/four-way particles, C8 = full coupled at D3Q27 fidelity (MF-ζ acceptance).
Only C8 is the M-Star-comparable full number claims-ledger.md:26 gates on.
Stage B rides MF-γ/ε/ζ (QA-sweep-owned); not schedulable until W-VOF unblocks.

Measurement protocol (both stages) — per the loaded-window trap
(lbmflow-whitepaper-benchmark memory; claims-ledger.md:22 1353→2791 artifact):
(1) quiet window (verify machine idle); (2) warmup ≥10 excluded; (3) A/B/A
interleave, two A's must agree within a few % or discard+re-run; (4) 5-run
median, report min/max band never a single number; (5) fixed step budget so
each run ≥ few hundred ms; (6) WALE-update and apply_rotating_ibm MUST be
inside the timed loop for C2/C3/C4 (they are part of the physics-step cost —
do not hoist out).

## 5. Ownership

ME-4 tracks M-F (PLAN.md:165, MF-ζ acceptance); the blocking track W-VOF is
QA-sweep-owned (HANDOFF §4). Stage A touches only landed core APIs + a new
read-only bench example → dispatchable independently NOW without waiting on
W-VOF. Stage B (C5–C8) sequences behind MF-γ/ε/ζ. Recommended: dispatch
Stage A now so the harness + quiet-window protocol are proven on landed
physics; when MF-γ/ε/ζ land, Stage B adds two toggle rows to an already-proven
harness and the full M-Star-comparable ratio drops out of the same protocol.
