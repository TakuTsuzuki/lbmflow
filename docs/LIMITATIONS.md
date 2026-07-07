# LBMFlow Limitations — Bioprocess CFD

Lifecycle: living (release-facing trust boundary; kept in sync with the
capability registry defined in BCFD-002).

This file is machine-readable in spirit — the capability registry in
`crates/lbm-cli` should be able to regenerate it, and drift between the
two is a merge-queue defect. The registry statuses used below map to the
credibility tiers in
[SPEC_BIOPROCESS_CORE.md](SPEC_BIOPROCESS_CORE.md) §3.

**Status legend:**

- `Unsupported` — the combination is rejected with a structured error.
- `Experimental` — code path runs, but not validated to a bioprocess band.
  Results are demos / screening; do not use in decisions.
- `Engineering` — bioprocess-band validation green; usable for
  design-of-experiments and internal decisions; no calibration / holdout
  pinned.
- `EvidenceBlocked` — Engineering ceiling reached; needs
  calibration + holdout + UQ + sensitivity for the specific QOI.
- `EvidenceReady` — BCFD-091 gate passed for this QOI.

Every row cites either a code location (source of truth) or a doc entry.

## 1. Bioprocess capabilities (BCFD-002 registry snapshot)

The table below is keyed by capability id and must match
`crates/lbm-cli/src/capabilities.rs`. `crates/lbm-cli/tests/capabilities_drift_guard.rs`
fails if a status or tier ceiling differs from the registry.

| Capability id | Capability | Status | Tier ceiling | Reason | Evidence / release note |
|---|---|---|---|---|---|
| `single_phase_stirred_tank` | Single-phase stirred tank | Experimental | Screening | BCFD-030 runner path implemented; VB-01 not yet green. | `docs/PLAN.md#bcfd-030`; closed 2026-07-08 via BCFD-030. |
| `rotating_ibm` | Rotating IBM impeller | Experimental | Screening | BCFD-030 integration implemented; stirred-tank validation pending. | `docs/PLAN.md#bcfd-021`; closed 2026-07-08 via BCFD-021 / BCFD-030. |
| `passive_scalar` | Passive scalar transport | Experimental | Screening | BCFD-034 ADE path implemented; VB-02 reducer tests are green, but registry promotion still waits on integrated product artefacts. | `docs/PLAN.md#bcfd-034`; closed 2026-07-08 via BCFD-034. |
| `phase_field_vof` | Phase-field VOF | Unsupported | Screening | Conservative Allen-Cahn path is not yet implemented in the registry. | `docs/PLAN.md#bcfd-040`. |
| `oxygen_kla` | Oxygen transport and kLa | Experimental | Screening | BCFD-050..053 oxygen scalar, Henry flux, kLa fit, and OUR hooks implemented; VB-06 synthetic reducer tests are green, but calibrated kL and holdout artefacts are absent. | `docs/PLAN.md#bcfd-050`; closed 2026-07-08 via BCFD-050..053. |
| `point_bubbles` | Point bubbles | Unsupported | Screening | Point-bubble entity store is not yet implemented. | `docs/PLAN.md#bcfd-070`. |
| `pbm` | Population balance model | Unsupported | Screening | PBM bins and kernels are not yet implemented. | `docs/PLAN.md#bcfd-073`. |
| `cell_exposure` | Cell and microcarrier exposure | Unsupported | Screening | Cell tracer and exposure QOIs are not yet implemented in the registry. | `docs/PLAN.md#bcfd-060`. |
| `evidence_tier_report` | Evidence-tier report | Unsupported | Evidence | Evidence gate is implemented, but required artefacts are absent: validation matrix pass, calibration/holdout separation, mesh/time-step sensitivity, QOI uncertainty interval, and limitation report. | `docs/PLAN.md#bcfd-091`. |

Promotion rule: a capability is `Experimental` when its code path exists
but lacks the full bioprocess validation package; `Engineering` requires
registry promotion after the relevant VB group and numerical sensitivity
records are green; `EvidenceReady` is per QOI and only after BCFD-091
returns ready. As of 2026-07-08, no QOI is `EvidenceReady`.

## 2. Legacy LBM capabilities preserved (with demo warning)

Legacy scenario paths continue to run and their tests remain green.
`lbm presets run <name>` emits `legacy LBM demo preset; not bioprocess
decision-grade` to stderr.

### 2.1 Lattice and boundaries

| Area | Status | Notes |
|---|---|---|
| D3Q27 open faces | Engineering (technical); Unsupported (bioprocess) | CPU full boundary coverage; GPU rejects. `crates/lbm-core/tests/d3q27_open_bc.rs`. |
| Curved walls (Bouzidi) | Engineering (technical); analytic-only geometries | `crates/lbm-core/src/bouzidi.rs`. Non-analytic geometry needs BCFD-023 STL import. |
| Geometry import | Unsupported | BCFD-023 optional feature; not on M0 critical path. |

### 2.2 Backends

| Area | Status | Notes |
|---|---|---|
| Scenario GPU dispatch | Unsupported for bioprocess coupled physics | Rejects multiphase, rotor, particles, non-rest initialization, force probes — the entire bioprocess coupled path. 2D f32 GPU scenarios run for demo purposes. |
| Localized GPU sources / face patches | Unsupported | Backend trait marks optional; wgpu returns false. |
| Gravity | Engineering (technical) | Backend-side body-force composition landed pre-pivot. |
| GPU availability in tests | Optional feature | `--features gpu`; not covered by workspace default. |

### 2.3 Precision and collision

| Area | Status | Notes |
|---|---|---|
| FP16 storage | Unsupported (validation-grade) | Capacity / throughput mode only. Long transients accumulate storage rounding. Never used for a QOI reported in a bioprocess decision. |
| Scenario schema (lattice, collision, storage) | Engineering (technical) for legacy | Narrow honored paths: D3Q27 CPU-only; cumulant 3D CPU-only; f16 only 2D D2Q9 GPU. Otherwise rejected with explicit errors — no silent fallback. |

### 2.4 Particles

| Area | Status | Notes |
|---|---|---|
| Coupling model | Engineering (technical, one-way) | `ParticleSet` is deterministic one-way Lagrangian. No reaction force to fluid. Two-way (BCFD-063), four-way, Saffman, Basset, Faxen, collision, LES stochastic dispersion are not implemented. |
| Schiller-Naumann range | Enforced Re_p ≤ 800 | Violations return `ParticleError` (particle index + offending Re); no silent extrapolation. |
| Near-wall sampling | Clamped grid samples | Not extrapolated wall-resolved velocities; interpret as clamped. |

### 2.5 LES

| Area | Status | Notes |
|---|---|---|
| WALE eddy viscosity | Engineering (technical, channel Re_τ=180) | Solver-level driver; one-step lag. `crates/lbm-core/src/les.rs`. |
| ν_t clipping | Diagnostic guard only | Default off. When enabled, `clipped_fraction` and `max_nu_t_before_clipping` are required in every reported validation. |
| Wall treatment | Unsupported | y+ wall-function / wall-fitted near-wall handling not implemented; design spec in `docs/archive/2026-07-07-pivot/proposals/LES_WALL_TREATMENT_SPEC.md`. |

### 2.6 Multiphase (legacy)

| Area | Status | Notes |
|---|---|---|
| Shan-Chen SCMP flat interface, Laplace law, contact angle | Demo | Documented density-ratio and spurious-current weaknesses. NOT production gas-liquid — see BCFD-040..048. |
| Shan-Chen MCMP Rayleigh-Taylor | Demo | — |
| Free surface / high-density-ratio gas-liquid / stirred-tank aeration | Unsupported | Belongs to BCFD-045..048. |

### 2.7 Checkpoint and restart

| Area | Status | Notes |
|---|---|---|
| Solver-state checkpoint (v3 format) | Engineering (technical) | Populations + stale buffer + moments + solid mask + optional force field + scalar distributions + phase field + QOI accumulators; per-rank MPI files with rank-0 manifest; strict layout / version guards. |
| Serialized state coverage | Partial | Scalar, phase, and QOI-stat sections are serialized when present; future cell-tracer / bubble / RNG sections have scaffolded traits but remain absent until BCFD-060 / BCFD-070 producers land. |
| Large-scale resilience | Partial | Parallel field slabs are supported for MPI source fields; failure recovery, atomic publish, and partial-write repair across ranks are still not implemented. |

### 2.8 MPI and scale

| Area | Status | Notes |
|---|---|---|
| Multi-rank single-node | Engineering (technical) | `crates/lbm-core/src/dist.rs`. |
| Multi-node weak scaling ≥80% @ 64 rank | Unsupported | RED pre-pivot; deferred until BCFD-100 + cluster access. |
| Memory scaling | Engineering (technical) | `MpiSolver::new_local` builds masks/material samples from owned-cell callbacks; legacy `MpiSolver::new` remains small-scale only. |
| Parallel I/O | Engineering (technical) | MPI ranks can write per-rank binary field slabs plus a manifest for velocity, phase `phi`, oxygen/scalar `C`, shear rate, and gas holdup; legacy rank-0 gather remains for small validation cases. |

### 2.9 Moving bodies

| Area | Status | Notes |
|---|---|---|
| Prescribed rigid rotation | Engineering (technical) | Marker set with fixed centre and ω; target U = ω × r. `crates/lbm-core/src/rotating_ibm.rs`. |
| Structural FSI / deformation / added mass DOF | Unsupported | Diagnostics (slip, reaction torque, force, momentum spreading) only. MRF and overset are out of scope. |

## 3. Explicit non-goals

The following are outside the product mission and will not be implemented
even if the underlying code makes them feasible:

- General-purpose OpenFOAM parity, arbitrary FVM numerics port.
- General CAD mesher (BCFD-023 is a minimal STL voxeliser for
  screening tier; not a mesher).
- Compressible CFD, combustion, general solid mechanics.
- Arbitrary chemistry kinetics engine (reaction hooks are for OUR /
  simple source terms only, BCFD-053).
- Fully general non-structured mesh solver.
- Public web GUI ahead of BCFD-081 report generator.
- GMP / CMC readiness claim without the BCFD-091 evidence gate.
