# LBMFlow Limitations

This file is a release-facing trust boundary for the current worktree. It states
what LBMFlow supports today and where users should not rely on it yet. For
planned work, see [PLAN.md](PLAN.md); this manifest does not make roadmap
promises.

## 1. Lattice and Boundaries

| Area | Current limitation | Evidence |
|---|---|---|
| D3Q27 open faces | D3Q27 velocity inlet, pressure outlet, outflow, and convective open-face closures are not supported. Open-face validation accepts only 3-unknown D2Q9 and 5-unknown D3Q19 closures; other lattices return `UnsupportedOpenFaceLattice`. D3Q27 has 9 unknown populations per face. | `crates/lbm-core/src/solver.rs:233-240`, `crates/lbm-core/src/solver.rs:349-353`, `crates/lbm-core/src/solver.rs:586-600`, `crates/lbm-core/src/lattice.rs:446-451` |
| Curved walls | Curved-wall helpers are analytic Bouzidi circle and sphere records, plus an explicit low-level link hook for validation. They are not a general geometry importer. | `crates/lbm-core/src/bouzidi.rs:47-56`, `crates/lbm-core/src/bouzidi.rs:105-114`, `crates/lbm-core/src/solver.rs:1865-1896`, `crates/lbm-core/src/solver.rs:1919-1923` |
| Geometry import | Scenario obstacles are built-in primitives (`circle`, `rect`, `sphere`). STL/CAD import is outside the current solver spec; no voxelization import path is exposed by the scenario schema. | `crates/lbm-scenario/src/lib.rs:306-320`, `docs/REQ_STIRRED_REACTOR.md:568-571` |

## 2. Backend Coverage

| Area | Current limitation | Evidence |
|---|---|---|
| Scenario GPU dispatch | Scenario-level GPU dispatch rejects `f64`; for 3D scenarios it rejects multiphase, rotor, particles, non-rest initialization, and force probes. The error string also names non-f32 storage as unsupported in this path. | `crates/lbm-scenario/src/lib.rs:110-112`, `crates/lbm-scenario/src/lib.rs:148-164` |
| Localized sources and face patches on GPU | The backend trait marks localized volume sources and masked face patches as optional; `WgpuBackend` returns `false`, and solver construction rejects non-empty sources or patches with `UnsupportedOnGpu`. | `crates/lbm-core/src/backend.rs:165-170`, `crates/lbm-core/src/gpu/backend.rs:1885-1890`, `crates/lbm-core/src/solver.rs:1378-1388` |
| Gravity performance | Runs with gravity cannot use the chunked fast path. `run()` falls back to `step()` per iteration, and each staged step performs host stage-out, overlays `rho*g`, uploads, runs one backend step, reads back, removes the overlay, and uploads again. | `crates/lbm-core/src/solver.rs:1513-1565`, `crates/lbm-core/src/solver.rs:1587-1605`, `crates/lbm-core/src/solver.rs:1634-1641` |
| GPU availability in tests | GPU tests are optional feature tests and require a native adapter. T14/T16 either require or skip based on adapter availability; default workspace tests do not cover GPU hosts. | `README.md:162-164`, `docs/VALIDATION.md:232-242`, `crates/lbm-core/tests/t14_backend_equiv.rs:40-43`, `crates/lbm-core/tests/t16_fp16_storage.rs:25-33` |

## 3. Precision and Collision Exposure

| Area | Current limitation | Evidence |
|---|---|---|
| FP16 storage | FP16 is a capacity/throughput mode, not a validation-grade long-transient reference mode. Distribution storage narrows to f16 while arithmetic remains f32; steady flows re-converge, but long transients accumulate storage rounding. | `crates/lbm-core/src/gpu/backend.rs:142-148`, `crates/lbm-core/tests/t16_fp16_storage.rs:15-23`, `docs/PHYSICS.md:313-324`, `docs/PERFORMANCE.md:59-62` |
| Scenario schema | Scenario JSON exposes `collision: bgk | trt | cumulant` and `compute.storage: f32 | f16`, but with narrow honored paths: cumulant only on 3D D3Q19 CPU scenarios; f16 only on 2D D2Q9 GPU scenarios (SHADER_F16 adapter required). All other combinations are rejected with explicit errors — no silent fallback. Scenario-level lattice selection (D3Q27) and MPI hints remain unexposed. | `crates/lbm-scenario/src/lib.rs` (`CollisionSpec`, `StorageSpec`, `CUMULANT_*`/`GPU_F16_*` errors), manifest `provenance` in `crates/lbm-cli/src/runner.rs` |

## 4. Particles

| Area | Current limitation | Evidence |
|---|---|---|
| Coupling model | Particles are one-way Lagrangian particles. They feel sampled velocity, buoyancy-reduced gravity, and Schiller-Naumann drag, but they do not apply reaction forces to the fluid. Two-way/four-way coupling, Saffman, Basset, Faxen, collision models, and stochastic LES dispersion are not implemented. | `crates/lbm-core/src/particles.rs:1-15` |
| Schiller-Naumann range | Drag correction is valid for `Re_p <= 800`; exceeding it returns a `ParticleError` (particle index + offending Re) — runs do not silently continue outside the correlation's validity domain. | `crates/lbm-core/src/particles.rs` (`SCHILLER_NAUMANN_RE_MAX`), PHYSICS.md 2026-07-07 entry |
| Near-wall sampling | `sample_grid` clamps sample positions to grid bounds, uses solid-neighbor velocity as zero, and returns the solid flag of the clamped lower node. Near-wall and out-of-domain particle samples therefore need interpretation as clamped grid samples, not extrapolated wall-resolved velocities. | `crates/lbm-core/src/particles.rs:276-285`, `crates/lbm-core/src/particles.rs:302-321`, `crates/lbm-core/src/particles.rs:324-331` |

## 5. LES

| Area | Current limitation | Evidence |
|---|---|---|
| WALE scope | WALE is landed as a solver-level eddy-viscosity driver. It computes `nu_t`, converts it directly to an `omega_plus` field, and installs that field for the next collision, giving a one-step lag. | `crates/lbm-core/src/les.rs:1-7`, `crates/lbm-core/src/les.rs:49-57`, `crates/lbm-core/src/les.rs:104-113`, `crates/lbm-core/src/solver.rs:2483-2518` |
| Remaining LES product treatment | tau_eff upper clipping with mandatory diagnostics landed 2026-07-07 (explicit config, default off — see PHYSICS.md entry). Still missing: y+ wall-function or wall-fitted near-wall handling (design spec in docs/proposals/LES_WALL_TREATMENT_SPEC.md) and turbulence-predictive acceptance beyond the channel Re_tau=180 characterization. | `docs/REQ_STIRRED_REACTOR.md:202-210`, `crates/lbm-core/src/les.rs` (`WaleLes` clipping + diagnostics) |

## 6. Multiphase

| Area | Current limitation | Evidence |
|---|---|---|
| Validated scope | Current multiphase validation covers Shan-Chen SCMP flat interface, Laplace law, contact angle, and two-component MCMP Rayleigh-Taylor growth within the documented bands. | `docs/VALIDATION.md:167-200`, `docs/MULTIPHASE_DESIGN.md:3-8`, `docs/MULTIPHASE_DESIGN.md:58-74`, `README.md:46-48` |
| Method limits | The Shan-Chen design documents density-ratio and spurious-current weaknesses and coupling between surface tension and density ratio. | `docs/MULTIPHASE_DESIGN.md:19-23`, `crates/lbm-core/src/compat/multiphase.rs:42-49` |
| Free surface and aeration | Conservative Allen-Cahn free surface (`W-VOF`) is pending. High-density-ratio gas-liquid and stirred-tank aeration acceptance remains T17 work with bands frozen after implementation, not a validated release claim today. | `docs/REQ_STIRRED_REACTOR.md:30-36`, `docs/REQ_STIRRED_REACTOR.md:229-233`, `docs/REQ_STIRRED_REACTOR.md:537-542`, `docs/REQ_STIRRED_REACTOR.md:590-594`, `docs/VALIDATION.md:286-311` |

## 7. Checkpoint and Restart

| Area | Current limitation | Evidence |
|---|---|---|
| Rank scope | Checkpoints are single-rank only. Saving with more than one local rank/part returns `CKPT_UNSUPPORTED`; loading also requires exactly one rank and one local part. | `crates/lbm-core/src/solver.rs:2556-2567`, `crates/lbm-core/src/solver.rs:2786-2790` |
| Serialized state | The checkpoint writes populations, stale stash, moments, solid mask, and optional force field. The manifest explicitly reserves `rng`, `particles`, and `stats` as `false`, so RNG state, particle state, and statistics accumulators are not serialized today. | `crates/lbm-core/src/solver.rs:2575-2600`, `crates/lbm-core/src/solver.rs:2648-2665`, `crates/lbm-core/src/solver.rs:2853-2888` |
| Distributed restart | Distributed checkpoint/restart is still a solver-improvement item, not the landed checkpoint path. | `docs/SOLVER_IMPROVEMENT_SPEC.md:190-194` |

## 8. MPI and Scale

| Area | Current limitation | Evidence |
|---|---|---|
| Functional coverage | MPI is verified for multi-rank execution within a single node; true multi-node weak scaling awaits cluster measurement. | `docs/MPI_GUIDE.md:7-17`, `docs/ARCHITECTURE_V2.md:157-160`, `docs/VALIDATION.md:220-227` |
| Multi-node performance claim | The 64-rank multi-node weak-scaling claim is RED in the release status table. | `README.md:138-146`, `docs/PERFORMANCE.md:70-71`, `docs/PLAN.md:124` |
| Memory scaling | `MpiSolver::new` still requires global compact arrays on all ranks, and the solver-improvement spec identifies this global-array replication as an OOM blocker at large grids. | `docs/SOLVER_IMPROVEMENT_SPEC.md:121-134` |
| Parallel I/O | MPI output is rank-0 gather only; parallel VTK/HDF5 is not started. | `docs/MPI_GUIDE.md:18-26`, `docs/SOLVER_IMPROVEMENT_SPEC.md:149-156` |

## 9. Moving Bodies

| Area | Current limitation | Evidence |
|---|---|---|
| Prescribed rigid rotation | The native IBM body is a marker set with fixed center and angular velocity; target velocity is prescribed as `U = Omega x r`. The compat rotor similarly prescribes solid-body target velocity and a ramped angular speed. | `crates/lbm-core/src/rotating_ibm.rs:1-8`, `crates/lbm-core/src/rotating_ibm.rs:22-32`, `crates/lbm-core/src/rotating_ibm.rs:89-101`, `crates/lbm-core/src/compat/rotor.rs:192-210` |
| Diagnostics, not structural FSI | The current moving-body implementation reports slip, reaction torque, force, and momentum-spreading diagnostics. It does not expose structural degrees of freedom, deformation state, or added-mass coupling analysis; current M-F requirements describe this landed path as IBM-inertial rotating-body forcing, with MRF and overset as separate reference/relaxation tracks. | `crates/lbm-core/src/rotating_ibm.rs:132-153`, `crates/lbm-core/src/solver.rs:2003-2015`, `crates/lbm-core/src/solver.rs:2148-2173`, `docs/REQ_STIRRED_REACTOR.md:52-55`, `docs/REQ_STIRRED_REACTOR.md:213-222` |

