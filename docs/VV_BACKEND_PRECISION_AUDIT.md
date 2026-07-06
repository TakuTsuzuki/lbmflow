# V&V Backend / Precision Audit

Suborder: G, V&V-BACKEND  
Date: 2026-07-06  
Worktree: `/Users/taku/projects/lbmflow-wt-cx-vv-backend`  
Base commit audited: `4eac49e2aebccfa52283c29cd7f22aabba3120f5`

## Scope and Rule

This is an audit of backend, precision, GPU, MPI, and partition equivalence
coverage. It does not add new physics, new constants, new tolerances, or new
acceptance bands.

Equivalence evidence is not the same as absolute physical validation. A
CPU-vs-GPU or partition-vs-monolithic match can prove that two paths implement
the same computation, but it cannot prove that the shared computation is
physically correct. Absolute physics gates remain T1-T12, T15, T17/T18, and
future GPU-direct sentinels where called out below.

## Evidence From This Session

| Command | Result | Evidence notes |
|---|---:|---|
| `cargo build --workspace --release` | PASS | Finished release profile in 39.02 s. |
| `cargo test --workspace --release` | PASS | Full default workspace suite exited 0. Notable backend evidence inside: `backend_simd_equiv` 21 passed; `t13_split_invariance` 9 passed; T14/T16 GPU files compiled as 0 tests without `gpu`; T15 f32 3D 3 passed; T15.5 default sanity passed after 275.32 s; Ghia Re100/Re400 passed. |
| `cargo test --release -p lbm-core --test backend_simd_equiv` | PASS | 21 passed, 0 failed. Covers CpuScalar vs CpuSimd in f64/f32 over D2Q9 plus D3Q19 and D3Q27 TGV/duct cases and CpuSimd split cases. |
| `cargo test --release -p lbm-core t13` | PASS | T13 filter exited 0. `t13_adversarial`: 6 passed, 1 ignored long-run. `t13_split_invariance`: 9 passed, including D3Q19 and D3Q27 2x2x2 split invariance. |
| `cargo build -p lbm-core --release --features gpu` | PASS | GPU feature compiled successfully. |
| `cargo test --release -p lbm-core --features gpu --test t14_adversarial near -- --nocapture` | PASS | 1 passed. Runtime GPU adapter was available; near-MAX_SPEED TGV field deltas stayed <= 1.644e-6. |
| `cargo test --release -p lbm-core --features gpu --test t14_backend_equiv -- --nocapture` | PASS | 8 passed. D2Q9 CPU-vs-wgpu runtime equivalence passed; pressure channel used the documented 1e-4 relaxed line with CPU-vs-CPU 1-ulp control. |
| `cargo test --release -p lbm-core --features gpu --test t14_3d_backend_equiv -- --nocapture` | PASS | 3 passed. D3Q19 GPU runtime equivalence passed for TGV3D, cavity3D, and open3D+force. |
| `cargo test --release -p lbm-core --features gpu cumulant_gpu_matches_cpu_measured_tgv3d_tolerance -- --nocapture` | PASS | D3Q19 delta 1.5497207641601563e-6; D3Q27 delta 3.4570693969726563e-6 over 200-step f32 cumulant TGV3D. |
| `PATH="$HOME/.local/openmpi/bin:$PATH" cargo build -p lbm-core --release --features mpi` | PASS | MPI feature compiled with native arm64 Open MPI 5.0.9. |
| `./scripts/test_mpi.sh` | BENCH-PENDING | Built the MPI example, then every `mpirun` failed before executing the test because PRTE could not bind sockets: `Operation not permitted`; sandbox/runtime environment issue, not an MPI equivalence result. |

All commands emitted the existing `compat/domain.rs:is_wall` dead-code warning.
GPU test builds also emitted deprecation warnings for the `GpuSolver` wrapper.

## Coverage Inventory

| Axis | Existing coverage | Status | Audit judgment |
|---|---|---|---|
| CpuScalar vs CpuSimd | `tests/backend_simd_equiv.rs`: f64/f32 D2Q9 scenarios; D3Q19 and D3Q27 3D TGV; 3D duct; mass, momentum, probed_force diagnostics; CpuSimd split-vs-CpuScalar cases. | VERIFIED | Strong equivalence coverage, including diagnostics beyond fields. Not an independent physics proof. |
| Monolithic vs in-process partition | `tests/t13_split_invariance.rs`, `tests/t13_adversarial.rs`: D2Q9, D3Q19, D3Q27, two-pass, seams, probes, Shan-Chen, open/moving boundaries. | VERIFIED | Strong field bit-match and diagnostic-tolerance evidence. The ignored 20k cavity long-run remains deeper-horizon pending. |
| MPI vs monolithic | `examples/mpi_t13.rs` + `scripts/test_mpi.sh`: ranks 1/2/3/4/5/6/8, 2D cases, D3Q19 3D at n=6/8, mismatch-nu negative test. | BENCH-PENDING | Feature build passed. Runtime could not start in this sandbox because Open MPI socket bind was denied. Do not claim MPI runtime validation from this session. |
| CPU vs wgpu D2Q9 | `tests/t14_backend_equiv.rs` and `tests/t14_adversarial.rs`: TGV, cavity, channel, cylinder+probe, force field, convective outflow, pressure BC sensitivity, near-MAX_SPEED attack. | VERIFIED | Runtime GPU evidence exists in this session. Mostly relative equivalence; pressure case has an explicit 1-ulp sensitivity control. |
| CPU vs wgpu D3Q19 | `tests/t14_3d_backend_equiv.rs`: TGV3D, cavity3D, open faces with body force. `t14_wale_gpu_equiv.rs` exists but was not run in this session. | VERIFIED-PARTIAL | Basic D3Q19 runtime equivalence passed. WALE GPU equivalence was not run here. |
| CPU vs wgpu D3Q27 | GPU unit test `cumulant_gpu_matches_cpu_measured_tgv3d_tolerance` covers D3Q27 cumulant TGV3D. D3Q27 open faces are rejected by construction. | VERIFIED-ONLY | Useful kernel-level/runtime equivalence for D3Q27 cumulant, but not a public T14-style D3Q27 scenario matrix. |
| f64 vs f32 | Core absolute tests include T6 f32 conservation and `tests/t15_3d_f32.rs`; `backend_simd_equiv` covers f32 backend equivalence. | VALIDATED/PARTIAL | D2Q9 f32 is well covered by T6 and standard physics validations; D3Q19 f32 has T15 f32 product-path sentinels. D3Q27 f32 absolute physics coverage is thinner. |
| FP16 storage | `tests/t16_fp16_storage.rs` has ignored GPU-heavy TGV and cavity degradation tests with frozen bands recorded in `PHYSICS.md`. | BENCH-PENDING | Implemented and documented, but heavy ignored T16 tests were not run in this session. Default suite does not validate FP16 runtime. |
| D2Q9 | T1-T14, T16/T18 pieces, CPU/GPU, partition. | VALIDATED | Best-covered lattice. |
| D3Q19 | T15/T15.5, T15 f32, T13, T14 3D GPU. | VALIDATED/PARTIAL | Strong CPU physics and basic GPU equivalence. GPU direct absolute physics sentinels are still needed. |
| D3Q27 | Lattice invariants, CPU smoke, T13, T15.4 D3Q27 TGV3D, CpuSimd equivalence, D3Q27 cumulant GPU unit. Open faces rejected. | VERIFIED-ONLY/PARTIAL | Good early-stage basis, but no broad D3Q27 physical validation matrix and no D3Q27 open-boundary physics. |

## Relative-Only Cases Lacking Absolute Physics Validation

1. GPU D2Q9 T14 verifies CPU-vs-wgpu agreement over representative scenarios, but
   direct GPU TGV convergence order, direct GPU Poiseuille profile, direct GPU
   cavity RMS vs Ghia, and direct GPU conservation sentinels are not first-class
   absolute acceptance tests.
2. GPU D3Q19 T14 verifies CPU-vs-wgpu agreement for short 3D cases, but it does
   not independently prove T15.2 duct, T15.3 sphere drag, T15.4 convergence, or
   T15.5 cavity on the GPU path.
3. D3Q27 GPU coverage is a D3Q27 cumulant TGV3D CPU-vs-GPU delta test, not an
   absolute physics suite.
4. FP16 storage has frozen bands in ignored tests and `PHYSICS.md`; because the
   tests are ignored and GPU-heavy, the default suite cannot be used to claim
   FP16 validation.
5. MPI runtime evidence depends on `mpirun`. In this session, the runtime never
   entered the solver due sandbox socket restrictions.

## Proposed Absolute GPU Sentinels

These should be added as GPU-feature tests with `LBM_REQUIRE_GPU=1` converting
adapter absence from skip to fail:

1. GPU TGV order: run D2Q9 f32 directly on Wgpu at N=32/64 using the T1
   pressure-consistent initialization; require order >= 1.7 and viscosity
   within the existing T1 band.
2. GPU Poiseuille profile: run D2Q9 TRT body-force Poiseuille directly on Wgpu;
   compare to the half-way-wall analytic profile and check mass/momentum drift.
3. GPU cavity RMS: run a shorter Re=100 cavity sentinel directly on Wgpu and
   compare centerlines to the existing Ghia table band adjusted only by recorded
   f32 measurements, not post-hoc tuning.
4. GPU conservation: periodic box mass conservation and uniform-force momentum
   growth directly on Wgpu, using the T6 f32 lines.
5. GPU D3Q19 direct T15 slice: at minimum TGV3D decay rate plus closed-box mass
   drift directly on Wgpu.
6. FP16 storage smoke in the non-ignored GPU suite: a smaller steady-flow case
   that finishes quickly and proves the f16 path is exercised, with the heavy
   T16 cases remaining ignored for characterization.

## Diagnostic Invariance Beyond Fields

Current tests already include several non-field diagnostics:

- `backend_simd_equiv`: total mass, total momentum, probed force.
- `t13_split_invariance`: total mass, momentum, probed force; population planes.
- `mpi_t13`: total mass, momentum, probed force, nonfinite count, all f planes.
- `t14_backend_equiv`: total mass, momentum, probed force in D2Q9 GPU cases.
- `t15_3d_f32`: f32 D3Q19 mass drift.
- T18/T13 deposition tests: partition-stable particle/deposition records.

Gaps:

- Torque invariance is covered in rotating/IBM-focused tests, not in the backend
  equivalence matrices.
- Scalar totals are not a first-class backend-equivalence diagnostic outside the
  Shan-Chen/force-field paths and future T17/T18 scopes.
- GPU D3Q19 T14 does not currently compare diagnostics as comprehensively as the
  D2Q9 T14 helper.

## BENCH-PENDING

- MPI runtime: blocked by sandbox socket bind denial in `mpirun`, despite native
  Open MPI and successful feature build.
- FP16 storage characterization: implemented and documented, but ignored heavy
  tests were not run here.
- GPU performance benchmarks: not run here; this audit only ran correctness
  tests. Use `cargo run -p lbm-core --release --features gpu --example
  bench_gpu` for MLUPS evidence.
- Full heavy validation: `cargo test --release -- --include-ignored` was not run
  in this suborder.

## Merge Recommendation

Merge recommendation: docs-only audit can merge after review. It changes no
runtime behavior and adds no physics. Do not use this document to claim MPI
runtime validation or FP16 heavy validation from this session. GPU runtime
equivalence for the executed T14/T14-3D/D3Q27-cumulant cases is supported by
session evidence.
