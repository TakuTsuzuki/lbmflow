# Solver Improvement Specification (Full Review 2026-07-05)

> **main integration note (2026-07-05 PM)**: This document was written on the review branch
> `claude/amazing-mirzakhani-4060d3` (based on commit 84abaa3, before V1 retirement).
> Since main has already retired V1, read it with the following path substitutions:
> `crates/lbm-core2` → `crates/lbm-core` (renamed), old V1's `crates/lbm-core/src/*` →
> `crates/lbm-core/src/compat/*` (consolidated into the facade; standalone V1 has been deleted).
> **A-1 (S0) is already resolved on main** (sync-tests.sh was ported to perl and then deleted upon V1 retirement;
> the duplicated suite has been promoted to the canonical tests that reference compat directly). B-4's "compat migration" was also
> carried out with V1 retirement. Experiments can be re-run in `scripts/spec-experiments/` (paths already translated) —
> E2/E7 have been confirmed to match the spec's numbers on renamed main.

> **Status: v1 (experimentally validated, actionable version)**. All v0 claims have been
> validated in §3 by experiments E1–E10. **Zero items were refuted**; 2 descriptive corrections (the symptom character of A-3, the numbers of A-6).
> **A-1 (S0) was carried out on this branch concurrently with validation** (perl port of sync-tests.sh plus
> a postcondition guard; the regenerated suite was confirmed fully green).
> Experiment code can be re-run in `scripts/spec-experiments/` (e.g. `cargo run --release e2`).

- Target commit: 84abaa3 (equivalent to main, branch claude/amazing-mirzakhani-4060d3)
- Scope: `crates/lbm-core` (V1, frozen reference implementation), `crates/lbm-core2` (V2:
  kernels / lattice / fields / solver / backend / halo / subdomain / dist / gpu / compat),
  validation and benchmark infrastructure (tests / VALIDATION.md / scripts / CI)
- Review method: 6 independent review streams (V2 physics kernels / V2 structure / GPU / MPI / V1 core /
  V&V infrastructure) were conducted in parallel, and the PM backed S0/S1-class findings with real code. Line numbers are as of the target commit.
- Severity: **S0** = correctness error (e.g. tests give false assurance) / **S1** = high risk
  (latent bug, scaling failure, significant design defect) / **S2** = improvement opportunity / **S3** = minor.
- Effort: S = a few hours / M = 1 day / L = several days.

---

## 0. Executive Summary

**There is no S0 in the correctness of the physics kernels.** D3Q19 Zou–He has been carried through
literature cross-check and hand-calculation verification against Hecht & Harting (2010), and matches. The TRT decomposition of Guo forcing, the constant folding of deviation storage,
the sufficiency of halo exchange (exact match with the unknown set), and the single-writer guarantee of GPU push-type fusion have all
been confirmed. **The issues are concentrated in "entry, structure, operations, and validation infrastructure"**:

1. **1 S0 in the validation infrastructure**: sync-tests.sh's sed is not BSD-sed compatible (`\b`) and fails silently, so
   the 57 tests of the "V2 validation suite" are actually re-running V1 (the compat layer has never received physics validation
   = R5 not demonstrated). The fix exists on the unmerged `v1-retirement` branch.
2. **Absence of entry guards**: The V2 native path (the main path for 3D/MPI) has no configuration-validation layer, and
   "uncovered faces," "periodic×open coaxial," "ν=0," etc. become non-physical computations without symptoms. On the V1 side too,
   NaN velocities pass through validation (confirmed by measurement).
3. **M-E's type-level blocker**: Due to the `Fields = SoaFields<T>` binding of `Solver`/`MpiSolver`, the
   GPU backend does not mount on the orchestrator, and the step sequence is duplicated in GpuSolver.
   Multi-GPU / MPI+GPU cannot be started with the current structure (independently identified by the structure review and the GPU review).
4. **Triple implementation of physics**: The Shan–Chen force runs in parallel in 3 places — V1 / compat / V2 native —
   and the feature matrix is asymmetric per axis (contact angle and two-component remain 2D/CPU-only with fixation risk).
5. **Scaling failures**: MPI setup replicates global arrays on all ranks (at 10⁹ grid,
   wall_u ≈ 24 GB/rank), no communication overlap, only rank-0 serial gather,
   no checkpointing. In addition, rayon may launch while `MPI_THREAD_SINGLE` is still declared.
6. **Validation gaps**: No CI (all quality claims are manual snapshots), the MPI path executes 0 lines under cargo test,
   GPU has only CPU-relative equivalence (zero absolute physics validation), and f32×3D is a product path yet unvalidated.

---

## 1. Reviewed Scope and Confirmed Items (No Problems)

The following are foundations confirmed not to need modification. They are also the "must not break" list during modifications.

- **D2Q9/D3Q19 tables**: C/W/OPP/TRT pairs/FACE_UNKNOWNS are all const-fn derived from C, so
  drift is structurally impossible. 0th–4th order moment isotropy 1e-15, bit-lock inspected against the V1 tables.
- **D3Q19 face Zou–He**: The closure ρ(1−u·n) = S0+2S⁻, NEBB reconstruction, and the signs of the 6 terms of the tangential correction N_k were
  cross-checked with the literature and match. Exact satisfaction of mass and normal momentum was confirmed by derivation.
- **Guo forcing**: TRT-consistent symmetric/antisymmetric decomposition (cp=1−ω⁺/2, cm=1−ω⁻/2), the F/2 correction of moments and
  Reduction, identical in definition to V1.
- **Constant folding of deviation storage**: half-way BB invariance, +2w physicalization of the probe, constant cancellation of convective mass pinning,
  weight neutrality of outflow — consistent on all paths.
- **Halo exchange**: the directions read from the halo in pull = exact match with the unknown set (has an exhaustive-inspection test).
  Corner satisfaction of the x→y→z two-phase forward, no intra-phase hazards.
- **MPI protocol**: All-Irecv-first-post → Isend → wait for deadlock freedom,
  no wildcard receives, unique tags, communicator separation, POD raw-byte transfer.
  There is field bit-match demonstration in T13-MPI (2D corner, 3D 2×2×2, including ψ exchange).
- **GPU push-type fusion**: single writer for every (q, cell) slot, consistency of parity management on all paths,
  term-by-term match of WGSL and kernels.rs (including combination order). f64 is rejected at compile time, so no silent degrade.
- **V1 fusion path**: exclusivity of band boundaries, the unsafe row-access contract, periodic wrap of copy_span /
  open-end slot preservation confirmed for all cx cases. No substantive drift in the compat layer's copy
  (differences are only doc/import/unused_mut).
- **Test determinism**: no randomness or environment dependence, no flake source. T13 is at the field-level
  `assert_eq!(d, 0.0)` bit-match level. T14 pressure-BC relaxation has a 1-ulp control test.

---

## 2. Improvement Items

### WP-A: Correctness / Entry Guards (immediate; premise is that all items pass existing tests green without modification)

#### A-1 [S0] Substitution fix of sync-tests.sh and recurrence guard
- Target: `scripts/sync-tests.sh:42`, generation target `crates/lbm-core2/tests/` (17 files)
- Current state: `\b` is invalid under BSD sed → substitution-free copy. All 17 generated files remain `use lbm_core::`,
  0 compat references (confirmed by PM measurement). The fix commit 622bbb2 (perl port) exists only on the
  unmerged `v1-retirement`.
- Improvement: Apply the equivalent of 622bbb2 (perl substitution) to the main line, and add a postcondition to the script
  (exit 1 if `use lbm_core::` remains in the output or the substitution count is 0). Add a static guard test that inspects that
  files with the generation header import `lbm_core2::compat`.
- Acceptance: `grep -rl "AUTO-GENERATED" crates/lbm-core2/tests | xargs grep -L "lbm_core2::compat"`
  is empty. The regenerated suite is fully green via compat (→ experiment E1 pre-validates this).
- Effort: S / Experiment: **E1 = confirmed, already carried out** (applied the perl port + postcondition guard, regenerated 16 files
  → `cargo test -p lbm-core2 --release` **whole suite green (exit 0)**.
  Physics suites including contact-angle frozen values, RT, and open-BC systems ran and passed for the first time via compat =
  **first demonstration of R5**. Remaining work is only adding the static guard test)

#### A-2 [S1] Make the configuration validation of V1+compat NaN-safe
- Target: `crates/lbm-core/src/domain.rs:306-310`, `crates/lbm-core2/src/compat/domain.rs` (same code)
- Current state: `if s > MAX_SPEED` is false for s=NaN (passes through). A NaN inlet makes the field NaN in 3 steps;
  a NaN MovingWall **silently degrades into a static wall** due to the comparison failure of rim-velocity selection (measured in review).
  `Trt { magic }` and `force` are also unvalidated. The rho check directly below is the NaN-safe form `!(x > 0.0)`, so
  the safe idiom itself already exists in the codebase.
- Improvement: Invert the velocity check to `if !(s <= MAX_SPEED)`. Add `magic > 0` and `is_finite()` of `force`/`u` to
  `validate()`. The same patch for V1 and compat (maintaining text identity).
- Acceptance: `build()` with NaN/inf u, force, or magic≤0 all return Err. Legal configurations are bit-invariant.
- Effort: S / Experiment: **E6 = confirmed** (NaN inlet: build()=Ok remains, and 42 cells have non-finite rho after 3 steps.
  NaN MovingWall: build()=Ok, and the field bit-matches a static wall — rigorously demonstrating silent static-wall degradation)

#### A-3 [S1] Reject configurations of Outflow/ConvectiveOutflow × solid adjacency's silent mass leak
- Target: `crates/lbm-core/src/sim.rs:765-774,667-669`, `crates/lbm-wasm/src/lib.rs:240-242`
- Current state: If the inner neighbor of an open-end cell is solid, the BC is skipped (confirmed in code by
  `solid[i] || solid[j]`), and the unknown slot is frozen forever at its initial value. It is a production path reachable from
  GUI painting (which rejects only the outermost 1 cell).
- Experimental result (mechanism demonstrated by E5/E5b, but the symptom character is shape-dependent, correcting the v0 description):
  the sign and scale of mass drift depend on shape and are not a discriminating metric on their own. The decisive one is
  **the steady non-physical velocity in a static system**: in a stationary box + right Outflow + a pocket (initial rho=2.0),
  after 2000 steps, the edge cell of the buggy path keeps a huge steady velocity of **ux = −0.115**
  (control = at 1 cell inside the plug, ux = +0.00000 with complete rest, as physics dictates).
  Since all cells remain finite, it cannot be caught by NaN monitoring — it is "silent" as v0 claimed.
- Improvement: Following the freeze policy, the minimal fix = add to `set_solid` (V1 `sim.rs:792` / compat / wasm)
  an assert that rejects "placing solid at the immediate inner neighbor of an open-end edge cell."
  The BC fallback implementation is a V2 issue (recorded in B-8's design note).
- Acceptance: panic on the E5b shape (GUI rejects the paint). Legal shapes are probe_state_hash bit-invariant.
  Turn E5b into a regression test (confirm rejection).
- Effort: M / Experiment: **E5/E5b = mechanism confirmed, description corrected**

#### A-4 [S1] Establish the V2 native configuration validation layer `GlobalSpec::validate`
- Target: `crates/lbm-core2/src/solver.rs:286-357`, `crates/lbm-core2/src/params.rs:24-35`,
  duplicated validation on the consumer side `crates/lbm-scenario/src/lib.rs:582-654`
- Current state: `Solver::build` has only dimension and array-length asserts. (1) On an **uncovered face** (a face that is neither periodic, nor open,
  nor a wall rim), stale values are mixed in as real data every step and become non-physical without symptoms
  (reachable just by using D3Q19 with `GlobalSpec::default()`). (2) ν=0 makes omega_m=0 in
  `omegas()` and proceeds with non-physical relaxation. (3) Double application of periodic×open coaxial.
  (4) MAX_SPEED / rho>0 / u_conv range / cross-axis open faces (V1 rejects with
  `AdjacentOpenEdges`, V2 passes through = the Zou–He premise breaks on 3D edges) are also unchecked.
  Equivalent validation is **duplicately implemented** in compat and lbm-scenario.
- Improvement: Establish `GlobalSpec::validate() -> Result<(), SpecError>` in core2 and enforce it at
  the head of `Solver::build`. Check items: ν>0 / coverage of all non-periodic faces (open BC or all-face solid rim) /
  periodic×open exclusivity / cross-axis open rejection / MAX_SPEED (NaN-safe) / rho_bc>0 / u_conv∈(0,1] /
  force[2]==0 in 2D / open-face axis extent≥3. MAX_SPEED check also in `set_inlet_profile`.
  Replace scenario's hand-written checks with a validate call + error conversion.
- Acceptance: A unit test where all of the above invalid configurations become Err. Existing T13/T14/T15/v1_match green without modification.
  scenario's 3D validation test green without modification, and the duplicated check code disappears.
- Effort: M / Experiment: **E2 = confirmed** (D3Q19, z-face uncovered, z-uniform initial condition, 100 steps:
  mass drift 2.7e-3, z-invariance breakage 1.9e-4, spurious uz 2.6e-3 while nonfinite=0.
  The covered control is 0.0 on all metrics — quantitatively demonstrating "quietly non-physical without emitting NaN"),
  **E3 = confirmed** (`omegas(nu=0)` → TRT (ω₊,ω₋)=(2,0), BGK (2,2). Solver has
  no error in either construction or 10 steps)

#### A-5 [S1] Build-time rejection of halo-exchange scope misuse
- Target: `crates/lbm-core2/src/halo.rs:47-62,246-268`, `crates/lbm-core2/src/solver.rs:262-284`
- Current state: `Subdomain::neighbors` is a global part id but `exchange_f_generic` resolves it as a local
  index. Misuse of `new_local_part` + `LocalPeriodic`/`InProcess` is only a doc note, and
  when the neighbor id is 0 it becomes **silently wrong physics** as a self-wrap (OOB panic if id≥1).
- Improvement: Add `const SCOPE: ExchangeScope { Local, Remote }` to `HaloExchange`, and
  require `SCOPE == Remote` in `Solver::build(only=Some(part))` (mismatch is a build-time panic).
- Acceptance: A regression test where the misuse configuration becomes a build-time error. T13 / T13-MPI green without modification.
- Effort: S / Experiment: **E4 = confirmed** (part=1 of [2,1,1] periodic x + LocalPeriodic ran without panicking,
  and compared with the correct 2-part InProcess run, the rho of the owned block deviated by up to 7.7e-2
  — demonstrating silently wrong physics)

#### A-6 [S2] Reject the normal component of MovingWall
- Target: `crates/lbm-core/src/domain.rs:284-327`, `crates/lbm-core/src/sim.rs:1421-1425`
- Current state: The momentum injection of half-way BB is consistent only for the tangential wall velocity. The normal component does not diverge and keeps
  silently injecting/draining mass (the sign depends on the direction).
- Improvement: In `validate()`, reject a MovingWall with an edge normal component with `InvalidParameter`
  (V1+compat simultaneously). Add the reason to the doc.
- Acceptance: One with a normal component returns Err. Existing cavity-system tests green, bit-invariant.
- Effort: S / Experiment: **E7 = confirmed** (32×32 closed box, 500 steps: tangential u=[0.05,0] is
  drift +1.1e-13 (exact conservation), normal u=[0,−0.05] is mass 900→395.5, **−56.1%**.
  No error — since it does not diverge either, you cannot notice)

#### A-7 [S2] Input validation of `init_with` (V1+compat)
- Target: `crates/lbm-core/src/sim.rs:830-910`
- Current state: NaN immediately at rho=0 (0×inf), no MAX_SPEED check on velocity (`set_inlet_profile` has
  a check, so it is asymmetric within the API). The GUI's Droplet initialization pours JSON values in without validation.
- Improvement: `assert!(r > 0 && r.is_finite())` + MAX_SPEED check (message with coordinates), a Panics section in
  the doc. Note in compat's `init_with` that "the closure must be pure
  (re-evaluated up to 5 times near neighbors)."
- Acceptance: Coordinate-tagged panic on an invalid closure. probe_state_hash bit-invariant.
- Effort: S

#### A-8 [S2] D3Q19-only guard of `zou_he_face_3d` and ConvectiveOutflow contract test
- Target: `crates/lbm-core2/src/kernels.rs:494-508` (unknown 5 hardcoded),
  `kernels.rs:589-594` + `fields.rs:114-118` (stale-slot implicit contract)
- Current state: (1) The unknown set is fixed at 5. When D3Q27 is added (Q_MAX=27, which is in plan), it passes both
  compilation and execution while leaving 4 in place, wrongly making it undetectable. (2) ConvectiveOutflow's
  memory term depends on the implicit contract that "streaming does not write the unknown slot," spanning 4 modules
  (GPU re-implements it independently with an edge stash). It will surely break with in-place streaming (an M-E candidate).
- Improvement: (1) A `assert_eq!(L::unknowns(face).len(), 5)` guard at the head of the function (in the mid-term, resolve it together with the panic path of
  `dir_index` by const-ing the face-direction table). (2) Make it explicit in the `Backend::stream` contract that
  "the open-face unknown slot is invariant," and add a contract test that inspects the bit match with the pre-stream value in
  both CPU/GPU.
- Acceptance: A unit test of the guard assert. The contract test is green on both backends.
- Effort: S+S

#### A-9 [S2] V2 runtime NaN watchdog
- Target: `crates/lbm-core2/src/solver.rs:877-891` (`local_nonfinite_count` is manual-call only)
- Current state: Divergence detection is scattered across 4 ways in CLI 2D / CLI 3D / MPI / GPU, and the GPU path has no means.
  The kernel is unguarded to maintain V1 equivalence (this is correct).
- Improvement: Standardize `Solver::run_guarded(steps, check_every) -> Result<(), Diverged{step}>` using the finite check of
  the existing f64 aggregate `local_mass_partials` (NaN propagates to the total sum).
  The CLI/MPI drivers just call this. GPU is provisionally handled via the readback path of the same API.
- Acceptance: Injecting NaN into 1 cell → detected with a step number within N steps. Overhead <1% (512²).
- Effort: S–M

#### A-10 [S3] Small-grain consistency bundle
- (a) Remove the unused_mut of `crates/lbm-core/src/multiphase.rs:310` (resolves diff noise with compat,
  binary-invariant). (b) Fix the misleading comment of `sim.rs:801-806` (multiphase does not read solid rho).
- (c) Resolve the ±25%/±15% notation inconsistency of `t15_3d.rs:455-470` (delete the old numbers on the VALIDATION.md side).
- (d) Update the "bit-for-bit" claim in the doc at the head of kernels.rs to reality (starts bit-matching pre-fusion V1,
  the current one is a ≤1e-11 constraint, measured ~1.6e-14/50steps).
  (e) Add a solid/periodic-consistency debug_assert to MCMP `update_forces`. (f) A bit-match
  property test of `equilibrium()` and feq inside collide.
- Effort: S (bulk)

### WP-B: Structural Modifications (premise preparation to be completed before starting M-E)

#### B-1 [S1] Backend `Fields` generalization and GpuSolver integration (M-E most important)
- Target: `crates/lbm-core2/src/solver.rs:218,242-243`, `dist.rs:280-284`,
  `gpu/solver.rs:40-50,142-169`, `gpu/backend.rs:810-814`, `halo.rs:33-40`
- Current state: `Solver`/`MpiSolver` are fixed to `Fields = SoaFields<T>`, and `WgpuBackend`
  (`Fields = GpuFields`) does not mount. GpuSolver duplicates the step sequence (with known deviation).
  The GPU-side `stream` asserts `CellRange::full` and rejects the two_pass split.
  Multi-GPU / MPI+GPU cannot be composed at the type level.
- Improvement (staged commissioning):
  1. Formalize `stage_in/stage_out` (host⇔device transcription; `WgpuBackend::upload` is the existing substance) in
     `Backend`, and have `Solver` hold `SoaFields` as host staging, transcribing only at edit boundaries
     (generalization of GpuSolver's `host_dirty`/`device_ahead` mechanism).
  2. Unify gather/diagnostics via `read_moments`/`reduce`, establish `Solver<D2Q9, f32, WgpuBackend,
     LocalPeriodic>`, and delete GpuSolver's own step sequence.
  3. Add band dispatch (y-range specification via uniform) to the fused kernel and withdraw the
     `stream(range)` assert (the premise for overlap and future multi-GPU).
  4. Make `HaloExchange` `Backend::Fields`-generic at the pack/unpack buffer boundary
     (the GPU edge stash and dist.rs's face protocol are the templates).
- Acceptance: T14 is green via the same orchestrator. GpuSolver deleted. The MLUPS of bench_gpu
  regresses ≤3% (same measurement procedure).
- Effort: L / Experiment: **E10** (measurement of the current submit granularity = the performance baseline at modification time)

#### B-2 [S1] Organization of the Backend synchronization-point contract (probe / moments / end_step)
- Target: `crates/lbm-core2/src/backend.rs:94-110`, `gpu/backend.rs:803-831,855-868`
- Current state: Against the contract that `stream` returns the probe force synchronously, GPU returns zero (a trap where, once mounted on Solver,
  `probed_force()` silently returns 0), and `update_moments` is meaning-repurposed as a submit hook.
- Improvement: Remove the probe force from the return value of `stream` and formalize it as `read_probed_force` (explicit readback).
  Add an `end_step` hook to the trait and separate submit. Make `update_moments` explicit as a lazy contract.
  Declare two-pass non-support with a capability method (a transitional measure until it is resolved by B-1's item 3).
  Make V2's probed_force a fixed-order fold of band partial sums, deterministic in bits across runs
  (resolving rayon-reduce non-determinism, enabling incorporation into the state hash).
- Acceptance: No residual zero-return or meaning-repurposing. Adding a T14 case with a probe, CPU/GPU match on the same API.
  probed_force bit-matches across 2 runs with the same thread count.
- Effort: M

#### B-3 [S1] Unification of the Shan–Chen implementation and V2 native multiphase
- Target: `crates/lbm-core2/src/solver.rs:679-748`, `compat/multiphase.rs:134-374`,
  V1 `multiphase.rs` (scenario imports it in production)
- Current state: The SC force stencil is 3 systems, 5 loops (PM confirmed: V1:195,351 / compat:196,352 /
  solver.rs:736). Wall adhesion, virtual wall density, and two-component are compat/V1 form only (= 2D/CPU-limited),
  MPI/3D is neutral single-phase only, GPU has no multiphase.
- Improvement: Absorb the wall term (g_wall, wall_rho, same accumulation order as compat version 347-365) into
  `Solver::update_shan_chen_force` → replace compat `ShanChen` with a thin delegation → port `MultiComponent`
  to V2 native as "2 Solvers + exchange_scalar." GPU multiphase is M-F consideration
  (out of scope in this spec, recorded in the B-9 note).
- Acceptance: The SC stencil loop is in 1 place within lbm-core2. validation_contact_angle /
  multiphase / rt green without modification. Add 1 T13 extension for contact angle with wall_rho on the MPI path.
- Effort: M–L

#### B-4 [S2] Compat migration of the 2D production path and true freezing of V1
- Target: `crates/lbm-scenario/src/lib.rs:7-8`, `crates/lbm-cli/src/runner.rs:6-7`,
  `crates/lbm-wasm/Cargo.toml:12`
- Current state: There is 0 production code that imports compat. 2D scenario, CLI, and wasm/GUI run on V1, and
  V1 is load-bearing while claiming to be "frozen." The V1 review confirmed that all of the APIs used by wasm are
  covered by compat, zero-copy `*_ptr` compatibility, and buildability with `default-features = false`.
- Improvement: Replace scenario / CLI / wasm's `lbm_core::` with `lbm_core2::compat::`. Demote V1 to a
  dev-dependency (equivalence-test only). Judge the complete deletion of the `v1-retirement` branch
  (33e130a) separately after this migration accumulates a track record.
- Acceptance: The whole workspace green, wasm-pack build succeeds, 0 `lbm_core::` references in production crates,
  the field hash of GUI presets matches the V1 version (add 1 wasm smoke = D-11).
- Effort: M / Dependency: A-1 (the compat suite being genuine)

#### B-5 [S2] Restart/state-injection API (snapshot)
- Target: `crates/lbm-core2/src/solver.rs:948-965` (gather only, no load)
- Improvement: Implement `Solver::snapshot() -> StateV1 { f[Q], solid, wall_u, force_field, time }` /
  `restore()` as a pair (internally halo fill + update_moments). MpiSolver uses the same API via rank0
  gather/scatter (the foundation for C-8's distributed checkpoint).
- Acceptance: "N steps → snapshot → restore → M steps" = "N+M steps continuously" is an f64 bit match.
- Effort: M

#### B-6 [S2] Groundwork for per-cell relaxation rate (LES premise, avoiding the three-way modification right before M-F)
- Target: `crates/lbm-core2/src/params.rs:76-105`, `kernels.rs:91-173`
- Improvement: Add `omega_field: Option<Vec<T>>` to `SoaFields`, and have `collide_row` use per-cell omega only when Some.
  The None path is bit-identical (assured by probe_state_hash). GPU is flag-controlled with 1
  storage buffer (the premise is GPU-8's limit raise). The LES body is M-F.
- Acceptance: Bit-identical to all existing tests at None. A test where the uniform value = the scalar specification match.
- Effort: M

#### B-7 [S2] Sealing of the public backdoor and f64 unification of diagnostics
- (a) Make `fields_mut` pub(crate), and route necessary operations to dedicated methods with automatic dirty management
  (the source of the worst failure mode where a single-rank edit mistake becomes a silent hang under MPI).
- (b) Add `total_mass_f64()` to the facade (resolving diagnostic quantization ~0.06/10⁶ cells at f32).
- (c) A dirty-flag consistency Allreduce (1 byte) only in debug at the head of MpiSolver::step for
  fail-fast (zero cost in release).
- Acceptance: Automation of dirty on the mask-edit path, the single-rank-edit debug test fails an assert rather than hanging.
- Effort: S–M

#### B-8 [S2] Kernel extension-point design note (no implementation, 1 docs sheet)
- Per-cell omega passing convention (B-6), placement of MRT/cumulant kernels (the CollisionKind branch and
  the location of transformation matrices), the per-link wall-distance sparse structure of curved boundaries (Bouzidi), the
  ConvectiveOutflow alternative at in-place streaming (generalization of the GPU edge-stash scheme), the
  BC fallback of Outflow × solid adjacency (the permanent solution of A-3). Fix "existing-configuration bit invariance" in each extension's DoD.
- Effort: M (until review approval)

### WP-C: Scale / Operations (MPI / GPU)

#### C-1 [S1] Localization of MPI setup (resolution of global-array replication)
- Target: `crates/lbm-core2/src/dist.rs:305-312`, `solver.rs:298-333,76-125,628-654`
- Current state: `MpiSolver::new` requires global compact arrays (solid/wall_u) on all ranks,
  `build_wall_rims` also generates the whole domain, and `set_solid` is a global all-cells×all-ranks call.
  At 10⁹ grid, wall_u ≈ 24 GB/rank → guaranteed OOM in the weak-scaling configuration.
  **A structural blocker that manifests with grid size, not rank count** (the current test is undetected because n≤8, small grid).
- Improvement: A closure-taking `MpiSolver::new_with(solid: impl Fn(x,y,z)->bool, …)` (local evaluation of the same
  form as init_with) + a local version of `build_wall_rims` + batch `set_solids_where(pred)`.
  The existing API is left in place for small scale.
- Acceptance: Peak RSS per rank is O(N/P)+constant (measured). Field bit match via the T13-MPI new API.
- Effort: M

#### C-2 [S2] Exchange overlap (post/finish split and two_pass connection)
- Target: `dist.rs:171-189`, `solver.rs:372-429`
- Improvement: Add `post_f()/finish_f()` to `HaloExchange` (InProcess completes immediately).
  Rewire the step to collide → post → interior stream (the interior range of the existing two_pass) → finish →
  boundary shell. The first stage can overlap the x phase only.
  **Premise fix**: two_pass's boundary_shells overlap at 1-width axes and the probe is double-counted
  (`solver.rs:972-1007`; the field is idempotent and harmless, invisible in T13). Fix the shells to be mutually disjoint
  before connecting.
- Acceptance: probed_force match of two_pass on/off at a 1-width axis (→ experiment E8 demonstrates the current double-count).
  T13-MPI all PASS maintained. Measurable reduction of exchange-wait occupancy.
- Effort: L / Experiment: **E8 = confirmed** ([64,1,1], obstacle 1 cell, probe, 20 steps:
  the on/off probed_force ratio = **2.000** (exact double-count), total_mass matches for both
  — the field is intact and only the probe breaks = undetectable in T13's field comparison, as claimed)

#### C-3 [S2] Parallel I/O (per-rank raw + manifest, eliminating the rank0 whole-domain buffer)
- Target: `dist.rs:504-574`
- Improvement: Short-term = each rank writes its own block to an individual file, rank0 writes only the manifest.
  Mid-term = MPI-IO subarray (rsmpi's File support needs confirmation). gather_* is left, marked as validation-use.
- Acceptance: At output, rank0 peak RSS is O(N/P). Output time is non-increasing with rank count.
- Effort: M–L

#### C-4 [S2] Deferral of probe Allreduce
- Target: `dist.rs:369-379` (a 3-double Allreduce every step when probe is enabled)
- Improvement: Make `probed_force()` collective and cache it with the `time` key (reduce only when queried).
- Acceptance: In a probe-enabled benchmark, the per-step collective disappears from the profile. mpi_t13 PASS maintained.
- Effort: S

#### C-5 [S2] Persistence of exchange buffers (a stepping stone to GPU-aware MPI)
- Target: `dist.rs:125-126,175-181,199-229`, `solver.rs:691-711` (ψ plane allocated every step)
- Improvement: Have `MpiExchange` hold send/receive buffers per face×kind and reuse them. Make ψ plane / staging
  fields too. Concentrate buffer ownership in MpiExchange (make the swap point to GPUDirect a single place).
- Acceptance: Zero heap allocation during steady steps (measured). bench_mpi non-regressing. T13-MPI bit match maintained.
- Effort: S–M

#### C-6 [S2] Inter-rank spec consistency check
- Target: `dist.rs:305-333`
- Current state: A mismatch in nu, faces, or mask content is **not detected because the message length matches**, and
  emits a "plausible" field discontinuous at seams (a job-script-accident class).
- Improvement: At build time, compare the spec-normalized byte sequence + mask FNV hash by Allreduce(min/max) and
  abort with the item name on mismatch.
- Acceptance: A 2-rank injection test with only nu changed immediately gives an explicit error. Normal-case cost immeasurable.
- Effort: S

#### C-7 [S2] Funneled-ization of the MPI thread level
- Target: `examples/mpi_t13.rs:391`, `examples/bench_mpi.rs:27`, `backend.rs:36,133-147`
- Current state: rsmpi 0.8.1's `initialize()` = `Threading::Single` (confirmed by PM from a registry source).
  On the other hand, the default feature `parallel` launches rayon above 16,384 cells → at real sizes,
  multithreaded execution under an MPI_THREAD_SINGLE declaration (an MPI-spec violation. Can corrupt/hang on UCX/OFI systems.
  The current test is coincidentally serial at ≤6,144 cells/rank).
- Improvement: Add `dist::init_mpi() -> Universe` (requests Funneled, explicit error if provided is insufficient AND parallel is enabled)
  and migrate the examples/guide.
- Acceptance: With 2 ranks × forced rayon (parallel_min_cells lowered), the T13-MPI equivalent PASSes,
  and provided ≥ Funneled is confirmed in the log.
- Effort: S / Experiment: **E9 = confirmed (source level, see §3)**

#### C-8 [S2] Distributed checkpoint/restart
- On top of B-5, a collective `MpiSolver::save(dir)/load(world, dir, backend)` (per-rank raw +
  rank0 manifest, spec-hash and decomp consistency validation). Bit-match on resume via raw storage of the deviation-storage f.
- Acceptance: "50 steps → save → load → 50 steps" = "100 steps continuously" is a field bit match. A manifest mismatch is an explicit error.
- Effort: M / Dependency: B-5, C-6

#### C-9 [S1] Time calibration of the GPU submit chunk and Result-ification of device-lost
- Target: `gpu/backend.rs:187-188,226` (`submit_chunk: 200` fixed), `:94-98,301-314` (expect panic)
- Current state: 200 steps (up to ~1000 dispatch) in a single submit. Mechanically time-unbounded →
  on a slow GPU, exceeding Windows TDR (default 2 s) → device removed → **process panic**.
  No recovery path.
- Improvement: Auto-calibrate to a target of 1 submit ≈ 100–250 ms from the first-chunk measurement (keeping the upper limit 200).
  Result-ify `wait_idle`/`map_staging` to `Result<_, GpuError>` and propagate device lost.
  Capture the reason with `set_device_lost_callback`.
- Acceptance: The MLUPS of bench_gpu regresses ≤3%. A unit test of the calibration logic. A test where a poll failure returns Err.
- Effort: M / Experiment: **E10 = confirmed** (measured on this machine M5 Max/Metal: 2048² TGV 5,719 MLUPS →
  200-step chunk = **147 ms/submit**. Since this value is on a top-class consumer GPU,
  the extrapolation holds that on a ~15× slower GPU (a several-hundred-MLUPS-class iGPU) the same grid exceeds TDR 2 s.
  For reference: 512²=11,509 / 1024²=6,607 MLUPS, −5.3 to −18.0% vs proto (within the ±20% acceptance line))

#### C-10 [S2] Pre-validation of GPU resource limits
- Target: `gpu/backend.rs:650-692,66-91,236-241`
- Improvement: At the head of `alloc`, check required bytes vs `device.limits()`, `Q*n ≤ u32::MAX` (D3Q19 overflows at
  226 million cells), and dispatch count ≤65,535, and Err with a reason.
- Acceptance: `GpuSolver::new` gives an explicit error on an over-limit grid. T14 unchanged green.
- Effort: S

#### C-11 [S2] Efficiency improvement of the GPU diagnostic path
- Target: `gpu/backend.rs:318-343,870-911`, `gpu/solver.rs:184-207`
- Improvement: Arc-ify `f_cache` (eliminate cloning), make FluidCells readback unnecessary (complete with host_solid),
  consolidate the 3 readbacks of sync into 1 encoder/1 wait. Add the GPU-side 2-stage reduction as a fast mode of
  M-E (keep the host f64 path for T14).
- Acceptance: The sync+diagnostics triple is ≥3× faster at 2048². T14 diagnostic values bit-identical.
- Effort: S (+M)

#### C-12 [S2] FP16 plumbing (a premise of the M-E body)
- Target: `gpu/backend.rs:80` (Features::empty fixed), `gpu/wgsl.rs:216`, element size `*4` scattered
- Improvement: Conditional request of `SHADER_F16`, `generate::<L>(cfg: KernelCfg { storage: F32|F16 })`-ization
  (confine the change points to the 2 places "buffer declaration + load/store wrapper," keeping the arithmetic f32),
  hold element size in `GpuFields`. Non-supporting adapters give an explicit error (do not silent-fallback).
- Acceptance: Establish T16 (quantify the degradation of f16 storage in a frozen band). MLUPS ≥1.5× vs the f32 version at 2048² TGV.
- Effort: M / Dependency: B-1 (least rework after orchestrator integration)

#### C-13 [S2] Contract-layer connection of GPU (scenario `backend: "gpu" | "auto"`)
- Target: `crates/lbm-scenario/src/lib.rs:575-577`, `crates/lbm-cli/src/main.rs:190`
- Improvement: Unlock in the feature-gpu build, implement a capability check (reject f64/3D/unsupported BC with a reason) in
  validate, and `auto` selects by a measured threshold (e.g. n≥256²) and logs the result explicitly.
- Acceptance: The cavity/cylinder presets complete with `backend:"gpu"` and the field difference from CPU is within T14 tolerance.
  `f64`+`gpu` is rejected with a reason.
- Effort: M / Dependency: B-1
- The public benchmark (M-E) has publishing of the reproduction procedure as its main body — this item is the premise.

#### C-14 [S2] Reduction of the GPU memory footprint
- Target: `gpu/backend.rs:677-678,691,133`
- Improvement: Dummy-buffer-ize force_field/wall_u when absent, lazy-allocate and right-size staging.
  The current state is +284 MB (about 1.9×) at 2048² TGV — 3D-ization directly hits the maximum grid of an 8–16 GB-class GPU.
- Acceptance: Without a force field or solid, the allocation is ≤ 2×f + moments + mask + O(perimeter). T14 green.
- Effort: S

#### C-15 [S3] GPU small-grain bundle
- (a) Raise limits such as `max_storage_buffers_per_shader_stage` (currently the step kernel is
  exactly at the default upper limit of 8 — a landmine that panics at runtime with 1 added). (b) A parse+validate unit test of the generated WGSL by naga
  (no GPU required, with a golden file). (c) Single-table generation or a collation test of the Rust index ⇔
  WGSL field of BcParams. (d) Result-ify `GpuContext::new` to
  `Result<_, GpuInitError>` (with adapter_info; makes the auto fallback diagnosable).
  (e) Note in one line in GPU_EVALUATION.md the non-determinism of the probe-force f32 CAS addition.
- Effort: S×5

#### C-16 [S3] MPI small-grain bundle
- (a) `choose_decomp` (auto split with minimum surface area) + arbitrary-n support of mpi_t13 (add cases for n=3,5,6 and
  non-divisible dims). (b) 3D/D3Q19 and strong/weak-mode-ization of bench_mpi
  (the current 2D-strip-fixed cannot be used for R3 main measurement). (c) A deterministic total sum `total_mass_deterministic()`
  (fixed-order composition of global row-wise partial sums), added as an option.
- Effort: S+S+M

### WP-D: Validation / Process

#### D-1 [S1] Three-stage setup of CI
- Currently `.github/` is absent (PM confirmed). All quality claims are manual snapshots.
- Improvement: (1) `cargo test --workspace --release` on push/PR, (2) nightly `--include-ignored`,
  (3) GPU/MPI on a self-hosted runner (this machine M5 Max + Open MPI in ~/.local) with tagged execution.
  While remote is undecided, start with a pre-merge local hook + appending execution logs to `docs/CI_LOG.md`.
- Acceptance: A default suite green is mechanically enforced for merge. Weekly full+gpu+mpi record.
- Effort: M

#### D-2 [S1] Cargo-testing of MPI logic
- Target: `dist.rs` (0 `#[test]` — PM confirmed)
- Improvement: Unit-test pack/unpack and phase plan as pure functions (same buffer content as InProcess,
  tolerance ==0.0). Turn `test_mpi.sh` into a nightly job. Bring the 1-rank self-exchange smoke into the cargo sphere.
- Acceptance: The main logic of dist.rs is tested without mpirun. T13-MPI weekly PASS log.
- Effort: M

#### D-3 [S1] Absolute physics validation of GPU and skip mechanism
- Target: `tests/t14_backend_equiv.rs` (CPU-relative equivalence, only 8, expect panic without an adapter)
- Current state: A bug where CPU and GPU **break in the same direction** (a shared spec-interpretation mistake) is undetectable by relative equivalence.
- Improvement: 2 GPU-direct absolute tests (TGV convergence order ≥1.7, cavity Ghia RMS ≤0.02U — calibrated and frozen with
  f32 measurement). Adapter absence is skip (`LBM_REQUIRE_GPU=1` promotes to fail).
  Make non-support of 3D GPU / GPU multiphase explicit in VALIDATION.md's known-limitations section.
- Acceptance: The 2 absolute tests are green with `--features gpu`. Skip exit on a GPU-less host.
- Effort: M

#### D-4 [S1] Addition of f32×3D validation
- Target: `crates/lbm-scenario/src/lib.rs:490-492` (`Sim3Handle::F32` = product path),
  `tests/t15_3d.rs` (0 f32 occurrences — PM confirmed)
- Improvement: f32 versions of t15-1 (z-invariant 2D degeneration, f32 relative ≤1e-5) and t15-4 (TGV3D decay rate ±2%) +
  mass drift ≤1e-5/10³step. Freeze the measured values as an f32 row in VALIDATION.md T15.
- Effort: S

#### D-5 [S1] Horizon extension of equivalence validation and native absolute validation
- Current state: V2's absolute physics validation is only the transitivity of "≤1e-11 with V1 at 500 steps → V1 is validated"
  (Ghia steady state is ~99k steps, 200× the horizon). The V2 native API (3D/scenario path) is outside the chain.
- Improvement: After A-1 completes, (1) 1 long-time equivalence (cavity Re=100 for 20k steps, ≤1e-9, may be #[ignore]),
  (2) 1 TGV-convergence-order test of the native `Solver` directly (also makes the deviation between facade and core detectable).
- Effort: S / Dependency: A-1

#### D-6 [S2] Consistency recovery of the acceptance-criteria source
- Target: `docs/COMPETITIVE_SPEC.md:57-59` (R1: sphere ±5% etc. / R3: weak scaling ≥85%) vs the implementation
  (sphere ±10%, D_h normalization / R3 is 73.2% at n=8 → shrunk to "n≤4 local" in a commit)
- Improvement: Update R1/R3 with revision history (a link to the basis of ±10% and the D_h definition, and an explicit note of the
  n≤4 local line and the cluster-condition non-attainment). Correct PLAN's "R1 achieved" notation to "(excluding 3D cavity, ±10% revised version)."
  Give T15.5 (3D cavity Ghia table) a backlog position.
- Acceptance: Every "achieved" declaration corresponds 1:1 to the current version of the source.
- Effort: S

#### D-7 [S2] Addition of T13/T14 sections to VALIDATION.md
- Centralize T14's "6 configurations 1e-5 / pressure 1e-4 + 1-ulp control" and T13's "field ==0.0 / diagnostics 1e-12" into
  the spec (in a form that can be commissioned to codex). Make the test headers thin, referencing the spec.
- Effort: S

#### D-8 [S2] codex adversarial-test order #7 (following up on T14/T15)
- Current state: Adversarial commissioning goes up to T13 (order #6). GPU equivalence and 3D physics are implementation-side self-made only.
- Improvement: Commission from D-7's spec revision: T14 attacks (initial discontinuity right above a boundary face, a probe touching a face,
  u→near MAX_SPEED), T15 attacks (a perturbation that breaks z degeneration, an extreme aspect ratio, an off-center sphere).
- Acceptance: 1 round + a triage record remains in TESTING_NOTES.
- Effort: M (including commissioning and triage)

#### D-9 [S2] Performance regression detection
- Improvement: A core2 version of bench_mlups (resolving that V2 CPU's performance claim is still a citation of V1 Phase 9 numbers) +
  a `--check` mode (fail on deviation ±25% from the frozen value) + nightly execution appending a JSON history to `docs/bench_history/`.
  probe_state_hash is also enforced nightly.
- Acceptance: A 25% regression manifests as a fail within 24h.
- Effort: M / Experiment: **E10** provides the GPU-side baseline value

#### D-10 [S2] `lbm verify` and LIMITATIONS.md
- Improvement: A `lbm verify` that runs a validation subset (a-few-minutes class) and outputs a comparison table against the acceptance bands.
  Centralize known limitations (GPU=2D/f32, no CP/restart, long-time validation limit, τ→0.5 guideline, etc.) into
  `docs/LIMITATIONS.md` and ship it paired with releases.
- Effort: L (verify) + S (LIMITATIONS)

#### D-11 [S2] wasm smoke test
- TGV 100-step mass conservation + native f64 match ≤1e-12 with `wasm-pack test --node`.
  The web side is minimal, only a schema round-trip of the written-out JSON.
- Effort: M / Dependency: efficient to carry out simultaneously with B-4

#### D-12 [S2] Explicitation of the CpuScalar performance gap (handover to M-E)
- Current state: V2 CPU is phase-separated 3-pass + per-cell branching, and the single-core gap vs the V1 fused kernel is
  unmeasured (the denominator of the 97-99% weak scaling may be lax).
- Improvement: A q-major port of V1 step_band in M-E's CpuSimd (planned). Until then, in the public benchmark
  co-note the V1 single-core ratio. Acceptance: CpuSimd bit-matches CpuScalar and single-core ≥ 0.9x of V1.
- Effort: L (the M-E body)

---

## 3. Results of Validation Experiments (carried out 2026-07-05, this machine M5 Max. Re-run: `scripts/spec-experiments/`)

**All 10 experiments carried out. Zero refutations, 2 descriptive corrections** (E5→the symptom character of A-3,
E7→the numbers and sign of A-6). The raw output of the experiment logs is transcribed into the table below as each experiment's RESULT line.

| ID | Validation target | Measured result | Judgment |
|---|---|---|---|
| E1 | A-1/B-4/D-5 | perl port + guard applied → 16 files regenerated → `cargo test -p lbm-core2 --release` **whole green (exit 0)** | **confirmed, carried out**. First demonstration of R5. Resolves B-4's blocker |
| E2 | A-4 | Uncovered z-face: mass drift 2.7e-3, z-invariance breakage 1.9e-4, spurious uz 2.6e-3 while nonfinite=0 (the covered control is 0.0 on all metrics) | **confirmed** (quantitatively demonstrates silent non-physics) |
| E3 | A-4 | `omegas(0)` = TRT (2, 0), BGK (2, 2). No error in Solver construction or 10 steps | **confirmed** |
| E4 | A-5 | part=1+[2,1,1]+LocalPeriodic: no panic, rho max deviation vs the correct answer **7.7e-2** | **confirmed** (silent wrong physics) |
| E5/E5b | A-3 | Mass drift is shape-dependent and unfit for judgment (corrected the v0 description). Decisive: the edge cell of a pocket in a static box retains a steady velocity of **ux=−0.115** even after 2000 steps (the control is +0.00000) | **mechanism confirmed, description corrected** |
| E6 | A-2 | NaN inlet: build()=Ok→42 non-finite cells in 3 steps. NaN MovingWall: build()=Ok→**bit-matches** a static wall | **confirmed** |
| E7 | A-6 | Tangential: drift +1.1e-13 (exact conservation). Normal u=[0,−0.05]: mass 900→395.5 (**−56.1%**, no error). v0's "+115%" turned out to be a sign difference dependent on direction | **confirmed, numbers corrected** |
| E8 | C-2 | [64,1,1] 1-width axis: the two_pass on/off probed_force ratio = **2.000**, total_mass matches (the field is intact) | **confirmed** |
| E9 | C-7 | Directly confirmed in rsmpi 0.8.1 `environment.rs:268-270` that `initialize()` = `initialize_with_threading(Threading::Single)`. Against PARALLEL_MIN_CELLS=16,384, rayon launch at a real-operation grid (512²/rank=262,144 cells) is arithmetically certain. Runtime demo omitted (symptoms are hard to reproduce on the shared-memory BTL and are unnecessary for the judgment) | **confirmed (source level)** |
| E10 | C-9/B-1/D-9 | 2048² TGV **5,719 MLUPS** → 200-step chunk = **147 ms/submit** (512²=11,509 / 1024²=6,607. −5.3 to −18.0% vs proto, within the ±20% line) | **confirmed** (TDR extrapolation holds, the basis of the calibration initial value) |

Auxiliary (already confirmed by code inspection, no experiment needed): the `Fields=SoaFields` binding (B-1), GPU probe zero return
(B-2), SC 5 loops (B-3), `MpiSolver::new`'s global-array signature and the 24 GB/rank arithmetic
(C-1), `submit_chunk: 200` (C-9), 0 `#[test]` in dist.rs (D-2), absence of `.github` (D-1),
0 f32 occurrences in t15 (D-4), rsmpi `initialize()`→`Threading::Single` (C-7, confirmed from a registry source).

---

## 4. Order of Execution (confirmed)

All experiment dependencies are resolved. Each item is described at a granularity that can be commissioned independently.

- **Carried out (this branch)**: The body of A-1 (script fix + guard + regeneration + all-green confirmation).
  Remaining: 1 static guard test (include in R-Phase 1).
- **R-Phase 1 (immediate, ~2 days)**: A-2–A-10 + A-1 remaining work + D-6/D-7 (document consistency).
  All S–M. Passing existing tests green without modification and legal-configuration bit invariance are the common DoD.
- **R-Phase 2 (M-E premise preparation, ~1.5 weeks)**: B-1 → B-2 → (in parallel) B-3, B-5–B-8 /
  C-9–C-11 / D-1–D-5. B-4 (compat migration) can be started anytime since it is demonstrated by E1
  (recommended to carry out simultaneously with D-11's wasm smoke).
- **R-Phase 3 (scale, can run in parallel with M-E)**: C-1, C-4–C-7, C-16 → (C-2, C-3, C-8) /
  C-12–C-15 after B-1 / D-8–D-10.
- The M-E body (FP16 implementation, multi-GPU, public benchmark, CpuSimd = D-12) takes this spec's
  B-1/B-2/C-9/C-12/C-13/D-9 as premises.

Organization: Implementation is commissioned to Opus/Sonnet subagents / codex per WP. Validation tests
(D-3/D-4/D-8 and each acceptance test) are created adversarially by codex from the VALIDATION.md revised version (D-7) and
separated from the implementation (maintaining the conventional protocol).
