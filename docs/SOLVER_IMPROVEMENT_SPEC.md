# Solver Improvement Specification — active residuals

Origin: full-review 2026-07-05 of the V2 core plus supporting infrastructure.
For the historical review record (E1–E10 experiments, A-1 sync-tests fix, WP-A
correctness/entry guards), see git history around commit 84abaa3 and its
predecessors.

**Audit note (2026-07-07):** This pass checked each numbered item against the
current worktree by grepping the named symbols/behaviors in `crates/`,
`scripts/`, `.github/`, and relevant tests/docs. Status markers below correct
the live/open inventory; item bodies are retained as the historical record.

**R-Phase 1 landed.** All WP-A items (A-1..A-10) plus D-6/D-7 doc-consistency
work were bundled into R-Phase 1 and are on main. V1 was retired in the same
window; `crates/lbm-core2` was renamed to `crates/lbm-core`, and legacy V1 lives
in `crates/lbm-core/src/compat/`. Below is the surviving spec: R-Phase 2 (M-E
premise) and R-Phase 3 (scale, MPI/GPU) items still not landed. Cross-refs to
HANDOFF-PM-2026-07-07 §4 name the currently open bundle.

Severity: S0 correctness / S1 high risk / S2 improvement / S3 minor.
Effort: S = hours, M = ~1 day, L = several days.

Every item's DoD includes: existing tests green without modification;
legal-configuration bit invariance (probe_state_hash unchanged where noted).

---

## WP-B — Structural (M-E premise prep)

### B-1 [S1] Backend `Fields` generalization and GpuSolver integration (M-E's most important)
**STATUS (2026-07-07 audit): PARTIALLY RESOLVED — `Solver` now owns `B::Fields` and can mount `WgpuBackend` (`crates/lbm-core/src/solver.rs:1250`, `crates/lbm-core/src/gpu/backend.rs:1885`, `crates/lbm-core/src/gpu/backend.rs:2459`); remaining: monolithic GPU only, no generic GPU halo/multi-part path, and `GpuSolver` still exists as a wrapper (`crates/lbm-core/src/gpu/backend.rs:1933`, `crates/lbm-core/src/gpu/solver.rs:17`).**
Target: `crates/lbm-core/src/{solver.rs, dist.rs, halo.rs, gpu/{solver,backend}.rs}`.
`Solver`/`MpiSolver` are fixed to `Fields = SoaFields<T>`, so `WgpuBackend`
(`Fields = GpuFields`) does not mount; GpuSolver duplicates the step sequence,
`stream` asserts `CellRange::full`, and multi-GPU / MPI+GPU cannot compose at
the type level.

Staged plan:
1. Formalize `stage_in/stage_out` (host⇔device) on `Backend`; keep `SoaFields`
   as host staging, transcribe at edit boundaries (generalize GpuSolver's
   `host_dirty`/`device_ahead`).
2. Unify gather/diagnostics via `read_moments`/`reduce`; establish
   `Solver<D2Q9, f32, WgpuBackend, LocalPeriodic>` and delete GpuSolver's step.
3. Add band dispatch (y-range via uniform) to the fused kernel; withdraw the
   `stream(range)` assert (premise for overlap and multi-GPU).
4. Make `HaloExchange` `Backend::Fields`-generic at the pack/unpack boundary.

Acceptance: T14 green on the unified orchestrator; GpuSolver deleted; bench_gpu
MLUPS regresses ≤3%. Effort: L. Handoff status: **HANDOFF §4 lists B-1 as owner-
scheduled into the M-E/M-F integration campaign; not currently dispatched.**

### B-2 [S1] Backend synchronization-point contract
**STATUS (2026-07-07 audit): PARTIALLY RESOLVED — explicit GPU probe readback exists and `GpuSolver::probed_force()` uses it (`crates/lbm-core/src/gpu/backend.rs:960`, `crates/lbm-core/src/gpu/solver.rs:162`); remaining: `Backend::stream` still returns probe force and the unified `WgpuBackend` slot is still zeroed (`crates/lbm-core/src/backend.rs:189`, `crates/lbm-core/src/gpu/backend.rs:2100`).**
Target: `crates/lbm-core/src/backend.rs`, `gpu/backend.rs`.
`stream` is contracted to return probe force synchronously; GPU returns zero
(silent trap once mounted on Solver). `update_moments` is meaning-repurposed as
a submit hook.

Fix: remove probe force from `stream`; formalize `read_probed_force` (explicit
readback). Add `end_step` hook and split submit. Make `update_moments` an
explicit lazy contract. Declare two-pass non-support via a capability method
(transitional until B-1 item 3 resolves it). Make V2's probed_force a
fixed-order fold of band partial sums (deterministic in bits, resolves rayon-
reduce non-determinism, enables incorporation into the state hash).

Acceptance: no zero-return or meaning-repurposing; adding a T14 case with a
probe, CPU/GPU match on the same API; probed_force bit-matches across 2 runs
with the same thread count. Effort: M.

### B-3 [S1] Unify Shan-Chen and V2 native multiphase
**STATUS (2026-07-07 audit): PARTIALLY RESOLVED — native single-component Shan-Chen now has wall adhesion / virtual wall density and MPI scalar exchange (`crates/lbm-core/src/solver.rs:2278`, `crates/lbm-core/src/dist.rs:744`); remaining: compat `ShanChen` and `MultiComponent` still carry separate 2D facade implementations (`crates/lbm-core/src/compat/multiphase.rs:307`, `crates/lbm-core/src/compat/multiphase.rs:138`).**
Target: `crates/lbm-core/src/solver.rs`, `compat/multiphase.rs`.
SC force stencil currently exists in 3 places (V1, compat, V2 native).
Wall adhesion, virtual wall density, and two-component are compat/V1 only (2D/
CPU-limited); MPI/3D is neutral single-phase only; GPU has no multiphase.

Fix: absorb the wall term (`g_wall`, `wall_rho`, accumulation order per compat
lines 347-365) into `Solver::update_shan_chen_force`. Replace compat `ShanChen`
with a thin delegation. Port `MultiComponent` to V2 native as "2 Solvers +
`exchange_scalar`". GPU multiphase is M-F (out of scope, tracked in B-8).

Acceptance: SC stencil loop lives in one place inside `lbm-core`;
validation_contact_angle / multiphase / rt green; add 1 T13 extension for
contact angle with wall_rho on the MPI path. Effort: M-L. **HANDOFF §4 lists
B-3 among the owner-scheduled R-Phase 2 residuals.**

### B-6 [S2] Per-cell relaxation-rate groundwork (LES premise)
**STATUS (2026-07-07 audit): PARTIALLY RESOLVED — `SoaFields::omega_field`, `Solver::set_omega_field`, and CPU collide consumption landed with tests (`crates/lbm-core/src/fields.rs:195`, `crates/lbm-core/src/solver.rs:2483`, `crates/lbm-core/src/kernels.rs:191`, `crates/lbm-core/tests/wale_les.rs:48`); remaining: no generic staged `omega_field` upload path for `WgpuBackend`, only WALE-specific GPU omega buffers (`crates/lbm-core/src/gpu/backend.rs:1458`).**
Target: `crates/lbm-core/src/params.rs`, `kernels.rs`.
Add `omega_field: Option<Vec<T>>` to `SoaFields`; `collide_row` uses per-cell
omega only when `Some`. The `None` path is bit-identical (probe_state_hash
guarded). GPU flag-controlled with 1 storage buffer (premise: GPU-8's limit
raise). LES body itself is M-F.

Acceptance: bit-identical to all existing tests at `None`; a uniform-value case
matches scalar spec. Effort: M.

### B-7 [S2] Seal the public backdoor; f64 diagnostics
**STATUS (2026-07-07 audit): RESOLVED — `fields_mut` is `pub(crate)`, `total_mass_f64()` exists on solver/compat/MPI, and the MPI dirty debug Allreduce is covered by `dirty-mismatch` (`crates/lbm-core/src/solver.rs:2462`, `crates/lbm-core/src/solver.rs:2930`, `crates/lbm-core/src/dist.rs:642`, `scripts/test_mpi.sh:63`).**
- (a) Make `fields_mut` pub(crate); route mutations through dedicated methods
  with automatic dirty management (single-rank edit mistake under MPI is the
  worst failure mode — a silent hang).
- (b) Add `total_mass_f64()` to the facade (resolves diagnostic quantization
  ~0.06/10⁶ cells at f32).
- (c) Debug-only dirty-flag consistency Allreduce (1 byte) at head of
  `MpiSolver::step` for fail-fast (zero cost in release).

Acceptance: dirty automated on mask-edit path; the single-rank-edit debug test
asserts rather than hangs. Effort: S-M.

### B-8 [S2] Kernel extension-point design note
**STATUS (2026-07-07 audit): RESOLVED — `docs/KERNEL_EXTENSION_POINTS.md` covers per-cell omega, MRT/cumulant placement, Bouzidi records, ConvectiveOutflow, and Outflow x solid fallback (`docs/KERNEL_EXTENSION_POINTS.md:43`, `docs/KERNEL_EXTENSION_POINTS.md:73`, `docs/KERNEL_EXTENSION_POINTS.md:124`, `docs/KERNEL_EXTENSION_POINTS.md:171`, `docs/KERNEL_EXTENSION_POINTS.md:221`).**
One docs sheet; no implementation. Covers: per-cell omega passing convention
(B-6); placement of MRT/cumulant kernels (CollisionKind branching and transform
matrix location); per-link wall-distance sparse structure for curved boundaries
(Bouzidi); ConvectiveOutflow alternative under in-place streaming (generalize
the GPU edge-stash scheme); BC fallback for Outflow × solid adjacency (the
permanent solution of A-3). Each extension's DoD includes existing-config bit
invariance. Effort: M (through review approval).

Handoff note: a design note landed 2026-07-06 (`Add kernel extension-point
design note`, commit 1758814). Confirm coverage of the above five topics
before closing B-8.

---

## WP-C — Scale (MPI / GPU), still open

R-Phase 3 items owner-scheduled to run alongside M-E. C-2, C-3, C-8 depend on
B-1; C-12, C-13 also depend on B-1.

### C-1 [S1] Localize MPI setup (kill global-array replication)
**STATUS (2026-07-07 audit): STILL OPEN**
Target: `crates/lbm-core/src/{dist.rs, solver.rs}`.
`MpiSolver::new` requires global compact arrays (solid/wall_u) on all ranks;
`build_wall_rims` also generates the whole domain; `set_solid` is a
global-all-cells×all-ranks call. At 10⁹ grid, `wall_u ≈ 24 GB/rank` = guaranteed
OOM in the weak-scaling configuration. Structural blocker that manifests with
grid size, not rank count (current tests undetected because n≤8, small grid).

Fix: closure-taking `MpiSolver::new_with(solid: impl Fn(x,y,z)->bool, …)`
(local eval, mirroring `init_with`) + local `build_wall_rims` + batch
`set_solids_where(pred)`. Legacy API stays for small scale.

Acceptance: per-rank peak RSS = O(N/P)+constant (measured); T13-MPI bit match
via the new API. Effort: M.

### C-2 [S2] Exchange overlap + boundary-shell disjointness fix
**STATUS (2026-07-07 audit): PARTIALLY RESOLVED — disjoint boundary shells and probe split-invariance tests landed (`crates/lbm-core/src/backend.rs:79`, `crates/lbm-core/tests/t13_split_invariance.rs:286`); remaining: no `post_f()`/`finish_f()` overlap API and `WgpuBackend` still rejects two-pass streaming (`crates/lbm-core/src/gpu/backend.rs:2044`).**
Target: `dist.rs`, `solver.rs`.
Add `post_f()/finish_f()` to `HaloExchange` (InProcess completes immediately).
Rewire step to collide → post → interior stream → finish → boundary shell.
**Premise fix**: two_pass's `boundary_shells` overlap on 1-width axes; the
probe is double-counted (field is idempotent so invisible in T13). E8 confirmed
the on/off probed_force ratio = 2.000 on `[64,1,1]`. Fix shells disjoint before
connecting.

Acceptance: two_pass on/off probed_force match on 1-width axis; T13-MPI PASS;
exchange-wait occupancy measurably reduced. Effort: L. Depends: B-1 item 3
(withdrawal of the `stream(range)` assert).

### C-3 [S2] Parallel I/O (per-rank raw + manifest)
**STATUS (2026-07-07 audit): STILL OPEN**
Target: `dist.rs`.
Short-term: each rank writes its block; rank0 writes only the manifest. Mid-
term: MPI-IO subarray (rsmpi `File` support to confirm). Keep gather_* as
validation-only.

Acceptance: rank0 peak RSS at output = O(N/P); output time non-increasing with
rank count. Effort: M-L.

### C-4 [S2] Defer probe Allreduce
**STATUS (2026-07-07 audit): STILL OPEN**
Target: `dist.rs` (3-double Allreduce every step when probe enabled).
Fix: make `probed_force()` collective and cache with `time` key (reduce only
when queried). Acceptance: per-step collective disappears from profile; mpi_t13
PASS. Effort: S.

### C-5 [S2] Persistent exchange buffers (stepping stone to GPU-aware MPI)
**STATUS (2026-07-07 audit): PARTIALLY RESOLVED — `MpiExchange` now owns reusable send/recv buffers and `Solver` owns persistent psi planes (`crates/lbm-core/src/dist.rs:266`, `crates/lbm-core/src/solver.rs:1274`); remaining: steady exchange still creates request vectors in `transfer_axis`, so zero-allocation acceptance is not proven (`crates/lbm-core/src/dist.rs:337`).**
Target: `dist.rs`, `solver.rs` (ψ plane allocated every step).
Fix: `MpiExchange` holds send/recv buffers per face×kind and reuses them; ψ
plane / staging fields likewise. Consolidate buffer ownership in MpiExchange
(single swap point to GPUDirect later).
Acceptance: zero heap allocation during steady steps; bench_mpi non-regressing;
T13-MPI bit match. Effort: S-M.

### C-6 [S2] Inter-rank spec consistency check
**STATUS (2026-07-07 audit): RESOLVED — rank spec hashes are Allreduce-compared with named mismatch failures and tests/scripts cover the changed-item path (`crates/lbm-core/src/dist.rs:244`, `crates/lbm-core/src/dist.rs:1176`, `scripts/test_mpi.sh:50`).**
Target: `dist.rs`.
A mismatch in nu, faces, or mask content is undetected because message length
matches — emits a plausible field discontinuous at seams. Fix: at build time,
Allreduce-compare spec-normalized byte sequence + mask FNV hash (min/max);
abort with item name on mismatch. Effort: S.

### C-7 [S2] MPI thread level = Funneled
**STATUS (2026-07-07 audit): PARTIALLY RESOLVED — `bench_mpi` requests `Threading::Funneled` and logs the provided level (`crates/lbm-core/examples/bench_mpi.rs:67`); remaining: no shared `dist::init_mpi()` and `mpi_t13` still calls `mpi::initialize()` (`crates/lbm-core/examples/mpi_t13.rs:385`).**
Target: `examples/mpi_t13.rs`, `bench_mpi.rs`, `backend.rs`.
rsmpi 0.8.1 `initialize()` = `Threading::Single`; default feature `parallel`
launches rayon above 16,384 cells. Real sizes get multi-threaded execution
under a Single declaration (MPI-spec violation; UCX/OFI can corrupt/hang).
Current tests coincidentally serial at ≤6,144 cells/rank.
Fix: add `dist::init_mpi() -> Universe` (requests Funneled; explicit error if
provided is insufficient AND parallel is enabled); migrate examples/guide.
Acceptance: 2 ranks × forced rayon (parallel_min_cells lowered) → T13-MPI-
equivalent PASS; provided ≥ Funneled logged. Effort: S.

### C-8 [S2] Distributed checkpoint/restart
**STATUS (2026-07-07 audit): STILL OPEN**
On top of B-5's snapshot API: collective `MpiSolver::save(dir)/load(world, dir,
backend)` — per-rank raw + rank0 manifest, spec-hash and decomp consistency
validation. Deviation-storage f raw for bit-match on resume.
Acceptance: 50 steps → save → load → 50 steps = 100 steps continuously (bit
match); manifest mismatch is an explicit error. Effort: M. Depends: B-5, C-6.

### C-9 [S1] GPU submit-chunk time calibration + device-lost Result-ification
**STATUS (2026-07-07 audit): RESOLVED — submit chunk calibration, poll/device-lost error paths, and tests are present (`crates/lbm-core/src/gpu/backend.rs:249`, `crates/lbm-core/src/gpu/backend.rs:447`, `crates/lbm-core/src/gpu/backend.rs:470`, `crates/lbm-core/src/gpu/backend.rs:2346`).**
Target: `gpu/backend.rs` (`submit_chunk: 200` fixed; `expect` panic on device
loss).
200 steps ≈ ~1000 dispatch in a single submit → on a slow GPU, exceeds Windows
TDR (default 2 s) → device removed → process panic; no recovery.
Fix: auto-calibrate to target 1 submit ≈ 100–250 ms from first-chunk
measurement (upper limit 200 preserved). Result-ify `wait_idle`/`map_staging`
to `Result<_, GpuError>`; propagate device lost. Capture reason via
`set_device_lost_callback`.
Acceptance: bench_gpu MLUPS regresses ≤3%; calibration unit test; poll failure
returns Err. Effort: M. E10 baseline (M5 Max): 2048² TGV 5,719 MLUPS →
147 ms/submit; TDR extrapolation for ~15× slower GPUs holds.

### C-10 [S2] Pre-validate GPU resource limits
**STATUS (2026-07-07 audit): RESOLVED — GPU resource planning validates `Q*n`, buffer/storage limits, and workgroup counts with unit tests (`crates/lbm-core/src/gpu/backend.rs:283`, `crates/lbm-core/src/gpu/backend.rs:341`, `crates/lbm-core/src/gpu/backend.rs:2369`).**
Target: `gpu/backend.rs`.
At head of `alloc`, check required bytes vs `device.limits()`, `Q*n ≤ u32::MAX`
(D3Q19 overflows at 226M cells), and dispatch count ≤65,535. Err with reason.
Acceptance: over-limit grid → `GpuSolver::new` explicit error; T14 green.
Effort: S.

### C-11 [S2] GPU diagnostic-path efficiency
**STATUS (2026-07-07 audit): PARTIALLY RESOLVED — `f_cache` is `Arc`, `FluidCells` uses `host_solid`, and sync combines populations/moments/probe in one readback (`crates/lbm-core/src/gpu/backend.rs:547`, `crates/lbm-core/src/gpu/backend.rs:2154`, `crates/lbm-core/src/gpu/backend.rs:2269`); remaining: reductions still read populations back to the host, so GPU-side 2-stage reduction is not landed (`crates/lbm-core/src/gpu/backend.rs:2137`).**
Target: `gpu/{backend,solver}.rs`.
Arc-ify `f_cache` (eliminate cloning). Make FluidCells readback unnecessary
(host_solid suffices). Consolidate sync's 3 readbacks into 1 encoder/1 wait.
Add GPU-side 2-stage reduction as an M-E fast mode (keep host f64 path for
T14).
Acceptance: sync+diagnostics triple ≥3× faster at 2048²; T14 diagnostic values
bit-identical. Effort: S (+M).

### C-12 [S2] FP16 plumbing
**STATUS (2026-07-07 audit): RESOLVED — `GpuStorage::F16`, `SHADER_F16` context requirement, f16 WGSL load/store, and T16 frozen-band tests are present (`crates/lbm-core/src/gpu/backend.rs:144`, `crates/lbm-core/src/gpu/backend.rs:398`, `crates/lbm-core/src/gpu/wgsl.rs:228`, `crates/lbm-core/tests/t16_fp16_storage.rs:15`).**
Target: `gpu/backend.rs`, `gpu/wgsl.rs`.
Conditional `SHADER_F16` request; `generate::<L>(cfg: KernelCfg { storage:
F32|F16 })` (confine change to buffer decl + load/store wrapper, keep arithmetic
in f32); hold element size in `GpuFields`. Non-supporting adapters return an
explicit error (no silent fallback).
Acceptance: T16 established (band-freeze the f16 storage degradation); MLUPS
≥1.5× vs f32 at 2048² TGV. Effort: M. Depends: B-1.

**Landed 2026-07-06** (per HANDOFF §3 claims-ledger snapshot): FP16 storage ×2
capacity GREEN. Bands frozen (TGV transient 1.401e-1 vs band 2e-1; cavity
steady 2.579e-3 vs band 5e-3). ~2.0× MLUPS @2048²; D3Q19 f16 >5 GLUPS. C-12
retained here for historical spec; body is done.

### C-13 [S2] scenario `backend: "gpu" | "auto"`
**STATUS (2026-07-07 audit): RESOLVED — scenario schema has `auto|cpu|gpu`, auto thresholding/capability rejection, CLI GPU dispatch, and `build_gpu2d` (`crates/lbm-scenario/src/lib.rs:95`, `crates/lbm-scenario/src/lib.rs:126`, `crates/lbm-cli/src/runner.rs:63`, `crates/lbm-scenario/src/lib.rs:1304`).**
Target: `crates/lbm-scenario/src/lib.rs`, `crates/lbm-cli/src/main.rs`.
Unlock in feature-gpu build; validate capability (reject f64/3D/unsupported
BC with reason); `auto` selects by measured threshold (e.g. n≥256²) and logs.
Acceptance: cavity/cylinder presets complete with `backend:"gpu"`; field diff
from CPU within T14; `f64`+`gpu` rejected with reason. Effort: M. Depends: B-1.

**Landed 2026-07-06** (HANDOFF §3, commit 1a14d90): explicit `backend:"gpu"`
GREEN. Retained for historical spec.

### C-14 [S2] GPU memory-footprint reduction
**STATUS (2026-07-07 audit): PARTIALLY RESOLVED — optional dummy wall/force allocation exists and upload expands on demand (`crates/lbm-core/src/gpu/backend.rs:1675`, `crates/lbm-core/src/gpu/backend.rs:1712`, `crates/lbm-core/src/gpu/backend.rs:1596`); remaining: default backend `alloc()` still requests full wall/force buffers (`crates/lbm-core/src/gpu/backend.rs:1671`, `crates/lbm-core/src/gpu/backend.rs:1892`).**
Target: `gpu/backend.rs`.
Dummy-buffer-ize force_field/wall_u when absent; lazy-allocate and right-size
staging. Current +284 MB (~1.9×) at 2048² TGV — 3D scaling directly hits the
max grid of an 8-16 GB card.
Acceptance: without force field or solid, allocation ≤ 2×f + moments + mask +
O(perimeter); T14 green. Effort: S.

### C-15 [S3] GPU small-grain bundle
**STATUS (2026-07-07 audit): PARTIALLY RESOLVED — raised adapter limits, naga validation tests, BcParams table test, and `GpuContext::new -> Result` landed (`crates/lbm-core/src/gpu/backend.rs:426`, `crates/lbm-core/src/gpu/wgsl.rs:1410`, `crates/lbm-core/src/gpu/wgsl.rs:1443`, `crates/lbm-core/src/gpu/backend.rs:392`); remaining: probe CAS nondeterminism is documented in WGSL comments, not in `docs/GPU_EVALUATION.md` (`crates/lbm-core/src/gpu/wgsl.rs:786`).**
- (a) Raise limits (`max_storage_buffers_per_shader_stage` — step kernel at
  the default upper limit 8, adding 1 storage buffer panics at runtime).
- (b) `naga` parse+validate unit test of the generated WGSL (no GPU required,
  golden file).
- (c) Single-table generation or collation test between Rust index and WGSL
  field of `BcParams`.
- (d) Result-ify `GpuContext::new` to `Result<_, GpuInitError>` (with
  `adapter_info`; makes auto fallback diagnosable).
- (e) One-line note in GPU_EVALUATION.md on the non-determinism of the probe-
  force f32 CAS addition.
Effort: S×5.

### C-16 [S3] MPI small-grain bundle
**STATUS (2026-07-07 audit): RESOLVED — surface-minimizing `choose_decomp`, non-divisible rank coverage in `test_mpi.sh`, 3D/weak/strong `bench_mpi` modes, and deterministic mass are present (`crates/lbm-core/src/dist.rs:104`, `scripts/test_mpi.sh:43`, `crates/lbm-core/examples/bench_mpi.rs:7`, `crates/lbm-core/src/dist.rs:810`).**
- (a) `choose_decomp` (minimum surface area) + arbitrary-n support of mpi_t13
  (cases n=3,5,6 and non-divisible dims).
- (b) 3D/D3Q19 and strong/weak-mode-ization of bench_mpi (current 2D-strip
  fixed cannot serve R3 main measurement).
- (c) Deterministic total sum `total_mass_deterministic()` (fixed-order
  composition of global row-wise partial sums), optional.
Effort: S+S+M.

**Landed 2026-07-06** for (b): bench_mpi 3D + weak modes ready per HANDOFF §4;
consumed by CLUSTER_OPTIONS. (a) and (c) still open.

---

## WP-D — Validation / Process, still open

### D-1..D-5, D-8..D-12
R-Phase 1 completed D-6 (competitive spec / PLAN reconciliation) and D-7
(VALIDATION.md T13/T14 sections). Remaining:

- **D-1 [S1]** CI three-stage setup (push/PR default; nightly `--include-
  ignored`; self-hosted GPU/MPI).
  **STATUS (2026-07-07 audit): PARTIALLY RESOLVED — `.github/workflows/ci.yml` now runs release tests, GPU compile, MPI, web build, and defines a heavy-validation job (`.github/workflows/ci.yml:21`, `.github/workflows/ci.yml:23`, `.github/workflows/ci.yml:25`, `.github/workflows/ci.yml:43`); remaining: no `schedule:` trigger is present and no self-hosted GPU/MPI runner is configured (`.github/workflows/ci.yml:4`).**
  Currently `.github/` absent; all quality claims are manual snapshots. Effort:
  M.
- **D-2 [S1]** Cargo-test MPI logic (`dist.rs` has 0 `#[test]`).
  **STATUS (2026-07-07 audit): PARTIALLY RESOLVED — dist/halo pure tests and one-rank MPI smoke exist (`crates/lbm-core/src/dist.rs:1099`, `crates/lbm-core/src/dist.rs:1191`); remaining: `test_mpi.sh` is in CI but not an active nightly because the workflow has no schedule trigger (`.github/workflows/ci.yml:4`, `.github/workflows/ci.yml:43`).**
  Pack/unpack and phase-plan as pure functions (InProcess reference, tolerance
  ==0.0). Turn `test_mpi.sh` into a nightly job; bring 1-rank self-exchange
  smoke into cargo. Effort: M.
- **D-3 [S1]** Absolute physics validation of GPU + skip mechanism.
  **STATUS (2026-07-07 audit): STILL OPEN**
  TGV convergence order ≥1.7 and cavity Ghia RMS ≤0.02U — GPU-direct, f32-
  calibrated bands, frozen. Adapter absence = skip; `LBM_REQUIRE_GPU=1`
  promotes to fail. Non-support of 3D GPU / multiphase explicit in
  VALIDATION.md limitations. Effort: M.
- **D-4 [S1]** f32×3D validation.
  **STATUS (2026-07-07 audit): RESOLVED — `t15_3d_f32.rs` covers z-invariant degeneration, TGV3D decay, and mass drift for the f32 D3Q19 product path (`crates/lbm-core/tests/t15_3d_f32.rs:1`, `crates/lbm-core/tests/t15_3d_f32.rs:82`, `crates/lbm-core/tests/t15_3d_f32.rs:161`, `crates/lbm-core/tests/t15_3d_f32.rs:195`).**
  `Sim3Handle::F32` is a product path with 0 f32 tests. Add: t15-1 z-invariant
  2D degeneration (f32 rel ≤1e-5), t15-4 TGV3D decay rate ±2%, mass drift
  ≤1e-5/10³ step. Freeze in VALIDATION T15. Effort: S.
- **D-5 [S1]** Extend equivalence horizon + native absolute validation.
  **STATUS (2026-07-07 audit): RESOLVED — `d5_long_horizon.rs` contains the ignored Re=100 20k-step compat/native cavity equivalence and native `Solver` TGV convergence test (`crates/lbm-core/tests/d5_long_horizon.rs:75`, `crates/lbm-core/tests/d5_long_horizon.rs:149`).**
  Currently ≤1e-11 vs V1 at 500 steps only; Ghia steady is ~99k steps (200×
  the horizon); V2 native (3D / scenario path) outside the chain. Add: cavity
  Re=100 20k-step equivalence (≤1e-9, may be `#[ignore]`); TGV convergence
  order directly on native `Solver` (detects facade/core drift). Effort: S.
- **D-8 [S2]** codex adversarial order #7 for T14/T15 (initial discontinuity
  above a boundary face; probe touching a face; u near MAX_SPEED; z-degeneration
  break; extreme aspect ratio; off-center sphere).
  **STATUS (2026-07-07 audit): RESOLVED — adversarial T14/T15 tests cover the named cases (`crates/lbm-core/tests/t14_adversarial.rs:5`, `crates/lbm-core/tests/t14_adversarial.rs:122`, `crates/lbm-core/tests/t14_adversarial.rs:189`, `crates/lbm-core/tests/t14_adversarial.rs:275`, `crates/lbm-core/tests/t15_adversarial.rs:173`, `crates/lbm-core/tests/t15_adversarial.rs:290`).**
  Effort: M.
- **D-9 [S2]** Performance regression detection.
  **STATUS (2026-07-07 audit): PARTIALLY RESOLVED — core benchmark coverage exists in `bench_backends` (`crates/lbm-core/examples/bench_backends.rs:1`); remaining: no `--check`, no `docs/bench_history/` JSON history, and no nightly `probe_state_hash` enforcement found by grep.**
  Core version of `bench_mlups` (V2 CPU perf still cites V1 Phase 9 numbers) +
  `--check` mode (fail on ±25% deviation from frozen value) + nightly JSON
  history in `docs/bench_history/`. Nightly probe_state_hash enforcement.
  Effort: M.
- **D-10 [S2]** `lbm verify` (few-minutes validation subset comparing against
  acceptance bands) + `docs/LIMITATIONS.md` shipped paired with releases.
  **STATUS (2026-07-07 audit): STILL OPEN**
  Effort: L (verify) + S (limitations).
- **D-11 [S2]** wasm smoke test (`wasm-pack test --node`: TGV 100-step mass
  conservation + native f64 match ≤1e-12; schema round-trip on the web side).
  **STATUS (2026-07-07 audit): PARTIALLY RESOLVED — a wasm-bindgen TGV smoke compares wasm f32 fields bitwise to native compat f32 (`crates/lbm-wasm/src/lib.rs:360`, `crates/lbm-wasm/src/lib.rs:385`); remaining: no f64 ≤1e-12 native match or web-side schema round-trip found.**
  Effort: M.
- **D-12 [S2]** CpuScalar performance-gap explicitation → M-E.
  **STATUS (2026-07-07 audit): RESOLVED — `CpuSimd` equivalence gates and V1 ratio notes are present (`crates/lbm-core/tests/backend_simd_equiv.rs:1`, `crates/lbm-core/examples/bench_backends.rs:12`).**
  Q-major port of V1 `step_band` in `CpuSimd`; until then, publish co-noted V1
  single-core ratio. Acceptance: `CpuSimd` bit-matches `CpuScalar`; single-core
  ≥ 0.9x of V1. Effort: L (M-E body).

---

## Sequencing

Per HANDOFF §4:

- **R-Phase 2 residuals owner-scheduled into M-E/M-F**: B-2, B-3, B-6..B-8,
  R2-C (mechanical TRT port with ANOM-P2-001 fix — order text staged in
  `scratchpad/order-r2c.txt`; PM-owned). B-1 (highest-value structural
  prerequisite) tracked separately as the M-E gate.
- **R-Phase 3 (scale)**: C-1, C-4..C-7, C-16(a,c) can run in parallel with
  M-E once B-1 lands. C-2/C-3/C-8 wait on B-1. C-12/C-13 body done; C-14/C-15
  remain.
- **WP-D**: D-1/D-2/D-9 are quickest wins for closing the "all quality claims
  are manual snapshots" gap; D-3/D-4/D-5 gate GPU/f32 claims.

Validation tests (D-3/D-4/D-8 and each acceptance) are authored adversarially
by codex from the VALIDATION.md T13/T14 spec; test worktrees never share with
implementation worktrees.
