> **STATUS: IMPLEMENTED 2026-07-06 via 20d0e10 (central-moment collision on SIMD + GPU, MF-alpha stage 3)**

# MF-α implementation plan: D3Q27 lattice + cumulant/central-moment collision

Produced 2026-07-06 by a PM-commissioned design survey (Plan agent) of
crates/lbm-core, for the M-F table row MF-α (docs/PLAN.md) and FR-CORE-01/02
(docs/REQ_STIRRED_REACTOR.md). Line references are as of main 004d77b.

## Q1 — Lattice trait and every Q-hardwired site

**What a D3Q27 impl must provide** (trait at `src/lattice.rs:116-151`): `D=3`,
`Q=27`, `C` (27 velocities), `W` (8/27, 2/27 ×6, 1/54 ×12, 1/216 ×8), `OPP`,
`PAIRS` (13 pairs), `FACE_UNKNOWNS` (9 per face: `n`, `n±t1`, `n±t2`, and the 4
corners `n±t1±t2`). All derived tables come from the existing generic
`const fn`s (`opp_table` lattice.rs:157, `pairs_table` :174, `face_unknowns`
:192) — zero new machinery. `Q_MAX` is already 27 (lattice.rs:366), so kernel
scratch arrays (`[T; Q_MAX]` kernels.rs:110,172-173, backend_simd.rs:1861) and
`KParams` (params.rs:108-110) already fit. Recommended ordering: extend the
D3Q19 convention (lattice.rs:279-291) — rest, 6 axes, 12 edge diagonals, then
8 corners in adjacent opposite pairs, making D3Q19 a strict prefix and
preserving the `OPP[q]=q±1` property (lattice.rs:544-550).

**Q=9/Q=19 hardwired sites (exhaustive):**

| Site | Nature |
|---|---|
| `src/kernels.rs:531-544` | `zou_he_face_3d` asserts exactly 5 unknowns (guarded panic anticipating D3Q27). The 3D Zou–He derivation (kernels.rs:509-619) reconstructs only 5 directions; Q27 needs 9 incl. corner terms — new derivation, not a loop change. |
| `src/gpu/wgsl.rs:77-142` | `BcParams` uniform has exactly `unk0..unk4` + `unk_count` (:104-109) and D3Q19-specific slots (`q_pp/q_pm/...`); the generated `bc` kernel branches on `unk_count == 3u/5u` (:785, :819, :939) and loops `for k in 0..5` (:956, :981). |
| `src/gpu/backend.rs:1057-1060` | `assert!(unk.len() == 3 || unk.len() == 5)` for GPU open faces. |
| `src/lattice.rs:473` | Test helper `check_face_unknowns` hardcodes `expected_count = if L::D == 2 { 3 } else { 5 }` — must become 9 for Q27. |
| `src/solver.rs:576-590` | `lattice_name()`/`lattice_id()` map only (2,9)/(3,19); D3Q27 checkpoints would silently write `"unknown"`/id 0. |
| `crates/lbm-scenario/src/lib.rs:743-746` | 3D engine type alias pinned to `D3Q19` (FR-CORE-01 selectability is a facade concern, deferrable). |

**Generic already (verified, no change needed):** streaming `stream_row`
(kernels.rs:230-305; all |c|≤1 so corner directions work); halo plane packing
(halo.rs:224-260) and the x→y→z extended-layer corner forwarding
(halo.rs:154-194 handles (±1,±1,±1) by two-hop forwarding); wall rims
(solver.rs:798, mask-only); SIMD `step_band` (ring sized `L::Q` planes,
backend_simd.rs:1100; blocked collide auto-selected by `L::Q > 9`,
backend_simd.rs:183-186); WGSL `step`/`moments`/edge-stash generation
(wgsl.rs:243-536, `stash_len` :147); reductions, Bouzidi, moments,
`outflow_face`/`convective_face` (kernels.rs:623-699).

## Q2 — Collision architecture

`CollisionKind` (src/params.rs:15-23, Bgk | Trt{magic}) is **erased at solver
build**: solver.rs:1035 collapses it to `(omega_p, omega_m)`; `StepParams`
(params.rs:82-91) carries only the two rates; `collide_row`
(kernels.rs:127-216) is a single TRT-pair loop (BGK ≡ ω⁻=ω⁺). Dispatch points
a `Cumulant` variant must reach:

- **CpuScalar**: backend.rs:351-391 → sibling `collide_row_cumulant`.
- **CpuSimd**: single funnel `collide_span_dispatch` (backend_simd.rs:508-524),
  called from the fused band pass AND from `capture_stale`
  (backend_simd.rs:1861-1893, which re-collides open-face cells — cumulant must
  be used there too or convective-BC memory drifts).
- **GPU**: collide section of `emit_step_entry` (wgsl.rs:374-421) + `Params`
  uniform (wgsl.rs:561-567 has only `omega_p/omega_m/cp/cm`).

**Layout/pass-structure fit:** yes. SoA q-major already gathers all Q values
per cell per collision; cumulant collision is cell-local, so the invariant
`collide → halo → stream → BC → moments` and the fused ring are untouched.
Exhaustive matches on `CollisionKind` (solver.rs:510-514 spec hash,
dist.rs:184-189 checkpoint hash, solver.rs:317-321 validate) fail to compile
until extended — good. Real design item: **Guo forcing** is per-direction
inside collide (kernels.rs:181-186, 199-213 with `cp/cm` prefactors); the
cumulant operator needs forcing transformed to (central-)moment space with the
same half-force velocity convention (`moments_row` kernels.rs:352-356).

## Q3 — omega/relaxation field

`set_omega_field` (solver.rs:2105-2135) → per-part compact
`fields.omega_field` (fields.rs:197) → consumed only inside collide kernels:
kernels.rs:191-196 (`op = omega[x]`, `cp = 1 - op/2`; **ω⁻ stays uniform**)
and identically in backend_simd.rs (:289-294 / blocked). WALE installs it
(les.rs:113), bounded to `Fields = SoaFields<T>` — CPU-only today; the GPU
shader has no omega binding (being added by the cx/gpu-wale order). Cumulant
composition is clean: the per-cell field maps to the second-order (shear)
cumulant relaxation rate ω₁ per cell (the "collision kernels only replace the
local omega_plus fetch" contract, solver.rs:2102-2104), higher-order rates
uniform. REQ §4.6 anticipates this.

## Q4 — WGSL generator at Q=27

Generic: prologue, collide, push-stream, edge stash, moments, `stash_len` all
generate from lattice tables — a Q27 `step`/`moments` shader would emit and
validate today. Breaks/risks:

1. **`bc` entry point + `BcParams`** — the 5-unknown ceiling. Needs
   `unk0..unk8`, corner slots, Q27 Zou–He formula emission.
2. No workgroup memory used (all registers) — nothing hard breaks; the
   unrolled step kernel grows 19→27 loads + push blocks (~1.4-1.5× code);
   register pressure / Metal compile time is a perf risk, not correctness.
3. Buffer sizes: f buffers scale as `Q*n` with a `u32::MAX` element guard
   (`GpuResourcePlan::for_grid`, gpu/backend.rs:289-300) — ceiling drops from
   ~226M to ~159M cells (~540³), fine for the ≤256³ dev line. Memory 27/19 ≈
   1.42× (REQ budget table already assumes D3Q27).
4. f16 storage path (`storage_load/store`) is orthogonal and carries over.

## Q5 — Test extensibility

Cheap (generic helpers; add an instantiation or genericize a
`type S3 = Solver<D3Q19,...>` alias):
- lattice.rs invariant suite (:372-495) — one new `d3q27_invariants` test; fix
  :473. 4th-order isotropy passes for Q27; ADD the third-order-diagonal /
  6th-order checks that DISTINGUISH Q27 from Q19 (D3Q19's known deficit).
- tests/t15_3d.rs `tgv3d_short` + `t15_4_tgv3d_diffusive_convergence`
  (:597-693) — the natural stage-1 TGV3D order gate.
- tests/t13_split_invariance.rs `t13_tgv3d_2x2x2_split_invariant` (:306).
- tests/d3q19_smoke.rs, tests/backend_simd_equiv.rs case 7 (3D TGV),
  kernels.rs bit-exact feq fixed-point test (:774-784) — generic bodies.
- tests/t14_3d_backend_equiv.rs — GPU gate, reusable at stage 3.

Not cheap: t15_1b/1c Zou–He degeneracy and any open-face case — blocked on the
9-unknown Zou–He (stage 4).

## Staged order plan

Survey correction to the naive staging: **CpuSimd is already fully generic for
BGK/TRT at Q=27**, so stage 1 gets SIMD for free and gates it with
backend_simd_equiv; "SIMD in stage 3" applies only to the cumulant operator.

**Stage 1 — D3Q27 lattice, BGK/TRT reuse (CpuScalar + CpuSimd).** ~400-550 lines.
- src/lattice.rs: D3Q27 tables + const derivations + invariant tests (incl.
  :473 fix, Q19-vs-Q27 discriminating moment checks).
- src/solver.rs: lattice_name/lattice_id entries; build-time SpecError guard
  rejecting open FaceBC when the lattice face-unknown count isn't 3/5 (turn
  the kernels.rs:538 runtime panic into a construction error).
- src/lib.rs export + docs; src/halo.rs:7, src/lattice.rs:1 doc updates.
- SEPARATE VALIDATION-TEST ORDER (adversarial, own worktree, after stage-1
  impl lands): genericize tests/t15_3d.rs, tests/t13_split_invariance.rs,
  tests/backend_simd_equiv.rs 3D case; new tests/d3q27_smoke.rs.
- Gate: d3q27_invariants + TGV3D order ≥ ~1.7 + T13 2×2×2 split invariance +
  SIMD equivalence.

**Stage 2 — Cumulant collision, CPU scalar reference.** ~700-1000 lines.
- src/params.rs: `CollisionKind::Cumulant { rates }` (central-moment cascaded
  form is acceptable first per FR-CORE-02; cumulant shares the transform
  scaffolding); StepParams carries the operator; KParams extension.
- src/kernels.rs: collide_row_cumulant (forward transform, per-rate
  relaxation, back-transform, moment-space Guo forcing), ω-field consuming the
  shear rate per Q3.
- src/solver.rs (:317 validate — rate positivity per FR-CORE-02; :510 hash),
  src/dist.rs (:184 hash), src/backend.rs (:351 dispatch). Stages 1+2 both
  touch solver.rs/kernels.rs — sequential ordering satisfies the same-file rule.
- Tests (separate adversarial order): Galilean-invariance band (advected
  TGV/shear wave at nonzero mean velocity, cumulant vs BGK error ratio),
  rotational isotropy of decay rates, TGV3D under cumulant, rate-range bench,
  ω-field composition extending tests/wale_les.rs.

**Stage 3 — Cumulant on SIMD + GPU.** ~700-900 lines.
- src/backend_simd.rs: blocked cumulant span kernel behind
  collide_span_dispatch (:508) INCLUDING capture_stale (:1861); extend
  backend_simd_equiv with a measured cumulant tolerance (no pair-reassociation
  equivalence argument exists for cumulant — the 1e-11 TRT contract does not
  carry; measure a fresh bound).
- src/gpu/wgsl.rs cumulant emission + Params rates; gpu/backend.rs plumbing;
  extend t14_3d_backend_equiv.

**Stage 4 (parallel-capable after stage 1) — D3Q27 open-face BCs + facade.**
~550-700 lines.
- src/kernels.rs: 9-unknown NEBB/Zou–He derivation (corner unknowns); remove
  the stage-1 guard.
- src/gpu/wgsl.rs + gpu/backend.rs: BcParams → unk0..8 + corner slots, bc
  kernel regeneration, :1057 assert.
- crates/lbm-scenario/src/lib.rs: lattice selection knob (FR-CORE-01).
- Tests: t15_3d Zou–He degeneracy for Q27, duct case.

## Top-3 hardwiring risks

1. **Open-face BC stack is structurally D3Q19-shaped end to end** — CPU
   derivation (kernels.rs:509-619), GPU uniform layout + generated kernel
   (wgsl.rs:77-142, :717-996), GPU word-packing (gpu/backend.rs:1011-1110).
   The 9-unknown closure is new physics derivation work; the GPU uniform is a
   binary layout change touching three files at once.
2. **Collision-kind erasure**: everything downstream of solver.rs:1035 assumes
   the operator ≡ (omega_p, omega_m, cp, cm) — StepParams, KParams, WGSL
   Params, SIMD stale-capture re-collide. Keep BGK/TRT paths byte-identical
   and add cumulant as a parallel branch (the V1-port bit-exactness contracts,
   kernels.rs:3-14, must not perturb).
3. **Identity/equivalence infrastructure**: checkpoint lattice_name/lattice_id
   silently degrade at Q27; collision hashes need new arms; the
   backend-equivalence tolerance philosophy has no cumulant analogue — an
   unmeasured tolerance either masks bugs or flakes CI.
