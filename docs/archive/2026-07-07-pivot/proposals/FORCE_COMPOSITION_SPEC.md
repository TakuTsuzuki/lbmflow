# Force-Composition Graph — Freeze Spec (pre-W-VOF prerequisite)

Status: **FROZEN** (reviewer directive, 2026-07-07). This document freezes the
force-composition graph *before* W-VOF implementation begins, so that
device-resident force composition cannot silently regress to host staging when
new force terms (surface tension, scalar buoyancy, particle reaction) land.

It is a design contract, not an implementation order. Every "composed today"
claim below was verified against the code at the cited `file:line`. Consistency
target: `docs/proposals/WVOF_IMPL_SPEC.md` (decisions D4/D5/D7, §4.4) and
`docs/PHYSICS.md` (2026-07-07 backend-gravity entry, §"Well-balanced gravity
composition point").

---

## 0. The graph

Every body force in the engine enters collision through **one** Guo forcing
term. There is exactly one per-cell accumulator, `SoaFields::force_field`
(compact-core `Vec<[T;3]>`, optional; `None` = no per-cell force), plus the
scalar uniform force `StepParams::force` and the per-mass gravity
`StepParams::gravity`. The frozen total, per fluid cell `x`, per step:

```
F_total(x) =   F_user_uniform                     [StepParams::force,  scalar, backend-side]
             + F_cell_field(x)                     [force_field[x],     backend-side]
             + (rho_or_rho_phi(x) - rho_ref) * g   [StepParams::gravity, backend-side; rho_ref=0 today]
             + F_surface_tension(x)                [W-VOF: into force_field, solver pre-pass]
             + F_scalar_buoyancy(x)                [W-SCAL: into force_field, solver pre-pass]
             + F_particle_reaction(x)              [RESERVED (two-way); not implemented]
             + F_rotating_ibm(x)                   [into force_field, solver pre-pass]
```

The single arithmetic point at which `F_user_uniform`, `F_cell_field` and
gravity are summed and handed to Guo forcing is `KParams::force_at`
(`crates/lbm-core/src/params.rs:232-251`) on the CPU backends, and the
WGSL-emitted collide preamble (`crates/lbm-core/src/gpu/wgsl.rs:269-286`) on the
GPU backend. Every other term (surface tension, scalar buoyancy, IBM) is a
**solver-orchestration-level pre-pass** that writes into `force_field` *before*
the backend step, so it reaches Guo forcing through `F_cell_field` and never
touches a collision or BC kernel. This is the invariant W-VOF must preserve.

The Guo half-force velocity correction is applied inside `force_at`'s consumer:
`u = (m + F_total/2)/rho` (CPU `moments_row`; GPU `wgsl.rs:302-305`). **The full
`F_total` sum must be formed before that half-force division** — see §1.

---

## 1. Per-term specification

Legend for (a): **backend-side** = composed inside the backend step from data
already resident in `Backend::Fields`, no per-step host round-trip;
**host-side (pre-pass)** = composed by the solver into `force_field` on host
staging before the backend step; **host-staged overlay** = written into
`force_field`, uploaded, then unwound after the step; **not implemented**.

### T1. `F_user_uniform` — uniform body force
- **(a) Today:** backend-side. `StepParams::force: [T;3]`
  (`params.rs:150`), summed first in `force_at` (`params.rs:235,240,245,249`);
  GPU seeds `fvx = P.fx` (`wgsl.rs:270-272`). Never staged.
- **(b) After W-VOF:** unchanged, backend-side. It is a scalar uniform, carries
  no `phi` dependence.
- **(c) GPU:** supported today (`wgsl.rs:270-272`); no change.
- **(d) MPI/halo:** none. Cell-local constant; no neighbor data.
- **(e) Pass position:** read inside `collide`, before the half-force
  correction. First summand of `F_total`.

### T2. `F_cell_field` — caller-owned per-cell force field
- **(a) Today:** backend-side. `force_field` is uploaded once and read in
  `collide`/`moments` via `force_at` (`params.rs:236,241`); GPU reads
  `force_field[i]` under `FLAG_FF` (`wgsl.rs:273-278`). It is the shared
  accumulator that IBM / Shan-Chen / (future) surface-tension write into. The
  *field itself* is resident; the pre-passes that fill it run host-side (see
  T4/T6/T7).
- **(b) After W-VOF:** unchanged as the accumulator. W-VOF adds *contributors*
  (T4/T5), not a new consumption path.
- **(c) GPU:** supported today (`gpu/backend.rs:1602-1642`, `FLAG_FORCE_FIELD`).
- **(d) MPI/halo:** the *field read* is cell-local (no halo). Contributors that
  compute it from a gradient stencil own their own halo (T4/T6).
- **(e) Pass position:** read in `collide`, before the half-force correction.
  Second summand.

### T3. `rho * g` — per-mass gravity (single-phase today)
- **(a) Today:** **backend-side on all three backends**, with a host-staged
  fallback retained. `StepParams::gravity: Option<[T;3]>` (`params.rs:159`);
  `force_at` composes `field[idx][a] + rho*g[a]` grouped so
  `F_uniform + (F_cell + rho*g)` (`params.rs:234-248`). GPU composes
  `fv += rho * P.g` under `FLAG_GRAVITY` (`wgsl.rs:282-286`), gated by
  `supports_gravity_body_force()==true` for `CpuScalar` (`backend.rs:396-398`),
  `CpuSimd` (`backend_simd.rs:1496-1498`), and `WgpuBackend`
  (`gpu/backend.rs:1904-1906`). The solver takes the backend path when
  `gravity.is_some() && backend.supports_gravity_body_force()`
  (`solver.rs:1690-1700, 1721-1734, 3238-3241`). If a backend returns `false`,
  `run_staged_step` (`solver.rs:1651-1671`) falls back to the host-staged
  overlay `stage_gravity`/`unstage_gravity` (`solver.rs:1591-1649`), which
  writes `rho*g` into `force_field`, uploads, runs, then unwinds — bit-identical
  to a raw per-cell field (frozen contract, `solver.rs:1588-1590`).
- **(b) After W-VOF:** backend-side, with the density factor replaced. Per
  WVOF_IMPL_SPEC D5 and the `stage_gravity` doc-comment
  (`solver.rs:1582-1590`), the composed term becomes `(rho(phi) - rho_ref)*g`
  with `rho_ref` = reference (default `rho_g`, lighter phase); `rho_ref = 0`
  keeps single-phase bit-identical. **The composition point does not move**:
  W-VOF supplies `rho(phi)` into the same slot (`force_at` on CPU, the
  `rho * P.g` line on GPU). W-VOF phase-1 lands the `rho(phi)` edit on the
  host-staged / `force_field` path (WVOF §4.4); whenever a device-side
  `rho(phi)*g` composition is written it must consume the identical `rho(phi)`
  and `rho_ref` inputs (WVOF §4.4 hand-off requirement).
- **(c) GPU:** supported today for `rho*g` (`wgsl.rs:282-286`). GPU `rho(phi)*g`
  is a follow-on gated on the `phi` field being device-resident (B-1); the CPU
  reference is the oracle first.
- **(d) MPI/halo:** none for single-phase (`rho` is cell-local). For
  `rho(phi)`, `phi` itself is halo-exchanged by the W-VOF phase-field pre-pass,
  but the gravity multiply is cell-local once `phi(x)` is known.
- **(e) Pass position:** in `collide`, before the half-force correction. Third
  summand. Note the well-balanced requirement (VR-STR-06): `F_s` and
  `(rho(phi)-rho_ref)*g` must be summed *together* into `F_total` so their
  gradient parts cancel against the pressure gradient at the same discrete
  precision — they cannot be composed in separate passes.

### T4. `F_surface_tension = mu_phi * grad(phi)` — chemical-potential form
- **(a) Today:** **not implemented** (single-phase core has no `phi`). The
  Shan-Chen SCMP/MCMP cohesion force (T6) is the closest existing analog and is
  the code path W-VOF reuses.
- **(b) After W-VOF:** host-side pre-pass into `force_field`. Adopt the
  chemical-potential form (WVOF D4): `mu_phi = 4 beta phi(phi-1)(phi-1/2) -
  kappa grad^2(phi)`, `F_s = mu_phi grad(phi)` (WVOF eq. 13-14). Composed at the
  solver level in the phase-field pre-pass (WVOF §4.2 step 1c), structurally
  identical to `update_shan_chen_force` (`solver.rs:2399-2499`): exchange one
  padded `phi` plane, run a q-stencil with halo, write `force_field`. **CSF
  `sigma kappa n_hat delta_s` is validation-only** (WVOF D4), never the
  production force.
- **(c) GPU:** deferred (WVOF D8/phase-3). Needs a `phi` gradient shader + the
  `g`-distribution buffer; must not block CPU phase-1. `Backend::Fields`
  reserves the storage extension (`backend.rs:130-135`).
- **(d) MPI/halo:** **yes — needs neighbor data.** `grad(phi)` and
  `grad^2(phi)` are finite-difference stencils; the pre-pass must exchange the
  padded `phi` plane (`HaloExchange::exchange_scalar`, the `psi_planes` reuse,
  `solver.rs:2429-2440`) before evaluating the stencil, exactly as Shan-Chen
  does.
- **(e) Pass position:** solver pre-pass **before** `collide`; the value lands
  in `force_field` and is read in `collide` before the half-force correction.
  Fourth summand. Must be summed jointly with T3 (well-balanced).

### T5. `F_scalar_buoyancy` — active-scalar (thermal/species) buoyancy
- **(a) Today:** **not implemented** in the hydrodynamic force path. (An
  active-scalar-feedback proposal exists at
  `docs/proposals/active-scalar-feedback.md`; the buoyancy coupling into
  `force_field` is not landed.)
- **(b) After W-VOF / W-SCAL:** host-side pre-pass into `force_field`, e.g.
  Boussinesq `F_b = rho * beta_s (s - s_ref) g_dir`, composed in the same
  accumulation as T3/T4 (WVOF §4.4 lists `F_b^scalar` explicitly). Any concrete
  closure lands with its own PHYSICS.md provenance entry + validation test —
  this spec only reserves its slot and composition point.
- **(c) GPU:** deferred with T4 (needs the scalar field device-resident).
- **(d) MPI/halo:** the buoyancy multiply is cell-local; the *scalar transport*
  that produces `s(x)` owns its own halo. The force term itself needs no
  neighbor data.
- **(e) Pass position:** solver pre-pass before `collide`; read in `collide`
  before the half-force correction. Fifth summand.

### T6. `F_shan_chen` (the existing multiphase cohesion, folded under T2)
- **(a) Today:** host-side pre-pass, **overwrite semantics.**
  `update_shan_chen_force_with_walls` (`solver.rs:2399-2499`) exchanges the
  padded `psi` plane, runs the q-stencil, and **writes** `ff[c] = -(psi_i(g s +
  ...))` (`solver.rs:2488-2492`) — it also zeroes non-participating cells
  `ff[c] = [0;3]` (`solver.rs:2457`). It does **not** accumulate.
- **(b) After W-VOF:** unchanged for the SCMP/MCMP path; W-VOF's surface-tension
  pre-pass (T4) is a *separate* AC-model contributor and must **accumulate**,
  not overwrite (see §2, ordering rule R3).
- **(c) GPU:** not device-resident; host pre-pass only.
- **(d) MPI/halo:** **yes.** `exchange_scalar` on the `psi` plane
  (`solver.rs:2429-2440`); solid neighbors feed the cohesion sum via the
  exchanged mask (`solver.rs:2470-2479`).
- **(e) Pass position:** solver pre-pass before `step`; read in `collide`.

### T7. `F_rotating_ibm` — rotating rigid-body direct-forcing IBM
- **(a) Today:** host-side pre-pass, **accumulate semantics.**
  `apply_rotating_ibm` (`solver.rs:2113-2311`) interpolates fluid velocity to
  markers, computes the direct-forcing increment, spreads it to the Eulerian
  grid, and **accumulates** into `force_field`:
  `ff[c][a] = ff[c][a] + add[a]` (`solver.rs:2301-2304`). The volume-penalization
  path likewise feeds `force_field`. Both enter collision only through Guo
  forcing (PHYSICS.md §"Rotating bodies").
- **(b) After W-VOF:** unchanged, host-side pre-pass, accumulate.
- **(c) GPU:** not device-resident (marker interpolation/spreading is host-side);
  the accumulated `force_field` is then uploaded like any per-cell field.
- **(d) MPI/halo:** the marker interpolation/spreading stencil
  (`marker_stencil`, `rotating_ibm.rs:163`) reads/writes a neighborhood of
  Eulerian cells around each marker. In a decomposed run a marker near a seam
  needs cells on the neighbor part; the current path gathers to a global buffer
  (`solver.rs:2131-2184`, `gather_rho`) — a full distributed IBM spread is
  out of scope here but the halo dependency is **yes, needs neighbor data.**
- **(e) Pass position:** solver pre-pass before `step`; read in `collide`.

### T8. `F_particle_reaction` — two-way particle→fluid reaction (RESERVED)
- **(a) Today:** **not implemented.** `ParticleSet` is one-way only:
  `particles.rs:5-7` header — "It implements only the one-way FR-PART-01 …
  drag, but they do not apply any reaction force back to" the fluid. No
  `force_field` write exists in `particles.rs`.
- **(b) After the two-way target lands:** host-side pre-pass into `force_field`,
  **accumulate** (Newton's third law: `-sum_p F_drag,p` spread to the grid),
  composed in the same accumulation as T7. Reserved slot; the concrete closure
  lands with its own provenance + test.
- **(c) GPU:** deferred.
- **(d) MPI/halo:** **yes** — same spreading-stencil seam dependency as T7.
- **(e) Pass position:** solver pre-pass before `step`; read in `collide`.

---

## 2. Composition-order invariants

The step-invariant order is fixed (CLAUDE.md; verified `backend.rs:258-320`
`run_span`): `collide → exchange_f (halo) → stream → apply_bouzidi → swap →
apply_open_faces → apply_volume_sources → update_moments`. Force composition
lives entirely inside `collide` (the read of `F_total`) and in the solver
pre-passes that run before `run_span`.

**R1 — Full sum before the half-force.** The complete `F_total` (all summands
T1–T8 present in a given configuration) must be formed *before* the Guo
half-force velocity correction `u = (m + F_total/2)/rho`. Splitting the
half-force across partial sums is forbidden: `(m + F_a/2)/rho` then adding
`F_b/2` later changes the equilibrium `u` seen by the *rest* of the collision
and breaks Galilean consistency. `force_at` (`params.rs:232-251`) and the WGSL
preamble (`wgsl.rs:269-305`) already satisfy this; any new backend-side term
must be added to the `fv` accumulator *above* the `u = (m + 0.5*fv)*inv` line
(`wgsl.rs:303-305`), never after.

**R2 — Fixed summation order (bit gates).** Floating-point addition is
non-associative; T13 (partition invariance) and T14 (backend equivalence) gate
on bit/threshold parity of the composed result. The frozen order is the one in
`force_at`: `F_uniform + (F_cell_field + gravity)`, i.e. the per-cell terms are
summed into `force_field` first (in pre-pass write order T4→T5→T7→T8, all
accumulating), then gravity is added to that per-cell sum, then the uniform is
added outermost. The grouping parentheses in `params.rs:234-248` are
**load-bearing** and must be reproduced identically on every backend (the GPU
does `fv = P.f; if FF fv += ff; if GRAV fv += rho*g` — same result because
uniform+cell are exact-order and gravity is added last to the running `fv`).
Rule the implementation must follow: **a term's contribution to
`force_field` is committed by its pre-pass in a fixed, documented order; the
backend then reads `force_field` and adds uniform + gravity in the exact
`force_at` grouping.** Any reordering must be justified against T13/T14 and
land a new frozen recipe.

**R3 — Accumulate vs overwrite.** Contributors that share `force_field` must
**accumulate** (`ff[c] += ...`), matching IBM (`solver.rs:2303`). The one
exception is Shan-Chen SCMP/MCMP, which **overwrites** and zeroes
(`solver.rs:2457,2488`) because it is the *sole* multiphase force for that
model. W-VOF surface tension (T4) and scalar buoyancy (T5) are **additive**
contributors and must accumulate; they must never be layered on top of an
overwriting Shan-Chen pass in the same configuration (the two multiphase models
are mutually exclusive). A pre-pass that overwrites when another contributor is
active is a composition bug.

**R4 — Commutativity assumption.** The per-cell additive terms (T2/T4/T5/T7/T8)
are mathematically commutative, but *not* bit-commutative. The frozen write
order (R2) is the contract; do not rely on commutativity to reorder pre-passes.

---

## 3. Capability-negotiation rule

A backend advertises whether it can compose a term device-side. Today the only
such capability is `Backend::supports_gravity_body_force()`
(`backend.rs:189-191`; `true` for both CPU backends and GPU). The rule, frozen:

> A solver configuration that requires term `X` on a backend that cannot
> compose `X` device-side must **either** (i) fall back to an **explicit**
> host-staging path with a recorded performance caveat, **or** (ii) reject the
> configuration with a `SpecError`. It must **never silently change the
> composition point per step**, and never silently drop the term.

- **Precedent (gravity):** the solver checks
  `gravity.is_some() && supports_gravity_body_force()` once
  (`solver.rs:1690,1721,3238`) and takes *one* path for the whole run — backend
  composition when capable, else `run_staged_step` host overlay
  (`solver.rs:1651-1671`). The choice is stable across steps, not re-decided
  per step.
- **Extension for future terms:** each new backend-composable term gets its own
  `supports_*` trait method with default `false` (so a new backend is
  fail-safe: it does not claim a capability it lacks). W-VOF surface
  tension/gravity-`rho(phi)` follow the *host pre-pass* path on phase-1 (always
  available, backend-agnostic), so no negotiation is needed until a device-side
  `rho(phi)*g` shader lands; when it does, it adds `supports_phase_gravity()`
  and the solver picks one path per run, exactly like gravity.
- **The rejection arm** uses `SpecError::UnsupportedOnGpu { feature }`
  (`solver.rs:301-305`) — already the pattern for localized sources/face
  patches that GPU cannot do.

Forbidden: a backend that composes a term device-side on some steps and via
host staging on others within a single run; a `supports_*` method that returns
`true` but silently no-ops; a fallback that stages without recording the
caveat.

---

## 4. Conformance-test plan

For **each composable term**, two gates.

**(A) Host-composed vs backend-composed equivalence.** For any term a backend
can compose device-side, a test must run the identical scenario twice — once
forcing the host-staged / pre-pass path, once the backend path — and assert
bit/threshold parity of the velocity field. Precedent, gravity:
`gravity.rs::gravity_channel_is_bit_identical_to_raw_rho_g_force_field`
(`tests/gravity.rs:270`) proves the backend `rho*g` equals a raw per-cell
`force_field` filled with `rho*g`; and
`t14_backend_equiv::t14_gravity_body_force_device_resident`
(`tests/t14_backend_equiv.rs:460`, `#[ignore]`, BENCH-PENDING on a native GPU
adapter) proves the GPU path matches. Required coverage:

| Term | Host-vs-backend equivalence test (required) |
|------|---------------------------------------------|
| T1 F_user_uniform | covered by existing backend-equiv (uniform force in T14) |
| T2 F_cell_field | covered by `t14_backend_equiv` force-field cases |
| T3 rho*g gravity | `gravity.rs` bit-identical + `t14_gravity_body_force_device_resident` |
| T3' (rho(phi)-rho_ref)*g | **NEW (W-VOF):** CPU-ref vs future GPU shader must reuse the same rho(phi)/rho_ref inputs; add a `t14_phase_gravity_device_resident` when the shader lands |
| T4 F_surface_tension | **NEW (W-VOF):** potential-form vs CSF cross-check (VR-STR-03) + host-pre-pass determinism |
| T5 F_scalar_buoyancy | **NEW (W-SCAL):** analytic Boussinesq balance test |
| T6 F_shan_chen | existing multiphase validation (T11c and buoyancy-sign tests) |
| T7 F_rotating_ibm | existing IBM validation |
| T8 F_particle_reaction | **NEW (two-way):** momentum-conservation (sum of reaction = -sum of drag) |

**(B) T13/T14 coverage (mandatory for any composition change).** Every change
to how a term is composed must keep green, in `--release`:
- `tests/backend_simd_equiv.rs` — CpuScalar vs CpuSimd bit/threshold parity.
- `tests/t13_split_invariance.rs` + `tests/t13_adversarial.rs` — partition
  invariance (the composed force must be identical under any decomposition;
  this is why per-cell terms are compact-core and gravity is cell-local).
- `tests/t14_backend_equiv.rs` + `tests/t14_3d_backend_equiv.rs` +
  `tests/t14_adversarial.rs` — backend equivalence (fused vs reference).
- GPU build gate `cargo build -p lbm-core --release --features gpu`; runtime GPU
  parity via the `#[ignore]` T14 device-resident tests on a GPU host.

A new backend-composable term is not "landed" until its (A) equivalence test
exists **and** the full (B) T13/T14 suite is green with the term active.

---

## 5. Definition of regression

A term is marked **backend-side** in §0/§1 (T1, T2, T3, and any future
`supports_*`-gated term). The following are regressions and **must fail a named
test**:

1. **Re-introducing per-step host staging for a backend-side term.** Any future
   PR that makes a term marked backend-side round-trip through host staging on
   *every step* when the backend advertises the capability is a regression. It
   must fail `tests/t14_backend_equiv.rs`
   (`t14_gravity_body_force_device_resident` and its successors): the device
   path would no longer be exercised, or would diverge from the host path, and
   the equivalence assertion breaks. Additionally, a performance-guard test on
   the gravity path (staging is O(cells) host work per step) must flag the
   reappearance of a per-step `stage_out_all → stage_gravity → stage_in` when
   `supports_gravity_body_force()` is `true`.

2. **Moving the composition point per step.** A PR that composes a term
   device-side on some steps and host-side on others within one run violates §3
   and must fail T13/T14 (the composed `force_field` would differ step-to-step
   under the same decomposition/backend).

3. **Splitting the half-force.** A PR that applies the Guo half-force before the
   full `F_total` is summed (violates R1) must fail `tests/backend_simd_equiv.rs`
   and the T14 backend-equivalence suite (the recovered velocity moment changes).

4. **Reordering the frozen summation without a new recipe.** A PR that changes
   the `force_at` grouping (`F_uniform + (F_cell + gravity)`) or the pre-pass
   write order without landing a new frozen T13/T14 recipe must fail the
   bit-parity assertions in `tests/t13_split_invariance.rs` /
   `tests/backend_simd_equiv.rs`.

5. **Overwriting a shared accumulator.** A PR that makes an additive contributor
   (T4/T5/T7/T8) overwrite `force_field` instead of accumulating (violates R3)
   must fail the multi-source composition test (e.g. gravity + IBM active
   together): the earlier contributor's force vanishes and the momentum budget
   breaks.

If any of tests 1–5 does not yet exist for a term, writing it is part of that
term's landing order — a backend-composable term ships with its regression
guard or it does not ship.
