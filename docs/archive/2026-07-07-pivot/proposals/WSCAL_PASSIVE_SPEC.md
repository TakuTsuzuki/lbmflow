# W-SCAL Phase 1 Implementation Specification — Passive Scalar Advection–Diffusion (ADE)

**Document ID**: SPEC-W-SCAL-PASSIVE (rev.1, 2026-07-07).
**Scope**: the M-F item `W-SCAL passive scalar ADE` of
`docs/REQ_STIRRED_REACTOR.md` (§11 DAG; `⊂ MF-ε`) — phase 1, the **passive**
scalar transport distribution `h`. Passive = the scalar is advected and
diffused by the *resolved* velocity field `u` with **NO feedback** to density,
viscosity, or surface tension. Active feedback (`ρ(C)`, `μ(C)`, `σ(C)`,
`F_b^scalar`) is a later phase (REQ §1 "Scalar" row: `active` fidelity default
is PENDING; `passive` relaxation is the phase-1 delivery here).
**Target core**: `crates/lbm-core` (D3Q7 ADE distribution `h`, carried
alongside the D3Q19/D3Q27 hydrodynamic `f` and — when W-VOF lands — the D3Q19
phase-field `g`).
**Acceptance**: VALIDATION.md **T17** row **VR-STR-04** (scalar/reaction:
Taylor–Aris dispersion is the named MF-ε gate) plus VR-STR-05 (scalar
total-mass conservation) and the REQ §8 mandatory "phase-wise scalar
total-mass conservation" consistency test. Provisional MVP bands in §5.

This spec is **executable**: every literature choice is decided and justified,
every code touchpoint is cited against the current worktree, and every gate is
mapped to a T17 row with a provisional band. A follow-on codex implementation
order should not need to re-derive any design decision here.

> **Discipline note (CLAUDE.md prime directive; PHYSICS.md;
> `.claude/skills/lbmflow-physics-discipline`).** Every closure below is
> resolved from the governing ADE or a literature-backed closure with a
> recorded derivation, validity domain, and a dedicated validation test. No
> band-calibrated constant, no case-keyed branch, no transport-absorbing clamp
> appears anywhere in this design (the only clamp discussed — negative-`C`
> limiting — is explicitly banned in phase 1 as it would mask a transport bug;
> see §6.5). The mandatory PHYSICS.md entry text is in §7.4; the
> behavior-validity review checklist is in §5.6.

---

## 0. Summary of decisions (read this first)

| # | Decision | Justification (short) |
|---|---|---|
| P1 | **BGK-relaxed ADE distribution `h` with a TRT option**; equilibrium **linear in velocity** (advection-only second moment not required for the passive ADE). Diffusivity mapping `τ_s = D/cs_s² + 1/2`. | Standard ADE-LBM (Krüger et al. 2017 §8.3; Chopard, Falcone & Latt 2009). BGK is the reference; TRT (Ginzburg 2005) is provided because a well-chosen magic parameter `Λ = 1/4` removes the BGK numerical-diffusion anisotropy and makes the wall closure exactly located — the accurate-diffusion path the REQ's Taylor–Aris gate needs. |
| P2 | **Scalar lattice = D3Q7** (new `Lattice` impl), NOT D3Q19. | D3Q7 is sufficient for isotropic ADE at O(Ma²) (the diffusion tensor of D3Q7 is isotropic; Suga et al. 2015; Krüger §8.3.3), at **2.7× less memory per component** than D3Q19. The REQ NFR-01 budget row is written for **D3Q7 × 2 × f32 = 56 B/cell** — this spec matches the budget the REQ already assumed. Adding one `Lattice` impl is ~40 lines; reusing D3Q19 would triple the per-component cost the REQ explicitly budgeted for `O(10⁸–10⁹)` cells with potentially many scalar components. |
| P3 | **`h` carried in a new `SoaFields` slot** (`h: Option<Vec<T>>`, `htmp`, `conc: Option<Vec<T>>`), all `Option` so `None` is bit-identical to the single-phase / scalar-free path (the B-6 invariance discipline). | `Backend::Fields` explicitly reserves the `h` extension (`backend.rs:130-135`, verbatim names `g`, `h`). Mirrors the W-VOF `g` slot structure so the two mount without structural conflict (§8). |
| P4 | **Scalar collide/stream is a solver-orchestration-level sub-step AFTER the hydrodynamic `f` step** (reads the just-updated `u`), analogous to where `update_shan_chen_force` and WALE's `set_omega_field` sit. The invariant `f` pass order and the `Backend` trait are untouched. | REQ §5 FR-COUP-01 passive dataflow: `… fused collide-stream-moments → boundary → scalar ADE → reaction (split)`. The scalar reads the resolved `u`; it must run after `update_moments`. |
| P5 | **Advection velocity is the resolved hydrodynamic `u` read from `SoaFields::{ux,uy,uz}`** (which already carry the Guo `F/2` correction — CLAUDE.md invariant "`sim.ux()` return physical velocity"). No separate velocity storage. | The physical (F/2-corrected) velocity is the transport velocity for the scalar; using the raw moment would inject a spurious drift. |
| P6 | **Phase-1 BC set**: zero-flux (bounce-back) walls, fixed-concentration (anti-bounce-back) Dirichlet inlet, zero-gradient (Neumann) outflow. Reactive/adsorption walls (FR-BC-04), SGS scalar flux (FR-LES-04 `Sc_t`), active feedback + `F_b^scalar`, and interfacial `S^if` are all OUT of phase 1 (hooks noted). | REQ §4.7 FR-BC-04 lists no-flux / adsorption / reactive as the wall menu; phase 1 lands no-flux + Dirichlet + zero-gradient, the minimal set the Taylor–Aris + conservation gates need. |
| P7 | **`h` halo reuses `exchange_f_generic::<D3Q7, T>`** (`halo.rs:308`); the scalar plane needs no separate scalar-plane exchange (unlike the phase-field `∇φ` stencil, the ADE stencil is the streaming step itself). | The `h` distribution is a full LBE population set, so its neighbour transfer is the ordinary population halo exchange, parameterized on `L = D3Q7`. |
| P8 | **CPU-first (CpuScalar reference → CpuSimd fused), GPU deferred** behind the same `Backend::Fields` extension (`backend.rs:130-135`). | B-1 (generic staged multi-set GPU upload) is only PARTIALLY RESOLVED; forcing GPU into phase 1 would stall on B-1. Identical staging posture to W-VOF D8. |
| P9 | **No force composition; the `F_b^scalar` slot in FORCE_COMPOSITION_SPEC T5 is reserved but written by NOBODY in phase 1.** | Passive ⇒ zero feedback. The hook is documented (§4.4) so the active phase can add the Boussinesq contributor into `force_field` without a structural change. |
| P10 | **Multi-component ready but phase 1 lands ONE component.** `h`/`conc` are single-component `Vec<T>`; a `Vec<ScalarField>` generalization is a phase-2 refactor, API-reserved but not built. | Minimal scope (CLAUDE.md). The Taylor–Aris and conservation gates are single-component; nothing in phase 1 needs a second scalar. |

---

## 1. Governing equation + LBE

### 1.1 The passive advection–diffusion equation

For a single passive scalar concentration `C(x,t)` transported by the resolved
incompressible velocity `u(x,t)` with molecular diffusivity `D` (constant,
uniform — no dependence on `C` in phase 1):

```
∂C/∂t + ∇·(C u) = ∇·(D ∇C)                                              (1)
```

With `∇·u = 0` (low-Mach LBM incompressibility) this is the equivalent
non-conservative form `∂C/∂t + u·∇C = D ∇²C`. REQ §3 gives the phase-1 target
verbatim (single-phase passive, `ρ`, `α` uniform), *with the SGS term and
reaction dropped for passive phase 1*:

```
∂C/∂t + u·∇C = ∇·[(D + ν_t/Sc_t) ∇C] + R(C) + Ṡ^if      ← REQ full form
∂C/∂t + u·∇C = D ∇²C                                     ← W-SCAL phase 1     (2)
```

The dropped terms are **hooks, not omissions** (§4.5): `ν_t/Sc_t` is FR-LES-04
(SGS scalar flux; OUT of phase 1 per REQ §11 "SGS part waits on W-LES"; W-LES
is landed so this is the *first* phase-2 add), `R(C)` is W-REACT, `Ṡ^if` is
interfacial mass transfer (waits on W-VOF).

### 1.2 ADE-LBM: the `h` distribution (decision P1, P2)

Carry a **D3Q7** distribution `h_i(x,t)`, `i = 0..6`, relaxed toward an
equilibrium that is **linear in velocity** (the passive ADE needs only the
first velocity moment; the O(u²) terms of the full hydrodynamic equilibrium are
not required and are omitted, which is the standard ADE-LBM equilibrium —
Krüger et al. 2017 *The Lattice Boltzmann Method* §8.3.2, Eq. 8.20; Chopard,
Falcone & Latt 2009):

```
h_i^eq = w_i^s C [ 1 + (c_i · u)/cs_s² ]                                  (3)
```

`w_i^s` and `cs_s²` are the D3Q7 quadrature weights and lattice sound speed
(see §1.4). `C = Σ_i h_i` is the macroscopic concentration (zeroth moment); the
equilibrium's zeroth moment is `C` and its first moment is `C u`, so the
recovered macroscopic equation via Chapman–Enskog is exactly (1).

**BGK collision (the reference path):**

```
h_i(x + c_i, t+1) = h_i(x,t) − ω_s [ h_i(x,t) − h_i^eq(x,t) ]            (4)
ω_s = 1/τ_s
```

**TRT collision (the accurate-diffusion path, default for the Taylor–Aris
gate):** split into symmetric (`+`) and antisymmetric (`−`) parts with two
rates `ω_s^+`, `ω_s^−`:

```
h_i^± = ½( h_i ± h_{ī} ),   h_i^{eq,±} = ½( h_i^eq ± h_{ī}^eq )
h_i(x+c_i,t+1) = h_i − ω_s^+ (h_i^+ − h_i^{eq,+}) − ω_s^- (h_i^- − h_i^{eq,-})   (5)
```

`ω_s^+ = 1/τ_s` sets the diffusivity (see §1.3); `ω_s^-` is a free
antisymmetric rate fixed by the **magic parameter**
`Λ = (1/ω_s^+ − ½)(1/ω_s^- − ½)`. Choosing `Λ = 1/4` (Ginzburg 2005; Ginzburg,
Verhaeghe & d'Humières 2008) makes the anti-bounce-back Dirichlet wall lie
*exactly* half-way (matching the CLAUDE.md wall invariant) and removes the
BGK numerical-diffusion dependence on `τ_s`, i.e. it removes the spurious
anisotropic/`τ`-dependent diffusion that would otherwise contaminate the
Taylor–Aris effective-diffusivity measurement. This is the **decided default**
for VR-STR-04; BGK remains selectable and is the bit-exact reference the TRT
path is validated against on the pure-diffusion test (§5, V1).

### 1.3 Diffusivity → relaxation mapping (decision P1)

Chapman–Enskog on (4)/(5) with equilibrium (3) recovers (1) with diffusivity

```
D = cs_s² ( τ_s − ½ )     ⇒     τ_s = D / cs_s² + ½                       (6)
```

expressed in the **scalar lattice's** `cs_s²` (NOT the hydrodynamic `cs² = 1/3`
— D3Q7's sound speed differs; see §1.4). This mirrors the hydrodynamic
`τ = 3ν + 0.5` convention (CLAUDE.md invariant) with `ν → D` and
`3 = 1/cs²`, generalized to `1/cs_s²`. Stability / accuracy domain:
`τ_s ∈ (0.5, ~1.0]` for BGK (near 0.5 the collision is under-relaxed and the
scheme loses positivity for advection-dominated cells; above ~1 numerical
diffusion grows). Under TRT with `Λ = 1/4` the accurate window is wider because
the magic parameter decouples the wall/diffusion error from `τ_s`. The chosen
`τ_s` and `Λ` are **reported per run and frozen in PHYSICS.md** after the
characterization sweep — they are physical/closure parameters with a validity
domain, not band-fit constants.

### 1.4 D3Q7 lattice definition (decision P2 — the new `Lattice` impl)

D3Q7 velocity set (rest + 6 axis neighbours), matching the existing
`Lattice` trait shape (`lattice.rs:117-152`; `C`, `W`, `OPP`, `PAIRS`,
`FACE_UNKNOWNS`, `CS2` overridable):

```
c_0=(0,0,0)  c_1=(1,0,0)  c_2=(-1,0,0)  c_3=(0,1,0)
c_4=(0,-1,0) c_5=(0,0,1)  c_6=(0,0,-1)
```

Standard D3Q7 weights and sound speed (Krüger §8.3.3; the isotropic choice):

```
w_0 = 1 − 6λ,   w_{1..6} = λ,   with λ = 1/8  (⇒ w_0 = 1/4, w_i = 1/8)
cs_s² = Σ_i w_i c_{i,x}² = 2λ = 1/4                                       (7)
```

so **`cs_s² = 1/4`** (contrast the hydrodynamic `cs² = 1/3`), and (6) becomes
`τ_s = 4D + ½`. `λ = 1/8` is the value giving an isotropic second-order
diffusion tensor with the D3Q7 velocity set; it is a lattice-derivation
constant (fixed by requiring `Σ w_i c_iα c_iβ = cs_s² δ_αβ` with a chosen
`cs_s²`), NOT a tuned parameter. The `OPP`/`PAIRS`/`FACE_UNKNOWNS` tables are
generated by the existing `const fn` derivations (`lattice.rs:158-212`,
`opp_table`/`pairs_table`/`face_unknowns`) — the D3Q7 impl only supplies `C`,
`W`, and `CS2 = 1/4`; the closure identity `Σ_{c·n=0} w + 2 Σ_{c·n<0} w == 1`
(the trait doc-comment `lattice.rs:114-116`) holds for D3Q7 by construction and
must be asserted in a unit test.

> **Direction-ordering note.** The hydrodynamic single-source-of-truth ordering
> (CLAUDE.md core invariant) is the D2Q9/D3Q19 ordering; the D3Q7 `h` lattice is
> an **independent** distribution with its own ordering as above (rest, then
> ±x, ±y, ±z). It does not touch the hydrodynamic ordering. The `h`-set halo
> uses `exchange_f_generic::<D3Q7, T>` (§3.3) which reads `D3Q7::C`/`OPP`, so
> the two orderings never mix.

---

## 2. Boundary conditions (decision P6)

All three phase-1 BCs are applied to `h` **after `h` streaming**, on the faces
this subdomain owns that lie on a global boundary — the same structural slot as
the hydrodynamic `apply_open_faces` pass (`backend.rs:233`), but for the scalar
sub-step. The 1-cell solid rim and half-way wall placement (CLAUDE.md
invariant) apply identically to the scalar.

### 2.1 Zero-flux (no-flux / adiabatic) walls — bounce-back

A solid wall with **zero normal scalar flux** (`∂C/∂n = 0`, `D ∂C/∂n = 0`) is
imposed by **standard bounce-back** of the unknown `h_i` populations at the
wall (the same half-way bounce-back the fluid `f` uses for no-slip). Bounce-back
of the ADE distribution reflects the diffusive flux back, giving exactly the
Neumann zero-flux condition at the half-way wall location (Krüger §8.4.1;
Zhang, Bengough & al. review of ADE-LBM BCs). Closure: for each wall-adjacent
fluid cell, the incoming unknown `h_i` (direction `i` pointing from wall into
fluid) is set to the post-stream `h_{ī}` that would have streamed into the
wall. This reuses the exact bounce-back mechanics of the hydrodynamic wall
(the `apply_bouzidi`/half-way BB path, `backend.rs:229`) applied to the `h`
set; no new geometry.

### 2.2 Fixed-concentration inlet (Dirichlet) — anti-bounce-back

A **fixed concentration** `C = C_in` on an inlet face is imposed by
**anti-bounce-back** of the unknown `h_i`:

```
h_i(x_f, t+1) = − h_{ī}(x_f, t+ ) + 2 w_i^s C_in [ 1 + (c_i·u_w)²/(2cs_s⁴) − u_w·u_w/(2cs_s²) ]
```

For the passive linear-equilibrium ADE the target reduces to
`h_i = − h_{ī} + 2 w_i^s C_in` (the velocity-correction term is O(u²) and is
dropped consistent with the linear equilibrium (3)); this places the Dirichlet
value exactly at the half-way wall with second-order accuracy when combined
with the TRT `Λ = 1/4` choice (Ginzburg 2005; Zhang et al.). Anti-bounce-back
(sign flip) is the decided closure because a naive "set `h_i = h_i^eq(C_in)`"
Zou–He-style Dirichlet is only first-order and location-ambiguous for ADE.

### 2.3 Zero-gradient outflow (Neumann) — copy / convective

A **zero streamwise gradient** (`∂C/∂n = 0`) outflow is imposed by copying the
unknown `h_i` from the nearest interior neighbour along the face normal
(`h_i(x_face) = h_i(x_face − n)`), the ADE analog of the hydrodynamic
convective/outflow face (`outflow_face_selected`, `kernels.rs`; the existing
`ConvectiveOutflow` BC path `backend.rs`). This lets the scalar plume leave the
domain without reflecting. Decided over a fixed-`C` outflow because the passive
plume's downstream concentration is unknown a priori; the Neumann copy is the
non-reflecting choice validated by the Taylor–Aris channel (§5, V3) where the
outlet must not clamp the dispersing profile.

### 2.4 Deferred (hooks, not phase 1)

- **Reactive / adsorption walls (FR-BC-04)** — belong to W-REACT; the wall
  pass is the same slot, with a flux source added. Not phase 1.
- **Interfacial `S^if`** (Henry partition / `k_L a`) — waits on W-VOF.
- **SGS wall scalar flux** — waits on the FR-LES-04 `Sc_t` add (phase 2).

---

## 3. Data-structure mapping (decision P3, P7)

### 3.1 What the landed machinery already supports (verified in code)

- **`Backend::Fields` is an open composite storage boundary**, reserving the
  `h` set by name: `backend.rs:130-135` (verbatim) — *"Future multiphase/scalar
  work can add additional distribution sets (`g`, `h`), per-cell properties, and
  Lagrangian buffers to this associated type…"*.
- **`SoaFields<T>`** (`fields.rs:168-210`) holds the q-major padded layout
  `f[q*n_padded + cell]`, compact-core moments (`rho, ux, uy, uz`,
  `fields.rs:179-186`), and the `Option`-gated extension fields `force_field`
  (`:196`), `omega_field` (`:199`) — the precedent that a new distribution/field
  is added as an `Option` with `None` = bit-identical legacy path.
- The padded index formula is `cell = z·(pnx·pny) + y·pnx + x`
  (`LocalGeom::pidx`, `fields.rs:78-95`); `n_core()`/`n_padded()` give buffer
  sizes.
- **Generic population halo exchange** `exchange_f_generic::<L, T>`
  (`halo.rs:308`) is parameterized on the lattice `L`; it will exchange the
  D3Q7 `h` set when instantiated `::<D3Q7, T>`. (The `exchange_scalar` plane
  path `halo.rs:71,371` is for single scalar planes like Shan–Chen `ψ`; `h` is
  a full population set and uses `exchange_f_generic` instead — decision P7.)

### 3.2 The W-SCAL additions to `SoaFields<T>`

Add to `SoaFields<T>` (`fields.rs:168`), all `Option<…>` so `None` is
bit-identical to today's scalar-free path (B-6 invariance discipline):

```rust
/// Passive-scalar ADE distribution set (D3Q7), q-major padded planes.
/// `None` ⇒ no scalar (no allocation, bit-identical legacy path).
pub h: Option<Vec<T>>,
/// Ping-pong partner of `h`. Scalar streaming writes here, then swapped.
pub htmp: Option<Vec<T>>,
/// Macroscopic concentration C = Σ_i h_i, compact core. C ≥ 0 physically
/// (positivity is a diagnostic, NOT a clamp — see §6.5).
pub conc: Option<Vec<T>>,
```

Rationale for placement / naming: these mirror the W-VOF `g`/`gtmp`/`phi`
triplet exactly (WVOF_IMPL_SPEC §3.2), so `f` (hydrodynamic), `g` (phase field),
`h` (scalar) form a uniform pattern and the two orders (W-VOF O1, W-SCAL O1)
touch **disjoint fields** of the same struct — no structural conflict (§8).
The `conc` field is the only persisted macroscopic scalar quantity; there is no
scratch gradient field (the ADE stencil is the streaming step itself, unlike
the phase-field `∇φ` reconstruction).

### 3.3 Memory cost per cell (D3Q7 h, matches NFR-01)

Per the REQ NFR-01 budget table (§7 of REQ), one scalar component with
ping-pong ×2:

| Component | Layout | B/cell (f32) |
|---|---|---|
| Scalar `h` (D3Q7 × 2) | 7 × 2 × f32 | **56** |
| `conc` (C, compact core) | 1 × f32 | 4 |

**≈ 60 B/cell per scalar component**, exactly the REQ's "Scalar h (per
component) — D3Q7 × 2 × f32 = 56" budget row. Had we reused D3Q19 (decision P2
rejected alternative), this would be 152 + 4 B/cell — 2.6× worse and a
violation of the budget the REQ sized `O(10⁸–10⁹)` cells against. At 1e8 cells
one D3Q7 scalar is ≈ 6 GB (f32); development/validation at ≤256³ (1.7e7 cells)
is ≈ 1 GB — negligible on the M5 Max dev box.

### 3.4 Halo plan (decision P7)

The `h` set is exchanged by the ordinary population halo:
`exchange_f_generic::<D3Q7, T>(subs, /* h planes */)`. Because the solver holds
`f` and `h` in the *same* `SoaFields` per part, the scalar sub-step calls the
halo exchange on the `h` buffer between `h`-collide and `h`-stream — the same
interleave the hydrodynamic step uses (`backend.rs:272`, `exchange_f` between
`collide` and `stream`). The D3Q7 stencil transfers only the 6 axis
face-neighbours (`FACE_UNKNOWNS` has 1 unknown per face for D3Q7 — the single
axis population entering through each face), a strictly smaller halo than
D3Q19's 5-per-face. T13 partition invariance (§5, V4b) gates that the exchanged
`h` halo makes the multi-part result bit-identical to single-part.

---

## 4. Pass structure (decision P4, P5, P8, P9)

### 4.1 The invariant step order (verified `backend.rs:258-300`, `run_span`)

The landed hydrodynamic per-step order is:

```
collide → exchange_f (halo) → stream (interior, then boundary shells)
        → apply_bouzidi → swap → apply_open_faces → update_moments
```

(CLAUDE.md invariant: collide → halo → stream → open BCs → boundary moments.)
CpuSimd fuses collide+stream+moments in `step_band`. This `f` order and the
`Backend` trait are **untouched** by W-SCAL.

### 4.2 Where the scalar sub-step slots in (decision P4)

W-SCAL adds a **scalar ADE sub-step composed at the solver-orchestration
level** (the same level as `update_shan_chen_force` `solver.rs:2381` and WALE's
`set_omega_field` `solver.rs:2596`), running **after** the hydrodynamic step so
it reads the just-updated physical velocity. Per solver step (REQ §5
FR-COUP-01 passive dataflow):

```
1. HYDRODYNAMIC f STEP (unchanged run_span): collide → halo → stream →
   open BCs → update_moments.  Produces ρ, u = (ux,uy,uz) [F/2-corrected].
2. SCALAR ADE SUB-STEP (new, solver-level, after moments):
   a. h collide (4)/(5): read C = Σ h_i and the resolved u from
      SoaFields::{ux,uy,uz} (decision P5); relax toward h^eq (3) with ω_s.
   b. exchange h halo: exchange_f_generic::<D3Q7,T> on the h buffer.
   c. h stream: pull-stream h → htmp; swap h/htmp.
   d. scalar BCs (§2): bounce-back walls / anti-bounce-back Dirichlet inlet /
      zero-gradient outflow, on this part's global faces.
   e. C = Σ_i h_i  (scalar moment update into `conc`).
3. (phase 2 hooks, NOT phase 1) reaction split-step R(C); active property
   update ρ(C)/μ(C); F_b^scalar composition. — none run in passive phase 1.
```

**Ordering rationale (decided, physical):** the scalar is transported by the
*current* resolved velocity, so step 2 reads the `u` produced by step 1's
`update_moments` in the same solver step — this is a **within-step, not
lagged**, coupling for the passive one-way direction (there is no back-coupling
to lag). This is exactly the REQ §5 FR-COUP-01 passive sequence
(`… → boundary → scalar ADE → reaction (split) → …`). Contrast W-VOF's
phase-field pre-pass which runs *before* the `f` step (it feeds `ρ`, `μ`, `F_s`
into the `f` collision); W-SCAL runs *after* because passive scalar consumes
`u` but produces nothing the `f` step needs. The two sub-steps therefore occupy
**different slots in the step** and do not contend (§8).

### 4.3 CPU-first / GPU staging plan (decision P8)

- **Phase 1 (CpuScalar reference).** Implement step 2 against `SoaFields` on
  `CpuScalar` — the reference backend (`backend.rs:125`). All validation (§5)
  runs here first; CpuScalar is the bit-exact oracle.
- **Phase 2 (CpuSimd fused).** Fold the `h` collide+stream into a scalar analog
  of the fused `step_band` (`FusedScratch`, `backend_simd.rs`). Gate:
  `tests/backend_simd_equiv.rs` bit/threshold parity + T13 (CLAUDE.md
  invariant: any pass-structure change must pass these before landing).
- **Phase 3 (GPU, deferred to a separate order).** A `WgpuBackend` `h`-set
  buffer + the ADE collide/stream shaders are a follow-on gated on B-1
  (PARTIALLY RESOLVED — no generic staged multi-distribution upload). Must not
  block phase 1. `Backend::Fields` reserves the storage (`backend.rs:130-135`).

### 4.4 Force composition: the reserved `F_b^scalar` hook (decision P9)

**Passive phase 1 composes NO force.** FORCE_COMPOSITION_SPEC T5
(`F_scalar_buoyancy`) reserves the slot and its composition point (a host-side
pre-pass into `force_field`, e.g. Boussinesq `F_b = ρ β_s (C − C_ref) g_dir`,
"composed in the same accumulation as T3/T4"). W-SCAL phase 1:

1. writes **nothing** into `force_field` (the passive scalar exerts no body
   force);
2. the `conc` field it produces is exactly the `s(x)` the future active-phase
   T5 contributor will read;
3. records that when the active phase lands, the buoyancy pre-pass slots into
   the phase-field/force pre-pass ordering (FORCE_COMPOSITION_SPEC §2 R1–R4:
   accumulate, before the Guo half-force, in the frozen summation order), and
   must carry its own PHYSICS.md provenance entry + validation test
   (VR-STR-06+: `C ≡ C_0 ⇒ F_b^scalar = 0` exact-zero degeneration).

No force-path code is touched in phase 1; this bullet is the hand-off contract,
not an implementation item.

### 4.5 Interaction with LES (SGS scalar flux — OUT of phase 1, hook noted)

FR-LES-04 requires the SGS scalar flux be reflected in the ADE relaxation time
via a turbulent Schmidt number `Sc_t` (default 0.7): the effective diffusivity
becomes `D_eff = D + ν_t/Sc_t`, hence `τ_s = D_eff/cs_s² + ½`. **This is OUT of
W-SCAL phase 1** (REQ §11: "SGS part waits on W-LES"; W-LES is LANDED, so this
is the *first* phase-2 add). The hook: the scalar collide reads `τ_s` from a
per-cell scalar-relaxation field the same way the hydrodynamic collide reads
`omega_field` (`fields.rs:199`, `set_omega_field` `solver.rs:2596`, consumed by
`collide_row` `kernels.rs`). Phase 1 uses a **uniform `τ_s`** (constant `D`); a
per-cell `omega_s_field` following the `omega_field` precedent is the phase-2
add that lets WALE's `ν_t` (already produced on-device per REQ §0) feed
`ν_t/Sc_t` into `D_eff`. Phase 1 must **not** silently include `ν_t`; if LES is
active, phase 1 transports with molecular `D` only and logs that SGS scalar
flux is off (a documented model limitation, not a hidden approximation —
PHYSICS.md §7.4).

---

## 5. Validation plan mapped to T17 (decision P6)

Tests are **authored adversarially by codex/Opus from this spec**, in a test
worktree that never shares with the implementation worktree (CLAUDE.md team
convention; REQ §8). Each row = metric / reference / band / grid / steps /
backend / pass-fail. Bands are provisional MVP gates (T17 "Band governance
rev.4"): tightening always allowed, loosening requires a recorded PHYSICS.md
rationale.

| ID | Test | Metric & band | Grid / steps / backend | T17 row |
|---|---|---|---|---|
| **V1** | **1D diffusion of a Gaussian vs analytic** (order-of-accuracy). Initial `C(x,0) = exp(−x²/2σ₀²)`, `u = 0`, pure diffusion. Compare to the analytic spreading Gaussian `σ²(t) = σ₀² + 2Dt`. | (a) L2 profile error vs analytic **< 1%** at `t` where `σ² = 2σ₀²`; (b) **order ≥ 1.9** (grid refinement `Δx → Δx/2` at fixed physical time, error ratio ≥ ~3.5). TRT `Λ=1/4` and BGK both run; TRT L2 ≤ BGK L2 (accuracy check). | 1D-in-3D: `256×4×4` and `512×4×4`, periodic transverse, `D` s.t. `τ_s∈[0.6,0.9]`, CpuScalar | VR-STR-04 |
| **V2** | **Advected Gaussian in uniform flow** (Galilean / translation invariance). Same blob, uniform `u=(U,0,0)`, `U=0.05`, periodic; advect one domain period. | (a) blob **centroid** returns to start within **< 0.5 Δx**; (b) profile L2 vs the pure-diffusion V1 profile at the same elapsed time **< 1%** (advection must not add spurious diffusion or distortion — Galilean check); (c) **no negative-C undershoot** beyond round-off ahead of the blob. | periodic `256×4×4`, one period, `U=0.05`, CpuScalar | VR-STR-04 |
| **V3** | **Taylor–Aris dispersion in a channel vs analytic effective diffusivity** — the **named MF-ε gate**. Point/line source of scalar in fully-developed plane-Poiseuille flow; measure the long-time streamwise variance growth rate → effective dispersion coefficient. | Compare measured `D_eff` to the **Taylor–Aris analytic** `D_eff = D (1 + Pe²/210)` for plane-Poiseuille (Aris 1956; the `/210` factor is the plane-channel coefficient — `/48` is the pipe coefficient, do NOT confuse), Pe = `Ū·H/D`. Band: **±10%** on `D_eff` (provisional MVP; tighten after characterization). Measurement window: after the initial-transient time `t > 0.5 H²/D` (cross-channel diffusion equilibrated), variance growth linear-fit R² ≥ 0.999. | 3D channel `H=64` cross, `L≥8H` streamwise, no-flux walls (§2.1), parabolic `u` frozen (or a converged Poiseuille run), Pe ∈ {10, 50, 100}, `D` s.t. `τ_s∈[0.6,0.9]`, CpuScalar | **VR-STR-04 (MF-ε, named)** |
| **V4** | **Scalar total-mass conservation to round-off with zero-flux walls.** No inlet/outlet, all walls no-flux (§2.1), arbitrary initial `C` and internal `u`. | (a) `\|Σ_cell C(t) − Σ_cell C(0)\| / Σ C(0)` **< 1e-12** (f64) / **< 1e-6** (f32) at every step to end — round-off only, NOT a band; the no-flux BB wall conserves the zeroth moment exactly. (b) **V4b partition invariance (T13):** the same run under any decomposition is bit-identical (the `h` halo, §3.4). | closed box `64×64×64`, a stirred `u` field, 20k steps, CpuScalar (+ CpuSimd for V4b) | VR-STR-05, REQ §8 "phase-wise scalar total-mass conservation" |
| **V5** | **`h=None` bit-identity DoD.** Any scenario with the scalar disabled (`h=None`, `conc=None`) produces a `probe_state_hash` **bit-identical** to the pre-W-SCAL engine on the same scenario. | exact `probe_state_hash` equality (single-backend regression, VR-STR-05 semantics) | cavity + cylinder presets, CpuScalar | VR-STR-05 (B-6 invariance) |

**Mandatory negative / consistency tests (REQ §8):**

- **Scalar total-mass (V4)** is the REQ §8 "phase-wise scalar total-mass
  conservation" gate for the passive single-phase special case (REQ §3:
  "non-conservative single-phase is a special case" of the conservative form;
  passive conservation is the zeroth-moment sum).
- **Dirichlet-vs-Neumann sign (V3/V1 negative arm):** a mutant that uses a
  first-order "set `h_i=h_i^eq(C_in)`" Dirichlet instead of anti-bounce-back
  (§2.2) must FAIL the Taylor–Aris band (proves the BC closure is load-bearing
  for accuracy).
- **No-flux vs leaky wall (V4 negative arm):** a mutant that applies a
  zero-gradient (copy) wall instead of bounce-back on the closed box must FAIL
  the mass-conservation gate (proves the no-flux closure actually conserves).
- **`ν_t` leak guard:** with LES active, the scalar must transport with
  molecular `D` only (phase 1); a test asserts `D_eff == D` (no silent SGS
  contribution) — the FR-LES-04 hook is off until phase 2.

### 5.6 Behavior-validity review (mandatory, REQ / CLAUDE.md)

After each validation run, before reporting: review the *observed* pattern, not
just the gated metric. Specifically: (a) the diffusing Gaussian stays Gaussian
and symmetric (no lattice-aligned anisotropy — the D3Q7 isotropy + TRT `Λ=1/4`
claim; a diamond/star artifact would indicate wrong weights or BGK anisotropy);
(b) the advected blob does not develop a leading/trailing asymmetry or a
grid-locked wobble (Galilean defect); (c) the Taylor–Aris variance grows
*linearly* in the long-time regime (a super/sub-linear growth means the
measurement window started before cross-channel equilibration, or the outflow
BC is reflecting); (d) `C` stays non-negative except for round-off-scale
undershoot at sharp fronts — a growing negative region signals an
advection-dominated instability (`τ_s` too near 0.5), which must be fixed by
resolving/`τ_s`, **not** by clamping (§6.5). Record the review in PHYSICS.md or
the track's findings file. A metric passing its band does **not** validate a
pattern no band covers.

---

## 6. Stability & parameter domain

### 6.1 Relaxation window

`τ_s = D/cs_s² + ½ = 4D + ½` (D3Q7, `cs_s²=1/4`). Operating band
`τ_s ∈ (0.5, ~1.0]`. Near 0.5 (`D → 0`, advection-dominated) BGK loses
positivity; TRT with `Λ=1/4` widens the usable window but does not remove the
grid-Péclet limit. Report `τ_s` and `Λ` per run.

### 6.2 Grid-Péclet (advection-diffusion CFL)

The cell Péclet `Pe_Δ = |u| Δx / D` governs positivity. For BGK, `Pe_Δ ≲ 2`
keeps the scheme monotone; above it, spurious oscillations (negative `C`)
appear at sharp fronts. TRT relaxes this bound. The unit-conversion feasibility
check (REQ §2, "diffusion-number / CFL violation") must warn on `Pe_Δ` out of
range for the chosen collision. This is a **resolution requirement**, not a
tunable — a scenario that violates it is refined or switched to TRT, never
clamped.

### 6.3 Mach / low-Mach consistency

The linear equilibrium (3) is consistent to O(Ma²); the same `Ma_lattice ≤ 0.1`
bound as the hydrodynamic field applies to the advecting `u` (REQ §2). No extra
constraint beyond the hydrodynamic one.

### 6.5 Positivity is a diagnostic, NOT a clamp (discipline)

`C ≥ 0` is physical, but phase 1 **must not** clamp negative `C` to zero.
A clamp would silently absorb transport error (a banned "position clamp /
transport-absorbing cap", CLAUDE.md prime directive;
`.claude/skills/lbmflow-physics-discipline` ban list). Negative `C` beyond
round-off is a *symptom* of a grid-Péclet / `τ_s` violation (§6.2) and must be
surfaced as a diagnostic (V2c, V5.6d) and fixed by resolution/collision, not
masked. The `conc` doc-comment (§3.2) states this.

---

## 7. Phased landing plan (decision P8, P10) + conflict surface

Three orders, file-conflict-aware. One order = one bundle = one dedicated
worktree (CLAUDE.md team convention). Implementation and adversarial-test
orders never share a worktree.

| Order | Scope | Primary files (conflict boundary) | Gate |
|---|---|---|---|
| **O1 — ADE transport core (CpuScalar)** | D3Q7 `Lattice` impl (§1.4); `h`/`htmp`/`conc` slots in `SoaFields` (§3.2); `h` BGK+TRT collide (4/5) with linear-`u` equilibrium (3); `h` halo via `exchange_f_generic::<D3Q7>`; `h` stream+swap; `C=Σh_i`; the scalar sub-step wiring at solver level (§4.2); uniform `τ_s` from `D` (6); scalar BCs (§2) bounce-back / anti-bounce-back / zero-gradient. | `lattice.rs` (add D3Q7), `fields.rs` (add h/htmp/conc), `solver.rs` (scalar sub-step orchestration), `kernels.rs` (h collide/stream/BC rows) | V1 (diffusion order), V2 (advection Galilean), V4 (mass conservation), V5 (h=None bit-identity) green on CpuScalar. |
| **O2 — Taylor–Aris + scenario/CLI plumbing** | scenario schema for a passive scalar (initial `C` field, `D`, `inlet_C`, wall type per face), config validation (grid-Péclet warning §6.2); CLI/output of `conc` field (VTI, `manifest.json`); the Taylor–Aris channel scenario. | `crates/lbm-scenario/src/lib.rs` (schema + validation), `crates/lbm-cli` (conc output) — **disjoint files from O1** | V3 (Taylor–Aris D_eff, the MF-ε gate). Depends: O1. |
| **O3 — Validation authorship (codex adversarial, separate worktree)** | All of §5 (V1–V5) + the negative/consistency tests (Dirichlet-sign, leaky-wall, ν_t-leak). Authored from THIS spec, not from the impl. | `crates/lbm-core/tests/wscal_*.rs` + `crates/lbm-scenario/tests/*` (new files only — no impl-file conflict) | Tests compile red against a stub, go green against O1/O2 as they land; freeze bands in VALIDATION.md T17 VR-STR-04/05. Runs alongside O1/O2. |

**Critical-path ordering:** O1 → O2. O3 runs concurrently from the start (test
worktree). CpuSimd fused (phase 2) and GPU (phase 3) are follow-on orders, out
of this plan's scope.

### 7.1 In-flight conflict surface with `cx/wvof-o1` (W-VOF O1)

The W-VOF O1 order (branch `cx/wvof-o1` per WVOF_IMPL_SPEC §8) implements the
D3Q19 phase-field `g` transport with the **same machinery** W-SCAL uses (a
second distribution set in `SoaFields`, a solver-level sub-step, per-cell
`Option` fields). The two orders **share three files** — `fields.rs`,
`solver.rs`, `kernels.rs` — but touch **disjoint regions** of each, by design:

| File | W-VOF O1 touches | W-SCAL O1 touches | Conflict? |
|---|---|---|---|
| `fields.rs` | adds `g`, `gtmp`, `phi` fields to `SoaFields` struct + `new()` init | adds `h`, `htmp`, `conc` fields to the same struct + `new()` init | **Adjacent, not overlapping** — both append `Option` fields to the same struct body and both add an init line in `new()`. A textual merge conflict is likely at the *struct field list* and the *`new()` initializer block*; it is a trivial both-add resolution (keep both field groups). **Merge order rule: land whichever lands first, rebase the second — both edits are additive `Option` fields with `None` default; there is no semantic conflict.** |
| `solver.rs` | phase-field **pre-pass BEFORE** the `f` step (feeds ρ/μ/F_s into collide) | scalar **sub-step AFTER** the `f` step (reads u) | **No semantic conflict** — different slots in the step (§4.2). Textual proximity only if both insert into the same `step()` orchestration method; resolve by placing the W-VOF pre-pass call before `run_span` and the W-SCAL sub-step call after `update_moments`, per each spec's §4.2. |
| `kernels.rs` | `g` collide row (D3Q19, phase-field eq + AC source) | `h` collide/stream/BC rows (D3Q7, linear-u eq) | **No conflict** — distinct functions (`collide_g_row` vs `collide_h_row`), distinct lattices. Both are new `pub(crate) fn`s appended to the module. |

**Required merge order (PM dependency queue):**
1. Neither O1 hard-depends on the other; they are **parallel** (REQ §11 DAG:
   W-SCAL and W-VOF are both `after W0`, independent edges).
2. Land order is **first-ready-first**, then rebase the second. Because both
   `SoaFields` edits are additive `None`-default `Option` fields (B-6
   invariance), the second-to-land rebase is mechanical: re-apply the struct
   field group + the `new()` init line + the module fn, no logic merge.
3. The **frozen invariant both must preserve:** with *its own* set `None`, each
   path is bit-identical to pre-feature (V5 for W-SCAL; the analogous
   `g=None`/`phi=None` DoD for W-VOF, WVOF §8 per-order DoD). This is what makes
   the two independently mergeable in either order.
4. **Shared-machinery hand-off:** if W-VOF lands first and introduces a
   `scalar_planes`/`exchange_scalar` reuse helper, W-SCAL does **not** need it
   (W-SCAL uses `exchange_f_generic::<D3Q7>` for a population set, not a scalar
   plane — decision P7); the two do not share the halo path, removing that
   contention entirely.

### 7.2 Per-order DoD (all orders)

Existing tests green *without modification*; `h=None`/`conc=None` path
bit-identical to today (`probe_state_hash` unchanged where applicable — B-6
invariance, V5); the phase-1 PHYSICS.md entry (§7.4) landed with O1;
behavior-validity review (§5.6) recorded for every validation run;
`backend_simd_equiv.rs` + T13 green (they exercise `h=None` and must stay
bit-identical after the O1 orchestration change).

### 7.4 The PHYSICS.md validity-domain statement (mandatory entry text)

The O1 order must add, on landing, a PHYSICS.md §1 stack entry and a §2
decision entry containing:

> **Passive scalar transport — ADE-LBM (Krüger et al. 2017 §8.3), D3Q7,
> BGK/TRT.** A passive concentration `C` is advected by the resolved
> (F/2-corrected) velocity `u` and diffused with constant molecular `D`, **no
> feedback** to ρ/μ/σ (active feedback = later phase). Distribution `h_i`
> (D3Q7, `cs_s²=1/4`, weights `w_0=1/4, w_{1..6}=1/8`), linear-in-`u`
> equilibrium `h_i^eq = w_i^s C[1 + c_i·u/cs_s²]`, `C=Σh_i`. Diffusivity mapping
> `τ_s = D/cs_s² + ½ = 4D + ½`. TRT with magic parameter `Λ=1/4` is the default
> (removes BGK `τ`-dependent numerical-diffusion anisotropy + places the
> anti-bounce-back Dirichlet wall exactly half-way); BGK is the reference. BCs:
> no-flux = bounce-back, fixed-C inlet = anti-bounce-back, zero-gradient outflow
> = neighbour copy. **Validity domain**: `τ_s ∈ (0.5,~1.0]`, cell-Péclet
> `Pe_Δ=|u|Δx/D ≲ 2` (BGK; wider under TRT) — outside → refine or switch
> collision, NEVER clamp `C`. `Ma_lattice ≤ 0.1`. SGS scalar flux
> (`ν_t/Sc_t`, FR-LES-04) is OFF in phase 1 (transports molecular `D` only even
> when LES is active — documented limitation, first phase-2 add). **Why here
> (not derivable from code)**: D3Q7 (not D3Q19) chosen for the ADE set to match
> the NFR-01 56 B/cell budget at 2.7× less memory (D3Q7 diffusion tensor is
> isotropic — sufficient for ADE); record the measured Taylor–Aris `D_eff` band
> and the frozen `τ_s`/`Λ`.

---

## 8. Coexistence with W-VOF (structural summary)

W-SCAL (`h`, D3Q7, post-`f` sub-step) and W-VOF (`g`, D3Q19, pre-`f` pre-pass)
are designed to **mount simultaneously without structural conflict**:

- **Storage:** `f` / `g` / `h` are three `Option` distribution sets in one
  `SoaFields`, each with a `None`-default bit-identity guarantee. Adding both is
  purely additive (§7.1).
- **Step slots:** phase-field pre-pass (before `f` collide) ≠ hydrodynamic `f`
  step ≠ scalar sub-step (after `f` moments). Three disjoint slots (§4.2).
- **Halo:** `g` reuses a scalar-plane exchange for `∇φ`; `h` uses
  `exchange_f_generic::<D3Q7>` for its population set. No shared halo buffer.
- **Force:** W-VOF writes `F_s`/`(ρ(φ)−ρ_ref)g` into `force_field`; W-SCAL
  passive writes **nothing** (P9). The active-scalar `F_b^scalar` (later phase)
  will accumulate into the same `force_field` in the FORCE_COMPOSITION_SPEC
  frozen order — that is when the two force paths first interact, and it is a
  later phase, not phase 1.
- **Active-scalar future coupling with W-VOF:** the REQ two-phase phase-wise
  conservative scalar form (`∂(α_q C_{k,q})/∂t + …`, REQ §3) reads `φ` (=`α_liq`)
  from the W-VOF `g` path. W-SCAL phase 1 is single-phase passive (uniform `ρ`,
  `α`); the phase-wise form is a phase-2 add gated on W-VOF, API-reserved by the
  single-component `conc` slot (P10) and the `φ`-availability the `g` path
  provides.

---

## 9. Load-bearing code references (grounding index)

| Claim | File:line |
|---|---|
| `Backend::Fields` reserves `g`/`h` distribution sets (verbatim) | `crates/lbm-core/src/backend.rs:130-135` |
| Invariant step order in `run_span` | `crates/lbm-core/src/backend.rs:258-300` |
| `Lattice` trait shape (C/W/OPP/PAIRS/CS2/FACE_UNKNOWNS) | `crates/lbm-core/src/lattice.rs:117-152` |
| `const fn` OPP/PAIRS/FACE_UNKNOWNS derivations (reused for D3Q7) | `crates/lbm-core/src/lattice.rs:158-212` |
| face-closure identity D3Q7 must assert | `crates/lbm-core/src/lattice.rs:114-116` |
| existing D3Q19/D3Q27 impls (no D3Q7 today — new impl needed) | `crates/lbm-core/src/lattice.rs:293-465` |
| `SoaFields` struct, `force_field`/`omega_field` `Option` precedent | `crates/lbm-core/src/fields.rs:168-210`, `:196-199` |
| q-major padded index formula | `crates/lbm-core/src/fields.rs:78-95` |
| generic population halo (parameterize `::<D3Q7,T>`) | `crates/lbm-core/src/halo.rs:308` |
| scalar-plane halo (used by Shan–Chen, NOT by `h`) | `crates/lbm-core/src/halo.rs:71`, `:371` |
| `set_omega_field` (per-cell relaxation; τ_s hook precedent) | `crates/lbm-core/src/solver.rs:2596` |
| Shan–Chen force pre-pass (solver-level sub-step precedent) | `crates/lbm-core/src/solver.rs:2381`, `:2399-2499` |
| gravity host-staged force composition (F_b^scalar slot neighbourhood) | `crates/lbm-core/src/solver.rs:1591` |
| T17 VR-STR-04 scalar/Taylor–Aris row; VR-STR-05 conservation | `docs/VALIDATION.md:348-349` |
| REQ passive dataflow (FR-COUP-01) | `docs/REQ_STIRRED_REACTOR.md` §5 |
| REQ scalar governing forms (single-phase passive special case) | `docs/REQ_STIRRED_REACTOR.md` §3, §4.2 FR-LES-04 |
| REQ NFR-01 D3Q7 56 B/cell budget row | `docs/REQ_STIRRED_REACTOR.md` §7 |
| W-VOF `g` slot / in-flight conflict surface | `docs/proposals/WVOF_IMPL_SPEC.md` §3.2, §8 |
| F_b^scalar reserved slot / composition contract | `docs/proposals/FORCE_COMPOSITION_SPEC.md` T5, §2 R1–R4 |

**Literature (decided references):**
Krüger, Kusumaatmaja, Kuzmin, Shardt, Silva & Viggen 2017, *The Lattice
Boltzmann Method* §8.3 (ADE-LBM: equilibrium, `τ_s` mapping, BCs) — **adopted
governing discretization.**
Chopard, Falcone & Latt 2009 (Eur. Phys. J. ST 171:245) — ADE-LBM
Chapman–Enskog + numerical-diffusion analysis.
Ginzburg 2005 (Adv. Water Resour. 28:1171); Ginzburg, Verhaeghe & d'Humières
2008 — TRT for ADE, magic parameter `Λ=1/4` (wall location + isotropy).
Suga, Kuwata, Takashima & Chikasue 2015 — D3Q7 ADE isotropy/accuracy.
Aris 1956 (Proc. R. Soc. A 235:67), Taylor 1953 — Taylor–Aris dispersion
`D_eff = D(1 + Pe²/210)` plane channel (the VR-STR-04 MF-ε analytic reference).
```
