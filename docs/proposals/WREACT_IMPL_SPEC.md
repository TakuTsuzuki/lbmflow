# W-REACT Implementation Specification — Homogeneous Reaction Networks (T5 tier)

**Document ID**: SPEC-W-REACT (rev.1, 2026-07-07).
**Scope**: the M-F item `W-REACT reaction / active feedback` of
`docs/REQ_STIRRED_REACTOR.md` (§0 landed-vs-pending table row "W-REACT"; §11 DAG
`after W-SCAL`) — the **T5 tier**: a *general homogeneous reaction network*
that adds a **split-step source term `R_k(C)`** to the `N`-species concentration
fields `C_k` produced by W-SCAL. Homogeneous = reactions occur in the bulk
fluid; heterogeneous (reactive/adsorption **wall**) is design-only here (§4.6,
FR-BC-04, deferred). Active feedback (`ρ(C)`, `μ(C)`, `σ(C)`, `F_b^scalar`) is a
sibling later phase and is NOT part of W-REACT T5; the reaction source is
*passive* on the flow (it moves atoms between the transported scalar fields but
does not itself write `force_field` — see §8, W-SCAL P9).
**Target core**: `crates/lbm-core` — a new `reaction` module composing on the
D3Q7 scalar `h`/`conc` machinery that W-SCAL defines.
**Dependency (hard)**: **W-REACT depends on W-SCAL multicomponent.** The single
scalar of `WSCAL_PASSIVE_SPEC` must first be generalized to `N` species
(`WSCAL_MULTICOMPONENT_SPEC`, **which does not exist yet as of this rev**). This
spec therefore **designs against the WSCAL_PASSIVE data structures generalized
per the REQ §3 multicomponent form** — every place a per-species field or loop
is assumed, the assumption is stated and flagged as a WSCAL-MC dependency (§0,
§2, §8, STOP-RULE flag S-1). If WSCAL-MC lands with a different `C_k` container
than assumed here, §2 is a mechanical re-type, not a redesign.
**Acceptance**: VALIDATION.md **T17** rows **VR-STR-04** (scalar/reaction:
"reaction-diffusion at known `Da`" is the named reaction gate) and **VR-STR-05 /
MJ-011** (scalar/element total conservation), plus the REQ §5 FR-COUP-01 /
FR-COUP-02 split-step + dt-halving convergence gate (**MJ-007**). Provisional
bands in §5.

This spec is **executable**: every kinetics closure is decided and cited with a
derivation, validity domain, units (physical vs lattice, with the explicit
`R_k` unit conversion), and its own validation test; every code touchpoint is
cited against the current worktree and the W-SCAL spec; every gate is mapped to
a T17 row with a provisional band and a behavior anchor. A follow-on codex
implementation order should not need to re-derive any design decision here.

> **Discipline note (CLAUDE.md prime directive; PHYSICS.md;
> `.claude/skills/lbmflow-physics-discipline`).** Every rate law below is a
> **component-agnostic** closure resolved from mass-action / enzyme-kinetics
> theory, referenced by **registry index only** — there is **no branch on a
> species name or a case identity** anywhere (ban-list item 1). No rate constant
> is calibrated to pass a band: `k`, `μ_max`, `K_s`, `E_a`, `A` are
> **user/scenario inputs** with a stated validity domain, reported per run and
> frozen in PHYSICS.md (ban-list item 2). **Positivity of `C_k` is preserved by
> the integrator (§4), NOT by a clamp** — a `.clamp`/`.max(0.0)` on a
> concentration is explicitly banned here because it would silently create mass
> (ban-list item 3; §4.5). The mandatory PHYSICS.md entry text is in §7.4; the
> behavior-validity review checklist is §5.6. Rule 1 artifacts for every
> kinetics term are §1.3.

---

## 0. Summary of decisions (read this first)

| # | Decision | Justification (short) |
|---|---|---|
| R1 | **Stoichiometry matrix `ν` (species×reactions) + general rate assembly `R_k = Σ_j ν_{kj} r_j`.** One reaction `j` has a rate `r_j ≥ 0` (extent velocity) and a column `ν_{·j}` of integer/rational stoichiometric coefficients (negative = reactant, positive = product). `R_k` is the net production of species `k`. | This is *the* general homogeneous network form (Aris 1965; standard chemical-reaction-engineering). It makes the network **data**, not code: any mechanism is a `(ν, {r_j})` pair. Element conservation (§5 V-CONS) is `E ν = 0` for the element-composition matrix `E` — a checkable property of the data, not of the code. |
| R2 | **Component-agnostic kinetics library, referenced by registry index.** Four decided rate-law kinds (§1.3): mass-action power law, Monod, Michaelis–Menten, Arrhenius-modulated. A reaction stores *indices* into the species registry (`reactant: SpeciesIdx`, `substrate: SpeciesIdx`, …), never a name. | CLAUDE.md ban: "branches keyed to sample/case identity". A reaction referencing species by name would be a case-identity branch; by index it is pure data. The four kinds cover batch decay, enzyme/microbial growth, and temperature dependence — the VR-STR-04 reaction menu. |
| R3 | **`R_k` is applied as an operator-split ODE sub-step AFTER the ADE transport sub-step**, per REQ §5 FR-COUP-01 (`… → scalar ADE → reaction (split) → …`). The reaction sub-step advances the local ODE `dC_k/dt = R_k(C, [T])` over `Δt` **cell-by-cell, independently** (homogeneous ⇒ no spatial coupling inside the reaction step). | Strang/Godunov operator splitting is the standard, GPU-embarrassingly-parallel treatment (LeVeque 2002; Krüger §8; the reaction is a pointwise source). Placing it after ADE means it acts on the just-transported `C_k`. |
| R4 | **Two integrator paths, selected by a computed stiffness criterion, NOT by case.** Explicit path: **SSP-RK2 (Heun)** on the local ODE, valid when the reaction Fourier/Damköhler timescale ratio `Δt/τ_rxn ≤ θ_expl` (§6.1). Stiff path: **Rosenbrock–W (ROS2, 2nd-order, L-stable)** on the local ODE with the analytic `2N`-or-smaller Jacobian, or backward-Euler as its 1st-order fallback (§4.2). | REQ §5 FR-COUP-02: "switches explicit / implicit / Rosenbrock-BDF by stiffness". The switch key is the **measured local stiffness** (a computed number), never a species/case identity. ROS2 is linearly-implicit (one Jacobian factorization per step, no Newton iteration) → deterministic, bit-reproducible, GPU-portable (FR-EXT-01). |
| R5 | **Positivity is an integrator property, never a clamp.** The explicit path is sub-cycled (§4.1) so its step respects the positivity CFL `Δt_r ≤ min_k C_k/|R_k^-|` (production-destruction split, Sandu 2001); the stiff path (backward-Euler / ROS2 on a Patankar-linearized destruction term) is **unconditionally positive** (Burchard, Deleersnijder & Meister 2003, MPRK). A negative `C_k` is a *symptom* surfaced as a diagnostic, fixed by sub-cycling / stiff-path selection, **never** masked (CLAUDE.md; W-SCAL §6.5 already bans the scalar clamp). | Ban-list item 3 (transport-absorbing clamp). Positive-preserving integrators exist in the literature, so the clamp is unnecessary AND banned. |
| R6 | **Arrhenius temperature dependence reads `T` from the thermal scalar axis when present, else an isothermal constant.** `k(T) = A exp(−E_a/(R T))`. `T` is *a species field in the same `C_k` registry* (or a dedicated thermal field the WSCAL-MC container carries); if no thermal axis is configured, `k` is the frozen isothermal `k = A exp(−E_a/(R T_ref))` with `T_ref` a reported scenario input. | No case-identity branch: "thermal present?" is a **registry query** (is a `ThermalIdx` configured), not a per-case switch. The isothermal degeneration is exact and testable (§5 negative test N-3). |
| R7 | **`Option`-gated for the B-6 bit-identity invariant.** The reaction network lives behind `SoaFields`-adjacent `reactions: Option<ReactionNetwork>` (or on the solver config). `None` ⇒ **no reaction sub-step runs, bit-identical to the pure-W-SCAL transport path** — the ABLATION guard (§5 V-ABL) asserts exactly this. | W-SCAL P3/P10 established the `Option`-`None`-is-bit-identical discipline; W-REACT extends it: `R_k = 0` (either `reactions = None`, or an empty network, or a network with all `r_j ≡ 0`) must reproduce W-SCAL's `probe_state_hash`. |
| R8 | **Units: kinetics constants are declared in PHYSICAL units, converted to LATTICE units once at scenario-build time; the ODE integrates in lattice units.** The explicit `R_k` conversion factor is stated in §1.4 (`R_k^lattice = R_k^phys · Δt_phys · (C_ref^lattice / C_ref^phys)` per species, with the first-order-rate special case `k^lattice = k^phys · Δt_phys`). NaN/±Inf detection at the contract boundary (FR-EXT-01) rejects a diverged `R_k`. | REQ FR-EXT-01: "I/O signatures with physical vs lattice units explicit; NaN/divergence detection at the contract boundary". Converting once (not per step) keeps the hot loop unit-free and deterministic. |
| R9 | **Reactive/adsorption WALL (FR-BC-04) is DESIGN-ONLY here (deferred).** The homogeneous network is the T5 deliverable; the heterogeneous surface reaction is a wall-flux source `J_wall = f(C_wall)` added to the scalar BC pass, sketched in §4.6 as a hook so it mounts later without a structural change. | Minimal scope (CLAUDE.md). The VR-STR-04 reaction gate (reaction-diffusion front) is a *homogeneous* front; no phase-1 T5 gate needs the wall reaction. |
| R10 | **CPU-first (CpuScalar reference → CpuSimd), GPU deferred**; the reaction sub-step is a pointwise per-cell kernel with no halo (homogeneous), so it is the *most* GPU-portable piece — but it lands after W-SCAL's own GPU staging (W-SCAL P8), gated on B-1. | Identical staging posture to W-SCAL P8 / W-VOF D8. The reaction kernel needs no neighbour data, so once the `C_k` fields are on-device it is a trivial elementwise pass. |

---

## 1. Governing source `R_k` + kinetics library

### 1.1 The reactive ADE and the split

W-SCAL delivers the transport of `N` passive species (REQ §3 single-phase form,
SGS + reaction dropped for passive):

```
∂C_k/∂t + u·∇C_k = ∇·(D_k ∇C_k)                              (W-SCAL, per species k)   (1)
```

W-REACT restores the reaction source, giving the reactive ADE (REQ §3
single-phase passive with `R_k`, `Ṡ^if` still deferred to W-VOF):

```
∂C_k/∂t + u·∇C_k = ∇·(D_k ∇C_k) + R_k(C[,T])                 ← W-REACT T5           (2)
```

**Operator split (decision R3, REQ §5 FR-COUP-01).** Over one solver step `Δt`
(= 1 lattice time unit), advance (2) as the composition of the W-SCAL transport
operator `𝒯` and the reaction operator `ℛ`:

```
Godunov (1st order, default): C^{n+1} = ℛ_Δt ∘ 𝒯_Δt (C^n)
Strang (2nd order, optional): C^{n+1} = 𝒯_{Δt/2} ∘ ℛ_Δt ∘ 𝒯_{Δt/2} (C^n)             (3)
```

`𝒯_Δt` is exactly the W-SCAL `h`-collide/stream/BC/moment sub-step (unchanged).
`ℛ_Δt` is the **new** sub-step: solve, *independently in every fluid cell*, the
local reaction ODE

```
dC_k/dt = R_k(C(x,·)[, T(x,·)]),   C(x,0) = C^{post-transport}(x),   over t∈[0,Δt]    (4)
```

Homogeneous ⇒ (4) has **no spatial derivative** — every cell is an independent
`N`-dimensional ODE. This is why `ℛ` needs no halo (§3.4) and is GPU-trivial
(R10). Default is Godunov (matches the FR-COUP-01 passive dataflow ordering
literally); Strang is an opt-in scenario flag for the front-speed accuracy gate
(§5 V-RDF) where 2nd-order split error matters. The split-error itself is a
gated quantity (§5 V-SPLIT, REQ FR-COUP-02 "split-error acceptance defined").

### 1.2 Stoichiometry matrix and rate assembly (decision R1)

A network of `M` reactions over `N` species is the pair `(ν, r)`:

- `ν ∈ ℤ^{N×M}` (or rational) — **stoichiometry matrix**. `ν_{kj}` = net moles
  of species `k` produced by one unit extent of reaction `j` (negative =
  consumed). Example `A + 2B → C`: column `= (−1, −2, +1)ᵀ` for `(A,B,C)`.
- `r_j(C[,T]) ≥ 0` — the **rate** (extent velocity) of reaction `j`, a kinetics
  closure from §1.3.

The net production of species `k`:

```
R_k(C) = Σ_{j=1}^{M} ν_{kj} r_j(C[,T])                                              (5)
```

**Element conservation is a property of `ν`, not of the code.** Let
`E ∈ ℤ^{P×N}` be the element-composition matrix (`E_{pk}` = atoms of element `p`
in a molecule of species `k`). A well-posed reaction conserves every element,
i.e. **`E ν = 0`** (each reaction column is in the null space of `E`). The
scenario schema (O3) *validates* `E ν = 0` at build time and rejects a network
that violates it (a mass-balance error is a scenario error, not a runtime
clamp). The runtime conservation gate (§5 V-CONS) is then the discrete check
that `Σ_cell Σ_k E_{pk} C_k` is invariant under `ℛ` to round-off, for every
element `p`.

### 1.3 Kinetics library — the four decided rate laws (Rule 1 artifacts)

All four are **component-agnostic** and reference species by **registry index**.
Each carries: form, derivation, validity domain, units, and its own validation
test (§5). The library is closed in T5 (adding a fifth kind is a later,
separately-specced order).

#### (K1) Mass-action power law

```
r = k · Π_{i∈reactants} C_i^{a_i}                                                   (6)
```

- **Form / derivation.** Elementary and empirical mass-action kinetics: the rate
  is proportional to the product of reactant concentrations raised to their
  partial orders `a_i` (for an elementary reaction `a_i = |ν_{ij}|`; general
  power-law allows fitted orders). Derives from collision theory / law of mass
  action (Guldberg–Waage 1864; Aris 1965). `a_i` are **reaction data**, not
  tuned to a band.
- **Validity domain.** Well-mixed at the cell scale (LBM cell ≫ molecular mean
  free path — always true here); dilute enough that activities ≈ concentrations.
  Order `a_i ≥ 0`. First-order (`Σa=1`) and second-order (`Σa=2`) are the
  VR-STR-04 analytic-decay special cases (§5 V-BATCH).
- **Units.** `[k]` depends on overall order `n = Σ a_i`: `[k] = (conc)^{1−n}
  time^{−1}`. First-order `[k] = s^{−1}`; second-order `[k] = m³ mol^{−1} s^{−1}`.
  Lattice conversion §1.4.
- **Test.** §5 V-BATCH (1st- and 2nd-order decay vs analytic), V-RDF (front).

#### (K2) Monod growth

```
r = μ_max · C_S / (K_s + C_S)                                                       (7)
```

- **Form / derivation.** Saturating microbial-growth rate: `μ_max` the maximum
  specific rate, `K_s` the half-saturation constant (`r = μ_max/2` at `C_S =
  K_s`), `C_S` the limiting-substrate concentration (by registry index).
  Empirical saturation kinetics (Monod 1949); mathematically identical form to
  Michaelis–Menten (K3) but named for biomass growth. Used with a biomass
  species via `ν` (e.g. `r_growth` produces biomass, consumes substrate).
- **Validity domain.** Single limiting substrate, non-inhibitory (`C_S ≥ 0`);
  `K_s > 0` strictly (avoids the `0/0` at the origin — a division-by-zero guard
  is a numerical-validity check at the contract boundary, NOT a physics clamp).
  Reduces to first-order (`r ≈ (μ_max/K_s) C_S`) when `C_S ≪ K_s` and to
  zeroth-order (`r ≈ μ_max`) when `C_S ≫ K_s` — both limits are testable (§5
  V-MONOD checks both asymptotes).
- **Units.** `[μ_max] = time^{−1}`, `[K_s] = conc`. `μ_max`/`K_s` are
  **reaction-system-dependent, reported per run, frozen in PHYSICS.md**
  (envelope §Envelope).
- **Test.** §5 V-CSTR (Monod chemostat steady state vs analytic), V-MONOD
  (batch growth curve + both asymptotes).

#### (K3) Michaelis–Menten enzyme kinetics

```
r = V_max · C_S / (K_m + C_S)                                                       (8)
```

- **Form / derivation.** Quasi-steady-state approximation for the
  enzyme–substrate complex `E + S ⇌ ES → E + P` (Michaelis & Menten 1913;
  Briggs–Haldane 1925 QSSA): `V_max = k_cat · E_tot`, `K_m = (k_{−1}+k_cat)/k_1`.
  Same functional form as Monod (7) — implemented by the *same* saturating
  evaluator with a `(V_max, K_m)` parameterization; kept as a **named alias**,
  not a duplicated branch, so there is one code path (avoids a decorative
  duplicate — ban-list item 4).
- **Validity domain.** QSSA valid when `E_tot ≪ K_m + S_0` (enzyme far below
  substrate); `K_m > 0`. Product-inhibition / reversible forms are out of T5.
- **Units.** `[V_max] = conc·time^{−1}`, `[K_m] = conc`.
- **Test.** Shares V-MONOD/V-CSTR machinery with a `(V_max,K_m)` parameter set;
  the two-parameterizations-one-evaluator identity is a unit test (V-MM-ALIAS).

#### (K4) Arrhenius temperature dependence (decision R6)

```
k(T) = A · exp( −E_a / (R T) )                                                      (9)
```

- **Form / derivation.** Temperature dependence of a rate constant (Arrhenius
  1889): `A` the pre-exponential factor, `E_a` the activation energy, `R` the
  gas constant, `T` absolute temperature. Multiplies into K1/K2/K3 by replacing
  the constant `k`/`μ_max`/`V_max` with `k(T)`. **`T` is read from the thermal
  axis when present** (a `ThermalIdx` into the species/field registry), **else
  the isothermal constant** `k(T_ref) = A exp(−E_a/(R T_ref))` with `T_ref` a
  scenario input (decision R6).
- **Validity domain.** `T > 0` K (physical); `E_a ≥ 0`; Arrhenius (not
  modified-Arrhenius `A T^b exp(...)`) — the `T^b` prefactor is out of T5. The
  isothermal degeneration must be **exact** when no thermal axis exists (§5 N-3).
- **Units.** `[A] = ` same as the rate constant it modulates; `[E_a] = J
  mol^{−1}`; `R = 8.314 J mol^{−1} K^{−1}`; `[T] = K`. `E_a`/`A` are
  **reaction-system-dependent, per-run reported, PHYSICS.md frozen** (envelope).
- **Coupling to the thermal field.** In T5 the thermal field, if present, is
  itself a W-SCAL scalar (a temperature ADE with its own `D_T = α` thermal
  diffusivity); the reaction *heat release* `Σ_j (−ΔH_j) r_j` feeding back into
  `T` is an **active-feedback item deferred to the active phase** (§8) — T5
  reads `T` but does not write reaction heat into it (that would be a `ρ/μ/T`
  feedback, outside T5's passive-on-flow scope). Documented as a limitation
  (§7.4), not a hidden term.
- **Test.** §5 V-ARRH (two-temperature rate-ratio vs `exp(−E_a/R·(1/T₁−1/T₂))`;
  isothermal degeneration N-3).

### 1.4 Units: physical ↔ lattice `R_k` conversion (decision R8, FR-EXT-01)

The ADE integrates in **lattice units** (`Δx = Δt = 1`; concentration carried as
a lattice value `C^lat`). Kinetics constants are declared **physical** and
converted **once** at scenario-build time. The source term in the lattice ODE:

```
dC_k^lat/dt^lat = R_k^lat,   with   R_k^lat = R_k^phys · Δt_phys · (C_ref^lat / C_ref^phys)   (10)
```

where `Δt_phys` is the physical time per lattice step (from the W-SCAL /
scenario unit system) and `C_ref^lat/C_ref^phys` is the concentration scale.
Special cases made explicit (so O3 validation can check them):

| Rate kind | Physical constant | Lattice constant |
|---|---|---|
| 1st-order mass-action (K1, `n=1`) | `k^phys [s⁻¹]` | `k^lat = k^phys · Δt_phys` (dimensionless per step) |
| 2nd-order mass-action (K1, `n=2`) | `k^phys [m³mol⁻¹s⁻¹]` | `k^lat = k^phys · Δt_phys · (C_ref^phys)` (scaled by conc unit) |
| Monod/MM rate (K2/K3) | `μ_max/V_max [s⁻¹ or conc·s⁻¹]`, `K_s/K_m [conc]` | `μ_max^lat = μ_max^phys·Δt_phys`; `K_s^lat = K_s^phys·(C_ref^lat/C_ref^phys)` |
| Arrhenius (K4) | `A` (units of modulated `k`), `E_a [J/mol]` | `A^lat` converted as its host `k`; the `exp(−E_a/RT)` factor is dimensionless (no conversion; `T` in K) |

**Contract-boundary NaN/divergence detection (FR-EXT-01).** After building each
`r_j` and after each `ℛ` sub-step, the reaction module checks `R_k`/`C_k` for
`NaN`/`±Inf` and rejects (returns a diagnosed error, does not silently
continue). A diverged user closure is an error at the boundary, not a default
(W-SCAL ban on silent fallbacks). The closure contract (§2.3) requires
determinism (same input → bit-identical output) and state-freedom (no interior
mutable state between cells) so the pointwise kernel is GPU-portable and the
`probe_state_hash` single-backend equivalence (FR-COUP-04) holds.

---

## 2. Reaction-network & kinetics data structures (Rust API sketch)

Designed against the W-SCAL `SoaFields` `conc`/`h` machinery **generalized to
`N` species** (the WSCAL-MC dependency, decision R7, STOP-RULE flag S-1). Where
this spec assumes a multi-species container it says so; if WSCAL-MC lands a
different `C_k` type the change is a re-type of the `species` slice below.

### 2.1 Registry-index-based network (decision R1, R2)

```rust
/// Index of a species (or the thermal field) in the W-SCAL multicomponent
/// registry. NEVER a name at runtime — names live only in the scenario schema
/// and are resolved to indices at build time (no case-identity branch, R2).
pub struct SpeciesIdx(pub u16);

/// A single reaction: a rate law + a stoichiometry column.
pub struct Reaction {
    /// Kinetics closure for r_j ≥ 0 (component-agnostic, index-based).
    pub rate: RateLaw,
    /// Sparse stoichiometry column ν_{·j}: (species, coefficient). Signed:
    /// negative = reactant, positive = product. Rational stored as f64 but
    /// element-balance (E ν = 0) is checked exactly at build time.
    pub stoich: Vec<(SpeciesIdx, f64)>,
}

/// The four decided rate laws (K1–K4), all index-based, all state-free.
pub enum RateLaw {
    /// K1: r = k · Π C_i^{a_i}. `k` may be temperature-modulated (Arrhenius).
    MassAction { k: RateConst, powers: Vec<(SpeciesIdx, f64)> },
    /// K2/K3 (one evaluator; alias distinguishes reporting only):
    ///   r = vmax · C_s / (ksat + C_s)
    Saturating { vmax: RateConst, ksat: f64, substrate: SpeciesIdx, alias: SatKind },
}

pub enum SatKind { Monod, MichaelisMenten }

/// A rate constant: either isothermal, or Arrhenius-modulated reading a thermal
/// axis by index (K4, decision R6). Values are ALREADY in lattice units (§1.4
/// conversion done at build time); `a_pre`/`e_a` too.
pub enum RateConst {
    Const(f64),
    Arrhenius { a_pre: f64, e_a_over_r: f64, thermal: Option<SpeciesIdx>, t_ref: f64 },
}

/// The whole homogeneous network. `Option<ReactionNetwork>` on the solver
/// config; `None` ⇒ no reaction sub-step, bit-identical to W-SCAL (R7, V-ABL).
pub struct ReactionNetwork {
    pub reactions: Vec<Reaction>,       // M reactions
    pub n_species: usize,               // N (must match the W-SCAL registry)
    /// Element-composition matrix E (P×N), for the runtime conservation gate
    /// and the build-time E·ν = 0 check (§1.2). Stored, not recomputed.
    pub elements: ElementMatrix,
    /// Integrator selection + stiffness threshold (§4). NOT per-case; a single
    /// solver-config value.
    pub integrator: ReactionIntegrator,
}
```

Note `Option`-gating (R7): the network hangs off the solver config as
`Option<ReactionNetwork>`; `None` skips §3 step 2f entirely.

### 2.2 The reaction sub-step interface (decision R3)

```rust
/// Advance the local reaction ODE (4) over one lattice step in EVERY fluid
/// cell, in place on the multicomponent concentration fields. No halo, no
/// neighbour reads (homogeneous). Pointwise ⇒ GPU-portable (R10, FR-EXT-01).
///
/// `conc`: the WSCAL-MC concentration container (assumed `&mut [SpeciesField]`
/// or the MC generalization of W-SCAL `conc: Option<Vec<T>>`; S-1 dependency).
pub(crate) fn react_substep(
    net: &ReactionNetwork,
    conc: &mut MultiConc,   // WSCAL-MC container (S-1)
    geom: &LocalGeom,       // fluid-cell mask + indexing (fields.rs)
    dt: f64,                // = 1 lattice step (or Δt/2 under Strang)
) -> Result<(), ReactionError>;   // Err on NaN/±Inf at the boundary (§1.4)
```

### 2.3 Closure contract (FR-EXT-01) — enforced properties

Every `RateLaw::eval(&C_local, &T_local) -> f64` must be:
1. **Deterministic / bit-identical**: same inputs → same bits (no RNG, no
   floating-nondeterministic reduction order). Enables FR-COUP-04
   `probe_state_hash` single-backend equivalence.
2. **State-free**: no mutable state carried between cells or steps (the local
   ODE state IS `C_local`). Enables the pointwise GPU kernel.
3. **GPU-evaluable**: expressible in wgsl (the four kinds are `mul`/`div`/`exp`
   only — all wgsl intrinsics).
4. **NaN-guarded at the boundary**: caller checks output (§1.4); the closure
   itself must not `unwrap_or` a physical default (ban-list item 5).

---

## 3. Mapping to the solver step (the split-step slot)

### 3.1 The invariant + W-SCAL step order (verified `backend.rs`; W-SCAL §4)

Hydrodynamic `f` step (CLAUDE.md invariant, unchanged): `collide → halo →
stream → open BCs → moments`. W-SCAL adds the scalar ADE sub-step **after**
moments (W-SCAL §4.2). W-REACT adds the reaction sub-step **after** the ADE
sub-step, exactly the FR-COUP-01 ordering `… → scalar ADE → reaction (split)`.

### 3.2 Where the reaction sub-step slots in (decision R3)

Per solver step (extending W-SCAL §4.2 step 2 with the new step 2f):

```
1. HYDRODYNAMIC f STEP (unchanged): produces ρ, u [F/2-corrected].
2. SCALAR ADE SUB-STEP (W-SCAL, per species k, UNCHANGED):
   a. h_k collide   b. exchange h_k halo   c. h_k stream+swap
   d. scalar BCs    e. C_k = Σ_i h_{k,i}   (→ post-transport C_k)
   f. REACTION SUB-STEP (NEW, W-REACT; only if reactions = Some):
      for each fluid cell (independent, no halo — decision R3):
        solve dC_k/dt = R_k(C[,T]) over Δt via the selected integrator (§4);
        write updated C_k back in place (NO clamp — §4.5).
      Under Strang (3): step 2f runs with Δt/2 BEFORE 2a and Δt/2 AFTER 2e
      (the transport is the middle half-steps). Godunov default: 2f once, after 2e.
3. (active-phase hooks, NOT T5) property update ρ(C)/μ(C)/σ(C); F_b^scalar;
   reaction-heat feedback into T. — none run in W-REACT T5 (§8).
```

**Ordering rationale (decided, physical).** The reaction acts on the
concentration field *after* it has been transported to its current location in
this step — the split (3). The reaction produces **nothing the transport or
`f` step needs within the same step** (T5 is passive on the flow), so it is the
last scalar operation; there is no back-coupling to lag. This is the literal
REQ §5 FR-COUP-01 passive sequence. (The active phase inserts step 3 and, for
strong coupling, the predictor–corrector iteration — out of T5.)

### 3.3 CPU-first / GPU staging (decision R10)

- **Phase 1 (CpuScalar reference).** `react_substep` over `SoaFields` on
  CpuScalar — the bit-exact oracle. All §5 validation runs here first.
- **Phase 2 (CpuSimd).** The pointwise ODE vectorizes trivially across cells
  (no neighbour data). Gate: `backend_simd_equiv.rs` + T13 stay bit-identical
  with `reactions = None` (must, since 2f is skipped) and within the SIMD
  threshold with reactions on.
- **Phase 3 (GPU, deferred, gated on B-1).** The reaction kernel is the *most*
  GPU-friendly pass (elementwise, no halo — R10); it lands after W-SCAL GPU.

### 3.4 No halo (decision R3)

`react_substep` reads and writes only the local cell's `C_k` (and local `T`).
It is **purely pointwise** — no `exchange_*` call, no neighbour access. T13
partition invariance is therefore automatic for the reaction sub-step (it
cannot introduce a decomposition dependence); the §5 V-CONS/T13 gate confirms
it does not.

---

## 4. Stiff-integration design (decision R4, R5; REQ FR-COUP-02)

The local ODE (4) can be **stiff**: fast reactions (large `Da`, §6) relax on a
timescale `τ_rxn ≪ Δt`, so an explicit step of size `Δt` overshoots and can go
negative. FR-COUP-02 requires the solver "switch explicit / implicit /
Rosenbrock-BDF by stiffness". Two paths, switched by a **computed stiffness
number** (R4), never by case.

### 4.1 Explicit path — sub-cycled SSP-RK2 (Heun), positivity by sub-cycling

For non-stiff cells (stiffness `< θ_stiff`, §6.1): integrate (4) with an
**SSP-RK2 (Heun / explicit trapezoid)** scheme, which is 2nd-order and
strong-stability-preserving (Gottlieb, Shu & Tadmor 2001), sub-cycled into `s`
equal micro-steps `δ = Δt/s` where `s` is chosen so the **positivity CFL** holds:

```
δ ≤ CFL_pos · min_{k: R_k<0} ( C_k / |R_k^-| )                                     (11)
```

with `R_k^- = ` the destruction (negative) part of `R_k` (production–destruction
split, Sandu 2001; `CFL_pos ≈ 0.5` for RK2), and `s = ⌈Δt/δ⌉`. This guarantees
`C_k` cannot cross zero within a micro-step — **positivity from the step size,
not a clamp** (R5). `s` is computed per cell from the local state (a number),
not chosen per case. If `s` would exceed a cap `s_max` (§6.1), the cell is
**re-routed to the stiff path** (that is the stiffness criterion firing).

### 4.2 Stiff path — Rosenbrock–W (ROS2), unconditionally positive

For stiff cells (`stiffness ≥ θ_stiff`): a **2-stage Rosenbrock–W method (ROS2)**
— linearly-implicit, 2nd-order, L-stable (Verwer, Spee, Blom & Hundsdorfer
1999) — on the local ODE with the analytic Jacobian `J_{kl} = ∂R_k/∂C_l` (the
four rate laws have closed-form derivatives; §4.3). One `N×N` linear solve per
stage (`N` is small — the network species count, ≤ tens), **no Newton
iteration** ⇒ deterministic and bit-reproducible (FR-EXT-01). Positivity is
secured by applying ROS2 to the **Patankar-linearized** destruction term (the
Modified-Patankar–Runge–Kutta construction, Burchard, Deleersnijder & Meister
2003, is unconditionally positive AND conservative); backward-Euler
(1st-order) is the fallback when the Jacobian factorization is ill-conditioned.
The stiff path takes a **single** `Δt` step (its stability does not need
sub-cycling); accuracy at large `Δt` is the split-error's concern (§5 V-SPLIT),
not the integrator's stability.

### 4.3 The analytic Jacobian (state-free, GPU-portable)

`J = ν · (∂r/∂C)`, where `∂r_j/∂C_l` is closed-form per kind:
- K1 mass-action: `∂r/∂C_l = a_l · r / C_l` (for `C_l > 0`; the `C_l→0` limit is
  handled by evaluating `k a_l Π_{i≠l} C_i^{a_i} · C_l^{a_l−1}` directly — no
  division).
- K2/K3 saturating: `∂r/∂C_S = vmax · ksat / (ksat + C_S)²`, zero w.r.t. others.
- K4 Arrhenius factor: multiplies the host law's derivative; `∂/∂T` enters only
  when `T` is itself an integrated species (active-feedback, out of T5 — in T5
  `T` is frozen over the reaction sub-step, so `∂r/∂T` is not assembled).

Because `J` is assembled from the same index-based data as `R_k`, it is
component-agnostic and state-free (no per-species-name branch).

### 4.4 The stiffness criterion (decision R4 — a computed number)

Stiffness is the ratio of the step to the fastest local reaction timescale,
estimated from the Jacobian spectral radius / diagonal dominance:

```
stiffness(x) = Δt · max_k | ∂R_k/∂C_k |  ≈  Δt · ρ(J(x))                            (12)
```

`stiffness < θ_stiff` (e.g. `θ_stiff ≈ 1`, characterize-and-freeze) → explicit
path (§4.1); else stiff path (§4.2). This is the local Damköhler-per-step
(§6.2). The threshold `θ_stiff` and `s_max` are **solver-config numbers with a
stated validity domain, reported per run and frozen in PHYSICS.md** — they are
integrator-accuracy parameters, not band-fit physics constants (the physics is
the same either path; the switch only trades cost for stability).

### 4.5 Positivity is NEVER a clamp (decision R5, discipline)

`C_k ≥ 0` is physical. W-REACT **must not** `.clamp(0.0, _)` / `.max(0.0)` a
concentration. Positivity is delivered by construction: sub-cycling (11) on the
explicit path, MPRK/backward-Euler on the stiff path — both are
literature-backed positive-preserving integrators. A negative `C_k` emerging is
a **diagnostic** (the explicit CFL was violated and the cell should have taken
the stiff path, or `Δt` is too large): surface it (§5.6 review, §5 V-POS), fix
by integrator selection, **never mask**. This is the exact CLAUDE.md /
physics-discipline ban (transport-absorbing clamp) that W-SCAL §6.5 already
applied to the scalar transport, extended to the reaction step. Element
conservation (§1.2, §5 V-CONS) is the second guard: a clamp would break `Eν=0`
conservation, so the conservation gate would catch a smuggled clamp anyway.

### 4.6 Reactive/adsorption wall (FR-BC-04) — DESIGN-ONLY hook (decision R9)

Deferred from T5; sketched so it mounts without a structural change. A
heterogeneous surface reaction adds a **wall-flux source** to the W-SCAL scalar
BC pass (W-SCAL §2): instead of zero-flux bounce-back, the wall imposes a
Robin/flux BC `D_k ∂C_k/∂n |_wall = J_{k,wall}(C_wall)`, where `J_{k,wall}` is a
surface rate law (e.g. first-order adsorption `J = −k_s C_wall`, or
Langmuir–Hinshelwood surface kinetics). Implementation slot: the same
wall-adjacent cell pass as W-SCAL §2.1 bounce-back, with the reflected
population biased by the flux (a flux-augmented anti-bounce-back). It carries its
own Rule 1 artifacts + a validation test (adsorbing-wall depletion vs analytic)
when specced. **No wall-reaction code is written in T5** — this bullet is the
hand-off contract.

---

## 5. Validation plan mapped to T17 (bands + behavior anchors + negative tests)

Tests are **authored adversarially by codex/Opus from THIS spec**, in a test
worktree that never shares with the implementation worktree (CLAUDE.md; REQ §8).
Each row = metric / reference / band (with its **denominator** stated) / grid /
backend / T17 row. Bands are provisional MVP gates (tightening always allowed;
loosening needs a recorded PHYSICS.md rationale). Every row has a **behavior
anchor** (Rule 3 layer 2), not just a scalar band.

| ID | Test | Metric & band (denominator) | Behavior anchor | Grid / backend | T17 row |
|---|---|---|---|---|---|
| **V-BATCH** | **Batch reactor vs analytic decay.** Well-mixed closed box (`u=0`, no transport), single species, K1. 1st-order `A→∅`: `C(t)=C₀e^{−kt}`. 2nd-order `2A→∅`: `C(t)=C₀/(1+2kC₀t)`. | L2rel of `C(t)` vs analytic **< 0.5%**; **denominator = C₀** (initial concentration). Order-of-convergence in `Δt`: halving `Δt` reduces split+integrator error by **≥ 3.5×** (order ≥ 1.9) — explicit RK2 path. | `C(t)` **monotone decreasing**, strictly positive at every step (no undershoot); the 2nd-order curve decays **slower at late time** than the 1st-order at matched initial slope (sign/curvature check). | `8×8×8` (or 1 cell), CpuScalar | VR-STR-04 |
| **V-MONOD** | **Monod batch growth curve + asymptotes.** Biomass `X`, substrate `S`, `dX/dt=μ_max X S/(K_s+S)`, `dS/dt=−(1/Y)dX/dt`. Compare to the analytic implicit Monod batch solution. | L2rel vs reference **< 1%** (denom = `X₀+S₀·Y` = max attainable biomass). Both asymptotes: at `S≫K_s` the initial slope matches zeroth-order `μ_max X`; at `S≪K_s` the tail matches first-order `(μ_max/K_s)XS`, each within **5%** (denom = the respective analytic slope). | `X` **monotone increasing** to the yield plateau; `S` monotone decreasing to 0; the growth-rate peak occurs while `S>K_s` then decays — extremum-location anchor. | `8×8×8`, CpuScalar | VR-STR-04 |
| **V-CSTR** | **CSTR steady-state conversion vs analytic.** Continuous stirred-tank: inflow `C_in` at rate `Q`, outflow at `Q`, volume `V`, mean residence `τ=V/Q`, 1st-order reaction. Steady conversion `X = kτ/(1+kτ)`. Realized as a well-mixed cell with W-SCAL Dirichlet inflow + zero-gradient outflow + K1 reaction, run to steady state. | steady `C_out` vs `C_in/(1+kτ)` within **±3%** (denom = `C_in`). | `X` **monotone increasing in `Da=kτ`** across a `Da` sweep (`assert X(Da_hi)>X(Da_lo)`) — the FR-COUP-01 coupling-with-transport anchor. | small tank, CpuScalar | VR-STR-04 |
| **V-RDF** | **Reaction–diffusion traveling front vs analytic (the named `Da` gate).** Fisher–KPP `∂C/∂t = D∂²C/∂x² + rC(1−C)` (logistic K1 network `A→2A` style source): the pulled front travels at the analytic **minimum speed `c* = 2√(rD)`** with a `sech²`-family profile. Measure the front position vs time → speed; fit the front profile. | measured front speed vs `c*=2√(rD)` within **±10%** (denom = `c*`); profile L2rel vs the analytic KPP traveling wave **< 5%** (denom = front amplitude). Reported `Da = rL²/D` (or `r·τ_diff`) per run (envelope). Strang split ON for this row (2nd-order split needed at front). | front is a **single monotone advancing** transition (no spurious back-front, no oscillatory wake); speed is **constant** in the long-time pulled regime (linear position-vs-time, R²≥0.999); speed **increases with √(rD)** across a small sweep (monotonic anchor). | 1D-in-3D `512×4×4`, CpuScalar | **VR-STR-04 (named reaction-diffusion `Da` gate)** |
| **V-ARRH** | **Arrhenius rate-ratio vs analytic.** Batch K1 decay at two frozen temperatures `T₁,T₂` (thermal axis present, held constant). The ratio of measured effective rate constants must equal `exp(−(E_a/R)(1/T₂−1/T₁))`. | rate-ratio vs analytic within **±2%** (denom = analytic ratio). | rate **increases with T** (Arrhenius sign: `k(T₂)>k(T₁)` for `T₂>T₁`) — sign anchor. | `8×8×8`×2 runs, CpuScalar | VR-STR-04 |
| **V-SPLIT** | **Operator-split-error dt-halving convergence (MJ-007).** Coupled transport+reaction (advected + reacting blob). Halve `Δt` (via sub-stepping the whole step) and measure the change in the converged field. | error vs the `Δt→0` Richardson-extrapolated solution converges at **order ≥ 1** (Godunov) / **≥ 2** (Strang); the halving-ratio matches the claimed order within **±0.2** (denom = the order). | the error **shrinks monotonically** under successive halving (no plateau = no hidden `Δt`-independent bias such as a clamp) — the MJ-007 anchor. | advected reacting blob, `128×64×4`, CpuScalar | VR-STR-04 / MJ-007 |
| **V-CONS** | **Element/mass conservation under reaction (MJ-011).** Closed box (no-flux walls, no in/out), a multi-species network with `Eν=0`, arbitrary internal `u` and initial `C_k`. For every element `p`: `Σ_cell Σ_k E_{pk} C_k(t)` invariant. | per-element `\|total_p(t) − total_p(0)\|/total_p(0)` **< 1e-12** (f64) / **< 1e-6** (f32) at every step — **round-off, NOT a band**; denom = `total_p(0)`. Species totals individually change (reaction moves atoms) but every **element** total is invariant. **V-CONS-b partition invariance (T13):** identical under any decomposition (the reaction step is pointwise, §3.4). | at least one species total **decreases** and another **increases** over the run (proves the network is actually reacting, not a no-op) while every element total holds flat — the "moves atoms, conserves elements" anchor. | closed `64×64×64`, stirred `u`, 20k steps, CpuScalar (+CpuSimd for -b) | **VR-STR-05 / MJ-011** |
| **V-ABL** | **ABLATION guard — `R_k=0` ⇒ pure transport, bit-identical to W-SCAL.** Run a W-SCAL scenario with (i) `reactions=None`, (ii) an empty network, (iii) a network with all `r_j≡0` (e.g. `k=0`). Each must produce `probe_state_hash` **bit-identical to the pure-W-SCAL** engine on the same scenario. | exact `probe_state_hash` equality (single-backend, FR-COUP-04). | with a *nonzero* network the hash **differs** (ablation-guard converse: the reaction step is actually load-bearing — Rule 3 ablation template). | W-SCAL Taylor–Aris + closed-box scenarios, CpuScalar | VR-STR-05 (B-6 invariance) |

**Mandatory negative / consistency tests (REQ §8):**

- **N-1 Positivity-under-stiffness (V-POS negative arm).** A stiff network run on
  the **explicit path with sub-cycling disabled** (`s=1` forced) must produce a
  negative `C_k` (proving the sub-cycling / stiff-path selection is
  load-bearing for positivity); the *correct* path (auto-selected integrator)
  must keep `C_k ≥ −round-off` on the same case. Guards R5 (no clamp masking).
- **N-2 Element-balance rejection.** A scenario network with `Eν ≠ 0` (a
  deliberately unbalanced reaction) must be **rejected at scenario-build
  validation** (O3), not run — proves the build-time `Eν=0` check (§1.2) is
  live and that conservation is not silently "fixed" at runtime.
- **N-3 Isothermal Arrhenius degeneration.** With no thermal axis configured,
  the Arrhenius law must reduce **exactly** to `k=A exp(−E_a/(R T_ref))` and
  match a plain `Const(k)` network to round-off (proves R6's registry-query
  degeneration is exact, not an approximation).
- **N-4 Split-error is not a clamp bias (V-SPLIT converse).** V-SPLIT's monotone
  error-shrink under halving must NOT plateau at a nonzero floor; a plateau
  would indicate a `Δt`-independent bias (a hidden clamp or a conservation
  leak). This is the MJ-007 anti-Goodhart anchor.
- **N-5 MM/Monod one-evaluator identity (V-MM-ALIAS).** A `Monod(μ_max,K_s)` and
  a `MichaelisMenten(V_max=μ_max, K_m=K_s)` network must produce bit-identical
  `C(t)` (proves K3 is a reporting alias of K2, not a duplicated/divergent code
  path — guards ban-list "decorative duplicate").

### 5.6 Behavior-validity review (mandatory, REQ / CLAUDE.md)

After each validation run, before reporting, execute the six-step review
(physics-discipline SKILL). Reaction-specific anchors: (a) the batch decay curve
has the **correct convexity** (1st-order exponential vs 2nd-order slower tail) —
a wrong sign/curvature means the stoichiometry column or the rate assembly (5)
is transposed; (b) the Monod growth shows the substrate-limited **inflection**,
not unbounded exponential growth (would mean the saturation denominator is
inactive); (c) the reaction-diffusion front is a **single advancing monotone
transition** at the pulled speed — an oscillatory wake or a stalled front means
the split-error/positivity path is wrong; (d) `C_k` **non-negative** except
round-off — a growing negative region is the stiffness-criterion mis-routing a
cell to the explicit path (fix by integrator selection §4.4, **not** a clamp
§4.5); (e) every **element** total flat while **species** totals move — a
drifting element total means `Eν≠0` slipped through, or a clamp is creating
mass. Every run leaves a visual artifact (concentration field PNG / `C_k(t)`
time-series / front-position plot). Record the review in PHYSICS.md or the
track findings file. A metric passing its band does NOT validate a pattern no
band covers.

---

## 6. Stability & stiffness domain

### 6.1 Explicit-path validity domain (Damköhler / reaction-vs-diffusion)

The explicit RK2 path (§4.1) is valid (positive, accurate) when the
**reaction Damköhler number per step** is bounded. Define the reaction
timescale `τ_rxn = 1/max_k|∂R_k/∂C_k|` and the two relevant Damköhler numbers:

```
Da_I  = τ_flow / τ_rxn   (reaction vs advection)
Da_II = τ_diff / τ_rxn = L²/(D τ_rxn)   (reaction vs diffusion — the front regime)
```

The explicit path per-step stiffness (12) `= Δt/τ_rxn`. The domain:

```
Δt/τ_rxn ≤ θ_stiff  (≈ 1, characterize-freeze)  →  explicit path, sub-cycled to (11)
Δt/τ_rxn > θ_stiff                               →  stiff path (§4.2, unconditional)
```

Sub-cycle count `s = ⌈(Δt/τ_rxn)/CFL_pos⌉` capped at `s_max`; beyond `s_max` the
cell routes to the stiff path (that cap IS the practical stiffness switch). All
of `θ_stiff`, `CFL_pos`, `s_max` are integrator-accuracy config numbers,
reported per run, PHYSICS.md-frozen — they do not change the recovered physics
(both paths integrate the same ODE), only cost/stability.

### 6.2 Coupling with the transport CFL

The reaction sub-step does not touch the W-SCAL cell-Péclet / diffusion CFL
(§W-SCAL 6.2) — the split (3) means transport and reaction each obey their own
stability bound. The reaction bound is (11)/(12); the transport bound is
W-SCAL's. A stiff reaction does NOT force the transport `Δt` down (that is the
whole point of the operator split + stiff reaction integrator).

### 6.3 Determinism / GPU portability envelope

The pointwise ODE integrators (RK2, ROS2, backward-Euler) use only
`+,−,×,÷,exp` and a fixed-size `N×N` solve — all bit-reproducible with a fixed
evaluation order and wgsl-portable (FR-EXT-01, R8). NaN/±Inf at the boundary
→ diagnosed error (§1.4). This keeps FR-COUP-04 `probe_state_hash` single-backend
equivalence intact with reactions on.

### 6.4 Envelope (per REQ / task)

- Monod `μ_max`/`K_s` and Arrhenius `E_a`/`A` are **reaction-system-dependent**:
  per-run reported, PHYSICS.md frozen (never a code constant).
- `k_L a` interfacial-transfer domain 1–100 h⁻¹ is a **W-VOF / active-phase**
  interfacial-source concern (`Ṡ^if`), NOT a T5 homogeneous-reaction item — noted
  as an envelope boundary; T5 does not implement `k_L a` (deferred with `Ṡ^if`,
  §8).
- Grid ≤ 256³ for validation (matches W-SCAL dev envelope).

---

## 7. Phased landing plan + CODEX ORDER BREAKDOWN

Four orders, file-conflict-aware. One order = one bundle = one dedicated
worktree (CLAUDE.md team convention). Implementation and adversarial-test orders
**never** share a worktree.

### 7.1 CODEX ORDER BREAKDOWN

| Order | Scope | Primary files (conflict boundary) | Depends | Gate (DoD) |
|---|---|---|---|---|
| **O1 — kinetics + network core (CpuScalar)** | The `reaction` module: `SpeciesIdx`/`Reaction`/`RateLaw`/`RateConst`/`ReactionNetwork`/`ElementMatrix` (§2.1); the four rate laws K1–K4 (§1.3) index-based + state-free; `R_k = Σν r` assembly (5); the analytic Jacobian (§4.3); unit conversion (§1.4); the `react_substep` pointwise interface (§2.2) wired into the solver step AFTER the W-SCAL ADE sub-step (§3.2 step 2f), `Option`-gated (R7). Explicit RK2 path only (§4.1) in O1. | **NEW file** `crates/lbm-core/src/reaction.rs` (module, all types + rate laws + assembly + RK2); `solver.rs` (add the step-2f call after the W-SCAL sub-step — a single insertion point) | **W-SCAL O1** (needs the multicomponent `C_k` container — S-1) | V-BATCH, V-MONOD, V-ARRH, V-ABL green on CpuScalar; N-3, N-5 green; `reactions=None` bit-identical to W-SCAL (`probe_state_hash`). |
| **O2 — stiff integrator** | The ROS2 Rosenbrock–W path + backward-Euler fallback (§4.2), MPRK positivity (§4.2), the stiffness criterion + auto-selection (§4.4), sub-cycling (§4.1 (11)). Extends O1's `react_substep` integrator dispatch. | `crates/lbm-core/src/reaction.rs` — **the integrator submodule only** (`reaction::integrate`); **disjoint from O1's rate-law/network code** by living in a separate `mod integrate` within the file, or a sibling `reaction/integrate.rs`. Land after O1 to avoid same-region churn. | **O1** | V-SPLIT (dt-halving order), V-POS / N-1 (positivity, stiff path); stiff-network cases from V-BATCH/V-MONOD converge; N-4 (no plateau). |
| **O3 — scenario schema for reactions** | Scenario JSON schema for a reaction network: species registry (names→indices), per-reaction `ν` column + kinetics-kind declaration + physical constants + Arrhenius params + thermal-axis reference; the **build-time `Eν=0` element-balance validation** (§1.2, N-2), the physical→lattice unit conversion at build (§1.4), stiffness-config parse; CLI/output of per-species `C_k` fields (extends W-SCAL O2 `conc` output to `N` species). | `crates/lbm-scenario/src/lib.rs` (schema + `Eν=0` + unit conversion), `crates/lbm-cli` (multi-species output) — **disjoint files from O1/O2** | **O1** (network types), **W-SCAL O2** (scalar scenario plumbing) | Round-trips a reaction scenario to a `ReactionNetwork`; N-2 (unbalanced network rejected at build); CSTR scenario (V-CSTR) builds + runs; V-RDF front scenario builds. |
| **O4 — adversarial-test authorship (separate worktree)** | All of §5 (V-BATCH, V-MONOD, V-CSTR, V-RDF, V-ARRH, V-SPLIT, V-CONS, V-ABL) + negatives N-1…N-5. Authored from **THIS spec**, not the impl. Analytic references derived independently (decay, Monod batch, CSTR conversion, Fisher–KPP `c*=2√(rD)`, Arrhenius ratio, element-null-space). | **NEW files only** `crates/lbm-core/tests/wreact_*.rs` + `crates/lbm-scenario/tests/wreact_scenario_*.rs` — **no impl-file conflict** | reads spec; goes green against O1/O2/O3 as they land | Tests compile red against a stub, green against O1→O2→O3; bands frozen into VALIDATION.md T17 VR-STR-04/05. Runs concurrently from day one (test worktree). |

**Critical-path ordering:** W-SCAL-MC → **O1 → O2**; **O3** after O1 (needs the
network types); **O4** concurrent from the start (test worktree). CpuSimd
(phase 2) and GPU (phase 3) are follow-on orders, out of this plan's scope.

### 7.2 File-conflict boundaries (impl orders)

- **O1 owns `reaction.rs` outright** (new file — no conflict with anything) and
  makes a **single-line insertion** in `solver.rs` step 2f. The only shared file
  is `solver.rs`; the insertion point is *after* the W-SCAL sub-step call, so it
  does not collide with the W-SCAL O1 `solver.rs` edit (which adds the sub-step
  *before* this line). Merge order rule: land W-SCAL first, then O1 rebases the
  one insertion after it.
- **O2 extends `reaction.rs` in a disjoint region** (the integrator submodule);
  land after O1, mechanical rebase.
- **O3 and O4 touch disjoint files** from O1/O2 (`lbm-scenario`, `lbm-cli`, and
  new `tests/` files) — no impl-file conflict.

### 7.3 Per-order DoD (all orders)

Existing tests green *without modification*; `reactions=None` path bit-identical
to W-SCAL (`probe_state_hash` unchanged — V-ABL, B-6 invariance); the PHYSICS.md
entry (§7.4) landed with O1 (and updated by O2 with the integrator entry);
behavior-validity review (§5.6) recorded for every validation run;
`backend_simd_equiv.rs` + T13 green (they exercise `reactions=None` and must stay
bit-identical after the O1 step-2f insertion); ban-list grep over the diff clean
(no `.clamp`/`.max(0.0)` on a concentration, no species-name branch, no
band-calibrated float literal in a rate expression — §4.5, R2, R8).

### 7.4 The PHYSICS.md validity-domain statement (mandatory entry text)

O1 adds a PHYSICS.md §1 stack entry + §2 decision entry:

> **Homogeneous reaction network — split-step source `R_k` on the W-SCAL scalar
> fields (T5 tier).** A network of `M` reactions over `N` species is the pair
> `(ν, r)`: stoichiometry matrix `ν` (species×reaction, signed) and per-reaction
> rates `r_j ≥ 0`; net production `R_k = Σ_j ν_{kj} r_j`. Kinetics library
> (component-agnostic, referenced by **registry index**, never species name):
> mass-action power law `r=kΠC_i^{a_i}` (Guldberg–Waage; Aris 1965); Monod
> `r=μ_max C_S/(K_s+C_S)` (Monod 1949); Michaelis–Menten `r=V_max C_S/(K_m+C_S)`
> (Briggs–Haldane QSSA 1925, same evaluator as Monod); Arrhenius
> `k(T)=A exp(−E_a/RT)` reading `T` from the thermal scalar axis when present,
> else isothermal `k(T_ref)`. Applied as an operator split (Godunov default,
> Strang optional) AFTER the W-SCAL ADE sub-step (REQ §5 FR-COUP-01): the local
> ODE `dC_k/dt=R_k` is integrated cell-by-cell (homogeneous ⇒ no halo). Integrator
> selected by a **computed stiffness** `Δt·ρ(J)` (never by case): explicit
> SSP-RK2 sub-cycled to the positivity CFL for non-stiff cells; Rosenbrock-W
> (ROS2, L-stable) / MPRK backward-Euler for stiff cells (REQ FR-COUP-02).
> **Positivity is delivered by the integrator (sub-cycling / MPRK), NEVER by a
> clamp.** Constants (`k, μ_max, K_s, E_a, A`) are scenario inputs in physical
> units, converted once to lattice units (`k^lat=k^phys·Δt_phys`, §1.4);
> NaN/±Inf checked at the closure contract boundary (FR-EXT-01). **Element
> conservation** `Eν=0` validated at scenario build; runtime `Σ_cell Σ_k E_{pk}
> C_k` invariant to round-off. **Validity domain**: explicit path
> `Δt/τ_rxn ≤ θ_stiff` (else stiff path); reaction Damköhler `Da_I=τ_flow/τ_rxn`,
> `Da_II=L²/(Dτ_rxn)` reported per run; grid ≤256³. Reactive/adsorption **wall**
> (FR-BC-04), active feedback (`ρ/μ/σ(C)`, `F_b^scalar`, reaction-heat into `T`),
> and interfacial `k_L a`/`Ṡ^if` are **OUT of T5** (documented limitations, later
> phases). **Why here (not derivable from code)**: the four kinetics kinds, the
> stiff-integrator choice (ROS2/MPRK for unconditional positivity without a
> clamp), and the stiffness threshold are literature-backed closures; record the
> measured reaction-diffusion front speed band, the frozen `θ_stiff`/`CFL_pos`/
> `s_max`, and per-run `μ_max/K_s/E_a/A/Da`.

---

## 8. Coexistence with W-SCAL / W-VOF (structural summary)

W-REACT (`R_k` split-step, pointwise, after the W-SCAL ADE sub-step) mounts on
W-SCAL and coexists with W-VOF without structural conflict:

- **Storage.** W-REACT adds **no new distribution set** — it reads and writes the
  W-SCAL(-MC) `C_k` concentration fields in place. The `ReactionNetwork` is an
  `Option` on the solver config (R7); `None` = bit-identical to W-SCAL (V-ABL).
  No `SoaFields` field is added (contrast W-SCAL's `h`/`conc`, W-VOF's `g`/`phi`)
  → **zero struct-field conflict** with either sibling order.
- **Step slots.** phase-field pre-pass (before `f`, W-VOF) ≠ `f` step ≠ scalar
  ADE sub-step (after `f`, W-SCAL) ≠ **reaction sub-step (after ADE, W-REACT)**.
  Four disjoint slots; W-REACT is strictly downstream of W-SCAL's slot (§3.2).
- **Halo.** W-REACT is **pointwise, no halo** (§3.4) — it adds no halo contention
  with W-SCAL's `exchange_f_generic::<D3Q7>` or W-VOF's scalar-plane exchange.
- **Force.** W-REACT T5 writes **nothing** into `force_field` (passive on the
  flow, like W-SCAL passive P9). The `F_b^scalar` active-scalar contributor
  (FORCE_COMPOSITION_SPEC T5) is the **active phase**, not W-REACT T5; when it
  lands it reads the same `C_k` fields W-REACT updates and composes into
  `force_field` in the frozen summation order — that is the *only* place the
  reaction and force paths interact, and it is a later phase.
- **Two-phase future coupling with W-VOF.** The REQ phase-wise conservative
  reactive form `∂(α_q C_{k,q})/∂t + … = … + α_q R_{k,q} + S_{k,q}^{if}` (REQ §3)
  applies the *same* `R_k` per phase `q` weighted by `α_q=φ` from W-VOF, plus the
  interfacial `Ṡ^if`/`k_L a` transfer between phases. W-REACT T5 is single-phase
  homogeneous; the phase-wise `α_q R_{k,q}` weighting and `k_L a` are a later
  phase gated on W-VOF, API-reachable because the `R_k = Σν r` assembly (5) is
  already per-cell and phase-agnostic (multiply by `α_q(x)` when W-VOF's `φ` is
  available — a scalar prefactor, no structural change).

---

## 9. STOP-RULE flags carried by this spec

- **S-1 (dependency, not a physics stop).** WSCAL_MULTICOMPONENT_SPEC does not
  exist yet. W-REACT O1 **cannot land** until W-SCAL provides an `N`-species
  `C_k` container. This spec designs against the WSCAL_PASSIVE structures
  generalized per REQ §3; if WSCAL-MC lands a different container, §2 is a
  re-type. **Action: WSCAL-MC must be specced+landed before W-REACT O1 is
  dispatched.** (This is a sequencing flag for the PM, not a physics stop-rule.)
- No Rule-1 physics stop-rule fires: every kinetics closure and every integrator
  is literature-backed with a derivation, validity domain, and a validation test
  (§1.3, §4); no band is reachable only via a banned pattern. Positivity is
  delivered by a positive-preserving integrator (a literature capability), so the
  banned clamp is genuinely unnecessary, not merely prohibited.

---

## 10. Load-bearing references (grounding index)

| Claim | Source |
|---|---|
| Reactive ADE single-phase passive form + `R_k` | `docs/REQ_STIRRED_REACTOR.md` §3 (Scalars/reactions block, lines 158–171) |
| FR-COUP-01 split-step placement (`… → scalar ADE → reaction (split)`) | REQ §5 FR-COUP-01 |
| FR-COUP-02 explicit/implicit/Rosenbrock-BDF by stiffness; negative-conc limiting, element conservation, split-error acceptance | REQ §5 FR-COUP-02 |
| FR-EXT-01 closure contract (deterministic/bit-identical/GPU-evaluable/state-free/NaN at boundary/units) | REQ §4.8 FR-EXT-01 |
| VR-STR-04 scalar/reaction (Taylor–Aris + reaction-diffusion at known `Da` + `k_L a`) | REQ §8 VR-STR-04; `docs/VALIDATION.md` T17 VR-STR-04 |
| VR-STR-05 conservation (scalar totals, `probe_state_hash` single-backend) | REQ §8 VR-STR-05; VALIDATION.md T17 VR-STR-05 |
| MJ-011 scalar total-mass / element conservation | REQ §8 "scalar total-mass conservation in phase-wise form (MJ-011)" |
| MJ-007 active-scalar / split dt-halving convergence | REQ §8 "active-scalar dt-halving convergence (MJ-007)"; §5 FR-COUP-01 acceptance |
| W-SCAL `C_k`/`h`/`conc` machinery, `Option`-`None`-bit-identity discipline, positivity-not-a-clamp | `docs/proposals/WSCAL_PASSIVE_SPEC.md` §3.2, §6.5, P3/P10 |
| W-VOF `g`/`phi` slot, step-slot coexistence | `docs/proposals/WVOF_IMPL_SPEC.md` §3.2, §8 |
| `F_b^scalar` reserved slot / composition contract | `docs/proposals/FORCE_COMPOSITION_SPEC.md` T5 |
| `SoaFields` `Option` field precedent; solver-level sub-step precedent | `crates/lbm-core/src/fields.rs`, `crates/lbm-core/src/solver.rs` (Shan–Chen force pre-pass, `set_omega_field`) |
| physics-discipline ban list (clamp / case-branch / calibrated constant) | `.claude/skills/lbmflow-physics-discipline/SKILL.md` Rule 2 |

**Literature (decided references):**
- Aris, R. 1965, *Introduction to the Analysis of Chemical Reactors* — stoichiometry matrix `ν`, general rate assembly `R_k=Σν_{kj}r_j`, element-null-space `Eν=0`.
- Guldberg & Waage 1864 — law of mass action (K1).
- Monod, J. 1949 (Annu. Rev. Microbiol. 3:371) — saturating microbial growth (K2).
- Michaelis & Menten 1913; Briggs & Haldane 1925 — enzyme QSSA (K3).
- Arrhenius, S. 1889 (Z. Phys. Chem. 4:226) — `k(T)=A exp(−E_a/RT)` (K4).
- LeVeque, R. 2002, *Finite Volume Methods for Hyperbolic Problems* — operator (Godunov/Strang) splitting (§1.1).
- Gottlieb, Shu & Tadmor 2001 (SIAM Rev. 43:89) — SSP Runge–Kutta (explicit path, §4.1).
- Verwer, Spee, Blom & Hundsdorfer 1999 (SIAM J. Sci. Comput. 20:1456) — Rosenbrock-W (ROS2) for stiff reaction ODEs (§4.2).
- Burchard, Deleersnijder & Meister 2003 (Appl. Numer. Math. 47:1) — Modified Patankar (MPRK): unconditionally positive AND conservative reaction integration (§4.2, §4.5).
- Sandu, A. 2001 (J. Comput. Phys. 170:589) — positive numerical integration of chemical kinetics, production–destruction split (§4.1 positivity CFL).
- Fisher 1937 / Kolmogorov–Petrovsky–Piskunov 1937 — reaction-diffusion traveling wave, pulled-front speed `c*=2√(rD)` (§5 V-RDF analytic reference).
