# T6 Thermal / Energy Axis Implementation Specification — Temperature ADE, Reaction Enthalpy, Arrhenius Coupling, Thermal Property Feedback

**Document ID**: SPEC-THERMAL-T6 (rev.1, 2026-07-07).
**Status**: draft, executable spec. Contains a **proposed REQ amendment
(FR-THERM block, §0.5) awaiting PM ratification** — flagged inline; nothing
in this document is an authoritative REQ requirement until that block is
appended to `REQ_STIRRED_REACTOR.md`.
**Scope**: the **thermal / energy axis (T6 tier)** of the reaction-engineering
goal — temperature field `T`, its transport, thermal property feedback, and
its coupling to reaction kinetics via Arrhenius `k(T)`. PM decision (2026-07-07):
the thermal axis is **upgraded from API-reserved to in-scope**, because
Arrhenius temperature dependence of reaction rates is required for the
reaction-engineering deliverable (`W-REACT`). This supersedes the "deferred
extension point + API hook only" recommendation of
`active-scalar-feedback.md` §4.2 — that recommendation is explicitly
overridden by the PM directive; the *formulas* it fixed (§4.1) are adopted here
verbatim.
**Target core**: `crates/lbm-core`. Temperature `T` is carried on the **same
D3Q7 ADE-LBM distribution + relaxation machinery** as the passive scalar `C_k`
(`WSCAL_PASSIVE_SPEC.md`) — `T` is a *special scalar with feedback targets and
a reaction-enthalpy source*, **not a new lattice type**
(`active-scalar-feedback.md` §4.2 isomorphism).
**Depends on**: `W-SCAL` (D3Q7 `h`-machinery — LANDED per this spec's
assumption that WSCAL O1 lands first; `T` reuses its `Lattice` impl, halo,
sub-step slot) and `W-REACT` (reaction rates `r_j` and enthalpies `ΔH_{r,j}`;
Arrhenius `k(T)` closes the two-way loop). The active-feedback path
(`ρ(T)`, `μ(T)`, `σ(T)`) shares the `active-scalar-feedback.md` closures.
**Acceptance**: VALIDATION.md **T17** — new rows **VR-STR-12..15** (this spec
proposes them, §5) plus the reuse of **VR-STR-08** (thermal-capillary
migration, already sketched in `active-scalar-feedback.md` §6) as the `σ(T)`
gate. Provisional bands + two-layer gates + denominators in §5.

This spec is **executable** (mirrors `WSCAL_PASSIVE_SPEC.md` format): every
closure is decided and cited, every code touchpoint references the current
worktree or the WSCAL/WREACT contract, every gate is a T17 row with a
provisional band and an explicit denominator, and each closure carries a
mandatory PHYSICS.md entry (§7) and its own validation test (§5).

> **Discipline note (CLAUDE.md prime directive; PHYSICS.md;
> `.claude/skills/lbmflow-physics-discipline`).** Every closure below is
> resolved from the governing energy equation or a literature-backed closure
> with a recorded derivation, validity domain, and a dedicated validation test.
> **No band-calibrated constant, no case-keyed branch, no transport-absorbing
> clamp** appears anywhere in this design. The three parameters that look like
> constants — `Pr_t = 0.85`, `cs_s² = 1/4` (D3Q7), `Λ = 1/4` (TRT magic) — are
> each a literature/lattice-derivation value with a fixed validity domain, not
> a fit (§0, §6). The Arrhenius `k(T)` and the property functions `μ(T)/ρ(T)/
> σ(T)` are physical constitutive laws whose coefficients (`E_a`, `β_T`, `σ_T`,
> `B`) are **user/material inputs from the scenario**, never chosen to pass a
> band; the validation anchors (§5) test that the *implemented functional form*
> matches the analytic law, not that a tuned number lands in a window. The only
> clamp discussed anywhere — a temperature floor — is **banned** in §6.4
> exactly as WSCAL bans the negative-`C` clamp.

---

## 0. Summary of decisions (read this first)

| # | Decision | Justification (short) |
|---|---|---|
| **T1** | **`T` reuses the WSCAL D3Q7 `h`-machinery unchanged** — same `Lattice` impl (`cs_s²=1/4`, weights `w_0=1/4, w_{1..6}=1/8`), same BGK/TRT collide, same `exchange_f_generic::<D3Q7,T>` halo, same post-`f` sub-step slot. The temperature distribution is a **second `h`-style set** `hT` (D3Q7 × 2), with a macroscopic-`T` compact field `temp`. | `active-scalar-feedback.md` §4.2: temperature is *mathematically isomorphic* to `C_k` on the ADE infrastructure; the only differences are (a) feedback targets `σ_T/β_T/μ(T)` and (b) an enthalpy source `Ṡ_T`. Reusing WSCAL's landed machinery is ~zero new lattice code and inherits its VR-STR-04 accuracy gates (Taylor–Aris, mass conservation, D3Q7 isotropy) for free. |
| **T2** | **Diffusivity mapping uses thermal diffusivity `α`** (not molecular `D`): `τ_T = α/cs_s² + ½ = 4α + ½` (D3Q7). Effective thermal diffusivity under LES adds the SGS heat flux `α_t = ν_t/Pr_t`: `α_eff = α + ν_t/Pr_t`, `τ_T = α_eff/cs_s² + ½`. | Energy ADE (§1); FR-LES-04 already names `Pr_t` as the SGS **heat** flux closure paired with `Sc_t` for scalars. `α = k_c/(ρ c_p)` (thermal conductivity / (ρ·specific heat)). `Pr_t = 0.85` is the standard turbulent-Prandtl value (Kays 1994; not a fit — §6.5). This is the *exact same code path* as WSCAL's `Sc_t`, with `Sc_t → Pr_t` and `D → α`. |
| **T3** | **Reaction enthalpy source** `Ṡ_T = −Σ_j ΔH_{r,j} r_j / (ρ c_p)` is a **split-step source** added into `hT` after the `hT` collide-stream (operator splitting), read from the WREACT reaction rates `r_j` in the SAME solver step. Exothermic `ΔH_r<0` ⇒ `Ṡ_T>0` (temperature rises). | Energy balance derivation §1.2. Operator-split source is the standard ADE-LBM reaction/source treatment (Krüger §8.3.5): the transport LBE is unchanged, `Ṡ_T·Δt` is added to `T` (equivalently distributed onto `hT` via the equilibrium weights). Uses WREACT's `r_j` — the same rates that feed the species `R_k(C)` — so heat release and species consumption are consistent by construction. |
| **T4** | **Arrhenius `k(T) = A_j exp(−E_{a,j}/(𝓡 T))`** consumed by WREACT. This is the reason thermal is in scope: WREACT's rate `r_j = k_j(T)·∏ C^{ν}` reads the `temp` field per cell. `T` is **Kelvin (absolute)** throughout the thermal axis (non-negotiable — Arrhenius is undefined for `T≤0`). | The Arrhenius law is the governing constitutive relation for reaction-rate temperature dependence (Arrhenius 1889; Levenspiel *Chemical Reaction Engineering* 1999). `A_j` (pre-exponential), `E_{a,j}` (activation energy), `𝓡` (gas constant) are per-reaction material inputs from the WREACT scenario block. Two-way loop: `T→k(T)→r_j→Ṡ_T→T`. |
| **T5** | **Thermal property feedback** (`μ(T)`, `ρ(T)/β_T`, `σ(T)/σ_T`) reuses the `active-scalar-feedback.md` §4.1 closures **verbatim**, composed **additively with the solutal feedback** (double-diffusive): thermal buoyancy `F_b^T = ρ_0 β_T(T−T_0)g` adds to `F_b^{scalar}`; thermal Marangoni `∇σ=(∂σ/∂T)∇T` superposes linearly on the surfactant `∇σ=(∂σ/∂C)∇C` (§4.1 linear-superposition rule) inside the **single Convention-D1 variable-σ force** (never double-adding a normal capillary term). `μ(T)` multiplies the FR-LES-02 viscosity composition as one more factor `f_T(T)=μ(T)/μ_0`. | `active-scalar-feedback.md` §3 (additive Boussinesq), §4.1 (thermal forms), §1.3 Convention D1 (double-counting avoidance). Reusing the same force-composition and viscosity-composition points means thermal feedback introduces **no new force path** — it feeds the existing `force_field` and `τ_eff` accumulation. |
| **T6** | **Thermal wall BCs**: Dirichlet (fixed `T_w`) via **anti-bounce-back** (lands first — needed for the conduction MMS and Rayleigh–Bénard hot/cold walls); Neumann (fixed flux `q_w`) via a **bounce-back + source correction**; Robin (convective `−k_c ∂T/∂n = h_c(T−T_∞)`) as a linear combination of the two (**deferred to a second wave** — needs the Neumann flux closure validated first). Zero-flux (adiabatic) is the WSCAL bounce-back reused. | Mirrors WSCAL §2 BC menu (bounce-back / anti-bounce-back / zero-gradient) with the thermal-specific flux/convective closures added. Dirichlet + adiabatic are the minimal set the conduction MMS (VR-STR-12) and natural-convection (VR-STR-14) gates need; Robin lands with its own convective-cooling validation. |
| **T7** | **Option-gated for B-6 bit-identity**: `hT`, `hTtmp`, `temp` are `Option<Vec<T>>` in `SoaFields`, `None` ⇒ **bit-identical to the isothermal engine** (the ablation gate B-6, §5 VR-STR-15). Thermal-off is a hard zero-cost path — no allocation, no sub-step call, `probe_state_hash` unchanged. | Same discipline as WSCAL P3 (`h`/`htmp`/`conc` all `Option`) and the `force_field`/`omega_field` precedent (`fields.rs:196-199`). The PM's mandatory ablation gate (thermal off → bit-identical) is a structural `Option`-None guarantee, not a runtime branch. |
| **T8** | **CPU-first (CpuScalar reference → CpuSimd fused), GPU deferred**, identical staging posture to WSCAL P8 (gated on B-1). | Do not block the thermal axis on GPU multi-distribution upload (B-1 PARTIALLY RESOLVED). `Backend::Fields` reserves the storage (`backend.rs:130-135`). |
| **T9** | **Coupling mode = predictor-corrector for the active thermal path** (REQ FR-COUP-01 active dataflow), degrading to within-step explicit for the passive-transport-only case. The `k(T)→r_j→Ṡ_T→T` loop and the `T→ρ/μ/σ` property loop are **both** in the active predictor-corrector sequence; time-lagged is a flagged relaxation (`active_scalar_lagged`) with a dt-halving convergence gate. | REQ FR-COUP-01: temperature "is another active scalar feeding property updates" (§5). The thermal reaction-heat loop is stiff for high `E_a` (exotherm runaway); the predictor-corrector + subcycling + the `Δt_r`/thermal-CFL constraint (§6) is the REQ-mandated treatment. dt-halving convergence is the MJ-007 negative test extended to `T`. |
| **T10** | **Single temperature field, phase-wise thermal deferred.** Phase 1 lands **one** `T` (single-phase or a single mixture temperature). Two-phase phase-wise energy (`α_q ρ_q c_{p,q} T_q`, latent heat of phase change) is API-reserved, gated on W-VOF, NOT built. | Minimal scope (CLAUDE.md). The reaction-engineering goal needs single-mixture `T` for Arrhenius; latent heat / phase-change enthalpy blurs the focus and waits on W-VOF exactly as `active-scalar-feedback.md` §4.2 noted. `Ṡ_T` in phase 1 carries reaction heat only (no latent-heat term). |

---

## 0.5. FR-THERM REQ amendment block (PROPOSED — for PM ratification)

> **⚠ PROPOSED REQ AMENDMENT — NOT YET AUTHORITATIVE.** The block below is
> written in the exact functional-requirement style of
> `REQ_STIRRED_REACTOR.md` §4 so that the PM can append it verbatim to REQ §4
> as a new subsection **§4.9 (W-THERM thermal / energy axis)**, and add the
> corresponding rows to the §1 runtime matrix, §11 DAG, and §8 VR list. Until
> the PM appends it, this is a proposal only. The PM directive of 2026-07-07
> (thermal upgraded from API-reserved to in-scope for the reaction-engineering
> goal) is the authorization to *draft* this; ratification = appending it.
>
> **Companion edits the PM must make on ratification** (called out so they are
> not forgotten):
> - **§1 runtime matrix**: change the *Scalar* row's `[T]` bracket and the
>   *Phase-2 (API-reserved)* line's "thermal axis" entry — thermal moves from
>   "Phase 2 (API-reserved)" to the initial-delivery / in-scope set.
> - **§10**: strike "Thermal axis recommended as API-reserved extension" and
>   replace with a pointer to FR-THERM (this block).
> - **§11 DAG**: add the `W-THERM` row (deps below) and a parallel-wave entry.
> - **§8**: add VR-STR-12..15 (and note VR-STR-08 doubles as the `σ(T)` gate).

### §4.9 W-THERM — thermal / energy axis (PROPOSED)

- **FR-THERM-01 (temperature transport)**: The solver shall evolve an absolute
  temperature field `T` [K] by the advection–diffusion–reaction-heat equation
  ```
  ∂T/∂t + u·∇T = ∇·[(α + ν_t/Pr_t)∇T] + Ṡ_T
  ```
  on the **shared D3Q7 ADE-LBM distribution and relaxation machinery** used for
  species scalars `C_k` (REQ §3; `WSCAL_PASSIVE_SPEC.md`) — temperature is a
  scalar with feedback targets and a heat source, NOT a distinct lattice type.
  Molecular thermal diffusivity `α = k_c/(ρ c_p)`; the SGS heat flux is closed
  with the turbulent Prandtl number `Pr_t` (default `0.85`), paired with the
  scalar `Sc_t` under FR-LES-04. The diffusivity→relaxation mapping is
  `τ_T = (α + ν_t/Pr_t)/cs_s² + ½`. `T` is measured in Kelvin throughout;
  configuration and I/O accept a unit tag but the core is Kelvin-only.

- **FR-THERM-02 (reaction enthalpy source)**: The heat source shall be the
  reaction-enthalpy release
  ```
  Ṡ_T = −(1/(ρ c_p)) Σ_j ΔH_{r,j} r_j
  ```
  where `r_j` are the WREACT reaction rates [mol·m⁻³·s⁻¹] and `ΔH_{r,j}` the
  molar reaction enthalpies [J·mol⁻¹] (`ΔH_r<0` exothermic ⇒ `Ṡ_T>0`).
  `Ṡ_T` shall be evaluated from the SAME `r_j` that drive the species source
  `R_k(C)` in the same solver step (heat/species consistency by construction),
  and added as an operator-split source (FR-COUP-01). Phase 1 carries reaction
  heat only; latent heat of phase change is API-reserved (waits on W-VOF).
  The adiabatic-batch limit `ΔT_ad = −ΔH_r ΔC/(ρ c_p)` shall be a validation
  anchor (VR-STR-13).

- **FR-THERM-03 (Arrhenius rate dependence)**: WREACT reaction rate constants
  shall carry the Arrhenius temperature dependence
  ```
  k_j(T) = A_j exp(−E_{a,j} / (𝓡 T))
  ```
  read per cell from the `T` field (`A_j` pre-exponential, `E_{a,j}` activation
  energy [J·mol⁻¹], `𝓡` gas constant). This closes the two-way thermal–kinetic
  loop `T→k(T)→r_j→Ṡ_T→T`. An optional modified-Arrhenius form
  `k_j(T)=A_j T^{n_j} exp(−E_{a,j}/(𝓡T))` shall be schema-reserved (a single
  power-law factor; NOT built in phase 1 unless a reaction needs it). Rate
  evaluation shall guard `T>0` at the WREACT contract boundary (FR-EXT-01
  NaN/divergence detection), NOT by clamping `T`.

- **FR-THERM-04 (thermal property feedback)**: When the active thermal path is
  enabled, property fields shall depend on `T`, composed **additively with the
  solutal feedback** (double-diffusive), reusing the `active-scalar-feedback.md`
  §4.1 closures:
  - **Viscosity** `μ(T)=μ_0 exp[B(1/T − 1/T_0)]` (Arrhenius, liquid-viscosity
    standard) or linear `μ(T)=μ_0[1−b(T−T_0)]`; enters the FR-LES-02/FR-STRESS-04
    viscosity composition as one multiplicative factor `f_T(T)=μ(T)/μ_0`.
  - **Density / thermal buoyancy** `F_b^T = ρ_0 β_T(T−T_0)g` (Boussinesq),
    added to the momentum force **additively with the solutal `F_b^{scalar}`**;
    exactly `0` at `T≡T_0` (so VR-STR-06 static stratification is preserved with
    thermal ON — a mandatory degeneracy gate). Non-Boussinesq is deferred
    (`|β_T ΔT|≲0.1` validity, per active-scalar-feedback §3.3).
  - **Surface tension / thermal Marangoni** `σ(T)=σ_0+σ_T(T−T_0)` (`σ_T=∂σ/∂T`,
    usually `<0`), evaluated inside the **single Convention-D1 variable-σ
    interfacial force** with `∇σ=(∂σ/∂T)∇T` superposed **linearly** on the
    surfactant `∇σ=(∂σ/∂C)∇C`. The `∇σ=0` degeneration MUST match the
    constant-σ reference (Convention D2; VR-STR-06+). Surface-tension feedback
    is meaningful only with W-VOF present.

- **FR-THERM-05 (thermal wall boundary conditions)**: The scalar wall menu
  (FR-BC-04) shall be extended for temperature with: **Dirichlet** (fixed wall
  temperature `T_w`, anti-bounce-back), **Neumann** (fixed wall heat flux `q_w`,
  bounce-back + source correction), **Robin/convective**
  (`−k_c ∂T/∂n = h_c(T−T_∞)`, linear combination), and **adiabatic** (zero-flux,
  bounce-back). Dirichlet and adiabatic are the initial-delivery set; Robin
  lands with its own convective-cooling validation. The half-way wall placement
  and 1-cell solid rim (CLAUDE.md invariant) apply to the thermal BCs
  identically.

- **FR-THERM-06 (coupling / time integration)**: The thermal axis shall follow
  the FR-COUP-01 **active predictor-corrector** dataflow (temperature is an
  active scalar feeding property updates and reaction rates). The stiff
  reaction-heat loop shall respect a thermal-diffusion stability bound
  (`τ_T>½`, grid-Péclet), a reaction-heat timestep `Δt_r`, and — when thermal
  buoyancy is active — the buoyancy `Δt_b` of active-scalar-feedback §5.1 with
  `β_C→β_T`. Time-lagged explicit thermal feedback is permitted only as the
  flagged relaxation `active_scalar_lagged=true` with a dt-halving convergence
  benchmark (MJ-007 extended to `T`).

- **FR-THERM-07 (thermal-off bit-identity — ablation gate)**: With the thermal
  axis disabled the engine shall be **bit-identical** to the isothermal engine
  (`probe_state_hash` unchanged) — a hard zero-cost `Option`-None path, no
  allocation and no sub-step. This is validation gate **B-6** (VR-STR-15).

**Runtime-matrix / DAG deltas (for §1, §11)**: *Thermal* axis — fidelity
default `resolved-energy` (D3Q7 ADE + reaction heat + property feedback);
relaxation `passive-transport-only` (transport with no feedback / no reaction
heat). DAG: `W-THERM` hard-deps `W-SCAL` (h-machinery) and `W-REACT` (rates +
enthalpies); active feedback path additionally needs `W-VOF` for `σ(T)`/two-phase.
Parallel with W-REACT once W-SCAL lands.

*(End of proposed FR-THERM block.)*

---

## 1. Governing equation + energy balance + enthalpy source

### 1.1 The temperature advection–diffusion–reaction-heat equation

For an absolute temperature `T(x,t)` [K] transported by the resolved
incompressible velocity `u`, with thermal diffusivity `α` and a volumetric heat
source `Ṡ_T`:

```
∂T/∂t + u·∇T = ∇·[(α + ν_t/Pr_t)∇T] + Ṡ_T                               (1)
```

This is the constant-property (single-mixture) energy equation written in
temperature form. It is **term-for-term the WSCAL passive ADE (2)** with three
substitutions and one addition:

| WSCAL passive scalar | Thermal axis | Note |
|---|---|---|
| `C` (concentration) | `T` (temperature, K) | macroscopic zeroth moment `T=Σ_i hT_i` |
| `D` (molecular diffusivity) | `α` (thermal diffusivity `k_c/(ρc_p)`) | §2 mapping `τ_T=4α+½` |
| `ν_t/Sc_t` (SGS scalar flux) | `ν_t/Pr_t` (SGS heat flux) | `Pr_t=0.85`; FR-LES-04 |
| — | `+ Ṡ_T` (reaction heat) | §1.2; split-step source (T3) |

**Derivation of (1) from the energy balance.** Start from the constant-property
thermal-energy equation (incompressible, negligible viscous dissipation and
pressure work, single mixture):

```
ρ c_p (∂T/∂t + u·∇T) = ∇·(k_c ∇T) + q̇                                  (1a)
```

`q̇` [W·m⁻³] = volumetric heat generation. Dividing by `ρ c_p` (constant) and
defining thermal diffusivity `α ≡ k_c/(ρ c_p)` [m²·s⁻¹] gives
`∂T/∂t + u·∇T = α∇²T + q̇/(ρ c_p)`. Adding the LES SGS heat flux as a gradient
diffusion `−⟨u'T'⟩ = (ν_t/Pr_t)∇T` (FR-LES-04) and identifying
`Ṡ_T ≡ q̇/(ρ c_p)` [K·s⁻¹] yields (1). Viscous dissipation `Φ/(ρc_p)` and
compressibility (pressure work `βT Dp/Dt`) are dropped: valid in the low-Mach,
low-Eckert regime `Ec = U²/(c_p ΔT) ≪ 1` (the LBM operating regime,
`Ma_lattice≤0.1`) — recorded as a validity-domain limitation in §6.3 and
PHYSICS.md, NOT a silent omission.

### 1.2 The reaction-enthalpy source `Ṡ_T` (decision T3, FR-THERM-02)

For a set of reactions `j` with molar rates `r_j` [mol·m⁻³·s⁻¹] and molar
reaction enthalpies `ΔH_{r,j}` [J·mol⁻¹]:

```
q̇ = −Σ_j ΔH_{r,j} r_j        [W·m⁻³]        (exothermic ΔH_r<0 ⇒ q̇>0)
Ṡ_T = q̇/(ρ c_p) = −(1/(ρ c_p)) Σ_j ΔH_{r,j} r_j     [K·s⁻¹]           (2)
```

**Sign / units check.** Exothermic reaction: `ΔH_r<0`, `r_j>0`, so
`−ΔH_r r_j>0` ⇒ `q̇>0` ⇒ heat released ⇒ `T` rises. Units:
`[J·mol⁻¹]·[mol·m⁻³·s⁻¹] = J·m⁻³·s⁻¹ = W·m⁻³` ✓; divided by
`[kg·m⁻³]·[J·kg⁻¹·K⁻¹] = J·m⁻³·K⁻¹` gives `K·s⁻¹` ✓ (a temperature rate,
matching the LHS `∂T/∂t`).

**Consistency with the species source (load-bearing).** The `r_j` in (2) are
**the same rates** WREACT uses for the species production `R_k = Σ_j ν_{kj} r_j`
(`ν_{kj}` = stoichiometric coefficient). Heat and species are therefore
consistent by construction: a mole reacted removes species and releases exactly
its enthalpy. This is why `Ṡ_T` reads WREACT's `r_j` directly rather than
recomputing — a recomputation could drift from the species path and is banned.

**Adiabatic batch limit (the VR-STR-13 anchor).** For a well-mixed adiabatic
batch (`u=0`, no diffusion loss, single reaction A→products consuming `ΔC` of A),
integrating (2) with `r = −dC_A/dt` gives the closed-form adiabatic temperature
rise
```
ΔT_ad = −ΔH_r ΔC_A / (ρ c_p)                                            (3)
```
independent of the rate law — a pure enthalpy-balance identity. VR-STR-13 gates
the implemented `Ṡ_T` against (3).

### 1.3 ADE-LBM for `T` — the `hT` distribution (decision T1)

Carry a **second D3Q7 set** `hT_i(x,t)`, `i=0..6`, identical in structure to the
WSCAL `h` set (same `Lattice` impl, same weights `w_0=1/4, w_{1..6}=1/8`, same
`cs_s²=1/4`), relaxed toward the linear-in-velocity equilibrium

```
hT_i^eq = w_i^s T [ 1 + (c_i·u)/cs_s² ]                                  (4)
T = Σ_i hT_i     (zeroth moment = macroscopic temperature)
```

BGK (4→5 in WSCAL) and TRT with magic parameter `Λ=1/4` (default for the
accurate-conduction gate) are inherited **unchanged** from WSCAL §1.2 — the only
difference is the relaxation rate is set from `α` not `D` (§2) and the
reaction-heat source is applied post-stream (§3). Chapman–Enskog on this scheme
recovers (1) minus `Ṡ_T` (the transport part); the source is added by operator
splitting (T3).

**Source application (operator splitting, decision T3).** After the `hT`
collide-stream produces the transported `T*`, apply the source over one step:
```
T^{n+1} = T* + Ṡ_T Δt                                                    (5)
```
distributed back onto `hT_i` via `hT_i += w_i^s Ṡ_T Δt` (adds `Ṡ_T Δt` to the
zeroth moment while leaving higher moments of the increment isotropic —
the standard ADE-LBM zeroth-moment source, Krüger §8.3.5). This is a
**forward-Euler split source**; its split error is `O(Δt)` and is bounded by the
`Δt_r` reaction-heat timestep (§6.2) — for stiff exotherms the source substep is
subcycled or the predictor-corrector (T9) is used, per FR-COUP-01.

---

## 2. Data structures (decision T1, T2, T7)

### 2.1 Reuse of the landed WSCAL machinery (no new lattice)

`T` mounts on the **same D3Q7 `Lattice` impl** WSCAL O1 adds (`lattice.rs`);
the **same** `exchange_f_generic::<D3Q7,T>` halo (`halo.rs:308`); the **same**
post-`f` solver sub-step slot (`solver.rs`, after `update_moments`); the **same**
collide/stream/BC kernel rows parameterized by lattice + relaxation rate. No new
lattice type, no new halo path, no new step-order change beyond what WSCAL
already introduced. This is the T1 isomorphism made concrete.

### 2.2 The W-THERM additions to `SoaFields<T>` (decision T7)

Add to `SoaFields<T>` (`fields.rs:168`), all `Option<…>` so `None` is
bit-identical to the isothermal / thermal-free path (B-6 invariance):

```rust
/// Temperature ADE distribution set (D3Q7), q-major padded planes.
/// `None` ⇒ isothermal (no allocation, bit-identical to the isothermal engine).
/// Shares the D3Q7 Lattice impl, halo, and sub-step machinery of the WSCAL `h` set.
pub h_t: Option<Vec<T>>,
/// Ping-pong partner of `h_t`. Temperature streaming writes here, then swapped.
pub h_t_tmp: Option<Vec<T>>,
/// Macroscopic temperature T = Σ_i h_t_i, compact core, Kelvin (absolute).
/// T > 0 physically (positivity is a diagnostic + a hard WREACT-contract guard,
/// NOT a transport-absorbing clamp — see §6.4).
pub temp: Option<Vec<T>>,
```

Optional per-cell thermal-relaxation field (the LES `α_t = ν_t/Pr_t` hook,
mirroring WSCAL's phase-2 `omega_s_field` and the landed `omega_field`
precedent `fields.rs:199`):

```rust
/// Per-cell thermal relaxation ω_T = 1/τ_T when LES SGS heat flux is active
/// (τ_T from α_eff = α + ν_t/Pr_t). `None` ⇒ uniform molecular α. Phase-2 add,
/// exactly analogous to WSCAL's Sc_t hook; phase 1 uses a uniform τ_T.
pub omega_t_field: Option<Vec<T>>,
```

Placement/naming mirror the WSCAL `h`/`htmp`/`conc` triplet and the W-VOF
`g`/`gtmp`/`phi` triplet exactly, so `f` / `g` / `h` (scalar) / `h_t`
(temperature) form a uniform additive-`Option` pattern; the four orders touch
**disjoint fields** of the same struct (§8 coexistence).

### 2.3 Feedback hooks (Option-gated; consumed by existing composition points)

Thermal feedback writes into the **existing** force / viscosity accumulation —
no new force path (decision T5):

- **`F_b^T` (thermal buoyancy)** — accumulated into `force_field` (`fields.rs:196`)
  alongside the solutal `F_b^{scalar}` and `F_s`/`ρg`, in the frozen
  FORCE_COMPOSITION_SPEC summation order (§4.2). A `thermal_buoyancy` closure
  object (`β_T`, `T_0`, `g`) gates it; `None` ⇒ no contribution.
- **`f_T(T)=μ(T)/μ_0`** — one multiplicative factor in the FR-LES-02 /
  FR-STRESS-04 viscosity composition (§4.1), applied at the `τ_eff` assembly
  point, composed with the solutal `f_C(C)` multiplicatively.
- **`σ(T)`** — feeds the Convention-D1 variable-σ interfacial force via
  `∇σ += (∂σ/∂T)∇T`; a `thermal_sigma` closure (`σ_T`) gates it. Meaningful only
  with W-VOF present; `None` ⇒ no thermal-Marangoni contribution.

### 2.4 Memory cost per cell (D3Q7 `h_t`, matches NFR-01)

Per the REQ NFR-01 budget row "Scalar h (per component) — D3Q7 × 2 × f32 = 56":
temperature is exactly one more scalar component.

| Component | Layout | B/cell (f32) |
|---|---|---|
| Temperature `h_t` (D3Q7 × 2) | 7 × 2 × f32 | **56** |
| `temp` (T, compact core) | 1 × f32 | 4 |

**≈ 60 B/cell for the thermal axis** — one scalar-component slot in the NFR-01
budget. At 1e8 cells ≈ 6 GB (f32); ≤256³ dev/validation ≈ 1 GB. Under the
mixed-safe precision profile (NFR-02), `T` inside the interface band and all
thermal reductions (mean `T`, `ΔT_ad` counters) promote to f64.

---

## 3. Solver-step slot + coupling to WREACT and property updates (decision T3, T9)

### 3.1 Where the thermal sub-step slots in (FR-COUP-01 active dataflow)

The thermal axis extends the WSCAL post-`f` sub-step. Per solver step, active
predictor-corrector (REQ FR-COUP-01 active branch, T9):

```
1. (if active) scalar/reaction PREDICTOR + property update ρ(C,T), μ(C,T), σ(C,T):
   a. evaluate Arrhenius k_j(T) from the CURRENT temp field (FR-THERM-03).
   b. evaluate reaction rates r_j = k_j(T)·∏C^ν (WREACT).
   c. property update: f_T(T), F_b^T, σ(T) written into the feedback hooks (§2.3).
2. FORCE COMPOSITION: F_s + ρg + F_b^{scalar} + F_b^T + …  (frozen order, §4.2).
3. HYDRODYNAMIC f STEP (unchanged run_span): collide → halo → stream → open BCs
   → update_moments.  Produces ρ, u (F/2-corrected).
4. SCALAR ADE SUB-STEP (WSCAL): h collide/stream/BC → C = Σ h_i.
5. TEMPERATURE ADE SUB-STEP (new, W-THERM), reusing the WSCAL h-machinery:
   a. hT collide (BGK/TRT, rate from τ_T): read T=Σ hT_i and resolved u; relax to (4).
   b. exchange hT halo: exchange_f_generic::<D3Q7,T>.
   c. hT stream (pull) → h_t_tmp; swap.
   d. thermal BCs (§ T6 / FR-THERM-05): Dirichlet / Neumann / adiabatic / Robin.
   e. T* = Σ_i hT_i.
   f. reaction-heat source (5): Ṡ_T = −Σ_j ΔH_{r,j} r_j /(ρ c_p) from step 1b's r_j;
      T^{n+1} = T* + Ṡ_T Δt   (subcycled/predictor-corrector if stiff, T9).
6. (if active) CORRECTOR: re-evaluate k_j(T^{n+1}), r_j, R_k, properties;
   → optional flow-scalar-thermal iteration for stiff coupling (FR-COUP-01).
```

**Ordering rationale (decided, physical).** `k(T)` is evaluated from the `T`
available at step 1 (predictor) and re-evaluated at step 6 (corrector) — this is
the FR-COUP-01 active-scalar predictor-corrector, extended so that temperature
(an active scalar feeding both properties AND reaction rates) participates. The
reaction-heat source (5f) uses the `r_j` computed in the same step (1b), so heat
and species consumption are synchronized. The temperature sub-step runs **after**
the species scalar sub-step (both read the resolved `u`); the property/force
feedback it produces is consumed by the **next** step's force composition (or,
under strong coupling, the within-step iteration). For the passive-transport-only
relaxation mode (no reaction, no feedback), steps 1/2/6 collapse and `T` is a
one-way within-step coupling exactly like WSCAL.

### 3.2 Coupling to WREACT: the `k(T)` contract

WREACT owns the rate law; W-THERM owns the `T` field. The contract (FR-EXT-01
signature discipline):

- **W-THERM provides**: the per-cell `temp` field (Kelvin), guaranteed `>0` at
  the read boundary (the WREACT rate evaluator receives `T` and asserts
  `T>0`/finite — a contract-boundary NaN/divergence guard, NOT a `T` clamp).
- **WREACT provides**: `r_j(T, C)` [mol·m⁻³·s⁻¹] via Arrhenius `k_j(T)=A_j
  exp(−E_{a,j}/(𝓡T))`, and the enthalpies `ΔH_{r,j}` from the scenario.
- **W-THERM consumes**: `r_j` and `ΔH_{r,j}` to form `Ṡ_T` (step 5f), reading the
  SAME `r_j` array WREACT used for `R_k` (step 4/1b) — no recomputation.

If WREACT is absent (thermal transport only, no reactions), `Ṡ_T≡0` and the
thermal axis is a pure ADE (VR-STR-12/14 need no reactions). If W-THERM is absent
(isothermal), WREACT uses a fixed reference-temperature rate `k_j(T_0)` — the
Arrhenius factor collapses to a constant; this is the degenerate isothermal-
kinetics path.

### 3.3 Property updates (decision T5) — see §4

The `T→ρ/μ/σ` property feedback (step 1c) is detailed in §4; it reuses the
`active-scalar-feedback.md` closures and composition points verbatim, additively
with the solutal feedback.

---

## 4. Property-feedback closures μ(T)/ρ(T)/σ(T) (decision T5, FR-THERM-04)

All three reuse `active-scalar-feedback.md` §4.1 **verbatim** and compose with
the solutal feedback per the rules there. This section states the composition;
it does **not** re-derive the closures (they are fixed in that document, which is
the single source of truth for the functional forms and their literature).

### 4.1 μ(T) — thermal viscosity (additive-in-log with solutal)

```
μ(T) = μ_0 exp[B(1/T − 1/T_0)]        (Arrhenius, liquid-viscosity standard)
       or  μ_0 [1 − b(T−T_0)]          (linear, weak dependence)
f_T(T) = μ(T)/μ_0
```

Composed into the FR-LES-02 / FR-STRESS-04 viscosity assembly (active-scalar-
feedback §2.2) as **one more multiplicative molecular-viscosity factor**:
```
ν_mol(γ̇, C, T) = [μ(γ̇) · f_C(C) · f_T(T)] / ρ(C,T)
ν_eff = ν_mol + ν_t,   τ_eff = ½ + ν_eff/(cs²Δt)
```
Multiplicative because `μ(γ̇)`, `f_C(C)`, `f_T(T)` are all factors of molecular
viscosity (eddy `ν_t` is additive). The τ upper/lower clip (REQ FR-LES-03
`[τ_min,τ_max]`) applies to the final `τ_eff`; a clip activation is a
"property-out-of-range" diagnostic (active-scalar-feedback §5.2), never a hidden
transport cap.

### 4.2 ρ(T)/β_T — thermal buoyancy (additive Boussinesq → double-diffusive)

```
F_b^T = ρ_0 β_T (T − T_0) g,   β_T = −(1/ρ_0)(∂ρ/∂T)|_p  (thermal expansion)
```

Added into `force_field` **additively with the solutal** `F_b^{scalar} = ρ_0
β_C(C−C_0)g` — the double-diffusive form (active-scalar-feedback §3, §4.1). The
Boussinesq consistency conditions C-B1/C-B2 (active-scalar-feedback §3.2) apply
verbatim with `β_C→β_T`: `F_b^T≡0` at `T≡T_0` (exact-zero degeneration —
**mandatory** VR-STR-06 preservation with thermal ON), kept as an independent
perturbation force separate from the `ρ(φ)g` well-balanced cancellation.
Boussinesq validity `|β_T ΔT|≲0.1`; beyond that, warn (feasibility check) and
either narrow ΔT or defer to a non-Boussinesq path (API-reserved). Sign: `β_T>0`
(usual), hot fluid rises.

### 4.3 σ(T)/σ_T — thermal Marangoni (linear superposition, single Convention-D1 force)

```
σ(T) = σ_0 + σ_T (T − T_0),   σ_T = ∂σ/∂T  (usually < 0 for pure substances)
```

Evaluated inside the **single** Convention-D1 variable-σ interfacial force
(active-scalar-feedback §1.3), with the temperature contribution to the
tension gradient **superposed linearly** on the surfactant contribution:
```
∇σ = (∂σ/∂C)∇C + (∂σ/∂T)∇T                                              (6)
```
This is the linear-superposition rule of active-scalar-feedback §4.1/§7 (when
multiple σ models act on the same interface, the first-version rule is
additive = linear superposition). **Convention D1 double-counting rule is
load-bearing**: when *any* σ-feedback (solutal or thermal) is active, the whole
interfacial force uses the Convention-D1 chemical-potential-combined form — the
normal capillary term is computed **once** there; a separate `F_s=μ_φ∇φ` is
NOT added on top (that would double-count the normal contribution — the exact
error Convention D1 exists to prevent). The `∇σ=0` degeneration (both `∇C=0`
and `∇T=0`) must reproduce the constant-σ reference (Convention D2; VR-STR-06+
degeneracy gate). `σ(T)` feedback requires W-VOF (an interface must exist) — in
its absence it is rejected/inert, exactly as active-scalar-feedback §7 states
for the surfactant σ path.

**Open derivation inherited (active-scalar-feedback §8.1):** the Convention-D1
coefficient consistency (`W²↔κ` for the Liu et al. 2306.11320 combined form vs
REQ's `σ=√(2κβ)/6`, `W=4√(κ/(2β))`) is **derivation-required before implementing
the σ path** — this thermal spec does not resolve it; it inherits the same
open point and the same VR-STR-06+ degeneracy verification. The σ(T) order (if
dispatched) STOPS and reports if the derivation is unresolved (stop-rule flag,
§7 orders).

---

## 5. Validation plan → T17 + NEW VR-STR rows (decision T6)

Tests are **authored adversarially by codex/Opus from THIS spec**, in a test
worktree that never shares with the implementation worktree (CLAUDE.md; REQ §8).
Each row = metric / reference / **two-layer gate** (a metric band AND a behavior/
pattern gate) / grid / steps / backend / **explicit denominator** for every
relative error. Bands are **provisional MVP gates** (T17 band governance:
tightening always allowed; loosening requires a recorded PHYSICS.md rationale).

**New VR-STR rows proposed** (numbering continues the series; VR-STR-08..11 are
reserved by active-scalar-feedback §6 for the solutal/surfactant/thermal-capillary
items — **VR-STR-08 (thermal-capillary migration) is claimed here as the `σ(T)`
gate**, see below; the thermal-transport/reaction items are the new
VR-STR-12..15):

| ID | Test | Two-layer gate: metric band (denominator) + behavior gate | Grid / steps / backend | T17 row |
|---|---|---|---|---|
| **VR-STR-12** | **Thermal conduction/convection MMS** (manufactured solution). Impose an analytic `T(x,t)` (e.g. decaying sinusoid `T=T_0+A e^{−α k² t}cos(kx)` for pure conduction; add a uniform `u` for advection–diffusion), inject the corresponding source, measure grid-convergence order. | **Metric**: (a) L2 error vs analytic **< 1%** at the reference time, denominator = `‖T_analytic−T_0‖_2` (the analytic deviation-from-reference L2 norm, so the relative error is not inflated by the `T_0` offset); (b) **grid-convergence order ≥ 1.9** (Δx→Δx/2 at fixed physical time, error ratio ≥ ~3.5). TRT `Λ=1/4` and BGK both run; TRT L2 ≤ BGK L2. **Behavior**: the conducted field stays symmetric/isotropic (no D3Q7 lattice-aligned diamond artifact — the isotropy claim); the advected sinusoid shows no leading/trailing asymmetry (Galilean). | 1D/2D-in-3D `256×…` & `512×…`, periodic transverse, `α` s.t. `τ_T∈[0.6,0.9]`, CpuScalar | **VR-STR-12 (new)** → VR-STR-04 family |
| **VR-STR-13** | **Adiabatic reaction temperature rise vs analytic** (batch exotherm). Well-mixed adiabatic box (`u=0`, adiabatic walls), single reaction A→P with `ΔH_r<0`, run to completion. Measure the final `ΔT`. | **Metric**: `\|ΔT_measured − ΔT_ad\| / ΔT_ad` **< 2%**, `ΔT_ad = −ΔH_r ΔC_A/(ρ c_p)` (eq. 3), denominator = the analytic `ΔT_ad`. **Behavior**: the `T(t)` curve is monotone-increasing and its *shape* tracks the conversion `ΔC_A(t)` (heat release synchronized with species consumption — the T3 consistency claim); no overshoot beyond `ΔT_ad` (an overshoot signals a heat/species desync or a split-error blowup). | closed box `32³` (well-mixed), run to `>99%` conversion, CpuScalar | **VR-STR-13 (new)** → VR-STR-04 (reaction) |
| **VR-STR-14** | **Thermal-buoyancy natural convection vs Rayleigh correlation.** Rayleigh–Bénard (or side-heated cavity): hot/cold Dirichlet walls (FR-THERM-05), adiabatic sidewalls, Boussinesq `F_b^T`. Measure the Nusselt number vs Rayleigh number. | **Metric**: `Nu` vs the reference `Nu–Ra` correlation/benchmark within **±10%** (denominator = reference `Nu`); at least 3 `Ra` values spanning ~10⁴–10⁶, and the **onset** `Ra_c≈1708` (Rayleigh–Bénard) reproduced within band. Reference: de Vahl Davis (1983) side-heated cavity benchmark, or the classical `Nu∝Ra^{1/4}` scaling. **Behavior**: convection cells appear only above `Ra_c` (subcritical → conduction-only, `Nu≈1`); cell count/orientation physically plausible; `F_b^T≡0` and quiescence at `ΔT=0` (degeneracy). | 2D/3D cavity `128²`/`64³`, steady, CpuScalar | **VR-STR-14 (new)** → VR-STR-04/06 |
| **VR-STR-15** | **Arrhenius rate-vs-T anchor.** Isothermal batches at a sweep of fixed `T`; measure the initial reaction rate `r_0(T)`. | **Metric**: `ln r_0` vs `1/T` is linear with **slope = −E_a/𝓡** within **±5%** (denominator = the analytic `E_a/𝓡`); i.e. the implemented `k(T)` reproduces the Arrhenius slope, and rate **increases with T** (correct sign). Anchored over ≥4 temperatures spanning a decade in rate. **Behavior**: the Arrhenius plot is a straight line (no curvature — confirms the exponential form, not a linearized surrogate); `r` strictly monotone-increasing in `T`. | well-mixed `16³`, fixed `T` each run, CpuScalar | **VR-STR-15 (new)** → VR-STR-04 (reaction) |
| **VR-STR-08** (claimed here) | **Thermal-capillary droplet migration** (Young–Goldstein–Block). The `σ(T)` gate, already sketched in active-scalar-feedback §6. Requires W-VOF. | **Metric**: terminal velocity `V` vs `V_YGB = 2U/[(2+k̃)(2+3μ̃)]`, `U=−σ_T G_T R/μ_B`, **L2rel(V) < 5%** at Ma→0 (denominator = analytic `V_YGB`), + 2nd-order gradient convergence. **Behavior**: droplet migrates toward the **hot** side for `σ_T<0` (correct thermocapillary direction); `∇σ=0` (uniform T) degeneration reproduces the constant-σ static droplet (Convention D2). | droplet R=20–40 l.u., `μ̃=k̃=1`, `G_T` const, Ma=Re=O(0.1), CpuScalar; **needs W-VOF** | **VR-STR-08** (active-scalar-feedback §6) |

**Mandatory negative / consistency / ablation tests (the two-layer gate's
second layer + the PM-mandated ablation):**

- **B-6 ablation (VR-STR-15 companion, FR-THERM-07)** — **thermal OFF ⇒
  bit-identical to the isothermal engine.** With `h_t=None`/`temp=None`, the
  `probe_state_hash` on any isothermal scenario (cavity, cylinder presets) is
  **bit-identical** to the pre-W-THERM engine. This is the PM-mandated ablation
  gate: thermal is a zero-cost `Option`-None path. (Distinct row; call it
  **B-6 / VR-STR-15-ablation**.)
- **Reaction-heat/species consistency (VR-STR-13 negative arm):** a mutant that
  computes `Ṡ_T` from an independently recomputed rate (not the WREACT `r_j`
  array) must FAIL a stricter conversion-vs-ΔT synchronization check — proves the
  T3 "same `r_j`" contract is load-bearing.
- **Buoyancy exact-zero degeneracy (VR-STR-06 extension):** with thermal buoyancy
  ON and `T≡T_0`, `max|u|<ε` static stratification is preserved (`F_b^T` exact
  zero) — the FR-THERM-04 mandatory degeneracy; a mutant with a non-zero-at-`T_0`
  buoyancy must FAIL.
- **σ(T) ∇σ=0 degeneracy (VR-STR-06+ extension):** with `σ(T)` ON and uniform `T`
  (and uniform `C`), the variable-σ Convention-D1 force must reproduce the
  constant-σ reference (no double-counted normal capillary term) — Convention D2.
- **Arrhenius sign / `T>0` guard:** a rate law that decreases with `T`, or a `T≤0`
  fed to Arrhenius without the contract guard tripping, must FAIL (proves the
  sign and the boundary guard).
- **dt-halving thermal-feedback convergence (MJ-007 extended):** on an
  active-thermal bench (natural convection or exotherm), the feedback coupling
  error converges under `Δt→Δt/2` — the FR-COUP-01 active-scalar convergence gate
  extended to the thermal loop.
- **`ν_t/Pr_t` SGS-flux guard:** in phase 1 (uniform `τ_T`), with LES active the
  thermal field transports with molecular `α` only; a test asserts `α_eff==α`
  (no silent SGS heat flux) until the `omega_t_field` phase-2 hook lands.

### 5.6 Behavior-validity review (mandatory, REQ / CLAUDE.md)

After each thermal validation run, before reporting: review the *observed*
pattern, not just the gated metric (the second layer above formalizes this per
test). Specifically: (a) conduction stays isotropic (no D3Q7 diamond); (b) the
exotherm `T(t)` tracks conversion and never overshoots `ΔT_ad`; (c) natural
convection cells appear only above `Ra_c` and quiescence holds at `ΔT=0`;
(d) the Arrhenius plot is straight and monotone; (e) thermocapillary migration
is toward the hot side for `σ_T<0`. Record the review in PHYSICS.md or the
track's findings file. A metric passing its band does **not** validate a pattern
no band covers (the natural-convection cell structure and the exotherm curve
shape are pattern gates no single scalar band covers — origin-of-directive case).

---

## 6. Stability & validity domain

### 6.1 Thermal-diffusion relaxation window (thermal CFL)

`τ_T = α/cs_s² + ½ = 4α + ½` (D3Q7, `cs_s²=1/4`), under LES `τ_T=4α_eff+½` with
`α_eff=α+ν_t/Pr_t`. Operating band `τ_T∈(0.5,~1.0]`, identical to WSCAL §6.1.
Near `0.5` (`α→0`, advection-dominated heat transport) BGK loses positivity; TRT
`Λ=1/4` widens the usable window. The diffusive stability limit
`α Δt/Δx² ≤ ½` (von Neumann) is the continuous-form statement of `τ_T>½`; the
LBM realizes it through the relaxation-time floor. Report `τ_T` and `Λ` per run;
freeze in PHYSICS.md after characterization.

### 6.2 Grid-Péclet + reaction-heat timestep

- **Thermal grid-Péclet** `Pe_Δ^T = |u|Δx/α` governs positivity of the advected
  temperature exactly as WSCAL §6.2 (`Pe_Δ≲2` BGK, wider TRT). A scenario
  violating it is refined or switched to TRT, **never clamped**.
- **Reaction-heat timestep** `Δt_r`: the split source (5) is forward-Euler; for a
  stiff exotherm the temperature rate `Ṡ_T` can be large, so the step must
  resolve the reaction-heat timescale `τ_heat ~ (ρc_p ΔT_ad)/(|ΔH_r| r_max)` — if
  `Ṡ_T Δt` would exceed a fraction of the local `T`, the source is subcycled or
  the predictor-corrector iteration (T9) is invoked (FR-COUP-01). This is a
  resolution/integration requirement, not a tunable band.
- **Thermal buoyancy** `Δt_b ≤ C_b √(Δx/|β_T ΔT g|)` — the active-scalar-feedback
  §5.1 buoyancy constraint with `β_C→β_T` (`C_b` inherited, still needs
  calibration there; back-calculated via VR-STR-14).

### 6.3 Dropped-term validity (Eckert / low-Mach)

Viscous dissipation and pressure work are dropped from (1) — valid for
`Ec = U²/(c_p ΔT) ≪ 1` and low-Mach (`Ma_lattice≤0.1`), the LBM regime. This is
a **documented model limitation** (PHYSICS.md §7), not a hidden approximation:
a scenario with large `Ec` (e.g. high-speed viscous heating) is out of the
thermal axis's validity domain and must be flagged, not silently run.

### 6.4 Temperature positivity is a diagnostic + a contract guard, NOT a clamp

`T>0` (Kelvin) is physical and is a **hard guard at the WREACT contract
boundary** (Arrhenius is undefined for `T≤0`; the rate evaluator asserts finite
`T>0` — FR-EXT-01 NaN/divergence detection). But the transport step **must not
clamp `T`**: a `T`-floor that silently absorbs a negative excursion is a banned
transport-absorbing clamp (CLAUDE.md prime directive;
`.claude/skills/lbmflow-physics-discipline` ban list). A `T<0` beyond round-off
is a *symptom* of a grid-Péclet / `τ_T` / split-source-stiffness violation
(§6.1–6.2) and must be surfaced as a diagnostic and fixed by resolution /
collision / subcycling, **never masked**. The `temp` doc-comment (§2.2) states
this. (Direct analog of WSCAL §6.5 negative-`C`.)

### 6.5 `Pr_t = 0.85` is a literature value, not a fit

`Pr_t=0.85` is the standard turbulent Prandtl number for the SGS heat flux
(Kays 1994; the value paired with `Sc_t=0.7` in FR-LES-04). It is a closure
parameter with a validity domain (near-unity `Pr_t` for most turbulent flows),
overridable per scenario, reported in metadata — NOT calibrated to pass a band.
Likewise `cs_s²=1/4` and `w_0=1/4, w_{1..6}=1/8` are D3Q7 lattice-derivation
constants (WSCAL §1.4) and `Λ=1/4` is the TRT magic value (Ginzburg 2005).

---

## 7. CODEX ORDER BREAKDOWN

Four orders, file-conflict-aware. One order = one bundle = one dedicated
worktree (CLAUDE.md). Implementation and adversarial-test orders never share a
worktree. **Every physics-affecting order embeds the physics-discipline clauses
(lbmflow-codex-dispatch Step 1.5): citation + derivation + validity + own test +
PHYSICS.md entry; NO calibrated constants / case branches / transport-absorbing
clamps; STOP and report if a gate can't be met without a hack.**

**Prerequisite (hard dependency for O1/O2):** `W-SCAL` (D3Q7 h-machinery) MUST
be landed — O1 reuses its `Lattice` impl, halo, and sub-step slot. `W-REACT`
(rate law + `r_j`/`ΔH_r`) MUST be landed or co-developed — O2 consumes its rate
contract. If either is absent, O1/O2 STOP and report (stop-rule flag).

| Order | Scope | Primary files (conflict boundary) | DoD / gate |
|---|---|---|---|
| **O1 — Temperature ADE transport + reaction-enthalpy source (CpuScalar)** | `h_t`/`h_t_tmp`/`temp` `Option` slots in `SoaFields` (§2.2); `hT` BGK+TRT collide reusing the WSCAL D3Q7 `Lattice` (rate from `τ_T=4α+½`, eq. 4); `hT` halo via `exchange_f_generic::<D3Q7>`; `hT` stream+swap; `T=Σ hT_i`; the temperature sub-step wiring at solver level (§3.1 step 5, after the WSCAL scalar sub-step); **operator-split reaction-heat source** `Ṡ_T=−Σ ΔH_r r_j/(ρc_p)` reading a passed-in `r_j` array (eq. 5); thermal wall BCs Dirichlet(anti-BB)+adiabatic(BB)+Neumann(BB+source); uniform `τ_T`; PHYSICS.md entry (§7.entry). **NOT in O1: Robin BC (second wave), Arrhenius k(T), property feedback, LES α_t hook.** | `fields.rs` (add h_t/h_t_tmp/temp/omega_t_field), `solver.rs` (temperature sub-step orchestration + source apply), `kernels.rs` (hT collide/stream/BC + source rows) — **same three files WSCAL O1 touched, but disjoint new fns/fields (append `collide_ht_row`, `apply_thermal_bc`, etc.); textual both-add merge only**. | VR-STR-12 (conduction/convection MMS order≥1.9), VR-STR-13 (adiabatic ΔT_ad with a stubbed constant `r_j`), B-6 ablation (h_t=None bit-identity) green on CpuScalar. |
| **O2 — Arrhenius k(T) coupling + thermal property feedback** | Arrhenius `k_j(T)=A_j exp(−E_a/(𝓡T))` in the WREACT rate path (FR-THERM-03), reading `temp`; the `T>0` contract guard; the two-way `T→k(T)→r_j→Ṡ_T` wiring (feed O1's source with WREACT's real `r_j`); property feedback `f_T(T)` (μ composition factor, §4.1), `F_b^T` (thermal-buoyancy contributor into `force_field`, additive with solutal, §4.2), `σ(T)` contributor into the Convention-D1 variable-σ force (§4.3) — **`σ(T)` STOPS and reports if the Convention-D1 `W²↔κ` derivation (§4.3 / active-scalar-feedback §8.1) is unresolved and W-VOF absent**; predictor-corrector coupling hook (T9). | WREACT rate module (Arrhenius), `solver.rs` (property update step 1c + force-composition contributor — **coordinate with the active-scalar-feedback solutal contributor at the same accumulation point; additive**), viscosity-composition point, Convention-D1 force point (shared with W-VOF/active-scalar — land after those). Depends: O1 + W-REACT + (for σ) W-VOF. | VR-STR-15 (Arrhenius slope ±5%), VR-STR-13 with real k(T) exotherm, VR-STR-14 (natural convection Nu–Ra) green; VR-STR-08 (σ(T) migration) green **only when W-VOF present**; buoyancy/σ degeneracy gates green. |
| **O3 — Scenario schema thermal block + T / reaction-rate outputs** | Scenario JSON `temperature` block (enabled, T_0, α or (k_c, c_p), Pr_t, BCs per face with `type∈{dirichlet,neumann,robin,adiabatic}` + values in Kelvin, unit tag) and per-reaction `arrhenius{A,E_a}` + `dH_r` in the WREACT block; config validation (Kelvin>0, thermal grid-Péclet §6.2 warning, Boussinesq `|β_T ΔT|≲0.1` warning, `σ(T)` requires W-VOF reject, Robin-not-yet reject); CLI/output of `temp` field (VTI, manifest.json), reaction-rate `r_j` and `Ṡ_T` field outputs; `ΔT_ad` / mean-T reduction outputs. | `crates/lbm-scenario/src/lib.rs` (schema + validation), `crates/lbm-cli` (temp/r_j/Ṡ_T output) — **disjoint files from O1/O2**. Depends: O1 (field to output), O2 (Arrhenius params to parse). | schema round-trips; validation rejects/warns per the rules; `temp`/`r_j`/`Ṡ_T` appear in a run's manifest; VR-STR-14 scenario (Rayleigh–Bénard) authored via the schema. |
| **O4 — Adversarial validation authorship (codex, separate worktree)** | All of §5 (VR-STR-12/13/14/15 + VR-STR-08 σ-gate) + the negative/consistency/ablation tests (B-6 bit-identity, heat/species-consistency mutant, buoyancy exact-zero, σ ∇σ=0 degeneracy, Arrhenius sign + T>0 guard, dt-halving thermal convergence, ν_t/Pr_t SGS-leak guard). Authored from **THIS spec**, not the impl. Freeze VR-STR-12..15 bands in VALIDATION.md T17 on landing. | `crates/lbm-core/tests/wtherm_*.rs` + `crates/lbm-scenario/tests/*` (new files only — no impl-file conflict). Runs concurrently from the start (test worktree); tests compile red against stubs, go green as O1/O2/O3 land. | tests green against O1/O2/O3; bands frozen in T17; behavior-validity review (§5.6) recorded for each run. |

**Critical-path ordering:** `W-SCAL` (prereq) → **O1** → **O2**; **O3** depends on
O1 (field) and O2 (params); **O4** runs concurrently from the start (test
worktree). `W-REACT` is a hard dep of O2 (rate contract) — if W-REACT is not yet
landed, O1 can land standalone (transport + source with a stubbed `r_j`), and O2
lands when W-REACT is ready. The `σ(T)` sub-item of O2 is gated on W-VOF and on
the Convention-D1 derivation (STOP-RULE flag). CpuSimd fused + GPU are follow-on
orders out of this plan's scope (T8, gated on `backend_simd_equiv.rs`/T13 and
B-1 respectively).

### 7.entry — The PHYSICS.md validity-domain statements (mandatory)

**O1 must add** a PHYSICS.md §1 stack entry + §2 decision entry:

> **Temperature transport — ADE-LBM (Krüger 2017 §8.3), D3Q7, BGK/TRT, on the
> shared WSCAL scalar machinery.** Absolute temperature `T` [K] is advected by
> the resolved (F/2-corrected) `u` and diffused with thermal diffusivity
> `α=k_c/(ρc_p)`; distribution `hT_i` (D3Q7, `cs_s²=1/4`, `w_0=1/4,w_{1..6}=1/8`),
> linear-in-`u` equilibrium, `T=Σ hT_i`; `τ_T=α/cs_s²+½=4α+½` (`4α_eff+½` under
> LES with `α_eff=α+ν_t/Pr_t`, `Pr_t=0.85` Kays 1994). Reaction-enthalpy source
> `Ṡ_T=−Σ_j ΔH_{r,j}r_j/(ρc_p)`, operator-split forward-Euler, read from the SAME
> WREACT `r_j` as the species source (heat/species consistency). BCs:
> Dirichlet=anti-BB, adiabatic=BB, Neumann=BB+source (Robin second wave).
> **Validity**: `τ_T∈(0.5,~1.0]`, thermal grid-Péclet `≲2` (BGK, wider TRT),
> `Ma_lattice≤0.1`, low-Eckert `Ec≪1` (viscous dissipation + pressure work
> dropped — documented limitation). `T>0` is a diagnostic + WREACT-contract guard,
> NEVER a clamp (analog of WSCAL negative-C). SGS heat flux OFF in phase 1
> (molecular `α` only even under LES — first phase-2 add). Record the measured
> conduction MMS order, ΔT_ad match, and the frozen `τ_T`/`Λ`.

**O2 must add** a §2 decision entry for the Arrhenius + feedback closures:

> **Arrhenius rate + thermal property feedback.** `k_j(T)=A_j exp(−E_a/(𝓡T))`
> (Arrhenius 1889; Levenspiel 1999); closes the two-way `T→k(T)→r_j→Ṡ_T→T` loop.
> `A_j,E_{a,j}` are material inputs (never band-fit). Property feedback reuses
> `active-scalar-feedback.md` §4.1 verbatim, additive with solutal:
> `μ(T)=μ_0 exp[B(1/T−1/T_0)]` (one multiplicative viscosity factor `f_T`);
> `F_b^T=ρ_0 β_T(T−T_0)g` (Boussinesq, additive with `F_b^{scalar}`, exact-zero at
> `T≡T_0`, validity `|β_T ΔT|≲0.1`); `σ(T)=σ_0+σ_T(T−T_0)` with
> `∇σ=(∂σ/∂C)∇C+(∂σ/∂T)∇T` linearly superposed inside the SINGLE Convention-D1
> variable-σ force (normal capillary counted once — Convention D1/D2 double-count
> rule; `∇σ=0`→constant-σ reference). Convention-D1 `W²↔κ` coefficient
> consistency is derivation-required before the σ path lands (inherited open
> point, active-scalar-feedback §8.1). Record the Arrhenius slope match, Nu–Ra
> match, and any characterized coefficients.

---

## 8. Coexistence with WSCAL / WREACT / WVOF (structural summary)

W-THERM (`h_t`, D3Q7, post-`f` sub-step, active feedback + reaction heat) mounts
alongside the other distribution axes without structural conflict:

- **Storage:** `f` (D3Q19/27 hydro) / `g` (D3Q19 phase-field, W-VOF) / `h`
  (D3Q7 scalar, W-SCAL) / `h_t` (D3Q7 temperature, W-THERM) are **four `Option`
  distribution sets** in one `SoaFields`, each with a `None`-default bit-identity
  guarantee (B-6). Adding `h_t` is purely additive — the second-to-land order
  rebases the struct field group + `new()` init + module fns mechanically (same
  both-add resolution WSCAL §7.1 documents for W-VOF).
- **Step slots:** phase-field pre-pass (before `f`) ≠ `f` step ≠ scalar sub-step
  (after `f`) ≠ **temperature sub-step (after the scalar sub-step)**. Four
  disjoint slots (§3.1). W-THERM's temperature sub-step is a sibling of WSCAL's
  scalar sub-step, running immediately after it (both read the same resolved `u`).
- **Halo:** `h_t` reuses `exchange_f_generic::<D3Q7,T>` — the **same** path
  WSCAL's `h` uses (a D3Q7 population set), instantiated on a different buffer.
  No shared buffer, no new halo path.
- **WSCAL reuse (not conflict):** W-THERM depends on WSCAL's D3Q7 `Lattice` impl,
  halo, and sub-step scaffolding being landed — it *reuses* them. If WSCAL lands
  first (expected), W-THERM adds a second `h`-style set with zero new lattice
  code. The two are the SAME machinery on two buffers with two relaxation rates
  (`D→τ_s`, `α→τ_T`) and two feedback-target sets.
- **WREACT coupling:** the `k(T)→r_j→Ṡ_T` loop (§3.2) is where W-THERM and
  W-REACT interlock — WREACT owns the rate law and `r_j`; W-THERM owns `T` and
  consumes `r_j`/`ΔH_r` for `Ṡ_T`. Single shared `r_j` array (no recompute).
- **Force / feedback:** thermal buoyancy `F_b^T` accumulates into the SAME
  `force_field` as W-VOF's `F_s`/`ρg` and active-scalar's `F_b^{scalar}`, in the
  frozen FORCE_COMPOSITION_SPEC order — **additive** (double-diffusive). The
  thermal-Marangoni `∇σ` term superposes linearly inside the SINGLE Convention-D1
  variable-σ force W-VOF+active-scalar establish (never a second normal capillary
  term). This is the only place the four axes' physics interact, and it is the
  active-feedback phase (O2), gated on W-VOF for the σ path.
- **Two-phase phase-wise energy** (`α_q ρ_q c_{p,q} T_q`, latent heat) is
  API-reserved, gated on W-VOF, NOT built (T10) — the single-`T` `temp` slot and
  the `φ`-availability from the `g` path reserve it, exactly as WSCAL P10 reserves
  the phase-wise scalar form.

---

## 9. Load-bearing references (grounding index)

| Claim | Source |
|---|---|
| `T` reuses D3Q7 h-machinery; thermal ≅ solutal ADE isomorphism | `active-scalar-feedback.md` §4.2; `WSCAL_PASSIVE_SPEC.md` §1.4, P1/P2 |
| Thermal property closures μ(T)/ρ(T)/σ(T), additive with solutal | `active-scalar-feedback.md` §4.1, §3 (Boussinesq), §1.3 Convention D1/D2 |
| Convention-D1 `W²↔κ` derivation-required (inherited open point) | `active-scalar-feedback.md` §8.1 |
| VR-STR-08 thermal-capillary migration sketch (the σ(T) gate) | `active-scalar-feedback.md` §6 |
| `Pr_t` SGS heat flux closure paired with `Sc_t` | REQ FR-LES-04; `active-scalar-feedback.md` §4.1 |
| Temperature is an active scalar feeding property updates (dataflow) | REQ FR-COUP-01 (§5) active branch |
| D3Q7 `Lattice`, halo, sub-step slot, Option-field/B-6 discipline | `WSCAL_PASSIVE_SPEC.md` §1.4/§3/§4, P3/P8; `fields.rs:168-210`, `:196-199`; `backend.rs:130-135`, `:258-300`; `halo.rs:308` |
| Force composition point / frozen order for `F_b` contributors | `FORCE_COMPOSITION_SPEC.md` T5, §2 R1–R4; REQ §3 `F_b^{scalar}` |
| Scalar wall BC menu (extended for thermal) | REQ FR-BC-04; `WSCAL_PASSIVE_SPEC.md` §2 |
| NFR-01 D3Q7 56 B/cell budget (temperature = one scalar component) | REQ §7 NFR-01 budget table |
| Thermal-off bit-identity ablation (B-6) | `WSCAL_PASSIVE_SPEC.md` V5 (h=None); REQ NFR-04 / FR-COUP-04 probe_state_hash |
| MJ-007 dt-halving convergence (extended to thermal) | REQ §8 mandatory negative tests; FR-COUP-01 |
| Positivity-as-diagnostic-not-clamp discipline (analog of neg-C) | `WSCAL_PASSIVE_SPEC.md` §6.5; CLAUDE.md prime directive; `.claude/skills/lbmflow-physics-discipline` ban list |

**Literature (decided references):**
- Arrhenius, S. 1889 (*Z. Phys. Chem.* 4:226) — the rate–temperature law
  `k=A exp(−E_a/𝓡T)`. Levenspiel, *Chemical Reaction Engineering* (3rd ed., 1999)
  — reaction-engineering standard for `k(T)`, `ΔH_r`, adiabatic temperature rise.
- Krüger, Kusumaatmaja, Kuzmin, Shardt, Silva & Viggen 2017, *The Lattice
  Boltzmann Method* §8.3 — ADE-LBM equilibrium, `τ` mapping, source term (§8.3.5),
  BCs (adopted governing discretization, shared with WSCAL).
- Kays, W.M. 1994 (*J. Heat Transfer* 116:284) — turbulent Prandtl number
  `Pr_t≈0.85` (the SGS-heat-flux closure value).
- Ginzburg 2005 (*Adv. Water Resour.* 28:1171) — TRT for ADE, magic `Λ=1/4`
  (wall location + isotropy); inherited from WSCAL.
- de Vahl Davis, G. 1983 (*Int. J. Numer. Methods Fluids* 3:249) — differentially
  heated square-cavity natural-convection benchmark (VR-STR-14 `Nu`–`Ra`
  reference). Rayleigh–Bénard onset `Ra_c≈1708` (Chandrasekhar 1961).
- Young, N.O., Goldstein, J.S. & Block, M.J. 1959 (*J. Fluid Mech.* 6:350) —
  thermocapillary droplet migration terminal velocity (VR-STR-08 analytic).
- Liu, Wu, Ba, Xi et al. 2023 (*Phys. Rev. E* 108:055306) — well-balanced
  variable-σ (incl. linear σ–T thermocapillary) phase-field LBM (the
  Convention-D1 combined force; inherited from active-scalar-feedback §1.3/§4.1).
