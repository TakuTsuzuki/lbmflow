# T3-Tier Implementation Specification — Multicomponent Interfacial Mass Transfer

**Document ID**: SPEC-XFER-T3 (rev.1, 2026-07-07).
**Scope**: the shared **closure library** for gas–liquid interfacial mass
transfer required by `docs/REQ_STIRRED_REACTOR.md` **FR-VOF-05** ("Interfacial
mass transfer separated: resolved interface (normal flux + Henry partition +
phase-wise diffusion) vs point-bubble (`k_L a(C*−C)`). Henry and Sherwood
applicability explicit."). This spec owns the **Henry partition**, the
**partial-pressure output**, and the **Sherwood-correlation film-coefficient
closure** consumed by BOTH branches; it does not own W-VOF's phase-field
transport (`WVOF_IMPL_SPEC.md`) or W-BUB's population-balance / bubble-motion
model (`WBUB_PBM_IMPL_SPEC.md`, **not yet written** — see §8 dependency note).
**Target core**: `crates/lbm-core` (species-registry-consuming closure
functions; no new distribution set of its own — this tier calls into the
W-SCAL `h`/`conc` machinery per component and, on the resolved branch, into
the W-VOF `φ` field).
**Acceptance**: VALIDATION.md **T17** row **VR-STR-04** (scalar/reaction,
"`k_L a` (formula = interface integral or correlation, explicit)") and
**VR-STR-02a** (single-bubble, extended here to multi-gas absorption/
stripping). New rows are proposed in §5 as **T17 additions** pending PM
freeze — this spec does not itself edit VALIDATION.md.

This spec is **executable**: every closure is decided and justified, every
Rust API is written against the registry design in
`WSCAL_MULTICOMPONENT_SPEC.md` — **which does not exist yet at the time of
writing** (see §8) — with the minimal registry surface this spec needs stated
explicitly so the two specs can be reconciled without a redesign. No band is
calibrated; every constant carries a citation.

> **Discipline note (CLAUDE.md prime directive; PHYSICS.md;
> `.claude/skills/lbmflow-physics-discipline`).** Every closure below is
> either resolved from a diffusion/partition governing equation or a
> literature-backed correlation with a recorded derivation, stated validity
> domain (`Re_b`, `Sc`, `Eo`), and its own validation test (§5). **No
> species-name branch appears anywhere** — `pO2`/`pCO2`/`pH2`/`pN2` are
> printed *instances* of one generic per-species output, never special-cased
> in code. Non-volatile species get zero flux from a registry **flag**, never
> from a name comparison. No transport-absorbing clamp. The mandatory
> PHYSICS.md entry text is in §7.2; the behavior-validity review checklist is
> in §5.6.

---

## 0. Summary of decisions (read this first)

| # | Decision | Justification (short) |
|---|---|---|
| T1 | **Henry's law partition is per-species, driven by a registry `volatile: bool` flag — never a name comparison.** Non-volatile species (glucose, salts, biomass) get an interfacial flux of **exactly 0** by construction (the flux closure is not evaluated for them, not clamped to zero after evaluation). | REQ FR-VOF-05 + CLAUDE.md ban list ("case-identity branch... model must not know which test case it is running"). A `match species.name { "O2" => ... }` branch is exactly the banned pattern; a registry boolean carried per species is the resolved-physics equivalent (the species either partitions across a gas-liquid interface or it does not — that is a property of the molecule, encoded once at registration, not re-derived per call site). |
| T2 | **Two independent transfer branches sharing one closure library**, gated by which multiphase mode is active (REQ §1 Interface axis): **resolved-interface** (W-VOF `φ`-field present) uses normal diffusive flux + Henry jump condition; **point-bubble** (W-BUB Lagrangian/Eulerian bubbles, W-VOF absent or bubbles below `d_b/Δx` threshold) uses `k_L a(C*−C)`. A run configured with neither active multiphase mode present gets **no interfacial transfer term** (single-phase passive scalar, W-SCAL only) — not an error, just the FR-VOF-05 term evaluating over an empty interface/bubble set. | REQ §1 interface axis: `resolved-phasefield` (fidelity) vs `point-bubble` (relaxation) are mutually exclusive per FR-VOF-04 switching criteria (`d_b/Δx, d_b/W, Eo, Re_b, α_g, We_b`); FR-VOF-05 restates the same split specifically for mass transfer. |
| T3 | **Partial pressure `p_Xi = C_i / H_i` is the ONE general output function**, called identically for every volatile species in the registry; O2/CO2/H2/N2 in the envelope (§6) are registry *entries*, not hardcoded cases. Units and the dimensionless-vs-dimensional Henry convention are fixed once (§1.1) so `p_Xi` is unambiguous. | Task requirement: "do NOT hardcode any [gas]." Matches T1's registry-driven design — the same per-species loop produces flux AND partial pressure. |
| T4 | **Henry coefficient convention: `H_i` in `[Pa·m³/mol]` (or the SI-consistent lattice equivalent), defined by `p_i = H_i C_i^{liquid}` (concentration-based Henry's law, the "Henry volatility" `H^cp` inverse convention — see Sander 2015 compilation).** Equilibrium liquid concentration at a given gas partial pressure is `C*_i = p_i / H_i`. | Sander (2015, *Atmos. Chem. Phys.* 15:4399) catalogs Henry constants in multiple conventions; picking ONE (concentration-based, dimensional) and stating it once removes the single most common bug class in gas-liquid transfer code (unit-convention mismatches). This spec's `H_i` is **always** this convention; the registry stores it pre-converted (§2.1) so no call site re-derives units. |
| T5 | **Resolved-interface flux is the phase-wise diffusion equation's own diffusive flux, jump-matched by Henry's law at the interface — not a separate source term bolted onto the bulk ADE.** The interfacial condition is a **boundary condition on the phase-wise scalar system** (REQ §3 "Two-phase phase-wise conservative" form), not an extra `S^if` volumetric term computed independently of the diffusion solve. | REQ §3 explicitly frames `S^if` as appearing in the phase-wise conservative ADE (`∂(α_q C_{k,q})/∂t + ... = ... + S_{k,q}^{if}`) with `S_{k,liq}^{if} = −S_{k,gas}^{if}`. Modeling it as an interface-band source (§3.2) that enforces the Henry jump self-consistently is the only way to guarantee `Σ_q∫α_q C_{k,q}dV` conservation (REQ §3, VR-STR-05) without double-counting the diffusive flux the ADE solver already computes across the diffuse interface. |
| T6 | **Point-bubble film coefficient `k_L` from the Higbie penetration theory as the primary closure, Frössling as a cross-check / high-`Re_b` alternative — both stated, Higbie decided as default** because it needs no empirical exponents beyond `Sc` and is the standard choice for freely-circulating bubbles in the Grace-diagram regime this engine's VR-STR-02a already targets. | Higbie (1935) penetration model: `k_L = 2√(D_i/(π t_c))`, contact time `t_c = d_b/U_t` (Kawase & Moo-Young 1990 review; Kantarci, Borak & Ulgen 2005 bioreactor mass-transfer review). Frössling (1938) / Clift-Grace-Weber (1978) boundary-layer correlation `Sh = 2 + 0.6 Re_b^{1/2} Sc^{1/3}` is the literature alternative for rigid/contaminated spheres; stated as the validity-domain cross-check (§6) because real bioreactor bubbles are often contaminated (surfactant-immobilized surface), which breaks Higbie's free-circulation assumption — this is flagged as risk R1 in §6. |
| T7 | **Interfacial area `a` from Sauter mean diameter `d_32` and gas holdup `ε_g`: `a = 6 ε_g / d_32`.** `d_32` and `ε_g` are NOT computed by this spec — they are consumed from whichever upstream source is active (W-VOF interface segmentation per REQ FR-VOF-04, or W-BUB's PBM `d_32`, or a user-specified single bubble diameter for VR-STR-02a). | Standard bubble-column/stirred-tank interfacial-area relation (Kantarci et al. 2005; Garcia-Ochoa & Gomez 2009 review, their Eq. for `a`). This spec's job is to consume `d_32, ε_g` and emit `k_L a`, not to produce `d_32` — that is FR-VOF-04's job (already assigned to W-VOF/W-BUB), avoiding a duplicate computation this task explicitly warns against ("this spec OWNS the Henry+Sherwood closure definitions; WBUB consumes them" — not the reverse). |
| T8 | **Multicomponent point-bubble gas-side balance: each bubble carries a per-species gas-phase inventory `n_{i,gas}` (mol); each species partitions independently via its own `H_i, D_i, k_{L,i}`,** with the *only* cross-species coupling being the shared bubble volume/pressure (ideal-gas or given-EOS closure relating `Σ_i n_{i,gas}` to bubble volume at the local hydrostatic + Laplace pressure). No species-pair interaction term (Henry partition is dilute-solution, additive across species by construction — Sander 2015 §2). | REQ FR-VOF-05 + task requirement "Multicomponent gas-side balance: each bubble's gas composition evolves as species partition in/out." Dilute multicomponent Henry's law has no cross terms (each solute's activity coefficient is referenced to infinite dilution in the *solvent*, not to other solutes) — this is a resolved consequence of the dilute-solution assumption already implicit in Henry's law itself, not an added approximation. |
| T9 | **Ablation is a first-class configuration state, not a test-only hack**: `transfer: Option<TransferConfig>` at the top of both branches; `None` ⇒ species are inert across the interface (bit-identical to W-SCAL/W-VOF/W-BUB running with no transfer term at all). | Mirrors the `Option`-gated `h`/`g` precedent (WSCAL_PASSIVE_SPEC P3, WVOF_IMPL_SPEC) — B-6 invariance discipline. Required negative test in §5 (ablation: transfer off → inert). |
| T10 | **This spec produces no new distribution set and no new `Lattice` impl.** It is a closure/BC layer that (a) on the resolved branch, modifies the phase-wise scalar BC/source at the interface band (consuming W-SCAL's per-phase `h_liq`/`h_gas` — itself a W-SCAL phase-2 add, flagged as a dependency in §8); (b) on the point-bubble branch, is a per-bubble ODE right-hand-side term consumed by W-BUB's bubble integrator. | Minimal scope (CLAUDE.md). The task brief is explicit that this is "the SHARED closure library both... paths use" — a library, not a new transport mechanism. |

---

## 1. Governing transfer equations

### 1.1 Henry's law partition (the interface equilibrium)

For a volatile species `i` at a gas–liquid interface, local thermodynamic
equilibrium (fast compared to bulk transport — the standard two-film/
penetration assumption, Lewis & Whitman 1924) sets the **interfacial**
liquid-side concentration equal to the value in equilibrium with the local
gas-phase partial pressure `p_i`:

```
C*_i = p_i / H_i                                                        (1)
```

with the **decided convention (T4)**: `H_i` is the concentration-based Henry
volatility, `[H_i] = Pa·m³·mol⁻¹` (equivalently `atm·L·mol⁻¹` etc. — the
registry stores SI-consistent `H_i` after unit conversion at load time, §2.1),
defined by

```
p_i = H_i · C_i^{liquid, equilibrium}                                    (1')
```

`H_i` is temperature-dependent (van't Hoff form, `H_i(T) = H_i(T_ref)
exp[-d(ln H_i)/d(1/T) · (1/T - 1/T_ref)]`) — **phase 1 of this spec takes
`H_i` as a registry constant at the run's fixed temperature** (isothermal
reactor, consistent with REQ's current scope which has no landed thermal
axis); the `T`-dependence hook is noted in §6 as a validity-domain bound, not
implemented.

**Non-volatile species** (registry `volatile = false`): equation (1) is never
evaluated; the species has no gas-phase partial pressure and no interfacial
flux (T1). This is the resolved consequence of "non-volatile" — the species
literally has no equilibrium vapor pressure in the model, not a flux forced
to zero after being computed.

### 1.2 Partial-pressure output (general, for every volatile species)

The inverse of (1) is the **general diagnostic output** required by the task
("partial pressure `p_Xi = C_i/H_i` as the general output for ALL volatile
species"):

```
p_{X_i} = C_i · H_i                                                      (2)
```

**Sign/direction note**: equation (1) computes the *equilibrium liquid
concentration* from a *given* gas partial pressure (used inside the flux
closures, §1.3–1.4, where the gas-side state is the input). Equation (2)
computes the *partial pressure implied by* a *given* liquid concentration
(used as the reporting/diagnostic output — "what pO2 does this measured
dissolved-O2 concentration correspond to"). Both are the same equilibrium
relation `p_i = H_i C_i` read in the two directions the two use-sites need;
implemented as ONE function `partial_pressure(species, c_liquid) -> p_i =
c_liquid * H_i` (§2.2) and ONE function `equilibrium_conc(species, p_gas) ->
c_star = p_gas / H_i` — inverses of each other by construction, sharing the
same `H_i` lookup, so they cannot drift apart.

### 1.3 Resolved-interface branch: normal flux + Henry jump (W-VOF-gated)

At a diffuse gas–liquid interface (W-VOF `φ`, `φ=1` liquid / `φ=0` gas per
CLAUDE.md / REQ §3 convention), each phase diffuses species `i` with its own
phase diffusivity `D_{i,liq}`, `D_{i,gas}` (registry-supplied per species per
phase — a species generally diffuses faster in gas than in liquid, e.g. `D_gas
~ 10⁴×D_liq` for typical dissolved gases). The **phase-wise conservative
scalar form** (REQ §3, restated per-species):

```
∂(α_gas C_{i,gas})/∂t + ∇·(α_gas u C_{i,gas})
    = ∇·[α_gas D_{i,gas} ∇C_{i,gas}] + S_{i,gas}^{if}
∂(α_liq C_{i,liq})/∂t + ∇·(α_liq u C_{i,liq})
    = ∇·[α_liq D_{i,liq} ∇C_{i,liq}] + S_{i,liq}^{if}
α_liq = φ,  α_gas = 1 − φ
```

with the **interfacial coupling condition** (T5): continuity of normal
diffusive flux across the interface, with the *concentration* jump set by
Henry's law rather than continuity of `C_i` itself (the physical content of
"two immiscible phases in local equilibrium" — the flux, not the
concentration, is what a diffuse interface actually must match):

```
Flux continuity:      D_{i,gas} ∇C_{i,gas}·n̂  =  D_{i,liq} ∇C_{i,liq}·n̂     (3)
Henry jump:            C_{i,liq}|_interface = p_{i,gas}|_interface / H_i     (4)
                       (p_{i,gas} = C_{i,gas} R T  via ideal-gas law
                        relating the local gas-phase concentration to
                        partial pressure at temperature T)
Mass-conserving source: S_{i,liq}^{if} = −S_{i,gas}^{if}                    (5)
```

**Implementation as an interface-band source (T5, the practical closure)**:
directly enforcing (3)-(4) as a sharp jump condition inside a diffuse-interface
LBM is the classical VOF/phase-field interfacial-BC problem; the decided
closure (consistent with the diffuse-interface Allen-Cahn formalism REQ §3
already adopts for `ρ(φ), μ(φ)`) is a **smoothed interfacial source term**
confined to the interface band (the same `interface_band = max(3W, 6Δx)` band
NFR-02 already reserves for f64 promotion), analogous to a continuous surface
reaction:

```
S_{i,liq}^{if}(x) = k_int,i · |∇φ| · ( C*_i(x) − C_{i,liq}(x) )             (6)
S_{i,gas}^{if}(x) = −S_{i,liq}^{if}(x) · (α_liq/α_gas)  [mass-conserving
                     partition of the same molar flux across the two phase
                     volume fractions in the interface band cell]
```

`|∇φ|` localizes the source to the interface (zero away from it, by
construction — φ is uniform in bulk phases); `C*_i` is the local Henry
equilibrium value (1) evaluated from the *local* gas-side `C_{i,gas}` via the
ideal-gas relation in (4). `k_int,i` is an **interfacial exchange rate, not a
free/calibrated constant** — it is fixed by requiring (6) to recover the sharp
flux-continuity/Henry-jump limit (3)-(4) as `Δx/W → 0` (a resolved-limit
requirement, matched via the same asymptotic analysis used for the Allen-Cahn
mobility `M` in WVOF_IMPL_SPEC — the exact `k_int,i(D_{i,gas}, D_{i,liq}, W)`
closed form is an **open derivation item**, flagged as risk R2 in §6, blocking
implementation until derived and validated against test V3 in §5). This is
the honest state: the *equilibrium* physics (1)/(4) and the *conservation*
requirement (5) are fully resolved; the *numerical* realization of the jump
inside a diffuse interface is a stated open derivation, not a filled-in
constant — per the physics-discipline Rule 1 provenance table, this is
explicitly flagged rather than shipped with a guessed `k_int,i`.

### 1.4 Point-bubble branch: `k_L a(C*−C)` film model (W-BUB-gated)

For an unresolved (point) bubble, the classical two-film/penetration model
gives the molar transfer rate of species `i` into the liquid per unit
dispersion volume:

```
Ṡ_i = k_{L,i} · a · ( C*_i − C_{i,liq} )                                  (7)
```

`C*_i` is the Henry equilibrium (1) evaluated at the **bubble's own local
gas-phase partial pressure** `p_{i,bubble} = x_i · P_bubble` (`x_i` = mole
fraction of species `i` in the bubble's gas inventory, `P_bubble` = local
total pressure = hydrostatic + Laplace `2σ/r_b`). This is the resolved-physics
consequence of Dalton's law applied to the bubble's own multicomponent gas
mixture — no cross-species term (T8).

**Film coefficient `k_{L,i}` (Higbie penetration, decided default, T6):**

```
k_{L,i} = 2 √( D_{i,liq} / (π t_c) ),      t_c = d_b / U_t                (8)
```

`d_b` = bubble diameter (from `d_32`, T7), `U_t` = bubble terminal rise
velocity (already computed by the existing VR-STR-02a Grace correlation path
— REQ §3 point-bubble force balance / this engine's Grace-diagram terminal
velocity). `t_c` is the surface-renewal contact time for a bubble rising a
distance ~ its own diameter (Higbie 1935; Kawase & Moo-Young 1990).

**Cross-check / contaminated-surface alternative (Frössling, stated per T6,
not the default):**

```
Sh_i = k_{L,i} d_b / D_{i,liq} = 2 + 0.6 Re_b^{1/2} Sc_i^{1/3}             (9)
Re_b = ρ_l U_t d_b / μ_l,           Sc_i = μ_l / (ρ_l D_{i,liq})
```

(Frössling 1938; Clift, Grace & Weber 1978 Ch. 3 — the rigid-sphere
boundary-layer correlation.) Equation (8) (Higbie) applies to **clean,
freely-circulating** bubbles (internal circulation continuously renews the
surface); equation (9) (Frössling) applies to **rigid or contaminated**
bubbles (surfactants immobilize the surface, suppressing internal
circulation). The two differ by up to ~2× in the intermediate `Re_b` range —
this is validity-domain risk **R1** (§6), not resolved by picking one formula;
this spec's API (§2.3) exposes BOTH so the caller/scenario states which
surface-mobility regime applies, defaulting to Higbie per T6.

**Interfacial area (T7):**

```
a = 6 ε_g / d_32                                                         (10)
```

(surface-area-to-volume ratio of a Sauter-mean-diameter population — standard
bubble-column relation, Kantarci et al. 2005 Eq. 2; consumed, not derived,
per T7.)

### 1.5 Multicomponent gas-side balance (point-bubble, T8)

Each bubble's gas-phase mole inventory `n_{i,gas}` (mol per bubble) evolves by
the same flux (7), integrated over the bubble's own surface area `A_b = π
d_b²`, with sign flipped (species leaving the liquid enters the gas):

```
dn_{i,gas}/dt = − Ṡ_i · V_b   [mol/s per bubble; Ṡ_i in mol/(m³·s) is a
                               volumetric rate over the *dispersion*, so the
                               per-bubble form scales by the bubble's own
                               volume V_b = π d_b³/6 and by 1/ε_g of the local
                               cell to convert dispersion-volumetric to
                               per-bubble — see §2.3 for the exact per-bubble
                               vs per-dispersion-volume unit reconciliation]
x_i = n_{i,gas} / Σ_j n_{j,gas}          (mole fraction, for Dalton's law)
P_bubble = P_hydrostatic(z) + 2σ/r_b     (local total pressure)
p_{i,gas} = x_i · P_bubble               (Dalton's law, per-species)
```

`Σ_i n_{i,gas}` relates to bubble volume via the ideal gas law (dilute gas
at reactor conditions — `P_bubble V_b = (Σ_i n_{i,gas}) R T`), closing the
loop: as species partition out of the bubble, `Σn_{i,gas}` drops, `V_b` at
fixed `P_bubble` drops (bubble shrinks — e.g. a stripping N2 bubble absorbing
CO2 while losing O2 changes size), which is passed back to W-BUB's PBM/motion
model as an external volume-forcing term (an explicit **hand-off contract to
W-BUB**, not implemented here — §3.2).

---

## 2. Closure data structures (Rust API)

### 2.1 Minimal species-registry surface this spec needs

**Dependency note (§8): `WSCAL_MULTICOMPONENT_SPEC.md` does not exist at
time of writing.** This section states the *minimal* registry fields this
spec's closures require, so that whichever spec lands the registry first
implements a superset compatible with the other. If the multicomponent
registry spec is written before this one lands, its `Species` type MUST
contain (at least) these fields under these names or an explicitly documented
mapping:

```rust
/// One species entry in the multicomponent registry (minimal surface
/// this spec's closures require — see WSCAL_MULTICOMPONENT_SPEC.md for
/// the full registry, which is a superset: molar mass, reaction
/// stoichiometry, etc. are NOT needed here).
pub struct SpeciesTransferProps {
    /// Registry index / stable id — used as the array index into all
    /// per-species field arrays (conc, flux, etc.), never a name match.
    pub id: SpeciesId,
    /// Whether this species partitions across a gas-liquid interface at
    /// all. `false` ⇒ every function in this module is a documented no-op
    /// for this species (T1) — checked once at registry validation time,
    /// not re-checked per cell/per step.
    pub volatile: bool,
    /// Henry volatility, concentration-based convention (T4), SI units
    /// `Pa*m^3/mol`, at the run's fixed reference temperature. Required
    /// (construction error) when `volatile == true`; absent/unused when
    /// `volatile == false`.
    pub henry_h: Option<f64>,
    /// Liquid-phase molecular diffusivity, m^2/s (SI; converted to
    /// lattice units by the existing unit-conversion layer, REQ §2).
    pub diffusivity_liquid: f64,
    /// Gas-phase molecular diffusivity, m^2/s. Required only on the
    /// resolved-interface branch (W-VOF); point-bubble treats the gas
    /// phase as a well-mixed 0-D inventory per bubble (§1.5) and does not
    /// need a gas-phase spatial diffusivity.
    pub diffusivity_gas: Option<f64>,
}
```

This spec treats `SpeciesTransferProps` as **read-only input** — it is
populated by the registry (whichever spec lands it), validated once at
scenario-construction time (Henry required iff volatile, per T1 — a
**construction-time error**, not a runtime branch, if a volatile species is
missing `henry_h`), and indexed by `SpeciesId` everywhere below. No function
in this module ever matches on a species name or string.

### 2.2 The shared closure library (both branches call these)

```rust
/// Equilibrium liquid concentration in contact with gas at partial
/// pressure `p_gas_pa` (Eq. 1). Returns `None` for a non-volatile
/// species — callers must handle `None` as "no interfacial term", never
/// coerce to 0.0 and proceed (that would silently compute a flux against
/// a fabricated equilibrium of 0, which is a different physical claim
/// from "this species has no interface physics at all").
pub fn equilibrium_conc(species: &SpeciesTransferProps, p_gas_pa: f64) -> Option<f64> {
    species.henry_h.map(|h| p_gas_pa / h)   // Eq. 1: C* = p / H
}

/// Partial pressure implied by a liquid-phase concentration (Eq. 2, the
/// general diagnostic output — called identically for every volatile
/// species; O2/CO2/H2/N2 are registry entries, never named here).
pub fn partial_pressure(species: &SpeciesTransferProps, c_liquid: f64) -> Option<f64> {
    species.henry_h.map(|h| c_liquid * h)   // Eq. 2: p = C * H
}

/// Film mass-transfer coefficient, Higbie penetration (Eq. 8, T6 default).
/// `terminal_velocity_m_s`, `bubble_diameter_m` come from the existing
/// Grace-correlation / d_32 pipeline (T7) — not computed here.
pub fn film_coefficient_higbie(
    diffusivity_liquid_m2_s: f64,
    bubble_diameter_m: f64,
    terminal_velocity_m_s: f64,
) -> f64 {
    let contact_time_s = bubble_diameter_m / terminal_velocity_m_s;
    2.0 * (diffusivity_liquid_m2_s / (std::f64::consts::PI * contact_time_s)).sqrt()
}

/// Film mass-transfer coefficient, Frössling correlation (Eq. 9, T6
/// contaminated-surface cross-check — validity domain §6 risk R1).
pub fn film_coefficient_frossling(
    diffusivity_liquid_m2_s: f64,
    bubble_diameter_m: f64,
    terminal_velocity_m_s: f64,
    kinematic_viscosity_liquid_m2_s: f64,
) -> f64 {
    let re_b = terminal_velocity_m_s * bubble_diameter_m / kinematic_viscosity_liquid_m2_s;
    let sc = kinematic_viscosity_liquid_m2_s / diffusivity_liquid_m2_s;
    let sh = 2.0 + 0.6 * re_b.sqrt() * sc.cbrt();
    sh * diffusivity_liquid_m2_s / bubble_diameter_m
}

/// Interfacial area from Sauter mean diameter + gas holdup (Eq. 10, T7).
/// `d32_m`, `gas_holdup` are consumed from the active interface-tracking
/// mode (W-VOF segmentation or W-BUB PBM) — not produced here.
pub fn interfacial_area_per_volume(gas_holdup: f64, d32_m: f64) -> f64 {
    6.0 * gas_holdup / d32_m
}
```

All four functions are **pure, stateless, unit-explicit** (SI in, SI out;
lattice-unit conversion is the caller's job per the existing unit layer,
REQ §2) — satisfying FR-EXT-01's GPU-evaluability / determinism contract for
extension-point closures without depending on it structurally (this module
has no `Backend` trait dependency; both branches call it as ordinary Rust
functions from whatever context they run in, host or per-cell kernel).

### 2.3 Branch-selection enum (T2) and per-bubble gas inventory (T8)

```rust
/// Which interfacial-transfer branch is active for a run. Mirrors the
/// REQ §1 Interface axis mutual exclusion (resolved-phasefield vs
/// point-bubble); `None` at the top level = no transfer at all (T9).
pub enum TransferBranch {
    /// W-VOF-gated: Eq. 3-6, interface-band source on the phase-wise
    /// scalar system. Requires `phi: Some(_)` (W-VOF active) and
    /// per-phase scalar fields (W-SCAL phase-2 multicomponent + phase-wise
    /// split — see §8 dependency).
    ResolvedInterface,
    /// W-BUB-gated: Eq. 7-9, k_L a film model against each bubble's own
    /// gas inventory. Requires W-BUB's bubble population (Lagrangian
    /// point-bubble set or Eulerian number-density field).
    PointBubble { sherwood_correlation: SherwoodCorrelation },
}

pub enum SherwoodCorrelation { Higbie, Frossling }

/// Per-bubble multicomponent gas inventory (Eq. in §1.5). One instance
/// per tracked bubble in the W-BUB population; owned/stepped by W-BUB's
/// integrator, which calls this module's `bubble_species_flux` each
/// substep — this struct crosses the W-BUB API boundary and is therefore
/// specified here as the hand-off contract (§3.2), not implemented here.
pub struct BubbleGasInventory {
    /// mol of each registry species currently in this bubble's gas core,
    /// indexed by `SpeciesId` (dense array sized to the registry's
    /// volatile-species count; non-volatile species never appear here by
    /// construction — T1 — so this array is never indexed by a
    /// non-volatile SpeciesId).
    pub moles: Vec<f64>,
}

/// The per-bubble, per-species molar transfer rate (Eq. 7 evaluated at
/// this bubble's own partial pressure via Eq. in §1.5), to be
/// integrated by W-BUB's own ODE stepper alongside bubble motion. Pure
/// function of the inputs; does not mutate `inventory` (caller applies
/// dn/dt).
pub fn bubble_species_flux(
    species: &SpeciesTransferProps,
    inventory: &BubbleGasInventory,
    species_id: SpeciesId,
    c_liquid_local: f64,
    bubble_total_pressure_pa: f64,
    k_l: f64,           // from film_coefficient_higbie/frossling (Eq. 8/9)
    bubble_area_m2: f64,
) -> Option<f64> {
    if !species.volatile { return None; }   // T1: no branch, no evaluation
    let x_i = inventory.moles[species_id.index()] / inventory.moles.iter().sum::<f64>();
    let p_i = x_i * bubble_total_pressure_pa;               // Dalton's law
    let c_star = equilibrium_conc(species, p_i)?;           // Eq. 1
    Some(k_l * bubble_area_m2 * (c_star - c_liquid_local))  // Eq. 7, per-bubble scaled
}
```

---

## 3. The two branches: gating and hand-off

### 3.1 Resolved-interface branch — W-VOF-gated

**Preconditions (construction-time check, not a runtime branch):** `phi:
Some(_)` (W-VOF landed and active on this run) AND the phase-wise scalar
split exists (W-SCAL carrying separate `h_gas`/`h_liq` per volatile species —
**this is a W-SCAL phase-2 capability that does not exist yet**; WSCAL_PASSIVE_
SPEC phase 1 is explicitly single-phase, single-component, REQ §3's "phase-wise
conservative" form is future work per that spec's §8). Until that phase-wise
split lands, `TransferBranch::ResolvedInterface` is **rejected at scenario
construction** with an explicit error naming the missing capability (the
FR-VOF-04-style config-guard pattern already used elsewhere in this codebase,
e.g. D-track's α_p threshold guard) — never silently degraded to point-bubble.

**What this spec contributes once the precondition is met**: the interface-band
source terms (6), the `k_int,i` derivation (flagged open, §1.3, §6 risk R2),
and the flux-continuity test (§5 V3). **What W-VOF/W-SCAL contribute**: the
`φ` field, `|∇φ|`, the phase-wise `C_{i,gas}/C_{i,liq}` transport itself.

### 3.2 Point-bubble branch — W-BUB-gated

**Preconditions:** W-BUB's bubble population (Lagrangian point-bubbles or an
Eulerian number-density + `d_32` field) is active. **This spec does not
implement W-BUB** — `WBUB_PBM_IMPL_SPEC.md` does not exist yet (§8). The
hand-off contract:

- W-BUB owns: bubble tracking/PBM (`d_32`, breakup/coalescence per
  Luo-Svendsen/Prince-Blanch, FR-VOF-04), bubble motion (buoyancy/drag/lift
  per REQ §3 point-bubble force balance), the `BubbleGasInventory` **storage**
  (one per tracked bubble/bin) and its **integration** (calling
  `bubble_species_flux` §2.3 each substep, applying `dn_i/dt` to update
  `moles`, applying the resulting bubble-volume change from §1.5's EOS
  closure to its own size/PBM state).
- This spec owns: `equilibrium_conc`, `partial_pressure`,
  `film_coefficient_higbie/frossling`, `interfacial_area_per_volume`,
  `bubble_species_flux` (§2.2-2.3) — the **pure closure functions** W-BUB's
  integrator calls. This spec does NOT own the bubble ODE stepper, the PBM
  kernels, or the liquid-side sink term's placement in the Eulerian scalar
  field's source (`S_{i,liq}^{if} = -Ṡ_i` deposited into the liquid ADE at the
  bubble's location/kernel footprint — W-BUB's masked-source or kernel-
  deposition mechanism, structurally analogous to the D-track particle-
  deposition source pattern already landed, `DISPERSED_DEPOSITION.md`).

**Liquid-side sink placement note**: `Ṡ_i` deposits into the liquid ADE
exactly like the D-track's landed CR-1 interior volumetric source (T18.1) —
this spec's contribution ends at producing the scalar rate; the deposition
kernel (regularized point-source, mass-conserving) is a W-BUB implementation
detail reusing already-landed infrastructure, not a new mechanism this spec
must design.

---

## 4. Outputs

### 4.1 Partial pressure (per volatile species, general — T3)

`p_{X_i}` (§1.2, Eq. 2) is emitted for **every** volatile species in the
registry as a per-cell (resolved branch) or per-bubble/per-cell-averaged
(point-bubble branch) field, named generically `partial_pressure[species_id]`
in the output schema — never `pO2`/`pCO2` as distinct hardcoded fields. The
envelope species (§6) are registry entries whose `partial_pressure` output
happens to be labeled `pO2` etc. **in the registry's own display name**, not
in this module's code.

### 4.2 Composition evolution (point-bubble)

Per-bubble (or per-bin, in an Eulerian PBM binning) `BubbleGasInventory.moles`
time series — the "each bubble's gas composition evolves" output the task
requires — is exposed as a diagnostic time series per bubble/bin, mole
fraction `x_i` and implied bubble `P_bubble, V_b` derived per §1.5. This
reuses the existing per-particle/per-bubble diagnostic output pattern (FR-IO
particle histogram precedent, REQ §4.5) rather than inventing a new output
category.

### 4.3 Mass-transfer rate diagnostics

`k_L a` per species (resolved: interface-integrated Eq. 6 flux ÷ driving
force, reported both ways — as the fitted `k_L a` AND the raw flux integral,
so the two are cross-checkable; point-bubble: `k_{L,i} a` directly from Eq.
8-9/10) is the primary VR-STR-02c/04 comparison quantity and must be emitted
per species, per the same "no hardcoded species" rule.

---

## 5. Validation plan mapped to T17 (VR-STR-04 / VR-STR-02a)

Tests are **authored adversarially by codex/Opus from this spec**, in a
worktree that never shares with the implementation worktree (CLAUDE.md; REQ
§8). All bands below are **provisional MVP gates** (REQ §8 "Band governance"
— tightening always allowed, loosening needs a PHYSICS.md rationale). Every
row states its denominator explicitly (physics-discipline Rule 3 layer 1) and
at least one behavior anchor (layer 2).

| ID | Test | Metric & band (denominator stated) | Grid/setup | T17 row |
|---|---|---|---|---|
| **V1** | **Single-bubble multi-gas absorption/stripping vs film-model prediction.** A single rising bubble (VR-STR-02a geometry: `d_b/Δx≥20`, ρ ratio 100, Grace-regime `U_t` already validated) initialized with a registry gas mixture at non-equilibrium composition (e.g. air-like N2/O2 bubble rising through O2-depleted, CO2-supersaturated liquid — absorbs CO2, strips both O2 and N2 simultaneously). Point-bubble branch. | Per-species: measured cumulative moles transferred vs the **closed-form film-model prediction** `∫k_{L,i}a(C*_i-C_i)dt` integrated with the SAME `k_L,i,a` this spec computes (this is a **self-consistency** check of the numerical integrator against its own closure, denominator = closed-form total transferred moles) — band **< 2%** (numerical-integration-only error, not a physics band). **Behavior anchor**: sign check — the species with `C_i > C*_i` (supersaturated CO2) transfers gas→liquid direction NEGATIVE net for gas inventory... i.e., CO2 net INTO liquid (bubble CO2 moles decrease... unless bubble is a source), O2/N2 net OUT of liquid into bubble (stripping) — assert the SIGN of `dn_i/dt` matches the sign of `(C*_i - C_i)` per species independently. | 3D, single bubble, Grace regime, multi-species registry (envelope §6), CpuScalar | **VR-STR-02a extension (multi-gas)** |
| **V2** | **Henry equilibrium approached at long time.** Same single-bubble (or a well-mixed batch, zero relative velocity, to isolate the film model from hydrodynamics) run to `t → ∞` (many film time constants `1/(k_L a)`). | `\|C_i(t) - C*_i(∞)\| / \|C_i(0) - C*_i(∞)\|` (denominator = initial departure from equilibrium) decays **monotonically** and reaches **< 1%** by `t > 5/(k_L a)` (5 time constants of first-order relaxation — standard criterion) for every volatile species independently. **Behavior anchor**: exponential-decay shape — fit `ln\|C_i(t)-C*_i\|` vs `t`, slope should match `-k_{L,i}a` to **within 5%** (proves the *dynamics*, not just the endpoint, match the film model). | batch / single-bubble, all envelope species, CpuScalar | VR-STR-04 (k_L a formula-vs-measurement gate) |
| **V3** | **Total species mass conserved across the interface, gas loss = liquid gain, per species AND per element.** Closed system (no inflow/outflow), resolved-interface branch (once W-VOF phase-wise split lands, §3.1) OR point-bubble (bubble inventory + liquid field). | (a) **Per-species**: `\|Δn_{i,gas} + Δn_{i,liq}\| / n_{i,total}(0)` (denominator = total initial moles of species i) **< 1e-10** (f64) / **< 1e-6** (f32) — round-off only, the interfacial source (6)/(7) is antisymmetric BY CONSTRUCTION (T5 Eq. 5), so any measured drift is a bug, not a band. (b) **Per-element** (if the registry tracks elemental composition — e.g. carbon in CO2 — a WSCAL_MULTICOMPONENT_SPEC dependency, flagged if absent): same round-off-only band summed over every species containing that element. **Behavior anchor**: gas-phase and liquid-phase inventories move in EXACT lockstep (`Δn_gas(t) = -Δn_liq(t)` at every sampled step, not just at the end). | closed batch, multi-species, CpuScalar (+CpuSimd cross-check per B-6) | VR-STR-05 (species/element conservation) |
| **V4** | **Volatile vs non-volatile behavior anchor.** Registry includes a mix of volatile (envelope §6) and at least one explicitly non-volatile species (glucose, per the task's own example) in the SAME run. | Non-volatile species: interfacial flux **exactly 0.0** (bit-exact, not a small-band check — T1 says the closure is never evaluated) at every step. Volatile species: nonzero flux (sanity — the branch is actually wired). **Behavior anchor**: this is itself the anchor — "zero flux for non-volatile, nonzero for volatile in the identical run" directly tests that the registry `volatile` flag (not a name) gates the physics. | single-bubble or batch, mixed registry, CpuScalar | VR-STR-04 negative/consistency (REQ §8 "sparger phase unit test" analog for species) |
| **V5** | **Partial-pressure output correctness.** For every volatile species, independently perturb `C_i` and confirm `p_{X_i} = C_i · H_i` (Eq. 2) exactly, and confirm `equilibrium_conc(partial_pressure(C)) == C` (round-trip, Eq. 1/2 are exact inverses) to f64 round-off. | exact equality to float round-off (**not a physics band** — this is a unit-test-grade algebraic identity check on the closure functions themselves, §2.2). **Behavior anchor**: none needed (this is a pure-function algebra test, not a simulation behavior test). | unit test, no simulation, all envelope species + at least 2 synthetic species with extreme `H_i` (very soluble / very insoluble) to check no overflow/underflow at the envelope's stated `H_i` range | VR-STR-04 (partial-pressure formula gate) |
| **V6 (ablation, T9)** | **Transfer off → species inert across interface.** `TransferBranch: None` (or `TransferConfig: None` at the top level) on an otherwise-identical single-bubble or batch scenario. | Every volatile species' concentration in BOTH phases is **bit-identical** to a scalar-transport-only run with the SAME species registry but no transfer term wired (the W-SCAL/W-VOF/W-BUB machinery runs, this module contributes literally nothing) — `probe_state_hash` equality, single-backend (VR-STR-05 semantics, B-6 invariance). | single-bubble, CpuScalar | VR-STR-05 (ablation / B-6 invariance) |

**Mandatory negative / consistency tests (REQ §8 pattern, restated for this
spec):**

- **V4** (above) is the non-volatile-zero-flux negative anchor.
- **Sherwood-correlation swap negative test**: substituting Frössling (9) for
  Higbie (8) at a `Re_b` where the two differ by >20% (the R1 risk band, §6)
  must produce a measurably different `k_L a` in V1/V2 — proving the
  correlation choice is load-bearing, not decorative (physics-discipline "if
  it does nothing measurable, delete it").
- **Henry-jump sign negative test**: flipping the sign in `C*_i - C_i` (Eq. 7)
  must make V2's equilibrium-approach test FAIL (diverge from equilibrium
  instead of approaching it) — proves the closure's sign convention is
  load-bearing (mirrors REQ CR-004's forcing-moment negative-test pattern).
- **V6 ablation** (above) is the "transfer off → inert" mandatory negative
  test the task explicitly names.

### 5.6 Behavior-validity review (mandatory, REQ / CLAUDE.md)

After every validation run, before reporting: (a) confirm the **direction** of
each species' net transfer matches its supersaturation sign (`C_i` vs `C*_i`)
independently — a multicomponent run where one species transfers the wrong
direction while the aggregate mass balance still closes (V3) is a
**silently-cancelling-error** pattern (physics-discipline "boundary artifact
sweep") that only V1's per-species sign anchor catches, not V3's total-mass
band; (b) confirm the equilibrium-approach (V2) is monotonic and exponential,
not oscillatory (oscillation would indicate `k_int,i` in the resolved branch
is over the numerical-stability limit — a resolution issue, never to be fixed
by clamping the flux, per the physics-discipline ban list); (c) confirm the
bubble-shrinkage/growth (§1.5, point-bubble) direction matches the net
mole-inventory change sign (a stripping-dominated bubble should grow, an
absorption-dominated bubble should shrink, holding `T,P_hydrostatic` fixed);
(d) record the review in PHYSICS.md per the template in §7.2.

---

## 6. Validity domains + risks

| Closure | Validity domain | Source |
|---|---|---|
| Henry's law (Eq. 1) | Dilute solution (`C_i` well below saturation for all species — Henry's law is the low-concentration limit of the full activity-coefficient equilibrium); isothermal (registry `H_i` fixed at run temperature — no van't Hoff `T`-dependence implemented, §1.1). | Sander 2015 |
| Higbie penetration (Eq. 8) | Clean, freely-circulating (internal-circulation) bubbles; small-to-moderate `Re_b` where surface renewal dominates. Bubble must be in the Grace-diagram spherical/ellipsoidal regime this engine already validates (VR-STR-02a) — Higbie assumes a well-defined, periodically-renewed contact time, which breaks down for cap-regime/highly-deformed bubbles. | Higbie 1935; Kawase & Moo-Young 1990 |
| Frössling (Eq. 9) | Rigid or surfactant-immobilized (contaminated) bubbles, `Re_b` in the boundary-layer range the correlation was fit over (Clift-Grace-Weber cite `Re_b` up to ~1000s for rigid spheres). | Frössling 1938; Clift, Grace & Weber 1978 |
| Sauter-mean area (Eq. 10) | Requires a meaningful `d_32` — i.e., either W-VOF interface segmentation resolving individual bubbles, or W-BUB's PBM producing a converged size distribution. Meaningless (and this spec must reject, not silently proceed) for a single unresolved bubble with no PBM (use the bubble's own `A_b = πd_b²` directly instead, per §1.5's per-bubble form — Eq. 10 is for a *population*, not a single tracked bubble). | Kantarci et al. 2005 |
| Interface-band source `k_int,i` (Eq. 6) | **OPEN — not yet derived** (§1.3). Blocks resolved-branch implementation until a closed form is derived and validated against the sharp-interface limit. | — (flagged, not filled) |

**RISKS:**

- **R1 (high-Re_b Sherwood applicability — the task's named risk).** Higbie
  and Frössling diverge substantially outside their fitted ranges, and real
  bioreactor bubbles occupy a continuum between clean and fully contaminated
  depending on surfactant/protein content (a bioreactor-specific concern —
  cell-culture media is surfactant-rich). Neither correlation has a validated
  closed-form criterion in this codebase for *which* regime a given run is in;
  §2.3's `SherwoodCorrelation` enum makes the choice an explicit scenario
  input (not auto-detected) precisely because auto-detection would require an
  additional un-derived closure (surface-contamination criterion) that this
  spec does not have license to invent (physics-discipline Rule 1: "you would
  have to invent a constant... STOP"). **Mitigation**: both formulas
  implemented and exposed; VR-STR-02a/02c validation runs report which was
  used; a follow-up spec is needed if auto-detection is ever required — not
  attempted here.
- **R2 (interface-band `k_int,i` closed form, resolved branch — open
  derivation, §1.3).** Without it, `TransferBranch::ResolvedInterface` cannot
  be implemented to more than a stub. **This is the STOP-RULE-flagged item of
  this spec** — see the final summary.
- **R3 (isothermal Henry, §1.1).** No `T`-dependence is implemented; a run
  spanning a temperature gradient (not currently a REQ scope item — no landed
  thermal axis) would need the van't Hoff extension before this closure
  applies. Flagged as future work, not a gap in current scope.
- **R4 (W-SCAL phase-wise multicomponent dependency, §3.1).** The resolved
  branch is unimplementable until W-SCAL's phase-wise, multi-species split
  lands — currently WSCAL_PASSIVE_SPEC is single-phase, single-component.
  This is an explicit external dependency, not something this spec can
  resolve.

---

## 7. CODEX ORDER BREAKDOWN

One order = one bundle = one dedicated worktree (CLAUDE.md team convention).
Implementation and adversarial-test orders never share a worktree. Every
physics-affecting order embeds the four lbmflow-codex-dispatch Step 1.5
clauses (reading/ban/stop-rule/two-layer) — restated per-order below as the
concrete band+anchor pair from §5.

| Order | Scope | Primary files (conflict boundary) | Depends on | DoD |
|---|---|---|---|---|
| **O1 — Henry + Sherwood closure library + `p_Xi` output** | `SpeciesTransferProps` (§2.1, against whatever registry lands — see dependency note), `equilibrium_conc`, `partial_pressure`, `film_coefficient_higbie`, `film_coefficient_frossling`, `interfacial_area_per_volume` (§2.2), all pure/stateless/unit-explicit. No simulation wiring. | new module `crates/lbm-core/src/transfer.rs` (new file — no existing-file conflict) | **WSCAL_MULTICOMPONENT_SPEC's registry type** (if unlanded when this order runs, O1 defines `SpeciesTransferProps` standalone per §2.1 and documents the reconciliation TODO) | V5 (partial-pressure algebra) green; unit tests for Eq. 8/9 against literature-tabulated `Sh` numbers at 2-3 reference `Re_b,Sc` points; ban-list grep clean (no species-name branches — the module must compile and pass tests using ONLY synthetic species ids, proving no O2/CO2 string appears in code). |
| **O2 — Point-bubble branch (`k_L a`, multicomponent gas-side balance)** | `TransferBranch::PointBubble`, `BubbleGasInventory`, `bubble_species_flux` (§2.3), the per-bubble EOS volume-coupling (§1.5). Wiring into W-BUB's bubble integrator (call site only — does not implement W-BUB's own ODE stepper). | depends on files W-BUB creates (bubble population storage/integrator) — **coordinate with WBUB_PBM_IMPL_SPEC's file list once it exists**; until then, this order implements `bubble_species_flux` as a standalone pure function (no integration call site) and stops there. | **W-BUB** (hard — no bubble population to attach to without it); O1 (uses its closures) | V1, V2, V4, V6 (§5) green once W-BUB's integrator exists to call into; standalone: O1-style unit tests on `bubble_species_flux` algebra pass. **If W-BUB is not yet landed when this order is dispatched, its DoD is "closure function implemented + unit-tested standalone" — integration DoD deferred, explicitly logged as blocked, not silently skipped.** |
| **O3 — Resolved-interface branch (interface-band source, Henry jump)** | `TransferBranch::ResolvedInterface`, the interface-band source (6), the flux-continuity accounting (3)-(5). **BLOCKED on R2 (§6) — the `k_int,i` closed form is undelivered.** | W-VOF (phase field `φ`) + W-SCAL phase-wise multicomponent split (R4, §6) — **neither landed as of this spec's writing** | O1 | **STOP-RULE applies to this order as currently scoped** — see final summary. Do not dispatch O3 until R2's derivation is delivered by a design session (not a codex implementation order — deriving `k_int,i`'s closed form is a physics-derivation task, not a coding task) AND W-VOF + the W-SCAL phase-wise split have landed. |
| **O4 — Adversarial validation tests** | All of §5 (V1-V6) + the three named negative/consistency tests. Authored from THIS spec, not from O1/O2/O3's implementation. | `crates/lbm-core/tests/xfer_*.rs` (new files only) | O1 (compiles against its closures immediately); O2/O3 (their tests go red→green as those land) | Tests compile and run red against stubs from day one (proves they test the spec, not the implementation); go green incrementally as O1→O2 land; V1-V6 bands frozen into VALIDATION.md T17 as new VR-STR-04/02a sub-rows once O1/O2 are green. O3's tests (V3 resolved-branch variant) stay red/ignored until O3 unblocks. |

**Critical-path ordering:** O1 first (no dependencies, unblocks everything).
O2 can implement its closure algebra immediately (parallel with O1 if bundled
carefully, but simpler to run O1→O2 sequentially given O2 imports O1's
functions) but its integration DoD is genuinely blocked on W-BUB landing. O3
is **not dispatchable yet** — R2 (interface-band closed form) and R4 (W-SCAL
phase-wise split) are both open. O4 runs concurrently from the start in its
own worktree, red until O1/O2 land.

---

## 8. Coexistence / dependency map

- **WSCAL_MULTICOMPONENT_SPEC.md — does not exist yet.** This spec designs
  against the minimal registry surface in §2.1 (`SpeciesTransferProps`:
  `id`, `volatile`, `henry_h`, `diffusivity_liquid/gas`) and asks that
  whichever spec lands the full multicomponent registry (molar mass, reaction
  stoichiometry, elemental composition for V3's per-element conservation
  check, etc.) either uses these exact field names/types or documents an
  explicit mapping. No code in this spec should need to change if the full
  registry is a strict superset of §2.1.
- **WBUB_PBM_IMPL_SPEC.md — does not exist yet.** This spec is written to be
  the closure library WBUB **consumes**: §2.2's four pure functions and §2.3's
  `bubble_species_flux` are the complete surface WBUB needs to call. WBUB owns
  bubble tracking/PBM/motion and the `BubbleGasInventory` storage/integration
  loop (§3.2 hand-off contract). When WBUB_PBM_IMPL_SPEC is written, it should
  cite this spec's §2.2-2.3 as its mass-transfer dependency rather than
  re-deriving Henry/Sherwood closures.
- **WSCAL_PASSIVE_SPEC.md** (landed as a design, not yet implemented per
  repo state at time of writing) — this spec's resolved-interface branch (§3.1)
  needs WSCAL's **phase-2** multicomponent + phase-wise split, which is
  explicitly out of WSCAL's phase-1 scope (that spec's §8, P10: "Multi-component
  ready but phase 1 lands ONE component"). No conflict, just a hard sequencing
  dependency (R4).
- **WVOF_IMPL_SPEC.md** — the resolved-interface branch's `φ`, `|∇φ|` inputs
  come from here; this spec adds no new field to `SoaFields`, so there is no
  structural file-conflict surface with WVOF's O1-order file list (unlike
  WSCAL vs WVOF's shared `fields.rs`/`solver.rs`/`kernels.rs` touches, T3-tier
  transfer is a pure closure module plus BC/source terms that read existing
  fields — see O1/O3 file boundaries above, both new files).
- **DISPERSED_DEPOSITION.md** — the point-bubble liquid-side sink placement
  (§3.2) reuses the landed CR-1 interior-source infrastructure (T18.1) by
  structural analogy; no code sharing required, just a design precedent this
  spec points to instead of re-inventing a deposition mechanism.

---

## PHYSICS.md entries (mandatory, land with O1/O2/O3 respectively)

### O1 entry template

> **Henry's law partition + Sherwood-correlation film coefficients —
> multicomponent interfacial mass-transfer closure library
> (`crates/lbm-core/src/transfer.rs`).** Per-species equilibrium
> `C*_i = p_i/H_i` (concentration-based Henry convention, Sander 2015);
> partial-pressure diagnostic `p_i = C_i H_i` (exact inverse). Film
> coefficient `k_L` from Higbie penetration (default, Eq. 8, Kawase &
> Moo-Young 1990) or Frössling (Eq. 9, Clift-Grace-Weber 1978), caller-selected
> per the surface-mobility regime (no auto-detection — R1). Interfacial area
> `a = 6ε_g/d_32` (Kantarci et al. 2005), consumed from W-VOF/W-BUB, not
> computed here. **Volatile flag from the species registry gates every
> function — non-volatile species never enter a flux computation (not
> clamped to zero after one).** Validity domain: dilute solution, isothermal
> (no van't Hoff `T` dependence — R3); Higbie needs clean/circulating
> bubbles, Frössling needs rigid/contaminated bubbles, and this codebase does
> not auto-select between them (R1, open). Validation:
> `crates/lbm-core/tests/xfer_v5_partial_pressure.rs` (algebra), literature
> `Sh` cross-check unit tests. Replaces/interacts with: nothing prior (new
> closure); feeds W-BUB (§O2) and the resolved-interface branch (§O3, blocked).

### O2 entry template (fill at O2 landing)

> **Point-bubble multicomponent gas-side balance
> (`crates/lbm-core/src/transfer.rs::bubble_species_flux` + W-BUB integrator
> call site).** Each bubble's per-species gas inventory evolves by Eq. 7
> evaluated at the bubble's own Dalton's-law partial pressure (Eq. §1.5); no
> cross-species term (dilute-solution Henry's law has none — resolved
> consequence, not an added approximation). Bubble volume couples to
> `Σn_{i,gas}` via ideal-gas EOS at local hydrostatic+Laplace pressure,
> handed to W-BUB's PBM as an external volume-forcing term. Validation:
> VR-STR-02a multi-gas extension (V1, V2, V4, V6). [fill measured `k_L a`
> values + frozen bands at landing.]

### O3 entry template (do NOT land until R2 is derived — write at that time)

> [BLOCKED — see STOP-RULE below. Entry to be written when the interface-band
> `k_int,i` closed form is derived and the resolved branch is implemented.]

**Literature (decided references):**
Sander 2015 (*Atmos. Chem. Phys.* 15:4399) — compiled Henry's law constants,
convention taxonomy (adopted: concentration-based `H^cp`-inverse). Lewis &
Whitman 1924 (*Ind. Eng. Chem.* 16:1215) — two-film theory (equilibrium at the
interface, resistance in the bulk films — the assumption behind Eq. 7).
Higbie 1935 (*Trans. AIChE* 31:365) — penetration theory, Eq. 8. Frössling
1938 (*Gerlands Beitr. Geophys.* 52:170); Clift, Grace & Weber 1978, *Bubbles,
Drops, and Particles* Ch. 3 — Eq. 9. Kawase & Moo-Young 1990 (*Chem. Eng. Sci.*
45:1435) — penetration-model review for bubble columns/bioreactors. Kantarci,
Borak & Ulgen 2005 (*Process Biochem.* 40:2263) — bubble-column mass-transfer
review, interfacial-area relation Eq. 10. Garcia-Ochoa & Gomez 2009
(*Biotechnol. Adv.* 27:153) — stirred-tank `k_L a` correlation review (context
for VR-STR-02c). Grace (Grace, Wairegi & Nguyen 1976) — terminal-velocity
diagram already adopted by VR-STR-02a, reused here for `U_t` in Eq. 8/9.

---

**STOP-RULE (this spec's honest disposition on the resolved-interface
branch):**

```
STOP-RULE: T17 gate for TransferBranch::ResolvedInterface (§O3) is
unreachable without either an undelivered closed-form derivation or a
banned hack.
Attempted: derived the sharp-interface flux-continuity + Henry-jump
conditions (Eq. 3-5) fully from first principles; attempted to write the
diffuse-interface realization (Eq. 6) but the interfacial exchange rate
k_int,i has no closed form yet relating it to (D_gas, D_liq, W) that
recovers the sharp limit as Delta x / W -> 0 — filling it with a fitted
number would be exactly the banned "calibrated constant" pattern.
Blocking mechanism: a diffuse-interface (phase-field) representation of a
sharp Henry jump condition requires an asymptotic matched-inner-outer
derivation (the same class of derivation Allen-Cahn's mobility M required)
that has not been done for the two-species-diffusivity jump case.
Options for the PM: (a) commission a dedicated derivation task (physics
design session, not a codex coding order) before dispatching O3; (b) accept
point-bubble-only interfacial transfer as the initial delivery and defer
resolved-interface transfer to a later phase (matches REQ's own "Phase 2
API-reserved" framing for point-bubble machinery, inverted here since
point-bubble transfer is actually the SIMPLER path to land first); (c) no
literature closure found yet for the diffuse-interface Henry-jump rate
constant specifically (adjacent Allen-Cahn mobility literature does not
directly transfer because Henry's law is a concentration JUMP, not a
continuous field, unlike phase-field density/viscosity blending).
```

This report is a SUCCESS outcome of this spec's O3 section, not a failure —
per `.claude/skills/lbmflow-physics-discipline` Rule 4.
