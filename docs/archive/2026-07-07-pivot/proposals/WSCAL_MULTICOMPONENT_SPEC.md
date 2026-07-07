# W-SCAL Phase 2 Implementation Specification — Multi-Species ADE with SGS Scalar Flux and Two-Phase Partition

**Document ID**: SPEC-W-SCAL-MULTICOMPONENT (rev.1, 2026-07-07).
**Scope**: the phase-2 generalization named by `docs/proposals/WSCAL_PASSIVE_SPEC.md`
decision **P10** ("Multi-component ready but phase 1 lands ONE component… a
`Vec<ScalarField>` generalization is a phase-2 refactor, API-reserved but not
built") and by that spec's decision **P9**/§4.5 hooks (`F_b^scalar` reserved,
`ν_t/Sc_t` SGS scalar flux OUT of phase 1). This spec is the **T2 tier of the
reaction-engineering goal**: N independently-diffusing, independently-phase-
partitioning chemical species (O2, CO2, H2, N2, substrate, product, ionic
species, …) transported through the same D3Q7 ADE machinery WSCAL_PASSIVE
landed for one component, with the REQ §3 phase-wise conservative form and the
FR-LES-04 `Sc_t` hook now implementable (W-LES/WALE `ν_t` is landed).
**Target core**: `crates/lbm-core` — extends `SoaFields` with an `N`-slot
scalar-field vector; touches `lattice.rs` (reuses D3Q7, no new lattice),
`fields.rs`, `solver.rs`, `kernels.rs`.
**Acceptance**: VALIDATION.md **T17** row **VR-STR-04** (scalar/reaction,
multi-species Taylor–Aris extension) and **VR-STR-05** (total-species-mass +
element-balance conservation), REQ §3 "two-phase phase-wise conservative"
governing form, REQ §4.2 FR-LES-04 (`Sc_t`), REQ §7 NFR-01 (memory budget:
`N × 56 B/cell` is additive per the existing budget row — no new row needed).

This spec is **executable**: it does not re-derive any WSCAL_PASSIVE decision
(D3Q7 lattice, `τ_s` mapping, BGK/TRT collision, halo path, wall closures) —
it generalizes their *scope* from one scalar to `N`, adds the SGS flux term
the parent left OUT, and adds the two-phase partition form the parent's
single-phase equation is the special case of. Every new decision below is
either resolved from the governing equations (§1), a literature-backed closure
with citation + derivation + validity domain + its own test (Rule 1 of
`.claude/skills/lbmflow-physics-discipline`), or flagged STOP-RULE (none are,
in this document — see §5 negative-test list and §7 for the one deferred
item that is out of scope rather than blocked).

> **Discipline note (CLAUDE.md prime directive; PHYSICS.md;
> `.claude/skills/lbmflow-physics-discipline`).** No species name is ever
> compiled into a branch (O2/CO2/H2/N2/substrate/product/ion are configuration
> instances of one generic `Species` record, §2). No band-calibrated constant
> appears; `Sc_t = 0.7` is the REQ §3-cited literature default and is treated
> as a closure parameter, not a fitted constant (§1.4 gives its provenance and
> validity domain). The mandatory PHYSICS.md entry text is in §7.4; the
> behavior-validity review checklist is in §5.6, extending WSCAL_PASSIVE §5.6
> item-for-item to N species and to the Sc_t ablation.

---

## 0. Summary of decisions (read this first)

| # | Decision | Justification (short) |
|---|---|---|
| M1 | **`SoaFields` carries `scalars: Vec<ScalarField<T>>`** (length `N`, `N=0` legal and bit-identical to today), replacing the single `Option<h>/Option<htmp>/Option<conc>` triplet WSCAL_PASSIVE reserved. Each `ScalarField<T>` owns its own `h: Vec<T>`, `htmp: Vec<T>`, `conc: Vec<T>`, `tau_s: T` (uniform) or `Option<Vec<T>> omega_s_field` (per-cell), and a `SpeciesId` back-reference into the registry (§2). | A `Vec<ScalarField>` is the exact generalization WSCAL_PASSIVE P10 named. Keeping `N=0 ⇒ empty Vec` preserves the B-6 `None`-is-bit-identical discipline losslessly (an empty `Vec` allocates nothing, iterates zero times, byte-identical to the pre-feature struct layout modulo the `Vec`'s own header — verified by the V5 bit-identity gate, §5). |
| M2 | **Species registry is a host-side, `Copy`-free struct `Species { name: String, henry_h: Option<f64>, diffusivity_d: f64, charge_z: f64, molecular_weight_mw: f64, volatile: bool }`**, stored once per scenario (`Vec<Species>`, index = `SpeciesId(usize)`), NOT duplicated per cell. Core code reads only `diffusivity_d` (→ `τ_s`) and `volatile`/`henry_h` (→ interface BC selection, phase-2-of-phase-2, hook only). `charge_z`/`molecular_weight_mw` are carried for future electroneutrality / element-balance bookkeeping (§1.6) and are NOT read by any collision/streaming kernel in this spec — pure metadata passthrough. | Component-agnostic by construction: the registry is data, the kernels are generic over `N` and read only per-species *numbers* (`D_i`, `τ_{s,i}`), never `name`. This is the literal ban-list check for "case-identity branch" (Rule 2 of the physics-discipline skill) — grepping the diff for `if name ==` / `match name` must return zero hits in any kernel path. |
| M3 | **One D3Q7 `h`-set per species, independently relaxed.** `τ_{s,i} = D_i/cs_s² + ½` per WSCAL_PASSIVE (6), unchanged formula, applied `N` times with `N` different `D_i`. No new lattice, no new equilibrium form — §1.2/§1.3/§1.4 of WSCAL_PASSIVE apply verbatim per-species. | Species do not interact through the ADE operator itself in this spec (no cross-diffusion / Stefan-Maxwell coupling — REQ does not ask for it, and adding it now would be an uncited closure; §6.6 records this as an explicit non-goal, not a silent omission). Reaction coupling (`R_k(C)`, W-REACT) is a separate DAG node and out of scope. |
| M4 | **Two-phase phase-wise conservative form is the governing equation when a phase field `φ` is present** (REQ §3 verbatim); **single-phase passive is the `φ≡1` (or W-VOF absent) special case**, matching WSCAL_PASSIVE's own framing of its equation (2) as "the equivalent non-conservative form" of the general one. The phase-wise form is implemented as a **`φ`-weighted equilibrium and a `φ`-weighted zeroth moment**, NOT a separate code path — §1.5 gives the single formula that degenerates correctly. | REQ §3: "non-conservative single-phase is a special case" of the two-phase conservative form (REV-CFD-MJ-011) — the spec must implement ONE equation, not two, or the single-phase path silently diverges from the reference when W-VOF later lands (a coexistence bug, not a physics bug, but one this spec must not create). |
| M5 | **Per-species, per-phase diffusivity: `D_{i,q}` (species × phase), NOT a single `D_i`.** When `φ` is present, `τ_{s,i}` becomes **per-cell** via the existing `omega_s_field` mechanism (WSCAL_PASSIVE §4.5 hook, generalized to per-species): `D_i(x) = D_{i,gas} + φ(x)(D_{i,liq} − D_{i,gas})` (linear-in-φ interpolation of the *diffusivity itself* — see §1.5 for why this, not harmonic, is the correct interpolation for a Fickian flux, contrasting WVOF's harmonic-in-μ choice for viscosity). | O2/CO2/H2/N2 have markedly different `D_liq`/`D_gas` ratios (gas-phase diffusivities are ~10⁴× liquid-phase for typical small molecules — Cussler 2009 Table 3.1-1); a single `D_i` cannot represent aeration transport. This is the concrete instance of REQ §3's "per-species diffusivity may differ by phase" (task envelope). |
| M6 | **SGS scalar flux: `D_eff,i = D_i + ν_t/Sc_t`, per species, `Sc_t = 0.7` default**, read from the ALREADY-LANDED `WaleLes::nu_t()` field (`les.rs:103`, "global compact order") the same way `omega_field` (hydrodynamic) is installed by `set_omega_field` (`solver.rs:2663`). A new `set_omega_s_field(species_idx, omega_s)` follows the identical precedent, per-species. `Sc_t` is a **per-species scenario parameter** (registry-adjacent, not in `Species` itself — it is a closure parameter of the *turbulence-scalar coupling*, not a species property; see §2.3), defaulting to 0.7. | REQ §4.2 FR-LES-04 mandate, WSCAL_PASSIVE §4.5's named "first phase-2 add" now that W-LES is landed. `Sc_t = 0.7` is the standard literature default for gas/liquid turbulent mass transfer (Reynolds analogy value; Launder & Spalding 1974; Tominaga & Stathopoulos 2007 review report the 0.2–1.3 measured range with 0.7–0.9 the common CFD default) — cited, not fitted; validity domain and the mandatory ablation guard test are in §5. |
| M7 | **Total-species-mass conservation is per-species-per-phase** (the REQ §3 conservation statement `Σ_q ∫ α_q C_{k,q} dV` changes only by boundary flux/reaction), reducing to WSCAL_PASSIVE V4's single-phase zeroth-moment sum when `φ≡1`. **Element-balance conservation is a host-side linear-combination check** (`Σ_k ν_{k,e} × (species k total mass)` for a user-declared stoichiometric element map `e`), NOT a new core invariant — no reaction is implemented in this spec (R(C)=0), so element balance reduces to "every species individually conserves its own mass," a corollary of M7's per-species gate, not a new mechanism. | No reaction path exists yet (W-REACT is a separate, later DAG node per REQ §11); claiming a "real" element-balance gate before a reaction term exists would be Rule-2-banned decorative physics (a test that always trivially passes with `R=0` and validates nothing new). §5 V9 states this precisely and is written to FAIL if a future reaction term breaks per-species conservation without a compensating element map — a forward-guard, not a current-value gate. |
| M8 | **`N=0` (empty `Vec<ScalarField>`) is bit-identical to the pre-multicomponent engine, including bit-identical to WSCAL_PASSIVE's own `N=1` machinery when `N=1` and the phase-1 code path is exercised** — i.e. this spec's generalization must not perturb WSCAL_PASSIVE's own V1–V5 gates. Verified by re-running WSCAL_PASSIVE's exact V1/V2/V4/V5 tests unmodified against the generalized `Vec`-backed storage. | B-6 invariance, applied twice: once at `N=0` (no scalars at all) and once at `N=1` (this spec must be a strict superset of WSCAL_PASSIVE's behavior, not a rewrite that happens to also work for one species). |
| M9 | **No cross-species interaction in the collision operator** (M3) and **no interfacial mass-transfer source `S^if`** (Henry-partition flux term) in this spec — both are explicitly deferred to W-BUB / interfacial-transfer (REQ §11: "waits on W-VOF" for the resolved-interface flux, separate DAG node from this generalization). The registry's `henry_h: Option<f64>` field is carried **as inert metadata** for that later phase; nothing reads it in this spec's kernels. | Scope discipline (CLAUDE.md "minimal scope" + Rule 1: `S^if` is a literature closure requiring the Henry-partition jump condition at the resolved interface, which needs `φ` gradient/normal reconstruction machinery this spec does not touch). Grep-checkable: `henry_h` appears in `Species` struct and nowhere else in `crates/lbm-core/src`. |

---

## 1. Governing equations + generalized LBE

### 1.1 What WSCAL_PASSIVE already resolved (unchanged, cited not re-derived)

WSCAL_PASSIVE §1.1–§1.4 fixed, for ONE scalar:
- the ADE `∂C/∂t + ∇·(Cu) = ∇·(D∇C)` (its Eq. 1/2);
- the D3Q7 `h_i` distribution, linear-in-`u` equilibrium `h_i^eq = w_i^s C[1 + c_i·u/cs_s²]` (its Eq. 3);
- BGK (Eq. 4) and TRT with `Λ=1/4` (Eq. 5) collision;
- the diffusivity mapping `τ_s = D/cs_s² + ½ = 4D + ½` (its Eq. 6, `cs_s²=1/4`);
- the D3Q7 lattice definition, weights `w_0=1/4, w_{1..6}=1/8` (its §1.4);
- BCs: bounce-back no-flux, anti-bounce-back Dirichlet, zero-gradient outflow (its §2).

**Every one of these is reused verbatim, per species `i = 0..N-1`, with
`D → D_i`, `τ_s → τ_{s,i}`, `C → C_i`, `h → h^{(i)}`.** This spec adds nothing
to the per-species ADE-LBM mechanics; it adds (a) the `N`-indexing, (b) the
SGS flux term, (c) the two-phase partition weighting.

### 1.2 REQ §3 governing forms — the three-tier hierarchy this spec must implement as ONE code path

REQ §3 gives three forms (verbatim, reordered here from general to the
special case WSCAL_PASSIVE already implements):

```
Density-based active (out of scope — W-REACT/active feedback, not this spec):
   ∂(ρY_k)/∂t + ∇·(ρ u Y_k + J_k) = R_k + S_k^if

Two-phase phase-wise conservative (q ∈ {gas, liquid}, α_liq=φ, α_gas=1−φ):
   ∂(α_q C_{k,q})/∂t + ∇·(α_q u C_{k,q})
       = ∇·[α_q(D_{k,q} + ν_t/Sc_t)∇C_{k,q}] + α_q R_{k,q} + S_{k,q}^if      (REQ-3)

Single-phase passive (ρ, α uniform; α_q ≡ 1, R_k ≡ 0, S^if ≡ 0):
   ∂C_k/∂t + u·∇C_k = ∇·[(D_k + ν_t/Sc_t)∇C_k]                              (REQ-3′)
```

`REQ-3′` is exactly WSCAL_PASSIVE's Eq. (2) *plus* the `ν_t/Sc_t` term that
spec explicitly left OUT (its §4.5). This spec's job is: (a) add the `ν_t/Sc_t`
term to the single-phase form (§1.4 below — applies even with `N` species and
no `φ`), and (b) implement `REQ-3` for when `φ` is present, such that setting
`α_q ≡ 1` in the `REQ-3` implementation recovers `REQ-3′` bit-for-bit (not
approximately — see §1.5's degeneration argument).

### 1.3 Per-species D3Q7 LBE (mechanical generalization of WSCAL_PASSIVE Eq. 3–5)

For species `i`, phase weight `α(x,t)` (defined in §1.5; `α≡1` if no `φ`
field exists), local diffusivity `D_i(x,t)` (§1.5), local SGS-augmented
diffusivity `D_{eff,i}(x,t) = D_i(x,t) + ν_t(x,t)/Sc_{t,i}` (§1.4):

```
h_i^{eq,(k)} = w_i^s α(x) C_k [ 1 + (c_i · u)/cs_s² ]                       (1)
τ_{s,k}(x)   = D_{eff,k}(x) / cs_s² + ½ = 4 D_{eff,k}(x) + ½                (2)
```

BGK/TRT collision (WSCAL_PASSIVE Eq. 4/5) applied per-`k`, with `τ_{s,k}(x)`
now **per-cell** (via `omega_s_field`, §1.4) whenever `ν_t≠0` or `φ` varies,
and **uniform** (the WSCAL_PASSIVE fast path, unchanged) when neither varies —
these are the same code path with the field either `None` (uniform scalar
read) or `Some` (per-cell gather), exactly the `omega_field` precedent
(`fields.rs:199`). `C_k = Σ_i h_i^{(k)} / α(x)` when `α` is carried in the
equilibrium (§1.5 explains why the zeroth moment must be divided by `α`, not
multiplied into the transported quantity another way).

### 1.4 SGS scalar flux — the `ν_t/Sc_t` closure (decision M6)

**Closure**: `D_{eff,k}(x) = D_k(x) + ν_t(x)/Sc_{t,k}`, `Sc_{t,k}` a per-species
scenario scalar defaulting to **0.7**.

**Source / derivation**: this is the standard gradient-diffusion SGS scalar
flux closure, the scalar analog of the eddy-viscosity momentum closure already
landed for WALE (`les.rs`): filtered scalar transport produces an unresolved
subgrid flux `⟨u'C'⟩` which is modeled as `−(ν_t/Sc_t)∇C̄` (Reynolds-analogy
gradient closure; Pope 2000 *Turbulent Flows* §10.4; the `Sc_t` value itself
is an empirical closure constant with a well-characterized range — Tominaga &
Stathopoulos (2007) review measured `Sc_t ∈ [0.2, 1.3]` across
building/environmental-flow experiments, with **0.7–0.9 the standard CFD
default** and REQ §3 naming **0.7** explicitly as "default." This spec adopts
the REQ-mandated value; it is a **cited literature default, not a fitted
constant** — Rule 1 row 2 (literature-backed closure) applies, and the four
required artifacts are: (1) this citation+derivation, (2) validity domain
below, (3) the ablation-guard test V8 (§5), (4) the PHYSICS.md entry (§7.4).

**Validity domain**: valid wherever the underlying WALE `ν_t` is valid (LES
resolved-turbulence regime, `y+/cell` resolution caveats already recorded in
PHYSICS.md for WALE); `Sc_t` is treated as constant-in-space (no dynamic
`Sc_t` model) — a stated model limitation, not a hidden approximation. If LES
is inactive (`ν_t≡0` everywhere, no `WaleLes` driver installed), `D_eff,k ≡
D_k` exactly — this is WSCAL_PASSIVE's own required behavior (its §4.5 "phase
1 must not silently include `ν_t`") and is now the **general** case's `ν_t=0`
degeneration, verified by the same negative test WSCAL_PASSIVE specified,
generalized to N species (§5 V8b).

**Timing note (one-step lag, inherited from WALE, not new)**: `les.rs:5-7`
documents that WALE's `nu_t`/`omega_plus` field "has a one-step lag" (computed
from the current velocity gradient, applied to the *next* collision). The
scalar `D_eff,k` inherits the identical lag — `set_omega_s_field` reads the
same `nu_t()` snapshot `set_omega_field` used for the hydrodynamic step in the
same solver iteration, so hydrodynamic and scalar relaxation see the *same*
lagged `ν_t`, preserving internal consistency (no new lag is introduced
relative to the existing hydrodynamic LES lag).

**Field mechanics**: extends the `omega_field` precedent (`fields.rs:199`,
`set_omega_field` `solver.rs:2663`) to a **per-species** `omega_s_field`
living inside each `ScalarField<T>` (§3.2), `None` when `D_{eff,k}` is uniform
(no LES, no `φ`), `Some(Vec<T>)` (compact-core) otherwise. The solver-level
call sequence, once per step, per species with LES active:

```
nu_t = wale.nu_t();                              // from this step's velocity grad
for each species k:
    for each cell x:  tau_s_k[x] = 4.0*(D_k(x) + nu_t[x]/Sc_t_k) + 0.5;
    solver.set_omega_s_field(k, Some(&tau_s_k));  // installs before h-collide(k)
```

### 1.5 Two-phase phase-wise weighting (decision M4, M5) — the single-formula degeneration argument

REQ-3's `α_q C_{k,q}` product is the transported "phase-partial concentration"
(REQ's own term for `α_liq × C_{k,liq}` or `α_gas × C_{k,gas}`, summed
implicitly by the fact that a species k crossing the interface is represented
by ONE `C_k(x)` field weighted by whichever phase occupies `x`). This spec
implements it as: the transported LBM moment is `α(x) C_k(x)`, i.e.

```
h_i^{(k),eq} = w_i^s [α(x) C_k(x)] [ 1 + c_i·u/cs_s² ],   α(x) = φ(x)         (3)
Σ_i h_i^{(k)} = α(x) C_k(x)   ⇒   C_k(x) = (Σ_i h_i^{(k)}) / α(x)             (4)
```

**Degeneration to REQ-3′ (single-phase, decision M4 requirement)**: when no
`φ` field exists (W-VOF not mounted, or the scenario declares single-phase),
`α(x) ≡ 1` identically (a compile-time/config-time constant, not a runtime
`φ=1` field read — no division-by-`α` cost or hazard in the single-phase
path), and (3)/(4) reduce to WSCAL_PASSIVE's Eq. (3)/`C=Σh_i` **exactly,
bit-for-bit**, because `α≡1` makes the multiply/divide no-ops that the
compiler/kernel author elides via the SAME `Option`-gated branch `omega_field`
already uses (`α_field: Option<Vec<T>>`, `None` ⇒ the uniform-`α=1` fast path,
mirroring M1's storage decision). This is the concrete mechanism that
satisfies M8 (bit-identical to WSCAL_PASSIVE at `N=1`, no-`φ`).

**Per-phase diffusivity interpolation — linear-in-D, not harmonic (decision
M5)**: `D_k(x) = D_{k,gas} + φ(x)(D_{k,liq} − D_{k,gas})`. **Why linear, not
harmonic** (contrasting WVOF_IMPL_SPEC §1.3's harmonic-in-μ default for
viscosity): viscosity's harmonic mixing rule is derived from continuity of
*shear stress* across a diffuse interface (momentum-flux matching — the
correct closure for a mechanical modulus in series). Fickian diffusive flux
`J = −D∇C` has no analogous series-resistance structure *within* the diffuse
band when the model already treats `C_k` as one continuous field advected by
one continuous `α` weight (the resistances-in-series picture applies to a
*sharp* two-film interface model, which is the point-bubble `k_L a` closure
REQ §3 lists SEPARATELY as the `S^if` term this spec does NOT implement — see
M9). For the resolved-diffuse-interface, phase-wise-conservative form REQ-3
prescribes, the diffusivity that appears inside `∇·[α_q D_{k,q} ∇C_{k,q}]` is
the *local* material diffusivity of whichever phase occupies each point, and
volume-fraction-linear interpolation is the standard diffuse-interface mixing
rule for a transport coefficient that is NOT a series modulus (Kim 2012,
*Phase-Field Models for Multi-Component Fluid Flows*, reviews linear mixing
for scalar diffusivity in Cahn-Hilliard-coupled transport; contrast his
harmonic treatment of viscosity in the same reference, which WVOF already
cites). **Validity domain**: interface band width `W` (WVOF's parameter) must
resolve the transition — this spec inherits WVOF's `W` sizing, adds no new
resolution requirement. If W-VOF is not mounted, this paragraph is inert
(§1.5's first paragraph's degeneration applies).

### 1.6 Element balance (decision M7) — explicitly NOT a new mechanism

REQ does not define a reaction operator in this spec's scope (`R_k ≡ 0`
throughout — M9). Therefore "element balance" reduces to: for a user-declared
map `e: SpeciesId → stoichiometric coefficient` (e.g. carbon-count for a
CO2/glucose/product set), `Σ_k e(k) × (total mass of species k)` is conserved
**iff** each species individually conserves its own total mass (M7's
per-species gate) — linearity, not a new physical mechanism. This spec
implements the per-species total-mass gate (§5 V6/V7) and a **host-side
utility** (no new kernel) that evaluates the linear combination for a
caller-supplied `e` map, purely as a convenience/regression tool for the
future reaction phase — it asserts nothing that the per-species gate doesn't
already imply, and its test (§5 V9) is written to catch the day a reaction
term is added without updating this bookkeeping (a forward guard, per M7's
justification column).

---

## 2. Species registry data structure (decision M1, M2, M6)

### 2.1 The `Species` record (host-side, scenario/config layer)

```rust
/// One chemical species in the multi-component scalar registry. Pure
/// configuration data — no kernel reads `name`; kernels read only the
/// numeric fields resolved into per-species LBM parameters (see
/// `ScalarField::tau_s` / `omega_s_field`, §3.2).
#[derive(Clone, Debug, PartialEq)]
pub struct Species {
    /// Human-readable label (diagnostics/output only — e.g. "O2", "CO2",
    /// "substrate"). NEVER matched on by any transport/collision code path
    /// (banned case-identity branch, physics-discipline Rule 2).
    pub name: String,
    /// Henry's law constant (dimensionless or [Pa·m³/mol] per the scenario's
    /// declared convention — unit convention frozen at the O2 scenario-schema
    /// order, not here). `None` ⇒ non-volatile (e.g. glucose, ionic species
    /// with no gas-phase presence). Inert metadata in THIS spec (M9) —
    /// reserved for the interfacial-transfer phase.
    pub henry_h: Option<f64>,
    /// Molecular diffusivity in the LIQUID phase, physical units [m^2/s].
    /// Always required (every species has a liquid-phase diffusivity even if
    /// it never appears in the gas phase).
    pub diffusivity_d_liquid: f64,
    /// Molecular diffusivity in the GAS phase, physical units [m^2/s].
    /// Required when `henry_h.is_some()` (the species can occupy gas);
    /// otherwise ignored (never read — `volatile=false` species use
    /// `diffusivity_d_liquid` unconditionally via `α≡1`-equivalent gas-phase
    /// masking, §3.2).
    pub diffusivity_d_gas: f64,
    /// Ionic charge (dimensionless, signed). Metadata only in this spec — no
    /// electromigration/electroneutrality term is implemented (would be an
    /// uncited Nernst-Planck closure; out of scope, not a hook-in-waiting for
    /// this document since REQ does not request it — recorded here only
    /// because the task envelope names "ion" as an example species).
    pub charge_z: f64,
    /// Molecular weight [g/mol]. Metadata only — used by the host-side
    /// element-balance utility (§1.6) and by scenario-layer unit conversions
    /// (e.g. mol/L ↔ mass/L for `conc`), never by a core kernel.
    pub molecular_weight_mw: f64,
    /// Whether this species can occupy the gas phase (drives whether a
    /// gas-phase diffusivity / Henry partition is meaningful). `false` ⇒
    /// `diffusivity_d_gas` is never read; `D_k(x) ≡ diffusivity_d_liquid`
    /// unconditionally (the species is liquid-confined by definition, not by
    /// a runtime branch on `name`).
    pub volatile: bool,
}
```

**Registry storage**: `Vec<Species>` held once per scenario (scenario/config
layer, `crates/lbm-scenario`), NOT per-cell, NOT in `SoaFields`. Core
(`lbm-core`) never sees a `Species`; it sees only the resolved numeric
products (`D_i` per phase, `τ_{s,i}` or `omega_s_field`) that the scenario
layer computes ONCE at setup and hands to `Solver` — identical in spirit to
how WSCAL_PASSIVE's single `D` becomes `τ_s` at config time, generalized to a
per-index array. `SpeciesId(usize)` is the index into both the registry
`Vec<Species>` and the `SoaFields.scalars: Vec<ScalarField<T>>` — the two
`Vec`s are **parallel arrays**, same length `N`, same order; core code never
needs the registry to run, only `N` and the per-index `D`/`τ_s` it was
configured with.

### 2.2 Rust API sketch — solver-facing (core, generic over N, no species names)

```rust
// crates/lbm-core/src/fields.rs — SoaFields<T> addition (replaces the
// WSCAL_PASSIVE single h/htmp/conc Option triplet with a Vec):

/// One species' D3Q7 ADE state + relaxation parameters, compact/padded per
/// the same layout convention as `f`/`ftmp` (§3.1). `SoaFields::scalars` is
/// `Vec::new()` (length 0) on the legacy/no-scalar path — bit-identical to
/// pre-W-SCAL (decision M1, M8).
#[derive(Clone, Debug)]
pub struct ScalarField<T: Real> {
    /// D3Q7 populations, q-major padded planes (current state). Same layout
    /// as `f` (§3.1 of WSCAL_PASSIVE), lattice fixed to D3Q7 regardless of
    /// the hydrodynamic lattice.
    pub h: Vec<T>,
    /// Ping-pong partner of `h`.
    pub htmp: Vec<T>,
    /// Macroscopic concentration C_k = (Sum_i h_i) / alpha, compact core.
    /// alpha = 1 (no phase field) is the fast path; see `fields.rs` doc on
    /// `SoaFields::phase_field` for where alpha is read from when present.
    pub conc: Vec<T>,
    /// Uniform relaxation rate `tau_s` for this species when both `D_k` and
    /// `Sc_t`-scaled `nu_t` are spatially uniform (i.e. `omega_s_field` is
    /// `None`). Always populated (even when `omega_s_field` is `Some`, this
    /// holds the last-uniform value for diagnostics) — mirrors how the
    /// hydrodynamic uniform `omega` and `omega_field` coexist.
    pub tau_s: T,
    /// Optional per-cell relaxation rate override (`omega_s = 1/tau_s(x)`),
    /// compact core. `None` preserves the uniform-`tau_s` collide path
    /// exactly (same precedent as `SoaFields::omega_field`, `fields.rs:199`).
    /// Populated whenever LES SGS flux is active for this species (`D_eff,k`
    /// varies in space via `nu_t(x)`) or a two-phase `D_k(x)` varies (M5).
    pub omega_s_field: Option<Vec<T>>,
    /// Registry back-reference (index into the scenario's `Vec<Species>`,
    /// §2.1). Used only for output/diagnostics labeling — never for kernel
    /// branching (no kernel reads this field's target `Species.name`).
    pub species_id: usize,
}

// SoaFields<T> gains:
pub struct SoaFields<T: Real> {
    // ...unchanged fields (f, ftmp, rho, ux, uy, uz, solid, ...)...
    /// Multi-species passive-scalar ADE state, one entry per registered
    /// species, `N = scalars.len()`. Empty on the scalar-free path — no
    /// allocation, bit-identical to pre-W-SCAL (decision M1).
    pub scalars: Vec<ScalarField<T>>,
    /// Phase-fraction field `alpha(x) = phi(x)` read (not owned) from the
    /// W-VOF `g`-path's `phi` field when present; `None` when W-VOF is not
    /// mounted or the scenario is declared single-phase (decision M4 fast
    /// path — `alpha == 1` identically, no per-cell read). This is a
    /// *reference*-shaped hook (the actual `phi` storage lives in the W-VOF
    /// `Option<Vec<T>>` slot per WVOF_IMPL_SPEC §3.2); W-SCAL-multi does not
    /// duplicate it.
    pub scalar_phase_weight: Option<()>,  // placeholder marker; see O2 note in §7 — the
                                            // concrete wiring (borrow vs copy of `phi`)
                                            // is an O2-time decision, not re-litigated here.
}
```

```rust
// crates/lbm-core/src/solver.rs — Solver<L,B,T> additions (mirrors
// set_omega_field, solver.rs:2663, one call per species):

impl<L: Lattice, B: Backend<L, T>, T: Real> Solver<L, B, T> {
    /// Register a new species' D3Q7 scalar field, allocated quiescent
    /// (C=0). Returns its `SpeciesId` (== index into `scalars`). Must be
    /// called before the first `step()`/`run_span()` — the scalar sub-step
    /// (§4) iterates `0..scalars.len()` every step, so adding a species
    /// mid-run is unsupported (a config-time operation, not a runtime one).
    pub fn add_species(&mut self, diffusivity_d: T, /* + phase-aware variant, §2.3 */) -> usize;

    /// Set/replace the uniform relaxation rate for species `k` (no LES, no
    /// phase weighting) — the fast path.
    pub fn set_species_tau(&mut self, k: usize, tau_s: T);

    /// Install/clear a per-cell relaxation-rate override for species `k`
    /// (`omega_s_field`), mirroring `set_omega_field`'s `Option`-swap
    /// mechanics exactly (`solver.rs:2663-2690`). Called once per step from
    /// the LES/phase-weighting pre-pass (§1.4) when either varies spatially.
    pub fn set_omega_s_field(&mut self, k: usize, omega_s: Option<&[T]>);

    /// Set a species' boundary condition per global face (bounce-back /
    /// anti-bounce-back Dirichlet with C_in / zero-gradient) — the
    /// WSCAL_PASSIVE §2 BC menu, now indexed by species (a wall can be
    /// no-flux for O2 and simultaneously Dirichlet for a tracer dye, e.g.).
    pub fn set_species_face_bc(&mut self, k: usize, face: Face, bc: ScalarFaceBc<T>);
}
```

### 2.3 `Sc_t` storage — a per-species scenario parameter, not a `Species` field (decision M6)

`Sc_{t,k}` is stored alongside the per-species LES wiring at the
scenario/solver-config layer (a `Vec<T>` parallel to `scalars`, or a single
scalar applied uniformly across species when the scenario declares one `Sc_t`
for all — both are legal; the registry `Species` struct itself does NOT carry
`Sc_t` because it is a turbulence-scalar *coupling* parameter, not an
intrinsic species property, matching REQ §3's phrasing "`SGS flux closed with
Sc_t (default 0.7)`" as a closure-level default, not a per-molecule constant).

---

## 3. Data-structure mapping to `SoaFields` — code references

### 3.1 What the landed machinery already supports (verified in code, this session)

- **No W-SCAL code has landed yet.** `grep -n "D3Q7\|htmp\|conc\|scalar" crates/lbm-core/src/{fields.rs,lattice.rs,solver.rs,kernels.rs}`
  (run 2026-07-07 against `main` at commit `19c1b57`) returns **zero** hits
  for `D3Q7`/`htmp`/scalar-ADE symbols — WSCAL_PASSIVE's O1 order has not been
  dispatched/merged. This spec's §3.2 `ScalarField<T>`/`scalars: Vec<...>`
  therefore **replaces** WSCAL_PASSIVE §3.2's `h`/`htmp`/`conc` Option triplet
  in the same PR/order that implements WSCAL_PASSIVE O1 — see §7's revised
  order plan; there is no already-landed single-component code to migrate.
- **`Backend::Fields` reserves `g`/`h` distribution sets** (verbatim,
  `crates/lbm-core/src/backend.rs:130-135`) — unchanged claim, still true,
  still the structural anchor both WSCAL_PASSIVE and this spec build on.
- **`SoaFields<T>`** (`crates/lbm-core/src/fields.rs:168-210`) holds `f`,
  `ftmp`, compact-core moments, and the `Option`-gated `force_field` (`:196`),
  `omega_field` (`:199`) — the precedent `omega_s_field` (§2.2) copies exactly.
- **`WaleLes::nu_t()`** (`crates/lbm-core/src/les.rs:102-105`) returns `&[T]`
  "in global compact order" — the exact input §1.4's `D_eff,k` formula reads.
  **One-step lag documented at `les.rs:5-7`** — carried into §1.4's timing note.
- **`set_omega_field`** (`crates/lbm-core/src/solver.rs:2663-2690`) is the
  `Option`-swap precedent §2.2's `set_omega_s_field` copies verbatim, once per
  species.
- **Generic population halo** `exchange_f_generic::<L, T>`
  (`crates/lbm-core/src/halo.rs:308`) — unchanged; called once per species per
  step with `L = D3Q7` (§4.3), exactly as WSCAL_PASSIVE §3.4 specified for its
  single `h`. No new halo mechanism for `N>1` — it is the same call, N times.
- **`cx/wvof-o1`** exists as a git branch but its diff against `main` (this
  session, `git diff main cx/wvof-o1 --stat`) shows only **deletions** of
  files that exist on `main` (including `WSCAL_PASSIVE_SPEC.md` itself),
  confirming it predates both specs and is a **stale/abandoned branch**, not
  an in-flight W-VOF landing. **W-VOF is PENDING per REQ §0**, matching
  REQ_STIRRED_REACTOR.md's own landed/pending table. This spec's §1.5
  two-phase form is therefore written against an **absent** `φ` today; its
  `α≡1` degeneration path (§1.5) is what actually ships in O1/O2 below, and
  the `φ`-weighted path activates automatically once W-VOF lands (no
  W-SCAL-multi code change needed — this is the coexistence contract, §8).

### 3.2 The W-SCAL-multi additions to `SoaFields<T>`

Replaces the WSCAL_PASSIVE §3.2 triplet (`h: Option<Vec<T>>`, `htmp`, `conc`)
with:

```rust
pub scalars: Vec<ScalarField<T>>,   // see §2.2 for ScalarField<T> definition
```

`scalars.len() == 0` is the `None`-equivalent legacy path (B-6 invariance,
decision M1/M8) — no allocation beyond the empty `Vec`'s header, and every
scalar sub-step loop (`for k in 0..scalars.len()`) executes zero iterations.

### 3.3 Memory cost per cell (N species, D3Q7 h × N, matches NFR-01 additively)

Per the REQ NFR-01 budget table row "Scalar h (per component) — D3Q7 × 2 ×
f32 = 56", this spec's cost is **strictly additive per species**, no
per-species overhead beyond WSCAL_PASSIVE's own row:

| Component | Layout | B/cell (f32), per species |
|---|---|---|
| Scalar `h` (D3Q7 × 2) | 7 × 2 × f32 | 56 |
| `conc` (compact core) | 1 × f32 | 4 |
| `omega_s_field` (compact core, only when LES/phase-weighted active) | 1 × f32 | +4 (conditional) |

**N species ⇒ ≈ 60 N B/cell** (64N with per-cell `omega_s_field` active for
every species) — e.g. the envelope's 5-species set (O2, CO2, H2, N2, glucose)
is **≈ 300–320 B/cell**, well inside the REQ §7 NFR-01 "development /
validation ≤256³ (1.7e7 cells ≈ 10 GB)" envelope (5 species × 1.7e7 cells ×
64 B ≈ 5.4 GB, additive to the base hydrodynamic budget, still single-digit-GB
on the M5 Max dev box). At `N=8` (adding a couple of ionic tracers) the cost
is ≈ 8.6 GB — still ≤256³-feasible. No new NFR-01 budget row is needed; the
existing row is correct as an **N=1** instance and this spec documents its
linear scaling rather than replacing it.

### 3.4 Halo plan (unchanged mechanism, applied N times)

`exchange_f_generic::<D3Q7, T>` is called once per species per step
(§4.3) — no batching/fusing across species in phase 1 of this spec (a
per-species halo call is the WSCAL_PASSIVE-established mechanism; batching N
species' halo traffic into one call is a legitimate future performance
optimization but is NOT required for correctness and is explicitly deferred —
see §7 O2 scope note — to avoid over-scoping this spec with an unrequested
optimization, consistent with CLAUDE.md "minimal scope").

---

## 4. Pass structure — where the N-species sub-step slots in

### 4.1 Unchanged invariant order (verified `backend.rs:258-300`, `run_span`)

Identical to WSCAL_PASSIVE §4.1 — this spec touches nothing about the
hydrodynamic `f` pass order. Restated for completeness:

```
collide → exchange_f (halo) → stream (interior, then boundary shells)
        → apply_bouzidi → swap → apply_open_faces → update_moments
```

### 4.2 Where the multi-species scalar sub-step slots in (generalizes WSCAL_PASSIVE §4.2)

Unchanged **position** — after the hydrodynamic step, at solver-orchestration
level (same slot as `update_shan_chen_force` `solver.rs:2381` and WALE's
`set_omega_field` `solver.rs:2596`) — generalized to iterate over species and
to insert the SGS/phase-weighting pre-pass WSCAL_PASSIVE explicitly deferred:

```
1. HYDRODYNAMIC f STEP (unchanged run_span): collide → halo → stream →
   open BCs → update_moments. Produces rho, u, AND (if W-LES active)
   nu_t via WaleLes::update() [wherever WALE already runs it — unchanged].

2. [NEW, this spec] PER-SPECIES tau_s / D_eff REFRESH (host-side, before the
   scalar collide loop, once per step, O(N x n_cells)):
   a. If W-VOF phi is present: read alpha(x) = phi(x) (borrowed, not copied
      — see the §2.2 SoaFields.scalar_phase_weight wiring note).
   b. For each species k: D_k(x) = D_gas_k + alpha(x)*(D_liq_k - D_gas_k)
      [degenerates to the scalar D_k when no phase field, decision M4/M5].
   c. If W-LES active: D_eff_k(x) = D_k(x) + nu_t(x)/Sc_t_k  [decision M6].
   d. If D_eff_k(x) is spatially uniform (no LES, no phi): skip (b)/(c),
      keep the existing uniform tau_s fast path — no field write, no cost.
      Otherwise: tau_s_k(x) = 4*D_eff_k(x) + 0.5; solver.set_omega_s_field(k, Some(&tau_s_k)).

3. SCALAR ADE SUB-STEP (generalizes WSCAL_PASSIVE step 2, now `for k in 0..N`):
   for k in 0..scalars.len() {
       a. h collide (WSCAL_PASSIVE Eq 4/5, with alpha-weighted eq per Eq 3
          of THIS spec): read C_k = Sum h_i^(k) / alpha(x) and the resolved
          u from SoaFields::{ux,uy,uz}; relax toward h^eq with
          omega_s_k(x) [uniform or per-cell, per step 2].
       b. exchange h^(k) halo: exchange_f_generic::<D3Q7,T> on scalars[k].h.
       c. h^(k) stream: pull-stream -> scalars[k].htmp; swap.
       d. scalar BCs (WSCAL_PASSIVE §2, per-species-configurable):
          bounce-back / anti-bounce-back Dirichlet / zero-gradient, on this
          part's global faces, using species k's own BC config (§2.2
          set_species_face_bc).
       e. C_k = alpha(x) * Sum h_i^(k)  ... wait, see Eq (4): C_k = Sum h_i^(k) / alpha(x).
   }

4. (out of scope, this spec — later DAG nodes) reaction split-step R_k(C);
   interfacial S^if flux; active-feedback rho(C)/mu(C); F_b^scalar
   composition. NONE run in this spec (decision M9).
```

**Ordering rationale for step 2 (new, decided here)**: the `D_eff,k` refresh
must run strictly between the hydrodynamic `update_moments` (which produces
this step's `u` AND, when LES is active, is the point after which `WaleLes`
recomputes `nu_t` for the *next* step per its documented one-step-lag
convention, `les.rs:5-7`) and the scalar collide loop, so that every species'
`τ_{s,k}` in a given step uses the **same** `nu_t` snapshot the hydrodynamic
collide of the **next** step will also use — internal consistency inherited
from WALE's existing convention (§1.4's timing note), not a new lag.

**Per-species independence**: step 3's loop body has no data dependency
across `k` (species do not interact — decision M3/M9), so the `N` iterations
are trivially parallelizable (a future SIMD/GPU optimization, not required by
this spec — CpuScalar reference runs them sequentially, §7 O1 scope).

### 4.3 CPU-first / GPU staging (unchanged posture, generalizes WSCAL_PASSIVE §4.3)

- **Phase 1 (this spec's O1, CpuScalar reference).** The `N`-loop of §4.2 step
  3 against `SoaFields::scalars`, on `CpuScalar` (`backend.rs:125`).
- **Phase 2 (CpuSimd fused).** Same gate as WSCAL_PASSIVE: `backend_simd_equiv.rs`
  + T13 must stay green; deferred to a follow-on order, not blocking this
  spec's landing (CLAUDE.md invariant, restated).
- **Phase 3 (GPU).** Deferred behind B-1 (PARTIALLY RESOLVED), unchanged from
  WSCAL_PASSIVE's posture — `N` independent per-species buffers make the GPU
  staging problem `N`× WSCAL_PASSIVE's single-buffer problem, not a new kind
  of problem; no new blocking dependency introduced.

### 4.4 Force composition — unchanged, still nobody writes (decision M9, restates WSCAL_PASSIVE P9)

This spec composes **no force**, for any species. The reserved
`F_b^scalar` slot (FORCE_COMPOSITION_SPEC T5) is unaffected — T5(b) already
anticipates "After W-VOF / W-SCAL" as the trigger for implementing the
Boussinesq contributor; **this spec is not that trigger** (T5 remains
`(a) not implemented` after this spec lands, same as after WSCAL_PASSIVE
alone). Any of the N species' `conc` fields is available as a candidate
`s(x)` input to a future T5 closure — no different from WSCAL_PASSIVE's single
`conc`, just N of them now.

---

## 5. Validation plan mapped to T17 (extends WSCAL_PASSIVE §5 to N species + SGS + partition)

All bands are **provisional MVP**, each states its denominator/normalization
per `.claude/skills/lbmflow-physics-discipline` Rule 3 (band + behavior
anchor, both required). Tests authored adversarially by codex/Opus from this
spec, in a worktree that never shares with the implementation worktree.

| ID | Test | Metric & band (denominator stated) | Behavior anchor | Grid / steps / backend | T17 row |
|---|---|---|---|---|---|
| **V6** | **Multi-species Taylor–Aris, independent dispersion.** N=5 species (envelope §"ENVELOPE" values: O2 D=2.5e-9, CO2 D=1.9e-9, H2 D=4.5e-9, N2 D≈2.0e-9 [Cussler Table for N2 in water at 37C — record exact value used in PHYSICS.md at impl time], glucose D=6.7e-10 m²/s) co-injected in the SAME plane-Poiseuille channel as WSCAL_PASSIVE V3, same `u`. | Per-species `D_eff,k` vs Taylor–Aris analytic `D_eff,k = D_k(1+Pe_k²/210)`, `Pe_k = Ū·H/D_k` (denominator = the analytic `D_eff,k` itself, i.e. relative error). Band: **±10%** per species (same provisional band as WSCAL_PASSIVE V3, now checked N times independently). | **Each species' measured `D_eff,k` ranking must match its `D_k` ranking** (H2 disperses fastest, glucose slowest — monotonic in `D_k` at fixed `u`, since `Pe_k` ∝ `1/D_k` and the Taylor term dominates at these `D_k` values) — a cross-species monotonicity check, catching a species-index mixup that a per-species-only band could miss (e.g. species 2's buffer accidentally reading species 3's `D`). | 3D channel `H=64`, `L≥8H`, no-flux walls, Pe∈{10,50,100} per species (5 species × 3 Pe = 15 runs, or one run with 5 simultaneous species at their physical `D_k` ratios — both required, see negative-test list below), CpuScalar | VR-STR-04 |
| **V7** | **Total-species-mass conservation, N species, closed box, no reaction.** Generalizes WSCAL_PASSIVE V4 to N independently-initialized species in the same closed no-flux box. | For each species k independently: `\|Σ_cell C_k(t) − Σ_cell C_k(0)\| / Σ_cell C_k(0)` **< 1e-12** (f64) / **< 1e-6** (f32) at every step — round-off only (denominator = each species' own initial total mass; NOT normalized by a cross-species total, which would let one species' error hide in another's larger mass). | **No species' error correlates with another's initial concentration magnitude** (a shared-buffer aliasing bug would show a fixed absolute error scaling with a DIFFERENT species' mass) — checked by initializing species at deliberately mismatched magnitudes (1e-6, 1e0, 1e6 in lattice units) and confirming each relative error stays at its own round-off floor independent of the others' magnitude. | closed box `64×64×64`, N=5, mismatched initial magnitudes (above), stirred `u`, 20k steps, CpuScalar (+CpuSimd for partition invariance below) | VR-STR-05 |
| **V7b** | **Partition invariance (T13), N species.** Same V7 run under any subdomain decomposition. | Bit-identical `C_k(t)` for every k, every decomposition (same T13 mechanism WSCAL_PASSIVE V4b used, now over N independent halo exchanges). | Bit-identity is itself the strongest possible anchor — no separate behavior check needed beyond confirming the exchange loop (§3.4) correctly indexes species k's OWN halo buffer, not a shared/aliased one. | same as V7, 1/2/4-part decompositions, CpuScalar+CpuSimd | VR-STR-05 (T13) |
| **V8** | **`Sc_t` ablation guard (SGS scalar flux, the FR-LES-04 gate).** Same channel as V6, W-LES active (`ν_t>0` measured, nonzero), run once WITH the SGS term (`D_eff,k=D_k+ν_t/Sc_t`) and once with it artificially disabled (`D_eff,k≡D_k`, i.e. `Sc_t→∞` equivalent). | `\|D_eff,k(Sc_t=0.7) − D_eff,k(no-SGS)\| / D_eff,k(no-SGS)` **> band-width of V6's ±10% band** (i.e. **> 10%**, proving the term is not decorative — physics-discipline Rule 3 "ablation guard" template, applied with V6's own band as the required-detectable-effect threshold). Requires a flow/resolution combination where `ν_t` is not itself negligible (reuse a WALE-characterized turbulent channel, e.g. the landed Re_τ=178 MKM channel geometry, PHYSICS.md 2026-07-07 WALE entries) rather than the laminar V6 channel (laminar `ν_t≡0` would make this test vacuous — a **known trap**, explicitly avoided here). | The `D_eff,k` INCREASE from disabling-to-enabling SGS must be **larger for smaller molecular `D_k`** in relative terms (SGS flux is a fixed additive `ν_t/Sc_t` term — its relative contribution to `D_eff` is larger when the molecular `D_k` baseline is smaller, e.g. glucose's SGS contribution dominates its `D_eff` far more than H2's does) — a cross-species monotonicity anchor analogous to V6's. | turbulent channel (WALE-characterized geometry, `ν_t` confirmed nonzero via `WaleLesDiagnostics`), N≥2 species spanning the envelope's `D` range, CpuScalar+GPU(if available) | VR-STR-04 (FR-LES-04) |
| **V8b** | **`ν_t` leak guard, LES-inactive.** Generalizes WSCAL_PASSIVE's negative arm: with NO `WaleLes` driver installed (or `nu_t≡0`), assert `D_eff,k == D_k` exactly (bit-level, not approximate) for every species — no silent SGS contribution when LES is off. | n/a (exact-zero degeneration, not a statistical band) | This is itself the behavior anchor (WSCAL_PASSIVE's own template, §5.6d): a nonzero `D_eff,k − D_k` here is unconditionally a bug, not a tunable. | any scenario, LES off, N≥1, CpuScalar | VR-STR-04 |
| **V9** | **Per-phase partition consistency.** With a W-VOF `phi` field present (mocked/injected directly into `SoaFields` for this test if W-VOF's own scenario plumbing isn't landed yet — see the O-order dependency note in §7), verify: (a) `α≡1` degeneration recovers WSCAL_PASSIVE's own V1 pure-diffusion profile bit-for-bit (decision M4's degeneration claim, empirically checked, not just argued); (b) with `φ` a step function (sharp liquid/gas half-box, no dynamics — a static `φ` injected for this unit test only) and `D_liq ≠ D_gas` for one species, the measured local diffusive spreading rate on each side matches `D_liq`/`D_gas` respectively within the V6 band (±10%), i.e. the interpolation (§1.5) actually selects the right per-phase `D`, not an average. | (a) exact bit-identity vs the WSCAL_PASSIVE V1 baseline; (b) ±10% per-side `D_eff` vs the phase's own `D` (same denominator convention as V6). | (b)'s anchor: the spreading-rate DISCONTINUITY at the `φ` step must be visible and located AT the step (not smeared across the whole domain) — a diffuse-interface artifact (over-wide numerical smearing of the property jump) would show as the transition band exceeding the configured `φ` interface width `W`, caught by the same interface-width diagnostic WVOF already specifies. | 1D-in-3D diffusion, static half-box `φ` injection, CpuScalar | VR-STR-04, REQ §3 phase-wise form |
| **V10** | **Element-balance forward-guard (decision M7).** With `R_k≡0` (this spec's scope), assert the host-side element-balance utility (§1.6) evaluated on ANY user-declared stoichiometric map `e` reproduces exactly the per-species V7 conservation (i.e. `Σ_k e(k) × drift_k` is bounded by `Σ_k |e(k)| × (V7's own per-species round-off bound)` — a linear-combination identity, not a new physical claim). | Exact linear-algebra identity check (round-off only, same bound as V7). | This test is explicitly a **forward guard**: it is written so that if a future reaction order adds `R_k≠0` WITHOUT wiring the corresponding element map, THIS test starts failing (element balance breaks while V7's raw per-species mass check would legitimately also change, but by an amount inconsistent with any valid stoichiometry) — documented in the test's own doc-comment as "this test's PASS today is a corollary of V7, not new physics; its job is to fail loudly when reaction lands without bookkeeping." | same as V7, CpuScalar | VR-STR-05 (REQ §3 conservation statement) |
| **V11 (bit-identity DoD)** | **`N=0` and `N=1`-matches-WSCAL_PASSIVE bit-identity.** (a) Any scenario with `scalars=Vec::new()` (N=0) produces a `probe_state_hash`-equivalent regression match to the pre-W-SCAL(-multi) engine (decision M1/M8, first half). (b) Any scenario with exactly N=1 species, no `φ`, no LES, reproduces WSCAL_PASSIVE's own V1/V2/V4/V5 test outputs bit-for-bit under the generalized `Vec`-backed storage (decision M8, second half — this spec must not silently change WSCAL_PASSIVE's single-component behavior). | exact bit-identity (regression-style match, VR-STR-05 semantics, "single-backend regression only" per REQ FR-COUP-04) | n/a — bit-identity is the anchor | (a) cavity + cylinder presets; (b) re-run WSCAL_PASSIVE's own V1/V2/V4/V5 fixtures verbatim against this spec's storage, CpuScalar | VR-STR-05 (B-6 invariance) |

**Mandatory negative / consistency tests (extends WSCAL_PASSIVE §5's list):**

- **Species-index aliasing (V6/V7 negative arm):** a mutant that swaps two
  species' `D_k`/`τ_{s,k}` mid-registry (e.g. an off-by-one in the
  `SpeciesId` → `scalars` index map) must FAIL both V6's monotonicity anchor
  AND V7's mismatched-magnitude anchor — this is the concrete test that
  exercises the "no code branches on species identity, only on index" claim:
  if indices are handled correctly, WHICH species sits at index 3 is
  irrelevant to correctness; a mutant that assumes a specific index means a
  specific species (e.g. hardcoding "index 0 is always O2") must be caught.
- **Sc_t decorativeness (V8, restated as the Rule-3 requirement):** V8 IS the
  mandatory ablation guard for the SGS term — its band is written to FAIL
  (i.e. the guard is satisfied) precisely because disabling the term changes
  the metric by more than V6's band-width; a version of V8 that shows no
  detectable difference between SGS-on and SGS-off would mean `Sc_t`'s
  implementation is decorative (Rule 2 ban) and must be fixed, not the test
  loosened.
- **Phase-interpolation-direction sign (V9 negative arm):** a mutant that
  swaps `D_liq`/`D_gas` in the §1.5 interpolation formula (i.e. treats `φ=1`
  as gas instead of liquid, inconsistent with the CLAUDE.md/REQ `φ=1 liquid`
  convention) must FAIL V9(b)'s per-side band (the wrong side gets the wrong
  `D`) — this is the scalar-transport analog of REQ FR-VOF-03's sparger
  phase-convention test.
- **Harmonic-vs-linear misapplication guard (V9, documentation-level):** if a
  future implementer copies WVOF's harmonic-in-μ formula for `D_k(x)` instead
  of this spec's linear-in-D formula (§1.5), V9(b)'s measured per-side
  `D_eff` will still individually match each pure-phase limit (harmonic and
  linear interpolation agree at `φ∈{0,1}` by construction) — **V9(b) alone
  cannot catch this substitution**. A dedicated **V9c** (not yet in the table
  above — flagged here as a required addition at O3 authorship time) must
  test an INTERMEDIATE `φ` value (e.g. `φ=0.5`) against the linear-formula
  prediction specifically, where harmonic and linear diverge, closing this
  gap. (Recorded as a spec requirement on the O3 test author, not deferred
  silently.)

### 5.6 Behavior-validity review (mandatory, extends WSCAL_PASSIVE §5.6 to N species / SGS / partition)

After each validation run, before reporting, in addition to WSCAL_PASSIVE's
(a)-(d) (Gaussian symmetry, Galilean invariance, linear variance growth,
non-negativity):

(e) **Cross-species consistency**: with N species sharing the same `u` field,
each species' plume shape must differ ONLY by its diffusivity's effect on
spreading rate — any species showing a shape anomaly the others don't (e.g. a
lattice-aligned artifact appearing in ONE species but not others at similar
`τ_s`) indicates a per-species buffer/indexing bug, not a physics effect
(species do not interact, M3 — any species-specific shape difference beyond
what its own `D_k`/`τ_{s,k}` predicts is a bug).
(f) **SGS contribution spatial correlation**: where `ν_t(x)` is large (near
resolved turbulent structures, per the existing WALE diagnostics), `D_eff,k`
should visibly track it (a species' effective spreading should broaden in
high-`ν_t` regions and narrow in low-`ν_t` regions within the same run) — a
`D_eff,k` that is spatially uniform despite a highly non-uniform `nu_t()`
output signals the `omega_s_field` per-cell wiring (§1.4) did not actually
activate (fell through to the uniform fast path silently).
(g) **Phase-boundary sharpness** (when `φ` is present, V9): the diffusive
spreading-rate transition must be located AT the `φ` interface, within `W`
(WVOF's interface-width parameter) — not smeared across the domain (a
too-wide transition indicates either an interpolation bug or, if `φ` itself is
diffuse per WVOF's own model, a legitimate physical smear that must be
distinguished from a bug by checking against WVOF's own characterized `W`).

Record the review in PHYSICS.md or the track's findings file, per the
existing template (physics-discipline skill, "Post-run behavior-validity
review").

---

## 6. Stability & parameter domain

### 6.1 Relaxation window — unchanged per species, independently

`τ_{s,k} = D_{eff,k}/cs_s² + ½ = 4 D_{eff,k} + ½` (WSCAL_PASSIVE Eq. 6,
unchanged formula, applied per `k`). Operating band `τ_{s,k} ∈ (0.5, ~1.0]`
for BGK, wider under TRT `Λ=1/4` — identical to WSCAL_PASSIVE §6.1, checked
**independently for every species** (a fast-diffusing species like H2,
`D≈4.5e-9` m²/s, and a slow one like glucose, `D≈6.7e-10` m²/s, will generally
require DIFFERENT lattice `D` values at a shared `Δx`/`Δt`, hence different
`τ_{s,k}` — the unit-conversion/feasibility layer must check EACH species'
`τ_{s,k}` against the window, not just one representative value. This is a
config-validation requirement, not a new physics decision — restated from
WSCAL_PASSIVE §6.1, generalized to N independent checks.)

### 6.2 Grid-Péclet — per species, per phase

`Pe_{Δ,k}(x) = |u| Δx / D_{eff,k}(x)` — WSCAL_PASSIVE §6.2's bound (`≲2` BGK,
wider TRT) applies per species, evaluated at the SMALLEST `D_{eff,k}(x)`
across both phases and across the domain for that species (the slowest-
diffusing phase/region is the binding constraint) — a resolution requirement
per REQ §2, never a clamp (restated, unchanged principle).

### 6.3 Mach / low-Mach consistency — unchanged

Identical to WSCAL_PASSIVE §6.3; the linear equilibrium (this spec's Eq. 1) is
consistent to O(Ma²) exactly as WSCAL_PASSIVE's Eq. 3 was — no new Mach
constraint from adding `N` or the `α`/`ν_t` weighting (both are cell-local
multiplicative factors on the SAME linear-in-`u` equilibrium form).

### 6.4 `Sc_t` domain — restates §1.4's validity domain

`Sc_t ∈ [0.2, 1.3]` is the literature-measured range (Tominaga & Stathopoulos
2007); this spec's default 0.7 (REQ-mandated) is a point within that range,
not the full domain — a scenario that sets `Sc_t` outside `[0.2,1.3]` should
warn (config-validation-time), matching the discipline WSCAL_PASSIVE §6.2
applies to grid-Péclet (a resolution/parameter-domain warning, not a clamp on
the transported quantity itself — `Sc_t` is a scenario INPUT parameter, not a
transported field, so bounding it is not the banned "transport-absorbing
clamp" pattern; it is ordinary input validation, same class as the existing
`Ma_lattice≤0.1` check).

### 6.5 Positivity — unchanged, restated

Identical discipline to WSCAL_PASSIVE §6.5, per species independently: `C_k≥0`
is a diagnostic, never a clamp, for every `k`.

### 6.6 Explicit non-goal: no cross-species diffusion coupling

Restates decision M3 as a stability-domain boundary: this spec's validity
domain EXCLUDES scenarios requiring Stefan-Maxwell / Fickian cross-diffusion
coupling between species (relevant at high concentration or strongly
non-ideal mixtures — Cussler 2009 ch. 3). The dilute-solution approximation
(each species obeys its own independent Fickian ADE, REQ §3's form) is assumed
throughout. This is a stated model limitation (PHYSICS.md entry, §7.4), not a
silent gap — flagged here so a future high-concentration reactor scenario does
not silently misuse this spec outside its domain.

---

## 7. Phased landing plan — CODEX ORDER BREAKDOWN

Since **no WSCAL_PASSIVE code has landed** (§3.1), this spec's orders
**supersede and absorb** WSCAL_PASSIVE §7's O1/O2/O3 rather than layering on
top of already-merged code. The PM should dispatch the orders below INSTEAD
OF WSCAL_PASSIVE's original three; WSCAL_PASSIVE remains the correct
per-species derivation reference (§1.1 cites it verbatim) but its own §7 order
plan is superseded by this table.

| Order | Scope | Primary files (conflict boundary) | Gate |
|---|---|---|---|
| **O1 — D3Q7 ADE core + N-species storage (CpuScalar)** | D3Q7 `Lattice` impl (WSCAL_PASSIVE §1.4, unchanged); `ScalarField<T>` struct + `SoaFields.scalars: Vec<ScalarField<T>>` (§2.2, §3.2); per-species BGK+TRT collide (§1.3 Eq 1/2) with `α≡1` fast path (no `φ` yet — W-VOF absent, §3.1); `h^{(k)}` halo via `exchange_f_generic::<D3Q7>` called per species (§3.4); `h^{(k)}` stream+swap; `C_k=Σh_i^{(k)}` (`α=1` case); scalar sub-step wiring at solver level iterating `0..N` (§4.2 steps 1+3, step 2's `D_eff` refresh WITHOUT the LES/phase terms yet — i.e. `D_eff,k≡D_k` uniform, deferred to O2); scalar BCs (WSCAL_PASSIVE §2) generalized to per-species-configurable (§2.2 `set_species_face_bc`); species registry `Species`/`SpeciesId` (§2.1, scenario-adjacent but the core-facing numeric-only API lives here). | `lattice.rs` (add D3Q7), `fields.rs` (add `ScalarField`/`scalars` Vec), `solver.rs` (scalar sub-step orchestration, `add_species`/`set_species_tau`/`set_species_face_bc`), `kernels.rs` (h collide/stream/BC rows, generic over species index) | V7, V7b, V11(a,b) green on CpuScalar. (V6/V8/V9 need O2's `D_eff` refresh — not gated here.) |
| **O2 — SGS scalar flux (`Sc_t`) + per-phase `D` interpolation** | The §4.2 step-2 `D_eff,k(x)` refresh (reads `WaleLes::nu_t()`, writes `omega_s_field` per species via `set_omega_s_field`, §1.4/§2.2); per-species `Sc_t` scenario parameter (§2.3); the §1.5 `α(x)`-weighted equilibrium (Eq 3/4) AND the linear-in-D per-phase interpolation (§1.5), wired to read `φ` from the W-VOF `phi` slot **when present**, `α≡1` fast path **when absent** (the `SoaFields.scalar_phase_weight` borrow-vs-copy wiring, §2.2's placeholder note — resolve the concrete borrow mechanism here, against whatever W-VOF's `phi` storage looks like at the time this order runs; if W-VOF has NOT landed yet, implement ONLY the `α≡1` path and leave the `φ`-read as a `todo!()`-free `None`-only stub gated by a scenario config flag that the schema validator rejects as "unsupported: two-phase scalar requires W-VOF" — do NOT block O2 on W-VOF landing). | `solver.rs` (D_eff refresh pre-pass, distinct function from O1's sub-step), `les.rs` (read-only use of `nu_t()`, no changes expected), `fields.rs` (`omega_s_field` per `ScalarField`, `scalar_phase_weight` wiring) — **mostly disjoint from O1's kernel-row work, same files but different functions; sequence after O1, do not parallelize against O1 in the same files** | V6, V8, V8b green (single-phase, LES on/off). V9 green ONLY if W-VOF's `phi` is available by then (mocked/injected `φ` for the unit test is acceptable per §5 V9's own note — this does not require W-VOF's full scenario plumbing, only a `phi: Vec<T>` buffer the test constructs directly against `SoaFields`). |
| **O3 — Scenario/CLI plumbing: registry + multi-species schema** | `Species`/registry schema (name, Henry, D_liq, D_gas, charge, MW, volatile — §2.1) in scenario JSON; config validation (per-species grid-Péclet warning generalizing WSCAL_PASSIVE O2's single-species check, §6.1/§6.2, N independent checks); `Sc_t` scenario field (§2.3); CLI/output of each species' `conc` field (VTI, `manifest.json`, one output channel per species, labeled by `Species.name` — the ONE place `name` is used, for human-readable output labeling, never for kernel branching per M2); the multi-species Taylor–Aris channel scenario (V6's fixture). | `crates/lbm-scenario/src/lib.rs` (schema + validation), `crates/lbm-cli` (per-species `conc` output) — **disjoint from O1/O2's core files** | V6's scenario-level fixture constructible and running end-to-end via CLI. Depends: O1 (needs `add_species` API); parallelizable against O2 (scenario schema doesn't need `D_eff` internals, only the `Species` record shape). |
| **O4 — Validation authorship (codex adversarial, separate worktree)** | All of §5 (V6–V11) + the negative/consistency tests (species-index aliasing, Sc_t decorativeness / the V8 ablation itself, phase-interpolation sign, harmonic-vs-linear V9c per the flagged gap) + the §5.6 behavior-validity review record for each. Authored from THIS spec, not from the impl. | `crates/lbm-core/tests/wscal_multi_*.rs` + `crates/lbm-scenario/tests/*` (new files only — no impl-file conflict) | Tests compile red against a stub, go green against O1/O2/O3 as they land; freeze bands in VALIDATION.md T17 VR-STR-04/05. Runs alongside O1-O3. |

**Critical-path ordering:** O1 → O2 (same files, sequenced). O3 depends only
on O1's `add_species`/`ScalarField` API (parallelizable against O2). O4 runs
concurrently from the start (test worktree, per CLAUDE.md team convention:
"a test order and an implementation order never share a worktree").

### 7.1 In-flight conflict surface with `cx/wvof-o1`

**Resolved finding (this session)**: `cx/wvof-o1` is a **stale branch** — its
diff against `main` shows only deletions of files (including
`WSCAL_PASSIVE_SPEC.md`) that exist on current `main`, meaning the branch was
forked BEFORE those files were added and has not been rebased since; it does
NOT represent in-flight W-VOF O1 work competing for `fields.rs`/`solver.rs`/
`kernels.rs` edits right now. **W-VOF is PENDING** per REQ §0's own
landed/pending table — there is no live conflict to resolve today. If/when a
FRESH W-VOF O1 order is dispatched (superseding the stale branch), the same
conflict-boundary analysis WSCAL_PASSIVE §7.1 already worked out applies
unchanged (both additions are `Option`/`Vec`-typed additive struct fields;
mechanical rebase, no semantic merge) — this spec's `scalars: Vec<ScalarField>`
addition is exactly as additive as WSCAL_PASSIVE's original `h`/`htmp`/`conc`
triplet was, so WSCAL_PASSIVE §7.1's merge-order rule carries over verbatim
with `scalars` substituted for `h`/`htmp`/`conc`.

### 7.2 Per-order DoD (all orders)

Existing tests green *without modification*; `scalars=Vec::new()` path
bit-identical to today (V11a); the phase-2 PHYSICS.md entry (§7.4) landed with
O1/O2 (split if convenient — one entry for the D3Q7-multi-species mechanics
at O1, one for the SGS/`Sc_t` closure at O2, matching Rule 1's "one entry per
new term" granularity); behavior-validity review (§5.6) recorded for every
validation run; `backend_simd_equiv.rs` + T13 green (exercised at `N=0`,
unaffected by this spec's `CpuScalar`-only phase 1).

### 7.3 STOP-RULE check (explicit, per the task instruction)

No gate in §5 requires a banned pattern to reach its band:
- V6/V7/V8/V9's bands are all either round-off-floor (V7/V11, not a tunable
  band at all) or the SAME ±10% Taylor–Aris band WSCAL_PASSIVE already
  established as reachable (its own V3 gate, unchanged physics, just N
  independent instances) — no new calibration risk.
- V8's ablation guard is satisfied BY CONSTRUCTION once `Sc_t` is wired
  correctly (the guard measures whether the term does something, and the term
  is a straightforward additive `ν_t/Sc_t`, not a fitted closure that might
  fail to produce a detectable effect) — no STOP-RULE risk identified.
- The one flagged gap (§5's "Harmonic-vs-linear misapplication guard," V9c)
  is a **test-completeness** note for the O4 author, not a physics
  blocker — the correct closure (§1.5) is already decided and derived; V9c
  just needs to be written to actually distinguish it from the alternative,
  which is a test-design task, not a STOP-RULE physics gap.

**No STOP-RULE is raised by this spec.**

### 7.4 The PHYSICS.md validity-domain statement (mandatory entry text)

O1/O2 must add, on landing, PHYSICS.md entries following the Rule 1 template.
Two entries (mechanics + closure), per §7.2's granularity note:

> **Multi-species passive scalar transport — N-species ADE-LBM generalization
> (crates/lbm-core/src/fields.rs:ScalarField, solver.rs: per-species
> sub-step).** N independent D3Q7 `h^{(k)}` distributions, one per registered
> species, each obeying WSCAL_PASSIVE's ADE-LBM mechanics
> (Krüger et al. 2017 §8.3) unchanged — linear-in-`u` equilibrium, BGK/TRT
> collision, `τ_{s,k}=D_{eff,k}/cs_s²+½`. Species do not interact in the
> collision operator (dilute-solution approximation, no Stefan-Maxwell
> coupling — §6.6). **Validity domain**: identical per-species
> `τ_{s,k}∈(0.5,~1.0]`/cell-Péclet bound as WSCAL_PASSIVE, checked
> independently for every species (the slowest-diffusing species/phase is the
> binding resolution constraint, §6.2). **Why here**: `Vec<ScalarField>`
> (not `[ScalarField; N_MAX]` or per-species-hardcoded fields) chosen so `N`
> is a runtime scenario parameter, matching REQ's "component-agnostic
> registry" requirement — no species name ever appears in a kernel branch
> (grep-verified at landing time).
>
> **SGS turbulent scalar flux — `D_eff,k = D_k + ν_t/Sc_t,k` (crates/lbm-core/
> src/solver.rs: `set_omega_s_field` pre-pass; les.rs: `WaleLes::nu_t()`
> read-only consumer).** Gradient-diffusion closure for the unresolved
> turbulent scalar flux, the scalar analog of the already-landed WALE
> eddy-viscosity closure. `Sc_t,k` defaults to 0.7 (REQ §3-mandated;
> Tominaga & Stathopoulos 2007 measured range [0.2,1.3]). **Validity domain**:
> wherever the underlying WALE `ν_t` is valid (inherits WALE's documented
> resolution caveats); constant-in-space `Sc_t` (no dynamic model) is a stated
> limitation. `ν_t=0`/no-LES ⇒ `D_eff,k≡D_k` exactly (V8b). **Why here**: the
> one-step lag on `nu_t()` (les.rs:5-7) is inherited, not introduced, by
> reading the SAME snapshot the hydrodynamic collide uses.
>
> **Two-phase per-species diffusivity — linear-in-D interpolation (contrast
> WVOF's harmonic-in-μ).** `D_k(x) = D_{k,gas} + φ(x)(D_{k,liq}−D_{k,gas})`,
> linear (not harmonic) because Fickian flux inside a resolved diffuse
> interface has no series-modulus structure (Kim 2012). **Validity domain**:
> requires W-VOF's interface band width `W` to resolve the transition (V9's
> sharpness anchor); inert (`α≡1`) when W-VOF is absent.

---

## 8. Coexistence with W-VOF (`g`) and W-SCAL-passive (single `h`)

- **W-SCAL-passive is subsumed, not parallel-mounted.** Because no
  WSCAL_PASSIVE code has landed (§3.1), there is no separate "single-`h`"
  code path to keep alive alongside this spec's `Vec<ScalarField>` — O1 above
  IS the WSCAL_PASSIVE implementation, generalized from the start. `N=1` is
  simply the smallest legal `scalars.len()`, verified bit-identical to what
  WSCAL_PASSIVE itself would have produced (V11b) — there is no dual-
  maintenance burden.
- **W-VOF (`g`, D3Q19, pre-`f` pre-pass) mounts independently**, exactly per
  WSCAL_PASSIVE §8's original analysis: `f`/`g`/`scalars` are disjoint storage
  (a `Vec` instead of an `Option` triplet does not change the additivity
  argument — `Vec::new()` and `None` are both zero-cost, zero-allocation
  legacy states). Step slots remain disjoint (phase-field pre-pass ≠
  hydrodynamic `f` step ≠ scalar sub-step, §4.2 unchanged from WSCAL_PASSIVE
  §4.2's rationale). Halo remains disjoint (`g`'s scalar-plane exchange vs
  `scalars[k]`'s `exchange_f_generic::<D3Q7>`, called N times).
- **The `φ`-read wiring (§2.2's `scalar_phase_weight`, §7 O2's resolution
  note) is the ONE genuinely new coupling point** this spec introduces beyond
  what WSCAL_PASSIVE already established: O2 must read W-VOF's `phi` buffer
  (whatever its landed storage shape turns out to be) without copying it or
  creating a second source of truth. This is flagged explicitly as an
  O2-time integration decision (not resolved here, because W-VOF's own
  storage is PENDING/unlanded — resolving it now would be speculative) but
  is NOT a blocking dependency: O2's `α≡1` path ships independently, and the
  `φ`-aware path activates additively whenever W-VOF lands, per the
  degeneration argument (§1.5) that guarantees no regression either way.
- **Force**: unchanged from WSCAL_PASSIVE §8 — this spec still writes NO
  force (M9), the `F_b^scalar` slot remains reserved for a later, separate
  active-feedback phase, now with N candidate `conc_k` inputs instead of one.

---

## 9. Load-bearing code reference index

| Claim | File:line |
|---|---|
| No W-SCAL code landed yet (grep, this session, commit `19c1b57`) | `crates/lbm-core/src/{fields.rs,lattice.rs,solver.rs,kernels.rs}` (absence of `D3Q7`/`htmp`/scalar-ADE symbols) |
| `Backend::Fields` reserves `g`/`h` distribution sets (verbatim) | `crates/lbm-core/src/backend.rs:130-135` |
| Invariant step order in `run_span` | `crates/lbm-core/src/backend.rs:258-300` |
| `Lattice` trait shape (reused unchanged for D3Q7) | `crates/lbm-core/src/lattice.rs:117-152` |
| existing D3Q19/D3Q27 impls (no D3Q7 today — new impl needed, per WSCAL_PASSIVE §1.4) | `crates/lbm-core/src/lattice.rs:259-465` |
| `SoaFields` struct, `force_field`/`omega_field` `Option` precedent for `omega_s_field` | `crates/lbm-core/src/fields.rs:168-210` (struct), `:196` (`force_field`), `:199` (`omega_field`) |
| `set_omega_field` (per-cell relaxation; `set_omega_s_field` precedent) | `crates/lbm-core/src/solver.rs:2663-2690` |
| generic population halo (parameterize `::<D3Q7,T>`, called N times) | `crates/lbm-core/src/halo.rs:308` |
| scalar-plane halo (NOT used by `h^{(k)}`, used by Shan–Chen `ψ` / candidate `φ` carrier) | `crates/lbm-core/src/halo.rs:71`, `:371` |
| `WaleLes::nu_t()` — "global compact order," the `D_eff,k` input | `crates/lbm-core/src/les.rs:102-105` |
| WALE one-step lag on `nu_t`/`omega_plus` (inherited timing, not new) | `crates/lbm-core/src/les.rs:5-7` |
| Shan–Chen force pre-pass (solver-level sub-step precedent) | `crates/lbm-core/src/solver.rs:2381`, `:2399-2499` |
| `cx/wvof-o1` is a stale/abandoned branch (diff-only-deletions vs `main`, verified this session) | git branch `cx/wvof-o1` vs `main` at `19c1b57` |
| T17 VR-STR-04 scalar/Taylor–Aris row (now multi-species); VR-STR-05 conservation | `docs/VALIDATION.md:348-349` |
| REQ three-tier scalar governing forms | `docs/REQ_STIRRED_REACTOR.md` §3 |
| REQ FR-LES-04 (`Sc_t` hook) | `docs/REQ_STIRRED_REACTOR.md` §4.2 |
| REQ dimensionless conventions (`Sc`, `Pe_N`, `Pe_tip`) | `docs/REQ_STIRRED_REACTOR.md` §2 |
| REQ NFR-01 memory budget (D3Q7 56 B/cell row, additive per species) | `docs/REQ_STIRRED_REACTOR.md` §7 |
| REQ §11 DAG: W-SCAL passive, parallel wave 1; W-REACT depends on W-SCAL (active: W-VOF) | `docs/REQ_STIRRED_REACTOR.md` §11 |
| WSCAL_PASSIVE §1–§6 (per-species ADE-LBM mechanics, cited not re-derived) | `docs/proposals/WSCAL_PASSIVE_SPEC.md` |
| WVOF harmonic-in-μ mixing rule (contrasted, not reused, for `D_k`) | `docs/proposals/WVOF_IMPL_SPEC.md` §1.3 |
| WVOF `phi`/`g` slot, `SoaFields` additions (the future `φ`-read coupling point) | `docs/proposals/WVOF_IMPL_SPEC.md` §3.2 |
| FORCE_COMPOSITION_SPEC T5 `F_b^scalar` reserved slot (unaffected by this spec) | `docs/proposals/FORCE_COMPOSITION_SPEC.md` T5 |
| capability-gap note on MCMP per-component sources (orthogonal V&V finding, not this spec's mechanism) | `docs/proposals/capability_gap_mcmp_sources.md` |

**Literature (decided references, extends WSCAL_PASSIVE's list):**
Krüger, Kusumaatmaja, Kuzmin, Shardt, Silva & Viggen 2017, *The Lattice
Boltzmann Method* §8.3 — unchanged per-species ADE-LBM base (cited, not
re-derived, in this document).
Pope 2000, *Turbulent Flows* §10.4 — gradient-diffusion SGS scalar flux
closure (the `ν_t/Sc_t` form).
Launder & Spalding 1974 — Reynolds-analogy `Sc_t` origin.
Tominaga & Stathopoulos 2007 (Atmos. Environ. 41:8091) — `Sc_t` measured
range [0.2,1.3], CFD default 0.7–0.9 — the cited validity domain for M6/§1.4.
Cussler 2009, *Diffusion: Mass Transfer in Fluid Systems* 3rd ed., Table
3.1-1 — gas/liquid molecular diffusivity ratios (motivates M5's per-phase `D`
requirement) and ch. 3 (Stefan-Maxwell multicomponent coupling, the explicit
non-goal of §6.6).
Kim 2012, *Phase-Field Models for Multi-Component Fluid Flows* (World
Scientific) — diffuse-interface scalar-diffusivity mixing-rule review
(motivates the linear-in-D choice, §1.5, contrasted against WVOF's harmonic-
in-μ).
Nicoud & Ducros 1999 — WALE (unchanged, cited via `les.rs`'s own header, the
source of the `ν_t` this spec consumes read-only).

**Envelope reference values (validation-relevant physical values at 37°C
water, task-specified; exact literature source for each to be pinned in the
PHYSICS.md entry at O1/O4 landing time, not fixed here to avoid a stale
citation if the implementer finds a more precise source at scenario-authoring
time):** O2 `D≈2.5e-9` m²/s, `C*≈0.21` mM; CO2 `D≈1.9e-9` m²/s; H2
`D≈4.5e-9` m²/s; N2 (value to be pinned, same order as O2/CO2 per Cussler);
glucose (non-volatile) `D≈6.7e-10` m²/s. Grid `≤256³` (§3.3's memory-budget
envelope).

