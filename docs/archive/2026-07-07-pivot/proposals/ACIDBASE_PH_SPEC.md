# T5b Implementation Specification — Acid–Base Equilibrium & pH as a Derived Field

**Document ID**: SPEC-T5b-ACIDBASE (rev.1, 2026-07-07).
**Scope**: the tier-5b extension of `docs/REQ_STIRRED_REACTOR.md` §3
"Scalars/reactions" — fast-equilibrium acid–base speciation and pH as a
**derived diagnostic field**, layered on top of the W-SCAL transported total-
concentration species. This is a DAE (differential-algebraic) closure, NOT a
new transport equation: it consumes the totals `C_k` that W-SCAL already
advects/diffuses and solves a local algebraic system each step to recover the
speciation and `[H+]`.
**Target core**: `crates/lbm-core`, as a per-cell post-transport pass reading
`SoaFields::conc`-family totals and writing a new derived-field slot (pH,
speciation fractions) — no new distribution, no new lattice.
**Acceptance**: VALIDATION.md **T17** row **VR-STR-04** extension (a new sub-
row, "VR-STR-04b acid–base equilibrium"), analogous in status to VR-STR-06+
(an "extension of an existing VR item," not a new VR item number) — see §5.
**Depends on**: `docs/proposals/WSCAL_PASSIVE_SPEC.md` (single-component
passive ADE, **landed spec**, phase-1 lands exactly one transported total;
this spec's multi-species need is flagged as a dependency on the phase-2
`Vec<ScalarField>` generalization reserved by WSCAL decision P10 — see §0
D0/D9)); `docs/proposals/WSCAL_MULTICOMPONENT_SPEC.md` (**not yet written**
at the time of this spec — see §0 D9 for the interim path); a `WREACT_IMPL_
SPEC.md` for the kinetic-reaction tier (**not yet written**; §8 states the
reconciliation contract this spec commits to once it exists).

This spec is **executable**: every closure choice is decided and justified,
every code touchpoint is cited against the current worktree, and every gate
is mapped to a VALIDATION.md row with a provisional band. A follow-on codex
implementation order should not need to re-derive any design decision here.

> **Discipline note (CLAUDE.md prime directive; PHYSICS.md;
> `.claude/skills/lbmflow-physics-discipline`).** The proton-condition solve
> below is resolved from mass-action equilibrium constants (`pKa`) and exact
> charge/mass balance — it is algebra, not a fitted closure. The only
> configuration inputs are literature `pKa` values and charge numbers per
> species (§2), each with its own citation. No band-calibrated constant, no
> case-keyed branch (species names never appear in solver code, only in
> configuration data — §2.1), and no result-absorbing clamp appears anywhere
> in this design (the Newton iterate is bounded by a *provable* bracketing
> interval from the physics, not an arbitrary clamp — see §3.3). Ideal-
> solution activity is the phase-1 default with a stated validity domain
> (§6); Davies-law activity correction is an explicit, separately validated
> opt-in, never silently blended in.

---

## 0. Summary of decisions (read this first)

| # | Decision | Justification (short) |
|---|---|
| D0 | **Fast-equilibrium acid–base speciation is an algebraic DAE closure applied AFTER the W-SCAL transport sub-step**, operating on the just-transported total concentrations `C_T,k` (e.g. total carbonate, total ammonia). It is NOT a new LBM distribution and NOT a new transport equation. | REQ §3 treats `R_k(C)` (reaction) as a source term added to the transported `C_k` ADE; a *fast* equilibrium is the `Da → ∞` (infinitely fast reaction) limit of that same term, which degenerates the ODE reaction into an algebraic constraint (§8, §1.5) — the standard operator-split DAE treatment (Seeburger & Fahrenholz-style geochemical speciation codes; Steefel & MacQuarrie 1996 review). |
| D1 | **The equilibrium system is species-agnostic**: a configuration object supplies (a) a list of *total* (conservative, transported) components, (b) a list of *acid/base equilibria* each relating a subset of species by a `pKa` and a stoichiometric proton count, (c) a charge number `z_i` per species, (d) `[H+]`/`[OH-]` via `Kw`. Carbonate (`pKa1≈6.35`, `pKa2≈10.33`) is loaded as a **data instance** of this generic system, not hardcoded. | CLAUDE.md ban list: "case-identity branch… must not know which test case it is running." The solver function takes `EquilibriumSystem` data, never a species name. |
| D2 | **The unknown solved for is `h ≡ [H+]` (or equivalently `pH`), one scalar Newton unknown per cell**, via the **proton condition** (charge-balance polynomial in `h`, §1.3) built from the current cell's total concentrations. All other species concentrations and `[OH-]` are then evaluated in closed form from `h` (no additional unknowns) — this is the standard single-unknown reformulation of pH-speciation problems (Morel & Hering 1993, *Principles and Applications of Aquatic Chemistry*, ch. 3; Stumm & Morgan 1996 ch. 3). | Reduces an N-species nonlinear system to ONE Newton unknown regardless of how many acid–base families are configured — the system stays agnostic (D1) and the solver stays cheap (one Newton loop per cell per step, not an N×N system). |
| D3 | **Newton's method on `log h` (not `h`) is the iteration variable**, i.e. Newton on `pH = -log10(h)` internally via `x = ln h`. | `h` spans many orders of magnitude (`10^0` to `10^-14`) inside a single titration; Newton in `h` directly is ill-conditioned near the true root when `h` is small. Newton in `x = ln h` (equivalently pH-space) is the standard well-conditioned reformulation (Van Der Vorst-style geochemical solvers; also the reformulation used by PHREEQC's Newton–Raphson core, Parkhurst & Appelo 2013 ch. 2). |
| D4 | **Operator-split placement: equilibrium closure runs after W-SCAL transport, and — when W-REACT kinetic reactions are also active — after the kinetic reaction sub-step**, each solver step. It is the *last* algebraic correction applied to the totals before the step ends. | The DAE structure is `d(totals)/dt = transport + kinetics`, `algebraic constraint = equilibrium`; standard sequential operator-splitting for index-1 DAEs (Hairer & Wanner 1996 ch. VII) applies the differential part first, then projects onto the algebraic manifold. §8 states the full reconciliation with kinetic W-REACT. |
| D5 | **The closure does not change the transported totals `C_T,k` themselves** — `C_T,k` (e.g. total carbonate) is invariant under the equilibrium solve by construction (mass action redistributes protons among species of a fixed total; it does not create or destroy the total). Only the **derived speciation fields** (per-species concentration fractions, `[H+]`, pH) are written. | This is what makes the closure "fast-equilibrium, not reaction": total-carbonate mass conservation (already gated by W-SCAL V4/VR-STR-05) is untouched; the ABLATION test (§5, A1) checks exactly this — totals identical with equilibrium on/off. |
| D6 | **Ideal-solution activity is the phase-1 default** (`a_i = [i]`, activity coefficient `γ_i = 1`); the Davies equation ionic-strength correction (`log γ_i = -A z_i² (√I/(1+√I) - 0.3 I)`) is a **separate, opt-in, independently validated** extension, never silently active. | REQ has no existing ionic-strength machinery; inventing one implicitly would be an undocumented closure. Ideal behavior is exact in the dilute limit (I → 0) and its validity domain is stated (§6.4); Davies (1962) is the standard bridge to moderate ionic strength (I ≲ 0.5 M) and is spec'd as an explicit follow-on (§6.5), never merged into the phase-1 default silently. |
| D7 | **The Newton tolerance and the DAE split-error convergence test are the two required numerical-accuracy gates** — a Newton absolute tolerance on the charge-balance residual (`|g(h)| < tol_residual`, §3.2) and a Δt-halving split-error convergence order test (§5, V3) analogous to WSCAL's grid-refinement order test. | REQ §8 mandatory-test menu includes "active-scalar dt-halving convergence (MJ-007)" — this spec is the acid–base-specific instance of that same requirement, extended to the algebraic-constraint projection. |
| D8 | **pH is purely a derived/output field** (`pub ph: Option<Vec<T>>` alongside speciation fractions), never fed back into transport, viscosity, or density in phase 1 (mirrors WSCAL's "passive" framing — no active feedback yet). | Matches WSCAL P9/REQ's phase separation between passive scalars and the later active-feedback tier; a pH-dependent rate `k(pH)` is explicitly a W-REACT-phase concern (§8), not this spec's. |
| D9 | **Interim multi-species path**: because `docs/proposals/WSCAL_MULTICOMPONENT_SPEC.md` does not exist yet, this spec designs the equilibrium closure against a `Vec<TotalSpecies>` abstraction that is satisfied EITHER by (a) `N` independent single-component W-SCAL `h`/`conc` slots run in parallel (mechanically legal today — WSCAL P10 phase 1 lands one component but nothing prevents instantiating the D3Q7 sub-step `N` times with `N` separate `Option<Vec<T>>` slots, at `N×60 B/cell`), or (b) the future `Vec<ScalarField>` generalization once multicomponent lands. The equilibrium closure's Rust API (§2.2) takes a slice of concentration values and is blind to which storage strategy produced them — **no code in this spec's O1/O2 orders depends on the multicomponent spec existing**; only the scenario-authoring convenience of "one named list of totals" benefits from it later. | Minimal scope: do not block T5b on an unwritten sibling spec. The carbonate reference case needs exactly ONE transported total (`C_T` = total dissolved inorganic carbon) plus the universal `[H+]`/`[OH-]` pair, which needs no additional transported species at all (§1.4) — so the MVP validation (§5) runs entirely on the ALREADY-LANDED WSCAL single-component path, with zero dependency on multicomponent. Systems needing ≥2 independent transported totals (e.g. coupled carbonate + ammonia) DO need (a) or (b) above; flagged as a scale-out note in §2.3, not blocking. |
| D10 | **B-6 bit-identity discipline**: the equilibrium closure is `Option`-gated (`EquilibriumSystem: Option<...>` on the solver/scenario side); when absent, zero code executes in the per-step pass and all new fields (`ph`, speciation fractions) are `None`/unallocated, bit-identical to pre-T5b. | Same discipline as WSCAL P3 / WVOF D7 — every prior landed test stays green with the feature off. |

---

## 1. Governing system: charge balance + dissociation equilibria + the proton condition

### 1.1 Mass-action equilibria (per configured acid/base family)

For a weak polyprotic acid `H_nA` with `n` dissociation steps, each step is an
independent mass-action equilibrium:

```
H_nA          <=> H_{n-1}A^-   + H+     K_a1 = [H_{n-1}A^-][H+] / [H_nA]
H_{n-1}A^-    <=> H_{n-2}A^2-  + H+     K_a2 = [H_{n-2}A^2-][H+] / [H_{n-1}A^-]
...
```

`pKa_j = -log10(K_aj)`. For carbonate (`n=2`, the reference instance, D1):

```
H2CO3*  <=> HCO3-  + H+      pKa1 = 6.35   (Plummer & Busenberg 1982, 25 C, I=0)
HCO3-   <=> CO3^2- + H+      pKa2 = 10.33  (Plummer & Busenberg 1982, 25 C, I=0)
```

(`H2CO3*` is the conventional lumped "CO2(aq) + true H2CO3" species used in
all standard carbonate-system treatments — Zeebe & Wolf-Gladrow 2001, *CO2 in
Seawater* §1.1 — a notational convention, not a modeling approximation.)

Water autoionization is always active (universal, not configured per system):

```
H2O <=> H+ + OH-      Kw = [H+][OH-],   pKw = 14.00 at 25 C (CODATA/Harned & Owen 1958)
```

### 1.2 Speciation as a function of `h = [H+]` (closed form, no extra unknowns)

For an `n`-step acid with total concentration `C_T` (the W-SCAL-transported
conservative sum of all protonation states), each species fraction is a
rational function of `h` alone (standard alpha-fraction decomposition, Morel
& Hering ch. 3; Stumm & Morgan ch. 3). Define, for carbonate (`n=2`):

```
D(h) = h^2 + K_a1 h + K_a1 K_a2
alpha_0(h) = h^2        / D(h)     [H2CO3*]  = alpha_0(h) C_T
alpha_1(h) = K_a1 h      / D(h)     [HCO3-]   = alpha_1(h) C_T
alpha_2(h) = K_a1 K_a2   / D(h)     [CO3^2-]  = alpha_2(h) C_T
```

with `alpha_0 + alpha_1 + alpha_2 = 1` identically (a closed-form algebraic
identity, asserted as a unit test, §5 A3). This generalizes to `n` steps by
the standard recursion `alpha_j(h) = (Prod_{i<=j} K_ai) h^{n-j} / D(h)`,
`D(h) = Sum_{j=0}^{n} (Prod_{i<=j} K_ai) h^{n-j}` (`K_a0 ≡ 1`). This is pure
algebra given `{K_aj}` and `h` — no iteration needed once `h` is known.

### 1.3 The proton condition (charge balance) — the equation actually solved

The unknown `h` is fixed by requiring **exact electroneutrality** across all
species present in the cell (D2): background inert ions (charge `z_i`,
concentration `C_i`, from the species registry, REQ §3/§8) plus the
acid–base family's charged species plus `H+`/`OH-`:

```
g(h) = h - Kw/h + Sum_{acid families} Sum_{j=1}^{n} (-j) * alpha_j(h) * C_T
         + Sum_{inert ions i} z_i * C_i  [+ strong-acid/base titrant charge]
       = 0                                                              (1)
```

`(-j)` is the charge of the `j`-th deprotonated species of an acid family
(e.g. `HCO3-` carries charge `-1` at `j=1`, `CO3^2-` carries `-2` at `j=2` —
consistent with the general charge-vector convention `z_i` the species
registry already carries for transported species, REQ §3 registry). Equation
(1) is the **proton condition** / TOTH (total hydrogen ion) equation in the
sense of Morel & Hering ch. 3 — it is exactly charge balance, reorganized so
that the acid-family terms are written via the closed-form `alpha_j(h)` of
§1.2 rather than as separate unknowns. `Sum z_i C_i` includes any strong acid
or base added as a titrant (a signed concentration contributing `+C_HCl` or
`-C_NaOH` to the residual) — this is how the titration curve (§5, V1) is
generated: sweep the titrant term, resolve `h` at each point.

**This is the single nonlinear equation solved per cell, per step** (D2).
`g` is strictly monotonically DECREASING in `h > 0` for any physically valid
input (each `alpha_j(h)` term is monotone in `h` by construction, and `h` is
monotone increasing, `-Kw/h` is monotone increasing — sum of monotone
non-decreasing minus a non-decreasing acid term is not obviously monotone in
general polyprotic systems, so the model-specific monotonicity proof is
deferred to the Newton bracketing argument in §3.3, which relies only on
sign-consistency of `g` at the bracket endpoints, not global monotonicity).

### 1.4 The minimal carbonate-only case needs NO extra transported species (ties to D9)

For the reference validation system (§5), the "inert ion" and "titrant"
terms in (1) are the ONLY extra input besides the single transported total
`C_T` (dissolved inorganic carbon). A strong-acid/base titration adds a
scalar (not even a transported field — a scenario-level scalar ramp or a
second W-SCAL total if spatial titrant transport is wanted) to the residual.
**Zero new transported species are required for V1–V4 (§5)** — confirming
D9's claim that this spec's MVP is unblocked by the multicomponent spec.

---

## 2. Equilibrium-system data structure (decision D1, D9, D10)

### 2.1 Configuration object — species-agnostic (Rust API)

```rust
/// One acid/base dissociation family (e.g. carbonate, ammonia, phosphate).
/// Agnostic: holds only pKa values, a name for I/O labeling, and the index
/// of the transported total concentration it consumes. No species-specific
/// branching anywhere downstream of this struct — see the Rule-2 ban-list
/// grep in the DoD (§7).
pub struct AcidFamily {
    /// Human-readable label for I/O only (e.g. "carbonate"). Never matched
    /// on in solver logic — see ban list, CLAUDE.md/physics-discipline Rule 2.
    pub label: String,
    /// pKa_1..pKa_n, protonation steps most-protonated-first.
    pub pka: Vec<f64>,
    /// Index into the scenario's transported-total list (§2.3) supplying
    /// C_T for this family.
    pub total_index: usize,
}

/// A single inert (non-equilibrating) charged species contributing to
/// charge balance but not itself part of any AcidFamily (e.g. Na+, Cl-).
/// Concentration may be a scenario constant or (future, D9-b) a transported
/// total; phase 1 supports the constant form only (§2.3 note).
pub struct InertIon {
    pub label: String,
    pub charge: i32,       // z_i, signed
    pub concentration: f64, // mol/L, scenario-constant in phase 1
}

/// Full equilibrium-system configuration for one cell-local solve.
/// `None` at the scenario/solver level => zero equilibrium code runs,
/// bit-identical to pre-T5b (D10, B-6 discipline).
pub struct EquilibriumSystem {
    pub families: Vec<AcidFamily>,
    pub inert_ions: Vec<InertIon>,
    /// pKw, temperature-dependent; default 14.00 @ 25 C (configurable per
    /// REQ §2 dimensionless/physical-unit conventions — CODATA value with
    /// its own citation if changed).
    pub pkw: f64,
    /// Newton convergence tolerance on the charge-balance residual g(h),
    /// same units as concentration (mol/L). See §3.2 for the default and
    /// its derivation.
    pub newton_tol: f64,
    /// Hard iteration cap (diagnostic trip, not a silent fallback — a cell
    /// that fails to converge within this cap is a validation ERROR, not a
    /// defaulted pH; see §3.4).
    pub newton_max_iter: u32,
    /// Ideal (phase 1 default) or Davies ionic-strength-corrected activity
    /// (§6.5, opt-in, independently validated — never the silent default).
    pub activity_model: ActivityModel,
}

pub enum ActivityModel {
    Ideal,
    Davies { a_debye: f64 }, // A ~ 0.509 (25 C, water) — citation in §6.5
}
```

`total_index` and `InertIon::concentration` are the ONLY points where a
scenario supplies numbers; the solver code (§3) never contains a carbonate-
or ammonia-specific branch — it only ever iterates `families: &[AcidFamily]`
and `inert_ions: &[InertIon]` generically. This satisfies D1 and the Rule-2
ban-list requirement mechanically (grep for `"carbonate"` or `"ammonia"` in
`crates/lbm-core/src/*.rs` must return zero hits outside doc comments/tests).

### 2.2 Per-cell solve entry point (Rust API)

```rust
/// Solve the proton condition (1) for one cell given its transported totals
/// and (constant) inert-ion charges. Pure function: no lattice/backend
/// dependency, no I/O — this is what makes it embeddable identically in
/// CpuScalar and (later) a GPU per-cell kernel (§6.6).
///
/// `totals[i]` corresponds to `families[i].total_index`'s C_T for THIS cell
/// (already resolved by the caller from SoaFields conc-family storage,
/// per the D9 storage-agnostic contract).
///
/// Returns the converged `h = [H+]` plus per-family alpha fractions, or an
/// error variant (never a silently defaulted value — D10/Rule-2 "no
/// silent fallback").
pub fn solve_proton_condition<T: Real>(
    system: &EquilibriumSystem,
    totals: &[T],
) -> Result<Speciation<T>, EquilibriumError>;

pub struct Speciation<T> {
    pub h: T,                       // [H+]
    pub ph: T,                      // -log10(h)  (activity-based if Davies active, §6.5)
    pub oh: T,                      // Kw / h
    pub alpha: Vec<Vec<T>>,         // alpha[family_idx][step_idx]
    pub newton_iters: u32,          // diagnostic
    pub residual: T,                // final |g(h)|, diagnostic
}

pub enum EquilibriumError {
    /// Newton did not converge within newton_max_iter — a hard error
    /// surfaced to the caller (validation failure / diagnostic field),
    /// never masked (Rule 2 "no silent fallback").
    NotConverged { iters: u32, residual: f64 },
    /// A supplied total or ion concentration was negative beyond
    /// round-off — an upstream data error, not something this function
    /// papers over.
    NegativeInput { index: usize, value: f64 },
}
```

### 2.3 `SoaFields` mount point (decision D8, D10)

```rust
/// Derived pH field, compact core. `None` <=> no equilibrium system
/// configured — zero allocation, bit-identical legacy path (B-6).
pub ph: Option<Vec<T>>,
/// Per-family, per-step speciation fractions (alpha), compact core,
/// flattened [family][step][cell]. Diagnostic/output only — never read
/// back into transport, viscosity, or density in phase 1 (D8).
pub speciation: Option<Vec<Vec<T>>>,
```

Mirrors the `Option`-gated placement precedent of `force_field`/`omega_field`
(`fields.rs:196,199`) and the WSCAL `h`/`htmp`/`conc` triad (WSCAL §3.2).
**Scale-out note (ties D9):** if/when ≥2 independent transported totals are
needed (e.g. carbonate + ammonia sharing one charge balance), each total is
either a separate WSCAL single-component `conc` slot (today, mechanically)
or a `Vec<ScalarField>` element (once multicomponent lands) — `solve_proton_
condition`'s `totals: &[T]` signature is unchanged either way (§2.2 is
storage-blind by design).

---

## 3. The algebraic solver (Newton) + operator-split placement

### 3.1 Newton iteration in log-space (decision D3)

Let `x = ln h`. Rewrite (1) as `G(x) = g(e^x) = 0`. Newton update:

```
x_{k+1} = x_k - G(x_k) / G'(x_k),     G'(x) = g'(h) * h    (chain rule, h = e^x)
```

`g'(h)` is available in closed form by differentiating (1) term-by-term
(each `alpha_j(h)` is a rational function of `h`, differentiable in closed
form — no numerical differencing, which would be both slower and a source
of avoidable error). This is standard Newton–Raphson on the log-transformed
proton condition (Morel & Hering ch. 3.5; PHREEQC's core solver, Parkhurst &
Appelo 2013 ch. 2, uses the equivalent "-log(activity)" internal variables
for exactly this conditioning reason).

### 3.2 Convergence criterion (decision D7)

Converged when `|g(h_k)| < newton_tol` (absolute residual on the ORIGINAL
charge-balance equation (1), in mol/L — not on `x`, because the physically
meaningful quantity is charge-balance closure, §5 V2). **Default `newton_tol
= 1e-12` mol/L** for f64 (matches WSCAL V4's `1e-12` round-off mass-
conservation gate — the same numerical floor rationale: charge balance
should close to round-off, not to a band). For f32 cell-local solves
(matching a CpuSimd/GPU f32 build), `newton_tol = 1e-6` mol/L (matching
WSCAL V4's f32 floor). **This tolerance is a numerical-convergence
parameter, not a physical band** — it is reported per run (§5), never
tuned to pass a downstream metric (Rule 2, "calibrated constant" ban does
not apply to a solver stopping criterion, but the SAME discipline of "no
band-fitting" applies: tightening the tolerance must only improve accuracy
monotonically, verified by the order test V3).

### 3.3 Initial guess and the bracketing safeguard (NOT a clamp — see distinction)

Newton on a log-transformed monotone-ish function converges quadratically
from a good start but can overshoot for a bad one. The safeguard is a
**bracketing bisection fallback**, standard globalized Newton (Press et al.,
*Numerical Recipes* §9.4, "Newton-Raphson with bisection backup"):

1. Establish a bracket `[h_lo, h_hi]` where `g(h_lo) > 0` and `g(h_hi) < 0`
   (guaranteed to exist because `g(h) -> +inf` as `h -> inf` — the `h` term
   dominates — and `g(h) -> -inf` as `h -> 0` — the `-Kw/h` term dominates;
   this is a mathematical fact about (1)'s asymptotics, not a tuned bound).
   `h_lo = 1e-14 mol/L` (below pH 14), `h_hi = 1 mol/L` (pH 0) are safe
   universal starting brackets for aqueous systems (`Kw` and all `Ka` values
   are constrained to the aqueous range by construction — a system outside
   pH [0,14] is outside this spec's validity domain, §6.4, and is a
   configuration error, not silently handled).
2. Each Newton step that would leave `[h_lo, h_hi]` is replaced by a
   bisection step instead (standard safeguard, not a clamp on the SOLUTION —
   it bounds the ITERATE search interval using a mathematically guaranteed
   sign-change bracket, shrinking monotonically as convergence narrows it;
   contrast the banned "position clamp" pattern, which caps a *physical*
   transported quantity at an arbitrary bound to hide transport error. Here
   the bound is exact bracketing of a root whose existence and uniqueness
   in `(0, infinity)` is guaranteed by the intermediate value theorem plus
   the asymptotic argument above — no error is being hidden, and the
   bracket does not change the converged answer, only the path to it).
3. **Warm start**: the previous step's converged `h` for this cell is the
   default initial guess (Newton converges in 2-4 iterations typically when
   `h` changes smoothly step-to-step); the bracket midpoint is the fallback
   if no previous value exists (first step) or if the warm-started Newton
   diverges outside the bracket.

### 3.4 No silent fallback on non-convergence (Rule 2)

If Newton (with bisection backup) does not converge within `newton_max_iter`
(default 50 — bisection alone halves the bracket each step and the bracket
width `1 mol/L` needs ~47 halvings to reach `1e-14` mol/L resolution, so 50
is a generous cap derived from the bracket width, not an arbitrary number),
`solve_proton_condition` returns `EquilibriumError::NotConverged`. The
solver-level caller surfaces this as a per-run diagnostic count (cells not
converged) and — for phase 1 — treats ANY non-convergence as a hard run
error (abort with the diagnostic), never a defaulted pH. This is deliberately
strict (matching WSCAL §6.5's "positivity is a diagnostic, not a clamp"
posture) because a silently-defaulted pH would be a physical-integrity
violation invisible to every downstream gate.

### 3.5 Operator-split placement relative to transport and reaction (decision D4)

Per solver step, when `EquilibriumSystem` is configured:

```
1. HYDRODYNAMIC f STEP (unchanged, WSCAL/REQ invariant order).
2. SCALAR ADE SUB-STEP (WSCAL §4.2, unchanged): transport each configured
   total C_T,k (e.g. total carbonate) via its own h_k/htmp_k/conc_k D3Q7
   distribution (or the future multicomponent equivalent, D9).
3. (if W-REACT kinetic reactions are also configured, §8): apply the
   kinetic reaction increment to the SAME totals — kinetics and fast
   equilibrium act on the same conservative pool.
4. EQUILIBRIUM PROJECTION (this spec, NEW, runs LAST):
   for each fluid cell: totals[] <- read conc_k for each configured family
                        (h, ph, alpha, oh) <- solve_proton_condition(system, totals)
                        write ph, speciation into SoaFields
   (C_T,k itself is NOT modified — decision D5; only derived fields change.)
```

**Split-error argument**: this is a first-order Lie/Godunov operator split
between the "transport+kinetics" differential part and the "equilibrium"
algebraic projection — standard for index-1 DAEs (Hairer & Wanner 1996 ch.
VII.2; Strang splitting is a natural second-order refinement, noted as a
phase-2 option in §6.6, not required for the MVP band, §5 V3). The split
error vanishes as `dt -> 0` at the expected first order for the chosen
(non-Strang) splitting; this is exactly what the dt-halving test (§5, V3)
verifies.

---

## 4. pH output derivation

### 4.1 Definition

```
pH = -log10(a_H+)
```

`a_H+ = h` (ideal, phase 1 default, D6) or `a_H+ = gamma_H+ * h` (Davies
activity correction, §6.5, opt-in). The `Speciation::ph` field (§2.2) always
reports the ACTIVITY-based pH consistent with whichever `ActivityModel` is
configured — this is the standard definition of "pH" in aqueous chemistry
(IUPAC; Stumm & Morgan ch. 3.1) and is why the ideal/Davies choice must be
explicit rather than defaulted silently (D6): the two produce numerically
different pH for the identical `h` once ionic strength is nonzero.

### 4.2 Output wiring

`ph: Option<Vec<T>>` mounts into the same field-export path as `conc`
(WSCAL P3 precedent) — VTI/manifest output (per
`lbmflow-user-postprocess` conventions), no new export mechanism.
`speciation` (per-family, per-step fractions) exports as additional named
scalar fields, one per `(family, step)` pair, labeled from `AcidFamily::
label` (I/O-only use of the label — §2.1 ban-list note).

---

## 5. Validation plan → VR-STR-04 extension (T5b acid–base sub-row)

Tests are authored adversarially by codex/Opus from THIS spec, in a worktree
disjoint from the implementation worktree (CLAUDE.md convention). Bands are
provisional MVP gates (tightening always allowed; loosening needs a recorded
PHYSICS.md rationale, per T17 band governance).

| ID | Test | Metric & band | Grid / config | VALIDATION.md row |
|---|---|---|---|---|
| **V1** | **Titration curve vs analytic Bjerrum/carbonate solution.** Single well-mixed cell (0D, `u=0`, no transport needed — isolates the algebraic closure from ADE error), sweep added strong acid (HCl) concentration `C_a` from 0 to `2*C_T` at fixed `C_T` (total DIC = 2 mmol/L, a typical freshwater value — Stumm & Morgan Table 4.1). Compute `pH(C_a)` from `solve_proton_condition` at each point; compare to the closed-form analytic titration curve (charge balance (1) solved to machine precision by an independent reference implementation — e.g. `scipy.optimize.brentq` on the same equation, computed OUT of this codebase as the reference oracle). | (a) `max\|pH_computed - pH_reference\|` over the sweep **< 1e-6 pH units** (f64) / **< 1e-3** (f32) — this is an algebra-implementation-correctness gate, not a physical-model band, so it is tight; (b) **behavior anchor**: `pH` is **strictly monotonically non-increasing** in `C_a` (physical: adding acid never raises pH) — assert `pH[i+1] <= pH[i] + eps` across the sweep. | 0D single-cell harness, carbonate `pKa1=6.35, pKa2=10.33`, `C_T=2e-3 mol/L`, `pKw=14.0`, ideal activity | **VR-STR-04b** (new sub-row) |
| **V2** | **Charge balance residual -> round-off.** For every converged cell in V1 and V3, evaluate `g(h)` (equation (1)) directly from the converged `h` and the same inputs. | `\|g(h)\|` **< newton_tol** by construction (this is the convergence criterion itself, §3.2) AND **< 1e-12 mol/L** independent of the configured `newton_tol` when `newton_tol` is set to its f64 default — i.e. verify the DEFAULT closes to round-off, not just to whatever tolerance was configured (a tautology-avoidance check: confirms `newton_tol=1e-12` is actually achievable, not just requested). | Same as V1/V3 | **VR-STR-04b** |
| **V3** | **DAE/operator-split convergence (dt-halving).** A spatially-transported case: 1D channel, uniform `u`, a step change in total-carbonate inlet concentration, WITH a slow kinetic perturbation active (a synthetic first-order relaxation toward a shifted equilibrium, standing in for W-REACT until it lands — §8) so the split error is nonzero. Run at `dt`, `dt/2`, `dt/4`; measure `pH` profile L2 error vs a `dt/8` reference (Richardson extrapolation, same convention as WSCAL V1's grid-refinement order test). | **Order >= 0.9** (first-order Lie split expected order = 1; MVP band allows some slack for the algebraic-projection nonlinearity, tighten after characterization) — error ratio between successive halvings >= ~1.8 (vs the ideal 2.0 for order 1). | 1D-in-3D channel `256x4x4`, periodic transverse, carbonate system, synthetic kinetic term rate `k` s.t. Da = O(1) (deliberately NOT infinite, so the split error is measurable — the fast-equilibrium limit Da->infinity is validated separately by ABLATION A2) | **VR-STR-04b** |
| **V4 (ablation)** | **Equilibrium OFF => totals inert, no speciation.** Same channel as V3 but `EquilibriumSystem = None`. | `ph` and `speciation` fields are unallocated (`None`) — a compile/API-level check, not a numeric one; AND the transported `conc_k` totals are **bit-identical** to a pre-T5b WSCAL-only run of the same scenario (the equilibrium projection must not perturb the totals it reads, decision D5, verified by literal bit-identity, not a tolerance). | Same grid as V3, `EquilibriumSystem=None` | **VR-STR-04b** (ablation guard, physics-discipline Rule 3 requirement) |
| **V5 (Bjerrum structure)** | **Speciation fractions cross at the pKa points** — the defining structural feature of a Bjerrum plot (Stumm & Morgan Fig. 4.1-style). At `pH = pKa1`, assert `alpha_0(h) == alpha_1(h)` (both `H2CO3*` and `HCO3-` fractions equal, `= 0.5` each, `alpha_2` negligible) to solver tolerance; at `pH = pKa2`, assert `alpha_1(h) == alpha_2(h)`. | `\|alpha_0 - alpha_1\|` at `h = 10^-pKa1` **< 1e-10** (algebraic identity, tight); same for the `pKa2` crossing. **This is the mandatory "speciation fractions cross at the pKa points" behavior anchor** named in the task brief — a metric passing V1's L2 band does NOT by itself prove the crossing structure is right (physics-discipline Rule 3, two-layer gate). | 0D single-cell harness, same carbonate system as V1 | **VR-STR-04b** |

**Mandatory negative/consistency tests** (physics-discipline Rule 3 +
REQ §8 pattern):

- **Wrong-sign charge test**: a mutant that flips the sign of one `alpha_j`
  charge contribution in (1) must FAIL V1's monotonicity anchor (proves the
  charge-balance signs are load-bearing, not decorative).
- **Non-convergence surfacing test**: a deliberately pathological input
  (e.g. `pKa` values outside [0,14], §3.3 validity bound) must return
  `EquilibriumError`, not a silently-produced `pH` — asserted directly
  against the public API (§2.2), not inferred from a downstream symptom.
- **`total_index`/`InertIon` agnosticism test**: running the SAME solver
  code path with the carbonate `pka` values swapped for a synthetic
  2-step acid with different `pKa`s (no other code change) reproduces V1's
  analytic-vs-computed agreement — proves D1's agnosticism is real, not
  aspirational (grep-based Rule-2 check is necessary but not sufficient;
  this test is the behavioral proof).

### 5.6 Behavior-validity review (mandatory, REQ / CLAUDE.md)

After each validation run, before reporting: (a) the titration curve (V1) is
S-shaped with inflections AT the configured `pKa`s, not offset — an offset
inflection would indicate a sign or stoichiometry error in (1) even if the
L2 band happens to pass elsewhere on the curve; (b) the Bjerrum crossing
(V5) fractions sum to 1 at every `pH` sampled, not just at the crossing
points (spot-check `alpha_0+alpha_1+alpha_2==1` across the full V1 sweep,
§1.2's algebraic identity); (c) the dt-halving profile (V3) shows the error
concentrated where the kinetic perturbation and equilibrium projection are
BOTH active (spatially, near the concentration front) — a uniform-everywhere
error would indicate the split error is dominated by the ADE discretization
itself (already gated by WSCAL) rather than the equilibrium projection this
spec adds, meaning V3 would not actually be testing what it claims; (d) the
ablation (V4) totals-bit-identity check is inspected cell-by-cell, not just
as an aggregate norm (a norm can hide a canceling sign error). Record the
review in PHYSICS.md (entry template below, §7.4) or the track's findings
file.

---

## 6. Solver convergence / validity domain

### 6.1 Existence and uniqueness of the proton-condition root

For any physically valid input (`C_T,k >= 0`, `Kw > 0`, `K_aj > 0`), `g(h)`
in (1) is continuous on `(0, infinity)`, `g(h) -> +infinity` as `h ->
infinity` (the `h` term dominates all bounded `alpha_j C_T` and `-Kw/h`
terms) and `g(h) -> -infinity` as `h -> 0+` (the `-Kw/h` term dominates).
By the intermediate value theorem a root exists; for the single-unknown
proton-condition reformulation this root is the unique physical pH (Morel &
Hering ch. 3, the standard existence argument for TOTH-form equations).
Uniqueness in general polyprotic multi-family systems is not a universal
theorem (pathological configured `pKa` sets could in principle admit
multiple roots) — phase 1 relies on the bracketing bisection (§3.3) to find
A root and additionally asserts monotonic charge balance as a validity
check (any configured system failing global monotonicity is flagged, not
silently accepted); this is the stated boundary of phase-1 validity, not a
gap papered over.

### 6.2 Newton/bisection convergence order

Newton (unsafeguarded) converges quadratically near the root (standard
Newton theory, given `G` is `C^2` — true here, all `alpha_j` terms are
smooth rational functions of `h` away from `h=0`). The bisection fallback
converges linearly (halves the bracket each step) but only engages when
Newton would overshoot — expected behavior is 1-3 bisection steps (if any)
followed by quadratic Newton convergence, empirically characterized and
reported per configured system in PHYSICS.md (§7.4) at implementation time.

### 6.3 Grid/timestep coupling to WSCAL

The equilibrium projection's OWN error (Newton residual, §3.2) is
independent of `dt`/`dx` — it is a per-cell algebraic solve to a fixed
numerical tolerance. The COUPLED error (§5 V3) comes entirely from the
operator split between transport/kinetics and the projection, which DOES
scale with `dt` (§3.5). The WSCAL transport error (grid-Péclet, `tau_s`
window, WSCAL §6.1-6.2) is unaffected by this spec — the equilibrium
projection reads `conc_k` but does not feed back into the ADE (D8), so
WSCAL's own stability bounds apply unchanged.

### 6.4 Ideal-activity validity domain (decision D6)

`a_i = [i]` (activity = concentration) is exact only in the infinite-
dilution limit `I -> 0` (ionic strength). For typical freshwater/mild
process-water conditions (`I <~ 0.01-0.05 M`), the ideal approximation's
error in `pH` is small (a few percent in `K_a`-shifted equilibrium
constants, Davies 1962 characterization) but grows systematically with `I`.
**Stated validity domain: `I <= 0.01 M`** for the ideal default to be within
the V1 band's implicit accuracy (this is an activity-model validity bound,
analogous to WSCAL's `Pe_Delta` bound — outside it, switch `ActivityModel`
to `Davies`, never silently accept degraded accuracy). `pH` range
`[0, 14]` (aqueous window, §3.3); temperature fixed at the configured `pKw`
(no explicit temperature-dependence model for `K_aj` in phase 1 — each
configured `pKa` is assumed measured at the run's nominal temperature; a
temperature-dependent `pKa(T)` correlation is an explicit future extension,
not silently assumed constant across a non-isothermal run — if W-SCAL's
active-scalar thermal axis (REQ §10 "Thermal axis recommended as API-
reserved extension") lands with a temperature field, this spec's `pKa`
inputs must be revisited, flagged here as a forward dependency, not solved).

### 6.5 Davies ionic-strength correction (opt-in extension, decision D6)

```
log10(gamma_i) = -A * z_i^2 * ( sqrt(I)/(1+sqrt(I)) - 0.3*I ),   I = 0.5 * Sum_i z_i^2 [i]
```

`A ~ 0.509` (25 C, water; Davies 1962, *Ion Association*, extending Debye-
Hückel to `I` up to ~0.5 M with the empirical `-0.3I` term). This is a
LITERATURE-BACKED CLOSURE per Rule 1 and requires, before it may be enabled
in any run whose results are reported: (1) this citation + the derivation
note above (done here); (2) the stated validity domain `I <~ 0.5 M` (Davies'
own stated range); (3) its OWN validation test — a Davies-vs-ideal
divergence test at a known `I` matching a published activity-coefficient
table (e.g. Stumm & Morgan Table 3.4) to within the published table's
precision; (4) a PHYSICS.md entry (template §7.4) separate from the ideal-
default entry. **Phase 1 ships the `ActivityModel::Davies` code path
API-complete but its OWN validation test is a follow-on order (§7, O3
scope note)** — it must not be reported as validated until that test lands
and passes; until then any run using it must disclose that it is running
an unvalidated activity path.

### 6.6 GPU / performance note (deferred, not blocking)

The per-cell Newton solve (§3) is branchy and iteration-count-variable —
unfriendly to naive GPU SIMT (divergent iteration counts across a warp).
Phase 1 is CPU-only (`CpuScalar` reference, mirroring WSCAL P8's staging
posture); a GPU path would need either a fixed-iteration-count Newton
(pad to `newton_max_iter` always, wasting some cycles) or a warp-uniform
bisection-only variant — both are follow-on performance work, explicitly
out of this spec's scope, and must not silently change the converged
answer (any GPU variant needs its own equivalence gate against the CPU
oracle, analogous to `backend_simd_equiv.rs`).

---

## 7. Codex order breakdown

Three orders, file-conflict-aware, one order = one bundle = one dedicated
worktree (CLAUDE.md convention). Implementation and adversarial-test orders
never share a worktree.

| Order | Scope | Primary files (conflict boundary) | DoD |
|---|---|---|---|
| **O1 — Equilibrium closure + Newton solver** | `EquilibriumSystem`/`AcidFamily`/`InertIon`/`ActivityModel` config types (§2.1); `Speciation`/`EquilibriumError` (§2.2); `solve_proton_condition` (charge balance (1), alpha-fraction closed forms §1.2, Newton-in-log-space + bracketing bisection §3.1/3.3); ideal activity (D6) with the `Davies` variant API-complete but its OWN validation deferred to a follow-on order per §6.5; the `ph`/`speciation` `SoaFields` slots (§2.3) and the solver-level projection pass wired AFTER the WSCAL scalar sub-step and any W-REACT kinetic sub-step (§3.5) — `Option`-gated, `None` = zero code executes (D10). | NEW module (e.g. `crates/lbm-core/src/equilibrium.rs`) for the pure-function solver + config types; small additive edits to `fields.rs` (two `Option` fields, same pattern as WSCAL §3.2) and `solver.rs` (one new pass call site, placed last in the per-step order). **No edits to `lattice.rs`, `kernels.rs`, or `backend.rs`** — this closure needs no new distribution, no new BC, no new lattice. | `solve_proton_condition` unit tests (V1 numerics, V5 Bjerrum-crossing identity, V2 residual round-off) green on `CpuScalar`; `EquilibriumSystem=None` bit-identical to pre-T5b on every existing preset (ablation half of V4, the compile/API-level `None` check); Rule-2 ban-list grep clean (no species-name branch in `equilibrium.rs`/`solver.rs`/`fields.rs`); PHYSICS.md entry (§7.4) landed with O1. |
| **O2 — Scenario schema for equilibrium systems + pH output** | `crates/lbm-scenario` JSON schema additions: an `equilibrium` block (families with `pka`/`label`/`total_index` referencing the scenario's existing scalar-total declarations from WSCAL's O2 schema work — reuse, not duplicate, WSCAL's `inlet_C`/`D` scalar block; §0 D9 storage-agnostic contract makes this a schema-level reference, not a new transport declaration), `inert_ions`, `pkw`, `newton_tol`/`newton_max_iter`, `activity_model`; validation warnings for `I > 0.01` under `Ideal` (§6.4) and for configured `pKa` outside `[0,14]` (§3.3); CLI/output wiring for `ph`/`speciation` fields (VTI/manifest, mirrors WSCAL O2's `conc` output). The synthetic first-order kinetic perturbation needed ONLY for V3's split-error test is scenario-expressible here too (a minimal stand-in, explicitly labeled as a test fixture, not a W-REACT feature — §8 note). | `crates/lbm-scenario/src/lib.rs` (schema + validation), `crates/lbm-cli` (ph/speciation output) — **disjoint files from O1**. | V1/V3's channel-scenario JSON validates and runs end-to-end through the CLI producing `ph` output files; grid-Péclet-style warnings fire on the two new validity-domain checks (§6.4, §3.3) without false-positiving on the V1-V5 reference configs. Depends: O1 (needs the `EquilibriumSystem` types to (de)serialize into). |
| **O3 — Validation authorship (codex adversarial, separate worktree)** | All of §5 (V1-V5) + the three negative/consistency tests (wrong-sign, non-convergence surfacing, agnosticism). Authored from THIS spec, not from the O1/O2 implementation. Includes the external reference-oracle computation for V1 (an independent script, e.g. Python/`scipy.optimize.brentq`, checked into the test fixture directory as a documented oracle-generation script, NOT re-derived by hand each time). | `crates/lbm-core/tests/acidbase_*.rs` + `crates/lbm-scenario/tests/*` (new files only — no impl-file conflict) + a `scripts/` or `tests/fixtures/` oracle-generation script. | Tests compile red against a stub, go green against O1/O2 as they land; freeze the VR-STR-04b bands in VALIDATION.md (this spec's provisional bands are MVP; O3 characterizes and the PM freezes per T17 band governance). Runs concurrently with O1/O2 from the start. |

**Critical-path ordering:** O1 -> O2 (schema needs the O1 types). O3 runs
concurrently from the start (test worktree, stubs against the spec's public
API signatures in §2.1/§2.2, which are frozen by this document). No
dependency on a not-yet-written `WSCAL_MULTICOMPONENT_SPEC.md` or `WREACT_
IMPL_SPEC.md` blocks O1/O2/O3 (§0 D9, §8) — the reference validation system
(carbonate, §5) needs exactly one already-landed WSCAL total.

### 7.4 The PHYSICS.md validity-domain statement (mandatory entry text)

The O1 order must add, on landing, a PHYSICS.md entry containing:

> **Acid–base fast-equilibrium speciation & pH — algebraic DAE closure
> (Morel & Hering 1993 ch. 3; Stumm & Morgan 1996 ch. 3).** Given transported
> total concentrations `C_T,k` (produced by W-SCAL, unmodified by this
> closure — decision D5), solve the charge-balance proton condition `g(h)=0`
> (equation (1)) per cell via Newton-in-log(h) with a bracketing-bisection
> safeguard (§3.1/§3.3) to recover `[H+]`, `pH=-log10(a_H+)`, `[OH-]=Kw/h`,
> and closed-form speciation fractions `alpha_j(h)` (§1.2) for each
> configured acid/base family (`pKa` set + charge vector, agnostic — the
> carbonate system `pKa1=6.35/pKa2=10.33`, Plummer & Busenberg 1982, is a
> configuration instance, not a hardcoded species). Ideal activity
> (`a_i=[i]`) is the phase-1 default; **validity domain `I<=0.01 M`,
> `pH in [0,14]`, isothermal (pKa fixed at run temperature)**; outside `I<=
> 0.01M`, switch to the Davies correction (§6.5, opt-in, independently
> validated separately from this entry — do not report Davies results
> under this entry's validation). Newton convergence tolerance
> `newton_tol=1e-12 mol/L` (f64) / `1e-6` (f32), reported per run, is a
> numerical-solver parameter, not a physical band. **Non-convergence is a
> hard error (`EquilibriumError::NotConverged`), never a defaulted pH.**
> Operator-split placement: LAST pass in the per-step order, after WSCAL
> transport and any W-REACT kinetics (§3.5) — a first-order Lie/Godunov
> split (Hairer & Wanner 1996 ch. VII.2), split-error convergence order
> characterized and reported (target >= 0.9, VR-STR-04b V3). **Validation**:
> `crates/lbm-core/tests/acidbase_*.rs` V1 (titration vs analytic, L2
> < 1e-6 pH f64), V2 (charge-balance residual round-off), V3 (dt-halving
> split-error order), V4 (equilibrium-off ablation: totals bit-identical),
> V5 (Bjerrum pKa-crossing identity). **Replaces/interacts with**: reads
> W-SCAL's transported totals (WSCAL PHYSICS.md entry) without modifying
> them; is superseded (for infinitely-fast reactions) or coexists (for
> finite-rate reactions, §8) with W-REACT once that spec lands.

---

## 8. Coexistence with W-REACT (equilibrium vs kinetic reactions)

`docs/proposals/WREACT_IMPL_SPEC.md` does not exist at the time of this
spec (checked: absent from `docs/proposals/`). This section states the
reconciliation contract T5b commits to so that when W-REACT is specified,
neither side needs to be redesigned:

- **Physical distinction**: fast-equilibrium (this spec) is the `Da ->
  infinity` limit of a general reversible reaction `R_k(C)` (REQ §3's
  reaction source term) — the forward/backward rates are so much faster
  than transport/mixing that the reaction is always at local equilibrium,
  collapsing the ODE `dC/dt = R(C)` to the algebraic constraint `R(C) = 0`
  solved this spec's way (§1.3). W-REACT's finite-rate kinetics are the
  general case where `Da` is O(1) or smaller and the ODE must be integrated
  (not collapsed) — e.g. an explicit or implicit sub-step advancing `C_k`
  by `dt * R_k(C)` between the transport sub-step and (if ALSO configured)
  this spec's equilibrium projection.
- **Operator-split ordering (commitment, §3.5 already states the mechanics)**:
  transport -> kinetics (W-REACT, if configured) -> equilibrium projection
  (this spec, if configured) -> end of step. Rationale: kinetics changes the
  totals (mass-action source/sink between DIFFERENT total pools, e.g. a
  reaction consuming ammonia total while producing nitrate total); the
  fast-equilibrium projection then re-partitions EACH pool's own internal
  speciation (e.g. total-ammonia's `NH4+`/`NH3` split) at fixed post-
  kinetics totals. The two operate on different axes (inter-pool for
  kinetics, intra-pool for equilibrium) and therefore commute in the sense
  that swapping their order changes only the split error (§3.5), not the
  physical fixed point — but the stated order (kinetics-then-equilibrium)
  is the one this spec validates (§5 V3) and is the one to implement first.
- **A reaction rate that itself depends on `pH`** (e.g. an acid-catalyzed
  kinetic step) is explicitly a W-REACT-phase concern: W-REACT would READ
  this spec's `ph` output field (§4) as an input to its rate law, the same
  read-only relationship WSCAL's `conc` has to this spec's totals. This
  spec's `ph` field must therefore be computed and available BEFORE
  W-REACT's kinetic sub-step can consume it in a future coupled scheme —
  flagged here as a probable reordering need once W-REACT lands (kinetics
  needing `pH` would want equilibrium computed on the PREVIOUS step's
  totals, a one-step-lagged coupling, analogous to how WSCAL itself reads
  velocity from the current step rather than lagging — the exact resolution
  is deferred to WREACT_IMPL_SPEC's own §4-equivalent pass-structure
  section, not decided here).
- **No code overlap**: W-REACT's kinetic sub-step and this spec's
  equilibrium projection are, by the ordering above, DISJOINT pass slots
  (like WSCAL/WVOF's disjoint pre-pass/post-pass slots, WSCAL §7.1) — when
  WREACT_IMPL_SPEC is written, its file-conflict table should list
  `equilibrium.rs`/the projection call site in `solver.rs` as "adjacent, not
  overlapping" with its own kinetics module, by the same reasoning WSCAL
  used against WVOF.

---

## 9. Load-bearing code references (grounding index)

| Claim | File:line |
|---|---|
| `SoaFields` struct, `force_field`/`omega_field` `Option` precedent (mount pattern for `ph`/`speciation`) | `crates/lbm-core/src/fields.rs:168-210`, `:196-199` |
| WSCAL `h`/`htmp`/`conc` `Option` triad (the transported totals this spec consumes) | `docs/proposals/WSCAL_PASSIVE_SPEC.md` §3.2 |
| WSCAL scalar sub-step slot (solver-orchestration level, after `update_moments`) — this spec's projection runs immediately after it | `docs/proposals/WSCAL_PASSIVE_SPEC.md` §4.2 |
| REQ scalar/reaction governing forms (`R_k(C)`, the reaction source this spec's equilibrium is the `Da->infinity` limit of) | `docs/REQ_STIRRED_REACTOR.md` §3 |
| REQ §8 mandatory test menu ("active-scalar dt-halving convergence (MJ-007)" — the template for this spec's V3) | `docs/REQ_STIRRED_REACTOR.md` §8 |
| REQ §11 DAG: `W-REACT reaction / active feedback` depends on `W-SCAL` | `docs/REQ_STIRRED_REACTOR.md` §11 |
| VR-STR-04 scalar/reaction row (this spec's VR-STR-04b sub-row extends it) | `docs/VALIDATION.md:348` |
| `F_scalar_buoyancy` reserved slot (future active-feedback coupling point; this spec's `ph`/speciation stay OUT of it in phase 1, D8) | `docs/proposals/FORCE_COMPOSITION_SPEC.md` T5 |
| `docs/proposals/WSCAL_MULTICOMPONENT_SPEC.md` | **absent** — confirmed via `docs/proposals/` listing at spec time; this spec's D9 states the unblocked interim path. |
| `docs/proposals/WREACT_IMPL_SPEC.md` | **absent** — confirmed via `docs/proposals/` listing at spec time; §8 states the reconciliation contract to honor once it exists. |

**Literature (decided references):**
Morel, F.M.M. & Hering, J.G. 1993, *Principles and Applications of Aquatic
Chemistry*, Wiley, ch. 3 — the proton-condition (TOTH) single-unknown
reformulation this spec's §1.3/§3 is built on.
Stumm, W. & Morgan, J.J. 1996, *Aquatic Chemistry*, 3rd ed., Wiley, ch. 3-4
— alpha-fraction speciation, Bjerrum-plot structure (§1.2, §5 V5), typical
DIC concentrations (§5 V1).
Plummer, L.N. & Busenberg, E. 1982, *Geochim. Cosmochim. Acta* 46:1011 —
carbonate system `pKa1=6.35`, `pKa2=10.33` at 25 C, `I=0` (the reference
instance, §1.1).
Zeebe, R.E. & Wolf-Gladrow, D. 2001, *CO2 in Seawater: Equilibrium,
Kinetics, Isotopes*, Elsevier, §1.1 — `H2CO3*` lumped-species convention.
Davies, C.W. 1962, *Ion Association*, Butterworths — the Davies
ionic-strength activity correction (§6.5), `A~0.509` at 25 C.
Parkhurst, D.L. & Appelo, C.A.J. 2013, *PHREEQC Version 3* (USGS TWRI Book
6, ch. A43), ch. 2 — Newton–Raphson on log-transformed variables in
geochemical speciation solvers (the precedent for D3's log-space Newton).
Hairer, E. & Wanner, G. 1996, *Solving Ordinary Differential Equations II:
Stiff and Differential-Algebraic Problems*, 2nd ed., Springer, ch. VII —
operator splitting for index-1 DAEs (§3.5, §8's split-ordering argument).
Press, W.H. et al., *Numerical Recipes*, 3rd ed., §9.4 — Newton-Raphson
with bisection safeguard (§3.3).
