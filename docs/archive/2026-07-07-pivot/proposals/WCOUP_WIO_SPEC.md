# W-COUP / W-IO Implementation Specification — Coupling Orchestration + Output/Analysis API

**Document ID**: SPEC-W-COUP-WIO (rev.1, 2026-07-07).
**Scope**: the M-F items `W-COUP coupling loop (FR-COUP)` and
`W-IO I/O & analysis (FR-IO)` of `docs/REQ_STIRRED_REACTOR.md` §11 DAG
(`⊂ MF-ζ`, cross-cutting, incremental across producing subsystems). This is
the **T7 tier**: it does not implement any subsystem's physics — it
*orchestrates* the subcycling/predictor-corrector schedule across whichever
subsystems are landed, and it defines the **reduction/output layer** that
turns their fields into the reactor-engineering quantities REQ §6 names
(concentration, partial pressure, pH, T, gas holdup, `d_32`, `k_L a`,
reaction rates, conversion/yield/selectivity, blend time, RTD).
**Target core**: `crates/lbm-core` (orchestration in `solver.rs`, reductions
as a new `analysis.rs`/`reactor_outputs.rs` module), `crates/lbm-scenario`
(schema extension), `crates/lbm-cli` (output plumbing).
**Acceptance**: VALIDATION.md **T17** rows **VR-STR-05** (coupled regression
+ conservation drift) and **VR-STR-07** (initialization independence), plus
the REQ §8 mandatory active-scalar dt-halving convergence test and the
element-balance/RTD/conversion behavior anchors defined below.

This spec is **executable**: every orchestration decision cites the concrete
sub-step slot in the subsystem spec that owns it, every QOI is a formula +
units + source field, and every gate is mapped to a T17 row with a
provisional band. A follow-on codex implementation order should not need to
re-derive any design decision here.

> **Discipline note (CLAUDE.md prime directive; PHYSICS.md;
> `.claude/skills/lbmflow-physics-discipline`).** This spec adds **zero new
> physics terms** — it is pure orchestration (call order, dt-min, subcycle
> counts) and pure reduction (post-hoc arithmetic on already-computed
> fields). Every formula in §2 is either an exact definition (partial
> pressure from Henry's law, conversion from stoichiometry, RTD from a
> tracer response) or a literature-backed correlation identical to the one
> the owning subsystem spec (WVOF/WSCAL/etc.) already cites — **no new
> citation is introduced here that is not already load-bearing in a sibling
> spec**. The one place this spec could accidentally introduce ad-hoc
> physics — clamping a QOI to a "sane" range (negative rate, `ε_g > 1`,
> yield > 1) — is explicitly banned in §4.5: out-of-range QOIs are
> diagnostics that point at an upstream bug, never clamped.

---

## 0. Summary of decisions (read this first)

| # | Decision | Justification (short) |
|---|---|---|
| C1 | **The coupling loop is a fixed, ordered sequence of solver-level sub-step *slots* — not a new physics model.** Each slot is a call-out to a subsystem's own sub-step (already specified in WVOF_IMPL_SPEC §4.2, WSCAL_PASSIVE_SPEC §4.2, or reserved for WREACT/WBUB/ACIDBASE/THERMAL, none of which are landed yet). W-COUP's job is exclusively: (a) fix the slot order, (b) fix the dt-min rule across slots, (c) fix the predictor-corrector iteration for the active path, (d) provide the `active_scalar_lagged` relaxation flag. | REQ §5 FR-COUP-01 specifies the dataflow order textually; no subsystem spec claims ownership of the *sequencing* — this is the missing orchestration layer the DAG (REQ §11) calls `W-COUP`, dependent on "active set" (whichever of W-VOF/W-SCAL/W-REACT/W-BUB/pH/thermal are landed at integration time). |
| C2 | **Disabled subsystems are no-op slots, not branches.** Each slot is implemented as `if let Some(state) = &self.vof_state { ... }` (mirrors the `Option`-gated `g`/`h` field pattern already frozen in WVOF §3.2 / WSCAL §3.2). With every optional subsystem `None`, `W-COUP`'s orchestration loop degenerates to exactly today's `run_staged_step` call sequence — bit-identical (B-6 invariance). | This is the SAME invariance discipline every sibling spec already carries (WVOF O1 DoD, WSCAL V5); W-COUP inherits it structurally rather than re-deriving it, and it is the acceptance criterion for VR-STR-05's B-6 arm here (§5, test W1). |
| C3 | **dt is a single scalar per solver step, `dt = min(dt_capillary, dt_reaction, dt_buoyancy_Marangoni, dt_particle, dt_hydrodynamic)`** — the REQ FR-COUP-01 constraint list plus the FR-COUP-02 reaction dt plus the active-scalar-feedback.md §5.1 Marangoni/buoyancy dt, taking the **minimum**, never an average or a per-subsystem independent dt (no local time-stepping in phase 1 — REQ does not ask for it and AMR, which would need it, is Phase 2/advanced per REQ §1). | REQ §5 FR-COUP-01 lists the three named constraints "capillary `Δt_σ`, particle `Δt_p`, reaction `Δt_r`" as a set to respect; active-scalar-feedback.md §5.1 is the only place a Marangoni/buoyancy dt formula already exists, decided-not-reinvented here (§1.3). A single global dt (not per-subsystem local stepping) is consistent with the landed engine, which has one `Δt=1` lattice step; W-COUP's dt-min governs how *coarse* that lattice step's physical mapping may be, plus the *subcycle count* for faster-than-hydrodynamic processes (C4), not a second temporal grid. |
| C4 | **Fast processes (reaction, equilibrium/pH) subcycle *within* one hydrodynamic step; slow processes (flow, phase-field, species ADE) do not.** Subcycle count `n_sub = ceil(dt_hydro / dt_fast)`, `n_sub ≥ 1`. Reaction/equilibrium run `n_sub` internal sub-steps of size `dt_hydro/n_sub` while `u`, `φ`, `ρ` are frozen (operator-split, standard Strang-type splitting for stiff chemistry coupled to slow transport). | REQ §5 FR-COUP-01: "For strong coupling / stiff reactions / surface-tension waves: operator-splitting error, subcycling, iterative strong coupling required." FR-COUP-02 requires the reaction solver to switch explicit/implicit/Rosenbrock-BDF *by stiffness* — that switch is WREACT's job (reserved slot, §1.4); W-COUP only owns the subcycle-count/freeze-transport contract around whatever solver WREACT plugs in. |
| C5 | **The active-feedback loop is predictor-corrector by default; `active_scalar_lagged=true` is an explicit relaxation flag, not a silent shortcut.** Predictor: scalar/reaction half-step → property re-evaluation `ρ(C),μ(C),σ(C)[,T]` → force composition → flow step. Corrector: scalar ADE + reaction corrector using the **just-computed** `u,φ` → property re-evaluation → (optional) iterate if the correction changes properties by more than a convergence threshold (`strong coupling`, REQ FR-COUP-01 last clause). Lagged mode skips the corrector (properties feed the *next* step's flow, one step stale) and is gated on a documented stability condition (§1.5) and a lag-error benchmark test (§5, test C1). | REQ §5 FR-COUP-01 verbatim: "Active (fidelity default = predictor-corrector): scalar/reaction predictor → property update → force composition → flow step → scalar ADE corrector → reaction corrector → property re-evaluation (→ optional flow-scalar iteration for stiff coupling). Time-lagged explicit feedback allowed only as flagged relaxation `active_scalar_lagged=true`, with stability conditions + lag-error benchmark." This spec fixes the two loop bodies to the letter of that sentence — no invention. |
| C6 | **The output layer is a pure reduction pass, run on a cadence independent of the step cadence** (`output_every` per QOI, as already the pattern for `OutputSpec.every`), reading only already-computed fields (`SoaFields`, subsystem-owned per-cell state) — it writes NOTHING back into the simulation state. | Read-only reduction cannot introduce a feedback path or a Goodharted metric into the physics; it is exactly the same posture as the existing `FieldKind`→PNG/CSV/VTK pipeline (`crates/lbm-cli/src/runner.rs`), extended with new `FieldKind` variants and a new class of **scalar** (not field) outputs (blend time, RTD, `Np`, conversion) that reduce over the whole domain or a probe set rather than producing a per-cell field. |
| C7 | **New scalar QOIs are a new `ScenarioMetric`/`ReactorOutput` schema class, distinct from `FieldKind` (per-cell) and `ProbeSpec` (per-point time series).** `FieldKind` gets new **per-cell** variants (`Concentration{species}`, `PartialPressure{species}`, `Ph`, `Temperature`, `GasHoldup{definition}`) for visualization/export; whole-domain scalars (blend time, RTD, `Np`, `k_L a`, conversion/yield/selectivity, drift diagnostics) get a new `ReactorOutput` enum written to `manifest.json` under a `reactor_outputs` object, analogous to how `Manifest.diagnostics` already carries `tau`/etc. | `FieldKind` (`crates/lbm-scenario/src/lib.rs:570-589`) is declared `#[non_exhaustive]`-free but every consumer (`runner.rs` match arms at :229, :565, :637, :791) is an exhaustive match — adding variants is mechanical and matches the existing extension pattern (ShearRate/DissipationRate/VorticityMag/QCriterion were added this way). Scalars need a *different* shape (they are not `Vec<f64>` per cell) — reusing `FieldKind` for them would break every exhaustive-match consumer's assumption that a `FieldKind` produces a per-cell `Vec<f64>` (`export_field_f64`, `runner.rs:637`). |
| C8 | **`inlet_phase: gas\|liquid` (already FR-VOF-03/WVOF §5.1) is the ONLY inlet phase declaration; this spec adds NO new raw-φ surface** — species/reaction/equilibrium/thermal inlets are concentration/temperature values riding on top of an `inlet_phase`-typed velocity face, never a φ value. Config validation extends the existing FR-VOF-03 rejection rule to the new fields (a `species_inlet.c` or `thermal.t_inlet` value never implies or requires a raw φ). | REQ FR-VOF-03 "Schema never exposes raw φ for inlets" is already load-bearing (WVOF §5.1, O3); this spec's schema additions (§3) must not create a second door. |
| C9 | **Conservation drift is measured, never enforced.** Every drift monitor in §4 is a read-only accumulator compared against a frozen threshold in a test; nothing in the solver ever rescales a field to force a monitored total back to its initial value (that would be a transport-absorbing clamp on a conserved quantity — banned outright by the physics-discipline ban list). | `.claude/skills/lbmflow-physics-discipline` Rule 2 ban list: "Transport-absorbing clamp... Silently converts transport into accumulation at the bound." A mass-rescaling corrector is exactly this pattern relocated to a global sum instead of a per-cell position; REQ §8's "drift thresholds" language is itself a diagnostic framing (VR-STR-05 "Energy-like = monitoring only, not exactly conserved"), generalized here to every monitored quantity. |
| C10 | **This tier lands LAST**, after at least one of {W-SCAL, W-VOF, W-REACT, pH, W-BUB, thermal} is present, because the orchestration slots and QOI reductions have nothing to orchestrate/reduce over until a producing subsystem exists. Phase 1 of *this* spec's own implementation (§7) targets the **currently-landed set** (W-SCAL passive if/when it lands, W-VOF if/when it lands) and API-reserves the rest; each new subsystem that lands later adds its slot/QOI without restructuring W-COUP/W-IO. | REQ §11 DAG: `W-COUP coupling loop (FR-COUP) | active set | incremental`, `W-IO I/O & analysis (FR-IO) | each producing subsystem | incremental`. Per CLAUDE.md working discipline this is explicit in the codex order breakdown (§7) and the final summary (STOP-RULE flag on "no active subsystem landed yet"). |

---

## 1. The coupling schedule (master `step()` orchestration)

### 1.1 The slot sequence (extends the CLAUDE.md invariant `f`-step; touches nothing inside it)

The CLAUDE.md core invariant step order — `collide → halo exchange → streaming
→ open-boundary BCs → boundary moments correction` — is the **hydrodynamic `f`
step** and is untouched by this spec (it is `Backend::run_span`,
`backend.rs:258-300`, already frozen by WVOF §4.1 and WSCAL §4.1). W-COUP adds
solver-orchestration-level slots **around** it, extending the pattern both
WVOF (`solver.rs` pre-pass before `run_span`, WVOF §4.2) and WSCAL (`solver.rs`
sub-step after `update_moments`, WSCAL §4.2) already use — this spec is the
**union schedule** those two specs each show one slice of:

```
PER SOLVER STEP (physical time dt, this step's value — see §1.2):

 0. [predictor half] SCALAR/REACTION PREDICTOR  (active-feedback only, C5)
    a. species ADE half-step (WSCAL slot, using u,phi from the END of the
       PREVIOUS step) — n_sub reaction sub-steps interleaved (C4) using
       whichever solver WREACT selects for the current stiffness (FR-COUP-02;
       reserved call, no default solver decided by W-COUP).
    b. equilibrium/pH solve to convergence at each reaction sub-step
       (ACIDBASE reserved slot; algebraic, not integrated in time — see §1.4).
    c. property predictor: rho(C), mu(C), sigma(C)[, T] from the half-stepped
       C (WSCAL active-feedback formulas, active-scalar-feedback.md §1-§4;
       NOT re-derived here — W-COUP only calls the property function).

 1. FORCE COMPOSITION (existing FORCE_COMPOSITION_SPEC pre-pass; T4 surface
    tension, T5 scalar/thermal buoyancy read the predictor's rho/mu/sigma).

 2. PHASE-FIELD SUB-STEP (WVOF slot, if W-VOF active): g collide -> halo ->
    stream -> phi = sum g_i  (WVOF_IMPL_SPEC §4.2 steps 1-2, verbatim, unchanged).

 3. HYDRODYNAMIC f STEP (Backend::run_span, UNCHANGED): collide -> halo ->
    stream -> open BCs -> update_moments. Produces rho, u = (ux,uy,uz)
    [F/2-corrected] for THIS step.

 4. [corrector] SCALAR/REACTION CORRECTOR  (active-feedback only, C5)
    a. species ADE corrector using the u,phi just produced by step 3
       (WSCAL slot, WSCAL_PASSIVE_SPEC §4.2 step 2, generalized to read the
       active velocity/phase state).
    b. reaction corrector, n_sub sub-steps (C4), consistent with 0.a/0.b.
    c. equilibrium/pH re-solve.
    d. property re-evaluation rho(C), mu(C), sigma(C)[, T] -> feeds NEXT
       step's force composition (step 1) and NEXT step's omega_field.
    e. IF strong-coupling iteration requested (stiff Da or stiff capillary,
       REQ FR-COUP-01 last clause): compare this corrector's property delta
       to the predictor's; if above a configured tolerance, re-run steps
       0.c-4.d once more (bounded iteration count, logged). This is the
       "optional flow-scalar iteration for stiff coupling."

 5. INTERFACIAL TRANSFER + BUBBLE/PBM UPDATE (WBUB reserved slot, if active):
    k_L a mass exchange between phases (INTERFACIAL_TRANSFER_SPEC reserved),
    then bubble/PBM number-density and diameter update (WBUB_PBM reserved) —
    both READ the corrector's converged C/phi/u and WRITE only their own
    state (bubble population, per-species interfacial flux term consumed as
    a source S^if by the NEXT step's species ADE, per REQ §3 S^if convention).

 6. OUTPUT/ANALYSIS REDUCTION (W-IO, this spec §2; cadence-gated, read-only,
    §1.6) + CONSERVATION DRIFT ACCUMULATION (§4; every step, O(1) cost).
```

**Ordering rationale (decided, physical):** slot 0 (predictor) must precede
force composition (slot 1) because the REQ sentence names `predictor →
property update → force composition → flow step` in that order — the flow
step needs `rho(φ,C)`/`μ(φ,C)`/`σ(φ,C)` *before* it collides. Slot 2
(phase-field) sits between the property predictor and the `f` step for the
same reason WVOF places it there (§4.2 of that spec: `φ` feeds `ρ,μ,F_s` into
collide). Slot 4 (corrector) runs after the `f` step because it needs the
*updated* `u,φ` (this is what makes it a corrector rather than a second
predictor — REQ's literal `scalar ADE corrector → reaction corrector →
property re-evaluation`). Slot 5 (interfacial transfer/PBM) runs last among
physics because it is REQ's `active(...) → property update → ...` composed
with the two-phase interfacial exchange, which by REQ §3's `S^if` convention
is consumed by "the next step's" ADE — it is correct for it to be the final
physics slot, writing state the *next* iteration's slot 0 reads. Slot 6 is
last unconditionally (pure reduction, no physics feedback, C6).

### 1.2 The dt-min rule (decision C3)

One scalar `dt` (in physical units; the lattice step itself is always
`Δt_lattice = 1`, unaffected — `dt` here is the **physical-time mapping** the
unit-conversion layer uses, REQ §2) is computed **once per run** from the
static configuration (not re-evaluated per step in phase 1 — a per-step
re-evaluation would require an inner iteration akin to adaptive time-stepping,
which REQ does not ask for and which First-version W-COUP defers):

```
dt = min(
       dt_capillary,     if W-VOF active:  Δt_sigma <= sqrt(rho_bar * dx^3 / (2*pi*sigma_max))
                          (REQ FR-COUP-01; sigma_max per active-scalar-feedback.md
                          §5.1 "when sigma is variable, use sigma_max")
       dt_marangoni,     if sigma feedback active: Δt_Ma <= C_Ma * mu_bar * dx
                          / (|grad_sigma|_max * W)     (active-scalar-feedback.md
                          §5.1; C_Ma is that spec's own derivation-required
                          constant — W-COUP consumes it, does not re-derive it)
       dt_buoyancy,      if scalar/thermal buoyancy active: Δt_b <= C_b *
                          sqrt(dx / |beta_C * Delta_C * g|)  (ibid.)
       dt_reaction,      if W-REACT active AND explicit reaction integration
                          selected: the explicit stability dt of whichever
                          scheme FR-COUP-02 selects (reserved; if
                          implicit/Rosenbrock-BDF is selected instead, this
                          term drops out because implicit steps are
                          unconditionally stable at the price of subcycling,
                          C4, not dt-limiting)
       dt_particle,      if two-way/four-way particles active: FR-PART-06's
                          Δt_p <~ T_col/10 (reserved; particles remain
                          one-way D-track today, no dt constraint from this
                          term at present)
       dt_hydrodynamic   the existing low-Mach unit-conversion feasibility
                          check (REQ §2: Ma_lattice <= 0.1, CFL/diffusion-number)
     )
```

`dt` is reported in `manifest.json` alongside which constraint was binding
(the `argmin`) — this is a **diagnostic**, not a tunable; a scenario whose
binding constraint is `dt_marangoni` at an unreasonably small value is a
signal to refine `W`/mesh or revisit `σ` variability, never a signal to
loosen the formula. **All terms whose subsystem is inactive drop out of the
`min` entirely** (not substituted with infinity via a magic sentinel — the
`min` is computed over a `Vec` built by pushing only active terms, so an
inactive subsystem cannot accidentally become the binding constraint through
a stale default).

### 1.3 Subcycling of fast processes (decision C4)

```
n_sub = ceil(dt_hydrodynamic_step / dt_fast_process)     n_sub >= 1
```

applied independently to **reaction** (`dt_fast_process = dt_reaction_intrinsic`,
the reaction system's own stiffness-limited or user-requested internal step,
FR-COUP-02) and **equilibrium/pH** (`dt_fast_process` is not a time step at
all — equilibrium is solved algebraically to convergence at *every* reaction
sub-step, so its "subcycle count" is really "solve once per reaction
sub-step," §1.4). During all `n_sub` reaction sub-steps, **transport state is
frozen** (`u`, `φ`, `ρ` from the current predictor/corrector snapshot do not
change) — this is the standard operator-splitting assumption (Strang-type)
that makes subcycling valid: reaction is locally an ODE in a frozen
advective/diffusive background over one hydrodynamic `dt`. The **splitting
error** this introduces is exactly the object REQ FR-COUP-01's "operator-
splitting error ... required" language flags, and it is measured (not
bounded a priori) by the dt-halving convergence test (§5, test C2): halving
`dt_hydrodynamic` while holding the reaction stiffness fixed must shrink the
corrector's property delta at the expected first-order (Strang: second-order
if symmetrized) rate.

### 1.4 Equilibrium (pH) as an algebraic constraint, not a transported field

The equilibrium solve (ACIDBASE_PH_SPEC, reserved) is **not** a rate process —
it is the statement "at every point and every reaction sub-step, the
acid-base species concentrations satisfy their equilibrium constants exactly"
(fast relative to transport and to the *other* finite-rate reactions in the
network). W-COUP's contract with that (unwritten) spec is only: **the
equilibrium solve runs after every reaction sub-step, using that sub-step's
just-updated total concentrations as the equilibrium system's conserved
totals** (charge balance + mass balance per acid-base pair), and its solved
pH/species-split feeds back into the *next* reaction sub-step's rate law
evaluation (since many rate laws are pH-dependent) within the same `n_sub`
loop. This is a slot contract (order + input/output), not a physics
definition — the equilibrium system itself (which pairs, which K values) is
100% ACIDBASE_PH_SPEC's scope.

### 1.5 The `active_scalar_lagged` relaxation flag (decision C5)

`active_scalar_lagged: bool` (default `false` = predictor-corrector fidelity
default). When `true`:

- Slot 0 (predictor) still runs (species ADE + reaction advance one full
  step using the **previous** step's `u,φ`).
- Slot 4 (corrector) is **skipped entirely** — no re-solve, no property
  re-evaluation against the new `u,φ`.
- Property fields (`ρ(C), μ(C), σ(C)[,T]`) feed the **current** step's force
  composition (slot 1) computed from the predictor only, and therefore are
  exactly one step stale relative to the `u,φ` they will next influence
  (explicit lag).

**Stability condition (required before this flag may be set `true` — a
config-validation warning, not a runtime clamp):** the lag is stable only
when the property-feedback timescale `τ_feedback` (the fastest of the
buoyancy/Marangoni/capillary response times, i.e. `1/max(dt_marangoni,
dt_buoyancy, dt_capillary)⁻¹`... concretely, the **explicit-feedback CFL**:
`dt_hydrodynamic / τ_feedback < 1`, where `τ_feedback` is estimated from the
same coefficients as §1.2's dt terms) is satisfied; this is the same class of
condition WALE's one-step SGS lag already documents (`les.rs:5-7`,
"one-step lag" precedent cited by WVOF §4.2) generalized to the active-scalar
lag. Config validation computes this ratio and warns (following the existing
`validate()` warning pattern, `lib.rs:622-690`) when it exceeds a threshold;
it does not reject, because the lag error is a *quantified*, benchmarked
approximation (measured by test C1, §5), not an undefined one.

**Lag-error benchmark (mandatory, REQ):** on the Marangoni or
concentration-dependent-viscosity bench already named by REQ FR-COUP-01
("on active-scalar bench (Marangoni or concentration-dependent viscosity),
feedback error converges under dt-halving"), run BOTH modes
(`active_scalar_lagged=false` reference vs `=true` relaxation) at a sequence
of halved `dt`; the lagged mode's error relative to the predictor-corrector
reference must **shrink monotonically** as `dt→0` (first-order lag error is
expected and acceptable; growing or non-monotonic error fails the test — see
§5 test C1).

### 1.6 Output cadence is independent of the physics cadence (decision C6)

Every `ReactorOutput` (§2) and every new `FieldKind` (§2) carries its own
`every` (steps) exactly like the existing `OutputSpec.every`
(`lib.rs:593-600`); the reduction pass (slot 6) only executes a given
output's computation on the steps where `step % every == 0` (or at run end).
This is unchanged from the existing `OutputSpec` cadence model — W-IO adds
new *kinds* of output, not a new cadence mechanism.

---

## 2. Output/analysis API — the QOI definitions

Every QOI below is: **formula, units, source field(s), definition metadata
(if any parameter/kernel choice exists), and its schema surface** (new
`FieldKind` per-cell variant, or new `ReactorOutput` whole-domain/probe
scalar). None of these are new physics — each formula is either an exact
algebraic reduction of already-computed conserved/derived quantities, or (for
`k_L a`, blend time, RTD) the same literature definition already named in
REQ §6/§2.

### 2.1 Per-species dissolved concentration `C_i` — `FieldKind::Concentration{ species }`

```
C_i(x,t) = conc field of species i's own ADE distribution (WSCAL h-set,
           WSCAL_PASSIVE_SPEC §3.2 "conc"; multi-component generalizes P10's
           single Vec<T> to Vec<ScalarField>, one per species — the
           WSCAL_MULTICOMPONENT_SPEC's scope, referenced not redefined here)
Units: physical concentration units per the scenario's unit-conversion layer
       (REQ §2) — mol/m^3 or kg/m^3 as configured; lattice-unit C is
       converted at the output boundary only (never inside the transport).
Source: SoaFields::conc[i] (per-species slot, WSCAL_MULTICOMPONENT_SPEC scope)
```

Schema: `FieldKind::Concentration { species: String }` — the per-cell field,
exported exactly like `Rho`/`Ux` today (`runner.rs` match arm addition).

### 2.2 Partial pressure `p_Xi` (general — pO2/pCO2/pH2/pN2 as instances) — `FieldKind::PartialPressure{ species }` + `ReactorOutput::PartialPressureBulk{ species }`

```
p_Xi(x,t) = C_i(x,t) / H_i          (Henry's law, dimensional form;
                                     H_i = Henry's constant for species i,
                                     [pressure]/[concentration], temperature-
                                     dependent if THERMAL_AXIS is active —
                                     H_i(T) supplied by INTERFACIAL_TRANSFER_SPEC,
                                     not redefined here)
Units: Pa (or bar/atm per scenario unit config)
Source: FieldKind::Concentration{species} (per-cell) -> per-cell division by
        H_i(T(x)) [T(x)=T_ref if thermal axis inactive]
Bulk (whole-liquid-phase) reduction:
  p_Xi_bulk = <C_i>_V,liquid / H_i(<T>_V)     (volume-average concentration
                                              over the liquid phase, weighted
                                              by phi if W-VOF active, else
                                              plain volume average)
```

`p_Xi` is a *general* definition parameterized by species — `pO2`, `pCO2`,
`pH2`, `pN2` are simply `species ∈ {O2, CO2, H2, N2}` instances of the SAME
`FieldKind::PartialPressure{species}` / `ReactorOutput::PartialPressureBulk`,
never separate variants (no per-gas special-casing, consistent with the
"general, all volatile species" requirement and the physics-discipline ban on
case-identity branches).

### 2.3 pH — `FieldKind::Ph` + `ReactorOutput::PhBulk`

```
pH(x,t) = -log10( a_H+(x,t) )       (a_H+ = activity of H+; phase 1 assumes
                                     ideal solution, a_H+ = C_H+ in mol/L, so
                                     pH = -log10(C_H+); non-ideal activity
                                     coefficients are an ACIDBASE_PH_SPEC
                                     extension, not redefined here)
Units: dimensionless (the standard pH scale; C_H+ itself in mol/L before the log)
Source: C_H+(x,t), the equilibrium-solved hydrogen-ion concentration
        (ACIDBASE_PH_SPEC's per-cell output, §1.4's algebraic equilibrium slot)
Bulk: pH_bulk = -log10( <C_H+>_V )   (volume-averaged BEFORE the log — a
      pH computed from an averaged-then-logged concentration is the
      physically meaningful "bulk pH"; averaging pH values directly would be
      averaging in log-space and is explicitly NOT the definition used here)
```

### 2.4 Temperature `T` — `FieldKind::Temperature` + `ReactorOutput::TemperatureBulk`

```
T(x,t) = the thermal-axis ADE distribution's macroscopic moment
         (THERMAL_AXIS_SPEC's own "conc"-equivalent for the temperature
         field, using the same D3Q7 ADE machinery as WSCAL with D -> alpha
         [thermal diffusivity] and Sc_t -> Pr_t; NOT redefined here —
         THERMAL_AXIS_SPEC owns the governing equation)
Units: K (or the scenario's configured temperature unit)
Source: SoaFields' thermal scalar slot (THERMAL_AXIS_SPEC scope)
Bulk: <T>_V, volume-averaged (phi-weighted over the liquid phase if W-VOF
      active, matching the p_Xi_bulk convention §2.2)
```

If the thermal axis is not active (REQ §1 lists it as "API-reserved
extension", not landed), `FieldKind::Temperature` / `ReactorOutput::
TemperatureBulk` are schema-valid but config validation rejects a scenario
that requests them without `thermal.enabled=true` — same rejection pattern
as any output requesting a field its producing subsystem does not have
active (§3.4).

### 2.5 Gas holdup `ε_g` (with definition metadata) — `FieldKind::GasHoldup{ definition }` + `ReactorOutput::GasHoldupBulk{ definition }`

Per REQ FR-IO-01 (load-bearing, verbatim carried forward — this spec does not
alter the definitions, only wires them into the output schema):

```
resolved-phasefield:
  eps_g_raw(x,t)             = <1 - phi>_V_filter        (local smeared void
                                                          fraction; V_filter =
                                                          the reduction/
                                                          averaging kernel,
                                                          metadata below)
  eps_g_thresholded(x,t; phi_c) = volume(phi < phi_c) / V_filter, default
                                   phi_c = 0.5    (both variants ALWAYS
                                   output together per REQ FR-IO-01: "both
                                   output")
point-bubble (WBUB, reserved):
  eps_g_bubble(x,t) = sum_b [ V_b * W_kernel(x - x_b) ] / V_filter
hybrid:
  eps_g_total = eps_g_resolved + eps_g_bubble  MINUS the double-counted
                overlap region (the resolved region's own bubble content is
                excluded from the point-bubble sum by construction — WBUB's
                switching criterion, FR-VOF-04, guarantees a bubble is EITHER
                resolved OR point, never both, so "minus double-count" is
                zero by the switching invariant, stated for completeness)
```

**Definition metadata (mandatory on every ε_g output — REQ FR-IO-01
"LOAD-BEARING")**: `filter_width` (`V_filter`'s linear scale, e.g. `max(3W,
6Δx)` if reusing the WVOF interface-band convention, or a user-set averaging
box), `averaging_volume` (whole-domain / a named sub-region / per-cell local
box), `time_window` (instantaneous / a stated moving-average window),
`kernel` (for point-bubble: the smearing kernel name, e.g. Gaussian/top-hat,
and its width). These four fields are **serialized alongside every ε_g
value**, never as a separate lookup — "any ε_g must be recomputable from a
snapshot" (REQ FR-IO-01).

Schema: `FieldKind::GasHoldup { definition: GasHoldupDefinition }` where
`GasHoldupDefinition` is `Raw | Thresholded{ phi_c: f64 } | Bubble |
Total`, each producing a per-cell field; `ReactorOutput::GasHoldupBulk{
definition, filter_width, averaging_volume, time_window }` for the
whole-domain scalar with metadata attached in the output record itself.

### 2.6 Sauter mean diameter `d_32` — `ReactorOutput::SauterMeanDiameter`

```
d_32 = sum_b(N_b * d_b^3) / sum_b(N_b * d_b^2)
```

Two source paths, matching REQ FR-VOF-04 / VR-STR-02:

```
resolved-phasefield (default): d_32 measured by INTERFACE SEGMENTATION —
   connected-component labeling of the phi < phi_c iso-region (phi_c=0.5,
   the SAME threshold as eps_g_thresholded, §2.5, for consistency), each
   component's volume V_k converted to an equivalent sphere diameter
   d_k = (6 V_k / pi)^(1/3), then d_32 = sum(d_k^3)/sum(d_k^2) over
   segmented components. (REQ §4.4 FR-VOF-04: "resolved-phasefield default
   measures d_32 by interface segmentation" — exact wording adopted.)
point-bubble + PBM (WBUB, reserved): d_32 from the PBM's own number-density
   distribution N(d) directly (no segmentation needed — PBM tracks a size
   distribution natively; Luo-Svendsen / Prince-Blanch kernels per REQ
   FR-VOF-04, owned by WBUB_PBM_IMPL_SPEC).
```

Segmentation is a pure post-processing reduction over `phi` (already-computed
field) — no new physics, a connected-component algorithm.

### 2.7 Per-species mass-transfer coefficient `k_L a` — `ReactorOutput::KlaPerSpecies{ species }`

Two formulas depending on interface mode (REQ §6 FR-IO-05 lineage / VR-STR-04):

```
Formula A (interface-flux integral, resolved-phasefield):
  k_L a = ( Integral_over_interface  S_if_i  dA )  /  ( V_liquid * (C_i,sat - C_i,bulk) )
  where S_if_i is the resolved normal Henry-partition flux (REQ §3
  "S_{k,liq}^{if} = -S_{k,gas}^{if}, flux positive into liquid",
  INTERFACIAL_TRANSFER_SPEC's own per-cell output at the interface) and
  C_i,sat = H_i * p_Xi_bulk,gas (the Henry-law saturation concentration, §2.2)

Formula B (correlation / point-bubble k_L a(C*-C) model):
  k_L a = the POINT-BUBBLE model's own coefficient (REQ §3
  "point-bubble = k_L a(C*-C)"), read directly from WBUB_PBM_IMPL_SPEC's
  per-bubble-class parameter (already a model input in that regime, not a
  post-hoc reduction) OR back-calculated from the bulk mass-balance rate:
  k_L a = ( d<C_i>_V,liquid/dt )  /  ( C_i,sat - <C_i>_V,liquid )
  (the "formula = interface integral or correlation, explicit" REQ VR-STR-04
  requirement — both paths must be reported with which one was used)
```

Units: 1/s (`k_L` has units of velocity [m/s], `a` is interfacial area per
liquid volume [1/m]; the product is what this QOI reports — never reported
as `k_L` or `a` separately without the product, since only the product is
uniquely defined without an independent area measurement in the resolved
path).

### 2.8 Reaction rates — `ReactorOutput::ReactionRate{ reaction_id }` + `FieldKind::ReactionRate{ reaction_id }`

```
r_j(x,t) = the WREACT rate law's own evaluated value for reaction j
           (mol/(m^3 s) or the configured rate units) — this is WREACT's own
           computed quantity (its rate law IS the physics; §1.4/§1.3's
           subcycle loop calls it every reaction sub-step). W-IO's
           contribution is exclusively EXPOSING the already-computed r_j(x,t)
           as an output field/reduction (spatial average <r_j>_V, or a
           per-cell field for visualization) — no new formula.
```

Schema: `FieldKind::ReactionRate { reaction_id: String }` (per-cell, for
visualization of hot spots), `ReactorOutput::ReactionRate { reaction_id }`
(volume-averaged scalar time series).

### 2.9 Conversion, yield, selectivity — `ReactorOutput::{Conversion,Yield,Selectivity}`

Defined from stoichiometry + inlet/outlet concentrations (standard reaction-
engineering definitions — an exact algebraic reduction, not a new closure).
For a reaction network with stoichiometric matrix `ν` (already WREACT's own
input; REQ FR-COUP-02 references element conservation against the same `ν`)
and a designated limiting reactant `A` and desired product `P`:

```
Conversion   X_A  = (C_A,in - C_A,out) / C_A,in            (dimensionless, [0,1])
Yield        Y_P  = ( (C_P,out - C_P,in) / |nu_P| )
                     / ( (C_A,in - C_A,out) / |nu_A| )      (moles P formed per
                                                             mole A consumed,
                                                             stoichiometry-
                                                             normalized;
                                                             dimensionless, [0,1]
                                                             for a single-pass
                                                             product-forming step)
Selectivity  S_P  = Y_P / X_A                               (= yield per unit
                                                             conversion;
                                                             dimensionless,
                                                             [0,1] by
                                                             construction when
                                                             Y_P <= X_A, which
                                                             holds whenever P
                                                             is not also
                                                             consumed by a
                                                             side reaction that
                                                             REGENERATES A —
                                                             an edge case the
                                                             element-balance
                                                             diagnostic, §4,
                                                             flags rather than
                                                             the definition
                                                             clamping)
```

`C_*,in`/`C_*,out` are the inlet-face and outlet-face flow-weighted average
concentrations (`Integral(C*u.n dA) / Integral(u.n dA)` over the respective
open face) — the same face-integration convention REQ FR-ROT-04 already uses
for `N_Q` ("integration surface, velocity components ... defined"), reused
here rather than re-invented. **Bounded [0,1] is a consequence of the
definition for a well-posed single-limiting-reactant network, not an
enforced clamp** — if a run produces `X_A` or `Y_P` outside `[0,1]`, that is
a diagnostic signal (mass-balance violation, wrong limiting-reactant choice,
or a genuinely multi-pass/recycle topology this single-pass formula does not
model), reported as-is (§4.5 discipline), never clamped into range.

### 2.10 Blend time — `ReactorOutput::BlendTime`

```
Definition (REQ FR-IO-05, adopted verbatim): the elapsed time from a tracer
injection event until the domain (or a named detection region/probe set)'s
coefficient of variation of the tracer concentration falls below a
configured threshold:
  CoV(t) = std_dev( C_tracer(x,t) over detection region ) / mean( C_tracer(x,t) )
  t_blend = min{ t : CoV(t) < CoV_threshold  for all t' >= t within the
                 window }     (the "for all t' >= t" clause guards against a
                 transient dip below threshold that is not yet true mixedness
                 — a single instantaneous crossing is not sufficient)
Units: s (physical time, via the unit-conversion layer)
Source: FieldKind::Concentration{tracer} sampled over the injection/detection
        surfaces the scenario declares (REQ FR-IO-05: "injection/detection
        surfaces explicit per scenario" — new schema fields, §3.5)
CoV_threshold: scenario-declared (no engine default baked in — REQ says
        "Thresholds ... explicit per scenario", so this is a REQUIRED
        scenario field, not a silent default; omitting it is a config
        validation error per the "no silent fallback for a physical
        parameter" ban-list rule)
```

### 2.11 Residence-time distribution (RTD) — `ReactorOutput::Rtd`

```
Definition (standard tracer-response RTD, e.g. Levenspiel): impulse or step
tracer injection at the inlet; E(t) = the normalized exit-age distribution
measured from the OUTLET flow-weighted tracer concentration time series
C_out(t):
  Impulse-response form:  E(t) = C_out(t) / Integral_0^inf C_out(t') dt'
  Step-response form:     F(t) = C_out(t) / C_out(infinity),  E(t) = dF/dt
Mean residence time: tau_bar = Integral_0^inf t*E(t) dt
Variance:            sigma_t^2 = Integral_0^inf (t-tau_bar)^2 * E(t) dt
Normalization check (behavior anchor, NOT a band): Integral_0^inf E(t) dt = 1
        (exactly, up to the truncation error of the finite run window and the
        numerical-integration quadrature — see §5, test C4)
```

Units: `E(t)` has units of inverse time (1/s); `τ̄`, `σ_t` in s.
Source: `FieldKind::Concentration{tracer}` flow-weighted at the outlet face
(the same face-integration convention as §2.9's inlet/outlet averages), time-
series-accumulated exactly like an existing `ProbeSpec::Point` but face-
integrated instead of point-sampled — a new `ProbeSpec::FaceFlux{ face,
species, every }` probe kind (§3.5) feeding the RTD reduction.

---

## 3. Scenario schema extension

### 3.1 Species registry (top-level)

```jsonc
"species": [
  { "name": "O2", "diffusivity": 2.1e-9, "sct": 0.7,
    "henry": { "H": <Pa per (mol/m^3)>, "H_temperature_dependence": null },
    "initial": { "type": "uniform", "value": 0.0 } }
  // WSCAL_MULTICOMPONENT_SPEC owns per-species transport fidelity; W-COUP/
  // W-IO only requires the registry entries needed to LABEL and UNIT-CONVERT
  // each species' C_i, p_Xi, and reaction/conversion bookkeeping (§2.1-2.2,
  // §2.9). No transport physics field is redefined here.
]
```

### 3.2 Reaction network (top-level, WREACT scope; registry shape only)

```jsonc
"reactions": [
  { "id": "r1", "stoichiometry": { "A": -1, "P": 1 },   // nu, signed
    "kinetics": { "type": "power_law", "k": <rate const>, "order": { "A": 1 } },
    "stiffness_hint": "auto" }                            // FR-COUP-02 solver
                                                          // switch input;
                                                          // WREACT-owned
]
```

`stoichiometry` is the `ν` matrix row this spec's §2.9 conversion/yield/
selectivity/element-balance (§4.3) formulas consume — **the only field of
this block W-COUP/W-IO reads**; the kinetics block's contents are entirely
WREACT_IMPL_SPEC's scope, registered here only so the schema has one place
species+reactions live (REQ §5 names them together in the coupling
dataflow).

### 3.3 Equilibrium system (top-level, ACIDBASE_PH_SPEC scope; registry shape only)

```jsonc
"equilibrium": {
  "acid_base_pairs": [
    { "acid": "HA", "base": "A-", "pKa": <value> }
  ],
  "charge_balance_species": ["H+", "OH-", "A-", "Na+"]
}
```

Read by §1.4's slot contract only for "which totals are conserved across the
solve" — the solver algorithm is ACIDBASE_PH_SPEC's scope.

### 3.4 Thermal block (top-level, THERMAL_AXIS_SPEC scope; registry + enable flag)

```jsonc
"thermal": {
  "enabled": false,                    // default false = axis inactive
  "alpha": <thermal diffusivity>,
  "pr_t": 0.85,
  "initial": { "type": "uniform", "value": 293.15 },
  "feedback": { "sigma_T": null, "beta_T": null, "viscosity": null }
}
```

Config validation: any scenario output requesting `FieldKind::Temperature`
or `ReactorOutput::TemperatureBulk` with `thermal.enabled=false` is a
**hard config-validation error** (not a warning) — "requested output has no
producing subsystem," the same class of rejection FR-VOF-03 already uses for
raw-φ (C8).

### 3.5 Inlet phase + species/mixing probes (extends existing `edges`/`probes`)

```jsonc
"edges": [
  { "name": "inlet1", "kind": "velocity", "inlet_phase": "liquid",   // FR-VOF-03,
                                                                      // unchanged
    "species_inlet": { "O2": 0.0, "tracer": 1.0 },                   // NEW: per-
                                                                      // species C
                                                                      // at this
                                                                      // face; NOT
                                                                      // a phi value
    "temperature_inlet": 293.15 }                                   // NEW, only
                                                                      // valid if
                                                                      // thermal.enabled
]
"probes": [
  { "type": "faceFlux", "face": "outlet1", "species": "tracer", "every": 10 }
  // NEW ProbeSpec variant feeding RTD (§2.11) and inlet/outlet conversion
  // averages (§2.9) — a face-integrated flow-weighted C(t), not a point sample.
]
"mixing": {                                                          // NEW top-
                                                                      // level block
  "blend_time": { "tracer": "tracer", "injection_face": "inlet1",
                   "detection_region": "domain", "cov_threshold": 0.05 },
  "rtd": { "tracer": "tracer", "injection": "impulse",
           "injection_face": "inlet1", "outlet_probe": "faceFlux@outlet1" }
}
```

### 3.6 Config validation rules (new, additive to the existing `validate()` pass)

- **Raw-φ rejection extended (C8)**: `species_inlet`, `temperature_inlet`
  never accept or imply a `phi` value; the existing FR-VOF-03 validator
  (`lib.rs`, the `inlet_phase` check) is untouched — these new fields simply
  cannot express φ (they are keyed by species name / are a plain
  temperature), so the rejection is structural (no φ-shaped field exists to
  reject) rather than a new runtime check. A schema-level test (§5, test W3)
  asserts no new field of this spec's schema additions accepts a bare `phi`
  key.
- **Output-without-producer rejection (§3.4 pattern, generalized)**: any
  `OutputSpec`/`ReactorOutput` naming a species not in the `species` registry,
  a `reaction_id` not in `reactions`, or a thermal output with
  `thermal.enabled=false`, is a hard config-validation error, not a silent
  no-op and not a runtime panic.
- **`active_scalar_lagged` stability warning (§1.5)**: computed and warned
  per the existing `warn()` pattern, never rejected.
- **`mixing.blend_time.cov_threshold` required, no default** (§2.10, "no
  silent fallback for a physical parameter").
- **Disabled-subsystem bit-identity (B-6, C2)**: this is a *test* obligation
  (§5, test W1), not a schema validation rule — noted here for completeness
  of the schema/behavior contract.

---

## 4. Conservation drift diagnostics

Every quantity below is accumulated **every step** (O(1) reduction cost,
piggybacked on the existing global-reduction pass already used for
`probe_state_hash` and the mass/momentum bookkeeping WVOF's V4/V7 and WSCAL's
V4 already perform per-subsystem) and reported as a **time series** in
`manifest.json`'s `reactor_outputs.drift` object plus a final summary
(`max` and `final` drift). None of these values are ever fed back into the
simulation (C9) — accumulate, compare to threshold, report; never rescale.

| Quantity | Formula | Expected behavior | Threshold class |
|---|---|---|---|
| **Mass** | `Σ_cell ρ(x) dV` (single-phase) or `Σ_cell [φρ_l+(1-φ)ρ_g] dV` (two-phase, W-VOF active) | conserved except through open boundary mass flux (accounted separately as an explicit boundary-flux term, not lumped into "drift") | VR-STR-05 mass drift band |
| **Momentum** | `Σ_cell ρ(x) u(x) dV` | conserved except through boundary flux + net body force (gravity, IBM reaction) — the drift-minus-accounted-forces residual is the diagnostic, not the raw momentum change | VR-STR-05 momentum drift band |
| **Per-species scalar totals** | `Σ_cell C_i(x) dV` (or phase-wise `Σ_cell α_q C_{i,q} dV`, REQ §3 conservative form) | conserved except boundary flux + net reaction source `Σ_cell R_i(C) dV` (accounted, subtracted before calling the residual "drift") | VR-STR-05 scalar-total drift band (WSCAL V4 precedent, generalized per-species) |
| **Gas-phase volume** | `Σ_cell (1-φ) dV` (resolved) or `Σ_b V_b` (point-bubble) or the hybrid sum (§2.5's `eps_g_total`, integrated over volume) | conserved except sparger inflow + degassing outflow (accounted) — the FR-VOF-03 gas-volume-balance test (WVOF V6) generalized to a continuous drift monitor | VR-STR-05 gas-volume drift band |
| **Element balance** | for each atomic element `E` conserved by the reaction network: `Σ_cell Σ_i (n_{E,i} C_i(x)) dV`, where `n_{E,i}` = atoms of `E` per mole of species `i` (from the species registry's molecular formula — a REQUIRED registry field for any species participating in `reactions`, not inferred) | conserved except boundary flux (accounted) — reactions by construction cannot change element totals (`Σ_i n_{E,i} ν_{i,j} = 0` for every reaction `j`, an algebraic identity of a balanced stoichiometry, checked at config-validation time, not at runtime) | new: VR-STR-05 extension, element-balance closure (§5, test C3) |
| **Particle / bubble count** | `N_particles(t)` / `N_bubbles(t)` (population sizes) | conserved except explicit injection/removal/breakup-coalescence events (each logged and accounted so the residual isolates unaccounted loss/gain — e.g. a particle silently vanishing off-grid) | VR-STR-05 particle-count drift band (D-track precedent, generalized) |
| **Energy-like quantities** | kinetic energy `Σ_cell ½ρ|u|²dV`, interfacial free energy `Σ_cell [β φ²(1-φ)² + κ|∇φ|²/2] dV` (W-VOF), particle kinetic energy | **monitoring only, NOT exactly conserved** (REQ §8 VR-STR-05 explicit: "Energy-like = monitoring only") — no threshold failure on drift alone; reported as a time series for behavior-validity review (dissipation should be monotone-decreasing in an unforced decay, growing energy without an active energy source is the diagnostic signal) | monitoring, no pass/fail band |

**Accounted vs. unaccounted drift.** Every non-energy row above computes
"drift" as the **residual after subtracting every physically-expected
change already accounted for by a logged mechanism** (boundary flux
integral, reaction source integral, sparger injection volume, particle
injection/removal event log). This is what makes the diagnostic
non-trivial: a monitor that just reports "total changed" would fire on every
run with an inlet/outlet; the monitor that matters is "total changed by
MORE than the sum of the mechanisms we already know about." Each accounted
mechanism is itself computed from already-existing reduction machinery
(face-integrated flux, §2.9's inlet/outlet convention; reaction rate
integral, §2.8) — no new physics, purely bookkeeping arithmetic.

### 4.5 Diagnostics are never clamps (discipline, restates C9)

No drift monitor, no QOI (conversion/yield/selectivity bounded-by-
construction, §2.9; ε_g; d_32) is ever used to rescale, clamp, or correct a
transported field. An out-of-band drift or an out-of-[0,1] QOI is reported
exactly as measured and routed per the physics-discipline escalation table
(a core-engine defect if the drift traces to a code path, a spec defect if
the definition itself is inapplicable to the scenario topology, §2.9's
recycle-topology caveat).

---

## 5. Validation → T17 VR-STR-05 / VR-STR-07

Tests are **authored adversarially by codex/Opus from this spec**, in a test
worktree that never shares with the implementation worktree (CLAUDE.md team
convention; REQ §8). Because this tier lands last (C10), every test below is
written against **whichever subsystems are actually landed at test-authorship
time** — the table's "Grid/backend" column notes the minimal active set each
test needs; a test that needs an unlanded subsystem is authored now (compiles
red against a stub) and goes green when that subsystem lands (the same
"compiles red, goes green incrementally" pattern WVOF §7/WSCAL §5 already
use for their own O1→O2 sequencing).

| ID | Test | Metric & band | Grid / steps / backend / active-set needed | T17 row |
|---|---|---|---|---|
| **C1** | **Active-scalar lag-error dt-halving convergence** (mandatory per REQ §5 FR-COUP-01). Marangoni or concentration-dependent-viscosity bench (whichever active-scalar-feedback.md feedback lands first), run at `active_scalar_lagged=false` (reference) and `=true` (relaxation) across a halved-dt sequence. | Lagged-mode error vs the `=false` reference **shrinks monotonically** as `dt→0` (first-order lag error expected; growing/non-monotonic FAILS). | active set: W-VOF + at least one active-scalar feedback (σ or μ). Until then: compiles red against a stub scenario. | VR-STR-05 (REQ §5 mandatory dt-halving) |
| **C2** | **Reaction-subcycle operator-splitting convergence.** A stiff reaction (`Da >> 1`) coupled to slow transport (`Pe` moderate); run at halved `n_sub` (equivalently halved `dt_hydrodynamic` at fixed reaction stiffness). | Corrector property delta (or a probed species concentration) converges at the expected splitting order (first-order Strang, or the order WREACT's actual splitting delivers — measured, then frozen) as `dt→0`. | active set: W-SCAL + W-REACT (reserved). Compiles red until W-REACT lands. | VR-STR-05 |
| **C3** | **Element-balance closure.** A reaction network with 2+ elements tracked (e.g. a simple `A -> P` with a labeled atom), run to steady/quasi-steady state with known inlet/outlet/reaction-consistent stoichiometry. | `\|Σ_element(t) − Σ_element(0) − accounted_boundary_flux(t)\| / Σ_element(0) < 1e-6` (f64) / `< 1e-4` (f32) — round-off-scale, not a band (reaction cannot change element totals by construction, §4). | active set: W-SCAL + W-REACT (reserved). Compiles red until W-REACT lands. | VR-STR-05 (new element-balance row) |
| **C4** | **RTD normalization + conversion/yield/selectivity behavior anchors.** Full mini-reactor smoke (see C6 below) impulse tracer at inlet, exit-age `E(t)` measured. | (a) `\|∫E(t)dt − 1\| < 1%` (finite-window truncation tolerance, tightened as the run window grows); (b) **conversion increases monotonically with residence time** (behavior anchor — sweep mean residence time via flow rate, assert `X_A` is monotone increasing, `.claude/skills/lbmflow-physics-discipline` Rule 3 "Monotone trend" template); (c) `Y_P ∈ [0,1]`, `S_P ∈ [0,1]` for the well-posed single-limiting-reactant network (§2.9). | active set: full mini-reactor (C6). Compiles red until W-VOF+W-SCAL+W-REACT land. | VR-STR-05, behavior anchors |
| **C5** | **`probe_state_hash` coupled regression, single-backend.** With a fixed active subsystem set and fixed scenario, `probe_state_hash` is bit-identical across repeated runs on the SAME backend (FR-COUP-04: "single-backend only"). | exact bit equality, single backend (CpuScalar) | whatever active set exists at test-authorship time; runs incrementally as each subsystem lands. | VR-STR-05 (FR-COUP-04) |
| **W1** | **Disabled-subsystem bit-identity (B-6, decision C2).** A scenario with EVERY W-COUP-orchestrated subsystem `None`/disabled produces a `probe_state_hash` bit-identical to the pre-W-COUP engine on the same scenario (the same class of test as WSCAL's V5 and WVOF's per-order DoD, generalized to the ORCHESTRATION layer itself — proves W-COUP's slot insertion adds no-op branches only, never a hidden cost or hidden state mutation when nothing is active). | exact `probe_state_hash` equality | cavity + cylinder presets, CpuScalar | VR-STR-05 (B-6 invariance) |
| **VR-STR-07-1** | **Initialization independence — spin-up time.** Same scenario run with varying impeller/gas ramp duration (FR-INIT-01) before the statistics-sampling window starts. | Quasi-steady statistics (e.g. mean `<u>`, `Np`, blend time) computed over the post-spin-up window agree within a frozen threshold across ramp durations. | active set: whatever produces a quasi-steady stat (minimally W-ROT, already LANDED — this test can be authored and run TODAY against the landed W-ROT/stirring path, independent of W-COUP/W-IO landing). | VR-STR-07 |
| **VR-STR-07-2** | **Initialization independence — statistics-start time.** Same scenario, varying WHEN the statistics accumulator starts sampling (holding total run length fixed). | Quasi-steady statistics agree within the same frozen threshold as VR-STR-07-1 regardless of stats-start time, once past the physically-required spin-up. | Same active-set note as VR-STR-07-1. | VR-STR-07 |
| **C6** | **Full mini-reactor integration smoke** (sparger + species + reaction + pH + transfer) — exercises the WHOLE loop (§1.1's slots 0-6) and reports every §2 output. Not a tight band test; a smoke gate that the loop runs to completion, every output is finite/well-formed, and every drift monitor (§4) stays within its band. | pass/fail = (a) run completes without divergence; (b) every §2 output (C_i, p_Xi, pH, ε_g, d_32, k_L a, r_j, X_A/Y_P/S_P, blend time, RTD) is produced, finite, and carries its required metadata (ε_g's filter_width/averaging_volume/time_window, §2.5); (c) every §4 drift monitor within its band; (d) the behavior anchors of C4 hold. | active set: ALL of {W-VOF, W-SCAL, W-REACT, ACIDBASE, INTERFACIAL_TRANSFER, W-BUB} — the terminal integration gate for this entire tier; compiles red until every dependency lands, goes green last. | VR-STR-05 (coupled regression, the REQ-named "mini-reactor integration smoke") |

**Mandatory negative / consistency tests (REQ §8, this tier's share):**

- **Raw-φ rejection extended (C8, §3.6)**: a scenario JSON that attempts to
  smuggle a `phi` value through `species_inlet` or `temperature_inlet` (e.g.
  a malformed key) must be REJECTED by config validation — since the schema
  structurally cannot express φ in these fields, this test asserts the
  *type* of `species_inlet`/`temperature_inlet` never deserializes a
  `{"phi": ...}` shape (a schema-level negative test, not a runtime check).
- **Disabled-path bit-identity (W1)**: already listed above; restated here
  because REQ frames it as part of the mandatory negative-test set for any
  B-6-governed feature.
- **dt-min argmin sanity**: a scenario where the reaction dt term is
  artificially made binding (very stiff kinetics, explicit integration
  forced) must show `dt`'s reported binding constraint as `dt_reaction`,
  not silently falling back to `dt_hydrodynamic` — proves the `min()` is
  computed over the live `Vec` of active terms (C3), not a stale default.
- **Drift-vs-accounted-mechanism sign test**: disabling the boundary-flux
  accounting (a test-only mutant) on an open-boundary scenario must make the
  mass/scalar drift monitor FAIL its band (proves the "accounted" subtraction
  in §4 is load-bearing, not decorative).

### 5.6 Behavior-validity review (mandatory, REQ / CLAUDE.md)

After every validation run, before reporting: review the *observed* pattern,
not just the gated metric. Specifically for this tier: (a) conversion vs.
residence-time curve has the expected concave-toward-saturation shape for a
first-order-ish reaction (not a step function or a non-monotone wiggle
outside noise — a wiggle would indicate an RTD/probe timing bug, not a real
effect); (b) RTD `E(t)` is a single-humped, non-negative curve for a
well-mixed or plug-flow-like reactor topology (a negative excursion or a
long-time non-decay signals a leaking/reflecting outlet BC, not real
long-tail mixing, unless the geometry genuinely has dead zones — check the
flow field before concluding "dead zone"); (c) drift-monitor time series are
flat (not slowly ramping) once the accounted-mechanism subtraction (§4) is
applied — a slow ramp in the RESIDUAL (post-subtraction) means an unaccounted
leak, not a real physical effect; (d) the lagged-vs-corrector comparison
(test C1) shows the lag error concentrated where `∇σ`/`∇μ` is largest (near
the interface/mixing front), not spread uniformly (uniform lag error would
suggest a bug unrelated to the feedback term itself, e.g. a dt-mapping
error). Record the review in PHYSICS.md or the track's findings file. A
metric passing its band does **not** validate a pattern no band covers.

---

## 6. Stability — the composite dt

Restated from §1.2/§1.3 as the single stability contract this tier owns:

```
dt_physical = min( dt_capillary, dt_marangoni, dt_buoyancy, dt_reaction_explicit,
                    dt_particle, dt_hydrodynamic )   [only over ACTIVE terms]

n_sub_reaction    = ceil( dt_physical / dt_reaction_intrinsic )
n_sub_equilibrium = solved algebraically every reaction sub-step (not a
                    time-stepped quantity — no CFL-style dt of its own)

active_scalar_lagged stability (config-validation warning threshold):
  dt_physical / tau_feedback < 1,   tau_feedback estimated from the same
  coefficients as dt_marangoni / dt_buoyancy / dt_capillary (§1.5)
```

No term in this composite is invented by this spec — each is either the
literal REQ FR-COUP-01 formula, the active-scalar-feedback.md §5.1 formula
(both cited, not re-derived), or the existing low-Mach feasibility check
(REQ §2). W-COUP's sole addition is the **`min()` over the live active set**
and the **subcycle-count formula** — pure scheduling arithmetic.

---

## 7. Codex order breakdown

Four orders, file-conflict-aware. **This tier lands LAST** (decision C10):
none of these orders should be dispatched until at least one producing
subsystem (W-SCAL passive and/or W-VOF phase 1) has landed enough that the
orchestration loop and the reduction layer have real state to operate on.
One order = one bundle = one dedicated worktree (CLAUDE.md team convention).
Implementation and adversarial-test orders never share a worktree.

| Order | Scope | Primary files (conflict boundary) | Gate | Depends on |
|---|---|---|---|---|
| **O1 — Coupling orchestration** | The slot sequence (§1.1) as a new solver-level orchestration method (`Solver::run_coupled_step`, alongside the existing `run_staged_step`); dt-min computation (§1.2) + reporting into `manifest.json`; reaction/equilibrium subcycle loop driver (§1.3-1.4) — calling INTO whichever subsystem sub-step functions exist, with every not-yet-landed slot compiled as a `None`-gated no-op (C2); `active_scalar_lagged` flag + stability-warning validation (§1.5). | `crates/lbm-core/src/solver.rs` (new orchestration method; does NOT touch `run_staged_step`/`run_span` bodies, only adds a caller above them), `crates/lbm-scenario/src/lib.rs` (`active_scalar_lagged` field + its validation warning) | W1 (disabled-subsystem bit-identity) green immediately (trivially, since every slot is `None`); C1/C2 compile red against stubs, go green as W-VOF/W-SCAL/W-REACT land. | At least the `Option`-gated state shape of ONE producing subsystem (W-SCAL and/or W-VOF) landed, so there is a real (non-stub) slot to call. |
| **O2 — Output-API / QOI reductions** | The full §2 reduction layer as a new module (`crates/lbm-core/src/reactor_outputs.rs`): per-species `C_i`/`p_Xi`/pH/T per-cell reductions (§2.1-2.4), `ε_g` with metadata (§2.5), `d_32` segmentation (§2.6), `k_L a` both formulas (§2.7), reaction-rate exposure (§2.8), conversion/yield/selectivity face-integration (§2.9), blend time (§2.10), RTD (§2.11). Pure read-only functions over `SoaFields`/subsystem state — no solver-state mutation. | `crates/lbm-core/src/reactor_outputs.rs` (new file — no conflict with O1's `solver.rs` edit beyond a call-site addition at slot 6, §1.1) | Each reduction unit-tested against a hand-computed value on a synthetic field (e.g. a known-analytic `phi` blob for `d_32` segmentation, a known-analytic Henry constant for `p_Xi`); full integration deferred to O4/C6. | Same producing-subsystem prerequisite as O1; disjoint file from O1, so **parallelizable with O1** once that prerequisite holds. |
| **O3 — Scenario schema + CLI outputs** | §3's schema additions (`species`, `reactions` registry shape, `equilibrium` registry shape, `thermal` block, `inlet_phase`-adjacent `species_inlet`/`temperature_inlet`, `mixing` block, `ProbeSpec::FaceFlux`); config validation rules (§3.6); new `FieldKind` variants (`Concentration`, `PartialPressure`, `Ph`, `Temperature`, `GasHoldup`, `ReactionRate`) wired into the existing exhaustive-match consumers (`crates/lbm-cli/src/runner.rs` PNG/CSV/VTK export arms, mirroring how `ShearRate`/`QCriterion` were added); `manifest.json`'s new `reactor_outputs` object (§4's drift time series + §2's whole-domain `ReactorOutput` scalars). | `crates/lbm-scenario/src/lib.rs` (schema + validation — same file as O1's `active_scalar_lagged` field; both are additive struct-field/enum-variant edits, mechanical merge like the WVOF/WSCAL `fields.rs` precedent, WSCAL §7.1), `crates/lbm-cli/src/runner.rs` (new `FieldKind` match arms + `Manifest.reactor_outputs`) | New `FieldKind` variants exhaustive-match-compile against O2's reduction functions; `manifest.json` schema test asserts every §2/§4 QOI round-trips through serde. | Depends on O2 (needs the reduction functions to call); schema-only edits can start in parallel with O2's implementation and integrate at the end. |
| **O4 — Adversarial integration tests (codex adversarial, separate worktree)** | All of §5 (C1-C6, W1, VR-STR-07-1/2) + the mandatory negative/consistency tests (raw-φ-rejection extension, disabled-path bit-identity, dt-min argmin sanity, drift-accounting sign test). Authored from THIS spec, not from the impl. VR-STR-07-1/2 can be authored and start running IMMEDIATELY (they need only the already-landed W-ROT stirring path, no dependency on O1-O3). | `crates/lbm-core/tests/wcoup_*.rs`, `crates/lbm-scenario/tests/*` (new files only — no impl-file conflict) | Tests compile red against stubs, go green against O1-O3 as they land; freeze bands in VALIDATION.md T17 VR-STR-05/07. VR-STR-07-1/2 run standalone today. | None for VR-STR-07-1/2 (landed W-ROT suffices); O1-O3 for C1-C6/W1. |

**Critical-path ordering within this tier:** {O1 ∥ O2} → O3 → C6 (the
terminal mini-reactor smoke, which needs every schema/output wired). O4 runs
concurrently from the start (test worktree), with its VR-STR-07 rows
executable today independent of this tier's own landing.

**Per-order DoD (all orders):** existing tests green *without modification*;
every orchestration slot bit-identical to today when its subsystem is
`None`/disabled (W1, B-6 invariance); every new `FieldKind`/`ReactorOutput`
serializes/deserializes round-trip; behavior-validity review (§5.6) recorded
for every validation run that produces one.

### 7.4 The PHYSICS.md entry (mandatory, but content-light — this tier adds no physics)

Because this spec introduces **zero new physical models** (pure
orchestration + pure reduction, per the discipline note in the preamble),
the PHYSICS.md entry O1 lands is a **decision-record entry**, not a
model-stack entry:

> **W-COUP coupling orchestration + W-IO output/analysis layer — pure
> scheduling and reduction, no new physics.** The predictor-corrector active-
> feedback loop, dt-min rule, and reaction/equilibrium subcycling implement
> REQ §5 FR-COUP-01/02 literally (no new formula); the QOI reductions (§2)
> are exact algebraic definitions (Henry's law, stoichiometric conversion/
> yield/selectivity, RTD moments) or literature definitions already cited by
> the producing subsystem's own spec (WVOF's ε_g/d_32, INTERFACIAL_TRANSFER's
> k_L a). **Why here (not derivable from code)**: this tier lands LAST in the
> DAG (REQ §11) because every slot/QOI is a call-out to a producing
> subsystem; record here which subsystems were active when each T17
> VR-STR-05/07 band was frozen, since a band frozen with only {W-SCAL, W-VOF}
> active does not automatically extend to the full mini-reactor set (C6) —
> re-characterize when W-REACT/ACIDBASE/W-BUB land.

---

## 8. Coexistence / full dependency ordering across all tiers

This is the **top of the DAG** (REQ §11: `W-COUP | active set | incremental`,
`W-IO | each producing subsystem | incremental`). The full ordering, restated
from REQ §11 with this spec's internal order (§7) inserted at the tip:

```
W0 (LANDED: D3Q19/27, cumulant, Guo)
 |
 +--> W-GRAV (LANDED) --> W-VOF (PENDING, critical path) --+
 |                          |                              |
 +--> W-SCAL passive (SPEC'd, WSCAL_PASSIVE_SPEC) ----------+--> W-SCAL active
 |                          |                              |    feedback
 |                          v                              |    (active-scalar-
 +--> W-STRESS (PARTIAL) -> W-LES (LANDED) -> W-PART (LANDED  |    feedback.md)
                                               D-track P2)     |
                          W-VOF --> W-BCTOP (PENDING)          |
                          W-VOF, W-SCAL --> W-REACT             |
                          (WREACT_IMPL_SPEC, not yet written)  |
                          W-SCAL --> ACIDBASE (pH)               |
                          (ACIDBASE_PH_SPEC, not yet written)  |
                          W-VOF, W-SCAL, W-EXT --> W-BUB (PBM)   |
                          (WBUB_PBM_IMPL_SPEC, not yet written) |
                          W-VOF --> INTERFACIAL_TRANSFER          |
                          (INTERFACIAL_TRANSFER_SPEC, not yet    |
                          written)                                |
                          W-SCAL --> THERMAL_AXIS (optional)       |
                          (THERMAL_AXIS_SPEC, not yet written)    |
                                                                    v
                          ============================================
                          THIS SPEC: W-COUP (orchestrates ALL of the
                          above's sub-steps) + W-IO (reduces ALL of the
                          above's fields into REQ §6 QOIs)
                          ============================================
                                     |
                                     v
                          T17 VR-STR-05 (coupled regression + drift)
                          T17 VR-STR-07 (init independence)
                          C6 mini-reactor integration smoke (terminal gate)
```

**What this spec's own §7 orders (O1-O4) can start on TODAY**: O1's `None`-
gated orchestration skeleton and O4's VR-STR-07-1/2 rows need only the
LANDED set (W0, W-GRAV, W-ROT, W-LES, W-PART D-track). **What waits on a
sibling spec landing**: every C1-C6/W1 test that names a specific subsystem
(W-VOF, W-SCAL active, W-REACT, ACIDBASE, W-BUB, INTERFACIAL_TRANSFER) is
authored now, compiles red, and goes green incrementally — this is the
explicit "incremental across producing subsystems" delivery mode REQ §11
names for both W-COUP and W-IO, not a reason to delay authorship.

**Coexistence with each sibling spec's own O-numbering**: this tier's O1
(`solver.rs` orchestration caller) sits ABOVE every sibling spec's own
solver-level sub-step (WVOF's phase-field pre-pass, WSCAL's scalar sub-step)
— it calls them, it does not modify their internals, so there is no file-
region overlap beyond the same additive-edit-to-a-shared-file pattern
WVOF §7.1 already documents and resolves mechanically (both-add, no logic
merge). This spec's O3 schema edits interleave with WSCAL_MULTICOMPONENT/
WREACT/ACIDBASE/THERMAL's own eventual schema edits to `crates/lbm-scenario/
src/lib.rs` the same way — additive top-level JSON blocks, mechanical rebase.

---

## 9. Load-bearing code references (grounding index)

| Claim | File:line |
|---|---|
| Invariant step order in `run_span` (untouched by this spec) | `crates/lbm-core/src/backend.rs:258-300` |
| `run_staged_step` (gravity pre-pass ordering precedent for O1's new orchestration method) | `crates/lbm-core/src/solver.rs:1587-1606` |
| `Backend::Fields` open composite storage boundary (`g`/`h` precedent this spec's slots call into) | `crates/lbm-core/src/backend.rs:130-135` |
| `OutputSpec` / cadence model (`every`) reused for every new QOI's cadence | `crates/lbm-scenario/src/lib.rs:591-600` |
| `FieldKind` enum + exhaustive-match consumers (extension pattern for §2's new per-cell variants) | `crates/lbm-scenario/src/lib.rs:568-589`; `crates/lbm-cli/src/runner.rs:229-232,565-671,791-816` |
| `ProbeSpec` enum (extension pattern for `FaceFlux`) | `crates/lbm-scenario/src/lib.rs:552-566` |
| `validate()` warning pattern (extension pattern for §3.6's new validation rules) | `crates/lbm-scenario/src/lib.rs:620-690` |
| `Manifest` struct / `manifest.json` plumbing (extension point for `reactor_outputs`) | `crates/lbm-cli/src/runner.rs:187` (Manifest construction), `crates/lbm-cli/src/main.rs:59-138` |
| WVOF phase-field pre-pass slot (this spec's §1.1 step 2, called not redefined) | `docs/proposals/WVOF_IMPL_SPEC.md` §4.2 |
| WSCAL scalar sub-step slot (this spec's §1.1 steps 0/4, called not redefined) | `docs/proposals/WSCAL_PASSIVE_SPEC.md` §4.2 |
| FORCE_COMPOSITION_SPEC T4/T5 + R1-R4 composition-order invariants (this spec's §1.1 step 1, called not redefined) | `docs/proposals/FORCE_COMPOSITION_SPEC.md` §1, §2 |
| active-scalar-feedback.md §5.1 dt formulas (cited verbatim in §1.2/§6, not re-derived) | `docs/proposals/active-scalar-feedback.md` §5.1 |
| REQ FR-COUP-01/02/03/04/05 (this spec's normative source for §1) | `docs/REQ_STIRRED_REACTOR.md` §5 |
| REQ FR-IO-01..06 (this spec's normative source for §2, esp. ε_g metadata) | `docs/REQ_STIRRED_REACTOR.md` §6 |
| REQ FR-VOF-03 raw-φ rejection (this spec's §3.6 extends, does not re-decide) | `docs/REQ_STIRRED_REACTOR.md` §4.4 |
| REQ §8 VR-STR-05/07 acceptance criteria | `docs/REQ_STIRRED_REACTOR.md` §8 |
| REQ §11 DAG placement of W-COUP/W-IO ("incremental", "active set") | `docs/REQ_STIRRED_REACTOR.md` §11 |
| VALIDATION.md T17 VR-STR-05/07 rows | `docs/VALIDATION.md:349-351` |
| WVOF §7.1 shared-file mechanical-merge precedent (this spec's §8 coexistence rule) | `docs/proposals/WVOF_IMPL_SPEC.md` §7.1 |
| physics-discipline ban list (transport-absorbing clamp — governs §4.5/C9) | `.claude/skills/lbmflow-physics-discipline/SKILL.md` Rule 2 |

**Literature (decided references — all already cited by a sibling spec; no
new citation introduced by this tier, per the preamble discipline note):**
Fakhari, Mitchell, Leonardi & Bolster 2017 (PRE 96:053301) — ε_g/d_32
resolved-phasefield definitions, cited via WVOF_IMPL_SPEC.
Krüger et al. 2017 §8.3 — species ADE machinery underlying `C_i`, cited via
WSCAL_PASSIVE_SPEC.
Levenspiel, *Chemical Reaction Engineering* — RTD `E(t)`/`F(t)` definitions
and moments (§2.11), the standard reaction-engineering reference (no
alternative formulation needed; this is the textbook definition REQ FR-IO-05
names as "RTD").
Brackbill, Kothe & Zemach 1992 (JCP 100:335) — capillary dt lineage, cited
via active-scalar-feedback.md §5.1 (itself citing this for the `Δt_σ`
formula this spec's §1.2/§6 reuse verbatim).

---

## 10-line summary

- **O1** (`solver.rs` new orchestration method + `active_scalar_lagged` flag/
  validation) and **O2** (`reactor_outputs.rs` new module, pure reductions)
  are parallelizable, disjoint files — start together once one producing
  subsystem (W-SCAL and/or W-VOF) lands enough state to call into.
- **O3** (schema + CLI wiring: `species`/`reactions`/`equilibrium`/`thermal`/
  `mixing` blocks, new `FieldKind` variants, `manifest.json` `reactor_outputs`)
  depends on O2's reduction functions; schema authorship can start early,
  integration waits for O2.
- **O4** (adversarial tests, separate worktree) authors ALL of §5 now; its
  VR-STR-07-1/2 rows run TODAY against the already-landed W-ROT stirring
  path with zero dependency on O1-O3; C1-C6/W1 compile red and go green
  incrementally as O1-O3 and each sibling subsystem land.
- **Cross-tier dependency**: this tier is the DAG tip — `W0/W-GRAV/W-ROT/
  W-LES/W-PART (LANDED)` feed VR-STR-07 today; `W-VOF (pending, critical
  path) → W-SCAL active / W-BCTOP / W-REACT / ACIDBASE / W-BUB / thermal
  (all pending, none yet spec'd except W-VOF/WSCAL-passive)` feed C1-C6; the
  terminal gate C6 (full mini-reactor smoke) needs the ENTIRE active set.
- **STOP-RULE flags**: (1) C1 (active-scalar dt-halving) and C2/C3 (reaction
  subcycle/element-balance) cannot be run until W-REACT and at least one
  active-scalar feedback exist — authored now, explicitly red, not a
  physics gap in THIS spec; (2) this spec introduces no new physics term, so
  no Rule-1-row-3 stop condition applies to W-COUP/W-IO itself — any future
  stop-rule on this tier would originate in a SIBLING spec's physics, routed
  there, not fixed here by a shortcut in the orchestration/reduction layer.
