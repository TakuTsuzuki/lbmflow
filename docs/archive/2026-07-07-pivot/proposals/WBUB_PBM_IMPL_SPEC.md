# W-BUB Implementation Specification — Euler–Lagrange Point Bubbles + Population Balance (PBM) + Interfacial Mass Transfer

**Document ID**: SPEC-W-BUB-PB (rev.1, 2026-07-07).
**Scope**: the M-F item `W-BUB point bubbles + PBM + interfacial transfer` of
`docs/REQ_STIRRED_REACTOR.md` (§11 DAG; hard deps `W0, W-SCAL, W-EXT`;
Phase 2 / API-reserved). This is the **T4 tier** of the interface axis (REQ §1:
interface relaxation = `point-bubble`; the resolved `resolved-phasefield`
Allen–Cahn W-VOF is the T-fidelity reference the point-bubble path is
degradation-benchmarked against — VR-STR-RELAX). **This is the PRIMARY
full-tank aeration route** per the stirred-reactor goal: resolving every bubble
in an `O(10⁸)`-cell aerated tank is infeasible, so full-tank `ε_g`/`d_32`/`k_L a`
is delivered by tracking `O(10⁶–10⁷)` Lagrangian point bubbles two-way coupled
to the carrier, with a population balance closing the size distribution.
**Target core**: `crates/lbm-core` (D3Q19/D3Q27 carrier `f`; the D-track
one-way `particles.rs` Lagrangian machinery is the reused precedent;
`W-SCAL` `h`/`conc` fields carry the dissolved species the bubbles exchange
with).
**Acceptance**: VALIDATION.md **T17** rows **VR-STR-02b** (bubble swarm: `ε_g`
distribution, hindered rise, `d_32`), **VR-STR-02c** (aeration: `ε_g, d_32,
k_L a` vs correlations), **VR-STR-04** (the `k_L a` scalar row), and
**VR-STR-RELAX** (point-bubble-vs-resolved-phasefield degradation). Provisional
bands in §6.

This spec is **executable**: every force law, PBM kernel, and mass-transfer
correlation is decided, cited, derived, given a validity domain, and mapped to a
dedicated validation test with a provisional band. A follow-on codex order
should not need to re-derive any design decision here.

> **Discipline note (CLAUDE.md prime directive; PHYSICS.md;
> `.claude/skills/lbmflow-physics-discipline`).** Every force law, PBM kernel,
> and transfer correlation below is a **literature-backed closure** and
> therefore carries, without exception, all four Rule-1 artifacts: citation +
> derivation, stated validity domain (Eo–Mo–Re_b / α_g / Re_b ranges), its own
> validation test vs an analytic/correlation reference, and a PHYSICS.md entry
> (§8). There are **no band-calibrated constants** (every coefficient is a
> tabulated correlation value with a source), **no case-identity branches** (the
> Tomiyama contaminated/pure/distilled selector is a *system-property* input,
> not a per-test switch — §1.2), and **no transport-absorbing position clamps**:
> a bubble reaching a wall is handled by the wall-lubrication *force law* (§1.5)
> and the staircase contact reflection reused from `particles.rs`
> (`resolve_solid_contact`), which is a **documented wall model**, NOT a
> `.clamp()` that silently absorbs momentum. Gas holdup `ε_g` is built by a
> **smoothing kernel** (§1.7), not a bound. The high-Re_b Sherwood validity
> domain is flagged as a technical risk (§7). The mandatory PHYSICS.md entries
> are in §8; behavior anchors are in §6.5.

> **Multicomponent dependency note.** This spec references a species registry /
> multicomponent `C_k` ADE infrastructure. `WSCAL_MULTICOMPONENT_SPEC.md` does
> **not exist in the tree as of 2026-07-07** (only `WSCAL_PASSIVE_SPEC.md`, the
> single-component `h`/`conc` phase-1 delivery, is landed-designed). This spec
> is therefore designed against **WSCAL_PASSIVE_SPEC + the REQ §3 multicomponent
> form** (`C_k`, per-component `D_k`, `S_k^if`, Henry partition) and states the
> dependency explicitly (§0 B7, §5.1): W-BUB's per-species transfer needs a
> `Vec<ScalarField>`-shaped registry that WSCAL_PASSIVE §P10 reserved but did
> not build. The **single-component special case** (one gas species, e.g. O₂
> into water) is fully specified here and is the MVP; the multicomponent
> generalization is gated on the species-registry order and is API-reserved, not
> built by W-BUB itself.

---

## 0. Summary of decisions (read this first)

| # | Decision | Justification (short) |
|---|---|---|
| B1 | **Euler–Lagrange point bubbles** are the delivered aeration model: each bubble is a Lagrangian marker with position, velocity, per-bubble diameter, and a per-bubble multicomponent gas inventory, advanced by an explicit force balance (§1.1). The carrier `f` sees the bubbles only through (a) a two-way momentum scatter and (b) a gas-holdup field `ε_g`. | REQ §1 interface-relaxation `point-bubble`; REQ §4.4 FR-VOF-04. Full-tank aeration at `O(10⁸)` cells cannot resolve every bubble; the point-bubble route is the only tractable full-tank `ε_g/d_32/k_L a` producer. Resolved W-VOF is the RELAX reference, not the production aeration path. |
| B2 | **Force balance** = buoyancy + **Tomiyama drag** + **Tomiyama lift** (with the Eo-dependent sign change) + **added mass** + **wall lubrication** (§1.2–1.5). Turbulent-dispersion `F_TD` is API-reserved (needs SGS `k` — deferred with W-LES scalar coupling, §1.6). | REQ §3 point-bubble momentum equation verbatim: `m_b dv_b/dt = F_buoy + F_drag(Tomiyama) + F_lift + F_addedmass + F_walllub + F_TD`. Each force is a cited closure (§8). |
| B3 | **Two-way momentum coupling via a conservative scatter kernel** reusing the `particles.rs` Lagrangian-buffer + grid-sampling machinery (`sample_grid` trilinear, `ParticleSet` container pattern). The reaction force `−F_bubble→fluid` is scattered onto the Eulerian `force_field` with the **same interpolation weights** used to sample `u` at the bubble (the standard PSIC/PCM consistency requirement), guaranteeing discrete momentum conservation. | REQ §3 particles row "regularized reaction-force scatter; momentum conservation validated"; §4.5 FR-PART "reaction-force scatter kernel + momentum-conservation validation". The one-way `particles.rs` is the code precedent (`sample_grid`, trilinear weights, staircase wall); W-BUB adds the *scatter* (transpose of the *gather*). |
| B4 | **PBM = class method (discrete bins of fixed pivot sizes), Hounslow-style fixed-pivot with mass-conservative redistribution.** Kernels: **coalescence = Prince & Blanch (1990)** (turbulent + buoyancy + laminar-shear collision + film-drainage efficiency); **breakup = Luo & Svendsen (1996)** (turbulent eddy collision, energy-density criterion). `d_32 = Σ n_i d_i³ / Σ n_i d_i²` from the bin populations. | REQ §4.4 "Population balance (PBM) mandatory … Luo–Svendsen / Prince–Blanch kernels; mono-disperse point-bubble cannot support `d_32` acceptance." Class method chosen over QMOM (§4.1 justification): it carries the *actual* size distribution (needed for the spatially-resolved `d_32` field and for per-bin rise-velocity, which QMOM's abstract moments cannot give the Lagrangian tracker directly). |
| B5 | **`ε_g` via kernel smearing** of point bubbles onto the Eulerian grid: `ε_g(x) = Σ_b V_b W(x − x_b; h) / V_cell`, `W` a compact normalized smoothing kernel (Gaussian-like tri-weight or the trilinear volume kernel), width `h = 2Δx` (§1.7). This is a **smoothing kernel, NOT a clamp** — `ε_g` is bounded above by physics (packing), and the model *warns/switches* to resolved above a stated `α_g`, it does not clip transport. | REQ §6 FR-IO-01 `ε_g_bubble = Σ_b V_b W_kernel(x − x_b)/V_filter` verbatim (the point-bubble ε_g definition; carries filter width as metadata). |
| B6 | **Per-species mass transfer `ṁ_k = k_L a (C_k* − C_k)`** per bubble: `k_L` from a **bubble Sherwood correlation** (Higbie penetration `Sh = (2/√π)√(Pe_b)` for clean small bubbles, with the Frössling/Ranz–Marshall `Sh = 2 + 0.6 Re_b^{1/2} Sc^{1/3}` alternative for the rigid/contaminated regime — selector = the same interface-condition input as the drag law, §5.2); `a` = interfacial area from the bubble surface; `C_k*` = **Henry partition** of the per-bubble gas partial pressure. Gas composition EVOLVES: the per-bubble inventory `n_{k,b}` changes by `−ṁ_k`, and the dissolved `C_k` field gains `+ṁ_k` via a source scattered like B3. | REQ §3 `S^if` point-bubble = `k_L a(C*−C)`; §4.5 FR-VOF-05 "point-bubble = `k_L a(C*−C)`; Henry and Sherwood applicability explicit." Shares the Henry/Sherwood closures with the T3 resolved-interface transfer spec — see §5.4 coordination. |
| B7 | **Multicomponent gas per bubble** stored as a small per-bubble `Vec`/array of species moles `n_{k,b}`; couples to the W-SCAL dissolved-`C_k` fields through the species registry. **Single component is the MVP** (one gas species; e.g. O₂→water); multicomponent is API-reserved, gated on the (not-yet-existing) `WSCAL_MULTICOMPONENT_SPEC` species registry. | REQ §3 conservative multicomponent `C_k`; WSCAL_PASSIVE §P10 reserved the `Vec<ScalarField>`. The dependency is explicit; W-BUB does not build the registry. |
| B8 | **Hindered (swarm) rise** via a void-fraction correction to the drag: `C_D,swarm = C_D(1 − α_g)^{−n}` (Richardson–Zaki-type exponent from the bubble-swarm literature) — reduces terminal rise velocity as `α_g` grows. This is a **closure with a stated exponent and validity domain**, not a fitted constant (§1.8). | REQ §8 VR-STR-02b "hindered rise"; the swarm effect on rise velocity is a named acceptance behavior anchor. |
| B9 | **Lagrangian sub-step composed at solver-orchestration level, AFTER the carrier `f` step** (reads the just-updated physical `u` = F/2-corrected velocity — the same slot W-SCAL uses, `solver.rs` sub-step level). The invariant `f` pass order and the `Backend` trait are **untouched**; the two-way force is applied on the *next* carrier step through `force_field` (a one-step-lagged explicit two-way coupling, standard for point-particle LBM; §3). | REQ §5 FR-COUP-01 dataflow: `… → boundary → scalar ADE → reaction → particle integration`. Bubbles read the resolved `u`; the reaction force is composed into `force_field` for the next Guo-forced collide (FORCE_COMPOSITION_SPEC slot). |
| B10 | **Engine-agnostic bubble module** mirroring `particles.rs`: the bubble set is advanced from caller-supplied sampler + scatter closures and holds no solver reference. New file `crates/lbm-core/src/bubbles.rs` (force balance + swarm + wall), new file `crates/lbm-core/src/pbm.rs` (class-method coalescence/breakup), keeping W-BUB's core touch to two *new* files plus a thin solver-level orchestration hook — disjoint from `particles.rs`, `fields.rs` distribution slots, and the W-VOF/W-SCAL fields (§9). | CLAUDE.md minimal-scope + `particles.rs` precedent (deliberately engine-agnostic, no solver reference). Two new files minimize the conflict surface with in-flight W-VOF/W-SCAL orders. |
| B11 | **CPU-first (CpuScalar reference), GPU deferred.** The Lagrangian sub-step, scatter, PBM, and transfer are host-side over the staging `SoaFields`; GPU bubble tracking is a follow-on gated on B-1 (staged multi-buffer upload) exactly as W-SCAL/W-VOF defer GPU. | B-1 PARTIALLY RESOLVED; forcing GPU into Phase 2 W-BUB would stall. Lagrangian counts (`O(10⁶–10⁷)`) are ~1 GB (REQ NFR-01 "Particles 10⁷ × ~100 B = 1 GB (negligible)") — host-side is affordable at dev/validation scale. |

---

## 1. Governing force balance + PBM equations

Notation follows REQ §2 (`Eo, Mo, Re_b, We_b`), all in **lattice units** inside
the core (the unit-conversion layer maps SI ↔ lattice; §7.3). A bubble `b` has
position `x_b`, velocity `v_b`, diameter `d_b`, volume `V_b = π d_b³/6`, gas
density `ρ_g`, gas mass `m_g = ρ_g V_b`. The carrier is the *liquid* with
density `ρ_l`, kinematic viscosity `ν_l`, dynamic `μ_l = ρ_l ν_l`, surface
tension `σ`, gravity `g`. Slip velocity `w = u(x_b) − v_b` (liquid minus bubble;
`u` is the resolved F/2-corrected velocity sampled at `x_b`). Bubble Reynolds
`Re_b = |w| d_b / ν_l`.

### 1.1 The point-bubble momentum equation (decision B2)

The bubble's translational motion (REQ §3, verbatim structure):

```
(ρ_g + C_A ρ_l) V_b  dv_b/dt
   = F_buoy + F_drag + F_lift + F_addedmass* + F_walllub  [+ F_TD]          (1)
```

where the added-mass term is written on the LHS as an effective inertia
(`C_A ρ_l V_b`, decision §1.4) and `F_addedmass*` on the RHS carries the
*fluid-acceleration* part (`C_A ρ_l V_b Du/Dt`) so the two together are the
complete added-mass closure. Because `ρ_g ≪ ρ_l`, the bubble is
*inertia-light*: the effective mass is dominated by the added mass `C_A ρ_l V_b`
(the reason a bubble tracks the surrounding acceleration far more tightly than a
heavy particle — this is physically load-bearing, not a stiffness hack).

Each RHS force is defined below with citation + derivation + validity.

### 1.2 Drag — Tomiyama (1998) (decision B2)

The drag force:

```
F_drag = ½ C_D ρ_l A_b |w| w,     A_b = π d_b²/4                            (2)
```

**Drag coefficient** — Tomiyama, Kataoka, Zun & Sakaguchi (1998), the standard
single-bubble correlation parameterized by the **system contamination state**,
returning `C_D` as the max of a Reynolds-dependent Stokes/Schiller–Naumann
branch and an Eötvös-dependent distorted-bubble branch:

```
Contaminated (dirty water, industrial default):
  C_D = max[ (24/Re_b)(1 + 0.15 Re_b^0.687),  (8/3) Eo/(Eo + 4) ]           (3a)

Pure / slightly contaminated:
  C_D = max[ min( (24/Re_b)(1+0.15 Re_b^0.687), (72/Re_b) ), (8/3)Eo/(Eo+4) ] (3b)

Distilled / hyper-clean (mobile interface):
  C_D = max[ (16/Re_b)(1+0.15 Re_b^0.687), (8/3)Eo/(Eo+4) ]                  (3c)
```

with `Eo = Δρ g d_b² / σ` (REQ §2), `Δρ = ρ_l − ρ_g`.

**Derivation / basis.** The first argument is the Schiller–Naumann viscous drag
(the same rigid-sphere law already in `particles.rs::schiller_naumann_drag_correction`;
reused verbatim), valid at low–moderate `Re_b`; the second is the
Mendelson/Tomiyama distorted-cap-bubble limit where surface tension (`Eo`)
governs the terminal shape and drag. `max(·)` selects the governing regime —
this is Tomiyama's published correlation form, **not** a case switch: the branch
is chosen by the *physical* `Re_b, Eo` of the bubble, and the
contaminated/pure/distilled variant is a **fluid-system property** (an input
field of the scenario, e.g. `interface_condition: contaminated`), the same input
that selects the Sherwood correlation (§5.2). It is decided at scenario level
per the physical liquid, never per test case.

**Validity domain (frozen, PHYSICS.md §8):** `10⁻² ≲ Re_b ≲ 10⁵`,
`10⁻³ ≲ Eo ≲ 40`, `Mo` spanning air–water to viscous systems
(`Mo ∈ [10⁻¹¹, 10³]`). The industrial default is (3a) contaminated. Outside the
`Eo`/`Mo` box the model must warn (a validity flag, not a silent extrapolation).

### 1.3 Lift — Tomiyama et al. (2002) lift coefficient with Eo sign change (decision B2)

The shear-lift force on a bubble in a velocity-gradient field:

```
F_lift = −C_L ρ_l V_b ( w × ω ),   ω = ∇ × u(x_b)                          (4)
```

**Lift coefficient** — Tomiyama, Tamai, Zun & Hosokawa (2002), the correlation
that captures the experimentally observed **sign reversal** of the lateral
migration for large deformed bubbles:

```
        ⎧ min[ 0.288 tanh(0.121 Re_b),  f(Eo_d) ]      Eo_d < 4
  C_L = ⎨ f(Eo_d)                                       4 ≤ Eo_d ≤ 10
        ⎩ −0.27                                         Eo_d > 10                (5)
  f(Eo_d) = 0.00105 Eo_d³ − 0.0159 Eo_d² − 0.0204 Eo_d + 0.474
  Eo_d = Δρ g d_H² / σ ,  d_H = d_b (1 + 0.163 Eo^0.757)^{1/3}  (Wellek aspect ratio)
```

`Eo_d` uses the horizontal long-axis diameter `d_H` via the Wellek et al. (1966)
aspect-ratio correlation.

**Derivation / basis.** For small bubbles the lift is the classical Saffman/Auton
positive lift (`C_L > 0`, bubbles migrate toward the low-pressure/high-velocity
side); as the bubble deforms (large `Eo_d`) the wake asymmetry reverses the
lateral force (`C_L < 0`) — Tomiyama's tanh + cubic fit reproduces the measured
crossover near `Eo_d ≈ 4` and saturates at `−0.27`. The **sign change is the
physics**, and the correlation coefficients are Tomiyama's published values (no
recalibration). This sign reversal is a mandatory behavior anchor (§6.5): small
bubbles must migrate toward the wall in a rising swarm near a wall, large bubbles
toward the core.

**Validity domain:** `1.39 ≤ log₁₀ Mo ≤ ... ` per Tomiyama (2002) air–water /
Glycerol systems; `Eo_d ∈ [1.2, 5.7]` is the fitted range, extrapolated with a
warning beyond it. `Re_b ∈ [1.4, ...]`.

### 1.4 Added mass (decision B2)

```
F_addedmass = C_A ρ_l V_b ( Du/Dt |_{x_b} − dv_b/dt ),   C_A = ½ (sphere)     (6)
```

The `−C_A ρ_l V_b dv_b/dt` part is moved to the LHS effective inertia (eqn 1);
the `+C_A ρ_l V_b Du/Dt` part (material derivative of the liquid velocity at the
bubble, computed from stored `u` history / local advection) stays on the RHS.
`C_A = 1/2` is the exact potential-flow added-mass coefficient for a sphere
(Lamb, *Hydrodynamics*); for strongly deformed bubbles it is a mild
over-estimate but is the standard closure — validity `Eo ≲ 4` for the
sphere value, flagged otherwise. **No tuning.**

### 1.5 Wall lubrication — Antal / Tomiyama (decision B2; the wall is a MODEL not a clamp)

Near a wall a rising bubble feels a repulsive lateral force preventing
overlap; this is a **force law**, the physically correct way to keep bubbles off
the wall — categorically different from a `.clamp()` on position:

```
F_walllub = −C_W (ρ_l/2) |w_∥|² (d_b/2) ( 1/y_w² ... ) n_w                   (7)
```

using the **Tomiyama (1998) wall force** form
`C_W(Eo)` with the wall-distance decay
`F_WL = C_WL (d_b/2)( 1/y_w² − 1/(D−y_w)²)` between two walls (or the single-wall
`1/y_w²` term), `y_w` = distance to the wall, `n_w` = wall-normal, `w_∥` = slip
component parallel to the wall. `C_WL(Eo)` is Tomiyama's tabulated
`Eo`-dependent coefficient.

**Derivation / basis.** The lubrication film between the bubble and the wall
supports a pressure that pushes the bubble away, scaling as `1/y_w²`; Antal,
Lahey & Flaherty (1991) and Tomiyama (1998) give the closed forms. This force,
combined with lift, sets the near-wall void-fraction profile (the peak/coring
behavior). **The wall is resolved by this force, not by clamping `x_b`.** The
only position operation is the staircase-contact reflection reused from
`particles.rs::resolve_solid_contact` for the degenerate case of a bubble that
still reaches a solid cell (a documented wall reflection model, restitution
input, subdivided to prevent tunneling — NOT a transport-absorbing clamp; the
lubrication force is what normally keeps bubbles off the wall).

**Validity:** `y_w/d_b ≳ 0.5`; below that the point-bubble abstraction breaks and
the model must warn (the bubble should be resolved, RELAX regime).

### 1.6 Turbulent dispersion (API-reserved, deferred)

`F_TD = −C_TD ρ_l k ∇α_g` (Lopez de Bertodano) needs the SGS turbulent kinetic
energy `k` from W-LES. W-LES is landed but exposes `ν_t`, not `k`; deriving `k`
from `ν_t` is a separate closure. `F_TD` is **API-reserved** (a slot in the
force accumulation, written by nobody in this phase — the same discipline as the
W-SCAL reserved `F_b^scalar`). Its absence is a documented model limitation
(PHYSICS.md §8), not a hidden approximation.

### 1.7 Gas holdup `ε_g` by kernel smearing (decision B5)

```
ε_g(x) = Σ_b V_b W(x − x_b; h) / V_cell ,   Σ_cells W(x−x_b;h) V_cell = 1     (8)
```

`W` is a **compact, normalized smoothing kernel** — the delivered choice is the
**trilinear volume kernel** (the exact transpose of the `sample_grid` trilinear
gather, so scatter and gather share weights → discrete consistency with the
momentum scatter B3), optionally widened to a `h = 2Δx` tri-weight tent kernel
to reduce grid noise at high bubble counts. The normalization `Σ W V_cell = 1`
guarantees `∫ ε_g dV = Σ_b V_b` (total injected gas volume is exactly
represented — a mass-consistency test, §6). This is a **smoothing kernel with a
declared width `h`, carried as `ε_g` output metadata** (REQ §6 FR-IO-01: filter
width + averaging volume + time window mandatory). It is **not** a clamp: `ε_g`
is a diagnostic field bounded by physics, and above a stated `α_g` the config
*warns / recommends resolved-phasefield* (the point-bubble validity limit), it
does not clip.

### 1.8 Hindered (swarm) rise (decision B8)

The swarm reduces each bubble's rise velocity relative to a single bubble via a
void-fraction drag correction:

```
C_D,swarm = C_D,single · (1 − α_g)^{−n}                                      (9)
```

`α_g = ε_g(x_b)` sampled at the bubble; the exponent `n` is the
Richardson–Zaki-type bubble-swarm value (`n ≈ 2` for the distorted-bubble
regime per Garnier, Lance & Marié 2002 / Ishii–Zuber drift-flux; the exact `n`
is a **stated closure value with a validity domain**, frozen in PHYSICS.md after
the swarm characterization sweep — it is a physical correlation exponent, not a
band-fit constant). Higher `α_g` ⇒ larger `C_D` ⇒ lower terminal `v_b`
(hindering). This is the mechanism the VR-STR-02b hindered-rise anchor checks
(§6.5): mean swarm rise velocity must decrease monotonically with `α_g`.

### 1.9 Population balance (PBM) — coalescence + breakup (decision B4)

The number density `n(x, d, t)` of bubbles of diameter `d` evolves by
birth/death from coalescence and breakup (Ramkrishna 2000 PBE, spatially local
per Eulerian cell — the source terms operate on the bubble population *within a
cell/filter volume*):

```
∂n/∂t + advection = B_coal − D_coal + B_break − D_break                     (10)
```

**Coalescence — Prince & Blanch (1990).** Coalescence rate =
collision frequency × coalescence efficiency:

```
Γ_c(d_i, d_j) = [θ_ij^T + θ_ij^B + θ_ij^LS] · λ_ij                           (11)
  Turbulent collision:  θ_ij^T = C_1 (d_i+d_j)² (u_i'² + u_j'²)^{1/2},
                        u_i'  = 1.4 ε^{1/3} d_i^{1/3}   (inertial-range turbulent velocity)
  Buoyancy collision:   θ_ij^B = (π/4)(d_i+d_j)² |u_ri − u_rj|   (differential rise)
  Laminar-shear:        θ_ij^LS = (1/6)(d_i+d_j)³ γ̇   (mean shear rate γ̇)
  Efficiency:           λ_ij = exp( − t_ij / τ_ij )
       t_ij = { r_ij³ ρ_l / (16 σ) }^{1/2} ln(h_0/h_f)  (film-drainage time)
       τ_ij = r_ij^{2/3} / ε^{1/3}                       (contact time)
       r_ij = ½ (2/d_i + 2/d_j)^{-1}                     (equivalent radius)
```

`ε` = turbulent dissipation rate (from the resolved strain + SGS; §5.3),
`γ̇ = √(2 S:S)` (the FR-STRESS shear-rate field), `h_0, h_f` = initial/critical
film thickness (physical constants of the liquid, tabulated — air–water
`h_0 ≈ 10⁻⁴ m`, `h_f ≈ 10⁻⁸ m`, cited not fitted). `C_1 ≈ 0.089` is Prince &
Blanch's published turbulent-collision constant.

**Breakup — Luo & Svendsen (1996).** Breakup by turbulent-eddy collision, with
an energy-density criterion (an eddy breaks a bubble only if it carries enough
energy to create the new surface):

```
Ω_B(d:d_i) = ∫_{ξ_min}^{1}  0.923 (1−α_g) n_i (ε/d_i²)^{1/3}
             · (1+ξ)²/ξ^{11/3} · exp( − 12 c_f σ / (β ρ_l ε^{2/3} d_i^{5/3} ξ^{11/3}) ) dξ  (12)
  ξ = λ/d_i (eddy-size ratio), c_f = increase in surface area on breakup,
  β = 2.047 (Luo–Svendsen constant), ξ_min = 11.4 η/d_i (η = Kolmogorov scale)
```

All constants (`0.923`, `β = 2.047`, `ξ_min` coefficient) are Luo & Svendsen's
published values — **no recalibration**.

**Sauter mean diameter** from the discretized population (§4):

```
d_32 = Σ_i n_i d_i³ / Σ_i n_i d_i²                                          (13)
```

**Validity domains (PHYSICS.md §8):** both kernels assume **isotropic
inertial-range turbulence** (breaking eddies in the inertial subrange) — valid
for `d_b` between the Kolmogorov scale `η` and the integral scale, `α_g ≲ 0.2`
(dilute-to-moderate, the `(1−α_g)` factor is the leading dilute correction).
Outside, warn. The kernels are the accepted stirred-tank/bubble-column standard
(Prince–Blanch, Luo–Svendsen are the two REQ §4.4 names).

---

## 2. Data structures (Rust API, reusing particles.rs) (decisions B1, B3, B7, B10)

New engine-agnostic modules, mirroring `particles.rs` (no solver reference;
advanced from caller-supplied sampler + scatter closures).

### 2.1 `crates/lbm-core/src/bubbles.rs`

```rust
/// One Lagrangian point bubble, lattice units. Mirrors `particles::Particle`
/// but carries a per-bubble multicomponent gas inventory and the deformation
/// diameters the force laws need.
#[derive(Clone, Debug, PartialEq)]
pub struct Bubble {
    pub pos: [f64; 3],
    pub vel: [f64; 3],
    /// Equivalent spherical diameter d_b (volume-equivalent).
    pub d: f64,
    /// Per-species gas inventory, moles (or lattice mass units): n_{k,b}.
    /// Length == species registry length. Single-component MVP => len 1.
    pub gas: Vec<f64>,
    /// PBM size-class index this bubble currently belongs to (fixed-pivot).
    pub class: usize,
}

impl Bubble {
    pub fn volume(&self) -> f64 { std::f64::consts::PI * self.d.powi(3) / 6.0 }
    /// Wellek horizontal diameter for Eo_d (lift), eqn (5).
    pub fn d_horizontal(&self, eo: f64) -> f64 {
        self.d * (1.0 + 0.163 * eo.powf(0.757)).cbrt()
    }
}

/// Liquid-carrier + system parameters and the interface-condition selector.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum InterfaceCondition { Contaminated, Pure, Distilled }

/// Bubble container + carrier parameters, lattice units. Mirrors ParticleSet.
#[derive(Clone, Debug)]
pub struct BubbleSet {
    pub bubbles: Vec<Bubble>,
    pub rho_l: f64,
    pub rho_g: f64,
    pub nu_l: f64,
    pub sigma: f64,
    pub g: [f64; 3],
    pub interface: InterfaceCondition,
    pub c_added: f64,          // = 0.5 (sphere); explicit, not hidden
    pub swarm_exponent: f64,   // n in (9); frozen closure value
    pub restitution: f64,      // reused wall-contact reflection
}

/// Sampled carrier state at a bubble position (extends particles::Sample with
/// what the bubble force balance needs).
#[derive(Clone, Copy, Debug)]
pub struct CarrierSample {
    pub u: [f64; 3],            // resolved F/2-corrected velocity
    pub dudt: [f64; 3],         // material derivative Du/Dt (added mass)
    pub vort: [f64; 3],         // vorticity ω = ∇×u (lift)
    pub shear_rate: f64,        // γ̇ = √(2 S:S) (PBM laminar-shear collision)
    pub eps_turb: f64,          // turbulent dissipation ε (PBM)
    pub alpha_g: f64,           // local ε_g (hindered rise)
    pub y_wall: f64,            // distance to nearest wall (lubrication)
    pub n_wall: [f64; 3],       // wall normal
    pub solid: bool,
}

/// The reaction the bubble applies back to the carrier at its position:
/// -F_bubble→fluid (momentum) and +ṁ_k (species source into dissolved C_k).
#[derive(Clone, Debug)]
pub struct BubbleReaction {
    pub pos: [f64; 3],
    pub force: [f64; 3],        // scattered into force_field (B3)
    pub species_src: Vec<f64>,  // per-k source into conc/C_k fields (B6)
    pub gas_holdup: f64,        // V_b, scattered into ε_g field (B5)
}

impl BubbleSet {
    /// Advance all bubbles one lattice step; return the per-bubble reactions
    /// the caller scatters (transpose of the sampler's gather weights).
    /// `sample` supplies CarrierSample at arbitrary positions (trilinear,
    /// reusing particles::sample_grid machinery). Errors when a force-law
    /// validity domain is exceeded (Eo/Re_b/y_wall out of range) — a
    /// validation error, never a silent clamp.
    pub fn step<F>(&mut self, sample: F) -> Result<Vec<BubbleReaction>, BubbleError>
    where F: Fn([f64; 3]) -> CarrierSample { /* §3 */ unimplemented!() }
}

/// Validity-domain violation (mirrors particles::ParticleError shape).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BubbleError {
    pub bubble_index: usize,
    pub kind: BubbleErrorKind,   // DragEo, LiftEo, WallTooClose, ...
    pub value: f64,
    pub domain_max: f64,
}
```

### 2.2 `crates/lbm-core/src/pbm.rs`

```rust
/// Fixed-pivot size classes (Hounslow). Pivots d_i are geometric (d_{i+1}=2^{1/3}d_i
/// => volume-doubling), so coalescence of two class-i bubbles lands in class i+1
/// exactly (mass-conservative pivot). Diameters are lattice units.
#[derive(Clone, Debug)]
pub struct SizeGrid {
    pub pivots: Vec<f64>,        // d_i, ascending
}

/// Per-cell (or per-filter-volume) bubble-number populations by class.
#[derive(Clone, Debug)]
pub struct ClassPopulation {
    pub n: Vec<f64>,             // number per class, len == pivots.len()
}

impl ClassPopulation {
    /// Sauter mean diameter d_32 = Σ n_i d_i³ / Σ n_i d_i²  (eqn 13).
    pub fn sauter(&self, grid: &SizeGrid) -> f64 { /* … */ 0.0 }
}

/// Coalescence (Prince–Blanch, eqn 11) + breakup (Luo–Svendsen, eqn 12) kernels.
/// Pure functions of (d_i, d_j, ε, γ̇, α_g, σ, ρ_l, ν_l): no solver state,
/// GPU-portable (FR-EXT-01).
pub fn coalescence_rate(di: f64, dj: f64, eps: f64, shear: f64,
                        rho_l: f64, sigma: f64, /*film consts*/ h0: f64, hf: f64) -> f64;
pub fn breakup_rate(di: f64, eps: f64, alpha_g: f64,
                    rho_l: f64, nu_l: f64, sigma: f64) -> f64;

/// Advance the population one PBM step (fixed-pivot birth/death with
/// mass-conservative redistribution, eqn 10). Returns the updated populations;
/// total gas volume Σ n_i (π d_i³/6) is conserved to round-off (a test, §6).
pub fn advance_population(pop: &mut ClassPopulation, grid: &SizeGrid,
                          fields: &CellPbmFields, dt: f64);
```

### 2.3 Eulerian side fields (in `SoaFields<T>`, all `Option` — B6 invariance)

The bubble reactions land in Eulerian fields; only **one new distribution-free
`Option` field** is added, plus reuse of the existing `force_field` slot:

```rust
// crates/lbm-core/src/fields.rs, appended to SoaFields<T> (all Option => None
// is bit-identical to the bubble-free path, the B-6 invariance discipline):
/// Gas holdup ε_g (point-bubble kernel-smeared, eqn 8), compact core.
pub eps_g: Option<Vec<T>>,
/// Sauter mean diameter d_32 field (PBM, eqn 13), compact core. Diagnostic.
pub d32: Option<Vec<T>>,
```

The two-way momentum reaction reuses the **existing `force_field: Option<Vec<[T;3]>>`**
(`fields.rs:196`) — the same slot Guo forcing / gravity / (future) `F_b^scalar`
compose into (FORCE_COMPOSITION_SPEC). The per-species dissolved source lands in
the **W-SCAL `conc`/`h`** fields (single component) or the multicomponent `C_k`
registry (reserved). No new distribution set is introduced — bubbles are
Lagrangian, `ε_g`/`d32` are compact scalar diagnostics, so W-BUB adds **zero new
LBE populations** (contrast W-SCAL's D3Q7 `h`).

---

## 3. Solver-step slot + carrier coupling (decisions B3, B9)

### 3.1 The invariant step order (unchanged)

The landed carrier per-step order (CLAUDE.md invariant, `backend.rs` `run_span`):

```
collide → halo → stream → open BCs → boundary moments (update_moments)
```

is **untouched**. The Guo-forced collide already reads `force_field`
(`fields.rs:196`), so the two-way bubble force is applied by *writing*
`force_field` before the next collide — no new pass, no `Backend`-trait change.

### 3.2 Where the bubble sub-step slots in (decision B9)

Per solver step (REQ §5 FR-COUP-01 dataflow `… → boundary → scalar ADE →
reaction → particle integration`):

```
1. CARRIER f STEP (unchanged run_span). Produces ρ, u=(ux,uy,uz) [F/2-corrected].
2. W-SCAL scalar ADE sub-step (if present): dissolved C_k transported by u.
3. W-BUB LAGRANGIAN + PBM SUB-STEP (new, solver-level, after moments+scalar):
   a. For each bubble: sample CarrierSample at x_b (trilinear gather, reuse
      particles::sample_grid): u, Du/Dt, ω, γ̇, ε, α_g=ε_g(x_b), y_wall, n_wall,
      and per-species dissolved C_k at x_b.
   b. Force balance (1): buoyancy + Tomiyama drag (2,3) + Tomiyama lift (4,5) +
      added mass (6) + wall lubrication (7); hindered-rise drag correction (9).
      Integrate v_b, x_b explicitly (semi-implicit in the linear drag term,
      exactly the `particles.rs` (τ_p v + u + τ_p g)/(τ_p+1) update generalized
      to the full force set). Staircase wall contact reused (resolve_solid_contact).
   c. Per-species mass transfer (§5): ṁ_k = k_L a (C_k* − C_k); update gas
      inventory n_{k,b} -= ṁ_k dt; accumulate species source +ṁ_k for scatter.
   d. Emit BubbleReaction{ force = −(sum of hydrodynamic forces the fluid does on
      the bubble, i.e. drag+lift+added-mass reaction), species_src, gas_holdup=V_b }.
4. SCATTER (conservative, B3): for each reaction, scatter `force` into force_field,
   `species_src` into conc/C_k, `gas_holdup` into eps_g — ALL using the SAME
   trilinear weights the gather (3a) used at x_b (PSIC consistency => discrete
   momentum & species conservation).
5. PBM sub-step (per cell/filter): advance_population with cell ε, γ̇, α_g;
   reconcile per-bubble diameters d_b with the updated class populations
   (bubbles split/merge => Lagrangian markers added/removed, deterministic order);
   recompute d32 field (13).
```

### 3.3 Two-way coupling: one-step-lagged explicit (decision B9)

The bubble reads `u` from step-1 moments and writes `force_field`, which the
**next** carrier collide consumes via Guo forcing. This is a **one-step-lagged
explicit two-way coupling** — the standard point-particle-in-LBM scheme
(Nguyen–Ladd; Feng–Michaelides). The lag is `O(Δt)` and its error converges
under `Δt`-halving (a coupling-convergence test is the RELAX/COUP acceptance,
§6). Momentum conservation is *exact per step* because the reaction scattered
onto the fluid is the negation of the force the fluid exerted on the bubbles,
with identical interpolation weights (the scatter is the transpose of the
gather) — VR-STR-05 momentum-budget gate.

### 3.4 CPU-first / GPU staging (decision B11)

Phase-2 W-BUB is host-side over the staging `SoaFields` (CpuScalar reference:
the bit-exact oracle). GPU bubble tracking + on-device scatter is a follow-on
gated on B-1, exactly as W-SCAL D8 / W-VOF D8 defer GPU. The PBM kernels
(`pbm.rs`) are written state-free / GPU-portable (FR-EXT-01) so the follow-on can
lift them onto the device unchanged.

---

## 4. PBM discretization (decision B4)

### 4.1 Method choice — fixed-pivot class method (justification)

**Chosen: the Hounslow (1988) fixed-pivot class method** with geometric
volume-doubling pivots `V_{i+1} = 2 V_i` (`d_{i+1} = 2^{1/3} d_i`).

Rejected alternatives and why:
- **QMOM / DQMOM (moment methods)** track only a few moments (`m_0..m_5`) and
  reconstruct abstract quadrature nodes. They are cheaper but **cannot hand the
  Lagrangian tracker a per-size rise velocity or a resolved size distribution** —
  W-BUB needs the *actual* distribution because each Lagrangian bubble has a real
  diameter feeding its own force balance (rise velocity, lift sign, drag). A
  moment method would force a synthetic re-sampling of bubbles from moments every
  step, adding an unvalidated closure. The class method carries the distribution
  the tracker already needs.
- **Fully continuous PBE PDE in `d`** is accurate but expensive and overkill for
  the `~10–20` classes a stirred tank needs.

### 4.2 Fixed-pivot mechanics

- **Pivots**: `d_i = d_min · 2^{i/3}`, `i = 0..N_c` (`N_c ≈ 12–16`), so
  coalescence of two class-`i` bubbles produces exactly a class-`i+1` volume
  (no inter-pivot mass smearing for equal-size coalescence). Unequal-size events
  land between pivots and are **redistributed to the two bounding pivots
  conserving both number and mass** (Hounslow's two-property-conservative
  redistribution — this is the mass-conservation mechanism, tested to round-off
  §6, NOT a clamp).
- **Coalescence** (11): for each pivot pair `(i,j)`, `Γ_c(d_i,d_j) n_i n_j`
  removes from `i,j`, adds to the pivot(s) bounding `V_i+V_j`.
- **Breakup** (12): `Ω_B(d_i)` removes from `i`, adds the daughter distribution
  to lower pivots (binary breakup default: two equal daughters at
  `V_i/2` → the pivot bounding `d_i·2^{-1/3}`).
- **Coupling to the Lagrangian markers**: the per-cell `ClassPopulation` is the
  *aggregate*; the Lagrangian bubbles are the *carriers*. A coalescence event
  removes two markers and inserts one at their mass-weighted centroid with summed
  gas inventory; a breakup event removes one marker and inserts two. Insertions
  and removals are in **deterministic bubble-index order** (bit-reproducibility,
  FR-COUP-04). Marker count is bounded by capping the resolved marker density and
  representing sub-grid populations by parcel weights (a "computational parcel"
  weight `w_b` — standard for `O(10⁷)` bubbles; the weight is a *representation*
  count, carried explicitly, not a physics fudge).

### 4.3 `d_32` field

`d_32(x)` (eqn 13) is computed per filter volume from the local
`ClassPopulation` (or directly from the kernel-smeared marker diameters) and
written to the `d32` diagnostic field (§2.3). It is the VR-STR-02b/02c acceptance
quantity and the behavior anchor "`d_32` responds to breakup rate" target (§6.5).

### 4.4 Turbulent dissipation `ε` input

The kernels need `ε` (turbulent dissipation). From the resolved+SGS field:
`ε = 2(ν_l + ν_t) S:S` (the total viscous dissipation of the resolved strain
plus the SGS eddy viscosity from W-LES `ν_t`). This reuses the FR-STRESS shear
tensor `S` and the landed WALE `ν_t` — no new closure, a documented composition
(PHYSICS.md §8). Where W-LES is off, `ε` is the resolved-only `2ν_l S:S` with a
logged limitation (under-estimates sub-grid breakup, same posture as W-SCAL's
molecular-`D` fallback).

---

## 5. Per-species mass transfer `k_L a (C* − C)` (decision B6)

This is the **point-bubble branch of FR-VOF-05**. The resolved-interface branch
(normal jump + Henry + phase-wise diffusion) lives in the T3 interfacial-transfer
spec; the two **share the Henry and Sherwood closures** — see §5.4 coordination.

### 5.1 The per-bubble, per-species balance

Each bubble carries gas moles `n_{k,b}`; the dissolved liquid field is `C_k`
(W-SCAL `conc` for the single-component MVP; the multicomponent `C_k` registry
otherwise — the B7 dependency). Transfer per species:

```
ṁ_k = k_L,k · a_b · ( C_k* − C_k(x_b) )                                     (14)
  a_b = π d_b²                     (bubble interfacial area, one bubble)
  dn_{k,b}/dt = − ṁ_k             (gas inventory shrinks/grows as species leave/enter)
  species source into liquid: + ṁ_k  scattered to C_k at x_b (§3.4)
```

The **aggregate `k_L a`** reported for validation (VR-STR-02c/04) is
`k_L a = Σ_b k_L,k a_b / V_filter`, equivalently
`a = 6 α_g / d_32` (the standard interfacial-area density from holdup and Sauter
diameter — a derived identity, cited, that the test cross-checks against the
direct sum). This is the FR-IO / VR-STR-04 "`k_L a` (formula stated)"
requirement: `a = 6 ε_g / d_32`, `k_L` from §5.2.

### 5.2 `k_L` from a bubble Sherwood correlation (selector = interface condition)

```
k_L,k = Sh_k · D_k,liq / d_b ,    Sh_k = Sherwood number                    (15)

Mobile / clean interface (Distilled/Pure): Higbie penetration
   Sh = (2/√π) √(Pe_b) ,   Pe_b = Re_b · Sc_k ,   Sc_k = ν_l / D_k,liq       (15a)
Rigid / contaminated interface: Frössling / Ranz–Marshall
   Sh = 2 + 0.6 Re_b^{1/2} Sc_k^{1/3}                                        (15b)
```

**Derivation / basis.** Higbie (1935) surface-renewal / penetration theory:
for a mobile bubble interface the contact time is `d_b/|w|`, giving
`k_L = 2√(D/π t_c)` ⇒ the `√Pe_b` scaling. The Frössling/Ranz–Marshall form is
the rigid-sphere (contaminated, immobile interface) boundary-layer result with
the `Sh→2` diffusive floor. **The selector is the SAME `InterfaceCondition`
input as the drag law** (§1.2, §2.1) — a physical fluid-system property, not a
per-case switch: a contaminated liquid gets rigid drag (3a) *and* rigid Sherwood
(15b), consistently.

**Validity domain (technical risk, §7):** Higbie (15a) is derived for
low–moderate `Re_b, Sc`; at **high `Re_b`** (large vigorously-wobbling bubbles)
the penetration assumption degrades and (15a) can over-predict — the high-`Re_b`
Sherwood domain is a flagged risk (§7), warned when `Re_b` exceeds the
correlation's stated range. `Sc_k ≳ 1` (liquids). No extrapolation without a
warning.

### 5.3 Henry partition for `C_k*` (shared closure, §5.4)

The interface equilibrium liquid concentration:

```
C_k* = H_k · p_{k,b} ,   p_{k,b} = (n_{k,b}/Σ_j n_{j,b}) · P_b               (16)
```

`H_k` = Henry solubility coefficient (species+liquid property, tabulated —
e.g. O₂ in water), `p_{k,b}` = partial pressure of species `k` in bubble `b`
(mole fraction × bubble pressure `P_b`; `P_b` = ambient + hydrostatic +
Laplace `4σ/d_b`). **Gas composition evolves**: as `ṁ_k` depletes `n_{k,b}`
differentially per species (O₂ transfers faster than N₂), the mole fractions
and hence `p_{k,b}`, `C_k*` change — the multicomponent gas-side balance is the
`dn_{k,b}/dt = −ṁ_k` set (14), closed by (16). Total moles set `P_b V_b = (Σn) R T`
(ideal gas; isothermal MVP) so the bubble **shrinks as gas dissolves** —
`d_b` updates from `V_b = Σn_{k,b}/(P_b/RT)`, feeding back into every force law
and the PBM (a validated coupling, not a clamp).

### 5.4 Coordination with the T3 resolved-interface transfer spec

FR-VOF-05 splits interfacial transfer into resolved (T3) and point-bubble
(here). The **shared closures** are the **Henry coefficients `H_k`** (§5.3) and
the **Sherwood/Sc definitions** (§5.2) — these live in a single
species-properties table (the species registry, B7) referenced by both specs.
The **T3 spec owns** the resolved normal-jump flux + phase-wise diffusion
(`S_{k,liq}^if = −S_{k,gas}^if`); **this spec owns** the point-bubble
`k_L a (C*−C)` lumped form. Where they overlap (a scenario running hybrid
resolved+point regions), the ε_g double-count exclusion (REQ §6 FR-IO-01
`ε_g_total = ε_g_resolved + ε_g_bubble` with exclusion) applies — reserved,
not built here.

---

## 6. Validation plan mapped to T17 (decisions all)

Tests **authored adversarially by codex/Opus from this spec**, in a test
worktree that never shares with the implementation worktree (CLAUDE.md; REQ §8).
Each row = metric / reference / band / grid / steps / backend. Bands are
provisional MVP gates (T17 band governance: tightening always allowed, loosening
needs a recorded PHYSICS.md rationale). **Two-layer** per Rule 3: a scalar band
**and** a behavior anchor.

### 6.1 Bands + denominators (layer 1)

| ID | Test | Metric & band (denominator stated) | Grid / steps / backend | T17 row |
|---|---|---|---|---|
| **VB1** | **Single-bubble terminal velocity vs Grace** (Eo–Mo–Re map). One bubble in quiescent liquid, released, rise to terminal `U_t`; sweep `(Eo, Mo)` across the Grace (1976) regime map (spherical / ellipsoidal / spherical-cap). | `\|U_t,measured − U_t,Grace\| / U_t,Grace` **≤ 10%** per `(Eo,Mo)` point (denominator = Grace-chart value). Run all three `InterfaceCondition` variants; contaminated is the primary gate. | quiescent box `128×128×256`, single bubble, ~50k steps to terminal, CpuScalar | VR-STR-02a/02b |
| **VB2** | **Bubble-swarm hindered rise + ε_g distribution.** `N_b` bubbles released from a bottom patch into a column; measure mean swarm rise velocity `⟨v_b,z⟩` vs `α_g` and the vertical `ε_g(z)` profile. | (a) `⟨v_b,z⟩(α_g) / U_t,single` follows `(1−α_g)^{n}` within **±15%** (denominator = single-bubble `U_t` from VB1); (b) `ε_g` profile shape matches drift-flux `ε_g = j_g/(v_slip + j)` within **±20%** (denominator = drift-flux value). | column `64×64×256`, `N_b∈{10³,10⁴}`, `α_g∈{0.02,0.05,0.10}`, CpuScalar | VR-STR-02b |
| **VB3** | **`d_32` from PBM vs correlation.** Aerated stirred cell at fixed `ε` (dissipation) and `α_g`; PBM reaches steady `d_32`; compare to a published `d_32 ∝ We^{-0.6}` / Calderbank-type correlation for the same `ε, σ, ρ_l, α_g`. | `\|d_32,PBM − d_32,corr\| / d_32,corr` **≤ 25%** (denominator = correlation value). | stirred cell `128³`, PBM `N_c=14` classes, to steady `d_32`, CpuScalar | VR-STR-02c |
| **VB4** | **Aeration `k_L a` vs published correlation.** Sparged tank; dissolve gas into initially-degassed liquid; fit the `C_k(t)` rise to `dC/dt = k_L a (C*−C)` ⇒ measured `k_L a`; compare to a Van't Riet-type `k_L a = c (P/V)^a (v_s)^b` correlation. | `\|k_L a_measured − k_L a_corr\| / k_L a_corr` **≤ 25%** (denominator = correlation; matches REQ §8 "02c ±25%"). Cross-check `a = 6ε_g/d_32` against the direct area sum within **±10%**. | tank `128³`, sparger patch, ~100k steps, single-species O₂-in-water params, CpuScalar | VR-STR-02c / VR-STR-04 |
| **VB5** | **Point-bubble-vs-resolved RELAX degradation.** Same single-bubble and small-swarm cases run BOTH point-bubble (this spec) and resolved W-VOF (fidelity ref). | Relative degradation `\|q_PB − q_VOF\|/\|q_VOF\|` reported (NOT a hard band day-one — RELAX bands freeze when the extension lands) for `q ∈ {U_t, ε_g, d_32}`; provisional target ≤ **20%** on `U_t`. Denominator = the resolved-W-VOF value. | single bubble + `10²` swarm, matched `(Eo,Mo)`, CpuScalar | **VR-STR-RELAX** |
| **VB6** | **Two-way momentum conservation (VR-STR-05).** Closed periodic box, bubbles + carrier, no external forcing except gravity+buoyancy (which cancel for the total). | Total momentum drift `\|Σ(ρ u)+Σ m_b v_b at t − at 0\| / initial\|` **< 1e-10** (f64) — the scatter-is-transpose-of-gather guarantee; round-off only, NOT a band. Also FR-COUP-04 `probe_state_hash` bit-identity single-backend. | periodic `64³`, `10³` bubbles, 20k steps, CpuScalar | VR-STR-05 |
| **VB7** | **Coupling `Δt`-convergence (FR-COUP-01).** The one-step-lagged two-way coupling error under `Δt`-halving on a bubble-in-shear case. | terminal migration position converges at **order ≥ ~1** under `Δt→Δt/2` (denominator = Richardson-extrapolated limit). | shear cell, CpuScalar | VR-STR-05 / COUP |

### 6.2 Mandatory negative / consistency tests (REQ §8; Rule 3 ablation guards)

- **NEG-1 (sparger off ⇒ zero holdup — THE named negative test):** with the gas
  inlet `φ`/bubble-injection **off** (no bubbles injected), the `ε_g` field is
  **exactly zero everywhere** and `d_32` is undefined/zero — no phantom holdup.
  A mutant that leaks a nonzero `ε_g` with no bubbles FAILS. (This is the
  `sparger φ off → zero holdup` negative test the order requires.)
- **NEG-2 (lift sign-change ablation):** a mutant that drops the Eo-dependent
  sign reversal (forces `C_L > 0` always) must FAIL the near-wall void-profile
  anchor (§6.5) for large bubbles — proves the (5) sign change is load-bearing.
- **NEG-3 (drag interface-condition consistency):** a scenario declared
  `contaminated` must use (3a)+(15b) *together*; a test asserts the selector is
  a single system input, and that mixing (3a) drag with (15a) Higbie Sherwood is
  rejected/impossible (no per-quantity case switch).
- **NEG-4 (mass conservation, PBM):** total gas volume `Σ_i n_i (π d_i³/6)`
  under coalescence+breakup with **zero transfer** is conserved to round-off
  (the fixed-pivot two-property redistribution) — a mutant with naive
  nearest-pivot assignment (mass-losing) FAILS.
- **NEG-5 (transfer mass balance):** `Σ_k Δn_{k,b}` (gas lost) `= −∫ Σ_k ṁ_k dV`
  (dissolved gained), to round-off — species conservation across the interface
  (VR-STR-05 scalar total-mass, phase-wise).
- **NEG-6 (`ν_t`/`ε` leak guard):** with W-LES off, the PBM `ε` uses molecular
  `2ν_l S:S` only (no silent SGS contribution) — a test asserts no `ν_t` term
  leaks in (same posture as W-SCAL's molecular-`D` guard).

### 6.5 Behavior-validity review (mandatory, layer 2 — Rule 3 anchors)

After each run, before reporting, review the *observed pattern*, not just the
band (record in PHYSICS.md / track findings):

- **A1 — rise velocity monotonic in `d_b`:** in the small-bubble (viscous) branch,
  larger bubbles rise faster (`U_t` increases with `d_b`) up to the distorted-cap
  plateau — `assert!(U_t(d_hi) > U_t(d_lo))` in the ellipsoidal regime; the Grace
  map's non-monotonicity at the cap transition is itself an anchor.
- **A2 — `ε_g` peaks above the sparger:** the holdup field must be maximal in the
  plume above the injection patch and decay laterally — assert the vertical `ε_g`
  profile has its interior maximum in the sparger column, not at a wall/seam
  (boundary-artifact sweep: `ε_g` accumulation exactly at a wall is guilty until
  proven physical — the kernel smear must not pile at the domain edge).
- **A3 — `d_32` responds to breakup rate:** raising the dissipation `ε` (more
  vigorous stirring) must DECREASE steady `d_32` (more breakup) —
  `assert!(d_32(eps_hi) < d_32(eps_lo))`. This is the breakup-kernel ablation
  guard: disabling breakup must raise `d_32` beyond band.
- **A4 — lift sign reversal:** small bubbles migrate toward the wall in a rising
  near-wall swarm (`C_L>0`), large bubbles toward the core (`C_L<0`) — assert the
  cross-stream `ε_g` peak location flips across `Eo_d ≈ 4`.
- **A5 — hindered rise monotone:** `⟨v_b,z⟩` decreases monotonically with `α_g`
  (VB2) — the swarm slows itself.
- **Boundary-artifact sweep:** inspect `ε_g`, momentum, and species at every
  wall, sparger patch, top outlet, and MPI seam — accumulation exactly at a bound
  is an artifact until proven physical (the deposition edge-ring lesson).

Every run leaves a **visual artifact** (ε_g cross-section, bubble scatter,
`d_32` map — lbmflow-qa-viewer); codex generates the artifact and lists its path,
the reviewing session does the looking (physics-discipline Skill).

---

## 7. Stability, validity domains, and RISKS

### 7.1 Force-balance stability

- The added-mass effective inertia `(ρ_g + C_A ρ_l)V_b` keeps the bubble ODE
  well-conditioned even for `ρ_g ≪ ρ_l` (without it the bare `ρ_g V_b` inertia
  would make the drag term stiff). The linear drag term is integrated
  semi-implicitly (the `particles.rs` `(τ_p v + …)/(τ_p+1)` structure), stable
  for any `τ_p > 0`.
- Sub-cycling: if the bubble response time `τ_b` is `≪ Δt`, sub-step the bubble
  integration within one carrier step (FR-COUP-01 particle `Δt_p`). Stated, not
  hidden.

### 7.2 Point-bubble validity domain (when the abstraction holds)

The point-bubble model is valid only when `d_b/Δx ≲ 1` (sub-grid bubble) **and**
`y_w/d_b ≳ 0.5` (bubble not touching the wall) **and** `α_g ≲ 0.2`
(dilute-to-moderate, the PBM `(1−α_g)` and swarm `(1−α_g)^{-n}` leading
corrections). Outside — large resolved bubbles, dense packing — the config
**warns and recommends resolved-phasefield** (the FR-VOF-04 switching criteria
`d_b/W, d_b/Δx, Eo, Re_b, α_g, We_b`); it does not silently extrapolate.

### 7.3 RISKS (technical)

| # | Risk | Why | Mitigation / flag |
|---|---|---|---|
| **R1** | **PBM kernel → `d_32` uncertainty.** | Prince–Blanch and Luo–Svendsen kernels have order-of-magnitude scatter between published implementations (film-drainage `h_0/h_f`, `C_1`, `β` sensitivities); `d_32` can be off by more than the 25% band in some regimes. | The band VB3 is **provisional** and frozen only after the characterization sweep records the achieved value + rationale in PHYSICS.md (band governance). The kernels are the REQ-mandated pair; if 25% is unreachable with published constants, STOP-RULE (do NOT recalibrate constants to pass — that is a banned band-fit). |
| **R2** | **High-`Re_b` Sherwood validity.** | Higbie penetration (15a) is derived for moderate `Re_b`; at high `Re_b` (large wobbling bubbles) it over-predicts `k_L`, biasing `k_L a` high. | `k_L` correlation warns above its stated `Re_b` range; VB4 restricts the primary gate to the correlation's validity box; the high-`Re_b` domain is documented as a limitation (PHYSICS.md), NOT patched with a fudge factor. Flagged as the order's STOP-RULE candidate if the aeration band needs high-`Re_b` bubbles. |
| **R3** | **Marker-count / parcel-weight representation.** | `O(10⁷)` bubbles with coalescence/breakup inserting/removing markers can drift the represented number vs the PBM aggregate. | Parcel weight `w_b` is explicit and its conservation (Σ weighted volume = PBM Σ) is a test (NEG-4 extension); deterministic insertion order (FR-COUP-04). |
| **R4** | **One-step-lag two-way coupling** at high `α_g`. | The explicit lag can destabilize dense two-way coupling. | VB7 `Δt`-convergence; sub-cycling (§7.1); the `α_g ≲ 0.2` validity cap (§7.2). |
| **R5** | **`ε` (dissipation) input quality.** | PBM kernels depend strongly on `ε`; resolved-only `ε` under-estimates sub-grid breakup when W-LES is off. | `ε = 2(ν_l+ν_t)S:S` with W-LES on (landed); molecular-only fallback logged (NEG-6). |

**STOP-RULE readiness:** if VB3 (`d_32`) or VB4 (`k_L a`) cannot be met with the
**published** kernel/correlation constants, the correct outcome is the
Rule-4 stop-rule report (spec/band revision or a documented resolved-W-VOF
fallback for that regime) — **never** a recalibrated constant or a case branch.

---

## 8. Codex order breakdown (decisions B10, B11)

Five orders. **One order = one bundle = one dedicated git worktree**
(CLAUDE.md). Implementation and adversarial-test orders **never share a
worktree**. File boundaries are chosen so the four impl orders touch **disjoint
files** wherever possible; the two shared files (`fields.rs`, `solver.rs`) are
touched additively (`Option` fields / a sub-step hook) exactly like the
W-SCAL/W-VOF coexistence discipline (§9).

| Order | Scope | Primary files (conflict boundary) | Deps | DoD |
|---|---|---|---|---|
| **O1 — Point-bubble force balance** | `bubbles.rs` (NEW): `Bubble`, `BubbleSet`, `CarrierSample`, `BubbleReaction`, `InterfaceCondition`; buoyancy + Tomiyama drag (2,3) + Tomiyama lift (4,5) + added mass (6) + wall lubrication (7) + hindered-rise (9); explicit semi-implicit integrator reusing the `particles.rs` update structure + `sample_grid` + `resolve_solid_contact`; conservative two-way scatter (B3, transpose-of-gather) into `force_field`; `eps_g` slot + kernel smear (8). | `crates/lbm-core/src/bubbles.rs` (NEW); `fields.rs` (+`eps_g` Option field, +init line); `solver.rs` (bubble sub-step hook §3.2 + scatter). | W0 (landed), particles.rs (landed) | VB1 (Grace `U_t` ≤10%), VB6 (momentum conservation to round-off), NEG-1 (sparger-off ⇒ ε_g≡0), NEG-2 (lift-sign ablation), A1/A4/A5 anchors, green on CpuScalar; `eps_g=None` bit-identity (B-6). PHYSICS.md §8 force-law entries landed. |
| **O2 — PBM (coalescence + breakup + d_32)** | `pbm.rs` (NEW): `SizeGrid`, `ClassPopulation`, `coalescence_rate` (Prince–Blanch, 11), `breakup_rate` (Luo–Svendsen, 12), `advance_population` (fixed-pivot Hounslow, mass-conservative), `sauter` (13); Lagrangian marker split/merge reconciliation (§4.2); `d32` field slot + population→`d_32` (13); `ε` composition `2(ν_l+ν_t)S:S` (§4.4). | `crates/lbm-core/src/pbm.rs` (NEW); `fields.rs` (+`d32` Option field, +init line); `bubbles.rs` (marker split/merge hook — coordinate with O1 via a shared trait, land O1 first). | O1 | VB3 (`d_32` ≤25%), NEG-4 (PBM mass conservation to round-off), A3 anchor (`d_32` decreases with `ε`); `d32=None` bit-identity. PHYSICS.md §8 kernel entries + frozen `d_32` band rationale. |
| **O3 — Per-species `k_L a` transfer** | Transfer (14) `ṁ_k = k_L a(C*−C)`; `k_L` Sherwood (15,15a/b) selector shared with drag; Henry `C*` (16); evolving gas inventory `n_{k,b}` + bubble-shrink `d_b` update; species source scatter into W-SCAL `conc` (single-component MVP) / `C_k` registry (reserved); aggregate `k_L a = 6ε_g/d_32 · k_L` + direct-sum cross-check. | `bubbles.rs` (transfer methods on `BubbleSet` — appended fns; coordinate with O1); `solver.rs` (species-source scatter into `conc` sub-step — after W-SCAL scalar step). | O1, W-SCAL (`conc` field) | VB4 (`k_L a` ≤25%), NEG-5 (transfer mass balance to round-off), R2 high-`Re_b` warning path tested; single-component green on CpuScalar. PHYSICS.md §8 Sherwood/Henry entries + coordination note w/ T3 transfer spec. |
| **O4 — Scenario / CLI / outputs (`ε_g`, `d_32`, `k_L a`)** | Scenario schema: sparger/bubble-injection config (gas volumetric flow, `InterfaceCondition`, species+`H_k`/`D_k`), `inlet_phase: gas\|liquid` (never raw φ, FR-VOF-03), size-grid config, validity-warning wiring (`d_b/Δx`, `y_w/d_b`, `α_g`, `Re_b` — §7.2); CLI/VTI output of `ε_g` (with filter-width/window metadata, FR-IO-01), `d_32`, aggregate `k_L a`, bubble scatter dump. | `crates/lbm-scenario/src/lib.rs` (schema + validation); `crates/lbm-cli` (outputs) — **disjoint from O1–O3 core files**. | O1, O2, O3 | Scenario round-trips; `ε_g`/`d_32`/`k_L a` written with metadata; NEG-1 reproducible from CLI (sparger-off scenario ⇒ zero holdup); validity warnings fire. |
| **O5 — Adversarial validation authorship** | ALL of §6: VB1–VB7 + NEG-1..6 + behavior anchors A1–A5. Authored **from this spec**, not from the impl. Grace map, drift-flux, `d_32`/`k_L a` correlation references encoded; freeze provisional bands in VALIDATION.md T17 VR-STR-02b/02c/04/RELAX. | `crates/lbm-core/tests/wbub_*.rs`, `crates/lbm-core/tests/wpbm_*.rs`, `crates/lbm-scenario/tests/*` (NEW files only — no impl-file conflict). | (spec only) | Tests compile red vs stub, go green as O1–O4 land; bands frozen with recorded rationale; every run lists a visual-artifact path (physics-discipline Skill). Runs concurrently from the start in its own worktree. |

**Critical-path ordering:** `O1 → O2`, `O1 → O3` (O2, O3 parallel after O1);
`O4` after O1–O3; **O5 runs concurrently from the start** (separate test
worktree). GPU (host→device bubble tracking + scatter) and `F_TD`
turbulent-dispersion are follow-on orders, out of this plan's scope.

**STOP-RULE flags carried into orders:** O2 (`d_32` band VB3) and O3 (`k_L a`
band VB4) each embed the Rule-4 clause verbatim — if the published-constant
kernel/correlation cannot meet the provisional band, the order emits the
stop-rule report (spec/band revision), it does **not** recalibrate a constant or
add a case branch.

---

## 9. Coexistence with W-SCAL / W-VOF / particles.rs (structural summary)

W-BUB is designed to mount alongside the in-flight and landed subsystems
without structural conflict:

- **New files:** `bubbles.rs`, `pbm.rs` are wholly new — zero conflict with
  `particles.rs` (which W-BUB *reuses* by pattern, not by editing: `sample_grid`,
  the `(τ_p v + …)/(τ_p+1)` integrator shape, `resolve_solid_contact`, the
  `Error`-on-validity-domain posture are all re-instantiated for bubbles, leaving
  the one-way particle code untouched).
- **`fields.rs`:** W-BUB adds two `Option` compact scalar fields (`eps_g`,
  `d32`) — **additive, `None`-default, bit-identical** to the bubble-free path
  (B-6 invariance). This is the same additive pattern as W-SCAL's
  `h`/`htmp`/`conc` and W-VOF's `g`/`gtmp`/`phi`; a textual merge at the struct
  field list / `new()` initializer is trivial both-add. W-BUB adds **no new LBE
  distribution set** (bubbles are Lagrangian; `ε_g`/`d32` are diagnostics), so it
  does not contend with the D3Q7 `h` or D3Q19 `g` distribution slots at all.
- **`solver.rs`:** W-BUB's bubble sub-step slots **after** the carrier `f` step
  and **after** the W-SCAL scalar ADE sub-step (§3.2) — a distinct
  orchestration slot from W-VOF's phase-field *pre*-pass (before `f`) and the
  carrier step itself. Three+ disjoint slots (pre-pass / f-step / scalar / bubble)
  → no semantic contention; textual proximity only, resolved by ordering the
  calls per each spec's §3/§4.
- **`force_field`:** W-BUB scatters the two-way reaction into the **existing**
  `force_field` (Guo/gravity/`F_b^scalar` slot), composing in the
  FORCE_COMPOSITION_SPEC frozen order. This is where W-BUB's force first
  interacts with W-VOF's `F_s` and gravity — an *additive accumulation* into the
  same buffer, in the frozen summation order (before the Guo half-force), each
  contributor with its own provenance. No structural change to the force path.
- **W-VOF dependency:** the RELAX reference (VB5) needs resolved W-VOF; W-BUB's
  *production* path does not depend on W-VOF landing (it is the sub-grid
  alternative). The hybrid `ε_g_total` double-count exclusion (REQ §6) is
  reserved for when both run.
- **W-SCAL dependency:** the transfer source (O3) writes into W-SCAL `conc`;
  single-component MVP needs only the landed-designed W-SCAL `conc` field. The
  multicomponent `C_k` registry (B7) is the one **hard external dependency** and
  is API-reserved (the `Vec<f64> gas` per bubble + `Vec` species source already
  shaped for it), gated on the not-yet-existing `WSCAL_MULTICOMPONENT_SPEC`.

---

## PHYSICS.md entries (mandatory, landed with the respective orders)

Copy into PHYSICS.md §1 (stack) and §2 (decisions), one per closure, no field
optional (physics-discipline Rule 1 template):

> **Point-bubble drag — Tomiyama (1998) (`bubbles.rs::drag`, O1).**
> Form: `F_drag = ½ C_D ρ_l (π d_b²/4)|w|w`, `C_D` = eqns (3a/b/c) by
> `InterfaceCondition`. Source: Tomiyama, Kataoka, Zun & Sakaguchi, JSME Int. J.
> 41(2) 1998. Validity: `10⁻²≲Re_b≲10⁵`, `10⁻³≲Eo≲40`, `Mo∈[10⁻¹¹,10³]`;
> contaminated (3a) industrial default. Validation: `wbub_terminal.rs`
> (VB1, Grace map ≤10%). Interacts: reused Schiller–Naumann branch from
> `particles.rs`; selector shared with Sherwood (O3).

> **Point-bubble lift — Tomiyama (2002) (`bubbles.rs::lift`, O1).**
> Form: `F_lift = −C_L ρ_l V_b (w×ω)`, `C_L` = eqn (5) with Eo-sign change and
> Wellek `d_H`. Source: Tomiyama, Tamai, Zun & Hosokawa, Chem. Eng. Sci. 57 2002;
> Wellek et al. 1966. Validity: `Eo_d∈[1.2,5.7]` fitted, warn beyond; sign
> reversal at `Eo_d≈4`. Validation: `wbub_lift.rs` (A4 anchor, NEG-2 ablation).
> Interacts: needs resolved vorticity `ω=∇×u`.

> **Added mass (`bubbles.rs::added_mass`, O1).** Form: eqn (6), `C_A=½` sphere
> (Lamb). Validity: `Eo≲4` for sphere value, warn beyond. Validation:
> `wbub_terminal.rs` transient. Interacts: LHS effective inertia (stiffness).

> **Wall lubrication — Tomiyama/Antal (`bubbles.rs::wall_lub`, O1).** Form:
> eqn (7), `1/y_w²` decay, `C_WL(Eo)`. Source: Antal, Lahey & Flaherty 1991;
> Tomiyama 1998. Validity: `y_w/d_b≳0.5`, warn/switch-to-resolved below.
> Validation: near-wall void-profile (A4). Interacts: this is the wall MODEL;
> `resolve_solid_contact` reflection is the degenerate backstop, NOT a clamp.

> **Hindered swarm rise (`bubbles.rs::swarm`, O1).** Form: eqn (9),
> `C_D,swarm=C_D(1−α_g)^{-n}`. Source: Richardson–Zaki / Garnier, Lance & Marié
> 2002; Ishii–Zuber drift-flux. Validity: `α_g≲0.2`; `n` frozen after VB2 sweep.
> Validation: `wbub_swarm.rs` (VB2, A5). Interacts: `α_g=ε_g(x_b)`.

> **Gas holdup smear (`bubbles.rs::eps_g`, O1).** Form: eqn (8), normalized
> trilinear/tent kernel width `h=2Δx`, `∫ε_g=Σ_b V_b`. Source: REQ §6 FR-IO-01
> point-bubble ε_g; SPH-type conservative smoothing. Validity: `α_g≲0.2`.
> Validation: NEG-1 (sparger-off⇒0), A2 (peaks above sparger), mass-consistency.
> Interacts: smoothing kernel, NOT a clamp; carries filter-width metadata.

> **PBM coalescence — Prince & Blanch (1990) (`pbm.rs::coalescence_rate`, O2).**
> Form: eqn (11), turbulent+buoyancy+shear collision × film-drainage efficiency;
> `C_1≈0.089`, `h_0,h_f` tabulated liquid constants. Source: Prince & Blanch,
> AIChE J. 36(10) 1990. Validity: isotropic inertial-range turbulence,
> `η<d_b<L_int`, `α_g≲0.2`. Validation: `wpbm_d32.rs` (VB3, A3, NEG-4).

> **PBM breakup — Luo & Svendsen (1996) (`pbm.rs::breakup_rate`, O2).**
> Form: eqn (12), turbulent-eddy collision energy-density criterion; `β=2.047`,
> `0.923`, `ξ_min=11.4η/d_i`. Source: Luo & Svendsen, AIChE J. 42(5) 1996.
> Validity: as coalescence. Validation: `wpbm_d32.rs` (VB3, A3). Interacts:
> `ε=2(ν_l+ν_t)S:S` (W-LES `ν_t`); `d_32=Σn_i d_i³/Σn_i d_i²`.

> **Bubble mass transfer — Higbie/Frössling Sherwood + Henry
> (`bubbles.rs::transfer`, O3).** Form: eqns (14–16), `k_L=Sh D/d_b`,
> `Sh`=(15a) Higbie mobile / (15b) Frössling rigid by `InterfaceCondition`,
> `C*=H_k p_{k,b}`, evolving `n_{k,b}`, shrinking `d_b`. Source: Higbie 1935;
> Ranz & Marshall 1952; Henry's law. Validity: `Sc≳1`; **Higbie high-`Re_b`
> over-predicts — RISK R2**, warn beyond correlation range. Validation:
> `wbub_kla.rs` (VB4 ≤25%, NEG-5). Interacts: point-bubble branch of FR-VOF-05;
> shares Henry/Sherwood with T3 resolved-transfer spec; writes W-SCAL `conc`/`C_k`.

---

## Literature (decided references)

- **Grace 1976** (Trans. IChemE 54:167) — single-bubble terminal-velocity /
  shape regime map (Eo–Mo–Re), the VR-STR-02a `U_t` reference.
- **Tomiyama, Kataoka, Zun & Sakaguchi 1998** (JSME Int. J. 41:472) — bubble
  drag correlation (contaminated/pure/distilled), wall force.
- **Tomiyama, Tamai, Zun & Hosokawa 2002** (Chem. Eng. Sci. 57:1849) — lift
  coefficient with Eo-dependent sign reversal; **Wellek, Agrawal & Skelland
  1966** (aspect ratio `d_H`).
- **Antal, Lahey & Flaherty 1991** (Int. J. Multiphase Flow 17:635) — wall
  lubrication force.
- **Prince & Blanch 1990** (AIChE J. 36:1485) — bubble coalescence kernel.
- **Luo & Svendsen 1996** (AIChE J. 42:1225) — bubble breakup kernel.
- **Hounslow, Ryall & Marshall 1988** (AIChE J. 34:1821) — fixed-pivot
  discretized population balance (mass-conservative class method).
- **Ramkrishna 2000**, *Population Balance* — the PBE and discretization theory.
- **Higbie 1935** (Trans. AIChE 31:365) — penetration-theory `k_L`;
  **Ranz & Marshall 1952** — Frössling/Ranz–Marshall Sherwood.
- **Richardson & Zaki 1954**; **Garnier, Lance & Marié 2002** (Exp. Therm.
  Fluid Sci. 26:811); **Ishii & Zuber 1979** — hindered/swarm rise, drift-flux.
- **Lamb 1932**, *Hydrodynamics* — added-mass `C_A=½`.
- **Nguyen & Ladd 2002**; **Feng & Michaelides 2004** — point-particle two-way
  coupling in LBM (scatter-as-gather-transpose momentum conservation).
- REQ_STIRRED_REACTOR §2 (nondim), §3 (point-bubble momentum, `S^if`), §4.4
  (FR-VOF-04 PBM mandatory), §4.5 (FR-PART), §6 FR-IO-01 (`ε_g` defn), §8
  (VR-STR-02b/c, RELAX); WSCAL_PASSIVE_SPEC (the `conc`/`h` this transfers into;
  the format mirrored here); FORCE_COMPOSITION_SPEC (the `force_field` slot);
  `crates/lbm-core/src/particles.rs` (the reused Lagrangian precedent).
```
