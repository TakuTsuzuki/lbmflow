# Requirements Specification (final): rotating-boundary, high-density-ratio two-phase, LES-coupled 3D multiphysics LBM solver

**Document ID**: REQ-M-F-STR / **Revision**: rev.3 (rev.1 = all 48 items of the codex adversarial review reflected / rev.1a = reflects the PM decision "default is fidelity-first, relaxation is a bolt-on extension point" / rev.1b = PM integration: domain-neutralization of the title, follow-up to the core rename, addition of the §7 memory-budget table, VALIDATION T17 wiring / rev.2 = all 11 items of the codex second review adopted: new relaxation-equivalence validation VR-STR-RELAX, clarification of scope wording, variable-σ surface-tension convention, addition of F_b^scalar, budget-table arithmetic fixes, etc. — see docs/proposals/req-round2-findings.md / **rev.3** = competitive-review triage diff merged (authored as "rev.1c" against rev.1b, layered here onto rev.2): P1 population balance, P2 §4.8 extension contracts, P3 FR-IO-05/06, P4 reference datasets, P5 product-layer scope note, §11 implementation dependency DAG. New content in English per the 2026-07-05 language directive. / **rev.4** = external numerical-physics review (REV-CFD-*, filed against rev.1a) triaged and merged: all 4 Critical fixes (sparger phase inversion, Allen–Cahn/continuity mass-flux consistency, non-equilibrium stress stage convention, forcing second-moment single-definition), Ca_spurious dimensional fix, Pe_N/Pe_tip split, active-scalar predictor–corrector dataflow, precision-profile enum, conservative scalar forms, four-way contact contract (extension-gated), viscosity-interpolation & σ(κ,β) freeze, ε_g processing definitions, and **provisional numeric acceptance bands** under a band-governance rule. Disposition table in TESTING_NOTES.md; MJ-008 was already fixed in rev.2. This file is PM-owned for translation — excluded from the bulk translation session.)
**Positioning**: M-F (vertical feature) of `docs/PLAN.md` / upstream requirement onto `ARCHITECTURE_V2.md`.
Validation acceptance is [VALIDATION.md](VALIDATION.md) **T17** (wiring VR-STR-01 to 07).
**Target core**: `lbm-core` (formerly lbm-core2. D3Q19/D3Q27, CpuScalar/CpuSimd/wgpu, MPI partitioning)
**Representative application problem**: stirred-tank reactor (functional requirements are defined domain-neutrally; the validation benches of §2 and §8
concretize this application. A region that squarely overlaps the core use case of M-Star CFD)

## 0. Review-reflection policy

- All Critical/Major equation, sign, and definition bugs (non-equilibrium stress evaluation, `τ_eff`, MRF apparent forces, surface tension, conservative Allen-Cahn, Np/N_Q, dimensionless numbers) were fixed and their coefficients frozen.
- Missing physics (top-side degassing, hydrostatic well-balanced, turbulent-Schmidt-number SGS scalar, particle SGS dispersion, bubble-induced turbulence BIT, contact angle, initialization/spin-up, memory budget) was added as new requirements.
- **Delivery scope (clarified in rev.2, enumeration finalized in rev.4)**: **the fidelity-default subsystem group (the "default" column of the §1 table) is implemented together simultaneously**. Relaxation extensions (MRF, point-bubble, one-way, block-AMR, aggressive f32) in the first release **only reserve the trait boundary, config schema, and validation items (VR-STR-RELAX)**, with implementation bolted on later. codex's "mode split" is dropped, not as a delivery phase, but as a **runtime-mode mutual-exclusion constraint** into the configuration matrix of §1. Modes that physically conflict within the same computation (MRF+IBM in the same zone, phase-averaging+MRF, etc.) are rejected at config validation.
  **Initial delivery (release gate)**: IBM-inertial, resolved-phasefield, active scalar
  (predictor–corrector coupling, §5), two-way particles, uniform grid, fidelity
  precision profile (`mixed_safe`, §7). **Phase 2 (API-reserved now, implemented
  later, accepted via VR-STR-RELAX)**: MRF-frozen-rotor, point-bubble (+PBM),
  one-way particles, four-way particle contact (FR-PART-04..06), block-AMR,
  aggressive f32 (`mixed_fast`), hybrid interface, thermal axis.
  "All subsystems simultaneously" means all physics axes exist in the initial
  delivery — not all modes of each axis. (REV-CFD-MJ-008; same substance as the
  rev.2 fix, list made explicit.)
- Unverifiable wording ("guarantee", "naturally tolerated", "stably integrated") was replaced with measurable errors, conservation laws, and applicability ranges.
- AMR was demoted to "the initial version is based on a uniform grid; AMR is an advanced option (enabled after defining coarse-fine conservative interpolation, the time-step ratio, and validation problems)" (#29).

---

## 1. Runtime-mode configuration matrix (mutual exclusion) / design principles

**Design principle (fidelity default, relaxation is a bolt-on extension point)**: each mode axis is abstracted with a strategy/trait. **The default for every axis is the fidelity-first implementation (= reference solution)**, and low-cost approximations (MRF, point-bubble, one-way, AMR, aggressive f32) are added as **bolt-on extensions** behind the same trait. The structure allows swapping without changing the core coupling loop (§5), and relaxation modes are validated against the corresponding fidelity reference solution by tolerance (thresholds defined in §8 VR). The fidelity default has the maximum computational cost (the §7 memory budget is sized on the fidelity configuration).

All modes are implemented, but exactly one is selected exclusively from each axis per computation. The config-validation layer rejects inconsistent combinations.

| Axis | Default (fidelity-first) | Relaxation extension (later) / reference tier | Exclusion constraints / notes |
|---|---|---|---|
| Rotation | `IBM-inertial` (unsteady, time-accurate) | relaxation: `MRF-frozen-rotor` (steady approximation) / reference tier: `sliding-overset` | MRF cannot be combined with IBM moving blades (#6). Phase-averaged statistics only for IBM/overset (#37). |
| Interface | `resolved-phasefield` (conservative Allen-Cahn. Interface / mass-transfer fidelity first) | relaxation: `point-bubble` (Euler-Lagrange) / `hybrid` | Switching decision is `d_b/Δx, d_b/W, Eo, Re_b, α_g, We_b` (#12). hybrid defines interphase mass / momentum / scalar conservation laws (§5). |
| Scalar | `active` (feedback to σ, viscosity, density, [temperature] enabled) | relaxation: `passive` (feedback opt-out) | active's feedback targets and stabilization made explicit (#13). |
| Particle coupling | `two-way` (`four-way` at high `α_p`) | relaxation: `one-way` | `α_p` / mass-loading threshold (#16). Accompanied by a reaction-force scatter kernel and momentum-conservation validation. |
| Precision | Fidelity profile: `f64` for the interface neighborhood, conserved quantities, torque, interface curvature, and reductions; `f32` only for the far-field bulk | relaxation: aggressive `f32` / reference tier: all `f64` | #32, §7. |
| Grid | `uniform` (fully resolved at the required resolution) | relaxation: `block-AMR` (when it satisfies coarse-fine conservative interpolation, the time-step ratio, and validation) | #29. AMR is bolted on later owing to implementation risk. |

---

## 2. Representative application problem, representative quantities, dimensionless numbers (stirred-tank reactor)

This section concretizes the **representative application** of the §8 validation benches; the functional requirements of §4
apply generally to rotating boundaries, high-density-ratio two-phase, LES, and scalar/particle coupling (consistent with the neutralization of the title).

A 3D cylindrical (or prismatic) vessel. A continuous phase (Newtonian / non-Newtonian liquid), a dispersed gas phase from a bottom sparger (`ρ_l/ρ_g ≈ 10³`, `μ_l/μ_g ≈ 10²`), a rigid-body rotating impeller at constant angular velocity `Ω`, suspended particles near neutral buoyancy, and multiple scalars with interfacial mass transfer and liquid-phase reactions. Target observables: time / phase-averaged 3D velocity field, shear-stress field based on the second invariant of the strain rate, particle Lagrangian cumulative shear exposure, gas holdup `ε_g`, and dissolved-scalar concentration field.

### 2.1 Dimensionless-number definitions (representative quantities fixed, #26, #23, #24)

Representative rotation rate `N = Ω/(2π)` [rev/s], impeller diameter `D`, tank diameter `T`, liquid depth `H`, gravity `g`, bubble diameter `d_b`, particle diameter `d_p`, molecular diffusivity `D_m`, surface tension `σ`, `Δρ = ρ_l − ρ_g`.

```
Re   = ρ_l N D² / μ_l                 (stirring Reynolds, representative velocity U_tip = πND)
Fr   = N² D / g
We   = ρ_l N² D³ / σ
Eo   = Δρ g d_b² / σ                  (=Bond, bubble scale)
Mo   = g μ_l⁴ Δρ / (ρ_l² σ³)
Ca   = μ_l U / σ
Sc   = ν_l / D_m
Pe_N   = N D² / D_m = Re·Sc          (impeller velocity scale ND — REV-CFD-MJ-006)
Pe_tip = U_tip D / D_m = π·Re·Sc     (tip speed U_tip = πND; each use site must
                                      state which Pe it means — no bare "Pe")
Da_n = k C_ref^{n-1} · (L/U)          (reaction of order n, k is the rate constant; noted separately per order)
St   = τ_p / τ_f,   τ_p = ρ_p d_p² / (18 μ_l)
Np   = P / (ρ_l N³ D⁵),  P = 2π N T_q = Ω T_q   (T_q=torque; N is in rev/s, ρ is on the liquid-phase basis)
N_Q  = Q / (N D³)                     (Q=net volumetric flow through the impeller discharge surface)
```

Lattice-side constraints: `Ma_lattice = U_tip/c_s ≤ 0.1`, Cahn number `Cn = W/L`, interface Péclet `Pe_φ = U W / M`, relaxation time `τ ∈ [τ_min, τ_max]`.

### 2.2 Matching priority (when simultaneous matching is impossible, #25)

The degrees of freedom of the physical→lattice conversion are finite. When all dimensionless numbers cannot be matched simultaneously, the priority is fixed:
**(1) Re → (2) density ratio, viscosity ratio + We/Eo (interface dynamics) → (3) Fr (when the free surface / buoyancy dominates) → (4) Sc/Pe, Da (scalar, reaction) → (5) St (particles)**.
The unit-conversion layer must run a feasibility check, and on `Ma>0.1` / `τ∉[τ_min,τ_max]` / excessive `Cn` / diffusion-number or CFL violation it warns with the compromised dimensionless numbers and the error made explicit.

---

## 3. Governing equation system (revised)

```
Continuous phase (recovered by low-Mach LBM, well-balanced gravity with phase-wise density.
rev.4 / REV-CFD-CR-002: mass-flux consistency with the phase-field diffusion —
with ρ=ρ(φ) and a diffusive phase flux J_φ, the naive ∂ρ/∂t+∇·(ρu)=0 cannot hold;
the density flux J_ρ = (ρ_l−ρ_g) J_φ must appear in BOTH the continuity identity
and the momentum advection (consistent/AGG-type formulation — mandatory at
ρ_l/ρ_g ≈ 10³)):
  ∂ρ/∂t + ∇·(ρu + J_ρ) = 0,        J_ρ = (ρ_l−ρ_g) J_φ
  ∂(ρu)/∂t + ∇·[(ρu + J_ρ) u] = -∇p + ∇·[ (μ(γ̇)+μ_t)(∇u+∇uᵀ) ]
                        + F_s + ρ g + F_b^{scalar} + F_g^{disp} + F_p + F_rot
  - The SAME discrete J_ρ is used in both equations (single code path — verified
    by code review and the advected-droplet conservation test, §8 VR-STR-03/05).
  - If a quasi-incompressible / pressure-evolution formulation is adopted instead,
    its continuity statement, divergence condition, and conservation checks must
    replace the above explicitly — silence is not an option.
  - Gravity imposes ρg on all phases and discretizes the hydrostatic ∇p_hydro = ρg in a well-balanced way (#34).
  - F_b^{scalar} (rev.2, active density feedback): the Boussinesq perturbation force of the
    solute buoyancy F_b = ρ_0 β_C (C−C_0) g. It is exactly 0 at C≡C_0, and is **not mixed**
    with the well-balanced hydrostatic cancellation of ρ(φ)g (composed as an independent
    force source. For details see docs/proposals/active-scalar-feedback.md).
  - F_rot is only for the MRF mode (§4.3). The composition of μ(γ̇) and μ_t follows the implicit consistency of §4.7.

Two-phase interface (conservative Allen–Cahn phase field, fixed to the Fakhari 2017 family, #8.
rev.4: written in explicit conservative-flux form so J_φ is a first-class object):
  ∂φ/∂t + ∇·(φu + J_φ) = 0,   J_φ = −M [ ∇φ − (4/W) φ(1−φ) n̂ ]
  n̂ = ∇φ / (|∇φ| + ε),  φ∈[0,1] (φ=1: liquid, φ=0: gas),  M[length²/time]
  Density interpolation: ρ(φ) = ρ_g + φ(ρ_l−ρ_g).
  Viscosity interpolation (rev.4 / REV-CFD-MJ-013 — default frozen):
    **harmonic in μ**:  1/μ(φ) = φ/μ_l + (1−φ)/μ_g
  Alternatives (linear-in-μ, linear-in-ν) are explicit config options, logged in
  run metadata, and are NOT covered by the default validation bands.

Surface tension (based on the chemical-potential form, with a convention branch when σ is variable. rev.2, #7):
  μ_φ = 4β φ(φ−1)(φ−1/2) − κ ∇²φ
  σ = √(2κβ)/6,   W = 4√(κ/(2β))    ← these ARE the definitions for the adopted
  double-well free energy (rev.4 / REV-CFD-MJ-013: the "coefficients are
  model-defined" hedge is removed; internal consistency of {σ, W, κ, β, μ_φ} was
  verified in the codex round-2 review; the (σ,W) ↔ (κ,β) inversion is unique)
  - When σ is constant (reference form): F_s = μ_φ ∇φ (the CSF-equivalent σκn̂δ_s goes to the validation items)
  - When active with σ depending on C_k / temperature: F_s = μ_φ∇φ is not used directly; it is
    unified into a well-balanced CSF / chemical-potential combined form that avoids double-counting
    with the Marangoni tangential force (docs/proposals/active-scalar-feedback.md conventions D1/D2).
    The coefficients are frozen after being derived to the (κ,β,W,σ) convention of this section
    (derivation required — mandatory before implementation). A degeneracy test that agrees with the
    σ-constant reference form under degeneration to ∇σ=0 is placed in §8.

Dispersed gas phase (point-bubble mode, Euler-Lagrange, #12):
  m_b dv_b/dt = F_buoy + F_drag(Tomiyama) + F_lift + F_addedmass + F_walllub + F_TD
  The BIT (bubble-induced turbulence) production term is added to the LES (§4.2, #46)

Dispersed particles (Euler-Lagrange, #16):
  m_p dv_p/dt = F_drag(Schiller-Naumann, Re_p range made explicit) + F_buoy
               + [in the high-accuracy case] F_Saffman + F_Basset + F_Faxen
  In the two/four-way case, reaction forces are scattered with a regularized kernel and momentum conservation is validated

Scalars / reactions (component k, active/passive made explicit, #13, #14.
rev.4 / REV-CFD-MJ-011: conservative forms are normative for two-phase and
active-density cases — the non-conservative single-phase form is a special case):
  Single-phase passive (simplified form, valid only when ρ, α uniform):
    ∂C_k/∂t + u·∇C_k = ∇·[ (D_k + ν_t/Sc_t) ∇C_k ] + R_k(C) + Ṡ_k^{if}
  Two-phase, phase-wise conservative (normative for gas–liquid scalars,
  q ∈ {gas, liquid}, α_liq = φ, α_gas = 1−φ):
    ∂(α_q C_{k,q})/∂t + ∇·(α_q u C_{k,q})
      = ∇·[ α_q (D_{k,q} + ν_t/Sc_t) ∇C_{k,q} ] + α_q R_{k,q}(C) + S_{k,q}^{if}
  Density-based active scalar (when the scalar feeds back into ρ):
    ∂(ρY_k)/∂t + ∇·(ρ u Y_k + J_k) = R_k + S_k^{if}
  - The SGS scalar flux is closed with the turbulent Schmidt number Sc_t (default Sc_t=0.7, variable).
  - S^{if}: for the resolved interface, normal jump + partition coefficient (Henry partition — the sign convention
    is: S_{k,liq}^{if} = −S_{k,gas}^{if}, interfacial flux positive into liquid);
    for point-bubble, k_L a(C*−C) (#35).
  - Conservation statement: Σ_q ∫ α_q C_{k,q} dV changes only by boundary fluxes
    and reactions — this is the quantity tested in VR-STR-05 scalar drift.

Eddy viscosity (SGS, Smagorinsky and WALE separated, #4):
  Smagorinsky: ν_t = (C_s Δ)² |S̄|,  |S̄|=√(2 S̄:S̄)
  WALE:        ν_t = (C_w Δ)² (S^d:S^d)^{3/2} / [ (S̄:S̄)^{5/2} + (S^d:S^d)^{5/4} ]
               S^d is the deviatoric symmetric part of the square of the velocity-gradient tensor (local gradient reconstruction is required)
```

---

## 4. Functional requirements (numerical methods, coefficient-frozen version)

### 4.1 Foundation LBM core
- **FR-CORE-01**: D3Q19/D3Q27 selectable. **The condition for D3Q27 to be the default is limited to "multiphase or strong forcing or cumulant use"** (#30). The Hermite order of the equilibrium distribution held on each lattice and the recovery accuracy are defined separately (D3Q19 is limited in third-order isotropy). **The M-F fidelity-default scenario falls under the multiphase / strong-forcing condition, so it is always D3Q27**. In single-phase / weak-forcing derived scenarios (e.g. VR-STR-01 single-phase stirring) D3Q19 is permitted (rev.2, finding 11).
- **FR-CORE-02**: Central moments (cascaded) / cumulant are implemented. Stability is stipulated not by "guarantee" but by **the allowable relaxation-rate range on the target bench, positivity, and the presence/absence of regularization/filtering/entropic limiter** (#31).
- **FR-CORE-03**: Guo forcing. The velocity moment is `ρu = Σ c_i f_i + Δt F/2`.
  Stress evaluation uses the forcing second-moment correction **as defined by the
  single equation in FR-STRESS-01** — prose words like "subtract"/"add" are banned
  from this topic; the equation is the only definition (rev.4 / REV-CFD-CR-004).
- **FR-CORE-04**: `Ma_lattice ≤ 0.1`, control of the compressibility error `O(Ma²)`. Include the consistency of acoustic scaling and incompressibility in the unit-conversion feasibility (#25).

### 4.2 Turbulence model (LES-LBM)
- **FR-LES-01**: Smagorinsky (including dynamic Germano) and WALE are implemented **as separate equations**. WALE is the default (`ν_t→0` near walls). **Because WALE requires the full velocity gradient, the "no finite differences" requirement is withdrawn**, and the local-gradient reconstruction method (moments or compact differences) is made explicit (#4).
- **FR-LES-02**: The eddy-viscosity relaxation-time reflection is **`τ_eff = 1/2 + (ν_0+ν_t)/(c_s²Δt)`** (general form). The lattice-unit `c_s²=1/3, Δt=1` simplified form `Δτ_t = 3ν_t` is noted separately (#3 fix).
- **FR-LES-03**: The wall-shear-dominated region uses a `y⁺` wall function or a wall-fitted interpolated boundary. `τ_eff` has not only the lower bound `>1/2` but also **upper clipping and diagnostics** (to avoid over-diffusion and boundary-accuracy degradation, #27).
- **FR-LES-04**: The SGS scalar flux (turbulent Schmidt number `Sc_t`) and the SGS heat flux (turbulent Prandtl number) are reflected into the ADE-LBM relaxation time (#14).

### 4.3 Rotating impeller (mode exclusion, #5, #6, #21, #22)
- **FR-ROT-01** (IBM-inertial): direct-forcing IBM (Uhlmann type). Target rigid-body velocity `U=Ω×r`. **"Galilean invariance guaranteed" is removed**, and **thresholds on slip velocity, torque error, and momentum-conservation error** are set for Taylor-Couette, rotating cylinder, and moving-wall Couette (including the adoption conditions for multi-direct-forcing / implicit IBM).
- **FR-ROT-02** (MRF-frozen-rotor): inside the rotating zone, **solve the relative velocity `u_rel = u_abs − Ω×r`**, and impose on the body force the Coriolis `−2ρ Ω×u_rel` and centrifugal `−ρ Ω×(Ω×r)`. **MRF is not applied to stationary tank walls or baffles**. The velocity-consistency condition at the rotating-zone boundary is defined. Cannot be started simultaneously with IBM moving blades.
- **FR-ROT-03**: Stationary walls / baffles = interpolated bounce-back (Bouzidi/Ginzburg), **moving blades = IBM or moving-wall interpolated BB**, clearly separated (#22). The update frequency of the STL distance field and the geometric error during rotation are defined.
- **FR-ROT-04**: `Np = P/(ρ_l N³ D⁵)`, `P = Ω T_q`, `N = Ω/(2π)` are fixed (2π double-counting prohibited). During gas aeration, **the ungassed `Np_0`, the gassed `Np_g`, and the gassed power-drop ratio** are output separately. `N_Q = Q/(ND³)`, and the integration surface, velocity components, time/phase averaging, and backflow handling of `Q` are defined (#23, #24).
- **FR-ROT-05** (sliding-overset, advanced): overset-grid halo interpolation compatible with MPI. For reference-tier validation.

### 4.4 High-density-ratio two-phase flow
- **FR-VOF-01**: Conservative Allen-Cahn (fixed in §3). The mass-conservation error is stipulated per bench — for a **closed static droplet / rising single bubble / sparger open boundary**, tolerances including time, grid resolution, and outflow/inflow amounts are set (#9). Shan-Chen is not adopted for this use.
- **FR-VOF-02** (rev.4 / REV-CFD-MJ-005 — dimensional fix): spurious currents on a
  static droplet are bounded by the (dimensionless) capillary number
  `Ca_spurious = μ_l |u|_spurious / σ < 10⁻³` (target We→0, resolution stated).
  The old `|u|·L/(σ/μ)` form carried a stray length dimension and is void. A
  length-bearing indicator, if wanted, is `Re_spurious = |u|_spurious L/ν_l` —
  a separate metric, never called Ca. well-balanced chemical-potential form (implements the coefficient relations of #7).
- **FR-VOF-03** (sparger, rev.4 / REV-CFD-CR-001 — **phase-inversion fix**):
  the sparger injects GAS; under the §3 definition (φ=1: liquid, φ=0: gas) the
  injected phase value is **φ=0**. The rev.1 text banned "plain `φ=1` + velocity
  Dirichlet" — that read as a liquid-injection ban and inverted the phase BC.
  Corrected requirements:
  - Choose from gas-phase volumetric-flow boundary / stochastic bubble injection / resolved orifice. A plain
    **`φ=0` + velocity Dirichlet alone is banned** — the injection model must
    simultaneously satisfy gas volumetric-flow conservation, pressure consistency,
    contact angle, and the `d_b/W`, `d_b/Δx` lower bounds (#10).
  - **The scenario schema/API never exposes raw φ for inlets**: config says
    `inlet_phase: gas | liquid` and the core maps it (gas→φ=0, liquid→φ=1) —
    enforced by config validation (A-4 style). Outputs report `φ_liquid` and
    `α_g = 1−φ` with explicit names.
  - Acceptance: gas-inlet setting injects φ=0 (unit test); a sparger-only case
    balances injected gas volume vs. domain gas-volume increase within tolerance
    (VR-STR-02c precursor); no schema field accepts a raw φ boundary value.
  Breakup / coalescence is weakened to "numerically tolerated", and it is stated explicitly that real thin-film drainage is not resolved (#11).
- **FR-VOF-04** (point-bubble): the switching condition includes `d_b/W, Eo, Re_b, α_g, We_b, mass-transfer consistency`. Interphase mass / momentum / scalar conservation laws in the hybrid-mixed case are defined (#12).
  **(rev.3, P1)** Population balance modelling (PBM) of the bubble-size distribution is
  required on the point-bubble path (breakup/coalescence kernels, e.g. Luo–Svendsen /
  Prince–Blanch): a mono-disperse point-bubble model cannot support the `d_32`
  acceptance of VR-STR-02 (internal consistency). Per-bubble gas-phase composition
  bookkeeping (component inventory and interfacial transfer budgets) must reconcile
  with FR-VOF-05. *Scope alignment (rev.2/§0)*: point-bubble is a relaxation extension
  (API-reserved in v1); this PBM requirement binds when that extension is implemented —
  in the resolved-phasefield default, `d_32` is measured from the resolved interface.
- **FR-VOF-05**: Interfacial mass transfer is **separated into the resolved interface (normal flux / partition coefficient / phase-wise diffusion) and point-bubble (`k_L a(C*−C)`)** (#35). The applicability ranges of Henry's law and the Sherwood number are made explicit.

### 4.5 Dispersed particles
- **FR-PART-01**: one/two/four-way switching by `α_p`/mass-loading (thresholds made explicit). The `Re_p` applicability range of Schiller-Naumann, the reaction-force scatter kernel, and momentum-conservation validation are required (#16). For neutrally buoyant fine particles, the need for Saffman/Basset/Faxen is decided by `d_p/Δx` and `St`.
- **FR-PART-02**: Switchable to a resolved-particle method (PSM/Noble-Torczynski, Ladd/Aidun-Lu).
- **FR-PART-03**: Record `∫γ̇dt` and `max γ̇` along trajectories. **When tracking under LES, enable SGS turbulent dispersion (stochastic dispersion)**, or state resolved-only explicitly (to avoid grid dependence of the exposure PDF/CDF, #17).
- **FR-PART-04 (rev.4 / REV-CFD-MJ-012 — four-way contact contract; Phase-2
  extension, API-reserved in v1)**: four-way coupling requires a soft-sphere
  normal-collision model with explicit parameters: restitution `e_n`, collision
  time `T_col`, spring `k`, dashpot `η`, max overlap `δ_max`, particle substep
  `Δt_p` (with `Δt_p ≲ T_col/10`).
- **FR-PART-05 (rev.4)**: when `d_p/Δx` does not resolve the lubrication gap, a
  lubrication correction (or calibrated implicit lubrication) is required; the
  applicability condition is stated with the model.
- **FR-PART-06 (rev.4 — config guard, initial delivery)**: while four-way is
  unimplemented/unvalidated, runs exceeding the `α_p` / mass-loading threshold of
  the two-way regime are **rejected at config validation** (A-4 style), with the
  threshold and its source stated in the error message. Initial delivery ships
  two-way + this guard; contact benches (particle–particle, particle–wall,
  settling, sheared suspension, overlap ≤ δ_max) gate the Phase-2 extension.

### 4.6 Stress-field evaluation (convention fixed, #1, #2, #18, #19, #20)
- **FR-STRESS-01** (rev.4 / REV-CFD-CR-003, CR-004 — stage convention and forcing
  correction fixed by equations, not prose): strain rate is evaluated locally from
  non-equilibrium distributions. **The default stage is pre-collision /
  post-streaming** (the distribution as it arrives, before collide — the stage the
  standard coefficient below is derived for):
  ```
  f_i^{neq,pre} = f_i^{pre} − f_i^{eq}(ρ, u)        (u includes the F/2 correction)
  Π_neq_raw     = Σ_i c_iα c_iβ f_i^{neq,pre}
  Π_force       = −(Δt/2)(u_α F_β + u_β F_α)         (Guo forcing second moment,
                                                      for THIS engine's u/f_eq defs)
  Π_neq_corr    = Π_neq_raw − Π_force  =  Π_neq_raw + (Δt/2)(uF + Fu)
  S_αβ          = − Π_neq_corr / (2 ρ c_s² τ_eff Δt)
  ```
  `Π_neq_corr` is the ONLY normative definition; natural-language sign words are
  non-normative. The exact sign of Π_force is derivation-frozen against this
  engine's Guo discretisation **before implementation** and locked by a negative
  test (body-force Poiseuille must FAIL with the sign flipped — §8 VR-STR-03).
  **If the post-collision / pre-streaming stage is used instead** (e.g. inside a
  fused kernel where it is cheaper), the stage transform is mandatory:
  BGK: `Π_neq,pre = Π_neq,post / (1 − 1/τ_eff)`; MRT/cumulant: apply the inverse
  shear-moment relaxation `R(τ_shear)⁻¹` — then proceed with the equations above.
  The stress-evaluation API takes a required `neq_stage` enum
  (`PreCollision | PostCollision`) — no default-by-silence, misuse is a compile-
  or construct-time error (same philosophy as A-4/A-5 guards).
  For cumulant/MRT, the coefficient is corrected by the shear-moment relaxation rate. The circular dependence of the Smagorinsky closure
  is solved in **algebraic closed form** (`τ_eff` obtained explicitly from `|Q|`; Hou et al.-type quadratic).
- **FR-STRESS-02**: The output stress is defined separately into **`resolved viscous` / `SGS` / `capillary` / `particle`**. For `γ̇=√(2S:S)`, the second invariant `II_S`, and von Mises, the source tensor is restricted (#19).
- **FR-STRESS-03**: Wall shear is defined per mode (**tangential velocity-gradient reconstruction / IBM forcing integral / MEM**). The handling when the non-equilibrium quantity near the interpolated boundary does not represent the wall gradient is made explicit. Validation includes a curved moving wall (#20).
- **FR-STRESS-04**: The composition rule, iteration procedure, convergence criterion, `τ_min/τ_max`, and LES applicability range of non-Newtonian `μ(γ̇)` (Carreau-Yasuda/Casson/power-law) and `μ_t` are made explicit (to avoid double-counting and divergence, #18).

### 4.7 Boundaries, gravity, initialization (new, #33, #34, #45, #47)
- **FR-BC-01** (top boundary, mandatory specification): choose from `closed` / `free-surface` / `degassing-outlet`. During sparging, a gas-phase exhaust outlet is required (a closed tank + gas-phase inflow only is unphysical due to gas accumulation, #33). Headspace pressure, free-surface deformation, and free-surface contact angle are defined.
- **FR-BC-02** (gravity): `ρg` on all phases, dynamic-pressure/hydrostatic decomposition, and a well-balanced hydrostatic test (`|u|<ε` in static stratification) are required (#34).
- **FR-BC-03** (wettability): per-wall contact-angle boundary condition, slip/no-slip, and phase-field flux condition are defined (#47).
- **FR-BC-04** (scalar wall): choose from no-flux/adsorption/reactive wall (#35).
- **FR-INIT-01**: Initial velocity/pressure/phase-field/scalar/particle placement, impeller ramp-up, gas-flow-rate ramp, statistics-sampling start time, and quasi-steady decision criterion are required (#45).

### 4.8 Extension & closure contracts (rev.3, P2)

- **FR-EXT-01**: Define explicit contracts for the trait/strategy extension points of
  §1 and for user-supplied closures — reaction rates `R_k`, non-Newtonian viscosity
  `μ(γ̇)`, body-force sources, and the relaxation-mode implementations
  (MRF / point-bubble / one-way / AMR):
  - input/output signatures with explicit physical vs. lattice units;
  - determinism (identical inputs → bit-identical outputs);
  - GPU evaluability (state-free, portable to wgpu);
  - error handling (NaN/divergence detection at the contract boundary);
  - schema versioning and backward compatibility.
  The primary boundary is Rust traits; foreign-language ABI/SDK is deferred to a
  separate API specification (see §10 product-layer note). The fidelity-default
  implementation is the default of each trait; relaxation implementations swap in
  under the same contract and are accepted via VR-STR-RELAX.
  *Implementation note*: this contract work is co-designed with the R-Phase 2 / B-1
  trait-boundary design (SOLVER_IMPROVEMENT_SPEC WP-B) — one design, two consumers.

---

## 5. Coupling / time integration (#28)

- **FR-COUP-01** (rev.4 / REV-CFD-MJ-007 — the dataflow is split by scalar mode so
  "active" is not silently one step lagged):
  **passive scalar**: phase-field update → ρ/μ field update → force-source composition (`F_s+ρg+F_g+F_p+F_rot`) →
  fused collide-stream-moments → boundary → scalar ADE → reaction (split) → particle integration.
  **active scalar (fidelity default — predictor–corrector)**:
  scalar/reaction predictor → property update `ρ(C), μ(C), σ(C)[, T]` →
  force-source composition (incl. `F_b^{scalar}`, Marangoni) → flow step → scalar ADE corrector →
  reaction corrector → property re-evaluation (→ optional flow–scalar iteration
  for stiff coupling). Time-lagged explicit feedback is allowed only as the
  flagged relaxation `active_scalar_lagged=true`, with stated applicability
  (weak feedback, non-stiff), stability conditions, and a lag-error benchmark —
  accepted via VR-STR-RELAX. Mode (coupled/lagged) is logged in run metadata.
  **For strong coupling, stiff reactions, and surface-tension waves, operator-splitting error, subcycling, and iterative strong coupling are required**.
  The respective constraints of the capillary time step `Δt_σ ≤ √(ρ̄ Δx³/(2πσ))`, the particle `Δt_p`, and the reaction ODE `Δt_r` are imposed.
  Acceptance: on the active-scalar standard bench (Marangoni or
  concentration-dependent viscosity), the feedback error converges under
  time-step halving (§8 VR-STR-06+/RELAX).
- **FR-COUP-02**: The reaction solver switches among explicit/implicit/Rosenbrock-BDF by stiffness detection. **Acceptance criteria for negative-concentration limiting, element-conservation error, and split error** are defined (#15).
- **FR-COUP-03**: Dimensionless matching follows the priority of §2.2 + feasibility check (#25).
- **FR-COUP-04**: `probe_state_hash` bit-equivalence is **limited to implementation regression on a single backend**. Physical validity and conservation laws are on separate criteria (§8, #28, #42).
- **FR-COUP-05**: AMR is an advanced option. When enabled, coarse-fine conservative interpolation, the time-step ratio, and dedicated validation are required (#29).

---

## 6. I/O / visualization

- **FR-IO-01**: 3D field output is **uniform grid = VTI, structured-curved = VTS, unstructured/AMR = VTU/AMR** (#43). `φ` is a diffuse-interface indicator, not the void fraction (#36).
  **ε_g processing definitions (rev.4 / REV-CFD-MN-014)** — every ε_g output
  carries filter width, averaging volume, and time window as metadata:
  - resolved-phasefield: `ε_g_raw = ⟨1−φ⟩_V` and
    `ε_g_thresholded(φ_c) = volume(φ<φ_c)/V`, default `φ_c = 0.5` — both output.
  - point-bubble: `ε_g_bubble = Σ_b V_b W_kernel(x−x_b) / V_filter`
    (kernel-smoothed void fraction).
  - hybrid: `ε_g_total = ε_g_resolved + ε_g_bubble` with double-count exclusion
    over the resolved region.
  Any ε_g indicator must be recomputable from a snapshot; experiment comparisons
  state which definition was used.
- **FR-IO-02**: Time-averaged / phase-averaged statistics (mean field, RMS, Reynolds stress). **Phase averaging is only for IBM/overset unsteady modes**. MRF is output separately as a rotating-frame average / quasi-steady (#37).
- **FR-IO-03**: 3D display in the Web GUI (slices, isosurfaces, shear heatmap, time-series probes). The existing 2D canvas is extended to WebGL/WebGPU.
- **FR-IO-04**: Histogram/CDF of particle cumulative shear exposure (presence/absence of SGS dispersion made explicit).
- **FR-IO-05 (rev.3, P3 — mixing metrics)**: Derived outputs for **blend time**
  (time until the coefficient of variation CoV of a tracer falls below a stated
  threshold) and **RTD** (tracer response `E(t)`, mean residence time, variance).
  The homogenisation threshold and the tracer injection/detection surfaces must be
  explicitly defined per scenario.
- **FR-IO-06 (rev.3, P3 — large-scale I/O & resilience)**: Full-field dumps are
  impractical at target scales (§7 budget); require **parallel I/O**
  (HDF5/ADIOS2-class) + compression + in-situ statistics / downsampling.
  **Deterministic checkpoint/restart with crash recovery** (bit-reproducible resume
  including RNG state, particle state, and statistics accumulators) is mandatory.
  Formats are sized against the §7 budget table.
  *Convergence note*: builds on SOLVER_IMPROVEMENT_SPEC B-5 (snapshot API),
  C-3 (per-rank parallel I/O), C-8 (distributed checkpoint) — reuse, don't duplicate.

---

## 7. Non-functional requirements

- **NFR-01 (scale / memory budget, #44)**: `O(10⁸–10⁹)` cells. **A budget table for bytes-per-cell, number of distributions (D3Q27 × phase field × scalars), particle count, GPU memory, I/O volume, and checkpoint frequency** is mandatory. Since 1e9 cells × many distributions is on the order of **0.6 TB** at the fidelity default, and on the order of **TB to several TB** with all f64, multiple scalars, simultaneous checkpoint retention, and I/O buffers included (rev.2 fix), an estimate for wgpu multi-GPU + MPI partitioning is attached.

  **Budget table (rev.1b, fidelity-default configuration, deviation storage, ping-pong ×2. Per cell)**:

  | Component | Lattice/type | bytes/cell |
  |---|---|---|
  | Fluid distribution f | D3Q27 × 2 × f32 | 216 |
  | Phase-field distribution g (conservative Allen-Cahn) | D3Q19 × 2 × f32 | 152 |
  | Scalar distribution h (per component) | D3Q7 × 2 × f32 | 56 |
  | moments / property fields (ρ, u×3, φ, μ_φ, ∇φ×3, ν_t, γ̇, τ_eff) | 12 × f32 | 48 |
  | Mask / flags | u8×2 | 2 |
  | Statistics accumulators (mean u×3, RMS×3, Reynolds stress 6, etc.) | ~13 × f32 to f64 | 52–104 |
  | Interface-band f64 promotion (band width ~2W, 5–10% of all cells, +368 B/band-cell of f+g amortized) | amortized | +18–37 |
  | Additional overhead when including the interface-band curvature / reduction workspace (rev.2 recheck) | amortized | up to ~+40 |
  | **Total (1 scalar component)** | | **≈ 540–620 B/cell** |

  Conversion: **1e8 cells ≈ 56–62 GB** (feasible on a single node with this machine's M5 Max 128 GB, upper limit ~1.5e8) /
  **1e9 cells ≈ 0.56–0.62 TB** (f32 bulk), **≈ 1.1–1.2 TB** at the all-f64 reference tier.
  Particles 10⁷ × ~100 B = 1 GB (negligible). A checkpoint, as raw storage of the distributions,
  is the same amount as the field data per dump (~0.5 TB/dump at 1e9) → the frequency is back-calculated from the I/O bandwidth,
  defaulting to 2–5 dumps per job. GPU: 8–16 GB/card → 1.3–2.6e7 cells/card in the f32 configuration,
  and 1e9 cells **require 40–80 cards of multi-GPU or CPU-cluster MPI** (not feasible on a single GPU).
  Conclusion: 1e9 at the fidelity default is cluster-only. Development / validation use ≤256³ (1.7e7 cells ≈ 10 GB) as
  the standard, and scale measurements are integrated into the R3 cluster plan (CLUSTER_OPTIONS.md).
- **NFR-02 (precision policy, #32, rev.2 vocabulary cleanup, rev.4 / REV-CFD-MJ-009 —
  enumerated so array design / GPU kernels / memory budget can bind to it)**:
  `precision_profile ∈ { full_f64, mixed_safe (default), mixed_fast }`:
  - **full_f64** (reference tier): all distributions, phase field, scalars,
    particle statistics, reductions in f64. High-density-ratio reference
    validations also run here.
  - **mixed_safe** (fidelity default = §1 profile): bulk distributions f32;
    **f64 fixed for**: `φ, ∇φ, κ(curvature), μ_φ, F_s, ρ(φ), μ(φ)`,
    distributions inside the interface band, all global reductions, torque,
    `Np`, `N_Q`, mass/volume counters, particle cumulative exposure.
    **interface_band = max(3W, 6Δx)** — provisional default; the band width is
    re-frozen by the W-VOF characterization (§10) and recorded in PHYSICS.md.
  - **mixed_fast** (relaxation extension): single-phase / weak-coupling only;
    permitted only when density ratio ≤ stated limit AND the Ca_spurious and
    mass-drift validations pass; config validation rejects out-of-range use.
    Accepted via VR-STR-RELAX-f32.
  Each profile has an array-type table and memory-budget column (§7).
  Consistent with the `ρ_l/ρ_g≈10³`, `Ca_spurious<10⁻³`, and mass-conservation requirements.
- **NFR-03 (performance)**: Integrate phase field, scalars, and forcing into the fused `step_band`, and maintain the 3D extension of ring double-buffering and SoA plane-major.
- **NFR-04 (determinism)**: Reductions in deterministic order. GPU/MPI use tolerance-based regression (bit-equivalence limited to a single backend, #42).

---

## 8. Validation / acceptance criteria (quantified, already wired as VALIDATION.md **T17**. With thresholds, #38–#42)

The validation tests are authored adversarially from this specification by codex/Opus and separated from the implementation.

**Band governance (rev.4 / REV-CFD-MJ-010 — reconciling "numbers now" with the
experiment-driven freeze protocol)**: every VR-STR item carries a **provisional
numeric band from day one** (table below — these are the MVP gate). Bands are
finalized by the established protocol (implement → characterize → record rationale
in PHYSICS.md → freeze in VALIDATION.md T17) under one asymmetric rule:
**tightening a band is always allowed; loosening a provisional band requires a
recorded physical rationale in PHYSICS.md** (reference uncertainty, method order,
resolution limit — as exercised for T15.5). This removes both failure modes:
un-testable placeholder specs AND post-hoc self-serving thresholds.
Each test is specified with: metric / target·reference / tolerance / resolution /
time window / backend / pass-fail rule (the T17 row format).

**Provisional bands (MVP gate; supersede the "±tolerance%" placeholders)**:
- Rushton `Np` vs experimental correlation: **±10%**
- PIV/LDA velocity profiles (VR-STR-01): **L2_rel < 15%, L∞_rel < 30%** per line
- static droplet mass drift: **< 0.1% / 1000 steps**; advected droplet (one period,
  periodic box): total mass drift **< 0.1%** (CR-002 acceptance)
- single-bubble terminal velocity vs Grace (02a): **±10%**
- `k_L a` vs correlation (02c): **±25%**
- well-balanced static stratification (VR-STR-06): **max|u| < 10⁻⁶ (lattice units)**
  at ρ ratio 10³ — provisional; retighten after discretisation freeze
- GPU/MPI cross-backend drift (VR-STR-05): mean quantities **< 2%**, higher-order
  statistics **< 5%** (bit-equality stays single-backend-only)
- `Ca_spurious < 10⁻³` (already fixed, dimensionally corrected FR-VOF-02)

**Mandatory negative/consistency tests added by rev.4**:
- forcing-moment sign negative test: body-force Poiseuille FAILS with Π_force sign
  flipped (CR-004);
- stress stage-convention cross-check: pre-collision evaluation vs post-collision+
  transform agree within tolerance on Couette/Poiseuille/Taylor–Couette (CR-003);
- J_ρ consistency code-path check + droplet advection conservation (CR-002);
- sparger phase unit test: gas inlet injects φ=0, gas-volume balance closes
  (CR-001);
- scalar total-mass conservation in the phase-wise form (MJ-011);
- active-scalar dt-halving convergence (MJ-007).

- **VR-STR-01 (single-phase stirring)**: standard baffled tank (`D/T`, `C/T`, blade geometry, number of baffles fixed), specified Re range, ungassed. Rushton `Np` = experimental correlation ±tolerance%; the impeller discharge-velocity profile is compared against threshold `L2/L∞rel` on PIV/LDA reference lines (#38).
  **(rev.3, P4) Reference datasets**: Wu & Patterson (1989) LDA; Deen et al. (2002)
  PIV (standard Rushton, D/T=1/3, 4 baffles); standard `Np` correlations. Numeric
  bands are frozen via the T17 experiment-driven protocol — not hardcoded here.
- **VR-STR-02 (gas-liquid, split into 02a/b/c in rev.2)**: **02a single bubble** = `U_t` compared by relative error against the Grace diagram Eo-Mo-Re. **02b bubble swarm** = `ε_g` spatial distribution, swarm rise velocity (hindered rise), `d_32` when breakup/coalescence is allowed, and turbulence intensity when BIT is used (`ν_t` response). **02c stirred-tank aeration** = experimental-correlation comparison of `ε_g, d_32, k_L a` (#39).
  **(rev.3, P4) References**: single bubble = Grace diagram (Eo-Mo-Re); aerated tank =
  published `ε_g`/`d_32`/`k_L a` data and correlations. In point-bubble / RELAX-PB
  evaluations, `d_32` presupposes the FR-VOF-04 population balance (P1); in the
  resolved-phasefield default it is measured by interface segmentation.
- **VR-STR-03 (shear / stress)**: separate the method of manufactured solutions (MMS) single-phase, curved Couette, rotating cylinder, non-Newtonian Poiseuille, and multiphase static droplet, and set the **grid convergence order** and `L2/L∞`. Line design accounting for the divergent severity of `L∞` near walls (#40).
- **VR-STR-04 (scalar / reaction)**: Taylor-Aris dispersion, a reaction-diffusion front with known `Da`, `k_L a` (the computation formula = interface integral or correlation, made explicit). The tolerance, target `Pe/Da/Sc`, and boundary conditions of each are specified (#41).
- **VR-STR-05 (coupled regression / conservation)**: `probe_state_hash` is limited to single-backend regression. **Drift thresholds for mass, momentum, scalar totals, gas-phase volume, particle count, and energy-like quantities are set individually**. Energy-like quantities (kinetic energy, interfacial free energy, particle kinetic energy) are treated **as monitoring quantities for unphysical drift, not as exactly conserved** (rev.2). GPU/MPI are tolerance-based (#42).
- **VR-STR-06 (well-balanced)**: `|u|<ε` in static stratification (#34). **06+ (rev.2)**: with the active scalar ON and `C≡C_0`, the same quiescence is satisfied (exact-zero degeneration of `F_b^{scalar}`). The `∇σ=0` degeneration of the variable-σ form (agreement with the σ-constant reference form) is also placed in this group.
- **VR-STR-07 (initialization independence)**: varying the spin-up / statistics-start conditions, the quasi-steady statistics agree within threshold (#45).
- **VR-STR-RELAX (relaxation-mode equivalence, newly added in rev.2 — finding 1)**: relaxation extensions are accepted by relative degradation against the **corresponding fidelity reference solution**. The comparison target and measured quantity of each axis (tolerances via characterize→freeze):
  - **RELAX-MRF**: against the IBM-inertial (or sliding-overset) reference at the same geometry and same Re, the tolerance of `Np`, the impeller discharge-velocity line, the mean velocity field, and torque. Application is limited to configurations where the steady approximation holds.
  - **RELAX-PB (point-bubble)**: against the resolved-phasefield reference, the tolerance of `ε_g`, `d_32`, `k_L a`, and the momentum/scalar budget. The applicability range is limited by `d_b/Δx, d_b/W, α_g` (same as the switching condition of FR-VOF-04).
  - **RELAX-1W (one-way)**: against the two-way reference, the tolerance of particle statistics and the mass-loading upper limit for which neglecting the reaction force is permitted.
  - **RELAX-AMR**: against the uniform reference, conserved quantities, interface position, torque, velocity-field norm, and the budget error when crossing coarse-fine boundaries.
  - **RELAX-f32 (aggressive f32)**: against the fidelity-profile (or all-f64) reference, the tolerance of conserved-quantity drift, `Ca_spurious`, `Np`, interface curvature, and reduction quantities. Application is limited to single-phase / weak coupling (NFR-02).

---

## 9. Major technical risks

| # | Difficulty | Risk | Mitigation |
|---|---|---|---|
| 1 | High-density-ratio two-phase (`10³`) | Interface instability, spurious currents, f32 rounding | well-balanced phase field + D3Q27 + f64 reduction, point-bubble alternative (§1) |
| 2 | High-Re stability | Divergence, hyperviscosity, positivity | cumulant + WALE, algebraic closure `τ_eff`, limiter (§4.1,4.6) |
| 3 | Rotating-boundary conservation | IBM slip, torque error | multi-direct-forcing, overset reference validation, thresholding (§4.3) |
| 4 | Coupling stiffness | Time-scale divergence of reaction / interface / rotation | operator-splitting-error evaluation, subcycling, capillary dt (§5) |
| 5 | Computational cost / memory | 1e9 cells × many distributions = TB scale | memory budget table, multi-GPU + MPI, AMR advanced option (§7) |
| 6 | Mass transfer / BIT | Turbulence of unresolved aeration, insufficient transfer | Sc_t SGS, BIT production term, resolved/point separation (§4.2,4.4) |

---

## 10. Design decisions (finalized items and remaining implementation details)

**Finalized (rev.1a, PM decision)**: the default of every axis is set to **fidelity-first**, and low-cost approximations are implemented as bolt-on extension points (§1 design principle). The four previously unresolved items are resolved as follows:

- Interface default = `resolved-phasefield` (interface / mass-transfer fidelity first). `point-bubble` is a relaxation extension.
- Scalar default = `active` (with property feedback). `passive` is a relaxation extension.
- Grid default = `uniform` (fully resolved). `block-AMR` is a relaxation extension.
- Precision default = fidelity profile (interface neighborhood, conserved quantities, reductions in f64; far-field bulk in f32). All-f64 is the reference tier, aggressive f32 is a relaxation extension.

**Remaining implementation details (spec refinement, not decisions)** — status as of rev.3:

- Concrete equations and stabilization (including Marangoni) of the feedback targets (σ, viscosity, density, [temperature]) of the `active` scalar.
  → researched: docs/proposals/active-scalar-feedback.md. **One derivation is mandatory
  before implementation** (Marangoni coefficient consistency with the (κ,β) convention,
  §3). Thermal axis recommended as API-reserved extension.
- Tolerance thresholds against the fidelity reference solution for each relaxation extension (appended to §8 VR).
  → structure defined as VR-STR-RELAX (rev.2); numeric bands frozen at relaxation
  implementation time.
- The f64/f32 boundary of the fidelity profile (band width near the interface, reduction range).
  → frozen experimentally during W-VOF implementation (characterize→freeze).
- API definition of the trait boundary (strategy swap point) of each mode axis.
  → contract requirements fixed as FR-EXT-01 (§4.8); concrete Rust API co-designed
  with R-Phase 2 / B-1.

**Product-layer scope note (rev.3, P5)**: GUI/CAD & STL import, materials DB,
Python/CLI SDK, parameter sweep & optimizer, cloud/cluster/queue integration,
packaged validation assets, and competitive benchmark tables are **out of scope for
this solver specification**. They are version-managed in separate volumes
(Product Requirements / API Specification / Validation Pack / Performance Benchmark).

---

## 11. Implementation dependency graph (rev.3 — priority + dependency DAG, not stage gates)

Items with no dependency edge between them are implemented **concurrently**
(parallel-agent worktrees, per the standing parallelization directive).
Mapping to the PLAN.md M-F delegation tracks is noted per row
(MF-α…ζ are the delegation bundles; W-items are the fine-grained DAG nodes).

| Item | Hard deps (must precede) | Parallel | Notes / PLAN track |
|---|---|---|---|
| W0 core basis (D3Q19/27, cumulant, Guo forcing) | — (strengthens M-C 3D basis) | — | prerequisite for all; = MF-α |
| W-EXT trait contracts (FR-EXT-01) | W0 | yes | early definition = prerequisite of all relaxation modes; low cost, high leverage; co-designed with R-Phase 2 B-1 |
| W-UNIT unit/nondimensional feasibility (§2.2) | W0 | yes | independent, early |
| W-STRESS stress fields (FR-STRESS) | W0 | yes | top priority (primary output + prerequisite of LES & particle exposure); ⊂ MF-β |
| W-ROT rotating IBM (FR-ROT-01) | W0 | yes | prerequisite of Np/N_Q; MRF/overset live behind W-EXT as relaxation/reference tiers; = MF-δ |
| W-GRAV well-balanced gravity (FR-BC-02) | W0 | yes | prerequisite of the interface track; ⊂ MF-γ |
| W-SCAL passive scalar ADE (§3 scalar eq.; SGS flux part waits on W-LES) | W0 | yes | ⊂ MF-ε |
| W-LES turbulence SGS (FR-LES) | W-STRESS | conditional | \|S\| closure needs the stress evaluation; ⊂ MF-β |
| W-VOF resolved interface (FR-VOF-01/02) | W-GRAV | conditional | fidelity default; hardest item; **critical path**; ⊂ MF-γ |
| W-PART particles + cumulative exposure (FR-PART) | W-STRESS (SGS dispersion: W-LES) | conditional | exposure integral needs the γ̇ field; ⊂ MF-ε |
| W-REACT reaction / active feedback (§3, FR-COUP-02; active feedback needs W-VOF) | W-SCAL | conditional | ⊂ MF-ε |
| W-BUB point bubbles + PBM + interfacial transfer (FR-VOF-03/04/05) | W0, W-SCAL, W-EXT | conditional | relaxation extension (API-reserved in v1, per §0) |
| W-BCTOP top boundary / degassing / contact angle (FR-BC-01/03) | W-VOF | conditional | ⊂ MF-γ |
| W-COUP coupling loop (FR-COUP) | active subsystem set | incremental | grows as tracks land; ⊂ MF-ζ |
| W-IO I/O & analysis (FR-IO incl. -05/-06) | each producing subsystem | incremental | Np←ROT, blend/RTD←SCAL, exposure←PART; ⊂ MF-ζ |
| W-VAL validation T17 (VR-STR-01–07, RELAX) | each subsystem | yes | codex adversarial authorship, separated from implementation |

**Parallel waves** (sets that start together):
1. After W0, mutually independent: **W-EXT / W-UNIT / W-STRESS / W-ROT / W-GRAV / W-SCAL** (6-way parallel).
2. After their deps: W-LES (←STRESS) / W-VOF (←GRAV) / W-PART (←STRESS) / W-REACT (←SCAL).
3. Later: W-BCTOP (←VOF) / W-BUB (←SCAL,EXT) / active feedback & interfacial transfer (←VOF).
4. Cross-cutting throughout: W-COUP / W-IO / W-VAL.

**Critical paths** (staff first):
`W0 → W-GRAV → W-VOF → W-BCTOP/interfacial transfer` (interface chain — longest, hardest) and
`W0 → W-STRESS → W-LES → W-PART` (stress/exposure chain).

*Boundary decisions upheld (rev.3)*: throughput/scaling KPIs stay delegated to
CLUSTER_OPTIONS.md (R3) — not duplicated here; no hardcoded numeric thresholds
(P4 adds dataset names only — bands freeze via the T17 protocol); the product
ecosystem (P5 list) lives in separate volumes.
