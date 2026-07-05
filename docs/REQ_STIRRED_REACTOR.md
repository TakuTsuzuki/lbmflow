# Requirements Specification (finalized): rotating-boundary / high-density-ratio two-phase / LES-coupled 3D multiphysics LBM solver

**Document ID**: REQ-M-F-STR / **Revision**: rev.3 (rev.1 = all 48 items from the codex adversarial review reflected / rev.1a = PM decision "default is fidelity-first, relaxation is an add-on extension point" reflected / rev.1b = PM integration: domain-neutralization of the title, follow-up to the core rename, addition of the §7 memory budget table, VALIDATION T17 wiring / rev.2 = all 11 items of the codex second review adopted: new relaxation-equivalence validation VR-STR-RELAX, clarification of scope wording, variable-σ surface-tension convention, F_b^scalar added, budget-table arithmetic fixes, and others — see docs/proposals/req-round2-findings.md / **rev.3** = competitive-review triage diff merged (authored as "rev.1c" against rev.1b, layered here onto rev.2): P1 population balance, P2 §4.8 extension contracts, P3 FR-IO-05/06, P4 reference datasets, P5 product-layer scope note, §11 implementation dependency DAG. New content in English per the 2026-07-05 language directive; full translation of this document is delegated to the translation session.)
**Positioning**: M-F (vertical feature) of `docs/PLAN.md` / upper-level requirements for `ARCHITECTURE_V2.md`.
Validation acceptance is [VALIDATION.md](VALIDATION.md) **T17** (wiring VR-STR-01 to 07).
**Target core**: `lbm-core` (formerly lbm-core2. D3Q19/D3Q27, CpuScalar/CpuSimd/wgpu, MPI partitioning)
**Representative application problem**: stirred-tank reactor (functional requirements are defined domain-neutrally, and the validation benches of §2 and §8
concretize this application. A region that directly overlaps the core use cases of M-Star CFD)

## 0. Review-reflection policy

- All Critical/Major formula, sign, and definition bugs (non-equilibrium stress evaluation, `τ_eff`, MRF apparent forces, surface tension, conservative Allen-Cahn, Np/N_Q, dimensionless numbers) were fixed and their coefficients frozen.
- Missing physics (top degassing, hydrostatic well-balanced, turbulent Schmidt-number SGS scalar, particle SGS dispersion, bubble-induced turbulence BIT, contact angle, initialization/spin-up, memory budget) was added as new requirements.
- **Delivery scope (clarified in rev.2)**: the **fidelity-default subsystem group (the "default" column of the §1 table) is implemented all at once, concurrently**. Relaxation extensions (MRF, point-bubble, one-way, block-AMR, aggressive f32) in the first version **only reserve the trait boundary, configuration schema, and validation items (VR-STR-RELAX)**; implementation is an add-on. The codex "mode split" is not a delivery phase but is mapped into the §1 configuration matrix as a **runtime-mode mutual-exclusion constraint**. Modes that physically conflict within the same computation (MRF + IBM in the same zone, phase-averaging + MRF, etc.) are rejected by config validation.
- Non-verifiable wording ("guaranteed," "naturally tolerated," "stably integrated") was replaced with measurable errors, conservation laws, and ranges of applicability.
- AMR was demoted to "the initial version is based on a uniform grid; AMR is an upper-level option (enabled once coarse-fine conservative interpolation, time-step ratio, and validation problems are defined)" (#29).

---

## 1. Runtime-mode configuration matrix (mutually exclusive) / design principles

**Design principle (fidelity default, relaxation is an add-on extension point)**: each mode axis is abstracted with a strategy/trait. **The default for every axis is the fidelity-first implementation (= reference solution)**, and low-cost approximations (MRF, point-bubble, one-way, AMR, aggressive f32) are added as **add-on extensions** behind the same trait. The structure is swappable without changing the core coupling loop (§5), and relaxation modes are validated against the tolerance to the corresponding fidelity reference solution (thresholds defined in §8 VR). The fidelity default has the maximum computational cost (the §7 memory budget is sized on the fidelity configuration).

All modes are implemented, but per computation exactly one is selected from each axis, exclusively. The config-validation layer rejects inconsistent combinations.

| Axis | Default (fidelity-first) | Relaxation extension (later) / reference tier | Exclusion constraints / notes |
|---|---|---|---|
| Rotation | `IBM-inertial` (unsteady, time-accurate) | Relaxation: `MRF-frozen-rotor` (steady approximation) / reference tier: `sliding-overset` | MRF cannot be combined with IBM moving blades (#6). Phase-averaged statistics only with IBM/overset (#37). |
| Interface | `resolved-phasefield` (conservative Allen-Cahn. Interface and mass-transfer fidelity prioritized) | Relaxation: `point-bubble` (Euler-Lagrange) / `hybrid` | Switching decision is `d_b/Δx, d_b/W, Eo, Re_b, α_g, We_b` (#12). hybrid defines interphase mass, momentum, and scalar conservation laws (§5). |
| Scalar | `active` (feedback to σ, viscosity, density, [temperature] enabled) | Relaxation: `passive` (feedback opt-out) | Feedback targets and stabilization of active are made explicit (#13). |
| Particle coupling | `two-way` (`four-way` at high `α_p`) | Relaxation: `one-way` | `α_p` / mass-loading threshold (#16). Accompanied by a reaction-force scattering kernel and momentum-conservation validation. |
| Precision | Fidelity profile: near-interface, conserved quantities, torque, interface curvature, and reductions are `f64`; only the far bulk is `f32` | Relaxation: aggressive `f32` / reference tier: all `f64` | #32, §7. |
| Grid | `uniform` (full resolution at the required resolution) | Relaxation: `block-AMR` (when coarse-fine conservative interpolation, time-step ratio, and validation are satisfied) | #29. AMR is an add-on due to implementation risk. |

---

## 2. Representative application problem / representative quantities / dimensionless numbers (stirred-tank reactor)

This section concretizes the **representative application** of the §8 validation benches; the functional requirements of §4
apply generally to rotating boundaries, high-density-ratio two-phase, LES, and scalar/particle coupling (consistent with the neutralization of the title).

3D cylindrical (or prismatic) vessel. Continuous phase (Newtonian/non-Newtonian liquid), dispersed gas phase from a bottom sparger (`ρ_l/ρ_g ≈ 10³`, `μ_l/μ_g ≈ 10²`), rigid-body rotating blades at constant angular velocity `Ω`, suspended particles near neutral buoyancy, multiple scalars with interfacial mass transfer and liquid-phase reactions. Target observables: time-/phase-averaged 3D velocity field, shear-stress field based on the second invariant of the strain rate, particle Lagrangian cumulative shear exposure, gas holdup `ε_g`, and dissolved-scalar concentration field.

### 2.1 Dimensionless-number definitions (representative quantities fixed, #26, #23, #24)

Representative rotational speed `N = Ω/(2π)` [rev/s], impeller diameter `D`, tank diameter `T`, liquid depth `H`, gravity `g`, bubble diameter `d_b`, particle diameter `d_p`, molecular diffusivity `D_m`, surface tension `σ`, `Δρ = ρ_l − ρ_g`.

```
Re   = ρ_l N D² / μ_l                 (stirring Reynolds, representative velocity U_tip = πND)
Fr   = N² D / g
We   = ρ_l N² D³ / σ
Eo   = Δρ g d_b² / σ                  (=Bond, bubble scale)
Mo   = g μ_l⁴ Δρ / (ρ_l² σ³)
Ca   = μ_l U / σ
Sc   = ν_l / D_m,   Pe = U L / D_m = Re·Sc
Da_n = k C_ref^{n-1} · (L/U)          (reaction order n, k is the rate constant; noted separately per order)
St   = τ_p / τ_f,   τ_p = ρ_p d_p² / (18 μ_l)
Np   = P / (ρ_l N³ D⁵),  P = 2π N T_q = Ω T_q   (T_q=torque; N in rev/s, ρ referenced to the liquid phase)
N_Q  = Q / (N D³)                     (Q=net volumetric flow rate at the blade discharge surface)
```

Lattice-side constraints: `Ma_lattice = U_tip/c_s ≤ 0.1`, Cahn number `Cn = W/L`, interface Péclet `Pe_φ = U W / M`, relaxation time `τ ∈ [τ_min, τ_max]`.

### 2.2 Matching priority order (when simultaneous matching is impossible, #25)

The physical→lattice conversion has finite degrees of freedom. When not all dimensionless numbers can be matched simultaneously, the priority order is fixed:
**(1) Re → (2) density ratio / viscosity ratio + We/Eo (interface dynamics) → (3) Fr (when the free surface / buoyancy dominates) → (4) Sc/Pe / Da (scalar / reaction) → (5) St (particles)**.
The unit-conversion layer must run a feasibility check, and when `Ma>0.1` / `τ∉[τ_min,τ_max]` / `Cn` is too large / diffusion-number or CFL violations occur, it warns while explicitly stating the compromised dimensionless numbers and the errors.

---

## 3. Governing equation system (revised)

```
Continuous phase (recovered by low-Mach LBM, well-balanced gravity with phase-wise density):
  ∂ρ/∂t + ∇·(ρu) = 0
  ∂(ρu)/∂t + ∇·(ρuu) = -∇p + ∇·[ (μ(γ̇)+μ_t)(∇u+∇uᵀ) ]
                        + F_s + ρ g + F_b^{scalar} + F_g^{disp} + F_p + F_rot
  - Gravity applies ρg to all phases and discretizes the hydrostatic pressure ∇p_hydro = ρg in a well-balanced manner (#34).
  - F_b^{scalar} (rev.2, active density feedback): the Boussinesq perturbation force of
    the solute buoyancy F_b = ρ_0 β_C (C−C_0) g. It is exactly 0 at C≡C_0 and is
    **not mixed** with the well-balanced hydrostatic cancellation of ρ(φ)g (composed as an independent force source. See docs/proposals/active-scalar-feedback.md for details).
  - F_rot is only for the MRF mode (§4.3). The composition of μ(γ̇) and μ_t is the implicit consistency of §4.7.

Two-phase interface (conservative Allen–Cahn phase field, fixed to the Fakhari 2017 family, #8):
  ∂φ/∂t + ∇·(φu) = ∇·[ M ( ∇φ − (4/W) φ(1−φ) n̂ ) ]
  n̂ = ∇φ / (|∇φ| + ε),  φ∈[0,1] (φ=1:liquid, φ=0:gas),  M[length²/time]
  Density/viscosity interpolation: ρ(φ)=ρ_g+φ(ρ_l−ρ_g),  1/μ or μ is fixed by a specified interpolation rule.

Surface tension (chemical-potential form as the baseline, with a convention branch when σ is variable. rev.2, #7):
  μ_φ = 4β φ(φ−1)(φ−1/2) − κ ∇²φ
  σ = √(2κβ)/6·(coefficient is model-defined), W = 4√(κ/(2β)); the relation fixes (σ,W)↔(κ,β)
  - When σ is constant (baseline form): F_s = μ_φ ∇φ (the CSF-equivalent σκn̂δ_s goes to the validation items)
  - When active with σ depending on C_k / temperature: F_s = μ_φ∇φ is not used directly; it is unified into a
    well-balanced CSF / chemical-potential combined form that avoids double-counting with the Marangoni tangential force
    (docs/proposals/active-scalar-feedback.md conventions D1/D2). The coefficients are frozen after being derived
    into the (κ,β,W,σ) convention of this section (derivation required — mandatory before implementation). A degeneracy
    test that agrees with the σ-constant baseline form under degeneracy to ∇σ=0 is placed in §8.

Dispersed gas phase (point-bubble mode, Euler-Lagrange, #12):
  m_b dv_b/dt = F_buoy + F_drag(Tomiyama) + F_lift + F_addedmass + F_walllub + F_TD
  Add the BIT (bubble-induced turbulence) generation term to the LES (§4.2, #46)

Dispersed particles (Euler-Lagrange, #16):
  m_p dv_p/dt = F_drag(Schiller-Naumann, Re_p range explicit) + F_buoy
               + [high-accuracy] F_Saffman + F_Basset + F_Faxen
  For two/four-way, scatter the reaction force with a regularization kernel and verify momentum conservation

Scalar / reaction (component k, active/passive explicit, #13, #14):
  ∂C_k/∂t + u·∇C_k = ∇·[ (D_k + ν_t/Sc_t) ∇C_k ] + R_k(C) + Ṡ_k^{if}
  - The SGS scalar flux is closed with the turbulent Schmidt number Sc_t (default Sc_t=0.7, variable).
  - Ṡ_k^{if}: the resolved interface uses a normal jump + partition coefficient; point-bubble uses k_L a(C*−C) (#35).

Eddy viscosity (SGS, Smagorinsky and WALE separated, #4):
  Smagorinsky: ν_t = (C_s Δ)² |S̄|,  |S̄|=√(2 S̄:S̄)
  WALE:        ν_t = (C_w Δ)² (S^d:S^d)^{3/2} / [ (S̄:S̄)^{5/2} + (S^d:S^d)^{5/4} ]
               S^d is the deviatoric symmetric part of the squared velocity-gradient tensor (local gradient recovery required)
```

---

## 4. Functional requirements (numerical methods, coefficients-fixed version)

### 4.1 Foundational LBM core
- **FR-CORE-01**: D3Q19/D3Q27 selectable. **The condition under which D3Q27 is the default is limited to "multiphase or strong forcing or cumulant use"** (#30). The Hermite order of the equilibrium distribution retained on each lattice and the recovery accuracy are defined separately (D3Q19 has restrictions on 3rd-order isotropy). **The M-F fidelity-default scenario falls under the multiphase / strong-forcing condition and is therefore always D3Q27**. Single-phase / weak-forcing derived scenarios (e.g., VR-STR-01 single-phase stirring) permit D3Q19 (rev.2, finding 11).
- **FR-CORE-02**: Implement central moments (cascaded) / cumulant. Stability is not "guaranteed" but specified by the **allowable relaxation-rate range on the target bench, positivity, and the presence/absence of regularization/filtering/entropic limiter** (#31).
- **FR-CORE-03**: Guo forcing. The velocity moment is `ρu = Σ c_i f_i + Δt F/2`. **In stress evaluation, subtract the forcing second-moment correction** (§4.6, #2).
- **FR-CORE-04**: `Ma_lattice ≤ 0.1`, compressibility error `O(Ma²)` controlled. Include the consistency of acoustic scaling and incompressibility in the unit-conversion feasibility (#25).

### 4.2 Turbulence model (LES-LBM)
- **FR-LES-01**: Implement Smagorinsky (including dynamic Germano) and WALE **as separate formulas**. WALE default (`ν_t→0` near the wall). **Because WALE requires the full velocity gradient, the "no finite differences" requirement is withdrawn**, and the local-gradient recovery method (moment or compact difference) is made explicit (#4).
- **FR-LES-02**: The relaxation-time reflection of eddy viscosity is **`τ_eff = 1/2 + (ν_0+ν_t)/(c_s²Δt)`** (general form). The simplified form `Δτ_t = 3ν_t` for lattice units `c_s²=1/3, Δt=1` is noted separately (#3 fix).
- **FR-LES-03**: The wall-shear-dominated region uses a `y⁺` wall function or wall-fitted interpolated boundary. `τ_eff` has not only the lower bound `>1/2` but also **upper-bound clipping and diagnostics** (avoiding over-diffusion / boundary-accuracy degradation, #27).
- **FR-LES-04**: Reflect the SGS scalar flux (turbulent Schmidt number `Sc_t`) and the SGS heat flux (turbulent Prandtl number) into the ADE-LBM relaxation time (#14).

### 4.3 Rotating impeller (mode-exclusive, #5, #6, #21, #22)
- **FR-ROT-01** (IBM-inertial): direct-forcing IBM (Uhlmann type). Target rigid-body velocity `U=Ω×r`. **"Guarantee Galilean invariance" is removed**, and **thresholds are set on slip velocity, torque error, and momentum-conservation error** for Taylor-Couette, rotating cylinder, and moving-wall Couette (including multi-direct-forcing / implicit IBM adoption conditions).
- **FR-ROT-02** (MRF-frozen-rotor): within the rotating zone, **solve the relative velocity `u_rel = u_abs − Ω×r`** and apply to the body force the Coriolis `−2ρ Ω×u_rel` and centrifugal `−ρ Ω×(Ω×r)`. **Do not apply MRF to stationary tank walls / baffles**. Define the velocity-matching condition at the rotating-zone boundary. Cannot be started simultaneously with IBM moving blades.
- **FR-ROT-03**: Clearly separate stationary walls/baffles = interpolated bounce-back (Bouzidi/Ginzburg), and **moving blades = IBM or moving-wall interpolated BB** (#22). Define the update frequency of the STL distance field and the geometric error during rotation.
- **FR-ROT-04**: Fix `Np = P/(ρ_l N³ D⁵)`, `P = Ω T_q`, `N = Ω/(2π)` (2π double-counting prohibited). During gas sparging, separately output **the ungassed `Np_0`, the gassed `Np_g`, and the gassed power-reduction ratio**. Define `N_Q = Q/(ND³)`, and the integration surface, velocity component, time/phase averaging, and backflow handling of `Q` (#23, #24).
- **FR-ROT-05** (sliding-overset, upper level): make overlapping-grid halo interpolation compatible with MPI. For reference-tier validation.

### 4.4 High-density-ratio two-phase flow
- **FR-VOF-01**: Conservative Allen-Cahn (fixed in §3). The mass-conservation error is specified per bench — set tolerances including time, grid resolution, and inflow/outflow amounts for a **closed static droplet / rising single bubble / sparger open boundary** (#9). Shan-Chen is not adopted for this use.
- **FR-VOF-02**: Spurious currents on a static droplet are `|u|_spurious·L/(σ/μ) = Ca_spurious < 10⁻³` (target We→0, resolution stated). Well-balanced chemical-potential form (implement the coefficient relation of #7).
- **FR-VOF-03** (sparger): Select from "gas-phase volumetric-flow-rate boundary / stochastic bubble injection / resolved orifice" (a simple `φ=1` + velocity Dirichlet alone is prohibited). State flow-rate conservation, pressure boundary, contact angle, and the lower bounds `d_b/W` and `d_b/Δx` (#10). Weaken breakup/coalescence to "numerically tolerated" and state that real thin-film drainage is not solved (#11).
- **FR-VOF-04** (point-bubble): Include in the switching conditions `d_b/W, Eo, Re_b, α_g, We_b, mass-transfer consistency`. Define interphase mass, momentum, and scalar conservation laws when hybrid is mixed (#12).
  **(rev.3, P1)** Population balance modelling (PBM) of the bubble-size distribution is
  required on the point-bubble path (breakup/coalescence kernels, e.g. Luo–Svendsen /
  Prince–Blanch): a mono-disperse point-bubble model cannot support the `d_32`
  acceptance of VR-STR-02 (internal consistency). Per-bubble gas-phase composition
  bookkeeping (component inventory and interfacial transfer budgets) must reconcile
  with FR-VOF-05. *Scope alignment (rev.2/§0)*: point-bubble is a relaxation extension
  (API-reserved in v1); this PBM requirement binds when that extension is implemented —
  in the resolved-phasefield default, `d_32` is measured from the resolved interface.
- **FR-VOF-05**: **Separate interfacial mass transfer into the resolved interface (normal flux, partition coefficient, phase-wise diffusion) and point-bubble (`k_L a(C*−C)`)** (#35). State the ranges of applicability of Henry's law and the Sherwood number.

### 4.5 Dispersed particles
- **FR-PART-01**: Switch one/two/four-way by `α_p`/mass-loading (thresholds explicit). Require the `Re_p` range of applicability of Schiller-Naumann, the reaction-force scattering kernel, and momentum-conservation validation (#16). For neutrally buoyant fine particles, decide the necessity of Saffman/Basset/Faxen by `d_p/Δx` and `St`.
- **FR-PART-02**: Can switch to a resolved-particle method (PSM/Noble-Torczynski, Ladd/Aidun-Lu).
- **FR-PART-03**: Record `∫γ̇dt` and `max γ̇` along the trajectory. **During LES tracking, enable SGS turbulent dispersion (stochastic dispersion)**, or state resolved-only (avoiding grid dependence of the exposure PDF/CDF, #17).

### 4.6 Stress-field evaluation (convention fixed, #1, #2, #18, #19, #20)
- **FR-STRESS-01**: The strain rate is evaluated locally from the non-equilibrium distribution. **The convention is fixed**:
  ```
  f_i^{neq} = f_i − f_i^{eq} (post-collision, pre-streaming; u includes the F/2 correction)
  Q_αβ = Σ_i c_iα c_iβ f_i^{neq} + (Δt/2)(u_α F_β + u_β F_α)   ← forcing second-moment correction
  S_αβ = − Q_αβ / (2 ρ c_s² τ_eff Δt)
  ```
  For cumulant/MRT, apply a coefficient correction with the shear-moment relaxation rate. Solve the circular dependence of the Smagorinsky closure in **algebraic closed form** (obtain `τ_eff` explicitly from `|Q|`; a Hou et al.-type quadratic).
- **FR-STRESS-02**: Define the output stress separately into **`resolved viscous` / `SGS` / `capillary` / `particle`**. Limit the source tensor of `γ̇=√(2S:S)`, the second invariant `II_S`, and von Mises (#19).
- **FR-STRESS-03**: Define wall shear per mode (**tangential velocity-gradient reconstruction / IBM forcing integration / MEM**). State the handling when the non-equilibrium quantity near the interpolated boundary does not represent the wall gradient. Validation includes a curved moving wall (#20).
- **FR-STRESS-04**: State the composition rule, iteration procedure, convergence criteria, `τ_min/τ_max`, and LES range of applicability for non-Newtonian `μ(γ̇)` (Carreau-Yasuda/Casson/power-law) and `μ_t` (avoiding double-counting / divergence, #18).

### 4.7 Boundaries, gravity, and initialization (new, #33, #34, #45, #47)
- **FR-BC-01** (top boundary, mandatory specification): Select from `closed` / `free-surface` / `degassing-outlet`. During sparging, require a gas-phase discharge outlet (a closed tank + gas-phase inflow only is non-physical due to gas accumulation, #33). Define the headspace pressure, liquid-surface deformation, and liquid-surface contact angle.
- **FR-BC-02** (gravity): Require `ρg` for all phases, dynamic/hydrostatic pressure decomposition, and a well-balanced hydrostatic test (`|u|<ε` in stationary stratification) (#34).
- **FR-BC-03** (wettability): Define the contact-angle boundary condition, slip/no-slip, and phase-field flux condition per wall (#47).
- **FR-BC-04** (scalar wall): Select no-flux / adsorption / reactive wall (#35).
- **FR-INIT-01**: Require the initial velocity/pressure/phase-field/scalar/particle placement, impeller ramp-up, gas-flow-rate ramp, statistics-sampling start time, and quasi-steady decision criterion (#45).

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

- **FR-COUP-01**: Default data flow "phase-field update → ρ/μ field update → force-source composition (`F_s+ρg+F_b^{scalar}+F_g+F_p+F_rot`) → fused collide-stream-moments → boundary → scalar ADE → reaction (split) → particle integration". **For strong coupling, stiff reactions, and surface-tension waves, require operator-splitting-error, subcycling, and iterative strong coupling**. Impose the constraints of the capillary time step `Δt_σ ≤ √(ρ̄ Δx³/(2πσ))`, particle `Δt_p`, and reaction ODE `Δt_r`.
- **FR-COUP-02**: The reaction solver switches explicit/implicit/Rosenbrock-BDF by stiffness detection. Define the **acceptance criteria for the negative-concentration limit, element-conservation error, and split error** (#15).
- **FR-COUP-03**: Dimensionless matching is the priority order of §2.2 + feasibility check (#25).
- **FR-COUP-04**: `probe_state_hash` bit-equivalence is **limited to single-backend implementation regression**. Physical validity / conservation laws are separate criteria (§8, #28, #42).
- **FR-COUP-05**: AMR is an upper-level option. When enabled, require coarse-fine conservative interpolation, time-step ratio, and dedicated validation (#29).

---

## 6. Input/output / visualization

- **FR-IO-01**: 3D field output is **uniform grid = VTI, structured curvilinear = VTS, unstructured/AMR = VTU/AMR** (#43). `ε_g` is defined per resolved-phase-field / point-bubble / hybrid (volume-average window and filter width stated; `φ` is a diffuse-interface indicator, not the void fraction, #36).
- **FR-IO-02**: Time-averaged / phase-averaged statistics (mean field, RMS, Reynolds stress). **Phase averaging only in IBM/overset unsteady modes**. MRF is output separately as a rotating-coordinate average / pseudo-steady (#37).
- **FR-IO-03**: 3D display in the Web GUI (slices, isosurfaces, shear heatmap, time-series probes). Extend the existing 2D canvas to WebGL/WebGPU.
- **FR-IO-04**: Histogram/CDF of particle cumulative shear exposure (presence/absence of SGS dispersion stated).
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

- **NFR-01 (scale / memory budget, #44)**: `O(10⁸–10⁹)` cells. **A budget table for bytes per cell, number of distributions (D3Q27 × phase field × scalar), particle count, GPU memory, I/O volume, and checkpoint frequency** is mandatory. A 1e9-cell × multi-distribution case is **0.6 TB class** at the fidelity default, and **TB to several-TB class** including all-f64, multiple scalars, simultaneous checkpoint retention, and I/O buffers (rev.2 fix), so attach the estimate for wgpu multi-GPU + MPI partitioning.

  **Budget table (rev.1b, fidelity-default configuration, deviation storage, ping-pong ×2. Per cell)**:

  | Component | Lattice/type | bytes/cell |
  |---|---|---|
  | Fluid distribution f | D3Q27 × 2 × f32 | 216 |
  | Phase-field distribution g (conservative Allen-Cahn) | D3Q19 × 2 × f32 | 152 |
  | Scalar distribution h (per component) | D3Q7 × 2 × f32 | 56 |
  | moments / property fields (ρ, u×3, φ, μ_φ, ∇φ×3, ν_t, γ̇, τ_eff) | 12 × f32 | 48 |
  | Mask / flags | u8×2 | 2 |
  | Statistics accumulators (mean u×3, RMS×3, Reynolds stress 6, etc.) | ~13 × f32 to f64 | 52–104 |
  | Interface-band f64 promotion (band width ~2W, 5–10% of all cells, amortizing +368 B/band-cell of f+g) | amortized | +18–37 |
  | Additional amount when including the curvature / reduction work area of the interface band (rev.2 recheck) | amortized | up to ~+40 |
  | **Total (1 scalar component)** | | **≈ 540–620 B/cell** |

  Conversion: **1e8 cells ≈ 56–62 GB** (single-node feasible on this machine's M5 Max 128 GB, upper limit ~1.5e8) /
  **1e9 cells ≈ 0.56–0.62 TB** (f32 bulk), and **≈ 1.1–1.2 TB** at the all-f64 reference tier.
  10⁷ particles × ~100 B = 1 GB (negligible). A checkpoint, saving the raw distributions,
  is the same amount as the field itself per instance (~0.5 TB/instance at 1e9) → the frequency is back-calculated
  from I/O bandwidth, with 2–5 per job as the default. GPU: 8–16 GB/card → 1.3–2.6e7 cells/card in the f32 configuration,
  and a 1e9-cell case **requires 40–80 cards of multi-GPU or a CPU-cluster MPI** (single GPU is infeasible).
  Conclusion: 1e9 at the fidelity default is cluster-only. Development/validation take ≤256³ (1.7e7 cells ≈ 10 GB) as the
  standard, and scale measurement is integrated into the R3 cluster plan (CLUSTER_OPTIONS.md).
- **NFR-02 (precision policy, #32, rev.2 vocabulary tidy-up)**: **The default is the fidelity profile** = `f32-bulk` + **near-interface, conserved quantities, torque, interface curvature, and reductions are f64** (same as §1). **The range of applicability of aggressive f32 (relaxation extension) is limited to "single-phase / weak coupling"**, and the relative degradation is validated with VR-STR-RELAX-f32. Make it consistent with the `ρ_l/ρ_g≈10³`, `Ca_spurious<10⁻³`, and mass-conservation requirements.
- **NFR-03 (performance)**: Integrate phase field, scalar, and forcing into the fused `step_band`; maintain the 3D extension of ring double-buffering and SoA plane-major.
- **NFR-04 (determinism)**: Reductions are in deterministic order. GPU/MPI is tolerance-based regression (bit-equivalence is single-backend only, #42).

---

## 8. Validation / acceptance criteria (quantified, wired as VALIDATION.md **T17**. With thresholds, #38–#42)

The validation tests are authored adversarially by codex/Opus from this specification and separated from the implementation. Unfixed tolerance bands
("±tolerance %" notation) are quantified by the existing protocol (implementation → characterization measurement → record the basis in
PHYSICS.md → freeze), and the frozen values are recorded in VALIDATION.md T17.

- **VR-STR-01 (single-phase stirring)**: Standard baffled tank (`D/T`, `C/T`, blade geometry, and baffle count fixed), specified Re range, ungassed. Rushton `Np` = experimental correlation ±tolerance %, and blade-discharge velocity profile matched against a `L2/L∞rel` threshold on PIV/LDA reference survey lines (#38).
  **(rev.3, P4) Reference datasets**: Wu & Patterson (1989) LDA; Deen et al. (2002)
  PIV (standard Rushton, D/T=1/3, 4 baffles); standard `Np` correlations. Numeric
  bands are frozen via the T17 experiment-driven protocol — not hardcoded here.
- **VR-STR-02 (gas-liquid, split into 02a/b/c in rev.2)**: **02a single bubble** = match `U_t` against the Grace diagram Eo-Mo-Re by relative error. **02b bubble swarm** = `ε_g` spatial distribution, swarm rise velocity (hindered rise), `d_32` when breakup/coalescence is allowed, and turbulence intensity (`ν_t` response) when BIT is used. **02c stirred-tank aeration** = experimental-correlation ratio of `ε_g, d_32, k_L a` (#39).
  **(rev.3, P4) References**: single bubble = Grace diagram (Eo-Mo-Re); aerated tank =
  published `ε_g`/`d_32`/`k_L a` data and correlations. In point-bubble / RELAX-PB
  evaluations, `d_32` presupposes the FR-VOF-04 population balance (P1); in the
  resolved-phasefield default it is measured by interface segmentation.
- **VR-STR-03 (shear / stress)**: Separate single-phase method of manufactured solutions (MMS), curved-surface Couette, rotating cylinder, non-Newtonian Poiseuille, and multiphase static droplet, and set the **grid convergence order** and `L2/L∞`. Survey-line design that accounts for the divergent severity of near-wall `L∞` (#40).
- **VR-STR-04 (scalar / reaction)**: Taylor-Aris dispersion, a reaction-diffusion front with known `Da`, and `k_L a` (state whether the computation formula is an interface integral or a correlation). Specify the tolerance, target `Pe/Da/Sc`, and boundary conditions for each (#41).
- **VR-STR-05 (coupled regression / conservation)**: `probe_state_hash` is limited to single-backend regression. **Set drift thresholds individually for the totals of mass, momentum, scalar, gas-phase volume, particle count, and energy-like quantities**. Energy-like quantities (kinetic energy, interfacial free energy, particle kinetic energy) are treated **as monitoring quantities for non-physical drift, not as exactly conserved** (rev.2). GPU/MPI is tolerance-based (#42).
- **VR-STR-06 (well-balanced)**: `|u|<ε` in stationary stratification (#34). **06+ (rev.2)**: with the active scalar ON and `C≡C_0`, satisfy the same stationarity (exact-zero degeneracy of `F_b^{scalar}`). The `∇σ=0` degeneracy of the variable-σ form (agreement with the σ-constant baseline form) is also placed in this group.
- **VR-STR-07 (initialization independence)**: Vary the spin-up / statistics-start conditions and the quasi-steady statistics agree within the threshold (#45).
- **VR-STR-RELAX (relaxation-mode equivalence, newly added in rev.2 — finding 1)**: relaxation extensions are accepted on the basis of relative degradation against the **corresponding fidelity baseline solution**. Comparison target and measured quantity for each axis (tolerances are characterize→freeze):
  - **RELAX-MRF**: against an IBM-inertial (or sliding-overset) baseline at the same geometry and same Re, the tolerances on `Np`, the impeller-discharge velocity profile line, the mean velocity field, and torque. Applicability is limited to configurations where the steady approximation holds.
  - **RELAX-PB (point-bubble)**: against a resolved-phasefield baseline, the tolerances on `ε_g`, `d_32`, `k_L a`, and the momentum/scalar balance. The applicable range is limited by `d_b/Δx, d_b/W, α_g` (identical to the switch-over condition of FR-VOF-04).
  - **RELAX-1W (one-way)**: against a two-way baseline, the tolerance on particle statistics and the upper mass-loading limit at which neglecting the reaction is permissible.
  - **RELAX-AMR**: against a uniform baseline, the conserved quantities, interface position, torque, velocity-field norm, and the balance error when crossing the coarse-fine boundary.
  - **RELAX-f32 (aggressive f32)**: against a fidelity-profile (or all-f64) baseline, the tolerances on conserved-quantity drift, `Ca_spurious`, `Np`, interface curvature, and reduced quantities. Applicability is limited to single-phase / weak coupling (NFR-02).

---

## 9. Major technical risks

| # | Hard part | Risk | Mitigation |
|---|---|---|---|
| 1 | High-density-ratio two-phase (`10³`) | Interface instability, parasitic currents, f32 rounding | well-balanced phase field + D3Q27 + f64 reduction, point-bubble alternative (§1) |
| 2 | High-Re stability | Divergence, hyperviscosity, positivity | cumulant + WALE, algebraic closure `τ_eff`, limiter (§4.1,4.6) |
| 3 | Rotating-boundary conservation | IBM slip, torque error | multi-direct-forcing, overset baseline validation, thresholding (§4.3) |
| 4 | Coupling stiffness | Time-scale divergence between reaction, interface, and rotation | operator-splitting error evaluation, subcycling, capillary dt (§5) |
| 5 | Compute cost / memory | 1e9 lattice × many distributions = TB scale | memory budget table, multi-GPU + MPI, AMR upper option (§7) |
| 6 | Mass transfer / BIT | Turbulence of unresolved aeration, insufficient transfer | Sc_t SGS, BIT source term, resolved/point separation (§4.2,4.4) |

---

## 10. Design decisions (settled items and remaining implementation details)

**Settled (rev.1a, PM decision)**: set the default of every axis to **fidelity first**, and implement low-cost approximations as bolt-on extension points (§1 design principles). The 4 previously unresolved items are resolved as follows:

- Interface default = `resolved-phasefield` (interface / mass-transfer fidelity first). `point-bubble` is a relaxation extension.
- Scalar default = `active` (with property feedback). `passive` is a relaxation extension.
- Lattice default = `uniform` (fully resolved). `block-AMR` is a relaxation extension.
- Precision default = fidelity profile (f64 near the interface, for conserved quantities, and for reductions; f32 in the far bulk). All-f64 is the baseline grade, aggressive f32 is a relaxation extension.

**Remaining implementation details (spec refinement, not decisions)** — status as of rev.3:

- Concrete formulas and stabilization (including Marangoni) for the feedback targets of the `active` scalar (σ / viscosity / density / [temperature]).
  → researched: docs/proposals/active-scalar-feedback.md. **One derivation is mandatory
  before implementation** (Marangoni coefficient consistency with the (κ,β) convention,
  §3). Thermal axis recommended as API-reserved extension.
- Allowed error thresholds for each relaxation extension against its fidelity baseline solution (added to the §8 VR).
  → structure defined as VR-STR-RELAX (rev.2); numeric bands frozen at relaxation
  implementation time.
- The f64/f32 boundary of the fidelity profile (band width near the interface, reduction range).
  → frozen experimentally during W-VOF implementation (characterize→freeze).
- API definition of the trait boundary (strategy swap point) for each mode axis.
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
