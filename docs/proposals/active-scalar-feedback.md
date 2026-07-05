# Proposal: Active Scalar Property Feedback and Stabilization Spec

**Document ID**: PROP-ACTIVE-SCALAR / **Status**: draft (spec drafted based on literature survey. Implementation not started)
**Parent**: [REQ_STIRRED_REACTOR.md](../REQ_STIRRED_REACTOR.md) §10 remaining implementation details, "concrete formulas and stabilization (including Marangoni) for the `active` scalar's feedback targets (σ, viscosity, density, [temperature])"
**Notation**: Follows REQ §3. φ∈[0,1] (φ=1: liquid, φ=0: gas), σ=√(2κβ)/6, W=4√(κ/(2β)), μ_φ=4β φ(φ−1)(φ−1/2)−κ∇²φ, Sc_t default 0.7. Component concentration C_k, default reference concentration C_0.
**Legend**: "**derivation required**" = places where the literature does not substantiate the formula, and the coefficients must be determined in-house from Chapman-Enskog / equilibrium profile integration before implementation.

---

## 0. Scope and Premises

REQ §1 fixes the scalar default to `active` (feedback to σ, viscosity, density, [temperature] enabled). This document drafts the **concrete formulas, double-counting avoidance, stabilization, validation, and API** for each feedback. The phase field follows the conservative Allen-Cahn form (Fakhari 2017 lineage, fixed in REQ §3), and surface tension is based on the chemical-potential form F_s = μ_φ∇φ. Feedback means "making property fields depend on C_k" — it does not change the evolution law of the phase field φ itself (so as not to break the interface-tracking equivalence-freeze values).

**Confirmation of REQ's coefficient relations** (making explicit that there is no contradiction): For the double-well f_b = β φ²(1−φ)², μ_φ = df_b/dφ − κ∇²φ = 2β φ(1−φ)(1−2φ) − κ∇²φ = 4β φ(φ−1)(φ−1/2) − κ∇²φ (consistent with REQ §3). The equilibrium profile φ_eq(x) = ½[1+tanh(2x/W)], W = 4√(κ/(2β)) = √(8κ/β), σ = ∫ κ(dφ/dx)² dx = √(2κβ)/6. Hereafter this (σ,W)↔(κ,β) correspondence is taken as the sole source of truth.

---

## 1. σ(C): Surfactant Dependence and Marangoni Force

### 1.1 Choice of Adsorption Isotherm

For the surfactant coverage Γ (dimensionless Γ=Γ*/Γ_∞*, Γ_∞*=saturation adsorption amount), two **equations of state** are provided, with Langmuir as the default:

- **Linear (default sub-tier, low-coverage regime Γ≲0.3)**: σ(Γ) = σ_0 (1 − E·Γ). Elasticity number E = 𝓡T Γ_∞*/σ_0 (the first-order term of a Langmuir expansion at Γ→0).
- **Langmuir (default, full coverage range)**: σ*(Γ*) = σ_0 [1 + (𝓡T Γ_∞*/σ_0) ln(1 − Γ*/Γ_∞*)], dimensionless form σ(Γ) = 1 + E ln(1−Γ) (Kwakkel/Prosperetti lineage, corresponds to the equation of state in van der Sman & van der Graaf 2006. Sources: arXiv:2409.19374, Springer Rheol. Acta 44:365 (2005)). Range of applicability: Γ<1 (since σ→−∞ as Γ→1, the lower-bound clip in §1.4 is mandatory).

### 1.2 Handling of Bulk C and Interfacial Adsorption Γ (Simplification for the M-F First Version)

Strictly speaking, Γ is an independent field constrained to the interface, coupled to bulk C via adsorption/desorption (Langmuir kinetics dΓ/dt = k_a C_s(Γ_∞−Γ) − k_d Γ) (Teigen et al. 2011, van der Sman 2006). **The M-F first version proposes the following two tiers**:

- **Tier-A (first-version default, simplified)**: Assume equilibrium adsorption, and algebraically determine Γ from the local interfacial bulk concentration C via the Langmuir equilibrium isotherm Γ_eq(C) = Γ_∞ K C/(1+K C) (K=k_a/k_d), using σ(Γ_eq(C)). **This allows direct dependence on bulk C** (an approximation for soluble, quasi-equilibrium surfactants). No additional scalar field is needed — the existing C_k ADE infrastructure (REQ §3) can be reused as-is. Explicit range of applicability: valid only when the interfacial adsorption timescale τ_ads=1/(k_a C+k_d) is sufficiently faster than the interface deformation timescale.
- **Tier-B (deferred extension point)**: Solve the interface-constrained field Γ̂=δ_Σ Γ via conservative transport (∂Γ̂/∂t + ∇·(Γ̂u) = (1/Pe_Γ)∇·[∇Γ̂ − (1−2φ)Γ̂/(√2 Cn)·∇φ/|∇φ|], profile-preserving, arXiv:2409.19374), to handle non-equilibrium adsorption (soluble/insoluble). Explicitly marked as unsupported in the first version.

### 1.3 Phase-Field-Consistent Form of the Marangoni Force (Double-Counting Avoidance)

**Baseline (well-balanced, default)**: The interfacial force for variable σ decomposes into normal (capillary) and tangential (Marangoni) components:

  F_s = [ −σ κ_φ n̂ + (I − n̂⊗n̂)·∇σ ] δ_Σ,  n̂ = ∇φ/(|∇φ|+ε),  δ_Σ = (3√2/4) W |∇φ|²

(Source: Liu, Wu, Ba, Xi et al., Phys. Rev. E 108, 055306 (2023) = arXiv:2306.11320. The coefficient of δ_Σ corresponds to the well-balanced phase-field discretization in that paper).

REQ's chemical-potential form F_s = μ_φ∇φ is a **normal-component-only** equivalent expression that presumes σ=constant. When σ is variable, simply adding a tangential term to this double-counts the normal contribution. **The double-counting avoidance convention is fixed as follows**:

- **Convention D1 (adopted)**: Evaluate the entire interfacial force as a field of σ(x) in the CSF form above. That is, **when σ is variable, do not use F_s = μ_φ∇φ**; instead, as a substitute in potential form, use the chemical-potential-combined form from Liu et al. (2306.11320):
  F_s = (3√2/4) W [ |∇φ|² ∇σ − (∇σ·∇φ) ∇φ + (σ/W²) M_φ ∇φ ]
  The 1st and 2nd terms are the tangential Marangoni term (the |∇φ|²-weighted expansion of (I−n̂n̂)∇σ), and the 3rd term is the normal capillary term that does not explicitly compute curvature (M_φ is the chemical potential, corresponding to μ_φ in this paper's notation). **Caution**: this paper's definition of μ_φ uses W² for the interface-width parameter, and whether the coefficients match REQ's (κ,β) convention is **derivation required** (the W²↔κ correspondence must be confirmed via the equilibrium profile — do not casually reuse it).
- **Convention D2 (consistency when degenerating to σ=constant)**: When ∇σ=0, the 1st and 2nd terms vanish, leaving only the 3rd term, which must match the existing F_s = μ_φ∇φ (normal capillary) form (this degeneracy test is included as a validation item; §6 VR-STR-06 lineage).
- **For σ=constant passive/non-surfactant cases, continue using F_s = μ_φ∇φ as before**, and switch to the variable-σ path (Convention D1) only when feedback.sigma is enabled (limiting the branch to a single location protects the equivalence freeze).

### 1.4 Stabilization

- **σ lower-bound clip**: Impose σ_min = c_σ·σ_0 (default c_σ=0.05, needs tuning) to prevent the Langmuir σ→−∞ divergence and negative surface tension (interface numerical collapse). Add the clip activation rate to the monitored quantities.
- **Lattice gradient evaluation of ∇σ**: Since σ is a composite function of φ and C, evaluate ∇σ = (∂σ/∂C)∇C using an isotropically weighted central difference (the lattice gradient ∇ψ ≈ (1/c_s²Δt)Σ_i w_i c_i ψ(x+c_i) for D2Q9/D3Q19). Gradient reconstruction from non-equilibrium moments (∇C ∝ Σ(g_i−g_i^eq)) is also an option, but for active scalars where SGS diffusion is present, the default is lattice-gradient evaluation of C (consistent with the gradient reconstruction in REQ FR-STRESS, **validation required**).
- **Timescale constraint**: Balance the CFL for Marangoni convective velocity U_Ma ~ |∇σ|·W/μ against the capillary dt (§5). In systems with large E (elasticity number), the Marangoni dt can become rate-limiting.

**Key literature**: Teigen, Song, Lowengrub, Voigt, *J. Comput. Phys.* 230:375 (2011) (soluble surfactant phase-field, independent Γ field); van der Sman & van der Graaf, *Rheol. Acta* 44:365 (2005) / *Comput. Phys. Commun.* (2006) (diffuse-interface free energy with Langmuir isotherm); Liu et al., Phys. Rev. E 108, 055306 (2023) (well-balanced variable-σ force for phase-field LBM); Kwakkel et al., arXiv:2409.19374 (2024) (profile-preserving surfactant + Marangoni, equation of state σ=1+E ln(1−Γ)); Farhat/Lee lineage, Zheng-Shu-Chew, *J. Comput. Phys.* 218:353 (2006) (classical implementation of phase-field LBM Marangoni).

---

## 2. μ(C): Viscosity Dependence

### 2.1 Typical Forms (3 types provided, default = exponential)

- **Linear (dilute, weak dependence)**: μ(C) = μ_0 [1 + k_μ (C−C_0)]. Range of applicability: k_μ(C−C_0)≪1.
- **Exponential (Arrhenius type, default)**: μ(C) = μ_0 exp[A(C−C_0)]. Of the form η=η_0 exp(A C_w) (Source: Arrhenius-type concentration dependence for liquids/colloidal dispersions. ScienceDirect "Arrhenius Equation" overview, ternary mixture systems ResearchGate 321081670). Guarantees positivity, stable over a wide range.
- **Krieger-Dougherty (suspensions/high loading, when the solute is akin to a solid volume fraction φ_p)**: μ(C) = μ_0 (1 − C/C_max)^(−[η]C_max), [η]=2.5 (spheres), C_max=maximum packing fraction (monodisperse ~0.64, polydisperse higher). Source: Krieger & Dougherty (1959), Anton-Paar wiki, Springer Rheol. Acta. Range of applicability: only suspensions where C can be regarded as a volume fraction. Since it diverges as C→C_max, the τ upper-bound clip in §5 is mandatory.

### 2.2 Composition Order into τ_eff (Consistency with REQ FR-LES-02)

REQ FR-LES-02's general form: τ_eff = 1/2 + (ν_0+ν_t)/(c_s²Δt). Viscosity feedback enters here as a **multiplicative correction on molecular kinematic viscosity**. Composition rule (**compose multiplicatively on the molecular side, then reflect into τ exactly once**, to avoid double-counting):

  ν_mol(γ̇, C) = [μ(γ̇) · f_C(C)] / ρ(C),   f_C(C) = μ(C)/μ_0

  ν_eff = ν_mol(γ̇, C) + ν_t,   τ_eff = 1/2 + ν_eff/(c_s²Δt)

- **Composition order (fixed proposal)**: (1) non-Newtonian shear dependence μ(γ̇) (determined from γ̇ via the iteration in REQ FR-STRESS-04) → (2) multiply by the concentration factor f_C(C) → (3) divide by ρ(C) to convert to kinematic viscosity → (4) add the SGS ν_t → (5) τ_eff. Since both μ(γ̇) and μ(C) are factors of molecular viscosity, **multiplicative composition** is physically appropriate (μ_t, the eddy viscosity, is additive). The non-Newtonian iterative convergence criterion (FR-STRESS-04) runs in an inner loop with f_C(C) held fixed, while C is frozen explicitly for one step (operator splitting).
- **Relationship to τ upper/lower-bound clipping**: REQ FR-LES-03's [τ_min, τ_max] is applied to the final τ_eff. When μ(C) is of KD type and swings toward divergence, the clip takes effect first, so monitor clip activation as a diagnostic for "property outside its valid range" (§5). The same applies on the τ_min side (low viscosity / high-C dilution).

---

## 3. ρ(C): Density Dependence (Solutal Buoyancy)

### 3.1 Boussinesq Approximation Form

Buoyancy due to solute concentration differences is added to the momentum equation as a **Boussinesq body force**:

  F_b = ρ_0 β_C (C − C_0) g,   β_C = −(1/ρ_0)(∂ρ/∂C)|_{T,p} (solutal expansion coefficient)

(As required by REQ. Distinct from the two-phase phase density ρ(φ), this treats solute density modulation within the liquid phase as a perturbation. The standard form for double-diffusive LBM. Source: LBM double-diffusive natural convection review, Springer JTAC 10973-022-11354-z, ResearchGate 6331163). In dimensionless form, characterized by the solutal Grashof/Rayleigh number Ra_C = β_C ΔC g L³/(ν D_m).

### 3.2 Consistency with Well-Balanced Gravity for High-Density-Ratio Two-Phase Flow (Condition for Not Breaking VR-STR-06)

REQ FR-BC-02's well-balanced hydrostatic discretization requires "static stratification with |u|<ε" (VR-STR-06). **Consistency conditions** when adding solutal buoyancy:

- **Condition C-B1**: F_b is exactly 0 at the reference state (C=C_0). Well-balancing is performed against the hydrostatic pressure ∇p_hydro=ρ(φ)g of the reference density ρ(φ), and F_b is added separately as a **perturbation force on top of that**. In a uniform field with C≡C_0, F_b≡0, so the static stratification test of VR-STR-06 must remain unchanged even with C feedback ON (this is added as a degeneracy test to the VR, §6).
- **Condition C-B2**: Do not fold F_b and ρ(φ)g together into a single well-balanced force (since the Boussinesq perturbation is a non-equilibrium driving force, treating it as a target of well-balanced cancellation would break the stratification). Implementation-wise, retain F_b as an independent term at the force-source composition stage (F_s+ρg+F_b+… , REQ FR-COUP-01).
- **Sign/direction**: g points downward. If β_C>0 (heavy solute), fluid parcels with C>C_0 sink. Verify the sign with solutal RT (§6).

### 3.3 Determining the Concentration Range Where Non-Boussinesq Treatment Is Needed

The Boussinesq approximation is valid for density modulation |Δρ/ρ_0| = |β_C ΔC| ≲ 0.1 (conventional value). For concentration ranges exceeding this:

- **Criterion**: Upon detecting |β_C(C_max−C_min)| > 0.1, warn via the feasibility check (REQ §2.2), and either (a) switch to a low-Mach variable-density path that directly reflects phase-density interpolation ρ(φ,C) into the momentum/continuity equations, or (b) narrow the target concentration range. The M-F first version defaults to the Boussinesq path, and non-Boussinesq is explicitly a deferred extension point (reason: rigorously reflecting variable density into the continuous-phase LBM requires additional pressure/mass-conservation validation, causing scope expansion).

---

## 4. [Temperature] (Optional Axis)

### 4.1 Temperature ADE and Property Feedback Formulas

When the temperature field T is added on the scalar ADE infrastructure (isomorphic to REQ §3's C_k):

  ∂T/∂t + u·∇T = ∇·[(α + ν_t/Pr_t) ∇T] + Ṡ_T,   Pr_t=turbulent Prandtl number (paired with REQ FR-LES-04)

Thermal property feedback:

- **ν(T)**: μ(T) = μ_0 exp[B(1/T − 1/T_0)] (Arrhenius, standard for liquid viscosity), or linear μ(T)=μ_0[1−b(T−T_0)]. Add f_T(T)=μ(T)/μ_0 as an additional multiplicative factor to §2's composition rule.
- **σ(T) (thermal Marangoni)**: σ(T) = σ_0 + σ_T(T−T_0), σ_T=∂σ/∂T (usually <0 for pure substances). Source: Liu et al., Phys. Rev. E 108, 055306 (2023) (linear σ-T and phase-field LBM thermocapillary). The Marangoni force is evaluated via the same mechanism as §1.3's CSF form, with ∇σ=(∂σ/∂T)∇T (additive with surfactant ∇σ=(∂σ/∂C)∇C, linear superposition).
- **ρ(T)**: Thermal buoyancy F_b^T = ρ_0 β_T(T−T_0)g, β_T=thermal expansion coefficient (additive with §3's solutal Boussinesq → double-diffusive).

### 4.2 Sharing the Scalar ADE Infrastructure and Implementation Recommendation

Temperature can share the **same ADE-LBM distribution and relaxation-time mechanism** as C_k (the only differences being that the feedback targets are σ_T/β_T/μ(T), and that Ṡ_T may include reaction heat/phase-change latent heat). Sharing potential is high.

**Recommendation (material for PM decision)**: The temperature axis should be treated as a **"deferred extension point" in the M-F first version, with only design hooks provided, excluding [temperature] from the first-version scope**. Reasons: (1) The formulas for thermal Marangoni/thermal buoyancy are mathematically isomorphic to the solutal versions, and can share the ADE infrastructure, CSF force, and Boussinesq force, so the cost of adding it later is low (little necessity to force it into the first version). (2) Latent heat, reaction heat, wall thermal boundary conditions (Neumann/Robin), and Pr_t validation would additionally be required, blurring the focus of the first version's stirred-reactor validation (VR-STR-01–07). (3) REQ §1 also brackets [temperature] as a conditional axis. → Reserving a temperature hook in the API (§7) while fixing the formulas in this document is the minimum-risk approach.

---

## 5. Cross-Cutting Stabilization Policy

### 5.1 Stability Conditions for Explicit Time Integration

Since all feedback forces are added explicitly, adopt the **minimum** of the following dt constraints (extending REQ FR-COUP-01's capillary dt):

- **Capillary dt**: Δt_σ ≤ √(ρ̄ Δx³/(2π σ)) (as stated in REQ. Resolves the shortest capillary wave. Source: generally Δt ≤ √(ρ_m Δ³/(π σ_max)), Brackbill lineage). For σ, **when σ is variable, use σ_max (not the minimum σ)** (the wave speed √(σk³/ρ) increases with σ).
- **Marangoni dt**: Δt_Ma ≤ C_Ma · μ̄ Δx / (|∇σ|_max W) (the advective CFL for Marangoni convective velocity U_Ma~|∇σ|W/μ. C_Ma≈0.5, **coefficient is derivation required/needs tuning** — no formulated phase-field Marangoni CFL was found in the literature; this is an estimate from advective scaling).
- **Buoyancy dt**: Δt_b ≤ C_b · √(Δx / |β_C ΔC g|) (the √(Δx/a) constraint for buoyancy acceleration a=|β_C ΔC|g. Thermal buoyancy β_T is analogous. C_b≈0.5, needs tuning).

### 5.2 Clipping/Relaxation of Feedback Quantities

- **σ clip**: §1.4's σ_min.
- **τ clip**: §2.2, REQ's [τ_min,τ_max].
- **Under-relaxation**: When updating property fields from C, in rapidly changing situations use χ^(n+1) = (1−ω)χ^(n) + ω χ_target(C) (χ∈{σ,μ,ρ_pert}, ω∈(0,1], default ω=1, ω<1 under stiff conditions). Used together with strong coupling/iteration (REQ FR-COUP-01).
- **Startup ramp**: Linearly ramp feedback strength 0→1 during startup (synchronized with REQ FR-INIT-01's impeller/gas ramp). Suppresses unphysical Marangoni spikes during the initial transient. The ramp time completes before the quasi-steady-state judgment.

### 5.3 Monitored Quantities to Add for Divergence Detection

In addition to REQ's existing monitoring: σ_min hit rate, τ clip activation rate, |β_C ΔC| (Boussinesq validity), effective Marangoni Reynolds/CFL number, maximum ∇σ, effective elasticity number E, (for Tier-B) Γ's excursion outside [0,1], and drift of the sum of Σσ over the interface. All trigger warnings above threshold, detected before divergence.

---

## 6. Validation Proposals (Additions to the VR-STR Series, 3–5 items each)

| ID (proposed) | Validation problem | Measured quantity | Tolerance (initial proposal) | Grid/parameters (initial proposal) |
|---|---|---|---|---|
| **VR-STR-08** | Thermal capillary droplet migration (Young-Goldstein-Block, small Ma, Re) | Terminal velocity V vs V_YGB = 2U/[(2+k̃)(2+3μ̃)], U=−σ_T G_T R/μ_B | L2rel(V) < 5% (Ma→0), 2nd-order gradient convergence | Droplet R=20–40 lattice units, μ̃=1, k̃=1, G_T constant, Ma=Re=O(0.1). Source: Young, Goldstein & Block (1959), J. Fluid Mech. 6:350 |
| **VR-STR-09** | Surfactant Marangoni droplet migration / steady surface-tension-gradient driving | Steady Marangoni convective velocity vs elasticity number E scaling | L2rel < 10% in the linear-E regime of velocity | Coverage-gradient droplet in quiescent fluid. E=0.1–0.5, Langmuir σ=1+E ln(1−Γ). Source: arXiv:2409.19374 |
| **VR-STR-10** | Viscosity-stratified Poiseuille flow (μ(C) two-layer) | Velocity profile vs analytical solution (two-layer viscosity Poiseuille) | L2rel(u) < 2%, interfacial velocity continuity | Parallel plates, C step upper/lower → viscosity ratio 2–10, steady state. μ(C) exponential form. Analytical: two-layer Couette/Poiseuille continuous-stress solution |
| **VR-STR-11** | Solutal Rayleigh-Taylor / solutal fingering instability | Growth rate, mixing-layer thickness vs linear stability/literature | Growth rate rel < 15% (early linear regime) | Specified Ra_C, constant β_C, 2D/3D. Boussinesq F_b. Sign check included as the non-stratified variant of VR-STR-06 |
| **VR-STR-06+ (degeneracy)** | Static stratification with active feedback ON | \|u\|_max | < ε (equal to VR-STR-06) | Uniform C≡C_0 giving F_b≡0, satisfying Convention D2 with σ=constant. Degeneracy test for §3.2 C-B1 |

Recommended adoption: as a minimum configuration, VR-STR-08 (thermal Marangoni, strongest since an analytical solution exists), VR-STR-10 (viscosity feedback, analytical solution exists), and VR-STR-06+ (degeneracy, regression for double-counting/stratification breakage) are mandatory. VR-STR-09 (surfactant) and VR-STR-11 (solutal buoyancy) are added when the corresponding feedback is implemented.

---

## 7. Draft API / Scenario Representation

Add a `feedback` object to `scalars[k]` in the scenario JSON (default is all null = passive equivalent; under REQ's active default, at least one must be enabled):

```jsonc
"scalars": [{
  "name": "O2",
  "feedback": {
    "sigma":     { "model": "langmuir", "sigma0": <σ_0>, "E": <elasticity number>,
                   "adsorption": { "K": <k_a/k_d>, "gamma_inf": <Γ_∞> },
                   "sigma_min_ratio": 0.05, "tier": "A" },   // A=bulk-C equilibrium, B=interfacial-Γ field (deferred)
    "viscosity": { "model": "exponential", "A": <A>, "C0": <C_0>,
                   "clip": "inherit_tau_bounds" },           // linear|exponential|krieger_dougherty
    "buoyancy":  { "model": "boussinesq", "beta_C": <β_C>, "C0": <C_0>,
                   "g": [0,0,-9.81] }                        // non_boussinesq is deferred
  }
}],
"temperature": {                                             // optional axis (first version is disable-hook only)
  "enabled": false,
  "feedback": { "sigma_T": <∂σ/∂T>, "beta_T": <β_T>, "viscosity": {"model":"arrhenius","B":<B>} },
  "pr_t": 0.85
}
```

**Range of applicability / mutual-exclusion validation items (additions to the config validation layer, REQ §1)**:

- `sigma.model=langmuir` with `sigma.tier=A` presumes "interfacial adsorption is quasi-equilibrium." Display a warning (do not reject) based on the ratio of τ_ads to the interface deformation timescale.
- `viscosity.model=krieger_dougherty` interprets C as a volume fraction, so when the same C_k is used together with `buoyancy=boussinesq` (which presumes a dilute solute), issue a consistency warning (dilute and high-loading are mutually exclusive physics).
- `sigma` feedback is only meaningful for phase-field scalars in modes where the interface (`resolved-phasefield`) is enabled (in `point-bubble`-only mode, the Marangoni force cannot be defined via an interfacial delta → reject, or delegate to the k_L a side).
- `buoyancy` presumes separation from well-balanced gravity (REQ FR-BC-02) per §3.2 C-B2. When enabled, make the VR-STR-06+ degeneracy test mandatory in CI.
- For multiple scalars' σ feedback (surfactant + temperature), ∇σ is linearly superposed (§4.1). The composition rule when multiple σ models are given to the same interface (additive or minimum) must be made explicit (first version: additive = linear superposition).
- Sc_t inherits the REQ default of 0.7 (overridable per `scalars[k]`). Temperature uses Pr_t (default 0.85 proposed, needs literature confirmation).

---

## 8. Remaining Open Points (Require PM Decision)

1. **Coefficient consistency for Convention D1 (most important)**: The chemical-potential-combined Marangoni form from Liu et al. (2306.11320) uses W² for the interface width, and whether the coefficients match REQ's (κ,β) (σ=√(2κβ)/6, W=4√(κ/(2β))) is **unconfirmed (derivation required)**. Before implementation, cross-check W²↔κ via the equilibrium profile, and numerically verify via the degeneracy test (Convention D2, VR-STR-06+).
2. **Whether to allow the Tier-A simplification for surfactants (direct bulk-C dependence) in the first version**: Since it is a quasi-equilibrium adsorption approximation, it introduces error in non-equilibrium, strongly deforming interfaces. Whether to include Tier-B (independent Γ field) in the first version or defer it is in tension with mass-transfer fidelity requirements (REQ §1's interface default = fidelity priority). This document proposes Tier-A for the first version, Tier-B deferred.
3. **Whether to include the temperature axis in scope**: §4.2 recommends "deferred extension point + API hook reservation." If included, additional VR (thermal capillary convection + thermal boundary) and memory budget (additional distributions) entries are needed.
4. **Coefficients for Marangoni/buoyancy dt** (C_Ma, C_b) are estimates, with no formula in the literature. Calibration is needed after implementation (back-calculated via VR-STR-08/11).
5. **Physical exclusivity of KD-type viscosity and Boussinesq buoyancy**: A single scalar simultaneously satisfying "high-loading suspension" and "dilute solutal buoyancy" is not typical. The question is how strongly the validation layer should enforce physical consistency of the feedback model per C_k.
