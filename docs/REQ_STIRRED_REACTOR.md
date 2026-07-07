# M-F Requirements: rotating boundary + high-density-ratio two-phase + LES + scalars/particles

**Document ID**: REQ-M-F-STR (active through rev.4; git history holds the
rev.1..rev.4 triage record).
**Target core**: `crates/lbm-core` (D3Q19/D3Q27, CpuScalar/CpuSimd/wgpu, MPI).
**Positioning**: M-F vertical of `docs/PLAN.md`. Acceptance = VALIDATION.md
**T17** (VR-STR-01..07 + RELAX).
**Representative application**: stirred-tank reactor. Requirements below are
domain-neutral; §2 concretizes the application.

## 0. Landed vs. pending (2026-07-07)

**Landed** — retained as one-line pointers only:
- **W-ROT** rotating-body direct-forcing IBM (bab1cae, merged 99bb32a) —
  IBM-inertial fidelity default. Test-side contract in ANOM-P4-009 (db50b24);
  8-probe accuracy audit f08aad0.
- **W-GRAV** well-balanced gravity (860cd7b, merged 79a539f).
- **W-LES** WALE core (1d3d692) + on-device omega pass (68feb27); frozen
  characterizations TGV64 `ν_eff` multimode (7e1e9f7) and Re_tau=178 vs MKM
  1999 DNS (6656089). WALE is the default; Smagorinsky implemented as a
  separate equation.
- **MF-alpha stage 3** central-moment collision on SIMD + GPU (20d0e10);
  cumulant offset triaged as resolution-point calibration (ANOM-P4-008
  disposition C, commits e569fb7 / 5eae598).
- **D3Q27 stage 1** validation (add5d5b, cx/d3q27-val).
- **D-track P2** one-way dispersed deposition — CR-1 interior sources, CR-2
  masked face patches, CR-3 deposition-aware stepping (3b1bcdc, 76b5071);
  T18.1-3 adversarial tests landed. Two-way/four-way particle coupling remains
  pending.

**Pending — the critical path**:
- **W-VOF** conservative Allen-Cahn free surface (fidelity-default interface).
  Owned by QA-sweep session per HANDOFF §4. Gates VR-STR-02 (single bubble /
  swarm / aeration) and blocks W-BCTOP (top surface / degassing / contact
  angle) and interfacial mass transfer.
- **W-BCTOP**, **W-BUB** (point-bubble + PBM), **active-feedback scalar**,
  **W-REACT** — all wait on W-VOF.

Sections retained below cover (a) W-VOF spec detail (unlanded), (b)
load-bearing conventions the codebase and its validation depend on
(dimensionless conventions §2, forcing-moment equation FR-STRESS-01, sparger
phase convention FR-VOF-03, ε_g output metadata FR-IO-01, precision profile
NFR-02), and (c) the DAG + parallel-wave plan.

---

## 1. Runtime-mode configuration matrix

Fidelity default per axis; relaxation approximations are bolt-on extensions
behind the same trait. Config validation rejects inconsistent combinations
(runtime mutual exclusion, not release-phase gating).

| Axis | Fidelity default | Relaxation / reference | Notes |
|---|---|---|---|
| Rotation | `IBM-inertial` (LANDED) | rel: `MRF-frozen-rotor` (PENDING) / ref: `sliding-overset` (PENDING) | MRF cannot combine with IBM moving blades. Phase-averaged stats IBM/overset only. |
| Interface | `resolved-phasefield` Allen-Cahn (**W-VOF PENDING**) | rel: `point-bubble` (PENDING) / `hybrid` (PENDING) | Switching by `d_b/Δx, d_b/W, Eo, Re_b, α_g, We_b`. |
| Scalar | `active` (feedback to σ, μ, ρ, [T]) (PENDING) | rel: `passive` (PENDING) | Feedback targets and stabilization explicit. |
| Particle | `two-way` (PENDING, FR-PART target); `four-way` at high `α_p` (PENDING) | rel: `one-way` (LANDED D-track) | Four-way = Phase 2 (FR-PART-04..06 contract). |
| Precision | `mixed_safe` (NFR-02 profile PENDING; explicit `f32`/`f64` LANDED) | ref: `full_f64` (profile PENDING) / rel: `mixed_fast` (PENDING) | ρ-ratio limit gates `mixed_fast`. |
| Grid | `uniform` (LANDED) | rel: `block-AMR` (PENDING) | AMR needs conservative coarse-fine + dt-ratio + validation. |

**Initial delivery (release gate target, not current landed state)**:
IBM-inertial + resolved-phasefield + active scalar (predictor-corrector) +
two-way particles + uniform + `mixed_safe`. **Current landed relaxation
subset**: IBM-inertial + one-way particles + uniform grids + explicit
`f32`/`f64`; current multiphase is Shan-Chen, not W-VOF; scalar ADE/feedback
and precision profiles remain pending. **Phase 2 (API-reserved)**: MRF,
point-bubble + PBM, four-way contact, block-AMR, `mixed_fast`, hybrid
interface, thermal axis. Accepted via VR-STR-RELAX.

---

## 2. Dimensionless conventions (load-bearing)

Representative quantities: `N = Ω/(2π)` [rev/s], `D` (impeller), `T` (tank),
`H` (depth), `g`, `d_b`, `d_p`, `D_m`, `σ`, `Δρ = ρ_l - ρ_g`.

```
Re     = ρ_l N D² / μ_l                (U_tip = πND)
Fr     = N² D / g
We     = ρ_l N² D³ / σ
Eo     = Δρ g d_b² / σ
Mo     = g μ_l⁴ Δρ / (ρ_l² σ³)
Ca     = μ_l U / σ
Sc     = ν_l / D_m
Pe_N   = N D² / D_m = Re·Sc            (impeller velocity scale ND)
Pe_tip = U_tip D / D_m = π·Re·Sc       (U_tip = πND) — every use site MUST
                                       state which Pe; no bare "Pe"
Da_n   = k C_ref^{n-1} · (L/U)         (reaction order n)
St     = τ_p / τ_f, τ_p = ρ_p d_p² / (18 μ_l)
Np     = P / (ρ_l N³ D⁵), P = Ω T_q    (T_q = torque; N in rev/s; liquid basis)
N_Q    = Q / (N D³)                    (Q = net vol flow through discharge)
```

Lattice: `Ma_lattice = U_tip/c_s ≤ 0.1`, `Cn = W/L`, `Pe_φ = U W / M`,
`τ ∈ [τ_min, τ_max]`.

**Matching priority when simultaneous matching impossible**:
(1) Re → (2) ρ ratio, μ ratio + We/Eo → (3) Fr → (4) Sc/Pe/Da → (5) St.

Unit-conversion layer runs a feasibility check; on `Ma>0.1`, `τ∉[τ_min,τ_max]`,
excessive `Cn`, or diffusion-number / CFL violation, warn with compromised
dimensionless numbers explicit.

---

## 3. Governing equations (retained conventions)

**Continuous phase** (low-Mach LBM, well-balanced gravity). Density-flux
consistency with phase-field diffusion (REV-CFD-CR-002) — `J_ρ = (ρ_l−ρ_g) J_φ`
**must appear identically in both** continuity and momentum advection
(consistent/AGG-type; mandatory at ρ ratio 10³):

```
∂ρ/∂t + ∇·(ρu + J_ρ) = 0,        J_ρ = (ρ_l − ρ_g) J_φ
∂(ρu)/∂t + ∇·[(ρu + J_ρ) u]
    = -∇p + ∇·[(μ(γ̇) + μ_t)(∇u + ∇uᵀ)]
      + F_s + ρg + F_b^{scalar} + F_g^{disp} + F_p + F_rot
```

Same discrete `J_ρ` in both equations (single code path; VR-STR-03/05
verifies). `F_b^{scalar} = ρ_0 β_C (C−C_0) g`, exactly 0 at `C≡C_0`,
independent of the `ρ(φ)g` well-balanced cancellation. `F_rot` MRF-only.

**Two-phase interface — W-VOF PENDING**. Conservative Allen-Cahn (Fakhari
2017), explicit conservative-flux form:

```
∂φ/∂t + ∇·(φu + J_φ) = 0,   J_φ = −M [∇φ − (4/W) φ(1−φ) n̂]
n̂ = ∇φ / (|∇φ| + ε), φ∈[0,1] (φ=1 liquid, φ=0 gas)
ρ(φ) = ρ_g + φ(ρ_l − ρ_g)
1/μ(φ) = φ/μ_l + (1−φ)/μ_g          (harmonic-in-μ, default frozen
                                     REV-CFD-MJ-013; linear-in-μ/ν are
                                     opt-in, logged, not gate-covered)
```

**Surface tension** (chemical-potential form):

```
μ_φ = 4β φ(φ−1)(φ−1/2) − κ ∇²φ
σ = √(2κβ)/6,   W = 4√(κ/(2β))    ← these ARE the definitions
```

σ constant (reference): `F_s = μ_φ ∇φ` (CSF-equivalent `σκn̂δ_s` for
validation). σ active: unified well-balanced CSF / chemical-potential form,
avoids double-counting with Marangoni; `∇σ=0` degeneration must agree with
σ-constant reference (VR-STR-06+).

**Point-bubble (Phase 2)**:
`m_b dv_b/dt = F_buoy + F_drag(Tomiyama) + F_lift + F_addedmass + F_walllub + F_TD`.
BIT production term added to LES.

**Particles (one-way D-track LANDED; two/four-way PENDING)**:
`m_p dv_p/dt = F_drag(Schiller-Naumann) + F_buoy [+ F_Saffman/F_Basset/F_Faxen]`.
Two/four-way: regularized reaction-force scatter; momentum conservation
validated.

**Scalars/reactions** (REV-CFD-MJ-011 — conservative forms normative for
two-phase and active-ρ; non-conservative single-phase is a special case):

```
Single-phase passive (ρ, α uniform):
   ∂C_k/∂t + u·∇C_k = ∇·[(D_k + ν_t/Sc_t)∇C_k] + R_k(C) + Ṡ_k^{if}

Two-phase phase-wise conservative (q∈{gas,liquid}, α_liq=φ, α_gas=1−φ):
   ∂(α_q C_{k,q})/∂t + ∇·(α_q u C_{k,q})
       = ∇·[α_q(D_{k,q} + ν_t/Sc_t)∇C_{k,q}] + α_q R_{k,q} + S_{k,q}^{if}

Density-based active:
   ∂(ρY_k)/∂t + ∇·(ρ u Y_k + J_k) = R_k + S_k^{if}
```

SGS flux closed with `Sc_t` (default 0.7). `S^{if}`: resolved interface =
normal jump + Henry partition (`S_{k,liq}^{if} = −S_{k,gas}^{if}`, flux
positive into liquid); point-bubble = `k_L a(C* − C)`. Conservation:
`Σ_q ∫ α_q C_{k,q} dV` changes only by boundary fluxes and reactions
(VR-STR-05).

**Eddy viscosity (LANDED W-LES)** — WALE default, Smagorinsky (incl.
dynamic Germano) as separate equation:

```
Smagorinsky: ν_t = (C_s Δ)² |S̄|,        |S̄| = √(2 S̄:S̄)
WALE:        ν_t = (C_w Δ)² (S^d:S^d)^{3/2}
                    / [(S̄:S̄)^{5/2} + (S^d:S^d)^{5/4}]
```

`S^d` = deviatoric symmetric part of the square of the velocity-gradient
tensor; local gradient reconstruction required.

---

## 4. Functional requirements

### 4.1 Core LBM
- **FR-CORE-01**: D3Q19/D3Q27 selectable. D3Q27 default when multiphase or
  strong forcing or cumulant; **M-F fidelity default is always D3Q27**.
- **FR-CORE-02**: Cumulant / central-moment collision **LANDED** (MF-alpha
  stage 3).
- **FR-CORE-03**: Guo forcing. `ρu = Σ c_i f_i + Δt F/2`. Stress evaluation
  uses the single equation in FR-STRESS-01 — prose ("add"/"subtract") is
  banned; the equation is the only definition.
- **FR-CORE-04**: `Ma_lattice ≤ 0.1`, `O(Ma²)` compressibility control.

### 4.2 W-LES turbulence — LANDED
- **FR-LES-01**: WALE (default) and Smagorinsky (incl. dynamic Germano) as
  separate equations. WALE needs the full velocity gradient — local gradient
  reconstruction (moments or compact differences) explicit.
- **FR-LES-02**: `τ_eff = 1/2 + (ν_0 + ν_t)/(c_s² Δt)`; lattice-unit
  simplification `Δτ_t = 3 ν_t`.
- **FR-LES-03**: Wall region uses `y⁺` wall function or wall-fitted
  interpolated boundary. `τ_eff` has both lower bound `>1/2` and **upper
  clipping + diagnostics** (avoid over-diffusion + boundary degradation).
- **FR-LES-04**: SGS scalar flux (`Sc_t`) and SGS heat flux (`Pr_t`)
  reflected in ADE-LBM relaxation time.

### 4.3 W-ROT rotating impeller — IBM-inertial LANDED; MRF Phase 2
- **FR-ROT-01** (IBM-inertial): direct-forcing IBM (Uhlmann type). Target
  `U = Ω × r`. Thresholds on slip velocity, torque error, momentum-conservation
  error set for Taylor-Couette, rotating cylinder, moving-wall Couette.
- **FR-ROT-02** (MRF): inside rotating zone, solve `u_rel = u_abs − Ω × r`;
  impose Coriolis `−2ρ Ω × u_rel` and centrifugal `−ρ Ω × (Ω × r)`. **MRF
  not applied to stationary walls / baffles.** Cannot start with IBM
  moving blades. Phase 2.
- **FR-ROT-03**: Stationary walls / baffles = interpolated bounce-back
  (Bouzidi/Ginzburg); moving blades = IBM or moving-wall interpolated BB.
- **FR-ROT-04**: `Np = P/(ρ_l N³ D⁵)`, `P = Ω T_q`, `N = Ω/(2π)` (2π double-
  counting prohibited). During aeration: `Np_0`, `Np_g`, gassed power-drop
  ratio output separately. `N_Q = Q/(ND³)`; integration surface, velocity
  components, time/phase averaging, backflow handling of `Q` defined.
- **FR-ROT-05**: sliding-overset (advanced reference-tier).

### 4.4 W-VOF two-phase — PENDING (critical path)
- **FR-VOF-01**: Conservative Allen-Cahn (§3). Mass-conservation error
  stipulated per bench (closed static droplet / rising single bubble /
  sparger open boundary). **Shan-Chen not adopted here.**
- **FR-VOF-02** (REV-CFD-MJ-005, dimensional fix): `Ca_spurious = μ_l
  |u|_spurious / σ < 10⁻³` (target We→0, resolution stated). The old
  `|u|·L/(σ/μ)` form was dimensionally void. Length-bearing indicator
  (if wanted): `Re_spurious = |u|_spurious L / ν_l` — never called Ca.
  Well-balanced chemical-potential form.
- **FR-VOF-03 (sparger; REV-CFD-CR-001 — phase-inversion fix — LOAD-BEARING)**:
  the sparger injects **GAS**. Under §3 (φ=1 liquid, φ=0 gas), the injected
  phase value is **φ=0**. Requirements:
  - Choose from gas-phase volumetric-flow / stochastic bubble injection /
    resolved orifice. **Plain `φ=0` + velocity Dirichlet alone is banned**
    — injection must simultaneously satisfy gas volumetric-flow
    conservation, pressure consistency, contact angle, and `d_b/W, d_b/Δx`
    lower bounds.
  - **Schema never exposes raw φ for inlets**: config uses
    `inlet_phase: gas | liquid`; core maps it (gas→φ=0, liquid→φ=1);
    enforced by config validation. Outputs report `φ_liquid` and
    `α_g = 1 − φ`.
  - Acceptance: gas-inlet unit test injects φ=0; sparger-only case balances
    injected gas volume vs domain gas-volume increase (VR-STR-02c precursor);
    no schema field accepts a raw φ boundary value.
- **FR-VOF-04** (point-bubble, Phase 2): switching by `d_b/W, Eo, Re_b, α_g,
  We_b, mass-transfer consistency`. Interphase mass/momentum/scalar
  conservation in hybrid case defined. **Population balance (PBM) mandatory**
  on this path (Luo-Svendsen / Prince-Blanch kernels) — mono-disperse
  point-bubble cannot support `d_32` acceptance. In resolved-phasefield
  default, `d_32` is measured by interface segmentation.
- **FR-VOF-05**: Interfacial mass transfer separated: resolved interface
  (normal flux + Henry partition + phase-wise diffusion) vs point-bubble
  (`k_L a(C* − C)`). Henry and Sherwood applicability explicit.

### 4.5 Dispersed particles — one-way D-track P2 LANDED; two/four-way PENDING
- **FR-PART-01..06**: one/two/four-way switching by `α_p`/mass-loading
  (thresholds explicit); Schiller-Naumann `Re_p` range explicit; reaction-
  force scatter kernel + momentum-conservation validation. LES tracking
  needs SGS turbulent dispersion (stochastic) OR resolved-only statement.
  Four-way (Phase 2) requires a soft-sphere normal-collision model with
  explicit `e_n, T_col, k, η, δ_max, Δt_p (≲ T_col/10)`. Lubrication
  correction when `d_p/Δx` doesn't resolve the gap. **Config guard
  (initial delivery)**: runs exceeding two-way regime `α_p`/mass-loading
  threshold are rejected at config validation (A-4 style) with threshold and
  source in error message.

### 4.6 Stress fields — FR-STRESS-01 LOAD-BEARING
- **FR-STRESS-01 (REV-CFD-CR-003, CR-004 — stage convention + forcing
  correction fixed by equations, not prose)**. Default stage is
  **pre-collision / post-streaming** (the coefficient below is derived for
  this stage):

  ```
  f_i^{neq,pre} = f_i^{pre} − f_i^{eq}(ρ, u)        (u includes F/2)
  Π_neq_raw     = Σ_i c_iα c_iβ f_i^{neq,pre}
  Π_force       = −(Δt/2)(u_α F_β + u_β F_α)        (Guo second moment,
                                                     for THIS engine's u/f_eq)
  Π_neq_corr    = Π_neq_raw − Π_force  =  Π_neq_raw + (Δt/2)(uF + Fu)
  S_αβ          = − Π_neq_corr / (2 ρ c_s² τ_eff Δt)
  ```

  `Π_neq_corr` is the **ONLY** normative definition. The exact sign of
  `Π_force` is derivation-frozen against this engine's Guo discretisation
  **before implementation** and locked by a negative test (body-force
  Poiseuille FAILS with the sign flipped — VR-STR-03).

  **If post-collision / pre-streaming stage is used instead**, the stage
  transform is mandatory:
  BGK: `Π_neq,pre = Π_neq,post / (1 − 1/τ_eff)`;
  MRT/cumulant: apply inverse shear-moment relaxation `R(τ_shear)⁻¹`.
  Then proceed with the equations above.

  The stress-evaluation API takes a required `neq_stage` enum
  (`PreCollision | PostCollision`) — no default-by-silence; misuse is a
  compile- or construct-time error.

  For cumulant/MRT, coefficient corrected by shear-moment relaxation rate.
  Smagorinsky closure's circular dependence solved in **algebraic closed
  form** (`τ_eff` from `|Q|`; Hou et al. quadratic).
- **FR-STRESS-02**: Output stress separated into **resolved viscous / SGS /
  capillary / particle**. For `γ̇ = √(2S:S)`, `II_S`, von Mises — source
  tensor restricted.
- **FR-STRESS-03**: Wall shear defined per mode (tangential-velocity-
  gradient reconstruction / IBM forcing integral / MEM). Handling when
  non-eq quantity near interpolated boundary doesn't represent the wall
  gradient explicit. Validation includes curved moving wall.
- **FR-STRESS-04**: Composition, iteration, convergence, `τ_min/τ_max`, and
  LES applicability of `μ(γ̇)` (Carreau-Yasuda / Casson / power-law) and
  `μ_t` explicit.

### 4.7 Boundaries, gravity, initialization
- **FR-BC-01** (top boundary, **W-BCTOP waits on W-VOF**): choose from
  `closed` / `free-surface` / `degassing-outlet`. During sparging, a gas-
  phase exhaust outlet is required. Headspace pressure, free-surface
  deformation, contact angle defined.
- **FR-BC-02 (gravity)** — **LANDED W-GRAV**: `ρg` on all phases; dynamic-
  pressure/hydrostatic decomposition; well-balanced hydrostatic test
  (`|u| < ε` in static stratification).
- **FR-BC-03** (wettability): per-wall contact angle, slip/no-slip,
  phase-field flux condition.
- **FR-BC-04** (scalar wall): no-flux / adsorption / reactive wall.
- **FR-INIT-01**: Initial fields, impeller ramp-up, gas ramp, statistics-
  sampling start time, quasi-steady decision criterion required.

### 4.8 Extension & closure contracts
- **FR-EXT-01**: Explicit contracts for §1 trait/strategy extension points
  and user-supplied closures (`R_k`, `μ(γ̇)`, body-force sources,
  relaxation-mode implementations). Requirements: I/O signatures with
  physical vs lattice units explicit; determinism (bit-identical outputs);
  GPU evaluability (state-free, wgpu-portable); NaN/divergence detection at
  the contract boundary; schema versioning + backward compatibility. Primary
  boundary = Rust traits; foreign-language ABI/SDK is a separate API spec.
  Co-designed with SOLVER_IMPROVEMENT_SPEC B-1. **First design note landed
  2026-07-06 (commit 1758814, kernel extension points).**

---

## 5. Coupling / time integration

- **FR-COUP-01 (REV-CFD-MJ-007 — dataflow split by scalar mode)**:
  - **Passive**: phase-field → ρ/μ → force composition → fused collide-
    stream-moments → boundary → scalar ADE → reaction (split) → particle
    integration.
  - **Active (fidelity default = predictor-corrector)**: scalar/reaction
    predictor → property update `ρ(C), μ(C), σ(C)[, T]` → force composition
    (incl. `F_b^{scalar}`, Marangoni) → flow step → scalar ADE corrector →
    reaction corrector → property re-evaluation (→ optional flow-scalar
    iteration for stiff coupling).
  Time-lagged explicit feedback allowed only as flagged relaxation
  `active_scalar_lagged=true`, with stability conditions + lag-error
  benchmark. Mode logged in metadata. For strong coupling / stiff reactions
  / surface-tension waves: operator-splitting error, subcycling, iterative
  strong coupling required. Constraints: capillary `Δt_σ ≤ √(ρ̄ Δx³ /
  (2πσ))`, particle `Δt_p`, reaction `Δt_r`. Acceptance: on active-scalar
  bench (Marangoni or concentration-dependent viscosity), feedback error
  converges under dt-halving.
- **FR-COUP-02**: Reaction solver switches explicit / implicit / Rosenbrock-
  BDF by stiffness. Negative-concentration limiting, element conservation,
  split-error acceptance defined.
- **FR-COUP-03**: Dimensionless matching per §2 priority + feasibility.
- **FR-COUP-04**: `probe_state_hash` bit-equivalence single-backend only.
- **FR-COUP-05**: AMR advanced; coarse-fine conservative interpolation,
  dt ratio, and dedicated validation required when enabled.

---

## 6. I/O / visualization

- **FR-IO-01**: 3D field output — uniform=VTI, structured-curved=VTS,
  unstructured/AMR=VTU/AMR. `φ` is the diffuse-interface indicator, not
  the void fraction.
  **ε_g processing definitions (LOAD-BEARING — every ε_g output carries
  filter width, averaging volume, time window as metadata)**:
  - resolved-phasefield: `ε_g_raw = ⟨1 − φ⟩_V` **and**
    `ε_g_thresholded(φ_c) = volume(φ < φ_c) / V`, default `φ_c = 0.5` —
    both output.
  - point-bubble: `ε_g_bubble = Σ_b V_b W_kernel(x − x_b) / V_filter`.
  - hybrid: `ε_g_total = ε_g_resolved + ε_g_bubble` with double-count
    exclusion over the resolved region.
  Any ε_g must be recomputable from a snapshot; comparisons state which
  definition was used.
- **FR-IO-02**: Time / phase-averaged statistics. Phase averaging IBM/overset
  only; MRF output as rotating-frame average / quasi-steady.
- **FR-IO-03**: 3D display in web GUI (slices, isosurfaces, shear heatmap,
  time-series probes). Extend existing 2D canvas to WebGL/WebGPU.
- **FR-IO-04**: Histogram/CDF of particle cumulative shear exposure (SGS
  dispersion presence/absence explicit).
- **FR-IO-05** (mixing metrics): blend time (until tracer CoV < threshold)
  and RTD (`E(t)`, mean residence time, variance). Thresholds and injection/
  detection surfaces explicit per scenario.
- **FR-IO-06** (large-scale I/O + resilience): Full dumps impractical at
  target scale (§7); require parallel I/O (HDF5/ADIOS2) + compression +
  in-situ statistics. **Deterministic checkpoint/restart with crash recovery**
  (bit-reproducible resume incl. RNG, particle state, statistics
  accumulators) mandatory. Builds on SOLVER_IMPROVEMENT_SPEC B-5 (snapshot),
  C-3 (per-rank parallel I/O), C-8 (distributed checkpoint) — reuse, don't
  duplicate.

---

## 7. Non-functional

- **NFR-01 (scale / memory budget)**: `O(10⁸–10⁹)` cells.

  **Budget table (fidelity default, deviation storage, ping-pong ×2)**:

  | Component | Lattice/type | B/cell |
  |---|---|---|
  | Fluid distribution f | D3Q27 × 2 × f32 | 216 |
  | Phase-field g (Allen-Cahn) | D3Q19 × 2 × f32 | 152 |
  | Scalar h (per component) | D3Q7 × 2 × f32 | 56 |
  | Moments/properties (ρ, u×3, φ, μ_φ, ∇φ×3, ν_t, γ̇, τ_eff) | 12 × f32 | 48 |
  | Mask / flags | u8 × 2 | 2 |
  | Statistics accumulators | ~13 × f32-f64 | 52-104 |
  | Interface-band f64 promotion (~2W band, 5-10% cells) | amortized | +18-37 |
  | Curvature / reduction workspace overhead in band | amortized | up to +40 |
  | **Total (1 scalar component)** | | **≈ 540-620 B/cell** |

  Conversion: **1e8 cells ≈ 56-62 GB** (single-node feasible on M5 Max
  128 GB, upper limit ~1.5e8) / **1e9 cells ≈ 0.56-0.62 TB** (f32 bulk) /
  ~1.1-1.2 TB (all-f64 reference). Particles 10⁷ × ~100 B = 1 GB (negligible).
  Checkpoint dump ≈ field-data size; frequency back-calculated from I/O
  bandwidth, default 2-5 dumps/job. GPU 8-16 GB/card → 1.3-2.6e7 cells/card
  in f32 — **1e9 cells require 40-80 cards multi-GPU or CPU-cluster MPI**.
  **Conclusion: 1e9 fidelity-default is cluster-only**. Development /
  validation ≤256³ (1.7e7 cells ≈ 10 GB); scale measurements integrated with
  the R3 cluster plan (CLUSTER_OPTIONS.md).
- **NFR-02 (precision profile — LOAD-BEARING, REV-CFD-MJ-009 — bound by
  array design / GPU kernels / memory budget)**:
  `precision_profile ∈ { full_f64, mixed_safe (default), mixed_fast }`.
  - **full_f64** (reference): all distributions, phase field, scalars,
    particle statistics, reductions in f64.
  - **mixed_safe** (fidelity default): bulk distributions f32; **f64 fixed
    for** `φ, ∇φ, κ, μ_φ, F_s, ρ(φ), μ(φ)`, distributions inside the
    interface band, all global reductions, torque, `Np`, `N_Q`, mass/volume
    counters, particle cumulative exposure.
    `interface_band = max(3W, 6Δx)` — provisional; re-frozen by W-VOF
    characterization and recorded in PHYSICS.md.
  - **mixed_fast** (relaxation): single-phase / weak-coupling only; permitted
    only when ρ ratio ≤ stated limit AND `Ca_spurious` + mass-drift
    validations pass; config validation rejects out-of-range use.
- **NFR-03 (performance)**: Integrate phase field, scalars, forcing into
  fused `step_band`; maintain 3D ring double-buffering and SoA plane-major.
- **NFR-04 (determinism)**: Reductions deterministic-order; GPU/MPI use
  tolerance-based regression (bit-equivalence single-backend only).

---

## 8. Validation / T17 acceptance

Tests authored adversarially from this spec by codex/Opus, separated from
implementation. Each test = metric / target·reference / tolerance /
resolution / time window / backend / pass-fail rule (T17 row format).

**Band governance**: every VR-STR item carries a provisional band from day
one (below = MVP gate). Bands finalized by the standard protocol
(implement → characterize → record rationale in PHYSICS.md → freeze in
VALIDATION.md T17). **Asymmetric rule**: tightening always allowed;
loosening a provisional band requires recorded physical rationale in
PHYSICS.md.

**Provisional bands (MVP gate)**:
- Rushton `Np` vs correlation: ±10%
- PIV/LDA (VR-STR-01): L2_rel < 15%, L∞_rel < 30% per line
- static droplet mass drift: < 0.1% / 1000 steps; advected droplet
  (periodic, one period): < 0.1% (CR-002)
- single-bubble terminal velocity vs Grace (02a): ±10%
- `k_L a` vs correlation (02c): ±25%
- well-balanced static stratification (VR-STR-06): max|u| < 10⁻⁶ (LU) at
  ρ ratio 10³ — retighten after discretisation freeze
- GPU/MPI cross-backend drift (VR-STR-05): mean < 2%, higher-order < 5%
- `Ca_spurious < 10⁻³`

**Mandatory negative / consistency tests**:
- forcing-moment sign negative test: body-force Poiseuille FAILS with
  `Π_force` sign flipped (CR-004).
- stress stage-convention cross-check: pre-collision vs post-collision
  + transform on Couette/Poiseuille/Taylor-Couette (CR-003).
- `J_ρ` consistency code-path + droplet advection conservation (CR-002).
- sparger phase unit test: gas inlet injects φ=0; gas-volume balance
  closes (CR-001).
- scalar total-mass conservation in phase-wise form (MJ-011).
- active-scalar dt-halving convergence (MJ-007).

**VR items**:
- **VR-STR-01** single-phase stirring: baffled tank, specified Re, ungassed.
  Rushton `Np` vs correlation; discharge-velocity vs PIV/LDA. References:
  Wu & Patterson (1989) LDA; Deen et al. (2002) PIV (Rushton D/T=1/3, 4
  baffles); standard `Np` correlations. Bands frozen via T17.
- **VR-STR-02** gas-liquid (02a/b/c): 02a single bubble `U_t` vs Grace
  Eo-Mo-Re; 02b bubble swarm (`ε_g` distribution, hindered rise, `d_32`
  when breakup/coalescence allowed, turbulence intensity when BIT used);
  02c aeration (`ε_g, d_32, k_L a` vs correlations). Refs: 02a = Grace;
  02c = published. In RELAX-PB, `d_32` presupposes PBM; resolved-phasefield
  default measures `d_32` by interface segmentation.
- **VR-STR-03** shear/stress: MMS, curved Couette, rotating cylinder, non-
  Newtonian Poiseuille, multiphase static droplet. Grid convergence order +
  `L2/L∞`. Line design accounts for `L∞` severity near walls.
- **VR-STR-04** scalar/reaction: Taylor-Aris dispersion; reaction-diffusion
  front with known `Da`; `k_L a` (formula = interface integral or correlation,
  explicit). Tolerance, `Pe/Da/Sc`, BCs specified.
- **VR-STR-05** coupled regression / conservation: `probe_state_hash`
  single-backend only. Drift thresholds individually for mass, momentum,
  scalar totals, gas volume, particle count, energy-like quantities.
  Energy-like = monitoring only, not exactly conserved. GPU/MPI tolerance-
  based.
- **VR-STR-06** well-balanced (gravity axis reference landed): `|u| < ε`
  static stratification. **06+**: active scalar ON, `C≡C_0` → same
  quiescence (`F_b^{scalar}` exact-zero degeneration); `∇σ=0` degeneration
  of variable-σ form (agreement with σ-constant reference).
- **VR-STR-07** initialization independence: varying spin-up / statistics-
  start, quasi-steady statistics agree within threshold.
- **VR-STR-RELAX**: relaxation extensions accepted by relative degradation
  vs corresponding fidelity reference. Tolerances via characterize→freeze:
  MRF vs IBM-inertial/overset (`Np`, discharge, mean velocity, torque; steady
  regime only); PB vs resolved-phasefield (`ε_g, d_32, k_L a`, momentum/
  scalar budget; FR-VOF-04 applicability); 1W vs two-way (particle statistics,
  mass-loading limit for neglecting reaction force); AMR vs uniform
  (conserved quantities, interface position, torque, `L2` norm, coarse-fine
  budget); f32 vs fidelity/full_f64 (drift, `Ca_spurious`, `Np`, curvature,
  reductions; single-phase / weak coupling only).

---

## 9. Major technical risks

| # | Difficulty | Risk | Mitigation / status |
|---|---|---|---|
| 1 | High-ρ-ratio two-phase (10³) | Interface instability, spurious currents, f32 rounding | well-balanced phase field + D3Q27 + f64 reductions, point-bubble alternative. **Open — W-VOF pending.** |
| 2 | High-Re stability | Divergence, hyperviscosity, positivity | cumulant + WALE, algebraic `τ_eff`, limiter. **W-LES landed; cumulant landed.** |
| 3 | Rotating-boundary conservation | IBM slip, torque error | multi-direct-forcing, overset reference, thresholding. **W-ROT landed.** |
| 4 | Coupling stiffness | Time-scale divergence | operator-splitting, subcycling, capillary Δt. |
| 5 | Compute cost / memory | 1e9 cells × many distributions = TB scale | memory budget + multi-GPU + MPI + AMR advanced. **Blocked on ME-3 cluster campaign.** |
| 6 | Mass transfer / BIT | Unresolved aeration turbulence, insufficient transfer | Sc_t SGS, BIT production, resolved/point separation. |

---

## 10. Design decisions (finalized)

Each axis default = fidelity-first; low-cost approximations = bolt-on
extensions. Interface default = `resolved-phasefield`; scalar default =
`active`; grid default = `uniform`; precision default = `mixed_safe`;
all-f64 is the reference tier; aggressive f32 (`mixed_fast`) is relaxation.

**Remaining spec-refinement work** (not decisions):
- Active-scalar feedback: concrete equations + stabilization (incl.
  Marangoni). Researched: `docs/proposals/active-scalar-feedback.md`. **One
  derivation is mandatory before implementation** (Marangoni coefficient
  consistency with (κ,β) convention). Thermal axis recommended as API-
  reserved extension.
- Tolerance thresholds vs fidelity reference for each relaxation extension
  (§8 VR). Structure defined as VR-STR-RELAX; numeric bands frozen at
  relaxation implementation time.
- f64/f32 boundary of fidelity profile (band width near interface,
  reduction range). Frozen experimentally during **W-VOF** implementation.
- API definition of trait boundary of each mode axis. Contract requirements
  fixed as FR-EXT-01; concrete Rust API co-designed with SOLVER_IMPROVEMENT_
  SPEC B-1.

**Product-layer scope note**: GUI/CAD, STL import, materials DB, Python/CLI
SDK, sweeps/optimizer, cloud/cluster integration, packaged validation
assets, competitive benchmark tables — **out of scope for this solver
spec**. Managed in separate volumes.

---

## 11. Implementation DAG

Items with no dependency edge implement **concurrently** (parallel-agent
worktrees). Mapping to PLAN.md M-F tracks (MF-α..ζ) noted per row.

| Item | Hard deps | Parallel | Status / notes |
|---|---|---|---|
| W0 core basis (D3Q19/27, cumulant, Guo) | — | — | = MF-α. **LANDED** |
| W-EXT trait contracts (FR-EXT-01) | W0 | yes | co-designed with B-1. Design note landed 2026-07-06 |
| W-UNIT unit/nondim feasibility (§2) | W0 | yes | independent, early |
| W-STRESS stress fields (FR-STRESS) | W0 | yes | ⊂ MF-β. Top priority |
| W-ROT rotating IBM (FR-ROT-01) | W0 | yes | = MF-δ. **LANDED** |
| W-GRAV well-balanced gravity (FR-BC-02) | W0 | yes | ⊂ MF-γ. **LANDED** |
| W-SCAL passive scalar ADE | W0 | yes | ⊂ MF-ε. SGS part waits on W-LES |
| W-LES turbulence SGS (FR-LES) | W-STRESS | conditional | ⊂ MF-β. **LANDED (WALE + GPU)** |
| **W-VOF resolved interface (FR-VOF-01/02)** | **W-GRAV** | conditional | fidelity default; hardest item; **CRITICAL PATH; PENDING**; ⊂ MF-γ |
| W-PART particles + exposure | W-STRESS (SGS: W-LES) | conditional | ⊂ MF-ε. **D-track P2 LANDED** |
| W-REACT reaction / active feedback | W-SCAL (active: W-VOF) | conditional | ⊂ MF-ε |
| W-BUB point bubbles + PBM + interfacial transfer | W0, W-SCAL, W-EXT | conditional | Phase 2 (API-reserved) |
| W-BCTOP top boundary / degassing / contact angle | W-VOF | conditional | ⊂ MF-γ. **Waits on W-VOF** |
| W-COUP coupling loop (FR-COUP) | active set | incremental | ⊂ MF-ζ |
| W-IO I/O & analysis (FR-IO) | each producing subsystem | incremental | ⊂ MF-ζ |
| W-VAL validation T17 | each subsystem | yes | codex adversarial authorship, separated from implementation |

**Parallel waves**:
1. After W0 (LANDED): W-EXT / W-UNIT / W-STRESS / W-ROT / W-GRAV / W-SCAL —
   6-way parallel. Of these, W-ROT + W-GRAV LANDED.
2. After deps: W-LES (←STRESS) / **W-VOF (←GRAV) — currently blocking** /
   W-PART (←STRESS) / W-REACT (←SCAL).
3. Later: W-BCTOP (←VOF) / W-BUB (←SCAL, EXT) / active feedback +
   interfacial transfer (←VOF).
4. Cross-cutting: W-COUP / W-IO / W-VAL.

**Critical path**:
- `W0 → W-GRAV → W-VOF → W-BCTOP / interfacial transfer` — **stuck at
  W-VOF**.
- `W0 → W-STRESS → W-LES → W-PART` — **LANDED**.

Throughput/scaling KPIs stay delegated to CLUSTER_OPTIONS.md (R3) — not
duplicated here; no hardcoded numeric thresholds (bands freeze via T17);
product ecosystem in separate volumes.
