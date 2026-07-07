# docs/qa/multiphase-hard-cases.md — Hard multiphase validation cases with known answers

Status: research deliverable (2026-07-07). Scope: cases **beyond** the in-flight set (T11 statics, T12 RT growth, capillary-wave dispersion, Lamb droplet, two-layer Couette, Lucas–Washburn, Hysing/Grace gated on phase-field).

## Available knob sets (repo-grounded)

| ID | Knobs | Where |
|----|-------|-------|
| **K1** — SCMP 2D | `G`, `Psi` (Classic/Exponential), `g_wall`, `wall_rho`, per-cell `force_field_mut`, `set_gravity` | `crates/lbm-core/src/compat/multiphase.rs` (`ShanChen`), `compat/sim.rs` |
| **K2** — SCMP 3D | `Solver::update_shan_chen_force_with_walls(g, g_wall, psi_wall, psi)` — dimension-generic (D3Q19), halo-aware; body-force field + gravity; moving-wall BC exists (T15.5 cavity) | `crates/lbm-core/src/solver.rs` ~L2264 |
| **K3** — MCMP 2D | `MultiComponent`: `G_ab` (ψ=ρ), per-component gravity `g_a`/`g_b`, per-component wall affinity `g_wall_a/b`; **per-component ν** (each component is its own `Simulation`) | `compat/multiphase.rs` |
| Missing | 3D MCMP (no cross-component force in V2 `Solver`); variable σ / thermal field; tunable-EOS ψ beyond the two forms | — |

Global constraints that shape every acceptance band below: coexistence at G=−5 is ρ_l/ρ_v = 1.888/0.1194 ≈ **15.8** (T11 frozen), σ ≈ 3.32e-2 (T11), spurious currents ≤ 5e-3 (velocity noise floor), Ma limit |u| ≤ ~0.1 (ceiling). Signal velocities must live in the ~1-decade window between the two. SCMP shares ν across phases ⇒ dynamic viscosity ratio = density ratio ≈ 16 (bubble-like, λ ≈ 0.063). σ and ρ_l,v are **emergent** — always use the T11-measured values, never nominal inputs (zero-free-parameter cross-prediction is the whole point of the top-ranked cases).

---

## Ranked cases (information-per-compute among TODAY-feasible)

### 1. Jurin capillary rise statics — **YES** (K1: `wall_rho` + gravity)
- **Hard because**: it is the first test that couples all three separately-frozen calibrations — σ (T11 Laplace), θ (T11c wall_rho), and gravity — with **no free parameters**. SC failure modes it exposes: vapor condensation in narrow slots, barometric vapor stratification under g, contact-line spurious currents that never settle, vapor weight neglected in naive h formulas.
- **Known answer**: parallel-plate slot of gap w: h = 2σ cos θ / (Δρ g w), Δρ = ρ_l − ρ_v. Textbook (Jurin 1718; de Gennes, Brochard-Wyart & Quéré, *Capillarity and Wetting Phenomena*, ch. 2).
- **Acceptance**: with σ from T11 and θ from T11c taken as inputs, measured meniscus height within ±10% for ≥3 gaps w ∈ {16, 24, 32}; linearity h vs 1/w with R² ≥ 0.99; run two wall_rho values (θ ≈ 63°, ~107° → capillary depression!) — sign flip of h is a free adversarial check.
- **Compute**: tiny (≤ 256×256, 2D). Highest info/compute on the list.

### 2. Prosperetti damped interfacial standing wave + gravity–capillary dispersion crossover — **YES** (K3: G_ab + per-component gravity; T12 rig, stable side)
- **Hard because**: exact **viscous, transient** solution — tests frequency AND damping AND the initial-value transient. Diffuse-interface excess dissipation and interface-thickness artifacts show up here before anywhere else. This is the standard hard benchmark used by Gerris/Basilisk.
- **Known answer**: [Prosperetti, "Motion of two superposed viscous fluids", Phys. Fluids **24**, 1217 (1981)](https://pubs.aip.org/aip/pfl/article-abstract/24/7/1217/437190/Motion-of-two-superposed-viscous-fluids) — exact IVP solution for the interface amplitude a(t), valid when the two fluids have **equal kinematic viscosity** — exactly the T12 configuration (both bulk ρ=1, shared ν). Inviscid frequency ω₀² = (g_eff Δρ_eff k + σ_AB k³)/(ρ₁+ρ₂).
- **Acceptance**: (a) fitted ω within 5% and damping rate within 10% of the exact solution at 2–3 k values; (b) **k-scan across the crossover** k_σ = √(g Δρ/σ): ω(k) tracks gk-dominated → σk³-dominated branches within 10% (this absorbs the "gravity-capillary crossover" candidate); (c) waveform L2 error vs analytic ≤ 10% over 3 periods.
- **Compute**: small 2D (256², short runs). Reuses T12 measurement code (Fourier mode projection).

### 3. Rayleigh–Taylor surface-tension cutoff scan — **YES** (K3; direct extension of T12)
- **Hard because**: T12 checks one growth rate; this checks the **stability boundary** — the sign structure of γ²(k) = A g k − σ_AB k³/(ρ₁+ρ₂) (+ viscous correction already in the T12 γ_th formula). A wrong σ discretization or a spurious-current-driven mixing layer moves k_c.
- **Known answer**: cutoff k_c = √(g Δρ/σ); most-unstable k_max = k_c/√3 (inviscid). Chandrasekhar, *Hydrodynamic and Hydromagnetic Stability* (1961), ch. X (textbook).
- **Acceptance**: scan mode number in a fixed box (or box width at fixed mode): growth for k < k_c, damped oscillation for k > k_c, transition bracketed within ±1 mode of prediction; γ(k) within the T12 band [0.75, 1.25] on the unstable side.
- **Compute**: ~6 T12-sized runs. Nearly free given T12 infrastructure.

### 4. Taylor–Culick film retraction — **YES** (K1/K2, SCMP)
- **Hard because**: exact momentum-balance speed with a singular rim; tests whether surface energy → kinetic energy conversion is quantitatively right (any spurious interface dissipation gives systematically slow rims). Almost never used in LBM validation despite being exact — high adversarial value.
- **Known answer**: v_TC = √(2σ/(ρ_l h)) for film thickness h (Taylor, Proc. R. Soc. A 253, 313 (1959); Culick, J. Appl. Phys. 31, 1128 (1960)). Textbook-established.
- **Acceptance**: 2D liquid strip in vapor, one edge released: steady rim speed within ±10% of v_TC for h ∈ {16, 24, 32}; fitted exponent of v vs h = −0.5 ± 0.05. Caveats to control: SC evaporation mass loss from the thin film (monitor film mass), vapor drag (small at ratio 16 but measurable — report).
- **Compute**: small 2D. 

### 5. Kelvin–Helmholtz marginal stability with surface tension — **PARTIAL** (K3: G_ab + per-component gravity + initial shear profile)
- **Hard because**: analytic **threshold**, not just a rate — bisecting an instability onset is a merciless test of both σ and the stratification force coupling. Main honesty caveat: the analytic result is for a vortex sheet; the initialized shear layer has finite thickness δ growing as √(νt), so the measurement must complete while kδ ≪ 1.
- **Known answer**: for fixed k: instability iff (ρ₁ρ₂/(ρ₁+ρ₂)) k (ΔU)² > g Δρ + σk²; minimizing over k: (ΔU_c)² = 2(ρ₁+ρ₂)/(ρ₁ρ₂) · √(g Δρ σ) at k_c = √(g Δρ/σ). Chandrasekhar 1961 §101 (textbook; also [MIT 1.63 notes](https://web.mit.edu/1.63/www/Lec-notes/chap5_instability/5-2KHdiscont.pdf)).
- **Acceptance**: at fixed box k, bisect ΔU: onset (mode grows ≥3× before kδ > 0.3) brackets analytic ΔU_c within ±15%; below threshold, mode oscillates at the Doppler-shifted frequency.
- **Compute**: ~8–10 small 2D runs per k. 

### 6. Droplet coalescence neck growth — inertial regime **YES**, viscous regime **GATED-resolution** (K1 2D exponent / K2 3D quantitative)
- **Hard because**: near-singular initial dynamics; the neck must traverse decades of scale; diffuse interface sets a hard floor (neck meaningful only once ≫ 4–5 cell interface width). Discriminates interface-topology handling quality.
- **Known answer**: inertial regime r_neck = C (σR/ρ)^{1/4} t^{1/2}, C ≈ 1.4–1.6 ([Duchemin, Eggers & Josserand, JFM 487, 167 (2003)](https://arxiv.org/abs/physics/0212075); crossover physics: [Paulsen, Burton & Nagel, PRL 106, 114501 (2011)](https://arxiv.org/abs/1012.1298), Paulsen et al. PNAS 109, 6859 (2012)). Viscous regime r ~ (σ/πμ) t with log correction (Eggers, Lister & Stone, JFM 401, 293 (1999)) — needs neck resolution below interface width at our sizes → gated. 2D inertial exponent 1/2 also expected but prefactor is **community folklore** — exponent-only acceptance in 2D.
- **Acceptance**: 3D SCMP (two touching R=32 droplets, 128×128×256): fitted exponent 0.50 ± 0.05 over ≥1 decade of t once r > 6 interface widths; prefactor reported (informational band 1.2–1.8, not gating — outer vapor at ratio 16 is not vacuum).
- **Compute**: one medium 3D run (+1 2D cheap screen).

### 7. Plateau–Rayleigh 3D thread breakup — **YES** single-fluid viscous theory (K2); Tomotika two-fluid version **GATED (needs 3D MCMP)**
- **Hard because**: full 3D interface pinch-off; mode selection tests the σk³ physics anisotropy on the lattice; satellite formation tests topology change.
- **Known answer**: inviscid: ω² = (σ/ρR³) · x I₁(x)/I₀(x) · (1−x²), x = kR; max at kR ≈ 0.697, ω_max ≈ 0.34 √(σ/ρR³) (Rayleigh 1878; see [Eggers & Villermaux, Rep. Prog. Phys. 71, 036601 (2008)](https://arxiv.org/pdf/1701.06157) for the viscous generalization). At τ=1, Oh ≈ 0.25 for R≈24 — **use the viscous single-fluid dispersion (Rayleigh 1892/Weber 1931 form in Eggers & Villermaux), not the inviscid one**; outer vapor treated as negligible (μ_v/μ_l ≈ 0.06 — report as systematic). [Tomotika 1935 (Proc. R. Soc. A 150, 322)](https://makingscience.royalsociety.org/items/rr_57_29/referees-report-by-leonard-bairstow-on-a-paper-on-the-instability-of-a-cylindrical-thread-of-a-viscous-liquid-surrounded-by-another-viscous-fluid-by-s-tomotika) gives the arbitrary-viscosity-ratio dispersion with dominant-wavelength selection — the right target the day 3D MCMP lands.
- **Acceptance**: seeded single mode: γ_fit/γ_th(viscous) ∈ [0.85, 1.15] at 3 values of kR straddling 0.697; white-noise seed: dominant emergent kR within ±15% of the viscous-theory maximum; breakup produces main + satellite drops (qualitative gate).
- **Compute**: 64×64×512-ish 3D, moderate.

### 8. Lock-exchange gravity current front speed — **PARTIAL** (K3: per-component gravity; σ must be pushed toward miscible limit)
- **Hard because**: large-scale buoyancy-inertia balance, front condition is a nonlinear integral constraint — tests the per-component gravity coupling at finite amplitude, far from the linear-instability cases above. Caveat: theory is for miscible fluids; MCMP has σ_AB → require Bo = Δρ_eff g H²/σ_AB ≳ 10² (large H, G_ab near separation threshold), report σ correction.
- **Known answer**: energy-conserving front: u_f = ½ √(g′H) (Fr_H = 1/2, equivalently 1/√2 on current depth h = H/2) — Benjamin, JFM 31, 209 (1968); DNS reference values: [Härtel, Meiburg & Necker, JFM 418, 189 (2000)](https://www.researchgate.net/publication/30011744_High_resolution_numerical_simulations_of_lock-exchange_gravity-driven_flows) (no-slip fronts a few % below ½).
- **Acceptance**: constant-velocity slumping-phase front speed within ±10% of ½√(g′H) at the largest feasible Bo and Re ≥ ~10³; monotone approach to the bound with increasing Bo.
- **Compute**: one long-box 2D run per (Bo, Re) — moderate.

### 9. Bubble deformation in shear (Taylor 1934, λ ≈ 0.06 branch) — **PARTIAL** (K2: 3D SCMP + moving-wall Couette)
- **Hard because**: quantitative viscous-stress vs σ balance in 3D with an imposed outer flow; orientation (45°) and the λ-dependence are free extra checks. SCMP fixes λ = μ_v/μ_l ≈ 0.063 (vapor bubble in sheared liquid — use the **liquid** viscosity in Ca). PARTIAL because wall confinement corrections (Shapira & Haber 1990) must be applied or walls kept ≥ 4a away, and Ca window is squeezed by the spurious-velocity floor.
- **Known answer**: D ≡ (L−B)/(L+B) = Ca · (19λ+16)/(16λ+16) ≈ 1.01·Ca at λ=0.063, small Ca ([Taylor, Proc. R. Soc. A 146, 501 (1934)](https://arxiv.org/html/2407.10880); formula ubiquitously reproduced). Do **not** run this in 2D — Taylor's formula is 3D; 2D counterparts in LBM papers are apples-to-oranges.
- **Acceptance**: D vs Ca slope within ±15% for Ca ∈ {0.05, 0.1, 0.2}; major-axis angle 45° ± 5° at smallest Ca.
- **Compute**: ~128³ × 3 runs — significant but standard.

### 10. Bretherton film deposition — **PARTIAL** (K3: 2D channel, body-force-driven; resolution-limited)
- **Hard because**: THE thin-film/lubrication benchmark — the deposited film is set by a delicate front-region curvature matching. Diffuse interface must be ≪ film ≪ gap: at 1.34 Ca^{2/3} with film ≥ 6 cells and ≥ 3 interface widths you need half-gap ≥ ~200 cells at Ca ≈ 0.05. Honest expectation: exponent yes, coefficient strained.
- **Known answer**: h∞/r_meniscus = 1.34 Ca^{2/3}, Ca = μ_l U_b/σ, small Ca ([Bretherton, JFM 10, 166 (1961)](https://pubs.aip.org/aip/pof/article/33/12/123303/1062080/On-the-extension-of-Bretherton-theory-for-thin)); planar-channel form and moderate-Ca extensions collected in [de Lózar et al. / Aussillous–Quéré lineage and planar-vs-axisymmetric comparison](https://arxiv.org/pdf/1711.10447).
- **Acceptance**: fitted exponent 0.67 ± 0.07 over Ca ∈ [0.02, 0.2]; coefficient within ±25% (informational at the low end); bubble speed vs mean flow relation (U_b/U_mean − 1 ~ Ca^{2/3}) as a second observable.
- **Compute**: long 2D channel at high resolution — expensive for 2D.

### 11. Taylor bubble rise / Davies–Taylor Froude number — **PARTIAL** (K2: 3D SCMP + gravity + voxel tube)
- **Hard because**: integrates gravity, curvature-dominated nose, draining wall film, and wall wetting in one steady observable. Insensitive to gas density (good for ratio 16) but the wall film needs resolution; σ and μ corrections push you off the ideal 0.351 unless Eo ≳ 40, and reaching Eo ≳ 40 at fixed σ ≈ 0.033 means large D and small g — long transients.
- **Known answer**: Fr = U/√(gD) = 0.351 (Dumitrescu 1943 — regarded as the accurate value; Davies & Taylor 1950 got 0.328) for inertial regime; full (Eo, Mo) correlation: [Viana et al., JFM 494, 379 (2003)](https://dept.aem.umn.edu/~./faculty/joseph/archive/docs/332_taylor-bubble-jfm-mar17-nov25.pdf); regime boundaries White & Beardmore (1962).
- **Acceptance**: Fr within ±10% of the Viana correlation evaluated at the **simulated** Eo, Mo (not the ideal 0.351); nose shape spherical-cap fit radius vs Dumitrescu profile (qualitative).
- **Compute**: D ≥ 64 cells, L ≥ 6D, 3D, long — the most expensive YES-track item so far.

### 12. Hadamard–Rybczynski terminal velocity (+ Hasimoto periodic correction) — **PARTIAL** (K2: 3D SCMP + gravity)
- **Hard because**: uniquely tests **interfacial momentum transfer / internal circulation** — a clean rigid-sphere-vs-fluid-sphere 3/2 velocity discriminant. Killers: Re < 1 needs tiny g ⇒ terminal velocity dangerously close to the spurious-current floor; periodic-image hindrance is large.
- **Known answer**: U = (2Δρ g a²)/(9μ_c) · 3(1+λ)/(2+3λ); λ→0 bubble limit: 1.5× the rigid Stokes velocity (Hadamard 1911; Rybczynski 1911; Clift, Grace & Weber, *Bubbles, Drops and Particles*, 1978, ch. 3). Periodic-box correction: Hasimoto, JFM 5, 317 (1959) — apply, don't hand-wave.
- **Acceptance**: Hasimoto-corrected U within ±15% at Re ≤ 0.5 for two box sizes (consistency between the two after correction is the real gate); internal circulation pattern present (streamlines in the drop frame — behavior-validity review item).
- **Compute**: 3D, slow settling — expensive per data point.

### 13. Spinodal decomposition domain-growth exponents — **PARTIAL** (K1/K3 2D screen; 3D for the real thing)
- **Hard because**: statistical, multi-decade-in-time scaling with slow crossovers; tests coarsening hydrodynamics no other case touches. But: prefactors non-universal, crossovers broad — low information density per CPU-hour, and 2D scaling is known to **break** (Wagner & Yeomans, PRL 80, 1429 (1998): 2D binary-fluid scale-invariance breakdown), so 2D is a qualitative screen only.
- **Known answer**: 3D symmetric binary fluid: diffusive L ~ t^{1/3} (Lifshitz–Slyozov), viscous hydrodynamic L ~ t (Siggia, PRA 20, 595 (1979)), inertial L ~ t^{2/3} (Furukawa); definitive LB study with regime windows: [Kendon, Cates, Pagonabarraga, Desplat & Bladon, JFM 440, 147 (2001)](https://arxiv.org/abs/cond-mat/0006026).
- **Acceptance**: local exponent d ln L/d ln t within ±0.1 of the regime value over ≥1 decade, regime identified a priori from (L/L₀, t/t₀) reduced units per Kendon et al. — not post-hoc regime shopping.
- **Compute**: long runs; 3D version properly expensive. Rank low despite good theory.

### 14. Bhaga–Weber rising-bubble shape/Re points — **PARTIAL** (K2 3D; density-ratio caveat)
- **Hard because**: finite-deformation gravity-σ-viscosity balance; the closed-wake regimes are steady and photogenic (good behavior-review material). PARTIAL: experiments are at density ratio ~10⁵; at 16 the internal density is not negligible — stick to Mo > 4e-3, Re < 110 where drag is outer-controlled.
- **Known answer**: [Bhaga & Weber, JFM 105, 61 (1981)](https://www.semanticscholar.org/paper/Bubbles-in-viscous-liquids:-shapes,-wakes-and-Bhaga-Weber/ea0ab083c9d7ddb889b42e84cb0d2ce17874fcc3) — Re(Eo, Mo) correlations + shape photographs; complements the gated Hysing/Grace track with 3 concrete (Eo, Mo) points.
- **Acceptance**: terminal Re within ±15% of the B-W correlation; aspect ratio within ±10%; shape class matches the photograph (oblate ellipsoid / dimpled cap).
- **Compute**: 3D, 2–3 points. Expensive; partially redundant with the phase-field-gated Grace work — run at most as a bridge.

### 15. Liquid-bridge Plateau stability limit — **PARTIAL** (K2: 3D SCMP between solid disks)
- **Hard because**: bifurcation-point detection with contact-line pinning at disk edges (Gibbs pinning via bounce-back geometry is only approximately sharp on a voxel lattice).
- **Known answer**: cylindrical bridge between coaxial disks unstable iff L > 2πR (Plateau limit; slenderness analysis Gillette & Dyson, Chem. Eng. J. 2, 44 (1971); textbook in e.g. Langbein, *Capillary Surfaces*). 
- **Acceptance**: bisect L/2πR ∈ [0.8, 1.2]: breakup onset at 1.0 ± 0.1; below limit, bridge relaxes to stable equal-volume shape.
- **Compute**: 3D moderate. Pinning fidelity risk → rank low.

### 16. Cox–Voinov dynamic contact angle / Tanner spreading — **PARTIAL** (K1: wall_rho + forced spreading)
- **Hard because**: contact-line motion is regularized by the diffuse interface (effective slip length is an emergent property, not an input) — this case *characterizes* a closure rather than validating first-principles physics. Per physics-discipline: the fitted log constant must be recorded as a closure with its measured value ~ interface width.
- **Known answer**: θ_d³ = θ_e³ + 9 Ca ln(x/L_s) (Voinov 1976; Cox, JFM 168, 169 (1986)); complete-wetting spreading R ~ t^{1/10} (Tanner 1979). Both textbook (de Gennes et al., *Capillarity and Wetting*).
- **Acceptance**: θ_d³ − θ_e³ linear in Ca with slope 9·ln(x/L_s), L_s fitted once and required to be O(interface width) and constant across θ_e; Tanner exponent 0.10 ± 0.02. Directly de-risks the in-flight Lucas–Washburn case (same contact-line closure).
- **Compute**: small 2D. Ranked low only because one constant is fitted.

### 17. Saffman–Taylor fingering — **PARTIAL, needs a harness-level Hele-Shaw closure** (K3 + per-cell force API)
- **Hard because / honest caveat**: the classic results (finger width → ½ channel at low σ-parameter, McLean & Saffman JFM 102, 455 (1981); Chuoke linear dispersion) live in **Darcy/Hele-Shaw**, not 2D Navier–Stokes — a direct 2D MCMP run does NOT converge to them. Feasible route: add a depth-averaged Brinkman drag −(12ν/b²)u via the existing per-cell force field each step (literature-backed closure — must go through lbmflow-physics-discipline provenance, with its own Poiseuille-in-Hele-Shaw sub-validation). Viscosity contrast itself is available (per-component ν in K3).
- **Acceptance (if closure lands)**: Chuoke most-unstable wavelength within ±15%; single-finger relative width vs the McLean–Saffman curve at 2 values of the control parameter.
- Rank last among feasible: real physics value, but gated on a new (if small and legitimate) closure.

---

## Assessed and NOT recommended today

| Case | Verdict | Reason |
|------|---------|--------|
| **Tomotika two-fluid thread breakup** | GATED — 3D MCMP | The dispersion relation (viscosity-ratio-dependent dominant wavelength, Tomotika 1935) is the best-in-class target, but `MultiComponent` is 2D-only. Queue behind a V2 cross-component force. |
| **Rayleigh–Plesset / Minnaert bubble oscillation** | NOT TODAY (honest) | The SC vapor obeys the SC EOS near its critical point, not a polytropic ideal gas; Minnaert's ω² = 3γp₀/(ρR²) does not apply. One could derive an SC-EOS analog, but that validates the code against itself (no external truth). Revisit with a tunable-EOS ψ (e.g. Carnahan–Starling) or free-energy model. |
| **Sessile-drop evaporation (d²-law, constant-angle modes)** | N/A | Needs thermal/vapor-diffusion control; isothermal SC evaporation has no clean external reference. |
| **Marangoni (thermocapillary migration, Young et al. 1959)** | N/A | No temperature field and no spatially variable G/σ anywhere in the stack (checked: `Psi` fixed forms, scalar `g`). Note for the phase-field track: variable-σ knob would unlock the exact YGB migration velocity — best-in-class when available. |
| **Grace-diagram sweep / Hysing benchmark** | Already gated in-flight | Do not duplicate; case 14 (Bhaga–Weber points) is the SCMP-feasible bridge. |
| **Partial-coalescence cascade (daughter ratio ~0.5)** | Watch list | Blanchette & Bigioni (Nat. Phys. 2006) reference exists but bands are loose — community folklore territory for acceptance purposes. |

## Cross-cutting design rules
1. Every case uses **measured** σ (T11 Laplace) and θ (T11c) as inputs — cases 1–5 then have zero free parameters; a failure is unambiguous.
2. Signal velocity window: 5e-3 (spurious floor) ≪ u ≪ 0.1 (Ma ceiling). Cases that can't fit a decade of signal in that window (very-low-Re HR, very-low-Ca Bretherton) get PARTIAL for that reason alone.
3. Per the behavior-validity directive: each case's acceptance includes at least one **pattern** check (orientation angle, shape class, internal circulation, satellite drops), not just the scalar band.

---

**Report summary**: Deliverable above, grounded in the repo (SCMP is available in 3D via `Solver::update_shan_chen_force_with_walls` in `/Users/taku/projects/流体シミュレータ/crates/lbm-core/src/solver.rs`; MCMP with per-component ν/gravity/wall-affinity is 2D-only in `/Users/taku/projects/流体シミュレータ/crates/lbm-core/src/compat/multiphase.rs`; frozen density ratio is ≈15.8 per `docs/VALIDATION.md` T11). 17 new cases ranked; top-4 (Jurin statics, Prosperetti standing wave, RT σ-cutoff, Taylor–Culick) are cheap 2D runs with exact answers and zero free parameters. Two honest rejections (Rayleigh–Plesset, Marangoni) with the exact missing knob named. Not written to disk (read-only grounding specified) — commit body is the markdown between the `---` markers.

Sources: [Duchemin, Eggers & Josserand — Inviscid coalescence of drops](https://arxiv.org/abs/physics/0212075) · [Paulsen, Burton & Nagel — Viscous-to-inertial crossover](https://arxiv.org/abs/1012.1298) · [Paulsen et al. PNAS 2012](https://pnas.org/content/early/2012/04/16/1120775109) · [Eggers & Villermaux review (jet/thread dispersion)](https://arxiv.org/pdf/1701.06157) · [Tomotika 1935 (Royal Society record)](https://makingscience.royalsociety.org/items/rr_57_29/referees-report-by-leonard-bairstow-on-a-paper-on-the-instability-of-a-cylindrical-thread-of-a-viscous-liquid-surrounded-by-another-viscous-fluid-by-s-tomotika) · [Prosperetti 1981, Phys. Fluids 24, 1217](https://pubs.aip.org/aip/pfl/article-abstract/24/7/1217/437190/Motion-of-two-superposed-viscous-fluids) · [MIT 1.63 KH notes](https://web.mit.edu/1.63/www/Lec-notes/chap5_instability/5-2KHdiscont.pdf) · [Bretherton extensions, Phys. Fluids 33, 123303](https://pubs.aip.org/aip/pof/article/33/12/123303/1062080/On-the-extension-of-Bretherton-theory-for-thin) · [Planar vs axisymmetric Taylor films](https://arxiv.org/pdf/1711.10447) · [Härtel, Meiburg & Necker 2000 lock-exchange DNS](https://www.researchgate.net/publication/30011744_High_resolution_numerical_simulations_of_lock-exchange_gravity-driven_flows) · [Kendon et al. JFM 440 (2001) spinodal LB study](https://arxiv.org/abs/cond-mat/0006026) · [Wagner & Yeomans 2D spinodal](https://www.researchgate.net/publication/11150360_Lattice_Boltzmann_study_of_spinodal_decomposition_in_two_dimensions) · [Bhaga & Weber JFM 105 (1981)](https://www.semanticscholar.org/paper/Bubbles-in-viscous-liquids:-shapes,-wakes-and-Bhaga-Weber/ea0ab083c9d7ddb889b42e84cb0d2ce17874fcc3) · [Viana et al. Taylor-bubble correlation](https://dept.aem.umn.edu/~./faculty/joseph/archive/docs/332_taylor-bubble-jfm-mar17-nov25.pdf) · [Taylor 1934 deformation formula (review)](https://arxiv.org/html/2407.10880) · [Rayleigh–Plateau parameters](https://www.pnas.org/doi/10.1073/pnas.2306088120)