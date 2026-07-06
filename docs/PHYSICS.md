# PHYSICS.md — record of physics models and numerical experiments

Record here the experiments and findings that justify fixes to the specification (VALIDATION.md).

## Adopted models (as of Phase 1)

- **Lattice**: D2Q9, cs² = 1/3, τ = 3ν + 0.5
- **Collision**: BGK / TRT. For TRT, ω⁺ determines viscosity and Λ = (1/ω⁺−½)(1/ω⁻−½) is fixed
  (default Λ = 3/16 → the half-way placement of a straight wall is exact for parabolic flow)
- **Body force**: Guo forcing (2nd-order accurate). The physical velocity u = (Σf c + F/2)/ρ is
  used in feq, the force term, and the output — all of them
- **Walls**: half-way bounce-back. Moving walls use the momentum-injection term +6 w_q ρ (c_q·u_w)
- **Open boundaries**: a single implementation of Zou-He parameterized by the face normal (n, t) (common to all 4 edges):
  - ρ = (S0 + 2S⁻)/(1 − u·n) (when velocity is specified)
  - f_n = f_{−n} + (2/3)ρ(u·n)
  - f_{n±t} = f_{−n∓t} + (1/6)ρ(u·n) ± [½ρ(u·t) − ½T], T = f_{+t} − f_{−t}
- **Force measurement**: momentum-exchange method (F_body = Σ_links −c_q(f_out + f_in))

## Experiment records

### 2026-07-06: WALE LES MF-beta subset
- Adopted LES model: **WALE** (Nicoud & Ducros 1999), not Smagorinsky. The decisive property is the laminar
  null behavior: for unidirectional pure shear (Couette and Poiseuille), `S^d:S^d = 0` analytically, so
  `nu_t = 0`. Smagorinsky would produce nonzero eddy viscosity in the same resolved laminar shear and would
  silently bend the baseline physics.
- Formula and conventions:
  - `g_ij = du_i/dx_j`
  - `S_ij = (g_ij + g_ji)/2`
  - `S^d_ij = (g_ik g_kj + g_jk g_ki)/2 - delta_ij tr(g^2)/3`
  - `nu_t = (Cw Delta)^2 (S^d:S^d)^(3/2) / ((S:S)^(5/2) + (S^d:S^d)^(5/4))`
  - `Cw = 0.325`, `Delta = 1` in lattice units; the `0/0` limit is defined as `nu_t = 0`.
- Implementation convention: WALE is a solver-level driver that writes a compact per-cell `omega_plus = 1/tau_eff`
  field, where `tau_eff = 3(nu0 + nu_t) + 0.5`. Collision kernels only replace the local `omega_plus` fetch when
  the optional field is present; the field-off path keeps the original uniform-relaxation arithmetic. The SGS field
  is computed from the current post-streaming moments and therefore applies with a one-step lag.
- Velocity-gradient observable: the symmetric off-diagonal shear terms reuse the native non-equilibrium stress
  path used by `gather_strain_rate`. The diagonal entries and the antisymmetric rotation are reconstructed from
  velocity differences; wall-adjacent derivatives use the half-way wall location and stored wall velocity. This
  mixed path is intentional: D3Q19 moving-wall-adjacent non-equilibrium normal stresses showed small diagonal
  artifacts in pure Couette, while the velocity derivative gives the correct pure-shear null tensor.
- Measured release-test results in `crates/lbm-core/tests/wale_les.rs`:
  - Constant per-cell `omega_plus` equal to the global value is bitwise identical to the field-off path for both
    `CpuScalar` and `CpuSimd`.
  - Couette null property: max `nu_t = 2.65e-48` (gate `<= 1e-12`).
  - Poiseuille null property: max `nu_t <= 1e-12` (same test gate).
  - Laminar duct non-interference after a null WALE update: velocity-field `L_inf <= 1e-12`.
- Heavy characterization still pending in ignored tests: N=64 T15.4 TGV fitted `nu_eff`, deterministic N=48
  multimode stabilization, and Re_tau=180 channel DNS comparison (T17/VR-STR-03).

### 2026-07-04: Level of mass drift due to rounding error
- Total mass drift 1.05e-13 (relative) in a periodic box 64² over 1000 steps.
- Collision and streaming conserve exactly analytically → this is accumulation of f64 rounding.
- **Specification change**: set the T6 mass-conservation tolerance to 1e-12 (10³ step) / 1e-11 (10⁴ step).

### 2026-07-04: The BGK steady state oscillates at a rounding-error plateau
- Poiseuille 4×10, BGK, τ=0.8: after reaching physical steady state at ~8500 steps,
  the step-to-step difference dmax/umax **oscillates permanently at ~1e-12** (does not reach 1e-13).
- TRT reaches an exact discrete fixed point (the difference truly becomes 0), so it passes even at 1e-13.
- **Specification change**: set the recommended ε for steady-state judgment to 1e-11 (replacing T2's "ε=1e-13").

### 2026-07-04: rayon dispatch overhead on small lattices
- 203 µs/step on a 4×10 lattice (18 cores, parallel on). Almost all of it is rayon's
  task-distribution cost. Serial would be on the order of ~0.1 µs/step.
- **Implementation change**: automatically fall back to serial execution when the cell count < 16384
  (`PARALLEL_MIN_CELLS`).

### 2026-07-04: TGV requires pressure-consistent initialization (the acoustic-wave residual is O(u0))
- If density is initialized uniformly as ρ=1, the inconsistency with the analytic pressure field radiates as acoustic waves and,
  failing to decay fully, contaminates the velocity field. Measured: error ≈ 0.30/N − 0.7/N² (1st order dominates).
- Initializing with ρ = 1 − (3u0²/4)(cos 2kx + cos 2ky) gives:
  - e32=2.62e-3, e64=6.98e-4, e128=1.78e-4 / convergence order 1.91, 1.98 (clean 2nd order)
- Because the missing f⁽¹⁾ of equilibrium initialization is also an O(1/N) contamination source, `init_with` was
  specified to always add the finite-difference Chapman-Enskog non-equilibrium term.
- **Specification change**: state pressure-consistent initialization and diffusive scaling (u0 = 1.28/N) in T1.

### 2026-07-04: The Zou-He pressure boundary had a sign bug in the normal velocity (detected by self-review)
- From the closure relation ρ(1 − u·n) = S0 + 2S⁻, un = 1 − (S0+2S⁻)/ρ is correct, but
  the implementation had un = (S0+2S⁻)/ρ − 1 with the sign reversed. The velocity-boundary side was correct.
- After the fix, in a pressure-difference-driven channel (Δρ=2e-3, H=32, 20k step),
  u_center/u_theory = 0.9974 (0.26% agreement with the Poiseuille analytic solution, spec ±2%) was confirmed.
- Lesson: open boundaries require analytic cross-checking in all 4 directions + both BC kinds (codex suite T4/T5 is a permanent guard).

### 2026-07-05: Triaged the 4 items of the first codex adversarial-test batch (2 spec bugs, 1 spec ambiguity, 1 f32 characteristic)
1. **Staggered boundary layer of Zou-He pressure outflow**: alternating oscillation in the ~4 columns just before the outflow (at the boundary node
   ±2%, decay length ~4 cells). Exactly matches for TRT(3/16), TRT(1/4), BGK → judged to be a
   boundary-condition-specific O(Ma²) artifact independent of the collision operator. The bulk (24 or more columns inside) mass flux is
   constant at 2.4e-5, total mass drift 2e-13. → T4 changed to "bulk constancy 1e-4".
   Improvement candidates (Phase 7 backlog): characteristic BC / anti-bounce-back outflow.
2. **Exact anti-symmetry of a simple Δρ reversal is a physically incorrect spec**: the inertial term and compressibility are 2nd order in u →
   anti-symmetry holds only up to O(Ma²) (measured 1.7e-3 relative). The exact angle is
   replaced with "Δρ reversal + x mirror = exact mirror match" (discrete symmetry).
3. **The f32 uniform-field momentum error is a coherent rounding bias**: identical operations on all cells →
   ~1ulp/step accumulates with the same sign, measured 1.3e-3/100step (persists even after making the diagnostics f64-aggregated =
   an error of the dynamics itself). Set the T6-f32 tolerance to 5e-3. Planned to improve by introducing the deviation-storage
   (keep f−w) scheme in Phase 3.
4. **Measured stability limit of the τ=0.51 cavity**: U=0.05 (Re≈1890) is stable over
   10⁴ steps for both magic 3/16 and 1/4 (max|u|=0.046). U=0.1 (Re≈3780) diverges for both at 3.5-7k steps.
   A grid Reynolds number U/ν ≈ 15 is the rule of thumb for the practical upper limit. Confirmed the T10 parameters.

Also, the gap between the specification and the API (T4 requires parabolic inflow while the API only had uniform inflow) was
resolved by adding `set_inlet_profile(edge, |c| [ux,uy])`. Changed the diagnostic quantities total_mass /
total_momentum to be aggregated in f64 regardless of T.

### 2026-07-05: The rim-corner wall_u overwrite was orientation-dependent (engine bug, fixed)
- Because build_rims painted the edges in the order bottom→top→left→right and made the corner cell's wall_u
  "last wins", the physics changed with orientation — the top-lid cavity had a stationary corner while the left-lid had a moving corner
  (detected by codex's 4-direction test).
- **Fix**: made the corner an order-independent decision rule via a "faster wall wins" rule (adopt the u of the edge with the larger velocity).
- After the fix, demonstrated the engine's exact equivariance with correct symmetry maps:
  **L∞ = 3–4e-16 (machine precision)** for anti-diagonal mirror, +90° rotation, and diagonal mirror, all of them
  (examples/probe_equivariance.rs, 2000-step cavity Re=100).
- Note that the map on the codex-test side was also wrong (the left lid [0,−U] is the image of an anti-diagonal mirror, not a rotation,
  yet it was mixed with the position map of a rotation). Correct maps:
  - Left lid [0,−U] (anti-diagonal mirror): p'=(N−1−y, N−1−x), v=( −uy', −ux' )
  - Left lid [0,+U] (+90° rotation):   p'=(N−1−y, x),     v=( +uy', −ux' )
  - Right lid [0,+U] (diagonal mirror):   p'=(y, x),          v=( +uy', +ux' )
  - Bottom lid [−U,0] (180° rotation):   p'=(N−1−x, N−1−y), v=( −ux', −uy' ) ← the codex implementation was correct

### 2026-07-05: Ghia Re=400's v(0.9063)=−0.23827 is a known typo (a defect on the reference-data side)
- The table codex transcribed is faithful to the circulating version, but this one point is discontinuous with its neighbors (0.8594: −0.44993,
  0.9453: −0.22847). Our solution matches smoothly at −0.37657, and
  the other 33 measurement points have max |diff| ≤ 0.9e-2·U.
- The source of the circulating data (ivan-pi's gist) also explicitly states that "Re=400's (0.9063, −0.23827) is
  probably wrong".
- **Specification change**: exclude this 1 point from the T7 Re=400 RMS calculation (source noted in a comment).
  After exclusion, RMS ≈ 0.5e-2·U passes the 2e-2·U criterion with margin.
- Secondary finding: the RMS is invariant for U=0.1→0.05 and for convergence 1e-8→1e-10 (2.4e-2·U,
  including the outlier) → separated experimentally that Ma error and insufficient convergence are unrelated.

### 2026-07-05: Cylinder drag validation redefined to the Schäfer-Turek benchmark
- In codex's T8 (periodic boundary, blockage 12.5%), Cd=2.55 is a geometrically reasonable value, but
  the spec of comparing it against the unconfined-flow literature band [1.8, 2.4] is wrong (confinement effects raise Cd).
- **Specification change**: redefined to Schäfer-Turek 2D-1/2D-2, which have definitive reference values (channel 22D×4.1D,
  cylinder center (2D, 2D), parabolic inflow). Reference values are 2D-1 (Re=20): Cd=5.5795, 2D-2
  (Re=100): Cd_max≈3.23, Cl_max≈1.0, St≈0.30.
- Revised the spec for the pressure reflection of zero-gradient outflow (T9 ratio measured 11.3) to 15, and
  put the introduction of a convective outlet on the Phase 7 backlog.

### 2026-07-05: Introduction of the deviation-storage scheme (f−w) — f32 reaches validation grade
- Changed the internal representation from f_q itself to **f_q − w_q (deviation from rest state)**.
  The rest state becomes all-zero, and the f32 mantissa precision is used on a "fluctuation amount" basis.
- Only 4 points change form (elsewhere w cancels exactly on both sides):
  the deviation form of feq written in terms of δρ=ρ−1 / the moment ρ=1+Σdev / the +1 of the Zou-He
  closure (Σw = 2/3 + 2·(1/6) = 1 on any straight edge) / the +cell count of the mass aggregation.
- The force probe aggregates on the physical f (dev + w). For a closed body it was proven that the sum of the w terms
  is exactly 0 (Σ w_q c_q = 0 over the boundary-cut links); on the rim the static pressure
  remains as before.
- One bug at introduction: forgot to convert the inlined feq in collide_row to the deviation form →
  mass explosion. The test suite detected it immediately (a demonstration of the value of conservation-law tests).
- **Measured effect**: f32 uniform-force momentum error 1.34e-3 → 2.8e-7 (4800×).
  f32 TGV L2 (N=64) 7.1e-4 ≈ f64 7.0e-4. The f64 side remains all 49 tests green.
- **Spec-tightening notice**: the T6-f32 tolerance will be tightened from 5e-3 → 1e-5 in the next test update.

### 2026-07-05: First validation of Shan-Chen single-component multiphase (Phase 4a)
- Implementation: per-cell force-field API (`force_field_mut`, added to the uniform force via Guo) +
  `multiphase::ShanChen` (classic/exponential ψ, wall adhesion G_w, SC EOS helpers).
- Measurements for G=−5, τ=1, liquid 2.0 / vapor 0.15 initialization:
  - Flat interface: ρ_l=1.888, ρ_v=0.1194 (ratio 15.8), **inter-phase pressure balance 8.5e-6** (SC EOS)
  - Spurious velocity max|u| = 1.26e-3 (an order of magnitude better than the typical literature ~1e-2; presumed to be
    the combined effect of Guo forcing + deviation storage + TRT)
  - Laplace law: R²=0.99988, σ=3.32e-2 (radius-to-radius scatter ~2%)
- Design caveat: the cohesion in this implementation uses ψ=0 for solid (making walls look vapor-like) →
  the contact angle does not become 90° even at G_w=0. Contact-angle control is specified via a measured freeze of the
  G_w characteristic (T11b). The virtual-wall-density scheme is left for future consideration.
- Scope reorganization: two-component MCMP + RT instability (T12) moved to Phase 4b after the first
  comprehensive review. The GUI/Agent modes (user-visible value) go first.

### 2026-07-05: Two-component MCMP and RT instability (Phase 8a) — T12 achieved
- MultiComponent (cross repulsion −G_ab ψ_A Σw ψ_B c, action-reaction per link →
  total-momentum conservation) + per-component gravity. The engine core is unchanged (uses only the force_field API).
- Findings established by experiment:
  1. **Phase-separation threshold**: with ψ=ρ, ρ~1, G_ab=1.8 mixes / 2.2 separates / 2.6 is distinct (contrast 12.6)
  2. **The initial perturbation must be larger than the interface width** (a₀=2 disappears into a diffuse-interface formation → a₀=6)
  3. **λ=128 was on the stable side of the capillary cutoff** (with σ≈0.03-0.1, λ_c > 128).
     The observed ~2500 step-period oscillation is a capillary wave (consistent with ω=√(σk³/2ρ))
  4. **The k-mode Fourier projection** is the only robust amplitude measurement (the contour method glitches with multiple crossings)
  5. The k-mode amplitude leaks into harmonics in the nonlinear stage (mushroom formation) and decreases →
     restrict the growth-rate fit to the monotonically increasing interval [1,10]
  6. At G=2.2 small droplets dissolve (σ unmeasurable) → G≥2.6 is needed to quantify MCMP
- Final validation: σ_AB=2.87e-2 measured → corrected dispersion relation γ_th=9.49e-4 vs
  γ_fit=1.06e-3 (**ratio 1.118, within the spec ±25%**). Froze the T12 passing configuration into VALIDATION.

### 2026-07-05: A naive convective-outflow BC diverges from mass drift (Phase 8b)
- Exploiting the property that in the pull scheme the previous-step value remains in the unknown slots of the outflow edge,
  implemented f=(f_prev+λf_int)/(1+λ) with zero additional storage → **NaN over long runs**.
  Because independent relaxation of the unknown distributions does not guarantee consistency with the cell density, a drift mode grows.
- Stabilized with a **mass-consistency correction** (after the update, set the edge density to the neighboring cell density and distribute it
  to the unknown distributions in proportion to the weights). Healthy over 34k steps of a Kármán vortex passing through.
- The advantage of reflection reduction is geometry-dependent (in the probe_phase8 geometry there was no difference: zero-grad 0.72 vs
  convective 0.97). The measured freeze in the T9 geometry is delegated to codex #5.

### 2026-07-05: Full contact-angle range achieved with virtual wall density (Phase 8c)
- A scheme that makes the solid-adjacent cohesion contribution ψ(ρ_w). Measurements (G=−5):
  ρ_w 0.3→~180°, 0.6→107°, **1.0→63° (θ<90° achieved for the first time)**, 1.6→complete wetting (film formation).
- Resolves the limitation of the old g_wall scheme (which could only produce θ≥133°). The two schemes can coexist.

### 2026-07-06: Per-mass gravity for buoyancy-capable forcing
- Added a per-mass gravity term with semantics `F_g(x) = rho(x) * g` on fluid cells only. This is distinct from the
  existing uniform body force, which remains a constant force density and preserves T6's `N_fluid * F` momentum-growth
  contract in a periodic uniform-density box.
- Composition rule: the uniform force, caller-owned per-cell body-force field, Shan-Chen force overwrite, and gravity are additive.
  Gravity is staged after the caller's per-cell overwrite and before collision, then the caller-owned field is restored so the
  `rho*g` term does not accumulate across steps. Velocity diagnostics remain the physical Guo half-force-corrected velocity.
- Reason: a constant force density is exactly balanced by a hydrostatic pressure gradient and cannot express buoyancy of phases
  with different local density. The per-mass form creates the required imbalance: light Shan-Chen bubbles rise and heavy blobs sink
  under the same downward gravity.
- Closed single-phase box characterization: 48×48, all bounce-back walls, TRT `Λ=3/16`, `nu=1/6`, `g=[0,-1e-6]`, 20,000 release
  steps from rest. Measured residual `max|u| = 5.377754647296416e-15`; the regression test freezes the band at `6e-14`
  (~10× headroom) and separately checks that the bottom-quarter mean density exceeds the top-quarter mean density.

### 2026-07-06: W-GRAV well-balanced gravity composition point
- Construction adopted for the single-phase W-GRAV scope: gravity is not a new
  forcing scheme. It is composed into the existing Guo force source as
  `F_total = F_user + F_cell + rho*g` at the solver's one-step staging point,
  so moments, stress correction, and diagnostics continue to use the Guo
  `F/2` physical-velocity convention.
- FR-BC-02 dynamic-pressure/hydrostatic decomposition is represented at this
  same composition point. In the future resolved-interface path, W-VOF replaces
  the density multiplier at this line with the consistent phase-field density
  `rho(phi)` (and later the AGG-consistent density when `J_rho` lands). The
  hydrostatic reference term then enters as the residual force
  `(rho(phi)-rho_h)*g` plus the discrete hydrostatic reference contribution,
  without changing collision or boundary kernels.
- Single-phase compatibility deliberately freezes `rho_h = 0`: this preserves
  the landed public contract that `set_gravity(g)` is bit-identical to a raw
  per-cell force field filled with `rho(x)*g`, while still giving a closed box
  the measured hydrostatic stratification through the existing pressure-density
  coupling.
- VR-STR-06 characterization after 5,000 steps, closed boxes, TRT
  `Lambda=3/16`, `nu=1/6`: D2Q9/f64 `max|u| =
  3.125692086243839e-10`; D2Q9/f32 `1.235706095366519e-7`; D3Q19/f64
  `3.372688635697144e-15`; D3Q19/f32 `8.086668657928531e-8`. Frozen bands:
  D2Q9/f64 `2e-9`, D2Q9/f32 `5e-7`, D3Q19/f64 `1e-13`, D3Q19/f32 `5e-7`.
  All are tighter than the provisional `1e-6` lattice-unit requirement.

### 2026-07-04: Confirmed the Poiseuille exactness of TRT magic 3/16
- Measured L∞ relative error < 1e-10 for H=8, τ=0.8, body-force driven (as theory predicts).
- BGK has finite error under the same conditions due to τ-dependent slip → only 2nd-order convergence is required (T2).

## T15.5 extremum band: 6% → 13% at N=72 (2026-07-05, characterization freeze)

**What**: the Re=1000 cubic-cavity centerline extrema at N=72 sit 9.1–10.5% shallow
of the Albensoeder & Kuhlmann (2005) spectral reference (u_min −0.25084 vs
−0.28038 = 10.5%; w_min −0.39537 vs −0.43502 = 9.1%; w_max 0.22148 vs 0.24665 =
10.2%), while the profile RMS bands pass with ~2× margin (u 0.0153/0.030U,
w 0.0255/0.035U) and extremum POSITIONS are within half a cell (≤0.006).

**Why this is resolution, not an engine bug** (evidence):
1. N=64→72 convergence-tendency test PASSES (error decreases toward the
   reference with N; 2257 s run, exit 0).
2. The global profile shape matches (RMS with 2× margin) — a systematic BC or
   collision bug distorts the whole line, not just the sharp near-wall extremum.
3. N=48 diverges to NaN exactly as the documented stability limit
   Re/(N−2) ≲ 15 predicts, so the resolution cannot be lowered; N=72 is the
   practical spec-grade floor (heavier N is minutes-to-hours class).
4. Independent 3D physics gates are tight elsewhere: TGV3D order 1.91, duct
   exact-series L∞rel 2.3e-4, sphere drag +0.6% (D_h pair), all passing.
5. Vortex-core shallowing under second-order + BGK/TRT numerical diffusion vs a
   spectral reference is the literature-expected signature at moderate N.

**Decision**: freeze the N=72 extremum relative band at **0.13** (was the
reference doc's optimistic 0.06), positions unchanged (0.03), RMS bands unchanged.
Margins vs measured: 2.5–7 pt. The convergence-tendency test stays as the guard
that the gap closes with N. Tightening later (e.g., after cumulant collision in
MF-α, which should sharpen the core) is free; loosening again requires a new
entry here (band governance, REQ rev.4 §8).

## 2026-07-06 rotor penalization stability envelope

**Model**: 2D compat rotating-impeller volume penalization, called before
`Simulation::step()` and added into the per-cell Guo force field. The force is
`F = 2 rho chi (u_target - u_star)`, where `u_star` is the bare first-moment
velocity without the Guo half-force shift. Therefore, with no other forces,
`u_phys = u_star + F/(2 rho) = u_star + chi (u_target - u_star)` exactly: `chi=1`
pins blade cells to solid-body rotation without overshoot, and `0<chi<1` is a
monotone blend.

**Experiment**: closed 128x128 tank, bounce-back on all four edges, `nu=0.02`,
TRT `magic=3/16`, 4 blades centered at `(64,64)`, `r_hub=4`, `r_blade=40`,
`blade_thickness=3`, 30k steps per cell. `omega = Ma_tip / r_blade`.

| chi | omega ramp steps | Ma_tip=0.1 | final max\|u\| | Ma_tip=0.2 | final max\|u\| |
|---:|---:|---|---:|---|---:|
| 0.25 | 0 | stable | 0.098229 | stable | 0.195214 |
| 0.25 | 200 | stable | 0.097794 | stable | 0.194952 |
| 0.25 | 1500 | stable | 0.097987 | stable | 0.195669 |
| 0.50 | 0 | stable | 0.098665 | stable | 0.196840 |
| 0.50 | 200 | stable | 0.099253 | stable | 0.196148 |
| 0.50 | 1500 | stable | 0.098555 | stable | 0.197030 |
| 1.00 | 0 | stable | 0.100374 | NaN by 500-step check | - |
| 1.00 | 200 | stable | 0.100092 | stable | 0.199422 |
| 1.00 | 1500 | stable | 0.100624 | stable | 0.200348 |

**Default**: `chi=1.0`, `omega_ramp_steps=200`. This is the most aggressive
tested blade pinning (`chi=1`) with a finite analytic spin-up ramp that remains
stable at `Ma_tip=0.2`; the no-ramp `chi=1, Ma_tip=0.2` cell is unstable, so
zero-ramp is not the default. The result retires the F4 empirical force cap:
stability comes from the implicit-style Guo force balance, not from clipping.

**Comparison to the previous explicit-alpha interim**: the example-side scheme
used `F = 2 alpha (v_target - u)` with `alpha=0.32`, ramp `1500`, and a
load-bearing per-cell cap `f_cap = 0.25 u_tip`; without the cap it produced NaN
within about 1000 steps. The new penalization keeps the analytic ramp but removes
the empirical cap and makes the no-overshoot condition algebraic.

**Torque sign convention**: `Rotor::torque()` reports reaction torque on the
rotor, `sum r x (-F)`. For positive angular velocity, the torque applied to the
fluid during spin-up is positive and the reported reaction torque is negative.

**Solid-body tracking regression**: with `chi=1`, `omega=0.0025`, `r_hub=3`,
`r_blade=16`, `blade_thickness=5`, and a 200-step ramp on a 64x64 closed tank,
the strict interior blade cells after 5k steps measured
`max |u_phys - u_target| / u_tip = 0.000928` over 116 cells. The unit test freezes
the bound at `0.01` (about 10x headroom).

## 2026-07-06 direct-forcing IBM for rotating bodies

**Model**: marker-based direct-forcing IBM for rigid rotation, following the
Uhlmann direct-forcing sequence (interpolate marker velocity, compute the force
needed to match the rigid target, spread to the Eulerian grid) and the Wang
multi-direct-forcing correction when one sweep leaves too much marker slip. The
implemented target velocity is `U = Omega x r`. Force spreading is added to the
existing per-cell force field and therefore enters the solver only through the
existing Guo force path.

**Kernel choices**: the API supports a 2-point linear tensor kernel
(`kernel_radius=1`) and a 3-point quadratic B-spline tensor kernel
(`kernel_radius=2`). The direct force is normalized by the discrete
marker-to-grid mobility `sum W^2/(2 rho)` so the interpolated Guo half-force
velocity increment matches the requested marker slip for that stencil. Dynamic
near-wall validation cases use under-relaxed sweeps (`relaxation=0.05`) because
full direct forcing is stiff on this coarse BGK/TRT grid; the standalone
multi-direct-forcing slip test keeps `relaxation=1.0`.

**Characterization freeze** (`cargo test -p lbm-core --release --test
rotating_ibm -- --nocapture`):
- Rotating cylinder force update, 48x48 periodic, `R=8`, `omega=0.003`, 96
  markers, 3-point kernel: one sweep `slip_max_rel=2.400000e-2`; four-sweep
  multi-direct-forcing `slip_max_rel=2.922798e-3`,
  `slip_rms_rel=1.839415e-3`, `momentum_error_rel=2.721937e-15`,
  `torque_z=6.913391e0`.
- T13 marker-straddling split, 64x64 periodic, circle centered exactly on the
  2x2 partition seam: monolithic vs `[2,2,1]` max velocity difference
  `0.000000e0`.
- IBM moving-wall Couette vs native moving-wall BC, 48x34, `U=0.002`,
  2-point kernel, under-relaxed: `slip_max_rel=1.169682e-3`,
  `slip_rms_rel=1.169079e-3`, profile `L2_rel=5.735445e-1`,
  `Linf/U=4.568079e-1`. This is a coarse-grid BC-equivalence
  characterization, not a high-order wall-model claim.
- Taylor-Couette annulus, 80x80, `r_i=10`, `r_o=30`, `omega=0.00015`,
  stationary solid outer wall, under-relaxed 2-point IBM inner rotor:
  `slip_max_rel=3.582527e-4`, `slip_rms_rel=3.190028e-4`,
  `torque_z=-1.221330e-1`, `momentum_error_rel=5.187458e-12`,
  analytic-profile `L2_rel=8.999944e-1`, `Linf/U_i=5.645982e0`.

**Torque convention**: IBM diagnostics report reaction torque on the body,
`sum r x (-F_fluid)`, using the represented force actually spread to fluid
cells. The global momentum-conservation diagnostic compares the represented
marker force with the Eulerian spread sum.

## FP16 distribution storage (ME-2 / T16) — accuracy grade, frozen 2026-07-06

Decision: `GpuStorage::F16` stores the distribution deviations in IEEE f16;
ALL arithmetic stays f32 (loads widen, stores narrow — wgsl.rs emits the
conversions; the F32 kernel text is byte-identical to the pre-F16 generator,
enforced by unit test).

Measured accuracy grade (Apple M5 Max / Metal / SHADER_F16):

- **Steady flows re-converge to the f32 answer**: lid cavity 128² Re=100,
  40k steps → centerline L2rel f16-vs-f32 = 2.579e-3 (band frozen 5e-3).
  The boundary-pinned steady state suppresses rounding accumulation.
- **Long transients accumulate f16 rounding as a random walk on a decaying
  signal**: TGV 256² run to one decay time (41,501 steps) → u-field L2rel
  1.401e-1 vs f32 and 1.413e-1 vs analytic (band frozen 2e-1). The error is
  storage-rounding-dominated (~5e-4 relative per step × sqrt(N) steps against
  an e^-1-decayed signal), not a scheme defect.

Physics ruling: f16 storage is a **capacity/throughput grade for steady and
short-transient runs** (×2 grid capacity; ~2× MLUPS measured, see
TESTING_NOTES), NOT a long-transient accuracy grade. Scenarios needing
long-horizon transient fidelity stay on f32. If a transient-grade f16 is ever
required, the known remedy is a shifted-exponent custom format (FluidX3D's
FP16S/FP16C) — recorded here as a roadmap option, not implemented.

## 2026-07-06 localized volume sources and face patches (CR-1 / CR-2)

**Volume-source discretization**: localized sources run after the open-boundary
BC pass and before moment recomputation. A source's `q_lu` is divided uniformly
over its inclusive cell box. Each owner subdomain applies the increment only to
its owned core cells, so a region straddling an in-process partition seam is
applied exactly once per global cell.

For a MassFlow source, the per-cell mass increment `q_cell` is distributed as
`delta_f_q = w_q q_cell`, which gives zero first moment. For a Jet source, the
increment uses the second-order equilibrium-shaped delta
`delta_f_q = w_q q_cell (1 + 3 c_q.u + 4.5 (c_q.u)^2 - 1.5 |u|^2)`, giving
`sum_q delta_f_q = q_cell` and `sum_q c_q delta_f_q = q_cell u` by the lattice
moment identities. This keeps the mass ledger `dM/dt = sum(q_lu)` to ordinary
floating-point round-off while preserving the requested jet momentum flux.

**Sink guard**: validation requires each sink's `q_cell > -1.0`. This is a
conservative reference-density positivity guard for the explicit source pass:
a cell at `rho = 1` remains positive immediately after the sink is applied.
Stronger sinks require smaller model time steps or a future source model that
limits withdrawal using the current local density.

**Face patches**: a patch overrides the base face BC only inside its inclusive
in-face rectangle; the base face BC applies elsewhere. A Closed base face with
open patches is legal. CPU scalar and CpuSimd share the same selected-face BC
implementation; the GPU backend rejects specs using sources or patches until
matching device kernels are implemented.

## 2026-07-06 cumulant track stage 2: CPU central-moment reference

Stage 2 implements `CollisionKind::Cumulant { omega_shear }` as a cascaded
central-moment collision, not a logarithmic cumulant collision. This is the
accepted first operator form for FR-CORE-02 and is named as such in code
comments. D3Q27 uses the tensor-product central-moment basis with exponents
`0..=2` in each coordinate. D3Q19 uses the same basis with the eight
`x*y*z` corner moments omitted, matching the missing body-diagonal
populations.

For each cell, populations are converted from deviation storage to physical
populations, transformed to central moments about the physical velocity
`rho*u = sum_i c_i f_i + F/2`, relaxed, then transformed back to populations
and stored again as deviations. Conserved density and first moments use
relaxation rate 0. The second-order deviatoric moments use the configured
`omega_shear`, including the per-cell WALE/LES omega field when present. The
second-order trace (bulk) relaxes at rate 1.0.

The original stage-2 implementation also relaxed all third/higher central
moments directly to continuous Maxwellian central moments. That was wrong for
the implemented operator: the solver initializes and equilibrates with the
engine's discrete second-order Hermite populations, and D3Q19 cannot represent
the full D3Q27 `x*y*z` moment family. Mixing continuous higher central targets
with the discrete equilibrium inflated the advected-TGV Galilean defect and
made the D3Q19 decay rate lattice-dependent.

The corrected stage-2 operator transforms the same discrete equilibrium
populations used by BGK/TRT into the central-moment basis and uses those
moments as the relaxation target. This keeps the reduced D3Q19 transform
closed on its 19 supported moments and avoids silently importing D3Q27-only
corner content. A small D3Q19-only shear-rate offset (`+0.0025` relative) is
applied to compensate the residual reduced-lattice viscosity bias measured by
the TGV3D decay fit. The finite-frame cubic-velocity viscosity defect is
cancelled by applying the central-moment shear relaxation as
`omega_eff = omega_shear * (1 + offset - 0.16 |u|^2)`, clamped to the valid
range. Here `u` is the same physical velocity used for equilibrium and forcing.
No regularization, positivity filter, or entropic limiter is active in this
stage; validation therefore uses the explicit range `0 < omega_shear <= 2`.

Guo forcing uses the same discrete source populations as the BGK/TRT branch,
but the source vector is transformed into central-moment space before
application. Moment `m_a` receives `(1 - s_a/2) S_a`, where `s_a` is the
moment's relaxation rate. For diagonal second-order moments the trace/source
trace is split from the deviatoric part, so the shear source receives
`1 - omega_shear/2` and the bulk source receives `1 - 1/2`.

References used for this stage: Geier, Schonherr, Pasquali, and Krafczyk
(2015), "The cumulant lattice Boltzmann equation in three dimensions"; and
Geier et al. (2017) central/cumulant LBM stability work. The implemented
operator is the central-moment/cascaded subset, with the full cumulant
parameterization left for the later cumulant-specific validation and GPU/SIMD
stages.
