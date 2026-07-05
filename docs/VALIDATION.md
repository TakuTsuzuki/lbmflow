# VALIDATION.md — Boundary Condition / Physics Validation Test Spec Matrix

This file is the **order spec** for the validation test suite. The test author (codex)
shall write tests based only on this spec and the public API (do not copy the engine's
internal implementation). The acceptance criteria assume f64 and `--release` execution.

## Public API (used from tests)

```rust
use lbm_core::prelude::*;

let mut sim: Simulation<f64> = SimConfig {
    nx: 64, ny: 64,
    nu: 0.02,                        // kinematic viscosity (lattice units). tau = 3*nu + 0.5
    collision: Collision::Trt { magic: 0.1875 },   // or Collision::Bgk
    edges: Edges {
        left:   EdgeBC::Periodic,
        right:  EdgeBC::Periodic,
        bottom: EdgeBC::BounceBack,               // stationary wall (1-cell solid rim)
        top:    EdgeBC::MovingWall { u: [0.1, 0.0] },
    },
    force: [1e-6, 0.0],              // uniform body force (Guo)
    ..Default::default()
}.build().unwrap();                  // an invalid config yields Err(ConfigError)

sim.set_solid(x, y);                       // internal obstacle (any cell outside the rim)
sim.set_solid_region(|x, y| bool);         // bulk specification via predicate
sim.set_inlet_profile(Edge::Left, |c| [ux, uy]); // per-node profile for a VelocityInlet edge
                                           // c is the along-edge coordinate (left/right edge=y, top/bottom=x)
sim.init_with(|x, y| (rho, ux, uy));       // initialize with f = feq(rho,u) (rho=1, u=0 if not called)
sim.step();                                // 1 time step
sim.run(n);                                // n steps

sim.nx(); sim.ny(); sim.time();            // shape / elapsed steps
sim.rho(x, y); sim.ux(x, y); sim.uy(x, y); // macroscopic quantities (physical velocity, incl. force half-correction)
sim.rho_field(); sim.ux_field(); sim.uy_field(); // &[T] (cell = y*nx + x)
sim.is_solid(x, y);
sim.total_mass(); sim.total_momentum();    // Σρ, [Σρux, Σρuy] (fluid cells only)
sim.set_force_probe(|x, y| bool);          // target solid set for momentum-exchange force measurement
sim.probed_force();                        // [Fx, Fy] of the most recent step
sim.fluid_cell_count();                    // number of non-solid cells (for momentum tests)
```

- Edge kinds: `Periodic` / `BounceBack` / `MovingWall{u}` / `VelocityInlet{u}` (Zou-He) /
  `PressureOutlet{rho}` (Zou-He) / `Outflow` (zero gradient)
- A wall edge is realized as a solid rim (1 cell). **The wall surface lies half-way between
  the rim cell center and the adjacent fluid cell center** (half-way). Example: for top/bottom
  walls with grid height `Ny`, the rim is at y=0 and y=Ny-1, the fluid rows are y=1..=Ny-2, and
  the wall surfaces are at y=0.5 and y=Ny-1.5. Therefore
  **channel width H = Ny-2 (= number of fluid cell rows)**, and the distance of a fluid cell
  center from the wall is y_w = j - 0.5 (j = 1..H). The Poiseuille maximum velocity is g·H²/(8ν).
- Build-time errors (`ConfigError`): tau ≤ 0.5 (nu ≤ 0), unpaired Periodic,
  an orthogonal edge of a Zou-He/Outflow edge that is neither wall nor Periodic, nx/ny < 3, etc.

## Notation

- Error norm: `L2rel(u, u_ref) = sqrt(Σ|u-u_ref|²) / sqrt(Σ|u_ref|²)` (fluid cells only)
- Convergence order: `order = log2(err(N) / err(2N))`
- Steady-state criterion: `max|u^{t+Δ} - u^t| / max|u| < ε` (Δ=500 step). **ε = 1e-11 is recommended**:
  BGK oscillates permanently at a round-off error plateau of ~1e-12, so 1e-13 is unreachable
  (see the experiment record in docs/PHYSICS.md).

---

## Test Matrix

### T1. Taylor–Green vortex (periodic boundaries, viscous decay, convergence order)
- Setup: Periodic on all edges, N×N (N=32, 64), ν=0.02, **diffusive scaling u0 = 1.28/N**,
  k=2π/N. Analytic solution `ux = -u0 cos(kx) sin(ky) e^{-2νk²t}`, `uy = +u0 sin(kx) cos(ky) e^{-2νk²t}`.
  Initialization must pass, via `init_with`, the **pressure-consistent density ρ = 1 − (3u0²/4)(cos 2kx + cos 2ky)**
  (a uniform ρ=1 causes O(u0) contamination from residual acoustic waves; see docs/PHYSICS.md).
- Acceptance criteria:
  - At t = 1/(2νk²), velocity-field L2rel ≤ 1.5e-3 (N=64, TRT; measured 7.0e-4)
  - Convergence order order ≥ 1.7 (N=32→64; measured 1.91)
  - The effective viscosity ν_eff from the decay-rate fit is within ±2% of the nominal ν (N=64)
- Angle: equivalent for BGK and TRT; for an initial field rotated by 90° the result is rotationally symmetric (L∞ ≤ 1e-12).

### T2. Body-force-driven Poiseuille flow (exactness of half-way BB)
- Setup: BounceBack top/bottom, Periodic left/right, force F=[g,0] (e.g. g=1e-6), arbitrary ny (H=ny-2).
  Analytic solution: `ux(y) = g/(2ν) * y_w (H - y_w)`, y_w = (distance of cell center from wall) = j-0.5.
- Acceptance criteria:
  - TRT (Λ=3/16): at steady state (ε=1e-11), L∞rel ≤ 1e-10 (**exact**; holds even at H=8)
  - BGK: convergence order ≥ 1.7 for H=8→16 (BGK has a τ-dependent slip error, so exactness is not required)
  - Top/bottom symmetry of the profile: |ux(j) - ux(H+1-j)| ≤ 1e-13
- Angle: the same profile results even when the same setup is rotated by 90° (left/right walls, F=[0,g]).

### T3. Couette flow (moving wall)
- Setup: top MovingWall{u:[U,0]} (U=0.1), bottom BounceBack, left/right Periodic.
  Analytic solution: `ux(y_w) = U * y_w / H` (wall position referenced to half-way).
- Acceptance criteria: at steady state L∞rel ≤ 1e-10 (all of BGK/TRT and τ∈{0.6, 1.0, 1.4}).
- Angle: equivalent when the bottom wall is moved or a vertical Couette is set up with left/right walls. Mass conservation of the moving wall (total-mass drift ≤ 1e-12 relative / 10⁴ step).

### T4. Zou-He velocity inlet + pressure outlet channel
- Setup: left VelocityInlet + a parabolic profile via `set_inlet_profile`
  (u_max=0.05, the analytic form from T2, rim coordinate is [0,0]), right PressureOutlet{rho:1},
  BounceBack top/bottom. 96×34 (H=32), TRT.
- Acceptance criteria (after steady state):
  - The mass flux Q(x)=Σ_y ρ·ux is constant in the **bulk region** (a cross-section at least 24 columns
    away from the outflow boundary): max|Q−Q̄|/Q̄ ≤ 1e-4 (measured 2.4e-5)
  - Central-cross-section profile L2rel ≤ 2e-3 vs the parabola
  - Total-mass drift ≤ 1e-11 over 10⁴ step after steady state (measured 2e-13)
- **Known artifact (by spec)**: an O(Ma²) staggered oscillation appears in the ~4 columns immediately
  before the pressure-outlet boundary (about ±2%, independent of the collision operator, decay length ~4 cells).
  This is an intrinsic characteristic of the Zou-He pressure boundary and is not a failure condition (see PHYSICS.md).
- Angle: equivalent results in all 4 directions (left→right, right→left, bottom→top, top→bottom).

### T5. Pressure-difference-driven channel (Zou-He pressure-pressure)
- Setup: left PressureOutlet{rho_in}, right PressureOutlet{rho_out} (small Δρ, e.g. 2e-3),
  BounceBack top/bottom. Analytic: Poiseuille with dp/dx = cs²Δρ/L.
  **L = nx−1** (the pressure-specified nodes lie on the boundary columns, and their spacing is the effective channel length).
- Acceptance criteria: steady-state flow rate within ±2% of the analytic value (TRT, H=32; measured 0.26%).
  Linearity of the pressure field p(x)=cs²ρ(x) (bulk, excluding 8 columns at each end) R² ≥ 0.999.
- Angle:
  - **Exact**: with a Δρ sign reversal + x mirror, the field matches exactly under mirroring (L∞ ≤ 1e-12;
    by the x-inversion symmetry of the discrete system)
  - **Approximate**: a plain sign reversal without mirroring breaks the symmetry at O(Ma²) via the inertial term
    and compressibility, so relative L∞ ≤ 5e-3 (measured 1.7e-3). Requiring an exact 1e-12 here is physically wrong.

### T6. Conservation laws / consistency
- Periodic box + arbitrary initial field: total mass constant to within relative 1e-11 over 10⁴ step
  (because round-off accumulation of ~1e-13/10³step has been measured; physically it is exactly conserved).
- BB box (walls on all edges): same as above.
- Periodic box with uniform force F: total momentum grows by `N_fluid * F` per step (relative 1e-10).
- The 0th/1st/2nd moment identities of feq (unit test): Σfeq=ρ, Σfeq c=ρu, Σfeq cc = ρ(cs²I+uu)
  (1e-14 at several points with |u|≤0.1).
- Angle (f32): mass drift ≤ 1e-5 (10³step), relative error of the force-driven momentum growth
  ≤ 1e-5 (10²step). **After the deviation-storage scheme (introduced 2026-07-05), the measured value is 2.8e-7**
  (before introduction it was 1.3e-3 due to a coherent round-off bias on the uniform field; see PHYSICS.md).

### T7. Lid-driven cavity (comparison with Ghia et al. 1982)
- Setup: walls on all edges, top edge MovingWall{[U,0]}, Re = U*L/ν ∈ {100, 400, 1000}, N=129
  (L=N-2), U=0.1, TRT. Until steady state (ε=1e-8 or a 300k step cap).
- Acceptance criteria: compare the geometric centerlines u(y) / v(x) with the 17 points of the Ghia table,
  RMS error ≤ 0.02·U (Re=100/400), ≤ 0.03·U (Re=1000). The primary vortex center position is within ±0.02L of the reference.
  **Known typo**: Re=400 v(x=0.9063)=−0.23827 is a known error in the circulated data
  (discontinuous with neighboring points; see PHYSICS.md 2026-07-05). Exclude this one point from the RMS.
- Angle: same solution when the lid direction is rotated in 4 directions. Use the **correct symmetry map**
  (described in PHYSICS.md; e.g. the left lid [0,−U] is the anti-diagonal mirror p'=(N−1−y,N−1−x), v=(−uy',−ux')).
  Acceptance criterion L∞ ≤ 1e-10 (measured at machine precision ~4e-16). Allowed at Re=100, 2000 step.

### T8. Flow around a cylinder — Schäfer–Turek benchmark (force measurement, vortex shedding)
Adopt a standard benchmark with established reference values (Schäfer & Turek 1996, "Benchmark computations
of laminar flow around a cylinder"). Geometry (strict ratios):
channel 22D × 4.1D, cylinder center at 2D from inflow and 2D from the bottom wall (**slightly asymmetric**,
which is the trigger for vortex shedding; center y/H = 0.4878). BounceBack top/bottom,
left VelocityInlet + `set_inlet_profile` parabola u(y) = 4 u_max y_w(H−y_w)/H²,
right PressureOutlet{1.0}. U_mean = (2/3) u_max, Re = U_mean·D/ν.
Cd = 2Fx/(ρ U_mean² D), Cl = 2Fy/(ρ U_mean² D).

- **2D-1 (Re=20, steady)** reference values: Cd = 5.5795, Cl = 0.0106, Δp* = Δp/(ρU_mean²) = 2.9375
  - D=20 (grid 440×82, u_max=0.075, ν=0.05 → Re = 0.05·20/0.05 = 20)
    default suite: Cd ∈ [5.2, 6.0], Cl ∈ [−0.05, 0.08] (staircase coarse-grid band)
  - D=40 (grid 880×164, u_max=0.075, ν=0.1) #[ignore]:
    Cd ∈ [5.35, 5.85], convergence trend (|Cd(40)−5.5795| < |Cd(20)−5.5795|)
- **2D-2 (Re=100, unsteady)** reference values: Cd_max ≈ 3.22–3.24, Cl_max ≈ 0.99–1.01,
  St ≈ 0.295–0.305. D=40, u_max=0.15 (U_mean=0.1, ν=0.04) #[ignore]:
  - St ∈ [0.28, 0.32] (measured from the zero crossings of Cl)
  - Cd_max ∈ [3.0, 3.5], Cl_max ∈ [0.8, 1.2]
  - Periodicity of vortex shedding: the variation in the length of consecutive Cl periods ≤ 2%
- Note: because of the staircase approximation, the bands are wider than the reference values. They will be
  tightened when curved boundaries (a Phase 7 candidate) are introduced. The old spec (comparison against periodic
  boundaries / unconstrained bands) is retired due to geometric inconsistency (PHYSICS.md 2026-07-05).

### T9. Soundness of Outflow (zero gradient)
- Setup: in a channel equivalent to T8-2D-2, replace the right edge with Outflow. It does not diverge even when
  a vortex passes through the outflow face.
- Acceptance criteria: no NaN/Inf over 10⁵ step, reverse-flow mass flux at most 5% of the total inflow, and the
  rms of pressure oscillations near the outflow face (x>0.9L) is **within 15×** that of the central region (measured 11.3;
  zero-gradient outflow partially reflects pressure waves, an intrinsic characteristic. An improvement, a convective
  outlet, is considered in the Phase 7 backlog; PHYSICS.md 2026-07-05).

### T10. Robustness / error paths
- τ ≤ 0.5, unpaired Periodic, Zou-He orthogonal-edge violation, nx<3 → `ConfigError`.
- **Stability-limit case (parameters fixed)**: a cavity with τ=0.51, N=128, U=0.05 (Re≈1890),
  TRT Λ=3/16 has no NaN/Inf over 10⁴ step (measured stable at max|u|=0.046).
  U=0.1 (Re≈3780) diverges in ~3.5-7k step with both Λ=3/16 and 1/4 (a known limit,
  the grid Reynolds number U/ν ≈ 30 is exceeded. Guideline: for τ→0.5, U/ν ≤ 15).
- Placing `set_solid` on an open-boundary edge panics (by spec).
- If the moving-wall / inflow velocity is |u| > MAX_SPEED(=0.3), `ConfigError::VelocityTooHigh`.
  A velocity overflow in `set_inlet_profile` panics.

### T11. Shan-Chen single-component multiphase (Phase 4a; measurement-calibrated 2026-07-05)
Common setup: `ShanChen::new(-5.0)` (classic ψ = 1−e^{−ρ}), τ=1 (nu=1/6),
initialize liquid ρ=2.0 / vapor ρ=0.15, each step `sc.update_force(&mut sim); sim.step()`.
Compare pressure **always with the SC EOS** (`sc.pressure(rho)` = cs²ρ + (G cs²/2)ψ²).

- **Flat interface** (64×128 periodic, 30k step):
  - Coexistence densities ρ_l = 1.888 ± 2%, ρ_v = 0.1194 ± 3% (measured regression values)
  - Inter-phase pressure equilibrium: |p_l − p_v|/p ≤ 1e-4 (measured 8.5e-6)
  - Spurious velocity max|u| ≤ 5e-3 (measured 1.26e-3)
  - Total-mass drift ≤ 1e-10 relative (the SC force has no zeroth mass moment)
- **Laplace law** (128², R₀ ∈ {12,16,20,24}, 40k step):
  - Linearity of Δp vs 1/R_fit R² ≥ 0.999 (measured 0.99988)
  - Slope σ = 3.32e-2 ± 10% (measured regression value); each droplet's σ=Δp·R is within ±5% of the slope
  - Radius measured from the isocontour area at the median density (square root of area/π)
- **f32 angle**: the flat-interface case is stable even in f32 (no NaN, coexistence densities ±5%)

### T11b. Contact angle (freezing the G_w characteristic)
- With a wall-attached droplet (left/right Periodic, top/bottom BounceBack, the top wall far enough from the droplet),
  measure at G_w ∈ {−1.5, 0, +1.5}.
  Because this implementation uses ψ=0 for solid (excluded from cohesion) + a separate term −G_w ψ Σw s c,
  **G_w=0 does not give 90°** (it leans toward the non-wetting side). The test is:
  - θ(G_w) is monotonic (the more negative G_w, the more wetting = smaller θ)
  - Regression-freeze the 3 measured angles (±8°): **G_w=−1.5: 133.2°, 0: 160.4°, +1.5: 163.7°**
    (measured 2026-07-05). The measurement method is a spherical-cap fit (θ = 2·atan(2h/w)).
- Known constraint: the current scheme has a narrow range on the wetting side (it is hard to produce θ < 90°).
  A switch to the virtual-wall-density scheme (ψ(ρ_w) for solid) is placed in the Phase 7 backlog.

### T12. Two-component MCMP: Rayleigh–Taylor growth rate (measurement-calibrated 2026-07-05)
Self-consistent two-stage validation (`multiphase::MultiComponent`, ψ=ρ, G_ab=2.6, trace 0.05):
1. **Measurement of σ_AB**: A droplet in B (128², r₀=24, ν=0.1, 20k step),
   with p = cs²(ρ_A + ρ_B + G_ab ρ_A ρ_B), Δp·R_fit = σ_AB (measured ≈ 2.87e-2).
2. **RT growth rate**: 256×256 (left/right Periodic, top/bottom BounceBack), both components bulk ρ=1,
   at the interface y₀=128 a perturbation a₀=6·cos(kx) (k=2π/256), gravity on the heavy component only
   g_a=[0,−1e-4] (equivalent to an effective Atwood 0.5), ν=0.1.
   - Amplitude measured by the **k-mode Fourier projection of the column mass** (the contour method glitches)
   - ln(amp) regression over the monotonically increasing interval (amp ∈ [1,10]) → γ_fit
   - Reference: **γ_th = sqrt(gk/2 − σ_AB k³/2 + ν²k⁴) − νk²** (incl. tension / viscosity corrections)
   - Pass: γ_fit/γ_th ∈ [0.75, 1.25] (measured 1.118), amp reaches 10 or more
     (confirming the instability actually exists), total-mass drift ≤ 1e-10
- Separation smoke: with a half-and-half initialization (96², 5k step), G_ab=2.2 separates (contrast ≥ 3),
  G_ab=1.8 mixes (≤ 1.5) — confirming the existence of a phase-separation threshold.

### T11c. Full-range contact angle via virtual wall density (measured 2026-07-05)
`ShanChen::with_wall_rho(ρ_w)` (G=−5, liquid 2.0/vapor 0.15, 160×100, 30k step):
- Monotonicity: ρ_w ↑ → θ ↓. Regression-freeze: **ρ_w=0.3: ~180° (non-wetting), 0.6: 107°,
  1.0: 63° (θ<90° achieved)** each ±8°
- ρ_w=1.6 is complete wetting (film over the entire wall, contact width = full width) → as a qualitative case,
  assert that it "forms a film"

### T9b. Convective outflow (ConvectiveOutflow)
- Implementation with mass-consistency correction (pin the edge density to the adjacent cell). u_conv ∈ (0,1] is
  validated at build time (InvalidParameter).
- Compare against Outflow with the same geometry and same metrics as T9 and freeze at the measured values
  (in the PM's probe_phase8 geometry it is 0.97 vs 0.72, so the advantage is geometry-dependent.
  At minimum, require stability, non-divergence, and reverse flow ≤5%).

---

### T15. (M-C: when introducing 3D/D3Q19) 3D physics validation
Acceptance criteria for the core V2 D3Q19 (COMPETITIVE_SPEC R1).
Test body: `crates/lbm-core/tests/t15_3d.rs` (measured values are from the 2026-07-05 M-C implementation):
1. **Degeneracy match of the z-invariant 2D-TGV**: initialize a z-invariant 2D TGV on a 3D grid (N×N×4, z periodic)
   and **match field-by-field** with the identical D2Q9 scenario (f64 ≤1e-12; measured 8.9e-16/648step —
   the z-invariant projection is preserved at essentially bit precision).
   The most important smoke to first flush out bugs in D3Q19 streaming/weights/symmetry.
   Angle: confirm the same degeneracy also via 3D Zou-He faces (parabolic inflow + pressure-outlet channel)
   (measured 1.1e-15/500step). Also confirm as a unit the exactness of the prescribed moments at face nodes
   (|u−u_bc| ≤ 1e-14, measured 6.9e-18).
2. **Rectangular-duct Poiseuille (exact series solution)**: a rectangular duct of cross-section a×b, body-force-driven,
   4 walls half-way BB. u(y,z) = (g/2ν)·[analytic series] (Fourier series,
   u = (16a²g)/(νπ³) Σ_{n odd} (1/n³)[1 − cosh(nπz/2a)/cosh(nπb/2a)] sin(nπy/2a)).
   With TRT, L∞rel ≤ 1e-3 (the series is truncated at n≤99, incl. a convergence check; measured: for a 32² cross-section
   L∞rel 2.3e-4, truncation error ≤1e-4·umax confirmed).
   Flow rate Q matches the analytic value within ±0.5% (Q = (64a³g)/(νπ⁴) Σ (1/n⁴)[2b − (4a/nπ)tanh(nπb/2a)],
   measured 0.094%).
3. **Drag on a sphere**: a sphere in uniform flow (staircase), Re ∈ {20, 100}.
   Measure Cd by momentum-exchange and compare with the **Schiller-Naumann correlation**
   Cd = (24/Re)(1 + 0.15 Re^0.687) within ±10% (formula values Re=20: 2.6095, Re=100: 1.0917.
   The old wording "Re=20: ≈2.09" was a misstatement of the value at Re≈28 — see TESTING_NOTES.md 2026-07-05.
   The tolerance ±10% is unchanged).
   Blockage ≤ 3% (domain ≥ 8D), D ≥ 24 lattice (D=24, 192×128×128, periodic side faces —
   heavy, so #[ignore]. The default suite has a lightweight D=12 version, band ±15% with D_h normalization,
   measured +2.3%. The old wording "±25% / +14.2%" was the pre-D_h-normalization value — removed).
   Cd is the average over a 500-step window (to counter the weakly-damped acoustic ripple between inflow↔outflow, see TESTING_NOTES).
   **Normalization (triage confirmed 2026-07-05)**: use the hydrodynamic pair —
   Cd_h = F/(½ρU²π(r+½)²), reference SN(Re_h), Re_h = U(D+1)/ν. The half-way BB wall lies a half-link outside
   the solid cells (Ladd's staircase-sphere calibration), and the nominal D normalization was incompatible with
   the D=24 band because of the half-link bias (~+2/D).
   Measured (D_h basis): Re=20/D=24 **+7.1%**, Re=100/D=24 **+0.6%**,
   Re=20/D=12 lightweight **+2.3%** (band ±15%) — all pass.
4. **3D-TGV (true 3D, low-Re convergence order)**: run the classic TGV u=(sin x cos y cos z, −cos x sin y cos z, 0)
   only for a short time (t = 0.1/(νk²), where vortex stretching is weak), and the initial decay rate matches the
   diffusive limit 2ν(3k²) (±2%; measured 0.11% at N=64), and order ≥1.7 for N=32→64
   (diffusive scaling u0∝1/N; measured 1.91).
   **Calibration of the u0 coefficient (measured 2026-07-05)**: because the classic 3D TGV is not an exact NS solution,
   the nonlinear deviation from the diffusive-limit reference becomes a **resolution-independent relative offset**
   under diffusive scaling (measured L2rel ≈ 0.13·u0/(νk); with u0=1.28/N it is 0.165, which destroys the convergence order).
   Freeze the coefficient at u0 = 1.28e-4/N, which is sufficiently smaller than the spatial-error floor
   (e32=1.29e-3, e64=3.44e-4) (offset ~2e-5, more than 6 digits of margin in f64).
5. **3D cavity = T15.5**: the reference data is [T15_5_CAVITY3D_REFERENCE.md](T15_5_CAVITY3D_REFERENCE.md)
   (Albensoeder & Kuhlmann 2005, Re=1000, provenance-verified, 7-digit cross-checked, smoothness-audited.
   **The lesson from the Ghia typo incident has been applied**). Acceptance band: centerline RMS ≤ 0.030U and other bands in that document.
   By the stability constraint Re/(N−2) ≲ 15, N ≥ 72. The test body `t15_5_cavity3d.rs` was adversarially authored
   under codex order #7 (ordered 2026-07-05). Re=100/400 are on hold because the original source has not been obtained.

---

### T13. Partition invariance (core V2 equivalence; adversarially verified under codex order #6)
Bodies: `t13_split_invariance.rs` / `t13_adversarial.rs` / `examples/mpi_t13.rs`.
- **InProcess (in-thread) partitioning**: for 1×1 vs 2×2 / 4×1 / 1×4 / (3D) 2×2×2, the
  fields (rho/u/all f planes) are a **bit match** with `assert_eq!(d, 0.0)`.
  Adversarial angles (order #6): obstacles/probes/Zou-He faces on the seam, a cavity whose lid straddles the seam,
  Shan-Chen ψ exchange, a corner droplet, etc. — 8 kinds, all withstood.
- **T13-MPI** (mpirun -n {1,2,4,8}): the rank-0 gathered field has max|Δ| = 0.0 vs the single-rank reference.
  For diagnostics (mass/momentum/probed_force/NaN count), only the f64 recombination difference of the rank partial-sums → Allreduce
  is tolerated: atol+rtol each 1e-12 (field) / 1e-11 (diagnostics). Reproduce: `./scripts/test_mpi.sh`.
- **Known blind spot**: the probe double-counting of the two_pass boundary shell (improvement spec E8/C-2) does not
  appear in the fields, so it cannot be detected by the T13 field comparison. After the C-2 shell fix, a
  two_pass on/off probe-match test on a width-1 axis will be added to this section's acceptance.

### T14. Backend equivalence (CPU vs Wgpu, `--features gpu`)
Body: `t14_backend_equiv.rs` (6 configurations: TGV / cavity / cylinder+probe / Zou-He inflow-outflow /
force field / moving wall, f32).
- Relative field difference ≤ 1e-5. Only the pressure-boundary case ≤ 1e-4 (fixed via a **1-ulp control test** that
  "the relaxation is a round-off-order difference, not a physics difference"). Diagnostic values (mass/momentum) match.
- f64+GPU is rejected at compile time (no silent degrade).
- **Known gap (improvement spec D-3, resolved in R-Phase 2)**: currently only CPU relative equivalence, so a
  shared-spec-interpretation bug where CPU/GPU break in the same direction cannot be detected. Two absolute-physics
  tests directly on GPU (TGV convergence order ≥1.7 / cavity Ghia RMS ≤ 0.02U — calibrated and frozen against the f32
  measurements) will be added. Absence of an adapter is a skip; `LBM_REQUIRE_GPU=1` promotes it to a fail.

### T16. (M-E: when introducing FP16 storage) Precision-mode validation — **not yet implemented**
- Quantify the degradation of f16 storage (f32 arithmetic) on TGV / cavity and freeze the tolerance band.
  Deviation storage (f−w) is a prerequisite. Implementation is improvement spec C-12 → the M-E body.

### T17. (M-F: coupled multiphysics) VR-STR acceptance matrix — **spec wired, awaiting implementation**
Source: [REQ_STIRRED_REACTOR.md](REQ_STIRRED_REACTOR.md) §8 (rev.1b).
The tests are adversarially authored by codex/Opus from the REQ and separated from the implementation (the conventional protocol).
"Frozen after implementation" = implement → characterization measurement → PHYSICS.md record → list the frozen values in this table.

| ID | Target | Acceptance criteria (the part already fixed in REQ) | Band status |
|---|---|---|---|
| VR-STR-01 | Single-phase stirring (standard baffled tank, non-aerated) | Rushton Np = experimental correlation ratio; verify impeller discharge velocity against PIV/LDA reference survey lines with L2/L∞rel. Np = P/(ρ_l N³D⁵), P = Ω T_q (no 2π double-counting) | tolerance % frozen after implementation |
| VR-STR-02 | Gas-liquid (split into 02a/02b/02c) | **02a single bubble**: verify U_t against the Grace diagram (Eo-Mo-Re). **02b bubble swarm**: ε_g spatial distribution, swarm rise velocity, (when coalescence/breakup allowed) d_32, ν_t response under BIT. **02c aerated stirring**: experimental correlation ratio of ε_g, d_32, k_L a | relative-error band frozen after implementation |
| VR-STR-03 | Shear / stress field | Separate MMS single-phase, curved Couette, rotating cylinder, non-Newtonian Poiseuille, multiphase static droplet. **Spurious velocity Ca_spurious < 10⁻³ (fixed)**. Near-wall L∞ survey-line design required | convergence order / L2/L∞ bands frozen after implementation |
| VR-STR-04 | Scalar / reaction | Taylor-Aris dispersion, reaction-diffusion front at a known Da, k_L a (formula stated explicitly). State the target Pe/Da/Sc for each test. SGS scalar uses Sc_t (default 0.7) | tolerance frozen after implementation |
| VR-STR-05 | Coupled regression / conservation | probe_state_hash bit equivalence is **limited to single-backend regression**. Set individual drift thresholds for mass, momentum, total scalar, gas-phase volume, particle count, and **energy-like quantities (treated as monitored quantities: kinetic E, interfacial free E, particle kinetic E)**. GPU/MPI are tolerance-based | thresholds frozen after implementation |
| VR-STR-06 | well-balanced hydrostatics | maintain \|u\| < ε in a static stratification (high density ratio 10³). **06+**: with active ON and C≡C_0, the same static behavior (F_b^scalar degenerates to exactly zero); with ∇σ=0, matches the constant-σ reference form | ε frozen after the discretization is decided |
| VR-STR-07 | Initialization independence | quasi-steady statistics match within threshold when the run-up / statistics start are varied | statistics window / threshold frozen after implementation |
| VR-STR-RELAX | Relaxation-mode equivalence (new in rev.2) | compare each relaxation extension with the corresponding fidelity reference solution: MRF→IBM reference (Np, survey lines, torque) / point-bubble→resolved reference (ε_g, d_32, k_La, budget) / one-way→two-way (particle statistics, mass-loading cap) / AMR→uniform (conserved quantities, interface position, coarse-fine budget) / aggressive f32→fidelity profile (drift, Ca_spurious, Np, curvature) | frozen when the relaxation extension is implemented (the initial version only reserves the trait/schema/validation items) |

## Test implementation conventions (for codex)

- Location: `crates/lbm-core/tests/validation_*.rs` (one theme per file).
- Shared helpers in `crates/lbm-core/tests/common/mod.rs`.
- Attach `#[ignore]` to heavy computations (T7 Re=1000, T8, T9); the CI-equivalent is
  `cargo test --release`, the full one is `cargo test --release -- --include-ignored`.
- No randomness (deterministic). Attach a message containing the measured value to each `assert!`
  (e.g. `assert!(err < 5e-3, "L2rel = {err}")`).
- Adding external crates is allowed only for `approx` (anything else requires consultation).
- Embed the Ghia reference data as a constant table inside the test file (with a source comment).
