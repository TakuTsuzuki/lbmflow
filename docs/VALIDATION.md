# VALIDATION.md — Boundary Condition / Physics Validation Test Spec Matrix

This file is the **order spec** for the validation test suite. Tests are written from
this spec + the public API only (not from engine internals). Acceptance assumes f64
and `--release`. Historical calibration prose lives in PHYSICS.md; this file states
the resulting bands.

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
}.build().unwrap();                  // invalid config → Err(ConfigError)

sim.set_solid(x, y);                       // internal obstacle (any cell outside the rim)
sim.set_solid_region(|x, y| bool);         // bulk specification via predicate
sim.set_inlet_profile(Edge::Left, |c| [ux, uy]); // per-node profile for a VelocityInlet edge
                                           // c is the along-edge coordinate (left/right=y, top/bottom=x)
sim.init_with(|x, y| (rho, ux, uy));       // initialize with f = feq(rho,u) (rho=1, u=0 if not called)
sim.step();
sim.run(n);

sim.nx(); sim.ny(); sim.time();
sim.rho(x, y); sim.ux(x, y); sim.uy(x, y); // macroscopic (physical velocity, incl. force half-correction)
sim.rho_field(); sim.ux_field(); sim.uy_field(); // &[T] (cell = y*nx + x)
sim.is_solid(x, y);
sim.total_mass(); sim.total_momentum();
sim.set_force_probe(|x, y| bool);
sim.probed_force();                        // [Fx, Fy] of the most recent step
sim.fluid_cell_count();
```

- Edge kinds: `Periodic` / `BounceBack` / `MovingWall{u}` / `VelocityInlet{u}` (Zou-He) /
  `PressureOutlet{rho}` (Zou-He) / `Outflow` (zero gradient).
- Wall edges are a 1-cell solid rim. **The wall surface is half-way between the rim cell
  center and the adjacent fluid cell center**. Channel width **H = Ny − 2** (fluid rows).
  Distance from wall of fluid cell j is y_w = j − 0.5. Poiseuille peak = g·H²/(8ν).
- Build-time errors (`ConfigError`): tau ≤ 0.5 (nu ≤ 0), unpaired Periodic, an orthogonal
  edge of a Zou-He/Outflow edge that is neither wall nor Periodic, nx/ny < 3.

## Notation

- `L2rel(u, u_ref) = sqrt(Σ|u−u_ref|²) / sqrt(Σ|u_ref|²)` over fluid cells.
- Convergence order: `order = log2(err(N) / err(2N))`.
- Steady-state criterion: `max|u^{t+Δ} − u^t| / max|u| < ε` (Δ=500). **Use ε = 1e-11**
  (BGK round-off plateau is ~1e-12; see PHYSICS.md).

---

## Test Matrix

### T1. Taylor–Green vortex (periodic, viscous decay, convergence order)
- Setup: Periodic all edges, N×N (N=32, 64), ν=0.02, **diffusive u0 = 1.28/N**, k=2π/N.
  Analytic `ux = −u0 cos(kx) sin(ky) e^{−2νk²t}`, `uy = +u0 sin(kx) cos(ky) e^{−2νk²t}`.
  Initialize via `init_with` with **ρ = 1 − (3u0²/4)(cos 2kx + cos 2ky)** (uniform ρ=1 seeds
  O(u0) acoustic contamination).
- Accept: at t = 1/(2νk²), L2rel ≤ 1.5e-3 (N=64, TRT); order ≥ 1.7 (N=32→64);
  ν_eff from the decay-rate fit within ±2% of ν (N=64).
- Angle: BGK/TRT equivalent; 90°-rotated initial field is rotationally symmetric (L∞ ≤ 1e-12).

### T2. Body-force Poiseuille (half-way BB exactness)
- Setup: BB top/bottom, Periodic left/right, F=[g,0] (g=1e-6). Analytic
  `ux(y) = g/(2ν) y_w (H − y_w)`.
- Accept: TRT (Λ=3/16) at ε=1e-11, L∞rel ≤ 1e-10 (**exact**, even at H=8);
  BGK order ≥ 1.7 (H=8→16, τ-dependent slip error means not exact);
  top/bottom symmetry |ux(j) − ux(H+1−j)| ≤ 1e-13.
- Angle: same profile under 90° rotation (walls left/right, F=[0,g]).

### T3. Couette (moving wall)
- Setup: top MovingWall{[U,0]} U=0.1, bottom BB, left/right Periodic.
  Analytic `ux(y_w) = U y_w / H`.
- Accept: L∞rel ≤ 1e-10 at steady state, for BGK/TRT and τ ∈ {0.6, 1.0, 1.4}.
- Angle: bottom moved / vertical Couette equivalent. Moving-wall mass drift ≤ 1e-12 rel / 10⁴ step.

### T4. Zou-He velocity inlet + pressure outlet channel
- Setup: left VelocityInlet + parabolic `set_inlet_profile` (u_max=0.05, rim origin [0,0]),
  right PressureOutlet{1}, BB top/bottom. 96×34 (H=32), TRT.
- Accept (after steady state):
  - Bulk mass flux Q(x)=Σ_y ρ ux is constant on cross-sections ≥ 24 cols from the outlet:
    max|Q−Q̄|/Q̄ ≤ 1e-4.
  - Central profile L2rel ≤ 2e-3 vs the parabola.
  - Total-mass drift ≤ 1e-11 over 10⁴ step post-steady.
- **Known artifact (by spec)**: staggered O(Ma²) ripple in the ~4 cols immediately before
  the pressure outlet (±2%, decay length ~4 cells), intrinsic to Zou-He (see PHYSICS.md).
- Angle: equivalent in all 4 directions.

### T5. Pressure-difference channel (Zou-He pressure-pressure)
- Setup: left PressureOutlet{ρ_in}, right PressureOutlet{ρ_out} (Δρ=2e-3), BB top/bottom.
  Analytic Poiseuille with dp/dx = cs²Δρ/L, **L = nx−1** (specified nodes on the boundary columns).
- Accept: flow rate ±2% of analytic (TRT, H=32); bulk pressure linearity R² ≥ 0.999 (excluding
  8 cols at each end).
- Angle:
  - **Exact** mirror: Δρ sign reversal + x-mirror gives L∞ ≤ 1e-12 (discrete symmetry).
  - **Approximate** plain sign flip: relative L∞ ≤ 5e-3 (inertial + compressibility O(Ma²)
    break the mirror — asking for 1e-12 here is physically wrong).

### T6. Conservation laws / consistency
- Periodic box, arbitrary initial field: total mass constant to 1e-11 rel over 10⁴ step.
- BB box: same.
- Periodic box + uniform F: total momentum grows by `N_fluid · F` per step (1e-10 rel).
- feq moment identities (unit test): Σfeq = ρ, Σfeq c = ρu, Σfeq cc = ρ(cs²I + uu) to 1e-14
  at several |u| ≤ 0.1 points.
- Angle (f32): mass drift ≤ 1e-5 (10³ step); force-driven momentum-growth rel error ≤ 1e-5
  (10² step). Requires the deviation-storage scheme (PHYSICS.md 2026-07-05).

### T7. Lid-driven cavity (Ghia et al. 1982)
- Setup: BB all edges, top MovingWall{[U,0]}, Re ∈ {100, 400, 1000}, N=129 (L=N−2), U=0.1, TRT.
  Until steady (ε=1e-8 or 300k step cap).
- Accept: centerlines u(y) / v(x) vs Ghia's 17 tabulated points, RMS ≤ 0.02·U (Re=100/400),
  ≤ 0.03·U (Re=1000). Primary-vortex center within ±0.02L.
  **Known typo**: Re=400 v(x=0.9063)=−0.23827 (discontinuous with neighbors; PHYSICS.md).
  Exclude this point from the RMS.
- Angle: same solution for the 4 lid rotations under the anti-diagonal symmetry map in
  PHYSICS.md. L∞ ≤ 1e-10.

### T8. Cylinder — Schäfer–Turek 1996 benchmark
Geometry (strict): channel 22D × 4.1D, cylinder center 2D from inflow, 2D from bottom
(asymmetry triggers shedding; center y/H = 0.4878). BB top/bottom, left VelocityInlet +
parabolic `set_inlet_profile` u(y) = 4 u_max y_w(H−y_w)/H², right PressureOutlet{1}.
U_mean = (2/3) u_max, Re = U_mean·D/ν. Cd = 2Fx/(ρ U_mean² D), Cl = 2Fy/(ρ U_mean² D).

- **2D-1 (Re=20, steady)**: ref Cd=5.5795, Cl=0.0106, Δp*=2.9375.
  - D=20 (440×82, u_max=0.075, ν=0.05): Cd ∈ [5.2, 6.0], Cl ∈ [−0.05, 0.08] (staircase band).
  - D=40 (880×164, u_max=0.075, ν=0.1) #[ignore]: Cd ∈ [5.35, 5.85]; convergence trend
    |Cd(40)−ref| < |Cd(20)−ref|.
- **2D-2 (Re=100, unsteady)**: ref Cd_max ≈ 3.22–3.24, Cl_max ≈ 0.99–1.01, St ≈ 0.295–0.305.
  D=40, u_max=0.15 #[ignore]:
  - St ∈ [0.28, 0.32] (from Cl zero crossings).
  - Cd_max ∈ [3.0, 3.5], Cl_max ∈ [0.8, 1.2].
  - Consecutive Cl-period variation ≤ 2%.
- Bands are widened for staircase geometry; tighten when curved boundaries land (Phase 7).

### T9. Outflow (zero gradient) robustness
- Setup: T8-2D-2 channel with right edge = Outflow.
- Accept: no NaN/Inf over 10⁵ step; reverse-flow mass flux ≤ 5% of inflow; rms of pressure
  oscillations at x > 0.9L within **15×** of central region (zero-gradient partially reflects
  pressure waves — intrinsic; convective outlet candidate on Phase 7 backlog).

### T9b. Convective outflow (`ConvectiveOutflow`)
- Mass-consistency correction: pin edge density to adjacent cell. u_conv ∈ (0, 1] validated at
  build time (`InvalidParameter`).
- Compare to Outflow with the same geometry and T9 metrics; freeze at measured values.
  Advantage is geometry-dependent — at minimum require stability, non-divergence, reverse
  flow ≤ 5%.

### T10. Robustness / error paths
- τ ≤ 0.5, unpaired Periodic, Zou-He orthogonal-edge violation, nx < 3 → `ConfigError`.
- **Stability-limit fixed case**: cavity τ=0.51, N=128, U=0.05 (Re≈1890), TRT Λ=3/16 —
  no NaN/Inf over 10⁴ step. U=0.1 (Re≈3780) diverges (known: grid Reynolds U/ν ≈ 30
  exceeds the τ→0.5 guideline U/ν ≤ 15).
- `set_solid` on an open-boundary edge panics (by spec).
- Moving-wall / inflow |u| > MAX_SPEED (=0.3) → `ConfigError::VelocityTooHigh`. Overflow in
  `set_inlet_profile` panics.

### T11. Shan-Chen single-component multiphase
Common: `ShanChen::new(-5.0)`, ψ = 1−e^{−ρ}, τ=1 (ν=1/6), initial ρ_l=2.0 / ρ_v=0.15, each step
`sc.update_force(&mut sim); sim.step()`. Pressure compared **always with the SC EOS** (`sc.pressure`).

- **Flat interface** (64×128 periodic, 30k step):
  - Coexistence: ρ_l = 1.888 ± 2%, ρ_v = 0.1194 ± 3%.
  - |p_l − p_v| / p ≤ 1e-4.
  - Spurious max|u| ≤ 5e-3.
  - Total-mass drift ≤ 1e-10 rel.
- **Laplace** (128², R₀ ∈ {12,16,20,24}, 40k step):
  - Linearity Δp vs 1/R_fit R² ≥ 0.999.
  - Slope σ = 3.32e-2 ± 10%; each droplet's σ = Δp·R within ±5% of the slope.
  - R measured from median-density isocontour area (√(area/π)).
- **f32 angle**: flat-interface stable (no NaN, coexistence ±5%).

Frozen bands asserted in `crates/lbm-core/tests/validation_multiphase.rs`.

### T11b. Contact angle (frozen G_w regression)
- Wall-attached droplet, Periodic left/right, BB top/bottom, top wall far from droplet.
  Because this implementation uses ψ=0 for solid + a separate −G_w ψ Σw s c cohesion, **G_w=0
  is not 90°** (leans non-wetting).
- Accept: θ(G_w) monotonic; regression freeze the 3 measured angles (±8°) —
  **G_w=−1.5: 133.2°, 0: 160.4°, +1.5: 163.7°**. Method: spherical-cap fit θ = 2·atan(2h/w).
- Wetting-side narrow range known; virtual-wall-density switch (T11c) is the workaround.

Frozen bands asserted in `crates/lbm-core/tests/validation_contact_angle.rs`.

### T11c. Full-range contact angle via virtual wall density
- `ShanChen::with_wall_rho(ρ_w)`, G=−5, liquid 2.0 / vapor 0.15, 160×100, 30k step.
- Accept: monotone ρ_w ↑ → θ ↓. Regression freeze (±8°):
  **ρ_w=0.3: ~180°, 0.6: 107°, 1.0: 63° (θ<90° achieved)**.
- ρ_w=1.6 = complete wetting (film over full wall) → assert "forms a film" qualitatively.

### T12. Two-component MCMP: Rayleigh–Taylor growth
Two-stage self-consistent (`multiphase::MultiComponent`, ψ=ρ, G_ab=2.6, trace 0.05):

1. **σ_AB measurement**: droplet in B (128², r₀=24, ν=0.1, 20k step);
   p = cs²(ρ_A + ρ_B + G_ab ρ_A ρ_B); Δp · R_fit = σ_AB (≈ 2.87e-2).
2. **RT growth**: 256×256, Periodic left/right, BB top/bottom, both components bulk ρ=1;
   interface at y₀=128, perturbation a₀=6·cos(kx) k=2π/256; gravity on the heavy component
   only g_a=[0,−1e-4] (effective Atwood 0.5), ν=0.1.
   - Amplitude via **k-mode Fourier projection of column mass** (contour method glitches).
   - ln(amp) regression over amp ∈ [1, 10] → γ_fit.
   - Reference **γ_th = √(gk/2 − σ_AB k³/2 + ν²k⁴) − νk²** (tension + viscosity corrected).
   - Accept: γ_fit/γ_th ∈ [0.75, 1.25]; amp reaches ≥ 10; total-mass drift ≤ 1e-10.
- Separation smoke: half-and-half init (96², 5k step) — G_ab=2.2 separates (contrast ≥ 3),
  G_ab=1.8 mixes (≤ 1.5).

Asserted in `crates/lbm-core/tests/validation_rt.rs`.

---

### T13. Partition invariance (core V2 equivalence)
Files: `t13_split_invariance.rs` / `t13_adversarial.rs` / `examples/mpi_t13.rs`.

- **InProcess partitioning**: 1×1 vs 2×2 / 4×1 / 1×4 / (3D) 2×2×2 — fields (rho/u/all f
  planes) are a **bit match** with `assert_eq!(d, 0.0)`. Adversarial angles: obstacles /
  probes / Zou-He faces on seams, lid straddling seam, Shan-Chen ψ exchange, corner droplet.
- **T13-MPI** (mpirun -n {1, 2, 4, 8}): rank-0 gathered field has max|Δ| = 0.0 vs single-rank.
  Diagnostics (mass/momentum/probed_force/NaN count) tolerated at atol+rtol 1e-12 (field) /
  1e-11 (diagnostics) — only f64 rank-partial-sum recombination difference. Repro: `./scripts/test_mpi.sh`.
- **Known blind spot**: probe double-counting on the two_pass boundary shell (spec E8/C-2)
  is invisible to field diff; after the C-2 shell fix, a two_pass on/off probe-match test on
  a width-1 axis will be added here.

### T14. Backend equivalence (CPU vs Wgpu, `--features gpu`)
File: `t14_backend_equiv.rs` (6 configs: TGV / cavity / cylinder+probe / Zou-He /
force field / moving wall, f32).

- Rel field diff ≤ 1e-5; pressure-boundary case ≤ 1e-4 (justified by a 1-ulp control test —
  "round-off, not physics"). Diagnostic values match.
- f64 + GPU rejected at compile time (no silent degrade).
- **Known gap (spec D-3, resolved in R-Phase 2)**: CPU-relative equivalence only, so a
  spec-interpretation bug where CPU/GPU break the same way is invisible. Adds two GPU
  absolute-physics tests (TGV order ≥ 1.7, cavity Ghia RMS ≤ 0.02U — calibrated to the f32
  measurements). No adapter → skip; `LBM_REQUIRE_GPU=1` promotes to fail.

### T15. (M-C) 3D physics validation (D3Q19)
File: `crates/lbm-core/tests/t15_3d.rs`.

1. **z-invariant 2D-TGV degeneracy**: init a z-invariant 2D TGV on N×N×4 (z Periodic) and
   **field-by-field match** the identical D2Q9 scenario, f64 ≤ 1e-12. Angle: same
   degeneracy via 3D Zou-He faces (parabolic inflow + pressure outlet); face-node
   prescribed-moment exactness |u − u_bc| ≤ 1e-14.
2. **Rectangular-duct Poiseuille (exact series)**: a×b cross-section, body-force driven,
   4 half-way BB walls. u(y,z) = (g/2ν)·[series]. TRT: L∞rel ≤ 1e-3 (series truncated at
   n ≤ 99 with convergence check). Flow rate Q vs analytic within ±0.5%.
3. **Sphere drag** (staircase), Re ∈ {20, 100}. Momentum-exchange Cd vs **Schiller–Naumann**
   Cd = (24/Re)(1 + 0.15 Re^0.687) within ±10% (Re=20: 2.6095, Re=100: 1.0917).
   Blockage ≤ 3% (domain ≥ 8D), D ≥ 24 (D=24, 192×128×128, periodic sides; #[ignore]).
   Default suite: lightweight D=12 with D_h normalization, band ±15%.
   Cd averaged over a 500-step window (acoustic ripple mitigation; see TESTING_NOTES).
   **Normalization (frozen 2026-07-05)**: hydrodynamic pair Cd_h = F/(½ρU²π(r+½)²),
   Re_h = U(D+1)/ν — half-way BB wall lies a half-link outside solid cells (Ladd staircase),
   nominal-D normalization was incompatible with the D=24 band via ~+2/D bias.
4. **3D-TGV (short-time convergence order)**: classic u=(sin x cos y cos z, −cos x sin y cos z, 0),
   short time t = 0.1/(νk²). Initial decay rate matches 2ν(3k²) within ±2%; order ≥ 1.7 (N=32→64).
   **u0 coefficient frozen at u0 = 1.28e-4/N** (classic 3D TGV isn't an exact NS solution → a
   resolution-independent relative offset ≈ 0.13·u0/(νk) under diffusive scaling would
   destroy the order; PHYSICS.md).
   D3Q27 coverage is currently periodic / closed-wall only: any open face
   (velocity inlet, pressure outlet, outflow, convective outflow) is rejected
   with `UnsupportedOpenFaceLattice` because D3Q27 has 9 unknown populations
   per open face. Open-face D3Q27 validation gates are future work; current
   D3Q27 tests cover periodic TGV / split invariance and closed-wall specs
   (`crates/lbm-core/src/solver.rs::validate_lattice`,
   `crates/lbm-core/tests/t15_3d.rs`,
   `crates/lbm-core/tests/t13_split_invariance.rs`).
5. **T15.5 3D cavity** (Albensoeder & Kuhlmann 2005, Re=1000): reference in
   [T15_5_CAVITY3D_REFERENCE.md](T15_5_CAVITY3D_REFERENCE.md); file `t15_5_cavity3d.rs`.
   - Default N=64 qualitative sentinel: mass drift, symmetry-plane |v|/U, extrema
     signs/locations qualitative.
   - `#[ignore]` N=72 spec profiles: centerline u/w RMS ≤ 0.030U / 0.035U; midplane
     symmetry ≤ 1e-14; anti-2D guard. Extremum band **frozen at 0.13** (convergence-tendency
     test 2026-07-05 confirmed error decreases toward A&K with N; PHYSICS.md
     "T15.5 extremum band"). Tightening candidate after MF-α cumulant.
   - N=48 diverges (Re/(N−2) ≲ 15); N ≥ 72 for spec-grade profiles.
   - Re=100/400 on hold (originals not yet obtained).
   - **Endpoint rule**: table endpoints are boundary values themselves — sampler returns
     u(z=0)=0, u(z=1)=U, w(x=0)=w(x=1)=0 directly. Using adjacent fluid cells conflates
     the half-way moving-wall layer with reference endpoints (u-line RMS degrades to 0.037
     at N=72 otherwise).

### T16. (M-E, FP16 storage) — **implemented, gated for capacity/throughput mode**
File: `crates/lbm-core/tests/t16_fp16_storage.rs` (`#[ignore]`, feature `gpu`, requires
SHADER_F16 adapter). Implementation surface: `GpuStorage::F16` on the wgpu backend stores
distribution buffers as IEEE f16 and keeps arithmetic in f32; generated WGSL validates f16
storage for D2Q9, D3Q19, and D3Q27. The frozen accuracy bands are currently enforced on
D2Q9 wgpu f32-vs-f16 scenarios:

- **TGV 256² transient over one decay time (41,501 steps)**: f16-vs-f32 velocity L2_rel
  must be ≤ 2.0e-1; frozen measured value 1.401e-1. f16-vs-analytic velocity L2_rel
  must be ≤ 2.0e-1; frozen measured value 1.413e-1.
- **Lid cavity 128², Re=100, 40k steps**: f16-vs-f32 centerline L2_rel must be ≤
  5.0e-3; frozen measured value 2.579e-3.

Use: FP16 storage is a capacity/throughput mode, not a validation-grade long-time reference
mode. It doubles on-device distribution capacity and performance characterization records
~2.0x MLUPS at 2048² plus D3Q19 f16 > 5 GLUPS (`docs/PERFORMANCE.md`), but long transients
accumulate f16 store-rounding as a random walk on a decaying signal (`docs/PHYSICS.md`).
Use f32/f64 storage for validation-grade long-horizon references.

### T17. (M-F, VR-STR coupled multiphysics) — **mixed: landed subsystems plus pending critical path**
Source: [REQ_STIRRED_REACTOR.md](REQ_STIRRED_REACTOR.md) §8 rev.4. Tests are adversarially
authored from the REQ by orders separate from implementation.

Subsystem status synced to `REQ_STIRRED_REACTOR.md` "Landed vs. pending" and spot-checked
against code/tests where feasible:

| Subsystem | Spec status | Implementation status | Validation status |
|---|---|---|---|
| W0 / MF-alpha core basis (D3Q19/D3Q27, cumulant, Guo) | Fidelity basis | LANDED | VALIDATED for current scope (`cumulant_acceptance.rs`, T13/T15 D3Q27 periodic/closed-wall gates) |
| W-ROT rotating IBM | IBM-inertial fidelity default; MRF Phase 2 | LANDED | VALIDATED at subsystem level (`rotating_ibm.rs`, `mf_interim.rs`); full VR-STR-01 tank/PIV gate NOT YET |
| W-GRAV well-balanced gravity | Gravity axis reference | LANDED | VALIDATED for single-phase static stratification / force equivalence (`gravity.rs`, VR-STR-06 gravity axis); active-scalar 06+ NOT YET |
| W-LES turbulence SGS | WALE default; Smagorinsky separate equation | LANDED | VALIDATED/CHARACTERIZED for WALE TGV, stabilization, CPU/GPU equivalence (`wale_les.rs`, `t14_wale_gpu_equiv.rs`); Re_tau DNS acceptance NOT YET |
| W-STRESS stress fields | FR-STRESS convention fixed | PARTIAL/PENDING | NOT YET complete T17 VR-STR-03 |
| W-VOF resolved interface | Conservative Allen-Cahn free surface, fidelity default | PENDING | NOT YET; blocks VR-STR-02 and W-BCTOP/interfacial transfer |
| W-BCTOP top boundary / degassing / contact angle | Waits on W-VOF | PENDING | NOT YET |
| W-SCAL scalar ADE | Passive scalar path specified; active scalar fidelity default waits on coupling/VOF | PENDING | NOT YET; scalar total-mass and active dt-halving negative/consistency gates remain pending |
| W-REACT reaction / active feedback | Waits on W-SCAL; active mode waits on W-VOF | PENDING | NOT YET |
| W-PART / D-track particles + deposition | D-track P2 landed; higher-density four-way Phase 2 | LANDED for current D-track P2 scope | VALIDATED for current particle/deposition scope (`t18_*.rs`, `particles_deposition_smoke.rs`, `accuracy_audit_particles.rs`); full VR-STR particle coupling still incremental |
| W-BUB point bubbles + PBM | Phase 2 / API-reserved | PENDING | NOT YET |
| W-COUP / W-IO coupled loop and reactor outputs | Incremental across producing subsystems | PENDING | NOT YET full coupled VR-STR-05/07 |

**Band governance (rev.4)**: each row has a provisional numeric band (Np ±10%; PIV L2<15% /
L∞<30%; droplet mass drift <0.1%/1000 step; U_t ±10%; k_La ±25%; stratification |u|<1e-6 l.u.;
GPU/MPI mean <2% / higher-order <5%; Ca_spurious = μ_l|u|/σ < 1e-3). Tighten freely;
**loosen only with PHYSICS.md rationale** (reference uncertainty / method order / resolution
limit — the T15.5 precedent).

Mandatory negative/consistency tests (rev.4): forcing-moment sign negative test; stress
stage cross-check; J_ρ consistency + advected-droplet conservation; sparger gas-phase (φ=0)
unit test; phase-wise scalar total-mass conservation; active-scalar dt-halving convergence.
Each test row records: metric / target·reference / tolerance / resolution / time window /
backend / pass-fail rule.

| ID | Target | Acceptance criteria (fixed in REQ) | Band status |
|---|---|---|---|
| VR-STR-01 | Single-phase stirring (baffled tank, non-aerated) | Rushton Np = experimental correlation ratio; impeller-discharge PIV/LDA reference lines by L2/L∞rel. Np = P/(ρ_l N³D⁵), P = Ω T_q (no 2π double-count) | tolerance frozen after impl |
| VR-STR-02 | Gas-liquid (02a/02b/02c) | **02a single bubble**: U_t vs Grace (Eo-Mo-Re). **02b swarm**: ε_g spatial dist, swarm rise velocity, d_32, ν_t under BIT. **02c aerated stirring**: ε_g, d_32, k_L a ratios | frozen after impl |
| VR-STR-03 | Shear / stress field | MMS single-phase, curved Couette, rotating cylinder, non-Newtonian Poiseuille, multiphase static droplet. **Ca_spurious < 10⁻³ (fixed)**. Near-wall L∞ survey lines | order / L2/L∞ frozen after impl |
| VR-STR-04 | Scalar / reaction | Taylor–Aris dispersion, reaction-diffusion at known Da, k_L a (formula stated). Per-test target Pe/Da/Sc. SGS scalar uses Sc_t=0.7 | frozen after impl |
| VR-STR-05 | Coupled regression / conservation | probe_state_hash bit-equiv **is single-backend regression only**. Individual drift thresholds for mass, momentum, total scalar, gas-phase volume, particle count, and energy-like monitored quantities (kinetic E, interfacial free E, particle kinetic E). GPU/MPI tolerance-based | frozen after impl |
| VR-STR-06 | Well-balanced hydrostatics | |u| < ε in static stratification (density ratio 10³). **06+**: active ON with C≡C_0 → same static (F_b^scalar → 0); with ∇σ=0 → matches constant-σ reference | frozen after discretization decided |
| VR-STR-07 | Initialization independence | quasi-steady statistics match within threshold under varied run-up / stats start | window/threshold frozen after impl |
| VR-STR-RELAX | Relaxation-mode equivalence | each relaxation vs its fidelity reference — MRF→IBM (Np, lines, torque) / point-bubble→resolved (ε_g, d_32, k_La, budget) / one-way→two-way (particle stats, mass-loading cap) / AMR→uniform (conserved, interface, coarse-fine budget) / aggressive f32→fidelity profile (drift, Ca_spurious, Np, curvature) | frozen when the relaxation extension lands |

### T18. (D-track, dispersed-phase deposition) — **spec wired, P2 in progress**
Source: [DISPERSED_DEPOSITION.md](DISPERSED_DEPOSITION.md) (frozen; T18.1–.3 public surface
frozen in §5). Vocabulary stays domain-neutral (particle-laden multiphase flow). Tests are
adversarial and dispatched separately from implementations.

- **T18.1 Localized interior volume source/sink (CR-1)** — `tests/t18_1_interior_source_sink.rs`:
  - Point-like sink in a filled closed box (D3Q19, all faces Closed) reproduces
    incompressible sink far-field `u_r(r) = q/(4πr²)` on shell-averaged radial profiles for
    r away from region and walls; band ±10%.
  - Mass ledger: d(total_mass)/step = Σ q_lu to summation round-off. **Band 1e-6 rel**
    (achievable error is bounded by N·ε·M/|Σq|; PHYSICS.md 2026-07-06 T18 reconciliation).
  - `Jet` source delivers prescribed momentum flux within ±2% (measured inside the
    pre-wall-contact window — after the acoustic front hits walls, BB absorbs momentum by
    construction).
  - Errors: region touching a face, overlapping solids/sources, |u| > MAX_SPEED, sink strong
    enough to drive local ρ ≤ 0.
  - T13 bit match with the source region on a seam.
  - GPU: a spec with `sources` is rejected with `SpecError` (no silent physics).
- **T18.2 Per-cell masked face BC (CR-2)** — `tests/t18_2_masked_face.rs`:
  - Impinging jet: central Velocity patch (down) + coaxial Pressure annulus (4 rects) on the
    same top face; closed floor and sides. Spin-up 2400 steps, then mass drift ≤ 5e-9 rel /
    200 steps. Floor wall jet is radial — stagnation on axis (≤ 0.25× peak), off-axis peak
    > 1e-5, monotone decay past the peak.
  - Frozen semantics (PHYSICS.md 2026-07-06): patch rects are global in-face coords
    (seam-safe per subdomain); non-patch cells of a bare Closed base face with patches form
    a zero-velocity lid; a Closed patch on an open base is a lid on its rect.
  - Errors: patch out of face bounds, overlapping patches, one-open-axis rule violated by
    base ∪ patch union.
  - T13 with a patch straddling a seam; GPU rejection.
- **T18.3 Dispersed-phase deposition layer (CR-3)** — `tests/t18_3_particle_deposition.rs`:
  - Terminal velocity: single particle under gravity matches Schiller–Naumann (reducing to
    Stokes v_s = (ρ_p/ρ_f − 1) g d²/(18ν) at Re_p < 0.1) across a St sweep ≥ 2 decades; ±2%.
  - Deposition-map determinism: bit-identical deposit records (positions and order) between
    unpartitioned vs T13-partitioned runs of the same scenario.
  - Floor-crossing capture: straight-line trajectory recorded at exact interpolated
    crossing point; particle leaves suspended set (n_deposited + n_suspended conserved).
- **T18.4 Forward-model monotone anchors** — example-level, `#[ignore]` (runs
  `--example dispersed_seeding`): (1) ejection rate ↑ → CV ↑; (2) agitation → CV ↑ vs
  quiescent (quantitative after P3; stub until then); (3) fill volume ↑ → CV ↑; (4) repeated
  passes ↑ → CV ↓. Regression: gentle CV ∈ [1.05, 1.30], empty_bin_fraction ≤ 0.15;
  Ma ≤ 0.3 / τ ≥ 0.51 aborts.
- **T18.5 Inverse recovery (P4)** — **stub**: inverse solver recovers a known-good recipe on
  a synthetic n\* generated by the forward model itself (parameter recovery within tol,
  objective ≤ frozen band).

## Test implementation conventions (for codex)

- Location: `crates/lbm-core/tests/validation_*.rs` (one theme per file).
- Shared helpers in `crates/lbm-core/tests/common/mod.rs`.
- `#[ignore]` on heavy computations (T7 Re=1000, T8, T9). CI-equivalent: `cargo test --release`;
  full: `cargo test --release -- --include-ignored`.
- No randomness (deterministic). Each `assert!` carries the measured value
  (e.g. `assert!(err < 5e-3, "L2rel = {err}")`).
- External crates allowed only for `approx` (anything else → consult).
- Embed the Ghia reference as a constant table inside the test file (with source comment).

## Construction-rejection semantics (R-Phase 1, A-2..A-7 — merged 2026-07-05)

Invalid configurations fail loudly at build/placement time instead of running silently
non-physically (E2/E4/E5b/E6/E7 pathologies are all construction errors). New surface:
`SpecError` (native `GlobalSpec::validate`, enforced in `Solver::build`; scenario 3D routes
through it), `ExchangeScope` (HaloExchange scope guard), `Diverged` + `run_guarded` (runtime
NaN watchdog, CPU/GPU/MPI), `Simulation::set_solid_allowed` (non-panicking placement probe;
wasm paint uses it), `params::MAX_SPEED` relocation (compat re-exports). Native
`Solver::init_with` / GPU seed validation is deliberately still open (S2) — scheduled into
R-Phase 2.
