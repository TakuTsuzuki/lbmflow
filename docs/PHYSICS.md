# PHYSICS.md — physics models and load-bearing decisions

Purpose: index the physics live in `crates/lbm-core` today, and record ONLY
the decisions whose reasoning is not derivable from the code + test suite in
front of you. This is not a change log. If a "lesson" is now an assertion in
a test or an obvious LBM fact, it does not belong here.

**When to add an entry**: you changed a physical behavior AND a future reader
staring at the code alone could not reconstruct why. Otherwise, put the
rationale in the test's assert message or a source comment and move on.

---

## 1. Current physics stack

Each item is grep-checked against `crates/lbm-core`. One-line "why" per item.

- **Lattices**: D2Q9, D3Q19, D3Q27. `cs² = 1/3`, `τ = 3ν + 0.5`. Single
  source of truth for direction ordering: the `Lattice` trait impl in
  `lattice.rs` (see CLAUDE.md invariant).
- **Collision — default TRT, magic `Λ = 3/16`**. `Λ = (1/ω⁺−½)(1/ω⁻−½)`.
  Why default: half-way straight-wall placement is exact for Poiseuille
  (`smoke_poiseuille.rs::trt_magic_is_exact`), which BGK cannot achieve.
- **Collision — BGK** available as a lighter alternative. Same convergence
  order as TRT away from walls; loses the Poiseuille exactness.
- **Collision — cascaded central-moment** for D3Q19/D3Q27,
  `CollisionKind::CentralMoment { omega_shear }`. Implemented as a
  cascaded central-moment operator (not logarithmic cumulants). D3Q19 uses
  the D3Q27 tensor-product basis with the eight `x·y·z` corner moments
  dropped, matching the missing body-diagonal populations. Targets are the
  DISCRETE second-order Hermite equilibria used by BGK/TRT, not continuous
  Maxwellian moments (see §2 "Cumulant target choice").
- **Distribution storage — deviation form `f_q − w_q`**. Rest state is the
  zero vector; f32 mantissa is spent on the fluctuation, not on the rest.
  This is what lets f32 hit validation grade (see §2 "Deviation storage").
- **Guo forcing with `F/2` correction in `u`**. Physical velocity
  `u = (Σf c + F/2)/ρ` is used in `feq`, in the source term, and in every
  velocity accessor. Every backend enforces this; CLAUDE.md invariant.
- **Walls — half-way bounce-back**. 1-cell solid rim; wall surface midway
  between rim center and fluid center.
- **Moving walls — momentum injection** `+6 w_q ρ (c_q · u_w)` on the
  bounced link. Corner cells: "faster wall wins" (see §2 "Rim-corner
  orientation").
- **Virtual wall density for Shan-Chen contact angle** (`ShanChen::with_wall_rho`).
  Cohesion sum treats solid neighbours as if they carried ψ(ρ_w).
  Why this scheme: the older `g_wall` variant could only reach θ ≥ 133°;
  virtual wall density unlocks the full 0–180° range on the same solver.
- **Open BC — Zou-He, single implementation** parameterized by face normal
  `(n, t)`. Handles velocity and pressure faces on all four (2D) / six (3D)
  faces via one formula; with body force, it closes on the Guo raw momentum
  `rho u - F/2` so velocity faces prescribe physical velocity. See doc
  comments in `params.rs` / `kernels.rs`.
- **Open BC — convective outflow with mass-drift correction** (`params.rs`,
  `kernels.rs::convective_outflow`). The naive `f = (f_prev + λ f_int)/(1+λ)`
  drifts to NaN; the live implementation re-pins edge density to the
  neighbouring cell and redistributes to unknown populations by weight (see
  §2 "Convective outflow mass pinning").
- **Multiphase SCMP — Shan-Chen single component**, classic and exponential
  ψ, Carnahan-Starling (CS) EOS helpers, wall adhesion via `g_wall` OR
  virtual wall density.
- **Multiphase MCMP — Shan-Chen two component** with cross repulsion
  `−G_ab ψ_A Σ w ψ_B c` applied action-reaction per link (total momentum
  conserved); per-component gravity.
- **Density-ratio ceiling from CS-EOS**. The SC EOS helpers set the
  physically achievable density ratio; higher ratios need a different EOS,
  not tuning of `G`.
- **LES — WALE default** (not Smagorinsky). Nicoud & Ducros (1999), `Cw=0.325`.
  Why default: `S^d:S^d ≡ 0` analytically in pure shear, so `ν_t = 0` in
  laminar Couette/Poiseuille. Smagorinsky would silently bend the baseline
  physics on any resolved laminar shear.
- **Per-mass gravity** `F_g(x) = ρ(x) g` on fluid cells, distinct from the
  uniform force-density body force. Why: a constant force density is
  cancelled exactly by hydrostatic pressure and cannot express buoyancy;
  `ρ·g` creates the required imbalance for Shan-Chen bubbles/blobs.
- **Well-balanced gravity composition point**. Gravity is not a separate
  forcing scheme — it is composed at the solver's single Guo source point
  as `F_total = F_user + F_cell + ρ·g`. This is the line W-VOF will later
  swap `ρ` for `ρ(φ)` without touching collision or BC kernels (see §2
  "W-GRAV composition point").
- **Bouzidi 2nd-order interpolated bounce-back** for curved walls
  (`bouzidi.rs`), from Bouzidi-Firdaouss-Lallemand 2001. Moving walls add
  the Ladd momentum source at every interpolated reflected point; for the
  qd < 1/2 two-fluid-node branch this means the first and second interpolation
  points together contribute `2 w_q rho (c_opp·u_w)/cs²`.
- **Rotating bodies** — two paths: volume penalization (`F = 2ρχ(u_target − u*)`,
  algebraic no-overshoot at χ=1 with a finite spin-up ramp) and marker-based
  direct-forcing IBM (Uhlmann sequence + Wang multi-direct-forcing
  correction). Both feed the existing per-cell force field, so they only
  enter the solver through the Guo path.
- **Localized volume sources and face patches**. Volume sources apply
  after the open-BC pass and before moment recomputation. MassFlow uses
  `Δf_q = w_q q_cell` (zero first moment); Jet uses the equilibrium-shaped
  `Δf_q = w_q q_cell · (1 + 3 c·u + 4.5 (c·u)² − 1.5 |u|²)` (`Σ = q_cell`,
  `Σ c = q_cell u`). Face patches override the base face BC only inside
  their in-face rectangle; a Closed base face with a patched region is
  legal (bare cells become a zero-velocity Zou-He lid — see §2 "Patched
  Closed face is a lid").
- **Dispersed-phase deposition closures** (D-track): adhesion capture and
  resuspension per `docs/DISPERSED_DEPOSITION.md`. Both are literature-backed
  closures with recorded derivations and their own validation tests (T18
  family); revisions to either must land through PHYSICS.md.

---

## 2. Load-bearing decisions

### 2026-07-07 Bouzidi qd<1/2 moving-wall second-point correction (`bouzidi.rs::apply_bouzidi_impl`)
- Form: for `qd < 1/2` with a second fluid node,
  `f_opp = sigma*(f_q(x_f)+W) + (1-sigma)*(f_q(x_f-c_q)+W)`,
  where `sigma = 2 qd` and
  `W = 2 w_q rho (c_opp·u_w)/cs² = 6 w_q rho (c_opp·u_w)` for D2Q9/D3Q19/D3Q27.
- Source: Bouzidi-Firdaouss-Lallemand (2001) interpolated bounce-back with
  the same moving-wall equilibrium shift as the existing half-way moving-wall
  rule. The old implementation applied `W` only at the first interpolation
  point, imposing `sigma*u_w` instead of `u_w`.
- Validity domain: Bouzidi links with `0 < qd < 1`, local density from the
  adjacent fluid cell, and wall speeds inside the existing `MAX_SPEED` limit.
  The correction is active only for moving walls; static-wall behavior is
  byte-identical because `W = 0`.
- Validation:
  `crates/lbm-core/tests/accuracy_audit_bouzidi_moving.rs::qd_sweep_moving_wall_couette_should_match_offgrid_linear_profile_all_qd`
  checks `qd={0.25,0.5,0.75}` against the off-grid Couette line with
  `max_rel_dev <= 2e-3`, and
  `::qd_half_moving_wall_is_bitwise_half_way_moving_wall` preserves the
  `qd=0.5` half-way moving-wall bitwise degeneracy.
- Replaces / interacts with: fixes ANOM-P4-025 without changing collision,
  static Bouzidi walls, or the half-way moving-wall branch.

### 2026-07-07 behavior review — `vv_sedim_2d` rev 2 attempt
Pattern: with `d=6`, `|g|=1e-4`, and `30_000` steps, the run deposits all
500 particles but reports `mean_deposition_x=5.005724e2`, outside the 128-cell
basin; the histogram piles 444 counts into the final clamped bin.
Mechanism: particles that pass the pressure outlet continue to be sampled from
the clamped boundary state and later cross the deposition plane out of domain,
so the deposition map is dominated by post-outlet transport rather than in-basin
settling.
Resolved vs closure: Stokes/SN particle settling remains the active validated
closure (`v_s=1.2e-3`, `Re_p=4.32e-2`), but the headline spatial pattern is an
example/domain bookkeeping artifact from combining an open outlet with clamped
particle sampling.
Artifacts checked: `out/vv_sedim_2d/density_00000.png`,
`out/vv_sedim_2d/density_15000.png`,
`out/vv_sedim_2d/density_30000.png`, and
`out/vv_sedim_2d/deposition_map.png`.
Verdict: ARTIFACT.
Routing: PM/spec decision required before claiming the rev-2 behavior anchor:
add a physically stated particle outlet/escape rule, revise the mean-x anchor,
or change the protocol so particles deposit inside the basin without relying on
post-outlet clamped sampling.

### 2026-07-07 behavior review — `vv_sedim_2d` rev 3 closed basin
Pattern: in the 128x64 all-bounce-back box, 500 particles seeded uniformly on
the horizontal line `x=10..118, z=60` settle to the bottom without lateral drift:
`deposition_fraction=1.000000e0`, `mean_deposition_x=6.400000e1`,
raw `deposit_x_std=3.123933e1` from the seed-line span, and seed-relative
`lateral_scatter_std=0.000000e0`.
Mechanism: with no fluid crossflow and no pressure outlet, the sampled fluid
velocity stays zero and the only active acceleration is vertical particle
gravity, so each seed column follows the same low-Re settling trajectory.
Resolved vs closure: the closed fluid box and zero velocity field are resolved
LBM behavior; particle motion uses the existing one-way Schiller-Naumann/Stokes
low-Re closure at `v_stokes=1.200000e-3` and `Re_p=4.320000e-2`.
Artifacts checked: the rev-2 outlet/clamped-boundary artifact is removed by the
closed box. The deposition counting plane remains `z=1.0`, the first fluid row,
so floor crossings are recorded before the staircase solid-contact model pins a
particle to the bottom rim; this affects event bookkeeping only, not lateral
transport. Visual artifacts:
`out/vv_sedim_2d/density_00000.png`,
`out/vv_sedim_2d/density_30000.png`,
`out/vv_sedim_2d/density_60000.png`, and
`out/vv_sedim_2d/deposition_map.png`.
Verdict: PHYSICAL.
Routing: none.

### 2026-07-06 dispersed seeding closure removal (`examples/dispersed_seeding`)
- Form: removed the example-local harshness switch, analytic jet/wall-jet
  superposition, stochastic lateral dispersion, side-wall particle clamps,
  direct agitation kicks, reservoir scoring heuristic, and mystery reservoir
  force. Tray particles now advance one step after each tray LBM step and
  trilinearly sample the live resolved velocity field. Translational agitation
  uses the non-inertial-frame pseudo-acceleration `-A omega^2 sin(omega t)` on
  the fluid via the core per-mass Guo forcing path; particles receive the
  matching pseudo-force through `ParticleSet::g`, divided by
  `(1-rho_f/rho_p)` because the core applies the same buoyancy weighting as
  gravity.
- Source: governing-equation frame transformation for a translating frame plus
  existing T18.3 Stokes/SN particle settling validation. Reservoir withdrawal
  samples the 1D concentration at `depth_frac` by backtracing
  `z0 + v_s(d)t` into the initially filled column.
- Validity domain: one-way Lagrangian particles; Stokes/SN particle model
  validity from T18.3; resolved tray flow only, with no unresolved turbulent or
  wall-jet closure.
- Validation: `cargo build --release -p lbm-cli --example dispersed_seeding`
  passed. Gentle resolved-only run completed with `Ma=0.093`, `tau=0.536`,
  `n_deposited=0`, `n_suspended=10000`, `n_extracted=10000`, wall time
  `1607.25 s`.
- Replaces / interacts with: replaces the closure-driven P1.1/T18.4 example
  path. The old CV/empty-bin bands are invalid after this removal.

### 2026-07-06 behavior review — dispersed seeding gentle resolved-only
Pattern: no deposition map forms; all extracted particles remain suspended.
Mechanism: the resolved tray field and protocol duration do not carry particles
from the nozzle region to the floor after the closure layer is removed.
Resolved vs closure: the reported pattern is resolved-only; no example-local
transport closure, clamp, or kick remains live.
Artifacts checked: `out/dispersed_seeding/gentle/density.csv`,
`out/dispersed_seeding/gentle/density.png`,
`out/dispersed_seeding/gentle/tray_velocity.vtk`, and
`out/dispersed_seeding/gentle/near_floor_radial_velocity.csv`. The former edge
ring does not survive because there are no floor deposits. The near-floor
radial velocity profile is weak, order `1e-6 m/s`, and does not support a
resolved edge-ring deposition mechanism over this sample duration.
Verdict: UNKNOWN / CAPABILITY GAP.
Routing: PM/core decision required: revise the demo budget/protocol, improve
resolved jet/free-surface capability, or introduce a literature-backed closure
with the full provenance and validation package.

## 2026-07-06 cumulant track stage 2: CPU central-moment reference

Stage 2 implements `CollisionKind::CentralMoment { omega_shear }` as a cascaded
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

On 2026-07-07 the central-moment transform was factored algebraically without
changing the operator. The previous implementation built the shifted
central-moment matrix `M(u)` and inverted it per cell. The live CPU and GPU
paths now compute fixed raw moments `R f`, apply the binomial raw-to-central
shift by `u`, keep the same relaxation/source schedule, inverse-shift
central-to-raw moments, and multiply by the fixed `R^-1` matrix for the
lattice. This is the identity
`(c - u)^e = sum_{k <= e} binom(e, k) c^k (-u)^(e-k)` and its inverse; it
adds no model term, limiter, calibration constant, or changed validity domain.

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
corner content. ANOM-P4-008 later showed that the D3Q19-only `+0.0025`
shear-rate offset was a banned finite-resolution calibration and removed it.
The finite-frame cubic-velocity term remains pending as a compile-visible
ablation target:
`omega_eff = omega_shear * (1 - 0.16 |u|^2)`, clamped to the valid range,
unless `CENTRAL_MOMENT_DISABLE_VELOCITY_CORRECTION_FOR_ABLATION` is set for
the E1 rerun. Here `u` is the same physical velocity used for equilibrium
and forcing. No regularization, positivity filter, or entropic limiter is
active in this stage; validation therefore uses the explicit range
`0 < omega_shear <= 2`.

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


Chronological. Each entry states what future readers of the code alone
would get wrong without it.

### 2026-07-04 — BGK steady-state ε floor is ~1e-11, not 1e-13
Poiseuille `4×10`, BGK, τ=0.8: after physical steady state at ~8500 steps,
`dmax/umax` oscillates permanently at ~1e-12 (never reaches 1e-13). TRT
reaches an exact discrete fixed point. The T-tests use ε = 1e-11 for BGK
steady-state gates for this reason (`validation_conservation.rs`,
`accuracy_audit_probe.rs`). Without this note, someone tightening a BGK
steady-state gate to 1e-13 would chase a rounding plateau forever.

### 2026-07-04 — TRT magic `Λ = 3/16` Poiseuille exactness
Measured `L∞_rel < 1e-10` for H=8, τ=0.8, body-force driven — matches the
analytic magic property. This is WHY TRT is the default collision. BGK has
finite τ-dependent slip under the same conditions; T2 only asks BGK for
2nd-order convergence. Pinned in `smoke_poiseuille.rs::trt_magic_is_exact`.

### 2026-07-04 — TGV requires pressure-consistent initialization
Uniform ρ=1 injects an O(u₀) acoustic residual that decays too slowly and
contaminates the velocity field (measured error ≈ 0.30/N − 0.7/N²,
first-order dominant). Initializing with
`ρ = 1 − (3u₀²/4)(cos 2kx + cos 2ky)` restores clean 2nd order
(e32=2.62e-3, e64=6.98e-4, e128=1.78e-4). `init_with` also adds the
Chapman-Enskog non-equilibrium term for the same reason. T1 spec-locks
this initialization and the diffusive scaling `u₀ = 1.28/N`. Without this
note, someone "simplifying" the TGV init to uniform ρ silently breaks the
convergence test.

### 2026-07-05 — Rim-corner orientation: "faster wall wins"
`build_rims` used to paint edges in bottom→top→left→right order, so the
corner cell's `wall_u` was last-write-wins and the physics rotated with
orientation (top-lid cavity had a stationary corner; left-lid had a moving
one). Fix: pick the edge with larger |u|. Post-fix the engine is exact
lattice-equivariant (`L∞ ~ 3–4e-16` for anti-diagonal mirror, +90°
rotation, and diagonal mirror on a 2000-step Re=100 cavity). Correct
symmetry maps are recorded in the cavity/equivariance tests; the wrong
maps used by an early adversarial batch cost hours of "engine bug" hunt.

### 2026-07-05 — Ghia Re=400 `v(x=0.9063) = −0.23827` is a typo
The neighbouring points are −0.44993 (0.8594) and −0.22847 (0.9453); our
solution matches smoothly at −0.37657. The upstream gist maintainer also
notes the entry as suspect. T7 excludes this one point from the Re=400 RMS
(`validation_cavity.rs` line ~133, with the point held in `GHIA_X` so the
exclusion is transparent). Without this note, a well-meaning cleanup that
re-includes the point pushes T7 RMS ~5× over the band.

### 2026-07-05 — Deviation storage `f − w` unlocks f32 validation grade
Internal representation stores `f_q − w_q`; the rest state becomes exactly
zero and f32 mantissa is spent on the fluctuation. Only four points in the
code need to know: the deviation-form `feq` written in `δρ = ρ − 1`; the
`ρ = 1 + Σdev` moment; the `+1` in the Zou-He closure (`Σw = 1` on any
straight edge); and the `+cell_count` in the mass aggregator. Measured
effect: f32 uniform-force momentum error 1.34e-3 → 2.8e-7 (~4800×); f32
TGV L2 (N=64) 7.1e-4 ≈ f64 7.0e-4. Rationale for the whole `real.rs` /
storage layout is here; there is no cheaper way to reach f32 validation
grade on any lattice this size.

### 2026-07-05 — Convective outflow needs mass-consistency pinning
Naive `f = (f_prev + λ f_int)/(1+λ)` (previous-step values in the pull
scheme's unknown slots) NaNs over long runs, because independent
relaxation of unknown distributions grows a drift mode. Live implementation
sets the edge density to the neighbouring cell density each step and
distributes it to the unknown populations in proportion to weights. Healthy
over 34k steps of a Kármán wake. Reflection reduction vs zero-gradient
outflow is geometry-dependent (probe_phase8: 0.72 vs 0.97), so the choice
between them is a scenario decision, not "convective is strictly better."

### 2026-07-05 — Full contact-angle range via virtual wall density
Measured with `G = −5`: `ρ_w` 0.3 → θ ≈ 180°; 0.6 → 107°; **1.0 → 63°
(θ < 90° achieved)**; 1.6 → complete wetting. The old `g_wall` scheme only
reaches θ ≥ 133°. This is the reason `ShanChen::with_wall_rho` exists as
a separate constructor and why T11c is the acceptance gate.

### 2026-07-06 — Per-mass gravity + W-GRAV well-balanced composition point
Two entries share a rationale:
1. **Per-mass form** `F_g = ρ·g` (not constant force density): a constant
   force density is exactly cancelled by hydrostatic pressure and cannot
   express buoyancy. Per-mass creates the imbalance so Shan-Chen bubbles
   rise and heavy blobs sink under the same downward `g`.
2. **Composition point**: gravity is added into the existing Guo source
   as `F_total = F_user + F_cell + ρ·g` at the solver's single one-step
   staging line. This is deliberate: W-VOF later swaps `ρ` → `ρ(φ)` here,
   and the hydrostatic residual `(ρ(φ) − ρ_h) g + hydrostatic_ref` enters
   at the same line — collision and BC kernels never change. Single-phase
   freezes `ρ_h = 0` so `set_gravity(g)` is bit-identical to a raw
   per-cell force field filled with `ρ(x) g`.

Without this note, someone splitting gravity into its own forcing scheme
breaks the W-VOF handoff plan silently.

### 2026-07-06 — T15.5 3D-cavity extremum band 6% → 13% at N=72
Re=1000 cubic-cavity centerline extrema at N=72 sit 9.1–10.5% shallow of
Albensoeder & Kuhlmann (2005) spectral, while profile RMS bands pass with
2× margin and extremum positions are within half a cell. Evidence this is
resolution (not an engine bug): N=64→72 convergence-tendency test PASSES
(error decreases with N); global profile shape matches (RMS with 2× margin);
N=48 diverges to NaN exactly where `Re/(N−2) ≲ 15` predicts (so N cannot
be lowered); independent 3D gates are tight (TGV3D order 1.91, duct
`L∞_rel` 2.3e-4, sphere drag +0.6%). Decision: freeze the N=72 extremum
band at 0.13; positions and RMS unchanged; the convergence-tendency test
is the guard that the gap closes with N. Loosening this band again
requires a fresh entry.

### 2026-07-06 — Rotor penalization: default `χ = 1`, ramp = 200 steps
The most aggressive stable configuration from the stability envelope
(`χ = 1`, ramp = 200) is the default because it algebraically pins blade
cells to solid-body rotation (`u_phys = u_star + χ (u_target − u_star)`)
with no overshoot. `χ = 1` with no ramp is unstable at Ma_tip=0.2; hence
the ramp is not zero. This retires the earlier F4 empirical force cap —
stability comes from the implicit-style Guo force balance, not from
clipping. Torque convention: `Rotor::torque()` is REACTION torque on the
rotor (`Σ r × (−F)`).

### 2026-07-06 — FP16 storage is a capacity/throughput grade, NOT a
### long-transient accuracy grade
`GpuStorage::F16` stores deviations in IEEE f16; arithmetic stays f32
(loads widen, stores narrow). Measured: steady flows re-converge to the
f32 answer (cavity 128² Re=100 centerline L2_rel = 2.579e-3, band 5e-3);
long transients accumulate rounding as a random walk on a decaying signal
(TGV 256² over one decay time: L2_rel 1.401e-1 vs f32, 1.413e-1 vs
analytic, band 2e-1). The transient bound is
`~5e-4 rel/step · √N_steps` against an e⁻¹-decayed signal — this is
storage-rounding-dominated, not a scheme defect. Roadmap option if a
long-transient f16 is ever needed: FluidX3D's shifted-exponent FP16S/FP16C;
NOT currently implemented.

### 2026-07-06 — Cumulant target choice: DISCRETE second-order Hermite
An earlier stage-2 draft relaxed third/higher central moments toward
CONTINUOUS Maxwellian central moments. That was wrong: the solver
initializes and equilibrates with the engine's discrete second-order
Hermite populations, and D3Q19 cannot represent the D3Q27 `x·y·z` moment
family at all. Mixing continuous higher targets with the discrete
equilibrium inflated the advected-TGV Galilean defect and made the D3Q19
decay rate lattice-dependent. Live implementation transforms the SAME
discrete equilibrium used by BGK/TRT into the central-moment basis and
uses those as relaxation targets. The D3Q19-only shear-rate offset
(+0.0025 rel) once recorded here was removed by ANOM-P4-008 as a banned
finite-resolution calibration; the remaining finite-frame cubic-velocity
term `ω_eff = ω_shear · (1 − 0.16 |u|²)` is pending ablation.
Anyone "cleaning up" this to continuous Maxwellian will regress the T15
family; this note is the reason not to.

### 2026-07-06 — Patched Closed face becomes a zero-velocity Zou-He lid
When a base Closed face has open patches, cells outside every patch were
initially left with no boundary treatment — inbound populations undefined
after streaming, T18.2 impinging-jet diverges at ~1.7k steps. Frozen
semantics: those cells get a zero-velocity Zou-He (impermeable no-slip
lid). Same rule for a Closed patch on an open base face. Rim-covered cells
are unaffected (kernel skips solids). Without this note, someone
"simplifying" the patch pass to skip non-patched cells silently NaNs the
impinging-jet class of scenarios.

### 2026-07-06 — T18.1/T18.2 relative bands are cancellation-bounded
- **Mass ledger**: 1e-12 → 1e-6 rel/step. `dm` is the difference of two
  `O(N_cells)` sums; cancellation bounds the achievable error near
  `N · ε · M / |Σq|` ≈ 3e-7 for the 18³ case; measured 3.5e-8. The 1e-12
  provisional was dimensionally impossible for f64.
- **Jet momentum** is only measurable in the pre-wall-contact window
  (~8 steps for the 24³ box); after the acoustic front hits, bounce-back
  absorbs momentum by construction. Within the clean window, injection
  delivers `q·u` within 2%.
- **Impinging-jet wall-jet peak floor** 1e-4 → 1e-5 (measured 2.4e-5 at
  `Re_jet=5`): the jet diffuses over the 22-cell drop; the floor only
  needs to separate a coherent radial outflow from round-off noise, which
  sits 7 orders lower.

### 2026-07-07 D3Q27 open-face velocity/pressure closure (`kernels.rs::zou_he_face_d3q27`)
- Form: D3Q27 velocity inlet and pressure outlet use non-equilibrium
  bounce-back on the nine incoming face links. For a face with inward normal
  `n` and tangent axes `t1,t2`, unknown links are
  `U = {q: c_q·n = 1}`. The density/normal-velocity closure is
  `rho (1 - u_n) = S0 + 2 S-`, with `S0 = sum_{c·n=0} f_q`,
  `S- = sum_{c·n=-1} f_q`; deviation storage adds the same analytic `+1`
  constant used by the D2Q9/D3Q19 closure. Velocity faces prescribe
  `u_n,u_t1,u_t2` and solve `rho`; pressure faces prescribe `rho`, solve
  `u_n`, and set `u_t1 = u_t2 = 0`. Each unknown is reconstructed as
  `f_q = f_opp(q) + 6 w_q rho (c_q·u) + delta_q`, with
  `delta_q = C1 w_q c_q,t1 / S1 + C2 w_q c_q,t2 / S2`,
  `Sk = sum_{q in U} w_q c_q,tk^2`, and
  `Ck = rho u_tk - Q0_k - 6 rho u_tk Sk`,
  `Q0_k = sum_{c·n=0} c_tk f_q`. D3Q27 tensor-product symmetry makes
  `sum_U w c_t1 c_t2 = 0`, so these two correction components do not
  cross-couple. The correction has zero mass and zero normal moment because
  `sum_U w c_tk = 0`; it supplies exactly `Ck` to the corresponding tangent
  moment. Therefore the post-closure node satisfies `rho` and all three
  velocity components to rounding.
- Source: Zou and He, "On pressure and velocity flow boundary conditions for
  the lattice Boltzmann BGK model" (Physics of Fluids, 1997; arXiv
  `comp-gas/9611001`, https://arxiv.org/abs/comp-gas/9611001) introduced
  the non-equilibrium bounce-back pressure / velocity construction. Hecht
  and Harting, "Implementation of on-site velocity boundary conditions for
  D3Q19 lattice Boltzmann" (J. Stat. Mech. P01018, 2010; arXiv `0811.4593`,
  https://arxiv.org/abs/0811.4593) give the D3Q19 on-site moment closure.
  The D3Q27 addition is the algebra above: the D3Q19 tangent deficit
  correction is distributed over the full nine-link D3Q27 incoming plane by
  the lattice weights, including the four body-diagonal corner links.
- Validity domain: planar axis-aligned D3Q27 velocity and pressure faces on
  the CPU backends, under the existing one-open-axis rule and low-Mach
  prescribed-speed guard (`MAX_SPEED = 0.3`). It does not implement D3Q27
  `Outflow` or `Convective` faces, and GPU/WGSL still rejects D3Q27 open
  faces explicitly.
- Validation: `crates/lbm-core/tests/d3q27_open_bc.rs`:
  `d3q27_open_faces_enforce_velocity_and_pressure_moments_all_orientations`
  measured max velocity-face moment error `6.939e-18`, pressure density error
  `0`, and pressure transverse velocity error `0` across all six faces;
  `d3q27_open_duct_matches_series_shape_and_d3q19` measured duct
  flux-scaled profile L2rel `2.143e-4`, unscaled L2rel `7.416e-3`
  (compressible mass-flux scaling `1.007413`), D3Q27-vs-D3Q19 L2rel
  `3.421e-4`, mass-flux imbalance `4.212e-5`, cross-flow ratio
  `3.588e-7`, and monotone pressure drop; `t13_d3q27_open_duct_split_invariant_with_bc_seams`
  passed bit-exact split invariance with seams crossing inlet/outlet cells;
  `d3q27_unimplemented_open_face_kinds_are_rejected` preserves explicit
  rejection for unimplemented D3Q27 open kinds.
- Replaces / interacts with: lifts the previous `UnsupportedOpenFaceLattice`
  restriction only for D3Q27 velocity inlet / pressure outlet. D2Q9 and
  D3Q19 keep their existing code paths. The D3Q27 outflow and convective
  closures remain unimplemented rather than inheriting a new model silently.

### 2026-07-07 Guo-force-consistent Zou-He closure (`kernels.rs::zou_he_face_selected`)
- Form: velocity and pressure Zou-He faces reconstruct unknown populations
  from the raw boundary velocity
  `v = u_phys - F/(2 rho)`, where `F = F0 + rho g` is the same composed Guo
  force used by collision and `moments_row`. For velocity faces,
  `rho (1 - v_n) = S0 + 2 S- + 1` gives
  `rho = (S0 + 2 S- + 1 - F0_n/2) / (1 - u_n + g_n/2)`;
  reconstruction then uses `rho v_n` and `rho v_t`. For pressure faces,
  `rho` is prescribed, `v_n = 1 - (S0 + 2 S- + 1)/rho`, and zero physical
  tangential velocity is enforced by `v_t = -F_t/(2 rho)`. This is the
  pre-force-population form equivalent to adding the Guo source correction
  to the NEBB closure, but it preserves the existing no-force arithmetic
  exactly.
- Source: resolved from the engine's Guo moment convention
  `rho u = sum_q c_q f_q + F/2` and the existing Zou-He/Hecht-Harting
  moment closure. No empirical coefficient is introduced.
- Validity domain: planar axis-aligned D2Q9, D3Q19, and D3Q27 velocity and
  pressure faces, including masked face patches, under the existing
  one-open-axis and low-Mach prescribed-speed guards. The correction uses the
  backend-composed single-phase force `F0 + rho g`; other closures must feed
  the same Guo force path before this boundary pass.
- Validation: `crates/lbm-core/tests/zou_he_force.rs` passes for D2Q9,
  D3Q19, and D3Q27 with uniform force plus per-mass gravity
  (max velocity-face error `1.388e-17` or lower). Kernel unit
  `zou_he_d2q9_zero_force_matches_legacy_formula_bitwise` pins the zero-force
  D2Q9 branch against the legacy formula. The interaction matrix
  `feature_interaction_conservation_matrix` flips the uniform-force ×
  face-patch and gravity × face-patch cells green while the other feature
  pairs remain green/skip as documented.
- Replaces / interacts with: fixes ANOM-P4-021. Whole-face Zou-He and T18.2
  face patches share `zou_he_face_selected`, so both paths receive the same
  correction. The generated GPU `bc` kernel mirrors the same raw-velocity
  closure; the F32 collision byte-identity scope is unchanged.

### 2026-07-07 behavior review — D3Q27 open rectangular duct
Pattern: the duct run is predominantly unidirectional in `x`; transverse
velocity is negligible (`cross_rel = 3.588e-7`), mass flux is balanced
inlet-to-outlet within `4.212e-5`, and the plane-averaged density decreases
from inlet to outlet.
Mechanism: the prescribed inlet flux plus fixed outlet density establish an
axial pressure gradient; half-way wall rims impose the rectangular-duct shear
profile, and the D3Q27 closure supplies only the missing incoming face
populations needed to satisfy the boundary moments.
Resolved vs closure: the duct shear and pressure-gradient response are
resolved LBM dynamics with half-way walls; the only non-resolved term active
is the D3Q27 open-face moment closure recorded above.
Artifacts checked: `crates/lbm-core/target/d3q27_open_duct_profile.csv`.
No clamps, caps, walls-as-absorbers beyond the documented wall rim, or
partition seam artifacts were observed; the T13 open-duct split case was
bit-exact with seams crossing the boundary cells.
Verdict: PHYSICAL.
Routing: none.

Without these notes, someone re-tightening the bands from "1e-12 looks
tight" chases physically impossible numbers.

### 2026-07-07 — ANOM-P4-008: D3Q19 central-moment offset removed
- Form removed: the D3Q19 central-moment path previously adjusted the
  configured second-order shear relaxation as
  `omega_eff = omega_shear * (1 + 0.0025 - 0.16 |u|^2)`, clamped to `<= 2`.
  The `+0.0025` term was removed from CPU scalar/SIMD and generated GPU WGSL
  paths. The remaining velocity term is explicitly ablatable via
  `CENTRAL_MOMENT_DISABLE_VELOCITY_CORRECTION_FOR_ABLATION` and remains
  pending E1 verdict; normal builds keep it active so the Galilean-defect
  acceptance stays covered.
- Source: the removed offset was empirical calibration, not a first-principles
  closure. The decisive audit in
  `crates/lbm-core/tests/accuracy_audit_cumulant.rs` separates
  finite-resolution `O(h^2)` error from continuum bias with
  `nu_eff/nu - 1 = a + b/N^2`. With the offset removed, the light canary
  measured D3Q19 defects `(24, 4.009218693e-2)` and
  `(32, 2.256135694e-2)`, giving `a = 2.171836972e-5`, below the
  `|a| <= 4e-3` light band and consistent with the heavy acceptance target
  `|a| <= 2e-3`. The offset's own tau-space footprint matched the poisoned
  continuum intercept, so it was a banned calibration.
- Training / holdout split: training = TGV3D decay at the calibration settings
  used to select `+0.0025` and `0.16`; holdout =
  `crates/lbm-core/tests/cumulant_holdout.rs` with (1) advected TGV3D at mean
  frames `u_frame = {0, 0.05, 0.1}`, (2) off-calibration-Re TGV3D with
  `nu = 0.04`, and (3) D3Q19-vs-D3Q27 TGV3D at the off-calibration Re.
- Validation results from
  `cargo test -p lbm-core --release --test cumulant_holdout -- --include-ignored --nocapture`
  on 2026-07-07:
  - Advected TGV3D, D3Q19 cumulant, `N=32`, `nu=0.02`, `u0=0.012`,
    160 steps: rates were `4.634861882e-3` (`u_frame=0`),
    `4.639833752e-3` (`u_frame=0.05`), and `4.654339667e-3`
    (`u_frame=0.1`) vs analytic `4.626377063e-3`; relative errors were
    `1.834009427e-3`, `2.908688268e-3`, and `6.044168839e-3`. The
    frame-spread was `4.195075506e-3`, exceeding the derived
    `Ma_frame,max^2 * (k dx)^2 = 1.156594266e-3` band. **Finding:** the
    calibrated correction does not establish Galilean-invariant viscous decay
    on this holdout.
  - Off-calibration Re, D3Q19 cumulant, `N=32`, `nu=0.04`, `u0=0.012`,
    160 steps: rate `9.223420342e-3` vs analytic `9.252754126e-3`,
    `nu_eff = 3.987318896e-2`, relative error `3.170276018e-3`, passing the
    T15 decay-rate class band `2e-2`.
  - D3Q19 vs D3Q27 cross-check at the off-calibration Re: D3Q19 relative
    error `3.170276018e-3`; D3Q27 rate `9.299938977e-3`, analytic
    `9.252754126e-3`, relative error `5.099546655e-3`. D3Q19 is not an
    outlier against the D3Q27 error plus the T15 band (`2.509954666e-2`).
- Corrected acceptance: the viscosity gate is resolution-aware h² intercept,
  not a single N=32 value. The finite-N smoke in
  `crates/lbm-core/tests/cumulant_acceptance.rs` is re-frozen to the
  uncorrected D3Q19 N=32 measurement
  `nu_eff = 2.0454550535750255e-2`, relative error
  `2.2727526787512733e-2`, with a narrow `2.4e-2` band. D3Q27 remains under
  the existing `2.0e-2` smoke band with `nu_eff = 2.0187636744471944e-2`,
  relative error `9.381837223597193e-3`.
- Validity domain: no D3Q19 lattice-offset closure is live. The remaining
  `-0.16 |u|^2` central-moment velocity term is not validated as a viscosity
  correction; E1 remains SPEC-GAP until rerun with the ablation flag. It
  still has a measured Galilean-defect effect in the current acceptance:
  D3Q19 BGK `2.570488585e-3` vs CentralMoment `9.996091795e-4`, D3Q27 BGK
  `2.533062587e-3` vs CentralMoment `1.161446988e-3`.
- Replaces / interacts with: removes the previous D3Q19 empirical
  viscosity-offset calibration from the central-moment collision. It does
  not replace resolved LBM viscosity (`tau = 3 nu + 0.5`) and should not be
  reused as a generic LES, wall, or stability limiter coefficient.

### 2026-07-07 behavior review — cumulant holdout integral runs
Pattern: all reported TGV3D runs had positive decay rates and monotonically
decreasing fluctuation kinetic energy; the advected-frame sequence showed a
systematic increase in measured decay rate as `u_frame` increased.
Mechanism: the Fourier-mode velocity field decays by viscous diffusion, while
the frame trend indicates residual frame-dependent numerical viscosity after
the empirical D3Q19 correction.
Resolved vs closure: viscous decay and periodic streaming are resolved by the
LBM update. The reviewed historical run had both the now-removed D3Q19-only
`+0.0025` offset and the pending `-0.16 |u|^2` shear-rate adjustment active;
new runs after ANOM-P4-008 retain only the ablatable velocity term.
Artifacts checked: these are fully periodic integral-metric tests with no
walls, outlets, clamps, or seams. Per REV-6, no field-visualization artifact
was produced or expected for these integral metrics; the behavior anchor is
the positive-rate and monotone-energy check embedded in every run.
Verdict: CLOSURE-DRIVEN finding for finite-frame Galilean invariance; passing
off-Re and D3Q27 cross-checks do not validate the failed frame-dependence
behavior.
Routing: core cumulant follow-up / claim narrowing; do not loosen the
holdout band without a new physics derivation recorded here.

### 2026-07-07 WALE `tau_eff` upper clipping (`crates/lbm-core/src/les.rs:WaleLes`)
- Form: when explicitly configured, the WALE driver limits the effective
  symmetric relaxation time to `tau_eff <= tau_eff_max`, where
  `tau_eff = 1/2 + (nu_0 + nu_t)/(c_s^2 Delta t) = 1/2 + 3(nu_0 + nu_t)`
  in lattice units. The equivalent eddy-viscosity cap is
  `nu_t <= (tau_eff_max - (1/2 + 3nu_0)) / 3`. `None` leaves the raw WALE
  `nu_t` bit-identical to the unclipped implementation.
- Source: LBM BGK/TRT relaxation uses `omega_plus = 1/tau_eff`; very large
  `tau_eff` drives `omega_plus` toward zero and over-diffuses the resolved
  field. FR-LES-02 defines the effective relaxation relation and FR-LES-03
  requires upper clipping with diagnostics. The limiter is therefore a
  configured numerical-stability bound on the collision relaxation, not a
  calibration term for any validation band.
- Validity domain: applies only to WALE SGS viscosity in lattice units after
  the Nicoud-Ducros WALE closure has computed raw `nu_t`; the configured
  `tau_eff_max` must be finite, greater than `1/2`, and at least the laminar
  `tau_0 = 1/2 + 3nu_0` for the solver. Default is off; no silent physical
  default is installed.
- Validation: `crates/lbm-core/tests/wale_les.rs::wale_unset_clipping_matches_raw_wale_bitwise_on_sheared_field`,
  `::wale_tau_eff_clipping_diagnostics_match_reference`, and
  `::wale_tau_eff_clipping_count_is_monotone_with_bound`.
- Replaces / interacts with: augments the WALE SGS relaxation-field driver.
  Diagnostics are mandatory per update: clipped-cell count, clipped-cell
  fraction, maximum raw `nu_t` before clipping, configured `tau_eff_max`, and
  the equivalent active `nu_t` bound.

### 2026-07-07 behavior review — REV-4 WALE `tau_eff` clipping tests
Pattern: with clipping unset, the sheared multimode field preserves every raw
WALE `nu_t` bit; with clipping engaged, only cells whose raw `nu_t` exceeds
the active bound are clipped to that bound.
Mechanism: the pattern follows directly from limiting
`tau_eff = tau_0 + 3nu_t`, so clipping only lowers the SGS contribution in
cells that would otherwise exceed the explicitly configured relaxation-time
ceiling.
Resolved vs closure: velocity gradients and collision relaxation remain the
resolved LBM path; WALE is the active SGS closure; the new limiter is a
diagnosed numerical-stability limiter on that closure output.
Artifacts checked: no field-visualization artifact was produced because this
was a code+unit-test order, not an experiment/demo run; the behavior anchor is
the exact cell-by-cell clipped-set assertion.

### 2026-07-07 Backend-side gravity body-force composition (`StepParams::gravity`)
- Form: `F_total(x,t) = F_uniform + F_cell(x,t) + rho(x,t) * g`, composed in
  the backend collide/moment force path and entering the existing Guo source
  term. For a caller-owned per-cell force field, arithmetic grouping preserves
  the old staged overlay: `F_uniform + (F_cell + rho*g)`.
- Source: resolved from the existing per-mass gravity model and Guo forcing
  invariant recorded above; no new physics term, closure, or constant.
- Validity domain: same as the existing single-phase gravity path and Guo
  forcing path. Solid cells are skipped through the existing collide/moment
  masks. Future VOF/AGG work must replace only the density factor at this
  composition point, not the forcing scheme.
- Validation: `cargo test -p lbm-core gravity --release` passed
  `gravity.rs::closed_box_gravity_forms_stable_hydrostatic_stratification`,
  `gravity.rs::gravity_channel_is_bit_identical_to_raw_rho_g_force_field`,
  `gravity.rs::vr_str_06_static_stratification_quiescent_all_lattices_and_precisions`,
  and `gravity.rs::shan_chen_gravity_composes_with_additive_force_field_and_creates_buoyancy`.
  `cargo test -p lbm-core --test backend_simd_equiv --release`,
  `cargo test -p lbm-core --test t13_split_invariance --release`,
  `cargo test -p lbm-core --test t13_adversarial --release`, and
  `cargo test --workspace --release` passed. GPU build gate
  `cargo build -p lbm-core --release --features gpu` passed. Runtime GPU

### 2026-07-07 Shan-Chen additive force-field composition (`compat::multiphase`)
- Form: `F_cell(x,t) <- F_cell(x,t) + F_SC(x,t)` for SCMP and
  `F_cell_sigma(x,t) <- F_cell_sigma(x,t) + F_MCMP_sigma(x,t)` for MCMP.
  Callers clear/reset the caller-owned per-cell field once per step before
  composing transient sources; all composed sources enter the existing Guo
  source point together.
- Source: resolved from the existing Shan-Chen interaction force and Guo
  forcing composition invariant. This changes staging from overwrite to
  addition; it adds no new physical force term, closure, branch, or constant.
- Validity domain: same as the existing SCMP and MCMP Shan-Chen closures and
  the existing per-cell Guo force-field path. Repeated calls without a
  per-step clear intentionally accumulate source terms and are caller error,
  matching the rotor source contract.
- Validation:
  `validation_multiphase.rs::t11_shan_chen_adds_to_existing_force_field_anom_p4_022`
  asserts cell-by-cell `gravity + SC` composition and checks the composed field
  is neither contribution alone.
- Replaces / interacts with: replaces the old `copy_from_slice` staging in
  `ShanChen::update_force` and `MultiComponent::update_forces`; composes with
  rotor, raw per-cell sources, gravity, and uniform force through the existing
  Guo path.
  equivalence is covered by ignored T14 test
  `t14_backend_equiv::t14_gravity_body_force_device_resident` and is
  BENCH-PENDING on a native GPU adapter.
- Replaces / interacts with: replaces the host staging overlay used by
  `Solver::stage_gravity` on capable backends. The staged overlay remains as
  fallback for backends that do not advertise backend-side gravity. Shan-Chen,
  IBM, uniform force, and explicit per-cell force fields still combine through
  the same Guo force path.

### 2026-07-07 behavior review — backend-side gravity composition tests
Pattern: hydrostatic closed boxes stayed quiescent within the existing
machine-level/static-stratification bands; dense phase moved in the `-g`
direction and light phase rose in the Shan-Chen buoyancy sign test.
Mechanism: the backend composes the same `rho*g` body-force density into the
same Guo source term, so pressure/gravity balance and buoyancy signs are
unchanged while avoiding per-step host staging.
Resolved vs closure: gravity and Guo forcing are resolved engine terms; the
Shan-Chen cohesion used by the buoyancy sign test is the existing documented
multiphase closure and was not changed here.
Artifacts checked: no clamp, wall, seam, or outlet accumulation artifact was
introduced; T13 split invariance stayed green and the bit-identical
gravity-vs-raw-force-field test stayed green. No field-visualization artifact
is expected for this code-and-test order; the evidence is scalar assertions
from the named validation tests.
Verdict: PHYSICAL.
Routing: none.

---

## 3. Prohibited patterns

Ad-hoc physics is BANNED (owner directive, CLAUDE.md §Working discipline).
Every physical behavior must be resolved from the governing equations or a
literature-backed closure with a recorded derivation, validity domain, and
its own validation test.

Explicit prohibitions:
- Constants calibrated to pass a specific acceptance band.
- Branches keyed to sample or case identity ("harshness" switches).
- Position clamps or caps that silently absorb transport.
- Decorative physics terms with no derivation.

If a gate cannot be met without such a hack, **STOP and report** — the
spec gets revised, not the physics faked.

Executable procedure — grep-able ban-list smells, provenance decision
table, two-layer gate template, stop-rule report format, escalation table
— lives in `.claude/skills/lbmflow-physics-discipline/SKILL.md`. Every
developer agent follows it mechanically; every physics-affecting codex
order embeds its clauses. Do not duplicate that content here.

---

### 2026-07-07 Schiller-Naumann particle drag validity domain (`crates/lbm-core/src/particles.rs:particle_velocity`)
- Form: one-way particle drag uses
  `f(Re_p) = 1 + 0.15 Re_p^0.687`,
  `tau_p = rho_p d^2 / (18 rho_f nu f(Re_p))`, and the semi-implicit update
  `v_{n+1,a} = (tau_p v_{n,a} + u_a + tau_p g_a(1 - rho_f/rho_p)) / (tau_p + 1)`.
  The implemented validity boundary is `SCHILLER_NAUMANN_RE_MAX = 800.0`;
  `Re_p > 800` returns `ParticleError` with the particle index and offending
  Reynolds number.
- Source: Schiller and Naumann (1935), standard drag correction for isolated
  spherical particles; retained here as the existing T18.3 particle closure.
- Validity domain: one-way dilute spherical particles in the Schiller-Naumann
  drag-correction range, enforced as `0 <= Re_p <= 800` in this code path.
- Validation: `crates/lbm-core/src/particles.rs::tests::schiller_naumann_in_domain_matches_formula_and_is_monotone`
  checks bit-identical in-domain formula evaluation and monotonicity on
  `[0, 800]`; `crates/lbm-core/src/particles.rs::tests::schiller_naumann_out_of_domain_reports_particle_index_and_re`
  checks the hard error for `Re_p > 800`; existing T18.3 settling/deposition
  tests continue to cover in-domain particle transport.
- Replaces / interacts with: replaces the previous silent release-build
  `Re_p.min(800)` clipping in particle drag evaluation. The clip was a banned
  transport-absorbing behavior because it silently switched the drag law
  outside the closure's validity domain instead of making the invalid state
  visible to the caller.

---

## T17/VR-STR-03 — Re_tau=178.12 turbulent channel vs MKM DNS (first measurement, frozen 2026-07-07)

Setup: minimal-flow-unit body-force channel (delta=48: 128x98x72, Lx+=475,
Lz+=267, u_tau=0.008, TRT magic 3/16 + WALE), deterministic multimode init,
20 Te warmup + 30 Te statistics. Primary characterization runs on the GPU
(f32, on-device WALE; 280 s wall vs ~5.5 h CPU); the f64 CPU variant is the
reference-precision cross-check. GPU f32 reproduces the CPU harness smoke
values to 4 significant digits.

Measured (delta=48, equilibrium verified):
- Sustained turbulence: -<u'v'>+ at y+~30 over the last 10 Te = 0.729
  (MKM DNS ~0.75); peak -<u'v'>+ = 0.726. The resolved Reynolds shear stress
  MATCHES DNS.
- Total-stress force balance: nu dU/dy - <u'v'> vs u_tau^2(1 - y/delta),
  L2rel = 0.0535 over y/delta in [0.2, 0.8] — statistical equilibrium holds;
  the residual is the (unaccounted) SGS mean contribution + finite window.
- Mean profile: U+ L2rel vs MKM = 0.2328 over y+ in [5, 150]; centerline U+
  = 22.0 vs DNS 18.30 (+20%).

Behavior-validity review: the error pattern (correct Reynolds stress, correct
force balance, over-predicted mean gradient in the buffer/log region) is the
documented signature of wall-UNRESOLVED LES without a wall model at
y+/cell = 3.7. It is a resolution grade, not a model defect: WALE's null
gates (laminar shear, TGV small-strain) remain exact, and the equilibrium
diagnostics above rule out forcing/BC artifacts. Bands frozen accordingly:
mean U+ L2rel <= 0.30 (coarse-LES grade, measured 0.2328), stress balance
<= 0.10 (measured 0.0535), turbulence guard -<u'v'>+ > 0.4 (measured 0.729).
Improving the mean-profile grade needs either finer near-wall resolution
(delta >= 96, y+ <= 1.9 — GPU cost ~30 min, planned as a follow-up
characterization) or a wall model (roadmap item, not implemented).

Resolution ladder COMPLETE (2026-07-07, upsample-restart protocol for
delta>=80): delta=96 (256x194x144, y+/cell = 1.9, stage-A delta=48
equilibrium field upsampled trilinearly, 10 Te warmup, 30 Te statistics,
EXIT:0): mean U+ L2rel = 0.0549, centerline U+ = 19.08 (DNS 18.30, +4.3%),
total-stress balance L2rel 0.0674, -<u'v'>+ (y+~30, last 10 Te) = 0.740
(DNS ~0.75), peak 0.803. The ladder 0.233 (y+3.7) -> 0.155 (y+2.8) ->
0.055 (y+1.9) converges cleanly to the DNS profile: at near-wall-resolved
spacing the WALE channel reproduces MKM within ~5% mean-profile error with
DNS-grade Reynolds shear stress. Behavior review: all three equilibrium
diagnostics (stress linearity, sustained -u'v'+, achieved Re_tau) hold at
every rung; no closure or boundary artifact remains visible at delta=96.
The frozen delta=48 routine gate stands (coarse-grade band 0.30); the
ladder is the accuracy statement.

Resolution trend (delta=64, 171x130x96, y+/cell = 2.8, 50 Te warmup,
equilibrated — measured 2026-07-07, 1341 s GPU): mean U+ L2rel 0.1549
(from 0.2328), centerline U+ 20.58 (from 22.02; DNS 18.30), total-stress
L2rel 0.0192, -<u'v'>+ (y+~30, last 10 Te) 0.755 (DNS ~0.75), peak 0.733.
Every metric converges toward DNS with near-wall resolution, confirming the
resolution-grade classification. The frozen delta=48 bands stand; a
delta>=96 point and/or a wall model are the paths to a tighter mean band.

Harness notes: (i) the total-stress fold initially applied the half-channel
sign only to the viscous term — the Reynolds term needs the same fold
(upper-half error ~2x line; fixed before freezing); (ii) delta=64 at 20 Te
warmup produced a non-equilibrated stats window (peak -<u'v'>+ 0.975 above
the equilibrium ceiling, stress residual 34%) — equilibration is delta-
dependent; 50 Te equilibrates it (measurements above).

Fine-grid restart note (2026-07-07): PM-measured delta=80/96 GPU runs stayed
laminar from t=0 under both the original 3-mode seed and the deterministic
1/|k| multiscale ladder, then accelerated under the constant body force until
the Mach ceiling was exceeded. This is a transition-tripping failure on clean
fine grids, not a new turbulence model. The GPU validation path now uses the
standard channel-DNS remedy for delta>=80: equilibrate the frozen delta=48
case first, gather rho/u, trilinearly interpolate the turbulent field onto the
target grid with periodic x/z wrapping and one-sided wall-normal interpolation,
then run a short 10-Te relaxation warmup before collecting statistics.
Velocity is carried over unscaled because u_tau is held fixed by construction;
wall-rim velocities remain zero. The sustained-turbulence guard and the
force-balance/mean-profile bands remain the seed-independence checks for the
statistics window.

---

### 2026-07-07 WALE clipping validation-claim disclosure follow-up

WALE `tau_eff` clipping remains a diagnosed numerical-stability guard, not a
turbulence-model calibration knob. Any validation claim made with clipping
active must disclose that clipping was active and report the corresponding
`clipped_fraction` and `max_nu_t_before_clipping` diagnostics from
`WaleLesDiagnostics`.

---

### 2026-07-07 Guo force-field ingestion timing (`crates/lbm-core/src/solver.rs:refresh_moments_after_force_change`)
- Form: after an exactly uniform per-cell force field is installed or cleared
  (and after gravity is set on all-fluid domains), stored moments are refreshed
  with the same physical velocity definition used by uniform forcing:
  `u = (m + F/2) / rho`, where `F` is the composed Guo body force.
- Source: resolved Guo forcing contract already used by `moments_row` and the
  uniform-force initialization path; this removes ANOM-P2-001's inconsistent
  late force-field staging rather than adding a new closure.
- Validity domain: BGK/TRT Guo forcing paths in this engine for equivalent
  uniform force-field ingestion. Nonuniform model force fields and wall-rim
  gravity/hydrostatic cases keep their validated pre-existing staging.
- Validation: `crates/lbm-core/tests/accuracy_audit.rs::uniform_force_impulse_matches_force_field_anom_p2_001`
  asserts uniform force, equivalent per-cell force field, and equivalent
  gravity have identical step-1 momentum impulse to `1e-14`; T13 split
  invariance and `backend_simd_equiv` were rerun unchanged.
- Replaces / interacts with: replaces the previous equivalent uniform
  force-field path that collided step 1 from stale moments and lost the
  transient `1/(2 tau_minus) * F` contribution while steady slopes remained
  correct.

### 2026-07-07 D3Q27 open-face completion (`kernels.rs::{outflow_face_selected,convective_face_selected}`, `gpu/wgsl.rs`)
- Form: D3Q27 `Outflow` uses the existing zero-gradient open-face closure
  `f_q(x_face) = f_q(x_face + n)` for every incoming link
  `q in {c_q·n = 1}`. For D3Q27 this set has nine links, including the four
  body-diagonal corner links; no special corner coefficient is introduced.
  D3Q27 `Convective` uses the existing radiation update
  `f_q^{new}(face) = (f_q^{prev}(face) + U_c f_q^{new}(interior)) / (1 + U_c)`
  on the same incoming set, followed by the existing mass-consistency pin
  `rho(face) := rho(interior)` with the density deficit distributed as
  `Delta f_q = Delta rho * w_q / sum_{incoming} w_q`. The WGSL path now emits
  the D3Q27 velocity/pressure NEBB branch from the lattice tables and applies
  the same generic nine-slot outflow/convective loops.
- Source: zero-gradient extrapolation and the convective/radiation outlet are
  the same closures documented for the D2Q9/D3Q19 paths in this file,
  including the 2026-07-05 convective mass-consistency pinning entry. The
  D3Q27 change is only the direction-set extension: `L::unknowns(face)` is the
  analytic incoming plane for the lattice, and the formulas are per-link
  scalar extrapolation/advection equations, so corner links use the identical
  equation as axial and face-diagonal links. Velocity/pressure GPU NEBB uses
  the D3Q27 derivation recorded in the 2026-07-07 D3Q27 velocity/pressure
  entry above.
- Validity domain: planar axis-aligned D3Q27 open faces under the existing
  one-open-axis rule, low-Mach velocity guard, positive outlet density guard,
  and `0 < U_c <= 1` convective-speed guard. Zero-gradient and convective
  outlets are extrapolating/radiation closures, not pressure-prescribing
  outlets; wall-bounded duct tests therefore pin outlet-local distortion
  against the inherited D3Q19-equivalent envelope rather than claiming the
  pressure-outlet mass-flux behavior.
- Validation: `crates/lbm-core/tests/d3q27_open_bc.rs`:
  `d3q27_outflow_duct_matches_profile_and_d3q19` measured D3Q27-vs-D3Q19
  L2rel `3.504e-4`, flux-scaled profile L2rel `9.844e-4`, outlet-local
  flux envelope `2.301e-1`, cross-flow ratio `4.203e-3`, and monotone
  pressure drop. `d3q27_convective_duct_matches_profile_and_d3q19` measured
  D3Q27-vs-D3Q19 L2rel `4.152e-4`, flux-scaled profile L2rel `2.279e-3`,
  outlet-local flux envelope `1.982`, cross-flow ratio `3.630e-2`, and
  monotone pressure drop. `d3q27_open_outlets_balance_uniform_through_flow_mass_flux`
  measured uniform-flow mass-flux imbalance `0.000e0` for Outflow and
  `1.735e-16` for Convective. `t13_d3q27_open_duct_split_invariant_with_bc_seams`
  now covers pressure, outflow, and convective outlets with seams crossing BC
  cells. `crates/lbm-core/tests/d3q27_open_metamorphic.rs` converts the old
  rejection pins to all-kind CPU acceptance and adds a T14-style D3Q27 GPU
  CPU-vs-GPU duct equivalence test; the sandbox run compiled the GPU feature
  path and skipped execution because no adapter was available
  (PENDING-NATIVE-RUN).
- Replaces / interacts with: removes the CPU `UnsupportedOpenFaceKind` guard
  for D3Q27 `Outflow` and `Convective`, and removes the GPU
  `UnsupportedOnGpu { feature: "D3Q27 open faces" }` guard for implemented
  D3Q27 open faces. Unimplemented GPU localized features remain rejected.

### 2026-07-07 behavior review — D3Q27 outflow/convective ducts
Pattern: pressure-outlet D3Q27 remains a clean unidirectional rectangular
duct; D3Q27 Outflow and Convective match D3Q19 closely but show larger
outlet-local flux and transverse-velocity distortion in the wall-bounded duct.
Mechanism: zero-gradient and radiation outlets extrapolate incoming
populations rather than prescribing pressure, so the immediate outlet adjusts
to evacuate the imposed inlet profile while the D3Q27 corner links follow the
same per-link extrapolation/advection law as D3Q19's incoming set.
Resolved vs closure: the duct core, wall friction, and NEBB inlet are resolved
LBM behavior; the outlet-local distortion is closure-driven by the
zero-gradient/radiation outlet model and the convective mass pin. Uniform
through-flow confirms the closures balance mass when their extrapolation
assumption is exactly satisfied.
Artifacts checked: no clamps or case-identity branches were introduced;
T13 seams crossing pressure/outflow/convective BC cells are bit-exact; the
wall-bounded outlet artifact is pinned separately from the upstream profile
and D3Q19 consistency metrics. Visual artifact for PM review:
`crates/lbm-core/target/d3q27_pressure_duct_profile.csv`,
`crates/lbm-core/target/d3q27_outflow_duct_profile.csv`, and
`crates/lbm-core/target/d3q27_convective_duct_profile.csv`.
Verdict: CLOSURE-DRIVEN (validated against inherited D3Q19-equivalent behavior
and exact uniform-flow mass balance).
Routing: none; native GPU execution of the new T14-style D3Q27 open-face
equivalence test remains PENDING-NATIVE-RUN on a host with a GPU adapter.

### 2026-07-07 W-VOF O1 conservative Allen-Cahn phase-field transport (`crates/lbm-core/src/phase_field.rs`)
- Form: D3Q19 phase-field distribution `g_i` in ordinary q-major padded SoA
  form, `phi = sum_i g_i`, prescribed-velocity equilibrium
  `g_i^eq = phi w_i [1 + 3 c_i.u + 4.5(c_i.u)^2 - 1.5 u.u]`,
  `tau_phi = 3M + 0.5`, and collision-stream transport recovering
  `d_t phi + div(phi u) = div(M[grad(phi) - (4/W)phi(1-phi)n])`.
  The shared interface-flux helper is
  `J_phi = -M[grad(phi) - (4/W)phi(1-phi)n]`, with D3Q19 isotropic stencils
  `grad(phi) = 3 sum_{i>0} w_i c_i phi(x+c_i)` and
  `lap(phi) = 6 sum_{i>0} w_i [phi(x+c_i)-phi(x)]`.
  O1 implements transport only: no density feedback, viscosity feedback,
  surface-tension force, gravity edit, wetting boundary, sparger inlet, or
  hydrodynamic momentum `J_rho` coupling is active.
- Source: Chiu and Lin 2011 conservative Allen-Cahn counter-term, adopted via
  the Fakhari, Mitchell, Leonardi and Bolster 2017 velocity-based LBE form
  frozen in `docs/proposals/WVOF_IMPL_SPEC.md`. The source first moment is
  mobility-scaled so the relaxation-provided diffusive flux and the
  counter-flux share the same `M`, as required by the governing equation.
  This resolves the spec Eq. (4)/(6) prefactor ambiguity by enforcing Eq. (1)
  rather than introducing an extra fitted coefficient.
- Validity domain: O1 enforces `W in [4,5]`, `M in (0,1/6]`
  (`tau_phi in (0.5,1.0]`). Validation runs here used `W=4`, `M=0.04`,
  `Ma_lattice=|u|=0.055` for diagonal advection, periodic D3Q19 domains,
  and prescribed velocity only. Resolved droplet test uses `d/W=4`, the
  W-VOF lower bound; smaller bubbles remain a point-bubble follow-up, not O1.
- Validation: `crates/lbm-core/tests/wvof_o1_phase_field.rs`:
  `flat_interface_at_rest_holds_tanh_profile` passed the stated
  second-order `W=4` profile band (`L2_rel < 0.08`);
  `diagonal_periodic_droplet_advection_conserves_mass_and_profile` measured
  mass drift `1.9475e-14` over 1000 steps (band `<0.1%/1000 steps`),
  returned-profile `L2_rel = 1.2477e-1` (band `<0.14`), and interface width
  `7.193` (band `[0.5W,2W] = [2,8]`). `g_none_keeps_existing_hydrodynamic_path_bit_identical`
  proved the disabled path leaves D3Q19 hydrodynamic fields and all `f`
  planes bit-identical. Artifact:
  `target/wvof_o1/droplet_profile_before_after.csv`.
- Replaces / interacts with: adds the reserved optional `g`, `gtmp`, and
  compact `phi` fields to `SoaFields`; all are `None` by default. Existing
  hydrodynamic `f` pass order remains unchanged. Future O2 surface tension,
  density/viscosity feedback, gravity composition, and `J_rho` momentum
  correction must consume this same `J_phi` path.

### 2026-07-07 behavior review — W-VOF O1 advected droplet
Pattern: a diffuse spherical `phi=1` droplet translated diagonally by one
period in a fully periodic D3Q19 box under uniform prescribed velocity and
returned to its initial profile within the stated O1 profile band.
Mechanism: periodic pull-streaming advects the `g` distribution while the
conservative Allen-Cahn counter-flux balances relaxation diffusion at the
interface.
Resolved vs closure: the pattern comes from the resolved D3Q19 LBE transport
plus the literature-backed conservative Allen-Cahn counter-term; no
hydrodynamic coupling, surface-tension force, density feedback, or limiter was
active.
Artifacts checked: no position or `phi` clamp is present; periodic seams are
handled through the same ordered halo-layer exchange as `f`; exported line
profile before/after is `target/wvof_o1/droplet_profile_before_after.csv`.
Verdict: PHYSICAL for the O1 prescribed-velocity, matched-density transport
scope.

### 2026-07-07 wall-metrics observable for y+ / u_tau (`crates/lbm-core/src/wall_model.rs`, `Solver::gather_wall_metrics`)
- Form: read-only wall-adjacent diagnostics report `y_w`, tangential speed
  `u_parallel`, `u_tau`, `y+ = y_w u_tau / nu`, and `tau_w / rho = u_tau^2`
  in global compact cell order. Half-way rim walls use `y_w = 0.5 dx`;
  Bouzidi records use `y_w = qd |c_q|`. The turbulent branch solves
  `u_parallel / u_tau = ln(y_w u_tau / nu) / kappa + B` by Newton iteration
  with `kappa = 0.41`, `B = 5.2`; if the resulting `y+ < 11.6`, the observable
  reports the viscous branch `u_tau = sqrt(nu u_parallel / y_w)`.
- Source: half-way wall placement is the existing LBM wall convention in this
  codebase; Bouzidi `qd` is the geometric wall-link distance from
  Bouzidi-Firdaouss-Lallemand (2001). The log-law constants and branch switch
  are the standard smooth-wall equilibrium law used by the Malaspinas and
  Sagaut LBM wall-model class, frozen in
  `docs/proposals/LES_WALL_TREATMENT_SPEC.md` for W1/W2.
- Validity domain: W1 is instrumentation only. The log-law diagnostic is
  meaningful for attached smooth-wall equilibrium layers, with the future wall
  model's intended log-region domain `30 <= y+ <= 300`; sub-buffer cells are
  surfaced by the laminar branch and their reported `y+`, not hidden.
- Validation: `crates/lbm-core/tests/wall_metrics.rs` checks the laminar
  Poiseuille wall-shear band from the first-node one-sided discretization,
  quiescent `u_tau = 0`, no metrics in a fully periodic box, and Bouzidi
  distance-controlled y+ scaling for fixed `u_tau`. Focused gate:
  `cargo test -p lbm-core --test wall_metrics --release` passed on
  2026-07-07.
- Replaces / interacts with: this adds no closure and does not call
  `set_omega_field`; W2 may consume the observable, but W1 has zero effect on
  computed density, velocity, populations, or collision rates.

### 2026-07-07 behavior review — wall-metrics W1 code/test order
Pattern: no flow field pattern was generated or modified; the order adds a
read-only diagnostic over existing wall geometry and velocity moments.
Mechanism: reported wall metrics follow directly from geometric wall distance,
tangential velocity projection, and the specified wall-law/viscous formulas.
Resolved vs closure: half-way/Bouzidi wall distance is resolved geometry;
`kappa`, `B`, and the `y+ = 11.6` switch are diagnostic wall-law closure
constants, not active simulation terms in W1.
Artifacts checked: no clamps, outlets, seams, or wall-population writes are
introduced; no field visualization artifact is expected for this code+test
diagnostic-only order.
Verdict: PHYSICAL for diagnostics; no computed-field behavior claim is made.
Routing: none.

### 2026-07-07 falsification record — cumulant |u|^2 term removal experiment
- Experiment (branch cx/galilean-fix, not merged): removed the D3Q19
  `omega_eff = omega_shear * (1 - 0.16 |u|^2)` factor per
  docs/proposals/CUMULANT_GALILEAN_FIX.md path (c), which PREDICTED the
  advected-TGV3D frame spread would drop below the derived band 1.156594266e-3.
- Measured: spread WITH the term 4.195075506e-3 (prior holdout record);
  spread WITHOUT the term 1.051034711e-2 — 2.5x WORSE. The falsifier fired;
  the term was restored and no change was merged.
- Interpretation: the scalar |u|^2 factor does compensate a real |u|^2-scaling
  portion of the D3Q19 defect (presumably the isotropic/trace part); the
  proposal's diagnosis that removal would land sub-band is falsified. The
  residual anisotropic (tensorial) part remains uncorrected either way.
- Standing state: the calibrated term is RETAINED with the narrowed claim
  (valid for non-advected/weakly-framed decay; Galilean invariance at finite
  frame velocity NOT established on D3Q19 — holdout test stays
  `#[ignore = FINDING]` at its derived band). Any future fix must reproduce
  BOTH measured spreads (with-term 4.20e-3, without-term 1.05e-2 at
  u_frame <= 0.1, N=32, nu=0.02, u0=0.012) before its correction claim is
  trusted.

### 2026-07-07 addendum — Galilean-fix anchors invalidated by r2-c collision seam
- PM bisect: the merge d35faf4 (r2-c scalar TRT collision seam + ANOM-P2-001
  impulse fix) changed the D3Q19 cumulant TGV3D decay observable: advected
  u_frame=0 rel_err 1.834009427e-3 (at c272909) -> 2.506837906e-2 (at d35faf4
  and after). D3Q27 cross-check bit-unchanged -> D3Q19-specific. The default
  bands did not catch the shift (off-Re band 2e-2-class; band vacuity).
- Consequence: the falsification-record anchor numbers above and the
  CUMULANT_GALILEAN_FIX round-2 error model are valid for the PRE-r2-c tree
  only. The cubic-correction implementation order's STEP-0 anchor gate
  correctly detected the mismatch and halted without touching code.
- Status: cumulant Galilean work is BLOCKED on the r2-c triage (routed to the
  r2-c owner: intentional physics change to be recorded + holdout values
  re-frozen, or regression to be fixed). Do not rebuild the error model until
  the baseline is adjudicated.

### 2026-07-07 W-VOF O1 counter-term interface-maintenance fix (`crates/lbm-core/src/phase_field.rs`, `Solver::phase_field_step_prescribed_velocity`)
- Form: the D3Q19 conservative Allen-Cahn source remains the Fakhari
  velocity-form counter-term
  `S_i = w_i (c_i . [(4/W) phi(1-phi)n]) / cs^2`, but the discrete collision
  update now applies it as `(M/tau_phi) S_i = omega_phi M S_i`. Since the
  recovered source flux carries the BGK factor `tau_phi`, this gives the
  governing counter-flux `M(4/W) phi(1-phi)n` and balances the relaxation
  diffusion `M grad(phi)` at the tanh profile. The D3Q19 gradient is still
  `grad(phi) = 3 sum_i w_i c_i phi(x+c_i)`; its edge sums are evaluated in a
  permutation-invariant order so coordinate rotations do not seed round-off
  differences into the nonlinear normal.
- Source: Chiu and Lin 2011 conservative Allen-Cahn equation as discretized by
  Fakhari, Mitchell, Leonardi and Bolster 2017, PRE 96, 053301; equations and
  validity domain frozen in `docs/proposals/WVOF_IMPL_SPEC.md` section 1. The
  prior implementation multiplied `S_i` by `M` directly, leaving the recovered
  counter-flux short by `1/tau_phi` and causing relaxation diffusion to broaden
  the interface.
- Validity domain: unchanged W-VOF O1 domain, `W in [4,5]`,
  `M in (0,1/6]`, D3Q19 prescribed-velocity transport with no density,
  surface-tension, wetting, or hydrodynamic momentum coupling.
- Validation: `cargo test --release -p lbm-core --test wvof_o1_adversarial
  -- --nocapture` passed all six tests after un-ignoring the four FINDING
  gates. Measured gate values: width transit `u=0.02` final width
  `4.116276031`, phi drift `5.35e-13`; width transit `u=0.08` final width
  `4.103498626`, phi drift `2.00e-13`; two-droplet one-period
  `L2rel=1.311946682e-2`, phi drift `2.38e-13`; rotation metamorphic
  `Linf=8.881784197e-16`; counter-term sign anchor width
  `5.000000000 -> 4.089462844`; mobility edges `M=1e-4` width
  `4.177072435`, `M=1/6` width `4.105201912`.
- Replaces / interacts with: replaces the O1 phase-field source coefficient
  only. The public flux helper `J_phi = -M[grad(phi) -
  (4/W)phi(1-phi)n]` is unchanged; `g=None` hydrodynamic behavior remains
  bit-identical.

### 2026-07-07 behavior review — W-VOF O1 counter-term sign anchor
Pattern: a deliberately diffused `W=5` spherical interface at zero flow
monotonically sharpened toward the configured `W=4` profile, with the fitted
width decreasing from `5.000000000` to `4.089462844` over 500 phase-field
steps and total `phi` conserved at round-off in the adversarial suite.
Mechanism: with `u=0`, the conservative Allen-Cahn relaxation diffusion
`M grad(phi)` is opposed by the recovered counter-flux
`M(4/W)phi(1-phi)n`; because the initial interface is wider than the target,
the anti-diffusive term dominates until the tanh balance is approached.
Resolved vs closure: the pattern comes from the literature-backed CAC
counter-term and D3Q19 isotropic stencil; no limiter, position clamp, tuned
constant, surface-tension closure, density feedback, or hydrodynamic coupling
was active.
Artifacts checked: no boundary accumulation is present in the fully periodic
box; the rotation gate confirms the source stencil follows the repo's
orientation convention. Width-vs-time CSV:
`target/wvof_o1/counter_term_width_vs_time.csv`; midplane PGM:
`/tmp/lbmflow_wvof_o1_adversarial/counter_term_sharpening.pgm`.
Verdict: PHYSICAL for W-VOF O1 prescribed-velocity, matched-density transport.
Routing: none.

---

### 2026-07-07 direct-forcing IBM full-step impulse and overlap mobility (`crates/lbm-core/src/solver.rs:apply_rotating_ibm`)
- Form: each IBM sweep solves the marker-space correction
  `M q = U_marker - I[u]`, with marker impulse unknowns `q_k` spread as
  `q_k W_k(x)` and mobility
  `M_jk = sum_x W_j(x) W_k(x) / rho(x)`. The implemented Richardson sweep
  uses the row-sum preconditioner
  `G_j = sum_k M_jk = sum_x W_j(x) sum_k W_k(x) / rho(x)` and applies
  `q_j += relaxation * slip_j / G_j`; the cell predictor and force field both
  use the realized full Guo impulse `delta u = F / rho`.
- Source: resolved Guo forcing contract after the R2-C force-field impulse
  fix. The old IBM sizing targeted the half-force diagnostic increment
  `sum_x W F/(2 rho) = slip`, which made the realized post-step impulse twice
  the requested slip. For a single marker, `G_j = M_jj`, so the correction is
  the exact full-step direct-forcing impulse. For overlapping markers,
  `M` is the symmetric positive regularized-delta Gram matrix used in
  Uhlmann/Wang multi-direct-forcing analyses; `G^-1 M` has eigenvalues in
  `[0, 1]` under the row-sum bound, so relaxation `1.0` is non-amplifying for
  the represented marker modes instead of exciting the dense-marker collective
  gain.
- Validity domain: marker-based rotating IBM with finite positive density,
  kernel radius 1 or 2, and audited circular marker spacings `ds/h` in
  `[0.39, 1.0]`. The row-sum preconditioner treats marker-overlap stability;
  it is not a volume-penalization cure. Compat rotor penalization applies a
  per-cell correction directly in every solid cell with no interpolation-
  spreading Gram operator, so the coherent `chi=1` disc can still realize the
  explicit reflection `u -> 2 U_target - u` and remains ANOM-P4-010 until it
  gets its own implicit or derived relaxation treatment.
- Validation: `crates/lbm-core/tests/accuracy_audit_ibm.rs` default audit
  passed. Measured B1 torque ratio `1.06598` (relative error `6.598e-2`);
  B2 sub-cell torque spread `1.641e-3`; B3 near-edge conservation diagnostic
  `2.313e-5` on a near-zero net force; B4 slip refinement `5.802e-5 ->
  1.112e-5`; B5 kernel/relaxation torques `3.685e-2`, `3.738e-2`,
  `3.703e-2`; B6 `Omega -> -Omega` torque antisymmetry `1.041e-16` and
  mapped field difference `9.015e-17`; B8 Taylor-Couette profile
  `L2_rel=4.522e-2`, `Linf/U_i=3.542e-2`. `rotating_ibm.rs` now pins the
  corrected unit profiles with Taylor-Couette `L2_rel=2.887e-2`,
  `Linf/U_i=2.778e-2`, and torque within 10% of the annular-Couette value;
  the one-cell-off-wall Couette case is marked as smoke-only and superseded by
  the audit for quantitative IBM accuracy.
- Replaces / interacts with: replaces the old `2 * slip / M_jj` half-force
  sizing and immediate marker-order Gauss-Seidel update. Diagnostics now
  accumulate all sweep impulses and report slip after the applied sweep; net
  force conservation uses a `1e-12` relative-scale floor so a nearly zero net
  force is not reported as an arbitrary huge relative error.

### 2026-07-07 behavior review — ANOM-P4-001 IBM audit
Pattern: the corrected rotating IBM runs remain finite at the default
`relaxation=1.0`; torque is linear/stable, `Omega -> -Omega` produces the
mirrored velocity field to round-off, and the Taylor-Couette profile has the
expected annular monotone swirl between the rotating IBM circle and the outer
stationary wall.
Mechanism: Guo forcing supplies the full-step marker impulse, while the
row-sum marker mobility damps the collective overlap mode that previously
doubled and amplified the applied slip.
Resolved vs closure: the Guo impulse, TRT collision, and half-way outer wall
are resolved engine behavior; the IBM regularized-delta marker coupling is a
documented direct-forcing discretization with the validity domain above.
Artifacts checked: no clamps, ramps, or case-identity branches were added.
The B8 audit no longer places a stationary Eulerian solid core inside the IBM
kernel support, avoiding a separate half-way wall artifact in the profile
gate. No visual artifact was generated in this coding session; the PM/V&V
viewer pass should use the B8 velocity field if spatial inspection is needed.
Verdict: PHYSICAL within the marker-IBM validity domain; ANOM-P4-010 remains
separate for volume penalization.
Routing: none for ANOM-P4-001; ANOM-P4-010 stays routed separately.

---

### 2026-07-07 behavior review — ANOM-P4-016 MCMP i3 RT cutoff canary
Pattern: after the ANOM-P4-022 additive-force fix, the i3 canary still fails.
Mode 3 develops large downward bulk momentum (`p_total_y=-1.44e2` at step 10,
`-1.11e3` at step 100), then a lower-wall velocity/density failure
(`max|u|=4.013e3` at `heavy:(10,12)`, `rho_min=-6.009` at
`heavy:(160,10)` by step 400). Mode 7 remains finite through step 400 but
shows the same bulk-momentum signature and lower-wall high-speed locus.
Mechanism: the hard i3 protocol applies gravity to the heavy component only
inside a closed box, injecting a large nonzero net vertical body force; the
closed-wall return flow becomes wall-adjacent and violates the low-Mach/density
stability envelope before the intended mid-height RT cutoff behavior is
measurable.
Resolved vs closure: SC/MCMP interaction forces and Guo force ingestion are
the existing documented closures/source path; no new term was introduced. The
failure is driven by the test's per-component body-force protocol interacting
with closed walls, not by the ANOM-P4-022 staging bug.
Artifacts checked: lower-wall loci in diagnostic output; density maps
`target/vv_rt_i3/rt_mode3_step400_heavy.pgm`,
`target/vv_rt_i3/rt_mode3_step400_light.pgm`,
`target/vv_rt_i3/rt_mode7_step400_heavy.pgm`, and
`target/vv_rt_i3/rt_mode7_step400_light.pgm`.
Verdict: UNKNOWN/SPEC-COUPLING. A physical fix requires a derived MCMP
buoyancy forcing protocol (for example, a validated pressure-balanced or
zero-net-force formulation) and a validation update; no mean-force
subtraction or tuning was applied in this order.
Routing: PM/spec decision for ANOM-P4-016; core implementation only after the
forcing model is derived and accepted.

---

### 2026-07-07 checkerboard ghost-mode decay envelope (`crates/lbm-core/tests/accuracy_audit_modes.rs`)
- Decision: the `(pi, pi)` odd-even density-mode audit asserts monotone
  non-increase of a three-sample absolute-amplitude envelope, not strict
  adjacent-sample monotonicity of `|amp(pi,pi)|`.
- Rationale: the checkerboard perturbation is a Brillouin-corner ghost mode.
  It is strongly damped by the BGK/TRT collision spectrum, but the signed
  coefficient can oscillate after the mode reaches round-off and for TRT near
  `tau=0.5`. Adjacent absolute values may therefore have tiny rebounds even
  though the physical envelope is decaying. The load-bearing physics anchor is
  envelope damping plus no growth and no transfer into orthogonal staggered
  modes.
- Validation: `cargo test -p lbm-core --release --test accuracy_audit_modes -- --nocapture`
  passed on 2026-07-07. Measured BGK `tau=0.6`: raw monotonicity `0.300`,
  envelope monotonicity `1.000`, final/initial `1.942077e-12`,
  max growth `1.000000`, max leakage sum `1.409463e-18`. Measured TRT
  `tau=0.51`, `Lambda=3/16`: raw monotonicity `0.800`, envelope monotonicity
  `1.000`, final/initial `2.941800e-5`, max growth `1.000000`, max leakage
  sum `1.843144e-18`.
  The frozen gates are final/initial `<1e-8` for BGK and `<1e-3` for TRT,
  leakage `<1e-12`, and max growth `<2`.

### 2026-07-07 — Volume-penalization validity domain (ANOM-P4-010 / P4-027, V&V-ruled)

Volume penalization of a rotating body (`compat::rotor`, `F = ρ·(implicit
sizing)·(u_target − u*)`; the half-force-pinning `2ρχ` form was replaced with
the physical-implicit sizing whose one-cell recurrence has |gain| ≤ 1) has a
bounded validity domain, established by the rescoped rotor audit (F1–F6) and
V&V cross-referee against the direct-forcing IBM:

- **Domain = thin/porous structures AND short transients.** Coherent solid
  interiors are out of domain: the F6 witness shows a full disc grows
  detectably (crosses |u|=0.3 at step ~698 with the fixed sizing, vs ~96 with
  the old sizing) rather than holding rigid — distributed Darcy drag cannot
  enforce a rigid interior. Route coherent solids to rotating IBM (validated,
  slip ~9e-4) or Bouzidi.
- **Horizon caveat (ANOM-P4-027, S3).** Even a thin blade accumulates
  Darcy-drag dissipation error at long horizons: with the fixed sizing the
  thin-blade case is stable through ~6k steps (torque −6.19e-2, matching IBM
  −4.17e-2 to cross_rel 0.39) but diverges to NaN at ~26–28k steps at moderate
  resolution (was NaN at ~6k with the old sizing — the sizing fix is a strict
  improvement, not a cure). Validity window ≈ < 20k steps at this resolution;
  route long-time thin-blade studies to rotating IBM (stable indefinitely).
- The IBM–penalization cross-model torque difference (~0.39 at moderate
  resolution / before true steady state) is within the O(20–40%) distributed-
  drag-vs-sharp-interface range of the Angot–Bruneau–Fabrie penalization
  literature; not a defect.

Ruling (V&V concur 2026-07-07): the sizing fix landed as a strict improvement;
the audit F1–F3 are re-scoped as domain-boundary witnesses (#[ignore]'d,
runnable), F5/F6 green, F4 thin-blade referee documented-red pending an
acceptance criterion that states the short-horizon window. All rotor-audit
items are mf-interim-gated and do not affect default landing gates.
