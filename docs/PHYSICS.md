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
- **Collision — cumulant (cascaded central-moment)** for D3Q19/D3Q27,
  `CollisionKind::Cumulant { omega_shear }`. Implemented as a cascaded
  central-moment operator (not logarithmic cumulants); the code comment in
  `solver.rs` names this honestly. D3Q19 uses the D3Q27 tensor-product
  basis with the eight `x·y·z` corner moments dropped, matching the missing
  body-diagonal populations. Targets are the DISCRETE second-order Hermite
  equilibria used by BGK/TRT, not continuous Maxwellian moments (see §2
  "Cumulant target choice").
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
  faces via one formula; see doc comments in `params.rs` / `kernels.rs`.
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
  (`bouzidi.rs`), from Bouzidi-Firdaouss-Lallemand 2001.
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
uses those as relaxation targets, plus a D3Q19-only shear-rate offset
(+0.0025 rel) and a finite-frame cubic-velocity correction
`ω_eff = ω_shear · (1 + offset − 0.16 |u|²)` calibrated against TGV3D.
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

Without these notes, someone re-tightening the bands from "1e-12 looks
tight" chases physically impossible numbers.

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
