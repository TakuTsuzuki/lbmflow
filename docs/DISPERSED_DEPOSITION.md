# DISPERSED_DEPOSITION.md — Dispersed-Phase Deposition-Design Tool (D-track)

Status: P0/P1/P1.1 complete and verified · P2 (core promotion) in progress.
Acceptance criteria live in [VALIDATION.md](VALIDATION.md) **T18**.
Decisions in §3 are **frozen** — implement, do not re-litigate.

## 1. Goal

An AI-agent-native design tool for depositing a rigid dispersed phase
(suspended in a carrier fluid) onto an M×N target partition, driven by a
protocol over four operation primitives (withdraw / eject / agitate / settle):

- **Forward**: predict the deposited number-density field n(x,y) over the
  target partition, plus its statistics (CV, Max/Mean, empty_bin_fraction).
- **Inverse** (P4): find protocol parameters producing n(x,y) closest to a
  prescribed TARGET field n\*(x,y). Uniform (n\* = const) is only the special
  case; gradients, rings, and arbitrary patterns are first-class. The solver
  returns the closest reachable solution within the set spanned by the control
  degrees of freedom.

**Vocabulary policy**: ALL vocabulary stays domain-neutral (particle-laden
multiphase flow). No application-domain terms anywhere — code, docs, JSON
keys, metrics, commit messages.

## 2. Current state (verified 2026-07-06, PM-reproduced)

Isolated 3D (D3Q19) demo at `crates/lbm-cli/examples/dispersed_seeding/`
(no core changes; raw findings in its `SPEC_FINDINGS.md`). Runs
withdraw→eject→agitate→settle over a tall reservoir (visualization-only LBM)
plus a shallow target tray; deposits particles; emits density CSV +
`metrics.json` + qa-viewer-compatible STRUCTURED_POINTS VTK.

Verified numbers (low-Mach band): gentle Ma 0.093 / harsh Ma 0.099, τ = 0.536;
gentle empty_bin_fraction = 0.0; trend CV(gentle) = 1.16 < CV(harsh) = 4.15.
The gentle case is a PHYSICAL center-heavy radial gradient from a single jet —
uniformity requires the multi-jet + agitation + repeated-pass levers, which is
what motivates the inverse solver.

## 3. Frozen spec

### 3.1 Units and numerical guards

- SI units are authoritative; the grid is derived from `dx_m` (a mismatch of
  more than 1 cell between supplied grid counts and SI/dx_m is a validation
  error). The reservoir may keep a coarse visualization spacing (`grid.dx_m`)
  while the tray refines via `grid.tray_dx_m`.
- Abort (validation error, not warning) if Ma > 0.3 or τ < 0.51.
- Nondimensionalize via **diffusive scaling** (refine dx, dt ∝ dx² at fixed
  physical ν) to stay in the accurate low-Mach band (target Ma ≤ 0.1).

### 3.2 Operation primitives

- `settle(duration_s)` — gravity-only evolution.
- `withdraw(depth_frac, volume_frac, rate_uLs)` — `depth_frac` 0 = filled
  surface, 1 = floor; extraction is a statistical 1D settling-column
  concentration model.
- `eject(points_xy_frac[], rate_uLs, height_m, nozzle_diameter_m)` — each
  point is a jet disk; `u_jet = Q / (π (nozzle/2)²)`; `nozzle_diameter_m` is
  required; multi-point supported.
- `agitate(pattern, count, speed_mms, amplitude_mm)` — translational only;
  an unknown `pattern` is a validation error; Fr = Aω²/g.

### 3.3 Particle model

One-way Lagrangian: semi-implicit Stokes drag + gravity + agitation body
force. Deposit on floor crossing. Suspended particles are reported as
`n_suspended` and are NEVER projected onto the floor. `max_particle_steps`
aborts rather than silently changing the physics.

### 3.4 Metrics (`metrics.json`)

`CV`, `max_over_mean`, `empty_bin_fraction`, `n_extracted`, `n_deposited`,
`n_suspended`, `Re_jet`, `St`, `Fr`, `Ma`, `tau`. Every run logs the
non-dimensional REGIME line.

### 3.5 Frozen bands

- Gentle single-jet CV band: **1.05 ≤ CV ≤ 1.30** (frozen at P1.1 low-Mach
  measurement; the old 0.95–1.40 band was measured at Ma ≈ 0.25 in the
  compressibility-error regime and is retired).
- Gentle empty_bin_fraction target: ≤ 0.15 (measured 0.0).

## 4. Phasing

| Phase | Content | Status / gate |
|---|---|---|
| P0 | Isolated skeleton; trend gate; 10 spec findings filed | ✅ done |
| P1 | Fidelity: no force-deposit; resolved tray wall-jet; multi-jet; gentle empty = 0.0 | ✅ done |
| P1.1 | Low-Mach: diffusive scaling; Ma ≤ 0.1, τ ≥ 0.51; CV band re-frozen 1.05–1.30 | ✅ done |
| P2 | Core promotion: land CR-1/CR-2/CR-3 (§5); example switches from substitutions to real BCs + the core particle layer, reproducing P1.1 numbers within band | Gate: T18.1–.3 green, then example parity |
| P3 | Free surface: single-phase mass-tracking (VOF-on-LBM), NOT a switch to FVM/FEM. Build ONLY if P1/P2 evidence shows the agitation-worsens-uniformity trend requires the interface — do not build speculatively | Gate: reproduce that trend from interface sloshing |
| P4 | Inverse solver: ship discrete-recipe comparison FIRST; then CMA-ES / Bayesian-opt + a response-surface surrogate | Gate: recover a known-good recipe on a synthetic n\* (T18.5) |

## 5. P2 core requirements (CR-1/2/3)

Each CR gets an adversarial acceptance test authored from this spec +
VALIDATION.md T18 by a SEPARATE order/worktree from its implementation.
Scope for all three: CPU scalar + CpuSimd backends with T13 partition
invariance preserved; the GPU backend must **reject** specs using these
features with a `SpecError` (honest failure, no silent wrong physics) until a
GPU follow-up is scheduled.

### CR-1 Localized interior volume source/sink

Motivation: withdraw needs a suction outlet at an interior region; eject needs
jet inflow not tied to a whole lattice face. `GlobalSpec` today exposes only
whole-face BCs.

Frozen public API surface (signatures may be refined at review, semantics may
not):

```rust
// params.rs
pub struct SourceRegion { pub lo: [usize; 3], pub hi: [usize; 3] } // inclusive cell box
pub enum SourceKind<T: Real> {
    MassFlow { q_lu: T },              // mass per lattice step, uniform over region; negative = sink
    Jet      { q_lu: T, u: [T; 3] },   // mass injection carrying prescribed velocity u
}
pub struct VolumeSource<T: Real> { pub region: SourceRegion, pub kind: SourceKind<T> }
// GlobalSpec gains: pub sources: Vec<VolumeSource<T>>,
```

Validation rules: region strictly interior (≥ 1 cell from every face), inside
the domain, no overlap with solids or another source, |u| ≤ MAX_SPEED, and a
sink strength bound such that local ρ stays positive (frozen: per-cell drain
> −1.0 per step). Mass ledger contract: d(total_mass)/step = Σ q_lu to
summation round-off — the achievable relative error is bounded by the
cancellation of two O(N_cells) sums (≈ N·ε·M/|Σq|), NOT by 1e-12; see
PHYSICS.md "T18 first-measurement reconciliation". Jet momentum contract:
the equilibrium-shaped injection delivers dP/step = q_lu·u, measurable only
inside the pre-wall-contact window of a closed box.

Acceptance = **T18.1**: a point-like sink in a filled closed box reproduces
the analytic incompressible sink far-field u_r(r) = q/(4πr²) within band on
radial shells away from walls; the global mass ledger holds exactly; a Jet
source delivers the prescribed momentum flux within band.

### CR-2 Per-cell (masked) boundary condition on a face

Motivation: inlet patches and an outlet on the SAME top face. Today one face
is one `FaceBC` and open faces are limited to one axis.

Frozen public API surface:

```rust
// params.rs
pub struct FacePatch<T: Real> {
    pub face: usize,        // 0..6, same indexing as GlobalSpec::faces
    pub lo: [usize; 2],     // inclusive in-face coords (remaining axes, ascending order)
    pub hi: [usize; 2],
    pub bc: FaceBC<T>,      // override inside the rect; base face BC applies outside
}
// GlobalSpec gains: pub face_patches: Vec<FacePatch<T>>,
```

Validation rules: patch inside face bounds; patches must not overlap; the
existing one-open-axis rule applies to the union of base faces and patches;
velocity/pressure parameter limits as for whole faces.

Semantics frozen at first measurement (2026-07-06, PHYSICS.md "T18
first-measurement reconciliation"): patch rects are GLOBAL in-face
coordinates (translated per subdomain, seam-safe); the non-patch cells of a
bare Closed base face carrying patches are an impermeable zero-velocity lid
(no-BC-at-all diverges); a Closed patch on an open base face is likewise a
lid on its rectangle.

Acceptance = **T18.2**: impinging jet (central Velocity patch, downward) with
a coaxial outlet patch on the SAME top face over a closed floor conserves
global mass at steady state and produces the expected radial wall-jet profile
(stagnation at the axis, off-axis peak in u_r, monotone decay beyond the
peak).

### CR-3 In-core one-way Lagrangian dispersed-phase deposition layer

Motivation: promote the example's particle integrator so it is reusable
(CLI/scenario/inverse loop) and validated. **Extend the existing
`crates/lbm-core/src/particles.rs`** (one-way Schiller-Naumann `ParticleSet`,
engine-agnostic sampler closure) — do not create a parallel module.

Frozen semantics: a deposition-aware step in which any particle whose sub-step
segment crosses the floor plane is removed from the set and recorded at the
interpolated crossing point (deterministic index-order iteration, no
data-race-prone accumulation); suspended = remaining particles; drag stays the
existing semi-implicit Schiller-Naumann form (which reduces to Stokes at
Re_p → 0); step limits are caller-controlled and abort, never truncate
silently.

Acceptance = **T18.3**: (a) single-particle terminal velocity matches the
analytic Stokes/SN settling velocity across a St sweep; (b) deposition-map
determinism under partition invariance — bit-identical deposit records when
the sampled fluid field comes from an unpartitioned vs a partitioned (T13)
run; (c) floor-crossing capture is exact for a straight-line crossing.

Payload convention (frozen 2026-07-06, resolving the ambiguity filed in
SPEC_NOTES_T18_3): `DepositEvent.particle` carries the particle state AT
DEPOSITION — `pos` = the interpolated crossing point (duplicated in
`DepositEvent.pos`), `vel` = the impact velocity used for the crossing step
(the post-drag v_new; physically useful for future resuspension/bounce
criteria, unlike the pre-step state), `exposure` = as accumulated through the
depositing step's start sample.

## 6. Validation anchors (forward-model monotone signs) = T18.4

The forward model must reproduce these monotone trends:

1. ejection rate ↑ → CV ↑
2. agitation present → CV ↑ vs quiescent (quantitative only after P3)
3. fill volume ↑ → CV ↑
4. repeated ejection passes ↑ → CV ↓

Plus regression of the frozen gentle band (§3.5) and the Ma/τ guards.

## 7. Inverse solver (P4) = T18.5

Ship order: (1) discrete-recipe comparison (enumerate a recipe library, score
against n\*), (2) CMA-ES / Bayesian optimization over protocol parameters,
(3) response-surface surrogate to cut forward-run cost. Objective: field
distance to n\* on the M×N partition (L2 by default) + constraint penalties
(Ma/τ guards, suspended-fraction cap). Gate: recover a known-good recipe on a
synthetic n\* generated by the forward model itself.

## 8. File map

- Example: `crates/lbm-cli/examples/dispersed_seeding/` (`main.rs`,
  `protocol.rs`, `particles.rs`, `reservoir.rs`, `readout.rs`,
  `sample_gentle.json`, `sample_harsh.json`, `README.md`,
  `SPEC_FINDINGS.md` = raw findings, order records)
- Acceptance: [VALIDATION.md](VALIDATION.md) §T18
- Core landing sites (P2): `crates/lbm-core/src/params.rs`, `solver.rs`,
  backends (CR-1/2); `crates/lbm-core/src/particles.rs` (CR-3)
