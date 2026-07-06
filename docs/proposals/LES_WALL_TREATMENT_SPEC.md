# LES Wall Treatment + Turbulence-Validation Track (MF-β completion)

Status: proposal (design spec). Owner: turbulence-modeling architect.
Scope: the `FR-LES-03` wall-region clause of
[REQ_STIRRED_REACTOR.md](../REQ_STIRRED_REACTOR.md) §4.2 and the Re_τ DNS
acceptance that closes W-LES from "WALE core landed" (T17 row: *"Re_tau DNS
acceptance NOT YET"*) to predictive wall-bounded turbulent CFD.

This document decides the phase-1 wall model, fixes its formulation and its
interaction with the landed solver, and defines a validation ladder mapped to
[VALIDATION.md](../VALIDATION.md)-style bands. It creates no code; it is the
order-generating spec for the codex/adversarial-test landing plan in §5.

## 0. Grounding — what is actually landed (cite before you build)

- **WALE driver** (`crates/lbm-core/src/les.rs`). `WaleLes::update`
  (les.rs:57–114) gathers the full velocity-gradient tensor via
  `solver.gather_velocity_gradient` (les.rs:63), forms the Nicoud–Ducros 1999
  operator (les.rs:71–109), and installs a **one-step-lag** relaxation-rate
  field `omega[i] = 1/(3(ν₀+ν_t)+1/2)` (les.rs:111) through
  `solver.set_omega_field(Some(&self.omega))` (les.rs:113). `Cw = 0.325`
  (les.rs:16). The zero-gradient `0/0` limit is defined as `ν_t = 0`
  (les.rs:104–108) — this is what makes WALE null in laminar shear.
- **omega-field plumbing** (`solver.rs:2489–2519`). `set_omega_field` takes a
  **global compact-order** slice (`(z·ny+y)·nx+x`), asserts length = cell count
  (solver.rs:2492), slices it into each owned subdomain's `omega_field`
  (solver.rs:2498–2514). Collision kernels replace only the local `omega_plus`
  fetch when the field is present; `None` restores the uniform-rate path
  (solver.rs:2486–2488, 2516). This is the sole legal entry point for a
  per-cell relaxation modification — the wall model MUST also write through it.
- **Step pass order** (`solver.rs:1608–1632`, doc-comment 1608):
  `collide → stream → Bouzidi → swap → open faces → moments`. WALE runs
  **outside** `step()`: the test harness pattern is `les.update(&mut solver);
  solver.run(1);` (validation_channel_dns.rs:345–348) — i.e. `update` reads the
  post-moments state of step *n* and installs the field consumed by the collide
  of step *n+1*. The wall model inherits this same "pre-step field install"
  contract; it must not require a new pass inside the fused
  collide+stream+moments kernel (that would re-open the
  `backend_simd_equiv.rs` / T13 gates — see CLAUDE.md core invariants).
- **`tau_eff` used in stress** (`solver.rs:3112`): `tau_eff = 1/omega_p`, the
  *scalar* global relaxation time, structured (comment solver.rs:3133–3135) so
  a per-cell `omega_plus` field replaces the scalar. The parallel `cx/les-clip`
  order adds the FR-LES-03 **upper clip on `τ_eff`** plus diagnostics; this spec
  **assumes that clip lands** and specifies its interaction (§3.4), but does not
  itself implement it.
- **Wall representation** (three coexisting, all relevant):
  1. **1-cell solid rim + half-way bounce-back** — CLAUDE.md invariant: "Wall
     edges are a 1-cell solid rim. Wall surfaces are half-way." The first fluid
     node sits `y_w = 0.5` lattice units off the wall
     (validation_channel_dns.rs:173–177, `yw_from_y`).
  2. **Bouzidi interpolated bounce-back** (`crates/lbm-core/src/bouzidi.rs`).
     Each `BouzidiLink` (bouzidi.rs:13–26) carries the fluid `cell`, direction
     `q`, fractional wall distance `qd ∈ (0,1)` (bouzidi.rs:19), and a
     `wall_ref` neighbour used for wall velocity (bouzidi.rs:25). `qd` is the
     **exact geometric wall distance along the link** — this is the primitive a
     wall-fitted treatment reads (§2). The pass runs after streaming, before
     swap (bouzidi.rs:256, apply_bouzidi_impl), consistent with the step order
     above.
  3. **Moving-wall Guo/Ladd term** already threaded through both BB paths
     (bouzidi.rs:284, `wall_term = 6·w_q·ρ·(c_q·u_w)`).
- **Existing turbulent-channel infrastructure** (unmerged, `cx/chan180`,
  `crates/lbm-core/tests/validation_channel_dns.rs`, 504 lines):
  MKM 1999 Re_τ=178.12 DNS reference table (mean U⁺, validation_channel_dns.rs:14),
  `full_case` = 128×62×72, δ=30, u_τ=0.008, ν=u_τδ/Re_τ, force_x=u_τ²/δ
  (validation_channel_dns.rs:77–94), stat protocol ~50 eddy-turnovers
  (validation_channel_dns.rs:115–123), metrics: mean-profile L2rel over
  5≤y⁺≤150 (validation_channel_dns.rs:322–336), total-stress L2rel vs the
  analytic force-balance line u_τ²(1−y/δ) (validation_channel_dns.rs:—),
  sustained-turbulence guard −⟨u'v'⟩⁺(y⁺≈30) > 0.4. Bands `mean=0.15`,
  `stress=0.10` are **BAND-FREEZE-PENDING(PM)**. This is a **wall-resolved**
  characterization (no wall model): first-node y⁺ = 0.5·u_τ/ν ≈ **2.4**. It is
  the foundation this spec extends and re-uses; the wall-model ladder does not
  re-derive the harness.

## 1. Wall-treatment decision

### 1.1 The M-F use case that drives the decision

Target: single-phase baffled stirred tank, VR-STR-01
(REQ_STIRRED_REACTOR.md:498), matching Reynolds first
(REQ_STIRRED_REACTOR.md:98 priority list). Impeller Re = ND²/ν in the fully
turbulent regime (Re ≳ 10⁴), Ma_lattice = U_tip/c_s ≤ 0.1
(REQ_STIRRED_REACTOR.md:94). The wall region of interest is **the tank wall and
baffles** (stationary, Bouzidi/half-way BB) and **the impeller blade boundary
layer** (moving, IBM/moving-wall BB, FR-ROT-03). Blade boundary layers are thin
and high-shear; tank walls carry the returning circulation.

### 1.2 Option (a): wall-resolved LES, y⁺ ≈ 1 — grid-cost analysis

Wall-resolved LES requires the first off-wall node at y⁺ ≲ 1 and ≳ 15–20 nodes
across the buffer/log region. In LBM with half-way BB the first node is fixed at
y_w = 0.5Δx, so **y⁺ ≈ 1 ⟺ Δx = 2ν/u_τ** at the wall. For the stirred tank take
u_τ ≈ 0.05·U_tip (flat-plate estimate) and U_tip fixed by Ma ≤ 0.1:

- Channel calibration point (already built): Re_τ=178 wall-resolved needs
  δ=30 lattice units for first-node y⁺≈2.4 (validation_channel_dns.rs:80–83).
  Pushing to y⁺≈1 at the same Re_τ needs δ≈72 → ~2.4× per wall-normal direction.
- Stirred tank at Re=2×10⁴, D_tank ≈ 3D_impeller: a wall-resolved boundary
  layer on the tank wall scales the near-wall Δx as ν/u_τ. Resolving y⁺≈1 on all
  wetted walls of a D_tank domain drives the isotropic lattice to
  **O(500³)–O(800³)** cells (1.2×10⁸–5×10⁸). At the CpuSimd throughput implied by
  the channel runtime (§4.6), one flow-through would be **days**, not the
  40–60 min of the Re_τ=178 channel. LBM cannot grade the mesh (uniform lattice
  is a core invariant); wall-normal-only refinement is unavailable without
  multiblock/grid-refinement, which is out of MF-β scope.

**Verdict:** wall-resolved LES is the *validation reference* (channel, §4),
not the *production stirred-tank path*. Reject as the phase-1 stirred-tank
treatment on cost grounds.

### 1.3 Option (b): equilibrium algebraic log-law wall function (LBM-adapted)

Wall-modeled LES (WMLES): place the first node in the log layer
(30 ≲ y⁺ ≲ 150), do **not** resolve the viscous sublayer, and supply the wall
shear stress τ_w from an equilibrium law-of-the-wall solved at that node. In
LBM the established construction (Malaspinas & Sagaut, *J. Fluid Mech.* 700
(2012) 514–542, "Wall model for large-eddy simulation based on the lattice
Boltzmann method", and the consistency-corrected treatment in their 2014
follow-up work) does **not** try to impose τ_w by modifying the bounce-back
populations directly (that over-determines the near-wall moment). Instead it
**sets the effective near-wall relaxation** so that the modeled molecular+SGS
viscosity reproduces the log-law slope: solve

    u⁺ = (1/κ) ln(y⁺) + B,   κ=0.41, B=5.2   (equilibrium log law)

for u_τ from the first-node tangential speed u_∥ and its wall distance y_w
(Newton iteration on u_τ), then form the **wall eddy viscosity** ν_t,wall that,
added at the wall-adjacent node, yields the modeled stress τ_w = ρu_τ². This is
a *closure on ν_t at the wall-adjacent cell*, which composes cleanly with the
WALE ν_t field and the existing `set_omega_field` channel (§3). No population
hacking, no calibrated constant beyond the literature κ, B.

Cost: the first node at y⁺≈30 means Δx ≈ 30ν/u_τ — an O(30×) coarsening of the
wall-normal spacing versus option (a). This is what makes the stirred tank
tractable at O(200³)–O(300³) for Re=2×10⁴.

**Verdict:** this is the phase-1 stirred-tank wall treatment. It matches the
FR-LES-03 "y⁺ wall function" branch verbatim, is literature-backed
(Malaspinas & Sagaut class), introduces no calibrated-to-pass constant
(PHYSICS.md discipline satisfied: derivation = log law; validity domain =
attached equilibrium boundary layer, 30≲y⁺≲300), and reuses the landed
`set_omega_field` path with zero new collision pass.

### 1.4 Option (c): wall-fitted interpolated bounce-back + van Driest damping

Use Bouzidi links (`qd`, bouzidi.rs:19) to place the wall at its true
sub-cell position and apply near-wall SGS damping via a van Driest factor
`(1 − exp(−y⁺/A⁺))², A⁺=26` on ν_t. This *resolves* the sublayer (so it is a
wall-resolving refinement of option a, not a wall model) and its damping only
corrects the WALE ν_t over-prediction very near the wall. It does **not**
relax the y⁺≈1 grid requirement — van Driest is a wall-resolved-LES tool.

**Verdict:** defer to phase-2 as the *blade/curved-wall* refinement (where
Bouzidi is already the FR-ROT-03 boundary and the geometry is curved), and as
the accuracy upgrade for the channel reference. Not the phase-1 stirred-tank
answer. The van Driest damping is, however, adopted in phase-1 as the **blend
that switches off the log-law model in the near-wall/separated cells** (§3.3),
so the two options are not mutually exclusive — (b) is the phase-1 default, (c)
supplies the low-y⁺ limiter and the phase-2 curved-wall path.

### 1.5 Decision

**Phase-1 = equilibrium algebraic log-law wall function (option b)**, applied as
a wall-adjacent eddy-viscosity closure through the existing `set_omega_field`
channel, blended with WALE, with a van Driest-style near-wall cut so the model
deactivates where it is invalid. Wall-resolved LES (option a) is retained
strictly as the DNS-validation reference (channel ladder, §4). Wall-fitted
Bouzidi + van Driest (option c) is the phase-2 curved-wall/blade refinement.

## 2. y⁺ diagnostics design

Two quantities per near-wall fluid cell: **wall distance** y_w and **friction
velocity** u_τ (hence y⁺ = y_w u_τ/ν and u_τ² = τ_w/ρ). Both are computed from
data structures that already exist.

### 2.1 Wall distance y_w

- **Half-way BB rim cells:** y_w = 0.5Δx along the wall-normal link, exactly the
  `yw_from_y` convention already used (validation_channel_dns.rs:173–177). A
  wall-adjacent cell is any fluid cell with ≥1 solid neighbour among the lattice
  directions; the normal is the (averaged) direction(s) toward the solid
  neighbour(s). This is read from the solid mask (`fields.solid`) that the
  gradient gather already consults (solver.rs:3178).
- **Bouzidi links:** y_w = qd·|c_q| along link q (bouzidi.rs:19). This is the
  exact sub-cell distance and is preferred wherever a `BouzidiLink` exists for
  the cell. `wall_ref` (bouzidi.rs:25) gives the wall velocity for the
  tangential-speed projection.

### 2.2 Friction velocity u_τ

Solved at the wall-adjacent node from the tangential fluid speed u_∥ (physical,
Guo-corrected velocity from `sim.ux()` etc. — CLAUDE.md invariant) and y_w by
Newton iteration on the log law u_∥/u_τ = (1/κ)ln(y_w u_τ/ν)+B. Fallback to the
viscous branch u_τ = √(ν u_∥/y_w) when the Newton y⁺ lands below the buffer
(y⁺<11.6, matched-layer switch); this is a **derived** switch (continuity of
the two-layer profile), not a case-keyed branch.

### 2.3 Output contract

A new solver observable, mirroring `gather_strain_rate` / `gather_ux`
(solver.rs:3140, global compact order, solid cells → 0):

```
pub fn gather_wall_metrics(&self) -> Vec<WallCellMetric>;
// WallCellMetric { cell_index: usize (global compact),
//                  y_w: T, y_plus: T, u_tau: T, tau_w: T,
//                  source: WallSource /* HalfwayRim | Bouzidi */ }
```

Only wall-adjacent cells are populated (others absent / zeroed). The wall-model
driver (§3) consumes this internally; the public surface exposes it for the
FR-LES-03 diagnostics requirement and for the validation harness to report the
achieved first-node y⁺ (the number that decides whether the model is even in its
validity domain). This observable is **read-only** and adds no step pass — same
discipline as `gather_velocity_gradient`.

## 3. Formulation of the chosen wall model

### 3.1 Governing equations

At each wall-adjacent fluid cell c with wall distance y_w and tangential speed
u_∥ (relative to the local wall velocity u_w from `wall_ref`/`wall_u`):

1. Solve for u_τ (§2.2), equilibrium log law, κ=0.41, B=5.2.
2. Modeled wall shear stress τ_w = ρ u_τ².
3. Required **total** near-wall eddy viscosity to carry τ_w at that node:
   ν_tot,wall = τ_w / (ρ |∂u_∥/∂n|), evaluated with the one-sided normal
   gradient available from the gathered velocity field.
4. Wall SGS viscosity to add on top of molecular: ν_t,wall = max(0,
   ν_tot,wall − ν₀).
5. **Blend with WALE** (§3.3) → ν_t,eff at cell c.
6. Install through the existing field: `omega[c] = 1/(3(ν₀ + ν_t,eff) + 1/2)`
   — identical algebra to les.rs:111, same `set_omega_field` call.

Interior cells keep the pure WALE ν_t (les.rs:110). The wall model only
overwrites the wall-adjacent-cell entries of the same `omega` vector before it
is installed.

### 3.2 Where it enters the pass order

Exactly at the WALE install point: the wall model is a **post-processing step on
the `omega` vector inside/after `WaleLes::update`**, before
`set_omega_field` (les.rs:113). No change to `step()` (solver.rs:1610), the
fused SIMD kernel, or the Bouzidi pass. Consequently it inherits WALE's
one-step lag (u_τ is computed from step-n velocity, applied to step-n+1
collide) and does **not** touch the bit/threshold gates
(`backend_simd_equiv.rs`, T13) — the only field written is the same per-cell
`omega_plus` those gates already tolerate (verified constant-field bit-identity,
wale_les.rs:48–104).

### 3.3 Blend and near-wall cut (van Driest limiter)

To avoid double-counting SGS stress and to deactivate the model where it is
invalid:

    f_vd(y⁺) = (1 − exp(−y⁺/A⁺))²,  A⁺ = 26
    ν_t,eff = f_vd · ν_t,wall(log-law) + (1 − f_vd) · ν_t,WALE

- At y⁺ ≳ 30 (log region): f_vd → 1, model dominates (WMLES behaviour).
- At y⁺ ≲ 5: f_vd → 0, the model switches off and pure WALE (which is already
  correctly null in laminar sublayer shear, les.rs:104–108) carries the cell.
  This makes the treatment **safe when the grid is accidentally wall-resolved**
  (a wall-resolved run silently reduces to WALE — this is the mechanism behind
  ablation gate §4(iv) and the laminar-regression gate §4(i)).
- κ, B, A⁺ are the standard literature constants (log law; van Driest 1956) —
  no calibration to any acceptance band (PHYSICS.md discipline).

### 3.4 Interaction with the τ_eff upper clip (cx/les-clip)

The clip (assumed landed) bounds `τ_eff` from above with diagnostics. The wall
model **feeds the same `omega`/`τ_eff` field the clip operates on**, so the clip
sees the wall-model contribution automatically. Ordering: wall model writes
`omega` → clip bounds it → `set_omega_field`. The clip's diagnostic counter
should therefore attribute clipped cells; a wall-adjacent cell hitting the upper
clip is a **failure-mode signal** (grid too coarse for even the log law, or a
separated/impinging region where the equilibrium assumption is void — see §3.5).
The `gather_wall_metrics` output (§2.3) and the clip diagnostics together are
the required FR-LES-03 diagnostics.

### 3.5 Validity domain and failure modes

- **Valid:** attached, quasi-equilibrium turbulent boundary layer with the first
  node at 30 ≲ y⁺ ≲ 300; smooth wall; Ma_lattice ≤ 0.1.
- **Failure modes (must be surfaced, never silently absorbed):**
  - *Separation / impingement / stagnation* (u_∥→0): log law gives u_τ→0,
    ν_t,wall→0; the van Driest blend hands the cell back to WALE. Acceptable
    degradation, flagged when a wall cell's y⁺ leaves the valid band.
  - *First node in the buffer/sublayer* (y⁺<30): model over-relies on the
    matched viscous branch; f_vd shrinks the model weight. Report the achieved
    y⁺ distribution — if the median wall y⁺ is outside 30–300 the run is out of
    validity domain and the result is CHARACTERIZED-not-VALIDATED.
  - *Upper-clip saturation* (§3.4): banned as a silent pass; if wall cells
    saturate the clip, STOP and report per the PHYSICS.md stop-rule — do not
    widen the clip to pass a band.
  - *Curved/rotating blade:* equilibrium log law is weakest here; phase-2
    Bouzidi + van Driest (option c) supersedes. Phase-1 uses it only on the
    tank/baffle walls and reports blade-region y⁺ separately.

## 4. Validation ladder (bands, grids, step counts, M5 Max runtime)

All runs D3Q19, CpuSimd f64, WALE default. Reuse the `validation_channel_dns.rs`
harness on `cx/chan180`; the wall-model rungs add a WMLES case and a wall-metric
report. Runtimes estimated from the anchor: the existing full Re_τ=178 channel
(128×62×72 ≈ 5.7×10⁵ cells, ~2.7×10⁶ steps at ~50 T_e) is annotated
"~40–60 min CPU" (validation_channel_dns.rs, `#[ignore]` tag). That fixes the
throughput anchor used below; PM re-measures on landing.

| Rung | Purpose | Grid / steps | Metric → band | Est. M5 Max |
|---|---|---|---|---|
| (i) laminar channel regression | wall model + WALE must NOT perturb Poiseuille | 8×62×4, 12k steps | profile L2rel ≤ **0.02**; max ν_t ≤ **1e-12**; wall-model ν_t,eff ≤ **1e-12** (log-law null in laminar) | < 30 s |
| (ii) turbulent channel Re_τ=180 WMLES vs MKM DNS | primary predictive gate | WMLES coarse: 64×22×36, δ=10, first-node y⁺≈**30**, ~50 T_e (~1.3×10⁶ steps) | mean U⁺ L2rel(5≤y⁺≤150) ≤ **0.18**; total-stress L2rel(0.2≤y/δ≤0.8) ≤ **0.12**; −⟨u'v'⟩⁺ peak ∈ **[0.7,1.0]**; sustained −⟨u'v'⟩⁺(y⁺≈30) > **0.4**; median wall y⁺ ∈ **[30,120]** | ~15–25 min |
| (ii-ref) wall-RESOLVED Re_τ=178 (reference) | anchors the WMLES band; already built | 128×62×72 (existing full_case) | existing: mean ≤ 0.15, stress ≤ 0.10 (BAND-FREEZE-PENDING) | 40–60 min |
| (iii) Re_τ=395 stretch | second Re point, band-transfer check | WMLES: 96×30×54, δ=14, y⁺≈30, ~40 T_e | mean U⁺ L2rel ≤ **0.20**; log-slope 1/κ within **±12%**; sustained-turbulence guard | ~45–75 min |
| (iv) WALE on/off ablation | prove the SGS closure is load-bearing, not decorative | rung (ii) grid, WALE-off vs WALE-on | WALE-off must FAIL (diverge under run_guarded, or centerline U⁺ error > 2× the WALE-on error); WALE-on passes (ii) | 2× rung (ii) |
| (v) dt-halving (one-step-lag) study | bound the WALE/wall-model lag error | rung (ii) at Δt and Δt/2 (u_τ, force, ν rescaled to fix Re_τ, Ma) | mean-profile L2rel difference between Δt and Δt/2 ≤ **0.03** (lag error is O(Δt), must be small vs the DNS band) | ~2× rung (ii) |

Band justifications (VALIDATION.md governance style — provisional numeric band,
frozen after first landing, tightening candidate recorded):

- **(ii) mean 0.18 vs reference 0.15:** WMLES on a coarse grid legitimately
  loses the buffer layer, so a *wider* band than the wall-resolved reference is
  correct, not laxer physics. The band is anchored to the reference rung
  (ii-ref) so the two cannot drift independently.
- **(ii) −⟨u'v'⟩⁺ peak [0.7,1.0]:** DNS peak ≈0.72 near y⁺≈30; a coarse WMLES
  under-resolves the peak location, so the band is one-sided-generous upward
  and floored to reject a laminarized (no-stress) solution.
- **(iii) log-slope ±12%:** the physically meaningful WMLES output is the log-law
  slope 1/κ; testing it directly (not just L2rel) guards against a profile that
  matches in L2 while getting κ wrong.
- **(v) 0.03:** the one-step lag is an O(Δt) explicit-in-time error on the SGS
  and wall closure; halving Δt must move the mean profile by less than the
  gate width or the lag is contaminating the result (this is the concrete test
  the external review's "one-step-lag error" concern demands).

All heavy rungs `#[ignore]`, `LBM_CHAN180_SMOKE`-style env smoke variants for CI
(validation_channel_dns.rs:—, smoke_case). PM schedules the heavy sweep; the
smoke variants run in the normal `--release` suite.

## 5. Phased landing plan (conflict-aware orders)

Convention (CLAUDE.md): implementation orders and adversarial-test orders never
share a worktree; validation tests are written from this spec by codex/Opul,
not by the implementer. Each order = one bundle = one worktree = one background
`codex exec`. `cx/chan180` (harness + reference + MKM data) is a **prerequisite
merge** — land it first so all wall-model rungs build on the existing harness.

- **Order W1 (impl) — wall diagnostics + metric observable.** Files:
  `crates/lbm-core/src/solver.rs` (add `gather_wall_metrics`, wall-adjacency
  detection reusing the solid mask), a new `crates/lbm-core/src/wall_model.rs`
  (WallCellMetric, y_w/u_τ solver, §2). No behaviour change to the step; pure
  read-only observable + Newton u_τ solver. Gate: builds green;
  `gather_wall_metrics` returns 0.5·u_τ/ν first-node y⁺ matching the channel
  hand-computation. **Conflict:** touches solver.rs — must merge before/after
  W2 serially, not in parallel with it.
- **Order W2 (impl) — log-law wall closure + van Driest blend.** Files:
  `crates/lbm-core/src/les.rs` (extend `WaleLes` with an optional
  `WallModel` that post-processes the `omega` vector before `set_omega_field`,
  §3), `wall_model.rs`. Depends on W1's observable. Assumes `cx/les-clip`
  merged (feeds its clip). Gate: laminar-null (rung i) passes inside the impl
  worktree; WMLES coarse channel runs without divergence under `run_guarded`.
  **Conflict:** les.rs + solver.rs — serialize after W1.
- **Order T1 (adversarial test) — wall-model validation ladder.** Files:
  extend `crates/lbm-core/tests/validation_channel_dns.rs` (rungs i, ii, iv, v)
  and add the Re_τ=395 case (rung iii). Written from §4 by codex/Opus in a
  **separate worktree** from W1/W2. Gate: rungs (i) green; heavy rungs
  `#[ignore]`, smoke variants green; bands as §4 with `BAND-FREEZE-PENDING(PM)`
  markers.
- **Order T2 (adversarial test, optional phase-2 seed) — Bouzidi + van Driest
  curved-wall rung.** Files: a Taylor–Couette or rotating-cylinder wall-model
  test exercising the Bouzidi `qd` path (bouzidi.rs) with van Driest damping.
  Only after W1/W2/T1 land; seeds the phase-2 curved-wall treatment.

Merge order: `cx/chan180` → `cx/les-clip` → W1 → W2 → T1 (→ T2). PM runs the
build-verify ritual and the smoke rungs at each merge; the heavy DNS sweep runs
once after T1 to freeze the §4 bands.

---

## Decision summary

1. Phase-1 wall treatment = **equilibrium algebraic log-law wall function
   (Malaspinas & Sagaut class)**, applied as a wall-adjacent eddy-viscosity
   closure through the existing `set_omega_field` path — no new step pass.
2. Wall-resolved LES (y⁺≈1) is rejected for the stirred tank on grid cost
   (O(500³+), days/run) and kept only as the channel DNS reference.
3. Wall-fitted Bouzidi + van Driest is deferred to phase-2 for curved
   walls/blades; van Driest's damping is adopted now only as the near-wall cut
   that hands low-y⁺ cells back to WALE.
4. y⁺ diagnostics come from existing structures: y_w = 0.5Δx (half-way rim) or
   qd·|c_q| (Bouzidi links, bouzidi.rs:19); u_τ by Newton on the log law;
   exposed via a read-only `gather_wall_metrics` observable.
5. The model post-processes the `omega` vector inside `WaleLes::update` before
   `set_omega_field` (les.rs:113), inherits the one-step lag, feeds the
   `cx/les-clip` τ_eff clip, and touches no bit/T13 gate.
6. Validity domain 30≲y⁺≲300 attached equilibrium BL; failure modes
   (separation, sub-buffer node, clip saturation) are surfaced via diagnostics,
   never silently absorbed (PHYSICS.md stop-rule).
7. Validation ladder = laminar regression, Re_τ=180 WMLES vs MKM DNS (anchored
   to the landed wall-resolved reference), Re_τ=395 stretch, WALE on/off
   ablation, Δt-halving lag study — bands and M5-Max runtimes tabulated in §4.
8. Landing = 4 conflict-aware orders (W1 diagnostics, W2 closure, T1 tests, T2
   phase-2 seed), impl/test worktrees separated, merged after
   `cx/chan180` + `cx/les-clip`.

File: `docs/proposals/LES_WALL_TREATMENT_SPEC.md`
