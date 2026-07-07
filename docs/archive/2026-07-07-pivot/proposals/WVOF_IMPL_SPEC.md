# W-VOF Implementation Specification — Conservative Allen–Cahn Two-Phase Model

**Document ID**: SPEC-W-VOF (rev.1, 2026-07-07).
**Owner**: QA-sweep session (per HANDOFF-PM-2026-07-07 §4).
**Scope**: the M-F critical-path item `W-VOF` — the `resolved-phasefield`
interface fidelity default of `docs/REQ_STIRRED_REACTOR.md` §1. Delivers
FR-VOF-01/02/03, gates VR-STR-02 and unblocks W-BCTOP / W-BUB /
active-scalar / interfacial mass transfer.
**Target core**: `crates/lbm-core` (D3Q19 for the phase-field distribution `g`,
carried alongside the D3Q19/D3Q27 hydrodynamic `f`).
**Acceptance**: VALIDATION.md **T17** rows VR-STR-02/03/05/06, plus the
mandatory negative/consistency tests (J_ρ code-path, advected-droplet
conservation, sparger φ=0, well-balanced stratification).

This spec is **executable**: every literature choice is decided and justified,
every code touchpoint is cited against the current worktree, and every gate is
mapped to a T17 row with a provisional band. A follow-on implementation order
should not need to re-derive any design decision here.

> **Discipline note (CLAUDE.md prime directive; PHYSICS.md).** Every closure in
> §1 is resolved from the governing equations or a literature-backed closure
> with a recorded derivation, validity domain, and a dedicated validation test.
> No band-calibrated constant, no case-keyed branch, no transport-absorbing
> clamp appears anywhere in this design. The mandatory PHYSICS.md entry text is
> in §6.5; the behavior-validity review checklist is in §7.6.

---

## 0. Summary of decisions (read this first)

| # | Decision | Justification (short) |
|---|---|---|
| D1 | **Conservative Allen–Cahn** phase field, **Fakhari et al. 2017 (PRE 96, 053301) velocity-based LBE** form of the Geier/Chiu–Lin (2011) conservative AC equation. | Second-order, mass-conserving, no reinitialization, cheaper than a Cahn–Hilliard 4th-order stencil; the REQ already names Fakhari 2017 (§3). |
| D2 | Phase field carried by a **second D3Q19 distribution set `g`** (`GRAV`-parallel, deviation-free), collided with a single-relaxation model at `τ_φ = 3M + 0.5`. | Matches the NFR-01 budget row (D3Q19×2×f32 = 152 B/cell); D3Q19 is isotropy-sufficient for the AC surface term. |
| D3 | **Hydrodynamic coupling via the pressure/velocity-form incompressible two-phase LBE** (Fakhari–Bolster), density `ρ(φ)` linear, viscosity **harmonic-in-μ** (REQ REV-CFD-MJ-013 frozen default). | Well-balanced with the landed `set_gravity` path; harmonic-μ is the REQ-mandated default. |
| D4 | **Surface tension = chemical-potential form** `F_s = μ_φ ∇φ`, with `μ_φ = 4β φ(φ−1)(φ−1/2) − κ ∇²φ`. **CSF `σκn̂δ_s` is validation-only.** | REQ §3 makes the chemical-potential form normative; it is the well-balanced form that cancels against the pressure gradient exactly at equilibrium (kills a class of parasitic currents CSF leaves). |
| D5 | **Gravity stays on the landed host-staged `stage_gravity` path** (`solver.rs:1527`); W-VOF replaces the density factor `rho` with `ρ(φ)` and adds the well-balanced residual `(ρ(φ) − ρ_ref)·g`. | The code comment at `solver.rs:1518-1525` was written to be edited at exactly this line; no new gravity path. |
| D6 | **J_ρ enters continuity and momentum through ONE function** `phase_flux_Jphi(cell) -> [T;3]`, consumed once by the `g` collision (continuity) and once by the momentum source (`J_ρ = (ρ_l−ρ_g) J_φ`). No second discretization exists. | REQ REV-CFD-CR-002 mandate; VR-STR-03/05 negative test proves single-path. |
| D7 | **Per-step force + gradient composition reuses the `update_shan_chen_force` machinery** (`solver.rs:2274-2392`): the `psi_planes` padded-scalar halo exchange, the neighbor stencil, the `force_field` write. `g`-set halo reuses `exchange_f`. | This path already does exactly "exchange one padded scalar plane, run a q-stencil with halo, write force_field" — the AC gradient/curvature/force stencil is structurally identical. |
| D8 | **CPU-first (CpuScalar reference → CpuSimd fused), GPU deferred to a later order** behind the same `Backend::Fields` extension the trait already anticipates (`backend.rs:130-135`). | B-1 is only PARTIALLY RESOLVED (monolithic GPU only, no generic staged multi-set upload); forcing GPU into phase 1 would stall on B-1. |
| D9 | **Phase 1 lands: walls (wetting/contact angle via geometric normal BC), static + advected + rising-bubble validation, sparger φ=0 gas inlet.** Deferred: degassing top boundary (W-BCTOP), swarm/PBM `d_32`, GPU. | FR-VOF-03 sparger is load-bearing and cheap on the AC path; W-BCTOP is a separate DAG node (`REQ §11`). |
| D10 | **Density ratio path staged 1 → 10 → 100 → 10³**, `mixed_safe` precision with `φ, ∇φ, κ, μ_φ, F_s, ρ(φ), μ(φ)` in f64 (NFR-02). Interface width `W = 4`–`5` lattice units default, mobility `M` from `Pe_φ = UW/M`. | Fakhari 2017 demonstrates 10³ at these settings; f64 interface band is the REQ NFR-02 fixed set. |

---

## 1. Governing equations

### 1.1 Conservative Allen–Cahn phase field (decision D1)

We solve the conservative Allen–Cahn (CAC) equation in the **explicit
conservative-flux form** already written into REQ §3:

```
∂φ/∂t + ∇·(φ u) = ∇·( M [ ∇φ − (4/W) φ(1−φ) n̂ ] )      (1)
n̂ = ∇φ / (|∇φ| + ε),   φ ∈ [0,1]   (φ = 1 liquid, φ = 0 gas)
```

Equivalently, with the interface-normal flux collected as `J_φ`:

```
∂φ/∂t + ∇·(φ u + J_φ) = 0,   J_φ = −M [ ∇φ − (4/W) φ(1−φ) n̂ ]    (2)
```

**Why this variant, precisely.** The original Allen–Cahn equation is *not*
mass-conserving (it drives mean-curvature motion). Chiu & Lin (2011, JCP 230,
185) and independently Geier et al. (2015) add the counter-term
`(4/W) φ(1−φ) n̂` so the sharpening flux exactly balances the diffusive flux at
the equilibrium tanh profile, making the equation **conservative** while
retaining the second-order surface term (no 4th-order `∇²(∇²φ)` Cahn–Hilliard
stencil). We adopt the **Fakhari, Mitchell, Leonardi & Bolster (2017), PRE 96,
053301** LBE discretization of (1) — the "velocity-based" formulation — because:

1. It is a *velocity-form* LBE (the equilibrium is built on the hydrodynamic
   `u`, not on a separate chemical-potential relaxation), which couples cleanly
   to our existing Guo-forced `f` velocity moments (`u` already carries the F/2
   term everywhere — CLAUDE.md invariant).
2. It reaches **density ratio 10³** with a diffuse interface of `W ≈ 4–5`
   lattice cells at second-order convergence — the REQ risk-1 target.
3. The REQ (§3, FR-VOF-01) already names "Fakhari 2017" as the governing
   reference; this spec is not introducing a new model.

The equilibrium `tanh` interface profile that (1) preserves is:

```
φ(ξ) = ½ [ 1 + tanh( 2ξ / W ) ]      (ξ = signed distance across interface)   (3)
```

### 1.2 Velocity-based LBE for the phase-field distribution `g` (decision D2)

Carry a D3Q19 distribution `g_i` with a single relaxation rate
`ω_φ = 1/τ_φ`, `τ_φ = 3M + 0.5` (mirrors the hydrodynamic `τ = 3ν + 0.5`
convention, CLAUDE.md invariant; `cs² = 1/3`). The collision is:

```
g_i(x + c_i, t+1) = g_i(x,t) − ω_φ [ g_i − g_i^eq ] + (1 − ω_φ/2) F_i^φ      (4)
```

with the phase-field equilibrium (Fakhari 2017 Eq. 12; velocity-based):

```
g_i^eq = φ · w_i [ 1 + (c_i·u)/cs² + (c_i·u)²/(2 cs⁴) − u·u/(2 cs²) ]         (5)
```

and the **interface-sharpening source** (the AC counter-term, injected as a
forcing so the recovered macroscopic equation is exactly (1)):

```
F_i^φ = (1 − ω_φ/2) w_i · (c_i · [ (4/W) φ(1−φ) n̂ ]) / cs²                   (6)
```

Macroscopic phase field: `φ = Σ_i g_i`. Because `g^eq` sums to `φ` and the
source has zero zeroth moment, `Σ_i g_i` is conserved up to the divergence of
the advective + sharpening flux — i.e. exactly (1). **`g` is stored in ordinary
(non-deviation) form**: unlike `f` (deviation `f_q − w_q` for f32 mantissa
economy, `fields.rs:1-6`), `g^eq` is proportional to `φ ∈ [0,1]` and has no
quiescent rest state to subtract; the interface band is f64 under `mixed_safe`
anyway (NFR-02).

**Gradient / normal / curvature reconstruction.** `∇φ`, `∇²φ`, and hence `n̂`
and `κ = −∇·n̂`, are reconstructed with the isotropic lattice-weight stencils
(Fakhari 2017 Eq. 22–23), reusing the neighbor-with-halo iteration already
proven in `update_shan_chen_force` (`solver.rs:2355-2378`) and the
finite-difference/periodic neighbor logic in `gather_velocity_gradient`
(`solver.rs:3185-3198`):

```
∇φ(x)  = (1/cs²) Σ_{i≠0} w_i c_i φ(x + c_i)                                   (7)
∇²φ(x) = (2/cs²) Σ_{i≠0} w_i [ φ(x + c_i) − φ(x) ]                           (8)
```

These are the same D3Q19 weight-moment isotropic operators the Shan–Chen force
stencil uses (`F = −G ψ Σ_q w_q ψ(x+c_q) c_q`, `solver.rs:2283`), so the
implementation lifts that loop directly with `ψ → φ`.

### 1.3 Hydrodynamic coupling: density / viscosity interpolation (decision D3)

Density and kinematic viscosity are pointwise functions of `φ`:

```
ρ(φ) = ρ_g + φ (ρ_l − ρ_g)                                                   (9)
1/μ(φ) = φ/μ_l + (1−φ)/μ_g       (harmonic-in-μ; REV-CFD-MJ-013 frozen)      (10)
ν(φ) = μ(φ)/ρ(φ)  →  τ(φ) = 3 ν(φ) + 0.5  →  ω(φ) = 1/τ(φ)                  (11)
```

The per-cell `ω(φ)` is installed through the **already-landed per-cell
relaxation-rate field** `omega_field` (`fields.rs:195-197`,
`Solver::set_omega_field` `solver.rs:2489`, consumed by `collide_row`
`kernels.rs:191`). This is exactly the B-6 hook WALE uses (`les.rs:113`
`solver.set_omega_field(...)`). W-VOF composes `ω(φ)` the same way; when LES is
also active the two contributions combine as
`ν_eff = ν(φ) + ν_t`, `ω = 1/(3 ν_eff + 0.5)` before the single
`set_omega_field` call (documented composition order in §4.4).

The hydrodynamic `f` uses the **incompressible pressure-evolution two-phase LBE
(Fakhari–Bolster)**: the `f` equilibrium is built on `(p, u)` with density
`ρ(φ)`, so the recovered momentum equation is REQ §3's

```
∂(ρu)/∂t + ∇·[(ρu + J_ρ) u] = −∇p + ∇·[μ(φ)(∇u+∇uᵀ)] + F_s + ρg + …         (12)
```

In practice, at phase-1 density ratios ≤ 100 we retain the existing
density-weighted `f` moments (`u = (Σ f c + F/2)/ρ`) with `ρ = ρ(φ)` supplied
by the phase field, and add the two-phase forces below as Guo forcing. The
pressure-form refinement (splitting hydrodynamic pressure from `φ`) is required
only on the 10³ leg (§6.4) and lands in the density-ratio order (phase 3).

### 1.4 Surface-tension force (decision D4 — chemical-potential form)

We use the **chemical-potential (potential) form**, which is the REQ-normative
choice (§3):

```
μ_φ = 4β φ(φ−1)(φ−1/2) − κ ∇²φ                                              (13)
F_s = μ_φ ∇φ                                                                (14)
```

with the free-energy parameters *defined* (not calibrated) by:

```
σ = √(2κβ)/6,        W = 4 √(κ / (2β))                                       (15)
```

Given a target physical surface tension `σ` and interface width `W`, (15)
inverts uniquely to `β = 3σ/W`, `κ = 3σW/2` (no free constant — this is a
derivation, not a fit). **Justification for chemical-potential over CSF**: the
potential form `F_s = μ_φ ∇φ` is *well-balanced* — at a static equilibrium
interface `μ_φ` is spatially constant, so `F_s = μ_φ ∇φ = ∇(μ_φ φ) − φ∇μ_φ`
is a pure gradient that the pressure term absorbs exactly, driving the
parasitic (spurious) current toward zero. The CSF form `σ κ n̂ δ_s` requires a
curvature estimate `κ` whose discretization error injects a residual tangential
force that CSF cannot cancel. The CSF form is implemented **only** as an
independent validation cross-check (VR-STR-03 multiphase static droplet must
agree with the potential form), never as the production force.

### 1.5 Well-balanced gravity interaction (decision D5)

The landed gravity path composes `ρ(x)·g` into the per-cell `force_field`
before collision, on fluid cells only, in `stage_gravity`
(`solver.rs:1527-1565`). The method's own doc-comment (`solver.rs:1518-1525`)
states the W-VOF edit verbatim:

> "W-VOF must replace the density factor at this exact line with `rho(phi)` …
> In dynamic-pressure notation the future well-balanced residual is composed
> here as `F_s + (rho(phi) − rho_h) * g + F_b^scalar + …`; single-phase
> compatibility currently uses `rho_h = 0`."

W-VOF therefore:

1. Replaces `let rho = fields.rho[c];` (`solver.rs:1553`) with the phase-field
   density `ρ(φ) = ρ_g + φ(ρ_l − ρ_g)`.
2. Sets the reference density `ρ_ref` (the REQ's `rho_h`) to a chosen constant
   (default `ρ_ref = ρ_g`, i.e. the lighter phase, so the gas headspace is
   force-free), giving the composed gravity term `(ρ(φ) − ρ_ref)·g`.
3. Leaves the surface-tension force `F_s` (14) and any scalar buoyancy
   `F_b^scalar` composed into the same `force_field` accumulation, in the order
   documented in §4.4.

This makes VR-STR-06 (max|u| < 1e-6 in static stratification at ρ ratio 10³)
achievable because the hydrostatic pressure gradient balances `(ρ(φ)−ρ_ref)g`
term-by-term with the potential-form `F_s` (well-balanced by construction).
**No new gravity code path is created**; single-phase runs (`φ` field absent)
keep `ρ_ref = 0` and are bit-identical to the landed contract
(`solver.rs:1525`).

---

## 2. J_ρ consistency requirement (decision D6 — REV-CFD-CR-002)

### 2.1 The mandate

REQ §3 requires the density flux `J_ρ = (ρ_l − ρ_g) J_φ` to appear
**identically** in both continuity and momentum advection:

```
∂ρ/∂t + ∇·(ρu + J_ρ) = 0
∂(ρu)/∂t + ∇·[(ρu + J_ρ) u] = −∇p + ∇·[μ(∇u+∇uᵀ)] + F_s + ρg + …
```

i.e. the *same discrete* `J_φ` that transports the phase field must also
advect momentum. If the two are discretized independently, at ρ ratio 10³ the
inconsistency produces spurious interfacial momentum ("AGG-type" error) and the
droplet-advection mass/momentum budget diverges.

### 2.2 The one shared code path

There is exactly one function that computes the discrete interface flux:

```rust
// solver.rs, new — the SINGLE source of the interface flux.
// Reuses the ∇φ / n̂ stencil of update_shan_chen_force (solver.rs:2355-2378).
fn phase_flux_Jphi(&self, part: usize, x: usize, y: usize, z: usize) -> [T; 3];
//   J_φ = −M [ ∇φ − (4/W) φ(1−φ) n̂ ]      (eq. 2)
```

Two — and only two — consumers call it:

1. **Continuity (phase transport).** The `g`-collision source `F_i^φ` (eq. 6)
   is *derived from* the same `(4/W)φ(1−φ)n̂` term and the same `∇φ`; the AC
   sharpening flux in (2) is `J_φ`. The `g` LBE recovers `∂φ/∂t + ∇·(φu + J_φ)`.
2. **Momentum advection.** `J_ρ = (ρ_l − ρ_g) · J_φ` is added to the momentum
   flux. Concretely: the convective correction to the Guo body force carries
   `∇·(J_ρ u)`; it is composed into `force_field` using the value returned by
   the *identical* `phase_flux_Jphi` call — never a re-derived gradient.

Implementation rule (frozen): the momentum path **must not** recompute `∇φ`,
`n̂`, or `M`; it receives `J_φ` from `phase_flux_Jphi`. A code-review grep
gate (`.claude/skills/lbmflow-physics-discipline` ban-list style) forbids a
second `(4/W)` or second `∇φ` occurrence in the momentum composition.

### 2.3 The test that proves it (VR-STR-03/05, CR-002)

**Advected-droplet conservation (adversarial, authored separately).**
A circular/spherical droplet (`ρ_l/ρ_g = 100`, `W = 4`, `R₀ = 20`) in a fully
periodic box under a **uniform background velocity** `U = (0.05, 0, 0)` is
advected exactly one domain period.

- **Mass drift**: `|∫φ dV(t=T) − ∫φ dV(0)| / ∫φ dV(0) < 0.1%` after one period
  (REQ §8 provisional band; CR-002).
- **Shape/position**: the droplet returns to its start position with L2 profile
  error below the frozen band (no numerical drift from the advective frame).
- **Code-path negative test**: a mutant that uses a *second, independent*
  `∇φ` discretization for the momentum flux (breaking the single-path rule)
  must FAIL this test — this is the "J_ρ consistency code-path" mandatory test
  in REQ §8. The test asserts the production path calls `phase_flux_Jphi`
  exactly once per consumer (instrumented call-count or a snapshot equivalence
  between the flux fed to `g` and the flux fed to momentum).

Galilean invariance under the moving frame is the physical content: only a
consistent `J_ρ` in both equations keeps a translating droplet stationary in
its own frame.

---

## 3. Data-structure mapping (decision D2, D7)

### 3.1 What the landed machinery already supports (verified in code)

- **`Backend::Fields` is an open composite storage boundary**, explicitly
  designed for additional distribution sets. `backend.rs:130-135` (verbatim):
  *"Future multiphase/scalar work can add additional distribution sets (`g`,
  `h`), per-cell properties, and Lagrangian buffers to this associated type
  while the solver continues to transfer through the host staging object."*
- **`SoaFields<T>`** (`fields.rs:167-208`) already holds the q-major padded
  layout `f[q*n_padded + cell]` (halo-padded), compact core moment buffers
  (`rho, ux, uy, uz`), optional `force_field` (compact core,
  `fields.rs:193-194`), optional `omega_field` (compact core,
  `fields.rs:195-197`), `solid`/`wall_u`/`probe` (padded). The layout formula
  is `cell = z·(pnx·pny) + y·pnx + x` (`fields.rs:78-95`, `LocalGeom::pidx`).
- **Padded-scalar halo exchange** exists: `HaloExchange::exchange_scalar(subs,
  planes: &mut [&mut [T]])` (`halo.rs:71`, generic impl `halo.rs:371`), already
  used to exchange the Shan–Chen `ψ` plane (`solver.rs:2322-2333`). The solver
  owns reusable `psi_planes: Vec<Vec<T>>` (`solver.rs:1274`, allocated
  `solver.rs:1460`).
- **Per-cell property fields** (`force_field`, `omega_field`) are compact-core
  and set through dedicated dirty-managed methods (`set_body_force_field`
  `solver.rs:1969`, `set_omega_field` `solver.rs:2489`).

### 3.2 The W-VOF additions to `SoaFields<T>`

Add to `SoaFields<T>` (`fields.rs:167`), all `Option<…>` so `None` is
bit-identical to today's single-phase path (the B-6 invariance discipline):

```rust
/// Phase-field distribution set (D3Q19), q-major padded planes, non-deviation
/// form. `None` ⇒ single-phase (no allocation, bit-identical legacy path).
pub g: Option<Vec<T>>,
/// Ping-pong partner of `g`.
pub gtmp: Option<Vec<T>>,
/// Phase field φ ∈ [0,1], compact core (f64 in the interface band under
/// mixed_safe; see NFR-02). φ = Σ_i g_i.
pub phi: Option<Vec<T>>,
```

Derived per-cell quantities (`∇φ`, `∇²φ`, `n̂`, `κ`, `μ_φ`, `ρ(φ)`, `μ(φ)`)
are **transient per-step scratch**, computed in the force-composition pass and
consumed immediately — they are *not* stored fields (mirrors how Shan–Chen
computes the force into `force_field` each step without persisting `ψ`
gradients). The persistent state is `g`, `gtmp`, `phi`.

The `phi` field is exchanged as a padded scalar plane through the existing
`psi_planes` buffers (rename the concept to `scalar_planes` or add a second
`phi_planes: Vec<Vec<T>>` — a phase-1 order decision; the exchange call is
`solver.rs:2322-2333` verbatim with `ψ→φ`). The `g` distribution halo reuses
`exchange_f`-style pack/unpack; the g-set is D3Q19 so it can reuse
`exchange_f_generic` (`halo.rs:308`) parameterized on `L = D3Q19` regardless of
whether `f` is D3Q19 or D3Q27.

### 3.3 Memory cost per cell (D3Q19 f + g, matches NFR-01)

Per the REQ NFR-01 budget table (§7), fidelity default with ping-pong ×2:

| Component | Layout | B/cell (f32) |
|---|---|---|
| Hydrodynamic `f` (D3Q19 default here; D3Q27 = 216) | 19 × 2 × f32 | 152 |
| Phase-field `g` (D3Q19 × 2) | 19 × 2 × f32 | 152 |
| `φ` + `ρ(φ)` + `μ(φ)`/`ω(φ)` (compact) | 3 × f32 | 12 |
| `∇φ` (transient scratch, band only — amortized) | ~3 × f32 | ~2–12 |
| Interface-band f64 promotion (φ, ∇φ, κ, μ_φ, ρ, μ; ~5–10 % cells) | amortized | +18–37 |

**≈ 320–360 B/cell** for a D3Q19+g phase-1 configuration (the REQ's 540–620
B/cell figure includes D3Q27 f (216), a scalar `h` set, and full statistics
accumulators not in W-VOF phase 1). At 1e8 cells: ≈ 32–36 GB — single-node
feasible on the M5 Max 128 GB dev box for validation at ≤256³ (1.7e7 cells
≈ 6 GB).

---

## 4. Pass structure (decision D5, D7, D8)

### 4.1 The invariant step order (verified `backend.rs:225-286`, `run_span`)

The landed `Backend::run_span` per step is:

```
collide → exchange_f (halo) → stream (interior, then boundary shells)
        → apply_bouzidi → swap → apply_open_faces → apply_volume_sources
        → update_moments
```

(`backend.rs:237-284`; CLAUDE.md invariant: collide → halo → stream → open BCs
→ boundary moments). CpuSimd fuses collide+stream+moments in `step_band`; the
`fused` scratch is `FusedScratch<T>` (`fields.rs:126-164`).

Gravity is composed **before** the backend step by host-staging into
`force_field`: `run_staged_step` (`solver.rs:1587-1606`) does
`stage_out_all → stage_gravity → stage_in_if_dirty → run_span → stage_out_all →
unstage_gravity`.

### 4.2 Where the phase-field passes slot in

W-VOF adds a **phase-field sub-step and a force-composition pre-pass**, both
composed at the solver-orchestration level (the same level as
`update_shan_chen_force` and WALE's `set_omega_field`), so the `Backend` trait
and the invariant `f` pass order are untouched. Per solver step:

```
0. update_moments state from previous step provides ρ, u, φ.
1. PHASE-FIELD PRE-PASS (new, solver-level, before the f collide):
   a. exchange φ padded plane (halo.rs:71 exchange_scalar; ψ_planes reuse).
   b. for each core cell: reconstruct ∇φ (7), ∇²φ (8), n̂, κ; compute
      J_φ = phase_flux_Jphi (§2.2, THE single flux); μ_φ (13).
   c. compose force_field += F_s (14) + (ρ(φ)−ρ_ref)·g  [gravity edit, D5]
      + J_ρ momentum-advection correction [same J_φ, §2.2].
   d. set_omega_field(ω(φ) [⊕ ν_t if LES])  (11), via solver.rs:2489.
2. PHASE-FIELD LBE (new): g collide (4) with source (6) → exchange g halo
   (exchange_f_generic on D3Q19) → g stream → swap g → φ = Σ_i g_i.
3. HYDRODYNAMIC f STEP (unchanged run_span): collide → halo → stream →
   open BCs → moments, now reading the composed force_field and ω(φ).
```

Steps 1–2 are the phase-field analog of `update_shan_chen_force` (which today
runs before `step`); step 3 is the untouched `run_span`. Ordering rationale:
`φ` (hence `ρ`, `μ`, `F_s`) is evaluated from the *current* state and lags the
hydrodynamic step by one sub-step — the same explicit-lag convention WALE uses
(`les.rs:5-7`: "one-step lag"). This lag is a documented model property, not a
hidden approximation (PHYSICS.md entry §6.5).

### 4.3 CPU-first / GPU staging plan (decision D8)

- **Phase 1 (CpuScalar reference).** Implement steps 1–2 against `SoaFields`
  on `CpuScalar` — the reference backend (`backend.rs:1`). All validation
  (§7) runs here first; CpuScalar is the bit-exact oracle.
- **Phase 2 (CpuSimd fused).** Fold the `g` collide+stream into the fused
  `step_band` (`FusedScratch`, `backend_simd.rs`), and the φ-gradient/force
  composition into the pre-pass. Gate: `backend_simd_equiv.rs` bit/threshold
  parity (CLAUDE.md invariant: any pass-structure change must pass
  `tests/backend_simd_equiv.rs` and T13).
- **Phase 3 (GPU, deferred to a separate order).** The `WgpuBackend` mounts a
  `GpuFields` composite (`gpu/backend.rs`), and B-1 is only PARTIALLY RESOLVED
  (monolithic GPU, no generic staged multi-distribution upload —
  SOLVER_IMPROVEMENT_SPEC B-1 STATUS line). A `g`-set GPU buffer + the φ
  gradient shader are a follow-on; they must not block phase 1. The
  `Backend::Fields` associated type already reserves this
  (`backend.rs:130-135`).

### 4.4 Force-composition order (frozen; interaction with cx/gravity-device)

The single `force_field` accumulation, per cell, in this fixed order:

```
force_field[c] = F_s(φ)                       // (14) surface tension
               + (ρ(φ) − ρ_ref) · g           // (D5) well-balanced gravity
               + J_ρ advection correction      // (§2.2) shared flux
               + F_user / F_b^scalar / …        // existing sources, unchanged
```

`ω(φ)` (11), optionally combined with LES `ν_t` as
`ω = 1/(3(ν(φ)+ν_t)+0.5)`, is installed via the single `set_omega_field` call.

**Interaction with `cx/gravity-device`.** As of 2026-07-07 the branch
`cx/gravity-device` points at the same commit as `main`
(`git rev-parse main cx/gravity-device` → identical SHA `29c6304…`); no
divergent device-side gravity composition has landed, and the GPU backend does
no gravity today (grep of `gpu/backend.rs` finds only `rho` copy-out, no
`rho*g` composition). The current gravity path is entirely host-staged
(`stage_gravity`, `solver.rs:1527`). **Therefore W-VOF phase 1 edits the
host-staged path only** (D5), and is orthogonal to any future device-side
gravity work: whenever `cx/gravity-device` lands a device composition, it must
consume the *same* `ρ(φ)` and `ρ_ref` inputs that §4.4 defines — the W-VOF
order records this as a hand-off requirement, not a phase-1 dependency.

---

## 5. Boundary conditions (decision D9)

### 5.1 Phase 1 — lands

- **Walls: wetting / contact angle (FR-BC-03).** The equilibrium contact angle
  `θ` is imposed as a *geometric normal boundary condition* on `∇φ` at the
  solid rim: at a wall cell the interface normal is rotated so that
  `n̂·n_wall = −cos θ` (Ding & Spelt 2007 geometric wetting BC). Concretely, the
  `∇φ` stencil (7) at a fluid cell adjacent to a solid uses a *ghost* φ in the
  solid set by the wetting condition, exactly as `update_shan_chen_force`
  handles solid neighbors today with a virtual value (`solver.rs:2363-2372`,
  the `psi_wall` virtual-wall-density path). This mirrors the landed
  full-range contact-angle mechanism (T11c, `MULTIPHASE_DESIGN.md:62-66`;
  `PHYSICS.md` T11c) but derives θ from geometry, not a calibrated `G_w`. The
  1-cell solid rim and half-way wall placement (CLAUDE.md invariant) are
  unchanged.
- **Sparger gas inlet (FR-VOF-03 — LOAD-BEARING).** The sparger injects **gas**
  (`φ = 0`). Per FR-VOF-03 the schema **never exposes raw φ**: config uses
  `inlet_phase: gas | liquid`; the core maps `gas → φ = 0`, `liquid → φ = 1`,
  enforced by config validation. A gas inlet sets the phase-field boundary to
  `φ = 0` on the inlet face and a velocity Dirichlet on `f` (existing
  `set_inlet_profile` `solver.rs:2218` / Zou–He velocity face). **Plain `φ=0` +
  velocity alone is banned** (FR-VOF-03): the injection must simultaneously
  satisfy gas volumetric-flow conservation and the `d_b/W`, `d_b/Δx` lower
  bounds. Phase 1 implements the **gas-phase volumetric-flow inlet** variant
  (the simplest FR-VOF-03-compliant option): a velocity Dirichlet with `φ=0`
  whose face-integrated `∫(1−φ)u·n dA` equals the prescribed gas volumetric
  flow, validated by the gas-volume balance test (§7, VR-STR-02c precursor).
  Outputs report `φ_liquid` and `α_g = 1 − φ` (FR-VOF-03, FR-IO-01), never raw
  φ at inlets.

### 5.2 Deferred to W-BCTOP (separate DAG node, REQ §11)

- **Degassing / free-surface top boundary (FR-BC-01, W-BCTOP).** `closed` /
  `free-surface` / `degassing-outlet` top faces, headspace pressure, and
  free-surface deformation **wait on W-VOF** by the DAG (`REQ §11`: W-BCTOP hard
  dep = W-VOF). Phase 1 uses a `closed` top (no-slip wall) with gas headspace,
  which is sufficient for the closed static/advected/rising-bubble validation
  set. Degassing is the *first* follow-on order once W-VOF phase 1 is green.
- **Reactive / adsorption scalar walls (FR-BC-04)** — belong to W-SCAL/W-REACT.

---

## 6. Stability & parameter domain (decision D10)

### 6.1 Interface width `W`

Default `W = 4` lattice cells (Fakhari 2017 operating point); `W = 5` for the
10³ ratio leg to widen the resolved band. `Cn = W/L` (Cahn number) kept small
(interface thin relative to droplet: `d_b/W ≥ 4` — the resolvability lower
bound FR-VOF-03/04 references). Below `W = 3` the tanh profile (3) is
under-resolved and mass conservation degrades; below `d_b/W = 4` the droplet is
not a resolved interface and must switch to the point-bubble path (Phase 2).

### 6.2 Mobility `M`

`M = τ_φ-derived`, chosen from the phase Péclet number `Pe_φ = U W / M`
(REQ §2 lattice conventions). Too-large `M` over-diffuses the interface (mass
smearing); too-small `M` under-sharpens (the counter-term cannot hold the
profile) and can make `τ_φ = 3M + 0.5` approach 0.5 (instability). Operating
band: `τ_φ ∈ [0.5, 1.0]` → `M ∈ (0, 1/6]`; default `M = 0.02–0.1` (Fakhari
2017 range). `M` is a physical closure parameter with a validity domain, not a
band-fit constant — its value is reported per run and frozen in PHYSICS.md
after the characterization sweep.

### 6.3 Spurious currents (parasitic) bound

The chemical-potential well-balanced form (D4) holds parasitic currents to
`Ca_spurious = μ_l |u|_spurious / σ < 10⁻³` (FR-VOF-02, dimensionally-correct
form; the old `|u|·L/(σ/μ)` form is void). Literature (Fakhari 2017, Fig. 3–4)
reports `Ca_spurious ~ O(10⁻⁴–10⁻³)` at these `W`/`M` for a static droplet —
consistent with the target band. This is the acceptance bar, not a tuned
result.

### 6.4 Density-ratio path to 10³

Staged, each leg its own gate before proceeding:

1. **ρ ratio 1** (matched density): validates the AC transport + advected-drop
   conservation with no density coupling (pure kinematics).
2. **ρ ratio 10**: turns on `ρ(φ)`, `μ(φ)`, gravity coupling; rising-bubble.
3. **ρ ratio 100**: Laplace + parasitic-current gate at production coupling.
4. **ρ ratio 10³**: requires the pressure-form `f` refinement (§1.3), f64
   interface band (NFR-02), and `W = 5`. VR-STR-06 well-balanced at 10³ is the
   terminal gate. This is REQ risk-1 ("Open — W-VOF pending").

### 6.5 The PHYSICS.md validity-domain statement (mandatory entry text)

The implementation order must add, on landing, a PHYSICS.md §1 stack entry and
a §2 decision entry containing:

> **Two-phase interface — conservative Allen–Cahn (Fakhari 2017), potential-form
> surface tension.** `g` D3Q19 velocity-based LBE, `τ_φ = 3M + 0.5`; `φ = Σ g_i`,
> `φ=1` liquid / `φ=0` gas. Surface tension `F_s = μ_φ ∇φ`,
> `μ_φ = 4βφ(φ−1)(φ−½) − κ∇²φ`, with `β = 3σ/W`, `κ = 3σW/2` (from `σ=√(2κβ)/6`,
> `W=4√(κ/2β)` — derivation, not calibration). Density `ρ(φ)` linear,
> viscosity harmonic-in-μ (REV-CFD-MJ-013). Gravity `(ρ(φ)−ρ_ref)g`, `ρ_ref=ρ_g`
> default, composed on the landed host-staged force path. `J_ρ=(ρ_l−ρ_g)J_φ`
> from a single `phase_flux_Jphi` fed to both continuity (g source) and momentum.
> **Validity domain**: `W∈[4,5]`, `M∈(0,1/6]` (`τ_φ∈[0.5,1.0]`), `d_b/W≥4`
> resolvability (below → point-bubble Phase 2), density ratio ≤10³ demonstrated,
> `Ma_lattice≤0.1`. **Why here (not derivable from code)**: the potential form is
> chosen over CSF for well-balancedness (kills a parasitic-current class);
> record the measured `Ca_spurious` and the frozen `mixed_safe` interface-band
> width (`max(3W,6Δx)` provisional, re-frozen from the W-VOF characterization).

---

## 7. Validation plan mapped to T17 (decision D9)

Tests are **authored adversarially by codex/Opus from this spec**, in a test
worktree that never shares with the implementation worktree (CLAUDE.md team
convention; REQ §8). Each row = metric / reference / band / grid / steps /
backend / pass-fail. Bands are provisional MVP gates (REQ §8 "Provisional
bands"); tightening is always allowed, loosening requires a PHYSICS.md
rationale.

| ID | Test | Metric & band | Grid / steps / backend | T17 row |
|---|---|---|---|---|
| V1 | **Laplace law** | `Δp = σ/R` linear fit R² ≥ 0.999; slope `σ` within ±10 % of the (β,κ)-set target; each droplet `σ=Δp·R` within ±5 % of slope | 128², R₀∈{12,16,20,24}, ρ ratio 100, 40k steps, CpuScalar | VR-STR-03 |
| V2 | **Parasitic currents** | `Ca_spurious = μ_l|u|_spurious/σ < 10⁻³` (We→0, static droplet) | 128² (R₀=32), ρ ratio 100, steady, CpuScalar | VR-STR-03 (`Ca_spurious` fixed) |
| V3 | **Single rising bubble vs Grace** | terminal `U_t` vs Grace Eo–Mo–Re diagram, ±10 % | 3D `d_b/Δx≥20`, ρ ratio 100, until terminal, CpuScalar | VR-STR-02a |
| V4 | **Advected-droplet mass drift** | `|∫φ dV(T)−∫φ dV(0)|/∫φ dV(0) < 0.1 %` over one period; + J_ρ single-path negative test (§2.3) | periodic 128², U=(0.05,0,0), ρ ratio 100, one period | VR-STR-03/05, CR-002 |
| V5 | **Static stratification stillness (well-balanced)** | `max|u| < 10⁻⁶` (LU) in static two-layer stratification at ρ ratio 10³ | 64×128, ρ ratio 10³, W=5, 30k steps, f64 band | VR-STR-06 |
| V6 | **Sparger gas balance (φ=0)** | gas-inlet unit test injects φ=0; `∫(injected gas vol) = Δ(domain gas vol)` closes; no schema field accepts raw φ | small 3D sparger-only, CpuScalar | VR-STR-02c precursor, CR-001 |
| V7 | **Static droplet mass drift** | `< 0.1 % / 1000 steps` | 128² R₀=25, ρ ratio 100, 10k steps | VR-STR-05 |
| V8 | **CSF cross-check** | potential-form static droplet agrees with CSF `σκn̂δ_s` reference within frozen band | 128², CpuScalar | VR-STR-03 |

**Mandatory negative / consistency tests (REQ §8):**

- **J_ρ single-path** (V4 negative arm): a mutant with an independent momentum
  `∇φ` must FAIL — proves §2.2.
- **Sparger phase** (V6): a schema that accepts a raw φ inlet value must be
  rejected by config validation.
- **Well-balanced sign** (V5): flipping the `(ρ(φ)−ρ_ref)g` reference (using
  `ρ_ref=ρ_l` instead of `ρ_g` incorrectly, or dropping `F_s`) must break the
  `max|u|<10⁻⁶` gate — proves the balance is real, not coincidental.

### 7.6 Behavior-validity review (mandatory, REQ / CLAUDE.md)

After each validation run, before reporting: review the *observed* pattern, not
just the gated metric. Specifically: (a) the interface stays a `tanh` profile
of width `W` (not smeared, not oscillating); (b) the rising bubble's wake and
shape (spherical / ellipsoidal / cap) matches the Grace-diagram regime for its
Eo–Mo, not merely the terminal-velocity number; (c) parasitic currents are
localized at the interface and decay, not growing; (d) the static droplet does
not drift or lose volume in a preferred direction (would indicate a
non-conservative flux or a clamp). Record the review in PHYSICS.md or the
track's findings file. A metric passing its band does **not** validate a
pattern no band covers.

---

## 8. Phased landing plan (decision D8, D9)

Five orders, file-conflict-aware. One order = one bundle = one dedicated
worktree (CLAUDE.md team convention). Implementation and adversarial-test
orders never share a worktree.

| Order | Scope | Primary files (conflict boundary) | Gate |
|---|---|---|---|
| **O1 — AC transport core (CpuScalar)** | `g` D3Q19 set + `phi` field in `SoaFields`; `g` collide (4)+source (6); φ=Σg; ∇φ/∇²φ/n̂ stencils (7,8); `phase_flux_Jphi` (§2.2); φ padded-plane halo (reuse `exchange_scalar`). NO density coupling yet (ρ ratio 1). | `fields.rs` (add g/gtmp/phi), `solver.rs` (phase-field pre-pass + g LBE), `halo.rs` (φ plane reuse), `kernels.rs` (g collide row) | V4 (mass drift, ρ ratio 1) green on CpuScalar; single-phase bit-identical with g=None (B-6-style invariance). |
| **O2 — Two-phase coupling + surface tension** | `ρ(φ)`(9), harmonic `μ(φ)`(10)→`ω(φ)`(11) via `set_omega_field`; `μ_φ`(13), `F_s=μ_φ∇φ`(14); gravity edit (D5) in `stage_gravity`; force-composition order (§4.4); J_ρ momentum correction (shared flux). | `solver.rs` (`stage_gravity` edit `:1527`, force composition), reuses O1 stencils | V1 (Laplace), V2 (Ca_spurious), V7 (static drift), V5 (well-balanced, ρ ratio ≤100 leg) on CpuScalar. Depends: O1. |
| **O3 — Boundaries: wetting + sparger** | geometric contact-angle BC on ∇φ (§5.1, virtual-φ ghost like `solver.rs:2363`); sparger `inlet_phase: gas\|liquid` schema + core φ mapping + gas volumetric-flow inlet; `α_g=1−φ` output (FR-IO-01). | `crates/lbm-scenario/src/lib.rs` (schema `inlet_phase`, config validation), `solver.rs` (wetting ghost, gas inlet), `crates/lbm-cli` (α_g output) | V6 (sparger gas balance, CR-001), contact-angle monotonicity. Depends: O2. Separate scenario/CLI files → parallelizable with O4. |
| **O4 — CpuSimd fused + density-ratio 10³** | fold g collide+stream into `step_band` (FusedScratch); pressure-form `f` refinement (§1.3) for 10³; f64 interface band (NFR-02). | `backend_simd.rs`, `fields.rs` (FusedScratch g rings), `solver.rs` (precision band) | `backend_simd_equiv.rs` bit/threshold parity + T13 partition invariance; V5 at ρ ratio 10³. Depends: O2. |
| **O5 — Validation authorship (codex adversarial, separate worktree)** | All of §7 (V1–V8) + the three negative/consistency tests, authored from THIS spec, not from the impl. | `crates/lbm-core/tests/wvof_*.rs` (new files only — no impl-file conflict) | Tests compile red against a stub, go green against O1–O4 as they land; freeze bands in VALIDATION.md T17. Runs alongside O1–O4. |

**Critical-path ordering**: O1 → O2 → {O3 ∥ O4}. O5 runs concurrently from the
start (test worktree). GPU (Phase 3) is a post-W-VOF order gated on B-1 and is
out of this plan's scope.

**Per-order DoD (all orders):** existing tests green *without modification*;
`g=None`/`phi=None` path bit-identical to today (probe_state_hash unchanged
where applicable — B-6 invariance discipline); the phase-1 PHYSICS.md entry
(§6.5) landed with O2; behavior-validity review (§7.6) recorded for every
validation run.

---

## 9. Load-bearing code references (grounding index)

| Claim | File:line |
|---|---|
| `Backend::Fields` reserves g/h distribution sets | `crates/lbm-core/src/backend.rs:130-135` |
| Invariant step order in `run_span` | `crates/lbm-core/src/backend.rs:225-286` |
| `SoaFields` layout, force_field, omega_field | `crates/lbm-core/src/fields.rs:167-208`, `:193-197` |
| q-major padded index formula | `crates/lbm-core/src/fields.rs:78-95` |
| Gravity host-staged composition + W-VOF edit point (doc-comment) | `crates/lbm-core/src/solver.rs:1518-1565` |
| `run_staged_step` (gravity ordering) | `crates/lbm-core/src/solver.rs:1587-1606` |
| Shan–Chen force stencil = template for ∇φ/F_s stencil | `crates/lbm-core/src/solver.rs:2274-2392` |
| `set_omega_field` (per-cell ω, B-6 hook) | `crates/lbm-core/src/solver.rs:2489` |
| `exchange_scalar` padded-plane halo | `crates/lbm-core/src/halo.rs:71`, `:371` |
| `gather_velocity_gradient` (neighbor/periodic reconstruction) | `crates/lbm-core/src/solver.rs:3166-3258` |
| WALE `set_omega_field` composition + one-step lag | `crates/lbm-core/src/les.rs:57-114`, `:5-7` |
| B-1 GPU multi-set limitation (why GPU deferred) | `docs/SOLVER_IMPROVEMENT_SPEC.md` B-1 STATUS |
| `cx/gravity-device` == main (no divergent device gravity) | `git rev-parse main cx/gravity-device` → `29c6304…` |
| Landed contact-angle mechanism (virtual wall density) | `MULTIPHASE_DESIGN.md:62-66`, `docs/PHYSICS.md` T11c |

**Literature (decided references):**
Chiu & Lin 2011 (JCP 230:185) — conservative AC counter-term.
Geier et al. 2015 — conservative phase-field LBM.
Fakhari, Mitchell, Leonardi & Bolster 2017 (PRE 96:053301) — velocity-based
conservative-AC LBE, ρ ratio 10³ (**adopted governing discretization**).
Ding & Spelt 2007 — geometric wetting BC.
Nicoud & Ducros 1999 — WALE (already landed; ω-field composition precedent).
