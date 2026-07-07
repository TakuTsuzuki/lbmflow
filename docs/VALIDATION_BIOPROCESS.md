# Bioprocess Validation Plan (VB-01 … VB-08)

Lifecycle: living (owning doc for bioprocess validation groups; supersedes
the T1..T18 matrix which moved to `archive/2026-07-07-pivot/`).

Each VB group specifies **setup**, **QOI**, **acceptance band**, **tier
this validates**, and **current status**. Tests are written adversarially
from this spec (not from the engine), per the historical LBMFlow
"validation-driven development" convention.

Current status legend:

- **Not started** — no ticket in progress.
- **Ticket X open** — implementation ticket in progress or dispatched.
- **Landed / Engineering** — code + quick test green.
- **Engineering GREEN** — full validation band met (includes heavy tests).
- **Evidence GREEN** — Engineering + calibration + holdout + UQ + sensitivity
  records (see [CREDIBILITY_BIOPROCESS.md](CREDIBILITY_BIOPROCESS.md)).

Heavy validation runs are behind `cargo test --release --include-ignored`;
quick smoke tests run by default.

---

## VB-01 — Single-phase stirred tank Np

**Setup.** 3D stirred tank template (BCFD-020) with Rushton or pitched-blade
impeller (BCFD-021) and 4 baffles. Water-like fluid. Two reference
operating points from published correlations (e.g. Rushton Np ≈ 5.0 at
turbulent Re > 10⁴ with standard T:D=3; PBT-45 Np ≈ 1.3 at Re > 10⁴).

**QOI.** Time-averaged Np in the statistically stationary window.

**Acceptance.** ±15% vs published correlations at two independent operating
points, and ±5% mesh convergence between the two finest grids (three-grid
GCI or equivalent).

**Tier this validates.** Engineering (Tier 1) for `single_phase_stirred_tank`
and `power` QOIs (BCFD-031).

**Current status.** Not started. Depends on BCFD-030 / BCFD-031.

---

## VB-02 — Passive-scalar mixing

**Setup.** Same tank as VB-01, at the same rotational speed. Point-pulse
scalar injection near the surface; scalar diffusivity chosen for high Sc.

**QOI.** t95 (time to CV ≤ 0.05·CV₀) and t99 (≤ 0.01·CV₀).

**Acceptance.** Time-step invariance within 5% (halved Δt), monotonic CV
decay after the initial mixing transient, and dimensionally-consistent
Nθ = N · t95 correlation with a published band for the geometry.

**Tier this validates.** Engineering for `passive_scalar` and `mixing`
QOIs (BCFD-034, BCFD-035).

**Current status.** Not started. Depends on BCFD-034 / BCFD-035.

---

## VB-03 — Wall-shear and shear-rate fields

**Setup.** Two configurations exercised in parallel:

- Plane Couette (analytical `du/dy = U/H`).
- Plane Poiseuille (analytical `du/dy = -6·U_max·y/H²` sign convention).

**QOI.** Cell-wise `gamma_dot` field vs analytic; wall-shear-adjacent
tangential-gradient proxy vs analytic; L2 error and L∞ error.

**Acceptance.** L2 ≤ 1e-3 at resolution N ≥ 64, second-order convergence
across N ∈ {32, 64, 128}.

**Tier this validates.** Engineering for `stress` / `shear` fields
(BCFD-032, BCFD-033).

**Current status.** Not started. Depends on BCFD-032 / BCFD-033.

---

## VB-04 — Phase-field droplet and Laplace law

**Setup.** Static 2D and 3D droplets of radius R in periodic domain, with
surface tension σ and no gravity. Sweep R ∈ {8, 16, 32} lattice units.

**QOI.** Laplace pressure jump ΔP measured across the interface vs
analytical σ/R (2D) or 2σ/R (3D). Total-phi drift over 10⁴ steps.

**Acceptance.** Laplace slope within 10% of σ; slope through origin (no
constant offset > 5%); total-phi drift ≤ 0.1% over the run window.

**Tier this validates.** Engineering for `resolved_phase_field` and
`surface_tension` (BCFD-040..043).

**Current status.** Not started. Depends on BCFD-040..043 and BCFD-048.

---

## VB-05 — Sparger gas ledger

**Setup.** 3D tank with resting fluid, ring sparger injecting gas at fixed
volumetric flow Q. Closed lid (no gas escape). Track integrated gas volume.

**QOI.** `Σ (1 − φ) dV` growth rate vs Q.

**Acceptance.** Ledger balances within 2% of expected `Q·t` after the initial
transient. No negative φ. Rejects liquid injection through gas sparger.

**Tier this validates.** Engineering for `sparger` boundary and
`resolved_phase_field` gas injection (BCFD-046, BCFD-047).

**Current status.** Not started. Depends on BCFD-046 / BCFD-047.

---

## VB-06 — Oxygen kLa synthetic

**Setup.** Uniform gas-liquid mixture with fixed interfacial area a and
fixed kL; oxygen scalar initialised to C₀ = 0 with C* > 0. Dynamic gassing
transient.

**QOI.** Fit `dC/dt = kLa·(C* − C)` on the transient; recover the input kLa.

**Acceptance.** Recovered kLa within 5% of input. Fit R² ≥ 0.99.
Equilibrium (C = C*) case fits `kLa ≈ 0` within tolerance.

**Tier this validates.** Engineering for `oxygen` scalar and `kla` QOI
(BCFD-050..052).

**Current status.** Not started. Depends on BCFD-050..052.

---

## VB-07 — Cell shear-exposure integral

**Setup.** Analytical Couette between plates with prescribed uniform
`gamma_dot`. Tracers seeded uniformly. Exposure integral `E = ∫ max(0,
τ − τ_c)^m dt` with a threshold above and below the constant τ.

**QOI.** Percentile distribution P50/P90/P95/P99 of E across tracers; max
E; fraction above threshold; residence time above threshold.

**Acceptance.** Above-threshold case: E per tracer matches analytical
integral within 5% at coarse Δt and within 1% at halved Δt (second-order
time). Below-threshold case: E = 0 exactly for every tracer. Percentile
reducer verified against a synthetic distribution.

**Tier this validates.** Engineering for `cell_tracer`, `shear_exposure`,
and the damage integral (BCFD-060, BCFD-061).

**Current status.** Not started. Depends on BCFD-060 / BCFD-061.

---

## VB-08 — Synthetic scale-up decision

**Setup.** Two synthetic operating maps (small tank and large tank) with
known Np(N), kLa(N, Q_g), P95-shear(N), mixing time(N). Feed the maps into
the scale-up evaluator (BCFD-084) with constraint set: `kLa ≥ target`,
`P/V ≤ limit`, `P95_shear ≤ limit`, `mixing_time ≤ limit`.

**QOI.** Feasible operating window for the large tank; explicit conflict
table if empty.

**Acceptance.** Evaluator recovers the analytic feasible set within
tolerance; correctly reports the tightest constraint when the set is
empty; ranks constraints in the documented priority (constant kLa → P/V
→ tip speed → mixing time, unless user weights override).

**Tier this validates.** Engineering for `scale-up window` (BCFD-084).

**Current status.** Not started. Depends on BCFD-084 and BCFD-083.

---

## Cross-cutting: capability registry snapshot

Every landed VB group updates the capability registry (BCFD-002) so
`lbm capabilities --json` reports its status without doc-reading. Reports
generated by `lbm report` cite the VB status by ID.

## Cross-cutting: heavy tests

Bands `≤ 1e-3`, three-grid convergence, and finite-domain lift/drag tests
run under `cargo test --release --include-ignored`. Quick smoke tests
(setup, unit checks, one grid) run in the default `cargo test --release`
suite.
