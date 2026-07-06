---
name: lbmflow-user-scaleup-advisor
description: >-
  Pick the governing scale-up / scale-down criterion for a stirred / mixing
  process and translate it into the impeller-speed and geometry relations that
  carry lab conditions to production. Use when the user asks "how do I scale up /
  scale down this reactor", "which scale-up criterion", "will my lab conditions
  transfer to production", "constant P/V or constant tip speed?", "what N for the
  big tank", "why won't bench selectivity reproduce at scale", or describes a
  process (regime, phases, shear, thermal, gas-liquid) wanting a transfer rule.
  Owns the rate-limiting-regime → primary-criterion decision
  table, the primary-match-plus-guardrails output object, the geometric-
  similarity inverse-match relations (P/V, tip speed, Re, Fr), and an honest
  TODAY-vs-roadmap statement of what LBMFlow computes. In SI, not lattice
  units. Do NOT use it to build/run a scenario or pick collision/stability
  (author-scenario / tune-stability); it does not yet fingerprint a flow or match
  geometries (future S-Fingerprint / S-Match Skills).
---

# LBMFlow — scale-up / scale-down criterion advisor

Scale-up fails when the wrong quantity is held constant. You cannot hold
everything constant at once: for geometric similarity, matching power-per-volume,
tip speed, Reynolds number, and Froude number all demand **different** impeller
speeds (see §4). The engineering job is to identify the **rate-limiting regime**,
match the ONE quantity that governs it, and then quantify what the other regimes
give up so the residuals can be checked rather than silently violated.

This Skill turns a process description into a single structured recommendation:
one **primary match** (the governing criterion), a set of **guardrail
constraints** (the other active regimes, expressed as bounds not co-equal
matches), the **residuals to report**, the **items LBMFlow cannot compute today**,
and the **next runs**. It is advisory dimensional engineering — it does NOT
author or run an LBM scenario, and it reasons in SI/process quantities, not
lattice units.

**Output is always the structured object of §2 — never a single group, never
"use your judgment".** Every branch of the decision table is enumerated; a
composite regime (e.g. aerobic + shear-sensitive + exothermic) yields exactly one
primary plus the others as guardrails.

---

## 1. Definitions (v1.1 — dimensionally corrected; do NOT regress)

These are the exact relations the advice is built on. A dimensional bug in v1
(writing `P/V = Np N³ D⁵ / V`, which is missing `ρ_l` and wrong in units) is
**fixed** here. Never regress to that form.

Symbols: `N` = impeller speed [rev/s], `D` = impeller diameter [m], `ρ_l` =
liquid density [kg/m³], `V` = liquid volume [m³], `T_q` = shaft torque [N·m],
`Ω = 2πN` = angular speed [rad/s], `Np` = power number (dimensionless,
geometry/Re-dependent), `μ` = dynamic viscosity [Pa·s], `ν = μ/ρ_l` = kinematic
viscosity [m²/s], `g` = 9.81 m/s².

| Quantity | Correct definition | Units |
|---|---|---|
| Power number | `Np = P / (ρ_l N³ D⁵)` (dimensionless) | — |
| Power from torque | `P = Ω · T_q = 2πN · T_q` | W |
| Power from Np | `P = ρ_l · Np · N³ · D⁵` | W |
| **Power per volume** | `P/V = ρ_l · Np · N³ · D⁵ / V` | **W/m³** |
| **Mean dissipation** | `ε̄ = P / (ρ_l · V) = Np · N³ · D⁵ / V` | **W/kg** |
| Impeller Reynolds | `Re = ρ_l N D² / μ = N D² / ν` | — |
| Impeller Froude | `Fr = N² D / g` | — |
| Tip speed | `v_tip = π N D` | m/s |

**Dimensional guard — never write `P/V = Np N³ D⁵ / V`.** That expression is
missing `ρ_l` and is dimensionally wrong (it yields W/(m³·kg⁻¹)-nonsense).
`P/V` carries `ρ_l` (→ W/m³); `ε̄` does NOT (the `ρ_l` cancels → W/kg). Keeping
`P/V` and `ε̄` distinct — one with `ρ_l`, one without — is the load-bearing check.

`ε̄` is the **tank-mean** dissipation. The **local** dissipation at the impeller,
`ε_local`, is 20–100× `ε̄` and is what governs shear damage and micromixing; do
not use the mean where the local peak is the physics (see the shear and
micromixing rows in §3).

---

## 2. The mandatory output shape — primary match + guardrails

The advice is ALWAYS returned as this structured object. Producing a single
criterion with no guardrails, or telling the user to "use engineering judgment",
is a defect.

```
{
  primary_match:               <ONE governing criterion + the key quantity to match, from §3>,
  guardrail_constraints:       [ <each OTHER active regime as a bound / must-not-violate, not a co-equal match> ],
  residuals_to_report:         [ <what the primary match sacrifices — the quantities that will DRIFT at the new scale and must be checked> ],
  unsupported_or_spec_only_items: [ <quantities the recommendation needs that LBMFlow cannot compute today — routed per §5> ],
  run_next:                    [ <concrete next action: an LBM run that IS supported, a lab measurement, a correlation to apply, or an inverse-match number from §4> ]
}
```

Rules:
- **Exactly one `primary_match`.** If two regimes seem co-primary, the more
  selectivity-/quality-/safety-critical one wins the primary slot; the other
  becomes the tightest guardrail. Name why in one line.
- **Composite regimes are the norm, not the exception.** A real aerobic
  fermentation of a shear-sensitive cell line in an exothermic medium has
  primary = kLa (O₂ is rate-limiting), guardrails = {shear peak ε_local / tip
  speed under the cell threshold, A/V heat-removal adequacy}, residuals =
  {blend time θ_mix grows, Fr/surface behavior changes}.
- **Every guardrail is a bound**, e.g. "v_tip ≤ 2.5 m/s (cell lysis threshold)",
  not "also match tip speed".
- **`unsupported_or_spec_only_items` is mandatory when present** — never silently
  assume LBMFlow computes kLa, torque→Np, θ95, or SGS dissipation (it does not
  today; §5).

---

## 3. Governing-criterion decision table (§1.3 of the order)

Identify the **rate-limiting regime** (the step that gates the outcome), read its
row for the primary match and the key quantity. Rows are mutually exhaustive; the
DEFAULT row catches anything not otherwise classified. In a composite process,
the row whose regime is rate-limiting supplies `primary_match`; every other
matching row contributes a `guardrail_constraint`.

| Rate-limiting regime (trigger) | Primary match | Key quantity to hold / compute |
|---|---|---|
| **Fast homogeneous reaction, selectivity-sensitive** (competing/consecutive kinetics, mixing-sensitive yield; e.g. fast precipitation, nitration, diazo coupling) | **Micromixing** at the feed point | Damköhler `Da = τ_E / τ_rxn` with **engulfment time** `τ_E ≈ 17.2 (ν/ε)^½` evaluated at **ε at the FEED point** (local, NOT tank-mean), plus the **feed-pipe Reynolds** number `Re_feed`. Keep `Da < 1` (mixing faster than reaction). Match feed-zone ε and feed-pipe Re, not P/V. |
| **Slow / kinetically-limited reaction** (reaction is the slow step; mixing is "good enough") | Bulk concentration + temperature + residence time | Match bulk `c`, `T`, and mean residence time `τ_res = V/Q`. Impeller detail is second-order — do NOT over-specify P/V. |
| **O₂-limited gas–liquid** (aerobic fermentation, oxidation, hydrogenation where O₂/gas transfer gates rate) | **Volumetric mass-transfer coefficient `kLa`** | Match `kLa` (⇒ match OTR = kLa·(C*−C_L)). `kLa ∝ (P_g/V)^a · v_s^b` — hold gassed power-per-volume and superficial gas velocity `v_s` together. |
| **Shear-sensitive** (cells, flocs, crystals, emulsions, live tissue) | **Peak local shear** | Peak `ε_local` / **tip speed `v_tip = πND`** / max shear rate `γ̇`, checked against the **particle/aggregate size** and its lysis/breakage threshold. Cap the peak, don't match a mean. |
| **Solid suspension** (keep particles off the bottom / homogeneously suspended) | **Just-suspended speed `N_js`** | Zwietering `N_js = S · ν^0.1 · (g Δρ/ρ_l)^0.45 · X^0.13 · d_p^0.2 / D^0.85`. Match the **suspension state** (all solids in motion), i.e. run at/above `N_js` at each scale — this does NOT follow a simple power law, recompute per scale. |
| **Exothermic** (heat removal gates safe operation / runaway risk) | **Heat-removal capacity** | Surface-area-to-volume `A/V` (falls as `S⁻¹` with scale) + jacket/coil film coefficient `h`. `Q_removable = h·A·ΔT` must exceed reaction heat. A/V shrinking on scale-up is often the binding constraint. |
| **Free-surface / vortex-sensitive** (surface entrainment, vortex formation, gas ingestion at surface) | **Froude number `Fr = N²D/g`** | Match `Fr` to preserve surface/vortex behavior; consider baffling. Governs whether the free surface deforms similarly. |
| **DEFAULT — turbulent, mixing-sensitive but not micromixing-limited** (blend time / general homogenization matters; nothing above dominates) | **Power per volume `P/V`** + **geometric similarity** | Hold `P/V = ρ_l Np N³ D⁵ / V` constant with geometric similarity. **Residual is mandatory:** quantify how **blend time `θ_mix`** grows and how the **feed-zone** conditions change with scale (θ_mix ↑ on scale-up even at constant P/V) — never present constant-P/V as "everything is matched". |

Regime-identification hint: ask (a) what gates the rate — mixing, kinetics, gas
transfer, heat, suspension? (b) is anything fragile to shear? (c) is it
exothermic? (d) does the free surface matter? Each "yes" is a matching row: the
gating one is primary, the rest are guardrails.

---

## 4. Inverse-match relations — the impeller-speed transfer laws (§1.5)

Geometric similarity, scale ratio **`S = D_high / D_low`** (large ÷ small; `S>1`
scale-up, `S<1` scale-down). Given the criterion you chose as primary, the large
scale's impeller speed follows from the small scale's:

| Criterion held constant | Speed relation | Note |
|---|---|---|
| **P/V or ε̄** (constant power per volume/mass) | `N_high = N_low · S^(−2/3)` | The workhorse for the DEFAULT regime. Equivalent to `N_low · (D_low/D_high)^(2/3)`. |
| **Tip speed** `v_tip = πND` | `N_high = N_low · S^(−1)` | For shear-limited. Constant tip speed ⇒ N falls as 1/S. |
| **Reynolds** `Re = ND²/ν` | `N_high = N_low · S^(−2)` | Rarely the right target on scale-up (N drops steeply); mainly a diagnostic. |
| **Froude** `Fr = N²D/g` | `N_high = N_low · S^(−1/2)` | For free-surface/vortex similarity. |

(The exponents are written as `S^(−k)`; the order's §1.5 states the equivalent
inverse form `N_low = N_high · S^(+k)`. They are the same relation read in
opposite directions — confirm which scale is "high" before applying.)

**Conflict statement (mandatory):** these targets are **mutually exclusive under
geometric similarity** — you cannot hold P/V, tip speed, Re, and Fr constant
simultaneously, because each demands a different `N_high`. Scaling up at constant
P/V (`S^(−2/3)`) *raises* tip speed (∝ `S^(+1/3)`) and *lowers* Re-based mixing
intensity; scaling at constant tip speed under-powers the tank per unit volume.
State the trade explicitly whenever you emit an inverse-match number: "matching
X, quantity Y drifts by factor Z."

**Recycle-rig / loop-reactor flag:** when **circulation time (the time to pump the
whole volume once through the high-shear zone) dominates** the outcome — i.e. the
process is governed by how often fluid revisits the impeller/feed zone rather than
by the instantaneous local intensity (large tanks, jet-loop / recycle rigs,
crystallizers with external loops) — the geometric-similarity speed laws above are
**not sufficient**. Circulation-time gradients (well-mixed impeller zone vs.
poorly-mixed bulk) break simple similarity. Flag it: route the circulation-time /
blend-time quantification to the residuals and note it needs a resolved-flow
computation (see §5), not a one-line speed law.

---

## 5. Capability honesty — what LBMFlow computes TODAY vs roadmap

Be plain about this. The advice above is dimensional engineering that stands on
its own, but several downstream quantities it references are **not computable in
LBMFlow's current user surface**. Route any such quantity to
`unsupported_or_spec_only_items` and tell the user how to get it (correlation,
lab measurement, or a future track). Aligns with `docs/skills/b1-capability-map.md`
and the PM annotations in `docs/skills/scaleup/ORDER-v1.1.md`.

| Downstream quantity | Status TODAY | Route |
|---|---|---|
| **Resolved viscous dissipation `ε = 2ν S:S`** (local ε field) | **YELLOW.** Strain-rate/shear-rate gathers landed on the native Rust `Solver` (machine-precision Couette verification); computable at the Rust API level, but the **CLI/scenario output channel is NOT wired yet** (order C queued). Not reachable from a user scenario/CLI run today. | `unsupported_or_spec_only_items`. For feed-point ε / peak ε_local today: use a lab/CFD correlation or the local/mean multiplier (ε_local ≈ 20–100× ε̄). Note the Rust-API path exists for a developer-run. |
| **Rotating-impeller torque → `Np`** | **RED.** Momentum-exchange gives force/torque for **STATIC solids only**; the volume-penalization impeller is a body-force emulation with **no torque**. `Np` from a simulated stirred tank is not obtainable. | `unsupported_or_spec_only_items`. Get `Np` from the impeller vendor curve / standard correlation for the impeller type (Rushton ~5, pitched-blade ~1.3, etc.), not from LBMFlow. |
| **`kLa` (gas–liquid mass transfer)** | **RED.** No gas–liquid interfacial-transfer model in the user surface; 2D Shan–Chen multiphase exists but is not a kLa predictor. | `unsupported_or_spec_only_items`. Use an empirical `kLa ∝ (P_g/V)^a v_s^b` correlation from lab data. |
| **Lagrangian lifelines / circulation-time distribution** | **RED** until the M-F Lagrangian track lands. | `unsupported_or_spec_only_items` (this is what the recycle-rig flag needs). |
| **Scalar-transport blend/mixing time `θ95` (ADE)** | **RED** until the scalar-ADE track lands. No user-facing tracer-transport `θ95`. | `unsupported_or_spec_only_items`. Estimate θ_mix from a correlation (e.g. `θ_mix ∝ N⁻¹` at fixed geometry) for the residual. |
| **LES sub-grid dissipation `ε_sgs`** | **RED** until the LES track lands. Resolved ε only (and that YELLOW per above). | `unsupported_or_spec_only_items`. |
| Steady/transient **velocity field, wall force on a static obstacle**, PNG/CSV/VTK export | **GREEN** (per B1 capability map). | `run_next` — a supported LBM run (author + run via the sibling Skills). |

**Never claim** LBMFlow gives you kLa, an impeller Np from torque, θ95, or an
SGS/local-ε field from a user CLI run today. If the recommendation needs one, it
goes to `unsupported_or_spec_only_items` with the real-world route.

---

## 6. Worked example (composite regime, end-to-end)

**Task:** "We have a 10 L bench aerobic fermentation of a shear-sensitive
mammalian cell line; the broth is mildly exothermic. Rushton impeller,
D_low = 0.06 m, N_low = 300 rpm (5 rev/s). We want to scale to a 1000 L
production tank (geometrically similar). What criterion, and what impeller speed?"

**Step 1 — regimes present.** (a) O₂ transfer gates growth → gas–liquid,
rate-limiting. (b) Cells are shear-fragile → shear-sensitive. (c) Mildly
exothermic → heat removal. Three matching rows in §3; O₂ is rate-limiting ⇒
primary = kLa.

**Step 2 — scale ratio.** Volume ratio 1000/10 = 100. Geometric similarity ⇒
`S = D_high/D_low = 100^(1/3) = 4.64`. So `D_high ≈ 0.28 m`.

**Step 3 — primary (kLa) and its inverse-match.** kLa is held via gassed
`P_g/V` + superficial gas velocity. The `P/V`-constant speed law (§4) gives the
starting impeller speed: `N_high = N_low · S^(−2/3) = 5 · 4.64^(−2/3) = 5 · 0.357
= 1.79 rev/s ≈ 107 rpm`. Gas flow scaled to hold `v_s` (same superficial
velocity ⇒ VVM drops with scale). **kLa itself LBMFlow cannot compute today**
→ `unsupported_or_spec_only_items` (use the `kLa ∝ (P_g/V)^a v_s^b` correlation
from bench data to confirm the P/V that hits the target OTR).

**Step 4 — guardrails.**
- Shear: tip speed at the primary speed = `π·N_high·D_high = π·1.79·0.28 =
  1.57 m/s`. Compare to the cell threshold (constant-P/V scale-up RAISES tip
  speed by `S^(+1/3)=1.67×` vs bench: bench `π·5·0.06 = 0.94 m/s` → 1.57 m/s).
  Guardrail: `v_tip ≤` the line's lysis limit (often ~1.5–2 m/s for CHO). If
  1.57 m/s violates it, the shear guardrail **caps** the speed and you accept a
  lower kLa (or add spargers / go to a lower-shear impeller). Peak `ε_local` is
  the true damage driver — YELLOW in LBMFlow (§5).
- Heat: A/V falls as `S⁻¹` = 1/4.64, i.e. specific heat-transfer area drops
  ~4.6×. Guardrail: `h·A·ΔT ≥ Q_rxn` at 1000 L — check jacket adequacy; a coil
  may be needed. A/V shrinkage is the classic exotherm scale-up trap.

**Step 5 — residuals.** Blend time `θ_mix` grows on scale-up even at constant
P/V (feed-zone homogenization slower); Fr changes (`Fr_high/Fr_low = S^(?)` at
constant-P/V speed → surface behavior shifts). Both go to `residuals_to_report`.

**Step 6 — assemble the object.**
```
primary_match: "kLa (O₂ transfer is rate-limiting) — hold gassed P_g/V + superficial gas velocity v_s constant; starting N_high ≈ 1.79 rev/s (107 rpm) at D_high ≈ 0.28 m via the constant-P/V law N·S^(−2/3)."
guardrail_constraints: [
  "v_tip = πN_high D_high ≈ 1.57 m/s MUST stay ≤ the cell-line lysis threshold (~1.5–2 m/s CHO); if exceeded, cap speed and accept lower kLa or re-engineer sparging/impeller.",
  "Heat removal: A/V drops ~4.6× (S⁻¹); require h·A·ΔT ≥ Q_rxn at 1000 L — verify jacket, add coil if short."
]
residuals_to_report: [
  "Blend time θ_mix increases at constant P/V (slower feed-zone homogenization).",
  "Froude/free-surface behavior shifts vs bench under the constant-P/V speed."
]
unsupported_or_spec_only_items: [
  "kLa — RED in LBMFlow today; confirm target P_g/V from a kLa∝(P_g/V)^a v_s^b correlation on bench data.",
  "Impeller Np (Rushton) — from vendor/standard curve, NOT from LBMFlow torque (RED: impeller is body-force, no torque).",
  "Peak local ε_local for the shear check — YELLOW (resolved ε at Rust-API level only, not user CLI); use ε_local≈20–100×ε̄ or lab/CFD."
]
run_next: [
  "Apply the constant-P/V inverse-match: N_high = 5·4.64^(−2/3) = 1.79 rev/s.",
  "Confirm target P_g/V against the bench kLa correlation for the OTR you need.",
  "A LBMFlow-supported run (velocity field / static-obstacle force / VTK) can characterize the geometry, but does NOT yield kLa/Np/θ95 today (author+run via the sibling Skills)."
]
```

---

## 7. Top failure modes (and the fix)

- **Returned a single criterion with no guardrails / said "use your judgment".**
  The output is defective without the structured object. Fix: always emit
  `{primary_match, guardrail_constraints, residuals_to_report,
  unsupported_or_spec_only_items, run_next}`, every regime enumerated.
- **Wrote `P/V = Np N³ D⁵ / V` (missing ρ_l).** Dimensionally wrong (v1 bug).
  Fix: `P/V = ρ_l Np N³ D⁵ / V` [W/m³]; `ε̄ = Np N³ D⁵ / V` [W/kg] (no ρ_l).
  Keep the two distinct.
- **Used tank-mean ε̄ where local peak governs.** Micromixing and shear damage
  are set by ε at the feed point / at the impeller (`ε_local ≈ 20–100× ε̄`), not
  the mean. Fix: use the local quantity for those rows; τ_E uses feed-point ε.
- **Held constant P/V and called everything matched.** Constant P/V still lets
  blend time and feed-zone conditions drift on scale-up. Fix: the DEFAULT row's
  residual (θ_mix, feed zone) is mandatory.
- **Applied one speed law as if it satisfied all criteria.** P/V, tip speed, Re,
  Fr each need a different N_high. Fix: state the conflict and the drift factor
  whenever an inverse-match number is emitted.
- **Claimed LBMFlow gives kLa / impeller Np / θ95 / SGS-ε from a user run.** It
  does not today. Fix: route to `unsupported_or_spec_only_items` with the real
  route (correlation / vendor curve / lab). Resolved ε is YELLOW (Rust-API only).
- **Confused scale-up vs scale-down direction of S.** `S = D_high/D_low`; verify
  which tank is "high" before applying `S^(−k)`. Fix: state S numerically.
- **Missed the recycle-rig / circulation-time case.** When revisit frequency
  governs (loop reactors, large tanks), the similarity speed laws are
  insufficient. Fix: raise the recycle-rig flag; route circulation-time to
  residuals + roadmap (Lagrangian, RED).
- **Mixed lattice units into the advice.** This Skill is SI/process engineering.
  Fix: keep it dimensional; scenario/lattice work routes to author-scenario /
  tune-stability.

---

## 8. Inventory row (ownership / trigger / non-overlap)

| Field | Value |
|---|---|
| **Ownership boundary** | The rate-limiting-regime → primary-criterion decision table (§3); the primary-match-plus-guardrails structured object (§2); the corrected Np/P/V/ε̄ definitions (§1); the geometric-similarity inverse-match speed laws + conflict/recycle-rig flag (§4); the honest TODAY-vs-roadmap capability routing for scale-up quantities (§5). Reasons in SI/process regime. |
| **Trigger phrases** | "how do I scale up / scale down this reactor", "which scale-up criterion", "constant P/V or constant tip speed?", "will my lab conditions transfer to production", "what N for the big tank", "why won't bench selectivity reproduce at scale". |
| **Explicit non-overlap** | Emits an engineering criterion + speed relations, NOT a scenario. Does NOT choose collision/stability knobs (→ `lbmflow-user-tune-stability`) or build/run a scenario (→ `lbmflow-user-author-scenario`). Does NOT fingerprint a measured flow into a regime, nor numerically match two geometries — those are the future S-Fingerprint / S-Match Skills. Any lattice-unit question routes out (this Skill is SI-only). |
| **Do NOT use when** | The user wants to author/validate/run a simulation (→ author-scenario / run-preset / run-monitor-mcp), fix a divergence or pick bgk/trt (→ tune-stability), extract/plot output files (→ postprocess), OR wants a numerical fingerprint/match of an already-simulated flow (future Skills — say they are not yet available rather than improvising a number). |

---

## 9. Verification gate — the done check

A scale-up recommendation is done when ALL hold:
1. The output is the **complete structured object** of §2 (all five keys
   present; `primary_match` is exactly one; guardrails are bounds, not co-matches).
2. **Every active regime** in the description appears as either the primary or a
   guardrail (nothing dropped); the DEFAULT row's θ_mix/feed-zone residual is
   present when P/V is primary.
3. Any `P/V` or `ε̄` written matches §1 **with `ρ_l` in P/V and not in ε̄**
   (dimensional guard passes).
4. Every inverse-match number cites `S = D_high/D_low` numerically and states the
   **conflict / drift** it implies; recycle-rig flag raised when circulation time
   governs.
5. Every downstream quantity LBMFlow **cannot compute today** (kLa, torque→Np,
   θ95, SGS-ε; resolved-ε YELLOW) is in `unsupported_or_spec_only_items` with a
   real-world route — **no false capability claim**.
