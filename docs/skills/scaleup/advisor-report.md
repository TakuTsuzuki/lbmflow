# SU-Advisor — Authoring Record + Eval Contract

**Session: SU-advisor-author** (LBMFlow SCALEUP v1.1 workstream).
Worktree: `/Users/taku/projects/lbmflow-wt-su-advisor`, branch `skills/su-advisor`.
Date: 2026-07-05. All artifacts English. PM (Fable) owns gates and merges.
Deliverable Skill: `.claude/skills/lbmflow-user-scaleup-advisor/SKILL.md`.

This Skill is **core-independent and GREEN today**: it is dimensional scale-up
engineering (regime → criterion → impeller-speed law), with **no solver
dependency**. It does not author, validate, or run any LBM scenario; it only
*references* which downstream quantities the solver can/cannot compute (§5 of the
Skill) so it never over-promises.

---

## 0. skill-creator / validator record

| Item | Finding |
|---|---|
| skill-creator available? | Yes (`anthropic-skills:skill-creator`), no generator script — hand-authored `SKILL.md` per documented anatomy. |
| Skill output path | `.claude/skills/lbmflow-user-scaleup-advisor/SKILL.md` |
| Validate command | `python3 -m scripts.quick_validate <skill-dir>` run from the skill-creator dir (pyyaml present locally this session). |
| Validator result | **`Skill is valid!`** (exit 0). Note: the validator enforces `description ≤ 1024 chars`; the description was tightened to fit while keeping all trigger phrases + the non-overlap/roadmap clauses. |
| Structure | Single `SKILL.md` (no `references/` needed — the decision table, definitions, and relations are all load-bearing and belong in the main body for a smaller model). |

---

## 1. Authoring decisions (traceability to the order)

- **Definitions block (§1) is v1.1-corrected and guarded.** `P/V = ρ_l Np N³ D⁵ / V`
  [W/m³] carries `ρ_l`; `ε̄ = Np N³ D⁵ / V` [W/kg] does not. An explicit
  dimensional guard forbids the v1 bug `P/V = Np N³ D⁵ / V`. Cross-checked against
  the PM annotation in `ORDER-v1.1.md` (§1.4 verified matching REQ_STIRRED_REACTOR
  rev.4: `Np = P/(ρ N³ D⁵)`, `P = Ω T_q`).
- **Local vs mean ε** is called out because micromixing (τ_E at the feed point)
  and shear damage (peak ε_local) are set by the *local* value, `ε_local ≈
  20–100× ε̄`. The Skill never lets the mean stand in for the local peak.
- **Mandatory output object (§2).** `{primary_match, guardrail_constraints[],
  residuals_to_report[], unsupported_or_spec_only_items[], run_next[]}` — enforced
  as the ONLY acceptable output shape; a single group or "use your judgment" is a
  defect. Composite regimes → exactly one primary + others as guardrail bounds.
- **Decision table (§3)** enumerates all eight rows from order §1.3: fast-homogeneous
  → micromixing (Da, τ_E≈17.2(ν/ε)^½ at feed point + Re_feed); slow → bulk c/T/τ_res;
  O₂-limited → kLa; shear → peak ε_local/tip speed/γ̇ vs particle size; suspension →
  N_js (Zwietering); exothermic → A/V + film h; free-surface → Fr; DEFAULT → P/V +
  geometric similarity with **mandatory** θ_mix/feed-zone residual.
- **Inverse-match (§4)** carries all four laws with `S = D_high/D_low`: P/V→S^(−2/3),
  tip→S^(−1), Re→S^(−2), Fr→S^(−1/2); the mutual-exclusivity conflict statement;
  and the recycle-rig / circulation-time-dominance flag. The Skill notes the
  order's §1.5 inverse form (`N_low = N_high·S^(+k)`) is the same relation read the
  other direction, and requires stating the direction of S numerically.
- **Capability honesty (§5)** aligned to `b1-capability-map.md` + PM annotations:
  resolved ε = 2νS:S is **YELLOW** (Rust-API strain-rate gathers landed; CLI/scenario
  output channel not wired — order C queued); torque→Np is **RED** (static-solid
  momentum-exchange only; impeller is body-force, no torque); kLa, Lagrangian
  lifelines/circulation-time, scalar-ADE θ95, LES ε_sgs all **RED** until M-F tracks
  land. All routed to `unsupported_or_spec_only_items` with a real-world alternative.
- **Sonnet-parity checklist:** explicit decision table (no judgment calls); one
  end-to-end **composite** worked example (aerobic + shear + exothermic, §6);
  top failure modes with fixes (§7); inventory row with ownership/trigger/non-overlap
  and do-NOT-use vs `lbmflow-user-tune-stability` and the future S-Fingerprint/S-Match
  Skills (§8); a verification gate (§9).

---

## 2. Public examples (should-trigger; illustrative — NOT held-out)

These illustrate correct behavior. The **held-out** eval tasks are authored
separately and author-blind by codex from the contract in §3 below. I did not
write, see, or run any held-out task.

### PE-1 — DEFAULT regime (blend-time-sensitive, single-phase)

**Task:** "Scale a blending tank from 50 L (D=0.1 m, N=200 rpm) to 5000 L,
geometrically similar. Nothing fragile, single liquid phase, just need it well
mixed. What speed?"

**Expected:** DEFAULT row → primary = **P/V constant + geometric similarity**.
Volume ratio 5000/50 = 100 ⇒ `S = 100^(1/3) = 4.64`, D_high ≈ 0.46 m.
`N_high = N_low·S^(−2/3) = (200/60 rev/s)·4.64^(−2/3) = 3.33·0.357 = 1.19 rev/s
≈ 71 rpm`. **Mandatory residual:** θ_mix grows on scale-up even at
constant P/V; feed-zone homogenization slower — report it, do NOT declare
"everything matched". `unsupported_or_spec_only_items`: θ95/blend-time (RED,
scalar-ADE) — estimate from θ_mix∝N⁻¹ correlation. Structured object with all
five keys.

### PE-2 — Micromixing (fast homogeneous, selectivity-sensitive)

**Task:** "A fast competitive reaction (semi-batch, feed at the impeller) gives
95% selectivity in the 2 L lab reactor but we're worried about the 2000 L plant.
Which criterion?"

**Expected:** fast-homogeneous/selectivity row → primary = **micromixing at the
feed point**: `Da = τ_E/τ_rxn`, `τ_E ≈ 17.2(ν/ε)^½` at **feed-point ε** (local,
NOT tank-mean), plus **feed-pipe Re**. Keep `Da < 1`. Guardrails: feed location /
feed-pipe design; local ε at feed. Residuals: bulk P/V not the target — matching
tank-mean P/V would MISS the feed-zone micromixing. `unsupported_or_spec_only_items`:
resolved feed-point ε (**YELLOW** — Rust-API only, no CLI channel) → use
correlation/CFD or ε_local≈20–100×ε̄. Correctly rejects "just hold P/V".

### PE-3 — Composite: aerobic + shear-sensitive + exothermic

**Task:** the §6 worked example (10 L → 1000 L mammalian-cell aerobic
fermentation, Rushton, shear-fragile, mildly exothermic).

**Expected:** primary = **kLa** (O₂ rate-limiting); guardrails = {v_tip ≤ cell
lysis threshold — note constant-P/V scale-up RAISES tip speed by S^(+1/3); heat
removal h·A·ΔT ≥ Q_rxn with A/V falling as S⁻¹}; residuals = {θ_mix grows, Fr/
surface shifts}; `unsupported_or_spec_only_items` = {kLa RED→correlation, Np
RED→vendor curve, peak ε_local YELLOW}; `run_next` = {apply N_high = 5·4.64^(−2/3)
= 1.79 rev/s, confirm P_g/V vs bench kLa correlation, a supported LBM velocity/
force run that does NOT yield kLa/Np/θ95}. Exactly one primary; the other two
regimes are guardrail bounds.

### PE-4 — Solid suspension (Njs, non-power-law)

**Task:** "Keep catalyst particles suspended when I go from bench to 10× volume."

**Expected:** suspension row → primary = **N_js (Zwietering)**; key point is that
just-suspended speed does NOT follow a simple S^(−k) law — **recompute N_js per
scale** from the correlation (d_p, Δρ, X, D). Guardrails: any shear cap on the
particles; power draw. Correctly does NOT apply the P/V or tip-speed inverse law
as if it were the suspension criterion.

---

## 3. Eval-harness CONTRACT (hard gates — no held-out tasks here)

**Conventions (per A-pilot §5 / B2 §4):** Baseline = same prompt, no Skill.
Held-out tasks authored separately/adversarially and author-blind (codex). Each
assertion is objectively checkable against the model's produced structured object
/ text (string / field-presence / numeric-equality within tolerance). A
**non-discriminating-assertion guard** names the assertions that MUST fail
baseline to prove the Skill carries the expertise.

| ID | Assertion | Threshold / artifact |
|---|---|---|
| SA-1 | **Exact `primary_match` equality**: the criterion named equals the reference criterion for the task's rate-limiting regime (micromixing / bulk-kinetic / kLa / peak-shear / Njs / A-V-heat / Fr / P/V-default). | string/enum match to the reference primary; exactly ONE primary emitted |
| SA-2 | **`guardrail_constraints ⊇ reference required set`**: every other active regime in the task appears as a guardrail **bound** (not a co-equal match, not dropped). | superset check; each guardrail phrased as a bound (≤/≥/cap), 0 co-primary duplicates |
| SA-3 | **Unsupported items correctly flagged**: every quantity the recommendation needs that is RED/YELLOW today (kLa, torque→Np, θ95, ε_sgs, resolved-ε) appears in `unsupported_or_spec_only_items` with a real-world route; 0 false "LBMFlow computes X" claims. | field-presence + no-false-capability; resolved-ε flagged YELLOW (not GREEN, not RED) |
| SA-4 | **Structured-tuple semantics**: output is the full object `{primary_match, guardrail_constraints[], residuals_to_report[], unsupported_or_spec_only_items[], run_next[]}` — all five keys present; never a single group; never "use your judgment". | 5/5 keys present; 0 "use your judgment" / bare-criterion outputs |
| SA-5 | **Dimensional correctness**: any `P/V` written includes `ρ_l` (→ W/m³) and any `ε̄` excludes it (→ W/kg); the string `P/V = Np N³ D⁵ / V` (no ρ_l) never appears. | regex/units check; 0 occurrences of the v1 bug form |
| SA-6 | **Inverse-match numeric correctness**: when a speed is emitted, it uses `S = D_high/D_low` (direction stated) and the correct exponent per criterion (P/V:−2/3, tip:−1, Re:−2, Fr:−1/2); the conflict/drift is stated; N_js tasks recompute rather than apply a power law. | numeric within ±5%; exponent matches criterion; conflict statement present |
| SA-7 | **Non-discriminating-assertion guard**: on a task where circulation time governs (loop/recycle rig or large tank), the recycle-rig flag is raised and circulation-time routed to residuals + roadmap. | flag present when applicable; else N/A |

**Non-discriminating-assertion guard (MUST fail baseline):** **SA-1, SA-2, SA-4,
and SA-5** must fail the no-Skill baseline. A baseline model characteristically
(a) names a plausible-but-wrong primary or gives several co-equal criteria (fails
SA-1/SA-2), (b) returns a prose paragraph ending in "use engineering judgment"
rather than the structured object (fails SA-4), and (c) reproduces the dimensional
P/V bug or conflates P/V with ε̄ (fails SA-5). These four are the discriminators;
if the baseline passes any of them the task is not diagnostic and must be
re-authored.

### Scoring (per A-pilot §5)

Pass rate = fraction of held-out tasks where ALL applicable hard gates pass.
Target: with-Skill ≥ 0.9 and strictly > baseline, with the guard assertions
(SA-1/SA-2/SA-4/SA-5) failing baseline. Report mean ± stddev over 3 runs/task
(`aggregate_benchmark`), plus time/token deltas vs baseline.

---

## 4. Open questions for PM

1. **Inverse-form direction convention.** The Skill uses `S = D_high/D_low` and
   writes laws as `N_high = N_low·S^(−k)`; the order §1.5 states the equivalent
   `N_low = N_high·S^(+k)`. Both are correct (same relation, opposite reading).
   Confirm the house convention so S-Match reuses one direction consistently.
2. **Resolved-ε YELLOW → GREEN timing.** When order C wires the CLI/scenario ε
   output channel, §5 of the Skill should move resolved-ε from
   `unsupported_or_spec_only_items` to `run_next` (a supported ε field for
   feed-point / peak-shear checks). No Skill change needed until then; flag when C
   lands.
3. **Np source.** The Skill routes Np to vendor/standard correlations (impeller is
   body-force, no torque). If the rotating-solid torque track (RED) ever lands,
   revisit §5.
4. **Overlap with future S-Fingerprint/S-Match.** This Skill deliberately stops at
   "which criterion + which speed law". The numerical fingerprint of a measured
   flow and the geometry-to-geometry match are the future Skills; §8's do-NOT-use
   clause steers there. Confirm the boundary before those land.

---

## 5. Handoff

- **Deliverables:** `.claude/skills/lbmflow-user-scaleup-advisor/SKILL.md`
  (validator: `Skill is valid!`) + this report.
- **Out of scope by design:** held-out eval tasks (codex authors author-blind from
  §3), description-trigger optimization (`scripts/run_loop.py`) — after PM accepts.
- **No solver dependency:** this Skill is GREEN today and needs no core change; it
  only tracks capability status in §5, which the PM annotations already fixed.
