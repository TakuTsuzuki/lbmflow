# ORDER — Physics-rigor risk sweep (V&V session)

Commissioned by the D-track PM, user directive 2026-07-06 (refined same day:
the goal is to IDENTIFY every place that risks compromising physical rigor
and hedge each with MINIMAL EFFORT — prioritized triage, not blind mass
deletion). Addressee: the "Fluid-structure coupling simulator V&V" session.
This file is the canonical order. The executable discipline all fixes must
follow is `.claude/skills/lbmflow-physics-discipline/SKILL.md`.

## Standing principle (read first)

CLAUDE.md / AGENTS.md "Working discipline" now carry the prime directive:
**physical rigor is everything; ad-hoc physics is BANNED** everywhere in the
stack. Every physical behavior must be either resolved from the governing
equations or a literature-backed closure with derivation + validity domain +
its own validation test (PHYSICS.md entry mandatory). Prohibited outright:
constants calibrated to pass an acceptance band; branches keyed to
sample/case identity; position clamps/caps that silently absorb transport;
decorative physics terms.

## Task 1 — Inventory (FIRST; do not modify code yet)

Sweep `crates/` (core, scenario, cli incl. examples) and `web/` for every
physics-affecting term. Classify each:

- **(A) resolved physics** — OK.
- **(B) documented closure** — verify a PHYSICS.md entry + a validation test
  exist; flag if either is missing.
- **(C) AD-HOC — to kill.**

Known seed offenders (all in
`crates/lbm-cli/examples/dispersed_seeding/particles.rs` — `sample_tray` and
the step loop):

1. The `harshness` branch: `wall_jet_len` switches 0.42 vs 0.10 of tray width
   BY SAMPLE IDENTITY.
2. The analytic jet-gaussian / wall-jet closure added ON TOP of the resolved
   LBM tray field.
3. The lateral dispersion constant `2.5e-5 m^2/s` (SPEC_FINDINGS admits it is
   calibrated to the acceptance gate, not validated).
4. The side-wall position clamps that trap particles at the rim — this is
   what makes the gentle deposition map ring at the edges/corners
   (behavior-validity finding; origin note in CLAUDE.md).
5. The deterministic `sin()` pseudo-agitation kick in the harsh branch.

## Task 2 — Minimal-effort hedge plan

For each (C) item, assign a risk tier (does it shape a REPORTED result? a
frozen band? only an unreported internal detail?) and pick the CHEAPEST
sufficient hedge from this ladder, in order:

1. **Delete** — if the term is decorative (removal changes no gated metric),
   remove it outright.
2. **Validate & document** — if a literature basis exists, add the four
   Rule-1 artifacts (citation/derivation, validity domain, validation test,
   PHYSICS.md entry) and keep the term. Often the cheapest real fix.
3. **Ablation-guard** — if the term's influence on reported results is the
   risk, add an ablation behavior anchor (term off ⇒ gated metric must move)
   so its effect is at least visible and tracked, and file the follow-up.
4. **Replace with resolved physics** — when the resolved field/solver can
   already provide the behavior (post-CR-1/2/3 this is true for more of the
   deposition demo than when the closures were written).
5. **CORE CAPABILITY GAP** — when none of the above is possible at current
   resolution, prepare the routing package (scenario JSON, exported fields,
   metrics, repro command) for the core-engine session per the routing rule
   in CLAUDE.md.

Bands frozen on top of a (C)-class term are NOT authoritative: flag them for
re-freezing from resolved physics (PM decision 2026-07-06 — e.g. the gentle
CV band 1.05–1.30 was calibrated on the wall-jet closure).

## Task 3 — Continuous V&V loop (after PM triage of the inventory)

Experiment matrix → visualize (lbmflow-qa-viewer; a tray-demo variant exists
in its assets) → behavior-validity review (spatial patterns / trends / signs,
not just band metrics) → findings ledger → route per the rule (core defect →
core session with data; demo defect → PM dispatches codex; spec defect →
spec revision + PHYSICS.md rationale).

## Deliverable

`docs/proposals/adhoc-inventory-2026-07-06.md`: per item — location,
classification, evidence, kill plan, priority. Report to the D-track PM
session when ready; the PM triages before any code change is ordered.

---

## PM ack log

- 2026-07-06 (D-track PM): V&V physics-validity note on near-neutral BBO
  terms ACCEPTED and frozen into DISPERSED_DEPOSITION.md §5.5 (validity
  floor |ρ_p/ρ_f − 1| ≥ 0.05 provisional; below it the model + validation
  target is the tracer + settling-slip limit). The offered tracer-limit
  audit row (v → u relaxation in nonuniform flow) is ACCEPTED — please add
  it when the Phase B spec lands. Order E confound handling (ANOM-P4-007,
  revised tau-fingerprint/h²-intercept probes) endorsed; verdict routes to
  the core session, cc D-track PM. (Message channel was down; this file is
  the ack.)
