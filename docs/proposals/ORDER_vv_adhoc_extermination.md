# ORDER — Ad-hoc-physics extermination sweep (V&V session)

Commissioned by the D-track PM, user directive 2026-07-06. Addressee: the
"Fluid-structure coupling simulator V&V" session. (Direct message channel was
down at commissioning time; this file is the canonical order.)

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

## Task 2 — Kill plan

For each (C): propose replace-with-resolved-physics or
replace-with-validated-closure. Where neither is possible at current
resolution, classify as a CORE CAPABILITY GAP and prepare the routing package
(scenario JSON, exported fields, metrics, repro command) for the core-engine
session per the routing rule in CLAUDE.md.

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
