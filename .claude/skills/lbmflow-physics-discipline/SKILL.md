---
name: lbmflow-physics-discipline
description: >-
  The mandatory discipline for ANY change that affects computed physical
  behavior in LBMFlow — core, scenario, CLI, examples, demos, or GUI. Use
  BEFORE writing a new model term, constant, branch, clamp, boundary
  treatment, example closure, or validation gate, and AFTER every
  experiment/demo run before reporting its results (behavior-validity
  review). It carries the provenance decision table, the ban list with
  grep-able smells, the two-layer gate template, the stop-rule report
  template, and the escalation table. Designed so a Sonnet-level agent can
  participate in development safely by following it mechanically. Do NOT use
  it to run build gates (lbmflow-build-verify) or dispatch codex orders
  (lbmflow-codex-dispatch).
---

# LBMFlow physics discipline — the developer's contract

**Prime directive (user directive 2026-07-06): physical rigor is everything;
ad-hoc physics is banned.** This Skill is the mechanical procedure that
enforces it. If you follow every checklist here you cannot silently damage
the physics; when a checklist tells you to STOP, stopping and reporting IS
the correct, expected outcome — never work around it.

## Step 0 — does this Skill apply?

Ask one question: **can my change alter any number a simulation computes or
any pattern it produces?** (New term, constant, branch, clamp, BC treatment,
initial condition, particle rule, example closure, output binning — yes.
Pure refactors with bit-identical output, docs, build scripts — no, but the
bit-identity claim must be backed by the equivalence gates in
`lbmflow-build-verify`.)

## Rule 1 — Provenance gate (before writing ANY new physics term)

For **every** new constant, formula, branch, or limiter, walk this table
top-down and take the FIRST matching row:

| The term is... | Then... |
|---|---|
| Resolved from the governing equations by the engine (LBM populations, Guo forcing, half-way walls, …) | OK. No extra artifact needed. |
| A literature-backed closure (drag law, SGS model, settling correlation, …) | Required artifacts, ALL FOUR, before merge: (1) citation + derivation note; (2) stated validity domain (Re/St/resolution range); (3) its own validation test vs an analytic/reference solution; (4) a PHYSICS.md entry (template below). |
| Neither (you would have to invent a constant, tune a coefficient, or add a case-specific branch to make it work) | **STOP. Do not write it.** Emit the stop-rule report (Rule 4). The spec gets revised, not the physics faked. |

PHYSICS.md entry template (copy, fill, no field optional):

```markdown
### <date> <term name> (<file:symbol>)
- Form: <equation as implemented>
- Source: <citation / derivation>
- Validity domain: <parameter ranges>
- Validation: <test file::test name, measured vs reference>
- Replaces / interacts with: <what resolved physics it augments>
```

## Rule 2 — Ban list (grep-able smells)

These patterns are prohibited in physics paths. Before handoff, grep your
diff for each; any hit must be either removed or proven to be outside a
physics path.

| Smell | Grep | Why it is banned |
|---|---|---|
| Case-identity branch | `harshness`, `if.*sample`, thresholds switching model constants per case | The model must not know which test case it is running (Goodhart). |
| Calibrated constant | any float literal in a physics expression without a provenance comment | "It makes the gate pass" is not a derivation. |
| Transport-absorbing clamp | `.clamp(` / `.min(` / `.max(` applied to positions or transported quantities | Silently converts transport into accumulation at the bound (the deposition edge-ring incident). Domain-boundary clamps are only legal as part of a documented wall model. |
| Decorative physics | terms whose removal does not change any gated metric | If it does nothing measurable, delete it; if it does, it needs Rule 1. |
| Silent fallback | `unwrap_or(<physical value>)`, defaulting a physical parameter | Missing physics input is a validation error, not a default. |

## Rule 3 — Two-layer acceptance (when authoring any gate or test)

Every acceptance criterion has TWO layers; a gate with only layer 1 is
incomplete and will be Goodharted:

1. **Band**: scalar metric within a frozen range (CV, L2rel, drift, …).
2. **Behavior anchor**: at least one assertion about the pattern itself —
   sign, monotonicity, spatial structure, symmetry. Templates:
   - Monotone trend: `assert!(cv(rate_hi) > cv(rate_lo), "...")`
   - Spatial structure: assert location of extrema / stagnation / peak
     ordering along a profile (see `t18_2_masked_face.rs` wall-jet profile:
     stagnation minimum at axis → off-axis peak → monotone decay).
   - Symmetry: mirrored input ⇒ mirrored field within stated band.
   - Ablation guard: with closure term X disabled, the gated metric must
     CHANGE by more than band-width (proves the gate actually sees X).

Band governance: tightening is always allowed; loosening requires a recorded
rationale in PHYSICS.md (measured value, why the old band was wrong — the
"T18 first-measurement reconciliation" entries are the worked example).
**A band that was calibrated on top of an unvalidated closure is not
authoritative: killing the closure invalidates the band, and the band gets
re-frozen from resolved physics (PM decision 2026-07-06).**

## Rule 4 — Stop rule (exact report format)

When Rule 1 row 3 fires, or a gate cannot be met without a banned pattern,
stop coding and emit exactly this in your final report:

```
STOP-RULE: gate <T-id / band> is unreachable without unphysical terms.
Attempted: <what physical approaches you tried>
Blocking mechanism: <one-sentence physics of why it cannot pass>
Options for the PM: (a) spec/band revision, (b) resolved-physics capability
needed in core: <what>, (c) validated closure exists in literature: <ref>.
```

This report is a SUCCESS outcome of your order, not a failure.

## Post-run behavior-validity review (MANDATORY before reporting results)

After every experiment/demo run, execute these six steps and attach the
record to your report. A metric passing its band does NOT validate a pattern
no band covers.

1. **Look**: visualize the result (lbmflow-qa-viewer; for the deposition
   demo use its tray variant). Never review from scalars alone.
   **Every experiment run MUST leave a visual artifact** (field PNG,
   cross-section, density map, or a viewer dashboard) alongside its
   metrics — a run that only emits scalars is unreviewable and its results
   may not be reported. Division of labor: agents without a browser (codex)
   GENERATE the artifact and list its path in the report; the reviewing
   Claude session (PM / V&V) does the LOOKING and writes the record.
2. **Name the mechanism**: one sentence — "the pattern is X because Y".
   If you cannot fill Y with physics, escalate (table below).
3. **List active non-resolved terms**: every closure/limiter that was live
   in this run (from PHYSICS.md + your Rule 1 artifacts).
4. **Attribution check**: for each listed term — would the headline pattern
   survive its removal? If unknown and the term is example-local, run the
   ablation (disable the term, rerun, compare); if the ablation is
   expensive, flag it as UNVERIFIED ATTRIBUTION in the report.
5. **Boundary artifact sweep**: inspect the pattern at every clamp, wall,
   seam, and outlet: accumulation exactly at a bound is guilty until proven
   physical.
6. **Record**: append to PHYSICS.md or the track findings file:

```markdown
### <date> behavior review — <run id / sample>
Pattern: <what is observed>
Mechanism: <step-2 sentence>
Resolved vs closure: <which parts of the pattern come from which>
Artifacts checked: <clamps/walls/seams — findings>
Verdict: PHYSICAL | CLOSURE-DRIVEN (needs validation) | ARTIFACT (fix)
Routing: none | core session (data package) | codex order | spec revision
```

## Escalation table (never decide these yourself)

| Situation | Action |
|---|---|
| A gate needs a term Rule 1 forbids | Stop-rule report → PM |
| Observed pattern has no physics explanation you can defend | Behavior review with Verdict=UNKNOWN → PM |
| The fix seems to require core capabilities (resolution, coupling, new BC) | Prepare data package (scenario JSON, exported fields, metrics, repro command) → PM routes to core session |
| A frozen band blocks a physically-correct result | Report measured value + rationale draft → PM (band governance) |
| Two rules of this Skill appear to conflict | Quote both, ask the PM; do not pick one silently |

## Handoff checklist (final, all boxes required)

- [ ] Rule 1 artifacts exist for every new term (or no new terms).
- [ ] Ban-list grep over the diff is clean (or hits justified in the report).
- [ ] Gates green per `lbmflow-build-verify` tier for the touched files.
- [ ] New/changed acceptance has both layers (band + behavior anchor).
- [ ] Behavior-validity review record attached (if anything was run).
- [ ] Every reported run lists its visual artifact path (PNG / VTK +
      viewer / dashboard) — no scalar-only runs.
- [ ] Report is evidence-based: every claim maps to a tool output from this
      session; unverified things are labeled unverified.
