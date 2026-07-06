# B2 Skill Specs — Track B (Skills for User)

**Session B2, LBMFlow Skills Initiative.** Turns the B1 capability map's GREEN
items into user-facing Skill specs + scaffolds, writes implementation orders for
the one YELLOW item, and spec-only notes for RED items. Single source of
capability truth = `docs/skills/b1-capability-map.md`.
All artifacts English. Branch `skills/b2`.

Base grounding for accepted options (verified against source, not memory):
`crates/lbm-scenario/src/lib.rs` (enums, `validate`, 3D build rejections),
`crates/lbm-cli/src/main.rs` + `mcp.rs` (command/tool surface),
`crates/lbm-cli/src/runner.rs` (output filenames).

---

## 1. Granularity decisions (justified against the B1 stage table)

The candidate list in the brief maps to B1 stages. Decision rule: one Skill per
**user goal with a distinct output artifact and a distinct verification gate**;
fold a stage into another Skill when it is that Skill's done-check rather than a
standalone goal; never exceed green evidence.

| # | Proposed Skill | B1 stage(s) covered | Green evidence | Decision & justification |
|---|---|---|---|---|
| 1 | **lbmflow-user-run-preset** | Presets run; gallery | B1 §2.2, §2.8 | **Ship.** Distinct zero-authoring goal ("run the demo, fetch results"). Output artifact = manifest + files. Fixed configs, so it stays clear of authoring/tuning. |
| 2 | **lbmflow-user-author-scenario** | Scheme *spec*, BC/property spec, mesh/lattice (2D/3D), **YELLOW obstacle composition**, config-validation | B1 §2.3, §2.4, §2.5, §3.1, §2.7 | **Ship (one Skill).** All produce the SAME artifact — a scenario JSON that passes `validate`. Folding obstacle composition in avoids two Skills fighting over `obstacles[]`. Config-validation is this Skill's done-check, not a separate goal. Owns the lattice-units disclaimer + routes unit conversion (RED) out. |
| 3 | **lbmflow-user-tune-stability** | Scheme *selection* (bgk\|trt), stability levers | B1 §2.3, §2.7 | **Ship.** Distinct goal: "which scheme / why did it diverge / clear these warnings". Operates on an existing scenario; splitting it from authoring keeps the author Skill about structure and this one about numbers. Gate = clean `validate`. |
| 4 | **lbmflow-user-run-monitor-mcp** | Async run/monitor via MCP | B1 §5 | **Ship.** Distinct subsystem (7 MCP tools, async loop, 4-cap). Separate from the blocking CLI runner. Gate = terminal `run_status`. |
| 5 | **lbmflow-user-postprocess** | Post/visualize (PNG/CSV/VTK/gallery, probes) | B1 §2.6, §2.8 | **Ship.** Distinct goal: declare + consume outputs. Owns the "3D PNG/CSV is a z-mid slice; VTK is the volume" nuance. |

**Folded, NOT separate Skills (justification):**
- *Config-validation* — folded into #2 and #3 as their verification gate. A
  standalone "validate" Skill would duplicate their done-criteria and overlap on
  `lbm validate` (the exact stacking hazard the brief warns against).
- *Obstacle composition* (B1 YELLOW §3.1) — folded into #2. B1 shows the
  primitives already execute (sphere ran, §2.5), so composition needs only a Skill
  wrapper, **no source change** → it is green-executable and belongs with the
  artifact it emits. (It is listed in §5 as a would-be yellow order ONLY as a
  contingency if PM rejects folding; the default is: no order needed.)
- *Blocking CLI `lbm run` of an authored scenario* — the mechanics live in the
  run-preset Skill's runner conventions and the MCP Skill's sync `run_scenario`;
  not worth a sixth Skill.

**Not shipped (RED — see §4):** GPU execution, CAD/mesh import, non-bgk/trt
schemes, 3D multiphase, quantitative validation-vs-reference, SI→lattice unit
conversion.

---

## 2. I/O contract per Skill

| Skill | Input (from user / upstream) | Output (artifact) | Verification gate (done-check) |
|---|---|---|---|
| run-preset | A preset name (one of 4) or "list/gallery" intent | A results dir: `manifest.json` + listed files; gallery `index.html` | Command line `status=completed` **and** `manifest.status=="completed"`. |
| author-scenario | NL flow description (+ optional preset template) in LATTICE units | A scenario JSON conforming to `references/schema.md` | `lbm validate` returns `ok:true` (warnings → route to tune-stability). |
| tune-stability | An existing scenario + a divergence/warning symptom | The scenario with adjusted `collision`/`nu`/`u`/resolution | `lbm validate` → `ok:true`, `warnings:[]` (or user-accepted). |
| run-monitor-mcp | A validated scenario OBJECT | `runId`(s) + retrieved manifest(s) | `run_status.state=="completed"` **and** `manifest.status=="completed"`. |
| postprocess | A completed run's `--out` dir (+ desired field/format) | Routed files (PNG/CSV/VTK/force/point) + reading guidance; or an `outputs`/`probes` block to add before running | Files exist in `manifest.files`, `manifest.status=="completed"`, format matches need (3D volume ⇒ VTK). |

---

## 3. Inventory table (ownership / triggers / non-overlap / do-NOT-use)

Skills deliberately touch the same repo surface (scenario JSON, `lbm validate`,
`manifest.json`), so each row fixes an ownership boundary so a smaller model can
route correctly.

| Skill | Ownership boundary (what it owns) | Trigger phrases | Explicit non-overlap rule | Do NOT use when… |
|---|---|---|---|---|
| **run-preset** | Preset discovery (`presets list/show`), running one preset (`presets run`), the whole-gallery pass (`gallery`), reading the resulting manifest/files. Fixed configs only. | "run a preset", "run the cavity/karman/droplet demo", "show me a simulation", "generate the gallery", "what presets exist" | Runs pre-tuned configs as-is. Any *edit* to the config → author-scenario. Does not author, tune, or use MCP. | The user wants a NEW or modified config (→ author-scenario / tune-stability) or async/sweep runs (→ run-monitor-mcp). |
| **author-scenario** | The scenario schema + accepted enums, 2D/3D choice, BC/property spec, obstacle composition from circle/rect/sphere, and the `validate` gate. Emits the JSON. | "set up a simulation", "make a scenario", "simulate flow past…", "add a cylinder/sphere/box", "write the scenario JSON" | Produces a validated scenario. Does NOT pick collision or fix stability warnings (→ tune-stability). Does NOT run it. Inputs are lattice units; routes unit conversion to the red note. | Running a fixed preset (→ run-preset), choosing scheme/fixing divergence (→ tune-stability), or running/monitoring (→ run-monitor-mcp / run-preset runner). |
| **tune-stability** | The two accepted collision options (bgk\|trt) and the three stability levers (tau floor, low-Mach cap, grid-Re). Adjusts numbers on an existing scenario. | "bgk or trt", "which collision", "why did it diverge / go NaN", "make it stable", "fix this tau/Mach/grid-Re warning", "what nu/velocity" | Tunes an already-drafted scenario's numeric knobs. Does NOT build structure/geometry (→ author-scenario) or run it. | Building a scenario from scratch (→ author-scenario) or running a fixed preset (→ run-preset). |
| **run-monitor-mcp** | The `lbm mcp` stdio server, its 7 tools, and the async loop (`start_run`→`run_status`→`list_runs`, 4-cap, result retrieval). | "run in the background", "kick off a long run", "run a sweep", "poll the run", "use the MCP server", "check run status", "list runs" | Owns async/parallel execution + monitoring. Does not author/validate the scenario (upstream) or post-process files (→ postprocess). | A quick blocking single run (→ run-preset runner) or authoring the JSON (→ author-scenario). |
| **postprocess** | Declaring `outputs`/`probes` so the run writes what's needed, and reading/routing PNG/CSV/VTK/force/point/gallery files (incl. the 3D-slice-vs-volume rule). | "visualize results", "get a PNG/CSV/VTK", "plot the drag over time", "open in ParaView", "make a gallery", "where are my files", "export the field" | Specifies + consumes outputs. Does NOT author physics (→ author-scenario) or run (→ run-preset / run-monitor-mcp). | Building the physics of the scenario (→ author-scenario) or executing it (→ run-preset / run-monitor-mcp). |

**Cross-Skill routing summary (the intended pipeline):**
author-scenario → tune-stability (if warned) → run-preset runner *or*
run-monitor-mcp (execute) → postprocess (consume). run-preset is the shortcut for
the fixed built-in configs.

---

## 4. Eval-harness CONTRACT per Skill (hard gates only — no held-out tasks)

Conventions (per A-pilot §5): **Baseline** = same prompt, no Skill. Held-out tasks
authored separately/adversarially. Each assertion is objectively checkable
(string / exit-code / file-presence / JSON-field). May run against a
fixture/dry-run harness asserting on the *commands/JSON the model would produce*
where live execution is costly. A **non-discriminating-assertion guard** names the
assertions that MUST fail the baseline to prove the Skill carries the expertise.

### 4.1 lbmflow-user-run-preset

| ID | Assertion | Threshold / artifact |
|---|---|---|
| RP-1 | Only one of the 4 real preset names is used (`cavity`, `cylinder-karman`, `two-phase-droplet`, `droplet-on-wall`) | 0 invented preset names |
| RP-2 | Uses `presets run <name>` (single) or `gallery` (all), matching the ask | correct command for the intent |
| RP-3 | Passes `--out <dir>` | present in every run command |
| RP-4 | "Done" asserted only on `manifest.status=="completed"`, not merely command exit | no success claim on `diverged`/`failed` |
| RP-5 | A config-edit request is routed to author-scenario, not run as a preset | 0 preset edits attempted |
| RP-6 | No scenario JSON authored inline (presets are fixed) | 0 hand-written scenarios |

Non-discriminating guard: **RP-1 and RP-4** must fail baseline (baseline invents
preset names / claims success on a diverged manifest).

### 4.2 lbmflow-user-author-scenario

| ID | Assertion | Threshold / artifact |
|---|---|---|
| AS-1 | Emitted JSON validates: `lbm validate` → `ok:true` (or model states the exact fix path) | pass on the fixture |
| AS-2 | Only accepted enum values used (`collision∈{bgk,trt}`, edges∈7-set, shapes∈{circle,rect,sphere}, field/format from the lists) | 0 out-of-schema tokens |
| AS-3 | 2D-vs-3D chosen correctly; 3D never emits multiphase / non-rest init / gpu | 0 build-rejected 3D features |
| AS-4 | Physical-unit input is NOT pasted as lattice values; conversion routed to the red note | 0 SI-as-lattice substitutions |
| AS-5 | Obstacle composition emits one `obstacles[]` element per requested shape | count matches the request |
| AS-6 | "Done" only after a `validate` pass (not on bare JSON emission) | no premature done |

Non-discriminating guard: **AS-2, AS-3, AS-4** must fail baseline (baseline
recommends MRT/gpu/3D-multiphase or treats SI numbers as lattice units).

### 4.3 lbmflow-user-tune-stability

| ID | Assertion | Threshold / artifact |
|---|---|---|
| TS-1 | Recommends only `bgk` or `trt` | 0 non-existent schemes |
| TS-2 | Fix is derived from the actual `validate` warnings for THIS scenario | model runs/reads validate before tuning |
| TS-3 | tau warning → `nu` raised so `tau≥0.55` (i.e. `nu≥0.0167`) | correct threshold applied |
| TS-4 | Mach: keeps `max|u|≤0.15` advisory, never above the hard 0.3 cap | 0 configs with `|u|>0.3` |
| TS-5 | grid-Re fix does not silently change the modeled Reynolds number when Re is the target (raises resolution instead of nu) | correct lever chosen per intent |
| TS-6 | "Done" only when `validate` → `warnings:[]` (or user-accepted) | no done with an open warning |

Non-discriminating guard: **TS-1, TS-3, TS-4** must fail baseline (baseline
suggests MRT, mis-sets the tau floor, or exceeds the Mach cap).

### 4.4 lbmflow-user-run-monitor-mcp

| ID | Assertion | Threshold / artifact |
|---|---|---|
| MC-1 | Only the 7 real tools are called (no `cancel_run`/`run_progress`/`stop`) | 0 nonexistent tool calls |
| MC-2 | Long/parallel work uses `start_run` (async), not blocking `run_scenario` | correct tool for the intent |
| MC-3 | ≤ 4 concurrent `start_run`; excess queued | max concurrent ≤ 4 |
| MC-4 | `state:"running"` is treated as in-progress, not stuck/failed | correct RUNNING verdict |
| MC-5 | "Done" requires terminal `state=="completed"` AND `manifest.status=="completed"` | no success on diverged manifest |
| MC-6 | Does not promise live %/step-progress/mid-run snapshot/cancel (not available) | 0 false-capability claims |

Non-discriminating guard: **MC-1, MC-3, MC-6** must fail baseline (baseline
invents a cancel/progress tool, over-fans past 4, or promises a live progress bar).

### 4.5 lbmflow-user-postprocess

| ID | Assertion | Threshold / artifact |
|---|---|---|
| PP-1 | Only accepted field/format used (`speed/ux/uy/rho/vorticity` × `png/csv/vtk`) | 0 out-of-list values |
| PP-2 | Outputs are declared in the scenario BEFORE running when the field wasn't written | adds `outputs`/`probes`, not a post-hoc recovery |
| PP-3 | 3D PNG/CSV described as a z-mid slice; true 3D volume routed to VTK | correct on a 3D task |
| PP-4 | Filenames sourced from `manifest.files`, not guessed | model reads the manifest |
| PP-5 | Force/point time series routed to `force.csv` / `point_<x>_<y>[_<z>].csv` for "plot over time" asks | correct file chosen |
| PP-6 | Outputs from a non-`completed` manifest are NOT presented as final | 0 diverged-as-final |

Non-discriminating guard: **PP-2, PP-3** must fail baseline (baseline tries to
recover an undeclared field, or presents a 3D PNG as the full volume).

### Scoring (all Skills, per A-pilot §5)

Pass rate = fraction of held-out tasks where ALL applicable hard gates pass.
Target: with-Skill ≥ 0.9 and strictly > baseline, with the guard assertions
failing baseline. Report mean ± stddev over 3 runs/task
(`aggregate_benchmark`), plus time/token deltas vs baseline.

---

## 5. Yellow-item implementation orders (codex-dispatchable)

Per the initiative rule, a Skill for a yellow item comes *after* an order lands.
The only B1 yellow item is **obstacle composition (§3.1)** — and B1's own "small"
check shows it needs **no source change** (primitives already execute). So the
**default recommendation is: NO order — fold it into author-scenario** (done; §1).

The order below is a **contingency** only if PM rejects folding and wants
composition as a first-class, tested capability with a smoke fixture. It is
written to be dispatchable as-is.

### ORDER Y1 (contingency) — obstacle-composition smoke fixture + doc

- **Scope:** Add a validated example scenario demonstrating multi-primitive
  obstacle composition (≥3 primitives incl. a 3D `sphere`) and a smoke test that
  `lbm validate` + a short `lbm run` both succeed on it. No engine/schema change —
  primitives already exist; this only locks the composition path with a test and a
  reference example the author-scenario Skill can point to.
- **Files (≤5):** a new example under `examples/` or `crates/lbm-cli/tests/`
  (composition fixture JSON), a smoke test in `crates/lbm-cli/tests/`, and a
  one-paragraph note in `docs/` referencing it. No changes to
  `crates/lbm-scenario/src/lib.rs`.
- **Acceptance:** (1) `lbm validate <fixture>` → `ok:true`; (2) a short
  `lbm run <fixture>` → `status:"completed"` with obstacle cells present (force
  probe non-trivial); (3) `cargo test --workspace --release` green; (4) the
  fixture uses ONLY `circle`/`rect`/`sphere` (no schema extension).
- **Effort:** ≤ 8h, ≤ 5 files, no new data model / numerical method / async
  subsystem — satisfies the "small" gate.
- **Dispatch note:** worktree-isolated codex order; verify with
  `lbmflow-build-verify` (Core tier + G8 CLI smoke) before merge.

**No other yellow orders** — B1 §3.2 confirms every other capability is either
already green (Skill wrapper only, all shipped in §1) or red (§4/§6, no Skill).

---

## 6. Red-item spec notes (research/spec only — NO Skill)

Each red item gets a short note so a Skill never tries to satisfy it, and so a
future order has a starting point. Negative evidence = B1 §4.

### R1 — SI→lattice unit conversion (B1 §4, surprise #4)

- **Status:** No conversion exists anywhere; all inputs are raw lattice units.
- **Skill routing:** author-scenario and tune-stability MUST NOT convert; they
  state "inputs are lattice units" and point physical-unit requests here.
- **Spec seed (for a future order):** a *prose/helper* non-dimensionalization
  guide — given a target Reynolds number and a chosen grid size + characteristic
  length in cells, derive a stable `nu` and inlet `u` in lattice units (respecting
  tau≥0.55, |u|≤0.15). This is dimensional analysis, not an engine feature; it
  could ship as a reference doc or a tiny helper script, NOT as a Skill that calls
  the engine (there is nothing to call). Adversarial risk: a naive Skill would
  silently paste SI numbers → wrong physics. Keep it a note until scoped.

### R2 — Quantitative validation vs. reference (B1 §4.5, surprise #3)

- **Status:** Analytic/reference checks (Ghia, Poiseuille/Couette/TGV, cylinder
  Cd/St, contact angle, conservation) exist ONLY as developer `cargo test` in
  `crates/lbm-core/tests/` — **no user-facing CLI/MCP compare command.**
- **Skill routing:** no user Skill claims accuracy-vs-reference. postprocess reads
  fields; it does not assert correctness against an analytic solution.
- **Spec seed:** an order could add a `lbm compare` subcommand (or MCP tool) that
  runs a scenario and reports L2/RMSE vs a named analytic case. This is a new
  subsystem (new CLI surface + reference data plumbing) → NOT small, NOT a Skill
  yet. Note only.

### R3 — GPU execution (B1 §4.4, surprise #2)

- **Status:** No GPU numerical backend. 3D rejects `backend:"gpu"` at build; 2D
  currently *silently falls back to CPU* (B1 surprise #2).
- **In-flight fix dependency (per brief):** explicit 2D `backend:"gpu"` will
  become a **validate-time error** instead of silent CPU fallback. **All B2 Skills
  are authored against that fixed behavior:** none recommend `gpu`; author-scenario
  §schema and tune-stability say gpu does not run and to use `cpu`/`auto`. When the
  fix lands, no Skill change is needed (they already never emit gpu). If the fix is
  NOT yet merged when these Skills ship, the only residual risk is a user manually
  writing `backend:"gpu"` in 2D and getting a silent CPU run — the Skills steer
  away from that but cannot suppress a hand-written value.

### R4 — CAD / mesh / STL import & unstructured/refined meshes (B1 §4.2)

- **Status:** No geometry-file import; lattice is uniform Cartesian only. Geometry
  = the three analytic primitives (covered by author-scenario).
- **Skill routing:** author-scenario offers ONLY primitive composition and says
  mesh/CAD import does not exist. Note only; a real importer is a large subsystem.

### R5 — 3D multiphase / 3D non-rest init (B1 §4.3)

- **Status:** Build-rejected. 3D is single-phase, `init:rest` only.
- **Skill routing:** author-scenario's Step-0 and schema §12 hard-steer 3D away
  from multiphase/non-rest init. Note only.

---

## 7. Deliverable manifest

- **This spec:** `docs/skills/b2-skill-specs.md`.
- **Shipped user Skills** (all validate clean via `python3 -m scripts.quick_validate`):
  - `.claude/skills/lbmflow-user-run-preset/SKILL.md`
  - `.claude/skills/lbmflow-user-author-scenario/SKILL.md` (+ `references/schema.md`)
  - `.claude/skills/lbmflow-user-tune-stability/SKILL.md`
  - `.claude/skills/lbmflow-user-run-monitor-mcp/SKILL.md`
  - `.claude/skills/lbmflow-user-postprocess/SKILL.md`
- **Yellow orders ready for dispatch:** Y1 (contingency only — default is fold,
  no order needed).
- **Red spec notes:** R1–R5 (no Skills).

## 8. Open questions for PM

1. **Fold vs. order for obstacle composition** — B2 folded it into author-scenario
   (no source change needed). Confirm, or dispatch contingency order Y1 for a
   first-class tested fixture.
2. **GPU fix timing** — Skills assume the 2D `backend:"gpu"` → validate-error fix
   lands. If it will not land soon, do we want a defensive line in author-scenario
   explicitly warning that a hand-written 2D `gpu` silently runs on CPU on current
   `main`? (Currently the Skills simply never emit gpu.)
3. **Unit-conversion note (R1)** — ship as a standalone reference doc / helper
   script next wave, or leave as prose routing inside author-scenario? B2 chose
   routing-only for now.
4. **run_scenario vs runner overlap** — blocking single-run mechanics are split
   between run-preset (runner conventions) and run-monitor-mcp (sync
   `run_scenario`). Acceptable, or promote a thin "run an authored scenario (CLI)"
   note? B2 judged a 6th Skill unwarranted.
