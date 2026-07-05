# Track A — Skills for Developer: Pilot Report

Session: **A-pilot-author** (LBMFlow Skills Initiative, Track A).
Worktree: `/Users/taku/projects/lbmflow-wt-skills-a`, branch `skills/a-pilot`.
Date: 2026-07-05. All artifacts English. PM (Fable) owns gates and merges.

---

## 0. skill-creator discovery record (mandatory gate)

| Item | Finding |
|---|---|
| skill-creator available? | **Yes** — via the Skill tool (`anthropic-skills:skill-creator`). |
| Docs path | `/Users/taku/Library/Application Support/Claude/local-agent-mode-sessions/skills-plugin/ea9d5f86-0ac5-4abc-9c2b-8fccb400f67d/eaa2f29c-a524-404b-a25d-81cb50844cb4/skills/skill-creator/SKILL.md` |
| Scaffold command | **None** — skill-creator has no generator script. A Skill is a hand-authored directory: `<name>/SKILL.md` (YAML frontmatter `name`+`description` required) plus optional `references/`, `scripts/`, `assets/`. I created the dirs + `SKILL.md` by hand per the documented anatomy. |
| Skill output path (this repo) | `.claude/skills/<name>/SKILL.md` |
| Validate command | `python3 -m scripts.quick_validate <skill-path>` (run from the skill-creator dir). **Blocked locally**: the interpreter lacks `pyyaml` (`ModuleNotFoundError: No module named 'yaml'`). Fell back to a manual frontmatter check (valid `---` block, `name`==dir, `description` present, block scalar parses) — both Skills pass. PM should `pip install pyyaml` and re-run for the canonical check. |
| Eval / benchmark commands | `python3 -m scripts.aggregate_benchmark <workspace>/iteration-N --skill-name <name>` (benchmark); `eval-viewer/generate_review.py` (review viewer); `python3 -m scripts.run_loop --eval-set <json> --skill-path <path> --model <id> --max-iterations 5` (description-trigger optimization). Held-out eval authoring is explicitly out of my scope (Step 6). |
| Package command | `python3 -m scripts.package_skill <skill-folder>` → `.skill` bundle. |

No hand-rolled Skill format was used — the standard `SKILL.md` anatomy was followed.

---

## 1. Track A inventory

Each candidate below is a repeated LBMFlow dev operation currently living as prose
in CLAUDE.md / AGENTS.md / SOLVER_IMPROVEMENT_SPEC.md §4 / TESTING_NOTES.md. The
candidates deliberately overlap on tests/gates/docs, so each row states an
ownership boundary, trigger phrases, an explicit non-overlap rule, and a
"do NOT use when…" clause — so a smaller model cannot pick the wrong Skill or
stack conflicting done-criteria.

| # | Candidate | Ownership boundary (what it owns) | Trigger phrases | Explicit non-overlap rule | Do NOT use when… |
|---|---|---|---|---|---|
| A | **codex fan-out dispatch** | Bundling items→orders, the exact `codex exec … < /dev/null -C <wt>` invocation, worktree-per-order, ~5–8 concurrency ceiling, rollout-jsonl progress/done/fail detection. Stops at "order done". | "dispatch to codex", "fan out", "parallelize N orders", "kick off codex", "monitor the codex run", "read the rollout log" | Never runs the gate suite itself (→ B) and never merges branches (→ C). Hands the landed branch to B/C. | You are the one building/testing the tree (→ B), authoring the validation tests (→ D), or merging (→ C). |
| B | **build & verify ritual** | The exact cargo/wasm/web/CLI gate command set, the gate-tier decision table, and the per-gate observable done-check. | "build", "run the tests", "verify", "is it green?", "run validation", "make sure nothing broke", "before commit/merge" | Owns *running* gates on an already-written tree. Does not dispatch codex (→ A), does not decide merge order (→ C), does not write tests (→ D). | Dispatching parallel work (→ A) or authoring tests (→ D). |
| C | **PM merge queue** | Ordering landed codex/subagent branches by dependency, merging, reverting out-of-scope formatting diffs, keeping trunk green. | "merge the branches", "integrate the orders", "in what order do I land these", "PM integration" | Consumes branches produced by A and verified by B. Does not build (calls B) and does not dispatch (calls A). | Only one branch / no ordering decision, or the branch is not yet verified (→ B first). |
| D | **adversarial validation-test authoring** | Writing validation tests from VALIDATION.md **against** the implementation, kept in a separate worktree from impl; recording findings in TESTING_NOTES.md. | "write the T7 test", "adversarial test", "author validation from the spec", "acceptance test for TXX" | A *test* order is dispatched via A but never shares a worktree with an impl order. Runs gates via B. Distinct from C (no merging). | Implementing the feature itself, or merging (a test order and impl order never share a worktree). |
| E | **physics-change protocol** | The rule set for changing physics (moments/BC/tau/forcing): require full `--include-ignored` validation, record rationale+experiment in PHYSICS.md, keep CLAUDE.md/AGENTS.md invariants in sync. | "change the collision", "adjust tau/forcing/BC", "modify the physics", "update the moments" | Triggers B's G4 full-validation tier as its gate; adds the PHYSICS.md documentation obligation on top. Not a build Skill (delegates gates to B). | A pure build/refactor with no physics-spec change (→ B). |
| F | **invariant guardians** (D2Q9 ordering / q-major SoA / one-step pass structure) | The frozen core invariants and the specific gate that protects each: any pass-structure/storage-order change MUST pass `backend_simd_equiv` + T13 before landing. | "change fields.rs layout", "reorder D2Q9 dirs", "refactor collide/stream/step_band", "fuse the kernel" | Names the invariant + its required gate; the gate itself is *run* by B (this Skill supplies the "which gate is mandatory" judgment B's tier table encodes). | No core layout/pass-structure change is involved. |

---

## 2. Classification — Skill-ify vs leave-as-prose

| Candidate | Decision | Reason |
|---|---|---|
| **A codex fan-out dispatch** | **Skill-ify** ✅ (pilot) | Highest frequency; failure modes are silent and expensive (`< /dev/null` hang, worktree collisions, misread rollout). Carries a precise, copy-runnable procedure + a real jsonl schema — exactly what a Skill should own so the model needn't rediscover it. |
| **B build & verify ritual** | **Skill-ify** ✅ (pilot) | Highest frequency; the `--release` rule and the wasm-pack `.gitignore` gotcha are load-bearing and easy to forget. Objectively verifiable done-checks make it an ideal Skill. |
| **D adversarial test authoring** | **Skill-ify** (next wave) | Repeated and rule-heavy (separate-worktree discipline, TESTING_NOTES findings format), but by design authored by codex/Opus adversarially *separate from* implementation. Worth a Skill later; not a pilot to avoid coupling to Track A's dev-ops focus. |
| **E physics-change protocol** | **Skill-ify** (next wave) | Clear decision procedure (mandatory G4 + PHYSICS.md record). Lower frequency than A/B; layers on top of B, so best authored after B stabilizes to reuse its gate contract. |
| **C PM merge queue** | **Leave-as-prose (for now)** | Ownership sits with the PM (Fable), who owns gates/merges. A developer-facing Skill would blur that boundary. Re-evaluate if merge work is delegated. |
| **F invariant guardians** | **Leave-as-prose — folded into B** | The invariants are already encoded as B's "which gate is mandatory" tier rows (backend change → `backend_simd_equiv` + T13). A standalone Skill would duplicate B's done-criteria and risk conflicting gates for a smaller model — exactly the stacking hazard the brief warns against. The invariant *facts* stay in CLAUDE.md/AGENTS.md as the single source of truth. |

---

## 3. Pilot Skills authored

| Skill | Path | Version |
|---|---|---|
| **lbmflow-codex-dispatch** | `.claude/skills/lbmflow-codex-dispatch/SKILL.md` (+ `references/rollout-schema.md`) | v1 (initial) |
| **lbmflow-build-verify** | `.claude/skills/lbmflow-build-verify/SKILL.md` | v1 (initial) |

Both meet the five Sonnet-parity requirements:
1. **Explicit decision procedures** — dispatch has bundling rules + a concurrency
   table + a rollout-state table; build has a gate-tier table + a per-gate
   done-check table. No "use your judgment".
2. **Exact, copy-runnable commands** — the corrected command forms only
   (`< /dev/null` mandatory; wasm-pack `.gitignore` removal baked in as G6).
3. **Verification gate baked in** — dispatch names `task_complete` as the terminal
   artifact with a copy-runnable probe; build names each gate's observable green
   condition (exit codes, `test result: ok`, regenerated `pkg`).
4. **Worked example + failure modes** — each has one end-to-end example and a
   "top failure modes (and the fix)" section.
5. **Narrow scope** — dispatch stops at "order done"; build only checks an
   already-written tree. Cross-references route the rest to the correct Skill.

The rollout-jsonl documentation is sourced from **actual** files under
`~/.codex/sessions/2026/07/05/rollout-*.jsonl` (codex CLI 0.142.4): record types
(`session_meta` / `turn_context` / `event_msg` / `response_item`), the terminal
marker (`event_msg` → `task_complete`, carrying `last_agent_message` +
`duration_ms`), progress signals (`agent_message`, `patch_apply_end.success`,
`token_count`), and how to match a rollout to its worktree via
`session_meta.payload.cwd`.

---

## 4. Public examples (for the eval harness)

These are illustrative should-trigger tasks the pilot Skills must handle. They are
**public examples only** — the held-out evaluation tasks are authored separately
and adversarially by another session from the contract in §5. I did not write,
see, or run any held-out task.

### lbmflow-codex-dispatch

- **PE-A1**: "I've got R-Phase 1 items A-2, A-3, and A-5 ready plus the wasm smoke
  test. A-2 and A-5 both edit `guards.rs`. Fan these out to codex and tell me when
  each is done." → Expects: A-2+A-5 bundled into ONE order (same file); A-3 and
  the wasm test as separate orders (test never shares a worktree with impl); each
  dispatched with `codex exec … -C <wt> "…" < /dev/null &`; ≤ concurrency ceiling;
  monitored to `task_complete` per order.
- **PE-A2**: "codex order in `../lbmflow-wt-cx-d4` has been quiet for 3 minutes —
  is it stuck or still running?" → Expects: locate the rollout by `cwd`, inspect
  for `task_complete`; if absent but `write_stdin`/`token_count`/`function_call_output`
  records still append (e.g. a long `--include-ignored` bench), report RUNNING not
  stuck; only call it hung if no new records AND no `task_complete`.
- **PE-A3**: "Kick off a codex order to implement B-3 in its own worktree." →
  Expects: worktree created, exact invocation incl. `< /dev/null`, background,
  then rollout-based monitoring; hand-off to verify/merge noted, not performed.

### lbmflow-build-verify

- **PE-B1**: "A codex branch changed `simd.rs` step_band fusion — verify it before
  I merge." → Expects: tier = Core + G3(`backend_simd_equiv` + T13) + G4; all run
  with `--release`; done only when every required gate is green; report which
  gates ran.
- **PE-B2**: "I regenerated the WASM engine and the GUI. Make sure it's committable."
  → Expects: G5 wasm-pack build, **G6** remove `pkg/.gitignore` + confirm
  `git status --short` lists regenerated `pkg`, G7 `(cd web && npm run build)`.
- **PE-B3**: "Just tweaked a preset in lbm-cli — quick check nothing's broken." →
  Expects: Core (G1+G2, `--release`) + G8 `./target/release/lbm presets run cavity`
  exits 0 with step output; does NOT claim GPU/MPI coverage.

---

## 5. Eval-harness contract (CONTRACT ONLY)

How held-out tasks (authored separately) will be judged. Hard gates, artifacts,
and thresholds — **not** the tasks themselves.

### Shared conventions

- **Baseline** = same prompt, no Skill (new-Skill baseline per skill-creator).
- **Environment**: held-out tasks may run against a fixture/dry-run harness so no
  real codex order or ~5-min bench is required; graders may assert on the
  *commands the model would run* (transcript) rather than live execution. The
  rollout probes are checked against the real sample files in
  `~/.codex/sessions/2026/07/05/`.
- Each assertion is objectively checkable (string/exit-code/file-presence).

### lbmflow-codex-dispatch — hard gates (all must pass)

| ID | Assertion | Threshold / artifact |
|---|---|---|
| CD-1 | Every dispatch command includes `< /dev/null` | 100% of `codex exec` invocations in the transcript |
| CD-2 | Every order is dispatched with `-C <worktree>` and one worktree per order | no two orders share a worktree path |
| CD-3 | Same-file items are bundled into ONE order; disjoint-file items are separate | matches the task's file map exactly |
| CD-4 | A test order never shares a worktree with an implementation order | 0 violations |
| CD-5 | Compile-bound concurrency ≤ 8 simultaneous | max concurrent ≤ 8 |
| CD-6 | Completion is determined by presence of an `event_msg`/`task_complete` record | model cites `task_complete` (not "no output" heuristics) |
| CD-7 | A still-appending log with no `task_complete` is reported RUNNING, not stuck | correct RUNNING verdict on the fixture |
| CD-8 | Rollout is matched to its order by `session_meta.payload.cwd` | correct file selected among ≥2 concurrent logs |
| CD-9 | No merge/gate work performed inline (routed to build-verify / PM) | 0 gate-suite or merge commands issued by this Skill |

Non-discriminating-assertion guard: CD-6/CD-7/CD-8 must fail the baseline
(no-Skill) run to prove the Skill carries the expertise.

### lbmflow-build-verify — hard gates (all must pass)

| ID | Assertion | Threshold / artifact |
|---|---|---|
| BV-1 | Every `cargo test` invocation includes `--release` | 100% |
| BV-2 | Correct gate tier selected for the change's file set | matches the Step-0 table for the task |
| BV-3 | Backend pass-structure change runs G3 (`backend_simd_equiv` **and** T13) | both present when applicable; else N/A |
| BV-4 | Physics change runs G4 (`--include-ignored`) | present when applicable |
| BV-5 | After a WASM change, G6 runs: removes `pkg/.gitignore` + checks `git status --short` | both steps present when a wasm build occurs |
| BV-6 | GUI change runs `(cd web && npm run build)` | present when applicable |
| BV-7 | "Done" is declared only when every required gate's observable green condition holds | no premature "done" with an unrun/red required gate |
| BV-8 | Does not claim GPU/MPI coverage from the standard suite | 0 false coverage claims |
| BV-9 | No parallel-codex dispatch or merge performed inline | routed out |

Non-discriminating-assertion guard: BV-1, BV-3, and BV-5 must fail the baseline
run — they encode the exact gotchas the Skill exists to prevent.

### Scoring

- **Pass rate** per Skill = fraction of held-out tasks where ALL applicable hard
  gates pass. Target: with-Skill pass rate ≥ 0.9 and strictly > baseline, with
  the discriminating assertions (CD-6/7/8, BV-1/3/5) failing baseline.
- Report mean ± stddev over 3 runs/task (skill-creator `aggregate_benchmark`),
  plus time/token deltas vs baseline.

---

## 6. Handoffs / blockers

- **Blocker (minor)**: `quick_validate.py` needs `pyyaml`, absent in the local
  python3. Manual frontmatter validation passed for both Skills; PM should
  `pip install pyyaml` and re-run `python3 -m scripts.quick_validate` for the
  canonical check before merge.
- **Out of scope by design**: held-out eval tasks (authored adversarially by a
  separate codex session from §5), and the description-trigger optimization
  (`scripts/run_loop.py`) — run after PM accepts the pilots.
- **Next wave (recommended)**: Skill-ify D (adversarial test authoring) and E
  (physics-change protocol), reusing this build-verify gate contract.
