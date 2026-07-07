---
name: lbmflow-codex-dispatch
description: >-
  Fan out LBMFlow implementation work to the codex CLI as parallel, worktree-
  isolated background orders, then monitor each order's rollout jsonl to
  completion. Use whenever the plan is to "dispatch to codex", "fan out",
  "parallelize", "run N codex orders", "kick off codex", or hand a bundle of
  implementation items to codex; also use when asked how to check on / monitor a
  running codex order or read its rollout log. This Skill owns order bundling,
  the exact `codex exec` invocation, worktree assignment, safe concurrency, and
  rollout-jsonl progress/completion/failure detection. Do NOT use it to actually
  build or run the test gates yourself (that is lbmflow-build-verify) or to merge
  landed branches (that is the PM merge-queue ritual). One order = one bundle =
  one dedicated worktree.
---

# LBMFlow codex fan-out dispatch

The LBMFlow implementation workhorse is the **codex CLI** (model gpt-5.5), fanned
out aggressively: many focused orders run in parallel, each in its own git
worktree, as background processes. This Skill carries the exact invocation, the
bundling rules, the concurrency ceiling, and — critically — how to read codex's
own session log (`rollout-*.jsonl`) to know whether an order is progressing,
done, or failed. Improvising any of these silently wastes hours (a wrong stdin
hangs codex forever; an over-bundled order serializes work that should parallel).

## The one invocation that always works (copy-runnable)

```bash
codex exec --sandbox workspace-write --skip-git-repo-check -C <worktree> "<order>" < /dev/null
```

- **`< /dev/null` is MANDATORY.** codex reads stdin; if stdin is a live pipe it
  waits for EOF and **hangs forever**. Redirecting from `/dev/null` gives an
  immediate EOF so codex runs the order and exits. This is the single most
  common way to lose a dispatch.
- **`-C <worktree>`** sets the working directory to the order's dedicated
  worktree. One order operates in exactly one worktree.
- Run it in the **background** so orders fan out concurrently (see Step 3).

## Step 1 — Bundle items into orders

The bundling rule follows one physical constraint: **git worktrees isolate by
directory, and two orders editing the same file will conflict.**

Decision procedure:

1. **Group candidate items by the file(s) they edit.**
2. **Items that touch the SAME file → bundle into ONE order.** Parallelize
   *across* files, never *within* a file. Two orders both rewriting `fields.rs`
   in separate worktrees produce conflicting branches that cannot both merge.
3. **Items on disjoint files → separate orders** (this is where the parallelism
   comes from).
4. **A test order and an implementation order NEVER share a worktree.**
   Validation tests are authored adversarially, kept separate from the
   implementation. Dispatch them as their own orders in their own worktrees.
5. One order should be a **focused item bundle** — coherent enough that a single
   codex run can finish and self-verify it, not a grab-bag of unrelated work.

## Step 1.5 — Physics-affecting orders: mandatory template clauses

If the order can alter any computed physical behavior (core, scenario,
examples, demos), its text MUST embed all four clauses — an order missing one
is malformed and gets Goodharted (2026-07-06 incident: closures calibrated to
CV bands produced an unreviewed edge-ring deposition pattern):

1. **Reading clause**: "Read and follow
   `.Codex/skills/lbmflow-physics-discipline/SKILL.md` — the provenance
   gate, ban list, and behavior-validity review are requirements of this
   order."
2. **Ban clause**: "No new constants without derivation + PHYSICS.md entry;
   no branches keyed to sample/case identity; no clamps/caps that absorb
   transport; no silent physical defaults."
3. **Stop-rule clause**: "If the gate cannot be met without violating the ban
   clause, STOP and emit the stop-rule report (a success outcome) — do not
   recalibrate constants to force green; that is a PM decision."
4. **Two-layer acceptance**: state the band AND at least one behavior anchor
   (sign / monotonicity / spatial structure) the result must satisfy, plus
   "attach the behavior-validity review record for every run you report,
   and list the visual artifact path (PNG / VTK / dashboard) for every run —
   scalar-only runs are unreportable; the PM does the looking".

The stop-rule clause has a proven track record: the P2 integration v1 order
stopped exactly right on a deposition-parity failure instead of recalibrating
the dispersion constant, which is what surfaced the real wiring bug.

## Step 2 — Assign one dedicated worktree per order

Each order gets its own worktree so the background runs never collide on the
working tree or the index.

```bash
# From the main repo (adjust branch/naming to your convention, e.g. cx-<topic>):
git worktree add ../lbmflow-wt-cx-<topic> -b cx/<topic>
```

Naming seen in practice: `lbmflow-wt-cx-<topic>` (e.g. `-cx-wasm-smoke`,
`-cx-d4`). Keep the worktree path stable — you pass it to `-C` and you locate its
rollout log by its `cwd` (Step 4).

## Step 3 — Fan out, respecting the concurrency ceiling

Concurrency is **machine-bounded**, not codex-bounded. `cargo` is CPU/RAM heavy.

| Order workload | Safe simultaneous count (this M5 Max) |
|---|---|
| Compile-bound (any order that runs `cargo build`/`cargo test`) | **~5–8** at once |
| Read/write-light (docs, small edits, no full compile) | more is fine; add on top |

Dispatch each order as a background process:

```bash
codex exec --sandbox workspace-write --skip-git-repo-check \
  -C ../lbmflow-wt-cx-<topic> "<order>" < /dev/null &
```

Do not exceed ~5–8 compile-bound orders concurrently — oversubscribing thrashes
RAM and makes *every* order slower, defeating the fan-out. Queue the rest.

## Step 4 — Monitor via the rollout jsonl

codex writes one session log per order at
`~/.codex/sessions/<YYYY>/<MM>/<DD>/rollout-*.jsonl`. Each line is one JSON
record. This is the authoritative progress/completion/failure signal.

**Locate the current order's log** (match by the worktree it runs in):

```bash
DAY=$(date +%Y/%m/%d)
WT=../lbmflow-wt-cx-<topic>              # the -C worktree you dispatched to
ABS=$(cd "$WT" && pwd)
# newest rollout whose session cwd == this worktree:
for f in $(ls -t ~/.codex/sessions/$DAY/rollout-*.jsonl); do
  head -1 "$f" | grep -q "\"cwd\": \"$ABS\"" && { echo "$f"; break; }
done
```

The first record is `type: "session_meta"` and its `payload.cwd` is the
worktree — that is how you disambiguate concurrent orders. (`turn_context`
records also carry `cwd`.)

**Progress / completion / failure — the record fields that matter:**

| You want to know | Look for | Field(s) that carry the signal |
|---|---|---|
| Did the order start? | `event_msg` / `task_started` | `payload.started_at`, `model_context_window` |
| Is it making progress? | `event_msg` / `agent_message` (narration) and `response_item` / `function_call` (name `exec_command` / `write_stdin`) | new records appended = alive; `agent_message.message` is human-readable status |
| Did an edit land? | `event_msg` / `patch_apply_end` | `payload.success` (bool), `payload.stdout` lists updated files |
| **Is it DONE?** | `event_msg` / `task_complete` | **terminal marker.** `payload.last_agent_message` = the final summary; `duration_ms`; `completed_at` |
| Token/context usage | `event_msg` / `token_count` | `payload.info.total_token_usage.total_tokens`, `model_context_window` |

**Terminal-state rule:** an order is finished **iff** a `task_complete` record
exists. `task_complete.last_agent_message` is codex's own end-of-run report —
read it to learn what changed, what was verified, and (crucially) which tests
codex left **red on purpose**. Absence of `task_complete` while records keep
being appended = still running (a long `cargo test --release -- --include-ignored`
can append `write_stdin`/`token_count` records for minutes — that is alive, not
stuck).

Full field reference and copy-runnable jq/python probes:
[references/rollout-schema.md](references/rollout-schema.md).

**Quick "is it done?" probe:**

```bash
python3 - "$f" <<'PY'
import json,sys
recs=[json.loads(l) for l in open(sys.argv[1]) if l.strip()]
done=[r for r in recs if r.get("type")=="event_msg" and r["payload"].get("type")=="task_complete"]
if done:
    p=done[-1]["payload"]
    print("DONE in %.1fs"%(p["duration_ms"]/1000))
    print(p["last_agent_message"][:1200])
else:
    last=recs[-1]["payload"]
    print("RUNNING — last record:", last.get("type"))
PY
```

## Step 5 — Hand off (do NOT do these here)

This Skill stops at "order is done". Two things happen next, owned by other
rituals — do not fold them in:

- **Verify the branch** with the build-and-verify gates → use
  `lbmflow-build-verify`. codex self-verifies, but the PM re-runs the gate tier
  before trusting a branch. Read `task_complete.last_agent_message` for which
  tests codex reports red.
- **Merge landed branches** in dependency order → the PM merge-queue ritual.

## Worked example (end-to-end)

Goal: dispatch R-Phase 1 items A-2…A-4 (all edit different files) plus a WASM
smoke test, in parallel.

1. **Bundle (Step 1):** A-2, A-3, A-4 touch disjoint files → 3 orders. The wasm
   smoke test is a *test* → its own order, own worktree (never shares with impl).
   → 4 orders total.
2. **Worktrees (Step 2):** `git worktree add ../lbmflow-wt-cx-a2 -b cx/a2`
   (…-a3, …-a4, …-wasm-smoke).
3. **Fan out (Step 3):** 4 orders, all compile-bound, under the 5–8 ceiling —
   dispatch all with the exact invocation, each `< /dev/null &`.
4. **Monitor (Step 4):** for each worktree, find its rollout by `cwd`, run the
   "is it done?" probe until `task_complete` appears. `agent_message` lines give
   live status; long `include-ignored` benches keep the log growing = alive.
5. **Hand off (Step 5):** for each `task_complete`, read `last_agent_message`,
   then verify with `lbmflow-build-verify` and queue for merge.

## Top failure modes (and the fix)

- **Forgot `< /dev/null`.** codex hangs with no output forever. Kill it and
  re-dispatch with the redirect. This is the #1 dispatch failure.
- **Two orders share a worktree (or edit the same file in separate worktrees).**
  Branches conflict and cannot both merge. Fix: re-bundle by file (Step 1); one
  worktree per order.
- **Over-fanned compile-bound orders (>8).** Everything slows via RAM thrash.
  Fix: cap at ~5–8 compile-bound; queue the rest.
- **Mistook a long bench for a hang.** No `task_complete` yet but `write_stdin` /
  `token_count` records keep appending → it is running a heavy
  `--include-ignored` suite. Do not kill it.
- **Read the wrong order's log.** With many concurrent orders, always match the
  rollout by `session_meta.payload.cwd == <worktree abs path>`, not by "newest".
- **Trusted `task_complete` as "all green".** It only means the run ended —
  `last_agent_message` may report intentionally red tests. Verify with the gate
  Skill before merging.
- **Sandbox cannot `git commit` in a worktree (`index.lock` EPERM on the shared
  `.git`).** Intermittent; 2/4 orders hit it on 2026-07-06. The order text must
  always include the committed-ready fallback clause; PM commits on codex's
  behalf at merge time. Do not treat it as a failed order.
- **Sandbox cannot acquire the Metal GPU adapter.** `bench_gpu`/GPU tests report
  `no usable GPU adapter was found` inside codex while the same binary works
  outside. GPU benches and GPU suites are ALWAYS PM-run landing evidence; orders
  must say "build only, report BENCH-PENDING (sandbox adapter)".
