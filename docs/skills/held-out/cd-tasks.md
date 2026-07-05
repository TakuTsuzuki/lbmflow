# Held-Out Tasks: lbmflow-codex-dispatch

These prompts evaluate the codex fan-out dispatch Skill only. They are
adversarial held-out tasks derived from `a-pilot-contract.md` and must not be
shown during Skill repair.

## CD-HO-01: Shared Guard File Plus Separate Test

Prompt text:

> Fan out these R-Phase guard items. A-2 edits `crates/lbm-core/src/guards.rs`,
> A-5 also edits `crates/lbm-core/src/guards.rs`, A-8 edits
> `crates/lbm-scenario/src/schema.rs`, and the adversarial A-8 test edits
> `crates/lbm-scenario/tests/schema_ranges.rs`. Use separate worktrees where
> appropriate, start the codex orders, and tell me how you will know each order
> is complete. Do not run the verification suite or merge anything.

Fixture requirements:

- `fixtures/repo-maps/dispatch-file-maps.md#map-cd-01`

Applicable gates:

- CD-1, CD-2, CD-3, CD-4, CD-5, CD-6, CD-9

Machine-checkable pass criteria:

- Transcript contains at least three `codex exec` dispatches, each with
  `< /dev/null` and `-C`.
- Transcript contains one combined order naming both `A-2` and `A-5`.
- Transcript contains separate orders for `A-8` and `T-A8` or the phrase
  `schema_ranges`.
- Transcript does not place the implementation item `A-8` and the test item
  `T-A8` in the same worktree path.
- Transcript mentions `task_complete` as the completion criterion.
- Transcript does not contain merge commands, `git merge`, or full gate-suite
  commands such as `cargo test --workspace --release`.

## CD-HO-02: Nine Compile-Bound Orders

Prompt text:

> Dispatch the M-F wave in parallel from the synthetic map: W-EXT, W-UNIT,
> W-STRESS, W-ROT, W-GRAV, W-SCAL, W-VOF, W-PART, and W-IO. They are all
> compile-bound except W-IO, and I care more about not melting the machine than
> shaving a minute. Show the launch shape and the exact command form. Leave
> verification and merge for later.

Fixture requirements:

- `fixtures/repo-maps/dispatch-file-maps.md#map-cd-02`

Applicable gates:

- CD-1, CD-2, CD-4, CD-5, CD-6, CD-9

Machine-checkable pass criteria:

- Every `codex exec` command contains `< /dev/null` and `-C`.
- Transcript states a maximum concurrent compile-bound launch count of `8` or
  fewer, and indicates that the ninth order waits or is queued.
- Distinct worktree paths are used for distinct orders.
- Transcript mentions `task_complete`.
- Transcript does not run cargo gate suites or merge commands.

## CD-HO-03: Running Log Must Not Be Called Stuck

Prompt text:

> The order in `/tmp/lbmflow-eval/wt-cx-mf-alpha` has no final answer yet. Check
> the fixture rollouts and tell me whether it is done, running, or stuck. Do not
> infer from quiet terminal output; use the rollout file.

Fixture requirements:

- `fixtures/rollouts/running-cx-mf-alpha.jsonl`

Applicable gates:

- CD-6, CD-7, CD-8, CD-9

Machine-checkable pass criteria:

- Transcript names `running-cx-mf-alpha.jsonl` or the cwd
  `/tmp/lbmflow-eval/wt-cx-mf-alpha`.
- Transcript references `session_meta.payload.cwd` or an equivalent cwd-based
  rollout match.
- Transcript states `RUNNING` or `running`.
- Transcript explains that no `task_complete` record is present.
- Transcript does not state `stuck`, `hung`, `complete`, or `done` as the final
  verdict.
- Transcript does not run verification or merge commands.

## CD-HO-04: Completed Log With Distractor

Prompt text:

> Two rollout files started at the same minute. I need the status for the order
> whose worktree is `/tmp/lbmflow-eval/wt-cx-rphase-a2-a5`, not the A-3 order.
> Pick the right log by cwd and report whether the order is complete.

Fixture requirements:

- `fixtures/rollouts/concurrent-a-correct.jsonl`
- `fixtures/rollouts/concurrent-b-distractor.jsonl`

Applicable gates:

- CD-6, CD-8, CD-9

Machine-checkable pass criteria:

- Transcript names `concurrent-a-correct.jsonl` or the cwd
  `/tmp/lbmflow-eval/wt-cx-rphase-a2-a5`.
- Transcript references `session_meta.payload.cwd` or equivalent cwd matching.
- Transcript states the order is complete because `task_complete` is present.
- Transcript does not select `concurrent-b-distractor.jsonl` as the target
  status file.
- Transcript does not run verification or merge commands.

## CD-HO-05: Stalled Log Is Not Complete

Prompt text:

> Check `/tmp/lbmflow-eval/wt-cx-d9`. The last order may have failed before a
> final response. Classify it from the rollout fixture and tell me what evidence
> is missing for completion.

Fixture requirements:

- `fixtures/rollouts/stalled-cx-d9.jsonl`

Applicable gates:

- CD-6, CD-7, CD-8

Machine-checkable pass criteria:

- Transcript names `stalled-cx-d9.jsonl` or the cwd
  `/tmp/lbmflow-eval/wt-cx-d9`.
- Transcript references cwd matching through `session_meta.payload.cwd`.
- Transcript states it is not complete because `task_complete` is absent.
- Transcript classifies the order as stalled, failed, or needs intervention,
  not `RUNNING`.
- Transcript cites the command failure or lack of new appended records as the
  reason for not calling it still running.

## CD-HO-06: Dispatch Only Despite User Asking "Make It Green"

Prompt text:

> Kick off codex for the GPU cleanup and the WASM smoke test from map-cd-03.
> When the orders finish, make sure they are green and ready to merge.

Fixture requirements:

- `fixtures/repo-maps/dispatch-file-maps.md#map-cd-03`

Applicable gates:

- CD-1, CD-2, CD-3, CD-4, CD-6, CD-9

Machine-checkable pass criteria:

- Transcript launches two separate `codex exec` orders with `-C` and
  `< /dev/null`.
- Transcript keeps GPU cleanup and WASM smoke test in separate worktrees.
- Transcript mentions `task_complete` for dispatch completion.
- Transcript explicitly routes verification to build-verify or says it will not
  run gate suites inline.
- Transcript does not contain `cargo test`, `wasm-pack build`, `npm run build`,
  `git merge`, or merge claims.

## CD-HO-07: Completed Single Order Is Done By Event Marker

Prompt text:

> Is `/tmp/lbmflow-eval/wt-cx-b7` finished? I only want the rollout-based
> answer, not a guess from the branch name.

Fixture requirements:

- `fixtures/rollouts/completed-cx-b7.jsonl`

Applicable gates:

- CD-6, CD-8

Machine-checkable pass criteria:

- Transcript names `completed-cx-b7.jsonl` or the cwd
  `/tmp/lbmflow-eval/wt-cx-b7`.
- Transcript references cwd matching through `session_meta.payload.cwd`.
- Transcript cites `event_msg` and `task_complete`, or at minimum
  `task_complete`, as the completion evidence.
- Transcript states the order is complete/done.

## CD-HO-08: Worktree Collision Trap

Prompt text:

> Start three codex orders: `docs English cleanup` in `docs/PLAN.md`,
> `perf table update` in `docs/PERFORMANCE.md`, and `GPU note update` in
> `docs/GPU_EVALUATION.md`. These are light, so you can use one worktree named
> `/tmp/lbmflow-eval/wt-docs` for all three, right? Show the exact commands.

Fixture requirements:

- None.

Applicable gates:

- CD-1, CD-2, CD-3, CD-5, CD-6, CD-9

Machine-checkable pass criteria:

- Transcript refuses or corrects the shared-worktree suggestion.
- Each `codex exec` command has a distinct `-C` path and `< /dev/null`.
- Transcript either dispatches three separate orders or explicitly bundles only
  if it describes a single combined docs order with one worktree; it must not
  show three independent orders sharing `/tmp/lbmflow-eval/wt-docs`.
- Transcript mentions `task_complete`.
- Transcript does not run verification or merge commands.

## CD-HO-09: Test Authoring Must Be Separate From Implementation

Prompt text:

> I want one codex order to implement `W-BCTOP` in
> `crates/lbm-core/src/boundary/topology.rs` and also write the adversarial T17
> topology tests in `crates/lbm-core/tests/t17_topology.rs`. It is faster if
> they share a branch. Dispatch it.

Fixture requirements:

- None.

Applicable gates:

- CD-1, CD-2, CD-4, CD-6, CD-9

Machine-checkable pass criteria:

- Transcript refuses to put implementation and adversarial test authoring in
  the same worktree/order.
- Transcript shows or describes two separate `codex exec` orders, each with a
  distinct `-C` path and `< /dev/null`.
- Transcript mentions `task_complete`.
- Transcript does not run gate suites or merge commands.

