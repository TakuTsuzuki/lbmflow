# A-Pilot Held-Out Evaluation Protocol

This runbook evaluates the two Track A developer Skills from the frozen
held-out task set in `docs/skills/held-out/`. The evaluator must not reveal
held-out task text, held-out-unique fixture names, expected patches, or
held-out command output to Skill authors during repair.

## Scope

Skills under test:

- `lbmflow-codex-dispatch`
- `lbmflow-build-verify`

Task files:

- `docs/skills/held-out/cd-tasks.md`
- `docs/skills/held-out/bv-tasks.md`

Fixtures:

- `docs/skills/held-out/fixtures/rollouts/`
- `docs/skills/held-out/fixtures/repo-maps/`

Grader:

```bash
python3 docs/skills/held-out/grade.py <TASK_ID> <captured-transcript.txt>
```

## Frozen Models And Runs

Pinned models:

- Opus: `claude-opus-4-8`
- Sonnet: `claude-sonnet-5`

Run each task `N=3` per model. Use the same prompt text, fixture availability,
tool permissions, command timeout, and decoding settings for both models. If a
provider supports explicit decoding controls, record them; otherwise record
`defaults`.

A run succeeds only if every applicable hard gate for that task passes in
`grade.py`. A tool timeout, command failure, missing artifact, missing transcript,
or unverifiable output is a failed run. An infrastructure outage may be excluded
only if both models are rerun from scratch for the affected task set.

## Baseline Leg

Run a no-Skill baseline for every held-out task using the same model pool and
fixture access but without loading the candidate Skill.

Purpose:

- Verify the non-discriminating-assertion guard from the contract.
- CD-6, CD-7, and CD-8 must fail the baseline for dispatch tasks that exercise
  rollout forensics.
- BV-1, BV-3, and BV-5 must fail the baseline for build-verify tasks that
  exercise release-mode cargo tests, backend equivalence/T13, and the WASM
  `.gitignore` artifact gotcha.

If the baseline routinely passes these discriminating assertions, the held-out
task or grader is non-discriminating and must be revised before judging Sonnet.

## Per-Run Metadata

Each run record must include:

- Task id.
- Skill name and Skill version/path.
- Exact model id.
- Provider/version.
- Decoding settings or `defaults`.
- Tool permissions and sandbox/approval settings.
- Repo commit SHA.
- Command timeout.
- Machine/OS.
- Fixture bundle version or commit SHA.
- Captured transcript path.
- Grader JSON output.

Recommended transcript naming:

```text
eval-runs/<round>/<skill>/<model>/<task-id>/run-<1..3>.txt
eval-runs/<round>/<skill>/<model>/<task-id>/run-<1..3>.grade.json
```

## Scoring

For each Skill and model:

```text
pass_rate = successful_runs / (num_tasks * 3)
```

Acceptance:

1. Opus pass rate must be at least 80%. If Opus is below 80%, the Skill or eval
   is underspecified and must be revised before judging Sonnet.
2. Sonnet pass rate must be within 10 percentage points of Opus.
3. For every task/gate where Opus passes at least 2 of 3 runs, Sonnet must pass
   that same gate at least 2 of 3 runs.

Report per-Skill and per-model:

- Overall pass rate.
- Per-task pass counts.
- Per-gate pass counts.
- Baseline pass/fail summary for the discriminating assertions.
- Mean and standard deviation of time and tokens when available.

## Repair Protocol

If a frozen round fails, the evaluator may give generalized feedback only:

- Missing decision branch.
- Missing command or flag.
- Gate ambiguity.
- Artifact-path ambiguity.
- Uncovered failure mode.

The evaluator must not disclose:

- Held-out task text.
- Held-out-unique fixture filenames.
- Expected patch or exact expected answer.
- Held-out command output.
- Which hidden distractor caused the failure, beyond the generalized category.

After a Skill revision, start a fresh frozen round from scratch. Do not carry
partial successes forward. Record the old and new Skill versions/paths.

## Failure Classification

Use these categories in the round report:

- `skill-miss`: the transcript failed an unambiguous gate.
- `eval-ambiguous`: the task or grader could not objectively distinguish a
  correct answer from an incorrect one.
- `infra-outage`: provider/tooling failure that affected both models and was
  rerun from scratch.
- `fixture-error`: fixture missing, malformed, or inconsistent with the task
  spec.

`eval-ambiguous` and `fixture-error` findings require an eval revision and a new
frozen round before comparing Opus and Sonnet.

## GATE-AMBIGUITY Findings

The held-out grader freezes these assumptions because the contract references
some procedures without fully specifying their table entries:

- `BV-2` refers to a Step-0 gate-tier table that is not included verbatim in the
  contract. This eval reconstructs tiers from the hard gates, public examples,
  and repository instructions.
- `BV-3` explicitly names `backend_simd_equiv` and T13, while public example
  PE-B1 also expects full validation for pass-structure changes. Held-out tasks
  that alter backend pass structure may therefore require `BV-4` as well.
- `BV-5` says "when a wasm build occurs." This eval treats committed regenerated
  WASM package artifacts as enough to require the `.gitignore` removal and
  `git status --short` artifact check even if the prompt says no TypeScript
  source changed.
- `CD-5` is difficult to prove from a static transcript in dry-run mode. The
  grader accepts an explicit launch cap/queue plan rather than live process
  accounting.
