# PM Session Handoff — 2026-07-06

Handoff document for the next PM thread. Written by the outgoing PM (this Claude
Code session) covering everything a successor needs to pick up without
re-derivation. English throughout per owner directive; conversation with the
user may remain Japanese.

Author-session sessionId: `926ff0b4-2122-4cb4-8fd1-f9cd79d2786f`
Owner: Taku Tsuzuki (tsuzuki@epistra.jp).
Product: LBMFlow (commercial-grade LBM simulator, Rust core + TypeScript GUI +
Agent mode).

---

## 1. Standing owner directives (all binding)

These override any default behavior. Enforce mechanically.

1. **Product mission (Stop-hook /goal)**: "run continuously until it becomes the
   best conceivable LBM simulator." User sleeps/works elsewhere; PM decides
   everything within these directives.
2. **Full-power delegation with codex-max** (2026-07-05): the implementation
   workhorse is `codex exec`, fanned out aggressively (up to ~100 parallel;
   effective concurrency 5-8 compile-bound on M5 Max). One order = one bundle =
   one worktree = one codex. **`< /dev/null` is mandatory** or codex hangs.
3. **English only** (2026-07-05): all code, identifiers, commits, docs, GUI/CLI
   strings, error messages. Conversation with the user may remain Japanese.
4. **Physical rigor is the prime directive; ad-hoc physics is BANNED**
   (2026-07-06). Every physical behavior anywhere in the stack must be either
   resolved from the governing equations or a literature-backed closure with a
   recorded derivation, validity domain, and its own validation test
   (PHYSICS.md mandatory). Prohibited: constants calibrated to pass an
   acceptance band; branches keyed to sample/case identity; position clamps or
   caps that silently absorb transport; decorative physics terms. If a gate
   cannot be met without a hack, STOP and report. The executable procedure is
   `.claude/skills/lbmflow-physics-discipline`.
5. **Behavior-validity review** (2026-07-06): after every experiment/demo run,
   review whether the OBSERVED pattern (spatial, sign, trend) is physically
   plausible — a metric passing its band does NOT validate a pattern no band
   covers.
6. **Accuracy-first, turbulence-inclusive** (2026-07-06): turbulence-inclusive
   accuracy is the product's edge; comprehensively dig out any approximation
   that bends physics. This spawned the `lbmflow-accuracy-audit` Skill.
7. **Owner paper strategy overrides PM ruling B**: technical paper written in
   present tense as ideal state (paper = product spec). **RELEASE GATE** =
   paper stays internal draft until every present-tense claim is true; the
   `docs/paper/claims-ledger.md` is the gate instrument.
8. **Evidence-based progress**: before reporting anything as done, match each
   claim against a tool result from THIS session (test output, file diff, run
   log). Report unverified as unverified. "Done" = the Build & test gate
   passed here — a codex order finishing is not evidence its branch is green.
9. **Finish, don't announce**: never end a turn on a plan/checklist/"I'll now
   do X" — do X first, then end. Stop only when task complete or blocked on
   input only the user can provide.
10. **Minimal scope**: no drive-by refactors, no helpers for one-off ops, no
    abstractions for hypothetical futures. Validate at system boundaries only.

## 2. Team & session topology

Multiple long-running sessions collaborate on this repo. Communicate via
`mcp__ccd_session_mgmt__send_message` (cross-session messaging). Names below
appear in cross-session pings.

- **PM (this session)**: dispatches codex fleet, runs merge queue, arbitrates
  ownership splits, holds the claims ledger. `sessionId 926ff0b4-2122-...`.
- **QA sweep**: automated Physics Anomaly Sweep loop. Runs V&V matrices,
  behavior-reviews outputs, files findings, dispatches its own codex for
  MF-interim / W-LES / W-VOF / accuracy audits. Session
  `local_fbae3513-8550-41ee-a89c-bcb86a54991b`. Very active.
- **Skill-builder**: builds developer Skills (accuracy-audit, others).
  `local_a33c532c-b966-4c9c-8747-556e0edce493`.
- **Sales/paper**: writes the technical paper against the claims ledger.
- **Auto-approver** (background daemon): auto-approves `send_message`
  prompts. Location + restart in memory `cfd-auto-approver`.
- **Translation session** (completed): translated legacy Japanese to English.

Ownership deconfliction was negotiated with QA sweep on 2026-07-06 and is
recorded in cross-session messages — the split you inherit:
- **PM owns**: R-Phase 2 (B-1..B-8), M-E (ME-1..ME-4), W-GRAV proper,
  W-ROT proper, R2-C mechanical TRT port with ANOM-P2-001 fix.
- **QA sweep owns**: MF-interim (mf-grav/mf-rotor/mf-particles), scenario/runner
  integration, resuspension harness, W-LES, W-VOF, W-BCTOP, accuracy-audit
  dispatches, physics anomaly log.

## 3. Trunk state (as of handoff)

```
origin/main = 4eac49e qa: benchmark backlog (lane 4.2)
```

Local worktree (`/Users/taku/projects/流体シミュレータ`) is on branch **main**,
tip 4eac49e, matches origin. Working tree clean (only untracked `.agents/`,
`.codex/` metadata dirs).

Trunk contains all Wave 2 work: B-1 (Backend Fields generalization + run_span),
B-5 (TLV checkpoint/restart), ME-1 (D3Q19 WGSL + open-face BC + C-13 CLI +
FP16 scaffold — later resolved GREEN by QA-sweep quiet-window A/B), W-GRAV
proper, W-ROT proper (Uhlmann direct-forcing IBM), MF-interim (gravity, rotor,
particles), cx/wles (WALE LES with heavy-characterization freeze), cx/acc
(accuracy audit adversarial suite), qa/skill-accuracy-audit, MF-α stage 3
(central-moment collision on SIMD + GPU), D3Q27 track, VV master plan, cold
review triage, benchmark backlog. See `git log --oneline -50` for full history.

## 4. Claims ledger — release-gate status

From `docs/paper/claims-ledger.md`:

| Claim | Status |
|---|---|
| 3D GPU D3Q19 acceleration | **GREEN** (T14-3D + quiet-window 192³ 2791-2813, 128³ 2778-2880 MLUPS on unmodified main) |
| explicit `backend:"gpu"` | **GREEN** (C-13 landed 1a14d90) |
| FP16 storage ×2 grid capacity | **GREEN** (bands frozen; ~2.0× MLUPS @2048²; D3Q19 f16 >5 GLUPS) |
| Multi-node scaling ≥80% weak @64 rank | RED — needs cluster access confirm |
| Full-physics stirred workload | RED — longest lead |
| 2D GPU/CPU/T13/n≤4/wasm/MCP+Skills | GREEN |

The paper stays internal draft until Multi-node + Full-physics rows are true or
the paper is edited at release true-up.

## 5. Live codex processes and scratchpad orders

**No codex processes currently running** (`ps aux | grep codex exec` = 0).

Order files staged in
`/private/tmp/claude-501/-Users-taku-projects---------/926ff0b4-2122-4cb4-8fd1-f9cd79d2786f/scratchpad/`
— many are historical (already served). Ones a successor should evaluate before
re-using: `order-r2c.txt` (R2-C mechanical TRT port + ANOM-P2-001 single
forcing definition fix). Historical/served: order-b1*, order-me1*,
order-bouzidi*, order-g2, order-units, order-wgrav, order-wrot, order-b5,
order-d8, order-r2d, order-strain — all have landed content on trunk.

Codex logs `/tmp/codex-*.log` are transient (previous session tail-only).

## 6. Worktree inventory

`git worktree list` shows ~60 worktrees. Most are stale (branches long merged).
The primary checkout is on **main** at 4eac49e (path
`/Users/taku/projects/流体シミュレータ`). Do NOT prune indiscriminately — some
sessions may still hold references. Only prune worktrees whose branch is fully
merged AND whose parent session is confirmed done. When in doubt, leave them.

Notable worktrees a successor may need:
- `/Users/taku/projects/lbmflow-wt-r2-b1` — historic r2-b1 (B-1 landed via
  b1-land merge; branch superseded).
- `/Users/taku/projects/lbmflow-wt-me1{,-perf}` — ME-1 tracks (landed +
  perf-rejected; superseded).
- `.claude/worktrees/nostalgic-hellman-8b2924` [qa/anomaly-sweep] — QA sweep
  session's own worktree; DO NOT touch.
- `.claude/worktrees/wizardly-gauss-a20290` [qa/skill-accuracy-audit] —
  Skill-builder session's worktree; DO NOT touch.
- `.claude/worktrees/lucid-fermi-b7bd08` — translation session (done);
  `.claude/worktrees/youthful-germain-fd7ac4` — misc claude worktree.

## 7. Skills (developer + user)

`.claude/skills/`:

Developer Skills (agent-facing):
- **lbmflow-codex-dispatch** — order bundling, `codex exec` invocation,
  worktree assignment, safe concurrency, rollout-jsonl monitoring. Encodes the
  sandbox commit + Metal adapter traps.
- **lbmflow-build-verify** — exact cargo/wasm/web/CLI gate sequence before
  commit or merge.
- **lbmflow-physics-discipline** — mandatory for any change affecting computed
  physics. Provenance decision table, ban list with grep-able smells,
  two-layer gate template, stop-rule report, escalation table.
- **lbmflow-accuracy-audit** — bent-physics audit loop (5 attack axes A1-A5,
  P1-P4 phase/model-tier split, shared metrics library).
- **lbmflow-qa-viewer** — self-contained interactive 3D viewer of exported
  field volumes.

User Skills: `lbmflow-user-{author-scenario, postprocess, run-monitor-mcp,
run-preset, scaleup-advisor, tune-stability}`.

## 8. Merge-queue rules (hard requirements, encoded in accuracy-audit Skill)

Encode in every P4 fix order and landing gate:

(a) Landing gate = `cargo test --workspace --release --no-fail-fast`
    with **UNPIPED exit code** (`; echo EXIT:$?`). Never pipe a gate through
    tail/grep — pipe eats the exit code (bit twice this session; earlier
    backend_equiv regression was masked ~4h by fail-fast).
(b) GPU landing evidence is ALWAYS PM-run outside the codex sandbox. Adapter
    access in-sandbox is intermittent (2 of 5 orders had it in one wave).
    Codex GPU orders say "build + CPU gates only; report BENCH/GPU-PENDING".
(c) Audit/adversarial test orders NEVER share a worktree with implementation
    orders.
(d) Band freezes = measured value + stated headroom, measured values printed
    in assert messages; any loosening needs PHYSICS.md rationale + PM
    sign-off; the DENOMINATOR/normalization of every relative tolerance must
    be stated in the assert (per-component-relative vs scale-relative flipped
    pass/fail on identical physics this session).
(e) Current-wrong-value pins carry their ANOM id in the assert message; the
    fix order flips the pin to the correct value IN THE SAME COMMIT as the
    fix. Pin retighten never deferred.
(f) Anomaly-log entry (docs/qa/anomaly-log.md) required before merge for
    every P3-confirmed finding, including test-side dispositions.

## 9. Known traps (learned this window)

- **Backticks in inline codex order strings** die in zsh command substitution
  (strain-order incident). Pass orders via file:
  `codex exec ... "$(cat <order-file>)" < /dev/null > /tmp/codex-<tag>.log 2>&1 &`.
- **Sandbox `git commit` fails** intermittently with `index.lock` EPERM on the
  shared .git in worktrees (2/4 orders in one wave). Order text must include
  the "committed-ready fallback" clause — PM commits on codex's behalf at
  merge time.
- **Metal GPU adapter denial** in-sandbox is intermittent. `bench_gpu` /
  `--features gpu` tests need PM verification outside the sandbox.
- **`cargo test` fail-fast** masks regressions after the first failing binary.
  Landing gates MUST use `--no-fail-fast`.
- **Piping a gate through `| tail`** eats the exit code (reported green when
  red). Always separate: run gate raw, then `echo EXIT:$?`.
- **Keep-both merge resolution via naive regex on `<<<<<<< HEAD` /
  `=======` / `>>>>>>>` markers** can drop the closing `}` of the HEAD block
  when the two sides are two full functions concatenated. When you use the
  simple python one-liner, ALWAYS `cargo build --workspace --release` before
  committing the merge and hand-fix any unclosed delimiter (bit twice this
  session in solver.rs and scenario/lib.rs).
- **B-1's async run() contract**: `run(n)` on the GPU backend submits recorded
  chunks per C-9 calibration (auto-timed to ~100-250 ms) and returns after
  submission — completion is guaranteed at the next explicit `sync()` /
  `gather_*`. Tests that assert "run() blocks" are wrong; assert submissions
  counter + sync fence instead.
- **Cargo test binary ordering**: alphabetically-earlier tests run first —
  hiding later regressions when fail-fast is on. `t14_adversarial` runs
  before `t14_backend_equiv`; a bad startup-moment fix on the latter was
  masked all night by this.

## 10. Communication protocol with other sessions

Use `mcp__ccd_session_mgmt__send_message` with sessionId. The QA sweep and
Skill-builder sessions send narrative pings automatically at merge/finding
milestones; reply with:
- Ack + trunk tip they should rebase against.
- Ownership deconfliction (who owns W-ROT, W-LES etc.).
- Sequencing guidance ("hold your rebase until my B-1 merges").
- Rulings on tolerance/framing edge cases (denominator, headroom).

Recent binding rulings sent this session (all hold):
- 6 merge-queue rules to QA sweep (encoded in accuracy-audit Skill).
- Bouzidi as accuracy-audit Skill dry-run target.
- Metrics library placement: `crates/lbm-core/tests/common/metrics.rs` +
  `scripts/qa/metrics.py`, drift-guarded to 1e-12.
- WALE (not Smagorinsky) as W-LES default.
- ANOM-P2-001 folded into R2-C order.
- Probe-force tolerance denominator: force-vector L∞ scale (not per-component).

## 11. Immediate next moves

If picked up from an autonomous loop tick (no user prompt), the next moves are
QA-driven (QA sweep is the active driver of trunk motion this window):

1. Poll for new cross-session pings (they arrive automatically; you'll see
   them as `<cross-session-message>` on the next turn).
2. Check trunk (`git fetch && git log --oneline -10`) for QA-sweep merges you
   don't recognize — they land continuously.
3. If everything is quiet, do NOT proactively start new codex orders — the
   R-Phase 2 residuals (B-2, B-3, B-6..B-8, R2-C) are pending owner scheduling
   into the M-E/M-F integration campaign. R2-C is the highest-value next
   dispatch (ANOM-P2-001 physics fix bundled) if the owner or QA sweep asks
   for it.

If picked up with a user prompt, act on the prompt directly.

## 12. Files a successor should read first

1. `CLAUDE.md` — invariants, working discipline, physics rigor directive.
2. `docs/PLAN.md` — phase plan, current queue.
3. `docs/paper/claims-ledger.md` — release gate.
4. `docs/VALIDATION.md` — acceptance criteria.
5. `docs/PHYSICS.md` — physics decisions + experiment log.
6. `docs/qa/anomaly-log.md` — findings triage.
7. `.claude/skills/lbmflow-codex-dispatch/SKILL.md` — dispatch mechanics.
8. `.claude/skills/lbmflow-physics-discipline/SKILL.md` — physics ban list.
9. `.claude/skills/lbmflow-accuracy-audit/SKILL.md` — audit loop.
10. `TESTING_NOTES.md` — measured values, band freezes, bench data.

## 13. Memory pointers

Session-persistent memory index at
`/Users/taku/.claude/projects/-Users-taku-projects---------/memory/MEMORY.md`.
Key entries covering this project: `lbmflow-project-setup`,
`lbmflow-night-progress`, `lbmflow-solver-review-2026-07`,
`lbmflow-whitepaper-benchmark`, `lbmflow-research-agent-role`,
`lbmflow-qa-anomaly-sweep`, `lbmflow-dispersed-deposition-track`,
`lbmflow-physical-rigor-directive`.

Update `lbmflow-night-progress.md` after any large window of work.
