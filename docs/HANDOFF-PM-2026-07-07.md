# PM Session Handoff — 2026-07-07 (refresh of 2026-07-06)

Handoff document for the next PM thread. Refreshed 2026-07-07 with actual
current trunk state (many QA-sweep-driven merges landed between the previous
handoff and this one). Reflects owner directive to remove paper-as-truth
framing from all operational docs. English throughout; conversation with the
user may remain Japanese.

Author-session sessionId: `926ff0b4-2122-4cb4-8fd1-f9cd79d2786f`
Owner: Taku Tsuzuki.
Product: LBMFlow (commercial-grade LBM simulator, Rust core + TypeScript GUI +
Agent mode).

---

## 1. Standing owner directives (all binding)

Enforce mechanically. These override any default behavior.

1. **Product mission**: "run continuously until it becomes the best conceivable
   LBM simulator" (Stop-hook /goal). User sleeps/works elsewhere; PM decides
   everything within these directives. The goal is the *simulator*, not any
   artifact about it.
2. **The paper is a living draft, NOT a goal** (2026-07-07 clarification,
   overrides the 2026-07-05 "paper = product spec" ruling). Do NOT hold
   implementation to a paper claim ahead of its measurement. Implementation
   converges on physics/spec targets (docs/PLAN.md, docs/VALIDATION.md,
   docs/PHYSICS.md). Update the paper when measurements change, not the other
   way around. `docs/paper/claims-ledger.md` is a working measurement-status
   snapshot, NOT a release gate — treat it as read-mostly context, not as a
   scheduler.
3. **Full-power delegation via codex-max** (2026-07-05): the implementation
   workhorse is `codex exec`, fanned out aggressively (~100 in principle;
   effective concurrency 5-8 compile-bound on M5 Max). One order = one bundle =
   one worktree = one codex. `< /dev/null` is mandatory or codex hangs.
4. **English only** (2026-07-05): all code, identifiers, commits, docs,
   GUI/CLI strings, error messages. Conversation with the user may remain
   Japanese.
5. **Physical rigor is the prime directive; ad-hoc physics is BANNED**
   (2026-07-06). Every physical behavior anywhere in the stack must be either
   resolved from the governing equations or a literature-backed closure with
   a recorded derivation, validity domain, and its own validation test (a
   PHYSICS.md entry is mandatory when physics changes). Prohibited: constants
   calibrated to pass a band; branches keyed to sample/case identity; position
   clamps or caps that silently absorb transport; decorative physics terms.
   If a gate cannot be met without a hack, STOP and report — the spec gets
   revised, not the physics faked. The executable procedure is
   `.claude/skills/lbmflow-physics-discipline`.
6. **Behavior-validity review** (2026-07-06): after every experiment/demo run,
   before reporting: review whether the OBSERVED pattern (spatial, sign,
   trend) is physically plausible. A metric passing its band does NOT
   validate a pattern no band covers.
7. **Accuracy-first, turbulence-inclusive** (2026-07-06): turbulence-inclusive
   accuracy is the product's edge; comprehensively dig out any approximation
   that bends physics. Spawned `lbmflow-accuracy-audit` Skill.
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
`mcp__ccd_session_mgmt__send_message` with sessionId. Names below appear in
cross-session pings.

- **PM (this session)**: dispatches codex fleet, runs merge queue, arbitrates
  ownership splits, holds the release-tracking view.
  sessionId `926ff0b4-2122-4cb4-8fd1-f9cd79d2786f`.
- **QA sweep**: the automated Physics Anomaly Sweep loop. Runs V&V matrices,
  behavior-reviews outputs, files findings, dispatches its own codex for
  MF-interim / W-LES / W-VOF / accuracy audits / D-track / MF-alpha /
  D3Q27 / VV master plan. Session
  `local_fbae3513-8550-41ee-a89c-bcb86a54991b`. Very active — most of the
  post-2026-07-06 landings are QA-sweep-driven.
- **Skill-builder**: builds developer Skills (accuracy-audit, others).
  Session `local_a33c532c-b966-4c9c-8747-556e0edce493`.
- **Auto-approver** (background daemon): auto-approves `send_message`
  prompts. Location + restart in memory `cfd-auto-approver`.
- **Translation session** (completed): translated legacy Japanese to English.

Ownership deconfliction negotiated with QA sweep on 2026-07-06 and applied
since:
- **PM owns**: R-Phase 2 (B-1..B-8) and M-E (ME-1..ME-4). W-GRAV proper and
  W-ROT proper were mine; both landed. R2-C (mechanical TRT port with
  ANOM-P2-001 fix) is staged in `scratchpad/order-r2c.txt` — still owned by
  PM if it dispatches.
- **QA sweep owns**: MF-interim, scenario/runner integration, resuspension,
  W-LES, W-VOF, W-BCTOP, accuracy-audit dispatches, physics anomaly log,
  D-track P0..P4, MF-alpha (D3Q27 + cumulant), T18 series, benchmark backlog.
  Since 2026-07-06 QA sweep has been the primary driver of trunk motion.

## 3. Trunk state — snapshot 2026-07-07

`origin/main = 62c964f` (contains the previous handoff commit; local main is
FF-caught up).

Recent post-2026-07-06 merges (all QA-sweep-driven unless noted):

- `4eac49e` qa: benchmark backlog (lane 4.2) — Kovasznay/Womersley/
  Sangani-Acrivos ready-now; MF-eta rigid-body proposal.
- `714da6a` qa: cold-review triage — 3 independent confirmations, F19/F20
  refuted (ANOM-P4-011), cumulant naming routed.
- `2e121c8` qa: band-vacuity scan (lane 2.2) — PM-triaged retighten queue.
- `fbdfd5b` qa: V&V master plan — 8 axes, wave-scheduled, sub-agent
  execution.
- `db50b24` rotor: document ANOM-P4-009 contracts (caller-must-clear, hub
  hole).
- `0caf4db` qa: ANOM-P4-009 (rotor audit test-side) + ANOM-P4-010 (S1
  penalization solid-disc divergence).
- `e569fb7` inventory: 3.1 cumulant offset reclassified (C) — confirmed by
  Order E measurement.
- `20d0e10` cx/cumulant-s3 merged — central-moment collision on SIMD + GPU
  (MF-alpha stage 3).
- `797ad14` cx/mock-engine-warning — runtime warn + MOCK badge.
- `3dde2bf` research: adhesion-capture + resuspension closure proposal
  (Phase B).
- `2b6510f` **discipline: physical-rigor prime directive + V&V routing +
  behavior-validity review** — this is where directives 5–7 above land in
  the operational docs.
- `97faaba` skill: **lbmflow-physics-discipline** — mechanical rigor
  procedure for developer agents.
- `5063452` skill: every experiment run must leave a visual artifact.
- `76b5071` docs: mark D-track P2 complete (CR-1/2/3 landed, example parity
  0.4%).
- Merges of cx/t18-{1,2,3}, cx/cr{1,2,3}-impl, cx/d3q27-{s1,val}, cx/cumulant-val,
  cx/gpu-wale, cx/dispersed-seeding — all part of the M-F / MF-alpha /
  D-track / M-E integration campaign.

M-E claims (measurement snapshot, `docs/paper/claims-ledger.md`):
- 3D GPU D3Q19: **GREEN**. Quiet-window MLUPS A/B/A: 192³ 2791-2813,
  128³ 2778-2880 (target ≥1,500). T14-3D GREEN.
- explicit `backend:"gpu"`: **GREEN** (C-13, commit 1a14d90).
- FP16 storage ×2 capacity: **GREEN**. Bands frozen (TGV transient 1.401e-1
  vs band 2e-1; cavity steady 2.579e-3 vs band 5e-3). ~2.0× MLUPS @2048²,
  D3Q19 f16 >5 GLUPS.
- Multi-node scaling ≥80% weak @64 rank: **RED** — needs cluster access
  (AWS hpc7g×8 ~¥13k on standby, or Fugaku trial; awaits owner spend
  confirm).
- Full-physics stirred workload: **RED** — longest lead, tracks M-F.

## 4. What is NOT done

The two RED rows above and the M-F integration wave. Concretely:

- **ME-3 cluster campaign** — bench_mpi 3D + weak modes ready in code; blocked
  on cluster spend confirm.
- **ME-4 full-physics benchmark** — waits for M-F integration completion
  (interface + LES + particles + scalar; a large fraction of the QA-sweep
  wave 3 work is inside this).
- **R2-C mechanical TRT port** (SPEC_COLLISION_COMPOSITION first PR) with
  ANOM-P2-001 (uniform-force vs force-field step-1 impulse mismatch) fix.
  Order staged in `scratchpad/order-r2c.txt`. Not yet dispatched — sequence
  it after any collision-adjacent QA-sweep landings settle.
- **Remaining R-Phase 2 items**: B-2 (transitional two-pass non-support
  capability), B-3 (Shan-Chen unification per SPEC_COLLISION_COMPOSITION),
  B-6..B-8 (per SOLVER_IMPROVEMENT_SPEC). Owner-scheduled into the M-E/M-F
  integration campaign; not blocking today.
- **W-VOF proper** (Allen-Cahn free surface, MF-γ, gated on W-GRAV proper
  landing which is DONE). QA sweep owns.
- **W-LES heavy characterizations** on the post-B-1 tree (Re_tau=180 vs DNS)
  — MKM 1999 reference profiles landed (7ae2baf); QA sweep planned to
  freeze after ME-1 stabilized (which it has).
- **Physics anomaly triage** on the pass-4 log — cadence QA-sweep-driven,
  new entries land continuously.

## 5. Live codex processes and scratchpad orders

**No codex processes currently running** (last check).

Order files in
`/private/tmp/claude-501/-Users-taku-projects---------/926ff0b4-2122-4cb4-8fd1-f9cd79d2786f/scratchpad/`
are historical (all named `order-*.txt` from earlier waves). The only one a
successor might dispatch as-is:
- `order-r2c.txt` — R2-C mechanical TRT port + ANOM-P2-001 fix. Read it,
  refresh against current trunk before dispatching.

All others (`order-b1*`, `order-me1*`, `order-bouzidi*`, `order-g2`,
`order-units`, `order-wgrav`, `order-wrot`, `order-b5`, `order-d8`,
`order-r2d`, `order-strain`) served their landings.

## 6. Worktree inventory

`git worktree list` shows ~60+ worktrees. Most are stale — branches long
merged. Do NOT prune indiscriminately; some sessions still hold references.
Prune only worktrees whose branch is fully merged AND whose parent session is
confirmed done.

Never touch these (owned by other sessions):
- `.claude/worktrees/nostalgic-hellman-8b2924` [qa/anomaly-sweep] — QA sweep.
- `.claude/worktrees/wizardly-gauss-a20290` [qa/skill-accuracy-audit] —
  Skill-builder.
- `.claude/worktrees/lucid-fermi-b7bd08` — translation session (done).
- `.claude/worktrees/youthful-germain-fd7ac4` — misc claude worktree.

## 7. Skills

`.claude/skills/`:

Developer Skills (agent-facing):
- **lbmflow-codex-dispatch** — order bundling, `codex exec` invocation,
  worktree assignment, safe concurrency, rollout-jsonl monitoring. Encodes
  the sandbox commit + Metal adapter traps.
- **lbmflow-build-verify** — exact cargo/wasm/web/CLI gate sequence before
  commit or merge.
- **lbmflow-physics-discipline** — mandatory for any change affecting
  computed physics. Provenance decision table, ban list with grep-able
  smells, two-layer gate template, stop-rule report, escalation table.
- **lbmflow-accuracy-audit** — bent-physics audit loop (5 attack axes A1–A5,
  P1–P4 phase/model-tier split, shared metrics library).
- **lbmflow-qa-viewer** — self-contained interactive 3D viewer of exported
  field volumes.

User Skills: `lbmflow-user-{author-scenario, postprocess, run-monitor-mcp,
run-preset, scaleup-advisor, tune-stability}`.

## 8. Merge-queue rules (hard requirements)

Encoded in the accuracy-audit Skill. Encode in every P4 fix order and landing
gate:

(a) Landing gate = `cargo test --workspace --release --no-fail-fast` with
    **UNPIPED exit code** (`; echo EXIT:$?`). Never pipe a gate through
    tail/grep — pipe eats the exit code. Both failure modes bit this session.
(b) GPU landing evidence is ALWAYS PM-run outside the codex sandbox.
    In-sandbox Metal adapter access is intermittent (2 of 5 orders had it
    tonight; also see the loaded-window MLUPS trap in
    `lbmflow-whitepaper-benchmark`). Codex GPU orders say "build + CPU gates
    only; report BENCH/GPU-PENDING".
(c) Audit/adversarial test orders NEVER share a worktree with implementation
    orders.
(d) Band freezes = measured value + stated headroom, measured values printed
    in assert messages; loosening needs PHYSICS.md rationale + PM sign-off;
    the DENOMINATOR/normalization of every relative tolerance must be stated
    in the assert (per-component vs scale-relative flipped pass/fail on
    identical physics this session).
(e) Current-wrong-value pins carry their ANOM id in the assert message; the
    fix order flips the pin IN THE SAME COMMIT as the fix. Pin retighten
    never deferred.
(f) Anomaly-log entry (docs/qa/anomaly-log.md) required before merge for
    every P3-confirmed finding, including test-side dispositions.

## 9. Known traps (learned this window)

- **Backticks in inline codex order strings** die in zsh command substitution.
  Pass orders via file: `codex exec ... "$(cat <order-file>)" < /dev/null >
  /tmp/codex-<tag>.log 2>&1 &`.
- **Sandbox `git commit` fails intermittently** with `index.lock` EPERM on
  the shared .git in worktrees. Order text must include the
  "committed-ready fallback" clause — PM commits on codex's behalf at merge
  time.
- **Metal GPU adapter denial** in-sandbox is intermittent. GPU tests / bench
  are PM-run.
- **Loaded-window MLUPS false-negative trap**: on unified memory, background
  cargo/codex load halves-to-thirds GPU MLUPS (D3Q19 192³: 1353 loaded vs
  2791-2813 quiet — cost a false ME-1 RED for hours). NEVER flip a perf gate
  RED (or dispatch kernel-opt orders) from a loaded-window number; re-measure
  quiet with A/B/A interleave.
- **`cargo test` fail-fast** masks regressions after the first failing
  binary. Landing gates MUST use `--no-fail-fast`.
- **Piping a gate through `| tail`** eats the exit code. Always separate:
  run gate raw, then `echo EXIT:$?`.
- **Keep-both merge resolution** via naive regex on `<<<<<<< HEAD` /
  `=======` / `>>>>>>>` markers can drop the closing `}` of the HEAD block
  when both sides are two full functions concatenated. Always
  `cargo build --workspace --release` before committing the merge and
  hand-fix any unclosed delimiter (bit twice this session in solver.rs and
  scenario/lib.rs).
- **B-1's async run() contract**: `run(n)` on the GPU backend submits
  recorded chunks per C-9 calibration and returns after submission;
  completion is guaranteed at the next explicit `sync()` / `gather_*`. Tests
  that assert "run() blocks" are wrong; assert submissions counter + sync
  fence instead.
- **Cargo test binary ordering**: alphabetically-earlier tests run first —
  hiding later regressions when fail-fast is on. `t14_adversarial` runs
  before `t14_backend_equiv`; a bad startup-moment fix on the latter was
  masked all night by this.

## 10. Communication with other sessions

Use `mcp__ccd_session_mgmt__send_message` with sessionId. QA sweep and
Skill-builder send narrative pings automatically at merge/finding milestones;
reply with:
- Ack + trunk tip they should rebase against.
- Ownership deconfliction.
- Sequencing guidance ("hold your rebase until my B-1 merges").
- Rulings on tolerance/framing edge cases (denominator, headroom).

Binding rulings sent this session (all hold):
- 6 merge-queue rules to QA sweep (encoded in accuracy-audit Skill).
- Bouzidi as accuracy-audit Skill dry-run target.
- Metrics library placement: `crates/lbm-core/tests/common/metrics.rs` +
  `scripts/qa/metrics.py`, drift-guarded to 1e-12.
- WALE (not Smagorinsky) as W-LES default.
- ANOM-P2-001 folded into R2-C order.
- Probe-force tolerance denominator: force-vector L∞ scale (not
  per-component).

## 11. Immediate next moves

If picked up from an autonomous loop tick (no user prompt), the next moves
are QA-driven since QA sweep is the active driver of trunk motion:

1. Poll for cross-session pings (`<cross-session-message>` on next turn).
2. `git fetch && git log --oneline -10` — QA-sweep merges continuously.
3. If everything is quiet, do NOT proactively start new codex orders. The
   R-Phase 2 residuals (B-2, B-3, B-6..B-8, R2-C) are owner-scheduled into
   the M-E/M-F integration campaign. R2-C is the highest-value discretionary
   dispatch (ANOM-P2-001 physics fix bundled) if the owner or QA sweep asks
   for it.

If picked up with a user prompt, act on the prompt directly.

## 12. Files a successor should read first

1. `CLAUDE.md` — invariants, working discipline, physics rigor directive.
2. `docs/PLAN.md` — phase plan, current queue.
3. `docs/VALIDATION.md` — acceptance criteria (T1..T18.x).
4. `docs/PHYSICS.md` — physics decisions + experiment log.
5. `docs/qa/anomaly-log.md` — findings triage.
6. `.claude/skills/lbmflow-codex-dispatch/SKILL.md` — dispatch mechanics.
7. `.claude/skills/lbmflow-physics-discipline/SKILL.md` — physics ban list.
8. `.claude/skills/lbmflow-accuracy-audit/SKILL.md` — audit loop.
9. `docs/paper/claims-ledger.md` — measurement-status snapshot (read-mostly).
10. `TESTING_NOTES.md` — measured values, band freezes, bench data.

## 13. Memory pointers

Session-persistent memory index at
`/Users/taku/.claude/projects/-Users-taku-projects---------/memory/MEMORY.md`.
Key entries: `lbmflow-project-setup`, `lbmflow-night-progress`,
`lbmflow-solver-review-2026-07`, `lbmflow-whitepaper-benchmark` (updated
2026-07-07 to reflect paper-as-living-draft framing),
`lbmflow-research-agent-role`, `lbmflow-qa-anomaly-sweep`,
`lbmflow-dispersed-deposition-track`, `lbmflow-physical-rigor-directive`.

Update `lbmflow-night-progress.md` after any large window of work.

## Appendix — what was reframed on 2026-07-07

The 2026-07-05 owner strategy that "the technical paper is the product spec
and implementation converges to the paper" was retracted 2026-07-07: it
misleads goals. The paper is now a living draft describing what has been
measured; the claims-ledger is a status snapshot, not a release gate;
implementation converges to physics/spec targets (PLAN.md, VALIDATION.md,
PHYSICS.md). Files edited to reflect this: `docs/paper/claims-ledger.md`,
`docs/paper/LBMFlow-whitepaper.md`, `docs/PLAN.md`, and the two memory files
(`lbmflow-whitepaper-benchmark`, `lbmflow-night-progress`, `MEMORY.md`).
