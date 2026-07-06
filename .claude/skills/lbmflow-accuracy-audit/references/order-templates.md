# P2 / P4 codex order templates

Fill the `{{...}}` slots and dispatch per lbmflow-codex-dispatch (one order =
one dedicated worktree; `codex exec --sandbox workspace-write
--skip-git-repo-check -C <worktree> "<order>" < /dev/null`, backgrounded;
monitor the rollout jsonl). Both templates are written to be runnable
verbatim by codex/Sonnet-class agents: every judgment call is made by the
audit list / triage disposition BEFORE the order is cut.

Known dispatch trap: codex may fail to `git commit` in shared-.git worktrees
(index.lock EPERM). If the rollout shows the work done but no commit, commit
for it yourself, **scoped to the ordered paths only** (`git add <paths from
the order>` — never `git add -A`).

---

## P2 ENCODE order template

```
Write adversarial accuracy-audit tests for {{subsystem}} in
{{test_file, e.g. crates/lbm-core/tests/accuracy_audit_bouzidi.rs}}.
Implement EXACTLY the audit list below - one test (or #[ignore] SPEC-GAP
stub) per row, no extra tests, no changes outside {{test_file}}.

AUDIT LIST (from P1; each row: approximation / analytic reference /
expected order-band / axis / cost):
{{audit_list_rows}}

HARD CONVENTIONS (all mandatory):
1. Derive every reference formula analytically in comments in the test file,
   from the physics, so the test reviews standalone. Do not cite the
   implementation as the reference.
2. Bands come from theory with ~10x headroom over float noise. Never tune a
   band to make a failing test pass - a failing test is a P3 finding, leave
   it failing and report it.
3. If the public API cannot express a probe, write the test body as the
   analytic derivation in comments, mark it #[ignore = "SPEC-GAP: ..."],
   and list it in your final report. Never silently drop a row.
4. Known-anomaly pins (rows marked PIN in the audit list) assert the CURRENT
   WRONG value stated in the row, with a comment naming the anomaly id and
   the correct value, and #[ignore = "expected failure until {{fix_ref}}
   fixes {{anomaly_id}}; this current-wrong-value pin must then fail loudly
   and be retightened"].
5. CPU reference backend only (compat Simulation<f64> or Solver with
   CpuScalar). No SIMD/GPU/MPI variants - T14 gates those transitively.
6. Use the shared metrics library: mod common; use common::metrics::*;
   (l2_rel, linf_rel, order_fit, envelope_fit, phase_fit, monotonicity,
   curve_agreement - see crates/lbm-core/tests/common/metrics.rs). Do NOT
   reimplement any metric inline.
7. Cost tags: rows marked "light" must run in <~1 s each; rows marked
   "heavy" get #[ignore = "heavy ACC-AUDIT {{description}}"]. Every heavy
   row with a light-canary note also gets the coarse light variant.
8. Each test prints its measured values (println! "ACC {{...}}: ...") so
   triage can read numbers without re-running.

VERIFY before reporting done:
- cargo test -p lbm-core --release --test {{test_file_stem}} compiles and
  runs; report the pass/fail status of EVERY test by name.
- cargo test -p lbm-core --release --test {{test_file_stem}} --
  --include-ignored --list shows the SPEC-GAP and heavy stubs.
- Failures are EXPECTED output at this stage: report them verbatim
  (assertion message + printed measured values). Do not fix the engine, do
  not loosen bands.

Commit to the current branch with message
"tests: accuracy audit for {{subsystem}} ({{n}} probes, {{m}} spec-gaps)".
```

### Slot-filling rules (for the dispatcher, not codex)

- `{{audit_list_rows}}`: paste P1 rows verbatim, including the derivation
  sketch — codex must not have to invent physics. A row whose derivation
  sketch you cannot write is not ready to dispatch.
- One P2 order per test file. Multiple subsystems never share a test file.
- The P2 worktree is a TEST worktree: it must never be reused for an
  implementation (P4) order.

---

## P4 FIX order template

```
Fix {{anomaly_id}} in {{owned_files, exact list}}.

FINDING (from docs/qa/anomaly-log.md - do not re-triage):
{{anomaly_log_entry_verbatim}}

ROOT CAUSE HYPOTHESIS (from P3):
{{mechanism, e.g. "per-cell force-field path applies the source term with
weight w_q only, missing the 1/(2 tau_minus) TRT correction that the
uniform path applies in ..."}}

SCOPE:
- Touch ONLY {{owned_files}}. If the fix genuinely needs another file, STOP
  and report instead of editing it (that file belongs to another order).
- No drive-by refactors, no new helpers beyond the fix.

ACCEPTANCE (all must pass, in this order):
1. RETIGHTEN THE PIN IN THE SAME COMMIT AS THE FIX: {{pin_test_name}} in
   {{pin_test_file}} currently asserts the wrong value {{wrong_value}} with
   an assert message containing "{{anomaly_id}}". Change it to assert the
   correct value {{correct_value}} with band {{band}} (denominator/
   normalization named explicitly in the assert message), remove its
   #[ignore], and rename it from *_current_wrong_value_pin_* to
   *_regression_pin_*. The pin test MUST remain in the file - it is the
   regression witness. Deleting it is an automatic reject. Deferring the
   pin retighten to a follow-up PR is also an automatic reject.
2. LANDING GATE (UNPIPED exit code, never pipe to tail/grep — the pipe eats
   the exit code and can report red as green):
     cargo test --workspace --release --no-fail-fast; echo EXIT:$?
   The gate passes iff EXIT:0 AND every test line is `ok`. --no-fail-fast
   is mandatory: without it, a later regression can be hidden behind an
   earlier failure.
3. GPU: this fix touches {{gpu_touched, y/n}}. If y, the order does BUILD +
   CPU gates only and reports "BENCH/GPU-PENDING for PM" - do NOT attempt
   GPU benches in the codex sandbox (Metal adapter access is intermittent
   in-sandbox; PM runs GPU evidence outside the sandbox).
4. PHYSICS: if this is a physics change (moments/BC/tau/forcing/collision):
   the full validation tier
     cargo test --workspace --release --no-fail-fast -- --include-ignored;
     echo EXIT:$?
   also passes, AND a docs/PHYSICS.md entry is added with rationale and
   measured before/after values.
5. If you find you must LOOSEN any band to make a test pass, STOP and
   report - loosening needs a PHYSICS.md rationale AND PM sign-off; do it
   yourself in the fix order and the order is rejected.

Commit with message "fix: {{anomaly_id}} {{one-line}}" and report the
before/after measured values printed by the pin test verbatim.
```

### Slot-filling rules

- `{{owned_files}}` comes from the PM/QA ownership map — confirm with the PM
  that no in-flight order is churning the same files before dispatch
  (collision-kernel and backend files are chronically contended).
- Cut one P4 order per anomaly unless two anomalies share a root cause in
  the same files — then one order fixes both and retightens both pins.
- P4 orders never edit the audit test file except for the pin retightening
  named in the acceptance — the test stays adversarial.
