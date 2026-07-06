# P3 triage — decision table and logging format

A failing adversarial test is a *claim*. Calibration (2026-07-06): the first
six adversarial failures of the accuracy-audit program were ALL test-side
physics-design errors — zero engine bugs. Triage exists to keep that ratio
from generating six bogus fix orders.

## Step 1 — Derive before blaming the engine

Re-derive the reference independently of both the test and the engine.
Checklist of the recurring test-side traps (each has burned a real test):

| Trap | What the test got wrong |
|---|---|
| Wall placement | Half-way convention: the wall surface is midway between rim center and fluid center — analytic profiles must be evaluated at `y_wall = y_rim + 0.5`, not at the rim center. |
| Forcing half-step | `sim.ux()` etc. already include the Guo F/2 correction; adding it again (or subtracting it) shifts every force-driven reference by F/2. |
| Initialization transient | `init_with` equilibrium init injects a discrete transient; a reference assuming smooth startup needs either warmup steps or a startup-aware band. |
| Asymptotic regime | Order fits sampled outside the asymptotic range (too coarse, or error at float floor) give garbage slopes — check `r2` and the raw errors before reading the slope. |
| Compressibility pollution | Fixed u while refining h changes Ma per level; O(Ma²) error contaminates the order fit. Use diffusive scaling (u ∝ h). |
| Discrete time alignment | Step-count observables: is the observable after streaming or after collision? One half-step of ambiguity is first-order visible in single-step probes. |

Only when the independent derivation confirms the test does the failure
become a finding.

## Step 2 — Disposition

| Evidence pattern | Disposition | Action |
|---|---|---|
| Independent derivation contradicts the test's reference | **test-fix** | Fix the test yourself in the P2 worktree; log the trap (so the next audit's P1 avoids it); no order. |
| Derivation confirms test; error has a coherent mechanism in the engine (you can name the missing/wrong term) | **engine-fix order** | Log the anomaly; add a current-wrong-value PIN row to the test file (P2 convention 4); cut a P4 order. |
| Derivation confirms test; mechanism spans a design decision or an API that cannot express the correct behavior | **SPEC-GAP pin** | `#[ignore]` + `SPEC-GAP:` comment carrying the derivation; log it; raise to PM for the contract decision — do not cut a P4 order. |
| Test cannot be confirmed or refuted (no closed-form reference after honest effort) | **downgrade to cross-path (A5)** | Rewrite the probe as a consistency check between two paths; note the lost absolute anchor in the log. |

## Step 3 — Severity (house taxonomy)

| Severity | Definition | Example |
|---|---|---|
| **S0** | Silently-wrong physics: bounded, plausible output, no runtime signal | Out-of-envelope Ma 0.52 run reporting STABLE (pass-1 A1) |
| **S1** | Divergence/NaN leak in a supported configuration | grid-Re ~1160 case reaching NaN with no abort |
| **S2** | Below-expected accuracy or wrong transient, steady-state invisible | ANOM-P2-001 impulse deficit 4/7·F; rotor odd-blade mirror arms |
| **S3** | Minor: contract ambiguity, doc gap, cosmetic | "rotor chi=0: rejected vs no-op" contract gaps |

S0/S1 findings interrupt: message the PM immediately rather than batching to
the end of the audit. S2/S3 batch into the log.

## Step 4 — Log (mandatory, append-only)

Every finding — INCLUDING test-side errors — is appended to
`docs/qa/anomaly-log.md`:

```
**ANOM-<pass>-<nnn> — <one-line title>** — S<0-3>, disposition: <test-fix |
engine-fix order (owner/ref) | SPEC-GAP pin>.
- Scenario+config: <exact reproducer: file, sizes, tau, collision, steps>
- Expected: <value/curve WITH its source (derivation or citation)>
- Observed: <measured excerpt, verbatim from the test's println output>
- Impact: <what downstream measurement this bends>
```

Done-check: `grep` the log for every failing test name from the P2 report —
zero failing tests without a logged disposition.
