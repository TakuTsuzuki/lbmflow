# SPEC — Observer / function-object framework

> Research-agent deliverable (OpenFOAM's product moat: a composable, scheduled, restart-safe
> runtime + post-hoc measurement framework). Ready-to-dispatch, same depth/format as the four prior
> specs. Grounded in firsthand reads of `lbm-scenario/src/lib.rs` (`ProbeSpec`, `OutputSpec`,
> `FieldKind` L296–341), `lbm-cli/src/runner.rs` (`field_values`, `field_values_3d` L206–382), and
> `lbm-core/src/solver.rs` (`run_guarded` / `local_mass_partials` / `Diverged` L784–806).
>
> Maps COMPETITOR_ANALYSIS OpenFOAM §B3 (function objects) + OpenLB functor framework onto the
> JSON-scenario / CLI / MCP idiom. License: generic; no competitor code.

---

## 1. Scope & the generalization (additive)

Today measurement is split across two flat lists on the scenario: `outputs: Vec<OutputSpec>`
(a `FieldKind` → Png/Csv/Vtk every N steps) and `probes: Vec<ProbeSpec>` (`Force{every}` /
`Point{x,y,z,every}`). This covers snapshots + point series but not: time-averaging, line/plane
sampling, volume reductions, non-dimensional force coefficients, or convergence monitoring — the
checklist industrial users measure a CFD tool by (COMPETITOR OpenFOAM §B3).

Generalize both into **one scheduled `observers` list**, each observer = `kind` + `schedule` +
`region` + `sink`. **Additive & backward-compatible:** `OutputSpec`/`ProbeSpec` remain as
deserialization sugar that desugars into observers (a `FieldKind` snapshot and a Force/Point
observer), so every existing scenario and the 200+ suite stay green unmodified (R5). New scenarios
use `observers` directly.

---

## 2. The observer model

```jsonc
"observers": [
  { "id": "wake",     "kind": "snapshot",     "field": "vorticity", "sink": {"vtk": {}},
    "schedule": {"writeEvery": 500} },
  { "id": "drag",     "kind": "forceCoeffs",  "of": "cylinder",     "sink": {"csv": {}},
    "schedule": {"sampleEvery": 10} },
  { "id": "umean",    "kind": "fieldAverage", "field": "speed",     "stats": ["mean","rms"],
    "schedule": {"window": {"start": 20000, "sampleEvery": 5}, "writeEvery": 1000} },
  { "id": "profile",  "kind": "lineSample",   "field": "ux",
    "line": {"from": [40,0,0], "to": [40,128,0], "n": 128}, "sink": {"csv": {}},
    "schedule": {"writeEvery": 1000} },
  { "id": "converge", "kind": "residual",     "of": ["mass","speed"], "tol": 1e-5,
    "schedule": {"sampleEvery": 100} }
]
```

- **`id`** — stable key; names the output file(s), the checkpoint STATS blob (§5), and the MCP
  handle. Required and unique.
- **`kind`** — the catalog (§4).
- **`region`** — one of `field` (whole grid), `of` (a named obstacle/material group or a builtin
  signal), `line`/`plane`/`box` (sampling geometry), or a point. What's valid depends on `kind`.
- **`sink`** — `png`/`csv`/`vtk` (existing `OutputFormat`, extended) or `series` (in-memory /
  manifest time series, the current probe behaviour); MCP-readable regardless (§9).
- **`schedule`** — §3.

---

## 3. Scheduling model (executeControl vs writeControl)

Two independent clocks, mirroring OpenFOAM's execute/write split:

- **`sampleEvery` (execute clock)** — how often the observer *computes/accumulates* (a probe point,
  a statistic sample, a residual check). Default = every step for cheap observers, else required.
- **`writeEvery` (write clock)** — how often the observer *emits* to its sink. For snapshots this is
  the only clock (`OutputSpec::every` maps here; `0` = at end only, preserving today's semantics,
  lib.rs L328–330).
- **`window: {start, end?, sampleEvery}`** — for statistics: accumulate only within `[start, end]`
  (default `end = run end`) at `sampleEvery` stride. Decouples "average over the developed regime"
  from "write a snapshot of the running average every writeEvery."

A pure snapshot uses `writeEvery`; a probe uses `sampleEvery` (+ optional `writeEvery` to flush);
`fieldAverage` uses `window.sampleEvery` to accumulate and `writeEvery` to emit — the three
clocks are orthogonal and validated at load (e.g. `writeEvery` must be a multiple of the run's
output cadence is NOT required; any positive value is legal).

---

## 4. Observer catalog

| kind | consumes | emits | notes |
|---|---|---|---|
| `snapshot` | a `FieldKind` (via the single field-compute site, §6) | Png/Csv/Vtk | = today's `OutputSpec`; adds line/plane not — that's `lineSample`/`planeSample` |
| `probe` | point (x,y,z) | (ux,uy,uz,rho) series | = today's `ProbeSpec::Point` |
| `forceCoeffs` | `probed_force()` on a named group + the UnitConverter | Cd, Cl, (St from shedding freq) series | T8's `drag_lift` made user-facing (validation_cylinder.rs); **must read the active wall model's force**, §12-F2 |
| `fieldAverage` | a `FieldKind`, over `window` | mean and/or rms field → Vtk/Csv; restart-safe | accumulators = checkpoint STATS section, §5 |
| `lineSample` | a `FieldKind` along a poly-line (`from`,`to`,`n`) | Csv/Vtk 1D | linear interp between cells |
| `planeSample` | a `FieldKind` on an axis plane / cut | Csv/Vtk 2D | for a slice through a 3D field |
| `volReduction` | a `FieldKind` over `box`/`field`/material | min/max/sum/mean/integral scalar series | e.g. total enstrophy, mean speed |
| `residual` | `local_mass_partials()` + optional field-L2 delta | convergence status + value series | steady-state detection, §4.1 |

### 4.1 Residual / convergence monitor (cheap signal reuse)

Reuse `run_guarded`'s `local_mass_partials()` (solver.rs L802) — already computed for the divergence
guard — as the always-cheap mass residual. Optionally add a field-L2 steady-state delta:
`r = ‖φ(t) − φ(t−Δ)‖₂ / ‖φ(t)‖₂` for a chosen `FieldKind` sampled every `sampleEvery`. Emit a
machine-readable status `{value, tol, converged: bool}`; when `converged` and a scenario opts in
(`stopOnConverge: true`), the run may terminate early — directly serving agent-driven sweeps
(Pillar 1 / R4): an agent stops a converged case without a fixed step budget. `Diverged` (solver.rs
L181) remains the hard failure; `residual` is the soft convergence complement.

---

## 5. Restart-safe statistics (the checkpoint STATS section, made concrete)

`fieldAverage` (and `volReduction` time-integrals) keep **running accumulators**, not a buffer of
snapshots:

```
per observer id, per accumulated field:
  count : u64                    // samples taken
  sum   : [T; n_core]            // Σ φ            → mean = sum / count
  sumsq : [T; n_core]            // Σ φ²           → rms  = sqrt(sumsq/count − mean²)
```

**This IS the checkpoint `STATS` section reserved in SPEC_CHECKPOINT_RESTART §2/§3.** The blob
layout above is the concrete definition of that reserved section: keyed by observer `id`, versioned
under the same TLV table. On resume, the accumulators reload and the average continues seamlessly
across the restart — so a time-mean over `[20000, 100000]` is identical whether or not the run was
checkpointed at step 60000. (Numerical note: use the `sum`/`sumsq` form in the working precision; a
Welford variant is optional if f32 rms cancellation is measured to matter — flag as a follow-up, not
a v1 requirement.)

**Cross-reference resolved (finding §12-F1):** the checkpoint spec labelled STATS "reserved (M-F)".
`fieldAverage` is a P1/P2 observer that may ship **before** M-F, so STATS must be activated when the
first statistics observer lands — **not gated on M-F**. Correction filed as §12-F1.

---

## 6. Consuming FieldKinds (single compute site — do not recompute)

`FieldKind` today = `{Speed, Ux, Uy, Rho, Vorticity}` (lib.rs L314–320), computed by
`field_values` / `field_values_3d` (runner.rs L206–382). The in-flight ε-channel order is adding
**3D vorticity, Q-criterion (and λ2/enstrophy)** as `FieldKind`s. Per your directive, the observer
layer **consumes `FieldKind` values from one canonical compute site — it never re-derives them.**

Requirement (finding §12-F3): there must be exactly **one** `field_value(FieldKind, region)`
provider, consumed by `snapshot`, `lineSample`, `planeSample`, `fieldAverage`, and `volReduction`
alike. The ε-channel's Q/λ2/vorticity implementation is that site; the observer framework calls it.
Two vorticity implementations (one in ε-channel, one in observers) would be a defect — the spec
mandates the shared provider.

---

## 7. SI-valued output (via the resolved UnitConverter)

Observers that emit **physical** quantities (`forceCoeffs` Cd/Cl/St, any `sink` requesting SI units)
pull their conversion factors from the **resolved UnitConverter** of SPEC_UNIT_CONVERTER — they do
**not** re-derive them:

- `Cd = 2 F_x / (ρ U² A)`, `Cl = 2 F_y / (ρ U² A)`, `St = f_shed · L / U`, with `ρ, U, L, A` the SI
  characteristic values from the `units` block; `F` from `probed_force()` converted by
  `C_force = ρ·dx⁴/dt²` (SPEC_UNIT_CONVERTER §2).
- A CSV/manifest emitting SI velocity/pressure applies `C_velocity`, `C_pressure` from the same block.

**Cross-reference / placement (finding §12-F2):** because the core stays lattice-unit-only
(SPEC_UNIT_CONVERTER §6 invariant), the physical-value computation lives in the **scenario/runner
layer**, which owns the resolved UnitConverter and calls core for the raw lattice
`probed_force()`/fields. Observers emitting SI must therefore be a runner-layer concern (or be handed
the precomputed non-dim factors) — never inside `lbm-core`. Filed as §12-F2. If no `units` block is
present (raw-lattice scenario), `forceCoeffs` emits the lattice-unit coefficients and says so in the
output header.

---

## 8. `lbm postprocess --func` (post-hoc over saved artifacts)

OpenFOAM `foamPostProcess -func`: run observers over **already-saved** results without re-simulating.

- `lbm postprocess --func <observerId|kind> --from <ckpt_dir | field.vtk | run_dir>` applies an
  observer to a saved checkpoint (SPEC_CHECKPOINT_RESTART `ckpt_<step>/`) or a snapshot artifact.
- For fields present in the artifact it reads them directly; for a `FieldKind` derivable from the
  saved populations (a checkpoint carries `f` + moments) it recomputes via the §6 provider.
- Emits to the same sinks. This makes every observer **both** a runtime and a post-hoc tool — the
  OpenFOAM property that lets users add a measurement after a long run without repeating it.

Cross-reference: postprocess reading `ckpt_*/` binds this spec to SPEC_CHECKPOINT_RESTART §3 (it
consumes that layout); no contradiction — postprocess is a read-only consumer of the checkpoint.

---

## 9. MCP catalog + run_status (R4, agent-native)

- **Catalog discovery:** `lbm schema` and an MCP tool emit the full observer catalog — kinds, their
  params, and the currently-available `FieldKind`s — so an agent attaches diagnostics without
  hand-editing JSON and gets "valid choices" the way OpenFOAM's runTimeSelection errors do
  (COMPETITOR OpenFOAM §B6).
- **Live results:** the MCP `run_status` tool surfaces each observer's latest values (force
  coefficients, residual status, probe series tail) as structured JSON — same `diagnostics[]`/series
  shape the async job API already uses — so an agent can watch convergence and stop/branch a sweep
  (Pillar 1). `residual.converged` is the machine signal for early-stop.

---

## 10. Invariants & guardrails

- **Additive:** no `observers` ⇒ `outputs`/`probes` behave exactly as today (desugared); the 200+
  suite stays green (R5). `deny_unknown_fields` on the new structs (as `OutputSpec` already uses,
  lib.rs L323) keeps schema strictness.
- **No trajectory perturbation:** observers are read-only on simulation state; computing/emitting a
  measurement must not touch live buffers (same rule as checkpointing, SPEC_CHECKPOINT_RESTART §5).
- **Determinism (R4):** observer values are pure functions of the sampled state; MPI reductions
  accumulate rank partials in f64 like the existing diagnostics (dist.rs), so results are
  rank-count-independent to the diagnostic tolerance.
- **Single field-compute site** (§6) and **single UnitConverter** (§7) — no duplicated physics.

---

## 11. Acceptance / adversarial test matrix (for codex)

1. **Desugar equivalence:** a scenario using `outputs`/`probes` and the equivalent `observers` list
   produce byte-identical artifacts (proves additive back-compat).
2. **forceCoeffs vs T8:** the `forceCoeffs` observer reproduces `validation_cylinder.rs`'s
   `drag_lift` Cd/Cl within tolerance on the Schäfer-Turek case (same probed_force path, now
   user-facing) — and consumes the Bouzidi-corrected force when `wall: bouzidi` is set (§12-F2).
3. **fieldAverage correctness:** a steady Poiseuille `fieldAverage(mean)` equals the analytic
   profile; `rms` → 0 in steady flow and matches a reference in a known unsteady case.
4. **Restart-safe stats (ties to checkpoint):** a `fieldAverage` over `[N, N+M]` is **bit-identical**
   whether or not the run is checkpointed and resumed mid-window — the STATS blob round-trips
   (shared acceptance with SPEC_CHECKPOINT_RESTART §7).
5. **Line/plane sampling:** `lineSample` of `ux` across a channel matches direct cell reads (interp
   correctness); `planeSample` of a 3D field matches the corresponding VTK slice.
6. **Residual monitor:** `residual` reports decreasing mass/field deltas and flips `converged` when a
   steady case settles; a diverging case still trips `Diverged` (guard unchanged).
7. **postprocess parity:** `lbm postprocess --func` over a saved checkpoint yields the same observer
   output as the equivalent runtime observer (proves the §6 provider is shared runtime/post-hoc).
8. **MCP catalog:** the emitted observer catalog lists every kind + params; `run_status` returns
   live observer values in the structured shape (schema test).
9. **No-perturbation:** a run with observers attached has a bit-identical population trajectory to
   the same run without them.

---

## 12. Cross-spec bindings

- **STATS section** (SPEC_CHECKPOINT_RESTART reserved slot): activated by the
  first statistics observer (`fieldAverage`); blob layout in §5. Not gated on M-F.
- **SI-emitting observers** (`forceCoeffs`): live in runner/scenario layer using the
  UnitConverter (core stays lattice-only); must read the active wall model's
  force output (Bouzidi-consistent MEM when `wall: bouzidi`).
- **Single `field_value(FieldKind, …)` provider** — the ε-channel Q/λ2/vorticity
  site is that provider; observers consume, never re-derive.
