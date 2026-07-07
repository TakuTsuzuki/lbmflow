# V&V Visual Output and Anomaly Guide

Lifecycle: living — V&V visual-evidence requirements, updated in place.

This guide makes scalar-only V&V reports non-compliant. Any run that produces a
spatial field or spatial behavior must leave a visual artifact path next to its
scalar metrics, and the behavior review must cite that artifact before the run
is reported as evidence.

## Current Output Formats

The existing CLI and QA tools already provide the required raw material:

| Format | Current producer / parser | Use in V&V |
|---|---|---|
| `manifest.json` | `lbm run ... --out <dir>` / `scripts/qa/run_sweep.py` | Run status, diagnostics, warnings, generated file list, units. |
| PNG field images | CLI outputs with `format: "png"` | Human inspection of 2D fields and 3D mid-plane slices. |
| Field CSV | CLI outputs with `format: "csv"`; parsed by `scripts/qa/qa_checks.py` | Numeric postprocessing, profile checks, contour/slice plotting. For 3D, CLI CSV is a mid-plane slice. |
| Legacy ASCII VTK | CLI outputs with `format: "vtk"`; parsed by `scripts/qa/qa_checks.py` | Full 3D volume inspection in ParaView or qa-viewer-compatible tooling. |
| Probe CSV | `force.csv`, `point_<x>_<y>.csv`, `torque.csv`, `particles_<step>.csv` | Time-series plots and sign/trend checks. |
| Example artifacts | e.g. dispersed seeding `density.csv`, `metrics.json`, `*_velocity.vtk` | Track-specific evidence, still copied into the standard V&V layout. |

PNG/CSV/VTK are not interchangeable in claims. A 3D field claim needs a VTK
volume or explicitly named slices; a 2D pattern claim may use PNG plus raw CSV.

## Standard V&V Run Layout

Every V&V run pack SHOULD use this layout:

```text
out/vv/<run-id>/
  manifest.json
  scenario.json
  command.txt
  metrics.csv
  fields/
  plots/
  behavior_review.md
```

Required semantics:

- `manifest.json`: the exact CLI manifest, or a wrapper manifest with the CLI
  manifest embedded under `sourceManifest`. It must include status, steps,
  diagnostics, warnings, and all generated artifact paths.
- `scenario.json`: the exact scenario input. If the run is from a Rust test or
  example instead of the scenario CLI, include the equivalent test command and
  test parameters in `command.txt`.
- `command.txt`: one copy-runnable command per line, including environment
  variables that affect physics, backend, precision, or decomposition.
- `metrics.csv`: one row per sampled step. The first column is `step`. Include
  the scalar gates actually used for this run, such as `total_mass`,
  `mass_drift_rel`, `max_speed`, `max_mach`, `force_x`, `force_y`, `torque_z`,
  `l2_rel`, `linf_rel`, `cv`, or `spurious_u_max`.
- `fields/`: raw or directly copied field artifacts. Use names that preserve
  field, step, slice, and format, for example `speed_20000.png`,
  `rho_20000.csv`, `velocity_20000.vtk`, `midz_vorticity_20000.png`.
- `plots/`: derived visualizations of metrics or field reductions, preferably
  SVG/PNG. Every plotted quantity must be reproducible from `metrics.csv`,
  field CSV/VTK, or probe CSV.
- `behavior_review.md`: the behavior-validity review record. It is required
  before the run may be cited as validation evidence.

Compatibility note: existing `lbm run --out <dir>` output is still valid raw
evidence. For V&V campaigns, copy or symlink those files into `fields/` and write
the run pack files above rather than changing solver output code.

## Lightweight Plotting Utility

Use `scripts/vv/plot_metrics.py` to turn headered V&V metric/probe CSV files
into SVG plots without adding dependencies:

```bash
python3 scripts/vv/plot_metrics.py out/vv/<run-id>/metrics.csv \
  --out out/vv/<run-id>/plots
```

The script expects a header row, uses `step` as the x-axis when present, and
writes one SVG per numeric y-column. It is intended for `metrics.csv`,
`force.csv`, `torque.csv`, and point-probe CSVs, not raw field CSV files.

## Scalar-Plus-Visual Rule

A V&V run is reportable only when all applicable rows are satisfied:

| Run type | Minimum scalar evidence | Minimum visual artifact |
|---|---|---|
| Field-producing 2D run | Gate metrics in `metrics.csv` or test output | At least one PNG field image, plus raw CSV when a profile/symmetry claim is made. |
| Field-producing 3D run | Gate metrics in `metrics.csv` or test output | VTK volume, or named orthogonal slices if the claim is explicitly slice-local. |
| Probe-only run | Probe CSV plus derived scalar metric | Plot of the probe time series in `plots/`. |
| Particle/deposition run | Conservation and deposition metrics | Density map or VTK/tray field that shows the spatial distribution. |
| Backend/partition equivalence run | Difference norms and diagnostics | Difference field image/VTK when any nonzero spatial difference is observed. |

If a run emits only scalars while spatial behavior exists, mark it
`UNREVIEWABLE-SCALAR-ONLY` and do not use it as validation evidence.

## Behavior Review Template

```markdown
# Behavior Review: <run-id>

- Command: `<copy-runnable command>`
- Scenario: `scenario.json`
- Manifest: `manifest.json`
- Metrics: `metrics.csv`
- Visual artifacts reviewed: `fields/...`, `plots/...`

Pattern:
<what is observed spatially or temporally>

Mechanism:
<one sentence linking the pattern to resolved physics or a named validated closure>

Resolved vs closure:
<what comes from governing-equation physics, what comes from model closures or example-local rules>

Boundary / seam / outlet sweep:
<walls, clamps, partition seams, open boundaries, source/sink regions checked>

Verdict:
PHYSICAL | CLOSURE-DRIVEN | ARTIFACT | UNKNOWN | UNREVIEWABLE-SCALAR-ONLY

Routing:
none | core-engine defect | scenario/CLI defect | demo/example defect | spec revision | BENCH-PENDING
```

## Anomaly Templates

Each anomaly report must name a scalar trigger and a visual artifact. If either
is missing, the finding is incomplete.

### Symmetry Break

- Trigger: mirror/rotation difference exceeds the relevant band, or the visual
  field shows asymmetric extrema in a symmetric setup.
- Required artifacts: paired fields, difference field, scalar norm table.
- Inspect: wall-corner rules, moving-wall orientation, lattice direction maps,
  partition boundaries, asymmetric initialization.
- Report:

```markdown
## Symmetry break: <run-id>
- Expected symmetry: <mirror/rotation/translation>
- Scalar trigger: <metric and band>
- Visual artifacts: <paths>
- Locality: <bulk | wall | corner | seam | outlet>
- Physical mechanism found: <yes/no; describe>
- Disposition: core-engine defect | scenario defect | spec bug | needs more evidence
```

### Boundary Accumulation

- Trigger: mass, particles, scalar, density, or speed accumulates exactly at a
  wall, clamp, source boundary, or domain limit without a documented wall model.
- Required artifacts: field image/VTK covering the boundary, profile normal to
  the boundary, mass or particle ledger.
- Inspect: position clamps, sink/source guards, wall adhesion terms, deposition
  crossing logic, bounce-back/open-face adjacency.
- Disposition default: ARTIFACT until a wall model and validation test explain it.

### Outlet Reflection

- Trigger: pressure or velocity oscillations near an outlet exceed the accepted
  reflection band, or vortices visibly reflect upstream.
- Required artifacts: near-outlet field snapshots over time, central-region vs
  outlet-region pressure/speed RMS plot.
- Inspect: Outflow vs ConvectiveOutflow choice, reverse-flow fraction, pressure
  boundary staggered layer, stale-slot preservation.
- Note: the T4 pressure-outlet four-column staggered layer and T9 zero-gradient
  reflection are known artifacts with documented bands; do not call them
  validation failures unless they exceed those bands or appear outside their
  validity domain.

### Seam Artifact

- Trigger: discontinuity or extrema aligned with a partition, MPI rank boundary,
  GPU tile, slice boundary, or output stitching line.
- Required artifacts: field/difference image with seam overlay or coordinates,
  monolithic-vs-partition scalar comparison, manifest/decomposition metadata.
- Inspect: halo exchange scope, patch coordinate translation, source ownership,
  force/probe double counting, rank-local output stitching.
- Disposition default: core-engine or I/O defect unless reproduced in the
  monolithic reference.

### Spurious Current

- Trigger: nonzero velocity around a static interface, closed hydrostatic state,
  or equilibrium droplet exceeds the relevant band.
- Required artifacts: speed image, interface/density field, max-speed time
  series, location of max velocity relative to interface/walls.
- Inspect: force balance, pressure/EOS consistency, surface-tension closure,
  gravity composition, contact-angle wall term.
- Report whether the maximum is interface-local, wall-local, or bulk-local.

### Mass Leak

- Trigger: total mass, phase mass, scalar total, particle count, or gas volume
  drifts beyond its gate after accounting for explicit sources/sinks.
- Required artifacts: mass-ledger plot, source/sink schedule, density or phase
  field at first and last sampled steps.
- Inspect: open boundaries, source/sink ledger sign, density positivity guards,
  face patches, particle deposition/removal accounting.
- A scalar drift alone is not enough; include the spatial field that shows where
  mass is gained or lost when fields exist.

### Force / Torque Sign Error

- Trigger: force, lift, torque, or power sign contradicts the configured motion,
  symmetry, or reference benchmark.
- Required artifacts: force/torque time-series plot, velocity/vorticity field,
  geometry/probe mask image or coordinates.
- Inspect: momentum-exchange sign convention, reaction-vs-applied force naming,
  rotation direction, normal direction, reference-frame transforms.
- Report both conventions explicitly, for example "force on body" vs "force on
  fluid" or "reaction torque" vs "applied torque".

## Reporting Status Vocabulary

- `VALIDATED`: scalar band and behavior anchor pass, visual artifact reviewed,
  and validity domain matches the claim.
- `VERIFIED-ONLY`: implementation equivalence or regression passed, but no
  independent physical reference or behavior review exists.
- `SPEC-ONLY`: requirement exists, implementation or validation evidence does
  not.
- `BENCH-PENDING`: build/static evidence exists, but required hardware or heavy
  runtime evidence is unavailable.
- `UNSAFE-CLAIM`: current evidence does not support the claim or contradicts it.
- `UNREVIEWABLE-SCALAR-ONLY`: spatial behavior exists but no visual artifact was
  generated or reviewed.
