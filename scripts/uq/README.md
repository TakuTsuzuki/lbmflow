# Minimal UQ Sweep Harness

This directory contains a small input-sensitivity harness for the existing
`lbm` CLI. It runs parameter combinations, records raw quantities of interest
(QOIs), and computes one-at-a-time sensitivities plus bootstrap confidence
intervals where exact repeats exist.

This is an input-sensitivity harness, not a full ASME V&V 20 uncertainty
decomposition. It does not separate numerical error, model-form error,
discretization error, validation-data uncertainty, or calibration uncertainty.
Those decompositions are future work.

## Requirements

- Python 3.
- The already-built CLI binary at `target/release/lbm`, or pass another path
  with `--lbm`.
- `numpy` is optional. If available, bootstrap resampling uses it; otherwise
  the scripts fall back to the Python standard library.
- `matplotlib` is optional and used only by `analyze.py --plots`.

## Sweep Spec

The sweep spec is JSON. The basic form maps scenario parameter paths to a list
of values:

```json
{
  "parameters": {
    "physics.nu": [0.01, 0.02, 0.04],
    "grid.nx,grid.ny": [[33, 33], [49, 49]]
  },
  "repeats": 1
}
```

Path syntax is dot-separated JSON object access, with optional list indexes
such as `obstacles[0].r`. Comma-separated paths are a grouped parameter: each
listed value must be a list with the same length, and the paths are varied
together. The example above runs a 3 viscosity by 2 resolution sweep, not all
four combinations of `grid.nx` and `grid.ny`.

## Running Sweeps

From the repository root:

```bash
python3 scripts/uq/sweep.py \
  --preset cavity \
  --spec scripts/uq/example/sweep_spec.json \
  --out-csv scripts/uq/example/sweep.csv \
  --work-dir scripts/uq/example/runs \
  --lbm target/release/lbm
```

For a scenario JSON instead of a preset:

```bash
python3 scripts/uq/sweep.py \
  --base path/to/scenario.json \
  --spec path/to/sweep_spec.json \
  --out-csv out/uq/sweep.csv
```

The script reads presets through `lbm presets show <name>`, materializes a
modified scenario JSON for each combination, and runs `lbm run <scenario.json>
--out <run-dir> --json`. This keeps preset and JSON workflows on the same
runner path and allows preset parameters to be swept.

## QOI Definitions

The CSV is tidy: one row per run. Parameter columns are named `param.*`; QOI
columns are named `qoi.*`; execution metadata columns are named `run.*`.

Manifest QOIs come from `<out>/manifest.json`:

- `qoi.manifest.status`: manifest `status`.
- `qoi.manifest.stepsRun`: manifest `stepsRun`.
- `qoi.manifest.wallSeconds`: manifest `wallSeconds`.
- `qoi.manifest.mlups`: manifest `mlups`.
- `qoi.manifest.totalMass`: manifest `diagnostics.totalMass`.
- `qoi.manifest.maxSpeed`: manifest `diagnostics.maxSpeed`.
- `qoi.manifest.tau`: manifest `diagnostics.tau`.
- `qoi.manifest.warningCount`: length of manifest `warnings`.

Probe QOIs come from CSV files listed in manifest `files`:

- `force.csv` columns `fx`, `fy` become
  `qoi.probe.force.<component>.<last|mean|std|min|max>`.
- `point_X_Y.csv` columns `ux`, `uy`, `rho` become
  `qoi.probe.point_X_Y.<component>.<last|mean|std|min|max>`.
- Each probe also records `.sampleCount`.

Field snapshot QOIs come from CSV field outputs listed in manifest `files`:

- Scalar matrices such as `speed_200.csv` become
  `qoi.field.speed.<last|mean|std|min|max|rms|sum|cellCount>`.
- For `speed` fields, `qoi.field.speed.kineticEnergy` is computed as
  `0.5 * sum(speed^2)` over the exported field snapshot. This is a lattice-unit
  diagnostic over the exported cells, not a dimensional energy.

If a scenario writes only PNG or VTK outputs, no field-derived QOIs are computed
from those files. The harness does not infer hidden fields from images.

## Analysis

```bash
python3 scripts/uq/analyze.py \
  scripts/uq/example/sweep.csv \
  --out-dir scripts/uq/example
```

Outputs:

- `sensitivity.csv`: per-parameter one-at-a-time sensitivity. The normalized
  slope is `(dQOI/dparam) * mean(abs(param)) / mean(abs(QOI))`, averaged across
  groups where all other parameters are fixed.
- `bootstrap_ci.csv`: 95% bootstrap confidence intervals for the mean QOI at
  exact repeated parameter points with `n >= 2`.
- `summary.txt` and `summary.md`: plain-text and Markdown summaries.

The analysis does not drop outliers, smooth data, or recalibrate QOIs. Raw sweep
rows remain in the sweep CSV; statistics are clearly labeled derived outputs.
