# Visual Anomaly Trawl Guide

Lifecycle: living — V&V master plan lane 5.3 harness guide, updated in place.

This guide documents `scripts/qa/visual_trawl.py`, the dependency-free Python
scanner for existing LBMFlow PNG and legacy ASCII VTK field outputs. It does not
run simulations. Operator execution of the scenario matrix is a separate task.

## Purpose

The trawl catches spatial failures that scalar gates can miss:

- velocity fields exceeding the low-Mach hard ceiling;
- density or mass fields drifting from the first checkpoint;
- PNG fields carrying a grid-scale checkerboard mode;
- scalar or magnitude fields violating a declared mirror/rotation symmetry.

The script prints a flagged file list and a per-file summary. A flagged file is
not automatically a physics defect; it is a behavior-review queue item for the
PM or V&V operator.

## Default Thresholds

The built-in profile is intentionally conservative:

| Threshold | Default | Justification |
|---|---:|---|
| `max_speed_hard` | `0.3` | LBMFlow's low-Mach hard limit. Prescribed speeds above this are invalid, and observed speeds above it are treated as an anomaly until explained. |
| `mass_drift_rel` | `1e-6` per checkpoint from step 0 | A postprocessing sentinel for conservation drift in output checkpoints. Tighter physics gates exist in validation tests; this lane is a broad trawl, so `1e-6` catches visible or accumulating drift without replacing test-specific bands. |
| `checkerboard_factor` | `0.5` | PNG checkerboard-mode power is flagged when it exceeds `0.5 * u0^2`, matching the init-ringing measurement rule used for lane-5 visual screening. |

`u0` defaults to `0.1` only as a fallback for ad hoc scans. Scenario campaigns
must provide `u0` per scenario type or via `--u0`; TGV-style diffusive setups
must not use the fallback silently.

The symmetry band defaults to `symmetry_rel = 1e-6`. It has no universal physics
meaning; it is a numerical trawl threshold for scenarios that explicitly declare
an exact visual symmetry.

## Basic Usage

Scan one output directory:

```bash
python3 scripts/qa/visual_trawl.py out/run-dir --scenario-type cavity --u0 0.1
```

Scan recursively and write a report:

```bash
python3 scripts/qa/visual_trawl.py out/qa/visual-trawl --recursive \
  --thresholds scripts/qa/visual-thresholds.json \
  --report out/qa/visual-trawl/anomaly-scan.txt
```

Emit machine-readable JSON:

```bash
python3 scripts/qa/visual_trawl.py out/run-dir --json \
  --report out/run-dir/anomaly-scan.json
```

Run only the synthetic harness check:

```bash
python3 scripts/qa/visual_trawl.py --self-test
```

Expected verification output:

```text
visual_trawl.py self-test PASS
```

## Gallery Mode

Gallery mode is for already-produced gallery output. It scans the current
gallery VTK naming convention:

```bash
python3 scripts/qa/visual_trawl.py --gallery gallery \
  --report out/qa/gallery-anomaly-scan.txt
```

This iterates:

```text
gallery/*/vtk_field_*.vtk
```

It does not generate the gallery and does not invoke `lbm gallery`.

## Threshold Profiles

Threshold magnitudes are scenario-specific. Use JSON profiles for campaign
runs:

```json
{
  "default": {
    "max_speed_hard": 0.3,
    "mass_drift_rel": 1e-6,
    "checkerboard_factor": 0.5,
    "u0": 0.1,
    "symmetry_rel": 1e-6
  },
  "tgv": {
    "u0": 0.02,
    "symmetry_rel": 1e-8
  },
  "droplet": {
    "u0": 0.005,
    "mass_drift_rel": 1e-7
  }
}
```

Run a named profile:

```bash
python3 scripts/qa/visual_trawl.py out/tgv --scenario-type tgv \
  --thresholds scripts/qa/visual-thresholds.json
```

Command-line flags override profile values:

```bash
python3 scripts/qa/visual_trawl.py out/tgv --scenario-type tgv \
  --thresholds scripts/qa/visual-thresholds.json --u0 0.015
```

## Scenario Hints

If a scenario JSON includes a QA block, the scanner can read it with
`--scenario`:

```json
{
  "qa": {
    "visualTrawl": {
      "scenarioType": "cavity",
      "symmetry": "mirror-x",
      "thresholds": {
        "u0": 0.1,
        "symmetry_rel": 1e-6
      }
    }
  }
}
```

Supported symmetry declarations are `mirror-x`, `mirror-y`, and `rot180`.
Symmetry is computed on scalar fields or velocity magnitudes; 3D volume symmetry
is not implemented in this first harness.

## Field Semantics

VTK:

- `VECTORS` arrays are reduced to velocity magnitude for `max|u|` and symmetry.
- scalar fields named like `speed`, `velocity`, `vel`, `u_mag`, or `umag` are
  treated as speed fields for the hard speed threshold.
- scalar fields named like `rho`, `density`, or `mass` are summed for mass drift
  relative to the first checkpoint in the same directory and field name.

PNG:

- The script decodes non-interlaced 8-bit or 16-bit PNGs using only the Python
  standard library.
- PNG luminance is mapped to field units with `--png-value-scale`; the default
  is `1.0`.
- Checkerboard mode is the normalized Nyquist FFT coefficient power,
  `|sum(field[x,y] * (-1)^(x+y)) / (nx*ny)|^2`.

Color-map PNGs are therefore screening artifacts, not quantitative validation
evidence, unless the campaign records the color-to-field scale used to generate
them.

## Operator Workflow

1. Produce the scenario/gallery outputs separately.
2. Run `visual_trawl.py` with a scenario-specific threshold profile.
3. Attach the report to the V&V run pack under `plots/` or the run root.
4. Behavior-review every flagged file using `docs/qa/VV_VISUAL_ANOMALY_GUIDE.md`.
5. Route confirmed defects to `docs/qa/anomaly-log.md`; do not treat the scan
   report alone as a final anomaly disposition.

As of 2026-07-07, only the harness and self-test are landed. No end-to-end
visual trawl campaign is claimed by this guide.
