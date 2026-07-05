---
name: lbmflow-user-postprocess
description: >-
  Get LBMFlow simulation results out as images, data, or 3D volumes and know
  which file is which — PNG fields, CSV grids, VTK volumes (ParaView), probe time
  series, and the self-contained HTML gallery. Use whenever the user wants to
  "visualize the results", "get a PNG/CSV/VTK", "plot the drag/force over time",
  "open it in ParaView", "make a gallery", "export the velocity/vorticity field",
  "where are my output files", or asks how to read `manifest.json`, `force.csv`,
  or a step-numbered `.vtk` field file. This Skill owns declaring
  `outputs`/`probes` in a scenario
  so the run WRITES what you need, and reading/routing the files it produced. Do
  NOT use it to author the physics of the scenario (that is
  lbmflow-user-author-scenario) or to run it (run-preset / run-monitor-mcp) — this
  Skill is about specifying and consuming outputs, not computing them.
---

# LBMFlow — post-process & visualize outputs

A run only writes what its scenario asks for. This Skill has two jobs: (1) declare
the right `outputs`/`probes` so the simulation emits the fields and formats the
user wants, and (2) read back and route the files it produced (PNG/CSV/VTK,
`force.csv`, point series, `manifest.json`, the gallery `index.html`).

## Step 1 — Ask for the right output up front

Outputs must be declared in the scenario BEFORE the run — you cannot recover a
field that was never written. Add to the scenario:

```json
"outputs": [
  { "field": "speed",     "format": "png", "every": 0 },
  { "field": "vorticity", "format": "png", "every": 1000 },
  { "field": "rho",       "format": "vtk", "every": 0 }
],
"probes": [
  { "type": "force", "every": 100 },
  { "type": "point", "x": 128, "y": 40, "every": 50 }
]
```

Accepted values (recommend ONLY these):

| Slot | Accepted values |
|---|---|
| `outputs[].field` | `speed`, `ux`, `uy`, `rho`, `vorticity` |
| `outputs[].format` | `png`, `csv`, `vtk` |
| `outputs[].every` | N = snapshot every N steps; **`0` = only at the end** |
| `probes[].type` | `force` (obstacle force → `force.csv`), `point` (`x,y[,z],every` → point series) |

## Step 2 — Pick the format by what the user wants

| User wants… | Format | Why |
|---|---|---|
| A quick look / a figure to paste | `png` | Rendered field image. 2D = full field; **3D = z-mid slice**. |
| Raw numbers to plot/analyze yourself | `csv` | Grid of values. 2D = full field; **3D = z-mid slice**. |
| A real 3D volume in ParaView/VisIt | `vtk` | Legacy STRUCTURED_POINTS. 2D = `DIMENSIONS nx ny 1`; **3D = full volume `nx ny nz`**. |
| Drag/lift or a force time series | `probes: force` | `force.csv`: header `step,fx,fy` (2D) / `step,fx,fy,fz` (3D). |
| A value at one location over time | `probes: point` | `point_<x>_<y>[_<z>].csv` time series. |
| A browsable summary of ALL presets | `gallery` command | Self-contained `index.html` with embedded PNGs. |

**3D caveat to state clearly:** PNG and CSV in 3D are a **single z-mid slice**,
not the volume. For true 3D data the user must use `vtk` (the full volume) — do
not imply a 3D PNG shows the whole field.

## Step 3 — Find the files the run wrote

The run writes into its `--out <dir>` (or `outDir`). The authoritative list is
`manifest.json`'s `files` array — read it rather than guessing filenames:

```bash
cat <out>/manifest.json        # -> .files lists every output written
```

Filename conventions (from the runner):

| File | Pattern | Content |
|---|---|---|
| Field snapshot | `<field>_<step>.<png\|csv\|vtk>` | e.g. `speed_20000.png`, `rho_200.vtk` |
| Obstacle force | `force.csv` | `step,fx,fy[,fz]` time series |
| Point probe | `point_<x>_<y>[_<z>].csv` | single-point `(ux,uy,rho[,uz])` over time |
| Manifest | `manifest.json` | run status, diagnostics, and the `files` list |
| Gallery | `index.html` | self-contained; open in a browser directly |

## Step 4 — Route the file to the right tool

- **PNG** → show/attach directly; it is a finished raster image.
- **CSV** (field grid) → row-major `values[y*nx + x]`, one row per y. Load into
  pandas/numpy/a plotter. For 3D remember it is the z-mid slice.
- **VTK** → open in ParaView/VisIt (`File → Open`). It is ASCII legacy VTK
  STRUCTURED_POINTS; no conversion needed. This is the path for genuine 3D volume
  rendering / iso-surfaces / slicing.
- **force.csv / point CSV** → plot the column(s) vs `step` (e.g. drag `fx` vs
  step to see the Kármán shedding oscillation).
- **gallery `index.html`** → open in a browser; everything is embedded (no server,
  no external assets).

## Verification gate — the done check

Post-processing is done when:

1. `manifest.json` `status` is `completed` (a diverged/failed run's outputs are
   partial/meaningless — say so, don't present them as final), AND
2. every file you route to the user is actually listed in `manifest.files` and
   exists on disk (`ls <out>`), AND
3. the format you handed over matches the user's need per the Step-2 table (e.g.
   you did NOT hand a 3D PNG slice when they asked for the 3D volume — that needs
   VTK).

## Worked example (end-to-end)

Task: "I ran the cylinder case — give me the drag curve and a field image, and I
want the 3D version openable in ParaView."

1. **Declare outputs (Step 1):** for the drag curve add
   `probes:[{type:force,every:100}]`; for the image
   `outputs:[{field:vorticity,format:png,every:0}]`; for ParaView (3D) add
   `{field:speed,format:vtk,every:0}` and make it a 3D scenario (`nz>1`).
2. **Run** (via run-preset runner or MCP) → `manifest.status:"completed"`.
3. **Find files (Step 3):** `cat out/.../manifest.json` → `files` =
   `force.csv, vorticity_40000.png, speed_40000.vtk`.
4. **Route (Step 4):** plot `force.csv` `fx` vs `step` for the drag curve; attach
   `vorticity_40000.png`; tell them to open `speed_40000.vtk` in ParaView for the
   3D volume.
5. **Report** the three files and what each is.

## Top failure modes (and the fix)

- **Field/probe never declared, then "where's my CSV?"** Outputs are opt-in and
  must be in the scenario before the run. Fix: add the `outputs`/`probes` entry
  and re-run — data not written cannot be recovered.
- **Recommended a field/format that doesn't exist.** Only `speed/ux/uy/rho/
  vorticity` × `png/csv/vtk`. Fix: use the accepted lists.
- **Presented a 3D PNG/CSV as the whole field.** In 3D those are a z-mid slice.
  Fix: for the volume use `vtk`; say "slice" for PNG/CSV.
- **Guessed filenames.** Fix: read `manifest.files`; names follow
  `<field>_<step>.<ext>` / `force.csv` / `point_<x>_<y>[_<z>].csv`.
- **Presented outputs from a diverged run as final.** Fix: check
  `manifest.status == "completed"` first.
- **`every: 0` surprise.** `0` means "only at the end", not "every step". For
  time-lapse frames set a positive `every`.
