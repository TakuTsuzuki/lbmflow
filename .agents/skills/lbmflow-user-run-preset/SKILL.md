---
name: lbmflow-user-run-preset
description: >-
  Run a built-in LBMFlow preset simulation from the CLI and fetch its results —
  the fastest way to get a working simulation and its output files without
  authoring any JSON. Use whenever the user wants to "run a preset", "run the
  cavity/cylinder/karman/droplet demo", "show me a simulation", "run the
  built-in example", "generate the gallery", or asks "what presets are
  available". This Skill owns preset discovery (`presets list`/`show`), running
  one preset (`presets run`), the whole-gallery pass (`gallery`), and reading the
  resulting `manifest.json` + output files. Do NOT use it to author a NEW
  scenario from a natural-language description (that is
  lbmflow-user-author-scenario) or to choose collision/stability settings (that
  is lbmflow-user-tune-stability) — presets are fixed, pre-tuned configs you run
  as-is.
---

# LBMFlow — run a built-in preset

LBMFlow ships four pre-tuned, physics-validated presets. Running one is the
zero-authoring path to a real simulation: no scenario JSON, no unit choices, no
stability tuning — the preset already encodes all of that. This Skill covers
discovering presets, running one (or all of them as a gallery), and reading back
the results the run wrote to disk.

**Prerequisite:** the `lbm` binary must exist at `./target/release/lbm`. If it is
missing, build it once with `cargo build --workspace --release` (release is
required — LBM is ~50x slower in debug). Everything below is run from the repo
root.

## The four presets (exact names — nothing else is a preset)

| Name | What it is | Rough wall time |
|---|---|---|
| `cavity` | Lid-driven cavity, steady-state detection | ~60–70 s |
| `cylinder-karman` | Kármán vortex street past a cylinder + drag probe | ~160 s |
| `two-phase-droplet` | Shan-Chen two-phase droplet equilibration (2D) | ~55 s |
| `droplet-on-wall` | Contact-angle demo on a wall (`wallRho=1.0`, θ≈63°) | (completes) |

If the user names anything else ("run the airfoil preset"), it does not exist —
say so and offer the author-scenario Skill instead. Do not invent preset names.

## Decision procedure — which command

| User intent | Command |
|---|---|
| "What presets are there?" / list them | `presets list` |
| "Show me the cavity config / JSON" | `presets show cavity` |
| Run ONE named preset | `presets run <name> --out <dir>` |
| Run ALL presets + build an HTML gallery | `gallery --out <dir>` |

## The commands (copy-runnable, exact)

Run from the repo root. Always pass `--out <dir>` so you know where results land.

```bash
# List the available presets (name + one-line description)
./target/release/lbm presets list

# Inspect one preset's full scenario JSON (useful as an authoring template)
./target/release/lbm presets show cavity

# Run ONE preset; results go to the given directory
./target/release/lbm presets run cavity --out out/cavity

# Run ALL presets and build a self-contained HTML gallery (~5 min total)
./target/release/lbm gallery --out out/gallery
```

## Verification gate — the done check

A preset run is done only when BOTH hold:

1. The command's final line reports `status=completed` (NOT `diverged` or
   `failed`). Example of a good line:
   `status=completed steps=20000 wall=69.3s mlups=5 out=.../out/cavity`
2. `<out>/manifest.json` exists and its `"status"` field is `"completed"`.

Read the manifest to fetch results — this is the machine-readable record of the
run:

```bash
cat out/cavity/manifest.json
```

A completed manifest looks like this (fields you report to the user):

```json
{ "scenario": "cavity", "status": "completed", "stepsRun": 20000,
  "wallSeconds": 69.34, "mlups": 4.725,
  "diagnostics": { "totalMass": 15876.0, "maxSpeed": 0.09531, "tau": 0.56 },
  "warnings": [], "files": [ "speed_20000.png" ] }
```

The `files` array lists every output file written into `<out>/` — report those
paths so the user can open them. For the gallery, the entry point is
`<out>/index.html` (a self-contained HTML file with embedded PNGs — openable
directly in a browser, no server needed).

If `status` is `diverged` or `failed`, the preset did not complete: report the
status and the `diagnostics` verbatim. Do not claim success. (Presets are tuned
to complete, so a divergence usually means a broken build or environment, not a
bad config — the presets themselves are fixed.)

## Worked example (end-to-end)

Task: "Run the Kármán vortex demo and tell me where the images are."

1. **Which command (decision table):** run one named preset → `presets run`.
2. **Run it:**
   `./target/release/lbm presets run cylinder-karman --out out/karman`
3. **Watch the final line:** `status=completed steps=40000 wall=163.0s …`. Good.
4. **Read the manifest:** `cat out/karman/manifest.json` → `"status":"completed"`,
   `files` lists the PNG(s) and `force.csv`.
5. **Report:** "Completed in ~163 s. Outputs in `out/karman/`:
   `speed_40000.png`, `force.csv` (drag/lift time series), `manifest.json`."

## Top failure modes (and the fix)

- **`lbm` binary not found.** It is a release artifact. Fix: run
  `cargo build --workspace --release` first, then re-run the preset.
- **Invented a preset name.** Only the four names in the table above exist. Fix:
  run `presets list` to confirm; if the user wants something else, hand off to
  `lbmflow-user-author-scenario`.
- **No `--out`, then "where are the files?"** Without `--out`, results default to
  `out/<scenario name>`. Always pass `--out <dir>` so the location is explicit.
- **Reported success on a `diverged` run.** Check the manifest `status`, not just
  that the command exited. `diverged`/`failed` are NOT success.
- **Ran the gallery expecting one preset.** `gallery` runs ALL four (~5 min). For
  a single preset use `presets run <name>`.
- **User wants to tweak a preset (change nu, collision, grid).** That is no
  longer "run a preset" — presets are fixed. Route the tweak to
  `lbmflow-user-author-scenario` (use `presets show <name>` as the starting
  template) and, for stability of the tweak, `lbmflow-user-tune-stability`.
