---
name: lbmflow-qa-viewer
description: >-
  Build a self-contained, interactive 3D viewer for an exported LBMFlow field
  volume — a linked 3D scene (Three.js/WebGL, orbit + zoom) plus horizontal and
  vertical cross-sections that scan with sliders, a field selector (speed, shear,
  vorticity, Q-criterion, axial/radial/swirl velocity, kinetic energy), an
  optional multi-dataset switch (e.g. stir-speed levels), a shear-stress
  threshold check, and Stokes-settling particles by diameter — all in physical
  units. Use whenever you want to LOOK at a 3D result interactively, "make a
  viewer / dashboard", "visualize the volume in 3D", "compare fields / cases
  side by side", or produce a shareable Artifact of a run. Owns the
  volume→self-contained-HTML pipeline (inline Three.js + base64 data, so it opens
  offline by double-click AND passes the Artifact CSP). Do NOT use it to compute
  the fields (run/example step) or to JUDGE them for defects (that is the
  core-owned anomaly-scan, which encodes VALIDATION-band semantics) — this Skill
  only turns exported data into an interactive visualization.
---

# LBMFlow — interactive 3D QA viewer

The visualization step of the physics-QA loop: **run → export → scan →
[visualize] → report**. Turns an exported field volume into one self-contained
HTML page — no server, no external requests — that you can open by double-click,
screenshot for a report, or publish as an Artifact.

## Input format

`volume.bin` (f32 LE, `(k*vn+j)*vn+i`, comps `[vx,vy,vz,shear]`) + `volume.json`
(`{vn,n,cx,cy,r_tank,zc,tip_r,disk_r,hub_r,shaft_r,blade_hh,disk_hh,n_blades,
blade_hw,baffle_len,baffle_hw,omega,u_tip,speed_max,shear_max}`). Same format the
stirred-tank example writes and the core-owned anomaly-scan reads. The viewer
derives vorticity / Q-criterion in-page from the velocity field; shear is read
from comp 3 (prefer the core `gather_shear_rate` value at export time).

## One-time: fetch the (inlined) libraries

The page inlines Three.js r128 + OrbitControls so it has zero external requests.
Fetch once into your working dir:

```bash
curl -sL https://unpkg.com/three@0.128.0/build/three.min.js -o three.min.js
curl -sL https://unpkg.com/three@0.128.0/examples/js/controls/OrbitControls.js -o OrbitControls.js
```

## Build a viewer

Assets live in `assets/` next to this file:
- `dashboard_template.html` — the linked multi-field / multi-dataset dashboard.
- `build_dashboard.py` — assembles a BUNDLE (geometry + one or more datasets) and
  inlines Three.js/OrbitControls/data into a full local HTML.
- `inline_single.py` — the simpler single-volume inliner (tokens
  `/*__THREE__*/ /*__ORBIT__*/ /*__META__*/ /*__DATA__*/`).

Local (full-res), one or more datasets — label:dir pairs:

```bash
python3 assets/build_dashboard.py assets/dashboard_template.html \
  three.min.js OrbitControls.js  out/dashboard.html \
  Slow:runs/slow  Medium:runs/med  Fast:runs/fast
```

Single dataset also works: pass one `Label:dir`.

## Publish as an Artifact (optional)

The template is a full HTML doc; the Artifact host supplies `<head>`/`<body>`, so
strip the wrapper and prepend a charset meta (keeps the µ/ρ/τ glyphs correct).
Downsample the volume first (≈40³) to keep the Artifact light (~2–5 MB):

```python
# strip wrapper -> artifact-ready inner HTML
import re, pathlib
h = pathlib.Path("out/dashboard.html").read_text()
title = re.search(r"<title>.*?</title>", h, re.S).group(0)
style = re.search(r"<style>.*?</style>", h, re.S).group(0)
body  = re.search(r"<body>(.*)</body>", h, re.S).group(1)
pathlib.Path("out/dashboard_artifact.html").write_text(
    f'<meta charset="utf-8">\n{title}\n{style}\n{body}')
```

Then publish `out/dashboard_artifact.html` with the Artifact tool (favicon 🌀).

## Verify before you trust a screenshot

Serve locally and check in the preview harness — WebGL bugs hide between source
and render:

```json
// .claude/launch.json
{ "version":"0.0.1", "configurations":[
  { "name":"static","runtimeExecutable":"python3","runtimeArgs":["-m","http.server","8899"],"port":8899 }]}
```

`preview_start` → navigate to `/out/dashboard.html` → `preview_console_logs`
(errors) → `preview_screenshot`. Known WebGL gotchas already handled in the
template but worth re-checking if you fork it:
- transparent vessel meshes must set `depthWrite:false` or they occlude the
  interior point cloud;
- `renderer.setSize(w,h)` (NOT `false`) so the canvas CSS size matches the host;
- `Object3D.position` is read-only in r128 — set `mesh.position.set(...)`, never
  `Object.assign(mesh,{position:...})`;
- keep all UI copy ASCII/entity-encoded (or prepend the charset meta) so glyphs
  survive a head-less host.

## Boundary

Physical scales (Δx, Δt, effective viscosity, Pa, m/s, rpm, Stokes v_s) are
computed IN the template for display. Per the QA hand-off these belong in a core
UnitConverter long-term; until it lands, the viewer's dimensionalization is a
display convenience with the effective-viscosity caveat shown in-page. Field
math the viewer still does client-side (vorticity, Q) should be swapped for the
core FieldKinds when they land.
