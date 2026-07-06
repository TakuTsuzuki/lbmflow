---
name: lbmflow-user-author-scenario
description: >-
  Author a valid LBMFlow scenario JSON (2D D2Q9 or 3D D3Q19) from a natural-
  language description, compose obstacle geometry from primitives, and validate
  it before running. Use whenever the user wants to "set up a simulation", "make
  a scenario", "simulate flow past a cylinder / around obstacles / in a channel",
  "write the scenario JSON", "add a sphere/cylinder/box obstacle", "put three
  staggered cylinders", or describes a flow they want to model in words. This
  Skill owns the scenario schema, the accepted enum values, obstacle composition
  from circle/rect/sphere primitives, and the `lbm validate` gate. Inputs are RAW
  LATTICE UNITS — this Skill does NOT convert SI/physical units (see the routing
  rule for unit-conversion requests). Do NOT use it to RUN a built-in preset
  as-is (that is lbmflow-user-run-preset), to pick collision scheme or fix a
  stability warning (that is lbmflow-user-tune-stability), or to run/monitor the
  finished scenario (CLI run = the run-preset Skill's runner; async = run-monitor-
  mcp).
---

# LBMFlow — author & validate a scenario

A scenario is a single JSON document that fully specifies one simulation. This
Skill turns a natural-language flow description into a scenario that **passes
`lbm validate`** — that validation pass is the done-criterion. It covers the 2D
(D2Q9) and 3D (D3Q19) paths and composing obstacle geometry from the three
primitive shapes.

**The one rule that prevents wrong physics: every number is in LATTICE UNITS.**
There is no SI/physical-unit conversion anywhere in LBMFlow. `nu`, all velocities
(`u`), `rho`, and `force` are raw lattice quantities. If the user gives physical
units ("water at 20°C", "1 m/s inlet", "Reynolds 200 in a 5 cm pipe"), see
**"Unit conversion is out of scope"** below — do NOT silently treat physical
numbers as lattice numbers.

The complete accepted schema — every field, every enum value, with the exact
serde names — is in [references/schema.md](references/schema.md). Read it before
emitting JSON. Recommend ONLY values that appear there; anything else is rejected
by `deny_unknown_fields` / unknown-variant errors.

## Step 0 — 2D or 3D?

| Signal in the request | Path | Consequence |
|---|---|---|
| Planar flow, "2D", cavity, channel, vortex street, ANY multiphase/droplet | **2D** — omit `grid.nz` (or set `nz:1`) | D2Q9; multiphase allowed |
| "3D", a volume, a sphere/ball, flow "through" a duct, z-extent named | **3D** — set `grid.nz > 1` (needs `nz >= 3`) | D3Q19, CPU only, single-phase, `init:rest` only |

3D is real and runs, but it is narrow: **no multiphase, `init` must be `rest`,
backend must be cpu/auto.** If the user asks for a 3D droplet or 3D multiphase,
that is rejected at build time — say so (route to the red spec note in
references/schema.md §"3D limits"), do not emit it.

## Step 1 — Fill the required skeleton

Every scenario needs `name`, `grid`, `physics`, `edges`, `run`. Minimal 2D
skeleton (fill the values from the request; keep everything in lattice units):

```json
{
  "name": "channel-flow",
  "grid": { "nx": 200, "ny": 80 },
  "physics": { "nu": 0.02, "collision": { "type": "trt" }, "precision": "f64" },
  "edges": {
    "left":   { "type": "velocityInlet", "u": [0.05, 0.0] },
    "right":  { "type": "pressureOutlet", "rho": 1.0 },
    "top":    { "type": "bounceBack" },
    "bottom": { "type": "bounceBack" }
  },
  "run": { "steps": 20000 },
  "outputs": [ { "field": "speed", "format": "png", "every": 0 } ]
}
```

For 3D, add `"nz"` to `grid` and OPTIONAL `front`/`back` z-faces to `edges`
(omitted = periodic). Defaults you can rely on (do not restate them unless the
user asks): `collision` = `trt`, `precision` = `f64`, `force` = `[0,0]`,
`init` = `rest`, `compute.backend` = `auto`.

## Step 2 — Boundary conditions (7 accepted edge types)

Set each of `left/right/top/bottom` (+ `front/back` in 3D) to one of exactly
these — the full list, with required fields:

| `type` | Required fields | Use for |
|---|---|---|
| `periodic` | — | wrap-around (must pair on the same axis) |
| `bounceBack` | — | stationary solid wall |
| `movingWall` | `u:[ux,uy]` | driven wall (e.g. cavity lid) |
| `velocityInlet` | `u:[ux,uy]` | prescribed-velocity inflow |
| `pressureOutlet` | `rho` | fixed-density (pressure) outflow |
| `outflow` | — | zero-gradient outflow |
| `convectiveOutflow` | `uConv` (in (0,1]) | low-reflection outflow; `uConv` ≈ mean outflow speed |

`inletProfile` (optional) turns a flat inlet into a Poiseuille parabola:
`{ "edge": "left", "kind": "parabolic", "umax": 0.1 }`. Only `parabolic` exists.

## Step 3 — Compose obstacle geometry from primitives

Obstacles are an array of primitive shapes voxelized onto the uniform lattice.
**Only three shapes exist** — compose them; there is NO mesh/CAD/STL import.

| `shape` | Fields | Dimensionality |
|---|---|---|
| `circle` | `cx, cy, r` | 2D disk; in 3D = a cylinder extruded along z |
| `rect` | `x0, y0, x1, y1` (integer cells) | 2D rectangle; in 3D = a box extruded along z |
| `sphere` | `cx, cy, cz, r` | 3D only — a solid ball |

Composition = emit one array element per shape. Examples:

- "three staggered cylinders in a 2D channel" → three `circle` entries at
  different `(cx,cy)`.
- "a backward-facing step" → one `rect` blocking part of the inlet height.
- "a ball in a 3D duct" → one `sphere` (requires `grid.nz > 1`).

```json
"obstacles": [
  { "shape": "circle", "cx": 60, "cy": 40, "r": 8 },
  { "shape": "circle", "cx": 90, "cy": 24, "r": 8 },
  { "shape": "circle", "cx": 90, "cy": 56, "r": 8 }
]
```

To measure force on obstacles, add a probe:
`"probes": [ { "type": "force", "every": 100 } ]` (writes `force.csv`).

## Step 4 — Validate (the done gate)

The scenario is NOT done until `lbm validate` accepts it. Validate reports two
things: hard build errors (`ok:false` + `error`) and non-fatal stability
warnings. Run:

```bash
./target/release/lbm validate scenario.json
# or pipe JSON on stdin:
cat scenario.json | ./target/release/lbm validate -
```

| Validate output | Meaning | Your action |
|---|---|---|
| `{ "error": null, "ok": true, "warnings": [] }` | Fully valid, no concerns | **Done.** |
| `ok:true` with `warnings:[…]` | Valid but physically risky (tau/Mach/grid-Re) | Valid to run, but the warnings are a stability matter → route to `lbmflow-user-tune-stability`. |
| `ok:false` + `error:"invalid-scenario-json"` | Schema violation (unknown field/variant) | Fix against references/schema.md and re-validate. |
| `ok:false` + a build error (e.g. `3D … 未対応: multiphase`) | Feature not supported in this configuration | Remove the unsupported feature; explain the limit. Do not retry with the same JSON. |

A schema error is precise — it names the bad token, e.g.
`unknown variant \`mrt\`, expected \`bgk\` or \`trt\``. Trust it and fix that token.

## Unit conversion is out of scope (route it, never fake it)

LBMFlow has NO SI→lattice conversion. If the request contains physical units
(m/s, Pa, °C, cm, real Reynolds numbers tied to physical size):

1. **Do not** paste physical numbers into `nu`/`u`/`rho`/`force` as if they were
   lattice units — that produces silently wrong physics.
2. State plainly that scenario inputs are lattice units and that dimensional
   conversion is a separate concern this Skill does not perform.
3. You MAY still author a scenario using directly-given LATTICE quantities, or
   ask the user for lattice-unit values. Non-dimensionalization guidance (picking
   `nu`, `u`, grid size to hit a target Reynolds number in lattice units) belongs
   to the unit-conversion spec note — see `docs/skills/b2-skill-specs.md`
   §"Red spec notes → unit conversion". Point the user there rather than
   improvising a conversion inside the scenario.

## Worked example (end-to-end)

Task: "Simulate 2D flow past a cylinder in a channel, parabolic inlet, and
measure the drag."

1. **2D or 3D (Step 0):** planar, "2D" → 2D, omit `nz`.
2. **Skeleton (Step 1):** `grid 260×80`, `nu 0.02`, `trt`, `f64`.
3. **BCs (Step 2):** `left velocityInlet u:[0.08,0]`, `right pressureOutlet
   rho:1.0`, `top`/`bottom bounceBack`; add `inletProfile` parabolic on `left`,
   `umax:0.1`.
4. **Geometry (Step 3):** one `circle` `cx:60,cy:40,r:9`; add
   `probes:[{type:force,every:100}]` for drag.
5. **Validate (Step 4):** `lbm validate scenario.json` → `ok:true`. If a grid-Re
   warning appears, hand the stability question to `lbmflow-user-tune-stability`.
6. **Report:** the JSON + "validated clean; run it with the runner or async MCP."

## Top failure modes (and the fix)

- **Physical units pasted as lattice units.** Symptom: `nu`/`u` look like SI
  numbers; results are nonsense or diverge. Fix: keep everything lattice; route
  conversion to the spec note (see above).
- **Recommended a scheme/feature not in the schema.** MRT, cumulant, STL import,
  GPU, 3D multiphase — none exist. Fix: recommend only values in
  references/schema.md; `validate` will reject the rest.
- **3D multiphase / 3D non-rest init.** Build-time rejected. Fix: 3D is
  single-phase `init:rest` only; do 2D for multiphase.
- **Unpaired periodic edges.** `periodic` on `left` but not `right` errors at
  build. Fix: periodic faces must pair on the same axis.
- **`convectiveOutflow` with `uConv` outside (0,1].** Build error. Fix: set
  `uConv` near the expected mean outflow speed (e.g. 0.05–0.15).
- **`rect` coordinates as floats.** `x0..y1` are integer cell indices. Fix: use
  whole numbers; `circle`/`sphere` centers/radii may be fractional.
- **Declared done without validating.** Emitting JSON is not done. Fix: `done`
  requires an `lbm validate` pass (Step 4).
