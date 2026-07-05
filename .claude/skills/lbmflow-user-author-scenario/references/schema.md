# LBMFlow scenario JSON (v0) — accepted schema reference

Authoritative accepted-value list, transcribed from
`crates/lbm-scenario/src/lib.rs` on `main`. JSON keys are **camelCase**. Top-level
objects and most nested ones use `deny_unknown_fields` — an unknown key is a hard
`invalid-scenario-json` error. Recommend ONLY what is listed here.

## Table of contents
1. Top-level shape
2. grid
3. physics (nu, collision, force, precision)
4. compute (backend)
5. edges (7 BC types) + inletProfile
6. obstacles (3 shapes)
7. init
8. multiphase (2D only)
9. run / steady
10. probes
11. outputs (fields × formats)
12. 3D limits (build-time rejections)

## 1. Top-level shape

Required: `name` (string), `grid`, `physics`, `edges`, `run`.
Optional: `version` (u32, default 0), `compute`, `inletProfile`, `obstacles` (array),
`init` (default `rest`), `multiphase`, `probes` (array), `outputs` (array).

## 2. grid

`{ "nx": <usize>, "ny": <usize>, "nz": <usize?> }`
- `nz` omitted or `1` → 2D (D2Q9). `nz > 1` → 3D (D3Q19); 3D requires `nz >= 3`
  (and `nx,ny >= 3`).

## 3. physics

```json
"physics": {
  "nu": 0.02,                       // kinematic viscosity, LATTICE units; tau = 3*nu + 0.5
  "collision": { "type": "trt" },   // "trt" (default, recommended) | "bgk"  — NOTHING ELSE
  "force": [0.0, 0.0],              // uniform body force (lattice units). 2-vector; z is 0 in 3D
  "precision": "f64"                // "f64" (default) | "f32"
}
```
- `collision.type` accepts exactly `bgk` and `trt`. MRT/cumulant/regularized do
  NOT exist (`unknown variant` error). Scheme choice → `lbmflow-user-tune-stability`.

## 4. compute (optional)

```json
"compute": { "backend": "auto" }    // "auto" (default) | "cpu" | "gpu"
```
- **`gpu` does not run.** 3D rejects it at build time. (In-flight fix: 2D will
  also reject `gpu` at validate time instead of silently falling back to CPU —
  author against that fixed behavior: never recommend `gpu`.) Use `auto` or `cpu`.

## 5. edges

Required faces: `left`, `right`, `bottom`, `top`. 3D-only optional: `front`
(z=0), `back` (z=nz-1); omitted = periodic. Each face is one of:

| `type` | Fields | Notes |
|---|---|---|
| `periodic` | — | must pair on the same axis |
| `bounceBack` | — | stationary no-slip wall |
| `movingWall` | `u: [ux, uy]` | driven wall |
| `velocityInlet` | `u: [ux, uy]` | prescribed-velocity inflow |
| `pressureOutlet` | `rho: <f64>` | fixed density/pressure outflow |
| `outflow` | — | zero-gradient |
| `convectiveOutflow` | `uConv: <f64 in (0,1]>` | low-reflection; ≈ mean outflow speed |

`inletProfile` (optional, one object):
```json
"inletProfile": { "edge": "left", "kind": "parabolic", "umax": 0.1 }
```
- `edge`: `left|right|bottom|top`. `kind`: `parabolic` (only value).

## 6. obstacles (array; each element is one shape)

| `shape` | Fields | Dimensionality |
|---|---|---|
| `circle` | `cx, cy, r` (f64) | 2D disk; 3D = z-extruded cylinder |
| `rect` | `x0, y0, x1, y1` (usize cells) | 2D rectangle; 3D = z-extruded box |
| `sphere` | `cx, cy, cz, r` (f64) | 3D only — solid ball |

There is NO mesh/CAD/STL/OBJ/STEP import. Geometry = these primitives only,
staircase-voxelized on the uniform lattice.

## 7. init

```json
"init": { "kind": "rest" }          // default
```
- `rest` (default). `droplet { cx, cy, r, rhoLiquid, rhoVapor }` and
  `pool { heightFrac, rhoLiquid, rhoVapor }` are 2D-multiphase only. **3D: `rest`
  only** (anything else is build-rejected).

## 8. multiphase (2D ONLY)

```json
"multiphase": { "g": -5.0, "gWall": 0.0, "wallRho": 1.0 }
```
- `g`: Shan-Chen cohesion (negative; `-5.0` validated; `g > -4` won't phase-
  separate — warned). `wallRho` (optional) controls contact angle (preferred over
  `gWall`). **Rejected in 3D** (`3D … 未対応: multiphase`).

## 9. run

```json
"run": { "steps": 20000, "stopWhenSteady": { "epsilon": 1e-6, "checkEvery": 500 } }
```
- `steps` required. `stopWhenSteady` optional (early-stop on steady state).

## 10. probes (array)

| `type` | Fields | Output |
|---|---|---|
| `force` | `every` | `force.csv` (obstacle momentum-exchange force) |
| `point` | `x, y, z?, every` | `point_<x>_<y>[_<z>].csv` time series |

## 11. outputs (array)

```json
"outputs": [ { "field": "speed", "format": "png", "every": 0 } ]
```
- `field`: `speed | ux | uy | rho | vorticity`.
- `format`: `png | csv | vtk`.
- `every`: snapshot every N steps; `0` = only at the end.
- Post-processing / reading these files → `lbmflow-user-postprocess`.

## 12. 3D limits (build-time rejections — do NOT emit these in 3D)

When `grid.nz > 1`:
- `multiphase` present → error `3D (nz > 1) では未対応: multiphase（多相流）`.
- `init` other than `rest` → error `init は rest のみ`.
- `compute.backend: "gpu"` → error `compute.backend "gpu"（cpu / auto を指定…）`.

3D supported set: single-phase, `init:rest`, cpu/auto backend, the 7 edge BCs,
`circle`/`rect`/`sphere` obstacles, force/point probes, all fields/formats
(PNG/CSV are z-mid slices; VTK is the full volume).
