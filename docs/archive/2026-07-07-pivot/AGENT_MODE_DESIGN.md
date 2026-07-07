# Agent Mode Design (Phase 6)

**Status (2026-07-07)**: **Landed**: `lbm` CLI, scenario runner, schema printing, presets, and stdio MCP server are implemented.
**Landed**: MCP exposes the designed seven tools, including async `start_run`/`run_status`/`list_runs` with a 4-run cap.
**Superseded/current intent**: the JSON example below is historical; authoritative schema is `crates/lbm-scenario/src/lib.rs`, which now includes 3D, units, compute backend, rotor, particles, and related fields.

Enable agents (Claude/Codex, etc.) and CLI users to run simulations through the
same entry point. **The scenario JSON is the single execution contract**, and the
GUI presets also internally generate the same schema (unifying the 3 modes).

## Crate structure

(landed with additions 2026-07-07 — CLI subcommands exist; `gallery`, checkpoint flags, and JSON manifest output were added beyond this sketch)

- `crates/lbm-cli` → binary name `lbm`
  - `lbm run <scenario.json> [--out DIR]`: executes and outputs artifacts + manifest.json
  - `lbm validate <scenario.json>`: schema/physics validity check only (does not execute)
  - `lbm presets [list|show NAME|run NAME]`: built-in presets
  - `lbm schema`: outputs the scenario JSON Schema to stdout (for agent self-discovery)
  - `lbm mcp`: starts an MCP server over stdio
- Dependencies: serde/serde_json, clap, png (visualization output), (MCP is either rmcp or hand-written JSON-RPC)

## Scenario JSON schema (v0)

(partially superseded 2026-07-07 — use `lbm schema`/`lbm-scenario::Scenario`; this example omits newer fields and includes `collision.magic`, which is not in the current `CollisionSpec`)

```jsonc
{
  "version": 0,
  "name": "cylinder-re100",                 // also used as the output directory name
  "grid": { "nx": 440, "ny": 160 },
  "physics": {
    "nu": 0.01,                              // lattice units
    "collision": { "type": "trt", "magic": 0.1875 },   // {"type":"bgk"} also allowed
    "force": [0.0, 0.0],
    "precision": "f64"                       // "f32" | "f64" (precision/speed trade-off)
  },
  "edges": {                                 // 1:1 with lbm-core's EdgeBC
    "left":   { "type": "velocityInlet", "u": [0.05, 0.0] },
    "right":  { "type": "outflow" },
    "bottom": { "type": "bounceBack" },
    "top":    { "type": "bounceBack" }
  },
  "inletProfile": { "edge": "left", "kind": "parabolic", "umax": 0.05 },  // optional
  "obstacles": [                             // expanded via set_solid_region
    { "shape": "circle", "cx": 110, "cy": 80, "r": 12 },
    { "shape": "rect", "x0": 0, "y0": 0, "x1": 10, "y1": 5 }
  ],
  "init": { "kind": "rest" },                // rest | taylorGreen | custom (future)
  "run": {
    "steps": 100000,
    "stopWhenSteady": { "epsilon": 1e-11, "checkEvery": 500 }   // optional
  },
  "probes": [                                // time-series recording
    { "type": "force", "target": "obstacles", "every": 10 },
    { "type": "point", "x": 220, "y": 80, "fields": ["ux","uy","rho"], "every": 100 }
  ],
  "outputs": [
    { "field": "speed", "at": "end", "format": "png", "colormap": "viridis" },
    { "field": "ux",    "at": "end", "format": "csv" },
    { "snapshotEvery": 10000, "field": "vorticity", "format": "png" }
  ]
}
```

### Execution results (out directory)

(landed with manifest drift 2026-07-07 — actual manifest uses camelCase `stepsRun`, status `completed|steady|diverged`, and errors are returned by the runner/MCP call rather than persisted as `status="error"`)

- `manifest.json`: machine-readable summary
  ```jsonc
  {
    "scenario": "cylinder-re100",
    "status": "completed" | "diverged" | "error",
    "steps": 100000, "wallSeconds": 42.1, "mlups": 168.3,
    "diagnostics": { "totalMass": ..., "maxSpeed": ..., "reynolds": ... },
    "warnings": ["..."],
    "files": [ {"path": "speed_100000.png", "kind": "field", ...} ]
  }
  ```
- probes are CSV (`force.csv`: step,fx,fy)
- If NaN is detected mid-run → status="diverged" and it terminates immediately, retaining
  the last diagnostics (includes a cause-hint string so the agent can fix parameters and retry)

### Design principles (agent UX)

(landed with scope drift 2026-07-07 — self-description and deterministic runs exist; validation reports JSON errors/warnings/units, not the exact field-level repair schema sketched here)

1. **Self-describing**: usage can be discovered with just `lbm schema` and `lbm presets list`
2. **Structured failures**: validation errors return JSON stating "which field, why,
   and how to fix it" (e.g. "nu must be > 0; tau = 3*nu + 0.5 must exceed 0.5")
3. **Numerical guardrails**: at the validate stage, stability heuristics (|u|>0.15 warning,
   grid Reynolds number U/ν > 15 warning, etc.) are emitted into warnings
4. **Determinism**: same scenario → same result (maintains the seed-free design)

## MCP server (`lbm mcp`)

(landed 2026-07-07 — all seven designed tools are registered in `tools_list` and dispatched in `tools_call`)

The stdio MCP server exposes 7 tools. The synchronous `run_scenario` blocks the
client; the async trio (`start_run` / `run_status` / `list_runs`) is the way to
launch long or multiple simulations without blocking (4-concurrent cap).

| Tool | Purpose |
|---|---|
| `run_scenario(scenario, outDir?)` | Blocking run; returns manifest on completion. |
| `start_run(scenario, outDir?) -> runId` | Non-blocking; returns a run handle. |
| `run_status(runId) -> {state, progress, manifest?}` | Poll a running/finished run. |
| `list_runs() -> [{runId, state, ...}]` | Enumerate active + recent runs. |
| `validate_scenario(scenario) -> {ok, errors[], warnings[]}` | Schema + stability heuristics only. |
| `list_presets() -> [{name, description, scenario}]` | Built-in preset catalog. |
| `get_schema() -> JSON Schema` | Self-discovery of the scenario schema. |

## Relationship with the GUI

(partially landed 2026-07-07 — GUI scenario export exists; presets still use GUI `EngineConfig` and are converted when exported)

- GUI presets = scenario JSON (generated from web/src/presets.ts or shared JSON)
- Future: an "export scenario" button in the GUI → reproducible in Agent Mode
