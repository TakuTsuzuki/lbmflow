# Agent Mode Design (Phase 6)

Enable agents (Claude/Codex, etc.) and CLI users to run simulations through the
same entry point. **The scenario JSON is the single execution contract**, and the
GUI presets also internally generate the same schema (unifying the 3 modes).

## Crate structure

- `crates/lbm-cli` → binary name `lbm`
  - `lbm run <scenario.json> [--out DIR]`: executes and outputs artifacts + manifest.json
  - `lbm validate <scenario.json>`: schema/physics validity check only (does not execute)
  - `lbm presets [list|show NAME|run NAME]`: built-in presets
  - `lbm schema`: outputs the scenario JSON Schema to stdout (for agent self-discovery)
  - `lbm mcp`: starts an MCP server over stdio
- Dependencies: serde/serde_json, clap, png (visualization output), (MCP is either rmcp or hand-written JSON-RPC)

## Scenario JSON schema (v0)

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

1. **Self-describing**: usage can be discovered with just `lbm schema` and `lbm presets list`
2. **Structured failures**: validation errors return JSON stating "which field, why,
   and how to fix it" (e.g. "nu must be > 0; tau = 3*nu + 0.5 must exceed 0.5")
3. **Numerical guardrails**: at the validate stage, stability heuristics (|u|>0.15 warning,
   grid Reynolds number U/ν > 15 warning, etc.) are emitted into warnings
4. **Determinism**: same scenario → same result (maintains the seed-free design)

## MCP server (`lbm mcp`)

Tools:
- `run_scenario(scenario: object, outDir?: string) -> manifest` (supports progress notifications)
- `validate_scenario(scenario: object) -> {ok, errors[], warnings[]}`
- `list_presets() -> [{name, description, scenario}]`
- `get_schema() -> JSON Schema`
- `read_field(runDir, field, format="csv") -> data` (retrieval of completed run results)

## Relationship with the GUI

- GUI presets = scenario JSON (generated from web/src/presets.ts or shared JSON)
- Future: an "export scenario" button in the GUI → reproducible in Agent Mode
