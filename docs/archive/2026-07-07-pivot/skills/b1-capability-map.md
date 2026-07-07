# B1 Capability Map — LBMFlow (current `main`)

> **Status note (2026-07-07, PM):** this map is a dated snapshot (2026-07-05).
> Its GPU findings are superseded: 2D `backend:"gpu"` scenario dispatch is now
> landed (feature `gpu`, f32 only) and unsupported cases raise an explicit
> error instead of silently falling back to CPU
> (`crates/lbm-scenario/src/lib.rs:110-112`, `crates/lbm-cli/src/runner.rs`).
> See the README capability matrix for current status. The historical trap
> documented in §4.4 is retained on purpose.

**Track B (Skills for User), session B1.** This report classifies what a USER can do
on current `main`, per workflow stage, color-coded with hard evidence. Source of truth is
**executable CLI/MCP/schema on `main` only** — not docs claims, not branches, not planned
features. No Skill is authored here; this is a capability map for PM review before B2.

---

## 0. Base SHA, environment, and worktree note

| Item | Value |
|------|-------|
| Base SHA (worktree HEAD) | `e2442ea61856b7cb0d0bd7ee763a45be4943e55d` |
| Branch | `skills/b1-capability-map` |
| `main` HEAD at checkout | `e2442ea` (identical — worktree created from local main) |
| Toolchain | `rustc 1.93.0 (254b59607 2026-01-19)`, `cargo 1.93.0` |
| Platform | Darwin arm64, macOS 26.4 (build 25E246) |
| Build | `cargo build --workspace --release` → exit 0; binary `target/release/lbm` (1.94 MB) |

**Worktree/remote deviation (recorded per brief):** The brief said the repo has NO remote
configured and instructed creating the worktree from local main HEAD. **This is now inaccurate:**
`git remote -v` in the source repo shows a real remote —
`origin https://github.com/TakuTsuzuki/lbmflow.git` (fetch+push). The worktree already existed at
`/Users/taku/projects/lbmflow-b1` on branch `skills/b1-capability-map`, based on local main HEAD
`e2442ea` (matching main). No `git fetch`/`git worktree add` was re-run. The branch contains only
this report file; no source/docs/tests/schema were modified. Untracked `target/` (build output)
is present as expected.

---

## 1. Per-stage color table (summary)

| Stage | Color | One-line verdict |
|-------|-------|------------------|
| Geometry / modeling | 🟡 Yellow | Only parametric primitives (circle/rect/sphere) via scenario JSON. No CAD/mesh import. |
| Mesh / lattice | 🟢 Green | Uniform Cartesian lattice via `grid{nx,ny[,nz]}`; 2D D2Q9 / 3D D3Q19 both execute. No refinement/unstructured. |
| Scheme selection | 🟢 Green | `physics.collision` = `bgk`\|`trt` both execute; precision f32/f64; backend cpu/auto. |
| BCs & properties — BCs | 🟢 Green | 7 edge BC types + inlet profiles + obstacles all execute. |
| BCs & properties — unit→lattice | 🔴 Red | No SI→lattice conversion anywhere; all inputs are raw lattice units. |
| Run / monitor (async MCP) | 🟢 Green | `lbm mcp` stdio; async `start_run`/`run_status`/`list_runs` verified working. |
| Post / visualize | 🟢 Green | PNG/CSV/VTK outputs + self-contained HTML gallery, all executed. |
| Validate — config | 🟢 Green | `lbm validate` returns errors + stability warnings (tau/Mach/grid-Re). |
| Validate — compare vs. reference | 🔴 Red | No user-facing comparison/analytic-check command; only developer `cargo test`. |

---

## 2. Green — evidence blocks (command + pasted output)

### 2.1 CLI surface (self-discovery)

```
$ ./target/release/lbm --help
Commands:
  run       シナリオ JSON を実行し、結果を出力ディレクトリに書き出す
  validate  シナリオを実行せずに検証する（エラー/警告を JSON で報告）
  presets   組み込みプリセットの一覧・表示・実行
  gallery   全プリセットを順に実行し、自己完結 HTML ギャラリー（index.html）を生成する
  schema    シナリオ JSON の書式説明を出力する（エージェントの自己発見用）
  mcp       MCP サーバーとして stdio で待ち受ける（AI エージェント連携）
```

### 2.2 Presets — list + run one (GREEN)

```
$ ./target/release/lbm presets list
cavity               リッド駆動キャビティ（定常判定つき）
cylinder-karman      円柱まわりのカルマン渦列 + 抗力プローブ
two-phase-droplet    Shan-Chen 二相液滴の平衡化
droplet-on-wall      壁上液滴の接触角デモ（仮想壁密度 wallRho=1.0 → θ≈63°）

$ ./target/release/lbm presets run cavity --out <SP>/cavity_out
status=completed steps=20000 wall=69.3s mlups=5 out=.../cavity_out
  speed_20000.png
```

`manifest.json` produced:
```json
{ "scenario": "cavity", "status": "completed", "stepsRun": 20000,
  "wallSeconds": 69.34, "mlups": 4.725,
  "diagnostics": { "totalMass": 15876.0, "maxSpeed": 0.09531, "tau": 0.56 },
  "warnings": [], "files": [ "speed_20000.png" ] }
```

### 2.3 Scheme selection — `lbm schema` (GREEN; fixed ruling confirmed)

`lbm schema` output (excerpt — this is the authoritative accepted-options list a Skill may recommend):

```
  "physics": {
    "nu": 0.02,                            // 動粘性係数（格子単位）。tau = 3*nu + 0.5
    "collision": { "type": "trt" },        // "trt"（推奨） | "bgk"
    "force": [0.0, 0.0],                   // 一様体積力（重力など。3D では z 成分 0）
    "precision": "f64"                     // "f32" | "f64"
  },
  "compute": { "backend": "auto" },        // 省略可: "auto" | "cpu" | "gpu"（gpu は未提供）
  "edges": {
    "left":   { "type": "velocityInlet", "u": [0.1, 0.0] },
    "right":  { "type": "pressureOutlet", "rho": 1.0 },
    "bottom": { "type": "bounceBack" }, "top": { "type": "bounceBack" }
    // 他: periodic, movingWall{u}, outflow, convectiveOutflow{uConv}
  },
  ...
  // field: speed | ux | uy | rho | vorticity     format: png | csv | vtk
```

**Executable collision options on `main`: exactly `bgk` and `trt`.** A Skill may recommend only
these two (NL selection among existing executable options = GREEN per fixed ruling). Any other
scheme (MRT, cumulant, etc.) is not in the enum — see Red §4.1 for negative evidence.

### 2.4 BCs & properties — BC coverage (GREEN)

Source enum `EdgeSpec` (crates/lbm-scenario/src/lib.rs) exposes 7 edge types, all wired to core:
`Periodic, BounceBack, MovingWall{u}, VelocityInlet{u}, PressureOutlet{rho}, Outflow,
ConvectiveOutflow{uConv}`. Plus optional `inletProfile` (parabolic), `obstacles`, `multiphase`
(Shan-Chen, 2D only), `init` (rest/droplet/pool). All exercised by presets that run to completion
(§2.2, §2.9). This is the *boundary/property specification* half — GREEN.

### 2.5 3D scenario round-trip (GREEN) — write minimal 3D scenario with nz, run it

Wrote a 24×24×24 scenario (nz=24 → D3Q19), velocityInlet + pressureOutlet + bounceBack walls,
a `sphere` obstacle, force probe, and speed/rho/ux outputs.

```
$ ./target/release/lbm validate mini3d.json
{ "error": null, "ok": true, "warnings": [] }

$ ./target/release/lbm run mini3d.json --out <SP>/mini3d_out --json
{ "scenario": "mini3d", "status": "completed", "stepsRun": 200,
  "wallSeconds": 0.237, "mlups": 11.66,
  "diagnostics": { "totalMass": 13061.9, "maxSpeed": 0.0718, "tau": 0.65 },
  "warnings": [], "files": [ "force.csv", "speed_200.png", "rho_200.vtk", "ux_200.csv" ] }
```

Files on disk: `force.csv, manifest.json, rho_200.vtk (229 KB), speed_200.png, ux_200.csv`.

### 2.6 Output formats — VTK / PNG / CSV (GREEN; 2D and 3D)

VTK from the 3D run is a full 3D volume (not a slice):
```
$ head rho_200.vtk
# vtk DataFile Version 3.0
LBMFlow rho step=200
ASCII
DATASET STRUCTURED_POINTS
DIMENSIONS 24 24 24
ORIGIN 0 0 0
SPACING 1 1 1
POINT_DATA 13824
SCALARS rho double 1
LOOKUP_TABLE default
```
Format matrix (from crates/lbm-cli/src/runner.rs + observed runs):
- **PNG**: 2D full field; 3D = z-mid slice. Fields: speed, ux, uy, rho, vorticity.
- **CSV**: 2D full field; 3D = z-mid slice.
- **VTK**: 2D `DIMENSIONS nx ny 1`; 3D full volume `DIMENSIONS nx ny nz` (ParaView-openable).
- **Probe CSV**: `force.csv` (2D header `step,fx,fy`; 3D header `step,fx,fy,fz` — both observed),
  and point time-series `point_<x>_<y>[_<z>].csv`.

### 2.7 Config validation + error/warning reporting (GREEN — divergence/error UX)

**Invalid JSON (unknown collision variant) → structured error, exit 1:**
```
$ ./target/release/lbm validate bad.json     # collision.type = "mrt"
Error: {
  "error": "invalid-scenario-json",
  "hint": "lbm schema で書式を、lbm presets show <name> で実例を確認できます",
  "message": "unknown variant `mrt`, expected `bgk` or `trt` at line 1 column 97"
}
```

**Physically unstable but parseable → ok:false + graded warnings (tau / Mach / grid-Re):**
```
$ ./target/release/lbm validate unstable.json   # nu=1e-4, inlet u=0.5
{ "error": "prescribed speed 0.5 exceeds the low-Mach limit 0.3 (lattice units)",
  "ok": false,
  "warnings": [
    { "field": "physics.nu", "message": "tau = 0.500 は安定限界に近い（0.55 未満）..." },
    { "field": "edges", "message": "流入/壁速度 0.500 は圧縮性誤差が目立つ水準（0.15 超）" },
    { "field": "physics", "message": "グリッドレイノルズ数 U/ν = 5000.0 > 15: 発散の恐れ..." } ] }
```

**Runtime divergence → `status:"diverged"`, NaN detected, run halted early:**
```
$ ./target/release/lbm run diverge2.json --json   # nu=1e-7, inlet u=0.25, 5000 steps
{ "scenario": "diverge2", "status": "diverged", "stepsRun": 1000,
  "diagnostics": { "totalMass": null, "maxSpeed": 0.0, "tau": 0.5000003 }, ... }
```
(Requested 5000 steps; engine detected NaN and stopped at step 1000, reporting `diverged`.)

### 2.8 Gallery (GREEN) — all presets + self-contained HTML

```
$ ./target/release/lbm gallery --out <SP>/gallery_out
[gallery] cavity: status=completed steps=20000 wall=60.2s
[gallery] cylinder-karman: status=completed steps=40000 wall=163.0s
[gallery] two-phase-droplet: status=completed steps=20000 wall=56.9s
[gallery] droplet-on-wall: ... (completed)
```
Output: `index.html` (384 KB, `<!DOCTYPE html>` self-contained, base64-embedded PNGs) plus per-preset
subdirectories with manifests, PNGs, VTK, and force.csv. Full-pipeline (~5 min) runs green end to end.

### 2.9 Multiphase (GREEN, 2D only)

`two-phase-droplet` and `droplet-on-wall` presets (Shan-Chen single-component, contact-angle demo)
run to completion in the gallery pass. This is a real executable capability in 2D.
(3D multiphase is rejected — see Red §4.3.)

---

## 3. Yellow — small items (bounded notes)

### 3.1 Geometry / modeling — 🟡

**Entry point exists:** scenario `obstacles[]` accepts `circle{cx,cy,r}`, `rect{x0,y0,x1,y1}` (2D +
z-extruded in 3D), and `sphere{cx,cy,cz,r}` (3D). These execute (§2.5 sphere ran). So *parametric
primitive geometry* is actually GREEN. The **yellow** classification is for the natural next step a
user expects — **compositing/parameterizing primitive geometry from an NL request** (e.g. "put three
staggered cylinders", "a channel with a backward-facing step built from rects"). That is authorable
as a Skill that emits `obstacles[]` arrays.

**"small" check (all must hold):** compose existing primitive shapes into `obstacles[]` JSON.
- ≤8 h: yes (JSON generation from geometry description; primitives already exist).
- ≤5 files: yes (Skill file only — no source change needed; primitives already accepted).
- No new persistent data model: yes (uses existing `obstacles[]`).
- No new numerical method / geometry kernel: yes (reuses circle/rect/sphere staircase kernels).
- No new async subsystem / external service: yes.
- Validatable by existing smoke test: yes (`lbm validate` + a run).

> NOTE: CAD import, arbitrary meshes, and true 3D modeling are **NOT** yellow — see Red §4.2.
> Yellow here is strictly "generate `obstacles[]` from existing primitives".

### 3.2 (No other yellow items.)

Every other executable capability is either already GREEN (needs only a Skill wrapper, not new
engineering) or RED (requires a new subsystem and cannot be made to "pass" within the small budget).

---

## 4. Red — negative evidence (no Skill to be authored)

### 4.1 Collision/scheme selection beyond bgk/trt — 🔴

`lbm validate` with `collision.type = "mrt"` returns
`unknown variant \`mrt\`, expected \`bgk\` or \`trt\`` (§2.7). The `CollisionSpec` enum has exactly
two variants. MRT/cumulant/regularized are spec-only. A Skill must recommend ONLY bgk/trt.

### 4.2 CAD import / 3D modeling / mesh generation — 🔴

Negative evidence:
- No CLI subcommand accepts a geometry file. `lbm --help` exposes only `run/validate/presets/
  gallery/schema/mcp`; `lbm run` takes a scenario JSON (or `-` stdin), nothing else.
- No STL/STEP/OBJ/IGES/mesh parser anywhere: `grep` for `stl|step|obj|iges|gmsh|mesh import` in
  `crates/lbm-scenario` and `crates/lbm-cli` finds no import path. Geometry is limited to the three
  analytic primitives in §3.1 (staircase-voxelized on the uniform lattice).
- Lattice is uniform Cartesian only (`grid{nx,ny,nz}`); no refinement, no unstructured/body-fitted
  mesh. No executable pipeline + supported format exists on `main` → research/spec-only.

### 4.3 3D multiphase (and 3D non-rest init) — 🔴

```
$ ./target/release/lbm validate mp3d.json     # nz=16 + multiphase
{ "error": "3D (nz > 1) では未対応: multiphase（多相流）", "ok": false, ... }
```
Source builder rejects 3D multiphase and 3D non-`rest` init at build time. 3D is single-phase,
`init:rest` only.

### 4.4 GPU execution — 🔴

Schema states `gpu は未提供` ("not provided"). Empirically:
- 3D path rejects `backend:"gpu"` at build time (source: `crates/lbm-scenario/src/lib.rs:575`,
  message "compute.backend \"gpu\"（cpu / auto を指定してください）").
- 2D path **silently ignores** `backend:"gpu"` and runs on CPU — a 2D `backend:"gpu"` scenario
  returned `status:"completed"` with normal MLUPS (i.e. no GPU ran; no error surfaced).

  **Surprise / minor inconsistency worth flagging to PM:** `lbm validate` on the 2D gpu scenario
  returned `{ "ok": true, "warnings": [] }` (no warning), yet the schema says gpu is unavailable.
  There is no user-visible signal that gpu silently fell back to CPU in 2D. A user could believe a
  GPU run happened. (No source change made — reported only.)

No GPU numerical execution exists → 🔴. Any "run on GPU" user request is not satisfiable.
*(Superseded 2026-07-07: 2D GPU scenario dispatch landed with explicit errors — see status note at top.)*

### 4.5 Quantitative validation / comparison against reference — 🔴

There is **no user-facing** compare/analytic-check command: `grep` for
`compare|reference|analytic|ghia|golden|rmse|l2.?error` in `crates/lbm-cli/src` finds nothing in the
CLI/MCP surface. Reference/analytic validation (Ghia cavity, Poiseuille/Couette/TGV, cylinder Cd/St,
contact-angle, conservation) exists ONLY as developer `cargo test` targets in
`crates/lbm-core/tests/` (e.g. `validation_cavity.rs`, `validation_channel.rs`, `validation_tgv.rs`,
`validation_cylinder.rs`) — not reachable by a user via CLI or MCP. So the "validate/compare" stage
splits: **config validation = GREEN** (`lbm validate`, §2.7); **compare-to-reference = RED** for users.

### 4.6 Run monitoring beyond MCP async — 🔴 (partial)

Synchronous CLI `lbm run` blocks with no progress stream (only a final line). Live monitoring
(step-by-step progress, mid-run field snapshots on demand, cancellation) is available ONLY inside the
MCP async model (poll `run_status`), and even there status is coarse (running/completed/failed — the
`run_status` note carries elapsed seconds but not a step counter/percentage). No CLI progress bar,
no cancel command. The async *subsystem* is GREEN (§5); finer-grained monitoring is not present.

---

## 5. MCP inspection appendix (empirical — do not trust memory)

**Exact command driving the server** (throwaway script outside the repo tree, `/private/tmp/mcp_drive.sh`),
piping newline-delimited JSON-RPC into stdio:

```
{ printf '{"jsonrpc":"2.0","id":1,"method":"initialize","params":{...}}\n'
  printf '{"jsonrpc":"2.0","id":2,"method":"tools/list","params":{}}\n'
  printf '{"jsonrpc":"2.0","id":3,"method":"tools/call","params":{"name":"list_presets","arguments":{}}}\n'
  sleep 0.5
} | /Users/taku/projects/lbmflow-b1/target/release/lbm mcp
```

**`initialize` reply:**
```json
{"id":1,"jsonrpc":"2.0","result":{"capabilities":{"tools":{}},"protocolVersion":"2024-11-05",
 "serverInfo":{"name":"lbmflow","version":"0.1.0"}}}
```

**Tool count found empirically: 7** (CLAUDE.md's "7 tools incl. async start_run/run_status/list_runs"
is CONFIRMED). Full `tools/list` (names + input schemas):

| # | Tool | Required input | Purpose |
|---|------|----------------|---------|
| 1 | `run_scenario` | `scenario` (obj; req: name,grid,physics,edges,run) `[+outDir]` | Synchronous run; blocks; returns manifest. |
| 2 | `start_run` | `scenario` `[+outDir]` | **Async**: spawns background run, returns `runId` immediately. Cap 4 concurrent. |
| 3 | `run_status` | `runId` (string) | Returns `state: running\|completed\|failed`; manifest on completion. |
| 4 | `list_runs` | (none) | Enumerate all runs of this server in start order. |
| 5 | `validate_scenario` | `scenario` | Validate without running; config errors + stability warnings. |
| 6 | `list_presets` | (none) | All built-in presets with full scenario JSON. |
| 7 | `get_schema` | (none) | Full scenario JSON (v0) format reference (= `lbm schema`). |

Raw `tools/list` excerpt (input schema for the two async-critical tools):
```json
{"name":"start_run","inputSchema":{"type":"object",
  "properties":{"scenario":{"type":"object","required":["name","grid","physics","edges","run"]},
                "outDir":{"type":"string"}},"required":["scenario"]}}
{"name":"run_status","inputSchema":{"type":"object",
  "properties":{"runId":{"type":"string"}},"required":["runId"]}}
```

**Real tool call — async round-trip on a small scenario (start_run → list_runs → run_status):**

`start_run` reply (immediate, non-blocking, deterministic runId):
```json
{"id":2,"result":{"content":[{"text":"{ \"note\": \"バックグラウンドで実行中...\",
  \"outDir\": \".../async_out\", \"runId\": \"run-1-async-smoke\" }"}]}}
```
`list_runs` reply:
```json
{"id":3,"result":{"content":[{"text":"[ { \"outDir\": \".../async_out\",
  \"runId\": \"run-1-async-smoke\", \"scenarioName\": \"async-smoke\", \"state\": \"completed\" } ]"}]}}
```
`run_status` reply (carries manifest):
```json
{"id":4,"result":{"content":[{"text":"{ \"manifest\": { \"status\": \"completed\",
  \"stepsRun\": 300, \"mlups\": 31.2, \"diagnostics\": {...} },
  \"runId\": \"run-1-async-smoke\", \"state\": \"completed\" }"}]}}
```

**Verdict: the async MCP claim HOLDS.** `lbm mcp` starts cleanly over stdio, exposes exactly 7 tools,
and the async subsystem (`start_run`/`run_status`/`list_runs`, deterministic `run-<seq>-<name>` ids,
4-run concurrency cap in source) works end to end. MCP did NOT fail to start.

---

## 6. Top surprises (docs/memory claims vs. empirical `main`)

1. **Remote IS configured.** The brief and MEMORY claim no remote exists; `git remote -v` shows a
   live `origin` (github.com/TakuTsuzuki/lbmflow.git). Recorded in §0.
2. **GPU silently falls back to CPU in 2D.** Schema says gpu "未提供", 3D rejects it at build, but a
   2D `backend:"gpu"` scenario runs to `completed` on CPU with `validate` reporting `ok:true` and no
   warning — a user could wrongly believe a GPU run occurred (§4.4).
3. **"Validation" is developer-only.** Extensive analytic/reference validation exists but ONLY as
   `cargo test`; there is no user-facing compare/analytic command. The user-facing "validate" is
   config/stability checking, not quantitative accuracy comparison (§4.5).
4. **No unit→lattice conversion at all.** Every quantity (nu, u, rho, force) is raw lattice units;
   there is no SI/dimensional helper. The "BCs & properties (unit→lattice)" stage is BC-green but
   unit-conversion-red — a genuine gap where a future Skill would have to do dimensional analysis in
   prose, not by calling any engine feature (§1, §4 context).
5. **3D is more capable than expected but narrow.** Full D3Q19 3D runs (validate/run/VTK-volume/
   PNG-slice/CSV/force-fz) all work — but only single-phase, `init:rest`, CPU. 3D multiphase and 3D
   non-rest init are hard-rejected (§4.3).
6. **Tool count confirmed at 7** (not a surprise, but empirically verified against memory rather than
   trusted — CLAUDE.md was correct).

---

## 7. Handoff to PM (do NOT proceed to B2)

- Skill-eligible GREEN stages: scheme selection (bgk/trt only), BC/property specification, run/monitor
  via MCP async, post/visualize (PNG/CSV/VTK/gallery), config-validation, mesh/lattice (uniform
  Cartesian 2D/3D), and — as the one YELLOW — parametric-primitive geometry composition (`obstacles[]`).
- RED stages that must NOT get a Skill: CAD/mesh import, non-bgk/trt schemes, GPU execution, 3D
  multiphase, and quantitative validation-against-reference.
- Flag #2 and #4 above are the items most likely to mislead a user-facing Skill; B2 should decide how
  Skills communicate the CPU-only reality and the manual unit-conversion burden.
