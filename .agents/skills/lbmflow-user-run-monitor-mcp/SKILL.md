---
name: lbmflow-user-run-monitor-mcp
description: >-
  Run LBMFlow scenarios asynchronously via its MCP server and monitor them to
  completion — the way to launch long or multiple simulations without blocking.
  Use whenever the user wants to "run this in the background", "kick off a long
  run and check on it", "run a sweep / several scenarios at once", "start a run
  and poll it", "use the MCP server", "connect an agent to LBMFlow", or asks how
  to check a run's status / list running simulations. This Skill owns the `lbm
  mcp` stdio server, its 7 tools, and the async loop (`start_run` → `run_status`
  → `list_runs`) including the 4-concurrent cap and result retrieval. Do NOT use
  it for a quick blocking CLI run of one preset (that is lbmflow-user-run-preset)
  or to author/validate the scenario JSON first (author-scenario / tune-stability
  produce the scenario this Skill runs).
---

# LBMFlow — async run & monitor via MCP

For a long simulation (tens of thousands of steps, a big grid) or several at
once, the blocking CLI `lbm run` is the wrong tool — it blocks until the run
finishes with no progress. LBMFlow's MCP server exposes an **async** model:
`start_run` returns a `runId` immediately and the run proceeds in the background;
you poll `run_status` for the result. This Skill covers starting the server,
launching async runs, and polling them to completion.

**Prerequisite:** the validated scenario JSON (from
`lbmflow-user-author-scenario` / `lbmflow-user-tune-stability`). MCP tools take
the scenario as an OBJECT, and `start_run`/`run_scenario` require it to contain
`name, grid, physics, edges, run`.

## Start the server

```bash
./target/release/lbm mcp        # speaks MCP over stdio (JSON-RPC, newline-delimited)
```

It advertises `serverInfo.name = "lbmflow"` and exposes exactly **7 tools**. In a
normal agent setup you register this command as an MCP server and call the tools;
the tool set is what matters.

## The 7 tools (exact — nothing else exists)

| Tool | Required input | What it does |
|---|---|---|
| `run_scenario` | `scenario` (obj) `[+ outDir]` | **Synchronous** run; blocks; returns the manifest. Use only for short runs. |
| `start_run` | `scenario` (obj) `[+ outDir]` | **Async** — spawns a background run, returns `{ runId }` immediately. Cap **4 concurrent**. |
| `run_status` | `runId` (string) | Returns `state: running \| completed \| failed`; the manifest once completed. |
| `list_runs` | (none) | Lists all runs of this server in start order, with `state`. |
| `validate_scenario` | `scenario` (obj) | Validate without running (config errors + stability warnings). |
| `list_presets` | (none) | All built-in presets with full scenario JSON. |
| `get_schema` | (none) | Full scenario JSON (v0) format reference (= `lbm schema`). |

## Decision procedure — sync or async

| Situation | Tool |
|---|---|
| One short run, you want the result inline now | `run_scenario` (blocks) |
| A long run, or you want to keep working while it runs | `start_run` → poll `run_status` |
| Several runs / a parameter sweep | multiple `start_run` (≤ 4 in flight), poll each |
| Just checking a config, not running | `validate_scenario` |

## The async loop (start → poll → fetch)

1. **(Optional) `validate_scenario`** first so a bad config fails fast, before you
   spend a background run on it.
2. **`start_run { scenario, outDir? }`** → returns immediately with a
   deterministic id, e.g. `runId: "run-1-<name>"`. The reply text notes the run
   is executing in the background and gives its `outDir`.
3. **Poll `run_status { runId }`** until `state` is `completed` or `failed`. While
   it is `running`, wait and poll again — do NOT treat `running` as stuck. The
   status note carries elapsed seconds (coarse: there is no step counter or
   percentage — see limits below).
4. **On `completed`**, `run_status` carries the full `manifest` (status, stepsRun,
   mlups, diagnostics, files). That is the result. On `failed`, report it.
5. **`list_runs`** any time to see every run and its state at a glance.

**Concurrency cap:** `start_run` allows **at most 4 concurrent** background runs.
For a sweep of more than 4, launch 4, poll for completions, then launch the next
as slots free. Do not fire all N at once.

## Verification gate — the done check

An async run is done only when `run_status` reports:

```json
{ "state": "completed",
  "runId": "run-1-<name>",
  "manifest": { "status": "completed", "stepsRun": <N>, "mlups": <..>,
                "diagnostics": { ... }, "files": [ ... ] } }
```

Both `state == "completed"` AND `manifest.status == "completed"` must hold. If
`state == "failed"`, or the manifest `status` is `"diverged"`/`"failed"`, the run
did NOT succeed — report the state and diagnostics, do not claim success. A run is
only finished when its state is terminal (`completed`/`failed`); `running` means
keep polling.

## Known limits (be honest about these)

- **Monitoring is coarse.** `run_status` gives `running/completed/failed` plus
  elapsed seconds — **no step counter, no percentage, no mid-run field snapshot,
  no cancel.** Do not promise a live progress bar or the ability to peek at fields
  mid-run; those do not exist. Report elapsed time and state only.
- **4 concurrent max** on `start_run`. Larger sweeps queue.
- GPU is not available; runs execute on CPU regardless of any `backend` hint.

## Worked example (end-to-end)

Task: "Run these three cavity variants in the background and tell me when each
finishes."

1. **Validate** each with `validate_scenario` → all `ok:true`.
2. **Launch:** three `start_run` calls (≤ 4 cap, fine) → `run-1-cavA`,
   `run-2-cavB`, `run-3-cavC`, each returned immediately.
3. **Poll:** `list_runs` shows all three `running`. Poll each `run_status` until
   `completed`. A long one stays `running` for a while — that is alive, not stuck.
4. **Fetch:** on each `completed`, read `manifest.files` and `diagnostics` from
   the `run_status` reply. Report per-run: status, stepsRun, mlups, output files.

## Top failure modes (and the fix)

- **Used `run_scenario` for a huge run and blocked everything.** Fix: use
  `start_run` for long runs; reserve `run_scenario` for short ones.
- **Called a tool that doesn't exist** (`cancel_run`, `run_progress`, `stop`).
  Only the 7 tools above exist. Fix: there is no cancel/progress tool; poll
  `run_status`.
- **Treated `state:"running"` as failure/stuck.** It means in progress. Fix: keep
  polling until terminal; report elapsed seconds meanwhile.
- **Fired more than 4 `start_run` at once for a sweep.** Fix: cap at 4 in flight;
  launch the rest as runs complete.
- **Claimed a live progress %/mid-run snapshot.** Not available. Fix: report only
  state + elapsed seconds.
- **Passed the scenario as a JSON string instead of an object.** MCP tools take
  `scenario` as an object with `name/grid/physics/edges/run`. Fix: pass the parsed
  object.
- **Reported `completed` state but ignored a `diverged` manifest.** Check BOTH
  `state` and `manifest.status`. Fix: a diverged manifest is not success.
