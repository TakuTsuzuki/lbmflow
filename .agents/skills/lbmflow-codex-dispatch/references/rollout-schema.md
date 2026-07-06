# codex rollout jsonl — record schema reference

Sourced from real files under `~/.codex/sessions/2026/07/05/rollout-*.jsonl`
(codex CLI 0.142.4, model provider openai / gpt-5.5, `originator: codex_exec`).
Each line is one JSON object. Two top-level keys always present: `timestamp`
(ISO-8601 UTC) and `type`. The payload shape depends on `type`.

## Top-level record types (observed counts in one ~350-record order)

| `type` | count | role |
|---|---|---|
| `session_meta` | 1 (first line) | session header: `cwd`, `session_id`, `cli_version`, `base_instructions` |
| `turn_context` | 1 | per-turn context: `cwd`, `workspace_roots`, `approval_policy`, `sandbox_policy`, `current_date`, `timezone` |
| `event_msg` | many | high-level lifecycle events (see below) — the layer you monitor |
| `response_item` | many | raw model turn items (messages, reasoning, function_call, tool calls) |

## `event_msg` payload types (the monitoring layer)

`record.payload.type` is one of:

| payload.type | meaning | key fields |
|---|---|---|
| `task_started` | order began | `started_at` (epoch s), `model_context_window`, `collaboration_mode_kind` |
| `user_message` | the order text you dispatched | `message` (the full order prompt) |
| `agent_message` | codex narration / status | `message` (human-readable), `phase` (e.g. `"commentary"`) |
| `patch_apply_end` | a file edit was applied | `success` (bool), `stdout` (lists `M`/`A` files), `stderr`, `changes` |
| `token_count` | usage checkpoint | `info.total_token_usage.total_tokens`, `info.last_token_usage.*`, `info.model_context_window`, `rate_limits` |
| `task_complete` | **TERMINAL** — order finished | `last_agent_message` (final summary), `duration_ms`, `completed_at`, `time_to_first_token_ms`, `turn_id` |

### session_meta (first record)

```json
{
  "timestamp": "2026-07-05T11:52:28.223Z",
  "type": "session_meta",
  "payload": {
    "session_id": "019f321f-...",
    "cwd": "/Users/taku/projects/lbmflow-wt-cx-d4",
    "originator": "codex_exec",
    "cli_version": "0.142.4",
    "source": "exec",
    "model_provider": "openai",
    "base_instructions": { "text": "You are Codex, a coding agent..." }
  }
}
```

`payload.cwd` is the dispatch worktree — use it to match a rollout file to the
order you dispatched with `-C <worktree>`.

### task_started

```json
{ "type": "task_started", "turn_id": "019f...", "started_at": 1783177777,
  "model_context_window": 258400, "collaboration_mode_kind": "default" }
```

### patch_apply_end

```json
{ "type": "patch_apply_end", "call_id": "call_...", "turn_id": "019f...",
  "stdout": "Success. Updated the following files:\nM crates/lbm-core/tests/validation_conservation.rs\n",
  "stderr": "", "success": true, "changes": { "<abs path>": { ... } } }
```

### task_complete (terminal marker)

```json
{ "type": "task_complete", "turn_id": "019f...",
  "last_agent_message": "Implemented the revised tests ... Verification:\n- `cargo test --release`: failed as intended on red validation tests. ...",
  "completed_at": 1783179057, "duration_ms": 1279456, "time_to_first_token_ms": 3025 }
```

`last_agent_message` is codex's own end-of-run report. It typically lists changed
files, the verification commands it ran, and any tests it deliberately left red
(with the failing values). Read it before trusting the branch.

### token_count

```json
{ "type": "token_count",
  "info": {
    "total_token_usage": { "input_tokens": 2861552, "output_tokens": 26053, "total_tokens": 2887605 },
    "last_token_usage": { "total_tokens": 87829 },
    "model_context_window": 258400 },
  "rate_limits": { "primary": { "used_percent": 3.0, "window_minutes": 300 }, "secondary": { ... } } }
```

## `response_item` payload types (raw turn items — deeper detail)

| payload.type | notes |
|---|---|
| `message` | assistant/user message; `role`, `content[].text` |
| `reasoning` | model reasoning trace |
| `function_call` | tool invocation; `name` is `exec_command` or `write_stdin`, `arguments` is a JSON string (`cmd`, `workdir`, `yield_time_ms`, `max_output_tokens`) |
| `function_call_output` | the tool result; `output` (stdout chunk, wall time, session id) |
| `custom_tool_call` / `custom_tool_call_output` | non-shell tool calls (e.g. apply_patch) |

`exec_command` / `write_stdin` function_calls are how codex runs shell commands
and feeds long-running processes. A stream of these with growing
`function_call_output` (e.g. `test tXX ... ok`) is the signature of a live,
long-running test suite — **not** a hang.

## Copy-runnable probes

**Locate the rollout for a given worktree:**

```bash
DAY=$(date +%Y/%m/%d); ABS=$(cd ../lbmflow-wt-cx-<topic> && pwd)
for f in $(ls -t ~/.codex/sessions/$DAY/rollout-*.jsonl); do
  head -1 "$f" | grep -q "\"cwd\": \"$ABS\"" && { echo "$f"; break; }
done
```

**Done / running check + last summary:**

```bash
python3 - "$f" <<'PY'
import json,sys
recs=[json.loads(l) for l in open(sys.argv[1]) if l.strip()]
done=[r for r in recs if r.get("type")=="event_msg" and r["payload"].get("type")=="task_complete"]
if done:
    p=done[-1]["payload"]; print("DONE %.1fs"%(p["duration_ms"]/1000)); print(p["last_agent_message"])
else:
    print("RUNNING — last:", recs[-1]["payload"].get("type"))
PY
```

**Latest status narration + token usage:**

```bash
python3 - "$f" <<'PY'
import json,sys
recs=[json.loads(l) for l in open(sys.argv[1]) if l.strip()]
msgs=[r["payload"]["message"] for r in recs if r.get("type")=="event_msg" and r["payload"].get("type")=="agent_message"]
tc=[r for r in recs if r.get("type")=="event_msg" and r["payload"].get("type")=="token_count"]
if msgs: print("last status:", msgs[-1])
if tc: print("tokens:", tc[-1]["payload"]["info"]["total_token_usage"]["total_tokens"])
PY
```
