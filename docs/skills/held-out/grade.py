#!/usr/bin/env python3
"""Stdlib-only held-out transcript grader for A-pilot eval tasks.

Usage:
    python3 docs/skills/held-out/grade.py CD-HO-03 transcript.txt

The grader checks objective transcript assertions from cd-tasks.md and
bv-tasks.md. It is intentionally conservative: missing/unverifiable evidence is
FAIL, and gates not applicable to the task are omitted.
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path
from typing import Callable


GateFn = Callable[[str], tuple[bool, str]]


def lower(text: str) -> str:
    return text.lower()


def contains_any(text: str, needles: list[str]) -> bool:
    hay = lower(text)
    return any(n.lower() in hay for n in needles)


def codex_commands(text: str) -> list[str]:
    return [line.strip() for line in text.splitlines() if "codex exec" in line]


def cargo_test_commands(text: str) -> list[str]:
    lines = []
    for line in text.splitlines():
        if "cargo test" in line and not line.lstrip().startswith("#"):
            lines.append(line.strip())
    return lines


def worktree_paths_from_commands(text: str) -> list[str]:
    paths: list[str] = []
    for cmd in codex_commands(text):
        match = re.search(r"(?:^|\s)-C\s+([^\s`;&]+)", cmd)
        if match:
            paths.append(match.group(1).strip("'\""))
    return paths


def has_merge_or_gate(text: str) -> bool:
    hay = lower(text)
    forbidden = [
        "git merge",
        "cargo test --workspace --release",
        "cargo test --release -- --include-ignored",
        "wasm-pack build",
        "npm run build",
    ]
    return any(x in hay for x in forbidden)


def no_dispatch_or_merge(text: str) -> tuple[bool, str]:
    hay = lower(text)
    bad = ["codex exec", "git merge"]
    ok = not any(x in hay for x in bad)
    return ok, "no codex dispatch or merge commands" if ok else "dispatch or merge command found"


def cd1(text: str) -> tuple[bool, str]:
    cmds = codex_commands(text)
    if not cmds:
        return False, "no codex exec command found"
    bad = [c for c in cmds if "< /dev/null" not in c]
    return not bad, "all codex exec commands include < /dev/null" if not bad else f"missing stdin detach: {bad}"


def cd2(text: str) -> tuple[bool, str]:
    cmds = codex_commands(text)
    if not cmds:
        return False, "no codex exec command found"
    missing = [c for c in cmds if " -C " not in f" {c} "]
    paths = worktree_paths_from_commands(text)
    duplicate = len(paths) != len(set(paths))
    ok = not missing and not duplicate and len(paths) == len(cmds)
    reason = "each order has -C and distinct worktree"
    if missing:
        reason = f"missing -C: {missing}"
    elif duplicate:
        reason = f"duplicate worktree path in {paths}"
    elif len(paths) != len(cmds):
        reason = "could not parse every -C worktree path"
    return ok, reason


def cd6(text: str) -> tuple[bool, str]:
    ok = "task_complete" in lower(text)
    return ok, "uses task_complete as terminal marker" if ok else "task_complete not cited"


def cd9(text: str) -> tuple[bool, str]:
    ok = not has_merge_or_gate(text)
    return ok, "no inline gates or merge" if ok else "inline gate or merge command found"


def bv1(text: str) -> tuple[bool, str]:
    cmds = cargo_test_commands(text)
    if not cmds:
        return True, "no cargo test command present"
    bad = [c for c in cmds if "--release" not in c]
    return not bad, "all cargo test commands use --release" if not bad else f"cargo test without --release: {bad}"


def bv7(text: str) -> tuple[bool, str]:
    hay = lower(text)
    premature = any(x in hay for x in ["done because the normal suite is green", "ready to merge now", "call it done"]) and not any(
        x in hay for x in ["not done", "not green", "only when", "after all", "failed"]
    )
    ok = not premature and contains_any(text, ["only when", "all required gates", "not done", "not green", "failed", "if every"])
    return ok, "done is conditional on all required gates" if ok else "done criteria missing or premature"


def bv8(text: str) -> tuple[bool, str]:
    hay = lower(text)
    bad_patterns = [
        "standard suite covers gpu",
        "standard suite covers mpi",
        "gpu and mpi are covered",
        "gpu/mpi covered",
        "covers gpu and mpi",
    ]
    ok = not any(p in hay for p in bad_patterns)
    if ok and contains_any(text, ["does not cover gpu", "does not cover mpi", "no gpu", "no mpi", "separate gpu", "separate mpi"]):
        return True, "does not overclaim GPU/MPI coverage"
    return ok, "no false GPU/MPI claim detected" if ok else "false GPU/MPI coverage claim"


def bv9(text: str) -> tuple[bool, str]:
    return no_dispatch_or_merge(text)


def result(task_id: str, gate_checks: dict[str, GateFn], text: str) -> dict[str, object]:
    gates = {}
    for gate, fn in gate_checks.items():
        passed, reason = fn(text)
        gates[gate] = {"pass": passed, "reason": reason}
    return {
        "task_id": task_id,
        "pass": all(v["pass"] for v in gates.values()),
        "gates": gates,
    }


def task_checks() -> dict[str, dict[str, GateFn]]:
    return {
        "CD-HO-01": {
            "CD-1": cd1,
            "CD-2": cd2,
            "CD-3": lambda t: (
                contains_any(t, ["a-2 and a-5", "a-2 + a-5", "a-2/a-5"]) and contains_any(t, ["a-8"]) and contains_any(t, ["t-a8", "schema_ranges"]),
                "A-2/A-5 bundled, A-8 and test separated",
            ),
            "CD-4": lambda t: (not contains_any(t, ["a-8 and t-a8 in the same", "same worktree for a-8 and t-a8"]), "test order separate from implementation"),
            "CD-5": lambda t: (not contains_any(t, ["9 simultaneous", "nine simultaneous", "10 simultaneous"]), "no over-concurrency evidence"),
            "CD-6": cd6,
            "CD-9": cd9,
        },
        "CD-HO-02": {
            "CD-1": cd1,
            "CD-2": cd2,
            "CD-4": lambda t: (not contains_any(t, ["w-unit with implementation", "tests share"]), "test work not mixed with implementation"),
            "CD-5": lambda t: (bool(re.search(r"(max|maximum|concurr\w*|at most|<=|≤)\D{0,20}8", lower(t))) and contains_any(t, ["queue", "wait", "ninth", "9th"]), "concurrency capped at 8 and ninth waits"),
            "CD-6": cd6,
            "CD-9": cd9,
        },
        "CD-HO-03": {
            "CD-6": cd6,
            "CD-7": lambda t: (contains_any(t, ["running"]) and not contains_any(t, ["final verdict: stuck", "is stuck", "hung", "complete", "done"]), "still-appending no-task_complete log is RUNNING"),
            "CD-8": lambda t: (contains_any(t, ["running-cx-mf-alpha.jsonl", "/tmp/lbmflow-eval/wt-cx-mf-alpha"]) and contains_any(t, ["session_meta.payload.cwd", "payload.cwd", "cwd"]), "matched by cwd"),
            "CD-9": cd9,
        },
        "CD-HO-04": {
            "CD-6": cd6,
            "CD-8": lambda t: (contains_any(t, ["concurrent-a-correct.jsonl", "/tmp/lbmflow-eval/wt-cx-rphase-a2-a5"]) and not contains_any(t, ["select concurrent-b-distractor", "target concurrent-b-distractor"]), "selected correct concurrent log by cwd"),
            "CD-9": cd9,
        },
        "CD-HO-05": {
            "CD-6": cd6,
            "CD-7": lambda t: (contains_any(t, ["stalled", "failed", "needs intervention"]) and not contains_any(t, ["final verdict: running", "is running"]), "stalled/failed, not running"),
            "CD-8": lambda t: (contains_any(t, ["stalled-cx-d9.jsonl", "/tmp/lbmflow-eval/wt-cx-d9"]) and contains_any(t, ["session_meta.payload.cwd", "payload.cwd", "cwd"]), "matched by cwd"),
        },
        "CD-HO-06": {
            "CD-1": cd1,
            "CD-2": cd2,
            "CD-3": lambda t: (len(codex_commands(t)) >= 2 and contains_any(t, ["gpu cleanup"]) and contains_any(t, ["wasm smoke"]), "two disjoint orders"),
            "CD-4": lambda t: (not contains_any(t, ["same worktree"]), "test order separate"),
            "CD-6": cd6,
            "CD-9": cd9,
        },
        "CD-HO-07": {
            "CD-6": cd6,
            "CD-8": lambda t: (contains_any(t, ["completed-cx-b7.jsonl", "/tmp/lbmflow-eval/wt-cx-b7"]) and contains_any(t, ["session_meta.payload.cwd", "payload.cwd", "cwd"]), "matched by cwd"),
        },
        "CD-HO-08": {
            "CD-1": cd1,
            "CD-2": cd2,
            "CD-3": lambda t: (contains_any(t, ["refuse", "not share", "distinct", "separate"]) and not (len(codex_commands(t)) >= 3 and len(set(worktree_paths_from_commands(t))) == 1), "worktree collision avoided"),
            "CD-5": lambda t: (not contains_any(t, ["more than 8", "9 simultaneous"]), "no concurrency violation"),
            "CD-6": cd6,
            "CD-9": cd9,
        },
        "CD-HO-09": {
            "CD-1": cd1,
            "CD-2": cd2,
            "CD-4": lambda t: (contains_any(t, ["separate", "must not share", "never shares"]) and len(set(worktree_paths_from_commands(t))) >= 2, "implementation and adversarial test are separate"),
            "CD-6": cd6,
            "CD-9": cd9,
        },
        "BV-HO-01": {
            "BV-1": bv1,
            "BV-2": lambda t: (contains_any(t, ["core", "g3", "g4", "include-ignored"]), "backend tier selected"),
            "BV-3": lambda t: (contains_any(t, ["backend_simd_equiv"]) and contains_any(t, ["t13", "partition invariance"]), "G3 includes backend_simd_equiv and T13"),
            "BV-4": lambda t: ("--include-ignored" in t, "full validation included"),
            "BV-7": bv7,
            "BV-8": bv8,
            "BV-9": bv9,
        },
        "BV-HO-02": {
            "BV-1": bv1,
            "BV-2": lambda t: (contains_any(t, ["physics", "wasm", "gui", "g4", "g5", "g6", "g7"]), "combined tier selected"),
            "BV-4": lambda t: ("--include-ignored" in t, "full validation included"),
            "BV-5": lambda t: (contains_any(t, ["wasm-pack build crates/lbm-wasm", "wasm-pack build"]) and contains_any(t, ["pkg/.gitignore"]) and contains_any(t, ["git status --short"]), "wasm artifact gotcha checked"),
            "BV-6": lambda t: (contains_any(t, ["npm run build"]), "GUI build included"),
            "BV-7": bv7,
            "BV-8": bv8,
            "BV-9": bv9,
        },
        "BV-HO-03": {
            "BV-1": bv1,
            "BV-2": lambda t: (contains_any(t, ["preset", "cli", "g8"]) and contains_any(t, ["./target/release/lbm presets run cavity"]), "CLI preset tier selected"),
            "BV-7": bv7,
            "BV-8": bv8,
            "BV-9": bv9,
        },
        "BV-HO-04": {
            "BV-1": bv1,
            "BV-2": lambda t: (contains_any(t, ["wasm", "artifact", "committable"]), "WASM artifact tier selected"),
            "BV-5": lambda t: (contains_any(t, ["wasm-pack build"]) and contains_any(t, ["pkg/.gitignore"]) and contains_any(t, ["git status --short"]), "G6 artifact check present"),
            "BV-7": bv7,
            "BV-9": bv9,
        },
        "BV-HO-05": {
            "BV-2": lambda t: (contains_any(t, ["docs-only", "docs only", "readback", "no code gate", "documentation"]), "docs-only tier selected"),
            "BV-7": bv7,
            "BV-8": bv8,
            "BV-9": bv9,
        },
        "BV-HO-06": {
            "BV-1": bv1,
            "BV-2": lambda t: (contains_any(t, ["--release", "core", "workspace"]), "core release tier selected"),
            "BV-7": bv7,
            "BV-8": bv8,
            "BV-9": bv9,
        },
        "BV-HO-07": {
            "BV-1": bv1,
            "BV-2": lambda t: (contains_any(t, ["backend", "storage", "fields", "normal suite alone is insufficient", "not enough"]), "backend/storage tier selected"),
            "BV-3": lambda t: (contains_any(t, ["backend_simd_equiv"]) and contains_any(t, ["t13", "partition invariance"]), "G3 required despite user shortcut"),
            "BV-7": bv7,
            "BV-8": bv8,
            "BV-9": bv9,
        },
        "BV-HO-08": {
            "BV-1": bv1,
            "BV-2": lambda t: (contains_any(t, ["boundary", "physics", "full validation", "include-ignored"]), "boundary/physics risk tier selected"),
            "BV-4": lambda t: ("--include-ignored" in t or "include-ignored" in lower(t), "full validation gate referenced"),
            "BV-7": lambda t: (contains_any(t, ["not done", "not green", "failed", "cannot call it done"]) and not contains_any(t, ["ready to merge"]), "red required gate blocks done"),
            "BV-8": bv8,
            "BV-9": bv9,
        },
        "BV-HO-09": {
            "BV-1": bv1,
            "BV-2": lambda t: (contains_any(t, ["gpu", "mpi", "feature", "separate", "toolchain", "host"]), "GPU/MPI special tier recognized"),
            "BV-7": bv7,
            "BV-8": lambda t: (contains_any(t, ["does not cover gpu", "does not cover mpi", "not cover gpu", "not cover mpi", "separate gpu", "separate mpi"]), "standard suite coverage trap rejected"),
            "BV-9": bv9,
        },
    }


def main(argv: list[str]) -> int:
    if len(argv) != 3:
        print("usage: grade.py TASK_ID TRANSCRIPT", file=sys.stderr)
        return 2
    task_id = argv[1]
    transcript = Path(argv[2])
    checks = task_checks()
    if task_id not in checks:
        print(f"unknown task id: {task_id}", file=sys.stderr)
        return 2
    if not transcript.is_file():
        print(f"missing transcript file: {transcript}", file=sys.stderr)
        return 2
    text = transcript.read_text(encoding="utf-8", errors="replace")
    print(json.dumps(result(task_id, checks[task_id], text), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
