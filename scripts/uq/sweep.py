#!/usr/bin/env python3
"""Run parameter sweeps through the LBMFlow CLI and collect tidy QOI rows."""

from __future__ import annotations

import argparse
import csv
import itertools
import json
import math
import os
import re
import subprocess
import sys
import time
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_LBM = ROOT / "target" / "release" / "lbm"
STATS = ("last", "mean", "std", "min", "max")


def load_json(path: Path) -> Any:
    with path.open("r", encoding="utf-8") as f:
        return json.load(f)


def write_json(path: Path, data: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, sort_keys=True)
        f.write("\n")


def path_parts(path: str) -> list[str | int]:
    parts: list[str | int] = []
    for token in path.split("."):
        if not token:
            raise ValueError(f"empty path segment in {path!r}")
        pos = 0
        m = re.match(r"^[^\[]+", token)
        if m:
            parts.append(m.group(0))
            pos = m.end()
        while pos < len(token):
            m = re.match(r"\[(\d+)\]", token[pos:])
            if not m:
                raise ValueError(f"cannot parse path segment {token!r} in {path!r}")
            parts.append(int(m.group(1)))
            pos += m.end()
    return parts


def set_path(obj: Any, path: str, value: Any) -> None:
    cur = obj
    parts = path_parts(path)
    for part in parts[:-1]:
        cur = cur[part]
    cur[parts[-1]] = value


def grouped_paths(key: str) -> list[str]:
    return [p.strip() for p in key.split(",") if p.strip()]


def spec_parameters(spec: dict[str, Any]) -> tuple[dict[str, list[Any]], int]:
    params = spec.get("parameters", spec)
    if not isinstance(params, dict):
        raise ValueError("sweep spec must be an object or contain a 'parameters' object")
    repeats = int(spec.get("repeats", 1)) if "parameters" in spec else 1
    if repeats < 1:
        raise ValueError("repeats must be >= 1")
    out: dict[str, list[Any]] = {}
    for key, values in params.items():
        if key in {"repeats", "metadata"} and "parameters" not in spec:
            continue
        if not isinstance(values, list) or not values:
            raise ValueError(f"parameter {key!r} must map to a non-empty list")
        out[key] = values
    return out, repeats


def expand_combinations(params: dict[str, list[Any]]) -> list[dict[str, Any]]:
    keys = list(params)
    combos: list[dict[str, Any]] = []
    for values in itertools.product(*(params[k] for k in keys)):
        combo: dict[str, Any] = {}
        for key, value in zip(keys, values):
            paths = grouped_paths(key)
            if len(paths) == 1:
                combo[paths[0]] = value
                continue
            if not isinstance(value, list) or len(value) != len(paths):
                raise ValueError(
                    f"grouped parameter {key!r} needs list values of length {len(paths)}"
                )
            for path, item in zip(paths, value):
                combo[path] = item
        combos.append(combo)
    return combos


def scenario_from_preset(lbm: Path, preset: str) -> dict[str, Any]:
    proc = subprocess.run(
        [str(lbm), "presets", "show", preset],
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    if proc.returncode != 0:
        raise RuntimeError(
            f"failed to read preset {preset!r}: {proc.stderr.strip() or proc.stdout.strip()}"
        )
    return json.loads(proc.stdout)


def flatten_manifest(manifest: dict[str, Any]) -> dict[str, Any]:
    diagnostics = manifest.get("diagnostics") or {}
    return {
        "qoi.manifest.status": manifest.get("status"),
        "qoi.manifest.stepsRun": manifest.get("stepsRun"),
        "qoi.manifest.wallSeconds": manifest.get("wallSeconds"),
        "qoi.manifest.mlups": manifest.get("mlups"),
        "qoi.manifest.totalMass": diagnostics.get("totalMass"),
        "qoi.manifest.maxSpeed": diagnostics.get("maxSpeed"),
        "qoi.manifest.tau": diagnostics.get("tau"),
        "qoi.manifest.warningCount": len(manifest.get("warnings") or []),
    }


def numeric_stats(values: list[float]) -> dict[str, float]:
    if not values:
        return {}
    n = len(values)
    mean = sum(values) / n
    var = sum((v - mean) ** 2 for v in values) / (n - 1) if n > 1 else 0.0
    return {
        "last": values[-1],
        "mean": mean,
        "std": math.sqrt(var),
        "min": min(values),
        "max": max(values),
    }


def read_probe_csv(path: Path) -> dict[str, Any]:
    qois: dict[str, Any] = {}
    with path.open("r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f)
        rows = list(reader)
    if not rows:
        return qois
    stem = path.stem
    prefix = "qoi.probe.force" if stem == "force" else f"qoi.probe.{stem}"
    qois[f"{prefix}.sampleCount"] = len(rows)
    for col in reader.fieldnames or []:
        if col == "step":
            continue
        vals: list[float] = []
        for row in rows:
            try:
                vals.append(float(row[col]))
            except (TypeError, ValueError):
                pass
        for stat, val in numeric_stats(vals).items():
            qois[f"{prefix}.{col}.{stat}"] = val
    return qois


def read_field_csv(path: Path) -> dict[str, Any]:
    values: list[float] = []
    with path.open("r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            for cell in line.split(","):
                try:
                    values.append(float(cell))
                except ValueError:
                    pass
    if not values:
        return {}
    field = path.stem.rsplit("_", 1)[0]
    prefix = f"qoi.field.{field}"
    stats = numeric_stats(values)
    rms = math.sqrt(sum(v * v for v in values) / len(values))
    qois: dict[str, Any] = {
        f"{prefix}.cellCount": len(values),
        f"{prefix}.rms": rms,
        f"{prefix}.sum": sum(values),
    }
    for stat, val in stats.items():
        qois[f"{prefix}.{stat}"] = val
    if field == "speed":
        qois[f"{prefix}.kineticEnergy"] = 0.5 * sum(v * v for v in values)
    return qois


def collect_csv_qois(out_dir: Path, manifest: dict[str, Any]) -> dict[str, Any]:
    qois: dict[str, Any] = {}
    for rel in manifest.get("files") or []:
        path = out_dir / rel
        if path.suffix.lower() != ".csv" or not path.is_file():
            continue
        if path.name == "force.csv" or path.name.startswith("point_"):
            qois.update(read_probe_csv(path))
        else:
            qois.update(read_field_csv(path))
    return qois


def stringify(value: Any) -> str:
    if isinstance(value, (dict, list)):
        return json.dumps(value, separators=(",", ":"), sort_keys=True)
    return "" if value is None else str(value)


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    keys: list[str] = []
    seen: set[str] = set()
    preferred_prefixes = ("run.", "param.", "qoi.")
    for prefix in preferred_prefixes:
        for row in rows:
            for key in row:
                if key.startswith(prefix) and key not in seen:
                    keys.append(key)
                    seen.add(key)
    for row in rows:
        for key in row:
            if key not in seen:
                keys.append(key)
                seen.add(key)
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=keys, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: stringify(row.get(k)) for k in keys})


def run_one(
    lbm: Path,
    scenario: dict[str, Any],
    params: dict[str, Any],
    work_dir: Path,
    combo_index: int,
    repeat: int,
    timeout: float | None,
) -> dict[str, Any]:
    scenario = json.loads(json.dumps(scenario))
    for path, value in params.items():
        set_path(scenario, path, value)
    safe_parts = [f"{p.replace('.', '_').replace(',', '_')}={stringify(v)}" for p, v in params.items()]
    run_id = f"run_{combo_index:04d}_rep_{repeat:02d}"
    scenario["name"] = f"{scenario.get('name', 'scenario')}_{run_id}"
    run_dir = work_dir / run_id
    scenario_path = run_dir / "scenario.json"
    write_json(scenario_path, scenario)
    row: dict[str, Any] = {
        "run.comboIndex": combo_index,
        "run.repeat": repeat,
        "run.id": run_id,
        "run.outputDir": str(run_dir),
        "run.scenarioPath": str(scenario_path),
        "run.parameterLabel": ";".join(safe_parts),
    }
    for path, value in params.items():
        row[f"param.{path}"] = value

    cmd = [str(lbm), "run", str(scenario_path), "--out", str(run_dir), "--json"]
    t0 = time.perf_counter()
    try:
        proc = subprocess.run(
            cmd,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=timeout,
            check=False,
        )
        row["run.returnCode"] = proc.returncode
        row["run.stdout"] = proc.stdout.strip()
        row["run.stderr"] = proc.stderr.strip()
    except subprocess.TimeoutExpired as e:
        row["run.returnCode"] = "timeout"
        row["run.stdout"] = (e.stdout or "").strip() if isinstance(e.stdout, str) else ""
        row["run.stderr"] = f"timeout after {timeout} seconds"
    row["run.wallSeconds"] = time.perf_counter() - t0

    manifest_path = run_dir / "manifest.json"
    if manifest_path.is_file():
        manifest = load_json(manifest_path)
        row.update(flatten_manifest(manifest))
        row.update(collect_csv_qois(run_dir, manifest))
    elif row.get("run.returnCode") == 0 and row.get("run.stdout"):
        try:
            manifest = json.loads(str(row["run.stdout"]))
            row.update(flatten_manifest(manifest))
            row.update(collect_csv_qois(run_dir, manifest))
        except json.JSONDecodeError:
            row["qoi.manifest.status"] = "missing-manifest"
    else:
        row["qoi.manifest.status"] = "run-failed"
    return row


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    source = ap.add_mutually_exclusive_group(required=True)
    source.add_argument("--base", type=Path, help="base scenario JSON")
    source.add_argument("--preset", help="built-in preset name, read via 'lbm presets show'")
    ap.add_argument("--spec", type=Path, required=True, help="sweep spec JSON")
    ap.add_argument("--out-csv", type=Path, required=True, help="tidy output CSV")
    ap.add_argument("--work-dir", type=Path, default=Path("out/uq"), help="run output root")
    ap.add_argument("--lbm", type=Path, default=DEFAULT_LBM, help="path to lbm binary")
    ap.add_argument("--timeout", type=float, default=None, help="per-run timeout in seconds")
    ap.add_argument("--fail-fast", action="store_true", help="stop after the first failed run")
    args = ap.parse_args()

    lbm = args.lbm
    if not lbm.is_file() or not os.access(lbm, os.X_OK):
        print(f"error: lbm binary is not executable: {lbm}", file=sys.stderr)
        return 2

    if args.base:
        base = load_json(args.base)
    else:
        base = scenario_from_preset(lbm, args.preset)

    params, repeats = spec_parameters(load_json(args.spec))
    combos = expand_combinations(params)
    args.work_dir.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, Any]] = []
    total = len(combos) * repeats
    ordinal = 0
    for combo_index, params_for_run in enumerate(combos):
        for repeat in range(repeats):
            ordinal += 1
            print(f"[{ordinal}/{total}] {params_for_run}", flush=True)
            row = run_one(
                lbm,
                base,
                params_for_run,
                args.work_dir,
                combo_index,
                repeat,
                args.timeout,
            )
            rows.append(row)
            write_csv(args.out_csv, rows)
            if args.fail_fast and row.get("qoi.manifest.status") in {
                "run-failed",
                "missing-manifest",
            }:
                return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
