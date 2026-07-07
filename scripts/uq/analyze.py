#!/usr/bin/env python3
"""Analyze LBMFlow UQ sweep CSVs: OAT sensitivity and bootstrap CIs."""

from __future__ import annotations

import argparse
import csv
import math
import random
import re
from collections import defaultdict
from pathlib import Path
from typing import Any

try:
    import numpy as np  # type: ignore
except Exception:  # pragma: no cover - exercised on minimal Python installs
    np = None


EXCLUDE_QOI_SUFFIXES = {".status"}


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f))


def as_float(value: Any) -> float | None:
    try:
        if value is None or value == "":
            return None
        v = float(value)
        return v if math.isfinite(v) else None
    except (TypeError, ValueError):
        return None


def percentile(values: list[float], pct: float) -> float:
    if not values:
        return math.nan
    values = sorted(values)
    rank = (len(values) - 1) * pct / 100.0
    lo = math.floor(rank)
    hi = math.ceil(rank)
    if lo == hi:
        return values[lo]
    return values[lo] * (hi - rank) + values[hi] * (rank - lo)


def mean(values: list[float]) -> float:
    return sum(values) / len(values) if values else math.nan


def param_cols(rows: list[dict[str, str]]) -> list[str]:
    cols: list[str] = []
    for row in rows:
        for key in row:
            if key.startswith("param.") and key not in cols:
                cols.append(key)
    return cols


def qoi_cols(rows: list[dict[str, str]]) -> list[str]:
    cols: list[str] = []
    for row in rows:
        for key, value in row.items():
            if not key.startswith("qoi.") or key.endswith(tuple(EXCLUDE_QOI_SUFFIXES)):
                continue
            if key not in cols and as_float(value) is not None:
                cols.append(key)
    return cols


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def oat_sensitivity(
    rows: list[dict[str, str]], pcols: list[str], qcols: list[str]
) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    for pcol in pcols:
        other_params = [p for p in pcols if p != pcol]
        groups: dict[tuple[str, ...], list[dict[str, str]]] = defaultdict(list)
        for row in rows:
            groups[tuple(row.get(p, "") for p in other_params)].append(row)
        for qcol in qcols:
            slopes: list[float] = []
            normalized: list[float] = []
            group_count = 0
            for group_rows in groups.values():
                points: dict[float, list[float]] = defaultdict(list)
                for row in group_rows:
                    p = as_float(row.get(pcol))
                    q = as_float(row.get(qcol))
                    if p is None or q is None:
                        continue
                    points[p].append(q)
                if len(points) < 2:
                    continue
                xs = sorted(points)
                ys = [mean(points[x]) for x in xs]
                dx = xs[-1] - xs[0]
                if dx == 0:
                    continue
                dy = ys[-1] - ys[0]
                slope = dy / dx
                p_scale = max(abs(mean(xs)), 1e-300)
                q_scale = max(abs(mean(ys)), 1e-300)
                slopes.append(slope)
                normalized.append(slope * p_scale / q_scale)
                group_count += 1
            if slopes:
                out.append(
                    {
                        "parameter": pcol,
                        "qoi": qcol,
                        "groups": group_count,
                        "rawSlopeMean": mean(slopes),
                        "rawSlopeMin": min(slopes),
                        "rawSlopeMax": max(slopes),
                        "normalizedSlopeMean": mean(normalized),
                        "normalizedSlopeMin": min(normalized),
                        "normalizedSlopeMax": max(normalized),
                    }
                )
    return out


def bootstrap_means(values: list[float], samples: int, seed: int) -> list[float]:
    if np is not None:
        rng = np.random.default_rng(seed)
        arr = np.array(values, dtype=float)
        draws = rng.choice(arr, size=(samples, len(arr)), replace=True)
        return [float(x) for x in draws.mean(axis=1)]
    rng = random.Random(seed)
    return [mean([rng.choice(values) for _ in values]) for _ in range(samples)]


def bootstrap_ci(
    rows: list[dict[str, str]],
    pcols: list[str],
    qcols: list[str],
    samples: int,
    seed: int,
) -> list[dict[str, Any]]:
    grouped: dict[tuple[str, ...], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[tuple(row.get(p, "") for p in pcols)].append(row)
    out: list[dict[str, Any]] = []
    for key, group_rows in grouped.items():
        if len(group_rows) < 2:
            continue
        params = {p: key[i] for i, p in enumerate(pcols)}
        for qcol in qcols:
            values = [v for v in (as_float(r.get(qcol)) for r in group_rows) if v is not None]
            if len(values) < 2:
                continue
            draws = bootstrap_means(values, samples, seed)
            row: dict[str, Any] = {
                **params,
                "qoi": qcol,
                "n": len(values),
                "mean": mean(values),
                "ci95Low": percentile(draws, 2.5),
                "ci95High": percentile(draws, 97.5),
            }
            out.append(row)
    return out


def render_reports(
    rows: list[dict[str, str]],
    sensitivity: list[dict[str, Any]],
    ci_rows: list[dict[str, Any]],
    pcols: list[str],
    qcols: list[str],
    out_dir: Path,
) -> None:
    lines = [
        "LBMFlow UQ Sweep Summary",
        "========================",
        "",
        f"Runs: {len(rows)}",
        f"Parameters: {', '.join(pcols) if pcols else '(none)'}",
        f"Numeric QOIs: {len(qcols)}",
        f"Bootstrap CI rows with repeats: {len(ci_rows)}",
        "",
        "Top normalized OAT sensitivities by absolute mean:",
    ]
    top = sorted(
        sensitivity,
        key=lambda r: abs(float(r.get("normalizedSlopeMean") or 0.0)),
        reverse=True,
    )[:20]
    if top:
        for row in top:
            lines.append(
                "- {parameter} -> {qoi}: normalized={normalizedSlopeMean:.6g}, "
                "raw={rawSlopeMean:.6g}, groups={groups}".format(**row)
            )
    else:
        lines.append("- No one-at-a-time sensitivities could be computed.")
    lines += [
        "",
        "Bootstrap confidence intervals:",
        "- Computed only for exact repeated parameter points with n >= 2.",
        "- Rows are written to bootstrap_ci.csv; raw sweep rows are unchanged.",
    ]
    (out_dir / "summary.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")

    md = [
        "# LBMFlow UQ Sweep Summary",
        "",
        f"- Runs: {len(rows)}",
        f"- Parameters: {', '.join(f'`{p}`' for p in pcols) if pcols else '(none)'}",
        f"- Numeric QOIs: {len(qcols)}",
        f"- Bootstrap CI rows with repeats: {len(ci_rows)}",
        "",
        "## Top Normalized OAT Sensitivities",
        "",
    ]
    if top:
        md += ["| Parameter | QOI | Normalized slope | Raw slope | Groups |", "|---|---:|---:|---:|---:|"]
        for row in top:
            md.append(
                "| `{parameter}` | `{qoi}` | {normalizedSlopeMean:.6g} | "
                "{rawSlopeMean:.6g} | {groups} |".format(**row)
            )
    else:
        md.append("No one-at-a-time sensitivities could be computed.")
    md += [
        "",
        "## Bootstrap Confidence Intervals",
        "",
        "Computed only for exact repeated parameter points with `n >= 2`. "
        "See `bootstrap_ci.csv` for the full table.",
    ]
    (out_dir / "summary.md").write_text("\n".join(md) + "\n", encoding="utf-8")


def maybe_plot(
    rows: list[dict[str, str]], pcols: list[str], qcols: list[str], out_dir: Path
) -> None:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt  # type: ignore
    except Exception:
        return
    plot_dir = out_dir / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True)
    for pcol in pcols:
        xs = [as_float(r.get(pcol)) for r in rows]
        if not any(x is not None for x in xs):
            continue
        for qcol in qcols[:12]:
            pts = [
                (as_float(r.get(pcol)), as_float(r.get(qcol)))
                for r in rows
                if as_float(r.get(pcol)) is not None and as_float(r.get(qcol)) is not None
            ]
            if len(pts) < 2:
                continue
            fig, ax = plt.subplots()
            ax.scatter([p[0] for p in pts], [p[1] for p in pts])
            ax.set_xlabel(pcol)
            ax.set_ylabel(qcol)
            safe = re_safe(f"{pcol}_{qcol}")
            fig.tight_layout()
            fig.savefig(plot_dir / f"{safe}.png", dpi=160)
            plt.close(fig)


def re_safe(name: str) -> str:
    return "".join(c if c.isalnum() or c in "._-" else "_" for c in name)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("csv", type=Path, help="sweep CSV from sweep.py")
    ap.add_argument("--out-dir", type=Path, required=True, help="analysis output directory")
    ap.add_argument("--bootstrap-samples", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=12345)
    ap.add_argument("--plots", action="store_true", help="write scatter plots if matplotlib imports")
    args = ap.parse_args()

    rows = read_rows(args.csv)
    pcols = param_cols(rows)
    qcols = qoi_cols(rows)
    sens = oat_sensitivity(rows, pcols, qcols)
    cis = bootstrap_ci(rows, pcols, qcols, args.bootstrap_samples, args.seed)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_csv(
        args.out_dir / "sensitivity.csv",
        sens,
        [
            "parameter",
            "qoi",
            "groups",
            "rawSlopeMean",
            "rawSlopeMin",
            "rawSlopeMax",
            "normalizedSlopeMean",
            "normalizedSlopeMin",
            "normalizedSlopeMax",
        ],
    )
    write_csv(args.out_dir / "bootstrap_ci.csv", cis, [*pcols, "qoi", "n", "mean", "ci95Low", "ci95High"])
    render_reports(rows, sens, cis, pcols, qcols, args.out_dir)
    if args.plots:
        maybe_plot(rows, pcols, qcols, args.out_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
