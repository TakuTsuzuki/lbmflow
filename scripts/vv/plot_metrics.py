#!/usr/bin/env python3
"""Create lightweight SVG line plots from V&V CSV metrics.

This script intentionally uses only the Python standard library. It is meant for
V&V run packs such as out/vv/<run-id>/metrics.csv, but it also works with probe
CSVs like force.csv or torque.csv.
"""

import argparse
import csv
import math
from pathlib import Path


def _safe_name(name):
    return "".join(c if c.isalnum() or c in "._-" else "_" for c in name)


def _read_csv(path):
    if not Path(path).exists():
        raise SystemExit(f"{path}: file not found")
    with Path(path).open(newline="") as f:
        reader = csv.DictReader(f)
        if reader.fieldnames is None:
            raise SystemExit(f"{path}: missing CSV header")
        rows = list(reader)
    if not rows:
        raise SystemExit(f"{path}: no data rows")
    return reader.fieldnames, rows


def _as_float(value):
    try:
        v = float(value)
    except (TypeError, ValueError):
        return None
    return v if math.isfinite(v) else None


def _numeric_columns(fieldnames, rows):
    cols = []
    for name in fieldnames:
        values = [_as_float(row.get(name)) for row in rows]
        if any(v is not None for v in values):
            cols.append(name)
    return cols


def _series(rows, x_col, y_col):
    pts = []
    for i, row in enumerate(rows):
        x = _as_float(row.get(x_col)) if x_col else float(i)
        y = _as_float(row.get(y_col))
        if x is not None and y is not None:
            pts.append((x, y))
    return pts


def _scale(value, lo, hi, out_lo, out_hi):
    if hi == lo:
        return 0.5 * (out_lo + out_hi)
    return out_lo + (value - lo) * (out_hi - out_lo) / (hi - lo)


def _svg_polyline(points, title, x_label, y_label):
    width, height = 900, 520
    left, right, top, bottom = 80, 25, 45, 70
    xs = [p[0] for p in points]
    ys = [p[1] for p in points]
    xmin, xmax = min(xs), max(xs)
    ymin, ymax = min(ys), max(ys)
    if ymin == ymax:
        pad = abs(ymin) * 0.05 + 1.0
        ymin -= pad
        ymax += pad
    coords = []
    for x, y in points:
        px = _scale(x, xmin, xmax, left, width - right)
        py = _scale(y, ymin, ymax, height - bottom, top)
        coords.append(f"{px:.2f},{py:.2f}")
    x0, x1 = left, width - right
    y0, y1 = height - bottom, top
    return f"""<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">
  <rect width="100%" height="100%" fill="white"/>
  <text x="{left}" y="26" font-family="sans-serif" font-size="18">{title}</text>
  <line x1="{x0}" y1="{y0}" x2="{x1}" y2="{y0}" stroke="#222"/>
  <line x1="{x0}" y1="{y0}" x2="{x0}" y2="{y1}" stroke="#222"/>
  <polyline points="{' '.join(coords)}" fill="none" stroke="#2364aa" stroke-width="2"/>
  <text x="{left}" y="{height - 25}" font-family="sans-serif" font-size="13">{x_label}: {xmin:.6g} to {xmax:.6g}</text>
  <text x="{left}" y="{height - 8}" font-family="sans-serif" font-size="13">{y_label}: {ymin:.6g} to {ymax:.6g}</text>
</svg>
"""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("csv_path", help="metrics.csv, force.csv, torque.csv, or another headered CSV")
    ap.add_argument("--out", default=None, help="output directory for SVG plots")
    ap.add_argument("--x", default=None, help="x column; defaults to step if present, else first numeric column")
    ap.add_argument("--columns", default="", help="comma-separated y columns; defaults to all numeric columns except x")
    args = ap.parse_args()

    csv_path = Path(args.csv_path)
    fieldnames, rows = _read_csv(csv_path)
    numeric = _numeric_columns(fieldnames, rows)
    if not numeric:
        raise SystemExit(f"{csv_path}: no numeric columns")
    x_col = args.x or ("step" if "step" in numeric else numeric[0])
    if x_col not in numeric:
        raise SystemExit(f"{csv_path}: x column is not numeric: {x_col}")
    if args.columns:
        y_cols = [c.strip() for c in args.columns.split(",") if c.strip()]
    else:
        y_cols = [c for c in numeric if c != x_col]
    missing = [c for c in y_cols if c not in numeric]
    if missing:
        raise SystemExit(f"{csv_path}: non-numeric or missing columns: {', '.join(missing)}")
    out_dir = Path(args.out) if args.out else csv_path.parent / "plots"
    out_dir.mkdir(parents=True, exist_ok=True)

    written = []
    for y_col in y_cols:
        pts = _series(rows, x_col, y_col)
        if len(pts) < 2:
            continue
        title = f"{csv_path.name}: {y_col} vs {x_col}"
        svg = _svg_polyline(pts, title, x_col, y_col)
        out_path = out_dir / f"{csv_path.stem}_{_safe_name(y_col)}.svg"
        out_path.write_text(svg)
        written.append(out_path)
    if not written:
        raise SystemExit(f"{csv_path}: no plottable series with at least two points")
    for path in written:
        print(path)


if __name__ == "__main__":
    main()
