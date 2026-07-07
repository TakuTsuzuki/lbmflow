#!/usr/bin/env python3
"""Visual anomaly trawl for LBMFlow field artifacts.

Scans legacy ASCII VTK and PNG field outputs, computes file-level summary
statistics, and reports files whose metrics violate configured anomaly bands.
This is a postprocessing harness only; it does not run simulations.
"""

from __future__ import annotations

import argparse
import binascii
import json
import math
import re
import struct
import sys
import tempfile
import zlib
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable


DEFAULT_THRESHOLDS: dict[str, Any] = {
    "max_speed_hard": 0.3,
    "mass_drift_rel": 1e-6,
    "checkerboard_factor": 0.5,
    "u0": 0.1,
    "symmetry_rel": 1e-6,
    "symmetry_abs_floor": 1e-12,
    "png_value_scale": 1.0,
    "mass_fields": ["rho", "density", "mass"],
}

SUPPORTED_SYMMETRIES = {"mirror-x", "mirror-y", "rot180", "none", ""}


@dataclass
class FieldRecord:
    path: Path
    source_kind: str
    field_name: str
    step: int | None
    nx: int
    ny: int
    nz: int
    values: list[float]
    is_velocity: bool = False
    metrics: dict[str, float | None] = field(default_factory=dict)
    flags: list[dict[str, str]] = field(default_factory=list)


def parse_step(path: Path) -> int | None:
    match = re.search(r"_(\d+)(?:\.[^.]+)$", path.name)
    return int(match.group(1)) if match else None


def infer_field_from_name(path: Path) -> str:
    stem = path.stem
    match = re.match(r"(.+?)_(\d+)$", stem)
    return match.group(1) if match else stem


def parse_vtk(path: Path) -> FieldRecord:
    lines = path.read_text().splitlines()
    nx = ny = nz = None
    point_count = None
    for line in lines:
        parts = line.split()
        if len(parts) == 4 and parts[0] == "DIMENSIONS":
            nx, ny, nz = int(parts[1]), int(parts[2]), int(parts[3])
        elif len(parts) == 2 and parts[0] == "POINT_DATA":
            point_count = int(parts[1])

    if nx is None or ny is None or nz is None:
        raise ValueError(f"{path}: missing DIMENSIONS")
    expected = nx * ny * nz
    if point_count is not None and point_count != expected:
        raise ValueError(f"{path}: POINT_DATA {point_count} != {expected}")

    for idx, line in enumerate(lines):
        parts = line.split()
        if not parts:
            continue
        if parts[0] == "VECTORS" and len(parts) >= 3:
            field_name = parts[1]
            raw = _float_tokens(lines[idx + 1 :], expected * 3)
            magnitudes = [
                math.sqrt(raw[i] * raw[i] + raw[i + 1] * raw[i + 1] + raw[i + 2] * raw[i + 2])
                for i in range(0, len(raw), 3)
            ]
            return FieldRecord(
                path=path,
                source_kind="vtk-vector",
                field_name=field_name,
                step=parse_step(path),
                nx=nx,
                ny=ny,
                nz=nz,
                values=magnitudes,
                is_velocity=True,
            )
        if parts[0] == "SCALARS" and len(parts) >= 3:
            field_name = parts[1]
            components = int(parts[3]) if len(parts) >= 4 else 1
            data_start = None
            for j in range(idx + 1, len(lines)):
                if lines[j].startswith("LOOKUP_TABLE"):
                    data_start = j + 1
                    break
            if data_start is None:
                raise ValueError(f"{path}: SCALARS {field_name} missing LOOKUP_TABLE")
            raw = _float_tokens(lines[data_start:], expected * components)
            if components == 1:
                values = raw
            else:
                values = []
                for i in range(0, len(raw), components):
                    values.append(math.sqrt(sum(v * v for v in raw[i : i + components])))
            return FieldRecord(
                path=path,
                source_kind="vtk-scalar",
                field_name=field_name,
                step=parse_step(path),
                nx=nx,
                ny=ny,
                nz=nz,
                values=values,
                is_velocity=_is_speed_field(field_name, path),
            )

    raise ValueError(f"{path}: no VECTORS or SCALARS data array found")


def _float_tokens(lines: Iterable[str], limit: int) -> list[float]:
    vals: list[float] = []
    for line in lines:
        stripped = line.strip()
        if not stripped:
            continue
        vals.extend(float(tok) for tok in stripped.split())
        if len(vals) >= limit:
            break
    if len(vals) != limit:
        raise ValueError(f"expected {limit} values, got {len(vals)}")
    return vals


def parse_png(path: Path, png_value_scale: float) -> FieldRecord:
    width, height, gray = read_png_luminance(path)
    values = [v * png_value_scale for v in gray]
    return FieldRecord(
        path=path,
        source_kind="png",
        field_name=infer_field_from_name(path),
        step=parse_step(path),
        nx=width,
        ny=height,
        nz=1,
        values=values,
        is_velocity=_is_speed_field(infer_field_from_name(path), path),
    )


def read_png_luminance(path: Path) -> tuple[int, int, list[float]]:
    data = path.read_bytes()
    if not data.startswith(b"\x89PNG\r\n\x1a\n"):
        raise ValueError(f"{path}: not a PNG")
    pos = 8
    width = height = bit_depth = color_type = interlace = None
    palette: list[tuple[int, int, int]] = []
    idat_parts: list[bytes] = []
    while pos < len(data):
        if pos + 8 > len(data):
            raise ValueError(f"{path}: truncated PNG chunk")
        length = struct.unpack(">I", data[pos : pos + 4])[0]
        ctype = data[pos + 4 : pos + 8]
        payload = data[pos + 8 : pos + 8 + length]
        pos += 12 + length
        if ctype == b"IHDR":
            width, height, bit_depth, color_type, _comp, _filter, interlace = struct.unpack(
                ">IIBBBBB", payload
            )
        elif ctype == b"PLTE":
            palette = [
                (payload[i], payload[i + 1], payload[i + 2])
                for i in range(0, len(payload), 3)
            ]
        elif ctype == b"IDAT":
            idat_parts.append(payload)
        elif ctype == b"IEND":
            break

    if width is None or height is None or bit_depth is None or color_type is None:
        raise ValueError(f"{path}: missing IHDR")
    if interlace != 0:
        raise ValueError(f"{path}: interlaced PNG is not supported")
    if bit_depth not in (8, 16):
        raise ValueError(f"{path}: only 8-bit and 16-bit PNGs are supported")

    channels = {0: 1, 2: 3, 3: 1, 4: 2, 6: 4}.get(color_type)
    if channels is None:
        raise ValueError(f"{path}: unsupported PNG color type {color_type}")
    if color_type == 3 and bit_depth != 8:
        raise ValueError(f"{path}: paletted PNGs require 8-bit indices")
    bpp = max(1, channels * bit_depth // 8)
    row_bytes = width * bpp
    raw = zlib.decompress(b"".join(idat_parts))
    rows = _png_unfilter(raw, width, height, row_bytes, bpp)

    scale = float((1 << bit_depth) - 1)
    out: list[float] = []
    for row in rows:
        for x in range(width):
            off = x * bpp
            if color_type == 0:
                lum = _sample(row, off, bit_depth) / scale
            elif color_type == 3:
                idx = row[off]
                if idx >= len(palette):
                    raise ValueError(f"{path}: palette index {idx} out of range")
                r, g, b = palette[idx]
                lum = (0.2126 * r + 0.7152 * g + 0.0722 * b) / 255.0
            elif color_type == 4:
                lum = _sample(row, off, bit_depth) / scale
            else:
                r = _sample(row, off, bit_depth) / scale
                g = _sample(row, off + bit_depth // 8, bit_depth) / scale
                b = _sample(row, off + 2 * (bit_depth // 8), bit_depth) / scale
                lum = 0.2126 * r + 0.7152 * g + 0.0722 * b
            out.append(lum)
    return width, height, out


def _sample(row: bytes | bytearray, off: int, bit_depth: int) -> int:
    if bit_depth == 8:
        return row[off]
    return (row[off] << 8) | row[off + 1]


def _png_unfilter(raw: bytes, width: int, height: int, row_bytes: int, bpp: int) -> list[bytearray]:
    rows: list[bytearray] = []
    pos = 0
    prev = bytearray(row_bytes)
    for _y in range(height):
        if pos >= len(raw):
            raise ValueError("truncated PNG image data")
        filter_type = raw[pos]
        pos += 1
        row = bytearray(raw[pos : pos + row_bytes])
        pos += row_bytes
        if len(row) != row_bytes:
            raise ValueError("truncated PNG row")
        for i in range(row_bytes):
            left = row[i - bpp] if i >= bpp else 0
            up = prev[i]
            up_left = prev[i - bpp] if i >= bpp else 0
            if filter_type == 0:
                recon = row[i]
            elif filter_type == 1:
                recon = row[i] + left
            elif filter_type == 2:
                recon = row[i] + up
            elif filter_type == 3:
                recon = row[i] + ((left + up) // 2)
            elif filter_type == 4:
                recon = row[i] + _paeth(left, up, up_left)
            else:
                raise ValueError(f"unsupported PNG filter {filter_type}")
            row[i] = recon & 0xFF
        rows.append(row)
        prev = row
    if width * height == 0:
        raise ValueError("empty PNG")
    return rows


def _paeth(a: int, b: int, c: int) -> int:
    p = a + b - c
    pa = abs(p - a)
    pb = abs(p - b)
    pc = abs(p - c)
    if pa <= pb and pa <= pc:
        return a
    if pb <= pc:
        return b
    return c


def _is_speed_field(field_name: str, path: Path) -> bool:
    text = f"{field_name} {path.stem}".lower()
    return any(tok in text for tok in ("speed", "velocity", "vel", "u_mag", "umag"))


def _is_mass_field(field_name: str, thresholds: dict[str, Any]) -> bool:
    name = field_name.lower()
    return any(tok.lower() in name for tok in thresholds["mass_fields"])


def load_thresholds(path: Path | None, scenario_type: str) -> dict[str, Any]:
    thresholds = dict(DEFAULT_THRESHOLDS)
    if path is None:
        return thresholds
    loaded = json.loads(path.read_text())
    if not isinstance(loaded, dict):
        raise ValueError("threshold JSON must be an object")
    if any(k in loaded for k in DEFAULT_THRESHOLDS):
        thresholds.update(loaded)
    else:
        thresholds.update(loaded.get("default", {}))
        thresholds.update(loaded.get(scenario_type, {}))
    return thresholds


def apply_overrides(thresholds: dict[str, Any], args: argparse.Namespace) -> dict[str, Any]:
    out = dict(thresholds)
    for name in (
        "max_speed_hard",
        "mass_drift_rel",
        "checkerboard_factor",
        "u0",
        "symmetry_rel",
        "png_value_scale",
    ):
        value = getattr(args, name)
        if value is not None:
            out[name] = value
    return out


def scenario_hints(path: Path | None) -> tuple[str | None, str | None, dict[str, Any]]:
    if path is None:
        return None, None, {}
    data = json.loads(path.read_text())
    qa = data.get("qa", {}) if isinstance(data, dict) else {}
    vt = qa.get("visualTrawl", {}) if isinstance(qa, dict) else {}
    scenario_type = vt.get("scenarioType") or qa.get("scenarioType") or data.get("type")
    symmetry = vt.get("symmetry") or qa.get("symmetry") or data.get("symmetry")
    if isinstance(symmetry, dict):
        symmetry = symmetry.get("kind")
    thresholds = vt.get("thresholds", {}) if isinstance(vt, dict) else {}
    return scenario_type, symmetry, thresholds


def discover_files(input_dir: Path | None, recursive: bool) -> list[Path]:
    if input_dir is None:
        return []
    if input_dir.is_file():
        return [input_dir]
    pattern = "**/*" if recursive else "*"
    return sorted(
        p
        for p in input_dir.glob(pattern)
        if p.is_file() and p.suffix.lower() in {".vtk", ".png"}
    )


def discover_gallery(root: Path) -> list[Path]:
    return sorted(root.glob("*/vtk_field_*.vtk"))


def parse_records(files: Iterable[Path], thresholds: dict[str, Any]) -> tuple[list[FieldRecord], list[str]]:
    records: list[FieldRecord] = []
    errors: list[str] = []
    for path in files:
        try:
            if path.suffix.lower() == ".vtk":
                records.append(parse_vtk(path))
            elif path.suffix.lower() == ".png":
                records.append(parse_png(path, float(thresholds["png_value_scale"])))
        except Exception as exc:
            errors.append(f"{path}: {type(exc).__name__}: {exc}")
    return records, errors


def analyze_records(
    records: list[FieldRecord], thresholds: dict[str, Any], symmetry: str | None
) -> None:
    for rec in records:
        max_abs = max(abs(v) for v in rec.values)
        rec.metrics["max_abs_value"] = max_abs
        if rec.is_velocity:
            rec.metrics["max_speed"] = max_abs
            if max_abs > thresholds["max_speed_hard"]:
                _flag(
                    rec,
                    "max_speed_hard",
                    max_abs,
                    thresholds["max_speed_hard"],
                    "max|u| exceeds low-Mach hard limit",
                )
        else:
            rec.metrics["max_speed"] = None

        if rec.source_kind == "png":
            cb = checkerboard_mode_power(rec.values, rec.nx, rec.ny)
            rec.metrics["checkerboard_amplitude"] = cb
            cb_threshold = thresholds["checkerboard_factor"] * thresholds["u0"] * thresholds["u0"]
            if cb > cb_threshold:
                _flag(
                    rec,
                    "checkerboard_mode",
                    cb,
                    cb_threshold,
                    "Nyquist checkerboard-mode power exceeds init-ringing band",
                )
        else:
            rec.metrics["checkerboard_amplitude"] = None

    _apply_mass_drift(records, thresholds)

    if symmetry and symmetry not in {"", "none"}:
        if symmetry not in SUPPORTED_SYMMETRIES:
            raise ValueError(f"unsupported symmetry {symmetry!r}; expected {sorted(SUPPORTED_SYMMETRIES)}")
        for rec in records:
            if rec.nz != 1:
                rec.metrics["symmetry_violation"] = None
                continue
            violation = symmetry_violation(
                rec.values,
                rec.nx,
                rec.ny,
                symmetry,
                float(thresholds["symmetry_abs_floor"]),
            )
            rec.metrics["symmetry_violation"] = violation
            if violation > thresholds["symmetry_rel"]:
                _flag(
                    rec,
                    "symmetry_violation",
                    violation,
                    thresholds["symmetry_rel"],
                    f"{symmetry} relative violation",
                )
    else:
        for rec in records:
            rec.metrics["symmetry_violation"] = None


def _apply_mass_drift(records: list[FieldRecord], thresholds: dict[str, Any]) -> None:
    groups: dict[tuple[Path, str], list[FieldRecord]] = {}
    for rec in records:
        if _is_mass_field(rec.field_name, thresholds):
            groups.setdefault((rec.path.parent, rec.field_name), []).append(rec)
        else:
            rec.metrics["mass_drift_rel"] = None
    for group in groups.values():
        group.sort(key=lambda r: (r.step is None, r.step if r.step is not None else 0, r.path.name))
        baseline = sum(group[0].values)
        for rec in group:
            total = sum(rec.values)
            if baseline == 0.0:
                drift = 0.0 if total == 0.0 else math.inf
            else:
                drift = abs(total - baseline) / abs(baseline)
            rec.metrics["mass"] = total
            rec.metrics["mass_drift_rel"] = drift
            if rec is not group[0] and drift > thresholds["mass_drift_rel"]:
                _flag(
                    rec,
                    "mass_drift_rel",
                    drift,
                    thresholds["mass_drift_rel"],
                    "mass drift from first checkpoint exceeds band",
                )


def _flag(rec: FieldRecord, check: str, observed: float, threshold: float, detail: str) -> None:
    rec.flags.append(
        {
            "check": check,
            "observed": f"{observed:.6g}",
            "threshold": f"{threshold:.6g}",
            "detail": detail,
        }
    )


def checkerboard_mode_power(values: list[float], nx: int, ny: int) -> float:
    total = 0.0
    for y in range(ny):
        row = y * nx
        for x in range(nx):
            total += values[row + x] * (1.0 if ((x + y) % 2 == 0) else -1.0)
    coeff = abs(total) / float(nx * ny)
    return coeff * coeff


def symmetry_violation(values: list[float], nx: int, ny: int, symmetry: str, floor: float) -> float:
    max_abs = max(max(abs(v) for v in values), floor)
    worst = 0.0
    for y in range(ny):
        for x in range(nx):
            if symmetry == "mirror-x":
                xx, yy = nx - 1 - x, y
            elif symmetry == "mirror-y":
                xx, yy = x, ny - 1 - y
            elif symmetry == "rot180":
                xx, yy = nx - 1 - x, ny - 1 - y
            else:
                return 0.0
            worst = max(worst, abs(values[y * nx + x] - values[yy * nx + xx]))
    return worst / max_abs


def report(records: list[FieldRecord], errors: list[str], as_json: bool) -> str:
    flagged = [rec for rec in records if rec.flags]
    if as_json:
        payload = {
            "scanned": len(records),
            "parseErrors": errors,
            "flagged": [
                {
                    "path": str(rec.path),
                    "field": rec.field_name,
                    "step": rec.step,
                    "metrics": rec.metrics,
                    "flags": rec.flags,
                }
                for rec in flagged
            ],
            "files": [
                {
                    "path": str(rec.path),
                    "field": rec.field_name,
                    "step": rec.step,
                    "metrics": rec.metrics,
                    "flags": rec.flags,
                }
                for rec in records
            ],
        }
        return json.dumps(payload, indent=2, sort_keys=True)

    lines = [f"visual_trawl: scanned {len(records)} file(s), flagged {len(flagged)} file(s)"]
    if errors:
        lines.append("")
        lines.append("PARSE ERRORS:")
        lines.extend(f"  {err}" for err in errors)
    lines.append("")
    lines.append("FLAGGED FILES:")
    if not flagged:
        lines.append("  none")
    for rec in flagged:
        for flag in rec.flags:
            lines.append(
                f"  {rec.path}: {flag['check']} observed={flag['observed']} "
                f"threshold={flag['threshold']} ({flag['detail']})"
            )
    lines.append("")
    lines.append("SUMMARY:")
    for rec in records:
        metrics = ", ".join(
            f"{k}={_fmt_metric(v)}"
            for k, v in sorted(rec.metrics.items())
            if v is not None
        )
        flag_names = ",".join(flag["check"] for flag in rec.flags) or "-"
        lines.append(
            f"  {rec.path}: field={rec.field_name} step={rec.step if rec.step is not None else '-'} "
            f"kind={rec.source_kind} {metrics} flags={flag_names}"
        )
    return "\n".join(lines)


def _fmt_metric(value: float | None) -> str:
    if value is None:
        return "-"
    if not math.isfinite(value):
        return str(value)
    return f"{value:.6g}"


def run_scan(args: argparse.Namespace) -> tuple[int, str]:
    scenario_type_hint, symmetry_hint, scenario_thresholds = scenario_hints(args.scenario)
    scenario_type = args.scenario_type or scenario_type_hint or "default"
    thresholds = load_thresholds(args.thresholds, scenario_type)
    thresholds.update(scenario_thresholds)
    thresholds = apply_overrides(thresholds, args)

    symmetry = args.symmetry if args.symmetry is not None else symmetry_hint
    if symmetry is None:
        symmetry = "none"

    files = discover_files(args.input, args.recursive)
    if args.gallery is not None:
        files.extend(discover_gallery(args.gallery))
    files = sorted(set(files))
    records, errors = parse_records(files, thresholds)
    analyze_records(records, thresholds, symmetry)
    text = report(records, errors, args.json)
    if args.report is not None:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(text + "\n")
    return sum(len(rec.flags) for rec in records) + len(errors), text


def self_test() -> None:
    with tempfile.TemporaryDirectory(prefix="visual_trawl_selftest_") as td:
        root = Path(td)
        _write_scalar_vtk(root / "speed_0.vtk", "speed", 4, 4, [0.1] * 16)
        speed = [0.1] * 16
        speed[5] = 0.31
        _write_scalar_vtk(root / "speed_10.vtk", "speed", 4, 4, speed)
        _write_scalar_vtk(root / "rho_0.vtk", "rho", 4, 4, [1.0] * 16)
        rho = [1.0] * 16
        rho[0] = 1.001
        _write_scalar_vtk(root / "rho_10.vtk", "rho", 4, 4, rho)
        _write_png_gray8(root / "speed_20.png", 8, 8, [255 if (x + y) % 2 == 0 else 0 for y in range(8) for x in range(8)])

        args = argparse.Namespace(
            input=root,
            gallery=None,
            recursive=False,
            scenario=None,
            scenario_type="default",
            thresholds=None,
            max_speed_hard=None,
            mass_drift_rel=None,
            checkerboard_factor=None,
            u0=None,
            symmetry_rel=None,
            png_value_scale=None,
            symmetry="none",
            json=True,
            report=None,
        )
        code, text = run_scan(args)
        payload = json.loads(text)
        checks = {flag["check"] for rec in payload["flagged"] for flag in rec["flags"]}
        assert "max_speed_hard" in checks, checks
        assert "mass_drift_rel" in checks, checks
        assert "checkerboard_mode" in checks, checks
        assert code >= 3, code
    print("visual_trawl.py self-test PASS")


def _write_scalar_vtk(path: Path, field_name: str, nx: int, ny: int, values: list[float]) -> None:
    lines = [
        "# vtk DataFile Version 3.0",
        f"synthetic {field_name}",
        "ASCII",
        "DATASET STRUCTURED_POINTS",
        f"DIMENSIONS {nx} {ny} 1",
        "ORIGIN 0 0 0",
        "SPACING 1 1 1",
        f"POINT_DATA {nx * ny}",
        f"SCALARS {field_name} double 1",
        "LOOKUP_TABLE default",
    ]
    lines.extend(f"{v:.17g}" for v in values)
    path.write_text("\n".join(lines) + "\n")


def _write_png_gray8(path: Path, width: int, height: int, values: list[int]) -> None:
    rows = []
    for y in range(height):
        rows.append(bytes([0]) + bytes(values[y * width : (y + 1) * width]))
    raw = b"".join(rows)
    payload = struct.pack(">IIBBBBB", width, height, 8, 0, 0, 0, 0)
    chunks = [
        _png_chunk(b"IHDR", payload),
        _png_chunk(b"IDAT", zlib.compress(raw)),
        _png_chunk(b"IEND", b""),
    ]
    path.write_bytes(b"\x89PNG\r\n\x1a\n" + b"".join(chunks))


def _png_chunk(kind: bytes, payload: bytes) -> bytes:
    crc = binascii.crc32(kind)
    crc = binascii.crc32(payload, crc) & 0xFFFFFFFF
    return struct.pack(">I", len(payload)) + kind + payload + struct.pack(">I", crc)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Scan LBMFlow VTK/PNG field outputs for visual anomaly triggers."
    )
    parser.add_argument("input", nargs="?", type=Path, help="output directory or single .vtk/.png file")
    parser.add_argument("--gallery", type=Path, help="gallery root; scans gallery/*/vtk_field_*.vtk")
    parser.add_argument("--recursive", action="store_true", help="recurse under input directory")
    parser.add_argument("--scenario", type=Path, help="scenario JSON with optional qa.visualTrawl hints")
    parser.add_argument("--scenario-type", default=None, help="threshold profile name from --thresholds")
    parser.add_argument("--thresholds", type=Path, help="JSON threshold profile file")
    parser.add_argument("--max-speed-hard", type=float, default=None)
    parser.add_argument("--mass-drift-rel", type=float, default=None)
    parser.add_argument("--checkerboard-factor", type=float, default=None)
    parser.add_argument("--u0", type=float, default=None, help="driving speed for 0.5*u0^2 checkerboard band")
    parser.add_argument("--symmetry-rel", type=float, default=None)
    parser.add_argument("--png-value-scale", type=float, default=None, help="maps PNG luminance 0..1 to field units")
    parser.add_argument("--symmetry", choices=sorted(SUPPORTED_SYMMETRIES), default=None)
    parser.add_argument("--json", action="store_true", help="emit machine-readable JSON")
    parser.add_argument("--report", type=Path, help="write the same report to this path")
    parser.add_argument("--self-test", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.self_test:
        self_test()
        return 0
    if args.input is None and args.gallery is None:
        parser.error("provide an input directory/file, --gallery, or --self-test")
    failures, text = run_scan(args)
    print(text)
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
