#!/usr/bin/env python3
"""Compare stored LBMFlow and OpenLB benchmark fields.

This is a roadmap/operator harness for V&V lane 6.1. It does not launch
OpenLB. It reads already-produced fields, maps both sides to a common grid by
nearest neighbor, reports L2rel and Linf, and exits nonzero when bands fail.

Usage:
  python3 scripts/qa/openlb_compare.py manifest.json
  python3 scripts/qa/openlb_compare.py --self-test
  python3 scripts/qa/openlb_compare.py --check-openlb-build
"""

from __future__ import annotations

import argparse
import json
import math
import re
import struct
import subprocess
import sys
import tempfile
import zlib
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable
from xml.etree import ElementTree


OPENLB_ROOT = Path("/Users/taku/projects/cfd-bench")
REQUIRED_KEYS = {
    "benchmark",
    "lbmflow_output_dir",
    "openlb_output_dir",
    "grid_size",
    "band_L2rel",
    "band_linf",
}
BENCHMARKS = {"cavity", "cylinder", "tgv"}
EXTENSIONS = {".vtk", ".vti", ".csv", ".png"}


@dataclass(frozen=True)
class Field:
    values: list[float]
    dims: tuple[int, ...]
    source: Path
    label: str


def fail(message: str) -> None:
    raise SystemExit(f"error: {message}")


def parse_grid_size(raw: object, benchmark: str) -> tuple[int, ...]:
    if isinstance(raw, int):
        if raw <= 0:
            fail("grid_size must be positive")
        return (raw, raw)
    if isinstance(raw, list) and len(raw) in (2, 3):
        dims = tuple(int(v) for v in raw)
        if any(v <= 0 for v in dims):
            fail("grid_size entries must be positive")
        return dims
    fail(f"grid_size for {benchmark!r} must be an integer or [nx, ny[, nz]]")


def load_manifest(path: Path) -> dict:
    manifest = json.loads(path.read_text())
    missing = sorted(REQUIRED_KEYS - set(manifest))
    if missing:
        fail(f"{path}: missing required key(s): {', '.join(missing)}")
    benchmark = manifest["benchmark"]
    if benchmark not in BENCHMARKS:
        fail(f"{path}: benchmark must be one of {sorted(BENCHMARKS)}, got {benchmark!r}")
    manifest["grid_size"] = parse_grid_size(manifest["grid_size"], benchmark)
    for key in ("band_L2rel", "band_linf"):
        try:
            manifest[key] = float(manifest[key])
        except (TypeError, ValueError):
            fail(f"{path}: {key} must be numeric")
        if manifest[key] < 0.0:
            fail(f"{path}: {key} must be non-negative")
    return manifest


def normalize_name(name: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", name.lower())


def is_number_token(token: str) -> bool:
    try:
        float(token)
        return True
    except ValueError:
        return False


def file_score(path: Path, preferred_field: str | None) -> tuple[int, int, int, str]:
    name = path.stem.lower()
    field = (preferred_field or "speed").lower()
    field_score = 0
    if field and field in name:
        field_score = 4
    elif any(s in name for s in ("speed", "velocity", "vel", "u_")):
        field_score = 3
    elif any(s in name for s in ("ux", "uy", "uz", "rho")):
        field_score = 2
    step_match = re.search(r"_(\d+)$", path.stem)
    step = int(step_match.group(1)) if step_match else -1
    ext_score = {".vtk": 4, ".vti": 3, ".csv": 2, ".png": 1}.get(path.suffix.lower(), 0)
    return field_score, ext_score, step, path.name


def discover_field(output_dir: Path, solver: str, preferred_field: str | None) -> Field:
    if not output_dir.is_dir():
        fail(f"{solver} output dir does not exist: {output_dir}")

    manifest_files: list[Path] = []
    run_manifest = output_dir / "manifest.json"
    if run_manifest.is_file():
        try:
            files = json.loads(run_manifest.read_text()).get("files", [])
            manifest_files = [output_dir / f for f in files if (output_dir / f).suffix.lower() in EXTENSIONS]
        except (json.JSONDecodeError, AttributeError):
            manifest_files = []

    candidates = [p for p in manifest_files if p.is_file()]
    if not candidates:
        candidates = [p for p in output_dir.rglob("*") if p.is_file() and p.suffix.lower() in EXTENSIONS]
    if not candidates:
        fail(f"{solver} output dir has no supported field files: {output_dir}")

    errors: list[str] = []
    for path in sorted(candidates, key=lambda p: file_score(p, preferred_field), reverse=True):
        try:
            return read_field(path, preferred_field)
        except Exception as exc:
            errors.append(f"{path.name}: {type(exc).__name__}: {exc}")
    fail(f"no readable {solver} field in {output_dir}; tried: {'; '.join(errors[:8])}")


def read_field(path: Path, preferred_field: str | None) -> Field:
    ext = path.suffix.lower()
    if ext == ".vtk":
        return read_legacy_vtk(path, preferred_field)
    if ext == ".vti":
        return read_vti(path, preferred_field)
    if ext == ".csv":
        return read_csv_field(path, preferred_field)
    if ext == ".png":
        return read_png_luminance(path)
    raise ValueError(f"unsupported extension {ext}")


def read_legacy_vtk(path: Path, preferred_field: str | None) -> Field:
    lines = path.read_text(errors="replace").splitlines()
    dims: tuple[int, ...] | None = None
    point_data = None
    for i, line in enumerate(lines):
        parts = line.split()
        if len(parts) >= 4 and parts[0].upper() == "DIMENSIONS":
            raw = tuple(int(v) for v in parts[1:4])
            dims = raw if raw[2] != 1 else raw[:2]
        elif len(parts) >= 2 and parts[0].upper() == "POINT_DATA":
            point_data = i + 1
            break
    if dims is None:
        raise ValueError("missing DIMENSIONS")
    if point_data is None:
        raise ValueError("missing POINT_DATA")
    count = math.prod(dims)
    want = normalize_name(preferred_field or "")
    blocks: list[Field] = []
    i = point_data
    while i < len(lines):
        parts = lines[i].split()
        if not parts:
            i += 1
            continue
        kind = parts[0].upper()
        if kind == "SCALARS" and len(parts) >= 2:
            name = parts[1]
            comps = int(parts[3]) if len(parts) >= 4 and parts[3].isdigit() else 1
            i += 1
            if i < len(lines) and lines[i].strip().upper().startswith("LOOKUP_TABLE"):
                i += 1
            raw, i = collect_vtk_numbers(lines, i, count * comps)
            values = collapse_components(raw, comps)
            blocks.append(Field(values, dims, path, name))
            continue
        if kind == "VECTORS" and len(parts) >= 2:
            name = parts[1]
            raw, i = collect_vtk_numbers(lines, i + 1, count * 3)
            values = collapse_components(raw, 3)
            blocks.append(Field(values, dims, path, name))
            continue
        i += 1
    if not blocks:
        raise ValueError("no SCALARS or VECTORS point-data block")
    return choose_block(blocks, want)


def collect_vtk_numbers(lines: list[str], start: int, needed: int) -> tuple[list[float], int]:
    values: list[float] = []
    i = start
    block_headers = {"SCALARS", "VECTORS", "FIELD", "CELL_DATA", "POINT_DATA"}
    while i < len(lines) and len(values) < needed:
        parts = lines[i].split()
        if parts and parts[0].upper() in block_headers and values:
            break
        for token in parts:
            if is_number_token(token):
                values.append(float(token))
        i += 1
    if len(values) < needed:
        raise ValueError(f"VTK data ended after {len(values)} value(s), expected {needed}")
    return values[:needed], i


def collapse_components(raw: list[float], comps: int) -> list[float]:
    if comps == 1:
        return raw
    return [
        math.sqrt(sum(raw[i + c] * raw[i + c] for c in range(comps)))
        for i in range(0, len(raw), comps)
    ]


def choose_block(blocks: list[Field], want: str) -> Field:
    if want:
        for block in blocks:
            if normalize_name(block.label) == want or want in normalize_name(block.label):
                return block
    for fallback in ("speed", "velocity", "vel", "u"):
        for block in blocks:
            if fallback in normalize_name(block.label):
                return block
    return blocks[0]


def read_vti(path: Path, preferred_field: str | None) -> Field:
    root = ElementTree.parse(path).getroot()
    image = root.find(".//ImageData")
    piece = root.find(".//Piece")
    extent_text = None
    if piece is not None:
        extent_text = piece.attrib.get("Extent")
    if extent_text is None and image is not None:
        extent_text = image.attrib.get("WholeExtent")
    if extent_text is None:
        raise ValueError("missing ImageData/Piece extent")
    extent = [int(v) for v in extent_text.split()]
    raw_dims = (extent[1] - extent[0] + 1, extent[3] - extent[2] + 1, extent[5] - extent[4] + 1)
    dims = raw_dims if raw_dims[2] != 1 else raw_dims[:2]
    count = math.prod(dims)
    arrays = root.findall(".//PointData/DataArray") or root.findall(".//CellData/DataArray")
    blocks: list[Field] = []
    for arr in arrays:
        fmt = arr.attrib.get("format", "ascii").lower()
        if fmt != "ascii":
            raise ValueError(f"XML VTK DataArray format {fmt!r} is not supported")
        name = arr.attrib.get("Name", "field")
        comps = int(arr.attrib.get("NumberOfComponents", "1"))
        text = arr.text or ""
        raw = [float(t) for t in text.split()]
        if len(raw) < count * comps:
            raise ValueError(f"DataArray {name!r} has {len(raw)} value(s), expected {count * comps}")
        blocks.append(Field(collapse_components(raw[: count * comps], comps), dims, path, name))
    if not blocks:
        raise ValueError("no PointData/CellData DataArray")
    return choose_block(blocks, normalize_name(preferred_field or ""))


def read_csv_field(path: Path, preferred_field: str | None) -> Field:
    lines = [ln.strip() for ln in path.read_text(errors="replace").splitlines() if ln.strip()]
    if not lines:
        raise ValueError("empty CSV")

    first = lines[0]
    if first.startswith("#"):
        m = re.search(r"nx=(\d+),\s*ny=(\d+)(?:,\s*nz=(\d+))?", first)
        if m:
            dims = (int(m.group(1)), int(m.group(2)))
            if m.group(3):
                dims = dims + (int(m.group(3)),)
            values = [float(tok) for ln in lines[1:] for tok in split_row(ln)]
            if len(values) != math.prod(dims):
                raise ValueError(f"CSV has {len(values)} value(s), expected {math.prod(dims)}")
            return Field(values, dims, path, path.stem)
        lines = [ln for ln in lines if not ln.startswith("#")]

    rows = [split_row(ln) for ln in lines]
    if not rows:
        raise ValueError("CSV has no data rows")

    header = None
    if any(not is_number_token(tok) for tok in rows[0]):
        header = [normalize_name(tok) for tok in rows[0]]
        rows = rows[1:]
    numeric = [[float(tok) for tok in row] for row in rows if row]
    if not numeric:
        raise ValueError("CSV has no numeric rows")

    if header:
        return read_table_csv(path, header, numeric, preferred_field)
    if all(len(row) == len(numeric[0]) for row in numeric) and len(numeric) > 1:
        return Field([v for row in numeric for v in row], (len(numeric[0]), len(numeric)), path, path.stem)
    raise ValueError("headerless CSV must be a dense matrix")


def split_row(line: str) -> list[str]:
    return [tok for tok in re.split(r"[,\s]+", line.strip()) if tok]


def read_table_csv(path: Path, header: list[str], rows: list[list[float]], preferred_field: str | None) -> Field:
    coord_names = {"x", "ix", "i", "y", "iy", "j", "z", "iz", "k"}
    coord_idx = [i for i, name in enumerate(header) if name in coord_names]
    if len(coord_idx) < 2:
        raise ValueError("table CSV needs x/y coordinate columns")
    coord_idx = coord_idx[:3]
    value_idx = choose_csv_value_columns(header, preferred_field, coord_idx)
    coords = [[row[i] for i in coord_idx] for row in rows]
    axes = [sorted({int(c[d]) for c in coords}) for d in range(len(coord_idx))]
    dims = tuple(len(axis) for axis in axes)
    axis_pos = [{v: i for i, v in enumerate(axis)} for axis in axes]
    values = [math.nan] * math.prod(dims)
    for row, coord in zip(rows, coords):
        idx = tuple(axis_pos[d][int(coord[d])] for d in range(len(dims)))
        linear = idx[0] + dims[0] * idx[1] if len(dims) == 2 else idx[0] + dims[0] * (idx[1] + dims[1] * idx[2])
        comps = [row[i] for i in value_idx]
        values[linear] = comps[0] if len(comps) == 1 else math.sqrt(sum(v * v for v in comps))
    if any(math.isnan(v) for v in values):
        raise ValueError("table CSV does not cover a dense coordinate grid")
    label = ",".join(header[i] for i in value_idx)
    return Field(values, dims, path, label)


def choose_csv_value_columns(header: list[str], preferred_field: str | None, coord_idx: list[int]) -> list[int]:
    want = normalize_name(preferred_field or "")
    if want:
        exact = [i for i, name in enumerate(header) if i not in coord_idx and name == want]
        if exact:
            return exact[:1]
    vector_sets = [
        ("ux", "uy", "uz"),
        ("u0", "u1", "u2"),
        ("velocityx", "velocityy", "velocityz"),
        ("velx", "vely", "velz"),
    ]
    for names in vector_sets:
        idx = [header.index(n) for n in names if n in header and header.index(n) not in coord_idx]
        if len(idx) >= 2:
            return idx
    candidates = [i for i, _ in enumerate(header) if i not in coord_idx]
    if not candidates:
        raise ValueError("table CSV has no value columns")
    return candidates[-1:]


def read_png_luminance(path: Path) -> Field:
    data = path.read_bytes()
    if data[:8] != b"\x89PNG\r\n\x1a\n":
        raise ValueError("not a PNG")
    pos = 8
    width = height = color_type = bit_depth = interlace = None
    payload = b""
    while pos < len(data):
        length = struct.unpack(">I", data[pos : pos + 4])[0]
        kind = data[pos + 4 : pos + 8]
        chunk = data[pos + 8 : pos + 8 + length]
        pos += 12 + length
        if kind == b"IHDR":
            width, height, bit_depth, color_type, _compression, _filter, interlace = struct.unpack(">IIBBBBB", chunk)
        elif kind == b"IDAT":
            payload += chunk
        elif kind == b"IEND":
            break
    if width is None or height is None or color_type is None:
        raise ValueError("missing IHDR")
    if bit_depth != 8 or interlace != 0 or color_type not in (0, 2, 6):
        raise ValueError("only 8-bit non-interlaced grayscale/RGB/RGBA PNG is supported")
    channels = {0: 1, 2: 3, 6: 4}[color_type]
    raw = zlib.decompress(payload)
    stride = width * channels
    rows: list[bytes] = []
    prev = bytes(stride)
    p = 0
    for _ in range(height):
        filt = raw[p]
        p += 1
        row = bytearray(raw[p : p + stride])
        p += stride
        unfilter(row, prev, channels, filt)
        rows.append(bytes(row))
        prev = bytes(row)
    values = [0.0] * (width * height)
    for image_y, row in enumerate(rows):
        y = height - 1 - image_y
        for x in range(width):
            px = x * channels
            if channels == 1:
                lum = row[px] / 255.0
            else:
                lum = (0.2126 * row[px] + 0.7152 * row[px + 1] + 0.0722 * row[px + 2]) / 255.0
            values[y * width + x] = lum
    return Field(values, (width, height), path, "png-luminance")


def unfilter(row: bytearray, prev: bytes, channels: int, filt: int) -> None:
    for i in range(len(row)):
        left = row[i - channels] if i >= channels else 0
        up = prev[i]
        up_left = prev[i - channels] if i >= channels else 0
        if filt == 0:
            add = 0
        elif filt == 1:
            add = left
        elif filt == 2:
            add = up
        elif filt == 3:
            add = (left + up) // 2
        elif filt == 4:
            add = paeth(left, up, up_left)
        else:
            raise ValueError(f"bad PNG filter {filt}")
        row[i] = (row[i] + add) & 0xFF


def paeth(a: int, b: int, c: int) -> int:
    p = a + b - c
    pa = abs(p - a)
    pb = abs(p - b)
    pc = abs(p - c)
    if pa <= pb and pa <= pc:
        return a
    if pb <= pc:
        return b
    return c


def resample_nearest(field: Field, dims: tuple[int, ...]) -> list[float]:
    if field.dims == dims:
        return field.values
    if len(field.dims) != len(dims):
        raise ValueError(f"cannot compare {field.dims} field on {dims} common grid")
    src = field.dims
    out: list[float] = []
    for idx in iter_indices(dims):
        src_idx = tuple(nearest_index(idx[d], dims[d], src[d]) for d in range(len(dims)))
        out.append(field.values[linear_index(src_idx, src)])
    return out


def iter_indices(dims: tuple[int, ...]) -> Iterable[tuple[int, ...]]:
    if len(dims) == 2:
        nx, ny = dims
        for y in range(ny):
            for x in range(nx):
                yield (x, y)
    else:
        nx, ny, nz = dims
        for z in range(nz):
            for y in range(ny):
                for x in range(nx):
                    yield (x, y, z)


def nearest_index(i: int, n_dst: int, n_src: int) -> int:
    if n_dst == 1 or n_src == 1:
        return 0
    return min(n_src - 1, max(0, round(i * (n_src - 1) / (n_dst - 1))))


def linear_index(idx: tuple[int, ...], dims: tuple[int, ...]) -> int:
    if len(dims) == 2:
        x, y = idx
        return y * dims[0] + x
    x, y, z = idx
    return x + dims[0] * (y + dims[1] * z)


def l2rel(actual: list[float], reference: list[float]) -> float:
    num = sum((a - r) ** 2 for a, r in zip(actual, reference))
    den = sum(r * r for r in reference)
    if den == 0.0:
        return 0.0 if num == 0.0 else math.inf
    return math.sqrt(num / den)


def linf(actual: list[float], reference: list[float]) -> float:
    return max((abs(a - r) for a, r in zip(actual, reference)), default=0.0)


def compare(manifest: dict) -> int:
    preferred_field = manifest.get("field")
    common_dims = manifest["grid_size"]
    lbm = discover_field(Path(manifest["lbmflow_output_dir"]), "LBMFlow", preferred_field)
    openlb = discover_field(Path(manifest["openlb_output_dir"]), "OpenLB", preferred_field)
    lbm_values = resample_nearest(lbm, common_dims)
    openlb_values = resample_nearest(openlb, common_dims)
    l2 = l2rel(lbm_values, openlb_values)
    li = linf(lbm_values, openlb_values)
    ok_l2 = l2 <= manifest["band_L2rel"]
    ok_linf = li <= manifest["band_linf"]
    status = "PASS" if ok_l2 and ok_linf else "FAIL"
    print(f"benchmark={manifest['benchmark']} field={preferred_field or 'auto'} status={status}")
    print(f"common_grid={common_dims}")
    print(f"lbmflow={lbm.source} label={lbm.label} dims={lbm.dims}")
    print(f"openlb={openlb.source} label={openlb.label} dims={openlb.dims}")
    print(f"L2rel={l2:.17g} band={manifest['band_L2rel']:.17g} ok={ok_l2}")
    print(f"Linf={li:.17g} band={manifest['band_linf']:.17g} ok={ok_linf}")
    return 0 if status == "PASS" else 1


def write_synthetic_vtk(path: Path, dims: tuple[int, int], values: list[float]) -> None:
    nx, ny = dims
    lines = [
        "# vtk DataFile Version 3.0",
        "synthetic speed",
        "ASCII",
        "DATASET STRUCTURED_POINTS",
        f"DIMENSIONS {nx} {ny} 1",
        "ORIGIN 0 0 0",
        "SPACING 1 1 1",
        f"POINT_DATA {nx * ny}",
        "SCALARS speed double 1",
        "LOOKUP_TABLE default",
    ]
    for i in range(0, len(values), 8):
        lines.append(" ".join(str(v) for v in values[i : i + 8]))
    path.write_text("\n".join(lines) + "\n")


def self_test() -> int:
    with tempfile.TemporaryDirectory(prefix="openlb_compare_") as tmp:
        root = Path(tmp)
        lbm_dir = root / "lbm"
        olb_dir = root / "openlb"
        lbm_dir.mkdir()
        olb_dir.mkdir()
        dims = (4, 3)
        values = [float(x + 10 * y) for y in range(dims[1]) for x in range(dims[0])]
        write_synthetic_vtk(lbm_dir / "speed_0.vtk", dims, values)
        write_synthetic_vtk(olb_dir / "speed_0.vtk", dims, values)
        manifest_path = root / "compare.json"
        manifest_path.write_text(json.dumps({
            "benchmark": "tgv",
            "lbmflow_output_dir": str(lbm_dir),
            "openlb_output_dir": str(olb_dir),
            "grid_size": list(dims),
            "band_L2rel": 0.0,
            "band_linf": 0.0,
            "field": "speed",
        }) + "\n")
        rc = compare(load_manifest(manifest_path))
        if rc != 0:
            return rc
    print("self-test PASS")
    return 0


def check_openlb_build() -> int:
    ls = subprocess.run(["ls", str(OPENLB_ROOT)], text=True, capture_output=True)
    test = subprocess.run(["test", "-d", str(OPENLB_ROOT / "build")])
    print(f"ls {OPENLB_ROOT}: {'OK' if ls.returncode == 0 else 'FAIL'}")
    if ls.stdout.strip():
        print(ls.stdout.rstrip())
    if ls.stderr.strip():
        print(ls.stderr.rstrip(), file=sys.stderr)
    print(f"test -d {OPENLB_ROOT / 'build'}: {'OK' if test.returncode == 0 else 'FAIL'}")
    return 0 if ls.returncode == 0 and test.returncode == 0 else 1


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", nargs="?", help="comparison manifest JSON")
    parser.add_argument("--self-test", action="store_true", help="run synthetic identical-field test")
    parser.add_argument("--check-openlb-build", action="store_true", help="check cfd-bench root/build dirs only")
    args = parser.parse_args()

    modes = sum(bool(v) for v in (args.self_test, args.check_openlb_build, args.manifest))
    if modes != 1:
        parser.error("choose exactly one: MANIFEST, --self-test, or --check-openlb-build")
    if args.self_test:
        return self_test()
    if args.check_openlb_build:
        return check_openlb_build()
    return compare(load_manifest(Path(args.manifest)))


if __name__ == "__main__":
    sys.exit(main())
