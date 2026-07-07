#!/usr/bin/env python3
"""Compare native and WASM LBMFlow field snapshots.

Default band: 1e-4 L2rel for f32 WASM macroscopic fields compared with the
native f64 snapshot. Expected nonzero differences come from f32 arithmetic in
the browser/WASM path, tau/omega precision, and reconstruction of macroscopic
fields from the deviation-storage population representation.
"""

from __future__ import annotations

import argparse
import json
import math
import tempfile
from pathlib import Path
from typing import Any, Iterable, List, Sequence

DEFAULT_BAND = 1.0e-4
DEFAULT_FIELDS = "rho,velocity"


def _load_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as f:
        data = json.load(f)
    if not isinstance(data, dict):
        raise ValueError(f"{path}: top-level JSON value must be an object")
    return data


def _as_float_list(data: Any, label: str) -> List[float]:
    if not isinstance(data, list):
        raise ValueError(f"{label}: expected a JSON array")
    out = [float(v) for v in data]
    bad = [i for i, v in enumerate(out) if not math.isfinite(v)]
    if bad:
        raise ValueError(f"{label}: non-finite values at indices {bad[:5]}")
    return out


def _field(snapshot: dict[str, Any], name: str) -> List[float]:
    fields = snapshot.get("fields")
    if isinstance(fields, dict) and name in fields:
        return _as_float_list(fields[name], f"fields.{name}")
    if name in snapshot:
        return _as_float_list(snapshot[name], name)
    raise ValueError(f"snapshot is missing field '{name}'")


def _has_field(snapshot: dict[str, Any], name: str) -> bool:
    fields = snapshot.get("fields")
    return name in snapshot or (isinstance(fields, dict) and name in fields)


def _velocity(snapshot: dict[str, Any]) -> List[float]:
    ux = _field(snapshot, "ux")
    uy = _field(snapshot, "uy")
    if len(ux) != len(uy):
        raise ValueError(f"velocity component length mismatch: ux={len(ux)} uy={len(uy)}")
    uz: List[float] = []
    try:
        uz = _field(snapshot, "uz")
    except ValueError:
        pass
    if uz and len(uz) != len(ux):
        raise ValueError(f"velocity component length mismatch: ux={len(ux)} uz={len(uz)}")
    out: List[float] = []
    for i, (u, v) in enumerate(zip(ux, uy)):
        out.extend([u, v])
        if uz:
            out.append(uz[i])
    return out


def _mask(snapshot: dict[str, Any], include_solid: bool) -> List[bool] | None:
    if include_solid or "solid" not in snapshot:
        return None
    solid = snapshot["solid"]
    if not isinstance(solid, list):
        raise ValueError("solid: expected a JSON array")
    return [not bool(v) for v in solid]


def _apply_mask(values: Sequence[float], mask: Sequence[bool] | None, components: int) -> List[float]:
    if mask is None:
        return list(values)
    if len(values) != len(mask) * components:
        raise ValueError(
            f"mask length {len(mask)} is incompatible with value length {len(values)} "
            f"and component count {components}"
        )
    out: List[float] = []
    for cell, keep in enumerate(mask):
        if keep:
            start = cell * components
            out.extend(values[start:start + components])
    return out


def l2rel(actual: Sequence[float], reference: Sequence[float]) -> float:
    if len(actual) != len(reference):
        raise ValueError(f"length mismatch: actual={len(actual)} reference={len(reference)}")
    num = sum((a - r) * (a - r) for a, r in zip(actual, reference))
    den = sum(r * r for r in reference)
    if den == 0.0:
        return 0.0 if num == 0.0 else math.inf
    return math.sqrt(num / den)


def _compare_field(
    name: str,
    native: dict[str, Any],
    wasm: dict[str, Any],
    mask: Sequence[bool] | None,
) -> float:
    if name == "velocity":
        native_components = 3 if _has_field(native, "uz") else 2
        wasm_components = 3 if _has_field(wasm, "uz") else 2
        if native_components != wasm_components:
            raise ValueError(
                f"velocity dimensionality mismatch: native={native_components}D "
                f"wasm={wasm_components}D"
            )
        reference = _apply_mask(_velocity(native), mask, components=native_components)
        actual = _apply_mask(_velocity(wasm), mask, components=wasm_components)
    else:
        reference = _apply_mask(_field(native, name), mask, components=1)
        actual = _apply_mask(_field(wasm, name), mask, components=1)
    return l2rel(actual, reference)


def compare(
    native_path: Path,
    wasm_path: Path,
    fields: Iterable[str],
    band: float,
    include_solid: bool,
) -> int:
    native = _load_json(native_path)
    wasm = _load_json(wasm_path)

    for key in ("scenario_id", "nx", "ny", "steps"):
        if key in native and key in wasm and native[key] != wasm[key]:
            raise ValueError(f"{key} mismatch: native={native[key]!r} wasm={wasm[key]!r}")

    native_mask = _mask(native, include_solid)
    wasm_mask = _mask(wasm, include_solid)
    if native_mask is not None and wasm_mask is not None and native_mask != wasm_mask:
        raise ValueError("solid mask mismatch between native and WASM snapshots")

    print(f"native={native_path}")
    print(f"wasm={wasm_path}")
    print(f"band={band:.3e}")

    failed = False
    for field in fields:
        value = _compare_field(field, native, wasm, native_mask)
        status = "PASS" if value <= band else "FAIL"
        print(f"{status} {field}: L2rel={value:.9e}")
        failed = failed or value > band
    return 1 if failed else 0


def _self_test() -> None:
    native = {
        "schema_version": 1,
        "scenario_id": "self-test",
        "nx": 2,
        "ny": 2,
        "rho": [1.0, 1.0, 1.0, 1.0],
        "ux": [0.0, 0.1, 0.2, 0.3],
        "uy": [0.0, -0.1, -0.2, -0.3],
        "solid": [True, False, False, False],
    }
    wasm = {
        **native,
        "rho": [1.0, 1.0 + 2.0e-7, 1.0 - 1.0e-7, 1.0],
        "ux": [0.0, 0.1000002, 0.1999999, 0.3000001],
        "uy": [0.0, -0.0999999, -0.2000001, -0.2999999],
    }
    with tempfile.TemporaryDirectory() as td:
        native_path = Path(td) / "native.json"
        wasm_path = Path(td) / "wasm.json"
        native_path.write_text(json.dumps(native), encoding="utf-8")
        wasm_path.write_text(json.dumps(wasm), encoding="utf-8")
        rc = compare(native_path, wasm_path, ["rho", "velocity"], 1.0e-4, include_solid=False)
    if rc != 0:
        raise AssertionError("self-test comparison unexpectedly failed")
    print("wasm_parity_check.py self-test PASS")


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compare native and WASM JSON field snapshots with L2rel."
    )
    parser.add_argument("paths", nargs="*", help="optional positional paths: NATIVE_JSON WASM_JSON")
    parser.add_argument("--native", type=Path, help="native f64 JSON snapshot")
    parser.add_argument("--wasm", type=Path, help="WASM f32 JSON snapshot")
    parser.add_argument("--band", type=float, default=DEFAULT_BAND, help="maximum allowed L2rel")
    parser.add_argument(
        "--fields",
        default=DEFAULT_FIELDS,
        help="comma-separated fields to compare; use velocity for ux/uy[/uz] combined",
    )
    parser.add_argument(
        "--include-solid",
        action="store_true",
        help="compare all cells instead of masking out solid cells from native solid[]",
    )
    parser.add_argument("--self-test", action="store_true", help="run built-in parser/metric test")
    return parser.parse_args()


def main() -> int:
    args = _parse_args()
    if args.self_test:
        _self_test()
        return 0

    native = args.native
    wasm = args.wasm
    if args.paths:
        if len(args.paths) != 2:
            raise SystemExit("expected exactly two positional paths: NATIVE_JSON WASM_JSON")
        if native is not None or wasm is not None:
            raise SystemExit("use either positional paths or --native/--wasm, not both")
        native, wasm = Path(args.paths[0]), Path(args.paths[1])
    if native is None or wasm is None:
        raise SystemExit("provide NATIVE_JSON and WASM_JSON, or use --native and --wasm")

    fields = [f.strip() for f in args.fields.split(",") if f.strip()]
    if not fields:
        raise SystemExit("--fields must name at least one field")
    return compare(native, wasm, fields, args.band, args.include_solid)


if __name__ == "__main__":
    raise SystemExit(main())
