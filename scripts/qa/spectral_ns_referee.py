#!/usr/bin/env python3
"""Independent 2D spectral Navier-Stokes referee for Taylor-Green vortex.

This is intentionally separate from the Rust solver.  It solves periodic
incompressible Navier-Stokes on a square with a Fourier pseudo-spectral method:

  du/dt + P[(u.grad)u] = nu Laplacian(u),  div u = 0,

where P is the Fourier-space Leray projector.  The quadratic convolution is
evaluated with 3/2 zero-padding and the resulting nonlinear term is filtered
with the standard 2/3 cutoff (Orszag dealiasing).  The viscous linear term is
advanced exactly by using an integrating factor exp(-nu |K|^2 dt).

CLI comparison schema:
  [
    {"t": 0, "ux": [[...], ...], "uy": [[...], ...]},
    {"t": 100, "ux": [[...], ...], "uy": [[...], ...]}
  ]

`t` is interpreted in the same time units as nu.  Use --time-scale if the JSON
stores integer step counts that need conversion to physical time.

TODO(lane-6.2-follow-up): add the Rust-side route that exports LBMFlow TGV
snapshots in this schema.  This order deliberately does not touch Rust.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple

import numpy as np


ArrayPair = Tuple[np.ndarray, np.ndarray]


def _validate_inputs(N: int, nu: float, u0: float, k: float, t: float) -> None:
    if N < 8 or N % 2 != 0:
        raise ValueError("N must be an even integer >= 8")
    if nu <= 0.0:
        raise ValueError("nu must be positive")
    if u0 < 0.0:
        raise ValueError("u0 must be non-negative")
    if k <= 0.0:
        raise ValueError("k must be positive")
    if t < 0.0:
        raise ValueError("t must be non-negative")


def _wave_numbers(N: int, k: float) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    # The TGV parameter k is the magnitude of the fundamental two-component
    # wavevector.  Each square-axis component is therefore alpha = k/sqrt(2).
    alpha = k / math.sqrt(2.0)
    n = np.fft.fftfreq(N, d=1.0 / N)
    kx, ky = np.meshgrid(alpha * n, alpha * n, indexing="xy")
    k2 = kx * kx + ky * ky

    cutoff = math.floor(N / 3)
    nx, ny = np.meshgrid(n, n, indexing="xy")
    two_thirds = (np.abs(nx) <= cutoff) & (np.abs(ny) <= cutoff)
    return kx, ky, k2, two_thirds


def _initial_tgv(N: int, u0: float) -> ArrayPair:
    theta = 2.0 * math.pi * np.arange(N, dtype=np.float64) / N
    x, y = np.meshgrid(theta, theta, indexing="xy")

    # 2D Taylor-Green/Kolmogorov projected mode from streamfunction
    # psi = -u0 sin(x) sin(y).  It is divergence-free and has |K|^2 = k^2
    # under the physical scaling in _wave_numbers, so its kinetic energy
    # decays as E(t)/E(0) = exp(-2 nu k^2 t).
    ux = -u0 * np.cos(x) * np.sin(y)
    uy = u0 * np.sin(x) * np.cos(y)
    return ux, uy


def _pad_spectrum(a: np.ndarray, M: int) -> np.ndarray:
    N = a.shape[0]
    if a.shape != (N, N):
        raise ValueError("expected a square spectrum")
    if M < N:
        raise ValueError("padded size must be >= N")
    out = np.zeros((M, M), dtype=np.complex128)
    shifted = np.fft.fftshift(a)
    start = (M - N) // 2
    out[start:start + N, start:start + N] = shifted
    # numpy's ifft2 divides by grid_size^2.  Scaling preserves the same
    # physical field amplitude when evaluating the N-mode spectrum on M^2
    # collocation points.
    return np.fft.ifftshift(out) * ((M / N) ** 2)


def _truncate_spectrum(a: np.ndarray, N: int) -> np.ndarray:
    M = a.shape[0]
    if a.shape != (M, M):
        raise ValueError("expected a square spectrum")
    shifted = np.fft.fftshift(a)
    start = (M - N) // 2
    cropped = shifted[start:start + N, start:start + N]
    # Convert M-grid FFT normalization back to the N-grid normalization.
    return np.fft.ifftshift(cropped) * ((N / M) ** 2)


def _project(
    ax: np.ndarray,
    ay: np.ndarray,
    kx: np.ndarray,
    ky: np.ndarray,
    k2: np.ndarray,
) -> ArrayPair:
    out_x = ax.copy()
    out_y = ay.copy()
    nonzero = k2 > 0.0
    dot = kx[nonzero] * ax[nonzero] + ky[nonzero] * ay[nonzero]
    out_x[nonzero] = ax[nonzero] - kx[nonzero] * dot / k2[nonzero]
    out_y[nonzero] = ay[nonzero] - ky[nonzero] * dot / k2[nonzero]
    return out_x, out_y


def _nonlinear_rhs(
    ux_hat: np.ndarray,
    uy_hat: np.ndarray,
    kx: np.ndarray,
    ky: np.ndarray,
    k2: np.ndarray,
    two_thirds: np.ndarray,
) -> ArrayPair:
    """Return -P[(u.grad)u] in Fourier space."""
    N = ux_hat.shape[0]
    M = 3 * N // 2

    ux_m = np.fft.ifft2(_pad_spectrum(ux_hat, M)).real
    uy_m = np.fft.ifft2(_pad_spectrum(uy_hat, M)).real
    dux_dx = np.fft.ifft2(_pad_spectrum(1j * kx * ux_hat, M)).real
    dux_dy = np.fft.ifft2(_pad_spectrum(1j * ky * ux_hat, M)).real
    duy_dx = np.fft.ifft2(_pad_spectrum(1j * kx * uy_hat, M)).real
    duy_dy = np.fft.ifft2(_pad_spectrum(1j * ky * uy_hat, M)).real

    adv_x = ux_m * dux_dx + uy_m * dux_dy
    adv_y = ux_m * duy_dx + uy_m * duy_dy
    adv_x_hat = _truncate_spectrum(np.fft.fft2(adv_x), N)
    adv_y_hat = _truncate_spectrum(np.fft.fft2(adv_y), N)

    adv_x_hat = np.where(two_thirds, adv_x_hat, 0.0)
    adv_y_hat = np.where(two_thirds, adv_y_hat, 0.0)
    proj_x, proj_y = _project(adv_x_hat, adv_y_hat, kx, ky, k2)

    # The spatial integral of (u.grad)u is zero for periodic incompressible
    # flow.  Removing round-off in the zero Fourier mode prevents artificial
    # mean acceleration without changing the governing equations.
    proj_x[0, 0] = 0.0
    proj_y[0, 0] = 0.0
    return -proj_x, -proj_y


def _default_dt(t: float, nu: float, u0: float, k: float) -> float:
    if t == 0.0:
        return 1.0
    advective = 0.05 / max(u0 * k, 1e-30)
    viscous = 0.20 / max(nu * k * k, 1e-30)
    return min(t, advective, viscous)


def _advance(
    ux_hat: np.ndarray,
    uy_hat: np.ndarray,
    nu: float,
    kx: np.ndarray,
    ky: np.ndarray,
    k2: np.ndarray,
    two_thirds: np.ndarray,
    t_end: float,
    dt: float,
) -> ArrayPair:
    if t_end == 0.0:
        return ux_hat, uy_hat

    step = min(dt, t_end)
    if step <= 0.0:
        raise ValueError("dt must be positive")

    # Integrating-factor RK4: w = exp(nu |K|^2 t) u_hat, so viscosity is exact
    # and RK4 is applied only to the projected nonlinear term.
    def physical_hat(w_x: np.ndarray, w_y: np.ndarray, at: float) -> ArrayPair:
        damping = np.exp(-nu * k2 * at)
        return damping * w_x, damping * w_y

    def w_rhs(w_x: np.ndarray, w_y: np.ndarray, at: float) -> ArrayPair:
        px, py = physical_hat(w_x, w_y, at)
        rx, ry = _nonlinear_rhs(px, py, kx, ky, k2, two_thirds)
        growth = np.exp(nu * k2 * at)
        return growth * rx, growth * ry

    elapsed = 0.0
    w_x = ux_hat.copy()
    w_y = uy_hat.copy()
    while elapsed < t_end:
        h = min(step, t_end - elapsed)
        k1x, k1y = w_rhs(w_x, w_y, elapsed)
        k2x, k2y = w_rhs(w_x + 0.5 * h * k1x, w_y + 0.5 * h * k1y, elapsed + 0.5 * h)
        k3x, k3y = w_rhs(w_x + 0.5 * h * k2x, w_y + 0.5 * h * k2y, elapsed + 0.5 * h)
        k4x, k4y = w_rhs(w_x + h * k3x, w_y + h * k3y, elapsed + h)
        w_x = w_x + (h / 6.0) * (k1x + 2.0 * k2x + 2.0 * k3x + k4x)
        w_y = w_y + (h / 6.0) * (k1y + 2.0 * k2y + 2.0 * k3y + k4y)
        elapsed += h

    return physical_hat(w_x, w_y, t_end)


def _snapshot(t: float, ux: np.ndarray, uy: np.ndarray) -> Dict[str, Any]:
    speed = np.sqrt(ux * ux + uy * uy)
    return {
        "t": t,
        "ux": ux.copy(),
        "uy": uy.copy(),
        "kinetic_energy": float(0.5 * np.mean(ux * ux + uy * uy)),
        "max_speed": float(np.max(speed)),
    }


def tgv_reference(
    N: int,
    nu: float,
    u0: float,
    k: float,
    t: float,
    *,
    field_times: Optional[Iterable[float]] = None,
    dt: Optional[float] = None,
) -> Dict[str, Any]:
    """Solve 2D periodic TGV and return final metrics plus snapshots.

    Args:
        N: even grid size in each direction.
        nu: kinematic viscosity.
        u0: initial peak speed.
        k: magnitude of the TGV wavevector.  The square-axis components are
           k/sqrt(2), so E(t)/E(0) = exp(-2*nu*k*k*t).
        t: final time.
        field_times: optional snapshot times.  Defaults to [t].
        dt: optional RK step for the nonlinear integrating-factor solve.
    """
    _validate_inputs(N, nu, u0, k, t)
    times = sorted({float(v) for v in (field_times if field_times is not None else [t])})
    if not times or times[-1] > t or times[0] < 0.0:
        raise ValueError("field_times must be within [0, t]")

    kx, ky, k2, two_thirds = _wave_numbers(N, k)
    ux, uy = _initial_tgv(N, u0)
    ux_hat = np.fft.fft2(ux)
    uy_hat = np.fft.fft2(uy)
    step_dt = dt if dt is not None else _default_dt(t, nu, u0, k)

    snapshots: List[Dict[str, Any]] = []
    now = 0.0
    for target in times:
        ux_hat, uy_hat = _advance(
            ux_hat, uy_hat, nu, kx, ky, k2, two_thirds, target - now, step_dt
        )
        now = target
        ux_real = np.fft.ifft2(ux_hat).real
        uy_real = np.fft.ifft2(uy_hat).real
        snapshots.append(_snapshot(now, ux_real, uy_real))

    final = snapshots[-1]
    return {
        "t": t,
        "kinetic_energy": final["kinetic_energy"],
        "max_speed": final["max_speed"],
        "field_snapshots": snapshots,
    }


def _l2rel_pair(ux: np.ndarray, uy: np.ndarray, ref_ux: np.ndarray, ref_uy: np.ndarray) -> float:
    num = np.sum((ux - ref_ux) ** 2 + (uy - ref_uy) ** 2)
    den = np.sum(ref_ux ** 2 + ref_uy ** 2)
    return float(math.sqrt(num / den))


def _read_snapshot_json(path: Path) -> List[Dict[str, Any]]:
    raw = json.loads(path.read_text())
    rows = raw["snapshots"] if isinstance(raw, dict) and "snapshots" in raw else raw
    if not isinstance(rows, list) or not rows:
        raise ValueError("snapshot JSON must be a non-empty list or {'snapshots': [...]}")
    for i, row in enumerate(rows):
        if not all(key in row for key in ("t", "ux", "uy")):
            raise ValueError(f"snapshot {i} must contain t, ux, and uy")
    return rows


def _compare(args: argparse.Namespace) -> int:
    rows = _read_snapshot_json(Path(args.json_file))
    first_ux = np.asarray(rows[0]["ux"], dtype=np.float64)
    if first_ux.ndim != 2 or first_ux.shape[0] != first_ux.shape[1]:
        raise ValueError("ux/uy snapshots must be square 2D arrays")
    N = int(first_ux.shape[0])

    print("t,time,l2rel,max_lbm,max_ref,ke_ref")
    failed = 0
    for row in rows:
        ux = np.asarray(row["ux"], dtype=np.float64)
        uy = np.asarray(row["uy"], dtype=np.float64)
        if ux.shape != (N, N) or uy.shape != (N, N):
            raise ValueError("all ux/uy snapshots must share the first square shape")
        t_raw = float(row["t"])
        time = t_raw * args.time_scale
        ref = tgv_reference(N, args.nu, args.u0, args.k, time, dt=args.dt)
        snap = ref["field_snapshots"][-1]
        l2 = _l2rel_pair(ux, uy, snap["ux"], snap["uy"])
        max_lbm = float(np.max(np.sqrt(ux * ux + uy * uy)))
        print(
            f"{row['t']},{time:.12g},{l2:.12e},"
            f"{max_lbm:.12e},{snap['max_speed']:.12e},{snap['kinetic_energy']:.12e}"
        )
        if args.fail_above is not None and l2 > args.fail_above:
            failed += 1
    return 1 if failed else 0


def _self_verify() -> int:
    N = 32
    nu = 0.01
    u0 = 0.05
    k = math.sqrt(2.0)
    t = 1.0 / (nu * k * k)
    ref0 = tgv_reference(N, nu, u0, k, 0.0)
    ref1 = tgv_reference(N, nu, u0, k, t)
    observed = ref1["kinetic_energy"] / ref0["kinetic_energy"]
    expected = math.exp(-2.0)
    err = abs(observed - expected)
    assert err < 1e-3, f"energy decay {observed:.12e} != {expected:.12e}; err={err:.3e}"
    print(
        "PASS spectral_ns_referee self-verification: "
        f"E(t)/E(0)={observed:.12e}, expected={expected:.12e}, abs_err={err:.3e}"
    )
    return 0


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="2D Fourier pseudo-spectral Taylor-Green vortex referee"
    )
    sub = parser.add_subparsers(dest="cmd")

    self_p = sub.add_parser("self-test", help="run analytic TGV decay verification")
    self_p.set_defaults(func=lambda _args: _self_verify())

    cmp_p = sub.add_parser("compare", help="compare LBMFlow TGV JSON snapshots")
    cmp_p.add_argument("json_file")
    cmp_p.add_argument("--nu", type=float, required=True)
    cmp_p.add_argument("--u0", type=float, required=True)
    cmp_p.add_argument("--k", type=float, required=True,
                       help="TGV wavevector magnitude; use sqrt(2)*2*pi/N for mode-1 LBM grids")
    cmp_p.add_argument("--time-scale", type=float, default=1.0,
                       help="multiply JSON t by this factor before spectral comparison")
    cmp_p.add_argument("--dt", type=float, default=None,
                       help="optional nonlinear RK step")
    cmp_p.add_argument("--fail-above", type=float, default=None,
                       help="exit nonzero if any L2rel exceeds this threshold")
    cmp_p.set_defaults(func=_compare)
    return parser


def main(argv: Optional[List[str]] = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)
    if args.cmd is None:
        return _self_verify()
    try:
        return int(args.func(args))
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
