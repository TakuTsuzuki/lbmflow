"""Anomaly-detection checks for the physics-QA sweep (stdlib only).

Every check returns a finding dict:
  {check, ok, expected, observed, detail}
run_sweep.py attaches severity/disposition when a check fails.

Parsers understand the CLI collection surface exactly as written by
crates/lbm-cli/src/runner.rs: field CSV (one `#` header line + ny rows of nx
values), legacy-ASCII VTK structured points (x fastest, then y, then z),
force.csv / point CSVs, manifest.json.
"""

import json
import math
import re
from pathlib import Path


# ------------------------------------------------------------------ parsers

def read_manifest(out_dir):
    return json.loads((Path(out_dir) / "manifest.json").read_text())


def read_grid_csv(path):
    """Field CSV -> (values row-major y*nx+x, nx, ny)."""
    lines = Path(path).read_text().splitlines()
    header = lines[0]
    m = re.search(r"nx=(\d+), ny=(\d+)", header)
    nx, ny = int(m.group(1)), int(m.group(2))
    vals = []
    for line in lines[1:]:
        if not line.strip():
            continue
        vals.extend(float(t) for t in line.split(","))
    assert len(vals) == nx * ny, f"{path}: {len(vals)} != {nx}*{ny}"
    return vals, nx, ny


def read_vtk(path):
    """Legacy ASCII structured-points VTK -> (values, nx, ny, nz).

    Layout: x fastest, then y, then z (cell = (z*ny + y)*nx + x).
    """
    lines = Path(path).read_text().splitlines()
    nx = ny = nz = None
    data_at = None
    for i, line in enumerate(lines):
        if line.startswith("DIMENSIONS"):
            _, a, b, c = line.split()
            nx, ny, nz = int(a), int(b), int(c)
        if line.startswith("LOOKUP_TABLE"):
            data_at = i + 1
            break
    vals = []
    for line in lines[data_at:]:
        if line.strip():
            vals.extend(float(t) for t in line.split())
    assert len(vals) == nx * ny * nz, f"{path}: {len(vals)} != {nx}*{ny}*{nz}"
    return vals, nx, ny, nz


def read_force_csv(path):
    """force.csv -> list of (step, fx, fy[, fz])."""
    rows = []
    for line in Path(path).read_text().splitlines()[1:]:
        if line.strip():
            parts = line.split(",")
            rows.append(tuple(float(p) for p in parts))
    return rows


def snapshots(out_dir, kind):
    """All `<kind>_<step>.csv|vtk` files sorted by step -> [(step, path)]."""
    hits = []
    for p in Path(out_dir).iterdir():
        m = re.fullmatch(rf"{kind}_(\d+)\.(csv|vtk)", p.name)
        if m:
            hits.append((int(m.group(1)), p))
    return sorted(hits)


def read_field_any(path):
    """CSV or VTK -> (values, nx, ny, nz)."""
    if str(path).endswith(".vtk"):
        return read_vtk(path)
    v, nx, ny = read_grid_csv(path)
    return v, nx, ny, 1


def _region_sum(vals, nx, ny, nz, region):
    if region == "all":
        return sum(vals)
    if region == "interior":  # 2D: strip the 1-cell rim
        return sum(vals[y * nx + x] for y in range(1, ny - 1) for x in range(1, nx - 1))
    if region == "interior3d":
        return sum(vals[(z * ny + y) * nx + x]
                   for z in range(1, nz - 1) for y in range(1, ny - 1)
                   for x in range(1, nx - 1))
    raise ValueError(region)


def _finding(check, ok, expected, observed, detail=""):
    return {"check": check, "ok": bool(ok), "expected": expected,
            "observed": observed, "detail": detail}


# ------------------------------------------------------------------- checks
# Each check(out_dir, cfg, args) -> finding (or list of findings).

def check_finite_and_status(out_dir, cfg, args=None):
    man = read_manifest(out_dir)
    bad = man["status"] == "diverged" or not math.isfinite(man["diagnostics"]["maxSpeed"]) \
        or not math.isfinite(man["diagnostics"]["totalMass"])
    return _finding(
        "finite/status", not bad,
        "status in {completed, steady}, finite diagnostics",
        f"status={man['status']}, maxSpeed={man['diagnostics']['maxSpeed']:.4g}, "
        f"totalMass={man['diagnostics']['totalMass']:.6g}")


def check_speed_ceiling(out_dir, cfg, args=None):
    man = read_manifest(out_dir)
    ms = man["diagnostics"]["maxSpeed"]
    return _finding("|u| hard ceiling", ms <= 0.3,
                    "max|u| <= 0.3 (MAX_SPEED, low-Mach hard limit)",
                    f"max|u| = {ms:.4g}")


def check_speed_scale(out_dir, cfg, args):
    man = read_manifest(out_dir)
    ms = man["diagnostics"]["maxSpeed"]
    lim = args["factor"] * args["u_ref"]
    return _finding("|u| scale sanity", ms <= lim,
                    f"max|u| <= {lim:.3g} ({args['factor']}x driving speed)",
                    f"max|u| = {ms:.4g}")


def check_mass_drift(out_dir, cfg, args):
    snaps = snapshots(out_dir, "rho")
    if len(snaps) < 2:
        return _finding("mass drift", False, ">=2 rho snapshots",
                        f"{len(snaps)} snapshot(s)", "collection gap")
    sums, steps = [], []
    for step, path in snaps:
        v, nx, ny, nz = read_field_any(path)
        sums.append(_region_sum(v, nx, ny, nz, args["region"]))
        steps.append(step)
    m0 = sums[0]
    drift = max(abs(s - m0) / abs(m0) for s in sums[1:])
    span = max(steps) - min(steps)
    band = args["band_per_1e4"] * max(1.0, span / 1e4)
    return _finding(
        "mass drift", drift <= band,
        f"rel drift <= {band:.3g} over {span} steps "
        f"({args['band_per_1e4']:.0e}/1e4 steps, region={args['region']})",
        f"rel drift = {drift:.3g}",
        f"mass[first]={m0:.12g}, mass[last]={sums[-1]:.12g}")


def check_momentum_growth(out_dir, cfg, args):
    """Uniform force on a periodic box: du/dt = F/rho exactly (T6)."""
    snaps = snapshots(out_dir, "ux")
    vals = []
    for step, path in snaps:
        v, nx, ny, nz = read_field_any(path)
        vals.append((step, sum(v) / len(v)))
    worst = 0.0
    for (s0, u0), (s1, u1) in zip(vals, vals[1:]):
        slope = (u1 - u0) / (s1 - s0)
        worst = max(worst, abs(slope - args["fx"]) / args["fx"])
    return _finding("momentum growth", worst <= args["band"],
                    f"du/dstep = F = {args['fx']:.3g} (rel err <= {args['band']:.0e})",
                    f"max rel err = {worst:.3g}",
                    f"{len(vals)} snapshots")


def check_field_uniform(out_dir, cfg, args):
    snaps = snapshots(out_dir, args["field"])
    worst = 0.0
    for step, path in snaps:
        v, *_ = read_field_any(path)
        worst = max(worst, max(v) - min(v))
    return _finding("field uniformity", worst <= args["band"],
                    f"max-min <= {args['band']:.0e} ({args['field']}, all snapshots)",
                    f"max-min = {worst:.3g}")


def _final_field(out_dir, kind):
    snaps = snapshots(out_dir, kind)
    assert snaps, f"no {kind} output in {out_dir}"
    return read_field_any(snaps[-1][1])


def check_poiseuille_exact(out_dir, cfg, args):
    v, nx, ny, _ = _final_field(out_dir, "ux")
    h = ny - 2
    g, nu = args["g"], args["nu"]
    ref = [g / (2 * nu) * (j - 0.5) * (h - (j - 0.5)) for j in range(1, ny - 1)]
    x = nx // 2
    prof = [v[j * nx + x] for j in range(1, ny - 1)]
    umax = max(abs(r) for r in ref)
    linf = max(abs(p - r) for p, r in zip(prof, ref)) / umax
    return _finding("Poiseuille exactness", linf <= args["band"],
                    f"LInf_rel <= {args['band']:.0e} vs u(y)=g/(2nu) y_w(H-y_w)",
                    f"LInf_rel = {linf:.3g}",
                    f"u_max sim={max(prof):.6g} ref={max(ref):.6g}")


def check_poiseuille_l2(out_dir, cfg, args):
    """L2rel vs analytic (order input for the BGK pair; always 'ok')."""
    v, nx, ny, _ = _final_field(out_dir, "ux")
    h = ny - 2
    g, nu = args["g"], args["nu"]
    ref = [g / (2 * nu) * (j - 0.5) * (h - (j - 0.5)) for j in range(1, ny - 1)]
    x = nx // 2
    prof = [v[j * nx + x] for j in range(1, ny - 1)]
    num = math.sqrt(sum((p - r) ** 2 for p, r in zip(prof, ref)))
    den = math.sqrt(sum(r ** 2 for r in ref))
    l2 = num / den
    f = _finding("Poiseuille L2rel (order input)", True, "recorded", f"L2rel = {l2:.4g}")
    f["l2rel"] = l2
    return f


def check_profile_symmetry(out_dir, cfg, args):
    v, nx, ny, _ = _final_field(out_dir, "ux")
    x = nx // 2
    h = ny - 2
    worst = max(abs(v[j * nx + x] - v[(ny - 1 - j) * nx + x])
                for j in range(1, ny // 2))
    umax = max(abs(v[j * nx + x]) for j in range(1, ny - 1))
    return _finding("top/bottom symmetry", worst <= args["band"],
                    f"|ux(j) - ux(H+1-j)| <= {args['band']:.0e} (T2 absolute)",
                    f"max asym = {worst:.3g} (umax={umax:.3g})")


def check_couette_exact(out_dir, cfg, args):
    v, nx, ny, _ = _final_field(out_dir, "ux")
    h = ny - 2
    u = args["u_wall"]
    ref = [u * (j - 0.5) / h for j in range(1, ny - 1)]
    x = nx // 2
    prof = [v[j * nx + x] for j in range(1, ny - 1)]
    linf = max(abs(p - r) for p, r in zip(prof, ref)) / u
    return _finding("Couette exactness", linf <= args["band"],
                    f"LInf_rel <= {args['band']:.0e} vs u(y)=U y_w/H",
                    f"LInf_rel = {linf:.3g}")


def check_t4_flow_rate(out_dir, cfg, args):
    ux, nx, ny, _ = _final_field(out_dir, "ux")
    rho, *_ = _final_field(out_dir, "rho")
    x_hi = nx - 1 - args["outlet_margin"]
    qs = []
    for x in range(1, x_hi + 1):
        qs.append(sum(rho[j * nx + x] * ux[j * nx + x] for j in range(1, ny - 1)))
    qbar = sum(qs) / len(qs)
    dev = max(abs(q - qbar) for q in qs) / abs(qbar)
    return _finding("T4 bulk flow-rate constancy", dev <= args["band"],
                    f"max|Q-Qbar|/Qbar <= {args['band']:.0e} "
                    f"(bulk x in [1,{x_hi}]; measured 2.4e-5 in spec)",
                    f"max dev = {dev:.3g}", f"Qbar = {qbar:.6g}")


def check_t4_profile(out_dir, cfg, args):
    ux, nx, ny, _ = _final_field(out_dir, "ux")
    h = ny - 2
    um = args["umax"]
    ref = [4 * um * (j - 0.5) * (h - (j - 0.5)) / (h * h) for j in range(1, ny - 1)]
    x = nx // 2
    prof = [ux[j * nx + x] for j in range(1, ny - 1)]
    num = math.sqrt(sum((p - r) ** 2 for p, r in zip(prof, ref)))
    den = math.sqrt(sum(r ** 2 for r in ref))
    return _finding("T4 central profile", num / den <= args["band"],
                    f"L2rel <= {args['band']:.0e} vs parabola",
                    f"L2rel = {num / den:.3g}")


def check_ghia_rms(out_dir, cfg, args):
    import matrix as mx
    u_ref, v_ref, tol_u, exclude_typo = mx.GHIA[args["re"]]
    ux, nx, ny, _ = _final_field(out_dir, "ux")
    uy, *_ = _final_field(out_dir, "uy")
    n = nx
    big_l = float(n - 2)
    u_lid = 0.1

    def sample(vals, frac, along_y):
        pos = 0.5 + frac * big_l
        c0 = int(min(max(math.floor(pos), 1), n - 2))
        c1 = min(c0 + 1, n - 2)
        t = pos - c0
        mid = n // 2
        if along_y:
            return (1 - t) * vals[c0 * nx + mid] + t * vals[c1 * nx + mid]
        return (1 - t) * vals[mid * nx + c0] + t * vals[mid * nx + c1]

    total, cnt = 0.0, 0
    for i in range(17):
        du = sample(ux, mx.GHIA_Y[i], True) - u_lid * u_ref[i]
        total += du * du
        cnt += 1
        if exclude_typo and abs(mx.GHIA_X[i] - 0.9063) < 1e-12:
            continue
        dv = sample(uy, mx.GHIA_X[i], False) - u_lid * v_ref[i]
        total += dv * dv
        cnt += 1
    rms = math.sqrt(total / cnt)
    return _finding(f"Ghia RMS (Re={args['re']})", rms <= tol_u * u_lid,
                    f"RMS <= {tol_u} * U = {tol_u * u_lid:.4g}",
                    f"RMS = {rms:.4g} ({rms / u_lid:.4f} U)")


def check_cd_cl_steady(out_dir, cfg, args):
    rows = read_force_csv(Path(out_dir) / "force.csv")
    scale = 2.0 / (args["u_mean"] ** 2 * args["d"])
    sel = [(r[1] * scale, r[2] * scale) for r in rows if r[0] >= args["sample_start"]]
    cd = sum(a for a, _ in sel) / len(sel)
    cl = sum(b for _, b in sel) / len(sel)
    lo, hi = args["cd_band"]
    clo, chi = args["cl_band"]
    out = [
        _finding("T8 Cd band", lo <= cd <= hi,
                 f"Cd in [{lo}, {hi}] (ref 5.5795)", f"Cd = {cd:.4f}",
                 f"{len(sel)} samples from step {args['sample_start']}"),
        _finding("T8 Cl band", clo <= cl <= chi,
                 f"Cl in [{clo}, {chi}] (ref 0.0106)", f"Cl = {cl:.4f}"),
    ]
    return out


def check_karman(out_dir, cfg, args):
    rows = read_force_csv(Path(out_dir) / "force.csv")
    scale = 2.0 / (args["u_mean"] ** 2 * args["d"])
    win = [(r[0], r[1] * scale, r[2] * scale) for r in rows if r[0] >= args["window_start"]]
    steps = [w[0] for w in win]
    cds = [w[1] for w in win]
    cls = [w[2] for w in win]
    half = len(cls) // 2
    mean_cl = sum(cls) / len(cls)
    amp1 = max(abs(c - mean_cl) for c in cls[:half])
    amp2 = max(abs(c - mean_cl) for c in cls[half:])
    saturated = amp2 > 0 and abs(amp2 - amp1) / amp2 <= 0.15
    # upward zero crossings of Cl - mean
    xs = []
    for (s0, c0), (s1, c1) in zip(zip(steps, cls), zip(steps[1:], cls[1:])):
        a0, a1 = c0 - mean_cl, c1 - mean_cl
        if a0 < 0.0 <= a1:
            xs.append(s0 + (s1 - s0) * (-a0) / (a1 - a0))
    periods = [b - a for a, b in zip(xs, xs[1:])]
    out = [_finding("Karman saturation guard", saturated,
                    "Cl amplitude change <= 15% between window halves",
                    f"amp1={amp1:.3f}, amp2={amp2:.3f}",
                    "if false: shedding not yet saturated -> St/Cd inconclusive")]
    if not periods:
        out.append(_finding("Strouhal", False, "detectable Cl oscillation",
                            "no zero crossings in window",
                            "no shedding detected"))
        return out
    pm = sum(periods) / len(periods)
    st = args["d"] / (args["u_mean"] * pm)
    pv = max(abs(p - pm) for p in periods) / pm
    lo, hi = args["st_band"]
    out.append(_finding("Strouhal", lo <= st <= hi,
                        f"St in [{lo}, {hi}] (ref 0.295-0.305)",
                        f"St = {st:.4f} (period {pm:.1f} steps, n={len(periods)})"))
    out.append(_finding("period regularity", pv <= args["period_var"],
                        f"consecutive-period variation <= {args['period_var']:.0%}",
                        f"max variation = {pv:.2%}"))
    clo, chi = args["cdmax_band"]
    out.append(_finding("Cd_max band", clo <= max(cds) <= chi,
                        f"Cd_max in [{clo}, {chi}]", f"Cd_max = {max(cds):.4f}"))
    llo, lhi = args["clmax_band"]
    clmax = max(abs(c) for c in cls)
    out.append(_finding("Cl_max band", llo <= clmax <= lhi,
                        f"|Cl|_max in [{llo}, {lhi}]", f"|Cl|_max = {clmax:.4f}"))
    return out


def check_reverse_flow(out_dir, cfg, args):
    out = []
    for step, path in snapshots(out_dir, "ux"):
        ux, nx, ny, _ = read_field_any(path)
        rho_path = Path(out_dir) / path.name.replace("ux", "rho")
        if not rho_path.exists():
            continue
        rho, *_ = read_field_any(rho_path)
        x_out, x_in = nx - 1, 0
        rev = -sum(min(ux[j * nx + x_out], 0.0) * rho[j * nx + x_out]
                   for j in range(1, ny - 1))
        inflow = sum(ux[j * nx + x_in] * rho[j * nx + x_in] for j in range(1, ny - 1))
        frac = rev / inflow if inflow > 0 else float("inf")
        out.append(_finding(f"T9 reverse flow @step {step}", frac <= args["band"],
                            f"reverse mass flux <= {args['band']:.0%} of inflow",
                            f"{frac:.2%}"))
    return out or _finding("T9 reverse flow", False, "ux+rho snapshots", "none found")


def check_duct_ref(y_w, z_c, h, g, nu, nmax=99):
    """Rectangular-duct series, VALIDATION T15.2 notation: 2a = 2b = H."""
    a = h / 2.0
    b = h / 2.0
    s = 0.0
    for n in range(1, nmax + 1, 2):
        s += (1.0 / n ** 3) * (1 - math.cosh(n * math.pi * z_c / (2 * a))
                               / math.cosh(n * math.pi * b / (2 * a))) \
            * math.sin(n * math.pi * y_w / (2 * a))
    return 16 * a * a * g / (nu * math.pi ** 3) * s


def check_duct_exact(out_dir, cfg, args):
    v, nx, ny, nz, = _final_field(out_dir, "ux")
    h = ny - 2
    g, nu = args["g"], args["nu"]
    x = nx // 2
    worst, umax = 0.0, 0.0
    for k in range(1, nz - 1):
        for j in range(1, ny - 1):
            y_w = j - 0.5
            z_c = (k - 0.5) - h / 2.0
            ref = check_duct_ref(y_w, z_c, h, g, nu)
            sim = v[(k * ny + j) * nx + x]
            worst = max(worst, abs(sim - ref))
            umax = max(umax, abs(ref))
    return _finding("duct series exactness", worst / umax <= args["band"],
                    f"LInf_rel <= {args['band']:.0e} vs Fourier series (n<=99)",
                    f"LInf_rel = {worst / umax:.4g}")


def check_duct_flow_rate(out_dir, cfg, args):
    v, nx, ny, nz = _final_field(out_dir, "ux")
    h = ny - 2
    g, nu = args["g"], args["nu"]
    a = b = h / 2.0
    q_ref = 0.0
    for n in range(1, 100, 2):
        q_ref += (1.0 / n ** 4) * (2 * b - (4 * a / (n * math.pi))
                                   * math.tanh(n * math.pi * b / (2 * a)))
    q_ref *= 64 * a ** 3 * g / (nu * math.pi ** 4)
    x = nx // 2
    q = sum(v[(k * ny + j) * nx + x]
            for k in range(1, nz - 1) for j in range(1, ny - 1))
    dev = abs(q - q_ref) / q_ref
    return _finding("duct flow rate", dev <= args["band"],
                    f"|Q-Q_ref|/Q_ref <= {args['band']:.1%} (measured 0.094%)",
                    f"dev = {dev:.3%}", f"Q={q:.6g} Q_ref={q_ref:.6g}")


def check_duct_yz_symmetry(out_dir, cfg, args):
    v, nx, ny, nz = _final_field(out_dir, "ux")
    assert ny == nz
    x = nx // 2
    worst = max(abs(v[(k * ny + j) * nx + x] - v[(j * ny + k) * nx + x])
                for k in range(1, nz - 1) for j in range(1, ny - 1))
    umax = max(abs(v[(k * ny + j) * nx + x])
               for k in range(1, nz - 1) for j in range(1, ny - 1))
    return _finding("duct y<->z symmetry", worst <= args["band"] * umax,
                    f"|u(y,z)-u(z,y)| <= {args['band']:.0e} * umax",
                    f"max asym = {worst:.3g} (umax={umax:.3g})")


def check_z_mirror_symmetry(out_dir, cfg, args):
    v, nx, ny, nz = _final_field(out_dir, args["field"])
    worst = 0.0
    for z in range(nz // 2):
        zm = nz - 1 - z
        for y in range(1, ny - 1):
            for x in range(1, nx - 1):
                worst = max(worst, abs(v[(z * ny + y) * nx + x]
                                       - v[(zm * ny + y) * nx + x]))
    rel = worst / args["u_ref"]
    return _finding("3D cavity z-mirror symmetry", rel <= args["band"],
                    f"max|f(z)-f(mirror)|/U <= {args['band']:.0e} "
                    "(T15.5 sentinel measured ~2e-15)",
                    f"rel asym = {rel:.3g}")


def check_xy_mirror_symmetry(out_dir, cfg, args):
    """Mirror about the droplet center with periodic wrap.

    The droplet sits at integer (cx, cy), so the IC's symmetry map on the
    periodic grid is x -> (2*cx - x) mod nx (NOT nx-1-x, which is the mirror
    about the half-integer domain center and does not fix the IC).
    """
    v, nx, ny, _ = _final_field(out_dir, args["field"])
    cx = int(cfg["scenario"]["init"]["cx"])
    cy = int(cfg["scenario"]["init"]["cy"])
    worst = 0.0
    for y in range(ny):
        for x in range(nx):
            xm = (2 * cx - x) % nx
            ym = (2 * cy - y) % ny
            worst = max(worst, abs(v[y * nx + x] - v[y * nx + xm]),
                        abs(v[y * nx + x] - v[ym * nx + x]))
    return _finding("droplet x/y mirror symmetry", worst <= args["band"],
                    f"max mirror asym <= {args['band']:.0e} (observation band)",
                    f"max asym = {worst:.3g}")


def _sc_pressure(rho, g):
    psi = 1.0 - math.exp(-rho)
    return rho / 3.0 + (g / 6.0) * psi * psi


def check_spurious_current(out_dir, cfg, args):
    man = read_manifest(out_dir)
    ms = man["diagnostics"]["maxSpeed"]
    return _finding("spurious currents", ms <= args["band"],
                    f"max|u| <= {args['band']:.0e} at equilibrium "
                    "(T11; flat-interface measured 1.26e-3)",
                    f"max|u| = {ms:.4g}")


def check_laplace_sigma(out_dir, cfg, args):
    rho, nx, ny, _ = _final_field(out_dir, "rho")
    g = args["g"]
    lo, hi = min(rho), max(rho)
    thr = 0.5 * (lo + hi)
    area = sum(1 for r in rho if r > thr)
    r_fit = math.sqrt(area / math.pi)
    cx, cy = (nx - 1) / 2.0, (ny - 1) / 2.0
    p_in, n_in, p_out, n_out = 0.0, 0, 0.0, 0
    for y in range(ny):
        for x in range(nx):
            d = math.hypot(x - cx, y - cy)
            if d <= 0.5 * r_fit:
                p_in += _sc_pressure(rho[y * nx + x], g)
                n_in += 1
            elif d >= r_fit + 10:
                p_out += _sc_pressure(rho[y * nx + x], g)
                n_out += 1
    sigma = (p_in / n_in - p_out / n_out) * r_fit
    ref = args["sigma_ref"]
    dev = abs(sigma - ref) / ref
    return _finding("Laplace sigma", dev <= args["band_rel"],
                    f"sigma = {ref:.3g} +- {args['band_rel']:.0%} (T11 slope band)",
                    f"sigma = {sigma:.4g} (dev {dev:.1%})",
                    f"R_fit = {r_fit:.2f}, dp = {sigma / r_fit:.4g}")


def check_contact_angle(out_dir, cfg, args):
    rho, nx, ny, _ = _final_field(out_dir, "rho")
    lo, hi = min(rho), max(rho)
    thr = 0.5 * (lo + hi)
    blob = [(x, y) for y in range(1, ny - 1) for x in range(nx)
            if rho[y * nx + x] > thr]
    if not blob:
        return _finding("contact angle", False, "a liquid blob", "none found")
    y_max = max(y for _, y in blob)
    h = y_max - 0.5  # cell-center height above the half-way wall surface (y=0.5)
    row1 = [x for x, y in blob if y == 1]
    w = max(row1) - min(row1) + 1 if row1 else 0
    if w == 0:
        return _finding("contact angle", False,
                        "droplet touching the wall (contact row occupied)",
                        "no contact-line cells",
                        "droplet may have detached (non-wetting)")
    theta = math.degrees(2 * math.atan2(2 * h, w))
    ref, tol = args["theta_ref"], args["tol"]
    return _finding("contact angle (spherical cap)", abs(theta - ref) <= tol,
                    f"theta = {ref} +- {tol} deg (T11c, wallRho=1.0)",
                    f"theta = {theta:.1f} deg", f"h={h}, w={w}")


CHECKS = {
    "mass_drift": check_mass_drift,
    "momentum_growth": check_momentum_growth,
    "field_uniform": check_field_uniform,
    "poiseuille_exact": check_poiseuille_exact,
    "poiseuille_l2": check_poiseuille_l2,
    "profile_symmetry": check_profile_symmetry,
    "couette_exact": check_couette_exact,
    "t4_flow_rate": check_t4_flow_rate,
    "t4_profile": check_t4_profile,
    "ghia_rms": check_ghia_rms,
    "cd_cl_steady": check_cd_cl_steady,
    "karman": check_karman,
    "reverse_flow": check_reverse_flow,
    "duct_exact": check_duct_exact,
    "duct_flow_rate": check_duct_flow_rate,
    "duct_yz_symmetry": check_duct_yz_symmetry,
    "z_mirror_symmetry": check_z_mirror_symmetry,
    "xy_mirror_symmetry": check_xy_mirror_symmetry,
    "spurious_current": check_spurious_current,
    "laplace_sigma": check_laplace_sigma,
    "contact_angle": check_contact_angle,
    "speed_scale": check_speed_scale,
}
