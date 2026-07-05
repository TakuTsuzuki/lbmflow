"""Agreement metrics for accuracy-audit analysis (stdlib only).

1:1 mirror (names + semantics) of the Rust library
`crates/lbm-core/tests/common/metrics.rs`. Rust is the source of truth for
in-suite assertions; this module serves post-hoc analysis of CLI/scenario
output (compose with the parsers in qa_checks.py). Keep the two in sync —
if you change one, change the other in the same commit.

Run self-tests: python3 scripts/qa/metrics.py
"""

import math
from typing import Callable, List, NamedTuple, Sequence, Tuple


def l2_rel(actual: Sequence[float], reference: Sequence[float]) -> float:
    """Relative L2 error ||a - r||_2 / ||r||_2."""
    assert len(actual) == len(reference)
    num = sum((a - r) ** 2 for a, r in zip(actual, reference))
    den = sum(r * r for r in reference)
    return math.sqrt(num / den)


def linf_rel(actual: Sequence[float], reference: Sequence[float],
             floor: float) -> float:
    """Relative L-infinity error: max|a-r| / max(max|r|, floor)."""
    assert len(actual) == len(reference)
    dmax = max(abs(a - r) for a, r in zip(actual, reference))
    rmax = max(max(abs(r) for r in reference), floor)
    return dmax / rmax


class LinFit(NamedTuple):
    """y ~ slope*x + intercept; r2 rejects sloppy fits (assert it too)."""
    slope: float
    intercept: float
    r2: float


def linear_fit(x: Sequence[float], y: Sequence[float]) -> LinFit:
    assert len(x) == len(y) and len(x) >= 2
    n = float(len(x))
    mx, my = sum(x) / n, sum(y) / n
    sxx = sum((xi - mx) ** 2 for xi in x)
    sxy = sum((xi - mx) * (yi - my) for xi, yi in zip(x, y))
    syy = sum((yi - my) ** 2 for yi in y)
    slope = sxy / sxx
    intercept = my - slope * mx
    r2 = 1.0 if syy == 0.0 else (sxy * sxy) / (sxx * syy)
    return LinFit(slope, intercept, r2)


def order_fit(h: Sequence[float], err: Sequence[float]) -> LinFit:
    """Observed convergence order (log-log slope). Assert slope AND r2."""
    return linear_fit([math.log(v) for v in h], [math.log(v) for v in err])


def envelope_fit(y: Sequence[float], amp: Sequence[float]) -> LinFit:
    """amp ~ A*exp(-k*y): fit of ln(amp) on y; A = exp(intercept), k = -slope."""
    return linear_fit(list(y), [math.log(v) for v in amp])


def phase_fit(t: Sequence[float], signal: Sequence[float],
              omega: float) -> Tuple[float, float]:
    """(amplitude, phase) with signal ~ amplitude*sin(omega*t + phase).

    Quadrature projection; sample an integer number of periods.
    """
    assert len(t) == len(signal) and t
    n = float(len(t))
    s = sum(si * math.sin(omega * ti) for ti, si in zip(t, signal))
    c = sum(si * math.cos(omega * ti) for ti, si in zip(t, signal))
    a_sin, a_cos = 2.0 * s / n, 2.0 * c / n
    return math.hypot(a_sin, a_cos), math.atan2(a_cos, a_sin)


def monotonicity(xs: Sequence[float]) -> float:
    """Fraction of adjacent pairs strictly decreasing (1.0 = monotone decay)."""
    assert len(xs) >= 2
    dec = sum(1 for a, b in zip(xs, xs[1:]) if b < a)
    return dec / (len(xs) - 1)


class CurveAgreement(NamedTuple):
    max_rel_dev: float
    worst_x: float
    frac_in_band: float


def curve_agreement(theory: Callable[[float], float],
                    samples: Sequence[Tuple[float, float]],
                    rel_band: float, floor: float) -> CurveAgreement:
    """A3 primitive: each sample must lie ON the theory curve, point by point."""
    assert samples
    max_rel_dev, worst_x, in_band = 0.0, samples[0][0], 0
    for x, measured in samples:
        th = theory(x)
        dev = abs(measured - th) / max(abs(th), floor)
        if dev > max_rel_dev:
            max_rel_dev, worst_x = dev, x
        if dev <= rel_band:
            in_band += 1
    return CurveAgreement(max_rel_dev, worst_x, in_band / len(samples))


def _selftest() -> None:
    # Mirrors crates/lbm-core/tests/metrics_selftest.rs fixtures.
    assert abs(l2_rel([1, 2.3, 3, 4], [1, 2, 3, 4]) - 0.3 / math.sqrt(30)) < 1e-15
    assert abs(linf_rel([1, 2.3, 3, 4], [1, 2, 3, 4], 0.0) - 0.3 / 4) < 1e-15
    assert abs(linf_rel([0.1], [0.0], 1.0) - 0.1) < 1e-15

    h: List[float] = [0.1, 0.05, 0.025, 0.0125]
    fit = order_fit(h, [3 * x * x for x in h])
    assert abs(fit.slope - 2.0) < 1e-12 and fit.r2 > 1 - 1e-12
    assert order_fit(h, [3e-2, 8e-3, 5e-3, 4.8e-3]).r2 < 0.98

    y = [0.5, 1.0, 2.0, 3.5, 5.0]
    fit = envelope_fit(y, [0.7 * math.exp(-0.35 * v) for v in y])
    assert abs(fit.slope + 0.35) < 1e-12
    assert abs(math.exp(fit.intercept) - 0.7) < 1e-12

    omega = 2 * math.pi / 100
    t = list(range(400))
    amp, phase = phase_fit(t, [0.02 * math.sin(omega * ti + 0.6) for ti in t], omega)
    assert abs(amp - 0.02) < 1e-6 and abs(phase - 0.6) < 1e-6

    assert monotonicity([4, 3, 2, 1]) == 1.0
    assert abs(monotonicity([4, 3, 3.5, 1]) - 2 / 3) < 1e-15

    agree = curve_agreement(lambda x: x * x,
                            [(1, 1), (2, 4), (3, 9.9), (4, 16)], 0.05, 0.0)
    assert abs(agree.max_rel_dev - 0.1) < 1e-12
    assert agree.worst_x == 3 and abs(agree.frac_in_band - 0.75) < 1e-12
    print("metrics.py selftest OK")


if __name__ == "__main__":
    _selftest()
