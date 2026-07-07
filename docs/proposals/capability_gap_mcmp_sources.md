# Capability request — per-component (MCMP) volume sources

Filed 2026-07-06 by the D-track PM, addressed to the core-engine session.
Origin: V&V Axis-9.1 sparger experiment (ANOM-P4-020, main `5a2833d`,
harness on `cx/vv-sparger`).

## Finding

SCMP (single-component Shan-Chen) **cannot** express gas injection: phase
identity is local density, and a MassFlow source inside a liquid pool just
densifies the pool at the source (measured ρ_max = 2.385 super-liquid at
the injection cell across a 5× rate sweep, zero bubbles). Mechanism-level
finding, not a tuning problem.

## Ask

Add **per-component volume sources for MCMP** (`MultiComponent`, the two-
distribution Shan-Chen path). Injecting only the gas component's `f_g`
inside the liquid domain gives a real low-density-ratio sparger analog
before MF-γ (conservative Allen-Cahn phase-field) lands. Native `Solver`'s
CR-1 sources are single-distribution today.

## Downstream

REQ VR-STR-02 sparger unit test's home is now unambiguous: either MCMP
per-component sources (interim, low-ratio) or MF-γ phase-field gas inflow
(fidelity path). The V&V harness on `cx/vv-sparger` — detachment
detector, rise-velocity anchor, Laplace check, honesty clauses — transfers
directly to whichever lands first.

## D-track relevance

None on the deposition path (dispersed_seeding uses single-distribution
particle Lagrangian + CR-1 mass/jet sources for the LIQUID tray flow —
this capability gap is orthogonal to D-track's kill-closures work). Filed
here for cross-track visibility and to record the routing.
