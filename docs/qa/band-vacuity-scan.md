# Band-vacuity scan — pending retighten queue (V&V master plan lane 2.2)

Original scan 2026-07-06 (commit 2e121c8). Rows below are the surviving
retighten candidates not yet landed. Delete a row when its retighten commits.

Retighten discipline: every retightened assert prints its measured value;
tightening is always allowed; a test failing a proposed tighter band is a
finding (reported, not "fixed"). No band may be LOOSENED by the retighten
order.

## Retighten queue

| Target | Current band | Raw-scan estimate | Action |
|---|---|---|---|
| `g2_analytic_strain.rs` (4 asserts: Couette + Poiseuille dissipation linf/mean) | 2e-4 / 1e-4 / 1e-6 / 1.1e-3 | measured ~1e-12..1e-15 | measure, retighten to measured×20, print measured values |
| `strain_rate.rs:91,120` (Couette Sxy, Poiseuille shear) | 2e-12 | ~1e-14 | same |
| `smoke_conservation.rs:108` + `validation_conservation.rs:86` (uniform-force momentum) | 1e-10 rel | ~1e-14 | same |
| `validation_channel.rs:126` (rotated Poiseuille) | 1e-10 | ~1e-13 | same |
| `accuracy_audit.rs:210` (acoustic damping gamma_rel) | 3e-1 | ~0.1-0.15 | measure; retighten with a stated physical error model (finite-k damping correction), not a bare ratio |
| `accuracy_audit.rs:273` (Galilean TRT defect) | 1.5e-1 | ~0.05-0.08 | same |
| `d3q19_smoke.rs:93` (TGV decay rate 15%) | 1.5e-1 | ~0.05 | same |
| `d3q19/d3q27` smoke lid gates (`>1e-4` existence-only) | lower-bound only | — | add upper bound + magnitude vs lid speed (two-layer rule); keep smoke label |
| `t15_3d.rs:686` (TGV3D decay rate 2%) | 2e-2 | ~0.5-1e-2 | measure; retighten if ≥5x slack confirmed |

Caveat: raw-scan estimates are estimates. The retighten order MUST measure
before tightening (band governance).

## Excluded (spec-frozen — do NOT retighten without band governance)

Ghia cavity RMS ≤ 0.02U (T7); duct flow ±0.5% (T15.2); sphere drag ±10/15%
(T15.3 triage); f32 tolerances (T6); cylinder Cd/Cl staircase bands (T8);
Stokes-I audit bands (first-measurement + stated headroom); effective-
viscosity ±2% (T1). `rotating_ibm.rs` loose profile gates (0.65/0.95/5.8)
are the only green witness of the weak-config path — retire as part of the
ANOM-P4-001 fix landing, not before.

## Systemic actions (fold into the same order)

1. Every retightened assert prints its measured value.
2. Two-layer rule retrofit for existence-only gates (upper bound).
3. No loosening; failure of a proposed tighter band = triage finding.
