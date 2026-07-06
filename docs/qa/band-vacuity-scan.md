# Band-vacuity scan — 2026-07-06 (V&V master plan lane 2.2)

Read-only sweep of numeric tolerances across `crates/lbm-core/tests`. Raw
scan by a subagent (48 asserts; 22 flagged VACUOUS-RISK), then PM-triaged
below. Caveat on the raw scan: several "measured" values were estimated,
not extracted — the retighten order must MEASURE before tightening (band
governance: tightening is always allowed; every retightened assert prints
its measured value).

## PM triage — surviving retighten candidates (queue: one codex order, W2)

| Target | Current band | Raw-scan estimate | Action |
|---|---|---|---|
| g2_analytic_strain.rs (4 asserts: Couette + Poiseuille dissipation linf/mean) | 2e-4 / 1e-4 / 1e-6 / 1.1e-3 | measured ~1e-12..1e-15 | measure, retighten to measured×20, add printed measured values |
| strain_rate.rs:91,120 (Couette Sxy, Poiseuille shear) | 2e-12 | ~1e-14 | same |
| smoke_conservation.rs:108 + validation_conservation.rs:86 (uniform-force momentum) | 1e-10 rel | ~1e-14 | same |
| validation_channel.rs:126 (rotated Poiseuille) | 1e-10 | ~1e-13 | same |
| accuracy_audit.rs:210 (acoustic damping gamma_rel) | 3e-1 | ~0.1-0.15 | measure; retighten with a stated physical error model (finite-k damping correction), not a bare ratio |
| accuracy_audit.rs:273 (Galilean TRT defect) | 1.5e-1 | ~0.05-0.08 | same |
| d3q19_smoke.rs:93 (TGV decay rate 15%) | 1.5e-1 | ~0.05 | same |
| d3q19/d3q27 smoke lid gates (`>1e-4` existence-only) | lower-bound only | — | add upper bound + magnitude expectation vs lid speed (two-layer rule); keep the smoke label |
| t15_3d.rs:686 (TGV3D decay rate 2%) | 2e-2 | ~0.5-1e-2 | measure; retighten if ≥5x slack confirmed |

## Superseded / no action

- rotating_ibm.rs loose profile gates (0.65/0.95/5.8): superseded by the
  adversarial cx/audit-ibm suite; these get retightened/retired as part of
  the ANOM-P4-001 fix landing, not before (they are the only green witness
  of the weak-config path until then).
- Spec-frozen bands correctly excluded from retightening without band
  governance: Ghia cavity RMS ≤ 0.02U (VALIDATION.md T7), duct flow ±0.5%
  (T15.2), sphere drag ±10/15% (T15.3 triage), f32 tolerances (T6),
  cylinder Cd/Cl staircase bands (T8), Stokes-I audit bands (frozen at
  first measurement with stated headroom), effective-viscosity ±2% (T1).
  The raw scan over-flagged these; they carry documented rationale.

## Systemic actions (fold into the same order)

1. Every retightened assert prints its measured value (assert-message
   convention already standard in audit files; retrofit these older files).
2. Two-layer rule retrofit for existence-only gates (upper bounds).
3. No band may be LOOSENED by this order; any test that fails a proposed
   tighter band is a triage finding, reported not "fixed".
