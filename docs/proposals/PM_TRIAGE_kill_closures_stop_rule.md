# PM triage — kill-closures STOP-RULE outcome

Date: 2026-07-06. Author: D-track PM (Fable).

## What happened

Order `cx/kill-deposition-closures` completed successfully by triggering
its **STOP-RULE** (per the physics-discipline skill). The dispersed_seeding
example was rebuilt on resolved physics only — no harshness switch, no
analytic jet/wall-jet superposition, no calibrated dispersion constant, no
side-wall clamps, no direct agitation kicks, no reservoir scoring
heuristics. Particles co-evolve with the tray LBM and trilinearly sample
the live velocity field; agitation is a non-inertial-frame pseudo-force on
both fluid and particles; reservoir extraction is the analytic 1D
settling-column backtrace.

Result: **gentle sample deposits 0/10000 particles** (all suspended);
**density map is uniformly black** — the honest visual of "no ad-hoc
transport, no deposition". Metrics: `n_deposited=0, n_suspended=10000,
Ma=0.093, tau=0.536`. Wall time 1607 s.

## Why (mechanism, evidence)

The exported near-floor radial velocity profile
(`out/dispersed_seeding/gentle/near_floor_radial_velocity.csv`) is the
smoking gun:

- Mean |u_r| in the near-floor band across radii ~6–119 cells:
  **~1e-6 to 3e-6 m/s**.
- 100 µm bead Stokes settling: **v_s = 2.72e-4 m/s** — **2 orders of
  magnitude larger** than the near-floor flow.
- Time for a bead released at nozzle height (~8 mm) to sediment purely
  ballistically: **~30 s**.
- Protocol time (12000 steps × dt = 4.22e-4 s): **~5.1 s**.

The old CV=1.163 gentle result was produced by two closures conspiring:
(a) the analytic wall-jet term supplied radial transport 2 orders larger
than the resolved field actually carries at this tray resolution, and
(b) the side-wall clamp accumulated particles at the rim once transported
there. Neither was physics.

## The four honest options

The STOP-RULE report listed three; I refine and rank them for PM decision:

**A. Extend the protocol time budget** to physically-required duration.
- Cost: ~30 s protocol → step count grows ×6 → wall time ~2.5 h per sample
  (from measured 1607 s at 5 s). CpuSimd may cut that in half.
- Physical validity: NONE ADDED. It simulates a longer static-tray settle
  after ejection, which is honest and reveals what the resolved flow
  actually does. If the answer is "the ring is a resolved wall-jet
  transient", we will see it appear.
- Recommended action: run gentle at ×6 protocol duration on CpuSimd as a
  diagnostic BEFORE deciding B or C. Cheap on the ordering side (one
  parameter tweak); the wall time is the price of finding out.

**B. Adopt the adhesion-capture closure stack** already researched
(`docs/proposals/adhesion-capture-closure.md`).
- C1 (contact capture within one radius with named α, default 1.0) +
  C4 (JKR rolling-detachment threshold for silicone/PDMS) — both come
  with derivation, validity domain, validation test, and draft
  PHYSICS.md entries per Rule 1.
- This is the smallest defensible model addition to make deposition
  physically real without reintroducing the banned patterns.
- Requires **one user input** (already flagged): sign and rough magnitude
  of Δρ. If Δρ < 0 (creaming), floor deposition is the wrong framing
  entirely for the real 20 µm material — do NOT build B before this
  question is answered.

**C. Change the protocol / geometry** so the resolved field actually
carries particles to the floor in-budget: shorter drop height (thinner
tray), stronger jet, agitation-driven convection cell. This is a spec
change, honest, and probably necessary in combination with B for the real
20 µm near-neutral case.

**D. Rejected: add a "turbulent dispersion" fudge to compensate.** This is
the same shape of closure the discipline skill bans; it would resurrect
the edge-ring artifact under a different name.

## Recommendation

1. **Immediate diagnostic (A-only, no code changes)**: rerun gentle with
   protocol duration ×6 to observe what the resolved field actually does
   over the physical timescale. Wall time cost is acceptable for a
   one-shot answer. Deliverable: 2D density map + comparison to the
   black-map above.
2. **Simultaneously, block Phase B on the user's Δρ answer** (sign +
   rough magnitude, and shape: solid / capsule / flake). Adhesion-capture
   closure adoption is ready-to-order but conditional on that answer.
3. **After A's result**: if the resolved flow shows a physical mechanism
   for a deposition pattern (even if quantitative details are off), then
   B adds the wall interaction and Phase B reparameterizes to 20 µm. If
   A shows the resolved flow is essentially quiescent for the settling
   time, the honest conclusion is that this example needs a different
   protocol/geometry (C) — and that is a spec change to propose to the
   user, not a knob to add.

## What the discipline achieved (record)

The STOP-RULE fired exactly as designed. A codex order did NOT
recalibrate constants to force a plausible-looking density map when the
resolved physics could not produce one. A pre-discipline order would
have hit exactly this gap and quietly restored some fraction of the
closure layer to "make it work"; the black density map is the visible
proof of the difference. The V&V-side inventory + Order E cumulant
verdict + kill-closures STOP-RULE landing on the same day validates the
whole design-review loop we adopted at 2026-07-06.
