# SPEC: Bouzidi T8 drag-coefficient recovery — close the Cd band or prove the band wrong

Status: dispatch-ready (2026-07-06). Owner: PM. Implements in the existing `r2-bouzidi`
worktree (branch `r2-bouzidi`, commit `d3cd72d` present).

## 1. Context and evidence

The Bouzidi interpolated bounce-back (post-stream pass per SPEC_BOUZIDI_STL.md) is
implemented with the qd=1/2 bitwise BB degeneracy and Wen momentum-exchange force. The
first pass could not meet the T8 acceptance band and was frozen as characterization only:

- Measured (Re=20, cylinder D=20, 10,000 samples): **Cd = 5.83340474**, Cl = 0.00867670.
- Target band: **Cd ∈ [5.41, 5.75]** (T8, VALIDATION.md — the ±3% band around the
  half-way-BB-converged reference for the confined-cylinder configuration).
- Miss: +1.4% above the upper edge; +8% above band center.
- Deferred in pass 1: convergence slope D={10,20,40}, off-grid Poiseuille, GPU Bouzidi
  entry point.

A steady +5-8% Cd bias at D=20 is the classic signature of an **effective-diameter error
of O(1) lattice spacing** (Cd scales ≈ linearly with D at fixed Re here), so geometry
definitions are the first suspects, before the force formula.

## 2. Diagnosis matrix (execute in this order, record every number)

1. **Radius/qd convention audit.** The wall surface must be the *analytic circle*
   r = D/2 measured from the cylinder center in lattice units, with qd = fractional
   distance from the FLUID node to the surface along each link (0 < qd ≤ 1). Off-by-half
   traps: (a) treating cell centers as surface (adds ~0.5Δx to the radius), (b) computing
   qd from the solid node instead of the fluid node (qd ↦ 1−qd), (c) inconsistent center
   placement (on-node vs inter-node; the reference config places the center OFF-lattice
   at a half-integer coordinate — check what the test uses vs what the reference assumed).
   Cross-check: at qd = 1/2 exactly, the scheme must degenerate BITWISE to half-way BB
   (the existing degeneracy test pins this — but that test uses a flat wall; verify the
   *cylinder map generator* produces qd values consistent with the analytic circle by
   dumping min/max/mean qd and the count of links per boundary cell).
2. **Re definition.** nu from tau, U from the actual inflow profile (parabolic mean vs
   max — a 2/3 factor error in U moves Cd massively at Re=20). Confirm Re = U_mean·D/nu
   matches the reference's definition (Schäfer-Turek uses U_mean for 2D-1/2D-2).
3. **Blockage and domain.** The band was frozen for a specific H/D and inlet distance.
   Confirm the test grid matches; at Re=20 confined flow, H/D=4 vs 5 shifts Cd by several %.
4. **Force evaluation.** Wen (2014) Galilean-invariant momentum exchange:
   F = Σ_links [f_i⁺(c_i − u_w) − f_ī(c_ī − u_w)] over boundary links crossing the surface,
   with populations taken at the correct time (post-collision incoming, post-stream
   returning). Audit: u_w = 0 here so Galilean terms vanish — if pass 1 used plain ME
   (f c_i), that is NOT the source of an 8% error at rest. Check instead for double
   counting of links at cells with multiple cut links and for missing diagonal links.
5. **Convergence study (mandatory regardless of the above).** D = {10, 20, 40} at fixed
   Re=20 and fixed blockage (scale the domain with D). Expect ≈2nd-order convergence of
   Cd toward a limit. If the D→∞ extrapolation lands INSIDE [5.41, 5.75], the residual
   at D=20 is discretization and the acceptance moves to D=40 (band unchanged); if it
   lands OUTSIDE, there is a code bug — return to items 1–4.

## 3. Required changes

- Fix whatever the matrix convicts. Every fix must keep the qd=1/2 flat-wall bitwise
  degeneracy test green unmodified.
- Complete the pass-1 deferrals:
  - Off-grid Poiseuille: channel walls at non-integer positions via Bouzidi; parabolic
    profile error must beat half-way BB on the same (deliberately misaligned) grid.
  - Convergence test: the D={10,20,40} study above, committed as an `#[ignore]`d heavy
    test with the slope asserted (order ≥ 1.7 measured in the assert message).
- GPU Bouzidi remains OUT of scope (it lands with the ME-1 kernel work); make sure the
  CPU pass structure keeps the seam (post-stream pass over a link list) that
  SPEC_BOUZIDI_STL.md defines for the GPU port.

## 4. Acceptance

1. T8 band met: Cd ∈ [5.41, 5.75] at the acceptance resolution (D=20 if a bug is found
   and fixed; D=40 with the documented convergence slope if the matrix proves pure
   discretization). Cl bound per VALIDATION.md. Never loosen the band; if your evidence
   says the band itself encodes a different configuration (blockage/Re definition),
   STOP and report — do not re-freeze unilaterally (band governance).
2. `cargo test --workspace --release` green; degeneracy + off-grid Poiseuille +
   convergence tests in; measured values printed in asserts.
3. TESTING_NOTES.md append: the full diagnosis matrix results (all numbers), root cause,
   before/after Cd.

## 5. Rules

Work only in `r2-bouzidi`. English. Never commit red. No tolerance/band edits. Sandbox
git-commit failure → committed-ready note.
