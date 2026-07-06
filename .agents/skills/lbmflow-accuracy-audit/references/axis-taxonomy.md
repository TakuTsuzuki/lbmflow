# The five attack axes (P1 excavation scaffold)

Each axis is a *way an approximation's error becomes first-order visible*.
For every approximation in the target subsystem, walk all five axes and ask
the axis question; most approximations are killable on 2-3 of them. Every
axis below carries a worked example that found a real defect (or ruled a
design) in this repo on 2026-07-05/06 — cite-able calibration that the axis
works.

## A1 — Convergence order vs theoretical order

**Question**: the scheme claims order p; does the measured error scale as
h^p (with a clean log-log fit), and does it *degrade* to the expected lower
order exactly where theory says it must?

- Probe: run the same analytic case at 3-4 resolutions in the asymptotic
  regime; `order_fit(h, err)` from the metrics library; assert BOTH the slope
  band (e.g. `1.8 <= slope <= 2.3` for a second-order claim) and
  `fit.r2 >= 0.98` (a slope through non-asymptotic points is noise).
- Order *degradation* is a finding even when absolute error is tiny: a
  second-order BC that measures first-order is mis-implemented.
- Diffusive scaling: refine u together with h (u ∝ h) to hold Re constant,
  otherwise compressibility error (O(Ma²)) pollutes the fit — see A2.

**Worked example (design-ruling)**: half-way bounce-back walls are
second-order at the midpoint wall; a staircase-approximated curved wall
degrades to first order. That degradation — not any absolute error — is the
analytic case that motivates Bouzidi interpolated BB, and the audit case that
verifies Bouzidi actually restores it (see worked-example-bouzidi.md).

## A2 — Invariances

**Question**: the continuous equations are invariant under X; the
discretization is not exactly — is the breakage at the expected order, or is
something asymmetric leaking?

Field-proven sub-axes:

- **Rotation/reflection**: run the mirrored/rotated configuration and map the
  fields back; discrepancy above float noise means direction-dependent code.
  *Worked example*: **ANOM-P2-002** — the rotor blade indicator produced
  mirror arms for odd blade counts (3 blades rendered as 6 half-thickness
  arms); an along-blade sign check was missing. All fields looked plausible;
  found by cross-reading the adversarial suite's independently-written
  geometry against the implementation. Independent geometry derivation in
  tests is itself an A2 instrument.
- **Sub-cell translation**: translate the geometry by a non-integer offset
  (e.g. cylinder center +0.3 cells); physical observables (drag, flux, slip)
  must move smoothly and stay in band. Staircase BCs jump; interpolated BCs
  must not.
- **Galilean**: superpose a uniform advection velocity U0 on a known solution
  in a periodic domain; the co-moving solution must be preserved to the
  scheme's Galilean-defect order (O(Ma³) cubic lattice defect for D2Q9).
- **Diffusive scaling / same-Re-different-Ma**: two runs at the same Re but
  different Ma (u and nu scaled together) must converge to the same
  nondimensional solution as Ma → 0, with the gap shrinking as O(Ma²). A gap
  that does not shrink is a modeling error, not compressibility.

## A3 — Functional-form agreement with theory

**Question**: theory predicts not just a magnitude but a *curve* (error or
observable as a function of tau, k, distance…). Does the measurement lie ON
the curve, point by point? **"Small" is not a pass** — a small error with the
wrong functional form is a different (wrong) scheme.

- Probe: sample the parameter, `curve_agreement(theory_fn, samples, band)`.
- Canonical LBM case: the bounce-back slip law — the effective wall position
  of BB depends on tau in a known closed form. Asserting only "slip is small"
  at one tau passes a broken implementation; asserting the slip *vs tau*
  curve pins the scheme.

**Worked examples**:
- **FD shear under-report (TESTING_NOTES 2026-07-06 F5)**: post-hoc
  finite-difference shear reconstruction systematically under-reports peak
  shear where gradients are sharpest: peak −13% at u_tip 0.045, −22% at
  0.080, −35% at 0.120. Each individual number could be waved off as "a few
  percent mean error"; the *systematic trend vs gradient sharpness* is the
  functional-form finding that ruled LBM-native f_neq stress evaluation
  (FR-STRESS-01) as the reference.
- **Smagorinsky spurious nu_t (W-LES design ruling)**: Smagorinsky produces
  nonzero eddy viscosity in *pure shear* where theory says nu_t must vanish;
  WALE's operator does vanish there (nu_t → 0 near walls). The functional
  form of nu_t vs the strain invariants — not any magnitude — pre-empted the
  bug class and ruled WALE the default (REQ FR-LES-01).

## A4 — Transient fidelity

**Question**: steady states are over-constrained (everything relaxes to
them); transients expose the *bookkeeping* — half-step forcing corrections,
initialization consistency, boundary-condition time alignment.

- **Stokes' first problem** (impulsively started plate): u/U = erfc(y/2√(νt));
  L2 error must decay in time at the diffusive rate (`monotonicity`, ratio
  bands).
- **Stokes' second problem** (oscillating plate): amplitude envelope
  U·exp(−ky), phase lag k·y with k = √(ω/2ν) — `envelope_fit` + `phase_fit`.
  (Currently a SPEC-GAP in the compat API: no runtime MovingWall velocity
  setter — the ignored test carries the derivation.)
- **Acoustics**: a standing density mode must oscillate at cs = 1/√3 and damp
  at the known viscous rate — `phase_fit` for cs, `envelope_fit` for damping.
- **Startup impulse bookkeeping**: after exactly ONE step under force F, the
  momentum gain has an exact expected value (Guo: Δp = F with the F/2 moment
  correction). Single-step probes are the sharpest tool in this axis.

**Worked example**: **ANOM-P2-001** — at tau=1 TRT (Λ=3/16), the uniform
`SimConfig::force` path gains u(1) = 1.5·F (exact Guo) after step 1, but the
per-cell force-field path gains 0.9286·F — a one-time impulse deficit of
1/(2·tau_minus)·F = 4/7·F. Growth is F/step on both paths afterwards, so
EVERY steady gate (T2/T6/T11) is blind to it; the offset seeds slowly
diverging trajectories near obstacles. Found by a single-step A4 probe
crossed with A5. Pinned at the wrong ratio 7/3 in
`crates/lbm-core/tests/accuracy_audit.rs` (branch cx/acc) so the R2-C fix
must retighten it.

## A5 — Cross-path consistency

**Question**: two code paths claim to implement the same physics — do they
agree to the documented tolerance on an observable *designed to differ if
either is wrong*?

Pairs that exist in this repo: uniform force vs per-cell force field vs
gravity path; CPU scalar vs SIMD vs GPU vs MPI-partitioned (owned by
T13/T14 — do not duplicate); half-way BB vs Bouzidi at qd=0.5 (must
degenerate exactly); 2D D2Q9 vs 3D D3Q19 on a z-invariant configuration;
BC variants imposing the same physical wall.

- Cross-path probes need no analytic reference — the paths referee each
  other — which makes them cheap to write and a good *first* sweep. But pair
  them with one absolute (A1-A4) probe per family, or two consistently-wrong
  paths pass together.
- Degeneracy probes (Bouzidi qd=0.5 ≡ half-way BB; 3D with nz=1 ≡ 2D) are
  the strongest form: the expected difference is *exactly zero* (or float
  noise), so any signal is a finding.

**Worked example**: ANOM-P2-001 above was surfaced as an A5 disagreement
(uniform vs field path) and diagnosed to first order by the A4 single-step
probe. The pairing matters: A5 said "these differ", A4 said "this one is
wrong and by exactly 4/7·F".

## Cost tagging

Tag each audit row `light` (runs <~1 s, lands in the default suite) or
`heavy` (`#[ignore]`, runs under `--include-ignored` in the ~5-min full
validation tier). Convergence sweeps at 3-4 resolutions are usually heavy;
single-step probes and degeneracy checks are light. A finding that can only
be shown heavy should still get a light "canary" version (coarsest
resolution, loosest band) so the default suite carries some witness.
