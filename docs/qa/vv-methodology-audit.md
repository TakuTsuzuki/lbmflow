# V&V Methodology Audit — ASME V&V 20 / Oberkampf-Roy vs. LBMFlow Practice

**Lane 4.3 of the V&V master plan.** Read-only repo audit against the two
canonical scientific-computing V&V frameworks. Compiled 2026-07-07.

- ASME V&V 20-2009 (reaffirmed 2016), *Standard for Verification and Validation
  in Computational Fluid Dynamics and Heat Transfer* (ASME PTC/V&V Committee,
  H. W. Coleman & F. Stern principal author).
- Oberkampf & Roy, *Verification and Validation in Scientific Computing*,
  Cambridge University Press (2010).
- P. J. Roache, *Verification and Validation in Computational Science and
  Engineering*, Hermosa (1998); "Quantification of uncertainty in
  computational fluid dynamics", ARFM 29:123 (1997) — original GCI paper.
- C. J. Roy, "Review of code and solution verification procedures for
  computational simulation", JCP 205:131 (2005).

Section numbers are not quoted from the primary sources beyond what the web
search returned; only definitions are used, not paginated citations.

---

## 1. Framework summary (what we are auditing against)

### 1.1 Vocabulary (Oberkampf-Roy chapter 2; ASME V&V 20 Annex)

- **Code verification** — evidence that the discretized equations, as coded,
  solve the intended continuous mathematical model correctly, i.e. that the
  numerical scheme achieves its formal order of accuracy as *h*→0. The
  primary evidence is **observed order of accuracy** on manufactured or
  exact analytic solutions, and the presence/absence of bugs is decided by
  whether *p_obs* → *p_formal*.
- **Solution verification** — for a specific problem, estimate the
  **numerical uncertainty** *U_num* of the reported quantity of interest.
  Standard machinery: Richardson extrapolation on a three-grid family with a
  fixed refinement ratio *r*, safety-factor-guarded Grid Convergence Index.
- **Model validation** — quantitative comparison of the simulation *S*
  against experimental data *D* at a validation point, expressed as the
  **validation comparison error** *E = S − D* and the **validation
  uncertainty** *U_val*. In V&V 20's formulation
  *U_val² = U_num² + U_input² + U_D²*, where *U_input* is uncertainty of the
  simulation inputs (properties, geometry, BCs) and *U_D* is the experimental
  measurement uncertainty. The **model form uncertainty** *δ_model* is
  bracketed by |*E*| ± *U_val*.
- **Calibration** — adjustment of model parameters against data. In V&V 20 /
  Oberkampf-Roy this is emphatically **not** validation: parameters fit
  against a dataset cannot be re-validated against the same dataset. The
  discipline is a training/holdout split (or hold-out predictions for
  data never seen), and any calibrated parameter that leaks into a
  validation gate poisons the gate.

### 1.2 The GCI machinery (Roache; Roy 2005)

For three grids with cell sizes *h_1 < h_2 < h_3* and constant ratio
*r = h_{i+1}/h_i* (typically 2), and a quantity *φ_i*:

1. Observed order of accuracy: *p_obs = ln((φ_3 − φ_2)/(φ_2 − φ_1)) / ln r*
   (or a nonlinear solve for non-constant *r*).
2. Richardson-extrapolated exact value: *φ_h=0 ≈ φ_1 + (φ_1 − φ_2)/(r^{p_obs} − 1)*.
3. Grid Convergence Index (fine-grid discretization error band):
   *GCI_fine = F_s · |φ_1 − φ_2| / (r^{p_obs} − 1) / |φ_1|*
   with safety factor *F_s = 1.25* when three grids are used and
   *p_obs ≈ p_formal* (asymptotic range confirmed), *F_s = 3.0* otherwise
   (two grids, or *p_obs* materially off).
4. Asymptotic range check: the three-grid consistency *(φ_3 − φ_2)/(φ_2 − φ_1) ≈ r^{p_obs}*
   must actually hold; if it doesn't, the mesh isn't in the asymptotic range
   and the GCI is not an uncertainty band, it's a stronger warning.

*U_num* on a headline result is then reported as, e.g., *Cd = 5.5795 ±
GCI_fine%*.

### 1.3 What V&V 20 does *not* say

It doesn't require any particular test set. It requires that headline
quantities carry stated uncertainty, that validation compares like-for-like
(*E* vs. *U_val*), and that calibration is not laundered into validation.

---

## 2. LBMFlow practice — how it maps to the framework

### 2.1 Where practice CONFORMS to V&V 20 / Oberkampf-Roy

| # | Practice | Framework equivalent | Evidence |
|---|---|---|---|
| C1 | Analytic/exact-solution gates T1 (Taylor–Green), T2 (Poiseuille), T3 (Couette), T15.1 (z-invariant TGV bit-match), T15.2 (rectangular-duct exact series). | **Code verification** — order-of-accuracy against known analytic solutions. Best-practice channel. | `VALIDATION.md` T1/T2/T3/T15; test files `crates/lbm-core/tests/validation_*.rs`, `t15_3d.rs`. |
| C2 | Observed-order-of-accuracy asserts, e.g. T1 `order = log2(err32/err64) ≥ 1.7`, T2 BGK `H=8→16 ≥ 1.7`, T15.4 3D TGV `N=32→64 ≥ 1.7`, T8 Bouzidi `≥ 1.7` from three D∈{10,20,40} grids. | **Observed order of accuracy** in the Roy/Roache sense. This is the correct verification metric. | `smoke_poiseuille.rs:85`, `validation_channel.rs:113`, `t15_3d.rs:683`, `d5_long_horizon.rs:153`, `validation_cylinder.rs:424`. |
| C3 | T8 Bouzidi three-grid extrapolation: `extrapolated_limit = cd40 + (cd40 − cd20)/(2^p_obs − 1)`, then band-check the extrapolated limit against the Schäfer-Turek reference. | **Richardson extrapolation** in exactly Roache's form, with observed order (not assumed formal order). | `validation_cylinder.rs:425`. |
| C4 | Cumulant `h²-intercept` audit: fit `nu_eff/nu − 1 = a + b/N²` across N ∈ {24, 32, 48} to separate the resolution-independent closure offset *a* from the O(*h²*) spatial-error floor *b/N²*. Band on *a* alone. | **Solution verification separating model error from discretization error** — the classical error-budget decomposition. This is textbook V&V 20 thinking. | `accuracy_audit_cumulant.rs:285-347`; ANOM-P4-008 in `anomaly-log.md`. |
| C5 | T15.5 3D cavity band derivation: explicit budget of reference uncertainty (~1e-5, negligible) + grid bias (0.005–0.010 U from Zhang series) + interpolation error (~0.003 U) + O(Ma²) (~0.003 U), summed to an expected ~0.008–0.015 U, then a 2–3x margin. | **Informal *U_val* budget**: separates *U_num* (grid+Ma²+interp) from *U_D* (reference uncertainty). One of the few places we actually decompose in the V&V 20 style. | `T15_5_CAVITY3D_REFERENCE.md` §6.3. |
| C6 | The +0.0025 D3Q19 cumulant "viscosity offset" was **removed** on the finding that it was a calibration to N=32, and a training/holdout split was constructed: training = TGV3D at calibration settings; holdouts = advected TGV3D (Galilean frame sweep), off-Re TGV3D, D3Q19-vs-D3Q27 cross-check. | **Correct handling of a calibration-vs-validation conflation.** This is exactly the "don't validate against the training set" doctrine, executed after the fact with a real holdout suite. | `PHYSICS.md` 2026-07-07 ANOM-P4-008 entry; `crates/lbm-core/tests/cumulant_holdout.rs`. |
| C7 | Physics discipline ban list explicitly forbids (a) constants calibrated to pass an acceptance band, (b) case-identity branches, (c) transported-quantity clamps that silently absorb error, (d) decorative physics terms. Enforced by grep-sweep (V&V lane 7.1). | **Prohibition of the calibration-hidden-in-validation pattern**, at the source level. Nothing in V&V 20 asks for this, but the intent — "don't quietly fit and then call it validated" — is the same. | `CLAUDE.md` "Working discipline"; `.claude/skills/lbmflow-physics-discipline`. |
| C8 | Adversarial-test principle: validation tests written by orders separate from the implementation, from the spec and public API only. Includes mutation testing (V&V lane 2.1) already killing 3 seeded mutants. | **Independent verification** in the Oberkampf-Roy sense (the coder does not write the tests). Not required by V&V 20 but consonant with it. | `VV_MASTER_PLAN.md` axes 2 and 3. |
| C9 | Backend equivalence (T13 partition invariance bit-exact; T14 CPU-vs-GPU bit/threshold; T14 GPU-absolute-physics after R-Phase 2). | **Cross-platform code verification** — same code, different execution paths, must agree. | `t13_split_invariance.rs`, `t14_backend_equiv.rs`, `gpu_absolute.rs`. |
| C10 | T17/VR-STR band governance rule: "tighten freely; loosen only with PHYSICS.md rationale (reference uncertainty / method order / resolution limit)". Every loosening carries an error-source justification. | **Error-source attribution before band change** — a working substitute for a formal *U_num*/*U_input*/*U_D* form. | `VALIDATION.md` §T17 "Band governance". |

### 2.2 Where practice DEVIATES from V&V 20 / Oberkampf-Roy

| # | Gap | What V&V 20 or Oberkampf-Roy prescribes | Severity | Cheapest fix |
|---|---|---|---|---|
| G1 | **No headline quantity carries a stated numerical-uncertainty band *U_num*.** Bands are asserted (e.g. Cd ∈ [5.35, 5.85], k_La ±25%, Np ±10%) but reported as *pass/fail thresholds*, not as *simulation result ± U_num*. | ASME V&V 20 requires that any validation comparison report the simulation-side numerical uncertainty *U_num* as a number, together with *U_input* and *U_D*, so that *U_val = √(U_num² + U_input² + U_D²)* is meaningful and *E vs. U_val* is the actual verdict. Pass/fail against a fixed band is a weaker artifact. | **Medium — the concepts are used inline (C4, C5, T17 governance) but never assembled at the level of "the number we report".** | Add a per-headline-result JSON line to the manifest emitted by presets/gallery: `{quantity, value, U_num, U_num_method, U_num_grids?}` for Cd (T8), Np (VR-STR-01), k_La (VR-STR-04), extrema (T15.5). Populate for cases that already have three grids or a defensible one-off error model; leave `U_num: null, U_num_method: "none"` and *say so* for cases that don't. This is Oberkampf-Roy's rule: "compute it or say you didn't". No new tests required; it's a reporting-layer change. |
| G2 | **No Grid Convergence Index reported anywhere.** T8 Bouzidi does the three-grid Richardson step and prints `extrapolated_limit`, but stops short of computing `GCI_fine = F_s · |φ_1 − φ_2| / (r^{p_obs} − 1) / |φ_1|` with Roache's safety factor. Everywhere else that uses `log2(errN/err2N)` is two-grid, so the observed order is measured but the asymptotic-range check `(φ_3 − φ_2)/(φ_2 − φ_1) ≈ r^{p_obs}` isn't performed. | Roache / Roy prescribe three grids, computed *p_obs*, asymptotic-range consistency check, then *F_s = 1.25* if *p_obs ≈ p_formal* else *F_s = 3.0*, then quote the resulting GCI band as the numerical uncertainty. This is the V&V 20 default machinery for *U_num*. | **Medium.** The audit's own conclusions rest on three-grid extrapolation (C4) and error-source decomposition (C5), but they aren't standardised. Where three grids exist we're one arithmetic step away from a proper GCI; where only two grids exist we can only quote *F_s = 3.0* which is honest. | Add a shared helper `common/gci.rs` with `pub fn gci(fine, medium, coarse, r) -> Gci { p_obs, phi_extrap, gci_fine, asymptotic_range_ok }` (asymptotic check = `(medium−coarse)/(fine−medium)` within 20% of `r^p_obs`). Wire it into (a) T8 Bouzidi (already three grids), (b) T15.4 3D TGV (add N=48 to the existing 32/64 pair to make three grids), (c) T15.2 duct where three widths already exist. Anywhere we can't afford three grids, document `F_s = 3.0` two-grid fallback per Roache. |
| G3 | **`order = log2(errN/err2N)` is used with only two grids in T1, T2, T15.4, and D5.** Roache/Roy: two-grid observed order is valid only if you *assume* the formal order and check consistency; otherwise you need three grids to *measure* *p_obs*. Our current form measures a slope on two points, which coincides with *p_obs* only under the constant-ratio assumption. | The two-grid form is fine as a screen; it isn't a solution-verification statement. To claim "the scheme is second-order in this regime" you need three grids. | Low (screening use is honest as long as it's labelled as such). | Either add a third resolution (usually the cheapest, and adds real information), or rename the metric in test names/output to `two_grid_slope_screen` and reserve `observed_order_of_accuracy` for the three-grid variant. |
| G4 | **The V&V 20 uncertainty framework (*U_num*, *U_input*, *U_D*, *U_val*) is not the vocabulary in use.** Local vocabulary is "SPEC-GAP", "STOP-RULE", "frozen band", "PHYSICS.md rationale", "training/holdout". The concepts map — SPEC-GAP ≈ documented model-form gap; STOP-RULE ≈ hitting `E > U_val` and refusing to report a false pass; frozen band ≈ threshold, not uncertainty; PHYSICS.md rationale ≈ error-source attribution — but the mapping isn't stated anywhere. | Not a defect, a translation gap. Downstream users familiar with V&V 20 don't know we're doing much of it under different names. | Low. | Add a one-page "V&V vocabulary crosswalk" section to `PHYSICS.md` or `VALIDATION.md` §0 that maps our terms to V&V 20: SPEC-GAP → known-*δ_model* domain restriction; STOP-RULE → |E| > U_val decision; frozen-band → engineering acceptance threshold, not an uncertainty band; training/holdout → calibration-validation separation. |
| G5 | **`U_input` is not itself audited.** Simulation inputs — ν, τ, u_max, geometry-to-lattice mapping, D (nominal vs. hydrodynamic), Ma — carry real uncertainty when we compare to Schäfer-Turek, Ghia, Albensoeder-Kuhlmann. The T15.3 sphere-drag "hydrodynamic pair (D+1)/2" normalization is an explicit acknowledgement that the geometric input has O(1/D) uncertainty in staircase, but this is not propagated. | V&V 20's *U_input* is meant to be estimated (typically Monte-Carlo, or a linearised sensitivity study) and combined into *U_val*. | Low today, higher as we approach VR-STR-01 (aerated stirred tank) where the geometric inputs (impeller clearance, blade thickness) are physically uncertain, not just discretisation-uncertain. | For the current headline cases, add a sensitivity ledger: for T8, ±1 cell in *D* and ±5% in inlet-profile amplitude → measured ±Δ*Cd*; for VR-STR-01, ±1 cell impeller clearance → ±Δ*Np*. This is a bench measurement, one small helper, no framework overhaul. |
| G6 | **`U_D` (reference uncertainty) is inconsistently treated.** T15.5 does it right (§6.3 explicit ~1e-5 reference band, negligible). Ghia (T7), Schäfer-Turek (T8), Schiller-Naumann (T15.3), MKM DNS (T17/VR-STR-03) are treated as truth without a reported *U_D*. Ghia's tabulated points have several-percent scatter with later references; Schäfer-Turek gives error bars in the original paper. | V&V 20 requires *U_D* on the experimental (or reference) side, and specifically warns against treating a canonical reference as if *U_D = 0*. | Low-to-medium. Our bands are wide enough that this is usually not load-bearing, but the T15.5 tightening candidate and VR-STR-03 DNS comparison will hit this quickly. | For each canonical reference in `VALIDATION.md`, add a one-line "reference uncertainty" note (source scatter, published error bar, or "not stated in the original — treat as 0"). Where a validation gate ends up assuming *U_D = 0* against a reference that has stated error bars, promote to a note and consider widening the band by *U_D* on next retighten. |
| G7 | **No standing check for silent calibration in the T11 (Shan-Chen), T11b/c (contact angle), and T17 (VR-STR "frozen after impl") tracks.** T11 explicitly says "coexistence bands frozen to measurement" and T11b/c freeze angles at ±8° to what the current implementation happens to produce. That is calibration-as-regression-pin. It is *labelled* honest (band-vacuity scan §2.2 and T17 rev.4 band governance), so it is not the +0.0025 pattern — but there is no positive holdout for these subsystems. | Oberkampf-Roy: a regression pin against your own measurement is not validation, it's a change-detector. Fine as long as no gate uses it to claim the physics is right. | **Medium** at the process level (not the code level). The ban is soft; a new subsystem could import the pattern under "frozen to measurement" and not attract review. | Extend `.claude/skills/lbmflow-physics-discipline` with an explicit "regression-pin vs. validation-gate" distinction. Every `frozen at measurement` band gets a one-line label in the test — `// REGRESSION_PIN` vs. `// VALIDATION_GATE` — and the physics-discipline grep sweep counts regression pins that lack a linked holdout. Also: the T11 track should acquire at least one *external* validation (Laplace slope σ vs. any published SC parameter fit) as a positive step. |
| G8 | **The band-vacuity scan (V&V lane 2.2, `band-vacuity-scan.md`) is not run per landing.** The 2026-07-06 sweep found asserts loose by ~10× to ~10¹²× against measured values. This is the exact failure mode Oberkampf-Roy call out: a "test" whose band is so loose that any implementation passes. | The framework needs *effective* tests; a band 10¹² wider than the measured error is a token. | Low today (the surviving retighten queue is short), but drifts up as the codebase grows. | Fold band-vacuity into CI: for every test with a numeric assert, dump `(measured, band)` and require `band < 50 × measured` unless the assert carries a `// LOOSE_BY_DESIGN reason=...` tag. Cheap grep sweep. |
| G9 | **No formal "asymptotic range confirmed" statement per convergence claim.** T1/T2/T15.4 assert *order ≥ 1.7*; we don't check that we're inside the asymptotic range for the reported grid family. If a grid family shows *p_obs* = 1.4 (below asymptotic-range expectation for second-order LBM), we currently just fail the band; V&V 20 would require flagging "not in asymptotic range → GCI unreliable → *F_s = 3.0*". | Roache: a lone *p_obs* below the formal order is not "the code is broken", it's "we're not asymptotic yet; grow the grid or trust nothing tighter than *F_s = 3.0*". | Low. | Once G2's `gci` helper exists, its `asymptotic_range_ok` flag settles this. Message on failure: distinguish "wrong order" from "not asymptotic". |
| G10 | **Model-form uncertainty *δ_model* is never reported as a number.** Oberkampf-Roy: after all *U*s are collected, |E| − *U_val* (when positive) is a lower bound on the model-form uncertainty. LBMFlow reports pass/fail and PHYSICS.md prose; the residual error attributable to the choice of collision operator, forcing scheme, subgrid model, etc. is never quoted. | Not required by V&V 20 for a pass/fail gate; expected in a validation report suitable for a customer. | Low today; **medium the day the paper ships or a customer asks "how good is your k_La?"** — that's precisely the number they want. | Add a "model-form gap" column to the T17 acceptance table (currently only *target/tolerance*). Populate it from the observed |E| minus the current U_num budget on VR-STR-01/02/04, once the impl lands. |

---

## 3. The +0.0025 case in V&V 20 terms

The removed D3Q19 cumulant offset is the audit's cleanest specimen and worth
naming clearly because it will keep coming up.

- **What it was.** A constant added to the shear-relaxation rate that
  cancelled the ordinary O(*h²*) spatial-discretisation error of TGV3D
  decay at exactly N=32, chosen so the acceptance band passed. Under a
  V&V 20 accounting: it moved a residual from *U_num* (which V&V 20
  requires to shrink under refinement) into *δ_model* (which does not).
  Because the acceptance gate compared simulation to the analytic decay
  rate at a single N, the calibration was invisible to the gate.
- **How it was caught.** The `nu_eff/nu − 1 = a + b/h²` fit across three
  grids (C4) — a Richardson-style resolution separation. The intercept *a*
  measured the operator's continuum bias; the slope *b* measured the
  spatial-error floor. The offset's own footprint −0.0025·2/(2−ω) matched
  the poisoned intercept to 99.8%, exposing that the "correction" was
  a resolution-point fit.
- **What replaced it.** (i) The offset was removed. (ii) The acceptance
  gate was rewritten as an *h²*-intercept band on *a* (`accuracy_audit_cumulant.rs`
  E2), which is invariant across grids. (iii) A training/holdout split
  was constructed (`cumulant_holdout.rs`): advected TGV3D (Galilean-frame
  sweep), off-Re TGV3D, D3Q19-vs-D3Q27 cross-check. The advected-frame
  holdout in turn found a frame-dependence signal that the calibrated
  form had hidden — a real *δ_model* finding disclosed in PHYSICS.md as a
  currently-open closure gap.

That is a textbook Oberkampf-Roy episode. It also demonstrates that the
missing framework machinery (formal GCI, stated *U_num* on the headline
result) would have caught this at first landing, not on audit.

## 4. Is anything similar still lurking?

Systematic search for the same pattern — constants calibrated against an
acceptance observable, used inside a gate against the same observable.

| Candidate | Verdict | Reasoning |
|---|---|---|
| Cumulant `−0.16 |u|²` velocity term | **Open (labelled).** PHYSICS.md 2026-07-07 explicitly states "not validated as a viscosity correction; E1 remains SPEC-GAP until rerun with the ablation flag", and there is a `CENTRAL_MOMENT_DISABLE_VELOCITY_CORRECTION_FOR_ABLATION` toggle. Not laundered; disclosed as a `δ_model` gap. Ablation study queued. | Same family as +0.0025; hasn't been fully audited yet because it needs the toggle. Not a hidden hack, an open one. |
| T11 flat-interface coexistence ρ_l = 1.888 ± 2%, ρ_v = 0.1194 ± 3% | **Regression pin, not validation.** Frozen to measurement; SC coexistence is known to depend on the forcing scheme (Li/Luo/Li 2012), so treating measured coexistence as a spec is legitimate as *regression* but not as *validation*. | The band-vacuity scan already labels this correctly. **G7** applies. Positive holdout candidate: Maxwell-rule Δρ predicted vs. measured under Shan's F/ρ form (published SC parameter). |
| T11b contact-angle frozen G_w=−1.5:133.2°, 0:160.4°, +1.5:163.7° (±8°) | **Regression pin, not validation.** The doc says as much: "Because this implementation uses ψ=0 for solid + a separate −G_w ψ Σ w s c cohesion, G_w=0 is not 90°". | Same as T11 — pinned to what THIS implementation gives, no independent reference. **G7** applies. |
| T15.5 extremum band widened 6% → 13% at N=72 | **Not a calibration** — the widening carries a PHYSICS.md rationale in the V&V 20 style: reference scatter + measured convergence direction. Loosening was done under band governance with a stated error-source justification. | Textbook error-source attribution. |
| T9b convective outflow "frozen pressure ratio 7.96" | **Regression pin, correctly labelled.** T9b's spec: "Advantage is geometry-dependent — at minimum require stability, non-divergence, reverse flow ≤ 5%. Compare to Outflow with the same geometry and T9 metrics; freeze at measured values." | Pinned to detect regressions in the mass-consistency pin, not to claim outflow physics is validated. Correctly scoped. |
| T7 Ghia Re=400 v(0.9063)=−0.23827 excluded from RMS | **Documented reference error**, not a LBMFlow calibration. Ghia's published typo, excluded per PHYSICS.md 2026-07-05. | Not a lurking issue; exactly how Oberkampf-Roy say to handle a reference-side defect. |
| T17 rev.4 "loosen only with PHYSICS.md rationale" plus a track of frozen bands "after impl" | **Process risk, not yet a specific defect.** Every "frozen after impl" is a future calibration opportunity if the impl author writes their own band without a holdout. | **G7** applies; the fix is process (`REGRESSION_PIN` vs. `VALIDATION_GATE` labels), not code. |
| T18.1 mass ledger 1e-6 rel, T18.2 mass drift 5e-9 rel | **First-measurement band with derivation** — the T18.1 band is derived analytically from summation round-off (N·ε·M/|Σq|), not fit. | Not a calibration; a bound. |
| ANOM-P4-011 (F19/F20) refuted by derivation | **Correctly refuted.** Cold-review triage derived the expected artifact and confirmed it was model-form, not a bug. | This is Oberkampf-Roy's "derive before blaming the engine" rule, executed. |

**Net.** One open closure gap (−0.16 |u|² velocity term) is being handled
under disclosure. The T11 and T11b/c frozen bands are process-hygiene risks
(G7) rather than active calibration hacks. Nothing else looks like the
+0.0025 pattern under this pass. The dominant residual is systemic (no
GCI, no U_num on the headline, no crosswalk vocabulary) rather than any
specific hidden fit.

## 5. Priority ranking of the fix proposals

Cheapest impactful first. Each has an anchor to §2.2 above.

1. **[G1 / G2] Add `common/gci.rs` and report *U_num* on 3 headlines** (Cd
   T8 already three-grid; sphere Cd T15.3; TGV3D order T15.4 with an added
   N=48). One helper, three call sites, changes the reporting layer only.
   Delivers the missing V&V 20 headline that the +0.0025 case exposed as a
   real gap.
2. **[G4] V&V vocabulary crosswalk** in `VALIDATION.md` or `PHYSICS.md`.
   Documentation only; removes the "we're doing V&V 20 in local dialect"
   friction.
3. **[G7] `REGRESSION_PIN` vs. `VALIDATION_GATE` labels + T11 external
   holdout.** Discipline update in the physics-discipline skill; adds one
   real T11 validation (Maxwell-rule Δρ prediction) to break the
   pinned-to-self pattern. Codex-scale order.
4. **[G3] Rename two-grid `order` to `two_grid_slope_screen` where used
   as a screen; add a third grid where it's used as a claim.** T1/T2/T15.4
   go from "screen" to "verification" for the cost of one extra resolution
   each.
5. **[G8] Fold band-vacuity into CI.** One script + a `LOOSE_BY_DESIGN`
   tag; prevents future band-erosion silently.
6. **[G6] Reference-uncertainty notes in `VALIDATION.md`.** Line-item
   documentation.
7. **[G5] Input-sensitivity ledger for T8, VR-STR-01.** Bench measurement,
   defers formal *U_input* propagation until VR-STR-01 lands.
8. **[G10] `δ_model` column on T17 acceptance table.** Populate as VR-STR
   subsystems land.
9. **[G9] "Not-asymptotic" flag** falls out of G1/G2's helper automatically.

## 6. Deliverable

Everything above is process/reporting. **No physics is being alleged
wrong by this audit** beyond what PHYSICS.md and anomaly-log.md already
report. The audit's job is to name the framework LBMFlow is already
partially executing, and to point at where it stops short of what a
customer with a V&V 20 procurement clause would ask for.

The +0.0025 episode was the framework catching itself: three-grid resolution
separation exposed a resolution-point calibration, holdouts were
constructed after the fact, PHYSICS.md recorded the ban rule. Formalising
the machinery (G1/G2) is the change that lets the framework catch the
next such episode before it lands, not on audit.
