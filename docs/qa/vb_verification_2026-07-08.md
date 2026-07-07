# VB-01..VB-08 verification pass (2026-07-08)

Lifecycle: snapshot (adversarial verification on branch
`cx/vb-verify-2026-07-08`).

Verification scope: landed implementation as of `git log --oneline`
showing BCFD-000..084 and BCFD-100..102 present. BCFD-046, BCFD-047,
BCFD-048, BCFD-090, BCFD-091, and BCFD-092 were not present in the
checked log.

Summary:
- Un-ignored tests: 8
- Still ignored: 14
- Anomalies filed: 5 (`BCFD-VV-001`..`BCFD-VV-005`)
- Acceptance bands weakened: no
- Test setup notes: VB-07 adds a `1.0e-14` roundoff floor only to the
  error-vs-error refinement comparison for an analytically constant
  shear integral; the 5% and 1% exposure bands are unchanged.

## VB-01 Single-phase stirred tank Np

Impl landed: yes (`BCFD-030`, `BCFD-031`).
Test count: 3 (un-ignored: 0, still ignored: 3, anomalies: 1)
Un-ignored tests + result: []
Anomalies:
- `BCFD-VV-004`: Test-side disposition: keep all three tests ignored because
  the adversarial skeleton requires real 3D stirred-tank reference runs at
  the published operating points, and no public validation driver/reference
  fixture is wired into the test. Impl-side disposition: BCFD-030/031 landed
  the runner/QOI pieces, but not a public reproducible VB-01 validation
  harness that yields statistically stationary Np for Rushton/PBT and the
  three-grid sequence.
Disposition: Not-runnable-until-BCFD-VV-004

## VB-02 Passive-scalar mixing

Impl landed: yes (`BCFD-034`, `BCFD-035`; depends on landed `BCFD-030`).
Test count: 2 (un-ignored: 0, still ignored: 2, anomalies: 1)
Un-ignored tests + result: []
Anomalies:
- `BCFD-VV-005`: Test-side disposition: keep both tests ignored because the
  skeleton requires a real stirred-tank point-pulse ADE run, not synthetic
  CV reducer data. Impl-side disposition: the scalar/mixing QOI API exists,
  but no public VB-02 stirred-tank scalar validation run is wired into the
  test surface.
Disposition: Not-runnable-until-BCFD-VV-005

## VB-03 Wall-shear and shear-rate fields

Impl landed: yes (`BCFD-032`, `BCFD-033`).
Test count: 3 (un-ignored: 3, still ignored: 0, anomalies: 0)
Un-ignored tests + result:
- `couette_gamma_dot_matches_u_over_h_at_n64`: PASS
- `poiseuille_gradient_matches_signed_analytic_profile_at_n64`: PASS
- `shear_fields_converge_second_order_on_32_64_128`: PASS
Anomalies: []
Disposition: Engineering-GREEN

## VB-04 Phase-field droplet and Laplace law

Impl landed: no for full VB-04 (`BCFD-040`..`BCFD-043` landed; `BCFD-048`
not landed; Bundle P not merged).
Test count: 3 (un-ignored: 0, still ignored: 3, anomalies: 0)
Un-ignored tests + result: []
Anomalies: []
Disposition: Not-runnable-until-BCFD-048

## VB-05 Sparger gas ledger

Impl landed: no (`BCFD-046`, `BCFD-047` not landed; Bundle P not merged).
Test count: 3 (un-ignored: 0, still ignored: 3, anomalies: 0)
Un-ignored tests + result: []
Anomalies: []
Disposition: Not-runnable-until-BCFD-046/047

## VB-06 Oxygen kLa synthetic

Impl landed: yes (`BCFD-050`, `BCFD-051`, `BCFD-052`).
Test count: 2 (un-ignored: 1, still ignored: 1, anomalies: 1)
Un-ignored tests + result:
- `dynamic_gassing_fit_recovers_input_kla_with_high_r2`: PASS
Anomalies:
- `BCFD-VV-001`: Test-side disposition: keep
  `equilibrium_case_fits_zero_kla` ignored to avoid committing a red gate.
  Impl-side disposition: `dynamic_gassing_kla_fit` returns a skipped QOI for
  steady `C = C*`, while VB-06 acceptance says the equilibrium case fits
  `kLa approximately 0`.
Disposition: Impl-anomaly

## VB-07 Cell shear-exposure integral

Impl landed: yes (`BCFD-060`, `BCFD-061`).
Test count: 3 (un-ignored: 2, still ignored: 1, anomalies: 1)
Un-ignored tests + result:
- `above_threshold_couette_exposure_matches_analytic_at_coarse_and_halved_dt`: PASS
- `below_threshold_couette_exposure_is_exactly_zero_for_every_tracer`: PASS
Anomalies:
- `BCFD-VV-002`: Test-side disposition: keep
  `percentile_reducer_matches_synthetic_distribution` ignored because the
  adversarial skeleton expects nearest-rank percentiles. Impl-side
  disposition: the landed reducer uses linearly interpolated percentiles.
  Follow-up must either freeze the interpolated method in the VB acceptance
  text or change the reducer to the adversarial nearest-rank contract.
Disposition: Impl-anomaly

## VB-08 Synthetic scale-up decision

Impl landed: yes (`BCFD-083`, `BCFD-084`).
Test count: 3 (un-ignored: 2, still ignored: 1, anomalies: 1)
Un-ignored tests + result:
- `evaluator_recovers_analytic_large_tank_feasible_set_from_synthetic_maps`: PASS
- `infeasible_case_emits_explicit_conflict_table`: PASS
Anomalies:
- `BCFD-VV-003`: Test-side disposition: keep
  `tightest_constraint_and_default_priority_are_reported` ignored to avoid a
  red gate. Impl-side disposition: `evaluate_operating_window` ranks
  constraints by violation magnitude, while VB-08 requires the documented
  default priority order `constant kLa -> P/V -> tip speed -> mixing time`
  unless user weights override. The public constraint set also has no tip
  speed constraint field.
Disposition: Impl-anomaly
