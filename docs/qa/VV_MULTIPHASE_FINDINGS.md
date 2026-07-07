# V&V-MULTIPHASE Findings

Lifecycle: snapshot (2026-07-06) — dated audit record, never edited.

Date: 2026-07-06
Branch: `cx/vv-multiphase`
Base commit audited: `4eac49e2aebccfa52283c29cd7f22aabba3120f5`

## Scope

Suborder F audited the existing Shan-Chen and two-component multiphase validation
surface:

- Shan-Chen coexistence density and pressure equilibrium.
- Laplace law and droplet pressure jump.
- Contact angle through `G_w` and virtual wall density `wall_rho`.
- Rayleigh-Taylor growth and two-component separation threshold.
- Behavior anchors: droplet radius/centering, pressure plateaus, spurious-current
  location, mass drift, and mirror symmetry.

No solver physics was changed. The only code edits add QA behavior checks that
consume already-exported fields.

## Coverage Status

| Area | Current evidence | Status | Gap / downgrade |
|---|---|---|---|
| T11 flat Shan-Chen coexistence | `validation_multiphase::t11_flat_interface_coexistence_pressure_currents_and_mass` passed in this session. Checks liquid/vapor densities, SC-EOS pressure balance, max spurious speed, and total mass drift. | VALIDATED for 2D CPU/f64 facade scalar bands | No visual artifact was produced for the flat-interface test in this session. Do not treat this as behavior-reviewed beyond the scalar assertions. |
| T11 f32 flat interface | `validation_multiphase::t11_f32_flat_interface_stays_finite_and_close_to_f64` passed. | VERIFIED/VALIDATED for stated f32 flat-interface smoke | Only flat-interface f32 is covered; no f32 Laplace/contact/RT atlas. |
| T11 Laplace droplet | Rust single-radius smoke passed. QA `droplet-t11` produced `sigma = 0.03488`, 5.1% above the frozen `0.0332`, within the 15% single-droplet band. | VALIDATED for representative single-radius 2D CPU case | Full four-radius linearity is still `#[ignore]`; run before making the full slope claim from this branch. |
| Droplet behavior anchors | QA `droplet-t11`: `R_fit=22.26` for `R0=20`, centroid `(64.00,64.00)`, pressure plateau spans below Laplace jump, max spurious current at `d/R_fit=1.009` on the half-density interface, mirror asymmetry `8.95e-12`. | VALIDATED for representative artifact run | This is not a full radius atlas. |
| T11b `G_w` contact angle | Rust default test passed with printed angles: `G_w=-1.5 -> 133.191 deg`, `0 -> 160.435 deg`, `+1.5 -> 163.740 deg`. | VALIDATED for the frozen `G_w` characteristic | The known wetting-range limitation remains; `G_w=0` is not a 90-degree neutral wall. |
| T11c `wall_rho` full range | Rust film case `wall_rho=1.6` passed. QA `droplet-on-wall-t11c` at `wall_rho=1.0` measured `theta=64.5 deg`, within the 63 +/- 8 deg band. | PARTIAL | The full `wall_rho={0.3,0.6,1.0}` 30k sweep is `#[ignore]` and was not run here. Mark full-range monotonicity BENCH-PENDING for this branch. |
| T12 two-component separation | `validation_rt::t12_mcmp_separation_threshold_smoke` passed for `G_ab=2.2` separates and `G_ab=1.8` mixes, with mass drift within band. | VALIDATED for default separation smoke | No visual artifact for separation in this session. |
| T12 sigma_AB | `validation_rt::t12_mcmp_sigma_ab_laplace_regression` passed and printed `sigma_AB = 2.86969302e-2`. | VALIDATED for scalar sigma regression | No artifact generated for the AB droplet. |
| T12 Rayleigh-Taylor | Default smoke passed: `max_amp=23.853410`, mass drift `4.6875e-13`. | VERIFIED behavior exists; full validation BENCH-PENDING | The full `256^2 x 12k` growth-rate fit is `#[ignore]` and was not run here. Do not claim the full `gamma_fit/gamma_th` gate from this session. |
| Backend/lattice coverage | Multiphase scenario validation warns/rejects unsupported GPU/3D paths; T13 has Shan-Chen split-invariance coverage. | VERIFIED-ONLY for partition equivalence; SPEC-LIMITED for GPU/3D | Multiphase physics validation here is 2D CPU/facade. No GPU or 3D multiphase physics claim is allowed. |

## Artifacts

Representative artifact run:

```bash
python3 scripts/qa/run_sweep.py --bin target/release/lbm --out out/vv/multiphase --only droplet-t11,droplet-on-wall-t11c
```

Machine-readable summary:

- `out/vv/multiphase/results.json`

Free droplet T11 artifacts:

- `out/vv/multiphase/droplet-t11/scenario.json`
- `out/vv/multiphase/droplet-t11/validate.json`
- `out/vv/multiphase/droplet-t11/manifest.json`
- `out/vv/multiphase/droplet-t11/rho_30000.png`
- `out/vv/multiphase/droplet-t11/rho_5000.csv` through `rho_30000.csv`
- `out/vv/multiphase/droplet-t11/speed_30000.csv`

Wall-contact T11c artifacts:

- `out/vv/multiphase/droplet-on-wall-t11c/scenario.json`
- `out/vv/multiphase/droplet-on-wall-t11c/validate.json`
- `out/vv/multiphase/droplet-on-wall-t11c/manifest.json`
- `out/vv/multiphase/droplet-on-wall-t11c/rho_30000.png`
- `out/vv/multiphase/droplet-on-wall-t11c/rho_30000.csv`

## Behavior Review

### 2026-07-06 behavior review - `droplet-t11`

Pattern: A circular dense droplet remains centered in vapor; SC-EOS pressure is
higher inside than outside; the largest spurious current is localized on the
diffuse interface.

Mechanism: Shan-Chen cohesion creates a surface-tension pressure jump balanced by
the curved diffuse interface; parasitic velocity appears where the force stencil
samples the largest density gradient.

Resolved vs closure: The LBM transport, Guo forcing, and mass conservation are
resolved by the solver. Surface tension and coexistence are Shan-Chen
pseudopotential closure behavior, validated only over the frozen T11/T12 domain.

Artifacts checked: `rho_30000.png`, `rho_30000.csv`, `speed_30000.csv`,
`manifest.json`, and `results.json`. No wall, outlet, or clamp was active in this
periodic droplet run.

Verdict: PHYSICAL within the validated Shan-Chen closure domain.

Routing: none.

### 2026-07-06 behavior review - `droplet-on-wall-t11c`

Pattern: The `wall_rho=1.0` wall-attached droplet wets the bottom wall and the
spherical-cap fit reports `theta=64.5 deg`, inside the frozen T11c band.

Mechanism: The virtual wall density contributes solid-neighbor pseudopotential
cohesion, increasing wetting as `wall_rho` approaches the liquid density.

Resolved vs closure: The wall is a half-way bounce-back boundary. The contact
angle is closure-driven by the virtual-wall-density Shan-Chen wall model and is
valid only over the measured frozen cases.

Artifacts checked: `rho_30000.png`, `rho_30000.csv`, `manifest.json`, and
`results.json`. The contact row is occupied; no edge-ring accumulation or outlet
reflection was present.

Verdict: PHYSICAL within the frozen T11c virtual-wall-density closure domain.

Routing: none.

## Findings

- **F-MP-01, resolved in this branch:** the QA multiphase droplet check had scalar
  spurious-current and Laplace checks but lacked behavior anchors for radius,
  pressure plateau, and current location. Added checks in `scripts/qa/qa_checks.py`
  and wired them into `droplet-t11`.
- **F-MP-02, validation gap:** full T11 Laplace linearity, T11c full-range
  wall-rho monotonicity, and T12 full RT growth-rate fit remain `#[ignore]`.
  They are BENCH-PENDING unless `cargo test --release -- --include-ignored` is run
  and completes green.
- **F-MP-03, claim downgrade:** multiphase validation here is 2D CPU/facade
  validation. GPU and 3D multiphase are unsupported product paths and must not be
  advertised as validated.
- **F-MP-04, scope boundary:** Shan-Chen validation does not validate the M-F
  high-density-ratio resolved phase-field requirements in `REQ_STIRRED_REACTOR.md`.
  Those remain SPEC-ONLY under T17 until W-VOF is implemented and validated.

## Merge Recommendation

Merge recommended: the normal workspace release test was green in this session.
The branch improves V&V observability without touching solver physics. It should
not be used to upgrade any GPU, 3D, or high-density-ratio phase-field claim.
