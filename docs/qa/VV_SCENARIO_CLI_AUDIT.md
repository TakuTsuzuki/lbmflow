# V&V Scenario / CLI / MCP Audit

Lifecycle: snapshot (2026-07-06) — dated audit record, never edited.

Date: 2026-07-06  
Base commit: `4eac49e2aebccfa52283c29cd7f22aabba3120f5`  
Scope: suborder H, scenario schema, built-in presets, CLI, gallery, and MCP product paths.

## Summary

The scenario layer is a strict JSON contract (`deny_unknown_fields`) shared by CLI
and MCP. Hard physical rejections are enforced by the same build path used for
`lbm validate` and MCP `validate_scenario`: `lbm_scenario::build_check`, which
dispatches to the 2D compat builder or 3D native `GlobalSpec::validate`.
`validate` is advisory only and must not be treated as the hard gate.

No new physical model, closure, constant, or acceptance band was added in this
audit. The only code change adds scenario-boundary rejection for invalid
`inletProfile` requests so product paths return a structured build error instead
of reaching a core assertion.

## Built-In Presets

| preset | physical claim | scenario parameters | Ma / tau / rho / velocity constraints | core validation coverage | artifacts |
|---|---|---|---|---|---|
| `cavity` | Lid-driven cavity smoke/demo with steady-state detection. Not a full Ghia validation run. | 128x128, TRT, `nu=0.02`, top moving wall `u=[0.1,0]`, closed walls, 20k steps or steady. | `tau=0.56`; wall speed `0.1 < MAX_SPEED=0.3`; single density near 1. | T7 validation exists in core tests for Ghia/reference behavior; preset is smaller/product smoke only. | `speed_<step>.png`, `manifest.json`. |
| `cylinder-karman` | Cylinder wake and force-probe demo; not the full ignored T8 benchmark. | 440x164, cylinder `r=20`, parabolic inlet `umax=0.15`, pressure outlet, TRT, `nu=0.04`, 40k steps. | `tau=0.62`; inlet max `0.15 < 0.3`; nominal mean Re approx 100 on D=40 if using mean `2/3 umax`. | T8/T9 core validation covers cylinder/outflow behavior with reference bands; this preset is a product run path. | `vorticity_*.png`, `speed_<step>.png`, `force.csv`, `manifest.json`. |
| `two-phase-droplet` | Shan-Chen droplet equilibration demo. | 128x128 periodic, TRT, `nu=1/6`, droplet `rho_liquid=2.0`, `rho_vapor=0.15`, `G=-5`. | `tau=1.0`; no imposed velocity; multiphase recommends f64 and uses f64. | T11 core validation covers coexistence, pressure balance, Laplace law, and f32 smoke. | `rho_<step>.png`, `manifest.json`. |
| `droplet-on-wall` | Contact-angle demo for virtual wall density `wallRho=1.0`, approximately 63 degrees from frozen T11c characterization. | 160x100, periodic x, top/bottom bounce-back, droplet on bottom wall, `nu=1/6`, `G=-5`, `wallRho=1.0`, 30k steps. | `tau=1.0`; no imposed velocity; multiphase f64. | T11c core validation covers monotonic wall-rho contact angle and film case. | `rho_<step>.png`, `rho_<step>.vtk`, `manifest.json`. |

## Rejection Coverage

| path | status | evidence / rule |
|---|---|---|
| invalid `nu` / `tau <= 0.5` | Rejected by `build_check`; covered by `build_check_rejects_hard_scenario_physics_errors`. | 2D compat `SimConfig::build`; 3D `GlobalSpec::validate`. |
| velocity too high | Rejected by `build_check`; covered for moving wall and parabolic inlet profile. | Uses `MAX_SPEED=0.3`; NaN-safe `!(speed <= MAX_SPEED)`. |
| illegal orthogonal open faces | Rejected by `build_check`; covered for 2D in the new test and 3D in existing `build3d_runs_and_guards`. | Maintains V1 corner rule: open boundaries must lie on one axis. |
| unsupported GPU physics | Explicit `backend:"gpu"` is rejected when unavailable or unsupported; existing test covers 2D, existing 3D guard covers 3D. | `build_check` returns the same error used by CLI/MCP validation; `auto` may fall back to CPU. |
| invalid source/sink | Scenario schema does not expose `sources` or `sinks`; unknown fields are rejected by serde and covered by `unsupported_source_sink_schema_fields_are_rejected`. | Core native source/sink validation exists under T18, but is not a scenario/MCP product surface yet. |
| invalid core native source/patch definitions | Covered below scenario by `GlobalSpec::validate` and T18 tests. | Scenario always builds `sources: Vec::new()` and `face_patches: Vec::new()`. |

## CLI / Gallery / MCP

- `lbm validate <scenario>` parses JSON, resolves optional units, runs advisory
  warnings, then runs `build_check`; it exits non-zero on hard build errors.
- `lbm run <scenario>` and `lbm presets run <name>` run through the same runner
  and write `manifest.json` plus declared output files.
- `lbm presets list/show/run` uses the same `lbm_scenario::presets()` table as
  the library tests.
- `lbm gallery` runs every preset and embeds produced PNGs into `index.html`.
  This is a useful visual product smoke, but it is long because it runs all
  presets at their configured step counts.
- MCP exposes seven tools: `run_scenario`, `start_run`, `run_status`,
  `list_runs`, `validate_scenario`, `list_presets`, and `get_schema`. The
  validation tool uses the same `build_check` hard gate as CLI validation.

## Findings

1. **S2 validation gap closed in this suborder**: before this audit, a scenario
   with `inletProfile` on a non-inlet edge or with `umax > MAX_SPEED` could reach
   the compat assertion path. The scenario builder now rejects those inputs with
   structured `ConfigError` before constructing the profile.
2. **S2 product-surface limitation**: scenario JSON has no `sources`, `sinks`,
   or `facePatches` surface yet. That is honest and safe because strict serde
   rejects those fields. When these become product features, their schema must
   route directly through `GlobalSpec::validate` and add CLI/MCP rejection tests
   for out-of-face patches, source/solid overlap, source overlap, and
   sink-too-strong errors.
3. **Product-path smoke completed**: `lbm presets run cavity` and `lbm gallery`
   both completed in this audit. The gallery is evidence that the fixed preset
   table, runner, PNG/VTK/CSV outputs, manifest writing, and HTML embedding path
   are wired end-to-end.

## Merge Recommendation

Merge candidate if the release build, workspace tests, and cavity CLI smoke pass.
The change is small, improves hard rejection behavior, and does not alter the
physics model for valid scenarios.

## Behavior Review Record

Run id: `vv-scenario-gallery-smoke-20260706`  
Artifacts:
`out/vv-scenario-cavity-smoke-20260706/speed_20000.png`,
`out/vv-scenario-gallery-smoke-20260706/index.html`,
`out/vv-scenario-gallery-smoke-20260706/cavity/speed_20000.png`,
`out/vv-scenario-gallery-smoke-20260706/cylinder-karman/vorticity_10000.png`,
`out/vv-scenario-gallery-smoke-20260706/cylinder-karman/vorticity_20000.png`,
`out/vv-scenario-gallery-smoke-20260706/cylinder-karman/vorticity_30000.png`,
`out/vv-scenario-gallery-smoke-20260706/cylinder-karman/vorticity_40000.png`,
`out/vv-scenario-gallery-smoke-20260706/cylinder-karman/speed_40000.png`,
`out/vv-scenario-gallery-smoke-20260706/two-phase-droplet/rho_20000.png`,
`out/vv-scenario-gallery-smoke-20260706/droplet-on-wall/rho_30000.png`,
`out/vv-scenario-gallery-smoke-20260706/droplet-on-wall/rho_30000.vtk`.

Pattern: All four built-in presets completed and produced visual field
artifacts plus manifests.  
Mechanism: These are product smoke runs through the scenario/CLI output paths;
their qualitative physical mechanisms are the documented preset mechanisms
(lid-driven recirculation, cylinder wake, Shan-Chen droplet equilibration, and
wall-density contact-angle relaxation).  
Resolved vs closure: cavity and cylinder use resolved single-phase LBM with
validated wall/open-boundary models; droplet cases use the Shan-Chen and
virtual-wall-density closures recorded in `PHYSICS.md` and validated by T11/T11c.  
Artifacts checked: This audit verified artifact creation and manifest
consistency from the run outputs. It did not perform an image-level human
physical-pattern review inside this session; treat the gallery as ready for PM
visual review, not as new physical validation.  
Verdict: VERIFIED product-path smoke, not additional physical validation.  
Routing: none for product path; no new anomaly found.
