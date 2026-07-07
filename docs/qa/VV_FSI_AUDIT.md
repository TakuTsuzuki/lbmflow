# V&V FSI Audit

Lifecycle: snapshot (2026-07-06) — dated audit record, never edited.

Date: 2026-07-06  
Suborder: E. V&V-FSI  
Scope: FSI / IBM / rotating-boundary / particle-coupling claims only. Fluid-only validation is not counted as FSI validation.

## Executive Status

Current FSI-grade capability is limited and must be described precisely:

- Rotating 2D compat impeller volume penalization is **VERIFIED-ONLY** against algebraic no-overshoot, torque sign, solid-body tracking, and stability-envelope characterization. It is not validated as a stirred-tank power-number model.
- Native marker-based rotating IBM is **VERIFIED-ONLY**. It has marker-slip, momentum-spreading conservation, partition-invariance, and coarse Taylor-Couette/Couette sentinels, but the current profile bands are characterization-grade and not a high-order curved-wall validation.
- One-way Lagrangian particles and T18 deposition are **VERIFIED-ONLY** for Schiller-Naumann settling, deterministic interpolation, floor-crossing deposition, and partition-invariant deposition maps. They are not two-way or resolved-particle FSI.
- T17 stirred-reactor, `Np`, `N_Q`, gas-liquid coupling, resolved high-density-ratio phase field, two-way particles, Turek-Hron, and beam FSI are **SPEC-ONLY**.

No current claim should say that the solver has validated full FSI. The correct claim is: selected rotating-boundary and one-way particle submodels have verification sentinels and characterization records; full coupled FSI validation remains future work.

## Inventory And Classification

| Area | Implementation / evidence | Classification | Finding route |
|---|---|---:|---|
| 2D compat rotor volume penalization | `crates/lbm-core/src/compat/rotor.rs`; unit tests for algebraic no-overshoot, ramped angle, torque sign, solid-body tracking; PHYSICS.md rotor stability envelope | VERIFIED-ONLY | Future validation gate for stirred tank `Np`, `N_Q`, discharge profile |
| Native rotating IBM | `crates/lbm-core/src/rotating_ibm.rs`, `Solver::apply_rotating_ibm`; `crates/lbm-core/tests/rotating_ibm.rs` covers slip reduction, torque sign symmetry, partition seam, coarse Couette/Taylor-Couette characterization | VERIFIED-ONLY | Core-model limitation: current curved-wall profile bands are loose characterization |
| Moving wall / rotating cylinder torque sanity | Compat rotor documents reaction torque on the body and has a sign test. Native IBM has exact torque sign symmetry under omega reversal, but the positive-omega sign observed in the sentinel should be reviewed against the reaction-torque prose before any `Np` validation claim | VERIFIED-ONLY | Core-model/documentation finding: settle native IBM torque sign convention before power-number validation |
| Particles | `crates/lbm-core/src/particles.rs`; one-way Schiller-Naumann drag + buoyancy-reduced gravity, deterministic sampler, solid contact, deposition | VERIFIED-ONLY | Model-validity limitation |
| Near-neutral particles | Current model has no added mass, lift, Basset history, Faxen correction, or resolved particle coupling; scenario validation now emits a warning for particle-model scope | SPEC-ONLY for full finite-size near-neutral FSI | Model-validity limitation |
| Added mass / reaction force | Not implemented for particles; no two-way reaction-force scatter from particles to fluid | SPEC-ONLY | Future validation gate; warn users now |
| T18 localized source/sink | `t18_1_interior_source_sink.rs`; analytic sink far-field, mass ledger, jet momentum, partition invariance, GPU rejection | VALIDATED for CR-1 scope only | Not FSI; relevant as forcing/source infrastructure |
| T18 masked faces | `t18_2_masked_face.rs`; impinging-jet behavior anchor, mass drift, validation errors, partition invariance, GPU rejection | VALIDATED for CR-2 scope only | Not FSI; relevant as boundary infrastructure |
| T18 deposition | `t18_3_particle_settling.rs`, `t18_3_particle_deposition.rs`; settling and deposition sentinels | VERIFIED-ONLY / VALIDATED only for one-way CR-3 pieces listed in T18 | Model-validity limitation |
| Stirred reactor T17 | VALIDATION.md and REQ_STIRRED_REACTOR.md wire VR-STR-01..07, but implementation is incomplete | SPEC-ONLY | Future validation gate |
| `Np`, `N_Q` stirred-tank claims | Definitions exist in REQ; no accepted benchmark run | SPEC-ONLY | Future validation gate |
| Turek-Hron / beam FSI | No implementation or tests found | SPEC-ONLY | Future validation gate |

## Unsafe Claims To Avoid

- "FSI validated" is an **UNSAFE-CLAIM** unless tied to a specific future benchmark such as Turek-Hron or beam FSI.
- "Stirred reactor validated" is an **UNSAFE-CLAIM** before VR-STR-01 power number and survey-line gates pass.
- "Two-way particles" or "resolved particle FSI" is an **UNSAFE-CLAIM** for the current `particles` module.
- "Near-neutral particle fidelity" is an **UNSAFE-CLAIM** unless added mass / lift / history / Faxen or resolved-particle coupling is implemented and validated.

## Sentinel Tests Added

- `rotating_cylinder_torque_flips_with_omega_sign`: asserts IBM torque sign symmetry under omega reversal, with identical slip magnitude as the behavior anchor. The positive-omega sign itself is not promoted to validation until the native IBM torque convention is reviewed.
- `validate_downgrades_one_way_particle_model_claims`: asserts scenario validation warns that particles are one-way Schiller-Naumann only and do not validate added mass, near-neutral finite-size behavior, or full FSI.

No new physical force law, closure constant, or acceptance band was introduced.

## Future Acceptance Gates

| Gate | Required evidence |
|---|---|
| Full FSI benchmark | Turek-Hron FSI2/FSI3 or equivalent beam/cylinder benchmark: displacement amplitude/frequency, drag/lift phase, mesh/time convergence, visual field artifacts |
| Rotating-boundary stirred tank | Rushton or equivalent geometry with `Np = P/(rho N^3 D^5)`, `N_Q`, discharge-line velocity profile, torque time series, and visual vortical/velocity artifacts |
| Native IBM curved wall | Taylor-Couette profile with tight L2/Linf bands after resolution study, torque vs analytic wall shear, marker slip, and wall-adjacent behavior artifact |
| Two-way particles | Reaction-force scatter with global momentum conservation, mass-loading limit, one-way vs two-way relaxation comparison, and particle/fluid visual artifacts |
| Near-neutral particle regime | Added-mass/lift/history/Faxen or resolved-particle method with stated validity domain and benchmark against literature/analytic reference |
| Added-mass red zone | Config validation rejects or warns when a requested particle regime needs added mass but only one-way drag is enabled |

## Behavior-Validity Record

The added tests are unit-level sign/scope sentinels and do not produce standalone spatial-field artifacts. A separate CLI rotor smoke run was used for a visual behavior anchor:

- Artifacts: `target/vv_fsi_artifacts/rotor_smoke/speed_120.png`, `target/vv_fsi_artifacts/rotor_smoke/vorticity_120.png`, `target/vv_fsi_artifacts/rotor_smoke/torque.csv`, `target/vv_fsi_artifacts/rotor_smoke/manifest.json`
- Pattern: speed is localized around the rotating blade region with a decaying induced flow in the closed tank.
- Mechanism: compat rotor volume penalization adds Guo body force inside blade cells, driving tangential motion during spin-up.
- Resolved vs closure: the LBM fluid response and Guo forcing are resolved; the blade-fluid coupling is a penalization model characterized in PHYSICS.md, not a fully validated impeller FSI closure.
- Boundary artifact sweep: closed bounce-back tank walls are present; this smoke run is not used as a wall-validation metric.
- Verdict: CLOSURE-DRIVEN, suitable only as a visual sentinel and not as validation of stirred-reactor `Np` or full FSI.
