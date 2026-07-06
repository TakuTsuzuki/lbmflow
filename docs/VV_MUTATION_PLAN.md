# V&V Mutation Plan

Scope: suborder C only. The purpose is to test whether existing validation fails when specific physics invariants are broken. Mutations are temporary and must never be committed as engine code.

## Safe Temporary Protocol

Runner: `scripts/vv/mutation_probe.sh`.

Protocol:
- Refuse to mutate any dirty engine target file.
- Copy every target file into a temporary backup directory under `/tmp`.
- Apply one exact string mutation.
- Run the configured sentinel command.
- Interpret a nonzero test result as `KILLED`; interpret exit 0 as `SURVIVED`.
- Restore all target files in a trap on normal exit, failure, interrupt, or termination.
- Verify `git diff --quiet -- <target-files>` after restore.

This keeps the committed branch limited to the mutation plan, runner, and findings.

## Mutation Inventory

| ID | Mutation | Primary target | Sentinel |
|---|---|---|---|
| `guo-f2-velocity-removed` | Remove the Guo `F/2` correction from physical velocity / momentum diagnostics | `backend.rs`, `kernels.rs` | `validation_conservation::t6_periodic_uniform_force_adds_exact_momentum` |
| `forcing-sign-flipped` | Flip the Guo source term sign | `kernels.rs` | `validation_conservation::t6_periodic_uniform_force_adds_exact_momentum` |
| `trt-relaxation-swapped` | Use the symmetric relaxation for the antisymmetric TRT branch | `kernels.rs` | `validation_channel::t2_trt_magic_poiseuille_is_exact_and_symmetric` |
| `d2q9-opposite-broken` | Break a D2Q9 opposite-direction entry | `lattice.rs` | `lattice::tests::d2q9_invariants` |
| `d3q19-opposite-broken` | Break a D3Q19 opposite-direction entry | `lattice.rs` | `lattice::tests::d3q19_invariants` |
| `d3q27-face-unknown-broken` | Break a D3Q27 face-unknown set | `lattice.rs` | `lattice::tests::d3q27_invariants` |
| `halfway-wall-shifted` | Read the reflected population from the solid cell, effectively shifting the half-way wall treatment | `kernels.rs` | `validation_channel::t2_trt_magic_poiseuille_is_exact_and_symmetric` |
| `moving-wall-sign-flipped` | Flip moving-wall momentum injection | `kernels.rs` | `validation_channel::t3_top_wall_couette_exact_for_bgk_and_trt_all_taus` |
| `zou-he-pressure-normal-sign-flipped` | Flip pressure-Zou-He normal velocity closure | `kernels.rs` | `t15_3d::t15_1c_zou_he_3d_enforces_prescribed_moments` |
| `pressure-outlet-correction-removed` | Remove tangential correction from Zou-He reconstruction | `kernels.rs` | `validation_open_bc::t4_velocity_inlet_pressure_outlet_channel_all_four_orientations` |
| `outflow-stale-slot-broken` | Let streaming overwrite open-face unknown slots instead of preserving stale slots | `kernels.rs` | `stream_contract::cpu_stream_preserves_open_face_unknowns_d2q9` |
| `mpi-halo-x-direction-swapped` | Perturb MPI halo message tagging in the x-exchange family | `dist.rs` | MPI feature unit sentinel if MPI is available |
| `probe-force-physicalization-removed` | Remove `+2w` physicalization from momentum-exchange force probes | `kernels.rs` | `accuracy_audit_probe::a2_steady_poiseuille_wall_friction_balance` |
| `shan-chen-force-sign-flipped` | Flip Shan-Chen cohesion force sign | `compat/multiphase.rs` | `validation_multiphase::t11_laplace_single_radius_smoke` |
| `contact-angle-wall-term-sign-flipped` | Flip legacy Shan-Chen wall-adhesion sign | `compat/multiphase.rs` | `validation_contact_angle::t11b_wall_adhesion_contact_angles_are_monotone_and_frozen` |
| `f32-deviation-storage-disabled` | Remove the deviation-storage `+1` density reconstruction term | `kernels.rs` | `validation_conservation::t6_f32_mass_and_momentum_hold_with_tightened_tolerance` |

## Run Set

Executed lightweight probes for this suborder:
- `moving-wall-sign-flipped`: killed by `validation_channel::t3_top_wall_couette_exact_for_bgk_and_trt_all_taus`.
- `zou-he-pressure-normal-sign-flipped`: survived `t15_3d::t15_1c_zou_he_3d_enforces_prescribed_moments`.
- `outflow-stale-slot-broken`: killed by `stream_contract::cpu_stream_preserves_open_face_unknowns_d2q9`.

The run set intentionally covers wall physics, open-boundary pressure closure, and open-face streaming memory. Longer multiphase, contact-angle, and MPI probes remain listed but were not executed in this lightweight pass.
