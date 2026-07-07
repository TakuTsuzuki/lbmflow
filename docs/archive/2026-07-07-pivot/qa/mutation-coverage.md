# Mutation Coverage Matrix

Lane: V&V master plan 2.1, validation-suite mutation testing.

Date: 2026-07-07.

Runner: `scripts/vv/mutation_probe.sh`.

Mechanism:

- `mutation_catalog` declares the available mutant IDs.
- `set_mutation` maps each ID to a description, target source files, and a sentinel command.
- Before applying a mutant, the runner refuses dirty target files.
- The runner copies target files to a temporary backup directory, applies exact text/regex replacements, runs the sentinel, and restores target files from backup through an `EXIT` trap.
- Sentinel failure means the mutant is killed and the runner exits `0`; sentinel success means the mutant survived and the runner exits `10`.
- After restore, the runner checks `git diff --quiet -- "${TARGETS[@]}"` and exits `99` if any temporary source edit remains.

Scope note: the harness catalog contains additional exploratory mutants. This matrix records the 10 lane-2.1 mutants requested for the current V&V pass: the three baseline mutants from `cx/vv-mutation` plus seven added physics mutants.

## Matrix

| Mutant ID | Mutation | Sentinel command | Result | First failing test / failure signal |
|---|---|---|---|---|
| `moving-wall-sign-flipped` | Flip moving-wall bounce-back momentum-injection sign. | `cargo test --release -p lbm-core --test validation_channel t3_top_wall_couette_exact_for_bgk_and_trt_all_taus -- --exact` | CAUGHT | `t3_top_wall_couette_exact_for_bgk_and_trt_all_taus`; T3 top Couette error `1.875` at `tau=0.6`, BGK. |
| `zou-he-pressure-normal-sign-flipped` | Flip pressure Zou-He normal-velocity closure sign. | `cargo test --release -p lbm-core --test t15_3d t15_1d_zou_he_pressure_faces_drive_from_high_density_to_low_density -- --exact` | CAUGHT | `t15_1d_zou_he_pressure_faces_drive_from_high_density_to_low_density`; FA-MUT-001 high-density XNeg pressure face drove negative `ux`. |
| `outflow-stale-slot-broken` | Let streaming overwrite open-face unknown slots instead of preserving stale slots. | `cargo test --release -p lbm-core --test stream_contract cpu_stream_preserves_open_face_unknowns_d2q9 -- --exact` | CAUGHT | `cpu_stream_preserves_open_face_unknowns_d2q9`; open-face unknown slot `65` was overwritten. |
| `guo-source-sign-flipped` | Flip the Guo source term sign in `kernels.rs::collide_row`. | `cargo test --release -p lbm-core --test validation_conservation t6_periodic_uniform_force_adds_exact_momentum -- --exact` | CAUGHT | `t6_periodic_uniform_force_adds_exact_momentum`; momentum gain had wrong sign/magnitude, axis 0 relative error `1.428571428571385`. |
| `moving-wall-factor-two` | Use factor `2.0` instead of `6.0` in the half-way moving-wall term. | `cargo test --release -p lbm-core --test validation_channel t3_top_wall_couette_exact_for_bgk_and_trt_all_taus -- --exact` | CAUGHT | `t3_top_wall_couette_exact_for_bgk_and_trt_all_taus`; T3 top Couette error `0.625` at `tau=0.6`, BGK. |
| `d2q9-weight-nonopposite-swapped` | Swap two non-opposite D2Q9 quadrature weights. | `cargo test --release -p lbm-core lattice::tests::d2q9_invariants -- --exact` | CAUGHT | `lattice::tests::d2q9_invariants`; opposite-direction weights no longer matched (`1/9` vs `1/36`). |
| `halfway-wall-source-cell-offbyone` | Read the reflected wall-link population from `si` instead of the adjacent fluid cell. | `cargo test --release -p lbm-core --test validation_channel t2_trt_magic_poiseuille_is_exact_and_symmetric -- --exact` | CAUGHT | `t2_trt_magic_poiseuille_is_exact_and_symmetric`; T2 TRT `Linf_rel = 0.15238095238094762`. |
| `probe-corner-links-dropped` | Drop diagonal/corner bounce-back links from momentum-exchange probe accumulation. | `cargo test --release -p lbm-core --test accuracy_audit_probe a2_steady_poiseuille_wall_friction_balance -- --exact` | CAUGHT | `a2_steady_poiseuille_wall_friction_balance`; wall-friction residual `-2.56e-4` versus band `2.56e-14`. |
| `trt-omega-minus-uses-tau` | Compute TRT `omega_minus` using `tau` instead of `lambda_plus = tau - 0.5`. | `cargo test --release -p lbm-core params::tests::omegas_match_v1_derivation -- --exact` | CAUGHT | `params::tests::omegas_match_v1_derivation`; `omega_m` was `1.1978609625668448` instead of `0.27586206896551746`. |
| `equilibrium-u2-coefficient-four` | Change `equilibrium()`'s `c_u^2` coefficient from `4.5` to `4.0`. | `cargo test --release -p lbm-core kernels::tests::collide_feq_matches_equilibrium_bitwise_d2q9 -- --exact` | CAUGHT | `collide_feq_matches_equilibrium_bitwise_d2q9`; bitwise fixed-point check failed at slot `64` (`q=1`, cell `0`). |

## Survivors

No survivors in this 10-mutant matrix. No new sentinel proposal is required from this pass.

## Restore Evidence

Each mutant cycle printed `git status --short crates/lbm-core/src` after the runner exited. The status output was empty after every cycle, confirming the temporary mutations were restored. The final verification also checks the pristine source tree before the full workspace test gate.
