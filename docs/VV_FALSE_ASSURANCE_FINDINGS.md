# V&V False-Assurance Findings

This file records mutation survivors from suborder C. A survivor means the test suite accepted a deliberately wrong physics mutation.

## 2026-07-06 Suborder C

### FA-MUT-001: pressure-Zou-He sign mutation survives a prescribed-moment sentinel

- Mutation: `zou-he-pressure-normal-sign-flipped`.
- Temporary edit: changed the pressure-boundary closure from `un = 1 - closure / rho_bc` to `un = closure / rho_bc - 1` in `crates/lbm-core/src/kernels.rs`.
- Sentinel command: `cargo test --release -p lbm-core --test t15_3d t15_1c_zou_he_3d_enforces_prescribed_moments -- --exact`.
- Evidence: the sentinel exited 0 and the runner reported `RESULT: SURVIVED`.
- Interpretation: `t15_1c_zou_he_3d_enforces_prescribed_moments` checks prescribed velocity on the velocity face and prescribed density on the pressure face, but it does not detect a flipped pressure-face normal velocity sign. The pressure-boundary mutation should be covered by a pressure-driven-flow or pressure-face normal-momentum sentinel, not only density enforcement.
- Resolution: added `t15_3d::t15_1d_zou_he_pressure_faces_drive_from_high_density_to_low_density`, which starts from rest with high density on `XNeg` and low density on `XPos` and asserts both pressure-face normal velocities and the interior mean velocity point high-to-low. The mutation runner now uses that sentinel, and the same mutation is killed.
