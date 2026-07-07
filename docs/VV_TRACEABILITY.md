# V&V Traceability Matrix

Generated for suborder B, V&V-TRACE, on 2026-07-06 from repository inspection only.

This file is a coverage matrix, not a fresh test report. I did not run the physics
validation suites in this session. A row marked `VALIDATED` means the repository
contains an executable or documented validation gate with a physical reference and
acceptance metric. It does not mean that gate passed in this session. Hardware,
manual, or ignored rows keep their execution class explicit.

Status vocabulary:

- `VALIDATED`: direct physical validation evidence exists in tests and/or recorded
  experiment docs.
- `VERIFIED-ONLY`: equivalence, construction, smoke, or implementation-regression
  evidence exists, but it is not an absolute physical validation.
- `SPEC-ONLY`: acceptance is specified but no executable validation evidence was
  found for the full claim.
- `MISSING`: neither executable validation nor a complete spec was found.
- `BENCH-PENDING`: validation depends on heavy, ignored, manual, GPU, MPI, or
  cluster execution not run here.
- `UNSAFE-CLAIM`: a claim would overstate the available evidence or conflicts with
  known findings.

The CSV companion is [vv_traceability.csv](vv_traceability.csv).

## Highest-Risk Gaps

1. T17 stirred-reactor validation is mostly `SPEC-ONLY`; production claims for
   full coupled stirred-reactor physics, aeration, active scalar, scalar reaction,
   and relaxation-mode accuracy would be unsafe.
2. T17 rotating/IBM evidence is not production-grade. `rotating_ibm.rs` contains
   sentinels, but `docs/qa/anomaly-log.md` records legal/default IBM divergence
   findings that require core routing before strong FSI or impeller claims.
3. T14 GPU is still primarily CPU-relative equivalence. That is useful
   verification, but it cannot catch shared CPU/GPU physics errors without the
   absolute GPU physics sentinels called out in `SOLVER_IMPROVEMENT_SPEC.md` D-3.
4. T13-MPI depends on `scripts/test_mpi.sh` and a native MPI toolchain. Standard
   `cargo test --workspace --release` does not cover MPI runtime.
5. T16 FP16 storage has feature-gated ignored tests and PHYSICS band records, but
   `VALIDATION.md` still says "not yet implemented"; the acceptance source is
   internally inconsistent.
6. T18.4 forward-model trend anchors remain example-level and partly stubbed; the
   ad-hoc inventory warns that earlier P1.1 bands measured through example-local
   closures must not be generalized.
7. T18.5 inverse recovery is explicitly a stub and not implementable yet.
8. T15 D3Q27 coverage is partial. D3Q27 has smoke/TGV/partition evidence, but not
   the full D3Q19 T15 physical suite.
9. T9/T9b full outflow pressure-reflection horizons are ignored/heavy; default
   tests cover stability/backflow but not every long-horizon metric.
10. Scenario/CLI reproduction of several validations is incomplete; the QA log
    says T1 analytic initialization is not expressible through scenario JSON.

## Matrix

| ID | Source | Acceptance Metric | Evidence | Command | Run Class | Backend | Precision | Lattice/Dimension | Type | Status |
|---|---|---|---|---|---|---|---|---|---|---|
| T1 | VALIDATION.md T1 | TGV L2rel <= 1.5e-3, order >= 1.7, nu_eff +/-2%, rotation symmetry | `crates/lbm-core/tests/validation_tgv.rs` | `cargo test --release -p lbm-core --test validation_tgv` | default | CPU compat | f64 | D2Q9 / 2D | analytic convergence | VALIDATED |
| T2 | VALIDATION.md T2 | TRT half-way Poiseuille exactness, BGK order, symmetry, rotated profile | `crates/lbm-core/tests/validation_channel.rs` | `cargo test --release -p lbm-core --test validation_channel t2` | default | CPU compat | f64 | D2Q9 / 2D | analytic wall validation | VALIDATED |
| T3 | VALIDATION.md T3 | Couette exactness for BGK/TRT/tau set, alternate wall/orientation, mass drift | `crates/lbm-core/tests/validation_channel.rs` | `cargo test --release -p lbm-core --test validation_channel t3` | default | CPU compat | f64 | D2Q9 / 2D | analytic moving-wall validation | VALIDATED |
| T4 | VALIDATION.md T4 | Inlet/outlet channel bulk flux, center profile, mass drift, four orientations | `crates/lbm-core/tests/validation_open_bc.rs` | `cargo test --release -p lbm-core --test validation_open_bc t4` | default | CPU compat | f64 | D2Q9 / 2D | open-boundary physical validation | VALIDATED |
| T5 | VALIDATION.md T5 | Pressure-pressure flow rate +/-2%, pressure linearity, mirror/sign symmetries | `crates/lbm-core/tests/validation_open_bc.rs` | `cargo test --release -p lbm-core --test validation_open_bc t5` | default | CPU compat | f64 | D2Q9 / 2D | open-boundary physical validation | VALIDATED |
| T6-f64 | VALIDATION.md T6 | Mass conservation, forced momentum growth, feq moments | `crates/lbm-core/tests/validation_conservation.rs` | `cargo test --release -p lbm-core --test validation_conservation` | default | CPU compat | f64 | D2Q9 / 2D | conservation / moment identity | VALIDATED |
| T6-f32 | VALIDATION.md T6 | f32 mass and momentum drift <= 1e-5 class after deviation storage | `crates/lbm-core/tests/validation_conservation.rs` | `cargo test --release -p lbm-core --test validation_conservation t6_f32` | default | CPU compat | f32 | D2Q9 / 2D | precision conservation | VALIDATED |
| T7 | VALIDATION.md T7 | Ghia Re=100/400 default RMS/vortex; Re=1000 ignored; four lid orientations | `crates/lbm-core/tests/validation_cavity.rs` | `cargo test --release -p lbm-core --test validation_cavity`; full: `cargo test --release -- --include-ignored` | default + ignored | CPU compat | f64 | D2Q9 / 2D | literature benchmark | VALIDATED |
| T8-2D-1 | VALIDATION.md T8 | Schaefer-Turek Re=20 Cd/Cl coarse D=20, D=40 convergence ignored | `crates/lbm-core/tests/validation_cylinder.rs` | `cargo test --release -p lbm-core --test validation_cylinder t8_2d1` | default + ignored | CPU compat | f64 | D2Q9 / 2D | literature benchmark / force | VALIDATED |
| T8-2D-2 | VALIDATION.md T8 | Re=100 shedding St/Cd/Cl/periodicity | `crates/lbm-core/tests/validation_cylinder.rs` | `cargo test --release -p lbm-core --test validation_cylinder t8_2d2 -- --ignored` | ignored | CPU compat | f64 | D2Q9 / 2D | unsteady literature benchmark | BENCH-PENDING |
| T9 | VALIDATION.md T9 | Outflow no NaN, backflow <=5%, pressure oscillation ratio <=15x | `crates/lbm-core/tests/validation_outflow.rs` | `cargo test --release -p lbm-core --test validation_outflow t9` | default + ignored | CPU compat | f64 | D2Q9 / 2D | outlet stability / artifact | VALIDATED |
| T9b | VALIDATION.md T9b | ConvectiveOutflow validation and comparison to Outflow | `crates/lbm-core/tests/validation_outflow.rs` | `cargo test --release -p lbm-core --test validation_outflow t9b`; full ignored for pressure ratio | default + ignored | CPU compat | f64 | D2Q9 / 2D | outlet stability / comparison | VALIDATED |
| T10 | VALIDATION.md T10 | Config errors, tau=0.51 stability point, open-edge solid panic | `crates/lbm-core/tests/validation_robustness.rs` | `cargo test --release -p lbm-core --test validation_robustness` | default | CPU compat | f64 | D2Q9 / 2D | construction and stability guard | VALIDATED |
| T11-flat | VALIDATION.md T11 | SC coexistence densities, EOS pressure equilibrium, spurious current, mass | `crates/lbm-core/tests/validation_multiphase.rs` | `cargo test --release -p lbm-core --test validation_multiphase t11_flat` | default | CPU compat | f64 | D2Q9 / 2D | multiphase physical validation | VALIDATED |
| T11-laplace | VALIDATION.md T11 | Laplace law slope/R2 and per-droplet sigma; full sweep ignored | `crates/lbm-core/tests/validation_multiphase.rs` | `cargo test --release -p lbm-core --test validation_multiphase t11_laplace`; full ignored sweep | default + ignored | CPU compat | f64 | D2Q9 / 2D | multiphase physical validation | VALIDATED |
| T11-f32 | VALIDATION.md T11 | f32 flat-interface stability and density bands | `crates/lbm-core/tests/validation_multiphase.rs` | `cargo test --release -p lbm-core --test validation_multiphase t11_f32` | default | CPU compat | f32 | D2Q9 / 2D | multiphase precision smoke | VALIDATED |
| T11b | VALIDATION.md T11b | Contact-angle G_w monotonicity and frozen angles | `crates/lbm-core/tests/validation_contact_angle.rs` | `cargo test --release -p lbm-core --test validation_contact_angle t11b` | default | CPU compat | f64 | D2Q9 / 2D | wall-wetting validation | VALIDATED |
| T11c | VALIDATION.md T11c | Wall-rho contact angle full range and complete-wetting qualitative case | `crates/lbm-core/tests/validation_contact_angle.rs` | `cargo test --release -p lbm-core --test validation_contact_angle t11c`; full angle sweep ignored | default + ignored | CPU compat | f64 | D2Q9 / 2D | wall-wetting validation | VALIDATED |
| T12 | VALIDATION.md T12 | MCMP separation, sigma_AB, RT growth fit; full 256^2 rate ignored | `crates/lbm-core/tests/validation_rt.rs` | `cargo test --release -p lbm-core --test validation_rt`; full ignored RT fit | default + ignored | CPU compat | f64 | D2Q9 / 2D | multiphase instability validation | VALIDATED |
| T13-inprocess | VALIDATION.md T13 | Monolithic vs split fields bit-match; diagnostics tolerance; adversarial seams | `crates/lbm-core/tests/t13_split_invariance.rs`, `t13_adversarial.rs` | `cargo test --release -p lbm-core t13` | default + one ignored long run | CpuScalar / decomposition | f64 | D2Q9, D3Q19, D3Q27 / 2D+3D | partition equivalence | VERIFIED-ONLY |
| T13-MPI | VALIDATION.md T13 | MPI rank counts gathered field match, diagnostics tolerance, mismatch-nu negative | `crates/lbm-core/examples/mpi_t13.rs`, `scripts/test_mpi.sh` | `./scripts/test_mpi.sh` | manual MPI | CPU MPI | f64 | D2Q9, D3Q19 / 2D+3D | distributed equivalence | BENCH-PENDING |
| T14-D2 | VALIDATION.md T14 | CPU vs wgpu f32 relative fields/diagnostics across 7 configs | `crates/lbm-core/tests/t14_backend_equiv.rs`, `t14_adversarial.rs` | `cargo test --release -p lbm-core --features gpu --test t14_backend_equiv` | feature gpu | CpuScalar vs Wgpu | f32 | D2Q9 / 2D | backend equivalence | VERIFIED-ONLY |
| T14-D3 | VALIDATION.md T14 and later tests | CPU vs wgpu D3Q19 relative equivalence | `crates/lbm-core/tests/t14_3d_backend_equiv.rs`, `t14_wale_gpu_equiv.rs` | `cargo test --release -p lbm-core --features gpu --test t14_3d_backend_equiv` | feature gpu | CpuScalar vs Wgpu | f32 | D3Q19 / 3D | backend equivalence | VERIFIED-ONLY |
| T15.1 | VALIDATION.md T15.1 | 3D z-invariant TGV and Zou-He degeneracy to D2Q9; prescribed moments | `crates/lbm-core/tests/t15_3d.rs` | `cargo test --release -p lbm-core --test t15_3d t15_1` | default | CPU native | f64 | D3Q19 / 3D | dimensional degeneracy | VALIDATED |
| T15.2 | VALIDATION.md T15.2 | Rectangular duct series L_inf and flow rate | `crates/lbm-core/tests/t15_3d.rs` | `cargo test --release -p lbm-core --test t15_3d t15_2` | default | CPU native | f64 | D3Q19 / 3D | analytic 3D validation | VALIDATED |
| T15.3 | VALIDATION.md T15.3 | Sphere drag vs Schiller-Naumann, light default and D=24 ignored cases | `crates/lbm-core/tests/t15_3d.rs`, `t15_adversarial.rs` | default light: `cargo test --release -p lbm-core --test t15_3d t15_3_sphere_drag_re20_light`; full ignored | default + ignored | CPU native | f64 | D3Q19 / 3D | literature/correlation validation | VALIDATED |
| T15.4 | VALIDATION.md T15.4 | True 3D TGV decay rate and order; D3Q27 TGV included | `crates/lbm-core/tests/t15_3d.rs` | `cargo test --release -p lbm-core --test t15_3d t15_4` | default | CPU native | f64 | D3Q19, D3Q27 / 3D | analytic convergence | VALIDATED |
| T15.5 | VALIDATION.md T15.5 | 3D cavity Re=1000 default sentinel; N=72 profiles and convergence ignored | `crates/lbm-core/tests/t15_5_cavity3d.rs`, `docs/T15_5_CAVITY3D_REFERENCE.md` | default: `cargo test --release -p lbm-core --test t15_5_cavity3d`; full ignored | default + ignored | CPU native | f64 | D3Q19 / 3D | literature benchmark | VALIDATED |
| T15-f32 | SOLVER D-4 / current tests | f32 3D degeneracy, TGV decay, mass drift | `crates/lbm-core/tests/t15_3d_f32.rs` | `cargo test --release -p lbm-core --test t15_3d_f32` | default | CPU native | f32 | D3Q19 / 3D | precision validation | VALIDATED |
| T16 | VALIDATION.md T16 vs PHYSICS T16 | FP16 storage degradation vs f32 for TGV/cavity; GPU SHADER_F16 required | `crates/lbm-core/tests/t16_fp16_storage.rs`, `docs/PHYSICS.md` FP16 section | `cargo test --release -p lbm-core --features gpu --test t16_fp16_storage -- --ignored` | feature gpu + ignored | Wgpu | f16 storage / f32 compute | D2Q9 / 2D | precision characterization | BENCH-PENDING |
| T17-01 | VALIDATION.md T17 / REQ section 8 | Stirred-tank Np and PIV/LDA velocity lines | no full executable T17 test found | none found | not implemented | intended CPU/GPU/MPI | mixed_safe planned | D3Q19/D3Q27 / 3D | coupled multiphysics benchmark | SPEC-ONLY |
| T17-02 | VALIDATION.md T17 / REQ section 8 | Single bubble, bubble swarm, aerated stirring, epsilon_g/d32/kLa | no full executable T17 test found | none found | not implemented | intended CPU/GPU/MPI | mixed_safe planned | D3Q27 / 3D | gas-liquid validation | SPEC-ONLY |
| T17-03 | VALIDATION.md T17 / REQ section 8 | Stress/shear MMS, curved Couette, rotating cylinder, non-Newtonian, droplet | partial: `wale_les.rs`, `strain_rate.rs`, `rotating_ibm.rs`; DNS skeleton ignored; QA anomaly log records IBM divergence | targeted tests exist, full T17 absent | partial default + ignored skeleton | CPU native | f64/f32 partial | D2Q9/D3Q19 / 2D+3D | stress / LES / IBM validation | UNSAFE-CLAIM |
| T17-04 | VALIDATION.md T17 / REQ section 8 | Scalar and reaction Taylor-Aris, reaction front, kLa | no executable validation found | none found | not implemented | intended CPU/GPU/MPI | mixed_safe planned | 3D | scalar/reaction validation | SPEC-ONLY |
| T17-05 | VALIDATION.md T17 / REQ section 8 | Coupled conservation/regression: mass, momentum, scalar, gas volume, particles, energy-like monitors | partial conservation tests exist outside full coupled system | none for full coupled T17 | not implemented | intended CPU/GPU/MPI | mixed_safe planned | 3D | coupled conservation | SPEC-ONLY |
| T17-06 | VALIDATION.md T17 / REQ section 8 | Well-balanced static stratification and active-scalar degeneracies | partial: `crates/lbm-core/tests/gravity.rs` has `vr_str_06_static_stratification...` | `cargo test --release -p lbm-core --test gravity vr_str_06` | default partial | CPU native | f64/f32 | D2Q9, D3Q19 / 2D+3D | well-balanced gravity | VERIFIED-ONLY |
| T17-07 | VALIDATION.md T17 / REQ section 8 | Initialization independence of quasi-steady statistics | no executable validation found | none found | not implemented | intended CPU/GPU/MPI | mixed_safe planned | 3D | statistical validation | SPEC-ONLY |
| T17-RELAX | VALIDATION.md T17 / REQ section 8 | MRF/PB/one-way/AMR/f32 relative degradation vs fidelity references | no full executable validation found | none found | not implemented | intended CPU/GPU/MPI | mixed profiles | 3D | relaxation equivalence | SPEC-ONLY |
| T18.1 | VALIDATION.md T18.1 / DISPERSED_DEPOSITION CR-1 | Source/sink far-field, mass ledger, jet momentum, validation errors, partition, GPU rejection | `crates/lbm-core/tests/t18_1_interior_source_sink.rs` | default: `cargo test --release -p lbm-core --test t18_1_interior_source_sink`; GPU rejection with `--features gpu` | default + ignored + feature gpu | CPU scalar/partition, GPU rejection | f64 | D3Q19 / 3D | source/sink validation | VALIDATED |
| T18.2 | VALIDATION.md T18.2 / DISPERSED_DEPOSITION CR-2 | Masked face impinging jet mass and wall-jet behavior, validation errors, seam, GPU rejection | `crates/lbm-core/tests/t18_2_masked_face.rs` | `cargo test --release -p lbm-core --test t18_2_masked_face`; GPU rejection with `--features gpu` | default + feature gpu | CPU scalar/partition, GPU rejection | f64 | D3Q19 / 3D | boundary patch validation | VALIDATED |
| T18.3 | VALIDATION.md T18.3 / DISPERSED_DEPOSITION CR-3 | Particle settling vs SN/Stokes, floor crossing, determinism, partition-invariant deposits | `crates/lbm-core/tests/t18_3_particle_settling.rs`, `t18_3_particle_deposition.rs` | `cargo test --release -p lbm-core --test t18_3_particle_settling --test t18_3_particle_deposition` | default | CPU sampler / partitioned field | f64 | particle model over 3D samples | particle/deposition validation | VALIDATED |
| T18.4 | VALIDATION.md T18.4 / DISPERSED_DEPOSITION section 6 | Forward-model monotone trends and gentle CV band | example-level evidence in `DISPERSED_DEPOSITION.md`; no committed CLI example path found under current `crates/lbm-cli/examples` | none found in current tree | example/manual | CLI/example intended | mixed | application-level deposition | forward-model behavior anchors | SPEC-ONLY |
| T18.5 | VALIDATION.md T18.5 / DISPERSED_DEPOSITION section 7 | Inverse recovery of known-good synthetic recipe | explicit stub | none | not implementable yet | not implemented | not implemented | not implemented | inverse validation | SPEC-ONLY |

## Notes

- Standard `cargo test --workspace --release` does not cover `--features gpu`,
  `--features mpi`, or ignored heavy validation.
- `docs/VALIDATION.md` ordering is historical: T15 appears before T13/T14/T16.
  The matrix keeps the source IDs rather than file order.
- `docs/PHYSICS.md` and `docs/DISPERSED_DEPOSITION.md` contain newer T16/T18
  records than the T16 heading in `docs/VALIDATION.md`; reconcile those docs
  before using T16 status in public claims.
- Relative backend equivalence is intentionally marked `VERIFIED-ONLY` unless a
  row has an absolute physical reference on that backend.
