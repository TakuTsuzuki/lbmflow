# Claims-Ledger Cross-Check — V&V Master Plan Lane 3.3

**Date:** 2026-07-07 · **Owner:** V&V worker (lane 3.3) · **Scope:** every
explicit physics / performance / capability claim currently visible in
`docs/PHYSICS.md` and `docs/paper/*.md` mapped to the test that proves it in
`crates/lbm-core/tests/` or `docs/qa/`. Read-only audit; no test or spec
changes made here.

## Sources scanned

- `docs/PHYSICS.md` (719 lines, physics model index + decision log).
- `docs/paper/LBMFlow-whitepaper.md` (246 lines, living technical paper).
- `docs/paper/claims-ledger.md` (34 lines, status-snapshot ledger).
- `docs/paper/benchmark-results.md` (69 lines, raw benchmark tables).
- `docs/qa/anomaly-log.md` (278 lines, live ANOM ledger).
- `docs/qa/VV_MASTER_PLAN.md` §3.3 (references `docs/qa/VV_TRACEABILITY.md`).

**Missing companion file (finding):** `docs/qa/VV_TRACEABILITY.md` is named in
`VV_MASTER_PLAN.md` §3.3 as landed by `cx/vv-trace` (27 VALIDATED / 4
VERIFIED-ONLY / 3 BENCH-PENDING / 8 SPEC-ONLY / 1 UNSAFE-CLAIM). No such file
exists in `docs/` or `docs/qa/`. Either the merge is pending or the file was
lost; this cross-check treats the anchor set as tests-and-ANOM-only until the
matrix lands.

## Verdict taxonomy

- **PROVEN** — claim exact wording is exercised by a named test whose
  assertion band matches the claim's magnitude and whose gate would trip on a
  regression to the wrong side of the claim.
- **PARTIAL** — a test exists but (a) covers a strict subset of the claim
  (fewer cases, coarser resolution, or a smoke variant), (b) its band is
  looser than the paper wording, or (c) claim is measured but not
  regression-pinned in the test.
- **UNPROVEN** — no test in-tree gates the claim (may be BENCH-PENDING,
  SPEC-ONLY, or asserted only in a run-log).
- **STALE-RISK** — claim references a constant / feature / band that a live
  ANOM row has retired, reclassified, or opened for referee. Claim needs a
  conditional qualifier or a follow-up edit.

Tests referenced are inside `crates/lbm-core/tests/` unless noted.

---

## Table — 30 claims covered

| # | Claim (verbatim + doc:line) | Source doc | Proving test file :: fn | Verdict | Notes |
|---|---|---|---|---|---|
| 1 | "**Taylor–Green vortex**: 2nd-order convergence (measured order 1.91), effective viscosity within ±2% of nominal." (whitepaper:78–79) | whitepaper §4 | `validation_tgv.rs :: t1_tgv_trt_accuracy_and_second_order_convergence`, `:: t1_tgv_effective_viscosity_within_two_percent`, `:: t1_tgv_bgk_accuracy_is_comparable` | **PROVEN** | Order and ±2% ν gate both live; also `smoke_tgv.rs` smoke variant. |
| 2 | "**Poiseuille flow (TRT)**: exact to ≤1×10⁻¹⁰ (half-way bounce-back is analytically exact)." (whitepaper:80–81) | whitepaper §4 | `smoke_poiseuille.rs :: trt_magic_is_exact`; `validation_channel.rs :: t2_trt_magic_poiseuille_is_exact_and_symmetric` | **PROVEN** | PHYSICS.md §2 2026-07-04 entry cites the same test. |
| 3 | "**Lid-driven cavity vs Ghia et al. (1982)**: centerline RMS within tolerance at Re = 100/400/1000" (whitepaper:82–83) | whitepaper §4 | `validation_cavity.rs :: t7_lid_driven_cavity_re100_matches_ghia`, `:: t7_lid_driven_cavity_re400_matches_ghia`, `:: t7_lid_driven_cavity_re1000_matches_ghia` (`#[ignore]` heavy) | **PARTIAL** | Re=1000 gate is ignored-by-default; only Re=100/400 run in `--release`. Re=1000 requires `--include-ignored`. |
| 4 | "3D cavity vs Albensoeder–Kuhlmann (2005) at Re = 1000." (whitepaper:83–84) | whitepaper §4 | `t15_5_cavity3d.rs :: t15_5_cavity3d_re1000_default_sanity`, `:: t15_5_cavity3d_re1000_profiles_n72` (`#[ignore]` heavy) | **PARTIAL** | Extremum band frozen at 0.13 (13%) per PHYSICS.md 2026-07-06; profile RMS is the tight gate; N=72 profile run is `#[ignore]`. |
| 5 | "**Cylinder vs Schäfer–Turek**: drag/lift coefficients and Strouhal number in the benchmark bands (2D-1 steady, 2D-2 vortex shedding)." (whitepaper:84–86) | whitepaper §4 | `validation_cylinder.rs :: t8_2d1_d20_cylinder_steady_drag_lift_are_in_reference_bands`, `:: t8_2d2_d40_cylinder_vortex_shedding_matches_reference_bands` (`#[ignore]`) | **PARTIAL** | 2D-2 shedding case is ignored-by-default; only 2D-1 steady runs on the default gate. |
| 6 | "**3D duct**: exact Fourier-series solution to L∞rel 2.3×10⁻⁴" (whitepaper:86–87) | whitepaper §4 | `t15_3d.rs :: t15_2_rectangular_duct_poiseuille_matches_series` | **PROVEN** | Same 2.3e-4 figure cited in PHYSICS.md 2026-07-06 T15.5 entry. |
| 7 | "sphere drag vs Schiller–Naumann within ±10%." (whitepaper:87–88) | whitepaper §4 | `t15_3d.rs :: t15_3_sphere_drag_re20_light` (default), `:: t15_3_sphere_drag_re20` (`#[ignore]`), `:: t15_3_sphere_drag_re100` (`#[ignore]`) | **PARTIAL** | Re=20 heavy and Re=100 both `#[ignore]`; default gate is the light Re=20 variant only. `SCHILLER_NAUMANN_RE_MAX = 800` domain enforced (PHYSICS.md 2026-07-07 entry). |
| 8 | "**Multiphase (Shan–Chen)**: Laplace law R² = 0.9999" (whitepaper:88–89) | whitepaper §4 | `validation_multiphase.rs :: t11_laplace_single_radius_smoke` (default), `:: t11_laplace_four_radius_sweep_is_linear` (`#[ignore]`) | **PARTIAL / STALE-RISK** | Full 4-radius R² sweep is ignored-by-default. Absolute σ value is under three-way referee per ANOM-P4-017 (Taylor-Culick prefactor 0.49×), ANOM-P4-014 CLOSED (Jurin 1.54×), ANOM-P4-018 (near-wall SC artifact), ANOM-P4-019 (contact-line immobility). Claim conditioning: R² is a linearity metric, not an absolute-σ claim, so R² itself survives; any *σ value* claim needs the SC-pressure-tensor lane-1.7 audit gate. |
| 9 | "contact-angle full range" (whitepaper:89) | whitepaper §4 | `validation_contact_angle.rs :: t11c_virtual_wall_density_contact_angles_are_monotone_and_frozen` (`#[ignore]`), `:: t11c_virtual_wall_density_1p6_completely_wets_wall`; `:: t11b_wall_adhesion_contact_angles_are_monotone_and_frozen` | **PARTIAL** | Full 0-180° sweep `#[ignore]`; only monotone-and-frozen 1.6 wetting case and the older `g_wall` T11b run on the default gate. PHYSICS.md 2026-07-05 entry documents the exact θ table. |
| 10 | "Rayleigh–Taylor growth rate γ within 12% of the tension-and-viscosity-corrected reference" (whitepaper:89–92) | whitepaper §4 | `validation_rt.rs :: t12_rt_growth_rate_matches_corrected_dispersion` (`#[ignore]`), `:: t12_rt_default_growth_smoke_reaches_amp_8`, `:: t12_mcmp_sigma_ab_laplace_regression`, `:: t12_mcmp_separation_threshold_smoke` | **PARTIAL** | The 12% γ gate is ignored-by-default (heavy 256²×12k). Default gate is smoke `reaches_amp_8`. |
| 11 | "Split a domain across four sub-blocks … **max\|Δ\| = 0.0** — bit-for-bit identical." (whitepaper:66–68; Appendix A) | whitepaper §4 | `t13_split_invariance.rs` (10 test fns incl. `t13_tgv_periodic_split_invariant`, `t13_cavity_split_invariant`, `t13_channel_profile_outflow_split_invariant`, `t13_cylinder_probe_split_invariant`, `t13_tgv3d_2x2x2_split_invariant`, `t13_tgv3d_d3q27_2x2x2_split_invariant`, `t13_shan_chen_droplet_native_split_invariant`, `t13_uneven_split_and_deeper_decomp`); adversarial: `t13_adversarial.rs` (7 fns) | **PROVEN** | Adversarial suite covers straddling obstacles, split-line probe, L-shape 3-subdomain, per-cell force droplet on 4-rank corner, uneven 3-1-1 splits. |
| 12 | "Rotate an entire problem 90° and the solution reproduces to **4×10⁻¹⁶**, machine precision." (whitepaper:68–69) | whitepaper §4 | `validation_tgv.rs :: t1_tgv_rotated_initial_field_stays_rotationally_symmetric`; `validation_cavity.rs :: t7_re100_cavity_is_exact_under_four_lid_orientations`; `validation_channel.rs :: t2_poiseuille_rotated_90deg_matches_horizontal_profile`, `:: t3_top_wall_couette_exact_for_bgk_and_trt_all_taus`, `:: t3_bottom_wall_driven_couette_exact`, `:: t3_vertical_couette_exact` | **PROVEN** | PHYSICS.md 2026-07-05 rim-corner-orientation entry pins the 3-4e-16 figure. |
| 13 | "The GPU backend reproduces the CPU trajectory to ≤1×10⁻⁵ relative across six 2D scenario classes; a one-ulp control test pins the residual to rounding order" (whitepaper:96–98) | whitepaper §4 | `t14_backend_equiv.rs :: t14_tgv_periodic`, `:: t14_lid_cavity`, `:: t14_channel_inlet_profile_outflow`, `:: t14_cylinder_probe`, `:: t14_cell_force_bgk`, `:: t14_convective_outflow`, `:: t14_pressure_driven_channel`, `:: t14_pressure_bc_ulp_sensitivity_control`; 3D: `t14_3d_backend_equiv.rs` (3 fns); `t14_backend_equiv.rs :: t14_gravity_body_force_device_resident` (`#[ignore]`) | **PROVEN** | `FIELD_TOL = 1e-5` in file. Absolute checks: `gpu_absolute.rs` covers analytic-vs-GPU T1/T2/T7/T15.4 (Ghia case `#[ignore]`). Ulp control test explicitly named. |
| 14 | "**Precision transparency.** Deviation-storage f32 is validation-grade: uniform-force momentum error 1.3×10⁻³ → 2.8×10⁻⁷; Taylor–Green f32 error 7.1×10⁻⁴ vs f64 7.0×10⁻⁴." (whitepaper:157–158) | whitepaper §5 | `backend_simd_equiv.rs :: tgv_2d_f32` etc.; `validation_multiphase.rs :: t11_f32_flat_interface_stays_finite_and_close_to_f64`; PHYSICS.md 2026-07-05 "Deviation storage" entry names 1.34e-3 → 2.8e-7 (4800×) | **PARTIAL** | No single named test asserts *both* specific numbers as regression pins; PHYSICS.md carries them; `t15_3d_f32.rs` and `smoke_tgv.rs` cover TGV f32 pathways. |
| 15 | "τ = 3ν + ½ (cs² = ⅓)." (whitepaper:52–53; PHYSICS.md:18) | whitepaper §3 | Enforced structurally in `crates/lbm-core/src/*.rs`; validated indirectly by every quantitative viscosity test (T1, T2, T15). | **PROVEN** | Structural invariant per CLAUDE.md. |
| 16 | "TRT with the magic parameter Λ=3/16 is the default: it makes the wall position exact and is as fast as BGK." (whitepaper:46–47) | whitepaper §3 | `smoke_poiseuille.rs :: trt_magic_is_exact`; PHYSICS.md 2026-07-04 entry pins `L∞_rel < 1e-10` at H=8, τ=0.8. | **PROVEN** | |
| 17 | "Body forces use second-order Guo forcing" + "Velocity moments include the Guo forcing F/2 correction" (whitepaper:48; CLAUDE.md invariant; PHYSICS.md §1) | whitepaper §3 | Multi-file: `body_force_field.rs`, `gravity.rs :: gravity_channel_is_bit_identical_to_raw_rho_g_force_field`, all TGV/Poiseuille tests; ANOM-P2-001 stub `mf_interim.rs:265` and `accuracy_audit.rs:471` (open) | **PARTIAL / STALE-RISK** | Steady-state Guo is proven. **ANOM-P2-001 OPEN**: uniform-force vs per-cell force-field transient impulse mismatch (Guo half-force delivery differs — 1/(2 τ_minus)·F impulse deficit measured). Whitepaper's "second-order Guo" claim is steady-state true but transient wrong for the per-cell-field path. |
| 18 | "solid walls are half-way bounce-back" (whitepaper:48; PHYSICS.md:39) | whitepaper §3 | `validation_channel.rs` Poiseuille exactness (T2) is the direct proof of half-way placement; `t15_2_rectangular_duct_poiseuille_matches_series` L∞rel 2.3e-4. | **PROVEN** | |
| 19 | "open boundaries use a normal-parameterized Zou–He implementation covering D2Q9 and D3Q19 velocity inlets and pressure outlets on open faces. D3Q27 open faces are not supported today." (whitepaper:48–52) | whitepaper §3 | `t15_3d.rs :: t15_1b_zou_he_channel_degenerates_to_d2q9`, `:: t15_1c_zou_he_3d_enforces_prescribed_moments`; `d3q27_open_bc.rs :: d3q27_open_faces_enforce_velocity_and_pressure_moments_all_orientations`, `:: d3q27_open_duct_matches_series_shape_and_d3q19`, `:: t13_d3q27_open_duct_split_invariant_with_bc_seams`, `:: d3q27_unimplemented_open_face_kinds_are_rejected` | **STALE (claim wording)** | PHYSICS.md 2026-07-07 entry landed D3Q27 velocity inlet + pressure outlet. Whitepaper §3 wording "D3Q27 open faces are not supported today" needs update — only D3Q27 `Outflow` and `Convective` remain unsupported. |
| 20 | "**Bit-exact partition invariance and backend equivalence.**" (whitepaper:94–96) | whitepaper §4 | `backend_simd_equiv.rs` (17 named fns incl. `split_2x2_cpusimd_matches_monolithic_cpuscalar`, `two_pass_streaming_cpusimd_matches_cpuscalar`, `split_2x2x1_duct_3d_cpusimd_matches_cpuscalar`); `t13_split_invariance.rs` bit-identical (0.0) split. | **PROVEN** | CpuSimd-vs-CpuScalar and split-vs-monolithic both explicitly bit-identical. |
| 21 | "2D D2Q9 1024²: **7,073** MLUPS" (whitepaper §5 GPU table; benchmark-results.md §2) | whitepaper §5 | Not a test-gated claim. Reproduction command: `cargo run --release --features gpu -p lbm-core --example bench_gpu` (Appendix A). Raw CSVs in `~/projects/cfd-bench/`. | **UNPROVEN (in-tree)** | Benchmark run, not regression-pinned; falls under claims-ledger "GREEN measured tonight or earlier" row. Machine-tied to Apple M5 Max. |
| 22 | "3D D3Q19 192³: **2,791–2,813** MLUPS (quiet-window A/B/A, target ≥1,500)" (whitepaper §5) | claims-ledger row 1 | Not a test; benchmark run. `example bench_gpu3d`. | **PARTIAL** | Claims-ledger explicitly marks as GREEN measured; no regression pin. The earlier 1353 MLUPS was a loaded-window artifact per claims-ledger note. |
| 23 | "**FP16 storage** … capacity/throughput mode: ~2× MLUPS at 2048², D3Q19 f16 >5 GLUPS; TGV transient 1.401×10⁻¹ (band 2×10⁻¹), cavity steady 2.579×10⁻³ (band 5×10⁻³)" (whitepaper §5; PHYSICS.md 2026-07-06 entry) | whitepaper §5 | `t16_fp16_storage.rs :: t16_tgv2d_f16_storage_degradation_vs_f32_gpu` (`#[ignore]`), `:: t16_cavity2d_f16_storage_degradation_vs_f32_gpu` (`#[ignore]`) | **PARTIAL** | Bands frozen, tests exist, but both are ignored-by-default (need SHADER_F16 GPU adapter). |
| 24 | "2D CPU MLUPS: 1,480 (2048²/18T/f32)" (whitepaper §5 CPU table; benchmark-results.md §1) | benchmark-results.md | Bench example: `bench_backends -- simd f32 2048 18 400`; no regression-pinned test. | **UNPROVEN (in-tree)** | Benchmark, machine-tied. |
| 25 | "3D CPU MLUPS: 302 (192³/18T/f32); 267 (128³)" (whitepaper §5) | benchmark-results.md | Bench example. | **UNPROVEN (in-tree)** | Same as #24. |
| 26 | "OpenLB head-to-head 128³ 18-thread: LBMFlow 266.6 vs OpenLB 298.8 MLUPS" (whitepaper §5) | benchmark-results.md | External. `~/projects/cfd-bench/run_openlb_sweep.sh`. | **UNPROVEN (in-tree)** | Machine-and-window-tied; whitepaper discloses the SIMD-vs-scalar fairness caveat. |
| 27 | "seven MCP tools including an asynchronous job lifecycle" / "JSON scenario schema" (whitepaper §1, §6) | whitepaper §6 | Structural: `crates/lbm-cli` MCP server. No lane-3.3 physics claim. | **PROVEN (structural)** | Verified by `lbm mcp` CLI feature listing per CLAUDE.md; not a physics gate. |
| 28 | "**Multi-node scaling**: 64-rank weak ≥80% measured" (claims-ledger row 4) | claims-ledger | `crates/lbm-core/examples/bench_mpi.rs`, `examples/mpi_t13.rs`, `scripts/bench_mpi.sh`, `scripts/qa/mpi_local_preflight.sh` | **UNPROVEN** | Claims-ledger explicitly **RED**. Whitepaper §7 discloses as "not claimed as measured today"; MPI is functional coverage / n≤4 weak scaling only. Consistent — no stale-risk. |
| 29 | "**Full-physics stirred workload**" (claims-ledger row 5) | claims-ledger | Subsystems: T17 W0/W-ROT/W-GRAV/W-LES landed per VALIDATION.md; W-VOF, W-BCTOP, W-SCAL, W-REACT, W-BUB, coupled W-COUP/W-IO pending. | **UNPROVEN** | Claims-ledger **RED**. Whitepaper §7 discloses. Consistent. |
| 30 | "T17/VR-STR-03 Reτ=178 turbulent channel vs MKM DNS: mean U+ L2rel ≤ 0.30 (measured 0.2328), stress balance ≤ 0.10 (measured 0.0535), turbulence guard -<u'v'>+ > 0.4 (measured 0.729)" (PHYSICS.md:665–706) | PHYSICS.md | `validation_channel_dns.rs :: channel_re_tau_180_wale_vs_mkm_dns` (`#[ignore]`, T17/VR-STR-03 heavy CPU ~40-60 min), `:: channel_re_tau_180_wale_vs_mkm_dns_gpu` (`#[ignore]`), `:: wale_channel_laminar_harness_smoke` (default) | **PARTIAL** | Bands frozen and gate exists; heavy gate ignored-by-default. PHYSICS.md discloses "wall-UNRESOLVED LES coarse-LES grade" — behavior-anchor honestly attached. |

## Additional PHYSICS.md-only claims (spot-check for coverage)

| # | Claim | Proving test | Verdict | Notes |
|---|---|---|---|---|
| 31 | "Rim-corner orientation: 'faster wall wins' … L∞ ~ 3–4e-16" (PHYSICS.md 2026-07-05) | `validation_cavity.rs :: t7_re100_cavity_is_exact_under_four_lid_orientations` | **PROVEN** | Explicit 4-lid orientation cross-check. |
| 32 | "Convective outflow needs mass-consistency pinning … healthy over 34k steps" (PHYSICS.md 2026-07-05) | `validation_outflow.rs`, `backend_simd_equiv.rs :: convective_outflow_2d_f64/f32`, `t14_backend_equiv.rs :: t14_convective_outflow`; PHYSICS.md pins the mechanism. | **PROVEN** | |
| 33 | "Well-balanced gravity composition point + backend-side gravity" (PHYSICS.md 2026-07-06, 2026-07-07) | `gravity.rs :: closed_box_gravity_forms_stable_hydrostatic_stratification`, `:: gravity_channel_is_bit_identical_to_raw_rho_g_force_field`, `:: vr_str_06_static_stratification_quiescent_all_lattices_and_precisions`, `:: shan_chen_gravity_composes_with_force_overwrite_and_creates_buoyancy` | **PROVEN** | |
| 34 | "WALE default … `S^d:S^d ≡ 0` in pure shear, ν_t = 0 in laminar" (PHYSICS.md §1) | `wale_les.rs :: wale_null_for_steady_couette_and_poiseuille`, `:: les_on_does_not_change_laminar_duct_after_null_update`, `:: constant_omega_field_is_bitwise_identical_to_field_off_scalar/simd`, `:: wale_tgv64_nu_eff_characterization` (`#[ignore]`), `:: wale_unset_clipping_matches_raw_wale_bitwise_on_sheared_field`, `:: wale_tau_eff_clipping_diagnostics_match_reference`, `:: wale_tau_eff_clipping_count_is_monotone_with_bound` | **PROVEN** | |
| 35 | "Rotating IBM: direct-forcing iterations reduce marker slip and conserve momentum; marker straddling partition seam is bit-identical" (PHYSICS.md §1 rotating bodies) | `rotating_ibm.rs :: direct_forcing_iterations_reduce_marker_slip_and_conserve_momentum`, `:: marker_straddling_partition_seam_is_bit_identical`, `:: ibm_moving_wall_couette_matches_native_moving_wall_profile`, `:: taylor_couette_marker_profile_matches_analytic_band` | **PROVEN / STALE-RISK** | Direct forcing at the *tested* stable point works; **ANOM-P4-001 OPEN**: the *default-config* time-stepped direct-forcing IBM diverges (n_markers ∈ {63, 160}, relaxation=1.0). Any whitepaper claim that IBM is broadly production-grade under default configuration must be conditioned on the P4-001 fix. |
| 36 | "T18.3 particle deposition Stokes/Schiller-Naumann settling" (PHYSICS.md §1, 2026-07-07 entry) | `t18_3_particle_settling.rs :: t18_3_single_particle_terminal_velocity_matches_schiller_naumann`, `:: t18_3_particle_step_is_bit_deterministic`; `t18_3_particle_deposition.rs`; `crates/lbm-core/src/particles.rs :: tests::schiller_naumann_in_domain_matches_formula_and_is_monotone`, `:: schiller_naumann_out_of_domain_reports_particle_index_and_re` | **PROVEN** | Explicit `SCHILLER_NAUMANN_RE_MAX = 800` domain enforced; ANOM-Dropped confirms the silent clip was removed. |
| 37 | "Rotor penalization default χ=1, ramp=200; torque = reaction torque" (PHYSICS.md 2026-07-06 rotor entry) | `rotating_ibm.rs`; compat rotor tests | **PROVEN / STALE-RISK** | Nominal path proven; **ANOM-P4-010 OPEN**: compat volume-penalization diverges for solid disc at Re=0.09 (thin-blade / 2-blade edge case). Same family as P4-001. |
| 38 | "Cumulant / CentralMoment: discrete second-order Hermite equilibrium as target" (PHYSICS.md 2026-07-06/07 entries) | `cumulant_holdout.rs`, `cumulant_acceptance.rs`, `accuracy_audit_cumulant.rs :: h²-intercept e2 canary` | **PROVEN** | ANOM-P4-008 CLOSED verifies the offset removal (a = 2.17e-5, |a| ≤ 4e-3 light band). |

---

## Stale-risk rows (require conditional wording, follow-up, or edit)

- **STALE (should edit).** Whitepaper §3 line "D3Q27 open faces are not
  supported today" is now false for velocity inlet + pressure outlet; only
  D3Q27 Outflow/Convective remain unsupported. PHYSICS.md 2026-07-07 entry
  landed the closure with `d3q27_open_bc.rs` and `d3q27_open_metamorphic.rs`
  green. Whitepaper needs a one-line correction.
- **STALE-RISK (needs conditioning).** All whitepaper §4 Shan-Chen absolute-σ
  claims should be qualified. The Laplace R² = 0.9999 claim is a *linearity*
  claim and survives, but the SC pressure tensor delivers three different σ
  values on curved menisci vs flat vs retracting rims (three-way referee open
  under ANOM-P4-017 / ANOM-P4-014-closed / ANOM-P4-018 / ANOM-P4-019). Any
  future whitepaper claim of an absolute σ or a specific-scenario dynamic SC
  behavior needs to await the lane-1.7 SC-pressure-tensor audit or be
  explicitly narrowed. **Physics.md does not yet flag this** — recommend
  adding a §2 "SC dynamic validity" entry citing the referee.
- **STALE-RISK (behavior anchor).** Whitepaper §4 "second-order Guo forcing"
  is steady-state correct but hides an OPEN transient defect
  (ANOM-P2-001: 4/7·F impulse deficit for per-cell force-field path;
  test-side pin lives in `mf_interim.rs:265` and `accuracy_audit.rs:471`
  with `#[ignore]`). No paper claim contradicts this today, but any
  future *transient-impulse* claim would need the R2-C fix first.
- **STALE-RISK (rotational bodies).** ANOM-P4-001 (time-stepped IBM
  divergence at default relaxation=1.0) and ANOM-P4-010 (compat
  penalization solid-disc divergence at Re=0.09) are both OPEN with
  test-side red gates on `cx/audit-ibm` B1-B8 and `cx/audit-rotor` F1-F5.
  Whitepaper §3 currently states IBM as a landed capability ("Prescribed
  rigid rotating-boundary IBM reports force, torque, slip, and
  momentum-spreading diagnostics; it is not a general structural FSI
  solver") — this narrow scope disclosure survives, but any claim of a
  broader default-config production-grade IBM would be stale.
- **STALE-RISK (banned resolution-point calibration).** ANOM-P4-008 CLOSED
  (2026-07-07): the D3Q19 central-moment `+0.0025` offset was a banned
  calibration and is removed from CPU/SIMD/WGSL. PHYSICS.md 2026-07-07
  entry documents the removal. Any lingering paper draft or doc that
  quotes cumulant behavior at the `+0.0025`-active operating point would
  be stale — **spot-check clean**: whitepaper does not cite the offset;
  PHYSICS.md 2026-07-06 cumulant entry has been superseded in place by
  the 2026-07-07 ANOM-P4-008 entry (both cited above). The remaining
  `-0.16|u|²` term is under ablation flag and validated as (B) per
  ANOM-P4-008 FULLY CLOSED addendum. **No live stale text**.
- **STALE (traceability infrastructure).** `docs/qa/VV_TRACEABILITY.md` is
  named as landed in VV_MASTER_PLAN §3.3 but does not exist in the tree.
  Either merge is pending (task #12) or the file was lost. This
  cross-check is the interim substitute.

## Traceability-only rows lacking a paper anchor (informational)

Not applicable this pass — `docs/qa/VV_TRACEABILITY.md` is absent, so no
traceability-only rows exist to reconcile against paper anchors. When the
master matrix lands, this section should be re-scanned for rows whose
`Claim` column does not tie to `PHYSICS.md` or the whitepaper text.

## PHYSICS.md rows without an obvious traceability anchor (interim, until VV_TRACEABILITY lands)

- PHYSICS.md 2026-07-06 "dispersed seeding closure removal" —
  behavior-review verdict UNKNOWN / CAPABILITY GAP, not a positive claim.
  Not paper-anchored; correctly not in the paper. No stale risk.
- PHYSICS.md 2026-07-07 "WALE `tau_eff` clipping" — diagnosed
  numerical-stability guard, explicitly *not* a validation-band knob;
  the follow-up entry mandates disclosure of `clipped_fraction` and
  `max_nu_t_before_clipping` on any validation claim made with clipping
  active. Test-anchored (`wale_les.rs`). No paper claim to align.
- PHYSICS.md 2026-07-06 "Patched Closed face is a zero-velocity Zou-He
  lid" — internal frozen semantics, no paper claim; anchored via
  `face_patch_smoke.rs`, T18.2 impinging-jet.

## Summary

- **30 paper/PHYSICS.md claims mapped** (rows 1–30) + **8 PHYSICS.md-only
  spot-checks** (rows 31–38) = 38 rows.
- Verdict distribution: **PROVEN 15**, **PROVEN (structural) 1**,
  **PROVEN / STALE-RISK 3**, **PARTIAL 10**, **PARTIAL / STALE-RISK 1**,
  **STALE 1 (D3Q27 open-face wording)**, **UNPROVEN (in-tree) 4**
  (all four are benchmark MLUPS numbers with disclosed reproduction
  commands, machine-tied), **UNPROVEN 2** (ME-3, ME-4 — both openly RED
  in claims-ledger, correctly disclosed in whitepaper §7).
- **Zero fabricated claims found**: every paper claim has either a
  test-backed gate or an explicit "not claimed" §7 disclosure.
- **Follow-up edits recommended**:
  1. Whitepaper §3 line "D3Q27 open faces are not supported today" →
     narrow to "D3Q27 open-face `Outflow` and `Convective` are not
     supported today; velocity inlet and pressure outlet landed
     2026-07-07 with tensor-product tangent-deficit distribution".
  2. Whitepaper §4 Shan-Chen block → add a footnote linking to the
     ANOM-P4-014/017/018/019 SC dynamic-limits referee; strengthen the
     existing "This is the validated Shan-Chen scope" hedge with an
     explicit "static Laplace σ = 1; dynamic mechanical σ is under
     characterization (three-way referee)".
  3. PHYSICS.md §2 → new entry summarizing the SC three-way referee
     status so the whitepaper footnote has a canonical PHYSICS-side
     anchor.
  4. Land `docs/qa/VV_TRACEABILITY.md` (VV_MASTER_PLAN §3.3 task #12) or
     mark the master-plan row as SUPERSEDED-BY-THIS-DOC.

**Ledger health, one line:** the paper is measurement-honest — every
non-RED claim has a test or a reproduction command, every RED claim is
disclosed in §7, and the four live stale-risks are catalogued above with
either an edit path or a follow-up gate.
