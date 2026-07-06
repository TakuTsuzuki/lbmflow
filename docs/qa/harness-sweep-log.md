# Physics Anomaly Sweep — Log

Automated physics-QA over runnable (scenario × config) pairs. Reference values
come ONLY from docs/VALIDATION.md frozen bands or analytic solutions — no
invented bands. Harness: `scripts/qa/` (matrix.py = the runnable matrix,
qa_checks.py = detectors, run_sweep.py = RUN→COLLECT→DETECT driver).

Severity taxonomy (house): **S0** correctness (silently wrong physics /
false-green), **S1** high risk (divergence/leak in a supported config),
**S2** improvement (accuracy below expected order/band, tooling gaps),
**S3** minor. Documented limitations (capability-map reds) are
"expected-limitation", not anomalies.

Entry format: `{id, scenario+config, expected (source), observed (data
excerpt), visualization, severity, disposition-proposal}`.

---

## Pass 1 — 2026-07-06, branch `qa/anomaly-sweep`

Matrix: 18 configs across 8 tracks (conservation / channel-analytic /
open-boundary / cavity / cylinder / 3D / multiphase / robustness), all via
the `lbm run` scenario path, f64: 2D D2Q9 (CpuSimd compat) + 3D D3Q19
(CpuScalar). Machine-readable evidence: [pass1/results-main.json]
(pass1/results-main.json) (18-config sweep) and
[pass1/results-rerun-fixed.json](pass1/results-rerun-fixed.json) (4-config
rerun after the harness fixes below). Field/PNG artifacts under
`out/qa-pass1*/<id>/` (gitignored, reproducible via
`python3 scripts/qa/run_sweep.py --bin target/release/lbm --out out/qa-pass1`).

**Headline: no engine physics anomaly found (no S0/S1).** All frozen-band and
analytic checks pass; all four initial check failures were harness or
matrix-fidelity defects, fixed and re-verified in the rerun. PM was NOT
messaged (S0/S1 only rule).

### Runnable-matrix evidence (per-config, main sweep + fixed rerun)

| id | grid | status | steps | key observed vs band |
|---|---|---|---|---|
| t6-momentum-periodic | 64² | completed | 2000 | du/dstep vs F rel err 8.2e-14 (≤1e-10, T6); ux uniform to 0; mass drift 0 |
| poiseuille-trt-h32 | 16×34 | steady | 26000 | L∞rel 1.35e-11 (≤1e-10 exact, T2); symmetry 1.4e-18 (≤1e-13); mass drift 0 |
| poiseuille-bgk-h8/h16 | 12×10/18 | steady | 2500/7500 | BGK convergence order 2.000 (≥1.7, T2) |
| couette-trt-tau08 | 16×34 | steady | 25500 | L∞rel 1.35e-11 (≤1e-10, T3); mass drift 0 |
| couette-bgk-tau06 | 16×34 | steady | 72000 | L∞rel 5.63e-11 (≤1e-10, T3) |
| channel-zouhe-t4 | 96×34 | steady | 67000 | bulk Q dev 6.4e-7 (≤1e-4); profile L2rel 1.08e-3 (≤2e-3, T4) |
| cavity-re100 | 129² | steady | 188000 | Ghia RMS 0.0053U (≤0.02U, T7); mass drift 2.3e-16 |
| cavity-re400 | 129² | 300k cap | 300000 | Ghia RMS 0.0088U (≤0.02U, typo pt excluded) |
| cavity-re1000 | 129² | 300k cap | 300000 | Ghia RMS 0.0133U (≤0.03U); τ=0.538 advisory present |
| cylinder-re20-t8 | 440×82 | completed | 30000 | Cd 5.851 ∈ [5.2,6.0], Cl 0.0092 ∈ [−0.05,0.08] (T8 2D-1) |
| cylinder-re100-karman-t8 | 880×164 | completed | 150000 | St 0.3066 ∈ [0.28,0.32]; Cd_max 3.447 ∈ [3.0,3.5]; |Cl|max 1.091 ∈ [0.8,1.2]; period var 0.06% ≤2% (T8 2D-2, from rest) |
| outflow-karman-t9 | 880×164 | completed | 60000 | reverse flow ≤0.42% (≤5%, T9); no NaN |
| duct3d-t15 | 12×34×34 | steady | 13500 | series L∞rel 2.31e-4 (≤1e-3; spec measured 2.3e-4); Q dev 0.094% (≤0.5%; spec measured 0.094%); y↔z asym 3.3e-19 |
| cavity3d-re1000-n64 | 64³ | 40k cap | 40000 | mass drift 1.2e-16; z-mirror asym 4.0e-15/U (T15.5 sentinel ~2e-15); τ+grid-Re advisories present |
| droplet-t11 | 128² | completed | 30000 | spurious max|u| 4.3e-3 (≤5e-3, T11); Laplace σ 0.0349 = ref+5.1% (±15%); mass drift 3.9e-12 (≤1e-10); mirror asym 7.1e-12 |
| droplet-on-wall-t11c | 160×100 | completed | 30000 | θ 64.5° vs 63±8° (T11c wallRho=1.0) |
| tau-margin-cavity-t10 | 128² | completed | 10000 | stable, max|u| 0.0459 (T10 frozen point measured 0.046); τ=0.51 advisory present |

### Anomalies / findings

**ANOM-P1-001 — conservation diagnostics are not first-class collection
surface** — S2 (tooling), disposition-proposal: engine/CLI improvement.
- Scenario+config: any; bitten twice in pass 1 (`poiseuille-trt-h32`,
  `duct3d-t15` — steady-stop at 26k/13.5k steps landed before the first
  periodic Rho snapshot, leaving mass drift uncomputable until the snapshot
  cadence was retuned per-config).
- Expected: VALIDATION T6 mass/momentum bands directly checkable from run
  outputs.
- Observed: manifest carries end-state `totalMass`/`maxSpeed` only; drift
  must be reconstructed from full-field Rho snapshots whose cadence must be
  hand-matched to an unpredictable steady-stop step. Works, but is fragile
  and costs one full-field write per sample (98 MB evidence dir for 18
  configs, dominated by snapshot fields).
- Proposal: manifest gains a periodic diagnostics series (step, totalMass,
  totalMomentum[, maxSpeed]) at `checkEvery` cadence — T6 then becomes a
  zero-I/O check for any agent. (Matches the observer/function-object
  framework direction, d66d0cb.)

**ANOM-P1-002 — T4 profile band is calibrated to the frozen ν and does not
transfer** — S3 (spec-documentation), disposition-proposal: VALIDATION
footnote.
- Scenario+config: `channel-zouhe-t4` variant at ν=0.05 (first authoring)
  vs frozen ν=0.02 (validation_open_bc.rs).
- Expected: central-profile L2rel ≤ 2e-3 (VALIDATION T4).
- Observed: 2.75e-3 at ν=0.05 (steady, 35k steps) vs 1.08e-3 at the frozen
  ν=0.02 — the O(Ma²) axial-pressure-gradient distortion grows with ν and
  leaves the band while the physics stays correct.
- Proposal: one-line footnote in VALIDATION T4 that the 2e-3 band holds at
  the frozen parameters (ν=0.02, u_max=0.05, H=32) — prevents false
  anomaly reports by future automated QA (this sweep now pins ν=0.02).

**ANOM-P1-003 — runtime |u| can exceed the compressibility advisory in
validated configs and nothing reports it** — S3 (monitoring gap),
disposition-proposal: runtime Ma watermark in manifest.
- Scenario+config: `cylinder-re100-karman-t8` (T8 2D-2 spec parameters,
  u_max=0.15 = exactly the validator's compressibility threshold).
- Expected: `lbm validate` warns above 0.15; the low-Mach hard limit is 0.3.
- Observed: runtime max|u| = 0.228 (1.52× the inlet peak, vortex
  acceleration) — validated silently, ran fine, T8 bands absorbed the
  compressibility error by construction. A user config validated at
  u=0.14 could reach ~0.21+ at runtime with no signal anywhere.
- Proposal: `maxSpeed` is already in the manifest — add the advisory
  cross-check at run end (warning entry when the runtime watermark crosses
  0.15/0.3), i.e. promote the validate thresholds to runtime. Cheap, catches
  the whole class. (A run that reaches 0.3+ currently just diverges or
  clips physics silently.)

### Harness errata (pass-1 self-findings, fixed in this commit)

- `xy_mirror_symmetry` initially mirrored about the half-integer domain
  center (x→nx−1−x); the droplet IC at integer (64,64) on a periodic grid is
  symmetric under x→(2cx−x) mod nx instead. Measured asym 0.564 (wrong axis)
  vs 7.1e-12 (correct axis) — a reminder that symmetry checks must use the
  IC's own symmetry group, not the array's.
- Snapshot cadence vs steady-stop (see ANOM-P1-001); `channel-zouhe-t4`
  re-pinned to the frozen T4 parameters.

### Coverage gaps (not anomalies — matrix items that cannot run today)

- **GAP-1 Taylor–Green vortex (T1)**: not expressible via scenario JSON — no
  analytic-init surface (`init` supports rest/droplet/pool only) and no
  pressure-consistent density init. T1 stays test-suite-only. Proposal: an
  `init: expression` or named-analytic-init scenario extension if agent-mode
  users are expected to reproduce validation cases.
- **GAP-2 buoyant droplet / bubble rise**: physically inexpressible.
  `physics.force` is a constant force *density* (VALIDATION T6: momentum grows
  N_fluid·F/step), so a uniform force is exactly balanced by a hydrostatic
  pressure gradient — no buoyancy without ρ-proportional gravity. The
  two-component path (`MultiComponent`, per-component gravity, used by T12) is
  not exposed through the scenario schema. Proposal: `force: {perMass: true}`
  or per-phase gravity in `MultiphaseSpec` when M-F needs it.
- **GAP-3 particles**: M-F spec only (FR-PART-*), nothing landed — track
  stubbed, per PM capability ruling.
- **GAP-4 3D multiphase, LES, rotating boundaries**: hard-rejected or absent
  on main (build3d rejects multiphase; no LES/rotating surface). Stubs only.
- **GAP-5 no Uz FieldKind**: 3D runs can output Speed/Ux/Uy/Rho/Vorticity but
  not Uz; z-velocity is only reachable via point probes or
  sqrt(Speed²−Ux²−Uy²) (sign lost). The strain/shear + 3D vorticity/Q
  FieldKind order in flight should include plain **Uz** for symmetry checks.
- **GAP-6 T9 horizon**: pass 1 ran the outflow-vortex case to 60k steps, not
  the spec's 1e5 (time cap). Full-horizon rerun queued for pass 2. Likewise
  the T9 pressure-oscillation-rms metric (≤15× central) and T9b
  ConvectiveOutflow comparison are pass-2 items.
- **GAP-7 2D-2 cylinder staircase parity**: the scenario `circle` obstacle is
  boundary-inclusive (d²≤r²) while the frozen T8 D=40 bands were measured with
  the exclusive staircase (d²<r²). Worked around with r=19.995 (lattice-exact
  equivalent); documented here so nobody "fixes" the radius in matrix.py.
- **GAP-8 vortex-center metric (T7)**: primary-vortex position (±0.02L) needs
  a streamfunction/ψ extremum; pass 1 checked centerline RMS only. Pass-2
  candidate (compute ψ by integrating ux over y from the Ux field).
- **GAP-9 f32 angles**: pass 1 is f64-only; T6 f32 drift bands and the T11
  f32 flat-interface angle are pass-2 items (precision switch is one line in
  matrix.py).

### Expected-limitation notes (validator advisories, working as intended)

- `cavity-re1000` (τ=0.538) and `cavity3d-re1000-n64` (τ=0.519, U/ν=16.1)
  draw the documented sub-0.55-τ / grid-Re advisories and still run stable —
  consistent with VALIDATION T7/T15.5 measured behaviour. `tau-margin-cavity-
  t10` (τ=0.51, U=0.05) is the frozen T10 stability-limit point; its
  stability (measured max|u| 0.0459 vs T10's frozen 0.046) is an acceptance
  criterion, not an anomaly.
- `cavity-re400`/`re1000` hit the 300k-step cap with the steady flag unset
  (spec's own cap; Ghia RMS passes regardless) — recorded for completeness.
