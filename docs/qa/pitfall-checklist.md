# LBM Implementation Pitfall Checklist (V&V lane 4.1)

Compendium of known LBM implementation pitfalls from literature and community
knowledge, each mapped against `crates/lbm-core` test coverage. Compiled
2026-07-06 from web research + read-only repo audit (research agent; PM
committed). Verdicts: **COVERED** (existing test would catch it), **PARTIAL**
(some exposure remains), **GAP** (no coverage; kill-case proposed), **N/A**
(feature absent). Tally: 16 COVERED / 8 PARTIAL / 1 GAP / 1 N/A / 1
covered-by-design across 27 entries.

Method note: verdicts grounded by two repo sweeps plus direct reads of
accuracy_audit.rs, kernels.rs (stream_row probe form, line ~536),
t16_fp16_storage.rs, VALIDATION.md. All references verified against search
results; items marked "community folklore" have no reliable citation.

## Summary matrix

| # | Pitfall | Verdict | Anchor / evidence |
|---|---------|---------|-------------------|
| 1 | Mixed forcing-scheme conventions (Guo vs SC-shift vs He/EDM) | **PARTIAL** | forcing_path one-step match + ANOM-P2-001 pin; MISSING: force-driven vs pressure-driven Poiseuille equivalence |
| 2 | Missing F/2 in velocity moments | COVERED | design invariant + trt_magic_is_exact (impossible with half-force bias) |
| 3 | Non-equilibrium initialization transients | **PARTIAL** | pressure-consistent init used in T1; no direct ringing assertion |
| 4 | Zou-He corner/edge underdetermination | COVERED (by design) | SpecError::OpenFacesOnMultipleAxes forbids the corner; 4-orientation T4/T5 |
| 5 | Checkerboard / odd-even decoupling (staggered modes) | **GAP** | kill-case: seeded (−1)^(x+y) density mode must decay per collision spectrum |
| 6 | tau→0.5 over-relaxation fragility | COVERED | t10 stable-at-0.51 + documented divergence at grid-Re≈30 |
| 7 | Viscosity-dependent BB slip (BGK) vs TRT Λ=3/16 | **PARTIAL** | TRT exact; BGK slip-law pin is #[ignore] SPEC-GAP (convention freeze pending) |
| 8 | Half-way wall off-by-one (effective width) | COVERED | t2 exact at two widths |
| 9 | Moving-wall BB momentum term | COVERED | Couette all-tau exact; cavity 4-orientation; MW mass drift |
| 10 | Momentum-exchange Galilean invariance (Wen 2014) | **PARTIAL** | conventional Ladd form; static probes safe; kill-case: co-moving walls must read zero force |
| 11 | Bouzidi mass non-conservation | **PARTIAL** | no total-mass assert with Bouzidi active; kill-case: closed channel + off-grid cylinder drift pin |
| 12 | Fresh-node refill at moving solids | N/A | no moving voxel walls; moving geometry is diffuse IBM |
| 13 | Staircase boundary artifacts | COVERED | rotor half-cell translation pins; T8 staircase-vs-Bouzidi bands |
| 14 | Outflow reflections / vortex crossing outlet | COVERED | t9/t9b frozen pressure-ratio pins; no non-reflecting BC = feature choice |
| 15 | Ma² compressibility contamination | COVERED | cavity same-Re half-Mach scaling (light+heavy); TGV diffusive scaling |
| 16 | O(u³) equilibrium-truncation Galilean defect | COVERED (pinned) | advected-TGV defect pinned ≤1.5e-1 (loose by design, regression detector) |
| 17 | Cumulant/central-moment GI claims vs measured | COVERED | cumulant_acceptance::measure_galilean; viscosity-offset separately audited (ANOM-P4-008) |
| 18 | D3Q19 rotational-symmetry deficit | **PARTIAL** | no D3Q19-vs-D3Q27 anisotropy discriminator; bouzidi g6 SPEC-GAP related |
| 19 | Shan-Chen spurious currents | COVERED | flat_interface max\|u\| pin + Laplace R² |
| 20 | SC thermodynamic inconsistency / forcing dependence of coexistence | **PARTIAL** | coexistence frozen to measurement (honest) but no tau-sweep pin of the dependence surface |
| 21 | LES near-wall eddy viscosity | **PARTIAL** | WALE-only; defining pure-shear nu_t≡0 property untested directly (wale_null covers Couette/Poiseuille — verify scope) |
| 22 | Acoustic dispersion/damping errors | COVERED | sound speed 1e-3 + damping envelope (band flagged loose in lane 2.2) |
| 23 | Single/half precision round-off | COVERED | deviation storage + t15_3d_f32 + d5_long_horizon + t16 (GPU-gated) |
| 24 | Halo / parallel synchronization | COVERED (exemplary) | T13 bit-exact splits + seam adversarial |
| 25 | Streaming / OPP[] table bugs | **PARTIAL** | stream_contract covers faces; no single-delta propagation unit or OPP[] table test |
| 26 | Silent clamps / NaN scrubbing | COVERED | no physics clamps found; run_guarded watchdog; ban-list discipline |
| 27 | Unit-conversion / parameter-range misuse | COVERED (boundary) | scenario validator (tau/Ma/grid-Re warnings), MAX_SPEED rejection |

## Known-defect ledger (pinned, intentional)

ANOM-P2-001 uniform-force impulse (#[ignore] pin until R2-C); BGK slip-law
SPEC-GAP; Stokes-II SPEC-GAP (no runtime MovingWall setter); Bouzidi g6
rotational-anisotropy SPEC-GAP; Galilean TGV band deliberately loose; T9
O(Ma²) outlet stagger documented intrinsic. Plus open ANOM-P4-001/008/010
(core routings, gates on cx/audit-ibm / cx/audit-cumulant / cx/audit-rotor).

## Proposed new kill-cases, priority order (feed lane 1.7 / W2 orders)

1. **P5 checkerboard-mode decay** (the one GAP): periodic box, seed
   rho = 1 + 1e-4·(−1)^(x+y), u=0; BGK and TRT at tau ∈ {0.6, 0.51};
   project rho onto the mode each step; assert monotone decay at the
   collision-spectrum rate, no growth/plateau above round-off, no leakage
   into (π,0)/(0,π) modes above 1e-12. 64², ~500 steps, sub-second.
   Suggested file: tests/accuracy_audit_modes.rs.
2. **P10 co-moving-frame zero wall force**: periodic box, uniform u0, both
   walls MovingWall at u0 (rigid co-moving, zero shear): physical wall force
   is exactly 0; conventional momentum exchange reports F ∝ rho·u0 per link
   if frame-dependent. One step. Decides the Wen-2014 question here.
3. **P21 WALE pure-shear nu_t ≡ 0**: laminar Couette WALE-on must be
   bit-identical (≤2 ULP) to WALE-off with nu_t ≤ round-off everywhere
   including first fluid cell (verify wale_null scope covers the field, not
   only integral behavior).
4. **P25 single-delta streaming + OPP[] table test**: per lattice, assert
   c[OPP[q]] == −c[q], w[OPP[q]] == w[q]; inject δ at one cell/direction,
   stream once without collision, assert arrival exactly at x+c_q only.
   Milliseconds; pins the storage-order contract SIMD/GPU fusions preserve.
5. **P11 Bouzidi mass-drift pin**: closed periodic-x channel + off-grid
   cylinder, 10⁴ steps, freeze |ΣΔrho|/Σrho band vs staircase reference
   (expected nonzero — bound and monitor).
6. **P20 SC coexistence tau-sweep**: flat interface at tau {0.6, 1.0, 1.4},
   assert coexistence-density shift within a frozen band (pins the
   forcing-scheme-dependence surface). Maxwell-rule deviation = doc task.
7. **P3 init-ringing comparison**: TGV flat-rho vs pressure-consistent init;
   acoustic k-mode amplitude ratio must match the O(u0)/O(u0²) theory.
8. **P1 force-driven vs pressure-driven Poiseuille equivalence** — cheapest
   cross-forcing-convention detector (schedule after ANOM-P2-001 fix).
9. **P18 D3Q19-vs-D3Q27 rotated-flow anisotropy discriminator** — heavier;
   schedule with the D3Q27 validation track.

## Key references

Guo/Zheng/Shi PRE 65 046308 (2002); Li/Luo/Li PRE 86 016709 (2012); Skordos
PRE 48 4823 (1993); Mei/Luo/Lallemand/d'Humières CF 35 855 (2006); Zou&He
PF 9 1591 (1997); Hecht&Harting JSTAT P01018 (2010); Lallemand&Luo PRE 61
6546 (2000); He/Zou/Luo/Dembo JSP 87 115 (1997); Ginzburg&d'Humières PRE 68
066614 (2003); d'Humières&Ginzburg CMA 58 823 (2009); Wen et al. JCP 266 161
(2014); Bouzidi/Firdaouss/Lallemand PF 13 3452 (2001); Lallemand&Luo JCP 184
406 (2003); Shan PRE 73 047701 (2006); Wagner PRE 74 056703 (2006);
Nicoud&Ducros FTC 62 183 (1999); Premnath et al. PRE 79 026703 (2009);
Geier et al. CMA 70 507 (2015), JCP 348 (2017); Holdych et al. JCP 193 595
(2004); Lehmann et al. PRE 106 015308 (2022); Krüger et al., The LBM:
Principles and Practice (Springer 2017); Zanetti (1989) staggered
invariants; OpenLB/Palabos forum unit-conversion threads (community).
