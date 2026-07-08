# External architecture & physics review — disposition record (2026-07-07)

**Lifecycle:** snapshot (dated, never edited).
**Reviewer input:** external CFD/numerics architect, static code+doc audit
(reviewer did not run the build/test gates). **Disposition by:** PM (Fable),
each finding verified against code and, where feasible, against a test run in
this session before accepting.

This record maps every finding to a **verdict** and a **disposition** so the
review is traceable and not re-discovered. Verdict classes follow the
reviewer's own taxonomy: (a) demonstrably wrong, (b) unproven claim,
(c) defensible but poorly recorded.

## Headline correction

The reviewer's top P0 (F-001, "D3Q27 open-face coverage claim is false") was
itself **backwards**, and the way it went backwards is the review's most
important structural lesson. README/LIMITATIONS were **correct**; the STALE
side was the internal spec (`VALIDATION.md`, `ARCHITECTURE_V2.md`), which still
said D3Q27 `Outflow`/`Convective` were rejected and D3Q27 was not
scenario-selectable. D3Q27 outflow/convective **landed 2026-07-07** (second
same-day landing after velocity/pressure) and pass their gates
(`d3q27_open_bc.rs` 6/6 green, verified this session). A hand-maintained
capability matrix drifted and made a competent reviewer emit a confident,
wrong P0 — corroborating the reviewer's own Architectural Verdict #1
("hand-written matrix will keep lying"). **Only GPU-path D3Q27 open faces
remain rejected.**

## Findings disposition

| ID | Sev | Verdict | Disposition |
|----|-----|---------|-------------|
| F-001 | P0→rejected | (a), reversed | Docs were the stale party, not the code. Fixed VALIDATION.md + ARCHITECTURE_V2.md (`cx/ext-review-docs`). D3Q27 outflow/convective are supported on CPU and tested. |
| F-002 | P0→resolved (retain) | (b) | Central-moment `-0.16\|u\|²` velocity correction. A concurrent session adjudicated this fully on `main` (PHYSICS.md 2026-07-07 falsification record + `docs/proposals/CUMULANT_GALILEAN_FIX.md`): **removing the term was tested and made the advected-TGV3D frame spread 2.5× worse** (1.05e-2 vs 4.20e-3), falsifying the removal hypothesis — so it compensates a real \|u\|²-scaling part of the D3Q19 defect and is retained. Standing verdict: **retained with a narrowed claim** — the residual anisotropic part is uncorrected and Galilean invariance at finite frame velocity is NOT established on D3Q19 (holdout stays `#[ignore = FINDING]`). This review's contribution: a `LIMITATIONS.md` collision row + README scope note surfacing that verdict to the trust boundary; the reviewer's underlying concern (coefficient not first-principles-derived) is real and remains the open question under adjudication in the fix proposal. No new PHYSICS.md verdict was added here — the concurrent record owns it. |
| F-003 | P0 | (a/c) | README overclaims fixed (`cx/ext-review-docs`): "physically rigorous" reframed as enforced policy with inventoried open items; WALE "by construction" → operator design intent, DNS/wall-treatment as validation-queue; "bit-reproducible" scoped to the covered matrix with the T13 two-pass probe blind spot named. |
| F-004 | P1 | (a) confirmed | Native moving-wall `wall_u` bypassed the low-Mach guard. **Fixed** (`cx/ext-review-boundary`): `Solver::try_new` now validates `wall_u` (finite + speed ≤ MAX_SPEED, NaN-safe, reusing `SpecError::VelocityTooHigh`/`NonFiniteParameter` — no new variant, no clamp). Regression test in `boundary_input_validation.rs`. |
| F-005 | P1 | (a) confirmed | `init_with` seeded raw (ρ,u) unvalidated. **Fixed** (`cx/ext-review-boundary`): fail-loud panic on ρ≤0/non-finite/\|u\|>MAX_SPEED with offending coordinates; `# Panics` documented. Regression test included. |
| F-006 | P1→P2 | (a) partial | `set_bouzidi_links` accepts out-of-range `qd`. Repro in the finding was wrong (qd=0 hits the qd<0.5 multiply branch, no div-by-zero), but NaN/negative/≥1 still extrapolate non-physically. Low-level validation hook, documented as such. **Not fixed here** — queue item; `cx/vv-bouzimoq` covers the mixed-qd force ledger but not qd range validation. |
| F-007 | P1 | (b) | TRT magic-Λ per-BC guidance. Doc/queue item. |
| F-008 | P1 | (b) | D3Q19 drops body-diagonal moments → isotropy limits. Corroborates making D3Q27 mandatory for anisotropy-sensitive 3D; `cx/vv-rotaniso` in flight. |
| F-009 | P1 | (b) | GPU absolute-physics gates too narrow. In flight: `cx/vv-gpuabs3d`. |
| F-010 | P1 | (b) | WALE not yet turbulence-predictive. Already bounded in LIMITATIONS §5; README now matches. `cx/vv-walex` in flight. |
| F-011 | P1 | (b/c) | Rotating IBM subsystem-level, not stirred-reactor fidelity. Already GATE in VV_MASTER_PLAN 1.2/1.6 (`cx/audit-ibm`, `cx/audit-rotor`, ANOM-P4-001/010). |
| F-012 | P1 | (b) | Multiphase dynamic fidelity under-validated. Already bounded in LIMITATIONS §6; `cx/mp-dynamics`, `cx/mp-hard`, `cx/vv-sparger` in flight. |
| F-013 | P1 | (b) | Outflow reflection weak (T9 15× band). Corroborates `cx/vv-convout`; reflection-coefficient metric is the genuinely new ask (queue). |
| F-014 | P1 | (a) confirmed | README claimed particle "adhesion-capture and resuspension closures" — **no such closures exist** (only Shan-Chen *wall* adhesion for contact angle; deposition is one-way Lagrangian floor-crossing bookkeeping). Fixed README (`cx/ext-review-docs`); matches DISPERSED_DEPOSITION.md §3 + LIMITATIONS §4. Upgraded from reviewer's (b/c) to (a). |
| F-015 | P2 | (c) | FP16 = capacity mode. Already correct in LIMITATIONS §3; README FP16 line states bands are frozen-to-measured. |
| F-016 | P2 | (b/c) | T13 blind spot (two-pass probe double-count). Now named in README bit-reproducible caveat. |
| F-017 | P2 | (c) | q-major SoA cross-backend ABI coupling — accepted architectural risk; treat layout as ABI. Queue: layout-contract test. |
| F-018 | P2 | (a) confirmed | ARCHITECTURE_V2.md stale D3Q27 rows. **Fixed** with F-001. |
| F-019 | P2 | (b) | MPI functional not scale-proven. Already RED in release table + LIMITATIONS §8. |
| F-020 | P2 | (c) | Step-ordering doc inconsistency. Queue: executable phase-diagram assertion. |
| F-021 | P2 | (b) | WASM/GUI parity unproven. Already W3 in VV_MASTER_PLAN. |
| F-022 | P3 | (c) | "cumulant" naming = central-moment. Now stated in LIMITATIONS §3 collision row. |
| F-023 | P3 | (c) | Some bands frozen-around-behavior. `cx/vv-gci` (T15.5 GCI) + benchmark-backlog address this. |

## Fixed in this session (landed on clean, main-based branches)

- `cx/ext-review-docs` — F-001, F-003, F-014, F-018, F-022, F-016 (README
  caveat), and the F-002 trust-boundary surfacing (LIMITATIONS row + README
  note; PHYSICS.md verdict left to the concurrent record). Docs only; verified
  against code + `d3q27_open_bc` 6/6.
- `cx/ext-review-boundary` — F-004, F-005. `solver.rs` guards +
  `boundary_input_validation.rs`. Gates run this session:
  `boundary_input_validation` 5/5, `backend_simd_equiv` 21/21 (no regression),
  plus init_with-heavy smoke/validation suites.

## Corroborated existing queue/in-flight (no new duplicate work created)

F-008→`cx/vv-rotaniso`, F-009→`cx/vv-gpuabs3d`, F-010→`cx/vv-walex`,
F-011→`cx/audit-ibm`/`cx/audit-rotor`, F-012→`cx/mp-dynamics`/`cx/mp-hard`,
F-013→`cx/vv-convout`, F-023→`cx/vv-gci`. Genuinely-new queue candidates:
T9 reflection-coefficient metric (F-013), Bouzidi qd range validation (F-006),
q-major layout-contract test (F-017), step-ordering phase-diagram assertion
(F-020).

## Reviewer's architectural verdict — PM position

1. **No generated capability ledger → hand-written matrices lie.** ACCEPTED as
   the top structural risk; F-001 is the proof. Candidate: generate the
   capability matrix from code (`LatticeSpec`, `FaceBC` support, feature gates)
   so docs cannot drift.
2. **CPU-relative equivalence ≠ absolute physics.** ACCEPTED; per-backend
   absolute gates are the `cx/vv-gpuabs3d` direction.
3. **cumulant naming + live `0.16` correction.** ACCEPTED; renamed in
   LIMITATIONS, `0.16` recorded as derive-or-remove.
