# V&V Claim Status Guide

This document defines the claim-status vocabulary used by the V&V traceability
matrix, evidence report, and final campaign report. It is deliberately stricter
than "tests passed": a claim is valid only to the scope covered by its evidence.

## Status Vocabulary

| Status | Meaning | Minimum evidence |
|---|---|---|
| `VALIDATED` | The physical behavior is compared with an analytic solution, literature reference, or accepted benchmark, and the observed pattern is physically reviewed. | Passing command log, scalar metrics against a frozen band, visual/spatial artifact when fields exist, behavior-review verdict, and provenance for any closure. |
| `VERIFIED-ONLY` | The implementation is internally consistent, equivalent, deterministic, or regression-safe, but not independently validated against physics. | Passing equivalence/regression tests with command logs and clear scope limits. |
| `SPEC-ONLY` | The requirement or acceptance criterion is written, but no implementation or executable validation evidence exists. | Source section reference and missing implementation/test path. |
| `MISSING` | The claim has no adequate spec, test, or evidence. | Search/audit evidence showing no current coverage. |
| `BENCH-PENDING` | The claim requires hardware, runtime, cluster access, or long-duration validation not available in the current session. | Build/static evidence if available, plus exact pending command/environment/artifact. |
| `UNSAFE-CLAIM` | The claim is stronger than the evidence, rests on scalar-only evidence, relies on a banned pattern, or extrapolates outside a validity domain. | Contradicting evidence, missing visual/behavior review, unsupported product path, or physics-discipline violation. |
| `STOP-RULE` | The gate cannot be met without ad-hoc physics or another banned term. | Stop-rule report with attempted physical approaches and PM options. |

## Downgrade Rules

Apply the first matching downgrade.

| Condition | Required status |
|---|---|
| No executable test or run evidence exists | `SPEC-ONLY` or `MISSING` |
| Only CPU/GPU, scalar/SIMD, partition, or regression equivalence exists | `VERIFIED-ONLY` |
| A physical run emitted only scalar metrics while fields/spatial behavior exist | `UNSAFE-CLAIM` until visual artifact and behavior review exist |
| GPU adapter, MPI toolchain, cluster, or long-duration environment was unavailable | `BENCH-PENDING` |
| Implementation exists but feature is explicitly rejected or unsupported on a product path | `SPEC-ONLY`, `MISSING`, or `UNSAFE-CLAIM`, depending on the public claim |
| Evidence covers fluid only but the claim mentions FSI, particles, deposition, or multiphase coupling | `UNSAFE-CLAIM` for the broader claim |
| A gate would require calibrated constants, case-identity branches, silent clamps, or physical fallbacks | `STOP-RULE` |

## Evidence Requirements by Claim Area

| Claim area | Validation evidence required before `VALIDATED` | Common downgrade traps |
|---|---|---|
| D2Q9 fluid core | T1-T10 with release command logs, scalar bands, symmetry/behavior anchors, and PHYSICS.md rationale for known artifacts. | Treating green unit tests as validation when the benchmark or behavior anchor was not run. |
| D3Q19 fluid core | T15 physics tests plus T13 partition invariance and relevant diagnostics. | Claiming all 3D physics from z-invariant or equivalence tests alone. |
| D3Q27 | D3Q27-specific physics and isotropy/face-boundary validation. | Inheriting D3Q19 validation without D3Q27 evidence. |
| GPU/wgpu | T14 equivalence plus absolute GPU physics sentinels and runtime evidence on a named adapter. | CPU-relative equivalence only; adapter unavailable; f64 or unsupported physics silently assumed. |
| MPI/partition | T13 in-process and MPI runtime evidence, diagnostics, and command logs. | Monolithic-only tests; field-only match without diagnostic invariance. |
| f32/f16 precision | T6/T15/T16 precision-specific bands and documented degradation/validity domain. | Assuming f64 validation transfers to f32/f16 or long transients. |
| Multiphase | T11/T11b/T11c/T12 bands plus behavior anchors for interface shape, pressure plateau, spurious-current location, mass drift, and contact-angle trend. | Calibrated bands without closure provenance or spatial review. |
| Particle/deposition | T18 terminal velocity, source/sink, face-patch, deposition determinism, and artifact-backed spatial review. | Example/demo trends without model validity domain or visual evidence. |
| FSI/IBM/rotating boundary | Torque/slip/momentum benchmarks, Couette/Taylor-Couette or accepted reference cases, and explicit validity domain. | Fluid-only validation presented as FSI validation; coarse IBM characterization overstated. |
| Scenario/CLI/MCP | Schema validation, preset audits, CLI smoke, artifact output, and honest rejection of unsupported physics. | Product path silently using defaults or unsupported backend physics. |

## Claim Table Template

| Claim | Area | Status | Evidence | Missing evidence | Allowed wording | Forbidden wording |
|---|---|---|---|---|---|---|
|  |  |  |  |  |  |  |

## Review Checklist

- [ ] Every claim maps to a command, diff, log, or artifact from the current session.
- [ ] Every physical run that produced fields lists a visual artifact path.
- [ ] Every behavior claim has a mechanism statement and boundary/seam/outlet sweep.
- [ ] Every closure has a source, derivation note, validity domain, validation test, and PHYSICS.md entry.
- [ ] GPU, MPI, cluster, f16, FSI, and long-duration claims are not inferred from default CPU tests.
- [ ] Claims outside evidence are explicitly downgraded or forbidden.

