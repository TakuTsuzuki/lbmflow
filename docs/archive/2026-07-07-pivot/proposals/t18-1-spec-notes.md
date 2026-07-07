> **STATUS: spec absorbed into VALIDATION.md T18.1**

# SPEC_NOTES_T18_1

Adversarial findings while authoring T18.1 from `docs/VALIDATION.md` and
`docs/DISPERSED_DEPOSITION.md` CR-1.

1. The frozen CR-1 API specifies `SourceRegion`, `SourceKind`, `VolumeSource`,
   and `GlobalSpec::sources`, but it does not name the fallible solver-builder
   entry point. The tests use `GlobalSpec::validate(3, solid)` for CPU
   validation errors and assume a future `GpuSolver::build(...) ->
   Result<_, SpecError>` for GPU rejection, because the acceptance criterion
   requires `SpecError`, not panic.
2. `SourceRegion::hi` is inclusive. The overlap test intentionally makes two
   regions share exactly one corner cell; an implementation that treats `hi` as
   exclusive will miss that error.
3. The sink positivity rule says local rho must stay positive, but does not
   define the exact bound for a multi-cell region. The test interprets
   `q_lu` as the total source strength uniformly distributed over the inclusive
   region and rejects a sink whose per-cell drain can exceed a unit-density
   local population in one step.
4. The analytic far-field statement does not define how to measure shell
   radius for a finite 2^3 source. The tests use the geometric center of the
   inclusive box at half-cell coordinates and shell-average `-u dot rhat` for a
   sink.
5. The far-field criterion says "after a quasi-steady window" but does not
   define convergence. The default test fixes 360 steps on 32^3 and the ignored
   heavy test fixes 1200 steps on 64^3; those are intentionally concrete so the
   first implementation has to expose whether the band is realistic.
6. The `Jet` momentum-flux contract does not specify whether the measured
   flux is instantaneous at the source, a global momentum ledger, or a plane
   integral downstream. The test uses the global total-momentum increment per
   step before wall feedback should dominate, expecting `dP/dt = q_lu * u`.
7. The CR-1 text says scope includes CpuScalar and CpuSimd, but T18.1 does not
   specify whether every numerical assertion must run on both. This file
   targets CpuScalar only except for partition invariance, matching the focused
   acceptance items and keeping the default suite light.
8. The GPU rejection requirement is clear semantically but not mechanically:
   if GPU construction currently requires a live adapter before spec validation,
   a headless `--features gpu` environment could fail for the wrong reason.
   The implementation should validate unsupported `sources` before any adapter
   work when practical.
