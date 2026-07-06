# SPEC_FINDINGS.md

This report was produced while implementing and running the minimal
`examples/dispersed_seeding` demo. The sample runs are experiments, not just
compile checks: the implementation was kept small so failures and substitutions
would expose gaps in the order specification.

## Findings

1. **Reservoir suction is not expressible through the public D3Q19 API.**
   The requested setup needs a localized suction outlet at an interior depth and
   open replenishment at the top. The public `GlobalSpec` exposes whole-face
   boundary conditions and rejects open faces on multiple axes. Proposed fix:
   specify either a new public localized source/sink or make the MVP extraction
   model explicitly statistical.

2. **The target tray asks for a top open outlet and one or more top velocity
   inlets on the same face.** The current face-BC model supports one boundary
   kind per face, although `set_inlet_profile_with` can localize velocity on a
   velocity face. Proposed fix: define a mixed top boundary mask, or specify
   that the MVP uses a top velocity face with zero velocity outside jet patches.

3. **Flow-rate to inlet velocity is underdetermined.** `rate_uLs` is provided
   but nozzle diameter or jet patch area is not. The example uses an effective
   patch radius derived from grid spacing and tray size. Proposed fix: add
   `nozzle_diameter_m` or `inlet_patch_radius_m` to each `eject` point.

4. **`depth_frac` origin is ambiguous.** It is not stated whether `0` is the
   floor, open top, or filled-liquid surface. The example treats `0` as the
   reservoir floor and `1` as `fill_height_m`. Proposed fix: define the origin
   and whether the coordinate is absolute reservoir height or filled height.

5. **Protocol physical durations are too expensive for an example-scale
   release run after Mach-safe nondimensionalization.** The harsh sample can
   imply many thousands of particle steps depending on the chosen inlet area.
   The example caps particle integration at 900 steps and force-deposits any
   survivors at the end. Proposed fix: add an example-mode step cap or provide
   smaller sample counts/durations for CI-style acceptance.

6. **Expected CV trend is qualitative but not quantitatively specified.**
   The order says gentle should have low CV and harsh higher CV, but no
   acceptance margin is stated. Proposed fix: freeze an expected band after the
   physical model is made unambiguous.

7. **Reservoir dimensions are duplicated by SI geometry and grid counts.**
   `height_m / dx_m` and `width_m / dx_m` do not have to match `res_n*`.
   Proposed fix: define whether grid counts or SI dimensions are authoritative
   when they disagree.

8. **Agitation `pattern` has only one named value.** The spec defines
   translational agitation but does not state whether unknown patterns should
   error, warn, or be ignored. Proposed fix: require validation failure for
   unsupported patterns.

9. **Deposition boundary behavior is underspecified after the requested settle
   duration.** Some particles may remain suspended in a reduced model. The
   example freezes survivors at the final floor projection to keep readout
   finite. Proposed fix: define whether `n_deposited < n_extracted` is valid or
   whether the run must continue until all particles deposit.

10. **VTK coordinate convention is not specified.** The order accepts VTK but
    does not say whether points or cells should be used. The example emits
    `STRUCTURED_POINTS` point vectors in compact grid order. Proposed fix:
    define the volume schema expected by downstream QA viewers.

## Experimental outcomes

Runs performed with:

```bash
cargo run --release -p lbm-cli --example dispersed_seeding -- \
  crates/lbm-cli/examples/dispersed_seeding/sample_gentle.json
cargo run --release -p lbm-cli --example dispersed_seeding -- \
  crates/lbm-cli/examples/dispersed_seeding/sample_harsh.json
```

Observed metrics:

| sample | CV | max/mean | deposited | extracted | Ma | tau |
|---|---:|---:|---:|---:|---:|---:|
| gentle | 2.0144649711523934 | 8.4384 | 10000 | 10000 | 0.017847199124801216 | 0.536 |
| harsh | 5.301738079787294 | 33.128 | 18000 | 18000 | 0.13856406460551018 | 0.5115838486046297 |

The implementation writes finite `metrics.json` for both sample inputs and
emits `reservoir_velocity.vtk` plus `tray_velocity.vtk` for each sample. The
harsh input produces a more concentrated deposition pattern than the gentle
input in this reduced model (`CV_gentle < CV_harsh`), so the qualitative trend
gate is exercised by a real end-to-end run rather than by editing sample
parameters.
