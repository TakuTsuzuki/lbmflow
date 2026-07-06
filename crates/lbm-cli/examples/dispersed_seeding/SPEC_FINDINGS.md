# SPEC_FINDINGS.md

This is the v5 record for `examples/dispersed_seeding`. The example remains a
small runnable protocol demo. CR-2/CR-3 substitutions have been replaced by core
APIs; the former example-local closure layer has been removed.

## Resolved v1 ambiguities

1. Reservoir suction is now explicit statistical extraction. The reservoir LBM
   domain is visualization-only; extraction samples a 1D settled concentration
   column at the requested depth.
2. Tray injection uses a top velocity face. Nozzle disks carry downward
   velocity; all other top-face nodes are set to zero velocity.
3. `eject.nozzle_diameter_m` is required and is used directly to compute
   `u_jet = Q / area`.
4. `depth_frac` is frozen as `0 = filled-liquid surface`, `1 = reservoir floor`.
5. Particles are integrated until floor crossing or protocol end. Suspended
   particles are reported as `n_suspended` and are excluded from the density
   map. `max_particle_steps` aborts rather than silently changing the physics.
6. SI dimensions are authoritative. Validation fails if supplied grid counts
   differ from `SI / dx_m` by more than one cell. The reservoir keeps the
   coarse visualization spacing `grid.dx_m`; the tray can specify
   `grid.tray_dx_m` so low-Mach tray injection does not force a costly reservoir
   refinement.
7. Unknown `agitate.pattern` values are validation errors.
8. VTK output is ASCII `STRUCTURED_POINTS` with vector data in grid order.
9. Multi-point ejection is supported through `points_xy_frac`.
10. The gentle spread band is now quantitative, not qualitative.

## Closure-removal outcome (2026-07-06)

Order `cx/kill-deposition-closures` removed the example-local transport
closures: no harshness branch, no analytic jet/wall-jet superposition, no
stochastic lateral dispersion coefficient, no side-wall particle clamps, no
direct agitation kicks, no reservoir scoring heuristic, and no mystery
reservoir body force. Particles are stepped together with the tray solver and
sample the live resolved velocity through trilinear interpolation of the
current solver moment fields. Translational agitation is applied as an
oscillating frame acceleration through the core per-mass Guo forcing path, with
the matching density-weighted particle pseudo-force routed through
`ParticleSet::g`.

Reservoir withdrawal is now a statistical one-dimensional settling-column
sample at `depth_frac`: for each diameter, the concentration at the withdrawal
height is nonzero when the Stokes/T18.3 settling backtrace lies inside the
initially filled column. No rate band, settled bonus, or score weight remains.
Nozzle exit positions are sampled uniformly over the physical nozzle disk.

First resolved-only gentle run:

```bash
cargo build --release -p lbm-cli --example dispersed_seeding
./target/release/examples/dispersed_seeding crates/lbm-cli/examples/dispersed_seeding/sample_gentle.json
```

Result:

- `RESULT CV=0.000000 max_over_mean=0.000000 empty_bin_fraction=1.000000 n_deposited=0 n_suspended=10000 n_extracted=10000 out=out/dispersed_seeding/gentle`
- Wall time: `1607.25 s`, above the order's ~15 min/sample redesign threshold.
- Artifacts:
  - `out/dispersed_seeding/gentle/density.csv`
  - `out/dispersed_seeding/gentle/density.png`
  - `out/dispersed_seeding/gentle/tray_velocity.vtk`
  - `out/dispersed_seeding/gentle/near_floor_radial_velocity.csv`

Behavior consequence: the old CV/empty-bin bands are invalidated. With only
resolved transport at the current tray resolution and protocol duration, the
gentle sample deposits no particles; all extracted particles remain suspended.
The former edge ring does not survive because there is no deposition map. The
near-floor radial profile is weak (order `1e-6 m/s`) and does not transport
particles to the floor over the sample duration. This is a PM/core capability
gap for the demo budget/protocol or resolved jet/free-surface capability, not a
reason to restore the deleted closure layer.

## Historical experimental outcomes (closure-driven, invalidated)

Commands run on 2026-07-06:

```bash
cargo build --release -p lbm-cli --example dispersed_seeding
./target/release/examples/dispersed_seeding crates/lbm-cli/examples/dispersed_seeding/sample_gentle.json
./target/release/examples/dispersed_seeding crates/lbm-cli/examples/dispersed_seeding/sample_harsh.json
```

Low-Mach scaling:

- `grid.tray_dx_m = 1.875e-4 m`, refined from the old shared `dx_m = 5e-4 m`.
- `dt = 4.21875e-4 s`, computed by diffusive scaling
  `dt = 0.012 * tray_dx_m^2 / nu_phys`.
- `nu* = 0.012`, so physical viscosity stays `1.0e-6 m^2/s` and
  `tau = 3*nu* + 0.5 = 0.536`.
- Physical inputs are unchanged: gentle jet velocity is
  `2.387324e-2 m/s`, harsh jet velocity is `2.546479e-2 m/s`, gravity remains
  `9.80665 m/s^2`, and Stokes settling remains `2.7241e-4 m/s`.
- Reservoir grid remains `32x32x200` at `dx_m = 5e-4 m`; tray grid is
  `256x256x64`.

Observed metrics:

| sample | CV | max/mean | empty bins | deposited | suspended | extracted | Re_jet | St | Fr | Ma | tau |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| gentle | 1.138149 | 5.053599 | 0.0 | 7836 | 2164 | 10000 | 19.1 | 0.01741 | 0.0 | 0.093 | 0.536 |
| harsh | 4.147195 | 25.345883 | 0.8333333333333334 | 17743 | 257 | 18000 | 127.3 | 0.002971 | 0.165 | 0.099 | 0.536 |

Acceptance consequences:

- Both samples are in the low-Mach band: gentle `Ma = 0.09304`, harsh
  `Ma = 0.09924`; both retain `tau = 0.536 >= 0.51`.
- Gentle spreading is achieved: `empty_bin_fraction = 0.0`, below the target
  of `0.15`.
- Trend holds: `CV_gentle = 1.138149 < CV_harsh = 4.147195`.
- Re-frozen gentle CV band for this low-Mach example and sample set:
  `1.05 <= CV <= 1.30`. The old `0.95 <= CV <= 1.40` band was measured at
  `Ma ~= 0.25` in the compressibility-error regime and is retired.
  This band is invalid after closure removal and must not be used for
  resolved-only acceptance.

## Core requirements represented by substitutions

- Localized source/sink or withdrawal boundary for a reservoir, including an
  internal suction plane and replenishment/open top behavior.
- CR-2 replacement is complete for this example: tray injection uses the core
  masked face-patch API instead of a local mixed-face boundary substitution.
- Public D3Q19 support for localized nozzle area/flow-rate boundary conditions
  that preserve the specified volumetric flow after lattice discretization.
- Resolved transport capability sufficient to carry particles from the nozzle
  region to the deposition floor within the protocol budget, or a PM-approved
  protocol/budget revision. No example-local closure is currently allowed.
- CR-3 replacement is complete for this example: suspended/deposited accounting
  and floor-crossing events come from the core particle layer. The sampler is a
  deterministic function of particle position and step; stochastic dispersion is
  intentionally not part of core sampling.

## New issues

- The resolved-only gentle sample currently produces no floor deposits within
  the protocol duration. Future acceptance specs should decide whether the
  protocol duration, resolved tray capability, or free-surface/jet model needs
  revision before new bands are frozen.
