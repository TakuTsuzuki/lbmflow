# Behavior-validity review — example sweep (2026-07-07, PM-run)

Case: lid-driven cavity spin-up transient (1000 steps, NOT steady state),
3 viscosities x 2 resolutions, one repeat. Runs executed natively by the PM
(the authoring codex sandbox had no built binary).

Observed pattern:
- Total kinetic energy (extensive, sum over cells) increases monotonically
  with `nu` at fixed N (N=33: 0.195 -> 0.250 -> 0.290; N=49: 0.349 -> 0.465 ->
  0.582).
- `maxSpeed` stays below the lid speed and increases slightly with `nu`.
- `tau = 3 nu + 0.5` reproduced exactly (0.53 / 0.56 / 0.62).

Dominant mechanism: during spin-up the flow is set in motion by momentum
diffusion from the moving lid; larger `nu` diffuses lid momentum into the
interior faster, so at a fixed early time more fluid is moving — total KE
rises with `nu`. This is resolved transient physics (no closures, no boundary
artifacts involved); signs, monotonicity, and bounds are consistent with the
mechanism. Note the trend is specific to the transient window: at steady
state the KE-vs-Re relation is not trivially monotone, so this example must
not be quoted as a steady-state sensitivity.
