# SPEC_NOTES_T18_3

## Ambiguities and traps

- `VALIDATION.md` lists the T18.3 terminal-velocity test under
  `t18_3_particle_deposition.rs`, but this order requires the runnable
  pre-CR-3 terminal-velocity test to live in
  `t18_3_particle_settling.rs` so it is not blocked by missing deposition API.
- The phrase "swept over particle diameters spanning >= 2 decades of Stokes
  number, each case at Re_p < 0.1 and at least one case at moderate Re_p
  (~1-10)" is internally tense. The runnable test interprets this as a low-Re
  Stokes-limit sweep spanning many decades of particle response time plus one
  additional moderate-Re Schiller-Naumann case.
- The CR-3 spec says "step limits are caller-controlled and abort, never
  truncate silently" but does not name a step-limit parameter or error type in
  the suggested API. The future deposition tests therefore avoid step-limit
  behavior and focus only on frozen capture/determinism semantics.
- "Floor plane" is frozen as `floor_z`, but the current solid-contact
  documentation describes generic solid sampling rather than a distinguished
  floor plane. The CR-3 tests treat deposition as crossing `z == floor_z`,
  independent of `Sample::solid`.
- The exact payload in `DepositEvent { pos, particle }` is not fully specified:
  whether `particle` should be the pre-step particle, the interpolated
  crossing-state particle, or the removed post-drag particle. The tests only
  require the pre-step particle for deterministic index-order checks; crossing
  position is asserted via `event.pos`.
- `ParticleSet::step` uses `dt = 1` implicitly. The CR-3 suggested API also
  has no explicit `dt`; the tests assume one call is one lattice time step and
  that interpolation is over that unit step segment.
- The public docs for existing `ParticleSet::step` promise solid handling via a
  staircase-wall model with reflection/resting contact. CR-3 deposition instead
  requires removal and record emission on floor-plane crossing. That is a
  semantic extension, not a behavior promised by the existing public docs.
- The current public docs promise exposure accumulation through an
  `exposure_rate` closure before ordinary movement. CR-3 keeps the same
  closure in the suggested signature, but does not state whether a particle
  that deposits mid-step receives full-step, partial-step, or no exposure
  accumulation. The deposition tests pass `None` for exposure to avoid
  overspecifying this missing rule.
- The existing `ParticleSet::step` docs do not promise stable event ordering
  because there are no events. CR-3 explicitly freezes deterministic
  index-order recording; the deposition tests assert this directly.
