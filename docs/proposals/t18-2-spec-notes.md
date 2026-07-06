> **STATUS: spec absorbed into VALIDATION.md T18.2**

# SPEC_NOTES_T18_2

Notes from authoring the adversarial T18.2 acceptance tests from the frozen CR-2
surface only.

- The CR-2 API freezes rectangular `FacePatch` regions only. The impinging-jet
  acceptance text asks for a "coaxial annular/offset Outflow or Pressure patch";
  a true annulus is not expressible as one rectangular patch, so the test uses
  four non-overlapping rectangular Pressure patches to form a rectangular
  annulus around the central Velocity patch.
- The frozen API says `lo`/`hi` use in-face coordinates in remaining-axis
  ascending order, but it does not explicitly state whether `hi` is inclusive in
  validation math beyond the code comment. The tests assume `hi` is inclusive.
- The specification does not name patch-specific `SpecError` variants for
  out-of-bounds patches, overlapping patches, or one-open-axis violations from
  the base+patch union. The tests assert `is_err()` for those paths and use an
  exact `SpecError::VelocityTooHigh` match only for the already-public velocity
  guard.
- The public docs require GPU construction with non-empty `face_patches` to
  yield a `SpecError`, but existing public GPU construction patterns expose a
  `GpuSolver::try_new(...) -> Result<_, GpuInitError>` shape. The GPU test
  therefore checks that the error text surfaces either `SpecError` or
  `face_patches`; the exact wrapper variant should be frozen when CR-2 lands.
- Interactions between `WallSpec` solid rims and masked open patches on an
  otherwise closed face are underspecified. The impinging-jet tests mark the
  floor and side faces as wall rims, but leave the top face to be controlled by
  `FaceBC::Closed` plus `face_patches`, so top patches are not pre-empted by
  solid rim nodes.
- The acceptance text asks for "quasi-steady state" but does not define a
  convergence criterion. The default-suite impinging-jet test uses a fixed
  deterministic run length and freezes the mass-drift band at `2.0e-10` pending
  first measurement after CR-2 implementation lands.
- The radial wall-jet profile band is feature-based rather than reference-data
  based: axis stagnation minimum, positive off-axis peak, and monotone decay
  after the peak with a small `5.0e-5` slack. This is intentional because no
  reference profile or Reynolds number target is specified for the CR-2 gate.
- The patch/base equivalence test freezes an absolute field band of `1.0e-14`
  for f64 after 80 steps. If CR-2 applies the same operation in a different but
  round-off-equivalent order, this band may need first-measurement adjustment.
