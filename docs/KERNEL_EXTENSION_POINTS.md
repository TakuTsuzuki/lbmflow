# Kernel Extension Points

This note records the R-Phase 2 B-8 extension contracts for kernel work that
must remain additive to the existing validated configurations. It is docs-only:
it describes where future kernels attach, which current code paths define the
contracts, and which gates are the minimum definition of done.

Primary references:

- `docs/SOLVER_IMPROVEMENT_SPEC.md` B-6 and B-8.
- `docs/reference/SPEC_COLLISION_COMPOSITION.md`.
- `docs/reference/SPEC_BOUZIDI_STL.md`.
- `docs/PHYSICS.md` entries for WALE LES and the 2026-07-06 cumulant track.
- Current code paths: `crates/lbm-core/src/fields.rs`,
  `crates/lbm-core/src/params.rs`, `crates/lbm-core/src/kernels.rs`,
  `crates/lbm-core/src/backend.rs`, `crates/lbm-core/src/backend_simd.rs`,
  `crates/lbm-core/src/bouzidi.rs`, `crates/lbm-core/src/gpu/wgsl.rs`, and
  `crates/lbm-core/src/compat/sim.rs`.

## 1. Per-cell omega passing convention

B-6 has landed as `SoaFields::omega_field: Option<Vec<T>>`, a compact-core
per-cell `omega_plus = 1 / tau_eff` field. The scalar collision path passes
the relevant row slice into `kernels::collide_row` and
`kernels::collide_row_central_moment`; the SIMD fused path carries the same
slice through its `FusedCtx`; the solver-level WALE driver writes it through
`Solver::set_omega_field`. `None` is a semantic and numerical contract: the
kernel must take the original uniform-omega arithmetic path exactly. Existing
configurations with no field must remain bit-identical to the pre-B-6 path.

The field layout is compact core order, not padded SoA order. For a global
field passed to `Solver::set_omega_field`, the solver scatters into each
subdomain's compact `g.cidx(x, y, z)` order. Consumers must not index this
field with padded population indices and must not expose it as part of the
public distribution layout. `omega_field[x]` replaces only the local symmetric
relaxation rate `omega_p`; the TRT antisymmetric rate `omega_m` remains the
configured value unless a future operator explicitly defines and validates a
different per-cell rule.

On GPU the equivalent storage is the `omega_out` storage buffer declared by
`gpu/wgsl.rs` at bind group binding 15. The current WALE path writes
`omega_out[i]` in the `wale_omega` entry point, and generated fused step
entry points such as `step_wale` and `step_wale_cumulant` read that buffer
when `FLAG_WALE` or the WALE-specific entry point is active. This is a
storage buffer of one `f32` value per device cell in the GPU cell order
`i = z * (nx * ny) + y * nx + x`, matching the generated shader's compact
device grid, not q-major population storage.

Any new relaxation-rate consumer, for example an LES-body order or a
Smagorinsky comparison, must preserve these invariants:

- `None` on CPU and disabled WALE/per-cell omega on GPU must be exactly the
  uniform-relaxation configuration, including Guo prefactors `cp` and `cm`.
- A uniform per-cell field equal to the configured `omega_p` must match the
  uniform configuration under the existing bit/tolerance gates for that
  backend path.
- The field value is an input to collision only. It must not mutate density,
  velocity moments, force fields, solid masks, wall velocities, or halo
  exchange state.
- The value must remain in the physical relaxation domain required by the
  selected collision operator. If the model needs clipping, the derivation,
  validity range, and validation gate belong in `docs/PHYSICS.md`; it must not
  be calibrated to an acceptance band.
- CPU partitioning and MPI decomposition must not change results: the field
  is compact-local after solver scatter, and halo exchange still exchanges
  populations, not relaxation fields.

Definition of done: existing field-off tests remain bit-identical, the uniform
field test remains equal to the scalar specification, and any new omega
producer has a physics entry plus a direct validation test for its null or
reference behavior.

## 2. Placement of MRT / cumulant kernels

The current selection point is `CollisionKind` in `params.rs`, lowered into
`KParams`. CPU scalar dispatch currently branches in `Backend::collide`: if
`KParams::cumulant` is true it calls
`kernels::collide_row_central_moment`; otherwise it calls `kernels::collide_row`
for BGK/TRT. GPU generation similarly has separate generated entry points for
TRT/BGK and central-moment variants in `gpu/wgsl.rs`.

Post-R2-C, the intended branch is the `CollisionKind` to `Collision` operator
selection at the abstract-arithmetic boundary from
`SPEC_COLLISION_COMPOSITION.md`. The hot loops should call one operator body
through an `Arith` interpreter: scalar CPU, block/SIMD CPU, and WGSL emission.
The branch belongs outside the per-cell arithmetic body, selecting the operator
implementation once for the run or generated entry point, not sprinkling
operator-specific arithmetic through streaming, moments, or boundary code.

MRT belongs as a moment-space `Collision` implementation behind the same
boundary. Its transformation matrices `M` and `M^-1` are lattice-owned static
tables or generated lattice tables, colocated with the collision operator or
with lattice-specific kernel support, never embedded independently in scalar,
SIMD, and WGSL copies. Relaxation vectors are operator parameters derived from
`CollisionKind`/`StepParams`, and conserved moments must keep relaxation rate
zero.

The post-MF-alpha stage-3 central-moment/cumulant path lands in the same
operator slot. The current CPU reference functions in `kernels.rs`
(`central_basis`, `central_phi`, `solve_moment_system`, and
`collide_row_central_moment`) document the stage-2 scalar reference, but they
are not the long-term placement for three backend copies. Stage 3 should move
the transform and relaxation sequence behind the Arith-compatible collision
body. Central-moment basis metadata and any `M`/`M^-1` or equivalent transform
generation must have one source of truth per lattice, including D3Q19's
reduced moment family and D3Q27's full tensor-product family.

The physics contract is the 2026-07-06 `PHYSICS.md` cumulant-track entry:
the implemented operator is currently a cascaded central-moment collision, not
a full logarithmic cumulant model. It transforms physical populations, relaxes
second-order deviatoric central moments with the selected shear rate, handles
trace/source prefactors as documented, and uses the same physical velocity
definition with Guo's half-force correction. A future full cumulant operator
must add its own derivation and validation rather than silently changing this
meaning under the existing name.

Definition of done: every future collision kernel must prove existing-
configuration bit invariance before it adds new physics. BGK/TRT default
scenarios must keep the current backend-equivalence tolerances
(`backend_simd_equiv`, T14 CPU/GPU where applicable, equilibrium fixed-point
tests). New MRT or cumulant configurations then add their own validation bands
without relaxing the existing ones.

## 3. Per-link wall-distance sparse structure of Bouzidi

Curved walls attach through a sparse record list, not by forking the streaming
kernel. The current data model is `BouzidiLinks<T> { records:
Vec<BouzidiLink<T>> }` in `bouzidi.rs`, stored as
`SoaFields::bouzidi: Option<BouzidiLinks<T>>`. Each record names a padded
fluid cell, a lattice direction from the fluid cell toward the wall, the
fractional wall distance `qd`, whether the second fluid node exists, and a
wall reference cell for wall velocity and probe membership. Records are sorted
by `(cell, q)`.

The SPEC_BOUZIDI_STL placement is:

```text
collide -> halo exchange -> stream -> Bouzidi correction -> swap
        -> open-face BC -> volume sources -> moments
```

`Backend::run_span` already places `apply_bouzidi` after all stream passes and
before `swap`. `bouzidi::apply_bouzidi_impl` reads post-collide populations
from `fields.f`, overwrites the incoming slots in `fields.ftmp`, and leaves
unrecorded links as the half-way bounce-back result written by streaming. This
keeps the hot pull-streaming path in `kernels::stream_row` and the fused CPU
path unchanged for the no-Bouzidi case.

The degeneracy contract is strict: `qd = 1/2` must reduce to the existing
half-way bounce-back formula bitwise. That means the record pass must compute
`f_in = f*_q + W_q` for that link, with the same moving-wall term convention as
`stream_row`, and must not introduce interpolation arithmetic on the exact
half-way branch. `set_bouzidi_half_way_links` exists as the harness for this
contract.

The STL-lattice extension point is record generation. Analytic builders
(`circle_links`, `sphere_links`) and a future STL voxelizer should both emit
the same `BouzidiLinks` structure. The voxelizer may use material tags, BVH
ray intersections, or triangle-derived wall groups, but the runtime kernel
must still see sorted link records with `qd` in `(0, 1)` and one writer per
`(q_bar, cell)` output slot. A ray miss or ambiguous geometry fallback must be
documented; using `qd = 1/2` is allowed only as a declared half-way fallback,
not as a hidden accuracy shortcut.

Definition of done: `bouzidi = None` keeps every existing staircase wall case
unchanged; all-`qd = 1/2` links match the half-way path bitwise; analytic or
STL Bouzidi cylinder validation lands T8 coarse drag in the documented
`+/-3%` Cd band and adds the corresponding force/momentum-exchange validation
for curved links.

## 4. ConvectiveOutflow alternative at in-place streaming

The current CPU scalar scheme is pull streaming. When a source lies beyond a
non-periodic, no-halo face, `stream_row` skips that population and leaves the
out-buffer slot untouched. After the population swap, the skipped face slots
therefore hold the previous step's post-collide values. `convective_face`
depends on this memory term:

```text
f(edge, t+1) = (f_prev_post_collide(edge) + Uc * f(interior, t+1)) / (1 + Uc)
```

`backend_simd.rs` preserves that same contract for fused streaming by
capturing open-face unknown slots into a per-face stale stash during the fused
stream pass, restoring the previous stash before `apply_open_faces_impl`, and
then rotating the stash generations. The GPU generated step uses the same idea
explicitly: `stash_in` supplies the skipped unknown slots in `f_out`, and
`stash_out` captures this step's post-collide unknowns for the next step. The
GPU `bc` pass then operates on the same logical post-swap state as CPU
`apply_open_faces`.

An in-place streaming or esoteric-pull scheme must generalize this edge-stash
contract, not reinterpret ConvectiveOutflow. The minimum abstraction is a
per-open-face memory source that provides the previous step's post-collide
unknown populations in canonical `Face::index()` and `L::unknowns(face)` order,
plus a capture point that records the current step's post-collide unknowns
before those slots are overwritten. Whether the backend obtains those values
from untouched in-place memory, a side buffer, or a generated shader stash is
an implementation detail; the BC pass must see the same values.

Ordering remains:

```text
collide -> halo exchange -> stream/capture skipped unknowns
        -> Bouzidi if present -> swap or logical parity flip
        -> restore previous unknowns -> open-face BC -> moments
```

For GPU-style push streaming, this is the existing edge-stash scheme. For
in-place or esoteric-pull variants, the same logical stash can be represented
as an alias-safe side band, a parity-indexed ring, or a face-local scratch
buffer, but it must be sized by the lattice unknown count and face extent, not
by ad-hoc D2Q9 assumptions.

Definition of done: enabling a new streaming backend with no open faces must
be covered by existing backend equivalence; enabling it with Outflow or
ConvectiveOutflow must keep existing T9/T9b results bit-invariant or within
the already established CPU/GPU tolerance. No T9/T9b acceptance band may be
relaxed to accommodate a changed stale-slot convention.

## 5. BC fallback of Outflow x solid adjacency

A-3's current permanent mitigation is runtime rejection. In the compat facade,
`Simulation::set_solid` rejects cells on an open edge and cells directly one
step inward from an open edge because open-face BCs read that interior
neighbour. The documented failure mode was a solid pocket directly inward from
an Outflow, causing the BC to skip reconstruction and freeze unknown
populations into a finite but non-physical steady velocity. V2 validation also
rejects invalid open-face geometry at build/patch boundaries.

The permanent solution should make this a documented boundary resolver rather
than an error. The resolver belongs between geometry/material classification
and `apply_open_faces`, where each open-face cell can classify the one-cell
inward neighbour used by Outflow or ConvectiveOutflow. If the neighbour is
fluid, the existing copy/convective formula is used. If the neighbour is
solid, the resolver must choose a physically documented fallback for that face
cell, for example treating the local segment as a wall closure or using a
declared one-sided extrapolation with a validity domain. The choice must be
encoded as face-cell metadata so CPU scalar, CPU SIMD, GPU BC, and future
patch/material paths all execute the same rule.

The resolver must not silently skip the cell, preserve stale unknowns, or
branch on scenario identity. It must account for corners and patches through
the same face-cell selection tables used by `for_face_cells_selected`, and it
must define what happens when only part of an outlet is blocked by a body. If
the fallback represents a wall, its wall location and moving-wall behavior
must be consistent with the half-way or Bouzidi wall model active on the
adjacent solid. If it represents an extrapolating outlet, its mass and momentum
effect must be measured and documented.

Definition of done:

- The current rejection tests remain until the resolver is implemented and
  validated.
- The replacement accepts the formerly rejected Outflow/ConvectiveOutflow x
  solid-adjacent geometry only when every affected face cell has an explicit
  resolver classification.
- The A-3 stationary-pocket reproducer no longer develops the non-physical
  frozen velocity, and mass/momentum behavior is covered by a dedicated
  validation case.
- Existing valid open-boundary cases remain bit-invariant, including T9/T9b
  and CPU/SIMD/GPU open-face equivalence.
- The fallback's physical rationale, validity domain, and behavior-review
  record are added to `docs/PHYSICS.md` before the runtime rejection is
  removed.
