# SPEC — Bouzidi curved-boundary record list + STL voxelization

**Status split:**
- **Record-list infrastructure + circle/sphere links: LANDED.**
  See `crates/lbm-core/src/bouzidi.rs` (`BouzidiLink`, `BouzidiLinks`,
  `circle_links`, `sphere_links`, `half_way_links`, `apply_bouzidi_impl`).
- **STL voxelization pipeline: NOT LANDED.** Forward-looking; kept below.

License discipline: implement from Bouzidi/Firdaouss/Lallemand 2001 and
Guo/Zheng/Shi 2002. **No code** from OpenLB (GPL-2) or Palabos (AGPL-3.0) —
re-derived from published papers only.

---

## 1. Design invariants (as landed)

- Curved-BC correction is a **separate pass over a precomputed record list**,
  after streaming. Hot streaming kernel unforked.
- `SoaFields` carries an optional Bouzidi record buffer; None ⇒ every wall is
  the half-way rim (today's default).
- Records are sorted by `cell` for coalesced access.
- At qd = 1/2 the scheme degenerates **bitwise** to half-way BB (pinned by
  degeneracy test; blocks any regression in the flat-wall path).
- Wen 2014 Galilean-invariant momentum-exchange force at moving walls.

## 2. STL voxelization — forward-looking

The remaining work is a scenario-side pipeline that consumes an STL surface
mesh and emits a `BouzidiLinks` record list against the lattice:

1. Voxelize the STL: inside/outside classification per cell (ray parity /
   winding number; robust against non-manifold input).
2. For each fluid cell adjacent to a solid cell, compute `qd` per lattice
   direction by ray-tracing the STL along that link.
3. Emit the sorted record list into `SoaFields`.

Acceptance = 2nd-order wall-force convergence on a sphere-drag test vs
Schiller–Naumann, gated against the `sphere_links` analytic baseline.

Where the record-list contract lives: `crates/lbm-core/src/bouzidi.rs`.
GPU port slot: post-stream pass over the link list (kernel signature already
compatible with the WGSL generator).
