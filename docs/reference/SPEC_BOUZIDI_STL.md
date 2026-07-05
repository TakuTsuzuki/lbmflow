# SPEC — Bouzidi curved-boundary record list + STL voxelization

> Research-agent deliverable, P1 #2 (biggest physics-accuracy leverage). Ready-to-dispatch,
> same depth/format as SPEC_UNIT_CONVERTER.md and SPEC_COLLISION_COMPOSITION.md. Grounded in
> firsthand reads of the live boundary code — `kernels.rs::stream_row` (half-way BB, L269–297),
> `gpu/wgsl.rs` push-BB (L374–408), `fields.rs::SoaFields` (L168–223), and the T8 acceptance
> test `tests/validation_cylinder.rs`.
>
> **License discipline (mandatory, embedded in the codex order):** implement from the ORIGINAL
> papers — Bouzidi, Firdaouss & Lallemand 2001 (Phys. Fluids 13, 3452, interpolated bounce-back);
> Guo, Zheng & Shi 2002 (Chin. Phys. 11, 366, for the moving-wall/force terms). **No code from
> OpenLB (GPL-2) or Palabos (AGPL-3.0)** — designs are re-derived from the literature only.

---

## 1. Scope & the problem being solved

Today all walls are a **half-way bounce-back off a 1-cell solid rim** (CLAUDE.md invariant; the
inline BB in `stream_row` L272–282 and the fused/GPU push). For a curved body voxelised to solid
cells, the true wall sits at a per-link fractional distance `q ∈ [0,1]`, but half-way BB forces
`q = 1/2` on every link → a **staircase** whose wall-force error is **O(Δx)** (first order). The
cost is visible in the T8 (Schäfer-Turek) suite:

- `t8_2d1_d20` (D=20 staircase) accepts a **wide** band `Cd ∈ 5.2..6.0` (±~7% around
  `CD_REF_2D1 = 5.5795`), `Cl ∈ −0.05..0.08` — the band is wide *because* the staircase is only
  first-order at the wall (validation_cylinder.rs L166–181).
- `t8_2d1_d40` only asserts drag **converges** toward the reference (`err40 < err20`,
  L107–120) — i.e. you must *refine* to shrink the error.

**Bouzidi interpolated bounce-back** restores the true wall position per link (**2nd order**),
so the coarse D=20 case can hit a tightened band without refinement, and the wall force converges
at ~2nd order. This is the single highest-leverage physics upgrade (COMPETITOR_ANALYSIS §B: OpenLB
flagship gap, Palabos ships a whole off-lattice family). It **degrades gracefully to today's
behaviour** at `q = 1/2` (OpenLB's q=0.5 fallback), so it is additive, not a rewrite.

The design keeps the hot streaming kernel **unforked**: the curved-BC correction is a **separate
pass over a precomputed record list**, mirroring how open-face BCs already run as post-passes.

---

## 2. Geometry representation (records, not per-cell objects)

Geometry enters as **staged precomputed data** (COMPETITOR §B1 — both OpenLB and Palabos do this;
the anti-pattern to avoid is Palabos's per-cell virtual `Dynamics*`, COMPETITOR §D1). Add to
`SoaFields` one optional field, allocated only when a curved body is present (zero cost otherwise,
like `force_field`/`probe` today, fields.rs L192–194):

```rust
/// Bouzidi wall-distance records for curved boundaries. None ⇒ every wall is
/// the half-way rim (today's behaviour). Sorted by `cell` (coalesced access),
/// then by `q`. Built at scenario load; read by the post-stream Bouzidi pass.
pub bouzidi: Option<BouzidiLinks<T>>,

pub struct BouzidiLinks<T> {
    /// One record per (fluid boundary cell, wall-crossing direction q).
    records: Vec<BouzidiLink<T>>,
    /// CSR-style row offsets by cell for the scatter kernels (optional; the
    /// flat sorted list already coalesces).
}

struct BouzidiLink<T> {
    cell: u32,      // padded index i of the FLUID node x_f (same indexing as f: q*np + i)
    q: u8,          // direction pointing from x_f toward the wall
    qd: T,          // wall distance fraction ∈ (0,1): |x_f → wall| / |x_f → x_solid| = |x_f→wall|/Δx
    has_second: bool, // whether x_ff (second fluid node, x_f − c_q) exists & is fluid (q<1/2 path)
    wall_ref: u32,  // index into wall_u for the moving-wall term (reuses existing wall_u storage)
}
```

`qd` semantics (Bouzidi 2001 convention): along the link from fluid node `x_f` in lattice direction
`q` toward the solid, `qd = δ/Δx` where `δ` is the distance from `x_f` to the wall intersection.
`qd = 1/2` ⇒ wall exactly half-way ⇒ the formula reduces to today's half-way BB.

**Optional material-number grid** (COMPETITOR §B6, OpenLB `rename` / Palabos voxel flags): an
`Option<Vec<u8>>` tagging fluid / solid / inlet / outlet / named-wall, so BCs and force groups bind
to a material number rather than ad-hoc masks. Additive; not required for the first cylinder cut.

**Two builders** (same record output; the pass is builder-agnostic):
- **Analytic primitives (do first):** circle / sphere / rect / capsule already in the scenario. For
  these `qd` is **closed-form exact** — a ray-vs-circle solve per crossing link — so the
  Schäfer-Turek cylinder (a circle) tests 2nd order with **no STL needed**. This is the fastest path
  to the T8 win.
- **STL (do second):** see §4.

---

## 3. The Bouzidi post-stream pass

**Placement in the step:** `collide → halo → stream → [Bouzidi correction] → open-face BC →
moments`. It runs **after streaming, before swap**, so it reads **post-collide** populations `f*`
from the in-buffer `f` (the pull scheme leaves `f` = post-collide until swap; this is the same
state `stream_row`'s BB reads via `f[opp(q)*np+i]`, L276) and **overwrites** the incoming-population
slots that streaming already filled with half-way BB in the out-buffer. All reads from `f*`, all
writes to `out[q̄*np + i]`; nothing else is touched.

**Linear (2-point) Bouzidi interpolation** (Bouzidi 2001), for the unknown incoming population at
`x_f` in direction `q̄ = OPP[q]` after streaming, with `f*` = post-collide populations:

```text
q̄-population at x_f, moving-wall term W_q added last (see below):

qd < 1/2 :  f_q̄(x_f) = 2·qd · f*_q(x_f) + (1 − 2·qd) · f*_q(x_ff)               + W_q
qd ≥ 1/2 :  f_q̄(x_f) = (1/(2·qd)) · f*_q(x_f) + ((2·qd − 1)/(2·qd)) · f*_q̄(x_f) + W_q
qd = 1/2 :  → f_q̄(x_f) = f*_q(x_f) + W_q      (≡ today's half-way BB, L280 with W_q = 6 w_q ρ c_q·u_w)
```

- `x_ff = x_f − c_q` is the **second fluid node** away from the wall (the `has_second` flag; if it
  is solid or outside, fall back to the `qd ≥ 1/2` branch, which needs only `x_f` — the standard
  degeneracy handling).
- **Moving-wall term** `W_q = 6 w_q ρ (c_q · u_wall)` (Guo 2002 / Ladd), identical to the existing
  half-way injection (`six * p.wr[q] * rho_row[x] * cu`, L280) — reuse it verbatim, scaled per
  branch per Bouzidi 2001 (the `2·qd` / `1/(2·qd)` weighting also multiplies the wall term in the
  full form; specify the exact per-branch coefficient in the order from the paper).
- **Deviation storage:** `w_q = w_q̄`, so the interpolation is weight-neutral and applies to
  deviation populations unchanged — exactly the argument that makes today's BB deviation-safe
  (kernels.rs L273 "In deviation storage the formula is unchanged").

The pass is a flat loop over `records`; each record does 2–3 reads from `f*` + 1 write to `out`.
Optionally offer the **quadratic (3-point) Bouzidi** as a second operator behind the same record
list (adds a third node `x_fff = x_f − 2c_q`, one more field per record) — but linear first (it is
already 2nd order and is the "simplest post-stream record pass").

---

## 4. STL voxelization pipeline (second builder)

Pure-Rust, additive, plugs into the same record output (COMPETITOR §B6 / Palabos B1 staged-data):

1. **Parse** STL (ascii/binary) via a permissively-licensed crate (`stl_io`, MIT/Apache) → triangle
   soup in lattice coordinates (scaled by the UnitConverter `dx`, SPEC_UNIT_CONVERTER §2).
2. **Inside/outside** per cell: ray casting with odd-even crossing (3-ray majority vote for
   robustness on grazing/edge hits — OpenLB's 4-mode idea, re-derived), or signed distance via a
   BVH over triangles. Writes the `solid` mask (fields.rs L188).
3. **Per-link intersection** `qd`: for each fluid cell adjacent to solid, for each direction `q`
   crossing the surface, ray-vs-triangle (Möller–Trumbore) from `x_f` along `c_q` → nearest hit
   distance `δ` → `qd = δ/(|c_q|·Δx)` (link length `|c_q|` = 1 for axis dirs, √2 for diagonals).
   `qd = 1/2` fallback when the ray misses (grazing) → graceful half-way BB on that link.
4. **Material tagging** (optional): flood/rename to tag inlet/outlet/wall groups (§2).
5. **Emit** the sorted `records` list.

A **BVH** over triangles keeps steps 2–3 near-linear in cells (the aneurysm-class geometry Palabos
runs). Sparse/weighted domains (COMPETITOR §B7) are a separate later item — not required here.

---

## 5. Integration into the existing backends (keep the hot kernel unforked)

- **CpuScalar / CpuSimd:** the Bouzidi pass is a new backend method `apply_bouzidi(f, out, records)`
  run between `stream` and `apply_open_faces`. The fused `CpuSimd` streaming kernel
  (`fused_band`, backend_simd.rs L1037) is **unchanged** — it still writes half-way BB into the
  out-buffer; the Bouzidi pass then overwrites only the recorded links from `f*` (still resident in
  the in-buffer pre-swap). No change to the ring/band machinery. Records are cell-sorted so the pass
  is row-coalesced.
- **WGSL/GPU:** a new small compute entry point `bouzidi` = a scatter kernel over the record buffer
  (one invocation per record), reading `f_in` (post-collide, the push kernel writes `f_out` and
  never `f_in`, wgsl.rs L36–39) and writing `f_out[q̄*n + cell]`. It sits alongside the existing
  `bc`/`clear_probe` post-passes (wgsl.rs L456–547) and is dispatched only when records exist. The
  push `step` kernel is **not** modified — same anti-fork property as SPEC_COLLISION_COMPOSITION:
  the curved BC is data + one generated post-pass, not a per-body shader.
- **Records upload:** one storage buffer of `BouzidiLink` (packed: `cell:u32, q:u32, qd:f32,
  flags:u32, wall_ref:u32`), built host-side once per scenario.

---

## 6. Momentum-exchange force with Bouzidi (drag/lift consistency)

The current probe force is the half-way MEM `ftot = fout + fin + 2 w_q` summed over probed solid
links (kernels.rs L288–291, fused L1237). With interpolated walls the **Galilean-invariant MEM**
(Wen, Zhao & Chen 2014, re-derived — not copied) must use the interpolated incoming/outgoing pairs
so the exchanged momentum is measured at the true wall position; otherwise the force keeps a
staircase bias even though the populations are 2nd order. Specify: the force pass iterates the same
`records`, using `(f*_q(x_f), f_q̄(x_f))` and the link's `qd`, and accumulates
`Σ c_q [f*_q + f_q̄]` with the Bouzidi-consistent correction term. `probed_force()`
(validation_cylinder.rs uses it for Cd/Cl) then reflects the 2nd-order wall. This is a required part
of the order — the T8 acceptance (§7) measures the force, so an uncorrected force would mask the win.

---

## 7. T8 acceptance path — staircase bands → 2nd-order bands

The acceptance anchor is the existing Schäfer-Turek suite, tightened:

1. **Coarse-grid accuracy (the headline):** with Bouzidi on the analytic circle, `t8_2d1_d20`
   (D=20) must land in a **tightened** band around `CD_REF_2D1 = 5.5795` — target **±3%**
   (`Cd ∈ ~5.41..5.75`) versus today's ±~7% staircase band `5.2..6.0`, at the **same resolution**.
   Add a matching tightened `Cl` band. Gate the new band behind a `wall: bouzidi` scenario flag so
   the staircase test stays as the regression baseline.
2. **Convergence order:** a D ∈ {10, 20, 40} sweep shows the drag error `|Cd − 5.5795|` falling at
   **~2nd order** (slope ≈ 2 in log–log) for Bouzidi vs **~1st order** for the staircase — replacing
   the current `err40 < err20` "it converges" assertion (L107–120) with a measured order ≥ ~1.7
   (same style as the TGV order gates).
3. **2D-2 unsteady (Re=100):** `t8_2d2_d40` vortex-shedding `Cl_max` and Strouhal fall in tightened
   bands; the wake is symmetric (no staircase-induced asymmetry).
4. **Degeneracy regression:** a scenario whose body is aligned so every link is `qd = 1/2` produces
   **bit-identical** results to the half-way-BB path (proves graceful degradation and that Bouzidi
   is a strict generalisation).

---

## 8. Invariants & guardrails

- **Additive & gated:** `bouzidi = None` ⇒ today's behaviour exactly; the existing 200+ suite stays
  green unmodified (R5). The staircase T8 tests remain as the baseline; Bouzidi bands are new tests
  behind a scenario flag.
- **Pass order & storage untouched:** deviation `f−w`, q-major SoA `f[q*np+i]`, the one-step pass
  order, and the D2Q9 direction ordering are unchanged; the pass is strictly a post-stream overwrite
  of recorded links (CLAUDE.md core invariants preserved).
- **Backend equivalence:** CpuScalar↔CpuSimd Bouzidi results within the existing tolerance
  (backend_simd_equiv style); CPU↔GPU within T14 (≤1e-5). The scatter over records is order-
  independent per link (each writes a distinct `(q̄, cell)` slot), so it is deterministic (R4).
- **License:** originals only (Bouzidi 2001, Guo 2002, Wen 2014); no GPL/AGPL source. State this in
  the codex order header verbatim, as with prior orders.

---

## 9. Acceptance / adversarial test matrix (for codex)

1. **Poiseuille off-grid wall:** channel walls placed at fractional offsets (`qd ≠ 1/2`, e.g. 0.3,
   0.7) reproduce the analytic parabola to ≤1e-3 — the classic Bouzidi unit test; half-way BB fails
   this at the same resolution.
2. **Cylinder coarse accuracy:** `t8_2d1_d20` + Bouzidi lands in the ±3% Cd band and tightened Cl
   band (§7.1).
3. **Convergence order:** D∈{10,20,40} drag-error slope ≥ ~1.7 (§7.2).
4. **Degeneracy bit-identity:** all-`qd=1/2` body ≡ half-way BB bitwise (§7.4).
5. **Moving wall:** a translating/rotating curved body — the `W_q` term reproduces the imposed
   surface velocity (tangential slip → 0) to 2nd order.
6. **Second-node fallback:** thin gaps where `x_ff` is solid correctly switch to the `qd≥1/2` branch;
   no out-of-bounds read (the `has_second` flag path).
7. **STL round-trip:** a sphere STL voxelised → `qd` records match the analytic-sphere builder to
   ≤1e-3, and sphere drag (Re∈{20,100}) lands in the existing T15 Schiller-Naumann band — proving
   the STL builder and analytic builder agree and that curved 3D BC works.
8. **Force consistency:** the Bouzidi MEM force (§6) equals a control-volume momentum-integral force
   to within the force tolerance (the force is measured at the true wall, not the staircase).
9. **GPU parity:** the `bouzidi` scatter kernel matches the CPU pass within T14 on the cylinder case.
