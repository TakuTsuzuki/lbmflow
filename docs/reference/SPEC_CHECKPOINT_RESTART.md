# SPEC — Checkpoint / restart contract (B-5 / C-8, P1 table-stakes)

> Research-agent deliverable for the R-Phase 2 wave (feeds B-5 serialization + C-8 per-rank I/O
> layout). Ready-to-dispatch, same depth/format as the prior three P1 specs. Grounded in firsthand
> reads of `fields.rs::SoaFields` (L168–223), `subdomain.rs::Subdomain`, `dist.rs` (MPI rank =
> part id; the existing cross-rank spec-consistency check at L254 = C-6's hashing), and the pass /
> storage invariants in kernels.rs / backend_simd.rs.
>
> Named in COMPETITOR_ANALYSIS as a table-stakes gap held by all three competitors and not us:
> OpenLB `Serializer`/`save()/load()` (B8), OpenFOAM `startFrom latestTime`, Palabos
> `saveBinaryBlock` (B4). License note: the design is generic serialization — no competitor code.

---

## 1. Design in one sentence

A checkpoint stores **only the state that cannot be re-derived** (the population field, the
open-face stale slots, the compact moments, the step counter, and — reserved — RNG/particle/stat
state), plus **hashes of everything that CAN be re-derived** (scenario, decomposition, masks); on
load, geometry is rebuilt from the same scenario and the hashes must match, so a checkpoint is a
small state blob guarded by a strong "is this the same run?" check.

**f−w is why f32 resume is bit-identical.** Populations are stored in **deviation form (f − w)**
(fields.rs L173; kernels.rs L98–100). The checkpoint payload is the raw deviation bytes; reload is a
byte-for-byte copy back into the buffer. There is **no reconstruction** (no re-adding `w`, which
would lose ulps in f32), so an f32 run resumes **bit-identically**, not just closely — a property
codes that store full `f` cannot offer. This is the load-bearing reason to serialize the deviation
buffer verbatim rather than any "physical `f`" form.

---

## 2. What is serialized vs re-derived

**Serialized (the state blob, per rank):**

| Section | Source | Size | Why it can't be re-derived |
|---|---|---|---|
| `F_PRIMARY` | `SoaFields::f` (deviation, q-major padded) | `Q · np` | the simulation state itself |
| `STALE_STASH` | open-face unknown slots of the ping-pong partner `ftmp` | `O(surface)` | ConvectiveOutflow reads "previous step's post-collide" from the partner buffer (fields.rs L176–178; backend_simd L34–46) — **not derivable**; a resume without it diverges at convective/open faces on the first step |
| `MOMENTS` | compact `rho, ux, uy, uz` | `4 · n_core` | solid-cell densities are **carried, not recomputed** (multiphase wall ρ; backend_simd L1292) — recomputing from `f` would drop them |
| `STEP` | step counter + physical time | 16 B | run position |
| `RNG` *(reserved, M-F)* | stochastic-collision / turbulent-inlet PRNG state | small | FR-IO-06 crash recovery of stochastic runs |
| `PARTICLES` *(reserved, M-F)* | Lagrangian marker / particle buffers | variable | resolved-particle & IBM state |
| `STATS` *(reserved, M-F)* | running time-average / transient-statistics accumulators | variable | observer fieldAverage windows must survive restart |

**Re-derived on load (NOT serialized — validated by hash):** `solid`, `wall_u`, `probe`,
`force_field`, `inlet_profiles`, the Bouzidi record list, the material grid, `LocalGeom`, and the
decomposition — all deterministic functions of the scenario JSON + rank count (fields.rs L188–223;
subdomain.rs). Also **not** serialized: the `CpuSimd` ring/band scratch (`FusedScratch`,
rebuilt lazily, fields.rs L194+).

Rationale: masks/geometry are large and fully determined by the scenario; storing them would bloat
the checkpoint and risk drift. Storing their **hash** instead both shrinks the file and turns "did
the user resume into a different geometry?" into a hard, detectable error (§4).

---

## 3. On-disk layout (C-8: per-rank raw + rank0 manifest)

Mirrors Palabos `saveBinaryBlock` + rank0 metadata and OpenFOAM's time directory; one directory per
checkpoint step:

```
ckpt_<step>/
  manifest.json         # written by rank 0 only
  rank_0000.bin         # rank r's state blob (raw, versioned header + TLV sections)
  rank_0001.bin
  ...
```

Single-rank runs = `nranks == 1`, one `rank_0000.bin`. MPI ranks each write their own file at a
precomputed path (no collective single-file offset write needed for v1 — one file per rank is
simpler and matches the Cartesian decomposition; a collated single-file variant is a later C-item,
COMPETITOR notes we don't need graph-partition I/O).

**Per-rank `.bin` versioned header** (forward-compatible via a TLV section table):

```
magic:        "LBMKPT\0"  (8 bytes)
format_ver:   u32         (bump only on incompatible layout change)
endianness:   u8          (0 = little; v1 rejects cross-endian, byte-swap is a later item)
dtype:        u8          (0 = f32, 1 = f64)
lattice_id:   u16         (D2Q9 / D3Q19 / future D3Q27)
q, d:         u16, u16
np, n_core:   u64, u64    (padded plane length, compact core count — cross-check vs rebuilt geom)
section_count:u32
section_table: [ { id: u32, offset: u64, byte_len: u64 } ; section_count ]   # TLV
payload:      concatenated sections (F_PRIMARY, STALE_STASH, MOMENTS, [RNG|PARTICLES|STATS])
payload_hash: u64         (of all sections; corruption check)
```

**Forward compatibility:** a loader skips any section `id` it does not recognise (the reserved
RNG/PARTICLES/STATS slots are defined now so an older reader gracefully ignores a newer checkpoint's
extra sections, and a newer reader tolerates their absence). A `format_ver` newer than the loader
supports → a clean "checkpoint too new" rejection, never a misparse.

**`manifest.json` (rank 0):**

```json
{
  "kind": "lbmflow-checkpoint",
  "format_version": 1,
  "step": 20000,
  "time": 1.0,
  "dtype": "f32",
  "lattice": "D3Q19",
  "global": [512, 128, 128],
  "scenario_hash": "blake3:…",     // canonicalised scenario JSON
  "decomp_hash":   "blake3:…",     // nranks + per-rank origin/core/neighbors + periodicity
  "nranks": 8,
  "ranks": [
    { "rank": 0, "file": "rank_0000.bin", "origin": [0,0,0], "core": [64,128,128],
      "bytes": 1234567, "payload_hash": "blake3:…", "mask_hash": "blake3:…" },
    …
  ],
  "reserved": { "rng": false, "particles": false, "stats": false }
}
```

---

## 4. Load-time validation (converges with C-6 rank-consistency hashes)

On resume, rebuild geometry/decomp from the (same) scenario, recompute hashes, and validate in
order. Each failure is a **structured, agent-readable diagnostic** (same `diagnostics[]` shape as
SPEC_UNIT_CONVERTER §5 and the MCP tools), with a specific code and a non-zero exit — **never a
partial load**:

| Code | Check | Detects |
|---|---|---|
| `CKPT_TOO_NEW` | `format_version ≤ supported` | a checkpoint from a newer build |
| `CKPT_BAD_MAGIC` / `CKPT_TRUNCATED` | magic + declared section lengths vs file size | corruption / wrong file |
| `CKPT_PAYLOAD_CORRUPT` | recomputed `payload_hash` == header | bit-rot / partial write |
| `CKPT_DTYPE_MISMATCH` | `dtype` == run precision | f32 ckpt into an f64 run (or vice-versa) |
| `CKPT_LATTICE_MISMATCH` | `lattice_id` == run lattice | D2Q9 ckpt into a D3Q19 run |
| `CKPT_SCENARIO_MISMATCH` | `scenario_hash` == rebuilt | user changed physics/BC/geometry |
| `CKPT_DECOMP_MISMATCH` | `decomp_hash` == current | different rank count / partition (repartition-on-load is a later item) |
| `CKPT_MASK_MISMATCH` | per-rank `mask_hash` == rebuilt masks (solid/wall_u/bouzidi/material) | scenario that parses identically but builds different masks — the real safety net; **shares the hashing C-6 already uses to assert spec agreement across ranks** (dist.rs L254) |
| `CKPT_GEOM_MISMATCH` | `np`/`n_core`/`origin`/`core` per rank == rebuilt geom | inconsistent geometry even when hashes are stale |

`scenario_hash` is over the **canonicalised** scenario (stable key order, normalised numbers) so
formatting-only edits don't spuriously reject. The `mask_hash` is the strong check: it catches the
case where the scenario JSON is byte-different-but-semantically-same (allowed) vs
same-looking-but-produces-different-solids (rejected).

---

## 5. API / agent surface

- **Scenario JSON:** a `checkpoint` block — `{ "interval": <steps>, "keepLast": <n>, "dir": "…" }`
  — the run writes `ckpt_<step>/` every `interval` steps, pruning to the last `n` (OpenFOAM
  `purgeWrite`).
- **Resume:** `startFrom` — `"latestTime" | "step:<N>" | 0` (OpenFOAM semantics). CLI
  `lbm run scenario.json --resume ckpt_20000/`. MCP `start_run` gains `resume_from`, so an agent can
  crash-recover or branch a parameter study from a converged state (Pillar 1; COMPETITOR §B2 mapping).
- **Determinism (R4):** save is a pure snapshot; load + resume is bit-exact (deviation bytes copied
  verbatim, §1). Checkpointing does not perturb the trajectory of a non-checkpointed run (writing a
  snapshot must not touch live buffers).

---

## 6. Invariants & guardrails

- **Additive:** no `checkpoint` block ⇒ today's behaviour; the 200+ suite stays green unmodified (R5).
- **Storage/pass invariants untouched:** deviation `f−w`, q-major SoA `f[q·np+i]`, one-step pass
  order, D2Q9 ordering (CLAUDE.md) are unaffected — the checkpoint reads/writes existing buffers at
  a **step boundary** (after moments, before the next collide), a quiescent point where `f` is the
  primary and fully consistent.
- **MPI (C-8):** each rank snapshots its own `Subdomain` (subdomain.rs) independently; no halo
  exchange at save. Rank 0 gathers only the small manifest metadata. The bit-identical-to-monolithic
  property (dist.rs L16) means a resumed MPI run must match a resumed single-rank run over the same
  covered region — an extension of the T13 partition-invariance gate.
- **Reserved sections are defined now** (RNG/PARTICLES/STATS) precisely so M-F's FR-IO-06 crash
  recovery slots in without a format-version bump (§3 forward-compat).

---

## 7. Acceptance / adversarial test matrix (for codex)

1. **Bit-identical resume (headline, f32 AND f64):** run `N+M` continuous; separately run `N` →
   save → load → `M`; assert `f`, `rho`, `ux/uy/uz` **bitwise identical** at step `N+M`. This is the
   whole contract; it must hold for **f32** (proving the deviation-byte-copy claim, §1) as well as f64.
2. **Convective-outflow resume:** a scenario with a `Convective` face — resume is bit-identical only
   if `STALE_STASH` is captured; a deliberately stash-omitted build must **fail** this test (proves
   the stash is load-bearing, §2).
3. **Multiphase solid-density resume:** a Shan-Chen case with wall densities — resume preserves solid
   ρ bitwise (proves `MOMENTS` is serialized, not recomputed, §2).
4. **Rejection matrix:** every row of §4 fires its specific code with non-zero exit and no partial
   state mutation — truncated file, one flipped payload byte, changed `nu`/geometry, changed nrank,
   f32↔f64, D2Q9↔D3Q19, and a hand-bumped `format_ver` → `CKPT_TOO_NEW`.
5. **MPI resume ≡ monolithic resume:** an 8-rank save→load→resume matches a single-rank save→load→
   resume over the covered region, bitwise per cell (T13-style).
6. **Forward-compat skip:** a checkpoint carrying an unknown reserved section id loads successfully
   with that section skipped (proves the TLV section table, §3).
7. **Non-perturbation:** a run with `checkpoint.interval = k` produces a trajectory bit-identical to
   the same run with checkpointing disabled (writing snapshots must not disturb the simulation).
8. **keepLast pruning:** with `keepLast = 2`, only the last two `ckpt_*` dirs survive; `latestTime`
   resolves to the newest.
