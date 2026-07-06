# SPEC — Checkpoint / restart contract — LANDED

Status: **LANDED**. See `crates/lbm-core/src/solver.rs`
(`Solver::save` / `Solver::load` / `Solver::restore`, `CheckpointError`,
`RankHeader`, on-disk parsing helpers `read_u{8,16,32,64}` / `take` /
`read_rank_file`).

## Design in one sentence

A checkpoint stores only state that cannot be re-derived (the population field
in deviation form, the open-face stale slots, the compact moments, the step
counter), plus hashes of everything that can be re-derived (scenario,
decomposition, masks). On load, geometry is rebuilt from the same scenario and
the hashes must match, so a checkpoint is a small state blob guarded by a
strong "is this the same run?" check.

**Why f32 resume is bit-identical.** Populations are stored in deviation form
(f − w) — the checkpoint payload is the raw deviation bytes; reload is a
byte-for-byte copy back into the buffer with no reconstruction (no re-adding
`w`, which would lose ulps in f32). Codes that store full `f` cannot offer this.

## On-disk layout

`ckpt_<step>/` directory:
- `manifest.json` — rank 0 only, carries `scenario_hash`, `decomp_hash`,
  `dtype`, `lattice`, per-rank {origin, core, file, mask_hash, payload_hash}.
- `rank_<NNNN>.bin` — per-rank state blob, versioned header
  (`magic=LBMKPT`, `format_ver`, `dtype`, `lattice_id`, `q`, `d`, `np`,
  `n_core`, TLV section table) + concatenated sections (`F_PRIMARY`,
  `STALE_STASH`, `MOMENTS`, reserved `RNG` / `PARTICLES` / `STATS`) + payload
  hash. Forward-compatible via TLV: newer sections are skipped by older
  loaders; unknown `format_ver` → clean rejection, never misparse.

## Reserved for later

`RNG`, `PARTICLES`, `STATS` (observer-framework `fieldAverage` accumulators)
slot IDs are defined so the file format is forward-compatible when M-F and
the observer framework land.
