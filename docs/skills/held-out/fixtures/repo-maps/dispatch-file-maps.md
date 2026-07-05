# Dispatch Fixture File Maps

These synthetic maps are for held-out prompts only. They describe the intended
bundling boundaries without requiring a real repository checkout.

## map-cd-01

Order items:
- `A-2 inlet guard`: `crates/lbm-core/src/guards.rs`
- `A-5 guard error enum`: `crates/lbm-core/src/guards.rs`
- `A-8 scenario schema range`: `crates/lbm-scenario/src/schema.rs`
- `T-A8 adversarial schema test`: `crates/lbm-scenario/tests/schema_ranges.rs`

Expected bundles:
- One implementation order for A-2 + A-5.
- One implementation order for A-8.
- One test order for T-A8.

## map-cd-02

Order items:
- `W-EXT fields extension`: `crates/lbm-core/src/fields.rs`, `crates/lbm-core/src/layout.rs`
- `W-UNIT lattice unit tests`: `crates/lbm-core/tests/d3q27_lattice.rs`
- `W-STRESS stress invariant`: `crates/lbm-core/src/stress.rs`
- `W-ROT rotating wall`: `crates/lbm-core/src/boundary/rotating.rs`
- `W-GRAV well-balanced gravity`: `crates/lbm-core/src/forcing/gravity.rs`
- `W-SCAL scalar ADE`: `crates/lbm-core/src/scalar.rs`
- `W-VOF phase field`: `crates/lbm-core/src/phase_field.rs`
- `W-PART particles`: `crates/lbm-core/src/particles.rs`
- `W-IO output stats`: `crates/lbm-cli/src/output.rs`

Expected compile-bound concurrency ceiling: launch at most 8 simultaneous codex
orders. The ninth order must wait for one of the first eight to complete.

## map-cd-03

Order items:
- `GPU cleanup`: `crates/lbm-core/src/backend/gpu/mod.rs`
- `WASM smoke test`: `crates/lbm-wasm/tests/smoke.rs`

Expected bundles:
- Separate worktrees because paths are disjoint and the second item is a test.
No validation suite or merge is part of dispatch.

