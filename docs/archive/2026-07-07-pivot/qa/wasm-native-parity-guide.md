# WASM/native parity guide

Lifecycle: living QA operator guide. Created 2026-07-07 for V&V master-plan
lane 8.2.

This lane compares the browser/WASM field output against the native Rust
reference for the same canonical scenario. The committed native test only
exports and re-imports the native snapshot. Building WASM and capturing the
WASM snapshot is an operator task, often run in a second terminal or on a
second machine.

## Canonical scenario

- Case: `VALIDATION.md` T2 body-force Poiseuille.
- Grid: D2Q9, `32 x 32`.
- Boundary conditions: periodic left/right, bounce-back bottom/top.
- Collision: TRT, magic `3/16`.
- Parameters: `nu = 0.1`, `tau = 0.8`, force `[1.0e-6, 0.0]`.
- Export path: `test/tmp/native_scenario_field.json`.
- Field order: row-major `y * nx + x`, full domain including the solid rim.
- Default comparison mask: fluid cells only, using `solid[]` from the native
  snapshot.

## Native snapshot

From the repository root:

```bash
cargo test -p lbm-core --release --test wasm_native_parity_export -- --nocapture
```

The test writes:

```text
test/tmp/native_scenario_field.json
```

It also validates the export/import round trip on the native side.

## WASM snapshot

Build the web engine in the normal operator environment. From the repository
root, the current project command is:

```bash
wasm-pack build crates/lbm-wasm --target web --release --out-dir ../../web/src/engine/pkg
```

Then run the WASM path with the same T2 scenario and export a JSON snapshot
with the same top-level shape:

```json
{
  "schema_version": 1,
  "scenario_id": "t2_poiseuille_periodic_x_bounceback_y_n32",
  "nx": 32,
  "ny": 32,
  "rho": [1.0],
  "ux": [0.0],
  "uy": [0.0],
  "solid": [true]
}
```

The arrays above are illustrative; real arrays must contain `nx * ny` values.
The comparator also accepts fields under a top-level `fields` object, for
example `{"fields": {"rho": [...], "ux": [...], "uy": [...]}}`.

## Compare

```bash
python3 scripts/qa/wasm_parity_check.py \
  test/tmp/native_scenario_field.json \
  test/tmp/wasm_scenario_field.json
```

Default assertions:

- `rho` L2rel <= `1e-4`.
- combined `velocity` L2rel <= `1e-4`.
- `scenario_id`, `nx`, `ny`, and solid masks must match when present.

Useful options:

```bash
python3 scripts/qa/wasm_parity_check.py --band 5e-5 --fields rho,velocity native.json wasm.json
python3 scripts/qa/wasm_parity_check.py --fields rho,ux,uy --include-solid native.json wasm.json
python3 scripts/qa/wasm_parity_check.py --self-test
```

## Difference budget

The default `1e-4` band is for f32 WASM fields compared with the native f64
snapshot. A passing nonzero difference is expected; it comes from f32
arithmetic in the WASM/browser path, tau/omega precision, operation ordering,
and macroscopic-field reconstruction from LBMFlow's deviation-storage
population representation. A failure is a routing signal for lane 8.2: first
confirm the snapshots came from the same scenario and step count, then inspect
whether the discrepancy is precision/order-only or a behavioral drift in the
WASM path.
