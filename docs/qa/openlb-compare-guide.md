# OpenLB Cross-Solver Comparison Harness

Status: roadmap/operator tool for V&V master plan lane 6.1. The harness does
not run OpenLB or LBMFlow; it compares stored output fields from runs the
operator has already produced.

## Manifest Schema

Required keys:

```json
{
  "benchmark": "cavity",
  "lbmflow_output_dir": "out/lbmflow/cavity",
  "openlb_output_dir": "out/openlb/cavity",
  "grid_size": [129, 129],
  "band_L2rel": 0.02,
  "band_linf": 0.005
}
```

- `benchmark`: one of `cavity`, `cylinder`, or `tgv`.
- `lbmflow_output_dir`: directory containing LBMFlow stored output. Numeric
  VTK/CSV is preferred; gallery PNG is accepted only as an approximate
  luminance fallback.
- `openlb_output_dir`: directory containing OpenLB stored output. Legacy ASCII
  VTK, ASCII `.vti`, and Palabos-style table CSV are supported.
- `grid_size`: common comparison grid. Use `[nx, ny]` or `[nx, ny, nz]`; an
  integer means `[N, N]`.
- `band_L2rel`: maximum relative L2 difference.
- `band_linf`: maximum absolute L-infinity difference.

Optional key:

- `field`: preferred field name, for example `speed`, `ux`, `rho`, or an
  OpenLB DataArray name. If omitted, the harness prefers speed/velocity-like
  fields and then falls back to the first readable field.

## Workflow

1. Build OpenLB in `/Users/taku/projects/cfd-bench` following
   `docs/BENCH_COMPARISON_DRAFT.md`.
2. Build LBMFlow separately and run the same benchmark scenario. Request VTK or
   CSV outputs when possible.
3. Run the OpenLB benchmark separately and store its VTK/CSV output under a
   stable directory.
4. Write the comparison manifest with the output directories and acceptance
   bands.
5. Run:

```bash
python3 scripts/qa/openlb_compare.py compare-openlb-cavity.json
```

The command prints the selected files, source dimensions, common grid, `L2rel`,
and `Linf`. It exits `0` only when both bands pass.

## Utility Modes

Check that the local OpenLB benchmark tree exists without launching OpenLB:

```bash
python3 scripts/qa/openlb_compare.py --check-openlb-build
```

Run the harness self-test on synthetic identical VTK fields:

```bash
python3 scripts/qa/openlb_compare.py --self-test
```

## Honesty Note

This is a roadmap comparison tool, not a completed benchmark result. Actual
OpenLB-vs-LBMFlow comparisons are operator tasks: run both solvers under the
same benchmark definition, preserve the raw output, record hardware and build
settings, then run this harness on the stored fields. PNG comparisons are only
for quick gallery sanity checks; claim-grade comparisons should use numeric
VTK or CSV fields.
