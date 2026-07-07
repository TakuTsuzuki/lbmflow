# Grid Convergence Index Usage

This note standardizes the Roache Grid Convergence Index (GCI) vocabulary used
by validation tests.

## Vocabulary

- `f_h` is the fine-grid quantity of interest.
- `f_2h` is the medium-grid value.
- `f_4h` is the coarse-grid value.
- `r_21 = h_medium / h_fine`; for a doubled resolution this is `2`.
- `r_32 = h_coarse / h_medium`.
- `p` is the observed order from the three-grid series.
- `GCI_21` is the fine-grid numerical uncertainty estimate.
- `GCI_32` is the medium-grid estimate.
- `asymptotic_ratio = GCI_32 / (r_21^p * GCI_21)`. Values near `1` indicate
  that the three grids are behaving as an asymptotic refinement family.

The helper reports GCI values as percentages in `GciResult`, but the primitive
`gci(...)` returns the fractional form.

## Adding GCI to a Convergence Test

Put the three quantity-of-interest values in coarse-to-fine order:

```rust
mod common;

use common::gci::gci_from_series;

let gci = gci_from_series([qoi_d20, qoi_d40, qoi_d80], 2.0);
println!(
    "GCI observed_order={:.4}; Richardson_limit={:.8}; GCI_21_pct={:.6}; asymptotic_range_ratio={:.6}",
    gci.observed_order,
    gci.richardson_limit,
    gci.gci_21_pct,
    gci.asymptotic_ratio
);
assert!((0.9..=1.1).contains(&gci.asymptotic_ratio));
```

For unequal ratios, call `gci_result(f_coarse, f_medium, f_fine, r_21, r_32,
1.25)`. Use the `1.25` safety factor only when the three-grid asymptotic-range
check passes. For two-grid estimates or non-asymptotic series, Roache's
conservative safety factor is `3.0`, and the result should be labelled as a
screen rather than a confirmed solution-verification uncertainty.
