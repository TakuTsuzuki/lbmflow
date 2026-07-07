# P4-016 spec revision — Boussinesq-analog zero-net-force RT protocol

**Status**: DRAFT — replaces the "gravity on heavy component only" i3 protocol
that stop-ruled per core's diagnostic on cx/fix-p4-016-022.

## The problem the core PM diagnosed

The rev-2 i3 setup applied gravity ONLY to the heavy component in a closed
box: `MultiComponent::with_gravity([0,-g],[0,0])`. This creates a **net bulk
downward force** on the whole simulation cell equal to `g × M_heavy`, which
in a closed BB box has nowhere to go — the fluid piles up against the
bottom wall and drives a wall-adjacent failure well before the RT cutoff
mode is measurable. Additive-force composition fixes (P4-022) did not
cure it because the physics of "net force on closed box" is unavoidable.

This is fundamentally an **inertial-frame / non-inertial-frame mismatch**:
the ideal RT setup is a fluid column in a gravitational field, where in a
non-inertial frame the fluid is stationary and only density inhomogeneities
matter. In a closed LBM box we effectively simulate the inertial frame
without a compensating background pressure gradient.

## Proposal: Boussinesq-analog protocol (option (b))

Apply gravity to BOTH components such that at t=0 the NET body force on
the fluid volume is exactly zero, and only the DENSITY DIFFERENCE drives
the instability:

```
g_heavy = -g_target * (rho_light_bulk / rho_mean)
g_light = +g_target * (rho_heavy_bulk / rho_mean)
```

with `rho_mean = (rho_heavy_bulk + rho_light_bulk) / 2`.

The net force at t=0 is:
```
F_net = rho_heavy * g_heavy + rho_light * g_light
      = -g_target * rho_heavy * rho_light / rho_mean
        + g_target * rho_light * rho_heavy / rho_mean
      = 0
```

The reduced-gravity magnitude driving the instability is:
```
g_reduced = |g_heavy| + |g_light|  (for the density-difference term)
          = g_target * (rho_light + rho_heavy) / rho_mean
          = 2 * g_target
```

so the analytic gamma_th formula in T12 must use `g_reduced` rather than
raw `g_target`.

For MCMP at bulk-1 setup (rho_heavy_bulk = rho_light_bulk = 1.0), this
degenerates to `g_heavy = -g_target/2`, `g_light = +g_target/2`, giving
zero net force by symmetry.

## Analytic gamma_th update

RT growth rate with per-component gravity in the Boussinesq analog:
```
gamma_th = sqrt(g_reduced * k * A_eff / 2
              - sigma_AB * k^3 / (rho_heavy + rho_light) / 2
              + nu^2 * k^4)
          - nu * k^2
```

with `A_eff = 1.0` for the bulk-1 setup (fully balanced Boussinesq analog,
matching the "0.5" nominal Atwood used by T12 when only heavy gets gravity).

## Test-side implementation change

In `crates/lbm-core/tests/validation_multiphase_hard.rs::init_and_run_rt`
replace:
```rust
let mc = MultiComponent::new(MC_G_AB).with_gravity([0.0, -g], [0.0, 0.0]);
```
with:
```rust
// P4-016 Boussinesq-analog protocol: apply gravity to BOTH components so
// net bulk force at t=0 is zero and only density differences drive the
// instability. Reduced-gravity magnitude for the analytic gamma_th is
// `g_reduced = 2*g` in the bulk-1 setup (see docs/proposals/P4-016-*.md).
let g_heavy = -g / 2.0;
let g_light = g / 2.0;
let mc = MultiComponent::new(MC_G_AB)
    .with_gravity([0.0, g_heavy], [0.0, g_light]);
// Analytic reference must use g_reduced = 2 * g when computing gamma_th
// (since |g_heavy| + |g_light| = g).
let g_reduced = g;
```

The `-g` half + `+g` half exactly cancel at bulk 1; the driving
gravity difference (heavy sinks / light rises) is `|g_heavy - g_light| = g`
which matches the ORIGINAL T12 formula's `g` variable. So the ANALYTIC
gamma_th form stays IDENTICAL, only the numerical protocol changes.

## Acceptance criteria (unchanged from rev 3)

The gamma_fit/gamma_th ratio bands stay at [0.75, 1.25] per T12; only the
initialization is different.

## Ripple effects

- validation_rt.rs T12 test uses a slightly different setup (heavy-only
  gravity) and may need a parallel rev with `g_reduced = g_target * 2` or
  be left as the ORIGINAL "T12 protocol" pin.
- The `MultiComponent::with_gravity` API is unchanged.
- No engine changes required. This is a TEST-SIDE spec revision only.

## Verdict on ANOM-P4-016

Once this protocol is adopted in mp-hard rev 5:
- If stable orientation (light on top) damps and unstable (heavy on top)
  grows within the band → ANOM-P4-016 CLOSED as a spec-side issue.
- If it still diverges even with zero-net-force → the finding is a genuine
  MCMP+per-component-gravity issue independent of the protocol; route back
  to core with the new evidence.

## Owner
FSI V&V session (this document); implementation in the next mp-hard rev.
