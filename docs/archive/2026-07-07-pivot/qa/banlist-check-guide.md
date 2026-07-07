# Ban-List Check Guide

Lifecycle: living — V&V master plan lane 7.1 grep-based regression guard,
updated in place.

This guide documents `scripts/qa/banlist_check.sh`, the dependency-free Bash
guard for the physics-discipline ban list. The script is intentionally a static
regression guard, not a final audit verdict: every unwhitelisted hit must be
reviewed by a developer who can decide whether it is outside a physics path,
documented physics, or a real violation.

## Purpose

The guard sweeps the physics-sensitive source areas for four grep-able smells:

- case-identity branches, such as `case_id`, `sample_name`, `harshness`, or
  `protocol_name`, where physics could depend on labels instead of physical
  parameters;
- clamps and bounds on transported quantities, using `.clamp(`, `.min(`, and
  `.max(` as broad sentinels;
- bare calibrated literals, using a small guardrail set such as `0.0025`,
  `0.15`, `0.16`, and `2.5e-5`;
- silent physical defaults, using `.unwrap_or(` as the sentinel for physical
  parameters defaulting instead of failing validation.

All four families are reported as HIGH severity. The script exits with status
`1` when any HIGH hit remains unwhitelisted.

## Basic Usage

Run the guard over the default source paths:

```bash
bash scripts/qa/banlist_check.sh
```

Default scan paths are:

```text
crates/lbm-core/src
crates/lbm-cli/examples
```

Run the synthetic positive control:

```bash
bash scripts/qa/banlist_check.sh --self-test
```

Expected verification output:

```text
banlist_check.sh self-test PASS
```

## Whitelist

The whitelist is `scripts/qa/banlist_whitelist.txt`. Each non-comment line is a
fixed source-line pattern. A match is excluded when either the printed
`path:line:source` record or its line-number-independent `path:source` form
contains that pattern.

Prefer `path:source` entries so unrelated line-number movement does not churn
the whitelist:

```text
crates/lbm-core/src/solver.rs:                        let inv_rho = 1.0 / rho[gi].as_f64().max(1.0e-30);
```

Use the whitelist only after review. Good whitelist reasons include:

- the line is an index, allocation, or loop-bound helper rather than physical
  transport;
- the clamp is an explicitly documented numerical underflow guard, such as
  `rho.max(1.0e-30)`;
- the source line has nearby provenance in code or docs, such as a named
  literature correlation or a PHYSICS.md entry;
- the default is CLI ergonomics or output formatting, not a physical parameter.

Do not whitelist a line merely because it is old. If a hit is a real violation,
fix or route it; if it is documented but still physics-affecting, add the
provenance note before whitelisting.

## Output

Each check prints unwhitelisted lines first, then a count:

```text
== CLAMPS ON TRANSPORTED QUANTITIES (HIGH) ==
No unwhitelisted matches.
Count: 0 flagged, 43 whitelisted, 43 raw code matches.
```

The three count fields mean:

- `flagged`: unwhitelisted HIGH matches that make the script fail;
- `whitelisted`: raw matches excluded by `banlist_whitelist.txt`;
- `raw code matches`: grep matches after full-line comments are ignored.

## Positive Control

`--self-test` creates a temporary Rust file containing a synthetic
case-identity branch and a transported-quantity `.max(` guard. It reruns the
scanner against that temporary path with an empty whitelist and requires the
scan to fail. This proves that the detector catches a new banned pattern.

## Limitations

This is a grep guard, so it is deliberately conservative. It cannot prove that
a line changes physics or that a nearby comment is sufficient provenance. It
only prevents new ban-list-shaped code from landing without a conscious review.

When in doubt, apply the physics-discipline stop rule: missing physical input is
a validation error, case-specific branches are not allowed, transport-absorbing
clamps need a documented wall or boundary model, and calibrated constants need
derivation, validity domain, validation, and a PHYSICS.md entry.
