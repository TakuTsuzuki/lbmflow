# LBMFlow UQ Sweep Summary

- Runs: 6
- Parameters: `param.physics.nu`, `param.grid.nx`, `param.grid.ny`
- Numeric QOIs: 32
- Bootstrap CI rows with repeats: 0

## Top Normalized OAT Sensitivities

| Parameter | QOI | Normalized slope | Raw slope | Groups |
|---|---:|---:|---:|---:|
| `param.physics.nu` | `qoi.manifest.warningCount` | -2.33333 | -33.3333 | 2 |
| `param.physics.nu` | `qoi.probe.point_16_16.uy.min` | 1.42261 | 0.00144462 | 2 |
| `param.physics.nu` | `qoi.probe.point_16_16.ux.max` | -1.34193 | -0.0601298 | 2 |
| `param.physics.nu` | `qoi.probe.point_16_16.rho.std` | 0.74939 | 0.00958586 | 2 |
| `param.physics.nu` | `qoi.probe.point_16_16.ux.last` | -0.613522 | -0.200087 | 2 |
| `param.physics.nu` | `qoi.probe.point_16_16.ux.mean` | -0.610879 | -0.175804 | 2 |
| `param.physics.nu` | `qoi.probe.point_16_16.ux.min` | -0.519786 | -0.212153 | 2 |
| `param.physics.nu` | `qoi.probe.point_16_16.ux.std` | 0.509465 | 0.0497985 | 2 |
| `param.physics.nu` | `qoi.probe.point_16_16.uy.mean` | 0.484375 | 0.041703 | 2 |
| `param.physics.nu` | `qoi.field.speed.kineticEnergy` | 0.346036 | 5.47513 | 2 |
| `param.physics.nu` | `qoi.probe.point_16_16.uy.std` | 0.337616 | 0.0128843 | 2 |
| `param.physics.nu` | `qoi.probe.point_16_16.uy.last` | 0.306539 | 0.0383793 | 2 |
| `param.physics.nu` | `qoi.field.speed.sum` | 0.296308 | 277.799 | 2 |
| `param.physics.nu` | `qoi.field.speed.mean` | 0.296308 | 0.154891 | 2 |
| `param.physics.nu` | `qoi.probe.point_16_16.uy.max` | 0.287475 | 0.0365955 | 2 |
| `param.physics.nu` | `qoi.field.speed.rms` | 0.175302 | 0.152272 | 2 |
| `param.physics.nu` | `qoi.manifest.tau` | 0.122807 | 3 | 2 |
| `param.physics.nu` | `qoi.field.speed.std` | 0.106948 | 0.073911 | 2 |
| `param.physics.nu` | `qoi.manifest.maxSpeed` | 0.0399765 | 0.15338 | 2 |
| `param.physics.nu` | `qoi.field.speed.max` | 0.0399765 | 0.15338 | 2 |

## Bootstrap Confidence Intervals

Computed only for exact repeated parameter points with `n >= 2`. See `bootstrap_ci.csv` for the full table.
