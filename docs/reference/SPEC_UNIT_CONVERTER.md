# SPEC — SI UnitConverter: echo-back & stability-diagnostic contract

> Research-agent deliverable for the R-Phase 2 wave order (W-UNIT node, REQ §11 DAG; closes
> the B1 audit #1 red and the user's stated #1 gap). Ready-to-dispatch spec section.
> Derived from firsthand reads of OpenLB `src/core/unitConverter.h` (3 named constructors,
> 9 conversion factors, `print()` stability warning) and Palabos `src/core/units.h`
> (`IncomprFlowParam<T>`). Formulas below are the **single source of truth** — codex
> adversarial tests are written from §2/§4/§7, implementation from §1–§6.
>
> **Core invariant preserved:** conversion happens ONLY at the scenario/validation boundary
> (in `lbm-scenario`). The core (`lbm-core`) stays lattice-unit-only. `f−w` deviation storage
> is unaffected. No physical units cross the public core API.

---

## 1. Input contract

The user supplies **physical (SI) anchors** + **two numerical knobs**. Three named constructors
match how engineers actually think (mirrors OpenLB); pick exactly one form per scenario:

| Constructor | User provides | The other knob is derived |
|---|---|---|
| `FromResolutionAndRelaxationTime` | `N`, `tau` | `u_lat` derived |
| `FromResolutionAndLatticeVelocity` | `N`, `u_lat` | `tau` derived |
| `FromRelaxationTimeAndLatticeVelocity` | `tau`, `u_lat` | `N` derived |

**Physical anchors (SI, all three forms require):**

| Symbol | JSON key | Unit | Meaning |
|---|---|---|---|
| `L_phys` | `characteristicLength` | m | characteristic length (defines Re and dx) |
| `U_phys` | `characteristicVelocity` | m/s | characteristic velocity (defines Re and dt) |
| `ν_phys` | `kinematicViscosity` | m²/s | kinematic viscosity |
| `ρ_phys` | `density` | kg/m³ | reference density (**required** for pressure/force conversion; no default — force explicit to avoid the silent-water-vs-air error) |

**Numerical knobs:**
- `N` = `resolution` — cells spanning `L_phys` (dimensionless integer). Sets `dx = L_phys/N`.
- `u_lat` = `latticeVelocity` — characteristic velocity in lattice units (dimensionless). The
  **Mach/compressibility knob**.
- `tau` = `relaxationTime` — BGK relaxation time (dimensionless).

**Domain extents** `lx, ly, lz` (m) map to cell counts `nx = round(lx / dx)` etc. (document the
`+1` rim convention already used by the wall-rim builder; off-lattice adds one more — see
SPEC_UNIT_CONVERTER consumers). `Re := U_phys·L_phys/ν_phys` is computed, never user-supplied.

**Optional physical inputs that also get converted** (present-if-declared):
- `g_phys` gravity / body acceleration (m/s²) → lattice via `C_accel`.
- `endTime` (s) or `endStepCount` → total lattice steps via `dt`.
- `p_ref` reference pressure (Pa), default 0.

---

## 2. Derivation formulas (canonical — cs² = 1/3, tau = 3ν+0.5)

Let `Re = U_phys·L_phys / ν_phys`.

```
dx            = L_phys / N                         # m per cell   (= C_length)
dt            = u_lat · dx / U_phys                # s per step
              = u_lat · L_phys / (N · U_phys)
ν_lat         = u_lat · N / Re                     # = ν_phys · dt / dx²   (identity)
tau           = 3 · ν_lat + 0.5                    # cs² = 1/3
omega         = 1 / tau
```

The three constructors invert the same relations (all consistent with the identity
`tau = 3·u_lat·N/Re + 0.5`):

```
FromResolutionAndLatticeVelocity (N, u_lat):  ν_lat = u_lat·N/Re ;  tau = 3ν_lat+0.5
FromResolutionAndRelaxationTime  (N, tau):    ν_lat = (tau−0.5)/3 ;  u_lat = ν_lat·Re/N
FromRelaxationTimeAndLatticeVelocity (tau,u_lat): ν_lat=(tau−0.5)/3 ; N = ν_lat·Re/u_lat  (round; warn on rounding drift)
```

**Conversion factors** (physical = factor × lattice), echoed for every quantity:

```
C_length      = dx                                 # m
C_time        = dt                                 # s
C_velocity    = dx / dt   (= U_phys / u_lat)       # m/s
C_viscosity   = dx² / dt                            # m²/s
C_density     = ρ_phys                              # kg/m³   (lattice ρ≈1 ↔ ρ_phys)
C_pressure    = ρ_phys · (dx/dt)²  = ρ_phys·C_velocity²   # Pa    (p_phys = C_pressure·(p_lat−p_ref_lat), p_lat = ρ_lat/3)
C_force       = ρ_phys · dx⁴ / dt²                  # N
C_accel       = dx / dt²                            # m/s²    (g_lat = g_phys / C_accel)
```

**Dimensionless outputs:**

```
Ma            = √3 · u_lat            # lattice Mach number (u_lat / cs, cs = 1/√3)
Re            = U_phys·L_phys/ν_phys  # recomputed for the manifest (sanity echo)
Re_grid       = u_lat / ν_lat  = Re / N     # per-cell (grid) Reynolds number
Kn            = Ma / Re               # informational (Knudsen ~ Ma/Re)
```

---

## 3. Echo-back: derived quantities to emit

Every resolved scenario emits a `units` block (§5) containing, at minimum:
`dx`, `dt`, `ν_lat`, `tau`, `omega`, `u_lat`, all 8 conversion factors (§2),
and the four dimensionless numbers `Re`, `Ma`, `Re_grid`, `Kn`. If `endTime`/`g_phys`/`p_ref`
were supplied, also echo `total_steps`, `g_lat`, `p_ref_lat`. This is the OpenLB `print()` table
made machine-readable and agent-consumable (Pillar 1).

---

## 4. Stability-diagnostic contract

Every diagnostic ties a **derived quantity** to a **knob the user can change**, so the remedy is
actionable. Thresholds are anchored to LBMFlow's existing validator (tune-stability Skill:
tau≥~0.55, |u_lat|≤0.15 hard 0.3, grid-Re ≤15) so the converter introduces **no contradictory
thresholds**. Severity ∈ {`error` (reject), `warn` (run, flag)}.

| id | Quantity | Condition | Severity | Formula / threshold | Remedy (echoed) |
|---|---|---|---|---|---|
| `TAU_UNSTABLE` | `tau` | `tau ≤ 0.5` | **error** | hard floor cs²>0 requires τ>0.5 | raise `N` or `u_lat` (both raise ν_lat); or lower Re |
| `TAU_LOW` | `tau` | `0.5 < tau < 0.55` | warn | over-relaxation → BGK instability near walls | raise `N` (preferred) or `u_lat`; consider TRT (magic Λ=3/16 stabilises) |
| `TAU_HIGH` | `tau` | `tau > 2.0` | warn | over-diffusive / cells wasted (accuracy fine, cost high) | lower `N` or `u_lat` to cut step count |
| `MACH_HARD` | `u_lat` (`Ma`) | `u_lat > 0.3` (`Ma > 0.52`) | **error** | compressibility error dominates; likely divergence | lower `u_lat`; compensate τ by raising `N` |
| `MACH_HIGH` | `u_lat` (`Ma`) | `0.15 < u_lat ≤ 0.3` (`0.26 < Ma ≤ 0.52`) | warn | non-negligible O(Ma²) compressibility error | lower `u_lat` toward ≤0.15; raise `N` to hold τ |
| `GRID_RE_HIGH` | `Re_grid` | `Re_grid > 15` | warn | under-resolved cell → dispersion/instability | raise `N` so `N ≥ Re/15` (since Re_grid = Re/N) |
| `RESOLUTION_ROUNDING` | `N` | `FromRelaxationTimeAndLatticeVelocity` gives non-integer N | warn | derived N rounded; actual τ or Re drifts from request | switch to an `N`-fixed constructor, or accept the reported drift |
| `DENSITY_MISSING` | `ρ_phys` | pressure/force output requested but `density` absent | **error** | conversion undefined | set `density` explicitly (e.g. 998.2 water, 1.204 air @20°C) |

**Feasible-window note (emit when any warn fires):** for fixed `Re`, the three knobs satisfy
`tau = 3·u_lat·N/Re + 0.5`, `Ma = √3·u_lat`, `Re_grid = Re/N`. The stable box is
`u_lat ≤ 0.15` **and** `N ≥ Re/15` **and** `tau ≥ 0.55`. The converter SHOULD emit one suggested
`(N, u_lat)` pair that lands inside the box (e.g. `N* = ceil(Re/15)`, then `u_lat* =
min(0.15, (0.55−0.5)/3 · Re/N*)`), so an agent can auto-correct without solving the algebra.

**Verdict** = `error` if any error fires, else `warn` if any warn fires, else `ok`. `lbm validate`
returns non-zero exit on `error`.

---

## 5. Machine-readable output schema (manifest.json + `lbm validate`)

```json
{
  "units": {
    "constructor": "FromResolutionAndLatticeVelocity",
    "inputs": { "characteristicLength": 0.1, "characteristicVelocity": 1.0,
                "kinematicViscosity": 1.0e-6, "density": 998.2,
                "resolution": 200, "latticeVelocity": 0.1 },
    "lattice": { "dx": 5.0e-4, "dt": 5.0e-5, "nu_lattice": 0.002,
                 "tau": 0.506, "omega": 1.976, "u_char_lattice": 0.1 },
    "conversion_factors": { "length_m": 5.0e-4, "time_s": 5.0e-5,
                 "velocity_m_s": 10.0, "viscosity_m2_s": 5.0e-3,
                 "density_kg_m3": 998.2, "pressure_Pa": 99820.0,
                 "force_N": 0.01247, "acceleration_m_s2": 200000.0 },
    "dimensionless": { "reynolds": 100000.0, "mach": 0.173,
                 "grid_reynolds": 50.0, "knudsen": 1.73e-6 },
    "verdict": "warn",
    "diagnostics": [
      { "id": "TAU_LOW", "severity": "warn", "quantity": "tau", "value": 0.506,
        "threshold": 0.55, "message": "tau below 0.55: BGK over-relaxation risk near walls.",
        "remedy": "raise resolution N or latticeVelocity; or use TRT collision." },
      { "id": "GRID_RE_HIGH", "severity": "warn", "quantity": "grid_reynolds", "value": 50.0,
        "threshold": 15.0, "message": "grid Reynolds > 15: cell under-resolved.",
        "remedy": "raise N so N >= Re/15 (>= 6667 here)." }
    ],
    "suggestion": { "resolution": 6667, "latticeVelocity": 0.15 }
  }
}
```

The `diagnostics[]` shape is the same one the MCP `run_status` / `validate` tools already surface,
so an agent gets structured, machine-readable stability guidance (Pillar 1: explainable failures).

---

## 6. Implementation & invariant notes

- Lives in `lbm-scenario` as a `FlowParams` struct + resolver; runs during `lbm validate` and at
  scenario load. Output = raw lattice `{tau (or omega), u_lat, nx/ny/nz, g_lat, total_steps}` fed to
  the core. **The core never sees SI.**
- The resolver is pure/deterministic (no floating-point env dependence beyond IEEE-754); same inputs
  → byte-identical derived values across backends (supports determinism guarantee / R4).
- Reuses the existing tune-stability thresholds verbatim — do not fork them. If a threshold changes,
  it changes in one shared constant consumed by both the validator and this converter.
- `f−w` deviation storage: unaffected — conversion is upstream of any distribution allocation.
- Backward compatibility: a scenario with no `units` block is interpreted as raw lattice units (today's
  behaviour), so the existing 200+ tests stay green without modification (R5).

---

## 7. Acceptance / adversarial test matrix (for codex, written from §2/§4)

1. **Round-trip identity:** for random valid `(L,U,ν,ρ,N,u_lat)`, `getPhysX(getLatticeX(q)) == q`
   for q ∈ {velocity, viscosity, pressure, force, accel} to ≤1e-12 relative.
2. **Constructor equivalence:** the three constructors, seeded to describe the same physical case,
   produce identical `(dx, dt, ν_lat, tau, u_lat, Re, Ma, Re_grid)` (≤1e-12).
3. **Known-case anchors:** Schäfer-Turek 2D-2 (Re=100) and a Poiseuille case — assert `tau`, `Ma`,
   `Re_grid` match hand-computed reference values; assert the momentum-exchange drag, once run,
   lands in the existing validation band (cross-check the converter didn't shift physics).
4. **Threshold boundaries:** cases placed exactly at τ=0.5, 0.55, 2.0; u_lat=0.15, 0.3; Re_grid=15
   fire exactly the expected diagnostic id + severity (off-by-epsilon on both sides).
5. **Density-missing:** requesting pressure/force output without `density` → `DENSITY_MISSING` error,
   non-zero exit.
6. **Rounding drift:** `FromRelaxationTimeAndLatticeVelocity` with inputs forcing non-integer N →
   `RESOLUTION_ROUNDING` warn + the reported actual τ/Re drift matches the rounded N.
7. **Suggestion validity:** for any case that fires a warn, the emitted `suggestion` `(N,u_lat)`,
   fed back through the converter, yields `verdict == "ok"`.
8. **Invariant guard:** a scenario without a `units` block still loads and runs bit-identically to
   the pre-converter baseline (R5 regression).
