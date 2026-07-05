# SU-S0 Capability Map — Scale-up/Scale-down descriptors (current `main`)

**SCALEUP v1.1 workstream, session SU-S0.** Evidence-graded capability map for the
reactor scale-up descriptors (Np, P, ε̄, U_tip, Re, Fr, N_Q, t_c, torque, resolved
dissipation, tracer θ95, Lagrangian lifelines). No Skill is authored here.

This map **builds on** [`docs/skills/b1-capability-map.md`](../b1-capability-map.md)
(B1) — the general CLI/MCP/schema/artifact/output surface was audited hours ago and is
**not repeated**. SU-S0 audits only the **scaleup-specific deltas**: rotating-impeller
expression, torque→Np, resolved strain/dissipation, scalar/tracer, lifelines,
free-surface/Fr, the Np/P arithmetic layer, and circulation time. Colors follow the
brief: **green** = runnable today with pasted output + no code change; **yellow** =
≤8 eng-h/≤5 files/no new physics or subsystem, with entry-point evidence; **red** =
new subsystem/research, with negative evidence.

---

## 0. Base SHA, environment

| Item | Value |
|------|-------|
| Base SHA (worktree HEAD = branch tip) | `52baafde538043cce41e466642a8fb3a96639c7e` |
| Branch | `skills/su-s0` |
| Toolchain | `rustc 1.93.0 (254b59607 2026-01-19)`, `cargo 1.93.0` |
| Platform | Darwin arm64 (matches B1 §0) |
| Build | `cargo build -p lbm-core --release` → exit 0; `cargo build -p lbm-cli --release` → exit 0; binary `target/release/lbm` |

Definitions used throughout are frozen in `docs/REQ_STIRRED_REACTOR.md` §2.1:
`Re = ρ N D²/μ` (U_tip = πND), `Fr = N²D/g`, `Np = P/(ρ N⁵… )` → `Np = P/(ρ_l N³ D⁵)`,
`P = Ω T_q`, `N = Ω/(2π)`, `N_Q = Q/(N D³)`. (ORDER-v1.1 PM annotation confirms these
match rev.4 §2.1 dimensionally.)

---

## 1. Per-primitive color summary

| # | Scaleup primitive | Color | One-line verdict |
|---|---|---|---|
| 1 | Rotating impeller (volume-penalization interim) | 🔴 Red (user surface) | Native `Solver::set_body_force_field` exists; **NO scenario/CLI surface** expresses it. MRF/IBM absent. |
| 2 | Torque → Np | 🔴 Red | Static-solid momentum-exchange **force** probe IS exposed (3D `force.csv` `step,fx,fy,fz`); but that is a **net linear force, not a torque**. No r×f, no rotating-solid torque. |
| 3 | Resolved strain / dissipation ε=2νS:S | 🟡 Yellow | `gather_strain_rate`/`gather_shear_rate` landed, machine-precision verified — but native-API-only; **no `FieldKind` output channel**. Wiring an output channel = the yellow bound. |
| 4 | Tracer / scalar ADE θ95 | 🔴 Red | No scalar/ADE/tracer field anywhere on `main` (only the `CpuScalar` backend name). |
| 5 | Lagrangian lifelines | 🔴 Red | No particle/tracer integrator; nothing to integrate along. |
| 6 | Fr / free surface | Fr 🟢 Green (arithmetic) / free-surface 🔴 Red | Fr = N²D/g is pure math on **user inputs** (N, D not run artifacts). Free-surface physics absent. |
| 7 | Np/P/V/ε̄/U_tip/Re arithmetic layer | 🟢/🔴 split | Computable from **user inputs alone**: U_tip, Re, Fr, V (geometry). Needs a **run artifact**: torque T_q → therefore Np, P, ε̄ are **RED** (torque never produced). |
| 8 | Circulation time t_c ≈ V/(N_Q N D³) | 🔴 Red | Needs N_Q = Q/(ND³), the integrated impeller discharge — no impeller, no Q. |

**Net:** the scaleup arithmetic layer is a **thin computable shell over a missing core**.
Everything that depends only on *user-supplied* N, D, T, H, g, ν (U_tip, Re, Fr, V, tip-Ma
lattice check) is green arithmetic. Everything that depends on a **torque** or a **discharge
flow** or a **scalar/particle field** is red because the impeller and the scalar/particle
subsystems do not exist on `main`. The one genuine near-term delta is resolved dissipation
(#3), which is yellow only because the physics already computes and the gap is an output
channel.

---

## 2. Green / Yellow — evidence blocks (pasted output)

### 2.1 Static-solid momentum-exchange force IS exposed for 3D scenarios (🟢 — but it is a FORCE, not torque)

The scenario `ProbeSpec::Force` is wired to the native `set_force_probe` on **both** the 2D
compat path and the **3D** build path (`crates/lbm-scenario/src/lib.rs:810-831`,
`s.set_force_probe(move |x, y, z| …)`). A 3D run with an obstacle + a force probe writes
`force.csv`. Real run (24³, sphere obstacle, TRT):

```
$ ./target/release/lbm run su3d.json --out .../su3d_out --json
{ "scenario": "su3d", "status": "completed", "stepsRun": 200, "wallSeconds": 0.276,
  "mlups": 10.0, "diagnostics": { "totalMass": 13224.4, "maxSpeed": 0.0731, "tau": 0.56 },
  "warnings": [], "files": [ "force.csv", "speed_200.vtk", "ux_200.csv" ] }

$ cat .../su3d_out/force.csv
step,fx,fy,fz
50,0.5088738670038795,0.025709954424429776,-4.857e-16
100,-0.5298132955648666,-0.025442445651948220,-1.179e-16
150,-0.8631815816341930,-0.015168118373320577,-4.704e-15
200,0.8851508567678479,-0.005936510305534402,1.339e-15
```

**This is the key positive-for-force / negative-for-torque finding.** `force.csv` is the
integrated momentum-exchange **force** `Σ f` over the *static* solid cells (rims excluded).
It is a **net linear force vector `(fx,fy,fz)`** — there is **no moment arm, no `r×f`, no
per-cell force position** retained, so **no torque `T_q` is computed or emitted**. Since
`Np = P/(ρ_l N³ D⁵)` with `P = Ω T_q` needs `T_q`, Np/P have **no artifact to consume**
(→ §3 primitives 2 & 7 RED). B1 §2.5/§2.6 already established the 3D `force.csv`
`step,fx,fy,fz` header; SU-S0 adds: it is drag-force data, torque-void.

### 2.2 Resolved strain / shear-rate gathers — landed, machine-precision, but native-API-only (🟡)

The tonight-landed gathers are real and correct. `cargo test -p lbm-core --release
--test strain_rate -- --nocapture`:

```
running 3 tests
test strain_rate_gather_is_bit_identical_across_inprocess_split ... ok
couette strain errors: Sxy=1.3877787807814457e-16, gamma consistency=0e0
test couette_strain_rate_matches_half_way_wall_gradient ... ok
poiseuille shear-rate error with Pi_force = -0.5(uF+Fu): 1.528708800171974e-14
test forced_poiseuille_shear_rate_uses_rev4_force_sign ... ok
test result: ok. 3 passed; 0 failed; 0 ignored; ...
```

API surface (`crates/lbm-core/src/solver.rs:1475`,`1497`; distributed facade
`crates/lbm-core/src/dist.rs:980`,`1004`; GPU `crates/lbm-core/src/gpu/solver.rs:337`,`342`):
- `gather_strain_rate() -> Vec<[T;6]>` — S tensor `[S_xx,S_yy,S_zz,S_xy,S_xz,S_yz]`,
  post-streaming/pre-collision stage, Guo `Π_force = −0.5(uF+Fu)` correction, TRT
  `τ_eff = 1/ω_plus`; solid cells zero (FR-STRESS-01 rev.4 compliant).
- `gather_shear_rate() -> Vec<T>` — `γ̇ = √(2 S:S)`.

**What is missing for a user to get an ε field from a scenario run** (the yellow bound):
1. `FieldKind` (`crates/lbm-scenario/src/lib.rs:314`) has exactly
   `Speed, Ux, Uy, Rho, Vorticity` — **no strain/shear/dissipation variant**.
2. The CLI runner (`crates/lbm-cli/src/runner.rs:206-215`, `332-342`) maps `FieldKind`
   only to `ux/uy/rho/speed/vorticity` gathers; it **never calls `gather_strain_rate`
   /`gather_shear_rate`**. (`grep 'strain\|shear\|dissipation' crates/lbm-cli/src`
   returns nothing but the unrelated steady-state `epsilon` tolerance.)

So resolved viscous dissipation `ε = 2ν S:S = ν·γ̇²` (since `γ̇² = 2 S:S`) is **computable
at the Rust API level today** but a scenario run cannot emit it. Wiring is bounded:
add a `FieldKind` variant + a runner arm calling the existing gather + PNG/CSV/VTK write
(reuses the existing scalar-field output path). ≤5 files, ≤8 h, **no new physics** →
YELLOW. (ORDER-v1.1: "order C queued" for exactly this.)

### 2.3 Fr and the input-only arithmetic (🟢 — from USER inputs, not run artifacts)

`Fr = N²D/g`, `Re = ρND²/μ` (U_tip = πND), `U_tip = πND`, `V` (from T, H geometry) are
pure arithmetic on quantities the user supplies. **Critical input-contract note for
S-Fingerprint:** the scenario schema has **no impeller** — `lbm schema | grep -i
'impeller|rpm|angular|diameter'` returns nothing. Therefore `N` (rev/s) and `D` (impeller
diameter) are **USER-declared inputs, never run artifacts**. A scale-up advisor computing
Fr/Re/U_tip must take N, D, g, ρ, μ as prose inputs; the engine neither stores nor emits
them. This is green *as arithmetic on inputs*, with zero dependence on any run.

---

## 3. Red — negative evidence

### 3.1 Rotating impeller — no user surface (🔴)

- Native interim exists: `Solver::set_body_force_field(impl Fn(x,y,z)->[T;3])`
  (`crates/lbm-core/src/solver.rs:995`) + `clear_body_force_field` (`:1019`), exercised by
  `crates/lbm-core/tests/body_force_field.rs`. A per-cell body force *can* emulate a
  penalization impeller in Rust.
- **But there is NO user surface.** `grep -i 'body_force_field|bodyForceField|set_body_force'
  crates/lbm-scenario/src crates/lbm-cli/src` → **empty**. The scenario `physics.force` is a
  single **uniform** `[fx,fy]` gravity vector (`crates/lbm-scenario/src/lib.rs:97`), not a
  per-cell field. No scenario/CLI/MCP path reaches `set_body_force_field`.
- MRF / IBM / sliding-mesh: `grep -i 'rotat|impeller|MRF|immersed|IBM|angular'
  crates/lbm-scenario/src crates/lbm-cli/src` → **empty**. Confirmed RED (matches
  ORDER-v1.1 expectation).

### 3.2 Torque → Np (🔴)

Static-solid force probing is exposed in 3D (§2.1) — but it emits a **linear force**, not a
torque. No `r×f`, no rotating-solid torque, no `T_q`. `Np`/`P` therefore have no producible
input. RED (matches ORDER-v1.1 PM annotation: "momentum-exchange exists for STATIC solids
only; the volume-penalization impeller is a body-force emulation with no torque").

### 3.3 Tracer / scalar ADE θ95 (🔴)

`grep -i 'scalar|tracer|ADE|concentration|species'` across `crates/lbm-scenario/src` and
`crates/lbm-core/src` returns **only** the `CpuScalar` backend type name (a compute backend,
`crates/lbm-core/src/backend.rs:131`) — **no transport scalar field, no advection-diffusion
equation**. `FieldKind` (§2.2) has no concentration channel. No θ95 blend-time is
computable. RED.

### 3.4 Lagrangian lifelines (🔴)

No Lagrangian particle/tracer integrator exists (the particle subsystem is FR-PART, an M-F
requirement not on `main`; `grep` finds no particle advance in scenario/CLI). With no
tracer field (§3.3) and no particle integrator, there is nothing to integrate a lifeline
along. RED.

### 3.5 Free-surface physics (🔴)

`grep -i 'free.surface|freesurface|froude|\bFr\b' crates/lbm-scenario/src crates/lbm-cli/src`
→ **empty**. No free-surface / VOF-height / degassing-top BC on `main`. 3D is single-phase,
`init:rest` only (B1 §4.3). Free-surface = RED. (Fr *arithmetic* is green — §2.3 — but that
is not free-surface physics.)

### 3.6 Circulation time t_c ≈ V/(N_Q N D³) (🔴)

`N_Q = Q/(N D³)` needs `Q`, the integrated impeller discharge across the blade-swept surface.
With no impeller (§3.1) there is no discharge to integrate; `Q` and hence `N_Q` and `t_c`
are not producible. RED.

### 3.7 Np/P/ε̄ from a completed run's artifacts (🔴 — the artifact audit)

What a completed **3D** run actually contains (real run, §2.1):
- `manifest.json`: `scenario, status, stepsRun, wallSeconds, mlups,
  diagnostics{totalMass, maxSpeed, tau}, warnings, files`. **No Np, P, torque, ε.**
- `force.csv`: `step,fx,fy,fz` — net linear force (drag), **no torque** (§2.1).
- `speed_200.vtk`: full 3D structured-points volume (B1 §2.6). **No strain/shear/ε field**
  (§2.2). `ux_200.csv`: z-mid slice.

So from **(b) a completed run's artifacts**, the *only* dynamical scaleup-relevant quantity
is the net drag force — insufficient for Np (needs torque), ε̄ (needs the un-wired dissipation
field or a volume-integrated torque·Ω), or N_Q (needs discharge). The computable-today set is
entirely **(a) user-input arithmetic** (§2.3): U_tip, Re, Fr, V, tip-Ma lattice feasibility.

---

## 4. Top surprises

1. **`force.csv` is torque-void, and it is the whole dynamical surface.** The 3D static-solid
   momentum-exchange probe *is* wired to scenarios (better than one might fear), but it emits
   a **net linear force**, not `r×f`. The single most load-bearing scale-up descriptor (Np,
   via torque) has **no artifact** on `main` — the drag force is the closest thing and it is
   the wrong quantity. This is the pivotal gap for S-Fingerprint.
2. **Resolved dissipation is one output-channel away, not a subsystem away.** `gather_strain_rate`
   /`gather_shear_rate` pass at 1e-16/1e-14. The ε capability is blocked purely by `FieldKind`
   having no variant and the runner not calling the gather — a genuine YELLOW, unusually close
   to green for a "resolved dissipation field" claim.
3. **N and D are inputs, never artifacts.** The schema has no impeller at all, so every
   Fr/Re/U_tip/Np formula consumes user-declared N, D — the engine cannot supply them. The
   arithmetic layer's inputs are *entirely* prose, which fixes S-Fingerprint's input contract:
   it must elicit N, D, T, H, g, ρ, μ from the user, not read them from a run.
4. **The scaleup arithmetic layer is a computable shell over a missing core.** Input-only
   descriptors (U_tip, Re, Fr, V) are green math; every artifact-dependent descriptor (Np, P,
   ε̄, N_Q, t_c, θ95, lifelines) is red because impeller/torque/scalar/particle subsystems are
   absent. The green surface is real but shallow.

---

## 5. Handoff to PM (STOP for review)

- **Green (arithmetic on user inputs, no run needed):** U_tip = πND, Re, Fr = N²D/g, V,
  tip-Ma lattice feasibility check. Reuse B1 §2 for the general run/validate/output surface.
- **Green (run artifact):** net drag force via `ProbeSpec::Force` in 3D (`force.csv`
  `step,fx,fy,fz`) — a *force*, usable for drag/thrust, **not** torque.
- **Yellow (one item):** resolved strain/dissipation ε=2νS:S — physics done and verified;
  bound = wire a `FieldKind` variant + runner arm to the existing gather (≤5 files, ≤8 h,
  no new physics). ORDER-v1.1 "order C".
- **Red (subsystem/research — no Skill):** rotating impeller user surface, torque→Np,
  scalar/tracer ADE θ95, Lagrangian lifelines, free-surface physics, N_Q, t_c.
- **S-Fingerprint input-contract flag:** N and D are USER inputs (no impeller in schema);
  the Np/P/ε̄ path is blocked on the missing torque artifact — an advisor may compute the
  input-only descriptors and must clearly report Np/P/ε̄/N_Q/t_c as not-computable-on-`main`.

Branch tip: `52baafde538043cce41e466642a8fb3a96639c7e` (this commit will move it forward).
