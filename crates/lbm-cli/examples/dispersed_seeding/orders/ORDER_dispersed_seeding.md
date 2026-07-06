# ORDER — `examples/dispersed_seeding` (3D dispersed-phase withdraw→eject→deposit demo)

## MISSION FRAMING (read first)
This order has **two deliverables of equal weight**:
1. A **small, self-contained, runnable 3D demo** (see spec below).
2. **`SPEC_FINDINGS.md`** — an adversarial report of every spec bug, ambiguity,
   physical inconsistency, unit error, or intractable step you discover *while
   implementing and running experiments*. Build the SMALLEST thing that runs
   end-to-end, run the two sample inputs as EXPERIMENTS, and let the experiments
   expose where this spec is wrong or underspecified. If a step proves
   intractable within a small implementation, implement the reduced form, RUN
   it, and document the gap + a proposed fix in `SPEC_FINDINGS.md` — do NOT
   silently drop it. Completeness of the spec-hardening report matters as much
   as the code compiling.

## 0. Intent & isolation (hard constraints)
Demonstrate: aspirate a dispersed rigid-particle phase from a **tall 3D
reservoir**, eject it as an impinging jet into a **shallow 3D target tray**, let
particles deposit on the tray floor, and score the deposited number-density
distribution over an M×N partition.

- **Do NOT modify any file under `crates/lbm-core/src`, `crates/lbm-scenario`,
  `crates/lbm-wasm`, or `web/`.** All new code lives under the example
  directory. The only edit allowed outside it is adding an `[[example]]` stanza
  to `crates/lbm-cli/Cargo.toml` (reuse `serde`/`serde_json` if already present).
- Use the **public** `lbm_core` / `lbm_core::compat` API only (create a **3D
  D3Q19** sim, set velocity inlet(s) + a suction outlet / open boundary,
  `step()`, read physical velocity `ux()/uy()/uz()` per cell). Discover exact
  signatures from the crate; never reach into private modules. If the public API
  cannot express something you need (e.g. a suction outlet BC), record it in
  `SPEC_FINDINGS.md` and use the closest available construct.
- **Domain-neutral vocabulary everywhere** (identifiers, comments, strings,
  output labels): "reservoir", "dispersed particles", "target tray",
  "deposition", "partition/bin". No application-domain terms.
- Language policy: all artifacts in English.

## 1. Files to create
```
crates/lbm-cli/examples/dispersed_seeding/
  main.rs          # CLI entry: parse protocol JSON, orchestrate phases, write outputs
  protocol.rs      # serde structs for the input schema + SI->lattice nondimensionalization
  reservoir.rs     # 3D reservoir phase: suction-driven withdraw, entrained-particle sampling
  particles.rs     # Lagrangian one-way integrator (Stokes drag + gravity + agitation + deposition)
  readout.rs       # M x N bin density, CV, Max/Mean, CSV + 3D volume export
  README.md        # physics, model assumptions, non-goals, how to run, how to read outputs
  SPEC_FINDINGS.md # the spec-hardening report (deliverable #2)
  sample_gentle.json
  sample_harsh.json
```

## 2. Physics specification (implement; deviate only with a SPEC_FINDINGS entry)

**World:** full 3D, `x,y` horizontal, `z` vertical (gravity −z). SI inputs,
nondimensionalized internally. **No free surface** — every domain stays
fluid-filled (this is a stated MVP simplification, not a bug; note it in README).

**Phase A — reservoir (3D, D3Q19):** a tall box `res_w × res_w × res_h`
(res_h ≫ res_w), fluid-filled, dispersed particles pre-settled inside. A
**suction outlet** at a chosen depth `depth_frac` drives an LBM flow; fluid mass
is replenished through an open top boundary (constant pressure) so the domain
stays filled (NO free surface). Particles are advected one-way in the live LBM
field (see §Particles); a particle that exits through the suction outlet is
**extracted** → it joins the batch handed to Phase B, carrying its diameter.
This makes "how much / what size distribution is extracted at a given depth"
a measured quantity (`n_extracted`, extracted-diameter histogram).
*If a suction-outlet BC is not expressible via the public API, fall back to
concentration-weighted sampling of pre-settled particles at `depth_frac`, RUN
it, and log the substitution in SPEC_FINDINGS.md.*

**Phase B — target tray (3D, D3Q19):** wide shallow box `tray_w × tray_d ×
tray_h` (tray_h small), fluid-filled, side+floor walls, top open (outlet). One
or more **downward velocity inlets** (jet points) at the top inject the
extracted batch. Particles advect one-way + settle; on crossing the floor plane
a particle is frozen and counted into bin `(i,j)`, `i=floor(x/(tray_w/M))`,
`j=floor(y/(tray_d/N))`.

**Dispersed particles (Lagrangian, in-example, one-way):** diameter `d_p` sampled
lognormal(mean `d_p`, cv `d_p_cv`); mass `m_p = ρ_p·(π/6)·d_p³`.
- EoM in container frame:
  `m_p dv/dt = 3π μ d_p (u_f − v) + m_p (1 − ρ_f/ρ_p) g ẑ + F_agit`, `μ = ρ_f·nu`.
- Settling velocity `v_s = (ρ_p−ρ_f) g d_p²/(18μ)`; relaxation `τ_p = ρ_p d_p²/(18μ)`;
  Stokes number `St = τ_p·U_jet/L`.
- Integrate drag **semi-implicitly** (exponential update) for stability at small τ_p.
- Walls/floor: elastic clamp, no penetration; simple exclusion only.

**`agitate` (translational):** during the window apply `F_agit = −m_p·a_c(t)`,
`a_c(t)=A ω² sin(ωt) x̂`; `count`→window = `count·(2π/ω)`; `speed`=peak wall
velocity `A·ω`; `amplitude`=`A`. Log `Fr = A ω²/g`.

**Nondimensionalization:** pick `dx` from grid, choose `dt` so max lattice
velocity ≤ 0.1 (Mach guard); `τ = 3·nu* + 0.5`. **Log**: `dx, dt, nu*, u_jet*,
Ma, Re_jet=u_jet·d_nozzle/nu, St, Fr, τ`. Abort with a clear message if
`Ma>0.3` or `τ<0.51`.

## 3. Protocol input schema (SI units)
```jsonc
{
  "grid":       { "res_nx": 32, "res_ny": 32, "res_nz": 128,
                  "tray_nx": 96, "tray_ny": 96, "tray_nz": 24, "dx_m": 5e-4 },
  "fluid":      { "nu_m2s": 1.0e-6, "rho_f_kgm3": 1000.0 },
  "particles":  { "rho_p_kgm3": 1050.0, "d_p_m": 100e-6, "d_p_cv": 0.1, "count": 20000, "seed": 1 },
  "reservoir":  { "height_m": 0.10, "width_m": 0.016, "fill_height_m": 0.06, "initial_conc": 1.0 },
  "target":     { "width_m": 0.048, "depth_m": 0.048, "height_m": 0.012, "partitions_x": 12, "partitions_y": 12 },
  "protocol": [
    { "op": "settle",   "duration_s": 30.0 },
    { "op": "withdraw", "volume_frac": 0.5, "rate_uLs": 300.0, "depth_frac": 0.8 },
    { "op": "eject",    "points_xy_frac": [[0.5,0.5]], "rate_uLs": 179.0, "height_m": 0.008 },
    { "op": "agitate",  "pattern": "translational", "count": 0, "speed_mms": 90.0, "amplitude_mm": 5.0 },
    { "op": "settle",   "duration_s": 20.0 }
  ],
  "output":     { "dir": "out/dispersed_seeding", "csv": true, "volume": true }
}
```

## 4. User-input samples (create verbatim, then RUN both)
- `sample_gentle.json`: the schema block above (single jet, `rate_uLs:179`,
  `agitate.count:0`, partial `withdraw`) → expected **low CV**.
- `sample_harsh.json`: `withdraw.rate_uLs:2000, depth_frac:0.2, volume_frac:0.9`;
  `eject.rate_uLs:4319`; `agitate.count:20`; pre-`settle.duration_s:5.0`
  → expected **higher CV**. (Same grid/fluid/particles/geometry blocks.)

## 5. Outputs
- `out/.../density.csv` — `bin_i, bin_j, x_center_m, y_center_m, count, normalized_density`.
- `out/.../metrics.json` — `{ CV, max_over_mean, n_deposited, n_extracted,
  Re_jet, St, Fr, Ma, tau }`.
- `out/.../*.vtk` (or the volume format the qa-viewer consumes) — the 3D reservoir
  and tray fields, so the run is viewable in 3D. If a public volume-export util
  exists, reuse it; otherwise emit legacy ASCII VTK. Do NOT add a new heavy dep.
- stdout: one-line non-dimensional regime summary + final CV / Max-Mean.

## 6. Run interface
```bash
cargo run --release -p lbm-cli --example dispersed_seeding -- \
  crates/lbm-cli/examples/dispersed_seeding/sample_gentle.json
cargo run --release -p lbm-cli --example dispersed_seeding -- \
  crates/lbm-cli/examples/dispersed_seeding/sample_harsh.json
```

## 7. Acceptance gates (report each with evidence in your final message)
1. `cargo build --release -p lbm-cli --example dispersed_seeding` compiles clean.
2. Both samples run to completion and write finite `metrics.json`.
3. **Trend check:** `CV(gentle) < CV(harsh)` — report BOTH numbers. This is the
   physical-credibility gate. If it fails, fix the MODEL, never tune the samples
   to force it; if you cannot make it pass, document the physics reason in
   SPEC_FINDINGS.md.
4. A 3D volume file is emitted for at least one sample.
5. `git status` shows changes ONLY under the example dir (+ the one `[[example]]`
   stanza in `crates/lbm-cli/Cargo.toml`).
6. `README.md` states model + non-goals; `SPEC_FINDINGS.md` lists discovered
   spec issues with proposed fixes.

## 8. Non-goals (do NOT implement)
Free-surface / gas-liquid interface; two-way or four-way coupling;
particle-particle collision dynamics; the inverse-design solver; any CLI
subcommand or core API change; GPU.

## 9. Commit / handoff
Try `git add -A && git commit`. If the sandbox cannot commit (shared-.git
`index.lock` EPERM is a known intermittent failure), **leave the tree staged and
say so explicitly in your final message** — the PM will commit on your behalf.
Do not treat an inability to commit as a failed order.
