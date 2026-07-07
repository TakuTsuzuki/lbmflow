---
name: lbmflow-user-tune-stability
description: >-
  Choose the LBMFlow collision scheme (bgk|trt) and fix numerical-stability
  settings — tau floor, low-Mach velocity limit, grid-Reynolds number — so a
  scenario runs without diverging. Use whenever the user asks "which collision
  should I use", "bgk or trt", "why did my run diverge / go NaN", "why is it
  unstable", "how do I fix this validate warning about tau/Mach/grid Reynolds",
  "make it stable", or "what nu/velocity/resolution should I pick". This Skill
  owns the two accepted collision options and the three stability levers the
  validator checks (tau≥~0.55, |u|≤0.15 Ma / hard 0.3, U/ν≤15). Do NOT use it to
  build the scenario skeleton or geometry from scratch (that is
  lbmflow-user-author-scenario) or to run a fixed preset (that is
  lbmflow-user-run-preset) — this Skill only adjusts collision + stability knobs
  on an already-drafted scenario.
---

# LBMFlow — choose collision scheme & stabilize a run

Divergence in LBM is almost always one of a few well-understood causes: the
relaxation time `tau` sitting too close to 0.5, the lattice velocity being too
high (compressibility error), or the grid-Reynolds number being too large for the
resolution. This Skill picks the collision scheme and turns those three
`lbm validate` warnings into concrete parameter fixes. It operates on a scenario
that already exists (drafted by `lbmflow-user-author-scenario` or a preset shown
via `presets show`).

**Everything here is in LATTICE UNITS.** There is no SI conversion; `nu`, `u`,
`rho` are lattice quantities.

## Collision scheme — pick bgk or trt (nothing else exists)

`physics.collision.type` accepts exactly two values. Recommend ONLY these:

| Scheme | When to choose it |
|---|---|
| `trt` (default, recommended) | Default for essentially all cases. Its viscosity-independent stability makes it robust at low `nu` / near the tau floor. Choose it unless there is a specific reason not to. |
| `bgk` | Simplest single-relaxation scheme; fine at moderate `nu` (tau comfortably > 0.6). Slightly cheaper. More sensitive near the stability limit. |

MRT, cumulant, regularized, etc. are **not implemented** — `validate` returns
`unknown variant \`mrt\`, expected \`bgk\` or \`trt\``. Never recommend them.
When in doubt, `trt`.

## The three stability levers (what the validator checks)

`lbm validate` emits a graded warning for each. The relationships you tune with:

- **`tau = 3*nu + 0.5`** (cs² = 1/3). Stability wants `tau` comfortably above 0.5.
- **Grid Reynolds** `= U / nu` where `U` is the max inlet/wall speed.
- **Lattice Mach** ∝ `U` (cs = 1/√3 ≈ 0.577), so `U` is the compressibility knob.

| Warning (from `validate`) | Trigger | Physical cause | Fixes (any one) |
|---|---|---|---|
| `tau = X は安定限界に近い（0.55 未満）` | `tau < 0.55`, i.e. `nu < ~0.0167` | Relaxation too fast; near the BGK/TRT stability floor | Raise `nu` so `tau ≥ 0.55` (i.e. `nu ≥ 0.0167`); or increase grid resolution and rescale; prefer `trt`. |
| `流入/壁速度 X …（0.15 超）` | `max |u| > 0.15` | Compressibility (Mach) error; hard low-Mach cap is **0.3** (build error above it) | Lower the inlet/wall speed below 0.15 (aim ≤ 0.1); keep it under 0.3 always. |
| `グリッドレイノルズ数 U/ν = X > 15` | `U/nu > 15` | Under-resolved for this Reynolds number → divergence risk | Lower `U`, raise `nu`, OR raise grid resolution (more cells for the same physical Re). |

**Note on precision:** for multiphase, `f32` triggers a warning (interfaces need
`f64` headroom). Prefer `precision: "f64"` (the default) for sharp-interface or
low-nu work.

## Decision procedure — stabilizing a diverged / warned run

1. **Run `lbm validate scenario.json`** and read the `warnings` array (and any
   `error`). This is the ground truth — tune to what it actually reports, not a
   guess.
2. **If `ok:false` with a hard `error`** (e.g. speed `0.5` exceeds low-Mach limit
   `0.3`): that is a build blocker, fix it first — bring `max |u|` under 0.3
   (ideally ≤ 0.15).
3. **For each warning**, apply a fix from the table. Prefer the fix that changes
   physics least:
   - tau warning → raise `nu` to hit `tau ≥ 0.55` is usually cheapest.
   - grid-Re warning → if you cannot lower `U` (it is the physics you want),
     raise resolution; raising `nu` changes the Reynolds number you are modeling.
   - Mach warning → lower `U`; if `U` is fixed by the target Re, raise resolution
     and lower `U` together.
4. **Choose the scheme:** `trt` unless the user specifically wants `bgk` and
   `tau` is comfortably > 0.6.
5. **Re-validate.** Iterate until warnings clear (or are consciously accepted).

## Verification gate — the done check

Stabilization is done when:

```bash
./target/release/lbm validate scenario.json
```

returns `ok:true` AND the `warnings` array is empty — **or** the remaining
warnings are explicitly acknowledged by the user as acceptable (e.g. a mild
grid-Re warning they accept for speed). A clean run:

```json
{ "error": null, "ok": true, "warnings": [] }
```

If you want to confirm no *runtime* divergence (NaN mid-run), a short run is the
harder check: a diverged run reports `status:"diverged"` and halts early —
`{ "status": "diverged", "stepsRun": 1000, ... }` means the fix did not hold, so
tune further. (Running is the runner's job — this Skill's own gate is the clean
`validate`.)

## Worked example (end-to-end)

Task: "My run goes NaN. nu=1e-4, inlet u=0.5, trt." Stabilize it.

1. **Validate:** reports `ok:false`, `error:"prescribed speed 0.5 exceeds the
   low-Mach limit 0.3"`, plus warnings: tau≈0.5003 near-limit, edge-speed 0.5
   compressible, grid-Re = 5000 ≫ 15.
2. **Hard error first (Step 2):** `u=0.5 > 0.3` blocks build. Lower inlet to
   `u:[0.08,0]` (≤ 0.15).
3. **tau warning:** `nu=1e-4` → `tau=0.5003`, far too low. Raise to `nu=0.02`
   → `tau=0.56` (≥ 0.55). Clears.
4. **grid-Re:** now `U/nu = 0.08/0.02 = 4 < 15`. Clears. (At the old `nu` it was
   thousands.)
5. **Scheme:** keep `trt` (robust at these settings).
6. **Re-validate:** `ok:true`, `warnings:[]`. Done. Report the three changed
   numbers and why (Mach cap, tau floor, grid-Re).

## Top failure modes (and the fix)

- **Recommended a scheme that doesn't exist.** Only `bgk`/`trt`. Fix: default to
  `trt`; never MRT/cumulant/etc.
- **Fixed only one warning and declared done.** Divergence usually needs ALL
  three levers in range. Fix: clear every warning (or get explicit acceptance).
- **Raised `nu` to fix grid-Re but changed the physics.** Raising `nu` lowers the
  Reynolds number you are simulating. Fix: if the Reynolds number matters, raise
  resolution instead of `nu`.
- **Ignored the hard 0.3 Mach cap.** Above 0.3, build errors — a warning at 0.15
  is advisory, 0.3 is a wall. Fix: keep `max |u| ≤ 0.15`, never above 0.3.
- **Treated physical velocity as lattice velocity.** "1 m/s" is not a lattice
  speed. Unit conversion is out of scope → route to
  `lbmflow-user-author-scenario`'s unit-conversion note; tune in lattice units.
- **Guessed the warnings instead of running validate.** Fix: always base the fix
  on the actual `validate` output for THIS scenario.
