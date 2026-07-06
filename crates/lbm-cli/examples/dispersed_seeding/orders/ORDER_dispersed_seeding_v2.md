# ORDER v2 — `examples/dispersed_seeding` (fix degenerate deposition + freeze the spec)

## CONTEXT (what v1 produced)
v1 built a runnable 3D skeleton but the deposition is **degenerate**: in the
gentle sample 63/144 bins are empty and one bin holds 586 of 10000 particles
(CV=2.01). Root cause = single central jet + one-way coupling + **force-deposit
of survivors at a step cap**, so particles fall straight down in a central pile
instead of spreading. v1 also filed 10 spec ambiguities in `SPEC_FINDINGS.md`.
This order resolves those ambiguities (decisions below are FINAL — implement them,
do not re-litigate) and fixes the model so the gentle case actually spreads.

Build on the existing files under `crates/lbm-cli/examples/dispersed_seeding/`.
Same isolation rules as v1: **no changes under `crates/lbm-core/src`,
`lbm-scenario`, `lbm-wasm`, `web`.** If a fix genuinely requires a core API that
does not exist, do NOT patch core — implement the closest in-example
substitution AND log it under "core requirements" in `SPEC_FINDINGS.md`.

## RESOLVED SPEC (implement exactly)
1. **Extraction is an explicit statistical model.** Model the reservoir as a 1D
   vertical settling column `c(z)` (gravitational settling of the size
   distribution over the pre-`settle` duration). `withdraw(depth_frac, volume_frac,
   rate)` samples the extracted batch by concentration at the withdraw depth. The
   3D reservoir LBM domain is run **for visualization only** (emit its VTK), it
   does not drive extraction. Report `n_extracted` and the extracted-diameter
   histogram.
2. **Tray top boundary = a velocity face; zero velocity outside jet patches.**
   Each jet patch is a disk of diameter `nozzle_diameter_m` centered at a
   `points_xy_frac` location, carrying downward velocity `u_jet = Q / (π
   (nozzle_diameter_m/2)^2)`.
3. **`eject` gains a required `nozzle_diameter_m` field** (per point or shared).
4. **`depth_frac`: 0 = filled-liquid surface, 1 = reservoir floor.**
5. **REMOVE force-deposit-at-cap.** Integrate particles until they cross the
   floor plane OR the protocol ends. Any still suspended at the end are counted
   as `n_suspended` and **excluded from the density map** (never projected).
   Add an explicit `max_particle_steps` guard that, if hit, ABORTS with a clear
   message (not a silent central dump).
6. **SI dimensions are authoritative.** Derive grid counts from `dx_m`; if
   provided grid counts disagree with `SI/dx` by more than 1 cell, fail validation
   with a clear message.
7. **Unknown `agitate.pattern` → validation error.**
8. **VTK = `STRUCTURED_POINTS`** vector field in grid order (qa-viewer compatible).

## MODEL FIX (the point of v2 — make gentle actually spread)
The gentle case must demonstrate lateral spreading, not a central pile:
- Resolve the tray **impinging-jet wall jet**: run the tray LBM long enough that
  the downward jet turns into an outward radial wall jet along the floor, and
  advect particles one-way in that live field so they are carried outward before
  settling. Choose tray resolution / step count so this radial transport is
  actually present (verify by inspecting that deposition is not single-binned).
- Support **multiple jet points** (`points_xy_frac` list) — the harsh/gentle
  contrast and multi-point coverage both rely on it.
- Keep one-way coupling, no free surface, no inter-particle collisions (unchanged
  non-goals).

## SAMPLES
Update `sample_gentle.json` / `sample_harsh.json` to the resolved schema (add
`nozzle_diameter_m`, e.g. 0.8e-3 gentle). Keep gentle = low rate / no agitation,
harsh = high rate / agitation, so the trend contrast persists.

## OUTPUTS (extend v1)
`metrics.json` now includes: `CV, max_over_mean, empty_bin_fraction,
n_extracted, n_deposited, n_suspended, Re_jet, St, Fr, Ma, tau`.

## ACCEPTANCE GATES (report each with the actual numbers)
1. `cargo build --release -p lbm-cli --example dispersed_seeding` clean.
2. Both samples run to completion; finite `metrics.json`.
3. **No force-deposit.** Grep your own code to confirm survivors are reported as
   `n_suspended`, never projected onto the floor.
4. **Spreading achieved:** report `empty_bin_fraction` for gentle. It MUST be
   markedly below v1's 0.44; target < 0.15. If you cannot get below 0.15,
   report the achieved value and the physical reason in `SPEC_FINDINGS.md`.
5. **Trend holds:** `CV(gentle) < CV(harsh)` — report both.
6. Isolation: `git status` shows changes only under the example dir (+ Cargo
   stanza already present).
7. Update `SPEC_FINDINGS.md`: mark the 10 v1 items resolved/deferred, add any NEW
   issues found, propose a **frozen CV band for gentle** based on the achieved
   number, and list concrete **core requirements** (e.g. localized source/sink
   BC, mixed-face BC mask) that the in-example substitutions are standing in for.

## COMMIT
`git add -A && git commit`. If the shared-worktree `index.lock` EPERM blocks it,
leave the tree staged and say so — PM commits on your behalf.
