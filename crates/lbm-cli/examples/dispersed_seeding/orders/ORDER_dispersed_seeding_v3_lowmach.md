# ORDER v3 (P1.1) — bring `examples/dispersed_seeding` into the low-Mach band

## CONTEXT
v2 is physically correct but both samples run at **Ma ≈ 0.25–0.27**, just under
the 0.3 hard guard — a compressibility-error regime. Bring both samples to
**Ma ≤ 0.1** while keeping the run **stable and the physics unchanged**, then
re-verify the deposition metrics and re-freeze the gentle CV band.

Same isolation as before: changes ONLY under
`crates/lbm-cli/examples/dispersed_seeding/` (no core, no other crates).

## THE NUMERICAL CONSTRAINT (do not naively shrink dt)
`Ma = u* / cs`, `cs² = 1/3`, `u* = U · dt/dx`, and `ν_phys = ν* · dx²/dt` with
`ν* = cs²(τ − 0.5)`. At **fixed dx**, reducing `dt` lowers `u*` (good) but ALSO
lowers `ν*`, driving **τ → 0.5 (unstable)**. So you cannot hit Ma ≤ 0.1 by
shrinking dt alone.

Use **diffusive scaling**: keep `τ` in a safe band (τ ≥ 0.51, ideally ~0.55) and
**refine dx** for the tray so that `u* = U·dt/dx ≤ 0.0577` (⇒ Ma ≤ 0.1) with
`dt ∝ dx²` chosen to hold physical ν fixed. Verify physical viscosity, physical
jet velocity, gravity, and settling velocity are unchanged (SI inputs are
authoritative). Keep the reservoir domain coarse (visualization-only) to control
cost. Respect `max_particle_steps` (raise it if the finer dt needs more steps,
but ABORT with a message rather than silently truncating).

## TASKS
1. Adjust the nondimensionalization so both `sample_gentle.json` and
   `sample_harsh.json` run at **Ma ≤ 0.1** with **τ ≥ 0.51**, physics preserved.
   Prefer refining `dx_m` / grid over changing physical inputs. If a sample
   cannot reach Ma ≤ 0.1 within a reasonable example-scale cost, report the best
   achieved Ma + the cost reason in `SPEC_FINDINGS.md` (do not silently exceed 0.1).
2. Re-run both samples. Report the new `REGIME` line and full `metrics.json`.
3. **Re-freeze the gentle CV band** based on the low-Mach run (the physics
   shouldn't move much; confirm and update the band in `SPEC_FINDINGS.md`,
   noting the old 0.95–1.40 was measured in the compressibility regime).
4. Update `README.md` regime notes.

## ACCEPTANCE GATES (report with actual numbers)
1. `cargo build --release -p lbm-cli --example dispersed_seeding` clean.
2. Both samples: **Ma ≤ 0.1** (or documented best-effort) and **τ ≥ 0.51**.
3. Stability preserved: finite metrics, `empty_bin_fraction(gentle)` still < 0.15.
4. Trend preserved: `CV(gentle) < CV(harsh)` — report both.
5. Isolation: `git status` shows changes only under the example dir.
6. `SPEC_FINDINGS.md` updated with the re-frozen gentle CV band + regime note.

## COMMIT
`git add -A && git commit`. If the shared-worktree `index.lock` EPERM blocks it,
leave the tree staged and say so — PM commits on your behalf.
