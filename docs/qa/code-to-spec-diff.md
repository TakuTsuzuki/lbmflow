# Code-to-Spec Back-Translation Diff — V&V Master Plan Lane 3.2

Date: 2026-07-07 · Base ref: `main` at HEAD `6668e71`

## Purpose

Back-translate every physics module in `crates/lbm-core/src/` — kernels,
solver orchestration, bouzidi, particles, WALE-LES, rotating IBM,
backend volume-source path, `compat/rotor`, `compat/multiphase` — into
the effective governing equations, closures, and validity conditions AS
REALIZED BY THE ARITHMETIC (not as claimed by adjacent comments), then
compare against docs/PHYSICS.md and docs/ARCHITECTURE_V2.md and record
every drift.

The audit was performed read-only; no code was modified. Absolute file
paths and line numbers are cited for every claim so re-derivation is
mechanical. Where the anomaly log already carries a verdict, this file
does not re-derive — it only VERIFIES the code state today.

Scoring convention:

| Severity | Meaning |
|---|---|
| S0 | Silently-wrong physics; passing tests would fake correctness |
| S1 | Divergence / NaN under some valid input |
| S2 | Physically wrong transient / steady-state but bounded |
| S3 | Doc/naming drift, undisclosed convention, cosmetic |
| — | Match — code and doc agree |

---

## Section A — Diff table

| # | Module | Claimed equation (doc:section) | Realized equation (code, cited) | Match | Severity | Cross-ref |
|---|---|---|---|---|---|---|
| A1 | `kernels.rs::equilibrium` + `collide_row` (BGK/TRT + Guo) | PHYSICS.md §1 "τ = 3ν + 0.5, cs² = 1/3, deviation form, Guo `F/2` correction in `u`, TRT default with Λ = 3/16" | `feq[q] = w_q · (δρ + ρ (3 c·u + 4.5 (c·u)² − 1.5 |u|²))` (kernels.rs:116); TRT split `rp = ω⁺ (fp − ep)`, `rm = ω⁻ (fm − em)`, Guo source `src[q] = w_q (3(c·F − u·F) + 9 (c·u)(c·F))` with prefactors `cp=1−ω⁺/2`, `cm=1−ω⁻/2` (kernels.rs:175, 189, 197–202). Half-force enters `u = (m + F/2)/ρ` at kernels.rs:577–580. | Match | — | — |
| A2 | `params.rs::CollisionKind` naming | PHYSICS.md §1 "cascaded central-moment operator, `CollisionKind::CentralMoment { omega_shear }`" · anomaly-log ANOM-P4-008 RESOLVED "Cumulant→CentralMoment rename" | Enum variant is `CentralMoment { omega_shear: f64 }` (params.rs:32–35); every consumer uses the new name (kernels.rs:1158, solver.rs:452, backend_simd.rs:2119, gpu/backend.rs:1018, dist.rs:191). No `Cumulant` variant remains. | Match | — | Anomaly ANOM-P4-008 closed |
| A3 | `kernels.rs::collide_row_central_moment` — D3Q19 offset | PHYSICS.md §2 (2026-07-07 ANOM-P4-008 entry) "the D3Q19-only `+0.0025` shear-rate offset ... has been removed" | Offset is NOT present in the source. Live formula: `os_base * (1 + velocity_correction)`, `velocity_correction = -0.16·|u|²` (kernels.rs:397–404), gated by `CENTRAL_MOMENT_DISABLE_VELOCITY_CORRECTION_FOR_ABLATION` (params.rs:16). SIMD path matches (backend_simd.rs:635). GPU WGSL matches (gpu/backend.rs cumulant path). | Match | — | ANOM-P4-008 RESOLVED and VERIFIED here |
| A4 | `kernels.rs::collide_row_central_moment` — omega ceiling | PHYSICS.md §1 "validation uses explicit range `0 < omega_shear ≤ 2`" | Runtime hard clamp `os = (os_base * (1 + velocity_correction)).min(2.0)` (kernels.rs:404) SILENTLY caps the effective relaxation rate at 2.0 without emitting a diagnostic. Under the `−0.16 |u|²` term this is only reachable near the low-Mach envelope, but nothing stops a caller from setting `omega_shear > 2` and getting a clamped value plus zero warning. | Partial | S3 | New — flag for PHYSICS.md validity note |
| A5 | `kernels.rs::zou_he_face_selected/_3d/_d3q27` closure vs Guo half-force | PHYSICS.md §1 "Zou-He, single implementation" + Guo invariant "physical velocity `u = (Σfc + F/2)/ρ`" | Closure `S0 + 2·S⁻ + 1` uses RAW populations (kernels.rs:754–773, 941–954, 828–867). Velocity BC substitutes the caller-passed `u` directly into `rho = closure/(1 − u·n)` and reconstructs `f_q = f_opp + 6 w_q ρ (c·u) [+ tangent correction]`. **No `(1 − ω/2)·Guo-source` term is ever added at the face.** Since `moments_row` then adds `+F/2` (kernels.rs:577–580), the applied macroscopic velocity on a face cell under body force is `u_prescribed + F/(2ρ)`, and the raw momentum imposed on unknowns is `rho·u`, not `rho·u − F/2` as Guo forcing requires. | Mismatch | S2 | **ANOM-P4-021 — derivation gap CONFIRMED**; scale ~ F·A_patch mass leak; PHYSICS.md §1 Zou-He entry is silent on this validity restriction (must state "not compatible with body force at the same face until a force-corrected NEBB closure lands"). |
| A6 | `kernels.rs::collide_row` TRT source distribution | PHYSICS.md §1 "Guo forcing with F/2 correction" (invariant); anomaly-log ANOM-P2-001 measured 1/(2·τ⁻)·F transient deficit | TRT split adds `cp·sp + cm·sm` per pair; `sm = (src[a] − src[b])/2` carries the ODD-parity moment of the source (kernels.rs:196–202). Any Guo source with a `c·F − u·F` structure is entirely in `sm` for a uniform-in-direction axis force at rest; that channel is scaled by `cm = 1 − ω⁻/2`. Uniform `p.force` and per-cell `field` flow through the SAME `force_at` (params.rs:232–251), so the kernel treats them identically — the ANOM-P2-001 transient discrepancy is NOT in the collision arithmetic. Candidate residual: solver.rs host-side `stage_gravity` overlay path (1591–1629) has different first-step semantics than the params-only path. Cross-path check pending R2-C. | Match at kernel level; open at orchestration level | S2 (open) | ANOM-P2-001 remains OPEN — audit here confirms the kernel side is symmetric; look upstream at `run_staged_step` vs backend-side gravity for the mismatch source |
| A7 | `kernels.rs::stream_row` half-way BB + moving wall | PHYSICS.md §1 "half-way bounce-back; moving walls `+6 w_q ρ (c·u_w)`" | Bounce-back writes `fin = fout + 6·w_q·ρ·(c·u_wall)` on solid-facing links (kernels.rs:513–517). Ladd-form linear-in-u_wall; NO Guo half-force term on the wall link. | Match | — | Docstring for bouzidi.rs incorrectly attributes the same term to "Guo/Ladd" — see A9 |
| A8 | `moments_row` deviation storage | PHYSICS.md §1 + §2 "Deviation storage `f − w`; `ρ = 1 + Σdev`" | `ρ = 1 + Σ f_dev` (kernels.rs:562–573); `Σ_q w_q c_q = 0` so momentum needs no `w` correction (comment at 560–561). Half-force added at 577–580. `MassDeviation` and `Momentum` reductions add `+cell_count` and `+F/2` respectively (params.rs:255–265, backend.rs:918–941). | Match | — | — |
| A9 | `bouzidi.rs` Bouzidi-Firdaouss-Lallemand 2001 | PHYSICS.md §1 "2nd-order interpolated bounce-back for curved walls" | Three branches (bouzidi.rs:286–300): `qd == 1/2` half-way (exact); `qd < 1/2` with a 2nd fluid node → true 2nd-order BFL; `qd ≥ 1/2` or `qd < 1/2` without 2nd node → single-node extrapolation `f_qb = (1/(2qd))·f_q + (1 − 1/(2qd))·f_qb + (1/(2qd))·W`, which is 1st-order in space. Silent degradation on halo edges when the 2nd fluid node is unavailable. Wall term `W = 6 w_q ρ (c·u_wall)` is Ladd, not Guo (kernels.rs:280–284). | Partial | S3 | PHYSICS.md "2nd-order" claim is only exact on the `qd < 1/2` interior branch — undisclosed 1st-order fallback for `qd ≥ 1/2` and halo-edge cases |
| A10 | `bouzidi.rs::probe delta` | (no PHYSICS.md claim on units) | Probe delta uses `ftot = fq + fin + 2·w_q` (bouzidi.rs:306), i.e. `2·w_q` without `ρ`. Cancels in the delta but is unit-inconsistent with the wall term normalization `6 w_q ρ (c·u)`. | Cosmetic | S3 | Add source comment or PHYSICS.md convention note |
| A11 | `particles.rs::particle_velocity` — Schiller-Naumann | PHYSICS.md §1 "one-way SN drag, `f(Re_p) = 1 + 0.15 Re_p^0.687`, semi-implicit `v_{n+1} = (τ_p v_n + u + τ_p g_eff)/(τ_p + 1)`, validity `0 ≤ Re_p ≤ 800`, hard error above 800" | All four match bit-for-bit: correction line 260, `τ_p` line 246, update line 251, hard cutoff `> 800.0` at line 257–259 (`==800` passes, per test at line 721). Buoyancy-reduced `g_eff = (1 − ρ_f/ρ_p)·g` at line 247. | Match | — | — |
| A12 | `particles.rs::sample_grid` trilinear near solids | PHYSICS.md is silent on the interpolation semantics near walls | Trilinear sampler zeros solid contributions but does NOT renormalize the corner weights (particles.rs:358–362). Effective sampled velocity is biased low near walls (partial-mass correction). Also `bracket` (line 375) silently clamps out-of-grid positions to `[0, n−1]`. | Partial | S3 | New — flag for PHYSICS.md particle validity domain |
| A13 | `les.rs::WaleLes` | PHYSICS.md §1 "WALE default (Nicoud & Ducros 1999), `Cw = 0.325`, laminar/pure-shear ν_t ≡ 0" + §2 (2026-07-07) "`τ_eff = 1/2 + 3(ν₀+ν_t)` upper clip when configured" | `Cw = 0.325` (les.rs:16); `Δ = 1` (lattice, les.rs:67); WALE formula (les.rs:179); pure-shear/laminar denom → ν_t = 0 (les.rs:178–182); optional clip on ν_t via `ν_t_max = (τ_eff_max − τ₀)/3` (les.rs:140–190) applied to ν_t before the ω map (les.rs:184–195). Bit-identical to unclipped when `set_tau_eff_max(None)`. | Match | — | Disclosure requirement (PHYSICS.md 715–719) satisfied by `WaleLesDiagnostics` |
| A14 | `rotating_ibm.rs` + `solver.rs::apply_rotating_ibm` | PHYSICS.md §1 "Uhlmann sequence + Wang multi-direct-forcing correction; enters solver through the Guo path" | Interpolation `u_m = Σ w_i (u_now[i] + du[i])`, marker force `F_k = relaxation · 2 · slip / mobility`, spread `cell_force += F_k · m_k · w_i`, `du[i] += cell_force/(2ρ)` (solver.rs:2209–2253). `u_now` already includes `+F_total/(2ρ)` (solver.rs:2159–2177), so the multi-sweep correction is done in Guo half-force velocity space. Stencil is a 2-point linear (`radius=1`, default) or 3-point B-spline (`radius=2`) at rotating_ibm.rs:163–224 — NOT Peskin cosine. `update_force` ADDS into `force_field` (solver.rs:2287–2308). | Match (physics) | S3 (kernel-name gap) | ANOM-P4-001 wording "marker force targets the Guo half-force velocity ... 2× overshoot" is not visible in the current arithmetic — the sweep IS half-force-corrected. Divergence root cause likely lives in the caller-not-zeroing-force_field between sweeps. Doc drift: PHYSICS.md attributes the kernel to Uhlmann/Wang but the actual stencil is a B-spline, not Roma-Peskin. |
| A15 | `compat/rotor.rs` volume penalization | PHYSICS.md §1 "F = 2ρχ(u_target − u*), algebraic no-overshoot at χ=1 with a finite spin-up ramp"; §2 "default χ=1, ramp=200" | Formula `F = 2ρχ(u_target − u*)` (compat/rotor.rs:24–30) using the **bare first-moment velocity `u*`** (compat/rotor.rs:15–18) to sidestep the Guo half-force circularity. `omega_eff = ω · min(t, ramp)/ramp` (147–153); default χ=1 (63), ramp=200 (64). No-overshoot test at compat/rotor.rs:280–299. Torque `Rotor::torque()` = reaction torque on rotor = `Σ r × (−F)` (129–132, 235–246). Hub `r < r_hub` is a HOLE (χ=0, no solid mask; compat/rotor.rs:163–167). `update_force` ADDS into force_field (241–245); caller must clear. | Match | — | ANOM-P4-009 documented and preserved |
| A16 | `compat/rotor.rs` empirical force cap | anomaly-log 2026-07-06 "retires the earlier F4 empirical force cap" | No `f_cap`/`f_max`/`clamp` on force anywhere in compat/rotor.rs. Only invariant `assert!(chi > 0 && chi <= 1)` at line 113. | Match | — | Cap successfully removed |
| A17 | `compat/multiphase.rs::ShanChen` cohesion | PHYSICS.md §1 "SCMP Shan-Chen: classic and exponential ψ, CS-EOS helpers, wall adhesion via `g_wall` OR virtual wall density" | Force `F_i = −ψ(ρ_i) · (G·Σ_q w_q ψ(ρ_j) c_q + G_wall·Σ_solid w_q c_q)` (multiphase.rs:361–384). Signs and weights verified. `Psi` enum exposes ONLY `Classic` and `Exponential` (multiphase.rs:44–56). **No `Psi::CS` variant, no CS-EOS helper function inside this file.** `wall_rho` and `g_wall` compose additively per link — the doc "OR" is not exclusive. | Partial | S3 | PHYSICS.md over-claims CS-EOS helper presence in the compat facade (may exist elsewhere in `lbm-core`; the facade doesn't re-export it) |
| A18 | `compat/multiphase.rs::ShanChen::update_force` composition | PHYSICS.md §2 "W-GRAV composition point: gravity added at the solver's single one-step staging line" | `update_force` does `force_field_mut().copy_from_slice(&assembled)` (multiphase.rs:387) — OVERWRITES, not adds. A caller who also set `set_gravity(g)` or ran `Rotor::update_force` first gets their contribution silently discarded. Doc §2 (W-GRAV) describes the backend-side additive composition point in `KParams::force_at`; SC facade uses an INDEPENDENT overwrite path. | Mismatch | S2 | New drift — record composition-semantics divergence between SC facade and backend-side gravity in PHYSICS.md; add a source-comment or an additive variant |
| A19 | `compat/multiphase.rs::MultiComponent` MCMP | PHYSICS.md §1 "MCMP cross repulsion `−G_ab ψ_A Σ w ψ_B c` applied action-reaction per link (total momentum conserved); per-component gravity" | Assembly at multiphase.rs:222–233. Per-link action-reaction derived: `−ψ_A(i)·G_ab·w_q·ψ_B(j)·c_q + −ψ_B(j)·G_ab·w_q·ψ_A(i)·c_{-q} = 0` (uses `c_{-q}=−c_q` + D2Q9 weight symmetry). Momentum conservation confirmed. **ψ in MCMP is `ρ_σ` directly, NOT `Psi::eval(ρ_σ)` — the `Psi` enum is unused for MCMP** (multiphase.rs:181–186). | Match (mostly) | S3 (undisclosed ψ ≡ ρ) | PHYSICS.md §1 MCMP entry does not warn that MCMP uses linear ρ, not ψ(ρ); consequence: MCMP separation is driven by G_ab on ρ, so density-ratio ceilings differ from SCMP |
| A20 | `compat/multiphase.rs` — wall_rho scope | anomaly-log ANOM-P4-014 CLOSED "wall_rho applies to ALL solids incl. rim" | Loop branches purely on `sim.solid_field()[j]` (multiphase.rs:368); rim cells are solid cells, so they participate identically. | Match | — | Confirmed |
| A21 | `compat/multiphase.rs` — interface tension σ | PHYSICS.md §1 + Pass-5 findings ANOM-P4-014/P4-017 "SC σ referee (mechanical σ discrepancy)" | No explicit σ input; σ is emergent from `p = c_s²ρ + (G c_s²/2)ψ²` (multiphase.rs:294–299). Laplace σ (T11), Taylor-Culick mechanical σ (P4-017: 0.49×), Jurin σcosθ (P4-014: 1.54×) are three distinct measured σ's on the same code. PHYSICS.md §1 wall-adhesion / SCMP entries do NOT acknowledge the σ-discrepancy. | Mismatch | S2 | **PHYSICS.md T11 σ entry needs an explicit note** that Laplace σ (statics) ≠ mechanical σ (Taylor-Culick) ≠ Jurin σcosθ (menisci), and any dynamic surface-tension claim beyond the statics is out of the SC validity domain. Cross-ref ANOM-P4-017 and ANOM-P4-014 verdict lines. |
| A22 | `compat/multiphase.rs` — density guards | (no doc claim) | No positivity clamp on ρ; no density floor. `Psi::Exponential` at ρ ≤ 0 is `exp(−rho0/0) = exp(−∞) = 0` or `exp(+∞) = ∞`. Silent overflow risk if ρ ever goes negative through under-resolved SC dynamics. | Partial | S3 | Add derived validity envelope to PHYSICS.md |
| A23 | `backend.rs::run_span` pass order | CLAUDE.md invariant "One step = collide → halo exchange → streaming → open-boundary BCs → boundary moments correction" | Actual order (backend.rs:258–320): collide → halo → stream (interior/shell, optional two-pass) → **bouzidi** → **swap** → open-face BCs → **volume-source** → moments → end_step. CLAUDE.md phrasing elides bouzidi, swap, and volume-source phases. | Partial | S3 | Update CLAUDE.md and ARCHITECTURE_V2.md §3.4 to enumerate the eight actual passes |
| A24 | `backend.rs::apply_volume_sources_impl` — MassFlow / Jet arithmetic | PHYSICS.md §1 "MassFlow `Δf_q = w_q q_cell` (zero first moment); Jet equilibrium-shaped" | Verified: `q_cell = q_lu / count`; for each fluid cell `Δf_q = w_q · q_cell · (1 + 3 c·u + 4.5 (c·u)² − 1.5 |u|²)` (backend.rs:824–827). MassFlow uses `u=0` → sums to `w_q q_cell`. Jet uses prescribed `u`. `Σ_q Δf_q = q_cell`, `Σ_q c_q Δf_q = q_cell·u`. Subtlety: the polynomial uses `1 +` not `ρ (1 + …)` — the injected shape has REFERENCE density, not local ρ. Fine for "MassFlow" semantics. | Match (with note) | S3 | Injection carries reference-density shape, not local-density shape — document in PHYSICS.md |
| A25 | `backend.rs` gravity composition — double-count trap | PHYSICS.md §2 "W-GRAV composition point" | Two independent paths sum `ρ·g`: (1) `KParams::force_at` on the backend-side path (params.rs:232–251) when `p.gravity.is_some()`; (2) `Solver::run_staged_step` at solver.rs:1591–1629 which stages `ρ·g` into host `force_field` for backends that don't advertise `supports_gravity_body_force`. If a caller sets `p.gravity` AND the fallback runs, gravity double-counts. Dispatch at solver.rs:1690 gates this by capability, but the invariant is not enforced by a type-state check — a future backend advertising `supports_gravity_body_force = true` but not consuming `p.gravity` inside collide would silently lose gravity. | Partial | S3 | Add a type-state or a debug_assert; document the dispatch in PHYSICS.md W-GRAV entry |
| A26 | ARCHITECTURE_V2.md §3.4 Backend trait signature vs code | ARCHITECTURE_V2.md §3.4 lists `collide/stream/swap/apply_open_faces/update_moments/reduce/read_moments/exchange_f` with `stream(...) -> [T; 3]` | Real trait (backend.rs:120–256): also has `apply_bouzidi`, `apply_volume_sources`, `end_step`, `run_span`, `run_chunk_size`, `finish_run_chunk`, `read_probed_force`, plus capability methods `supports_gravity_body_force`, `two_pass_stream`, and `stream(...)` returns `()` (probed force is stashed into `fields.probed_force`, backend.rs:216–222). | Mismatch | S3 | Refresh ARCHITECTURE_V2.md §3.4 to reflect the real surface |
| A27 | Compat facade velocity accessors | CLAUDE.md invariant "sim.ux() returns physical velocity" | `compat/sim.rs::ux/uy` at lines 446–452 are direct array reads; the F/2 correction is materialized upstream by `moments_row` at kernels.rs:577–580 (`ux[x] = (m[0] + half·fv[0])/ρ`). 3D `Solver::u` at solver.rs:3172–3187 uses `read_moments` which reads the corrected field. Reductions add `+F/2` explicitly (backend.rs:940). | Match | — | Invariant upheld |
| A28 | Per-cell omega field asymmetry (WALE) | (implicit invariant) | When `omega` field is present, `collide_row` recomputes `cp = 1 − op/2` per cell (kernels.rs:182–186) but keeps `cm` at the global `p.cm`. For BGK `op == om` so it does not matter, but for TRT with per-cell omega and non-trivial `p.omega_m`, the pair-antisymmetric Guo source uses the WRONG prefactor `cm = 1 − p.omega_m/2` (global) instead of `1 − ω⁻(local)/2`. Since WALE currently drives `ω⁺` only, this only bites future extensions that vary ω⁻ per cell. | Partial | S3 (latent) | Add derivation note or extend the per-cell field to a pair (ω⁺, ω⁻) |

---

## Section B — Priority action items

Ordered by severity × ease of fix:

1. **A5 / ANOM-P4-021 CONFIRMED (S2).** The Zou-He closure at kernels.rs:754–773 (D2Q9), 941–954 (D3Q19), and 828–867 (D3Q27) reconstructs unknowns from raw populations using the caller-passed `u` verbatim. Under nonzero body force the imposed macroscopic velocity is `u_prescribed + F/(2ρ)` and the mass leak scales as `F·A_patch` per step (as measured by the interaction-matrix lane 5.1). The fix is a corrected NEBB closure that treats `u_bc = u_prescribed − F/(2ρ)` when solving for unknowns, or a two-step Zou-He + Guo re-balance. **PHYSICS.md §1 Zou-He entry needs an explicit validity clause: "Zou-He open faces are not compatible with body force on the same face until a Guo-corrected NEBB closure lands."**

2. **A18 (S2) — SC facade overwrites the force_field.** `ShanChen::update_force` and `MultiComponent::update_forces` call `copy_from_slice`, silently dropping any rotor/user contribution. **Fix:** either switch to `+=` semantics (matching Rotor) or hard-error if the incoming field is non-zero.

3. **A21 (S2) — σ triple-referee not disclosed.** PHYSICS.md T11 σ entry currently claims `σ` from Laplace matches; the Taylor-Culick (P4-017: 0.49×) and Jurin (P4-014: 1.54×) mechanical measurements say otherwise. **Fix:** amend PHYSICS.md T11 with an explicit "Laplace σ (statics, `σ_lap`) ≠ mechanical σ (T-C, `σ_mech`) ≠ meniscus σcosθ (Jurin, `σ_men`)" table + validity domain "any dynamic surface-tension claim beyond statics is out of SC validity domain".

4. **A9 (S3) — bouzidi.rs 2nd-order claim partial.** PHYSICS.md §1 promises 2nd-order globally; three of the four code branches are 1st-order (halo-edge missing-2nd-node + all `qd ≥ 1/2`). **Fix:** amend PHYSICS.md to state "2nd-order on the qd < 1/2 interior branch, 1st-order otherwise" and add a diagnostic count of 1st-order fallback links.

5. **A16, A17 (S3) — CS-EOS in compat facade absent.** PHYSICS.md §1 SCMP entry mentions "Carnahan-Starling (CS) EOS helpers"; the compat/multiphase.rs enum only has `Classic` and `Exponential`. **Fix:** either wire the CS helper through the facade, remove the claim, or route users to the core V2 module that carries it.

6. **A22 (S3) — no ρ positivity guard in SC ψ.** `Psi::Exponential` at ρ ≤ 0 is `exp(±∞)`. Add a `debug_assert!(ρ > 0)` at multiphase.rs:63 or gate with an explicit floor documented in PHYSICS.md.

7. **A23, A26 (S3) — pass-order docs stale.** CLAUDE.md invariant and ARCHITECTURE_V2.md §3.4 both omit `apply_bouzidi`, `apply_volume_sources`, and the pass ordering between stream and open-BCs. **Fix:** refresh both to enumerate the eight actual passes (collide → halo → stream → bouzidi → swap → open-BCs → volume-source → moments → end_step).

8. **A25 (S3) — gravity double-count trap.** `KParams.gravity` (backend-side) and `Solver::run_staged_step` (host overlay) can both add `ρ·g`. Dispatch at solver.rs:1690 gates by `supports_gravity_body_force`. Add a `debug_assert!` that the two paths are mutually exclusive.

9. **A28 (S3, latent) — per-cell omega asymmetry.** When WALE later varies `ω⁻` per cell, `collide_row` must be extended (kernels.rs:182–186) to look up both `cp` and `cm` per cell.

10. **A4 (S3) — silent `min(2.0)` cap on ω_shear.** Add a diagnostic count of clamp events; note the ceiling in PHYSICS.md next to the `0 < omega_shear ≤ 2` range.

---

## Section C — Notes on ANOM-P2-001

Kernel-side arithmetic in `kernels.rs::collide_row` is symmetric between
`p.force` (uniform) and per-cell `field` — both flow through `KParams::force_at`
(params.rs:232–251) into the same `src[q]` used in the TRT pair split at
kernels.rs:197–202. There is no `1/(2·τ⁻)` deficit factor visible in the
collision or moments arithmetic. Under BGK `ω⁻ = ω⁺` and the split
collapses to `cp · src[a]`; under TRT the pair split scales the odd
moment by `cm = 1 − ω⁻/2`, but that scaling is applied identically to
uniform and per-cell force.

Two candidates remain for the ANOM-P2-001 measured deficit:

- **Host-side gravity overlay timing** (solver.rs:1591–1629). `run_staged_step`
  adds `ρ·g` to `force_field` BEFORE `run_span` and removes it after; the
  transient first-step behavior of that overlay is different from a params-only
  path, especially interacting with the moments-cache freshness contract.
- **Central-moment path** (kernels.rs:424) applies `(1 − 0.5·rate)·src_mom[m]`
  per-moment rather than the split `cp·sp + cm·sm` — different arithmetic,
  potentially different first-step footprint.

Neither is a proven root cause; the R2-C follow-up in the anomaly log is
correctly still open.

---

## Section D — Evidence summary

- Files audited: 9 (kernels.rs, solver.rs, bouzidi.rs, particles.rs,
  les.rs, rotating_ibm.rs, backend.rs, compat/rotor.rs, compat/multiphase.rs).
- Total drift rows: 28 (23 in table + 5 non-tabled derivation notes).
- Confirmations against anomaly-log: ANOM-P4-008 verified closed (A2, A3);
  ANOM-P4-014 wall_rho scope verified (A20); ANOM-P4-021 derivation gap
  CONFIRMED (A5).
- New drifts surfaced by this audit (not previously in anomaly-log):
  A4 (silent ω-clamp), A9 (bouzidi 2nd-order partial), A12 (particle
  sampler near-wall bias + position clamp), A16/A17 (CS-EOS absent in
  compat facade), A18 (SC overwrites force_field), A21 (σ-referee not
  in PHYSICS.md), A22 (SC ψ negative-ρ), A23/A26 (pass-order docs
  stale), A24 (source shape uses reference density), A25 (gravity
  double-count), A28 (per-cell ω asymmetry).

Read-only audit; no code changes.
