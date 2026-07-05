# REQ_STIRRED_REACTOR.md rev.1b Second-Round Adversarial Review Findings

Target: `docs/REQ_STIRRED_REACTOR.md` rev.1b
Scope: Document review only. No code changes, no cargo execution.

## Numbered Findings

1. **Critical — The "tolerance against the fidelity reference solution" for relaxed extensions is undefined in §8/T17**

   Relevant locations:

   > L22: Relaxed modes are validated by tolerance against the corresponding fidelity reference solution (thresholds defined in §8 VR).

   > L26-L33: MRF / point-bubble / passive / one-way / block-AMR / aggressive f32 are listed as relaxed extensions.

   > L221-L227: VR-STR-01–07 covers only physics benchmarks, conservation laws, and initialization independence.

   > `docs/VALIDATION.md` L322-L330: T17 has only the VR-STR-01–07 table.

   Problem: §1 declares that relaxed modes are validated by "tolerance against the fidelity reference solution," but §8/T17 has no comparison target, measured quantity, or tolerance defined per relaxed mode. In particular, the relative-degradation judgment for `MRF-frozen-rotor` vs `IBM-inertial/sliding-overset`, `point-bubble` vs `resolved-phasefield`, `one-way` vs `two-way`, `block-AMR` vs `uniform`, and aggressive `f32` vs the fidelity profile is undefined, making §1's acceptance criteria unexecutable.

   Concrete fix: Add a "relaxed-mode equivalence" group to §8 as VR-STR-08 onward, or as VR-STR-RELAX-*, fixing the following for each axis:

   - MRF: tolerance for `Np`, discharge velocity profile, mean velocity field, and torque against an IBM/overset reference at the same geometry and Re.
   - point-bubble: tolerance for `ε_g`, `d_32`, `k_La`, and momentum/scalar budgets against a resolved-phasefield reference. However, bound the range of applicability by `d_b/Δx`, `d_b/W`, and `α_g`.
   - one-way: allowable mass-loading range for particle statistics and neglect of fluid momentum reaction, against a two-way reference.
   - AMR: tolerance for conserved quantities, interface position, torque, velocity-field norms, and budget error at coarse-fine crossings, against a uniform reference.
   - Aggressive f32: tolerance for conserved-quantity drift, `Ca_spurious`, `Np`, interface curvature, and reduced quantities, against the fidelity profile or a full-f64 reference.

2. **Major — "batch full-subsystem implementation" and "relaxed extensions (later)" conflict in scope**

   Relevant locations:

   > L14: The delivery scope maintains batch (all subsystems implemented simultaneously).

   > L22: Low-cost approximations (MRF, point-bubble, one-way, AMR, aggressive f32) are added as deferred extensions behind the same trait.

   > L26: Relaxed extensions (later) / reference-grade

   Problem: Under rev.1a's fidelity-default policy, low-cost approximations are a deferred extension point. However, L14 maintains simultaneous implementation of all subsystems, which reads as if MRF, point-bubble, one-way, AMR, and aggressive f32 are also targets for the initial delivery. The delivery scope and the trait/API reservation scope are conflated.

   Concrete fix: Clarify L14 to state "fidelity-default subsystems are implemented in a batch. Relaxed extensions reserve the trait boundaries, config schema, and validation items in the first version; implementation is deferred." Also add a column to the §1 table's "relaxed extensions (later)" row and the §4 FR items indicating whether each is mandatory for the first version or merely API-reserved.

3. **Major — The surface-tension convention when active σ is variable contradicts REQ §3 and the active-scalar proposal**

   Relevant locations:

   > L84-L87: Surface tension (fixed to chemical-potential form) ... `F_s = μ_φ ∇φ`

   > L255: Concrete formulas and stabilization (including Marangoni) for the `active` scalar's feedback targets (σ, viscosity, density, [temperature]).

   > `docs/proposals/active-scalar-feedback.md` L42-L48: When σ is variable, do not use `F_s = μ_φ∇φ`; switch to the well-balanced CSF/chemical-potential-combined form instead.

   Problem: REQ §3 fixes surface tension to `F_s = μ_φ∇φ`, but for active σ, the proposal document adopts a convention that does not use this formula, in order to avoid double-counting the Marangoni tangential force and the normal capillary force. Since the REQ side does not condition this on "only when σ is constant," it conflicts with the active default.

   Concrete fix: Revise the surface-tension section of §3 to the following structure:

   - σ constant: as before, take `μ_φ` and `F_s = μ_φ∇φ` as the baseline form.
   - σ depends on `C_k` or temperature: do not directly use `F_s = μ_φ∇φ`; unify on the well-balanced CSF/chemical-potential-combined form.
   - Add a degeneracy test to §8/T17 confirming agreement with the existing form when degenerating to σ=constant.
   - State explicitly that the coefficients must be derived against REQ's `(κ,β,W,σ)` convention and then frozen.

4. **Major — Active density feedback / Boussinesq force is not reflected in REQ's governing equations or coupling flow**

   Relevant locations:

   > L75: `+ F_s + ρ g + F_g^{disp} + F_p + F_rot`

   > L98-L101: Scalar/reaction equations cover only ADE and interfacial mass transfer.

   > L167: Force-source composition (`F_s+ρg+F_g+F_p+F_rot`)

   > L255: `active` scalar's feedback targets (σ, viscosity, density, [temperature])

   > `docs/proposals/active-scalar-feedback.md` L85-L97: `F_b = ρ_0 β_C (C − C_0) g` is added as a Boussinesq perturbation force.

   Problem: §1 defaults to `active`, and §10 lists density feedback as a remaining implementation detail, yet §3's momentum equation and §5's force-source composition have no scalar-derived Boussinesq force. Since the active-scalar proposal treats this as a perturbation force separate from `ρ(φ)g`, the REQ side also needs a position for `F_b^{scalar}` and a convention separating it from well-balanced gravity.

   Concrete fix: Add `F_b^{scalar}` to §3's momentum equation and §5's FR-COUP-01 force-source composition, and state explicitly that it is "exactly 0 at `C=C_0`, and not mixed with the well-balanced hydrostatic cancellation of `ρ(φ)g`." Add a degeneracy validation to §8 confirming the same static behavior as VR-STR-06 when active is ON and `C≡C_0`.

5. **Major — VALIDATION T17 drops the "energy-like quantity" drift from REQ VR-STR-05**

   Relevant locations:

   > L225: Set individual drift thresholds for mass, momentum, total scalar amount, gas-phase volume, particle count, and energy-like quantities.

   > `docs/VALIDATION.md` L328: Set individual drift thresholds for mass, momentum, total scalar amount, gas-phase volume, and particle count.

   Problem: REQ §8 includes "energy-like quantities" among the conservation/regression diagnostics, but this is missing from the transcription into VALIDATION T17. This makes the acceptance table inconsistent with REQ, despite rev.1b's claim that T17 is already wired up.

   Concrete fix: Add "energy-like quantities" to the VR-STR-05 row of `docs/VALIDATION.md` T17, and on the REQ side, add candidate definitions for energy-like quantities — e.g., kinetic energy, interfacial free energy, particle kinetic energy — scoped as "monitored for unphysical drift, not strict conservation."

6. **Major — The granularity of bubble-swarm validation separation is weaker in the REQ §8 and T17 table**

   Relevant locations:

   > L222: Separate validations for single bubble (`U_t` ...) / bubble swarm / stirred-tank aeration (`ε_g, d_32, k_L a` ...).

   > `docs/VALIDATION.md` L325: Gas-liquid (3-way separation: single bubble / bubble swarm / aerated stirred tank) ... single bubble `U_t` ... `ε_g, d_32, k_L a` are experimental correlation ratios.

   Problem: T17's target column maintains the 3-way separation, but the acceptance-criteria column only has metrics for single bubble and aerated stirred tank, with none for bubble swarm alone. It is ambiguous which validation covers bubble-swarm coalescence, breakup, swarm rise velocity, holdup, and BIT.

   Concrete fix: Split VR-STR-02 into 02a/02b/02c, and for 02b (bubble swarm), set at least `ε_g` distribution, swarm rise velocity, `d_32` (if coalescence/breakup is included), and turbulence intensity or `ν_t` response (if BIT generation is included). When comparing point-bubble and resolved-phasefield, also link to the relaxed-mode-equivalence validation in Finding 1.

7. **Major — The precision default and NFR-02's "f32 default" wording retain the old policy**

   Relevant locations:

   > L32: Fidelity profile: near-interface, conserved quantities, torque, interface curvature, and reductions use `f64`; only the far-field bulk uses `f32`.

   > L209: `f32-bulk` + conserved quantities/torque/interface curvature/reductions use f64. The f32 default's range of applicability is limited to "single-phase/weakly coupled."

   Problem: L32's default is the "fidelity profile," where only the far-field bulk uses f32. Meanwhile, L209 retains the phrase "f32 default," which reads as if the old policy — where aggressive f32 was the default — still applies. This is inconsistent in terminology with rev.1a's PM decision that "the default is fidelity-first."

   Concrete fix: Rewrite L209 as "the aggressive-f32 relaxed mode's range of applicability is limited to single-phase/weakly coupled. The fidelity default is the `f32-bulk + critical f64` profile." Add a relative-degradation validation for aggressive f32 to §8.

8. **Minor — The amortized value for f64 promotion of the interface band in the §7 memory table does not fully match the assumed 5–10%**

   Relevant locations:

   > L192-L193: Fluid distribution f = 216 B/cell, phase-field distribution g = 152 B/cell.

   > L198: f64 promotion of the interface band (band width ~2W, assuming 5–10% of all cells, f+g equivalent) | +30–40

   Verification: The f32 ping-pong of `f+g` is `216+152=368 B/cell`. Promoting this to f64 adds the same amount, +368 B/interface-band cell, so amortized over 5–10% of all cells this is `+18.4–36.8 B/cell`. The table's `+30–40` corresponds to roughly 8–11% of the interface band, or to a scenario where something beyond `f+g` is also promoted.

   Concrete fix: If keeping the 5–10% assumption, revise to `+18–37 B/cell`. If keeping `+30–40`, either set the band fraction to `~8–11%`, or state explicitly that the promotion target also includes `μ_φ, ∇φ, curvature/reduction working arrays`, etc. The total `≈560–620 B/cell` is currently rounded toward the high side in the table, so align the breakdown and the rounding method.

9. **Minor — NFR-01's opening statement that "1e9-cell grids with multiple distributions are TB-class" is inconsistent with the conversion in the same section**

   Relevant locations:

   > L186: 1e9-cell grids with multiple distributions are TB-class.

   > L201-L202: 1e9 cells ≈ 0.56–0.62 TB; even at the full-f64 reference grade, ≈ 1.1–1.2 TB.

   Problem: The budget table in the same section has the fidelity default at 0.56–0.62 TB, and even full f64 at 1.1–1.2 TB, which is hard to describe as "TB-class" (plural). It could reach several TB when including multiple scalars, simultaneous checkpoint retention, working buffers, and statistics history, but L186 does not state these conditions.

   Concrete fix: Revise L186 to "1e9-cell grids with multiple distributions are 0.6-TB class; with full f64, multiple scalars, and I/O buffers included, this reaches TB to several-TB class."

10. **Minor — The §2 heading is out of alignment with its positioning after title neutralization**

    Relevant locations:

    > L7-L8: Representative application problem: stirred-tank reactor (functional requirements are defined domain-neutrally; §2 and §8's validation benches instantiate this application).

    > L37: Target problem, representative quantities, dimensionless numbers.

    > L39: 3D cylindrical (or rectangular) vessel ... bottom sparger ... constant-angular-velocity `Ω` rigid-body-rotation impeller ...

    Problem: The title has been neutralized to be domain-independent, and L7-L8 also frames the stirred tank as a "representative application problem," yet the §2 heading remains "target problem," and the body text begins with stirred-tank-specific conditions only. This makes it appear to the reader that the entire document's target is fixed to the stirred tank.

    Concrete fix: Change the §2 heading to "Representative application, representative quantities, and dimensionless numbers (stirred-tank reactor)," and add a sentence at the opening stating "the following is a representative application of the §8 validation benches; the functional requirements in §4 apply generally to rotating boundaries, high-density-ratio two-phase flow, and LES coupling."

11. **Minor — The relationship between the D3Q27 default condition and the fidelity default is ambiguous**

    Relevant locations:

    > L22: The default for each axis is the fidelity-first implementation (= the reference solution).

    > L114: The condition for defaulting to D3Q27 is limited to "multiphase, or strong forcing, or when using cumulant."

    Problem: Since the M-F fidelity-default configuration includes high-density-ratio two-phase flow, strong forcing, and cumulant candidates, it should effectively always be D3Q27, yet FR-CORE-01 makes this a conditional default as a general core requirement. Whether D3Q19 is permitted for single-phase, weak forcing, or whether the M-F fidelity reference solution is always D3Q27, is not made explicit between §1 and §4.

    Concrete fix: Add to §1 or FR-CORE-01 the statement "the M-F fidelity-default scenario falls under the multiphase/strong-forcing condition and therefore selects D3Q27. Derived scenarios with single-phase, weak forcing may use D3Q19."

## Perspective-by-Perspective Confirmation

- §1 configuration matrix vs §4 functional requirements vs §8 VR: Findings 1, 2, 7, 11. In particular, relative-error validation for relaxed extensions is not wired up.
- §7 memory-budget table arithmetic: the main rows `D3Q27×2×f32=216 B`, `D3Q19×2×f32=152 B`, `D3Q7×2×f32=56 B`, `12×f32=48 B`, `1e8=56–62 GB`, `1e9=0.56–0.62 TB`, and GPU `1.3–2.6e7` cells/card × `40–80` cards, are all confirmed. Only the interface-band f64 amortization and the "TB-class" (plural) wording are Findings 8, 9.
- Dimensions and coefficients of all formulas: the general `τ_eff` form, `Np=P/(ρN^3D^5)` and `P=ΩT_q`, Allen-Cahn's `M[length²/time]` and `Pe_φ=UW/M`, and the internal consistency of `σ=√(2κβ)/6`, `W=4√(κ/(2β))`, and `μ_φ` are confirmed with no findings. Boussinesq is a finding (Finding 4) as a gap in reflecting it into the REQ body.
- Consistency between §2 after title neutralization and the representative application problem: Finding 10.
- `docs/VALIDATION.md` T17 and REQ §8: Findings 5, 6.
- `docs/proposals/active-scalar-feedback.md` and REQ §3/§4: Findings 3, 4. In particular, the surface-tension convention when σ is variable requires a conditional branch on the REQ side.

## Findings Count Summary

- Critical: 1
- Major: 6
- Minor: 4
- Total: 11
