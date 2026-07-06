I have enough to compile the backlog. Here is the deliverable.

---

# LBMFlow V&V Benchmark Backlog — Lane 4.2

**Scope:** standard CFD / FSI / particle-laden / multiphase validation benchmarks the simulator does *not* yet run, ranked by immediate usability. Excludes everything already covered (TGV 2D/3D, Poiseuille/Couette, Ghia + Albensoeder-Kuhlmann cavities, Schäfer-Turek 2D-1/2D-2, duct series, SN sphere-drag staircase, Laplace/contact-angle/RT Shan-Chen, Stokes first problem, acoustic decay, annular-Couette torque, SN settling, sink/dipole).

**Capability-status legend**
- **READY NOW** — runs on shipped capabilities, no code changes.
- **READY (heavy)** — runs on shipped capabilities but 3D/turbulent, large compute.
- **IN-FLIGHT FIX** — blocked only on the marker-IBM moving/rotating fix already in flight.
- **NEEDS TRACK: rigid-body** — needs a free-flying rigid body with hydrodynamic force/torque feedback (1–6 DoF) + two-way coupling. **No such track exists yet — propose MF-η "resolved rigid body."**
- **NEEDS TRACK: two-way particle** — needs two-way point/finite-size particle coupling (momentum feedback), beyond the current one-way Lagrangian model.
- **NEEDS MF-γ** — needs the planned phase-field high-density-ratio two-phase track.
- **OUT OF SCOPE** — needs elastic solids, free surface, or an energy equation; none planned.

Reference-data provenance is stated explicitly per item. "Digitize" = values exist only as published *figures* and must be extracted; "table/analytic/correlation" = usable numbers exist directly.

---

## Part 1 — READY-NOW ranking (by information-per-compute)

| Rank | Benchmark | Compute | What it newly validates | Reference data | Proposed acceptance metric |
|------|-----------|---------|-------------------------|----------------|----------------------------|
| **1** | **Kovasznay flow** | Tiny: 64²→256² convergence ladder, ~1e4 steady steps | First **exact full nonlinear steady NS** check with a nontrivial **pressure field**; measures observed spatial order on the real solver (not a linear mode) | **Analytic** (closed form) | L2 velocity + pressure error vs exact; observed order ≥ 1.9 across grid ladder at Re=40 |
| **2** | **Periodic cylinder-array permeability (Sangani-Acrivos)** | Tiny: 64²–128² periodic unit cell, low-Re/Stokes | Stokes drag / porous-closure accuracy and half-way-BB curved-surface drag as a function of solid fraction — direct porous/Darcy anchor | **Series solution** (Sangani & Acrivos 1982, tabulated coefficients) | Dimensionless permeability k/a² vs series at φ ∈ {0.1,0.3,0.5} within a few % |
| **3** | **Womersley pulsatile channel** | Tiny–small: width 64–128, body-force/pressure drive, ~5–10 periods | First **unsteady exact** check: temporal accuracy, Stokes-layer resolution, oscillatory-BC fidelity | **Analytic** (Bessel-function profile) | L2 profile error at ≥4 phases; amplitude & phase-lag error at α ∈ {4,8,12} |
| **4** | **Impulsively started cylinder (Koumoutsakos-Leonard)** | Moderate: D≈60–100 cells, transient to T=Ut/D≈6, Re∈{550,1000,3000,9500} | Transient drag on a **stationary curved (Bouzidi) boundary**, early wake/secondary-vortex topology — exercises unsteady curved BC | **Digitize** C_D(t) & separation angle from KL 1995 JFM 296:1 (data also mirrored at CSE-Lab KL_JFM_1995 page) | C_D(t) within band vs KL at Re=550,1000; correct secondary-vortex onset time |
| **5** | **Backward-facing step (Armaly)** | Moderate: step S≈20 cells, channel ~40S long, steady, Re≤400 (2D-valid range) | Separation/reattachment length prediction — a classic separated-flow anchor absent from current suite | **Digitize** x_r/S vs Re from Armaly 1983 JFM 127:473 figures | Reattachment length x_r/S vs Armaly at Re∈{100,200,400}; **note: restrict to Re≲400** — 3D effects diverge from 2D above that |
| **6** | **Oscillating cylinder in still fluid (Dütsch)** | Moderate–high: Re=100, KC=5, several periods, prescribed in-line oscillation | Unsteady **moving-boundary** force history + Morison decomposition — first prescribed-motion moving-body validation | **Digitize** phase-resolved u,v profiles at 4 phase angles / 3 x-stations from Dütsch et al. 1998 JFM 360:249 | Velocity profiles vs Dütsch at phases {0°,96°,192°,288°}; C_D history |

**Notes on ranking rationale.** Items 1–3 are exact/series-referenced and nearly free — they should land first and become permanent regression gates (each is a convergence-order or L2-band test, not a figure-matching exercise). Item 2 doubles as the porous-media Darcy anchor (see Part 4). Items 4–5 cost more and rely on digitized figures, lowering info-per-compute. Item 6 sits at the boundary of Part 2: it is READY *only if* prescribed oscillatory motion is available via the volume-penalization path (which already supports prescribed rotor motion, so prescribed translation is a small extension) — if it must go through marker-IBM translation, reclassify as IN-FLIGHT FIX.

---

## Part 2 — READY but heavy (shipped capabilities, large compute; not information-per-compute competitive)

### Taylor-Couette wavy-vortex regime
- **Reference:** critical Taylor number for Taylor-vortex onset (DiPrima & Swinney review; Ta_c well established); Wendt/Lathrop torque scaling correlations (G vs Re) — **correlations/analytic onset**, torque digitizable.
- **Validates:** 3D rotating-boundary (penalization rotor) transition sequence, torque integration, secondary-flow wavenumber — beyond the in-flight annular-Couette *laminar torque* test.
- **Capability:** READY (heavy) — inner-cylinder rotation via volume penalization rotor (shipped).
- **Cost:** 3D annulus, O(10⁷) cells, long transient to reach wavy state. Expensive.
- **Metric:** critical Ta for vortex onset within a few %; torque G(Re) vs correlation; correct wavy azimuthal wavenumber.

### Rushton turbine power number + Costes-Couderc flow field
- **Reference:** turbulent-plateau power number Np ≈ 5.0–5.5 (**published correlation**, e.g. Rushton/Bates); radial-jet mean-velocity profiles from Costes & Couderc 1988 LDV (**digitize**) and Wu-Patterson.
- **Validates:** rotor penalization + WALE LES + torque→power integration in the M-F stirred-reactor context — the headline M-F acceptance case.
- **Capability:** READY (heavy) — rotor penalization + WALE LES shipped; this is squarely M-F / REQ_STIRRED_REACTOR territory.
- **Cost:** 3D tank, O(10⁷–10⁸) cells, turbulent averaging window. Very expensive.
- **Metric:** Np within ~15% of correlation at Re>10⁴; radial-jet velocity/turbulence profiles vs Costes-Couderc within band.

---

## Part 3 — Blocked on in-flight fix or a new rigid-body track

| Benchmark | Reference data | What it validates | Capability status | Cost | Proposed metric |
|-----------|----------------|-------------------|-------------------|------|-----------------|
| **Jeffery orbits (ellipsoid in shear)** | **Analytic** orbit period T=(2π/γ̇)(r+1/r) (Jeffery 1922) | Torque-free rotational DoF of a resolved body in shear | NEEDS TRACK: rigid-body (torque feedback on 1 rotational DoF). Marker-IBM in-flight fix is necessary but **not sufficient** (that fix is *prescribed* rotation; Jeffery needs *free* rotation) | Small once capability exists (single body, 2D/3D shear cell) | Orbit period vs analytic within few % across aspect ratios; drift-free over ≥3 orbits |
| **ten Cate single-sphere sedimentation** | **Digitize** trajectory z(t) & settling velocity from ten Cate et al. 2002 Phys. Fluids PIV; setup 100×160×100 mm, D=15 mm, Re∈{1.5,4.1,11.6,32.2} | 1-DoF free settling with two-way hydrodynamic force feedback + wall approach | NEEDS TRACK: rigid-body (1-DoF translation, two-way) | Moderate–high: ~200×320×200 lattice, settle-out transient | Settling velocity vs PIV at all 4 Re within band; approach-deceleration profile matched |
| **Drafting-kissing-tumbling (two spheres)** | **Digitize** from Fortes, Joseph & Lundgren 1987 JFM 177:467 (qualitative sequence + timing); widely reproduced quantitatively in later IBM papers | Wake-mediated multi-body interaction + collision/lubrication + 6-DoF | NEEDS TRACK: rigid-body (6-DoF, 2 bodies, near-contact model) | High (3D, long transient) | Reproduce D→K→T sequence; kissing time & post-tumble separation vs reference |
| **Segré-Silberberg inertial migration** | **Analytic/experimental** equilibrium radius ≈0.6R (Segré & Silberberg 1962), Re-dependent shift toward wall | Inertial lift on a neutrally buoyant freely-moving finite-size particle | NEEDS TRACK: two-way particle **or** rigid-body (neutrally buoyant, 6-DoF, two-way) | High: resolved 3D pipe, long migration transient | Equilibrium radius vs 0.6R at low Re_p; correct outward shift with increasing Re_p |

---

## Part 4 — Multiphase / porous / thermal disposition

| Benchmark | Reference data | Status | Note |
|-----------|----------------|--------|------|
| **Hysing rising bubble (2D)** | **Tabulated reference values** — rise velocity, circularity, center-of-mass vs time — published in Hysing et al. 2009 IJNMF 60:1259 and hosted at featflow.de benchmark repository (quantitative, directly usable) | NEEDS MF-γ | Shan-Chen cannot reach the density ratios (case 1: 10, case 2: 1000). Reserve as the **headline MF-γ acceptance case** — reference data quality is excellent, so it is the best-anchored two-phase target once MF-γ lands. |
| **Porous media — Darcy limit** | Sangani & Acrivos 1982 series (see Part 1 #2) | **READY NOW** (already ranked #2) | The Darcy/permeability limit is fully covered by the periodic-array item. |
| **Porous media — Ergun (packed bed)** | **Correlation** (Ergun equation) | PARTIAL / needs work | Geometrically resolvable (packed spheres) but inertial/turbulent and large; only the linear Darcy regime is a cheap ready target. Defer the full Ergun inertial branch. |
| **Turek-Hron FSI1-3 (flexible flag)** | Turek & Hron 2006 tabulated tip-displacement/force reference values | OUT OF SCOPE | Requires an elastic/deformable solid solver — not available, not planned. |
| **Dam break / sloshing** | Various experimental (Martin-Moyce, Koshizuka) | OUT OF SCOPE | Requires free-surface tracking — not available, not planned. |
| **Thermal benchmarks** (differentially heated cavity, Rayleigh-Bénard, natural convection Nu correlations) | de Vahl Davis tabulated; Nu correlations | OUT OF SCOPE | No energy equation / thermal LBM. Explicitly out of scope per task. |

---

## Recommended landing order

1. **Kovasznay, Sangani-Acrivos array, Womersley** — three cheap exact/series gates; land as permanent regression tests (convergence-order + L2 bands). Highest information-per-compute in the entire backlog.
2. **Impulsively started cylinder, backward-facing step** — moderate-cost separated/transient anchors; both figure-digitized, so budget a digitization step and keep bands generous.
3. **Oscillating cylinder** — once the prescribed-motion moving-boundary path is confirmed (penalization preferred over the in-flight IBM translation).
4. Open a **rigid-body track proposal (suggest MF-η)** whose acceptance ladder is Jeffery orbits (analytic, cheapest) → ten Cate (best PIV anchor) → DKT / Segré-Silberberg. This single track unlocks four benchmarks.
5. Reserve **Hysing** as the MF-γ headline case (best-quality two-phase reference data in the field).
6. Rushton/Costes-Couderc and Taylor-Couette wavy are READY but heavy — schedule against M-F compute budget, not the fast regression suite.

---

## Sources

- [Kovasznay flow — Wikipedia (closed-form solution + λ, Re relation)](https://en.wikipedia.org/wiki/Kovasznay_flow)
- [Womersley pulsatile LBM validation — Sci. Rep. (Nature)](https://www.nature.com/articles/s41598-022-05269-w) · [axisymmetric LBM pulsatile pipe (ResearchGate)](https://www.researchgate.net/publication/287726861_Simulation_of_pulsatile_flow_in_a_circular_pipe_using_an_axisymmetric_lattice_Boltzmann_method)
- [ten Cate et al. 2002 settling-sphere PIV + LBM (ResearchGate)](https://www.researchgate.net/publication/46014429_PIV_Experiments_and_Lattice-Boltzmann_Simulations_on_a_Single_Sphere_Settling_Under_Gravity)
- [Hysing et al. 2009 rising-bubble benchmark (Wiley IJNMF)](https://onlinelibrary.wiley.com/doi/10.1002/fld.1934) · [FeatFlow rising-bubble reference data](https://featflow.de/en/benchmarks/cfdbenchmarking/bubble.html) · [Hysing 2009 PDF (Karlin/Hron)](https://www.karlin.mff.cuni.cz/~hron/NMMO403/1934_ftp.pdf)
- [Armaly et al. 1983 backward-facing step (Scholars' Mine / JFM 127:473)](https://scholarsmine.mst.edu/mec_aereng_facwork/1047/)
- [Dütsch et al. 1998 oscillating cylinder KC (Journal of Hydrodynamics ref. context)](https://link.springer.com/article/10.1016/S1001-6058(11)60302-8)
- [Koumoutsakos & Leonard 1995 impulsively started cylinder — CSE-Lab data mirror](https://www.cse-lab.ethz.ch/kljfm1995/) · [JFM 296:1 (Cambridge Core)](https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/abs/highresolution-simulations-of-the-flow-around-an-impulsively-started-cylinder-using-vortex-methods/615FA46CB8A5BACF14FB4100A3CBC598)
- [Costes & Couderc 1988 Rushton turbine LDV — Reynolds-scaling validation context (PSU)](https://www.engr.psu.edu/ce/hydro/hill/publications/cespart1_05.pdf)
- [Segré & Silberberg 1962 inertial migration — DNS validation context (ResearchGate)](https://www.researchgate.net/publication/342925744_Direct_Numerical_Simulation_of_the_Segre-Silberberg_Effect_Using_Immersed_Boundary_Method)
- [Fortes, Joseph & Lundgren 1987 DKT — reproduction/validation context (ScienceDirect)](https://www.sciencedirect.com/science/article/abs/pii/S0045793014001042)
- [Sangani & Acrivos 1982 periodic cylinder-array permeability — validation context (ScienceDirect)](https://www.sciencedirect.com/science/article/pii/S0307904X14002868)

---

**Delivery note:** returned as final message per instructions (not written to disk). Ready to be committed as `docs/qa/benchmark-backlog.md`. All reference-data provenance is stated per item; no reference values were invented — where only figures exist, the item is flagged "digitize."
