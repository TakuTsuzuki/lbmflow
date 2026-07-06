# V&V benchmark backlog (V&V master plan lane 4.2)

Outstanding validation benchmarks the simulator does not yet run, with
reference-data provenance stated per item. Landed items dropped
(Kovasznay / Sangani-Acrivos periodic array / Womersley — commit 2c78d85).
Drop any row when its regression test lands.

**Capability legend**: READY (heavy) = shipped, large compute · IN-FLIGHT FIX
= blocked on marker-IBM moving/rotating fix (ANOM-P4-001 family) ·
NEEDS TRACK: rigid-body = new MF-η "resolved rigid body" (1–6 DoF two-way
coupling) proposed · NEEDS TRACK: two-way particle = two-way point/finite
particle coupling · NEEDS MF-γ = planned high-density-ratio phase-field
track · OUT OF SCOPE = elastic solids / free surface / energy equation.

## READY-NOW (moderate cost, digitized references)

| Rank | Benchmark | Compute | Validates | Reference | Metric |
|---|---|---|---|---|---|
| 1 | **Impulsively started cylinder (Koumoutsakos-Leonard)** | D≈60–100, T=Ut/D≈6, Re ∈ {550, 1000, 3000, 9500} | Transient drag on stationary curved (Bouzidi) boundary; wake / secondary-vortex topology | **Digitize** KL 1995 JFM 296:1 (also cse-lab KL_JFM_1995 mirror) | C_D(t) band vs KL at Re=550, 1000; secondary-vortex onset time |
| 2 | **Backward-facing step (Armaly)** | step S≈20, channel ~40S, steady, Re ≤ 400 | Separation / reattachment length | **Digitize** Armaly 1983 JFM 127:473 | x_r/S vs Armaly at Re ∈ {100, 200, 400}; **restrict Re ≲ 400** (3D above) |
| 3 | **Oscillating cylinder (Dütsch)** | Re=100, KC=5, several periods, prescribed in-line oscillation | Unsteady moving-boundary force + Morison decomposition — first prescribed-motion validation | **Digitize** Dütsch et al. 1998 JFM 360:249 | u,v profiles at phases {0°, 96°, 192°, 288°}; C_D history |

Route via volume penalization (already supports prescribed rotor motion → 
prescribed translation is a small extension). If forced through marker-IBM,
reclassify item 3 as IN-FLIGHT FIX.

## READY (heavy) — schedule against M-F budget

- **Taylor-Couette wavy-vortex regime** — critical Ta for onset (DiPrima &
  Swinney); Wendt/Lathrop G(Re) torque correlation. 3D annulus, O(10⁷)
  cells, long transient. Metric: Ta_c within few %; G(Re); azimuthal
  wavenumber.
- **Rushton turbine Np + Costes-Couderc field** — Np ≈ 5.0–5.5
  (correlation); radial-jet profiles from Costes & Couderc 1988 LDV
  (digitize). 3D tank, O(10⁷–10⁸), turbulent averaging. Metric: Np within
  ~15% at Re > 10⁴; radial-jet velocity/turbulence bands. Headline M-F
  acceptance case.

## Blocked on in-flight fix / new track

| Benchmark | Reference | Blocked on | Metric |
|---|---|---|---|
| Jeffery orbits (ellipsoid in shear) | **Analytic** T = (2π/γ̇)(r + 1/r) | rigid-body track (free rotational DoF; marker-IBM fix necessary but *not* sufficient — it fixes prescribed rotation) | Orbit period within few % across aspect ratios; drift-free ≥ 3 orbits |
| ten Cate single-sphere sedimentation | **Digitize** ten Cate 2002 Phys. Fluids PIV; 100×160×100 mm, D=15 mm, Re ∈ {1.5, 4.1, 11.6, 32.2} | rigid-body track (1-DoF two-way) | Settling velocity vs PIV all 4 Re; wall-approach deceleration |
| Drafting-kissing-tumbling | **Digitize** Fortes-Joseph-Lundgren 1987 JFM 177:467 | rigid-body track (6-DoF, 2 bodies, near-contact) | D→K→T sequence; kissing time & post-tumble separation |
| Segré-Silberberg inertial migration | **Analytic/experimental** eq. radius ≈ 0.6R (Segré-Silberberg 1962) | two-way particle **or** rigid-body (neutrally buoyant, 6-DoF, two-way) | Equilibrium radius vs 0.6R at low Re_p; outward shift with Re_p |

**Propose MF-η "resolved rigid body" track** — its acceptance ladder is
Jeffery (cheap analytic) → ten Cate (best PIV anchor) → DKT →
Segré-Silberberg. One track unlocks four benchmarks.

## Multiphase / porous / thermal disposition

- **Hysing rising bubble (2D)** — NEEDS MF-γ. Tabulated reference values
  (Hysing 2009 IJNMF 60:1259; featflow.de). Reserve as MF-γ headline case
  (best-quality two-phase data available; density ratios 10 and 1000
  unreachable with Shan-Chen).
- **Porous — Ergun** — geometrically resolvable but inertial/turbulent and
  large. Only the linear Darcy regime (covered by Sangani-Acrivos, already
  landed) is a cheap ready target.
- **Turek-Hron FSI1-3 (flexible flag)** — OUT OF SCOPE (elastic solid).
- **Dam break / sloshing** — OUT OF SCOPE (free surface).
- **Thermal (differentially heated cavity, Rayleigh-Bénard, Nu)** — OUT OF
  SCOPE (no energy equation / thermal LBM).

## Sources

- [Armaly 1983 backward-facing step (Scholars' Mine / JFM 127:473)](https://scholarsmine.mst.edu/mec_aereng_facwork/1047/)
- [Koumoutsakos-Leonard 1995 impulsively started cylinder (CSE-Lab mirror)](https://www.cse-lab.ethz.ch/kljfm1995/) · [JFM 296:1](https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/abs/highresolution-simulations-of-the-flow-around-an-impulsively-started-cylinder-using-vortex-methods/615FA46CB8A5BACF14FB4100A3CBC598)
- [Dütsch 1998 oscillating cylinder KC context](https://link.springer.com/article/10.1016/S1001-6058(11)60302-8)
- [Costes-Couderc 1988 Rushton LDV context (PSU)](https://www.engr.psu.edu/ce/hydro/hill/publications/cespart1_05.pdf)
- [ten Cate 2002 settling-sphere PIV + LBM](https://www.researchgate.net/publication/46014429_PIV_Experiments_and_Lattice-Boltzmann_Simulations_on_a_Single_Sphere_Settling_Under_Gravity)
- [Segré-Silberberg 1962 inertial migration — DNS context](https://www.researchgate.net/publication/342925744_Direct_Numerical_Simulation_of_the_Segre-Silberberg_Effect_Using_Immersed_Boundary_Method)
- [Fortes-Joseph-Lundgren 1987 DKT — reproduction context](https://www.sciencedirect.com/science/article/abs/pii/S0045793014001042)
- [Hysing 2009 rising bubble (Wiley IJNMF)](https://onlinelibrary.wiley.com/doi/10.1002/fld.1934) · [FeatFlow reference data](https://featflow.de/en/benchmarks/cfdbenchmarking/bubble.html)
