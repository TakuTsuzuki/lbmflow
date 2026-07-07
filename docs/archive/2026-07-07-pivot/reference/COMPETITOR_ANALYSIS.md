# Competitor Source Analysis — OpenFOAM / OpenLB / Palabos

> Reference material for LBMFlow ("流体sim"). Goal per user directive: LBMFlow must
> be a **practical strict superset ("上位互換")** of OpenFOAM, OpenLB, and Palabos.
> Produced by the research agent from **direct source reads** of shallow clones
> (OpenFOAM-dev @ 54064fa8, OpenLB 1.9.0 dev tree, Palabos). Each claim carries a
> file/class citation.
>
> **Status (2026-07-05):** OpenFOAM ✅ · OpenLB ✅ · Palabos ✅ — all three complete.

---

## Executive orientation (how to read this)

- **OpenLB** = the LBM feature-parity yardstick. Its moat is (1) multiplicative
  **dynamics composition** `Tuple<MOMENTA, EQUILIBRIUM, COLLISION, FORCING>`,
  (2) the **Bouzidi / STL / material-number geometry pipeline**, and (3) the
  **UnitConverter** SI UX. Closing OpenLB P1 items 1–8 makes the "strict superset
  of practically-used OpenLB" claim defensible.
- **OpenFOAM** = NOT an LBM code; do **not** port its FVM numerics. Its moat is a
  **product envelope**: dimensioned-SI I/O, checkpoint/restart, a composable
  runtime post-processing (function-object) framework, a time-varying BC semantic
  set, and a discoverable runtime-selectable model registry — all expressed as a
  self-describing, reproducible text-case bundle. LBMFlow must reproduce that
  envelope in its JSON-scenario / CLI / MCP idiom.
- **Where LBMFlow is already ahead** (keep and market): 200+ adversarial validation
  suite, MCP/agent-native async API, WASM web GUI, f32 **deviation storage (f−w)**,
  and raw per-node speed (1,480 MLUPS 2D / 6–12 GLUPS Apple GPU).
- **The single highest-leverage gap** flagged by both codes and by the user:
  **an SI unit system** at the API boundary (currently raw lattice units).

---

# OpenLB findings

Analyzed tree: OpenLB **1.9.0** dev tree (`CITATION.cff` version 1.9.0; `rules.mk`
`OLB_RELEASE := 1.8r1`), C++20, GPL-2. Scale check (verified directly): **139 example
Makefiles across 21 example categories** (`examples/{laminar,turbulence,thermal,
multiComponent,freeSurface,particles,fsi,optimization,porousMedia,adsorption,
advectionDiffusionReaction,radiativeTransport,electroChemistry,microfluidics,
gridRefinement,solid,noise,uncertaintyQuantification,pdeSolverEoc,showCase,web,
forBeginners}`). Platforms: `Platform::{CPU_SISD, CPU_SIMD, GPU_CUDA, GPU_HIP}`
(`src/core/platform/platform.h`). Descriptors: D2Q5/D2Q8/D2Q9, D3Q7/D3Q13/D3Q15/
D3Q19/D3Q27 (`src/descriptor/definition/common.h`). Even `examples/web/cylinder2d`
is an emscripten/WASM build of the full engine.

## A. Feature inventory table

| Feature | OpenLB status (source ref) | LBMFlow status | Gap? |
|---|---|---|---|
| **Lattices** | D2Q5/8/9, D3Q7/13/15/19/27 (`src/descriptor/definition/common.h`) + MRT/CUM/RTLBM descriptor extensions | D2Q9, D3Q19 | **Y** (D3Q27, ADE lattices D2Q5/D3Q7) |
| BGK / TRT | `collision::BGK`, `collision::TRT` with runtime `MAGIC` field (`src/dynamics/collision.h:382-411` — magic is tunable) | BGK, TRT (fixed Λ=3/16) | minor (tunable magic) |
| MRT | `collision::MRT` + per-moment `s<DESCRIPTOR>(i)` (`src/dynamics/collisionMRT.h`) | none | **Y** |
| Cumulant | `collision::CUM`, D3Q27 Chimera transform (`src/dynamics/collisionCUM.h`, `cumulantDynamics.h`) | planned | **Y** |
| KBC / entropic | `collision::KBC` — **D3Q19 only, static_assert** (`src/dynamics/collisionKBC.h:46`); `EntropicEqDynamics` | none | **Y** |
| Regularized | `collision::RLB/IRLB/CIRLB/ThirdOrderRLB`; HRR hybrid recursive regularized (`collisionHRR.h`) | none | **Y** |
| LES | Wrapper collisions in `src/dynamics/collisionLES.h`: `SmagorinskyEffectiveOmega<COLLISION>`, `ShearSmagorinsky…`, `ConStrainSmagorinsky…`, `WaleEffectiveOmega`, `LocalVanDriestSmagorinsky…`; + `stochasticSGSdynamics.h`, ADM functors | none (planned) | **Y** |
| RANS | `kEpsilonRANSDynamics.h` | none | Y (P3) |
| Wall functions | Spalding iterative + Van Driest damping (`src/boundary/setTurbulentWallModel.h`, `boundary/postprocessor/turbulentWallModelPostProcessor.h`) | none | **Y** |
| Turbulent inlet | Vortex method synthetic turbulence (`src/boundary/vortexMethod.h`) | none | Y |
| Non-Newtonian | Power-law + Herschel-Bulkley w/ yield stress via `powerlaw::OmegaFromCell` (`src/dynamics/powerLawBGKdynamics.h`), dedicated `PowerLawUnitConverter` | planned μ(γ̇) | **Y** |
| Forcing schemes | Guo, PlainGuo, AdeGuo, Liang, CorrectedLiang, LiangTRT, KupershtockhLiang, MCGuo, LaddVerberg (`src/dynamics/forcing.h`) | Guo only | Y (Liang family for phase-field) |
| Zou-He | `zouHeDynamics.h`, `zouHeVelocity/Pressure{2D,3D}.h` | 2D all edges + D3Q19 faces | N |
| Regularized local BC (Latt) | `localVelocity.h`/`localPressure.h` via `momenta::Regularized*BoundaryTuple` | none | **Y** |
| FD interpolated BC (Skordos) | `interpolatedVelocity/Pressure*.h` + `StraightFdBoundaryProcessor` | none | Y |
| **Curved boundary** | Bouzidi + Yu IBB + velocity variant + ADE variant, per-link `BOUZIDI_DISTANCE`, distance from analytic/STL, q=0.5 fallback (`src/boundary/setBouzidiBoundary.h`, `bouzidiFields.h`) — **dominant BC in examples** (aorta3d, cylinder2d, airfoil2d…) | none | **Y (flagship gap)** |
| Slip / partial slip | `slip{2D,3D}.h`, `partialSlip*.h` (TUNER α∈[0,1]) | none | **Y** |
| Characteristic (non-reflective) outlet | `characteristicBoundary.h` (Wissocq 2017) + sponge layer (`dynamics/spongeLayerDynamics.h`) | convective outflow only | **Y** |
| Convective outflow | `interpolatedConvection*.h` | yes | N |
| ADE/thermal BCs | ADE Dirichlet (Allen-Reis), `regularizedTemperature*.h`, `regularizedHeatFlux*.h`, Robin (`robin.h`) | none | **Y** |
| BC framework | `boundary::set<BC>(lattice, geometry, material)` + discrete-normal classification + momenta-swap unification (`src/boundary/setBoundary.h`) | per-edge BC enum | Y (architecture) |
| Shan-Chen SCMP/MCMP | SC93/94, Carnahan-Starling, Peng-Robinson EOS; N-component PR w/ Huron-Vidal mixing (`src/dynamics/interactionPotential.h`) | SC SCMP + 2-component (2D) | Y (EOS zoo, 3D) |
| Free energy | 2/3-component Swift-style: `collision::FreeEnergy`, wetting wall BCs (`freeEnergyDynamics.h`, `boundary/phaseFieldWall.h`) | none | Y |
| Phase field (Allen-Cahn) | Incompressible phase-field stack: `phaseFieldCoupling.h`, IRLB/CIRLB, Liang forcing, `phaseFieldInletOutlet.h` | planned (conservative AC 10³) | Y (planned) |
| **Free surface (VOF-style)** | Full module: mass tracking, interface reconstruction, cell conversion (`src/dynamics/freeSurfacePostProcessor{2D,3D}.h`), 7 examples | none | **Y** |
| Thermal / Boussinesq | Coupled ADE lattice + `SmagorinskyBoussinesqCoupling`, `TotalEnthalpyPhaseChangeCoupling` (melting), `ThermalUnitConverter` | none | **Y** |
| Reaction / adsorption | ADE reaction coupling, explicit-FD reaction, Langmuir/Freundlich isotherms (`src/reaction/`) | none | Y (M-F reactor relevance) |
| Radiative transport | RTLBM P1 (`rtlbmDynamics.h`, `radiativeUnitConverter.h`) | none | Y (P3) |
| **Grid refinement** | Lagrava 2012 vertex-centered + Rohde cell-centered (`src/refinement/algorithm/{lagrava,rohde}.h`), 2D+3D | uniform only | **Y** |
| Sparse domains | `CuboidDecomposition::{removeByWeight, shrink, splitByWeight}` (`src/geometry/cuboidDecomposition.h`) | dense grid | **Y** |
| **STL / CAD pipeline** | `STLreader<T> : IndicatorF3D` octree accel, 4 ray-cast modes, signed distance; CSG indicators; `SuperGeometry::rename()` (`src/io/stlReader.h`, `geometry/superGeometry.h`) | JSON primitives (circle/rect/sphere) | **Y** |
| **Units (SI)** | `UnitConverter` family (`src/core/unitConverter.h`) + thermal/ADE/multiphase/power-law/radiative converters | raw lattice units | **Y (flagship gap)** |
| MPI | Cuboid decomposition + `HeuristicLoadBalancer`, `HeterogeneousLoadBalancer` (CPU+GPU mix), `SuperCommunicator::requestField<FIELD>()` | halo exchange, bit-identical partition invariance, 97-99% weak scaling | partial (load balancing, field-selective comm) |
| GPU | CUDA + HIP/ROCm (`src/core/platform/gpu/{cuda,hip}/`) | wgpu (Metal/Vulkan/DX12 — arguably broader) | N (different trade) |
| SIMD CPU | AVX2/512 columns + `CyclicColumn` periodic-shift streaming | fused collide-stream 1,480 MLUPS | N |
| VTK output | Parallel `SuperVTMwriter2D/3D` (per-cuboid VTI + VTM + PVD) | legacy VTK, no parallel | Y |
| **Checkpoint/restart** | `Serializer`/`Serializable`; `SuperLattice` `save()/load()` (`src/core/serializer.h`) | none | **Y** |
| Functor framework | `GenericF` hierarchy: `SuperLatticePhysVelocity3D`, `SuperLatticePhysDrag`, `SuperPlaneIntegralFluxVelocity`, `SuperRelativeErrorLpNorm`, functor arithmetic | probes/force CSV, stress fields | Y (composability) |
| Wall shear stress | `SuperLatticePhysWallShearStress3D` | strain-rate field exists | minor |
| Solver/Case framework | `BaseSolver` (`src/solver/lbSolver.h`) + 2025 declarative `ConcreteCase<MAP>` (`src/case/case.h`) | scenario JSON runner (arguably ahead) | N |
| Particles (resolved/HLBM) | Momentum-exchange resolved particles (Wen 2014 / Ladd 1994) | none | Y |
| Particles (subgrid + DEM) | Verlet dynamics, overlap-volume contact, MPI particle migration | none | Y |
| FSI / IBM | `InterpolateVelocityO<WIDTH>`/`SpreadForceO<WIDTH>` + deformable membranes (`src/fsi/ibm/`) | none | Y (P3) |
| Adjoint / optimization | Dual-LBM adjoint, LBFGS/BFGS, topology opt, AD tape (`src/optimization/`, `utilities/aDiff.h`) | none | Y (P3) |
| UQ | gPC, MC, Sobol QMC, LHS (`src/uq/`) | none | Y (P3) |
| MCP / agent API, web GUI | none (only emscripten demo) | MCP server, TS GUI, HTML gallery | **N — LBMFlow ahead** |
| Adversarial validation suite | example-level EOC tests (`examples/pdeSolverEoc/`) | 200+ adversarial tests | **N — LBMFlow ahead** |
| Deviation storage (f−w) f32 | not present | yes | **N — LBMFlow ahead** |

## B. Top implementation techniques worth adopting

**B1. Composable dynamics: `Tuple<MOMENTA, EQUILIBRIUM, COLLISION, FORCING>`**
(`src/dynamics/interface.h`, `dynamics.h`, `momenta/aliases.h`). Every dynamics is a
4-slot composition; e.g. `SmagorinskyForcedBGKdynamics = Tuple<BulkTuple, SecondOrder,
SmagorinskyEffectiveOmega<BGK>, Guo<ForcedWithStress>>`. The collision-*wrapper* trick
is the key: LES, per-cell omega (`collision::OmegaFromCell`), porous, TRT-magic are all
decorators around an inner collision, so N collisions × M modifiers multiply without
duplication.
*Rust mapping:* `trait Collision<L: Lattice> { fn collide(&self, m: &Moments, feq: &[f32],
f: &mut [f32], p: &CellParams); }` with wrapper structs `Smagorinsky<C: Collision>`,
`OmegaFromField<C>`. Monomorphization gives zero-cost composition inside the existing
`step_band` fused kernel; for WGSL, mirror as shader-permutation generation keyed by the
same (collision, modifier, forcing) triple, gated by T14 CPU/GPU equivalence.

**B2. UnitConverter — adopt the API shape nearly verbatim** (`src/core/unitConverter.h`).
9 conversion factors (`_conversionLength/Time/Velocity/Density/Mass/Viscosity/Force/
Torque/Pressure`), char values (`_charPhysLength/Velocity, _physViscosity, _physDensity,
_charPhysPressure`), discretization (`_resolution, _latticeRelaxationTime,
_charLatticeVelocity`). Key derivations: `τ = ν_phys/(dx²/dt)·invCs2 + 0.5`;
`conversionForce = ρ·dx⁴/dt²`; `conversionPressure = conversionForce/dx²`;
`ρ_lat = (p−p_ref)·invCs2/conversionPressure + 1`. Pairwise `getPhysX/getLatticeX` for
every quantity; derived `getReynoldsNumber/getMachNumber/getKnudsenNumber`. **Three named
constructors** matching how users think: FromResolutionAndRelaxationTime (N, τ),
FromResolutionAndLatticeVelocity (N, u_lat), FromRelaxationTimeAndLatticeVelocity (τ,
u_lat). `print()` warns if `u_lat > 8(τ−0.5)`.
*Rust mapping:* one `UnitConverter` struct + three constructors, embedded as a `units`
block in scenario JSON; `lbm validate` reprints the table and reuses existing tau/Mach/
grid-Re warnings (tune-stability Skill levers map 1:1). Highest commercial value, low risk.

**B3. Bouzidi curved boundary as a per-link distance field**
(`src/boundary/setBouzidiBoundary.h`, `bouzidiFields.h`). Wall distance q∈[0,1] stored per
population direction (`BOUZIDI_DISTANCE`, −1 = inactive); setup ray-casts each link against
analytic/STL indicator, **falls back to q=0.5 = plain halfway BB on failure** (degrades
gracefully to LBMFlow's current BC). Post-stream two-branch interpolation; velocity variant
adds `veloCoeff·t_i·invCs2` for moving walls.
*Rust mapping:* optional `Option<Box<[f32]>>` q-major field with the same
`f[q*plane+cell]` layout as f; populated at scenario build (exact distance for circle/
sphere/rect is analytic — no STL needed initially); a post-stream pass over a boundary-cell
index list. Moves cylinder/sphere drag from staircase to 2nd order — highest-leverage
physics upgrade.

**B4. Momenta-swap BC unification + `boundary::set(lattice, geometry, material)`**
(`src/boundary/setBoundary.h`). One BC parameterized by a momenta bundle:
`RegularizedVelocityBoundaryTuple` (fix u, compute ρ) vs `RegularizedPressureBoundaryTuple`
(fix ρ, compute u) reuse the same code. The setter classifies each boundary cell's discrete
normal (flat/edge/corner) and installs dynamics per cell.
*Rust mapping:* refactor Zou-He to a `BoundaryMoments` trait + normal-classification pass;
3D edges/corners (current weak spot) fall out of the same path.

**B5. LES as effective-omega computed in-collide from f_neq** (`src/dynamics/collisionLES.h`).
`τ_eff = τ + ½(√(τ² + C·C_s²·‖Π_neq‖) − τ)` from the non-equilibrium stress already inside
the collision — no extra pass for plain Smagorinsky. WALE needs one velocity-gradient
stencil but is the wall-correct engineering default.
*Rust mapping:* wrapper collision (B1) reading Π_neq inside `step_band`; zero extra memory
traffic. Cheapest credible path to turbulence for the M-F stirred-reactor milestone.

**B6. Material-number geometry pipeline** (`src/geometry/superGeometry.h`,
`src/io/stlReader.h`). Recipe: voxelize STL (octree + ray casting) → integer material field
→ `rename(from, to, indicator/neighbor)` + `clean()/checkForErrors()` → each BC bound to a
material number. Geometry prep decoupled from physics.
*Rust mapping:* formalize an explicit `material: u8` grid with `rename`-style ops in the
scenario builder; STL import (pure-Rust `stl_io` + 3-ray majority vote) becomes additive.

**B7. Weighted sparse cuboid decomposition** (`src/geometry/cuboidDecomposition.h`,
`communication/heuristicLoadBalancer.h`). Domain = many small cuboids each carrying
`_weight` = fluid-cell count; `removeByWeight()` discards solid blocks; balancer cost =
`full + ratioFullEmpty·empty`.
*Rust mapping:* generalize the MPI partitioner from "slabs of one box" to "list of weighted
boxes"; T13 partition-invariance suite is the regression harness. Prerequisite for
competitive STL-geometry runs.

**B8. Serializer checkpoint framework** (`src/core/serializer.h`). Uniform `Serializable`
trait over columns/fields/geometry; `SuperLattice::save()/load()` per rank.
*Rust mapping:* serde/bincode over the SoA buffers + a manifest of (descriptor, step, RNG,
halo). Small effort, table-stakes commercial feature.

**B9. Lagrava vertex-centered refinement** (`src/refinement/algorithm/lagrava.h`).
Coarse↔fine coupling: interpolate ρ,u in space/time, rescale f_neq by `(τ_c − 0.25)/τ_c`,
with dx,dt both /2 (acoustic scaling). Prefer Lagrava over Rohde (better documented,
examples-backed).
*Rust mapping:* two lattices + a coupling operator on interface cell lists; between-step
coupler. Large effort — schedule after B1–B6.

## C. Prioritized gap list

**P1 — must-have to claim LBM superiority:**
1. **SI unit conversion (B2)** — S/M. Pure API work, no solver risk.
2. **Bouzidi interpolated bounce-back (B3)** — M. Analytic primitives first, STL later.
   Gate: Schäfer-Turek drag/lift vs current staircase.
3. **MRT + D3Q27 + cumulant collision** — M/L. MRT is S/M once B1 lands.
4. **Smagorinsky + WALE LES (B5)** — M. Validates on existing TGV3D reference.
5. **Slip/free-slip + symmetry BC** — S. Local; embarrassing to lack.
6. **Checkpoint/restart (B8)** — S/M. Commercial table stakes.
7. **Non-Newtonian power-law/Herschel-Bulkley** — S/M once per-cell-omega (B1) exists.
8. **Dynamics composition refactor (B1)** — M. The multiplier that makes 3,4,7 cheap.

**P2 — should-have:** STL import + material pipeline (B6, M/L); thermal ADE lattice +
Boussinesq (M); sparse weighted decomposition (B7, M); characteristic outlet + sponge (M);
grid refinement Lagrava (B9, L); free-surface VOF (L); regularized local BC + tunable TRT
magic (S each); parallel VTK + error-norm functors (M); wall functions after LES (M);
reaction/ADE coupling (M/L, REQ_STIRRED_REACTOR-driven); KBC/entropic (M — a D2Q9+D3Q19 KBC
would *exceed* OpenLB's D3Q19-only).

**P3 — optional/defer:** subgrid+resolved particles (L); IBM/FSI membranes (L); adjoint/
optimization+AD (L); UQ (M, niche); radiative/electrochem/elastodynamics/Knudsen (L each);
k-ε RANS (M); vortex-method turbulent inlet (S/M, with LES).

## D. Explicitly NOT worth adopting

- **CSE codegen pipeline** (`src/cse/`, SymPy/Mako-generated hash-named kernels) — exists to
  defeat C++ template-tree optimizer failures; rustc/LLVM already CSE monomorphized code, and
  measured MLUPS says kernels aren't the bottleneck.
- **2D/3D `.h`/`.hh` file duplication** (×~150 files + `legacy/` trees). LBMFlow's
  dimension-generic core is strictly better; OpenLB's own 2024-25 rewrites are them escaping it.
- **Runtime platform dispatch via virtual base + switch** (`callUsingConcretePlatform`) — Rust
  enums/generics + backend trait already do this without type-erasure hazards.
- **XML case files** (`src/io/xmlReader.h`) — self-describing JSON schema is a generation ahead;
  adopt the *converter*, not its XML front end.
- **CyclicColumn periodic-shift streaming** — conflicts with the fused collide-stream band
  kernel and wgpu coalescing; would re-open every bit-equivalence gate (T13/T14) for no gain.
- **GIF/Gnuplot writers**, **Boolakee elastodynamics BCs / OSM parser** (research one-offs),
  **OpenLB's OpenMP layer** (rayon covers it more safely).

**Bottom line:** OpenLB's moat = dynamics composition + Bouzidi/STL/material pipeline +
UnitConverter UX. LBMFlow is already ahead on validation rigor, agent/MCP UX, WASM GUI, f32
deviation storage, and per-node speed. Real OpenLB users overwhelmingly run Smagorinsky/
WALE-BGK + Bouzidi + interpolated velocity/pressure BC + UnitConverter + VTM output — so P1
items 1–8 are exactly the parity list.

---

# OpenFOAM findings

Source tree: OpenFOAM-dev commit `54064fa8` (2026-07-04), modular-solver layout (`foamRun`
+ solver modules under `applications/modules/`). 260 tutorial cases
(`find tutorials -path '*/system/controlDict' | wc -l`).

The single most important structural fact: OpenFOAM-dev is a **module-selectable solver**.
`applications/solvers/foamRun/foamRun.C` loads a solver module named by `solver <name>;` in
`controlDict` (`tutorials/incompressibleFluid/pitzDaily/system/controlDict`:
`solver incompressibleFluid;`). Everything else — schemes, models, BCs, function objects —
is runtime-selected from dictionaries. The moat is not the FVM numerics; it's that a case is
a **self-describing, reproducible, text-dictionary bundle** plus a **huge runtime-selectable
model/BC/post-processing library** with **dimensional-unit safety** in every field.

## A. Workflow-capability table

| Capability | OpenFOAM mechanism (source ref) | LBMFlow status | Gap? |
|---|---|---|---|
| Self-describing case bundle | `system/{controlDict,fvSchemes,fvSolution}` + `constant/` + `0/`; each field carries `FoamFile{class;object;}` + `dimensions` + `boundaryField` (`pitzDaily/0/U`) | JSON scenario schema (single file) | Partial |
| Runtime model selection | `runTimeSelectionTable` macros (`src/OpenFOAM/db/runTimeSelection/`); `type <name>;`; mismatch prints valid-choice list | CLI presets + scenario enums | Partial (narrower registry) |
| Dict language: includes/macros/coded | `functionEntries/`: `includeEntry`, `includeFuncEntry`, `codeStream`, `calcEntry`, `ifEntry/ifeqEntry`, `$var` expansion | None (flat JSON) | Y |
| **Dimensional units (SI everywhere)** | `dimensionSet` (7 base dims, `src/OpenFOAM/dimensionSet/dimensionSet.H:135`); arithmetic checks dims, `FatalError "LHS and RHS of + have different dimensions"` (`dimensionSet.C:575`); named dims `dimensions [velocity];` | **Lattice units only** | **Y — #1 gap** |
| Runtime post-processing (function objects) | `src/functionObjects/{field,forces,solvers,utilities}` + `src/sampling`; base contract `execute()/write()/end()` w/ per-object `executeControl/writeControl` | PNG/CSV/VTK, probes, force series | Partial (has primitives, lacks framework+breadth) |
| Post-hoc processing of saved results | `foamPostProcess -func <name>` over saved time dirs | Outputs at run time only | Y |
| Turbulence: LES + wall functions | `src/MomentumTransportModels/.../LES/` (Smagorinsky, WALE, kEqn, dynamicKEqn, dynamicLagrangian, Deardorff, SpalartAllmarasDES/DDES/IDDES, kOmegaSSTDES); wall functions `nut*/epsilon*/omega*/kqR*` | None (planned) | Y (roadmap) |
| BC semantics library | `src/finiteVolume/fields/fvPatchFields/{derived,constraint}` (~80 derived + 13 constraint) | BB + Zou-He + outflow | Y |
| Time-varying inputs | `Function1`: Table, CSV, Sine, Square, Ramp variants, Polynomial, Coded, Scale | None | Y |
| Adaptive dt (Courant) | `adjustTimeStep yes; maxCo; maxDeltaT` (pitzDaily: `maxCo 5;`) | Fixed lattice dt | N (different by design) |
| Checkpoint / restart | `startFrom latestTime;` — every time dir is a full checkpoint; `runTimeModifiable yes` re-reads dicts mid-run | **None** | **Y — P1** |
| Write control | `writeControl {timeStep\|runTime\|adjustableRunTime\|cpuTime\|clockTime}`, `writeInterval`, `purgeWrite`, `writeFormat`, `writeCompression` | Fixed cadence | Partial |
| Non-Newtonian viscosity | `generalisedNewtonian/.../strainRateViscosityModels/`: powerLaw, CrossPowerLaw, BirdCarreau, HerschelBulkley, Casson | Newtonian only | Y (P2, LBM-reachable) |
| Passive/active scalar transport | `functionObjects/solvers/scalarTransport`, `age`, `phaseScalarTransport` | None (planned) | Y (P1) |
| Buoyancy / heat | `fvModels/general/{buoyancyForce,buoyancyEnergy,heatSource}`; thermal solver modules | None | Y (P2) |
| Porous media | `src/fvModels/general/porosityForce/` (Darcy-Forchheimer) | None | Y (P2, LBM: IBM/greyscale) |
| MRF (rotating frames) | `src/finiteVolume/cfdTools/general/MRF/MRFZone*` | None | Different by design (LBM → moving/IBM) |
| Momentum source / channel driving | `fvConstraints/meanVelocityForce`; `fvModels/general/{acceleration,semiImplicitSource}` | Guo body force | Partial |
| Meshing: CAD→mesh | blockMesh + snappyHexMesh (STL→castellate→snap→layers) | Uniform Cartesian voxels | Different (LBM: voxelization) |
| Parallel decomposition | `decomposePar`/`reconstructPar`; scotch/ptscotch/metis/parMetis/zoltan/hierarchical/simple/manual/multiLevel/structured | MPI halo exchange | Partial |
| Collated parallel I/O | uncollated/collated/masterUncollated/hostCollated fileHandlers | N/A | N (Cartesian I/O simpler) |
| ParaView / VTK ecosystem | `bin/paraFoam`, `foamToVTK`, `foamToEnsight` | VTK output | N — covered |
| Dict-from-CLI / scripting | `foamDictionary`, `foamListTimes`, `foamMonitor`, `RunFunctions` | CLI + MCP (7 tools) | Partial — MCP arguably ahead |

## B. Product concepts LBMFlow must adopt

1. **A real dimensional-unit system (SI in, lattice internal).** Every input carries a unit;
   incompatible arithmetic is refused at read time — catches the #1 CFD setup error.
   Source: `src/OpenFOAM/dimensionSet/dimensionSet.H` (7-exponent array), `dimensionSet.C:575`,
   named dims in `dimensions.C:90`. *Mapping:* SI-input layer in scenario JSON
   (`{"value":1e-6,"unit":"m2/s"}` or a `dimensions:[...]` vector); a converter computes
   lattice tau/dt/dx/u from characteristic length/velocity/viscosity + resolution and **echoes
   back Re/Ma/tau + stability verdict** (aligns with tune-stability Skill). Enforce consistency
   in `lbm validate` with OpenFOAM-style mismatch errors. Highest-leverage single feature.
2. **Checkpoint/restart via "latest state is a resumable case."** `startFrom latestTime` makes
   any written time dir a complete restart point. *Mapping:* serialize full `f` (q-major SoA +
   halo) + step + scenario hash; JSON gains `startFrom` + `checkpointInterval`; MCP `start_run`
   gains `resume_from`. P1, mechanically simple on a fixed grid.
3. **Composable function-object framework (runtime + post-hoc).** Uniform, individually-
   scheduled, attachable observers; re-runnable over saved results via `foamPostProcess -func`.
   Source: `functionObject.H`; `functionObjects/field/` (fieldAverage, CourantNo, Q, Lambda2,
   vorticity, enstrophy, wallShearStress, yPlus, streamlines, vol/surfaceFieldValue);
   `functionObjects/forces/`; `src/sampling/` (probes; sampledSet line/arc/box; sampledSurface
   cutPlane/isoSurface; writers csv/vtk/ensight). *Mapping:* generalize `outputs`/`probes` into
   a scheduled **observers** list, each `type`+`writeInterval`+region. Priority observers:
   fieldAverage (mean+rms), forces/forceCoeffs, probes, line/plane sampling, vorticity/Q/λ2/
   enstrophy, wallShearStress+yPlus (post-LES), volume min/max/integrate, residual monitor.
   Add `lbm postprocess --func <name>` over an existing output dir; expose the catalog via MCP.
4. **Time-varying / tabulated inputs (`Function1`).** Source:
   `src/OpenFOAM/primitives/functions/Function1/` (Constant/Table/CSV/Sine/Square/Polynomial/
   Ramp/Scale/Coded). *Mapping:* any scalar/vector input may be a constant or
   `{"function":"sine"|"table"|"ramp",...}` evaluated per lattice step. Small, high-value.
5. **Richer BC semantics (not FVM math).** The 7 that matter: (1) time/space-varying Dirichlet
   → Zou-He/regularized; (2) **backflow-safe outlet (inletOutlet)** → mixed outflow BC; (3)
   total-pressure inlet; (4) symmetry + slip planes (periodic already present); (5) coded BC
   hook; (6) synthetic-turbulence inlet (co-req LES); (7) mapped/recycled inlet. Source:
   `fvPatchFields/derived` (~80), `constraint`. P1: inletOutlet, symmetry/slip, time-varying
   Dirichlet. P2: totalPressure, coded, synthetic-turbulence.
6. **Discoverable runtime-selectable registry (runTimeSelectionTable idea).** `type kOmegaSST;`
   with an error listing valid choices. *Mapping:* make `lbm schema` (and an MCP tool) emit the
   full catalog of collision schemes / BC types / observers / forcing models with params.
7. **Case/dict composition (`#include`/`#includeFunc`/`$var`).** For sweeps and DRY setup.
   *Mapping:* JSON `$include`, a `$params` block with `${var}` substitution, and
   `lbm run --set nu=1e-5,2e-5`; pairs with the MCP `start_run` sweep pattern.

## C. Prioritized gap list

**P1 — must-have for "practical superset of OpenFOAM (LBM-reachable flows)":**
- **SI unit system** (dimensioned inputs → lattice conversion + echoed Re/Ma/tau + validation
  errors). **L** — user's stated #1 gap.
- **Checkpoint/restart** (`startFrom latestTime`, checkpointInterval, MCP resume). **M**.
- **Scalar transport** (passive advection-diffusion; planned). **M**.
- **Observer framework** (fieldAverage, probes, line/plane sampling, vorticity/Q/λ2,
  forceCoeffs, residual monitor) + `lbm postprocess --func`. **L** (incremental).
- **BC semantics**: inletOutlet, symmetry/slip, time-varying Dirichlet (Function1 subset). **M**.

**P2:** LES + wall functions (**L**, co-req synthetic-turbulence inlet + yPlus/wallShearStress
observers); non-Newtonian viscosity (**M**); Boussinesq buoyancy + heat (**L**); porous media
(**M**); totalPressure/coded BC (**M**); Function1 full set + JSON include/macro for sweeps (**M**).

**P3:** Ensight export, live residual plotting (**S**, ~covered by web GUI); `foamDictionary`-
style CLI (**S**, ~MCP); synthetic-turbulence generators DFSEM-class (**M**, post-LES);
recycled/mapped inlet (**M**); user-facing decomposition method choice (**S**, low value).

## D. Explicitly out of scope for an LBM code (honest public one-liners)

- **Unstructured/body-fitted meshing (snappyHexMesh).** "LBMFlow uses a uniform Cartesian
  lattice by design; geometry enters by voxelization of STL/CAD, a deliberate trade of
  body-fitted accuracy for GPU-scalable regularity." (Equivalent bar it must still meet:
  STL → signed-distance voxelization with sub-voxel wall placement + local refinement.)
- **Adaptive Courant time stepping.** "LBM advances on a fixed lattice step tied to the grid;
  stability is governed by tau and Mach number, so LBMFlow tunes resolution/tau up front."
- **MRF for steady rotating machinery.** "LBMFlow models rotating geometry with moving/immersed
  boundaries (transient) rather than a frozen-rotor reference-frame source term."
- **Sliding-mesh / non-conformal AMI.** "On a fixed Cartesian lattice, relative motion is
  handled by immersed/moving boundaries, not mesh sliding."
- **Compressible/combustion/Eulerian-multiphase FVM solvers, Lagrangian spray.** "LBMFlow
  targets the incompressible-to-weakly-compressible, isothermal-to-buoyant flow classes LBM
  resolves well; high-Mach compressible, detailed combustion, and dense Euler-Euler are outside
  the LBM regime and not claimed."
- **Collated parallel I/O / graph-partition zoo.** "Cartesian domains partition trivially and
  near-optimally, so no graph-partitioning heuristics or collated reconstruction are needed."

**Framing for the design doc:** LBMFlow need not replicate OpenFOAM's *numerics library*; it
must replicate OpenFOAM's **product envelope** — dimensioned SI I/O, restart, a composable
observer framework, a time-varying BC semantic set, a discoverable model registry — in the
JSON-scenario/CLI/MCP idiom. The P1 list is the minimum bar for a defensible "上位互換" claim;
the unit system is both the largest single item and the most credibility-critical.

---

# Palabos findings

Source analyzed: shallow clone at `scratchpad/palabos` (paths relative to repo root). All claims
traced to actual source. **Licensing note: Palabos is AGPL-3.0 — techniques and designs may be
adopted, code must NOT be copied. Implement from the published papers (Guo 2002, Bouzidi 2001,
Coreixas 2019, Malaspinas LES, Inamuro IBM), not from Palabos sources.** (Verified firsthand that
`dco/` is a Developer-Certificate-of-Origin license file, not autodiff; `dsl2d/dsl3d` are
double-shear-layer collision benchmarks, not a DSL.)

## A. Feature inventory table

| Feature | Palabos status (source ref) | LBMFlow status | Gap? |
|---|---|---|---|
| **Lattices** | D2Q9, D3Q13/15/19/27 (`src/latticeBoltzmann/nearestNeighborLattices2D/3D.h`); extended D2Q37/D3Q39/D3Q121; ADE D2Q5/D3Q7/D3Q19/D3Q27 (`advectionDiffusionLattices.h`) | D2Q9, D3Q19 | **Y** (D3Q27, ADE lattices) |
| **BGK family** | BGK, IncBGK, QuasiIncBGK, CompleteBGK, Regularized/SecuredRegularized/CompleteRegularized BGK, ConstRho, Stochastic (`src/basicDynamics/isoThermalDynamics.h`) | BGK only | **Y** (incompressible + regularized variants) |
| **TRT** | TRT/IncTRT/Ma1TRT w/ runtime `setMagicParam()` (`src/complexDynamics/trtDynamics.h`); Complete/CompleteRegularized/Truncated TRT | TRT fixed Λ=3/16 | Partial (tunable magic, complete/reg variants) |
| **MRT** | MRT, IncMRT (`complexDynamics/mrtDynamics.h`); M/invM/S per descriptor (`mrtLattices.h`); Guo-forced MRT | None | **Y** |
| **Central-moment/cumulant ("comprehensive")** | RM, HM, CM, CHM, K (cumulant), GH, RR (`src/basicDynamics/comprehensiveIsoThermalDynamics.h`, `latticeBoltzmann/comprehensiveModelsTemplates*.h`, per Coreixas et al.); XML-selectable in `showCases/dsl2d,dsl3d`, GPU-exercised | Planned | **Y** (P1) |
| **Entropic / KBC** | EntropicDynamics, ForcedEntropic, VariableOmegaELBM (`entropicDynamics.h`); KBCDynamics (`kbcDynamics.h`) | None | Y (P3) |
| **Forcing** | Guo, Shan-Chen shift, He, naive (`basicDynamics/externalForceDynamics.h`) | Guo only, per-cell force field | Partial (He/SC-shift only for multiphase) |
| **LES** | ~20 Smagorinsky variants incl. Malaspinas **consistent** formulation (`complexDynamics/smagorinskyDynamics.h`, 694 lines), dynamic Smagorinsky (precomputed omega field), Smagorinsky-Boussinesq thermal LES. No WALE, no wall functions | None (planned) | **Y** (P1); WALE would exceed Palabos |
| **Non-Newtonian** | Carreau: CarreauDynamics/BGKCarreau/RegularizedBGKCarreau + `carreauUnits.h`; `showCases/carreauPoiseuille` | None (planned) | **Y** |
| **Straight-wall BCs** | Regularized (local), Skordos FD ("Interp"), Zou-He, Inamuro analytical, Equilibrium, generalized LS — factory `create*BoundaryCondition2D/3D()` with full corner/edge dynamics | Zou-He + zero-grad + convective | **Y** (regularized BC, 3D edge/corner completeness) |
| **Free-slip / symmetry** | SpecularReflection dynamics (`core/dynamics.h`); FreeSlipProfile3D for curved walls | None | **Y** (S) |
| **Outflow** | VirtualOutletDynamics + population-copy (`neumannCondition*.h`); NLD outflow | Zero-grad + convective | N |
| **Sponge / absorbing** | Viscosity/Smagorinsky SpongeZone2D/3D, masked, tanh/cos (`spongeZones*.h`); WaveAbsorptionDynamics | None | **Y** |
| **Curved boundaries (off-lattice)** | Guo extrapolation + Guo-FD, Bouzidi, Filippova-Hänel, Mei-Luo-Shyy, generalized (`offLattice/*OffLatticeModel3D.h`), moving-wall Ladd correction | None (planned) | **Y** (P1) |
| **Immersed boundary** | Multi-direct-forcing Inamuro, Favier/Inamuro delta kernels, force/torque reductions, 2D+ADE (`offLattice/immersedWalls*.h`); rotating impeller in `showCases/multiComponent3d/multiPhaseMixer3D.cpp` | None | **Y** (P2; impeller = M-F reactor relevant) |
| **STL → sim pipeline** | STLreader (`stlFileIO.h`), TriangleSet ops (`triangleSet.h`), DEFscaledMesh + TriangleBoundary3D, voxelizer 5-state + `revoxelize()`, auto inlet/outlet tagging, OffLatticeBoundaryCondition3D orchestrator; `showCases/aneurysm` | Primitives (circle/rect/sphere) | **Y** (P1 for CAD) |
| **Interior obstacles / MEM** | MomentumExchangeBounceBack (`bounceBackModels.h`); defineDynamics over voxel masks | Half-way BB rim + MEM force | N |
| **Porous media / Brinkman** | PartialBBdynamics (`partialBBdynamics.h`), `showCases/partialBounceBack`, GPU `multiComponentPorous`, `sandstone` | None | **Y** |
| **Shan-Chen** | Single + multi-component incl. external-field (`multiPhysics/shanChenProcessor*.h`); 4 psi functions (`interparticlePotential.h`) | SCMP + 2-component (2D) | Partial (3D SC, psi menu) |
| **Free surface (VOF)** | Full mass-tracking model w/ interface flags, contact angle, MPI bubble tracking (`freeSurfaceModel3D.*` 130KB, `bubbleHistory3D.h`); `damBreak3d`, `collidingBubbles3d` | None | **Y** |
| **He-Lee high-density-ratio** | heLeeProcessor3D.h | None; conservative Allen-Cahn planned (more modern) | Y (covered by planned ACE) |
| **Thermal/scalar (ADE)** | Full ADE stack: ~15 dynamics (BGK/TRT/RLB/Perko ± source; `advectionDiffusionDynamics.h`), ADE BCs, Boussinesq coupling incl. Smagorinsky (`boussinesqThermalProcessor*.h`); `boussinesqThermal2d/3d`, `settlingDrivenConvection` | None | **Y** (P1 — reactor mixing) |
| **Particles / tracers** | Particle3D hierarchy (Point/Rest/Verlet/NormedVelocity), pathline particles, Multi/Dense/Light fields, inject/count/advance (`src/particles/*.h`); `particlesInCone` | None | **Y** |
| **Deformable bodies (npFEM)** | ShapeOp RBC/platelet FEM (`coupledSimulators/npFEM/`; `bloodFlowDefoBodies`) | None | Y (P3, niche) |
| **Moving rigid bodies** | Prescribed motion via SurfaceVelocity + Guo/Ladd or IBM (`showCases/movingWall`, `multiPhaseMixer3D`); no 6-DOF FSI solver in core | None | **Y** (prescribed motion P1-adjacent for impellers) |
| **Grid refinement** | Two stacks: legacy MultiGridLattice2D/3D + ConvectiveRescaleEngine (`src/multiGrid/`), newer MultiLevelWrapper3D (`src/gridRefinement/`); `gridRefinement2d/dipole`, `gridRefinement3d/offLatticeExternalFlow` | Uniform only | **Y** |
| **Sparse domains** | SparseBlockStructure2D/3D — blocks only where fluid exists (`multiBlock/sparseBlockStructure*.h`); ~90% memory savings in aneurysm | Dense bounding box | **Y** |
| **Data processors (coupling)** | BoxProcessing/Dot/Bounded/reductive functionals, processor levels, modif-tracking, auto envelope sync (`atomicBlock/dataProcessingFunctional*.h`) | Fixed pass structure | Architectural difference (see B9) |
| **MPI** | MpiManager, ThreadAttribution block→rank, ParallelBlockCommunicator, static repartition only, **no dynamic load balancing** (`src/parallelism/`) | Halo exchange, bit-identical partition invariance | N (LBMFlow stronger on determinism; Palabos stronger on sparse decomposition) |
| **Threading** | None in core (MPI-only) | rayon + SIMD, 1,480 MLUPS 2D | N (LBMFlow superior) |
| **GPU** | acceleratedLattice: C++17 stdpar **via nvc++ only**, enum-dispatched kernels, hybrid CPU/GPU (`src/acceleratedLattice/*`); coProcessors stub D3Q19+BGK only | wgpu, 6–12 GLUPS, equivalence gates, vendor-portable | N (LBMFlow superior; note their GPU path runs cumulant) |
| **Checkpoint/restart** | saveBinaryBlock/loadBinaryBlock + dynamics-id registry (`io/multiBlockWriter*.*`, `mpiParallelIO.h` collective `writeRawData`; `core/dynamicsIdentifiers.h`) | None | **Y** (P1) |
| **Output** | VTK ImageData (base64), parallel VTK, structured grid, sparse VTK, **XDMF+HDF5** (`io/vtkDataOutput.h`, `xdmfDataOutput.h`, `hdfWrapper.h`), PPM/GIF via ImageMagick shell-out | PNG/CSV/VTK, manifest, HTML gallery | Partial (parallel/HDF5) |
| **Statistics/reductions** | BlockStatistics subscribe/gather auto every `collideAndStream()`; reductive functionals; transient windowed stats (`core/blockStatistics.h`, `transientStatistics3D.h`) | Probes/force series | Partial (free running + transient averaging) |
| **Units (SI)** | `IncomprFlowParam<T>`, `ComprFlowParam<T>`, general `Units3D` (`src/core/units.h`; exact API in B7) | Raw lattice units | **Y** (P1) |
| **Config** | XML via TinyXML (`TINYXML_xmlIO.h`, `plbInit`) | Self-describing JSON schema | N (LBMFlow superior) |
| **Inlet turbulence gen** | None (only initial random perturbations) | None | N (both lack — opportunity to exceed) |
| **Adjoint/autodiff** | None (`dco/` is a license file) | None | N |
| **Agent/async API, web GUI** | None | MCP server, WASM GUI | N (LBMFlow unique) |

## B. Top implementation techniques worth adopting

1. **Off-lattice boundary orchestration: geometry pipeline as staged data, not code.** STL →
   `TriangleSet` → `DEFscaledMesh` → `TriangleBoundary3D` → `voxelize()` (5-state flags,
   `offLattice/voxelizer.h`) → `OffLatticeBoundaryCondition3D` binds an interchangeable
   `OffLatticeModel3D` (Guo/Bouzidi/FH/MLS) to per-node precomputed dry-node lists with wall
   distances + `BoundaryProfile3D` per patch (auto-tagged inlets/outlets). *Rust mapping:*
   voxelization is preprocessing producing (a) a cell-flag field and (b) a compact list of
   `(cell, q, wall_distance, profile_id)` sorted by cell; the curved-BC pass is one more fixed
   pass after streaming iterating that list — fits `collide → halo → stream → open BCs → moments`
   without adopting data processors. Backend trait gets `apply_offlattice_bc(&records)`; GPU = a
   scatter kernel over the list.
2. **Comprehensive moment-space collision templates (RM/HM/CM/CHM/K/GH/RR).** One template family
   behind a uniform moments→relax→reconstruct interface, per lattice; GPU reuses via a
   `CollisionModel` enum. *Rust mapping:* a `CollisionKernel` trait on zero-sized types
   (Bgk/Trt/Mrt/Cm/Cumulant…) monomorphized into `step_band` + a runtime enum at the boundary;
   moment transforms are pure per-cell math, SIMD/WGSL-portable. Gets D3Q27 cumulant + recursive-
   regularized nearly for free once the transform machinery exists.
3. **Sparse block structure for irregular domains** (`SparseBlockStructure2D/3D`) — ~90% memory
   savings on vessel/reactor geometry. *Rust mapping:* domain = set of dense blocks + adjacency;
   the MPI halo machinery already handles inter-block exchange; T13 extends: sparse vs dense must
   agree bitwise on the covered region.
4. **Dynamics-id registry + decompose/recompose checkpointing** with MPI-collective offset write
   into one file. *Rust mapping:* far simpler for us (no per-cell virtual dynamics) — checkpoint =
   (scenario hash, step, f-deviation arrays per block, RNG/probe state) + versioned header; ranks
   write at precomputed offsets. Deviation storage (f−w) keeps f32 restarts validation-grade.
5. **Sponge zones as an omega/Smagorinsky ramp field** (`spongeZones3D.h`) — tanh/cos spatial ramp
   of relaxation near outflow kills reflections; a per-cell omega multiplier consumed in collide
   (same storage as our per-cell body force). No new pass.
6. **TRT runtime magic parameter** (`setMagicParam`) — Λ=3/16 suits half-way BB but Λ=1/4 suits
   other BC mixes. *Rust mapping:* optional `magic` in collision config, default 3/16. Effort: S.
7. **IncomprFlowParam-style unit converter (verified firsthand in `src/core/units.h`).** Inputs:
   `physicalU, latticeU, Re, physicalLength, resolution, lx,ly,lz`. Derived: `dx=L_phys/N`;
   `dt=dx·u_lat/u_phys`; `nu_lat=u_lat·N/Re`; `tau=3·nu_lat+0.5`; helpers `nCell(l)`, `nStep(t)`.
   Also `ComprFlowParam` (rho/T/Pe) and general `Units3D`. *Rust mapping:* a `FlowParams` struct in
   lbm-scenario resolved at validation into raw lattice units (core stays lattice-unit-only —
   invariant preserved); emit a derived-quantities block into manifest.json. **Corroborates OpenLB
   B2 and OpenFOAM #1 — the user fixes (u_phys, L_phys, Re) + two numerical knobs (u_lat=Mach, N).**
8. **Consistent-Smagorinsky formulation catalog** — distinguishes naive omega-rescale vs consistent
   (Malaspinas, stress-mode-isolated) LES per collision family. When we build LES, the consistent
   variant on TRT/central-moments is the correct target; slots into the moment-space kernel (B2).
9. **Data-processor pattern — adopt the concept selectively, not the machinery.** Their extensibility
   (and slowness/unfusability) comes from user-registered `BoxProcessingFunctional`s at integer
   levels after every `collideAndStream()`. LBMFlow's fixed pass structure is the right default
   (enables step_band fusion + bit gates), but multiphysics coupling (Boussinesq ADE↔NS, Shan-Chen
   force) needs *one* sanctioned extension point: a typed `CouplingPass` slot between collide and
   stream, registered per scenario, with declared read/write field sets so halo exchange can be
   scheduled — a static, enumerated version of their processor levels. Do NOT adopt arbitrary
   functionals, runtime level graphs, or envelope options.

## C. Prioritized gap list

**P1 — must-have to claim LBM-capability superiority:**
1. D3Q27 + central-moment/cumulant collision (planned; parity target = their K/CM/CHM/RR) — **M**
2. MRT + regularized BGK variants (incl. regularized velocity/pressure BC) — **M**
3. Curved-boundary BC: Bouzidi first (post-stream record pass), Guo extrapolation later; moving-wall
   Ladd correction for impellers — **L**
4. STL import + voxelization (triangle set, in/out flags, wall distances, per-patch profiles) — **L**
   (prerequisite of 3)
5. Advection-diffusion lattice (D2Q5/D3Q7) + Boussinesq + scalar mixing — **L** (M-F reactor: mixing
   time, temperature)
6. Smagorinsky LES (consistent, per-cell omega) — **M**
7. SI unit converter (IncomprFlowParam clone at scenario boundary) — **S**
8. Checkpoint/restart (versioned binary, MPI offset write) — **M**
9. Free-slip/symmetry BC (specular reflection) — **S**
10. Sponge zones (viscosity ramp near outlets) — **S**
11. Prescribed-motion moving boundaries (rotating impeller; Palabos `multiPhaseMixer3D`) — **M** given 3

**P2 — should-have:** free-surface VOF (**L**); sparse block domains (**M-L**); immersed boundary
(Inamuro multi-direct forcing) (**M-L**); non-Newtonian Carreau/power-law (**S-M**, small once
per-cell omega exists); grid refinement (2:1 convective-scaled) (**L**); tracer/pathline particles +
residence-time stats (reactor RTD) (**M**); parallel HDF5/XDMF + parallel VTK (**M**); Shan-Chen 3D +
psi menu (**M**); running/transient statistics (**S-M**); TRT magic exposure (**S**); partial
bounce-back porous media (**M**).

**P3 — optional:** entropic/KBC (**M**); He-Lee (skip unless ACE slips); Inamuro/Skordos/equilibrium
straight-wall BC menagerie (**S** each, low value); bubble tracking analytics (**M**); extended
lattices D3Q39/D3Q121 (**L**, no pull); npFEM deformable bodies (**XL**, out of scope); acoustic wave
absorption (**M**).

## D. Explicitly NOT worth adopting

1. **Per-cell virtual `Dynamics*` AoS storage** (`blockLattice2D.hh`: `Cell[nx*ny]` each with a
   dynamics pointer, virtual collide per cell) — kills SIMD/GPU, is the root cause Palabos needs a
   *separate* `acceleratedLattice` for GPU. Our q-major SoA + static dispatch is strictly better;
   keep BCs as index lists, not per-cell objects.
2. **The general data-processor runtime** — adopt only the single typed coupling slot (B9).
3. **stdpar/nvc++ GPU strategy** — vendor-locked; our wgpu backend exceeds it in portability with
   equivalence gates. Steal only the enum-kernel dispatch idea (already how we work).
4. **ImageMagick shell-out for GIF/PPM** — fragile external dep; our PNG writer is superior.
5. **TinyXML config + argv globals** — our JSON schema + MCP is a generation ahead.
6. **Dual grid-refinement stacks** (legacy `multiGrid/` + newer `gridRefinement/`, both alive) — a
   warning: design ONE refinement architecture and never ship two.
7. **coProcessors** — abandoned D3Q19+BGK offload stub; no content.
8. **Inamuro analytical BC** (D2Q9-only) and **extended lattices D2Q37/D3Q121** — academic surface
   with no commercial demand; note as deliberate exclusions in the comparison table.
9. **ThreadAttribution's unused thread dimension** — rayon covers it.
10. **AGPL code reuse of any kind** — re-derive from published methods; implement from papers.

---

## Consolidated cross-cutting priorities (all three codes agree)

1. **SI UnitConverter** — named by OpenLB (B2) and OpenFOAM (#1) and the user. Single highest
   leverage; low solver risk. Do first.
2. **Curved boundary (Bouzidi/off-lattice)** — OpenLB flagship gap; Palabos ships a whole
   off-lattice model family (Guo/Bouzidi/FH/MLS). Biggest physics-accuracy upgrade. Pair with an
   STL→voxelization pipeline (both codes treat geometry as staged precomputed data — a cell-flag
   field + a sorted `(cell, q, wall_distance, profile)` record list that fits our fixed passes).
3. **Dynamics composition refactor** — the multiplier that makes MRT/cumulant/LES/non-Newtonian
   cheap (OpenLB B1; Palabos's composite-dynamics is the same idea).
4. **Checkpoint/restart** — table stakes, named by both OpenLB (B8) and OpenFOAM. Low effort.
5. **Composable observer/function-object framework + SI-aware output** — OpenFOAM's product moat;
   OpenLB's functor framework is the same concept.
6. **LES (Smagorinsky/WALE) as effective-omega** — cheap once composition lands; unblocks M-F.
