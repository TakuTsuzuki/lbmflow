# Model Risk Matrix — Bioprocess CFD

Lifecycle: living (owning doc for per-model risk assessment and the
alternative-closure enumeration required for Evidence-tier model-form
uncertainty per [CREDIBILITY_BIOPROCESS.md](CREDIBILITY_BIOPROCESS.md) §3).

Every model row lists: **purpose**, **assumption(s) that could bite**,
**validity domain**, **failure mode**, **alternative closures** (used for
model-form UQ), and **current tier**. When a QOI depends on a model, the
QOI's evidence claim inherits the model's tier ceiling.

## 1. Fluid solver core

| Model | Assumption | Validity | Failure mode | Alternatives | Tier |
|---|---|---|---|---|---|
| LBM D3Q19 BGK | Isotropic collision, unity Pr | Low-to-moderate Re, weakly compressible | Anisotropy at high Re, Galilean-invariance drift under non-zero mean advection | TRT (Λ=3/16), D3Q27, cumulant (with Galilean holdout open) | Engineering |
| TRT (Λ=3/16) | Exact half-way wall for Poiseuille | Simple geometries; complex BCs need care | Loss of exactness under moving/non-orthogonal walls | BGK, cumulant, D3Q27 | Engineering |
| Cumulant D3Q19 | Central-moment collision improves high-Re stability + isotropy | High Re turbulent flows | Open Galilean-invariance holdout at finite frame velocity (PHYSICS.md 2026-07-07) — use TRT/f64 or D3Q27 references for validation-grade frame-shift studies | TRT/f64, D3Q27 | Engineering (with hold-out flagged) |
| Guo forcing | 2nd-order body force; F/2 in u | Any external volume force, up to |F| ≪ 1 lu units | Force can invalidate low-Ma assumption if large | Explicit Kupershtokh / He-Chen-Doolen (deprecated in this repo) | Engineering |

## 2. Turbulence

| Model | Assumption | Validity | Failure mode | Alternatives | Tier |
|---|---|---|---|---|---|
| WALE LES | Wall-resolved eddy viscosity, near-wall damping by construction | Structured grids; Re_τ not extremely high | Near-wall over/under-prediction outside characterisation range; requires y+ diagnostics (BCFD-033) | Smagorinsky (rejected as default), Dynamic Smag (not implemented), RANS (out of scope) | Engineering (channel Re_τ=180) |
| ν_t clipping | Diagnosed numerical-stability guard; default off | Only when instability observed | Silent turbulence-model calibration when misused | Off (raw WALE); report `clipped_fraction` if on | Engineering |

## 3. Multiphase / gas-liquid

| Model | Assumption | Validity | Failure mode | Alternatives | Tier |
|---|---|---|---|---|---|
| Conservative Allen-Cahn phase-field (BCFD-040..048) | Diffuse interface with mobility M; interface width W resolved | Moderate density ratio (guard until high-ratio path validated); We/Eo consistent with interface Cn | Mass drift, interface too thin (Cn feasibility) | Cahn-Hilliard (deferred), colour-gradient LBM (deferred), VOF-PLIC (out of scope) | Experimental → Engineering after BCFD-048 |
| Shan-Chen SCMP/MCMP | Pseudo-potential immiscibility | Low density ratio (~1:10); coupled σ-ratio; documented spurious currents | Density-ratio ceiling; spurious velocities at interface | Allen-Cahn (production path), colour-gradient | Demo only — NOT product |
| Sparger resolved injection (BCFD-046) | Gas-only injection at orifice; φ=0 boundary | Orifice resolved by ≥4 lattice cells | Under-resolved orifice → non-physical gas volume | Point-bubble mode (BCFD-070+) | Experimental until VB-05 green |
| Point-bubble (BCFD-070..074) | Sub-grid bubble entity with drag + buoyancy + added mass | d_bubble ≪ grid Δx; ρ_b/ρ_l ≪ 1; α_g moderate | Breaks at high α_g without continuum treatment; high-Re_p closures out of range | Resolved phase field (BCFD-040..048), hybrid (BCFD-075), Euler-Euler (out of scope) | Experimental until VB group green |
| PBM (BCFD-073) | Bin-conservative breakup / coalescence kernels | Kernel-domain-specific; validity declared per kernel | Coalescence dominance; wrong d32 tail | Luo-Svendsen, Prince-Blanch, constant | Experimental |

## 4. Mass transfer / oxygen / kLa

| Model | Assumption | Validity | Failure mode | Alternatives | Tier |
|---|---|---|---|---|---|
| Oxygen scalar ADE (BCFD-050) | Passive scalar transport; Sc set explicitly | Sc within resolved diffusion; boundary sources known | Sc too high → grid-Sc feasibility fail | Reject at BCFD-004 | Experimental |
| Henry equilibrium interfacial flux (BCFD-051) | Sharp local equilibrium at interface; kL constant or correlated | Constant kL for smoke; correlations validity-scoped | Wrong kL model dominates result | kL correlation options; calibrated kL | Experimental (constant) → Engineering (calibrated) |
| kLa from PBM (BCFD-074) | `a = 6 α_g / d32`; kL from selected model | d32 within PBM validity; kL model applicable | Wrong d32 → wrong a → wrong kLa | Resolved-interface kLa (BCFD-052) | Experimental |
| Dynamic gassing fit (BCFD-052) | First-order response `dC/dt = kLa (C*-C)` | Steady interfacial area during fit; homogeneous mixing assumption | Under-mixed tank → wrong kLa fit | Compartment-wise kLa (future) | Engineering when VB-06 green |

## 5. Particles / cells

| Model | Assumption | Validity | Failure mode | Alternatives | Tier |
|---|---|---|---|---|---|
| Cell tracers (BCFD-060) | Massless; sample flow only | ρ_cell ≈ ρ_liquid; deformation irrelevant to trajectory | Misses lift / near-wall effects for stiff cells | Microcarrier mode (BCFD-062) | Engineering when VB-07 green |
| Microcarrier one-way (BCFD-062) | Schiller-Naumann drag; d, ρ, buoyancy, restitution known | Re_p ≤ 800 (enforced) | Beyond validity → structured error | Two-way (BCFD-063), four-way (not implemented) | Engineering when validated |
| Shear damage integral (BCFD-061) | `E = ∫ max(0, τ-τ_c)^m dt`; m constant | m and τ_c from cell-line calibration | Wrong m / τ_c → wrong ranking | γ̇ threshold, ε threshold placeholder | Experimental until cell-line calibrated |

## 6. Boundaries

| Model | Assumption | Validity | Failure mode | Alternatives | Tier |
|---|---|---|---|---|---|
| Half-way bounce-back | Wall at midpoint between rim and fluid | Static and moving walls | 1-cell rim mandatory; edge cases if rim broken | Bouzidi (curved); Zou-He (velocity/pressure BC only) | Engineering |
| Bouzidi (curved) | 2nd-order interpolated bounce-back | Analytic geometries (circle, sphere) currently | Non-analytic geometry → requires STL import path (BCFD-023) | Half-way BB | Engineering |
| Zou-He inlet | Velocity or pressure prescribed at face | Uniform or profiled velocity; not corner | Corner cells require rim; not for curved walls | Guo-based inlet (not exposed) | Engineering |
| Rotating IBM (BCFD-021) | Prescribed rigid rotation; U = ω × r | Impeller as marker set; no elastic deformation | Not general FSI; no structural DOF | MRF (out of scope), overset (out of scope) | Engineering |
| Free surface (BCFD-045) | Top face reflects; degassing lets gas out | Small surface deformation; no wave physics | Cannot model spillover / vortex depression | Full VOF (out of scope) | Experimental |
| Contact angle (BCFD-044) | Wall imposes local ∇φ at contact line | Static contact angle only | Dynamic contact angle absent | Cox-Voinov (not implemented) | Experimental |

## 7. Numerical precision

| Precision | Assumption | Validity | Failure mode | Alternatives | Tier |
|---|---|---|---|---|---|
| f64 | Full double precision | Validation-grade | Slower / more memory | f32 (with deviation storage), f16 (capacity only) | Engineering / Evidence |
| f32 with deviation storage | Store `f-w`; recover full precision on bounce-back (linear invariant) | Most 2D/3D CPU/GPU runs | Zou-He needs constant term folded in; long transients still accumulate | f64 | Engineering (bands frozen) |
| f16 storage / f32 compute | Capacity/throughput mode | Steady flows that re-converge | Long transients accumulate storage rounding | f32/f64 | NOT validation-grade |

## 8. Backends

| Backend | Assumption | Validity | Failure mode | Alternatives | Tier |
|---|---|---|---|---|---|
| CpuScalar | Reference | Small grids; correctness gate | Slow | CpuSimd | Engineering |
| CpuSimd | Fused collide+stream+moments | Equivalence to CpuScalar to ≤1e-5 rel (T13/T14) | If storage / step order changed, must re-verify | CpuScalar | Engineering |
| Wgpu GPU | f32 compute; feature `gpu` | 2D scenarios end-to-end; 3D scenarios reject multiphase / rotor / particles / non-rest init / force probes — the entire bioprocess coupled physics | Silent perf regression from loaded-window measurement (documented trap) | CPU | Unsupported for bioprocess coupled runs |
| MPI | Halo exchange; per-rank ownership | Small-node verified; 64-rank weak scaling RED (needs cluster) | Global-array replication in `MpiSolver::new` (BCFD-100 fixes this) | Single-node | Engineering (small-scale) |

## 9. QOI-level risk

| QOI | Model risk | Alternative | Evidence path |
|---|---|---|---|
| Np | LBM core + rotating IBM torque | Alternative impeller marker density | VB-01 + calibration/holdout per geometry |
| P/V | Follows Np | — | Inherits Np |
| Mixing time | Passive scalar + WALE turbulence | Scalar with LES off (laminar limit) | VB-02 + Nθ correlation |
| Gas holdup (resolved) | Phase field + sparger | Point-bubble mode | VB-05 + calibration |
| Gas holdup (hybrid) | Resolved + PBM double-counting risk | Resolved-only, PBM-only bounds | Requires BCFD-075 validation |
| d32 | PBM kernels | Alternative kernel set | Requires calibrated kernel |
| kLa (resolved) | Interfacial area from δ-approximation | PBM-based kLa | VB-06 + calibrated kL |
| kLa (PBM) | d32 and kL model | Alternative kernel + alternative kL | Same |
| Shear exposure | γ̇ field + integral | γ̇ threshold vs stress threshold vs ε threshold | VB-07 + cell-line calibration |
| Scale-up window | Composite | Constant-P/V vs constant-tip vs constant-kLa | VB-08 + holdout at target scale |
