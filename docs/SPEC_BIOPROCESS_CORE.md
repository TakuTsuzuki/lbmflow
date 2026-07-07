# Bioprocess CFD Core — Specification

Lifecycle: living (owning doc for intended use, forbidden use, tiers, QOI
catalog).

## 1. Intended use

LBMFlow is a **bioprocess-specific CFD core**. Its intended use is the CFD
side of stirred-tank cell-culture / bioreactor process design:

- Single-phase stirred-tank hydrodynamics (Np, P/V, discharge Nq, mixing
  time)
- Shear-rate and viscous-stress fields, wall-shear diagnostics, and
  cell / microcarrier shear-exposure histories
- Passive-scalar transport and mixing-time QOI
- Resolved-interface gas-liquid (conservative Allen-Cahn) for sparger-driven
  aeration and free-surface handling
- Oxygen mass transport, Henry-equilibrium interfacial flux, kLa QOI
- Point-bubble / PBM engineering-mode gas-liquid with d32 and interfacial
  area
- Cell / microcarrier trajectories, exposure integrals, damage risk models
- Scale-up operating-window evaluation across sweep results

## 2. Forbidden use

The tool is not to be used to justify:

- Any **GMP / CMC filing claim** unless the associated QOIs cleared the
  evidence-tier gate (see [CREDIBILITY_BIOPROCESS.md](CREDIBILITY_BIOPROCESS.md)
  and BCFD-091).
- **Validated kLa** without a documented calibration + holdout dataset pair.
- Any **production gas-liquid decision** using the Shan-Chen path. The
  production gas-liquid model is conservative Allen-Cahn phase field
  (BCFD-040..048), and even that requires the phase-field validation suite
  (BCFD-048).
- Any decision derived from a QOI whose report is missing units, method,
  averaging window, or averaging region.
- Any decision derived from a `max` shear or `max` exposure alone.
  Percentile distributions (P50 / P90 / P95 / P99) and fraction-above-
  threshold are required.

## 3. Credibility tiers

Each capability, and each QOI derived from it, carries a tier:

- **Tier 0 – Screening.** For qualitative comparison across configurations.
  Runs may be under-resolved. Never used as an evidence artifact.
- **Tier 1 – Engineering.** Runs are mesh-and-time-step resolved to a
  documented convergence band; the physics model has its validation suite
  green; QOI is emitted with full provenance; but no calibration/holdout has
  been pinned. Suitable for design-of-experiments and internal decisions.
- **Tier 2 – Evidence.** Requires the Engineering tier plus a calibration
  dataset AND an independent holdout dataset (BCFD-082), an uncertainty
  interval on the reported QOI (BCFD-083), and mesh / time-step sensitivity
  records; enforced by the evidence gate (BCFD-091). Suitable for external
  submissions where the tool's role must be defensible.

Current status: **all bioprocess QOIs are pre-Tier-1** until their BCFD
ticket lands with its validation entry.

## 4. QOI catalog

Every QOI ships with `units`, `method`, `time_window`, `averaging_region`,
`source_fields`, and `validation_tier` metadata (BCFD-080). If the QOI
cannot be computed (e.g. required source field disabled), it is emitted as
`skipped` with a reason string, not silently omitted.

| QOI | Units | Method (default) | Depends on | Sensitivity guardrail |
|---|---|---|---|---|
| Np | dimensionless | Guo-force resolved torque `Tq → P = ω Tq → P/(ρ N³ D⁵)` | rotating IBM + torque probe | mesh convergence to ≤5% |
| P/V | W/m³ | `P / V_working` | Np | working volume declared explicitly |
| Nq | dimensionless | discharge-surface volumetric flow `Q / (N D³)` | discharge surface definition | skipped with reason if surface undefined |
| t95, t99 | s | scalar CV(t) fit: `CV ≤ 0.05·CV₀` / `≤ 0.01·CV₀` | passive scalar ADE + pulse injection | uniform initial scalar → skipped |
| gas holdup ε_g | dimensionless | resolved: `⟨1−φ⟩`; hybrid: resolved + Σ bubbles / V | phase field OR point bubbles | threshold, averaging region required |
| d32 | m | `Σ nᵢ dᵢ³ / Σ nᵢ dᵢ²` over PBM bins | PBM | bin conservation checked |
| kLa | 1/s (also 1/hr) | resolved: `∫ kL a (C* − C) dV`; PBM: `kL · a_PBM`; dynamic: fit `dC/dt = kLa (C*−C)` | oxygen scalar + interfacial model | fit R², window, method recorded |
| shear exposure | Pa·s | `∫ max(0, τ−τ_c)^m dt` over tracer | shear-rate field + cell tracers | percentile distribution + max, not max alone |
| oxygen exposure | history | tracer C_L(t) resampled onto trajectory | oxygen scalar + cell tracers | percentile distribution |
| microcarrier suspension | dimensionless | settled fraction, height distribution, residence-near-impeller | microcarrier particles | Schiller-Naumann validity range enforced |
| scale-up window | operating-condition set | feasibility of constraint set across sweep | Np, kLa, shear, mixing, gas-holdup | infeasibility → explicit conflict table |

## 5. Scenario schema (BCFD-003)

- Root type `BioprocessScenario`, distinct from the legacy `Scenario`.
- Required sections: `version`, `name`, `credibility_tier`, `reactor`,
  `fluids`, `operation`, `physics`, `qoi`, `run`, `outputs`. Optional:
  `cells`.
- `reactor.kind = StirredTank` in M0..M2; other reactor kinds are
  Unsupported.
- `physics` is a discriminated union of `single_phase`,
  `resolved_phase_field`, `point_bubble`, `hybrid`, `passive_scalar`,
  `oxygen`, `cell_tracer` — combinations are validated against the capability
  registry.
- Unknown fields are rejected.
- Impossible combinations are rejected (e.g. `kla` without oxygen; sparger
  without a gas model; evidence tier without a validation dataset registry
  reference).

Old `Scenario` inputs continue to parse for backward compatibility; the
legacy path runs with the demo warning.

## 6. Unit and dimensionless-feasibility layer (BCFD-004)

Every bioprocess scenario emits a `UnitReport` containing:

- SI inputs, working volume, geometry lengths
- Dimensionless groups: Re, Fr, We, Eo, Mo, Sc, Pe, St, Ma_lattice, Cn, Pe_φ
- Feasibility diagnostics:
  - `Ma_lattice > 0.1` → warning; `> configured max` → hard reject
  - `τ` too close to 0.5 → reject
  - interface width too thin → reject (Cn feasibility)
  - resolved bubble under-resolved → reject
  - scalar diffusion unstable → reject
- Matching priority (declared explicitly): Re → density/viscosity ratio +
  We/Eo → Fr → Sc/Pe/Da → St. Partial matching emits a warning.

## 7. Non-negotiables

The following are enforced mechanically, not by convention:

- **No silent fallback.** Unsupported combinations fail with a structured
  error carrying an `UnsupportedReason` and remediation.
- **QOI provenance.** Every QOI serialises with its metadata bundle; missing
  metadata fails serialisation.
- **Validation before Engineering-tier.** A physics model without its BCFD
  validation entry green is Experimental / Unsupported, regardless of what
  the code can technically run.
- **Calibration ≠ validation.** Datasets used for calibration cannot be
  reused as holdouts (BCFD-082).
- **`max` is never the only report.** Distributions and percentiles are
  mandatory for shear, exposure, and stress-derived QOIs.
- **Evidence gate is mechanical.** BCFD-091 must return `EvidenceReady`
  before any report is emitted with the evidence label.
