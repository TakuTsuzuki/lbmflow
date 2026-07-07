# Bioprocess Pivot (2026-07-07)

Lifecycle: living (pivot announcement, kept in sync with capability registry
and validation status).

## 1. What changed

LBMFlow's product mission is now a **bioprocess-specific CFD core**, not a
general-purpose LBM simulator. Development priorities are re-anchored on the
QOI (quantities of interest) that drive cell-culture and stirred-reactor
process decisions:

- Power number **Np**, power-per-volume **P/V**
- Mixing time (scalar CV → **t95 / t99**)
- **Gas holdup** (resolved and hybrid)
- Sauter mean diameter **d32**
- Volumetric mass-transfer coefficient **kLa** (resolved-interface and PBM
  variants)
- **Shear-exposure** distributions and cell-damage risk
- Oxygen-exposure history and OUR
- Cell / microcarrier damage risk and suspension metrics
- **Scale-up operating window** (P/V, tip speed, kLa, mixing time,
  P95-shear-limited feasibility set)

Everything else — GPU throughput, WASM GUI, generic CFD compatibility,
extended lattice choices, FP16 grid capacity — is a *means*, not a product
claim, and is deferred until the QOI pipeline validates end-to-end.

## 2. What is retracted as a product claim

The following prior claims are demoted to "engineering / demo" status until a
bioprocess-tier validation lands (see
[SPEC_BIOPROCESS_CORE.md](SPEC_BIOPROCESS_CORE.md) tiers):

| Prior product claim | New status | Reason |
|---|---|---|
| Commercial-grade general-purpose LBM simulator | Retracted | Product is now bioprocess-specific; general LBM benchmarks are demos, not validation. |
| Shan-Chen as production gas-liquid model | Unsupported | Insufficient density ratio; spurious currents; production path is conservative Allen-Cahn phase field (BCFD-040..048). |
| FP16 storage as validation-grade | Unsupported | Capacity/throughput mode only. Never used for a QOI reported to a bioprocess decision. |
| 3D GPU as coupled-bioprocess-physics path | Unsupported | Scenario GPU dispatch rejects multiphase, rotor, particles, non-rest init, and force probes — the *entire* bioprocess coupled physics. Product CPU path first. |
| GMP / CMC evidence claims | Unsupported | No bioprocess QOI has calibration + holdout + UQ + sensitivity records today. Evidence gate (BCFD-091) not implemented. |
| M-A … M-F milestones · T1 … T18 acceptance criteria | Superseded | Replaced by BCFD-000..110 tickets and VB-01..VB-08 validation groups. See `docs/archive/2026-07-07-pivot/`. |

## 3. What is preserved

- **LBM core code** (`crates/lbm-core`): D2Q9/D3Q19/D3Q27 lattices, CPU
  scalar/SIMD backends, wgpu GPU backend, MPI halo exchange, WALE LES,
  rotating IBM, Bouzidi walls, Guo forcing, Shan-Chen SCMP/MCMP. The
  underlying numerics and their test coverage remain. What changes is *what
  gets exposed as a product QOI* and *what gets validated to a bioprocess
  band*.
- **Legacy scenarios and CLI presets** run and their tests remain green,
  but the CLI emits a warning:
  `legacy LBM demo preset; not bioprocess decision-grade`
- **Physical rigor prime directive**: unchanged. Every physics term is
  either resolved or a literature-backed closure with derivation, validity
  domain, and its own validation test. See
  [.claude/skills/lbmflow-physics-discipline](../.claude/skills/lbmflow-physics-discipline).
- **English-only artifacts** and the codex fan-out development pattern
  remain.

## 4. Route to product

Four milestones, sequential. Each milestone's exit criteria are the QOIs it
adds becoming available and machine-readable with provenance:

- **M0** — single-phase stirred tank, Np/P/V, shear percentiles, mixing time,
  report scaffold.
- **M1** — resolved gas-liquid (Allen-Cahn), sparger, gas holdup, oxygen
  scalar, synthetic-kLa fit.
- **M2** — cell/microcarrier exposure, UQ/sweep, scale-up window, evidence
  gate, CLI/MCP surface.
- **M3** — point bubbles + PBM, kLa from interfacial area, hybrid gas
  bookkeeping.

MPI parallel work (BCFD-100..102), GPU bioprocess support, GUI, generic CAD
meshing beyond BCFD-023, and evidence-grade claims are behind hard cut lines
defined in [PLAN.md](PLAN.md).

## 5. Discipline unchanged

- Every new physical model must include validation tests before being marked
  Engineering or Evidence.
- Every unsupported combination must fail loudly with a structured error, not
  silently fall back.
- Every QOI must include units, method, time window, averaging region, and
  validation tier.
- Every experimental result gets a behavior-validity review before being
  reported — metric passing a band ≠ pattern is physically valid.
- Evidence-tier claims require calibration/holdout separation, UQ, and
  mesh/time-step sensitivity records.

## 6. Related documents

- [SPEC_BIOPROCESS_CORE.md](SPEC_BIOPROCESS_CORE.md) — intended use, tiers,
  QOI catalog.
- [VALIDATION_BIOPROCESS.md](VALIDATION_BIOPROCESS.md) — VB-01..VB-08 groups.
- [CREDIBILITY_BIOPROCESS.md](CREDIBILITY_BIOPROCESS.md) — calibration /
  holdout / UQ policy.
- [MODEL_RISK_MATRIX.md](MODEL_RISK_MATRIX.md) — per-model risk table.
- [PLAN.md](PLAN.md) — BCFD-000..110 tickets and M0–M3 milestones.
- [LIMITATIONS.md](LIMITATIONS.md) — machine-readable capability status.
- [archive/2026-07-07-pivot/](archive/2026-07-07-pivot/) — pre-pivot plan,
  T1–T18 validation matrix, R-Phase / M-A..M-F work, V&V ledger, whitepaper,
  claims-ledger.
