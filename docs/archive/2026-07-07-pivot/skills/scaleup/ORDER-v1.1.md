# SCALEUP v1.1 — Implementation Order (received from the skills-design session, 2026-07-05)

PM sanity-check annotations (against main at receipt):
- Strain-rate tensor gathers LANDED tonight (gather_strain_rate/gather_shear_rate on
  native Solver, machine-precision Couette verification) — resolved viscous
  epsilon = 2 nu S:S is computable at the Rust API level; CLI/scenario output channel
  NOT yet wired (order C queued) → expect YELLOW at the user surface.
- Rotating-solid torque → Np: RED (momentum-exchange exists for STATIC solids only;
  the volume-penalization impeller is a body-force emulation with no torque).
- Lagrangian lifelines, scalar ADE theta95, LES eps_sgs: RED until M-F tracks land.
- eps-bar vs P/V definitions in §1.4 verified dimensionally consistent and exactly
  matching REQ_STIRRED_REACTOR.md rev.4 §2.1 (Np = P/(rho N^3 D^5), P = Omega T_q).
- SU-S0 must BUILD ON docs/skills/b1-capability-map.md (CLI/MCP/schema/artifact
  surface already audited tonight) and audit only the scaleup-specific deltas.

[Order body follows — v1.1 verbatim from the skills-design session; see the
cross-session message of 2026-07-05. Key sections: Part 1 survey (do not commit
verbatim until source-pinned), Part 2 S-Advisor/S-Fingerprint/S-Match specs with
verification gates, Part 3 DAG (SU-S0 + SU-advisor parallel now; fingerprint/match
gated), fixed artifact paths under docs/skills/scaleup/.]
