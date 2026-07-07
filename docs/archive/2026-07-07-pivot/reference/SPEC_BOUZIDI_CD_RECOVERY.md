# SPEC: Bouzidi T8 Cd recovery — LANDED

Status: **LANDED**. The Bouzidi interpolated bounce-back closes T8 (Schäfer–Turek
Cd band). See `crates/lbm-core/src/bouzidi.rs`, tests
`crates/lbm-core/tests/bouzidi.rs` + `accuracy_audit_bouzidi.rs` +
`validation_cylinder.rs` (T8).

Retained diagnosis-matrix headline (kept for the ANOM record):
a steady +5–8% Cd bias at fixed Re on the pass-1 half-way-BB baseline is the
classical O(1) effective-diameter error → check radius/qd convention, Re
definition (U_mean vs U_max), blockage, and momentum-exchange force convention
(Wen 2014) before adjusting the band. Never loosen the band; if the evidence
indicates the band itself encodes a different configuration, STOP and report
(band governance).

Convergence: D={10,20,40} at fixed Re=20 → order ≥ 1.7 asserted.
qd=1/2 degenerates bitwise to half-way BB (pinned by degeneracy test).
