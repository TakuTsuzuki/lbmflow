use lbm_core::geometry::{
    build_stirred_tank_geometry, GridSpec, SpargerTemplate, TankBottom, TankSpec,
    SPARGER_ORIFICE_MIN_CELLS,
};
use lbm_core::sparger::{apply_resolved_gas_injection, ResolvedGasInjectionSpec, SpargerGasLedger};

fn ring_geometry() -> lbm_core::geometry::StirredTankGeometry {
    let dx_m = 0.01;
    build_stirred_tank_geometry(
        GridSpec {
            dims: [100, 100, 100],
            dx_m,
        },
        TankSpec {
            vessel_diameter_m: 0.9,
            liquid_height_m: 0.9,
            bottom: TankBottom::Flat,
        },
        &[],
        &[],
        &[SpargerTemplate::Ring {
            center_z_m: 0.15,
            outer_radius_m: 0.25,
            orifice_count: 8,
            orifice_diameter_m: SPARGER_ORIFICE_MIN_CELLS * dx_m,
            gas_volumetric_flow_m3_per_s: Some(1.0e-7),
            inlet_phase_gas: true,
        }],
    )
    .unwrap()
}

#[test]
fn ring_sparger_gas_injection_decreases_phi_only_on_orifice_cells() {
    let geom = ring_geometry();
    let mut phi = vec![1.0f64; geom.solid.len()];
    let mut ledger = SpargerGasLedger::default();
    apply_resolved_gas_injection(
        &mut phi,
        &geom.sparger_mask,
        &geom.solid,
        ResolvedGasInjectionSpec {
            gas_volumetric_flow_m3_per_s: 1.0e-7,
            dt_s: 0.01,
            dx_m: geom.dx_m,
            orifice_diameter_m: SPARGER_ORIFICE_MIN_CELLS * geom.dx_m,
        },
        &mut ledger,
    )
    .unwrap();

    let mut changed_on_sparger = 0usize;
    let mut changed_off_sparger = 0usize;
    for ((&p, &is_sparger), &is_solid) in phi.iter().zip(&geom.sparger_mask).zip(&geom.solid) {
        if is_solid {
            continue;
        }
        if p < 1.0 && is_sparger {
            changed_on_sparger += 1;
        }
        if p < 1.0 && !is_sparger {
            changed_off_sparger += 1;
        }
    }
    assert!(
        changed_on_sparger > 0,
        "gas injection must decrease phi at ring-sparger orifice cells"
    );
    assert_eq!(
        changed_off_sparger, 0,
        "resolved sparger injection must not alter non-orifice cells directly"
    );
}

#[test]
fn gas_volume_ledger_matches_q_times_t_within_two_percent() {
    let geom = ring_geometry();
    let mut phi = vec![1.0f64; geom.solid.len()];
    let mut ledger = SpargerGasLedger::default();
    let q = 5.0e-8;
    let dt = 0.01;
    let steps = 20;
    for _ in 0..steps {
        apply_resolved_gas_injection(
            &mut phi,
            &geom.sparger_mask,
            &geom.solid,
            ResolvedGasInjectionSpec {
                gas_volumetric_flow_m3_per_s: q,
                dt_s: dt,
                dx_m: geom.dx_m,
                orifice_diameter_m: SPARGER_ORIFICE_MIN_CELLS * geom.dx_m,
            },
            &mut ledger,
        )
        .unwrap();
    }
    let expected = q * dt * steps as f64;
    let rel = (ledger.injected_gas_volume_m3 - expected).abs() / expected;
    assert!(
        rel <= 0.02,
        "ledger injected volume {} must match Q*t {expected} within 2%, rel={rel}",
        ledger.injected_gas_volume_m3
    );
    assert!(
        ledger
            .pressure_diagnostics
            .iter()
            .all(|diag| diag.gas_volume_residual_m3.abs() <= 1.0e-18),
        "pressure consistency diagnostics must record negligible volume residuals"
    );
}

#[test]
fn resolved_injection_rejects_under_resolved_orifice() {
    let mut phi = vec![1.0f64; 2];
    let sparger = vec![true, false];
    let solid = vec![false, false];
    let mut ledger = SpargerGasLedger::default();
    let err = apply_resolved_gas_injection(
        &mut phi,
        &sparger,
        &solid,
        ResolvedGasInjectionSpec {
            gas_volumetric_flow_m3_per_s: 1.0e-9,
            dt_s: 1.0,
            dx_m: 1.0e-3,
            orifice_diameter_m: 2.9e-3,
        },
        &mut ledger,
    )
    .unwrap_err();
    assert!(err.message.contains("orifice_diameter_m / dx"));
}
