//! Oxygen scalar transport helpers and resolved-interface transfer terms.

use crate::kla::KlModel;
use crate::real::Real;
use crate::solver::UnsupportedReason;
use serde::{Deserialize, Serialize};

pub const OXYGEN_SCALAR_NAME: &str = "oxygen";
pub const DEFAULT_O2_PARTIAL_PRESSURE_PA: f64 = 21_000.0;

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct OxygenState {
    pub c_liquid: Vec<f64>,
    pub c_star: Vec<f64>,
}

#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct OxygenDiagnostics {
    pub clipped_cells: usize,
    pub clipped_fraction: f64,
    pub min_before_clip: f64,
}

#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct OxygenFluxLedger {
    pub cumulative_mol: f64,
    pub last_step_mol: f64,
    pub source_integral_mol_s: f64,
    pub active_interface_cells: usize,
}

#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct OxygenError {
    pub code: &'static str,
    pub message: String,
    pub reason: UnsupportedReason,
}

impl OxygenError {
    fn invalid(message: impl Into<String>, detail: impl Into<String>) -> Self {
        Self {
            code: "oxygen_out_of_validity_range",
            message: message.into(),
            reason: UnsupportedReason::OutOfValidityRange {
                detail: detail.into(),
            },
        }
    }

    fn evidence_rejected(message: impl Into<String>) -> Self {
        Self {
            code: "oxygen_evidence_gate_rejected",
            message: message.into(),
            reason: UnsupportedReason::EvidenceGateFailed {
                missing: vec!["calibrated_kL_table".to_string()],
            },
        }
    }
}

impl std::fmt::Display for OxygenError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{} ({:?})", self.message, self.reason)
    }
}

impl std::error::Error for OxygenError {}

pub type OxygenResult<T> = Result<T, OxygenError>;

pub fn henry_equilibrium(henry_constant: f64, partial_pressure_o2_pa: f64) -> OxygenResult<f64> {
    if !(henry_constant.is_finite() && henry_constant > 0.0) {
        return Err(OxygenError::invalid(
            "Henry constant must be finite and positive",
            "fluids.henry_constant must be provided and > 0",
        ));
    }
    if !(partial_pressure_o2_pa.is_finite() && partial_pressure_o2_pa >= 0.0) {
        return Err(OxygenError::invalid(
            "oxygen partial pressure must be finite and non-negative",
            "partial_pressure_o2_pa must be >= 0",
        ));
    }
    Ok(henry_constant * partial_pressure_o2_pa)
}

pub fn interfacial_area_density(
    phi: f64,
    grad_phi: [f64; 3],
    interface_width_m: f64,
) -> OxygenResult<f64> {
    if !(phi.is_finite() && (0.0..=1.0).contains(&phi)) {
        return Err(OxygenError::invalid(
            "phase fraction must be finite and in [0, 1]",
            "resolved-interface oxygen transfer requires bounded phi",
        ));
    }
    if !(interface_width_m.is_finite() && interface_width_m > 0.0) {
        return Err(OxygenError::invalid(
            "interface width must be finite and positive",
            "interface_width_m must be > 0",
        ));
    }
    let grad_abs =
        (grad_phi[0] * grad_phi[0] + grad_phi[1] * grad_phi[1] + grad_phi[2] * grad_phi[2]).sqrt();
    if !grad_abs.is_finite() {
        return Err(OxygenError::invalid(
            "phase gradient magnitude must be finite",
            "grad_phi components must be finite",
        ));
    }
    Ok(6.0 * grad_abs * phi * (1.0 - phi) * (4.0 / interface_width_m))
}

pub fn reject_uncalibrated_kl_for_evidence(model: &KlModel) -> OxygenResult<()> {
    match model {
        KlModel::Calibrated { table_ref, .. } if !table_ref.is_empty() => Ok(()),
        _ => Err(OxygenError::evidence_rejected(
            "resolved oxygen kLa evidence tier requires calibrated kL",
        )),
    }
}

pub fn oxygen_source_step(
    c_liquid: f64,
    c_star: f64,
    area_density_1_m: f64,
    model: &KlModel,
) -> OxygenResult<f64> {
    if !(c_liquid.is_finite()
        && c_star.is_finite()
        && area_density_1_m.is_finite()
        && area_density_1_m >= 0.0)
    {
        return Err(OxygenError::invalid(
            "oxygen transfer source inputs must be finite",
            "C_L and C_star must be finite; interfacial area density must be >= 0",
        ));
    }
    let kl = model.value_m_s().map_err(|e| OxygenError {
        code: e.code,
        message: e.message,
        reason: e.reason,
    })?;
    Ok(kl * area_density_1_m * (c_star - c_liquid))
}

pub fn apply_interfacial_flux_sources(
    c_liquid: &mut [f64],
    c_star: &[f64],
    area_density_1_m: &[f64],
    cell_volumes_m3: &[f64],
    dt_s: f64,
    model: &KlModel,
    ledger: &mut OxygenFluxLedger,
) -> OxygenResult<OxygenDiagnostics> {
    if c_liquid.len() != c_star.len()
        || c_liquid.len() != area_density_1_m.len()
        || c_liquid.len() != cell_volumes_m3.len()
    {
        return Err(OxygenError::invalid(
            "oxygen source arrays must have equal length",
            "C_L, C_star, a_local, and cell volumes must match",
        ));
    }
    if !(dt_s.is_finite() && dt_s >= 0.0) {
        return Err(OxygenError::invalid(
            "oxygen source time step must be finite and non-negative",
            "dt_s must be >= 0",
        ));
    }
    let mut source_integral = 0.0;
    let mut active = 0usize;
    for i in 0..c_liquid.len() {
        if !(cell_volumes_m3[i].is_finite() && cell_volumes_m3[i] > 0.0) {
            return Err(OxygenError::invalid(
                "cell volume must be finite and positive",
                "cell_volumes_m3 entries must be > 0",
            ));
        }
        let source = oxygen_source_step(c_liquid[i], c_star[i], area_density_1_m[i], model)?;
        if area_density_1_m[i] > 0.0 {
            active += 1;
        }
        c_liquid[i] += source * dt_s;
        source_integral += source * cell_volumes_m3[i];
    }
    let diagnostics = clip_negative_concentrations(c_liquid);
    ledger.source_integral_mol_s = source_integral;
    ledger.last_step_mol = source_integral * dt_s;
    ledger.cumulative_mol += ledger.last_step_mol;
    ledger.active_interface_cells = active;
    Ok(diagnostics)
}

pub fn clip_negative_concentrations<T: Real>(values: &mut [T]) -> OxygenDiagnostics {
    let mut clipped = 0usize;
    let mut min_before = f64::INFINITY;
    let len = values.len();
    for v in values.iter_mut() {
        let before = v.as_f64();
        if before < min_before {
            min_before = before;
        }
        if *v < T::zero() {
            *v = T::zero();
            clipped += 1;
        }
    }
    OxygenDiagnostics {
        clipped_cells: clipped,
        clipped_fraction: if len == 0 {
            0.0
        } else {
            clipped as f64 / len as f64
        },
        min_before_clip: if min_before.is_finite() {
            min_before
        } else {
            0.0
        },
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn oxygen_scalar_name_is_stable_for_ade_registration() {
        assert_eq!(OXYGEN_SCALAR_NAME, "oxygen");
    }

    #[test]
    fn uniform_oxygen_remains_uniform_without_sources() {
        let mut c = vec![0.25; 8];
        let c_star = vec![1.0; 8];
        let area = vec![0.0; 8];
        let volume = vec![1.0e-9; 8];
        let mut ledger = OxygenFluxLedger::default();
        let diag = apply_interfacial_flux_sources(
            &mut c,
            &c_star,
            &area,
            &volume,
            1.0,
            &KlModel::Constant { value_m_s: 1.0e-4 },
            &mut ledger,
        )
        .unwrap();
        assert_eq!(c, vec![0.25; 8]);
        assert_eq!(diag.clipped_cells, 0);
        assert_eq!(ledger.last_step_mol, 0.0);
    }

    #[test]
    fn pure_diffusive_flux_points_down_concentration_gradient() {
        let left = 2.0;
        let right = 1.0;
        let diffusivity = 2.0e-9;
        let dx = 1.0e-3;
        let flux_left_to_right = -diffusivity * (right - left) / dx;
        assert!(flux_left_to_right > 0.0);
    }

    #[test]
    fn henry_boundary_concentration_uses_scenario_inputs() {
        let c_star = henry_equilibrium(1.3e-5, DEFAULT_O2_PARTIAL_PRESSURE_PA).unwrap();
        assert!((c_star - 0.273).abs() < 1.0e-12);
    }

    #[test]
    fn negative_concentration_clip_reports_diagnostic() {
        let mut c = vec![1.0, -0.2, 0.0, -1.0];
        let diag = clip_negative_concentrations(&mut c);
        assert_eq!(c, vec![1.0, 0.0, 0.0, 0.0]);
        assert_eq!(diag.clipped_cells, 2);
        assert_eq!(diag.clipped_fraction, 0.5);
        assert_eq!(diag.min_before_clip, -1.0);
    }

    #[test]
    fn interfacial_flux_is_zero_when_area_or_driving_force_is_zero() {
        let model = KlModel::Constant { value_m_s: 2.0e-5 };
        assert_eq!(oxygen_source_step(0.1, 1.0, 0.0, &model).unwrap(), 0.0);
        assert_eq!(oxygen_source_step(0.5, 0.5, 20.0, &model).unwrap(), 0.0);
    }

    #[test]
    fn closed_system_ledger_matches_integral_source() {
        let mut c = vec![0.0, 0.0];
        let c_star = vec![1.0, 1.0];
        let area = vec![10.0, 20.0];
        let volume = vec![2.0, 3.0];
        let mut ledger = OxygenFluxLedger::default();
        apply_interfacial_flux_sources(
            &mut c,
            &c_star,
            &area,
            &volume,
            0.5,
            &KlModel::Constant { value_m_s: 0.01 },
            &mut ledger,
        )
        .unwrap();
        let expected = (0.01 * 10.0 * 2.0 + 0.01 * 20.0 * 3.0) * 0.5;
        assert!((ledger.last_step_mol - expected).abs() < 1.0e-12);
    }

    #[test]
    fn invalid_henry_is_rejected() {
        let err = henry_equilibrium(0.0, DEFAULT_O2_PARTIAL_PRESSURE_PA).unwrap_err();
        assert_eq!(err.code, "oxygen_out_of_validity_range");
    }

    #[test]
    fn evidence_tier_rejects_uncalibrated_kl() {
        let err = reject_uncalibrated_kl_for_evidence(&KlModel::Constant { value_m_s: 1.0e-5 })
            .unwrap_err();
        assert_eq!(err.code, "oxygen_evidence_gate_rejected");
    }
}
