//! Resolved phase-field sparger gas-injection helpers.

use crate::geometry::SPARGER_ORIFICE_MIN_CELLS;
use crate::phase_field::PhaseFieldError;
use crate::real::Real;
use serde::{Deserialize, Serialize};

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ResolvedGasInjectionSpec {
    pub gas_volumetric_flow_m3_per_s: f64,
    pub dt_s: f64,
    pub dx_m: f64,
    pub orifice_diameter_m: f64,
}

impl ResolvedGasInjectionSpec {
    pub fn validate(self) -> Result<Self, PhaseFieldError> {
        if !(self.gas_volumetric_flow_m3_per_s.is_finite()
            && self.gas_volumetric_flow_m3_per_s > 0.0)
        {
            return Err(PhaseFieldError {
                message: "sparger gas volumetric flow must be finite and > 0".to_string(),
            });
        }
        if !(self.dt_s.is_finite() && self.dt_s > 0.0) {
            return Err(PhaseFieldError {
                message: "sparger injection dt_s must be finite and > 0".to_string(),
            });
        }
        if !(self.dx_m.is_finite() && self.dx_m > 0.0) {
            return Err(PhaseFieldError {
                message: "sparger injection dx_m must be finite and > 0".to_string(),
            });
        }
        if !(self.orifice_diameter_m.is_finite() && self.orifice_diameter_m > 0.0) {
            return Err(PhaseFieldError {
                message: "sparger orifice diameter must be finite and > 0".to_string(),
            });
        }
        let cells = self.orifice_diameter_m / self.dx_m;
        if cells < SPARGER_ORIFICE_MIN_CELLS {
            return Err(PhaseFieldError {
                message: format!(
                    "orifice_diameter_m / dx must be >= {SPARGER_ORIFICE_MIN_CELLS} for resolved phase-field gas injection (got {cells:.3})"
                ),
            });
        }
        Ok(self)
    }
}

#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub struct SpargerPressureDiagnostic {
    pub requested_gas_volume_m3: f64,
    pub injected_gas_volume_m3: f64,
    pub gas_volume_residual_m3: f64,
    pub injection_cells: usize,
    pub mean_orifice_alpha_g: f64,
}

#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub struct SpargerGasLedger {
    pub requested_gas_volume_m3: f64,
    pub injected_gas_volume_m3: f64,
    pub gas_volume_residual_m3: f64,
    pub pressure_diagnostics: Vec<SpargerPressureDiagnostic>,
}

impl SpargerGasLedger {
    pub fn record(&mut self, diagnostic: SpargerPressureDiagnostic) {
        self.requested_gas_volume_m3 += diagnostic.requested_gas_volume_m3;
        self.injected_gas_volume_m3 += diagnostic.injected_gas_volume_m3;
        self.gas_volume_residual_m3 = self.requested_gas_volume_m3 - self.injected_gas_volume_m3;
        self.pressure_diagnostics.push(diagnostic);
    }
}

pub fn apply_resolved_gas_injection<T: Real>(
    phi: &mut [T],
    sparger_mask: &[bool],
    solid: &[bool],
    spec: ResolvedGasInjectionSpec,
    ledger: &mut SpargerGasLedger,
) -> Result<SpargerPressureDiagnostic, PhaseFieldError> {
    let spec = spec.validate()?;
    if phi.len() != sparger_mask.len() || phi.len() != solid.len() {
        return Err(PhaseFieldError {
            message: format!(
                "sparger injection arrays must have matching lengths (phi={}, sparger_mask={}, solid={})",
                phi.len(),
                sparger_mask.len(),
                solid.len()
            ),
        });
    }

    let injection_cells = sparger_mask
        .iter()
        .zip(solid)
        .filter(|&(&is_sparger, &is_solid)| is_sparger && !is_solid)
        .count();
    if injection_cells == 0 {
        return Err(PhaseFieldError {
            message: "sparger mask contains no fluid injection cells".to_string(),
        });
    }

    let requested = spec.gas_volumetric_flow_m3_per_s * spec.dt_s;
    let cell_volume = spec.dx_m * spec.dx_m * spec.dx_m;
    let alpha_increment = requested / (injection_cells as f64 * cell_volume);
    if !(alpha_increment.is_finite() && alpha_increment >= 0.0) {
        return Err(PhaseFieldError {
            message: "sparger gas fraction increment is non-finite".to_string(),
        });
    }
    for ((value, &is_sparger), &is_solid) in phi.iter().zip(sparger_mask).zip(solid) {
        if !is_sparger || is_solid {
            continue;
        }
        if value.as_f64() < alpha_increment {
            return Err(PhaseFieldError {
                message: format!(
                    "sparger injection volume exceeds local liquid capacity: required alpha increment {alpha_increment:.6e}, local phi {:.6e}",
                    value.as_f64()
                ),
            });
        }
    }

    let delta = T::r(alpha_increment);
    let mut alpha_sum_after = 0.0;
    for ((value, &is_sparger), &is_solid) in phi.iter_mut().zip(sparger_mask).zip(solid) {
        if !is_sparger || is_solid {
            continue;
        }
        *value = *value - delta;
        alpha_sum_after += 1.0 - value.as_f64();
    }
    let injected = alpha_increment * injection_cells as f64 * cell_volume;
    let diagnostic = SpargerPressureDiagnostic {
        requested_gas_volume_m3: requested,
        injected_gas_volume_m3: injected,
        gas_volume_residual_m3: requested - injected,
        injection_cells,
        mean_orifice_alpha_g: alpha_sum_after / injection_cells as f64,
    };
    ledger.record(diagnostic.clone());
    Ok(diagnostic)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rejects_under_resolved_orifice_using_geometry_floor() {
        let err = ResolvedGasInjectionSpec {
            gas_volumetric_flow_m3_per_s: 1.0e-9,
            dt_s: 1.0,
            dx_m: 1.0e-3,
            orifice_diameter_m: (SPARGER_ORIFICE_MIN_CELLS - 0.1) * 1.0e-3,
        }
        .validate()
        .unwrap_err();
        assert!(err.message.contains("orifice_diameter_m / dx"));
    }

    #[test]
    fn injection_decreases_phi_and_balances_ledger() {
        let mut phi = vec![1.0f64; 8];
        let sparger = vec![true, true, false, false, false, false, false, false];
        let solid = vec![false; 8];
        let mut ledger = SpargerGasLedger::default();
        let diag = apply_resolved_gas_injection(
            &mut phi,
            &sparger,
            &solid,
            ResolvedGasInjectionSpec {
                gas_volumetric_flow_m3_per_s: 2.0e-9,
                dt_s: 1.0,
                dx_m: 1.0e-3,
                orifice_diameter_m: SPARGER_ORIFICE_MIN_CELLS * 1.0e-3,
            },
            &mut ledger,
        )
        .unwrap();
        assert!(phi[0] < 1.0 && phi[1] < 1.0);
        assert_eq!(phi[2], 1.0);
        assert_eq!(diag.injection_cells, 2);
        assert!((ledger.injected_gas_volume_m3 - 2.0e-9).abs() < 1.0e-21);
    }
}
