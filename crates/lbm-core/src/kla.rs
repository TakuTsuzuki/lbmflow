//! kLa computation from point-bubble / PBM interfacial area.

use crate::pbm::{PbmBins, PbmError};
use crate::solver::UnsupportedReason;
use serde::{Deserialize, Serialize};

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum KlModel {
    Constant { value_m_s: f64 },
    PenetrationTheoryPlaceholder,
    Calibrated { table_ref: String, value_m_s: f64 },
}

impl KlModel {
    pub fn discriminant(&self) -> &'static str {
        match self {
            Self::Constant { .. } => "constant",
            Self::PenetrationTheoryPlaceholder => "penetration_theory_placeholder",
            Self::Calibrated { .. } => "calibrated",
        }
    }

    pub fn value_m_s(&self) -> KlaResult<f64> {
        match self {
            Self::Constant { value_m_s } | Self::Calibrated { value_m_s, .. } => {
                validate_positive("kL", *value_m_s)?;
                Ok(*value_m_s)
            }
            Self::PenetrationTheoryPlaceholder => Err(KlaError {
                code: "kla_kl_model_not_implemented",
                message: "penetration-theory kL is a placeholder and has no value".to_string(),
                reason: UnsupportedReason::NotImplemented,
            }),
        }
    }

    pub fn evidence_ready(&self) -> KlaResult<()> {
        match self {
            Self::Calibrated { table_ref, .. } if !table_ref.is_empty() => Ok(()),
            _ => Err(KlaError {
                code: "kla_evidence_gate_rejected",
                message: "PBM kLa evidence tier requires calibrated kL".to_string(),
                reason: UnsupportedReason::EvidenceGateFailed {
                    missing: vec!["calibrated_kL_table".to_string()],
                },
            }),
        }
    }
}

#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct KlaProvenance {
    pub method: &'static str,
    pub k_l_source: &'static str,
    pub d32_source: &'static str,
    pub units: &'static str,
}

#[derive(Clone, Debug, PartialEq)]
pub struct KlaCell {
    pub interfacial_area_density_1_m: f64,
    pub kla_1_s: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct KlaReport {
    pub cells: Vec<KlaCell>,
    pub volume_average_kla_1_s: f64,
    pub provenance: KlaProvenance,
}

pub fn interfacial_area_from_alpha_d32(alpha_g: f64, d32_m: f64) -> KlaResult<f64> {
    validate_fraction("alpha_g", alpha_g)?;
    validate_positive("d32_m", d32_m)?;
    Ok(6.0 * alpha_g / d32_m)
}

pub fn oxygen_transfer_rate_mol_m3_s(
    kla_1_s: f64,
    c_liquid_mol_m3: f64,
    c_star_mol_m3: f64,
) -> KlaResult<f64> {
    if !(kla_1_s.is_finite()
        && kla_1_s >= 0.0
        && c_liquid_mol_m3.is_finite()
        && c_star_mol_m3.is_finite())
    {
        return Err(KlaError::out_of_validity_range(
            "oxygen transfer requires finite kLa and concentrations",
            "kla_1_s must be >= 0; concentrations must be finite",
        ));
    }
    Ok(kla_1_s * (c_star_mol_m3 - c_liquid_mol_m3))
}

pub fn compute_kla_from_alpha_d32(
    alpha_g: &[f64],
    d32_m: &[f64],
    cell_volumes_m3: &[f64],
    model: &KlModel,
) -> KlaResult<KlaReport> {
    if alpha_g.len() != d32_m.len() || alpha_g.len() != cell_volumes_m3.len() {
        return Err(KlaError::out_of_validity_range(
            "kLa field lengths must match",
            "alpha_g, d32_m and cell_volumes_m3 must have equal length",
        ));
    }
    let kl = model.value_m_s()?;
    let mut cells = Vec::with_capacity(alpha_g.len());
    let mut weighted = 0.0;
    let mut volume = 0.0;
    for ((&alpha, &d32), &vol) in alpha_g.iter().zip(d32_m.iter()).zip(cell_volumes_m3.iter()) {
        validate_positive("cell_volume_m3", vol)?;
        let a = interfacial_area_from_alpha_d32(alpha, d32)?;
        let kla = kl * a;
        cells.push(KlaCell {
            interfacial_area_density_1_m: a,
            kla_1_s: kla,
        });
        weighted += kla * vol;
        volume += vol;
    }
    if volume <= 0.0 {
        return Err(KlaError::out_of_validity_range(
            "kLa averaging volume must be positive",
            "sum(cell_volumes_m3) must be > 0",
        ));
    }
    Ok(KlaReport {
        cells,
        volume_average_kla_1_s: weighted / volume,
        provenance: KlaProvenance {
            method: "point_bubble_pbm",
            k_l_source: model.discriminant(),
            d32_source: "pbm",
            units: "1/s",
        },
    })
}

pub fn compute_kla_from_pbm_bins(
    bins: &[PbmBins],
    cell_volumes_m3: &[f64],
    model: &KlModel,
) -> KlaResult<KlaReport> {
    if bins.len() != cell_volumes_m3.len() {
        return Err(KlaError::out_of_validity_range(
            "PBM and volume arrays must have equal length",
            "bins.len() must equal cell_volumes_m3.len()",
        ));
    }
    let kl = model.value_m_s()?;
    let mut cells = Vec::with_capacity(bins.len());
    let mut weighted = 0.0;
    let mut volume = 0.0;
    for (bin, &vol) in bins.iter().zip(cell_volumes_m3.iter()) {
        validate_positive("cell_volume_m3", vol)?;
        let a = bin.interfacial_area_density_1_m();
        let kla = kl * a;
        cells.push(KlaCell {
            interfacial_area_density_1_m: a,
            kla_1_s: kla,
        });
        weighted += kla * vol;
        volume += vol;
    }
    if volume <= 0.0 {
        return Err(KlaError::out_of_validity_range(
            "kLa averaging volume must be positive",
            "sum(cell_volumes_m3) must be > 0",
        ));
    }
    Ok(KlaReport {
        cells,
        volume_average_kla_1_s: weighted / volume,
        provenance: KlaProvenance {
            method: "point_bubble_pbm",
            k_l_source: model.discriminant(),
            d32_source: "pbm",
            units: "1/s",
        },
    })
}

#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct KlaError {
    pub code: &'static str,
    pub message: String,
    pub reason: UnsupportedReason,
}

impl KlaError {
    pub fn out_of_validity_range(message: impl Into<String>, detail: impl Into<String>) -> Self {
        Self {
            code: "kla_out_of_validity_range",
            message: message.into(),
            reason: UnsupportedReason::OutOfValidityRange {
                detail: detail.into(),
            },
        }
    }
}

impl From<PbmError> for KlaError {
    fn from(value: PbmError) -> Self {
        Self {
            code: "pbm_error",
            message: value.message,
            reason: value.reason,
        }
    }
}

impl std::fmt::Display for KlaError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{} ({:?})", self.message, self.reason)
    }
}

impl std::error::Error for KlaError {}

pub type KlaResult<T> = Result<T, KlaError>;

fn validate_positive(name: &'static str, value: f64) -> KlaResult<()> {
    if value.is_finite() && value > 0.0 {
        Ok(())
    } else {
        Err(KlaError::out_of_validity_range(
            format!("{name} must be finite and positive"),
            format!("{name} must be > 0"),
        ))
    }
}

fn validate_fraction(name: &'static str, value: f64) -> KlaResult<()> {
    if value.is_finite() && value >= 0.0 && value <= 1.0 {
        Ok(())
    } else {
        Err(KlaError::out_of_validity_range(
            format!("{name} must be a finite fraction"),
            format!("{name} must be in [0, 1]"),
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn synthetic_d32_alpha_gives_expected_area() {
        let a = interfacial_area_from_alpha_d32(0.05, 2.0e-3).unwrap();
        assert!((a - 150.0).abs() < 1.0e-12);
    }

    #[test]
    fn equilibrium_concentration_gives_zero_transfer() {
        let r = oxygen_transfer_rate_mol_m3_s(0.01, 0.2, 0.2).unwrap();
        assert_eq!(r, 0.0);
    }

    #[test]
    fn kla_units_and_volume_average_are_reported() {
        let report = compute_kla_from_alpha_d32(
            &[0.05, 0.10],
            &[2.0e-3, 2.0e-3],
            &[1.0, 3.0],
            &KlModel::Constant { value_m_s: 1.0e-4 },
        )
        .unwrap();
        assert_eq!(report.provenance.units, "1/s");
        assert_eq!(report.provenance.method, "point_bubble_pbm");
        assert!(report.volume_average_kla_1_s > report.cells[0].kla_1_s);
    }

    #[test]
    fn evidence_tier_rejects_constant_kl() {
        let err = KlModel::Constant { value_m_s: 1.0e-4 }
            .evidence_ready()
            .unwrap_err();
        assert_eq!(err.code, "kla_evidence_gate_rejected");
    }
}
