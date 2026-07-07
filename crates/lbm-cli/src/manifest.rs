use crate::capabilities::{CapabilityEntry, CapabilityRegistry};
use lbm_scenario::Scenario;
use serde::ser::{Error as SerError, SerializeStruct};
use serde::{Serialize, Serializer};
use sha2::{Digest, Sha256};

pub const MANIFEST_PATH: &str = "manifest.json";

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
pub struct Manifest {
    pub scenario: String,
    pub scenario_hash: String,
    pub manifest_path: String,
    pub bioprocess_schema_version: Option<String>,
    pub backend: BackendId,
    pub lattice: LatticeId,
    pub precision: PrecisionId,
    pub active_models: Vec<ActiveModelTag>,
    pub qoi_methods: Vec<QoiMethodDescriptor>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub unit_report: Option<lbm_scenario::UnitReport>,
    pub capability_report: Vec<CapabilityEntry>,
    pub status: String,
    pub steps_run: u64,
    pub wall_seconds: f64,
    pub mlups: f64,
    pub diagnostics: Diagnostics,
    pub provenance: Provenance,
    pub warnings: Vec<lbm_scenario::Warning>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub units: Option<lbm_scenario::LegacyUnitReport>,
    pub files: Vec<String>,
}

#[allow(dead_code)]
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum BackendId {
    Cpu,
    Gpu,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum LatticeId {
    D2q9,
    D3q19,
    D3q27,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum PrecisionId {
    F32,
    F64,
}

#[allow(dead_code)]
#[derive(Clone, Debug, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum ActiveModelTag {
    SinglePhase,
    ShanChen,
    RotatingIbm,
    Particles,
    PassiveScalar,
    PhaseFieldMixture,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct QoiMethodDescriptor {
    pub qoi: String,
    pub method: String,
    pub input_fields: Vec<String>,
    pub provenance: QoiProvenance,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct QoiProvenance {
    source_fields: Option<Vec<String>>,
    averaging_window: Option<String>,
    units: Option<String>,
    validation_tier: Option<String>,
}

impl QoiProvenance {
    #[allow(dead_code)]
    pub fn new(
        source_fields: Vec<String>,
        averaging_window: impl Into<String>,
        units: impl Into<String>,
        validation_tier: impl Into<String>,
    ) -> Self {
        Self {
            source_fields: Some(source_fields),
            averaging_window: Some(averaging_window.into()),
            units: Some(units.into()),
            validation_tier: Some(validation_tier.into()),
        }
    }

    #[cfg(test)]
    pub fn missing_units_for_test() -> Self {
        Self {
            source_fields: Some(vec!["ux".to_string()]),
            averaging_window: Some("not_applicable".to_string()),
            units: None,
            validation_tier: Some("screening".to_string()),
        }
    }

    #[cfg(test)]
    pub fn missing_averaging_window_for_test() -> Self {
        Self {
            source_fields: Some(vec!["ux".to_string()]),
            averaging_window: None,
            units: Some("m/s".to_string()),
            validation_tier: Some("screening".to_string()),
        }
    }
}

impl Serialize for QoiProvenance {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let source_fields = self
            .source_fields
            .as_ref()
            .ok_or_else(|| S::Error::custom("QoiProvenance.source_fields is mandatory"))?;
        let averaging_window = self
            .averaging_window
            .as_ref()
            .ok_or_else(|| S::Error::custom("QoiProvenance.averaging_window is mandatory"))?;
        let units = self
            .units
            .as_ref()
            .ok_or_else(|| S::Error::custom("QoiProvenance.units is mandatory"))?;
        let validation_tier = self
            .validation_tier
            .as_ref()
            .ok_or_else(|| S::Error::custom("QoiProvenance.validation_tier is mandatory"))?;
        let mut state = serializer.serialize_struct("QoiProvenance", 4)?;
        state.serialize_field("sourceFields", source_fields)?;
        state.serialize_field("averagingWindow", averaging_window)?;
        state.serialize_field("units", units)?;
        state.serialize_field("validationTier", validation_tier)?;
        state.end()
    }
}

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
pub struct Provenance {
    pub backend: lbm_scenario::BackendChoice,
    pub lattice: String,
    pub collision: CollisionProvenance,
    pub precision: lbm_scenario::Precision,
    pub storage: lbm_scenario::StorageSpec,
}

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
pub struct CollisionProvenance {
    pub kind: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub magic: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub omega_shear: Option<f64>,
}

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
pub struct Diagnostics {
    pub total_mass: f64,
    pub max_speed: f64,
    pub tau: f64,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub phase_field: Option<lbm_core::phase_field::PhaseFieldDiagnostics>,
}

pub fn scenario_hash<T: Serialize>(scenario: &T) -> Result<String, serde_json::Error> {
    let bytes = serde_json::to_vec(scenario)?;
    let digest = Sha256::digest(bytes);
    Ok(format!("sha256:{digest:x}"))
}

pub fn capability_report() -> Vec<CapabilityEntry> {
    CapabilityRegistry::new().iter().cloned().collect()
}

pub fn lattice_id(lattice: &str) -> LatticeId {
    match lattice {
        "D3Q19" => LatticeId::D3q19,
        "D3Q27" => LatticeId::D3q27,
        _ => LatticeId::D2q9,
    }
}

pub fn precision_id(precision: lbm_scenario::Precision) -> PrecisionId {
    match precision {
        lbm_scenario::Precision::F32 => PrecisionId::F32,
        lbm_scenario::Precision::F64 => PrecisionId::F64,
    }
}

pub fn active_models_for_legacy(sc: &Scenario) -> Vec<ActiveModelTag> {
    let mut out = Vec::new();
    if sc.multiphase.is_some() {
        out.push(ActiveModelTag::ShanChen);
    } else {
        out.push(ActiveModelTag::SinglePhase);
    }
    if sc.rotor.is_some() {
        out.push(ActiveModelTag::RotatingIbm);
    }
    if sc.particles.is_some() {
        out.push(ActiveModelTag::Particles);
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn diagnostics() -> Diagnostics {
        Diagnostics {
            total_mass: 1.0,
            max_speed: 0.0,
            tau: 1.0,
            phase_field: None,
        }
    }

    fn provenance() -> Provenance {
        Provenance {
            backend: lbm_scenario::BackendChoice::Cpu,
            lattice: "D3Q19".to_string(),
            collision: CollisionProvenance {
                kind: "trt".to_string(),
                magic: Some(lbm_core::params::CollisionKind::MAGIC_STD),
                omega_shear: None,
            },
            precision: lbm_scenario::Precision::F64,
            storage: lbm_scenario::StorageSpec::F32,
        }
    }

    #[test]
    fn manifest_snapshot_for_bioprocess_scenario_has_bioprocess_fields() {
        let manifest = Manifest {
            scenario: "bio-snapshot".to_string(),
            scenario_hash: "sha256:abc".to_string(),
            manifest_path: MANIFEST_PATH.to_string(),
            bioprocess_schema_version: Some("bioprocess-1.0".to_string()),
            backend: BackendId::Cpu,
            lattice: LatticeId::D3q19,
            precision: PrecisionId::F64,
            active_models: vec![ActiveModelTag::SinglePhase, ActiveModelTag::PassiveScalar],
            qoi_methods: vec![QoiMethodDescriptor {
                qoi: "mixing_time".to_string(),
                method: "scalar_cv_threshold".to_string(),
                input_fields: vec!["C:tracer".to_string()],
                provenance: QoiProvenance::new(
                    vec!["C:tracer".to_string()],
                    "t0..t_end",
                    "s",
                    "screening",
                ),
            }],
            unit_report: None,
            capability_report: capability_report(),
            status: "completed".to_string(),
            steps_run: 0,
            wall_seconds: 0.0,
            mlups: 0.0,
            diagnostics: diagnostics(),
            provenance: provenance(),
            warnings: Vec::new(),
            units: None,
            files: Vec::new(),
        };
        let value = serde_json::to_value(&manifest).unwrap();
        assert_eq!(value["bioprocessSchemaVersion"], "bioprocess-1.0");
        assert_eq!(value["activeModels"][1], "passive_scalar");
        assert_eq!(value["qoiMethods"][0]["provenance"]["units"], "s");
        assert!(value["capabilityReport"].as_array().unwrap().len() >= 8);
    }

    #[test]
    fn qoi_provenance_missing_units_fails_serialisation() {
        let err = serde_json::to_string(&QoiProvenance::missing_units_for_test()).unwrap_err();
        assert!(err.to_string().contains("units"));
    }

    #[test]
    fn qoi_provenance_missing_averaging_window_fails_serialisation() {
        let err =
            serde_json::to_string(&QoiProvenance::missing_averaging_window_for_test()).unwrap_err();
        assert!(err.to_string().contains("averaging_window"));
    }
}
