use crate::units::UnitReport;
use crate::UnsupportedReason;
use serde::de::{self, Deserializer, Visitor};
use serde::{Deserialize, Serialize, Serializer};
use std::fmt;
use std::path::PathBuf;

#[derive(Clone, Debug, PartialEq, Serialize)]
#[serde(rename_all = "snake_case")]
pub struct BioprocessScenario {
    pub version: String,
    pub name: String,
    pub credibility_tier: CredibilityTier,
    pub reactor: ReactorSpec,
    pub fluids: FluidsSpec,
    pub operation: OperationSpec,
    pub physics: PhysicsSpec,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub cells: Option<CellsSpec>,
    pub qoi: QoiSpec,
    pub run: RunSpec,
    pub outputs: OutputsSpec,
}

impl BioprocessScenario {
    pub fn from_json_str(text: &str) -> Result<Self, BioprocessScenarioError> {
        let fields: BioprocessScenarioFields =
            serde_json::from_str(text).map_err(BioprocessScenarioError::serde)?;
        Self::from_fields(fields)
    }

    pub fn validate(&self) -> Result<(), BioprocessScenarioError> {
        if self.qoi.kla.is_some() && !self.physics.has_oxygen() {
            return Err(BioprocessScenarioError::unsupported(
                "kLa QOI requires oxygen physics",
                UnsupportedReason::MissingDependency {
                    depends_on: "oxygen".to_string(),
                },
            ));
        }

        if self.reactor.has_spargers() && !self.physics.has_gas_model() {
            return Err(BioprocessScenarioError::unsupported(
                "spargers require a gas model in physics",
                UnsupportedReason::MissingDependency {
                    depends_on: "gas_model_in_physics".to_string(),
                },
            ));
        }

        if self.credibility_tier == CredibilityTier::Evidence
            && (self.qoi.calibration_dataset_id.is_none() || self.qoi.holdout_dataset_id.is_none())
        {
            return Err(BioprocessScenarioError::unsupported(
                "evidence tier requires calibration and holdout dataset registry references",
                UnsupportedReason::EvidenceGateFailed {
                    missing: vec!["calibration_and_holdout_dataset_registry".to_string()],
                },
            ));
        }

        if let Some(sparger) = self.reactor.non_gas_sparger() {
            return Err(BioprocessScenarioError::unsupported(
                format!(
                    "sparger inlet phase must be gas, got {}",
                    sparger.inlet_phase()
                ),
                UnsupportedReason::OutOfValidityRange {
                    detail: "sparger inlet_phase must be gas".to_string(),
                },
            ));
        }

        if self.run.lattice == Some(LatticeSpec::D2q9) {
            return Err(BioprocessScenarioError::unsupported(
                "bioprocess stirred-tank scenarios require a 3D lattice",
                UnsupportedReason::NotImplemented,
            ));
        }

        Ok(())
    }

    pub fn compute_unit_report(&self) -> Result<UnitReport, BioprocessScenarioError> {
        crate::units::bioprocess_unit_report(self)
    }

    pub fn unit_report_with_diagnostics(&self) -> Result<UnitReport, BioprocessScenarioError> {
        crate::units::bioprocess_unit_report_unchecked(self)
    }

    fn from_fields(fields: BioprocessScenarioFields) -> Result<Self, BioprocessScenarioError> {
        let scenario = Self {
            version: fields.version,
            name: fields.name,
            credibility_tier: fields.credibility_tier,
            reactor: fields.reactor,
            fluids: fields.fluids,
            operation: fields.operation,
            physics: fields.physics,
            cells: fields.cells,
            qoi: fields.qoi,
            run: fields.run,
            outputs: fields.outputs,
        };
        scenario.validate()?;
        Ok(scenario)
    }
}

impl<'de> Deserialize<'de> for BioprocessScenario {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        let fields = BioprocessScenarioFields::deserialize(deserializer)?;
        Self::from_fields(fields).map_err(de::Error::custom)
    }
}

#[derive(Clone, Debug, Deserialize)]
#[serde(rename_all = "snake_case", deny_unknown_fields)]
struct BioprocessScenarioFields {
    version: String,
    name: String,
    credibility_tier: CredibilityTier,
    reactor: ReactorSpec,
    fluids: FluidsSpec,
    operation: OperationSpec,
    physics: PhysicsSpec,
    cells: Option<CellsSpec>,
    qoi: QoiSpec,
    run: RunSpec,
    outputs: OutputsSpec,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub struct BioprocessScenarioError {
    pub message: String,
    pub reason: UnsupportedReason,
}

impl BioprocessScenarioError {
    fn serde(error: serde_json::Error) -> Self {
        Self {
            message: error.to_string(),
            reason: UnsupportedReason::OutOfValidityRange {
                detail: "invalid bioprocess scenario JSON".to_string(),
            },
        }
    }

    pub(crate) fn unsupported(message: impl Into<String>, reason: UnsupportedReason) -> Self {
        Self {
            message: message.into(),
            reason,
        }
    }
}

impl fmt::Display for BioprocessScenarioError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} ({:?})", self.message, self.reason)
    }
}

impl std::error::Error for BioprocessScenarioError {}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum CredibilityTier {
    Screening,
    Engineering,
    Evidence,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum ReactorSpec {
    StirredTank {
        vessel_diameter_m: f64,
        liquid_height_m: f64,
        working_volume_m3: f64,
        impellers: Vec<ImpellerSpec>,
        baffles: Vec<BaffleSpec>,
        spargers: Vec<SpargerSpec>,
    },
}

impl ReactorSpec {
    fn has_spargers(&self) -> bool {
        match self {
            ReactorSpec::StirredTank { spargers, .. } => !spargers.is_empty(),
        }
    }

    fn non_gas_sparger(&self) -> Option<&SpargerSpec> {
        match self {
            ReactorSpec::StirredTank { spargers, .. } => spargers
                .iter()
                .find(|sparger| !sparger.inlet_phase().is_gas()),
        }
    }
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case", deny_unknown_fields)]
pub struct ImpellerSpec {
    pub kind: ImpellerKind,
    pub diameter_m: f64,
    pub clearance_from_bottom_m: f64,
    pub rotational_speed_rpm: f64,
    pub blade_count: u32,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ImpellerKind {
    Rushton,
    PitchedBlade,
    Marine,
    CustomMarkerSet,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case", deny_unknown_fields)]
pub struct BaffleSpec {
    pub count: u32,
    pub width_m: f64,
    pub thickness_m: f64,
    pub wall_attached: bool,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum SpargerSpec {
    Ring {
        center_z_m: f64,
        outer_radius_m: f64,
        orifice_count: u32,
        orifice_diameter_m: f64,
        gas_volumetric_flow_m3_per_s: Option<f64>,
        vvm: Option<f64>,
        inlet_phase: InletPhase,
    },
    Pipe {
        center_z_m: f64,
        length_m: f64,
        diameter_m: f64,
        gas_volumetric_flow_m3_per_s: Option<f64>,
        vvm: Option<f64>,
        inlet_phase: InletPhase,
    },
    PointOrifices {
        center_z_m: f64,
        positions: Vec<[f64; 3]>,
        gas_volumetric_flow_m3_per_s: Option<f64>,
        vvm: Option<f64>,
        inlet_phase: InletPhase,
    },
}

impl SpargerSpec {
    fn inlet_phase(&self) -> &InletPhase {
        match self {
            SpargerSpec::Ring { inlet_phase, .. }
            | SpargerSpec::Pipe { inlet_phase, .. }
            | SpargerSpec::PointOrifices { inlet_phase, .. } => inlet_phase,
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum InletPhase {
    Gas,
    Other(String),
}

impl InletPhase {
    fn is_gas(&self) -> bool {
        matches!(self, InletPhase::Gas)
    }
}

impl fmt::Display for InletPhase {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            InletPhase::Gas => write!(f, "gas"),
            InletPhase::Other(value) => write!(f, "{value}"),
        }
    }
}

impl Serialize for InletPhase {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        serializer.serialize_str(&self.to_string())
    }
}

impl<'de> Deserialize<'de> for InletPhase {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        struct InletPhaseVisitor;

        impl<'de> Visitor<'de> for InletPhaseVisitor {
            type Value = InletPhase;

            fn expecting(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
                formatter.write_str("a sparger inlet phase string")
            }

            fn visit_str<E>(self, value: &str) -> Result<Self::Value, E>
            where
                E: de::Error,
            {
                Ok(if value == "gas" {
                    InletPhase::Gas
                } else {
                    InletPhase::Other(value.to_string())
                })
            }
        }

        deserializer.deserialize_str(InletPhaseVisitor)
    }
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case", deny_unknown_fields)]
pub struct FluidsSpec {
    pub liquid_density_kg_m3: f64,
    pub liquid_viscosity_pa_s: f64,
    pub gas_density_kg_m3: Option<f64>,
    pub gas_viscosity_pa_s: Option<f64>,
    pub surface_tension_n_m: Option<f64>,
    pub oxygen_diffusivity_m2_per_s: Option<f64>,
    pub henry_constant: Option<f64>,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case", deny_unknown_fields)]
pub struct OperationSpec {
    pub duration_s: f64,
    pub gas_inlet_temp_c: Option<f64>,
    pub initial_condition: InitialCondition,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum InitialCondition {
    Quiescent,
    ExistingCheckpoint { path: PathBuf },
}

#[derive(Clone, Debug, PartialEq, Serialize)]
#[serde(transparent)]
pub struct PhysicsSpec {
    pub models: Vec<PhysicsModel>,
}

impl PhysicsSpec {
    pub fn has_oxygen(&self) -> bool {
        self.models
            .iter()
            .any(|model| matches!(model, PhysicsModel::Oxygen { .. }))
    }

    pub fn has_gas_model(&self) -> bool {
        self.models.iter().any(|model| {
            matches!(
                model,
                PhysicsModel::ResolvedPhaseField { .. }
                    | PhysicsModel::PointBubble { .. }
                    | PhysicsModel::Hybrid { .. }
            )
        })
    }
}

impl<'de> Deserialize<'de> for PhysicsSpec {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        #[derive(Deserialize)]
        #[serde(untagged)]
        enum PhysicsInput {
            One(PhysicsModel),
            Many(Vec<PhysicsModel>),
        }

        match PhysicsInput::deserialize(deserializer)? {
            PhysicsInput::One(model) => Ok(Self {
                models: vec![model],
            }),
            PhysicsInput::Many(models) => Ok(Self { models }),
        }
    }
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum PhysicsModel {
    SinglePhase,
    ResolvedPhaseField {
        interface_width_m: f64,
        mobility_m2_per_s: f64,
        contact_angle_deg: Option<f64>,
    },
    PointBubble {
        max_bubble_count: u32,
    },
    Hybrid {
        phase_field: ResolvedPhaseFieldInner,
        point_bubble: PointBubbleInner,
    },
    PassiveScalar {
        diffusivity_m2_per_s: f64,
        initial_pulse: Option<PulseSpec>,
    },
    Oxygen {
        henry_constant: f64,
        interfacial_flux_model: OxygenFluxModel,
        our_model: OurModel,
    },
    CellTracer {
        count: u32,
        seed: u64,
        record_shear: bool,
        record_oxygen: bool,
        damage_model: Option<DamageModelSpec>,
    },
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case", deny_unknown_fields)]
pub struct ResolvedPhaseFieldInner {
    pub interface_width_m: f64,
    pub mobility_m2_per_s: f64,
    pub contact_angle_deg: Option<f64>,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case", deny_unknown_fields)]
pub struct PointBubbleInner {
    pub max_bubble_count: u32,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case", deny_unknown_fields)]
pub struct PulseSpec {
    pub center_m: [f64; 3],
    pub radius_m: f64,
    pub concentration: f64,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum OxygenFluxModel {
    HenryEquilibrium,
    ConstantKl,
    InterfacialArea,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum OurModel {
    None,
    Constant,
    Monod,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum DamageModelSpec {
    Threshold { threshold_pa: f64, exponent: f64 },
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case", deny_unknown_fields)]
pub struct CellsSpec {
    pub count: u32,
    pub seed: u64,
    pub record_shear: bool,
    pub record_oxygen: bool,
    pub damage_model: Option<DamageModelSpec>,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case", deny_unknown_fields)]
pub struct QoiSpec {
    pub power: Option<PowerQoiOpts>,
    pub mixing: Option<MixingQoiOpts>,
    pub gas_holdup: Option<GasHoldupQoiOpts>,
    pub bubble_size: Option<BubbleSizeQoiOpts>,
    pub kla: Option<KlaQoiOpts>,
    pub shear_exposure: Option<ShearExposureQoiOpts>,
    pub oxygen_exposure: Option<OxygenExposureQoiOpts>,
    pub calibration_dataset_id: Option<String>,
    pub holdout_dataset_id: Option<String>,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case", deny_unknown_fields)]
pub struct PowerQoiOpts {}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case", deny_unknown_fields)]
pub struct MixingQoiOpts {}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case", deny_unknown_fields)]
pub struct GasHoldupQoiOpts {}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case", deny_unknown_fields)]
pub struct BubbleSizeQoiOpts {}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case", deny_unknown_fields)]
pub struct KlaQoiOpts {}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case", deny_unknown_fields)]
pub struct ShearExposureQoiOpts {}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case", deny_unknown_fields)]
pub struct OxygenExposureQoiOpts {}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case", deny_unknown_fields)]
pub struct RunSpec {
    pub steps: u64,
    pub dt_s: f64,
    pub grid_nx: u32,
    pub grid_ny: u32,
    pub grid_nz: u32,
    pub backend: Option<BackendSpec>,
    pub precision: Option<Precision>,
    pub lattice: Option<LatticeSpec>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum BackendSpec {
    Auto,
    Cpu,
    Gpu,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Precision {
    F32,
    F64,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum LatticeSpec {
    D2q9,
    D3q19,
    D3q27,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case", deny_unknown_fields)]
pub struct OutputsSpec {
    pub manifest_path: PathBuf,
    pub fields_every_n_steps: Option<u64>,
    pub probes_every_n_steps: Option<u64>,
    pub emit_qoi_json: bool,
    pub emit_qoi_csv: bool,
}
