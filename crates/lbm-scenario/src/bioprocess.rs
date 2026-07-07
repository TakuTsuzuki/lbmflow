use crate::units::UnitReport;
use crate::UnsupportedReason;
use lbm_core::geometry as core_geometry;
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

        self.validate_reactor_geometry()?;

        Ok(())
    }

    pub fn compute_unit_report(&self) -> Result<UnitReport, BioprocessScenarioError> {
        crate::units::bioprocess_unit_report(self)
    }

    pub fn unit_report_with_diagnostics(&self) -> Result<UnitReport, BioprocessScenarioError> {
        crate::units::bioprocess_unit_report_unchecked(self)
    }

    pub fn build_geometry(
        &self,
    ) -> Result<core_geometry::StirredTankGeometry, BioprocessScenarioError> {
        let unit_report = self.unit_report_with_diagnostics()?;
        self.reactor
            .build_geometry(&self.run, unit_report.lattice.dx_m)
    }

    pub fn import_stl_geometry(
        &self,
    ) -> Result<core_geometry::StirredTankGeometry, BioprocessScenarioError> {
        self.reactor.import_stl_geometry(
            self.credibility_tier,
            &self.run,
            self.unit_report_with_diagnostics()?.lattice.dx_m,
        )
    }

    fn validate_reactor_geometry(&self) -> Result<(), BioprocessScenarioError> {
        self.reactor.validate_geometry(&self.run)
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

    fn from_geometry(error: core_geometry::GeometryError) -> Self {
        Self {
            message: error.message,
            reason: error.reason,
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
        #[serde(default)]
        bottom: TankBottomSpec,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        stl_import: Option<StlImportSpec>,
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

    fn validate_geometry(&self, run: &RunSpec) -> Result<(), BioprocessScenarioError> {
        match self {
            ReactorSpec::StirredTank {
                bottom, spargers, ..
            } => {
                if *bottom == TankBottomSpec::Dished {
                    return Err(BioprocessScenarioError::unsupported(
                        "dished-bottom stirred tanks are not implemented for M0",
                        UnsupportedReason::NotImplemented,
                    ));
                }
                if run.grid_nx < 3 || run.grid_ny < 3 || run.grid_nz < 3 {
                    return Err(BioprocessScenarioError::unsupported(
                        "stirred-tank grid dimensions must preserve a 1-cell solid rim",
                        UnsupportedReason::OutOfValidityRange {
                            detail: "run.grid_nx, run.grid_ny and run.grid_nz must each be >= 3"
                                .to_string(),
                        },
                    ));
                }
                for sparger in spargers {
                    if !sparger.raw_phi_boundary_fields().is_empty() {
                        return Err(BioprocessScenarioError::unsupported(
                            "raw phi boundary fields are not accepted for sparger geometry",
                            UnsupportedReason::OutOfValidityRange {
                                detail:
                                    "spargers declare geometry and gas metadata only; BCFD-046 owns gas injection"
                                        .to_string(),
                            },
                        ));
                    }
                }
            }
        }
        Ok(())
    }

    fn build_geometry(
        &self,
        run: &RunSpec,
        dx_m: f64,
    ) -> Result<core_geometry::StirredTankGeometry, BioprocessScenarioError> {
        match self {
            ReactorSpec::StirredTank {
                vessel_diameter_m,
                liquid_height_m,
                bottom,
                impellers,
                baffles,
                spargers,
                ..
            } => core_geometry::build_stirred_tank_geometry(
                grid_spec(run, dx_m),
                core_geometry::TankSpec {
                    vessel_diameter_m: *vessel_diameter_m,
                    liquid_height_m: *liquid_height_m,
                    bottom: (*bottom).into(),
                },
                &baffles.iter().map(Into::into).collect::<Vec<_>>(),
                &impellers.iter().map(Into::into).collect::<Vec<_>>(),
                &spargers.iter().map(Into::into).collect::<Vec<_>>(),
            )
            .map_err(BioprocessScenarioError::from_geometry),
        }
    }

    fn import_stl_geometry(
        &self,
        tier: CredibilityTier,
        run: &RunSpec,
        dx_m: f64,
    ) -> Result<core_geometry::StirredTankGeometry, BioprocessScenarioError> {
        match self {
            ReactorSpec::StirredTank {
                vessel_diameter_m,
                liquid_height_m,
                stl_import,
                ..
            } => {
                let Some(import) = stl_import else {
                    return self.build_geometry(run, dx_m);
                };
                #[cfg(not(feature = "geometry-import"))]
                {
                    let _ = (tier, run, dx_m, vessel_diameter_m, liquid_height_m, import);
                    Err(BioprocessScenarioError::unsupported(
                        "STL geometry import requires the geometry-import feature",
                        UnsupportedReason::NotImplemented,
                    ))
                }
                #[cfg(feature = "geometry-import")]
                {
                    let labels = import
                        .labels_path
                        .clone()
                        .unwrap_or_else(|| import.path.with_extension("json"));
                    let imported = lbm_core::voxel_import::import_binary_stl_with_labels(
                        &import.path,
                        &labels,
                        lbm_core::voxel_import::VoxelImportOptions {
                            dims: [
                                run.grid_nx as usize,
                                run.grid_ny as usize,
                                run.grid_nz as usize,
                            ],
                            dx_m,
                            credibility_tier: tier.into(),
                        },
                    )
                    .map_err(BioprocessScenarioError::from_geometry)?;
                    let n = (run.grid_nx * run.grid_ny * run.grid_nz) as usize;
                    Ok(core_geometry::StirredTankGeometry {
                        dims: [
                            run.grid_nx as usize,
                            run.grid_ny as usize,
                            run.grid_nz as usize,
                        ],
                        dx_m,
                        solid: imported.solid,
                        wall_velocity: vec![[0.0; 3]; n],
                        baffle_mask: imported.baffle_mask,
                        sparger_mask: imported.sparger_mask,
                        sparger_orifice_centers: Vec::new(),
                        impellers: Vec::new(),
                    })
                }
            }
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum TankBottomSpec {
    Flat,
    Dished,
}

impl Default for TankBottomSpec {
    fn default() -> Self {
        Self::Flat
    }
}

impl From<TankBottomSpec> for core_geometry::TankBottom {
    fn from(value: TankBottomSpec) -> Self {
        match value {
            TankBottomSpec::Flat => Self::Flat,
            TankBottomSpec::Dished => Self::Dished,
        }
    }
}

#[cfg(feature = "geometry-import")]
impl From<CredibilityTier> for lbm_core::voxel_import::CredibilityTier {
    fn from(value: CredibilityTier) -> Self {
        match value {
            CredibilityTier::Screening => Self::Screening,
            CredibilityTier::Engineering => Self::Engineering,
            CredibilityTier::Evidence => Self::Evidence,
        }
    }
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case", deny_unknown_fields)]
pub struct StlImportSpec {
    pub path: PathBuf,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub labels_path: Option<PathBuf>,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum ImpellerSpec {
    Rushton {
        diameter_m: f64,
        clearance_from_bottom_m: f64,
        rotational_speed_rpm: f64,
        blade_count: u32,
    },
    PitchedBlade {
        diameter_m: f64,
        clearance_from_bottom_m: f64,
        rotational_speed_rpm: f64,
        blade_count: u32,
    },
    Marine {
        diameter_m: f64,
        clearance_from_bottom_m: f64,
        rotational_speed_rpm: f64,
        blade_count: u32,
    },
    CustomMarkerSet {
        markers: Vec<[f64; 3]>,
        rotational_speed_rpm: f64,
    },
}

impl ImpellerSpec {
    pub fn diameter_m(&self) -> Option<f64> {
        match self {
            Self::Rushton { diameter_m, .. }
            | Self::PitchedBlade { diameter_m, .. }
            | Self::Marine { diameter_m, .. } => Some(*diameter_m),
            Self::CustomMarkerSet { .. } => None,
        }
    }

    pub fn rotational_speed_rpm(&self) -> f64 {
        match self {
            Self::Rushton {
                rotational_speed_rpm,
                ..
            }
            | Self::PitchedBlade {
                rotational_speed_rpm,
                ..
            }
            | Self::Marine {
                rotational_speed_rpm,
                ..
            }
            | Self::CustomMarkerSet {
                rotational_speed_rpm,
                ..
            } => *rotational_speed_rpm,
        }
    }
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
        #[serde(default, skip_serializing_if = "Vec::is_empty")]
        raw_phi_boundary_fields: Vec<String>,
        gas_volumetric_flow_m3_per_s: Option<f64>,
        vvm: Option<f64>,
        inlet_phase: InletPhase,
    },
    Pipe {
        center_z_m: f64,
        length_m: f64,
        diameter_m: f64,
        axis: PipeAxisSpec,
        orifice_count: u32,
        orifice_diameter_m: f64,
        #[serde(default, skip_serializing_if = "Vec::is_empty")]
        raw_phi_boundary_fields: Vec<String>,
        gas_volumetric_flow_m3_per_s: Option<f64>,
        vvm: Option<f64>,
        inlet_phase: InletPhase,
    },
    PointOrifices {
        center_z_m: f64,
        positions: Vec<[f64; 3]>,
        orifice_diameter_m: f64,
        #[serde(default, skip_serializing_if = "Vec::is_empty")]
        raw_phi_boundary_fields: Vec<String>,
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

    fn raw_phi_boundary_fields(&self) -> &[String] {
        match self {
            SpargerSpec::Ring {
                raw_phi_boundary_fields,
                ..
            }
            | SpargerSpec::Pipe {
                raw_phi_boundary_fields,
                ..
            }
            | SpargerSpec::PointOrifices {
                raw_phi_boundary_fields,
                ..
            } => raw_phi_boundary_fields,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum PipeAxisSpec {
    X,
    Y,
}

impl From<&BaffleSpec> for core_geometry::BaffleTemplate {
    fn from(value: &BaffleSpec) -> Self {
        Self {
            count: value.count,
            width_m: value.width_m,
            thickness_m: value.thickness_m,
            wall_attached: value.wall_attached,
        }
    }
}

impl From<&ImpellerSpec> for core_geometry::ImpellerTemplate {
    fn from(value: &ImpellerSpec) -> Self {
        match value {
            ImpellerSpec::Rushton {
                diameter_m,
                clearance_from_bottom_m,
                rotational_speed_rpm,
                blade_count,
            } => core_geometry::ImpellerTemplate::Parametric {
                kind: core_geometry::ImpellerKind::Rushton,
                diameter_m: *diameter_m,
                clearance_from_bottom_m: *clearance_from_bottom_m,
                rotational_speed_rpm: *rotational_speed_rpm,
                blade_count: *blade_count,
            },
            ImpellerSpec::PitchedBlade {
                diameter_m,
                clearance_from_bottom_m,
                rotational_speed_rpm,
                blade_count,
            } => core_geometry::ImpellerTemplate::Parametric {
                kind: core_geometry::ImpellerKind::PitchedBlade,
                diameter_m: *diameter_m,
                clearance_from_bottom_m: *clearance_from_bottom_m,
                rotational_speed_rpm: *rotational_speed_rpm,
                blade_count: *blade_count,
            },
            ImpellerSpec::Marine {
                diameter_m,
                clearance_from_bottom_m,
                rotational_speed_rpm,
                blade_count,
            } => core_geometry::ImpellerTemplate::Parametric {
                kind: core_geometry::ImpellerKind::Marine,
                diameter_m: *diameter_m,
                clearance_from_bottom_m: *clearance_from_bottom_m,
                rotational_speed_rpm: *rotational_speed_rpm,
                blade_count: *blade_count,
            },
            ImpellerSpec::CustomMarkerSet {
                markers,
                rotational_speed_rpm,
            } => core_geometry::ImpellerTemplate::CustomMarkerSet {
                markers_m: markers.clone(),
                rotational_speed_rpm: *rotational_speed_rpm,
            },
        }
    }
}

impl From<PipeAxisSpec> for core_geometry::PipeAxis {
    fn from(value: PipeAxisSpec) -> Self {
        match value {
            PipeAxisSpec::X => Self::X,
            PipeAxisSpec::Y => Self::Y,
        }
    }
}

impl From<&SpargerSpec> for core_geometry::SpargerTemplate {
    fn from(value: &SpargerSpec) -> Self {
        match value {
            SpargerSpec::Ring {
                center_z_m,
                outer_radius_m,
                orifice_count,
                orifice_diameter_m,
                gas_volumetric_flow_m3_per_s,
                inlet_phase,
                ..
            } => core_geometry::SpargerTemplate::Ring {
                center_z_m: *center_z_m,
                outer_radius_m: *outer_radius_m,
                orifice_count: *orifice_count,
                orifice_diameter_m: *orifice_diameter_m,
                gas_volumetric_flow_m3_per_s: *gas_volumetric_flow_m3_per_s,
                inlet_phase_gas: inlet_phase.is_gas(),
            },
            SpargerSpec::Pipe {
                center_z_m,
                length_m,
                diameter_m,
                axis,
                orifice_count,
                orifice_diameter_m,
                gas_volumetric_flow_m3_per_s,
                inlet_phase,
                ..
            } => core_geometry::SpargerTemplate::Pipe {
                center_z_m: *center_z_m,
                length_m: *length_m,
                diameter_m: *diameter_m,
                axis: (*axis).into(),
                orifice_count: *orifice_count,
                orifice_diameter_m: *orifice_diameter_m,
                gas_volumetric_flow_m3_per_s: *gas_volumetric_flow_m3_per_s,
                inlet_phase_gas: inlet_phase.is_gas(),
            },
            SpargerSpec::PointOrifices {
                center_z_m,
                positions,
                orifice_diameter_m,
                gas_volumetric_flow_m3_per_s,
                inlet_phase,
                ..
            } => core_geometry::SpargerTemplate::PointOrifices {
                center_z_m: *center_z_m,
                positions_m: positions.clone(),
                orifice_diameter_m: *orifice_diameter_m,
                gas_volumetric_flow_m3_per_s: *gas_volumetric_flow_m3_per_s,
                inlet_phase_gas: inlet_phase.is_gas(),
            },
        }
    }
}

fn grid_spec(run: &RunSpec, dx_m: f64) -> core_geometry::GridSpec {
    core_geometry::GridSpec {
        dims: [
            run.grid_nx as usize,
            run.grid_ny as usize,
            run.grid_nz as usize,
        ],
        dx_m,
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

#[cfg(test)]
mod tests {
    use super::*;

    fn base_scenario(reactor: ReactorSpec) -> BioprocessScenario {
        BioprocessScenario {
            version: "bioprocess-1.0".to_string(),
            name: "geometry-boundary-test".to_string(),
            credibility_tier: CredibilityTier::Screening,
            reactor,
            fluids: FluidsSpec {
                liquid_density_kg_m3: 1000.0,
                liquid_viscosity_pa_s: 0.001,
                gas_density_kg_m3: None,
                gas_viscosity_pa_s: None,
                surface_tension_n_m: None,
                oxygen_diffusivity_m2_per_s: None,
                henry_constant: None,
            },
            operation: OperationSpec {
                duration_s: 1.0,
                gas_inlet_temp_c: None,
                initial_condition: InitialCondition::Quiescent,
            },
            physics: PhysicsSpec {
                models: vec![PhysicsModel::SinglePhase],
            },
            cells: None,
            qoi: QoiSpec {
                power: None,
                mixing: None,
                gas_holdup: None,
                bubble_size: None,
                kla: None,
                shear_exposure: None,
                oxygen_exposure: None,
                calibration_dataset_id: None,
                holdout_dataset_id: None,
            },
            run: RunSpec {
                steps: 1,
                dt_s: 1.0e-4,
                grid_nx: 96,
                grid_ny: 96,
                grid_nz: 96,
                backend: Some(BackendSpec::Cpu),
                precision: Some(Precision::F64),
                lattice: Some(LatticeSpec::D3q19),
            },
            outputs: OutputsSpec {
                manifest_path: PathBuf::from("manifest.json"),
                fields_every_n_steps: None,
                probes_every_n_steps: None,
                emit_qoi_json: false,
                emit_qoi_csv: false,
            },
        }
    }

    fn stirred_tank(stl_import: Option<StlImportSpec>, spargers: Vec<SpargerSpec>) -> ReactorSpec {
        ReactorSpec::StirredTank {
            vessel_diameter_m: 1.0,
            liquid_height_m: 1.0,
            working_volume_m3: 0.78539816339,
            bottom: TankBottomSpec::Flat,
            stl_import,
            impellers: vec![ImpellerSpec::Rushton {
                diameter_m: 0.34,
                clearance_from_bottom_m: 0.35,
                rotational_speed_rpm: 120.0,
                blade_count: 6,
            }],
            baffles: vec![],
            spargers,
        }
    }

    #[cfg(not(feature = "geometry-import"))]
    #[test]
    fn feature_off_stl_import_returns_structured_unsupported() {
        let scenario = base_scenario(stirred_tank(
            Some(StlImportSpec {
                path: PathBuf::from("missing.stl"),
                labels_path: None,
            }),
            vec![],
        ));
        let err = scenario.import_stl_geometry().unwrap_err();
        assert_eq!(err.reason, UnsupportedReason::NotImplemented);
        assert!(err.message.contains("geometry-import"));
    }

    #[test]
    fn rejects_raw_phi_boundary_fields() {
        let mut scenario = base_scenario(stirred_tank(
            None,
            vec![SpargerSpec::Ring {
                center_z_m: 0.1,
                outer_radius_m: 0.2,
                orifice_count: 4,
                orifice_diameter_m: 0.04,
                raw_phi_boundary_fields: vec!["phi".to_string()],
                gas_volumetric_flow_m3_per_s: Some(1.0e-5),
                vvm: None,
                inlet_phase: InletPhase::Gas,
            }],
        ));
        scenario.physics = PhysicsSpec {
            models: vec![PhysicsModel::ResolvedPhaseField {
                interface_width_m: 0.05,
                mobility_m2_per_s: 1.0e-8,
                contact_angle_deg: None,
            }],
        };
        let err = scenario.validate().unwrap_err();
        assert!(matches!(
            err.reason,
            UnsupportedReason::OutOfValidityRange { .. }
        ));
    }
}
