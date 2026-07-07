//! Bioprocess stirred-tank geometry builders.
//!
//! The compact arrays here use the same global indexing convention as
//! `solver::Solver`: `cell = (z * ny + y) * nx + x`.

use crate::rotating_ibm::{IbmMarker, RotatingBody};
use crate::solver::UnsupportedReason;
use std::f64::consts::TAU;

pub const STIRRED_TANK_MIN_CELLS: f64 = 32.0;
pub const SPARGER_ORIFICE_MIN_CELLS: f64 = 3.0;
pub const RUSHTON_BLADE_THICKNESS_FRACTION: f64 = 0.04;
pub const RUSHTON_DISC_DIAMETER_FRACTION: f64 = 0.65;
pub const PITCHED_BLADE_ANGLE_DEG: f64 = 45.0;
pub const PITCHED_BLADE_THICKNESS_FRACTION: f64 = 0.035;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct GeometryError {
    pub message: String,
    pub reason: UnsupportedReason,
}

impl GeometryError {
    pub fn not_implemented(message: impl Into<String>) -> Self {
        Self {
            message: message.into(),
            reason: UnsupportedReason::NotImplemented,
        }
    }

    pub fn out_of_validity_range(message: impl Into<String>, detail: impl Into<String>) -> Self {
        Self {
            message: message.into(),
            reason: UnsupportedReason::OutOfValidityRange {
                detail: detail.into(),
            },
        }
    }

    pub fn evidence_gate_failed(message: impl Into<String>, missing: Vec<String>) -> Self {
        Self {
            message: message.into(),
            reason: UnsupportedReason::EvidenceGateFailed { missing },
        }
    }
}

impl std::fmt::Display for GeometryError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{} ({:?})", self.message, self.reason)
    }
}

impl std::error::Error for GeometryError {}

pub type GeometryResult<T> = Result<T, GeometryError>;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TankBottom {
    Flat,
    Dished,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GridSpec {
    pub dims: [usize; 3],
    pub dx_m: f64,
}

impl GridSpec {
    pub fn cell_count(self) -> usize {
        self.dims[0] * self.dims[1] * self.dims[2]
    }

    fn validate(self) -> GeometryResult<()> {
        if self.dims.iter().any(|&n| n < 3) {
            return Err(GeometryError::out_of_validity_range(
                "stirred-tank grid dimensions must preserve a 1-cell solid rim",
                "run.grid_nx, run.grid_ny and run.grid_nz must each be >= 3",
            ));
        }
        if !(self.dx_m.is_finite() && self.dx_m > 0.0) {
            return Err(GeometryError::out_of_validity_range(
                "grid spacing must be finite and positive",
                "dx_m must be finite and > 0",
            ));
        }
        Ok(())
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TankSpec {
    pub vessel_diameter_m: f64,
    pub liquid_height_m: f64,
    pub bottom: TankBottom,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BaffleTemplate {
    pub count: u32,
    pub width_m: f64,
    pub thickness_m: f64,
    pub wall_attached: bool,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum ImpellerKind {
    Rushton,
    PitchedBlade,
    Marine,
    CustomMarkerSet,
}

#[derive(Clone, Debug, PartialEq)]
pub enum ImpellerTemplate {
    Parametric {
        kind: ImpellerKind,
        diameter_m: f64,
        clearance_from_bottom_m: f64,
        rotational_speed_rpm: f64,
        blade_count: u32,
    },
    CustomMarkerSet {
        markers_m: Vec<[f64; 3]>,
        rotational_speed_rpm: f64,
    },
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PipeAxis {
    X,
    Y,
}

#[derive(Clone, Debug, PartialEq)]
pub enum SpargerTemplate {
    Ring {
        center_z_m: f64,
        outer_radius_m: f64,
        orifice_count: u32,
        orifice_diameter_m: f64,
        gas_volumetric_flow_m3_per_s: Option<f64>,
        inlet_phase_gas: bool,
    },
    Pipe {
        center_z_m: f64,
        length_m: f64,
        diameter_m: f64,
        axis: PipeAxis,
        orifice_count: u32,
        orifice_diameter_m: f64,
        gas_volumetric_flow_m3_per_s: Option<f64>,
        inlet_phase_gas: bool,
    },
    PointOrifices {
        center_z_m: f64,
        positions_m: Vec<[f64; 3]>,
        orifice_diameter_m: f64,
        gas_volumetric_flow_m3_per_s: Option<f64>,
        inlet_phase_gas: bool,
    },
}

#[derive(Clone, Debug, PartialEq)]
pub struct StirredTankGeometry {
    pub dims: [usize; 3],
    pub dx_m: f64,
    pub solid: Vec<bool>,
    pub wall_velocity: Vec<[f64; 3]>,
    pub baffle_mask: Vec<bool>,
    pub sparger_mask: Vec<bool>,
    pub sparger_orifice_centers: Vec<[f64; 3]>,
    pub impellers: Vec<ImpellerMarkerSet>,
}

impl StirredTankGeometry {
    pub fn solid_count(&self) -> usize {
        self.solid.iter().filter(|&&v| v).count()
    }

    pub fn fluid_count(&self) -> usize {
        self.solid.len() - self.solid_count()
    }

    pub fn index(&self, x: usize, y: usize, z: usize) -> usize {
        idx(self.dims, x, y, z)
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct ImpellerMarkerSet {
    pub kind: ImpellerKind,
    pub center: [f64; 3],
    pub omega: [f64; 3],
    pub markers: Vec<IbmMarker>,
    pub wall_velocity: Vec<[f64; 3]>,
}

impl ImpellerMarkerSet {
    pub fn rotating_body(&self) -> RotatingBody {
        RotatingBody::from_markers(self.center, self.omega, self.markers.clone())
    }
}

pub fn build_stirred_tank_geometry(
    grid: GridSpec,
    tank: TankSpec,
    baffles: &[BaffleTemplate],
    impellers: &[ImpellerTemplate],
    spargers: &[SpargerTemplate],
) -> GeometryResult<StirredTankGeometry> {
    grid.validate()?;
    validate_tank(tank)?;
    if tank.bottom == TankBottom::Dished {
        return Err(GeometryError::not_implemented(
            "dished-bottom stirred tanks are not implemented for M0",
        ));
    }
    validate_resolution(grid, tank, impellers)?;
    validate_impellers_inside_liquid(tank, impellers)?;
    validate_baffles_inside_tank(tank, baffles)?;
    validate_spargers(tank, grid, spargers)?;

    let n = grid.cell_count();
    let mut out = StirredTankGeometry {
        dims: grid.dims,
        dx_m: grid.dx_m,
        solid: vec![false; n],
        wall_velocity: vec![[0.0; 3]; n],
        baffle_mask: vec![false; n],
        sparger_mask: vec![false; n],
        sparger_orifice_centers: Vec::new(),
        impellers: Vec::new(),
    };
    apply_cylindrical_tank_mask(grid, tank, &mut out);
    for baffle in baffles {
        apply_baffle_mask(grid, tank, *baffle, &mut out)?;
    }
    for sparger in spargers {
        apply_sparger_mask(grid, tank, sparger, &mut out)?;
    }
    for impeller in impellers {
        let markers = generate_impeller_marker_set(grid, tank, impeller)?;
        stamp_impeller_wall_velocity(&markers, &mut out);
        out.impellers.push(markers);
    }
    Ok(out)
}

pub fn generate_impeller_marker_set(
    grid: GridSpec,
    tank: TankSpec,
    impeller: &ImpellerTemplate,
) -> GeometryResult<ImpellerMarkerSet> {
    grid.validate()?;
    validate_tank(tank)?;
    validate_impellers_inside_liquid(tank, std::slice::from_ref(impeller))?;
    let center_xy = tank_center_xy(grid);
    let omega_z = match impeller {
        ImpellerTemplate::Parametric {
            rotational_speed_rpm,
            ..
        }
        | ImpellerTemplate::CustomMarkerSet {
            rotational_speed_rpm,
            ..
        } => rpm_to_rad_s(*rotational_speed_rpm),
    };
    match impeller {
        ImpellerTemplate::Parametric {
            kind: ImpellerKind::Rushton,
            diameter_m,
            clearance_from_bottom_m,
            blade_count,
            ..
        } => build_radial_impeller(
            grid,
            tank,
            ImpellerKind::Rushton,
            center_xy,
            *diameter_m,
            *clearance_from_bottom_m,
            *blade_count,
            0.0,
            omega_z,
            RUSHTON_BLADE_THICKNESS_FRACTION,
        ),
        ImpellerTemplate::Parametric {
            kind: ImpellerKind::PitchedBlade,
            diameter_m,
            clearance_from_bottom_m,
            blade_count,
            ..
        } => build_radial_impeller(
            grid,
            tank,
            ImpellerKind::PitchedBlade,
            center_xy,
            *diameter_m,
            *clearance_from_bottom_m,
            *blade_count,
            PITCHED_BLADE_ANGLE_DEG,
            omega_z,
            PITCHED_BLADE_THICKNESS_FRACTION,
        ),
        ImpellerTemplate::Parametric {
            kind: ImpellerKind::Marine,
            ..
        } => Err(GeometryError::not_implemented(
            "marine impeller marker generation is not implemented",
        )),
        ImpellerTemplate::Parametric {
            kind: ImpellerKind::CustomMarkerSet,
            ..
        } => Err(GeometryError::not_implemented(
            "custom impeller marker generation requires ImpellerSpec::CustomMarkerSet markers",
        )),
        ImpellerTemplate::CustomMarkerSet { markers_m, .. } => {
            validate_custom_markers(tank, markers_m)?;
            Err(GeometryError::not_implemented(
                "custom impeller marker-set generation is not implemented",
            ))
        }
    }
}

pub fn impeller_target_velocity(center: [f64; 3], omega: [f64; 3], p: [f64; 3]) -> [f64; 3] {
    let r = [p[0] - center[0], p[1] - center[1], p[2] - center[2]];
    [
        omega[1] * r[2] - omega[2] * r[1],
        omega[2] * r[0] - omega[0] * r[2],
        omega[0] * r[1] - omega[1] * r[0],
    ]
}

fn validate_tank(tank: TankSpec) -> GeometryResult<()> {
    if !(tank.vessel_diameter_m.is_finite() && tank.vessel_diameter_m > 0.0) {
        return Err(GeometryError::out_of_validity_range(
            "vessel diameter must be finite and positive",
            "reactor.vessel_diameter_m must be > 0",
        ));
    }
    if !(tank.liquid_height_m.is_finite() && tank.liquid_height_m > 0.0) {
        return Err(GeometryError::out_of_validity_range(
            "liquid height must be finite and positive",
            "reactor.liquid_height_m must be > 0",
        ));
    }
    Ok(())
}

fn validate_resolution(
    grid: GridSpec,
    tank: TankSpec,
    impellers: &[ImpellerTemplate],
) -> GeometryResult<()> {
    let mut min_length = tank.vessel_diameter_m;
    for impeller in impellers {
        if let ImpellerTemplate::Parametric { diameter_m, .. } = impeller {
            min_length = min_length.min(*diameter_m);
        }
    }
    let cells = min_length / grid.dx_m;
    if cells < STIRRED_TANK_MIN_CELLS {
        return Err(GeometryError::out_of_validity_range(
            "stirred-tank geometry is under-resolved",
            format!(
                "min(impeller diameter, vessel diameter) / dx must be >= {STIRRED_TANK_MIN_CELLS}, got {cells:.3}"
            ),
        ));
    }
    Ok(())
}

fn validate_impellers_inside_liquid(
    tank: TankSpec,
    impellers: &[ImpellerTemplate],
) -> GeometryResult<()> {
    for impeller in impellers {
        match impeller {
            ImpellerTemplate::Parametric {
                diameter_m,
                clearance_from_bottom_m,
                rotational_speed_rpm,
                blade_count,
                ..
            } => {
                require_positive("impeller.diameter_m", *diameter_m)?;
                require_positive("impeller.rotational_speed_rpm", *rotational_speed_rpm)?;
                if *blade_count == 0 {
                    return Err(GeometryError::out_of_validity_range(
                        "impeller blade_count must be positive",
                        "impeller.blade_count must be > 0",
                    ));
                }
                let radius = 0.5 * diameter_m;
                if !(*clearance_from_bottom_m > 0.0
                    && *clearance_from_bottom_m + radius < tank.liquid_height_m)
                {
                    return Err(GeometryError::out_of_validity_range(
                        "impeller must be above the tank bottom and below the liquid surface",
                        "impeller clearance and radius must fit inside the liquid volume",
                    ));
                }
                if *diameter_m >= tank.vessel_diameter_m {
                    return Err(GeometryError::out_of_validity_range(
                        "impeller diameter must fit inside the vessel",
                        "impeller.diameter_m must be less than reactor.vessel_diameter_m",
                    ));
                }
            }
            ImpellerTemplate::CustomMarkerSet {
                markers_m,
                rotational_speed_rpm,
            } => {
                require_positive("impeller.rotational_speed_rpm", *rotational_speed_rpm)?;
                validate_custom_markers(tank, markers_m)?;
            }
        }
    }
    Ok(())
}

fn validate_custom_markers(tank: TankSpec, markers_m: &[[f64; 3]]) -> GeometryResult<()> {
    if markers_m.is_empty() {
        return Err(GeometryError::out_of_validity_range(
            "custom impeller marker set must not be empty",
            "impeller.markers must contain at least one marker",
        ));
    }
    for marker in markers_m {
        if !marker.iter().all(|v| v.is_finite()) {
            return Err(GeometryError::out_of_validity_range(
                "custom impeller markers must be finite",
                "impeller.markers entries must be finite",
            ));
        }
        let r = (marker[0] * marker[0] + marker[1] * marker[1]).sqrt();
        if !(r < 0.5 * tank.vessel_diameter_m
            && marker[2] > 0.0
            && marker[2] < tank.liquid_height_m)
        {
            return Err(GeometryError::out_of_validity_range(
                "custom impeller marker must lie inside the liquid volume",
                "custom marker radius and z must fit inside the tank liquid",
            ));
        }
    }
    Ok(())
}

fn validate_baffles_inside_tank(tank: TankSpec, baffles: &[BaffleTemplate]) -> GeometryResult<()> {
    for baffle in baffles {
        if baffle.count != 0 && baffle.count != 4 {
            return Err(GeometryError::not_implemented(
                "custom baffle placement is not implemented for M0; use count = 4",
            ));
        }
        if baffle.count == 0 {
            continue;
        }
        require_positive("baffle.width_m", baffle.width_m)?;
        require_positive("baffle.thickness_m", baffle.thickness_m)?;
        if !baffle.wall_attached {
            return Err(GeometryError::not_implemented(
                "detached baffle placement is not implemented for M0",
            ));
        }
        if baffle.width_m + baffle.thickness_m >= 0.5 * tank.vessel_diameter_m {
            return Err(GeometryError::out_of_validity_range(
                "baffles must fit inside the tank radius",
                "baffle width plus thickness must be less than the tank radius",
            ));
        }
    }
    Ok(())
}

fn validate_spargers(
    tank: TankSpec,
    grid: GridSpec,
    spargers: &[SpargerTemplate],
) -> GeometryResult<()> {
    for sparger in spargers {
        let (center_z_m, gas_flow, inlet_phase_gas, orifice_count, orifice_diameter) =
            sparger_validation_values(sparger);
        if !inlet_phase_gas {
            return Err(GeometryError::out_of_validity_range(
                "sparger inlet phase must be gas",
                "inlet_phase must be gas",
            ));
        }
        if !(center_z_m.is_finite() && center_z_m > 0.0 && center_z_m < tank.liquid_height_m) {
            return Err(GeometryError::out_of_validity_range(
                "sparger centre must be below the liquid surface and above the tank bottom",
                "sparger center_z_m must lie inside the liquid volume",
            ));
        }
        if gas_flow.map_or(true, |q| !(q.is_finite() && q > 0.0)) {
            return Err(GeometryError::out_of_validity_range(
                "sparger gas volumetric flow must be positive",
                "gas_volumetric_flow_m3_per_s must be finite and > 0",
            ));
        }
        if orifice_count == 0 {
            return Err(GeometryError::out_of_validity_range(
                "sparger orifice count must be positive",
                "orifice_count must be > 0",
            ));
        }
        require_positive("sparger.orifice_diameter_m", orifice_diameter)?;
        let cells = orifice_diameter / grid.dx_m;
        if cells < SPARGER_ORIFICE_MIN_CELLS {
            return Err(GeometryError::out_of_validity_range(
                "orifice under-resolved for resolved injection (BCFD-046)",
                format!(
                    "orifice_diameter_m / dx must be >= {SPARGER_ORIFICE_MIN_CELLS}, got {cells:.3}"
                ),
            ));
        }
        validate_sparger_inside_tank(tank, sparger)?;
    }
    Ok(())
}

fn sparger_validation_values(sparger: &SpargerTemplate) -> (f64, Option<f64>, bool, u32, f64) {
    match sparger {
        SpargerTemplate::Ring {
            center_z_m,
            orifice_count,
            orifice_diameter_m,
            gas_volumetric_flow_m3_per_s,
            inlet_phase_gas,
            ..
        } => (
            *center_z_m,
            *gas_volumetric_flow_m3_per_s,
            *inlet_phase_gas,
            *orifice_count,
            *orifice_diameter_m,
        ),
        SpargerTemplate::Pipe {
            center_z_m,
            orifice_count,
            orifice_diameter_m,
            gas_volumetric_flow_m3_per_s,
            inlet_phase_gas,
            ..
        } => (
            *center_z_m,
            *gas_volumetric_flow_m3_per_s,
            *inlet_phase_gas,
            *orifice_count,
            *orifice_diameter_m,
        ),
        SpargerTemplate::PointOrifices {
            center_z_m,
            positions_m,
            orifice_diameter_m,
            gas_volumetric_flow_m3_per_s,
            inlet_phase_gas,
        } => (
            *center_z_m,
            *gas_volumetric_flow_m3_per_s,
            *inlet_phase_gas,
            positions_m.len() as u32,
            *orifice_diameter_m,
        ),
    }
}

fn validate_sparger_inside_tank(tank: TankSpec, sparger: &SpargerTemplate) -> GeometryResult<()> {
    let radius = 0.5 * tank.vessel_diameter_m;
    match sparger {
        SpargerTemplate::Ring { outer_radius_m, .. } => {
            require_positive("sparger.outer_radius_m", *outer_radius_m)?;
            if *outer_radius_m >= radius {
                return Err(GeometryError::out_of_validity_range(
                    "ring sparger must be inside the tank",
                    "outer_radius_m must be less than the tank radius",
                ));
            }
        }
        SpargerTemplate::Pipe {
            length_m,
            diameter_m,
            ..
        } => {
            require_positive("sparger.length_m", *length_m)?;
            require_positive("sparger.diameter_m", *diameter_m)?;
            if 0.5 * *length_m + 0.5 * *diameter_m >= radius {
                return Err(GeometryError::out_of_validity_range(
                    "pipe sparger must be inside the tank",
                    "length_m/2 plus diameter_m/2 must be less than the tank radius",
                ));
            }
        }
        SpargerTemplate::PointOrifices { positions_m, .. } => {
            for p in positions_m {
                if !p.iter().all(|v| v.is_finite()) {
                    return Err(GeometryError::out_of_validity_range(
                        "point sparger positions must be finite",
                        "positions entries must be finite",
                    ));
                }
                let r = (p[0] * p[0] + p[1] * p[1]).sqrt();
                if !(r < radius && p[2] > 0.0 && p[2] < tank.liquid_height_m) {
                    return Err(GeometryError::out_of_validity_range(
                        "point sparger positions must lie inside the tank",
                        "each point-orifice position must be inside the liquid volume",
                    ));
                }
            }
        }
    }
    Ok(())
}

fn require_positive(name: &str, value: f64) -> GeometryResult<()> {
    if value.is_finite() && value > 0.0 {
        Ok(())
    } else {
        Err(GeometryError::out_of_validity_range(
            format!("{name} must be finite and positive"),
            format!("{name} must be finite and > 0"),
        ))
    }
}

fn apply_cylindrical_tank_mask(grid: GridSpec, tank: TankSpec, out: &mut StirredTankGeometry) {
    let radius = 0.5 * tank.vessel_diameter_m;
    let [cx, cy] = tank_center_xy(grid);
    for z in 0..grid.dims[2] {
        let z_m = cell_z_m(grid, z);
        for y in 0..grid.dims[1] {
            for x in 0..grid.dims[0] {
                let [px, py, _] = cell_center_m(grid, x, y, z);
                let r = ((px - cx) * (px - cx) + (py - cy) * (py - cy)).sqrt();
                let is_rim = x == 0
                    || y == 0
                    || z == 0
                    || x + 1 == grid.dims[0]
                    || y + 1 == grid.dims[1]
                    || z + 1 == grid.dims[2];
                let outside_liquid = r >= radius || z_m >= tank.liquid_height_m;
                if is_rim || outside_liquid {
                    out.solid[idx(grid.dims, x, y, z)] = true;
                }
            }
        }
    }
}

fn apply_baffle_mask(
    grid: GridSpec,
    tank: TankSpec,
    baffle: BaffleTemplate,
    out: &mut StirredTankGeometry,
) -> GeometryResult<()> {
    if baffle.count == 0 {
        return Ok(());
    }
    if baffle.count != 4 {
        return Err(GeometryError::not_implemented(
            "custom baffle placement is not implemented for M0; use count = 4",
        ));
    }
    let radius = 0.5 * tank.vessel_diameter_m;
    let [cx, cy] = tank_center_xy(grid);
    let radial_inner = radius - baffle.width_m;
    for z in 1..grid.dims[2].saturating_sub(1) {
        let z_m = cell_z_m(grid, z);
        if z_m >= tank.liquid_height_m {
            continue;
        }
        for y in 1..grid.dims[1].saturating_sub(1) {
            for x in 1..grid.dims[0].saturating_sub(1) {
                let [px, py, _] = cell_center_m(grid, x, y, z);
                let rx = px - cx;
                let ry = py - cy;
                let r = (rx * rx + ry * ry).sqrt();
                if r >= radius || r < radial_inner {
                    continue;
                }
                let near_x_wall = ry.abs() <= 0.5 * baffle.thickness_m;
                let near_y_wall = rx.abs() <= 0.5 * baffle.thickness_m;
                if near_x_wall || near_y_wall {
                    let i = idx(grid.dims, x, y, z);
                    out.solid[i] = true;
                    out.baffle_mask[i] = true;
                }
            }
        }
    }
    Ok(())
}

fn apply_sparger_mask(
    grid: GridSpec,
    tank: TankSpec,
    sparger: &SpargerTemplate,
    out: &mut StirredTankGeometry,
) -> GeometryResult<()> {
    match sparger {
        SpargerTemplate::Ring {
            center_z_m,
            outer_radius_m,
            orifice_count,
            orifice_diameter_m,
            ..
        } => {
            for k in 0..*orifice_count {
                let th = TAU * k as f64 / *orifice_count as f64;
                let center = [
                    *outer_radius_m * th.cos(),
                    *outer_radius_m * th.sin(),
                    *center_z_m,
                ];
                stamp_orifice(grid, tank, center, *orifice_diameter_m, out);
            }
        }
        SpargerTemplate::Pipe {
            center_z_m,
            length_m,
            axis,
            orifice_count,
            orifice_diameter_m,
            ..
        } => {
            let denom = (*orifice_count).max(2) - 1;
            for k in 0..*orifice_count {
                let s = if *orifice_count == 1 {
                    0.0
                } else {
                    -0.5 * length_m + *length_m * k as f64 / denom as f64
                };
                let center = match axis {
                    PipeAxis::X => [s, 0.0, *center_z_m],
                    PipeAxis::Y => [0.0, s, *center_z_m],
                };
                stamp_orifice(grid, tank, center, *orifice_diameter_m, out);
            }
        }
        SpargerTemplate::PointOrifices {
            positions_m,
            orifice_diameter_m,
            ..
        } => {
            for &center in positions_m {
                stamp_orifice(grid, tank, center, *orifice_diameter_m, out);
            }
        }
    }
    Ok(())
}

fn stamp_orifice(
    grid: GridSpec,
    tank: TankSpec,
    center_rel_m: [f64; 3],
    diameter_m: f64,
    out: &mut StirredTankGeometry,
) {
    let [cx, cy] = tank_center_xy(grid);
    let center = [cx + center_rel_m[0], cy + center_rel_m[1], center_rel_m[2]];
    let radius = 0.5 * diameter_m;
    out.sparger_orifice_centers.push(center);
    for z in 1..grid.dims[2].saturating_sub(1) {
        let z_m = cell_z_m(grid, z);
        if z_m >= tank.liquid_height_m {
            continue;
        }
        for y in 1..grid.dims[1].saturating_sub(1) {
            for x in 1..grid.dims[0].saturating_sub(1) {
                let p = cell_center_m(grid, x, y, z);
                let d = ((p[0] - center[0]).powi(2)
                    + (p[1] - center[1]).powi(2)
                    + (p[2] - center[2]).powi(2))
                .sqrt();
                if d <= radius {
                    out.sparger_mask[idx(grid.dims, x, y, z)] = true;
                }
            }
        }
    }
}

fn build_radial_impeller(
    grid: GridSpec,
    tank: TankSpec,
    kind: ImpellerKind,
    center_xy: [f64; 2],
    diameter_m: f64,
    clearance_from_bottom_m: f64,
    blade_count: u32,
    blade_angle_deg: f64,
    omega_z: f64,
    thickness_fraction: f64,
) -> GeometryResult<ImpellerMarkerSet> {
    if !(0.0..=90.0).contains(&blade_angle_deg) {
        return Err(GeometryError::out_of_validity_range(
            "impeller blade angle must be between 0 and 90 degrees",
            "blade_angle_deg must be in [0, 90]",
        ));
    }
    let radius = 0.5 * diameter_m;
    let center = [center_xy[0], center_xy[1], clearance_from_bottom_m];
    let radial_inner = match kind {
        ImpellerKind::Rushton => RUSHTON_DISC_DIAMETER_FRACTION * radius,
        ImpellerKind::PitchedBlade => 0.35 * radius,
        _ => 0.0,
    };
    let dr = grid.dx_m.max(radius / 20.0);
    let dz = grid.dx_m;
    let tangential_step = grid.dx_m;
    let thickness = (thickness_fraction * diameter_m).max(grid.dx_m);
    let mut markers = Vec::new();
    let mut wall_velocity = Vec::new();
    for b in 0..blade_count {
        let th = TAU * b as f64 / blade_count as f64;
        let er = [th.cos(), th.sin(), 0.0];
        let et = [-th.sin(), th.cos(), 0.0];
        let pitch = blade_angle_deg.to_radians();
        let z_span = if blade_angle_deg == 0.0 {
            thickness
        } else {
            thickness + (radius - radial_inner) * pitch.sin() * 0.18
        };
        let nr = ((radius - radial_inner) / dr).ceil().max(1.0) as usize;
        let nt = (thickness / tangential_step).ceil().max(1.0) as usize;
        let nz = (z_span / dz).ceil().max(1.0) as usize;
        for ir in 0..=nr {
            let rr = radial_inner + (radius - radial_inner) * ir as f64 / nr as f64;
            for it in 0..=nt {
                let offset_t = -0.5 * thickness + thickness * it as f64 / nt as f64;
                for iz in 0..=nz {
                    let offset_z = -0.5 * z_span + z_span * iz as f64 / nz as f64;
                    let pitched_z = offset_z + (rr - radial_inner) * pitch.sin() * 0.18;
                    let p = [
                        center[0] + rr * er[0] + offset_t * et[0],
                        center[1] + rr * er[1] + offset_t * et[1],
                        center[2] + pitched_z,
                    ];
                    let rel_r =
                        ((p[0] - center_xy[0]).powi(2) + (p[1] - center_xy[1]).powi(2)).sqrt();
                    if rel_r >= 0.5 * tank.vessel_diameter_m
                        || p[2] <= 0.0
                        || p[2] >= tank.liquid_height_m
                    {
                        continue;
                    }
                    let u = impeller_target_velocity(center, [0.0, 0.0, omega_z], p);
                    markers.push(IbmMarker {
                        position: [p[0] / grid.dx_m, p[1] / grid.dx_m, p[2] / grid.dx_m],
                        weight: grid.dx_m * grid.dx_m,
                    });
                    wall_velocity.push(u);
                }
            }
        }
    }
    if markers.is_empty() {
        return Err(GeometryError::out_of_validity_range(
            "impeller marker generation produced no markers",
            "increase grid resolution or impeller diameter",
        ));
    }
    Ok(ImpellerMarkerSet {
        kind,
        center: [
            center[0] / grid.dx_m,
            center[1] / grid.dx_m,
            center[2] / grid.dx_m,
        ],
        omega: [0.0, 0.0, omega_z],
        markers,
        wall_velocity,
    })
}

fn stamp_impeller_wall_velocity(markers: &ImpellerMarkerSet, out: &mut StirredTankGeometry) {
    for (marker, velocity) in markers.markers.iter().zip(markers.wall_velocity.iter()) {
        let x = nearest_index(marker.position[0], out.dims[0]);
        let y = nearest_index(marker.position[1], out.dims[1]);
        let z = nearest_index(marker.position[2], out.dims[2]);
        let i = idx(out.dims, x, y, z);
        out.wall_velocity[i] = *velocity;
    }
}

fn nearest_index(v: f64, n: usize) -> usize {
    v.round().max(0.0).min((n - 1) as f64) as usize
}

fn tank_center_xy(grid: GridSpec) -> [f64; 2] {
    [
        0.5 * grid.dims[0] as f64 * grid.dx_m,
        0.5 * grid.dims[1] as f64 * grid.dx_m,
    ]
}

fn cell_center_m(grid: GridSpec, x: usize, y: usize, z: usize) -> [f64; 3] {
    [
        (x as f64 + 0.5) * grid.dx_m,
        (y as f64 + 0.5) * grid.dx_m,
        cell_z_m(grid, z),
    ]
}

fn cell_z_m(grid: GridSpec, z: usize) -> f64 {
    (z as f64 + 0.5) * grid.dx_m
}

fn idx(dims: [usize; 3], x: usize, y: usize, z: usize) -> usize {
    (z * dims[1] + y) * dims[0] + x
}

fn rpm_to_rad_s(rpm: f64) -> f64 {
    rpm * TAU / 60.0
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    fn grid(n: usize) -> GridSpec {
        GridSpec {
            dims: [n, n, n],
            dx_m: 1.0 / n as f64,
        }
    }

    fn tank() -> TankSpec {
        TankSpec {
            vessel_diameter_m: 1.0,
            liquid_height_m: 1.0,
            bottom: TankBottom::Flat,
        }
    }

    fn rushton(diameter_m: f64) -> ImpellerTemplate {
        ImpellerTemplate::Parametric {
            kind: ImpellerKind::Rushton,
            diameter_m,
            clearance_from_bottom_m: 0.35,
            rotational_speed_rpm: 120.0,
            blade_count: 6,
        }
    }

    #[test]
    fn tank_mask_volume_within_1pct_of_analytic_cylinder() {
        let g = grid(160);
        let geom = build_stirred_tank_geometry(g, tank(), &[], &[rushton(0.34)], &[]).unwrap();
        let fluid_volume = geom.fluid_count() as f64 * g.dx_m.powi(3);
        let analytic = PI * 0.5f64.powi(2) * ((g.dims[2] - 2) as f64 * g.dx_m);
        let rel = ((fluid_volume - analytic) / analytic).abs();
        assert!(
            rel < 0.01,
            "voxel cylinder volume rel error {rel:.5} must be < 1%; fluid={fluid_volume:.6} analytic={analytic:.6}"
        );
    }

    #[test]
    fn baffle_count_and_placement_matches_spec() {
        let g = grid(96);
        let baffle = BaffleTemplate {
            count: 4,
            width_m: 0.08,
            thickness_m: 0.025,
            wall_attached: true,
        };
        let geom =
            build_stirred_tank_geometry(g, tank(), &[baffle], &[rushton(0.34)], &[]).unwrap();
        let [cx, cy] = tank_center_xy(g);
        let mut quadrants = [0usize; 4];
        for z in 0..g.dims[2] {
            for y in 0..g.dims[1] {
                for x in 0..g.dims[0] {
                    let i = geom.index(x, y, z);
                    if !geom.baffle_mask[i] {
                        continue;
                    }
                    let p = cell_center_m(g, x, y, z);
                    let dx = p[0] - cx;
                    let dy = p[1] - cy;
                    if dy.abs() < dx.abs() {
                        quadrants[usize::from(dx > 0.0)] += 1;
                    } else {
                        quadrants[2 + usize::from(dy > 0.0)] += 1;
                    }
                }
            }
        }
        assert!(
            quadrants.iter().all(|&count| count > 0),
            "expected baffle cells on all four cardinal placements: {quadrants:?}"
        );
    }

    #[test]
    fn rejects_impeller_outside_liquid_volume() {
        let bad = ImpellerTemplate::Parametric {
            kind: ImpellerKind::Rushton,
            diameter_m: 0.4,
            clearance_from_bottom_m: 0.9,
            rotational_speed_rpm: 120.0,
            blade_count: 6,
        };
        let err = build_stirred_tank_geometry(grid(96), tank(), &[], &[bad], &[]).unwrap_err();
        assert!(matches!(
            err.reason,
            UnsupportedReason::OutOfValidityRange { .. }
        ));
    }

    #[test]
    fn rejects_dished_bottom() {
        let mut t = tank();
        t.bottom = TankBottom::Dished;
        let err = build_stirred_tank_geometry(grid(96), t, &[], &[rushton(0.34)], &[]).unwrap_err();
        assert_eq!(err.reason, UnsupportedReason::NotImplemented);
    }

    #[test]
    fn rejects_under_resolved_tank() {
        let err =
            build_stirred_tank_geometry(grid(48), tank(), &[], &[rushton(0.2)], &[]).unwrap_err();
        assert!(matches!(
            err.reason,
            UnsupportedReason::OutOfValidityRange { .. }
        ));
    }

    #[test]
    fn rushton_marker_count_scales_with_resolution() {
        let low = generate_impeller_marker_set(grid(64), tank(), &rushton(0.5)).unwrap();
        let high = generate_impeller_marker_set(grid(96), tank(), &rushton(0.5)).unwrap();
        assert!(
            high.markers.len() > low.markers.len(),
            "marker count should increase with resolution: low={} high={}",
            low.markers.len(),
            high.markers.len()
        );
    }

    #[test]
    fn pitched_blade_marker_count_scales_with_resolution() {
        let imp = || ImpellerTemplate::Parametric {
            kind: ImpellerKind::PitchedBlade,
            diameter_m: 0.5,
            clearance_from_bottom_m: 0.35,
            rotational_speed_rpm: 120.0,
            blade_count: 4,
        };
        let low = generate_impeller_marker_set(grid(64), tank(), &imp()).unwrap();
        let high = generate_impeller_marker_set(grid(96), tank(), &imp()).unwrap();
        assert!(high.markers.len() > low.markers.len());
    }

    #[test]
    fn impeller_marker_bounding_box_inside_vessel() {
        let markers = generate_impeller_marker_set(grid(96), tank(), &rushton(0.5)).unwrap();
        for marker in &markers.markers {
            let x = marker.position[0] * grid(96).dx_m - 0.5;
            let y = marker.position[1] * grid(96).dx_m - 0.5;
            let z = marker.position[2] * grid(96).dx_m;
            assert!((x * x + y * y).sqrt() < 0.5);
            assert!(z > 0.0 && z < 1.0);
        }
    }

    #[test]
    fn impeller_rotation_target_velocity_matches_omega_cross_r() {
        let center = [10.0, 10.0, 5.0];
        let omega = [0.0, 0.0, 2.0];
        let p = [13.0, 14.0, 5.0];
        assert_eq!(impeller_target_velocity(center, omega, p), [-8.0, 6.0, 0.0]);
    }

    #[test]
    fn rejects_invalid_blade_geometry() {
        let err = build_radial_impeller(
            grid(64),
            tank(),
            ImpellerKind::PitchedBlade,
            tank_center_xy(grid(64)),
            0.4,
            0.35,
            4,
            120.0,
            1.0,
            PITCHED_BLADE_THICKNESS_FRACTION,
        )
        .unwrap_err();
        assert!(matches!(
            err.reason,
            UnsupportedReason::OutOfValidityRange { .. }
        ));
    }

    #[test]
    fn marine_impeller_rejects_with_structured_error() {
        let marine = ImpellerTemplate::Parametric {
            kind: ImpellerKind::Marine,
            diameter_m: 0.35,
            clearance_from_bottom_m: 0.3,
            rotational_speed_rpm: 120.0,
            blade_count: 3,
        };
        let err = generate_impeller_marker_set(grid(96), tank(), &marine).unwrap_err();
        assert_eq!(err.reason, UnsupportedReason::NotImplemented);
    }

    #[test]
    fn ring_sparger_produces_expected_orifice_count() {
        let sparger = SpargerTemplate::Ring {
            center_z_m: 0.12,
            outer_radius_m: 0.25,
            orifice_count: 12,
            orifice_diameter_m: 0.04,
            gas_volumetric_flow_m3_per_s: Some(1.0e-5),
            inlet_phase_gas: true,
        };
        let geom =
            build_stirred_tank_geometry(grid(120), tank(), &[], &[rushton(0.34)], &[sparger])
                .unwrap();
        assert_eq!(geom.sparger_orifice_centers.len(), 12);
        assert!(geom.sparger_mask.iter().any(|&v| v));
    }

    #[test]
    fn pipe_sparger_mask_inside_tank() {
        let sparger = SpargerTemplate::Pipe {
            center_z_m: 0.12,
            length_m: 0.45,
            diameter_m: 0.03,
            axis: PipeAxis::X,
            orifice_count: 5,
            orifice_diameter_m: 0.03,
            gas_volumetric_flow_m3_per_s: Some(1.0e-5),
            inlet_phase_gas: true,
        };
        let geom =
            build_stirred_tank_geometry(grid(120), tank(), &[], &[rushton(0.34)], &[sparger])
                .unwrap();
        for z in 0..geom.dims[2] {
            for y in 0..geom.dims[1] {
                for x in 0..geom.dims[0] {
                    let i = geom.index(x, y, z);
                    if geom.sparger_mask[i] {
                        let p = cell_center_m(grid(120), x, y, z);
                        let r = ((p[0] - 0.5).powi(2) + (p[1] - 0.5).powi(2)).sqrt();
                        assert!(r < 0.5);
                    }
                }
            }
        }
    }

    #[test]
    fn rejects_raw_phi_boundary_fields() {
        let sparger = SpargerTemplate::Ring {
            center_z_m: 0.12,
            outer_radius_m: 0.25,
            orifice_count: 12,
            orifice_diameter_m: 0.04,
            gas_volumetric_flow_m3_per_s: Some(1.0e-5),
            inlet_phase_gas: false,
        };
        let err = build_stirred_tank_geometry(grid(120), tank(), &[], &[rushton(0.34)], &[sparger])
            .unwrap_err();
        assert!(matches!(
            err.reason,
            UnsupportedReason::OutOfValidityRange { .. }
        ));
    }

    #[test]
    fn rejects_sparger_below_no_liquid_surface() {
        let sparger = SpargerTemplate::Ring {
            center_z_m: 1.2,
            outer_radius_m: 0.25,
            orifice_count: 12,
            orifice_diameter_m: 0.04,
            gas_volumetric_flow_m3_per_s: Some(1.0e-5),
            inlet_phase_gas: true,
        };
        let err = build_stirred_tank_geometry(grid(120), tank(), &[], &[rushton(0.34)], &[sparger])
            .unwrap_err();
        assert!(matches!(
            err.reason,
            UnsupportedReason::OutOfValidityRange { .. }
        ));
    }

    #[test]
    fn rejects_under_resolved_orifice() {
        let sparger = SpargerTemplate::Ring {
            center_z_m: 0.12,
            outer_radius_m: 0.25,
            orifice_count: 12,
            orifice_diameter_m: 0.01,
            gas_volumetric_flow_m3_per_s: Some(1.0e-5),
            inlet_phase_gas: true,
        };
        let err = build_stirred_tank_geometry(grid(120), tank(), &[], &[rushton(0.34)], &[sparger])
            .unwrap_err();
        assert!(matches!(
            err.reason,
            UnsupportedReason::OutOfValidityRange { .. }
        ));
    }
}
