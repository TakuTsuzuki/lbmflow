use crate::Scenario;
use lbm_core::compat::prelude::MAX_SPEED;
use serde::{Deserialize, Serialize};

pub const TAU_LOW_WARN_THRESHOLD: f64 = 0.55;
pub const TAU_HIGH_WARN_THRESHOLD: f64 = 2.0;
pub const LATTICE_SPEED_WARN_THRESHOLD: f64 = 0.15;
pub const GRID_RE_WARN_THRESHOLD: f64 = 15.0;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum UnitVerdict {
    #[serde(rename = "ok")]
    Ok,
    #[serde(rename = "warn")]
    Warn,
    #[serde(rename = "error")]
    Error,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum UnitConstructor {
    FromResolutionAndRelaxationTime,
    FromResolutionAndLatticeVelocity,
    FromRelaxationTimeAndLatticeVelocity,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(rename_all = "camelCase", deny_unknown_fields)]
pub struct FlowParams {
    pub constructor: UnitConstructor,
    pub characteristic_length: f64,
    pub characteristic_velocity: f64,
    pub kinematic_viscosity: f64,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub density: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub resolution: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub lattice_velocity: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub relaxation_time: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub end_time: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub end_step_count: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub gravity: Option<[f64; 2]>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub reference_pressure: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub re_physical: Option<f64>,
}

#[derive(Clone, Debug, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct UnitInputsEcho {
    pub characteristic_length: f64,
    pub characteristic_velocity: f64,
    pub kinematic_viscosity: f64,
    pub density: Option<f64>,
    pub resolution: Option<usize>,
    pub lattice_velocity: Option<f64>,
    pub relaxation_time: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub end_time: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub end_step_count: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub gravity: Option<[f64; 2]>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub reference_pressure: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub re_physical: Option<f64>,
}

#[derive(Clone, Debug, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct LatticeUnits {
    pub dx: f64,
    pub dt: f64,
    #[serde(rename = "nu_lattice")]
    pub nu_lattice: f64,
    pub tau: f64,
    pub omega: f64,
    #[serde(rename = "u_char_lattice")]
    pub u_char_lattice: f64,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub total_steps: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub g_lat: Option<[f64; 2]>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub p_ref_lat: Option<f64>,
}

#[derive(Clone, Debug, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct ConversionFactors {
    #[serde(rename = "length_m")]
    pub length_m: f64,
    #[serde(rename = "time_s")]
    pub time_s: f64,
    #[serde(rename = "velocity_m_s")]
    pub velocity_m_s: f64,
    #[serde(rename = "viscosity_m2_s")]
    pub viscosity_m2_s: f64,
    #[serde(rename = "density_kg_m3")]
    pub density_kg_m3: Option<f64>,
    #[serde(rename = "pressure_Pa")]
    pub pressure_pa: Option<f64>,
    #[serde(rename = "force_N")]
    pub force_n: Option<f64>,
    #[serde(rename = "acceleration_m_s2")]
    pub acceleration_m_s2: f64,
}

#[derive(Clone, Debug, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct DimensionlessNumbers {
    pub reynolds: f64,
    pub mach: f64,
    pub grid_reynolds: f64,
    pub knudsen: f64,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub re_declared: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub re_physical: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub matching_ratio: Option<f64>,
}

#[derive(Clone, Debug, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct UnitDiagnostic {
    pub id: String,
    pub severity: UnitVerdict,
    pub quantity: String,
    pub value: f64,
    pub threshold: f64,
    pub message: String,
    pub remedy: String,
}

#[derive(Clone, Copy, Debug, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct UnitSuggestion {
    pub resolution: usize,
    pub lattice_velocity: f64,
}

#[derive(Clone, Debug, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct UnitReport {
    pub constructor: UnitConstructor,
    pub inputs: UnitInputsEcho,
    pub lattice: LatticeUnits,
    pub conversion_factors: ConversionFactors,
    pub dimensionless: DimensionlessNumbers,
    pub verdict: UnitVerdict,
    pub diagnostics: Vec<UnitDiagnostic>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub suggestion: Option<UnitSuggestion>,
}

#[derive(Clone, Debug)]
pub struct ResolvedScenario {
    pub scenario: Scenario,
    pub report: UnitReport,
}

pub fn resolve(sc: &Scenario) -> Result<Option<ResolvedScenario>, String> {
    let Some(params) = &sc.units else {
        return Ok(None);
    };

    let report = report(params)?;
    if report.verdict == UnitVerdict::Error {
        let ids = report
            .diagnostics
            .iter()
            .filter(|d| d.severity == UnitVerdict::Error)
            .map(|d| d.id.as_str())
            .collect::<Vec<_>>()
            .join(", ");
        return Err(format!("unit conversion failed: {ids}"));
    }

    let mut scenario = sc.clone();
    scenario.physics.nu = report.lattice.nu_lattice;
    if let Some(g) = report.lattice.g_lat {
        scenario.physics.force = g;
    }
    if let Some(total_steps) = report.lattice.total_steps {
        scenario.run.steps = total_steps;
    }
    scenario.units = Some(params.clone());
    Ok(Some(ResolvedScenario { scenario, report }))
}

pub fn report(params: &FlowParams) -> Result<UnitReport, String> {
    let re =
        params.characteristic_velocity * params.characteristic_length / params.kinematic_viscosity;
    if !(params.characteristic_length > 0.0
        && params.characteristic_velocity > 0.0
        && params.kinematic_viscosity > 0.0
        && re.is_finite())
    {
        return Err("unit inputs must be finite and positive for L, U and nu".to_string());
    }

    let mut diagnostics = Vec::new();
    if params.density.is_none() {
        diagnostics.push(diag(
            "DENSITY_MISSING",
            UnitVerdict::Error,
            "density",
            0.0,
            0.0,
            "density is required for pressure and force conversion.",
            "set density explicitly, e.g. 998.2 for water or 1.204 for air at 20 C.",
        ));
    } else if !(params.density.unwrap() > 0.0) {
        diagnostics.push(diag(
            "DENSITY_MISSING",
            UnitVerdict::Error,
            "density",
            params.density.unwrap(),
            0.0,
            "density must be finite and positive.",
            "set density explicitly, e.g. 998.2 for water or 1.204 for air at 20 C.",
        ));
    }
    validate_constructor_knobs(params)?;

    let (resolution, u_lat, tau, rounded_from) = match params.constructor {
        UnitConstructor::FromResolutionAndLatticeVelocity => {
            let n = require_resolution(params)?;
            let u = require_lattice_velocity(params)?;
            let nu_lat = u * n as f64 / re;
            (n, u, 3.0 * nu_lat + 0.5, None)
        }
        UnitConstructor::FromResolutionAndRelaxationTime => {
            let n = require_resolution(params)?;
            let tau = require_relaxation_time(params)?;
            let nu_lat = (tau - 0.5) / 3.0;
            (n, nu_lat * re / n as f64, tau, None)
        }
        UnitConstructor::FromRelaxationTimeAndLatticeVelocity => {
            let tau = require_relaxation_time(params)?;
            let u = require_lattice_velocity(params)?;
            let nu_lat = (tau - 0.5) / 3.0;
            let n_raw = nu_lat * re / u;
            let n = n_raw.round().max(1.0) as usize;
            let rounded = if (n_raw - n as f64).abs() > 1e-12 {
                Some(n_raw)
            } else {
                None
            };
            let actual_nu = u * n as f64 / re;
            (n, u, 3.0 * actual_nu + 0.5, rounded)
        }
    };
    if !(u_lat >= 0.0 && tau.is_finite()) {
        return Err(
            "unit constructor produced negative lattice velocity or non-finite tau".to_string(),
        );
    }

    let dx = params.characteristic_length / resolution as f64;
    let dt = u_lat * dx / params.characteristic_velocity;
    let nu_lat = (tau - 0.5) / 3.0;
    let omega = 1.0 / tau;
    let c_velocity = dx / dt;
    let c_accel = dx / (dt * dt);
    let density = params.density;
    let pressure = density.map(|rho| rho * c_velocity * c_velocity);
    let force = density.map(|rho| rho * dx.powi(4) / (dt * dt));
    let total_steps = params
        .end_step_count
        .or_else(|| params.end_time.map(|t| (t / dt).round().max(0.0) as usize));
    let g_lat = params.gravity.map(|g| [g[0] / c_accel, g[1] / c_accel]);
    let p_ref_lat = match (params.reference_pressure, pressure) {
        (Some(p), Some(c)) => Some(p / c),
        _ => None,
    };

    let mach = 3.0f64.sqrt() * u_lat;
    let grid_re = re / resolution as f64;
    let knudsen = mach / re;

    if tau <= 0.5 {
        diagnostics.push(diag(
            "TAU_UNSTABLE",
            UnitVerdict::Error,
            "tau",
            tau,
            0.5,
            "tau must exceed 0.5 for positive lattice viscosity.",
            "raise resolution N or latticeVelocity; or lower Re.",
        ));
    } else if tau < TAU_LOW_WARN_THRESHOLD {
        diagnostics.push(diag(
            "TAU_LOW",
            UnitVerdict::Warn,
            "tau",
            tau,
            TAU_LOW_WARN_THRESHOLD,
            "tau below 0.55: BGK over-relaxation risk near walls.",
            "raise resolution N or latticeVelocity; or use TRT collision.",
        ));
    }
    if tau > TAU_HIGH_WARN_THRESHOLD {
        diagnostics.push(diag(
            "TAU_HIGH",
            UnitVerdict::Warn,
            "tau",
            tau,
            TAU_HIGH_WARN_THRESHOLD,
            "tau above 2.0: over-diffusive and likely wasting cells.",
            "lower resolution N or latticeVelocity to cut step count.",
        ));
    }
    if u_lat > MAX_SPEED {
        diagnostics.push(diag(
            "MACH_HARD",
            UnitVerdict::Error,
            "mach",
            mach,
            3.0f64.sqrt() * MAX_SPEED,
            "lattice velocity above 0.3: compressibility error dominates.",
            "lower latticeVelocity; compensate tau by raising resolution N.",
        ));
    } else if u_lat > LATTICE_SPEED_WARN_THRESHOLD {
        diagnostics.push(diag(
            "MACH_HIGH",
            UnitVerdict::Warn,
            "mach",
            mach,
            3.0f64.sqrt() * LATTICE_SPEED_WARN_THRESHOLD,
            "lattice velocity above 0.15: non-negligible compressibility error.",
            "lower latticeVelocity toward <= 0.15; raise N to hold tau.",
        ));
    }
    if grid_re > GRID_RE_WARN_THRESHOLD {
        diagnostics.push(diag(
            "GRID_RE_HIGH",
            UnitVerdict::Warn,
            "grid_reynolds",
            grid_re,
            GRID_RE_WARN_THRESHOLD,
            "grid Reynolds > 15: cell under-resolved.",
            &format!(
                "raise N so N >= Re/15 (>= {} here).",
                (re / GRID_RE_WARN_THRESHOLD).ceil() as usize
            ),
        ));
    }
    if let Some(n_raw) = rounded_from {
        diagnostics.push(diag(
            "RESOLUTION_ROUNDING",
            UnitVerdict::Warn,
            "resolution",
            resolution as f64,
            n_raw,
            "derived resolution was non-integer and was rounded.",
            "switch to an N-fixed constructor, or accept the reported drift.",
        ));
    }
    if let Some(re_physical) = params.re_physical {
        if (re_physical - re).abs() > 1e-12 * re_physical.abs().max(1.0) {
            diagnostics.push(diag(
                "EFFECTIVE_VISCOSITY_REGIME",
                UnitVerdict::Ok,
                "reynolds",
                re,
                re_physical,
                "declared Reynolds number is intentionally matched below the physical Reynolds number.",
                "interpret the simulation as an effective-viscosity regime.",
            ));
        }
    }

    let verdict = if diagnostics.iter().any(|d| d.severity == UnitVerdict::Error) {
        UnitVerdict::Error
    } else if diagnostics.iter().any(|d| d.severity == UnitVerdict::Warn) {
        UnitVerdict::Warn
    } else {
        UnitVerdict::Ok
    };
    let suggestion = if verdict == UnitVerdict::Warn {
        let mut n = (re / GRID_RE_WARN_THRESHOLD).ceil().max(1.0) as usize;
        let mut u = LATTICE_SPEED_WARN_THRESHOLD
            .min(((TAU_LOW_WARN_THRESHOLD - 0.5) / 3.0) * re / n as f64);
        while 3.0 * u * n as f64 / re + 0.5 + 1e-15 < TAU_LOW_WARN_THRESHOLD {
            n += 1;
            u = LATTICE_SPEED_WARN_THRESHOLD
                .min(((TAU_LOW_WARN_THRESHOLD - 0.5) / 3.0) * re / n as f64);
        }
        Some(UnitSuggestion {
            resolution: n,
            lattice_velocity: u,
        })
    } else {
        None
    };

    Ok(UnitReport {
        constructor: params.constructor,
        inputs: UnitInputsEcho {
            characteristic_length: params.characteristic_length,
            characteristic_velocity: params.characteristic_velocity,
            kinematic_viscosity: params.kinematic_viscosity,
            density,
            resolution: params.resolution,
            lattice_velocity: params.lattice_velocity,
            relaxation_time: params.relaxation_time,
            end_time: params.end_time,
            end_step_count: params.end_step_count,
            gravity: params.gravity,
            reference_pressure: params.reference_pressure,
            re_physical: params.re_physical,
        },
        lattice: LatticeUnits {
            dx,
            dt,
            nu_lattice: nu_lat,
            tau,
            omega,
            u_char_lattice: u_lat,
            total_steps,
            g_lat,
            p_ref_lat,
        },
        conversion_factors: ConversionFactors {
            length_m: dx,
            time_s: dt,
            velocity_m_s: c_velocity,
            viscosity_m2_s: dx * dx / dt,
            density_kg_m3: density,
            pressure_pa: pressure,
            force_n: force,
            acceleration_m_s2: c_accel,
        },
        dimensionless: DimensionlessNumbers {
            reynolds: re,
            mach,
            grid_reynolds: grid_re,
            knudsen,
            re_declared: Some(re),
            re_physical: params.re_physical,
            matching_ratio: params.re_physical.map(|rp| re / rp),
        },
        verdict,
        diagnostics,
        suggestion,
    })
}

fn validate_constructor_knobs(params: &FlowParams) -> Result<(), String> {
    let has_n = params.resolution.is_some();
    let has_u = params.lattice_velocity.is_some();
    let has_tau = params.relaxation_time.is_some();
    let ok = match params.constructor {
        UnitConstructor::FromResolutionAndLatticeVelocity => has_n && has_u && !has_tau,
        UnitConstructor::FromResolutionAndRelaxationTime => has_n && has_tau && !has_u,
        UnitConstructor::FromRelaxationTimeAndLatticeVelocity => has_tau && has_u && !has_n,
    };
    if ok {
        Ok(())
    } else {
        Err(format!(
            "{:?} requires exactly two matching numerical knobs: resolution={}, latticeVelocity={}, relaxationTime={}",
            params.constructor, has_n, has_u, has_tau
        ))
    }
}

fn require_resolution(params: &FlowParams) -> Result<usize, String> {
    let n = params
        .resolution
        .ok_or_else(|| "constructor requires resolution".to_string())?;
    if n == 0 {
        return Err("resolution must be positive".to_string());
    }
    Ok(n)
}

fn require_lattice_velocity(params: &FlowParams) -> Result<f64, String> {
    let u = params
        .lattice_velocity
        .ok_or_else(|| "constructor requires latticeVelocity".to_string())?;
    if !(u > 0.0) {
        return Err("latticeVelocity must be finite and positive".to_string());
    }
    Ok(u)
}

fn require_relaxation_time(params: &FlowParams) -> Result<f64, String> {
    let tau = params
        .relaxation_time
        .ok_or_else(|| "constructor requires relaxationTime".to_string())?;
    if !tau.is_finite() {
        return Err("relaxationTime must be finite".to_string());
    }
    Ok(tau)
}

fn diag(
    id: &str,
    severity: UnitVerdict,
    quantity: &str,
    value: f64,
    threshold: f64,
    message: &str,
    remedy: &str,
) -> UnitDiagnostic {
    UnitDiagnostic {
        id: id.to_string(),
        severity,
        quantity: quantity.to_string(),
        value,
        threshold,
        message: message.to_string(),
        remedy: remedy.to_string(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        unit_report_for, CollisionSpec, EdgeSpec, EdgesSpec, Grid, Physics, Precision, RunSpec,
    };

    fn base_params(constructor: UnitConstructor) -> FlowParams {
        FlowParams {
            constructor,
            characteristic_length: 0.1,
            characteristic_velocity: 0.2,
            kinematic_viscosity: 2.0e-4,
            density: Some(998.2),
            resolution: None,
            lattice_velocity: None,
            relaxation_time: None,
            end_time: None,
            end_step_count: None,
            gravity: None,
            reference_pressure: Some(0.0),
            re_physical: None,
        }
    }

    fn report_ids(r: &UnitReport) -> Vec<&str> {
        r.diagnostics.iter().map(|d| d.id.as_str()).collect()
    }

    fn assert_close(a: f64, b: f64) {
        let scale = a.abs().max(b.abs()).max(1.0);
        assert!((a - b).abs() <= 1e-12 * scale, "{a:e} != {b:e}");
    }

    #[test]
    fn conversion_factors_round_trip_to_tight_tolerance() {
        let mut p = base_params(UnitConstructor::FromResolutionAndLatticeVelocity);
        p.resolution = Some(160);
        p.lattice_velocity = Some(0.08);
        let r = report(&p).unwrap();
        let c = &r.conversion_factors;
        for (q, factor) in [
            (0.037, c.velocity_m_s),
            (1.7e-6, c.viscosity_m2_s),
            (23.0, c.pressure_pa.unwrap()),
            (4.2e-5, c.force_n.unwrap()),
            (9.81, c.acceleration_m_s2),
        ] {
            let lattice = q / factor;
            let physical = lattice * factor;
            assert_close(physical, q);
        }
    }

    #[test]
    fn three_constructors_are_equivalent() {
        let mut a = base_params(UnitConstructor::FromResolutionAndLatticeVelocity);
        a.resolution = Some(100);
        a.lattice_velocity = Some(0.1);
        let ra = report(&a).unwrap();

        let mut b = base_params(UnitConstructor::FromResolutionAndRelaxationTime);
        b.resolution = Some(100);
        b.relaxation_time = Some(ra.lattice.tau);
        let rb = report(&b).unwrap();

        let mut c = base_params(UnitConstructor::FromRelaxationTimeAndLatticeVelocity);
        c.relaxation_time = Some(ra.lattice.tau);
        c.lattice_velocity = Some(0.1);
        let rc = report(&c).unwrap();

        for r in [&rb, &rc] {
            assert_close(r.lattice.dx, ra.lattice.dx);
            assert_close(r.lattice.dt, ra.lattice.dt);
            assert_close(r.lattice.nu_lattice, ra.lattice.nu_lattice);
            assert_close(r.lattice.tau, ra.lattice.tau);
            assert_close(r.lattice.u_char_lattice, ra.lattice.u_char_lattice);
            assert_close(r.dimensionless.reynolds, ra.dimensionless.reynolds);
            assert_close(r.dimensionless.mach, ra.dimensionless.mach);
            assert_close(
                r.dimensionless.grid_reynolds,
                ra.dimensionless.grid_reynolds,
            );
        }
    }

    #[test]
    fn known_case_anchors_match_hand_computed_values() {
        let mut schaefer_turek = base_params(UnitConstructor::FromResolutionAndLatticeVelocity);
        schaefer_turek.resolution = Some(100);
        schaefer_turek.lattice_velocity = Some(0.1);
        let r = report(&schaefer_turek).unwrap();
        assert_close(r.dimensionless.reynolds, 100.0);
        assert_close(r.lattice.tau, 0.8);
        assert_close(r.dimensionless.mach, 3.0f64.sqrt() * 0.1);
        assert_close(r.dimensionless.grid_reynolds, 1.0);

        let mut poiseuille = FlowParams {
            constructor: UnitConstructor::FromResolutionAndLatticeVelocity,
            characteristic_length: 0.01,
            characteristic_velocity: 0.02,
            kinematic_viscosity: 1.0e-6,
            density: Some(998.2),
            resolution: Some(100),
            lattice_velocity: Some(0.05),
            relaxation_time: None,
            end_time: None,
            end_step_count: None,
            gravity: None,
            reference_pressure: None,
            re_physical: None,
        };
        let r = report(&poiseuille).unwrap();
        assert_close(r.dimensionless.reynolds, 200.0);
        assert_close(r.lattice.tau, 0.575);
        assert_close(r.dimensionless.mach, 3.0f64.sqrt() * 0.05);
        assert_close(r.dimensionless.grid_reynolds, 2.0);

        poiseuille.gravity = Some([0.0, -9.81]);
        poiseuille.end_time = Some(0.1);
        let r = report(&poiseuille).unwrap();
        assert!(r.lattice.g_lat.unwrap()[1] < 0.0);
        assert_eq!(r.lattice.total_steps, Some(400));
    }

    #[test]
    fn threshold_boundaries_fire_exactly() {
        let tau_case = |tau: f64| {
            let mut p = base_params(UnitConstructor::FromResolutionAndRelaxationTime);
            p.resolution = Some(100);
            p.relaxation_time = Some(tau);
            report(&p).unwrap()
        };
        assert!(report_ids(&tau_case(0.5)).contains(&"TAU_UNSTABLE"));
        assert!(report_ids(&tau_case(0.5 + 1e-13)).contains(&"TAU_LOW"));
        assert!(!report_ids(&tau_case(TAU_LOW_WARN_THRESHOLD)).contains(&"TAU_LOW"));
        assert!(!report_ids(&tau_case(TAU_HIGH_WARN_THRESHOLD)).contains(&"TAU_HIGH"));
        assert!(report_ids(&tau_case(TAU_HIGH_WARN_THRESHOLD + 1e-12)).contains(&"TAU_HIGH"));

        let speed_case = |u: f64| {
            let mut p = base_params(UnitConstructor::FromResolutionAndLatticeVelocity);
            p.resolution = Some(100);
            p.lattice_velocity = Some(u);
            report(&p).unwrap()
        };
        assert!(!report_ids(&speed_case(LATTICE_SPEED_WARN_THRESHOLD)).contains(&"MACH_HIGH"));
        assert!(
            report_ids(&speed_case(LATTICE_SPEED_WARN_THRESHOLD + 1e-12)).contains(&"MACH_HIGH")
        );
        assert!(report_ids(&speed_case(MAX_SPEED)).contains(&"MACH_HIGH"));
        assert!(!report_ids(&speed_case(MAX_SPEED)).contains(&"MACH_HARD"));
        assert!(report_ids(&speed_case(MAX_SPEED + 1e-12)).contains(&"MACH_HARD"));

        let mut grid = base_params(UnitConstructor::FromResolutionAndLatticeVelocity);
        grid.characteristic_velocity = 15.0;
        grid.kinematic_viscosity = 0.001;
        grid.resolution = Some(100);
        grid.lattice_velocity = Some(0.1);
        assert!(!report_ids(&report(&grid).unwrap()).contains(&"GRID_RE_HIGH"));
        grid.characteristic_velocity = 15.0000000001;
        assert!(report_ids(&report(&grid).unwrap()).contains(&"GRID_RE_HIGH"));
    }

    #[test]
    fn density_missing_rejects_units_block() {
        let mut p = base_params(UnitConstructor::FromResolutionAndLatticeVelocity);
        p.resolution = Some(100);
        p.lattice_velocity = Some(0.1);
        p.density = None;
        let r = report(&p).unwrap();
        assert_eq!(r.verdict, UnitVerdict::Error);
        assert!(report_ids(&r).contains(&"DENSITY_MISSING"));
    }

    #[test]
    fn rounding_drift_is_reported_for_tau_velocity_constructor() {
        let mut p = base_params(UnitConstructor::FromRelaxationTimeAndLatticeVelocity);
        p.relaxation_time = Some(0.61);
        p.lattice_velocity = Some(0.073);
        let r = report(&p).unwrap();
        assert!(report_ids(&r).contains(&"RESOLUTION_ROUNDING"));
        let expected_n = (((0.61_f64 - 0.5) / 3.0) * 100.0 / 0.073).round();
        assert_close(r.lattice.tau, 3.0 * 0.073 * expected_n / 100.0 + 0.5);
    }

    #[test]
    fn warning_suggestion_feeds_back_to_ok() {
        let mut p = base_params(UnitConstructor::FromResolutionAndLatticeVelocity);
        p.characteristic_velocity = 20.0;
        p.kinematic_viscosity = 2.0e-4;
        p.resolution = Some(100);
        p.lattice_velocity = Some(0.2);
        let r = report(&p).unwrap();
        assert_eq!(r.verdict, UnitVerdict::Warn);
        let suggestion = r.suggestion.unwrap();
        p.resolution = Some(suggestion.resolution);
        p.lattice_velocity = Some(suggestion.lattice_velocity);
        let fixed = report(&p).unwrap();
        assert_eq!(fixed.verdict, UnitVerdict::Ok, "{:?}", report_ids(&fixed));
    }

    #[test]
    fn effective_viscosity_regime_is_echoed() {
        let mut p = base_params(UnitConstructor::FromResolutionAndLatticeVelocity);
        p.resolution = Some(100);
        p.lattice_velocity = Some(0.1);
        p.re_physical = Some(1_000.0);
        let r = report(&p).unwrap();
        assert_close(r.dimensionless.re_declared.unwrap(), 100.0);
        assert_eq!(r.dimensionless.re_physical, Some(1_000.0));
        assert_close(r.dimensionless.matching_ratio.unwrap(), 0.1);
        assert!(report_ids(&r).contains(&"EFFECTIVE_VISCOSITY_REGIME"));
    }

    #[test]
    fn no_units_block_keeps_legacy_scenario_shape() {
        let sc = Scenario {
            version: 0,
            name: "legacy".to_string(),
            grid: Grid {
                nx: 16,
                ny: 16,
                nz: 1,
            },
            physics: Physics {
                nu: 0.02,
                collision: CollisionSpec::Trt,
                force: [0.0, 0.0],
                gravity: None,
                precision: Precision::F64,
            },
            units: None,
            compute: None,
            wall: None,
            rotor: None,
            particles: None,
            edges: EdgesSpec {
                left: EdgeSpec::BounceBack,
                right: EdgeSpec::BounceBack,
                bottom: EdgeSpec::BounceBack,
                top: EdgeSpec::MovingWall { u: [0.1, 0.0] },
                front: None,
                back: None,
            },
            inlet_profile: None,
            obstacles: vec![],
            init: Default::default(),
            multiphase: None,
            run: RunSpec {
                steps: 1,
                stop_when_steady: None,
            },
            probes: vec![],
            outputs: vec![],
        };
        assert!(unit_report_for(&sc).unwrap().is_none());
        assert!(resolve(&sc).unwrap().is_none());
        assert!(!serde_json::to_string(&sc).unwrap().contains("\"units\""));
    }
}
