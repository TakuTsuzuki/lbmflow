//! Quantity-of-interest formula helpers for bioprocess reports.

use serde::{Deserialize, Serialize};

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PowerQoiInput {
    pub torque_n_m: f64,
    pub omega_rad_s: f64,
    pub rho_kg_m3: f64,
    pub impeller_diameter_m: f64,
    pub working_volume_m3: f64,
    pub discharge_flow_m3_s: Option<f64>,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PowerQoiResult {
    pub torque_n_m: f64,
    pub power_w: f64,
    pub rotational_speed_hz: f64,
    pub np: f64,
    pub p_over_v_w_m3: f64,
    pub nq: Option<f64>,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct SkippedQoi {
    pub qoi: String,
    pub reason: String,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MixingTimeResult {
    pub cv0: f64,
    pub t95_s: Option<f64>,
    pub t99_s: Option<f64>,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct KlaFitWindow {
    pub start_s: f64,
    pub end_s: f64,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum KlaFitMethod {
    DynamicGassingFit,
    SteadyStateIntegral,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct KlaDynamicFitResult {
    pub kla_1_per_s: f64,
    pub kla_1_per_hr: f64,
    pub fit_r2: f64,
    pub fitting_window_start_s: f64,
    pub fitting_window_end_s: f64,
    pub method: KlaFitMethod,
    pub ci95_1_per_s: Option<[f64; 2]>,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct KlaDynamicFitOutcome {
    pub result: Option<KlaDynamicFitResult>,
    pub skipped: Option<SkippedQoi>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct CompartmentCv {
    pub name: String,
    pub cv: Option<f64>,
    pub cell_count: usize,
}

pub fn power_qois(input: PowerQoiInput) -> Result<PowerQoiResult, &'static str> {
    if !(input.omega_rad_s.is_finite()
        && input.rho_kg_m3.is_finite()
        && input.impeller_diameter_m.is_finite()
        && input.working_volume_m3.is_finite()
        && input.rho_kg_m3 > 0.0
        && input.impeller_diameter_m > 0.0
        && input.working_volume_m3 > 0.0)
    {
        return Err("power QOI inputs must be finite and positive where dimensional");
    }
    let n_hz = input.omega_rad_s / std::f64::consts::TAU;
    if !(n_hz > 0.0) {
        return Err("rotational speed must be positive for Np/Nq");
    }
    let power_w = input.omega_rad_s * input.torque_n_m;
    let denom = input.rho_kg_m3 * n_hz.powi(3) * input.impeller_diameter_m.powi(5);
    let nq = input
        .discharge_flow_m3_s
        .map(|q| q / (n_hz * input.impeller_diameter_m.powi(3)));
    Ok(PowerQoiResult {
        torque_n_m: input.torque_n_m,
        power_w,
        rotational_speed_hz: n_hz,
        np: power_w / denom,
        p_over_v_w_m3: power_w / input.working_volume_m3,
        nq,
    })
}

pub fn scalar_cv(values: &[f64], include: &[bool]) -> Option<f64> {
    assert_eq!(values.len(), include.len());
    let mut n = 0usize;
    let mut sum = 0.0;
    for (&v, &inside) in values.iter().zip(include) {
        if inside {
            n += 1;
            sum += v;
        }
    }
    if n == 0 {
        return None;
    }
    let mean = sum / n as f64;
    if mean == 0.0 {
        return None;
    }
    let mut var = 0.0;
    for (&v, &inside) in values.iter().zip(include) {
        if inside {
            let d = v - mean;
            var += d * d;
        }
    }
    Some((var / n as f64).sqrt() / mean.abs())
}

pub fn mixing_time_from_cv(series: &[(f64, f64)]) -> Option<MixingTimeResult> {
    let &(_, cv0) = series.iter().find(|(_, cv)| cv.is_finite() && *cv > 0.0)?;
    let t95_threshold = 0.05 * cv0;
    let t99_threshold = 0.01 * cv0;
    let mut t95 = None;
    let mut t99 = None;
    for &(t, cv) in series {
        if !cv.is_finite() {
            continue;
        }
        if t95.is_none() && cv <= t95_threshold {
            t95 = Some(t);
        }
        if t99.is_none() && cv <= t99_threshold {
            t99 = Some(t);
        }
    }
    Some(MixingTimeResult {
        cv0,
        t95_s: t95,
        t99_s: t99,
    })
}

pub fn dynamic_gassing_window_default(series: &[(f64, f64)]) -> Option<KlaFitWindow> {
    let first = series.first()?.0;
    let last = series.last()?.0;
    if !(first.is_finite() && last.is_finite() && last > first) {
        return None;
    }
    Some(KlaFitWindow {
        start_s: first + 0.4 * (last - first),
        end_s: last,
    })
}

pub fn dynamic_gassing_kla_fit(
    series: &[(f64, f64)],
    c_star: f64,
    window: Option<KlaFitWindow>,
    steady_epsilon: f64,
) -> Result<KlaDynamicFitOutcome, &'static str> {
    if !(c_star.is_finite() && steady_epsilon.is_finite() && steady_epsilon >= 0.0) {
        return Err("kLa fit inputs must be finite; steady_epsilon must be >= 0");
    }
    if series.len() < 3 {
        return Ok(kla_skipped(
            "kLa",
            "at least three concentration samples are required",
        ));
    }
    let fit_window = match window.or_else(|| dynamic_gassing_window_default(series)) {
        Some(w) => w,
        None => return Ok(kla_skipped("kLa", "fitting window could not be inferred")),
    };
    if !(fit_window.start_s.is_finite()
        && fit_window.end_s.is_finite()
        && fit_window.end_s > fit_window.start_s)
    {
        return Ok(kla_skipped("kLa", "invalid fitting window"));
    }

    let mut points = Vec::new();
    let mut max_abs_slope = 0.0;
    let mut prev: Option<(f64, f64)> = None;
    for &(t, c) in series {
        if !(t.is_finite() && c.is_finite()) {
            continue;
        }
        if let Some((pt, pc)) = prev {
            let dt = t - pt;
            if dt > 0.0 {
                let slope = (c - pc).abs() / dt;
                if slope > max_abs_slope {
                    max_abs_slope = slope;
                }
            }
        }
        prev = Some((t, c));
        if t < fit_window.start_s || t > fit_window.end_s {
            continue;
        }
        let gap = c_star - c;
        if gap <= 0.0 || !gap.is_finite() {
            continue;
        }
        points.push((t, gap.ln()));
    }
    if max_abs_slope <= steady_epsilon {
        return Ok(kla_skipped(
            "kLa",
            "oxygen concentration is steady within epsilon",
        ));
    }
    if points.len() < 3 {
        return Ok(kla_skipped(
            "kLa",
            "fitting window has fewer than three non-equilibrium samples",
        ));
    }
    let n = points.len() as f64;
    let mean_t = points.iter().map(|p| p.0).sum::<f64>() / n;
    let mean_y = points.iter().map(|p| p.1).sum::<f64>() / n;
    let mut sxx = 0.0;
    let mut sxy = 0.0;
    let mut syy = 0.0;
    for &(t, y) in &points {
        let dt = t - mean_t;
        let dy = y - mean_y;
        sxx += dt * dt;
        sxy += dt * dy;
        syy += dy * dy;
    }
    if sxx <= 0.0 || syy <= 0.0 {
        return Ok(kla_skipped("kLa", "fitting window has no dynamic range"));
    }
    let slope = sxy / sxx;
    let intercept = mean_y - slope * mean_t;
    let kla = -slope;
    if !(kla.is_finite() && kla >= 0.0) {
        return Ok(kla_skipped("kLa", "fitted kLa is negative or non-finite"));
    }
    let mut sse = 0.0;
    for &(t, y) in &points {
        let residual = y - (intercept + slope * t);
        sse += residual * residual;
    }
    let r2 = 1.0 - sse / syy;
    if r2 < 0.9 {
        return Ok(kla_skipped("kLa", "dynamic gassing fit R2 is below 0.9"));
    }
    let ci95 = if points.len() > 2 {
        let sigma2 = sse / (points.len() as f64 - 2.0);
        let se_slope = (sigma2 / sxx).sqrt();
        let half = 1.96 * se_slope;
        let lo = kla - half;
        let lower = if lo < 0.0 { 0.0 } else { lo };
        Some([lower, kla + half])
    } else {
        None
    };
    Ok(KlaDynamicFitOutcome {
        result: Some(KlaDynamicFitResult {
            kla_1_per_s: kla,
            kla_1_per_hr: kla * 3600.0,
            fit_r2: r2,
            fitting_window_start_s: fit_window.start_s,
            fitting_window_end_s: fit_window.end_s,
            method: KlaFitMethod::DynamicGassingFit,
            ci95_1_per_s: ci95,
        }),
        skipped: None,
    })
}

fn kla_skipped(qoi: &str, reason: &str) -> KlaDynamicFitOutcome {
    KlaDynamicFitOutcome {
        result: None,
        skipped: Some(SkippedQoi {
            qoi: qoi.to_string(),
            reason: reason.to_string(),
        }),
    }
}

pub fn compartment_cv(
    dims: [usize; 3],
    values: &[f64],
    solid: &[bool],
    impeller_center_z: Option<f64>,
) -> Vec<CompartmentCv> {
    let n = dims[0] * dims[1] * dims[2];
    assert_eq!(values.len(), n);
    assert_eq!(solid.len(), n);
    let mut top = vec![false; n];
    let mut bulk = vec![false; n];
    let mut near = vec![false; n];
    for z in 0..dims[2] {
        for y in 0..dims[1] {
            for x in 0..dims[0] {
                let i = idx(dims, x, y, z);
                if solid[i] {
                    continue;
                }
                let zn = (z as f64 + 0.5) / dims[2] as f64;
                if zn >= 2.0 / 3.0 {
                    top[i] = true;
                } else if impeller_center_z
                    .is_some_and(|z_near| ((z as f64 + 0.5) - z_near).abs() <= dims[2] as f64 / 6.0)
                {
                    near[i] = true;
                } else {
                    bulk[i] = true;
                }
            }
        }
    }
    [("top", top), ("bulk", bulk), ("near_impeller", near)]
        .into_iter()
        .map(|(name, mask)| CompartmentCv {
            name: name.to_string(),
            cv: scalar_cv(values, &mask),
            cell_count: mask.iter().filter(|&&v| v).count(),
        })
        .collect()
}

fn idx(dims: [usize; 3], x: usize, y: usize, z: usize) -> usize {
    (z * dims[1] + y) * dims[0] + x
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn zero_torque_gives_zero_power_and_np() {
        let q = power_qois(PowerQoiInput {
            torque_n_m: 0.0,
            omega_rad_s: std::f64::consts::TAU,
            rho_kg_m3: 1000.0,
            impeller_diameter_m: 0.1,
            working_volume_m3: 0.01,
            discharge_flow_m3_s: None,
        })
        .unwrap();
        assert_eq!(q.power_w, 0.0);
        assert_eq!(q.np, 0.0);
        assert!(q.nq.is_none());
    }

    #[test]
    fn power_formula_and_sign_convention_are_explicit() {
        let q = power_qois(PowerQoiInput {
            torque_n_m: 2.0,
            omega_rad_s: std::f64::consts::TAU * 3.0,
            rho_kg_m3: 1000.0,
            impeller_diameter_m: 0.5,
            working_volume_m3: 2.0,
            discharge_flow_m3_s: Some(0.75),
        })
        .unwrap();
        assert!((q.power_w - std::f64::consts::TAU * 6.0).abs() < 1.0e-12);
        assert!((q.p_over_v_w_m3 - q.power_w / 2.0).abs() < 1.0e-12);
        assert!((q.nq.unwrap() - 0.75 / (3.0 * 0.5f64.powi(3))).abs() < 1.0e-12);
        let reverse = power_qois(PowerQoiInput {
            torque_n_m: -2.0,
            ..q_input()
        })
        .unwrap();
        assert!(reverse.power_w < 0.0);
    }

    fn q_input() -> PowerQoiInput {
        PowerQoiInput {
            torque_n_m: 1.0,
            omega_rad_s: std::f64::consts::TAU,
            rho_kg_m3: 1000.0,
            impeller_diameter_m: 0.5,
            working_volume_m3: 1.0,
            discharge_flow_m3_s: None,
        }
    }

    #[test]
    fn synthetic_cv_series_finds_t95_t99() {
        let series = [(0.0, 1.0), (1.0, 0.2), (2.0, 0.05), (3.0, 0.009)];
        let m = mixing_time_from_cv(&series).unwrap();
        assert_eq!(m.cv0, 1.0);
        assert_eq!(m.t95_s, Some(2.0));
        assert_eq!(m.t99_s, Some(3.0));
    }

    #[test]
    fn uniform_scalar_has_zero_cv_and_no_mixing_time() {
        let values = vec![1.0; 8];
        let mask = vec![true; 8];
        assert_eq!(scalar_cv(&values, &mask), Some(0.0));
        assert!(mixing_time_from_cv(&[(0.0, 0.0), (1.0, 0.0)]).is_none());
    }

    #[test]
    fn synthetic_exponential_uptake_fit_recovers_kla_within_5_percent() {
        let c_star = 1.0;
        let kla = 0.04;
        let series: Vec<_> = (0..=100)
            .map(|i| {
                let t = i as f64;
                (t, c_star - (c_star - 0.1) * (-kla * t).exp())
            })
            .collect();
        let fit = dynamic_gassing_kla_fit(&series, c_star, None, 1.0e-12)
            .unwrap()
            .result
            .unwrap();
        let rel = (fit.kla_1_per_s - kla).abs() / kla;
        assert!(rel <= 0.05, "rel={rel}");
        assert!(fit.fit_r2 > 0.999);
        assert!((fit.kla_1_per_hr - kla * 3600.0).abs() < 1.0e-9);
    }

    #[test]
    fn bad_kla_fit_window_is_rejected_with_skip_reason() {
        let series = [(0.0, 0.0), (1.0, 0.1), (2.0, 0.2)];
        let outcome = dynamic_gassing_kla_fit(
            &series,
            1.0,
            Some(KlaFitWindow {
                start_s: 2.0,
                end_s: 1.0,
            }),
            0.0,
        )
        .unwrap();
        assert!(outcome.result.is_none());
        assert!(outcome
            .skipped
            .unwrap()
            .reason
            .contains("invalid fitting window"));
    }

    #[test]
    fn equilibrium_kla_data_is_rejected() {
        let series = [(0.0, 1.0), (1.0, 1.0), (2.0, 1.0), (3.0, 1.0)];
        let outcome = dynamic_gassing_kla_fit(&series, 1.0, None, 1.0e-12).unwrap();
        assert!(outcome.result.is_none());
        assert!(outcome.skipped.unwrap().reason.contains("steady"));
    }

    #[test]
    fn compartment_aggregation_counts_regions() {
        let dims = [4, 4, 6];
        let n = dims[0] * dims[1] * dims[2];
        let values: Vec<f64> = (0..n).map(|i| 1.0 + i as f64).collect();
        let solid = vec![false; n];
        let comps = compartment_cv(dims, &values, &solid, Some(2.5));
        assert_eq!(comps.len(), 3);
        assert!(comps.iter().all(|c| c.cell_count > 0));
        assert!(comps
            .iter()
            .any(|c| c.name == "near_impeller" && c.cv.is_some()));
    }
}

// ============================================================================
// Checkpoint scaffolding (BCFD-102, brought forward from Bundle Q merge)
// ============================================================================

/// Minimal deterministic accumulator snapshot for BCFD-102.
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct QoiAccumulatorSnapshot {
    pub name: String,
    pub count: u64,
    pub sum: f64,
    pub min: f64,
    pub max: f64,
    pub reservoir: Vec<f64>,
}

/// Container serialized by solver checkpoints when QOI statistics exist.
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct QoiCheckpointState {
    pub accumulators: Vec<QoiAccumulatorSnapshot>,
}

impl QoiCheckpointState {
    pub fn is_empty(&self) -> bool {
        self.accumulators.is_empty()
    }
}
