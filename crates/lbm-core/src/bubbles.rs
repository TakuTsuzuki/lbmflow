//! Point-bubble entities and conservative bookkeeping for engineering aeration.

use crate::solver::UnsupportedReason;
use serde::{Deserialize, Serialize};
use std::f64::consts::PI;

pub const POINT_BUBBLE_ALPHA_G_MAX: f64 = 0.3;

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Bubble {
    pub position: [f64; 3],
    pub velocity: [f64; 3],
    pub diameter_m: f64,
    pub gas_volume_m3: f64,
    pub age_s: f64,
    pub id: u64,
}

impl Bubble {
    pub fn new(
        position: [f64; 3],
        velocity: [f64; 3],
        diameter_m: f64,
        tank_diameter_m: f64,
        id: u64,
    ) -> BubbleResult<Self> {
        validate_bubble_diameter(diameter_m, tank_diameter_m)?;
        Ok(Self {
            position,
            velocity,
            diameter_m,
            gas_volume_m3: bubble_volume_from_diameter(diameter_m)?,
            age_s: 0.0,
            id,
        })
    }
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct BubbleSet {
    pub bubbles: Vec<Bubble>,
    pub next_id: u64,
}

impl BubbleSet {
    pub fn new() -> Self {
        Self {
            bubbles: Vec::new(),
            next_id: 0,
        }
    }

    pub fn inject_from_sparger(
        &mut self,
        injector: &mut SpargerBubbleInjector,
        t_start_s: f64,
        dt_s: f64,
    ) -> BubbleResult<usize> {
        injector.inject(self, t_start_s, dt_s)
    }

    pub fn point_holdup(
        &self,
        dims: [usize; 3],
        dx_m: f64,
        tank_diameter_m: f64,
    ) -> BubbleResult<Vec<f64>> {
        if !(dx_m.is_finite() && dx_m > 0.0) {
            return Err(BubbleError::out_of_validity_range(
                "grid spacing must be finite and positive",
                "dx_m must be > 0",
            ));
        }
        let n = dims[0] * dims[1] * dims[2];
        let mut alpha = vec![0.0; n];
        let cell_volume = dx_m * dx_m * dx_m;
        for bubble in &self.bubbles {
            validate_bubble_diameter(bubble.diameter_m, tank_diameter_m)?;
            if let Some(i) = cell_index_for_position(dims, dx_m, bubble.position) {
                alpha[i] += bubble.gas_volume_m3 / cell_volume;
                if alpha[i] > POINT_BUBBLE_ALPHA_G_MAX {
                    return Err(BubbleError::high_gas_holdup(alpha[i]));
                }
            }
        }
        Ok(alpha)
    }

    pub fn scatter_reaction_forces(
        &self,
        dims: [usize; 3],
        dx_m: f64,
        bubble_forces_n: &[[f64; 3]],
    ) -> BubbleResult<MomentumCouplingLedger> {
        if self.bubbles.len() != bubble_forces_n.len() {
            return Err(BubbleError::out_of_validity_range(
                "bubble force array length must match bubble count",
                "bubble_forces_n.len() must equal bubbles.len()",
            ));
        }
        if !(dx_m.is_finite() && dx_m > 0.0) {
            return Err(BubbleError::out_of_validity_range(
                "grid spacing must be finite and positive",
                "dx_m must be > 0",
            ));
        }
        let n = dims[0] * dims[1] * dims[2];
        let mut liquid_reaction_force_n = vec![[0.0; 3]; n];
        let mut total_bubble_force_n = [0.0; 3];
        let mut total_liquid_reaction_force_n = [0.0; 3];

        for (bubble, force) in self.bubbles.iter().zip(bubble_forces_n.iter()) {
            for a in 0..3 {
                total_bubble_force_n[a] += force[a];
            }
            let weights = regularized_weights(dims, dx_m, bubble.position);
            let mut weight_sum = 0.0;
            for (_, w) in &weights {
                weight_sum += *w;
            }
            if weight_sum == 0.0 {
                continue;
            }
            for (i, w) in weights {
                let normalized = w / weight_sum;
                for a in 0..3 {
                    let reaction = -force[a] * normalized;
                    liquid_reaction_force_n[i][a] += reaction;
                    total_liquid_reaction_force_n[a] += reaction;
                }
            }
        }
        for a in 0..3 {
            let residual = total_bubble_force_n[a] + total_liquid_reaction_force_n[a];
            let scale =
                total_bubble_force_n[a].abs() + total_liquid_reaction_force_n[a].abs() + 1.0;
            if residual.abs() > 1.0e-12 * scale {
                return Err(BubbleError::out_of_validity_range(
                    "bubble-liquid momentum ledger is not balanced",
                    format!("component {a}: F_bubble + F_liquid = {residual:e}, scale={scale:e}"),
                ));
            }
        }
        Ok(MomentumCouplingLedger {
            liquid_reaction_force_n,
            total_bubble_force_n,
            total_liquid_reaction_force_n,
        })
    }
}

impl Default for BubbleSet {
    fn default() -> Self {
        Self::new()
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct SpargerBubbleInjector {
    pub orifice_positions_m: Vec<[f64; 3]>,
    pub gas_volumetric_flow_m3_per_s: f64,
    pub bubble_diameter_m: f64,
    pub tank_diameter_m: f64,
    pub injected_gas_volume_m3: f64,
    next_injection_index: u64,
    next_orifice: usize,
}

impl SpargerBubbleInjector {
    pub fn new(
        orifice_positions_m: Vec<[f64; 3]>,
        gas_volumetric_flow_m3_per_s: f64,
        bubble_diameter_m: f64,
        tank_diameter_m: f64,
    ) -> BubbleResult<Self> {
        if orifice_positions_m.is_empty() {
            return Err(BubbleError::out_of_validity_range(
                "sparger must provide at least one orifice position",
                "orifice_positions_m must not be empty",
            ));
        }
        if !orifice_positions_m
            .iter()
            .all(|p| p.iter().all(|v| v.is_finite()))
        {
            return Err(BubbleError::out_of_validity_range(
                "sparger orifice positions must be finite",
                "all orifice coordinates must be finite",
            ));
        }
        if !(gas_volumetric_flow_m3_per_s.is_finite() && gas_volumetric_flow_m3_per_s > 0.0) {
            return Err(BubbleError::out_of_validity_range(
                "sparger gas flow must be finite and positive",
                "gas_volumetric_flow_m3_per_s must be > 0",
            ));
        }
        validate_bubble_diameter(bubble_diameter_m, tank_diameter_m)?;
        Ok(Self {
            orifice_positions_m,
            gas_volumetric_flow_m3_per_s,
            bubble_diameter_m,
            tank_diameter_m,
            injected_gas_volume_m3: 0.0,
            next_injection_index: 0,
            next_orifice: 0,
        })
    }

    pub fn bubble_volume_m3(&self) -> f64 {
        PI * self.bubble_diameter_m * self.bubble_diameter_m * self.bubble_diameter_m / 6.0
    }

    pub fn expected_injected_volume_m3(&self, elapsed_s: f64) -> BubbleResult<f64> {
        if !(elapsed_s.is_finite() && elapsed_s >= 0.0) {
            return Err(BubbleError::out_of_validity_range(
                "elapsed time must be finite and non-negative",
                "elapsed_s must be >= 0",
            ));
        }
        Ok(self.gas_volumetric_flow_m3_per_s * elapsed_s)
    }

    fn inject(&mut self, set: &mut BubbleSet, t_start_s: f64, dt_s: f64) -> BubbleResult<usize> {
        if !(t_start_s.is_finite() && t_start_s >= 0.0 && dt_s.is_finite() && dt_s >= 0.0) {
            return Err(BubbleError::out_of_validity_range(
                "injection time window must be finite and non-negative",
                "t_start_s and dt_s must be >= 0",
            ));
        }
        let vb = self.bubble_volume_m3();
        let t_end_s = t_start_s + dt_s;
        let mut inserted = 0usize;
        loop {
            let next_time_s =
                (self.next_injection_index as f64 + 1.0) * vb / self.gas_volumetric_flow_m3_per_s;
            if next_time_s > t_end_s {
                break;
            }
            if next_time_s > t_start_s {
                let position = self.orifice_positions_m[self.next_orifice];
                let bubble = Bubble::new(
                    position,
                    [0.0; 3],
                    self.bubble_diameter_m,
                    self.tank_diameter_m,
                    set.next_id,
                )?;
                set.next_id += 1;
                set.bubbles.push(bubble);
                self.injected_gas_volume_m3 += vb;
                inserted += 1;
                self.next_orifice += 1;
                if self.next_orifice == self.orifice_positions_m.len() {
                    self.next_orifice = 0;
                }
            }
            self.next_injection_index += 1;
        }
        Ok(inserted)
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct MomentumCouplingLedger {
    pub liquid_reaction_force_n: Vec<[f64; 3]>,
    pub total_bubble_force_n: [f64; 3],
    pub total_liquid_reaction_force_n: [f64; 3],
}

#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct BubbleError {
    pub code: &'static str,
    pub message: String,
    pub reason: UnsupportedReason,
}

impl BubbleError {
    pub fn out_of_validity_range(message: impl Into<String>, detail: impl Into<String>) -> Self {
        Self {
            code: "bubble_out_of_validity_range",
            message: message.into(),
            reason: UnsupportedReason::OutOfValidityRange {
                detail: detail.into(),
            },
        }
    }

    pub fn high_gas_holdup(alpha_g: f64) -> Self {
        Self {
            code: "high_gas_holdup_needs_continuum_mode",
            message: "point-bubble gas holdup exceeds the dilute-mode cap".to_string(),
            reason: UnsupportedReason::OutOfValidityRange {
                detail: format!(
                    "alpha_g_bubble={alpha_g:.6} exceeds {POINT_BUBBLE_ALPHA_G_MAX:.3}"
                ),
            },
        }
    }
}

impl std::fmt::Display for BubbleError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{} ({:?})", self.message, self.reason)
    }
}

impl std::error::Error for BubbleError {}

pub type BubbleResult<T> = Result<T, BubbleError>;

pub fn bubble_volume_from_diameter(diameter_m: f64) -> BubbleResult<f64> {
    if !(diameter_m.is_finite() && diameter_m > 0.0) {
        return Err(BubbleError::out_of_validity_range(
            "bubble diameter must be finite and positive",
            "diameter_m must be > 0",
        ));
    }
    Ok(PI * diameter_m * diameter_m * diameter_m / 6.0)
}

pub fn validate_bubble_diameter(diameter_m: f64, tank_diameter_m: f64) -> BubbleResult<()> {
    if !(tank_diameter_m.is_finite() && tank_diameter_m > 0.0) {
        return Err(BubbleError::out_of_validity_range(
            "tank diameter must be finite and positive",
            "tank_diameter_m must be > 0",
        ));
    }
    bubble_volume_from_diameter(diameter_m)?;
    if diameter_m > tank_diameter_m / 10.0 {
        return Err(BubbleError::out_of_validity_range(
            "point-bubble diameter exceeds sub-grid validity range",
            "diameter_m must be <= tank_diameter_m / 10",
        ));
    }
    Ok(())
}

pub fn cell_index_for_position(dims: [usize; 3], dx_m: f64, p: [f64; 3]) -> Option<usize> {
    if !p.iter().all(|v| v.is_finite()) {
        return None;
    }
    let x = (p[0] / dx_m).floor() as isize;
    let y = (p[1] / dx_m).floor() as isize;
    let z = (p[2] / dx_m).floor() as isize;
    if x < 0
        || y < 0
        || z < 0
        || x >= dims[0] as isize
        || y >= dims[1] as isize
        || z >= dims[2] as isize
    {
        return None;
    }
    Some(((z as usize * dims[1] + y as usize) * dims[0]) + x as usize)
}

fn regularized_weights(dims: [usize; 3], dx_m: f64, p: [f64; 3]) -> Vec<(usize, f64)> {
    let gx = p[0] / dx_m - 0.5;
    let gy = p[1] / dx_m - 0.5;
    let gz = p[2] / dx_m - 0.5;
    let x0 = gx.floor() as isize;
    let y0 = gy.floor() as isize;
    let z0 = gz.floor() as isize;
    let fx = gx - x0 as f64;
    let fy = gy - y0 as f64;
    let fz = gz - z0 as f64;
    let mut out = Vec::with_capacity(8);
    for dz in 0..=1 {
        for dy in 0..=1 {
            for dx in 0..=1 {
                let x = x0 + dx;
                let y = y0 + dy;
                let z = z0 + dz;
                if x < 0
                    || y < 0
                    || z < 0
                    || x >= dims[0] as isize
                    || y >= dims[1] as isize
                    || z >= dims[2] as isize
                {
                    continue;
                }
                let wx = if dx == 0 { 1.0 - fx } else { fx };
                let wy = if dy == 0 { 1.0 - fy } else { fy };
                let wz = if dz == 0 { 1.0 - fz } else { fz };
                let w = wx * wy * wz;
                if w > 0.0 {
                    out.push((
                        ((z as usize * dims[1] + y as usize) * dims[0]) + x as usize,
                        w,
                    ));
                }
            }
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bubble_volume_from_diameter_matches_pi_d3_over_6() {
        let d = 2.0e-3;
        let got = bubble_volume_from_diameter(d).unwrap();
        let expected = std::f64::consts::PI * d.powi(3) / 6.0;
        assert!((got - expected).abs() < 1.0e-18);
    }

    #[test]
    fn deterministic_injection_matches_schedule_and_ledger() {
        let mut set = BubbleSet::new();
        let d = 1.0e-3;
        let vb = bubble_volume_from_diameter(d).unwrap();
        let q = 10.0 * vb;
        let mut injector =
            SpargerBubbleInjector::new(vec![[0.01, 0.02, 0.03], [0.04, 0.05, 0.06]], q, d, 0.1)
                .unwrap();
        let first = set.inject_from_sparger(&mut injector, 0.0, 0.35).unwrap();
        let second = set.inject_from_sparger(&mut injector, 0.35, 0.65).unwrap();
        assert_eq!(first, 3);
        assert_eq!(second, 7);
        assert_eq!(set.bubbles.len(), 10);
        assert_eq!(set.bubbles[0].id, 0);
        assert_eq!(set.bubbles[1].position, [0.04, 0.05, 0.06]);
        let expected = injector.expected_injected_volume_m3(1.0).unwrap();
        let rel = ((injector.injected_gas_volume_m3 - expected) / expected).abs();
        assert!(rel <= 0.01, "ledger rel error {rel:e} must be <= 1%");
    }

    #[test]
    fn invalid_diameter_is_rejected() {
        assert!(validate_bubble_diameter(0.0, 1.0).is_err());
        assert!(validate_bubble_diameter(0.11, 1.0).is_err());
    }

    #[test]
    fn scatter_reaction_forces_balances_momentum() {
        let mut set = BubbleSet::new();
        set.bubbles
            .push(Bubble::new([0.15, 0.15, 0.15], [0.0; 3], 1.0e-3, 1.0, 0).unwrap());
        let ledger = set
            .scatter_reaction_forces([4, 4, 4], 0.1, &[[1.0, -2.0, 0.5]])
            .unwrap();
        for a in 0..3 {
            let residual = ledger.total_bubble_force_n[a] + ledger.total_liquid_reaction_force_n[a];
            assert!(residual.abs() < 1.0e-12);
        }
    }

    #[test]
    fn high_alpha_guard_rejects_point_bubble_holdup() {
        let mut set = BubbleSet::new();
        set.bubbles
            .push(Bubble::new([0.05, 0.05, 0.05], [0.0; 3], 0.09, 1.0, 0).unwrap());
        let err = set.point_holdup([2, 2, 2], 0.1, 1.0).unwrap_err();
        assert_eq!(err.code, "high_gas_holdup_needs_continuum_mode");
    }
}
