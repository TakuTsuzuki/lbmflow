//! Binned population balance model scaffolding for point-bubble aeration.

use crate::bubbles::BubbleError;
use crate::solver::UnsupportedReason;
use serde::{Deserialize, Serialize};

pub const DEFAULT_PBM_BIN_COUNT: usize = 20;
pub const PBM_ALPHA_G_MAX: f64 = 0.3;
pub const PBM_RE_BUBBLE_MAX: f64 = 800.0;
pub const PBM_KW_MIN: f64 = 0.0;
pub const PBM_KW_MAX: f64 = 1.0e6;

#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct PbmValidity {
    pub re_bubble_max: f64,
    pub alpha_g_max: f64,
    pub kolmogorov_weber_min: f64,
    pub kolmogorov_weber_max: f64,
}

impl Default for PbmValidity {
    fn default() -> Self {
        Self {
            re_bubble_max: PBM_RE_BUBBLE_MAX,
            alpha_g_max: PBM_ALPHA_G_MAX,
            kolmogorov_weber_min: PBM_KW_MIN,
            kolmogorov_weber_max: PBM_KW_MAX,
        }
    }
}

impl PbmValidity {
    pub fn validate(&self, state: PbmLocalState, kernel: &'static str) -> PbmResult<()> {
        if !(state.re_bubble.is_finite()
            && state.re_bubble >= 0.0
            && state.re_bubble <= self.re_bubble_max)
        {
            return Err(PbmError::out_of_validity_range(
                format!("{kernel} outside Re_bubble validity range"),
                format!(
                    "Re_bubble={:.6} must be in [0, {:.6}]",
                    state.re_bubble, self.re_bubble_max
                ),
            ));
        }
        if !(state.alpha_g.is_finite() && state.alpha_g >= 0.0 && state.alpha_g <= self.alpha_g_max)
        {
            return Err(PbmError::out_of_validity_range(
                format!("{kernel} outside alpha_g validity range"),
                format!(
                    "alpha_g={:.6} must be in [0, {:.6}]",
                    state.alpha_g, self.alpha_g_max
                ),
            ));
        }
        if !(state.kolmogorov_weber.is_finite()
            && state.kolmogorov_weber >= self.kolmogorov_weber_min
            && state.kolmogorov_weber <= self.kolmogorov_weber_max)
        {
            return Err(PbmError::out_of_validity_range(
                format!("{kernel} outside kW validity range"),
                format!(
                    "kW={:.6} must be in [{:.6}, {:.6}]",
                    state.kolmogorov_weber, self.kolmogorov_weber_min, self.kolmogorov_weber_max
                ),
            ));
        }
        Ok(())
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PbmLocalState {
    pub re_bubble: f64,
    pub alpha_g: f64,
    pub kolmogorov_weber: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct PbmBins {
    pub diameters_m: Vec<f64>,
    pub counts: Vec<f64>,
    pub cell_volume_m3: f64,
}

impl PbmBins {
    pub fn log_space(
        min_diameter_m: f64,
        max_diameter_m: f64,
        bin_count: Option<usize>,
        cell_volume_m3: f64,
    ) -> PbmResult<Self> {
        let n = match bin_count {
            Some(v) => v,
            None => DEFAULT_PBM_BIN_COUNT,
        };
        if n < 2 {
            return Err(PbmError::out_of_validity_range(
                "PBM requires at least two bins",
                "bin_count must be >= 2",
            ));
        }
        if !(min_diameter_m.is_finite()
            && max_diameter_m.is_finite()
            && min_diameter_m > 0.0
            && max_diameter_m > min_diameter_m
            && cell_volume_m3.is_finite()
            && cell_volume_m3 > 0.0)
        {
            return Err(PbmError::out_of_validity_range(
                "PBM bin construction requires positive diameters and cell volume",
                "0 < min_diameter_m < max_diameter_m and cell_volume_m3 > 0",
            ));
        }
        let log_min = min_diameter_m.ln();
        let log_max = max_diameter_m.ln();
        let mut diameters_m = Vec::with_capacity(n);
        for i in 0..n {
            let f = i as f64 / (n - 1) as f64;
            diameters_m.push((log_min + f * (log_max - log_min)).exp());
        }
        Ok(Self {
            diameters_m,
            counts: vec![0.0; n],
            cell_volume_m3,
        })
    }

    pub fn total_number(&self) -> f64 {
        self.counts.iter().sum()
    }

    pub fn d32_m(&self) -> PbmResult<f64> {
        let mut num = 0.0;
        let mut den = 0.0;
        for (&n, &d) in self.counts.iter().zip(self.diameters_m.iter()) {
            num += n * d * d * d;
            den += n * d * d;
        }
        if den <= 0.0 {
            return Err(PbmError::out_of_validity_range(
                "d32 is undefined for an empty PBM distribution",
                "sum(n_i d_i^2) must be > 0",
            ));
        }
        Ok(num / den)
    }

    pub fn interfacial_area_density_1_m(&self) -> f64 {
        let mut area = 0.0;
        for (&n, &d) in self.counts.iter().zip(self.diameters_m.iter()) {
            area += n * std::f64::consts::PI * d * d;
        }
        area / self.cell_volume_m3
    }

    pub fn apply_breakup<K: BreakupKernel>(
        &mut self,
        kernel: &K,
        dt_s: f64,
        state: PbmLocalState,
    ) -> PbmResult<()> {
        kernel.apply(self, dt_s, state)
    }

    pub fn apply_coalescence<K: CoalescenceKernel>(
        &mut self,
        kernel: &K,
        dt_s: f64,
        state: PbmLocalState,
    ) -> PbmResult<()> {
        kernel.apply(self, dt_s, state)
    }
}

pub trait BreakupKernel {
    fn name(&self) -> &'static str;
    fn validity(&self) -> PbmValidity;
    fn apply(&self, bins: &mut PbmBins, dt_s: f64, state: PbmLocalState) -> PbmResult<()>;
}

pub trait CoalescenceKernel {
    fn name(&self) -> &'static str;
    fn validity(&self) -> PbmValidity;
    fn apply(&self, bins: &mut PbmBins, dt_s: f64, state: PbmLocalState) -> PbmResult<()>;
}

#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub enum DisabledKernel {
    Disabled,
}

impl BreakupKernel for DisabledKernel {
    fn name(&self) -> &'static str {
        "disabled_breakup"
    }

    fn validity(&self) -> PbmValidity {
        PbmValidity::default()
    }

    fn apply(&self, _bins: &mut PbmBins, _dt_s: f64, state: PbmLocalState) -> PbmResult<()> {
        BreakupKernel::validity(self).validate(state, BreakupKernel::name(self))
    }
}

impl CoalescenceKernel for DisabledKernel {
    fn name(&self) -> &'static str {
        "disabled_coalescence"
    }

    fn validity(&self) -> PbmValidity {
        PbmValidity::default()
    }

    fn apply(&self, _bins: &mut PbmBins, _dt_s: f64, state: PbmLocalState) -> PbmResult<()> {
        CoalescenceKernel::validity(self).validate(state, CoalescenceKernel::name(self))
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct ConstantBreakup {
    pub rate_1_s: f64,
    pub validity: PbmValidity,
}

impl ConstantBreakup {
    pub fn new(rate_1_s: f64) -> PbmResult<Self> {
        if !(rate_1_s.is_finite() && rate_1_s >= 0.0) {
            return Err(PbmError::out_of_validity_range(
                "constant breakup rate must be finite and non-negative",
                "rate_1_s must be >= 0",
            ));
        }
        Ok(Self {
            rate_1_s,
            validity: PbmValidity::default(),
        })
    }
}

impl BreakupKernel for ConstantBreakup {
    fn name(&self) -> &'static str {
        "constant_breakup"
    }

    fn validity(&self) -> PbmValidity {
        self.validity
    }

    fn apply(&self, bins: &mut PbmBins, dt_s: f64, state: PbmLocalState) -> PbmResult<()> {
        self.validity.validate(state, self.name())?;
        validate_dt(dt_s)?;
        let mut delta = vec![0.0; bins.counts.len()];
        for i in 1..bins.counts.len() {
            let moved = bins.counts[i] * self.rate_1_s * dt_s;
            delta[i] -= moved;
            delta[i - 1] += moved;
        }
        for (n, d) in bins.counts.iter_mut().zip(delta.iter()) {
            *n += *d;
        }
        Ok(())
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct ConstantCoalescence {
    pub rate_1_s: f64,
    pub validity: PbmValidity,
}

impl ConstantCoalescence {
    pub fn new(rate_1_s: f64) -> PbmResult<Self> {
        if !(rate_1_s.is_finite() && rate_1_s >= 0.0) {
            return Err(PbmError::out_of_validity_range(
                "constant coalescence rate must be finite and non-negative",
                "rate_1_s must be >= 0",
            ));
        }
        Ok(Self {
            rate_1_s,
            validity: PbmValidity::default(),
        })
    }
}

impl CoalescenceKernel for ConstantCoalescence {
    fn name(&self) -> &'static str {
        "constant_coalescence"
    }

    fn validity(&self) -> PbmValidity {
        self.validity
    }

    fn apply(&self, bins: &mut PbmBins, dt_s: f64, state: PbmLocalState) -> PbmResult<()> {
        self.validity.validate(state, self.name())?;
        validate_dt(dt_s)?;
        let mut delta = vec![0.0; bins.counts.len()];
        let last = bins.counts.len() - 1;
        for i in 0..last {
            let moved = bins.counts[i] * self.rate_1_s * dt_s;
            delta[i] -= moved;
            delta[i + 1] += moved;
        }
        for (n, d) in bins.counts.iter_mut().zip(delta.iter()) {
            *n += *d;
        }
        Ok(())
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum FutureKernelHook {
    LuoSvendsenBreakup,
    PrinceBlanchCoalescence,
}

impl FutureKernelHook {
    pub fn reject(self) -> PbmError {
        PbmError {
            code: "pbm_kernel_not_implemented",
            message: format!("{self:?} is a design hook; implementation is deferred"),
            reason: UnsupportedReason::NotImplemented,
        }
    }
}

#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct PbmError {
    pub code: &'static str,
    pub message: String,
    pub reason: UnsupportedReason,
}

impl PbmError {
    pub fn out_of_validity_range(message: impl Into<String>, detail: impl Into<String>) -> Self {
        Self {
            code: "pbm_out_of_validity_range",
            message: message.into(),
            reason: UnsupportedReason::OutOfValidityRange {
                detail: detail.into(),
            },
        }
    }
}

impl From<BubbleError> for PbmError {
    fn from(value: BubbleError) -> Self {
        Self {
            code: "bubble_error",
            message: value.message,
            reason: value.reason,
        }
    }
}

impl std::fmt::Display for PbmError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{} ({:?})", self.message, self.reason)
    }
}

impl std::error::Error for PbmError {}

pub type PbmResult<T> = Result<T, PbmError>;

fn validate_dt(dt_s: f64) -> PbmResult<()> {
    if dt_s.is_finite() && dt_s >= 0.0 {
        Ok(())
    } else {
        Err(PbmError::out_of_validity_range(
            "PBM time step must be finite and non-negative",
            "dt_s must be >= 0",
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn state() -> PbmLocalState {
        PbmLocalState {
            re_bubble: 10.0,
            alpha_g: 0.01,
            kolmogorov_weber: 1.0,
        }
    }

    #[test]
    fn bin_total_number_is_conserved_by_placeholder_kernels() {
        let mut bins = PbmBins::log_space(1.0e-4, 4.0e-3, Some(8), 1.0e-6).unwrap();
        bins.counts[4] = 100.0;
        let before = bins.total_number();
        bins.apply_breakup(&ConstantBreakup::new(1.0).unwrap(), 0.1, state())
            .unwrap();
        bins.apply_coalescence(&ConstantCoalescence::new(1.0).unwrap(), 0.1, state())
            .unwrap();
        let after = bins.total_number();
        assert!((after - before).abs() < 1.0e-12);
    }

    #[test]
    fn d32_formula_matches_synthetic_distribution() {
        let mut bins = PbmBins {
            diameters_m: vec![1.0, 2.0],
            counts: vec![2.0, 1.0],
            cell_volume_m3: 1.0,
        };
        let got = bins.d32_m().unwrap();
        let expected = (2.0 * 1.0f64.powi(3) + 1.0 * 2.0f64.powi(3))
            / (2.0 * 1.0f64.powi(2) + 1.0 * 2.0f64.powi(2));
        assert!((got - expected).abs() < 1.0e-12);
        bins.counts = vec![0.0, 0.0];
        assert!(bins.d32_m().is_err());
    }

    #[test]
    fn breakup_shifts_distribution_smaller() {
        let mut bins = PbmBins::log_space(1.0e-4, 4.0e-3, Some(8), 1.0e-6).unwrap();
        bins.counts[6] = 10.0;
        bins.apply_breakup(&ConstantBreakup::new(1.0).unwrap(), 1.0, state())
            .unwrap();
        assert!(bins.counts[5] > 0.0);
        assert_eq!(bins.counts[6], 0.0);
    }

    #[test]
    fn coalescence_shifts_distribution_larger() {
        let mut bins = PbmBins::log_space(1.0e-4, 4.0e-3, Some(8), 1.0e-6).unwrap();
        bins.counts[1] = 10.0;
        bins.apply_coalescence(&ConstantCoalescence::new(1.0).unwrap(), 1.0, state())
            .unwrap();
        assert!(bins.counts[2] > 0.0);
        assert_eq!(bins.counts[1], 0.0);
    }

    #[test]
    fn disabled_kernel_preserves_bins() {
        let mut bins = PbmBins::log_space(1.0e-4, 4.0e-3, Some(8), 1.0e-6).unwrap();
        bins.counts[3] = 5.0;
        let before = bins.counts.clone();
        bins.apply_breakup(&DisabledKernel::Disabled, 1.0, state())
            .unwrap();
        bins.apply_coalescence(&DisabledKernel::Disabled, 1.0, state())
            .unwrap();
        assert_eq!(bins.counts, before);
    }
}
