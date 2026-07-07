//! Physical material-property fields for bioprocess model coupling.

use crate::params::{PhaseFieldMixtureParams, ViscosityInterpolation};
use crate::real::Real;

/// Compact-core material-property fields.
#[derive(Clone, Debug, PartialEq)]
pub struct MaterialFields {
    pub rho_phys: Vec<f64>,
    pub mu_phys: Vec<f64>,
    pub nu_phys: Vec<f64>,
    pub sigma: Vec<f64>,
    pub alpha_liquid: Vec<f64>,
    pub alpha_gas: Vec<f64>,
}

/// One compact-core material sample supplied by localized geometry builders.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaterialSample {
    pub rho_phys: f64,
    pub mu_phys: f64,
    pub nu_phys: f64,
    pub sigma: f64,
    pub alpha_liquid: f64,
    pub alpha_gas: f64,
}

impl MaterialFields {
    pub fn new(n: usize) -> Self {
        Self {
            rho_phys: vec![0.0; n],
            mu_phys: vec![0.0; n],
            nu_phys: vec![0.0; n],
            sigma: vec![0.0; n],
            alpha_liquid: vec![0.0; n],
            alpha_gas: vec![0.0; n],
        }
    }

    pub fn len(&self) -> usize {
        self.rho_phys.len()
    }

    pub fn is_empty(&self) -> bool {
        self.rho_phys.is_empty()
    }

    pub fn update_phase_field_mixture<T: Real>(
        &mut self,
        phi: &[T],
        params: &PhaseFieldMixtureParams,
    ) {
        assert_eq!(
            self.len(),
            phi.len(),
            "phi length must match material fields"
        );
        for (i, &p) in phi.iter().enumerate() {
            let phi = p.as_f64();
            let rho = rho_interpolation(phi, params.rho_gas, params.rho_liquid);
            let mu = match params.viscosity_interpolation {
                ViscosityInterpolation::Harmonic => {
                    harmonic_viscosity(phi, params.mu_gas, params.mu_liquid)
                }
                ViscosityInterpolation::Linear => {
                    linear_viscosity(phi, params.mu_gas, params.mu_liquid)
                }
            };
            self.rho_phys[i] = rho;
            self.mu_phys[i] = mu;
            self.nu_phys[i] = mu / rho;
            self.sigma[i] = params.sigma;
            self.alpha_liquid[i] = phi;
            self.alpha_gas[i] = 1.0 - phi;
        }
    }

    pub fn rho_phys(&self) -> &[f64] {
        &self.rho_phys
    }

    pub fn mu_phys(&self) -> &[f64] {
        &self.mu_phys
    }

    pub fn nu_phys(&self) -> &[f64] {
        &self.nu_phys
    }

    pub fn sigma(&self) -> &[f64] {
        &self.sigma
    }

    pub fn alpha_liquid(&self) -> &[f64] {
        &self.alpha_liquid
    }

    pub fn alpha_gas(&self) -> &[f64] {
        &self.alpha_gas
    }

    pub fn set_sample(&mut self, i: usize, sample: MaterialSample) {
        self.rho_phys[i] = sample.rho_phys;
        self.mu_phys[i] = sample.mu_phys;
        self.nu_phys[i] = sample.nu_phys;
        self.sigma[i] = sample.sigma;
        self.alpha_liquid[i] = sample.alpha_liquid;
        self.alpha_gas[i] = sample.alpha_gas;
    }
}

#[inline]
pub fn rho_interpolation(phi: f64, rho_gas: f64, rho_liquid: f64) -> f64 {
    rho_gas + phi * (rho_liquid - rho_gas)
}

#[inline]
pub fn harmonic_viscosity(phi: f64, mu_gas: f64, mu_liquid: f64) -> f64 {
    1.0 / ((1.0 - phi) / mu_gas + phi / mu_liquid)
}

#[inline]
pub fn linear_viscosity(phi: f64, mu_gas: f64, mu_liquid: f64) -> f64 {
    mu_gas + phi * (mu_liquid - mu_gas)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::params::MaterialModel;

    #[test]
    fn rho_interpolation_endpoints_match_gas_and_liquid() {
        assert_eq!(rho_interpolation(0.0, 1.2, 998.0), 1.2);
        assert_eq!(rho_interpolation(1.0, 1.2, 998.0), 998.0);
    }

    #[test]
    fn rho_mu_interpolation_endpoints_still_match() {
        let params =
            match MaterialModel::phase_field_mixture(1.2, 998.0, 1.8e-5, 1.0e-3, 0.072).unwrap() {
                MaterialModel::PhaseFieldMixture(params) => params,
                _ => unreachable!(),
            };
        let mut fields = MaterialFields::new(2);
        fields.update_phase_field_mixture(&[0.0f64, 1.0], &params);
        assert_eq!(fields.rho_phys(), &[1.2, 998.0]);
        assert_eq!(fields.mu_phys(), &[1.8e-5, 1.0e-3]);
    }

    #[test]
    fn harmonic_viscosity_matches_endpoints() {
        assert_eq!(harmonic_viscosity(0.0, 1.8e-5, 1.0e-3), 1.8e-5);
        assert_eq!(harmonic_viscosity(1.0, 1.8e-5, 1.0e-3), 1.0e-3);
    }

    #[test]
    fn linear_viscosity_opt_in_matches_endpoints() {
        assert_eq!(linear_viscosity(0.0, 1.8e-5, 1.0e-3), 1.8e-5);
        assert_eq!(linear_viscosity(1.0, 1.8e-5, 1.0e-3), 1.0e-3);
    }

    #[test]
    fn negative_density_or_viscosity_rejected() {
        assert!(MaterialModel::phase_field_mixture(-1.0, 998.0, 1.8e-5, 1.0e-3, 0.072).is_err());
        assert!(MaterialModel::phase_field_mixture(1.2, 998.0, -1.8e-5, 1.0e-3, 0.072).is_err());
    }

    #[test]
    fn rejects_density_ratio_over_1000_engineering() {
        let err =
            MaterialModel::phase_field_mixture(1.0, 1000.1, 1.0e-5, 1.0e-3, 0.01).unwrap_err();
        assert_eq!(err.parameter, "rho_liquid/rho_gas");
        assert!(err
            .message
            .contains("DENSITY_RATIO_BEYOND_ENGINEERING_TIER"));
    }
}
