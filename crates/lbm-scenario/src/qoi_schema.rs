//! Stable bioprocess QOI output-file contract.

use serde::{Deserialize, Serialize};

pub const QOI_BUNDLE_JSON: &str = "qoi.json";
pub const QOI_POWER_CSV: &str = "qoi_power.csv";
pub const QOI_MIXING_CSV: &str = "qoi_mixing.csv";
pub const QOI_GAS_CSV: &str = "qoi_gas.csv";
pub const QOI_KLA_CSV: &str = "qoi_kla.csv";
pub const QOI_SHEAR_CSV: &str = "qoi_shear.csv";
pub const QOI_CELLS_CSV: &str = "qoi_cells.csv";

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub struct QoiOutputSchema {
    pub bundle_json: &'static str,
    pub time_series_csv: Vec<&'static str>,
}

pub fn bioprocess_qoi_output_schema() -> QoiOutputSchema {
    QoiOutputSchema {
        bundle_json: QOI_BUNDLE_JSON,
        time_series_csv: vec![
            QOI_POWER_CSV,
            QOI_MIXING_CSV,
            QOI_GAS_CSV,
            QOI_KLA_CSV,
            QOI_SHEAR_CSV,
            QOI_CELLS_CSV,
        ],
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn qoi_schema_lists_required_files() {
        let schema = bioprocess_qoi_output_schema();
        assert_eq!(schema.bundle_json, "qoi.json");
        assert!(schema.time_series_csv.contains(&"qoi_power.csv"));
        assert!(schema.time_series_csv.contains(&"qoi_cells.csv"));
    }
}
