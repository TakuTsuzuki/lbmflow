use clap::ValueEnum;
use serde::Serialize;
use std::fs;
use std::path::Path;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, ValueEnum)]
#[serde(rename_all = "snake_case")]
pub enum OutputMode {
    Gather,
    PerRank,
}

pub fn select_output_mode(world_size: usize, requested: Option<OutputMode>) -> OutputMode {
    requested.unwrap_or(if world_size > 4 {
        OutputMode::PerRank
    } else {
        OutputMode::Gather
    })
}

#[derive(Clone, Debug, PartialEq, Serialize)]
#[serde(rename_all = "snake_case")]
#[allow(dead_code)]
pub struct KlaQoiOutput {
    pub kla_1_per_s: Option<f64>,
    pub kla_1_per_hr: Option<f64>,
    pub fit_r2: Option<f64>,
    pub fitting_window_start_s: Option<f64>,
    pub fitting_window_end_s: Option<f64>,
    pub method: String,
    pub ci95_1_per_s: Option<[f64; 2]>,
    pub skipped_reason: Option<String>,
}

impl From<lbm_core::qoi::KlaDynamicFitOutcome> for KlaQoiOutput {
    fn from(value: lbm_core::qoi::KlaDynamicFitOutcome) -> Self {
        if let Some(result) = value.result {
            Self {
                kla_1_per_s: Some(result.kla_1_per_s),
                kla_1_per_hr: Some(result.kla_1_per_hr),
                fit_r2: Some(result.fit_r2),
                fitting_window_start_s: Some(result.fitting_window_start_s),
                fitting_window_end_s: Some(result.fitting_window_end_s),
                method: "dynamic_gassing_fit".to_string(),
                ci95_1_per_s: result.ci95_1_per_s,
                skipped_reason: None,
            }
        } else {
            Self {
                kla_1_per_s: None,
                kla_1_per_hr: None,
                fit_r2: None,
                fitting_window_start_s: None,
                fitting_window_end_s: None,
                method: "dynamic_gassing_fit".to_string(),
                ci95_1_per_s: None,
                skipped_reason: value.skipped.map(|s| s.reason),
            }
        }
    }
}

#[allow(dead_code)]
pub fn write_kla_qoi_json(
    outcome: lbm_core::qoi::KlaDynamicFitOutcome,
    out_dir: &Path,
) -> std::io::Result<String> {
    let name = "kla_qoi.json";
    let output = KlaQoiOutput::from(outcome);
    let bytes = serde_json::to_vec_pretty(&output)
        .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
    fs::write(out_dir.join(name), bytes)?;
    Ok(name.to_string())
}
