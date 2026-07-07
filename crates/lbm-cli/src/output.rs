use clap::ValueEnum;
use serde::Serialize;

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
