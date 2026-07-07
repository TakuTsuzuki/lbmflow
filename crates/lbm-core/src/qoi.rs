//! Checkpoint scaffolding for running QOI statistics.

use serde::{Deserialize, Serialize};

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
