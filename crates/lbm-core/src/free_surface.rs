//! Explicit top-boundary modes for phase-field bioprocess runs.

use serde::{Deserialize, Serialize};

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum TopBoundaryMode {
    ClosedLid,
    FreeSurface {
        engineering: bool,
    },
    DegassingOutlet {
        engineering: bool,
        gas_threshold: f64,
    },
}

impl Default for TopBoundaryMode {
    fn default() -> Self {
        Self::ClosedLid
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct DegassingLedger {
    pub gas_outflow_kg: f64,
    pub liquid_retention_delta_kg: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct FreeSurfaceError {
    pub message: String,
}

impl std::fmt::Display for FreeSurfaceError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.message)
    }
}

impl std::error::Error for FreeSurfaceError {}
