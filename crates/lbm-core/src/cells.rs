//! Serialization seam for future cell and microcarrier tracer state.

/// Checkpoint payload producer for cell-tracer-like particle state.
pub trait CellCheckpointSection {
    fn checkpoint_bytes(&self) -> Option<Vec<u8>>;
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct NoCellTracers;

impl CellCheckpointSection for NoCellTracers {
    fn checkpoint_bytes(&self) -> Option<Vec<u8>> {
        None
    }
}
