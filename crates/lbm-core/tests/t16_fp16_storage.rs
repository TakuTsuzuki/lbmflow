#![cfg(feature = "gpu")]
//! T16: FP16 storage characterization.
//!
//! These tests are ignored in the sandbox continuation because the order
//! requires measured f16-vs-f32 degradation bands to be characterized before
//! freezing assertions. PM will unignore/freeze on an adapter that exposes
//! SHADER_F16 reliably.

use lbm_core::prelude::*;

#[test]
#[ignore = "requires SHADER_F16 adapter characterization; freeze measured bands before enabling"]
fn t16_tgv2d_f16_storage_degradation_vs_f32_gpu() {
    let _ctx = GpuContext::new_with_shader_f16(true)
        .expect("T16 requires a GPU adapter with SHADER_F16");
    panic!("T16 bands are BENCH-PENDING until PM adapter run freezes values");
}

#[test]
#[ignore = "requires SHADER_F16 adapter characterization; freeze measured bands before enabling"]
fn t16_cavity2d_f16_storage_degradation_vs_f32_gpu() {
    let _ctx = GpuContext::new_with_shader_f16(true)
        .expect("T16 requires a GPU adapter with SHADER_F16");
    panic!("T16 bands are BENCH-PENDING until PM adapter run freezes values");
}
