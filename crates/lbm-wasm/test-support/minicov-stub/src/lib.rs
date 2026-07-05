//! Stub used only to let Cargo resolve wasm-bindgen-test metadata offline.
//!
//! The real `minicov` crate is pulled by wasm-bindgen-test only under
//! `cfg(all(target_arch = "wasm32", wasm_bindgen_unstable_test_coverage))`.
//! LBMFlow does not enable that unstable coverage cfg for its wasm smoke.
