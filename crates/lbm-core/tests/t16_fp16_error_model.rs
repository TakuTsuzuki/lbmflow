#![cfg(feature = "gpu")]
//! T16 lane 1.7: FP16 deviation-storage error-model characterization.
//!
//! These tests are deliberately order-band tests, not validation-grade physics
//! references. `GpuStorage::F16` stores deviation populations (`f_q - w_q`) as
//! IEEE f16 and widens to f32 for arithmetic. A stored f16 value has relative
//! round-off `eps16 ~= 2^-10`; if the rest-state physical populations were
//! stored directly, a quiescent field would carry coherent O(eps16 * rho)
//! bias. In deviation storage the rest state is exactly zero, so the coherent
//! bias scale is set by the flow/deviation amplitude instead: O(eps16 *
//! u_scale) or O(eps16 * delta_rho), with per-step uncorrelated storage error
//! accumulating as a random walk, O(sqrt(N_steps) * eps16 * scale). This is the
//! model documented in `kernels.rs` and frozen for T16 in `docs/PHYSICS.md`.
//!
//! Consequences used below:
//! - A nearly uniform rho perturbation of amplitude 1e-6 should leave density
//!   and mass drift on the `sqrt(steps) * eps16 * 1e-6` scale, not
//!   `sqrt(steps) * eps16 * rho`.
//! - A short diffusive-TGV f16-vs-f32 velocity difference should be derived
//!   from the `sqrt(steps) * eps16 * u0` population scale. The asserted field
//!   RMS uses the D2Q9 moment projection factor: the velocity signal lives in
//!   weighted odd population pairs, and the RMS norm over the sinusoidal field
//!   is correspondingly smaller than the raw population-scale envelope.

use lbm_core::gpu::{GpuInitError, GpuStorage, KernelCfg};
use lbm_core::prelude::*;
use std::f64::consts::PI;
use std::sync::Arc;

type Gpu2 = Solver<D2Q9, f32, WgpuBackend<D2Q9>, LocalPeriodic>;

const EPS16: f64 = 1.0 / 1024.0;

fn shader_f16_ctx_or_skip() -> Option<Arc<GpuContext>> {
    match GpuContext::new_with_shader_f16(true) {
        Ok(ctx) => Some(ctx),
        Err(GpuInitError::NoAdapter) => {
            eprintln!("GPU-PENDING: skipping T16 FP16 error-model test; no usable GPU adapter");
            None
        }
        Err(e) => panic!("T16 FP16 error-model tests require SHADER_F16: {e}"),
    }
}

fn gpu_solver(ctx: &Arc<GpuContext>, spec: &GlobalSpec<f32>, storage: GpuStorage) -> Gpu2 {
    let backend = WgpuBackend::<D2Q9>::with_config(ctx.clone(), KernelCfg { storage });
    Solver::new(spec, &[], &[], [1, 1, 1], backend, LocalPeriodic)
}

fn velocity_l2_rms_diff(a: &mut Gpu2, b: &mut Gpu2) -> (f64, f64, f64) {
    let aux = a.gather_ux();
    let auy = a.gather_uy();
    let bux = b.gather_ux();
    let buy = b.gather_uy();
    assert_eq!(aux.len(), bux.len());

    let mut sx = 0.0f64;
    let mut sy = 0.0f64;
    for i in 0..aux.len() {
        let dx = aux[i] as f64 - bux[i] as f64;
        let dy = auy[i] as f64 - buy[i] as f64;
        sx += dx * dx;
        sy += dy * dy;
    }
    let n = aux.len() as f64;
    (
        (sx + sy).sqrt() / (2.0 * n).sqrt(),
        (sx / n).sqrt(),
        (sy / n).sqrt(),
    )
}

#[test]
#[ignore = "gpu-required"]
fn t16_fp16_uniform_density_roundoff_tracks_deviation_scale() {
    let Some(ctx) = shader_f16_ctx_or_skip() else {
        return;
    };

    let n = 128usize;
    let steps = 10_000usize;
    let rho_amp = 1.0e-6f64;
    let k = 2.0 * PI / n as f64;
    let spec = GlobalSpec::<f32> {
        dims: [n, n, 1],
        // High enough to damp the seeded acoustic density mode over 10^4
        // steps, leaving the f16 deviation-storage floor as the measured scale.
        nu: 0.2,
        periodic: [true, true, false],
        ..Default::default()
    };

    let mut sim = gpu_solver(&ctx, &spec, GpuStorage::F16);
    sim.init_with(|x, _y, _z| ((1.0 + rho_amp * (k * x as f64).sin()) as f32, [0.0; 3]));
    let initial_mass = sim.total_mass() as f64;
    sim.run(steps);

    let rho = sim.gather_rho();
    let final_mass = sim.total_mass() as f64;
    let cells = (n * n) as f64;
    let final_mean = final_mass / cells;
    let mass_drift_rel = (final_mass - initial_mass).abs() / initial_mass.abs().max(1.0);
    let max_centered_rho = rho
        .iter()
        .map(|&r| (r as f64 - final_mean).abs())
        .fold(0.0f64, f64::max);

    let predicted_cell_error = (steps as f64).sqrt() * EPS16 * rho_amp;
    let physical_population_bound = (steps as f64).sqrt() * EPS16;
    println!(
        "T16 FP16 uniform-density model: steps={steps}, predicted_cell={predicted_cell_error:.9e}, \
         mass_drift_rel={mass_drift_rel:.9e}, max_centered_rho={max_centered_rho:.9e}, \
         physical_population_bound={physical_population_bound:.9e}"
    );

    assert!(
        mass_drift_rel <= 10.0 * predicted_cell_error,
        "uniform f16 mass drift {mass_drift_rel:.9e} exceeds 10x deviation-storage prediction {predicted_cell_error:.9e}"
    );
    assert!(
        max_centered_rho <= 10.0 * predicted_cell_error,
        "uniform f16 density residual {max_centered_rho:.9e} exceeds 10x deviation-storage prediction {predicted_cell_error:.9e}"
    );
    assert!(
        max_centered_rho <= 1.0e-3 * physical_population_bound,
        "density residual {max_centered_rho:.9e} is not separated from rho-scale f16 storage bound {physical_population_bound:.9e}"
    );
}

#[test]
#[ignore = "gpu-required"]
fn t16_fp16_short_tgv_diff_matches_random_walk_model() {
    let Some(ctx) = shader_f16_ctx_or_skip() else {
        return;
    };

    let n = 128usize;
    let steps = 100usize;
    let nu = 0.02f64;
    let u0 = 1.28 / n as f64;
    let k = 2.0 * PI / n as f64;
    let spec = GlobalSpec::<f32> {
        dims: [n, n, 1],
        nu,
        periodic: [true, true, false],
        ..Default::default()
    };
    let init = |x: usize, y: usize, _z: usize| {
        let (xf, yf) = (k * x as f64, k * y as f64);
        let rho = 1.0 - 3.0 * u0 * u0 / 4.0 * ((2.0 * xf).cos() + (2.0 * yf).cos());
        (
            rho as f32,
            [
                (-u0 * xf.cos() * yf.sin()) as f32,
                (u0 * xf.sin() * yf.cos()) as f32,
                0.0,
            ],
        )
    };

    let mut f32s = gpu_solver(&ctx, &spec, GpuStorage::F32);
    let mut f16s = gpu_solver(&ctx, &spec, GpuStorage::F16);
    f32s.init_with(init);
    f16s.init_with(init);
    f32s.run(steps);
    f16s.run(steps);

    let (l2_rms, ux_rms, uy_rms) = velocity_l2_rms_diff(&mut f16s, &mut f32s);
    let raw_population_prediction = (steps as f64).sqrt() * EPS16 * u0;
    // D2Q9 projection from stored odd population deviations to field-RMS
    // velocity error. This is an order-model coefficient, not a fitted band:
    // axial/diagonal equilibrium deviations are weighted by 1/9 and 1/36,
    // the velocity moment subtracts opposite pairs, and the TGV field RMS
    // contributes another sinusoidal factor. The resulting coefficient is
    // O(1/10..1/100); 1/32 is the midpoint used for a factor-3 gate.
    let predicted_field_rms = raw_population_prediction / 32.0;
    let ratio = l2_rms / predicted_field_rms.max(1.0e-30);
    let component_ratio = ux_rms / uy_rms.max(1.0e-30);
    println!(
        "T16 FP16 short-TGV model: steps={steps}, u0={u0:.9e}, raw_population_prediction={raw_population_prediction:.9e}, \
         predicted_field_rms={predicted_field_rms:.9e}, \
         measured_l2_rms={l2_rms:.9e}, ratio={ratio:.6}, ux_rms={ux_rms:.9e}, uy_rms={uy_rms:.9e}"
    );

    assert!(
        (1.0 / 3.0..=3.0).contains(&ratio),
        "short-TGV f16-vs-f32 velocity L2 RMS {l2_rms:.9e} is not within factor 3 of projected random-walk prediction {predicted_field_rms:.9e} (raw population scale {raw_population_prediction:.9e})"
    );
    assert!(
        (0.5..=2.0).contains(&component_ratio),
        "short-TGV degradation should preserve x/y symmetry; ux_rms/uy_rms = {component_ratio:.6}"
    );
}
