//! Deprecated GPU convenience wrapper.
//!
//! `GpuSolver` used to own a separate GPU step sequence. B-1 routes GPU
//! stepping through the common [`Solver`] orchestrator; this type remains only
//! as a compatibility surface for the current examples/tests/CLI-facing code.

#![allow(deprecated)]

use std::sync::Arc;

use crate::halo::LocalPeriodic;
use crate::lattice::{Face, Lattice, D2Q9};
use crate::solver::{Diverged, GlobalSpec, Solver};

use super::backend::{GpuContext, GpuError, WgpuBackend};

/// Deprecated compatibility alias around the unified solver.
#[deprecated(
    since = "0.1.0",
    note = "use Solver<L, f32, WgpuBackend<L>, LocalPeriodic>; GpuSolver is a thin compatibility wrapper"
)]
pub struct GpuSolver<L: Lattice = D2Q9> {
    inner: Solver<L, f32, WgpuBackend<L>, LocalPeriodic>,
}

impl<L: Lattice> GpuSolver<L> {
    /// Build a monolithic GPU solver.
    pub fn new(
        spec: &GlobalSpec<f32>,
        solid: &[bool],
        wall_u: &[[f32; 3]],
        ctx: Arc<GpuContext>,
    ) -> Self {
        Self::try_new(spec, solid, wall_u, ctx).expect("GPU solver initialization failed")
    }

    /// Fallible variant of [`Self::new`].
    pub fn try_new(
        spec: &GlobalSpec<f32>,
        solid: &[bool],
        wall_u: &[[f32; 3]],
        ctx: Arc<GpuContext>,
    ) -> Result<Self, GpuError> {
        let backend = WgpuBackend::<L>::new(ctx);
        let inner = Solver::try_new(spec, solid, wall_u, [1, 1, 1], backend, LocalPeriodic)
            .map_err(|e| GpuError::Spec(e.to_string()))?;
        Ok(Self { inner })
    }

    /// Steps per queue submit during `run`.
    pub fn set_submit_chunk(&mut self, chunk: usize) {
        self.inner.backend_mut().set_submit_chunk(chunk);
    }

    /// Number of queue submissions issued by the underlying backend.
    pub fn submissions(&self) -> u64 {
        self.inner.backend().submissions()
    }

    /// Enable or disable on-device WALE LES omega generation.
    ///
    /// Default is disabled. When enabled, each GPU step refreshes the
    /// pre-collision moments, computes a per-cell `omega_plus` field on the
    /// device, then collides with that field.
    pub fn set_wale(&mut self, enabled: bool) {
        let base_omega = (1.0 / self.inner.tau()) as f32;
        let backend = self.inner.backend() as *const WgpuBackend<L>;
        let fields = self.inner.backend_fields_mut(0);
        // SAFETY: `backend` and `fields` are disjoint fields inside `Solver`.
        // The backend method flushes queued work before replacing WALE buffers.
        unsafe { (&*backend).set_wale_enabled(fields, enabled, base_omega) };
    }

    /// Second-order consistent initialisation.
    pub fn init_with(&mut self, init: impl Fn(usize, usize, usize) -> (f32, [f32; 3])) {
        self.inner.backend().clear_cached_moment_upload_marker();
        self.inner.init_with(init);
    }

    /// Mark a global cell solid.
    pub fn set_solid(&mut self, x: usize, y: usize, z: usize) {
        self.inner.set_solid(x, y, z);
    }

    /// Select the probed solid cells.
    pub fn set_force_probe(&mut self, pred: impl Fn(usize, usize, usize) -> bool) {
        self.inner.set_force_probe(pred);
    }

    /// Per-node inlet velocity profile on a `Velocity` face.
    pub fn set_inlet_profile(&mut self, face: Face, values: &[[f32; 3]]) {
        self.inner.set_inlet_profile(face, values);
    }

    /// Per-cell body force added to the uniform force (compact core layout).
    pub fn set_force_field(&mut self, field: Vec<[f32; 3]>) {
        let n = self.inner.dims().iter().product::<usize>();
        assert_eq!(field.len(), n, "force field must cover the whole grid");
        self.inner.set_body_force_field_values(&field);
        self.inner.backend().cache_next_upload_moments_once();
    }

    /// Advance `steps` steps on the GPU through the unified solver.
    pub fn run(&mut self, steps: usize) {
        self.inner.run(steps);
    }

    /// Fallible variant kept for compatibility. The unified solver panics on
    /// backend errors through the existing `Backend` trait, so this wrapper has
    /// no additional fallible path until B-2 separates submit/read errors.
    pub fn try_run(&mut self, steps: usize) -> Result<(), GpuError> {
        self.run(steps);
        Ok(())
    }

    /// Block until all submitted GPU work completed.
    pub fn wait_idle(&self) {
        self.try_wait_idle().expect("GPU wait failed");
    }

    /// Fallible variant of [`Self::wait_idle`].
    pub fn try_wait_idle(&self) -> Result<(), GpuError> {
        self.inner.backend().context().try_wait_idle()
    }

    /// Advance `steps` steps with a periodic non-finite watchdog.
    pub fn run_guarded(&mut self, steps: usize, check_every: usize) -> Result<(), Diverged> {
        self.inner.run_guarded(steps, check_every)
    }

    /// Read the device state back into the host staging mirror.
    pub fn sync(&mut self) {
        self.inner.sync_host();
    }

    /// Fallible variant kept for compatibility. See [`Self::try_run`].
    pub fn try_sync(&mut self) -> Result<(), GpuError> {
        self.sync();
        Ok(())
    }

    /// Completed time steps.
    pub fn time(&self) -> u64 {
        self.inner.time()
    }

    /// Global grid extents.
    pub fn dims(&self) -> [usize; 3] {
        self.inner.dims()
    }

    /// Density at a global cell.
    pub fn rho(&mut self, x: usize, y: usize, z: usize) -> f32 {
        self.inner.rho(x, y, z)
    }

    /// Velocity at a global cell.
    pub fn u(&mut self, x: usize, y: usize, z: usize) -> [f32; 3] {
        self.inner.u(x, y, z)
    }

    /// Momentum-exchange force on probed solids during the most recent step.
    pub fn probed_force(&mut self) -> [f32; 3] {
        self.inner.probed_force()
    }

    /// Explicit readback of the momentum-exchange force on probed solids
    /// during the most recent completed step.
    pub fn read_probed_force(&self) -> [f32; 3] {
        self.inner.read_probed_force()
    }

    /// Total mass over fluid cells.
    pub fn total_mass(&mut self) -> f32 {
        self.inner.total_mass()
    }

    /// Total physical momentum over fluid cells.
    pub fn total_momentum(&mut self) -> [f32; 3] {
        self.inner.total_momentum()
    }

    /// Global density field.
    pub fn gather_rho(&mut self) -> Vec<f32> {
        self.inner.gather_rho()
    }

    /// Global x-velocity field.
    pub fn gather_ux(&mut self) -> Vec<f32> {
        self.inner.gather_ux()
    }

    /// Global y-velocity field.
    pub fn gather_uy(&mut self) -> Vec<f32> {
        self.inner.gather_uy()
    }

    /// Global z-velocity field.
    pub fn gather_uz(&mut self) -> Vec<f32> {
        self.inner.gather_uz()
    }

    /// Global strain-rate tensor through the host staging mirror.
    pub fn gather_strain_rate(&mut self) -> Vec<[f32; 6]> {
        self.sync();
        self.inner.gather_strain_rate()
    }

    /// Global shear-rate invariant through the host staging mirror.
    pub fn gather_shear_rate(&mut self) -> Vec<f32> {
        self.sync();
        self.inner.gather_shear_rate()
    }

    /// Current on-device WALE `omega_plus` field in compact global order.
    pub fn gather_wale_omega(&mut self) -> Vec<f32> {
        let base_omega = (1.0 / self.inner.tau()) as f32;
        self.inner
            .backend()
            .try_read_wale_omega(self.inner.backend_fields(0), base_omega)
            .expect("GPU WALE omega readback failed")
    }

    /// Global deviation-population plane `q`.
    pub fn gather_f(&mut self, q: usize) -> Vec<f32> {
        self.sync();
        self.inner.gather_f(q)
    }

    /// Whether a global cell is solid.
    pub fn is_solid(&self, x: usize, y: usize, z: usize) -> bool {
        self.inner.is_solid(x, y, z)
    }
}
