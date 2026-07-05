//! `GpuSolver`: batch-first driver of the [`WgpuBackend`] with the V1 step
//! sequence and a host-side mirror for setup and accessors.
//!
//! Shape (GPU_EVALUATION.md §4 / adoption condition (b)):
//!
//! - **Setup is host-side**: an embedded monolithic `Solver<L, f32,
//!   CpuScalar, LocalPeriodic>` provides `init_with`, wall rims, obstacles,
//!   probes and inlet profiles with the exact V1 code paths; `run` uploads
//!   the staged state once.
//! - **`run(n)` is batched**: each step records `collide → stream(fused) →
//!   swap → apply_open_faces → update_moments` through the [`Backend`]
//!   trait; ops are submitted every `submit_chunk` steps *without waiting*
//!   (one submit per run for small `n`). No per-step CPU sync anywhere —
//!   the measured 9x cliff.
//! - **Readback is explicit**: [`GpuSolver::sync`] copies populations,
//!   moments and the probe force back into the host mirror; the accessors
//!   (`u`, `rho`, `total_mass`, `gather_*`, `probed_force`, …) take
//!   `&mut self` and sync lazily, so the only blocking points are visible
//!   in the API.
//!
//! Known divergence (documented, V1-mechanics related): editing geometry
//! *between* GPU runs re-uploads from the host mirror whose ping-pong
//! partner buffer (`ftmp`) is not reconstructed from the device, so the
//! ConvectiveOutflow previous-value state restarts from the mirror's stale
//! copy. Set up the scenario before stepping (the V1-suite pattern) and the
//! trajectories match T14-tight.

use std::sync::Arc;
use std::time::Instant;

use crate::backend::{Backend, CellRange, CpuScalar, HostMoments};
use crate::halo::LocalPeriodic;
use crate::lattice::{Face, Lattice, D2Q9};
use crate::params::StepParams;
use crate::solver::{GlobalSpec, Solver};
use crate::subdomain::Subdomain;

use super::backend::{GpuContext, GpuError, GpuFields, WgpuBackend};

/// GPU (wgpu) time-evolution driver over a monolithic domain, f32.
pub struct GpuSolver<L: Lattice = D2Q9> {
    inner: Solver<L, f32, CpuScalar, LocalPeriodic>,
    backend: WgpuBackend<L>,
    fields: GpuFields,
    sub: Subdomain,
    params: StepParams<f32>,
    time: u64,
    probed: [f32; 3],
    host_dirty: bool,
    device_ahead: bool,
}

impl<L: Lattice> GpuSolver<L> {
    /// Build a solver over the whole grid (monolithic decomposition), with
    /// compact global `solid` / `wall_u` arrays exactly like
    /// [`Solver::new`].
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
        let inner = Solver::new(
            spec,
            solid,
            wall_u,
            [1, 1, 1],
            CpuScalar::default(),
            LocalPeriodic,
        );
        let (omega_p, omega_m) = spec.collision.omegas(spec.nu);
        let params = StepParams {
            omega_p,
            omega_m,
            force: spec.force,
            faces: spec.faces,
        };
        let backend = WgpuBackend::<L>::new(ctx);
        let sub = inner.sub(0).clone();
        let fields = backend.try_alloc_with_options(&sub, solid.iter().any(|&s| s), false)?;
        Ok(Self {
            inner,
            backend,
            fields,
            sub,
            params,
            time: 0,
            probed: [0.0; 3],
            host_dirty: true,
            device_ahead: false,
        })
    }

    /// Steps per queue submit during `run` (default 200; anything on the
    /// ≥10 plateau of the proto's submit-granularity table is equivalent).
    pub fn set_submit_chunk(&mut self, chunk: usize) {
        self.backend.set_submit_chunk(chunk);
    }

    // ------------------------------------------------------------------
    // Setup (host mirror; re-uploaded on the next run)
    // ------------------------------------------------------------------

    fn edit_host(&mut self) -> &mut Solver<L, f32, CpuScalar, LocalPeriodic> {
        self.sync();
        self.host_dirty = true;
        &mut self.inner
    }

    /// Second-order consistent initialisation (V1 `init_with`).
    pub fn init_with(&mut self, init: impl Fn(usize, usize, usize) -> (f32, [f32; 3])) {
        self.edit_host().init_with(init);
    }

    /// Mark a global cell solid (half-way bounce-back obstacle).
    pub fn set_solid(&mut self, x: usize, y: usize, z: usize) {
        self.edit_host().set_solid(x, y, z);
    }

    /// Select the probed solid cells (momentum-exchange force).
    pub fn set_force_probe(&mut self, pred: impl Fn(usize, usize, usize) -> bool) {
        self.edit_host().set_force_probe(pred);
    }

    /// Per-node inlet velocity profile on a `Velocity` face.
    pub fn set_inlet_profile(&mut self, face: Face, values: &[[f32; 3]]) {
        self.edit_host().set_inlet_profile(face, values);
    }

    /// Per-cell body force added to the uniform force (compact core layout).
    pub fn set_force_field(&mut self, field: Vec<[f32; 3]>) {
        let n = self.inner.dims().iter().product::<usize>();
        assert_eq!(field.len(), n, "force field must cover the whole grid");
        self.edit_host().fields_mut(0).force_field = Some(field);
    }

    // ------------------------------------------------------------------
    // Time evolution
    // ------------------------------------------------------------------

    /// Advance `steps` steps on the GPU: encode everything, submit in
    /// chunks, **no wait** (readback APIs wait when asked).
    pub fn run(&mut self, steps: usize) {
        self.try_run(steps).expect("GPU run failed");
    }

    /// Fallible variant of [`Self::run`].
    pub fn try_run(&mut self, steps: usize) -> Result<(), GpuError> {
        if steps == 0 {
            return Ok(());
        }
        if self.host_dirty {
            self.backend
                .upload(&self.sub, &mut self.fields, self.inner.fields(0));
            self.host_dirty = false;
        }
        let mut submitted_steps = 0usize;
        for step in 0..steps {
            self.backend
                .collide(&self.sub, &mut self.fields, &self.params);
            let _ = self.backend.stream(
                &self.sub,
                &mut self.fields,
                &self.params,
                CellRange::full(&self.sub),
            );
            self.backend.swap(&mut self.fields);
            self.backend
                .apply_open_faces(&self.sub, &mut self.fields, &self.params);
            self.backend
                .update_moments(&self.sub, &mut self.fields, &self.params);
            let is_last = step + 1 == steps;
            let chunk = self.backend.submit_chunk();
            if submitted_steps + 1 >= chunk || is_last {
                let measured_steps = submitted_steps + 1;
                let calibrating = !self.backend.submit_chunk_calibrated();
                let t = Instant::now();
                self.backend.try_flush(&self.fields)?;
                if calibrating {
                    self.backend.context().try_wait_idle()?;
                    self.backend
                        .calibrate_submit_chunk(measured_steps, t.elapsed());
                }
                submitted_steps = 0;
            } else {
                submitted_steps += 1;
            }
        }
        self.time += steps as u64;
        self.device_ahead = true;
        Ok(())
    }

    /// Block until all submitted GPU work completed (benchmarking hook —
    /// `run` itself never waits).
    pub fn wait_idle(&self) {
        self.try_wait_idle().expect("GPU wait failed");
    }

    /// Fallible variant of [`Self::wait_idle`].
    pub fn try_wait_idle(&self) -> Result<(), GpuError> {
        self.backend.context().try_wait_idle()
    }

    // ------------------------------------------------------------------
    // Explicit readback
    // ------------------------------------------------------------------

    /// Read the device state back into the host mirror (populations,
    /// moments, probe force). All accessors call this lazily; it is a no-op
    /// when the mirror is current.
    pub fn sync(&mut self) {
        self.try_sync().expect("GPU sync failed");
    }

    /// Fallible variant of [`Self::sync`].
    pub fn try_sync(&mut self) -> Result<(), GpuError> {
        if !self.device_ahead {
            return Ok(());
        }
        let mut hm = HostMoments::default();
        let (f, probed) = self.backend.try_read_sync(&self.fields, &mut hm)?;
        self.probed = probed;
        let host = self.inner.fields_mut(0);
        let g = host.geom;
        let (nx, ny, np) = (g.core[0], g.core[1], g.n_padded());
        let n = nx * ny;
        for q in 0..L::Q {
            for y in 0..ny {
                for x in 0..nx {
                    host.f[q * np + g.pidx(x, y, 0)] = f[q * n + y * nx + x];
                }
            }
        }
        host.rho.copy_from_slice(&hm.rho);
        host.ux.copy_from_slice(&hm.ux);
        host.uy.copy_from_slice(&hm.uy);
        self.device_ahead = false;
        Ok(())
    }

    // ------------------------------------------------------------------
    // Accessors / diagnostics (auto-sync)
    // ------------------------------------------------------------------

    /// Completed time steps.
    pub fn time(&self) -> u64 {
        self.time
    }
    /// Global grid extents.
    pub fn dims(&self) -> [usize; 3] {
        self.inner.dims()
    }

    /// Density at a global cell.
    pub fn rho(&mut self, x: usize, y: usize, z: usize) -> f32 {
        self.sync();
        self.inner.rho(x, y, z)
    }
    /// Velocity at a global cell (physical, half-force corrected).
    pub fn u(&mut self, x: usize, y: usize, z: usize) -> [f32; 3] {
        self.sync();
        self.inner.u(x, y, z)
    }
    /// Momentum-exchange force on the probed solids during the most recent
    /// step.
    pub fn probed_force(&mut self) -> [f32; 3] {
        self.sync();
        self.probed
    }
    /// Total mass over fluid cells (V1 f64 accumulation).
    pub fn total_mass(&mut self) -> f32 {
        self.sync();
        self.inner.total_mass()
    }
    /// Total physical momentum over fluid cells.
    pub fn total_momentum(&mut self) -> [f32; 3] {
        self.sync();
        self.inner.total_momentum()
    }
    /// Global density field (compact layout).
    pub fn gather_rho(&mut self) -> Vec<f32> {
        self.sync();
        self.inner.gather_rho()
    }
    /// Global x-velocity field.
    pub fn gather_ux(&mut self) -> Vec<f32> {
        self.sync();
        self.inner.gather_ux()
    }
    /// Global y-velocity field.
    pub fn gather_uy(&mut self) -> Vec<f32> {
        self.sync();
        self.inner.gather_uy()
    }
    /// Global deviation-population plane `q` (compact layout).
    pub fn gather_f(&mut self, q: usize) -> Vec<f32> {
        self.sync();
        self.inner.gather_f(q)
    }
    /// Whether a global cell is solid.
    pub fn is_solid(&self, x: usize, y: usize, z: usize) -> bool {
        self.inner.is_solid(x, y, z)
    }
}
