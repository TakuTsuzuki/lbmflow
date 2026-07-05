//! `WgpuBackend`: the [`Backend`] trait on wgpu compute (D2Q9, f32).
//!
//! Execution model (GPU_EVALUATION.md §3.4 — the submit-granularity trap):
//! the trait's phase methods **record** dispatches into a per-fields op list
//! and never talk to the GPU synchronously. Ops are materialised into one
//! command encoder / one compute pass and submitted (a) every
//! `submit_chunk` steps, without waiting, and (b) before any readback. The
//! CPU blocks only inside the explicit readback APIs (`read_moments`,
//! `read_f`, `read_probed_force`, `reduce`).
//!
//! Phase mapping (same observable state evolution as `CpuScalar`, one fused
//! kernel instead of two passes):
//!
//! - `collide` — records nothing; arms the fused dispatch (the `step` kernel
//!   performs collide+push-stream in one pass, preserving the CPU's `S∘C`
//!   step order — see `wgsl.rs`).
//! - `stream` — records `clear_probe` (when probing) + the fused dispatch.
//!   Returns zeros: the probe force accumulates on-device and is fetched
//!   through [`WgpuBackend::read_probed_force`] (explicit readback), holding
//!   the most recent step's value like V1 `probed_force`.
//! - `swap` — flips the ping-pong index (bind groups pre-built per parity).
//! - `apply_open_faces` — records the per-face `bc` dispatches in
//!   `Face::ALL` order (CPU order).
//! - `update_moments` — lazy: the fused kernel re-derives the moments it
//!   needs in-kernel from the same pre-collide state the CPU caches, so the
//!   device moment buffers are only refreshed by `read_moments`. Doubles as
//!   the step-end hook for chunked submission.
//! - `reduce` — explicit readback of the populations plus the exact V1
//!   f64 accumulation loop on the host (wgpu has no f64; doing the loop
//!   host-side keeps the diagnostic convention bit-compatible with
//!   `CpuScalar::reduce` for identical fields).
//!
//! Field storage is compact (`cell = y*nx + x`, no halo ring): periodic
//! wrap happens in-kernel, which is the monolithic-subdomain configuration.
//! Multi-part GPU decompositions (device-side halo exchange) are out of
//! scope here (M-D).

use std::cell::RefCell;
use std::marker::PhantomData;
use std::sync::Arc;

use crate::backend::{Backend, CellRange, HostMoments};
use crate::fields::SoaFields;
use crate::lattice::{Face, Lattice};
use crate::params::{FaceBC, Reduction, StepParams};
use crate::subdomain::Subdomain;

use super::wgsl;

/// Shared wgpu instance/adapter/device/queue.
///
/// Requests the adapter's real buffer limits (a 2048² D2Q9 f32 plane set is
/// a 151 MB binding, past the 128 MiB WebGPU default).
pub struct GpuContext {
    /// wgpu device.
    pub device: wgpu::Device,
    /// Submission queue.
    pub queue: wgpu::Queue,
    /// Adapter description (diagnostics).
    pub adapter_info: wgpu::AdapterInfo,
}

impl GpuContext {
    /// Create a context on the highest-performance adapter, or `None` when
    /// no usable GPU exists.
    pub fn new() -> Option<Arc<Self>> {
        let instance = wgpu::Instance::default();
        let adapter = pollster::block_on(instance.request_adapter(&wgpu::RequestAdapterOptions {
            power_preference: wgpu::PowerPreference::HighPerformance,
            ..Default::default()
        }))
        .ok()?;
        let adapter_info = adapter.get_info();
        let al = adapter.limits();
        let mut limits = wgpu::Limits::default();
        limits.max_storage_buffer_binding_size = al.max_storage_buffer_binding_size;
        limits.max_buffer_size = al.max_buffer_size;
        let (device, queue) = pollster::block_on(adapter.request_device(&wgpu::DeviceDescriptor {
            label: Some("lbm-core-gpu"),
            required_features: wgpu::Features::empty(),
            required_limits: limits,
            memory_hints: wgpu::MemoryHints::Performance,
            trace: wgpu::Trace::Off,
        }))
        .ok()?;
        Some(Arc::new(Self {
            device,
            queue,
            adapter_info,
        }))
    }

    /// Block until all submitted work completed.
    pub fn wait_idle(&self) {
        self.device
            .poll(wgpu::PollType::Wait)
            .expect("device poll failed");
    }
}

struct Pipelines {
    step: wgpu::ComputePipeline,
    moments: wgpu::ComputePipeline,
    bc: wgpu::ComputePipeline,
    clear_probe: wgpu::ComputePipeline,
}

/// Recorded (not yet submitted) dispatch.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Op {
    ClearProbe,
    /// Fused collide+stream, reading buffer parity `bg`.
    Fused { bg: usize },
    /// Open-face BC on `face`, operating on buffer parity `bg`.
    Bc { face: usize, bg: usize },
    /// Moments refresh from buffer parity `bg`.
    Moments { bg: usize },
}

/// Interior-mutable recorder state (readback methods take `&Fields`).
struct RecState {
    cur: usize,
    ops: Vec<Op>,
    steps_recorded: usize,
    pending_collide: bool,
    /// Written Params uniform (asserts step-parameter stability per run).
    params_words: Option<[u32; 12]>,
    /// Written per-face BC uniforms.
    bc_words: Option<[[u32; 32]; 6]>,
    /// Bumped per fused dispatch; invalidates the readback cache.
    generation: u64,
    /// Cached population readback (generation, data).
    f_cache: Option<(u64, Vec<f32>)>,
}

/// Device-resident fields of one (monolithic) subdomain.
pub struct GpuFields {
    nx: u32,
    ny: u32,
    n: usize,
    f: [wgpu::Buffer; 2],
    stash: [wgpu::Buffer; 2],
    mask: wgpu::Buffer,
    wall_u: wgpu::Buffer,
    force_field: wgpu::Buffer,
    rho: wgpu::Buffer,
    ux: wgpu::Buffer,
    uy: wgpu::Buffer,
    probe_acc: wgpu::Buffer,
    params_ub: wgpu::Buffer,
    bc_ub: [wgpu::Buffer; 6],
    profiles: [wgpu::Buffer; 6],
    staging: wgpu::Buffer,
    staging_small: wgpu::Buffer,
    fused_bg: [wgpu::BindGroup; 2],
    moments_bg: [wgpu::BindGroup; 2],
    bc_bg: [[wgpu::BindGroup; 2]; 6],
    clear_bg: wgpu::BindGroup,
    state: RefCell<RecState>,
    // Host-side copies for `reduce` (set by the upload path).
    pub(crate) host_solid: Vec<bool>,
    pub(crate) host_ff: Option<Vec<[f32; 3]>>,
    pub(crate) has_probe: bool,
    pub(crate) profile_set: [bool; 6],
}

impl GpuFields {
    fn cur(&self) -> usize {
        self.state.borrow().cur
    }

    /// Compact plane length (`nx * ny`).
    pub fn n_cells(&self) -> usize {
        self.n
    }
}

/// The wgpu implementation of [`Backend`] for a 2D lattice, `T = f32`
/// (WGSL has no f64; f32 deviation storage is the validated GPU grade —
/// GPU_EVALUATION.md §2).
pub struct WgpuBackend<L: Lattice> {
    ctx: Arc<GpuContext>,
    pipelines: Arc<Pipelines>,
    /// Steps per queue submit during batched runs (no waits in between).
    /// The proto measured 7.3–7.4 GLUPS for 10–100 dispatches/submit vs
    /// 0.8 GLUPS with a wait per step; anything ≥ ~10 is on the plateau.
    pub submit_chunk: usize,
    _l: PhantomData<L>,
}

impl<L: Lattice> WgpuBackend<L> {
    /// Compile the generated WGSL and build the pipeline set.
    pub fn new(ctx: Arc<GpuContext>) -> Self {
        assert_eq!(
            L::D,
            2,
            "WgpuBackend currently generates 2D (D2Q9) kernels; D3Q19 lands with M-C+"
        );
        let source = wgsl::generate::<L>();
        let module = ctx
            .device
            .create_shader_module(wgpu::ShaderModuleDescriptor {
                label: Some("lbm-core-d2q9"),
                source: wgpu::ShaderSource::Wgsl(source.into()),
            });
        let mk = |entry: &str| {
            ctx.device
                .create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
                    label: Some(entry),
                    layout: None, // auto: each entry point binds exactly what it uses
                    module: &module,
                    entry_point: Some(entry),
                    compilation_options: Default::default(),
                    cache: None,
                })
        };
        let pipelines = Arc::new(Pipelines {
            step: mk("step"),
            moments: mk("moments"),
            bc: mk("bc"),
            clear_probe: mk("clear_probe"),
        });
        Self {
            ctx,
            pipelines,
            submit_chunk: 200,
            _l: PhantomData,
        }
    }

    /// The shared device context.
    pub fn context(&self) -> &Arc<GpuContext> {
        &self.ctx
    }

    fn workgroups(&self, fields: &GpuFields) -> (u32, u32) {
        (
            fields.nx.div_ceil(wgsl::WG.0),
            fields.ny.div_ceil(wgsl::WG.1),
        )
    }

    fn bc_extent(&self, fields: &GpuFields, face: usize) -> u32 {
        if Face::ALL[face].axis() == 0 {
            fields.ny
        } else {
            fields.nx
        }
    }

    /// Materialise recorded ops into one encoder/pass and submit — without
    /// waiting (waits only happen in the readback paths).
    pub fn flush(&self, fields: &GpuFields) {
        let mut st = fields.state.borrow_mut();
        if st.ops.is_empty() {
            return;
        }
        let (gx, gy) = self.workgroups(fields);
        let mut enc = self
            .ctx
            .device
            .create_command_encoder(&wgpu::CommandEncoderDescriptor {
                label: Some("lbm-steps"),
            });
        {
            let mut pass = enc.begin_compute_pass(&wgpu::ComputePassDescriptor {
                label: Some("lbm-steps"),
                timestamp_writes: None,
            });
            for op in &st.ops {
                match *op {
                    Op::ClearProbe => {
                        pass.set_pipeline(&self.pipelines.clear_probe);
                        pass.set_bind_group(0, &fields.clear_bg, &[]);
                        pass.dispatch_workgroups(1, 1, 1);
                    }
                    Op::Fused { bg } => {
                        pass.set_pipeline(&self.pipelines.step);
                        pass.set_bind_group(0, &fields.fused_bg[bg], &[]);
                        pass.dispatch_workgroups(gx, gy, 1);
                    }
                    Op::Bc { face, bg } => {
                        pass.set_pipeline(&self.pipelines.bc);
                        pass.set_bind_group(0, &fields.bc_bg[face][bg], &[]);
                        let ext = self.bc_extent(fields, face);
                        pass.dispatch_workgroups(ext.div_ceil(wgsl::WG_BC), 1, 1);
                    }
                    Op::Moments { bg } => {
                        pass.set_pipeline(&self.pipelines.moments);
                        pass.set_bind_group(0, &fields.moments_bg[bg], &[]);
                        pass.dispatch_workgroups(gx, gy, 1);
                    }
                }
            }
        }
        self.ctx.queue.submit(Some(enc.finish()));
        st.ops.clear();
        st.steps_recorded = 0;
    }

    fn map_staging(&self, staging: &wgpu::Buffer, bytes: u64) -> Vec<u8> {
        let slice = staging.slice(..bytes);
        let (tx, rx) = futures_intrusive::channel::shared::oneshot_channel();
        slice.map_async(wgpu::MapMode::Read, move |res| {
            let _ = tx.send(res);
        });
        self.ctx.wait_idle();
        pollster::block_on(rx.receive())
            .expect("map_async callback dropped")
            .expect("staging buffer map failed");
        let out = slice.get_mapped_range().to_vec();
        staging.unmap();
        out
    }

    /// Explicit readback of the current deviation populations, compact SoA
    /// layout `f[q*n + y*nx + x]`. Cached until the next recorded step.
    pub fn read_f(&self, fields: &GpuFields) -> Vec<f32> {
        {
            let st = fields.state.borrow();
            if let Some((gen, data)) = &st.f_cache {
                if *gen == st.generation && st.ops.is_empty() {
                    return data.clone();
                }
            }
        }
        self.flush(fields);
        let bytes = (L::Q * fields.n * 4) as u64;
        let mut enc = self
            .ctx
            .device
            .create_command_encoder(&wgpu::CommandEncoderDescriptor {
                label: Some("read-f"),
            });
        enc.copy_buffer_to_buffer(&fields.f[fields.cur()], 0, &fields.staging, 0, bytes);
        self.ctx.queue.submit(Some(enc.finish()));
        let raw = self.map_staging(&fields.staging, bytes);
        let data: Vec<f32> = bytemuck::cast_slice(&raw).to_vec();
        let mut st = fields.state.borrow_mut();
        let generation = st.generation;
        st.f_cache = Some((generation, data.clone()));
        data
    }

    /// Explicit readback of the momentum-exchange probe force accumulated
    /// during the most recent executed step (V1 `probed_force` semantics).
    pub fn read_probed_force(&self, fields: &GpuFields) -> [f32; 3] {
        self.flush(fields);
        let mut enc = self
            .ctx
            .device
            .create_command_encoder(&wgpu::CommandEncoderDescriptor {
                label: Some("read-probe"),
            });
        enc.copy_buffer_to_buffer(&fields.probe_acc, 0, &fields.staging_small, 0, 12);
        self.ctx.queue.submit(Some(enc.finish()));
        let raw = self.map_staging(&fields.staging_small, 12);
        let v: &[f32] = bytemuck::cast_slice(&raw);
        [v[0], v[1], v[2]]
    }

    /// Write the global Params uniform once; assert the step parameters do
    /// not change between steps of a run (V1 semantics: relaxation rates,
    /// force and face BCs are fixed per solver lifetime).
    fn ensure_params(&self, sub: &Subdomain, fields: &GpuFields, p: &StepParams<f32>) {
        // Relaxation constants exactly as KParams::new builds them: f64
        // parameters converted to f32 once.
        let words: [u32; 12] = [
            fields.nx,
            fields.ny,
            (p.omega_p as f32).to_bits(),
            (p.omega_m as f32).to_bits(),
            ((1.0 - p.omega_p / 2.0) as f32).to_bits(),
            ((1.0 - p.omega_m / 2.0) as f32).to_bits(),
            p.force[0].to_bits(),
            p.force[1].to_bits(),
            {
                let halo = sub.halo_flags();
                let mut flags = 0u32;
                for (i, &h) in halo.iter().take(4).enumerate() {
                    if h {
                        flags |= wgsl::FLAG_HALO[i];
                    }
                }
                if fields.host_ff.is_some() {
                    flags |= wgsl::FLAG_FORCE_FIELD;
                }
                if fields.has_probe {
                    flags |= wgsl::FLAG_PROBE;
                }
                flags
            },
            0,
            0,
            0,
        ];
        let mut st = fields.state.borrow_mut();
        match &st.params_words {
            Some(prev) => assert_eq!(
                *prev, words,
                "step parameters changed mid-run; re-upload the fields first"
            ),
            None => {
                self.ctx
                    .queue
                    .write_buffer(&fields.params_ub, 0, bytemuck::cast_slice(&words));
                st.params_words = Some(words);
            }
        }
    }

    /// Build and upload the per-face BC uniforms from the step parameters
    /// and the lattice face tables (kernels.rs `zou_he_face` /
    /// `outflow_face` / `convective_face` constants).
    fn ensure_bc(&self, sub: &Subdomain, fields: &GpuFields, p: &StepParams<f32>) {
        {
            let st = fields.state.borrow();
            if let Some(prev) = &st.bc_words {
                let words = self.bc_words(sub, fields, p);
                assert_eq!(*prev, words, "face BCs changed mid-run");
                return;
            }
        }
        let words = self.bc_words(sub, fields, p);
        for (face, w) in words.iter().enumerate() {
            if w[0] != 0 {
                self.ctx
                    .queue
                    .write_buffer(&fields.bc_ub[face], 0, bytemuck::cast_slice(w));
            }
        }
        fields.state.borrow_mut().bc_words = Some(words);
    }

    fn bc_words(&self, sub: &Subdomain, fields: &GpuFields, p: &StepParams<f32>) -> [[u32; 32]; 6] {
        let (nx, ny) = (fields.nx, fields.ny);
        let mut out = [[0u32; 32]; 6];
        for face in Face::ALL {
            let fi = face.index();
            if face.axis() >= L::D || !sub.touches_global_face(face) {
                continue;
            }
            let bc = &p.faces[fi];
            if !bc.is_open() {
                continue;
            }
            let (base, stride, extent, joff): (u32, u32, u32, i32) = match face {
                Face::XNeg => (0, nx, ny, 1),
                Face::XPos => (nx - 1, nx, ny, -1),
                Face::YNeg => (0, 1, nx, nx as i32),
                Face::YPos => ((ny - 1) * nx, 1, nx, -(nx as i32)),
                _ => unreachable!("2D lattice"),
            };
            // Zou–He direction indices, derived exactly like zou_he_face.
            let n_in = face.n_in();
            let (nxi, nyi) = (n_in[0] as i32, n_in[1] as i32);
            let (tx, ty) = (-nyi, nxi);
            let q_n = L::dir_index([nxi as i8, nyi as i8, 0]);
            let q_d1 = L::dir_index([(nxi + tx) as i8, (nyi + ty) as i8, 0]);
            let q_d2 = L::dir_index([(nxi - tx) as i8, (nyi - ty) as i8, 0]);
            let q_t = L::dir_index([tx as i8, ty as i8, 0]);
            let q_mt = L::dir_index([-tx as i8, -ty as i8, 0]);
            let unk = L::unknowns(face);
            assert_eq!(unk.len(), 3, "D2Q9 face unknown count");
            let (kind, p0, p1) = match *bc {
                FaceBC::Closed => unreachable!(),
                FaceBC::Velocity { u } => (wgsl::BC_VELOCITY, u[0], u[1]),
                FaceBC::Pressure { rho } => (wgsl::BC_PRESSURE, rho, 0.0),
                FaceBC::Outflow => (wgsl::BC_OUTFLOW, 0.0, 0.0),
                FaceBC::Convective { u_conv } => (wgsl::BC_CONVECTIVE, u_conv, 0.0),
            };
            // Mass-pinning constants (convective_face): f64 weight sum, one
            // conversion; per-link weights in f32.
            let mut ws = 0.0f64;
            for &q in unk {
                ws += L::W[q];
            }
            let cinv = 1.0f32 / (1.0f32 + p0); // meaningful for Convective only
            let w = &mut out[fi];
            w[0] = kind;
            w[1] = base;
            w[2] = stride;
            w[3] = extent;
            w[4] = joff as u32;
            w[5] = u32::from(fields.profile_set[fi]);
            w[6] = q_n as u32;
            w[7] = L::OPP[q_n] as u32;
            w[8] = q_d1 as u32;
            w[9] = L::OPP[q_d1] as u32;
            w[10] = q_d2 as u32;
            w[11] = L::OPP[q_d2] as u32;
            w[12] = q_t as u32;
            w[13] = q_mt as u32;
            w[14] = unk[0] as u32;
            w[15] = unk[1] as u32;
            w[16] = unk[2] as u32;
            w[17] = p0.to_bits();
            w[18] = p1.to_bits();
            w[19] = (nxi as f32).to_bits();
            w[20] = (nyi as f32).to_bits();
            w[21] = (tx as f32).to_bits();
            w[22] = (ty as f32).to_bits();
            w[23] = (L::W[unk[0]] as f32).to_bits();
            w[24] = (L::W[unk[1]] as f32).to_bits();
            w[25] = (L::W[unk[2]] as f32).to_bits();
            w[26] = (ws as f32).to_bits();
            w[27] = cinv.to_bits();
        }
        out
    }

    /// Upload host-staged fields (populations, masks, moments, profiles)
    /// into the device buffers. `fields_host` is the monolithic subdomain's
    /// SoA state; halo rings are stripped (the kernel wraps in-place).
    ///
    /// Also primes the edge stash from `ftmp` so the ConvectiveOutflow
    /// previous-value convention continues the host state exactly (both
    /// start all-zero on a fresh solver).
    pub fn upload(&self, sub: &Subdomain, fields: &mut GpuFields, host: &SoaFields<f32>) {
        let g = host.geom;
        assert_eq!(g.core[0] as u32, fields.nx);
        assert_eq!(g.core[1] as u32, fields.ny);
        let (nx, ny, n) = (g.core[0], g.core[1], fields.n);
        let np = g.n_padded();
        let q = self.ctx.queue.clone();
        {
            let st = fields.state.borrow();
            assert!(
                st.ops.is_empty() && !st.pending_collide,
                "upload with recorded but unsubmitted steps"
            );
        }
        // Populations: current -> f[cur], ping-pong partner -> f[1-cur].
        let cur = fields.cur();
        let mut buf = vec![0f32; L::Q * n];
        for (src, dst) in [(&host.f, cur), (&host.ftmp, 1 - cur)] {
            for qi in 0..L::Q {
                for y in 0..ny {
                    for x in 0..nx {
                        buf[qi * n + y * nx + x] = src[qi * np + g.pidx(x, y, 0)];
                    }
                }
            }
            q.write_buffer(&fields.f[dst], 0, bytemuck::cast_slice(&buf));
        }
        // Edge stash (stash_in of the next fused dispatch = stash[cur]):
        // the ftmp values V1 would leave in the skipped slots.
        let slen = wgsl::stash_len::<L>(nx, ny);
        let mut stash = vec![0f32; slen];
        let mut off = 0usize;
        for face in &Face::ALL[..4] {
            let unk = L::unknowns(*face);
            let ext = if face.axis() == 0 { ny } else { nx };
            if !sub.has_halo(*face) {
                for (k, &u) in unk.iter().enumerate() {
                    for t in 0..ext {
                        let (x, y) = match face {
                            Face::XNeg => (0, t),
                            Face::XPos => (nx - 1, t),
                            Face::YNeg => (t, 0),
                            Face::YPos => (t, ny - 1),
                            _ => unreachable!(),
                        };
                        stash[off + k * ext + t] = host.ftmp[u * np + g.pidx(x, y, 0)];
                    }
                }
            }
            off += unk.len() * ext;
        }
        q.write_buffer(&fields.stash[cur], 0, bytemuck::cast_slice(&stash));
        q.write_buffer(
            &fields.stash[1 - cur],
            0,
            bytemuck::cast_slice(&vec![0f32; slen]),
        );
        // Mask (bit0 solid, bit1 probe) + host copies for reduce().
        let mut mask = vec![0u32; n];
        let mut host_solid = vec![false; n];
        for y in 0..ny {
            for x in 0..nx {
                let pi = g.pidx(x, y, 0);
                let c = y * nx + x;
                if host.solid[pi] {
                    mask[c] |= 1;
                    host_solid[c] = true;
                }
                if let Some(pm) = &host.probe {
                    if pm[pi] {
                        mask[c] |= 2;
                    }
                }
            }
        }
        q.write_buffer(&fields.mask, 0, bytemuck::cast_slice(&mask));
        fields.host_solid = host_solid;
        fields.has_probe = host.probe.is_some();
        // Wall velocities (read only at solid neighbours).
        let mut wu = vec![0f32; 2 * n];
        for y in 0..ny {
            for x in 0..nx {
                let v = host.wall_u[g.pidx(x, y, 0)];
                wu[2 * (y * nx + x)] = v[0];
                wu[2 * (y * nx + x) + 1] = v[1];
            }
        }
        q.write_buffer(&fields.wall_u, 0, bytemuck::cast_slice(&wu));
        // Per-cell force field (compact already).
        fields.host_ff = host.force_field.clone();
        if let Some(ff) = &host.force_field {
            let mut fv = vec![0f32; 2 * n];
            for (c, v) in ff.iter().enumerate() {
                fv[2 * c] = v[0];
                fv[2 * c + 1] = v[1];
            }
            q.write_buffer(&fields.force_field, 0, bytemuck::cast_slice(&fv));
        }
        // Moments (compact already; carries V1's values at solid cells,
        // which the moments kernel never rewrites).
        q.write_buffer(&fields.rho, 0, bytemuck::cast_slice(&host.rho));
        q.write_buffer(&fields.ux, 0, bytemuck::cast_slice(&host.ux));
        q.write_buffer(&fields.uy, 0, bytemuck::cast_slice(&host.uy));
        // Inlet profiles.
        for face in Face::ALL {
            let fi = face.index();
            fields.profile_set[fi] = false;
            if let Some(prof) = &host.inlet_profiles[fi] {
                let mut pv = vec![0f32; 2 * prof.len()];
                for (t, u) in prof.iter().enumerate() {
                    pv[2 * t] = u[0];
                    pv[2 * t + 1] = u[1];
                }
                q.write_buffer(&fields.profiles[fi], 0, bytemuck::cast_slice(&pv));
                fields.profile_set[fi] = true;
            }
        }
        // Probe accumulator and cached uniforms reset (masks/probe/profile
        // presence may have changed the flags).
        q.write_buffer(&fields.probe_acc, 0, &[0u8; 12]);
        let mut st = fields.state.borrow_mut();
        st.params_words = None;
        st.bc_words = None;
        st.f_cache = None;
        st.generation += 1;
    }
}

impl<L: Lattice> Backend<L, f32> for WgpuBackend<L> {
    type Fields = GpuFields;

    fn alloc(&self, sub: &Subdomain) -> GpuFields {
        let g = sub.geom;
        assert_eq!(g.d, 2, "WgpuBackend fields are 2D (D2Q9)");
        assert_eq!(g.core[2], 1);
        let (nx, ny) = (g.core[0] as u32, g.core[1] as u32);
        let n = (nx as usize) * (ny as usize);
        let device = &self.ctx.device;
        use wgpu::BufferUsages as U;
        let buf = |label: &str, size: u64, usage: wgpu::BufferUsages| {
            device.create_buffer(&wgpu::BufferDescriptor {
                label: Some(label),
                size,
                usage,
                mapped_at_creation: false,
            })
        };
        let fbytes = (L::Q * n * 4) as u64;
        let f = [
            buf("f0", fbytes, U::STORAGE | U::COPY_DST | U::COPY_SRC),
            buf("f1", fbytes, U::STORAGE | U::COPY_DST | U::COPY_SRC),
        ];
        let slen = (wgsl::stash_len::<L>(nx as usize, ny as usize) * 4) as u64;
        let stash = [
            buf("stash0", slen, U::STORAGE | U::COPY_DST),
            buf("stash1", slen, U::STORAGE | U::COPY_DST),
        ];
        let mask = buf("mask", (n * 4) as u64, U::STORAGE | U::COPY_DST);
        let wall_u = buf("wall_u", (n * 8) as u64, U::STORAGE | U::COPY_DST);
        let force_field = buf("force_field", (n * 8) as u64, U::STORAGE | U::COPY_DST);
        let rho = buf("rho", (n * 4) as u64, U::STORAGE | U::COPY_DST | U::COPY_SRC);
        let ux = buf("ux", (n * 4) as u64, U::STORAGE | U::COPY_DST | U::COPY_SRC);
        let uy = buf("uy", (n * 4) as u64, U::STORAGE | U::COPY_DST | U::COPY_SRC);
        let probe_acc = buf("probe_acc", 12, U::STORAGE | U::COPY_DST | U::COPY_SRC);
        let params_ub = buf("params", 48, U::UNIFORM | U::COPY_DST);
        let bc_ub = std::array::from_fn(|i| {
            buf(&format!("bc{i}"), 128, U::UNIFORM | U::COPY_DST)
        });
        let profiles = std::array::from_fn(|i| {
            let ext = if i < 2 { ny } else { nx } as u64;
            buf(&format!("profile{i}"), (ext * 8).max(8), U::STORAGE | U::COPY_DST)
        });
        let staging = buf("staging", fbytes, U::MAP_READ | U::COPY_DST);
        let staging_small = buf("staging-small", 12, U::MAP_READ | U::COPY_DST);

        fn e(binding: u32, b: &wgpu::Buffer) -> wgpu::BindGroupEntry<'_> {
            wgpu::BindGroupEntry {
                binding,
                resource: b.as_entire_binding(),
            }
        }
        let step_layout = self.pipelines.step.get_bind_group_layout(0);
        let fused_bg = [0usize, 1].map(|p| {
            device.create_bind_group(&wgpu::BindGroupDescriptor {
                label: Some("fused"),
                layout: &step_layout,
                entries: &[
                    e(0, &params_ub),
                    e(1, &f[p]),
                    e(2, &f[1 - p]),
                    e(3, &mask),
                    e(4, &wall_u),
                    e(5, &force_field),
                    e(6, &stash[p]),
                    e(7, &stash[1 - p]),
                    e(8, &probe_acc),
                ],
            })
        });
        let moments_layout = self.pipelines.moments.get_bind_group_layout(0);
        let moments_bg = [0usize, 1].map(|p| {
            device.create_bind_group(&wgpu::BindGroupDescriptor {
                label: Some("moments"),
                layout: &moments_layout,
                entries: &[
                    e(0, &params_ub),
                    e(1, &f[p]),
                    e(3, &mask),
                    e(5, &force_field),
                    e(9, &rho),
                    e(10, &ux),
                    e(11, &uy),
                ],
            })
        });
        let bc_layout = self.pipelines.bc.get_bind_group_layout(0);
        let bc_bg = std::array::from_fn(|face| {
            [0usize, 1].map(|p| {
                device.create_bind_group(&wgpu::BindGroupDescriptor {
                    label: Some("bc"),
                    layout: &bc_layout,
                    entries: &[
                        e(0, &params_ub),
                        e(2, &f[p]),
                        e(3, &mask),
                        e(12, &bc_ub[face]),
                        e(13, &profiles[face]),
                    ],
                })
            })
        });
        let clear_layout = self.pipelines.clear_probe.get_bind_group_layout(0);
        let clear_bg = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("clear-probe"),
            layout: &clear_layout,
            entries: &[e(8, &probe_acc)],
        });

        GpuFields {
            nx,
            ny,
            n,
            f,
            stash,
            mask,
            wall_u,
            force_field,
            rho,
            ux,
            uy,
            probe_acc,
            params_ub,
            bc_ub,
            profiles,
            staging,
            staging_small,
            fused_bg,
            moments_bg,
            bc_bg,
            clear_bg,
            state: RefCell::new(RecState {
                cur: 0,
                ops: Vec::new(),
                steps_recorded: 0,
                pending_collide: false,
                params_words: None,
                bc_words: None,
                generation: 0,
                f_cache: None,
            }),
            host_solid: vec![false; n],
            host_ff: None,
            has_probe: false,
            profile_set: [false; 6],
        }
    }

    fn collide(&mut self, sub: &Subdomain, fields: &mut GpuFields, p: &StepParams<f32>) {
        self.ensure_params(sub, fields, p);
        let mut st = fields.state.borrow_mut();
        assert!(!st.pending_collide, "collide called twice without stream");
        st.pending_collide = true;
    }

    fn stream(
        &mut self,
        sub: &Subdomain,
        fields: &mut GpuFields,
        _p: &StepParams<f32>,
        range: CellRange,
    ) -> [f32; 3] {
        assert_eq!(
            range,
            CellRange::full(sub),
            "WgpuBackend streams the full grid in one fused dispatch (no two-pass split)"
        );
        let mut st = fields.state.borrow_mut();
        assert!(
            st.pending_collide,
            "stream without a preceding collide (the fused kernel does both)"
        );
        st.pending_collide = false;
        if fields.has_probe {
            st.ops.push(Op::ClearProbe);
        }
        let cur = st.cur;
        st.ops.push(Op::Fused { bg: cur });
        st.generation += 1;
        st.f_cache = None;
        // The probe force accumulates on-device; explicit readback via
        // read_probed_force (per-step CPU sync is the 9x trap).
        [0.0; 3]
    }

    fn swap(&mut self, fields: &mut GpuFields) {
        let mut st = fields.state.borrow_mut();
        st.cur ^= 1;
    }

    fn apply_open_faces(&mut self, sub: &Subdomain, fields: &mut GpuFields, p: &StepParams<f32>) {
        self.ensure_bc(sub, fields, p);
        let mut st = fields.state.borrow_mut();
        let cur = st.cur;
        for face in Face::ALL {
            if face.axis() >= L::D || !sub.touches_global_face(face) {
                continue;
            }
            if p.faces[face.index()].is_open() {
                st.ops.push(Op::Bc {
                    face: face.index(),
                    bg: cur,
                });
            }
        }
    }

    fn update_moments(&mut self, _sub: &Subdomain, fields: &mut GpuFields, _p: &StepParams<f32>) {
        // Lazy: the fused kernel re-derives (rho, u) from the identical
        // pre-collide state, so the device moment buffers are refreshed only
        // by read_moments. This call marks the end of a step — the chunked
        // submit hook.
        let flush_now = {
            let mut st = fields.state.borrow_mut();
            st.steps_recorded += 1;
            st.steps_recorded >= self.submit_chunk
        };
        if flush_now {
            self.flush(fields);
        }
    }

    fn reduce(
        &self,
        _sub: &Subdomain,
        fields: &GpuFields,
        p: &StepParams<f32>,
        kind: Reduction,
    ) -> f64 {
        // Host-side exact V1 loop over the read-back populations: compact
        // cell order (y, x ascending), q inner, f64 accumulation — the same
        // sequence as CpuScalar::reduce on a monolithic subdomain.
        let f = self.read_f(fields);
        let (nx, ny, n) = (fields.nx as usize, fields.ny as usize, fields.n);
        let mut acc = 0.0f64;
        for y in 0..ny {
            for x in 0..nx {
                let c = y * nx + x;
                if fields.host_solid[c] {
                    continue;
                }
                match kind {
                    Reduction::FluidCells => acc += 1.0,
                    Reduction::MassDeviation => {
                        for q in 0..L::Q {
                            acc += f[q * n + c] as f64;
                        }
                    }
                    Reduction::Momentum(a) => {
                        let mut m = 0.0f64;
                        for q in 0..L::Q {
                            m += L::C[q][a] as f64 * f[q * n + c] as f64;
                        }
                        let fa = match &fields.host_ff {
                            Some(field) => p.force[a] as f64 + field[c][a] as f64,
                            None => p.force[a] as f64,
                        };
                        acc += m + 0.5 * fa;
                    }
                }
            }
        }
        acc
    }

    fn read_moments(&self, fields: &GpuFields, out: &mut HostMoments<f32>) {
        {
            let mut st = fields.state.borrow_mut();
            let cur = st.cur;
            st.ops.push(Op::Moments { bg: cur });
        }
        self.flush(fields);
        let n = fields.n;
        let plane = (n * 4) as u64;
        let mut enc = self
            .ctx
            .device
            .create_command_encoder(&wgpu::CommandEncoderDescriptor {
                label: Some("read-moments"),
            });
        enc.copy_buffer_to_buffer(&fields.rho, 0, &fields.staging, 0, plane);
        enc.copy_buffer_to_buffer(&fields.ux, 0, &fields.staging, plane, plane);
        enc.copy_buffer_to_buffer(&fields.uy, 0, &fields.staging, 2 * plane, plane);
        self.ctx.queue.submit(Some(enc.finish()));
        let raw = self.map_staging(&fields.staging, 3 * plane);
        let v: &[f32] = bytemuck::cast_slice(&raw);
        out.rho.clear();
        out.rho.extend_from_slice(&v[..n]);
        out.ux.clear();
        out.ux.extend_from_slice(&v[n..2 * n]);
        out.uy.clear();
        out.uy.extend_from_slice(&v[2 * n..3 * n]);
        out.uz.clear();
        out.uz.resize(n, 0.0);
    }
}
