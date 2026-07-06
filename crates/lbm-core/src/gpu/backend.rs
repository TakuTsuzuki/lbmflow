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
//!   The probe force accumulates on-device and is fetched through
//!   [`WgpuBackend::read_probed_force`] (explicit readback), holding the most
//!   recent step's value like V1 `probed_force`.
//! - `swap` — flips the ping-pong index (bind groups pre-built per parity).
//! - `apply_open_faces` — records the per-face `bc` dispatches in
//!   `Face::ALL` order (CPU order).
//! - `update_moments` — lazy: the fused kernel re-derives the moments it
//!   needs in-kernel from the same pre-collide state the CPU caches, so the
//!   device moment buffers are only refreshed by `read_moments`.
//! - `end_step` — records the step boundary for chunked submission.
//! - `reduce` — explicit readback of the populations plus the exact V1
//!   f64 accumulation loop on the host (wgpu has no f64; doing the loop
//!   host-side keeps the diagnostic convention bit-compatible with
//!   `CpuScalar::reduce` for identical fields).
//!
//! Field storage is compact (`cell = y*nx + x`, no halo ring): periodic
//! wrap happens in-kernel, which is the monolithic-subdomain configuration.
//! Multi-part GPU decompositions (device-side halo exchange) are out of
//! scope here (M-D).

use std::cell::{Cell, RefCell};
use std::marker::PhantomData;
use std::sync::{Arc, Mutex};
use std::time::Duration;

use crate::backend::{write_host_moments, Backend, CellRange, HostMoments};
use crate::fields::SoaFields;
use crate::halo::HaloExchange;
use crate::lattice::{Face, Lattice};
use crate::params::{CollisionKind, FaceBC, Reduction, StepParams};
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
    device_lost: Arc<Mutex<Option<String>>>,
}

/// Runtime GPU failure reported by wgpu operations.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum GpuError {
    /// Device polling failed before queued work completed.
    Poll(String),
    /// Staging-buffer mapping failed or did not invoke its callback.
    Map(String),
    /// wgpu reported that the device was lost.
    DeviceLost(String),
    /// Requested grid exceeds the selected adapter/device resource limits.
    ResourceLimit(String),
    /// The spec was rejected before any device work was submitted (features
    /// the GPU backend does not implement — honest failure instead of silent
    /// wrong physics). Carries the `SpecError` message.
    Spec(String),
}

impl std::fmt::Display for GpuError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Poll(msg) => write!(f, "GPU device poll failed: {msg}"),
            Self::Map(msg) => write!(f, "GPU staging buffer map failed: {msg}"),
            Self::DeviceLost(msg) => write!(f, "GPU device was lost: {msg}"),
            Self::ResourceLimit(msg) => write!(f, "GPU resource limit exceeded: {msg}"),
            Self::Spec(msg) => write!(f, "GPU solver rejected spec: {msg}"),
        }
    }
}

impl std::error::Error for GpuError {}

/// GPU context initialization failure with enough adapter detail for fallback logs.
#[derive(Clone, Debug)]
pub enum GpuInitError {
    /// No adapter matched the requested power preference.
    NoAdapter,
    /// Adapter was found but device creation failed.
    RequestDevice {
        /// Adapter selected before the request failed.
        adapter_info: wgpu::AdapterInfo,
        /// wgpu's device-request failure text.
        message: String,
    },
    /// The selected adapter does not expose a required feature.
    MissingFeature {
        /// Adapter selected before the feature check failed.
        adapter_info: wgpu::AdapterInfo,
        /// Feature that was required.
        feature: &'static str,
    },
}

impl std::fmt::Display for GpuInitError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::NoAdapter => write!(f, "no usable GPU adapter was found"),
            Self::MissingFeature {
                adapter_info,
                feature,
            } => write!(
                f,
                "GPU adapter '{}' ({:?}) does not support required feature {feature}",
                adapter_info.name, adapter_info.backend
            ),
            Self::RequestDevice {
                adapter_info,
                message,
            } => write!(
                f,
                "failed to create GPU device for adapter '{}' ({:?}): {message}",
                adapter_info.name, adapter_info.backend
            ),
        }
    }
}

/// Distribution-buffer storage precision requested for generated kernels.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum GpuStorage {
    /// f32 distribution buffers.
    F32,
    /// f16 distribution buffers; shader arithmetic remains f32.
    F16,
}

impl GpuStorage {
    fn element_bytes(self) -> u64 {
        match self {
            Self::F32 => 4,
            Self::F16 => 2,
        }
    }
}

fn f32_to_f16_bits(value: f32) -> u16 {
    let x = value.to_bits();
    let sign = ((x >> 16) & 0x8000) as u16;
    let exp = ((x >> 23) & 0xff) as i32 - 127 + 15;
    let mant = x & 0x7f_ffff;
    if exp <= 0 {
        if exp < -10 {
            return sign;
        }
        let mant = mant | 0x80_0000;
        let shift = (14 - exp) as u32;
        let mut half = (mant >> shift) as u16;
        if ((mant >> (shift - 1)) & 1) != 0 {
            half = half.wrapping_add(1);
        }
        return sign | half;
    }
    if exp >= 31 {
        return sign | 0x7c00;
    }
    let mut half = sign | ((exp as u16) << 10) | ((mant >> 13) as u16);
    if (mant & 0x1000) != 0 {
        half = half.wrapping_add(1);
    }
    half
}

fn f16_bits_to_f32(bits: u16) -> f32 {
    let sign = ((bits & 0x8000) as u32) << 16;
    let exp = (bits >> 10) & 0x1f;
    let mant = (bits & 0x03ff) as u32;
    let out = if exp == 0 {
        if mant == 0 {
            sign
        } else {
            let mut mantissa = mant;
            let mut e = -14i32;
            while (mantissa & 0x400) == 0 {
                mantissa <<= 1;
                e -= 1;
            }
            mantissa &= 0x3ff;
            sign | (((e + 127) as u32) << 23) | (mantissa << 13)
        }
    } else if exp == 31 {
        sign | 0x7f80_0000 | (mant << 13)
    } else {
        sign | ((((exp as i32 - 15 + 127) as u32) << 23) | (mant << 13))
    };
    f32::from_bits(out)
}

fn encode_storage(storage: GpuStorage, values: &[f32]) -> Vec<u8> {
    match storage {
        GpuStorage::F32 => bytemuck::cast_slice(values).to_vec(),
        GpuStorage::F16 => values
            .iter()
            .flat_map(|&v| f32_to_f16_bits(v).to_le_bytes())
            .collect(),
    }
}

fn decode_storage(storage: GpuStorage, bytes: &[u8]) -> Vec<f32> {
    match storage {
        GpuStorage::F32 => bytemuck::cast_slice(bytes).to_vec(),
        GpuStorage::F16 => bytes
            .chunks_exact(2)
            .map(|b| f16_bits_to_f32(u16::from_le_bytes([b[0], b[1]])))
            .collect(),
    }
}

/// GPU kernel configuration.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct KernelCfg {
    /// Distribution-buffer storage format.
    pub storage: GpuStorage,
}

impl Default for KernelCfg {
    fn default() -> Self {
        Self {
            storage: GpuStorage::F32,
        }
    }
}

impl std::error::Error for GpuInitError {}

/// Pick the next submit size from one measured submit.
pub(crate) fn calibrate_submit_chunk(
    measured_steps: usize,
    elapsed: Duration,
    current: usize,
) -> usize {
    const MAX_CHUNK: usize = 200;
    const TARGET_MS: f64 = 175.0;
    if measured_steps == 0 {
        return current.clamp(1, MAX_CHUNK);
    }
    let elapsed_ms = elapsed.as_secs_f64() * 1000.0;
    if !(elapsed_ms.is_finite() && elapsed_ms > 0.0) {
        return current.clamp(1, MAX_CHUNK);
    }
    ((measured_steps as f64 * TARGET_MS / elapsed_ms).round() as usize).clamp(1, MAX_CHUNK)
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct GpuResourcePlan {
    n: usize,
    element_bytes: u64,
    fbytes: u64,
    stash_bytes: u64,
    mask_bytes: u64,
    vec2_bytes: u64,
    moments_bytes: u64,
    staging_bytes: u64,
    gx: u32,
    gy: u32,
    gz: u32,
    max_bc_groups: u32,
}

impl GpuResourcePlan {
    fn for_grid<L: Lattice>(
        nx: usize,
        ny: usize,
        nz: usize,
        element_bytes: u64,
    ) -> Result<Self, GpuError> {
        let n = nx
            .checked_mul(ny)
            .and_then(|xy| xy.checked_mul(nz))
            .ok_or_else(|| GpuError::ResourceLimit(format!("grid {nx}x{ny} overflows usize")))?;
        let qn = L::Q.checked_mul(n).ok_or_else(|| {
            GpuError::ResourceLimit(format!(
                "Q*n overflows usize for D{}Q{} on {nx}x{ny}",
                L::D,
                L::Q
            ))
        })?;
        if qn > u32::MAX as usize {
            return Err(GpuError::ResourceLimit(format!(
                "Q*n = {qn} exceeds u32::MAX; split the domain or reduce grid size"
            )));
        }
        let bytes = |items: usize, item_bytes: u64, what: &str| -> Result<u64, GpuError> {
            (items as u64)
                .checked_mul(item_bytes)
                .ok_or_else(|| GpuError::ResourceLimit(format!("{what} byte count overflows u64")))
        };
        let nx_u32 = u32::try_from(nx).map_err(|_| {
            GpuError::ResourceLimit(format!("nx = {nx} exceeds u32::MAX for WGSL indexing"))
        })?;
        let ny_u32 = u32::try_from(ny).map_err(|_| {
            GpuError::ResourceLimit(format!("ny = {ny} exceeds u32::MAX for WGSL indexing"))
        })?;
        let nz_u32 = u32::try_from(nz).map_err(|_| {
            GpuError::ResourceLimit(format!("nz = {nz} exceeds u32::MAX for WGSL indexing"))
        })?;
        let max_ext = nx_u32.max(ny_u32).max(nz_u32);
        Ok(Self {
            n,
            element_bytes,
            fbytes: bytes(qn, element_bytes, "population buffer")?,
            stash_bytes: bytes(
                wgsl::stash_len::<L>(nx, ny, nz),
                element_bytes,
                "edge stash",
            )?,
            mask_bytes: bytes(n, 4, "solid mask")?,
            vec2_bytes: bytes(n, 16, "vec3 field")?,
            moments_bytes: bytes(n, 4, "moment buffer")?,
            staging_bytes: bytes(qn, 4, "population staging")?,
            gx: nx_u32.div_ceil(wgsl::WG.0),
            gy: ny_u32.div_ceil(wgsl::WG.1),
            gz: nz_u32,
            max_bc_groups: max_ext.div_ceil(wgsl::WG_BC),
        })
    }

    fn validate(&self, limits: &wgpu::Limits) -> Result<(), GpuError> {
        let check_buffer = |label: &str, bytes: u64| -> Result<(), GpuError> {
            if bytes > limits.max_buffer_size {
                return Err(GpuError::ResourceLimit(format!(
                    "{label} requires {bytes} bytes, above max_buffer_size {}",
                    limits.max_buffer_size
                )));
            }
            Ok(())
        };
        let check_storage = |label: &str, bytes: u64| -> Result<(), GpuError> {
            check_buffer(label, bytes)?;
            if bytes > limits.max_storage_buffer_binding_size as u64 {
                return Err(GpuError::ResourceLimit(format!(
                    "{label} requires {bytes} bytes, above max_storage_buffer_binding_size {}",
                    limits.max_storage_buffer_binding_size
                )));
            }
            Ok(())
        };
        for (label, bytes) in [
            ("population f0/f1", self.fbytes),
            ("edge stash", self.stash_bytes),
            ("solid mask", self.mask_bytes),
            ("wall velocity", self.vec2_bytes),
            ("force field", self.vec2_bytes),
            ("rho moment", self.moments_bytes),
            ("ux moment", self.moments_bytes),
            ("uy moment", self.moments_bytes),
            ("uz moment", self.moments_bytes),
        ] {
            check_storage(label, bytes)?;
        }
        check_buffer("population staging", self.staging_bytes)?;
        let max_groups = limits.max_compute_workgroups_per_dimension.min(65_535);
        for (axis, groups) in [
            ("x", self.gx),
            ("y", self.gy),
            ("z", self.gz),
            ("boundary", self.max_bc_groups),
        ] {
            if groups > max_groups {
                return Err(GpuError::ResourceLimit(format!(
                    "dispatch {axis} workgroup count {groups} exceeds limit {max_groups}"
                )));
            }
        }
        Ok(())
    }
}

impl GpuContext {
    /// Create a context on the highest-performance adapter.
    pub fn new() -> Result<Arc<Self>, GpuInitError> {
        Self::new_with_shader_f16(false)
    }

    /// Create a context, optionally requiring `SHADER_F16`.
    pub fn new_with_shader_f16(require_shader_f16: bool) -> Result<Arc<Self>, GpuInitError> {
        let instance = wgpu::Instance::default();
        let mut adapter = None;
        for power_preference in [
            wgpu::PowerPreference::HighPerformance,
            wgpu::PowerPreference::LowPower,
            wgpu::PowerPreference::None,
        ] {
            if let Ok(found) =
                pollster::block_on(instance.request_adapter(&wgpu::RequestAdapterOptions {
                    power_preference,
                    ..Default::default()
                }))
            {
                adapter = Some(found);
                break;
            }
        }
        let adapter = adapter.ok_or(GpuInitError::NoAdapter)?;
        let adapter_info = adapter.get_info();
        let adapter_features = adapter.features();
        if require_shader_f16 && !adapter_features.contains(wgpu::Features::SHADER_F16) {
            return Err(GpuInitError::MissingFeature {
                adapter_info,
                feature: "SHADER_F16",
            });
        }
        let al = adapter.limits();
        let mut limits = wgpu::Limits::default();
        limits.max_storage_buffer_binding_size = al.max_storage_buffer_binding_size;
        limits.max_buffer_size = al.max_buffer_size;
        limits.max_storage_buffers_per_shader_stage = al.max_storage_buffers_per_shader_stage;
        let required_features = if require_shader_f16 {
            wgpu::Features::SHADER_F16
        } else {
            wgpu::Features::empty()
        };
        let (device, queue) = pollster::block_on(adapter.request_device(&wgpu::DeviceDescriptor {
            label: Some("lbm-core-gpu"),
            required_features,
            required_limits: limits,
            memory_hints: wgpu::MemoryHints::Performance,
            trace: wgpu::Trace::Off,
        }))
        .map_err(|e| GpuInitError::RequestDevice {
            adapter_info: adapter_info.clone(),
            message: e.to_string(),
        })?;
        let device_lost = Arc::new(Mutex::new(None));
        let lost_slot = Arc::clone(&device_lost);
        device.set_device_lost_callback(move |reason, message| {
            let text = if message.is_empty() {
                format!("{reason:?}")
            } else {
                format!("{reason:?}: {message}")
            };
            *lost_slot.lock().expect("device lost mutex poisoned") = Some(text);
        });
        Ok(Arc::new(Self {
            device,
            queue,
            adapter_info,
            device_lost,
        }))
    }

    /// Block until all submitted work completed.
    pub fn wait_idle(&self) {
        self.try_wait_idle().expect("device poll failed");
    }

    /// Block until all submitted work completed, propagating poll/device loss.
    pub fn try_wait_idle(&self) -> Result<(), GpuError> {
        self.device
            .poll(wgpu::PollType::Wait)
            .map(|_| ())
            .map_err(|e| GpuError::Poll(e.to_string()))?;
        self.check_device_lost()
    }

    fn check_device_lost(&self) -> Result<(), GpuError> {
        match self
            .device_lost
            .lock()
            .expect("device lost mutex poisoned")
            .clone()
        {
            Some(reason) => Err(GpuError::DeviceLost(reason)),
            None => Ok(()),
        }
    }
}

struct Pipelines {
    step: wgpu::ComputePipeline,
    step_cached: wgpu::ComputePipeline,
    step_wale: wgpu::ComputePipeline,
    step_cached_wale: wgpu::ComputePipeline,
    step_cumulant: wgpu::ComputePipeline,
    step_cached_cumulant: wgpu::ComputePipeline,
    step_wale_cumulant: wgpu::ComputePipeline,
    step_cached_wale_cumulant: wgpu::ComputePipeline,
    wale_omega: wgpu::ComputePipeline,
    moments: wgpu::ComputePipeline,
    bc: wgpu::ComputePipeline,
    clear_probe: wgpu::ComputePipeline,
}

/// Recorded (not yet submitted) dispatch.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Op {
    ClearProbe,
    /// Fused collide+stream, reading buffer parity `bg`.
    Fused {
        bg: usize,
        cached: bool,
        wale: bool,
        cumulant: bool,
    },
    /// WALE omega refresh from pre-step moments and populations.
    WaleOmega {
        bg: usize,
    },
    /// Open-face BC on `face`, operating on buffer parity `bg`.
    Bc {
        face: usize,
        bg: usize,
    },
    /// Moments refresh from buffer parity `bg`.
    Moments {
        bg: usize,
    },
}

/// Interior-mutable recorder state (readback methods take `&Fields`).
struct RecState {
    cur: usize,
    ops: Vec<Op>,
    steps_recorded: usize,
    pending_collide: bool,
    use_cached_moments_once: bool,
    /// Written Params uniform (asserts step-parameter stability per run).
    params_words: Option<[u32; 20]>,
    /// Written per-face BC uniforms.
    bc_words: Option<[[u32; 64]; 6]>,
    /// Bumped per fused dispatch; invalidates the readback cache.
    generation: u64,
    /// Cached population readback (generation, shared data).
    f_cache: Option<(u64, Arc<Vec<f32>>)>,
}

struct StagingBuffer {
    buffer: wgpu::Buffer,
    size: u64,
}

struct WaleBuffers {
    omega: wgpu::Buffer,
    fused_bg: [wgpu::BindGroup; 2],
    fused_cached_bg: [wgpu::BindGroup; 2],
    wale_bg: [wgpu::BindGroup; 2],
}

/// Device-resident fields of one (monolithic) subdomain.
pub struct GpuFields {
    nx: u32,
    ny: u32,
    nz: u32,
    n: usize,
    element_bytes: u64,
    f: [wgpu::Buffer; 2],
    stash: [wgpu::Buffer; 2],
    mask: wgpu::Buffer,
    wall_u: wgpu::Buffer,
    force_field: wgpu::Buffer,
    wall_u_full: bool,
    force_field_full: bool,
    rho: wgpu::Buffer,
    ux: wgpu::Buffer,
    uy: wgpu::Buffer,
    uz: wgpu::Buffer,
    probe_acc: wgpu::Buffer,
    params_ub: wgpu::Buffer,
    bc_ub: [wgpu::Buffer; 6],
    profiles: [wgpu::Buffer; 6],
    staging: RefCell<Option<StagingBuffer>>,
    fused_bg: [wgpu::BindGroup; 2],
    fused_cached_bg: [wgpu::BindGroup; 2],
    moments_bg: [wgpu::BindGroup; 2],
    bc_bg: [[wgpu::BindGroup; 2]; 6],
    clear_bg: wgpu::BindGroup,
    wale: Option<WaleBuffers>,
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

fn has_open_faces<L: Lattice>(sub: &Subdomain, p: &StepParams<f32>) -> bool {
    Face::ALL.iter().any(|&face| {
        face.axis() < L::D && sub.touches_global_face(face) && p.faces[face.index()].is_open()
    })
}

/// The wgpu implementation of [`Backend`] for a 2D lattice, `T = f32`
/// (WGSL has no f64; f32 deviation storage is the validated GPU grade —
/// GPU_EVALUATION.md §2).
pub struct WgpuBackend<L: Lattice> {
    ctx: Arc<GpuContext>,
    cfg: KernelCfg,
    pipelines: Arc<Pipelines>,
    submissions: Cell<u64>,
    cache_next_upload_moments_once: Cell<bool>,
    /// Steps per queue submit during batched runs (no waits in between).
    /// The proto measured 7.3–7.4 GLUPS for 10–100 dispatches/submit vs
    /// 0.8 GLUPS with a wait per step; anything ≥ ~10 is on the plateau.
    pub submit_chunk: usize,
    submit_chunk_calibrated: bool,
    _l: PhantomData<L>,
}

impl<L: Lattice> WgpuBackend<L> {
    /// Compile the generated WGSL and build the pipeline set.
    pub fn new(ctx: Arc<GpuContext>) -> Self {
        Self::with_config(ctx, KernelCfg::default())
    }

    /// Compile kernels with an explicit storage configuration.
    pub fn with_config(ctx: Arc<GpuContext>, cfg: KernelCfg) -> Self {
        let storage = match cfg.storage {
            GpuStorage::F32 => wgsl::Storage::F32,
            GpuStorage::F16 => wgsl::Storage::F16,
        };
        let source = wgsl::generate_with_storage::<L>(storage);
        let module = ctx
            .device
            .create_shader_module(wgpu::ShaderModuleDescriptor {
                label: Some("lbm-core-gpu"),
                source: wgpu::ShaderSource::Wgsl(source.into()),
            });
        let storage_ro = wgpu::BindingType::Buffer {
            ty: wgpu::BufferBindingType::Storage { read_only: true },
            has_dynamic_offset: false,
            min_binding_size: None,
        };
        let storage_rw = wgpu::BindingType::Buffer {
            ty: wgpu::BufferBindingType::Storage { read_only: false },
            has_dynamic_offset: false,
            min_binding_size: None,
        };
        let uniform = wgpu::BindingType::Buffer {
            ty: wgpu::BufferBindingType::Uniform,
            has_dynamic_offset: false,
            min_binding_size: None,
        };
        let fused_entries = [
            (0, uniform.clone()),
            (1, storage_ro.clone()),
            (2, storage_rw.clone()),
            (3, storage_ro.clone()),
            (4, storage_ro.clone()),
            (5, storage_ro.clone()),
            (6, storage_ro.clone()),
            (7, storage_rw.clone()),
            (8, storage_rw.clone()),
            (9, storage_rw.clone()),
            (10, storage_rw.clone()),
            (11, storage_rw.clone()),
            (14, storage_rw.clone()),
            (15, storage_rw.clone()),
        ]
        .map(|(binding, ty)| wgpu::BindGroupLayoutEntry {
            binding,
            visibility: wgpu::ShaderStages::COMPUTE,
            ty,
            count: None,
        });
        let fused_bgl = ctx
            .device
            .create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
                label: Some("fused-layout"),
                entries: &fused_entries,
            });
        let fused_layout = ctx
            .device
            .create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
                label: Some("fused-pipeline-layout"),
                bind_group_layouts: &[&fused_bgl],
                push_constant_ranges: &[],
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
        let mk_fused = |entry: &str| {
            ctx.device
                .create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
                    label: Some(entry),
                    layout: Some(&fused_layout),
                    module: &module,
                    entry_point: Some(entry),
                    compilation_options: Default::default(),
                    cache: None,
                })
        };
        let pipelines = Arc::new(Pipelines {
            step: mk_fused("step"),
            step_cached: mk_fused("step_cached"),
            step_wale: mk_fused("step_wale"),
            step_cached_wale: mk_fused("step_cached_wale"),
            step_cumulant: mk_fused("step_cumulant"),
            step_cached_cumulant: mk_fused("step_cached_cumulant"),
            step_wale_cumulant: mk_fused("step_wale_cumulant"),
            step_cached_wale_cumulant: mk_fused("step_cached_wale_cumulant"),
            wale_omega: mk("wale_omega"),
            moments: mk("moments"),
            bc: mk("bc"),
            clear_probe: mk("clear_probe"),
        });
        Self {
            ctx,
            cfg,
            pipelines,
            submissions: Cell::new(0),
            cache_next_upload_moments_once: Cell::new(false),
            submit_chunk: 200,
            submit_chunk_calibrated: false,
            _l: PhantomData,
        }
    }

    /// The shared device context.
    pub fn context(&self) -> &Arc<GpuContext> {
        &self.ctx
    }

    /// Kernel/storage configuration.
    pub fn config(&self) -> KernelCfg {
        self.cfg
    }

    /// Number of queue submissions issued by this backend instance.
    pub fn submissions(&self) -> u64 {
        self.submissions.get()
    }

    fn submit(&self, command_buffer: wgpu::CommandBuffer) {
        self.ctx.queue.submit(Some(command_buffer));
        self.submissions.set(self.submissions.get() + 1);
    }

    /// Make the next uploaded state use its host-staged moments for one
    /// collision. This is needed when a force field is edited after the last
    /// moment refresh: CPU collision consumes the already-staged moments once.
    pub fn cache_next_upload_moments_once(&self) {
        self.cache_next_upload_moments_once.set(true);
    }

    /// Cancel a pending one-shot cached-moment upload marker.
    pub fn clear_cached_moment_upload_marker(&self) {
        self.cache_next_upload_moments_once.set(false);
    }

    fn workgroups(&self, fields: &GpuFields) -> (u32, u32, u32) {
        (
            fields.nx.div_ceil(wgsl::WG.0),
            fields.ny.div_ceil(wgsl::WG.1),
            fields.nz,
        )
    }

    fn bc_extent(&self, fields: &GpuFields, face: usize) -> u32 {
        let (t1, t2) = Face::ALL[face].tangents();
        let dims = [fields.nx, fields.ny, fields.nz];
        dims[t1] * dims[t2]
    }

    /// Materialise recorded ops into one encoder/pass and submit — without
    /// waiting (waits only happen in the readback paths).
    pub fn flush(&self, fields: &GpuFields) {
        self.try_flush(fields).expect("GPU submit failed");
    }

    /// Materialise recorded ops and submit them, propagating prior device loss.
    pub fn try_flush(&self, fields: &GpuFields) -> Result<(), GpuError> {
        self.ctx.check_device_lost()?;
        let mut st = fields.state.borrow_mut();
        if st.ops.is_empty() {
            return Ok(());
        }
        let (gx, gy, gz) = self.workgroups(fields);
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
                    Op::Fused {
                        bg,
                        cached,
                        wale,
                        cumulant,
                    } => {
                        if wale {
                            let wb = fields.wale.as_ref().expect("WALE op without WALE buffers");
                            pass.set_pipeline(match (cached, cumulant) {
                                (true, true) => &self.pipelines.step_cached_wale_cumulant,
                                (false, true) => &self.pipelines.step_wale_cumulant,
                                (true, false) => &self.pipelines.step_cached_wale,
                                (false, false) => &self.pipelines.step_wale,
                            });
                            let bind_group = if cached {
                                &wb.fused_cached_bg[bg]
                            } else {
                                &wb.fused_bg[bg]
                            };
                            pass.set_bind_group(0, bind_group, &[]);
                        } else {
                            pass.set_pipeline(match (cached, cumulant) {
                                (true, true) => &self.pipelines.step_cached_cumulant,
                                (false, true) => &self.pipelines.step_cumulant,
                                (true, false) => &self.pipelines.step_cached,
                                (false, false) => &self.pipelines.step,
                            });
                            let bind_group = if cached {
                                &fields.fused_cached_bg[bg]
                            } else {
                                &fields.fused_bg[bg]
                            };
                            pass.set_bind_group(0, bind_group, &[]);
                        }
                        pass.dispatch_workgroups(gx, gy, gz);
                    }
                    Op::WaleOmega { bg } => {
                        let wb = fields.wale.as_ref().expect("WALE op without WALE buffers");
                        pass.set_pipeline(&self.pipelines.wale_omega);
                        pass.set_bind_group(0, &wb.wale_bg[bg], &[]);
                        pass.dispatch_workgroups(gx, gy, gz);
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
                        pass.dispatch_workgroups(gx, gy, gz);
                    }
                }
            }
        }
        self.submit(enc.finish());
        st.ops.clear();
        st.steps_recorded = 0;
        Ok(())
    }

    fn ensure_staging(&self, fields: &GpuFields, bytes: u64) {
        let mut staging = fields.staging.borrow_mut();
        if staging.as_ref().is_some_and(|s| s.size >= bytes) {
            return;
        }
        let size = bytes.max(4);
        let buffer = self.ctx.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("staging"),
            size,
            usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        *staging = Some(StagingBuffer { buffer, size });
    }

    fn map_staging(&self, staging: &wgpu::Buffer, bytes: u64) -> Result<Vec<u8>, GpuError> {
        let slice = staging.slice(..bytes);
        let (tx, rx) = futures_intrusive::channel::shared::oneshot_channel();
        slice.map_async(wgpu::MapMode::Read, move |res| {
            let _ = tx.send(res);
        });
        self.ctx.try_wait_idle()?;
        pollster::block_on(rx.receive())
            .ok_or_else(|| GpuError::Map("map_async callback dropped".to_string()))?
            .map_err(|e| GpuError::Map(e.to_string()))?;
        let out = slice.get_mapped_range().to_vec();
        staging.unmap();
        Ok(out)
    }

    /// Explicit readback of the current deviation populations, compact SoA
    /// layout `f[q*n + y*nx + x]`. Cached until the next recorded step.
    pub fn read_f(&self, fields: &GpuFields) -> Vec<f32> {
        self.try_read_f(fields)
            .expect("GPU population readback failed")
            .as_ref()
            .clone()
    }

    /// Fallible variant of [`Self::read_f`].
    pub fn try_read_f(&self, fields: &GpuFields) -> Result<Arc<Vec<f32>>, GpuError> {
        {
            let st = fields.state.borrow();
            if let Some((gen, data)) = &st.f_cache {
                if *gen == st.generation && st.ops.is_empty() {
                    return Ok(Arc::clone(data));
                }
            }
        }
        self.try_flush(fields)?;
        let bytes = L::Q as u64 * fields.n as u64 * fields.element_bytes;
        self.ensure_staging(fields, bytes);
        let staging_ref = fields.staging.borrow();
        let staging = &staging_ref.as_ref().expect("staging buffer exists").buffer;
        let mut enc = self
            .ctx
            .device
            .create_command_encoder(&wgpu::CommandEncoderDescriptor {
                label: Some("read-f"),
            });
        enc.copy_buffer_to_buffer(&fields.f[fields.cur()], 0, staging, 0, bytes);
        self.submit(enc.finish());
        let raw = self.map_staging(staging, bytes)?;
        let data = Arc::new(decode_storage(self.cfg.storage, &raw));
        let mut st = fields.state.borrow_mut();
        let generation = st.generation;
        st.f_cache = Some((generation, Arc::clone(&data)));
        Ok(data)
    }

    /// Explicit readback of the momentum-exchange probe force accumulated
    /// during the most recent executed step (V1 `probed_force` semantics).
    pub fn read_probed_force(&self, fields: &GpuFields) -> [f32; 3] {
        self.try_read_probed_force(fields)
            .expect("GPU probe-force readback failed")
    }

    /// Fallible variant of [`Self::read_probed_force`].
    pub fn try_read_probed_force(&self, fields: &GpuFields) -> Result<[f32; 3], GpuError> {
        self.try_flush(fields)?;
        self.ensure_staging(fields, 12);
        let staging_ref = fields.staging.borrow();
        let staging = &staging_ref.as_ref().expect("staging buffer exists").buffer;
        let mut enc = self
            .ctx
            .device
            .create_command_encoder(&wgpu::CommandEncoderDescriptor {
                label: Some("read-probe"),
            });
        enc.copy_buffer_to_buffer(&fields.probe_acc, 0, staging, 0, 12);
        self.submit(enc.finish());
        let raw = self.map_staging(staging, 12)?;
        let v: &[f32] = bytemuck::cast_slice(&raw);
        Ok([v[0], v[1], v[2]])
    }

    /// Current submit chunk after any runtime calibration.
    pub fn submit_chunk(&self) -> usize {
        self.submit_chunk
    }

    /// Force a submit chunk and mark it calibrated by the caller.
    pub fn set_submit_chunk(&mut self, chunk: usize) {
        self.submit_chunk = chunk.clamp(1, 200);
        self.submit_chunk_calibrated = true;
    }

    /// Whether the first measured submit has already calibrated chunk size.
    pub fn submit_chunk_calibrated(&self) -> bool {
        self.submit_chunk_calibrated
    }

    /// Update chunk size from one measured submit.
    pub fn calibrate_submit_chunk(&mut self, measured_steps: usize, elapsed: Duration) {
        self.submit_chunk = calibrate_submit_chunk(measured_steps, elapsed, self.submit_chunk);
        self.submit_chunk_calibrated = true;
    }

    /// Write the global Params uniform once; assert the step parameters do
    /// not change between steps of a run (V1 semantics: relaxation rates,
    /// force and face BCs are fixed per solver lifetime).
    fn ensure_params(&self, sub: &Subdomain, fields: &GpuFields, p: &StepParams<f32>) {
        // Relaxation constants exactly as KParams::new builds them: f64
        // parameters converted to f32 once.
        let use_cached_moments = fields.state.borrow().use_cached_moments_once;
        let collision_code = match p.collision {
            CollisionKind::Bgk => 0u32,
            CollisionKind::Trt { .. } => 1u32,
            CollisionKind::Cumulant { .. } => 2u32,
        };
        let gravity = p.gravity.unwrap_or([0.0; 3]);
        let words: [u32; 20] = [
            fields.nx,
            fields.ny,
            fields.nz,
            0,
            (p.omega_p as f32).to_bits(),
            (p.omega_m as f32).to_bits(),
            ((1.0 - p.omega_p / 2.0) as f32).to_bits(),
            ((1.0 - p.omega_m / 2.0) as f32).to_bits(),
            p.force[0].to_bits(),
            p.force[1].to_bits(),
            p.force[2].to_bits(),
            gravity[0].to_bits(),
            gravity[1].to_bits(),
            gravity[2].to_bits(),
            {
                let halo = sub.halo_flags();
                let mut flags = 0u32;
                for (i, &h) in halo.iter().take(2 * L::D).enumerate() {
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
                if use_cached_moments {
                    flags |= wgsl::FLAG_CACHED_MOMENTS;
                }
                if fields.wale.is_some() {
                    flags |= wgsl::FLAG_WALE;
                }
                if p.gravity.is_some() {
                    flags |= wgsl::FLAG_GRAVITY;
                }
                for face in Face::ALL {
                    if face.axis() >= L::D {
                        continue;
                    }
                    if p.faces[face.index()].is_open() {
                        flags |= wgsl::FLAG_OPEN_FACE[face.index()];
                    }
                }
                flags
            },
            (p.collision.omega_shear((1.0 / p.omega_p - 0.5) / 3.0) as f32).to_bits(),
            collision_code,
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

    fn bc_words(&self, sub: &Subdomain, fields: &GpuFields, p: &StepParams<f32>) -> [[u32; 64]; 6] {
        let dims = [fields.nx, fields.ny, fields.nz];
        let strides = [1u32, fields.nx, fields.nx * fields.ny];
        let mut out = [[0u32; 64]; 6];
        for face in Face::ALL {
            let fi = face.index();
            if face.axis() >= L::D || !sub.touches_global_face(face) {
                continue;
            }
            let bc = &p.faces[fi];
            if !bc.is_open() {
                continue;
            }
            let a = face.axis();
            let (t1, t2) = face.tangents();
            let base = if face.is_neg() {
                0
            } else {
                (dims[a] - 1) * strides[a]
            };
            let stride1 = strides[t1];
            let stride2 = strides[t2];
            let extent1 = dims[t1];
            let extent = dims[t1] * dims[t2];
            let joff = if face.is_neg() {
                strides[a] as i32
            } else {
                -(strides[a] as i32)
            };
            let n_in = face.n_in();
            let unit = |axis: usize, sign: i8| -> [i8; 3] {
                let mut v = [0i8; 3];
                v[axis] = sign;
                v
            };
            let add = |p: [i8; 3], q: [i8; 3]| [p[0] + q[0], p[1] + q[1], p[2] + q[2]];
            let q_n = L::dir_index(n_in);
            let q_p1 = L::dir_index(add(n_in, unit(t1, 1)));
            let q_m1 = L::dir_index(add(n_in, unit(t1, -1)));
            let q_t1 = L::dir_index(unit(t1, 1));
            let q_mt1 = L::dir_index(unit(t1, -1));
            let q_p2 = if L::D == 3 {
                L::dir_index(add(n_in, unit(t2, 1)))
            } else {
                q_p1
            };
            let q_m2 = if L::D == 3 {
                L::dir_index(add(n_in, unit(t2, -1)))
            } else {
                q_m1
            };
            let q_t2 = if L::D == 3 {
                L::dir_index(unit(t2, 1))
            } else {
                q_t1
            };
            let q_mt2 = if L::D == 3 {
                L::dir_index(unit(t2, -1))
            } else {
                q_mt1
            };
            let q_pp = if L::D == 3 {
                L::dir_index(add(unit(t1, 1), unit(t2, 1)))
            } else {
                q_t1
            };
            let q_pm = if L::D == 3 {
                L::dir_index(add(unit(t1, 1), unit(t2, -1)))
            } else {
                q_t1
            };
            let q_mp = if L::D == 3 {
                L::dir_index(add(unit(t1, -1), unit(t2, 1)))
            } else {
                q_mt1
            };
            let q_mm = if L::D == 3 {
                L::dir_index(add(unit(t1, -1), unit(t2, -1)))
            } else {
                q_mt1
            };
            let unk = L::unknowns(face);
            assert!(
                unk.len() == 3 || unk.len() == 5,
                "GPU open face unknown count must be 3 or 5"
            );
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
            w[2] = stride1;
            w[3] = extent;
            w[4] = joff as u32;
            w[5] = u32::from(fields.profile_set[fi]);
            w[6] = stride2;
            w[7] = extent1;
            w[8] = q_n as u32;
            w[9] = L::OPP[q_n] as u32;
            w[10] = q_p1 as u32;
            w[11] = L::OPP[q_p1] as u32;
            w[12] = q_m1 as u32;
            w[13] = L::OPP[q_m1] as u32;
            w[14] = q_p2 as u32;
            w[15] = L::OPP[q_p2] as u32;
            w[16] = q_m2 as u32;
            w[17] = L::OPP[q_m2] as u32;
            w[18] = q_t1 as u32;
            w[19] = q_mt1 as u32;
            w[20] = q_t2 as u32;
            w[21] = q_mt2 as u32;
            w[22] = q_pp as u32;
            w[23] = q_pm as u32;
            w[24] = q_mp as u32;
            w[25] = q_mm as u32;
            for (k, &q) in unk.iter().enumerate() {
                w[26 + k] = q as u32;
            }
            w[31] = unk.len() as u32;
            let p2 = match *bc {
                FaceBC::Velocity { u } => u[2],
                _ => 0.0,
            };
            w[32] = p0.to_bits();
            w[33] = p1.to_bits();
            w[34] = p2.to_bits();
            w[35] = (n_in[0] as f32).to_bits();
            w[36] = (n_in[1] as f32).to_bits();
            w[37] = (n_in[2] as f32).to_bits();
            w[38] = (unit(t1, 1)[0] as f32).to_bits();
            w[39] = (unit(t1, 1)[1] as f32).to_bits();
            w[40] = (unit(t1, 1)[2] as f32).to_bits();
            w[41] = (unit(t2, 1)[0] as f32).to_bits();
            w[42] = (unit(t2, 1)[1] as f32).to_bits();
            w[43] = (unit(t2, 1)[2] as f32).to_bits();
            for (k, &q) in unk.iter().enumerate() {
                w[44 + k] = (L::W[q] as f32).to_bits();
            }
            w[49] = (ws as f32).to_bits();
            w[50] = cinv.to_bits();
        }
        out
    }

    fn storage_buffer(&self, label: &str, size: u64, copy_dst: bool) -> wgpu::Buffer {
        let mut usage = wgpu::BufferUsages::STORAGE;
        if copy_dst {
            usage |= wgpu::BufferUsages::COPY_DST;
        }
        self.ctx.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some(label),
            size: size.max(8),
            usage,
            mapped_at_creation: false,
        })
    }

    fn rebuild_field_bind_groups(&self, fields: &mut GpuFields) {
        fn e(binding: u32, b: &wgpu::Buffer) -> wgpu::BindGroupEntry<'_> {
            wgpu::BindGroupEntry {
                binding,
                resource: b.as_entire_binding(),
            }
        }
        let device = &self.ctx.device;
        let step_layout = self.pipelines.step.get_bind_group_layout(0);
        let step_cached_layout = self.pipelines.step_cached.get_bind_group_layout(0);
        fields.fused_bg = [0usize, 1].map(|p| {
            device.create_bind_group(&wgpu::BindGroupDescriptor {
                label: Some("fused"),
                layout: &step_layout,
                entries: &[
                    e(0, &fields.params_ub),
                    e(1, &fields.f[p]),
                    e(2, &fields.f[1 - p]),
                    e(3, &fields.mask),
                    e(4, &fields.wall_u),
                    e(5, &fields.force_field),
                    e(6, &fields.stash[p]),
                    e(7, &fields.stash[1 - p]),
                    e(8, &fields.probe_acc),
                    e(9, &fields.rho),
                    e(10, &fields.ux),
                    e(11, &fields.uy),
                    e(14, &fields.uz),
                    e(15, fields.wale.as_ref().map_or(&fields.rho, |w| &w.omega)),
                ],
            })
        });
        fields.fused_cached_bg = [0usize, 1].map(|p| {
            device.create_bind_group(&wgpu::BindGroupDescriptor {
                label: Some("fused_cached"),
                layout: &step_cached_layout,
                entries: &[
                    e(0, &fields.params_ub),
                    e(1, &fields.f[p]),
                    e(2, &fields.f[1 - p]),
                    e(3, &fields.mask),
                    e(4, &fields.wall_u),
                    e(5, &fields.force_field),
                    e(6, &fields.stash[p]),
                    e(7, &fields.stash[1 - p]),
                    e(8, &fields.probe_acc),
                    e(9, &fields.rho),
                    e(10, &fields.ux),
                    e(11, &fields.uy),
                    e(14, &fields.uz),
                    e(15, fields.wale.as_ref().map_or(&fields.rho, |w| &w.omega)),
                ],
            })
        });
        let moments_layout = self.pipelines.moments.get_bind_group_layout(0);
        fields.moments_bg = [0usize, 1].map(|p| {
            device.create_bind_group(&wgpu::BindGroupDescriptor {
                label: Some("moments"),
                layout: &moments_layout,
                entries: &[
                    e(0, &fields.params_ub),
                    e(1, &fields.f[p]),
                    e(3, &fields.mask),
                    e(5, &fields.force_field),
                    e(9, &fields.rho),
                    e(10, &fields.ux),
                    e(11, &fields.uy),
                    e(14, &fields.uz),
                ],
            })
        });
        let bc_layout = self.pipelines.bc.get_bind_group_layout(0);
        fields.bc_bg = std::array::from_fn(|face| {
            [0usize, 1].map(|p| {
                device.create_bind_group(&wgpu::BindGroupDescriptor {
                    label: Some("bc"),
                    layout: &bc_layout,
                    entries: &[
                        e(0, &fields.params_ub),
                        e(2, &fields.f[p]),
                        e(3, &fields.mask),
                        e(5, &fields.force_field),
                        e(9, &fields.rho),
                        e(10, &fields.ux),
                        e(11, &fields.uy),
                        e(14, &fields.uz),
                        e(12, &fields.bc_ub[face]),
                        e(13, &fields.profiles[face]),
                    ],
                })
            })
        });
        if fields.wale.is_some() {
            let omega = fields
                .wale
                .as_ref()
                .expect("WALE buffer exists")
                .omega
                .clone();
            fields.wale = Some(self.create_wale_bind_groups(fields, &omega));
        }
    }

    fn create_wale_bind_groups(&self, fields: &GpuFields, omega: &wgpu::Buffer) -> WaleBuffers {
        fn e(binding: u32, b: &wgpu::Buffer) -> wgpu::BindGroupEntry<'_> {
            wgpu::BindGroupEntry {
                binding,
                resource: b.as_entire_binding(),
            }
        }
        let device = &self.ctx.device;
        let step_layout = self.pipelines.step_wale.get_bind_group_layout(0);
        let step_cached_layout = self.pipelines.step_cached_wale.get_bind_group_layout(0);
        let fused_bg = [0usize, 1].map(|p| {
            device.create_bind_group(&wgpu::BindGroupDescriptor {
                label: Some("fused-wale"),
                layout: &step_layout,
                entries: &[
                    e(0, &fields.params_ub),
                    e(1, &fields.f[p]),
                    e(2, &fields.f[1 - p]),
                    e(3, &fields.mask),
                    e(4, &fields.wall_u),
                    e(5, &fields.force_field),
                    e(6, &fields.stash[p]),
                    e(7, &fields.stash[1 - p]),
                    e(8, &fields.probe_acc),
                    e(9, &fields.rho),
                    e(10, &fields.ux),
                    e(11, &fields.uy),
                    e(14, &fields.uz),
                    e(15, omega),
                ],
            })
        });
        let fused_cached_bg = [0usize, 1].map(|p| {
            device.create_bind_group(&wgpu::BindGroupDescriptor {
                label: Some("fused-cached-wale"),
                layout: &step_cached_layout,
                entries: &[
                    e(0, &fields.params_ub),
                    e(1, &fields.f[p]),
                    e(2, &fields.f[1 - p]),
                    e(3, &fields.mask),
                    e(4, &fields.wall_u),
                    e(5, &fields.force_field),
                    e(6, &fields.stash[p]),
                    e(7, &fields.stash[1 - p]),
                    e(8, &fields.probe_acc),
                    e(9, &fields.rho),
                    e(10, &fields.ux),
                    e(11, &fields.uy),
                    e(14, &fields.uz),
                    e(15, omega),
                ],
            })
        });
        let wale_layout = self.pipelines.wale_omega.get_bind_group_layout(0);
        let wale_bg = [0usize, 1].map(|p| {
            device.create_bind_group(&wgpu::BindGroupDescriptor {
                label: Some("wale-omega"),
                layout: &wale_layout,
                entries: &[
                    e(0, &fields.params_ub),
                    e(1, &fields.f[p]),
                    e(3, &fields.mask),
                    e(4, &fields.wall_u),
                    e(5, &fields.force_field),
                    e(9, &fields.rho),
                    e(10, &fields.ux),
                    e(11, &fields.uy),
                    e(14, &fields.uz),
                    e(15, omega),
                ],
            })
        });
        WaleBuffers {
            omega: omega.clone(),
            fused_bg,
            fused_cached_bg,
            wale_bg,
        }
    }

    /// Enable or disable on-device WALE omega generation for these fields.
    pub fn set_wale_enabled(&self, fields: &mut GpuFields, enabled: bool, base_omega: f32) {
        self.flush(fields);
        if enabled {
            if fields.wale.is_none() {
                let omega = self.ctx.device.create_buffer(&wgpu::BufferDescriptor {
                    label: Some("wale_omega"),
                    size: ((fields.n * 4) as u64).max(8),
                    usage: wgpu::BufferUsages::STORAGE
                        | wgpu::BufferUsages::COPY_DST
                        | wgpu::BufferUsages::COPY_SRC,
                    mapped_at_creation: false,
                });
                self.ctx.queue.write_buffer(
                    &omega,
                    0,
                    bytemuck::cast_slice(&vec![base_omega; fields.n]),
                );
                fields.wale = Some(self.create_wale_bind_groups(fields, &omega));
            }
        } else {
            fields.wale = None;
        }
        let mut st = fields.state.borrow_mut();
        st.params_words = None;
        st.f_cache = None;
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
        assert_eq!(g.core[2] as u32, fields.nz);
        let (nx, ny, nz, n) = (g.core[0], g.core[1], g.core[2], fields.n);
        let xy = nx * ny;
        let np = g.n_padded();
        let q = self.ctx.queue.clone();
        {
            let st = fields.state.borrow();
            assert!(
                st.ops.is_empty() && !st.pending_collide,
                "upload with recorded but unsubmitted steps"
            );
        }
        let use_cached_moments_once = self.cache_next_upload_moments_once.replace(false);
        // Populations: current -> f[cur], ping-pong partner -> f[1-cur].
        let cur = fields.cur();
        let mut buf = vec![0f32; L::Q * n];
        for (src, dst) in [(&host.f, cur), (&host.ftmp, 1 - cur)] {
            for qi in 0..L::Q {
                for z in 0..nz {
                    for y in 0..ny {
                        for x in 0..nx {
                            let c = z * xy + y * nx + x;
                            buf[qi * n + c] = src[qi * np + g.pidx(x, y, z)];
                        }
                    }
                }
            }
            q.write_buffer(&fields.f[dst], 0, &encode_storage(self.cfg.storage, &buf));
        }
        // Edge stash (stash_in of the next fused dispatch = stash[cur]):
        // the ftmp values V1 would leave in the skipped slots.
        let slen = wgsl::stash_len::<L>(nx, ny, nz);
        let mut stash = vec![0f32; slen];
        let mut off = 0usize;
        for face in Face::ALL {
            if face.axis() >= L::D {
                continue;
            }
            let unk = L::unknowns(face);
            let (t1, t2) = face.tangents();
            let dims = [nx, ny, nz];
            let ext = dims[t1] * dims[t2];
            if !sub.has_halo(face) {
                for (k, &u) in unk.iter().enumerate() {
                    for c2 in 0..dims[t2] {
                        for c1 in 0..dims[t1] {
                            let t = c2 * dims[t1] + c1;
                            let mut pos = [0usize; 3];
                            pos[face.axis()] = if face.is_neg() {
                                0
                            } else {
                                dims[face.axis()] - 1
                            };
                            pos[t1] = c1;
                            pos[t2] = c2;
                            stash[off + k * ext + t] =
                                host.ftmp[u * np + g.pidx(pos[0], pos[1], pos[2])];
                        }
                    }
                }
            }
            off += unk.len() * ext;
        }
        q.write_buffer(
            &fields.stash[cur],
            0,
            &encode_storage(self.cfg.storage, &stash),
        );
        q.write_buffer(
            &fields.stash[1 - cur],
            0,
            &encode_storage(self.cfg.storage, &vec![0f32; slen]),
        );
        // Mask (bit0 solid, bit1 probe) + host copies for reduce().
        let mut mask = vec![0u32; n];
        let mut host_solid = vec![false; n];
        for z in 0..nz {
            for y in 0..ny {
                for x in 0..nx {
                    let pi = g.pidx(x, y, z);
                    let c = z * xy + y * nx + x;
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
        }
        q.write_buffer(&fields.mask, 0, bytemuck::cast_slice(&mask));
        fields.host_solid = host_solid;
        fields.has_probe = host.probe.is_some();
        let needs_wall_u = fields.host_solid.iter().any(|&solid| solid);
        let needs_force_field = host.force_field.is_some();
        let mut rebuild_bgs = false;
        if needs_wall_u && !fields.wall_u_full {
            fields.wall_u = self.storage_buffer("wall_u", (n * 16) as u64, true);
            fields.wall_u_full = true;
            rebuild_bgs = true;
        }
        if needs_force_field && !fields.force_field_full {
            fields.force_field = self.storage_buffer("force_field", (n * 16) as u64, true);
            fields.force_field_full = true;
            rebuild_bgs = true;
        }
        if rebuild_bgs {
            self.rebuild_field_bind_groups(fields);
        }
        // Wall velocities (read only at solid neighbours).
        if fields.wall_u_full {
            let mut wu = vec![0f32; 4 * n];
            for z in 0..nz {
                for y in 0..ny {
                    for x in 0..nx {
                        let c = z * xy + y * nx + x;
                        let v = host.wall_u[g.pidx(x, y, z)];
                        wu[4 * c] = v[0];
                        wu[4 * c + 1] = v[1];
                        wu[4 * c + 2] = v[2];
                    }
                }
            }
            q.write_buffer(&fields.wall_u, 0, bytemuck::cast_slice(&wu));
        }
        // Per-cell force field (compact already).
        fields.host_ff = host.force_field.clone();
        if let Some(ff) = &host.force_field {
            let mut fv = vec![0f32; 4 * n];
            for (c, v) in ff.iter().enumerate() {
                fv[4 * c] = v[0];
                fv[4 * c + 1] = v[1];
                fv[4 * c + 2] = v[2];
            }
            q.write_buffer(&fields.force_field, 0, bytemuck::cast_slice(&fv));
        }
        // Moments (compact already; carries V1's values at solid cells,
        // which the moments kernel never rewrites).
        q.write_buffer(&fields.rho, 0, bytemuck::cast_slice(&host.rho));
        q.write_buffer(&fields.ux, 0, bytemuck::cast_slice(&host.ux));
        q.write_buffer(&fields.uy, 0, bytemuck::cast_slice(&host.uy));
        q.write_buffer(&fields.uz, 0, bytemuck::cast_slice(&host.uz));
        // Inlet profiles.
        for face in Face::ALL {
            let fi = face.index();
            fields.profile_set[fi] = false;
            if let Some(prof) = &host.inlet_profiles[fi] {
                let mut pv = vec![0f32; 4 * prof.len()];
                for (t, u) in prof.iter().enumerate() {
                    pv[4 * t] = u[0];
                    pv[4 * t + 1] = u[1];
                    pv[4 * t + 2] = u[2];
                }
                q.write_buffer(&fields.profiles[fi], 0, bytemuck::cast_slice(&pv));
                fields.profile_set[fi] = true;
            }
        }
        // Probe accumulator and cached uniforms reset (masks/probe/profile
        // presence may have changed the flags).
        q.write_buffer(&fields.probe_acc, 0, &[0u8; 12]);
        let mut st = fields.state.borrow_mut();
        st.use_cached_moments_once = use_cached_moments_once;
        st.params_words = None;
        st.bc_words = None;
        st.f_cache = None;
        st.generation += 1;
    }
}

impl<L: Lattice> WgpuBackend<L> {
    /// Allocate device fields after validating adapter/device resource limits.
    pub fn try_alloc(&self, sub: &Subdomain) -> Result<GpuFields, GpuError> {
        self.try_alloc_with_options(sub, true, true)
    }

    /// Allocate fields with optional dummy wall/force buffers for known-absent features.
    pub fn try_alloc_with_options(
        &self,
        sub: &Subdomain,
        full_wall_u: bool,
        full_force_field: bool,
    ) -> Result<GpuFields, GpuError> {
        let g = sub.geom;
        assert_eq!(g.d, L::D, "WgpuBackend lattice/domain dimension mismatch");
        let (nx, ny, nz) = (g.core[0] as u32, g.core[1] as u32, g.core[2] as u32);
        let n = (nx as usize) * (ny as usize) * (nz as usize);
        let element_bytes = self.cfg.storage.element_bytes();
        let plan =
            GpuResourcePlan::for_grid::<L>(nx as usize, ny as usize, nz as usize, element_bytes)?;
        plan.validate(&self.ctx.device.limits())?;
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
        let fbytes = plan.fbytes;
        let f = [
            buf("f0", fbytes, U::STORAGE | U::COPY_DST | U::COPY_SRC),
            buf("f1", fbytes, U::STORAGE | U::COPY_DST | U::COPY_SRC),
        ];
        let slen =
            wgsl::stash_len::<L>(nx as usize, ny as usize, nz as usize) as u64 * element_bytes;
        let stash = [
            buf("stash0", slen, U::STORAGE | U::COPY_DST),
            buf("stash1", slen, U::STORAGE | U::COPY_DST),
        ];
        let mask = buf("mask", (n * 4) as u64, U::STORAGE | U::COPY_DST);
        let wall_u_size = if full_wall_u { (n * 16) as u64 } else { 16 };
        let force_field_size = if full_force_field {
            (n * 16) as u64
        } else {
            16
        };
        let wall_u = buf("wall_u", wall_u_size, U::STORAGE | U::COPY_DST);
        let force_field = buf("force_field", force_field_size, U::STORAGE | U::COPY_DST);
        let rho = buf(
            "rho",
            (n * 4) as u64,
            U::STORAGE | U::COPY_DST | U::COPY_SRC,
        );
        let ux = buf("ux", (n * 4) as u64, U::STORAGE | U::COPY_DST | U::COPY_SRC);
        let uy = buf("uy", (n * 4) as u64, U::STORAGE | U::COPY_DST | U::COPY_SRC);
        let uz = buf("uz", (n * 4) as u64, U::STORAGE | U::COPY_DST | U::COPY_SRC);
        let probe_acc = buf("probe_acc", 12, U::STORAGE | U::COPY_DST | U::COPY_SRC);
        let params_ub = buf("params", 80, U::UNIFORM | U::COPY_DST);
        let bc_ub = std::array::from_fn(|i| buf(&format!("bc{i}"), 256, U::UNIFORM | U::COPY_DST));
        let profiles = std::array::from_fn(|i| {
            let (t1, t2) = Face::ALL[i].tangents();
            let dims = [nx, ny, nz];
            let ext = (dims[t1] * dims[t2]) as u64;
            buf(
                &format!("profile{i}"),
                (ext * 16).max(16),
                U::STORAGE | U::COPY_DST,
            )
        });
        fn e(binding: u32, b: &wgpu::Buffer) -> wgpu::BindGroupEntry<'_> {
            wgpu::BindGroupEntry {
                binding,
                resource: b.as_entire_binding(),
            }
        }
        let step_layout = self.pipelines.step.get_bind_group_layout(0);
        let step_cached_layout = self.pipelines.step_cached.get_bind_group_layout(0);
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
                    e(9, &rho),
                    e(10, &ux),
                    e(11, &uy),
                    e(14, &uz),
                    e(15, &rho),
                ],
            })
        });
        let fused_cached_bg = [0usize, 1].map(|p| {
            device.create_bind_group(&wgpu::BindGroupDescriptor {
                label: Some("fused_cached"),
                layout: &step_cached_layout,
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
                    e(9, &rho),
                    e(10, &ux),
                    e(11, &uy),
                    e(14, &uz),
                    e(15, &rho),
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
                    e(14, &uz),
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
                        e(5, &force_field),
                        e(9, &rho),
                        e(10, &ux),
                        e(11, &uy),
                        e(14, &uz),
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

        Ok(GpuFields {
            nx,
            ny,
            nz,
            n,
            element_bytes,
            f,
            stash,
            mask,
            wall_u,
            force_field,
            wall_u_full: full_wall_u,
            force_field_full: full_force_field,
            rho,
            ux,
            uy,
            uz,
            probe_acc,
            params_ub,
            bc_ub,
            profiles,
            staging: RefCell::new(None),
            fused_bg,
            fused_cached_bg,
            moments_bg,
            bc_bg,
            clear_bg,
            wale: None,
            state: RefCell::new(RecState {
                cur: 0,
                ops: Vec::new(),
                steps_recorded: 0,
                pending_collide: false,
                use_cached_moments_once: false,
                params_words: None,
                bc_words: None,
                generation: 0,
                f_cache: None,
            }),
            host_solid: vec![false; n],
            host_ff: None,
            has_probe: false,
            profile_set: [false; 6],
        })
    }
}

impl<L: Lattice> Backend<L, f32> for WgpuBackend<L> {
    type Fields = GpuFields;

    fn supports_localized_features(&self) -> bool {
        false
    }

    fn supports_gravity_body_force(&self) -> bool {
        true
    }

    fn alloc(&self, sub: &Subdomain) -> GpuFields {
        self.try_alloc(sub).expect("GPU field allocation failed")
    }

    fn stage_in(&self, sub: &Subdomain, fields: &mut GpuFields, host: &SoaFields<f32>) {
        self.upload(sub, fields, host);
    }

    fn stage_out(&self, sub: &Subdomain, fields: &GpuFields, host: &mut SoaFields<f32>) {
        let mut hm = HostMoments::default();
        let (f, _) = self
            .try_read_sync(fields, &mut hm)
            .expect("GPU stage-out readback failed");
        let g = host.geom;
        debug_assert_eq!(sub.geom, g);
        let (nx, ny, nz, np) = (g.core[0], g.core[1], g.core[2], g.n_padded());
        let xy = nx * ny;
        let n = xy * nz;
        for q in 0..L::Q {
            for z in 0..nz {
                for y in 0..ny {
                    for x in 0..nx {
                        let c = z * xy + y * nx + x;
                        host.f[q * np + g.pidx(x, y, z)] = f[q * n + c];
                    }
                }
            }
        }
        write_host_moments(g, &hm, host);
    }

    fn handles_single_part_periodic_halo(&self) -> bool {
        true
    }

    fn supports_two_pass(&self) -> bool {
        false
    }

    fn exchange_f<H: HaloExchange<f32>>(
        &mut self,
        _exchange: &H,
        subs: &[Subdomain],
        fields: &mut [GpuFields],
    ) {
        assert_eq!(
            fields.len(),
            1,
            "WgpuBackend B-1 path supports one monolithic part"
        );
        assert_eq!(
            subs.len(),
            1,
            "WgpuBackend B-1 path supports one monolithic subdomain"
        );
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
    ) {
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
        let cached = st.use_cached_moments_once || has_open_faces::<L>(sub, _p);
        let wale = fields.wale.is_some();
        let cumulant = matches!(_p.collision, CollisionKind::Cumulant { .. });
        if wale {
            st.ops.push(Op::Moments { bg: cur });
            st.ops.push(Op::WaleOmega { bg: cur });
        }
        st.ops.push(Op::Fused {
            bg: cur,
            cached,
            wale,
            cumulant,
        });
        st.generation += 1;
        st.f_cache = None;
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

    fn update_moments(&mut self, sub: &Subdomain, fields: &mut GpuFields, p: &StepParams<f32>) {
        let has_uniform_force = p.force.iter().any(|&f| f != 0.0) || p.gravity.is_some();
        if has_uniform_force && fields.host_ff.is_none() {
            self.ensure_params(sub, fields, p);
            let mut st = fields.state.borrow_mut();
            let cur = st.cur;
            st.ops.push(Op::Moments { bg: cur });
        } else {
            // Lazy: the fused kernel re-derives (rho, u) from the identical
            // pre-collide state, so no device moment refresh is needed.
        }
    }

    fn end_step(&mut self, fields: &GpuFields) {
        fields.state.borrow_mut().steps_recorded += 1;
    }

    fn read_probed_force(&self, fields: &GpuFields) -> [f32; 3] {
        WgpuBackend::read_probed_force(self, fields)
    }

    fn run_span<H: HaloExchange<f32>>(
        &mut self,
        _exchange: &H,
        subs: &[Subdomain],
        fields: &mut [GpuFields],
        p: &StepParams<f32>,
        two_pass: bool,
        probed_force: &mut [f32; 3],
        steps: usize,
    ) {
        if steps == 0 {
            return;
        }
        assert!(
            !two_pass,
            "WgpuBackend streams the full grid in one fused dispatch (no two-pass split)"
        );
        assert_eq!(
            fields.len(),
            1,
            "WgpuBackend B-1 path supports one monolithic part"
        );
        assert_eq!(
            subs.len(),
            1,
            "WgpuBackend B-1 path supports one monolithic subdomain"
        );
        let sub = &subs[0];
        let field = &mut fields[0];
        let open_faces = std::array::from_fn::<_, 6, _>(|fi| {
            let face = Face::ALL[fi];
            face.axis() < L::D && sub.touches_global_face(face) && p.faces[fi].is_open()
        });
        self.ensure_params(sub, field, p);
        self.ensure_bc(sub, field, p);
        let mut st = field.state.borrow_mut();
        assert!(
            !st.pending_collide,
            "run_span called with an unfinished collide/stream pair"
        );
        for _ in 0..steps {
            if field.has_probe {
                st.ops.push(Op::ClearProbe);
            }
            let cur = st.cur;
            let cached = st.use_cached_moments_once || open_faces.iter().any(|&is_open| is_open);
            let wale = field.wale.is_some();
            let cumulant = matches!(p.collision, CollisionKind::Cumulant { .. });
            if wale {
                st.ops.push(Op::Moments { bg: cur });
                st.ops.push(Op::WaleOmega { bg: cur });
            }
            st.ops.push(Op::Fused {
                bg: cur,
                cached,
                wale,
                cumulant,
            });
            st.generation += 1;
            st.f_cache = None;
            st.cur ^= 1;
            let cur = st.cur;
            for (face, &is_open) in open_faces.iter().enumerate() {
                if is_open {
                    st.ops.push(Op::Bc { face, bg: cur });
                }
            }
            st.steps_recorded += 1;
        }
        let _ = probed_force;
    }

    fn run_chunk_size(&self, fields: &[GpuFields]) -> usize {
        if fields
            .iter()
            .any(|field| field.state.borrow().use_cached_moments_once)
        {
            return 1;
        }
        self.submit_chunk.max(1)
    }

    fn finish_run_chunk(&mut self, fields: &[GpuFields], steps: usize) {
        let start = std::time::Instant::now();
        let calibrating = !self.submit_chunk_calibrated;
        for field in fields {
            if field.state.borrow().ops.is_empty() {
                continue;
            }
            self.flush(field);
        }
        if calibrating {
            self.ctx.wait_idle();
        }
        for field in fields {
            let mut st = field.state.borrow_mut();
            if st.use_cached_moments_once {
                st.use_cached_moments_once = false;
                st.params_words = None;
            }
        }
        if calibrating {
            self.calibrate_submit_chunk(steps, start.elapsed());
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
        let (nx, ny, nz, n) = (
            fields.nx as usize,
            fields.ny as usize,
            fields.nz as usize,
            fields.n,
        );
        let xy = nx * ny;
        if kind == Reduction::FluidCells {
            return fields.host_solid.iter().filter(|&&solid| !solid).count() as f64;
        }
        let f = self
            .try_read_f(fields)
            .expect("GPU population readback failed");
        let mut acc = 0.0f64;
        for z in 0..nz {
            for y in 0..ny {
                for x in 0..nx {
                    let c = z * xy + y * nx + x;
                    if fields.host_solid[c] {
                        continue;
                    }
                    match kind {
                        Reduction::FluidCells => unreachable!("handled before population readback"),
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
                            let rho = 1.0 + (0..L::Q).map(|q| f[q * n + c] as f64).sum::<f64>();
                            let gravity_force = p.gravity.map_or(0.0, |g| rho * g[a] as f64);
                            let fa = match &fields.host_ff {
                                Some(field) => {
                                    p.force[a] as f64 + (field[c][a] as f64 + gravity_force)
                                }
                                None => p.force[a] as f64 + gravity_force,
                            };
                            acc += m + 0.5 * fa;
                        }
                    }
                }
            }
        }
        acc
    }

    fn read_moments(&self, fields: &GpuFields, out: &mut HostMoments<f32>) {
        self.try_read_moments(fields, out)
            .expect("GPU moments readback failed");
    }
}

impl<L: Lattice> WgpuBackend<L> {
    /// Fallible moment readback used by [`GpuSolver`](super::solver::GpuSolver).
    pub fn try_read_moments(
        &self,
        fields: &GpuFields,
        out: &mut HostMoments<f32>,
    ) -> Result<(), GpuError> {
        {
            let mut st = fields.state.borrow_mut();
            let cur = st.cur;
            st.ops.push(Op::Moments { bg: cur });
        }
        self.try_flush(fields)?;
        let n = fields.n;
        let plane = (n * 4) as u64;
        let bytes = 4 * plane;
        self.ensure_staging(fields, bytes);
        let staging_ref = fields.staging.borrow();
        let staging = &staging_ref.as_ref().expect("staging buffer exists").buffer;
        let mut enc = self
            .ctx
            .device
            .create_command_encoder(&wgpu::CommandEncoderDescriptor {
                label: Some("read-moments"),
            });
        enc.copy_buffer_to_buffer(&fields.rho, 0, staging, 0, plane);
        enc.copy_buffer_to_buffer(&fields.ux, 0, staging, plane, plane);
        enc.copy_buffer_to_buffer(&fields.uy, 0, staging, 2 * plane, plane);
        enc.copy_buffer_to_buffer(&fields.uz, 0, staging, 3 * plane, plane);
        self.submit(enc.finish());
        let raw = self.map_staging(staging, bytes)?;
        let v: &[f32] = bytemuck::cast_slice(&raw);
        out.rho.clear();
        out.rho.extend_from_slice(&v[..n]);
        out.ux.clear();
        out.ux.extend_from_slice(&v[n..2 * n]);
        out.uy.clear();
        out.uy.extend_from_slice(&v[2 * n..3 * n]);
        out.uz.clear();
        out.uz.extend_from_slice(&v[3 * n..4 * n]);
        Ok(())
    }

    /// Read the current WALE omega field. If WALE is disabled, returns the
    /// uniform base relaxation field for diagnostics.
    pub fn try_read_wale_omega(
        &self,
        fields: &GpuFields,
        base_omega: f32,
    ) -> Result<Vec<f32>, GpuError> {
        self.try_flush(fields)?;
        let Some(wale) = &fields.wale else {
            return Ok(vec![base_omega; fields.n]);
        };
        let bytes = (fields.n * 4) as u64;
        self.ensure_staging(fields, bytes);
        let staging_ref = fields.staging.borrow();
        let staging = &staging_ref.as_ref().expect("staging buffer exists").buffer;
        let mut enc = self
            .ctx
            .device
            .create_command_encoder(&wgpu::CommandEncoderDescriptor {
                label: Some("read-wale-omega"),
            });
        enc.copy_buffer_to_buffer(&wale.omega, 0, staging, 0, bytes);
        self.submit(enc.finish());
        let raw = self.map_staging(staging, bytes)?;
        Ok(bytemuck::cast_slice(&raw).to_vec())
    }

    /// Read populations, moments, and probe force through one copy encoder
    /// and one staging-buffer map for `GpuSolver::sync`.
    pub fn try_read_sync(
        &self,
        fields: &GpuFields,
        out: &mut HostMoments<f32>,
    ) -> Result<(Arc<Vec<f32>>, [f32; 3]), GpuError> {
        self.try_flush(fields)?;
        {
            let mut st = fields.state.borrow_mut();
            let cur = st.cur;
            st.ops.push(Op::Moments { bg: cur });
        }
        self.try_flush(fields)?;

        let n = fields.n;
        let fbytes = L::Q as u64 * n as u64 * fields.element_bytes;
        let plane = (n * 4) as u64;
        let moments_offset = fbytes;
        let probe_offset = moments_offset + 4 * plane;
        let bytes = probe_offset + 12;
        self.ensure_staging(fields, bytes);
        let staging_ref = fields.staging.borrow();
        let staging = &staging_ref.as_ref().expect("staging buffer exists").buffer;
        let mut enc = self
            .ctx
            .device
            .create_command_encoder(&wgpu::CommandEncoderDescriptor {
                label: Some("read-sync"),
            });
        enc.copy_buffer_to_buffer(&fields.f[fields.cur()], 0, staging, 0, fbytes);
        enc.copy_buffer_to_buffer(&fields.rho, 0, staging, moments_offset, plane);
        enc.copy_buffer_to_buffer(&fields.ux, 0, staging, moments_offset + plane, plane);
        enc.copy_buffer_to_buffer(&fields.uy, 0, staging, moments_offset + 2 * plane, plane);
        enc.copy_buffer_to_buffer(&fields.uz, 0, staging, moments_offset + 3 * plane, plane);
        enc.copy_buffer_to_buffer(&fields.probe_acc, 0, staging, probe_offset, 12);
        self.submit(enc.finish());
        let raw = self.map_staging(staging, bytes)?;

        let f_count = L::Q * n;
        let f = Arc::new(decode_storage(self.cfg.storage, &raw[..fbytes as usize]));
        let moments: &[f32] =
            bytemuck::cast_slice(&raw[moments_offset as usize..probe_offset as usize]);
        out.rho.clear();
        out.rho.extend_from_slice(&moments[..n]);
        out.ux.clear();
        out.ux.extend_from_slice(&moments[n..2 * n]);
        out.uy.clear();
        out.uy.extend_from_slice(&moments[2 * n..3 * n]);
        out.uz.clear();
        out.uz.extend_from_slice(&moments[3 * n..4 * n]);
        let probe: &[f32] = bytemuck::cast_slice(&raw[probe_offset as usize..bytes as usize]);

        let mut st = fields.state.borrow_mut();
        let generation = st.generation;
        debug_assert_eq!(f.len(), f_count);
        st.f_cache = Some((generation, Arc::clone(&f)));
        Ok((f, [probe[0], probe[1], probe[2]]))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::backend::CpuScalar;
    use crate::halo::LocalPeriodic;
    use crate::lattice::{D2Q9, D3Q19, D3Q27};
    use crate::solver::{GlobalSpec, Solver};
    use std::f64::consts::PI;
    use std::sync::OnceLock;

    fn ctx() -> Arc<GpuContext> {
        static CTX: OnceLock<Arc<GpuContext>> = OnceLock::new();
        CTX.get_or_init(|| GpuContext::new().expect("GPU backend tests require a GPU adapter"))
            .clone()
    }

    #[test]
    fn submit_chunk_calibration_targets_middle_of_window() {
        assert_eq!(
            calibrate_submit_chunk(200, Duration::from_millis(147), 200),
            200
        );
        assert_eq!(
            calibrate_submit_chunk(200, Duration::from_millis(2_000), 200),
            18
        );
        assert_eq!(
            calibrate_submit_chunk(10, Duration::from_millis(10), 200),
            175
        );
        assert_eq!(calibrate_submit_chunk(0, Duration::ZERO, 300), 200);
    }

    #[test]
    fn poll_failure_is_an_error_value() {
        let err = GpuError::Poll(wgpu::PollError::Timeout.to_string());
        assert!(err.to_string().contains("GPU device poll failed"));
    }

    #[test]
    fn resource_plan_rejects_qn_overflow() {
        let err = GpuResourcePlan::for_grid::<D3Q19>(15_050, 15_050, 1, 4).unwrap_err();
        assert!(err.to_string().contains("Q*n"));
        assert!(err.to_string().contains("u32::MAX"));
    }

    #[test]
    fn resource_plan_rejects_storage_binding_limit() {
        let plan = GpuResourcePlan::for_grid::<D2Q9>(64, 64, 1, 4).unwrap();
        let limits = wgpu::Limits {
            max_storage_buffer_binding_size: 1024,
            max_buffer_size: u64::MAX,
            ..wgpu::Limits::default()
        };
        let err = plan.validate(&limits).unwrap_err();
        assert!(err.to_string().contains("max_storage_buffer_binding_size"));
    }

    #[test]
    fn resource_plan_rejects_dispatch_limit() {
        let plan = GpuResourcePlan::for_grid::<D2Q9>(16_777_216, 1, 1, 4).unwrap();
        let limits = wgpu::Limits {
            max_storage_buffer_binding_size: u32::MAX,
            max_buffer_size: u64::MAX,
            max_compute_workgroups_per_dimension: 65_535,
            ..wgpu::Limits::default()
        };
        let err = plan.validate(&limits).unwrap_err();
        assert!(err.to_string().contains("workgroup count"));
    }

    fn omega_from_nu(nu: f64) -> f64 {
        1.0 / (3.0 * nu + 0.5)
    }

    fn cumulant_tgv_spec<L: Lattice>(n: usize) -> GlobalSpec<f32> {
        let nu = 0.02;
        let _ = L::Q;
        GlobalSpec {
            dims: [n, n, n],
            nu,
            collision: CollisionKind::Cumulant {
                omega_shear: omega_from_nu(nu),
            },
            periodic: [true, true, true],
            ..Default::default()
        }
    }

    fn init_cumulant_tgv<L, B>(s: &mut Solver<L, f32, B, LocalPeriodic>, n: usize)
    where
        L: Lattice,
        B: Backend<L, f32>,
    {
        let u0 = 1.28e-4 / n as f64;
        s.init_with(move |x, y, z| {
            let k = 2.0 * PI / n as f64;
            let (xf, yf, zf) = (k * x as f64, k * y as f64, k * z as f64);
            let p =
                u0 * u0 / 16.0 * (((2.0 * xf).cos() + (2.0 * yf).cos()) * ((2.0 * zf).cos() + 2.0));
            (
                (1.0 + 3.0 * p) as f32,
                [
                    (u0 * xf.sin() * yf.cos() * zf.cos()) as f32,
                    (-u0 * xf.cos() * yf.sin() * zf.cos()) as f32,
                    0.0,
                ],
            )
        });
    }

    fn max_delta(a: &[f32], b: &[f32]) -> f64 {
        a.iter()
            .zip(b)
            .map(|(x, y)| (*x as f64 - *y as f64).abs())
            .fold(0.0, f64::max)
    }

    fn cumulant_gpu_cpu_delta<L: Lattice>() -> f64 {
        let n = 8;
        let spec = cumulant_tgv_spec::<L>(n);
        let mut cpu: Solver<L, f32, CpuScalar, LocalPeriodic> = Solver::new(
            &spec,
            &[],
            &[],
            [1, 1, 1],
            CpuScalar::default(),
            LocalPeriodic,
        );
        let mut gpu: Solver<L, f32, WgpuBackend<L>, LocalPeriodic> = Solver::new(
            &spec,
            &[],
            &[],
            [1, 1, 1],
            WgpuBackend::<L>::new(ctx()),
            LocalPeriodic,
        );
        init_cumulant_tgv(&mut cpu, n);
        init_cumulant_tgv(&mut gpu, n);
        cpu.run(200);
        gpu.run(200);
        max_delta(&cpu.gather_rho(), &gpu.gather_rho())
            .max(max_delta(&cpu.gather_ux(), &gpu.gather_ux()))
            .max(max_delta(&cpu.gather_uy(), &gpu.gather_uy()))
            .max(max_delta(&cpu.gather_uz(), &gpu.gather_uz()))
    }

    #[test]
    fn cumulant_gpu_matches_cpu_measured_tgv3d_tolerance() {
        let d3q19 = cumulant_gpu_cpu_delta::<D3Q19>();
        let d3q27 = cumulant_gpu_cpu_delta::<D3Q27>();
        eprintln!("cumulant GPU vs CPU 200-step TGV3D: D3Q19={d3q19:e} D3Q27={d3q27:e}");
        // Measured 2026-07-06 on the Stage-3 WGSL path:
        // D3Q19 = 1.5497207641601563e-6, D3Q27 = 3.4570693969726563e-6
        // over a 200-step f32 TGV3D. Gates are measured * 10 headroom.
        assert!(d3q19 <= 1.6e-5, "D3Q19 cumulant GPU delta {d3q19:e}");
        assert!(d3q27 <= 3.5e-5, "D3Q27 cumulant GPU delta {d3q27:e}");
    }
}
