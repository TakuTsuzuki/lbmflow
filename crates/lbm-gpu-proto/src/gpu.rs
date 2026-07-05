//! Thin wgpu wrapper: device setup + the ping-pong D2Q9 pipeline.

use crate::hostinit;

/// Shared adapter/device/queue.
pub struct GpuContext {
    pub device: wgpu::Device,
    pub queue: wgpu::Queue,
    pub adapter_info: wgpu::AdapterInfo,
}

impl GpuContext {
    /// Returns `None` (with a printed reason) when no usable GPU exists —
    /// itself a valid evaluation result.
    pub fn new() -> Option<Self> {
        let instance = wgpu::Instance::default();
        let adapter = match pollster::block_on(instance.request_adapter(
            &wgpu::RequestAdapterOptions {
                power_preference: wgpu::PowerPreference::HighPerformance,
                ..Default::default()
            },
        )) {
            Ok(a) => a,
            Err(e) => {
                eprintln!("no GPU adapter available: {e}");
                return None;
            }
        };
        let adapter_info = adapter.get_info();
        // 2048^2 needs a single 151 MB storage binding (9 * n * 4 B); the
        // WebGPU *default* limits cap bindings at 128 MiB, so request the
        // adapter's actual limits for buffer size / binding size.
        let al = adapter.limits();
        let mut limits = wgpu::Limits::default();
        limits.max_storage_buffer_binding_size = al.max_storage_buffer_binding_size;
        limits.max_buffer_size = al.max_buffer_size;
        let (device, queue) = match pollster::block_on(adapter.request_device(
            &wgpu::DeviceDescriptor {
                label: Some("lbm-gpu-proto"),
                required_features: wgpu::Features::empty(),
                required_limits: limits,
                memory_hints: wgpu::MemoryHints::Performance,
                trace: wgpu::Trace::Off,
            },
        )) {
            Ok(dq) => dq,
            Err(e) => {
                eprintln!("request_device failed: {e}");
                return None;
            }
        };
        Some(Self {
            device,
            queue,
            adapter_info,
        })
    }
}

/// One D2Q9 simulation instance on the GPU: two population buffers
/// (ping-pong), the fused stream+collide pipeline and a moments/readback path.
pub struct GpuLbm<'a> {
    ctx: &'a GpuContext,
    nx: u32,
    ny: u32,
    wg: (u32, u32),
    step_pipeline: wgpu::ComputePipeline,
    moments_pipeline: wgpu::ComputePipeline,
    /// step_bg[p]: read f[p], write f[1-p].
    step_bg: [wgpu::BindGroup; 2],
    /// moments_bg[p]: read f[p], write vel.
    moments_bg: [wgpu::BindGroup; 2],
    f: [wgpu::Buffer; 2],
    vel: wgpu::Buffer,
    staging: wgpu::Buffer,
    /// Which f buffer holds the current state.
    cur: usize,
}

impl<'a> GpuLbm<'a> {
    pub fn new(ctx: &'a GpuContext, nx: u32, ny: u32, nu: f64, wg: (u32, u32)) -> Self {
        let device = &ctx.device;
        let n = (nx as u64) * (ny as u64);
        let f_bytes = n * 9 * 4;
        let vel_bytes = n * 2 * 4;

        let source = include_str!("shader.wgsl")
            .replace("__WGX__", &wg.0.to_string())
            .replace("__WGY__", &wg.1.to_string());
        let module = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("lbm-d2q9"),
            source: wgpu::ShaderSource::Wgsl(source.into()),
        });
        let mk_pipeline = |entry: &str| {
            device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
                label: Some(entry),
                layout: None, // auto layout: each entry point gets exactly its bindings
                module: &module,
                entry_point: Some(entry),
                compilation_options: Default::default(),
                cache: None,
            })
        };
        let step_pipeline = mk_pipeline("step");
        let moments_pipeline = mk_pipeline("moments");

        let buf = |label: &str, size: u64, usage: wgpu::BufferUsages| {
            device.create_buffer(&wgpu::BufferDescriptor {
                label: Some(label),
                size,
                usage,
                mapped_at_creation: false,
            })
        };
        use wgpu::BufferUsages as U;
        let f = [
            buf("f0", f_bytes, U::STORAGE | U::COPY_DST),
            buf("f1", f_bytes, U::STORAGE | U::COPY_DST),
        ];
        let vel = buf("vel", vel_bytes, U::STORAGE | U::COPY_SRC);
        let staging = buf("staging", vel_bytes, U::MAP_READ | U::COPY_DST);
        let params = buf("params", 16, U::UNIFORM | U::COPY_DST);

        let (omega_p, omega_m) = hostinit::omegas_f32(nu);
        let words: [u32; 4] = [nx, ny, omega_p.to_bits(), omega_m.to_bits()];
        ctx.queue
            .write_buffer(&params, 0, bytemuck::cast_slice(&words));

        let step_layout = step_pipeline.get_bind_group_layout(0);
        let step_bg = [0usize, 1].map(|p| {
            device.create_bind_group(&wgpu::BindGroupDescriptor {
                label: Some("step"),
                layout: &step_layout,
                entries: &[
                    wgpu::BindGroupEntry {
                        binding: 0,
                        resource: params.as_entire_binding(),
                    },
                    wgpu::BindGroupEntry {
                        binding: 1,
                        resource: f[p].as_entire_binding(),
                    },
                    wgpu::BindGroupEntry {
                        binding: 2,
                        resource: f[1 - p].as_entire_binding(),
                    },
                ],
            })
        });
        let moments_layout = moments_pipeline.get_bind_group_layout(0);
        let moments_bg = [0usize, 1].map(|p| {
            device.create_bind_group(&wgpu::BindGroupDescriptor {
                label: Some("moments"),
                layout: &moments_layout,
                entries: &[
                    wgpu::BindGroupEntry {
                        binding: 0,
                        resource: params.as_entire_binding(),
                    },
                    wgpu::BindGroupEntry {
                        binding: 1,
                        resource: f[p].as_entire_binding(),
                    },
                    wgpu::BindGroupEntry {
                        binding: 3,
                        resource: vel.as_entire_binding(),
                    },
                ],
            })
        });

        Self {
            ctx,
            nx,
            ny,
            wg,
            step_pipeline,
            moments_pipeline,
            step_bg,
            moments_bg,
            f,
            vel,
            staging,
            cur: 0,
        }
    }

    /// Upload an initial population state (deviation storage, SoA `f[q*n+i]`).
    pub fn upload(&self, f_soa: &[f32]) {
        assert_eq!(f_soa.len() as u64, self.n_cells() * 9);
        self.ctx
            .queue
            .write_buffer(&self.f[self.cur], 0, bytemuck::cast_slice(f_soa));
    }

    fn n_cells(&self) -> u64 {
        self.nx as u64 * self.ny as u64
    }

    fn workgroups(&self) -> (u32, u32) {
        (self.nx.div_ceil(self.wg.0), self.ny.div_ceil(self.wg.1))
    }

    /// Advance `steps` steps: one dispatch per step, `chunk` dispatches per
    /// command buffer, single wait at the end. When `wait_each_submit` is set
    /// the CPU blocks after every submit (worst-case interactive pattern).
    pub fn run_opts(&mut self, steps: usize, chunk: usize, wait_each_submit: bool) {
        let (gx, gy) = self.workgroups();
        let mut done = 0;
        while done < steps {
            let k = (steps - done).min(chunk.max(1));
            let mut enc = self
                .ctx
                .device
                .create_command_encoder(&wgpu::CommandEncoderDescriptor { label: Some("run") });
            {
                let mut pass = enc.begin_compute_pass(&wgpu::ComputePassDescriptor {
                    label: Some("lbm-steps"),
                    timestamp_writes: None,
                });
                pass.set_pipeline(&self.step_pipeline);
                for _ in 0..k {
                    pass.set_bind_group(0, &self.step_bg[self.cur], &[]);
                    pass.dispatch_workgroups(gx, gy, 1);
                    self.cur ^= 1;
                }
            }
            self.ctx.queue.submit(Some(enc.finish()));
            if wait_each_submit {
                self.wait();
            }
            done += k;
        }
        self.wait();
    }

    /// `run_opts` with a single final wait (the benchmark's submit→wait path).
    pub fn run(&mut self, steps: usize, chunk: usize) {
        self.run_opts(steps, chunk, false);
    }

    fn wait(&self) {
        self.ctx
            .device
            .poll(wgpu::PollType::Wait)
            .expect("device poll");
    }

    /// Compute moments on the GPU and read back interleaved `[ux, uy]` pairs
    /// (dispatch + copy + map, blocking).
    pub fn velocity(&self) -> Vec<f32> {
        let (gx, gy) = self.workgroups();
        let mut enc = self
            .ctx
            .device
            .create_command_encoder(&wgpu::CommandEncoderDescriptor {
                label: Some("moments"),
            });
        {
            let mut pass = enc.begin_compute_pass(&wgpu::ComputePassDescriptor {
                label: Some("moments"),
                timestamp_writes: None,
            });
            pass.set_pipeline(&self.moments_pipeline);
            pass.set_bind_group(0, &self.moments_bg[self.cur], &[]);
            pass.dispatch_workgroups(gx, gy, 1);
        }
        enc.copy_buffer_to_buffer(&self.vel, 0, &self.staging, 0, self.n_cells() * 2 * 4);
        self.ctx.queue.submit(Some(enc.finish()));

        let slice = self.staging.slice(..);
        let (tx, rx) = futures_intrusive::channel::shared::oneshot_channel();
        slice.map_async(wgpu::MapMode::Read, move |res| {
            let _ = tx.send(res);
        });
        self.wait();
        pollster::block_on(rx.receive())
            .expect("map_async dropped")
            .expect("buffer map failed");
        let out: Vec<f32> = bytemuck::cast_slice(&slice.get_mapped_range()).to_vec();
        self.staging.unmap();
        out
    }
}
