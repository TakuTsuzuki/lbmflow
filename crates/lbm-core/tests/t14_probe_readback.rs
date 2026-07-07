#![cfg(feature = "gpu")]
//! T14 B-2 probe readback contract: CPU and GPU expose the most recent
//! momentum-exchange probe through the same `Solver::read_probed_force` API.

use lbm_core::prelude::*;
use std::sync::{Arc, OnceLock};

type Cpu = Solver<D2Q9, f32, CpuScalar, LocalPeriodic>;
type Gpu = Solver<D2Q9, f32, WgpuBackend<D2Q9>, LocalPeriodic>;

const DIAG_TOL: f64 = 1e-4;

fn ctx() -> Arc<GpuContext> {
    static CTX: OnceLock<Arc<GpuContext>> = OnceLock::new();
    CTX.get_or_init(|| GpuContext::new().expect("T14 requires a GPU adapter"))
        .clone()
}

fn assert_probe_close(cpu: [f32; 3], gpu: [f32; 3], what: &str) {
    let scale = cpu
        .iter()
        .map(|v| (*v as f64).abs())
        .fold(1e-6_f64, f64::max);
    for c in 0..2 {
        let d = (cpu[c] as f64 - gpu[c] as f64).abs();
        let lim = DIAG_TOL * scale;
        assert!(
            d <= lim,
            "{what}: force[{c}] |delta|={d:.3e} > {lim:.3e} \
             (cpu={:.9e}, gpu={:.9e})",
            cpu[c] as f64,
            gpu[c] as f64
        );
    }
}

#[test]
fn t14_probe_force_explicit_readback_cpu_gpu() {
    let (nx, ny) = (112usize, 56usize);
    let mut walls = WallSpec::<f32>::default();
    walls.is_wall[Face::YNeg.index()] = true;
    walls.is_wall[Face::YPos.index()] = true;
    let mut faces = [FaceBC::Closed; 6];
    faces[Face::XNeg.index()] = FaceBC::Velocity {
        u: [0.055, 0.0, 0.0],
    };
    faces[Face::XPos.index()] = FaceBC::Outflow;
    let spec = GlobalSpec {
        dims: [nx, ny, 1],
        nu: 0.035,
        periodic: [false, false, false],
        faces,
        ..Default::default()
    };
    let (solid, wall_u) = build_wall_rims(2, spec.dims, &walls);
    let mut cpu = Cpu::new(
        &spec,
        &solid,
        &wall_u,
        [1, 1, 1],
        CpuScalar::default(),
        LocalPeriodic,
    );
    let mut gpu = Gpu::new(
        &spec,
        &solid,
        &wall_u,
        [1, 1, 1],
        WgpuBackend::new(ctx()),
        LocalPeriodic,
    );
    let inside = |x: usize, y: usize, _: usize| (30..42).contains(&x) && y <= 6;
    for y in 0..ny {
        for x in 0..nx {
            if inside(x, y, 0) {
                cpu.set_solid(x, y, 0);
                gpu.set_solid(x, y, 0);
            }
        }
    }
    cpu.set_force_probe(inside);
    gpu.set_force_probe(inside);

    for chunk in 1..=4 {
        cpu.run(50);
        gpu.run(50);
        let cpu_read = cpu.read_probed_force();
        let gpu_read = gpu.read_probed_force();
        eprintln!(
            "chunk {chunk} t={}: cpu=[{:.9e},{:.9e},{:.9e}] \
             gpu=[{:.9e},{:.9e},{:.9e}] delta=[{:.3e},{:.3e},{:.3e}]",
            cpu.time(),
            cpu_read[0] as f64,
            cpu_read[1] as f64,
            cpu_read[2] as f64,
            gpu_read[0] as f64,
            gpu_read[1] as f64,
            gpu_read[2] as f64,
            cpu_read[0] as f64 - gpu_read[0] as f64,
            cpu_read[1] as f64 - gpu_read[1] as f64,
            cpu_read[2] as f64 - gpu_read[2] as f64
        );
        assert_eq!(
            cpu.probed_force().map(f32::to_bits),
            cpu_read.map(f32::to_bits),
            "CPU cached/readback force differs after chunk {chunk}"
        );
        assert_eq!(
            gpu.probed_force().map(f32::to_bits),
            gpu_read.map(f32::to_bits),
            "GPU cached/readback force differs after chunk {chunk}"
        );
        assert_probe_close(cpu_read, gpu_read, "explicit probe readback");
    }
}
