#![cfg(feature = "gpu")]
//! T14-WALE: on-device GPU WALE omega pass vs CPU WALE driver.

use lbm_core::lattice::D3Q19;
use lbm_core::prelude::*;
use std::sync::{Arc, OnceLock};

type Cpu3 = Solver<D3Q19, f32, CpuScalar, LocalPeriodic>;

fn ctx() -> Arc<GpuContext> {
    static CTX: OnceLock<Arc<GpuContext>> = OnceLock::new();
    CTX.get_or_init(|| GpuContext::new().expect("T14-WALE requires a GPU adapter"))
        .clone()
}

fn velocity_l2_rel(cpu: &Cpu3, gpu: &mut GpuSolver<D3Q19>) -> f64 {
    let (ux, uy, uz) = (cpu.gather_ux(), cpu.gather_uy(), cpu.gather_uz());
    let (gx, gy, gz) = (gpu.gather_ux(), gpu.gather_uy(), gpu.gather_uz());
    let mut num = 0.0f64;
    let mut den = 0.0f64;
    for i in 0..ux.len() {
        let dx = ux[i] as f64 - gx[i] as f64;
        let dy = uy[i] as f64 - gy[i] as f64;
        let dz = uz[i] as f64 - gz[i] as f64;
        num += dx * dx + dy * dy + dz * dz;
        den += (ux[i] as f64).powi(2) + (uy[i] as f64).powi(2) + (uz[i] as f64).powi(2);
    }
    num.sqrt() / den.sqrt().max(1.0e-12)
}

fn omega_from_les(nu: f32, les: &WaleLes<f32>) -> Vec<f32> {
    les.nu_t()
        .iter()
        .map(|&nut| 1.0 / (3.0 * (nu + nut) + 0.5))
        .collect()
}

fn omega_max_abs_rel(cpu: &[f32], gpu: &[f32]) -> f64 {
    let mut d = 0.0f64;
    let mut scale = 0.0f64;
    for (&a, &b) in cpu.iter().zip(gpu) {
        d = d.max((a as f64 - b as f64).abs());
        scale = scale.max((a as f64).abs());
    }
    d / scale.max(1.0e-12)
}

fn run_wale_pair(cpu: &mut Cpu3, gpu: &mut GpuSolver<D3Q19>, steps: usize) -> Vec<f32> {
    let nu = (cpu.tau() as f32 - 0.5) / 3.0;
    let mut les = WaleLes::<f32>::new();
    gpu.set_wale(true);
    gpu.set_submit_chunk(1);
    for _ in 0..steps {
        les.update(cpu);
        cpu.run(1);
        gpu.run(1);
    }
    omega_from_les(nu, &les)
}

#[test]
fn t14_wale_tgv3d_gpu_equiv() {
    let n = 32usize;
    let nu = 0.02f64;
    let spec = GlobalSpec::<f32> {
        dims: [n, n, n],
        nu,
        collision: CollisionKind::Trt {
            magic: CollisionKind::MAGIC_STD,
        },
        periodic: [true, true, true],
        ..Default::default()
    };
    let mut cpu = Cpu3::new(
        &spec,
        &[],
        &[],
        [1, 1, 1],
        CpuScalar::default(),
        LocalPeriodic,
    );
    let mut gpu = GpuSolver::<D3Q19>::new(&spec, &[], &[], ctx());
    let u0 = 0.03f64;
    let k = 2.0 * std::f64::consts::PI / n as f64;
    let init = move |x: usize, y: usize, z: usize| {
        let (xf, yf, zf) = (k * x as f64, k * y as f64, k * z as f64);
        (
            1.0f32,
            [
                (u0 * xf.sin() * yf.cos() * zf.cos()) as f32,
                (-u0 * xf.cos() * yf.sin() * zf.cos()) as f32,
                0.0,
            ],
        )
    };
    cpu.init_with(init);
    gpu.init_with(init);

    let cpu_omega = run_wale_pair(&mut cpu, &mut gpu, 200);
    let u_rel = velocity_l2_rel(&cpu, &mut gpu);
    let omega_rel = omega_max_abs_rel(&cpu_omega, &gpu.gather_wale_omega());
    eprintln!("T14-WALE TGV3D: u_l2_rel={u_rel:.3e} omega_max_rel={omega_rel:.3e}");
    assert!(u_rel <= 1.0e-5, "TGV3D WALE u L2rel {u_rel:e}");
    assert!(
        omega_rel <= 1.0e-6,
        "TGV3D WALE omega max relative diff {omega_rel:e}"
    );
}

#[test]
fn t14_wale_cavity3d_gpu_equiv() {
    let n = 24usize;
    let mut walls = WallSpec::default();
    for face in Face::ALL {
        walls.is_wall[face.index()] = true;
    }
    walls.u[Face::YPos.index()] = [0.06, 0.0, 0.0];
    let spec = GlobalSpec::<f32> {
        dims: [n, n, n],
        nu: 0.03,
        collision: CollisionKind::Trt {
            magic: CollisionKind::MAGIC_STD,
        },
        periodic: [false, false, false],
        ..Default::default()
    };
    let (solid, wall_u) = build_wall_rims(3, spec.dims, &walls);
    let mut cpu = Cpu3::new(
        &spec,
        &solid,
        &wall_u,
        [1, 1, 1],
        CpuScalar::default(),
        LocalPeriodic,
    );
    let mut gpu = GpuSolver::<D3Q19>::new(&spec, &solid, &wall_u, ctx());

    let cpu_omega = run_wale_pair(&mut cpu, &mut gpu, 200);
    let u_rel = velocity_l2_rel(&cpu, &mut gpu);
    let omega_rel = omega_max_abs_rel(&cpu_omega, &gpu.gather_wale_omega());
    eprintln!("T14-WALE cavity3D: u_l2_rel={u_rel:.3e} omega_max_rel={omega_rel:.3e}");
    assert!(u_rel <= 1.0e-5, "cavity3D WALE u L2rel {u_rel:e}");
    assert!(
        omega_rel <= 1.0e-6,
        "cavity3D WALE omega max relative diff {omega_rel:e}"
    );
}

#[test]
fn t14_wale_uniform_periodic_keeps_base_omega() {
    let dims = [16usize, 12usize, 10usize];
    let nu = 0.025f64;
    let spec = GlobalSpec::<f32> {
        dims,
        nu,
        collision: CollisionKind::Trt {
            magic: CollisionKind::MAGIC_STD,
        },
        periodic: [true, true, true],
        ..Default::default()
    };
    let mut gpu = GpuSolver::<D3Q19>::new(&spec, &[], &[], ctx());
    gpu.init_with(|_, _, _| (1.0, [0.02, -0.01, 0.005]));
    gpu.set_wale(true);
    gpu.run(3);
    let base = (1.0 / (3.0 * nu + 0.5)) as f32;
    let omega = gpu.gather_wale_omega();
    let max_abs = omega
        .iter()
        .map(|&w| (w - base).abs())
        .fold(0.0f32, f32::max);
    assert!(
        max_abs == 0.0,
        "uniform periodic WALE should keep base omega exactly, max_abs={max_abs:e}"
    );
}
