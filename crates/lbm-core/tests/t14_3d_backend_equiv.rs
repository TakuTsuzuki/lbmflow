#![cfg(feature = "gpu")]
//! T14-3D: CpuScalar vs WgpuBackend on D3Q19.

use lbm_core::lattice::D3Q19;
use lbm_core::prelude::*;
use std::sync::{Arc, OnceLock};

fn ctx() -> Arc<GpuContext> {
    static CTX: OnceLock<Arc<GpuContext>> = OnceLock::new();
    CTX.get_or_init(|| GpuContext::new().expect("T14-3D requires a GPU adapter"))
        .clone()
}

type Cpu3 = Solver<D3Q19, f32, CpuScalar, LocalPeriodic>;

fn linf_rel(a: &[f32], b: &[f32], floor: f64) -> f64 {
    assert_eq!(a.len(), b.len());
    let mut d = 0.0f64;
    let mut m = 0.0f64;
    for (x, y) in a.iter().zip(b) {
        d = d.max((*x as f64 - *y as f64).abs());
        m = m.max((*x as f64).abs());
    }
    d / m.max(floor)
}

fn check(cpu: &mut Cpu3, gpu: &mut GpuSolver<D3Q19>, what: &str) {
    let t = cpu.time();
    let rho_rel = linf_rel(&cpu.gather_rho(), &gpu.gather_rho(), 1.0);
    let ux = cpu.gather_ux();
    let uy = cpu.gather_uy();
    let uz = cpu.gather_uz();
    let gx = gpu.gather_ux();
    let gy = gpu.gather_uy();
    let gz = gpu.gather_uz();
    let mut du = 0.0f64;
    let mut umax = 0.0f64;
    for i in 0..ux.len() {
        du = du
            .max((ux[i] as f64 - gx[i] as f64).abs())
            .max((uy[i] as f64 - gy[i] as f64).abs())
            .max((uz[i] as f64 - gz[i] as f64).abs());
        umax = umax.max((ux[i] as f64).hypot(uy[i] as f64).hypot(uz[i] as f64));
    }
    let u_rel = du / umax.max(1e-6);
    let mut df = 0.0f64;
    let mut fmax = 0.0f64;
    for q in 0..D3Q19::Q {
        let cf = cpu.gather_f(q);
        let gf = gpu.gather_f(q);
        for (a, b) in cf.iter().zip(&gf) {
            df = df.max((*a as f64 - *b as f64).abs());
            fmax = fmax.max((*a as f64).abs());
        }
    }
    let f_rel = df / fmax.max(1e-6);
    eprintln!("{what} t={t}: rho={rho_rel:.3e} u={u_rel:.3e} f={f_rel:.3e}");
    assert!(rho_rel <= 1e-5, "{what} rho rel {rho_rel:e}");
    assert!(u_rel <= 1e-5, "{what} u rel {u_rel:e}");
    assert!(f_rel <= 1e-4, "{what} f rel {f_rel:e}");
}

#[test]
fn t14_3d_tgv_periodic_d3q19() {
    let n = 32usize;
    let spec = GlobalSpec::<f32> {
        dims: [n, n, n],
        nu: 0.02,
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
    for _ in 0..3 {
        cpu.run(40);
        gpu.run(40);
        check(&mut cpu, &mut gpu, "TGV3D");
    }
}

#[test]
fn t14_3d_lid_cavity_d3q19() {
    let n = 24usize;
    let mut walls = WallSpec::default();
    for face in Face::ALL {
        walls.is_wall[face.index()] = true;
    }
    walls.u[Face::YPos.index()] = [0.06, 0.0, 0.0];
    let spec = GlobalSpec::<f32> {
        dims: [n, n, n],
        nu: 0.03,
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
    for _ in 0..3 {
        cpu.run(40);
        gpu.run(40);
        check(&mut cpu, &mut gpu, "cavity3D");
    }
}
