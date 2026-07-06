use lbm_core::backend::CpuScalar;
use lbm_core::compat::prelude::{EdgeBC, Edges, SimConfig, Simulation};
use lbm_core::halo::LocalPeriodic;
use lbm_core::lattice::D2Q9;
use lbm_core::prelude::{GlobalSpec, Solver};

const N: usize = 1024;
const EXTRA_MASS: f64 = 0.03125;

fn expected_mass() -> f64 {
    (N * N) as f64 + EXTRA_MASS
}

fn assert_f64_mass_preserves_sub_ulp_bump(precise: f64, source_compatible: f64) {
    let base = (N * N) as f64;
    let quantized = expected_mass() as f32 as f64;

    assert_eq!(quantized, base);
    assert_eq!(source_compatible, quantized);
    assert!(
        (precise - expected_mass()).abs() <= 1.0e-9,
        "precise mass {precise:.16e} should retain the sub-f32-ULP bump near {expected:.16e}",
        expected = expected_mass()
    );
    assert!(
        (precise - source_compatible - EXTRA_MASS).abs() <= 1.0e-9,
        "f64 diagnostic should preserve the bump that total_mass() quantizes away: precise={precise:.16e}, total_mass={source_compatible:.16e}"
    );
}

fn init_bump(x: usize, y: usize) -> (f32, [f32; 3]) {
    let rho = if (x, y) == (N / 2, N / 2) {
        1.0 + EXTRA_MASS as f32
    } else {
        1.0
    };
    (rho, [0.0, 0.0, 0.0])
}

#[test]
fn solver_total_mass_f64_avoids_f32_scalar_quantization_at_1024sq() {
    let spec = GlobalSpec::<f32> {
        dims: [N, N, 1],
        periodic: [true, true, false],
        ..Default::default()
    };
    let mut solver: Solver<D2Q9, f32, CpuScalar, LocalPeriodic> = Solver::new(
        &spec,
        &[],
        &[],
        [1, 1, 1],
        CpuScalar::default(),
        LocalPeriodic,
    );
    solver.init_with(|x, y, _| init_bump(x, y));

    let precise = solver.total_mass_f64();
    let source_compatible = solver.total_mass() as f64;

    assert_f64_mass_preserves_sub_ulp_bump(precise, source_compatible);
}

#[test]
fn compat_total_mass_f64_avoids_f32_scalar_quantization_at_1024sq() {
    let mut sim: Simulation<f32> = SimConfig {
        nx: N,
        ny: N,
        edges: Edges {
            left: EdgeBC::Periodic,
            right: EdgeBC::Periodic,
            bottom: EdgeBC::Periodic,
            top: EdgeBC::Periodic,
        },
        ..Default::default()
    }
    .build()
    .unwrap();
    sim.init_with(|x, y| {
        let (rho, u) = init_bump(x, y);
        (rho, u[0], u[1])
    });

    let precise = sim.total_mass_f64();
    let source_compatible = sim.total_mass() as f64;

    assert_f64_mass_preserves_sub_ulp_bump(precise, source_compatible);
}
