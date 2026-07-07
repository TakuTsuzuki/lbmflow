//! Radar #22 / pitfall #18: D3Q19 vs D3Q27 rotated-flow anisotropy.
//!
//! The continuum diffusion operator is rotationally invariant. A lattice
//! quadrature only inherits that invariance up to the tensor moments it
//! integrates exactly. D3Q19 and D3Q27 both reproduce the required second
//! moment and the fourth moment,
//!
//!   sum_q w_q c_qi c_qj = cs^2 delta_ij
//!   sum_q w_q c_qi c_qj c_qk c_ql
//!     = cs^4 (delta_ij delta_kl + delta_ik delta_jl + delta_il delta_jk),
//!
//! so the Navier-Stokes viscous stress is isotropic at the formal order. The
//! difference appears in the next truncation layer: D3Q19 omits the body
//! diagonals and cannot make every sixth-order tensor component isotropic
//! (for example the c_x^2 c_y^2 c_z^2 moment is identically zero), while
//! D3Q27 includes the full tensor-product stencil. In rotated shear, those
//! missing sixth-order components can show up as axis-aligned dissipation or
//! weak secondary flow. This is the standard D3Q19-vs-D3Q27 anisotropy
//! warning discussed in Krueger et al. and by White & Chong, JCP 230 (2011).
//!
//! A literal 45-degree rotation of the fundamental Fourier mode has wavevector
//! components k/sqrt(2), which is not periodic on an N x N x N lattice. To keep
//! the periodic native solver honest, this test uses the closest exact
//! diagonal Fourier construction: the rotated coordinates are
//!   xi   = k (x + y)
//!   eta  = k (-x + y)
//!   zeta = k z
//! with k = 2 pi / N. That is the standard Rz(45 deg) TGV coordinate frame
//! scaled by sqrt(2) in the rotated x-y directions so every phase advances by
//! an integer multiple of 2 pi across the periodic box. The analytic reference
//! is the linear diffusion limit: each product mode decays as
//! exp(-nu |K|^2 t), and here |K|^2 = |grad xi|^2 + |grad eta|^2
//! + |grad zeta|^2 = 5 k^2. The run time follows the T15.4 short-run
//! convention, t = 0.1 / (nu k^2), with diffusive-scaling velocity
//! u0 = 1.28e-4 / N so vortex stretching is below the spatial-error floor.

use lbm_core::prelude::*;
use std::f64::consts::{FRAC_1_SQRT_2, PI};

mod common;

const N: usize = 32;
const NU: f64 = 0.02;
const U0_COEF: f64 = 1.28e-4;

type CpuPeriodic<L> = Solver<L, f64, CpuScalar, LocalPeriodic>;

fn rotated_tgv_velocity(u0: f64, k: f64, x: usize, y: usize, z: usize) -> [f64; 3] {
    let xf = x as f64;
    let yf = y as f64;
    let zeta = k * z as f64;
    let xi = k * (xf + yf);
    let eta = k * (-xf + yf);

    let u_xi = u0 * xi.sin() * eta.cos() * zeta.cos();
    let u_eta = -u0 * xi.cos() * eta.sin() * zeta.cos();

    [
        FRAC_1_SQRT_2 * (u_xi - u_eta),
        FRAC_1_SQRT_2 * (u_xi + u_eta),
        0.0,
    ]
}

fn rotated_tgv_density(u0: f64, k: f64, x: usize, y: usize, z: usize) -> f64 {
    let xf = x as f64;
    let yf = y as f64;
    let zeta = k * z as f64;
    let xi = k * (xf + yf);
    let eta = k * (-xf + yf);

    // Classic TGV pressure in the rotated coordinate frame; rho = 1 + p/cs^2
    // and cs^2 = 1/3. The O(u0^2) pressure-consistent initialization removes
    // an avoidable acoustic transient without changing the low-Mach reference.
    let p = u0 * u0 / 16.0 * ((2.0 * xi).cos() + (2.0 * eta).cos()) * ((2.0 * zeta).cos() + 2.0);
    1.0 + 3.0 * p
}

fn run_rotated_tgv<L: Lattice>() -> (f64, f64) {
    let n = N;
    let u0 = U0_COEF / n as f64;
    let k = 2.0 * PI / n as f64;
    let spec = GlobalSpec::<f64> {
        dims: [n, n, n],
        nu: NU,
        periodic: [true, true, true],
        ..Default::default()
    };
    let mut s: CpuPeriodic<L> = Solver::new(
        &spec,
        &[],
        &[],
        [1, 1, 1],
        CpuScalar::default(),
        LocalPeriodic,
    );
    s.init_with(move |x, y, z| {
        (
            rotated_tgv_density(u0, k, x, y, z),
            rotated_tgv_velocity(u0, k, x, y, z),
        )
    });

    let steps = (0.1 / (NU * k * k)).round() as usize;
    s.run(steps);

    let decay = (-5.0 * NU * k * k * steps as f64).exp();
    let ux = s.gather_ux();
    let uy = s.gather_uy();
    let uz = s.gather_uz();
    let mut actual = Vec::with_capacity(3 * n * n * n);
    let mut reference = Vec::with_capacity(3 * n * n * n);
    let mut max_uz = 0.0f64;
    for z in 0..n {
        for y in 0..n {
            for x in 0..n {
                let i = (z * n + y) * n + x;
                let v = rotated_tgv_velocity(u0, k, x, y, z);
                actual.extend([ux[i], uy[i], uz[i]]);
                reference.extend([v[0] * decay, v[1] * decay, v[2] * decay]);
                max_uz = max_uz.max(uz[i].abs());
            }
        }
    }

    (common::metrics::l2_rel(&actual, &reference), max_uz / u0)
}

#[test]
fn d3q19_vs_d3q27_rotated_tgv_anisotropy_discriminator() {
    let (err19, uz19_rel) = run_rotated_tgv::<D3Q19>();
    let (err27, uz27_rel) = run_rotated_tgv::<D3Q27>();
    let ratio = err19 / err27;

    println!(
        "rotated TGV anisotropy N={N}: D3Q19 L2rel={err19:.9e}, D3Q27 L2rel={err27:.9e}, ratio={ratio:.6}, max|uz|/u0 D3Q19={uz19_rel:.3e}, D3Q27={uz27_rel:.3e}"
    );

    assert!(
        err27 <= 5.0e-3,
        "D3Q27 rotated TGV L2rel={err27:.9e} exceeds 5e-3; normalization=||u_D3Q27-u_diffusion||2/||u_diffusion||2"
    );
    assert!(
        err19 <= 5.0 * err27,
        "D3Q19 rotated TGV L2rel={err19:.9e} exceeds 5x D3Q27={err27:.9e}; ratio={ratio:.6}; normalization=||u-u_diffusion||2/||u_diffusion||2"
    );
    assert!(
        uz19_rel <= 1.0e-4 && uz27_rel <= 1.0e-4,
        "rotated in-plane TGV should not grow material out-of-plane flow above the low-Mach numerical floor: max|uz|/u0 D3Q19={uz19_rel:.3e}, D3Q27={uz27_rel:.3e}, band=1e-4"
    );
}
