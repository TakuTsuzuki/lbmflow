//! Dev probe: pressure-pressure driven channel sanity (flow must go from
//! high density to low density, magnitude near Poiseuille).

use lbm_core::prelude::*;

fn main() {
    let (nx, ny) = (64, 34); // H = 32
    let drho = 2e-3;
    let nu = 0.05;
    let mut sim: Simulation<f64> = SimConfig {
        nx,
        ny,
        nu,
        edges: Edges {
            left: EdgeBC::PressureOutlet {
                rho: 1.0 + drho / 2.0,
            },
            right: EdgeBC::PressureOutlet {
                rho: 1.0 - drho / 2.0,
            },
            bottom: EdgeBC::BounceBack,
            top: EdgeBC::BounceBack,
        },
        ..Default::default()
    }
    .build()
    .unwrap();
    sim.run(20_000);
    let h = (ny - 2) as f64;
    let l = (nx - 1) as f64;
    let cs2 = 1.0 / 3.0;
    let dpdx = cs2 * drho / l;
    let umax_theory = dpdx * h * h / (8.0 * nu);
    let mid = ny / 2;
    let umax_sim = sim.ux(nx / 2, mid);
    println!(
        "u_center(sim) = {umax_sim:.6e}  u_max(theory) = {umax_theory:.6e}  ratio = {:.4}",
        umax_sim / umax_theory
    );
    // flow direction check at quarter points
    println!(
        "ux at x=8: {:.3e}, x=32: {:.3e}, x=56: {:.3e}",
        sim.ux(8, mid),
        sim.ux(32, mid),
        sim.ux(56, mid)
    );
    let mass_ok = sim.total_mass().is_finite();
    println!("finite = {mass_ok}");
}
