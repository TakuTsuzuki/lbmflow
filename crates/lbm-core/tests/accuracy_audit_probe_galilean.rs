//! Radar row 16: Wen-2014 co-moving-frame wall-force Galilean decider.
//!
//! This is intentionally a one-step discriminator for the momentum-exchange
//! force probe used on moving walls. It does not use implementation internals
//! as an oracle; the reference is the continuum rigid-translation limit.

use lbm_core::compat::prelude::*;

const RHO0: f64 = 1.0;
const U0: f64 = 0.05;
const TAU: f64 = 0.8;

const TRT: Collision = Collision::Trt {
    magic: Collision::MAGIC_STD,
};

fn tau_to_nu(tau: f64) -> f64 {
    (tau - 0.5) / 3.0
}

#[test]
fn co_moving_moving_walls_have_zero_physical_x_force_wen_2014() {
    // Analytic reference:
    //
    // Put the whole channel in rigid translation: periodic x, uniform
    // velocity u0 everywhere, and both y-walls moving with exactly the same
    // tangential velocity u_w = u0. The velocity gradient is zero, so the
    // Newtonian wall shear is tau_xy = mu * d(ux)/dy = 0. With no acceleration
    // and no relative wall/fluid motion, the physical x-force on the walls is
    // therefore exactly zero.
    //
    // Momentum-exchange discriminator:
    //
    // The conventional Ladd probe uses
    //
    //   F_Ladd = - sum_links c_q * (f_out + f_in)
    //
    // on fluid-solid links. Wen et al. 2014 observed that this counts the
    // momentum of mass advected with a moving boundary, so the measured force
    // changes under a Galilean boost. The invariant wall-frame form subtracts
    // the wall velocity in every link contribution:
    //
    //   F_Wen = - sum_links (c_q - u_w) * (f_out + f_in).
    //
    // In this co-moving state, the conventional form leaves an O(rho*u0)
    // tangential force per wall link, i.e. about rho*u0*nx per wall. The
    // Wen-invariant form cancels the advected-mass momentum and reports zero.
    let nx = 32;
    let ny = 34; // H = ny - 2 = 32 fluid rows.
    let mut sim: Simulation<f64> = SimConfig {
        nx,
        ny,
        nu: tau_to_nu(TAU),
        collision: TRT,
        edges: Edges {
            left: EdgeBC::Periodic,
            right: EdgeBC::Periodic,
            bottom: EdgeBC::MovingWall { u: [U0, 0.0] },
            top: EdgeBC::MovingWall { u: [U0, 0.0] },
        },
        ..Default::default()
    }
    .build()
    .unwrap();

    sim.init_with(|_, _| (RHO0, U0, 0.0));
    sim.set_force_probe(move |_, y| y == 0 || y == ny - 1);

    sim.step();
    let f_probe = sim.probed_force();
    let ladd_per_wall = RHO0 * U0 * nx as f64;
    let wen = 0.0;
    let band = 1.0e-10 * ladd_per_wall;

    println!(
        "ACC PROBE GALILEAN WEN-2014 row16: F_probe=({:.12e},{:.12e}), Ladd conventional prediction per wall rho*u0*nx={:.12e}, Wen invariant prediction={:.12e}, band={:.12e}",
        f_probe[0], f_probe[1], ladd_per_wall, wen, band
    );

    assert!(
        f_probe[0].abs() <= band,
        "ACC PROBE GALILEAN WEN-2014 row16: measured F_probe_x={:.12e}, band={band:.12e}, denominator=rho*u0*nx per wall={ladd_per_wall:.12e}; if Ladd-conventional, F ~= rho*u0*nx per wall, while Wen-invariant predicts 0",
        f_probe[0]
    );
}
