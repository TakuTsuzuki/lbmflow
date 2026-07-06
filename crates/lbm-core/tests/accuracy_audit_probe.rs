//! ACC PROBE: adversarial accuracy audits for momentum-exchange force probes
//! and wall-coupling. References are derived from continuum or discrete
//! momentum balances in-place; engine internals are not used as an oracle.

mod common;

use common::metrics::*;
use common::run_to_steady;
use lbm_core::compat::prelude::*;
use std::f64::consts::PI;

const TRT: Collision = Collision::Trt {
    magic: Collision::MAGIC_STD,
};

fn tau_to_nu(tau: f64) -> f64 {
    (tau - 0.5) / 3.0
}

fn bounce_box(nx: usize, ny: usize, nu: f64, force: [f64; 2]) -> Simulation<f64> {
    SimConfig {
        nx,
        ny,
        nu,
        collision: TRT,
        edges: Edges {
            left: EdgeBC::BounceBack,
            right: EdgeBC::BounceBack,
            bottom: EdgeBC::BounceBack,
            top: EdgeBC::BounceBack,
        },
        force,
        ..Default::default()
    }
    .build()
    .unwrap()
}

fn poiseuille_sim(nx: usize, ny: usize, nu: f64, g: f64) -> Simulation<f64> {
    SimConfig {
        nx,
        ny,
        nu,
        collision: TRT,
        edges: Edges {
            left: EdgeBC::Periodic,
            right: EdgeBC::Periodic,
            bottom: EdgeBC::BounceBack,
            top: EdgeBC::BounceBack,
        },
        force: [g, 0.0],
        ..Default::default()
    }
    .build()
    .unwrap()
}

fn couette_sim(nx: usize, ny: usize, nu: f64, rho0: f64, u_wall: f64) -> Simulation<f64> {
    let mut sim: Simulation<f64> = SimConfig {
        nx,
        ny,
        nu,
        collision: TRT,
        edges: Edges {
            left: EdgeBC::Periodic,
            right: EdgeBC::Periodic,
            bottom: EdgeBC::BounceBack,
            top: EdgeBC::MovingWall { u: [u_wall, 0.0] },
        },
        ..Default::default()
    }
    .build()
    .unwrap();
    sim.init_with(|_, _| (rho0, 0.0, 0.0));
    sim
}

fn steady_or_panic(sim: &mut Simulation<f64>, label: &str) {
    assert!(
        run_to_steady(sim, 500, 1.0e-11, 200_000),
        "{label}: steady=false, time={}",
        sim.time()
    );
}

fn couette_profile(sim: &Simulation<f64>) -> Vec<f64> {
    (1..sim.ny() - 1).map(|y| sim.ux(0, y)).collect()
}

fn couette_reference(ny: usize, u_wall: f64) -> Vec<f64> {
    let h = (ny - 2) as f64;
    (1..ny - 1).map(|y| u_wall * (y as f64 - 0.5) / h).collect()
}

fn obstacle_l(x: usize, y: usize) -> bool {
    ((24..=29).contains(&x) && (12..=20).contains(&y))
        || ((24..=34).contains(&x) && (12..=14).contains(&y))
}

fn mirrored_obstacle_l(nx: usize, x: usize, y: usize) -> bool {
    let ox = nx - 1 - x;
    obstacle_l(ox, y)
}

#[test]
fn a1_flagship_per_step_global_momentum_ledger_closure() {
    // Discrete derivation: in a closed box, bounce-back exchanges momentum
    // pairwise across every fluid-solid link. Newton's third law makes the
    // impulse on all probed solid cells equal and opposite to the impulse lost
    // by the fluid at those links. The only remaining external impulse on the
    // fluid is the uniform Guo body force applied to every fluid cell. Thus
    // for one step, p_fluid(t+1)-p_fluid(t)=N_fluid*F-F_wall(t..t+1).
    //
    // Caveat: total_momentum() reports physical momentum with the Guo F/2
    // velocity shift. Since F is constant, that half-step offset is identical
    // at t and t+1 and cancels in the momentum difference; we sample after
    // step 5 so start-up transients cannot be mistaken for a timing issue.
    let nx = 40;
    let ny = 32;
    let force = [4.0e-6, 1.0e-6];
    let mut sim = bounce_box(nx, ny, tau_to_nu(0.8), force);
    sim.set_solid_region(|x, y| (12..=18).contains(&x) && (10..=16).contains(&y));
    sim.set_force_probe(|_, _| true);
    sim.init_with(|x, y| {
        let ux = 0.02 * (2.0 * PI * y as f64 / ny as f64).sin();
        let uy = 0.01 * (2.0 * PI * x as f64 / nx as f64).sin();
        (1.0, ux, uy)
    });
    sim.run(5);

    for sample in 0..10 {
        let p0 = sim.total_momentum();
        sim.step();
        let p1 = sim.total_momentum();
        let fp = sim.probed_force();
        let n_fluid = sim.fluid_cell_count() as f64;
        let expected = [n_fluid * force[0] - fp[0], n_fluid * force[1] - fp[1]];
        let measured = [p1[0] - p0[0], p1[1] - p0[1]];
        // Band model: the identity is exact in infinite precision; the floor
        // is the cancellation error of differencing two O(N_cells)-term sums
        // of magnitude ~|p|, plus the probe/force sums. Scale the band with
        // the sum of the ledger-term magnitudes (NOT only the net force,
        // which is orders smaller than |p| and under-floors the round-off).
        let den = [
            p0[0].abs() + n_fluid * force[0].abs() + fp[0].abs(),
            p0[1].abs() + n_fluid * force[1].abs() + fp[1].abs(),
        ];
        println!(
            "ACC PROBE A1: sample={} dp=({:.12e},{:.12e}) expected=({:.12e},{:.12e}) F_probe=({:.12e},{:.12e}) den=({:.12e},{:.12e})",
            sample, measured[0], measured[1], expected[0], expected[1], fp[0], fp[1], den[0], den[1]
        );
        for c in 0..2 {
            let err = (measured[c] - expected[c]).abs();
            let band = 1.0e-11 * den[c];
            assert!(
                err <= band,
                "ACC PROBE A1 component {c}: measured dp={:.12e}, expected={:.12e}, abs_err={:.12e}, band={:.12e}, denominator=|p_t|+N_fluid*|F|+|F_probe|={:.12e}; O(1e-2 relative) suggests probe-completeness/timing finding",
                measured[c],
                expected[c],
                err,
                band,
                den[c]
            );
        }
    }
}

#[test]
fn a2_steady_poiseuille_wall_friction_balance() {
    // Steady derivation: for periodic x and no acceleration, the summed
    // x-momentum equation over all fluid cells has zero left-hand side. The
    // uniform body force injects +g*N_fluid into the fluid each step, so the
    // walls must absorb exactly that: the force ON the probed walls
    // (probed_force sign convention: force exerted BY the fluid ON the solid,
    // cf. A1 ledger and A5 drag sign) sums to +g*N_fluid. The geometry and
    // forcing are symmetric about the centerline, so each wall carries half.
    // Normal (y) components are NOT zero: a resting wall in near-equilibrium
    // fluid receives the static-pressure push. Per column, the links entering
    // the top wall (c_y > 0: q2 w=1/9, q5 and q6 w=1/36) transfer
    // 2*rho*(1/9 + 1/36 + 1/36) = rho/3 = rho*cs^2; the O(u^2) equilibrium
    // corrections cancel exactly in this sum (+4.5(c.u)^2 terms of q5/q6
    // offset the -1.5u^2 terms), so F_top_y = +nx*rho*cs^2 and
    // F_bottom_y = -nx*rho*cs^2 to round-off in steady state.
    let nx = 8;
    let ny = 34;
    let g = 1.0e-6;
    let mut sim = poiseuille_sim(nx, ny, tau_to_nu(0.8), g);
    steady_or_panic(&mut sim, "ACC PROBE A2");

    // The public API has one active probe, so retarget it between two reads
    // and advance one steady step after each retarget to measure that set.
    sim.set_force_probe(|_, y| y == 0);
    sim.step();
    let bottom = sim.probed_force();
    sim.set_force_probe(move |_, y| y == ny - 1);
    sim.step();
    let top = sim.probed_force();

    let n_fluid = sim.fluid_cell_count() as f64;
    let den = g * n_fluid;
    let sum_x = top[0] + bottom[0] - den;
    let sym_x = top[0] - bottom[0];
    let p_push = nx as f64 / 3.0;
    println!(
        "ACC PROBE A2: F_top=({:.12e},{:.12e}) F_bottom=({:.12e},{:.12e}) gN={:.12e} sum_x_res={:.12e} sym_x={:.12e}",
        top[0], top[1], bottom[0], bottom[1], den, sum_x, sym_x
    );
    assert!(
        sum_x.abs() <= 1.0e-10 * den,
        "ACC PROBE A2 wall-friction balance: measured residual={sum_x:.12e}, band={:.12e}, denominator=g*N_fluid={den:.12e}, top_x={:.12e}, bottom_x={:.12e}",
        1.0e-10 * den,
        top[0],
        bottom[0]
    );
    assert!(
        sym_x.abs() <= 1.0e-12 * den,
        "ACC PROBE A2 top/bottom symmetry: measured F_top_x-F_bottom_x={sym_x:.12e}, band={:.12e}, denominator=g*N_fluid={den:.12e}",
        1.0e-12 * den
    );
    for (label, value, expected) in [("top_y", top[1], p_push), ("bottom_y", bottom[1], -p_push)] {
        let err = (value - expected).abs();
        assert!(
            err <= 1.0e-9 * p_push,
            "ACC PROBE A2 static-pressure normal force {label}: measured={value:.12e}, expected={expected:.12e}, abs_err={err:.12e}, band={:.12e}, denominator=nx*rho*cs^2={p_push:.12e}",
            1.0e-9 * p_push
        );
    }
    let y_cancel = (top[1] + bottom[1]).abs();
    assert!(
        y_cancel <= 1.0e-12 * p_push,
        "ACC PROBE A2 normal-force cancellation: measured |F_top_y+F_bottom_y|={y_cancel:.12e}, band={:.12e}, denominator=nx*rho*cs^2={p_push:.12e}",
        1.0e-12 * p_push
    );
}

#[test]
fn a3_moving_wall_momentum_term_uses_local_density() {
    // Moving-wall derivation: Ladd bounce-back adds
    // 2*w_q*rho*(c_q.u_w)/c_s^2 to the reflected population, where rho is the
    // adjacent local fluid density. If rho were hardcoded to 1, a uniform
    // rho0=1.05 run would inject only 1/rho0 of the required wall impulse,
    // making the Couette slope and effective wall speed smaller by that
    // factor. The incompressible Couette solution itself is independent of
    // uniform rho0: ux(y_w)=U*y_w/H.
    let nx = 8;
    let ny = 34;
    let u_wall = 0.1;
    let mut sim = couette_sim(nx, ny, tau_to_nu(0.8), 1.05, u_wall);
    steady_or_panic(&mut sim, "ACC PROBE A3");
    let actual = couette_profile(&sim);
    let reference = couette_reference(ny, u_wall);
    let err = linf_rel(&actual, &reference, u_wall);
    let scale = actual[actual.len() - 1] / reference[reference.len() - 1];
    println!(
        "ACC PROBE A3: linf_rel={err:.12e} scale_top_cell={scale:.12e} rho0=1.05 time={}",
        sim.time()
    );
    assert!(
        err <= 1.0e-8,
        "ACC PROBE A3 density-dependent moving-wall term: measured linf_rel={err:.12e}, band=1.000000000000e-08, normalization=floor U={u_wall:.12e}, scale_top_cell={scale:.12e} (rho-hardcoding would be ~0.952)"
    );
}

#[test]
fn a4_probe_mirror_equivariance() {
    // Symmetry derivation: D2Q9 velocities and weights are closed under
    // x-reflection. Mirroring geometry, initial data, and boundary driving
    // maps each population path to a reflected path, so the force on the
    // mirrored obstacle must transform as [Fx,Fy] -> [-Fx,Fy].
    //
    // A left VelocityInlet with u=[+0.05,0] mirrors to a right VelocityInlet
    // with u=[-0.05,0]; the pressure outlet swaps to the opposite side.
    let nx = 64;
    let ny = 34;
    let mut orig: Simulation<f64> = SimConfig {
        nx,
        ny,
        nu: tau_to_nu(0.8),
        collision: TRT,
        edges: Edges {
            left: EdgeBC::VelocityInlet { u: [0.05, 0.0] },
            right: EdgeBC::PressureOutlet { rho: 1.0 },
            bottom: EdgeBC::BounceBack,
            top: EdgeBC::BounceBack,
        },
        ..Default::default()
    }
    .build()
    .unwrap();
    orig.set_solid_region(obstacle_l);
    orig.set_force_probe(obstacle_l);
    orig.run(400);
    let fo = orig.probed_force();

    let mut mirror: Simulation<f64> = SimConfig {
        nx,
        ny,
        nu: tau_to_nu(0.8),
        collision: TRT,
        edges: Edges {
            left: EdgeBC::PressureOutlet { rho: 1.0 },
            right: EdgeBC::VelocityInlet { u: [-0.05, 0.0] },
            bottom: EdgeBC::BounceBack,
            top: EdgeBC::BounceBack,
        },
        ..Default::default()
    }
    .build()
    .unwrap();
    mirror.set_solid_region(|x, y| mirrored_obstacle_l(nx, x, y));
    mirror.set_force_probe(|x, y| mirrored_obstacle_l(nx, x, y));
    mirror.run(400);
    let fm = mirror.probed_force();

    let fx_err = (fm[0] + fo[0]).abs();
    let fx_band = 1.0e-12 * fo[0].abs();
    let fy_err = (fm[1] - fo[1]).abs();
    let fy_den = fo[1].abs().max(1.0e-30);
    let fy_band = (1.0e-12 * fy_den).max(1.0e-14);
    println!(
        "ACC PROBE A4: F_orig=({:.12e},{:.12e}) F_mirror=({:.12e},{:.12e}) fx_err={:.12e} fy_err={:.12e} fy_floor=1e-30 fy_abs_floor=1e-14",
        fo[0], fo[1], fm[0], fm[1], fx_err, fy_err
    );
    assert!(
        fx_err <= fx_band,
        "ACC PROBE A4 mirror Fx: measured |Fx_mirror+Fx_orig|={fx_err:.12e}, band={fx_band:.12e}, denominator=|Fx_orig|={:.12e}",
        fo[0].abs()
    );
    assert!(
        fy_err <= fy_band,
        "ACC PROBE A4 mirror Fy: measured |Fy_mirror-Fy_orig|={fy_err:.12e}, band={fy_band:.12e}, denominator=max(|Fy_orig|,1e-30)={fy_den:.12e} with abs floor 1e-14"
    );
}

#[test]
fn a5_moving_wall_probed_shear_matches_exact_couette_friction() {
    // Couette derivation: at steady state ux=U*y_w/H, hence du/dy=U/H.
    // Newtonian shear stress on a wall is tau_xy=mu*du/dy=rho*nu*U/H per
    // unit wall length. probed_force returns force on the probed solid. The
    // moving top wall is dragged backward by the fluid, so F_top_x<0 for
    // U>0; the static bottom wall feels the equal opposite force.
    let nx = 8;
    let ny = 34;
    let u_wall = 0.1;
    for tau in [0.6, 1.0] {
        let nu = tau_to_nu(tau);
        let mut sim = couette_sim(nx, ny, nu, 1.0, u_wall);
        steady_or_panic(&mut sim, "ACC PROBE A5");

        sim.set_force_probe(move |_, y| y == ny - 1);
        sim.step();
        let top = sim.probed_force();
        sim.set_force_probe(|_, y| y == 0);
        sim.step();
        let bottom = sim.probed_force();

        let h = (ny - 2) as f64;
        let expected_mag = nu * u_wall / h * nx as f64;
        let mag_err = (top[0].abs() - expected_mag).abs();
        let mag_band = 1.0e-8 * expected_mag;
        // The public API exposes one probe set at a time, so the two walls
        // are read on CONSECUTIVE steps; near the 1e-11 steady criterion the
        // field still drifts ~O(1e-10) relative per step, which bounds the
        // achievable cancellation between the two reads (measured 3.1e-10
        // relative at tau=0.6). Band 1e-8 keeps ~30x headroom while staying
        // far below any physical probe asymmetry (O(u^2/H) ~ 1e-4 relative).
        let cancel = top[0] + bottom[0];
        let cancel_band = 1.0e-8 * top[0].abs();
        println!(
            "ACC PROBE A5: tau={tau:.3} F_top=({:.12e},{:.12e}) F_bottom=({:.12e},{:.12e}) expected_mag={:.12e} sign_convention=force_on_solid_top_negative_bottom_positive",
            top[0], top[1], bottom[0], bottom[1], expected_mag
        );
        assert!(
            mag_err <= mag_band,
            "ACC PROBE A5 Couette top shear tau={tau:.3}: measured |F_top_x|={:.12e}, expected={expected_mag:.12e}, abs_err={mag_err:.12e}, band={mag_band:.12e}, denominator=rho*nu*U/H*nx={expected_mag:.12e}",
            top[0].abs()
        );
        assert!(
            top[0] < 0.0,
            "ACC PROBE A5 sign tau={tau:.3}: measured F_top_x={:.12e}, expected negative because probed_force is force on moving solid opposing +x wall motion",
            top[0]
        );
        assert!(
            cancel.abs() <= cancel_band,
            "ACC PROBE A5 wall-force cancellation tau={tau:.3}: measured F_top_x+F_bottom_x={cancel:.12e}, band={cancel_band:.12e}, denominator=|F_top_x|={:.12e}",
            top[0].abs()
        );
    }
}
