use lbm_core::particles::{Particle, ParticleSet, Sample};

fn fluid(u: [f64; 3]) -> impl Fn([f64; 3]) -> Sample {
    move |_| Sample { u, solid: false }
}

#[test]
fn straight_line_floor_crossing_is_interpolated_exactly() {
    let vel = [0.5, -0.25, -2.0];
    let mut set = ParticleSet::new(
        vec![Particle {
            pos: [1.0, 2.0, 1.0],
            vel,
            d: 0.1,
            rho_p: 1.0,
            exposure: 0.0,
        }],
        1.0,
        0.1,
        [0.0; 3],
    );
    let mut deposits = Vec::new();

    set.step_depositing(fluid(vel), None::<fn([f64; 3]) -> f64>, 0.0, &mut deposits);

    assert!(set.particles.is_empty());
    assert_eq!(deposits.len(), 1);
    assert_eq!(deposits[0].pos, [1.25, 1.875, 0.0]);
    assert_eq!(deposits[0].particle.pos, deposits[0].pos);
    assert_eq!(deposits[0].particle.vel, vel);
}

#[test]
fn deposited_and_suspended_counts_are_conserved_in_particle_order() {
    let vel = [0.0, 0.0, -1.0];
    let initial = vec![
        Particle {
            pos: [0.0, 0.0, 0.25],
            vel,
            d: 0.1,
            rho_p: 1.0,
            exposure: 10.0,
        },
        Particle {
            pos: [1.0, 0.0, 2.0],
            vel,
            d: 0.1,
            rho_p: 1.0,
            exposure: 20.0,
        },
        Particle {
            pos: [2.0, 0.0, 0.75],
            vel,
            d: 0.1,
            rho_p: 1.0,
            exposure: 30.0,
        },
    ];
    let n0 = initial.len();
    let mut set = ParticleSet::new(initial, 1.0, 0.1, [0.0; 3]);
    let mut deposits = Vec::new();

    set.step_depositing(fluid(vel), None::<fn([f64; 3]) -> f64>, 0.0, &mut deposits);

    assert_eq!(deposits.len() + set.particles.len(), n0);
    assert_eq!(deposits.len(), 2);
    assert_eq!(set.particles.len(), 1);
    assert_eq!(deposits[0].particle.exposure, 10.0);
    assert_eq!(deposits[1].particle.exposure, 30.0);
    assert_eq!(set.particles[0].exposure, 20.0);
    assert_eq!(deposits[0].pos, [0.0, 0.0, 0.0]);
    assert_eq!(deposits[1].pos, [2.0, 0.0, 0.0]);
}

#[test]
fn low_re_terminal_velocity_matches_analytic_stokes_settling() {
    let rho_f = 1.0;
    let rho_p = 1.01;
    let d = 0.1;
    let nu = 0.1;
    let g_mag = 1.0e-3;
    let mut set = ParticleSet::new(
        vec![Particle {
            pos: [0.0, 0.0, 10.0],
            vel: [0.0; 3],
            d,
            rho_p,
            exposure: 0.0,
        }],
        rho_f,
        nu,
        [0.0, 0.0, -g_mag],
    );
    let mut deposits = Vec::new();

    for _ in 0..2_000 {
        set.step_depositing(
            fluid([0.0; 3]),
            None::<fn([f64; 3]) -> f64>,
            -1.0,
            &mut deposits,
        );
    }

    assert!(deposits.is_empty());
    let got = -set.particles[0].vel[2];
    let want = (rho_p / rho_f - 1.0) * g_mag * d * d / (18.0 * nu);
    let re_p = got * d / nu;
    assert!(re_p < 0.1, "Re_p={re_p:e}");
    assert!(
        (got - want).abs() / want < 1.0e-3,
        "terminal velocity got {got:e}, want {want:e}"
    );
}
