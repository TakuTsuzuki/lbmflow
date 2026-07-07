//! Near-wall particle sampling behavior anchors.
//!
//! These tests pin conventions around `sample_grid` and particle-wall contact.
//! They are regression anchors for documented behavior, not claims that the
//! clamped sampler is a wall-resolved particle model.

use lbm_core::particles::{sample_grid, Particle, ParticleSet, Sample};

const RHO_F: f64 = 1.0;
const NU: f64 = 0.1;

fn assert_close(got: f64, want: f64, tol: f64, context: &str) {
    let err = (got - want).abs();
    assert!(
        err <= tol,
        "{context}: got={got:e}, want={want:e}, abs_err={err:e}, tol={tol:e}"
    );
}

fn assert_vec3_close(got: [f64; 3], want: [f64; 3], tol: f64, context: &str) {
    for axis in 0..3 {
        assert_close(
            got[axis],
            want[axis],
            tol,
            &format!("{context} axis {axis}"),
        );
    }
}

fn particle(pos: [f64; 3], vel: [f64; 3], d: f64, rho_p: f64) -> Particle {
    Particle {
        pos,
        vel,
        d,
        rho_p,
        exposure: 0.0,
    }
}

#[test]
fn out_of_domain_sample_equals_clamped_boundary_sample() {
    let dims = [4, 5, 3];
    let accessor = |x: usize, y: usize, z: usize| {
        let xf = x as f64;
        let yf = y as f64;
        let zf = z as f64;
        (
            [
                0.125 + 0.5 * xf - 0.25 * yf + 0.0625 * zf,
                -0.75 + 0.125 * xf + 0.5 * yf - 0.25 * zf,
                1.5 - 0.375 * xf + 0.25 * yf + 0.125 * zf,
            ],
            x == 0 && y == 2 && z == 2,
        )
    };

    let outside = [-2.0, 2.4, 99.0];
    let clamped = [0.0, 2.4, 2.0];
    let outside_sample = sample_grid(outside, dims, accessor);
    let clamped_sample = sample_grid(clamped, dims, accessor);

    assert_vec3_close(
        outside_sample.u,
        clamped_sample.u,
        0.0,
        "out-of-domain sample must be the clamped-position sample",
    );
    assert_eq!(
        outside_sample.solid,
        accessor(0, 2, 2).1,
        "solid flag must come from the clamped lower node"
    );
    assert_eq!(
        outside_sample.solid, clamped_sample.solid,
        "outside and explicit clamped samples must report the same solid flag"
    );
}

#[test]
fn wall_parallel_shear_does_not_induce_wall_normal_particle_drift() {
    let dims = [16, 8, 1];
    let shear = 2.0e-3;
    let y_wall = 0.5;
    let initial_y = 1.25;
    let steps = 4096;

    let sample = |pos: [f64; 3]| {
        sample_grid(pos, dims, |_, y, _| {
            if y == 0 {
                ([0.0; 3], true)
            } else {
                ([shear * (y as f64 - y_wall), 0.0, 0.0], false)
            }
        })
    };

    let mut set = ParticleSet::new(
        vec![particle([1.25, initial_y, 0.0], [0.0; 3], 0.2, RHO_F)],
        RHO_F,
        NU,
        [0.0; 3],
    );

    let mut min_y = initial_y;
    let mut max_y = initial_y;
    for step in 0..steps {
        set.step(sample, None::<fn([f64; 3]) -> f64>).unwrap();
        let p = &set.particles[0];
        min_y = min_y.min(p.pos[1]);
        max_y = max_y.max(p.pos[1]);
        assert!(
            p.pos[1] >= 1.0,
            "step {step}: particle entered the solid lower-node contact band: pos={:?}",
            p.pos
        );
    }

    // Derivation of the wall-normal floor:
    // Couette input is affine, u=(gamma*(y-y_wall), 0, 0). Trilinear
    // interpolation reproduces affine fields exactly, and every sampled node
    // has u_y=0. With rho_p=rho_f and g_y=0, the particle update maps v_y=0 to
    // v_y=0, then y_{n+1}=y_n+0. The exact arithmetic path adds signed zero to
    // y, so the derived drift floor for this case is 0.
    let drift = set.particles[0].pos[1] - initial_y;
    assert_eq!(
        drift.to_bits(),
        0.0f64.to_bits(),
        "wall-parallel shear produced wall-normal drift={drift:e}; min_y={min_y:e}, max_y={max_y:e}"
    );
    assert_eq!(min_y.to_bits(), initial_y.to_bits());
    assert_eq!(max_y.to_bits(), initial_y.to_bits());
}

#[test]
fn settling_particle_stops_at_solid_wall_without_entering_solid() {
    let wall = |pos: [f64; 3]| Sample {
        u: [0.0; 3],
        solid: pos[2] < 0.0,
    };
    let mut set = ParticleSet::new(
        vec![particle([0.0, 0.0, 0.025], [0.0, 0.0, -0.03], 0.4, 2.0)],
        RHO_F,
        NU,
        [0.0, 0.0, -2.0e-2],
    );

    let mut previous_speed = set.particles[0].vel[2].abs();
    let mut contact_pos = None;
    for step in 0..256 {
        set.step(wall, None::<fn([f64; 3]) -> f64>).unwrap();
        let p = &set.particles[0];
        assert!(
            p.pos[2] >= 0.0,
            "step {step}: particle passed into solid cells: pos={:?}",
            p.pos
        );
        let speed = p.vel[2].abs();
        assert!(
            speed <= previous_speed + 1.0e-15,
            "step {step}: wall-normal speed increased in drag-dominated approach: previous={previous_speed:e}, current={speed:e}"
        );
        previous_speed = speed;

        if p.vel[2].to_bits() == 0.0f64.to_bits() {
            contact_pos = Some(p.pos);
            break;
        }
    }

    let contact_pos = contact_pos.expect("particle never contacted the wall");
    assert_eq!(
        set.particles.len(),
        1,
        "staircase-wall contact is a stop/reflect convention, not deposition"
    );

    // Convention pin: with restitution 0, the blocked wall-normal component is
    // set to zero and the particle remains suspended at the last fluid point.
    for step in 0..16 {
        set.step(wall, None::<fn([f64; 3]) -> f64>).unwrap();
        let p = &set.particles[0];
        assert_vec3_close(
            p.pos,
            contact_pos,
            0.0,
            &format!("post-contact position at step {step}"),
        );
        assert_eq!(
            p.vel[2].to_bits(),
            0.0f64.to_bits(),
            "post-contact wall-normal velocity should remain stopped"
        );
    }
}

#[test]
fn restitution_extremes_pin_stop_and_elastic_reflection_conventions() {
    let wall = |pos: [f64; 3]| Sample {
        u: [0.0, 0.0, -0.5],
        solid: pos[2] < 0.0,
    };
    let initial = particle([0.0, 0.0, 0.2], [0.0, 0.0, -0.5], 0.4, RHO_F);

    let mut stopped =
        ParticleSet::new(vec![initial.clone()], RHO_F, NU, [0.0; 3]).with_restitution(0.0);
    stopped.step(wall, None::<fn([f64; 3]) -> f64>).unwrap();
    assert_eq!(
        stopped.particles.len(),
        1,
        "zero-restitution wall contact must not deposit/remove the particle"
    );
    assert_vec3_close(
        stopped.particles[0].pos,
        initial.pos,
        0.0,
        "zero-restitution contact position",
    );
    assert_eq!(
        stopped.particles[0].vel[2].to_bits(),
        0.0f64.to_bits(),
        "zero restitution must stop the blocked wall-normal component"
    );

    let mut elastic =
        ParticleSet::new(vec![initial.clone()], RHO_F, NU, [0.0; 3]).with_restitution(1.0);
    elastic.step(wall, None::<fn([f64; 3]) -> f64>).unwrap();
    assert_vec3_close(
        elastic.particles[0].pos,
        initial.pos,
        0.0,
        "unit-restitution contact position",
    );
    assert_close(
        elastic.particles[0].vel[2].abs(),
        initial.vel[2].abs(),
        0.0,
        "unit restitution must preserve wall-normal speed magnitude",
    );
    assert!(
        elastic.particles[0].vel[2] > 0.0,
        "unit restitution must reverse the wall-normal sign"
    );
}

#[test]
fn clamped_near_wall_sample_bias_is_frozen_characterization() {
    let dims = [8, 8, 1];
    let shear = 1.0e-2;
    let sample = |pos: [f64; 3]| {
        sample_grid(pos, dims, |_, y, _| {
            let solid = y == 0;
            ([shear * y as f64, 0.0, 0.0], solid)
        })
    };

    let near_wall = [3.25, -0.25, 0.0];
    let interior = [3.25, 2.25, 0.0];
    let near = sample(near_wall);
    let inside = sample(interior);
    let near_analytic = shear * near_wall[1];
    let interior_analytic = shear * interior[1];
    let near_error = near.u[0] - near_analytic;
    let interior_error = inside.u[0] - interior_analytic;

    // This is a convention characterization, not an accuracy claim: the
    // out-of-domain y=-0.25 sample is clamped to the resting boundary node
    // y=0, so the measured difference from the analytic linear continuation is
    // gamma*(0 - -0.25) = 2.5e-3. The interior affine sample remains exact.
    assert_close(near_error, 2.5e-3, 1.0e-15, "clamped near-wall shear bias");
    assert_close(
        interior_error,
        0.0,
        1.0e-15,
        "interior linear-shear interpolation error",
    );
    assert!(
        near.solid,
        "clamped near-wall sample must report the clamped lower-node solid flag"
    );
    assert!(!inside.solid, "interior sample should remain fluid");
}
