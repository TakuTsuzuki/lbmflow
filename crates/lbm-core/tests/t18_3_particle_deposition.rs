//! T18.3 CR-3 deposition acceptance tests.
//!
//! This file is intentionally written against the frozen future API from
//! `docs/DISPERSED_DEPOSITION.md`:
//! `DepositEvent { pos, particle }` and
//! `ParticleSet::step_depositing(sample, exposure_rate, floor_z, &mut events)`.
//! It is expected not to compile until CR-3 lands.

use approx::abs_diff_eq;
use lbm_core::particles::{sample_grid, DepositEvent, Particle, ParticleSet, Sample};
use lbm_core::prelude::*;

fn moving_particle(pos: [f64; 3], vel: [f64; 3], rho_p: f64) -> Particle {
    Particle {
        pos,
        vel,
        d: 1.0,
        rho_p,
        exposure: 0.0,
    }
}

fn assert_vec3_exact_bits(got: [f64; 3], want: [f64; 3], context: &str) {
    for axis in 0..3 {
        assert_eq!(
            got[axis].to_bits(),
            want[axis].to_bits(),
            "{context}: axis {axis} bit mismatch, got={:e}, want={:e}",
            got[axis],
            want[axis]
        );
    }
}

fn assert_particle_exact_bits(got: &Particle, want: &Particle, context: &str) {
    assert_vec3_exact_bits(got.pos, want.pos, &format!("{context} pos"));
    assert_vec3_exact_bits(got.vel, want.vel, &format!("{context} vel"));
    assert_eq!(
        got.d.to_bits(),
        want.d.to_bits(),
        "{context}: diameter bit mismatch, got={:e}, want={:e}",
        got.d,
        want.d
    );
    assert_eq!(
        got.rho_p.to_bits(),
        want.rho_p.to_bits(),
        "{context}: rho_p bit mismatch, got={:e}, want={:e}",
        got.rho_p,
        want.rho_p
    );
    assert_eq!(
        got.exposure.to_bits(),
        want.exposure.to_bits(),
        "{context}: exposure bit mismatch, got={:e}, want={:e}",
        got.exposure,
        want.exposure
    );
}

fn assert_deposit_records_bit_match(a: &[DepositEvent], b: &[DepositEvent], context: &str) {
    assert_eq!(
        a.len(),
        b.len(),
        "{context}: deposited record count mismatch, left={}, right={}",
        a.len(),
        b.len()
    );
    for (i, (ea, eb)) in a.iter().zip(b).enumerate() {
        assert_vec3_exact_bits(ea.pos, eb.pos, &format!("{context} event[{i}] crossing"));
        assert_particle_exact_bits(
            &ea.particle,
            &eb.particle,
            &format!("{context} event[{i}] particle"),
        );
    }
}

#[test]
fn t18_3_floor_crossing_is_interpolated_and_conserves_particle_count() {
    let floor_z = 0.0;
    let initial = vec![
        moving_particle([1.0, 2.0, 0.25], [0.25, -0.5, -0.5], 1.0),
        moving_particle([3.0, 4.0, 1.25], [-0.25, 0.125, -0.25], 1.0),
    ];
    let mut set = ParticleSet::new(initial.clone(), 1.0, 0.1, [0.0; 3]);
    let mut deposits = Vec::<DepositEvent>::new();
    let sample = |p: [f64; 3]| Sample {
        // rho_p == rho_f and g == 0, so u == v keeps the step segment straight.
        u: if p[0] < 2.0 {
            initial[0].vel
        } else {
            initial[1].vel
        },
        solid: false,
    };

    set.step_depositing(sample, None::<fn([f64; 3]) -> f64>, floor_z, &mut deposits)
        .unwrap();

    assert_eq!(
        deposits.len(),
        1,
        "first step deposited count = {}, want 1",
        deposits.len()
    );
    assert_eq!(
        set.particles.len(),
        1,
        "first step suspended count = {}, want 1",
        set.particles.len()
    );

    let crossing_t = (floor_z - initial[0].pos[2]) / initial[0].vel[2];
    let expected_crossing = [
        initial[0].pos[0] + crossing_t * initial[0].vel[0],
        initial[0].pos[1] + crossing_t * initial[0].vel[1],
        floor_z,
    ];
    for axis in 0..3 {
        assert!(
            abs_diff_eq!(
                deposits[0].pos[axis],
                expected_crossing[axis],
                epsilon = 0.0
            ),
            "axis {axis}: crossing got={:e}, want={:e}",
            deposits[0].pos[axis],
            expected_crossing[axis]
        );
    }

    for _ in 0..8 {
        let before_total = set.particles.len() + deposits.len();
        set.step_depositing(sample, None::<fn([f64; 3]) -> f64>, floor_z, &mut deposits)
            .unwrap();
        let after_total = set.particles.len() + deposits.len();
        assert_eq!(
            after_total, before_total,
            "suspended + deposited count changed from {before_total} to {after_total}"
        );
    }
    assert_eq!(
        deposits.len() + set.particles.len(),
        initial.len(),
        "global particle count not conserved: deposited={}, suspended={}, initial={}",
        deposits.len(),
        set.particles.len(),
        initial.len()
    );
    assert_eq!(
        set.particles.len(),
        0,
        "all particles should have crossed after the multi-step run; suspended={}",
        set.particles.len()
    );
}

#[test]
fn t18_3_deposit_records_are_in_particle_index_order_within_a_step() {
    let particles = vec![
        moving_particle([10.0, 0.0, 0.75], [0.0, 0.0, -1.0], 1.0),
        moving_particle([20.0, 0.0, 0.25], [0.0, 0.0, -1.0], 1.0),
        moving_particle([30.0, 0.0, 2.00], [0.0, 0.0, -0.1], 1.0),
    ];
    let mut set = ParticleSet::new(particles.clone(), 1.0, 0.1, [0.0; 3]);
    let mut deposits = Vec::<DepositEvent>::new();
    let sample = |p: [f64; 3]| Sample {
        u: if p[0] < 15.0 {
            particles[0].vel
        } else if p[0] < 25.0 {
            particles[1].vel
        } else {
            particles[2].vel
        },
        solid: false,
    };

    set.step_depositing(sample, None::<fn([f64; 3]) -> f64>, 0.0, &mut deposits)
        .unwrap();

    assert_eq!(
        deposits.len(),
        2,
        "same-step deposited count = {}, want 2",
        deposits.len()
    );
    // Frozen payload convention (DISPERSED_DEPOSITION.md §5 CR-3): the recorded
    // particle carries the deposition state — pos = interpolated crossing, vel =
    // impact velocity. Here slip = 0 and g = 0, so vel is unchanged and only pos
    // moves to the crossing point.
    let expected0 = Particle {
        pos: [10.0, 0.0, 0.0],
        ..particles[0].clone()
    };
    let expected1 = Particle {
        pos: [20.0, 0.0, 0.0],
        ..particles[1].clone()
    };
    assert_vec3_exact_bits(deposits[0].pos, expected0.pos, "event[0] crossing");
    assert_vec3_exact_bits(deposits[1].pos, expected1.pos, "event[1] crossing");
    assert_particle_exact_bits(&deposits[0].particle, &expected0, "event[0]");
    assert_particle_exact_bits(&deposits[1].particle, &expected1, "event[1]");
}

#[test]
fn t18_3_deposition_map_is_partition_invariant_and_order_stable() {
    type S3<H> = Solver<D3Q19, f64, CpuScalar, H>;

    let n = 16usize;
    let spec = GlobalSpec::<f64> {
        dims: [n, n, n],
        nu: 0.02,
        periodic: [true, true, true],
        ..Default::default()
    };
    let k = 2.0 * std::f64::consts::PI / n as f64;
    let init = move |x: usize, y: usize, z: usize| {
        let (xx, yy, zz) = (k * x as f64, k * y as f64, k * z as f64);
        (
            1.0,
            [
                0.03 * xx.sin() * yy.cos() * zz.cos(),
                -0.03 * xx.cos() * yy.sin() * zz.cos(),
                -0.015,
            ],
        )
    };

    let mut mono: S3<LocalPeriodic> = Solver::new(
        &spec,
        &[],
        &[],
        [1, 1, 1],
        CpuScalar::default(),
        LocalPeriodic,
    );
    let mut split: S3<InProcess> =
        Solver::new(&spec, &[], &[], [2, 2, 2], CpuScalar::default(), InProcess);
    mono.init_with(init);
    split.init_with(init);
    for _ in 0..24 {
        mono.step();
        split.step();
    }

    let mono_fields = (mono.gather_ux(), mono.gather_uy(), mono.gather_uz());
    let split_fields = (split.gather_ux(), split.gather_uy(), split.gather_uz());
    assert_eq!(
        mono_fields.0, split_fields.0,
        "ux field is not bit-identical"
    );
    assert_eq!(
        mono_fields.1, split_fields.1,
        "uy field is not bit-identical"
    );
    assert_eq!(
        mono_fields.2, split_fields.2,
        "uz field is not bit-identical"
    );

    let seeds = vec![
        moving_particle([2.25, 2.25, 1.10], [0.0; 3], 1.4),
        moving_particle([7.50, 5.25, 1.60], [0.0; 3], 1.4),
        moving_particle([11.75, 9.50, 2.10], [0.0; 3], 1.4),
        moving_particle([13.25, 12.0, 2.60], [0.0; 3], 1.4),
    ];
    let mut mono_particles = ParticleSet::new(seeds.clone(), 1.0, 0.05, [0.0, 0.0, -0.05]);
    let mut split_particles = ParticleSet::new(seeds, 1.0, 0.05, [0.0, 0.0, -0.05]);
    let mut mono_deposits = Vec::<DepositEvent>::new();
    let mut split_deposits = Vec::<DepositEvent>::new();

    let sample_from = |fields: &(Vec<f64>, Vec<f64>, Vec<f64>), pos: [f64; 3]| -> Sample {
        sample_grid(pos, spec.dims, |x, y, z| {
            let i = (z * spec.dims[1] + y) * spec.dims[0] + x;
            ([fields.0[i], fields.1[i], fields.2[i]], false)
        })
    };

    for step in 0..80 {
        mono_particles
            .step_depositing(
                |pos| sample_from(&mono_fields, pos),
                None::<fn([f64; 3]) -> f64>,
                0.0,
                &mut mono_deposits,
            )
            .unwrap();
        split_particles
            .step_depositing(
                |pos| sample_from(&split_fields, pos),
                None::<fn([f64; 3]) -> f64>,
                0.0,
                &mut split_deposits,
            )
            .unwrap();
        assert_deposit_records_bit_match(
            &mono_deposits,
            &split_deposits,
            &format!("after particle step {step}"),
        );
        assert_eq!(
            mono_particles.particles, split_particles.particles,
            "suspended particles differ after particle step {step}"
        );
    }
}
