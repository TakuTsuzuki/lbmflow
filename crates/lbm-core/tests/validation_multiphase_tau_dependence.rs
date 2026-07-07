//! V&V lane 1.7 / P20: Shan-Chen coexistence tau-dependence sweep.
//!
//! In the continuum pseudopotential model, flat-interface coexistence
//! densities are set by the equation of state and should not depend on the
//! BGK/TRT relaxation time. Li et al. 2012 document that discrete forcing
//! schemes can shift the mechanical-stability condition, producing a residual
//! tau-dependent coexistence curve. This test characterizes that pitfall for
//! the current Guo-forced LBMFlow Shan-Chen path.

use lbm_core::compat::multiphase::ShanChen;
use lbm_core::compat::prelude::*;

const NX: usize = 64;
const NY: usize = 128;
const G: f64 = -5.0;
const STEPS: usize = 30_000;
const TAUS: [f64; 3] = [0.6, 1.0, 1.4];
const MONOTONE_ABS_TOL: f64 = 1.0e-5;

#[derive(Clone, Copy, Debug)]
struct CoexistencePoint {
    tau: f64,
    rho_l: f64,
    rho_v: f64,
}

fn nu_from_tau(tau: f64) -> f64 {
    (tau - 0.5) / 3.0
}

fn rel_shift(value: f64, reference: f64) -> f64 {
    ((value - reference) / reference).abs()
}

fn bulk_mean(sim: &Simulation<f64>, ranges: &[(usize, usize)]) -> f64 {
    let mut sum = 0.0;
    let mut count = 0usize;
    for &(y0, y1) in ranges {
        for y in y0..y1 {
            for x in 0..NX {
                sum += sim.rho(x, y);
                count += 1;
            }
        }
    }
    sum / count as f64
}

fn run_flat_interface(tau: f64) -> CoexistencePoint {
    let mut sim: Simulation<f64> = SimConfig {
        nx: NX,
        ny: NY,
        nu: nu_from_tau(tau),
        ..Default::default()
    }
    .build()
    .unwrap();
    sim.init_with(|_, y| {
        let liquid = y >= NY / 4 && y < 3 * NY / 4;
        (if liquid { 2.0 } else { 0.15 }, 0.0, 0.0)
    });

    let sc = ShanChen::new(G);
    for _ in 0..STEPS {
        sim.force_field_mut().fill([0.0; 2]);
        sc.update_force(&mut sim);
        sim.step();
    }

    assert!(
        sim.rho_field().iter().all(|rho| rho.is_finite()),
        "SC tau={tau} produced non-finite density"
    );
    CoexistencePoint {
        tau,
        // Initial interfaces sit at y=32 and y=96; these windows stay well
        // inside the equilibrated bulk phases.
        rho_l: bulk_mean(&sim, &[(48, 80)]),
        rho_v: bulk_mean(&sim, &[(8, 24), (104, 120)]),
    }
}

fn assert_nondecreasing(values: [f64; 3], label: &str) {
    assert!(
        values[0] <= values[1] + MONOTONE_ABS_TOL
            && values[1] <= values[2] + MONOTONE_ABS_TOL,
        "{label} is not nondecreasing over tau={TAUS:?} within abs tol {MONOTONE_ABS_TOL:e}: values={values:?}"
    );
}

fn assert_nonincreasing(values: [f64; 3], label: &str) {
    assert!(
        values[0] + MONOTONE_ABS_TOL >= values[1]
            && values[1] + MONOTONE_ABS_TOL >= values[2],
        "{label} is not nonincreasing over tau={TAUS:?} within abs tol {MONOTONE_ABS_TOL:e}: values={values:?}"
    );
}

#[test]
fn shan_chen_coexistence_is_weakly_tau_dependent_and_monotone() {
    let points = TAUS.map(run_flat_interface);
    println!("SC coexistence tau sweep (G={G}, {NX}x{NY}, steps={STEPS})");
    println!("tau       rho_l          rho_v          rel_l_vs_tau1  rel_v_vs_tau1");
    let reference = points
        .iter()
        .find(|p| p.tau == 1.0)
        .expect("tau=1.0 reference missing");
    for p in &points {
        println!(
            "{:<8.3} {:<14.8} {:<14.8} {:<14.6e} {:<14.6e}",
            p.tau,
            p.rho_l,
            p.rho_v,
            rel_shift(p.rho_l, reference.rho_l),
            rel_shift(p.rho_v, reference.rho_v),
        );
    }

    for p in &points {
        let rho_l_shift = rel_shift(p.rho_l, reference.rho_l);
        let rho_v_shift = rel_shift(p.rho_v, reference.rho_v);
        assert!(
            rho_l_shift < 0.05,
            "P20 SC rho_l tau shift exceeded 5%: tau={:.3}, rho_l={:.8}, tau1={:.8}, rel={rho_l_shift:e}",
            p.tau,
            p.rho_l,
            reference.rho_l
        );
        assert!(
            rho_v_shift < 0.10,
            "P20 SC rho_v tau shift exceeded 10%: tau={:.3}, rho_v={:.8}, tau1={:.8}, rel={rho_v_shift:e}",
            p.tau,
            p.rho_v,
            reference.rho_v
        );
    }

    assert_nonincreasing(points.map(|p| p.rho_l), "rho_l");
    assert_nondecreasing(points.map(|p| p.rho_v), "rho_v");
}
