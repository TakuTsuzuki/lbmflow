//! Validation T8: channel flow past a circular cylinder.

use lbm_core::prelude::*;

#[derive(Clone, Copy, Debug)]
struct CylinderCase {
    nx: usize,
    ny: usize,
    d: f64,
    cx: f64,
    cy: f64,
    u: f64,
    re: f64,
    steps: usize,
    sample_start: usize,
}

fn build_case(case: CylinderCase, outflow: bool) -> Simulation<f64> {
    let nu = case.u * case.d / case.re;
    let right = if outflow {
        EdgeBC::Outflow
    } else {
        EdgeBC::PressureOutlet { rho: 1.0 }
    };
    let mut sim: Simulation<f64> = SimConfig {
        nx: case.nx,
        ny: case.ny,
        nu,
        collision: Collision::Trt { magic: 3.0 / 16.0 },
        edges: Edges {
            left: EdgeBC::VelocityInlet { u: [case.u, 0.0] },
            right,
            bottom: EdgeBC::Periodic,
            top: EdgeBC::Periodic,
        },
        ..Default::default()
    }
    .build()
    .unwrap();
    sim.set_inlet_profile(Edge::Left, |_| [case.u, 0.0]);
    let r = 0.5 * case.d;
    let is_cylinder = |x: usize, y: usize| {
        let dx = x as f64 - case.cx;
        let dy = y as f64 - case.cy;
        dx * dx + dy * dy <= r * r
    };
    sim.set_solid_region(is_cylinder);
    sim.set_force_probe(is_cylinder);
    sim.init_with(|x, y| {
        if is_cylinder(x, y) || x == 0 || x == case.nx - 1 {
            (1.0, 0.0, 0.0)
        } else {
            // Tiny deterministic asymmetry shortens Re=100 shedding spin-up.
            let dy = (y as f64 - case.cy) / case.d;
            (1.0, case.u, 1.0e-4 * case.u * dy)
        }
    });
    sim
}

fn drag_lift(force: [f64; 2], case: CylinderCase) -> (f64, f64) {
    let scale = 2.0 / (case.u * case.u * case.d);
    (scale * force[0], scale * force[1])
}

#[test]
fn t8_re20_cylinder_steady_drag_is_in_reference_band() {
    let case = CylinderCase {
        nx: 440,
        ny: 160,
        d: 20.0,
        cx: 110.0,
        cy: 80.0,
        u: 0.05,
        re: 20.0,
        steps: 10_000,
        sample_start: 7_000,
    };
    let mut sim = build_case(case, false);
    let mut cd_sum = 0.0;
    let mut n = 0usize;
    for step in 0..case.steps {
        sim.step();
        if step >= case.sample_start {
            let (cd, _) = drag_lift(sim.probed_force(), case);
            cd_sum += cd;
            n += 1;
        }
    }
    let cd = cd_sum / n as f64;
    assert!(
        (1.8..=2.4).contains(&cd),
        "T8 Re=20 Cd = {cd:e}, steps = {}, samples = {n}",
        case.steps
    );
}

fn zero_crossing_frequency(samples: &[(usize, f64)]) -> Option<f64> {
    let mut crossings = Vec::new();
    for w in samples.windows(2) {
        let (t0, y0) = w[0];
        let (t1, y1) = w[1];
        if y0 == 0.0 || (y0 < 0.0) != (y1 < 0.0) {
            let frac = y0.abs() / (y0.abs() + y1.abs());
            crossings.push(t0 as f64 + frac * (t1 - t0) as f64);
        }
    }
    if crossings.len() < 5 {
        return None;
    }
    let periods: Vec<f64> = crossings.windows(3).map(|w| w[2] - w[0]).collect();
    Some(1.0 / (periods.iter().sum::<f64>() / periods.len() as f64))
}

#[test]
#[ignore]
fn t8_re100_cylinder_vortex_shedding_has_expected_st_cd_cl() {
    let case = CylinderCase {
        nx: 440,
        ny: 160,
        d: 20.0,
        cx: 110.0,
        // Off-centre by one lattice cell to trigger shedding faster.
        cy: 81.0,
        u: 0.05,
        re: 100.0,
        steps: 80_000,
        sample_start: 40_000,
    };
    let mut sim = build_case(case, true);
    let mut cd_sum = 0.0;
    let mut cl_min = f64::INFINITY;
    let mut cl_max = f64::NEG_INFINITY;
    let mut cl_samples = Vec::new();
    for step in 0..case.steps {
        if step < 5_000 {
            let ramp = case.u * (step as f64 + 1.0) / 5_000.0;
            sim.set_inlet_profile(Edge::Left, |_| [ramp, 0.0]);
        }
        sim.step();
        if step >= case.sample_start {
            let (cd, cl) = drag_lift(sim.probed_force(), case);
            cd_sum += cd;
            cl_min = cl_min.min(cl);
            cl_max = cl_max.max(cl);
            cl_samples.push((step, cl));
        }
    }
    let cd = cd_sum / cl_samples.len() as f64;
    let amp = 0.5 * (cl_max - cl_min);
    let freq =
        zero_crossing_frequency(&cl_samples).expect("T8 Re=100 not enough Cl zero crossings");
    let st = freq * case.d / case.u;
    assert!(
        (0.15..=0.19).contains(&st),
        "T8 Re=100 St = {st:e}, freq = {freq:e}, zero_crossings_samples = {}",
        cl_samples.len()
    );
    assert!(
        (1.2..=1.5).contains(&cd),
        "T8 Re=100 mean Cd = {cd:e}, St = {st:e}, Cl amplitude = {amp:e}"
    );
    assert!(
        (0.2..=0.45).contains(&amp),
        "T8 Re=100 Cl amplitude = {amp:e}, St = {st:e}, mean Cd = {cd:e}, cl_min = {cl_min:e}, cl_max = {cl_max:e}"
    );
}
