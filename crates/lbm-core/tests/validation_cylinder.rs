// Inherited verbatim from the retired V1 suite at its retirement (2026-07-05,
// scripts/sync-tests.sh mechanical retarget); now the canonical facade tests.
//! Validation T8: Schaefer-Turek channel flow past a circular cylinder.

use lbm_core::compat::prelude::*;

const MAGIC: f64 = 3.0 / 16.0;
const RHO: f64 = 1.0;
const CD_REF_2D1: f64 = 5.5795;

#[derive(Clone, Copy, Debug)]
struct CylinderCase {
    nx: usize,
    ny: usize,
    d: f64,
    cx: f64,
    cy: f64,
    u_max: f64,
    nu: f64,
    steps: usize,
    sample_start: usize,
    include_radius_boundary: bool,
}

impl CylinderCase {
    fn u_mean(self) -> f64 {
        (2.0 / 3.0) * self.u_max
    }

    fn re(self) -> f64 {
        self.u_mean() * self.d / self.nu
    }

    fn height(self) -> f64 {
        (self.ny - 2) as f64
    }
}

fn schaefer_turek_2d1_d20() -> CylinderCase {
    CylinderCase {
        // Wall rims occupy y=0 and y=ny-1, so H=ny-2=80=4D.
        // The benchmark ratio is 4.1D; this lattice keeps the documented
        // grid 440x82 and places the center 2D from inlet and lower wall.
        // With U_mean=(2/3)u_max, nu=0.05 gives the specified Re=20.
        nx: 440,
        ny: 82,
        d: 20.0,
        cx: 40.0,
        cy: 40.0,
        u_max: 0.075,
        nu: 0.05,
        steps: 30_000,
        sample_start: 20_000,
        include_radius_boundary: true,
    }
}

fn schaefer_turek_2d1_d40() -> CylinderCase {
    CylinderCase {
        // H=162 approximates 4.1D, and center (80,80) preserves the
        // inlet offset 2D and lower-wall offset 2D from the benchmark.
        nx: 880,
        ny: 164,
        d: 40.0,
        cx: 80.0,
        cy: 80.0,
        u_max: 0.075,
        nu: 0.1,
        steps: 50_000,
        sample_start: 35_000,
        include_radius_boundary: false,
    }
}

fn schaefer_turek_2d2_d40() -> CylinderCase {
    CylinderCase {
        nx: 880,
        ny: 164,
        d: 40.0,
        cx: 80.0,
        cy: 80.0,
        u_max: 0.15,
        nu: 0.04,
        steps: 120_000,
        sample_start: 80_000,
        include_radius_boundary: false,
    }
}

fn inlet_velocity(case: CylinderCase, y: usize) -> [f64; 2] {
    if y == 0 || y == case.ny - 1 {
        return [0.0, 0.0];
    }
    let h = case.height();
    let y_w = y as f64 - 0.5;
    [4.0 * case.u_max * y_w * (h - y_w) / (h * h), 0.0]
}

fn build_case(case: CylinderCase) -> Simulation<f64> {
    build_case_with_wall(case, false)
}

fn build_case_with_wall(case: CylinderCase, bouzidi: bool) -> Simulation<f64> {
    let mut sim: Simulation<f64> = SimConfig {
        nx: case.nx,
        ny: case.ny,
        nu: case.nu,
        collision: Collision::Trt { magic: MAGIC },
        edges: Edges {
            left: EdgeBC::VelocityInlet {
                u: [case.u_max, 0.0],
            },
            right: EdgeBC::PressureOutlet { rho: RHO },
            bottom: EdgeBC::BounceBack,
            top: EdgeBC::BounceBack,
        },
        ..Default::default()
    }
    .build()
    .unwrap();
    sim.set_inlet_profile(Edge::Left, |y| inlet_velocity(case, y));
    let r = 0.5 * case.d;
    let is_cylinder = |x: usize, y: usize| {
        let dx = x as f64 - case.cx;
        let dy = y as f64 - case.cy;
        let r2 = r * r;
        let d2 = dx * dx + dy * dy;
        if case.include_radius_boundary {
            d2 <= r2
        } else {
            d2 < r2
        }
    };
    sim.set_solid_region(is_cylinder);
    if bouzidi {
        sim.set_bouzidi_circle(case.cx, case.cy, r);
    }
    sim.set_force_probe(is_cylinder);
    sim.init_with(|x, y| {
        if is_cylinder(x, y) || x == 0 || x == case.nx - 1 || y == 0 || y == case.ny - 1 {
            (RHO, 0.0, 0.0)
        } else {
            let u = inlet_velocity(case, y)[0];
            let dy = (y as f64 - case.cy) / case.d;
            (RHO, u, 1.0e-5 * case.u_max * dy)
        }
    });
    sim
}

fn drag_lift(force: [f64; 2], case: CylinderCase) -> (f64, f64) {
    let scale = 2.0 / (RHO * case.u_mean() * case.u_mean() * case.d);
    (scale * force[0], scale * force[1])
}

fn run_steady_cd_cl(case: CylinderCase) -> (f64, f64, usize) {
    run_steady_cd_cl_with_wall(case, false)
}

fn run_steady_cd_cl_with_wall(case: CylinderCase, bouzidi: bool) -> (f64, f64, usize) {
    let mut sim = build_case_with_wall(case, bouzidi);
    let mut cd_sum = 0.0;
    let mut cl_sum = 0.0;
    let mut n = 0usize;
    for step in 0..case.steps {
        sim.step();
        if step >= case.sample_start {
            let (cd, cl) = drag_lift(sim.probed_force(), case);
            cd_sum += cd;
            cl_sum += cl;
            n += 1;
        }
    }
    (cd_sum / n as f64, cl_sum / n as f64, n)
}

#[test]
fn t8_2d1_d20_cylinder_steady_drag_lift_are_in_reference_bands() {
    let case = schaefer_turek_2d1_d20();
    let (cd, cl, samples) = run_steady_cd_cl(case);
    assert!(
        (5.2..=6.0).contains(&cd),
        "T8 2D-1 D=20 Cd = {cd:e}, Cl = {cl:e}, Re = {:e}, steps = {}, samples = {samples}",
        case.re(),
        case.steps
    );
    assert!(
        (-0.05..=0.08).contains(&cl),
        "T8 2D-1 D=20 Cl = {cl:e}, Cd = {cd:e}, Re = {:e}, steps = {}, samples = {samples}",
        case.re(),
        case.steps
    );
}

#[test]
#[ignore = "Bouzidi D=20 characterization: explicit curved-wall run"]
fn t8_bouzidi_2d1_d20_cylinder_steady_drag_lift_characterization() {
    let case = schaefer_turek_2d1_d20();
    let (cd, cl, samples) = run_steady_cd_cl_with_wall(case, true);
    println!(
        "T8 Bouzidi 2D-1 D=20 measured Cd={cd:.8}, Cl={cl:.8}, Re={:.8}, samples={samples}",
        case.re()
    );
    assert!(
        (5.80..=5.87).contains(&cd),
        "T8 Bouzidi D=20 characterization Cd = {cd:e}, Cl = {cl:e}, Re = {:e}, steps = {}, samples = {samples}",
        case.re(),
        case.steps
    );
    assert!(
        (-0.03..=0.05).contains(&cl),
        "T8 Bouzidi D=20 characterization Cl = {cl:e}, Cd = {cd:e}, Re = {:e}, steps = {}, samples = {samples}",
        case.re(),
        case.steps
    );
}

#[test]
#[ignore]
fn t8_2d1_d40_cylinder_drag_converges_toward_reference() {
    let coarse = schaefer_turek_2d1_d20();
    let fine = schaefer_turek_2d1_d40();
    let (cd20, cl20, samples20) = run_steady_cd_cl(coarse);
    let (cd40, cl40, samples40) = run_steady_cd_cl(fine);
    assert!(
        (5.35..=5.85).contains(&cd40),
        "T8 2D-1 D=40 Cd = {cd40:e}, Cl = {cl40:e}, Re = {:e}, steps = {}, samples = {samples40}",
        fine.re(),
        fine.steps
    );
    let err20 = (cd20 - CD_REF_2D1).abs();
    let err40 = (cd40 - CD_REF_2D1).abs();
    assert!(
        err40 < err20,
        "T8 2D-1 convergence err40 = {err40:e}, err20 = {err20:e}, Cd40 = {cd40:e}, Cd20 = {cd20:e}, Cl40 = {cl40:e}, Cl20 = {cl20:e}, samples40 = {samples40}, samples20 = {samples20}"
    );
}

fn zero_crossing_periods(samples: &[(usize, f64)]) -> Vec<f64> {
    let mut crossings = Vec::new();
    for w in samples.windows(2) {
        let (t0, y0) = w[0];
        let (t1, y1) = w[1];
        if y0 == 0.0 || (y0 < 0.0) != (y1 < 0.0) {
            let frac = y0.abs() / (y0.abs() + y1.abs());
            crossings.push(t0 as f64 + frac * (t1 - t0) as f64);
        }
    }
    crossings.windows(3).map(|w| w[2] - w[0]).collect()
}

#[test]
#[ignore]
fn t8_2d2_d40_cylinder_vortex_shedding_matches_reference_bands() {
    let case = schaefer_turek_2d2_d40();
    let mut sim = build_case(case);
    let mut cd_max = f64::NEG_INFINITY;
    let mut cl_abs_max = 0.0f64;
    let mut cl_samples = Vec::new();
    for step in 0..case.steps {
        sim.step();
        if step >= case.sample_start {
            let (cd, cl) = drag_lift(sim.probed_force(), case);
            cd_max = cd_max.max(cd);
            cl_abs_max = cl_abs_max.max(cl.abs());
            cl_samples.push((step, cl));
        }
    }
    let periods = zero_crossing_periods(&cl_samples);
    assert!(
        periods.len() >= 3,
        "T8 2D-2 D=40 zero-crossing periods = {}, samples = {}, steps = {}",
        periods.len(),
        cl_samples.len(),
        case.steps
    );
    let mean_period = periods.iter().sum::<f64>() / periods.len() as f64;
    let period_spread = periods
        .iter()
        .map(|p| (p - mean_period).abs() / mean_period)
        .fold(0.0, f64::max);
    let st = case.d / (case.u_mean() * mean_period);
    assert!(
        (0.28..=0.32).contains(&st),
        "T8 2D-2 D=40 St = {st:e}, mean_period = {mean_period:e}, periods = {}, steps = {}",
        periods.len(),
        case.steps
    );
    assert!(
        (3.0..=3.5).contains(&cd_max),
        "T8 2D-2 D=40 Cd_max = {cd_max:e}, St = {st:e}, Cl_max = {cl_abs_max:e}, steps = {}",
        case.steps
    );
    assert!(
        (0.8..=1.2).contains(&cl_abs_max),
        "T8 2D-2 D=40 Cl_max = {cl_abs_max:e}, St = {st:e}, Cd_max = {cd_max:e}, steps = {}",
        case.steps
    );
    assert!(
        period_spread <= 0.02,
        "T8 2D-2 D=40 period spread = {period_spread:e}, mean_period = {mean_period:e}, periods = {:?}",
        periods
    );
}
