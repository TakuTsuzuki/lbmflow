//! ACC-AUDIT radar row 18 / pitfall #11: Bouzidi mass-drift pin.
//!
//! Derivation of the monitored quantity:
//! - The discrete total mass is M(t) = sum_{fluid cells} rho_i(t).
//! - Ordinary half-way bounce-back moves each outgoing wall-link population
//!   back to the opposite population of the same adjacent fluid node. That is
//!   a one-to-one permutation/replacement of populations on closed walls, so
//!   with periodic x and closed y walls it conserves M up to floating-point
//!   summation roundoff.
//! - Bouzidi-Firdaouss-Lallemand interpolated bounce-back instead reconstructs
//!   the reflected population from a linear interpolation along an off-grid
//!   link. For qd < 1/2 the update combines the current cell and its next
//!   fluid neighbor; for qd >= 1/2 it combines the streamed population and a
//!   wall reference term. Those interpolation weights improve curved-wall
//!   placement, but they are not a global conservative redistribution over the
//!   cut-cell boundary, so a curved Bouzidi wall can slowly leak or gain mass.
//! - This test therefore does not claim exact conservation for Bouzidi. It
//!   freezes the current relative drift
//!       |M(t) - M(0)| / M(0)
//!   over 10^4 steps, and keeps a staircase half-way bounce-back run as the
//!   conservative control.

use lbm_core::compat::prelude::*;

const N: usize = 64;
const RADIUS: f64 = 8.0;
const CENTER: f64 = 32.0;
const STEPS: usize = 10_000;
const SAMPLE_EVERY: usize = 1_000;
const FORCE: [f64; 2] = [1.0e-6, 0.0];
const NU_FOR_TAU_0_8: f64 = 0.1;
const STAIRCASE_DRIFT_BAND: f64 = 1.0e-11;
const BOUZIDI_DRIFT_PIN: f64 = 1.2e-8;
const BOUZIDI_FINDING_LINE: f64 = 1.0e-6;
const TRT: Collision = Collision::Trt {
    magic: Collision::MAGIC_STD,
};

#[derive(Clone, Copy)]
enum CylinderWall {
    Staircase,
    Bouzidi,
}

fn inside_cylinder(x: usize, y: usize) -> bool {
    let dx = x as f64 - CENTER;
    let dy = y as f64 - CENTER;
    dx * dx + dy * dy <= RADIUS * RADIUS
}

fn build_case(wall: CylinderWall) -> Simulation<f64> {
    let mut sim: Simulation<f64> = SimConfig {
        nx: N,
        ny: N,
        nu: NU_FOR_TAU_0_8,
        collision: TRT,
        force: FORCE,
        edges: Edges {
            left: EdgeBC::Periodic,
            right: EdgeBC::Periodic,
            bottom: EdgeBC::BounceBack,
            top: EdgeBC::BounceBack,
        },
        ..Default::default()
    }
    .build()
    .unwrap();

    sim.set_solid_region(inside_cylinder);
    if matches!(wall, CylinderWall::Bouzidi) {
        sim.set_bouzidi_circle(CENTER, CENTER, RADIUS);
    }
    sim
}

fn mass_drift_series(wall: CylinderWall) -> Vec<(usize, f64, f64)> {
    let mut sim = build_case(wall);
    let m0 = sim.total_mass_f64();
    let mut out = vec![(0, m0, 0.0)];
    for step in (SAMPLE_EVERY..=STEPS).step_by(SAMPLE_EVERY) {
        sim.run(SAMPLE_EVERY);
        let mass = sim.total_mass_f64();
        let drift = (mass - m0).abs() / m0.abs();
        out.push((step, mass, drift));
    }
    out
}

fn print_series(label: &str, series: &[(usize, f64, f64)]) {
    for &(step, mass, drift) in series {
        println!("{label} step={step:5} mass={mass:.16e} rel_drift={drift:.12e}");
    }
}

fn max_drift(series: &[(usize, f64, f64)]) -> f64 {
    series
        .iter()
        .map(|&(_, _, drift)| drift)
        .fold(0.0, f64::max)
}

#[test]
fn radar18_pitfall11_bouzidi_cylinder_mass_drift_is_pinned() {
    let staircase = mass_drift_series(CylinderWall::Staircase);
    let bouzidi = mass_drift_series(CylinderWall::Bouzidi);
    print_series("staircase", &staircase);
    print_series("bouzidi  ", &bouzidi);

    let staircase_max = max_drift(&staircase);
    let bouzidi_max = max_drift(&bouzidi);
    println!(
        "radar #18 pitfall #11 mass drift summary: staircase_max={staircase_max:.12e}, bouzidi_max={bouzidi_max:.12e}, force=[{:.1e},{:.1e}], tau=0.8, grid={N}x{N}, R={RADIUS}",
        FORCE[0], FORCE[1]
    );

    assert!(
        staircase_max <= STAIRCASE_DRIFT_BAND,
        "radar #18 pitfall #11 control: staircase half-way bounce-back should be mass-conservative; max relative drift={staircase_max:.12e}, band={STAIRCASE_DRIFT_BAND:.1e}, denominator=M(0)"
    );
    assert!(
        bouzidi_max <= BOUZIDI_DRIFT_PIN,
        "radar #18 pitfall #11 Bouzidi mass-drift pin: max relative drift={bouzidi_max:.12e}, band={BOUZIDI_DRIFT_PIN:.1e}, denominator=M(0). This monitors current curved-wall interpolation leakage; drift > {BOUZIDI_FINDING_LINE:.1e} would be a real Bouzidi mass-loss finding beyond the expected pin."
    );
}
