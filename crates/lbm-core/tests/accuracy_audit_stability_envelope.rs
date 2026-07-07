//! ACC-AUDIT lane 5.5: stability-envelope cartography.
//!
//! This file does not define new physics. It measures the documented
//! tune-stability envelope against short cavity sweeps:
//! - tau soft guideline: tau >= 0.55
//! - prescribed-speed/Mach proxy: U <= 0.15 soft, U <= 0.30 hard
//! - grid Reynolds guideline: U / nu <= 15
//!
//! The tests deliberately sample diagnostic lines instead of mapping the full
//! two-parameter divergence surface.

use lbm_core::compat::prelude::*;

const NX: usize = 64;
const NY: usize = 64;
const TRT: Collision = Collision::Trt {
    magic: Collision::MAGIC_STD,
};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Status {
    Pass,
    Diverged,
    OutOfEnvelope,
}

impl Status {
    fn label(self) -> &'static str {
        match self {
            Status::Pass => "PASS",
            Status::Diverged => "DIVERGED",
            Status::OutOfEnvelope => "OUT-OF-ENVELOPE",
        }
    }
}

#[derive(Clone, Copy, Debug)]
struct RunMeasurement {
    tau: f64,
    u_lid: f64,
    nu: f64,
    grid_re: f64,
    steps: usize,
    max_u: f64,
    finite: bool,
}

impl RunMeasurement {
    fn status_for(&self, in_documented_envelope: bool) -> Status {
        if !self.finite {
            Status::Diverged
        } else if in_documented_envelope {
            Status::Pass
        } else {
            Status::OutOfEnvelope
        }
    }
}

fn nu_from_tau(tau: f64) -> f64 {
    (tau - 0.5) / 3.0
}

fn tau_from_nu(nu: f64) -> f64 {
    3.0 * nu + 0.5
}

fn cavity_measurement(tau: f64, u_lid: f64, steps: usize) -> RunMeasurement {
    cavity_measurement_with_nu(tau, nu_from_tau(tau), u_lid, steps)
}

fn cavity_measurement_with_nu(
    reported_tau: f64,
    nu: f64,
    u_lid: f64,
    steps: usize,
) -> RunMeasurement {
    let mut sim: Simulation<f64> = SimConfig {
        nx: NX,
        ny: NY,
        nu,
        collision: TRT,
        edges: Edges {
            left: EdgeBC::BounceBack,
            right: EdgeBC::BounceBack,
            bottom: EdgeBC::BounceBack,
            top: EdgeBC::MovingWall { u: [u_lid, 0.0] },
        },
        ..Default::default()
    }
    .build()
    .unwrap();

    let mut max_u = 0.0f64;
    let mut finite = true;
    for _ in 0..steps {
        sim.step();
        let snapshot = scan_velocity(&sim);
        max_u = max_u.max(snapshot.0);
        finite &= snapshot.1;
        if !finite {
            break;
        }
    }

    RunMeasurement {
        tau: reported_tau,
        u_lid,
        nu,
        grid_re: u_lid / nu,
        steps,
        max_u,
        finite,
    }
}

fn scan_velocity(sim: &Simulation<f64>) -> (f64, bool) {
    let mut max_u = 0.0f64;
    let mut finite = true;
    for y in 0..sim.ny() {
        for x in 0..sim.nx() {
            let rho = sim.rho(x, y);
            let ux = sim.ux(x, y);
            let uy = sim.uy(x, y);
            finite &= rho.is_finite() && ux.is_finite() && uy.is_finite();
            max_u = max_u.max(ux.hypot(uy));
        }
    }
    (max_u, finite)
}

fn print_case(line: &str, measurement: RunMeasurement, status: Status) {
    println!(
        "ACC STAB {line}: tau={:.5}, Ma_proxy={:.5}, grid_Re={:.5}, nu={:.8e}, \
         steps={}, max|u|={:.8e}, label={}",
        measurement.tau,
        measurement.u_lid,
        measurement.grid_re,
        measurement.nu,
        measurement.steps,
        measurement.max_u,
        status.label()
    );
}

fn print_line_summary(line: &str, measurements: &[RunMeasurement]) {
    let n = measurements.len() as f64;
    let mean_tau = measurements.iter().map(|m| m.tau).sum::<f64>() / n;
    let mean_ma = measurements.iter().map(|m| m.u_lid).sum::<f64>() / n;
    let mean_grid_re = measurements.iter().map(|m| m.grid_re).sum::<f64>() / n;
    let max_u = measurements.iter().map(|m| m.max_u).fold(0.0f64, f64::max);
    println!(
        "ACC STAB {line} summary: mean_tau={mean_tau:.5}, mean_Ma_proxy={mean_ma:.5}, \
         mean_grid_Re={mean_grid_re:.5}, line_max|u|={max_u:.8e}"
    );
}

#[test]
fn g1_tau_line_fixed_ma_005_measures_soft_tau_boundary() {
    let taus = [0.51, 0.55, 0.60, 0.80, 1.00];
    let mut measurements = Vec::new();

    for tau in taus {
        let measurement = cavity_measurement(tau, 0.05, 1_000);
        let status = measurement.status_for(tau >= 0.55);
        print_case("G1 tau-line", measurement, status);
        measurements.push(measurement);

        if tau >= 0.55 {
            assert!(
                measurement.finite,
                "ACC STAB G1 tau-line: tau={tau:.5} is inside documented soft envelope \
                 tau>=0.55 but diverged; max|u|={:.8e}, grid_Re={:.5}",
                measurement.max_u, measurement.grid_re
            );
        }
    }

    let tau_051 = measurements
        .iter()
        .find(|m| (m.tau - 0.51).abs() < 1.0e-12)
        .unwrap();
    println!(
        "ACC STAB G1 tau=0.51 documented soft-limit finding: observed_label={} \
         max|u|={:.8e} grid_Re={:.5}; both finite and divergent outcomes are pinned",
        tau_051.status_for(false).label(),
        tau_051.max_u,
        tau_051.grid_re
    );
    print_line_summary("G1 tau-line", &measurements);
}

#[test]
fn g2_ma_line_fixed_tau_08_measures_soft_and_hard_speed_boundary() {
    let u_lids = [0.05, 0.10, 0.15, 0.20, 0.30];
    let mut measurements = Vec::new();

    for u_lid in u_lids {
        let measurement = cavity_measurement(0.80, u_lid, 500);
        let status = measurement.status_for(u_lid <= 0.15);
        print_case("G2 Ma-line", measurement, status);
        measurements.push(measurement);

        if u_lid <= 0.15 {
            assert!(
                measurement.finite,
                "ACC STAB G2 Ma-line: u_lid={u_lid:.5} is inside documented soft envelope \
                 U<=0.15 but diverged; max|u|={:.8e}, tau={:.5}",
                measurement.max_u, measurement.tau
            );
        }
    }

    let hard = measurements
        .iter()
        .find(|m| (m.u_lid - 0.30).abs() < 1.0e-12)
        .unwrap();
    assert!(
        hard.finite,
        "ACC STAB G2 Ma-line: hard-limit case u_lid=0.30 must remain finite; \
         max|u|={:.8e}, tau={:.5}, grid_Re={:.5}",
        hard.max_u, hard.tau, hard.grid_re
    );
    print_line_summary("G2 Ma-line", &measurements);
}

#[test]
fn g3_grid_re_line_low_speed_nu_sweep_stays_bounded() {
    let nus = [1.0 / 6.0, 1.0 / 12.0, 1.0 / 24.0, 1.0 / 48.0];
    let mut measurements = Vec::new();

    for nu in nus {
        let tau = tau_from_nu(nu);
        let measurement = cavity_measurement_with_nu(tau, nu, 0.05, 1_000);
        let status = measurement.status_for(measurement.grid_re <= 15.0);
        print_case("G3 grid-Re-line", measurement, status);
        measurements.push(measurement);

        assert!(
            measurement.finite,
            "ACC STAB G3 grid-Re-line: grid_Re={:.5} <= 15 should be stable; \
             tau={:.5}, nu={:.8e}, max|u|={:.8e}",
            measurement.grid_re, measurement.tau, measurement.nu, measurement.max_u
        );
    }

    let monotone_blow_up = measurements
        .windows(2)
        .all(|w| w[1].max_u > 1.25 * w[0].max_u)
        && measurements.last().unwrap().max_u > 2.0 * measurements.first().unwrap().max_u;
    assert!(
        !monotone_blow_up,
        "ACC STAB G3 grid-Re-line: max|u| shows monotone blow-up across all below-limit \
         grid_Re cases; values={:?}",
        measurements.iter().map(|m| m.max_u).collect::<Vec<_>>()
    );
    print_line_summary("G3 grid-Re-line", &measurements);
}

#[test]
#[ignore = "heavy lane 5.5 stability-envelope cartography: wider diagnostic sweeps"]
fn heavy_stability_envelope_cartography_extended_lines() {
    for tau in [0.505, 0.51, 0.525, 0.55, 0.575, 0.60, 0.70, 0.80, 1.00] {
        let measurement = cavity_measurement(tau, 0.05, 5_000);
        print_case(
            "HEAVY tau-line",
            measurement,
            measurement.status_for(tau >= 0.55),
        );
    }
    for u_lid in [0.05, 0.10, 0.15, 0.20, 0.25, 0.30] {
        let measurement = cavity_measurement(0.80, u_lid, 2_000);
        print_case(
            "HEAVY Ma-line",
            measurement,
            measurement.status_for(u_lid <= 0.15),
        );
    }
    for nu in [
        1.0 / 6.0,
        1.0 / 12.0,
        1.0 / 24.0,
        1.0 / 48.0,
        1.0 / 96.0,
        1.0 / 160.0,
    ] {
        let tau = tau_from_nu(nu);
        let measurement = cavity_measurement_with_nu(tau, nu, 0.05, 5_000);
        print_case(
            "HEAVY grid-Re-line",
            measurement,
            measurement.status_for(measurement.grid_re <= 15.0),
        );
    }
}
