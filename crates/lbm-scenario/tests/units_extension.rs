use lbm_core::compat::prelude::ConfigError;
use lbm_scenario::{
    build, build_check, CollisionSpec, EdgeSpec, EdgesSpec, FlowParams, Grid, InitSpec,
    LegacyUnitReport, Physics, Precision, RunSpec, Scenario, SimHandle, UnitConstructor,
};

const ROUND_TRIP_REL_TOL: f64 = 1.0e-6;
const GRAVITY_REL_TOL: f64 = 1.0e-10;
const FIELD_L2_REL_TOL: f64 = 1.0e-12;

#[derive(Clone, Copy, Debug)]
struct PhysicalTuple {
    name: &'static str,
    nu_m2s: f64,
    l_m: f64,
    u_m_s: f64,
    rho_kg_m3: f64,
    g_m_s2: f64,
}

#[derive(Clone, Debug)]
struct FieldSnapshot {
    rho: Vec<f64>,
    ux: Vec<f64>,
    uy: Vec<f64>,
    near_lid_ux: f64,
}

fn physical_tuples() -> [PhysicalTuple; 10] {
    [
        PhysicalTuple {
            name: "water_channel",
            nu_m2s: 1.0e-6,
            l_m: 1.0e-2,
            u_m_s: 1.0e-2,
            rho_kg_m3: 998.2,
            g_m_s2: 9.81,
        },
        PhysicalTuple {
            name: "room_air_duct",
            nu_m2s: 1.5e-5,
            l_m: 1.0e-1,
            u_m_s: 1.0,
            rho_kg_m3: 1.204,
            g_m_s2: 9.81,
        },
        PhysicalTuple {
            name: "glycerol_creep",
            nu_m2s: 1.0e-3,
            l_m: 5.0e-2,
            u_m_s: 2.0e-2,
            rho_kg_m3: 1_260.0,
            g_m_s2: 9.81,
        },
        PhysicalTuple {
            name: "microfluidic_water",
            nu_m2s: 1.0e-6,
            l_m: 1.0e-4,
            u_m_s: 1.0e-3,
            rho_kg_m3: 998.2,
            g_m_s2: 9.81,
        },
        PhysicalTuple {
            name: "light_oil",
            nu_m2s: 1.0e-5,
            l_m: 2.0e-2,
            u_m_s: 1.0e-1,
            rho_kg_m3: 850.0,
            g_m_s2: 9.81,
        },
        PhysicalTuple {
            name: "blood_scale",
            nu_m2s: 3.5e-6,
            l_m: 3.0e-3,
            u_m_s: 5.0e-2,
            rho_kg_m3: 1_060.0,
            g_m_s2: 9.81,
        },
        PhysicalTuple {
            name: "nanofluidic_low_nu",
            nu_m2s: 1.0e-9,
            l_m: 1.0e-6,
            u_m_s: 1.0e-2,
            rho_kg_m3: 997.0,
            g_m_s2: 9.81,
        },
        PhysicalTuple {
            name: "slow_mixer",
            nu_m2s: 5.0e-6,
            l_m: 5.0e-1,
            u_m_s: 2.0e-2,
            rho_kg_m3: 1_010.0,
            g_m_s2: 9.81,
        },
        PhysicalTuple {
            name: "warm_gas",
            nu_m2s: 4.0e-5,
            l_m: 2.0e-1,
            u_m_s: 5.0e-1,
            rho_kg_m3: 0.95,
            g_m_s2: 9.81,
        },
        PhysicalTuple {
            name: "viscous_syrup",
            nu_m2s: 5.0e-2,
            l_m: 1.0e-1,
            u_m_s: 1.0e-3,
            rho_kg_m3: 1_380.0,
            g_m_s2: 9.81,
        },
    ]
}

fn base_params(case: PhysicalTuple, constructor: UnitConstructor) -> FlowParams {
    FlowParams {
        constructor,
        characteristic_length: case.l_m,
        characteristic_velocity: case.u_m_s,
        kinematic_viscosity: case.nu_m2s,
        density: Some(case.rho_kg_m3),
        resolution: None,
        lattice_velocity: None,
        relaxation_time: None,
        end_time: None,
        end_step_count: None,
        gravity: Some([0.0, -case.g_m_s2]),
        reference_pressure: Some(0.0),
        re_physical: None,
    }
}

fn constructor_params(
    case: PhysicalTuple,
    constructor: UnitConstructor,
    resolution: usize,
    lattice_velocity: f64,
    tau: f64,
) -> FlowParams {
    let mut params = base_params(case, constructor);
    match constructor {
        UnitConstructor::FromResolutionAndLatticeVelocity => {
            params.resolution = Some(resolution);
            params.lattice_velocity = Some(lattice_velocity);
        }
        UnitConstructor::FromResolutionAndRelaxationTime => {
            params.resolution = Some(resolution);
            params.relaxation_time = Some(tau);
        }
        UnitConstructor::FromRelaxationTimeAndLatticeVelocity => {
            params.relaxation_time = Some(tau);
            params.lattice_velocity = Some(lattice_velocity);
        }
    }
    params
}

fn report_for(params: &FlowParams) -> LegacyUnitReport {
    lbm_scenario::unit_report(params).unwrap_or_else(|err| panic!("unit report failed: {err}"))
}

fn resolved_resolution(report: &LegacyUnitReport) -> usize {
    (report.inputs.characteristic_length / report.conversion_factors.length_m).round() as usize
}

fn assert_rel_close(label: &str, actual: f64, expected: f64, rel_tol: f64) {
    let scale = actual.abs().max(expected.abs()).max(1.0);
    let abs = (actual - expected).abs();
    assert!(
        abs <= rel_tol * scale,
        "{label}: actual={actual:e} expected={expected:e} abs={abs:e} rel_limit={:e}",
        rel_tol * scale
    );
}

fn assert_round_trips(case: PhysicalTuple, report: &LegacyUnitReport) {
    let n = resolved_resolution(report);
    let l_back = n as f64 * report.conversion_factors.length_m;
    let u_back = report.lattice.u_char_lattice * report.conversion_factors.velocity_m_s;
    let nu_back = report.lattice.nu_lattice * report.conversion_factors.viscosity_m2_s;
    let g_back = report.lattice.g_lat.unwrap()[1] * report.conversion_factors.acceleration_m_s2;

    assert_rel_close("L round-trip", l_back, case.l_m, ROUND_TRIP_REL_TOL);
    assert_rel_close("U round-trip", u_back, case.u_m_s, ROUND_TRIP_REL_TOL);
    assert_rel_close("nu round-trip", nu_back, case.nu_m2s, ROUND_TRIP_REL_TOL);
    assert_rel_close("g round-trip", g_back, -case.g_m_s2, ROUND_TRIP_REL_TOL);
}

fn canonical_scenario(units: FlowParams, lid_u_lattice: f64) -> Scenario {
    Scenario {
        version: 0,
        name: format!("units-consistency-{:?}", units.constructor),
        grid: Grid {
            nx: 32,
            ny: 32,
            nz: 1,
            lattice: None,
        },
        physics: Physics {
            nu: 0.1,
            collision: CollisionSpec::Trt,
            force: [0.0, 0.0],
            gravity: None,
            precision: Precision::F64,
        },
        units: Some(units),
        compute: None,
        wall: None,
        edges: EdgesSpec {
            left: EdgeSpec::BounceBack,
            right: EdgeSpec::BounceBack,
            bottom: EdgeSpec::BounceBack,
            top: EdgeSpec::MovingWall {
                u: [lid_u_lattice, 0.0],
            },
            front: None,
            back: None,
        },
        inlet_profile: None,
        obstacles: Vec::new(),
        init: InitSpec::Rest,
        multiphase: None,
        rotor: None,
        particles: None,
        run: RunSpec {
            steps: 100,
            stop_when_steady: None,
        },
        probes: Vec::new(),
        outputs: Vec::new(),
    }
}

fn run_snapshot(sc: &Scenario) -> FieldSnapshot {
    match build(sc).unwrap_or_else(|err| panic!("build failed for {}: {err}", sc.name)) {
        SimHandle::F64(mut sim, None) => {
            sim.run(sc.run.steps);
            FieldSnapshot {
                rho: sim.rho_field().to_vec(),
                ux: sim.ux_field().to_vec(),
                uy: sim.uy_field().to_vec(),
                near_lid_ux: sim.ux(16, 30),
            }
        }
        SimHandle::F64(_, Some(_)) => panic!("canonical scenario unexpectedly built multiphase"),
        SimHandle::F32(_, _) => panic!("canonical scenario should build as f64"),
    }
}

fn l2_rel(a: &FieldSnapshot, b: &FieldSnapshot) -> f64 {
    let mut numerator = 0.0;
    let mut denominator = 0.0;
    for (left, right) in a
        .rho
        .iter()
        .chain(a.ux.iter())
        .chain(a.uy.iter())
        .zip(b.rho.iter().chain(b.ux.iter()).chain(b.uy.iter()))
    {
        let diff = left - right;
        numerator += diff * diff;
        denominator += left * left;
    }
    numerator.sqrt() / denominator.sqrt().max(1.0e-300)
}

#[test]
fn si_lattice_si_round_trip_across_constructors() {
    let resolution = 192;
    let lattice_velocity = 0.05;
    for case in physical_tuples() {
        let re = case.u_m_s * case.l_m / case.nu_m2s;
        let tau = 3.0 * lattice_velocity * resolution as f64 / re + 0.5;
        for constructor in [
            UnitConstructor::FromResolutionAndLatticeVelocity,
            UnitConstructor::FromResolutionAndRelaxationTime,
            UnitConstructor::FromRelaxationTimeAndLatticeVelocity,
        ] {
            let params = constructor_params(case, constructor, resolution, lattice_velocity, tau);
            let report = report_for(&params);
            assert_round_trips(case, &report);
            println!(
                "round_trip case={} constructor={:?} Re={:.12e} N={} u_lat={:.12e} tau={:.12e} nu_lat={:.12e} dx={:.12e} dt={:.12e} g_lat_y={:.12e} verdict={:?}",
                case.name,
                constructor,
                report.dimensionless.reynolds,
                resolved_resolution(&report),
                report.lattice.u_char_lattice,
                report.lattice.tau,
                report.lattice.nu_lattice,
                report.lattice.dx,
                report.lattice.dt,
                report.lattice.g_lat.unwrap()[1],
                report.verdict,
            );
        }
    }
}

#[test]
fn parameter_range_edges_are_reported_without_silent_drift() {
    let low_nu = PhysicalTuple {
        name: "water_25c_low_nu_edge",
        nu_m2s: 1.0e-9,
        l_m: 1.0e-2,
        u_m_s: 1.0e-4,
        rho_kg_m3: 997.0,
        g_m_s2: 9.81,
    };
    let low_nu_report = report_for(&constructor_params(
        low_nu,
        UnitConstructor::FromResolutionAndLatticeVelocity,
        512,
        0.01,
        0.0,
    ));
    assert_round_trips(low_nu, &low_nu_report);
    assert!(
        low_nu_report.lattice.tau > 0.5,
        "low-nu edge must retain positive lattice viscosity"
    );
    println!(
        "edge low_nu case={} Re={:.12e} tau={:.12e} nu_lat={:.12e} verdict={:?}",
        low_nu.name,
        low_nu_report.dimensionless.reynolds,
        low_nu_report.lattice.tau,
        low_nu_report.lattice.nu_lattice,
        low_nu_report.verdict
    );

    let nanofluidics = PhysicalTuple {
        name: "one_micron_nanofluidics",
        nu_m2s: 1.0e-9,
        l_m: 1.0e-6,
        u_m_s: 1.0e-2,
        rho_kg_m3: 997.0,
        g_m_s2: 9.81,
    };
    let nanofluidics_report = report_for(&constructor_params(
        nanofluidics,
        UnitConstructor::FromResolutionAndLatticeVelocity,
        100,
        0.05,
        0.0,
    ));
    assert_round_trips(nanofluidics, &nanofluidics_report);
    println!(
        "edge nanofluidics case={} Re={:.12e} dx={:.12e} dt={:.12e} tau={:.12e} verdict={:?}",
        nanofluidics.name,
        nanofluidics_report.dimensionless.reynolds,
        nanofluidics_report.lattice.dx,
        nanofluidics_report.lattice.dt,
        nanofluidics_report.lattice.tau,
        nanofluidics_report.verdict
    );

    let aeroacoustic = PhysicalTuple {
        name: "aeroacoustic_velocity_limit",
        nu_m2s: 1.5e-5,
        l_m: 1.0,
        u_m_s: 1.0e2,
        rho_kg_m3: 1.204,
        g_m_s2: 9.81,
    };
    let mut too_fast = canonical_scenario(
        constructor_params(
            aeroacoustic,
            UnitConstructor::FromResolutionAndRelaxationTime,
            100,
            0.0,
            0.56,
        ),
        0.31,
    );
    too_fast.units = None;
    too_fast.name = "aeroacoustic-too-fast-wall".to_string();
    match build(&too_fast) {
        Err(lbm_scenario::BuildError::Core(ConfigError::VelocityTooHigh { speed })) => {
            println!(
                "edge aeroacoustic case={} U={:.12e} wall_u_lat={:.12e} error=VelocityTooHigh",
                aeroacoustic.name, aeroacoustic.u_m_s, speed
            );
            assert_rel_close("too-fast wall speed echo", speed, 0.31, 1.0e-15);
        }
        Ok(_) => panic!("expected VelocityTooHigh for aeroacoustic edge, got successful build"),
        Err(err) => panic!("expected VelocityTooHigh for aeroacoustic edge, got {err}"),
    }
    let check_err = build_check(&too_fast).expect_err("build_check must reject too-fast wall");
    assert!(
        check_err.contains("low-Mach limit"),
        "unexpected build_check error: {check_err}"
    );
}

#[test]
fn gravity_conversion_anchor_matches_derived_formula() {
    let case = PhysicalTuple {
        name: "gravity_anchor",
        nu_m2s: 1.0e-6,
        l_m: 1.0,
        u_m_s: 0.1,
        rho_kg_m3: 998.2,
        g_m_s2: 9.81,
    };
    let resolution = 100;
    let lattice_velocity = 0.05;
    let re = case.u_m_s * case.l_m / case.nu_m2s;
    let tau = 3.0 * lattice_velocity * resolution as f64 / re + 0.5;
    let report = report_for(&constructor_params(
        case,
        UnitConstructor::FromResolutionAndLatticeVelocity,
        resolution,
        lattice_velocity,
        tau,
    ));

    // Derivation:
    //   dx = L / N
    //   u_lat = U * dt / dx, so dt = u_lat * dx / U
    //   g_lat = g * dt^2 / dx
    //         = g * u_lat^2 * dx / U^2
    //         = g * u_lat^2 * L / (N * U^2)
    let expected_g_lat = case.g_m_s2 * lattice_velocity.powi(2) * case.l_m
        / (resolution as f64 * case.u_m_s.powi(2));
    let actual_g_lat = -report.lattice.g_lat.unwrap()[1];
    assert_rel_close(
        "gravity lattice acceleration",
        actual_g_lat,
        expected_g_lat,
        GRAVITY_REL_TOL,
    );
    println!(
        "gravity_anchor g={:.12e} L={:.12e} U={:.12e} N={} u_lat={:.12e} tau={:.12e} expected_g_lat={:.12e} actual_g_lat={:.12e}",
        case.g_m_s2,
        case.l_m,
        case.u_m_s,
        resolution,
        lattice_velocity,
        report.lattice.tau,
        expected_g_lat,
        actual_g_lat
    );
}

#[test]
fn constructors_produce_identical_hundred_step_behavior() {
    let case = PhysicalTuple {
        name: "canonical_cavity",
        nu_m2s: 2.0e-4,
        l_m: 0.1,
        u_m_s: 0.2,
        rho_kg_m3: 998.2,
        g_m_s2: 9.81,
    };
    let resolution = 100;
    let lattice_velocity = 0.08;
    let re = case.u_m_s * case.l_m / case.nu_m2s;
    let tau = 3.0 * lattice_velocity * resolution as f64 / re + 0.5;

    let scenarios = [
        UnitConstructor::FromResolutionAndLatticeVelocity,
        UnitConstructor::FromResolutionAndRelaxationTime,
        UnitConstructor::FromRelaxationTimeAndLatticeVelocity,
    ]
    .map(|constructor| {
        let params = constructor_params(case, constructor, resolution, lattice_velocity, tau);
        let report = report_for(&params);
        assert_eq!(
            resolved_resolution(&report),
            resolution,
            "constructor {:?} changed the canonical resolution",
            constructor
        );
        println!(
            "behavior_setup constructor={:?} N={} u_lat={:.12e} tau={:.12e} nu_lat={:.12e} verdict={:?}",
            constructor,
            resolved_resolution(&report),
            report.lattice.u_char_lattice,
            report.lattice.tau,
            report.lattice.nu_lattice,
            report.verdict
        );
        canonical_scenario(params, lattice_velocity)
    });

    let snapshots = scenarios.each_ref().map(run_snapshot);
    assert!(
        snapshots[0].near_lid_ux > 0.0,
        "moving lid should create positive near-lid x velocity, got {:.12e}",
        snapshots[0].near_lid_ux
    );
    for (idx, snapshot) in snapshots.iter().enumerate().skip(1) {
        let err = l2_rel(&snapshots[0], snapshot);
        println!(
            "behavior_compare reference=FromResolutionAndLatticeVelocity other_index={} l2rel={:.12e} near_lid_ux={:.12e}",
            idx, err, snapshot.near_lid_ux
        );
        assert!(
            err < FIELD_L2_REL_TOL,
            "constructor {idx} field L2rel {err:e} exceeds {FIELD_L2_REL_TOL:e}"
        );
    }
}
