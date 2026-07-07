use lbm_scenario::{
    build, build_check, presets, CollisionSpec, EdgeSpec, EdgesSpec, Grid, InitSpec, Obstacle,
    Physics, Precision, RunSpec, Scenario, SimHandle,
};

const FUZZ_SEED: u64 = 0x8100_0000_0000_0001;

#[test]
fn g1_presets_roundtrip_identity() {
    let preset_list = presets();
    assert!(
        !preset_list.is_empty(),
        "presets() must expose named scenarios"
    );

    for (name, _, original) in preset_list {
        let json = serde_json::to_string(&original).unwrap();
        let reloaded: Scenario = serde_json::from_str(&json)
            .unwrap_or_else(|e| panic!("G1 preset {name}: reload failed: {e}"));
        let json_again = serde_json::to_string(&reloaded).unwrap();
        let original_value = serde_json::to_value(&original).unwrap();
        let reloaded_value = serde_json::to_value(&reloaded).unwrap();

        println!("G1 preset {name}: round-trip bytes={} ok", json.len());
        assert_eq!(
            original_value, reloaded_value,
            "G1 preset {name}: serialized Scenario value changed"
        );
        assert_eq!(
            json, json_again,
            "G1 preset {name}: serialization is not stable after reload"
        );
    }
}

#[test]
fn g2_malformed_inputs_are_rejected_with_error_class() {
    let cases = [
        MalformedCase {
            name: "missing_required_grid",
            class: "serde::missing_field",
            json: r#"{
                "name": "missing-grid",
                "physics": { "nu": 0.05 },
                "edges": {
                    "left": { "type": "periodic" },
                    "right": { "type": "periodic" },
                    "bottom": { "type": "bounceBack" },
                    "top": { "type": "bounceBack" }
                },
                "run": { "steps": 1 }
            }"#,
            expected: "missing field `grid`",
        },
        MalformedCase {
            name: "wrong_type_grid_nx",
            class: "serde::wrong_type",
            json: r#"{
                "name": "wrong-type",
                "grid": { "nx": "wide", "ny": 16 },
                "physics": { "nu": 0.05 },
                "edges": {
                    "left": { "type": "periodic" },
                    "right": { "type": "periodic" },
                    "bottom": { "type": "bounceBack" },
                    "top": { "type": "bounceBack" }
                },
                "run": { "steps": 1 }
            }"#,
            expected: "invalid type",
        },
        MalformedCase {
            name: "out_of_range_nu",
            class: "build_check::NonPositiveViscosity",
            json: r#"{
                "name": "bad-nu",
                "grid": { "nx": 16, "ny": 16 },
                "physics": { "nu": 0.0 },
                "edges": {
                    "left": { "type": "periodic" },
                    "right": { "type": "periodic" },
                    "bottom": { "type": "bounceBack" },
                    "top": { "type": "bounceBack" }
                },
                "run": { "steps": 1 }
            }"#,
            expected: "kinematic viscosity must be > 0",
        },
        MalformedCase {
            name: "out_of_range_velocity",
            class: "build_check::VelocityTooHigh",
            json: r#"{
                "name": "too-fast-wall",
                "grid": { "nx": 16, "ny": 16 },
                "physics": { "nu": 0.05 },
                "edges": {
                    "left": { "type": "bounceBack" },
                    "right": { "type": "bounceBack" },
                    "bottom": { "type": "bounceBack" },
                    "top": { "type": "movingWall", "u": [0.31, 0.0] }
                },
                "run": { "steps": 1 }
            }"#,
            expected: "low-Mach limit",
        },
        MalformedCase {
            name: "illegal_orthogonal_open_faces",
            class: "build_check::AdjacentOpenEdges",
            json: r#"{
                "name": "orthogonal-open-faces",
                "grid": { "nx": 16, "ny": 16 },
                "physics": { "nu": 0.05 },
                "edges": {
                    "left": { "type": "velocityInlet", "u": [0.05, 0.0] },
                    "right": { "type": "pressureOutlet", "rho": 1.0 },
                    "bottom": { "type": "outflow" },
                    "top": { "type": "bounceBack" }
                },
                "run": { "steps": 1 }
            }"#,
            expected: "open edges",
        },
        MalformedCase {
            name: "invalid_inlet_profile_kind",
            class: "serde::unknown_variant",
            json: r#"{
                "name": "bad-profile-kind",
                "grid": { "nx": 16, "ny": 16 },
                "physics": { "nu": 0.05 },
                "edges": {
                    "left": { "type": "velocityInlet", "u": [0.05, 0.0] },
                    "right": { "type": "pressureOutlet", "rho": 1.0 },
                    "bottom": { "type": "bounceBack" },
                    "top": { "type": "bounceBack" }
                },
                "inletProfile": { "edge": "left", "kind": "plug", "umax": 0.05 },
                "run": { "steps": 1 }
            }"#,
            expected: "unknown variant",
        },
        MalformedCase {
            name: "invalid_inlet_profile_edge",
            class: "build_check::InvalidParameter",
            json: r#"{
                "name": "profile-not-on-inlet",
                "grid": { "nx": 16, "ny": 16 },
                "physics": { "nu": 0.05 },
                "edges": {
                    "left": { "type": "periodic" },
                    "right": { "type": "periodic" },
                    "bottom": { "type": "bounceBack" },
                    "top": { "type": "bounceBack" }
                },
                "inletProfile": { "edge": "top", "kind": "parabolic", "umax": 0.05 },
                "run": { "steps": 1 }
            }"#,
            expected: "inletProfile edge",
        },
        MalformedCase {
            name: "unknown_obstacle_shape",
            class: "serde::unknown_variant",
            json: r#"{
                "name": "unknown-obstacle",
                "grid": { "nx": 16, "ny": 16 },
                "physics": { "nu": 0.05 },
                "edges": {
                    "left": { "type": "periodic" },
                    "right": { "type": "periodic" },
                    "bottom": { "type": "bounceBack" },
                    "top": { "type": "bounceBack" }
                },
                "obstacles": [
                    { "shape": "triangle", "x0": 4, "y0": 4, "x1": 8, "y1": 8 }
                ],
                "run": { "steps": 1 }
            }"#,
            expected: "unknown variant",
        },
    ];

    for case in cases {
        let err = load_and_build_check(case.json).unwrap_err();
        println!("G2 {}: {}: {}", case.name, case.class, err);
        assert!(
            !err.trim().is_empty(),
            "G2 {}: error message must be non-empty",
            case.name
        );
        assert!(
            err.contains(case.expected),
            "G2 {}: expected error containing {:?}, got {err:?}",
            case.name,
            case.expected
        );
    }
}

#[test]
fn g3_seeded_legal_configs_build_and_step_without_nan() {
    let mut rng = Lcg::new(FUZZ_SEED);
    println!("G3 LCG seed: 0x{FUZZ_SEED:016x}");

    for case_index in 0..20 {
        let sc = legal_case(case_index, &mut rng);
        println!(
            "G3 seed=0x{FUZZ_SEED:016x} case={case_index} name={} nx={} ny={} nu={:.12}",
            sc.name, sc.grid.nx, sc.grid.ny, sc.physics.nu
        );
        let mut sim = match build(&sc) {
            Ok(SimHandle::F64(sim, None)) => sim,
            Ok(_) => panic!("G3 case {case_index}: expected f64 single-phase 2D build"),
            Err(e) => panic!("G3 case {case_index}: build failed: {e}"),
        };

        for step in 0..10 {
            sim.step();
            assert_finite_field(case_index, step + 1, "rho", sim.rho_field());
            assert_finite_field(case_index, step + 1, "ux", sim.ux_field());
            assert_finite_field(case_index, step + 1, "uy", sim.uy_field());
        }
    }
}

struct MalformedCase {
    name: &'static str,
    class: &'static str,
    json: &'static str,
    expected: &'static str,
}

fn load_and_build_check(json: &str) -> Result<Scenario, String> {
    let sc: Scenario = serde_json::from_str(json).map_err(|e| e.to_string())?;
    build_check(&sc).map_err(|e| e.to_string())?;
    Ok(sc)
}

fn legal_case(case_index: usize, rng: &mut Lcg) -> Scenario {
    let nx = rng.usize_inclusive(4, 128);
    let ny = rng.usize_inclusive(4, 128);
    let nu = rng.f64_between(1.0 / 12.0, 1.0 / 6.0);
    let edge_mode = rng.usize_inclusive(0, 3);
    let edges = match edge_mode {
        0 => EdgesSpec {
            left: EdgeSpec::Periodic,
            right: EdgeSpec::Periodic,
            bottom: EdgeSpec::Periodic,
            top: EdgeSpec::Periodic,
            front: None,
            back: None,
        },
        1 => EdgesSpec {
            left: EdgeSpec::Periodic,
            right: EdgeSpec::Periodic,
            bottom: EdgeSpec::BounceBack,
            top: EdgeSpec::MovingWall {
                u: [rng.f64_between(0.0, 0.1), 0.0],
            },
            front: None,
            back: None,
        },
        2 => EdgesSpec {
            left: EdgeSpec::VelocityInlet {
                u: [rng.f64_between(0.005, 0.08), 0.0],
            },
            right: if rng.bool() {
                EdgeSpec::PressureOutlet { rho: 1.0 }
            } else {
                EdgeSpec::Outflow
            },
            bottom: EdgeSpec::BounceBack,
            top: EdgeSpec::BounceBack,
            front: None,
            back: None,
        },
        _ => EdgesSpec {
            left: EdgeSpec::BounceBack,
            right: EdgeSpec::BounceBack,
            bottom: EdgeSpec::VelocityInlet {
                u: [0.0, rng.f64_between(0.005, 0.08)],
            },
            top: EdgeSpec::PressureOutlet { rho: 1.0 },
            front: None,
            back: None,
        },
    };

    Scenario {
        version: 0,
        name: format!("g3-legal-{case_index:02}"),
        grid: Grid {
            nx,
            ny,
            nz: 1,
            lattice: None,
        },
        physics: Physics {
            nu,
            collision: if rng.bool() {
                CollisionSpec::Trt
            } else {
                CollisionSpec::Bgk
            },
            force: [0.0, 0.0],
            gravity: None,
            precision: Precision::F64,
        },
        units: None,
        compute: None,
        wall: None,
        edges,
        inlet_profile: None,
        obstacles: obstacles(nx, ny, rng),
        init: InitSpec::Rest,
        multiphase: None,
        rotor: None,
        particles: None,
        run: RunSpec {
            steps: 10,
            stop_when_steady: None,
        },
        probes: Vec::new(),
        outputs: Vec::new(),
    }
}

fn obstacles(nx: usize, ny: usize, rng: &mut Lcg) -> Vec<Obstacle> {
    let mut obs = Vec::with_capacity(5);
    let x_min = 1usize;
    let x_max = nx.saturating_sub(2).max(1);
    let y_min = 1usize;
    let y_max = ny.saturating_sub(2).max(1);

    for _ in 0..5 {
        let x = rng.usize_inclusive(x_min, x_max);
        let y = rng.usize_inclusive(y_min, y_max);
        if rng.bool() {
            obs.push(Obstacle::Rect {
                x0: x,
                y0: y,
                x1: x,
                y1: y,
            });
        } else {
            obs.push(Obstacle::Circle {
                cx: x as f64,
                cy: y as f64,
                r: 0.5,
            });
        }
    }

    obs
}

fn assert_finite_field(case_index: usize, step: usize, field: &str, values: &[f64]) {
    for (idx, value) in values.iter().enumerate() {
        assert!(
            value.is_finite(),
            "G3 case {case_index} step {step}: {field}[{idx}] is {value}"
        );
    }
}

struct Lcg {
    state: u64,
}

impl Lcg {
    fn new(seed: u64) -> Self {
        Self { state: seed }
    }

    fn next_u64(&mut self) -> u64 {
        self.state = self
            .state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        self.state
    }

    fn bool(&mut self) -> bool {
        self.next_u64() & 1 == 0
    }

    fn usize_inclusive(&mut self, min: usize, max: usize) -> usize {
        assert!(min <= max);
        let width = max - min + 1;
        min + (self.next_u64() as usize % width)
    }

    fn f64_between(&mut self, min: f64, max: f64) -> f64 {
        assert!(min <= max);
        let unit = ((self.next_u64() >> 11) as f64) * (1.0 / ((1u64 << 53) as f64));
        min + unit * (max - min)
    }
}
