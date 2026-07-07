//! V&V master plan lane 5.2: randomized legal-configuration property sweep.
//!
//! The generator is intentionally small and deterministic: failures print the
//! seed, case index, complete `SimConfig`, and obstacle rectangles so the case
//! can be copied into a focused regression test.

use lbm_core::compat::prelude::*;

const DEFAULT_SEED: u64 = 0x5eed_5eed_c0de_5200;
const STEPS: usize = 500;
const MASS_REL_TOL: f64 = 1.0e-11;
const MAX_SPEED_BOUND: f64 = 0.5;

#[derive(Clone, Copy, Debug)]
struct Rect {
    x0: usize,
    y0: usize,
    w: usize,
    h: usize,
}

#[derive(Clone, Debug)]
struct FuzzCase {
    seed: u64,
    index: usize,
    config: SimConfig<f64>,
    obstacles: Vec<Rect>,
}

#[derive(Clone, Copy, Debug)]
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
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        self.state
    }

    fn usize_inclusive(&mut self, lo: usize, hi: usize) -> usize {
        assert!(lo <= hi);
        lo + (self.next_u64() as usize % (hi - lo + 1))
    }

    fn f64(&mut self) -> f64 {
        let mantissa = self.next_u64() >> 11;
        mantissa as f64 * (1.0 / ((1u64 << 53) as f64))
    }

    fn f64_range(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.f64()
    }

    fn sign(&mut self) -> f64 {
        if self.next_u64() & 1 == 0 {
            -1.0
        } else {
            1.0
        }
    }
}

fn configured_seed() -> u64 {
    std::env::var("LBM_PROPERTY_FUZZ_SEED")
        .ok()
        .and_then(|s| {
            u64::from_str_radix(s.trim_start_matches("0x"), 16)
                .or_else(|_| s.parse())
                .ok()
        })
        .unwrap_or(DEFAULT_SEED)
}

fn inlet_velocity(rng: &mut Lcg, ux: f64, uy: f64) -> [f64; 2] {
    let scale = rng.f64_range(0.0, 0.05);
    [ux * scale, uy * scale]
}

fn moving_wall_for(edge: Edge, rng: &mut Lcg) -> EdgeBC<f64> {
    match edge {
        Edge::Left | Edge::Right => EdgeBC::MovingWall {
            u: [0.0, rng.sign() * rng.f64_range(0.0, 0.15)],
        },
        Edge::Bottom | Edge::Top => EdgeBC::MovingWall {
            u: [rng.sign() * rng.f64_range(0.0, 0.15), 0.0],
        },
    }
}

fn wall_edge(edge: Edge, rng: &mut Lcg) -> EdgeBC<f64> {
    if rng.next_u64() % 4 == 0 {
        moving_wall_for(edge, rng)
    } else {
        EdgeBC::BounceBack
    }
}

fn legal_edges(rng: &mut Lcg, case_index: usize) -> Edges<f64> {
    match case_index % 8 {
        // Ensure the default sweep always includes closed all-BB cases for the
        // conservation invariant.
        0 => Edges {
            left: EdgeBC::BounceBack,
            right: EdgeBC::BounceBack,
            bottom: EdgeBC::BounceBack,
            top: EdgeBC::BounceBack,
        },
        1 => Edges {
            left: EdgeBC::Periodic,
            right: EdgeBC::Periodic,
            bottom: EdgeBC::Periodic,
            top: EdgeBC::Periodic,
        },
        2 => Edges {
            left: EdgeBC::Periodic,
            right: EdgeBC::Periodic,
            bottom: wall_edge(Edge::Bottom, rng),
            top: wall_edge(Edge::Top, rng),
        },
        3 => Edges {
            left: wall_edge(Edge::Left, rng),
            right: wall_edge(Edge::Right, rng),
            bottom: EdgeBC::Periodic,
            top: EdgeBC::Periodic,
        },
        4 => Edges {
            left: EdgeBC::VelocityInlet {
                u: inlet_velocity(rng, 1.0, 0.0),
            },
            right: EdgeBC::PressureOutlet {
                rho: rng.f64_range(0.9998, 1.0002),
            },
            bottom: EdgeBC::BounceBack,
            top: EdgeBC::BounceBack,
        },
        5 => Edges {
            left: EdgeBC::PressureOutlet {
                rho: rng.f64_range(0.9998, 1.0002),
            },
            right: EdgeBC::PressureOutlet {
                rho: rng.f64_range(0.9998, 1.0002),
            },
            bottom: EdgeBC::BounceBack,
            top: EdgeBC::BounceBack,
        },
        6 => Edges {
            left: EdgeBC::BounceBack,
            right: EdgeBC::BounceBack,
            bottom: EdgeBC::VelocityInlet {
                u: inlet_velocity(rng, 0.0, 1.0),
            },
            top: EdgeBC::PressureOutlet {
                rho: rng.f64_range(0.9998, 1.0002),
            },
        },
        _ => {
            if rng.next_u64() & 1 == 0 {
                Edges {
                    left: EdgeBC::VelocityInlet {
                        u: inlet_velocity(rng, 1.0, 0.0),
                    },
                    right: EdgeBC::Outflow,
                    bottom: EdgeBC::BounceBack,
                    top: EdgeBC::BounceBack,
                }
            } else {
                Edges {
                    left: EdgeBC::BounceBack,
                    right: EdgeBC::BounceBack,
                    bottom: EdgeBC::PressureOutlet {
                        rho: rng.f64_range(0.9998, 1.0002),
                    },
                    top: EdgeBC::Outflow,
                }
            }
        }
    }
}

fn interior_limits(edges: &Edges<f64>, nx: usize, ny: usize) -> (usize, usize, usize, usize) {
    let x0 = if matches!(
        edges.left,
        EdgeBC::VelocityInlet { .. }
            | EdgeBC::PressureOutlet { .. }
            | EdgeBC::Outflow
            | EdgeBC::ConvectiveOutflow { .. }
    ) {
        2
    } else {
        1
    };
    let x1 = if matches!(
        edges.right,
        EdgeBC::VelocityInlet { .. }
            | EdgeBC::PressureOutlet { .. }
            | EdgeBC::Outflow
            | EdgeBC::ConvectiveOutflow { .. }
    ) {
        nx - 3
    } else {
        nx - 2
    };
    let y0 = if matches!(
        edges.bottom,
        EdgeBC::VelocityInlet { .. }
            | EdgeBC::PressureOutlet { .. }
            | EdgeBC::Outflow
            | EdgeBC::ConvectiveOutflow { .. }
    ) {
        2
    } else {
        1
    };
    let y1 = if matches!(
        edges.top,
        EdgeBC::VelocityInlet { .. }
            | EdgeBC::PressureOutlet { .. }
            | EdgeBC::Outflow
            | EdgeBC::ConvectiveOutflow { .. }
    ) {
        ny - 3
    } else {
        ny - 2
    };
    (x0, x1, y0, y1)
}

fn random_obstacles(rng: &mut Lcg, edges: &Edges<f64>, nx: usize, ny: usize) -> Vec<Rect> {
    let (x0, x1, y0, y1) = interior_limits(edges, nx, ny);
    let obstacle_count = rng.usize_inclusive(0, 3);
    let mut rects = Vec::with_capacity(obstacle_count);
    for _ in 0..obstacle_count {
        let avail_w = x1 - x0 + 1;
        let avail_h = y1 - y0 + 1;
        let w = rng.usize_inclusive(1, avail_w.min(8));
        let h = rng.usize_inclusive(1, avail_h.min(8));
        let rx = rng.usize_inclusive(x0, x1 - w + 1);
        let ry = rng.usize_inclusive(y0, y1 - h + 1);
        rects.push(Rect {
            x0: rx,
            y0: ry,
            w,
            h,
        });
    }
    rects
}

fn generate_case(seed: u64, index: usize, rng: &mut Lcg) -> FuzzCase {
    let nx = rng.usize_inclusive(8, 96);
    let ny = rng.usize_inclusive(8, 96);
    let nu = rng.f64_range(1.0 / 12.0, 0.5);
    let edges = legal_edges(rng, index);
    let force_mag = rng.f64_range(0.0, 1.0e-5);
    let force_angle = rng.f64_range(0.0, std::f64::consts::TAU);
    let config = SimConfig {
        nx,
        ny,
        nu,
        collision: if rng.next_u64() & 1 == 0 {
            Collision::Bgk
        } else {
            Collision::Trt {
                magic: Collision::MAGIC_STD,
            }
        },
        edges,
        force: [force_mag * force_angle.cos(), force_mag * force_angle.sin()],
    };
    let obstacles = random_obstacles(rng, &config.edges, nx, ny);
    FuzzCase {
        seed,
        index,
        config,
        obstacles,
    }
}

fn case_label(case: &FuzzCase) -> String {
    format!(
        "property_fuzz seed=0x{:016x} case={} config={:#?} obstacles={:#?}",
        case.seed, case.index, case.config, case.obstacles
    )
}

fn build_case(case: &FuzzCase) -> Simulation<f64> {
    let label = case_label(case);
    let mut sim = case
        .config
        .clone()
        .build()
        .unwrap_or_else(|err| panic!("{label}\nfailed to build legal config: {err}"));
    for rect in &case.obstacles {
        for y in rect.y0..rect.y0 + rect.h {
            for x in rect.x0..rect.x0 + rect.w {
                assert!(
                    sim.set_solid_allowed(x, y),
                    "{label}\nillegal obstacle cell ({x},{y})"
                );
                sim.set_solid(x, y);
            }
        }
    }
    sim
}

fn is_all_bounce_back(edges: &Edges<f64>) -> bool {
    matches!(edges.left, EdgeBC::BounceBack)
        && matches!(edges.right, EdgeBC::BounceBack)
        && matches!(edges.bottom, EdgeBC::BounceBack)
        && matches!(edges.top, EdgeBC::BounceBack)
}

fn assert_finite_and_bounded(sim: &Simulation<f64>, label: &str) {
    for (i, &rho) in sim.rho_field().iter().enumerate() {
        assert!(
            rho.is_finite(),
            "{label}\nnon-finite rho at compact cell {i}: {rho:?}"
        );
    }
    let mut max_u = 0.0f64;
    for i in 0..sim.ux_field().len() {
        let ux = sim.ux_field()[i];
        let uy = sim.uy_field()[i];
        assert!(
            ux.is_finite(),
            "{label}\nnon-finite ux at compact cell {i}: {ux:?}"
        );
        assert!(
            uy.is_finite(),
            "{label}\nnon-finite uy at compact cell {i}: {uy:?}"
        );
        max_u = max_u.max(ux.hypot(uy));
    }
    assert!(
        max_u <= MAX_SPEED_BOUND,
        "{label}\nmax |u| = {max_u:e}, bound = {MAX_SPEED_BOUND:e}"
    );
}

fn assert_bit_identical_fields(a: &Simulation<f64>, b: &Simulation<f64>, label: &str) {
    for (name, fa, fb) in [
        ("rho", a.rho_field(), b.rho_field()),
        ("ux", a.ux_field(), b.ux_field()),
        ("uy", a.uy_field(), b.uy_field()),
    ] {
        assert_eq!(fa.len(), fb.len(), "{label}\n{name} length mismatch");
        for i in 0..fa.len() {
            assert_eq!(
                fa[i].to_bits(),
                fb[i].to_bits(),
                "{label}\n{name} bit mismatch at compact cell {i}: {:016x} != {:016x} ({} != {})",
                fa[i].to_bits(),
                fb[i].to_bits(),
                fa[i],
                fb[i]
            );
        }
    }
}

fn run_case(case: &FuzzCase) {
    println!(
        "property_fuzz seed=0x{:016x} case={}",
        case.seed, case.index
    );
    let label = case_label(case);
    let mut a = build_case(case);
    let mut b = build_case(case);
    let m0 = a.total_mass_f64();
    a.run(STEPS);
    b.run(STEPS);
    assert_finite_and_bounded(&a, &label);
    assert_bit_identical_fields(&a, &b, &label);
    if is_all_bounce_back(&case.config.edges) {
        let m1 = a.total_mass_f64();
        let drift = ((m1 - m0) / m0).abs();
        assert!(
            drift < MASS_REL_TOL,
            "{label}\nclosed all-BB mass drift = {drift:e}, bound = {MASS_REL_TOL:e}, m0={m0:.17e}, m1={m1:.17e}"
        );
    }
}

fn run_sweep(seed: u64, cases: usize) {
    println!("property_fuzz sweep seed=0x{seed:016x} cases={cases} steps={STEPS}");
    let mut rng = Lcg::new(seed);
    for index in 0..cases {
        let case = generate_case(seed, index, &mut rng);
        run_case(&case);
    }
}

#[test]
fn randomized_legal_config_invariant_sweep_50() {
    run_sweep(configured_seed(), 50);
}

#[test]
#[ignore = "deep randomized legal-config sweep; run with --include-ignored"]
fn randomized_legal_config_invariant_sweep_500() {
    run_sweep(configured_seed() ^ 0xfeed_face_d15e_a5e5, 500);
}
