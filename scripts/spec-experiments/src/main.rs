//! Small validation experiments for docs/SOLVER_IMPROVEMENT_SPEC.md (v0 -> v1).
//! Each subcommand prints machine-readable RESULT lines.

use lbm_core::lattice::{D2Q9, D3Q19};
use lbm_core::prelude::*;

fn main() {
    let which = std::env::args().nth(1).unwrap_or_default();
    match which.as_str() {
        "e2" => e2_uncovered_face(),
        "e3" => e3_nu_zero(),
        "e4" => e4_halo_scope(),
        "e5" => e5_outflow_plug(),
        "e6" => e6_nan_config(),
        "e7" => e7_normal_moving_wall(),
        "e8" => e8_two_pass_probe(),
        _ => eprintln!("usage: expt e2|e3|e4|e5|e6|e7|e8"),
    }
}

// ---------------------------------------------------------------- E2
// V2 native, D3Q19: z faces neither periodic nor open nor walled
// ("uncovered"). Claim: stale slots freeze at initial values and corrupt
// the solution silently (finite everywhere).
fn e2_uncovered_face() {
    let dims = [32usize, 16, 8];
    let n = dims[0] * dims[1] * dims[2];
    let build = |periodic_z: bool| {
        let spec = GlobalSpec::<f64> {
            dims,
            nu: 1.0 / 6.0,
            periodic: [true, true, periodic_z],
            ..Default::default() // faces: all Closed
        };
        let solid = vec![false; n];
        let wall_u = vec![[0.0f64; 3]; n];
        let mut s: Solver<D3Q19, f64, CpuScalar, LocalPeriodic> = Solver::new(
            &spec,
            &solid,
            &wall_u,
            [1, 1, 1],
            CpuScalar::default(),
            LocalPeriodic,
        );
        // z-uniform initial condition: solution must stay z-invariant.
        s.init_with(|x, _y, _z| {
            let k = 2.0 * std::f64::consts::PI * x as f64 / dims[0] as f64;
            (1.0 + 0.05 * k.sin(), [0.0; 3])
        });
        s
    };

    for (name, pz) in [("covered(periodic z)", true), ("UNCOVERED(z faces bare)", false)] {
        let mut s = build(pz);
        let m0 = s.total_mass();
        s.run(100);
        // z-invariance violation: max |rho(x,y,0) - rho(x,y,nz/2)|
        let rho = s.gather_rho();
        let uz = s.gather_uz();
        let plane = dims[0] * dims[1];
        let mut zdev = 0.0f64;
        for c in 0..plane {
            let d = (rho[c] - rho[c + plane * (dims[2] / 2)]).abs();
            if d > zdev {
                zdev = d;
            }
        }
        let uzmax = uz.iter().fold(0.0f64, |a, v| a.max(v.abs()));
        let drift = (s.total_mass() - m0).abs();
        println!(
            "RESULT e2 {name}: nonfinite={} mass_drift={drift:.3e} z_invariance_dev={zdev:.3e} max|uz|={uzmax:.3e}",
            s.local_nonfinite_count()
        );
    }
}

// ---------------------------------------------------------------- E3
// nu = 0 passes with no error; TRT omega_m collapses to 0.
fn e3_nu_zero() {
    let (op_t, om_t) = CollisionKind::Trt {
        magic: CollisionKind::MAGIC_STD,
    }
    .omegas(0.0);
    let (op_b, om_b) = CollisionKind::Bgk.omegas(0.0);
    println!("RESULT e3 trt: omega_p={op_t} omega_m={om_t}");
    println!("RESULT e3 bgk: omega_p={op_b} omega_m={om_b}");
    // And the solver constructs + steps without any error:
    let spec = GlobalSpec::<f64> {
        nu: 0.0,
        ..Default::default()
    };
    let n = 64 * 64;
    let mut s: Solver<D2Q9, f64, CpuScalar, LocalPeriodic> = Solver::new(
        &spec,
        &vec![false; n],
        &vec![[0.0; 3]; n],
        [1, 1, 1],
        CpuScalar::default(),
        LocalPeriodic,
    );
    s.init_with(|x, _, _| (1.0 + 0.01 * (x as f64 * 0.3).sin(), [0.0; 3]));
    s.run(10);
    println!(
        "RESULT e3 solver: built+stepped 10 with nu=0, nonfinite={} (no error raised)",
        s.local_nonfinite_count()
    );
}

// ---------------------------------------------------------------- E4
// new_local_part(part=1 of [2,1,1]) + LocalPeriodic: neighbor id 0 resolves
// to the local (only) part -> silent self-wrap instead of the real neighbor.
fn e4_halo_scope() {
    let dims = [64usize, 32, 1];
    let n = dims[0] * dims[1];
    let spec = GlobalSpec::<f64> {
        dims,
        nu: 1.0 / 6.0,
        periodic: [true, false, false],
        ..Default::default()
    };
    // y walls via rims so the y axis needs no exchange.
    let walls = WallSpec::<f64> {
        is_wall: [false, false, true, true, false, false],
        u: [[0.0; 3]; 6],
    };
    let (solid, wall_u) = build_wall_rims(2, dims, &walls);
    let init = |x: usize, _y: usize, _z: usize| {
        let k = 2.0 * std::f64::consts::PI * x as f64 / dims[0] as f64;
        (1.0 + 0.05 * k.sin(), [0.0; 3])
    };

    // Ground truth: correct 2-part in-process run.
    let mut good: Solver<D2Q9, f64, CpuScalar, InProcess> = Solver::new(
        &spec,
        &solid,
        &wall_u,
        [2, 1, 1],
        CpuScalar::default(),
        InProcess,
    );
    good.init_with(init);

    // Misuse: own only part 1, but wire a local exchange.
    let mut bad: Solver<D2Q9, f64, CpuScalar, LocalPeriodic> = Solver::new_local_part(
        &spec,
        &solid,
        &wall_u,
        [2, 1, 1],
        1,
        CpuScalar::default(),
        LocalPeriodic,
    );
    bad.init_with(init);

    good.run(50);
    bad.run(50);

    // Sample rho along the owned block (x in [32,64)) mid-height.
    let mut maxd = 0.0f64;
    for x in 32..64 {
        let d = (good.rho(x, 16, 0) - bad.rho(x, 16, 0)).abs();
        if d > maxd {
            maxd = d;
        }
    }
    println!(
        "RESULT e4: ran without panic; max|rho_good-rho_bad| over owned block = {maxd:.3e} \
         (0 would mean harmless; >>eps means silent wrong physics)"
    );
}

// ---------------------------------------------------------------- E5 (V1)
// Outflow edge cell whose inward neighbour is solid: BC skipped, unknown
// slots frozen -> silent steady flow into the wall + mass drift.
fn e5_outflow_plug() {
    use lbm_core::compat::prelude::*;
    let run = |plug_x_off: usize| -> (f64, f64, f64, f64) {
        let (nx, ny) = (32usize, 16usize);
        let mut sim: Simulation<f64> = SimConfig {
            nx,
            ny,
            nu: 1.0 / 6.0,
            edges: Edges {
                left: EdgeBC::VelocityInlet { u: [0.05, 0.0] },
                right: EdgeBC::Outflow,
                bottom: EdgeBC::BounceBack,
                top: EdgeBC::BounceBack,
            },
            ..Default::default()
        }
        .build()
        .unwrap();
        // Partial plug: solid column at x = nx-1-plug_x_off, lower half only.
        let px = nx - 1 - plug_x_off;
        for y in 1..=7 {
            sim.set_solid(px, y);
        }
        let m0 = sim.total_mass();
        sim.run(2000);
        let pocket_ux = sim.ux(nx - 1, 4); // edge cell behind the plug (pocket)
        let pocket_uy = sim.uy(nx - 1, 4);
        (m0, sim.total_mass(), pocket_ux, pocket_uy)
    };
    let (m0a, m1a, uxa, uya) = run(1); // plug at nx-2: inward neighbour of edge cell -> BUG PATH
    let (m0b, m1b, uxb, uyb) = run(2); // plug at nx-3: legal control
    println!(
        "RESULT e5 plug@nx-2 (bug path): mass {m0a:.3} -> {m1a:.3} (drift {:+.3}), pocket edge cell u=({uxa:+.4},{uya:+.4})",
        m1a - m0a
    );
    println!(
        "RESULT e5 plug@nx-3 (control):  mass {m0b:.3} -> {m1b:.3} (drift {:+.3}), pocket edge cell u=({uxb:+.4},{uyb:+.4})",
        m1b - m0b
    );

    // E5b: direct freeze proof. Seed the pocket with rho=2, watch whether the
    // edge cell relaxes to ambient (legal) or stays pinned by frozen unknown
    // slots (bug). No inlet: quiescent box with an Outflow on the right, so
    // relaxation to rho=1 is the only physics.
    use lbm_core::compat::prelude::*;
    let runb = |plug_x_off: usize| -> (f64, f64, f64) {
        let (nx, ny) = (32usize, 16usize);
        let mut sim: Simulation<f64> = SimConfig {
            nx,
            ny,
            nu: 1.0 / 6.0,
            edges: Edges {
                left: EdgeBC::BounceBack,
                right: EdgeBC::Outflow,
                bottom: EdgeBC::BounceBack,
                top: EdgeBC::BounceBack,
            },
            ..Default::default()
        }
        .build()
        .unwrap();
        let px = nx - 1 - plug_x_off;
        for y in 1..=7 {
            sim.set_solid(px, y);
        }
        sim.init_with(|x, y| {
            if x > px && (1..=7).contains(&y) {
                (2.0, 0.0, 0.0)
            } else {
                (1.0, 0.0, 0.0)
            }
        });
        let r0 = sim.rho(nx - 1, 4);
        sim.run(2000);
        (r0, sim.rho(nx - 1, 4), sim.ux(nx - 1, 4))
    };
    let (r0a, r1a, uxa2) = runb(1);
    let (r0b, r1b, uxb2) = runb(2);
    println!(
        "RESULT e5b plug@nx-2 (bug path): pocket rho {r0a:.3} -> {r1a:.4} after 2000 steps (ux={uxa2:+.5}); frozen if stays >>1"
    );
    println!(
        "RESULT e5b plug@nx-3 (control):  pocket rho {r0b:.3} -> {r1b:.4} after 2000 steps (ux={uxb2:+.5}); should relax to ~1"
    );
}

// ---------------------------------------------------------------- E6 (V1)
// NaN velocities pass SimConfig::validate.
fn e6_nan_config() {
    use lbm_core::compat::prelude::*;
    // (a) NaN inlet
    let r = SimConfig::<f64> {
        nx: 32,
        ny: 16,
        nu: 1.0 / 6.0,
        edges: Edges {
            left: EdgeBC::VelocityInlet {
                u: [f64::NAN, 0.0],
            },
            right: EdgeBC::Outflow,
            bottom: EdgeBC::BounceBack,
            top: EdgeBC::BounceBack,
        },
        ..Default::default()
    }
    .build();
    match r {
        Ok(mut sim) => {
            sim.run(3);
            let bad = sim.rho_field().iter().filter(|v| !v.is_finite()).count();
            println!("RESULT e6 NaN-inlet: build()=Ok (BUG), nonfinite rho cells after 3 steps = {bad}");
        }
        Err(e) => println!("RESULT e6 NaN-inlet: build()=Err({e:?}) (already fixed?)"),
    }
    // (b) NaN moving wall == silent stationary wall
    let mk = |top: EdgeBC<f64>| -> Simulation<f64> {
        let mut sim = SimConfig {
            nx: 32,
            ny: 32,
            nu: 1.0 / 6.0,
            edges: Edges {
                left: EdgeBC::BounceBack,
                right: EdgeBC::BounceBack,
                bottom: EdgeBC::BounceBack,
                top,
            },
            ..Default::default()
        }
        .build()
        .expect("build");
        sim.init_with(|x, y| (1.0, 0.01 * ((x + y) as f64 * 0.1).sin(), 0.0));
        sim.run(100);
        sim
    };
    let a = mk(EdgeBC::MovingWall {
        u: [f64::NAN, 0.0],
    });
    let b = mk(EdgeBC::BounceBack);
    let ident = a
        .rho_field()
        .iter()
        .zip(b.rho_field())
        .all(|(x, y)| x.to_bits() == y.to_bits())
        && a.ux_field()
            .iter()
            .zip(b.ux_field())
            .all(|(x, y)| x.to_bits() == y.to_bits());
    println!(
        "RESULT e6 NaN-MovingWall: build()=Ok (BUG), fields bit-identical to stationary wall = {ident} \
         (true = NaN wall silently degrades to BounceBack)"
    );
}

// ---------------------------------------------------------------- E7 (V1)
// Wall-normal MovingWall component: silent mass injection in a closed box.
fn e7_normal_moving_wall() {
    use lbm_core::compat::prelude::*;
    let run = |u: [f64; 2]| -> (f64, f64) {
        let mut sim: Simulation<f64> = SimConfig {
            nx: 32,
            ny: 32,
            nu: 1.0 / 6.0,
            edges: Edges {
                left: EdgeBC::BounceBack,
                right: EdgeBC::BounceBack,
                bottom: EdgeBC::MovingWall { u },
                top: EdgeBC::BounceBack,
            },
            ..Default::default()
        }
        .build()
        .unwrap();
        let m0 = sim.total_mass();
        sim.run(500);
        (m0, sim.total_mass())
    };
    let (t0, t1) = run([0.05, 0.0]); // tangential: legal
    let (n0, n1) = run([0.0, -0.05]); // normal: silently non-conservative
    println!(
        "RESULT e7 tangential u=[0.05,0]:  mass {t0:.6} -> {t1:.6} (drift {:+.3e})",
        t1 - t0
    );
    println!(
        "RESULT e7 normal     u=[0,-0.05]: mass {n0:.6} -> {n1:.6} (drift {:+.3e}, {:+.1}%)",
        n1 - n0,
        100.0 * (n1 - n0) / n0
    );
}

// ---------------------------------------------------------------- E8
// two_pass on a width-1 axis: boundary shells overlap -> probe force
// double-counted (fields stay correct).
fn e8_two_pass_probe() {
    let dims = [64usize, 1, 1];
    let n = dims[0];
    let run = |two_pass: bool| -> ([f64; 3], f64) {
        let spec = GlobalSpec::<f64> {
            dims,
            nu: 1.0 / 6.0,
            periodic: [true, true, false],
            force: [1e-5, 0.0, 0.0],
            ..Default::default()
        };
        let mut solid = vec![false; n];
        solid[32] = true; // obstacle in the 1-cell channel
        let wall_u = vec![[0.0f64; 3]; n];
        let mut s: Solver<D2Q9, f64, CpuScalar, LocalPeriodic> = Solver::new(
            &spec,
            &solid,
            &wall_u,
            [1, 1, 1],
            CpuScalar::default(),
            LocalPeriodic,
        );
        s.set_two_pass(two_pass);
        s.set_force_probe(|x, _, _| x == 32);
        s.run(20);
        (s.probed_force(), s.total_mass())
    };
    let (f_off, m_off) = run(false);
    let (f_on, m_on) = run(true);
    let ratio = if f_off[0] != 0.0 { f_on[0] / f_off[0] } else { f64::NAN };
    println!(
        "RESULT e8: probed_force x: off={:.6e} on={:.6e} ratio={ratio:.3} (2.0 = double count); \
         mass off={m_off:.9} on={m_on:.9} (equal = fields unharmed)",
        f_off[0], f_on[0]
    );
}
