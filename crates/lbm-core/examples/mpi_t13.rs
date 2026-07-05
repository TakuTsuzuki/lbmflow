//! T13-MPI: distributed-vs-monolithic equivalence verification
//! (docs/MPI_GUIDE.md; driven by scripts/test_mpi.sh).
//!
//! Every case runs the same scenario twice, in the same process group:
//!
//! - **distributed**: one part per rank ([`MpiSolver`] over `MpiExchange`),
//! - **baseline**: the monolithic single-part solver on rank 0
//!   (`LocalPeriodic`, the V1-equivalent configuration).
//!
//! At the T13 checkpoints the rank-0 gathered fields (rho, u, and every
//! deviation-population plane) must match the baseline to ≤ 1e-12 in f64
//! (bit-exact is expected and reported), and the Allreduce'd diagnostics
//! (total mass / momentum / probed force / NaN count) to ≤ 1e-11
//! (rank partial sums reassociate).
//!
//! Case sets by world size (optional `argv[1]` substring-filters cases):
//!
//! - `-n 1|2|4`: 2D TGV (periodic), lid-driven cavity (lid crossing the
//!   seam), cylinder + force probe on the seam, and a single-component
//!   Shan–Chen droplet (ψ halos over `exchange_scalar`; at `-n 4` the
//!   droplet sits exactly on the 2×2 corner).
//! - `-n 8`: 3D TGV on a 2×2×2 decomposition (D3Q19).

use lbm_core::dist::MpiSolver;
use lbm_core::lattice::{D2Q9, D3Q19};
use lbm_core::prelude::*;
use mpi::topology::SimpleCommunicator;
use mpi::traits::*;

const FIELD_TOL: f64 = 1e-12;
const DIAG_TOL: f64 = 1e-11;

type Init = Box<dyn Fn(usize, usize, usize) -> (f64, [f64; 3])>;
type Pred = Box<dyn Fn(usize, usize, usize) -> bool>;

/// One scenario, decomposition-agnostic.
struct Setup {
    name: &'static str,
    spec: GlobalSpec<f64>,
    walls: WallSpec<f64>,
    init: Option<Init>,
    profile: Option<(Face, Vec<[f64; 3]>)>,
    /// Solid + probe predicate (probe covers the same cells).
    obstacle: Option<Pred>,
    /// Shan–Chen cohesion strength (classic ψ), refreshed before every step.
    shan_chen: Option<f64>,
    steps: usize,
}

impl Setup {
    fn plain(name: &'static str, spec: GlobalSpec<f64>, steps: usize) -> Self {
        Self {
            name,
            spec,
            walls: WallSpec::default(),
            init: None,
            profile: None,
            obstacle: None,
            shan_chen: None,
            steps,
        }
    }
}

fn psi_classic(rho: f64) -> f64 {
    1.0 - (-rho).exp()
}

/// Max |a - b| over two equally sized fields.
fn max_absdev(a: &[f64], b: &[f64]) -> f64 {
    a.iter()
        .zip(b)
        .map(|(x, y)| (x - y).abs())
        .fold(0.0f64, f64::max)
}

/// Tracks the worst deviations and the pass/fail state of one case.
#[derive(Default)]
struct Score {
    max_field: f64,
    max_diag: f64,
    fail: bool,
}

impl Score {
    fn field(&mut self, name: &str, t: usize, dev: f64) {
        self.max_field = self.max_field.max(dev);
        if dev > FIELD_TOL {
            self.fail = true;
            eprintln!("  FAIL field {name} t={t}: |Δ| = {dev:e} > {FIELD_TOL:e}");
        }
    }

    fn diag(&mut self, name: &str, t: usize, got: f64, want: f64) {
        let dev = (got - want).abs();
        self.max_diag = self.max_diag.max(dev);
        if dev > DIAG_TOL + DIAG_TOL * want.abs() {
            self.fail = true;
            eprintln!("  FAIL diag {name} t={t}: {got:e} vs {want:e} (|Δ| = {dev:e})");
        }
    }
}

/// Run one case: distributed on all ranks, baseline on rank 0, lockstep.
/// Returns the case's pass flag (Allreduce'd, identical on every rank).
fn run_case<L: Lattice>(world: &SimpleCommunicator, setup: &Setup, decomp: [usize; 3]) -> bool {
    let rank = world.rank() as usize;
    let (solid, wall_u) = build_wall_rims(L::D, setup.spec.dims, &setup.walls);
    let mut dist: MpiSolver<L, f64, CpuScalar> = MpiSolver::new(
        world,
        &setup.spec,
        &solid,
        &wall_u,
        decomp,
        CpuScalar::default(),
    );
    let mut base: Option<Solver<L, f64, CpuScalar, LocalPeriodic>> = (rank == 0).then(|| {
        Solver::new(
            &setup.spec,
            &solid,
            &wall_u,
            [1, 1, 1],
            CpuScalar::default(),
            LocalPeriodic,
        )
    });
    let dims = setup.spec.dims;
    if let Some(pred) = &setup.obstacle {
        for z in 0..dims[2] {
            for y in 0..dims[1] {
                for x in 0..dims[0] {
                    if pred(x, y, z) {
                        dist.set_solid(x, y, z);
                        if let Some(b) = &mut base {
                            b.set_solid(x, y, z);
                        }
                    }
                }
            }
        }
        dist.set_force_probe(|x, y, z| pred(x, y, z));
        if let Some(b) = &mut base {
            b.set_force_probe(|x, y, z| pred(x, y, z));
        }
    }
    if let Some(init) = &setup.init {
        dist.init_with(init);
        if let Some(b) = &mut base {
            b.init_with(init);
        }
    }
    if let Some((face, values)) = &setup.profile {
        dist.set_inlet_profile(*face, values);
        if let Some(b) = &mut base {
            b.set_inlet_profile(*face, values);
        }
    }

    let steps = setup.steps;
    let mut checkpoints = vec![1, 2, 3, 5, steps / 2, steps];
    checkpoints.retain(|&c| c > 0);
    checkpoints.dedup();
    let mut score = Score::default();
    let mut t = 0usize;
    for &cp in &checkpoints {
        while t < cp {
            if let Some(g) = setup.shan_chen {
                dist.update_shan_chen_force(g, psi_classic);
                if let Some(b) = &mut base {
                    b.update_shan_chen_force(g, psi_classic);
                }
            }
            dist.step();
            if let Some(b) = &mut base {
                b.step();
            }
            t += 1;
            // Probed force: compare every step (Allreduce'd partial sums).
            if setup.obstacle.is_some() {
                let got = dist.probed_force();
                if let Some(b) = &base {
                    let want = b.probed_force();
                    for a in 0..L::D {
                        score.diag(&format!("probed_force[{a}]"), t, got[a], want[a]);
                    }
                }
            }
        }
        // Distributed NaN check + diagnostics (collective; every rank).
        let nonfinite = dist.nonfinite_count();
        let mass = dist.total_mass();
        let mom = dist.total_momentum();
        // Rank-0 gathers (collective).
        let rho = dist.gather_rho();
        let ux = dist.gather_ux();
        let uy = dist.gather_uy();
        let uz = dist.gather_uz();
        let fs: Vec<Option<Vec<f64>>> = (0..L::Q).map(|q| dist.gather_f(q)).collect();
        if let Some(b) = &base {
            if nonfinite != 0 {
                score.fail = true;
                eprintln!("  FAIL t={t}: {nonfinite} non-finite values");
            }
            score.diag("total_mass", t, mass, b.total_mass());
            let bm = b.total_momentum();
            for a in 0..L::D {
                score.diag(&format!("momentum[{a}]"), t, mom[a], bm[a]);
            }
            score.field("rho", t, max_absdev(rho.as_ref().unwrap(), &b.gather_rho()));
            score.field("ux", t, max_absdev(ux.as_ref().unwrap(), &b.gather_ux()));
            score.field("uy", t, max_absdev(uy.as_ref().unwrap(), &b.gather_uy()));
            if L::D == 3 {
                score.field("uz", t, max_absdev(uz.as_ref().unwrap(), &b.gather_uz()));
            }
            for (q, fq) in fs.iter().enumerate() {
                score.field(
                    &format!("f[{q}]"),
                    t,
                    max_absdev(fq.as_ref().unwrap(), &b.gather_f(q)),
                );
            }
        }
    }
    // Share rank 0's verdict so every rank reports the same exit code.
    let local = i32::from(score.fail);
    let mut global = 0i32;
    world.all_reduce_into(&local, &mut global, mpi::collective::SystemOperation::max());
    if rank == 0 {
        let verdict = if global == 0 { "PASS" } else { "FAIL" };
        println!(
            "{verdict} {name} decomp={decomp:?} steps={steps} \
             (max field |Δ| = {mf:.1e}, max diag |Δ| = {md:.1e})",
            name = setup.name,
            mf = score.max_field,
            md = score.max_diag,
        );
    }
    global == 0
}

// ---------------------------------------------------------------------------
// Cases
// ---------------------------------------------------------------------------

fn tgv2d() -> Setup {
    let (nx, ny) = (96usize, 64usize);
    let mut s = Setup::plain(
        "tgv2d",
        GlobalSpec {
            dims: [nx, ny, 1],
            nu: 0.02,
            periodic: [true, true, false],
            ..Default::default()
        },
        200,
    );
    s.init = Some(Box::new(move |x, y, _| {
        let kx = 2.0 * std::f64::consts::PI / nx as f64;
        let ky = 2.0 * std::f64::consts::PI / ny as f64;
        let u0 = 0.04;
        (
            1.0 + 0.01 * (kx * x as f64).cos() * (ky * y as f64).cos(),
            [
                -u0 * (kx * x as f64).cos() * (ky * y as f64).sin(),
                u0 * (kx * x as f64).sin() * (ky * y as f64).cos(),
                0.0,
            ],
        )
    }));
    s
}

/// Lid-driven cavity: the moving lid (YPos wall) spans the whole top edge,
/// so any x-split seam crosses it.
fn cavity() -> Setup {
    let mut walls = WallSpec::default();
    for f in [Face::XNeg, Face::XPos, Face::YNeg, Face::YPos] {
        walls.is_wall[f.index()] = true;
    }
    walls.u[Face::YPos.index()] = [0.1, 0.0, 0.0];
    let mut s = Setup::plain(
        "cavity",
        GlobalSpec {
            dims: [64, 64, 1],
            nu: 0.02,
            periodic: [false, false, false],
            ..Default::default()
        },
        200,
    );
    s.walls = walls;
    s
}

/// Channel with a cylinder + momentum-exchange probe centred on the seam
/// (at `-n 4` its cells span all four parts) and a parabolic inlet profile
/// distributed across parts.
fn cylinder() -> Setup {
    let (nx, ny) = (64usize, 40usize);
    let mut walls = WallSpec::default();
    walls.is_wall[Face::YNeg.index()] = true;
    walls.is_wall[Face::YPos.index()] = true;
    let mut faces = [FaceBC::Closed; 6];
    faces[Face::XNeg.index()] = FaceBC::Velocity {
        u: [0.05, 0.0, 0.0],
    };
    faces[Face::XPos.index()] = FaceBC::Outflow;
    let mut s = Setup::plain(
        "cylinder+probe",
        GlobalSpec {
            dims: [nx, ny, 1],
            nu: 0.02,
            periodic: [false, false, false],
            faces,
            ..Default::default()
        },
        250,
    );
    s.walls = walls;
    let (cx, cy, r) = (nx as f64 / 2.0, ny as f64 / 2.0 - 0.3, 5.4);
    s.obstacle = Some(Box::new(move |x, y, _| {
        let (dx, dy) = (x as f64 - cx, y as f64 - cy);
        dx * dx + dy * dy < r * r
    }));
    let profile: Vec<[f64; 3]> = (0..ny)
        .map(|y| {
            let yy = y as f64 / (ny - 1) as f64;
            [0.08 * 4.0 * yy * (1.0 - yy), 0.0, 0.0]
        })
        .collect();
    s.profile = Some((Face::XNeg, profile));
    s
}

/// Single-component Shan–Chen droplet centred on the domain midpoint —
/// exactly the 2×2 corner at `-n 4` — exercising `exchange_scalar` ψ halos.
fn droplet() -> Setup {
    let n = 48usize;
    let mut s = Setup::plain(
        "shan-chen-droplet",
        GlobalSpec {
            dims: [n, n, 1],
            nu: 1.0 / 6.0,
            periodic: [true, true, false],
            ..Default::default()
        },
        300,
    );
    s.init = Some(Box::new(move |x, y, _| {
        let (dx, dy) = (
            x as f64 + 0.5 - n as f64 / 2.0,
            y as f64 + 0.5 - n as f64 / 2.0,
        );
        let r = (dx * dx + dy * dy).sqrt();
        // Liquid ≈ 1.90 inside, vapour ≈ 0.16 outside (classic ψ, G = -5).
        let rho = 0.16 + (1.90 - 0.16) * 0.5 * (1.0 - ((r - 10.0) / 2.0).tanh());
        (rho, [0.0, 0.0, 0.0])
    }));
    s.shan_chen = Some(-5.0);
    s
}

fn tgv3d() -> Setup {
    let n = 24usize;
    let mut s = Setup::plain(
        "tgv3d",
        GlobalSpec {
            dims: [n, n, n],
            nu: 0.02,
            periodic: [true, true, true],
            ..Default::default()
        },
        100,
    );
    let k = 2.0 * std::f64::consts::PI / n as f64;
    s.init = Some(Box::new(move |x, y, z| {
        let (xx, yy, zz) = (k * x as f64, k * y as f64, k * z as f64);
        (
            1.0,
            [
                0.03 * xx.sin() * yy.cos() * zz.cos(),
                -0.03 * xx.cos() * yy.sin() * zz.cos(),
                0.0,
            ],
        )
    }));
    s
}

fn main() {
    let universe = mpi::initialize().expect("MPI initialize failed");
    let world = universe.world();
    let size = world.size() as usize;
    let rank = world.rank() as usize;
    let only = std::env::args().nth(1);

    let mut all_pass = true;
    let mut ran = 0usize;
    match size {
        1 | 2 | 4 => {
            let decomp = match size {
                1 => [1, 1, 1],
                2 => [2, 1, 1],
                _ => [2, 2, 1],
            };
            for setup in [tgv2d(), cavity(), cylinder(), droplet()] {
                if only.as_deref().is_some_and(|o| !setup.name.contains(o)) {
                    continue;
                }
                all_pass &= run_case::<D2Q9>(&world, &setup, decomp);
                ran += 1;
            }
        }
        8 => {
            let setup = tgv3d();
            if !only.as_deref().is_some_and(|o| !setup.name.contains(o)) {
                all_pass &= run_case::<D3Q19>(&world, &setup, [2, 2, 2]);
                ran += 1;
            }
        }
        n => {
            if rank == 0 {
                eprintln!("unsupported world size {n}: use -n 1, 2, 4 (2D) or 8 (3D)");
            }
            all_pass = false;
        }
    }
    if rank == 0 {
        println!(
            "mpi_t13 [n={size}]: {}",
            if all_pass {
                format!("ALL PASS ({ran} case(s))")
            } else {
                "FAILURES DETECTED".to_string()
            }
        );
    }
    drop(world);
    drop(universe);
    std::process::exit(i32::from(!all_pass));
}
