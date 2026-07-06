//! MPI benchmark driver for the cluster campaign checklist in docs/MPI_GUIDE.md.
//!
//! Backward-compatible default:
//! `mpirun -n <ranks> bench_mpi [local_edge] [steps]`
//! still runs the original 2D weak-scaling case.
//!
//! Cluster modes are dormant behind `--mode`:
//! - `weak2d`: D2Q9, fixed 2D block per rank.
//! - `weak3d`: D3Q19, fixed 3D block per rank.
//! - `strong3d`: D3Q19, fixed global 3D grid.
//! - `diagnostics`: D3Q19 workload plus timed Allreduce/gather diagnostics.
//! - `placement`: MPI thread level, environment, host, and affinity evidence.
//! - `correctness`: D3Q19 smoke that reports nonfinite/mass/momentum after steps.
//!
//! External MPI profiling remains external by design: run the same workload
//! under mpiP, MPI_T tooling, or launcher verbosity to capture MPI-stack ratios
//! and BTL/MTL thresholds without perturbing normal benchmark runs.

use lbm_core::dist::{choose_decomp, MpiSolver};
use lbm_core::lattice::{Lattice, D2Q9, D3Q19};
use lbm_core::prelude::*;
use mpi::traits::*;
use std::process::Command;
use std::time::Instant;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Mode {
    Weak2d,
    Weak3d,
    Strong3d,
    Diagnostics,
    Placement,
    Correctness,
}

#[derive(Clone, Debug)]
struct Config {
    mode: Mode,
    local_edge: usize,
    global_edge: usize,
    steps: usize,
    warmup: usize,
    diagnostics_every: usize,
    gather_rho: bool,
    parallel: bool,
    decomp: Option<[usize; 3]>,
}

impl Default for Config {
    fn default() -> Self {
        Self {
            mode: Mode::Weak2d,
            local_edge: 512,
            global_edge: 1024,
            steps: 200,
            warmup: 20,
            diagnostics_every: 10,
            gather_rho: false,
            parallel: false,
            decomp: None,
        }
    }
}

fn main() {
    let cfg = parse_args();
    let (universe, threading) =
        mpi::initialize_with_threading(mpi::Threading::Funneled).expect("MPI initialize failed");
    let world = universe.world();
    let rank = world.rank() as usize;

    if rank == 0 {
        println!("MPI_THREAD requested=Funneled provided={threading:?}");
    }

    match cfg.mode {
        Mode::Weak2d => run_d2::<D2Q9>(&world, &cfg),
        Mode::Weak3d | Mode::Strong3d | Mode::Diagnostics | Mode::Correctness => {
            run_d3::<D3Q19>(&world, &cfg)
        }
        Mode::Placement => report_placement(&world, &cfg),
    }

    drop(world);
    drop(universe);
}

fn parse_args() -> Config {
    let mut cfg = Config::default();
    let args: Vec<String> = std::env::args().skip(1).collect();

    if args.first().is_some_and(|a| !a.starts_with('-')) {
        cfg.local_edge = args.first().and_then(|a| a.parse().ok()).unwrap_or(512);
        cfg.steps = args.get(1).and_then(|a| a.parse().ok()).unwrap_or(200);
        return cfg;
    }

    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--mode" => {
                i += 1;
                cfg.mode = parse_mode(args.get(i).map(String::as_str));
            }
            "--local-edge" => {
                i += 1;
                cfg.local_edge = parse_usize(args.get(i), "--local-edge");
            }
            "--global-edge" => {
                i += 1;
                cfg.global_edge = parse_usize(args.get(i), "--global-edge");
            }
            "--steps" => {
                i += 1;
                cfg.steps = parse_usize(args.get(i), "--steps");
            }
            "--warmup" => {
                i += 1;
                cfg.warmup = parse_usize(args.get(i), "--warmup");
            }
            "--diagnostics-every" => {
                i += 1;
                cfg.diagnostics_every = parse_usize(args.get(i), "--diagnostics-every").max(1);
            }
            "--gather-rho" => cfg.gather_rho = true,
            "--parallel" => cfg.parallel = true,
            "--decomp" => {
                i += 1;
                cfg.decomp = Some(parse_decomp(args.get(i).map(String::as_str)));
            }
            "--help" | "-h" => {
                print_help();
                std::process::exit(0);
            }
            other => panic!("unknown argument: {other}"),
        }
        i += 1;
    }
    cfg
}

fn parse_mode(arg: Option<&str>) -> Mode {
    match arg {
        Some("weak2d") => Mode::Weak2d,
        Some("weak3d") => Mode::Weak3d,
        Some("strong3d") => Mode::Strong3d,
        Some("diagnostics") => Mode::Diagnostics,
        Some("placement") => Mode::Placement,
        Some("correctness") => Mode::Correctness,
        Some(other) => panic!("unknown --mode: {other}"),
        None => panic!("--mode requires a value"),
    }
}

fn parse_usize(arg: Option<&String>, name: &str) -> usize {
    arg.unwrap_or_else(|| panic!("{name} requires a value"))
        .parse()
        .unwrap_or_else(|_| panic!("{name} must be a positive integer"))
}

fn parse_decomp(arg: Option<&str>) -> [usize; 3] {
    let s = arg.expect("--decomp requires a value like 4x4x4");
    let parts: Vec<_> = s.split('x').collect();
    assert_eq!(parts.len(), 3, "--decomp must have form XxYxZ");
    [
        parts[0].parse().expect("bad --decomp x"),
        parts[1].parse().expect("bad --decomp y"),
        parts[2].parse().expect("bad --decomp z"),
    ]
}

fn print_help() {
    println!(
        "bench_mpi [local_edge] [steps]\n\
         bench_mpi --mode weak2d|weak3d|strong3d|diagnostics|placement|correctness \\\n\
         \t[--local-edge N] [--global-edge N] [--steps N] [--warmup N] \\\n\
         \t[--diagnostics-every N] [--gather-rho] [--parallel] [--decomp XxYxZ]"
    );
}

fn backend(parallel: bool) -> CpuScalar {
    CpuScalar {
        parallel_min_cells: if parallel { 0 } else { usize::MAX },
    }
}

fn run_d2<L: Lattice>(world: &mpi::topology::SimpleCommunicator, cfg: &Config) {
    let n = world.size() as usize;
    let decomp = cfg.decomp.unwrap_or([n, 1, 1]);
    assert_eq!(
        decomp[0] * decomp[1] * decomp[2],
        n,
        "decomp must match rank count"
    );
    let dims = [cfg.local_edge * decomp[0], cfg.local_edge * decomp[1], 1];
    let spec = GlobalSpec::<f64> {
        dims,
        nu: 0.02,
        periodic: [true, true, false],
        ..Default::default()
    };
    let mut s: MpiSolver<L, f64, CpuScalar> =
        MpiSolver::new(world, &spec, &[], &[], decomp, backend(cfg.parallel));
    let (kx, ky) = (
        2.0 * std::f64::consts::PI / dims[0] as f64,
        2.0 * std::f64::consts::PI / dims[1] as f64,
    );
    s.init_with(|x, y, _| {
        let u0 = 0.02;
        (
            1.0,
            [
                -u0 * (kx * x as f64).cos() * (ky * y as f64).sin(),
                u0 * (kx * x as f64).sin() * (ky * y as f64).cos(),
                0.0,
            ],
        )
    });
    run_workload(world, cfg, s, "weak2d", "D2Q9", dims, decomp);
}

fn run_d3<L: Lattice>(world: &mpi::topology::SimpleCommunicator, cfg: &Config) {
    let n = world.size() as usize;
    let (dims, decomp) = match cfg.mode {
        Mode::Strong3d => {
            let dims = [cfg.global_edge, cfg.global_edge, cfg.global_edge];
            let decomp = cfg.decomp.unwrap_or_else(|| choose_decomp(L::D, dims, n));
            (dims, decomp)
        }
        _ => {
            let decomp = cfg.decomp.unwrap_or_else(|| balanced_decomp(n, 3));
            let dims = [
                cfg.local_edge * decomp[0],
                cfg.local_edge * decomp[1],
                cfg.local_edge * decomp[2],
            ];
            (dims, decomp)
        }
    };
    assert_eq!(
        decomp[0] * decomp[1] * decomp[2],
        n,
        "decomp must match rank count"
    );
    let spec = GlobalSpec::<f64> {
        dims,
        nu: 0.02,
        periodic: [true, true, true],
        ..Default::default()
    };
    let mut s: MpiSolver<L, f64, CpuScalar> =
        MpiSolver::new(world, &spec, &[], &[], decomp, backend(cfg.parallel));
    let (kx, ky, kz) = (
        2.0 * std::f64::consts::PI / dims[0] as f64,
        2.0 * std::f64::consts::PI / dims[1] as f64,
        2.0 * std::f64::consts::PI / dims[2] as f64,
    );
    s.init_with(|x, y, z| {
        let u0 = 0.02;
        let sx = (kx * x as f64).sin();
        let cx = (kx * x as f64).cos();
        let sy = (ky * y as f64).sin();
        let cy = (ky * y as f64).cos();
        let sz = (kz * z as f64).sin();
        let cz = (kz * z as f64).cos();
        (
            1.0,
            [
                u0 * sx * cy * cz,
                -u0 * cx * sy * cz,
                0.5 * u0 * cx * cy * sz,
            ],
        )
    });
    let mode = match cfg.mode {
        Mode::Strong3d => "strong3d",
        Mode::Diagnostics => "diagnostics",
        Mode::Correctness => "correctness",
        _ => "weak3d",
    };
    run_workload(world, cfg, s, mode, "D3Q19", dims, decomp);
}

fn run_workload<L: Lattice>(
    world: &mpi::topology::SimpleCommunicator,
    cfg: &Config,
    mut s: MpiSolver<L, f64, CpuScalar>,
    mode: &str,
    lattice: &str,
    dims: [usize; 3],
    decomp: [usize; 3],
) {
    let rank = world.rank() as usize;
    let n = world.size() as usize;
    s.run(cfg.warmup);
    s.barrier();

    let t0 = Instant::now();
    let mut diag_time = 0.0;
    let mut diag_calls = 0usize;
    for step in 0..cfg.steps {
        s.step();
        if cfg.mode == Mode::Diagnostics && (step + 1) % cfg.diagnostics_every == 0 {
            let d0 = Instant::now();
            let _ = s.total_mass();
            let _ = s.total_momentum();
            let _ = s.nonfinite_count();
            diag_time += d0.elapsed().as_secs_f64();
            diag_calls += 1;
        }
    }
    s.barrier();
    let step_time = t0.elapsed().as_secs_f64();

    let gather_time = if cfg.gather_rho {
        s.barrier();
        let g0 = Instant::now();
        let gathered = s.gather_rho();
        s.barrier();
        if rank == 0 {
            assert_eq!(
                gathered.as_ref().map(Vec::len),
                Some(dims[0] * dims[1] * dims[2])
            );
        }
        g0.elapsed().as_secs_f64()
    } else {
        0.0
    };

    let nonfinite = s.nonfinite_count();
    assert_eq!(nonfinite, 0, "non-finite values after benchmark");
    let mass = s.total_mass();
    let momentum = s.total_momentum();

    print_rank_result(world, cfg, mode, lattice, dims);

    if rank == 0 {
        let cells = (dims[0] * dims[1] * dims[2]) as f64;
        let mlups = cells * cfg.steps as f64 / step_time / 1e6;
        let local_cells = cells / n as f64;
        println!(
            "RESULT mode={mode} lattice={lattice} ranks={n} decomp={}x{}x{} \
             local_cells={local_cells:.0} global={}x{}x{} steps={} warmup={} \
             time_s={step_time:.6} mlups_total={mlups:.3} mlups_per_rank={:.3} \
             diag_calls={diag_calls} diag_time_s={diag_time:.6} gather_rho={} \
             gather_time_s={gather_time:.6} nonfinite={nonfinite} mass={mass:.12e} \
             momentum_x={:.12e} momentum_y={:.12e} momentum_z={:.12e} \
             parallel={} rayon_threads={} hostname={}",
            decomp[0],
            decomp[1],
            decomp[2],
            dims[0],
            dims[1],
            dims[2],
            cfg.steps,
            cfg.warmup,
            mlups / n as f64,
            cfg.gather_rho,
            momentum[0],
            momentum[1],
            momentum[2],
            cfg.parallel,
            env_value("RAYON_NUM_THREADS"),
            processor_name(),
        );
    }
    drop(s);
}

fn balanced_decomp(ranks: usize, d: usize) -> [usize; 3] {
    assert!(d == 2 || d == 3);
    let mut best = [ranks, 1, 1];
    let mut best_score = usize::MAX;
    for dx in 1..=ranks {
        if ranks % dx != 0 {
            continue;
        }
        let rem = ranks / dx;
        for dy in 1..=rem {
            if rem % dy != 0 {
                continue;
            }
            let dz = rem / dy;
            if d == 2 && dz != 1 {
                continue;
            }
            let dims = if d == 2 { [dx, dy, 1] } else { [dx, dy, dz] };
            let max = dims.into_iter().max().unwrap();
            let min = dims.into_iter().filter(|v| *v > 0).min().unwrap();
            let score = max - min;
            if score < best_score {
                best_score = score;
                best = dims;
            }
        }
    }
    best
}

fn print_rank_result(
    world: &mpi::topology::SimpleCommunicator,
    cfg: &Config,
    mode: &str,
    lattice: &str,
    dims: [usize; 3],
) {
    let size = world.size() as usize;
    let rank = world.rank() as usize;
    for r in 0..size {
        world.barrier();
        if rank == r {
            println!(
                "RANK_RESULT mode={mode} lattice={lattice} rank={rank} size={size} \
                 global={}x{}x{} hostname={} affinity={} \
                 ompi_pml={} ompi_btl={} ompi_mtl={} fi_provider={} ucx_tls={} \
                 omp_num_threads={} rayon_threads={} parallel={}",
                dims[0],
                dims[1],
                dims[2],
                processor_name(),
                token(affinity()),
                env_value("OMPI_MCA_pml"),
                env_value("OMPI_MCA_btl"),
                env_value("OMPI_MCA_mtl"),
                env_value("FI_PROVIDER"),
                env_value("UCX_TLS"),
                env_value("OMP_NUM_THREADS"),
                env_value("RAYON_NUM_THREADS"),
                cfg.parallel,
            );
        }
    }
    world.barrier();
}

fn report_placement(world: &mpi::topology::SimpleCommunicator, cfg: &Config) {
    let dims = [
        cfg.local_edge.max(1),
        cfg.local_edge.max(1),
        cfg.local_edge.max(1),
    ];
    print_rank_result(world, cfg, "placement", "none", dims);
    if world.rank() == 0 {
        println!(
            "RESULT mode=placement ranks={} hostname={} mpi_vendor_hint=\"{}\" \
             ompi_pml={} ompi_btl={} ompi_mtl={} fi_provider={} ucx_tls={}",
            world.size(),
            processor_name(),
            command_output("mpirun", &["--version"]).replace('\n', " | "),
            env_value("OMPI_MCA_pml"),
            env_value("OMPI_MCA_btl"),
            env_value("OMPI_MCA_mtl"),
            env_value("FI_PROVIDER"),
            env_value("UCX_TLS"),
        );
    }
}

fn processor_name() -> String {
    mpi::environment::processor_name().unwrap_or_else(|_| {
        command_output("hostname", &[])
            .lines()
            .next()
            .unwrap_or("unknown")
            .to_string()
    })
}

fn affinity() -> String {
    let pid = std::process::id().to_string();
    let taskset = command_output("taskset", &["-pc", &pid]);
    if !taskset.is_empty() {
        return taskset.replace('\n', "; ");
    }
    let psrset = command_output("psrset", &["-q"]);
    if !psrset.is_empty() {
        return psrset.replace('\n', "; ");
    }
    "unavailable".to_string()
}

fn command_output(cmd: &str, args: &[&str]) -> String {
    match Command::new(cmd).args(args).output() {
        Ok(out) if out.status.success() => String::from_utf8_lossy(&out.stdout).trim().to_string(),
        _ => String::new(),
    }
}

fn env_value(key: &str) -> String {
    std::env::var(key).unwrap_or_else(|_| "-".to_string())
}

fn token(value: String) -> String {
    value
        .chars()
        .map(|c| {
            if c.is_ascii_alphanumeric() || matches!(c, '-' | '_' | '.' | ':' | '/' | '=') {
                c
            } else {
                '_'
            }
        })
        .collect()
}
