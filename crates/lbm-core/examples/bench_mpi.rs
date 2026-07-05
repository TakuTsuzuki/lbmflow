//! Weak-scaling benchmark: a fixed 512² D2Q9 f64 block **per rank**, strip
//! decomposition `[n, 1, 1]` (docs/MPI_GUIDE.md; driven by
//! scripts/bench_mpi.sh, which aggregates the RESULT lines into a table).
//!
//! Measurement hygiene:
//!
//! - The backend is forced serial per rank (`parallel_min_cells = MAX`), so
//!   ranks scale over cores 1:1 and the halo-exchange overhead is not hidden
//!   inside rayon scheduling noise. Production runs would use MPI × rayon
//!   hybrid; that mapping is a cluster-tuning question (docs/MPI_GUIDE.md).
//! - On a single machine every rank shares one memory system, so these
//!   numbers measure *shared-memory* MPI (Open MPI sm/vader BTL), not a real
//!   interconnect. Treat them as an upper bound on locality and a functional
//!   check of the overlap structure; true weak scaling needs a cluster
//!   (COMPETITIVE_SPEC §5, R3: ≥80% at 64 ranks).
//!
//! Usage: `mpirun -n <ranks> bench_mpi [local_edge] [steps]`
//! (defaults 512, 200; 20 warm-up steps are excluded from timing).

use lbm_core::dist::MpiSolver;
use lbm_core::lattice::D2Q9;
use lbm_core::prelude::*;
use mpi::traits::*;
use std::time::Instant;

fn main() {
    let universe = mpi::initialize().expect("MPI initialize failed");
    let world = universe.world();
    let n = world.size() as usize;
    let rank = world.rank() as usize;
    let mut args = std::env::args().skip(1);
    let local: usize = args.next().and_then(|a| a.parse().ok()).unwrap_or(512);
    let steps: usize = args.next().and_then(|a| a.parse().ok()).unwrap_or(200);
    let warmup = 20usize;

    let dims = [local * n, local, 1];
    let spec = GlobalSpec::<f64> {
        dims,
        nu: 0.02,
        periodic: [true, true, false],
        ..Default::default()
    };
    // Serial backend per rank: rank count = core count, no oversubscription.
    let backend = CpuScalar {
        parallel_min_cells: usize::MAX,
    };
    let mut s: MpiSolver<D2Q9, f64, CpuScalar> =
        MpiSolver::new(&world, &spec, &[], &[], [n, 1, 1], backend);
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

    s.run(warmup);
    s.barrier();
    let t0 = Instant::now();
    s.run(steps);
    s.barrier();
    let dt = t0.elapsed().as_secs_f64();

    // Functional sanity so a silently-corrupt run cannot "scale well".
    assert_eq!(s.nonfinite_count(), 0, "non-finite values after benchmark");
    let mass = s.total_mass();

    if rank == 0 {
        let cells = (dims[0] * dims[1]) as f64;
        let mlups = cells * steps as f64 / dt / 1e6;
        println!(
            "RESULT ranks={n} local={local}x{local} global={}x{} steps={steps} \
             time_s={dt:.3} mlups_total={mlups:.1} mass={mass:.6e}",
            dims[0], dims[1]
        );
    }
    // The solver owns duplicated communicators: free them before finalize.
    drop(s);
    drop(world);
    drop(universe);
}
