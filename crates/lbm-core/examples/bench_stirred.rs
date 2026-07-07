//! Landed-physics stirred-tank degradation benchmark (ME-4a Stage A).
//!
//! Run:
//! `cargo run --release -p lbm-core --example bench_stirred -- [n] [steps] [warmup] [--gpu]`
//!
//! Defaults are a quick smoke shape (`n=64`, heuristic step count). Published
//! numbers must use the protocol printed by this example: quiet window, A/B/A
//! interleave, and 5-run median.

use lbm_core::particles::{sample_grid, Particle, ParticleSet};
use lbm_core::prelude::*;
use std::f64::consts::PI;
use std::time::Instant;

#[cfg(feature = "gpu")]
use lbm_core::gpu::{GpuStorage, KernelCfg};

const NU: f64 = 0.02;
const U_TIP: f64 = 0.04;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Row {
    B0,
    C1,
    C2,
    C3,
    C4,
}

impl Row {
    fn all() -> [Self; 5] {
        [Self::B0, Self::C1, Self::C2, Self::C3, Self::C4]
    }

    fn collision(self) -> CollisionKind {
        match self {
            Self::B0 => CollisionKind::Trt {
                magic: CollisionKind::MAGIC_STD,
            },
            Self::C1 | Self::C2 | Self::C3 | Self::C4 => CollisionKind::CentralMoment {
                omega_shear: 1.0 / (3.0 * NU + 0.5),
            },
        }
    }

    fn label(self) -> &'static str {
        match self {
            Self::B0 => "B0",
            Self::C1 => "C1",
            Self::C2 => "C2",
            Self::C3 => "C3",
            Self::C4 => "C4",
        }
    }

    fn tag(self) -> &'static str {
        match self {
            Self::B0 => "ratio[baseline]",
            Self::C1 => "ratio[cm]",
            Self::C2 => "ratio[cm+wale]",
            Self::C3 => "ratio[cm+wale+ibm]",
            Self::C4 => "ratio[cm+wale+ibm+part]",
        }
    }

    fn wale(self) -> bool {
        matches!(self, Self::C2 | Self::C3 | Self::C4)
    }

    fn ibm(self) -> bool {
        matches!(self, Self::C3 | Self::C4)
    }

    fn particles(self) -> bool {
        matches!(self, Self::C4)
    }
}

#[derive(Clone, Copy, Debug)]
struct Options {
    n: usize,
    steps: usize,
    warmup: usize,
    particles: usize,
    run_gpu: bool,
}

impl Options {
    fn parse() -> Self {
        let mut positional = Vec::new();
        let mut particles = 1024usize;
        let mut run_gpu = false;
        let mut args = std::env::args().skip(1);
        while let Some(arg) = args.next() {
            match arg.as_str() {
                "--gpu" => run_gpu = true,
                "--particles" => {
                    particles = args
                        .next()
                        .and_then(|s| s.parse().ok())
                        .expect("--particles requires a positive integer");
                }
                "--help" | "-h" => {
                    print_help();
                    std::process::exit(0);
                }
                _ => positional.push(arg),
            }
        }

        let n = positional
            .first()
            .and_then(|s| s.parse().ok())
            .unwrap_or(64usize);
        let cells = n * n * n;
        let heuristic_steps = (5_000_000usize / cells.max(1)).clamp(20, 500);
        let steps = positional
            .get(1)
            .and_then(|s| s.parse().ok())
            .unwrap_or(heuristic_steps);
        let warmup = positional
            .get(2)
            .and_then(|s| s.parse().ok())
            .unwrap_or(10usize);
        Self {
            n,
            steps,
            warmup,
            particles,
            run_gpu,
        }
    }
}

fn print_help() {
    println!("bench_stirred [n] [steps] [warmup] [--gpu] [--particles count]");
    println!("  n: published grids are 96, 128, 192; default 64; n=32 is useful for smoke");
    println!("  steps: measured steps after warmup; default uses bench_backends-style scaling");
    println!("  warmup: excluded warmup steps; default 10");
    println!("  --gpu: run GPU f32 collision-only B0->C1 rows when built with --features gpu");
}

/// Rushton-turbine geometry lifted from crates/lbm-cli/examples/stirred_tank_3d.rs.
struct Geom {
    n: usize,
    cx: f64,
    cy: f64,
    r_tank: f64,
    zc: f64,
    tip_r: f64,
    disk_r: f64,
    hub_r: f64,
    shaft_r: f64,
    blade_hh: f64,
    disk_hh: f64,
    n_blades: usize,
    blade_hw: f64,
    baffle_len: f64,
    baffle_hw: f64,
}

impl Geom {
    fn new(n: usize) -> Self {
        let r_tank = n as f64 / 2.0 - 3.0;
        let tip_r = r_tank / 3.0;
        Self {
            n,
            cx: (n as f64 - 1.0) / 2.0,
            cy: (n as f64 - 1.0) / 2.0,
            r_tank,
            zc: n as f64 * 0.35,
            tip_r,
            disk_r: tip_r * 0.66,
            hub_r: tip_r * 0.22,
            shaft_r: (tip_r * 0.12).max(1.5),
            blade_hh: tip_r * 0.30,
            disk_hh: 1.2,
            n_blades: 6,
            blade_hw: 1.2,
            baffle_len: r_tank * 0.2,
            baffle_hw: 1.5,
        }
    }

    fn rad(&self, x: usize, y: usize) -> f64 {
        let dx = x as f64 - self.cx;
        let dy = y as f64 - self.cy;
        (dx * dx + dy * dy).sqrt()
    }

    fn is_solid(&self, x: usize, y: usize, z: usize) -> bool {
        let n = self.n;
        if x == 0 || x == n - 1 || y == 0 || y == n - 1 || z == 0 || z == n - 1 {
            return true;
        }
        let r = self.rad(x, y);
        if r > self.r_tank {
            return true;
        }
        let dx = x as f64 - self.cx;
        let dy = y as f64 - self.cy;
        for k in 0..4 {
            let beta = k as f64 * PI / 2.0;
            let s = dx * beta.cos() + dy * beta.sin();
            let p = -dx * beta.sin() + dy * beta.cos();
            if s > self.r_tank - self.baffle_len && s <= self.r_tank && p.abs() <= self.baffle_hw {
                return true;
            }
        }
        false
    }

    fn in_turbine(&self, x: usize, y: usize, z: usize, theta: f64) -> bool {
        let dx = x as f64 - self.cx;
        let dy = y as f64 - self.cy;
        let r = (dx * dx + dy * dy).sqrt();
        let zf = z as f64;
        if r <= self.shaft_r && zf >= self.zc {
            return true;
        }
        if (zf - self.zc).abs() > self.blade_hh {
            return false;
        }
        if r <= self.hub_r {
            return true;
        }
        if r <= self.disk_r && (zf - self.zc).abs() <= self.disk_hh {
            return true;
        }
        if r >= self.disk_r * 0.85 && r <= self.tip_r {
            let phi = dy.atan2(dx);
            for b in 0..self.n_blades {
                let beta = theta + b as f64 * 2.0 * PI / self.n_blades as f64;
                let d = (r * (phi - beta).sin()).abs();
                if d <= self.blade_hw && (phi - beta).cos() > 0.0 {
                    return true;
                }
            }
        }
        false
    }

    fn rotor_body(&self, omega: f64, theta: f64) -> RotatingBody {
        let mut markers = Vec::new();
        for z in 1..self.n - 1 {
            for y in 1..self.n - 1 {
                for x in 1..self.n - 1 {
                    if self.is_solid(x, y, z) || !self.in_turbine(x, y, z, theta) {
                        continue;
                    }
                    markers.push(IbmMarker {
                        position: [x as f64, y as f64, z as f64],
                        weight: 1.0,
                    });
                }
            }
        }
        RotatingBody::from_markers([self.cx, self.cy, self.zc], [0.0, 0.0, omega], markers)
    }
}

fn spec<T: Real>(n: usize, collision: CollisionKind) -> GlobalSpec<T> {
    GlobalSpec {
        dims: [n, n, n],
        nu: NU,
        periodic: [false, false, false],
        collision,
        ..Default::default()
    }
}

fn periodic_spec<T: Real, L: Lattice>(n: usize, collision: CollisionKind) -> GlobalSpec<T> {
    let dims = if L::D == 2 { [n, n, 1] } else { [n, n, n] };
    GlobalSpec {
        dims,
        nu: NU,
        periodic: [true, true, L::D == 3],
        collision,
        ..Default::default()
    }
}

fn init<T: Real>(n: usize) -> impl Fn(usize, usize, usize) -> (T, [T; 3]) + Copy {
    move |x, y, z| {
        let k = 2.0 * PI / n as f64;
        let (xf, yf, zf) = (k * x as f64, k * y as f64, k * z as f64);
        (
            T::one(),
            [
                T::r(0.02 * yf.sin()),
                T::r(0.02 * (2.0 * xf).sin()),
                T::r(0.01 * zf.sin()),
            ],
        )
    }
}

fn build_stirred_solver<T: Real>(row: Row, n: usize) -> Solver<D3Q19, T, CpuSimd, LocalPeriodic>
where
    CpuSimd: Backend<D3Q19, T, Fields = SoaFields<T>>,
{
    let g = Geom::new(n);
    let mut walls = WallSpec::<T>::default();
    for f in [
        Face::XNeg,
        Face::XPos,
        Face::YNeg,
        Face::YPos,
        Face::ZNeg,
        Face::ZPos,
    ] {
        walls.is_wall[f.index()] = true;
    }
    let spec = spec::<T>(n, row.collision());
    let (solid, wall_u) = build_wall_rims(3, spec.dims, &walls);
    let mut solver = Solver::new(
        &spec,
        &solid,
        &wall_u,
        [1, 1, 1],
        CpuSimd::default(),
        LocalPeriodic,
    );
    solver.init_with(init::<T>(n));
    for z in 0..n {
        for y in 0..n {
            for x in 0..n {
                let rim = x == 0 || x == n - 1 || y == 0 || y == n - 1 || z == 0 || z == n - 1;
                if !rim && g.is_solid(x, y, z) {
                    solver.set_solid(x, y, z);
                }
            }
        }
    }
    solver.mark_masks_dirty();
    solver
}

fn build_particles(g: &Geom, count: usize) -> ParticleSet {
    let mut particles = Vec::with_capacity(count);
    let mut ix = 0usize;
    let mut iy = 0usize;
    let z = (g.zc + g.tip_r * 0.7).min((g.n - 2) as f64);
    while particles.len() < count {
        let fx = (ix as f64 + 0.5) / 32.0;
        let fy = (iy as f64 + 0.5) / 32.0;
        let x = 1.0 + fx * (g.n as f64 - 2.0);
        let y = 1.0 + fy * (g.n as f64 - 2.0);
        let xi = x.round().clamp(1.0, (g.n - 2) as f64) as usize;
        let yi = y.round().clamp(1.0, (g.n - 2) as f64) as usize;
        if !g.is_solid(xi, yi, z as usize) {
            particles.push(Particle {
                pos: [x, y, z],
                vel: [0.0; 3],
                d: 0.2,
                rho_p: 1.0,
                exposure: 0.0,
            });
        }
        ix += 1;
        if ix == 32 {
            ix = 0;
            iy = (iy + 1) % 32;
        }
    }
    ParticleSet::new(particles, 1.0, NU, [0.0; 3])
}

fn bench_stirred<T: Real>(row: Row, options: Options) -> f64
where
    CpuSimd: Backend<D3Q19, T, Fields = SoaFields<T>>,
{
    let g = Geom::new(options.n);
    let omega = U_TIP / g.tip_r;
    let mut solver = build_stirred_solver::<T>(row, options.n);
    let mut wale = WaleLes::<T>::new();
    let ibm_cfg = DirectForcingConfig::default();
    let mut particles = row
        .particles()
        .then(|| build_particles(&g, options.particles));

    for step in 0..options.warmup {
        apply_row_step(
            row,
            &g,
            omega,
            step,
            &mut solver,
            &mut wale,
            ibm_cfg,
            &mut particles,
        );
    }

    let cells = options.n * options.n * options.n;
    let t0 = Instant::now();
    for step in options.warmup..options.warmup + options.steps {
        apply_row_step(
            row,
            &g,
            omega,
            step,
            &mut solver,
            &mut wale,
            ibm_cfg,
            &mut particles,
        );
    }
    let mlups = (cells * options.steps) as f64 / t0.elapsed().as_secs_f64() / 1e6;
    assert!(
        solver.total_mass_f64().is_finite(),
        "non-finite mass in stirred benchmark"
    );
    mlups
}

fn apply_row_step<T: Real>(
    row: Row,
    g: &Geom,
    omega: f64,
    step: usize,
    solver: &mut Solver<D3Q19, T, CpuSimd, LocalPeriodic>,
    wale: &mut WaleLes<T>,
    ibm_cfg: DirectForcingConfig,
    particles: &mut Option<ParticleSet>,
) where
    CpuSimd: Backend<D3Q19, T, Fields = SoaFields<T>>,
{
    if row.wale() {
        wale.update(solver);
    }
    if row.ibm() {
        let body = g.rotor_body(omega, omega * step as f64);
        let diag = solver.apply_rotating_ibm(&body, ibm_cfg);
        std::hint::black_box(diag);
    }
    solver.step();
    if let Some(particles) = particles {
        let ux = solver.gather_ux();
        let uy = solver.gather_uy();
        let uz = solver.gather_uz();
        let dims = [g.n, g.n, g.n];
        particles
            .step(
                |pos| {
                    sample_grid(pos, dims, |x, y, z| {
                        let i = (z * g.n + y) * g.n + x;
                        (
                            [ux[i].as_f64(), uy[i].as_f64(), uz[i].as_f64()],
                            g.is_solid(x, y, z),
                        )
                    })
                },
                None::<fn([f64; 3]) -> f64>,
            )
            .expect("particle step stayed inside Schiller-Naumann validity");
    }
}

fn bench_d3q27<T: Real>(row: Row, n: usize, warmup: usize, steps: usize) -> f64
where
    CpuSimd: Backend<D3Q27, T, Fields = SoaFields<T>>,
{
    let spec = periodic_spec::<T, D3Q27>(n, row.collision());
    let mut solver: Solver<D3Q27, T, CpuSimd, LocalPeriodic> = Solver::new(
        &spec,
        &[],
        &[],
        [1, 1, 1],
        CpuSimd::default(),
        LocalPeriodic,
    );
    solver.init_with(init::<T>(n));
    solver.run(warmup);
    let cells = n * n * n;
    let t0 = Instant::now();
    solver.run(steps);
    let mlups = (cells * steps) as f64 / t0.elapsed().as_secs_f64() / 1e6;
    assert!(
        solver.total_mass_f64().is_finite(),
        "non-finite mass in D3Q27 benchmark"
    );
    mlups
}

fn print_result(
    row: &str,
    lattice: &str,
    backend: &str,
    precision: &str,
    n: usize,
    tag: &str,
    mlups: f64,
    baseline: f64,
) {
    let ratio = mlups / baseline;
    println!(
        "row={row} grid={n}^3 lattice={lattice} backend={backend} precision={precision} \
         tag={tag} mlups={mlups:.3} ratio_vs_baseline={ratio:.4} slowdown={:.2}x",
        1.0 / ratio
    );
}

#[cfg(feature = "gpu")]
fn run_gpu_d3q19(options: Options) {
    if !options.run_gpu {
        println!(
            "BENCH-PENDING gpu: pass --gpu to run B0->C1 f32 collision-only rows; \
             no GPU C4 headline is emitted because WALE-update and IBM are CPU-side host ops."
        );
        return;
    }
    let ctx = match GpuContext::new() {
        Ok(ctx) => ctx,
        Err(e) => {
            println!("BENCH-PENDING gpu: no usable adapter ({e})");
            return;
        }
    };
    println!(
        "gpu_adapter={} backend={:?}",
        ctx.adapter_info.name, ctx.adapter_info.backend
    );
    let b0 = bench_gpu_stirred(Row::B0, options, &ctx);
    print_result(
        "B0",
        "D3Q19",
        "Wgpu",
        "f32",
        options.n,
        Row::B0.tag(),
        b0,
        b0,
    );
    let c1 = bench_gpu_stirred(Row::C1, options, &ctx);
    print_result(
        "C1",
        "D3Q19",
        "Wgpu",
        "f32",
        options.n,
        Row::C1.tag(),
        c1,
        b0,
    );
    println!(
        "gpu_caveat=collision-only: refusing GPU C4 ratio; WALE-update and apply_rotating_ibm \
         are not fused GPU kernels in landed code."
    );
}

#[cfg(feature = "gpu")]
fn bench_gpu_stirred(row: Row, options: Options, ctx: &std::sync::Arc<GpuContext>) -> f64 {
    let g = Geom::new(options.n);
    let spec = spec::<f32>(options.n, row.collision());
    let mut walls = WallSpec::<f32>::default();
    for f in [
        Face::XNeg,
        Face::XPos,
        Face::YNeg,
        Face::YPos,
        Face::ZNeg,
        Face::ZPos,
    ] {
        walls.is_wall[f.index()] = true;
    }
    let (solid, wall_u) = build_wall_rims(3, spec.dims, &walls);
    let backend = WgpuBackend::<D3Q19>::with_config(
        ctx.clone(),
        KernelCfg {
            storage: GpuStorage::F32,
        },
    );
    let mut solver: Solver<D3Q19, f32, WgpuBackend<D3Q19>, LocalPeriodic> =
        Solver::new(&spec, &solid, &wall_u, [1, 1, 1], backend, LocalPeriodic);
    solver.init_with(init::<f32>(options.n));
    for z in 0..options.n {
        for y in 0..options.n {
            for x in 0..options.n {
                let rim = x == 0
                    || x == options.n - 1
                    || y == 0
                    || y == options.n - 1
                    || z == 0
                    || z == options.n - 1;
                if !rim && g.is_solid(x, y, z) {
                    solver.set_solid(x, y, z);
                }
            }
        }
    }
    solver.mark_masks_dirty();
    solver.run(options.warmup);
    solver.backend().context().wait_idle();
    let cells = options.n * options.n * options.n;
    let t0 = Instant::now();
    solver.run(options.steps);
    solver.backend().context().wait_idle();
    let mlups = (cells * options.steps) as f64 / t0.elapsed().as_secs_f64() / 1e6;
    assert!(
        solver.total_mass_f64().is_finite(),
        "non-finite mass in GPU benchmark"
    );
    mlups
}

#[cfg(not(feature = "gpu"))]
fn run_gpu_d3q19(options: Options) {
    if options.run_gpu {
        println!(
            "BENCH-PENDING gpu: --gpu was requested, but this binary was built without \
             --features gpu"
        );
    } else {
        println!(
            "BENCH-PENDING gpu: rebuild with --features gpu and pass --gpu for B0->C1 f32 \
             collision-only rows; no GPU C4 headline is emitted because WALE-update and IBM \
             are CPU-side host ops."
        );
    }
}

fn main() {
    let options = Options::parse();
    println!("# ME-4a landed-physics stirred-tank degradation benchmark");
    println!(
        "protocol=quiet-window A/B/A interleave 5-run-median warmup-excluded; \
         single-shot example run only"
    );
    println!(
        "coverage=landed-physics-subset no resolved multiphase no scalar ADE no two-way particles"
    );
    println!(
        "gpu_policy=B0->C1 collision-only; GPU C4 headline ratio is refused because WALE-update \
         and apply_rotating_ibm are CPU-side host operations"
    );
    println!(
        "config grid={}^3 steps={} warmup={} particles={} nu={} u_tip={}",
        options.n, options.steps, options.warmup, options.particles, NU, U_TIP
    );
    if !matches!(options.n, 96 | 128 | 192) {
        println!(
            "grid_note={}^3 is a smoke/development grid; published grids are 96^3, 128^3, 192^3",
            options.n
        );
    }

    let mut baseline_f32 = None;
    for row in Row::all() {
        let mlups = bench_stirred::<f32>(row, options);
        let base = *baseline_f32.get_or_insert(mlups);
        print_result(
            row.label(),
            "D3Q19",
            "CpuSimd",
            "f32",
            options.n,
            row.tag(),
            mlups,
            base,
        );
    }

    let mut baseline_f64 = None;
    for row in Row::all() {
        let mlups = bench_stirred::<f64>(row, options);
        let base = *baseline_f64.get_or_insert(mlups);
        print_result(
            row.label(),
            "D3Q19",
            "CpuSimd",
            "f64",
            options.n,
            row.tag(),
            mlups,
            base,
        );
    }

    run_gpu_d3q19(options);

    let d27_b0_f32 = bench_d3q27::<f32>(Row::B0, options.n, options.warmup, options.steps);
    print_result(
        "D27-B0",
        "D3Q27",
        "CpuSimd",
        "f32",
        options.n,
        "ratio[d3q27-baseline]",
        d27_b0_f32,
        d27_b0_f32,
    );
    let d27_c1_f32 = bench_d3q27::<f32>(Row::C1, options.n, options.warmup, options.steps);
    print_result(
        "D27-C1",
        "D3Q27",
        "CpuSimd",
        "f32",
        options.n,
        "ratio[d3q27-cm]",
        d27_c1_f32,
        d27_b0_f32,
    );
    let d27_b0_f64 = bench_d3q27::<f64>(Row::B0, options.n, options.warmup, options.steps);
    print_result(
        "D27-B0",
        "D3Q27",
        "CpuSimd",
        "f64",
        options.n,
        "ratio[d3q27-baseline]",
        d27_b0_f64,
        d27_b0_f64,
    );
    let d27_c1_f64 = bench_d3q27::<f64>(Row::C1, options.n, options.warmup, options.steps);
    print_result(
        "D27-C1",
        "D3Q27",
        "CpuSimd",
        "f64",
        options.n,
        "ratio[d3q27-cm]",
        d27_c1_f64,
        d27_b0_f64,
    );
}
