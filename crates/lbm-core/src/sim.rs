//! The core simulation loop: collide → stream → open-edge BCs → moments.
//!
//! Design notes:
//! - Memory layout is plane-major SoA: `f[q*n + y*nx + x]` (one contiguous
//!   plane per direction). The row loops are decomposed into solid-free spans
//!   (via per-row solid runs) so the hot kernels are branch-free and
//!   auto-vectorize; streaming interior spans is a plain shifted copy.
//! - Streaming uses the pull scheme; wall edges are one-cell solid rims so
//!   half-way bounce-back handles them uniformly (no corner special cases).
//! - Macroscopic fields stored in `rho/ux/uy` always describe the *current*
//!   post-step state; velocities include the Guo half-force correction.

use crate::domain::{Collision, Edge, EdgeBC, Edges, SimConfig, MAX_SPEED};
use crate::lattice::{dir_index, CX, CY, OPP, Q, W};
use crate::real::Real;

#[cfg(feature = "parallel")]
use rayon::prelude::*;

/// Per-step constants handed to the inner loops.
#[derive(Clone, Copy)]
struct Params<T: Real> {
    omega_p: T,
    omega_m: T,
    /// `1 - omega_p/2` (Guo forcing prefactor, symmetric part).
    cp: T,
    /// `1 - omega_m/2` (Guo forcing prefactor, antisymmetric part).
    cm: T,
    fx: T,
    fy: T,
    cxr: [T; Q],
    cyr: [T; Q],
    wr: [T; Q],
}

#[derive(Clone, Copy)]
struct Geom {
    nx: usize,
    ny: usize,
    per_x: bool,
    per_y: bool,
}

enum ZouHe<T: Real> {
    Velocity([T; 2]),
    Pressure(T),
}

/// Below this many cells the row-parallel loops run serially: rayon's
/// dispatch overhead (~100 µs/step measured on an 18-core M-series) dwarfs
/// the actual work on small grids.
pub const PARALLEL_MIN_CELLS: usize = 16_384;

/// Row-sliced mutable access to a plane-major population buffer
/// (`buf[q*n + y*nx + x]`) shared across row-parallel tasks.
///
/// Soundness contract: parallel tasks partition the domain by `y`, and every
/// task only requests rows for its own `y`, so the `(q, y)` row slices handed
/// out are pairwise disjoint across live borrows.
struct PlaneRows<'a, T> {
    ptr: *mut T,
    n: usize,
    nx: usize,
    ny: usize,
    _marker: std::marker::PhantomData<&'a mut [T]>,
}

unsafe impl<T: Send> Send for PlaneRows<'_, T> {}
unsafe impl<T: Send> Sync for PlaneRows<'_, T> {}

impl<'a, T> PlaneRows<'a, T> {
    fn new(buf: &'a mut [T], nx: usize, ny: usize) -> Self {
        debug_assert_eq!(buf.len(), nx * ny * Q);
        Self {
            ptr: buf.as_mut_ptr(),
            n: nx * ny,
            nx,
            ny,
            _marker: std::marker::PhantomData,
        }
    }

    /// Row `y` of plane `q`.
    ///
    /// # Safety
    /// The caller must guarantee no other live slice overlaps `(q, y)`.
    #[inline]
    unsafe fn row(&self, q: usize, y: usize) -> &'a mut [T] {
        debug_assert!(q < Q && y < self.ny);
        std::slice::from_raw_parts_mut(self.ptr.add(q * self.n + y * self.nx), self.nx)
    }

}

/// Row-sliced mutable access to a single `[y][x]` field, shared across
/// row-parallel tasks under the same soundness contract as [`PlaneRows`]
/// (each task touches only its own `y`). A plain range `par_iter` over rows
/// with these wrappers parallelizes markedly better than rayon's three-way
/// `zip` of chunk iterators.
struct FieldRows<'a, T> {
    ptr: *mut T,
    nx: usize,
    ny: usize,
    _marker: std::marker::PhantomData<&'a mut [T]>,
}

unsafe impl<T: Send> Send for FieldRows<'_, T> {}
unsafe impl<T: Send> Sync for FieldRows<'_, T> {}

impl<'a, T> FieldRows<'a, T> {
    fn new(buf: &'a mut [T], nx: usize, ny: usize) -> Self {
        debug_assert_eq!(buf.len(), nx * ny);
        Self {
            ptr: buf.as_mut_ptr(),
            nx,
            ny,
            _marker: std::marker::PhantomData,
        }
    }

    /// Row `y`.
    ///
    /// # Safety
    /// The caller must guarantee no other live slice overlaps row `y`.
    #[inline]
    unsafe fn row(&self, y: usize) -> &'a mut [T] {
        debug_assert!(y < self.ny);
        std::slice::from_raw_parts_mut(self.ptr.add(y * self.nx), self.nx)
    }
}

/// D2Q9 lattice Boltzmann simulation on a rectangular grid.
pub struct Simulation<T: Real> {
    nx: usize,
    ny: usize,
    f: Vec<T>,
    ftmp: Vec<T>,
    rho: Vec<T>,
    ux: Vec<T>,
    uy: Vec<T>,
    /// Double buffers for the fused pass: the moments of step k are written
    /// here while the in-flight collide stage still reads step k-1's moments
    /// from `rho/ux/uy`; the pairs are swapped right after the pass. Solid
    /// cells are never written by the pass, so both buffers keep their
    /// values in sync via `init_with`/`set_solid`.
    rho2: Vec<T>,
    ux2: Vec<T>,
    uy2: Vec<T>,
    /// `ConvectiveOutflow` memory term: the post-collide unknown populations
    /// of the previous step at each edge cell (`[coord*3 + slot]`, slots
    /// ordered `[n, n+t, n-t]`). The fused pass no longer materialises
    /// post-collide populations, so this is captured explicitly
    /// (`capture_conv_stale`) into `conv_stale_next` and swapped in after
    /// the BC pass. Zero-initialised, matching the pre-fusion behaviour of
    /// reading the build-time all-zero deviation buffer on the first step.
    conv_stale: [Option<Vec<T>>; 4],
    conv_stale_next: [Option<Vec<T>>; 4],
    solid: Vec<bool>,
    /// Per-row maximal runs of consecutive solid cells, `[start, end)` in x.
    /// Lets the hot loops process the solid-free spans branch-free.
    solid_runs: Vec<Vec<(u32, u32)>>,
    /// Set by solid-mask mutations; `solid_runs` is rebuilt lazily.
    solid_dirty: bool,
    /// Wall velocity per cell; only meaningful for solid cells.
    wall_u: Vec<[T; 2]>,
    edges: Edges<T>,
    omega_p: f64,
    omega_m: f64,
    nu: f64,
    force: [T; 2],
    probe: Option<Vec<bool>>,
    probed_force: [T; 2],
    /// Per-edge inlet velocity profiles, indexed by [`Edge::index`]; overrides
    /// the uniform `VelocityInlet` velocity when set.
    inlet_profiles: [Option<Vec<[T; 2]>>; 4],
    /// Optional per-cell body force, added to the uniform `force` (used by
    /// multiphase models; rewritten every step by the caller).
    force_field: Option<Vec<[T; 2]>>,
    time: u64,
    /// Only consulted by the parallel dispatch (serial builds always run one
    /// band).
    #[cfg_attr(not(feature = "parallel"), allow(dead_code))]
    use_parallel: bool,
}

impl<T: Real> std::fmt::Debug for Simulation<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Simulation")
            .field("nx", &self.nx)
            .field("ny", &self.ny)
            .field("nu", &self.nu)
            .field("time", &self.time)
            .finish_non_exhaustive()
    }
}

impl<T: Real> Simulation<T> {
    pub(crate) fn from_config(cfg: SimConfig<T>) -> Self {
        let n = cfg.nx * cfg.ny;
        let tau = 3.0 * cfg.nu + 0.5;
        let omega_p = 1.0 / tau;
        let omega_m = match cfg.collision {
            Collision::Bgk => omega_p,
            Collision::Trt { magic } => {
                let lam_p = tau - 0.5;
                1.0 / (magic / lam_p + 0.5)
            }
        };
        // Deviation storage: `f` holds f_q - w_q, so the quiescent state
        // (rho = 1, u = 0) is exactly all-zero. All arithmetic then acts on
        // small fluctuations, which lifts the effective precision of f32
        // runs by orders of magnitude (docs/PHYSICS.md).
        let f = vec![T::zero(); n * Q];
        let conv_len = |edge: Edge| match edge {
            Edge::Left | Edge::Right => cfg.ny,
            Edge::Bottom | Edge::Top => cfg.nx,
        };
        let conv_buf = |bc: EdgeBC<T>, edge: Edge| {
            matches!(bc, EdgeBC::ConvectiveOutflow { .. })
                .then(|| vec![T::zero(); 3 * conv_len(edge)])
        };
        let conv_stale = [
            conv_buf(cfg.edges.left, Edge::Left),
            conv_buf(cfg.edges.right, Edge::Right),
            conv_buf(cfg.edges.bottom, Edge::Bottom),
            conv_buf(cfg.edges.top, Edge::Top),
        ];
        let mut sim = Self {
            nx: cfg.nx,
            ny: cfg.ny,
            ftmp: f.clone(),
            f,
            rho: vec![T::one(); n],
            ux: vec![T::zero(); n],
            uy: vec![T::zero(); n],
            rho2: vec![T::one(); n],
            ux2: vec![T::zero(); n],
            uy2: vec![T::zero(); n],
            conv_stale_next: conv_stale.clone(),
            conv_stale,
            solid: vec![false; n],
            solid_runs: Vec::new(),
            solid_dirty: true,
            wall_u: vec![[T::zero(); 2]; n],
            edges: cfg.edges,
            omega_p,
            omega_m,
            nu: cfg.nu,
            force: cfg.force,
            probe: None,
            probed_force: [T::zero(); 2],
            inlet_profiles: [None, None, None, None],
            force_field: None,
            time: 0,
            use_parallel: cfg!(feature = "parallel") && n >= PARALLEL_MIN_CELLS,
        };
        sim.build_rims();
        sim.update_moments();
        sim
    }

    /// Realise wall-type edges as one-cell solid rims.
    ///
    /// Where two rims share a corner cell, the faster wall's velocity wins,
    /// so the result does not depend on the order edges are applied and the
    /// whole setup stays exactly equivariant under rotations/mirrors (a
    /// lid-driven cavity gives identical flows for all four lid
    /// orientations). Equal speeds keep the first-applied edge
    /// (bottom, top, left, right order) — only relevant when two moving
    /// walls of identical speed meet at a corner.
    fn build_rims(&mut self) {
        let (nx, ny) = (self.nx, self.ny);
        let mut best_speed = vec![-1.0f64; nx * ny];
        let mut rim = |cells: Box<dyn Iterator<Item = usize>>, bc: EdgeBC<T>| {
            let u = match bc {
                EdgeBC::MovingWall { u } => u,
                _ => [T::zero(); 2],
            };
            let speed = u[0].as_f64().powi(2) + u[1].as_f64().powi(2);
            if bc.is_wall() {
                for i in cells {
                    self.solid[i] = true;
                    if speed > best_speed[i] {
                        best_speed[i] = speed;
                        self.wall_u[i] = u;
                    }
                }
            }
        };
        rim(Box::new(0..nx), self.edges.bottom);
        rim(
            Box::new((0..nx).map(move |x| (ny - 1) * nx + x)),
            self.edges.top,
        );
        rim(Box::new((0..ny).map(move |y| y * nx)), self.edges.left);
        rim(
            Box::new((0..ny).map(move |y| y * nx + nx - 1)),
            self.edges.right,
        );
        self.solid_dirty = true;
    }

    fn rebuild_solid_runs(&mut self) {
        let (nx, ny) = (self.nx, self.ny);
        self.solid_runs.resize(ny, Vec::new());
        for y in 0..ny {
            let row = &self.solid[y * nx..(y + 1) * nx];
            let runs = &mut self.solid_runs[y];
            runs.clear();
            let mut x = 0;
            while x < nx {
                if row[x] {
                    let start = x;
                    while x < nx && row[x] {
                        x += 1;
                    }
                    runs.push((start as u32, x as u32));
                } else {
                    x += 1;
                }
            }
        }
        self.solid_dirty = false;
    }

    fn params(&self) -> Params<T> {
        let mut cxr = [T::zero(); Q];
        let mut cyr = [T::zero(); Q];
        let mut wr = [T::zero(); Q];
        for q in 0..Q {
            cxr[q] = T::r(CX[q] as f64);
            cyr[q] = T::r(CY[q] as f64);
            wr[q] = T::r(W[q]);
        }
        Params {
            omega_p: T::r(self.omega_p),
            omega_m: T::r(self.omega_m),
            cp: T::r(1.0 - self.omega_p / 2.0),
            cm: T::r(1.0 - self.omega_m / 2.0),
            fx: self.force[0],
            fy: self.force[1],
            cxr,
            cyr,
            wr,
        }
    }

    fn geom(&self) -> Geom {
        Geom {
            nx: self.nx,
            ny: self.ny,
            per_x: self.edges.left.is_periodic(),
            per_y: self.edges.bottom.is_periodic(),
        }
    }

    // ------------------------------------------------------------------
    // Time stepping
    // ------------------------------------------------------------------

    /// Advance the simulation by one time step.
    ///
    /// The whole collide -> stream -> moments update runs as a single fused
    /// pass over the grid ([`step_band`]); the open-edge BC pass and the
    /// boundary-line moment fix then touch only edge cells. Results are
    /// step-for-step identical to the classic separate-pass formulation
    /// (see `probe_state_hash`).
    pub fn step(&mut self) {
        if self.solid_dirty {
            self.rebuild_solid_runs();
        }
        let pf = self.fused_pass();
        self.capture_conv_stale();
        std::mem::swap(&mut self.f, &mut self.ftmp);
        std::mem::swap(&mut self.rho, &mut self.rho2);
        std::mem::swap(&mut self.ux, &mut self.ux2);
        std::mem::swap(&mut self.uy, &mut self.uy2);
        self.probed_force = pf;
        self.apply_open_edges();
        std::mem::swap(&mut self.conv_stale, &mut self.conv_stale_next);
        self.fix_boundary_moments();
        self.time += 1;
    }

    /// Advance the simulation by `steps` time steps.
    pub fn run(&mut self, steps: usize) {
        for _ in 0..steps {
            self.step();
        }
    }

    /// The fused collide+stream+moments pass: one sweep over the grid per
    /// step. Rows are processed in bands; each band collides its source rows
    /// just-in-time into a cache-resident ring, streams destinations from
    /// the ring, and writes moments into the double buffers (`rho2/ux2/uy2`)
    /// while the collide stage keeps reading the previous step's moments.
    fn fused_pass(&mut self) -> [T; 2] {
        let p = self.params();
        let g = self.geom();
        let (nx, ny) = (g.nx, g.ny);
        let planes = PlaneRows::new(&mut self.ftmp, nx, ny);
        let rho_new = FieldRows::new(&mut self.rho2, nx, ny);
        let ux_new = FieldRows::new(&mut self.ux2, nx, ny);
        let uy_new = FieldRows::new(&mut self.uy2, nx, ny);
        let (f, solid, wall_u) = (&self.f, &self.solid, &self.wall_u);
        let (rho_old, ux_old, uy_old) = (&self.rho, &self.ux, &self.uy);
        let runs = &self.solid_runs;
        let probe = self.probe.as_deref();
        let ff = self.force_field.as_deref();
        let force_on = p.fx != T::zero() || p.fy != T::zero() || ff.is_some();
        #[cfg(feature = "parallel")]
        let nbands = if self.use_parallel {
            rayon::current_num_threads().clamp(1, (ny / 16).max(1))
        } else {
            1
        };
        #[cfg(not(feature = "parallel"))]
        let nbands = 1;
        let band_size = ny.div_ceil(nbands);
        let body = |band: usize| -> [T; 2] {
            let y0 = band * band_size;
            let y1 = ((band + 1) * band_size).min(ny);
            step_band(
                y0,
                y1,
                f,
                &planes,
                rho_old,
                ux_old,
                uy_old,
                &rho_new,
                &ux_new,
                &uy_new,
                solid,
                runs,
                wall_u,
                ff,
                force_on,
                probe,
                &g,
                &p,
            )
        };
        let add = |a: [T; 2], b: [T; 2]| [a[0] + b[0], a[1] + b[1]];
        #[cfg(feature = "parallel")]
        if self.use_parallel {
            return (0..nbands)
                .into_par_iter()
                .map(body)
                .reduce(|| [T::zero(); 2], add);
        }
        (0..nbands).map(body).fold([T::zero(); 2], add)
    }

    /// Save the post-collide unknown populations of every `ConvectiveOutflow`
    /// edge cell into `conv_stale_next` — the BC's "previous step" memory
    /// term for the *next* step. The fused pass discards post-collide
    /// populations, so the affected edge cells (a handful per step) are
    /// re-collided here from the still-untouched `f` and previous moments;
    /// this reproduces the ring's values bit-for-bit.
    fn capture_conv_stale(&mut self) {
        let p = self.params();
        let (nx, ny) = (self.nx, self.ny);
        let n = nx * ny;
        let ff = self.force_field.as_deref();
        let force_on = p.fx != T::zero() || p.fy != T::zero() || ff.is_some();
        for edge in Edge::ALL {
            if self.conv_stale_next[edge.index()].is_none() {
                continue;
            }
            let (nxi, nyi) = edge.n_in();
            let (tx, ty) = (-nyi, nxi);
            let unknowns = [
                dir_index(nxi, nyi),
                dir_index(nxi + tx, nyi + ty),
                dir_index(nxi - tx, nyi - ty),
            ];
            let cells = self.side_cells(edge);
            let stale = self.conv_stale_next[edge.index()]
                .as_mut()
                .expect("checked above");
            for (coord, (x, y)) in cells.into_iter().enumerate() {
                let i = y * nx + x;
                if self.solid[i] {
                    continue;
                }
                let mut cell = [T::zero(); Q];
                for (q, v) in cell.iter_mut().enumerate() {
                    *v = self.f[q * n + i];
                }
                let mut it = cell.chunks_mut(1);
                let mut fq: [&mut [T]; Q] =
                    std::array::from_fn(|_| it.next().expect("Q chunks"));
                let (r1, ux1, uy1) = ([self.rho[i]], [self.ux[i]], [self.uy[i]]);
                match (force_on, ff) {
                    (true, Some(field)) => collide_span::<T, true, true>(
                        &mut fq,
                        &r1,
                        &ux1,
                        &uy1,
                        &field[i..i + 1],
                        0,
                        1,
                        &p,
                    ),
                    (true, None) => {
                        collide_span::<T, true, false>(&mut fq, &r1, &ux1, &uy1, &[], 0, 1, &p)
                    }
                    (false, _) => {
                        collide_span::<T, false, false>(&mut fq, &r1, &ux1, &uy1, &[], 0, 1, &p)
                    }
                }
                for (s, &q) in unknowns.iter().enumerate() {
                    stale[coord * 3 + s] = cell[q];
                }
            }
        }
    }

    /// Recompute rho/u on the four boundary lines from the post-BC
    /// populations. The fused in-pass moments saw the pre-BC (stale) unknown
    /// slots on open edges; every other boundary cell recomputes to the
    /// identical value (the fix is idempotent there).
    fn fix_boundary_moments(&mut self) {
        let p = self.params();
        let (nx, ny) = (self.nx, self.ny);
        let n = nx * ny;
        let half = T::r(0.5);
        let f = &self.f;
        let solid = &self.solid;
        let ff = self.force_field.as_deref();
        let rho = &mut self.rho;
        let ux = &mut self.ux;
        let uy = &mut self.uy;
        let mut fix = |i: usize| {
            if solid[i] {
                return;
            }
            let fi: [T; Q] = std::array::from_fn(|q| f[q * n + i]);
            let dr = fi[0]
                + ((fi[1] + fi[3]) + (fi[2] + fi[4]))
                + ((fi[5] + fi[7]) + (fi[6] + fi[8]));
            let a = fi[5] - fi[7];
            let b = fi[8] - fi[6];
            let mx = (fi[1] - fi[3]) + (a + b);
            let my = (fi[2] - fi[4]) + (a - b);
            let (fx, fy) = match ff {
                Some(field) => (p.fx + field[i][0], p.fy + field[i][1]),
                None => (p.fx, p.fy),
            };
            let r = T::one() + dr;
            rho[i] = r;
            let inv = T::one() / r;
            ux[i] = (mx + half * fx) * inv;
            uy[i] = (my + half * fy) * inv;
        };
        for x in 0..nx {
            fix(x);
            fix((ny - 1) * nx + x);
        }
        for y in 1..ny - 1 {
            fix(y * nx);
            fix(y * nx + nx - 1);
        }
    }

    fn update_moments(&mut self) {
        if self.solid_dirty {
            self.rebuild_solid_runs();
        }
        let p = self.params();
        let nx = self.nx;
        let n = self.nx * self.ny;
        let (f, runs) = (&self.f, &self.solid_runs);
        let ff = self.force_field.as_deref();
        let body = |(y, ((rrow, uxrow), uyrow)): (usize, ((&mut [T], &mut [T]), &mut [T]))| {
            let base = y * nx;
            let fq: [&[T]; Q] = std::array::from_fn(|q| &f[q * n + base..][..nx]);
            moments_row(
                &fq,
                rrow,
                uxrow,
                uyrow,
                ff.map(|f| &f[base..base + nx]),
                &runs[y],
                &p,
            );
        };
        #[cfg(feature = "parallel")]
        if self.use_parallel {
            self.rho
                .par_chunks_mut(nx)
                .zip(self.ux.par_chunks_mut(nx))
                .zip(self.uy.par_chunks_mut(nx))
                .enumerate()
                .for_each(body);
            return;
        }
        self.rho
            .chunks_mut(nx)
            .zip(self.ux.chunks_mut(nx))
            .zip(self.uy.chunks_mut(nx))
            .enumerate()
            .for_each(body);
    }

    // ------------------------------------------------------------------
    // Open-edge boundary conditions (Zou–He, outflow)
    // ------------------------------------------------------------------

    fn edge_bc(&self, edge: Edge) -> EdgeBC<T> {
        match edge {
            Edge::Left => self.edges.left,
            Edge::Right => self.edges.right,
            Edge::Bottom => self.edges.bottom,
            Edge::Top => self.edges.top,
        }
    }

    /// Cells on an edge, ordered by the along-edge coordinate (`y` for
    /// left/right, `x` for bottom/top).
    fn side_cells(&self, edge: Edge) -> Vec<(usize, usize)> {
        match edge {
            Edge::Left => (0..self.ny).map(|y| (0, y)).collect(),
            Edge::Right => (0..self.ny).map(|y| (self.nx - 1, y)).collect(),
            Edge::Bottom => (0..self.nx).map(|x| (x, 0)).collect(),
            Edge::Top => (0..self.nx).map(|x| (x, self.ny - 1)).collect(),
        }
    }

    fn apply_open_edges(&mut self) {
        for edge in Edge::ALL {
            match self.edge_bc(edge) {
                EdgeBC::VelocityInlet { u } => self.zou_he(edge, ZouHe::Velocity(u)),
                EdgeBC::PressureOutlet { rho } => self.zou_he(edge, ZouHe::Pressure(rho)),
                EdgeBC::Outflow => self.outflow(edge),
                EdgeBC::ConvectiveOutflow { u_conv } => self.convective_outflow(edge, u_conv),
                _ => {}
            }
        }
    }

    /// Convective outflow: `f(edge,t+1) = (f(edge,t) + Uc f(interior,t+1)))
    /// / (1 + Uc)`. The memory term `f(edge,t)` — the previous step's
    /// post-collide populations, which the pre-fusion implementation read
    /// from the stale unknown slots left by streaming — is now supplied by
    /// the explicitly captured `conv_stale` buffer (bit-identical values).
    fn convective_outflow(&mut self, edge: Edge, u_conv: T) {
        let (nxi, nyi) = edge.n_in();
        let (tx, ty) = (-nyi, nxi);
        let unknowns = [
            dir_index(nxi, nyi),
            dir_index(nxi + tx, nyi + ty),
            dir_index(nxi - tx, nyi - ty),
        ];
        let lam = u_conv;
        let inv = T::one() / (T::one() + lam);
        let nx = self.nx;
        let n = self.nx * self.ny;
        let stale = self.conv_stale[edge.index()]
            .take()
            .expect("conv_stale allocated for ConvectiveOutflow edges at build");
        // weight share for the mass correction over the 3 unknown links
        let wsum = T::r(W[unknowns[0]] + W[unknowns[1]] + W[unknowns[2]]);
        for (coord, (x, y)) in self.side_cells(edge).into_iter().enumerate() {
            let i = y * nx + x;
            let j = ((y as i32 + nyi) as usize) * nx + (x as i32 + nxi) as usize;
            if self.solid[i] || self.solid[j] {
                continue;
            }
            for (s, &q) in unknowns.iter().enumerate() {
                let prev = stale[coord * 3 + s];
                self.f[q * n + i] = (prev + lam * self.f[q * n + j]) * inv;
            }
            // Mass-consistency correction: without it the independent
            // population relaxation lets the edge density drift away and the
            // run eventually diverges. Pin rho(edge) to rho(neighbour) by
            // distributing the deficit over the unknowns by weight.
            let mut di = T::zero();
            let mut dj = T::zero();
            for q in 0..Q {
                di = di + self.f[q * n + i];
                dj = dj + self.f[q * n + j];
            }
            let corr = dj - di;
            for q in unknowns {
                self.f[q * n + i] = self.f[q * n + i] + corr * T::r(W[q]) / wsum;
            }
        }
        self.conv_stale[edge.index()] = Some(stale);
    }

    /// Zou–He boundary parameterised by the face normal.
    ///
    /// With inward normal `n` and tangent `t = rot90(n)`, the three unknown
    /// populations after streaming are `n`, `n+t`, `n-t`. Writing `S0` for the
    /// sum of edge-parallel populations, `S-` for the sum of outgoing ones and
    /// `T` for the tangential flux `f_{+t} - f_{-t}`:
    ///
    /// ```text
    /// rho     = (S0 + 2 S-) / (1 - u.n)          (velocity BC)
    /// u.n     = (S0 + 2 S-) / rho - 1            (pressure BC, u.t = 0)
    /// f_n     = f_-n     + (2/3) rho (u.n)
    /// f_{n±t} = f_{-n∓t} + (1/6) rho (u.n) ± [ (1/2) rho (u.t) - (1/2) T ]
    /// ```
    fn zou_he(&mut self, edge: Edge, kind: ZouHe<T>) {
        let (nxi, nyi) = edge.n_in();
        let (tx, ty) = (-nyi, nxi);
        let q_n = dir_index(nxi, nyi);
        let q_d1 = dir_index(nxi + tx, nyi + ty);
        let q_d2 = dir_index(nxi - tx, nyi - ty);
        let q_t = dir_index(tx, ty);
        let q_mt = dir_index(-tx, -ty);
        let (half, c23, c16, two) = (T::r(0.5), T::r(2.0 / 3.0), T::r(1.0 / 6.0), T::r(2.0));
        let (nxr, nyr) = (T::r(nxi as f64), T::r(nyi as f64));
        let (txr, tyr) = (T::r(tx as f64), T::r(ty as f64));
        let nx = self.nx;
        let n = self.nx * self.ny;
        let profile = self.inlet_profiles[edge.index()].take();
        for (coord, (x, y)) in self.side_cells(edge).into_iter().enumerate() {
            let i = y * nx + x;
            if self.solid[i] {
                continue;
            }
            let f = &mut self.f;
            // Deviation storage: the physical S0 + 2 S- equals the deviation
            // sums plus sum(w) over those directions, which is exactly 1 for
            // any straight edge (3 edge-parallel + 2x3 outgoing weights).
            let s0 = f[i] + f[q_t * n + i] + f[q_mt * n + i];
            let sneg = f[OPP[q_n] * n + i] + f[OPP[q_d1] * n + i] + f[OPP[q_d2] * n + i];
            let closure = s0 + two * sneg + T::one();
            let (r, un, ut) = match kind {
                ZouHe::Velocity(u) => {
                    let u = profile.as_ref().map_or(u, |p| p[coord]);
                    let un = u[0] * nxr + u[1] * nyr;
                    let ut = u[0] * txr + u[1] * tyr;
                    (closure / (T::one() - un), un, ut)
                }
                ZouHe::Pressure(rho_bc) => {
                    // From the closure rho (1 - u.n) = S0 + 2 S-.
                    let un = T::one() - closure / rho_bc;
                    (rho_bc, un, T::zero())
                }
            };
            let tcorr = half * (r * ut - (f[q_t * n + i] - f[q_mt * n + i]));
            f[q_n * n + i] = f[OPP[q_n] * n + i] + c23 * r * un;
            f[q_d1 * n + i] = f[OPP[q_d1] * n + i] + c16 * r * un + tcorr;
            f[q_d2 * n + i] = f[OPP[q_d2] * n + i] + c16 * r * un - tcorr;
        }
        self.inlet_profiles[edge.index()] = profile;
    }

    /// Zero-gradient outflow: copy the unknown populations from the cell one
    /// step inward along the face normal.
    fn outflow(&mut self, edge: Edge) {
        let (nxi, nyi) = edge.n_in();
        let (tx, ty) = (-nyi, nxi);
        let unknowns = [
            dir_index(nxi, nyi),
            dir_index(nxi + tx, nyi + ty),
            dir_index(nxi - tx, nyi - ty),
        ];
        let nx = self.nx;
        let n = self.nx * self.ny;
        for (x, y) in self.side_cells(edge) {
            let i = y * nx + x;
            let j = ((y as i32 + nyi) as usize) * nx + (x as i32 + nxi) as usize;
            if self.solid[i] || self.solid[j] {
                continue;
            }
            for q in unknowns {
                self.f[q * n + i] = self.f[q * n + j];
            }
        }
    }

    // ------------------------------------------------------------------
    // Setup helpers
    // ------------------------------------------------------------------

    fn on_open_edge(&self, x: usize, y: usize) -> bool {
        (x == 0 && self.edges.left.is_open())
            || (x == self.nx - 1 && self.edges.right.is_open())
            || (y == 0 && self.edges.bottom.is_open())
            || (y == self.ny - 1 && self.edges.top.is_open())
    }

    /// Mark a cell as solid (half-way bounce-back obstacle).
    ///
    /// Panics if `(x, y)` lies on an open (inlet/outlet/outflow) edge, which
    /// is unsupported.
    pub fn set_solid(&mut self, x: usize, y: usize) {
        assert!(x < self.nx && y < self.ny, "cell ({x},{y}) out of bounds");
        assert!(
            !self.on_open_edge(x, y),
            "cannot place solid cells on an open (inlet/outlet/outflow) edge"
        );
        let i = y * self.nx + x;
        self.solid[i] = true;
        self.solid_dirty = true;
        // Freeze the cell's moments in both double buffers: solid cells are
        // never rewritten by the fused pass, and multiphase wall adhesion
        // reads rho at solids (virtual wall density).
        self.rho2[i] = self.rho[i];
        self.ux2[i] = self.ux[i];
        self.uy2[i] = self.uy[i];
    }

    /// Mark every cell for which `pred(x, y)` returns true as solid.
    ///
    /// Panics under the same conditions as [`Simulation::set_solid`].
    pub fn set_solid_region(&mut self, pred: impl Fn(usize, usize) -> bool) {
        for y in 0..self.ny {
            for x in 0..self.nx {
                if pred(x, y) {
                    self.set_solid(x, y);
                }
            }
        }
    }

    /// Initialise every cell from the given `(rho, ux, uy)` field.
    ///
    /// Second-order consistent: sets `f = feq + f_neq`, where the
    /// Chapman–Enskog non-equilibrium part
    /// `f_neq = -w rho tau [3 (cc):∇u - div u]` is evaluated with central
    /// finite differences (periodic wrap on periodic axes, one-sided at
    /// walls). Pure-equilibrium initialisation would inject an O(1/N)
    /// error into smooth flows (measured on Taylor–Green; docs/PHYSICS.md).
    pub fn init_with(&mut self, init: impl Fn(usize, usize) -> (T, T, T)) {
        let p = self.params();
        let (nx, ny) = (self.nx, self.ny);
        let n = nx * ny;
        // Pass 1: store the macroscopic fields.
        for y in 0..ny {
            for x in 0..nx {
                let (r, vx, vy) = init(x, y);
                let i = y * nx + x;
                self.rho[i] = r;
                self.ux[i] = vx;
                self.uy[i] = vy;
            }
        }
        // Pass 2: f = feq + f_neq(∇u).
        let g = self.geom();
        let tau = T::r(self.tau());
        let three = T::r(3.0);
        for y in 0..ny {
            for x in 0..nx {
                let i = y * nx + x;
                let feq = equilibrium(&p, self.rho[i], self.ux[i], self.uy[i]);
                for q in 0..Q {
                    self.f[q * n + i] = feq[q];
                }
                if self.solid[i] {
                    continue;
                }
                // Central differences with graceful fallback to one-sided
                // when the neighbour is missing (wall rim / domain edge).
                let sample = |xx: isize, yy: isize| -> Option<usize> {
                    let mut xx = xx;
                    let mut yy = yy;
                    if xx < 0 || xx >= nx as isize {
                        if g.per_x {
                            xx = (xx + nx as isize) % nx as isize;
                        } else {
                            return None;
                        }
                    }
                    if yy < 0 || yy >= ny as isize {
                        if g.per_y {
                            yy = (yy + ny as isize) % ny as isize;
                        } else {
                            return None;
                        }
                    }
                    let j = yy as usize * nx + xx as usize;
                    (!self.solid[j]).then_some(j)
                };
                let diff = |plus: Option<usize>, minus: Option<usize>, field: &[T], own: T| -> T {
                    match (plus, minus) {
                        (Some(pj), Some(mj)) => (field[pj] - field[mj]) * T::r(0.5),
                        (Some(pj), None) => field[pj] - own,
                        (None, Some(mj)) => own - field[mj],
                        (None, None) => T::zero(),
                    }
                };
                let (xi, yi) = (x as isize, y as isize);
                let (xp, xm) = (sample(xi + 1, yi), sample(xi - 1, yi));
                let (yp, ym) = (sample(xi, yi + 1), sample(xi, yi - 1));
                let duxdx = diff(xp, xm, &self.ux, self.ux[i]);
                let duydx = diff(xp, xm, &self.uy, self.uy[i]);
                let duxdy = diff(yp, ym, &self.ux, self.ux[i]);
                let duydy = diff(yp, ym, &self.uy, self.uy[i]);
                let div = duxdx + duydy;
                for q in 0..Q {
                    let (cx, cy) = (p.cxr[q], p.cyr[q]);
                    let ccgu = cx * cx * duxdx + cx * cy * (duydx + duxdy) + cy * cy * duydy;
                    let fneq = -p.wr[q] * self.rho[i] * tau * (three * ccgu - div);
                    self.f[q * n + i] = self.f[q * n + i] + fneq;
                }
            }
        }
        self.update_moments();
        // Keep the double buffers coherent at cells the fused pass never
        // writes (solids, incl. init-prescribed virtual wall densities).
        self.rho2.copy_from_slice(&self.rho);
        self.ux2.copy_from_slice(&self.ux);
        self.uy2.copy_from_slice(&self.uy);
    }

    /// Prescribe a per-node velocity profile for a `VelocityInlet` edge,
    /// overriding its uniform velocity. `profile(c)` receives the along-edge
    /// coordinate (`y` for left/right edges, `x` for bottom/top) and returns
    /// `[ux, uy]`.
    ///
    /// Panics if the edge is not `VelocityInlet` or any speed exceeds
    /// [`MAX_SPEED`].
    pub fn set_inlet_profile(&mut self, edge: Edge, profile: impl Fn(usize) -> [T; 2]) {
        assert!(
            matches!(self.edge_bc(edge), EdgeBC::VelocityInlet { .. }),
            "set_inlet_profile: {edge:?} is not a VelocityInlet edge"
        );
        let len = match edge {
            Edge::Left | Edge::Right => self.ny,
            Edge::Bottom | Edge::Top => self.nx,
        };
        let values: Vec<[T; 2]> = (0..len).map(&profile).collect();
        for (c, u) in values.iter().enumerate() {
            let s = (u[0].as_f64().powi(2) + u[1].as_f64().powi(2)).sqrt();
            assert!(
                s <= MAX_SPEED,
                "inlet profile speed {s} at coordinate {c} exceeds the low-Mach limit {MAX_SPEED}"
            );
        }
        self.inlet_profiles[edge.index()] = Some(values);
    }

    /// Mutable access to the per-cell force field (`[fx, fy]` per cell,
    /// indexed `y*nx + x`), allocating it zero-filled on first use. The field
    /// is *added* to the uniform `force` and is intended to be rewritten each
    /// step by multiphase models (see `multiphase::ShanChen::update_force`).
    pub fn force_field_mut(&mut self) -> &mut [[T; 2]] {
        let n = self.nx * self.ny;
        self.force_field
            .get_or_insert_with(|| vec![[T::zero(); 2]; n])
    }

    /// Remove the per-cell force field (reverts to the uniform force only).
    pub fn clear_force_field(&mut self) {
        self.force_field = None;
    }

    /// Select the set of solid cells whose momentum-exchange force is
    /// accumulated each step (e.g. an obstacle for drag/lift measurement).
    pub fn set_force_probe(&mut self, pred: impl Fn(usize, usize) -> bool) {
        let mut mask = vec![false; self.nx * self.ny];
        for y in 0..self.ny {
            for x in 0..self.nx {
                mask[y * self.nx + x] = pred(x, y);
            }
        }
        self.probe = Some(mask);
    }

    /// Momentum-exchange force `[Fx, Fy]` on the probed solids during the
    /// most recent [`Simulation::step`]. Zero if no probe is set.
    pub fn probed_force(&self) -> [T; 2] {
        self.probed_force
    }

    // ------------------------------------------------------------------
    // Accessors
    // ------------------------------------------------------------------

    /// Lattice width in cells.
    pub fn nx(&self) -> usize {
        self.nx
    }
    /// Lattice height in cells.
    pub fn ny(&self) -> usize {
        self.ny
    }
    /// Number of completed time steps.
    pub fn time(&self) -> u64 {
        self.time
    }
    /// Kinematic viscosity (lattice units).
    pub fn nu(&self) -> f64 {
        self.nu
    }
    /// Relaxation time `tau = 3 nu + 0.5`.
    pub fn tau(&self) -> f64 {
        3.0 * self.nu + 0.5
    }

    #[inline]
    fn idx(&self, x: usize, y: usize) -> usize {
        debug_assert!(x < self.nx && y < self.ny);
        y * self.nx + x
    }

    /// Density at a cell.
    pub fn rho(&self, x: usize, y: usize) -> T {
        self.rho[self.idx(x, y)]
    }
    /// x-velocity at a cell (physical: includes the Guo half-force term).
    pub fn ux(&self, x: usize, y: usize) -> T {
        self.ux[self.idx(x, y)]
    }
    /// y-velocity at a cell (physical: includes the Guo half-force term).
    pub fn uy(&self, x: usize, y: usize) -> T {
        self.uy[self.idx(x, y)]
    }
    /// Whether a cell is solid.
    pub fn is_solid(&self, x: usize, y: usize) -> bool {
        self.solid[self.idx(x, y)]
    }
    /// Solid mask, indexed `[y*nx + x]`.
    pub fn solid_field(&self) -> &[bool] {
        &self.solid
    }
    /// Whether the x axis wraps periodically.
    pub fn is_periodic_x(&self) -> bool {
        self.edges.left.is_periodic()
    }
    /// Whether the y axis wraps periodically.
    pub fn is_periodic_y(&self) -> bool {
        self.edges.bottom.is_periodic()
    }
    /// Density field, indexed `[y*nx + x]`.
    pub fn rho_field(&self) -> &[T] {
        &self.rho
    }
    /// x-velocity field, indexed `[y*nx + x]`.
    pub fn ux_field(&self) -> &[T] {
        &self.ux
    }
    /// y-velocity field, indexed `[y*nx + x]`.
    pub fn uy_field(&self) -> &[T] {
        &self.uy
    }

    /// Number of fluid (non-solid) cells.
    pub fn fluid_cell_count(&self) -> usize {
        self.solid.iter().filter(|&&s| !s).count()
    }

    /// Total mass over fluid cells, computed directly from the populations.
    /// Accumulated in `f64` regardless of `T` so the diagnostic itself does
    /// not drown in summation round-off on `f32` grids.
    pub fn total_mass(&self) -> T {
        // Deviation storage: physical mass = fluid_cell_count + sum(f_dev).
        let n = self.nx * self.ny;
        let mut m = 0.0f64;
        let mut fluid = 0usize;
        for i in 0..n {
            if self.solid[i] {
                continue;
            }
            fluid += 1;
            for q in 0..Q {
                m += self.f[q * n + i].as_f64();
            }
        }
        T::r(fluid as f64 + m)
    }

    /// Total momentum `[sum rho ux, sum rho uy]` over fluid cells (physical,
    /// includes the half-force correction). Accumulated in `f64` like
    /// [`Simulation::total_mass`].
    pub fn total_momentum(&self) -> [T; 2] {
        let n = self.nx * self.ny;
        let mut px = 0.0f64;
        let mut py = 0.0f64;
        let (ufx, ufy) = (self.force[0].as_f64(), self.force[1].as_f64());
        let ff = self.force_field.as_deref();
        for i in 0..n {
            if self.solid[i] {
                continue;
            }
            let mut mx = 0.0f64;
            let mut my = 0.0f64;
            for q in 0..Q {
                let fq = self.f[q * n + i].as_f64();
                mx += CX[q] as f64 * fq;
                my += CY[q] as f64 * fq;
            }
            let (fx, fy) = match ff {
                Some(field) => (ufx + field[i][0].as_f64(), ufy + field[i][1].as_f64()),
                None => (ufx, ufy),
            };
            px += mx + 0.5 * fx;
            py += my + 0.5 * fy;
        }
        [T::r(px), T::r(py)]
    }
}

/// Equilibrium distribution in deviation form: `feq_q - w_q`.
///
/// Written in terms of `drho = rho - 1` so no large-magnitude cancellation
/// occurs — essential for the f32 precision of the deviation storage.
fn equilibrium<T: Real>(p: &Params<T>, r: T, vx: T, vy: T) -> [T; Q] {
    let three = T::r(3.0);
    let f45 = T::r(4.5);
    let f15 = T::r(1.5);
    let usq = vx * vx + vy * vy;
    let drho = r - T::one();
    let mut feq = [T::zero(); Q];
    for q in 0..Q {
        let cu = p.cxr[q] * vx + p.cyr[q] * vy;
        feq[q] = p.wr[q] * (drho + r * (three * cu + f45 * cu * cu - f15 * usq));
    }
    feq
}

/// TRT collision (BGK when `omega_m == omega_p`) with Guo forcing over the
/// solid-free spans of one row. `ff` optionally supplies a per-cell force
/// added to the uniform one; the option is dispatched here, outside the hot
/// loops, so both force flavours stay branch-free and vectorizable.
fn collide_row<T: Real, const FORCE: bool>(
    fq: &mut [&mut [T]; Q],
    rho: &[T],
    ux: &[T],
    uy: &[T],
    ff: Option<&[[T; 2]]>,
    runs: &[(u32, u32)],
    p: &Params<T>,
) {
    let nx = rho.len();
    let mut cursor = 0usize;
    for &(a, b) in runs.iter().chain(std::iter::once(&(nx as u32, nx as u32))) {
        match ff {
            Some(field) => {
                collide_span::<T, FORCE, true>(fq, rho, ux, uy, field, cursor, a as usize, p)
            }
            None => collide_span::<T, FORCE, false>(fq, rho, ux, uy, &[], cursor, a as usize, p),
        }
        cursor = b as usize;
    }
}

/// Branch-free TRT collision over `[x0, x1)` of one row (deviation form,
/// pairwise formulation). The pair decomposition works directly on
/// `ep = (feq_a + feq_b)/2` and `em = (feq_a - feq_b)/2`, halving the
/// equilibrium arithmetic and keeping x/y expressions mirror-symmetric so
/// lattice equivariance is preserved bit-for-bit.
///
/// `FORCE` compiles Guo forcing in or out; `FF` selects the per-cell force
/// field (`field`, ignored and empty when false) over the uniform force.
/// Both are const so every flavour of the inner loop is branch-free.
#[allow(clippy::too_many_arguments)]
fn collide_span<T: Real, const FORCE: bool, const FF: bool>(
    fq: &mut [&mut [T]; Q],
    rho: &[T],
    ux: &[T],
    uy: &[T],
    field: &[[T; 2]],
    x0: usize,
    x1: usize,
    p: &Params<T>,
) {
    if x0 >= x1 {
        return;
    }
    let w0 = T::r(4.0 / 9.0);
    let w1 = T::r(1.0 / 9.0);
    let w2 = T::r(1.0 / 36.0);
    let three = T::r(3.0);
    let f45 = T::r(4.5);
    let f15 = T::r(1.5);
    let nine = T::r(9.0);
    let half = T::r(0.5);
    let (op, om, cp, cm) = (p.omega_p, p.omega_m, p.cp, p.cm);
    let [f0, f1, f2, f3, f4, f5, f6, f7, f8] = fq;
    let (f0, f1, f2) = (&mut f0[x0..x1], &mut f1[x0..x1], &mut f2[x0..x1]);
    let (f3, f4, f5) = (&mut f3[x0..x1], &mut f4[x0..x1], &mut f5[x0..x1]);
    let (f6, f7, f8) = (&mut f6[x0..x1], &mut f7[x0..x1], &mut f8[x0..x1]);
    let rho = &rho[x0..x1];
    let ux = &ux[x0..x1];
    let uy = &uy[x0..x1];
    let field = if FF { &field[x0..x1] } else { field };
    let len = x1 - x0;
    assert!(
        f0.len() == len
            && f1.len() == len
            && f2.len() == len
            && f3.len() == len
            && f4.len() == len
            && f5.len() == len
            && f6.len() == len
            && f7.len() == len
            && f8.len() == len
            && rho.len() == len
            && ux.len() == len
            && uy.len() == len
            && (!FF || field.len() == len)
    );
    for x in 0..len {
        let (r, vx, vy) = (rho[x], ux[x], uy[x]);
        let usq = vx * vx + vy * vy;
        let drho = r - T::one();
        let base = drho - f15 * r * usq;
        let r3 = three * r;
        let r45 = f45 * r;
        // Forcing pieces (compiled out when FORCE is false).
        let (fx, fy, uf3) = if FORCE {
            let (fx, fy) = if FF {
                (p.fx + field[x][0], p.fy + field[x][1])
            } else {
                (p.fx, p.fy)
            };
            (fx, fy, three * (vx * fx + vy * fy))
        } else {
            (T::zero(), T::zero(), T::zero())
        };
        // Rest population: feq0 = w0 * base, src0 = -w0 * 3(u.F).
        let v0 = f0[x];
        f0[x] = if FORCE {
            v0 - op * (v0 - w0 * base) + cp * (-w0 * uf3)
        } else {
            v0 - op * (v0 - w0 * base)
        };
        // Axis pairs (w = 1/9): (1,3) along x, (2,4) along y.
        // Diagonal pairs (w = 1/36): (5,7) along +x+y, (6,8) along -x+y.
        macro_rules! pair {
            ($fa:ident, $fb:ident, $w:ident, $cu:expr, $cf:expr) => {{
                let cu = $cu;
                let (fa, fb) = ($fa[x], $fb[x]);
                let ep = $w * (base + r45 * cu * cu);
                let em = $w * (r3 * cu);
                let fp = half * (fa + fb);
                let fm = half * (fa - fb);
                let rp = op * (fp - ep);
                let rm = om * (fm - em);
                if FORCE {
                    let cf = $cf;
                    let sp = $w * (nine * cu * cf - uf3);
                    let sm = $w * (three * cf);
                    $fa[x] = fa - rp - rm + cp * sp + cm * sm;
                    $fb[x] = fb - rp + rm + cp * sp - cm * sm;
                } else {
                    $fa[x] = fa - rp - rm;
                    $fb[x] = fb - rp + rm;
                }
            }};
        }
        pair!(f1, f3, w1, vx, fx);
        pair!(f2, f4, w1, vy, fy);
        pair!(f5, f7, w2, vx + vy, fx + fy);
        pair!(f6, f8, w2, vy - vx, fy - fx);
    }
}

/// Collide global source row `s` into a ring slot, unless already resident.
/// The evicted slot is one whose tag is not in `needed` (`usize::MAX`
/// entries mark unused sources). The ring stores rows as `Q` consecutive
/// `nx`-slices per slot.
#[allow(clippy::too_many_arguments)]
fn ensure_collided<T: Real>(
    s: usize,
    needed: &[usize; 3],
    ring: &mut [T],
    tags: &mut [usize; 3],
    nx: usize,
    n: usize,
    f: &[T],
    rho_old: &[T],
    ux_old: &[T],
    uy_old: &[T],
    ff: Option<&[[T; 2]]>,
    force_on: bool,
    runs: &[Vec<(u32, u32)>],
    p: &Params<T>,
) {
    if tags.contains(&s) {
        return;
    }
    // A slot is free if it is empty or holds a row not needed for the
    // current destination row. (`usize::MAX` doubles as the "missing row"
    // marker in `needed`, so empty slots must be tested first.)
    let slot = (0..3)
        .find(|&k| tags[k] == usize::MAX || !needed.contains(&tags[k]))
        .expect("three ring slots cover at most two other needed rows");
    tags[slot] = s;
    let region = &mut ring[slot * Q * nx..(slot + 1) * Q * nx];
    let row = s * nx;
    let mut it = region.chunks_mut(nx);
    let mut fq: [&mut [T]; Q] = std::array::from_fn(|_| it.next().expect("Q chunks"));
    for (q, dst) in fq.iter_mut().enumerate() {
        dst.copy_from_slice(&f[q * n + row..][..nx]);
    }
    let (r, ux, uy) = (
        &rho_old[row..row + nx],
        &ux_old[row..row + nx],
        &uy_old[row..row + nx],
    );
    let ffrow = ff.map(|field| &field[row..row + nx]);
    if force_on {
        collide_row::<T, true>(&mut fq, r, ux, uy, ffrow, &runs[s], p);
    } else {
        collide_row::<T, false>(&mut fq, r, ux, uy, ffrow, &runs[s], p);
    }
}

/// One fused collide+stream+moments band: the complete time-step kernel for
/// destination rows `[y0, y1)` in a single sweep. This is the portable
/// reference kernel for the V2 `CpuSimd` backend — it reads `f` and the old
/// moments once, writes `out` and the new moments once, and keeps all
/// intermediate (post-collide) state in a 3-row cache-resident ring.
///
/// Source rows are collided just-in-time: streaming destination row `y`
/// pulls from the collided rows `y-1, y, y+1` held in the ring. The band's
/// outermost halo rows are collided redundantly by the neighbouring band
/// (same inputs, same results, no synchronisation). Solid-free source spans
/// stream as plain shifted copies; solid runs bounce back with momentum
/// injection for moving walls. Returns the momentum-exchange force over
/// probed solid links.
#[allow(clippy::too_many_arguments)]
fn step_band<T: Real>(
    y0: usize,
    y1: usize,
    f: &[T],
    out: &PlaneRows<'_, T>,
    rho_old: &[T],
    ux_old: &[T],
    uy_old: &[T],
    rho_new: &FieldRows<'_, T>,
    ux_new: &FieldRows<'_, T>,
    uy_new: &FieldRows<'_, T>,
    solid: &[bool],
    runs: &[Vec<(u32, u32)>],
    wall_u: &[[T; 2]],
    ff: Option<&[[T; 2]]>,
    force_on: bool,
    probe: Option<&[bool]>,
    g: &Geom,
    p: &Params<T>,
) -> [T; 2] {
    let six = T::r(6.0);
    let two = T::r(2.0);
    let nx = g.nx;
    let ny = g.ny;
    let n = nx * ny;
    let mut pf = [T::zero(); 2];
    if y0 >= y1 {
        return pf;
    }
    // Ring of collided source rows, tagged with global row indices.
    let mut ring = vec![T::zero(); 3 * Q * nx];
    let mut tags = [usize::MAX; 3];
    for y in y0..y1 {
        // Per-row partial sum, added to the band total at the end of the
        // row: keeps the same floating-point grouping as the pre-fusion
        // row-parallel reduction.
        let mut pf_row = [T::zero(); 2];
        // Global source rows feeding this destination row: sy = y - cy.
        let src_y = |cy: i32| -> usize {
            let sy = y as isize - cy as isize;
            if sy < 0 || sy >= ny as isize {
                if g.per_y {
                    ((sy + ny as isize) % ny as isize) as usize
                } else {
                    usize::MAX // open/wall edge: unknown slots stay stale
                }
            } else {
                sy as usize
            }
        };
        let needed = [src_y(1), src_y(0), src_y(-1)];
        for &s in &needed {
            if s != usize::MAX {
                ensure_collided(
                    s, &needed, &mut ring, &mut tags, nx, n, f, rho_old, ux_old, uy_old, ff,
                    force_on, runs, p,
                );
            }
        }
        let ring_row = |s: usize, q: usize| -> &[T] {
            let slot = tags.iter().position(|&t| t == s).expect("resident row");
            &ring[(slot * Q + q) * nx..][..nx]
        };
        for q in 0..Q {
            let sy = src_y(CY[q]);
            if sy == usize::MAX {
                continue;
            }
            let cx = CX[q] as isize;
            let src = ring_row(sy, q);
            // SAFETY: this band is the only writer of rows in [y0, y1).
            let dst = unsafe { out.row(q, y) };
            let mut cursor = 0usize;
            for &(a, b) in runs[sy]
                .iter()
                .chain(std::iter::once(&(nx as u32, nx as u32)))
            {
                let (a, b) = (a as usize, b as usize);
                // Fluid span [cursor, a): shifted copy ring -> dst.
                copy_span(dst, src, cursor, a, cx, g.per_x);
                // Solid run [a, b): half-way bounce-back into the
                // destination cells (reflected populations come from the
                // destination cell's own collided row, which is resident).
                for sx in a..b {
                    let mut x = sx as isize + cx;
                    if x < 0 || x >= nx as isize {
                        if g.per_x {
                            x = (x + nx as isize) % nx as isize;
                        } else {
                            continue;
                        }
                    }
                    let x = x as usize;
                    let i = y * nx + x;
                    if solid[i] {
                        continue;
                    }
                    let s = sy * nx + sx;
                    // In deviation storage the formula is unchanged
                    // (w_q == w_opp(q)).
                    let fout = ring_row(y, OPP[q])[x];
                    let wu = wall_u[s];
                    let cu = p.cxr[q] * wu[0] + p.cyr[q] * wu[1];
                    let fin = fout + six * p.wr[q] * rho_old[i] * cu;
                    dst[x] = fin;
                    if let Some(mask) = probe {
                        if mask[s] {
                            // Momentum given to the wall through this link,
                            // using physical populations (deviation + weight)
                            // so the static-pressure contribution on open
                            // surfaces (e.g. rims) is retained; for closed
                            // bodies the weight terms sum to exactly zero.
                            let ftot = fout + fin + two * p.wr[q];
                            pf_row[0] = pf_row[0] - p.cxr[q] * ftot;
                            pf_row[1] = pf_row[1] - p.cyr[q] * ftot;
                        }
                    }
                }
                cursor = b;
            }
        }
        // Fused moments: the nine destination rows just written are still
        // cache-resident. Written to the double buffers so in-flight collides
        // of other rows keep reading the previous step's moments.
        let fq: [&[T]; Q] = std::array::from_fn(|q| {
            // SAFETY: this band is the only holder of row `y`, and the
            // mutable borrows from the streaming loop above have ended.
            let row: &mut [T] = unsafe { out.row(q, y) };
            &*row
        });
        // SAFETY: rows [y0, y1) of the new-moment buffers belong to this band.
        let (rrow, uxrow, uyrow) =
            unsafe { (rho_new.row(y), ux_new.row(y), uy_new.row(y)) };
        moments_row(
            &fq,
            rrow,
            uxrow,
            uyrow,
            ff.map(|field| &field[y * nx..(y + 1) * nx]),
            &runs[y],
            p,
        );
        pf[0] = pf[0] + pf_row[0];
        pf[1] = pf[1] + pf_row[1];
    }
    pf
}

/// Copy the source span `[s0, s1)` of one plane row into its destination
/// (shifted by `cx`), wrapping the boundary element when periodic in x and
/// dropping it otherwise (open-edge slots keep their previous contents and
/// are rewritten by the BC pass).
#[inline]
fn copy_span<T: Copy>(dst: &mut [T], src: &[T], s0: usize, s1: usize, cx: isize, per_x: bool) {
    if s0 >= s1 {
        return;
    }
    let nx = src.len();
    match cx {
        0 => dst[s0..s1].copy_from_slice(&src[s0..s1]),
        1 => {
            let e = s1.min(nx - 1);
            if s0 < e {
                dst[s0 + 1..e + 1].copy_from_slice(&src[s0..e]);
            }
            if s1 == nx && per_x {
                dst[0] = src[nx - 1];
            }
        }
        -1 => {
            let s = s0.max(1);
            if s < s1 {
                dst[s - 1..s1 - 1].copy_from_slice(&src[s..s1]);
            }
            if s0 == 0 && per_x {
                dst[nx - 1] = src[0];
            }
        }
        _ => unreachable!("D2Q9 x-shifts are -1, 0, 1"),
    }
}

/// Recompute macroscopic fields from the populations over the solid-free
/// spans of one row. The per-cell force-field option is dispatched here,
/// outside the hot loops.
fn moments_row<T: Real>(
    fq: &[&[T]; Q],
    rho: &mut [T],
    ux: &mut [T],
    uy: &mut [T],
    ff: Option<&[[T; 2]]>,
    runs: &[(u32, u32)],
    p: &Params<T>,
) {
    let nx = rho.len();
    let mut cursor = 0usize;
    for &(a, b) in runs.iter().chain(std::iter::once(&(nx as u32, nx as u32))) {
        match ff {
            Some(field) => moments_span::<T, true>(fq, rho, ux, uy, field, cursor, a as usize, p),
            None => moments_span::<T, false>(fq, rho, ux, uy, &[], cursor, a as usize, p),
        }
        cursor = b as usize;
    }
}

/// Branch-free moment update over `[x0, x1)` of one row. The signed sums are
/// grouped pairwise (axis, then diagonals) to stay mirror-symmetric in x/y.
/// `FF` selects the per-cell force field (`field`, empty when false).
#[allow(clippy::too_many_arguments)]
fn moments_span<T: Real, const FF: bool>(
    fq: &[&[T]; Q],
    rho: &mut [T],
    ux: &mut [T],
    uy: &mut [T],
    field: &[[T; 2]],
    x0: usize,
    x1: usize,
    p: &Params<T>,
) {
    if x0 >= x1 {
        return;
    }
    let half = T::r(0.5);
    let [f0, f1, f2, f3, f4, f5, f6, f7, f8] = fq;
    let (f0, f1, f2) = (&f0[x0..x1], &f1[x0..x1], &f2[x0..x1]);
    let (f3, f4, f5) = (&f3[x0..x1], &f4[x0..x1], &f5[x0..x1]);
    let (f6, f7, f8) = (&f6[x0..x1], &f7[x0..x1], &f8[x0..x1]);
    let rho = &mut rho[x0..x1];
    let ux = &mut ux[x0..x1];
    let uy = &mut uy[x0..x1];
    let field = if FF { &field[x0..x1] } else { field };
    let len = x1 - x0;
    assert!(
        f0.len() == len
            && f1.len() == len
            && f2.len() == len
            && f3.len() == len
            && f4.len() == len
            && f5.len() == len
            && f6.len() == len
            && f7.len() == len
            && f8.len() == len
            && rho.len() == len
            && ux.len() == len
            && uy.len() == len
            && (!FF || field.len() == len)
    );
    for x in 0..len {
        // Deviation storage: rho = 1 + sum(f_dev); sum(w c) = 0 so the
        // momentum needs no correction.
        let dr = f0[x] + ((f1[x] + f3[x]) + (f2[x] + f4[x])) + ((f5[x] + f7[x]) + (f6[x] + f8[x]));
        let a = f5[x] - f7[x];
        let b = f8[x] - f6[x];
        let mx = (f1[x] - f3[x]) + (a + b);
        let my = (f2[x] - f4[x]) + (a - b);
        let (fx, fy) = if FF {
            (p.fx + field[x][0], p.fy + field[x][1])
        } else {
            (p.fx, p.fy)
        };
        let r = T::one() + dr;
        rho[x] = r;
        let inv = T::one() / r;
        ux[x] = (mx + half * fx) * inv;
        uy[x] = (my + half * fy) * inv;
    }
}

#[cfg(test)]
mod tests {
    use crate::domain::{Collision, EdgeBC, Edges, SimConfig};

    #[test]
    fn equilibrium_moments_are_exact() {
        // Deviation form: sum feq_dev = rho - 1, sum feq_dev c = rho u,
        // sum feq_dev cc = rho (cs2 I + u u) - cs2 I.
        use crate::lattice::{CS2, CX, CY, Q, W};
        let sim = SimConfig::<f64>::default().build().unwrap();
        let p = sim.params();
        for &(r, vx, vy) in &[(1.0, 0.0, 0.0), (0.9, 0.08, -0.05), (1.1, -0.1, 0.02)] {
            let feq: Vec<f64> = super::equilibrium(&p, r, vx, vy)
                .iter()
                .zip(W.iter())
                .map(|(dev, w)| dev + w)
                .collect();
            let m0: f64 = feq.iter().sum();
            let mut m1 = [0.0; 2];
            let mut m2 = [[0.0; 2]; 2];
            for q in 0..Q {
                let c = [CX[q] as f64, CY[q] as f64];
                m1[0] += feq[q] * c[0];
                m1[1] += feq[q] * c[1];
                for a in 0..2 {
                    for b in 0..2 {
                        m2[a][b] += feq[q] * c[a] * c[b];
                    }
                }
            }
            assert!((m0 - r).abs() < 1e-14, "mass: {m0} vs {r}");
            assert!((m1[0] - r * vx).abs() < 1e-14);
            assert!((m1[1] - r * vy).abs() < 1e-14);
            let u = [vx, vy];
            for a in 0..2 {
                for b in 0..2 {
                    let expect = r * (if a == b { CS2 } else { 0.0 } + u[a] * u[b]);
                    assert!(
                        (m2[a][b] - expect).abs() < 1e-14,
                        "second moment [{a}][{b}]: {} vs {expect}",
                        m2[a][b]
                    );
                }
            }
        }
    }

    #[test]
    fn pairwise_collision_matches_generic_reference() {
        // The span kernel's pair decomposition must agree with the plain
        // per-direction TRT formula to round-off accuracy.
        use crate::lattice::{PAIRS, Q};
        let sim = SimConfig::<f64> {
            nu: 0.03,
            ..Default::default()
        }
        .build()
        .unwrap();
        let p = sim.params();
        let (r, vx, vy) = (1.05, 0.07, -0.04);
        let (fx, fy) = (1e-4, -2e-4);
        let fin: [f64; Q] = std::array::from_fn(|q| 1e-3 * (q as f64 - 4.0));
        // Reference: generic per-direction TRT with Guo forcing.
        let feq = super::equilibrium(&p, r, vx, vy);
        let uf = vx * fx + vy * fy;
        let mut src = [0.0; Q];
        for q in 0..Q {
            let cu = p.cxr[q] * vx + p.cyr[q] * vy;
            let cf = p.cxr[q] * fx + p.cyr[q] * fy;
            src[q] = p.wr[q] * (3.0 * (cf - uf) + 9.0 * cu * cf);
        }
        let mut want = fin;
        want[0] = fin[0] - p.omega_p * (fin[0] - feq[0]) + p.cp * src[0];
        for (a, b) in PAIRS {
            let fp = 0.5 * (fin[a] + fin[b]);
            let fm = 0.5 * (fin[a] - fin[b]);
            let ep = 0.5 * (feq[a] + feq[b]);
            let em = 0.5 * (feq[a] - feq[b]);
            let sp = 0.5 * (src[a] + src[b]);
            let sm = 0.5 * (src[a] - src[b]);
            let rp = p.omega_p * (fp - ep);
            let rm = p.omega_m * (fm - em);
            want[a] = fin[a] - rp - rm + p.cp * sp + p.cm * sm;
            want[b] = fin[b] - rp + rm + p.cp * sp - p.cm * sm;
        }
        // Kernel under test, on a one-cell row.
        let mut cell: [Vec<f64>; Q] = std::array::from_fn(|q| vec![fin[q]]);
        let mut fq: [&mut [f64]; Q] = {
            let mut it = cell.iter_mut();
            std::array::from_fn(|_| it.next().unwrap().as_mut_slice())
        };
        let ffield = [[fx - p.fx, fy - p.fy]; 1];
        super::collide_span::<f64, true, true>(&mut fq, &[r], &[vx], &[vy], &ffield, 0, 1, &p);
        for q in 0..Q {
            assert!(
                (cell[q][0] - want[q]).abs() < 1e-15,
                "direction {q}: {} vs {}",
                cell[q][0],
                want[q]
            );
        }
    }

    #[test]
    #[should_panic(expected = "open")]
    fn solid_on_open_edge_panics() {
        let mut sim = SimConfig::<f64> {
            nx: 16,
            ny: 8,
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
        sim.set_solid(0, 3);
    }

    #[test]
    fn trt_reduces_to_bgk_when_omegas_match() {
        // TRT with magic chosen so omega_m == omega_p must equal BGK exactly.
        let nu = 0.1;
        let tau = 3.0 * nu + 0.5;
        let lam = tau - 0.5;
        let magic = lam * lam; // (1/w+ - 1/2)(1/w- - 1/2) with w- = w+ => lam^2
        let mk = |collision: Collision| {
            let mut s = SimConfig::<f64> {
                nx: 16,
                ny: 16,
                nu,
                collision,
                ..Default::default()
            }
            .build()
            .unwrap();
            s.init_with(|x, y| {
                let k = 2.0 * std::f64::consts::PI / 16.0;
                (
                    1.0,
                    0.02 * (k * y as f64).sin(),
                    0.02 * (k * x as f64).sin(),
                )
            });
            s.run(50);
            s
        };
        let a = mk(Collision::Bgk);
        let b = mk(Collision::Trt { magic });
        for i in 0..16 * 16 {
            assert!((a.ux_field()[i] - b.ux_field()[i]).abs() < 1e-14);
        }
    }
}
