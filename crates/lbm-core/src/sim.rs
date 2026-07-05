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

    /// All nine plane rows at `y`.
    ///
    /// # Safety
    /// As [`PlaneRows::row`], for all nine `(q, y)` rows at once.
    #[inline]
    unsafe fn rows(&self, y: usize) -> [&'a mut [T]; Q] {
        std::array::from_fn(|q| self.row(q, y))
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
        let mut sim = Self {
            nx: cfg.nx,
            ny: cfg.ny,
            ftmp: f.clone(),
            f,
            rho: vec![T::one(); n],
            ux: vec![T::zero(); n],
            uy: vec![T::zero(); n],
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
    pub fn step(&mut self) {
        if self.solid_dirty {
            self.rebuild_solid_runs();
        }
        self.collide();
        let pf = self.stream();
        std::mem::swap(&mut self.f, &mut self.ftmp);
        self.probed_force = pf;
        self.apply_open_edges();
        self.fix_boundary_moments();
        self.time += 1;
    }

    /// Advance the simulation by `steps` time steps.
    pub fn run(&mut self, steps: usize) {
        for _ in 0..steps {
            self.step();
        }
    }

    fn collide(&mut self) {
        let p = self.params();
        let (nx, ny) = (self.nx, self.ny);
        let planes = PlaneRows::new(&mut self.f, nx, ny);
        let (rho, ux, uy) = (&self.rho, &self.ux, &self.uy);
        let runs = &self.solid_runs;
        let ff = self.force_field.as_deref();
        let force_on = p.fx != T::zero() || p.fy != T::zero() || ff.is_some();
        let body = |y: usize| {
            let mut fq = unsafe { planes.rows(y) };
            let r = y * nx..(y + 1) * nx;
            let (rho, ux, uy) = (&rho[r.clone()], &ux[r.clone()], &uy[r.clone()]);
            let ff = ff.map(|f| &f[r]);
            if force_on {
                collide_row::<T, true>(&mut fq, rho, ux, uy, ff, &runs[y], &p);
            } else {
                collide_row::<T, false>(&mut fq, rho, ux, uy, ff, &runs[y], &p);
            }
        };
        #[cfg(feature = "parallel")]
        if self.use_parallel {
            (0..ny).into_par_iter().for_each(body);
            return;
        }
        (0..ny).for_each(body);
    }

    /// Pull-stream `f` into `ftmp` and recompute each row's macroscopic
    /// moments from the freshly streamed (cache-resident) populations in the
    /// same pass. Boundary-line moments are provisional at this point —
    /// open-edge cells still hold stale unknown slots — and are re-fixed
    /// after the BC pass by [`Simulation::fix_boundary_moments`].
    fn stream(&mut self) -> [T; 2] {
        let p = self.params();
        let g = self.geom();
        let nx = g.nx;
        let planes = PlaneRows::new(&mut self.ftmp, nx, g.ny);
        let rho_rows = FieldRows::new(&mut self.rho, nx, g.ny);
        let ux_rows = FieldRows::new(&mut self.ux, nx, g.ny);
        let uy_rows = FieldRows::new(&mut self.uy, nx, g.ny);
        let (f, solid, wall_u) = (&self.f, &self.solid, &self.wall_u);
        let runs = &self.solid_runs;
        let probe = self.probe.as_deref();
        let ff = self.force_field.as_deref();
        let body = |y: usize| -> [T; 2] {
            // SAFETY: each task owns exactly the rows of its own `y`.
            let (rrow, uxrow, uyrow) =
                unsafe { (rho_rows.row(y), ux_rows.row(y), uy_rows.row(y)) };
            stream_row(
                y,
                &planes,
                f,
                solid,
                runs,
                wall_u,
                rrow,
                uxrow,
                uyrow,
                ff.map(|field| &field[y * nx..(y + 1) * nx]),
                probe,
                &g,
                &p,
            )
        };
        let add = |a: [T; 2], b: [T; 2]| [a[0] + b[0], a[1] + b[1]];
        #[cfg(feature = "parallel")]
        if self.use_parallel {
            return (0..g.ny)
                .into_par_iter()
                .map(body)
                .reduce(|| [T::zero(); 2], add);
        }
        (0..g.ny).map(body).fold([T::zero(); 2], add)
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

    /// Convective outflow. In the pull scheme the unknown slots at the edge
    /// still hold the *previous step's* populations after streaming, so
    /// `f(edge,t+1) = (f(edge,t) + Uc f(interior,t+1)) / (1 + Uc)` needs no
    /// extra storage.
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
        // weight share for the mass correction over the 3 unknown links
        let wsum = T::r(W[unknowns[0]] + W[unknowns[1]] + W[unknowns[2]]);
        for (x, y) in self.side_cells(edge) {
            let i = y * nx + x;
            let j = ((y as i32 + nyi) as usize) * nx + (x as i32 + nxi) as usize;
            if self.solid[i] || self.solid[j] {
                continue;
            }
            for q in unknowns {
                let prev = self.f[q * n + i];
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
        self.solid[y * self.nx + x] = true;
        self.solid_dirty = true;
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
/// added to the uniform one.
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
        collide_span::<T, FORCE>(fq, rho, ux, uy, ff, cursor, a as usize, p);
        cursor = b as usize;
    }
}

/// Branch-free TRT collision over `[x0, x1)` of one row (deviation form,
/// pairwise formulation). The pair decomposition works directly on
/// `ep = (feq_a + feq_b)/2` and `em = (feq_a - feq_b)/2`, halving the
/// equilibrium arithmetic and keeping x/y expressions mirror-symmetric so
/// lattice equivariance is preserved bit-for-bit.
#[allow(clippy::too_many_arguments)]
fn collide_span<T: Real, const FORCE: bool>(
    fq: &mut [&mut [T]; Q],
    rho: &[T],
    ux: &[T],
    uy: &[T],
    ff: Option<&[[T; 2]]>,
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
    let ff = ff.map(|f| &f[x0..x1]);
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
            let (fx, fy) = match ff {
                Some(field) => (p.fx + field[x][0], p.fy + field[x][1]),
                None => (p.fx, p.fy),
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

/// Pull-scheme streaming for one destination row, source-side decomposed:
/// solid-free source spans become plain shifted copies, solid runs become
/// half-way bounce-back. Ends by recomputing the row's macroscopic moments
/// from the freshly streamed populations (`rho/ux/uy` are this row's slices:
/// read for the moving-wall term, then overwritten). Returns the
/// momentum-exchange force accumulated over probed solid links.
#[allow(clippy::too_many_arguments)]
fn stream_row<T: Real>(
    y: usize,
    out: &PlaneRows<'_, T>,
    f: &[T],
    solid: &[bool],
    solid_runs: &[Vec<(u32, u32)>],
    wall_u: &[[T; 2]],
    rho: &mut [T],
    ux: &mut [T],
    uy: &mut [T],
    ff: Option<&[[T; 2]]>,
    probe: Option<&[bool]>,
    g: &Geom,
    p: &Params<T>,
) -> [T; 2] {
    let six = T::r(6.0);
    let two = T::r(2.0);
    let nx = g.nx;
    let n = nx * g.ny;
    let mut pf = [T::zero(); 2];
    for q in 0..Q {
        let mut sy = y as isize - CY[q] as isize;
        if sy < 0 || sy >= g.ny as isize {
            if g.per_y {
                sy = (sy + g.ny as isize) % g.ny as isize;
            } else {
                // Unknown populations on an open edge; filled by the
                // open-edge pass right after streaming.
                continue;
            }
        }
        let sy = sy as usize;
        let cx = CX[q] as isize;
        let src = &f[q * n + sy * nx..][..nx];
        // SAFETY: this task is the only writer of row `y`.
        let dst = unsafe { out.row(q, y) };
        let mut cursor = 0usize;
        for &(a, b) in solid_runs[sy]
            .iter()
            .chain(std::iter::once(&(nx as u32, nx as u32)))
        {
            let (a, b) = (a as usize, b as usize);
            // Fluid span [cursor, a): shifted copy src -> dst.
            copy_span(dst, src, cursor, a, cx, g.per_x);
            // Solid run [a, b): half-way bounce-back into the destination
            // cells, with momentum injection for moving walls.
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
                let fout = f[OPP[q] * n + i];
                let wu = wall_u[s];
                let cu = p.cxr[q] * wu[0] + p.cyr[q] * wu[1];
                let fin = fout + six * p.wr[q] * rho[x] * cu;
                dst[x] = fin;
                if let Some(mask) = probe {
                    if mask[s] {
                        // Momentum given to the wall through this link, using
                        // physical populations (deviation + weight) so the
                        // static-pressure contribution on open surfaces (e.g.
                        // rims) is retained; for closed bodies the weight
                        // terms sum to exactly zero.
                        let ftot = fout + fin + two * p.wr[q];
                        pf[0] = pf[0] - p.cxr[q] * ftot;
                        pf[1] = pf[1] - p.cyr[q] * ftot;
                    }
                }
            }
            cursor = b;
        }
    }
    // Fused moments: the nine destination rows just written are still
    // cache-resident, so recomputing rho/u here saves the separate
    // full-field moments pass.
    let fq: [&[T]; Q] = std::array::from_fn(|q| {
        // SAFETY: this task is the only holder of row `y`, and the mutable
        // borrows from the streaming loop above have all ended.
        let row: &mut [T] = unsafe { out.row(q, y) };
        &*row
    });
    moments_row(&fq, rho, ux, uy, ff, &solid_runs[y], p);
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
/// spans of one row.
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
        moments_span(fq, rho, ux, uy, ff, cursor, a as usize, p);
        cursor = b as usize;
    }
}

/// Branch-free moment update over `[x0, x1)` of one row. The signed sums are
/// grouped pairwise (axis, then diagonals) to stay mirror-symmetric in x/y.
#[allow(clippy::too_many_arguments)]
fn moments_span<T: Real>(
    fq: &[&[T]; Q],
    rho: &mut [T],
    ux: &mut [T],
    uy: &mut [T],
    ff: Option<&[[T; 2]]>,
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
    let ff = ff.map(|f| &f[x0..x1]);
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
    );
    for x in 0..len {
        // Deviation storage: rho = 1 + sum(f_dev); sum(w c) = 0 so the
        // momentum needs no correction.
        let dr = f0[x] + ((f1[x] + f3[x]) + (f2[x] + f4[x])) + ((f5[x] + f7[x]) + (f6[x] + f8[x]));
        let a = f5[x] - f7[x];
        let b = f8[x] - f6[x];
        let mx = (f1[x] - f3[x]) + (a + b);
        let my = (f2[x] - f4[x]) + (a - b);
        let (fx, fy) = match ff {
            Some(field) => (p.fx + field[x][0], p.fy + field[x][1]),
            None => (p.fx, p.fy),
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
        super::collide_span::<f64, true>(
            &mut fq,
            &[r],
            &[vx],
            &[vy],
            Some(&ffield),
            0,
            1,
            &p,
        );
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
