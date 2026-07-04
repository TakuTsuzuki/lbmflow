//! The core simulation loop: collide → stream → open-edge BCs → moments.
//!
//! Design notes:
//! - Memory layout is cell-major AoS: `f[(y*nx + x)*Q + q]`.
//! - Streaming uses the pull scheme; wall edges are one-cell solid rims so
//!   half-way bounce-back handles them uniformly (no corner special cases).
//! - Macroscopic fields stored in `rho/ux/uy` always describe the *current*
//!   post-step state; velocities include the Guo half-force correction.

use crate::domain::{Collision, Edge, EdgeBC, Edges, SimConfig, MAX_SPEED};
use crate::lattice::{dir_index, CX, CY, OPP, PAIRS, Q, W};
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
            wall_u: vec![[T::zero(); 2]; n],
            edges: cfg.edges,
            omega_p,
            omega_m,
            nu: cfg.nu,
            force: cfg.force,
            probe: None,
            probed_force: [T::zero(); 2],
            inlet_profiles: [None, None, None, None],
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
        self.collide();
        let pf = self.stream();
        std::mem::swap(&mut self.f, &mut self.ftmp);
        self.probed_force = pf;
        self.apply_open_edges();
        self.update_moments();
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
        let nx = self.nx;
        let (rho, ux, uy, solid) = (&self.rho, &self.ux, &self.uy, &self.solid);
        let body = |(y, frow): (usize, &mut [T])| {
            let r = y * nx..(y + 1) * nx;
            collide_row(frow, &rho[r.clone()], &ux[r.clone()], &uy[r.clone()], &solid[r], &p);
        };
        #[cfg(feature = "parallel")]
        if self.use_parallel {
            self.f.par_chunks_mut(nx * Q).enumerate().for_each(body);
            return;
        }
        self.f.chunks_mut(nx * Q).enumerate().for_each(body);
    }

    fn stream(&mut self) -> [T; 2] {
        let p = self.params();
        let g = self.geom();
        let (f, solid, wall_u, rho) = (&self.f, &self.solid, &self.wall_u, &self.rho);
        let probe = self.probe.as_deref();
        let body = move |(y, row): (usize, &mut [T])| -> [T; 2] {
            stream_row(y, row, f, solid, wall_u, rho, probe, &g, &p)
        };
        let add = |a: [T; 2], b: [T; 2]| [a[0] + b[0], a[1] + b[1]];
        #[cfg(feature = "parallel")]
        if self.use_parallel {
            return self
                .ftmp
                .par_chunks_mut(g.nx * Q)
                .enumerate()
                .map(body)
                .reduce(|| [T::zero(); 2], add);
        }
        self.ftmp
            .chunks_mut(g.nx * Q)
            .enumerate()
            .map(body)
            .fold([T::zero(); 2], add)
    }

    fn update_moments(&mut self) {
        let p = self.params();
        let nx = self.nx;
        let (f, solid) = (&self.f, &self.solid);
        let body = |(y, ((rrow, uxrow), uyrow)): (usize, ((&mut [T], &mut [T]), &mut [T]))| {
            moments_row(
                &f[y * nx * Q..(y + 1) * nx * Q],
                rrow,
                uxrow,
                uyrow,
                &solid[y * nx..(y + 1) * nx],
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
                _ => {}
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
        let profile = self.inlet_profiles[edge.index()].take();
        for (coord, (x, y)) in self.side_cells(edge).into_iter().enumerate() {
            let i = y * nx + x;
            if self.solid[i] {
                continue;
            }
            let o = i * Q;
            let f = &mut self.f;
            // Deviation storage: the physical S0 + 2 S- equals the deviation
            // sums plus sum(w) over those directions, which is exactly 1 for
            // any straight edge (3 edge-parallel + 2x3 outgoing weights).
            let s0 = f[o] + f[o + q_t] + f[o + q_mt];
            let sneg = f[o + OPP[q_n]] + f[o + OPP[q_d1]] + f[o + OPP[q_d2]];
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
            let tcorr = half * (r * ut - (f[o + q_t] - f[o + q_mt]));
            f[o + q_n] = f[o + OPP[q_n]] + c23 * r * un;
            f[o + q_d1] = f[o + OPP[q_d1]] + c16 * r * un + tcorr;
            f[o + q_d2] = f[o + OPP[q_d2]] + c16 * r * un - tcorr;
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
        for (x, y) in self.side_cells(edge) {
            let i = y * nx + x;
            let j = ((y as i32 + nyi) as usize) * nx + (x as i32 + nxi) as usize;
            if self.solid[i] || self.solid[j] {
                continue;
            }
            for q in unknowns {
                self.f[i * Q + q] = self.f[j * Q + q];
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
                let o = i * Q;
                let feq = equilibrium(&p, self.rho[i], self.ux[i], self.uy[i]);
                self.f[o..o + Q].copy_from_slice(&feq);
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
                    let ccgu = cx * cx * duxdx
                        + cx * cy * (duydx + duxdy)
                        + cy * cy * duydy;
                    let fneq = -p.wr[q] * self.rho[i] * tau * (three * ccgu - div);
                    self.f[o + q] = self.f[o + q] + fneq;
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
        let mut m = 0.0f64;
        let mut fluid = 0usize;
        for i in 0..self.nx * self.ny {
            if self.solid[i] {
                continue;
            }
            fluid += 1;
            for q in 0..Q {
                m += self.f[i * Q + q].as_f64();
            }
        }
        T::r(fluid as f64 + m)
    }

    /// Total momentum `[sum rho ux, sum rho uy]` over fluid cells (physical,
    /// includes the half-force correction). Accumulated in `f64` like
    /// [`Simulation::total_mass`].
    pub fn total_momentum(&self) -> [T; 2] {
        let mut px = 0.0f64;
        let mut py = 0.0f64;
        let (fx, fy) = (self.force[0].as_f64(), self.force[1].as_f64());
        for i in 0..self.nx * self.ny {
            if self.solid[i] {
                continue;
            }
            let o = i * Q;
            let mut mx = 0.0f64;
            let mut my = 0.0f64;
            for q in 0..Q {
                let fq = self.f[o + q].as_f64();
                mx += CX[q] as f64 * fq;
                my += CY[q] as f64 * fq;
            }
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

/// TRT collision (BGK when `omega_m == omega_p`) with Guo forcing, one row.
fn collide_row<T: Real>(
    f: &mut [T],
    rho: &[T],
    ux: &[T],
    uy: &[T],
    solid: &[bool],
    p: &Params<T>,
) {
    let three = T::r(3.0);
    let f45 = T::r(4.5);
    let f15 = T::r(1.5);
    let nine = T::r(9.0);
    let half = T::r(0.5);
    let force_on = p.fx != T::zero() || p.fy != T::zero();
    for x in 0..rho.len() {
        if solid[x] {
            continue;
        }
        let o = x * Q;
        let (r, vx, vy) = (rho[x], ux[x], uy[x]);
        let usq = vx * vx + vy * vy;
        let uf = vx * p.fx + vy * p.fy;
        let drho = r - T::one();
        let mut feq = [T::zero(); Q]; // deviation form: feq_q - w_q
        let mut src = [T::zero(); Q];
        for q in 0..Q {
            let cu = p.cxr[q] * vx + p.cyr[q] * vy;
            feq[q] = p.wr[q] * (drho + r * (three * cu + f45 * cu * cu - f15 * usq));
            if force_on {
                let cf = p.cxr[q] * p.fx + p.cyr[q] * p.fy;
                src[q] = p.wr[q] * (three * (cf - uf) + nine * cu * cf);
            }
        }
        f[o] = f[o] - p.omega_p * (f[o] - feq[0]) + p.cp * src[0];
        for (a, b) in PAIRS {
            let (fa, fb) = (f[o + a], f[o + b]);
            let fp = half * (fa + fb);
            let fm = half * (fa - fb);
            let ep = half * (feq[a] + feq[b]);
            let em = half * (feq[a] - feq[b]);
            let sp = half * (src[a] + src[b]);
            let sm = half * (src[a] - src[b]);
            let rp = p.omega_p * (fp - ep);
            let rm = p.omega_m * (fm - em);
            f[o + a] = fa - rp - rm + p.cp * sp + p.cm * sm;
            f[o + b] = fb - rp + rm + p.cp * sp - p.cm * sm;
        }
    }
}

/// Pull-scheme streaming for one destination row. Returns the
/// momentum-exchange force accumulated over probed solid links.
#[allow(clippy::too_many_arguments)]
fn stream_row<T: Real>(
    y: usize,
    out: &mut [T],
    f: &[T],
    solid: &[bool],
    wall_u: &[[T; 2]],
    rho: &[T],
    probe: Option<&[bool]>,
    g: &Geom,
    p: &Params<T>,
) -> [T; 2] {
    let six = T::r(6.0);
    let mut pf = [T::zero(); 2];
    let nx = g.nx;
    for x in 0..nx {
        let i = y * nx + x;
        if solid[i] {
            continue;
        }
        let o = x * Q;
        for q in 0..Q {
            let mut sx = x as isize - CX[q] as isize;
            let mut sy = y as isize - CY[q] as isize;
            if sx < 0 || sx >= nx as isize {
                if g.per_x {
                    sx = (sx + nx as isize) % nx as isize;
                } else {
                    // Unknown population on an open edge; filled by the
                    // open-edge pass right after streaming.
                    continue;
                }
            }
            if sy < 0 || sy >= g.ny as isize {
                if g.per_y {
                    sy = (sy + g.ny as isize) % g.ny as isize;
                } else {
                    continue;
                }
            }
            let s = sy as usize * nx + sx as usize;
            if solid[s] {
                // Half-way bounce-back off the wall between cells s and i,
                // with momentum injection for moving walls. In deviation
                // storage the formula is unchanged (w_q == w_opp(q)).
                let fout = f[i * Q + OPP[q]];
                let wu = wall_u[s];
                let cu = p.cxr[q] * wu[0] + p.cyr[q] * wu[1];
                let fin = fout + six * p.wr[q] * rho[i] * cu;
                out[o + q] = fin;
                if let Some(mask) = probe {
                    if mask[s] {
                        // Momentum given to the wall through this link, using
                        // physical populations (deviation + weight) so the
                        // static-pressure contribution on open surfaces (e.g.
                        // rims) is retained; for closed bodies the weight
                        // terms sum to exactly zero.
                        let ftot = fout + fin + T::r(2.0) * p.wr[q];
                        pf[0] = pf[0] - p.cxr[q] * ftot;
                        pf[1] = pf[1] - p.cyr[q] * ftot;
                    }
                }
            } else {
                out[o + q] = f[s * Q + q];
            }
        }
    }
    pf
}

/// Recompute macroscopic fields from the populations for one row.
fn moments_row<T: Real>(
    f: &[T],
    rho: &mut [T],
    ux: &mut [T],
    uy: &mut [T],
    solid: &[bool],
    p: &Params<T>,
) {
    let half = T::r(0.5);
    for x in 0..rho.len() {
        if solid[x] {
            continue;
        }
        let o = x * Q;
        // Deviation storage: rho = 1 + sum(f_dev); sum(w c) = 0 so the
        // momentum needs no correction.
        let mut dr = T::zero();
        let mut mx = T::zero();
        let mut my = T::zero();
        for q in 0..Q {
            let fq = f[o + q];
            dr = dr + fq;
            mx = mx + p.cxr[q] * fq;
            my = my + p.cyr[q] * fq;
        }
        let r = T::one() + dr;
        rho[x] = r;
        let inv = T::one() / r;
        ux[x] = (mx + half * p.fx) * inv;
        uy[x] = (my + half * p.fy) * inv;
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
