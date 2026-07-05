//! V1 `lbm_core::sim::Simulation` facade over the V2 solver.
//!
//! Public API, semantics and panic messages replicate V1; the physics runs
//! on `Solver<D2Q9, T, CpuScalar, LocalPeriodic>` (bit-exact against V1,
//! `tests/v1_match.rs`).
//!
//! ```
//! use lbm_core::compat::prelude::*;
//!
//! // Lid-driven cavity, Re = U*L/nu = 0.1*62/0.02 = 310
//! let mut sim: Simulation<f64> = SimConfig {
//!     nx: 64,
//!     ny: 64,
//!     nu: 0.02,
//!     edges: Edges {
//!         left: EdgeBC::BounceBack,
//!         right: EdgeBC::BounceBack,
//!         bottom: EdgeBC::BounceBack,
//!         top: EdgeBC::MovingWall { u: [0.1, 0.0] },
//!     },
//!     ..Default::default()
//! }
//! .build()
//! .unwrap();
//! sim.run(100);
//! assert!(sim.ux(32, 60) != 0.0); // fluid is being dragged by the lid
//! ```

use super::domain::{Collision, Edge, EdgeBC, Edges, SimConfig, MAX_SPEED};
use super::lattice::{CX, CY, Q};
use super::real::Real;
use crate::backend::CpuScalar;
use crate::fields::SoaFields;
use crate::halo::LocalPeriodic;
use crate::lattice::D2Q9;
use crate::params::{CollisionKind, FaceBC};
use crate::solver::{build_wall_rims, GlobalSpec, Solver, WallSpec};

pub use crate::backend::PARALLEL_MIN_CELLS;

type CoreSolver<T> = Solver<D2Q9, T, CpuScalar, LocalPeriodic>;

/// D2Q9 lattice Boltzmann simulation on a rectangular grid (V1 facade).
pub struct Simulation<T: Real> {
    solver: CoreSolver<T>,
    edges: Edges<T>,
    force: [T; 2],
    /// Compact global solid mask (source of truth for `solid_field`; the
    /// solver's padded copies are synced through its mask API).
    solid_mirror: Vec<bool>,
    /// V1-shaped per-cell force field, staged into the core (as `[fx,fy,0]`)
    /// at the start of every step. Diagnostics read this live copy directly,
    /// exactly like V1.
    force2: Option<Vec<[T; 2]>>,
}

impl<T: Real> std::fmt::Debug for Simulation<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Simulation")
            .field("nx", &self.nx())
            .field("ny", &self.ny())
            .field("nu", &self.nu())
            .field("time", &self.time())
            .finish_non_exhaustive()
    }
}

impl<T: Real> Simulation<T> {
    pub(crate) fn from_config(cfg: SimConfig<T>) -> Self {
        let edge_faces = [
            (Edge::Left, cfg.edges.left),
            (Edge::Right, cfg.edges.right),
            (Edge::Bottom, cfg.edges.bottom),
            (Edge::Top, cfg.edges.top),
        ];
        let mut walls = WallSpec::<T>::default();
        let mut faces = [FaceBC::Closed; 6];
        let mut periodic = [false, false, false];
        for (edge, bc) in edge_faces {
            let face = edge.face();
            match bc {
                EdgeBC::Periodic => periodic[face.axis()] = true,
                EdgeBC::BounceBack => walls.is_wall[face.index()] = true,
                EdgeBC::MovingWall { u } => {
                    walls.is_wall[face.index()] = true;
                    walls.u[face.index()] = [u[0], u[1], T::zero()];
                }
                EdgeBC::VelocityInlet { u } => {
                    faces[face.index()] = FaceBC::Velocity {
                        u: [u[0], u[1], T::zero()],
                    }
                }
                EdgeBC::PressureOutlet { rho } => faces[face.index()] = FaceBC::Pressure { rho },
                EdgeBC::Outflow => faces[face.index()] = FaceBC::Outflow,
                EdgeBC::ConvectiveOutflow { u_conv } => {
                    faces[face.index()] = FaceBC::Convective { u_conv }
                }
            }
        }
        let spec = GlobalSpec {
            dims: [cfg.nx, cfg.ny, 1],
            nu: cfg.nu,
            collision: match cfg.collision {
                Collision::Bgk => CollisionKind::Bgk,
                Collision::Trt { magic } => CollisionKind::Trt { magic },
            },
            periodic,
            faces,
            force: [cfg.force[0], cfg.force[1], T::zero()],
        };
        let (solid, wall_u) = build_wall_rims(2, spec.dims, &walls);
        let solver = Solver::new(
            &spec,
            &solid,
            &wall_u,
            [1, 1, 1],
            CpuScalar::default(),
            LocalPeriodic,
        );
        Self {
            solver,
            edges: cfg.edges,
            force: cfg.force,
            solid_mirror: solid,
            force2: None,
        }
    }

    #[inline]
    fn fields(&self) -> &SoaFields<T> {
        self.solver.fields(0)
    }

    /// Stage the V1-shaped per-cell force into the core fields.
    fn sync_force_field(&mut self) {
        match &self.force2 {
            Some(f2) => {
                let core = self.solver.fields_mut(0);
                let ff = core
                    .force_field
                    .get_or_insert_with(|| vec![[T::zero(); 3]; f2.len()]);
                for (dst, src) in ff.iter_mut().zip(f2.iter()) {
                    *dst = [src[0], src[1], T::zero()];
                }
            }
            None => self.solver.fields_mut(0).force_field = None,
        }
    }

    // ------------------------------------------------------------------
    // Time stepping
    // ------------------------------------------------------------------

    /// Advance the simulation by one time step.
    pub fn step(&mut self) {
        self.sync_force_field();
        self.solver.step();
    }

    /// Advance the simulation by `steps` time steps.
    pub fn run(&mut self, steps: usize) {
        for _ in 0..steps {
            self.step();
        }
    }

    // ------------------------------------------------------------------
    // Setup helpers
    // ------------------------------------------------------------------

    fn on_open_edge(&self, x: usize, y: usize) -> bool {
        (x == 0 && self.edges.left.is_open())
            || (x == self.nx() - 1 && self.edges.right.is_open())
            || (y == 0 && self.edges.bottom.is_open())
            || (y == self.ny() - 1 && self.edges.top.is_open())
    }

    fn edge_bc(&self, edge: Edge) -> EdgeBC<T> {
        match edge {
            Edge::Left => self.edges.left,
            Edge::Right => self.edges.right,
            Edge::Bottom => self.edges.bottom,
            Edge::Top => self.edges.top,
        }
    }

    /// Whether `(x, y)` is the cell directly one step inward from an open
    /// edge — the neighbour an open-face BC reads to fill its unknown slots.
    /// A solid here makes the BC skip that column/row (`solid[edge] ||
    /// solid[interior]`), freezing the unknown populations at their initial
    /// values (A-3 / E5b: a stationary box with a right-Outflow pocket held a
    /// permanent ux = -0.115 at the edge, silent — every cell stays finite).
    fn is_open_edge_interior(&self, x: usize, y: usize) -> bool {
        let (nx, ny) = (self.nx(), self.ny());
        (x == 1 && self.edges.left.is_open())
            || (x == nx - 2 && self.edges.right.is_open())
            || (y == 1 && self.edges.bottom.is_open())
            || (y == ny - 2 && self.edges.top.is_open())
    }

    /// Whether [`Simulation::set_solid`] would accept `(x, y)`: it is in
    /// bounds, not on an open edge, and not the interior neighbour of one.
    /// Callers that must not panic (e.g. the GUI's paint tool) check this
    /// first.
    pub fn set_solid_allowed(&self, x: usize, y: usize) -> bool {
        x < self.nx()
            && y < self.ny()
            && !self.on_open_edge(x, y)
            && !self.is_open_edge_interior(x, y)
    }

    /// Mark a cell as solid (half-way bounce-back obstacle).
    ///
    /// # Panics
    /// Panics if `(x, y)` lies on an open (inlet/outlet/outflow) edge, or is
    /// the cell directly inward from one: an open-face BC reads that interior
    /// neighbour to reconstruct its unknown populations, and a solid there
    /// makes it silently skip, freezing those populations (A-3).
    pub fn set_solid(&mut self, x: usize, y: usize) {
        assert!(
            x < self.nx() && y < self.ny(),
            "cell ({x},{y}) out of bounds"
        );
        assert!(
            !self.on_open_edge(x, y),
            "cannot place solid cells on an open (inlet/outlet/outflow) edge"
        );
        assert!(
            !self.is_open_edge_interior(x, y),
            "cannot place a solid cell ({x},{y}) directly inward from an open \
             (inlet/outlet/outflow) edge: the boundary condition reads this \
             neighbour and would silently freeze its unknown populations"
        );
        let nx = self.nx();
        self.solid_mirror[y * nx + x] = true;
        self.solver.set_solid(x, y, 0);
    }

    /// Mark every cell for which `pred(x, y)` returns true as solid.
    ///
    /// Panics under the same conditions as [`Simulation::set_solid`].
    pub fn set_solid_region(&mut self, pred: impl Fn(usize, usize) -> bool) {
        for y in 0..self.ny() {
            for x in 0..self.nx() {
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
    ///
    /// # Panics
    /// Panics (with the offending coordinate) if `init` returns a density
    /// that is not strictly positive and finite — `rho = 0` makes the
    /// equilibrium `0 × ∞ = NaN` immediately — or a speed exceeding
    /// [`MAX_SPEED`]. `init` must also be *pure*: the finite-difference
    /// stencil re-evaluates it at each cell's neighbours, so it is called up
    /// to five times per cell.
    pub fn init_with(&mut self, init: impl Fn(usize, usize) -> (T, T, T)) {
        self.sync_force_field();
        self.solver.init_with(|x, y, _| {
            let (r, vx, vy) = init(x, y);
            assert!(
                r.as_f64() > 0.0 && r.as_f64().is_finite(),
                "init_with: density at ({x},{y}) must be > 0 and finite, got {}",
                r.as_f64()
            );
            let s = (vx.as_f64().powi(2) + vy.as_f64().powi(2)).sqrt();
            assert!(
                s <= MAX_SPEED,
                "init_with: speed {s} at ({x},{y}) exceeds the low-Mach limit {MAX_SPEED}"
            );
            (r, [vx, vy, T::zero()])
        });
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
            Edge::Left | Edge::Right => self.ny(),
            Edge::Bottom | Edge::Top => self.nx(),
        };
        let values: Vec<[T; 2]> = (0..len).map(&profile).collect();
        for (c, u) in values.iter().enumerate() {
            let s = (u[0].as_f64().powi(2) + u[1].as_f64().powi(2)).sqrt();
            assert!(
                s <= MAX_SPEED,
                "inlet profile speed {s} at coordinate {c} exceeds the low-Mach limit {MAX_SPEED}"
            );
        }
        let v3: Vec<[T; 3]> = values.iter().map(|u| [u[0], u[1], T::zero()]).collect();
        self.solver.set_inlet_profile(edge.face(), &v3);
    }

    /// Mutable access to the per-cell force field (`[fx, fy]` per cell,
    /// indexed `y*nx + x`), allocating it zero-filled on first use. The field
    /// is *added* to the uniform `force` and is intended to be rewritten each
    /// step by multiphase models (see `multiphase::ShanChen::update_force`).
    pub fn force_field_mut(&mut self) -> &mut [[T; 2]] {
        let n = self.nx() * self.ny();
        self.force2.get_or_insert_with(|| vec![[T::zero(); 2]; n])
    }

    /// Remove the per-cell force field (reverts to the uniform force only).
    pub fn clear_force_field(&mut self) {
        self.force2 = None;
        self.solver.fields_mut(0).force_field = None;
    }

    /// Select the set of solid cells whose momentum-exchange force is
    /// accumulated each step (e.g. an obstacle for drag/lift measurement).
    pub fn set_force_probe(&mut self, pred: impl Fn(usize, usize) -> bool) {
        self.solver.set_force_probe(|x, y, _| pred(x, y));
    }

    /// Momentum-exchange force `[Fx, Fy]` on the probed solids during the
    /// most recent [`Simulation::step`]. Zero if no probe is set.
    pub fn probed_force(&self) -> [T; 2] {
        let f = self.solver.probed_force();
        [f[0], f[1]]
    }

    // ------------------------------------------------------------------
    // Accessors
    // ------------------------------------------------------------------

    /// Lattice width in cells.
    pub fn nx(&self) -> usize {
        self.solver.dims()[0]
    }
    /// Lattice height in cells.
    pub fn ny(&self) -> usize {
        self.solver.dims()[1]
    }
    /// Number of completed time steps.
    pub fn time(&self) -> u64 {
        self.solver.time()
    }
    /// Kinematic viscosity (lattice units).
    pub fn nu(&self) -> f64 {
        self.solver.nu()
    }
    /// Relaxation time `tau = 3 nu + 0.5`.
    pub fn tau(&self) -> f64 {
        self.solver.tau()
    }

    #[inline]
    fn idx(&self, x: usize, y: usize) -> usize {
        debug_assert!(x < self.nx() && y < self.ny());
        y * self.nx() + x
    }

    /// Density at a cell.
    pub fn rho(&self, x: usize, y: usize) -> T {
        self.fields().rho[self.idx(x, y)]
    }
    /// x-velocity at a cell (physical: includes the Guo half-force term).
    pub fn ux(&self, x: usize, y: usize) -> T {
        self.fields().ux[self.idx(x, y)]
    }
    /// y-velocity at a cell (physical: includes the Guo half-force term).
    pub fn uy(&self, x: usize, y: usize) -> T {
        self.fields().uy[self.idx(x, y)]
    }
    /// Whether a cell is solid.
    pub fn is_solid(&self, x: usize, y: usize) -> bool {
        self.solid_mirror[self.idx(x, y)]
    }
    /// Solid mask, indexed `[y*nx + x]`.
    pub fn solid_field(&self) -> &[bool] {
        &self.solid_mirror
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
        &self.fields().rho
    }
    /// x-velocity field, indexed `[y*nx + x]`.
    pub fn ux_field(&self) -> &[T] {
        &self.fields().ux
    }
    /// y-velocity field, indexed `[y*nx + x]`.
    pub fn uy_field(&self) -> &[T] {
        &self.fields().uy
    }

    /// Number of fluid (non-solid) cells.
    pub fn fluid_cell_count(&self) -> usize {
        self.solid_mirror.iter().filter(|&&s| !s).count()
    }

    /// Total mass over fluid cells, computed directly from the populations.
    /// Accumulated in `f64` regardless of `T` so the diagnostic itself does
    /// not drown in summation round-off on `f32` grids.
    pub fn total_mass(&self) -> T {
        self.solver.total_mass()
    }

    /// Total momentum `[sum rho ux, sum rho uy]` over fluid cells (physical,
    /// includes the half-force correction). Accumulated in `f64` like
    /// [`Simulation::total_mass`].
    ///
    /// Computed against the live V1-shaped force field (not the staged core
    /// copy), so mutations between steps are visible exactly as in V1.
    pub fn total_momentum(&self) -> [T; 2] {
        let f = self.fields();
        let g = f.geom;
        let np = g.n_padded();
        let (nx, ny) = (self.nx(), self.ny());
        let mut px = 0.0f64;
        let mut py = 0.0f64;
        let (ufx, ufy) = (self.force[0].as_f64(), self.force[1].as_f64());
        let ff = self.force2.as_deref();
        for y in 0..ny {
            for x in 0..nx {
                let i = y * nx + x;
                let pi = g.pidx(x, y, 0);
                if f.solid[pi] {
                    continue;
                }
                let mut mx = 0.0f64;
                let mut my = 0.0f64;
                for q in 0..Q {
                    let fq = f.f[q * np + pi].as_f64();
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
        }
        [T::r(px), T::r(py)]
    }
}

#[cfg(test)]
mod tests {
    use crate::compat::domain::{Collision, EdgeBC, Edges, SimConfig};
    use crate::compat::sim::Simulation;

    /// A-7: `init_with` rejects a zero density (immediate `0 × ∞` NaN) with a
    /// coordinate-bearing panic.
    #[test]
    #[should_panic(expected = "density at (3,3)")]
    fn init_with_rejects_zero_density() {
        let mut sim = SimConfig::<f64> {
            nx: 8,
            ny: 8,
            ..Default::default()
        }
        .build()
        .unwrap();
        sim.init_with(|x, y| {
            if (x, y) == (3, 3) {
                (0.0, 0.0, 0.0)
            } else {
                (1.0, 0.0, 0.0)
            }
        });
    }

    /// A-7: `init_with` rejects a non-finite density.
    #[test]
    #[should_panic(expected = "must be > 0 and finite")]
    fn init_with_rejects_nan_density() {
        let mut sim = SimConfig::<f64> {
            nx: 8,
            ny: 8,
            ..Default::default()
        }
        .build()
        .unwrap();
        sim.init_with(|_, _| (f64::NAN, 0.0, 0.0));
    }

    /// A-7: `init_with` rejects a super-sonic seed velocity (asymmetry with
    /// `set_inlet_profile`, which already checked, is closed).
    #[test]
    #[should_panic(expected = "exceeds the low-Mach limit")]
    fn init_with_rejects_too_fast() {
        let mut sim = SimConfig::<f64> {
            nx: 8,
            ny: 8,
            ..Default::default()
        }
        .build()
        .unwrap();
        sim.init_with(|_, _| (1.0, 0.9, 0.0));
    }

    /// A-7: a legal seed field initialises without panicking.
    #[test]
    fn init_with_accepts_valid_field() {
        let mut sim = SimConfig::<f64> {
            nx: 16,
            ny: 16,
            ..Default::default()
        }
        .build()
        .unwrap();
        sim.init_with(|x, _| (1.0 + 0.01 * x as f64, 0.02, 0.0));
        assert!(sim.rho(4, 4).is_finite());
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

    fn e5b_sim() -> Simulation<f64> {
        // E5b geometry: inlet left, Outflow right, walls top/bottom.
        SimConfig::<f64> {
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
        .unwrap()
    }

    /// A-3 (E5b): a solid on the cell directly inward from the right Outflow
    /// (x = nx-2) freezes the outflow's unknown populations (permanent
    /// non-physical ux, no NaN). It must be rejected.
    #[test]
    #[should_panic(expected = "directly inward from an open")]
    fn solid_inward_of_outflow_panics() {
        let mut sim = e5b_sim();
        sim.set_solid(14, 4); // nx-2 = 14, the pocket cell that broke E5b
    }

    /// A-3: the cell inward from the left inlet (x = 1) is likewise rejected.
    #[test]
    #[should_panic(expected = "directly inward from an open")]
    fn solid_inward_of_inlet_panics() {
        let mut sim = e5b_sim();
        sim.set_solid(1, 4);
    }

    /// A-3: `set_solid_allowed` mirrors the panic condition; a solid two cells
    /// inward from the outflow (the E5 "plug 1 cell further in" control) is
    /// legal, and placing it leaves the field finite (the physical case).
    #[test]
    fn set_solid_allowed_matches_and_legal_shape_runs() {
        let mut sim = e5b_sim();
        let (nx, ny) = (sim.nx(), sim.ny());
        // Open edges + their interior neighbours are disallowed.
        assert!(!sim.set_solid_allowed(0, 4)); // on the inlet edge
        assert!(!sim.set_solid_allowed(nx - 1, 4)); // on the outflow edge
        assert!(!sim.set_solid_allowed(1, 4)); // inward of inlet
        assert!(!sim.set_solid_allowed(nx - 2, 4)); // inward of outflow
                                                    // A cell two in from the outflow is fine.
        assert!(sim.set_solid_allowed(nx - 3, 4));
        sim.set_solid(nx - 3, 4);
        assert!(sim.is_solid(nx - 3, 4));
        sim.run(50);
        assert!(sim.rho_field().iter().all(|r| r.is_finite()));
        // y-open variant: inlet on the bottom edge disallows y == 1.
        let sim_y = SimConfig::<f64> {
            nx: 8,
            ny: 16,
            edges: Edges {
                left: EdgeBC::BounceBack,
                right: EdgeBC::BounceBack,
                bottom: EdgeBC::VelocityInlet { u: [0.0, 0.05] },
                top: EdgeBC::Outflow,
            },
            ..Default::default()
        }
        .build()
        .unwrap();
        let _ = ny;
        assert!(!sim_y.set_solid_allowed(4, 1)); // inward of bottom inlet
        assert!(!sim_y.set_solid_allowed(4, sim_y.ny() - 2)); // inward of top outflow
        assert!(sim_y.set_solid_allowed(4, 2)); // two in from bottom is fine
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
