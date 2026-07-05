//! Solver orchestrator: drives subdomains × backend × halo exchange through
//! the V1 step sequence (collide → stream → open faces → moments).
//!
//! The step-phase order, the diagnostics' f64 accumulation and the
//! initialisation paths reproduce V1 `Simulation` mechanics exactly; the
//! compat facade is a thin wrapper over this type with a monolithic (1×1×1)
//! decomposition.

use crate::backend::{Backend, HostMoments};
use crate::fields::SoaFields;
use crate::halo::{ExchangeScope, HaloExchange};
use crate::kernels::equilibrium;
use crate::lattice::{Face, Lattice};
use crate::params::{CollisionKind, FaceBC, KParams, Reduction, StepParams, MAX_SPEED};
use crate::real::Real;
use crate::subdomain::Subdomain;

/// Global scenario description, backend/decomposition agnostic.
#[derive(Clone, Debug)]
pub struct GlobalSpec<T: Real> {
    /// Global grid extents `[nx, ny, nz]` (`nz == 1` for 2D lattices).
    pub dims: [usize; 3],
    /// Kinematic viscosity (lattice units); `tau = 3 nu + 0.5`.
    pub nu: f64,
    /// Collision operator.
    pub collision: CollisionKind,
    /// Periodic wrap per axis.
    pub periodic: [bool; 3],
    /// Open BC per global face (`Closed` for periodic/wall faces).
    pub faces: [FaceBC<T>; 6],
    /// Uniform body force (Guo forcing).
    pub force: [T; 3],
}

impl<T: Real> Default for GlobalSpec<T> {
    fn default() -> Self {
        Self {
            dims: [64, 64, 1],
            nu: 1.0 / 6.0,
            collision: CollisionKind::Trt {
                magic: CollisionKind::MAGIC_STD,
            },
            periodic: [true, true, false],
            faces: [FaceBC::Closed; 6],
            force: [T::zero(); 3],
        }
    }
}

/// A rejected [`GlobalSpec`] (A-4): the V2-native counterpart of the compat
/// facade's `ConfigError`. `Solver::build` calls [`GlobalSpec::validate`]
/// before allocating, turning the previously-silent non-physical
/// configurations (stale data on an uncovered face, ν = 0, periodic × open
/// on one axis, …) into hard errors.
#[derive(Clone, Debug, PartialEq)]
pub enum SpecError {
    /// Kinematic viscosity must be finite and > 0 (`tau = 3ν + 0.5 > 0.5`);
    /// `ν = 0` gives `omega_m = 0` and a non-physical relaxation (E3).
    NonPositiveViscosity {
        /// Offending value.
        nu: f64,
    },
    /// A parameter that must be finite is NaN or infinite.
    NonFiniteParameter {
        /// Parameter name.
        what: &'static str,
    },
    /// TRT magic Λ must be finite and > 0.
    InvalidMagic {
        /// Offending value.
        magic: f64,
    },
    /// The domain must be at least 3 cells on every active axis.
    DomainTooSmall {
        /// Configured extents.
        dims: [usize; 3],
    },
    /// `periodic` must not be combined with an open BC on the same axis.
    PeriodicOpenConflict {
        /// Axis (0 = x, 1 = y, 2 = z).
        axis: usize,
    },
    /// Open faces may lie on at most one axis (a shared domain edge breaks the
    /// Zou–He face assumptions — the 3D lift of V1's corner rule).
    OpenFacesOnMultipleAxes,
    /// A non-periodic face is neither an open BC nor fully covered by a solid
    /// wall rim: its halo would feed stale data into the interior every step
    /// (E2 — silent non-physical drift, no NaN).
    UncoveredFace {
        /// The offending face index ([`Face::index`]).
        face: usize,
    },
    /// A prescribed velocity (inlet or z-normal component etc.) exceeds
    /// [`MAX_SPEED`] (NaN-safe: NaN is rejected here too).
    VelocityTooHigh {
        /// Offending speed magnitude.
        speed: f64,
    },
    /// A prescribed outlet density must be finite and > 0.
    NonPositiveDensity {
        /// Offending value.
        rho: f64,
    },
    /// A convective-outflow advection speed must lie in `(0, 1]`.
    InvalidConvectiveSpeed {
        /// Offending value.
        u_conv: f64,
    },
    /// A 2D lattice must have a zero z body-force component.
    NonZeroZForce2D {
        /// Offending value.
        fz: f64,
    },
    /// An open face's own axis must span at least 3 cells (the Zou–He /
    /// outflow stencil reads one cell inward).
    OpenFaceAxisTooShort {
        /// The offending face index.
        face: usize,
        /// The axis extent.
        extent: usize,
    },
}

impl std::fmt::Display for SpecError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SpecError::NonPositiveViscosity { nu } => write!(
                f,
                "kinematic viscosity must be > 0 (got {nu}); tau = 3*nu + 0.5 must exceed 0.5"
            ),
            SpecError::NonFiniteParameter { what } => {
                write!(f, "parameter {what} must be finite (got NaN or infinity)")
            }
            SpecError::InvalidMagic { magic } => {
                write!(f, "TRT magic must be finite and > 0 (got {magic})")
            }
            SpecError::DomainTooSmall { dims } => write!(
                f,
                "domain must be at least 3 cells on every active axis (got {dims:?})"
            ),
            SpecError::PeriodicOpenConflict { axis } => write!(
                f,
                "axis {axis} is periodic and also carries an open BC; a face is one or the other"
            ),
            SpecError::OpenFacesOnMultipleAxes => write!(
                f,
                "open faces (inlet/outlet/outflow) may lie on at most one axis; \
                 perpendicular faces must be walls or periodic"
            ),
            SpecError::UncoveredFace { face } => write!(
                f,
                "face {face} is neither periodic, an open BC, nor a full solid wall rim; \
                 its halo would feed stale values into the interior every step"
            ),
            SpecError::VelocityTooHigh { speed } => write!(
                f,
                "prescribed speed {speed} exceeds the low-Mach limit {MAX_SPEED} (lattice units)"
            ),
            SpecError::NonPositiveDensity { rho } => {
                write!(f, "prescribed density must be > 0 (got {rho})")
            }
            SpecError::InvalidConvectiveSpeed { u_conv } => {
                write!(f, "convective outflow u_conv = {u_conv} must lie in (0, 1]")
            }
            SpecError::NonZeroZForce2D { fz } => {
                write!(f, "2D lattice requires force[2] == 0 (got {fz})")
            }
            SpecError::OpenFaceAxisTooShort { face, extent } => write!(
                f,
                "open face {face} needs its own axis to span >= 3 cells (got {extent})"
            ),
        }
    }
}

impl std::error::Error for SpecError {}

/// A non-finite state detected by the run-time watchdog
/// ([`Solver::run_guarded`] and the GPU/MPI counterparts, A-9).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Diverged {
    /// Completed steps when the non-finite state was detected. The divergence
    /// itself occurred at most `check_every` steps earlier.
    pub step: u64,
}

impl std::fmt::Display for Diverged {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "simulation diverged: non-finite mass detected at step {}",
            self.step
        )
    }
}

impl std::error::Error for Diverged {}

impl<T: Real> GlobalSpec<T> {
    /// Validate the scenario before a solver is built (A-4). `d` is the
    /// lattice dimension (`L::D`); `solid` is the compact global solid mask
    /// (`z*(nx*ny) + y*nx + x`, empty = no solids) used to decide whether a
    /// non-open, non-periodic face is fully walled.
    ///
    /// Checks (all previously silent on the V2-native path): ν finite & > 0;
    /// TRT magic finite & > 0; body force finite and (2D) `force[2] == 0`;
    /// every active axis ≥ 3 cells; no axis both periodic and open; open faces
    /// confined to one axis; every non-periodic face open **or** a full solid
    /// rim; open-face velocity ≤ MAX_SPEED (NaN-safe); outlet ρ > 0; convective
    /// u_conv ∈ (0, 1]; each open face's own axis ≥ 3 cells.
    pub fn validate(&self, d: usize, solid: &[bool]) -> Result<(), SpecError> {
        // Viscosity (finite & positive: ν = 0 ⇒ omega_m = 0, E3).
        if !self.nu.is_finite() {
            return Err(SpecError::NonFiniteParameter { what: "nu" });
        }
        if !(self.nu > 0.0) {
            return Err(SpecError::NonPositiveViscosity { nu: self.nu });
        }
        // TRT magic (finite & positive).
        if let CollisionKind::Trt { magic } = self.collision {
            if !magic.is_finite() || !(magic > 0.0) {
                return Err(SpecError::InvalidMagic { magic });
            }
        }
        // Body force finiteness, plus the 2D z-component rule.
        for (a, comp) in self.force.iter().enumerate() {
            let v = comp.as_f64();
            if !v.is_finite() {
                return Err(SpecError::NonFiniteParameter {
                    what: match a {
                        0 => "force[0]",
                        1 => "force[1]",
                        _ => "force[2]",
                    },
                });
            }
            if a == 2 && d < 3 && v != 0.0 {
                return Err(SpecError::NonZeroZForce2D { fz: v });
            }
        }
        // Minimum extents on the active axes.
        for a in 0..d {
            if self.dims[a] < 3 {
                return Err(SpecError::DomainTooSmall { dims: self.dims });
            }
        }
        // Per-axis: periodic × open exclusivity, and gather open axes.
        let mut open_axes = 0usize;
        for a in 0..d {
            let (neg, pos) = (Face::ALL[2 * a], Face::ALL[2 * a + 1]);
            let axis_open = self.faces[neg.index()].is_open() || self.faces[pos.index()].is_open();
            if self.periodic[a] && axis_open {
                return Err(SpecError::PeriodicOpenConflict { axis: a });
            }
            if axis_open {
                open_axes += 1;
            }
        }
        if open_axes > 1 {
            return Err(SpecError::OpenFacesOnMultipleAxes);
        }
        // Per-face checks: coverage, BC parameter ranges, open-axis extent.
        for face in Face::ALL {
            let a = face.axis();
            if a >= d {
                continue;
            }
            let bc = &self.faces[face.index()];
            if self.periodic[a] {
                // Periodic axis: this face wraps, no coverage or BC needed
                // (periodic × open already rejected above).
                continue;
            }
            if bc.is_open() {
                // Open face: its own axis must span >= 3 cells (reads one cell
                // inward), and its BC parameters must be in range.
                if self.dims[a] < 3 {
                    return Err(SpecError::OpenFaceAxisTooShort {
                        face: face.index(),
                        extent: self.dims[a],
                    });
                }
                match bc {
                    FaceBC::Velocity { u } => {
                        let mut sq = 0.0f64;
                        for c in u.iter() {
                            let v = c.as_f64();
                            if !v.is_finite() {
                                return Err(SpecError::NonFiniteParameter {
                                    what: "inlet velocity",
                                });
                            }
                            sq += v * v;
                        }
                        let s = sq.sqrt();
                        if !(s <= MAX_SPEED) {
                            return Err(SpecError::VelocityTooHigh { speed: s });
                        }
                    }
                    FaceBC::Pressure { rho } => {
                        let r = rho.as_f64();
                        if !(r > 0.0) {
                            return Err(SpecError::NonPositiveDensity { rho: r });
                        }
                    }
                    FaceBC::Convective { u_conv } => {
                        let v = u_conv.as_f64();
                        if !(v > 0.0 && v <= 1.0) {
                            return Err(SpecError::InvalidConvectiveSpeed { u_conv: v });
                        }
                    }
                    FaceBC::Outflow | FaceBC::Closed => {}
                }
            } else {
                // Closed, non-periodic face: it must be a full solid wall rim,
                // else its halo feeds stale interior values (E2).
                if !face_is_full_solid_rim(face, self.dims, solid) {
                    return Err(SpecError::UncoveredFace { face: face.index() });
                }
            }
        }
        Ok(())
    }
}

/// Whether every cell on `face`'s plane is solid (a full wall rim). An empty
/// `solid` mask means no solids, so a bare non-periodic closed face is
/// uncovered.
fn face_is_full_solid_rim(face: Face, dims: [usize; 3], solid: &[bool]) -> bool {
    if solid.is_empty() {
        return false;
    }
    let a = face.axis();
    let fixed = if face.is_neg() { 0 } else { dims[a] - 1 };
    let (t1, t2) = face.tangents();
    for c2 in 0..dims[t2] {
        for c1 in 0..dims[t1] {
            let mut pos = [0usize; 3];
            pos[a] = fixed;
            pos[t1] = c1;
            pos[t2] = c2;
            let i = (pos[2] * dims[1] + pos[1]) * dims[0] + pos[0];
            if !solid[i] {
                return false;
            }
        }
    }
    true
}

/// Which global faces are walls, and their tangential velocities.
#[derive(Clone, Copy, Debug)]
pub struct WallSpec<T: Real> {
    /// Wall flag per face (`Face::index()` order).
    pub is_wall: [bool; 6],
    /// Wall velocity per face (used when the face is a wall).
    pub u: [[T; 3]; 6],
}

impl<T: Real> Default for WallSpec<T> {
    fn default() -> Self {
        Self {
            is_wall: [false; 6],
            u: [[T::zero(); 3]; 6],
        }
    }
}

/// Realise wall-type faces as one-cell solid rims over the global grid
/// (V1 `build_rims`). Where two rims share a corner cell the faster wall's
/// velocity wins (strict `>` on squared speed), so the result does not
/// depend on application order; equal speeds keep the first-applied face in
/// V1's order — bottom, top, left, right — i.e. YNeg, YPos, XNeg, XPos,
/// then ZNeg, ZPos.
///
/// Returns compact global `(solid, wall_u)` arrays.
pub fn build_wall_rims<T: Real>(
    d: usize,
    dims: [usize; 3],
    walls: &WallSpec<T>,
) -> (Vec<bool>, Vec<[T; 3]>) {
    let n = dims[0] * dims[1] * dims[2];
    let mut solid = vec![false; n];
    let mut wall_u = vec![[T::zero(); 3]; n];
    let mut best = vec![-1.0f64; n];
    const ORDER: [Face; 6] = [
        Face::YNeg,
        Face::YPos,
        Face::XNeg,
        Face::XPos,
        Face::ZNeg,
        Face::ZPos,
    ];
    for face in ORDER {
        let a = face.axis();
        if a >= d || !walls.is_wall[face.index()] {
            continue;
        }
        let u = walls.u[face.index()];
        let mut speed = u[0].as_f64().powi(2) + u[1].as_f64().powi(2);
        if d == 3 {
            speed += u[2].as_f64().powi(2);
        }
        let fixed = if face.is_neg() { 0 } else { dims[a] - 1 };
        let (t1, t2) = match a {
            0 => (1, 2),
            1 => (0, 2),
            _ => (0, 1),
        };
        for c2 in 0..dims[t2] {
            for c1 in 0..dims[t1] {
                let mut pos = [0usize; 3];
                pos[a] = fixed;
                pos[t1] = c1;
                pos[t2] = c2;
                let i = (pos[2] * dims[1] + pos[1]) * dims[0] + pos[0];
                solid[i] = true;
                if speed > best[i] {
                    best[i] = speed;
                    wall_u[i] = u;
                }
            }
        }
    }
    (solid, wall_u)
}

/// Cartesian decomposition of the global grid into `decomp[0] × decomp[1] ×
/// decomp[2]` subdomains. Remainder cells go to the lowest-index parts.
/// Part id = `(pz * decomp[1] + py) * decomp[0] + px`.
pub fn partition(
    d: usize,
    dims: [usize; 3],
    periodic: [bool; 3],
    decomp: [usize; 3],
) -> Vec<Subdomain> {
    for a in 0..3 {
        assert!(decomp[a] >= 1, "decomp must be >= 1 per axis");
        if a >= d {
            assert_eq!(decomp[a], 1, "cannot split inactive axis {a}");
        }
        assert!(decomp[a] <= dims[a], "more parts than cells on axis {a}");
    }
    // Per-axis part extents and origins.
    let mut extents: [Vec<usize>; 3] = [vec![], vec![], vec![]];
    let mut origins: [Vec<usize>; 3] = [vec![], vec![], vec![]];
    for a in 0..3 {
        let k = decomp[a];
        let base = dims[a] / k;
        let rem = dims[a] % k;
        let mut o = 0;
        for p in 0..k {
            let e = base + usize::from(p < rem);
            if k > 1 {
                assert!(
                    e >= 2,
                    "split parts must be at least 2 cells wide on axis {a}"
                );
            }
            extents[a].push(e);
            origins[a].push(o);
            o += e;
        }
        debug_assert_eq!(o, dims[a]);
    }
    let pid = |px: usize, py: usize, pz: usize| (pz * decomp[1] + py) * decomp[0] + px;
    let mut subs = Vec::with_capacity(decomp[0] * decomp[1] * decomp[2]);
    for pz in 0..decomp[2] {
        for py in 0..decomp[1] {
            for px in 0..decomp[0] {
                let pc = [px, py, pz];
                let mut neighbors = [None; 6];
                for face in Face::ALL {
                    let a = face.axis();
                    if a >= d {
                        continue;
                    }
                    let k = decomp[a];
                    let mut nb = pc;
                    let at_edge = if face.is_neg() {
                        pc[a] == 0
                    } else {
                        pc[a] == k - 1
                    };
                    if at_edge && !periodic[a] {
                        continue;
                    }
                    nb[a] = if face.is_neg() {
                        (pc[a] + k - 1) % k
                    } else {
                        (pc[a] + 1) % k
                    };
                    neighbors[face.index()] = Some(pid(nb[0], nb[1], nb[2]));
                }
                subs.push(Subdomain {
                    global: dims,
                    origin: [origins[0][px], origins[1][py], origins[2][pz]],
                    geom: crate::fields::LocalGeom::new(
                        d,
                        [extents[0][px], extents[1][py], extents[2][pz]],
                        1,
                    ),
                    neighbors,
                });
            }
        }
    }
    subs
}

/// Time-evolution driver over a decomposed grid.
pub struct Solver<L, T, B, H>
where
    L: Lattice,
    T: Real,
    B: Backend<L, T>,
    H: HaloExchange<T>,
{
    params: StepParams<T>,
    nu: f64,
    dims: [usize; 3],
    periodic: [bool; 3],
    subs: Vec<Subdomain>,
    /// Host staging fields used for setup edits and population readback.
    host_parts: Vec<SoaFields<T>>,
    /// Backend-owned compute fields.
    parts: Vec<B::Fields>,
    backend: B,
    exchange: H,
    time: u64,
    probed_force: [T; 3],
    masks_dirty: bool,
    host_dirty: bool,
    device_ahead: bool,
    psi_planes: Vec<Vec<T>>,
    gravity: Option<[T; 3]>,
    /// Split streaming into interior + boundary-shell passes (the overlap
    /// seam for asynchronous exchanges). Off by default: the single full
    /// pass reproduces V1's probe summation order bit-for-bit.
    two_pass: bool,
    _lattice: std::marker::PhantomData<L>,
}

impl<L, T, B, H> Solver<L, T, B, H>
where
    L: Lattice,
    T: Real,
    B: Backend<L, T>,
    H: HaloExchange<T>,
{
    /// Build a solver over `decomp` subdomains. `solid` / `wall_u` are
    /// compact global arrays (empty = no solids); see [`build_wall_rims`].
    ///
    /// Mirrors V1 `from_config`: quiescent deviation state, rims applied,
    /// then one `update_moments` (so `u(t=0)` includes the half-force term).
    pub fn new(
        spec: &GlobalSpec<T>,
        solid: &[bool],
        wall_u: &[[T; 3]],
        decomp: [usize; 3],
        backend: B,
        exchange: H,
    ) -> Self {
        Self::build(spec, solid, wall_u, decomp, None, backend, exchange)
    }

    /// Build a solver that *owns exactly one part* of the `decomp`
    /// decomposition (distributed-memory configuration: one process per
    /// part). Neighbour ids in the subdomain still refer to the global part
    /// numbering — the exchange implementation defines where those parts
    /// live (for MPI, part id = rank). `LocalPeriodic` / `InProcess` cannot
    /// serve such a solver (they index neighbours into the local part list).
    ///
    /// `solid` / `wall_u` are still the *global* compact arrays; every owner
    /// slices out its own core. Cell accessors (`rho`, `u`, `set_solid`, …)
    /// address global coordinates and must only be called for cells this
    /// part owns; `gather_*` fills only the owned block (the distributed
    /// gather assembles rank blocks on the root).
    pub fn new_local_part(
        spec: &GlobalSpec<T>,
        solid: &[bool],
        wall_u: &[[T; 3]],
        decomp: [usize; 3],
        part: usize,
        backend: B,
        exchange: H,
    ) -> Self {
        Self::build(spec, solid, wall_u, decomp, Some(part), backend, exchange)
    }

    fn build(
        spec: &GlobalSpec<T>,
        solid: &[bool],
        wall_u: &[[T; 3]],
        decomp: [usize; 3],
        only: Option<usize>,
        backend: B,
        exchange: H,
    ) -> Self {
        if L::D == 2 {
            assert_eq!(spec.dims[2], 1, "2D lattice requires nz == 1");
        }
        let n = spec.dims[0] * spec.dims[1] * spec.dims[2];
        assert!(solid.is_empty() || solid.len() == n);
        assert!(wall_u.is_empty() || wall_u.len() == n);
        // A-4: validate the scenario before allocating. Higher layers
        // (lbm-scenario, the compat facade) validate explicitly and surface a
        // typed error; this call is the last-line guard that turns an invalid
        // native `GlobalSpec` (uncovered face, ν = 0, periodic × open, …) into
        // a clear panic instead of silent non-physical output.
        spec.validate(L::D, solid)
            .unwrap_or_else(|e| panic!("invalid GlobalSpec: {e}"));
        let (omega_p, omega_m) = spec.collision.omegas(spec.nu);
        let params = StepParams {
            omega_p,
            omega_m,
            force: spec.force,
            faces: spec.faces,
        };
        let mut subs = partition(L::D, spec.dims, spec.periodic, decomp);
        if let Some(part) = only {
            assert!(part < subs.len(), "part {part} out of range for {decomp:?}");
            // A single-part owner keeps *global* neighbour ids in its
            // subdomain, so only a Remote exchange (MPI) can resolve them. A
            // Local exchange (LocalPeriodic/InProcess) would read a global id
            // as a local `parts` index — a silent self-wrap into part 0 when
            // the id is 0, or an out-of-bounds panic otherwise (A-5).
            assert_eq!(
                H::SCOPE,
                ExchangeScope::Remote,
                "new_local_part (single-part ownership of a {decomp:?} decomposition) requires a \
                 Remote halo exchange (e.g. MpiExchange); LocalPeriodic/InProcess resolve \
                 neighbour ids as local part indices and would silently wrap or panic"
            );
            subs = vec![subs[part].clone()];
        }
        let mut host_parts: Vec<SoaFields<T>> =
            subs.iter().map(|s| SoaFields::new(L::Q, s.geom)).collect();
        // Distribute the global masks into the parts' padded cores.
        for (sub, fields) in subs.iter().zip(host_parts.iter_mut()) {
            let g = sub.geom;
            for z in 0..g.core[2] {
                for y in 0..g.core[1] {
                    for x in 0..g.core[0] {
                        let gi = ((sub.origin[2] + z) * spec.dims[1] + (sub.origin[1] + y))
                            * spec.dims[0]
                            + (sub.origin[0] + x);
                        let pi = g.pidx(x, y, z);
                        if !solid.is_empty() {
                            fields.solid[pi] = solid[gi];
                        }
                        if !wall_u.is_empty() {
                            fields.wall_u[pi] = wall_u[gi];
                        }
                    }
                }
            }
        }
        let parts = subs.iter().map(|s| backend.alloc(s)).collect();
        let mut solver = Self {
            params,
            nu: spec.nu,
            dims: spec.dims,
            periodic: spec.periodic,
            subs,
            host_parts,
            parts,
            backend,
            exchange,
            time: 0,
            probed_force: [T::zero(); 3],
            masks_dirty: true,
            host_dirty: true,
            device_ahead: false,
            psi_planes: Vec::new(),
            gravity: None,
            two_pass: false,
            _lattice: std::marker::PhantomData,
        };
        solver.psi_planes = solver
            .subs
            .iter()
            .map(|sub| vec![T::zero(); sub.geom.n_padded()])
            .collect();
        solver.sync_masks();
        solver.stage_in_if_dirty();
        // V1 from_config ends with update_moments (u(t=0) = force/2 on fluid).
        for i in 0..solver.parts.len() {
            solver
                .backend
                .update_moments(&solver.subs[i], &mut solver.parts[i], &solver.params);
        }
        solver.device_ahead = true;
        solver
    }

    fn sync_masks(&mut self) {
        self.exchange
            .exchange_masks(&self.subs, &mut self.host_parts);
        self.host_dirty = true;
        self.masks_dirty = false;
    }

    fn sync_masks_if_dirty(&mut self) {
        if self.masks_dirty {
            self.sync_masks();
        }
    }

    fn stage_in_if_dirty(&mut self) {
        if !self.host_dirty {
            return;
        }
        for i in 0..self.parts.len() {
            self.backend
                .stage_in(&self.subs[i], &mut self.parts[i], &self.host_parts[i]);
        }
        self.host_dirty = false;
        self.device_ahead = false;
    }

    fn stage_out_all(&mut self) {
        if !self.device_ahead {
            return;
        }
        for i in 0..self.parts.len() {
            self.backend
                .stage_out(&self.subs[i], &self.parts[i], &mut self.host_parts[i]);
        }
        self.device_ahead = false;
    }

    fn stage_gravity(&mut self) -> Option<Vec<(bool, Vec<[T; 3]>)>> {
        let gvec = self.gravity?;
        let mut staged = Vec::with_capacity(self.host_parts.len());
        // Gravity is a transient host-staged overlay: stage_out_all() first
        // makes rho current, stage_in_if_dirty() uploads rho*g for the
        // backend step, and unstage_gravity() removes it from host storage.
        for (sub, fields) in self.subs.iter().zip(self.host_parts.iter_mut()) {
            let geo = sub.geom;
            let n_core = geo.n_core();
            let was_none = fields.force_field.is_none();
            let ff = fields
                .force_field
                .get_or_insert_with(|| vec![[T::zero(); 3]; n_core]);
            if ff.len() != n_core {
                ff.clear();
                ff.resize(n_core, [T::zero(); 3]);
            }
            let mut added = vec![[T::zero(); 3]; n_core];
            for z in 0..geo.core[2] {
                for y in 0..geo.core[1] {
                    for x in 0..geo.core[0] {
                        let pi = geo.pidx(x, y, z);
                        if fields.solid[pi] {
                            continue;
                        }
                        let c = geo.cidx(x, y, z);
                        let rho = fields.rho[c];
                        for a in 0..3 {
                            added[c][a] = rho * gvec[a];
                            ff[c][a] = ff[c][a] + added[c][a];
                        }
                    }
                }
            }
            staged.push((was_none, added));
        }
        self.host_dirty = true;
        Some(staged)
    }

    fn unstage_gravity(&mut self, staged: Option<Vec<(bool, Vec<[T; 3]>)>>) {
        let Some(staged) = staged else {
            return;
        };
        for ((was_none, added), fields) in staged.into_iter().zip(self.host_parts.iter_mut()) {
            let Some(ff) = fields.force_field.as_mut() else {
                continue;
            };
            for (dst, add) in ff.iter_mut().zip(added.iter()) {
                for a in 0..3 {
                    dst[a] = dst[a] - add[a];
                }
            }
            if was_none {
                fields.force_field = None;
            }
        }
        self.host_dirty = true;
    }

    fn run_staged_step(&mut self) {
        self.stage_out_all();
        let gravity_stage = self.stage_gravity();
        self.stage_in_if_dirty();
        self.backend.run_span(
            &self.exchange,
            &self.subs,
            &mut self.parts,
            &self.params,
            self.two_pass,
            &mut self.probed_force,
            1,
        );
        self.time += 1;
        self.device_ahead = true;
        self.backend.finish_run_chunk(&self.parts, 1);
        self.stage_out_all();
        self.unstage_gravity(gravity_stage);
        self.stage_in_if_dirty();
    }

    /// Advance one time step (V1 order: collide → stream → Bouzidi → swap →
    /// open faces → moments).
    pub fn step(&mut self) {
        self.sync_masks_if_dirty();
        if self.gravity.is_some() {
            self.run_staged_step();
            return;
        }
        self.stage_in_if_dirty();
        self.backend.run_span(
            &self.exchange,
            &self.subs,
            &mut self.parts,
            &self.params,
            self.two_pass,
            &mut self.probed_force,
            1,
        );
        self.time += 1;
        self.device_ahead = true;
        self.backend.finish_run_chunk(&self.parts, 1);
        if !self.backend.handles_single_part_periodic_halo() {
            self.stage_out_all();
        }
    }

    /// Advance `steps` time steps.
    pub fn run(&mut self, steps: usize) {
        if self.gravity.is_some() {
            for _ in 0..steps {
                self.step();
            }
            return;
        }
        self.sync_masks_if_dirty();
        self.stage_in_if_dirty();
        let mut remaining = steps;
        while remaining > 0 {
            let chunk = self
                .backend
                .run_chunk_size(&self.parts)
                .max(1)
                .min(remaining);
            self.backend.run_span(
                &self.exchange,
                &self.subs,
                &mut self.parts,
                &self.params,
                self.two_pass,
                &mut self.probed_force,
                chunk,
            );
            self.time += chunk as u64;
            self.device_ahead = true;
            self.backend.finish_run_chunk(&self.parts, chunk);
            if !self.backend.handles_single_part_periodic_halo() {
                self.stage_out_all();
            }
            remaining -= chunk;
        }
    }

    /// Advance `steps` steps with a periodic non-finite watchdog (A-9).
    ///
    /// Every `check_every` steps — and once more after the final step when
    /// `steps` is not a multiple — the f64 mass aggregation behind
    /// [`Solver::total_mass`] is inspected. A NaN or ±Inf anywhere in the
    /// fluid populations propagates into that sum, so a non-finite total
    /// detects the divergence **without touching the physics kernels** (they
    /// stay guard-free and V1-equivalent); the produced trajectory is
    /// bit-identical to [`Solver::run`]. `check_every == 0` is treated as 1.
    ///
    /// Cost: one extra O(N·Q) f64 reduction per check — measured < 1% of
    /// step cost at 512² with `check_every = 100` (`tests/run_guarded.rs`).
    ///
    /// On detection, returns the completed step count; the divergence
    /// occurred at most `check_every` steps earlier.
    pub fn run_guarded(&mut self, steps: usize, check_every: usize) -> Result<(), Diverged> {
        let check_every = check_every.max(1);
        let mut remaining = steps;
        while remaining > 0 {
            let chunk = remaining.min(check_every);
            self.run(chunk);
            remaining -= chunk;
            self.check_mass_finite()?;
        }
        Ok(())
    }

    fn check_mass_finite(&self) -> Result<(), Diverged> {
        let (fluid, m) = self.local_mass_partials();
        if (fluid + m).is_finite() {
            Ok(())
        } else {
            Err(Diverged { step: self.time })
        }
    }

    /// Toggle the interior/boundary two-pass streaming split.
    pub fn set_two_pass(&mut self, on: bool) {
        self.two_pass = on;
    }

    // ------------------------------------------------------------------
    // Setup (host-side staging)
    // ------------------------------------------------------------------

    /// Initialise every cell from `(rho, u) = init(x, y, z)` (global
    /// coordinates), second-order consistent: `f = feq + f_neq` with the
    /// Chapman–Enskog non-equilibrium part from central velocity
    /// differences (V1 `init_with`).
    ///
    /// V1 samples the *stored* pass-1 fields for the differences; since the
    /// stored values are exactly `init(...)`'s outputs, this implementation
    /// re-evaluates `init` at neighbour coordinates instead — bit-identical
    /// values, and it works across subdomain boundaries without a moment
    /// halo. Solid neighbours (looked up in the exchanged halo masks) fall
    /// back one-sided exactly like V1.
    pub fn init_with(&mut self, init: impl Fn(usize, usize, usize) -> (T, [T; 3])) {
        self.stage_out_all();
        self.sync_masks_if_dirty();
        let kp = KParams::new::<L>(&self.params);
        let tau = T::r(3.0 * self.nu + 0.5);
        let three = T::r(3.0);
        let half = T::r(0.5);
        let dims = self.dims;
        let periodic = self.periodic;
        for (sub, fields) in self.subs.iter().zip(self.host_parts.iter_mut()) {
            let g = sub.geom;
            let np = g.n_padded();
            // Pass 1: store the macroscopic fields (all core cells).
            for z in 0..g.core[2] {
                for y in 0..g.core[1] {
                    for x in 0..g.core[0] {
                        let (r, u) = init(sub.origin[0] + x, sub.origin[1] + y, sub.origin[2] + z);
                        let c = g.cidx(x, y, z);
                        fields.rho[c] = r;
                        fields.ux[c] = u[0];
                        fields.uy[c] = u[1];
                        fields.uz[c] = u[2];
                    }
                }
            }
            // Pass 2: f = feq + f_neq(grad u).
            for z in 0..g.core[2] {
                for y in 0..g.core[1] {
                    for x in 0..g.core[0] {
                        let c = g.cidx(x, y, z);
                        let pi = g.pidx(x, y, z);
                        let feq = equilibrium::<L, T>(
                            &kp,
                            fields.rho[c],
                            [fields.ux[c], fields.uy[c], fields.uz[c]],
                        );
                        for q in 0..L::Q {
                            fields.f[q * np + pi] = feq[q];
                        }
                        if fields.solid[pi] {
                            continue;
                        }
                        // Central differences with graceful fallback to
                        // one-sided when the neighbour is missing (wall rim /
                        // non-periodic domain edge) — V1 `sample`/`diff`.
                        let sample = |da: [isize; 3]| -> Option<[T; 3]> {
                            let mut gpos = [0isize; 3];
                            for a in 0..3 {
                                gpos[a] = [x, y, z][a] as isize + sub.origin[a] as isize + da[a];
                                if gpos[a] < 0 || gpos[a] >= dims[a] as isize {
                                    if a < L::D && periodic[a] {
                                        gpos[a] = (gpos[a] + dims[a] as isize) % dims[a] as isize;
                                    } else {
                                        return None;
                                    }
                                }
                            }
                            // Solid lookup via the local halo (exchanged).
                            let lp = g.pidx_i(
                                x as isize + da[0],
                                y as isize + da[1],
                                z as isize + da[2],
                            );
                            if fields.solid[lp] {
                                return None;
                            }
                            let (_, u) = init(gpos[0] as usize, gpos[1] as usize, gpos[2] as usize);
                            Some(u)
                        };
                        let own = [fields.ux[c], fields.uy[c], fields.uz[c]];
                        let diff = |plus: Option<[T; 3]>, minus: Option<[T; 3]>, b: usize| -> T {
                            match (plus, minus) {
                                (Some(pv), Some(mv)) => (pv[b] - mv[b]) * half,
                                (Some(pv), None) => pv[b] - own[b],
                                (None, Some(mv)) => own[b] - mv[b],
                                (None, None) => T::zero(),
                            }
                        };
                        // grad[a][b] = d u_b / d x_a.
                        let mut grad = [[T::zero(); 3]; 3];
                        for a in 0..L::D {
                            let mut dp = [0isize; 3];
                            dp[a] = 1;
                            let mut dm = [0isize; 3];
                            dm[a] = -1;
                            let (pv, mv) = (sample(dp), sample(dm));
                            for b in 0..L::D {
                                grad[a][b] = diff(pv, mv, b);
                            }
                        }
                        let mut div = grad[0][0];
                        for a in 1..L::D {
                            div = div + grad[a][a];
                        }
                        for q in 0..L::Q {
                            // ccgu = sum_ab c_a c_b (grad symmetrised),
                            // accumulated in V1's (0,0), (0,1), (1,1) order.
                            let cq = kp.cr[q];
                            let mut ccgu = cq[0] * cq[0] * grad[0][0];
                            for a in 0..L::D {
                                for b in a..L::D {
                                    if a == 0 && b == 0 {
                                        continue;
                                    }
                                    if a == b {
                                        ccgu = ccgu + cq[a] * cq[a] * grad[a][a];
                                    } else {
                                        ccgu = ccgu + cq[a] * cq[b] * (grad[a][b] + grad[b][a]);
                                    }
                                }
                            }
                            let fneq = -kp.wr[q] * fields.rho[c] * tau * (three * ccgu - div);
                            fields.f[q * np + pi] = fields.f[q * np + pi] + fneq;
                        }
                    }
                }
            }
        }
        self.host_dirty = true;
        self.stage_in_if_dirty();
        for i in 0..self.parts.len() {
            self.backend
                .update_moments(&self.subs[i], &mut self.parts[i], &self.params);
        }
        self.device_ahead = true;
    }

    /// Mark a global cell solid (half-way bounce-back obstacle). Open-face
    /// checks are the caller's (facade's) responsibility.
    pub fn set_solid(&mut self, x: usize, y: usize, z: usize) {
        let (i, lx, ly, lz) = self.locate(x, y, z);
        let pi = self.subs[i].geom.pidx(lx, ly, lz);
        self.stage_out_all();
        self.host_parts[i].solid[pi] = true;
        self.masks_dirty = true;
        self.host_dirty = true;
        self.stage_in_if_dirty();
    }

    /// Build analytic Bouzidi records for a circle. Solid cells must already
    /// be marked with the same geometry.
    pub fn set_bouzidi_circle(&mut self, cx: f64, cy: f64, r: f64) {
        self.stage_out_all();
        self.sync_masks_if_dirty();
        for (sub, fields) in self.subs.iter().zip(self.host_parts.iter_mut()) {
            let links =
                crate::bouzidi::circle_links(&fields.geom, sub.origin, &fields.solid, cx, cy, r);
            fields.bouzidi = (!links.is_empty()).then_some(links);
        }
        self.host_dirty = true;
    }

    /// Build analytic Bouzidi records for a sphere. Solid cells must already
    /// be marked with the same geometry.
    pub fn set_bouzidi_sphere(&mut self, cx: f64, cy: f64, cz: f64, r: f64) {
        self.stage_out_all();
        self.sync_masks_if_dirty();
        for (sub, fields) in self.subs.iter().zip(self.host_parts.iter_mut()) {
            let links = crate::bouzidi::sphere_links::<T, L>(
                &fields.geom,
                sub.origin,
                &fields.solid,
                cx,
                cy,
                cz,
                r,
            );
            fields.bouzidi = (!links.is_empty()).then_some(links);
        }
        self.host_dirty = true;
    }

    /// Install qd=1/2 records for every fluid-solid link. This is intended as
    /// a degeneracy regression for bitwise equivalence to half-way BB.
    pub fn set_bouzidi_half_way_links(&mut self) {
        self.stage_out_all();
        self.sync_masks_if_dirty();
        for fields in self.host_parts.iter_mut() {
            let links = crate::bouzidi::half_way_links::<T, L>(&fields.geom, &fields.solid);
            fields.bouzidi = (!links.is_empty()).then_some(links);
        }
        self.host_dirty = true;
    }

    /// Remove all Bouzidi records; subsequent steps use pure half-way BB.
    pub fn clear_bouzidi(&mut self) {
        self.stage_out_all();
        for fields in self.host_parts.iter_mut() {
            fields.bouzidi = None;
        }
        self.host_dirty = true;
    }

    /// Select the solid cells whose momentum-exchange force is accumulated
    /// each step (V1 `set_force_probe`).
    pub fn set_force_probe(&mut self, pred: impl Fn(usize, usize, usize) -> bool) {
        self.stage_out_all();
        for (sub, fields) in self.subs.iter().zip(self.host_parts.iter_mut()) {
            let g = sub.geom;
            let mut mask = vec![false; g.n_padded()];
            for z in 0..g.core[2] {
                for y in 0..g.core[1] {
                    for x in 0..g.core[0] {
                        mask[g.pidx(x, y, z)] =
                            pred(sub.origin[0] + x, sub.origin[1] + y, sub.origin[2] + z);
                    }
                }
            }
            fields.probe = Some(mask);
        }
        self.masks_dirty = true;
        self.host_dirty = true;
    }

    /// Prescribe the per-cell body force (Guo forcing) from a closure over
    /// global cell coordinates `(x, y, z)`. The closure is evaluated once per
    /// owned core cell and stored in the part's compact layout, so the result
    /// is decomposition-invariant (identical global field for any `decomp`).
    /// Existing allocations are reused; call it before [`Solver::step`] each
    /// time the field changes (e.g. a time-dependent force). The force enters
    /// collision with the usual Guo half-force velocity correction, so
    /// `u(x)` accessors keep returning the physical velocity.
    ///
    /// Unlike [`Solver::update_shan_chen_force`] this stencil is purely local
    /// (no neighbour reads, no halo exchange): it is the general hook for
    /// spatially/temporally varying forcing — uniform or linear forcing,
    /// sponge/absorbing layers, and volume-penalization (Brinkman) regions
    /// that relax the local velocity toward a prescribed target.
    pub fn set_body_force_field(&mut self, f: impl Fn(usize, usize, usize) -> [T; 3]) {
        self.stage_out_all();
        for (sub, fields) in self.subs.iter().zip(self.host_parts.iter_mut()) {
            let g = sub.geom;
            let n_core = g.n_core();
            let buf = fields
                .force_field
                .get_or_insert_with(|| vec![[T::zero(); 3]; n_core]);
            if buf.len() != n_core {
                buf.clear();
                buf.resize(n_core, [T::zero(); 3]);
            }
            for z in 0..g.core[2] {
                for y in 0..g.core[1] {
                    for x in 0..g.core[0] {
                        buf[g.cidx(x, y, z)] =
                            f(sub.origin[0] + x, sub.origin[1] + y, sub.origin[2] + z);
                    }
                }
            }
        }
        self.host_dirty = true;
    }

    /// Drop the per-cell body force field on every owned part (subsequent
    /// steps run force-free unless [`GlobalSpec::force`] is nonzero).
    pub fn clear_body_force_field(&mut self) {
        self.stage_out_all();
        for fields in self.host_parts.iter_mut() {
            fields.force_field = None;
        }
        self.host_dirty = true;
    }

    /// Set per-mass gravity `g`; at the start of each step, `rho(x) * g` is
    /// added to the per-cell force on fluid cells only.
    pub fn set_gravity(&mut self, g: [T; 3]) {
        self.stage_out_all();
        self.gravity = Some(g);
    }

    /// Prescribe a per-node inlet profile on a `Velocity` face, `values`
    /// indexed by the global along-face coordinate in canonical face order:
    /// with tangent axes `(t1, t2) = face.tangents()`, the index is
    /// `c2 * dims[t1] + c1` (`t1` fastest). For 2D lattices `dims[t2] == 1`,
    /// so this degenerates to the single tangent coordinate (V1 convention).
    pub fn set_inlet_profile(&mut self, face: Face, values: &[[T; 3]]) {
        assert!(
            matches!(self.params.faces[face.index()], FaceBC::Velocity { .. }),
            "set_inlet_profile: {face:?} is not a Velocity face"
        );
        let (t1, t2) = face.tangents();
        assert_eq!(
            values.len(),
            self.dims[t1] * self.dims[t2],
            "profile must cover the whole global face"
        );
        self.stage_out_all();
        for (sub, fields) in self.subs.iter().zip(self.host_parts.iter_mut()) {
            if !sub.touches_global_face(face) {
                fields.inlet_profiles[face.index()] = None;
                continue;
            }
            let (o1, o2) = (sub.origin[t1], sub.origin[t2]);
            let (e1, e2) = (sub.geom.core[t1], sub.geom.core[t2]);
            let mut local = Vec::with_capacity(e1 * e2);
            for c2 in 0..e2 {
                for c1 in 0..e1 {
                    local.push(values[(o2 + c2) * self.dims[t1] + (o1 + c1)]);
                }
            }
            fields.inlet_profiles[face.index()] = Some(local);
        }
        self.host_dirty = true;
    }

    /// Closure form of [`Solver::set_inlet_profile`]: `profile(c1, c2)` is
    /// evaluated at the global tangent coordinates of every face node
    /// (`(t1, t2) = face.tangents()`; 2D faces always pass `c2 = 0`).
    /// The natural way to build e.g. a rectangular-duct profile
    /// `u(y, z) = umax f(y) g(z)` on an X face.
    pub fn set_inlet_profile_with(&mut self, face: Face, profile: impl Fn(usize, usize) -> [T; 3]) {
        let (t1, t2) = face.tangents();
        let mut values = Vec::with_capacity(self.dims[t1] * self.dims[t2]);
        for c2 in 0..self.dims[t2] {
            for c1 in 0..self.dims[t1] {
                values.push(profile(c1, c2));
            }
        }
        self.set_inlet_profile(face, &values);
    }

    /// Single-component Shan–Chen cohesion: recompute the per-cell force
    /// field from the current density via the pseudopotential `psi`,
    /// exchanging one padded ψ plane per part (`HaloExchange::exchange_scalar`)
    /// so the force stencil sees remote neighbours — the decomposition-aware
    /// counterpart of `compat::ShanChen::update_force` (neutral walls: solid
    /// and out-of-domain neighbours contribute nothing).
    ///
    /// `F(x) = -G ψ(x) Σ_q w_q ψ(x + c_q) c_q`, accumulated in ascending-`q`
    /// order (V1 convention). Call before each [`Solver::step`]; collective
    /// over all owners of the decomposition.
    pub fn update_shan_chen_force(&mut self, g: T, psi: impl Fn(T) -> T) {
        self.update_shan_chen_force_with_walls(g, T::zero(), T::zero(), psi);
    }

    /// Wall-adhesion variant of [`Solver::update_shan_chen_force`] — the
    /// native port of V1 `ShanChen::with_wall` (`g_wall`) and
    /// `ShanChen::with_wall_rho` (virtual wall density; pass the
    /// pre-evaluated `psi_wall = ψ(wall_rho)`, or zero to disable):
    ///
    /// `F(x) = -ψ(x) [ G ( Σ_{q:fluid} w_q ψ(x+c_q) c_q + Σ_{q:solid} w_q ψ_wall c_q )
    ///                 + G_w Σ_{q:solid} w_q c_q ]`
    ///
    /// Solid neighbours feed the cohesion sum with `psi_wall` (contact-angle
    /// control) plus the legacy `g_wall` adhesion term; out-of-domain
    /// neighbours on non-periodic global edges contribute nothing to either
    /// sum (zero-gradient approximation). Operand order is V1-identical, so
    /// a monolithic run reproduces `compat::ShanChen::update_force`
    /// bit-exactly; halo solids are covered by the mask exchange.
    pub fn update_shan_chen_force_with_walls(
        &mut self,
        g: T,
        g_wall: T,
        psi_wall: T,
        psi: impl Fn(T) -> T,
    ) {
        self.sync_masks_if_dirty();
        // ψ planes, padded: core = ψ(rho) (0 on solids), halo = 0 until the
        // exchange fills it (stays 0 outside non-periodic global edges,
        // matching V1's "out-of-domain contributes nothing").
        for ((sub, fields), plane) in self
            .subs
            .iter()
            .zip(self.host_parts.iter())
            .zip(self.psi_planes.iter_mut())
        {
            let geo = sub.geom;
            plane.fill(T::zero());
            for z in 0..geo.core[2] {
                for y in 0..geo.core[1] {
                    for x in 0..geo.core[0] {
                        let pi = geo.pidx(x, y, z);
                        if !fields.solid[pi] {
                            plane[pi] = psi(fields.rho[geo.cidx(x, y, z)]);
                        }
                    }
                }
            }
        }
        if self.psi_planes.len() == 1 {
            let mut plane = self.psi_planes[0].as_mut_slice();
            self.exchange
                .exchange_scalar(&self.subs, std::slice::from_mut(&mut plane));
        } else {
            let mut refs: Vec<&mut [T]> = self
                .psi_planes
                .iter_mut()
                .map(|p| p.as_mut_slice())
                .collect();
            self.exchange.exchange_scalar(&self.subs, &mut refs);
        }
        self.stage_out_all();
        // Neutral walls keep the exact historical expression (no adhesion
        // term appended), so pre-walls callers stay bit-identical.
        let wet = g_wall != T::zero() || psi_wall != T::zero();
        for (i, (sub, fields)) in self.subs.iter().zip(self.host_parts.iter_mut()).enumerate() {
            let geo = sub.geom;
            let plane = &self.psi_planes[i];
            let ff = fields
                .force_field
                .get_or_insert_with(|| vec![[T::zero(); 3]; geo.n_core()]);
            for z in 0..geo.core[2] {
                for y in 0..geo.core[1] {
                    for x in 0..geo.core[0] {
                        let c = geo.cidx(x, y, z);
                        let psi_i = plane[geo.pidx(x, y, z)];
                        if fields.solid[geo.pidx(x, y, z)] || psi_i == T::zero() {
                            ff[c] = [T::zero(); 3];
                            continue;
                        }
                        let mut s = [T::zero(); 3];
                        let mut adh = [T::zero(); 3];
                        for q in 1..L::Q {
                            let cq = L::C[q];
                            let pi = geo.pidx_i(
                                x as isize + cq[0] as isize,
                                y as isize + cq[1] as isize,
                                z as isize + cq[2] as isize,
                            );
                            let w = T::r(L::W[q]);
                            if wet && fields.solid[pi] {
                                // V1: the virtual wall density feeds the
                                // cohesion sum; g_wall adds the legacy
                                // adhesion term on top. (Halo solids are
                                // synced; non-periodic out-of-domain halos
                                // are never solid, hence contribute nothing.)
                                for a in 0..L::D {
                                    s[a] = s[a] + w * psi_wall * T::r(cq[a] as f64);
                                    adh[a] = adh[a] + w * T::r(cq[a] as f64);
                                }
                            } else {
                                let pj = plane[pi];
                                for a in 0..L::D {
                                    s[a] = s[a] + w * pj * T::r(cq[a] as f64);
                                }
                            }
                        }
                        for a in 0..3 {
                            ff[c][a] = if wet {
                                -(psi_i * (g * s[a] + g_wall * adh[a]))
                            } else {
                                -(psi_i * (g * s[a]))
                            };
                        }
                    }
                }
            }
        }
        self.host_dirty = true;
    }

    // ------------------------------------------------------------------
    // Accessors / diagnostics
    // ------------------------------------------------------------------

    fn locate(&self, x: usize, y: usize, z: usize) -> (usize, usize, usize, usize) {
        debug_assert!(x < self.dims[0] && y < self.dims[1] && z < self.dims[2]);
        // Fast path only for a truly monolithic part (a single *local* part
        // of a wider decomposition has a non-trivial origin).
        if self.subs.len() == 1 && self.subs[0].geom.core == self.dims {
            return (0, x, y, z);
        }
        for (i, s) in self.subs.iter().enumerate() {
            let inside = (0..3).all(|a| {
                let c = [x, y, z][a];
                c >= s.origin[a] && c < s.origin[a] + s.geom.core[a]
            });
            if inside {
                return (i, x - s.origin[0], y - s.origin[1], z - s.origin[2]);
            }
        }
        unreachable!("cell ({x},{y},{z}) not covered by any subdomain")
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
    /// Global grid extents.
    pub fn dims(&self) -> [usize; 3] {
        self.dims
    }
    /// Number of subdomains.
    pub fn part_count(&self) -> usize {
        self.parts.len()
    }
    /// Backend reference (used by backend-specific compatibility shims).
    pub fn backend(&self) -> &B {
        &self.backend
    }
    /// Mutable backend reference (used by backend-specific compatibility shims).
    pub fn backend_mut(&mut self) -> &mut B {
        &mut self.backend
    }
    /// Backend-owned fields of part `i`.
    pub fn backend_fields(&self, i: usize) -> &B::Fields {
        &self.parts[i]
    }
    /// Subdomain descriptor `i`.
    pub fn sub(&self, i: usize) -> &Subdomain {
        &self.subs[i]
    }
    /// Fields of part `i` (host staging; padded mask edits must go through
    /// `set_solid` / `set_force_probe` so halos stay in sync).
    pub fn fields(&self, i: usize) -> &SoaFields<T> {
        &self.host_parts[i]
    }
    /// Mutable fields of part `i` (see [`Solver::fields`] caveat).
    pub fn fields_mut(&mut self, i: usize) -> &mut SoaFields<T> {
        self.stage_out_all();
        self.host_dirty = true;
        &mut self.host_parts[i]
    }

    /// Synchronize backend-owned populations and moments into host staging.
    /// Device backends use this only at explicit read/edit boundaries.
    pub fn sync_host(&mut self) {
        self.stage_out_all();
    }

    /// Set or clear the per-cell symmetric relaxation-rate field
    /// (`omega_plus = 1/tau`) in global compact order.
    ///
    /// The field is compact and solver-level by design: collision kernels only
    /// replace the local `omega_plus` fetch when this field is present. A
    /// `None` field uses the original uniform-rate path.
    pub fn set_omega_field(&mut self, omega: Option<&[T]>) {
        let n = self.dims[0] * self.dims[1] * self.dims[2];
        if let Some(values) = omega {
            assert_eq!(values.len(), n, "omega field length must match cell count");
        }
        for (sub, fields) in self.subs.iter().zip(self.parts.iter_mut()) {
            let g = sub.geom;
            match omega {
                Some(values) => {
                    let local = fields
                        .omega_field
                        .get_or_insert_with(|| vec![T::zero(); g.n_core()]);
                    if local.len() != g.n_core() {
                        local.resize(g.n_core(), T::zero());
                    }
                    for z in 0..g.core[2] {
                        for y in 0..g.core[1] {
                            for x in 0..g.core[0] {
                                let gi = ((sub.origin[2] + z) * self.dims[1] + (sub.origin[1] + y))
                                    * self.dims[0]
                                    + (sub.origin[0] + x);
                                local[g.cidx(x, y, z)] = values[gi];
                            }
                        }
                    }
                }
                None => fields.omega_field = None,
            }
        }
    }

    /// Momentum-exchange force on the probed solids during the most recent
    /// step (V1 `probed_force`).
    pub fn probed_force(&self) -> [T; 3] {
        self.probed_force
    }

    /// Density at a global cell.
    pub fn rho(&self, x: usize, y: usize, z: usize) -> T {
        let (i, lx, ly, lz) = self.locate(x, y, z);
        let g = self.subs[i].geom;
        let mut hm = HostMoments::default();
        self.backend.read_moments(&self.parts[i], &mut hm);
        hm.rho[g.cidx(lx, ly, lz)]
    }
    /// Velocity at a global cell (physical, half-force corrected).
    pub fn u(&self, x: usize, y: usize, z: usize) -> [T; 3] {
        let (i, lx, ly, lz) = self.locate(x, y, z);
        let g = self.subs[i].geom;
        let c = g.cidx(lx, ly, lz);
        let mut hm = HostMoments::default();
        self.backend.read_moments(&self.parts[i], &mut hm);
        [hm.ux[c], hm.uy[c], hm.uz[c]]
    }
    /// Whether a global cell is solid.
    pub fn is_solid(&self, x: usize, y: usize, z: usize) -> bool {
        let (i, lx, ly, lz) = self.locate(x, y, z);
        self.host_parts[i].solid[self.subs[i].geom.pidx(lx, ly, lz)]
    }

    /// Total mass over fluid cells (V1 `total_mass`: physical mass =
    /// fluid-cell count + deviation sum, both accumulated in `f64`).
    pub fn total_mass(&self) -> T {
        let (fluid, m) = self.local_mass_partials();
        T::r(fluid + m)
    }

    /// Local partial sums behind [`Solver::total_mass`]: `(fluid_cells,
    /// mass_deviation)` over the parts owned by this process, in `f64`.
    /// A distributed owner sums these across ranks (order-insensitive up to
    /// f64 reassociation) before forming `fluid + m`.
    pub fn local_mass_partials(&self) -> (f64, f64) {
        let mut fluid = 0.0f64;
        let mut m = 0.0f64;
        for (sub, fields) in self.subs.iter().zip(self.parts.iter()) {
            fluid += self
                .backend
                .reduce(sub, fields, &self.params, Reduction::FluidCells);
            m += self
                .backend
                .reduce(sub, fields, &self.params, Reduction::MassDeviation);
        }
        (fluid, m)
    }

    /// Total physical momentum over fluid cells (V1 `total_momentum`).
    pub fn total_momentum(&self) -> [T; 3] {
        let p = self.local_momentum_partials();
        [T::r(p[0]), T::r(p[1]), T::r(p[2])]
    }

    /// Local partial sums behind [`Solver::total_momentum`] (see
    /// [`Solver::local_mass_partials`] for the distributed contract).
    pub fn local_momentum_partials(&self) -> [f64; 3] {
        let mut p = [0.0f64; 3];
        for (sub, fields) in self.subs.iter().zip(self.parts.iter()) {
            for (a, pa) in p.iter_mut().enumerate() {
                *pa += self
                    .backend
                    .reduce(sub, fields, &self.params, Reduction::Momentum(a));
            }
        }
        if let Some(g) = self.gravity {
            for (sub, fields) in self.subs.iter().zip(self.host_parts.iter()) {
                let geo = sub.geom;
                for z in 0..geo.core[2] {
                    for y in 0..geo.core[1] {
                        for x in 0..geo.core[0] {
                            let pi = geo.pidx(x, y, z);
                            if fields.solid[pi] {
                                continue;
                            }
                            let rho = fields.rho[geo.cidx(x, y, z)].as_f64();
                            for a in 0..3 {
                                p[a] += 0.5 * rho * g[a].as_f64();
                            }
                        }
                    }
                }
            }
        }
        p
    }

    /// Number of non-finite (NaN/Inf) values in this process's parts, over
    /// the populations and the macroscopic moments. `0` on a healthy run;
    /// a distributed owner sums the counts across ranks.
    pub fn local_nonfinite_count(&self) -> u64 {
        let mut n = 0u64;
        for fields in &self.host_parts {
            let finite = |v: &[T]| v.iter().filter(|x| !x.is_finite()).count() as u64;
            n += finite(&fields.f);
            n += finite(&fields.rho);
            n += finite(&fields.ux);
            n += finite(&fields.uy);
            n += finite(&fields.uz);
        }
        n
    }

    /// Force a mask-halo refresh before the next step. Distributed owners
    /// call this on *every* rank when any rank edits masks: the refresh is a
    /// collective exchange, so the dirty flag must agree globally.
    pub fn mark_masks_dirty(&mut self) {
        self.masks_dirty = true;
    }

    /// Number of fluid (non-solid) cells.
    pub fn fluid_cell_count(&self) -> usize {
        let mut n = 0.0;
        for (sub, fields) in self.subs.iter().zip(self.parts.iter()) {
            n += self
                .backend
                .reduce(sub, fields, &self.params, Reduction::FluidCells);
        }
        n as usize
    }

    /// Assemble a global compact array from backend-read moment planes.
    fn gather_moment(&self, get: impl Fn(&HostMoments<T>, usize) -> T) -> Vec<T> {
        let mut out = vec![T::zero(); self.dims[0] * self.dims[1] * self.dims[2]];
        for (sub, fields) in self.subs.iter().zip(self.parts.iter()) {
            let g = sub.geom;
            let mut hm = HostMoments::default();
            self.backend.read_moments(fields, &mut hm);
            for z in 0..g.core[2] {
                for y in 0..g.core[1] {
                    for x in 0..g.core[0] {
                        let gi = ((sub.origin[2] + z) * self.dims[1] + (sub.origin[1] + y))
                            * self.dims[0]
                            + (sub.origin[0] + x);
                        out[gi] = get(&hm, g.cidx(x, y, z));
                    }
                }
            }
        }
        out
    }

    /// Global density field (compact layout).
    pub fn gather_rho(&self) -> Vec<T> {
        self.gather_moment(|m, c| m.rho[c])
    }
    /// Global x-velocity field.
    pub fn gather_ux(&self) -> Vec<T> {
        self.gather_moment(|m, c| m.ux[c])
    }
    /// Global y-velocity field.
    pub fn gather_uy(&self) -> Vec<T> {
        self.gather_moment(|m, c| m.uy[c])
    }
    /// Global z-velocity field.
    pub fn gather_uz(&self) -> Vec<T> {
        self.gather_moment(|m, c| m.uz[c])
    }

    fn strain_rate_at(&self, fields: &SoaFields<T>, x: usize, y: usize, z: usize) -> [T; 6] {
        let g = fields.geom;
        let pi = g.pidx(x, y, z);
        if fields.solid[pi] {
            return [T::zero(); 6];
        }
        let c = g.cidx(x, y, z);
        let r = fields.rho[c];
        let u = [fields.ux[c], fields.uy[c], fields.uz[c]];
        let kp = KParams::new::<L>(&self.params);
        let feq = equilibrium::<L, T>(&kp, r, u);
        let np = g.n_padded();
        let mut pi_neq = [T::zero(); 6];
        for (q, cq) in L::C.iter().enumerate().take(L::Q) {
            let fneq = fields.f[q * np + pi] - feq[q];
            let cx = T::r(cq[0] as f64);
            let cy = T::r(cq[1] as f64);
            let cz = T::r(cq[2] as f64);
            pi_neq[0] = pi_neq[0] + cx * cx * fneq;
            pi_neq[1] = pi_neq[1] + cy * cy * fneq;
            pi_neq[2] = pi_neq[2] + cz * cz * fneq;
            pi_neq[3] = pi_neq[3] + cx * cy * fneq;
            pi_neq[4] = pi_neq[4] + cx * cz * fneq;
            pi_neq[5] = pi_neq[5] + cy * cz * fneq;
        }
        let force = match fields.force_field.as_ref() {
            Some(ff) => [
                self.params.force[0] + ff[c][0],
                self.params.force[1] + ff[c][1],
                self.params.force[2] + ff[c][2],
            ],
            None => self.params.force,
        };
        let half = T::r(0.5);
        // FR-STRESS-01 rev.4: Pi_force = -(dt/2)(uF + Fu), dt=1, so
        // Pi_neq_corr = Pi_neq_raw - Pi_force = Pi_neq_raw + 0.5(uF + Fu).
        pi_neq[0] = pi_neq[0] + u[0] * force[0];
        pi_neq[1] = pi_neq[1] + u[1] * force[1];
        pi_neq[2] = pi_neq[2] + u[2] * force[2];
        pi_neq[3] = pi_neq[3] + half * (u[0] * force[1] + u[1] * force[0]);
        pi_neq[4] = pi_neq[4] + half * (u[0] * force[2] + u[2] * force[0]);
        pi_neq[5] = pi_neq[5] + half * (u[1] * force[2] + u[2] * force[1]);

        let tau_eff = T::r(1.0 / self.params.omega_p);
        let scale = -(T::one() / (T::r(2.0 * L::CS2) * r * tau_eff));
        for v in &mut pi_neq {
            *v = *v * scale;
        }
        if L::D == 2 {
            pi_neq[2] = T::zero();
            pi_neq[4] = T::zero();
            pi_neq[5] = T::zero();
        }
        pi_neq
    }

    /// Global strain-rate tensor in compact cell order.
    ///
    /// Components are `[S_xx, S_yy, S_zz, S_xy, S_xz, S_yz]`. The value is
    /// evaluated from the read-only post-streaming / pre-collision
    /// populations currently stored in the solver, using the physical
    /// half-force-corrected velocity for `f_eq`. Solid cells return zeros.
    ///
    /// For TRT, the viscous stress is carried by the even/symmetric modes;
    /// therefore `tau_eff = 1 / omega_plus` (`StepParams::omega_p`). This is
    /// currently the global relaxation time, structured so a future per-cell
    /// `omega_plus` field can replace the scalar in this denominator.
    ///
    /// The Guo force correction follows FR-STRESS-01 rev.4 for this engine's
    /// deviation-form `f_eq`: `Pi_force = -0.5 * (uF + Fu)`, so the corrected
    /// non-equilibrium moment is `Pi_neq_raw + 0.5 * (uF + Fu)`.
    pub fn gather_strain_rate(&self) -> Vec<[T; 6]> {
        let mut out = vec![[T::zero(); 6]; self.dims[0] * self.dims[1] * self.dims[2]];
        for (sub, fields) in self.subs.iter().zip(self.host_parts.iter()) {
            let g = sub.geom;
            for z in 0..g.core[2] {
                for y in 0..g.core[1] {
                    for x in 0..g.core[0] {
                        let gi = ((sub.origin[2] + z) * self.dims[1] + (sub.origin[1] + y))
                            * self.dims[0]
                            + (sub.origin[0] + x);
                        out[gi] = self.strain_rate_at(fields, x, y, z);
                    }
                }
            }
        }
        out
    }

    /// Global velocity-gradient tensor in compact cell order.
    ///
    /// Each entry is `g[i][j] = du_i/dx_j`. Off-diagonal symmetric shear is
    /// taken from [`Solver::gather_strain_rate`]'s non-equilibrium stress
    /// path. Diagonal entries and the antisymmetric rotation are reconstructed
    /// from velocity differences because the native stress observable contains
    /// no vorticity, and D3Q19 moving-wall-adjacent normal stresses can carry
    /// small pure-shear artifacts.
    pub fn gather_velocity_gradient(&self) -> Vec<[[T; 3]; 3]> {
        let n = self.dims[0] * self.dims[1] * self.dims[2];
        let strain = self.gather_strain_rate();
        let ux = self.gather_ux();
        let uy = self.gather_uy();
        let uz = self.gather_uz();
        let idx =
            |x: usize, y: usize, z: usize| -> usize { (z * self.dims[1] + y) * self.dims[0] + x };
        let cell_state = |x: usize, y: usize, z: usize| -> Option<([T; 3], bool)> {
            let (i, lx, ly, lz) = self.locate(x, y, z);
            let g = self.subs[i].geom;
            let pi = g.pidx(lx, ly, lz);
            if self.parts[i].solid[pi] {
                Some((self.parts[i].wall_u[pi], true))
            } else {
                let ci = idx(x, y, z);
                Some(([ux[ci], uy[ci], uz[ci]], false))
            }
        };
        let neighbor =
            |x: usize, y: usize, z: usize, a: usize, da: isize| -> Option<([T; 3], bool)> {
                let mut p = [x as isize, y as isize, z as isize];
                p[a] += da;
                if p[a] < 0 || p[a] >= self.dims[a] as isize {
                    if a < L::D && self.periodic[a] {
                        p[a] = (p[a] + self.dims[a] as isize) % self.dims[a] as isize;
                    } else {
                        return None;
                    }
                }
                let (xx, yy, zz) = (p[0] as usize, p[1] as usize, p[2] as usize);
                cell_state(xx, yy, zz)
            };
        let mut out = vec![[[T::zero(); 3]; 3]; n];
        for z in 0..self.dims[2] {
            for y in 0..self.dims[1] {
                for x in 0..self.dims[0] {
                    let i = idx(x, y, z);
                    if self.is_solid(x, y, z) {
                        continue;
                    }
                    let s = strain[i];
                    let mut sm = [[T::zero(); 3]; 3];
                    sm[0][0] = s[0];
                    sm[1][1] = s[1];
                    sm[2][2] = s[2];
                    sm[0][1] = s[3];
                    sm[1][0] = s[3];
                    sm[0][2] = s[4];
                    sm[2][0] = s[4];
                    sm[1][2] = s[5];
                    sm[2][1] = s[5];
                    let mut fd = [[T::zero(); 3]; 3];
                    let own = [ux[i], uy[i], uz[i]];
                    for comp in 0..L::D {
                        for a in 0..L::D {
                            let plus = neighbor(x, y, z, a, 1);
                            let minus = neighbor(x, y, z, a, -1);
                            fd[comp][a] = match (plus, minus) {
                                (Some((pv, false)), Some((mv, false))) => {
                                    (pv[comp] - mv[comp]) * T::r(0.5)
                                }
                                (Some((pv, false)), Some((wv, true))) => {
                                    -T::r(4.0 / 3.0) * wv[comp]
                                        + own[comp]
                                        + T::r(1.0 / 3.0) * pv[comp]
                                }
                                (Some((wv, true)), Some((mv, false))) => {
                                    T::r(4.0 / 3.0) * wv[comp]
                                        - own[comp]
                                        - T::r(1.0 / 3.0) * mv[comp]
                                }
                                (Some((pv, false)), None) => pv[comp] - own[comp],
                                (None, Some((mv, false))) => own[comp] - mv[comp],
                                _ => T::zero(),
                            };
                        }
                    }
                    for row in 0..L::D {
                        for col in 0..L::D {
                            if row == col {
                                out[i][row][col] = fd[row][col];
                            } else {
                                let w = T::r(0.5) * (fd[row][col] - fd[col][row]);
                                out[i][row][col] = sm[row][col] + w;
                            }
                        }
                    }
                }
            }
        }
        out
    }

    /// Global shear-rate invariant `gamma_dot = sqrt(2 S:S)`.
    ///
    /// Uses [`Solver::gather_strain_rate`]'s stage, force correction and
    /// solid-cell convention.
    pub fn gather_shear_rate(&self) -> Vec<T> {
        self.gather_strain_rate()
            .into_iter()
            .map(|s| {
                let ss = s[0] * s[0]
                    + s[1] * s[1]
                    + s[2] * s[2]
                    + T::r(2.0) * (s[3] * s[3] + s[4] * s[4] + s[5] * s[5]);
                (T::r(2.0) * ss).sqrt()
            })
            .collect()
    }

    /// Global deviation-population plane `q` (compact layout).
    pub fn gather_f(&self, q: usize) -> Vec<T> {
        let mut out = vec![T::zero(); self.dims[0] * self.dims[1] * self.dims[2]];
        for (sub, fields) in self.subs.iter().zip(self.host_parts.iter()) {
            let g = sub.geom;
            let np = g.n_padded();
            for z in 0..g.core[2] {
                for y in 0..g.core[1] {
                    for x in 0..g.core[0] {
                        let gi = ((sub.origin[2] + z) * self.dims[1] + (sub.origin[1] + y))
                            * self.dims[0]
                            + (sub.origin[0] + x);
                        out[gi] = fields.f[q * np + g.pidx(x, y, z)];
                    }
                }
            }
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::backend::CpuScalar;
    use crate::halo::{InProcess, LocalPeriodic};
    use crate::lattice::D2Q9;

    /// A-5 (E4): building a single-part owner of a wider decomposition with a
    /// Local exchange must fail at construction — such an owner keeps global
    /// neighbour ids that a Local exchange would resolve as local indices
    /// (E4: part=1 of [2,1,1] periodic-x + LocalPeriodic ran without panic
    /// and diverged from the correct 2-part result by up to 7.7e-2).
    #[test]
    #[should_panic(expected = "Remote halo exchange")]
    fn single_part_owner_rejects_local_periodic_exchange() {
        let spec = GlobalSpec::<f64> {
            dims: [8, 4, 1],
            // Both in-plane axes periodic so the config itself is valid (A-4);
            // the failure under test is the halo *scope*, not coverage.
            periodic: [true, true, false],
            ..Default::default()
        };
        // part=1 of a [2,1,1] decomposition, LocalPeriodic (a Local scope).
        let _s: Solver<D2Q9, f64, CpuScalar, LocalPeriodic> = Solver::new_local_part(
            &spec,
            &[],
            &[],
            [2, 1, 1],
            1,
            CpuScalar::default(),
            LocalPeriodic,
        );
    }

    /// The same misuse with `InProcess` is equally rejected (also Local).
    #[test]
    #[should_panic(expected = "Remote halo exchange")]
    fn single_part_owner_rejects_in_process_exchange() {
        let spec = GlobalSpec::<f64> {
            dims: [8, 4, 1],
            // Both in-plane axes periodic so the config itself is valid (A-4);
            // the failure under test is the halo *scope*, not coverage.
            periodic: [true, true, false],
            ..Default::default()
        };
        let _s: Solver<D2Q9, f64, CpuScalar, InProcess> = Solver::new_local_part(
            &spec,
            &[],
            &[],
            [2, 1, 1],
            0,
            CpuScalar::default(),
            InProcess,
        );
    }

    /// A full in-process decomposition (owns every part) is the legitimate
    /// Local use and must still build.
    #[test]
    fn full_in_process_decomposition_builds() {
        let spec = GlobalSpec::<f64> {
            dims: [8, 4, 1],
            // Both in-plane axes periodic so the config itself is valid (A-4);
            // the failure under test is the halo *scope*, not coverage.
            periodic: [true, true, false],
            ..Default::default()
        };
        let mut s: Solver<D2Q9, f64, CpuScalar, InProcess> =
            Solver::new(&spec, &[], &[], [2, 1, 1], CpuScalar::default(), InProcess);
        s.run(2);
        assert!(s.total_mass().is_finite());
    }

    // ----------------------------------------------------------------------
    // A-4: GlobalSpec::validate
    // ----------------------------------------------------------------------

    use crate::lattice::D3Q19;
    use crate::params::FaceBC;

    /// Full solid rims for a walled non-periodic D3Q19 box (so a "closed
    /// non-periodic face" is legitimately covered in the positive tests).
    fn walled_box_solid(dims: [usize; 3]) -> Vec<bool> {
        let mut walls = WallSpec::<f64>::default();
        for f in Face::ALL {
            walls.is_wall[f.index()] = true;
        }
        build_wall_rims(3, dims, &walls).0
    }

    /// E2: a non-periodic z-face that is neither open nor a solid rim is
    /// rejected (its halo would feed stale interior values every step —
    /// nonfinite=0 yet mass drift 2.7e-3, false uz 2.6e-3).
    #[test]
    fn validate_rejects_uncovered_face() {
        // z non-periodic, no z walls, no z open BC, no solids → uncovered.
        let spec = GlobalSpec::<f64> {
            dims: [6, 6, 6],
            periodic: [true, true, false],
            ..Default::default()
        };
        assert!(matches!(
            spec.validate(3, &[]),
            Err(SpecError::UncoveredFace { .. })
        ));
        // Covered by a full z-wall rim → OK.
        let mut walls = WallSpec::<f64>::default();
        walls.is_wall[Face::ZNeg.index()] = true;
        walls.is_wall[Face::ZPos.index()] = true;
        let (solid, _) = build_wall_rims(3, spec.dims, &walls);
        assert!(spec.validate(3, &solid).is_ok());
    }

    /// E3: ν = 0 (and non-finite ν) are rejected (omega_m collapses to 0).
    #[test]
    fn validate_rejects_bad_viscosity() {
        let dims = [6, 6, 6];
        let solid = walled_box_solid(dims);
        let zero_nu = GlobalSpec::<f64> {
            dims,
            nu: 0.0,
            periodic: [false, false, false],
            ..Default::default()
        };
        assert!(matches!(
            zero_nu.validate(3, &solid),
            Err(SpecError::NonPositiveViscosity { .. })
        ));
        let nan_nu = GlobalSpec::<f64> {
            dims,
            nu: f64::NAN,
            periodic: [false, false, false],
            ..Default::default()
        };
        assert!(matches!(
            nan_nu.validate(3, &solid),
            Err(SpecError::NonFiniteParameter { .. })
        ));
    }

    /// periodic × open on the same axis is rejected.
    #[test]
    fn validate_rejects_periodic_open_conflict() {
        let mut faces = [FaceBC::<f64>::Closed; 6];
        faces[Face::XNeg.index()] = FaceBC::Velocity {
            u: [0.05, 0.0, 0.0],
        };
        let spec = GlobalSpec::<f64> {
            dims: [6, 6, 6],
            periodic: [true, false, false], // x periodic AND x-open
            faces,
            ..Default::default()
        };
        assert!(matches!(
            spec.validate(3, &walled_box_solid([6, 6, 6])),
            Err(SpecError::PeriodicOpenConflict { axis: 0 })
        ));
    }

    /// Open faces on two different axes are rejected (Zou–He edge sharing).
    #[test]
    fn validate_rejects_open_on_multiple_axes() {
        let mut faces = [FaceBC::<f64>::Closed; 6];
        faces[Face::XNeg.index()] = FaceBC::Velocity {
            u: [0.05, 0.0, 0.0],
        };
        faces[Face::YNeg.index()] = FaceBC::Outflow;
        let spec = GlobalSpec::<f64> {
            dims: [6, 6, 6],
            periodic: [false, false, false],
            faces,
            ..Default::default()
        };
        // The remaining closed faces are walled so only the multi-axis rule
        // fires.
        assert!(matches!(
            spec.validate(3, &walled_box_solid([6, 6, 6])),
            Err(SpecError::OpenFacesOnMultipleAxes)
        ));
    }

    /// Out-of-range open-face BC parameters are rejected (NaN-safe speed,
    /// non-positive outlet ρ, convective u_conv ∉ (0,1]).
    #[test]
    fn validate_rejects_bad_face_bc_parameters() {
        let dims = [6, 6, 6];
        let base = |faces| GlobalSpec::<f64> {
            dims,
            periodic: [false, true, true],
            faces,
            ..Default::default()
        };
        // Only x is non-periodic here, so the x-faces carry the open BC and
        // the y/z axes are periodic (covered). Too-fast inlet:
        let mut f = [FaceBC::<f64>::Closed; 6];
        f[Face::XNeg.index()] = FaceBC::Velocity { u: [0.9, 0.0, 0.0] };
        f[Face::XPos.index()] = FaceBC::Outflow;
        assert!(matches!(
            base(f).validate(3, &[]),
            Err(SpecError::VelocityTooHigh { .. })
        ));
        // NaN inlet component (NaN-safe rejection).
        let mut f = [FaceBC::<f64>::Closed; 6];
        f[Face::XNeg.index()] = FaceBC::Velocity {
            u: [f64::NAN, 0.0, 0.0],
        };
        f[Face::XPos.index()] = FaceBC::Outflow;
        assert!(matches!(
            base(f).validate(3, &[]),
            Err(SpecError::NonFiniteParameter { .. })
        ));
        // Non-positive outlet density.
        let mut f = [FaceBC::<f64>::Closed; 6];
        f[Face::XNeg.index()] = FaceBC::Velocity {
            u: [0.05, 0.0, 0.0],
        };
        f[Face::XPos.index()] = FaceBC::Pressure { rho: 0.0 };
        assert!(matches!(
            base(f).validate(3, &[]),
            Err(SpecError::NonPositiveDensity { .. })
        ));
        // Convective u_conv out of (0, 1].
        let mut f = [FaceBC::<f64>::Closed; 6];
        f[Face::XNeg.index()] = FaceBC::Velocity {
            u: [0.05, 0.0, 0.0],
        };
        f[Face::XPos.index()] = FaceBC::Convective { u_conv: 1.5 };
        assert!(matches!(
            base(f).validate(3, &[]),
            Err(SpecError::InvalidConvectiveSpeed { .. })
        ));
    }

    /// A 2D lattice must have force[2] == 0; a too-small active axis and a bad
    /// TRT magic are rejected.
    #[test]
    fn validate_rejects_2d_zforce_small_dims_and_magic() {
        // force[2] != 0 on a 2D spec.
        let spec = GlobalSpec::<f64> {
            dims: [8, 8, 1],
            periodic: [true, true, false],
            force: [0.0, 0.0, 1e-6],
            ..Default::default()
        };
        assert!(matches!(
            spec.validate(2, &[]),
            Err(SpecError::NonZeroZForce2D { .. })
        ));
        // 2-cell active axis.
        let tiny = GlobalSpec::<f64> {
            dims: [2, 8, 1],
            periodic: [true, true, false],
            ..Default::default()
        };
        assert!(matches!(
            tiny.validate(2, &[]),
            Err(SpecError::DomainTooSmall { .. })
        ));
        // Non-positive TRT magic.
        let bad_magic = GlobalSpec::<f64> {
            dims: [8, 8, 1],
            periodic: [true, true, false],
            collision: CollisionKind::Trt { magic: -1.0 },
            ..Default::default()
        };
        assert!(matches!(
            bad_magic.validate(2, &[]),
            Err(SpecError::InvalidMagic { .. })
        ));
    }

    /// A fully-periodic box (no faces to cover) and a fully-walled box both
    /// validate — the legitimate configurations must not be rejected.
    #[test]
    fn validate_accepts_periodic_and_walled() {
        let periodic = GlobalSpec::<f64> {
            dims: [6, 6, 6],
            periodic: [true, true, true],
            ..Default::default()
        };
        assert!(periodic.validate(3, &[]).is_ok());

        let dims = [6, 6, 6];
        let walled = GlobalSpec::<f64> {
            dims,
            periodic: [false, false, false],
            ..Default::default()
        };
        assert!(walled.validate(3, &walled_box_solid(dims)).is_ok());
    }

    /// The internal build-time guard fires for an uncovered native spec even
    /// when a caller bypasses the scenario layer (defense in depth).
    #[test]
    #[should_panic(expected = "invalid GlobalSpec")]
    fn build_panics_on_uncovered_face() {
        let spec = GlobalSpec::<f64> {
            dims: [6, 6, 6],
            periodic: [true, true, false],
            ..Default::default()
        };
        let _s: Solver<D3Q19, f64, CpuScalar, LocalPeriodic> = Solver::new(
            &spec,
            &[],
            &[],
            [1, 1, 1],
            CpuScalar::default(),
            LocalPeriodic,
        );
    }
}
