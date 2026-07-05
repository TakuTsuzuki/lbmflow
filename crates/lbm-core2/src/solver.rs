//! Solver orchestrator: drives subdomains × backend × halo exchange through
//! the V1 step sequence (collide → stream → open faces → moments).
//!
//! The step-phase order, the diagnostics' f64 accumulation and the
//! initialisation paths reproduce V1 `Simulation` mechanics exactly; the
//! compat facade is a thin wrapper over this type with a monolithic (1×1×1)
//! decomposition.

use crate::backend::{Backend, CellRange};
use crate::fields::SoaFields;
use crate::halo::HaloExchange;
use crate::kernels::equilibrium;
use crate::lattice::{Face, Lattice};
use crate::params::{CollisionKind, FaceBC, KParams, Reduction, StepParams};
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
        assert!(
            decomp[a] <= dims[a],
            "more parts than cells on axis {a}"
        );
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
    B: Backend<L, T, Fields = SoaFields<T>>,
    H: HaloExchange<T>,
{
    params: StepParams<T>,
    nu: f64,
    dims: [usize; 3],
    periodic: [bool; 3],
    subs: Vec<Subdomain>,
    parts: Vec<SoaFields<T>>,
    backend: B,
    exchange: H,
    time: u64,
    probed_force: [T; 3],
    masks_dirty: bool,
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
    B: Backend<L, T, Fields = SoaFields<T>>,
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
            subs = vec![subs[part].clone()];
        }
        let mut parts: Vec<SoaFields<T>> = subs.iter().map(|s| backend.alloc(s)).collect();
        // Distribute the global masks into the parts' padded cores.
        for (sub, fields) in subs.iter().zip(parts.iter_mut()) {
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
        let mut solver = Self {
            params,
            nu: spec.nu,
            dims: spec.dims,
            periodic: spec.periodic,
            subs,
            parts,
            backend,
            exchange,
            time: 0,
            probed_force: [T::zero(); 3],
            masks_dirty: true,
            two_pass: false,
            _lattice: std::marker::PhantomData,
        };
        solver.sync_masks();
        // V1 from_config ends with update_moments (u(t=0) = force/2 on fluid).
        for i in 0..solver.parts.len() {
            solver
                .backend
                .update_moments(&solver.subs[i], &mut solver.parts[i], &solver.params);
        }
        solver
    }

    fn sync_masks(&mut self) {
        self.exchange.exchange_masks(&self.subs, &mut self.parts);
        self.masks_dirty = false;
    }

    fn sync_masks_if_dirty(&mut self) {
        if self.masks_dirty {
            self.sync_masks();
        }
    }

    /// Advance one time step (V1 `step` order: collide → stream → swap →
    /// open faces → moments).
    pub fn step(&mut self) {
        self.sync_masks_if_dirty();
        for i in 0..self.parts.len() {
            self.backend
                .collide(&self.subs[i], &mut self.parts[i], &self.params);
        }
        self.exchange.exchange_f::<L>(&self.subs, &mut self.parts);
        let mut pf = [T::zero(); 3];
        for i in 0..self.parts.len() {
            let part_pf = self.stream_part(i);
            pf = [pf[0] + part_pf[0], pf[1] + part_pf[1], pf[2] + part_pf[2]];
        }
        for i in 0..self.parts.len() {
            self.backend.swap(&mut self.parts[i]);
        }
        self.probed_force = pf;
        for i in 0..self.parts.len() {
            self.backend
                .apply_open_faces(&self.subs[i], &mut self.parts[i], &self.params);
        }
        for i in 0..self.parts.len() {
            self.backend
                .update_moments(&self.subs[i], &mut self.parts[i], &self.params);
        }
        self.time += 1;
    }

    fn stream_part(&mut self, i: usize) -> [T; 3] {
        let sub = &self.subs[i];
        if !self.two_pass {
            return self
                .backend
                .stream(sub, &mut self.parts[i], &self.params, CellRange::full(sub));
        }
        // Interior pass first (would overlap an async exchange), then the
        // one-cell boundary shell. Field results are identical to the full
        // pass; only the probe partials' summation order differs.
        let c = sub.geom.core;
        let interior = CellRange {
            lo: [1, 1, if sub.geom.d == 3 { 1 } else { 0 }],
            hi: [
                c[0].saturating_sub(1),
                c[1].saturating_sub(1),
                if sub.geom.d == 3 { c[2].saturating_sub(1) } else { c[2] },
            ],
        };
        let sub = sub.clone();
        let mut pf = self
            .backend
            .stream(&sub, &mut self.parts[i], &self.params, interior);
        for shell in boundary_shells(&sub, interior) {
            let p2 = self
                .backend
                .stream(&sub, &mut self.parts[i], &self.params, shell);
            pf = [pf[0] + p2[0], pf[1] + p2[1], pf[2] + p2[2]];
        }
        pf
    }

    /// Advance `steps` time steps.
    pub fn run(&mut self, steps: usize) {
        for _ in 0..steps {
            self.step();
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
        self.sync_masks_if_dirty();
        let kp = KParams::new::<L>(&self.params);
        let tau = T::r(3.0 * self.nu + 0.5);
        let three = T::r(3.0);
        let half = T::r(0.5);
        let dims = self.dims;
        let periodic = self.periodic;
        for (sub, fields) in self.subs.iter().zip(self.parts.iter_mut()) {
            let g = sub.geom;
            let np = g.n_padded();
            // Pass 1: store the macroscopic fields (all core cells).
            for z in 0..g.core[2] {
                for y in 0..g.core[1] {
                    for x in 0..g.core[0] {
                        let (r, u) = init(
                            sub.origin[0] + x,
                            sub.origin[1] + y,
                            sub.origin[2] + z,
                        );
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
                                        gpos[a] =
                                            (gpos[a] + dims[a] as isize) % dims[a] as isize;
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
                            let (_, u) =
                                init(gpos[0] as usize, gpos[1] as usize, gpos[2] as usize);
                            Some(u)
                        };
                        let own = [fields.ux[c], fields.uy[c], fields.uz[c]];
                        let diff = |plus: Option<[T; 3]>,
                                    minus: Option<[T; 3]>,
                                    b: usize|
                         -> T {
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
                                        ccgu = ccgu
                                            + cq[a] * cq[b] * (grad[a][b] + grad[b][a]);
                                    }
                                }
                            }
                            let fneq =
                                -kp.wr[q] * fields.rho[c] * tau * (three * ccgu - div);
                            fields.f[q * np + pi] = fields.f[q * np + pi] + fneq;
                        }
                    }
                }
            }
        }
        for i in 0..self.parts.len() {
            self.backend
                .update_moments(&self.subs[i], &mut self.parts[i], &self.params);
        }
    }

    /// Mark a global cell solid (half-way bounce-back obstacle). Open-face
    /// checks are the caller's (facade's) responsibility.
    pub fn set_solid(&mut self, x: usize, y: usize, z: usize) {
        let (i, lx, ly, lz) = self.locate(x, y, z);
        let pi = self.subs[i].geom.pidx(lx, ly, lz);
        self.parts[i].solid[pi] = true;
        self.masks_dirty = true;
    }

    /// Select the solid cells whose momentum-exchange force is accumulated
    /// each step (V1 `set_force_probe`).
    pub fn set_force_probe(&mut self, pred: impl Fn(usize, usize, usize) -> bool) {
        for (sub, fields) in self.subs.iter().zip(self.parts.iter_mut()) {
            let g = sub.geom;
            let mut mask = vec![false; g.n_padded()];
            for z in 0..g.core[2] {
                for y in 0..g.core[1] {
                    for x in 0..g.core[0] {
                        mask[g.pidx(x, y, z)] = pred(
                            sub.origin[0] + x,
                            sub.origin[1] + y,
                            sub.origin[2] + z,
                        );
                    }
                }
            }
            fields.probe = Some(mask);
        }
        self.masks_dirty = true;
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
        for (sub, fields) in self.subs.iter().zip(self.parts.iter_mut()) {
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
    }

    /// Closure form of [`Solver::set_inlet_profile`]: `profile(c1, c2)` is
    /// evaluated at the global tangent coordinates of every face node
    /// (`(t1, t2) = face.tangents()`; 2D faces always pass `c2 = 0`).
    /// The natural way to build e.g. a rectangular-duct profile
    /// `u(y, z) = umax f(y) g(z)` on an X face.
    pub fn set_inlet_profile_with(
        &mut self,
        face: Face,
        profile: impl Fn(usize, usize) -> [T; 3],
    ) {
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
        self.sync_masks_if_dirty();
        // ψ planes, padded: core = ψ(rho) (0 on solids), halo = 0 until the
        // exchange fills it (stays 0 outside non-periodic global edges,
        // matching V1's "out-of-domain contributes nothing").
        let mut planes: Vec<Vec<T>> = self
            .subs
            .iter()
            .zip(self.parts.iter())
            .map(|(sub, fields)| {
                let geo = sub.geom;
                let mut plane = vec![T::zero(); geo.n_padded()];
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
                plane
            })
            .collect();
        let mut refs: Vec<&mut [T]> = planes.iter_mut().map(|p| p.as_mut_slice()).collect();
        self.exchange.exchange_scalar(&self.subs, &mut refs);
        for (i, (sub, fields)) in self.subs.iter().zip(self.parts.iter_mut()).enumerate() {
            let geo = sub.geom;
            let plane = &planes[i];
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
                        for q in 1..L::Q {
                            let cq = L::C[q];
                            let pj = plane[geo.pidx_i(
                                x as isize + cq[0] as isize,
                                y as isize + cq[1] as isize,
                                z as isize + cq[2] as isize,
                            )];
                            let w = T::r(L::W[q]);
                            for a in 0..L::D {
                                s[a] = s[a] + w * pj * T::r(cq[a] as f64);
                            }
                        }
                        for a in 0..3 {
                            ff[c][a] = -(psi_i * (g * s[a]));
                        }
                    }
                }
            }
        }
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
    /// Subdomain descriptor `i`.
    pub fn sub(&self, i: usize) -> &Subdomain {
        &self.subs[i]
    }
    /// Fields of part `i` (host staging; padded mask edits must go through
    /// `set_solid` / `set_force_probe` so halos stay in sync).
    pub fn fields(&self, i: usize) -> &SoaFields<T> {
        &self.parts[i]
    }
    /// Mutable fields of part `i` (see [`Solver::fields`] caveat).
    pub fn fields_mut(&mut self, i: usize) -> &mut SoaFields<T> {
        &mut self.parts[i]
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
        self.parts[i].rho[g.cidx(lx, ly, lz)]
    }
    /// Velocity at a global cell (physical, half-force corrected).
    pub fn u(&self, x: usize, y: usize, z: usize) -> [T; 3] {
        let (i, lx, ly, lz) = self.locate(x, y, z);
        let g = self.subs[i].geom;
        let c = g.cidx(lx, ly, lz);
        [self.parts[i].ux[c], self.parts[i].uy[c], self.parts[i].uz[c]]
    }
    /// Whether a global cell is solid.
    pub fn is_solid(&self, x: usize, y: usize, z: usize) -> bool {
        let (i, lx, ly, lz) = self.locate(x, y, z);
        self.parts[i].solid[self.subs[i].geom.pidx(lx, ly, lz)]
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
        p
    }

    /// Number of non-finite (NaN/Inf) values in this process's parts, over
    /// the populations and the macroscopic moments. `0` on a healthy run;
    /// a distributed owner sums the counts across ranks.
    pub fn local_nonfinite_count(&self) -> u64 {
        let mut n = 0u64;
        for fields in &self.parts {
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

    /// Assemble a global compact array from a per-part compact getter
    /// (test/diagnostic helper).
    fn gather(&self, get: impl Fn(&SoaFields<T>, usize) -> T) -> Vec<T> {
        let mut out = vec![T::zero(); self.dims[0] * self.dims[1] * self.dims[2]];
        for (sub, fields) in self.subs.iter().zip(self.parts.iter()) {
            let g = sub.geom;
            for z in 0..g.core[2] {
                for y in 0..g.core[1] {
                    for x in 0..g.core[0] {
                        let gi = ((sub.origin[2] + z) * self.dims[1] + (sub.origin[1] + y))
                            * self.dims[0]
                            + (sub.origin[0] + x);
                        out[gi] = get(fields, g.cidx(x, y, z));
                    }
                }
            }
        }
        out
    }

    /// Global density field (compact layout).
    pub fn gather_rho(&self) -> Vec<T> {
        self.gather(|f, c| f.rho[c])
    }
    /// Global x-velocity field.
    pub fn gather_ux(&self) -> Vec<T> {
        self.gather(|f, c| f.ux[c])
    }
    /// Global y-velocity field.
    pub fn gather_uy(&self) -> Vec<T> {
        self.gather(|f, c| f.uy[c])
    }
    /// Global z-velocity field.
    pub fn gather_uz(&self) -> Vec<T> {
        self.gather(|f, c| f.uz[c])
    }
    /// Global deviation-population plane `q` (compact layout).
    pub fn gather_f(&self, q: usize) -> Vec<T> {
        let mut out = vec![T::zero(); self.dims[0] * self.dims[1] * self.dims[2]];
        for (sub, fields) in self.subs.iter().zip(self.parts.iter()) {
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

/// Boundary shells complementing the interior box: fixed order YNeg row,
/// YPos row, XNeg column, XPos column (minus corners already covered), then
/// Z planes for 3D. Only the probe partials' summation order depends on
/// this; field results do not.
fn boundary_shells(sub: &Subdomain, interior: CellRange) -> Vec<CellRange> {
    let c = sub.geom.core;
    let mut shells = Vec::new();
    let z_full = (0, c[2]);
    // y = 0 and y = ny-1 full-width rows.
    shells.push(CellRange {
        lo: [0, 0, z_full.0],
        hi: [c[0], interior.lo[1], z_full.1],
    });
    shells.push(CellRange {
        lo: [0, interior.hi[1], z_full.0],
        hi: [c[0], c[1], z_full.1],
    });
    // x columns between the y rows.
    shells.push(CellRange {
        lo: [0, interior.lo[1], z_full.0],
        hi: [interior.lo[0], interior.hi[1], z_full.1],
    });
    shells.push(CellRange {
        lo: [interior.hi[0], interior.lo[1], z_full.0],
        hi: [c[0], interior.hi[1], z_full.1],
    });
    if sub.geom.d == 3 {
        // z planes of the remaining interior-xy box.
        shells.push(CellRange {
            lo: [interior.lo[0], interior.lo[1], 0],
            hi: [interior.hi[0], interior.hi[1], interior.lo[2]],
        });
        shells.push(CellRange {
            lo: [interior.lo[0], interior.lo[1], interior.hi[2]],
            hi: [interior.hi[0], interior.hi[1], c[2]],
        });
    }
    shells.retain(|s| !s.is_empty());
    shells
}
