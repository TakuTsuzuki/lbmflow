//! SoA field storage over a halo-padded local box.
//!
//! Populations are stored in deviation form (`f_q - w_q`, V1 convention: the
//! quiescent state is exactly all-zero, which is what makes f32 storage
//! validation-grade) as q-major planes:
//!
//! ```text
//! f[q * n_padded + cell],   cell = z * (pnx*pny) + y * pnx + x   (padded coords)
//! ```
//!
//! For an unpadded box (`halo = 0`) this is exactly the storage contract of
//! docs/ARCHITECTURE_V2.md §2.2, `cell = z·(nx·ny) + y·nx + x`. Subdomains
//! (§2.3) pad every axis `< D` by `halo` (= 1) ghost cells on each side; the
//! same formula then runs over the padded extents. GPU kernels use the
//! identical layout, so host/device copies are plain memcpys per plane.
//!
//! Two categories of arrays:
//! - **padded** (`f`, `ftmp`, `solid`, `wall_u`, `probe`): read across
//!   subdomain boundaries by streaming, so they carry the halo.
//! - **compact** (`rho`, `ux/uy/uz`, `force_field`): only ever accessed at
//!   core cells; kept halo-free so they can be borrowed directly as V1-shaped
//!   `&[T]` fields (explicit-readback boundary for GPU backends later).

use crate::real::Real;

/// Geometry of one local box: core extents plus halo width.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct LocalGeom {
    /// Spatial dimension (2 or 3).
    pub d: usize,
    /// Core cell extents `[nx, ny, nz]`; `nz == 1` when `d == 2`.
    pub core: [usize; 3],
    /// Ghost-cell width added on both sides of every axis `< d`.
    pub halo: usize,
}

impl LocalGeom {
    /// New geometry. Axes `>= d` must have extent 1 and get no halo.
    pub fn new(d: usize, core: [usize; 3], halo: usize) -> Self {
        assert!(d == 2 || d == 3, "dimension must be 2 or 3");
        for (a, &n) in core.iter().enumerate() {
            if a < d {
                assert!(n >= 1, "axis {a} extent must be >= 1");
            } else {
                assert_eq!(n, 1, "axis {a} extent must be 1 for d = {d}");
            }
        }
        Self { d, core, halo }
    }

    /// Padded extents (core plus halo on active axes).
    #[inline]
    pub fn padded(&self) -> [usize; 3] {
        let mut p = self.core;
        for (a, e) in p.iter_mut().enumerate() {
            if a < self.d {
                *e += 2 * self.halo;
            }
        }
        p
    }

    /// Number of cells in the padded box (one storage plane).
    #[inline]
    pub fn n_padded(&self) -> usize {
        let p = self.padded();
        p[0] * p[1] * p[2]
    }

    /// Number of core cells.
    #[inline]
    pub fn n_core(&self) -> usize {
        self.core[0] * self.core[1] * self.core[2]
    }

    /// Padded-array index of a core cell given in core coordinates.
    #[inline(always)]
    pub fn pidx(&self, x: usize, y: usize, z: usize) -> usize {
        debug_assert!(x < self.core[0] && y < self.core[1] && z < self.core[2]);
        self.pidx_i(x as isize, y as isize, z as isize)
    }

    /// Padded-array index of a cell in core coordinates that may reach into
    /// the halo (each coordinate in `-halo .. core + halo`).
    #[inline(always)]
    pub fn pidx_i(&self, x: isize, y: isize, z: isize) -> usize {
        let h = self.halo as isize;
        let p = self.padded();
        let (hx, hy) = (h, h);
        let hz = if self.d == 3 { h } else { 0 };
        debug_assert!(x >= -hx && x < self.core[0] as isize + hx);
        debug_assert!(y >= -hy && y < self.core[1] as isize + hy);
        debug_assert!(z >= -hz && z < self.core[2] as isize + hz);
        ((z + hz) as usize * p[1] + (y + hy) as usize) * p[0] + (x + hx) as usize
    }

    /// Compact (halo-free) index of a core cell: `z*(nx*ny) + y*nx + x`.
    #[inline(always)]
    pub fn cidx(&self, x: usize, y: usize, z: usize) -> usize {
        debug_assert!(x < self.core[0] && y < self.core[1] && z < self.core[2]);
        (z * self.core[1] + y) * self.core[0] + x
    }
}

/// Scratch state of the fused `CpuSimd` backend, allocated lazily on the
/// first fused pass and ignored by `CpuScalar`.
///
/// - `rho2/ux2/uy2/uz2` are the moment double buffers (V1 `rho2` mechanics):
///   the fused pass writes step k's moments here while in-flight collides of
///   other rows still read step k-1's moments from the primary buffers; the
///   pairs are swapped in `Backend::swap`. Solid cells are refreshed from the
///   primary buffers every pass, so both buffers stay in sync at cells the
///   pass never computes (V1 `set_solid`/`init_with` invariant).
/// - `stale` reproduces V1's ping-pong stale-slot convention for open faces:
///   `CpuScalar` collides in place, so the unknown slots streaming skips
///   retain the *previous step's post-collide* populations, which the open
///   BCs (Convective in particular) read. The fused pass never materialises
///   post-collide state in `f`, so it captures those values per open face
///   (`stale[1]`, this step) and `apply_open_faces` writes the previous
///   step's capture (`stale[0]`) back into the unknown slots before the BC
///   pass, then swaps the pair. Zero-initialised: matches the all-zero
///   deviation ping-pong buffer a `CpuScalar` run reads on its first step.
/// - `fresh` marks that the fused pass has just written the double buffers,
///   so the next `update_moments` only needs the open-face boundary fix
///   (V1 `fix_boundary_moments`) instead of a full recompute.
#[derive(Clone, Debug)]
pub struct FusedScratch<T: Real> {
    /// Density double buffer, compact core layout.
    pub rho2: Vec<T>,
    /// x-velocity double buffer, compact core.
    pub ux2: Vec<T>,
    /// y-velocity double buffer, compact core.
    pub uy2: Vec<T>,
    /// z-velocity double buffer, compact core (untouched for 2D lattices).
    pub uz2: Vec<T>,
    /// Stale-slot stash per open face (`Face::index()` order), cell-major
    /// over the face's canonical cell order, unknown directions innermost.
    /// `stale[0]` feeds this step's BC pass, `stale[1]` is captured during
    /// this step's fused pass; swapped after the BC pass.
    pub stale: [[Option<Vec<T>>; 6]; 2],
    /// Set by the fused stream pass, consumed by `update_moments`.
    pub fresh: bool,
    /// Per-band source rings, reused across steps (contents are transient
    /// per-pass scratch; kept here so the buffers are neither reallocated
    /// nor re-zeroed every step).
    pub(crate) rings: Vec<crate::backend_simd::Ring<T>>,
}

impl<T: Real> FusedScratch<T> {
    /// Allocate for `n_core` cells, moments copied from the primary buffers'
    /// quiescent state semantics (all zero; the first fused pass rewrites
    /// every cell before the buffers are ever swapped in).
    pub fn new(n_core: usize) -> Self {
        Self {
            rho2: vec![T::zero(); n_core],
            ux2: vec![T::zero(); n_core],
            uy2: vec![T::zero(); n_core],
            uz2: vec![T::zero(); n_core],
            stale: Default::default(),
            fresh: false,
            rings: Vec::new(),
        }
    }
}

/// SoA fields of one subdomain (the `CpuScalar` backend's field storage).
#[derive(Clone, Debug)]
pub struct SoaFields<T: Real> {
    /// Local geometry (shared by every array below).
    pub geom: LocalGeom,
    /// Number of populations (`L::Q`).
    pub q: usize,
    /// Deviation populations, q-major padded planes (current state).
    pub f: Vec<T>,
    /// Ping-pong partner of `f`. Streaming writes here, then the pair is
    /// swapped. Unknown slots skipped by streaming retain this buffer's prior
    /// content — the ConvectiveOutflow BC depends on that (V1 mechanics).
    pub ftmp: Vec<T>,
    /// Conservative Allen-Cahn phase-field distribution set, q-major padded
    /// D3Q19 planes in ordinary (non-deviation) form. `None` means single
    /// phase; the hydrodynamic path must remain bit-identical in that state.
    pub g: Option<Vec<T>>,
    /// Ping-pong partner of `g`.
    pub gtmp: Option<Vec<T>>,
    /// Phase field `phi = sum_i g_i`, compact core layout.
    pub phi: Option<Vec<T>>,
    /// Density, compact core. `1` on quiescent build; moments skip solids.
    pub rho: Vec<T>,
    /// x-velocity (physical: includes the Guo half-force term), compact core.
    pub ux: Vec<T>,
    /// y-velocity, compact core.
    pub uy: Vec<T>,
    /// z-velocity, compact core (all zero for 2D lattices).
    pub uz: Vec<T>,
    /// Solid mask, padded.
    pub solid: Vec<bool>,
    /// Wall velocity per cell (meaningful on solids only), padded.
    pub wall_u: Vec<[T; 3]>,
    /// Momentum-exchange probe mask over solids, padded.
    pub probe: Option<Vec<bool>>,
    /// Momentum-exchange force accumulated during the most recent step.
    pub probed_force: [T; 3],
    /// Per-cell body force added to the uniform force, compact core.
    pub force_field: Option<Vec<[T; 3]>>,
    /// Optional per-cell symmetric relaxation rate (`omega_plus = 1/tau`),
    /// compact core. `None` preserves the uniform-omega collide path exactly.
    pub omega_field: Option<Vec<T>>,
    /// Optional Bouzidi curved-wall links. `None` preserves the half-way
    /// bounce-back path exactly and allocates no record storage.
    pub bouzidi: Option<crate::bouzidi::BouzidiLinks<T>>,
    /// Per-face inlet velocity profiles ([`crate::lattice::Face::index`]
    /// order), overriding a `Velocity` face's uniform value. Indexed by this
    /// part's local along-face coordinate (2D: the single tangent axis;
    /// 3D: `t2 * extent(t1) + t1`, tangent axes ascending).
    pub inlet_profiles: [Option<Vec<[T; 3]>>; 6],
    /// Fused-backend scratch (`CpuSimd`); `None` under `CpuScalar`.
    pub fused: Option<Box<FusedScratch<T>>>,
}

impl<T: Real> SoaFields<T> {
    /// Allocate a quiescent state: `f = 0` (deviation form ⇒ rho = 1, u = 0),
    /// no solids, no probe, no per-cell force.
    pub fn new(q: usize, geom: LocalGeom) -> Self {
        let np = geom.n_padded();
        let nc = geom.n_core();
        Self {
            geom,
            q,
            f: vec![T::zero(); q * np],
            ftmp: vec![T::zero(); q * np],
            g: None,
            gtmp: None,
            phi: None,
            rho: vec![T::one(); nc],
            ux: vec![T::zero(); nc],
            uy: vec![T::zero(); nc],
            uz: vec![T::zero(); nc],
            solid: vec![false; np],
            wall_u: vec![[T::zero(); 3]; np],
            probe: None,
            probed_force: [T::zero(); 3],
            force_field: None,
            omega_field: None,
            bouzidi: None,
            inlet_profiles: [None, None, None, None, None, None],
            fused: None,
        }
    }

    /// One storage plane's length (padded cell count).
    #[inline]
    pub fn plane_len(&self) -> usize {
        self.geom.n_padded()
    }

    /// Swap the population ping-pong pair (after streaming).
    #[inline]
    pub fn swap_f(&mut self) {
        std::mem::swap(&mut self.f, &mut self.ftmp);
    }

    /// Plane `q` of the current populations.
    #[inline]
    pub fn f_plane(&self, q: usize) -> &[T] {
        let n = self.plane_len();
        &self.f[q * n..(q + 1) * n]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn unpadded_index_matches_spec_formula() {
        // cell = z*(nx*ny) + y*nx + x per ARCHITECTURE_V2 §2.2.
        let g = LocalGeom::new(3, [4, 3, 2], 0);
        let (nx, ny) = (4, 3);
        for z in 0..2 {
            for y in 0..3 {
                for x in 0..4 {
                    assert_eq!(g.pidx(x, y, z), z * (nx * ny) + y * nx + x);
                    assert_eq!(g.cidx(x, y, z), z * (nx * ny) + y * nx + x);
                }
            }
        }
        assert_eq!(g.n_padded(), 24);
        assert_eq!(g.n_core(), 24);
    }

    #[test]
    fn halo_padding_2d() {
        let g = LocalGeom::new(2, [4, 3, 1], 1);
        assert_eq!(g.padded(), [6, 5, 1]);
        assert_eq!(g.n_padded(), 30);
        assert_eq!(g.n_core(), 12);
        // Core (0,0) sits one ring in.
        assert_eq!(g.pidx(0, 0, 0), 6 + 1);
        // Halo corner is index 0.
        assert_eq!(g.pidx_i(-1, -1, 0), 0);
        // z never pads in 2D.
        assert_eq!(g.pidx_i(4, 3, 0), 4 * 6 + 5);
    }

    #[test]
    fn halo_padding_3d() {
        let g = LocalGeom::new(3, [4, 3, 2], 1);
        assert_eq!(g.padded(), [6, 5, 4]);
        assert_eq!(g.pidx_i(-1, -1, -1), 0);
        assert_eq!(g.pidx(0, 0, 0), (6 * 5) + 6 + 1);
    }

    #[test]
    fn fields_allocate_quiescent() {
        let g = LocalGeom::new(2, [8, 4, 1], 1);
        let f: SoaFields<f64> = SoaFields::new(9, g);
        assert_eq!(f.f.len(), 9 * g.n_padded());
        assert_eq!(f.rho.len(), 32);
        assert!(f.rho.iter().all(|&r| r == 1.0));
        assert!(f.f.iter().all(|&v| v == 0.0));
        assert_eq!(f.f_plane(8).len(), g.n_padded());
    }

    #[test]
    #[should_panic(expected = "axis 2 extent must be 1")]
    fn rejects_3d_extent_for_2d() {
        LocalGeom::new(2, [4, 4, 2], 1);
    }
}
