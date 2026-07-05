//! Backend abstraction (docs/ARCHITECTURE_V2.md §2.4) and the `CpuScalar`
//! reference implementation.
//!
//! Design deviation from the §2.4 sketch (documented per the "improvements
//! allowed" clause): `step()` is split into phase methods (`collide` /
//! `stream` / `swap` / `apply_open_faces` / `update_moments`) because the
//! halo exchange interleaves between collision and streaming and is a
//! decomposition concern, not a backend one. The orchestrator (`Solver`)
//! owns the step sequence; `reduce` and `read_moments` keep the sketch's
//! shape (backend-side reductions, explicit readback).
//!
//! Fused-kernel note for future backends (GPU_EVALUATION.md §2): a backend
//! that fuses collide+stream into one pass computes `(C∘S)^k` per k steps,
//! while this reference computes `(S∘C)^k`. The operator identity
//! `(C∘S)^k ∘ C = C ∘ (S∘C)^k` plus collision invariance of (rho, momentum)
//! makes the two comparable: apply one extra collision to the initial state
//! of the fused backend and compare velocity fields 1:1 (frozen as the
//! T14 regression recipe).

use crate::fields::SoaFields;
use crate::kernels::{
    collide_row, convective_face, moments_row, outflow_face, stream_row, zou_he_face, RawSlice,
    ZhKind,
};
use crate::lattice::{Face, Lattice};
use crate::params::{FaceBC, KParams, Reduction, StepParams};
use crate::real::Real;
use crate::subdomain::Subdomain;

#[cfg(feature = "parallel")]
use rayon::prelude::*;

/// Below this many core cells the row-parallel loops run serially: rayon's
/// dispatch overhead dwarfs the work on small grids (V1 threshold, measured
/// on an 18-core M-series).
pub const PARALLEL_MIN_CELLS: usize = 16_384;

/// Host-side copy of the macroscopic fields (explicit readback target).
#[derive(Clone, Debug, Default)]
pub struct HostMoments<T: Real> {
    /// Density, compact core layout.
    pub rho: Vec<T>,
    /// x-velocity (physical: includes the Guo half-force term).
    pub ux: Vec<T>,
    /// y-velocity.
    pub uy: Vec<T>,
    /// z-velocity (zero for 2D lattices).
    pub uz: Vec<T>,
}

/// Inclusive-exclusive core-cell range for split streaming passes
/// (interior first, boundary shell after the halo arrives).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct CellRange {
    /// Low corner (inclusive), core coordinates.
    pub lo: [usize; 3],
    /// High corner (exclusive).
    pub hi: [usize; 3],
}

impl CellRange {
    /// The full core box of a subdomain.
    pub fn full(sub: &Subdomain) -> Self {
        Self {
            lo: [0, 0, 0],
            hi: sub.geom.core,
        }
    }

    /// Whether the range is empty on any axis.
    pub fn is_empty(&self) -> bool {
        (0..3).any(|a| self.lo[a] >= self.hi[a])
    }
}

/// A compute + storage target for one subdomain's fields.
///
/// One time step, orchestrated by the solver, is:
/// `collide` → halo exchange → `stream` (interior, then boundary) → `swap`
/// → `apply_open_faces` → `update_moments`.
pub trait Backend<L: Lattice, T: Real> {
    /// Device-resident field storage.
    type Fields;

    /// Allocate quiescent fields for a subdomain.
    fn alloc(&self, sub: &Subdomain) -> Self::Fields;

    /// TRT/BGK collision with Guo forcing over all core cells (in place).
    fn collide(&mut self, sub: &Subdomain, fields: &mut Self::Fields, p: &StepParams<T>);

    /// Pull-streaming into the out-buffer over `range`. Returns the
    /// momentum-exchange force accumulated over probed solid links in the
    /// range (deterministic row-order sum).
    fn stream(
        &mut self,
        sub: &Subdomain,
        fields: &mut Self::Fields,
        p: &StepParams<T>,
        range: CellRange,
    ) -> [T; 3];

    /// Swap the population ping-pong pair (after all stream ranges ran).
    fn swap(&mut self, fields: &mut Self::Fields);

    /// Open-face BC pass (Zou–He / outflow / convective) on the faces of
    /// this subdomain that lie on an open global face.
    fn apply_open_faces(&mut self, sub: &Subdomain, fields: &mut Self::Fields, p: &StepParams<T>);

    /// Recompute macroscopic moments from the populations.
    fn update_moments(&mut self, sub: &Subdomain, fields: &mut Self::Fields, p: &StepParams<T>);

    /// Backend-side reduction over fluid core cells, accumulated in `f64`
    /// in compact cell order (V1 diagnostic convention).
    fn reduce(
        &self,
        sub: &Subdomain,
        fields: &Self::Fields,
        p: &StepParams<T>,
        kind: Reduction,
    ) -> f64;

    /// Explicit readback of the macroscopic fields into host memory.
    fn read_moments(&self, fields: &Self::Fields, out: &mut HostMoments<T>);
}

/// Correctness-first CPU backend: scalar row kernels, optional rayon
/// row-parallel dispatch (identical per-cell arithmetic either way; the
/// probe reduction is an ordered row sum, so results are deterministic —
/// unlike V1's rayon `reduce`, whose tree shape may vary run to run).
#[derive(Clone, Copy, Debug)]
pub struct CpuScalar {
    /// Row-parallel dispatch threshold in core cells.
    pub parallel_min_cells: usize,
}

impl Default for CpuScalar {
    fn default() -> Self {
        Self {
            parallel_min_cells: PARALLEL_MIN_CELLS,
        }
    }
}

impl CpuScalar {
    fn use_parallel(&self, sub: &Subdomain) -> bool {
        cfg!(feature = "parallel") && sub.geom.n_core() >= self.parallel_min_cells
    }
}

impl<L: Lattice, T: Real> Backend<L, T> for CpuScalar {
    type Fields = SoaFields<T>;

    fn alloc(&self, sub: &Subdomain) -> SoaFields<T> {
        SoaFields::new(L::Q, sub.geom)
    }

    fn collide(&mut self, sub: &Subdomain, fields: &mut SoaFields<T>, p: &StepParams<T>) {
        let kp = KParams::new::<L>(p);
        let g = fields.geom;
        let np = g.n_padded();
        let nx = g.core[0];
        let rows = g.core[1] * g.core[2];
        let f = RawSlice::new(&mut fields.f);
        let (rho, ux, uy, uz) = (&fields.rho, &fields.ux, &fields.uy, &fields.uz);
        let solid = &fields.solid;
        let ff = fields.force_field.as_deref();
        let body = |r: usize| {
            let y = r % g.core[1];
            let z = r / g.core[1];
            let pb = g.pidx(0, y, z);
            let c0 = g.cidx(0, y, z);
            // SAFETY: each row index r is processed exactly once, and
            // collide_row writes only its own row's cells.
            unsafe {
                collide_row::<L, T>(
                    f,
                    np,
                    pb,
                    &rho[c0..c0 + nx],
                    &ux[c0..c0 + nx],
                    &uy[c0..c0 + nx],
                    &uz[c0..c0 + nx],
                    &solid[pb..pb + nx],
                    ff.map(|v| &v[c0..c0 + nx]),
                    &kp,
                )
            }
        };
        #[cfg(feature = "parallel")]
        if self.use_parallel(sub) {
            (0..rows).into_par_iter().for_each(body);
            return;
        }
        (0..rows).for_each(body);
    }

    fn stream(
        &mut self,
        sub: &Subdomain,
        fields: &mut SoaFields<T>,
        p: &StepParams<T>,
        range: CellRange,
    ) -> [T; 3] {
        if range.is_empty() {
            return [T::zero(); 3];
        }
        let kp = KParams::new::<L>(p);
        let g = fields.geom;
        let np = g.n_padded();
        let halo = sub.halo_flags();
        let out = RawSlice::new(&mut fields.ftmp);
        let f = &fields.f;
        let (solid, wall_u, rho) = (&fields.solid, &fields.wall_u, &fields.rho);
        let probe = fields.probe.as_deref();
        let (ny_r, nz_r) = (range.hi[1] - range.lo[1], range.hi[2] - range.lo[2]);
        let rows = ny_r * nz_r;
        let body = |r: usize| -> [T; 3] {
            let y = range.lo[1] + r % ny_r;
            let z = range.lo[2] + r / ny_r;
            let c0 = g.cidx(0, y, z);
            // SAFETY: each row index r is processed exactly once, and
            // stream_row writes only its own row segment's cells.
            unsafe {
                stream_row::<L, T>(
                    out,
                    f,
                    np,
                    &g,
                    halo,
                    y,
                    z,
                    range.lo[0],
                    range.hi[0],
                    solid,
                    wall_u,
                    &rho[c0..c0 + g.core[0]],
                    probe,
                    &kp,
                )
            }
        };
        // Deterministic probe total: fold row partials in row order.
        let fold = |partials: Vec<[T; 3]>| -> [T; 3] {
            partials.into_iter().fold([T::zero(); 3], |a, b| {
                [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
            })
        };
        #[cfg(feature = "parallel")]
        if self.use_parallel(sub) {
            let partials: Vec<[T; 3]> = (0..rows).into_par_iter().map(body).collect();
            return fold(partials);
        }
        fold((0..rows).map(body).collect())
    }

    fn swap(&mut self, fields: &mut SoaFields<T>) {
        fields.swap_f();
    }

    fn apply_open_faces(&mut self, sub: &Subdomain, fields: &mut SoaFields<T>, p: &StepParams<T>) {
        let g = fields.geom;
        let np = g.n_padded();
        for face in Face::ALL {
            if !sub.touches_global_face(face) {
                continue;
            }
            let bc = &p.faces[face.index()];
            if !bc.is_open() {
                continue;
            }
            let profiles = &fields.inlet_profiles;
            match bc {
                FaceBC::Closed => {}
                FaceBC::Velocity { u } => zou_he_face::<L, T>(
                    &mut fields.f,
                    np,
                    &g,
                    &fields.solid,
                    face,
                    &ZhKind::Velocity(*u),
                    profiles[face.index()].as_deref(),
                ),
                FaceBC::Pressure { rho } => zou_he_face::<L, T>(
                    &mut fields.f,
                    np,
                    &g,
                    &fields.solid,
                    face,
                    &ZhKind::Pressure(*rho),
                    None,
                ),
                FaceBC::Outflow => {
                    outflow_face::<L, T>(&mut fields.f, np, &g, &fields.solid, face)
                }
                FaceBC::Convective { u_conv } => {
                    convective_face::<L, T>(&mut fields.f, np, &g, &fields.solid, face, *u_conv)
                }
            }
        }
    }

    fn update_moments(&mut self, sub: &Subdomain, fields: &mut SoaFields<T>, p: &StepParams<T>) {
        let kp = KParams::new::<L>(p);
        let g = fields.geom;
        let np = g.n_padded();
        let nx = g.core[0];
        let f = &fields.f;
        let solid = &fields.solid;
        let ff = fields.force_field.as_deref();
        let body = |(r, ((rrow, uxrow), (uyrow, uzrow))): (
            usize,
            ((&mut [T], &mut [T]), (&mut [T], &mut [T])),
        )| {
            let y = r % g.core[1];
            let z = r / g.core[1];
            let pb = g.pidx(0, y, z);
            let c0 = g.cidx(0, y, z);
            moments_row::<L, T>(
                f,
                np,
                pb,
                rrow,
                uxrow,
                uyrow,
                uzrow,
                &solid[pb..pb + nx],
                ff.map(|v| &v[c0..c0 + nx]),
                &kp,
            );
        };
        #[cfg(feature = "parallel")]
        if self.use_parallel(sub) {
            fields
                .rho
                .par_chunks_mut(nx)
                .zip(fields.ux.par_chunks_mut(nx))
                .zip(
                    fields
                        .uy
                        .par_chunks_mut(nx)
                        .zip(fields.uz.par_chunks_mut(nx)),
                )
                .enumerate()
                .for_each(body);
            return;
        }
        fields
            .rho
            .chunks_mut(nx)
            .zip(fields.ux.chunks_mut(nx))
            .zip(fields.uy.chunks_mut(nx).zip(fields.uz.chunks_mut(nx)))
            .enumerate()
            .for_each(body);
    }

    fn reduce(
        &self,
        sub: &Subdomain,
        fields: &SoaFields<T>,
        p: &StepParams<T>,
        kind: Reduction,
    ) -> f64 {
        let g = fields.geom;
        let np = g.n_padded();
        let ff = fields.force_field.as_deref();
        let mut acc = 0.0f64;
        // Compact cell order (z, y, x ascending), q inner — V1's exact
        // f64 accumulation sequence.
        for z in 0..g.core[2] {
            for y in 0..g.core[1] {
                let pb = g.pidx(0, y, z);
                let c0 = g.cidx(0, y, z);
                for x in 0..g.core[0] {
                    if fields.solid[pb + x] {
                        continue;
                    }
                    match kind {
                        Reduction::FluidCells => acc += 1.0,
                        Reduction::MassDeviation => {
                            for q in 0..L::Q {
                                acc += fields.f[q * np + pb + x].as_f64();
                            }
                        }
                        Reduction::Momentum(a) => {
                            let mut m = 0.0f64;
                            for q in 0..L::Q {
                                m += L::C[q][a] as f64 * fields.f[q * np + pb + x].as_f64();
                            }
                            let fa = match ff {
                                Some(field) => {
                                    p.force[a].as_f64() + field[c0 + x][a].as_f64()
                                }
                                None => p.force[a].as_f64(),
                            };
                            acc += m + 0.5 * fa;
                        }
                    }
                }
            }
        }
        let _ = sub;
        acc
    }

    fn read_moments(&self, fields: &SoaFields<T>, out: &mut HostMoments<T>) {
        out.rho.clear();
        out.rho.extend_from_slice(&fields.rho);
        out.ux.clear();
        out.ux.extend_from_slice(&fields.ux);
        out.uy.clear();
        out.uy.extend_from_slice(&fields.uy);
        out.uz.clear();
        out.uz.extend_from_slice(&fields.uz);
    }
}
