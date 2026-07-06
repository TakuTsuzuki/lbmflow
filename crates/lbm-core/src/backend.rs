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

#[cfg(feature = "gpu")]
use crate::fields::LocalGeom;
use crate::fields::SoaFields;
use crate::halo::HaloExchange;
use crate::kernels::{
    collide_row, collide_row_central_moment, convective_face_selected, moments_row,
    outflow_face_selected, stream_row, zou_he_face_selected, FaceCellSelection, RawSlice, ZhKind,
};
use crate::lattice::{Face, Lattice};
use crate::params::{FaceBC, KParams, Reduction, SourceKind, StepParams};
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

/// Boundary shells complementing the interior box: fixed order YNeg row,
/// YPos row, XNeg column, XPos column (minus corners already covered), then
/// Z planes for 3D. Only the probe partials' summation order depends on
/// this; field results do not.
pub(crate) fn boundary_shells(sub: &Subdomain, interior: CellRange) -> Vec<CellRange> {
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

/// A compute + storage target for one subdomain's fields.
///
/// One time step, orchestrated by the solver, is:
/// `collide` → halo exchange → `stream` (interior, then boundary) → `swap`
/// → `apply_open_faces` → `update_moments` → `end_step`.
pub trait Backend<L: Lattice, T: Real> {
    /// Backend-owned field storage.
    ///
    /// This is intentionally a composite storage boundary owned by the
    /// backend, not an alias for the hydrodynamic `f` populations. Today the
    /// first member is the single distribution set (`f`, plus its ping-pong
    /// partner and moment/mask side fields). Future multiphase/scalar work can
    /// add additional distribution sets (`g`, `h`), per-cell properties, and
    /// Lagrangian buffers to this associated type while the solver continues
    /// to transfer through the host staging object at edit/read boundaries.
    type Fields;

    /// Allocate quiescent fields for a subdomain.
    fn alloc(&self, sub: &Subdomain) -> Self::Fields;

    /// Copy the host staging state into backend-owned storage.
    ///
    /// CPU backends use the same SoA layout on both sides, so this is a
    /// straight copy. Device backends strip/pack the host halo layout into
    /// their native composite storage.
    fn stage_in(&self, sub: &Subdomain, fields: &mut Self::Fields, host: &SoaFields<T>);

    /// Copy backend-owned storage back into host staging.
    ///
    /// This is the synchronization point used by setup shims and diagnostics
    /// that need population access. Moment-only diagnostics should prefer
    /// [`Backend::read_moments`], and scalar reductions should prefer
    /// [`Backend::reduce`], so device backends avoid unnecessary population
    /// readbacks.
    fn stage_out(&self, sub: &Subdomain, fields: &Self::Fields, host: &mut SoaFields<T>);

    /// Whether this backend's streaming kernel handles the single-part
    /// periodic halo itself. Transitional B-1 hook: stage 4 will make
    /// `HaloExchange` generic over `Backend::Fields`; until then this lets the
    /// monolithic GPU path use the common orchestrator without forcing a
    /// host-mediated halo copy between collide and stream.
    fn handles_single_part_periodic_halo(&self) -> bool {
        false
    }

    /// Whether this backend implements localized volume sources and masked
    /// face patches. CPU backends share the reference implementation; GPU
    /// rejects these features until device kernels are added.
    fn supports_localized_features(&self) -> bool {
        true
    }

    /// Whether this backend supports the orchestrator's interior/boundary
    /// streaming split. Backends that fuse whole-grid kernels can return
    /// `false` so callers can reject two-pass mode before a run is recorded.
    fn supports_two_pass(&self) -> bool {
        true
    }

    /// Whether this backend can compose per-mass gravity `rho(x) * g` inside
    /// the same backend-resident Guo force path used for uniform and per-cell
    /// force fields.
    fn supports_gravity_body_force(&self) -> bool {
        false
    }

    /// Exchange post-collision population halos for backend-owned fields.
    ///
    /// CPU backends delegate to the current `HaloExchange<SoaFields>`
    /// implementation. The monolithic GPU backend handles wrap/open-boundary
    /// inputs in-kernel. B-1 stage 4 will move the halo trait itself to
    /// `Backend::Fields`; this hook avoids host-mediated copies in the
    /// meantime.
    fn exchange_f<H: HaloExchange<T>>(
        &mut self,
        exchange: &H,
        subs: &[Subdomain],
        fields: &mut [Self::Fields],
    );

    /// TRT/BGK collision with Guo forcing over all core cells (in place).
    fn collide(&mut self, sub: &Subdomain, fields: &mut Self::Fields, p: &StepParams<T>);

    /// Pull-streaming into the out-buffer over `range`.
    ///
    /// Implementations that support force probes add their deterministic
    /// per-range contribution to backend-owned step state. Call
    /// [`Backend::read_probed_force`] after [`Backend::end_step`] to read the
    /// most recent step's force.
    fn stream(
        &mut self,
        sub: &Subdomain,
        fields: &mut Self::Fields,
        p: &StepParams<T>,
        range: CellRange,
    );

    /// Swap the population ping-pong pair (after all stream ranges ran).
    fn swap(&mut self, fields: &mut Self::Fields);

    /// Curved-wall post-stream correction. Default is a no-op for backends
    /// without Bouzidi support; CPU backends override it for `SoaFields`.
    fn apply_bouzidi(&mut self, _sub: &Subdomain, _fields: &mut Self::Fields, _p: &StepParams<T>) {}

    /// Open-face BC pass (Zou–He / outflow / convective) on the faces of
    /// this subdomain that lie on an open global face.
    fn apply_open_faces(&mut self, sub: &Subdomain, fields: &mut Self::Fields, p: &StepParams<T>);

    /// Recompute macroscopic moments from the populations.
    ///
    /// Lazy contract: backends may defer or partially elide the recompute when
    /// their step kernels have already produced moments equivalent to a full
    /// refresh. Any later [`Backend::read_moments`] must observe up-to-date
    /// moments for the completed step.
    fn update_moments(&mut self, sub: &Subdomain, fields: &mut Self::Fields, p: &StepParams<T>);

    /// Complete a solver step after all phases have been recorded.
    ///
    /// CPU backends execute eagerly and use the default no-op. Asynchronous
    /// backends use this hook for submit boundaries; readback APIs remain the
    /// only required blocking points.
    fn end_step(&mut self, _fields: &Self::Fields) {}

    /// Explicit readback of the momentum-exchange force accumulated over
    /// probed solid links during the most recent completed step.
    fn read_probed_force(&self, fields: &Self::Fields) -> [T; 3];

    /// Record/execute a span of whole steps. The default is exactly the
    /// generic orchestrator step loop; asynchronous backends may override it
    /// to keep their per-step hot loop inside backend-owned storage and
    /// recorder state.
    fn run_span<H: HaloExchange<T>>(
        &mut self,
        exchange: &H,
        subs: &[Subdomain],
        fields: &mut [Self::Fields],
        p: &StepParams<T>,
        two_pass: bool,
        probed_force: &mut [T; 3],
        steps: usize,
    ) {
        for _ in 0..steps {
            for i in 0..fields.len() {
                self.collide(&subs[i], &mut fields[i], p);
            }
            self.exchange_f(exchange, subs, fields);
            for i in 0..fields.len() {
                let sub = &subs[i];
                if !two_pass {
                    self.stream(sub, &mut fields[i], p, CellRange::full(sub));
                } else {
                    let c = sub.geom.core;
                    let interior = CellRange {
                        lo: [1, 1, if sub.geom.d == 3 { 1 } else { 0 }],
                        hi: [
                            c[0].saturating_sub(1),
                            c[1].saturating_sub(1),
                            if sub.geom.d == 3 {
                                c[2].saturating_sub(1)
                            } else {
                                c[2]
                            },
                        ],
                    };
                    self.stream(sub, &mut fields[i], p, interior);
                    for shell in boundary_shells(sub, interior) {
                        self.stream(sub, &mut fields[i], p, shell);
                    }
                };
            }
            for i in 0..fields.len() {
                self.apply_bouzidi(&subs[i], &mut fields[i], p);
            }
            for field in fields.iter_mut() {
                self.swap(field);
            }
            for i in 0..fields.len() {
                self.apply_open_faces(&subs[i], &mut fields[i], p);
            }
            for i in 0..fields.len() {
                self.apply_volume_sources(&subs[i], &mut fields[i], p);
            }
            for i in 0..fields.len() {
                self.update_moments(&subs[i], &mut fields[i], p);
            }
            for field in fields.iter() {
                self.end_step(field);
            }
            *probed_force = fields.iter().fold([T::zero(); 3], |a, field| {
                let b = self.read_probed_force(field);
                [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
            });
        }
    }

    /// Apply localized interior sources after boundary conditions and before
    /// moment recomputation.
    fn apply_volume_sources(
        &mut self,
        _sub: &Subdomain,
        _fields: &mut Self::Fields,
        _p: &StepParams<T>,
    ) {
    }

    /// Maximum number of steps the orchestrator should record before calling
    /// [`Backend::finish_run_chunk`]. CPU backends execute immediately, so the
    /// default runs the requested span as one chunk.
    fn run_chunk_size(&self, _fields: &[Self::Fields]) -> usize {
        usize::MAX
    }

    /// End-of-chunk hook used by asynchronous backends to submit recorded
    /// device work and, when needed, block until that chunk is complete.
    fn finish_run_chunk(&mut self, _fields: &[Self::Fields], _steps: usize) {}

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

    fn stage_in(&self, _sub: &Subdomain, fields: &mut SoaFields<T>, host: &SoaFields<T>) {
        *fields = host.clone();
    }

    fn stage_out(&self, _sub: &Subdomain, fields: &SoaFields<T>, host: &mut SoaFields<T>) {
        *host = fields.clone();
    }

    fn supports_gravity_body_force(&self) -> bool {
        true
    }

    fn exchange_f<H: HaloExchange<T>>(
        &mut self,
        exchange: &H,
        subs: &[Subdomain],
        fields: &mut [SoaFields<T>],
    ) {
        exchange.exchange_f::<L>(subs, fields);
    }

    fn collide(&mut self, sub: &Subdomain, fields: &mut SoaFields<T>, p: &StepParams<T>) {
        fields.probed_force = [T::zero(); 3];
        let kp = KParams::new::<L>(p);
        let g = fields.geom;
        let np = g.n_padded();
        let nx = g.core[0];
        let rows = g.core[1] * g.core[2];
        let f = RawSlice::new(&mut fields.f);
        let (rho, ux, uy, uz) = (&fields.rho, &fields.ux, &fields.uy, &fields.uz);
        let solid = &fields.solid;
        let ff = fields.force_field.as_deref();
        let omega = fields.omega_field.as_deref();
        let body = |r: usize| {
            let y = r % g.core[1];
            let z = r / g.core[1];
            let pb = g.pidx(0, y, z);
            let c0 = g.cidx(0, y, z);
            // SAFETY: each row index r is processed exactly once, and
            // collide_row writes only its own row's cells.
            unsafe {
                if kp.central_moment {
                    collide_row_central_moment::<L, T>(
                        f,
                        np,
                        pb,
                        &rho[c0..c0 + nx],
                        &ux[c0..c0 + nx],
                        &uy[c0..c0 + nx],
                        &uz[c0..c0 + nx],
                        &solid[pb..pb + nx],
                        ff.map(|v| &v[c0..c0 + nx]),
                        omega.map(|v| &v[c0..c0 + nx]),
                        &kp,
                    )
                } else {
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
                        omega.map(|v| &v[c0..c0 + nx]),
                        &kp,
                    )
                }
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
    ) {
        if range.is_empty() {
            return;
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
            let pf = fold(partials);
            fields.probed_force = [
                fields.probed_force[0] + pf[0],
                fields.probed_force[1] + pf[1],
                fields.probed_force[2] + pf[2],
            ];
            return;
        }
        let pf = fold((0..rows).map(body).collect());
        fields.probed_force = [
            fields.probed_force[0] + pf[0],
            fields.probed_force[1] + pf[1],
            fields.probed_force[2] + pf[2],
        ];
    }

    fn swap(&mut self, fields: &mut SoaFields<T>) {
        fields.swap_f();
    }

    fn apply_bouzidi(&mut self, _sub: &Subdomain, fields: &mut SoaFields<T>, p: &StepParams<T>) {
        let pf = crate::bouzidi::apply_bouzidi_impl::<L, T>(fields, p);
        fields.probed_force = [
            fields.probed_force[0] + pf[0],
            fields.probed_force[1] + pf[1],
            fields.probed_force[2] + pf[2],
        ];
    }

    fn apply_open_faces(&mut self, sub: &Subdomain, fields: &mut SoaFields<T>, p: &StepParams<T>) {
        apply_open_faces_impl::<L, T>(sub, fields, p);
    }

    fn apply_volume_sources(
        &mut self,
        sub: &Subdomain,
        fields: &mut SoaFields<T>,
        p: &StepParams<T>,
    ) {
        apply_volume_sources_impl::<L, T>(sub, fields, p);
    }

    fn update_moments(&mut self, sub: &Subdomain, fields: &mut SoaFields<T>, p: &StepParams<T>) {
        update_moments_impl::<L, T>(fields, p, self.use_parallel(sub));
    }

    fn reduce(
        &self,
        sub: &Subdomain,
        fields: &SoaFields<T>,
        p: &StepParams<T>,
        kind: Reduction,
    ) -> f64 {
        reduce_impl::<L, T>(sub, fields, p, kind)
    }

    fn read_moments(&self, fields: &SoaFields<T>, out: &mut HostMoments<T>) {
        read_moments_impl(fields, out);
    }

    fn read_probed_force(&self, fields: &SoaFields<T>) -> [T; 3] {
        fields.probed_force
    }
}

// ---------------------------------------------------------------------------
// Shared phase bodies (used verbatim by `CpuScalar` and `CpuSimd`, so the
// backends can never drift on the boundary-condition or diagnostic paths)
// ---------------------------------------------------------------------------

/// Open-face BC pass on the faces of this subdomain that lie on an open
/// global face (V1 `apply_open_edges` order: `Face::ALL`).
pub(crate) fn apply_open_faces_impl<L: Lattice, T: Real>(
    sub: &Subdomain,
    fields: &mut SoaFields<T>,
    p: &StepParams<T>,
) {
    let g = fields.geom;
    let np = g.n_padded();
    for face in Face::ALL {
        if !sub.touches_global_face(face) {
            continue;
        }
        let bc = &p.faces[face.index()];
        // Patch presence is a GLOBAL property of the face: a subdomain whose
        // local window contains no patch rect still owes the non-patch cells
        // their base treatment.
        let face_has_patches = p
            .face_patches
            .iter()
            .any(|patch| patch.face == face.index());
        if !bc.is_open() && !face_has_patches {
            continue;
        }
        // Patch rects are specified in GLOBAL in-face coordinates; the face
        // kernels iterate this subdomain's LOCAL core cells, so translate and
        // clip every rect by the subdomain origin (a patch straddling a seam
        // must land on the same global cells in every decomposition — T18.2).
        let patches: Vec<_> = p
            .face_patches
            .iter()
            .filter(|patch| patch.face == face.index())
            .filter_map(|patch| {
                patch_rect_local(sub, face, patch.lo, patch.hi).map(|(lo, hi)| (patch.bc, lo, hi))
            })
            .collect();
        let profiles = &fields.inlet_profiles;
        let excluded: Vec<_> = patches.iter().map(|(_, lo, hi)| (*lo, *hi)).collect();
        if !bc.is_open() && face_has_patches {
            // A Closed base face carrying patches is an impermeable lid on the
            // non-patch cells: prescribe u = 0 there (zero-velocity Zou-He).
            // Rim-covered cells are solid and skipped by the kernel, so faces
            // closed by a wall rim are unaffected; this arm exists for the
            // bare Closed-plus-patches face (the CR-2 motivating case). The
            // untreated alternative leaves those cells with no BC at all and
            // slowly diverges (T18.2 impinging jet, NaN at ~1.7k steps).
            zou_he_face_selected::<L, T>(
                &mut fields.f,
                np,
                &g,
                &fields.solid,
                face,
                &ZhKind::Velocity([T::zero(); 3]),
                None,
                FaceCellSelection::Excluding { rects: &excluded },
            );
        }
        if bc.is_open() {
            match bc {
                FaceBC::Closed => {}
                FaceBC::Velocity { u } => zou_he_face_selected::<L, T>(
                    &mut fields.f,
                    np,
                    &g,
                    &fields.solid,
                    face,
                    &ZhKind::Velocity(*u),
                    profiles[face.index()].as_deref(),
                    FaceCellSelection::Excluding { rects: &excluded },
                ),
                FaceBC::Pressure { rho } => zou_he_face_selected::<L, T>(
                    &mut fields.f,
                    np,
                    &g,
                    &fields.solid,
                    face,
                    &ZhKind::Pressure(*rho),
                    None,
                    FaceCellSelection::Excluding { rects: &excluded },
                ),
                FaceBC::Outflow => outflow_face_selected::<L, T>(
                    &mut fields.f,
                    np,
                    &g,
                    &fields.solid,
                    face,
                    FaceCellSelection::Excluding { rects: &excluded },
                ),
                FaceBC::Convective { u_conv } => convective_face_selected::<L, T>(
                    &mut fields.f,
                    np,
                    &g,
                    &fields.solid,
                    face,
                    *u_conv,
                    FaceCellSelection::Excluding { rects: &excluded },
                ),
            }
        }
        for (bc, lo, hi) in patches {
            match bc {
                // A Closed patch is an impermeable lid on its rect (same
                // reasoning as the Closed-base-with-patches arm above).
                FaceBC::Closed => zou_he_face_selected::<L, T>(
                    &mut fields.f,
                    np,
                    &g,
                    &fields.solid,
                    face,
                    &ZhKind::Velocity([T::zero(); 3]),
                    None,
                    FaceCellSelection::Rect { lo, hi },
                ),
                FaceBC::Velocity { u } => zou_he_face_selected::<L, T>(
                    &mut fields.f,
                    np,
                    &g,
                    &fields.solid,
                    face,
                    &ZhKind::Velocity(u),
                    None,
                    FaceCellSelection::Rect { lo, hi },
                ),
                FaceBC::Pressure { rho } => zou_he_face_selected::<L, T>(
                    &mut fields.f,
                    np,
                    &g,
                    &fields.solid,
                    face,
                    &ZhKind::Pressure(rho),
                    None,
                    FaceCellSelection::Rect { lo, hi },
                ),
                FaceBC::Outflow => outflow_face_selected::<L, T>(
                    &mut fields.f,
                    np,
                    &g,
                    &fields.solid,
                    face,
                    FaceCellSelection::Rect { lo, hi },
                ),
                FaceBC::Convective { u_conv } => convective_face_selected::<L, T>(
                    &mut fields.f,
                    np,
                    &g,
                    &fields.solid,
                    face,
                    u_conv,
                    FaceCellSelection::Rect { lo, hi },
                ),
            }
        }
    }
}

/// Translate a patch rect from global in-face coordinates to this subdomain's
/// local core coordinates, clipping to the owned tangent range. Returns `None`
/// when the rect does not intersect this subdomain's part of the face.
fn patch_rect_local(
    sub: &Subdomain,
    face: Face,
    lo: [usize; 2],
    hi: [usize; 2],
) -> Option<([usize; 2], [usize; 2])> {
    let (t1, t2) = face.tangents();
    let mut out_lo = [0usize; 2];
    let mut out_hi = [0usize; 2];
    for (i, t) in [t1, t2].into_iter().enumerate() {
        let o = sub.origin[t];
        let n = sub.geom.core[t];
        if hi[i] < o || lo[i] >= o + n {
            return None;
        }
        out_lo[i] = lo[i].saturating_sub(o);
        out_hi[i] = (hi[i] - o).min(n - 1);
    }
    Some((out_lo, out_hi))
}

/// Apply volume sources on owner core cells only. The source pass runs after
/// open-boundary BCs and before moment recomputation. Each source's `q_lu` is
/// divided uniformly over its inclusive region; the per-cell mass increment is
/// added as an equilibrium-shaped population delta, so `sum_q delta_f = q_cell`
/// exactly up to floating-point summation and Jet sources also carry first
/// moment `q_cell * u`.
pub(crate) fn apply_volume_sources_impl<L: Lattice, T: Real>(
    sub: &Subdomain,
    fields: &mut SoaFields<T>,
    p: &StepParams<T>,
) {
    if p.sources.is_empty() {
        return;
    }
    let g = fields.geom;
    let np = g.n_padded();
    for source in &p.sources {
        let lo = source.region.lo;
        let hi = source.region.hi;
        let count = (hi[0] - lo[0] + 1) * (hi[1] - lo[1] + 1) * (hi[2] - lo[2] + 1);
        let (q_lu, u) = match source.kind {
            SourceKind::MassFlow { q_lu } => (q_lu, [T::zero(); 3]),
            SourceKind::Jet { q_lu, u } => (q_lu, u),
        };
        let q_cell = q_lu / T::r(count as f64);
        let mut usq = u[0] * u[0];
        for a in 1..L::D {
            usq = usq + u[a] * u[a];
        }
        for gz in lo[2]..=hi[2] {
            if gz < sub.origin[2] || gz >= sub.origin[2] + g.core[2] {
                continue;
            }
            let z = gz - sub.origin[2];
            for gy in lo[1]..=hi[1] {
                if gy < sub.origin[1] || gy >= sub.origin[1] + g.core[1] {
                    continue;
                }
                let y = gy - sub.origin[1];
                for gx in lo[0]..=hi[0] {
                    if gx < sub.origin[0] || gx >= sub.origin[0] + g.core[0] {
                        continue;
                    }
                    let x = gx - sub.origin[0];
                    let i = g.pidx(x, y, z);
                    if fields.solid[i] {
                        continue;
                    }
                    for q in 0..L::Q {
                        let mut cu = T::r(L::C[q][0] as f64) * u[0];
                        for a in 1..L::D {
                            cu = cu + T::r(L::C[q][a] as f64) * u[a];
                        }
                        let delta = T::r(L::W[q])
                            * q_cell
                            * (T::one() + T::r(3.0) * cu + T::r(4.5) * cu * cu - T::r(1.5) * usq);
                        fields.f[q * np + i] = fields.f[q * np + i] + delta;
                    }
                }
            }
        }
    }
}

/// Full moment recompute over all core rows (V1 `update_moments`).
pub(crate) fn update_moments_impl<L: Lattice, T: Real>(
    fields: &mut SoaFields<T>,
    p: &StepParams<T>,
    parallel: bool,
) {
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
    if parallel {
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
    #[cfg(not(feature = "parallel"))]
    let _ = parallel;
    fields
        .rho
        .chunks_mut(nx)
        .zip(fields.ux.chunks_mut(nx))
        .zip(fields.uy.chunks_mut(nx).zip(fields.uz.chunks_mut(nx)))
        .enumerate()
        .for_each(body);
}

/// Backend-side reduction over fluid core cells in compact cell order
/// (z, y, x ascending, q inner) — V1's exact f64 accumulation sequence.
pub(crate) fn reduce_impl<L: Lattice, T: Real>(
    sub: &Subdomain,
    fields: &SoaFields<T>,
    p: &StepParams<T>,
    kind: Reduction,
) -> f64 {
    let g = fields.geom;
    let np = g.n_padded();
    let ff = fields.force_field.as_deref();
    let mut acc = 0.0f64;
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
                        let mut rho = 1.0f64;
                        for q in 0..L::Q {
                            rho += fields.f[q * np + pb + x].as_f64();
                        }
                        let field_force = ff.map(|field| &field[c0 + x]);
                        let gravity_force = p.gravity.map_or(0.0, |g| rho * g[a].as_f64());
                        let fa = match field_force {
                            Some(field) => {
                                p.force[a].as_f64() + (field[a].as_f64() + gravity_force)
                            }
                            None => p.force[a].as_f64() + gravity_force,
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

/// Explicit readback of the macroscopic fields into host memory.
pub(crate) fn read_moments_impl<T: Real>(fields: &SoaFields<T>, out: &mut HostMoments<T>) {
    out.rho.clear();
    out.rho.extend_from_slice(&fields.rho);
    out.ux.clear();
    out.ux.extend_from_slice(&fields.ux);
    out.uy.clear();
    out.uy.extend_from_slice(&fields.uy);
    out.uz.clear();
    out.uz.extend_from_slice(&fields.uz);
}

/// Copy compact-core host moments into the corresponding part of a padded
/// [`SoaFields`] host-staging object.
#[cfg(feature = "gpu")]
pub(crate) fn write_host_moments<T: Real>(
    geom: LocalGeom,
    moments: &HostMoments<T>,
    host: &mut SoaFields<T>,
) {
    let mut c = 0usize;
    for z in 0..geom.core[2] {
        for y in 0..geom.core[1] {
            for x in 0..geom.core[0] {
                let ci = geom.cidx(x, y, z);
                host.rho[ci] = moments.rho[c];
                host.ux[ci] = moments.ux[c];
                host.uy[ci] = moments.uy[c];
                host.uz[ci] = moments.uz[c];
                c += 1;
            }
        }
    }
}
