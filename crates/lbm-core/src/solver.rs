//! Solver orchestrator: drives subdomains × backend × halo exchange through
//! the V1 step sequence (collide → stream → open faces → moments).
//!
//! The step-phase order, the diagnostics' f64 accumulation and the
//! initialisation paths reproduce V1 `Simulation` mechanics exactly; the
//! compat facade is a thin wrapper over this type with a monolithic (1×1×1)
//! decomposition.

use crate::backend::{Backend, HostMoments};
use crate::fields::SoaFields;
use crate::halo::{exchange_g_generic, ExchangeScope, HaloExchange};
use crate::kernels::equilibrium;
use crate::lattice::{Face, Lattice, D3Q19};
use crate::params::{
    CollisionKind, FaceBC, FacePatch, KParams, Reduction, SourceKind, SourceRegion, StepParams,
    VolumeSource, MAX_SPEED,
};
use crate::real::Real;
use crate::rotating_ibm::{
    add3, cross, marker_stencil, norm, to_real3, DirectForcingConfig, IbmDiagnostics, RotatingBody,
};
use crate::subdomain::Subdomain;
use crate::wall_model::{friction_velocity, WallCellMetric, WallMetricSource};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::io::{Read, Write};
use std::path::Path;

const CKPT_FORMAT_VERSION: u32 = 2;
const CKPT_MAGIC: &[u8; 8] = b"LBMKPT\0\0";
const SEC_F_PRIMARY: u32 = 1;
const SEC_STALE_STASH: u32 = 2;
const SEC_MOMENTS: u32 = 3;
const SEC_SOLID: u32 = 4;
const SEC_FORCE_FIELD: u32 = 5;

/// Structured checkpoint/restart failure. The `code` field is intentionally
/// stable for CLI/MCP agents.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CheckpointError {
    pub code: &'static str,
    pub message: String,
}

impl CheckpointError {
    fn new(code: &'static str, message: impl Into<String>) -> Self {
        Self {
            code,
            message: message.into(),
        }
    }
}

impl std::fmt::Display for CheckpointError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}: {}", self.code, self.message)
    }
}

impl std::error::Error for CheckpointError {}

impl From<std::io::Error> for CheckpointError {
    fn from(e: std::io::Error) -> Self {
        Self::new("CKPT_IO", e.to_string())
    }
}

impl From<serde_json::Error> for CheckpointError {
    fn from(e: serde_json::Error) -> Self {
        Self::new("CKPT_MANIFEST_INVALID", e.to_string())
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub(crate) struct CheckpointManifest {
    kind: String,
    format_version: u32,
    step: u64,
    time: f64,
    dtype: String,
    lattice: String,
    global: [usize; 3],
    scenario_hash: String,
    decomp_hash: String,
    nranks: usize,
    ranks: Vec<CheckpointRank>,
    reserved: BTreeMap<String, bool>,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub(crate) struct CheckpointRank {
    pub(crate) rank: usize,
    pub(crate) part: usize,
    pub(crate) file: String,
    pub(crate) origin: [usize; 3],
    pub(crate) core: [usize; 3],
    pub(crate) bytes: u64,
    pub(crate) payload_hash: String,
    pub(crate) mask_hash: String,
}

#[derive(Clone, Copy, Debug)]
struct SectionEntry {
    id: u32,
    offset: u64,
    byte_len: u64,
}

#[derive(Clone, Debug)]
struct RankHeader {
    dtype: u8,
    lattice_id: u16,
    q: u16,
    d: u16,
    np: u64,
    n_core: u64,
    payload_hash: u64,
}

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
    /// Localized interior volume sources/sinks.
    pub sources: Vec<VolumeSource<T>>,
    /// Rectangular per-face boundary-condition overrides.
    pub face_patches: Vec<FacePatch<T>>,
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
            sources: Vec::new(),
            face_patches: Vec::new(),
        }
    }
}

/// A rejected [`GlobalSpec`] (A-4): the V2-native counterpart of the compat
/// facade's `ConfigError`. `Solver::build` calls [`GlobalSpec::validate`]
/// before allocating, turning the previously-silent non-physical
/// configurations (stale data on an uncovered face, ν = 0, periodic × open
/// on one axis, …) into hard errors.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum UnsupportedReason {
    NotImplemented,
    OutOfValidityRange { detail: String },
    MissingDependency { depends_on: String },
    EvidenceGateFailed { missing: Vec<String> },
    DemoOnly { rationale: String },
}

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
    /// Central-moment shear relaxation must be finite and in range.
    InvalidCentralMomentRate {
        /// Offending value.
        omega_shear: f64,
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
    /// Open-face kernels support D2Q9/D3Q19 open types, and D3Q27 velocity /
    /// pressure faces only.
    UnsupportedOpenFaceLattice {
        /// Lattice name.
        lattice: &'static str,
        /// Number of unknown populations on each straight face.
        unknowns: usize,
    },
    /// This lattice has an open-face kernel, but not for the requested face
    /// kind.
    UnsupportedOpenFaceKind {
        /// Lattice name.
        lattice: &'static str,
        /// The offending face index.
        face: usize,
        /// Boundary kind name.
        kind: &'static str,
    },
    /// A volume-source region is outside the domain or touches a global face.
    SourceRegionNotInterior {
        /// Offending source index.
        source: usize,
        /// Region low corner.
        lo: [usize; 3],
        /// Region high corner.
        hi: [usize; 3],
    },
    /// Two volume-source regions overlap.
    SourceOverlap {
        /// First source index.
        a: usize,
        /// Second source index.
        b: usize,
    },
    /// A volume-source region covers a solid cell.
    SourceOverlapsSolid {
        /// Offending source index.
        source: usize,
        /// First solid cell found.
        cell: [usize; 3],
    },
    /// A sink removes too much mass from each cell in one step.
    SourceSinkTooStrong {
        /// Offending source index.
        source: usize,
        /// Per-cell mass increment.
        q_cell: f64,
    },
    /// A face patch references an inactive/out-of-range face or lies outside
    /// the face's tangent-coordinate bounds.
    FacePatchOutOfBounds {
        /// Offending patch index.
        patch: usize,
        /// Face index.
        face: usize,
        /// Patch low coordinate.
        lo: [usize; 2],
        /// Patch high coordinate.
        hi: [usize; 2],
    },
    /// Two patches on the same face overlap.
    FacePatchOverlap {
        /// First patch index.
        a: usize,
        /// Second patch index.
        b: usize,
    },
    /// GPU backends do not yet implement localized sources or face patches.
    UnsupportedOnGpu {
        /// Feature name.
        feature: &'static str,
        /// Machine-readable reason for the rejection.
        reason: UnsupportedReason,
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
            SpecError::InvalidCentralMomentRate { omega_shear } => write!(
                f,
                "central_moment omega_shear must be finite and in (0, 2] (got {omega_shear})"
            ),
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
            SpecError::UnsupportedOpenFaceLattice { lattice, unknowns } => write!(
                f,
                "{lattice} has {unknowns} unknown populations per open face; no open-face \
                 closure is implemented for this lattice"
            ),
            SpecError::UnsupportedOpenFaceKind {
                lattice,
                face,
                kind,
            } => write!(
                f,
                "{lattice} open face {face} uses {kind}, but this lattice currently supports \
                 only velocity inlet and pressure outlet open faces"
            ),
            SpecError::SourceRegionNotInterior { source, lo, hi } => write!(
                f,
                "source {source} region {lo:?}..={hi:?} must be inside the domain and at least one cell from every face"
            ),
            SpecError::SourceOverlap { a, b } => {
                write!(f, "source regions {a} and {b} overlap")
            }
            SpecError::SourceOverlapsSolid { source, cell } => {
                write!(f, "source {source} overlaps solid cell {cell:?}")
            }
            SpecError::SourceSinkTooStrong { source, q_cell } => write!(
                f,
                "source {source} sink removes {q_cell} mass per cell per step; q_cell must be > -1.0 to keep reference-density cells positive"
            ),
            SpecError::FacePatchOutOfBounds { patch, face, lo, hi } => write!(
                f,
                "face patch {patch} on face {face} has out-of-bounds rectangle {lo:?}..={hi:?}"
            ),
            SpecError::FacePatchOverlap { a, b } => {
                write!(f, "face patches {a} and {b} overlap on the same face")
            }
            SpecError::UnsupportedOnGpu { feature, .. } => {
                write!(f, "GPU backend does not yet support {feature}")
            }
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
        // Collision-specific relaxation parameters.
        match self.collision {
            CollisionKind::Trt { magic } => {
                if !magic.is_finite() || !(magic > 0.0) {
                    return Err(SpecError::InvalidMagic { magic });
                }
            }
            CollisionKind::CentralMoment { omega_shear } => {
                if !omega_shear.is_finite() || !(omega_shear > 0.0 && omega_shear <= 2.0) {
                    return Err(SpecError::InvalidCentralMomentRate { omega_shear });
                }
            }
            CollisionKind::Bgk => {}
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
        // Per-axis: periodic × open exclusivity, and gather open axes from
        // both whole-face BCs and patch BCs.
        let mut open_axes = 0usize;
        for a in 0..d {
            let (neg, pos) = (Face::ALL[2 * a], Face::ALL[2 * a + 1]);
            let axis_open = self.faces[neg.index()].is_open()
                || self.faces[pos.index()].is_open()
                || self
                    .face_patches
                    .iter()
                    .any(|p| p.face < 6 && Face::ALL[p.face].axis() == a && p.bc.is_open());
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
                    FaceBC::Outflow => {}
                    _ => validate_face_bc(*bc)?,
                }
            } else {
                // Closed, non-periodic face: it must be a full solid wall rim,
                // unless open patches cover part of the face. A Closed base
                // with open patches is legal: closed cells are handled by the
                // solid rim and patched cells run the open-BC pass.
                if !face_is_full_solid_rim(face, self.dims, solid)
                    && !self.face_patches.iter().any(|p| p.face == face.index())
                {
                    return Err(SpecError::UncoveredFace { face: face.index() });
                }
            }
        }
        self.validate_sources(d, solid)?;
        self.validate_face_patches(d)?;
        Ok(())
    }

    fn validate_sources(&self, d: usize, solid: &[bool]) -> Result<(), SpecError> {
        for (i, source) in self.sources.iter().enumerate() {
            let SourceRegion { lo, hi } = source.region;
            let mut volume = 1usize;
            for a in 0..3 {
                let active = a < d;
                let expected = if active { self.dims[a] } else { 1 };
                if lo[a] > hi[a]
                    || hi[a] >= expected
                    || (active && (lo[a] == 0 || hi[a] + 1 >= expected))
                {
                    return Err(SpecError::SourceRegionNotInterior { source: i, lo, hi });
                }
                volume = volume.saturating_mul(hi[a] - lo[a] + 1);
            }
            let q_lu = match source.kind {
                SourceKind::MassFlow { q_lu } => q_lu,
                SourceKind::Jet { q_lu, u } => {
                    validate_velocity(u)?;
                    q_lu
                }
            };
            let q = q_lu.as_f64();
            if !q.is_finite() {
                return Err(SpecError::NonFiniteParameter {
                    what: "source q_lu",
                });
            }
            let q_cell = q / volume as f64;
            // Conservative positivity guard for deviation-form density:
            // a reference-density cell (rho = 1) remains positive after the
            // source pass. Stronger sinks need smaller time steps or a model
            // that couples sink strength to the current local rho.
            if !(q_cell > -1.0) {
                return Err(SpecError::SourceSinkTooStrong { source: i, q_cell });
            }
            if !solid.is_empty() {
                for z in lo[2]..=hi[2] {
                    for y in lo[1]..=hi[1] {
                        for x in lo[0]..=hi[0] {
                            let gi = (z * self.dims[1] + y) * self.dims[0] + x;
                            if solid[gi] {
                                return Err(SpecError::SourceOverlapsSolid {
                                    source: i,
                                    cell: [x, y, z],
                                });
                            }
                        }
                    }
                }
            }
        }
        for a in 0..self.sources.len() {
            for b in a + 1..self.sources.len() {
                if regions_overlap(self.sources[a].region, self.sources[b].region) {
                    return Err(SpecError::SourceOverlap { a, b });
                }
            }
        }
        Ok(())
    }

    /// Validate against the full lattice, including open-face kernel support.
    pub fn validate_lattice<L: Lattice>(&self, solid: &[bool]) -> Result<(), SpecError> {
        self.validate(L::D, solid)?;
        if (0..L::D).any(|a| {
            let neg = Face::ALL[2 * a];
            let pos = Face::ALL[2 * a + 1];
            self.faces[neg.index()].is_open() || self.faces[pos.index()].is_open()
        }) {
            let unknowns = L::unknowns(Face::ALL[0]).len();
            if unknowns != 3 && unknowns != 5 && unknowns != 9 {
                return Err(SpecError::UnsupportedOpenFaceLattice {
                    lattice: lattice_name::<L>(),
                    unknowns,
                });
            }
        }
        Ok(())
    }

    fn validate_face_patches(&self, d: usize) -> Result<(), SpecError> {
        for (i, patch) in self.face_patches.iter().enumerate() {
            if patch.face >= 6 {
                return Err(SpecError::FacePatchOutOfBounds {
                    patch: i,
                    face: patch.face,
                    lo: patch.lo,
                    hi: patch.hi,
                });
            }
            let face = Face::ALL[patch.face];
            if face.axis() >= d {
                return Err(SpecError::FacePatchOutOfBounds {
                    patch: i,
                    face: patch.face,
                    lo: patch.lo,
                    hi: patch.hi,
                });
            }
            let (t1, t2) = face.tangents();
            if patch.lo[0] > patch.hi[0]
                || patch.lo[1] > patch.hi[1]
                || patch.hi[0] >= self.dims[t1]
                || patch.hi[1] >= self.dims[t2]
            {
                return Err(SpecError::FacePatchOutOfBounds {
                    patch: i,
                    face: patch.face,
                    lo: patch.lo,
                    hi: patch.hi,
                });
            }
            validate_face_bc(patch.bc)?;
            if patch.bc.is_open() && self.dims[face.axis()] < 3 {
                return Err(SpecError::OpenFaceAxisTooShort {
                    face: patch.face,
                    extent: self.dims[face.axis()],
                });
            }
        }
        for a in 0..self.face_patches.len() {
            for b in a + 1..self.face_patches.len() {
                let pa = self.face_patches[a];
                let pb = self.face_patches[b];
                if pa.face == pb.face && rects_overlap(pa.lo, pa.hi, pb.lo, pb.hi) {
                    return Err(SpecError::FacePatchOverlap { a, b });
                }
            }
        }
        Ok(())
    }
}

fn validate_face_bc<T: Real>(bc: FaceBC<T>) -> Result<(), SpecError> {
    match bc {
        FaceBC::Velocity { u } => validate_velocity(u),
        FaceBC::Pressure { rho } => {
            let r = rho.as_f64();
            if !(r > 0.0) {
                return Err(SpecError::NonPositiveDensity { rho: r });
            }
            Ok(())
        }
        FaceBC::Convective { u_conv } => {
            let v = u_conv.as_f64();
            if !(v > 0.0 && v <= 1.0) {
                return Err(SpecError::InvalidConvectiveSpeed { u_conv: v });
            }
            Ok(())
        }
        FaceBC::Closed | FaceBC::Outflow => Ok(()),
    }
}

fn validate_velocity<T: Real>(u: [T; 3]) -> Result<(), SpecError> {
    let mut sq = 0.0f64;
    for c in u {
        let v = c.as_f64();
        if !v.is_finite() {
            return Err(SpecError::NonFiniteParameter { what: "velocity" });
        }
        sq += v * v;
    }
    let speed = sq.sqrt();
    if !(speed <= MAX_SPEED) {
        return Err(SpecError::VelocityTooHigh { speed });
    }
    Ok(())
}

fn regions_overlap(a: SourceRegion, b: SourceRegion) -> bool {
    (0..3).all(|i| a.lo[i] <= b.hi[i] && b.lo[i] <= a.hi[i])
}

fn rects_overlap(alo: [usize; 2], ahi: [usize; 2], blo: [usize; 2], bhi: [usize; 2]) -> bool {
    alo[0] <= bhi[0] && blo[0] <= ahi[0] && alo[1] <= bhi[1] && blo[1] <= ahi[1]
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

fn fnv1a64(bytes: &[u8]) -> u64 {
    let mut h = 0xcbf29ce484222325u64;
    hash_bytes(&mut h, bytes);
    h
}

fn hash_bytes(h: &mut u64, bytes: &[u8]) {
    for &b in bytes {
        *h ^= b as u64;
        *h = h.wrapping_mul(0x100000001b3);
    }
}

fn hash_u64(h: &mut u64, v: u64) {
    hash_bytes(h, &v.to_le_bytes());
}

fn hash_f64(h: &mut u64, v: f64) {
    hash_bytes(h, &v.to_bits().to_le_bytes());
}

fn hash_real<T: Real>(h: &mut u64, v: T) {
    if std::mem::size_of::<T>() == 4 {
        hash_bytes(h, &(v.as_f64() as f32).to_bits().to_le_bytes());
    } else {
        hash_bytes(h, &v.as_f64().to_bits().to_le_bytes());
    }
}

fn hash_face_bc<T: Real>(h: &mut u64, bc: FaceBC<T>) {
    match bc {
        FaceBC::Closed => hash_u64(h, 0),
        FaceBC::Velocity { u } => {
            hash_u64(h, 1);
            for v in u {
                hash_real(h, v);
            }
        }
        FaceBC::Pressure { rho } => {
            hash_u64(h, 2);
            hash_real(h, rho);
        }
        FaceBC::Outflow => hash_u64(h, 3),
        FaceBC::Convective { u_conv } => {
            hash_u64(h, 4);
            hash_real(h, u_conv);
        }
    }
}

fn spec_hash<T: Real, L: Lattice>(
    spec: &GlobalSpec<T>,
    _solid: &[bool],
    _wall_u: &[[T; 3]],
) -> u64 {
    let mut h = 0xcbf29ce484222325u64;
    hash_u64(&mut h, L::D as u64);
    hash_u64(&mut h, L::Q as u64);
    for v in spec.dims {
        hash_u64(&mut h, v as u64);
    }
    hash_f64(&mut h, spec.nu);
    match spec.collision {
        CollisionKind::Bgk => hash_u64(&mut h, 1),
        CollisionKind::Trt { magic } => {
            hash_u64(&mut h, 2);
            hash_f64(&mut h, magic);
        }
        CollisionKind::CentralMoment { omega_shear } => {
            hash_u64(&mut h, 3);
            hash_f64(&mut h, omega_shear);
        }
    }
    for v in spec.periodic {
        hash_u64(&mut h, u64::from(v));
    }
    for bc in spec.faces {
        hash_face_bc(&mut h, bc);
    }
    for v in spec.force {
        hash_real(&mut h, v);
    }
    for source in &spec.sources {
        for v in source.region.lo {
            hash_u64(&mut h, v as u64);
        }
        for v in source.region.hi {
            hash_u64(&mut h, v as u64);
        }
        match source.kind {
            SourceKind::MassFlow { q_lu } => {
                hash_u64(&mut h, 1);
                hash_real(&mut h, q_lu);
            }
            SourceKind::Jet { q_lu, u } => {
                hash_u64(&mut h, 2);
                hash_real(&mut h, q_lu);
                for v in u {
                    hash_real(&mut h, v);
                }
            }
        }
    }
    for patch in &spec.face_patches {
        hash_u64(&mut h, patch.face as u64);
        for v in patch.lo {
            hash_u64(&mut h, v as u64);
        }
        for v in patch.hi {
            hash_u64(&mut h, v as u64);
        }
        hash_face_bc(&mut h, patch.bc);
    }
    h
}

pub(crate) fn decomp_hash(subs: &[Subdomain], periodic: [bool; 3]) -> u64 {
    let mut h = 0xcbf29ce484222325u64;
    hash_u64(&mut h, subs.len() as u64);
    for v in periodic {
        hash_u64(&mut h, u64::from(v));
    }
    for sub in subs {
        for v in sub.origin {
            hash_u64(&mut h, v as u64);
        }
        for v in sub.geom.core {
            hash_u64(&mut h, v as u64);
        }
        for nb in sub.neighbors {
            hash_u64(&mut h, nb.map(|v| v as u64 + 1).unwrap_or(0));
        }
    }
    h
}

fn part_mask_hash<T: Real>(fields: &SoaFields<T>) -> u64 {
    let mut h = 0xcbf29ce484222325u64;
    hash_u64(&mut h, fields.solid.len() as u64);
    for &v in &fields.solid {
        hash_u64(&mut h, u64::from(v));
    }
    for u in &fields.wall_u {
        for v in *u {
            hash_real(&mut h, v);
        }
    }
    h
}

fn hash_string(h: u64) -> String {
    format!("fnv1a64:{h:016x}")
}

fn dtype_name<T: Real>() -> &'static str {
    if std::mem::size_of::<T>() == 4 {
        "f32"
    } else {
        "f64"
    }
}

fn lattice_name<L: Lattice>() -> &'static str {
    match (L::D, L::Q) {
        (2, 9) => "D2Q9",
        (3, 19) => "D3Q19",
        (3, 27) => "D3Q27",
        _ => "unknown",
    }
}

fn lattice_id<L: Lattice>() -> u16 {
    match (L::D, L::Q) {
        (2, 9) => 29,
        (3, 19) => 319,
        (3, 27) => 327,
        _ => 0,
    }
}

fn push_real_bytes<T: Real>(out: &mut Vec<u8>, values: &[T]) {
    if std::mem::size_of::<T>() == 4 {
        for &v in values {
            out.extend_from_slice(&(v.as_f64() as f32).to_bits().to_le_bytes());
        }
    } else {
        for &v in values {
            out.extend_from_slice(&v.as_f64().to_bits().to_le_bytes());
        }
    }
}

fn read_real_bytes<T: Real>(bytes: &[u8], expected: usize) -> Result<Vec<T>, CheckpointError> {
    let width = std::mem::size_of::<T>();
    if bytes.len() != expected * width {
        return Err(CheckpointError::new(
            "CKPT_GEOM_MISMATCH",
            format!(
                "section length {} does not match expected {} values of {width} bytes",
                bytes.len(),
                expected
            ),
        ));
    }
    let mut out = Vec::with_capacity(expected);
    if width == 4 {
        for chunk in bytes.chunks_exact(4) {
            let v = f32::from_bits(u32::from_le_bytes(chunk.try_into().unwrap()));
            out.push(T::r(v as f64));
        }
    } else {
        for chunk in bytes.chunks_exact(8) {
            let v = f64::from_bits(u64::from_le_bytes(chunk.try_into().unwrap()));
            out.push(T::r(v));
        }
    }
    Ok(out)
}

fn required_section<'a>(
    sections: &'a BTreeMap<u32, Vec<u8>>,
    id: u32,
    name: &'static str,
) -> Result<&'a [u8], CheckpointError> {
    sections.get(&id).map(Vec::as_slice).ok_or_else(|| {
        CheckpointError::new(
            "CKPT_TRUNCATED",
            format!("checkpoint is missing required section {name}"),
        )
    })
}

fn write_u16(out: &mut Vec<u8>, v: u16) {
    out.extend_from_slice(&v.to_le_bytes());
}

fn write_u32(out: &mut Vec<u8>, v: u32) {
    out.extend_from_slice(&v.to_le_bytes());
}

fn write_u64(out: &mut Vec<u8>, v: u64) {
    out.extend_from_slice(&v.to_le_bytes());
}

fn take<'a>(bytes: &'a [u8], pos: &mut usize, n: usize) -> Result<&'a [u8], CheckpointError> {
    if bytes.len().saturating_sub(*pos) < n {
        return Err(CheckpointError::new(
            "CKPT_TRUNCATED",
            "rank file ended before the declared header/table/payload",
        ));
    }
    let out = &bytes[*pos..*pos + n];
    *pos += n;
    Ok(out)
}

fn read_u8(bytes: &[u8], pos: &mut usize) -> Result<u8, CheckpointError> {
    Ok(take(bytes, pos, 1)?[0])
}

fn read_u16(bytes: &[u8], pos: &mut usize) -> Result<u16, CheckpointError> {
    Ok(u16::from_le_bytes(take(bytes, pos, 2)?.try_into().unwrap()))
}

fn read_u32(bytes: &[u8], pos: &mut usize) -> Result<u32, CheckpointError> {
    Ok(u32::from_le_bytes(take(bytes, pos, 4)?.try_into().unwrap()))
}

fn read_u64(bytes: &[u8], pos: &mut usize) -> Result<u64, CheckpointError> {
    Ok(u64::from_le_bytes(take(bytes, pos, 8)?.try_into().unwrap()))
}

fn read_rank_file(path: &Path) -> Result<(RankHeader, BTreeMap<u32, Vec<u8>>), CheckpointError> {
    let mut bytes = Vec::new();
    std::fs::File::open(path)?.read_to_end(&mut bytes)?;
    let mut pos = 0usize;
    if take(&bytes, &mut pos, 8)? != CKPT_MAGIC {
        return Err(CheckpointError::new(
            "CKPT_BAD_MAGIC",
            "rank file magic is not LBMKPT",
        ));
    }
    let format_ver = read_u32(&bytes, &mut pos)?;
    if format_ver != CKPT_FORMAT_VERSION {
        return Err(CheckpointError::new(
            "CKPT_VERSION_MISMATCH",
            format!("rank file format_version {format_ver} differs from supported {CKPT_FORMAT_VERSION}"),
        ));
    }
    let endian = read_u8(&bytes, &mut pos)?;
    if endian != 0 {
        return Err(CheckpointError::new(
            "CKPT_BAD_MAGIC",
            "cross-endian checkpoints are not supported by this checkpoint format",
        ));
    }
    let dtype = read_u8(&bytes, &mut pos)?;
    let lattice_id = read_u16(&bytes, &mut pos)?;
    let q = read_u16(&bytes, &mut pos)?;
    let d = read_u16(&bytes, &mut pos)?;
    let np = read_u64(&bytes, &mut pos)?;
    let n_core = read_u64(&bytes, &mut pos)?;
    let section_count = read_u32(&bytes, &mut pos)? as usize;
    let mut table = Vec::with_capacity(section_count);
    for _ in 0..section_count {
        table.push(SectionEntry {
            id: read_u32(&bytes, &mut pos)?,
            offset: read_u64(&bytes, &mut pos)?,
            byte_len: read_u64(&bytes, &mut pos)?,
        });
    }
    if bytes.len().saturating_sub(pos) < 8 {
        return Err(CheckpointError::new(
            "CKPT_TRUNCATED",
            "rank file is missing payload hash trailer",
        ));
    }
    let payload_end = bytes.len() - 8;
    let payload_hash = u64::from_le_bytes(bytes[payload_end..].try_into().unwrap());
    let payload = &bytes[pos..payload_end];
    if fnv1a64(payload) != payload_hash {
        return Err(CheckpointError::new(
            "CKPT_PAYLOAD_CORRUPT",
            "rank file payload hash does not match its contents",
        ));
    }
    let mut sections = BTreeMap::new();
    for entry in table {
        let start = entry.offset as usize;
        let len = entry.byte_len as usize;
        let end = start.checked_add(len).ok_or_else(|| {
            CheckpointError::new("CKPT_TRUNCATED", "section offset overflows usize")
        })?;
        if end > payload.len() {
            return Err(CheckpointError::new(
                "CKPT_TRUNCATED",
                format!(
                    "section {} range {}..{} exceeds payload length {}",
                    entry.id,
                    start,
                    end,
                    payload.len()
                ),
            ));
        }
        sections.insert(entry.id, payload[start..end].to_vec());
    }
    Ok((
        RankHeader {
            dtype,
            lattice_id,
            q,
            d,
            np,
            n_core,
            payload_hash,
        },
        sections,
    ))
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

fn tangential_speed<T: Real>(u: [T; 3], wall_u: [T; 3], normal: [T; 3]) -> T {
    let rel = [u[0] - wall_u[0], u[1] - wall_u[1], u[2] - wall_u[2]];
    let un = rel[0] * normal[0] + rel[1] * normal[1] + rel[2] * normal[2];
    let tan = [
        rel[0] - un * normal[0],
        rel[1] - un * normal[1],
        rel[2] - un * normal[2],
    ];
    (tan[0] * tan[0] + tan[1] * tan[1] + tan[2] * tan[2]).sqrt()
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
    collision: CollisionKind,
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
        Self::try_new(spec, solid, wall_u, decomp, backend, exchange)
            .unwrap_or_else(|e| panic!("invalid GlobalSpec: {e}"))
    }

    /// Fallible variant of [`Self::new`] that returns validation failures
    /// before allocating or entering backend kernels.
    pub fn try_new(
        spec: &GlobalSpec<T>,
        solid: &[bool],
        wall_u: &[[T; 3]],
        decomp: [usize; 3],
        backend: B,
        exchange: H,
    ) -> Result<Self, SpecError> {
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
        Self::try_new_local_part(spec, solid, wall_u, decomp, part, backend, exchange)
            .unwrap_or_else(|e| panic!("invalid GlobalSpec: {e}"))
    }

    /// Fallible variant of [`Self::new_local_part`].
    pub fn try_new_local_part(
        spec: &GlobalSpec<T>,
        solid: &[bool],
        wall_u: &[[T; 3]],
        decomp: [usize; 3],
        part: usize,
        backend: B,
        exchange: H,
    ) -> Result<Self, SpecError> {
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
    ) -> Result<Self, SpecError> {
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
        spec.validate_lattice::<L>(solid)?;
        if (!spec.sources.is_empty() || !spec.face_patches.is_empty())
            && !backend.supports_localized_features()
        {
            let feature = if !spec.sources.is_empty() {
                "localized volume sources"
            } else {
                "masked face patches"
            };
            return Err(SpecError::UnsupportedOnGpu {
                feature,
                reason: UnsupportedReason::NotImplemented,
            });
        }
        let has_open_faces = spec.faces.iter().any(FaceBC::is_open)
            || spec.face_patches.iter().any(|patch| patch.bc.is_open());
        if L::D == 3 && L::Q == 27 && has_open_faces && !backend.supports_d3q27_open_faces() {
            return Err(SpecError::UnsupportedOnGpu {
                feature: "D3Q27 open faces",
                reason: UnsupportedReason::NotImplemented,
            });
        }
        let (omega_p, omega_m) = spec.collision.omegas(spec.nu);
        let params = StepParams {
            collision: spec.collision,
            omega_p,
            omega_m,
            force: spec.force,
            gravity: None,
            faces: spec.faces,
            sources: spec.sources.clone(),
            face_patches: spec.face_patches.clone(),
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
            collision: spec.collision,
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
        Ok(solver)
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

    /// Compose gravity into the existing Guo force field for one step.
    ///
    /// This is the gravity source-composition point required by FR-BC-02:
    /// the caller-owned per-cell force field is `F_s + F_user + ...`, and
    /// this method adds the single-phase gravity term `rho(x) * g` on fluid
    /// cells before collision. W-VOF must replace the density factor at this
    /// exact line with `rho(phi)` (or the consistent AGG density field once
    /// the phase-flux correction lands), leaving the Guo forcing scheme and
    /// all other force sources unchanged. In dynamic-pressure notation the
    /// future well-balanced residual is composed here as
    /// `F_s + (rho(phi) - rho_h) * g + F_b^scalar + ...`; single-phase
    /// compatibility currently uses `rho_h = 0`, which preserves the landed
    /// public contract that `set_gravity(g)` is bit-identical to a raw
    /// per-cell force field filled with `rho(x) * g`.
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
        self.refresh_probed_force();
        self.stage_out_all();
        self.unstage_gravity(gravity_stage);
        self.stage_in_if_dirty();
    }

    fn refresh_probed_force(&mut self) {
        self.probed_force = self.parts.iter().fold([T::zero(); 3], |a, field| {
            let b = self.backend.read_probed_force(field);
            [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
        });
    }

    fn params_with_backend_gravity(&self) -> StepParams<T> {
        let mut params = self.params.clone();
        params.gravity = self.gravity;
        params
    }

    fn refresh_moments_after_force_change(&mut self) {
        self.stage_in_if_dirty();
        let params_with_gravity;
        let backend_gravity = self.gravity.is_some() && self.backend.supports_gravity_body_force();
        let params = if backend_gravity {
            params_with_gravity = self.params_with_backend_gravity();
            &params_with_gravity
        } else {
            &self.params
        };
        for i in 0..self.parts.len() {
            self.backend
                .update_moments(&self.subs[i], &mut self.parts[i], params);
        }
        self.device_ahead = true;
    }

    fn force_field_is_uniform(&self) -> bool {
        let mut first = None;
        for fields in &self.host_parts {
            let Some(ff) = fields.force_field.as_ref() else {
                return false;
            };
            for &v in ff {
                match first {
                    Some(base) if v != base => return false,
                    Some(_) => {}
                    None => first = Some(v),
                }
            }
        }
        first.is_some()
    }

    fn core_has_solids(&self) -> bool {
        for (sub, fields) in self.subs.iter().zip(self.host_parts.iter()) {
            let g = sub.geom;
            for z in 0..g.core[2] {
                for y in 0..g.core[1] {
                    for x in 0..g.core[0] {
                        if fields.solid[g.pidx(x, y, z)] {
                            return true;
                        }
                    }
                }
            }
        }
        false
    }

    /// Advance one time step (V1 order: collide → stream → Bouzidi → swap →
    /// open faces → moments).
    pub fn step(&mut self) {
        self.sync_masks_if_dirty();
        let backend_gravity = self.gravity.is_some() && self.backend.supports_gravity_body_force();
        if self.gravity.is_some() && !backend_gravity {
            self.run_staged_step();
            return;
        }
        self.stage_in_if_dirty();
        let params = if backend_gravity {
            self.params_with_backend_gravity()
        } else {
            self.params.clone()
        };
        self.backend.run_span(
            &self.exchange,
            &self.subs,
            &mut self.parts,
            &params,
            self.two_pass,
            &mut self.probed_force,
            1,
        );
        self.time += 1;
        self.device_ahead = true;
        self.backend.finish_run_chunk(&self.parts, 1);
        self.refresh_probed_force();
        if !self.backend.handles_single_part_periodic_halo() {
            self.stage_out_all();
        }
    }

    /// Advance `steps` time steps.
    pub fn run(&mut self, steps: usize) {
        let backend_gravity = self.gravity.is_some() && self.backend.supports_gravity_body_force();
        if self.gravity.is_some() && !backend_gravity {
            for _ in 0..steps {
                self.step();
            }
            return;
        }
        self.sync_masks_if_dirty();
        self.stage_in_if_dirty();
        let params = if backend_gravity {
            self.params_with_backend_gravity()
        } else {
            self.params.clone()
        };
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
                &params,
                self.two_pass,
                &mut self.probed_force,
                chunk,
            );
            self.time += chunk as u64;
            self.device_ahead = true;
            self.backend.finish_run_chunk(&self.parts, chunk);
            self.refresh_probed_force();
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
        assert!(
            !on || self.backend.supports_two_pass(),
            "selected backend does not support two-pass streaming"
        );
        self.two_pass = on;
    }

    /// Whether the selected backend supports the interior/boundary streaming
    /// split used by two-pass overlap experiments.
    pub fn supports_two_pass(&self) -> bool {
        self.backend.supports_two_pass()
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

    /// Install or remove explicit Bouzidi curved-wall records for one local
    /// part. This narrow hook is intended for validation cases that construct
    /// analytic link distances directly; higher-level geometry should prefer
    /// [`Solver::set_bouzidi_circle`], [`Solver::set_bouzidi_sphere`], or
    /// [`Solver::set_bouzidi_half_way_links`].
    pub fn set_bouzidi_links(
        &mut self,
        part: usize,
        links: Option<crate::bouzidi::BouzidiLinks<T>>,
    ) {
        self.stage_out_all();
        self.host_parts[part].bouzidi = links;
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
        if self.force_field_is_uniform() {
            self.refresh_moments_after_force_change();
        }
    }

    /// Drop the per-cell body force field on every owned part (subsequent
    /// steps run force-free unless [`GlobalSpec::force`] is nonzero).
    pub fn clear_body_force_field(&mut self) {
        self.stage_out_all();
        let had_force_field = self
            .host_parts
            .iter()
            .any(|fields| fields.force_field.is_some());
        if !had_force_field {
            return;
        }
        let refresh = self.force_field_is_uniform();
        for fields in self.host_parts.iter_mut() {
            fields.force_field = None;
        }
        self.host_dirty = true;
        if refresh {
            self.refresh_moments_after_force_change();
        }
    }

    /// Add a rotating rigid-body direct-forcing IBM source to the current
    /// per-cell body-force field and return marker-level diagnostics.
    ///
    /// This method does not advance the solver. Call it after all earlier
    /// force sources for the step have been written and before [`Solver::step`]
    /// so the IBM contribution enters collision through the existing Guo
    /// forcing path. The spread increment is accumulated into any existing
    /// force field, matching the documented source-composition order.
    pub fn apply_rotating_ibm(
        &mut self,
        body: &RotatingBody,
        cfg: DirectForcingConfig,
    ) -> IbmDiagnostics {
        assert!(L::D == 2 || L::D == 3, "IBM requires a 2D or 3D lattice");
        assert!(
            cfg.max_iterations > 0,
            "IBM max_iterations must be positive"
        );
        assert!(
            cfg.slip_tolerance >= 0.0 && cfg.slip_tolerance.is_finite(),
            "IBM slip_tolerance must be finite and non-negative"
        );
        assert!(
            cfg.relaxation > 0.0 && cfg.relaxation <= 1.0 && cfg.relaxation.is_finite(),
            "IBM relaxation must be finite and in (0, 1]"
        );
        self.stage_out_all();

        let n = self.dims[0] * self.dims[1] * self.dims[2];
        let rho = self.gather_rho();
        let mut u_now = vec![[0.0f64; 3]; n];
        let mut solid = vec![false; n];
        for (sub, fields) in self.subs.iter().zip(self.host_parts.iter()) {
            let g = sub.geom;
            let np = g.n_padded();
            for z in 0..g.core[2] {
                for y in 0..g.core[1] {
                    for x in 0..g.core[0] {
                        let gi = ((sub.origin[2] + z) * self.dims[1] + (sub.origin[1] + y))
                            * self.dims[0]
                            + (sub.origin[0] + x);
                        let pi = g.pidx(x, y, z);
                        let c = g.cidx(x, y, z);
                        solid[gi] = fields.solid[pi];
                        if solid[gi] {
                            continue;
                        }
                        let mut mom = [0.0f64; 3];
                        for q in 0..L::Q {
                            let fq = fields.f[q * np + pi].as_f64();
                            for (a, ma) in mom.iter_mut().enumerate().take(L::D) {
                                *ma += L::C[q][a] as f64 * fq;
                            }
                        }
                        let mut force = [
                            self.params.force[0].as_f64(),
                            self.params.force[1].as_f64(),
                            self.params.force[2].as_f64(),
                        ];
                        if let Some(ff) = fields.force_field.as_ref() {
                            for (a, fa) in force.iter_mut().enumerate() {
                                *fa += ff[c][a].as_f64();
                            }
                        }
                        if let Some(gvec) = self.gravity {
                            for (a, fa) in force.iter_mut().enumerate() {
                                *fa += rho[gi].as_f64() * gvec[a].as_f64();
                            }
                        }
                        let inv_rho = 1.0 / rho[gi].as_f64().max(1.0e-30);
                        for a in 0..L::D {
                            u_now[gi][a] = (mom[a] + 0.5 * force[a]) * inv_rho;
                        }
                    }
                }
            }
        }

        let mut du = vec![[0.0f64; 3]; n];
        let mut spread = vec![[0.0f64; 3]; n];
        let mut diag = IbmDiagnostics::default();
        let mut max_target = 0.0f64;
        for marker in body.markers() {
            max_target = max_target.max(norm(body.target_velocity(marker.position), L::D));
        }
        let scale = max_target.max(1.0);

        let marker_stencils: Vec<_> = body
            .markers()
            .iter()
            .map(|marker| {
                let stencil = marker_stencil(marker.position, self.dims, L::D, cfg.kernel_radius);
                assert!(
                    !stencil.is_empty(),
                    "IBM marker stencil is outside the domain"
                );
                stencil
            })
            .collect();
        let mut cell_kernel_sum = vec![0.0f64; n];
        for stencil in &marker_stencils {
            for sp in stencil {
                let gi = (sp.z * self.dims[1] + sp.y) * self.dims[0] + sp.x;
                if solid[gi] {
                    continue;
                }
                cell_kernel_sum[gi] += sp.w;
            }
        }

        for iter in 0..cfg.max_iterations {
            let mut marker_impulses = vec![[0.0f64; 3]; body.markers().len()];
            for (mi, (marker, stencil)) in body
                .markers()
                .iter()
                .zip(marker_stencils.iter())
                .enumerate()
            {
                let target = body.target_velocity(marker.position);
                let mut um = [0.0f64; 3];
                let mut mobility = 0.0f64;
                for sp in stencil {
                    let gi = (sp.z * self.dims[1] + sp.y) * self.dims[0] + sp.x;
                    if solid[gi] {
                        continue;
                    }
                    let rc = rho[gi].as_f64();
                    um[0] += sp.w * (u_now[gi][0] + du[gi][0]);
                    um[1] += sp.w * (u_now[gi][1] + du[gi][1]);
                    um[2] += sp.w * (u_now[gi][2] + du[gi][2]);
                    mobility += sp.w * cell_kernel_sum[gi] / rc.max(1.0e-30);
                }
                if mobility == 0.0 {
                    continue;
                }
                let slip = [target[0] - um[0], target[1] - um[1], target[2] - um[2]];

                // Direct-forcing IBM is the marker-space linear solve M f = s,
                // where s = U_marker - I[u] and
                // M_jk = sum_cells W_j(cell) W_k(cell) / rho(cell) for marker
                // impulse unknowns q_k spread as q_k W_k. Before the post-R2-C
                // force-field impulse fix this code targeted the Guo half-force
                // predictor, sum W F/(2 rho), and therefore used 2*s/M_jj. The
                // actual post-step momentum impulse of a Guo force field is
                // F/rho, so an isolated marker is exactly corrected by s/M_jj.
                //
                // Dense markers have positive off-diagonal overlap. Using only
                // M_jj then lets a collective mode receive roughly the sum of
                // neighbouring marker corrections. The denominator used here is
                // the row sum G_j = sum_k M_jk, computed as
                // sum_cells W_j(cell) * sum_k W_k(cell) / rho(cell). Gershgorin
                // gives every eigenvalue of G^-1 M in [0, 1] because M is a
                // symmetric positive regularized-delta Gram matrix and G is its
                // positive row-sum diagonal. Richardson sweeps with
                // relaxation=1 therefore have spectral radius <= 1 for the
                // residual operator I - G^-1 M, with the unit mobility mode
                // corrected in one sweep and the remaining represented modes
                // damped. This is the Uhlmann/Wang multi-direct-forcing overlap
                // correction; for one marker G_j=M_jj, so the full-step
                // correction remains exact.
                marker_impulses[mi] = [
                    cfg.relaxation * slip[0] / mobility,
                    cfg.relaxation * slip[1] / mobility,
                    cfg.relaxation * slip[2] / mobility,
                ];
            }
            for ((marker, stencil), marker_impulse) in body
                .markers()
                .iter()
                .zip(marker_stencils.iter())
                .zip(marker_impulses.iter())
            {
                let mut represented_marker_force = [0.0f64; 3];
                for sp in stencil {
                    let gi = (sp.z * self.dims[1] + sp.y) * self.dims[0] + sp.x;
                    if solid[gi] {
                        continue;
                    }
                    let cell_force = [
                        marker_impulse[0] * sp.w,
                        marker_impulse[1] * sp.w,
                        marker_impulse[2] * sp.w,
                    ];
                    add3(&mut spread[gi], cell_force);
                    add3(&mut diag.fluid_force, cell_force);
                    add3(&mut represented_marker_force, cell_force);
                    let inv_rho = 1.0 / rho[gi].as_f64().max(1.0e-30);
                    du[gi][0] += cell_force[0] * inv_rho;
                    du[gi][1] += cell_force[1] * inv_rho;
                    du[gi][2] += cell_force[2] * inv_rho;
                }
                add3(&mut diag.marker_force, represented_marker_force);
                let r = [
                    marker.position[0] - body.center()[0],
                    marker.position[1] - body.center()[1],
                    marker.position[2] - body.center()[2],
                ];
                let tq_fluid = cross(r, represented_marker_force);
                diag.torque[0] -= tq_fluid[0];
                diag.torque[1] -= tq_fluid[1];
                diag.torque[2] -= tq_fluid[2];
            }

            diag.iterations = iter + 1;
            diag.slip_max = 0.0;
            let mut slip_sq_weighted = 0.0;
            let mut weight_sum = 0.0;
            for (marker, stencil) in body.markers().iter().zip(marker_stencils.iter()) {
                let target = body.target_velocity(marker.position);
                let mut um = [0.0f64; 3];
                let mut active_weight = false;
                for sp in stencil {
                    let gi = (sp.z * self.dims[1] + sp.y) * self.dims[0] + sp.x;
                    if solid[gi] {
                        continue;
                    }
                    active_weight = true;
                    um[0] += sp.w * (u_now[gi][0] + du[gi][0]);
                    um[1] += sp.w * (u_now[gi][1] + du[gi][1]);
                    um[2] += sp.w * (u_now[gi][2] + du[gi][2]);
                }
                if !active_weight {
                    continue;
                }
                let slip = [target[0] - um[0], target[1] - um[1], target[2] - um[2]];
                let slip_mag = norm(slip, L::D);
                diag.slip_max = diag.slip_max.max(slip_mag);
                slip_sq_weighted += marker.weight * slip_mag * slip_mag;
                weight_sum += marker.weight;
            }
            diag.slip_rms = if weight_sum > 0.0 {
                (slip_sq_weighted / weight_sum).sqrt()
            } else {
                0.0
            };
            diag.slip_max_rel = diag.slip_max / scale;
            diag.slip_rms_rel = diag.slip_rms / scale;
            let err = [
                diag.fluid_force[0] - diag.marker_force[0],
                diag.fluid_force[1] - diag.marker_force[1],
                diag.fluid_force[2] - diag.marker_force[2],
            ];
            diag.momentum_error_rel = norm(err, L::D) / norm(diag.marker_force, L::D).max(1.0e-12);
            if diag.slip_max_rel <= cfg.slip_tolerance {
                break;
            }
        }

        for (sub, fields) in self.subs.iter().zip(self.host_parts.iter_mut()) {
            let g = sub.geom;
            let ff = fields
                .force_field
                .get_or_insert_with(|| vec![[T::zero(); 3]; g.n_core()]);
            if ff.len() != g.n_core() {
                ff.clear();
                ff.resize(g.n_core(), [T::zero(); 3]);
            }
            for z in 0..g.core[2] {
                for y in 0..g.core[1] {
                    for x in 0..g.core[0] {
                        let gi = ((sub.origin[2] + z) * self.dims[1] + (sub.origin[1] + y))
                            * self.dims[0]
                            + (sub.origin[0] + x);
                        let c = g.cidx(x, y, z);
                        let add = to_real3::<T>(spread[gi]);
                        for a in 0..3 {
                            ff[c][a] = ff[c][a] + add[a];
                        }
                    }
                }
            }
        }
        self.host_dirty = true;
        diag
    }

    /// Set per-mass gravity `g`; at the start of each step, `rho(x) * g` is
    /// added to the per-cell force on fluid cells only.
    pub fn set_gravity(&mut self, g: [T; 3]) {
        self.stage_out_all();
        self.gravity = Some(g);
        if !self.core_has_solids() {
            self.refresh_moments_after_force_change();
        }
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
    // W-VOF O1 conservative Allen-Cahn phase-field transport
    // ------------------------------------------------------------------

    fn require_phase_field_lattice() -> Result<(), crate::phase_field::PhaseFieldError> {
        if L::D != 3 || L::Q != D3Q19::Q {
            return Err(crate::phase_field::PhaseFieldError {
                message: "W-VOF O1 phase-field transport is implemented only for D3Q19".to_string(),
            });
        }
        Ok(())
    }

    fn fill_phase_planes(&mut self) -> Result<(), crate::phase_field::PhaseFieldError> {
        for ((sub, fields), plane) in self
            .subs
            .iter()
            .zip(self.host_parts.iter())
            .zip(self.psi_planes.iter_mut())
        {
            let Some(phi) = fields.phi.as_ref() else {
                return Err(crate::phase_field::PhaseFieldError {
                    message: "phase field is not enabled".to_string(),
                });
            };
            let geo = sub.geom;
            plane.fill(T::zero());
            for z in 0..geo.core[2] {
                for y in 0..geo.core[1] {
                    for x in 0..geo.core[0] {
                        let pi = geo.pidx(x, y, z);
                        if !fields.solid[pi] {
                            plane[pi] = phi[geo.cidx(x, y, z)];
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
        Ok(())
    }

    fn diagnose_phase_field(&self) -> crate::phase_field::PhaseFieldDiagnostics {
        let mut total_phi = 0.0;
        let mut min_phi = f64::INFINITY;
        let mut max_phi = f64::NEG_INFINITY;
        let mut seen = false;
        for (sub, fields) in self.subs.iter().zip(self.host_parts.iter()) {
            let Some(phi) = fields.phi.as_ref() else {
                continue;
            };
            let geo = sub.geom;
            for z in 0..geo.core[2] {
                for y in 0..geo.core[1] {
                    for x in 0..geo.core[0] {
                        if fields.solid[geo.pidx(x, y, z)] {
                            continue;
                        }
                        let v = phi[geo.cidx(x, y, z)].as_f64();
                        total_phi += v;
                        if v < min_phi {
                            min_phi = v;
                        }
                        if v > max_phi {
                            max_phi = v;
                        }
                        seen = true;
                    }
                }
            }
        }
        if !seen {
            min_phi = 0.0;
            max_phi = 0.0;
        }
        crate::phase_field::PhaseFieldDiagnostics {
            total_phi,
            min_phi,
            max_phi,
        }
    }

    /// Enable the W-VOF O1 phase-field distribution from a compact global
    /// `phi` field and a prescribed velocity used for the initial `g_eq`.
    pub fn enable_phase_field_prescribed_velocity(
        &mut self,
        params: crate::phase_field::PhaseFieldParams<T>,
        phi: &[T],
        velocity: impl Fn(usize, usize, usize) -> [T; 3],
    ) -> Result<(), crate::phase_field::PhaseFieldError> {
        Self::require_phase_field_lattice()?;
        params.validate()?;
        if H::SCOPE == ExchangeScope::Remote {
            return Err(crate::phase_field::PhaseFieldError {
                message: "W-VOF O1 phase-field transport currently supports local CPU decompositions only"
                    .to_string(),
            });
        }
        let n = self.dims[0] * self.dims[1] * self.dims[2];
        if phi.len() != n {
            return Err(crate::phase_field::PhaseFieldError {
                message: format!(
                    "phase field length {} does not match cell count {n}",
                    phi.len()
                ),
            });
        }
        self.stage_out_all();
        for (sub, fields) in self.subs.iter().zip(self.host_parts.iter_mut()) {
            let geo = sub.geom;
            let np = geo.n_padded();
            let nc = geo.n_core();
            let local_phi = fields.phi.get_or_insert_with(|| vec![T::zero(); nc]);
            if local_phi.len() != nc {
                local_phi.resize(nc, T::zero());
            }
            let gset = fields
                .g
                .get_or_insert_with(|| vec![T::zero(); D3Q19::Q * np]);
            if gset.len() != D3Q19::Q * np {
                gset.resize(D3Q19::Q * np, T::zero());
            }
            let gtmp = fields
                .gtmp
                .get_or_insert_with(|| vec![T::zero(); D3Q19::Q * np]);
            if gtmp.len() != D3Q19::Q * np {
                gtmp.resize(D3Q19::Q * np, T::zero());
            }
            gset.fill(T::zero());
            gtmp.fill(T::zero());
            for z in 0..geo.core[2] {
                for y in 0..geo.core[1] {
                    for x in 0..geo.core[0] {
                        let gi = ((sub.origin[2] + z) * self.dims[1] + (sub.origin[1] + y))
                            * self.dims[0]
                            + (sub.origin[0] + x);
                        let c = geo.cidx(x, y, z);
                        let pi = geo.pidx(x, y, z);
                        local_phi[c] = phi[gi];
                        if fields.solid[pi] {
                            continue;
                        }
                        let u = velocity(sub.origin[0] + x, sub.origin[1] + y, sub.origin[2] + z);
                        let geq = crate::phase_field::equilibrium(phi[gi], u);
                        for q in 0..D3Q19::Q {
                            gset[q * np + pi] = geq[q];
                        }
                    }
                }
            }
        }
        self.host_dirty = true;
        Ok(())
    }

    /// Advance only the W-VOF O1 phase-field distribution under a prescribed
    /// velocity field. The hydrodynamic `f` step is not executed.
    pub fn phase_field_step_prescribed_velocity(
        &mut self,
        params: crate::phase_field::PhaseFieldParams<T>,
        velocity: impl Fn(usize, usize, usize) -> [T; 3],
    ) -> Result<crate::phase_field::PhaseFieldDiagnostics, crate::phase_field::PhaseFieldError>
    {
        Self::require_phase_field_lattice()?;
        let params = params.validate()?;
        if H::SCOPE == ExchangeScope::Remote {
            return Err(crate::phase_field::PhaseFieldError {
                message: "W-VOF O1 phase-field transport currently supports local CPU decompositions only"
                    .to_string(),
            });
        }
        self.stage_out_all();
        self.fill_phase_planes()?;

        let omega = params.omega();
        for (part, (sub, fields)) in self.subs.iter().zip(self.host_parts.iter_mut()).enumerate() {
            let geo = sub.geom;
            let np = geo.n_padded();
            let Some(gset) = fields.g.as_mut() else {
                return Err(crate::phase_field::PhaseFieldError {
                    message: "phase-field g distribution is not enabled".to_string(),
                });
            };
            let Some(phi) = fields.phi.as_ref() else {
                return Err(crate::phase_field::PhaseFieldError {
                    message: "phase field is not enabled".to_string(),
                });
            };
            let phi_plane = &self.psi_planes[part];
            for z in 0..geo.core[2] {
                for y in 0..geo.core[1] {
                    for x in 0..geo.core[0] {
                        let pi = geo.pidx(x, y, z);
                        if fields.solid[pi] {
                            continue;
                        }
                        let c = geo.cidx(x, y, z);
                        let p = phi[c];
                        let u = velocity(sub.origin[0] + x, sub.origin[1] + y, sub.origin[2] + z);
                        let geq = crate::phase_field::equilibrium(p, u);
                        let (grad, _) = crate::phase_field::grad_lap(geo, phi_plane, x, y, z);
                        let mut source = [T::zero(); 19];
                        let mut source_sum = T::zero();
                        for q in 0..D3Q19::Q {
                            source[q] = crate::phase_field::collide_source(params, p, grad, q);
                            source_sum = source_sum + source[q];
                        }
                        source[D3Q19::REST] = source[D3Q19::REST] - source_sum;
                        for q in 0..D3Q19::Q {
                            let idx = q * np + pi;
                            let old = gset[idx];
                            // The source first moment is multiplied by
                            // tau_phi in the recovered flux, so the discrete
                            // update uses omega_phi*M to recover the
                            // governing M(4/W)phi(1-phi)n_hat counter-flux.
                            gset[idx] =
                                old - omega * (old - geq[q]) + omega * params.mobility * source[q];
                        }
                    }
                }
            }
        }

        exchange_g_generic::<D3Q19, T>(&self.subs, &mut self.host_parts);

        for (sub, fields) in self.subs.iter_mut().zip(self.host_parts.iter_mut()) {
            let geo = sub.geom;
            let np = geo.n_padded();
            {
                let gset = fields
                    .g
                    .as_ref()
                    .expect("phase-field g distribution must be enabled");
                let gout = fields
                    .gtmp
                    .as_mut()
                    .expect("phase-field g ping-pong distribution must be enabled");
                for z in 0..geo.core[2] {
                    for y in 0..geo.core[1] {
                        for x in 0..geo.core[0] {
                            let dst = geo.pidx(x, y, z);
                            if fields.solid[dst] {
                                continue;
                            }
                            for q in 0..D3Q19::Q {
                                let c = D3Q19::C[q];
                                let src = geo.pidx_i(
                                    x as isize - c[0] as isize,
                                    y as isize - c[1] as isize,
                                    z as isize - c[2] as isize,
                                );
                                gout[q * np + dst] = gset[q * np + src];
                            }
                        }
                    }
                }
            }
            std::mem::swap(
                fields.g.as_mut().expect("phase-field g distribution"),
                fields
                    .gtmp
                    .as_mut()
                    .expect("phase-field g ping-pong distribution"),
            );
            let gset = fields.g.as_ref().expect("phase-field g distribution");
            let phi = fields.phi.as_mut().expect("phase field");
            for z in 0..geo.core[2] {
                for y in 0..geo.core[1] {
                    for x in 0..geo.core[0] {
                        let pi = geo.pidx(x, y, z);
                        let c = geo.cidx(x, y, z);
                        if fields.solid[pi] {
                            phi[c] = T::zero();
                            continue;
                        }
                        let mut populations = [T::zero(); 19];
                        for q in 0..D3Q19::Q {
                            populations[q] = gset[q * np + pi];
                        }
                        phi[c] = crate::phase_field::sum_populations(populations);
                    }
                }
            }
        }
        self.host_dirty = true;
        Ok(self.diagnose_phase_field())
    }

    /// Run the phase-field prescribed-velocity pre-pass, then the existing
    /// hydrodynamic step. The hydrodynamic pass order is unchanged.
    pub fn step_with_phase_field_prescribed_velocity(
        &mut self,
        params: crate::phase_field::PhaseFieldParams<T>,
        velocity: impl Fn(usize, usize, usize) -> [T; 3],
    ) -> Result<crate::phase_field::PhaseFieldDiagnostics, crate::phase_field::PhaseFieldError>
    {
        let diag = self.phase_field_step_prescribed_velocity(params, velocity)?;
        self.step();
        Ok(diag)
    }

    /// Compute the conservative Allen-Cahn interface flux `J_phi` at one
    /// global cell from the current `phi` state.
    pub fn phase_flux_jphi(
        &mut self,
        params: crate::phase_field::PhaseFieldParams<T>,
        x: usize,
        y: usize,
        z: usize,
    ) -> Result<[T; 3], crate::phase_field::PhaseFieldError> {
        Self::require_phase_field_lattice()?;
        let params = params.validate()?;
        self.stage_out_all();
        self.fill_phase_planes()?;
        let (part, lx, ly, lz) = self.locate(x, y, z);
        let fields = &self.host_parts[part];
        let geo = fields.geom;
        let phi = fields
            .phi
            .as_ref()
            .ok_or_else(|| crate::phase_field::PhaseFieldError {
                message: "phase field is not enabled".to_string(),
            })?;
        let (grad, _) = crate::phase_field::grad_lap(geo, &self.psi_planes[part], lx, ly, lz);
        Ok(crate::phase_field::phase_flux_jphi(
            params,
            phi[geo.cidx(lx, ly, lz)],
            grad,
        ))
    }

    /// Gather the compact global `phi` field.
    pub fn gather_phi(&mut self) -> Result<Vec<T>, crate::phase_field::PhaseFieldError> {
        self.stage_out_all();
        let mut out = vec![T::zero(); self.dims[0] * self.dims[1] * self.dims[2]];
        for (sub, fields) in self.subs.iter().zip(self.host_parts.iter()) {
            let Some(phi) = fields.phi.as_ref() else {
                return Err(crate::phase_field::PhaseFieldError {
                    message: "phase field is not enabled".to_string(),
                });
            };
            let geo = sub.geom;
            for z in 0..geo.core[2] {
                for y in 0..geo.core[1] {
                    for x in 0..geo.core[0] {
                        let gi = ((sub.origin[2] + z) * self.dims[1] + (sub.origin[1] + y))
                            * self.dims[0]
                            + (sub.origin[0] + x);
                        out[gi] = phi[geo.cidx(x, y, z)];
                    }
                }
            }
        }
        Ok(out)
    }

    /// Export `phi` along an x-directed line as CSV.
    pub fn export_phi_x_profile_csv(
        &mut self,
        y: usize,
        z: usize,
        path: impl AsRef<Path>,
    ) -> Result<(), Box<dyn std::error::Error>> {
        let phi = self.gather_phi()?;
        let path = path.as_ref();
        if let Some(parent) = path.parent() {
            std::fs::create_dir_all(parent)?;
        }
        let mut file = std::fs::File::create(path)?;
        writeln!(file, "x,phi")?;
        for x in 0..self.dims[0] {
            let gi = (z * self.dims[1] + y) * self.dims[0] + x;
            writeln!(file, "{x},{}", phi[gi])?;
        }
        Ok(())
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
    /// Mutable backend-owned fields of part `i` for backend-specific shims.
    pub fn backend_fields_mut(&mut self, i: usize) -> &mut B::Fields {
        &mut self.parts[i]
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
    /// Mutable fields of part `i` for crate-internal setup and fault-injection
    /// tests. Public callers must use dedicated methods so host/mask dirty
    /// flags cannot be forgotten.
    #[allow(dead_code)]
    pub(crate) fn fields_mut(&mut self, i: usize) -> &mut SoaFields<T> {
        self.stage_out_all();
        self.host_dirty = true;
        &mut self.host_parts[i]
    }

    #[cfg(debug_assertions)]
    pub(crate) fn host_dirty_for_debug(&self) -> bool {
        self.host_dirty
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
        self.host_dirty = true;
        for (sub, fields) in self.subs.iter().zip(self.host_parts.iter_mut()) {
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

    /// Install a compact global per-cell body-force field.
    ///
    /// `values` is indexed as `((z * ny + y) * nx + x)` and is sliced into
    /// owned parts automatically. This is the slice-oriented counterpart of
    /// [`Solver::set_body_force_field`]; both mark the host staging as dirty
    /// so the next step uploads the edited force field before collision.
    pub fn set_body_force_field_values(&mut self, values: &[[T; 3]]) {
        let n = self.dims[0] * self.dims[1] * self.dims[2];
        assert_eq!(values.len(), n, "force field length must match cell count");
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
                        let gi = ((sub.origin[2] + z) * self.dims[1] + (sub.origin[1] + y))
                            * self.dims[0]
                            + (sub.origin[0] + x);
                        buf[g.cidx(x, y, z)] = values[gi];
                    }
                }
            }
        }
        self.host_dirty = true;
        if self.force_field_is_uniform() {
            self.refresh_moments_after_force_change();
        }
    }

    fn write_part_checkpoint(
        dir: &Path,
        rank: usize,
        part: usize,
        sub: &Subdomain,
        fields: &SoaFields<T>,
    ) -> Result<CheckpointRank, CheckpointError> {
        let np = fields.geom.n_padded();
        let nc = fields.geom.n_core();

        let mut sections: Vec<(u32, Vec<u8>)> = Vec::new();
        let mut f_primary = Vec::with_capacity(fields.f.len() * std::mem::size_of::<T>());
        push_real_bytes(&mut f_primary, &fields.f);
        sections.push((SEC_F_PRIMARY, f_primary));

        let mut ftmp = Vec::with_capacity(fields.ftmp.len() * std::mem::size_of::<T>());
        push_real_bytes(&mut ftmp, &fields.ftmp);
        sections.push((SEC_STALE_STASH, ftmp));

        let mut moments = Vec::with_capacity(4 * nc * std::mem::size_of::<T>());
        push_real_bytes(&mut moments, &fields.rho);
        push_real_bytes(&mut moments, &fields.ux);
        push_real_bytes(&mut moments, &fields.uy);
        push_real_bytes(&mut moments, &fields.uz);
        sections.push((SEC_MOMENTS, moments));

        let solid: Vec<u8> = fields.solid.iter().map(|&v| u8::from(v)).collect();
        sections.push((SEC_SOLID, solid));

        if let Some(force) = &fields.force_field {
            let mut bytes = Vec::with_capacity(3 * force.len() * std::mem::size_of::<T>());
            for v in force {
                push_real_bytes(&mut bytes, v);
            }
            sections.push((SEC_FORCE_FIELD, bytes));
        }

        let mut payload = Vec::new();
        let mut table = Vec::with_capacity(sections.len());
        for (id, bytes) in sections {
            let offset = payload.len() as u64;
            payload.extend_from_slice(&bytes);
            table.push(SectionEntry {
                id,
                offset,
                byte_len: bytes.len() as u64,
            });
        }
        let payload_hash = fnv1a64(&payload);

        let mut bin = Vec::new();
        bin.extend_from_slice(CKPT_MAGIC);
        write_u32(&mut bin, CKPT_FORMAT_VERSION);
        bin.push(0);
        bin.push(if dtype_name::<T>() == "f32" { 0 } else { 1 });
        write_u16(&mut bin, lattice_id::<L>());
        write_u16(&mut bin, L::Q as u16);
        write_u16(&mut bin, L::D as u16);
        write_u64(&mut bin, np as u64);
        write_u64(&mut bin, nc as u64);
        write_u32(&mut bin, table.len() as u32);
        for entry in &table {
            write_u32(&mut bin, entry.id);
            write_u64(&mut bin, entry.offset);
            write_u64(&mut bin, entry.byte_len);
        }
        bin.extend_from_slice(&payload);
        write_u64(&mut bin, payload_hash);

        let rank_file = format!("rank_{rank:04}_part_{part:04}.bin");
        let rank_path = dir.join(&rank_file);
        let mut file = std::fs::File::create(&rank_path)?;
        file.write_all(&bin)?;

        Ok(CheckpointRank {
            rank,
            part,
            file: rank_file,
            origin: sub.origin,
            core: sub.geom.core,
            bytes: bin.len() as u64,
            payload_hash: hash_string(payload_hash),
            mask_hash: hash_string(part_mask_hash(fields)),
        })
    }

    pub(crate) fn checkpoint_manifest(
        &self,
        nranks: usize,
        ranks: Vec<CheckpointRank>,
        decomp_hash: u64,
    ) -> CheckpointManifest {
        let mut reserved = BTreeMap::new();
        reserved.insert("rng".to_string(), false);
        reserved.insert("particles".to_string(), false);
        reserved.insert("stats".to_string(), false);
        CheckpointManifest {
            kind: "lbmflow-checkpoint".to_string(),
            format_version: CKPT_FORMAT_VERSION,
            step: self.time,
            time: self.time as f64,
            dtype: dtype_name::<T>().to_string(),
            lattice: lattice_name::<L>().to_string(),
            global: self.dims,
            scenario_hash: hash_string(self.current_spec_hash()),
            decomp_hash: hash_string(decomp_hash),
            nranks,
            ranks,
            reserved,
        }
    }

    pub(crate) fn write_checkpoint_manifest(
        dir: impl AsRef<Path>,
        manifest: &CheckpointManifest,
    ) -> Result<(), CheckpointError> {
        std::fs::write(
            dir.as_ref().join("manifest.json"),
            serde_json::to_string_pretty(&manifest)?,
        )?;
        Ok(())
    }

    pub(crate) fn save_owned_parts_for_checkpoint(
        &mut self,
        dir: impl AsRef<Path>,
        rank: usize,
        part_ids: &[usize],
    ) -> Result<Vec<CheckpointRank>, CheckpointError> {
        self.stage_out_all();
        if part_ids.len() != self.host_parts.len() {
            return Err(CheckpointError::new(
                "CKPT_DECOMP_MISMATCH",
                format!(
                    "checkpoint part id list has {} entries for {} owned parts",
                    part_ids.len(),
                    self.host_parts.len()
                ),
            ));
        }
        let dir = dir.as_ref();
        std::fs::create_dir_all(dir)?;
        let mut out = Vec::with_capacity(self.host_parts.len());
        for ((sub, fields), &part) in self
            .subs
            .iter()
            .zip(self.host_parts.iter())
            .zip(part_ids.iter())
        {
            out.push(Self::write_part_checkpoint(dir, rank, part, sub, fields)?);
        }
        Ok(out)
    }

    /// Save a checkpoint directory (`manifest.json` plus one payload per
    /// owned part). Single-process split runs keep all part payloads under
    /// rank 0 and require the same partition layout on restart.
    pub fn save(&mut self, dir: impl AsRef<Path>) -> Result<(), CheckpointError> {
        let dir = dir.as_ref();
        let part_ids: Vec<usize> = (0..self.host_parts.len()).collect();
        let ranks = self.save_owned_parts_for_checkpoint(dir, 0, &part_ids)?;
        let manifest = self.checkpoint_manifest(1, ranks, decomp_hash(&self.subs, self.periodic));
        Self::write_checkpoint_manifest(dir, &manifest)
    }

    /// Load a single-rank checkpoint into a freshly rebuilt solver. The caller
    /// supplies the same spec/masks/decomposition used to build the run; hashes
    /// are checked before any checkpoint bytes are installed.
    pub fn load(
        dir: impl AsRef<Path>,
        spec: &GlobalSpec<T>,
        solid: &[bool],
        wall_u: &[[T; 3]],
        decomp: [usize; 3],
        backend: B,
        exchange: H,
    ) -> Result<Self, CheckpointError> {
        let mut solver = Self::new(spec, solid, wall_u, decomp, backend, exchange);
        solver.load_into(dir, spec_hash::<T, L>(spec, solid, wall_u))?;
        Ok(solver)
    }

    /// Restore a checkpoint into this already-built solver. The current
    /// solver configuration and masks are used as the rebuilt spec/mask
    /// authority; mismatch errors are returned before checkpoint state is
    /// installed.
    pub fn restore(&mut self, dir: impl AsRef<Path>) -> Result<(), CheckpointError> {
        let expected = self.current_spec_hash();
        self.load_into(dir, expected)
    }

    fn current_spec_hash(&self) -> u64 {
        let spec = GlobalSpec {
            dims: self.dims,
            nu: self.nu,
            collision: self.collision,
            periodic: self.periodic,
            faces: self.params.faces,
            force: self.params.force,
            sources: self.params.sources.clone(),
            face_patches: self.params.face_patches.clone(),
        };
        let fields = &self.host_parts[0];
        spec_hash::<T, L>(&spec, &fields.solid, &fields.wall_u)
    }

    fn load_into(
        &mut self,
        dir: impl AsRef<Path>,
        expected_spec_hash: u64,
    ) -> Result<(), CheckpointError> {
        let dir = dir.as_ref();
        let manifest_text = std::fs::read_to_string(dir.join("manifest.json"))?;
        let manifest: CheckpointManifest = serde_json::from_str(&manifest_text)?;
        if manifest.kind != "lbmflow-checkpoint" {
            return Err(CheckpointError::new(
                "CKPT_BAD_MAGIC",
                format!("manifest kind is {}", manifest.kind),
            ));
        }
        if manifest.format_version != CKPT_FORMAT_VERSION {
            return Err(CheckpointError::new(
                "CKPT_VERSION_MISMATCH",
                format!(
                    "checkpoint format_version {} differs from supported {}",
                    manifest.format_version, CKPT_FORMAT_VERSION
                ),
            ));
        }
        if manifest.dtype != dtype_name::<T>() {
            return Err(CheckpointError::new(
                "CKPT_DTYPE_MISMATCH",
                format!(
                    "checkpoint dtype {} differs from run dtype {}",
                    manifest.dtype,
                    dtype_name::<T>()
                ),
            ));
        }
        if manifest.lattice != lattice_name::<L>() {
            return Err(CheckpointError::new(
                "CKPT_LATTICE_MISMATCH",
                format!(
                    "checkpoint lattice {} differs from run lattice {}",
                    manifest.lattice,
                    lattice_name::<L>()
                ),
            ));
        }
        if manifest.global != self.dims {
            return Err(CheckpointError::new(
                "CKPT_GEOM_MISMATCH",
                format!(
                    "checkpoint global {:?} differs from run {:?}",
                    manifest.global, self.dims
                ),
            ));
        }
        let expected = hash_string(expected_spec_hash);
        if manifest.scenario_hash != expected {
            return Err(CheckpointError::new(
                "CKPT_SCENARIO_MISMATCH",
                format!(
                    "scenario_hash differed: checkpoint={} current={expected}",
                    manifest.scenario_hash
                ),
            ));
        }
        let current_decomp = hash_string(decomp_hash(&self.subs, self.periodic));
        self.load_into_with_layout(dir, manifest, current_decomp, 1, 0)
    }

    #[cfg_attr(not(feature = "mpi"), allow(dead_code))]
    pub(crate) fn restore_distributed_checkpoint(
        &mut self,
        dir: impl AsRef<Path>,
        expected_decomp_hash: u64,
        nranks: usize,
        rank: usize,
    ) -> Result<(), CheckpointError> {
        let expected = self.current_spec_hash();
        let dir = dir.as_ref();
        let manifest_text = std::fs::read_to_string(dir.join("manifest.json"))?;
        let manifest: CheckpointManifest = serde_json::from_str(&manifest_text)?;
        if manifest.kind != "lbmflow-checkpoint" {
            return Err(CheckpointError::new(
                "CKPT_BAD_MAGIC",
                format!("manifest kind is {}", manifest.kind),
            ));
        }
        if manifest.format_version != CKPT_FORMAT_VERSION {
            return Err(CheckpointError::new(
                "CKPT_VERSION_MISMATCH",
                format!(
                    "checkpoint format_version {} differs from supported {}",
                    manifest.format_version, CKPT_FORMAT_VERSION
                ),
            ));
        }
        if manifest.dtype != dtype_name::<T>() {
            return Err(CheckpointError::new(
                "CKPT_DTYPE_MISMATCH",
                format!(
                    "checkpoint dtype {} differs from run dtype {}",
                    manifest.dtype,
                    dtype_name::<T>()
                ),
            ));
        }
        if manifest.lattice != lattice_name::<L>() {
            return Err(CheckpointError::new(
                "CKPT_LATTICE_MISMATCH",
                format!(
                    "checkpoint lattice {} differs from run lattice {}",
                    manifest.lattice,
                    lattice_name::<L>()
                ),
            ));
        }
        if manifest.global != self.dims {
            return Err(CheckpointError::new(
                "CKPT_GEOM_MISMATCH",
                format!(
                    "checkpoint global {:?} differs from run {:?}",
                    manifest.global, self.dims
                ),
            ));
        }
        let expected = hash_string(expected);
        if manifest.scenario_hash != expected {
            return Err(CheckpointError::new(
                "CKPT_SCENARIO_MISMATCH",
                format!(
                    "scenario_hash differed: checkpoint={} current={expected}",
                    manifest.scenario_hash
                ),
            ));
        }
        let current_decomp = hash_string(expected_decomp_hash);
        self.load_into_with_layout(dir, manifest, current_decomp, nranks, rank)
    }

    fn load_into_with_layout(
        &mut self,
        dir: &Path,
        manifest: CheckpointManifest,
        current_decomp: String,
        expected_nranks: usize,
        current_rank: usize,
    ) -> Result<(), CheckpointError> {
        if manifest.decomp_hash != current_decomp {
            return Err(CheckpointError::new(
                "CKPT_DECOMP_MISMATCH",
                format!(
                    "decomp_hash differed: checkpoint={} current={current_decomp}",
                    manifest.decomp_hash
                ),
            ));
        }
        if manifest.nranks != expected_nranks {
            return Err(CheckpointError::new(
                "CKPT_RANK_MISMATCH",
                format!(
                    "checkpoint rank count {} differs from current {}",
                    manifest.nranks, expected_nranks
                ),
            ));
        }
        let mut entries: Vec<&CheckpointRank> = manifest
            .ranks
            .iter()
            .filter(|entry| entry.rank == current_rank)
            .collect();
        entries.sort_by_key(|entry| entry.part);
        if entries.len() != self.host_parts.len() {
            return Err(CheckpointError::new(
                "CKPT_DECOMP_MISMATCH",
                format!(
                    "checkpoint has {} payload parts for rank {}, current owner has {}",
                    entries.len(),
                    current_rank,
                    self.host_parts.len()
                ),
            ));
        }

        for ((entry, sub), fields) in entries
            .into_iter()
            .zip(self.subs.iter())
            .zip(self.host_parts.iter_mut())
        {
            if entry.origin != sub.origin || entry.core != sub.geom.core {
                return Err(CheckpointError::new(
                    "CKPT_GEOM_MISMATCH",
                    format!(
                        "part geometry differed: checkpoint rank={} part={} origin={:?} core={:?}, current origin={:?} core={:?}",
                        entry.rank, entry.part, entry.origin, entry.core, sub.origin, sub.geom.core
                    ),
                ));
            }
            let current_mask = hash_string(part_mask_hash(fields));
            if entry.mask_hash != current_mask {
                return Err(CheckpointError::new(
                    "CKPT_MASK_MISMATCH",
                    format!(
                        "mask_hash differed for rank={} part={}: checkpoint={} current={current_mask}",
                        entry.rank, entry.part, entry.mask_hash
                    ),
                ));
            }

            let (header, payload) = read_rank_file(&dir.join(&entry.file))?;
            if header.dtype != (if dtype_name::<T>() == "f32" { 0 } else { 1 }) {
                return Err(CheckpointError::new(
                    "CKPT_DTYPE_MISMATCH",
                    "rank file dtype differs from requested run precision",
                ));
            }
            if header.lattice_id != lattice_id::<L>()
                || header.q != L::Q as u16
                || header.d != L::D as u16
            {
                return Err(CheckpointError::new(
                    "CKPT_LATTICE_MISMATCH",
                    "rank file lattice metadata differs from requested lattice",
                ));
            }
            let np = fields.geom.n_padded();
            let nc = fields.geom.n_core();
            if header.np != np as u64 || header.n_core != nc as u64 {
                return Err(CheckpointError::new(
                    "CKPT_GEOM_MISMATCH",
                    format!(
                        "rank file np/n_core {}/{} differs from current {np}/{nc}",
                        header.np, header.n_core
                    ),
                ));
            }
            if hash_string(header.payload_hash) != entry.payload_hash {
                return Err(CheckpointError::new(
                    "CKPT_PAYLOAD_CORRUPT",
                    format!(
                        "rank payload hash {} differs from manifest {}",
                        hash_string(header.payload_hash),
                        entry.payload_hash
                    ),
                ));
            }

            let f = required_section(&payload, SEC_F_PRIMARY, "F_PRIMARY")?;
            fields.f = read_real_bytes(f, L::Q * np)?;
            let ftmp = required_section(&payload, SEC_STALE_STASH, "STALE_STASH")?;
            fields.ftmp = read_real_bytes(ftmp, L::Q * np)?;
            let moments = required_section(&payload, SEC_MOMENTS, "MOMENTS")?;
            let vals = read_real_bytes(moments, 4 * nc)?;
            fields.rho.copy_from_slice(&vals[0..nc]);
            fields.ux.copy_from_slice(&vals[nc..2 * nc]);
            fields.uy.copy_from_slice(&vals[2 * nc..3 * nc]);
            fields.uz.copy_from_slice(&vals[3 * nc..4 * nc]);
            if let Some(solid) = payload.get(&SEC_SOLID) {
                if solid.len() != fields.solid.len() {
                    return Err(CheckpointError::new(
                        "CKPT_GEOM_MISMATCH",
                        "solid section length differs from current padded geometry",
                    ));
                }
                let loaded: Vec<bool> = solid.iter().map(|&v| v != 0).collect();
                if loaded != fields.solid {
                    return Err(CheckpointError::new(
                        "CKPT_MASK_MISMATCH",
                        "serialized solid mask differs from rebuilt mask",
                    ));
                }
            }
            if let Some(force) = payload.get(&SEC_FORCE_FIELD) {
                let vals = read_real_bytes(force, 3 * nc)?;
                let mut ff = vec![[T::zero(); 3]; nc];
                for c in 0..nc {
                    ff[c] = [vals[3 * c], vals[3 * c + 1], vals[3 * c + 2]];
                }
                fields.force_field = Some(ff);
            } else {
                fields.force_field = None;
            }
            fields.fused = None;
        }
        self.time = manifest.step;
        self.probed_force = [T::zero(); 3];
        self.host_dirty = true;
        self.device_ahead = false;
        self.masks_dirty = false;
        Ok(())
    }

    /// Momentum-exchange force on the probed solids during the most recent
    /// step (V1 `probed_force`).
    pub fn probed_force(&self) -> [T; 3] {
        self.probed_force
    }

    /// Explicit backend readback of the momentum-exchange force on probed
    /// solids during the most recent completed step.
    pub fn read_probed_force(&self) -> [T; 3] {
        self.parts.iter().fold([T::zero(); 3], |a, field| {
            let b = self.backend.read_probed_force(field);
            [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
        })
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
        T::r(self.total_mass_f64())
    }

    /// Total mass over fluid cells as an `f64` diagnostic.
    ///
    /// Population storage may be `f32`, but the reduction is still performed
    /// in `f64`; returning `f64` avoids quantizing million-cell diagnostics
    /// back to the current precision's scalar spacing.
    pub fn total_mass_f64(&self) -> f64 {
        let (fluid, m) = self.local_mass_partials();
        fluid + m
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
        let params_with_gravity;
        let backend_gravity = self.gravity.is_some() && self.backend.supports_gravity_body_force();
        let params = if backend_gravity {
            params_with_gravity = self.params_with_backend_gravity();
            &params_with_gravity
        } else {
            &self.params
        };
        let mut p = [0.0f64; 3];
        for (sub, fields) in self.subs.iter().zip(self.parts.iter()) {
            for (a, pa) in p.iter_mut().enumerate() {
                *pa += self
                    .backend
                    .reduce(sub, fields, params, Reduction::Momentum(a));
            }
        }
        if !backend_gravity {
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

    /// Wall metrics for wall-adjacent fluid cells in global compact order.
    ///
    /// This is a read-only W1 diagnostic. It reads the current velocity moments
    /// and wall geometry, but it does not write populations, moments, masks, or
    /// per-cell relaxation fields. Half-way rim cells use `y_w = 0.5`; Bouzidi
    /// cells use the shortest installed `qd * |c_q|` link for that fluid cell.
    pub fn gather_wall_metrics(&self) -> Vec<WallCellMetric<T>> {
        let mut out = Vec::new();
        for ((sub, host), fields) in self
            .subs
            .iter()
            .zip(self.host_parts.iter())
            .zip(self.parts.iter())
        {
            let g = sub.geom;
            let mut hm = HostMoments::default();
            self.backend.read_moments(fields, &mut hm);

            let mut bouzidi = vec![None; g.n_padded()];
            if let Some(links) = &host.bouzidi {
                for rec in &links.records {
                    let q = rec.q as usize;
                    assert!(q > 0 && q < L::Q, "invalid Bouzidi direction {q}");
                    let cq = L::C[q];
                    let len_sq = cq[0] as f64 * cq[0] as f64
                        + cq[1] as f64 * cq[1] as f64
                        + cq[2] as f64 * cq[2] as f64;
                    assert!(len_sq > 0.0, "Bouzidi wall link must be non-rest");
                    let len = T::r(len_sq).sqrt();
                    let y_w = rec.qd * len;
                    let inv_len = T::one() / len;
                    let normal = [
                        T::r(cq[0] as f64) * inv_len,
                        T::r(cq[1] as f64) * inv_len,
                        T::r(cq[2] as f64) * inv_len,
                    ];
                    let wall_u = host.wall_u[rec.wall_ref as usize];
                    let cell = rec.cell as usize;
                    let replace = match bouzidi[cell] {
                        Some((prev_y_w, _, _)) => y_w < prev_y_w,
                        None => true,
                    };
                    if replace {
                        bouzidi[cell] = Some((y_w, normal, wall_u));
                    }
                }
            }

            for z in 0..g.core[2] {
                for y in 0..g.core[1] {
                    for x in 0..g.core[0] {
                        let pi = g.pidx(x, y, z);
                        if host.solid[pi] {
                            continue;
                        }
                        let c = g.cidx(x, y, z);
                        let u = [hm.ux[c], hm.uy[c], hm.uz[c]];
                        let gi = ((sub.origin[2] + z) * self.dims[1] + (sub.origin[1] + y))
                            * self.dims[0]
                            + (sub.origin[0] + x);

                        let metric_input = if let Some((y_w, normal, wall_u)) = bouzidi[pi] {
                            Some((y_w, normal, wall_u, WallMetricSource::Bouzidi))
                        } else {
                            self.halfway_wall_metric_input(host, x, y, z)
                                .map(|(normal, wall_u)| {
                                    (T::r(0.5), normal, wall_u, WallMetricSource::HalfwayRim)
                                })
                        };

                        let Some((y_w, normal, wall_u, source)) = metric_input else {
                            continue;
                        };
                        let u_parallel = tangential_speed(u, wall_u, normal);
                        let u_tau = friction_velocity(u_parallel, y_w, T::r(self.nu));
                        out.push(WallCellMetric {
                            cell_index: gi,
                            y_w,
                            u_parallel,
                            u_tau,
                            y_plus: y_w * u_tau / T::r(self.nu),
                            tau_w: u_tau * u_tau,
                            source,
                        });
                    }
                }
            }
        }
        out.sort_by_key(|m| m.cell_index);
        out
    }

    fn halfway_wall_metric_input(
        &self,
        fields: &SoaFields<T>,
        x: usize,
        y: usize,
        z: usize,
    ) -> Option<([T; 3], [T; 3])> {
        let g = fields.geom;
        let mut normal = [T::zero(); 3];
        let mut wall_u = [T::zero(); 3];
        let mut n = 0usize;
        for q in 1..L::Q {
            let cq = L::C[q];
            let np = g.pidx_i(
                x as isize + cq[0] as isize,
                y as isize + cq[1] as isize,
                z as isize + cq[2] as isize,
            );
            if fields.solid[np] {
                for a in 0..3 {
                    normal[a] = normal[a] + T::r(cq[a] as f64);
                    wall_u[a] = wall_u[a] + fields.wall_u[np][a];
                }
                n += 1;
            }
        }
        if n == 0 {
            return None;
        }
        let norm_sq = normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2];
        assert!(
            norm_sq > T::zero(),
            "half-way wall normal is undefined because solid-link directions cancel"
        );
        let inv_norm = T::one() / norm_sq.sqrt();
        let inv_n = T::one() / T::r(n as f64);
        for a in 0..3 {
            normal[a] = normal[a] * inv_norm;
            wall_u[a] = wall_u[a] * inv_n;
        }
        Some((normal, wall_u))
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
        let params_with_gravity;
        let params = if self.gravity.is_some() {
            params_with_gravity = self.params_with_backend_gravity();
            &params_with_gravity
        } else {
            &self.params
        };
        let kp = KParams::new::<L>(params);
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
        let force = kp.force_at(fields.force_field.as_deref(), c, r);
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
            if self.host_parts[i].solid[pi] {
                Some((self.host_parts[i].wall_u[pi], true))
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
    use std::time::{SystemTime, UNIX_EPOCH};

    fn tmp_ckpt(name: &str) -> std::path::PathBuf {
        let n = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        std::env::temp_dir().join(format!("lbmflow-{name}-{n}"))
    }

    fn bits<T: Real>(v: T) -> u64 {
        if std::mem::size_of::<T>() == 4 {
            (v.as_f64() as f32).to_bits() as u64
        } else {
            v.as_f64().to_bits()
        }
    }

    fn assert_solver_bits_eq<L, T, HA, HB>(
        a: &Solver<L, T, CpuScalar, HA>,
        b: &Solver<L, T, CpuScalar, HB>,
    ) where
        L: Lattice,
        T: Real,
        HA: HaloExchange<T>,
        HB: HaloExchange<T>,
    {
        assert_eq!(a.time(), b.time());
        assert_eq!(a.dims(), b.dims());
        for q in 0..L::Q {
            let af = a.gather_f(q);
            let bf = b.gather_f(q);
            assert_eq!(af.len(), bf.len());
            for (i, (&x, &y)) in af.iter().zip(&bf).enumerate() {
                assert_eq!(bits(x), bits(y), "f[{q}][{i}] differs");
            }
        }
        for (name, av, bv) in [
            ("rho", a.gather_rho(), b.gather_rho()),
            ("ux", a.gather_ux(), b.gather_ux()),
            ("uy", a.gather_uy(), b.gather_uy()),
            ("uz", a.gather_uz(), b.gather_uz()),
        ] {
            for (i, (&x, &y)) in av.iter().zip(&bv).enumerate() {
                assert_eq!(bits(x), bits(y), "{name}[{i}] differs");
            }
        }
    }

    fn solid_box(dims: [usize; 3]) -> Vec<bool> {
        let mut solid = vec![false; dims[0] * dims[1] * dims[2]];
        let idx = |x: usize, y: usize, z: usize| (z * dims[1] + y) * dims[0] + x;
        solid[idx(dims[0] / 2, dims[1] / 2, dims[2] / 2)] = true;
        solid
    }

    fn run_resume_bits<L, T>(dims: [usize; 3])
    where
        L: Lattice,
        T: Real,
    {
        let spec = GlobalSpec::<T> {
            dims,
            nu: 0.08,
            periodic: [true, true, L::D == 3],
            force: [
                T::r(1.0e-6),
                T::r(-2.0e-6),
                T::r(if L::D == 3 { 3.0e-7 } else { 0.0 }),
            ],
            ..Default::default()
        };
        let solid = solid_box(dims);
        let wall_u = vec![[T::zero(); 3]; solid.len()];
        let mut continuous: Solver<L, T, CpuScalar, LocalPeriodic> = Solver::new(
            &spec,
            &solid,
            &wall_u,
            [1, 1, 1],
            CpuScalar::default(),
            LocalPeriodic,
        );
        let mut resumed: Solver<L, T, CpuScalar, LocalPeriodic> = Solver::new(
            &spec,
            &solid,
            &wall_u,
            [1, 1, 1],
            CpuScalar::default(),
            LocalPeriodic,
        );
        continuous.run(100);
        resumed.run(50);
        let dir = tmp_ckpt("bit-resume");
        resumed.save(&dir).unwrap();
        let mut loaded: Solver<L, T, CpuScalar, LocalPeriodic> = Solver::load(
            &dir,
            &spec,
            &solid,
            &wall_u,
            [1, 1, 1],
            CpuScalar::default(),
            LocalPeriodic,
        )
        .unwrap();
        loaded.run(50);
        assert_solver_bits_eq(&continuous, &loaded);
        let _ = std::fs::remove_dir_all(dir);
    }

    #[test]
    fn checkpoint_resume_bit_exact_f32_f64_d2q9_d3q19() {
        run_resume_bits::<D2Q9, f32>([16, 12, 1]);
        run_resume_bits::<D2Q9, f64>([16, 12, 1]);
        run_resume_bits::<D3Q19, f32>([8, 7, 6]);
        run_resume_bits::<D3Q19, f64>([8, 7, 6]);
    }

    #[test]
    fn checkpoint_multi_part_roundtrip_bit_exact() {
        let spec = GlobalSpec::<f64> {
            dims: [24, 18, 1],
            nu: 0.04,
            periodic: [true, true, false],
            force: [1.0e-6, -2.0e-6, 0.0],
            ..Default::default()
        };
        let mut solid = vec![false; spec.dims[0] * spec.dims[1]];
        solid[(spec.dims[1] / 2) * spec.dims[0] + spec.dims[0] / 2] = true;
        let wall_u = vec![[0.0f64; 3]; solid.len()];
        let decomp = [2, 2, 1];
        let mut continuous: Solver<D2Q9, f64, CpuScalar, InProcess> = Solver::new(
            &spec,
            &solid,
            &wall_u,
            decomp,
            CpuScalar::default(),
            InProcess,
        );
        let mut resumed: Solver<D2Q9, f64, CpuScalar, InProcess> = Solver::new(
            &spec,
            &solid,
            &wall_u,
            decomp,
            CpuScalar::default(),
            InProcess,
        );
        let init = |x: usize, y: usize, _z: usize| {
            let kx = 2.0 * std::f64::consts::PI / spec.dims[0] as f64;
            let ky = 2.0 * std::f64::consts::PI / spec.dims[1] as f64;
            (
                1.0 + 0.002 * (kx * x as f64).cos() * (ky * y as f64).sin(),
                [
                    0.02 * (kx * x as f64).sin(),
                    -0.015 * (ky * y as f64).cos(),
                    0.0,
                ],
            )
        };
        continuous.init_with(init);
        resumed.init_with(init);
        continuous.run(20);
        resumed.run(7);
        let dir = tmp_ckpt("split-roundtrip");
        resumed.save(&dir).unwrap();
        let mut loaded: Solver<D2Q9, f64, CpuScalar, InProcess> = Solver::load(
            &dir,
            &spec,
            &solid,
            &wall_u,
            decomp,
            CpuScalar::default(),
            InProcess,
        )
        .unwrap();
        loaded.run(13);
        assert_solver_bits_eq(&continuous, &loaded);
        let _ = std::fs::remove_dir_all(dir);
    }

    #[test]
    fn checkpoint_layout_and_version_mismatch_are_precise() {
        let spec = GlobalSpec::<f32> {
            dims: [24, 16, 1],
            nu: 0.05,
            periodic: [true, true, false],
            force: [5.0e-7, 0.0, 0.0],
            ..Default::default()
        };
        let mut s: Solver<D2Q9, f32, CpuScalar, InProcess> =
            Solver::new(&spec, &[], &[], [2, 2, 1], CpuScalar::default(), InProcess);
        s.run(4);
        let dir = tmp_ckpt("layout-version");
        s.save(&dir).unwrap();

        let err = match Solver::<D2Q9, f32, CpuScalar, InProcess>::load(
            &dir,
            &spec,
            &[],
            &[],
            [4, 1, 1],
            CpuScalar::default(),
            InProcess,
        ) {
            Ok(_) => panic!("changed decomposition must reject checkpoint load"),
            Err(e) => e,
        };
        assert_eq!(err.code, "CKPT_DECOMP_MISMATCH");
        assert!(err.message.contains("decomp_hash differed"));

        let manifest_path = dir.join("manifest.json");
        let mut manifest: serde_json::Value =
            serde_json::from_str(&std::fs::read_to_string(&manifest_path).unwrap()).unwrap();
        manifest["format_version"] = serde_json::Value::from(1);
        std::fs::write(
            &manifest_path,
            serde_json::to_string_pretty(&manifest).unwrap(),
        )
        .unwrap();
        let err = match Solver::<D2Q9, f32, CpuScalar, InProcess>::load(
            &dir,
            &spec,
            &[],
            &[],
            [2, 2, 1],
            CpuScalar::default(),
            InProcess,
        ) {
            Ok(_) => panic!("old checkpoint version must reject checkpoint load"),
            Err(e) => e,
        };
        assert_eq!(err.code, "CKPT_VERSION_MISMATCH");
        assert!(err.message.contains("format_version 1"));
        let _ = std::fs::remove_dir_all(dir);
    }

    #[test]
    fn checkpoint_reports_truncated_payload_and_spec_mismatch() {
        let spec = GlobalSpec::<f32> {
            dims: [12, 10, 1],
            periodic: [true, true, false],
            force: [1.0e-6, 0.0, 0.0],
            ..Default::default()
        };
        let mut s: Solver<D2Q9, f32, CpuScalar, LocalPeriodic> = Solver::new(
            &spec,
            &[],
            &[],
            [1, 1, 1],
            CpuScalar::default(),
            LocalPeriodic,
        );
        s.run(3);
        let dir = tmp_ckpt("errors");
        s.save(&dir).unwrap();

        let mut changed = spec.clone();
        changed.nu = 0.09;
        let err = match Solver::<D2Q9, f32, CpuScalar, LocalPeriodic>::load(
            &dir,
            &changed,
            &[],
            &[],
            [1, 1, 1],
            CpuScalar::default(),
            LocalPeriodic,
        ) {
            Ok(_) => panic!("changed spec must reject checkpoint load"),
            Err(e) => e,
        };
        assert_eq!(err.code, "CKPT_SCENARIO_MISMATCH");
        assert!(err.message.contains("scenario_hash differed"));

        let rank = dir.join("rank_0000_part_0000.bin");
        let clean_rank = std::fs::read(&rank).unwrap();
        let mut corrupt = clean_rank.clone();
        let flip_at = corrupt.len() / 2;
        corrupt[flip_at] ^= 0x01;
        std::fs::write(&rank, corrupt).unwrap();
        let err = match Solver::<D2Q9, f32, CpuScalar, LocalPeriodic>::load(
            &dir,
            &spec,
            &[],
            &[],
            [1, 1, 1],
            CpuScalar::default(),
            LocalPeriodic,
        ) {
            Ok(_) => panic!("corrupted rank file must reject checkpoint load"),
            Err(e) => e,
        };
        assert_eq!(err.code, "CKPT_PAYLOAD_CORRUPT");

        let mut bytes = clean_rank;
        bytes.truncate(bytes.len() - 5);
        std::fs::write(&rank, bytes).unwrap();
        let err = match Solver::<D2Q9, f32, CpuScalar, LocalPeriodic>::load(
            &dir,
            &spec,
            &[],
            &[],
            [1, 1, 1],
            CpuScalar::default(),
            LocalPeriodic,
        ) {
            Ok(_) => panic!("truncated rank file must reject checkpoint load"),
            Err(e) => e,
        };
        assert!(matches!(
            err.code,
            "CKPT_TRUNCATED" | "CKPT_PAYLOAD_CORRUPT"
        ));
        let _ = std::fs::remove_dir_all(dir);
    }

    #[test]
    fn checkpoint_preserves_stale_stash_and_carried_solid_moments() {
        let mut faces = [FaceBC::<f32>::Closed; 6];
        faces[Face::XNeg.index()] = FaceBC::Velocity {
            u: [0.03, 0.0, 0.0],
        };
        faces[Face::XPos.index()] = FaceBC::Convective { u_conv: 0.03 };
        let dims = [18, 10, 1];
        let mut walls = WallSpec::<f32>::default();
        walls.is_wall[Face::YNeg.index()] = true;
        walls.is_wall[Face::YPos.index()] = true;
        let (solid, wall_u) = build_wall_rims(2, dims, &walls);
        let spec = GlobalSpec::<f32> {
            dims,
            nu: 0.06,
            periodic: [false, false, false],
            faces,
            force: [2.0e-6, 0.0, 0.0],
            ..Default::default()
        };
        let mut s: Solver<D2Q9, f32, CpuScalar, LocalPeriodic> = Solver::new(
            &spec,
            &solid,
            &wall_u,
            [1, 1, 1],
            CpuScalar::default(),
            LocalPeriodic,
        );
        s.run(8);
        {
            let f = s.fields_mut(0);
            let c = f.geom.cidx(0, 0, 0);
            f.rho[c] = 1.2345;
        }
        let saved_ftmp = s.fields(0).ftmp.clone();
        let dir = tmp_ckpt("stash-moments");
        s.save(&dir).unwrap();
        let loaded: Solver<D2Q9, f32, CpuScalar, LocalPeriodic> = Solver::load(
            &dir,
            &spec,
            &solid,
            &wall_u,
            [1, 1, 1],
            CpuScalar::default(),
            LocalPeriodic,
        )
        .unwrap();
        assert_eq!(
            loaded
                .fields(0)
                .ftmp
                .iter()
                .map(|&v| bits(v))
                .collect::<Vec<_>>(),
            saved_ftmp.iter().map(|&v| bits(v)).collect::<Vec<_>>()
        );
        assert_eq!(
            bits(loaded.fields(0).rho[loaded.fields(0).geom.cidx(0, 0, 0)]),
            bits(1.2345f32)
        );
        let _ = std::fs::remove_dir_all(dir);
    }

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

    use crate::lattice::{D3Q19, D3Q27};
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

    #[test]
    fn d3q27_accepts_all_open_face_kinds_before_build() {
        let mut faces = [FaceBC::<f64>::Closed; 6];
        faces[Face::XNeg.index()] = FaceBC::Velocity {
            u: [0.02, 0.0, 0.0],
        };
        faces[Face::XPos.index()] = FaceBC::Outflow;
        let spec = GlobalSpec::<f64> {
            dims: [6, 6, 6],
            periodic: [false, true, true],
            faces,
            ..Default::default()
        };
        assert!(spec.validate_lattice::<D3Q27>(&[]).is_ok());
        assert!(Solver::<D3Q27, f64, CpuScalar, LocalPeriodic>::try_new(
            &spec,
            &[],
            &[],
            [1, 1, 1],
            CpuScalar::default(),
            LocalPeriodic,
        )
        .is_ok());
    }

    #[test]
    fn d3q27_periodic_and_walled_specs_validate() {
        let periodic = GlobalSpec::<f64> {
            dims: [6, 6, 6],
            periodic: [true, true, true],
            ..Default::default()
        };
        assert!(periodic.validate_lattice::<D3Q27>(&[]).is_ok());

        let dims = [6, 6, 6];
        let walled = GlobalSpec::<f64> {
            dims,
            periodic: [false, false, false],
            ..Default::default()
        };
        assert!(walled
            .validate_lattice::<D3Q27>(&walled_box_solid(dims))
            .is_ok());
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
