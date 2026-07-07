//! # lbm-core
//!
//! Lattice Boltzmann core (the V2 architecture, sole engine since the V1
//! retirement on 2026-07-05): the dimension / lattice / precision / backend /
//! decomposition axes are orthogonal (docs/ARCHITECTURE_V2.md). The physics
//! kernels are written once, generically over a [`lattice::Lattice`] and a
//! [`real::Real`], and specialise at compile time.
//!
//! Layer map (docs/ARCHITECTURE_V2.md §1):
//!
//! - [`lattice`] — compile-time velocity sets (D2Q9, D3Q19, D3Q27) with derived
//!   tables (TRT pairs, per-face unknown sets).
//! - [`fields`] — q-major SoA deviation storage over halo-padded local boxes.
//! - `kernels` (private) — the physics (collide/stream/moments/BCs), written
//!   once, generic over lattice and precision; V1-faithful arithmetic.
//! - [`backend`] — the compute-target trait and the `CpuScalar` reference.
//! - [`subdomain`] / [`halo`] — decomposition and halo exchange.
//! - [`solver`] — the orchestrator (V1 step sequence over parts).
//! - [`compat`] — the retired V1 engine's public API as a supported facade
//!   (scenario/CLI/wasm run through it); V1↔V2 equivalence was proven by
//!   `tests/v1_match.rs`, frozen and removed with V1 (see branch history).

pub mod backend;
pub mod backend_simd;
pub mod bench_regression;
pub mod bouzidi;
pub mod bubble_forces;
pub mod bubbles;
pub mod cells;
mod collision;
pub mod compat;
#[cfg(feature = "mpi")]
pub mod dist;
pub mod divergence;
pub mod fields;
pub mod free_surface;
pub mod geometry;
#[cfg(feature = "gpu")]
pub mod gpu;
pub mod halo;
pub mod hybrid_gas;
mod kernels;
pub mod kla;
pub mod lattice;
pub mod les;
pub mod materials;
pub mod params;
pub mod particles;
pub mod pbm;
pub mod phase_field;
pub mod qoi;
pub mod real;
pub mod rotating_ibm;
pub mod solver;
pub mod subdomain;
pub mod surface_tension;
#[cfg(feature = "geometry-import")]
pub mod voxel_import;
pub mod wall_model;
pub mod wetting;

/// Convenient glob import for the V2 API.
pub mod prelude {
    pub use crate::backend::{Backend, CellRange, CpuScalar, HostMoments, PARALLEL_MIN_CELLS};
    pub use crate::backend_simd::CpuSimd;
    pub use crate::bouzidi::{BouzidiLink, BouzidiLinks};
    pub use crate::bubble_forces::{
        added_mass_force, bubble_reynolds, buoyancy_force, lift_placeholder_force, rk4_substep,
        schiller_naumann_drag_force, turbulent_dispersion_placeholder_force,
        wall_lubrication_placeholder_force, BubbleForceContext, ClosureValidity,
        SCHILLER_NAUMANN_RE_MAX,
    };
    pub use crate::bubbles::{
        bubble_volume_from_diameter, validate_bubble_diameter, Bubble, BubbleError, BubbleSet,
        MomentumCouplingLedger, SpargerBubbleInjector, POINT_BUBBLE_ALPHA_G_MAX,
    };
    pub use crate::divergence::{DivergenceError, PhaseDiag};
    pub use crate::fields::{DistributionKind, LocalGeom, ScalarDistribution, SoaFields};
    pub use crate::free_surface::{DegassingLedger, FreeSurfaceError, TopBoundaryMode};
    pub use crate::geometry::{
        build_stirred_tank_geometry, generate_impeller_marker_set, BaffleTemplate, ContactAngleMap,
        GeometryError, GridSpec, ImpellerKind as GeometryImpellerKind, ImpellerMarkerSet,
        ImpellerTemplate, PipeAxis, SpargerTemplate, StirredTankGeometry, TankBottom, TankSpec,
        WallContactAngle, SPARGER_ORIFICE_MIN_CELLS, STIRRED_TANK_MIN_CELLS,
    };
    pub use crate::halo::{HaloExchange, InProcess, LocalPeriodic};
    pub use crate::hybrid_gas::{
        hybrid_gas_bookkeeping, reject_hybrid_evidence_tier, HybridGasReport,
    };
    pub use crate::kla::{
        compute_kla_from_alpha_d32, compute_kla_from_pbm_bins, interfacial_area_from_alpha_d32,
        oxygen_transfer_rate_mol_m3_s, KlModel, KlaProvenance, KlaReport,
    };
    pub use crate::lattice::{Face, Lattice, D2Q9, D3Q19, D3Q27};
    pub use crate::les::{WaleLes, WaleLesDiagnostics, WALE_CW};
    pub use crate::materials::{MaterialFields, MaterialSample};
    pub use crate::params::{
        CollisionKind, FaceBC, FacePatch, MaterialModel, MaterialParamError,
        PhaseFieldMixtureParams, Reduction, SourceKind, SourceRegion, StepParams,
        ViscosityInterpolation, VolumeSource,
    };
    pub use crate::pbm::{
        BreakupKernel, CoalescenceKernel, ConstantBreakup, ConstantCoalescence, DisabledKernel,
        FutureKernelHook, PbmBins, PbmLocalState, PbmValidity, DEFAULT_PBM_BIN_COUNT,
    };
    pub use crate::phase_field::{
        ClippingPolicy, PhaseFieldDiagnostics, PhaseFieldError, PhaseFieldParams,
    };
    pub use crate::qoi::{QoiAccumulatorSnapshot, QoiCheckpointState};
    pub use crate::real::Real;
    pub use crate::rotating_ibm::{DirectForcingConfig, IbmDiagnostics, IbmMarker, RotatingBody};
    pub use crate::solver::{
        build_wall_rims, partition, CheckpointError, Diverged, GlobalSpec, Solver,
        SolverFeatureError, SpecError, WallSpec,
    };
    pub use crate::subdomain::Subdomain;
    pub use crate::surface_tension::{
        SurfaceTensionDiagnostics, SurfaceTensionError, SurfaceTensionParams,
    };
    pub use crate::wall_model::{WallCellMetric, WallMetricSource};
    pub use crate::wetting::{ContactAngleParams, WettingError};

    #[cfg(feature = "gpu")]
    #[allow(deprecated)]
    pub use crate::gpu::{GpuContext, GpuSolver, WgpuBackend};

    #[cfg(feature = "mpi")]
    pub use crate::dist::{
        read_parallel_field, FieldDtype, GlobalIndex, MpiExchange, MpiSolver, ParallelFieldError,
        ParallelFieldManifest, RankFieldSlab, ScalarName,
    };
}
