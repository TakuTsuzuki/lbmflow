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
pub mod credibility;
pub mod damage;
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
pub mod microcarrier;
pub mod oxygen;
pub mod params;
pub mod particles;
pub mod pbm;
pub mod phase_field;
pub mod qoi;
pub mod reaction;
pub mod real;
pub mod rotating_ibm;
pub mod scalar;
pub mod scaleup;
pub mod solver;
pub mod stress;
pub mod subdomain;
pub mod surface_tension;
pub mod uq;
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
    pub use crate::cells::{
        CellCheckpointSection, CellFieldSample, CellTracer, CellTracerPopulation,
    };
    pub use crate::credibility::{
        CalibrationDataset, CredibilityError, DatasetRegistry, HoldoutDataset,
    };
    pub use crate::damage::{
        exposure_distribution, DamageIncrement, DamageModelError, DamageThreshold,
        ExposureDistribution, ShearDamageModel,
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
    pub use crate::microcarrier::{
        drag_force_on_particle, scatter_drag_reaction_forces, terminal_velocity_stokes,
        validate_mass_loading, validate_reynolds, MicrocarrierError, MicrocarrierPopulation,
        ParticleDragForce, SuspensionMetrics, TwoWayScatterReport, TWO_WAY_MASS_LOADING_MAX,
    };
    pub use crate::oxygen::{
        apply_interfacial_flux_sources, clip_negative_concentrations, henry_equilibrium,
        interfacial_area_density, oxygen_source_step, OxygenDiagnostics, OxygenError,
        OxygenFluxLedger, OxygenState, OXYGEN_SCALAR_NAME,
    };
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
    pub use crate::qoi::{
        compartment_cv, dynamic_gassing_kla_fit, dynamic_gassing_window_default,
        mixing_time_from_cv, power_qois, scalar_cv, CompartmentCv, KlaDynamicFitOutcome,
        KlaDynamicFitResult, KlaFitMethod, KlaFitWindow, MixingQoiSection, MixingTimeResult,
        PowerQoiInput, PowerQoiResult, PowerQoiSection, QoiAccumulatorSnapshot, QoiBundle,
        QoiCheckpointState, QoiInterval, QoiPercentiles, QoiProvenance, QoiScalar,
        QoiValidationStatus, ShearQoiSection, SkippedQoi, ValidationTier,
    };
    pub use crate::reaction::{
        apply_oxygen_reaction_source, oxygen_uptake_rate, OurLedger, OurModel, ReactionError,
    };
    pub use crate::real::Real;
    pub use crate::rotating_ibm::{DirectForcingConfig, IbmDiagnostics, IbmMarker, RotatingBody};
    pub use crate::scalar::scalar_equilibrium;
    pub use crate::scaleup::{
        evaluate_operating_window, ConstraintConflict, ConstraintSet, ConstraintTightness,
        OperatingPoint, ScaleUpEvaluation, ScaleUpMode, ScaleUpQois,
    };
    pub use crate::solver::{
        build_wall_rims, partition, CheckpointError, Diverged, GlobalSpec, Solver,
        SolverFeatureError, SpecError, WallSpec,
    };
    pub use crate::stress::{
        compute_stress_field, percentile_summary, wall_shear_proxy, PercentileSummary, StressCell,
        WallShearProxy,
    };
    pub use crate::subdomain::Subdomain;
    pub use crate::surface_tension::{
        SurfaceTensionDiagnostics, SurfaceTensionError, SurfaceTensionParams,
    };
    pub use crate::uq::{
        combine_interval, one_factor_local_sensitivity, UqComponent, UqComponentKind,
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
