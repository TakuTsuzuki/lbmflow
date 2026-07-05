//! Domain configuration (V1 `lbm_core::domain` facade): edge boundary
//! conditions, collision operator selection, and validated construction of a
//! [`crate::compat::sim::Simulation`]. Semantics identical to V1.

use super::real::Real;
use super::sim::Simulation;

/// Boundary condition attached to one edge of the rectangular domain.
///
/// Wall-type edges (`BounceBack`, `MovingWall`) are realised as a one-cell
/// solid rim; the physical wall sits half-way between the rim cell centre and
/// the adjacent fluid cell centre. Open edges (`VelocityInlet`,
/// `PressureOutlet`, `Outflow`) act on fluid cells lying on the edge itself.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum EdgeBC<T: Real> {
    /// Wrap-around. Must be paired with `Periodic` on the opposite edge.
    Periodic,
    /// Stationary no-slip wall (half-way bounce-back).
    BounceBack,
    /// No-slip wall moving tangentially with velocity `u` (lattice units).
    MovingWall {
        /// Wall velocity `[ux, uy]`.
        u: [T; 2],
    },
    /// Zou–He velocity boundary: prescribes velocity `u` on the edge.
    VelocityInlet {
        /// Prescribed velocity `[ux, uy]`.
        u: [T; 2],
    },
    /// Zou–He pressure boundary: prescribes density `rho` (p = cs^2 rho).
    PressureOutlet {
        /// Prescribed density.
        rho: T,
    },
    /// Zero-gradient outflow (copies unknown populations from the interior).
    Outflow,
    /// Convective (radiation) outflow: `df/dt + Uc df/dn = 0`, discretised as
    /// `f(edge) = (f_prev(edge) + Uc f(interior)) / (1 + Uc)`. Far less
    /// pressure-reflective than `Outflow`; set `u_conv` to the expected mean
    /// outflow speed (lattice units).
    ConvectiveOutflow {
        /// Advection speed of the outgoing characteristics.
        u_conv: T,
    },
}

impl<T: Real> EdgeBC<T> {
    pub(crate) fn is_periodic(&self) -> bool {
        matches!(self, EdgeBC::Periodic)
    }
    pub(crate) fn is_wall(&self) -> bool {
        matches!(self, EdgeBC::BounceBack | EdgeBC::MovingWall { .. })
    }
    pub(crate) fn is_open(&self) -> bool {
        matches!(
            self,
            EdgeBC::VelocityInlet { .. }
                | EdgeBC::PressureOutlet { .. }
                | EdgeBC::Outflow
                | EdgeBC::ConvectiveOutflow { .. }
        )
    }
    fn max_speed(&self) -> f64 {
        match self {
            EdgeBC::MovingWall { u } | EdgeBC::VelocityInlet { u } => {
                (u[0].as_f64().powi(2) + u[1].as_f64().powi(2)).sqrt()
            }
            _ => 0.0,
        }
    }
}

/// Identifies one edge of the rectangular domain.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Edge {
    /// `x = 0` column.
    Left,
    /// `x = nx - 1` column.
    Right,
    /// `y = 0` row.
    Bottom,
    /// `y = ny - 1` row.
    Top,
}

impl Edge {
    /// All four edges.
    pub const ALL: [Edge; 4] = [Edge::Left, Edge::Right, Edge::Bottom, Edge::Top];

    /// The corresponding V2 core face.
    pub(crate) fn face(self) -> crate::lattice::Face {
        match self {
            Edge::Left => crate::lattice::Face::XNeg,
            Edge::Right => crate::lattice::Face::XPos,
            Edge::Bottom => crate::lattice::Face::YNeg,
            Edge::Top => crate::lattice::Face::YPos,
        }
    }
}

/// Boundary conditions for the four domain edges.
///
/// `left` is the `x = 0` column, `right` is `x = nx-1`, `bottom` is `y = 0`,
/// `top` is `y = ny-1`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Edges<T: Real> {
    /// Edge at `x = 0`.
    pub left: EdgeBC<T>,
    /// Edge at `x = nx - 1`.
    pub right: EdgeBC<T>,
    /// Edge at `y = 0`.
    pub bottom: EdgeBC<T>,
    /// Edge at `y = ny - 1`.
    pub top: EdgeBC<T>,
}

impl<T: Real> Default for Edges<T> {
    fn default() -> Self {
        Self {
            left: EdgeBC::Periodic,
            right: EdgeBC::Periodic,
            bottom: EdgeBC::Periodic,
            top: EdgeBC::Periodic,
        }
    }
}

/// Collision operator selection (accuracy/stability trade-off axis).
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Collision {
    /// Single-relaxation-time BGK: fastest, least stable, tau-dependent wall slip.
    Bgk,
    /// Two-relaxation-time: `magic = 3/16` places straight walls exactly
    /// half-way for parabolic flows. Recommended default.
    Trt {
        /// Magic parameter Λ = (1/ω+ − 1/2)(1/ω− − 1/2).
        magic: f64,
    },
}

impl Collision {
    /// The standard "magic" value 3/16 (exact half-way walls for Poiseuille).
    pub const MAGIC_STD: f64 = 3.0 / 16.0;
}

impl Default for Collision {
    fn default() -> Self {
        Collision::Trt {
            magic: Self::MAGIC_STD,
        }
    }
}

/// Maximum allowed prescribed speed (lattice units) before construction fails.
pub const MAX_SPEED: f64 = 0.3;

/// Errors detected when validating a [`SimConfig`].
#[derive(Clone, Debug, PartialEq)]
pub enum ConfigError {
    /// Kinematic viscosity must be positive (tau > 0.5).
    NonPositiveViscosity {
        /// Offending value.
        nu: f64,
    },
    /// Domain must be at least 3x3 cells.
    DomainTooSmall {
        /// Configured width.
        nx: usize,
        /// Configured height.
        ny: usize,
    },
    /// `Periodic` must appear on both edges of an axis or on neither.
    UnpairedPeriodic {
        /// `"x"` or `"y"`.
        axis: &'static str,
    },
    /// An open edge (inlet/outlet/outflow) may not share a corner with
    /// another open edge; perpendicular edges must be walls or periodic.
    AdjacentOpenEdges,
    /// Prescribed wall/inlet speed exceeds [`MAX_SPEED`] lattice units.
    VelocityTooHigh {
        /// Offending speed magnitude.
        speed: f64,
    },
    /// Prescribed outlet density must be positive.
    NonPositiveDensity {
        /// Offending value.
        rho: f64,
    },
    /// A model parameter is outside its valid range.
    InvalidParameter {
        /// Parameter name.
        what: &'static str,
        /// Offending value.
        value: f64,
    },
    /// A parameter that must be finite is NaN or infinite (e.g. a NaN inlet
    /// velocity that would otherwise slip past a `>`-comparison and corrupt
    /// the field a few steps in, or silently degrade a `MovingWall` to a
    /// static wall).
    NonFiniteParameter {
        /// Parameter name.
        what: &'static str,
    },
}

impl std::fmt::Display for ConfigError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ConfigError::NonPositiveViscosity { nu } => {
                write!(
                    f,
                    "kinematic viscosity must be > 0 (got {nu}); tau = 3*nu + 0.5 must exceed 0.5"
                )
            }
            ConfigError::DomainTooSmall { nx, ny } => {
                write!(f, "domain must be at least 3x3 cells (got {nx}x{ny})")
            }
            ConfigError::UnpairedPeriodic { axis } => {
                write!(
                    f,
                    "periodic BC on the {axis} axis must be set on both opposing edges"
                )
            }
            ConfigError::AdjacentOpenEdges => write!(
                f,
                "open edges (inlet/outlet/outflow) may not meet at a corner; \
                 perpendicular edges must be walls or periodic"
            ),
            ConfigError::VelocityTooHigh { speed } => write!(
                f,
                "prescribed speed {speed} exceeds the low-Mach limit {MAX_SPEED} (lattice units)"
            ),
            ConfigError::NonPositiveDensity { rho } => {
                write!(f, "prescribed density must be > 0 (got {rho})")
            }
            ConfigError::InvalidParameter { what, value } => {
                write!(f, "parameter {what} = {value} is outside its valid range")
            }
            ConfigError::NonFiniteParameter { what } => {
                write!(f, "parameter {what} must be finite (got NaN or infinity)")
            }
        }
    }
}

impl std::error::Error for ConfigError {}

/// Simulation configuration; call [`SimConfig::build`] to obtain a validated
/// [`Simulation`].
#[derive(Clone, Debug)]
pub struct SimConfig<T: Real> {
    /// Lattice width in cells.
    pub nx: usize,
    /// Lattice height in cells.
    pub ny: usize,
    /// Kinematic viscosity in lattice units; `tau = 3*nu + 0.5`.
    pub nu: f64,
    /// Collision operator.
    pub collision: Collision,
    /// Edge boundary conditions.
    pub edges: Edges<T>,
    /// Uniform body force (Guo forcing, 2nd order).
    pub force: [T; 2],
}

impl<T: Real> Default for SimConfig<T> {
    fn default() -> Self {
        Self {
            nx: 64,
            ny: 64,
            nu: 1.0 / 6.0, // tau = 1
            collision: Collision::default(),
            edges: Edges::default(),
            force: [T::zero(), T::zero()],
        }
    }
}

impl<T: Real> SimConfig<T> {
    /// Validate the configuration and build the simulation.
    pub fn build(self) -> Result<Simulation<T>, ConfigError> {
        self.validate()?;
        Ok(Simulation::from_config(self))
    }

    fn validate(&self) -> Result<(), ConfigError> {
        if !(self.nu > 0.0) {
            return Err(ConfigError::NonPositiveViscosity { nu: self.nu });
        }
        // TRT magic must be finite and positive (Λ ≤ 0 or NaN gives a
        // non-physical or degenerate ω−; the NaN case would otherwise slip
        // past a bare `>` test).
        if let Collision::Trt { magic } = self.collision {
            if !magic.is_finite() {
                return Err(ConfigError::NonFiniteParameter { what: "magic" });
            }
            if !(magic > 0.0) {
                return Err(ConfigError::InvalidParameter {
                    what: "magic",
                    value: magic,
                });
            }
        }
        // Uniform body force must be finite (a NaN force poisons every
        // fluid cell on the first collide).
        for (a, comp) in self.force.iter().enumerate() {
            if !comp.as_f64().is_finite() {
                return Err(ConfigError::NonFiniteParameter {
                    what: if a == 0 { "force[0]" } else { "force[1]" },
                });
            }
        }
        if self.nx < 3 || self.ny < 3 {
            return Err(ConfigError::DomainTooSmall {
                nx: self.nx,
                ny: self.ny,
            });
        }
        let e = &self.edges;
        if e.left.is_periodic() != e.right.is_periodic() {
            return Err(ConfigError::UnpairedPeriodic { axis: "x" });
        }
        if e.bottom.is_periodic() != e.top.is_periodic() {
            return Err(ConfigError::UnpairedPeriodic { axis: "y" });
        }
        let x_open = e.left.is_open() || e.right.is_open();
        let y_open = e.bottom.is_open() || e.top.is_open();
        if x_open && y_open {
            return Err(ConfigError::AdjacentOpenEdges);
        }
        // `normal_axis`: 0 for the x-facing edges (Left/Right), 1 for the
        // y-facing edges (Bottom/Top). A MovingWall may only slide
        // tangentially, so its velocity component along this axis must be 0.
        for (normal_axis, bc) in [(0, &e.left), (0, &e.right), (1, &e.bottom), (1, &e.top)] {
            // Prescribed velocities must be finite *before* the magnitude
            // test: a NaN component makes `max_speed()` NaN, and a bare
            // `NaN > MAX_SPEED` is false — the exact silent pass-through that
            // corrupts an inlet field or degrades a MovingWall to a static
            // wall a few steps in.
            if let EdgeBC::MovingWall { u } | EdgeBC::VelocityInlet { u } = bc {
                if !u[0].as_f64().is_finite() || !u[1].as_f64().is_finite() {
                    return Err(ConfigError::NonFiniteParameter {
                        what: "edge velocity",
                    });
                }
            }
            // A-6: a wall-normal MovingWall component injects/removes mass at
            // the half-way bounce-back without diverging — a silent leak
            // (E7: −56% mass over 500 steps). Reject it; only tangential
            // motion is physically representable by BB momentum injection.
            if let EdgeBC::MovingWall { u } = bc {
                if u[normal_axis].as_f64() != 0.0 {
                    return Err(ConfigError::InvalidParameter {
                        what: "MovingWall normal velocity component (walls may only slide tangentially)",
                        value: u[normal_axis].as_f64(),
                    });
                }
            }
            let s = bc.max_speed();
            // NaN-safe: `!(s <= MAX_SPEED)` rejects NaN as well as too-fast.
            if !(s <= MAX_SPEED) {
                return Err(ConfigError::VelocityTooHigh { speed: s });
            }
            if let EdgeBC::PressureOutlet { rho } = bc {
                if !(rho.as_f64() > 0.0) {
                    return Err(ConfigError::NonPositiveDensity { rho: rho.as_f64() });
                }
            }
            if let EdgeBC::ConvectiveOutflow { u_conv } = bc {
                let v = u_conv.as_f64();
                if !(v > 0.0 && v <= 1.0) {
                    return Err(ConfigError::InvalidParameter {
                        what: "u_conv",
                        value: v,
                    });
                }
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_config_builds() {
        let sim = SimConfig::<f64>::default().build().unwrap();
        assert_eq!(sim.nx(), 64);
        assert_eq!(sim.ny(), 64);
    }

    #[test]
    fn rejects_bad_configs() {
        let bad_nu = SimConfig::<f64> {
            nu: 0.0,
            ..Default::default()
        };
        assert!(matches!(
            bad_nu.build().unwrap_err(),
            ConfigError::NonPositiveViscosity { .. }
        ));

        let unpaired = SimConfig::<f64> {
            edges: Edges {
                left: EdgeBC::Periodic,
                right: EdgeBC::BounceBack,
                bottom: EdgeBC::BounceBack,
                top: EdgeBC::BounceBack,
            },
            ..Default::default()
        };
        assert!(matches!(
            unpaired.build().unwrap_err(),
            ConfigError::UnpairedPeriodic { axis: "x" }
        ));

        let open_corner = SimConfig::<f64> {
            edges: Edges {
                left: EdgeBC::VelocityInlet { u: [0.05, 0.0] },
                right: EdgeBC::Outflow,
                bottom: EdgeBC::Outflow,
                top: EdgeBC::BounceBack,
            },
            ..Default::default()
        };
        assert!(matches!(
            open_corner.build().unwrap_err(),
            ConfigError::AdjacentOpenEdges
        ));

        let too_fast = SimConfig::<f64> {
            edges: Edges {
                left: EdgeBC::Periodic,
                right: EdgeBC::Periodic,
                bottom: EdgeBC::BounceBack,
                top: EdgeBC::MovingWall { u: [0.5, 0.0] },
            },
            ..Default::default()
        };
        assert!(matches!(
            too_fast.build().unwrap_err(),
            ConfigError::VelocityTooHigh { .. }
        ));
    }

    /// A-2 (E6): NaN/inf on any prescribed velocity, the body force, or the
    /// TRT magic parameter must be rejected at `build()` rather than slipping
    /// past a bare `>`-comparison (NaN inlet → field NaN in ~3 steps; NaN
    /// MovingWall → silent static wall).
    #[test]
    fn rejects_non_finite_parameters() {
        let nan = f64::NAN;
        let inf = f64::INFINITY;

        // NaN inlet velocity (would corrupt the field, not caught by `>`).
        let nan_inlet = SimConfig::<f64> {
            edges: Edges {
                left: EdgeBC::VelocityInlet { u: [nan, 0.0] },
                right: EdgeBC::Outflow,
                bottom: EdgeBC::BounceBack,
                top: EdgeBC::BounceBack,
            },
            ..Default::default()
        };
        assert!(matches!(
            nan_inlet.build().unwrap_err(),
            ConfigError::NonFiniteParameter { .. }
        ));

        // NaN MovingWall velocity (would degrade to a static wall, no error).
        let nan_wall = SimConfig::<f64> {
            edges: Edges {
                left: EdgeBC::Periodic,
                right: EdgeBC::Periodic,
                bottom: EdgeBC::BounceBack,
                top: EdgeBC::MovingWall { u: [nan, 0.0] },
            },
            ..Default::default()
        };
        assert!(matches!(
            nan_wall.build().unwrap_err(),
            ConfigError::NonFiniteParameter { .. }
        ));

        // inf inlet velocity.
        let inf_inlet = SimConfig::<f64> {
            edges: Edges {
                left: EdgeBC::VelocityInlet { u: [0.0, inf] },
                right: EdgeBC::Outflow,
                bottom: EdgeBC::BounceBack,
                top: EdgeBC::BounceBack,
            },
            ..Default::default()
        };
        assert!(matches!(
            inf_inlet.build().unwrap_err(),
            ConfigError::NonFiniteParameter { .. }
        ));

        // NaN body force component.
        let nan_force = SimConfig::<f64> {
            force: [0.0, nan],
            ..Default::default()
        };
        assert!(matches!(
            nan_force.build().unwrap_err(),
            ConfigError::NonFiniteParameter { .. }
        ));

        // NaN TRT magic.
        let nan_magic = SimConfig::<f64> {
            collision: Collision::Trt { magic: nan },
            ..Default::default()
        };
        assert!(matches!(
            nan_magic.build().unwrap_err(),
            ConfigError::NonFiniteParameter { .. }
        ));

        // Non-positive TRT magic (Λ ≤ 0).
        let bad_magic = SimConfig::<f64> {
            collision: Collision::Trt { magic: 0.0 },
            ..Default::default()
        };
        assert!(matches!(
            bad_magic.build().unwrap_err(),
            ConfigError::InvalidParameter { what: "magic", .. }
        ));

        // Sanity: the finite defaults still build.
        assert!(SimConfig::<f64>::default().build().is_ok());
    }

    /// A-6 (E7): a MovingWall velocity component normal to its own edge is a
    /// silent mass leak (bounce-back injects/removes mass without diverging).
    /// Tangential motion is fine; a normal component must be rejected.
    #[test]
    fn rejects_moving_wall_normal_component() {
        // Top wall (y-facing): tangential = x-component only.
        let tangential = SimConfig::<f64> {
            edges: Edges {
                left: EdgeBC::BounceBack,
                right: EdgeBC::BounceBack,
                bottom: EdgeBC::BounceBack,
                top: EdgeBC::MovingWall { u: [0.05, 0.0] },
            },
            ..Default::default()
        };
        assert!(tangential.build().is_ok(), "tangential lid must be allowed");

        // Top wall with a normal (y) component — the E7 leak.
        let normal = SimConfig::<f64> {
            edges: Edges {
                left: EdgeBC::BounceBack,
                right: EdgeBC::BounceBack,
                bottom: EdgeBC::BounceBack,
                top: EdgeBC::MovingWall { u: [0.0, -0.05] },
            },
            ..Default::default()
        };
        assert!(matches!(
            normal.build().unwrap_err(),
            ConfigError::InvalidParameter { .. }
        ));

        // Left wall (x-facing): normal = x-component.
        let left_normal = SimConfig::<f64> {
            edges: Edges {
                left: EdgeBC::MovingWall { u: [0.05, 0.0] },
                right: EdgeBC::BounceBack,
                bottom: EdgeBC::BounceBack,
                top: EdgeBC::BounceBack,
            },
            ..Default::default()
        };
        assert!(matches!(
            left_normal.build().unwrap_err(),
            ConfigError::InvalidParameter { .. }
        ));
        // Left wall sliding tangentially (y) is fine.
        let left_tangential = SimConfig::<f64> {
            edges: Edges {
                left: EdgeBC::MovingWall { u: [0.0, 0.05] },
                right: EdgeBC::BounceBack,
                bottom: EdgeBC::BounceBack,
                top: EdgeBC::BounceBack,
            },
            ..Default::default()
        };
        assert!(left_tangential.build().is_ok());
    }
}
