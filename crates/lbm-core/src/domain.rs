//! Domain configuration: edge boundary conditions, collision operator
//! selection, and validated construction of a [`crate::sim::Simulation`].

use crate::real::Real;
use crate::sim::Simulation;

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

    /// Inward-pointing unit normal.
    pub(crate) fn n_in(self) -> (i32, i32) {
        match self {
            Edge::Left => (1, 0),
            Edge::Right => (-1, 0),
            Edge::Bottom => (0, 1),
            Edge::Top => (0, -1),
        }
    }

    pub(crate) fn index(self) -> usize {
        match self {
            Edge::Left => 0,
            Edge::Right => 1,
            Edge::Bottom => 2,
            Edge::Top => 3,
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
        for bc in [&e.left, &e.right, &e.bottom, &e.top] {
            let s = bc.max_speed();
            if s > MAX_SPEED {
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
}
