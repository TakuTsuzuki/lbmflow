//! Marker-based direct-forcing immersed-boundary model for rigid rotation.
//!
//! The implementation follows the Uhlmann direct-forcing structure: interpolate
//! the fluid velocity to Lagrangian markers, compute the force needed to move
//! each marker to its rigid-body target velocity, then spread the force back to
//! the Eulerian Guo force field. Repeated sweeps are the Wang-style
//! multi-direct-forcing correction used when a single interpolation/spreading
//! sweep leaves too much marker slip.

use crate::real::Real;

/// One Lagrangian IBM marker in lattice coordinates.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IbmMarker {
    /// Marker position in cell-center lattice coordinates.
    pub position: [f64; 3],
    /// Quadrature weight in lattice cell units. For a 2D closed curve this is
    /// the marker arc length; for 3D surfaces it is the marker area.
    pub weight: f64,
}

/// Rigid body described by Lagrangian markers and an angular velocity.
#[derive(Clone, Debug)]
pub struct RotatingBody {
    center: [f64; 3],
    omega: [f64; 3],
    markers: Vec<IbmMarker>,
}

impl RotatingBody {
    /// Build from an explicit marker set.
    pub fn from_markers(center: [f64; 3], omega: [f64; 3], markers: Vec<IbmMarker>) -> Self {
        assert!(
            !markers.is_empty(),
            "IBM body must contain at least one marker"
        );
        for (i, m) in markers.iter().enumerate() {
            assert!(
                m.position.iter().all(|v| v.is_finite()) && m.weight.is_finite() && m.weight > 0.0,
                "IBM marker {i} must have finite position and positive finite weight"
            );
        }
        assert!(
            center.iter().all(|v| v.is_finite()) && omega.iter().all(|v| v.is_finite()),
            "IBM center and omega must be finite"
        );
        Self {
            center,
            omega,
            markers,
        }
    }

    /// Uniform marker set on a 2D rotating circular boundary.
    pub fn circle_2d(center: [f64; 2], radius: f64, omega_z: f64, n_markers: usize) -> Self {
        assert!(
            radius > 0.0 && radius.is_finite(),
            "radius must be finite and positive"
        );
        assert!(n_markers >= 8, "circle IBM needs at least 8 markers");
        let ds = std::f64::consts::TAU * radius / n_markers as f64;
        let mut markers = Vec::with_capacity(n_markers);
        for i in 0..n_markers {
            let th = std::f64::consts::TAU * i as f64 / n_markers as f64;
            markers.push(IbmMarker {
                position: [
                    center[0] + radius * th.cos(),
                    center[1] + radius * th.sin(),
                    0.0,
                ],
                weight: ds,
            });
        }
        Self::from_markers([center[0], center[1], 0.0], [0.0, 0.0, omega_z], markers)
    }

    pub fn center(&self) -> [f64; 3] {
        self.center
    }

    pub fn omega(&self) -> [f64; 3] {
        self.omega
    }

    pub fn markers(&self) -> &[IbmMarker] {
        &self.markers
    }

    /// Rigid target velocity `U = Omega x r`.
    pub fn target_velocity(&self, p: [f64; 3]) -> [f64; 3] {
        let r = [
            p[0] - self.center[0],
            p[1] - self.center[1],
            p[2] - self.center[2],
        ];
        [
            self.omega[1] * r[2] - self.omega[2] * r[1],
            self.omega[2] * r[0] - self.omega[0] * r[2],
            self.omega[0] * r[1] - self.omega[1] * r[0],
        ]
    }
}

/// Controls for one direct-forcing IBM update.
#[derive(Clone, Copy, Debug)]
pub struct DirectForcingConfig {
    /// Maximum number of direct-forcing sweeps. Values above one enable
    /// multi-direct-forcing.
    pub max_iterations: usize,
    /// Stop iterating when `max(|u_marker - U|) / max(|U|)` is at or below this
    /// threshold. The denominator is floored at 1 for stationary-marker cases.
    pub slip_tolerance: f64,
    /// Kernel support radius in cells. `1` is the 2-point linear kernel; `2`
    /// is the 3-point quadratic B-spline kernel.
    pub kernel_radius: usize,
    /// Under-relaxation for the force increment. `1` is the direct-forcing
    /// value; smaller values are useful for stiff marker sets near walls.
    pub relaxation: f64,
}

impl Default for DirectForcingConfig {
    fn default() -> Self {
        Self {
            max_iterations: 3,
            slip_tolerance: 1.0e-2,
            kernel_radius: 1,
            relaxation: 1.0,
        }
    }
}

/// Diagnostics returned by an IBM force update.
#[derive(Clone, Copy, Debug, Default)]
pub struct IbmDiagnostics {
    /// Number of direct-forcing sweeps actually performed.
    pub iterations: usize,
    /// Maximum marker slip magnitude after the final sweep.
    pub slip_max: f64,
    /// RMS marker slip magnitude after the final sweep.
    pub slip_rms: f64,
    /// Maximum marker slip divided by the maximum target speed.
    pub slip_max_rel: f64,
    /// RMS marker slip divided by the maximum target speed.
    pub slip_rms_rel: f64,
    /// Reaction torque on the body, `sum r x (-F_fluid)`.
    pub torque: [f64; 3],
    /// Net force applied to the fluid by the IBM spread.
    pub fluid_force: [f64; 3],
    /// Net force represented by the marker quadrature before spreading.
    pub marker_force: [f64; 3],
    /// Relative force-spreading conservation error.
    pub momentum_error_rel: f64,
}

#[derive(Clone, Copy, Debug)]
pub(crate) struct StencilPoint {
    pub x: usize,
    pub y: usize,
    pub z: usize,
    pub w: f64,
}

pub(crate) fn marker_stencil(
    p: [f64; 3],
    dims: [usize; 3],
    d: usize,
    radius: usize,
) -> Vec<StencilPoint> {
    assert!(
        radius == 1 || radius == 2,
        "IBM kernel radius must be 1 or 2"
    );
    let mut axes = [Vec::<(usize, f64)>::new(), Vec::new(), Vec::new()];
    for a in 0..3 {
        if a >= d {
            axes[a].push((0, 1.0));
            continue;
        }
        let lo = (p[a].floor() as isize) - radius as isize + 1;
        let hi = (p[a].floor() as isize) + radius as isize;
        for i in lo..=hi {
            if i < 0 || i >= dims[a] as isize {
                continue;
            }
            let r = (p[a] - i as f64).abs();
            let w = if radius == 1 {
                linear_kernel(r)
            } else {
                quadratic_kernel(r)
            };
            if w != 0.0 {
                axes[a].push((i as usize, w));
            }
        }
    }
    let mut out = Vec::new();
    for &(x, wx) in &axes[0] {
        for &(y, wy) in &axes[1] {
            for &(z, wz) in &axes[2] {
                out.push(StencilPoint {
                    x,
                    y,
                    z,
                    w: wx * wy * wz,
                });
            }
        }
    }
    out
}

fn linear_kernel(r: f64) -> f64 {
    (1.0 - r).max(0.0)
}

fn quadratic_kernel(r: f64) -> f64 {
    if r < 0.5 {
        0.75 - r * r
    } else if r < 1.5 {
        0.5 * (1.5 - r) * (1.5 - r)
    } else {
        0.0
    }
}

pub(crate) fn add3(a: &mut [f64; 3], b: [f64; 3]) {
    for i in 0..3 {
        a[i] += b[i];
    }
}

pub(crate) fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

pub(crate) fn norm(v: [f64; 3], d: usize) -> f64 {
    v.iter().take(d).map(|x| x * x).sum::<f64>().sqrt()
}

pub(crate) fn to_real3<T: Real>(v: [f64; 3]) -> [T; 3] {
    [T::r(v[0]), T::r(v[1]), T::r(v[2])]
}
