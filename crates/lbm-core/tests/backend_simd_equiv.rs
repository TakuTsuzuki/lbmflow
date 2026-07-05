//! Backend equivalence gate: `CpuScalar` vs `CpuSimd` (the fused V1
//! `step_band` port) must produce the same trajectories on identical
//! scenarios, in the v1_match style: **f64 max |Δ| ≤ 1e-11, f32 ≤ 1e-6**
//! per cell over fields (rho / u / fluid-cell populations).
//!
//! Diagnostics are gated per precision the way v1_match does: `f64`
//! compares total mass / momentum / probed force at the same absolute
//! 1e-11 (reassociation noise stays orders below it even after the f64
//! accumulation over the grid). For `f32`, mass and momentum are
//! *extensive* f64 sums of per-cell f32 state, so backend last-ulp drift
//! accumulates linearly with the fluid-cell count N — the dimensionally
//! consistent line is **|Δ| ≤ 1e-6 · N** (observed: ~1e-9 · N, i.e. three
//! orders inside; v1_match's f32 case gates fields only, sidestepping the
//! question). The probed force is a surface sum over O(10²) links of O(0.1)
//! physical populations; its f32 line is **1e-5 absolute** (observed
//! ≤ ~1.3e-6; any real link-accounting bug shifts it by ≥ 1e-3).
//!
//! The fused backend replicates the `kernels.rs` expressions for streaming /
//! bounce-back / moments / BCs and replays the probe fold in `CpuScalar`'s
//! order; its collision uses V1 `collide_span`'s pair-shared form (an exact
//! algebraic regrouping of `kernels::collide_row` — see the `backend_simd`
//! module docs for the measured 1T cost of the literal DAG), so the two
//! backends drift by TRT-pair last-ulp reassociation only: observed
//! ~1e-13 (f64) / ~1e-7 (f32) over the suite, far inside the lines above.
//! The measured worst deltas are printed for the record.
//!
//! Coverage (mission list, each in f64 and f32):
//!   1. 2D TGV, fully periodic, TRT       (fused kernel + periodic halo)
//!   2. 2D lid-driven cavity              (still + moving-wall bounce-back)
//!   3. 2D channel, inlet profile → Outflow
//!   4. 2D cylinder + momentum-exchange probe (probe force every step)
//!   5. 2D per-cell + uniform Guo force, BGK
//!   6. 2D channel → ConvectiveOutflow    (stale-slot capture path)
//!   7. 3D TGV (D3Q19, fully periodic)
//!   8. 3D duct: inlet profile → Outflow, four wall faces
//! plus the decomposition integration tests: InProcess 2×2 (and 3D 2×2×1)
//! over `CpuSimd` against the monolithic `CpuScalar` reference.

use lbm_core::lattice::{D2Q9, D3Q19};
use lbm_core::prelude::*;
use std::f64::consts::PI;

/// v1_match acceptance lines (per-cell fields).
fn tol<T: Real>() -> f64 {
    if std::mem::size_of::<T>() == 4 {
        1e-6
    } else {
        1e-11
    }
}

/// Extensive-diagnostic line (total mass / momentum): per-cell line scaled
/// by the fluid-cell count for f32 (see module docs); absolute for f64.
fn tol_extensive<T: Real>(n_fluid: usize) -> f64 {
    if std::mem::size_of::<T>() == 4 {
        1e-6 * n_fluid as f64
    } else {
        1e-11
    }
}

/// Probed-force line (surface sum; see module docs).
fn tol_probe<T: Real>() -> f64 {
    if std::mem::size_of::<T>() == 4 {
        1e-5
    } else {
        1e-11
    }
}

fn max_abs_diff<T: Real>(a: &[T], b: &[T]) -> f64 {
    assert_eq!(a.len(), b.len());
    a.iter()
        .zip(b)
        .map(|(x, y)| (x.as_f64() - y.as_f64()).abs())
        .fold(0.0, f64::max)
}

/// Compare fields + diagnostics of two solvers over identical scenarios.
/// Population planes are compared over fluid cells only (solid cells hold
/// dead ping-pong junk whose history legitimately differs between the
/// backends). Returns the worst difference seen.
fn compare_state<L, T, BA, BB, HA, HB>(
    a: &Solver<L, T, BA, HA>,
    b: &Solver<L, T, BB, HB>,
    what: &str,
) -> f64
where
    L: Lattice,
    T: Real,
    BA: Backend<L, T, Fields = SoaFields<T>>,
    BB: Backend<L, T, Fields = SoaFields<T>>,
    HA: HaloExchange<T>,
    HB: HaloExchange<T>,
{
    let lim = tol::<T>();
    let mut worst = 0.0f64;
    let mut chk = |name: &str, d: f64, lim: f64| {
        assert!(d <= lim, "{what}: {name} max|Δ| = {d:e} > {lim:e}");
        if d > worst {
            worst = d;
        }
    };
    chk("rho", max_abs_diff(&a.gather_rho(), &b.gather_rho()), lim);
    chk("ux", max_abs_diff(&a.gather_ux(), &b.gather_ux()), lim);
    chk("uy", max_abs_diff(&a.gather_uy(), &b.gather_uy()), lim);
    if L::D == 3 {
        chk("uz", max_abs_diff(&a.gather_uz(), &b.gather_uz()), lim);
    }
    // Populations over fluid cells.
    let dims = a.dims();
    let fluid: Vec<bool> = (0..dims[0] * dims[1] * dims[2])
        .map(|i| {
            let x = i % dims[0];
            let y = (i / dims[0]) % dims[1];
            let z = i / (dims[0] * dims[1]);
            !a.is_solid(x, y, z)
        })
        .collect();
    for q in 0..L::Q {
        let (fa, fb) = (a.gather_f(q), b.gather_f(q));
        let d = fa
            .iter()
            .zip(&fb)
            .zip(&fluid)
            .filter(|&(_, &fl)| fl)
            .map(|((x, y), _)| (x.as_f64() - y.as_f64()).abs())
            .fold(0.0, f64::max);
        chk(&format!("f[{q}]"), d, lim);
    }
    // f64-accumulated diagnostics (per-precision lines; see module docs).
    let n_fluid = fluid.iter().filter(|&&b| b).count();
    let ext = tol_extensive::<T>(n_fluid);
    chk(
        "total_mass",
        (a.total_mass().as_f64() - b.total_mass().as_f64()).abs(),
        ext,
    );
    let (pa, pb) = (a.total_momentum(), b.total_momentum());
    for c in 0..L::D {
        chk(
            &format!("momentum[{c}]"),
            (pa[c].as_f64() - pb[c].as_f64()).abs(),
            ext,
        );
    }
    let (fa, fb) = (a.probed_force(), b.probed_force());
    for c in 0..L::D {
        chk(
            &format!("probed_force[{c}]"),
            (fa[c].as_f64() - fb[c].as_f64()).abs(),
            tol_probe::<T>(),
        );
    }
    worst
}

/// A CpuScalar/CpuSimd pair over the identical monolithic scenario.
struct Pair<L: Lattice, T: Real> {
    a: Solver<L, T, CpuScalar, LocalPeriodic>,
    b: Solver<L, T, CpuSimd, LocalPeriodic>,
}

impl<L: Lattice, T: Real> Pair<L, T> {
    fn new(spec: &GlobalSpec<T>, walls: &WallSpec<T>) -> Self {
        let (solid, wall_u) = build_wall_rims(L::D, spec.dims, walls);
        Self {
            a: Solver::new(
                spec,
                &solid,
                &wall_u,
                [1, 1, 1],
                CpuScalar::default(),
                LocalPeriodic,
            ),
            b: Solver::new(
                spec,
                &solid,
                &wall_u,
                [1, 1, 1],
                CpuSimd::default(),
                LocalPeriodic,
            ),
        }
    }

    fn init(&mut self, f: impl Fn(usize, usize, usize) -> (T, [T; 3]) + Copy) {
        self.a.init_with(f);
        self.b.init_with(f);
    }

    /// Step both solvers, comparing every step for the first 5, every 50
    /// afterwards, and at the end; the probed force is compared every step.
    fn run_compare(&mut self, steps: usize, what: &str) {
        let lim = tol::<T>();
        let plim = tol_probe::<T>();
        let mut worst = compare_state(&self.a, &self.b, &format!("{what} t=0"));
        for s in 1..=steps {
            self.a.step();
            self.b.step();
            let (fa, fb) = (self.a.probed_force(), self.b.probed_force());
            for c in 0..L::D {
                let d = (fa[c].as_f64() - fb[c].as_f64()).abs();
                assert!(d <= plim, "{what} t={s}: probed_force[{c}] Δ = {d:e}");
                worst = worst.max(d);
            }
            if s <= 5 || s % 50 == 0 || s == steps {
                worst = worst.max(compare_state(&self.a, &self.b, &format!("{what} t={s}")));
            }
        }
        eprintln!("{what}: worst |Δ| over {steps} steps = {worst:e} (field tol {lim:e})");
    }
}

// ---------------------------------------------------------------------------
// 1. 2D TGV: fully periodic, TRT.
// ---------------------------------------------------------------------------

fn tgv_2d<T: Real>() {
    let (nx, ny) = (96usize, 64);
    let spec = GlobalSpec::<T> {
        dims: [nx, ny, 1],
        nu: 0.02,
        periodic: [true, true, false],
        ..Default::default()
    };
    let mut pair = Pair::<D2Q9, T>::new(&spec, &WallSpec::default());
    pair.init(move |x, y, _| {
        let kx = 2.0 * PI / nx as f64;
        let ky = 2.0 * PI / ny as f64;
        let u0 = 0.04;
        let (xf, yf) = (kx * x as f64, ky * y as f64);
        let rho = 1.0 + 0.01 * xf.cos() * yf.cos();
        (
            T::r(rho),
            [
                T::r(-u0 * xf.cos() * yf.sin()),
                T::r(u0 * xf.sin() * yf.cos()),
                T::zero(),
            ],
        )
    });
    pair.run_compare(400, "tgv-2d");
}

#[test]
fn tgv_2d_f64() {
    tgv_2d::<f64>();
}

#[test]
fn tgv_2d_f32() {
    tgv_2d::<f32>();
}

// ---------------------------------------------------------------------------
// 2. 2D lid-driven cavity: still walls + moving-wall bounce-back.
// ---------------------------------------------------------------------------

fn cavity_2d<T: Real>() {
    let n = 64usize;
    let mut walls = WallSpec::<T>::default();
    for f in [Face::XNeg, Face::XPos, Face::YNeg, Face::YPos] {
        walls.is_wall[f.index()] = true;
    }
    walls.u[Face::YPos.index()] = [T::r(0.1), T::zero(), T::zero()];
    let spec = GlobalSpec::<T> {
        dims: [n, n, 1],
        nu: 0.02,
        periodic: [false, false, false],
        ..Default::default()
    };
    let mut pair = Pair::<D2Q9, T>::new(&spec, &walls);
    pair.run_compare(400, "cavity-2d");
}

#[test]
fn cavity_2d_f64() {
    cavity_2d::<f64>();
}

#[test]
fn cavity_2d_f32() {
    cavity_2d::<f32>();
}

// ---------------------------------------------------------------------------
// 3. 2D channel: Zou–He velocity inlet with a per-node profile → Outflow.
// ---------------------------------------------------------------------------

fn channel_2d<T: Real>() {
    let (nx, ny) = (120usize, 48);
    let mut walls = WallSpec::<T>::default();
    walls.is_wall[Face::YNeg.index()] = true;
    walls.is_wall[Face::YPos.index()] = true;
    let mut faces = [FaceBC::Closed; 6];
    faces[Face::XNeg.index()] = FaceBC::Velocity {
        u: [T::r(0.05), T::zero(), T::zero()],
    };
    faces[Face::XPos.index()] = FaceBC::Outflow;
    let spec = GlobalSpec::<T> {
        dims: [nx, ny, 1],
        nu: 0.05,
        periodic: [false, false, false],
        faces,
        ..Default::default()
    };
    let mut pair = Pair::<D2Q9, T>::new(&spec, &walls);
    let profile: Vec<[T; 3]> = (0..ny)
        .map(|y| {
            let yy = y as f64 / (ny - 1) as f64;
            [T::r(0.08 * 4.0 * yy * (1.0 - yy)), T::zero(), T::zero()]
        })
        .collect();
    pair.a.set_inlet_profile(Face::XNeg, &profile);
    pair.b.set_inlet_profile(Face::XNeg, &profile);
    pair.run_compare(400, "channel-2d");
}

#[test]
fn channel_profile_outflow_2d_f64() {
    channel_2d::<f64>();
}

#[test]
fn channel_profile_outflow_2d_f32() {
    channel_2d::<f32>();
}

// ---------------------------------------------------------------------------
// 4. 2D cylinder + momentum-exchange force probe.
// ---------------------------------------------------------------------------

fn cylinder_probe_2d<T: Real>() {
    let (nx, ny) = (128usize, 64);
    let mut walls = WallSpec::<T>::default();
    walls.is_wall[Face::YNeg.index()] = true;
    walls.is_wall[Face::YPos.index()] = true;
    let mut faces = [FaceBC::Closed; 6];
    faces[Face::XNeg.index()] = FaceBC::Velocity {
        u: [T::r(0.05), T::zero(), T::zero()],
    };
    faces[Face::XPos.index()] = FaceBC::Outflow;
    let spec = GlobalSpec::<T> {
        dims: [nx, ny, 1],
        nu: 0.02,
        periodic: [false, false, false],
        faces,
        ..Default::default()
    };
    let mut pair = Pair::<D2Q9, T>::new(&spec, &walls);
    let (cx, cy, r) = (32.0f64, ny as f64 / 2.0 - 0.3, 7.3f64);
    let inside = move |x: usize, y: usize| {
        let (dx, dy) = (x as f64 - cx, y as f64 - cy);
        dx * dx + dy * dy < r * r
    };
    for y in 0..ny {
        for x in 0..nx {
            if inside(x, y) {
                pair.a.set_solid(x, y, 0);
                pair.b.set_solid(x, y, 0);
            }
        }
    }
    pair.a.set_force_probe(move |x, y, _| inside(x, y));
    pair.b.set_force_probe(move |x, y, _| inside(x, y));
    pair.run_compare(300, "cylinder-2d");
}

#[test]
fn cylinder_probe_2d_f64() {
    cylinder_probe_2d::<f64>();
}

#[test]
fn cylinder_probe_2d_f32() {
    cylinder_probe_2d::<f32>();
}

// ---------------------------------------------------------------------------
// 5. 2D Guo forcing: uniform + per-cell field, BGK, rewritten every step.
// ---------------------------------------------------------------------------

fn cell_force_2d<T: Real>() {
    let (nx, ny) = (64usize, 48);
    let spec = GlobalSpec::<T> {
        dims: [nx, ny, 1],
        nu: 0.03,
        collision: CollisionKind::Bgk,
        periodic: [true, true, false],
        force: [T::r(1e-5), T::r(-5e-6), T::zero()],
        ..Default::default()
    };
    let mut pair = Pair::<D2Q9, T>::new(&spec, &WallSpec::default());
    let k = 2.0 * PI / nx as f64;
    pair.init(move |x, y, _| {
        let (xf, yf) = (k * x as f64, k * y as f64);
        (
            T::one(),
            [T::r(0.02 * yf.sin()), T::r(0.02 * xf.sin()), T::zero()],
        )
    });
    let pat = move |x: usize, y: usize, t: usize| {
        let ph = t as f64 * 0.01;
        [
            T::r(1e-5 * ((k * x as f64) + ph).sin() * (k * y as f64).cos()),
            T::r(1e-5 * (k * x as f64).cos() * ((k * y as f64) - ph).sin()),
            T::zero(),
        ]
    };
    let lim = tol::<T>();
    let mut worst = 0.0f64;
    for s in 0..250usize {
        for solver_idx in 0..2 {
            let fields = if solver_idx == 0 {
                pair.a.fields_mut(0)
            } else {
                pair.b.fields_mut(0)
            };
            let ff = fields
                .force_field
                .get_or_insert_with(|| vec![[T::zero(); 3]; nx * ny]);
            for y in 0..ny {
                for x in 0..nx {
                    ff[y * nx + x] = pat(x, y, s);
                }
            }
        }
        pair.a.step();
        pair.b.step();
        if s % 25 == 0 || s == 249 {
            worst = worst.max(compare_state(&pair.a, &pair.b, &format!("cell-force t={s}")));
        }
    }
    eprintln!("cell-force-2d: worst |Δ| = {worst:e} (tol {lim:e})");
}

#[test]
fn cell_force_bgk_2d_f64() {
    cell_force_2d::<f64>();
}

#[test]
fn cell_force_bgk_2d_f32() {
    cell_force_2d::<f32>();
}

// ---------------------------------------------------------------------------
// 6. 2D ConvectiveOutflow: the stale-slot (previous post-collide) path.
// ---------------------------------------------------------------------------

fn convective_2d<T: Real>() {
    let (nx, ny) = (120usize, 48);
    let mut walls = WallSpec::<T>::default();
    walls.is_wall[Face::YNeg.index()] = true;
    walls.is_wall[Face::YPos.index()] = true;
    let mut faces = [FaceBC::Closed; 6];
    faces[Face::XNeg.index()] = FaceBC::Velocity {
        u: [T::r(0.06), T::zero(), T::zero()],
    };
    faces[Face::XPos.index()] = FaceBC::Convective { u_conv: T::r(0.06) };
    let spec = GlobalSpec::<T> {
        dims: [nx, ny, 1],
        nu: 0.03,
        periodic: [false, false, false],
        faces,
        ..Default::default()
    };
    let mut pair = Pair::<D2Q9, T>::new(&spec, &walls);
    pair.run_compare(400, "convective-2d");
}

#[test]
fn convective_outflow_2d_f64() {
    convective_2d::<f64>();
}

#[test]
fn convective_outflow_2d_f32() {
    convective_2d::<f32>();
}

// ---------------------------------------------------------------------------
// 7. 3D TGV (D3Q19, fully periodic, TRT).
// ---------------------------------------------------------------------------

fn tgv_3d<T: Real>() {
    let n = 32usize;
    let spec = GlobalSpec::<T> {
        dims: [n, n, n],
        nu: 0.02,
        periodic: [true, true, true],
        ..Default::default()
    };
    let mut pair = Pair::<D3Q19, T>::new(&spec, &WallSpec::default());
    let k = 2.0 * PI / n as f64;
    pair.init(move |x, y, z| {
        let (xf, yf, zf) = (k * x as f64, k * y as f64, k * z as f64);
        let u0 = 0.04;
        let rho = 1.0
            + 3.0 * u0 * u0 / 16.0 * ((2.0 * xf).cos() + (2.0 * yf).cos()) * ((2.0 * zf).cos() + 2.0);
        (
            T::r(rho),
            [
                T::r(u0 * xf.sin() * yf.cos() * zf.cos()),
                T::r(-u0 * xf.cos() * yf.sin() * zf.cos()),
                T::zero(),
            ],
        )
    });
    pair.run_compare(150, "tgv-3d");
}

#[test]
fn tgv_3d_f64() {
    tgv_3d::<f64>();
}

#[test]
fn tgv_3d_f32() {
    tgv_3d::<f32>();
}

// ---------------------------------------------------------------------------
// 8. 3D duct: inlet profile → Outflow, four wall faces (D3Q19 open-face +
//    bounce-back + Zou–He 5-unknown reconstruction all at once).
// ---------------------------------------------------------------------------

fn duct_3d<T: Real>() {
    let (nx, ny, nz) = (48usize, 20, 20);
    let mut walls = WallSpec::<T>::default();
    for f in [Face::YNeg, Face::YPos, Face::ZNeg, Face::ZPos] {
        walls.is_wall[f.index()] = true;
    }
    let mut faces = [FaceBC::Closed; 6];
    faces[Face::XNeg.index()] = FaceBC::Velocity {
        u: [T::r(0.05), T::zero(), T::zero()],
    };
    faces[Face::XPos.index()] = FaceBC::Outflow;
    let spec = GlobalSpec::<T> {
        dims: [nx, ny, nz],
        nu: 0.05,
        periodic: [false, false, false],
        faces,
        ..Default::default()
    };
    let mut pair = Pair::<D3Q19, T>::new(&spec, &walls);
    let prof = move |c1: usize, c2: usize| {
        let fy = (c1 as f64) / (ny - 1) as f64;
        let fz = (c2 as f64) / (nz - 1) as f64;
        [
            T::r(0.08 * 16.0 * fy * (1.0 - fy) * fz * (1.0 - fz)),
            T::zero(),
            T::zero(),
        ]
    };
    pair.a.set_inlet_profile_with(Face::XNeg, prof);
    pair.b.set_inlet_profile_with(Face::XNeg, prof);
    pair.run_compare(200, "duct-3d");
}

#[test]
fn duct_3d_f64() {
    duct_3d::<f64>();
}

#[test]
fn duct_3d_f32() {
    duct_3d::<f32>();
}

// ---------------------------------------------------------------------------
// Subdomain / halo integration: InProcess decomposition over CpuSimd must
// match the monolithic CpuScalar reference (band split is local to each
// part; halo exchange stays at the phase boundary).
// ---------------------------------------------------------------------------

/// 2D cylinder + probe straddling the 2×2 seam, inlet profile crossing the
/// y-seam, Outflow at the x-seam edge.
#[test]
fn split_2x2_cpusimd_matches_monolithic_cpuscalar() {
    type T = f64;
    let (nx, ny) = (96usize, 64);
    let mut walls = WallSpec::<T>::default();
    walls.is_wall[Face::YNeg.index()] = true;
    walls.is_wall[Face::YPos.index()] = true;
    let mut faces = [FaceBC::Closed; 6];
    faces[Face::XNeg.index()] = FaceBC::Velocity {
        u: [0.05, 0.0, 0.0],
    };
    faces[Face::XPos.index()] = FaceBC::Outflow;
    let spec = GlobalSpec::<T> {
        dims: [nx, ny, 1],
        nu: 0.02,
        periodic: [false, false, false],
        faces,
        ..Default::default()
    };
    let (solid, wall_u) = build_wall_rims(2, spec.dims, &walls);
    let mut mono: Solver<D2Q9, T, CpuScalar, LocalPeriodic> = Solver::new(
        &spec,
        &solid,
        &wall_u,
        [1, 1, 1],
        CpuScalar::default(),
        LocalPeriodic,
    );
    let mut split: Solver<D2Q9, T, CpuSimd, InProcess> = Solver::new(
        &spec,
        &solid,
        &wall_u,
        [2, 2, 1],
        CpuSimd::default(),
        InProcess,
    );
    // Cylinder centred on the 2×2 corner so bounce-back + probe links cross
    // both seams.
    let (cx, cy, r) = (nx as f64 / 2.0 - 0.2, ny as f64 / 2.0 + 0.4, 6.7f64);
    let inside = move |x: usize, y: usize| {
        let (dx, dy) = (x as f64 - cx, y as f64 - cy);
        dx * dx + dy * dy < r * r
    };
    for y in 0..ny {
        for x in 0..nx {
            if inside(x, y) {
                mono.set_solid(x, y, 0);
                split.set_solid(x, y, 0);
            }
        }
    }
    mono.set_force_probe(move |x, y, _| inside(x, y));
    split.set_force_probe(move |x, y, _| inside(x, y));
    let profile: Vec<[T; 3]> = (0..ny)
        .map(|y| {
            let yy = y as f64 / (ny - 1) as f64;
            [0.08 * 4.0 * yy * (1.0 - yy), 0.0, 0.0]
        })
        .collect();
    mono.set_inlet_profile(Face::XNeg, &profile);
    split.set_inlet_profile(Face::XNeg, &profile);
    let mut worst = 0.0f64;
    for s in 1..=300usize {
        mono.step();
        split.step();
        let (fa, fb) = (mono.probed_force(), split.probed_force());
        for c in 0..2 {
            let d = (fa[c] - fb[c]).abs();
            assert!(d <= 1e-11, "t={s}: probed_force[{c}] Δ = {d:e}");
            worst = worst.max(d);
        }
        if s <= 3 || s % 50 == 0 || s == 300 {
            worst = worst.max(compare_state(&mono, &split, &format!("2x2 t={s}")));
        }
    }
    eprintln!("split-2x2 cylinder: worst |Δ| = {worst:e}");
}

/// Fully periodic 2D TGV on a 2×2 CpuSimd decomposition (every band edge is
/// also a halo edge somewhere) vs monolithic CpuScalar.
#[test]
fn split_2x2_tgv_periodic_cpusimd_matches_cpuscalar() {
    type T = f64;
    let (nx, ny) = (64usize, 64);
    let spec = GlobalSpec::<T> {
        dims: [nx, ny, 1],
        nu: 0.02,
        periodic: [true, true, false],
        ..Default::default()
    };
    let mut mono: Solver<D2Q9, T, CpuScalar, LocalPeriodic> = Solver::new(
        &spec,
        &[],
        &[],
        [1, 1, 1],
        CpuScalar::default(),
        LocalPeriodic,
    );
    let mut split: Solver<D2Q9, T, CpuSimd, InProcess> =
        Solver::new(&spec, &[], &[], [2, 2, 1], CpuSimd::default(), InProcess);
    let init = move |x: usize, y: usize, _z: usize| {
        let kx = 2.0 * PI / nx as f64;
        let ky = 2.0 * PI / ny as f64;
        let (xf, yf) = (kx * x as f64, ky * y as f64);
        (
            1.0 + 0.01 * xf.cos() * yf.cos(),
            [-0.04 * xf.cos() * yf.sin(), 0.04 * xf.sin() * yf.cos(), 0.0],
        )
    };
    mono.init_with(init);
    split.init_with(init);
    let mut worst = compare_state(&mono, &split, "2x2 tgv t=0");
    for s in 1..=300usize {
        mono.step();
        split.step();
        if s <= 3 || s % 50 == 0 || s == 300 {
            worst = worst.max(compare_state(&mono, &split, &format!("2x2 tgv t={s}")));
        }
    }
    eprintln!("split-2x2 tgv: worst |Δ| = {worst:e}");
}

/// Two-pass streaming (the interior/boundary split that overlaps async
/// halo exchanges) must hold for the fused backend too: destination cells
/// are partitioned by the shell ranges, sources are ring-collided
/// redundantly per shell with identical inputs. Cylinder + probe exercises
/// the partial-range bounce-back and probe-fold paths.
#[test]
fn two_pass_streaming_cpusimd_matches_cpuscalar() {
    type T = f64;
    let (nx, ny) = (96usize, 48);
    let mut walls = WallSpec::<T>::default();
    walls.is_wall[Face::YNeg.index()] = true;
    walls.is_wall[Face::YPos.index()] = true;
    let mut faces = [FaceBC::Closed; 6];
    faces[Face::XNeg.index()] = FaceBC::Velocity {
        u: [0.05, 0.0, 0.0],
    };
    faces[Face::XPos.index()] = FaceBC::Outflow;
    let spec = GlobalSpec::<T> {
        dims: [nx, ny, 1],
        nu: 0.02,
        periodic: [false, false, false],
        faces,
        ..Default::default()
    };
    let (solid, wall_u) = build_wall_rims(2, spec.dims, &walls);
    let mut a: Solver<D2Q9, T, CpuScalar, LocalPeriodic> = Solver::new(
        &spec,
        &solid,
        &wall_u,
        [1, 1, 1],
        CpuScalar::default(),
        LocalPeriodic,
    );
    let mut b: Solver<D2Q9, T, CpuSimd, LocalPeriodic> = Solver::new(
        &spec,
        &solid,
        &wall_u,
        [1, 1, 1],
        CpuSimd::default(),
        LocalPeriodic,
    );
    b.set_two_pass(true);
    let (cx, cy, r) = (28.0f64, ny as f64 / 2.0 - 0.4, 6.1f64);
    let inside = move |x: usize, y: usize| {
        let (dx, dy) = (x as f64 - cx, y as f64 - cy);
        dx * dx + dy * dy < r * r
    };
    for y in 0..ny {
        for x in 0..nx {
            if inside(x, y) {
                a.set_solid(x, y, 0);
                b.set_solid(x, y, 0);
            }
        }
    }
    a.set_force_probe(move |x, y, _| inside(x, y));
    b.set_force_probe(move |x, y, _| inside(x, y));
    let mut worst = 0.0f64;
    for s in 1..=200usize {
        a.step();
        b.step();
        let (fa, fb) = (a.probed_force(), b.probed_force());
        for c in 0..2 {
            let d = (fa[c] - fb[c]).abs();
            assert!(d <= 1e-11, "t={s}: probed_force[{c}] Δ = {d:e}");
            worst = worst.max(d);
        }
        if s % 50 == 0 || s == 200 {
            worst = worst.max(compare_state(&a, &b, &format!("two-pass t={s}")));
        }
    }
    eprintln!("two-pass: worst |Δ| = {worst:e}");
}

/// 3D duct on a 2×2×1 CpuSimd decomposition vs monolithic CpuScalar: the
/// z-band split must stay local to each part while x/y halos cross seams.
#[test]
fn split_2x2x1_duct_3d_cpusimd_matches_cpuscalar() {
    type T = f64;
    let (nx, ny, nz) = (32usize, 16, 16);
    let mut walls = WallSpec::<T>::default();
    for f in [Face::YNeg, Face::YPos, Face::ZNeg, Face::ZPos] {
        walls.is_wall[f.index()] = true;
    }
    let mut faces = [FaceBC::Closed; 6];
    faces[Face::XNeg.index()] = FaceBC::Velocity {
        u: [0.05, 0.0, 0.0],
    };
    faces[Face::XPos.index()] = FaceBC::Outflow;
    let spec = GlobalSpec::<T> {
        dims: [nx, ny, nz],
        nu: 0.05,
        periodic: [false, false, false],
        faces,
        ..Default::default()
    };
    let (solid, wall_u) = build_wall_rims(3, spec.dims, &walls);
    let mut mono: Solver<D3Q19, T, CpuScalar, LocalPeriodic> = Solver::new(
        &spec,
        &solid,
        &wall_u,
        [1, 1, 1],
        CpuScalar::default(),
        LocalPeriodic,
    );
    let mut split: Solver<D3Q19, T, CpuSimd, InProcess> = Solver::new(
        &spec,
        &solid,
        &wall_u,
        [2, 2, 1],
        CpuSimd::default(),
        InProcess,
    );
    let mut worst = compare_state(&mono, &split, "2x2x1 duct t=0");
    for s in 1..=150usize {
        mono.step();
        split.step();
        if s <= 3 || s % 50 == 0 || s == 150 {
            worst = worst.max(compare_state(&mono, &split, &format!("2x2x1 duct t={s}")));
        }
    }
    eprintln!("split-2x2x1 duct: worst |Δ| = {worst:e}");
}
