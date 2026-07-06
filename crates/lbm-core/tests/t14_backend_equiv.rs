#![cfg(feature = "gpu")]
//! T14: backend equivalence, CpuScalar vs Wgpu (COMPETITIVE_SPEC.md §4, R2).
//!
//! The same scenario runs on the CPU reference backend and on the wgpu
//! backend from the same f32 initial state; fields must agree to
//! **L∞ ≤ 1e-5 relative** (per-field max-norm) and the f64-accumulated
//! diagnostics (mass / momentum / probe force) to **≤ 1e-4 relative**.
//!
//! Operator-order note (the GPU_EVALUATION.md §2 recipe): the proto's
//! pull-fused kernel computes `C∘S` per dispatch and therefore needed the
//! `(C∘S)^k ∘ C = C ∘ (S∘C)^k` identity (pre-collided upload) to compare
//! against the CPU's `S∘C` — and that identity breaks once a boundary pass
//! must run *between* S and C. The V2 backend instead fuses in **push**
//! form: one dispatch collides its own cell and scatters to the neighbours,
//! i.e. computes exactly `S∘C`, and the open-face BC kernels run after it
//! just like `CpuScalar::apply_open_faces`. Per-step operator order is
//! identical, so states are compared **directly at every checkpoint — no
//! recipe needed**. Measured residual (this suite, M5 Max / Metal): field
//! L∞ stays orders of magnitude inside the 1e-5 line over 300–500 steps —
//! pure f32 rounding from the Metal compiler's FMA/reassociation —
//! confirming the direct comparison is valid.
//!
//! Coverage (7 configurations × 300–500 steps; 1–6 hold the strict 1e-5
//! field line, satisfying the "≥6 configurations" acceptance):
//!   1. TGV, fully periodic, TRT              (fused kernel + wrap)
//!   2. Lid-driven cavity                     (still + moving bounce-back)
//!   3. Channel, Zou–He inlet profile → Outflow
//!   4. Cylinder + momentum-exchange probe    (probe force diagnostics)
//!   5. Per-cell + uniform Guo force, BGK     (force paths, ω− = ω+)
//!   6. Channel → ConvectiveOutflow           (previous-value convention +
//!                                             mass pinning via edge stash)
//!   7. Pressure-driven channel (Zou–He pressure both ends) — documented
//!      relaxed line 1e-4 + a CPU-vs-CPU 1-ulp control test proving the
//!      pressure closure's intrinsic sensitivity sets the floor

use lbm_core::lattice::D2Q9;
use lbm_core::prelude::*;
use std::sync::{Arc, OnceLock};

fn ctx() -> Arc<GpuContext> {
    static CTX: OnceLock<Arc<GpuContext>> = OnceLock::new();
    CTX.get_or_init(|| GpuContext::new().expect("T14 requires a GPU adapter (run on a GPU host)"))
        .clone()
}

type Cpu = Solver<D2Q9, f32, CpuScalar, LocalPeriodic>;

/// Field acceptance: L∞(Δ) / max(‖ref‖∞, floor) ≤ 1e-5 (the R2 line).
const FIELD_TOL: f64 = 1e-5;
/// Diagnostics acceptance: ≤ 1e-4 relative (per-quantity absolute floor).
const DIAG_TOL: f64 = 1e-4;

/// L∞(a-b) / max(‖a‖∞, floor), accumulated in f64.
fn linf_rel(a: &[f32], b: &[f32], floor: f64) -> f64 {
    assert_eq!(a.len(), b.len());
    let mut d = 0.0f64;
    let mut m = 0.0f64;
    for (x, y) in a.iter().zip(b) {
        d = d.max((*x as f64 - *y as f64).abs());
        m = m.max((*x as f64).abs());
    }
    d / m.max(floor)
}

fn assert_rel(a: f64, b: f64, floor: f64, what: &str) {
    let d = (a - b).abs();
    let lim = DIAG_TOL * a.abs().max(floor);
    assert!(
        d <= lim,
        "{what}: |Δ| = {d:e} > {lim:e} (cpu {a:e}, gpu {b:e})"
    );
}

struct Pair {
    cpu: Cpu,
    gpu: GpuSolver<D2Q9>,
    /// Characteristic velocity, used as the momentum-comparison floor
    /// (`n_cells * u_char`): total momentum can legitimately be ~0 (TGV).
    u_char: f64,
}

impl Pair {
    fn new(spec: &GlobalSpec<f32>, walls: &WallSpec<f32>, u_char: f64) -> Self {
        let (solid, wall_u) = build_wall_rims(2, spec.dims, walls);
        let cpu = Cpu::new(
            spec,
            &solid,
            &wall_u,
            [1, 1, 1],
            CpuScalar::default(),
            LocalPeriodic,
        );
        let gpu = GpuSolver::new(spec, &solid, &wall_u, ctx());
        Self { cpu, gpu, u_char }
    }

    fn init(&mut self, f: impl Fn(usize, usize) -> (f32, [f32; 3]) + Copy) {
        self.cpu.init_with(move |x, y, _| f(x, y));
        self.gpu.init_with(move |x, y, _| f(x, y));
    }

    fn run(&mut self, steps: usize) {
        self.cpu.run(steps);
        self.gpu.run(steps);
    }

    /// Field comparison with the R2/proto metric (GPU_EVALUATION.md §2):
    /// velocity `max(|Δux|, |Δuy|) / max‖u_cpu‖` and `|Δrho| / max|rho|`,
    /// both ≤ `tol` (1e-5 for every case except the pressure BC — see the
    /// pressure test for the measured justification). Population planes are
    /// checked as a supplementary, state-normalised statement,
    /// `max_q |Δf_q| / max_q ‖f_q‖∞ ≤ max(tol, 1e-4)` over fluid cells only
    /// (solid cells hold dead ping-pong junk on both backends, by V1
    /// mechanics): the drift floor of two independent f32 compilations
    /// (Metal fast-math FMA/reassociation) saturates around 2–4e-7 absolute
    /// here, which a *per-plane* 1e-5 relative line would sit below for the
    /// small-magnitude planes (measured injection ~3e-8/step with uniformly
    /// distributed argmax — noise, not a systematic defect).
    fn check_tol(&mut self, what: &str, tol: f64) {
        let t = self.cpu.time();
        // Velocity, proto metric.
        let (uxa, uxb) = (self.cpu.gather_ux(), self.gpu.gather_ux());
        let (uya, uyb) = (self.cpu.gather_uy(), self.gpu.gather_uy());
        let mut du = 0.0f64;
        let mut umax = 0.0f64;
        for i in 0..uxa.len() {
            du = du
                .max((uxa[i] as f64 - uxb[i] as f64).abs())
                .max((uya[i] as f64 - uyb[i] as f64).abs());
            let s = (uxa[i] as f64).hypot(uya[i] as f64);
            umax = umax.max(s);
        }
        let ru = du / umax.max(1e-6);
        eprintln!("{what} t={t}: u L∞/max|u| = {ru:.3e}");
        assert!(
            ru <= tol,
            "{what} t={t}: velocity L∞ rel = {ru:e} > {tol:e}"
        );
        // Density.
        let (ra, rb) = (self.cpu.gather_rho(), self.gpu.gather_rho());
        let rr = linf_rel(&ra, &rb, 1.0);
        eprintln!("{what} t={t}: rho L∞ rel = {rr:.3e}");
        assert!(rr <= tol, "{what} t={t}: rho L∞ rel = {rr:e} > {tol:e}");
        // Populations, state-normalised, fluid cells only.
        let dims = self.cpu.dims();
        let fluid: Vec<bool> = (0..dims[0] * dims[1])
            .map(|i| !self.cpu.is_solid(i % dims[0], i / dims[0], 0))
            .collect();
        let mut df = 0.0f64;
        let mut fmax = 0.0f64;
        for q in 0..9 {
            let (fa, fb) = (self.cpu.gather_f(q), self.gpu.gather_f(q));
            for ((x, y), &fl) in fa.iter().zip(&fb).zip(&fluid) {
                if !fl {
                    continue;
                }
                df = df.max((*x as f64 - *y as f64).abs());
                fmax = fmax.max((*x as f64).abs());
            }
        }
        let rf = df / fmax.max(1e-6);
        let ftol = tol.max(1e-4);
        eprintln!("{what} t={t}: f (state-norm) L∞ rel = {rf:.3e}");
        assert!(
            rf <= ftol,
            "{what} t={t}: population L∞ rel = {rf:e} > {ftol:e}"
        );
        // Diagnostics (f64 accumulation on both sides; the GPU reduce is
        // the identical host-side V1 loop over read-back populations).
        let n = (dims[0] * dims[1]) as f64;
        assert_rel(
            self.cpu.total_mass() as f64,
            self.gpu.total_mass() as f64,
            1.0,
            &format!("{what} t={t}: total_mass"),
        );
        let (pa, pb) = (self.cpu.total_momentum(), self.gpu.total_momentum());
        for c in 0..2 {
            assert_rel(
                pa[c] as f64,
                pb[c] as f64,
                n * self.u_char,
                &format!("{what} t={t}: momentum[{c}]"),
            );
        }
    }

    fn check(&mut self, what: &str) {
        self.check_tol(what, FIELD_TOL);
    }

    fn check_probe(&mut self, what: &str) {
        let t = self.cpu.time();
        let (fa, fb) = (self.cpu.probed_force(), self.gpu.probed_force());
        let scale = fa[0].abs().max(fa[1].abs()).max(1e-6) as f64;
        for c in 0..2 {
            assert_rel(
                fa[c] as f64,
                fb[c] as f64,
                scale,
                &format!("{what} t={t}: probed_force[{c}]"),
            );
        }
    }
}

// ---------------------------------------------------------------------------
// Shared pressure-channel scenario (case 7 + its sensitivity control):
// Poiseuille-initialised channel driven by Zou–He pressure faces, with
// delta_rho consistent with the profile (u0 = cs2 drho H^2 / (8 nu (nx-1)))
// and inside the BC's validated envelope (validation_open_bc.rs uses 2e-3).
// ---------------------------------------------------------------------------

struct PressureChannel {
    spec: GlobalSpec<f32>,
    walls: WallSpec<f32>,
    u0: f64,
    drho: f64,
}

fn pressure_channel(rho_out_bump_ulp: bool) -> PressureChannel {
    let (nx, ny) = (96usize, 48usize);
    let mut walls = WallSpec::default();
    walls.is_wall[Face::YNeg.index()] = true;
    walls.is_wall[Face::YPos.index()] = true;
    let nu = 0.02f64;
    let u0 = 0.1f64;
    let h = (ny - 2) as f64;
    let drho = u0 * 8.0 * nu * (nx - 1) as f64 / ((1.0 / 3.0) * h * h);
    let mut rho_out = (1.0 - 0.5 * drho) as f32;
    if rho_out_bump_ulp {
        rho_out = f32::from_bits(rho_out.to_bits() + 1);
    }
    let mut faces = [FaceBC::Closed; 6];
    faces[Face::XNeg.index()] = FaceBC::Pressure {
        rho: (1.0 + 0.5 * drho) as f32,
    };
    faces[Face::XPos.index()] = FaceBC::Pressure { rho: rho_out };
    PressureChannel {
        spec: GlobalSpec {
            dims: [nx, ny, 1],
            nu,
            periodic: [false, false, false],
            faces,
            ..Default::default()
        },
        walls,
        u0,
        drho,
    }
}

fn pressure_channel_init(
    nx: usize,
    ny: usize,
    u0: f64,
    drho: f64,
    x: usize,
    y: usize,
) -> (f32, [f32; 3]) {
    let rho = 1.0 + 0.5 * drho - drho * x as f64 / (nx - 1) as f64;
    let yy = (y as f64 - 1.0) / (ny as f64 - 3.0);
    let ux = if y == 0 || y == ny - 1 {
        0.0
    } else {
        4.0 * u0 * yy * (1.0 - yy)
    };
    (rho as f32, [ux as f32, 0.0, 0.0])
}

// ---------------------------------------------------------------------------
// 1. TGV (fully periodic, TRT): the fused kernel + periodic wrap.
// ---------------------------------------------------------------------------

#[test]
fn t14_tgv_periodic() {
    let n = 128usize;
    let spec = GlobalSpec::<f32> {
        dims: [n, n, 1],
        nu: 0.02,
        periodic: [true, true, false],
        ..Default::default()
    };
    let mut pair = Pair::new(&spec, &WallSpec::default(), 0.05);
    let u0 = 0.05f64;
    let k = 2.0 * std::f64::consts::PI / n as f64;
    pair.init(move |x, y| {
        let (xf, yf) = (k * x as f64, k * y as f64);
        let rho = 1.0 - 3.0 * u0 * u0 / 4.0 * ((2.0 * xf).cos() + (2.0 * yf).cos());
        (
            rho as f32,
            [
                (-u0 * xf.cos() * yf.sin()) as f32,
                (u0 * xf.sin() * yf.cos()) as f32,
                0.0,
            ],
        )
    });
    for _ in 0..4 {
        pair.run(100);
        pair.check("TGV");
    }
}

// ---------------------------------------------------------------------------
// 2. Lid-driven cavity: still walls + moving-wall bounce-back.
// ---------------------------------------------------------------------------

#[test]
fn t14_lid_cavity() {
    let n = 96usize;
    let mut walls = WallSpec::default();
    for f in [Face::XNeg, Face::XPos, Face::YNeg, Face::YPos] {
        walls.is_wall[f.index()] = true;
    }
    walls.u[Face::YPos.index()] = [0.08, 0.0, 0.0];
    let spec = GlobalSpec::<f32> {
        dims: [n, n, 1],
        nu: 0.02,
        periodic: [false, false, false],
        ..Default::default()
    };
    let mut pair = Pair::new(&spec, &walls, 0.08);
    pair.run(200);
    pair.check("cavity");
    pair.run(200);
    pair.check("cavity");
}

// ---------------------------------------------------------------------------
// 3. Channel: Zou–He velocity inlet with a per-node profile → Outflow.
// ---------------------------------------------------------------------------

#[test]
fn t14_channel_inlet_profile_outflow() {
    let (nx, ny) = (160usize, 64usize);
    let mut walls = WallSpec::default();
    walls.is_wall[Face::YNeg.index()] = true;
    walls.is_wall[Face::YPos.index()] = true;
    let mut faces = [FaceBC::Closed; 6];
    faces[Face::XNeg.index()] = FaceBC::Velocity {
        u: [0.05, 0.0, 0.0],
    };
    faces[Face::XPos.index()] = FaceBC::Outflow;
    let spec = GlobalSpec::<f32> {
        dims: [nx, ny, 1],
        nu: 0.05,
        periodic: [false, false, false],
        faces,
        ..Default::default()
    };
    let mut pair = Pair::new(&spec, &walls, 0.08);
    let profile: Vec<[f32; 3]> = (0..ny)
        .map(|y| {
            let yy = y as f32 / (ny - 1) as f32;
            [0.08 * 4.0 * yy * (1.0 - yy), 0.0, 0.0]
        })
        .collect();
    pair.cpu.set_inlet_profile(Face::XNeg, &profile);
    pair.gpu.set_inlet_profile(Face::XNeg, &profile);
    for _ in 0..4 {
        pair.run(100);
        pair.check("channel");
    }
}

// ---------------------------------------------------------------------------
// 4. Cylinder + momentum-exchange force probe.
// ---------------------------------------------------------------------------

#[test]
fn t14_cylinder_probe() {
    let (nx, ny) = (160usize, 80usize);
    let mut walls = WallSpec::default();
    walls.is_wall[Face::YNeg.index()] = true;
    walls.is_wall[Face::YPos.index()] = true;
    let mut faces = [FaceBC::Closed; 6];
    faces[Face::XNeg.index()] = FaceBC::Velocity {
        u: [0.05, 0.0, 0.0],
    };
    faces[Face::XPos.index()] = FaceBC::Outflow;
    let spec = GlobalSpec::<f32> {
        dims: [nx, ny, 1],
        nu: 0.02,
        periodic: [false, false, false],
        faces,
        ..Default::default()
    };
    let mut pair = Pair::new(&spec, &walls, 0.05);
    let (cx, cy, r) = (40.0f64, ny as f64 / 2.0 - 0.3, 8.2f64);
    let inside = move |x: usize, y: usize| {
        let (dx, dy) = (x as f64 - cx, y as f64 - cy);
        dx * dx + dy * dy < r * r
    };
    for y in 0..ny {
        for x in 0..nx {
            if inside(x, y) {
                pair.cpu.set_solid(x, y, 0);
                pair.gpu.set_solid(x, y, 0);
            }
        }
    }
    pair.cpu.set_force_probe(move |x, y, _| inside(x, y));
    pair.gpu.set_force_probe(move |x, y, _| inside(x, y));
    for _ in 0..4 {
        pair.run(100);
        pair.check_probe("cylinder");
    }
    pair.check("cylinder");
}

// ---------------------------------------------------------------------------
// 5. Guo forcing, uniform + per-cell field, BGK.
// ---------------------------------------------------------------------------

#[test]
fn t14_cell_force_bgk() {
    let n = 128usize;
    let spec = GlobalSpec::<f32> {
        dims: [n, n, 1],
        nu: 0.03,
        collision: CollisionKind::Bgk,
        periodic: [true, true, false],
        force: [1e-5, 0.0, 0.0],
        ..Default::default()
    };
    let mut pair = Pair::new(&spec, &WallSpec::default(), 0.05);
    // Per-cell swirl force on top of the uniform force.
    let k = 2.0 * std::f64::consts::PI / n as f64;
    let field: Vec<[f32; 3]> = (0..n * n)
        .map(|i| {
            let (x, y) = (i % n, i / n);
            let (xf, yf) = (k * x as f64, k * y as f64);
            [(2e-5 * yf.sin()) as f32, (-2e-5 * xf.sin()) as f32, 0.0]
        })
        .collect();
    pair.cpu.set_body_force_field_values(&field);
    pair.gpu.set_force_field(field);
    pair.init(move |x, y| {
        let (xf, yf) = (k * x as f64, k * y as f64);
        (
            1.0,
            [(0.02 * yf.sin()) as f32, (0.02 * xf.sin()) as f32, 0.0],
        )
    });
    for _ in 0..3 {
        pair.run(100);
        pair.check("cell-force");
    }
}

// ---------------------------------------------------------------------------
// Gravity body force: device-resident `rho(x) * g` composed with the existing
// Guo force path. Ignored here because the Codex sandbox has no Metal adapter;
// PM/GPU hosts run it with `--features gpu -- --ignored`.
// ---------------------------------------------------------------------------

#[test]
#[ignore = "BENCH-PENDING: requires a native GPU adapter"]
fn t14_gravity_body_force_device_resident() {
    let n = 96usize;
    let spec = GlobalSpec::<f32> {
        dims: [n, n, 1],
        nu: 0.04,
        periodic: [true, true, false],
        force: [2e-6, 0.0, 0.0],
        ..Default::default()
    };
    let mut pair = Pair::new(&spec, &WallSpec::default(), 0.03);
    let k = 2.0 * std::f64::consts::PI / n as f64;
    pair.init(move |x, y| {
        let (xf, yf) = (k * x as f64, k * y as f64);
        let rho = 1.0 + 0.02 * (xf.sin() * yf.cos());
        (
            rho as f32,
            [(0.01 * yf.sin()) as f32, (-0.01 * xf.cos()) as f32, 0.0],
        )
    });
    let field: Vec<[f32; 3]> = (0..n * n)
        .map(|i| {
            let (x, y) = (i % n, i / n);
            let (xf, yf) = (k * x as f64, k * y as f64);
            [(3e-6 * yf.cos()) as f32, (2e-6 * xf.sin()) as f32, 0.0]
        })
        .collect();
    pair.cpu.set_body_force_field_values(&field);
    pair.gpu.set_force_field(field);
    pair.cpu.set_gravity([0.0, -4e-6, 0.0]);
    pair.gpu.set_gravity([0.0, -4e-6, 0.0]);
    for _ in 0..3 {
        pair.run(100);
        pair.check("gravity-body-force");
    }
}

// ---------------------------------------------------------------------------
// 6. ConvectiveOutflow: previous-value convention + mass pinning (the edge
//    stash path — the one place the push-fused kernel must reproduce V1's
//    stale-slot mechanics; see gpu/wgsl.rs module docs).
// ---------------------------------------------------------------------------

#[test]
fn t14_convective_outflow() {
    let (nx, ny) = (160usize, 64usize);
    let mut walls = WallSpec::default();
    walls.is_wall[Face::YNeg.index()] = true;
    walls.is_wall[Face::YPos.index()] = true;
    let mut faces = [FaceBC::Closed; 6];
    faces[Face::XNeg.index()] = FaceBC::Velocity {
        u: [0.05, 0.0, 0.0],
    };
    faces[Face::XPos.index()] = FaceBC::Convective { u_conv: 0.05 };
    let spec = GlobalSpec::<f32> {
        dims: [nx, ny, 1],
        nu: 0.05,
        periodic: [false, false, false],
        faces,
        ..Default::default()
    };
    let mut pair = Pair::new(&spec, &walls, 0.05);
    for _ in 0..4 {
        pair.run(100);
        pair.check("convective");
    }
}

// ---------------------------------------------------------------------------
// 7. Pressure-driven channel: Zou–He pressure faces on both ends.
//
// This case runs against a *documented relaxed field line of 1e-4* (the six
// configurations above satisfy the strict R2 1e-5 line, which meets the
// "≥6 configurations" acceptance). Reason, established by measurement:
//
// The pressure closure `un = 1 - closure/rho_bc` maps ~ulp(1)-scale
// arithmetic differences of the O(1) closure (Metal fast-math reciprocal
// division / reassociation, deterministic per compiler) *directly into a
// velocity* at the face — unlike the velocity BC, where the corresponding
// division noise lands in rho and enters f multiplied by the small u_n. The
// scenario then amplifies it: the drift probe showed the CPU↔GPU difference
// pinned at the pressure faces, present from step 2 at ~2.2e-7, growing to
// ~2.5e-6 by t=100 (velocity-relative 2.5e-5 at u0 = 0.1). The control test
// below reproduces the same growth curve **CPU-vs-CPU** by perturbing
// rho_bc by exactly 1 ulp — proving the sensitivity is intrinsic to the
// f32 Zou–He pressure closure, not a Wgpu implementation defect.
// ---------------------------------------------------------------------------

#[test]
fn t14_pressure_driven_channel() {
    let pc = pressure_channel(false);
    let [nx, ny, _] = pc.spec.dims;
    let (u0, drho) = (pc.u0, pc.drho);
    let mut pair = Pair::new(&pc.spec, &pc.walls, u0);
    pair.init(move |x, y| pressure_channel_init(nx, ny, u0, drho, x, y));
    for _ in 0..3 {
        pair.run(100);
        pair.check_tol("pressure-channel", 1e-4);
    }
}

/// Sensitivity control for case 7: the same scenario run CPU-vs-CPU with
/// `rho_out` perturbed by **1 ulp** (7e-8 relative — the scale of one
/// fast-math rounding difference in the closure). If a single-ulp input
/// difference already produces a drift comparable to the measured CPU↔GPU
/// difference, the relaxed line of case 7 is justified by the BC's
/// conditioning, not by the backend. Measured on M5 Max (t=100):
/// CPU-vs-CPU(1ulp) ≈ 1.5e-6 abs vs CPU↔GPU ≈ 2.5e-6 abs — same decade,
/// same face-pinned argmax, same growth shape.
#[test]
fn t14_pressure_bc_ulp_sensitivity_control() {
    let build = |bump: bool| {
        let pc = pressure_channel(bump);
        let [nx, ny, _] = pc.spec.dims;
        let (u0, drho) = (pc.u0, pc.drho);
        let (solid, wall_u) = build_wall_rims(2, pc.spec.dims, &pc.walls);
        let mut s = Cpu::new(
            &pc.spec,
            &solid,
            &wall_u,
            [1, 1, 1],
            CpuScalar::default(),
            LocalPeriodic,
        );
        s.init_with(move |x, y, _| pressure_channel_init(nx, ny, u0, drho, x, y));
        s
    };
    let mut a = build(false);
    let mut b = build(true);
    a.run(100);
    b.run(100);
    let (ua, ub) = (a.gather_ux(), b.gather_ux());
    let mut d = 0.0f64;
    for (x, y) in ua.iter().zip(&ub) {
        d = d.max((*x as f64 - *y as f64).abs());
    }
    eprintln!("pressure control: CPU-vs-CPU(1ulp rho_bc) t=100 ux dLinf = {d:.3e}");
    // The 1-ulp input drift must itself exceed the strict field line at
    // u0 = 0.1 — i.e. no f32 implementation pair could meet 1e-5 here —
    // and stay within the relaxed 1e-4 line that case 7 asserts.
    assert!(
        d > 1e-5 * 0.1,
        "scenario no longer amplifies 1-ulp input differences (drift {d:e}); \
         the relaxed tolerance of t14_pressure_driven_channel can be tightened"
    );
    assert!(d < 1e-4 * 0.1, "1-ulp drift {d:e} outgrew the relaxed line");
}
