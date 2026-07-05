//! Host-side f32 replica of lbm-core's initialisation path.
//!
//! lbm-core does not expose its population array, so to hand the GPU an
//! initial state on the *same trajectory* as `Simulation<f32>` we re-derive
//! it here with the same arithmetic, in the same order, in f32:
//!
//! 1. `f0 = feq(rho, u) + f_neq(grad u)` — the second-order consistent init
//!    (`Simulation::init_with`, periodic branch),
//! 2. moments recomputed from `f0` (`update_moments`), then
//! 3. one TRT collision (`collide_row`) applied on the host.
//!
//! Step 3 exists because the fused GPU kernel computes
//! `f_new = Collide(Stream(f))` per dispatch while the CPU steps
//! `Stream(Collide(f))`. With `g0 = Collide(f0)` uploaded, k GPU dispatches
//! yield `(C∘S)^k C f0 = C((S∘C)^k f0) = C(cpu_k)`, and since collision
//! preserves density and momentum exactly, the macroscopic fields match the
//! CPU's step-k fields one for one.
//!
//! Everything below is a line-for-line port of `lbm-core/src/sim.rs`
//! specialised to: f32, fully periodic, no solids, no body force.

use lbm_core::compat::lattice::{CX, CY, PAIRS, Q, W};

/// Lattice constants pre-converted to f32 (mirrors `Simulation::params`).
struct Consts {
    wr: [f32; Q],
    cxr: [f32; Q],
    cyr: [f32; Q],
}

fn consts() -> Consts {
    let mut c = Consts {
        wr: [0.0; Q],
        cxr: [0.0; Q],
        cyr: [0.0; Q],
    };
    for q in 0..Q {
        c.wr[q] = W[q] as f32;
        c.cxr[q] = CX[q] as f32;
        c.cyr[q] = CY[q] as f32;
    }
    c
}

/// TRT relaxation rates (magic 3/16, lbm-core's default), converted to f32
/// the same way `Simulation::from_config` + `params` do.
pub fn omegas_f32(nu: f64) -> (f32, f32) {
    let tau = 3.0 * nu + 0.5;
    let omega_p = 1.0 / tau;
    let magic = 3.0 / 16.0; // Collision::MAGIC_STD
    let lam_p = tau - 0.5;
    let omega_m = 1.0 / (magic / lam_p + 0.5);
    (omega_p as f32, omega_m as f32)
}

/// Taylor–Green initial condition, shared verbatim by the CPU reference
/// (`sim.init_with`) and the GPU host init so both consume identical f32
/// triples. Includes the analytic pressure field (see smoke_tgv.rs).
pub fn tgv_ic(x: usize, y: usize, n: usize, u0: f64) -> (f32, f32, f32) {
    let k = 2.0 * std::f64::consts::PI / n as f64;
    let (xf, yf) = (k * x as f64, k * y as f64);
    let rho = 1.0 - 3.0 * u0 * u0 / 4.0 * ((2.0 * xf).cos() + (2.0 * yf).cos());
    (
        rho as f32,
        (-u0 * xf.cos() * yf.sin()) as f32,
        (u0 * xf.sin() * yf.cos()) as f32,
    )
}

/// Equilibrium in deviation form (`feq_q - w_q`), port of `sim::equilibrium`.
fn equilibrium(c: &Consts, r: f32, vx: f32, vy: f32) -> [f32; Q] {
    let usq = vx * vx + vy * vy;
    let drho = r - 1.0;
    let mut feq = [0.0f32; Q];
    for q in 0..Q {
        let cu = c.cxr[q] * vx + c.cyr[q] * vy;
        feq[q] = c.wr[q] * (drho + r * (3.0 * cu + 4.5 * cu * cu - 1.5 * usq));
    }
    feq
}

/// Populations (deviation storage, cell-major AoS `f[(y*nx+x)*9+q]`) plus the
/// moments recomputed from them.
pub struct HostState {
    pub f_aos: Vec<f32>,
    pub rho: Vec<f32>,
    pub ux: Vec<f32>,
    pub uy: Vec<f32>,
}

/// Port of `Simulation::init_with` (periodic, no solids) + `update_moments`
/// for an n x n Taylor–Green vortex.
pub fn init_tgv_f32(n: usize, u0: f64, nu: f64) -> HostState {
    let c = consts();
    let ncells = n * n;
    let mut rho = vec![0.0f32; ncells];
    let mut ux = vec![0.0f32; ncells];
    let mut uy = vec![0.0f32; ncells];
    for y in 0..n {
        for x in 0..n {
            let (r, vx, vy) = tgv_ic(x, y, n, u0);
            let i = y * n + x;
            rho[i] = r;
            ux[i] = vx;
            uy[i] = vy;
        }
    }

    // f = feq + f_neq, with f_neq from periodic central differences.
    let tau = (3.0 * nu + 0.5) as f32;
    let mut f = vec![0.0f32; ncells * Q];
    for y in 0..n {
        for x in 0..n {
            let i = y * n + x;
            let o = i * Q;
            let feq = equilibrium(&c, rho[i], ux[i], uy[i]);
            f[o..o + Q].copy_from_slice(&feq);
            let xp = (x + 1) % n;
            let xm = (x + n - 1) % n;
            let yp = (y + 1) % n;
            let ym = (y + n - 1) % n;
            let duxdx = (ux[y * n + xp] - ux[y * n + xm]) * 0.5;
            let duydx = (uy[y * n + xp] - uy[y * n + xm]) * 0.5;
            let duxdy = (ux[yp * n + x] - ux[ym * n + x]) * 0.5;
            let duydy = (uy[yp * n + x] - uy[ym * n + x]) * 0.5;
            let div = duxdx + duydy;
            for q in 0..Q {
                let (cx, cy) = (c.cxr[q], c.cyr[q]);
                let ccgu = cx * cx * duxdx + cx * cy * (duydx + duxdy) + cy * cy * duydy;
                let fneq = -c.wr[q] * rho[i] * tau * (3.0 * ccgu - div);
                f[o + q] += fneq;
            }
        }
    }

    // Recompute moments from f exactly like `update_moments` (the CPU's first
    // collision consumes these, not the analytic closure values).
    for i in 0..ncells {
        let o = i * Q;
        let mut dr = 0.0f32;
        let mut mx = 0.0f32;
        let mut my = 0.0f32;
        for q in 0..Q {
            let fq = f[o + q];
            dr += fq;
            mx += c.cxr[q] * fq;
            my += c.cyr[q] * fq;
        }
        let r = 1.0 + dr;
        rho[i] = r;
        let inv = 1.0 / r;
        ux[i] = mx * inv;
        uy[i] = my * inv;
    }

    HostState {
        f_aos: f,
        rho,
        ux,
        uy,
    }
}

/// One TRT collision over all cells, port of `sim::collide_row` (no force).
pub fn collide_trt_f32(f: &mut [f32], rho: &[f32], ux: &[f32], uy: &[f32], op: f32, om: f32) {
    let c = consts();
    for i in 0..rho.len() {
        let o = i * Q;
        let feq = equilibrium(&c, rho[i], ux[i], uy[i]);
        f[o] -= op * (f[o] - feq[0]);
        for (a, b) in PAIRS {
            let (fa, fb) = (f[o + a], f[o + b]);
            let fp = 0.5 * (fa + fb);
            let fm = 0.5 * (fa - fb);
            let ep = 0.5 * (feq[a] + feq[b]);
            let em = 0.5 * (feq[a] - feq[b]);
            let rp = op * (fp - ep);
            let rm = om * (fm - em);
            f[o + a] = fa - rp - rm;
            f[o + b] = fb - rp + rm;
        }
    }
}

/// Cell-major AoS -> direction-major SoA (`f[q*n + i]`), the GPU layout.
pub fn aos_to_soa(f_aos: &[f32], ncells: usize) -> Vec<f32> {
    let mut soa = vec![0.0f32; ncells * Q];
    for i in 0..ncells {
        for q in 0..Q {
            soa[q * ncells + i] = f_aos[i * Q + q];
        }
    }
    soa
}
