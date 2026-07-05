//! T15 — 3D (D3Q19) physics validation (VALIDATION.md T15, COMPETITIVE_SPEC R1).
//!
//! 1. z-invariant 2D-TGV degeneracy: D3Q19 on N×N×4 (z periodic) must match
//!    D2Q9 field-by-field (f64, ≤ 1e-12) — the sharpest smoke for streaming /
//!    weights / symmetry bugs. A second degeneracy angle drives the 3D Zou–He
//!    faces (velocity inlet with profile + pressure outlet) through the same
//!    projection.
//! 2. Rectangular-duct Poiseuille vs the exact Fourier series (n ≤ 99, with
//!    truncation-convergence check): TRT L∞rel ≤ 1e-3, flow rate ±0.5%.
//! 3. Sphere drag at Re ∈ {20, 100} vs Schiller–Naumann ±10% (D ≥ 24,
//!    blockage ≤ 3%, momentum-exchange probe). Heavy → #[ignore]; a light
//!    D=12 variant runs by default with a wider (resolution-limited) band.
//! 4. True 3D TGV, low-Re: short-run (t = 0.1/(νk²)) decay rate matches the
//!    diffusion limit and the field error converges at order ≥ 1.7 for
//!    N = 32→64 under diffusive scaling (u0 ∝ 1/N).

use lbm_core::prelude::*;
use std::f64::consts::PI;

type S2 = Solver<D2Q9, f64, CpuScalar, LocalPeriodic>;
type S3 = Solver<D3Q19, f64, CpuScalar, LocalPeriodic>;

fn max_abs_diff(a: &[f64], b: &[f64]) -> f64 {
    assert_eq!(a.len(), b.len());
    a.iter()
        .zip(b)
        .map(|(x, y)| (x - y).abs())
        .fold(0.0, f64::max)
}

/// Run until `max |u - u_prev| <= tol * max |u|` between checks (all three
/// velocity components), or `max_steps`. Returns true when steady.
fn run_to_steady3(s: &mut S3, check_every: usize, tol: f64, max_steps: usize) -> bool {
    let mut prev: Option<Vec<f64>> = None;
    let mut elapsed = 0;
    while elapsed < max_steps {
        s.run(check_every);
        elapsed += check_every;
        let mut cur = s.gather_ux();
        cur.extend(s.gather_uy());
        cur.extend(s.gather_uz());
        if let Some(p) = &prev {
            let dmax = max_abs_diff(&cur, p);
            let umax = cur.iter().fold(0.0f64, |m, v| m.max(v.abs()));
            if umax > 0.0 && dmax <= tol * umax {
                return true;
            }
        }
        prev = Some(cur);
    }
    false
}

// ===========================================================================
// 1. z-invariant degeneracy: D3Q19 (N×N×4, z periodic) vs D2Q9
// ===========================================================================

/// Slice a 3D compact field (z*(nx*ny) + y*nx + x) into its z = `z` plane.
fn zslice(f: &[f64], nx: usize, ny: usize, z: usize) -> &[f64] {
    &f[z * nx * ny..(z + 1) * nx * ny]
}

#[test]
fn t15_1_tgv_z_invariant_degenerates_to_d2q9() {
    let n = 32;
    let nz = 4;
    let nu = 0.02;
    let u0 = 1.28 / n as f64;
    let k = 2.0 * PI / n as f64;
    // Pressure-consistent 2D TGV initial state (docs/PHYSICS.md), z-invariant.
    let init = move |x: usize, y: usize| {
        let (xf, yf) = (k * x as f64, k * y as f64);
        let rho = 1.0 - 3.0 * u0 * u0 / 4.0 * ((2.0 * xf).cos() + (2.0 * yf).cos());
        (rho, [-u0 * xf.cos() * yf.sin(), u0 * xf.sin() * yf.cos(), 0.0])
    };
    let steps = (1.0 / (2.0 * nu * k * k)).round() as usize; // T1's t*

    let spec2 = GlobalSpec::<f64> {
        dims: [n, n, 1],
        nu,
        periodic: [true, true, false],
        ..Default::default()
    };
    let mut s2: S2 = Solver::new(&spec2, &[], &[], [1, 1, 1], CpuScalar::default(), LocalPeriodic);
    s2.init_with(move |x, y, _| init(x, y));
    s2.run(steps);

    let spec3 = GlobalSpec::<f64> {
        dims: [n, n, nz],
        nu,
        periodic: [true, true, true],
        ..Default::default()
    };
    let mut s3: S3 = Solver::new(&spec3, &[], &[], [1, 1, 1], CpuScalar::default(), LocalPeriodic);
    s3.init_with(move |x, y, _| init(x, y));
    s3.run(steps);

    let (r2, x2, y2) = (s2.gather_rho(), s2.gather_ux(), s2.gather_uy());
    let (r3, x3, y3, z3) = (s3.gather_rho(), s3.gather_ux(), s3.gather_uy(), s3.gather_uz());
    let mut dmax = 0.0f64;
    for z in 0..nz {
        dmax = dmax.max(max_abs_diff(zslice(&r3, n, n, z), &r2));
        dmax = dmax.max(max_abs_diff(zslice(&x3, n, n, z), &x2));
        dmax = dmax.max(max_abs_diff(zslice(&y3, n, n, z), &y2));
    }
    let uzmax = z3.iter().fold(0.0f64, |m, v| m.max(v.abs()));
    println!("TGV degeneracy over {steps} steps: max field diff {dmax:.3e}, max |uz| {uzmax:.3e}");
    assert!(dmax <= 1e-12, "D3Q19 z-invariant TGV drifted from D2Q9: {dmax:.3e}");
    assert!(uzmax <= 1e-13, "z-invariant flow grew uz = {uzmax:.3e}");
}

/// Same projection argument, but through the open-face BCs: a Zou–He channel
/// (parabolic velocity inlet + pressure outlet + y walls), z-periodic thin
/// slab vs the identical D2Q9 run. Pins the 3D face closure to the proven 2D
/// one (kernels.rs derivation note).
#[test]
fn t15_1b_zou_he_channel_degenerates_to_d2q9() {
    let (nx, ny, nz) = (64, 34, 4);
    let nu = 0.02;
    let umax = 0.05;
    let h = (ny - 2) as f64;
    let steps = 500;
    let parabola = move |y: usize| -> [f64; 3] {
        if y == 0 || y as f64 >= h + 1.0 {
            return [0.0; 3];
        }
        let yw = y as f64 - 0.5;
        [4.0 * umax * yw * (h - yw) / (h * h), 0.0, 0.0]
    };
    let mut walls = WallSpec::<f64>::default();
    walls.is_wall[Face::YNeg.index()] = true;
    walls.is_wall[Face::YPos.index()] = true;
    let mut faces = [FaceBC::Closed; 6];
    faces[Face::XNeg.index()] = FaceBC::Velocity { u: [0.0; 3] };
    faces[Face::XPos.index()] = FaceBC::Pressure { rho: 1.0 };

    let spec2 = GlobalSpec::<f64> {
        dims: [nx, ny, 1],
        nu,
        periodic: [false, false, false],
        faces,
        ..Default::default()
    };
    let (solid2, wu2) = build_wall_rims(2, spec2.dims, &walls);
    let mut s2: S2 = Solver::new(&spec2, &solid2, &wu2, [1, 1, 1], CpuScalar::default(), LocalPeriodic);
    s2.set_inlet_profile_with(Face::XNeg, |y, _| parabola(y));
    s2.run(steps);

    let spec3 = GlobalSpec::<f64> {
        dims: [nx, ny, nz],
        nu,
        periodic: [false, false, true],
        faces,
        ..Default::default()
    };
    let (solid3, wu3) = build_wall_rims(3, spec3.dims, &walls);
    let mut s3: S3 = Solver::new(&spec3, &solid3, &wu3, [1, 1, 1], CpuScalar::default(), LocalPeriodic);
    s3.set_inlet_profile_with(Face::XNeg, |y, _| parabola(y));
    s3.run(steps);

    let (r2, x2, y2) = (s2.gather_rho(), s2.gather_ux(), s2.gather_uy());
    let (r3, x3, y3, z3) = (s3.gather_rho(), s3.gather_ux(), s3.gather_uy(), s3.gather_uz());
    let mut dmax = 0.0f64;
    for z in 0..nz {
        dmax = dmax.max(max_abs_diff(zslice(&r3, nx, ny, z), &r2));
        dmax = dmax.max(max_abs_diff(zslice(&x3, nx, ny, z), &x2));
        dmax = dmax.max(max_abs_diff(zslice(&y3, nx, ny, z), &y2));
    }
    let uzmax = z3.iter().fold(0.0f64, |m, v| m.max(v.abs()));
    println!("Zou-He degeneracy over {steps} steps: max diff {dmax:.3e}, max |uz| {uzmax:.3e}");
    assert!(dmax <= 1e-12, "3D Zou-He channel drifted from D2Q9: {dmax:.3e}");
    assert!(uzmax <= 1e-13, "channel grew uz = {uzmax:.3e}");
}

/// The Zou–He reconstruction must satisfy all four moment constraints on the
/// face nodes exactly: after a step, u == u_bc on a velocity face (all three
/// components — the transverse corrections are what make uy/uz exact) and
/// rho == rho_bc on a pressure face, to rounding.
#[test]
fn t15_1c_zou_he_3d_enforces_prescribed_moments() {
    let (nx, ny, nz) = (10, 8, 6);
    let u_bc = [0.05, 0.01, -0.02];
    let mut faces = [FaceBC::Closed; 6];
    faces[Face::XNeg.index()] = FaceBC::Velocity { u: u_bc };
    faces[Face::XPos.index()] = FaceBC::Pressure { rho: 1.0 };
    let spec = GlobalSpec::<f64> {
        dims: [nx, ny, nz],
        nu: 0.05,
        periodic: [false, true, true],
        faces,
        ..Default::default()
    };
    let mut s: S3 = Solver::new(&spec, &[], &[], [1, 1, 1], CpuScalar::default(), LocalPeriodic);
    s.run(3);
    let mut du = 0.0f64;
    let mut drho = 0.0f64;
    for z in 0..nz {
        for y in 0..ny {
            let u = s.u(0, y, z);
            for a in 0..3 {
                du = du.max((u[a] - u_bc[a]).abs());
            }
            drho = drho.max((s.rho(nx - 1, y, z) - 1.0).abs());
        }
    }
    println!("Zou-He 3D moment enforcement: max |u - u_bc| = {du:.3e}, max |rho - rho_bc| = {drho:.3e}");
    assert!(du <= 1e-14, "velocity face violates prescribed u by {du:.3e}");
    assert!(drho <= 1e-14, "pressure face violates prescribed rho by {drho:.3e}");
}

// ===========================================================================
// 2. Rectangular duct Poiseuille vs exact series (T15.2)
// ===========================================================================

/// Exact series for the axial velocity in a rectangular duct
/// (−a ≤ ỹ ≤ a, −b ≤ z̃ ≤ b, no-slip, uniform body force g along x):
///
///   u(y, z̃) = (16 a² g)/(ν π³) Σ_{n odd} (1/n³)
///             [1 − cosh(nπz̃/2a)/cosh(nπb/2a)] sin(nπy/2a),   y = ỹ + a
///
/// (VALIDATION.md T15.2 form; satisfies ∇²u = −g/ν term-by-term since
/// Σ_{n odd} sin(nθ)/n = π/4 on (0, π).) Truncated at `nmax` (odd).
fn duct_series(y: f64, zt: f64, a: f64, b: f64, g: f64, nu: f64, nmax: usize) -> f64 {
    let pref = 16.0 * a * a * g / (nu * PI * PI * PI);
    let mut sum = 0.0;
    let mut n = 1;
    while n <= nmax {
        let nf = n as f64;
        let kn = nf * PI / (2.0 * a);
        // cosh ratio computed as exp differences to avoid overflow for large n.
        let ratio = ((kn * zt.abs()).exp() + (-kn * zt.abs()).exp())
            / ((kn * b).exp() + (-kn * b).exp());
        sum += (1.0 - ratio) * (kn * y).sin() / (nf * nf * nf);
        n += 2;
    }
    pref * sum
}

/// Exact volumetric flow rate of the truncated series:
/// Q = (64 a³ g)/(ν π⁴) Σ_{n odd} (1/n⁴) [2b − (4a/nπ) tanh(nπb/2a)].
fn duct_series_q(a: f64, b: f64, g: f64, nu: f64, nmax: usize) -> f64 {
    let pref = 64.0 * a * a * a * g / (nu * PI.powi(4));
    let mut sum = 0.0;
    let mut n = 1;
    while n <= nmax {
        let nf = n as f64;
        sum += (2.0 * b - 4.0 * a / (nf * PI) * (nf * PI * b / (2.0 * a)).tanh())
            / (nf * nf * nf * nf);
        n += 2;
    }
    pref * sum
}

#[test]
fn t15_2_rectangular_duct_poiseuille_matches_series() {
    // Cross-section (ny-2) × (nz-2) fluid cells, walls half-way behind the
    // rims; x periodic (the flow is x-invariant), driven by a uniform force.
    let (nx, ny, nz) = (4, 34, 34);
    let (hy, hz) = ((ny - 2) as f64, (nz - 2) as f64);
    let (a, b) = (hy / 2.0, hz / 2.0);
    let nu = 0.1;
    let g = 5e-6;
    let mut walls = WallSpec::<f64>::default();
    for f in [Face::YNeg, Face::YPos, Face::ZNeg, Face::ZPos] {
        walls.is_wall[f.index()] = true;
    }
    let spec = GlobalSpec::<f64> {
        dims: [nx, ny, nz],
        nu,
        periodic: [true, false, false],
        force: [g, 0.0, 0.0],
        ..Default::default()
    };
    let (solid, wall_u) = build_wall_rims(3, spec.dims, &walls);
    let mut s: S3 = Solver::new(&spec, &solid, &wall_u, [1, 1, 1], CpuScalar::default(), LocalPeriodic);
    let steady = run_to_steady3(&mut s, 500, 1e-11, 100_000);
    assert!(steady, "duct flow did not reach steady state");

    // Series truncation convergence: n ≤ 99 vs n ≤ 399 must agree far below
    // the acceptance tolerance (the tail is O(1/N²), ~2.5e-5 relative).
    let umax_ref = duct_series(a, 0.0, a, b, g, nu, 99);
    let mut trunc = 0.0f64;
    for j in 1..=(ny - 2) {
        for kz in 1..=(nz - 2) {
            let y = j as f64 - 0.5;
            let zt = kz as f64 - 0.5 - b;
            trunc = trunc
                .max((duct_series(y, zt, a, b, g, nu, 99) - duct_series(y, zt, a, b, g, nu, 399)).abs());
        }
    }
    let trunc_rel = trunc / umax_ref;
    assert!(
        trunc_rel < 1e-4,
        "series not converged at n<=99: {trunc_rel:.2e} of umax"
    );

    // L∞rel over the fluid cross-section (x = 0 slab; the field is
    // x-invariant), normalised by the analytic maximum.
    let ux = s.gather_ux();
    let uy = s.gather_uy();
    let uz = s.gather_uz();
    let idx = |x: usize, y: usize, z: usize| (z * ny + y) * nx + x;
    let mut err = 0.0f64;
    let mut q_meas = 0.0f64;
    let mut cross = 0.0f64;
    for kz in 1..=(nz - 2) {
        for j in 1..=(ny - 2) {
            let y = j as f64 - 0.5;
            let zt = kz as f64 - 0.5 - b;
            let ana = duct_series(y, zt, a, b, g, nu, 99);
            // Average the (identical to rounding) x slabs.
            let mut u = 0.0;
            for x in 0..nx {
                u += ux[idx(x, j, kz)];
            }
            u /= nx as f64;
            err = err.max((u - ana).abs());
            q_meas += u;
            cross = cross.max(uy[idx(0, j, kz)].abs()).max(uz[idx(0, j, kz)].abs());
        }
    }
    let linf_rel = err / umax_ref;
    let q_ana = duct_series_q(a, b, g, nu, 199);
    let q_rel = (q_meas - q_ana).abs() / q_ana;
    println!(
        "duct {}x{}: L∞rel = {linf_rel:.3e}, Q rel err = {q_rel:.3e} (meas {q_meas:.6e}, ana {q_ana:.6e}), max cross-flow {cross:.2e}, umax {umax_ref:.3e}",
        ny - 2,
        nz - 2
    );
    assert!(linf_rel <= 1e-3, "duct L∞rel = {linf_rel:.3e} > 1e-3");
    assert!(q_rel <= 5e-3, "duct flow rate off by {q_rel:.3e} (> 0.5%)");
    // Sanity, not in the spec: the numerical secondary flow (a staircase-
    // corner artifact) must stay far below the axial flow (measured 3.3e-6
    // of umax).
    assert!(
        cross < 1e-4 * umax_ref,
        "secondary flow {cross:.2e} is not small vs umax {umax_ref:.2e}"
    );
}

// ===========================================================================
// 3. Sphere drag vs Schiller–Naumann (T15.3)
// ===========================================================================

fn schiller_naumann(re: f64) -> f64 {
    24.0 / re * (1.0 + 0.15 * re.powf(0.687))
}

/// Schiller–Naumann evaluated at the hydrodynamic Reynolds number
/// Re_h = Re · (D+1)/D, pairing with the Cd normalisation by r_h = r + 0.5
/// (half-way bounce-back walls sit half a link outside the solid cells —
/// Ladd's staircase-sphere calibration). PM triage 2026-07-05: with the
/// (r_h, Re_h) pair the three sphere cases measure +0.6% / +7.1% / +2.3%
/// vs SN; with nominal (r, Re) the half-link bias (~ +2/D) consumed the
/// band at D = 24 (TESTING_NOTES #2 disposition).
fn sn_hydro(re: f64, d: usize) -> f64 {
    schiller_naumann(re * (d as f64 + 1.0) / d as f64)
}

struct SphereCase {
    d: usize,
    re: f64,
    u_in: f64,
    dims: [usize; 3],
}

/// Uniform flow past a staircase sphere: velocity inlet (XNeg), pressure
/// outlet (XPos), periodic lateral faces. Returns Cd from the
/// momentum-exchange probe once Cd stops changing (5e-4 relative between
/// 500-step windows).
///
/// Cd is the *mean* over each window, not the last-step sample: the
/// impulsive start rings a weakly damped inlet↔outlet standing acoustic
/// wave (decay time ~ 1/(νk²) ≈ 6e4 steps in the 192-long box) whose force
/// ripple is O(Ma) relative to the u² drag scale. At Re = 20 (u = 0.05)
/// instantaneous samples stall the convergence check for ~1e5 steps; the
/// window mean cancels the ripple and converges with the wake itself.
fn sphere_drag(c: &SphereCase) -> f64 {
    let [nx, ny, nz] = c.dims;
    let r = c.d as f64 / 2.0;
    let nu = c.u_in * c.d as f64 / c.re;
    let blockage = PI * r * r / ((ny * nz) as f64);
    assert!(blockage <= 0.03 + 1e-12, "blockage {blockage:.3} > 3%");
    let mut faces = [FaceBC::Closed; 6];
    faces[Face::XNeg.index()] = FaceBC::Velocity {
        u: [c.u_in, 0.0, 0.0],
    };
    faces[Face::XPos.index()] = FaceBC::Pressure { rho: 1.0 };
    let spec = GlobalSpec::<f64> {
        dims: c.dims,
        nu,
        periodic: [false, true, true],
        faces,
        ..Default::default()
    };
    let mut s: S3 = Solver::new(&spec, &[], &[], [1, 1, 1], CpuScalar::default(), LocalPeriodic);
    // Sphere centred 2.5 D downstream, laterally centred on the (odd)
    // half-integer grid centre so the staircase mask is mirror-symmetric.
    let cx = 2.5 * c.d as f64;
    let cy = (ny as f64 - 1.0) / 2.0;
    let cz = (nz as f64 - 1.0) / 2.0;
    let r2 = r * r;
    let inside = move |x: usize, y: usize, z: usize| {
        let (dx, dy, dz) = (x as f64 - cx, y as f64 - cy, z as f64 - cz);
        dx * dx + dy * dy + dz * dz <= r2
    };
    for z in 0..nz {
        for y in 0..ny {
            for x in 0..nx {
                if inside(x, y, z) {
                    s.set_solid(x, y, z);
                }
            }
        }
    }
    s.set_force_probe(inside);
    // Start from the free-stream state so the impulsive transient is short.
    s.init_with(move |x, y, z| {
        let u = if inside(x, y, z) { 0.0 } else { c.u_in };
        (1.0, [u, 0.0, 0.0])
    });
    // Cd is normalised with the hydrodynamic radius r_h = r + 0.5: half-way
    // bounce-back places the wall half a link outside the solid cells (Ladd's
    // staircase-sphere calibration). Measured 2026-07-05: with the nominal r
    // the D=24 runs read +13.2%/+7.2% vs Schiller-Naumann; with r_h all three
    // cases collapse to +0.6%/+7.1%/+2.3% (PM triage, TESTING_NOTES #2).
    let r_h = r + 0.5;
    let cd_of = |fx: f64| fx / (0.5 * c.u_in * c.u_in * PI * r_h * r_h);
    let mut cd_prev = f64::NAN;
    let mut cd = f64::NAN;
    for chunk in 0..200 {
        let mut fx_sum = 0.0;
        for _ in 0..500 {
            s.step();
            fx_sum += s.probed_force()[0];
        }
        cd = cd_of(fx_sum / 500.0);
        assert!(cd.is_finite(), "diverged at chunk {chunk}");
        if (cd - cd_prev).abs() <= 5e-4 * cd.abs() {
            break;
        }
        cd_prev = cd;
    }
    println!(
        "sphere D={} Re={} dims={:?} nu={nu:.4}: Cd = {cd:.4} (SN {:.4}) after {} steps",
        c.d,
        c.re,
        c.dims,
        schiller_naumann(c.re),
        s.time()
    );
    cd
}

/// Light default-suite variant: D = 12 (below the T15 D ≥ 24 resolution
/// floor, so the band is widened to ±25% — the spec-grade runs are the
/// #[ignore] tests below). Blockage π(D/2)²/(64·64) = 2.8% ≤ 3%, length 8D.
#[test]
fn t15_3_sphere_drag_re20_light() {
    let cd = sphere_drag(&SphereCase {
        d: 12,
        re: 20.0,
        u_in: 0.06,
        dims: [96, 64, 64],
    });
    let sn = sn_hydro(20.0, 12);
    let rel = (cd - sn).abs() / sn;
    assert!(rel <= 0.15, "Cd = {cd:.3} vs SN_h {sn:.3}: {rel:.3} > 15%");
}

/// T15.3 spec-grade: Re = 20, D = 24, blockage 2.8%, domain 8D.
/// Schiller–Naumann Cd = 2.6095 (the formula value; the spec's old "≈2.09"
/// parenthetical was a slip — TESTING_NOTES.md 2026-07-05), acceptance ±10%.
///
/// TRIAGED 2026-07-05: Cd and the SN reference both use the hydrodynamic
/// pair (r_h = r+0.5, Re_h = Re(D+1)/D) — see sn_hydro. Measured +7.1%
/// within the ±10% band. VALIDATION.md T15.3 updated accordingly.
#[test]
#[ignore = "heavy: ~3M cells, run with --include-ignored"]
fn t15_3_sphere_drag_re20() {
    let cd = sphere_drag(&SphereCase {
        d: 24,
        re: 20.0,
        u_in: 0.05,
        dims: [192, 128, 128],
    });
    let sn = sn_hydro(20.0, 24); // SN(20.833) = 2.544
    let rel = (cd - sn).abs() / sn;
    assert!(rel <= 0.10, "Cd = {cd:.4} vs SN_h {sn:.4}: {rel:.3} > 10%");
}

/// T15.3 spec-grade: Re = 100, D = 24 (steady axisymmetric wake regime).
/// Schiller–Naumann Cd ≈ 1.09, acceptance ±10%.
#[test]
#[ignore = "heavy: ~3M cells, run with --include-ignored"]
fn t15_3_sphere_drag_re100() {
    let cd = sphere_drag(&SphereCase {
        d: 24,
        re: 100.0,
        u_in: 0.1,
        dims: [192, 128, 128],
    });
    let sn = sn_hydro(100.0, 24); // SN(104.17) = 1.071
    let rel = (cd - sn).abs() / sn;
    assert!(rel <= 0.10, "Cd = {cd:.4} vs SN_h {sn:.4}: {rel:.3} > 10%");
}

// ===========================================================================
// 4. True 3D TGV: diffusion-limit decay rate + convergence order (T15.4)
// ===========================================================================

/// Run the classic 3D TGV for t* = 0.1/(νk²) steps under diffusive scaling
/// and return (L2rel error vs the diffusion-limit solution, relative decay-
/// rate error). Short run ⇒ vortex stretching stays negligible and the
/// linearised solution u = u_init e^{−3νk²t} is the reference.
fn tgv3d_short(n: usize, nu: f64, u0_coef: f64) -> (f64, f64) {
    let u0 = u0_coef / n as f64;
    let k = 2.0 * PI / n as f64;
    let spec = GlobalSpec::<f64> {
        dims: [n, n, n],
        nu,
        periodic: [true, true, true],
        ..Default::default()
    };
    let mut s: S3 = Solver::new(&spec, &[], &[], [1, 1, 1], CpuScalar::default(), LocalPeriodic);
    let vel = move |x: usize, y: usize, z: usize| -> [f64; 3] {
        let (xf, yf, zf) = (k * x as f64, k * y as f64, k * z as f64);
        [
            u0 * xf.sin() * yf.cos() * zf.cos(),
            -u0 * xf.cos() * yf.sin() * zf.cos(),
            0.0,
        ]
    };
    s.init_with(move |x, y, z| {
        let (xf, yf, zf) = (k * x as f64, k * y as f64, k * z as f64);
        // Classic TGV pressure field (Taylor & Green 1937); rho = 1 + p/cs².
        let p = u0 * u0 / 16.0 * (((2.0 * xf).cos() + (2.0 * yf).cos()) * ((2.0 * zf).cos() + 2.0));
        (1.0 + 3.0 * p, vel(x, y, z))
    });
    let ke = |s: &S3| -> f64 {
        let (ux, uy, uz) = (s.gather_ux(), s.gather_uy(), s.gather_uz());
        ux.iter()
            .zip(&uy)
            .zip(&uz)
            .map(|((a, b), c)| a * a + b * b + c * c)
            .sum()
    };
    let tstar = (0.1 / (nu * k * k)).round() as usize;
    let e0 = ke(&s);
    s.run(tstar);
    let e1 = ke(&s);
    // Decay rate vs the diffusion limit 2ν|k|² with |k|² = 3k².
    let rate = -(e1 / e0).ln() / tstar as f64;
    let rate_ref = 6.0 * nu * k * k;
    let rate_rel = (rate - rate_ref).abs() / rate_ref;
    // Field error vs the diffusion-limit solution.
    let decay = (-3.0 * nu * k * k * tstar as f64).exp();
    let (ux, uy, uz) = (s.gather_ux(), s.gather_uy(), s.gather_uz());
    let mut num = 0.0;
    let mut den = 0.0;
    for z in 0..n {
        for y in 0..n {
            for x in 0..n {
                let i = (z * n + y) * n + x;
                let v = vel(x, y, z);
                num += (ux[i] - v[0] * decay).powi(2)
                    + (uy[i] - v[1] * decay).powi(2)
                    + (uz[i] - v[2] * decay).powi(2);
                den += (v[0] * decay).powi(2) + (v[1] * decay).powi(2);
            }
        }
    }
    let l2 = (num / den).sqrt();
    println!(
        "TGV3D N={n}: t*={tstar}, L2rel = {l2:.4e}, rate rel err = {rate_rel:.4e} (u0 = {u0:.2e})"
    );
    (l2, rate_rel)
}

#[test]
fn t15_4_tgv3d_diffusive_convergence() {
    let nu = 0.02;
    // Low-Re constant: the classic 3D TGV is not an exact NS solution, so
    // its (nonlinear) deviation from the diffusion-limit reference is a
    // *resolution-independent relative* offset under diffusive scaling,
    // measured ∝ u0: L2rel ≈ 0.13 u0/(νk) (0.165 at coef 1.28, 0.0165 at
    // 0.128, …). The coefficient is chosen so that offset (~2e-5) sits far
    // below the spatial error floor (e32 = 1.29e-3, e64 = 3.44e-4), i.e.
    // the test measures spatial convergence, per the spec's low-Re /
    // weak-vortex-stretching intent. f64 keeps 6+ digits of headroom at
    // u0 = 2e-6.
    let u0_coef = 1.28e-4;
    let (e32, _) = tgv3d_short(32, nu, u0_coef);
    let (e64, rate64) = tgv3d_short(64, nu, u0_coef);
    let order = (e32 / e64).log2();
    println!("TGV3D convergence: e32 = {e32:.4e}, e64 = {e64:.4e}, order = {order:.3}");
    assert!(
        rate64 <= 0.02,
        "N=64 decay rate off the diffusion limit by {rate64:.3e} (> 2%)"
    );
    assert!(
        order >= 1.7,
        "order = {order:.3} (e32 = {e32:.3e}, e64 = {e64:.3e})"
    );
}
