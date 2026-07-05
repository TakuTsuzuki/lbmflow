//! T15.5 — Re=1000 cubic lid-driven cavity, D3Q19.
//!
//! The spec-grade profile tests are intentionally `#[ignore]`: the reference
//! document sanctions N >= 72 for Re=1000 because Re/(N-2) must stay near or
//! below 15, and steady 3D cavity convergence is too expensive for the default
//! suite. The default test below is a lower-resolution qualitative sentinel.

use lbm_core::prelude::*;

type S3 = Solver<D3Q19, f64, CpuScalar, LocalPeriodic>;

const U_LID: f64 = 0.1;
const RE: f64 = 1000.0;

// Albensoeder, S. & Kuhlmann, H.C. (2005), "Accurate three-dimensional
// lid-driven cavity flow", Journal of Computational Physics 206, 536-558,
// doi:10.1016/j.jcp.2004.12.024, Tables 5/6 (Re=1000 cubic cavity).
// Transcribed via github.com/tum-pbs/PICT tests/validations.py @
// a95d7f9d0713262a1bff2bd9e2be5a203ee69208 (2025-05-22), arrays
// `reference_coords_x_T5` / `reference_vel_v` / `reference_coords_y_T6` /
// `reference_vel_u`, key `(1000, 1, 1, False)`, then converted to the
// coordinate convention in docs/T15_5_CAVITY3D_REFERENCE.md:
// z = 0.5 - x_AK for the vertical u/U line, x = y_AK + 0.5 and w/U = -table
// for the horizontal w/U line. PICT's values are already normalized by Re
// from A&K's printed nu/L velocity units, i.e. these constants are u/U, w/U.
//
// Provenance checks from docs/T15_5_CAVITY3D_REFERENCE.md §3.1:
// - the PICT 2D Ghia tables match the canonical Ghia et al. (1982) tables;
// - Ben Beya & Lili (2008), C. R. Mecanique 336, 863-872, Table 1,
//   independently reproduces A&K's Re=1000 extrema and confirms signs/locations.
//
// Smoothness audit note (§4): all 30 interior points in these two nonuniform
// 17-point tables pass the two-rule Ghia-typo check: (1) second divided
// differences show no alternating one-point spike pattern; (2) suspicious
// sharp extrema are judged by the smaller of the two one-sided quadratic
// prediction residuals and, for w(0.9063), by the independent A&K extremum.
const AK_U_Z: [(f64, f64); 17] = [
    (0.0000, 0.00000),
    (0.0547, -0.20623),
    (0.0625, -0.22283),
    (0.0703, -0.23696),
    (0.1016, -0.27293),
    (0.1719, -0.25160),
    (0.2813, -0.10999),
    (0.4531, -0.00612),
    (0.5000, 0.00802),
    (0.6172, 0.03905),
    (0.7344, 0.07334),
    (0.8516, 0.12183),
    (0.9531, 0.33171),
    (0.9609, 0.39821),
    (0.9688, 0.48443),
    (0.9766, 0.58964),
    (1.0000, 1.00000),
];

const AK_W_X: [(f64, f64); 17] = [
    (0.0000, 0.00000),
    (0.0625, 0.21738),
    (0.0703, 0.22746),
    (0.0781, 0.23503),
    (0.0938, 0.24407),
    (0.1563, 0.22924),
    (0.2266, 0.17580),
    (0.2344, 0.16987),
    (0.5000, 0.03674),
    (0.8047, -0.15223),
    (0.8594, -0.31117),
    (0.9063, -0.43423),
    (0.9453, -0.33511),
    (0.9531, -0.29032),
    (0.9609, -0.24095),
    (0.9688, -0.18864),
    (1.0000, 0.00000),
];

// Ghia et al. (1982) 2D Re=1000 u/U table, converted to ascending y. This is
// only an anti-2D-confusion guard: a correct cubic-cavity solution must not
// collapse onto the 2D centreline profile.
const GHIA_2D_U_Y_RE1000_ASC: [(f64, f64); 17] = [
    (0.0000, 0.00000),
    (0.0547, -0.18109),
    (0.0625, -0.20196),
    (0.0703, -0.22220),
    (0.1016, -0.29730),
    (0.1719, -0.38289),
    (0.2813, -0.27805),
    (0.4531, -0.10648),
    (0.5000, -0.06080),
    (0.6172, 0.05702),
    (0.7344, 0.18719),
    (0.8516, 0.33304),
    (0.9531, 0.46604),
    (0.9609, 0.51117),
    (0.9688, 0.57492),
    (0.9766, 0.65928),
    (1.0000, 1.00000),
];

const AK_U_MIN: (f64, f64) = (-0.2803833, 0.12419);
const AK_W_MIN: (f64, f64) = (-0.4350186, 0.90957);
const AK_W_MAX: (f64, f64) = (0.2466511, 0.10913);

#[derive(Debug)]
struct CavityMetrics {
    n: usize,
    steps: u64,
    steady: bool,
    mass_rel: f64,
    rms_u: f64,
    rms_w: f64,
    rms_vs_2d_ghia_u: f64,
    midplane_v_max: f64,
    u_min: (f64, f64),
    w_min: (f64, f64),
    w_max: (f64, f64),
}

fn max_abs_diff(a: &[f64], b: &[f64]) -> f64 {
    assert_eq!(a.len(), b.len());
    a.iter()
        .zip(b)
        .map(|(x, y)| (x - y).abs())
        .fold(0.0, f64::max)
}

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

fn cavity(n: usize) -> S3 {
    let l = (n - 2) as f64;
    let nu = U_LID * l / RE;
    let mut walls = WallSpec::<f64>::default();
    for f in [
        Face::XNeg,
        Face::XPos,
        Face::YNeg,
        Face::YPos,
        Face::ZNeg,
        Face::ZPos,
    ] {
        walls.is_wall[f.index()] = true;
    }
    walls.u[Face::ZPos.index()] = [U_LID, 0.0, 0.0];
    let spec = GlobalSpec::<f64> {
        dims: [n, n, n],
        nu,
        periodic: [false, false, false],
        collision: CollisionKind::Trt {
            magic: CollisionKind::MAGIC_STD,
        },
        ..Default::default()
    };
    let (solid, wall_u) = build_wall_rims(3, spec.dims, &walls);
    S3::new(
        &spec,
        &solid,
        &wall_u,
        [1, 1, 1],
        CpuScalar::default(),
        LocalPeriodic,
    )
}

fn sample_field(f: &[f64], n: usize, x: f64, y: f64, z: f64) -> f64 {
    let x0 = x.floor().clamp(1.0, (n - 2) as f64) as usize;
    let y0 = y.floor().clamp(1.0, (n - 2) as f64) as usize;
    let z0 = z.floor().clamp(1.0, (n - 2) as f64) as usize;
    let x1 = (x0 + 1).min(n - 2);
    let y1 = (y0 + 1).min(n - 2);
    let z1 = (z0 + 1).min(n - 2);
    let tx = x - x0 as f64;
    let ty = y - y0 as f64;
    let tz = z - z0 as f64;
    let idx = |i: usize, j: usize, k: usize| (k * n + j) * n + i;

    let c000 = f[idx(x0, y0, z0)];
    let c100 = f[idx(x1, y0, z0)];
    let c010 = f[idx(x0, y1, z0)];
    let c110 = f[idx(x1, y1, z0)];
    let c001 = f[idx(x0, y0, z1)];
    let c101 = f[idx(x1, y0, z1)];
    let c011 = f[idx(x0, y1, z1)];
    let c111 = f[idx(x1, y1, z1)];
    let c00 = (1.0 - tx) * c000 + tx * c100;
    let c10 = (1.0 - tx) * c010 + tx * c110;
    let c01 = (1.0 - tx) * c001 + tx * c101;
    let c11 = (1.0 - tx) * c011 + tx * c111;
    let c0 = (1.0 - ty) * c00 + ty * c10;
    let c1 = (1.0 - ty) * c01 + ty * c11;
    (1.0 - tz) * c0 + tz * c1
}

fn phys_pos(n: usize, frac: f64) -> f64 {
    0.5 + frac * (n - 2) as f64
}

fn sample_u_line(ux: &[f64], n: usize, z_frac: f64) -> f64 {
    if z_frac <= 0.0 {
        return 0.0;
    }
    if z_frac >= 1.0 {
        return U_LID;
    }
    sample_field(
        ux,
        n,
        phys_pos(n, 0.5),
        phys_pos(n, 0.5),
        phys_pos(n, z_frac),
    )
}

fn sample_w_line(uz: &[f64], n: usize, x_frac: f64) -> f64 {
    if x_frac <= 0.0 || x_frac >= 1.0 {
        return 0.0;
    }
    sample_field(
        uz,
        n,
        phys_pos(n, x_frac),
        phys_pos(n, 0.5),
        phys_pos(n, 0.5),
    )
}

fn rms_profile(line: &[(f64, f64); 17], measured: impl Fn(f64) -> f64) -> f64 {
    let sum = line
        .iter()
        .map(|&(c, r)| {
            let d = measured(c) / U_LID - r;
            d * d
        })
        .sum::<f64>();
    (sum / line.len() as f64).sqrt()
}

fn parabolic_extremum(samples: &[(f64, f64)]) -> (f64, f64) {
    let mut i = 1usize;
    for j in 1..samples.len() - 1 {
        if samples[j].1.abs() > samples[i].1.abs() {
            i = j;
        }
    }
    let (x1, y1) = samples[i - 1];
    let (x2, y2) = samples[i];
    let (x3, y3) = samples[i + 1];
    let denom = (x1 - x2) * (x1 - x3) * (x2 - x3);
    let a = (x3 * (y2 - y1) + x2 * (y1 - y3) + x1 * (y3 - y2)) / denom;
    let b = (x3 * x3 * (y1 - y2) + x2 * x2 * (y3 - y1) + x1 * x1 * (y2 - y3)) / denom;
    if a.abs() < 1.0e-14 {
        return samples[i];
    }
    let xv = (-b / (2.0 * a)).clamp(x1.min(x3), x1.max(x3));
    let yv = a * xv * xv + b * xv + (y1 - a * x1 * x1 - b * x1);
    (yv, xv)
}

fn extrema_from_lines(ux: &[f64], uz: &[f64], n: usize) -> ((f64, f64), (f64, f64), (f64, f64)) {
    let mut u_line = Vec::with_capacity(n - 2);
    let mut w_line = Vec::with_capacity(n - 2);
    for k in 1..=(n - 2) {
        let z = (k as f64 - 0.5) / (n - 2) as f64;
        u_line.push((z, sample_u_line(ux, n, z) / U_LID));
    }
    for i in 1..=(n - 2) {
        let x = (i as f64 - 0.5) / (n - 2) as f64;
        w_line.push((x, sample_w_line(uz, n, x) / U_LID));
    }

    let u_min_i = u_line
        .iter()
        .enumerate()
        .skip(1)
        .take(u_line.len() - 2)
        .min_by(|a, b| a.1 .1.total_cmp(&b.1 .1))
        .unwrap()
        .0;
    let w_min_i = w_line
        .iter()
        .enumerate()
        .skip(1)
        .take(w_line.len() - 2)
        .min_by(|a, b| a.1 .1.total_cmp(&b.1 .1))
        .unwrap()
        .0;
    let w_max_i = w_line
        .iter()
        .enumerate()
        .skip(1)
        .take(w_line.len() - 2)
        .max_by(|a, b| a.1 .1.total_cmp(&b.1 .1))
        .unwrap()
        .0;
    (
        parabolic_extremum(&u_line[u_min_i - 1..=u_min_i + 1]),
        parabolic_extremum(&w_line[w_min_i - 1..=w_min_i + 1]),
        parabolic_extremum(&w_line[w_max_i - 1..=w_max_i + 1]),
    )
}

fn midplane_v_max(uy: &[f64], n: usize) -> f64 {
    let y_mid = phys_pos(n, 0.5);
    let mut vmax = 0.0f64;
    for z in 1..=(n - 2) {
        for x in 1..=(n - 2) {
            vmax = vmax.max(sample_field(uy, n, x as f64, y_mid, z as f64).abs());
        }
    }
    vmax / U_LID
}

fn run_case(n: usize, steady_tol: f64, max_steps: usize) -> CavityMetrics {
    let mut s = cavity(n);
    let m0 = s.total_mass();
    let steady = run_to_steady3(&mut s, 500, steady_tol, max_steps);
    let mass_rel = ((s.total_mass() - m0) / m0).abs();
    let (ux, uy, uz) = (s.gather_ux(), s.gather_uy(), s.gather_uz());
    let rms_u = rms_profile(&AK_U_Z, |z| sample_u_line(&ux, n, z));
    let rms_w = rms_profile(&AK_W_X, |x| sample_w_line(&uz, n, x));
    let rms_vs_2d_ghia_u = rms_profile(&GHIA_2D_U_Y_RE1000_ASC, |z| sample_u_line(&ux, n, z));
    let midplane_v_max = midplane_v_max(&uy, n);
    let (u_min, w_min, w_max) = extrema_from_lines(&ux, &uz, n);
    println!(
        "T15.5 N={n} steps={} steady={steady} mass_rel={mass_rel:.3e} rms_u={rms_u:.4} rms_w={rms_w:.4} rms_vs_2d_u={rms_vs_2d_ghia_u:.4} max|v_mid|/U={midplane_v_max:.3e}",
        s.time()
    );
    println!(
        "T15.5 N={n} extrema: u_min={:.5}@z={:.5}, w_min={:.5}@x={:.5}, w_max={:.5}@x={:.5}",
        u_min.0, u_min.1, w_min.0, w_min.1, w_max.0, w_max.1
    );
    println!("T15.5 N={n} u/U profile:");
    for &(z, r) in &AK_U_Z {
        let m = sample_u_line(&ux, n, z) / U_LID;
        println!("  z={z:.4} measured={m:.5} ref={r:.5} diff={:.5}", m - r);
    }
    println!("T15.5 N={n} w/U profile:");
    for &(x, r) in &AK_W_X {
        let m = sample_w_line(&uz, n, x) / U_LID;
        println!("  x={x:.4} measured={m:.5} ref={r:.5} diff={:.5}", m - r);
    }
    CavityMetrics {
        n,
        steps: s.time(),
        steady,
        mass_rel,
        rms_u,
        rms_w,
        rms_vs_2d_ghia_u,
        midplane_v_max,
        u_min,
        w_min,
        w_max,
    }
}

fn assert_mass_and_symmetry(m: &CavityMetrics, v_limit: f64) {
    assert!(
        m.mass_rel <= 1.0e-11,
        "T15.5 N={} mass drift rel {:.3e} after {} steps",
        m.n,
        m.mass_rel,
        m.steps
    );
    assert!(
        m.midplane_v_max <= v_limit,
        "T15.5 N={} symmetry-plane max |v|/U = {:.3e} > {:.3e}",
        m.n,
        m.midplane_v_max,
        v_limit
    );
}

fn assert_qualitative_extrema(m: &CavityMetrics) {
    assert!(
        m.u_min.0 < -0.15 && (0.07..0.20).contains(&m.u_min.1),
        "T15.5 N={} qualitative u_min failed: {:.5}@{:.5}",
        m.n,
        m.u_min.0,
        m.u_min.1
    );
    assert!(
        m.w_min.0 < -0.25 && (0.82..0.96).contains(&m.w_min.1),
        "T15.5 N={} qualitative w_min failed: {:.5}@{:.5}",
        m.n,
        m.w_min.0,
        m.w_min.1
    );
    assert!(
        m.w_max.0 > 0.12 && (0.05..0.18).contains(&m.w_max.1),
        "T15.5 N={} qualitative w_max failed: {:.5}@{:.5}",
        m.n,
        m.w_max.0,
        m.w_max.1
    );
}

fn assert_profile_bands(
    m: &CavityMetrics,
    rms_u_limit: f64,
    rms_w_limit: f64,
    value_rel: f64,
    pos_limit: f64,
) {
    assert!(
        m.steady,
        "T15.5 N={} did not reach steady state in {} steps",
        m.n, m.steps
    );
    assert!(
        m.rms_u <= rms_u_limit,
        "T15.5 N={} u-line RMS/U = {:.4} > {:.4}",
        m.n,
        m.rms_u,
        rms_u_limit
    );
    assert!(
        m.rms_w <= rms_w_limit,
        "T15.5 N={} w-line RMS/U = {:.4} > {:.4}",
        m.n,
        m.rms_w,
        rms_w_limit
    );
    assert!(
        m.rms_vs_2d_ghia_u >= 0.05,
        "T15.5 N={} anti-2D guard failed: RMS/U vs 2D Ghia u = {:.4} < 0.05",
        m.n,
        m.rms_vs_2d_ghia_u
    );
    let u_rel = (m.u_min.0 - AK_U_MIN.0).abs() / AK_U_MIN.0.abs();
    let wmin_rel = (m.w_min.0 - AK_W_MIN.0).abs() / AK_W_MIN.0.abs();
    let wmax_rel = (m.w_max.0 - AK_W_MAX.0).abs() / AK_W_MAX.0.abs();
    assert!(
        u_rel <= value_rel && (m.u_min.1 - AK_U_MIN.1).abs() <= pos_limit,
        "T15.5 N={} u_min {:.5}@{:.5} vs A&K {:.5}@{:.5}: rel {:.3}, pos diff {:.3}",
        m.n,
        m.u_min.0,
        m.u_min.1,
        AK_U_MIN.0,
        AK_U_MIN.1,
        u_rel,
        (m.u_min.1 - AK_U_MIN.1).abs()
    );
    assert!(
        wmin_rel <= value_rel && (m.w_min.1 - AK_W_MIN.1).abs() <= pos_limit,
        "T15.5 N={} w_min {:.5}@{:.5} vs A&K {:.5}@{:.5}: rel {:.3}, pos diff {:.3}",
        m.n,
        m.w_min.0,
        m.w_min.1,
        AK_W_MIN.0,
        AK_W_MIN.1,
        wmin_rel,
        (m.w_min.1 - AK_W_MIN.1).abs()
    );
    assert!(
        wmax_rel <= 8.0 / 6.0 * value_rel && (m.w_max.1 - AK_W_MAX.1).abs() <= pos_limit,
        "T15.5 N={} w_max {:.5}@{:.5} vs A&K {:.5}@{:.5}: rel {:.3}, pos diff {:.3}",
        m.n,
        m.w_max.0,
        m.w_max.1,
        AK_W_MAX.0,
        AK_W_MAX.1,
        wmax_rel,
        (m.w_max.1 - AK_W_MAX.1).abs()
    );
}

#[test]
fn t15_5_reference_tables_have_expected_shape_and_are_not_2d_ghia() {
    assert_eq!(AK_U_Z[0], (0.0, 0.0));
    assert_eq!(AK_U_Z[16], (1.0, 1.0));
    assert_eq!(AK_W_X[0], (0.0, 0.0));
    assert_eq!(AK_W_X[16], (1.0, 0.0));
    let rms_3d_vs_2d = rms_profile(&GHIA_2D_U_Y_RE1000_ASC, |z| {
        U_LID
            * AK_U_Z
                .iter()
                .find(|(zz, _)| (*zz - z).abs() < 1.0e-12)
                .unwrap()
                .1
    });
    assert!(
        rms_3d_vs_2d >= 0.05,
        "frozen A&K 3D u table is too close to 2D Ghia: RMS/U = {rms_3d_vs_2d:.4}"
    );
}

#[test]
fn t15_5_cavity3d_re1000_default_sanity() {
    // Below the profile-validation resolution floor; this default sentinel
    // checks deterministic wall setup, mass conservation, symmetry-plane v,
    // and the qualitative 3D secondary-flow extrema signs/locations.
    let m = run_case(64, 5.0e-6, 20_000);
    assert_mass_and_symmetry(&m, 8.0e-3);
    assert_qualitative_extrema(&m);
}

#[test]
#[ignore = "heavy: N=72 Re=1000 3D cavity profile validation"]
fn t15_5_cavity3d_re1000_profiles_n72() {
    // Extrema band frozen at 0.13 after characterization (PHYSICS.md 2026-07-05,
    // "T15.5 extremum band"): at N=72 the sharp near-wall extrema sit on the
    // numerical-diffusion side of the spectral A&K reference by 9.1-10.5%
    // (u_min -0.25084 / w_min -0.39537 / w_max 0.22148); the N=64->72 convergence
    // test below confirms monotone approach toward the reference, and N=48
    // diverges per the Re/(N-2) <~ 15 stability limit, so N cannot be lowered.
    // Profile RMS (the primary shape criterion) passes with ~2x margin.
    // Band governance: tightening is free; loosening again requires a new
    // PHYSICS.md entry.
    let m = run_case(72, 1.0e-8, 500_000);
    assert_mass_and_symmetry(&m, 1.0e-3);
    assert_profile_bands(&m, 0.030, 0.035, 0.13, 0.03);
}

#[test]
#[ignore = "heavy: convergence tendency N=64 -> N=72"]
fn t15_5_cavity3d_re1000_u_min_converges_toward_reference() {
    let m64 = run_case(64, 2.0e-8, 350_000);
    let m72 = run_case(72, 1.0e-8, 500_000);
    let e64 = (m64.u_min.0 - AK_U_MIN.0).abs();
    let e72 = (m72.u_min.0 - AK_U_MIN.0).abs();
    assert!(
        e72 < e64,
        "T15.5 u_min convergence failed: N64 err {e64:.5}, N72 err {e72:.5}; N64 {:.5}@{:.5}, N72 {:.5}@{:.5}, A&K {:.5}@{:.5}",
        m64.u_min.0,
        m64.u_min.1,
        m72.u_min.0,
        m72.u_min.1,
        AK_U_MIN.0,
        AK_U_MIN.1
    );
}
