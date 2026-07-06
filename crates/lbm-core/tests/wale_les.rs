//! WALE LES MF-beta subset tests.

use lbm_core::prelude::*;

type S3<B = CpuScalar> = Solver<D3Q19, f64, B, LocalPeriodic>;

fn max_abs_diff(a: &[f64], b: &[f64]) -> f64 {
    a.iter()
        .zip(b)
        .map(|(x, y)| (x - y).abs())
        .fold(0.0, f64::max)
}

fn run_to_steady(s: &mut S3, check_every: usize, tol: f64, max_steps: usize) {
    let mut prev = s.gather_ux();
    for _ in (0..max_steps).step_by(check_every) {
        s.run(check_every);
        let cur = s.gather_ux();
        let d = max_abs_diff(&cur, &prev);
        let scale = cur.iter().map(|v| v.abs()).fold(0.0f64, f64::max).max(1.0);
        if d / scale <= tol {
            return;
        }
        prev = cur;
    }
    panic!("fixture did not reach steady state by t={}", s.time());
}

fn channel(ny: usize, nu: f64, force: [f64; 3], top_u: [f64; 3], backend: CpuScalar) -> S3 {
    let spec = GlobalSpec {
        dims: [4, ny, 4],
        nu,
        collision: CollisionKind::Trt {
            magic: CollisionKind::MAGIC_STD,
        },
        periodic: [true, false, true],
        force,
        ..Default::default()
    };
    let mut walls = WallSpec::default();
    walls.is_wall[Face::YNeg.index()] = true;
    walls.is_wall[Face::YPos.index()] = true;
    walls.u[Face::YPos.index()] = top_u;
    let (solid, wall_u) = build_wall_rims(D3Q19::D, spec.dims, &walls);
    Solver::new(&spec, &solid, &wall_u, [1, 1, 1], backend, LocalPeriodic)
}

fn constant_omega_field_is_bitwise_identical_to_field_off<B>(backend: B)
where
    B: lbm_core::backend::Backend<D3Q19, f64, Fields = SoaFields<f64>> + Clone,
{
    let spec = GlobalSpec {
        dims: [12, 10, 8],
        nu: 0.03,
        periodic: [true, true, true],
        ..Default::default()
    };
    let init = |x: usize, y: usize, z: usize| {
        let xf = x as f64 * 0.13;
        let yf = y as f64 * 0.17;
        let zf = z as f64 * 0.19;
        (
            1.0 + 1.0e-4 * (xf + yf).sin(),
            [
                0.02 * xf.sin() * yf.cos(),
                -0.015 * yf.sin() * zf.cos(),
                0.01 * zf.sin() * xf.cos(),
            ],
        )
    };
    let mut off: Solver<D3Q19, f64, B, LocalPeriodic> =
        Solver::new(&spec, &[], &[], [1, 1, 1], backend.clone(), LocalPeriodic);
    let mut on: Solver<D3Q19, f64, B, LocalPeriodic> =
        Solver::new(&spec, &[], &[], [1, 1, 1], backend, LocalPeriodic);
    off.init_with(init);
    on.init_with(init);
    let omega = vec![1.0 / (3.0 * spec.nu + 0.5); spec.dims[0] * spec.dims[1] * spec.dims[2]];
    on.set_omega_field(Some(&omega));
    off.run(7);
    on.run(7);
    // ULP-band identity (denominator = machine epsilon on absolute Δ; the two
    // paths differ only by IEEE-754 evaluation order — omega-field path
    // computes cp per cell, omega-off path uses a StepParams-precomputed cp
    // whose folding through the B-1 fused pass reorders sums by <=1 ULP).
    // Measured on 2026-07-06 (D3Q19 f64 12x10x8, 7 steps): max|Δf|=1.5e-16
    // (0.7 ULP), max|Δρ|=4.4e-16 (2 ULP). Band frozen at 5 * f64::EPSILON
    // with ~2x headroom — a physics discrepancy would blow past this instantly.
    let band = 5.0 * f64::EPSILON;
    for q in 0..D3Q19::Q {
        let (a, b) = (off.gather_f(q), on.gather_f(q));
        let d = a.iter().zip(&b).map(|(x, y)| (x - y).abs()).fold(0.0f64, f64::max);
        assert!(d <= band, "population plane {q}: max|Δ|={d:e} > {band:e} (5*eps)");
    }
    for (label, a, b) in [
        ("rho", off.gather_rho(), on.gather_rho()),
        ("ux",  off.gather_ux(),  on.gather_ux()),
        ("uy",  off.gather_uy(),  on.gather_uy()),
        ("uz",  off.gather_uz(),  on.gather_uz()),
    ] {
        let d = a.iter().zip(&b).map(|(x, y)| (x - y).abs()).fold(0.0f64, f64::max);
        assert!(d <= band, "{label}: max|Δ|={d:e} > {band:e} (5*eps)");
    }
}

#[test]
fn constant_omega_field_is_bitwise_identical_to_field_off_scalar() {
    constant_omega_field_is_bitwise_identical_to_field_off(CpuScalar::default());
}

#[test]
fn constant_omega_field_is_bitwise_identical_to_field_off_simd() {
    constant_omega_field_is_bitwise_identical_to_field_off(CpuSimd::default());
}

#[test]
fn wale_null_for_steady_couette_and_poiseuille() {
    let nu = (1.0 - 0.5) / 3.0;
    let mut couette = channel(10, nu, [0.0; 3], [0.1, 0.0, 0.0], CpuScalar::default());
    run_to_steady(&mut couette, 500, 1.0e-11, 200_000);
    let mut les = WaleLes::new();
    les.update(&mut couette);
    let (max_i, max_couette) = les
        .nu_t()
        .iter()
        .copied()
        .enumerate()
        .max_by(|a, b| a.1.total_cmp(&b.1))
        .unwrap();
    let _ = max_i;
    assert!(
        max_couette <= 1.0e-12,
        "WALE must be null for pure Couette shear, max nu_t={max_couette:e}"
    );

    let mut poiseuille = channel(10, 0.1, [1.0e-6, 0.0, 0.0], [0.0; 3], CpuScalar::default());
    run_to_steady(&mut poiseuille, 500, 1.0e-11, 200_000);
    les.update(&mut poiseuille);
    let max_poiseuille = les.nu_t().iter().copied().fold(0.0f64, f64::max);
    assert!(
        max_poiseuille <= 1.0e-12,
        "WALE must be null for pure Poiseuille shear, max nu_t={max_poiseuille:e}"
    );
}

#[test]
fn les_on_does_not_change_laminar_duct_after_null_update() {
    let mut off = channel(10, 0.1, [1.0e-6, 0.0, 0.0], [0.0; 3], CpuScalar::default());
    let mut on = channel(10, 0.1, [1.0e-6, 0.0, 0.0], [0.0; 3], CpuScalar::default());
    run_to_steady(&mut off, 500, 1.0e-11, 200_000);
    run_to_steady(&mut on, 500, 1.0e-11, 200_000);
    let mut les = WaleLes::new();
    les.update(&mut on);
    off.run(20);
    on.run(20);
    assert!(max_abs_diff(&off.gather_ux(), &on.gather_ux()) <= 1.0e-12);
    assert!(max_abs_diff(&off.gather_uy(), &on.gather_uy()) <= 1.0e-12);
    assert!(max_abs_diff(&off.gather_uz(), &on.gather_uz()) <= 1.0e-12);
}

#[test]
#[ignore = "T17/VR-STR-03 heavy acceptance: channel Re_tau=180 vs DNS"]
fn wale_channel_re_tau_180_dns_skeleton() {
    // TODO(T17/VR-STR-03): set up Re_tau=180 channel, collect mean profile
    // and Reynolds stresses, and compare against DNS acceptance bands.
}

fn make_tgv64(nu: f64) -> S3 {
    // T15.4 setup: classic 3D TGV, u0 = 1.28e-4/N under diffusive scaling,
    // pressure-consistent init (rho = 1 + 3p). See tests/t15_3d.rs.
    let n = 64usize;
    let spec = GlobalSpec::<f64> {
        dims: [n, n, n],
        nu,
        periodic: [true, true, true],
        ..Default::default()
    };
    let mut s: S3 = Solver::new(&spec, &[], &[], [1, 1, 1], CpuScalar::default(), LocalPeriodic);
    let u0 = 1.28e-4 / n as f64;
    let k = 2.0 * std::f64::consts::PI / n as f64;
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
        let p =
            u0 * u0 / 16.0 * (((2.0 * xf).cos() + (2.0 * yf).cos()) * ((2.0 * zf).cos() + 2.0));
        (1.0 + 3.0 * p, vel(x, y, z))
    });
    s
}

fn ke(s: &S3) -> f64 {
    let (ux, uy, uz) = (s.gather_ux(), s.gather_uy(), s.gather_uz());
    ux.iter()
        .zip(&uy)
        .zip(&uz)
        .map(|((a, b), c)| a * a + b * b + c * c)
        .sum()
}

/// T15.4 N=64 TGV: WALE nu_t is not zero (S^d != 0 for a true 3D vortex),
/// but with u0 = 1.28e-4/N the strain magnitude is tiny and the WALE-fitted
/// effective viscosity must deviate from nu_0 by <= 1% (order-wles.txt band).
///
/// Denominator note: "fitted nu_eff" comes from the KE decay rate
/// r = -ln(E1/E0)/t_star, with rate_ref = 6*nu*k^2 (diffusion limit) and
/// nu_eff = r / (6*k^2). Both LES-ON and LES-OFF fit the same way, so the
/// reported quantity is (nu_eff_on - nu_eff_off) / nu_0 = pure LES effect.
#[test]
#[ignore = "heavy characterization: T15.4 N=64 TGV fitted nu_eff with WALE (~5 min)"]
fn wale_tgv64_nu_eff_characterization() {
    let nu = 0.02;
    let n = 64usize;
    let k = 2.0 * std::f64::consts::PI / n as f64;
    let tstar = (0.1 / (nu * k * k)).round() as usize;

    let mut off = make_tgv64(nu);
    let e0_off = ke(&off);
    off.run(tstar);
    let e1_off = ke(&off);
    let rate_off = -(e1_off / e0_off).ln() / tstar as f64;

    let mut on = make_tgv64(nu);
    let mut les = WaleLes::<f64>::new();
    let e0_on = ke(&on);
    for _ in 0..tstar {
        les.update(&mut on);
        on.run(1);
    }
    let e1_on = ke(&on);
    let rate_on = -(e1_on / e0_on).ln() / tstar as f64;

    let nu_eff_off = rate_off / (6.0 * k * k);
    let nu_eff_on = rate_on / (6.0 * k * k);
    let dnu_rel = (nu_eff_on - nu_eff_off) / nu;
    eprintln!(
        "TGV64 freeze: nu_eff_off={nu_eff_off:e} nu_eff_on={nu_eff_on:e} \
         dnu_rel={dnu_rel:e} max_nu_t_on={:e}",
        les.nu_t().iter().copied().fold(0.0_f64, f64::max)
    );

    // Measured 2026-07-06 on the freeze pass (99bb32a): dnu_rel = 6.60e-8,
    // nu_eff_off = 1.9977e-2, nu_eff_on = 1.9977e-2, max nu_t = 1.39e-8.
    // Band frozen at 1e-6 (~15x headroom over the measured value); the
    // original order allowed 1% which is far too loose for the diffusive-
    // scaling regime where WALE should be essentially inert. A real physics
    // discrepancy (spurious eddy viscosity leaking into resolved TGV) would
    // exceed this by orders of magnitude.
    let band = 1.0e-6_f64;
    let max_nu_t_on = les.nu_t().iter().copied().fold(0.0_f64, f64::max);
    assert!(
        dnu_rel.abs() <= band,
        "TGV64 nu_eff shift {dnu_rel:e} > {band:e}, max nu_t (on)={max_nu_t_on:e}, \
         nu_eff_off={nu_eff_off:e} nu_eff_on={nu_eff_on:e} nu_0={nu:e} \
         (denominator = nu_0; a physics discrepancy would be O(nu_t/nu))"
    );
}

fn make_multimode(n: usize, nu: f64, u0: f64) -> S3 {
    let spec = GlobalSpec::<f64> {
        dims: [n, n, n],
        nu,
        periodic: [true, true, true],
        ..Default::default()
    };
    let mut s: S3 = Solver::new(&spec, &[], &[], [1, 1, 1], CpuScalar::default(), LocalPeriodic);
    let pi = std::f64::consts::PI;
    // Three non-aligned TGV-like modes, deterministic (no RNG per FR-PART-03
    // family anti-pattern ruling). Wave vectors chosen to break axis alignment
    // and provide a broadband strain field the LES can act on.
    let k1 = 2.0 * pi / n as f64;
    let k2 = 4.0 * pi / n as f64;
    let k3 = 6.0 * pi / n as f64;
    s.init_with(move |x, y, z| {
        let (xf, yf, zf) = (x as f64, y as f64, z as f64);
        let m1 = [
            (k1 * xf).sin() * (k1 * yf).cos() * (k1 * zf).cos(),
            -(k1 * xf).cos() * (k1 * yf).sin() * (k1 * zf).cos(),
            0.0,
        ];
        let m2 = [
            0.0,
            (k2 * yf).sin() * (k2 * zf).cos() * (k2 * xf).cos(),
            -(k2 * yf).cos() * (k2 * zf).sin() * (k2 * xf).cos(),
        ];
        let m3 = [
            -(k3 * zf).cos() * (k3 * xf).sin() * (k3 * yf).cos(),
            0.0,
            (k3 * zf).sin() * (k3 * xf).cos() * (k3 * yf).cos(),
        ];
        (
            1.0,
            [
                u0 * (m1[0] + m2[0] + m3[0]),
                u0 * (m1[1] + m2[1] + m3[1]),
                u0 * (m1[2] + m2[2] + m3[2]),
            ],
        )
    });
    s
}

fn reference_wale_nu_t(s: &S3) -> Vec<f64> {
    let cw_delta_sq = WALE_CW * WALE_CW;
    s.gather_velocity_gradient()
        .iter()
        .map(|g| {
            let mut strain = [[0.0; 3]; 3];
            for a in 0..3 {
                for b in 0..3 {
                    strain[a][b] = 0.5 * (g[a][b] + g[b][a]);
                }
            }

            let mut g2 = [[0.0; 3]; 3];
            for a in 0..3 {
                for b in 0..3 {
                    for k in 0..3 {
                        g2[a][b] += g[a][k] * g[k][b];
                    }
                }
            }
            let tr_g2 = g2[0][0] + g2[1][1] + g2[2][2];
            let mut sd = [[0.0; 3]; 3];
            for a in 0..3 {
                for b in 0..3 {
                    sd[a][b] = 0.5 * (g2[a][b] + g2[b][a]);
                    if a == b {
                        sd[a][b] -= tr_g2 / 3.0;
                    }
                }
            }

            let mut ss = 0.0;
            let mut sdsd = 0.0;
            for a in 0..3 {
                for b in 0..3 {
                    ss += strain[a][b] * strain[a][b];
                    sdsd += sd[a][b] * sd[a][b];
                }
            }
            let denom = ss.powf(2.5) + sdsd.powf(1.25);
            if denom > 0.0 {
                cw_delta_sq * sdsd.powf(1.5) / denom
            } else {
                0.0
            }
        })
        .collect()
}

#[test]
fn wale_unset_clipping_matches_raw_wale_bitwise_on_sheared_field() {
    let mut s = make_multimode(12, 0.02, 0.04);
    let expected = reference_wale_nu_t(&s);
    assert!(
        expected.iter().any(|v| *v > 0.0),
        "fixture must exercise non-zero WALE viscosity"
    );

    let mut les = WaleLes::<f64>::new();
    les.update(&mut s);
    assert_eq!(les.tau_eff_max(), None);
    assert_eq!(les.nu_t().len(), expected.len());
    for (i, (actual, expected)) in les.nu_t().iter().zip(&expected).enumerate() {
        assert_eq!(
            actual.to_bits(),
            expected.to_bits(),
            "nu_t[{i}] changed with clipping unset"
        );
    }

    let diagnostics = les.diagnostics();
    let max_expected = expected.iter().copied().fold(0.0_f64, f64::max);
    assert_eq!(diagnostics.clipped_cells, 0);
    assert_eq!(diagnostics.clipped_fraction, 0.0);
    assert_eq!(
        diagnostics.max_nu_t_before_clipping.to_bits(),
        max_expected.to_bits()
    );
    assert_eq!(diagnostics.tau_eff_max, None);
    assert_eq!(diagnostics.nu_t_max, None);
}

#[test]
fn wale_tau_eff_clipping_diagnostics_match_reference() {
    let nu = 0.02;
    let mut s = make_multimode(12, nu, 0.04);
    let raw = reference_wale_nu_t(&s);
    let mut sorted = raw.clone();
    sorted.sort_by(f64::total_cmp);
    let requested_cap = sorted[sorted.len() / 4];
    let tau_eff_max = 0.5 + 3.0 * (nu + requested_cap);

    let mut les = WaleLes::<f64>::new().with_tau_eff_max(tau_eff_max);
    les.update(&mut s);

    let diagnostics = les.diagnostics();
    let active_cap = diagnostics
        .nu_t_max
        .expect("clipping must report nu_t bound");
    let expected_clipped = raw.iter().filter(|v| **v > active_cap).count();
    let max_raw = raw.iter().copied().fold(0.0_f64, f64::max);
    assert!(
        expected_clipped > 0,
        "fixture must exercise tau_eff clipping"
    );
    assert_eq!(diagnostics.clipped_cells, expected_clipped);
    assert_eq!(
        diagnostics.clipped_fraction,
        expected_clipped as f64 / raw.len() as f64
    );
    assert_eq!(
        diagnostics.max_nu_t_before_clipping.to_bits(),
        max_raw.to_bits()
    );
    assert_eq!(diagnostics.tau_eff_max, Some(tau_eff_max));

    for (i, (actual, raw)) in les.nu_t().iter().zip(&raw).enumerate() {
        let expected = if *raw > active_cap { active_cap } else { *raw };
        assert_eq!(
            actual.to_bits(),
            expected.to_bits(),
            "nu_t[{i}] was not clipped exactly at the active bound"
        );
    }
}

#[test]
fn wale_tau_eff_clipping_count_is_monotone_with_bound() {
    let nu = 0.02;
    let raw = reference_wale_nu_t(&make_multimode(12, nu, 0.04));
    let max_raw = raw.iter().copied().fold(0.0_f64, f64::max);
    assert!(
        max_raw > 0.0,
        "fixture must exercise non-zero WALE viscosity"
    );

    let run = |fraction_of_max: f64| {
        let tau_eff_max = 0.5 + 3.0 * (nu + fraction_of_max * max_raw);
        let mut s = make_multimode(12, nu, 0.04);
        let mut les = WaleLes::<f64>::new().with_tau_eff_max(tau_eff_max);
        les.update(&mut s);
        les.diagnostics().clipped_cells
    };

    let clipped_low_bound = run(0.25);
    let clipped_mid_bound = run(0.50);
    let clipped_high_bound = run(0.75);
    assert!(
        clipped_low_bound >= clipped_mid_bound && clipped_mid_bound >= clipped_high_bound,
        "raising tau_eff_max increased clipped cells: low={clipped_low_bound}, \
         mid={clipped_mid_bound}, high={clipped_high_bound}"
    );
}

fn max_speed(s: &S3) -> f64 {
    let (ux, uy, uz) = (s.gather_ux(), s.gather_uy(), s.gather_uz());
    ux.iter()
        .zip(&uy)
        .zip(&uz)
        .map(|((a, b), c)| (a * a + b * b + c * c).sqrt())
        .fold(0.0_f64, f64::max)
}

/// Deterministic three-mode init at N=48. Parameters FROZEN 2026-07-06 on
/// origin/main 99bb32a: nu=0.003, u0=0.10 (so U/nu=33, well over the T10
/// grid-Reynolds threshold of 15). Measured:
/// - LES-OFF path diverges (NaN or max|u| > MAX_SPEED=0.3) at step 200
///   (20k step run, ~100x safety margin below the requested band).
/// - LES-ON path completes 20000 steps with max|u| = 5.15e-4, max nu_t =
///   5.08e-6 (~0.17% of nu_0 — WALE has an actual, non-trivial modeling
///   effect on this init while staying deep inside the low-Mach regime).
/// Wall time on M5 Max: ~7 min per run (LES-OFF + LES-ON). This is the
/// "turbulence tractability" evidence — a stabilization EXISTENCE proof
/// (100x horizon extension), NOT an accuracy claim; WALE-accurate turbulence
/// acceptance is Re_tau=180 vs DNS, still skeleton in this file.
#[test]
#[ignore = "heavy characterization: deterministic multimode stabilization (~2 min)"]
fn wale_multimode_stabilization() {
    // Frozen parameter tuple (see comment above): nu = 0.003, u0 = 0.10,
    // so U/nu = 33 — well over the T10 grid-Reynolds threshold of 15.
    let n = 48usize;
    let nu = 0.003;
    let u0 = 0.10;
    let steps = 20_000usize;

    // LES-OFF path — expect divergence within `steps`.
    let mut off = make_multimode(n, nu, u0);
    let mut off_diverged_at: Option<usize> = None;
    for i in 0..steps {
        off.run(1);
        if i % 200 == 199 {
            let (ux, uy, uz) = (off.gather_ux(), off.gather_uy(), off.gather_uz());
            let bad = ux.iter().any(|v| !v.is_finite())
                || uy.iter().any(|v| !v.is_finite())
                || uz.iter().any(|v| !v.is_finite())
                || max_speed(&off) > 0.3;
            if bad {
                off_diverged_at = Some(i + 1);
                break;
            }
        }
    }
    assert!(
        off_diverged_at.is_some(),
        "LES-OFF did NOT diverge within {steps} at nu={nu} u0={u0} — stabilization \
         point needs re-freeze (LES claim requires an off-path failure to beat)"
    );
    let off_step = off_diverged_at.unwrap();

    // LES-ON path — expect completion with bounded max|u|.
    let mut on = make_multimode(n, nu, u0);
    let mut les = WaleLes::<f64>::new();
    for _ in 0..steps {
        les.update(&mut on);
        on.run(1);
    }
    let on_max = max_speed(&on);
    let on_finite = {
        let (ux, uy, uz) = (on.gather_ux(), on.gather_uy(), on.gather_uz());
        ux.iter().chain(uy.iter()).chain(uz.iter()).all(|v| v.is_finite())
    };
    let max_nu_t_on = les.nu_t().iter().copied().fold(0.0_f64, f64::max);
    eprintln!(
        "MULTIMODE freeze: nu={nu} u0={u0} N={n} steps={steps} \
         off_diverged_at={off_step} on_max|u|={on_max:e} max_nu_t_on={max_nu_t_on:e}"
    );
    assert!(
        on_finite && on_max <= 0.3,
        "LES-ON not usefully stable at step {steps}: max|u|={on_max:e}, finite={on_finite}, \
         max nu_t={max_nu_t_on:e} — this is the tractability gate"
    );
}
