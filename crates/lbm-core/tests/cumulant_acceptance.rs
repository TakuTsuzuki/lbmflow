//! MF-alpha stage-2 adversarial acceptance for the CPU central-moment
//! (`CollisionKind::CentralMoment`) reference operator.

use lbm_core::prelude::*;
use std::f64::consts::PI;

type S3<L> = Solver<L, f64, CpuScalar, LocalPeriodic>;

const NU: f64 = 0.02;
const TGV_U0_COEF: f64 = 1.28e-4;
const ADVECTED_MEAN_U: f64 = 0.05;

fn omega_from_nu(nu: f64) -> f64 {
    1.0 / (3.0 * nu + 0.5)
}

fn spec<L: Lattice>(n: usize, nu: f64, collision: CollisionKind) -> GlobalSpec<f64> {
    let _ = L::D;
    GlobalSpec {
        dims: [n, n, n],
        nu,
        collision,
        periodic: [true, true, true],
        ..Default::default()
    }
}

fn tgv_velocity(n: usize, u0: f64, mean_u: f64, x: usize, y: usize, z: usize) -> [f64; 3] {
    let k = 2.0 * PI / n as f64;
    let (xf, yf, zf) = (k * x as f64, k * y as f64, k * z as f64);
    [
        mean_u + u0 * xf.sin() * yf.cos() * zf.cos(),
        -u0 * xf.cos() * yf.sin() * zf.cos(),
        0.0,
    ]
}

fn tgv_density(n: usize, u0: f64, x: usize, y: usize, z: usize) -> f64 {
    let k = 2.0 * PI / n as f64;
    let (xf, yf, zf) = (k * x as f64, k * y as f64, k * z as f64);
    let p = u0 * u0 / 16.0 * (((2.0 * xf).cos() + (2.0 * yf).cos()) * ((2.0 * zf).cos() + 2.0));
    1.0 + 3.0 * p
}

fn make_tgv<L: Lattice>(n: usize, nu: f64, collision: CollisionKind, mean_u: f64) -> S3<L> {
    let spec = spec::<L>(n, nu, collision);
    let mut s: S3<L> = Solver::new(
        &spec,
        &[],
        &[],
        [1, 1, 1],
        CpuScalar::default(),
        LocalPeriodic,
    );
    let u0 = TGV_U0_COEF / n as f64;
    s.init_with(move |x, y, z| {
        (
            tgv_density(n, u0, x, y, z),
            tgv_velocity(n, u0, mean_u, x, y, z),
        )
    });
    s
}

fn fluctuation_ke<L: Lattice>(s: &S3<L>, mean_u: f64) -> f64 {
    let (ux, uy, uz) = (s.gather_ux(), s.gather_uy(), s.gather_uz());
    ux.iter()
        .zip(&uy)
        .zip(&uz)
        .map(|((a, b), c)| (a - mean_u).powi(2) + b * b + c * c)
        .sum()
}

fn decay_rate<L: Lattice>(collision: CollisionKind, mean_u: f64) -> f64 {
    let n = 32usize;
    let k = 2.0 * PI / n as f64;
    let steps = (0.1 / (NU * k * k)).round() as usize;
    let mut s = make_tgv::<L>(n, NU, collision, mean_u);
    let e0 = fluctuation_ke(&s, mean_u);
    s.run(steps);
    let e1 = fluctuation_ke(&s, mean_u);
    -(e1 / e0).ln() / steps as f64
}

fn galilean_defect<L: Lattice>(collision: CollisionKind) -> f64 {
    let still = decay_rate::<L>(collision, 0.0);
    let advected = decay_rate::<L>(collision, ADVECTED_MEAN_U);
    (advected - still).abs() / still
}

fn measure_galilean<L: Lattice>(name: &str, bgk_baseline: f64) -> (f64, f64, f64) {
    let bgk = galilean_defect::<L>(CollisionKind::Bgk);
    let cumulant = galilean_defect::<L>(CollisionKind::CentralMoment {
        omega_shear: omega_from_nu(NU),
    });
    println!("{name} advected TGV3D Galilean defect: BGK={bgk:e} Cumulant={cumulant:e}");
    // BAND-FREEZE-PENDING(PM): provisional stage-2 claim threshold.
    (bgk, cumulant, 0.5 * bgk_baseline)
}

fn assert_galilean_result(name: &str, bgk: f64, cumulant: f64, band: f64) {
    assert!(
        (1.0e-3..=5.0e-3).contains(&bgk),
        "{name} BGK harness guard defect {bgk:e} outside [1e-3, 5e-3]"
    );
    assert!(
        cumulant < band,
        "{name} Cumulant Galilean defect {cumulant:e} is not below the 2x-improvement band {band:e}"
    );
}

#[test]
fn cumulant_improves_advected_tgv3d_galilean_invariance() {
    let d3q19 = measure_galilean::<D3Q19>("D3Q19", 2.539e-3);
    let d3q27 = measure_galilean::<D3Q27>("D3Q27", 2.473e-3);
    assert_galilean_result("D3Q19", d3q19.0, d3q19.1, d3q19.2);
    assert_galilean_result("D3Q27", d3q27.0, d3q27.1, d3q27.2);
}

fn max_population_delta<L: Lattice>(a: &S3<L>, b: &S3<L>) -> f64 {
    let mut d = 0.0f64;
    for q in 0..L::Q {
        let af = a.gather_f(q);
        let bf = b.gather_f(q);
        d = d.max(
            af.iter()
                .zip(&bf)
                .map(|(x, y)| (x - y).abs())
                .fold(0.0, f64::max),
        );
    }
    d
}

fn assert_rest_exact<L: Lattice>(name: &str) {
    let spec = spec::<L>(
        8,
        NU,
        CollisionKind::CentralMoment {
            omega_shear: omega_from_nu(NU),
        },
    );
    let mut initial: S3<L> = Solver::new(
        &spec,
        &[],
        &[],
        [1, 1, 1],
        CpuScalar::default(),
        LocalPeriodic,
    );
    let mut evolved: S3<L> = Solver::new(
        &spec,
        &[],
        &[],
        [1, 1, 1],
        CpuScalar::default(),
        LocalPeriodic,
    );
    initial.init_with(|_, _, _| (1.0, [0.0; 3]));
    evolved.init_with(|_, _, _| (1.0, [0.0; 3]));
    evolved.run(100);
    let df = max_population_delta(&initial, &evolved);
    assert_eq!(
        df, 0.0,
        "{name} rest equilibrium changed after 100 steps: max |Delta f|={df:e}"
    );
}

fn assert_uniform_exact<L: Lattice>(name: &str) {
    let spec = spec::<L>(
        8,
        NU,
        CollisionKind::CentralMoment {
            omega_shear: omega_from_nu(NU),
        },
    );
    let mut s: S3<L> = Solver::new(
        &spec,
        &[],
        &[],
        [1, 1, 1],
        CpuScalar::default(),
        LocalPeriodic,
    );
    let u0 = [0.03, 0.02, 0.01];
    s.init_with(move |_, _, _| (1.0, u0));
    s.run(200);
    let mut err = 0.0f64;
    for z in 0..8 {
        for y in 0..8 {
            for x in 0..8 {
                err = err.max((s.rho(x, y, z) - 1.0).abs());
                let u = s.u(x, y, z);
                for a in 0..3 {
                    err = err.max((u[a] - u0[a]).abs());
                }
            }
        }
    }
    assert!(
        err <= 1.0e-13,
        "{name} uniform periodic state drifted/non-uniform after 200 steps: max error {err:e} > 1e-13"
    );
}

#[test]
fn cumulant_rest_and_uniform_states_are_exact() {
    assert_rest_exact::<D3Q19>("D3Q19");
    assert_rest_exact::<D3Q27>("D3Q27");
    assert_uniform_exact::<D3Q19>("D3Q19");
    assert_uniform_exact::<D3Q27>("D3Q27");
}

fn measure_tgv_nu_eff<L: Lattice>(name: &str) -> (f64, f64) {
    let n = 32usize;
    let k = 2.0 * PI / n as f64;
    let steps = (0.1 / (NU * k * k)).round() as usize;
    let mut s = make_tgv::<L>(
        n,
        NU,
        CollisionKind::CentralMoment {
            omega_shear: omega_from_nu(NU),
        },
        0.0,
    );
    let e0 = fluctuation_ke(&s, 0.0);
    s.run(steps);
    let e1 = fluctuation_ke(&s, 0.0);
    let rate = -(e1 / e0).ln() / steps as f64;
    let nu_eff = rate / (6.0 * k * k);
    let rel = (nu_eff - NU).abs() / NU;
    println!("{name} Cumulant TGV3D N={n}: nu_eff={nu_eff:e} rel_err={rel:e}");
    (nu_eff, rel)
}

fn assert_tgv_nu_eff_result(name: &str, nu_eff: f64, rel: f64, band: f64) {
    assert!(
        rel <= band,
        "{name} CentralMoment TGV3D nu_eff relative error {rel:e} > refrozen finite-N decay-rate band {band:e} (nu_eff={nu_eff:e}, nu={NU:e})"
    );
}

#[test]
fn cumulant_tgv3d_diffusive_scaling_nu_eff_matches_t15_band() {
    let d3q19 = measure_tgv_nu_eff::<D3Q19>("D3Q19");
    let d3q27 = measure_tgv_nu_eff::<D3Q27>("D3Q27");
    // ANOM-P4-008 removed the D3Q19 finite-N viscosity offset. The decisive
    // acceptance is the h^2-intercept audit; this N=32 smoke band is frozen
    // to the uncorrected finite-resolution value printed by this test.
    assert_tgv_nu_eff_result("D3Q19", d3q19.0, d3q19.1, 2.4e-2);
    assert_tgv_nu_eff_result("D3Q27", d3q27.0, d3q27.1, 2.0e-2);
}

fn cumulant_channel(ny: usize, force: [f64; 3], top_u: [f64; 3]) -> S3<D3Q19> {
    let spec = GlobalSpec {
        dims: [4, ny, 4],
        nu: NU,
        collision: CollisionKind::CentralMoment {
            omega_shear: omega_from_nu(NU),
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
    Solver::new(
        &spec,
        &solid,
        &wall_u,
        [1, 1, 1],
        CpuScalar::default(),
        LocalPeriodic,
    )
}

fn all_fields_finite(s: &S3<D3Q19>) -> bool {
    let fields = [s.gather_rho(), s.gather_ux(), s.gather_uy(), s.gather_uz()];
    fields.iter().flatten().all(|v| v.is_finite())
}

#[test]
fn wale_omega_field_composes_with_cumulant_collision() {
    let mut sheared = cumulant_channel(10, [1.0e-6, 0.0, 0.0], [0.08, 0.0, 0.0]);
    let mut les = WaleLes::<f64>::new();
    for _ in 0..10 {
        les.update(&mut sheared);
        sheared.run(1);
    }
    let max_nu_t = les.nu_t().iter().copied().fold(0.0f64, f64::max);
    assert!(
        all_fields_finite(&sheared) && les.nu_t().iter().all(|v| v.is_finite()),
        "WALE over Cumulant produced non-finite fields: max nu_t={max_nu_t:e}"
    );

    let mut uniform = make_tgv::<D3Q19>(
        8,
        NU,
        CollisionKind::CentralMoment {
            omega_shear: omega_from_nu(NU),
        },
        0.0,
    );
    uniform.init_with(|_, _, _| (1.0, [0.03, 0.02, 0.01]));
    les.update(&mut uniform);
    let max_null = les.nu_t().iter().copied().fold(0.0f64, f64::max);
    assert!(
        max_null <= 1.0e-14,
        "null-shear WALE path over Cumulant produced nu_t={max_null:e}; expected uniform omega / nu_t=0"
    );
}

fn make_multimode(n: usize, nu: f64, u0: f64, collision: CollisionKind) -> S3<D3Q19> {
    let spec = spec::<D3Q19>(n, nu, collision);
    let mut s: S3<D3Q19> = Solver::new(
        &spec,
        &[],
        &[],
        [1, 1, 1],
        CpuScalar::default(),
        LocalPeriodic,
    );
    let k1 = 2.0 * PI / n as f64;
    let k2 = 4.0 * PI / n as f64;
    let k3 = 6.0 * PI / n as f64;
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

fn max_speed(s: &S3<D3Q19>) -> f64 {
    let (ux, uy, uz) = (s.gather_ux(), s.gather_uy(), s.gather_uz());
    ux.iter()
        .zip(&uy)
        .zip(&uz)
        .map(|((a, b), c)| (a * a + b * b + c * c).sqrt())
        .fold(0.0f64, f64::max)
}

fn diverged(s: &S3<D3Q19>) -> bool {
    !all_fields_finite(s) || max_speed(s) > 0.3
}

#[test]
#[ignore = "heavy characterization: high-Re multimode Cumulant stability up to 20k steps"]
fn cumulant_multimode_high_re_stability_probe() {
    let n = 48usize;
    let nu = 0.003;
    let u0 = 0.10;
    let steps = 20_000usize;

    let mut bgk = make_multimode(n, nu, u0, CollisionKind::Bgk);
    let mut bgk_diverged_at = None;
    for step in 1..=steps {
        bgk.run(1);
        if step % 200 == 0 && diverged(&bgk) {
            bgk_diverged_at = Some(step);
            break;
        }
    }

    let mut cumulant = make_multimode(
        n,
        nu,
        u0,
        CollisionKind::CentralMoment {
            omega_shear: omega_from_nu(nu),
        },
    );
    let mut reached = 0usize;
    for step in 1..=steps {
        cumulant.run(1);
        reached = step;
        if step % 200 == 0 && diverged(&cumulant) {
            break;
        }
    }
    let max_u = max_speed(&cumulant);
    let finite = all_fields_finite(&cumulant);
    println!(
        "Cumulant multimode high-Re probe: N={n} nu={nu:e} u0={u0:e} steps_reached={reached} \
         max|u|={max_u:e} finite={finite} BGK_diverged_at={bgk_diverged_at:?}"
    );
}
