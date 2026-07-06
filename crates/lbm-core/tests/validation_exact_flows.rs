//! Lane 4.2 exact/series validation gates:
//! Kovasznay steady Navier-Stokes, Womersley pulsatile channel, and the
//! Sangani-Acrivos square cylinder-array permeability series.
//!
//! These are intentionally adversarial tests. The comments derive the analytic
//! target used by each assertion; acceptance bands are benchmark requirements,
//! not fitted constants.

mod common;

use common::metrics::*;
use common::run_to_steady;
use lbm_core::compat::prelude::*;
use std::f64::consts::PI;

const TRT: Collision = Collision::Trt {
    magic: Collision::MAGIC_STD,
};
const CS2: f64 = 1.0 / 3.0;

#[derive(Clone, Copy, Debug)]
struct Complex {
    re: f64,
    im: f64,
}

impl Complex {
    fn new(re: f64, im: f64) -> Self {
        Self { re, im }
    }

    fn sub(self, rhs: Self) -> Self {
        Self::new(self.re - rhs.re, self.im - rhs.im)
    }

    fn mul(self, rhs: Self) -> Self {
        Self::new(
            self.re * rhs.re - self.im * rhs.im,
            self.re * rhs.im + self.im * rhs.re,
        )
    }

    fn div(self, rhs: Self) -> Self {
        let den = rhs.re * rhs.re + rhs.im * rhs.im;
        Self::new(
            (self.re * rhs.re + self.im * rhs.im) / den,
            (self.im * rhs.re - self.re * rhs.im) / den,
        )
    }

    fn sqrt(self) -> Self {
        let r = (self.re * self.re + self.im * self.im).sqrt();
        let re = ((r + self.re) / 2.0).sqrt();
        let im = self.im.signum() * ((r - self.re) / 2.0).sqrt();
        Self::new(re, im)
    }

    fn cosh(self) -> Self {
        Self::new(
            self.re.cosh() * self.im.cos(),
            self.re.sinh() * self.im.sin(),
        )
    }
}

fn phase_delta(a: f64, b: f64) -> f64 {
    (a - b + PI).rem_euclid(2.0 * PI) - PI
}

// ---------------------------------------------------------------------------
// G1 Kovasznay flow
// ---------------------------------------------------------------------------

fn kovasznay_lambda(re: f64) -> f64 {
    re / 2.0 - (re * re / 4.0 + 4.0 * PI * PI).sqrt()
}

fn kovasznay_velocity(n: usize, u0: f64, re: f64, x: usize, y: usize) -> [f64; 2] {
    // The exact solution is written in nondimensional coordinates x', y':
    // u = U0 [1 - exp(lambda x') cos(2 pi y')]
    // v = U0 [lambda/(2 pi)] exp(lambda x') sin(2 pi y') .
    //
    // We use the mathematically valid one-period window x' in [0, 1],
    // y' in [-1/2, 1/2).  Left/right open boundaries live on the actual edge
    // columns, so x' = x/(N-1).  The y direction is periodic, so its lattice
    // period is N cells and y' = y/N - 1/2 maps y=0..N-1 to one full period.
    // No half-way-wall offset is involved because this setup has no y walls.
    let xp = x as f64 / (n - 1) as f64;
    let yp = y as f64 / n as f64 - 0.5;
    let lambda = kovasznay_lambda(re);
    let e = (lambda * xp).exp();
    [
        u0 * (1.0 - e * (2.0 * PI * yp).cos()),
        u0 * lambda / (2.0 * PI) * e * (2.0 * PI * yp).sin(),
    ]
}

fn kovasznay_pressure(n: usize, u0: f64, re: f64, x: usize) -> f64 {
    // p = p0 + U0^2/2 [1 - exp(2 lambda x')].  The right pressure outlet
    // prescribes rho=1, so choose p0 such that p(x'=1)=0 and rho=1+p/cs^2.
    // The validation fit subtracts means, so this arbitrary pressure gauge
    // cannot affect the measured pressure-field slope or R^2.
    let xp = x as f64 / (n - 1) as f64;
    let lambda = kovasznay_lambda(re);
    let p_raw = 0.5 * u0 * u0 * (1.0 - (2.0 * lambda * xp).exp());
    let p_right = 0.5 * u0 * u0 * (1.0 - (2.0 * lambda).exp());
    p_raw - p_right
}

fn build_kovasznay(n: usize, u0: f64, nu: f64) -> Simulation<f64> {
    let re = u0 * n as f64 / nu;
    let mut sim: Simulation<f64> = SimConfig {
        nx: n,
        ny: n,
        nu,
        collision: TRT,
        edges: Edges {
            left: EdgeBC::VelocityInlet { u: [0.0, 0.0] },
            right: EdgeBC::PressureOutlet { rho: 1.0 },
            bottom: EdgeBC::Periodic,
            top: EdgeBC::Periodic,
        },
        ..Default::default()
    }
    .build()
    .unwrap();
    sim.set_inlet_profile(Edge::Left, |y| kovasznay_velocity(n, u0, re, 0, y));
    sim.init_with(|x, y| {
        let [ux, uy] = kovasznay_velocity(n, u0, re, x, y);
        (1.0 + kovasznay_pressure(n, u0, re, x) / CS2, ux, uy)
    });
    sim
}

fn kovasznay_measure(n: usize, u0: f64, nu: f64) -> (f64, f64, f64, f64, bool, u64) {
    let re = u0 * n as f64 / nu;
    let mut sim = build_kovasznay(n, u0, nu);
    let steady = run_to_steady(&mut sim, 500, 1.0e-11, 40_000);

    let mut u_actual = Vec::new();
    let mut u_ref = Vec::new();
    let mut v_actual = Vec::new();
    let mut v_ref = Vec::new();
    let mut p_actual = Vec::new();
    let mut p_ref = Vec::new();
    for y in 0..n {
        for x in 0..(n - 4) {
            let [ux, uy] = kovasznay_velocity(n, u0, re, x, y);
            u_actual.push(sim.ux(x, y));
            u_ref.push(ux);
            v_actual.push(sim.uy(x, y));
            v_ref.push(uy);
            p_actual.push(CS2 * sim.rho(x, y));
            p_ref.push(kovasznay_pressure(n, u0, re, x));
        }
    }
    let mean_pa = p_actual.iter().sum::<f64>() / p_actual.len() as f64;
    let mean_pr = p_ref.iter().sum::<f64>() / p_ref.len() as f64;
    for p in &mut p_actual {
        *p -= mean_pa;
    }
    for p in &mut p_ref {
        *p -= mean_pr;
    }
    let pressure_fit = linear_fit(&p_ref, &p_actual);
    let eu = l2_rel(&u_actual, &u_ref);
    let ev = l2_rel(&v_actual, &v_ref);
    println!(
        "VAL EXACT KOVASZNAY: N={n} U0={u0:.9e} nu={nu:.9e} Re={re:.6} steps={} steady={steady} criterion=1e-11 L2rel_u={eu:.9e} L2rel_v={ev:.9e} pressure_slope={:.9e} pressure_r2={:.9e} bulk_excludes_outlet_cols=4",
        sim.time(),
        pressure_fit.slope,
        pressure_fit.r2
    );
    (
        eu,
        ev,
        pressure_fit.slope,
        pressure_fit.r2,
        steady,
        sim.time(),
    )
}

#[test]
fn g1_kovasznay_light_exact_velocity_and_pressure_field() {
    let (eu, ev, slope, r2, steady, steps) = kovasznay_measure(64, 0.05, 0.05 * 64.0 / 20.0);
    assert!(
        steady,
        "VAL EXACT KOVASZNAY light steady=false after {steps} steps, band=true, criterion=max|du|/max|u| <= 1e-11"
    );
    assert!(
        eu <= 2.0e-3,
        "VAL EXACT KOVASZNAY light L2rel(u)={eu:.9e}, band=2.0e-3, normalization=||u_exact||2 bulk excluding 4 outlet columns"
    );
    assert!(
        ev <= 2.0e-2,
        "VAL EXACT KOVASZNAY light L2rel(v)={ev:.9e}, band=2.0e-2, normalization=||v_exact||2 bulk excluding 4 outlet columns"
    );
    assert!(
        (0.9..=1.1).contains(&slope),
        "VAL EXACT KOVASZNAY pressure slope={slope:.9e}, band=[0.9,1.1], measured=cs^2*(rho-rho_mean), reference=p_exact-p_mean"
    );
    assert!(
        r2 >= 0.99,
        "VAL EXACT KOVASZNAY pressure r2={r2:.9e}, band>=0.99, fit=measured pressure vs exact pressure on bulk excluding 4 outlet columns"
    );
}

#[test]
#[ignore = "heavy VAL-EXACT Kovasznay N={48,96,192} convergence ladder"]
fn g1_kovasznay_heavy_second_order_ladder() {
    // Fixed Re requires U0*N/nu = const.  This ladder holds nu fixed at the
    // N=48 light-Mach value nu=0.05*48/Re=0.12 and scales U0=0.05*48/N, so
    // the nondimensional Kovasznay field is unchanged while Ma decreases with
    // refinement.  The spacing h is proportional to 1/N.
    let nu = 0.05 * 48.0 / 20.0;
    let ns = [48usize, 96, 192];
    let mut hs = Vec::new();
    let mut errs = Vec::new();
    for n in ns {
        let u0 = 0.05 * 48.0 / n as f64;
        let (eu, ev, slope, r2, steady, steps) = kovasznay_measure(n, u0, nu);
        println!(
            "VAL EXACT KOVASZNAY HEAVY row: N={n} h={:.9e} L2rel_u={eu:.9e} L2rel_v={ev:.9e} pressure_slope={slope:.9e} pressure_r2={r2:.9e}",
            1.0 / n as f64
        );
        assert!(
            steady,
            "VAL EXACT KOVASZNAY heavy N={n} steady=false after {steps} steps, band=true, criterion=max|du|/max|u| <= 1e-11"
        );
        hs.push(1.0 / n as f64);
        errs.push(eu);
    }
    let fit = order_fit(&hs, &errs);
    println!(
        "VAL EXACT KOVASZNAY HEAVY order: slope={:.9e} r2={:.9e} errors={errs:?}",
        fit.slope, fit.r2
    );
    assert!(
        (1.7..=2.3).contains(&fit.slope),
        "VAL EXACT KOVASZNAY heavy order slope={:.9e}, band=[1.7,2.3], normalization=L2rel(u)",
        fit.slope
    );
    assert!(
        fit.r2 >= 0.98,
        "VAL EXACT KOVASZNAY heavy order r2={:.9e}, band>=0.98, fit=log(error) vs log(h)",
        fit.r2
    );
}

// ---------------------------------------------------------------------------
// G2 Womersley pulsatile channel
// ---------------------------------------------------------------------------

fn womersley_velocity_coeff(y: f64, h: f64, nu: f64, omega: f64, f0: f64) -> Complex {
    // The x-momentum equation for a fully developed channel with an oscillating
    // uniform body force is du/dt = nu d2u/dy2 + F0 cos(omega t), u(+-h)=0.
    // With u = Re{A(y) exp(i omega t)} and F = Re{F0 exp(i omega t)}:
    // i omega A = nu A'' + F0.
    // A particular solution is F0/(i omega), and the wall-corrected solution is
    // A(y) = F0/(i omega) [1 - cosh(k y) / cosh(k h)],
    // k = sqrt(i omega / nu).  Alpha = h sqrt(omega/nu); alpha=3 is selected
    // here, then omega is rounded through the integer period T=round(2pi/omega)
    // so the phase projection samples an exact integer number of periods.
    let k = Complex::new(0.0, omega / nu).sqrt();
    let ratio = k
        .mul(Complex::new(y, 0.0))
        .cosh()
        .div(k.mul(Complex::new(h, 0.0)).cosh());
    Complex::new(1.0, 0.0)
        .sub(ratio)
        .mul(Complex::new(0.0, -f0 / omega))
}

fn womersley_amp_phase(y: f64, h: f64, nu: f64, omega: f64, f0: f64) -> (f64, f64) {
    let a = womersley_velocity_coeff(y, h, nu, omega, f0);
    // phase_fit uses signal ~= amplitude*sin(omega t + phase).  Since
    // Re{A exp(iwt)} = A.re*cos(wt) - A.im*sin(wt), the sin/cos coefficients
    // are (-A.im, A.re), hence phase=atan2(cos_coeff, sin_coeff).
    let amp = (a.re * a.re + a.im * a.im).sqrt();
    let phase = a.re.atan2(-a.im);
    (amp, phase)
}

fn womersley_measure(alpha: f64) -> (CurveAgreement, f64, f64, usize) {
    let ny = 34usize;
    let nx = 8usize;
    let h = (ny - 2) as f64 / 2.0;
    let nu = 0.02;
    let f0 = 1.0e-6;
    let omega_target = alpha * alpha * nu / (h * h);
    let period = (2.0 * PI / omega_target).round() as usize;
    let omega = 2.0 * PI / period as f64;
    let mut sim: Simulation<f64> = SimConfig {
        nx,
        ny,
        nu,
        collision: TRT,
        edges: Edges {
            left: EdgeBC::Periodic,
            right: EdgeBC::Periodic,
            bottom: EdgeBC::BounceBack,
            top: EdgeBC::BounceBack,
        },
        ..Default::default()
    }
    .build()
    .unwrap();
    sim.init_with(|_, _| (1.0, 0.0, 0.0));

    // ANOM-P2-001 interaction: the per-cell force-field path has a one-step
    // impulse deficit when the force changes.  This sinusoidal field-path test
    // intentionally quantifies the resulting amplitude/phase bias; if the
    // measured deviation has an F/2-per-step shape, that is a core finding, not
    // a reason to widen the band here.
    let warmup = 6 * period;
    let sample_steps = 4 * period;
    let mut times = Vec::with_capacity(sample_steps);
    let mut signals = vec![Vec::with_capacity(sample_steps); ny - 2];
    for step in 0..(warmup + sample_steps) {
        let force = f0 * (omega * step as f64).cos();
        sim.force_field_mut().fill([force, 0.0]);
        sim.step();
        if step >= warmup {
            times.push(sim.time() as f64);
            for y in 1..ny - 1 {
                signals[y - 1].push(sim.ux(nx / 2, y));
            }
        }
    }

    let mut amp_samples = Vec::new();
    let mut phase_center_err = 0.0f64;
    for y in 1..ny - 1 {
        let y_centered = y as f64 - 0.5 - h;
        let (amp, phase) = phase_fit(&times, &signals[y - 1], omega);
        let (amp_ref, phase_ref) = womersley_amp_phase(y_centered, h, nu, omega, f0);
        amp_samples.push((y_centered / h, amp));
        if y == ny / 2 || y + 1 == ny / 2 {
            phase_center_err = phase_center_err.max(phase_delta(phase, phase_ref).abs());
        }
        println!(
            "VAL EXACT WOMERSLEY row: alpha={alpha:.6} y={y} y_over_h={:.9e} amp={amp:.9e} amp_ref={amp_ref:.9e} phase={phase:.9e} phase_ref={phase_ref:.9e} phase_err={:.9e}",
            y_centered / h,
            phase_delta(phase, phase_ref).abs()
        );
    }
    let agreement = curve_agreement(
        |eta| womersley_amp_phase(eta * h, h, nu, omega, f0).0,
        &amp_samples,
        0.05,
        1.0e-14,
    );
    println!(
        "VAL EXACT WOMERSLEY: alpha={alpha:.6} h={h:.6} nu={nu:.9e} F0={f0:.9e} omega_target={omega_target:.9e} omega_rounded={omega:.9e} period={period} warmup_periods=6 sample_periods=4 amp_max_rel_dev={:.9e} amp_worst_y_over_h={:.9e} amp_frac_in_5pct={:.9e} center_phase_err={phase_center_err:.9e}",
        agreement.max_rel_dev, agreement.worst_x, agreement.frac_in_band
    );
    (agreement, phase_center_err, omega, period)
}

#[test]
fn g2_womersley_light_unsteady_amplitude_and_phase() {
    let (agreement, phase_center_err, _, _) = womersley_measure(3.0);
    assert!(
        agreement.max_rel_dev <= 0.05,
        "VAL EXACT WOMERSLEY light amplitude max_rel_dev={:.9e}, band=0.05, normalization=analytic amplitude with floor=1e-14, worst_y_over_h={:.9e}",
        agreement.max_rel_dev,
        agreement.worst_x
    );
    assert!(
        phase_center_err <= 0.05,
        "VAL EXACT WOMERSLEY light centerline phase-lag error={phase_center_err:.9e} rad, band=0.05 rad, phase convention=phase_fit sin(omega*t+phase); ANOM-P2-001 may contribute a force-change impulse-deficit bias"
    );
}

#[test]
#[ignore = "heavy VAL-EXACT Womersley alpha={2,4,8} amplitude/phase sweep"]
fn g2_womersley_heavy_alpha_sweep() {
    for alpha in [2.0, 4.0, 8.0] {
        let (agreement, phase_center_err, omega, period) = womersley_measure(alpha);
        println!(
            "VAL EXACT WOMERSLEY HEAVY row: alpha={alpha:.6} omega={omega:.9e} period={period} amp_max_rel_dev={:.9e} phase_center_err={phase_center_err:.9e}",
            agreement.max_rel_dev
        );
        assert!(
            agreement.max_rel_dev <= 0.05,
            "VAL EXACT WOMERSLEY heavy alpha={alpha} amplitude max_rel_dev={:.9e}, band=0.05, normalization=analytic amplitude profile",
            agreement.max_rel_dev
        );
        assert!(
            phase_center_err <= 0.05,
            "VAL EXACT WOMERSLEY heavy alpha={alpha} centerline phase-lag error={phase_center_err:.9e} rad, band=0.05 rad"
        );
    }
}

// ---------------------------------------------------------------------------
// G3 Sangani-Acrivos periodic cylinder array
// ---------------------------------------------------------------------------

fn sangani_series_s(phi: f64) -> f64 {
    -phi.ln() - 1.476_335_97 + 2.0 * phi - 1.774 * phi * phi + 4.076 * phi * phi * phi
}

#[derive(Clone, Copy, Debug)]
struct PermeabilityCase {
    phi_target: f64,
    n: usize,
    a_nominal: f64,
    bouzidi: bool,
}

fn case_for_phi(phi: f64, bouzidi: bool) -> PermeabilityCase {
    let a = 10.0;
    let n = (PI * a * a / phi).sqrt().round() as usize;
    PermeabilityCase {
        phi_target: phi,
        n,
        a_nominal: a,
        bouzidi,
    }
}

fn permeability_measure(case: PermeabilityCase) -> (f64, f64, f64, f64, f64, u64) {
    let nu = 0.1;
    let g = 1.0e-7;
    let n = case.n;
    let cx = (n as f64 - 1.0) / 2.0;
    let cy = (n as f64 - 1.0) / 2.0;
    let mut sim: Simulation<f64> = SimConfig {
        nx: n,
        ny: n,
        nu,
        collision: TRT,
        force: [g, 0.0],
        ..Default::default()
    }
    .build()
    .unwrap();
    sim.set_solid_region(|x, y| {
        let dx = x as f64 - cx;
        let dy = y as f64 - cy;
        dx * dx + dy * dy <= case.a_nominal * case.a_nominal
    });
    if case.bouzidi {
        sim.set_bouzidi_circle(cx, cy, case.a_nominal);
    }
    sim.init_with(|_, _| (1.0, 0.0, 0.0));
    let steady = run_to_steady(&mut sim, 1_000, 1.0e-11, 80_000);
    assert!(
        steady,
        "VAL EXACT SANGANI steady=false after {} steps, phi_target={:.6}, N={n}, bouzidi={}",
        sim.time(),
        case.phi_target,
        case.bouzidi
    );

    let solid = n * n - sim.fluid_cell_count();
    let phi_actual = solid as f64 / (n * n) as f64;
    let a_actual = (solid as f64 / PI).sqrt();
    let mean_u = sim.ux_field().iter().sum::<f64>() / (n * n) as f64;
    // Darcy's law in lattice units gives superficial U = k g / nu, hence
    // k = U nu/g.  The mean is over the total periodic cell area, not only
    // fluid cells, matching the superficial-velocity permeability convention.
    //
    // Sangani-Acrivos define S(phi)=4*pi/F* with
    // S=-ln(phi)-1.47633597+2phi-1.774phi^2+4.076phi^3+O(phi^4).
    // For the square 2-D array in the dilute Stokes limit this is equivalent
    // to k/a^2 = S/(8 phi).  To avoid freezing a disputed drag-normalization
    // constant, the assertion uses the ratio form S_measured = 8 phi k/a^2 and
    // compares S_measured directly against the published S(phi) polynomial.
    // phi_actual and a_actual come from the rasterized solid-cell count because
    // staircase geometry changes area by O(1/a).
    let k = mean_u * nu / g;
    let s_measured = 8.0 * phi_actual * k / (a_actual * a_actual);
    let s_ref = sangani_series_s(phi_actual);
    let re_a = mean_u.abs() * a_actual / nu;
    println!(
        "VAL EXACT SANGANI: phi_target={:.6} phi_actual={phi_actual:.9e} N={n} solid_cells={solid} a_actual={a_actual:.9e} bouzidi={} mean_u={mean_u:.9e} nu={nu:.9e} g={g:.9e} k={k:.9e} k_over_a2={:.9e} S_measured={s_measured:.9e} S_ref={s_ref:.9e} Re_a={re_a:.9e} steps={}",
        case.phi_target,
        case.bouzidi,
        k / (a_actual * a_actual),
        sim.time()
    );
    (phi_actual, k, s_measured, s_ref, re_a, sim.time())
}

#[test]
fn g3_sangani_acrivos_light_phi_010_canary() {
    let case = case_for_phi(0.10, false);
    let (phi, _k, s_measured, _s_ref, re_a, _) = permeability_measure(case);
    let agreement = curve_agreement(sangani_series_s, &[(phi, s_measured)], 0.10, 0.0);
    assert!(
        re_a < 0.05,
        "VAL EXACT SANGANI light Re_a={re_a:.9e}, band<0.05, normalization=mean_u*a_actual/nu"
    );
    assert!(
        agreement.max_rel_dev <= 0.10,
        "VAL EXACT SANGANI light S=4pi/F* ratio max_rel_dev={:.9e}, band=0.10, denominator=Sangani series S(phi_actual), phi_actual={phi:.9e}, S_measured={s_measured:.9e}",
        agreement.max_rel_dev
    );
}

#[test]
#[ignore = "heavy VAL-EXACT Sangani-Acrivos phi={0.05,0.10,0.20} permeability sweep"]
fn g3_sangani_acrivos_heavy_phi_sweep() {
    let mut samples = Vec::new();
    let mut ks = Vec::new();
    for phi in [0.05, 0.10, 0.20] {
        let (phi_actual, k, s_measured, s_ref, re_a, _) =
            permeability_measure(case_for_phi(phi, false));
        println!(
            "VAL EXACT SANGANI HEAVY row: phi_target={phi:.6} phi_actual={phi_actual:.9e} k={k:.9e} S_measured={s_measured:.9e} S_ref={s_ref:.9e} Re_a={re_a:.9e}"
        );
        assert!(
            re_a < 0.05,
            "VAL EXACT SANGANI heavy phi={phi} Re_a={re_a:.9e}, band<0.05, normalization=mean_u*a_actual/nu"
        );
        samples.push((phi_actual, s_measured));
        ks.push(k);
    }
    let agreement = curve_agreement(sangani_series_s, &samples, 0.10, 0.0);
    let k_monotone = monotonicity(&ks);
    println!(
        "VAL EXACT SANGANI HEAVY sweep: S_max_rel_dev={:.9e} worst_phi={:.9e} frac_in_10pct={:.9e} k_monotonicity_decreasing={k_monotone:.9e} ks={ks:?}",
        agreement.max_rel_dev, agreement.worst_x, agreement.frac_in_band
    );
    assert!(
        agreement.max_rel_dev <= 0.10,
        "VAL EXACT SANGANI heavy S=4pi/F* ratio max_rel_dev={:.9e}, band=0.10, denominator=Sangani series S(phi_actual)",
        agreement.max_rel_dev
    );
    assert!(
        k_monotone >= 1.0,
        "VAL EXACT SANGANI heavy permeability monotonicity={k_monotone:.9e}, band=1.0, sequence=k(phi=0.05,0.10,0.20) must strictly decrease"
    );
}

#[test]
#[ignore = "heavy VAL-EXACT Sangani-Acrivos A5 Bouzidi phi=0.10 cross-check"]
fn g3_sangani_acrivos_heavy_bouzidi_cross_check() {
    // Compat exposes set_bouzidi_circle(cx,cy,r), so this is an executable
    // cross-check rather than a SPEC-GAP stub.  Both runs use the same
    // rasterized solid-cell mask for phi_actual; Bouzidi only changes the
    // fluid-solid link intersection distances.  The physical expectation is
    // that the curved-wall result moves toward the Sangani-Acrivos series
    // relative to the staircase boundary.
    let stair = permeability_measure(case_for_phi(0.10, false));
    let bouzidi = permeability_measure(case_for_phi(0.10, true));
    let err_stair = (stair.2 - stair.3).abs() / stair.3.abs();
    let err_bouzidi = (bouzidi.2 - bouzidi.3).abs() / bouzidi.3.abs();
    println!(
        "VAL EXACT SANGANI BOUZIDI: phi_stair={:.9e} S_stair={:.9e} S_ref={:.9e} err_stair={err_stair:.9e} phi_bouzidi={:.9e} S_bouzidi={:.9e} S_ref_bouzidi={:.9e} err_bouzidi={err_bouzidi:.9e}",
        stair.0, stair.2, stair.3, bouzidi.0, bouzidi.2, bouzidi.3
    );
    assert!(
        err_bouzidi <= err_stair,
        "VAL EXACT SANGANI Bouzidi cross-check err_bouzidi={err_bouzidi:.9e}, err_stair={err_stair:.9e}, denominator=Sangani S(phi_actual); Bouzidi should approach/bracket the series value"
    );
}
