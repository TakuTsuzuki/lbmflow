//! T17/VR-STR-03 — Re_tau=178.12 turbulent channel WALE characterization.

use lbm_core::prelude::*;
use std::f64::consts::PI;

type ChanSolver<B = CpuSimd> = Solver<D3Q19, f64, B, LocalPeriodic>;

const TRT_MAGIC: f64 = CollisionKind::MAGIC_STD;
const RE_TAU_MKM180: f64 = 178.12;

// Downsampled from docs/reference/mkm180/chan180.means, columns y+ and Umean.
// Source: Moser, Kim & Mansour channel DNS, Re_tau = 178.12,
// http://www.tam.uiuc.edu/Faculty/Moser/channel (chan180/profiles/chan180.means).
const MKM180_MEAN_U_PLUS: &[(f64, f64)] = &[
    (2.6224, 2.5939),
    (5.3381, 5.1133),
    (8.9902, 7.9052),
    (13.5590, 10.3110),
    (19.0190, 12.0660),
    (25.3420, 13.2760),
    (30.0190, 13.8700),
    (37.7010, 14.5570),
    (46.1430, 15.1010),
    (55.3000, 15.5690),
    (65.1230, 15.9930),
    (75.5590, 16.3890),
    (86.5500, 16.7560),
    (98.0370, 17.0940),
    (109.9600, 17.4000),
    (122.2500, 17.6720),
    (134.8400, 17.9110),
    (147.6700, 18.1030),
    (160.6600, 18.2350),
    (173.7500, 18.2970),
];

#[derive(Clone, Copy, Debug)]
struct ChannelCase {
    delta: usize,
    nx: usize,
    ny: usize,
    nz: usize,
    u_tau: f64,
    nu: f64,
    force_x: f64,
}

#[derive(Clone, Copy, Debug)]
struct Protocol {
    warmup_steps: usize,
    stats_steps: usize,
    sample_every: usize,
    smoke: bool,
}

#[derive(Clone, Debug)]
struct Stats {
    samples: usize,
    mean_u: Vec<f64>,
    mean_v: Vec<f64>,
    mean_uv: Vec<f64>,
    last_window_uv_plus_at_30: f64,
}

#[derive(Clone, Debug)]
struct Report {
    case: ChannelCase,
    protocol: Protocol,
    achieved_re_tau: f64,
    centerline_u_plus: f64,
    mean_profile_l2rel: f64,
    total_stress_l2rel: f64,
    uv_plus_peak: f64,
    uv_plus_at_30_last_window: f64,
}

fn full_case() -> ChannelCase {
    let delta = 48usize;
    let u_tau = 0.008;
    let nx = 128usize; // Lx+ = 474.99, Lx/delta = 2.67.
    let nz = 72usize; // Lz+ = 267.18, Lz/delta = 1.50.
    let ny = 2 * delta + 2;
    let nu = u_tau * delta as f64 / RE_TAU_MKM180;
    let force_x = u_tau * u_tau / delta as f64;
    ChannelCase {
        delta,
        nx,
        ny,
        nz,
        u_tau,
        nu,
        force_x,
    }
}

fn smoke_case() -> ChannelCase {
    let delta = 24usize;
    let u_tau = 0.006;
    let nx = 32usize;
    let nz = 16usize;
    let ny = 2 * delta + 2;
    let nu = u_tau * delta as f64 / RE_TAU_MKM180;
    let force_x = u_tau * u_tau / delta as f64;
    ChannelCase {
        delta,
        nx,
        ny,
        nz,
        u_tau,
        nu,
        force_x,
    }
}

fn full_protocol(case: ChannelCase) -> Protocol {
    let eddy_turnover = case.delta as f64 / case.u_tau;
    Protocol {
        warmup_steps: (20.0 * eddy_turnover).ceil() as usize,
        stats_steps: (30.0 * eddy_turnover).ceil() as usize,
        sample_every: 40,
        smoke: false,
    }
}

fn smoke_protocol(case: ChannelCase) -> Protocol {
    let eddy_turnover = case.delta as f64 / case.u_tau;
    Protocol {
        warmup_steps: (2.0 * eddy_turnover).ceil() as usize,
        stats_steps: (2.0 * eddy_turnover).ceil() as usize,
        sample_every: 50,
        smoke: true,
    }
}

fn laminar_smoke_case() -> ChannelCase {
    let delta = 8usize;
    let u_tau = 0.004;
    let nx = 8usize;
    let nz = 4usize;
    let ny = 2 * delta + 2;
    let nu = 0.08;
    let force_x = u_tau * u_tau / delta as f64;
    ChannelCase {
        delta,
        nx,
        ny,
        nz,
        u_tau,
        nu,
        force_x,
    }
}

fn make_channel<B>(case: ChannelCase, backend: B) -> ChanSolver<B>
where
    B: Backend<D3Q19, f64, Fields = SoaFields<f64>>,
{
    let spec = GlobalSpec {
        dims: [case.nx, case.ny, case.nz],
        nu: case.nu,
        collision: CollisionKind::Trt { magic: TRT_MAGIC },
        periodic: [true, false, true],
        force: [case.force_x, 0.0, 0.0],
        ..Default::default()
    };
    let mut walls = WallSpec::default();
    walls.is_wall[Face::YNeg.index()] = true;
    walls.is_wall[Face::YPos.index()] = true;
    let (solid, wall_u) = build_wall_rims(D3Q19::D, spec.dims, &walls);
    Solver::new(&spec, &solid, &wall_u, [1, 1, 1], backend, LocalPeriodic)
}

fn yw_from_y(y: usize, ny: usize) -> f64 {
    let lower = y as f64 - 0.5;
    let upper = ny as f64 - 1.5 - y as f64;
    lower.min(upper)
}

fn init_turbulent_channel<B>(solver: &mut ChanSolver<B>, case: ChannelCase)
where
    B: Backend<D3Q19, f64, Fields = SoaFields<f64>>,
{
    let h = (case.ny - 2) as f64;
    let amp = 0.3 * case.u_tau;
    let u_max = 18.3 * case.u_tau;
    let nx = case.nx as f64;
    let nz = case.nz as f64;
    solver.init_with(move |x, y, z| {
        if y == 0 || y + 1 == case.ny {
            return (1.0, [0.0; 3]);
        }
        let yw = y as f64 - 0.5;
        let eta = yw / h;
        let wall_taper = (4.0 * eta * (1.0 - eta)).max(0.0);
        let base_u = 4.0 * u_max * eta * (1.0 - eta);
        let ax = 2.0 * PI * x as f64 / nx;
        let az = 2.0 * PI * z as f64 / nz;
        let ay = PI * eta;
        let ux =
            base_u + amp * wall_taper * (ax.sin() * (2.0 * az).cos() + 0.5 * (3.0 * ax + az).sin());
        let uy = amp
            * wall_taper
            * (ay.sin() * ax.cos() * az.sin() + 0.35 * (2.0 * ax - az).cos() * (2.0 * ay).sin());
        let uz = amp
            * wall_taper
            * (ay.sin() * ax.sin() * az.cos() - 0.25 * (ax + 2.0 * az).cos() * ay.sin());
        (1.0, [ux, uy, uz])
    });
}

fn sample_planes<B>(solver: &ChanSolver<B>, case: ChannelCase) -> (Vec<f64>, Vec<f64>, Vec<f64>)
where
    B: Backend<D3Q19, f64, Fields = SoaFields<f64>>,
{
    let ux = solver.gather_ux();
    let uy = solver.gather_uy();
    let plane_n = (case.nx * case.nz) as f64;
    let mut u = vec![0.0; case.ny];
    let mut v = vec![0.0; case.ny];
    let mut uv = vec![0.0; case.ny];
    for z in 0..case.nz {
        for y in 0..case.ny {
            for x in 0..case.nx {
                let i = z * case.nx * case.ny + y * case.nx + x;
                u[y] += ux[i];
                v[y] += uy[i];
                uv[y] += ux[i] * uy[i];
            }
        }
    }
    for y in 0..case.ny {
        u[y] /= plane_n;
        v[y] /= plane_n;
        uv[y] /= plane_n;
    }
    (u, v, uv)
}

fn add_sample(acc: &mut Stats, u: &[f64], v: &[f64], uv: &[f64]) {
    acc.samples += 1;
    for y in 0..u.len() {
        acc.mean_u[y] += u[y];
        acc.mean_v[y] += v[y];
        acc.mean_uv[y] += uv[y];
    }
}

fn finish_stats(mut stats: Stats) -> Stats {
    let inv = 1.0 / stats.samples as f64;
    for y in 0..stats.mean_u.len() {
        stats.mean_u[y] *= inv;
        stats.mean_v[y] *= inv;
        stats.mean_uv[y] *= inv;
    }
    stats
}

fn interp(xs: &[f64], ys: &[f64], x: f64) -> f64 {
    assert!(xs.len() == ys.len() && xs.len() >= 2);
    if x <= xs[0] {
        return ys[0];
    }
    for i in 0..xs.len() - 1 {
        if x <= xs[i + 1] {
            let t = (x - xs[i]) / (xs[i + 1] - xs[i]);
            return ys[i] * (1.0 - t) + ys[i + 1] * t;
        }
    }
    ys[ys.len() - 1]
}

fn folded_profile(case: ChannelCase, values: &[f64]) -> (Vec<f64>, Vec<f64>) {
    let mut y_plus = Vec::new();
    let mut out = Vec::new();
    for y in 1..=case.delta {
        let ym = case.ny - 1 - y;
        y_plus.push((y as f64 - 0.5) * case.u_tau / case.nu);
        out.push(0.5 * (values[y] + values[ym]));
    }
    (y_plus, out)
}

fn l2rel_mean_profile(case: ChannelCase, mean_u: &[f64]) -> f64 {
    let (y_plus, u_folded) = folded_profile(case, mean_u);
    let u_plus: Vec<f64> = u_folded.iter().map(|u| u / case.u_tau).collect();
    let mut err2 = 0.0;
    let mut ref2 = 0.0;
    for &(yp, dns_u) in MKM180_MEAN_U_PLUS {
        if (5.0..=150.0).contains(&yp) {
            let sim_u = interp(&y_plus, &u_plus, yp);
            err2 += (sim_u - dns_u).powi(2);
            ref2 += dns_u.powi(2);
        }
    }
    (err2 / ref2).sqrt()
}

fn reynolds_uv_plus(case: ChannelCase, stats: &Stats, y: usize) -> f64 {
    -(stats.mean_uv[y] - stats.mean_u[y] * stats.mean_v[y]) / (case.u_tau * case.u_tau)
}

fn uv_plus_peak(case: ChannelCase, stats: &Stats) -> f64 {
    (1..case.ny - 1)
        .map(|y| reynolds_uv_plus(case, stats, y))
        .fold(f64::NEG_INFINITY, f64::max)
}

fn y_index_near_plus(case: ChannelCase, target_y_plus: f64) -> usize {
    (1..case.ny - 1)
        .min_by(|&a, &b| {
            let ya = yw_from_y(a, case.ny) * case.u_tau / case.nu;
            let yb = yw_from_y(b, case.ny) * case.u_tau / case.nu;
            (ya - target_y_plus)
                .abs()
                .total_cmp(&(yb - target_y_plus).abs())
        })
        .unwrap()
}

fn total_stress_l2rel(case: ChannelCase, stats: &Stats) -> f64 {
    let mut err2 = 0.0;
    let mut ref2 = 0.0;
    for y in 2..case.ny - 2 {
        let y_delta = yw_from_y(y, case.ny) / case.delta as f64;
        if !(0.2..=0.8).contains(&y_delta) {
            continue;
        }
        let lower_side = y <= case.delta;
        let sign = if lower_side { 1.0 } else { -1.0 };
        let du_dy = (stats.mean_u[y + 1] - stats.mean_u[y - 1]) / 2.0;
        let viscous = case.nu * sign * du_dy;
        let resolved = reynolds_uv_plus(case, stats, y) * case.u_tau * case.u_tau;
        let total = viscous + resolved;
        let analytic = case.u_tau * case.u_tau * (1.0 - y_delta);
        err2 += (total - analytic).powi(2);
        ref2 += analytic.powi(2);
    }
    (err2 / ref2).sqrt()
}

fn run_channel<B>(mut solver: ChanSolver<B>, case: ChannelCase, protocol: Protocol) -> Report
where
    B: Backend<D3Q19, f64, Fields = SoaFields<f64>>,
{
    init_turbulent_channel(&mut solver, case);
    let mut les = WaleLes::new();
    for _ in 0..protocol.warmup_steps {
        les.update(&mut solver);
        solver.run(1);
    }

    let mut stats = Stats {
        samples: 0,
        mean_u: vec![0.0; case.ny],
        mean_v: vec![0.0; case.ny],
        mean_uv: vec![0.0; case.ny],
        last_window_uv_plus_at_30: 0.0,
    };
    let mut last = stats.clone();
    let last_window_steps =
        ((10.0 * case.delta as f64 / case.u_tau).ceil() as usize).min(protocol.stats_steps);
    let last_window_start = protocol.stats_steps.saturating_sub(last_window_steps);
    let y30 = y_index_near_plus(case, 30.0);

    for step in 0..protocol.stats_steps {
        les.update(&mut solver);
        solver.run(1);
        if (step + 1) % protocol.sample_every == 0 || step + 1 == protocol.stats_steps {
            let (u, v, uv) = sample_planes(&solver, case);
            add_sample(&mut stats, &u, &v, &uv);
            if step + 1 >= last_window_start {
                add_sample(&mut last, &u, &v, &uv);
            }
        }
    }

    let mut stats = finish_stats(stats);
    if last.samples > 0 {
        last = finish_stats(last);
        stats.last_window_uv_plus_at_30 = reynolds_uv_plus(case, &last, y30);
    }
    let mean_profile_l2rel = l2rel_mean_profile(case, &stats.mean_u);
    let total_stress_l2rel = total_stress_l2rel(case, &stats);
    let uv_plus_peak = uv_plus_peak(case, &stats);
    let centerline_u_plus =
        0.5 * (stats.mean_u[case.delta] + stats.mean_u[case.delta + 1]) / case.u_tau;
    Report {
        case,
        protocol,
        achieved_re_tau: case.u_tau * case.delta as f64 / case.nu,
        centerline_u_plus,
        mean_profile_l2rel,
        total_stress_l2rel,
        uv_plus_peak,
        uv_plus_at_30_last_window: stats.last_window_uv_plus_at_30,
    }
}

fn print_report(label: &str, r: &Report) {
    eprintln!(
        "{label}: dims={}x{}x{}, delta={}, u_tau={:.6}, nu={:.8e}, force_x={:.8e}, \
         warmup={}, stats={}, sample_every={}, samples~{}, achieved_Re_tau={:.4}, \
         centerline_U+={:.4}, mean_U+_L2rel={:.4}, total_stress_L2rel={:.4}, \
         peak_-u'v'+={:.4}, last10Te_-u'v'+(y+~30)={:.4}",
        r.case.nx,
        r.case.ny,
        r.case.nz,
        r.case.delta,
        r.case.u_tau,
        r.case.nu,
        r.case.force_x,
        r.protocol.warmup_steps,
        r.protocol.stats_steps,
        r.protocol.sample_every,
        r.protocol.stats_steps.div_ceil(r.protocol.sample_every),
        r.achieved_re_tau,
        r.centerline_u_plus,
        r.mean_profile_l2rel,
        r.total_stress_l2rel,
        r.uv_plus_peak,
        r.uv_plus_at_30_last_window,
    );
}

#[test]
fn wale_channel_laminar_harness_smoke() {
    let case = laminar_smoke_case();
    let mut solver = make_channel(case, CpuSimd::default());
    solver.init_with(|_, _, _| (1.0, [0.0; 3]));
    let mut les = WaleLes::new();
    let steps = 12_000usize;
    for _ in 0..steps {
        les.update(&mut solver);
        solver.run(1);
    }
    let (u, _, _) = sample_planes(&solver, case);
    let h = (case.ny - 2) as f64;
    let mut err2 = 0.0;
    let mut ref2 = 0.0;
    for y in 1..case.ny - 1 {
        let yw = y as f64 - 0.5;
        let want = case.force_x * yw * (h - yw) / (2.0 * case.nu);
        err2 += (u[y] - want).powi(2);
        ref2 += want.powi(2);
    }
    let l2rel = (err2 / ref2).sqrt();
    let max_nu_t = les.nu_t().iter().copied().fold(0.0_f64, f64::max);
    eprintln!(
        "laminar channel harness smoke: steps={steps}, L2rel={l2rel:.6e}, max_nu_t={max_nu_t:.6e}"
    );
    assert!(
        l2rel <= 0.02,
        "laminar WALE channel profile L2rel={l2rel:.6e} > 2.0e-2 after {steps} steps"
    );
    assert!(
        max_nu_t <= 1.0e-12,
        "WALE must stay null in laminar channel, max nu_t={max_nu_t:.6e}"
    );
}

#[test]
#[ignore = "T17/VR-STR-03 heavy: ~40-60 min CPU"]
fn channel_re_tau_180_wale_vs_mkm_dns() {
    let smoke = std::env::var_os("LBM_CHAN180_SMOKE").is_some();
    let case = if smoke { smoke_case() } else { full_case() };
    let protocol = if smoke {
        smoke_protocol(case)
    } else {
        full_protocol(case)
    };
    let solver = make_channel(case, CpuSimd::default());
    let report = run_channel(solver, case, protocol);
    if report.protocol.smoke {
        print_report(
            "NON-PHYSICAL SMOKE channel Re_tau=178.12 WALE vs MKM DNS harness",
            &report,
        );
        return;
    }

    print_report(
        "channel Re_tau=178.12 WALE vs MKM DNS characterization",
        &report,
    );
    let mean_band = 0.15; // BAND-FREEZE-PENDING(PM)
    assert!(
        report.mean_profile_l2rel <= mean_band,
        "mean U+ L2rel={:.6} > {:.6} vs MKM180 DNS",
        report.mean_profile_l2rel,
        mean_band
    );
    let stress_band = 0.10; // BAND-FREEZE-PENDING(PM)
    assert!(
        report.total_stress_l2rel <= stress_band,
        "total-stress L2rel={:.6} > {:.6} vs analytic force-balance line",
        report.total_stress_l2rel,
        stress_band
    );
    assert!(
        report.uv_plus_at_30_last_window > 0.4,
        "sustained turbulence guard failed: last-10Te -<u'v'>+ at y+~30 = {:.6} <= 0.4",
        report.uv_plus_at_30_last_window
    );
}
