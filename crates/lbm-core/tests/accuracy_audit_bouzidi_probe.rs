//! ACC BOUZIDI PROBE: ANOM-P4-011 kill-case for the Bouzidi force-probe
//! correction sign at `qd != 1/2`.
//!
//! The compat facade exposes `set_bouzidi_circle` and
//! `set_bouzidi_half_way_links`; arbitrary `qd` assignment is intentionally a
//! narrow native validation hook (`Solver::set_bouzidi_links`). These tests use
//! the native D2Q9 CPU scalar solver with the same 2D compat setup requested by
//! the audit: 64x32, TRT tau=0.8, periodic left/right, bounce-back top/bottom,
//! and a staircase disk of radius 6 centered at (32,16). The custom records are
//! installed only on disk links; the straight top/bottom walls stay on the
//! ordinary half-way bounce-back path.

use lbm_core::lattice::D2Q9;
use lbm_core::prelude::*;
use std::f64::consts::PI;

const NX: usize = 64;
const NY: usize = 32;
const CX: f64 = 32.0;
const CY: f64 = 16.0;
const R: f64 = 6.0;
const TAU: f64 = 0.8;
const NU: f64 = (TAU - 0.5) / 3.0;
const FORCE_X: f64 = 2.0e-6;
const TRT_MAGIC: f64 = 3.0 / 16.0;

type Sim = Solver<D2Q9, f64, CpuScalar, LocalPeriodic>;

#[derive(Clone, Copy)]
enum Probe {
    AllSolids,
    DiskOnly,
}

#[derive(Clone, Copy)]
enum Qd {
    Uniform(f64),
    Mixed,
}

fn disk_at(cx: f64, x: usize, y: usize) -> bool {
    let dx = x as f64 - cx;
    let dy = y as f64 - CY;
    dx * dx + dy * dy <= R * R
}

fn disk_at_i(cx: f64, x: isize, y: isize) -> bool {
    x >= 0
        && y >= 0
        && (x as usize) < NX
        && (y as usize) < NY
        && disk_at(cx, x as usize, y as usize)
}

fn spec(force_x: f64) -> (GlobalSpec<f64>, WallSpec<f64>) {
    let mut walls = WallSpec::default();
    walls.is_wall[Face::YNeg.index()] = true;
    walls.is_wall[Face::YPos.index()] = true;
    (
        GlobalSpec {
            dims: [NX, NY, 1],
            nu: NU,
            periodic: [true, false, false],
            force: [force_x, 0.0, 0.0],
            collision: CollisionKind::Trt { magic: TRT_MAGIC },
            ..Default::default()
        },
        walls,
    )
}

fn make_sim(cx: f64, force_x: f64, qd: Qd, probe: Probe, mirror_init: bool) -> Sim {
    let (spec, walls) = spec(force_x);
    let (solid, wall_u) = build_wall_rims(2, spec.dims, &walls);
    let mut sim = Solver::new(
        &spec,
        &solid,
        &wall_u,
        [1, 1, 1],
        CpuScalar::default(),
        LocalPeriodic,
    );
    for y in 0..NY {
        for x in 0..NX {
            if disk_at(cx, x, y) {
                sim.set_solid(x, y, 0);
            }
        }
    }
    match probe {
        Probe::AllSolids => sim.set_force_probe(|_, _, _| true),
        Probe::DiskOnly => sim.set_force_probe(move |x, y, _| disk_at(cx, x, y)),
    }
    install_disk_bouzidi_links(&mut sim, cx, qd);
    sim.init_with(|x, y, _| {
        if y == 0 || y == NY - 1 || disk_at(cx, x, y) {
            (1.0, [0.0, 0.0, 0.0])
        } else {
            let sx = if mirror_init { -1.0 } else { 1.0 };
            let ux_reflect_x = if mirror_init { NX - 1 - x } else { x };
            let ux = sx * 0.012 * (2.0 * PI * y as f64 / NY as f64).sin();
            let uy = 0.004 * (2.0 * PI * ux_reflect_x as f64 / NX as f64).sin();
            (1.0, [ux, uy, 0.0])
        }
    });
    sim
}

fn install_disk_bouzidi_links(sim: &mut Sim, cx: f64, qd: Qd) {
    let g = sim.sub(0).geom;
    let solid = &sim.fields(0).solid;
    let mut records = Vec::new();
    for y in 0..NY {
        for x in 0..NX {
            let cell = g.pidx(x, y, 0);
            if solid[cell] {
                continue;
            }
            for q in 1..D2Q9::Q {
                let c = D2Q9::C[q];
                let sx = x as isize + c[0] as isize;
                let sy = y as isize + c[1] as isize;
                if !disk_at_i(cx, sx, sy) {
                    continue;
                }
                let bx = x as isize - c[0] as isize;
                let by = y as isize - c[1] as isize;
                let has_second = bx >= 0
                    && by >= 0
                    && bx < NX as isize
                    && by < NY as isize
                    && !solid[g.pidx_i(bx, by, 0)];
                let wall_ref = g.pidx_i(sx, sy, 0);
                let qd_value = match qd {
                    Qd::Uniform(v) => v,
                    Qd::Mixed => {
                        let values = [0.2, 0.3, 0.4, 0.6, 0.7, 0.8];
                        values[(records.len() + q) % values.len()]
                    }
                };
                records.push(BouzidiLink {
                    cell: cell as u32,
                    q: q as u8,
                    qd: qd_value,
                    has_second,
                    wall_ref: wall_ref as u32,
                });
            }
        }
    }
    assert!(
        !records.is_empty(),
        "ANOM-P4-011 setup: disk Bouzidi record set must be non-empty"
    );
    if matches!(qd, Qd::Mixed) {
        let mut min_qd = f64::INFINITY;
        let mut max_qd = f64::NEG_INFINITY;
        for r in &records {
            min_qd = min_qd.min(r.qd);
            max_qd = max_qd.max(r.qd);
        }
        assert_eq!(min_qd, 0.2, "mixed qd lower endpoint");
        assert_eq!(max_qd, 0.8, "mixed qd upper endpoint");
    }
    sim.set_bouzidi_links(0, Some(BouzidiLinks::new(records)));
}

fn assert_momentum_ledger(label: &str, sim: &mut Sim, force_x: f64) {
    // Discrete derivation: over one full step, the fluid momentum changes by
    // the uniform Guo body-force impulse on every fluid cell minus the
    // momentum-exchange force reported on the probed solids:
    //
    //   p_fluid(t+1) - p_fluid(t) = N_fluid * F - F_probe.
    //
    // `total_momentum()` includes the Guo half-force physical-velocity shift,
    // but the shift is constant for a uniform body force and cancels in the
    // difference. Because the top and bottom bounce-back walls also absorb
    // x-momentum under body force, these ledger tests probe all solid cells;
    // the curved-wall sign risk remains isolated to the disk records.
    let n_fluid = sim.fluid_cell_count() as f64;
    for step in 0..10 {
        let p0 = sim.total_momentum();
        sim.step();
        let p1 = sim.total_momentum();
        let probe = sim.probed_force();
        let body = [n_fluid * force_x, 0.0, 0.0];
        let dp = [p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]];
        println!(
            "{label}: step={step} dp=({:.12e},{:.12e}) N_fluid*F=({:.12e},{:.12e}) F_probe=({:.12e},{:.12e})",
            dp[0], dp[1], body[0], body[1], probe[0], probe[1]
        );
        for c in 0..2 {
            let expected = body[c] - probe[c];
            let err = (dp[c] - expected).abs();
            let den = (p0[c].abs() + body[c].abs() + probe[c].abs()).max(1.0e-30);
            let band = (1.0e-11 * den).max(5.0e-14);
            assert!(
                err <= band,
                "{label} component {c}: measured dp={:.12e}, expected=N_fluid*F-F_probe={:.12e}, abs_err={:.12e}, band={:.12e}, denominator=|p_t|+|N_fluid*F|+|F_probe|={:.12e}; a Bouzidi sign error would produce an O(1) ledger break",
                dp[c],
                expected,
                err,
                band,
                den
            );
        }
    }
}

#[test]
fn b1_half_qd_global_momentum_ledger_closes() {
    let mut sim = make_sim(CX, FORCE_X, Qd::Uniform(0.5), Probe::AllSolids, false);
    sim.run(5);
    assert_momentum_ledger("ACC BOUZIDI B1 qd=0.5", &mut sim, FORCE_X);
}

#[test]
fn b2_mixed_qd_global_momentum_ledger_closes() {
    let mut sim = make_sim(CX, FORCE_X, Qd::Mixed, Probe::AllSolids, false);
    sim.run(5);
    assert_momentum_ledger("ACC BOUZIDI B2 qd=0.2..0.8", &mut sim, FORCE_X);
}

#[test]
fn b3_half_qd_mirror_obstacle_force_x_mirrors() {
    // Reflection derivation: under x -> (NX-1-x), D2Q9 maps q to its opposite
    // x-direction, the disk center maps from x=32 to x=31, and the driving
    // force plus initial ux change sign. Thus the force on the mirrored disk
    // must be [-Fx, Fy].
    let mut orig = make_sim(CX, FORCE_X, Qd::Uniform(0.5), Probe::DiskOnly, false);
    let mirror_cx = (NX - 1) as f64 - CX;
    let mut mirror = make_sim(mirror_cx, -FORCE_X, Qd::Uniform(0.5), Probe::DiskOnly, true);
    for _ in 0..80 {
        orig.step();
        mirror.step();
    }
    let fo = orig.probed_force();
    let fm = mirror.probed_force();
    let fx_err = (fo[0] + fm[0]).abs();
    let fy_err = (fo[1] - fm[1]).abs();
    let fx_band = (1.0e-12 * fo[0].abs().max(fm[0].abs())).max(5.0e-16);
    let fy_band = (1.0e-12 * fo[1].abs().max(fm[1].abs())).max(5.0e-16);
    println!(
        "ACC BOUZIDI B3: F_orig=({:.12e},{:.12e}) F_mirror=({:.12e},{:.12e}) fx_err={:.12e} fy_err={:.12e}",
        fo[0], fo[1], fm[0], fm[1], fx_err, fy_err
    );
    assert!(
        fx_err <= fx_band,
        "ACC BOUZIDI B3 mirror Fx: |Fx_orig+Fx_mirror|={fx_err:.12e}, band={fx_band:.12e}, denominator=max(|Fx_orig|,|Fx_mirror|,1e-30)"
    );
    assert!(
        fy_err <= fy_band,
        "ACC BOUZIDI B3 mirror Fy: |Fy_orig-Fy_mirror|={fy_err:.12e}, band={fy_band:.12e}, denominator=max(|Fy_orig|,|Fy_mirror|) with abs floor 5e-16"
    );
}

#[test]
fn b4_qd_sweep_drag_keeps_sign_and_changes_smoothly() {
    // A sign bug in the Bouzidi probe correction would flip the measured drag
    // as the interpolation crosses qd=1/2. The exact drag is geometry- and
    // transient-dependent here, so this is a qualitative kill-case: same sign
    // for all qd values and smooth same-order changes around the half-way
    // control. The two off-half samples are the same |qd-0.5| and therefore
    // must remain comparable rather than jumping by an order of magnitude.
    let mut samples = Vec::new();
    for qd in [0.25, 0.5, 0.75] {
        let mut sim = make_sim(CX, FORCE_X, Qd::Uniform(qd), Probe::DiskOnly, false);
        sim.run(500);
        let f = sim.probed_force();
        println!(
            "ACC BOUZIDI B4: qd={qd:.2} F_probe=({:.12e},{:.12e}) |Fx|={:.12e}",
            f[0],
            f[1],
            f[0].abs()
        );
        samples.push((qd, f[0]));
    }

    let sign = samples[0].1.signum();
    assert_ne!(
        sign, 0.0,
        "ACC BOUZIDI B4 qd=0.25 drag is exactly zero; cannot establish sign"
    );
    for (qd, fx) in &samples {
        assert_eq!(
            fx.signum(),
            sign,
            "ACC BOUZIDI B4 sign flip across qd sweep: qd={qd:.2}, Fx={fx:.12e}, baseline sign={sign}"
        );
    }

    let f025 = samples[0].1.abs();
    let f050 = samples[1].1.abs();
    let f075 = samples[2].1.abs();
    assert!(
        f025 >= 0.5 * f050 && f075 >= 0.5 * f050,
        "ACC BOUZIDI B4 monotonic-distance lower guard: |F(0.25)|={f025:.12e}, |F(0.5)|={f050:.12e}, |F(0.75)|={f075:.12e}; off-half samples should not collapse relative to half-way"
    );
    let ratio = f025.max(f075) / f025.min(f075).max(1.0e-30);
    assert!(
        ratio <= 10.0,
        "ACC BOUZIDI B4 smoothness guard at equal |qd-0.5|: |F(0.25)|={f025:.12e}, |F(0.75)|={f075:.12e}, ratio={ratio:.6e}, band=10"
    );
}
