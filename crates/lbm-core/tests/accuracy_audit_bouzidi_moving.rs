//! ACC-AUDIT lane 1.7: Bouzidi moving-wall coverage.
//!
//! API discovery result: native Bouzidi records support moving walls. The
//! native `Solver::set_bouzidi_links` records carry `wall_ref`, and
//! `crates/lbm-core/src/bouzidi.rs::apply_bouzidi_impl` reads
//! `fields.wall_u[wall_ref]` before applying the interpolated bounce-back
//! replacement. The compat facade exposes `set_bouzidi_half_way_links` and
//! `set_bouzidi_circle`, but arbitrary off-grid straight-wall qd placement is
//! a native validation hook, so these tests use the native D2Q9 CPU scalar
//! solver directly.
//!
//! Moving-wall derivation, following Bouzidi, Firdaouss, and Lallemand 2001
//! Section 4. For a wall moving with velocity u_w, bounce-back reflects the
//! nonequilibrium population while the equilibrium shift injects the wall
//! momentum. The half-way moving-wall link replacement is
//!
//!   f_opp(q)(x_f, t+1) = f_q(x_f, t) + 2 w_q rho (c_opp(q) . u_w) / c_s^2.
//!
//! D2Q9 has c_s^2 = 1/3, so the correction is
//!
//!   6 w_q rho (c_opp(q) . u_w).
//!
//! Bouzidi interpolation replaces the stationary reflected value by a linear
//! interpolation along the link. The same moving-wall momentum source must be
//! weighted by the interpolation coefficient sigma_i:
//!
//!   qd < 1/2:  f_opp(q) = 2 qd f_q(x_f)
//!                     + (1 - 2 qd) f_q(x_f - c_q)
//!                     + 2 qd * 2 w_q rho (c_opp(q) . u_w) / c_s^2
//!   qd >= 1/2: f_opp(q) = (1 / 2 qd) f_q(x_f)
//!                     + (1 - 1 / 2 qd) f_opp(q)(x_f)
//!                     + (1 / 2 qd) * 2 w_q rho (c_opp(q) . u_w) / c_s^2.
//!
//! Therefore sigma_i is {2 qd, 1/(2 qd)} and qd=1/2 must degenerate exactly
//! to the existing half-way moving-wall rule.
//!
//! Current triage note, 2026-07-07: the qd < 1/2 branch is supported by the
//! API but does not satisfy the off-grid Couette reference below. It currently
//! imposes about sigma_i * U_wall (0.5 U_wall at qd=0.25). The exact-reference
//! test is landed as ignored with ANOM-L1_7-001, and the default suite carries
//! an active current-wrong pin so the eventual fix fails this file until the
//! pin is retightened.

mod common;

use common::metrics::curve_agreement;
use lbm_core::prelude::*;

const NX: usize = 48;
const NY: usize = 34;
const TAU: f64 = 0.8;
const NU: f64 = (TAU - 0.5) / 3.0;
const U_WALL: f64 = 0.05;
const TRT_MAGIC: f64 = 3.0 / 16.0;

type Sim = Solver<D2Q9, f64, CpuScalar, LocalPeriodic>;

fn couette_wall_positions(qd: f64) -> (f64, f64) {
    let wall_lo = 1.0 - qd;
    let wall_hi = (NY - 2) as f64 + qd;
    (wall_lo, wall_hi)
}

fn build_couette(qd: Option<f64>) -> Sim {
    let mut walls = WallSpec::default();
    walls.is_wall[Face::YNeg.index()] = true;
    walls.is_wall[Face::YPos.index()] = true;
    walls.u[Face::YPos.index()] = [U_WALL, 0.0, 0.0];

    let spec = GlobalSpec {
        dims: [NX, NY, 1],
        nu: NU,
        periodic: [true, false, false],
        collision: CollisionKind::Trt { magic: TRT_MAGIC },
        ..Default::default()
    };
    let (solid, wall_u) = build_wall_rims(2, spec.dims, &walls);
    let mut sim = Solver::new(
        &spec,
        &solid,
        &wall_u,
        [1, 1, 1],
        CpuScalar::default(),
        LocalPeriodic,
    );

    let (wall_lo, wall_hi) = qd
        .map(couette_wall_positions)
        .unwrap_or((0.5, (NY - 1) as f64 - 0.5));
    let h = wall_hi - wall_lo;
    sim.init_with(|_, y, _| {
        if y == 0 || y == NY - 1 {
            (1.0, [0.0, 0.0, 0.0])
        } else {
            let y_w = y as f64 - wall_lo;
            (1.0, [U_WALL * y_w / h, 0.0, 0.0])
        }
    });

    if let Some(qd) = qd {
        install_straight_wall_bouzidi_links(&mut sim, qd);
    }

    sim
}

fn install_straight_wall_bouzidi_links(sim: &mut Sim, qd: f64) {
    assert!(qd > 0.0 && qd < 1.0, "Bouzidi qd must lie inside (0,1)");
    let g = sim.sub(0).geom;
    let mut records = Vec::new();
    for x in 0..NX {
        let bottom = g.pidx(x, 1, 0);
        let top = g.pidx(x, NY - 2, 0);
        for q in [4usize, 7, 8] {
            let c = D2Q9::C[q];
            let wall_ref = g.pidx_i(x as isize + c[0] as isize, 0, 0);
            records.push(BouzidiLink {
                cell: bottom as u32,
                q: q as u8,
                qd,
                has_second: true,
                wall_ref: wall_ref as u32,
            });
        }
        for q in [2usize, 5, 6] {
            let c = D2Q9::C[q];
            let wall_ref = g.pidx_i(x as isize + c[0] as isize, (NY - 1) as isize, 0);
            records.push(BouzidiLink {
                cell: top as u32,
                q: q as u8,
                qd,
                has_second: true,
                wall_ref: wall_ref as u32,
            });
        }
    }
    sim.set_bouzidi_links(0, Some(BouzidiLinks::new(records)));
}

fn mean_ux_profile(sim: &Sim) -> Vec<f64> {
    let mut profile = Vec::with_capacity(NY - 2);
    for y in 1..NY - 1 {
        let sum = (0..NX).map(|x| sim.u(x, y, 0)[0]).sum::<f64>();
        profile.push(sum / NX as f64);
    }
    profile
}

fn assert_monotone_profile(label: &str, profile: &[f64]) {
    for (i, pair) in profile.windows(2).enumerate() {
        assert!(
            pair[1] >= pair[0],
            "{label}: Couette behavior anchor failed at interior rows {}->{}: ux decreased from {:.12e} to {:.12e}",
            i + 1,
            i + 2,
            pair[0],
            pair[1]
        );
    }
}

fn couette_curve_agreement(
    qd: f64,
    scale: f64,
    rel_band: f64,
) -> (common::metrics::CurveAgreement, Vec<f64>) {
    let mut sim = build_couette(Some(qd));
    sim.run(4000);
    let profile = mean_ux_profile(&sim);
    let (wall_lo, wall_hi) = couette_wall_positions(qd);
    let h = wall_hi - wall_lo;
    let samples: Vec<(f64, f64)> = profile
        .iter()
        .enumerate()
        .map(|(i, &ux)| {
            let y = i + 1;
            let y_w = y as f64 - wall_lo;
            (y_w, ux)
        })
        .collect();
    (
        curve_agreement(|y_w| scale * U_WALL * y_w / h, &samples, rel_band, U_WALL),
        profile,
    )
}

fn assert_couette_curve(label: &str, qd: f64, scale: f64, rel_band: f64) {
    let (agreement, profile) = couette_curve_agreement(qd, scale, rel_band);
    let (wall_lo, wall_hi) = couette_wall_positions(qd);
    let h = wall_hi - wall_lo;
    let first = profile[0];
    let mid = profile[profile.len() / 2];
    let last = profile[profile.len() - 1];
    println!(
        "{label}: qd={qd:.2}, H={h:.6}, scale={scale:.6}, ux_first={first:.12e}, ux_mid={mid:.12e}, ux_last={last:.12e}, max_rel_dev={:.12e}, worst_y_w={:.6}, frac_in_band={:.3}, band={rel_band:.12e}",
        agreement.max_rel_dev, agreement.worst_x, agreement.frac_in_band
    );
    assert!(
        agreement.max_rel_dev <= rel_band,
        "{label}: qd={qd:.2} curve_agreement max_rel_dev={:.12e}, band={rel_band:.12e}, denominator=max(|scale*U_wall*y_w/H|, U_wall={U_WALL:.12e}); scale={scale:.6}",
        agreement.max_rel_dev,
    );
    assert!(
        agreement.frac_in_band == 1.0,
        "{label}: qd={qd:.2} curve_agreement frac_in_band={:.6}, expected 1.0; denominator=max(|scale*U_wall*y_w/H|, U_wall={U_WALL:.12e}); scale={scale:.6}",
        agreement.frac_in_band
    );
    assert_monotone_profile(label, &profile);
}

#[test]
fn qd_half_moving_wall_is_bitwise_half_way_moving_wall() {
    let mut half = build_couette(None);
    let mut bouzidi = build_couette(Some(0.5));

    for step in 0..40 {
        half.step();
        bouzidi.step();
        assert_eq!(
            half.gather_ux(),
            bouzidi.gather_ux(),
            "lane 1.7 qd=0.5 moving-wall ux bit-match failed after step {step}; U_wall={U_WALL}, tau={TAU}"
        );
        assert_eq!(
            half.gather_uy(),
            bouzidi.gather_uy(),
            "lane 1.7 qd=0.5 moving-wall uy bit-match failed after step {step}; U_wall={U_WALL}, tau={TAU}"
        );
        assert_eq!(
            half.gather_rho(),
            bouzidi.gather_rho(),
            "lane 1.7 qd=0.5 moving-wall rho bit-match failed after step {step}; U_wall={U_WALL}, tau={TAU}"
        );
        for q in 0..D2Q9::Q {
            assert_eq!(
                half.gather_f(q),
                bouzidi.gather_f(q),
                "lane 1.7 qd=0.5 moving-wall f[{q}] bit-match failed after step {step}; U_wall={U_WALL}, tau={TAU}"
            );
        }
    }

    let profile = mean_ux_profile(&bouzidi);
    assert_monotone_profile("lane 1.7 qd=0.5 moving wall", &profile);
}

#[test]
fn qd_sweep_moving_wall_couette_matches_offgrid_linear_profile_for_supported_branch() {
    for qd in [0.5, 0.75] {
        assert_couette_curve(
            "lane 1.7 Bouzidi moving Couette exact branch",
            qd,
            1.0,
            2.0e-3,
        );
    }
}

#[test]
fn qd_lt_half_current_wrong_pin_imposes_sigma_scaled_wall_speed() {
    let qd = 0.25;
    let sigma = 2.0 * qd;
    assert_couette_curve(
        "ANOM-L1_7-001 current wrong qd<0.5 Bouzidi moving-wall pin",
        qd,
        sigma,
        1.0e-2,
    );
}

#[test]
#[ignore = "ANOM-L1_7-001: qd<0.5 currently imposes about sigma_i*U_wall; enable after the Bouzidi moving-wall fix and delete/retighten the current-wrong pin"]
fn qd_sweep_moving_wall_couette_should_match_offgrid_linear_profile_all_qd() {
    for qd in [0.25, 0.5, 0.75] {
        assert_couette_curve(
            "lane 1.7 Bouzidi moving Couette exact all-qd",
            qd,
            1.0,
            2.0e-3,
        );
    }
}
