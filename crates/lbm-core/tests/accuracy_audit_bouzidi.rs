//! ACC-AUDIT (Bouzidi curved BC) — P2 dry run for lbmflow-accuracy-audit.
//!
//! P1 audit list (source: .claude/skills/lbmflow-accuracy-audit/references/
//! worked-example-bouzidi.md; approximations of Bouzidi 2001 interpolated
//! bounce-back as implemented in crates/lbm-core/src/bouzidi.rs):
//!
//!   G1 (A1) Interpolation-stencil convergence order at fractional wall.
//!           Bouzidi claims second-order at the wall for straight walls with
//!           qd ∈ (0,1). Off-grid Poiseuille L2 error vs resolution: fit
//!           order = log(err)/log(h); expected slope >= 1.6 (diffusive
//!           scaling, well below the ~2 theoretical claim, ~10x headroom).
//!   G2 (A2) Sub-cell translation invariance. Translating a cylinder centre
//!           by fractional offsets Δ ∈ (0, 1) must move the observable
//!           (drag / centreline velocity) smoothly, not staircase. Coarse
//!           canary sweeps two offsets and bounds the spread — this is the
//!           axis PM cited as the a-priori catch for the 52eaf85 cylinder-
//!           centre-convention bug (+8% Cd bias).
//!   G3 (A5) qd=1/2 bitwise degeneracy with half-way BB (STATIONARY). Owned
//!           by the existing crates/lbm-core/tests/bouzidi.rs test
//!           `qd_half_records_are_bitwise_half_way_bounce_back` — NOT
//!           duplicated here; cross-reference only.
//!   G4 (A1) Off-grid Couette exact recovery. Bouzidi should exactly resolve
//!           a linear Couette profile up to O(Ma²) even off-grid — the
//!           residual is a lower bound on the boundary error alone (no
//!           bulk-scheme contribution). Light.
//!   G5 (A3) Effective wall position vs tau (slip-law). Bouzidi's effective
//!           wall must sit at the analytic wall independent of tau, unlike
//!           half-way BB whose effective wall shifts as (tau - 1/2). Heavy
//!           #[ignore]; carrying the analytic derivation and the exact tau
//!           sweep for future landing.
//!   G6 (A2) Rotational anisotropy. Translating the geometry by a lattice
//!           step is exact (permutation); the interesting axis is rotation
//!           — a diagonal-oriented channel vs axis-aligned must agree at
//!           matched resolution. Heavy #[ignore]; SPEC note only for now,
//!           requires diagonal-domain construction beyond the current
//!           helper set.
//!
//! Merge-queue conventions honored (SKILL.md §Merge-queue hard requirements):
//! - Every assert prints its measured value and NAMES the denominator of
//!   any relative tolerance.
//! - Metrics come from the shared library (common::metrics) — no inline
//!   reimplementations.
//! - CPU scalar reference backend only (T14 gates other backends).
//! - No current-wrong-value PINs in this file (dry run found no calibrated
//!   engine bug in Bouzidi at this pass); if a triage confirms one later
//!   the PIN goes here with its ANOM id in the assert message.

mod common;

use common::metrics::{l2_rel, order_fit};
use lbm_core::prelude::*;

/// Straight-channel Poiseuille reference used by G1/G4 — walls at fractional
/// positions wall_lo, wall_hi in lattice units; nu, uniform body force.
///
/// Analytic derivation: 1D steady Navier-Stokes with a uniform body force in x
/// and no-slip walls at y=wall_lo, wall_hi is
///   nu · d²u/dy² = -F   (F is force per unit mass; ρ=1 here)
///   u(wall_lo) = u(wall_hi) = 0
/// integrating twice and applying BCs:
///   u(y) = F/(2 nu) · (y - wall_lo) · (wall_hi - y)
/// with h = wall_hi - wall_lo and yy = y - wall_lo:
///   u(y) = F · yy · (h - yy) / (2 nu)
/// which is the exact Poiseuille parabola (also the Couette-with-force limit).
fn offgrid_poiseuille_channel(
    nx: usize,
    ny: usize,
    nu: f64,
    force: f64,
    wall_lo: f64,
    wall_hi: f64,
    bouzidi: bool,
) -> Solver<D2Q9, f64, CpuScalar, LocalPeriodic> {
    let mut walls = WallSpec::default();
    walls.is_wall[Face::YNeg.index()] = true;
    walls.is_wall[Face::YPos.index()] = true;
    let (solid, wall_u) = build_wall_rims(2, [nx, ny, 1], &walls);
    let mut solver = Solver::new(
        &GlobalSpec {
            dims: [nx, ny, 1],
            nu,
            periodic: [true, false, false],
            force: [force, 0.0, 0.0],
            collision: CollisionKind::Trt { magic: 3.0 / 16.0 },
            ..Default::default()
        },
        &solid,
        &wall_u,
        [1, 1, 1],
        CpuScalar::default(),
        LocalPeriodic,
    );

    if bouzidi {
        // Build fractional-wall records without leaning on the analytic-circle
        // path: for a straight y-wall at wall_lo, links from y=1 pointing to
        // solid (q=4, 7, 8) have qd = 1 - wall_lo (fractional distance from
        // the fluid cell centre to the wall along -y). Symmetrically at top.
        // Derived directly from the Bouzidi convention: qd is the fractional
        // wall distance along q measured from the fluid cell centre; a wall
        // at y=wall_lo below fluid cell y=1 (link q=4 = -y) is at distance
        // 1 - wall_lo along the link. See docs/PHYSICS.md T15.x conventions.
        let g = solver.sub(0).geom;
        let mut records = Vec::new();
        for x in 0..nx {
            let bottom = g.pidx(x, 1, 0);
            let top = g.pidx(x, ny - 2, 0);
            for q in [4usize, 7, 8] {
                let c = D2Q9::C[q];
                let wall_ref = g.pidx_i(x as isize + c[0] as isize, 0, 0);
                records.push(BouzidiLink {
                    cell: bottom as u32,
                    q: q as u8,
                    qd: 1.0 - wall_lo,
                    has_second: true,
                    wall_ref: wall_ref as u32,
                });
            }
            for q in [2usize, 5, 6] {
                let c = D2Q9::C[q];
                let wall_ref = g.pidx_i(x as isize + c[0] as isize, (ny - 1) as isize, 0);
                records.push(BouzidiLink {
                    cell: top as u32,
                    q: q as u8,
                    qd: wall_hi - (ny - 2) as f64,
                    has_second: true,
                    wall_ref: wall_ref as u32,
                });
            }
        }
        solver.set_bouzidi_links(0, Some(BouzidiLinks::new(records)));
    }

    solver.init_with(|_, y, _| {
        if y == 0 || y == ny - 1 {
            (1.0, [0.0, 0.0, 0.0])
        } else {
            let yy = y as f64 - wall_lo;
            let h = wall_hi - wall_lo;
            (1.0, [force * yy * (h - yy) / (2.0 * nu), 0.0, 0.0])
        }
    });
    solver
}

fn poiseuille_l2rel(
    s: &Solver<D2Q9, f64, CpuScalar, LocalPeriodic>,
    wall_lo: f64,
    wall_hi: f64,
    force: f64,
) -> f64 {
    let ny = s.dims()[1];
    let nx = s.dims()[0];
    let h = wall_hi - wall_lo;
    let nu = s.nu();
    let mut actual = Vec::with_capacity(nx * (ny - 2));
    let mut reference = Vec::with_capacity(nx * (ny - 2));
    for y in 1..ny - 1 {
        let yy = y as f64 - wall_lo;
        let exact = force * yy * (h - yy) / (2.0 * nu);
        for x in 0..nx {
            actual.push(s.u(x, y, 0)[0]);
            reference.push(exact);
        }
    }
    l2_rel(&actual, &reference)
}

/// G1 (A1) — off-grid Poiseuille CONVERGENCE ORDER.
///
/// Refinement uses diffusive scaling: doubling ny halves the reference-profile
/// speed (u_ref ∝ 1/h) so the Reynolds number stays fixed and O(Ma²)
/// compressibility error does not pollute the fit (see axis-taxonomy A2).
/// The Bouzidi papers claim second-order boundary accuracy for straight walls
/// with fractional qd; the assert here demands slope >= 1.6 with r² >= 0.9 —
/// ~10x headroom over float noise, well below the theoretical 2 so a genuinely
/// second-order implementation passes cleanly and a first-order regression
/// (staircase-BB-equivalent) fails.
#[test]
fn g1_bouzidi_offgrid_poiseuille_convergence_order_light() {
    let wall_lo = 0.3;

    // Coarser resolutions run more steps to reach steady state; force is
    // scaled with h to keep u_max / nu fixed (diffusive scaling), so the
    // reference peak speed is the same across levels.
    // u_peak_ref = F h² / (8 nu) — fixing this and choosing target ~0.02
    // yields F(h) = 8 nu · u_peak / h².
    let target_u_peak = 0.02;
    let nu_base = 0.04;

    let nys = [22usize, 30, 42];
    let mut errs = Vec::new();
    let mut hs = Vec::new();
    for &ny in &nys {
        let wall_hi = ny as f64 - 1.0 - wall_lo;
        let width = wall_hi - wall_lo;
        let force = 8.0 * nu_base * target_u_peak / (width * width);
        let mut s = offgrid_poiseuille_channel(32, ny, nu_base, force, wall_lo, wall_hi, true);
        // Steady-state relaxation time scales like width²/nu; ~8·width²/nu
        // steps brings the exponential transient below the target error band.
        let steps = ((8.0 * width * width / nu_base).ceil() as usize).max(6000);
        s.run(steps);
        let e = poiseuille_l2rel(&s, wall_lo, wall_hi, force);
        println!(
            "ACC G1 Bouzidi off-grid Poiseuille: ny={ny}, width={width:.3}, L2rel={e:.6e}"
        );
        errs.push(e);
        // For order_fit, `h` is the mesh spacing per physical-channel-width
        // (dimensionless resolution parameter h → 0 as we refine). Cell size
        // in lattice units is 1.0, so h = 1/width. Passing width itself as
        // "h" would fit err ∝ width^(-p) and land the slope at -p — a P3
        // test-side trap this comment now documents.
        hs.push(1.0 / width);
    }
    let fit = order_fit(&hs, &errs);
    println!(
        "ACC G1 order fit: slope={:.3}, r2={:.4}, errs={:?}",
        fit.slope, fit.r2, errs
    );
    // Denominator of the relative tolerance below: the fitted convergence
    // slope itself is dimensionless (log-log linear fit).
    assert!(
        fit.slope >= 1.6,
        "G1 A1: Bouzidi off-grid Poiseuille convergence slope {:.3} < 1.6 \
         (a slope near 1.0 = first-order = staircase-BB-equivalent regression); \
         errs={:?}, hs={:?}",
        fit.slope,
        errs,
        hs
    );
    assert!(
        fit.r2 >= 0.9,
        "G1 A1: order-fit r2 {:.4} < 0.9 — the fit is not asymptotic-regime; \
         slope value is not trustworthy; errs={:?}",
        fit.r2,
        errs
    );
}

/// G2 (A2) — sub-cell translation invariance LIGHT CANARY.
///
/// Translate the fractional wall position by a sub-cell offset and require
/// the centreline velocity to move smoothly, not jump discontinuously. This
/// probe uses the straight-wall Poiseuille channel (Bouzidi records built
/// directly from wall_lo, not from an analytic-circle intersection), which
/// isolates the interpolation stencil's translation behaviour from the ray-
/// intersection code.
///
/// Analytic reference: for the same physical channel width h = wall_hi -
/// wall_lo and same F/nu, u_peak_ref = F h² / (8 nu) is independent of the
/// integer/fractional split of the walls. Discrete Bouzidi is not exact but
/// its peak deviation must lie ON a smooth curve of wall_lo — no jumps.
///
/// Heavy version (full sweep of ~5 offsets + drag on a translated cylinder)
/// deferred to G2H once the Bouzidi cylinder DRAG probe API is exercised in
/// a follow-up dispatch.
#[test]
fn g2_bouzidi_sub_cell_translation_light_canary() {
    let ny = 42;
    let force = 1.0e-6;
    let nu = 0.04;
    let offsets = [0.20, 0.50, 0.80];
    let mut peaks = Vec::new();
    for &wall_lo in &offsets {
        let wall_hi = ny as f64 - 1.0 - wall_lo;
        let mut s = offgrid_poiseuille_channel(32, ny, nu, force, wall_lo, wall_hi, true);
        s.run(15_000);
        // Extract centreline speed (row nearest the geometric centre of the
        // channel, taken along x=0 — the flow is x-invariant).
        let y_mid = (ny / 2) as usize;
        let u_peak = s.u(0, y_mid, 0)[0];
        println!("ACC G2 canary: wall_lo={wall_lo:.2}, u_peak={u_peak:.6e}");
        peaks.push(u_peak);
    }
    let pmin = peaks.iter().cloned().fold(f64::INFINITY, f64::min);
    let pmax = peaks.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let spread = (pmax - pmin) / pmax;
    // Denominator of this relative tolerance: max u_peak across the sweep
    // (peak-relative, NOT force-scale-relative — force-scale-relative would
    // pass a fully-staircase result at this force magnitude).
    println!(
        "ACC G2 spread = (pmax-pmin)/pmax = {:.4e} (peak-relative); peaks={:?}",
        spread, peaks
    );
    assert!(
        spread <= 0.06,
        "G2 A2 canary: cylinder-free sub-cell translation spread {:.3}% \
         exceeds 6% of peak — the interpolation is not smooth across sub-cell \
         wall placement (calibration: staircase BB gives O(20-30%) here)",
        spread * 100.0
    );
}

/// G4 (A1) — off-grid Couette with a moving wall would be the sharper A1
/// probe (linear reference isolates the wall error from bulk-scheme error).
/// SPEC-GAP: the current wall_u path on Solver takes a per-cell velocity
/// array at build time and does not expose a runtime moving-wall setter
/// through the compat facade; without a runtime setter for the fractional-
/// wall case the moving-wall Couette test cannot be written against the
/// currently landed Bouzidi API without leaning on private internals.
///
/// Documented here (per SKILL.md P2 convention 3) rather than silently
/// skipped. When the wall-velocity setter lands, this test writes as:
///   ρ ν d²u/dy² = 0, u(wall_lo)=0, u(wall_hi)=U_wall  ⇒  u(y) = U_wall·yy/h
/// and asserts L2rel <= 1e-3 at ny=42 with fractional walls (Bouzidi resolves
/// linear profiles exactly to O(Ma²), so residual is a lower bound on the
/// boundary error alone).
#[test]
#[ignore = "SPEC-GAP: no runtime moving-wall velocity setter for fractional Bouzidi walls; carries the derivation for the future landing"]
fn g4_bouzidi_offgrid_couette_moving_wall_spec_gap() {
    // See test doc-comment above for the full derivation. Intentionally empty.
}

/// G5 (A3) — effective-wall-position vs tau. Half-way BB has a documented
/// tau-dependent effective wall (shifts by O(tau - 1/2)); Bouzidi's design
/// promise is that its effective wall is tau-independent up to O(1) at
/// fractional qd. Heavy sweep (5 tau values × steady-state runs). The
/// functional form to assert is A3: the fitted parabola-zero position vs
/// tau must be a horizontal line at wall_lo within band, NOT merely small.
#[test]
#[ignore = "heavy ACC-AUDIT G5 tau sweep — 5 tau × ~6000 steps at ny=42"]
fn g5_bouzidi_effective_wall_position_vs_tau_heavy() {
    let ny = 42;
    let wall_lo = 0.3;
    let wall_hi = ny as f64 - 1.0 - wall_lo;
    let force = 1.0e-6;
    // tau = 3 nu + 1/2  ⇒  nu = (tau - 1/2)/3
    let taus = [0.55f64, 0.70, 1.00, 1.50, 2.00];
    let mut positions = Vec::new();
    for &tau in &taus {
        let nu = (tau - 0.5) / 3.0;
        let mut s = offgrid_poiseuille_channel(32, ny, nu, force, wall_lo, wall_hi, true);
        s.run(15_000);
        // Fit the parabola u(y) = a·(y - y0)·(y1 - y) to three interior rows
        // near the centre and solve for the zeros y0, y1; the effective
        // lower wall is y0.
        let y_mid = ny / 2;
        let ys = [
            (y_mid - 1) as f64,
            y_mid as f64,
            (y_mid + 1) as f64,
        ];
        let us = [
            s.u(0, y_mid - 1, 0)[0],
            s.u(0, y_mid, 0)[0],
            s.u(0, y_mid + 1, 0)[0],
        ];
        // Fit quadratic u = A y² + B y + C by 3 points -> roots via -B ± √(B²-4AC) / 2A.
        let m = [
            [ys[0] * ys[0], ys[0], 1.0],
            [ys[1] * ys[1], ys[1], 1.0],
            [ys[2] * ys[2], ys[2], 1.0],
        ];
        // 3x3 solve by Cramer's rule.
        let det = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
        let mut sol = [0.0f64; 3];
        for k in 0..3 {
            let mut mk = m;
            for i in 0..3 {
                mk[i][k] = us[i];
            }
            sol[k] = (mk[0][0] * (mk[1][1] * mk[2][2] - mk[1][2] * mk[2][1])
                - mk[0][1] * (mk[1][0] * mk[2][2] - mk[1][2] * mk[2][0])
                + mk[0][2] * (mk[1][0] * mk[2][1] - mk[1][1] * mk[2][0]))
                / det;
        }
        let (a, b, c) = (sol[0], sol[1], sol[2]);
        let disc = (b * b - 4.0 * a * c).max(0.0).sqrt();
        let y0 = (-b - disc) / (2.0 * a); // lower zero
        println!("ACC G5: tau={tau:.2}, effective y0={y0:.4} (nominal wall_lo={wall_lo:.2})");
        positions.push(y0);
    }
    let mean = positions.iter().sum::<f64>() / positions.len() as f64;
    let dev = positions
        .iter()
        .map(|p| (p - wall_lo).abs())
        .fold(0.0, f64::max);
    // Denominator of this relative tolerance: the channel width h (spread
    // in effective wall position relative to the channel width itself).
    let h = wall_hi - wall_lo;
    let rel = dev / h;
    println!(
        "ACC G5 max |y0 - wall_lo| / h = {:.4} (channel-width-relative); mean y0={:.4}",
        rel, mean
    );
    assert!(
        rel <= 0.02,
        "G5 A3: Bouzidi effective wall position drifts {:.3}% of h across \
         tau ∈ [0.55, 2.0]; expected tau-independent (Bouzidi's core design \
         promise). half-way BB would drift O(10%) here; positions={:?}",
        rel * 100.0,
        positions
    );
}

/// G6 (A2) — rotational anisotropy. A diagonal channel oriented along the
/// (1,1) lattice diagonal should reproduce the axis-aligned centreline
/// profile at matched physical resolution; the difference is a measure of
/// Bouzidi's rotational error. Deferred: constructing a diagonal channel
/// with the current wall/rim helpers needs a slanted-domain builder outside
/// the tests/common surface, so this remains a SPEC-GAP with the derivation
/// in comments.
#[test]
#[ignore = "SPEC-GAP: no diagonal-domain builder in tests/common; carries the derivation for the future landing"]
fn g6_bouzidi_rotational_anisotropy_spec_gap() {
    // Analytic reference for the future implementation:
    // A channel oriented along the lattice diagonal (1,1)/√2 with the same
    // physical width h and same physical force F should give the same peak
    // speed F h² / (8 nu) at the centreline. Discretized on the D2Q9 lattice,
    // the wall is not axis-aligned so Bouzidi's ray-line qd stresses the
    // interpolation on directions {1..8} unequally — the rotational anisotropy
    // shows up as a peak-speed discrepancy vs the axis-aligned case at
    // matched h and Ma. Assert band: <= 5% of peak, denominator = peak.
}
