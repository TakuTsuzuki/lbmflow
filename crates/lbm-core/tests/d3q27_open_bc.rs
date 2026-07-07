//! REV-5: D3Q27 open-face velocity/pressure closure.
//!
//! The D3Q27 kernel extends the existing Zou-He / Hecht-Harting
//! non-equilibrium bounce-back face closure from five to nine incoming links.
//! These gates pin the prescribed moments, duct behavior, D3Q19 consistency,
//! outflow/convective behavior, and split invariance when seams cross BC cells.

use lbm_core::prelude::*;
use std::f64::consts::PI;

type S19<H = LocalPeriodic> = Solver<D3Q19, f64, CpuScalar, H>;
type S27<H = LocalPeriodic> = Solver<D3Q27, f64, CpuScalar, H>;

fn opposite(face: Face) -> Face {
    face.opposite()
}

fn face_velocity(face: Face) -> [f64; 3] {
    let n = face.n_in();
    let (t1, t2) = face.tangents();
    let mut u = [0.0; 3];
    u[face.axis()] = 0.04 * n[face.axis()] as f64;
    u[t1] = 0.011;
    u[t2] = -0.013;
    u
}

fn face_cell_positions(face: Face, dims: [usize; 3]) -> Vec<[usize; 3]> {
    let a = face.axis();
    let fixed = if face.is_neg() { 0 } else { dims[a] - 1 };
    let (t1, t2) = face.tangents();
    let mut out = Vec::new();
    for c2 in 0..dims[t2] {
        for c1 in 0..dims[t1] {
            let mut p = [0usize; 3];
            p[a] = fixed;
            p[t1] = c1;
            p[t2] = c2;
            out.push(p);
        }
    }
    out
}

#[test]
fn d3q27_open_faces_enforce_velocity_and_pressure_moments_all_orientations() {
    for face in Face::ALL {
        let dims = [8, 7, 6];
        let out = opposite(face);
        let u_bc = face_velocity(face);
        let mut faces = [FaceBC::Closed; 6];
        faces[face.index()] = FaceBC::Velocity { u: u_bc };
        faces[out.index()] = FaceBC::Pressure { rho: 1.0 };
        let mut periodic = [true; 3];
        periodic[face.axis()] = false;
        let spec = GlobalSpec::<f64> {
            dims,
            nu: 0.05,
            periodic,
            faces,
            ..Default::default()
        };
        let mut s: S27 = Solver::new(
            &spec,
            &[],
            &[],
            [1, 1, 1],
            CpuScalar::default(),
            LocalPeriodic,
        );
        s.run(3);

        let mut du = 0.0f64;
        for p in face_cell_positions(face, dims) {
            let got = s.u(p[0], p[1], p[2]);
            for a in 0..3 {
                du = du.max((got[a] - u_bc[a]).abs());
            }
        }

        let mut drho = 0.0f64;
        let mut dut = 0.0f64;
        let (t1, t2) = out.tangents();
        for p in face_cell_positions(out, dims) {
            drho = drho.max((s.rho(p[0], p[1], p[2]) - 1.0).abs());
            let got = s.u(p[0], p[1], p[2]);
            dut = dut.max(got[t1].abs()).max(got[t2].abs());
        }
        println!(
            "D3Q27 {face:?}: max |u-u_bc|={du:.3e}, pressure max |rho-1|={drho:.3e}, max outlet transverse |u|={dut:.3e}"
        );
        assert!(
            du <= 2.0e-14,
            "{face:?} velocity inlet moment error {du:.3e}"
        );
        assert!(
            drho <= 2.0e-14,
            "{face:?} pressure outlet density error {drho:.3e}"
        );
        assert!(
            dut <= 2.0e-14,
            "{face:?} pressure outlet transverse velocity error {dut:.3e}"
        );
    }
}

fn duct_shape(y: usize, z: usize, ny: usize, nz: usize) -> f64 {
    if y == 0 || z == 0 || y + 1 == ny || z + 1 == nz {
        return 0.0;
    }
    let h = (ny - 2) as f64;
    let w = (nz - 2) as f64;
    let (a, b) = (h / 2.0, w / 2.0);
    let yy = y as f64 - 0.5;
    let zt = z as f64 - 0.5 - b;
    let pref = 16.0 * a * a / (PI * PI * PI);
    let mut sum = 0.0;
    let mut n = 1;
    while n <= 99 {
        let nf = n as f64;
        let kn = nf * PI / (2.0 * a);
        let ratio =
            ((kn * zt.abs()).exp() + (-kn * zt.abs()).exp()) / ((kn * b).exp() + (-kn * b).exp());
        sum += (1.0 - ratio) * (kn * yy).sin() / (nf * nf * nf);
        n += 2;
    }
    pref * sum
}

fn duct_profile(ny: usize, nz: usize, u_peak: f64) -> Vec<[f64; 3]> {
    let center = duct_shape(ny / 2, nz / 2, ny, nz);
    (0..nz)
        .flat_map(|z| {
            (0..ny).map(move |y| {
                let ux = u_peak * duct_shape(y, z, ny, nz) / center;
                [ux, 0.0, 0.0]
            })
        })
        .collect()
}

fn duct_spec_with_outlet(
    nx: usize,
    ny: usize,
    nz: usize,
    outlet: FaceBC<f64>,
) -> (GlobalSpec<f64>, Vec<bool>, Vec<[f64; 3]>) {
    let mut walls = WallSpec::<f64>::default();
    for face in [Face::YNeg, Face::YPos, Face::ZNeg, Face::ZPos] {
        walls.is_wall[face.index()] = true;
    }
    let mut faces = [FaceBC::Closed; 6];
    faces[Face::XNeg.index()] = FaceBC::Velocity { u: [0.0; 3] };
    faces[Face::XPos.index()] = outlet;
    let spec = GlobalSpec::<f64> {
        dims: [nx, ny, nz],
        nu: 0.05,
        periodic: [false, false, false],
        faces,
        ..Default::default()
    };
    let (solid, wall_u) = build_wall_rims(3, spec.dims, &walls);
    (spec, solid, wall_u)
}

fn duct_spec(nx: usize, ny: usize, nz: usize) -> (GlobalSpec<f64>, Vec<bool>, Vec<[f64; 3]>) {
    duct_spec_with_outlet(nx, ny, nz, FaceBC::Pressure { rho: 1.0 })
}

fn init_duct<L: Lattice, H: HaloExchange<f64>>(
    s: &mut Solver<L, f64, CpuScalar, H>,
    ny: usize,
    _nz: usize,
    profile: &[[f64; 3]],
) {
    s.set_inlet_profile(Face::XNeg, profile);
    s.init_with(|_, y, z| (1.0, profile[z * ny + y]));
}

fn run_to_steady<L: Lattice>(s: &mut Solver<L, f64, CpuScalar, LocalPeriodic>) -> (bool, f64) {
    let mut prev: Option<Vec<f64>> = None;
    let mut rel_last = f64::INFINITY;
    for _ in 0..40 {
        s.run(250);
        let mut cur = s.gather_ux();
        cur.extend(s.gather_uy());
        cur.extend(s.gather_uz());
        if let Some(p) = &prev {
            let dmax = cur
                .iter()
                .zip(p)
                .map(|(a, b)| (a - b).abs())
                .fold(0.0f64, f64::max);
            let umax = cur.iter().map(|v| v.abs()).fold(0.0f64, f64::max);
            if umax > 0.0 && dmax <= 1.0e-10 * umax {
                return (true, dmax / umax);
            }
            if umax > 0.0 {
                rel_last = dmax / umax;
            }
        }
        prev = Some(cur);
    }
    (false, rel_last)
}

struct DuctMetrics {
    l2_profile: f64,
    l2_profile_flux_scaled: f64,
    profile_flux_scale: f64,
    l2_vs_d3q19: f64,
    bulk_flux_rel: f64,
    flux_rel: f64,
    cross_rel: f64,
    pressure_drop_ok: bool,
}

fn analyze_duct(d27: &S27, d19: &S19, profile: &[[f64; 3]], artifact: Option<&str>) -> DuctMetrics {
    let [nx, ny, nz] = d27.dims();
    let x = nx / 2;
    let ux27 = d27.gather_ux();
    let uy27 = d27.gather_uy();
    let uz27 = d27.gather_uz();
    let rho27 = d27.gather_rho();
    let ux19 = d19.gather_ux();
    let idx = |x: usize, y: usize, z: usize| (z * ny + y) * nx + x;

    let mut e2 = 0.0;
    let mut r2 = 0.0;
    let mut dot_ref = 0.0;
    let mut e19 = 0.0;
    let mut r19 = 0.0;
    let mut cross = 0.0f64;
    let mut umax = 0.0f64;
    let mut csv = String::from("y,z,u_d3q27,u_d3q19,u_ref\n");
    for z in 1..nz - 1 {
        for y in 1..ny - 1 {
            let u = ux27[idx(x, y, z)];
            let r = profile[z * ny + y][0];
            let u19 = ux19[idx(x, y, z)];
            e2 += (u - r) * (u - r);
            r2 += r * r;
            dot_ref += u * r;
            e19 += (u - u19) * (u - u19);
            r19 += u19 * u19;
            cross = cross
                .max(uy27[idx(x, y, z)].abs())
                .max(uz27[idx(x, y, z)].abs());
            umax = umax.max(u.abs());
            if artifact.is_some() {
                csv.push_str(&format!("{y},{z},{u:.16e},{u19:.16e},{r:.16e}\n"));
            }
        }
    }
    if let Some(label) = artifact {
        std::fs::create_dir_all("target").expect("create target dir");
        std::fs::write(format!("target/d3q27_{label}_duct_profile.csv"), csv)
            .expect("write duct profile CSV");
    }

    let flux = |xx: usize| -> f64 {
        let mut q = 0.0;
        for z in 1..nz - 1 {
            for y in 1..ny - 1 {
                q += rho27[idx(xx, y, z)] * ux27[idx(xx, y, z)];
            }
        }
        q
    };
    let q_in = flux(1);
    let q_mid = flux(nx / 2);
    let q_out = flux(nx - 2);
    let bulk_flux_rel = (q_in - q_mid).abs() / q_mid.abs();
    let flux_rel = ((q_in - q_out).abs().max((q_mid - q_out).abs())) / q_mid.abs();
    let plane_rho = |xx: usize| -> f64 {
        let mut r = 0.0;
        let mut n = 0.0;
        for z in 1..nz - 1 {
            for y in 1..ny - 1 {
                r += rho27[idx(xx, y, z)];
                n += 1.0;
            }
        }
        r / n
    };
    let rho_in = plane_rho(1);
    let rho_mid = plane_rho(nx / 2);
    let rho_out = plane_rho(nx - 2);
    let scale = dot_ref / r2;
    let mut e2_scaled = 0.0;
    for z in 1..nz - 1 {
        for y in 1..ny - 1 {
            let u = ux27[idx(x, y, z)];
            let r = scale * profile[z * ny + y][0];
            e2_scaled += (u - r) * (u - r);
        }
    }
    DuctMetrics {
        l2_profile: (e2 / r2).sqrt(),
        l2_profile_flux_scaled: (e2_scaled / r2).sqrt(),
        profile_flux_scale: scale,
        l2_vs_d3q19: (e19 / r19).sqrt(),
        bulk_flux_rel,
        flux_rel,
        cross_rel: cross / umax,
        pressure_drop_ok: rho_in > rho_mid && rho_mid > rho_out,
    }
}

#[test]
fn d3q27_open_duct_matches_series_shape_and_d3q19() {
    let (nx, ny, nz) = (40, 16, 16);
    let profile = duct_profile(ny, nz, 0.035);
    let (spec, solid, wall_u) = duct_spec(nx, ny, nz);
    let mut d27: S27 = Solver::new(
        &spec,
        &solid,
        &wall_u,
        [1, 1, 1],
        CpuScalar::default(),
        LocalPeriodic,
    );
    let mut d19: S19 = Solver::new(
        &spec,
        &solid,
        &wall_u,
        [1, 1, 1],
        CpuScalar::default(),
        LocalPeriodic,
    );
    init_duct(&mut d27, ny, nz, &profile);
    init_duct(&mut d19, ny, nz, &profile);
    let (steady27, rel27) = run_to_steady(&mut d27);
    let (steady19, rel19) = run_to_steady(&mut d19);

    let m = analyze_duct(&d27, &d19, &profile, Some("pressure"));
    println!(
        "D3Q27 open duct: steady27={steady27} rel27={rel27:.3e}, steady19={steady19} rel19={rel19:.3e}, L2 profile={:.3e}, scaled L2={:.3e} (scale {:.6}), D3Q27-vs-D3Q19={:.3e}, flux_rel={:.3e}, cross_rel={:.3e}, pressure_drop_ok={}",
        m.l2_profile, m.l2_profile_flux_scaled, m.profile_flux_scale, m.l2_vs_d3q19, m.flux_rel, m.cross_rel, m.pressure_drop_ok
    );
    assert!(
        rel27 <= 1.0e-7,
        "D3Q27 open duct did not settle enough: rel={rel27:.3e}"
    );
    assert!(
        rel19 <= 1.0e-7,
        "D3Q19 open duct did not settle enough: rel={rel19:.3e}"
    );
    assert!(
        m.l2_profile <= 1.0e-2,
        "D3Q27 duct unscaled profile L2rel={:.3e}",
        m.l2_profile
    );
    assert!(
        m.l2_profile_flux_scaled <= 1.0e-3,
        "D3Q27 duct flux-scaled profile L2rel={:.3e}",
        m.l2_profile_flux_scaled
    );
    assert!(
        m.l2_vs_d3q19 <= 2.0e-3,
        "D3Q27 duct differs from D3Q19 by L2rel={:.3e}",
        m.l2_vs_d3q19
    );
    assert!(
        m.flux_rel <= 2.0e-4,
        "D3Q27 duct mass flux imbalance={:.3e}",
        m.flux_rel
    );
    assert!(
        m.cross_rel <= 1.0e-3,
        "D3Q27 duct cross-flow ratio={:.3e}",
        m.cross_rel
    );
    assert!(
        m.pressure_drop_ok,
        "D3Q27 duct pressure is not monotone inlet to outlet"
    );
}

fn run_d3q27_open_outlet_duct(outlet: FaceBC<f64>, label: &str) {
    let (nx, ny, nz) = (40, 16, 16);
    let profile = duct_profile(ny, nz, 0.01);
    let (spec, solid, wall_u) = duct_spec_with_outlet(nx, ny, nz, outlet);
    let mut d27: S27 = Solver::new(
        &spec,
        &solid,
        &wall_u,
        [1, 1, 1],
        CpuScalar::default(),
        LocalPeriodic,
    );
    let mut d19: S19 = Solver::new(
        &spec,
        &solid,
        &wall_u,
        [1, 1, 1],
        CpuScalar::default(),
        LocalPeriodic,
    );
    init_duct(&mut d27, ny, nz, &profile);
    init_duct(&mut d19, ny, nz, &profile);
    let (steady27, rel27) = run_to_steady(&mut d27);
    let (steady19, rel19) = run_to_steady(&mut d19);
    let m = analyze_duct(&d27, &d19, &profile, Some(label));
    println!(
        "D3Q27 {label} duct: steady27={steady27} rel27={rel27:.3e}, steady19={steady19} rel19={rel19:.3e}, L2 profile={:.3e}, scaled L2={:.3e} (scale {:.6}), D3Q27-vs-D3Q19={:.3e}, bulk_flux_rel={:.3e}, outlet_flux_rel={:.3e}, cross_rel={:.3e}, pressure_drop_ok={}",
        m.l2_profile,
        m.l2_profile_flux_scaled,
        m.profile_flux_scale,
        m.l2_vs_d3q19,
        m.bulk_flux_rel,
        m.flux_rel,
        m.cross_rel,
        m.pressure_drop_ok
    );
    assert!(
        rel27 <= 2.0e-7,
        "{label}: D3Q27 did not settle: rel={rel27:.3e}"
    );
    assert!(
        rel19 <= 2.0e-7,
        "{label}: D3Q19 did not settle: rel={rel19:.3e}"
    );
    assert!(
        m.l2_profile_flux_scaled
            <= if label == "convective" {
                3.0e-3
            } else {
                1.5e-3
            },
        "{label}: D3Q27 flux-scaled duct profile L2rel={:.3e}",
        m.l2_profile_flux_scaled
    );
    assert!(
        m.l2_vs_d3q19 <= 2.5e-3,
        "{label}: D3Q27 differs from D3Q19 by L2rel={:.3e}",
        m.l2_vs_d3q19
    );
    assert!(
        m.flux_rel <= if label == "convective" { 2.1 } else { 2.5e-1 },
        "{label}: D3Q27 outlet-local flux envelope={:.3e}",
        m.flux_rel
    );
    assert!(
        m.cross_rel
            <= if label == "convective" {
                4.0e-2
            } else {
                5.0e-3
            },
        "{label}: D3Q27 cross-flow ratio={:.3e}",
        m.cross_rel
    );
    assert!(
        m.pressure_drop_ok,
        "{label}: pressure is not monotone inlet to outlet"
    );
}

#[test]
fn d3q27_outflow_duct_matches_profile_and_d3q19() {
    run_d3q27_open_outlet_duct(FaceBC::Outflow, "outflow");
}

#[test]
fn d3q27_convective_duct_matches_profile_and_d3q19() {
    run_d3q27_open_outlet_duct(FaceBC::Convective { u_conv: 0.01 }, "convective");
}

#[test]
fn d3q27_open_outlets_balance_uniform_through_flow_mass_flux() {
    let dims = [24, 8, 8];
    for (outlet, label) in [
        (FaceBC::Outflow, "outflow"),
        (FaceBC::Convective { u_conv: 0.01 }, "convective"),
    ] {
        let mut faces = [FaceBC::Closed; 6];
        faces[Face::XNeg.index()] = FaceBC::Velocity {
            u: [0.01, 0.0, 0.0],
        };
        faces[Face::XPos.index()] = outlet;
        let spec = GlobalSpec::<f64> {
            dims,
            nu: 0.05,
            periodic: [false, true, true],
            faces,
            ..Default::default()
        };
        let mut s: S27 = Solver::new(
            &spec,
            &[],
            &[],
            [1, 1, 1],
            CpuScalar::default(),
            LocalPeriodic,
        );
        s.init_with(|_, _, _| (1.0, [0.01, 0.0, 0.0]));
        s.run(200);
        let [nx, ny, nz] = dims;
        let ux = s.gather_ux();
        let rho = s.gather_rho();
        let idx = |x: usize, y: usize, z: usize| (z * ny + y) * nx + x;
        let flux = |x: usize| -> f64 {
            let mut q = 0.0;
            for z in 0..nz {
                for y in 0..ny {
                    q += rho[idx(x, y, z)] * ux[idx(x, y, z)];
                }
            }
            q
        };
        let q1 = flux(1);
        let qm = flux(nx / 2);
        let qo = flux(nx - 2);
        let rel = (q1 - qm).abs().max((qm - qo).abs()) / qm.abs();
        println!("D3Q27 {label} uniform duct flux rel={rel:.3e}");
        assert!(
            rel <= 1.0e-10,
            "{label}: uniform mass-flux imbalance {rel:e}"
        );
    }
}

fn assert_d3q27_split_equal<HA: HaloExchange<f64>, HB: HaloExchange<f64>>(
    a: &S27<HA>,
    b: &S27<HB>,
    what: &str,
) {
    for (name, va, vb) in [
        ("rho", a.gather_rho(), b.gather_rho()),
        ("ux", a.gather_ux(), b.gather_ux()),
        ("uy", a.gather_uy(), b.gather_uy()),
        ("uz", a.gather_uz(), b.gather_uz()),
    ] {
        let d = va
            .iter()
            .zip(&vb)
            .map(|(x, y)| (x - y).abs())
            .fold(0.0f64, f64::max);
        assert_eq!(d, 0.0, "{what}: {name} differs by {d:e}");
    }
    for q in 0..D3Q27::Q {
        let (fa, fb) = (a.gather_f(q), b.gather_f(q));
        let d = fa
            .iter()
            .zip(&fb)
            .map(|(x, y)| (x - y).abs())
            .fold(0.0f64, f64::max);
        assert_eq!(d, 0.0, "{what}: f[{q}] differs by {d:e}");
    }
}

#[test]
fn t13_d3q27_open_duct_split_invariant_with_bc_seams() {
    let (nx, ny, nz) = (24, 12, 10);
    let profile = duct_profile(ny, nz, 0.03);

    for (outlet, label) in [
        (FaceBC::Pressure { rho: 1.0 }, "pressure"),
        (FaceBC::Outflow, "outflow"),
        (FaceBC::Convective { u_conv: 0.03 }, "convective"),
    ] {
        let (spec, solid, wall_u) = duct_spec_with_outlet(nx, ny, nz, outlet);
        for decomp in [[1, 2, 1], [1, 1, 2], [2, 2, 1]] {
            let mut base: S27 = Solver::new(
                &spec,
                &solid,
                &wall_u,
                [1, 1, 1],
                CpuScalar::default(),
                LocalPeriodic,
            );
            init_duct(&mut base, ny, nz, &profile);
            let mut split: S27<InProcess> = Solver::new(
                &spec,
                &solid,
                &wall_u,
                decomp,
                CpuScalar::default(),
                InProcess,
            );
            split.set_two_pass(true);
            init_duct(&mut split, ny, nz, &profile);
            for t in 1..=80 {
                base.step();
                split.step();
                if t <= 3 || t % 20 == 0 {
                    assert_d3q27_split_equal(
                        &base,
                        &split,
                        &format!("{label} decomp {decomp:?} t={t}"),
                    );
                }
            }
        }
    }
}
