//! 3D baffled stirred tank with a Rushton-turbine impeller — velocity and
//! shear-stress fields, plus a subsampled 3D volume for an interactive viewer.
//!
//! Standard stirred-tank proportions (tank diameter T, liquid height H ~ T,
//! impeller diameter D ~ T/3, off-bottom clearance C ~ T/3, four wall baffles
//! of width ~ T/10). The impeller is a six-blade Rushton turbine (hub + disk +
//! six flat blades) on a central shaft.
//!
//! The impeller is NOT a resolved moving solid (the core has no moving-boundary
//! / IBM API yet — docs/REQ_STIRRED_REACTOR.md FR-ROT-01, pending M-F). It is
//! modelled by *volume penalization*: fluid cells inside the rotating turbine
//! footprint get a Guo body force (via the public Solver::set_body_force_field
//! API) that drags their velocity toward the local rigid-body velocity
//! v = omega x r. Baffles are true no-slip solids (half-way bounce-back).
//!
//! Outputs into <outdir>:
//!   - top_speed_*.png / top_shear_*.png : impeller-plane slices (animation)
//!   - side_speed_*.png                  : vertical mid-plane slice
//!   - volume.bin + volume.json          : subsampled (vx,vy,vz,shear) volume
//!                                         + geometry meta for the 3D viewer
//!
//! Run:
//!   cargo run --release --example stirred_tank_3d -- <outdir> [n] [steps] [every]

use lbm_core::prelude::*;
use std::f64::consts::PI;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

type S3 = Solver<D3Q19, f64, CpuScalar, LocalPeriodic>;

const VIRIDIS: [[u8; 3]; 16] = [
    [68, 1, 84], [72, 26, 108], [71, 47, 125], [65, 68, 135],
    [57, 86, 140], [49, 104, 142], [42, 120, 142], [35, 136, 142],
    [31, 152, 139], [34, 168, 132], [53, 183, 121], [84, 197, 104],
    [122, 209, 81], [165, 219, 54], [210, 226, 27], [253, 231, 37],
];
const INFERNO: [[u8; 3]; 16] = [
    [0, 0, 4], [12, 8, 38], [36, 12, 79], [66, 10, 104],
    [93, 18, 110], [120, 28, 109], [147, 38, 103], [174, 48, 92],
    [199, 62, 76], [220, 81, 57], [237, 105, 37], [246, 133, 17],
    [251, 163, 12], [249, 195, 41], [240, 226, 96], [252, 255, 164],
];

fn lut(anchors: &[[u8; 3]], t: f64) -> [u8; 3] {
    let t = t.clamp(0.0, 1.0) * (anchors.len() - 1) as f64;
    let i = (t as usize).min(anchors.len() - 2);
    let f = t - i as f64;
    let (a, b) = (anchors[i], anchors[i + 1]);
    [
        (a[0] as f64 + (b[0] as f64 - a[0] as f64) * f) as u8,
        (a[1] as f64 + (b[1] as f64 - a[1] as f64) * f) as u8,
        (a[2] as f64 + (b[2] as f64 - a[2] as f64) * f) as u8,
    ]
}

fn write_png(path: &Path, field: &[f64], mask: &[bool], w: usize, h: usize,
             vmax: f64, anchors: &[[u8; 3]], sc: usize) {
    let (ow, oh) = (w * sc, h * sc);
    let mut buf = vec![0u8; ow * oh * 3];
    for oy in 0..oh {
        let y = oy / sc;
        for ox in 0..ow {
            let x = ox / sc;
            let i = y * w + x;
            let rgb = if mask[i] { [92u8, 96, 104] } else { lut(anchors, field[i] / vmax.max(1e-30)) };
            let px = ((oh - 1 - oy) * ow + ox) * 3;
            buf[px..px + 3].copy_from_slice(&rgb);
        }
    }
    let file = File::create(path).expect("create png");
    let mut enc = png::Encoder::new(BufWriter::new(file), ow as u32, oh as u32);
    enc.set_color(png::ColorType::Rgb);
    enc.set_depth(png::BitDepth::Eight);
    let mut wr = enc.write_header().expect("png header");
    wr.write_image_data(&buf).expect("png data");
}

/// Rushton-turbine geometry (all lengths in lattice cells), derived from n.
struct Geom {
    n: usize,
    cx: f64,
    cy: f64,
    r_tank: f64,
    zc: f64,       // impeller mid-plane height
    tip_r: f64,    // blade tip radius (= D/2)
    disk_r: f64,
    hub_r: f64,
    shaft_r: f64,
    blade_hh: f64, // blade half-height in z
    disk_hh: f64,
    n_blades: usize,
    blade_hw: f64, // blade half-thickness (perp)
    baffle_len: f64,
    baffle_hw: f64,
}

impl Geom {
    fn new(n: usize) -> Self {
        let r_tank = n as f64 / 2.0 - 3.0;
        let tip_r = r_tank / 3.0; // D = T/3 -> tip radius T/6 = R/3
        Geom {
            n,
            cx: (n as f64 - 1.0) / 2.0,
            cy: (n as f64 - 1.0) / 2.0,
            r_tank,
            zc: n as f64 * 0.35, // clearance ~ T/3 from the floor
            tip_r,
            disk_r: tip_r * 0.66,
            hub_r: tip_r * 0.22,
            shaft_r: (tip_r * 0.12).max(1.5),
            blade_hh: tip_r * 0.30,
            disk_hh: 1.2,
            n_blades: 6,
            blade_hw: 1.2,
            baffle_len: r_tank * 0.2, // width ~ T/10
            baffle_hw: 1.5,
        }
    }
    fn rad(&self, x: usize, y: usize) -> f64 {
        let dx = x as f64 - self.cx;
        let dy = y as f64 - self.cy;
        (dx * dx + dy * dy).sqrt()
    }
    /// True where a static solid sits: outer rim, outside the tank, or a baffle.
    fn is_solid(&self, x: usize, y: usize, z: usize) -> bool {
        let n = self.n;
        if x == 0 || x == n - 1 || y == 0 || y == n - 1 || z == 0 || z == n - 1 {
            return true;
        }
        let r = self.rad(x, y);
        if r > self.r_tank {
            return true;
        }
        // Four full-height wall baffles at 0/90/180/270 degrees.
        let dx = x as f64 - self.cx;
        let dy = y as f64 - self.cy;
        for k in 0..4 {
            let beta = k as f64 * PI / 2.0;
            let s = dx * beta.cos() + dy * beta.sin(); // radial coord along baffle
            let p = -dx * beta.sin() + dy * beta.cos(); // perpendicular
            if s > self.r_tank - self.baffle_len && s <= self.r_tank && p.abs() <= self.baffle_hw {
                return true;
            }
        }
        false
    }
    /// If (x,y,z) is inside the rotating turbine, return the tangential unit
    /// direction (rigid-body spin about +z). None otherwise.
    fn impeller_dir(&self, x: usize, y: usize, z: usize, theta: f64) -> Option<[f64; 2]> {
        let dx = x as f64 - self.cx;
        let dy = y as f64 - self.cy;
        let r = (dx * dx + dy * dy).sqrt();
        if r < 1e-6 {
            return None;
        }
        let t = [-dy / r, dx / r];
        let zf = z as f64;
        // Shaft: thin column from the impeller up to the top.
        if r <= self.shaft_r && zf >= self.zc {
            return Some(t);
        }
        if (zf - self.zc).abs() <= self.blade_hh {
            // Hub.
            if r <= self.hub_r {
                return Some(t);
            }
            // Disk (thin plate).
            if r <= self.disk_r && (zf - self.zc).abs() <= self.disk_hh {
                return Some(t);
            }
            // Six flat blades on the disk rim.
            let phi = dy.atan2(dx);
            if r >= self.disk_r * 0.85 && r <= self.tip_r {
                for b in 0..self.n_blades {
                    let beta = theta + b as f64 * 2.0 * PI / self.n_blades as f64;
                    let d = (r * (phi - beta).sin()).abs();
                    if d <= self.blade_hw && (phi - beta).cos() > 0.0 {
                        return Some(t);
                    }
                }
            }
        }
        None
    }
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let outdir = PathBuf::from(a.get(1).map(String::as_str).unwrap_or("stirred_tank_out"));
    let n: usize = a.get(2).and_then(|s| s.parse().ok()).unwrap_or(80);
    let steps: usize = a.get(3).and_then(|s| s.parse().ok()).unwrap_or(8000);
    let every: usize = a.get(4).and_then(|s| s.parse().ok()).unwrap_or(40);
    let u_tip: f64 = a.get(5).and_then(|s| s.parse().ok()).unwrap_or(0.08);
    let nu: f64 = a.get(6).and_then(|s| s.parse().ok()).unwrap_or(0.02);
    std::fs::create_dir_all(&outdir).expect("mkdir outdir");

    let g = Geom::new(n);
    let (nx, ny, nz) = (n, n, n);
    let idx = |x: usize, y: usize, z: usize| (z * ny + y) * nx + x;

    let omega = u_tip / g.tip_r;
    let alpha = 0.32_f64;
    let spin_up = 1500.0_f64;
    let f_cap = 0.25 * u_tip; // scale the penalization cap with tip speed (0.08 -> 0.02)
    let re = u_tip * (2.0 * g.tip_r) / nu;
    let cs = 1.0_f64 / 3.0_f64.sqrt();
    let ma_tip = u_tip / cs;
    let tau = 3.0 * nu + 0.5;

    // ---- build solver ------------------------------------------------------
    let mut walls = WallSpec::<f64>::default();
    for f in [Face::XNeg, Face::XPos, Face::YNeg, Face::YPos, Face::ZNeg, Face::ZPos] {
        walls.is_wall[f.index()] = true;
    }
    let spec = GlobalSpec::<f64> {
        dims: [nx, ny, nz],
        nu,
        periodic: [false, false, false],
        collision: CollisionKind::Trt { magic: CollisionKind::MAGIC_STD },
        ..Default::default()
    };
    let (solid, wall_u) = build_wall_rims(3, spec.dims, &walls);
    let mut sim = S3::new(&spec, &solid, &wall_u, [1, 1, 1], CpuScalar::default(), LocalPeriodic);

    // Carve the round tank wall + baffles (static solids).
    for z in 0..nz {
        for y in 0..ny {
            for x in 0..nx {
                let rim = x == 0 || x == nx - 1 || y == 0 || y == ny - 1;
                if !rim && g.is_solid(x, y, z) {
                    sim.set_solid(x, y, z);
                }
            }
        }
    }
    sim.mark_masks_dirty();

    let n3 = nx * ny * nz;
    let mut force = vec![[0.0f64; 3]; n3];

    // Impeller force bounding box (z-band tip box) + shaft column (thin, tall).
    let zb0 = (g.zc - g.blade_hh - 1.0).floor().max(1.0) as usize;
    let zb1 = (g.zc + g.blade_hh + 1.0).ceil().min(nz as f64 - 2.0) as usize;
    let bb = (g.tip_r + 2.0).ceil() as i64;
    let (bx0, bx1) = ((g.cx as i64 - bb).max(1) as usize, (g.cx as i64 + bb).min(nx as i64 - 2) as usize);
    let (by0, by1) = ((g.cy as i64 - bb).max(1) as usize, (g.cy as i64 + bb).min(ny as i64 - 2) as usize);
    let sb = (g.shaft_r + 2.0).ceil() as i64;
    let (sx0, sx1) = ((g.cx as i64 - sb).max(1) as usize, (g.cx as i64 + sb).min(nx as i64 - 2) as usize);
    let (sy0, sy1) = ((g.cy as i64 - sb).max(1) as usize, (g.cy as i64 + sb).min(ny as i64 - 2) as usize);

    println!(
        "Rushton stirred tank: {n}^3, tank r={:.1}, tip r={:.1}, {} blades, \
         omega={omega:.5}/step (period {:.0}), u_tip={u_tip}, nu={nu}, tau={tau:.3}, \
         Ma_tip={ma_tip:.3}, Re~{re:.0}, steps={steps}",
        g.r_tank, g.tip_r, g.n_blades, 2.0 * PI / omega
    );

    let v_speed = u_tip;
    let v_shear = 6.0e-4;
    let scale = if n <= 96 { 5 } else { 4 };
    let mut frame = 0usize;

    let stamp_region = |force: &mut [[f64; 3]], theta: f64, ramp: f64,
                        ux: &[f64], uy: &[f64], x: usize, y: usize, z: usize| {
        if g.is_solid(x, y, z) {
            return;
        }
        if let Some(t) = g.impeller_dir(x, y, z, theta) {
            let i = idx(x, y, z);
            let vb = omega * g.rad(x, y) * ramp;
            let mut fx = 2.0 * alpha * (t[0] * vb - ux[i]);
            let mut fy = 2.0 * alpha * (t[1] * vb - uy[i]);
            let fm = (fx * fx + fy * fy).sqrt();
            if fm > f_cap {
                fx *= f_cap / fm;
                fy *= f_cap / fm;
            }
            force[i] = [fx, fy, 0.0];
        }
    };

    // ---- time loop ---------------------------------------------------------
    for step in 0..steps {
        let theta = omega * step as f64;
        let ramp = (step as f64 / spin_up).min(1.0);
        let ux = sim.gather_ux();
        let uy = sim.gather_uy();

        // Zero the two force footprints, then stamp the turbine + shaft.
        for z in zb0..=zb1 {
            for y in by0..=by1 {
                for x in bx0..=bx1 {
                    force[idx(x, y, z)] = [0.0, 0.0, 0.0];
                }
            }
        }
        for z in (g.zc as usize)..nz - 1 {
            for y in sy0..=sy1 {
                for x in sx0..=sx1 {
                    force[idx(x, y, z)] = [0.0, 0.0, 0.0];
                }
            }
        }
        for z in zb0..=zb1 {
            for y in by0..=by1 {
                for x in bx0..=bx1 {
                    stamp_region(&mut force, theta, ramp, &ux, &uy, x, y, z);
                }
            }
        }
        for z in (g.zc as usize)..nz - 1 {
            for y in sy0..=sy1 {
                for x in sx0..=sx1 {
                    stamp_region(&mut force, theta, ramp, &ux, &uy, x, y, z);
                }
            }
        }

        sim.set_body_force_field(|x, y, z| force[idx(x, y, z)]);
        sim.step();

        if step % every == 0 {
            let ux = sim.gather_ux();
            let uy = sim.gather_uy();
            let uz = sim.gather_uz();
            let speed: Vec<f64> = (0..n3)
                .map(|i| (ux[i] * ux[i] + uy[i] * uy[i] + uz[i] * uz[i]).sqrt())
                .collect();
            let shear = shear_field(&ux, &uy, &uz, nx, ny, nz, nu);

            // top view (impeller plane)
            let zc = g.zc as usize;
            let mut ts = vec![0.0; nx * ny];
            let mut th = vec![0.0; nx * ny];
            let mut tm = vec![false; nx * ny];
            for y in 0..ny {
                for x in 0..nx {
                    let i2 = y * nx + x;
                    ts[i2] = speed[idx(x, y, zc)];
                    th[i2] = shear[idx(x, y, zc)];
                    tm[i2] = g.is_solid(x, y, zc) || g.impeller_dir(x, y, zc, theta).is_some();
                }
            }
            // side view (vertical mid-plane)
            let ymid = ny / 2;
            let mut ss = vec![0.0; nx * nz];
            let mut sm = vec![false; nx * nz];
            for z in 0..nz {
                for x in 0..nx {
                    let i2 = z * nx + x;
                    ss[i2] = speed[idx(x, ymid, z)];
                    sm[i2] = g.is_solid(x, ymid, z) || g.impeller_dir(x, ymid, z, theta).is_some();
                }
            }
            let name = |p: &str| outdir.join(format!("{p}_{frame:04}.png"));
            write_png(&name("top_speed"), &ts, &tm, nx, ny, v_speed, &VIRIDIS, scale);
            write_png(&name("top_shear"), &th, &tm, nx, ny, v_shear, &INFERNO, scale);
            write_png(&name("side_speed"), &ss, &sm, nx, nz, v_speed, &VIRIDIS, scale);

            let smax = speed.iter().cloned().fold(0.0f64, f64::max);
            if frame % 10 == 0 {
                println!("  step {step}/{steps} -> frame {frame}  max|u|={smax:.4}");
            }
            frame += 1;
            if !smax.is_finite() || smax > 0.5 {
                println!("  DIVERGED at step {step}, max|u|={smax:.4}");
                break;
            }
        }
    }

    // ---- export subsampled volume for the interactive viewer ---------------
    let ux = sim.gather_ux();
    let uy = sim.gather_uy();
    let uz = sim.gather_uz();
    // Shear stress rho*nu*gamma_dot. Prefer the core-native strain-rate
    // (f_neq / FR-STRESS-01, gather_shear_rate) over the finite-difference
    // reconstruction; report the delta between the two (native = reference).
    let shear_fd = shear_field(&ux, &uy, &uz, nx, ny, nz, nu);
    let shear: Vec<f64> = sim.gather_shear_rate().iter().map(|gd| nu * gd).collect();
    {
        let (mut num, mut den, mut maxd) = (0.0f64, 0.0f64, 0.0f64);
        for i in 0..n3 {
            let d = (shear_fd[i] - shear[i]).abs();
            num += d;
            den += shear[i].abs();
            maxd = maxd.max(d);
        }
        let fd_max = shear_fd.iter().cloned().fold(0.0f64, f64::max);
        let nat_max = shear.iter().cloned().fold(0.0f64, f64::max);
        println!(
            "SHEAR native(f_neq) vs FD: mean|Δ|/mean={:.1}%  max|Δ|={:.2e}  \
             native_max={nat_max:.5}  FD_max={fd_max:.5}",
            if den > 0.0 { 100.0 * num / den } else { 0.0 },
            maxd
        );
    }
    let final_max = (0..n3)
        .map(|i| (ux[i] * ux[i] + uy[i] * uy[i] + uz[i] * uz[i]).sqrt())
        .fold(0.0f64, f64::max);
    let shear_max = shear.iter().cloned().fold(0.0f64, f64::max);
    // Stability verdict against the documented tune-stability envelope, not just
    // "did it blow up": a bounded run can still be silently wrong (Ma>0.3 =
    // compressible; grid-Re U/nu>15 = under-resolved). See docs/qa/anomaly-log.md A1/A2.
    let ma_field = final_max / cs;
    let grid_re = final_max / nu;
    let status = if !final_max.is_finite() || final_max >= 0.5 {
        "DIVERGED"
    } else if ma_field > 0.3 || grid_re > 15.0 {
        "OUT-OF-ENVELOPE" // bounded but outside the validated regime
    } else {
        "STABLE"
    };
    println!(
        "SUMMARY u_tip={u_tip} nu={nu} tau={tau:.3} Ma_tip={ma_tip:.3} Re~{re:.0} \
         final_max|u|={final_max:.4} Ma_field={ma_field:.2} grid_Re={grid_re:.0} \
         shear_max={shear_max:.5} -> {status}"
    );
    export_volume(&outdir, &g, &ux, &uy, &uz, &shear, nx, ny, nz, omega, u_tip, nu);

    println!("\nWrote {frame} frames x3 + volume.bin/json to {}", outdir.display());
    for s in ["top_speed", "top_shear", "side_speed"] {
        println!("ffmpeg -y -framerate 25 -i {0}/{s}_%04d.png -c:v libx264 -pix_fmt yuv420p \
                  -crf 18 -vf scale=trunc(iw/2)*2:trunc(ih/2)*2 {0}/{s}.mp4", outdir.display());
    }
}

/// Shear stress rho*nu*gamma_dot, gamma_dot = sqrt(2 S:S), central differences.
fn shear_field(ux: &[f64], uy: &[f64], uz: &[f64], nx: usize, ny: usize, nz: usize, nu: f64) -> Vec<f64> {
    let idx = |x: usize, y: usize, z: usize| (z * ny + y) * nx + x;
    let mut out = vec![0.0; nx * ny * nz];
    for z in 1..nz - 1 {
        for y in 1..ny - 1 {
            for x in 1..nx - 1 {
                let d = |f: &[f64], ax: usize| -> f64 {
                    match ax {
                        0 => 0.5 * (f[idx(x + 1, y, z)] - f[idx(x - 1, y, z)]),
                        1 => 0.5 * (f[idx(x, y + 1, z)] - f[idx(x, y - 1, z)]),
                        _ => 0.5 * (f[idx(x, y, z + 1)] - f[idx(x, y, z - 1)]),
                    }
                };
                let (uxx, uxy, uxz) = (d(ux, 0), d(ux, 1), d(ux, 2));
                let (uyx, uyy, uyz) = (d(uy, 0), d(uy, 1), d(uy, 2));
                let (uzx, uzy, uzz) = (d(uz, 0), d(uz, 1), d(uz, 2));
                let sxy = 0.5 * (uxy + uyx);
                let sxz = 0.5 * (uxz + uzx);
                let syz = 0.5 * (uyz + uzy);
                let ss = uxx * uxx + uyy * uyy + uzz * uzz
                    + 2.0 * (sxy * sxy + sxz * sxz + syz * syz);
                out[idx(x, y, z)] = nu * (2.0 * ss).sqrt();
            }
        }
    }
    out
}

/// Subsample the domain to a `vn^3` grid and write vx,vy,vz,shear as f32 LE
/// plus a JSON meta block (all geometry in lattice cells).
#[allow(clippy::too_many_arguments)]
fn export_volume(outdir: &Path, g: &Geom, ux: &[f64], uy: &[f64], uz: &[f64], shear: &[f64],
                 nx: usize, ny: usize, nz: usize, omega: f64, u_tip: f64, nu: f64) {
    let idx = |x: usize, y: usize, z: usize| (z * ny + y) * nx + x;
    let vn = 60usize.min(nx);
    let map = |i: usize, m: usize| ((i as f64) * (m as f64 - 1.0) / (vn as f64 - 1.0)).round() as usize;
    let mut buf: Vec<f32> = Vec::with_capacity(vn * vn * vn * 4);
    let mut smax = 0.0f64;
    let mut shmax = 0.0f64;
    for k in 0..vn {
        let z = map(k, nz);
        for j in 0..vn {
            let y = map(j, ny);
            for i in 0..vn {
                let x = map(i, nx);
                let c = idx(x, y, z);
                buf.push(ux[c] as f32);
                buf.push(uy[c] as f32);
                buf.push(uz[c] as f32);
                buf.push(shear[c] as f32);
                smax = smax.max((ux[c] * ux[c] + uy[c] * uy[c] + uz[c] * uz[c]).sqrt());
                shmax = shmax.max(shear[c]);
            }
        }
    }
    let mut f = BufWriter::new(File::create(outdir.join("volume.bin")).expect("volume.bin"));
    for v in &buf {
        f.write_all(&v.to_le_bytes()).expect("write f32");
    }
    f.flush().ok();

    let meta = format!(
        concat!(
            "{{\"vn\": {}, \"n\": {}, \"layout\": \"x-fastest,[vx,vy,vz,shear] f32 LE\", ",
            "\"cx\": {}, \"cy\": {}, \"r_tank\": {}, \"zc\": {}, ",
            "\"tip_r\": {}, \"disk_r\": {}, \"hub_r\": {}, \"shaft_r\": {}, ",
            "\"blade_hh\": {}, \"disk_hh\": {}, \"n_blades\": {}, \"blade_hw\": {}, ",
            "\"baffle_len\": {}, \"baffle_hw\": {}, ",
            "\"omega\": {}, \"u_tip\": {}, \"nu\": {}, \"speed_max\": {}, \"shear_max\": {}}}\n"
        ),
        vn, g.n, g.cx, g.cy, g.r_tank, g.zc,
        g.tip_r, g.disk_r, g.hub_r, g.shaft_r,
        g.blade_hh, g.disk_hh, g.n_blades, g.blade_hw,
        g.baffle_len, g.baffle_hw, omega, u_tip, nu, smax, shmax
    );
    std::fs::write(outdir.join("volume.json"), meta).expect("volume.json");
    println!("volume: {vn}^3 (speed_max={smax:.4}, shear_max={shmax:.5})");
}
