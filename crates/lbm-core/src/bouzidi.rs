//! Bouzidi-Firdaouss-Lallemand interpolated bounce-back records and CPU pass.
//!
//! Implemented from the original Bouzidi et al. 2001 interpolated
//! bounce-back formula, with the same Guo/Ladd moving-wall term used by the
//! existing half-way bounce-back path. No GPL/AGPL implementation was used.

use crate::fields::{LocalGeom, SoaFields};
use crate::lattice::Lattice;
use crate::params::{KParams, StepParams};
use crate::real::Real;

/// One curved-wall link from a fluid cell toward a wall intersection.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BouzidiLink<T: Real> {
    /// Padded local index of the fluid cell.
    pub cell: u32,
    /// Direction from the fluid cell toward the wall.
    pub q: u8,
    /// Fractional wall distance along `q`, in `(0, 1)`.
    pub qd: T,
    /// Whether the second fluid node `x_f - c_q` exists and is fluid.
    pub has_second: bool,
    /// Padded local index of the wall/solid neighbour, used for wall velocity
    /// and probe membership.
    pub wall_ref: u32,
}

/// Sorted Bouzidi wall-distance records.
#[derive(Clone, Debug, PartialEq)]
pub struct BouzidiLinks<T: Real> {
    /// One record per fluid-cell/wall-crossing lattice link.
    pub records: Vec<BouzidiLink<T>>,
}

impl<T: Real> BouzidiLinks<T> {
    pub fn new(mut records: Vec<BouzidiLink<T>>) -> Self {
        records.sort_by_key(|r| (r.cell, r.q));
        records.dedup_by_key(|r| (r.cell, r.q));
        Self { records }
    }

    pub fn is_empty(&self) -> bool {
        self.records.is_empty()
    }
}

/// Build links for an analytic circle in a 2D local subdomain. The caller is
/// responsible for having already marked the corresponding solid cells.
pub(crate) fn circle_links<T: Real>(
    geom: &LocalGeom,
    origin: [usize; 3],
    solid: &[bool],
    cx: f64,
    cy: f64,
    r: f64,
) -> BouzidiLinks<T> {
    let mut out = Vec::new();
    let r2 = r * r;
    for y in 0..geom.core[1] {
        for x in 0..geom.core[0] {
            let cell = geom.pidx(x, y, 0);
            if solid[cell] {
                continue;
            }
            let gx = (origin[0] + x) as f64;
            let gy = (origin[1] + y) as f64;
            for q in 1..crate::lattice::D2Q9::Q {
                let c = crate::lattice::D2Q9::C[q];
                let sx = x as isize + c[0] as isize;
                let sy = y as isize + c[1] as isize;
                if sx < -(geom.halo as isize)
                    || sy < -(geom.halo as isize)
                    || sx >= geom.core[0] as isize + geom.halo as isize
                    || sy >= geom.core[1] as isize + geom.halo as isize
                {
                    continue;
                }
                let wall_ref = geom.pidx_i(sx, sy, 0);
                if !solid[wall_ref] {
                    continue;
                }
                let Some(qd) = ray_circle_qd(gx, gy, c[0] as f64, c[1] as f64, cx, cy, r2) else {
                    continue;
                };
                let bx = x as isize - c[0] as isize;
                let by = y as isize - c[1] as isize;
                let has_second = bx >= -(geom.halo as isize)
                    && by >= -(geom.halo as isize)
                    && bx < geom.core[0] as isize + geom.halo as isize
                    && by < geom.core[1] as isize + geom.halo as isize
                    && !solid[geom.pidx_i(bx, by, 0)];
                out.push(BouzidiLink {
                    cell: cell as u32,
                    q: q as u8,
                    qd: T::r(qd),
                    has_second,
                    wall_ref: wall_ref as u32,
                });
            }
        }
    }
    BouzidiLinks::new(out)
}

/// Build links for an analytic sphere in a 3D local subdomain.
pub(crate) fn sphere_links<T: Real, L: Lattice>(
    geom: &LocalGeom,
    origin: [usize; 3],
    solid: &[bool],
    cx: f64,
    cy: f64,
    cz: f64,
    r: f64,
) -> BouzidiLinks<T> {
    let mut out = Vec::new();
    let r2 = r * r;
    for z in 0..geom.core[2] {
        for y in 0..geom.core[1] {
            for x in 0..geom.core[0] {
                let cell = geom.pidx(x, y, z);
                if solid[cell] {
                    continue;
                }
                let gp = [
                    (origin[0] + x) as f64,
                    (origin[1] + y) as f64,
                    (origin[2] + z) as f64,
                ];
                for q in 1..L::Q {
                    let c = L::C[q];
                    let sp = [
                        x as isize + c[0] as isize,
                        y as isize + c[1] as isize,
                        z as isize + c[2] as isize,
                    ];
                    if (0..L::D).any(|a| {
                        sp[a] < -(geom.halo as isize)
                            || sp[a] >= geom.core[a] as isize + geom.halo as isize
                    }) {
                        continue;
                    }
                    let wall_ref = geom.pidx_i(sp[0], sp[1], sp[2]);
                    if !solid[wall_ref] {
                        continue;
                    }
                    let Some(qd) = ray_sphere_qd(
                        gp,
                        [c[0] as f64, c[1] as f64, c[2] as f64],
                        [cx, cy, cz],
                        r2,
                    ) else {
                        continue;
                    };
                    let bp = [
                        x as isize - c[0] as isize,
                        y as isize - c[1] as isize,
                        z as isize - c[2] as isize,
                    ];
                    let has_second = !(0..L::D).any(|a| {
                        bp[a] < -(geom.halo as isize)
                            || bp[a] >= geom.core[a] as isize + geom.halo as isize
                    }) && !solid[geom.pidx_i(bp[0], bp[1], bp[2])];
                    out.push(BouzidiLink {
                        cell: cell as u32,
                        q: q as u8,
                        qd: T::r(qd),
                        has_second,
                        wall_ref: wall_ref as u32,
                    });
                }
            }
        }
    }
    BouzidiLinks::new(out)
}

/// Build qd=1/2 records for every fluid-solid link. This is a strict
/// degeneracy harness: the pass must reproduce half-way bounce-back bitwise.
pub(crate) fn half_way_links<T: Real, L: Lattice>(
    geom: &LocalGeom,
    solid: &[bool],
) -> BouzidiLinks<T> {
    let mut out = Vec::new();
    for z in 0..geom.core[2] {
        for y in 0..geom.core[1] {
            for x in 0..geom.core[0] {
                let cell = geom.pidx(x, y, z);
                if solid[cell] {
                    continue;
                }
                for q in 1..L::Q {
                    let c = L::C[q];
                    let sp = [
                        x as isize + c[0] as isize,
                        y as isize + c[1] as isize,
                        z as isize + c[2] as isize,
                    ];
                    if (0..L::D).any(|a| {
                        sp[a] < -(geom.halo as isize)
                            || sp[a] >= geom.core[a] as isize + geom.halo as isize
                    }) {
                        continue;
                    }
                    let wall_ref = geom.pidx_i(sp[0], sp[1], sp[2]);
                    if solid[wall_ref] {
                        out.push(BouzidiLink {
                            cell: cell as u32,
                            q: q as u8,
                            qd: T::r(0.5),
                            has_second: false,
                            wall_ref: wall_ref as u32,
                        });
                    }
                }
            }
        }
    }
    BouzidiLinks::new(out)
}

fn ray_circle_qd(x: f64, y: f64, dx: f64, dy: f64, cx: f64, cy: f64, r2: f64) -> Option<f64> {
    let ox = x - cx;
    let oy = y - cy;
    let a = dx * dx + dy * dy;
    let b = 2.0 * (ox * dx + oy * dy);
    let c = ox * ox + oy * oy - r2;
    let disc = b * b - 4.0 * a * c;
    if disc < 0.0 {
        return None;
    }
    let s = disc.sqrt();
    [(-b - s) / (2.0 * a), (-b + s) / (2.0 * a)]
        .into_iter()
        .filter(|&t| t > 0.0 && t < 1.0)
        .min_by(|a, b| a.total_cmp(b))
}

fn ray_sphere_qd(p: [f64; 3], d: [f64; 3], c0: [f64; 3], r2: f64) -> Option<f64> {
    let o = [p[0] - c0[0], p[1] - c0[1], p[2] - c0[2]];
    let a = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
    let b = 2.0 * (o[0] * d[0] + o[1] * d[1] + o[2] * d[2]);
    let c = o[0] * o[0] + o[1] * o[1] + o[2] * o[2] - r2;
    let disc = b * b - 4.0 * a * c;
    if disc < 0.0 {
        return None;
    }
    let s = disc.sqrt();
    [(-b - s) / (2.0 * a), (-b + s) / (2.0 * a)]
        .into_iter()
        .filter(|&t| t > 0.0 && t < 1.0)
        .min_by(|a, b| a.total_cmp(b))
}

/// Apply Bouzidi records after streaming and before the population swap.
/// Returns the probe-force replacement delta for curved links.
pub(crate) fn apply_bouzidi_impl<L: Lattice, T: Real>(
    fields: &mut SoaFields<T>,
    p: &StepParams<T>,
) -> [T; 3] {
    let Some(links) = fields.bouzidi.as_ref() else {
        return [T::zero(); 3];
    };
    if links.records.is_empty() {
        return [T::zero(); 3];
    }
    let kp = KParams::new::<L>(p);
    let np = fields.geom.n_padded();
    let two = T::r(2.0);
    let one = T::one();
    let six = T::r(6.0);
    let half = T::r(0.5);
    let mut delta = [T::zero(); 3];
    for rec in &links.records {
        let cell = rec.cell as usize;
        let wall = rec.wall_ref as usize;
        let q = rec.q as usize;
        let qb = L::OPP[q];
        let rho = fields.rho[compact_from_padded(fields.geom, cell)];
        let wu = fields.wall_u[wall];
        let mut cu = kp.cr[qb][0] * wu[0];
        for d in 1..L::D {
            cu = cu + kp.cr[qb][d] * wu[d];
        }
        let wall_term = six * kp.wr[q] * rho * cu;
        let fq = fields.f[q * np + cell];
        let fin = if rec.qd == half {
            fq + wall_term
        } else if rec.qd < half && rec.has_second {
            let c = L::C[q];
            let ff = fields.geom.pidx_i(
                padded_coord(fields.geom, cell, 0) - c[0] as isize,
                padded_coord(fields.geom, cell, 1) - c[1] as isize,
                padded_coord(fields.geom, cell, 2) - c[2] as isize,
            );
            let a = two * rec.qd;
            // BFL qd<1/2 reconstructs the post-reflection value at x_f by
            // interpolating between the wall-reflected population at x_f and
            // the same reflected population at the second fluid node x_f-c_q:
            //
            //   sigma f_q(x_f) + (1-sigma) f_q(x_f-c_q)
            //     + [sigma + (1-sigma)] 2 w_q rho (c_opp . u_w) / cs^2
            //
            // with sigma = 2 qd and cs^2 = 1/3. Applying the moving-wall
            // correction only to the first interpolation point would impose
            // sigma*u_w. The second-point term supplies the missing
            // (1-sigma) share, so the wall velocity is u_w for every qd<1/2.
            a * (fq + wall_term) + (one - a) * (fields.f[q * np + ff] + wall_term)
        } else {
            let a = one / (two * rec.qd);
            a * fq + (one - a) * fields.f[qb * np + cell] + a * wall_term
        };
        fields.ftmp[qb * np + cell] = fin;

        if fields.probe.as_ref().is_some_and(|m| m[wall]) {
            let hw_fin = fq + wall_term;
            let hw_tot = fq + hw_fin + two * kp.wr[q];
            let bz_tot = fq + fin + two * kp.wr[q];
            let dft = bz_tot - hw_tot;
            for a in 0..L::D {
                delta[a] = delta[a] + kp.cr[q][a] * dft;
            }
        }
    }
    if let Some(s) = fields.fused.as_deref_mut() {
        s.fresh = false;
    }
    delta
}

#[inline]
fn padded_coord(g: LocalGeom, cell: usize, axis: usize) -> isize {
    let p = g.padded();
    let raw = match axis {
        0 => cell % p[0],
        1 => (cell / p[0]) % p[1],
        _ => cell / (p[0] * p[1]),
    };
    let h = if axis < g.d { g.halo as isize } else { 0 };
    raw as isize - h
}

#[inline]
fn compact_from_padded(g: LocalGeom, cell: usize) -> usize {
    let x = padded_coord(g, cell, 0) as usize;
    let y = padded_coord(g, cell, 1) as usize;
    let z = padded_coord(g, cell, 2) as usize;
    g.cidx(x, y, z)
}
