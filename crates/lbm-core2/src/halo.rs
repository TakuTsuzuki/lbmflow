//! Halo exchange (docs/ARCHITECTURE_V2.md §2.3).
//!
//! Streaming pulls from a one-cell halo ring; the exchange fills that ring
//! from the neighbouring parts after collision, before streaming:
//!
//! - **populations**: per face only the directions *entering* through it
//!   (`L::unknowns(face)`: 3 for D2Q9, 5 for D3Q19) need transferring;
//! - **masks** (`solid` / `wall_u` / `probe`): exchanged when geometry
//!   changes (bounce-back reads the wall data of halo cells);
//! - **scalar planes** (multiphase ψ): full-value exchange for force stencils.
//!
//! Corner/edge halo cells are *forwarded*, not exchanged diagonally: phases
//! run x → y → z, and each later phase transfers layers extended over the
//! earlier axes' halos (the standard two-phase trick). A corner value thus
//! hops through the face neighbour, and only 6 face links exist — the same
//! plan an MPI implementation uses.
//!
//! The transfer itself is pack → unpack through a contiguous buffer, i.e.
//! exactly message-shaped: the future `Mpi` implementation replaces the
//! buffer hand-off with send/recv and keeps the layer maths.

use crate::fields::{LocalGeom, SoaFields};
use crate::lattice::{Face, Lattice};
use crate::real::Real;
use crate::subdomain::Subdomain;

/// Fills halo rings from neighbouring parts. Implementations define where
/// neighbours live (same process, other threads, MPI ranks).
///
/// GPU note: backends with device-resident fields provide their own
/// implementation over device buffers; the layer geometry below is layout-
/// identical, so the plan (faces, phases, direction sets) is shared.
pub trait HaloExchange<T: Real> {
    /// Fill the population halos (post-collide, pre-stream).
    fn exchange_f<L: Lattice>(&self, subs: &[Subdomain], parts: &mut [SoaFields<T>]);
    /// Refresh halo copies of `solid` / `wall_u` / `probe` after edits.
    fn exchange_masks(&self, subs: &[Subdomain], parts: &mut [SoaFields<T>]);
    /// Exchange one padded scalar plane per part (e.g. multiphase ψ).
    fn exchange_scalar(&self, subs: &[Subdomain], planes: &mut [&mut [T]]);
}

/// Single-part exchange: periodic axes wrap onto the same part. This is the
/// V1-equivalent configuration (one monolithic domain).
#[derive(Clone, Copy, Debug, Default)]
pub struct LocalPeriodic;

impl<T: Real> HaloExchange<T> for LocalPeriodic {
    fn exchange_f<L: Lattice>(&self, subs: &[Subdomain], parts: &mut [SoaFields<T>]) {
        assert_eq!(parts.len(), 1, "LocalPeriodic serves a single part");
        exchange_f_generic::<L, T>(subs, parts);
    }

    fn exchange_masks(&self, subs: &[Subdomain], parts: &mut [SoaFields<T>]) {
        assert_eq!(parts.len(), 1, "LocalPeriodic serves a single part");
        exchange_masks_generic(subs, parts);
    }

    fn exchange_scalar(&self, subs: &[Subdomain], planes: &mut [&mut [T]]) {
        assert_eq!(planes.len(), 1, "LocalPeriodic serves a single part");
        exchange_scalar_generic(subs, planes);
    }
}

// ---------------------------------------------------------------------------
// Shared layer machinery (used by LocalPeriodic and InProcess)
// ---------------------------------------------------------------------------

/// Padded indices of the halo layer behind `recv_face` (unpack side), in
/// canonical order. Layers of phase `axis` extend over the halos of earlier
/// axes (< `axis`), which were exchanged in earlier phases.
fn layer_indices(geom: &LocalGeom, recv_face: Face, phase_axis: usize, unpack: bool) -> Vec<usize> {
    let a = recv_face.axis();
    debug_assert_eq!(a, phase_axis);
    let h = geom.halo as isize;
    // Unpack writes the receiver's halo layer; pack reads the sender's
    // opposite core boundary layer.
    let fixed: isize = match (unpack, recv_face.is_neg()) {
        (true, true) => -h,                              // receiver low halo
        (true, false) => geom.core[a] as isize,          // receiver high halo
        (false, true) => geom.core[a] as isize - 1,      // sender high core
        (false, false) => 0,                             // sender low core
    };
    let range = |t: usize| -> (isize, isize) {
        if t < phase_axis && t < geom.d {
            (-h, geom.core[t] as isize + h) // extended: forwards corners
        } else {
            (0, geom.core[t] as isize)
        }
    };
    let (t1, t2) = match a {
        0 => (1, 2),
        1 => (0, 2),
        _ => (0, 1),
    };
    let (r1, r2) = (range(t1), range(t2));
    let mut out = Vec::with_capacity(((r1.1 - r1.0) * (r2.1 - r2.0)) as usize);
    for c2 in r2.0..r2.1 {
        for c1 in r1.0..r1.1 {
            let mut pos = [0isize; 3];
            pos[a] = fixed;
            pos[t1] = c1;
            pos[t2] = c2;
            out.push(geom.pidx_i(pos[0], pos[1], pos[2]));
        }
    }
    out
}

/// Assert the Cartesian-decomposition invariant: sender and receiver share
/// tangent extents, so layer cells map 1:1.
fn check_tangent_match(a: usize, dst: &LocalGeom, src: &LocalGeom) {
    for t in 0..3 {
        if t != a {
            assert_eq!(
                dst.core[t], src.core[t],
                "non-Cartesian decomposition: tangent extents differ on axis {t}"
            );
        }
    }
}

/// Generic post-collide population exchange over any part set whose
/// neighbours live in `parts` (in-process or single-part periodic).
pub(crate) fn exchange_f_generic<L: Lattice, T: Real>(
    subs: &[Subdomain],
    parts: &mut [SoaFields<T>],
) {
    let d = subs[0].geom.d;
    let mut buf: Vec<T> = Vec::new();
    for axis in 0..d {
        for side in 0..2 {
            let recv_face = Face::ALL[2 * axis + side];
            let dirs = L::unknowns(recv_face);
            for di in 0..parts.len() {
                let Some(si) = subs[di].neighbors[recv_face.index()] else {
                    continue;
                };
                check_tangent_match(axis, &subs[di].geom, &subs[si].geom);
                // pack (sender's opposite core boundary layer)
                let src_idx = layer_indices(&parts[si].geom, recv_face, axis, false);
                let np_s = parts[si].plane_len();
                buf.clear();
                buf.reserve(src_idx.len() * dirs.len());
                for &cell in &src_idx {
                    for &q in dirs {
                        buf.push(parts[si].f[q * np_s + cell]);
                    }
                }
                // unpack (receiver's halo layer)
                let dst_idx = layer_indices(&parts[di].geom, recv_face, axis, true);
                debug_assert_eq!(dst_idx.len(), src_idx.len());
                let np_d = parts[di].plane_len();
                let mut k = 0;
                for &cell in &dst_idx {
                    for &q in dirs {
                        parts[di].f[q * np_d + cell] = buf[k];
                        k += 1;
                    }
                }
            }
        }
    }
}

/// Generic mask exchange (`solid`, `wall_u`, `probe`).
pub(crate) fn exchange_masks_generic<T: Real>(subs: &[Subdomain], parts: &mut [SoaFields<T>]) {
    let d = subs[0].geom.d;
    for axis in 0..d {
        for side in 0..2 {
            let recv_face = Face::ALL[2 * axis + side];
            for di in 0..parts.len() {
                let Some(si) = subs[di].neighbors[recv_face.index()] else {
                    continue;
                };
                check_tangent_match(axis, &subs[di].geom, &subs[si].geom);
                let src_idx = layer_indices(&parts[si].geom, recv_face, axis, false);
                let dst_idx = layer_indices(&parts[di].geom, recv_face, axis, true);
                debug_assert_eq!(dst_idx.len(), src_idx.len());
                let solid_buf: Vec<bool> = src_idx.iter().map(|&c| parts[si].solid[c]).collect();
                let wall_buf: Vec<[T; 3]> = src_idx.iter().map(|&c| parts[si].wall_u[c]).collect();
                let probe_buf: Option<Vec<bool>> = parts[si]
                    .probe
                    .as_ref()
                    .map(|m| src_idx.iter().map(|&c| m[c]).collect());
                for (k, &cell) in dst_idx.iter().enumerate() {
                    parts[di].solid[cell] = solid_buf[k];
                    parts[di].wall_u[cell] = wall_buf[k];
                }
                if let Some(pb) = probe_buf {
                    let dst = parts[di]
                        .probe
                        .as_mut()
                        .expect("probe mask must be materialised on every part");
                    for (k, &cell) in dst_idx.iter().enumerate() {
                        dst[cell] = pb[k];
                    }
                }
            }
        }
    }
}

/// Generic scalar-plane exchange (padded planes, full values).
pub(crate) fn exchange_scalar_generic<T: Real>(subs: &[Subdomain], planes: &mut [&mut [T]]) {
    let d = subs[0].geom.d;
    for axis in 0..d {
        for side in 0..2 {
            let recv_face = Face::ALL[2 * axis + side];
            for di in 0..planes.len() {
                let Some(si) = subs[di].neighbors[recv_face.index()] else {
                    continue;
                };
                check_tangent_match(axis, &subs[di].geom, &subs[si].geom);
                let src_idx = layer_indices(&subs[si].geom, recv_face, axis, false);
                let dst_idx = layer_indices(&subs[di].geom, recv_face, axis, true);
                debug_assert_eq!(dst_idx.len(), src_idx.len());
                let buf: Vec<T> = src_idx.iter().map(|&c| planes[si][c]).collect();
                for (k, &cell) in dst_idx.iter().enumerate() {
                    planes[di][cell] = buf[k];
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::lattice::D2Q9;

    /// Double-periodic single part: after exchange, every halo cell that
    /// streaming can read must hold the wrapped core value.
    #[test]
    fn local_periodic_fills_wrapped_halos() {
        let dims = [5usize, 4, 1];
        let sub = Subdomain::monolithic(2, dims, [true, true, false]);
        let mut fields: SoaFields<f64> = SoaFields::new(9, sub.geom);
        let np = fields.plane_len();
        // Tag every core cell with a unique value per direction.
        for q in 0..9 {
            for y in 0..4 {
                for x in 0..5 {
                    fields.f[q * np + sub.geom.pidx(x, y, 0)] =
                        (q * 100 + y * 10 + x) as f64;
                }
            }
        }
        let ex = LocalPeriodic;
        HaloExchange::<f64>::exchange_f::<D2Q9>(&ex, &[sub.clone()], std::slice::from_mut(&mut fields));
        let g = &sub.geom;
        let wrap = |v: isize, n: usize| ((v + n as isize) % n as isize) as usize;
        // Streaming reads halo cell s with direction q iff s + c_q is a core
        // cell, i.e. q enters through the crossed face(s).
        for q in 0..9 {
            let c = crate::lattice::D2Q9::C[q];
            for y in -1isize..5 {
                for x in -1isize..6 {
                    let in_halo = x < 0 || x >= 5 || y < 0 || y >= 4;
                    if !in_halo {
                        continue;
                    }
                    let (dx, dy) = (x + c[0] as isize, y + c[1] as isize);
                    let dest_is_core = dx >= 0 && dx < 5 && dy >= 0 && dy < 4;
                    if !dest_is_core {
                        continue; // this (cell, dir) is never pulled
                    }
                    let (wx, wy) = (wrap(x, 5), wrap(y, 4));
                    let got = fields.f[q * np + g.pidx_i(x, y, 0)];
                    let want = fields.f[q * np + g.pidx(wx, wy, 0)];
                    assert_eq!(got, want, "q={q} halo=({x},{y}) wrap=({wx},{wy})");
                }
            }
        }
    }

    #[test]
    fn mask_exchange_wraps_solids() {
        let dims = [4usize, 3, 1];
        let sub = Subdomain::monolithic(2, dims, [true, true, false]);
        let mut fields: SoaFields<f64> = SoaFields::new(9, sub.geom);
        fields.solid[sub.geom.pidx(3, 1, 0)] = true;
        fields.wall_u[sub.geom.pidx(3, 1, 0)] = [0.1, 0.0, 0.0];
        let ex = LocalPeriodic;
        HaloExchange::<f64>::exchange_masks(&ex, &[sub.clone()], std::slice::from_mut(&mut fields));
        assert!(fields.solid[sub.geom.pidx_i(-1, 1, 0)]);
        assert_eq!(fields.wall_u[sub.geom.pidx_i(-1, 1, 0)], [0.1, 0.0, 0.0]);
        assert!(!fields.solid[sub.geom.pidx_i(-1, 0, 0)]);
    }
}
