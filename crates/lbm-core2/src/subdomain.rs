//! Subdomain: one part of a decomposed global grid (docs/ARCHITECTURE_V2.md
//! §2.3).
//!
//! A subdomain owns a `core` box of global cells plus a one-cell halo ring on
//! every active axis. Face neighbours are part indices (later: MPI ranks);
//! `None` marks a non-periodic global boundary, where streaming leaves the
//! unknown slots untouched (V1 open-edge mechanics). Corner/edge halo data is
//! forwarded by the two-phase exchange (x, then y including x-halos, then z
//! including x/y-halos), so only face neighbours are ever stored.

use crate::fields::LocalGeom;
use crate::lattice::Face;

/// One part of the decomposed global grid.
#[derive(Clone, Debug)]
pub struct Subdomain {
    /// Global grid extents `[nx, ny, nz]` (`nz == 1` for 2D).
    pub global: [usize; 3],
    /// Global coordinate of this part's core cell (0, 0, 0).
    pub origin: [usize; 3],
    /// Local core extents + halo width.
    pub geom: LocalGeom,
    /// Neighbouring part index per face (`Face::index()` order). `Some(self)`
    /// encodes a periodic wrap onto the same part.
    pub neighbors: [Option<usize>; 6],
}

impl Subdomain {
    /// Single-part "decomposition": the whole global grid as one subdomain,
    /// with periodic axes wrapping onto itself. Reproduces V1 behaviour.
    pub fn monolithic(d: usize, dims: [usize; 3], periodic: [bool; 3]) -> Self {
        let mut neighbors = [None; 6];
        for face in Face::ALL {
            if face.axis() < d && periodic[face.axis()] {
                neighbors[face.index()] = Some(0);
            }
        }
        Self {
            global: dims,
            origin: [0, 0, 0],
            geom: LocalGeom::new(d, dims, 1),
            neighbors,
        }
    }

    /// Whether streaming may pull from the halo behind `face` (i.e. a
    /// neighbour fills it). `false` at non-periodic global boundaries.
    #[inline]
    pub fn has_halo(&self, face: Face) -> bool {
        self.neighbors[face.index()].is_some()
    }

    /// Per-face halo availability in `Face::index()` order (kernel input).
    #[inline]
    pub fn halo_flags(&self) -> [bool; 6] {
        let mut f = [false; 6];
        for face in Face::ALL {
            f[face.index()] = self.has_halo(face);
        }
        f
    }

    /// Whether this part's core touches the given global boundary face
    /// (where global open-face BCs apply).
    pub fn touches_global_face(&self, face: Face) -> bool {
        let a = face.axis();
        if a >= self.geom.d {
            return false;
        }
        if face.is_neg() {
            self.origin[a] == 0
        } else {
            self.origin[a] + self.geom.core[a] == self.global[a]
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn monolithic_periodicity_maps_to_halo_flags() {
        let s = Subdomain::monolithic(2, [8, 4, 1], [true, false, false]);
        assert!(s.has_halo(Face::XNeg) && s.has_halo(Face::XPos));
        assert!(!s.has_halo(Face::YNeg) && !s.has_halo(Face::YPos));
        assert!(!s.has_halo(Face::ZNeg));
        for f in Face::ALL {
            let touches = s.touches_global_face(f);
            assert_eq!(touches, f.axis() < 2, "{f:?}");
        }
    }
}
