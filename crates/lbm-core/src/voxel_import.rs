//! Minimal binary-STL voxel import for bioprocess geometry.

use crate::geometry::{GeometryError, GeometryResult};
use serde::Deserialize;
use std::fs;
use std::path::Path;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CredibilityTier {
    Screening,
    Engineering,
    Evidence,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum PatchLabel {
    Wall,
    Impeller,
    Baffle,
    Sparger,
    Unknown,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct VoxelImportOptions {
    pub dims: [usize; 3],
    pub dx_m: f64,
    pub credibility_tier: CredibilityTier,
}

#[derive(Clone, Debug, PartialEq)]
pub struct ImportedVoxelGeometry {
    pub solid: Vec<bool>,
    pub wall_mask: Vec<bool>,
    pub impeller_mask: Vec<bool>,
    pub baffle_mask: Vec<bool>,
    pub sparger_mask: Vec<bool>,
    pub labels: Vec<PatchLabel>,
}

#[derive(Clone, Copy, Debug)]
struct Triangle {
    v: [[f64; 3]; 3],
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct PatchFile {
    patches: Vec<PatchSpec>,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct PatchSpec {
    name: PatchLabel,
    faces: Vec<usize>,
}

pub fn import_binary_stl_with_labels(
    stl_path: &Path,
    labels_path: &Path,
    options: VoxelImportOptions,
) -> GeometryResult<ImportedVoxelGeometry> {
    validate_options(options)?;
    let triangles = read_binary_stl(stl_path)?;
    let labels = read_patch_labels(labels_path, triangles.len(), options.credibility_tier)?;
    Ok(voxelize(&triangles, &labels, options))
}

fn validate_options(options: VoxelImportOptions) -> GeometryResult<()> {
    if options.dims.iter().any(|&n| n == 0) || !(options.dx_m.is_finite() && options.dx_m > 0.0) {
        return Err(GeometryError::out_of_validity_range(
            "voxel import grid must be finite and positive",
            "dims must be positive and dx_m must be finite and > 0",
        ));
    }
    Ok(())
}

fn read_binary_stl(path: &Path) -> GeometryResult<Vec<Triangle>> {
    let bytes = fs::read(path).map_err(|e| {
        GeometryError::out_of_validity_range(
            "cannot read STL file",
            format!("{}: {e}", path.display()),
        )
    })?;
    if bytes.len() < 84 {
        return Err(GeometryError::out_of_validity_range(
            "STL file is too short",
            "binary STL requires an 80-byte header and 4-byte face count",
        ));
    }
    if bytes
        .get(0..5)
        .is_some_and(|prefix| prefix.eq_ignore_ascii_case(b"solid"))
    {
        return Err(GeometryError::not_implemented(
            "ASCII STL import is not implemented; provide binary STL",
        ));
    }
    let n_tri = u32::from_le_bytes(bytes[80..84].try_into().expect("slice length")) as usize;
    let expected = 84usize.saturating_add(n_tri.saturating_mul(50));
    if bytes.len() != expected {
        return Err(GeometryError::out_of_validity_range(
            "binary STL byte length does not match its face count",
            format!(
                "expected {expected} bytes for {n_tri} faces, got {}",
                bytes.len()
            ),
        ));
    }
    let mut triangles = Vec::with_capacity(n_tri);
    for face in 0..n_tri {
        let base = 84 + face * 50 + 12;
        let mut v = [[0.0; 3]; 3];
        for (corner, dst) in v.iter_mut().enumerate() {
            for (axis, comp) in dst.iter_mut().enumerate() {
                let off = base + corner * 12 + axis * 4;
                *comp = f32::from_le_bytes(bytes[off..off + 4].try_into().expect("slice length"))
                    as f64;
            }
        }
        triangles.push(Triangle { v });
    }
    Ok(triangles)
}

fn read_patch_labels(
    path: &Path,
    face_count: usize,
    tier: CredibilityTier,
) -> GeometryResult<Vec<PatchLabel>> {
    let text = fs::read_to_string(path).map_err(|e| {
        GeometryError::out_of_validity_range(
            "cannot read STL patch-label file",
            format!("{}: {e}", path.display()),
        )
    })?;
    let patches: PatchFile = serde_json::from_str(&text).map_err(|e| {
        GeometryError::out_of_validity_range(
            "invalid STL patch-label JSON",
            format!("{}: {e}", path.display()),
        )
    })?;
    let mut labels = vec![PatchLabel::Unknown; face_count];
    let mut labelled = vec![false; face_count];
    for patch in patches.patches {
        if patch.name == PatchLabel::Unknown && tier != CredibilityTier::Screening {
            return Err(GeometryError::out_of_validity_range(
                "unknown STL patch labels are allowed only for screening tier",
                "patch label unknown is not accepted for engineering/evidence tier",
            ));
        }
        for face in patch.faces {
            if face >= face_count {
                return Err(GeometryError::out_of_validity_range(
                    "STL patch label references a face outside the triangle list",
                    format!("face index {face} >= {face_count}"),
                ));
            }
            labels[face] = patch.name;
            labelled[face] = true;
        }
    }
    if tier == CredibilityTier::Evidence && labelled.iter().any(|&v| !v) {
        let missing = labelled
            .iter()
            .enumerate()
            .filter_map(|(i, labelled)| {
                (!labelled).then(|| format!("face {i} missing patch label"))
            })
            .collect();
        return Err(GeometryError::evidence_gate_failed(
            "evidence-tier STL import requires explicit labels for every face",
            missing,
        ));
    }
    Ok(labels)
}

fn voxelize(
    triangles: &[Triangle],
    labels: &[PatchLabel],
    options: VoxelImportOptions,
) -> ImportedVoxelGeometry {
    let n = options.dims[0] * options.dims[1] * options.dims[2];
    let mut out = ImportedVoxelGeometry {
        solid: vec![false; n],
        wall_mask: vec![false; n],
        impeller_mask: vec![false; n],
        baffle_mask: vec![false; n],
        sparger_mask: vec![false; n],
        labels: labels.to_vec(),
    };
    for z in 0..options.dims[2] {
        for y in 0..options.dims[1] {
            for x in 0..options.dims[0] {
                let p = [
                    (x as f64 + 0.5) * options.dx_m,
                    (y as f64 + 0.5) * options.dx_m,
                    (z as f64 + 0.5) * options.dx_m,
                ];
                let Some(label) = point_label_inside(p, triangles, labels) else {
                    continue;
                };
                let i = (z * options.dims[1] + y) * options.dims[0] + x;
                match label {
                    PatchLabel::Sparger => out.sparger_mask[i] = true,
                    PatchLabel::Impeller => out.impeller_mask[i] = true,
                    PatchLabel::Baffle => {
                        out.solid[i] = true;
                        out.baffle_mask[i] = true;
                    }
                    PatchLabel::Wall | PatchLabel::Unknown => {
                        out.solid[i] = true;
                        out.wall_mask[i] = true;
                    }
                }
            }
        }
    }
    out
}

fn point_label_inside(
    p: [f64; 3],
    triangles: &[Triangle],
    labels: &[PatchLabel],
) -> Option<PatchLabel> {
    let dir = [1.0, 0.0, 0.0];
    let mut hits: Vec<(f64, PatchLabel)> = triangles
        .iter()
        .zip(labels.iter())
        .filter_map(|(tri, &label)| ray_intersects_triangle(p, dir, tri).map(|t| (t, label)))
        .collect();
    hits.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
    hits.dedup_by(|a, b| (a.0 - b.0).abs() < 1.0e-9);
    if hits.len() % 2 == 1 {
        hits.first().map(|(_, label)| *label)
    } else {
        None
    }
}

fn ray_intersects_triangle(origin: [f64; 3], dir: [f64; 3], tri: &Triangle) -> Option<f64> {
    let eps = 1.0e-12;
    let e1 = sub(tri.v[1], tri.v[0]);
    let e2 = sub(tri.v[2], tri.v[0]);
    let h = cross(dir, e2);
    let a = dot(e1, h);
    if a.abs() < eps {
        return None;
    }
    let f = 1.0 / a;
    let s = sub(origin, tri.v[0]);
    let u = f * dot(s, h);
    if !(0.0..=1.0).contains(&u) {
        return None;
    }
    let q = cross(s, e1);
    let v = f * dot(dir, q);
    if v < 0.0 || u + v > 1.0 {
        return None;
    }
    let t = f * dot(e2, q);
    (t > eps).then_some(t)
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::solver::UnsupportedReason;
    use std::io::Write;

    #[test]
    fn imports_small_cube_stl_fixture() {
        let dir = std::env::temp_dir().join(format!("lbm-cube-{}", std::process::id()));
        let _ = fs::remove_dir_all(&dir);
        fs::create_dir_all(&dir).unwrap();
        let stl = dir.join("cube.stl");
        let labels = dir.join("cube.json");
        write_unit_cube_binary_stl(&stl);
        fs::write(
            &labels,
            r#"{ "patches": [ { "name": "wall", "faces": [0,1,2,3,4,5,6,7,8,9,10,11] } ] }"#,
        )
        .unwrap();
        let imported = import_binary_stl_with_labels(
            &stl,
            &labels,
            VoxelImportOptions {
                dims: [20, 20, 20],
                dx_m: 0.05,
                credibility_tier: CredibilityTier::Evidence,
            },
        )
        .unwrap();
        assert!(imported.solid.iter().any(|&v| v));
        let _ = fs::remove_dir_all(&dir);
    }

    #[test]
    fn voxel_volume_within_5pct_of_analytic_cube() {
        let dir = std::env::temp_dir().join(format!("lbm-cube-volume-{}", std::process::id()));
        let _ = fs::remove_dir_all(&dir);
        fs::create_dir_all(&dir).unwrap();
        let stl = dir.join("cube.stl");
        let labels = dir.join("cube.json");
        write_unit_cube_binary_stl(&stl);
        fs::write(
            &labels,
            r#"{ "patches": [ { "name": "wall", "faces": [0,1,2,3,4,5,6,7,8,9,10,11] } ] }"#,
        )
        .unwrap();
        let dx = 0.025;
        let imported = import_binary_stl_with_labels(
            &stl,
            &labels,
            VoxelImportOptions {
                dims: [40, 40, 40],
                dx_m: dx,
                credibility_tier: CredibilityTier::Evidence,
            },
        )
        .unwrap();
        let volume = imported.solid.iter().filter(|&&v| v).count() as f64 * dx.powi(3);
        let rel = (volume - 1.0).abs();
        assert!(rel <= 0.05, "cube voxel volume rel error {rel:.4}");
        let _ = fs::remove_dir_all(&dir);
    }

    #[test]
    fn rejects_evidence_tier_with_unlabelled_faces() {
        let dir = std::env::temp_dir().join(format!("lbm-cube-evidence-{}", std::process::id()));
        let _ = fs::remove_dir_all(&dir);
        fs::create_dir_all(&dir).unwrap();
        let stl = dir.join("cube.stl");
        let labels = dir.join("cube.json");
        write_unit_cube_binary_stl(&stl);
        fs::write(
            &labels,
            r#"{ "patches": [ { "name": "wall", "faces": [0,1] } ] }"#,
        )
        .unwrap();
        let err = import_binary_stl_with_labels(
            &stl,
            &labels,
            VoxelImportOptions {
                dims: [10, 10, 10],
                dx_m: 0.1,
                credibility_tier: CredibilityTier::Evidence,
            },
        )
        .unwrap_err();
        assert!(matches!(
            err.reason,
            UnsupportedReason::EvidenceGateFailed { .. }
        ));
        let _ = fs::remove_dir_all(&dir);
    }

    fn write_unit_cube_binary_stl(path: &Path) {
        let tris = [
            ([0.0, 0.0, 0.0], [1.0, 1.0, 0.0], [1.0, 0.0, 0.0]),
            ([0.0, 0.0, 0.0], [0.0, 1.0, 0.0], [1.0, 1.0, 0.0]),
            ([0.0, 0.0, 1.0], [1.0, 0.0, 1.0], [1.0, 1.0, 1.0]),
            ([0.0, 0.0, 1.0], [1.0, 1.0, 1.0], [0.0, 1.0, 1.0]),
            ([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 0.0, 1.0]),
            ([0.0, 0.0, 0.0], [1.0, 0.0, 1.0], [0.0, 0.0, 1.0]),
            ([0.0, 1.0, 0.0], [1.0, 1.0, 1.0], [1.0, 1.0, 0.0]),
            ([0.0, 1.0, 0.0], [0.0, 1.0, 1.0], [1.0, 1.0, 1.0]),
            ([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 1.0, 1.0]),
            ([0.0, 0.0, 0.0], [0.0, 1.0, 1.0], [0.0, 1.0, 0.0]),
            ([1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [1.0, 1.0, 1.0]),
            ([1.0, 0.0, 0.0], [1.0, 1.0, 1.0], [1.0, 0.0, 1.0]),
        ];
        let mut bytes = vec![0u8; 80];
        bytes.extend_from_slice(&(tris.len() as u32).to_le_bytes());
        for tri in tris {
            bytes.extend_from_slice(&[0u8; 12]);
            for vertex in [tri.0, tri.1, tri.2] {
                for component in vertex {
                    bytes.extend_from_slice(&(component as f32).to_le_bytes());
                }
            }
            bytes.extend_from_slice(&0u16.to_le_bytes());
        }
        let mut file = fs::File::create(path).unwrap();
        file.write_all(&bytes).unwrap();
    }
}
