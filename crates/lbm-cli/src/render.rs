//! Field → PNG rendering with compact colormap LUTs.

use anyhow::Result;
use std::fs::File;
use std::io::BufWriter;
use std::path::Path;

/// 16-anchor viridis; interpolated linearly.
const VIRIDIS: [[u8; 3]; 16] = [
    [68, 1, 84],
    [72, 26, 108],
    [71, 47, 125],
    [65, 68, 135],
    [57, 86, 140],
    [49, 104, 142],
    [42, 120, 142],
    [35, 136, 142],
    [31, 152, 139],
    [34, 168, 132],
    [53, 183, 121],
    [84, 197, 104],
    [122, 209, 81],
    [165, 219, 54],
    [210, 226, 27],
    [253, 231, 37],
];

/// 9-anchor RdBu (diverging, blue = negative).
const RDBU: [[u8; 3]; 9] = [
    [5, 48, 97],
    [33, 102, 172],
    [67, 147, 195],
    [146, 197, 222],
    [247, 247, 247],
    [244, 165, 130],
    [214, 96, 77],
    [178, 24, 43],
    [103, 0, 31],
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

/// Write a scalar field (row-major, y=0 at the bottom) as PNG.
/// `diverging` uses RdBu with symmetric range; otherwise viridis 0..max.
pub fn write_png(
    path: &Path,
    field: &[f64],
    solid: &[bool],
    nx: usize,
    ny: usize,
    diverging: bool,
) -> Result<()> {
    let mut lo = f64::INFINITY;
    let mut hi = f64::NEG_INFINITY;
    for (v, s) in field.iter().zip(solid) {
        if !s && v.is_finite() {
            lo = lo.min(*v);
            hi = hi.max(*v);
        }
    }
    if !lo.is_finite() || !hi.is_finite() {
        lo = 0.0;
        hi = 1.0;
    }
    let mut buf = vec![0u8; nx * ny * 3];
    for y in 0..ny {
        for x in 0..nx {
            let i = y * nx + x;
            let px = ((ny - 1 - y) * nx + x) * 3; // flip vertically for PNG
            let rgb = if solid[i] {
                [90u8, 94, 100]
            } else if diverging {
                let m = lo.abs().max(hi.abs()).max(1e-30);
                lut(&RDBU, 0.5 + 0.5 * field[i] / m)
            } else {
                let span = (hi - lo).max(1e-30);
                lut(&VIRIDIS, (field[i] - lo) / span)
            };
            buf[px..px + 3].copy_from_slice(&rgb);
        }
    }
    let file = File::create(path)?;
    let mut enc = png::Encoder::new(BufWriter::new(file), nx as u32, ny as u32);
    enc.set_color(png::ColorType::Rgb);
    enc.set_depth(png::BitDepth::Eight);
    let mut writer = enc.write_header()?;
    writer.write_image_data(&buf)?;
    Ok(())
}
