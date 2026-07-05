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

/// 16-anchor inferno (sequential, for magnitude fields like shear / dissipation).
const INFERNO: [[u8; 3]; 16] = [
    [0, 0, 4],
    [12, 8, 38],
    [36, 12, 79],
    [66, 10, 104],
    [93, 18, 110],
    [120, 28, 109],
    [147, 38, 103],
    [174, 48, 92],
    [199, 62, 76],
    [220, 81, 57],
    [237, 105, 37],
    [246, 133, 17],
    [251, 163, 12],
    [249, 195, 41],
    [240, 226, 96],
    [252, 255, 164],
];

/// Colormap selector for the shared PNG writer.
#[derive(Clone, Copy, Debug)]
pub enum Colormap {
    /// Sequential blue→yellow (`0..max`).
    Viridis,
    /// Sequential black→yellow (`0..max`) — magnitude fields (shear, speed).
    Inferno,
    /// Diverging blue↔red on a symmetric range — signed fields.
    RdBu,
}

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

/// Shared scalar-field PNG writer — the one raster surface for the runner and
/// examples (so examples stop re-implementing it). `cmap` picks the colormap;
/// `vmax` fixes the range (`Some(v)` → `[0,v]` sequential or `[-v,v]` diverging;
/// `None` → auto-range over finite non-solid cells); `scale` supersamples each
/// cell into a `scale × scale` block (`1` = native). Solids render grey; the
/// image is flipped vertically (y up) for PNG.
#[allow(clippy::too_many_arguments)]
pub fn write_png_scaled(
    path: &Path,
    field: &[f64],
    solid: &[bool],
    nx: usize,
    ny: usize,
    cmap: Colormap,
    vmax: Option<f64>,
    scale: usize,
) -> Result<()> {
    let anchors: &[[u8; 3]] = match cmap {
        Colormap::Viridis => &VIRIDIS,
        Colormap::Inferno => &INFERNO,
        Colormap::RdBu => &RDBU,
    };
    let diverging = matches!(cmap, Colormap::RdBu);
    let (lo, hi) = match vmax {
        Some(v) if diverging => (-v.abs(), v.abs()),
        Some(v) => (0.0, v.abs()),
        None => {
            let mut lo = f64::INFINITY;
            let mut hi = f64::NEG_INFINITY;
            for (v, s) in field.iter().zip(solid) {
                if !s && v.is_finite() {
                    lo = lo.min(*v);
                    hi = hi.max(*v);
                }
            }
            if !lo.is_finite() || !hi.is_finite() {
                (0.0, 1.0)
            } else {
                (lo, hi)
            }
        }
    };
    let sc = scale.max(1);
    let (ow, oh) = (nx * sc, ny * sc);
    let mut buf = vec![0u8; ow * oh * 3];
    for oy in 0..oh {
        let y = oy / sc;
        for ox in 0..ow {
            let x = ox / sc;
            let i = y * nx + x;
            let px = ((oh - 1 - oy) * ow + ox) * 3; // flip vertically for PNG
            let rgb = if solid[i] {
                [90u8, 94, 100]
            } else if diverging {
                let m = lo.abs().max(hi.abs()).max(1e-30);
                lut(anchors, 0.5 + 0.5 * field[i] / m)
            } else {
                let span = (hi - lo).max(1e-30);
                lut(anchors, (field[i] - lo) / span)
            };
            buf[px..px + 3].copy_from_slice(&rgb);
        }
    }
    let file = File::create(path)?;
    let mut enc = png::Encoder::new(BufWriter::new(file), ow as u32, oh as u32);
    enc.set_color(png::ColorType::Rgb);
    enc.set_depth(png::BitDepth::Eight);
    let mut writer = enc.write_header()?;
    writer.write_image_data(&buf)?;
    Ok(())
}
