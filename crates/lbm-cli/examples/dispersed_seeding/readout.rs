use crate::protocol::{ProtocolInput, Regime};
use lbm_core::particles::DepositEvent;
use serde::Serialize;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

#[derive(Clone, Debug, Serialize)]
pub struct Metrics {
    #[serde(rename = "CV")]
    pub cv: f64,
    pub max_over_mean: f64,
    pub empty_bin_fraction: f64,
    pub n_deposited: usize,
    pub n_suspended: usize,
    pub n_extracted: usize,
    #[serde(rename = "Re_jet")]
    pub re_jet: f64,
    #[serde(rename = "St")]
    pub st: f64,
    #[serde(rename = "Fr")]
    pub fr: f64,
    #[serde(rename = "Ma")]
    pub ma: f64,
    pub tau: f64,
}

pub fn write_outputs(
    input: &ProtocolInput,
    regime: &Regime,
    deposits: &[DepositEvent],
    n_suspended: usize,
    n_extracted: usize,
    reservoir_velocity: &[[f64; 3]],
    tray_velocity: &[[f64; 3]],
) -> anyhow::Result<Metrics> {
    let outdir = input.output_dir();
    std::fs::create_dir_all(&outdir)?;
    let bins = bin_counts(input, regime, deposits);
    if input.output.csv {
        write_density_csv(&outdir.join("density.csv"), input, &bins)?;
    }
    if input.output.volume {
        write_vtk(
            &outdir.join("reservoir_velocity.vtk"),
            [input.grid.res_nx, input.grid.res_ny, input.grid.res_nz],
            input.grid.dx_m,
            reservoir_velocity,
        )?;
        write_vtk(
            &outdir.join("tray_velocity.vtk"),
            [input.grid.tray_nx, input.grid.tray_ny, input.grid.tray_nz],
            regime.dx,
            tray_velocity,
        )?;
    }
    let n = bins.len() as f64;
    let sum: usize = bins.iter().sum();
    let mean = sum as f64 / n.max(1.0);
    let var = bins
        .iter()
        .map(|&c| {
            let d = c as f64 - mean;
            d * d
        })
        .sum::<f64>()
        / n.max(1.0);
    let cv = if mean > 0.0 { var.sqrt() / mean } else { 0.0 };
    let max = bins.iter().copied().max().unwrap_or(0) as f64;
    let empty_bin_fraction =
        bins.iter().filter(|&&c| c == 0).count() as f64 / bins.len().max(1) as f64;
    let metrics = Metrics {
        cv,
        max_over_mean: if mean > 0.0 { max / mean } else { 0.0 },
        empty_bin_fraction,
        n_deposited: sum,
        n_suspended,
        n_extracted,
        re_jet: regime.re_jet,
        st: regime.st,
        fr: regime.fr,
        ma: regime.ma,
        tau: regime.tau,
    };
    let file = File::create(outdir.join("metrics.json"))?;
    serde_json::to_writer_pretty(BufWriter::new(file), &metrics)?;
    Ok(metrics)
}

fn bin_counts(input: &ProtocolInput, regime: &Regime, deposits: &[DepositEvent]) -> Vec<usize> {
    let mx = input.target.partitions_x;
    let my = input.target.partitions_y;
    let mut bins = vec![0usize; mx * my];
    for event in deposits {
        let x_m = event.pos[0] * regime.dx;
        let y_m = event.pos[1] * regime.dx;
        let i = ((x_m / input.target.width_m) * mx as f64)
            .floor()
            .clamp(0.0, (mx - 1) as f64) as usize;
        let j = ((y_m / input.target.depth_m) * my as f64)
            .floor()
            .clamp(0.0, (my - 1) as f64) as usize;
        bins[j * mx + i] += 1;
    }
    bins
}

fn write_density_csv(path: &Path, input: &ProtocolInput, bins: &[usize]) -> anyhow::Result<()> {
    let mut w = BufWriter::new(File::create(path)?);
    writeln!(
        w,
        "bin_i,bin_j,x_center_m,y_center_m,count,normalized_density"
    )?;
    let mx = input.target.partitions_x;
    let my = input.target.partitions_y;
    let mean = bins.iter().sum::<usize>() as f64 / bins.len().max(1) as f64;
    for j in 0..my {
        for i in 0..mx {
            let c = bins[j * mx + i];
            let x = (i as f64 + 0.5) * input.target.width_m / mx as f64;
            let y = (j as f64 + 0.5) * input.target.depth_m / my as f64;
            let nd = if mean > 0.0 { c as f64 / mean } else { 0.0 };
            writeln!(w, "{i},{j},{x:.8},{y:.8},{c},{nd:.8}")?;
        }
    }
    Ok(())
}

fn write_vtk(path: &Path, dims: [usize; 3], dx: f64, velocity: &[[f64; 3]]) -> anyhow::Result<()> {
    let mut w = BufWriter::new(File::create(path)?);
    writeln!(w, "# vtk DataFile Version 3.0")?;
    writeln!(w, "LBMFlow dispersed seeding velocity")?;
    writeln!(w, "ASCII")?;
    writeln!(w, "DATASET STRUCTURED_POINTS")?;
    writeln!(w, "DIMENSIONS {} {} {}", dims[0], dims[1], dims[2])?;
    writeln!(w, "ORIGIN 0 0 0")?;
    writeln!(w, "SPACING {dx} {dx} {dx}")?;
    writeln!(w, "POINT_DATA {}", dims[0] * dims[1] * dims[2])?;
    writeln!(w, "VECTORS velocity float")?;
    for v in velocity {
        writeln!(w, "{} {} {}", v[0], v[1], v[2])?;
    }
    Ok(())
}
