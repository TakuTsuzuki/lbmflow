use crate::particles::Particle;
use crate::protocol::ProtocolInput;

#[derive(Clone, Debug)]
pub struct Extraction {
    pub batch: Vec<Particle>,
    pub histogram: Vec<(f64, usize)>,
}

pub fn extract_by_depth(input: &ProtocolInput, particles: &[Particle]) -> Extraction {
    let withdraw = input.op("withdraw").expect("withdraw operation");
    let depth = withdraw.depth_frac.unwrap_or(0.5).clamp(0.0, 1.0);
    let volume_frac = withdraw.volume_frac.unwrap_or(1.0).clamp(0.0, 1.0);
    let rate = withdraw.rate_ul_s.unwrap_or(0.0);
    let settle_s = input
        .protocol
        .iter()
        .take_while(|op| op.op != "withdraw")
        .filter(|op| op.op == "settle")
        .filter_map(|op| op.duration_s)
        .sum::<f64>();
    let z0 = (1.0 - depth) * input.reservoir.fill_height_m;
    let band = (0.08 + 0.10 * (rate / 2000.0).clamp(0.0, 1.5)) * input.reservoir.fill_height_m;
    let mut scored: Vec<(f64, usize)> = particles
        .iter()
        .enumerate()
        .map(|(i, p)| {
            let settled_z = settled_column_z(input, p.d_m, p.pos[2], settle_s);
            let dz = ((settled_z - z0) / band.max(1.0e-9)).abs();
            let size = p.d_m / input.particles.d_p_m;
            let settled_bonus = if depth > 0.5 {
                size
            } else {
                1.0 / size.max(0.2)
            };
            let score = dz - 0.18 * settled_bonus;
            (score, i)
        })
        .collect();
    scored.sort_by(|a, b| a.0.total_cmp(&b.0));
    let n = ((particles.len() as f64) * volume_frac).round() as usize;
    let mut batch = Vec::with_capacity(n.min(scored.len()));
    for (_, i) in scored.into_iter().take(n) {
        batch.push(particles[i].clone());
    }
    let histogram = diameter_histogram(&batch, input.particles.d_p_m, input.particles.d_p_cv);
    Extraction { batch, histogram }
}

fn settled_column_z(input: &ProtocolInput, d_m: f64, z_initial: f64, settle_s: f64) -> f64 {
    let mu = input.fluid.rho_f_kgm3 * input.fluid.nu_m2s;
    let v =
        (input.particles.rho_p_kgm3 - input.fluid.rho_f_kgm3) * 9.80665 * d_m * d_m / (18.0 * mu);
    (z_initial - v * settle_s).clamp(0.0, input.reservoir.fill_height_m)
}

fn diameter_histogram(batch: &[Particle], mean: f64, cv: f64) -> Vec<(f64, usize)> {
    let min = mean * (1.0 - 3.0 * cv).max(0.1);
    let max = mean * (1.0 + 3.0 * cv).max(1.1);
    let bins = 8usize;
    let mut counts = vec![0usize; bins];
    for p in batch {
        let t = ((p.d_m - min) / (max - min)).clamp(0.0, 0.999_999);
        counts[(t * bins as f64) as usize] += 1;
    }
    (0..bins)
        .map(|i| {
            let center = min + (i as f64 + 0.5) * (max - min) / bins as f64;
            (center, counts[i])
        })
        .collect()
}
