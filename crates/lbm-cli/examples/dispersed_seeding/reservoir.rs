use crate::particles::Particle;
use crate::protocol::ProtocolInput;

#[derive(Clone, Debug)]
pub struct Extraction {
    pub batch: Vec<Particle>,
    pub histogram: Vec<(f64, usize)>,
}

pub fn extract_by_depth(input: &ProtocolInput, particles: &[Particle]) -> Extraction {
    let withdraw = input.op("withdraw").expect("withdraw operation");
    let depth = withdraw.depth_frac.expect("withdraw.depth_frac validated");
    let volume_frac = withdraw
        .volume_frac
        .expect("withdraw.volume_frac validated");
    let settle_s = input
        .protocol
        .iter()
        .take_while(|op| op.op != "withdraw")
        .filter(|op| op.op == "settle")
        .filter_map(|op| op.duration_s)
        .sum::<f64>();
    let z0 = (1.0 - depth) * input.reservoir.fill_height_m;
    // Frozen §3.2 defines withdrawal as sampling a 1D settling column at the
    // requested depth. T18.3 validates the same Stokes/SN settling structure in
    // the core particle model; in a uniform column, concentration at z0 for a
    // diameter d is nonzero exactly when the back-traced origin
    // z0 + v_s(d) * t lies inside the initially filled column.
    let eligible = particles
        .iter()
        .enumerate()
        .filter_map(|(i, p)| concentration_at_depth(input, p.d_m, z0, settle_s).then_some(i))
        .collect::<Vec<_>>();
    let n = ((particles.len() as f64) * volume_frac).round() as usize;
    let mut batch = Vec::with_capacity(n.min(eligible.len()));
    for i in eligible.into_iter().take(n) {
        batch.push(particles[i].clone());
    }
    let histogram = diameter_histogram(&batch, input.particles.d_p_m, input.particles.d_p_cv);
    Extraction { batch, histogram }
}

fn concentration_at_depth(input: &ProtocolInput, d_m: f64, z0: f64, settle_s: f64) -> bool {
    let mu = input.fluid.rho_f_kgm3 * input.fluid.nu_m2s;
    let v =
        (input.particles.rho_p_kgm3 - input.fluid.rho_f_kgm3) * 9.80665 * d_m * d_m / (18.0 * mu);
    let origin = z0 + v * settle_s;
    (0.0..=input.reservoir.fill_height_m).contains(&origin)
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
