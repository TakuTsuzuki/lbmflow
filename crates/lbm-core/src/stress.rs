//! Velocity-gradient stress and wall-shear diagnostics for bioprocess QOIs.

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct StressCell {
    pub s: [f64; 6],
    pub gamma_dot: f64,
    pub viscous_stress_pa: f64,
    pub second_invariant_s2: f64,
    pub von_mises_proxy_pa: f64,
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct PercentileSummary {
    pub p50: f64,
    pub p90: f64,
    pub p95: f64,
    pub p99: f64,
    pub max: f64,
    pub fraction_above_threshold: Option<f64>,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WallShearProxy {
    pub cell_index: usize,
    pub tau_w_pa: f64,
    pub y_plus: Option<f64>,
    pub u_parallel_m_s: f64,
    pub y_m: f64,
    pub normal: [f64; 3],
}

pub fn compute_stress_field(
    dims: [usize; 3],
    ux: &[f64],
    uy: &[f64],
    uz: &[f64],
    solid: &[bool],
    dx_m: f64,
    mu_eff_pa_s: &[f64],
) -> Vec<StressCell> {
    let n = dims[0] * dims[1] * dims[2];
    assert_eq!(ux.len(), n);
    assert_eq!(uy.len(), n);
    assert_eq!(uz.len(), n);
    assert_eq!(solid.len(), n);
    assert_eq!(mu_eff_pa_s.len(), n);
    assert!(dx_m.is_finite() && dx_m > 0.0);
    let mut out = vec![StressCell::default(); n];
    for z in 0..dims[2] {
        for y in 0..dims[1] {
            for x in 0..dims[0] {
                let i = idx(dims, x, y, z);
                if solid[i] {
                    continue;
                }
                let grad = velocity_gradient(dims, ux, uy, uz, solid, dx_m, x, y, z);
                let s = [
                    grad[0][0],
                    grad[1][1],
                    grad[2][2],
                    0.5 * (grad[0][1] + grad[1][0]),
                    0.5 * (grad[0][2] + grad[2][0]),
                    0.5 * (grad[1][2] + grad[2][1]),
                ];
                let s_colon_s = s[0] * s[0]
                    + s[1] * s[1]
                    + s[2] * s[2]
                    + 2.0 * (s[3] * s[3] + s[4] * s[4] + s[5] * s[5]);
                let gamma_dot = (2.0 * s_colon_s).sqrt();
                let mu = mu_eff_pa_s[i];
                out[i] = StressCell {
                    s,
                    gamma_dot,
                    viscous_stress_pa: mu * gamma_dot,
                    second_invariant_s2: s_colon_s,
                    von_mises_proxy_pa: 3.0_f64.sqrt() * mu * gamma_dot,
                };
            }
        }
    }
    out
}

pub fn percentile_summary(values: &[f64], threshold: Option<f64>) -> Option<PercentileSummary> {
    let mut finite: Vec<f64> = values.iter().copied().filter(|v| v.is_finite()).collect();
    if finite.is_empty() {
        return None;
    }
    finite.sort_by(|a, b| a.total_cmp(b));
    let fraction_above_threshold =
        threshold.map(|t| finite.iter().filter(|&&v| v > t).count() as f64 / finite.len() as f64);
    Some(PercentileSummary {
        p50: percentile_sorted(&finite, 0.50),
        p90: percentile_sorted(&finite, 0.90),
        p95: percentile_sorted(&finite, 0.95),
        p99: percentile_sorted(&finite, 0.99),
        max: *finite.last().expect("non-empty"),
        fraction_above_threshold,
    })
}

pub fn wall_shear_proxy(
    dims: [usize; 3],
    ux: &[f64],
    uy: &[f64],
    uz: &[f64],
    solid: &[bool],
    wall_u: &[[f64; 3]],
    dx_m: f64,
    rho_kg_m3: f64,
    mu_pa_s: f64,
) -> Vec<WallShearProxy> {
    let n = dims[0] * dims[1] * dims[2];
    assert_eq!(ux.len(), n);
    assert_eq!(uy.len(), n);
    assert_eq!(uz.len(), n);
    assert_eq!(solid.len(), n);
    assert_eq!(wall_u.len(), n);
    assert!(dx_m.is_finite() && dx_m > 0.0);
    assert!(rho_kg_m3.is_finite() && rho_kg_m3 > 0.0);
    assert!(mu_pa_s.is_finite() && mu_pa_s > 0.0);
    let y_m = 0.5 * dx_m;
    let nu = mu_pa_s / rho_kg_m3;
    let mut out = Vec::new();
    for z in 0..dims[2] {
        for y in 0..dims[1] {
            for x in 0..dims[0] {
                let i = idx(dims, x, y, z);
                if solid[i] {
                    continue;
                }
                let Some((normal, wu)) =
                    adjacent_wall_normal_and_velocity(dims, solid, wall_u, x, y, z)
                else {
                    continue;
                };
                let u = [ux[i] - wu[0], uy[i] - wu[1], uz[i] - wu[2]];
                let un = dot(u, normal);
                let tan = [
                    u[0] - un * normal[0],
                    u[1] - un * normal[1],
                    u[2] - un * normal[2],
                ];
                let u_parallel = dot(tan, tan).sqrt();
                let tau = mu_pa_s * u_parallel / y_m;
                let u_tau = (tau / rho_kg_m3).sqrt();
                out.push(WallShearProxy {
                    cell_index: i,
                    tau_w_pa: tau,
                    y_plus: Some(y_m * u_tau / nu),
                    u_parallel_m_s: u_parallel,
                    y_m,
                    normal,
                });
            }
        }
    }
    out.sort_by_key(|v| v.cell_index);
    out
}

fn velocity_gradient(
    dims: [usize; 3],
    ux: &[f64],
    uy: &[f64],
    uz: &[f64],
    solid: &[bool],
    dx: f64,
    x: usize,
    y: usize,
    z: usize,
) -> [[f64; 3]; 3] {
    let fields = [ux, uy, uz];
    let pos = [x, y, z];
    let mut out = [[0.0; 3]; 3];
    for comp in 0..3 {
        for axis in 0..3 {
            if dims[axis] <= 1 {
                continue;
            }
            let plus = neighbor_value(dims, fields[comp], solid, pos, axis, 1);
            let minus = neighbor_value(dims, fields[comp], solid, pos, axis, -1);
            let own = fields[comp][idx(dims, x, y, z)];
            out[comp][axis] = match (plus, minus) {
                (Some(p), Some(m)) => (p - m) / (2.0 * dx),
                (Some(p), None) => (p - own) / dx,
                (None, Some(m)) => (own - m) / dx,
                (None, None) => 0.0,
            };
        }
    }
    out
}

fn neighbor_value(
    dims: [usize; 3],
    f: &[f64],
    solid: &[bool],
    mut pos: [usize; 3],
    axis: usize,
    delta: isize,
) -> Option<f64> {
    let p = pos[axis] as isize + delta;
    if p < 0 || p >= dims[axis] as isize {
        return None;
    }
    pos[axis] = p as usize;
    let i = idx(dims, pos[0], pos[1], pos[2]);
    (!solid[i]).then_some(f[i])
}

fn adjacent_wall_normal_and_velocity(
    dims: [usize; 3],
    solid: &[bool],
    wall_u: &[[f64; 3]],
    x: usize,
    y: usize,
    z: usize,
) -> Option<([f64; 3], [f64; 3])> {
    let mut normal = [0.0; 3];
    let mut wu = [0.0; 3];
    let mut count = 0.0;
    for dz in -1isize..=1 {
        for dy in -1isize..=1 {
            for dx in -1isize..=1 {
                if dx == 0 && dy == 0 && dz == 0 {
                    continue;
                }
                let nx = x as isize + dx;
                let ny = y as isize + dy;
                let nz = z as isize + dz;
                if nx < 0
                    || ny < 0
                    || nz < 0
                    || nx >= dims[0] as isize
                    || ny >= dims[1] as isize
                    || nz >= dims[2] as isize
                {
                    continue;
                }
                let j = idx(dims, nx as usize, ny as usize, nz as usize);
                if !solid[j] {
                    continue;
                }
                normal[0] += dx as f64;
                normal[1] += dy as f64;
                normal[2] += dz as f64;
                for a in 0..3 {
                    wu[a] += wall_u[j][a];
                }
                count += 1.0;
            }
        }
    }
    if count == 0.0 {
        return None;
    }
    let mag = dot(normal, normal).sqrt();
    if mag == 0.0 {
        return None;
    }
    for a in 0..3 {
        normal[a] /= mag;
        wu[a] /= count;
    }
    Some((normal, wu))
}

fn percentile_sorted(sorted: &[f64], p: f64) -> f64 {
    debug_assert!(!sorted.is_empty());
    if sorted.len() == 1 {
        return sorted[0];
    }
    let rank = p * (sorted.len() - 1) as f64;
    let lo = rank.floor() as usize;
    let hi = rank.ceil() as usize;
    if lo == hi {
        sorted[lo]
    } else {
        let t = rank - lo as f64;
        sorted[lo] * (1.0 - t) + sorted[hi] * t
    }
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn idx(dims: [usize; 3], x: usize, y: usize, z: usize) -> usize {
    (z * dims[1] + y) * dims[0] + x
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn couette_gamma_dot_matches_analytic_and_has_behavior_anchor() {
        let dims = [8, 10, 1];
        let n = dims[0] * dims[1];
        let mut ux = vec![0.0; n];
        let uy = vec![0.0; n];
        let uz = vec![0.0; n];
        let solid = vec![false; n];
        let h = (dims[1] - 1) as f64;
        for y in 0..dims[1] {
            for x in 0..dims[0] {
                ux[idx(dims, x, y, 0)] = y as f64 / h;
            }
        }
        let mu = vec![2.0; n];
        let stress = compute_stress_field(dims, &ux, &uy, &uz, &solid, 1.0, &mu);
        for y in 1..dims[1] - 1 {
            for x in 1..dims[0] - 1 {
                let g = stress[idx(dims, x, y, 0)].gamma_dot;
                assert!((g - 1.0 / h).abs() < 1.0e-12, "gamma_dot={g}");
            }
        }
        assert!(
            stress[idx(dims, 4, 8, 0)].viscous_stress_pa > stress[idx(dims, 4, 8, 0)].gamma_dot
        );
    }

    #[test]
    fn poiseuille_shear_profile_is_antisymmetric() {
        let dims = [8, 17, 1];
        let n = dims[0] * dims[1];
        let mut ux = vec![0.0; n];
        let uy = vec![0.0; n];
        let uz = vec![0.0; n];
        let solid = vec![false; n];
        let c = (dims[1] - 1) as f64 / 2.0;
        for y in 0..dims[1] {
            let yy = y as f64 - c;
            for x in 0..dims[0] {
                ux[idx(dims, x, y, 0)] = 1.0 - yy * yy / (c * c);
            }
        }
        let mu = vec![1.0; n];
        let stress = compute_stress_field(dims, &ux, &uy, &uz, &solid, 1.0, &mu);
        for off in 1..7 {
            let lo = stress[idx(dims, 4, (c as usize) - off, 0)].s[3];
            let hi = stress[idx(dims, 4, (c as usize) + off, 0)].s[3];
            assert!((lo + hi).abs() < 1.0e-12, "Sxy not antisymmetric");
            assert!(stress[idx(dims, 4, (c as usize) + off, 0)].gamma_dot > 0.0);
        }
    }

    #[test]
    fn percentile_reducer_matches_synthetic_distribution() {
        let values: Vec<f64> = (0..=100).map(|v| v as f64).collect();
        let s = percentile_summary(&values, Some(90.0)).unwrap();
        assert_eq!(s.p50, 50.0);
        assert_eq!(s.p90, 90.0);
        assert_eq!(s.p95, 95.0);
        assert_eq!(s.p99, 99.0);
        assert_eq!(s.max, 100.0);
        assert!((s.fraction_above_threshold.unwrap() - 10.0 / 101.0).abs() < 1.0e-12);
    }

    #[test]
    fn couette_wall_shear_proxy_matches_one_sided_analytic() {
        let dims = [8, 8, 1];
        let n = dims[0] * dims[1];
        let mut ux = vec![0.0; n];
        let uy = vec![0.0; n];
        let uz = vec![0.0; n];
        let mut solid = vec![false; n];
        let mut wall_u = vec![[0.0; 3]; n];
        for x in 0..dims[0] {
            solid[idx(dims, x, 0, 0)] = true;
            solid[idx(dims, x, dims[1] - 1, 0)] = true;
            wall_u[idx(dims, x, dims[1] - 1, 0)] = [1.0, 0.0, 0.0];
        }
        for y in 1..dims[1] - 1 {
            for x in 0..dims[0] {
                ux[idx(dims, x, y, 0)] = (y as f64 - 0.5) / (dims[1] as f64 - 2.0);
            }
        }
        let wall = wall_shear_proxy(dims, &ux, &uy, &uz, &solid, &wall_u, 1.0, 1000.0, 0.001);
        let bottom = wall
            .iter()
            .find(|m| m.cell_index == idx(dims, 4, 1, 0))
            .unwrap();
        assert!((bottom.tau_w_pa - 0.001 / 6.0).abs() < 1.0e-12);
        let top = wall
            .iter()
            .find(|m| m.cell_index == idx(dims, 4, 6, 0))
            .unwrap();
        assert!(
            top.tau_w_pa > 0.0,
            "moving-wall sign should report positive shear magnitude"
        );
        assert!(bottom.y_plus.unwrap().is_finite());
    }
}
