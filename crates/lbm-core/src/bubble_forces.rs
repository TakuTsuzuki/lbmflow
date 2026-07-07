//! Literature-backed point-bubble force closures with hard validity guards.

use crate::bubbles::{bubble_volume_from_diameter, Bubble, BubbleError, BubbleResult};

pub const SCHILLER_NAUMANN_RE_MAX: f64 = 800.0;
pub const POINT_BUBBLE_ALPHA_G_MAX: f64 = 0.3;
pub const PLACEHOLDER_KW_MIN: f64 = 0.0;
pub const PLACEHOLDER_KW_MAX: f64 = 1.0e6;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BubbleForceContext {
    pub re_bubble: f64,
    pub alpha_g: f64,
    pub kolmogorov_weber: f64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ClosureValidity {
    pub re_bubble_max: f64,
    pub alpha_g_max: f64,
    pub kolmogorov_weber_min: f64,
    pub kolmogorov_weber_max: f64,
}

impl ClosureValidity {
    pub fn point_bubble_default() -> Self {
        Self {
            re_bubble_max: SCHILLER_NAUMANN_RE_MAX,
            alpha_g_max: POINT_BUBBLE_ALPHA_G_MAX,
            kolmogorov_weber_min: PLACEHOLDER_KW_MIN,
            kolmogorov_weber_max: PLACEHOLDER_KW_MAX,
        }
    }

    pub fn validate(&self, ctx: BubbleForceContext, closure: &'static str) -> BubbleResult<()> {
        if !(ctx.re_bubble.is_finite() && ctx.re_bubble >= 0.0) {
            return Err(BubbleError::out_of_validity_range(
                format!("{closure} requires finite non-negative Re_bubble"),
                "Re_bubble must be finite and >= 0",
            ));
        }
        if ctx.re_bubble > self.re_bubble_max {
            return Err(BubbleError::out_of_validity_range(
                format!("{closure} is outside its Re_bubble validity range"),
                format!(
                    "Re_bubble={:.6} exceeds {:.6}",
                    ctx.re_bubble, self.re_bubble_max
                ),
            ));
        }
        if !(ctx.alpha_g.is_finite() && ctx.alpha_g >= 0.0 && ctx.alpha_g <= self.alpha_g_max) {
            return Err(BubbleError::out_of_validity_range(
                format!("{closure} requires dilute point-bubble holdup"),
                format!(
                    "alpha_g={:.6} must be in [0, {:.6}]",
                    ctx.alpha_g, self.alpha_g_max
                ),
            ));
        }
        if !(ctx.kolmogorov_weber.is_finite()
            && ctx.kolmogorov_weber >= self.kolmogorov_weber_min
            && ctx.kolmogorov_weber <= self.kolmogorov_weber_max)
        {
            return Err(BubbleError::out_of_validity_range(
                format!("{closure} is outside its kW validity range"),
                format!(
                    "kW={:.6} must be in [{:.6}, {:.6}]",
                    ctx.kolmogorov_weber, self.kolmogorov_weber_min, self.kolmogorov_weber_max
                ),
            ));
        }
        Ok(())
    }
}

pub fn bubble_reynolds(slip_velocity_m_s: [f64; 3], diameter_m: f64, nu_liquid_m2_s: f64) -> f64 {
    norm(slip_velocity_m_s) * diameter_m / nu_liquid_m2_s
}

pub fn buoyancy_force(
    rho_liquid_kg_m3: f64,
    gravity_m_s2: [f64; 3],
    bubble_volume_m3: f64,
    ctx: BubbleForceContext,
) -> BubbleResult<[f64; 3]> {
    ClosureValidity::point_bubble_default().validate(ctx, "buoyancy")?;
    if !(rho_liquid_kg_m3.is_finite() && rho_liquid_kg_m3 > 0.0 && bubble_volume_m3 > 0.0) {
        return Err(BubbleError::out_of_validity_range(
            "buoyancy requires positive liquid density and bubble volume",
            "rho_liquid_kg_m3 and bubble_volume_m3 must be > 0",
        ));
    }
    Ok([
        -rho_liquid_kg_m3 * bubble_volume_m3 * gravity_m_s2[0],
        -rho_liquid_kg_m3 * bubble_volume_m3 * gravity_m_s2[1],
        -rho_liquid_kg_m3 * bubble_volume_m3 * gravity_m_s2[2],
    ])
}

pub fn schiller_naumann_drag_force(
    rho_liquid_kg_m3: f64,
    nu_liquid_m2_s: f64,
    diameter_m: f64,
    slip_velocity_m_s: [f64; 3],
    ctx: BubbleForceContext,
) -> BubbleResult<[f64; 3]> {
    ClosureValidity::point_bubble_default().validate(ctx, "Schiller-Naumann drag")?;
    if !(rho_liquid_kg_m3.is_finite()
        && rho_liquid_kg_m3 > 0.0
        && nu_liquid_m2_s.is_finite()
        && nu_liquid_m2_s > 0.0
        && diameter_m.is_finite()
        && diameter_m > 0.0)
    {
        return Err(BubbleError::out_of_validity_range(
            "drag requires positive liquid density, viscosity, and bubble diameter",
            "rho_liquid_kg_m3, nu_liquid_m2_s and diameter_m must be > 0",
        ));
    }
    let speed = norm(slip_velocity_m_s);
    if speed == 0.0 {
        return Ok([0.0; 3]);
    }
    let re = bubble_reynolds(slip_velocity_m_s, diameter_m, nu_liquid_m2_s);
    if re > SCHILLER_NAUMANN_RE_MAX {
        return Err(BubbleError::out_of_validity_range(
            "Schiller-Naumann drag is outside its validity range",
            format!("Re_bubble={re:.6} exceeds {SCHILLER_NAUMANN_RE_MAX:.1}"),
        ));
    }
    let cd = if re == 0.0 {
        0.0
    } else {
        (24.0 / re) * (1.0 + 0.15 * re.powf(0.687))
    };
    let area = std::f64::consts::PI * diameter_m * diameter_m / 4.0;
    let mag = 0.5 * rho_liquid_kg_m3 * cd * area * speed * speed;
    Ok([
        mag * slip_velocity_m_s[0] / speed,
        mag * slip_velocity_m_s[1] / speed,
        mag * slip_velocity_m_s[2] / speed,
    ])
}

pub fn added_mass_force(
    rho_liquid_kg_m3: f64,
    bubble_volume_m3: f64,
    liquid_material_acceleration_m_s2: [f64; 3],
    bubble_acceleration_m_s2: [f64; 3],
    ctx: BubbleForceContext,
) -> BubbleResult<[f64; 3]> {
    ClosureValidity::point_bubble_default().validate(ctx, "added mass")?;
    if !(rho_liquid_kg_m3.is_finite() && rho_liquid_kg_m3 > 0.0 && bubble_volume_m3 > 0.0) {
        return Err(BubbleError::out_of_validity_range(
            "added mass requires positive liquid density and bubble volume",
            "rho_liquid_kg_m3 and bubble_volume_m3 must be > 0",
        ));
    }
    let c = 0.5 * rho_liquid_kg_m3 * bubble_volume_m3;
    Ok([
        c * (liquid_material_acceleration_m_s2[0] - bubble_acceleration_m_s2[0]),
        c * (liquid_material_acceleration_m_s2[1] - bubble_acceleration_m_s2[1]),
        c * (liquid_material_acceleration_m_s2[2] - bubble_acceleration_m_s2[2]),
    ])
}

pub fn lift_placeholder_force(
    rho_liquid_kg_m3: f64,
    bubble_volume_m3: f64,
    slip_velocity_m_s: [f64; 3],
    vorticity_1_s: [f64; 3],
    ctx: BubbleForceContext,
) -> BubbleResult<[f64; 3]> {
    ClosureValidity::point_bubble_default().validate(ctx, "constant-C_L lift placeholder")?;
    let coeff = 0.5 * rho_liquid_kg_m3 * bubble_volume_m3;
    Ok(scale(cross(slip_velocity_m_s, vorticity_1_s), coeff))
}

pub fn wall_lubrication_placeholder_force(
    rho_liquid_kg_m3: f64,
    bubble_volume_m3: f64,
    slip_velocity_m_s: [f64; 3],
    wall_normal_into_liquid: [f64; 3],
    wall_distance_m: f64,
    diameter_m: f64,
    ctx: BubbleForceContext,
) -> BubbleResult<[f64; 3]> {
    ClosureValidity::point_bubble_default()
        .validate(ctx, "Tomiyama-like wall lubrication placeholder")?;
    if !(wall_distance_m.is_finite() && wall_distance_m > 0.0 && diameter_m > 0.0) {
        return Err(BubbleError::out_of_validity_range(
            "wall lubrication requires positive wall distance and diameter",
            "wall_distance_m and diameter_m must be > 0",
        ));
    }
    if wall_distance_m > 2.0 * diameter_m {
        return Ok([0.0; 3]);
    }
    let slip2 = dot(slip_velocity_m_s, slip_velocity_m_s);
    let magnitude =
        0.5 * rho_liquid_kg_m3 * bubble_volume_m3 * slip2 / (wall_distance_m * wall_distance_m);
    Ok(scale(unit(wall_normal_into_liquid)?, magnitude))
}

pub fn turbulent_dispersion_placeholder_force(
    rho_liquid_kg_m3: f64,
    bubble_volume_m3: f64,
    alpha_g_gradient_1_m: [f64; 3],
    turbulent_kinetic_energy_m2_s2: f64,
    ctx: BubbleForceContext,
) -> BubbleResult<[f64; 3]> {
    ClosureValidity::point_bubble_default()
        .validate(ctx, "Lopez-de-Bertodano turbulent dispersion placeholder")?;
    if !(turbulent_kinetic_energy_m2_s2.is_finite() && turbulent_kinetic_energy_m2_s2 >= 0.0) {
        return Err(BubbleError::out_of_validity_range(
            "turbulent dispersion requires finite non-negative k",
            "turbulent_kinetic_energy_m2_s2 must be >= 0",
        ));
    }
    let coeff = -0.1 * rho_liquid_kg_m3 * bubble_volume_m3 * turbulent_kinetic_energy_m2_s2;
    Ok(scale(alpha_g_gradient_1_m, coeff))
}

pub fn rk4_substep<F>(
    bubble: &Bubble,
    dt_s: f64,
    mut acceleration: F,
) -> BubbleResult<([f64; 3], [f64; 3])>
where
    F: FnMut([f64; 3], [f64; 3]) -> BubbleResult<[f64; 3]>,
{
    if !(dt_s.is_finite() && dt_s >= 0.0) {
        return Err(BubbleError::out_of_validity_range(
            "RK4 time step must be finite and non-negative",
            "dt_s must be >= 0",
        ));
    }
    let y0 = [bubble.position, bubble.velocity];
    let k1 = ode_rhs(y0, &mut acceleration)?;
    let k2 = ode_rhs(add_state(y0, k1, 0.5 * dt_s), &mut acceleration)?;
    let k3 = ode_rhs(add_state(y0, k2, 0.5 * dt_s), &mut acceleration)?;
    let k4 = ode_rhs(add_state(y0, k3, dt_s), &mut acceleration)?;
    let mut pos = bubble.position;
    let mut vel = bubble.velocity;
    for a in 0..3 {
        pos[a] += dt_s * (k1[0][a] + 2.0 * k2[0][a] + 2.0 * k3[0][a] + k4[0][a]) / 6.0;
        vel[a] += dt_s * (k1[1][a] + 2.0 * k2[1][a] + 2.0 * k3[1][a] + k4[1][a]) / 6.0;
    }
    Ok((pos, vel))
}

pub fn terminal_velocity_schiller_naumann(
    diameter_m: f64,
    rho_liquid_kg_m3: f64,
    nu_liquid_m2_s: f64,
    gravity_m_s2: f64,
) -> BubbleResult<f64> {
    let volume = bubble_volume_from_diameter(diameter_m)?;
    let area = std::f64::consts::PI * diameter_m * diameter_m / 4.0;
    let buoyancy = rho_liquid_kg_m3 * gravity_m_s2 * volume;
    let mut v = gravity_m_s2 * diameter_m * diameter_m / (18.0 * nu_liquid_m2_s);
    for _ in 0..80 {
        let re = v * diameter_m / nu_liquid_m2_s;
        if re > SCHILLER_NAUMANN_RE_MAX {
            return Err(BubbleError::out_of_validity_range(
                "terminal velocity exceeds Schiller-Naumann validity range",
                format!("Re_bubble={re:.6} exceeds {SCHILLER_NAUMANN_RE_MAX:.1}"),
            ));
        }
        let cd = if re == 0.0 {
            0.0
        } else {
            (24.0 / re) * (1.0 + 0.15 * re.powf(0.687))
        };
        let next_v = (2.0 * buoyancy / (rho_liquid_kg_m3 * cd * area)).sqrt();
        if (next_v - v).abs() < 1.0e-14 {
            return Ok(next_v);
        }
        v = next_v;
    }
    Ok(v)
}

fn ode_rhs<F>(state: [[f64; 3]; 2], acceleration: &mut F) -> BubbleResult<[[f64; 3]; 2]>
where
    F: FnMut([f64; 3], [f64; 3]) -> BubbleResult<[f64; 3]>,
{
    Ok([state[1], acceleration(state[0], state[1])?])
}

fn add_state(state: [[f64; 3]; 2], k: [[f64; 3]; 2], h: f64) -> [[f64; 3]; 2] {
    let mut out = state;
    for row in 0..2 {
        for a in 0..3 {
            out[row][a] += h * k[row][a];
        }
    }
    out
}

fn norm(v: [f64; 3]) -> f64 {
    dot(v, v).sqrt()
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

fn scale(v: [f64; 3], s: f64) -> [f64; 3] {
    [s * v[0], s * v[1], s * v[2]]
}

fn unit(v: [f64; 3]) -> BubbleResult<[f64; 3]> {
    let n = norm(v);
    if n == 0.0 {
        return Err(BubbleError::out_of_validity_range(
            "normal vector must be non-zero",
            "wall_normal_into_liquid must have non-zero magnitude",
        ));
    }
    Ok([v[0] / n, v[1] / n, v[2] / n])
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bubbles::Bubble;

    fn ctx(re: f64) -> BubbleForceContext {
        BubbleForceContext {
            re_bubble: re,
            alpha_g: 0.01,
            kolmogorov_weber: 1.0,
        }
    }

    #[test]
    fn single_rising_bubble_terminal_velocity_matches_force_balance() {
        let d = 1.0e-3;
        let rho = 1000.0;
        let nu = 1.0e-6;
        let g = 9.81;
        let vt = terminal_velocity_schiller_naumann(d, rho, nu, g).unwrap();
        let re = vt * d / nu;
        let drag = schiller_naumann_drag_force(rho, nu, d, [0.0, 0.0, vt], ctx(re)).unwrap();
        let buoy = buoyancy_force(
            rho,
            [0.0, 0.0, -g],
            bubble_volume_from_diameter(d).unwrap(),
            ctx(re),
        )
        .unwrap();
        let rel = ((drag[2] - buoy[2]) / buoy[2]).abs();
        let dir = std::path::Path::new("target/bcfd_071");
        std::fs::create_dir_all(dir).unwrap();
        std::fs::write(
            dir.join("rising_bubble_terminal_velocity.csv"),
            format!("diameter_m,terminal_velocity_m_s,re_bubble,force_balance_rel_error\n{d:e},{vt:e},{re:e},{rel:e}\n"),
        )
        .unwrap();
        assert!(
            rel < 0.05,
            "terminal velocity force balance rel error {rel:e} must be < 5%; vt={vt:e}, Re={re:e}"
        );
        assert!(vt > 0.0, "rising terminal speed must be positive");
    }

    #[test]
    fn force_closures_return_finite_values() {
        let c = ctx(10.0);
        let vb = bubble_volume_from_diameter(1.0e-3).unwrap();
        let forces = [
            buoyancy_force(1000.0, [0.0, 0.0, -9.81], vb, c).unwrap(),
            schiller_naumann_drag_force(1000.0, 1.0e-6, 1.0e-3, [0.01, 0.0, 0.0], c).unwrap(),
            added_mass_force(1000.0, vb, [1.0, 0.0, 0.0], [0.0; 3], c).unwrap(),
            lift_placeholder_force(1000.0, vb, [0.01, 0.0, 0.0], [0.0, 0.0, 2.0], c).unwrap(),
            wall_lubrication_placeholder_force(
                1000.0,
                vb,
                [0.01, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                1.0e-3,
                1.0e-3,
                c,
            )
            .unwrap(),
            turbulent_dispersion_placeholder_force(1000.0, vb, [1.0, 0.0, 0.0], 1.0e-4, c).unwrap(),
        ];
        for f in forces {
            assert!(f.iter().all(|v| v.is_finite()));
        }
    }

    #[test]
    fn invalid_re_is_rejected() {
        let err = schiller_naumann_drag_force(1000.0, 1.0e-6, 1.0e-3, [1.0, 0.0, 0.0], ctx(1000.0))
            .unwrap_err();
        assert!(matches!(
            err.reason,
            crate::solver::UnsupportedReason::OutOfValidityRange { .. }
        ));
    }

    #[test]
    fn rk4_substep_is_deterministic() {
        let b = Bubble::new([0.0; 3], [1.0, 0.0, 0.0], 1.0e-3, 1.0, 1).unwrap();
        let acc = |_: [f64; 3], _: [f64; 3]| Ok([0.0, 0.0, -9.81]);
        let a = rk4_substep(&b, 0.01, acc).unwrap();
        let acc = |_: [f64; 3], _: [f64; 3]| Ok([0.0, 0.0, -9.81]);
        let b2 = rk4_substep(&b, 0.01, acc).unwrap();
        assert_eq!(a, b2);
        assert!(a.0[0] > 0.0);
        assert!(a.1[2] < 0.0);
    }
}
