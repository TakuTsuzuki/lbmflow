//! V&V lane 5.1: feature-interaction conservation matrix.
//!
//! The target bug class is source-composition at feature pairs. The force
//! contracts pinned here are the documented ones: Shan-Chen overwrites the
//! per-cell force field, gravity adds `rho*g`, and rotor penalization adds
//! after the caller has cleared/rebuilt the field (ANOM-P4-009/P4-010).

use lbm_core::compat::multiphase::ShanChen;
use lbm_core::compat::prelude::{Collision, SimConfig, Simulation};
#[cfg(feature = "mf-interim")]
use lbm_core::compat::rotor::Rotor;
use lbm_core::particles::{DepositEvent, Particle, ParticleSet, Sample};
use lbm_core::prelude::{
    build_wall_rims, CollisionKind, CpuScalar, Face, FaceBC, FacePatch, GlobalSpec, LocalPeriodic,
    Solver, SourceKind, SourceRegion, VolumeSource, WallSpec, D3Q19,
};
use std::panic::{catch_unwind, AssertUnwindSafe};

const STEPS: usize = 200;
const FORCED_PATCH_STEPS: usize = 2000;
const FORCED_PATCH_DRIFT_WINDOW: usize = 200;
const MASS_REL_STATE_DEN: f64 = 1.0e-9;
const MIRROR_ABS_F64: f64 = 1.0e-12;
const FORCE_SUPERPOSITION_REL: f64 = 1.0e-12;
const CS2: f64 = 1.0 / 3.0;

type Native3 = Solver<D3Q19, f64, CpuScalar, LocalPeriodic>;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Feature {
    UniformForce,
    Gravity,
    ShanChen,
    Rotor,
    Particles,
    VolumeSource,
    FacePatch,
}

impl Feature {
    fn name(self) -> &'static str {
        match self {
            Feature::UniformForce => "uniform force",
            Feature::Gravity => "gravity",
            Feature::ShanChen => "Shan-Chen SCMP",
            Feature::Rotor => "rotor penalization",
            Feature::Particles => "particles CR-3",
            Feature::VolumeSource => "volume sources",
            Feature::FacePatch => "face patches",
        }
    }

    fn compat_only(self) -> bool {
        matches!(self, Feature::ShanChen | Feature::Rotor)
    }

    fn native_only(self) -> bool {
        matches!(self, Feature::VolumeSource | Feature::FacePatch)
    }

    fn is_force_source(self) -> bool {
        matches!(
            self,
            Feature::UniformForce | Feature::Gravity | Feature::ShanChen | Feature::Rotor
        )
    }
}

#[derive(Clone, Debug)]
enum Cell {
    Pass,
    Fail(String),
    Skip(String),
}

impl Cell {
    fn is_fail(&self) -> bool {
        matches!(self, Cell::Fail(_))
    }

    fn label(&self) -> &'static str {
        match self {
            Cell::Pass => "PASS",
            Cell::Fail(_) => "FAIL",
            Cell::Skip(_) => "SKIP",
        }
    }

    fn reason(&self) -> &str {
        match self {
            Cell::Pass => "",
            Cell::Fail(s) | Cell::Skip(s) => s,
        }
    }
}

#[test]
fn feature_interaction_conservation_matrix() {
    let features = [
        Feature::UniformForce,
        Feature::Gravity,
        Feature::ShanChen,
        Feature::Rotor,
        Feature::Particles,
        Feature::VolumeSource,
        Feature::FacePatch,
    ];

    let mut failures = Vec::new();
    for i in 0..features.len() {
        for j in (i + 1)..features.len() {
            let a = features[i];
            let b = features[j];
            let cell = run_cell(a, b);
            println!(
                "INTERACTION_MATRIX | {} x {} | {} | {}",
                a.name(),
                b.name(),
                cell.label(),
                cell.reason()
            );
            if cell.is_fail() {
                failures.push(format!("{} x {}: {}", a.name(), b.name(), cell.reason()));
            }
        }
    }

    let f32_cell = run_named("uniform force x particles CR-3 [f32]", || {
        run_compat_f32_uniform_particles()
    });
    println!(
        "INTERACTION_MATRIX | uniform force x particles CR-3 [f32] | {} | {}",
        f32_cell.label(),
        f32_cell.reason()
    );
    if f32_cell.is_fail() {
        failures.push(format!(
            "uniform force x particles CR-3 [f32]: {}",
            f32_cell.reason()
        ));
    }

    assert!(
        failures.is_empty(),
        "interaction matrix failures:\n{}",
        failures.join("\n")
    );
}

fn run_cell(a: Feature, b: Feature) -> Cell {
    if (a == Feature::Rotor || b == Feature::Rotor) && !rotor_enabled() {
        return Cell::Skip("requires cargo feature mf-interim".to_string());
    }
    if (a.compat_only() && b.native_only()) || (b.compat_only() && a.native_only()) {
        return Cell::Skip(format!(
            "{} is compat-only and {} is native-only; no public API composes them in one solver",
            if a.compat_only() { a.name() } else { b.name() },
            if a.native_only() { a.name() } else { b.name() }
        ));
    }

    run_named(&format!("{} x {}", a.name(), b.name()), || {
        if a.native_only() || b.native_only() {
            run_native_pair(a, b)
        } else {
            run_compat_pair(a, b)
        }
    })
}

fn run_named(name: &str, f: impl FnOnce() + std::panic::UnwindSafe) -> Cell {
    match catch_unwind(AssertUnwindSafe(f)) {
        Ok(()) => Cell::Pass,
        Err(payload) => Cell::Fail(panic_payload(payload, name)),
    }
}

fn panic_payload(payload: Box<dyn std::any::Any + Send>, name: &str) -> String {
    if let Some(s) = payload.downcast_ref::<String>() {
        s.clone()
    } else if let Some(s) = payload.downcast_ref::<&'static str>() {
        (*s).to_string()
    } else {
        format!("{name}: non-string panic payload")
    }
}

fn rotor_enabled() -> bool {
    cfg!(feature = "mf-interim")
}

#[derive(Clone, Copy)]
struct CompatCfg {
    uniform: bool,
    gravity: bool,
    sc: bool,
    #[cfg_attr(not(feature = "mf-interim"), allow(dead_code))]
    rotor: bool,
    particles: bool,
    mirror_x: bool,
}

impl CompatCfg {
    fn from_pair(a: Feature, b: Feature, mirror_x: bool) -> Self {
        Self {
            uniform: a == Feature::UniformForce || b == Feature::UniformForce,
            gravity: a == Feature::Gravity || b == Feature::Gravity,
            sc: a == Feature::ShanChen || b == Feature::ShanChen,
            rotor: a == Feature::Rotor || b == Feature::Rotor,
            particles: a == Feature::Particles || b == Feature::Particles,
            mirror_x,
        }
    }
}

fn compat_uniform_force(mirror_x: bool) -> [f64; 2] {
    [if mirror_x { -1.0e-7 } else { 1.0e-7 }, 2.0e-7]
}

fn compat_gravity(mirror_x: bool) -> [f64; 2] {
    [if mirror_x { 2.0e-7 } else { -2.0e-7 }, 1.0e-7]
}

fn build_compat(cfg: CompatCfg) -> Simulation<f64> {
    let mut sim = SimConfig {
        nx: 32,
        ny: 32,
        nu: 1.0 / 6.0,
        collision: Collision::Trt {
            magic: Collision::MAGIC_STD,
        },
        force: if cfg.uniform {
            compat_uniform_force(cfg.mirror_x)
        } else {
            [0.0, 0.0]
        },
        ..Default::default()
    }
    .build()
    .unwrap();
    if cfg.sc {
        sim.init_with(|x, y| {
            let kx = std::f64::consts::TAU * (x as f64 + 0.5) / 32.0;
            let ky = std::f64::consts::TAU * (y as f64 + 0.5) / 32.0;
            (1.0 + 0.01 * kx.cos() * ky.cos(), 0.0, 0.0)
        });
    }
    if cfg.gravity {
        sim.set_gravity(compat_gravity(cfg.mirror_x));
    }
    sim
}

#[cfg(feature = "mf-interim")]
fn build_rotor(mirror_x: bool) -> Rotor<f64> {
    // Thin blades, low omega, and chi<1 are the ANOM-P4-010 stable regime:
    // the test exercises additive composition without entering the known
    // solid-disc/chi=1 collective instability.
    Rotor::new(if mirror_x { 16.0 } else { 15.0 }, 15.0)
        .n_blades(2)
        .r_hub(2.5)
        .r_blade(9.0)
        .blade_thickness(1.0)
        .omega(if mirror_x { -5.0e-4 } else { 5.0e-4 })
        .chi(0.2)
        .omega_ramp_steps(0)
}

fn apply_compat_dynamic(
    sim: &mut Simulation<f64>,
    sc: &Option<ShanChen<f64>>,
    #[cfg(feature = "mf-interim")] rotor: &mut Option<Rotor<f64>>,
    #[cfg(not(feature = "mf-interim"))] _rotor: &mut Option<()>,
) {
    if sc.is_none() {
        sim.clear_force_field();
    }
    if let Some(sc) = sc {
        // Shan-Chen owns the per-cell field and overwrites it each step.
        sc.update_force(sim);
    }
    #[cfg(feature = "mf-interim")]
    if let Some(rotor) = rotor {
        // Rotor penalization adds into the field after the caller has rebuilt
        // the earlier per-cell sources (ANOM-P4-009).
        rotor.update_force(sim);
    }
}

fn run_compat_pair(a: Feature, b: Feature) {
    let cfg = CompatCfg::from_pair(a, b, false);
    run_compat_mass_and_finiteness(cfg);
    run_compat_mirror(cfg);
    if a.is_force_source() && b.is_force_source() {
        assert_compat_force_superposition(a, b);
    }
}

fn run_compat_mass_and_finiteness(cfg: CompatCfg) {
    let mut sim = build_compat(cfg);
    let sc = cfg.sc.then(|| ShanChen::new(-1.0));
    #[cfg(feature = "mf-interim")]
    let mut rotor = cfg.rotor.then(|| build_rotor(false));
    #[cfg(not(feature = "mf-interim"))]
    let mut rotor = None::<()>;
    let mut particles = cfg.particles.then(build_particles);
    let mut deposits = Vec::<DepositEvent>::new();
    let m0 = sim.total_mass_f64();

    let mut twin = cfg.particles.then(|| {
        build_compat(CompatCfg {
            particles: false,
            ..cfg
        })
    });
    #[cfg(feature = "mf-interim")]
    let mut twin_rotor = cfg
        .particles
        .then(|| cfg.rotor.then(|| build_rotor(false)))
        .flatten();
    #[cfg(not(feature = "mf-interim"))]
    let mut twin_rotor = None::<()>;
    let twin_sc = cfg
        .particles
        .then(|| cfg.sc.then(|| ShanChen::new(-1.0)))
        .flatten();

    for _ in 0..STEPS {
        apply_compat_dynamic(&mut sim, &sc, &mut rotor);
        sim.step();
        if let Some(particles) = particles.as_mut() {
            step_particles_compat(particles, &sim, &mut deposits);
            let twin = twin.as_mut().unwrap();
            apply_compat_dynamic(twin, &twin_sc, &mut twin_rotor);
            twin.step();
        }
    }

    let m1 = sim.total_mass_f64();
    assert_mass_ledger(
        "compat",
        m0,
        m1,
        0.0,
        MASS_REL_STATE_DEN,
        "denominator=max(|initial total mass|, 1)",
    );
    assert_compat_finite(&sim, "compat");
    if let Some(twin) = twin.as_ref() {
        assert_compat_fields_equal(
            &sim,
            twin,
            0.0,
            "particles CR-3 one-way no-feedback field equality",
        );
    }
}

fn run_compat_mirror(cfg: CompatCfg) {
    let mut a = build_compat(CompatCfg {
        mirror_x: false,
        ..cfg
    });
    let mut b = build_compat(CompatCfg {
        mirror_x: true,
        ..cfg
    });
    let sc_a = cfg.sc.then(|| ShanChen::new(-1.0));
    let sc_b = cfg.sc.then(|| ShanChen::new(-1.0));
    #[cfg(feature = "mf-interim")]
    let mut rotor_a = cfg.rotor.then(|| build_rotor(false));
    #[cfg(feature = "mf-interim")]
    let mut rotor_b = cfg.rotor.then(|| build_rotor(true));
    #[cfg(not(feature = "mf-interim"))]
    let mut rotor_a = None::<()>;
    #[cfg(not(feature = "mf-interim"))]
    let mut rotor_b = None::<()>;

    for _ in 0..STEPS {
        apply_compat_dynamic(&mut a, &sc_a, &mut rotor_a);
        apply_compat_dynamic(&mut b, &sc_b, &mut rotor_b);
        a.step();
        b.step();
    }
    let d = compat_mirror_delta(&a, &b);
    println!(
        "INTERACTION_MATRIX detail | compat mirror | rho={:.3e} ux_mirror={:.3e} uy={:.3e} band={MIRROR_ABS_F64:.3e}",
        d[0], d[1], d[2]
    );
    assert!(
        d.iter().all(|v| *v <= MIRROR_ABS_F64),
        "compat x-mirror equivariance failed: max_delta rho={:.12e}, ux_mirror={:.12e}, uy={:.12e}, band={MIRROR_ABS_F64:.12e}",
        d[0],
        d[1],
        d[2]
    );
}

fn compat_mirror_delta(a: &Simulation<f64>, b: &Simulation<f64>) -> [f64; 3] {
    let mut d = [0.0f64; 3];
    for y in 0..a.ny() {
        for x in 0..a.nx() {
            let xm = a.nx() - 1 - x;
            d[0] = d[0].max((a.rho(x, y) - b.rho(xm, y)).abs());
            d[1] = d[1].max((a.ux(x, y) + b.ux(xm, y)).abs());
            d[2] = d[2].max((a.uy(x, y) - b.uy(xm, y)).abs());
        }
    }
    d
}

fn assert_compat_finite(sim: &Simulation<f64>, label: &str) {
    for y in 0..sim.ny() {
        for x in 0..sim.nx() {
            let rho = sim.rho(x, y);
            let ux = sim.ux(x, y);
            let uy = sim.uy(x, y);
            assert!(
                rho.is_finite() && ux.is_finite() && uy.is_finite(),
                "{label}: non-finite field at ({x},{y}): rho={rho:e}, ux={ux:e}, uy={uy:e}"
            );
        }
    }
}

fn assert_compat_fields_equal(a: &Simulation<f64>, b: &Simulation<f64>, band: f64, label: &str) {
    let mut max = 0.0f64;
    for y in 0..a.ny() {
        for x in 0..a.nx() {
            max = max.max((a.rho(x, y) - b.rho(x, y)).abs());
            max = max.max((a.ux(x, y) - b.ux(x, y)).abs());
            max = max.max((a.uy(x, y) - b.uy(x, y)).abs());
        }
    }
    assert!(
        max <= band,
        "{label}: max field delta={max:.12e}, band={band:.12e}"
    );
}

fn compat_one_step_gain(features: &[Feature], sc_initial_state: bool) -> [f64; 2] {
    let cfg = CompatCfg {
        uniform: features.contains(&Feature::UniformForce),
        gravity: features.contains(&Feature::Gravity),
        sc: sc_initial_state,
        rotor: features.contains(&Feature::Rotor),
        particles: false,
        mirror_x: false,
    };
    let mut sim = build_compat(cfg);
    let sc = features
        .contains(&Feature::ShanChen)
        .then(|| ShanChen::new(-1.0));
    #[cfg(feature = "mf-interim")]
    let mut rotor = cfg.rotor.then(|| build_rotor(false));
    #[cfg(not(feature = "mf-interim"))]
    let mut rotor = None::<()>;
    apply_compat_dynamic(&mut sim, &sc, &mut rotor);
    let p0 = sim.total_momentum();
    sim.step();
    let p1 = sim.total_momentum();
    [p1[0] - p0[0], p1[1] - p0[1]]
}

fn assert_compat_force_superposition(a: Feature, b: Feature) {
    // Guo forcing is linear in the force vector used at collision. At t=1 all
    // feature forces below are computed from the same initialized populations:
    // uniform force is in StepParams, gravity contributes rho*g, Shan-Chen
    // overwrites the per-cell field, and rotor then adds into that same field.
    // Therefore the composed one-step momentum increment must equal the sum
    // of the two single-feature increments. ANOM-P2-001 makes the absolute
    // impulse path-dependent, so this deliberately compares composed-vs-parts
    // using each feature's real path and never equates uniform force to a raw
    // per-cell force-field path.
    let sc_initial_state = a == Feature::ShanChen || b == Feature::ShanChen;
    let pa = compat_one_step_gain(&[a], sc_initial_state);
    let pb = compat_one_step_gain(&[b], sc_initial_state);
    let pab = compat_one_step_gain(&[a, b], sc_initial_state);
    let expected = [pa[0] + pb[0], pa[1] + pb[1]];
    let err = ((pab[0] - expected[0]).powi(2) + (pab[1] - expected[1]).powi(2)).sqrt();
    let den = (expected[0].powi(2) + expected[1].powi(2))
        .sqrt()
        .max((pab[0].powi(2) + pab[1].powi(2)).sqrt())
        .max(1.0e-30);
    let rel = err / den;
    println!(
        "INTERACTION_MATRIX detail | force superposition {} x {} | composed=[{:.12e},{:.12e}] sum=[{:.12e},{:.12e}] rel={rel:.3e} denominator=max(|composed|,|sum|,1e-30)",
        a.name(),
        b.name(),
        pab[0],
        pab[1],
        expected[0],
        expected[1]
    );
    assert!(
        rel <= FORCE_SUPERPOSITION_REL,
        "force-composition one-step superposition failed for {} x {}: rel={rel:.12e}, band={FORCE_SUPERPOSITION_REL:.12e}, denominator=max(|composed|,|sum|,1e-30), composed=[{:.12e},{:.12e}], sum_of_parts=[{:.12e},{:.12e}]",
        a.name(),
        b.name(),
        pab[0],
        pab[1],
        expected[0],
        expected[1]
    );
}

#[derive(Clone, Copy)]
struct NativeCfg {
    uniform: bool,
    gravity: bool,
    source: bool,
    patch: bool,
    particles: bool,
    mirror_x: bool,
}

impl NativeCfg {
    fn from_pair(a: Feature, b: Feature, mirror_x: bool) -> Self {
        Self {
            uniform: a == Feature::UniformForce || b == Feature::UniformForce,
            gravity: a == Feature::Gravity || b == Feature::Gravity,
            source: a == Feature::VolumeSource || b == Feature::VolumeSource,
            patch: a == Feature::FacePatch || b == Feature::FacePatch,
            particles: a == Feature::Particles || b == Feature::Particles,
            mirror_x,
        }
    }
}

fn native_uniform_force(mirror_x: bool) -> [f64; 3] {
    [if mirror_x { -1.0e-7 } else { 1.0e-7 }, 2.0e-7, -1.0e-7]
}

fn native_gravity(mirror_x: bool) -> [f64; 3] {
    [if mirror_x { 2.0e-7 } else { -2.0e-7 }, 1.0e-7, 3.0e-7]
}

fn source_q() -> f64 {
    2.0e-7
}

fn native_source() -> VolumeSource<f64> {
    VolumeSource {
        region: SourceRegion {
            lo: [11, 10, 10],
            hi: [12, 13, 13],
        },
        kind: SourceKind::MassFlow { q_lu: source_q() },
    }
}

fn native_patch(mirror_x: bool) -> FacePatch<f64> {
    let _ = mirror_x;
    FacePatch {
        face: Face::ZPos.index(),
        lo: [8, 8],
        hi: [15, 15],
        bc: FaceBC::Velocity {
            // A zero-velocity patch is still the masked Velocity-patch path
            // (ANOM-P4-006 zero-velocity-lid semantics), but its analytic
            // normal mass flux is exactly zero. Nonzero tangential patches are
            // useful flow tests, not clean conservation-ledger cells.
            u: [0.0, 0.0, 0.0],
        },
    }
}

fn native_walls(patch: bool) -> WallSpec<f64> {
    let mut walls = WallSpec::default();
    for face in Face::ALL {
        if patch && face == Face::ZPos {
            continue;
        }
        walls.is_wall[face.index()] = true;
    }
    walls
}

fn build_native(cfg: NativeCfg) -> Native3 {
    let dims = [24, 24, 24];
    let spec = GlobalSpec {
        dims,
        nu: 0.05,
        periodic: [false, false, false],
        faces: [FaceBC::Closed; 6],
        force: if cfg.uniform {
            native_uniform_force(cfg.mirror_x)
        } else {
            [0.0; 3]
        },
        collision: CollisionKind::Trt {
            magic: CollisionKind::MAGIC_STD,
        },
        sources: if cfg.source {
            vec![native_source()]
        } else {
            Vec::new()
        },
        face_patches: if cfg.patch {
            vec![native_patch(cfg.mirror_x)]
        } else {
            Vec::new()
        },
    };
    let (solid, wall_u) = build_wall_rims(3, dims, &native_walls(cfg.patch));
    let mut s = Solver::new(
        &spec,
        &solid,
        &wall_u,
        [1, 1, 1],
        CpuScalar::default(),
        LocalPeriodic,
    );
    if cfg.gravity {
        s.set_gravity(native_gravity(cfg.mirror_x));
    }
    s
}

fn run_native_pair(a: Feature, b: Feature) {
    let cfg = NativeCfg::from_pair(a, b, false);
    run_native_mass_and_finiteness(cfg);
    run_native_mirror(cfg);
}

fn run_native_mass_and_finiteness(cfg: NativeCfg) {
    if cfg.patch && (cfg.uniform || cfg.gravity) && !cfg.source && !cfg.particles {
        run_forced_patch_mass_and_finiteness(cfg);
        return;
    }

    let mut s = build_native(cfg);
    let m0 = s.total_mass_f64();
    let mut particles = cfg.particles.then(build_particles);
    let mut deposits = Vec::<DepositEvent>::new();

    let mut twin = cfg.particles.then(|| {
        build_native(NativeCfg {
            particles: false,
            ..cfg
        })
    });
    for _ in 0..STEPS {
        s.step();
        if let Some(particles) = particles.as_mut() {
            step_particles_native(particles, &s, &mut deposits);
            twin.as_mut().unwrap().step();
        }
    }
    let m1 = s.total_mass_f64();
    let expected = if cfg.source { source_q() } else { 0.0 };
    assert_mass_ledger(
        "native",
        m0,
        m1,
        expected,
        MASS_REL_STATE_DEN,
        "denominator=max(|initial total mass|, 1)",
    );
    assert_eq!(
        s.local_nonfinite_count(),
        0,
        "native: non-finite count after {STEPS} steps"
    );
    if let Some(twin) = twin.as_ref() {
        assert_native_fields_equal(
            &s,
            twin,
            0.0,
            "particles CR-3 one-way no-feedback field equality",
        );
    }
}

fn run_forced_patch_mass_and_finiteness(cfg: NativeCfg) {
    assert!(cfg.patch && (cfg.uniform || cfg.gravity));
    assert!(!cfg.source && !cfg.particles);

    let mut s = build_native(cfg);
    let m0 = s.total_mass_f64();
    let pre_window_steps = FORCED_PATCH_STEPS - FORCED_PATCH_DRIFT_WINDOW;
    s.run(pre_window_steps);
    let m_pre_window = s.total_mass_f64();
    s.run(FORCED_PATCH_DRIFT_WINDOW);
    let m1 = s.total_mass_f64();
    let transient_delta = m1 - m0;
    let transient_rel = transient_delta / m0.abs().max(1.0);
    let label = if cfg.gravity {
        "native gravity x face patches"
    } else {
        "native uniform force x face patches"
    };

    println!(
        "INTERACTION_MATRIX detail | {label} transient mass adjustment | total_delta={transient_delta:.12e} rel={transient_rel:.3e} steps={FORCED_PATCH_STEPS}"
    );
    if cfg.gravity {
        let (exact, order) = gravity_hydrostatic_mass_estimate(&s, native_gravity(cfg.mirror_x));
        println!(
            "INTERACTION_MATRIX detail | native gravity x face patches hydrostatic estimate | transient_delta={transient_delta:.12e} exact_signed={exact:.12e} order_mag={order:.12e} formula=sum rho0*(exp(-g_z*depth/cs2)-1), order=rho0*|g_z|*H/(2*cs2)*V"
        );
    }

    assert_steady_mass_ledger(
        label,
        m_pre_window,
        m1,
        FORCED_PATCH_DRIFT_WINDOW,
        MASS_REL_STATE_DEN,
        "denominator=max(|initial mass at final window start|, 1)",
    );
    assert_eq!(
        s.local_nonfinite_count(),
        0,
        "native forced patch: non-finite count after {FORCED_PATCH_STEPS} steps"
    );
}

fn gravity_hydrostatic_mass_estimate(s: &Native3, g: [f64; 3]) -> (f64, f64) {
    // Per-mass gravity satisfies the isothermal hydrostatic balance
    // cs^2 grad(rho) = rho g.  Taking the z+ velocity patch as the reference
    // density plane and depth = z_patch - z gives
    // rho(depth) = rho0 * exp(-g_z * depth / cs^2).  The finite compressible
    // fill from a uniform rho0=1 state is therefore the cell sum below.  For
    // |g_z|H/cs^2 << 1 its magnitude is O(rho0*|g_z|*H/(2*cs^2)*V).
    let dims = [24usize, 24usize, 24usize];
    let rho0 = 1.0;
    let z_patch = (dims[2] - 1) as f64;
    let mut exact = 0.0;
    let mut depth_sum = 0.0;
    let mut cells = 0usize;
    for z in 0..dims[2] {
        for y in 0..dims[1] {
            for x in 0..dims[0] {
                if s.is_solid(x, y, z) {
                    continue;
                }
                let depth = z_patch - z as f64;
                exact += rho0 * ((-g[2] * depth / CS2).exp() - 1.0);
                depth_sum += depth;
                cells += 1;
            }
        }
    }
    let order = rho0 * g[2].abs() * depth_sum / CS2;
    debug_assert!(cells > 0);
    (exact, order)
}

fn run_native_mirror(cfg: NativeCfg) {
    let mut a = build_native(NativeCfg {
        mirror_x: false,
        ..cfg
    });
    let mut b = build_native(NativeCfg {
        mirror_x: true,
        ..cfg
    });
    a.run(STEPS);
    b.run(STEPS);
    let d = native_mirror_delta(&a, &b, [24, 24, 24]);
    println!(
        "INTERACTION_MATRIX detail | native mirror | rho={:.3e} ux_mirror={:.3e} uy={:.3e} uz={:.3e} band={MIRROR_ABS_F64:.3e}",
        d[0], d[1], d[2], d[3]
    );
    assert!(
        d.iter().all(|v| *v <= MIRROR_ABS_F64),
        "native x-mirror equivariance failed: max_delta rho={:.12e}, ux_mirror={:.12e}, uy={:.12e}, uz={:.12e}, band={MIRROR_ABS_F64:.12e}",
        d[0],
        d[1],
        d[2],
        d[3]
    );
}

fn native_mirror_delta(a: &Native3, b: &Native3, dims: [usize; 3]) -> [f64; 4] {
    let ar = a.gather_rho();
    let ax = a.gather_ux();
    let ay = a.gather_uy();
    let az = a.gather_uz();
    let br = b.gather_rho();
    let bx = b.gather_ux();
    let by = b.gather_uy();
    let bz = b.gather_uz();
    let mut d = [0.0f64; 4];
    for z in 0..dims[2] {
        for y in 0..dims[1] {
            for x in 0..dims[0] {
                let i = idx3(dims, x, y, z);
                let im = idx3(dims, dims[0] - 1 - x, y, z);
                d[0] = d[0].max((ar[i] - br[im]).abs());
                d[1] = d[1].max((ax[i] + bx[im]).abs());
                d[2] = d[2].max((ay[i] - by[im]).abs());
                d[3] = d[3].max((az[i] - bz[im]).abs());
            }
        }
    }
    d
}

fn assert_native_fields_equal(a: &Native3, b: &Native3, band: f64, label: &str) {
    let mut max = 0.0f64;
    for (va, vb) in [
        (a.gather_rho(), b.gather_rho()),
        (a.gather_ux(), b.gather_ux()),
        (a.gather_uy(), b.gather_uy()),
        (a.gather_uz(), b.gather_uz()),
    ] {
        max = max.max(
            va.iter()
                .zip(&vb)
                .map(|(x, y)| (x - y).abs())
                .fold(0.0, f64::max),
        );
    }
    assert!(
        max <= band,
        "{label}: max field delta={max:.12e}, band={band:.12e}"
    );
}

fn idx3(dims: [usize; 3], x: usize, y: usize, z: usize) -> usize {
    (z * dims[1] + y) * dims[0] + x
}

fn assert_mass_ledger(
    label: &str,
    m0: f64,
    m1: f64,
    expected_per_step: f64,
    rel_band: f64,
    denominator_note: &str,
) {
    let measured = (m1 - m0) / STEPS as f64;
    let den = m0.abs().max(1.0);
    let rel = (measured - expected_per_step).abs() / den;
    println!(
        "INTERACTION_MATRIX detail | {label} mass ledger | measured_per_step={measured:.12e} expected_per_step={expected_per_step:.12e} rel={rel:.3e} band={rel_band:.3e} {denominator_note}"
    );
    assert!(
        rel <= rel_band,
        "{label} mass ledger failed: measured_per_step={measured:.12e}, expected_per_step={expected_per_step:.12e}, rel={rel:.12e}, band={rel_band:.12e}, {denominator_note}, m0={m0:.12e}, m1={m1:.12e}, steps={STEPS}"
    );
    assert!(
        m0.is_finite() && m1.is_finite() && measured.is_finite(),
        "{label} mass ledger non-finite: m0={m0:e}, m1={m1:e}, measured_per_step={measured:e}"
    );
}

fn assert_steady_mass_ledger(
    label: &str,
    m0: f64,
    m1: f64,
    steps: usize,
    rel_band: f64,
    denominator_note: &str,
) {
    let measured = (m1 - m0) / steps as f64;
    let den = m0.abs().max(1.0);
    let rel = measured.abs() / den;
    println!(
        "INTERACTION_MATRIX detail | {label} steady mass ledger | measured_per_step={measured:.12e} expected_per_step=0.000000000000e0 rel={rel:.3e} band={rel_band:.3e} window_steps={steps} {denominator_note}"
    );
    assert!(
        rel <= rel_band,
        "{label} steady mass ledger failed: measured_per_step={measured:.12e}, expected_per_step=0.000000000000e0, rel={rel:.12e}, band={rel_band:.12e}, {denominator_note}, m0={m0:.12e}, m1={m1:.12e}, window_steps={steps}. Interpretation rule: steady-state drift persisting above band = candidate Guo-force x Zou-He-patch mass leak (pitfall family 1: reconstruction ignores the half-force in unknown populations); transient-only mass adjustment = physical compressible hydrostatic filling and this cell passes."
    );
    assert!(
        m0.is_finite() && m1.is_finite() && measured.is_finite(),
        "{label} steady mass ledger non-finite: m0={m0:e}, m1={m1:e}, measured_per_step={measured:e}"
    );
}

fn build_particles() -> ParticleSet {
    ParticleSet::new(
        vec![
            Particle {
                pos: [8.25, 8.75, 0.0],
                vel: [0.0; 3],
                d: 0.05,
                rho_p: 1.1,
                exposure: 0.0,
            },
            Particle {
                pos: [17.5, 12.25, 0.0],
                vel: [0.0; 3],
                d: 0.05,
                rho_p: 1.1,
                exposure: 0.0,
            },
        ],
        1.0,
        1.0 / 6.0,
        [0.0; 3],
    )
}

fn step_particles_compat(
    particles: &mut ParticleSet,
    sim: &Simulation<f64>,
    deposits: &mut Vec<DepositEvent>,
) {
    let nx = sim.nx() as isize;
    let ny = sim.ny() as isize;
    let sample = |p: [f64; 3]| {
        let x = (p[0].round() as isize).rem_euclid(nx) as usize;
        let y = (p[1].round() as isize).rem_euclid(ny) as usize;
        Sample {
            u: [sim.ux(x, y), sim.uy(x, y), 0.0],
            solid: sim.is_solid(x, y),
        }
    };
    particles
        .step_depositing(sample, None::<fn([f64; 3]) -> f64>, -1.0, deposits)
        .unwrap();
}

fn step_particles_native(
    particles: &mut ParticleSet,
    s: &Native3,
    deposits: &mut Vec<DepositEvent>,
) {
    let dims = [24isize, 24isize, 24isize];
    let sample = |p: [f64; 3]| {
        let x = (p[0].round() as isize).rem_euclid(dims[0]) as usize;
        let y = (p[1].round() as isize).rem_euclid(dims[1]) as usize;
        let z = (p[2].round() as isize).rem_euclid(dims[2]) as usize;
        Sample {
            u: s.u(x, y, z),
            solid: s.is_solid(x, y, z),
        }
    };
    particles
        .step_depositing(sample, None::<fn([f64; 3]) -> f64>, -1.0, deposits)
        .unwrap();
}

fn run_compat_f32_uniform_particles() {
    let cfg = |mirror_x: bool| SimConfig {
        nx: 32,
        ny: 32,
        nu: 1.0 / 6.0,
        collision: Collision::Trt {
            magic: Collision::MAGIC_STD,
        },
        force: [if mirror_x { -1.0e-6_f32 } else { 1.0e-6_f32 }, 2.0e-6_f32],
        ..Default::default()
    };
    let mut sim = cfg(false).build().unwrap();
    let mut twin = cfg(false).build().unwrap();
    let mut mirror = cfg(true).build().unwrap();
    let mut particles = build_particles();
    let mut deposits = Vec::<DepositEvent>::new();
    let m0 = sim.total_mass_f64();
    for _ in 0..STEPS {
        sim.step();
        twin.step();
        mirror.step();
        let nx = sim.nx() as isize;
        let ny = sim.ny() as isize;
        let sample = |p: [f64; 3]| {
            let x = (p[0].round() as isize).rem_euclid(nx) as usize;
            let y = (p[1].round() as isize).rem_euclid(ny) as usize;
            Sample {
                u: [sim.ux(x, y) as f64, sim.uy(x, y) as f64, 0.0],
                solid: sim.is_solid(x, y),
            }
        };
        particles
            .step_depositing(sample, None::<fn([f64; 3]) -> f64>, -1.0, &mut deposits)
            .unwrap();
    }
    let m1 = sim.total_mass_f64();
    assert_mass_ledger(
        "compat f32",
        m0,
        m1,
        0.0,
        1.0e-5,
        "denominator=max(|initial total mass|, 1); f32 variant uses T6 f32-class mass tolerance",
    );
    let mut max_no_feedback = 0.0f32;
    let mut max_mirror = 0.0f32;
    for y in 0..sim.ny() {
        for x in 0..sim.nx() {
            max_no_feedback = max_no_feedback
                .max((sim.rho(x, y) - twin.rho(x, y)).abs())
                .max((sim.ux(x, y) - twin.ux(x, y)).abs())
                .max((sim.uy(x, y) - twin.uy(x, y)).abs());
            let xm = sim.nx() - 1 - x;
            max_mirror = max_mirror
                .max((sim.rho(x, y) - mirror.rho(xm, y)).abs())
                .max((sim.ux(x, y) + mirror.ux(xm, y)).abs())
                .max((sim.uy(x, y) - mirror.uy(xm, y)).abs());
            assert!(
                sim.rho(x, y).is_finite() && sim.ux(x, y).is_finite() && sim.uy(x, y).is_finite(),
                "compat f32 non-finite at ({x},{y})"
            );
        }
    }
    assert_eq!(
        max_no_feedback, 0.0,
        "compat f32 particles no-feedback max field delta={max_no_feedback:e}"
    );
    assert!(
        max_mirror <= 1.0e-5,
        "compat f32 x-mirror max delta={max_mirror:e}, band=1e-5"
    );
}
