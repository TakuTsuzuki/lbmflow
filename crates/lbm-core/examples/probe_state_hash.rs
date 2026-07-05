//! Deterministic state hash over a zoo of configurations.
//!
//! Prints one line per scenario: a hex hash of the exact bit patterns of
//! (rho, ux, uy) plus the probed force after a fixed number of steps.
//! Two builds that claim to be step-for-step equivalent must print identical
//! tables (bit-exact reproduction, far stricter than the tolerance-based
//! validation suite). Used to verify kernel rework (e.g. pass fusion) does
//! not change results at all.
//!
//! Run: `cargo run --release --example probe_state_hash`

use lbm_core::multiphase::ShanChen;
use lbm_core::prelude::*;

/// FNV-1a over the raw bits of a float slice (via f64 bit patterns; exact
/// for both f32 and f64 fields).
fn hash_fields<T: Real>(h: &mut u64, fields: &[&[T]]) {
    const PRIME: u64 = 0x100000001b3;
    for field in fields {
        for v in *field {
            let bits = v.as_f64().to_bits();
            for byte in bits.to_le_bytes() {
                *h ^= byte as u64;
                *h = h.wrapping_mul(PRIME);
            }
        }
    }
}

fn finish<T: Real>(name: &str, sim: &Simulation<T>, extra: [T; 2]) {
    // Fields must be bit-identical at any thread count (per-cell arithmetic
    // is deterministic). The probed force is a parallel reduction, so its
    // last-ulp value may depend on the partitioning; hash it separately.
    let mut h = 0xcbf29ce484222325u64;
    hash_fields(&mut h, &[sim.rho_field(), sim.ux_field(), sim.uy_field()]);
    let mut hp = 0xcbf29ce484222325u64;
    hash_fields(&mut hp, &[&extra[..]]);
    println!("{name}: fields={h:016x} pf={hp:016x}");
}

fn periodic_tgv<T: Real>(name: &str, n: usize, steps: usize) {
    let mut sim: Simulation<T> = SimConfig {
        nx: n,
        ny: n,
        nu: 0.02,
        ..Default::default()
    }
    .build()
    .unwrap();
    let k = 2.0 * std::f64::consts::PI / n as f64;
    sim.init_with(|x, y| {
        (
            T::one(),
            T::r(0.04 * (k * y as f64).sin()),
            T::r(0.04 * (2.0 * k * x as f64).sin()),
        )
    });
    sim.run(steps);
    finish(name, &sim, [T::zero(); 2]);
}

fn cavity<T: Real>(name: &str, steps: usize) {
    let mut sim: Simulation<T> = SimConfig {
        nx: 96,
        ny: 64,
        nu: 0.02,
        edges: Edges {
            left: EdgeBC::BounceBack,
            right: EdgeBC::BounceBack,
            bottom: EdgeBC::BounceBack,
            top: EdgeBC::MovingWall { u: [T::r(0.1), T::zero()] },
        },
        ..Default::default()
    }
    .build()
    .unwrap();
    sim.run(steps);
    finish(name, &sim, [T::zero(); 2]);
}

fn forced_channel<T: Real>(name: &str, steps: usize) {
    // Periodic-x channel driven by a body force, walls top/bottom.
    let mut sim: Simulation<T> = SimConfig {
        nx: 64,
        ny: 48,
        nu: 0.05,
        force: [T::r(1e-5), T::zero()],
        edges: Edges {
            left: EdgeBC::Periodic,
            right: EdgeBC::Periodic,
            bottom: EdgeBC::BounceBack,
            top: EdgeBC::BounceBack,
        },
        ..Default::default()
    }
    .build()
    .unwrap();
    sim.run(steps);
    finish(name, &sim, [T::zero(); 2]);
}

fn cylinder_probe<T: Real>(name: &str, steps: usize) {
    // Inlet -> outflow past a staircase cylinder, with the momentum probe on.
    let mut sim: Simulation<T> = SimConfig {
        nx: 128,
        ny: 64,
        nu: 0.02,
        edges: Edges {
            left: EdgeBC::VelocityInlet { u: [T::r(0.08), T::zero()] },
            right: EdgeBC::Outflow,
            bottom: EdgeBC::BounceBack,
            top: EdgeBC::BounceBack,
        },
        ..Default::default()
    }
    .build()
    .unwrap();
    let is_cyl = |x: usize, y: usize| {
        let (dx, dy) = (x as f64 - 32.0, y as f64 - 32.0);
        dx * dx + dy * dy <= 8.0 * 8.0
    };
    sim.set_solid_region(is_cyl);
    sim.set_force_probe(is_cyl);
    sim.init_with(|_, _| (T::one(), T::r(0.08), T::zero()));
    sim.run(steps);
    let pf = sim.probed_force();
    finish(name, &sim, pf);
}

fn convective_outlet<T: Real>(name: &str, steps: usize) {
    let mut sim: Simulation<T> = SimConfig {
        nx: 96,
        ny: 48,
        nu: 0.03,
        edges: Edges {
            left: EdgeBC::VelocityInlet { u: [T::r(0.06), T::zero()] },
            right: EdgeBC::ConvectiveOutflow { u_conv: T::r(0.06) },
            bottom: EdgeBC::BounceBack,
            top: EdgeBC::BounceBack,
        },
        ..Default::default()
    }
    .build()
    .unwrap();
    sim.init_with(|_, _| (T::one(), T::r(0.06), T::zero()));
    sim.run(steps);
    finish(name, &sim, [T::zero(); 2]);
}

fn pressure_driven<T: Real>(name: &str, steps: usize) {
    let mut sim: Simulation<T> = SimConfig {
        nx: 96,
        ny: 32,
        nu: 0.05,
        edges: Edges {
            left: EdgeBC::PressureOutlet { rho: T::r(1.02) },
            right: EdgeBC::PressureOutlet { rho: T::r(0.99) },
            bottom: EdgeBC::BounceBack,
            top: EdgeBC::BounceBack,
        },
        ..Default::default()
    }
    .build()
    .unwrap();
    sim.run(steps);
    finish(name, &sim, [T::zero(); 2]);
}

fn shan_chen_droplet<T: Real>(name: &str, steps: usize) {
    // Per-cell force field path (multiphase); frozen droplet relaxation.
    let mut sim: Simulation<T> = SimConfig {
        nx: 64,
        ny: 64,
        nu: 1.0 / 6.0,
        ..Default::default()
    }
    .build()
    .unwrap();
    sim.init_with(|x, y| {
        let (dx, dy) = (x as f64 - 32.0, y as f64 - 32.0);
        let r = (dx * dx + dy * dy).sqrt();
        (T::r(if r < 12.0 { 2.0 } else { 0.15 }), T::zero(), T::zero())
    });
    let sc = ShanChen::new(-5.0);
    for _ in 0..steps {
        sc.update_force(&mut sim);
        sim.step();
    }
    finish(name, &sim, [T::zero(); 2]);
}

fn main() {
    // Large enough to hit the parallel path (>= PARALLEL_MIN_CELLS).
    periodic_tgv::<f32>("tgv-f32-256", 256, 200);
    periodic_tgv::<f64>("tgv-f64-256", 256, 200);
    // Serial path + wall rims and moving lid.
    cavity::<f32>("cavity-f32", 300);
    cavity::<f64>("cavity-f64", 300);
    forced_channel::<f64>("channel-force-f64", 300);
    cylinder_probe::<f32>("cylinder-probe-f32", 200);
    cylinder_probe::<f64>("cylinder-probe-f64", 200);
    convective_outlet::<f64>("convective-f64", 250);
    pressure_driven::<f64>("zouhe-pressure-f64", 250);
    shan_chen_droplet::<f64>("shanchen-f64", 150);
}
