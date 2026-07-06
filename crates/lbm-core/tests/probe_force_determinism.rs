//! Probe-force repeatability for the B-2 backend synchronization contract.
//!
//! The CPU backends may compute row/band partials in parallel, but the public
//! probed force must be a fixed-order fold and bit-repeatable across runs at
//! the same Rayon thread count.

use lbm_core::prelude::*;

type S<B> = Solver<D2Q9, f32, B, LocalPeriodic>;

fn build<B>(backend: B) -> S<B>
where
    B: Backend<D2Q9, f32, Fields = SoaFields<f32>>,
{
    let (nx, ny) = (192usize, 112usize);
    let mut walls = WallSpec::<f32>::default();
    walls.is_wall[Face::YNeg.index()] = true;
    walls.is_wall[Face::YPos.index()] = true;
    let mut faces = [FaceBC::Closed; 6];
    faces[Face::XNeg.index()] = FaceBC::Velocity {
        u: [0.045, 0.0, 0.0],
    };
    faces[Face::XPos.index()] = FaceBC::Outflow;
    let spec = GlobalSpec {
        dims: [nx, ny, 1],
        nu: 0.03,
        periodic: [false, false, false],
        faces,
        ..Default::default()
    };
    let (solid, wall_u) = build_wall_rims(2, spec.dims, &walls);
    let mut s = Solver::new(&spec, &solid, &wall_u, [1, 1, 1], backend, LocalPeriodic);
    let (cx, cy, r) = (72.4f64, 55.7f64, 9.8f64);
    let inside = move |x: usize, y: usize, _: usize| {
        let dx = x as f64 - cx;
        let dy = y as f64 - cy;
        dx * dx + dy * dy < r * r
    };
    for y in 0..ny {
        for x in 0..nx {
            if inside(x, y, 0) {
                s.set_solid(x, y, 0);
            }
        }
    }
    s.set_force_probe(inside);
    s
}

fn bits(v: [f32; 3]) -> [u32; 3] {
    [v[0].to_bits(), v[1].to_bits(), v[2].to_bits()]
}

fn assert_backend_repeatable<B>(backend: B, name: &str)
where
    B: Backend<D2Q9, f32, Fields = SoaFields<f32>> + Copy,
{
    let mut a = build(backend);
    let mut b = build(backend);
    assert_eq!(
        bits(a.read_probed_force()),
        bits(b.read_probed_force()),
        "{name}: initial explicit readback differs"
    );
    for step in 1..=80 {
        a.step();
        b.step();
        let ac = bits(a.probed_force());
        let bc = bits(b.probed_force());
        assert_eq!(ac, bc, "{name}: cached probed_force differs at step {step}");
        assert_eq!(
            ac,
            bits(a.read_probed_force()),
            "{name}: cached/readback force differs at step {step}"
        );
    }
}

#[test]
fn cpu_scalar_probe_force_bit_repeatable() {
    assert_backend_repeatable(CpuScalar::default(), "CpuScalar");
}

#[test]
fn cpu_simd_probe_force_bit_repeatable() {
    assert_backend_repeatable(CpuSimd::default(), "CpuSimd");
}
