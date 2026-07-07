//! Accuracy audit for lattice opposite tables and collision-free streaming.
//!
//! The single-delta tests use the native fields/backend API directly. That is
//! the clean collision-free path: populate `SoaFields::f`, run the periodic
//! halo exchange, call `CpuScalar::stream` once, and inspect `SoaFields::ftmp`
//! before any swap or collision can change the injected population.

use lbm_core::prelude::*;

fn assert_opp_table<L: Lattice>(name: &str) {
    assert_eq!(L::C.len(), L::Q, "{name}: C length must equal Q");
    assert_eq!(L::W.len(), L::Q, "{name}: W length must equal Q");
    assert_eq!(L::OPP.len(), L::Q, "{name}: OPP length must equal Q");

    for q in 0..L::Q {
        let opp = L::OPP[q];
        assert!(opp < L::Q, "{name}: OPP[{q}]={opp} is outside 0..{}", L::Q);
        for a in 0..3 {
            assert_eq!(
                L::C[opp][a],
                -L::C[q][a],
                "{name}: c[OPP[{q}]][{a}] must equal -c[{q}][{a}]"
            );
        }
        assert_eq!(L::W[opp], L::W[q], "{name}: w[OPP[{q}]] must equal w[{q}]");
    }
}

fn stream_params() -> StepParams<f64> {
    StepParams {
        collision: CollisionKind::Bgk,
        omega_p: 1.0,
        omega_m: 1.0,
        force: [0.0; 3],
        gravity: None,
        faces: [FaceBC::Closed; 6],
        sources: Vec::new(),
        face_patches: Vec::new(),
    }
}

fn wrap_add(x: usize, dx: i8, n: usize) -> usize {
    (x as isize + dx as isize).rem_euclid(n as isize) as usize
}

fn assert_single_delta_streams<L: Lattice>(name: &str, dims: [usize; 3]) {
    let sub = Subdomain::monolithic(L::D, dims, [true, true, true]);
    let params = stream_params();
    let c0 = [0usize, 0usize, 0usize];

    for q in 0..L::Q {
        let mut fields: SoaFields<f64> = SoaFields::new(L::Q, sub.geom);
        let np = fields.plane_len();
        let source = fields.geom.pidx(c0[0], c0[1], c0[2]);
        fields.f[q * np + source] = 1.0;

        let mut backend = CpuScalar::default();
        <CpuScalar as Backend<L, f64>>::exchange_f(
            &mut backend,
            &LocalPeriodic,
            std::slice::from_ref(&sub),
            std::slice::from_mut(&mut fields),
        );
        <CpuScalar as Backend<L, f64>>::stream(
            &mut backend,
            &sub,
            &mut fields,
            &params,
            CellRange::full(&sub),
        );

        let c = L::C[q];
        let dst = [
            wrap_add(c0[0], c[0], dims[0]),
            wrap_add(c0[1], c[1], dims[1]),
            wrap_add(c0[2], c[2], dims[2]),
        ];
        let expected_slot = q * np + fields.geom.pidx(dst[0], dst[1], dst[2]);
        let mut nonzero_slots = 0usize;

        for (slot, &got) in fields.ftmp.iter().enumerate() {
            let want = if slot == expected_slot { 1.0 } else { 0.0 };
            let err = (got - want).abs();
            if got != 0.0 {
                nonzero_slots += 1;
            }
            assert!(
                err <= 1.0e-14,
                "{name}: q={q}, c={c:?}, source={c0:?}, expected_dst={dst:?}, \
                 slot={slot}, expected_slot={expected_slot}, got={got:.17e}, \
                 want={want:.17e}, abs_err={err:.3e}, band=1e-14"
            );
        }
        assert_eq!(
            nonzero_slots, 1,
            "{name}: q={q}, c={c:?} should leave exactly one nonzero ftmp slot"
        );
    }
}

#[test]
fn s1_opp_tables_are_structurally_consistent() {
    assert_opp_table::<D2Q9>("D2Q9");
    assert_opp_table::<D3Q19>("D3Q19");
    assert_opp_table::<D3Q27>("D3Q27");
}

#[test]
fn s2_single_delta_streams_once_d2q9() {
    assert_single_delta_streams::<D2Q9>("D2Q9", [8, 8, 1]);
}

#[test]
fn s2_single_delta_streams_once_d3q19() {
    assert_single_delta_streams::<D3Q19>("D3Q19", [8, 8, 8]);
}

#[test]
fn s2_single_delta_streams_once_d3q27() {
    assert_single_delta_streams::<D3Q27>("D3Q27", [6, 6, 6]);
}
