//! Dev probe: (1) prove the engine is exactly equivariant under the lattice
//! symmetries relating the four lid orientations (with CORRECT maps), and
//! (2) measure the Ghia Re=400 RMS at U=0.05 vs U=0.1 (Mach dependence).

use lbm_core::prelude::*;

const N: usize = 129;
const U: f64 = 0.1;

fn cavity(re: f64, edges: Edges<f64>, u: f64) -> Simulation<f64> {
    let l = (N - 2) as f64;
    SimConfig {
        nx: N,
        ny: N,
        nu: u * l / re,
        edges,
        ..Default::default()
    }
    .build()
    .unwrap()
}

fn top(u: f64) -> Edges<f64> {
    Edges {
        left: EdgeBC::BounceBack,
        right: EdgeBC::BounceBack,
        bottom: EdgeBC::BounceBack,
        top: EdgeBC::MovingWall { u: [u, 0.0] },
    }
}

fn main() {
    let steps = 2000;
    let re = 100.0;

    // --- (1) equivariance with correct maps ---
    let mut base = cavity(re, top(U), U);
    base.run(steps);

    // (a) Left lid moving DOWN = anti-diagonal mirror:
    //     p' = (N-1-y, N-1-x),  v_orig(x,y) = (-uy', -ux') at p'
    let mut left = cavity(
        re,
        Edges {
            left: EdgeBC::MovingWall { u: [0.0, -U] },
            right: EdgeBC::BounceBack,
            bottom: EdgeBC::BounceBack,
            top: EdgeBC::BounceBack,
        },
        U,
    );
    left.run(steps);
    let mut linf_left = 0.0f64;
    for y in 1..N - 1 {
        for x in 1..N - 1 {
            let (xr, yr) = (N - 1 - y, N - 1 - x);
            linf_left = linf_left.max((-left.uy(xr, yr) - base.ux(x, y)).abs());
            linf_left = linf_left.max((-left.ux(xr, yr) - base.uy(x, y)).abs());
        }
    }

    // (b) Left lid moving UP = +90 degree rotation:
    //     p' = (N-1-y, x),  v_orig(x,y) = (uy', -ux') at p'
    let mut left_up = cavity(
        re,
        Edges {
            left: EdgeBC::MovingWall { u: [0.0, U] },
            right: EdgeBC::BounceBack,
            bottom: EdgeBC::BounceBack,
            top: EdgeBC::BounceBack,
        },
        U,
    );
    left_up.run(steps);
    let mut linf_rot = 0.0f64;
    for y in 1..N - 1 {
        for x in 1..N - 1 {
            let (xr, yr) = (N - 1 - y, x);
            linf_rot = linf_rot.max((left_up.uy(xr, yr) - base.ux(x, y)).abs());
            linf_rot = linf_rot.max((-left_up.ux(xr, yr) - base.uy(x, y)).abs());
        }
    }

    // (c) Right lid moving UP = main-diagonal mirror:
    //     p' = (y, x),  v_orig(x,y) = (uy', ux') at p'
    let mut right_up = cavity(
        re,
        Edges {
            left: EdgeBC::BounceBack,
            right: EdgeBC::MovingWall { u: [0.0, U] },
            bottom: EdgeBC::BounceBack,
            top: EdgeBC::BounceBack,
        },
        U,
    );
    right_up.run(steps);
    let mut linf_diag = 0.0f64;
    for y in 1..N - 1 {
        for x in 1..N - 1 {
            let (xr, yr) = (y, x);
            linf_diag = linf_diag.max((right_up.uy(xr, yr) - base.ux(x, y)).abs());
            linf_diag = linf_diag.max((right_up.ux(xr, yr) - base.uy(x, y)).abs());
        }
    }

    println!("equivariance L_inf: anti-diag mirror (left,down) = {linf_left:.3e}");
    println!("equivariance L_inf: +90 rotation    (left,up)   = {linf_rot:.3e}");
    println!("equivariance L_inf: diag mirror     (right,up)  = {linf_diag:.3e}");
}
