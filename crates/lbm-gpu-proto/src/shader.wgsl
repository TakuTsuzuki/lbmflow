// D2Q9 LBM: fused stream(pull) + collide(TRT) kernel, periodic box, f32.
//
// Population storage is *deviation form* (f_q - w_q), exactly like lbm-core:
// the quiescent state is all-zero, so f32 round-off scales with the
// fluctuation amplitude instead of the O(1) background (docs/PHYSICS.md).
//
// Direction order (lbm-core lattice.rs is the single source of truth):
//   0:(0,0) 1:(1,0) 2:(0,1) 3:(-1,0) 4:(0,-1) 5:(1,1) 6:(-1,1) 7:(-1,-1) 8:(1,-1)
// TRT pairs: (1,3) (2,4) (5,7) (6,8).
//
// Memory layout is SoA, one plane per direction: f[q*n + y*nx + x]. Threads
// in a workgroup then touch consecutive addresses per plane (coalesced), the
// key layout difference vs the CPU's cache-friendly cell-major AoS.
//
// Operator-order note: one dispatch computes f_new = Collide(Stream(f_old)),
// whereas lbm-core's step is Stream(Collide(.)). The harness uploads a
// pre-collided initial state, so after k dispatches the GPU buffer holds
// Collide(cpu_state_k). Density and momentum are collision invariants, so
// velocity/density fields compare 1:1 against the CPU (see src/main.rs).
//
// __WGX__ / __WGY__ are substituted by the host before compilation.

struct Params {
    nx: u32,
    ny: u32,
    omega_p: f32,
    omega_m: f32,
}

@group(0) @binding(0) var<uniform> P: Params;
@group(0) @binding(1) var<storage, read> f_in: array<f32>;
@group(0) @binding(2) var<storage, read_write> f_out: array<f32>;
@group(0) @binding(3) var<storage, read_write> vel: array<vec2<f32>>;

// Equilibrium in deviation form: feq_q - w_q, written in terms of
// drho = rho - 1 to avoid large-magnitude cancellation (matches lbm-core).
fn feq_dev(w: f32, drho: f32, rho: f32, cu: f32, usq15: f32) -> f32 {
    return w * (drho + rho * (3.0 * cu + 4.5 * cu * cu - usq15));
}

@compute @workgroup_size(__WGX__, __WGY__, 1)
fn step(@builtin(global_invocation_id) gid: vec3<u32>) {
    let nx = P.nx;
    let ny = P.ny;
    let x = gid.x;
    let y = gid.y;
    if (x >= nx || y >= ny) {
        return;
    }
    let n = nx * ny;
    let row = y * nx;

    // Periodic neighbour coordinates.
    let xm = select(x - 1u, nx - 1u, x == 0u);
    let xp = select(x + 1u, 0u, x + 1u == nx);
    let ym = select(y - 1u, ny - 1u, y == 0u);
    let yp = select(y + 1u, 0u, y + 1u == ny);
    let rowm = ym * nx;
    let rowp = yp * nx;

    // Pull-stream: f_q comes from the cell at (x - cx_q, y - cy_q).
    let f0 = f_in[row + x];
    let f1 = f_in[n + row + xm];
    let f2 = f_in[2u * n + rowm + x];
    let f3 = f_in[3u * n + row + xp];
    let f4 = f_in[4u * n + rowp + x];
    let f5 = f_in[5u * n + rowm + xm];
    let f6 = f_in[6u * n + rowm + xp];
    let f7 = f_in[7u * n + rowp + xp];
    let f8 = f_in[8u * n + rowp + xm];

    // Moments (deviation form: rho = 1 + sum f_dev).
    let drho = f0 + f1 + f2 + f3 + f4 + f5 + f6 + f7 + f8;
    let rho = 1.0 + drho;
    let inv = 1.0 / rho;
    let ux = (f1 - f3 + f5 - f6 - f7 + f8) * inv;
    let uy = (f2 - f4 + f5 + f6 - f7 - f8) * inv;

    // Equilibria.
    let usq15 = 1.5 * (ux * ux + uy * uy);
    let w0 = 4.0 / 9.0;
    let ws = 1.0 / 9.0;
    let wd = 1.0 / 36.0;
    let cu5 = ux + uy;
    let cu6 = uy - ux;
    let e0 = feq_dev(w0, drho, rho, 0.0, usq15);
    let e1 = feq_dev(ws, drho, rho, ux, usq15);
    let e2 = feq_dev(ws, drho, rho, uy, usq15);
    let e3 = feq_dev(ws, drho, rho, -ux, usq15);
    let e4 = feq_dev(ws, drho, rho, -uy, usq15);
    let e5 = feq_dev(wd, drho, rho, cu5, usq15);
    let e6 = feq_dev(wd, drho, rho, cu6, usq15);
    let e7 = feq_dev(wd, drho, rho, -cu5, usq15);
    let e8 = feq_dev(wd, drho, rho, -cu6, usq15);

    // TRT relaxation (BGK falls out when omega_m == omega_p).
    let op = P.omega_p;
    let om = P.omega_m;
    let i = row + x;

    f_out[i] = f0 - op * (f0 - e0);

    var fp = 0.5 * (f1 + f3);
    var fm = 0.5 * (f1 - f3);
    var relp = op * (fp - 0.5 * (e1 + e3));
    var relm = om * (fm - 0.5 * (e1 - e3));
    f_out[n + i] = f1 - relp - relm;
    f_out[3u * n + i] = f3 - relp + relm;

    fp = 0.5 * (f2 + f4);
    fm = 0.5 * (f2 - f4);
    relp = op * (fp - 0.5 * (e2 + e4));
    relm = om * (fm - 0.5 * (e2 - e4));
    f_out[2u * n + i] = f2 - relp - relm;
    f_out[4u * n + i] = f4 - relp + relm;

    fp = 0.5 * (f5 + f7);
    fm = 0.5 * (f5 - f7);
    relp = op * (fp - 0.5 * (e5 + e7));
    relm = om * (fm - 0.5 * (e5 - e7));
    f_out[5u * n + i] = f5 - relp - relm;
    f_out[7u * n + i] = f7 - relp + relm;

    fp = 0.5 * (f6 + f8);
    fm = 0.5 * (f6 - f8);
    relp = op * (fp - 0.5 * (e6 + e8));
    relm = om * (fm - 0.5 * (e6 - e8));
    f_out[6u * n + i] = f6 - relp - relm;
    f_out[8u * n + i] = f8 - relp + relm;
}

// Writes (ux, uy) per cell from the *current* populations; used for readback
// (verification and, in a production backend, GUI field extraction).
@compute @workgroup_size(__WGX__, __WGY__, 1)
fn moments(@builtin(global_invocation_id) gid: vec3<u32>) {
    let nx = P.nx;
    let ny = P.ny;
    let x = gid.x;
    let y = gid.y;
    if (x >= nx || y >= ny) {
        return;
    }
    let n = nx * ny;
    let i = y * nx + x;
    let f0 = f_in[i];
    let f1 = f_in[n + i];
    let f2 = f_in[2u * n + i];
    let f3 = f_in[3u * n + i];
    let f4 = f_in[4u * n + i];
    let f5 = f_in[5u * n + i];
    let f6 = f_in[6u * n + i];
    let f7 = f_in[7u * n + i];
    let f8 = f_in[8u * n + i];
    let rho = 1.0 + f0 + f1 + f2 + f3 + f4 + f5 + f6 + f7 + f8;
    let inv = 1.0 / rho;
    vel[i] = vec2<f32>(
        (f1 - f3 + f5 - f6 - f7 + f8) * inv,
        (f2 - f4 + f5 + f6 - f7 - f8) * inv,
    );
}
