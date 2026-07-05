//! WGSL kernel generation from the [`Lattice`] tables.
//!
//! The shader source is *generated* (not hand-written) so the velocity set,
//! weights, opposites, TRT pairs and per-face unknown tables can never drift
//! from `lattice.rs` — the "transcribe the face tables" option of
//! docs/ARCHITECTURE_V2.md §2.4. The generation is table-driven and therefore
//! 3D-ready in shape; the cell indexing and face set emitted here are 2D
//! (D2Q9), which `WgpuBackend::new` asserts.
//!
//! Four entry points:
//!
//! - `step` — fused collide(TRT/BGK + Guo) + push-streaming with half-way
//!   bounce-back (still/moving walls), probe-force accumulation and the
//!   open-face "edge stash" (see below). One dispatch advances one time step
//!   with **exactly the CPU operator order** `S∘C`: the thread collides its
//!   own cell (recomputing the moments from `f_in` precisely as
//!   `moments_row` cached them) and pushes the post-collide populations to
//!   its neighbours — the scatter dual of `stream_row`'s gather, value-equal
//!   per link. Open-face BCs then run as small post-passes, exactly like the
//!   CPU's `apply_open_faces` — so no `(C∘S)^k ∘ C = C ∘ (S∘C)^k` recipe is
//!   needed (unlike the pull-fused lbm-gpu-proto kernel, which cannot host a
//!   BC pass between S and C).
//! - `bc` — one face's open-boundary pass (Zou–He velocity/pressure,
//!   outflow, convective), lattice indices delivered via a uniform.
//! - `moments` — recompute (rho, ux, uy) from the current populations
//!   (explicit-readback path; fluid cells only, matching `moments_row`).
//! - `clear_probe` — zero the probe-force accumulator (dispatched at the
//!   start of each step that has a probe, inside the same compute pass).
//!
//! **Edge stash** (V1 implicit spec): `stream_row` leaves the out-buffer
//! untouched for populations whose source lies beyond a no-halo face, so
//! after the swap those slots still hold the *previous step's post-collide*
//! values (the in-place CPU collide wrote them); `convective_face` reads
//! them as its `prev`. A push kernel writes the out-buffer, never the
//! in-buffer, so that stale content would be two steps old instead. The
//! stash restores V1 mechanics exactly: each face cell writes its own
//! post-collide unknown-slot values to `stash_out` and copies `stash_in`
//! (last step's values) into the skipped `f_out` slots. Cost: `O(perimeter)`
//! storage and traffic.
//!
//! Arithmetic fidelity: every generated expression keeps the operand order
//! and grouping of `kernels.rs` (which is V1's), with zero-coefficient terms
//! elided — IEEE-754 identities only (`0*x` additions and `x - y ≡ x + (-y)`
//! reassociations that are exact). Remaining CPU↔GPU drift comes from the
//! Metal compiler's FMA/reassociation, bounded by T14's tolerance.

use crate::lattice::{Face, Lattice};
use std::fmt::Write;

/// Workgroup size of the full-grid kernels (`step`, `moments`).
/// 256×1 was the fastest shape in the lbm-gpu-proto sweep (SoA row-direction
/// coalescing) and respects the WebGPU default 256-invocation limit.
pub(crate) const WG: (u32, u32) = (256, 1);
/// Workgroup size of the 1D face kernels.
pub(crate) const WG_BC: u32 = 64;

/// Flag bits of `Params.flags` (must match the generated WGSL constants).
pub(crate) const FLAG_HALO: [u32; 4] = [1, 2, 4, 8];
pub(crate) const FLAG_FORCE_FIELD: u32 = 16;
pub(crate) const FLAG_PROBE: u32 = 32;

/// BC kind codes of `BcParams.kind` (0 = inactive face).
pub(crate) const BC_VELOCITY: u32 = 1;
pub(crate) const BC_PRESSURE: u32 = 2;
pub(crate) const BC_OUTFLOW: u32 = 3;
pub(crate) const BC_CONVECTIVE: u32 = 4;

/// Stash slots per buffer: `sum_faces |unknowns(face)| * extent(face)`,
/// faces in `Face::index` order (the same order the offsets are generated
/// in). 2D only.
pub(crate) fn stash_len<L: Lattice>(nx: usize, ny: usize) -> usize {
    let mut n = 0;
    for face in &Face::ALL[..4] {
        let ext = if face.axis() == 0 { ny } else { nx };
        n += L::unknowns(*face).len() * ext;
    }
    n.max(1)
}

/// An f32 literal that WGSL parses back to exactly `v` (`f` suffix = typed
/// f32 literal, no double-rounding through AbstractFloat).
fn lit(v: f32) -> String {
    let mut s = format!("{v:?}");
    if !s.contains('.') && !s.contains('e') && !s.contains("inf") && !s.contains("NaN") {
        s.push_str(".0");
    }
    assert!(v.is_finite(), "non-finite shader constant {v}");
    s.push('f');
    s
}

/// `Σ coef[a] * var[a]` with V1's left-to-right association and elided zero
/// terms (`±1` coefficients fold into sign — exact in IEEE-754).
fn dot_expr(coef: [i8; 2], vars: [&str; 2]) -> String {
    let mut s = String::new();
    for (c, v) in coef.iter().zip(vars) {
        match (*c, s.is_empty()) {
            (0, _) => {}
            (1, true) => s.push_str(v),
            (-1, true) => {
                s.push('-');
                s.push_str(v);
            }
            (1, false) => {
                let _ = write!(s, " + {v}");
            }
            (-1, false) => {
                let _ = write!(s, " - {v}");
            }
            (c, _) => unreachable!("non-unit lattice velocity component {c}"),
        }
    }
    if s.is_empty() {
        s.push_str("0.0f");
    }
    s
}

/// `Σ_q coef(q) * name{q}` over ascending q with elided zero terms
/// (V1 accumulation order for `dr`, `mx`, `my`).
fn sum_expr<L: Lattice>(prefix: &str, suffix: &str, coef: impl Fn(usize) -> i8) -> String {
    let mut s = String::new();
    for q in 0..L::Q {
        match (coef(q), s.is_empty()) {
            (0, _) => {}
            (1, true) => {
                let _ = write!(s, "{prefix}{q}{suffix}");
            }
            (-1, true) => {
                let _ = write!(s, "-{prefix}{q}{suffix}");
            }
            (1, false) => {
                let _ = write!(s, " + {prefix}{q}{suffix}");
            }
            (-1, false) => {
                let _ = write!(s, " - {prefix}{q}{suffix}");
            }
            _ => unreachable!(),
        }
    }
    if s.is_empty() {
        s.push_str("0.0f");
    }
    s
}

/// Shared prologue of `step` / `moments`: bounds check, solid skip,
/// population loads (`f0..`), force vector and V1-order moments.
fn emit_cell_prologue<L: Lattice>(s: &mut String) {
    *s += "    let nx = P.nx;\n";
    *s += "    let ny = P.ny;\n";
    *s += "    let x = gid.x;\n";
    *s += "    let y = gid.y;\n";
    *s += "    if (x >= nx || y >= ny) { return; }\n";
    *s += "    let n = nx * ny;\n";
    *s += "    let i = y * nx + x;\n";
    *s += "    if ((mask[i] & 1u) != 0u) { return; }\n";
    for q in 0..L::Q {
        if q == 0 {
            let _ = writeln!(s, "    let f0 = f_in[i];");
        } else {
            let _ = writeln!(s, "    let f{q} = f_in[{q}u * n + i];");
        }
    }
    // fv = uniform force + optional per-cell field (V1 collide_row order).
    *s += "    var fvx = P.fx;\n";
    *s += "    var fvy = P.fy;\n";
    *s += "    if ((P.flags & FLAG_FF) != 0u) {\n";
    *s += "        let ffv = force_field[i];\n";
    *s += "        fvx = fvx + ffv.x;\n";
    *s += "        fvy = fvy + ffv.y;\n";
    *s += "    }\n";
    let _ = writeln!(s, "    let dr = {};", sum_expr::<L>("f", "", |_| 1));
    *s += "    let rho = 1.0f + dr;\n";
    *s += "    let inv = 1.0f / rho;\n";
    let _ = writeln!(s, "    let mx = {};", sum_expr::<L>("f", "", |q| L::C[q][0]));
    let _ = writeln!(s, "    let my = {};", sum_expr::<L>("f", "", |q| L::C[q][1]));
    // moments_row: u = (m + f/2) / rho — this is the value collide reads.
    *s += "    let ux = (mx + 0.5f * fvx) * inv;\n";
    *s += "    let uy = (my + 0.5f * fvy) * inv;\n";
}

/// Wrap-or-skip destination coordinate for one axis of a push (the scatter
/// dual of `stream_row`'s halo-or-skip source logic).
fn emit_push_coord(s: &mut String, axis: usize, c: i8, var: &str, extent: &str) {
    let (flag_neg, flag_pos) = (FLAG_HALO[2 * axis], FLAG_HALO[2 * axis + 1]);
    match c {
        0 => {
            let _ = writeln!(s, "        let d{var} = {var};");
        }
        1 => {
            let _ = writeln!(s, "        var d{var} = {var} + 1u;");
            let _ = writeln!(s, "        if (d{var} == {extent}) {{");
            let _ = writeln!(
                s,
                "            if ((P.flags & {flag_pos}u) != 0u) {{ d{var} = 0u; }} else {{ ok = false; }}"
            );
            let _ = writeln!(s, "        }}");
        }
        -1 => {
            let _ = writeln!(s, "        var d{var} = {var};");
            let _ = writeln!(s, "        if ({var} == 0u) {{");
            let _ = writeln!(
                s,
                "            if ((P.flags & {flag_neg}u) != 0u) {{ d{var} = {extent} - 1u; }} else {{ ok = false; }}"
            );
            let _ = writeln!(s, "        }} else {{ d{var} = {var} - 1u; }}");
        }
        _ => unreachable!(),
    }
}

/// Generate the complete shader module for lattice `L` (asserted 2D by the
/// backend). Everything below `binding(0..13)` is shared by all entry
/// points; auto pipeline layouts keep each entry point's bind group minimal.
pub(crate) fn generate<L: Lattice>() -> String {
    let (wgx, wgy) = WG;
    let mut s = String::with_capacity(32 * 1024);
    let _ = writeln!(
        s,
        "// Generated by lbm_core::gpu::wgsl from the Lattice tables (D{}Q{}).",
        L::D,
        L::Q
    );
    s += "// Deviation-form f32 SoA planes f[q*n + y*nx + x]; push-fused collide+stream.\n\n";
    s += "struct Params {\n";
    s += "    nx: u32,\n    ny: u32,\n";
    s += "    omega_p: f32,\n    omega_m: f32,\n";
    s += "    cp: f32,\n    cm: f32,\n";
    s += "    fx: f32,\n    fy: f32,\n";
    s += "    flags: u32,\n    pad0: u32,\n    pad1: u32,\n    pad2: u32,\n";
    s += "}\n\n";
    s += "struct BcParams {\n";
    s += "    kind: u32,\n    base: u32,\n    stride: u32,\n    extent: u32,\n";
    s += "    joff: i32,\n    has_profile: u32,\n";
    s += "    q_n: u32,\n    o_n: u32,\n    q_d1: u32,\n    o_d1: u32,\n";
    s += "    q_d2: u32,\n    o_d2: u32,\n    q_t: u32,\n    q_mt: u32,\n";
    s += "    unk0: u32,\n    unk1: u32,\n    unk2: u32,\n";
    s += "    p0: f32,\n    p1: f32,\n";
    s += "    nxr: f32,\n    nyr: f32,\n    txr: f32,\n    tyr: f32,\n";
    s += "    cw0: f32,\n    cw1: f32,\n    cw2: f32,\n    wsum: f32,\n    cinv: f32,\n";
    s += "    pad0: u32,\n    pad1: u32,\n    pad2: u32,\n    pad3: u32,\n";
    s += "}\n\n";
    s += "@group(0) @binding(0) var<uniform> P: Params;\n";
    s += "@group(0) @binding(1) var<storage, read> f_in: array<f32>;\n";
    s += "@group(0) @binding(2) var<storage, read_write> f_out: array<f32>;\n";
    s += "@group(0) @binding(3) var<storage, read> mask: array<u32>;\n";
    s += "@group(0) @binding(4) var<storage, read> wall_u: array<vec2<f32>>;\n";
    s += "@group(0) @binding(5) var<storage, read> force_field: array<vec2<f32>>;\n";
    s += "@group(0) @binding(6) var<storage, read> stash_in: array<f32>;\n";
    s += "@group(0) @binding(7) var<storage, read_write> stash_out: array<f32>;\n";
    s += "@group(0) @binding(8) var<storage, read_write> probe_acc: array<atomic<u32>, 3>;\n";
    s += "@group(0) @binding(9) var<storage, read_write> rho_out: array<f32>;\n";
    s += "@group(0) @binding(10) var<storage, read_write> ux_out: array<f32>;\n";
    s += "@group(0) @binding(11) var<storage, read_write> uy_out: array<f32>;\n";
    s += "@group(0) @binding(12) var<uniform> B: BcParams;\n";
    s += "@group(0) @binding(13) var<storage, read> profile: array<vec2<f32>>;\n\n";
    let _ = writeln!(s, "const FLAG_FF: u32 = {FLAG_FORCE_FIELD}u;");
    s += "\n";
    // f32 atomic add via compare-exchange (WGSL has no float atomics). The
    // accumulation order is nondeterministic; T14 compares the probe force
    // with a tolerance, like every diagnostic.
    s += "fn atomic_add_f32(idx: u32, v: f32) {\n";
    s += "    var old = atomicLoad(&probe_acc[idx]);\n";
    s += "    loop {\n";
    s += "        let nv = bitcast<u32>(bitcast<f32>(old) + v);\n";
    s += "        let r = atomicCompareExchangeWeak(&probe_acc[idx], old, nv);\n";
    s += "        if (r.exchanged) { break; }\n";
    s += "        old = r.old_value;\n";
    s += "    }\n";
    s += "}\n\n";

    // ------------------------------------------------------------- step
    let _ = writeln!(s, "@compute @workgroup_size({wgx}, {wgy}, 1)");
    s += "fn step(@builtin(global_invocation_id) gid: vec3<u32>) {\n";
    emit_cell_prologue::<L>(&mut s);
    // Collide (collide_row): equilibria + Guo sources per direction, then
    // TRT pair relaxation. cu/cf per q with V1's seeded-dot association.
    s += "    let usq = ux * ux + uy * uy;\n";
    s += "    let uf = ux * fvx + uy * fvy;\n";
    s += "    let drho = rho - 1.0f;\n";
    s += "    let op = P.omega_p;\n";
    s += "    let om = P.omega_m;\n";
    s += "    let cp = P.cp;\n";
    s += "    let cm = P.cm;\n";
    for q in 0..L::Q {
        let c = [L::C[q][0], L::C[q][1]];
        let w = lit(L::W[q] as f32);
        let cu = dot_expr(c, ["ux", "uy"]);
        let cf = dot_expr(c, ["fvx", "fvy"]);
        let _ = writeln!(s, "    let cu{q} = {cu};");
        let _ = writeln!(s, "    let cf{q} = {cf};");
        let _ = writeln!(
            s,
            "    let e{q} = {w} * (drho + rho * (3.0f * cu{q} + 4.5f * cu{q} * cu{q} - 1.5f * usq));"
        );
        let _ = writeln!(
            s,
            "    let s{q} = {w} * (3.0f * (cf{q} - uf) + 9.0f * cu{q} * cf{q});"
        );
    }
    let rest = L::REST;
    let _ = writeln!(
        s,
        "    let fc{rest} = f{rest} - op * (f{rest} - e{rest}) + cp * s{rest};"
    );
    for &(a, b) in L::PAIRS {
        let _ = writeln!(s, "    let fp{a} = 0.5f * (f{a} + f{b});");
        let _ = writeln!(s, "    let fm{a} = 0.5f * (f{a} - f{b});");
        let _ = writeln!(s, "    let ep{a} = 0.5f * (e{a} + e{b});");
        let _ = writeln!(s, "    let em{a} = 0.5f * (e{a} - e{b});");
        let _ = writeln!(s, "    let sp{a} = 0.5f * (s{a} + s{b});");
        let _ = writeln!(s, "    let sm{a} = 0.5f * (s{a} - s{b});");
        let _ = writeln!(s, "    let rp{a} = op * (fp{a} - ep{a});");
        let _ = writeln!(s, "    let rm{a} = om * (fm{a} - em{a});");
        let _ = writeln!(
            s,
            "    let fc{a} = f{a} - rp{a} - rm{a} + cp * sp{a} + cm * sm{a};"
        );
        let _ = writeln!(
            s,
            "    let fc{b} = f{b} - rp{a} + rm{a} + cp * sp{a} - cm * sm{a};"
        );
    }
    // Push (stream_row's scatter dual). Rest population stays home.
    let _ = writeln!(s, "    f_out[{rest}u * n + i] = fc{rest};");
    for q in 0..L::Q {
        if q == rest {
            continue;
        }
        let c = L::C[q];
        let o = L::OPP[q];
        let w32 = L::W[q] as f32;
        let sixw = lit(6.0f32 * w32);
        let twow = lit(2.0f32 * w32);
        // c_opp · wall_u — the pull-form bounce-back term of the receiving
        // (this) cell's opposite direction.
        let cub = dot_expr([L::C[o][0], L::C[o][1]], ["wu.x", "wu.y"]);
        let _ = writeln!(s, "    {{ // q = {q}, c = ({}, {}), opp = {o}", c[0], c[1]);
        s += "        var ok = true;\n";
        emit_push_coord(&mut s, 0, c[0], "x", "nx");
        emit_push_coord(&mut s, 1, c[1], "y", "ny");
        s += "        if (ok) {\n";
        s += "            let j = dy * nx + dx;\n";
        s += "            let mj = mask[j];\n";
        s += "            if ((mj & 1u) != 0u) {\n";
        s += "                let wu = wall_u[j];\n";
        let _ = writeln!(s, "                let cub = {cub};");
        let _ = writeln!(s, "                let fin = fc{q} + {sixw} * rho * cub;");
        let _ = writeln!(s, "                f_out[{o}u * n + i] = fin;");
        s += "                if ((mj & 2u) != 0u) {\n";
        let _ = writeln!(s, "                    let ftot = fc{q} + fin + {twow};");
        for (axis, &ca) in c.iter().take(2).enumerate() {
            match ca {
                1 => {
                    let _ = writeln!(s, "                    atomic_add_f32({axis}u, ftot);");
                }
                -1 => {
                    let _ = writeln!(s, "                    atomic_add_f32({axis}u, -ftot);");
                }
                _ => {}
            }
        }
        s += "                }\n";
        s += "            } else {\n";
        let _ = writeln!(s, "                f_out[{q}u * n + j] = fc{q};");
        s += "            }\n";
        s += "        }\n";
        s += "    }\n";
    }
    // Edge stash: carry last step's post-collide values into the skipped
    // unknown slots and stash this step's for the next (V1 stale-slot
    // mechanics; see module docs). Offsets in Face::index order.
    let mut off_ny = 0usize; // accumulated multiples of ny
    let mut off_nx = 0usize; // accumulated multiples of nx
    for face in &Face::ALL[..4] {
        let unk = L::unknowns(*face);
        let (cond, tvar, ext) = match face {
            Face::XNeg => ("x == 0u", "y", "ny"),
            Face::XPos => ("x == nx - 1u", "y", "ny"),
            Face::YNeg => ("y == 0u", "x", "nx"),
            Face::YPos => ("y == ny - 1u", "x", "nx"),
            _ => unreachable!(),
        };
        let flag = FLAG_HALO[face.index()];
        let off = match (off_ny, off_nx) {
            (0, 0) => "0u".to_string(),
            (a, 0) => format!("{a}u * ny"),
            (a, b) => format!("{a}u * ny + {b}u * nx"),
        };
        let _ = writeln!(
            s,
            "    if ({cond} && (P.flags & {flag}u) == 0u) {{ // {face:?} edge stash"
        );
        for (k, &u) in unk.iter().enumerate() {
            let _ = writeln!(s, "        let sl{k} = {off} + {k}u * {ext} + {tvar};");
            let _ = writeln!(s, "        f_out[{u}u * n + i] = stash_in[sl{k}];");
            let _ = writeln!(s, "        stash_out[sl{k}] = fc{u};");
        }
        s += "    }\n";
        if face.axis() == 0 {
            off_ny += unk.len();
        } else {
            off_nx += unk.len();
        }
    }
    s += "}\n\n";

    // ---------------------------------------------------------- moments
    let _ = writeln!(s, "@compute @workgroup_size({wgx}, {wgy}, 1)");
    s += "fn moments(@builtin(global_invocation_id) gid: vec3<u32>) {\n";
    emit_cell_prologue::<L>(&mut s);
    s += "    rho_out[i] = rho;\n";
    s += "    ux_out[i] = ux;\n";
    s += "    uy_out[i] = uy;\n";
    s += "}\n\n";

    // --------------------------------------------------------------- bc
    let c23 = lit((2.0f64 / 3.0) as f32);
    let c16 = lit((1.0f64 / 6.0) as f32);
    let _ = writeln!(s, "@compute @workgroup_size({WG_BC}, 1, 1)");
    s += "fn bc(@builtin(global_invocation_id) gid: vec3<u32>) {\n";
    s += "    let t = gid.x;\n";
    s += "    if (t >= B.extent) { return; }\n";
    s += "    let n = P.nx * P.ny;\n";
    s += "    let i = B.base + t * B.stride;\n";
    s += "    if ((mask[i] & 1u) != 0u) { return; }\n";
    let _ = writeln!(s, "    if (B.kind == {BC_VELOCITY}u || B.kind == {BC_PRESSURE}u) {{");
    s += "        let ft = f_out[B.q_t * n + i];\n";
    s += "        let fmt = f_out[B.q_mt * n + i];\n";
    let _ = writeln!(s, "        let s0 = f_out[{rest}u * n + i] + ft + fmt;");
    s += "        let sneg = f_out[B.o_n * n + i] + f_out[B.o_d1 * n + i] + f_out[B.o_d2 * n + i];\n";
    s += "        let closure = s0 + 2.0f * sneg + 1.0f;\n";
    s += "        var r = 0.0f;\n";
    s += "        var un = 0.0f;\n";
    s += "        var ut = 0.0f;\n";
    let _ = writeln!(s, "        if (B.kind == {BC_VELOCITY}u) {{");
    s += "            var ubx = B.p0;\n";
    s += "            var uby = B.p1;\n";
    s += "            if (B.has_profile != 0u) {\n";
    s += "                let pr = profile[t];\n";
    s += "                ubx = pr.x;\n";
    s += "                uby = pr.y;\n";
    s += "            }\n";
    s += "            un = ubx * B.nxr + uby * B.nyr;\n";
    s += "            ut = ubx * B.txr + uby * B.tyr;\n";
    s += "            r = closure / (1.0f - un);\n";
    s += "        } else {\n";
    s += "            r = B.p0;\n";
    s += "            un = 1.0f - closure / r;\n";
    s += "            ut = 0.0f;\n";
    s += "        }\n";
    s += "        let tcorr = 0.5f * (r * ut - (ft - fmt));\n";
    let _ = writeln!(s, "        f_out[B.q_n * n + i] = f_out[B.o_n * n + i] + {c23} * r * un;");
    let _ = writeln!(
        s,
        "        f_out[B.q_d1 * n + i] = f_out[B.o_d1 * n + i] + {c16} * r * un + tcorr;"
    );
    let _ = writeln!(
        s,
        "        f_out[B.q_d2 * n + i] = f_out[B.o_d2 * n + i] + {c16} * r * un - tcorr;"
    );
    s += "        return;\n";
    s += "    }\n";
    s += "    let j = u32(i32(i) + B.joff);\n";
    s += "    if ((mask[j] & 1u) != 0u) { return; }\n";
    let _ = writeln!(s, "    if (B.kind == {BC_OUTFLOW}u) {{");
    s += "        f_out[B.unk0 * n + i] = f_out[B.unk0 * n + j];\n";
    s += "        f_out[B.unk1 * n + i] = f_out[B.unk1 * n + j];\n";
    s += "        f_out[B.unk2 * n + i] = f_out[B.unk2 * n + j];\n";
    s += "        return;\n";
    s += "    }\n";
    let _ = writeln!(s, "    if (B.kind == {BC_CONVECTIVE}u) {{");
    s += "        let lam = B.p0;\n";
    for k in 0..3 {
        let _ = writeln!(
            s,
            "        f_out[B.unk{k} * n + i] = (f_out[B.unk{k} * n + i] + lam * f_out[B.unk{k} * n + j]) * B.cinv;"
        );
    }
    // Mass pinning: rho(edge) := rho(neighbour), deficit spread over the
    // unknowns by weight (convective_face, q-ascending sums).
    let di: Vec<String> = (0..L::Q)
        .map(|q| format!("f_out[{q}u * n + i]"))
        .collect();
    let dj: Vec<String> = (0..L::Q)
        .map(|q| format!("f_out[{q}u * n + j]"))
        .collect();
    let _ = writeln!(s, "        let di = {};", di.join(" + "));
    let _ = writeln!(s, "        let dj = {};", dj.join(" + "));
    s += "        let corr = dj - di;\n";
    for k in 0..3 {
        let _ = writeln!(
            s,
            "        f_out[B.unk{k} * n + i] = f_out[B.unk{k} * n + i] + corr * B.cw{k} / B.wsum;"
        );
    }
    s += "    }\n";
    s += "}\n\n";

    // ------------------------------------------------------ clear_probe
    s += "@compute @workgroup_size(1, 1, 1)\n";
    s += "fn clear_probe() {\n";
    s += "    atomicStore(&probe_acc[0], 0u);\n";
    s += "    atomicStore(&probe_acc[1], 0u);\n";
    s += "    atomicStore(&probe_acc[2], 0u);\n";
    s += "}\n";
    s
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::lattice::D2Q9;

    #[test]
    fn d2q9_source_contains_expected_pieces() {
        let src = generate::<D2Q9>();
        for needle in [
            "fn step(",
            "fn moments(",
            "fn bc(",
            "fn clear_probe(",
            "atomicCompareExchangeWeak",
            "0.11111111f", // 1/9 weight
            "0.027777778f", // 1/36 weight
        ] {
            assert!(src.contains(needle), "missing {needle:?} in generated WGSL");
        }
        // All four faces get stash blocks; offsets: 0, 3ny, 6ny, 6ny+3nx.
        assert!(src.contains("XNeg edge stash"));
        assert!(src.contains("6u * ny + 3u * nx"));
    }

    #[test]
    fn stash_len_counts_all_face_slots() {
        assert_eq!(stash_len::<D2Q9>(10, 7), 3 * 7 * 2 + 3 * 10 * 2);
    }

    #[test]
    fn literals_roundtrip() {
        for v in [4.5f32, 1.0 / 9.0, 1.0 / 36.0, 6.0 * (1.0f32 / 9.0), 2.0 / 3.0] {
            let l = lit(v);
            let parsed: f32 = l.trim_end_matches('f').parse().unwrap();
            assert_eq!(parsed, v, "literal {l} does not roundtrip");
        }
    }
}
