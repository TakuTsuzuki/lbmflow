//! WGSL kernel generation from the [`Lattice`] tables.
//!
//! The shader source is *generated* (not hand-written) so the velocity set,
//! weights, opposites, TRT pairs and per-face unknown tables can never drift
//! from `lattice.rs` — the "transcribe the face tables" option of
//! docs/ARCHITECTURE_V2.md §2.4.
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
//! - `moments` — recompute (rho, ux, uy, uz) from the current populations
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

/// Distribution-buffer storage precision. Arithmetic stays f32; this controls
/// only `f`/stash buffer declarations and load/store wrappers.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Storage {
    F32,
    F16,
}

/// Flag bits of `Params.flags` (must match the generated WGSL constants).
pub(crate) const FLAG_HALO: [u32; 6] = [1, 2, 4, 8, 16, 32];
pub(crate) const FLAG_FORCE_FIELD: u32 = 64;
pub(crate) const FLAG_PROBE: u32 = 128;
pub(crate) const FLAG_OPEN_FACE: [u32; 6] = [256, 512, 1024, 2048, 4096, 8192];
pub(crate) const FLAG_CACHED_MOMENTS: u32 = 16_384;

/// BC kind codes of `BcParams.kind` (0 = inactive face).
pub(crate) const BC_VELOCITY: u32 = 1;
pub(crate) const BC_PRESSURE: u32 = 2;
pub(crate) const BC_OUTFLOW: u32 = 3;
pub(crate) const BC_CONVECTIVE: u32 = 4;

/// `BcParams` field order. Rust writes `bc_words[face][index]` in this order.
pub(crate) const BC_PARAMS_FIELDS: [(&str, &str); 64] = [
    ("kind", "u32"),
    ("base", "u32"),
    ("stride1", "u32"),
    ("extent", "u32"),
    ("joff", "i32"),
    ("has_profile", "u32"),
    ("stride2", "u32"),
    ("extent1", "u32"),
    ("q_n", "u32"),
    ("o_n", "u32"),
    ("q_p1", "u32"),
    ("o_p1", "u32"),
    ("q_m1", "u32"),
    ("o_m1", "u32"),
    ("q_p2", "u32"),
    ("o_p2", "u32"),
    ("q_m2", "u32"),
    ("o_m2", "u32"),
    ("q_t1", "u32"),
    ("q_mt1", "u32"),
    ("q_t2", "u32"),
    ("q_mt2", "u32"),
    ("q_pp", "u32"),
    ("q_pm", "u32"),
    ("q_mp", "u32"),
    ("q_mm", "u32"),
    ("unk0", "u32"),
    ("unk1", "u32"),
    ("unk2", "u32"),
    ("unk3", "u32"),
    ("unk4", "u32"),
    ("unk_count", "u32"),
    ("p0", "f32"),
    ("p1", "f32"),
    ("p2", "f32"),
    ("nxr", "f32"),
    ("nyr", "f32"),
    ("nzr", "f32"),
    ("t1x", "f32"),
    ("t1y", "f32"),
    ("t1z", "f32"),
    ("t2x", "f32"),
    ("t2y", "f32"),
    ("t2z", "f32"),
    ("cw0", "f32"),
    ("cw1", "f32"),
    ("cw2", "f32"),
    ("cw3", "f32"),
    ("cw4", "f32"),
    ("wsum", "f32"),
    ("cinv", "f32"),
    ("pad0", "u32"),
    ("pad1", "u32"),
    ("pad2", "u32"),
    ("pad3", "u32"),
    ("pad4", "u32"),
    ("pad5", "u32"),
    ("pad6", "u32"),
    ("pad7", "u32"),
    ("pad8", "u32"),
    ("pad9", "u32"),
    ("pad10", "u32"),
    ("pad11", "u32"),
    ("pad12", "u32"),
];

/// Stash slots per buffer: `sum_faces |unknowns(face)| * extent(face)`,
/// faces in `Face::index` order (the same order the offsets are generated
/// in).
pub(crate) fn stash_len<L: Lattice>(nx: usize, ny: usize, nz: usize) -> usize {
    let mut n = 0;
    for face in Face::ALL {
        if face.axis() >= L::D {
            continue;
        }
        let (t1, t2) = face.tangents();
        let ext = [nx, ny, nz][t1] * [nx, ny, nz][t2];
        n += L::unknowns(face).len() * ext;
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
fn dot_expr(coef: [i8; 3], vars: [&str; 3]) -> String {
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
/// (V1 accumulation order for `dr`, `mx`, `my`, `mz`).
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
fn emit_cell_prologue<L: Lattice>(s: &mut String, allow_cached_moments: bool) {
    *s += "    let nx = P.nx;\n";
    *s += "    let ny = P.ny;\n";
    *s += "    let nz = P.nz;\n";
    *s += "    let x = gid.x;\n";
    *s += "    let y = gid.y;\n";
    *s += "    let z = gid.z;\n";
    *s += "    if (x >= nx || y >= ny || z >= nz) { return; }\n";
    *s += "    let xy = nx * ny;\n";
    *s += "    let n = xy * nz;\n";
    *s += "    let i = z * xy + y * nx + x;\n";
    *s += "    if ((mask[i] & 1u) != 0u) { return; }\n";
    for q in 0..L::Q {
        if q == 0 {
            let _ = writeln!(s, "    let fq0 = f_in[i];");
        } else {
            let _ = writeln!(s, "    let fq{q} = f_in[{q}u * n + i];");
        }
    }
    // fv = uniform force + optional per-cell field (V1 collide_row order).
    *s += "    var fvx = P.fx;\n";
    *s += "    var fvy = P.fy;\n";
    *s += "    var fvz = P.fz;\n";
    *s += "    if ((P.flags & FLAG_FF) != 0u) {\n";
    *s += "        let ffv = force_field[i];\n";
    *s += "        fvx = fvx + ffv.x;\n";
    *s += "        fvy = fvy + ffv.y;\n";
    *s += "        fvz = fvz + ffv.z;\n";
    *s += "    }\n";
    let _ = writeln!(s, "    let dr = {};", sum_expr::<L>("fq", "", |_| 1));
    *s += "    var rho = 1.0f + dr;\n";
    *s += "    var inv = 1.0f / rho;\n";
    let _ = writeln!(
        s,
        "    let mx = {};",
        sum_expr::<L>("fq", "", |q| L::C[q][0])
    );
    let _ = writeln!(
        s,
        "    let my = {};",
        sum_expr::<L>("fq", "", |q| L::C[q][1])
    );
    let _ = writeln!(
        s,
        "    let mz = {};",
        sum_expr::<L>("fq", "", |q| L::C[q][2])
    );
    // moments_row: u = (m + f/2) / rho — this is the value collide reads.
    *s += "    var ux = (mx + 0.5f * fvx) * inv;\n";
    *s += "    var uy = (my + 0.5f * fvy) * inv;\n";
    *s += "    var uz = (mz + 0.5f * fvz) * inv;\n";
    if allow_cached_moments {
        *s += "    var use_cached_moments = (P.flags & 16384u) != 0u;\n";
        for face in Face::ALL {
            if face.axis() >= L::D {
                continue;
            }
            let cond = match face {
                Face::XNeg => "x == 0u",
                Face::XPos => "x == nx - 1u",
                Face::YNeg => "y == 0u",
                Face::YPos => "y == ny - 1u",
                Face::ZNeg => "z == 0u",
                Face::ZPos => "z == nz - 1u",
            };
            let flag = FLAG_OPEN_FACE[face.index()];
            let _ = writeln!(
                s,
                "    if ({cond} && (P.flags & {flag}u) != 0u) {{ use_cached_moments = true; }}"
            );
        }
        *s += "    if (use_cached_moments) {\n";
        *s += "        rho = rho_out[i];\n";
        *s += "        ux = ux_out[i];\n";
        *s += "        uy = uy_out[i];\n";
        *s += "        uz = uz_out[i];\n";
        *s += "        inv = 1.0f / rho;\n";
        *s += "    }\n";
    }
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

fn emit_step_entry<L: Lattice>(s: &mut String, name: &str, allow_cached_moments: bool) {
    let (wgx, wgy) = WG;
    let _ = writeln!(s, "@compute @workgroup_size({wgx}, {wgy}, 1)");
    let _ = writeln!(
        s,
        "fn {name}(@builtin(global_invocation_id) gid: vec3<u32>) {{"
    );
    emit_cell_prologue::<L>(s, allow_cached_moments);
    if !allow_cached_moments {
        *s += "    let _keep_moment_bindings = arrayLength(&rho_out) + arrayLength(&ux_out) + arrayLength(&uy_out) + arrayLength(&uz_out);\n";
    }
    // Collide (collide_row): equilibria + Guo sources per direction, then
    // TRT pair relaxation. cu/cf per q with V1's seeded-dot association.
    *s += "    let usq = ux * ux + uy * uy + uz * uz;\n";
    *s += "    let uf = ux * fvx + uy * fvy + uz * fvz;\n";
    *s += "    let drho = rho - 1.0f;\n";
    *s += "    let op = P.omega_p;\n";
    *s += "    let om = P.omega_m;\n";
    *s += "    let cp = P.cp;\n";
    *s += "    let cm = P.cm;\n";
    for q in 0..L::Q {
        let c = L::C[q];
        let w = lit(L::W[q] as f32);
        let cu = dot_expr(c, ["ux", "uy", "uz"]);
        let cf = dot_expr(c, ["fvx", "fvy", "fvz"]);
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
        "    let fc{rest} = fq{rest} - op * (fq{rest} - e{rest}) + cp * s{rest};"
    );
    for &(a, b) in L::PAIRS {
        let _ = writeln!(s, "    let fp{a} = 0.5f * (fq{a} + fq{b});");
        let _ = writeln!(s, "    let fm{a} = 0.5f * (fq{a} - fq{b});");
        let _ = writeln!(s, "    let ep{a} = 0.5f * (e{a} + e{b});");
        let _ = writeln!(s, "    let em{a} = 0.5f * (e{a} - e{b});");
        let _ = writeln!(s, "    let sp{a} = 0.5f * (s{a} + s{b});");
        let _ = writeln!(s, "    let sm{a} = 0.5f * (s{a} - s{b});");
        let _ = writeln!(s, "    let rp{a} = op * (fp{a} - ep{a});");
        let _ = writeln!(s, "    let rm{a} = om * (fm{a} - em{a});");
        let _ = writeln!(
            s,
            "    let fc{a} = fq{a} - rp{a} - rm{a} + cp * sp{a} + cm * sm{a};"
        );
        let _ = writeln!(
            s,
            "    let fc{b} = fq{b} - rp{a} + rm{a} + cp * sp{a} - cm * sm{a};"
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
        let cub = dot_expr(L::C[o], ["wu.x", "wu.y", "wu.z"]);
        let _ = writeln!(
            s,
            "    {{ // q = {q}, c = ({}, {}, {}), opp = {o}",
            c[0], c[1], c[2]
        );
        *s += "        var ok = true;\n";
        emit_push_coord(s, 0, c[0], "x", "nx");
        emit_push_coord(s, 1, c[1], "y", "ny");
        emit_push_coord(s, 2, c[2], "z", "nz");
        *s += "        if (ok) {\n";
        *s += "            let j = dz * xy + dy * nx + dx;\n";
        *s += "            let mj = mask[j];\n";
        *s += "            if ((mj & 1u) != 0u) {\n";
        *s += "                let wu = wall_u[j];\n";
        let _ = writeln!(s, "                let cub = {cub};");
        let _ = writeln!(s, "                let fin = fc{q} + {sixw} * rho * cub;");
        let _ = writeln!(s, "                f_out[{o}u * n + i] = fin;");
        *s += "                if ((mj & 2u) != 0u) {\n";
        let _ = writeln!(s, "                    let ftot = fc{q} + fin + {twow};");
        for (axis, &ca) in c.iter().enumerate() {
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
        *s += "                }\n";
        *s += "            } else {\n";
        let _ = writeln!(s, "                f_out[{q}u * n + j] = fc{q};");
        *s += "            }\n";
        *s += "        }\n";
        *s += "    }\n";
    }
    // Edge stash: carry last step's post-collide values into the skipped
    // unknown slots and stash this step's for the next (V1 stale-slot
    // mechanics; see module docs). Offsets in Face::index order.
    let mut offset_terms: Vec<String> = Vec::new();
    for face in Face::ALL {
        if face.axis() >= L::D {
            continue;
        }
        let unk = L::unknowns(face);
        let (t1, t2) = face.tangents();
        let ext_names = ["nx", "ny", "nz"];
        let (cond, tvar, ext) = match face {
            Face::XNeg => ("x == 0u", format!("z * ny + y"), "ny * nz".to_string()),
            Face::XPos => ("x == nx - 1u", format!("z * ny + y"), "ny * nz".to_string()),
            Face::YNeg => ("y == 0u", format!("z * nx + x"), "nx * nz".to_string()),
            Face::YPos => ("y == ny - 1u", format!("z * nx + x"), "nx * nz".to_string()),
            Face::ZNeg => ("z == 0u", format!("y * nx + x"), "nx * ny".to_string()),
            Face::ZPos => ("z == nz - 1u", format!("y * nx + x"), "nx * ny".to_string()),
        };
        let flag = FLAG_HALO[face.index()];
        let off = if offset_terms.is_empty() {
            "0u".to_string()
        } else {
            offset_terms.join(" + ")
        };
        let _ = writeln!(
            s,
            "    if ({cond} && (P.flags & {flag}u) == 0u) {{ // {face:?} edge stash"
        );
        for (k, &u) in unk.iter().enumerate() {
            let _ = writeln!(s, "        let sl{k} = {off} + {k}u * ({ext}) + {tvar};");
            let _ = writeln!(s, "        f_out[{u}u * n + i] = stash_in[sl{k}];");
            let _ = writeln!(s, "        stash_out[sl{k}] = fc{u};");
        }
        *s += "    }\n";
        offset_terms.push(format!("{}u * {} * {}", unk.len(), ext_names[t1], ext_names[t2]));
    }
    *s += "}\n\n";
}

/// Generate the complete shader module for lattice `L` (asserted 2D by the
/// backend). Everything below `binding(0..13)` is shared by all entry
/// points; auto pipeline layouts keep each entry point's bind group minimal.
#[cfg(test)]
pub(crate) fn generate<L: Lattice>() -> String {
    generate_with_storage::<L>(Storage::F32)
}

/// Generate the complete shader module for lattice `L`.
pub(crate) fn generate_with_storage<L: Lattice>(storage: Storage) -> String {
    let (wgx, wgy) = WG;
    let rest = L::REST;
    let mut s = String::with_capacity(32 * 1024);
    let _ = writeln!(
        s,
        "// Generated by lbm_core::gpu::wgsl from the Lattice tables (D{}Q{}).",
        L::D,
        L::Q
    );
    if storage == Storage::F16 {
        s += "enable f16;\n";
    }
    s += "// Deviation-form SoA planes f[q*n + z*(nx*ny) + y*nx + x]; push-fused collide+stream.\n\n";
    s += "struct Params {\n";
    s += "    nx: u32,\n    ny: u32,\n    nz: u32,\n    pad_dim: u32,\n";
    s += "    omega_p: f32,\n    omega_m: f32,\n";
    s += "    cp: f32,\n    cm: f32,\n";
    s += "    fx: f32,\n    fy: f32,\n    fz: f32,\n";
    s += "    flags: u32,\n    pad0: u32,\n    pad1: u32,\n";
    s += "}\n\n";
    s += "struct BcParams {\n";
    for (name, ty) in BC_PARAMS_FIELDS {
        let _ = writeln!(s, "    {name}: {ty},");
    }
    s += "}\n\n";
    s += "@group(0) @binding(0) var<uniform> P: Params;\n";
    let f_ty = if storage == Storage::F16 { "f16" } else { "f32" };
    let _ = writeln!(s, "@group(0) @binding(1) var<storage, read> f_in: array<{f_ty}>;");
    let _ = writeln!(s, "@group(0) @binding(2) var<storage, read_write> f_out: array<{f_ty}>;");
    s += "@group(0) @binding(3) var<storage, read> mask: array<u32>;\n";
    s += "@group(0) @binding(4) var<storage, read> wall_u: array<vec3<f32>>;\n";
    s += "@group(0) @binding(5) var<storage, read> force_field: array<vec3<f32>>;\n";
    let _ = writeln!(s, "@group(0) @binding(6) var<storage, read> stash_in: array<{f_ty}>;");
    let _ = writeln!(s, "@group(0) @binding(7) var<storage, read_write> stash_out: array<{f_ty}>;");
    s += "@group(0) @binding(8) var<storage, read_write> probe_acc: array<atomic<u32>, 3>;\n";
    s += "@group(0) @binding(9) var<storage, read_write> rho_out: array<f32>;\n";
    s += "@group(0) @binding(10) var<storage, read_write> ux_out: array<f32>;\n";
    s += "@group(0) @binding(11) var<storage, read_write> uy_out: array<f32>;\n";
    s += "@group(0) @binding(14) var<storage, read_write> uz_out: array<f32>;\n";
    s += "@group(0) @binding(12) var<uniform> B: BcParams;\n";
    s += "@group(0) @binding(13) var<storage, read> profile: array<vec3<f32>>;\n\n";
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
    emit_step_entry::<L>(&mut s, "step", false);
    emit_step_entry::<L>(&mut s, "step_cached", true);

    // ---------------------------------------------------------- moments
    let _ = writeln!(s, "@compute @workgroup_size({wgx}, {wgy}, 1)");
    s += "fn moments(@builtin(global_invocation_id) gid: vec3<u32>) {\n";
    emit_cell_prologue::<L>(&mut s, false);
    s += "    rho_out[i] = rho;\n";
    s += "    ux_out[i] = ux;\n";
    s += "    uy_out[i] = uy;\n";
    s += "    uz_out[i] = uz;\n";
    s += "}\n\n";

    // --------------------------------------------- open-face moment fixup
    s += "fn fix_bc_moments(i: u32, n: u32) {\n";
    let terms: Vec<String> = (0..L::Q).map(|q| format!("f_out[{q}u * n + i]")).collect();
    let mx_terms: Vec<String> = (0..L::Q)
        .filter_map(|q| match L::C[q][0] {
            1 => Some(format!("f_out[{q}u * n + i]")),
            -1 => Some(format!("-f_out[{q}u * n + i]")),
            0 => None,
            _ => unreachable!(),
        })
        .collect();
    let my_terms: Vec<String> = (0..L::Q)
        .filter_map(|q| match L::C[q][1] {
            1 => Some(format!("f_out[{q}u * n + i]")),
            -1 => Some(format!("-f_out[{q}u * n + i]")),
            0 => None,
            _ => unreachable!(),
        })
        .collect();
    let mz_terms: Vec<String> = (0..L::Q)
        .filter_map(|q| match L::C[q][2] {
            1 => Some(format!("f_out[{q}u * n + i]")),
            -1 => Some(format!("-f_out[{q}u * n + i]")),
            0 => None,
            _ => unreachable!(),
        })
        .collect();
    let _ = writeln!(s, "    let dr = {};", terms.join(" + "));
    let _ = writeln!(
        s,
        "    let mx = {};",
        if mx_terms.is_empty() {
            "0.0f".to_string()
        } else {
            mx_terms.join(" + ")
        }
    );
    let _ = writeln!(
        s,
        "    let my = {};",
        if my_terms.is_empty() {
            "0.0f".to_string()
        } else {
            my_terms.join(" + ")
        }
    );
    let _ = writeln!(
        s,
        "    let mz = {};",
        if mz_terms.is_empty() {
            "0.0f".to_string()
        } else {
            mz_terms.join(" + ")
        }
    );
    s += "    var fvx = P.fx;\n";
    s += "    var fvy = P.fy;\n";
    s += "    var fvz = P.fz;\n";
    s += "    if ((P.flags & FLAG_FF) != 0u) {\n";
    s += "        let ffv = force_field[i];\n";
    s += "        fvx = fvx + ffv.x;\n";
    s += "        fvy = fvy + ffv.y;\n";
    s += "        fvz = fvz + ffv.z;\n";
    s += "    }\n";
    s += "    let r = 1.0f + dr;\n";
    s += "    let inv = 1.0f / r;\n";
    s += "    rho_out[i] = r;\n";
    s += "    ux_out[i] = (mx + 0.5f * fvx) * inv;\n";
    s += "    uy_out[i] = (my + 0.5f * fvy) * inv;\n";
    s += "    uz_out[i] = (mz + 0.5f * fvz) * inv;\n";
    s += "}\n\n";

    // --------------------------------------------------------------- bc
    let c23 = lit((2.0f64 / 3.0) as f32);
    let c16 = lit((1.0f64 / 6.0) as f32);
    let _ = writeln!(s, "@compute @workgroup_size({WG_BC}, 1, 1)");
    s += "fn bc(@builtin(global_invocation_id) gid: vec3<u32>) {\n";
    s += "    let t = gid.x;\n";
    s += "    if (t >= B.extent) { return; }\n";
    s += "    let n = P.nx * P.ny * P.nz;\n";
    s += "    let c1 = t % B.extent1;\n";
    s += "    let c2 = t / B.extent1;\n";
    s += "    let i = B.base + c1 * B.stride1 + c2 * B.stride2;\n";
    s += "    if ((mask[i] & 1u) != 0u) { return; }\n";
    let _ = writeln!(
        s,
        "    if (B.kind == {BC_VELOCITY}u || B.kind == {BC_PRESSURE}u) {{"
    );
    s += "        let ft1 = f_out[B.q_t1 * n + i];\n";
    s += "        let fmt1 = f_out[B.q_mt1 * n + i];\n";
    s += "        let ft2 = f_out[B.q_t2 * n + i];\n";
    s += "        let fmt2 = f_out[B.q_mt2 * n + i];\n";
    s += "        let fpp = f_out[B.q_pp * n + i];\n";
    s += "        let fpm = f_out[B.q_pm * n + i];\n";
    s += "        let fmp = f_out[B.q_mp * n + i];\n";
    s += "        let fmm = f_out[B.q_mm * n + i];\n";
    let _ = writeln!(s, "        var s0 = f_out[{rest}u * n + i] + ft1 + fmt1;");
    s += "        var sneg = f_out[B.o_n * n + i] + f_out[B.o_p1 * n + i] + f_out[B.o_m1 * n + i];\n";
    s += "        if (B.unk_count == 5u) {\n";
    s += "            s0 = s0 + ft2 + fmt2 + fpp + fpm + fmp + fmm;\n";
    s += "            sneg = sneg + f_out[B.o_p2 * n + i] + f_out[B.o_m2 * n + i];\n";
    s += "        }\n";
    s += "        let closure = s0 + 2.0f * sneg + 1.0f;\n";
    s += "        var r = 0.0f;\n";
    s += "        var un = 0.0f;\n";
    s += "        var ut1 = 0.0f;\n";
    s += "        var ut2 = 0.0f;\n";
    let _ = writeln!(s, "        if (B.kind == {BC_VELOCITY}u) {{");
    s += "            var ubx = B.p0;\n";
    s += "            var uby = B.p1;\n";
    s += "            var ubz = B.p2;\n";
    s += "            if (B.has_profile != 0u) {\n";
    s += "                let pr = profile[t];\n";
    s += "                ubx = pr.x;\n";
    s += "                uby = pr.y;\n";
    s += "                ubz = pr.z;\n";
    s += "            }\n";
    s += "            un = ubx * B.nxr + uby * B.nyr + ubz * B.nzr;\n";
    s += "            ut1 = ubx * B.t1x + uby * B.t1y + ubz * B.t1z;\n";
    s += "            ut2 = ubx * B.t2x + uby * B.t2y + ubz * B.t2z;\n";
    s += "            r = closure / (1.0f - un);\n";
    s += "        } else {\n";
    s += "            r = B.p0;\n";
    s += "            un = 1.0f - closure / r;\n";
    s += "            ut1 = 0.0f;\n";
    s += "            ut2 = 0.0f;\n";
    s += "        }\n";
    s += "        if (B.unk_count == 3u) {\n";
    s += "            let tcorr = 0.5f * (r * ut1 - (ft1 - fmt1));\n";
    let _ = writeln!(
        s,
        "        f_out[B.q_n * n + i] = f_out[B.o_n * n + i] + {c23} * r * un;"
    );
    let _ = writeln!(
        s,
        "            f_out[B.q_p1 * n + i] = f_out[B.o_p1 * n + i] + {c16} * r * un + tcorr;"
    );
    let _ = writeln!(
        s,
        "            f_out[B.q_m1 * n + i] = f_out[B.o_m1 * n + i] + {c16} * r * un - tcorr;"
    );
    s += "        } else {\n";
    s += "            let qt1 = ft1 - fmt1 + fpp + fpm - fmp - fmm;\n";
    s += "            let qt2 = ft2 - fmt2 + fpp - fpm + fmp - fmm;\n";
    s += "            let n1 = (1.0f / 3.0f) * r * ut1 - 0.5f * qt1;\n";
    s += "            let n2 = (1.0f / 3.0f) * r * ut2 - 0.5f * qt2;\n";
    s += "            f_out[B.q_n * n + i] = f_out[B.o_n * n + i] + (1.0f / 3.0f) * r * un;\n";
    s += "            f_out[B.q_p1 * n + i] = f_out[B.o_p1 * n + i] + (1.0f / 6.0f) * r * (un + ut1) + n1;\n";
    s += "            f_out[B.q_m1 * n + i] = f_out[B.o_m1 * n + i] + (1.0f / 6.0f) * r * (un - ut1) - n1;\n";
    s += "            f_out[B.q_p2 * n + i] = f_out[B.o_p2 * n + i] + (1.0f / 6.0f) * r * (un + ut2) + n2;\n";
    s += "            f_out[B.q_m2 * n + i] = f_out[B.o_m2 * n + i] + (1.0f / 6.0f) * r * (un - ut2) - n2;\n";
    s += "        }\n";
    s += "        fix_bc_moments(i, n);\n";
    s += "        return;\n";
    s += "    }\n";
    s += "    let j = u32(i32(i) + B.joff);\n";
    s += "    if ((mask[j] & 1u) != 0u) {\n";
    s += "        fix_bc_moments(i, n);\n";
    s += "        return;\n";
    s += "    }\n";
    let _ = writeln!(s, "    if (B.kind == {BC_OUTFLOW}u) {{");
    s += "        f_out[B.unk0 * n + i] = f_out[B.unk0 * n + j];\n";
    s += "        f_out[B.unk1 * n + i] = f_out[B.unk1 * n + j];\n";
    s += "        f_out[B.unk2 * n + i] = f_out[B.unk2 * n + j];\n";
    s += "        if (B.unk_count == 5u) {\n";
    s += "            f_out[B.unk3 * n + i] = f_out[B.unk3 * n + j];\n";
    s += "            f_out[B.unk4 * n + i] = f_out[B.unk4 * n + j];\n";
    s += "        }\n";
    s += "        fix_bc_moments(i, n);\n";
    s += "        return;\n";
    s += "    }\n";
    let _ = writeln!(s, "    if (B.kind == {BC_CONVECTIVE}u) {{");
    s += "        let lam = B.p0;\n";
    for k in 0..5 {
        let _ = writeln!(
            s,
            "        if (B.unk_count > {k}u) {{ f_out[B.unk{k} * n + i] = (f_out[B.unk{k} * n + i] + lam * f_out[B.unk{k} * n + j]) * B.cinv; }}"
        );
    }
    // Mass pinning: rho(edge) := rho(neighbour), deficit spread over the
    // unknowns by weight (convective_face, q-ascending sums).
    let di: Vec<String> = (0..L::Q).map(|q| format!("f_out[{q}u * n + i]")).collect();
    let dj: Vec<String> = (0..L::Q).map(|q| format!("f_out[{q}u * n + j]")).collect();
    let _ = writeln!(s, "        let di = {};", di.join(" + "));
    let _ = writeln!(s, "        let dj = {};", dj.join(" + "));
    s += "        let corr = dj - di;\n";
    for k in 0..5 {
        let _ = writeln!(
            s,
            "        if (B.unk_count > {k}u) {{ f_out[B.unk{k} * n + i] = f_out[B.unk{k} * n + i] + corr * B.cw{k} / B.wsum; }}"
        );
    }
    s += "        fix_bc_moments(i, n);\n";
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
    use crate::lattice::{D2Q9, D3Q19};

    #[test]
    fn generated_wgsl_parses_and_validates_with_naga() {
        let source = generate::<D2Q9>();
        let module = wgpu::naga::front::wgsl::parse_str(&source).expect("WGSL parse failed");
        let mut validator = wgpu::naga::valid::Validator::new(
            wgpu::naga::valid::ValidationFlags::all(),
            wgpu::naga::valid::Capabilities::empty(),
        );
        validator.validate(&module).expect("WGSL validation failed");
    }

    #[test]
    fn generated_d3q19_wgsl_parses_and_validates_with_naga() {
        let source = generate::<D3Q19>();
        let module = wgpu::naga::front::wgsl::parse_str(&source).expect("WGSL parse failed");
        let mut validator = wgpu::naga::valid::Validator::new(
            wgpu::naga::valid::ValidationFlags::all(),
            wgpu::naga::valid::Capabilities::empty(),
        );
        validator.validate(&module).expect("WGSL validation failed");
    }

    #[test]
    fn bc_params_field_table_matches_rust_word_indices() {
        let expected = [
            "kind",
            "base",
            "stride1",
            "extent",
            "joff",
            "has_profile",
            "stride2",
            "extent1",
            "q_n",
            "o_n",
            "q_p1",
            "o_p1",
            "q_m1",
            "o_m1",
            "q_p2",
            "o_p2",
            "q_m2",
            "o_m2",
            "q_t1",
            "q_mt1",
            "q_t2",
            "q_mt2",
            "q_pp",
            "q_pm",
            "q_mp",
            "q_mm",
            "unk0",
            "unk1",
            "unk2",
            "unk3",
            "unk4",
            "unk_count",
            "p0",
            "p1",
            "p2",
            "nxr",
            "nyr",
            "nzr",
            "t1x",
            "t1y",
            "t1z",
            "t2x",
            "t2y",
            "t2z",
            "cw0",
            "cw1",
            "cw2",
            "cw3",
            "cw4",
            "wsum",
            "cinv",
            "pad0",
            "pad1",
            "pad2",
            "pad3",
            "pad4",
            "pad5",
            "pad6",
            "pad7",
            "pad8",
            "pad9",
            "pad10",
            "pad11",
            "pad12",
        ];
        let actual: Vec<&str> = BC_PARAMS_FIELDS.iter().map(|(name, _)| *name).collect();
        assert_eq!(actual, expected);
    }

    #[test]
    fn d2q9_source_contains_expected_pieces() {
        let src = generate::<D2Q9>();
        for needle in [
            "fn step(",
            "fn moments(",
            "fn bc(",
            "fn clear_probe(",
            "atomicCompareExchangeWeak",
            "0.11111111f",  // 1/9 weight
            "0.027777778f", // 1/36 weight
        ] {
            assert!(src.contains(needle), "missing {needle:?} in generated WGSL");
        }
        // All four faces get stash blocks; exact offsets are generated from
        // the lattice face tables and covered by `stash_len_counts_all_face_slots`.
        assert!(src.contains("XNeg edge stash"));
    }

    #[test]
    fn stash_len_counts_all_face_slots() {
        assert_eq!(stash_len::<D2Q9>(10, 7, 1), 3 * 7 * 2 + 3 * 10 * 2);
    }

    #[test]
    fn literals_roundtrip() {
        for v in [
            4.5f32,
            1.0 / 9.0,
            1.0 / 36.0,
            6.0 * (1.0f32 / 9.0),
            2.0 / 3.0,
        ] {
            let l = lit(v);
            let parsed: f32 = l.trim_end_matches('f').parse().unwrap();
            assert_eq!(parsed, v, "literal {l} does not roundtrip");
        }
    }
}
