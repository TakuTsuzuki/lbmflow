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
pub(crate) const FLAG_WALE: u32 = 32_768;

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

fn storage_load(storage: Storage, expr: &str) -> String {
    match storage {
        Storage::F32 => expr.to_string(),
        Storage::F16 => format!("f32({expr})"),
    }
}

fn storage_store(storage: Storage, expr: &str) -> String {
    match storage {
        Storage::F32 => expr.to_string(),
        Storage::F16 => format!("f16({expr})"),
    }
}

/// Shared prologue of `step` / `moments`: bounds check, solid skip,
/// population loads (`f0..`), force vector and V1-order moments.
fn emit_cell_prologue<L: Lattice>(s: &mut String, allow_cached_moments: bool, storage: Storage) {
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
            let _ = writeln!(s, "    let fq0 = {};", storage_load(storage, "f_in[i]"));
        } else {
            let _ = writeln!(
                s,
                "    let fq{q} = {};",
                storage_load(storage, &format!("f_in[{q}u * n + i]"))
            );
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

fn emit_step_entry<L: Lattice>(
    s: &mut String,
    name: &str,
    allow_cached_moments: bool,
    wale_omega: bool,
    storage: Storage,
    central_moment: bool,
) {
    let (wgx, wgy) = WG;
    let _ = writeln!(s, "@compute @workgroup_size({wgx}, {wgy}, 1)");
    let _ = writeln!(
        s,
        "fn {name}(@builtin(global_invocation_id) gid: vec3<u32>) {{"
    );
    emit_cell_prologue::<L>(s, allow_cached_moments, storage);
    if !allow_cached_moments {
        *s += "    let _keep_moment_bindings = arrayLength(&rho_out) + arrayLength(&ux_out) + arrayLength(&uy_out) + arrayLength(&uz_out);\n";
    }
    let rest = L::REST;
    if central_moment {
        emit_central_moment_collide::<L>(s);
    } else {
    // Collide (collide_row): equilibria + Guo sources per direction, then
    // TRT pair relaxation. cu/cf per q with V1's seeded-dot association.
    *s += "    let usq = ux * ux + uy * uy + uz * uz;\n";
    *s += "    let uf = ux * fvx + uy * fvy + uz * fvz;\n";
    *s += "    let drho = rho - 1.0f;\n";
    if wale_omega {
        *s += "    let op = omega_out[i];\n";
        *s += "    let cp = 1.0f - 0.5f * op;\n";
        *s += "    let om = P.omega_m;\n";
        *s += "    let cm = P.cm;\n";
    } else {
        *s += "    let op = P.omega_p;\n";
        *s += "    let om = P.omega_m;\n";
        *s += "    let cp = P.cp;\n";
        *s += "    let cm = P.cm;\n";
    }
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
    }
    // Push (stream_row's scatter dual). Rest population stays home.
    let _ = writeln!(
        s,
        "    f_out[{rest}u * n + i] = {};",
        storage_store(storage, &format!("fc{rest}"))
    );
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
        let _ = writeln!(
            s,
            "                f_out[{o}u * n + i] = {};",
            storage_store(storage, "fin")
        );
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
        let _ = writeln!(
            s,
            "                f_out[{q}u * n + j] = {};",
            storage_store(storage, &format!("fc{q}"))
        );
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
            let _ = writeln!(
                s,
                "        f_out[{u}u * n + i] = {};",
                storage_store(storage, &storage_load(storage, &format!("stash_in[sl{k}]")))
            );
            let _ = writeln!(
                s,
                "        stash_out[sl{k}] = {};",
                storage_store(storage, &format!("fc{u}"))
            );
        }
        *s += "    }\n";
        offset_terms.push(format!(
            "{}u * {} * {}",
            unk.len(),
            ext_names[t1],
            ext_names[t2]
        ));
    }
    *s += "}\n\n";
}

fn central_basis_vec<L: Lattice>() -> Vec<[u8; 3]> {
    let mut basis = Vec::with_capacity(L::Q);
    for ax in 0..=2 {
        for ay in 0..=2 {
            for az in 0..=2 {
                if L::D == 2 && az != 0 {
                    continue;
                }
                if L::D == 3 && L::Q == 19 && ax > 0 && ay > 0 && az > 0 {
                    continue;
                }
                if basis.len() < L::Q {
                    basis.push([ax, ay, az]);
                }
            }
        }
    }
    basis
}

fn phi_expr(c: [i8; 3], exp: [u8; 3]) -> String {
    let term = |axis: usize, name: &str| -> String {
        let base = format!("({}.0f - {name})", c[axis]);
        match exp[axis] {
            0 => "1.0f".to_string(),
            1 => base,
            2 => format!("{base} * {base}"),
            _ => unreachable!(),
        }
    };
    let mut out = format!("{} * {}", term(0, "ux"), term(1, "uy"));
    if exp[2] != 0 {
        out = format!("{out} * {}", term(2, "uz"));
    }
    out
}

fn emit_central_moment_collide<L: Lattice>(s: &mut String) {
    let n = L::Q;
    let basis = central_basis_vec::<L>();
    *s += "    var phys: array<f32, ";
    let _ = writeln!(s, "{n}>;");
    *s += "    var source: array<f32, ";
    let _ = writeln!(s, "{n}>;");
    *s += "    var feq_phys: array<f32, ";
    let _ = writeln!(s, "{n}>;");
    *s += "    let usq = ux * ux + uy * uy + uz * uz;\n";
    *s += "    let uf = ux * fvx + uy * fvy + uz * fvz;\n";
    for q in 0..n {
        let w = lit(L::W[q] as f32);
        let cu = dot_expr(L::C[q], ["ux", "uy", "uz"]);
        let cf = dot_expr(L::C[q], ["fvx", "fvy", "fvz"]);
        let _ = writeln!(s, "    phys[{q}] = fq{q} + {w};");
        let _ = writeln!(
            s,
            "    source[{q}] = {w} * (3.0f * (({cf}) - uf) + 9.0f * ({cu}) * ({cf}));"
        );
        let _ = writeln!(
            s,
            "    feq_phys[{q}] = {w} * rho * (1.0f + 3.0f * ({cu}) + 4.5f * ({cu}) * ({cu}) - 1.5f * usq);"
        );
    }
    let _ = writeln!(s, "    var mom: array<f32, {n}>;");
    let _ = writeln!(s, "    var src_mom: array<f32, {n}>;");
    let _ = writeln!(s, "    var eq: array<f32, {n}>;");
    for m in 0..n {
        let _ = writeln!(s, "    mom[{m}] = 0.0f;");
        let _ = writeln!(s, "    src_mom[{m}] = 0.0f;");
        let _ = writeln!(s, "    eq[{m}] = 0.0f;");
        for q in 0..n {
            let phi = phi_expr(L::C[q], basis[m]);
            let _ = writeln!(s, "    mom[{m}] += ({phi}) * phys[{q}];");
            let _ = writeln!(s, "    src_mom[{m}] += ({phi}) * source[{q}];");
            let _ = writeln!(s, "    eq[{m}] += ({phi}) * feq_phys[{q}];");
        }
    }
    let offset = if L::D == 3 && L::Q == 19 { "0.0025f" } else { "0.0f" };
    if L::D == 3 {
        *s += "    let os_base = select(P.omega_shear, omega_out[i], (P.flags & FLAG_WALE) != 0u);\n";
    } else {
        *s += "    let os_base = P.omega_shear;\n";
    }
    let _ = writeln!(
        s,
        "    let os = min(2.0f, os_base * (1.0f + {offset} - 0.16f * usq));"
    );
    let _ = writeln!(s, "    var post: array<f32, {n}>;");
    for (m, e) in basis.iter().enumerate() {
        let order = e[0] as usize + e[1] as usize + e[2] as usize;
        let rate = match order {
            0 | 1 => "0.0f",
            2 => "os",
            _ => "1.0f",
        };
        let _ = writeln!(
            s,
            "    post[{m}] = mom[{m}] - ({rate}) * (mom[{m}] - eq[{m}]) + (1.0f - 0.5f * ({rate})) * src_mom[{m}];"
        );
    }
    let diag: Vec<usize> = (0..L::D)
        .map(|a| {
            basis
                .iter()
                .position(|e| {
                    let mut de = [0u8; 3];
                    de[a] = 2;
                    *e == de
                })
                .expect("missing diagonal central moment")
        })
        .collect();
    let inv_d = lit((1.0 / L::D as f32) as f32);
    let trace_neq = diag
        .iter()
        .map(|idx| format!("(mom[{idx}] - eq[{idx}])"))
        .collect::<Vec<_>>()
        .join(" + ");
    let trace_src = diag
        .iter()
        .map(|idx| format!("src_mom[{idx}]"))
        .collect::<Vec<_>>()
        .join(" + ");
    let _ = writeln!(s, "    let bulk_neq = ({trace_neq}) * {inv_d};");
    let _ = writeln!(s, "    let bulk_src = ({trace_src}) * {inv_d};");
    for idx in diag {
        let _ = writeln!(
            s,
            "    post[{idx}] = eq[{idx}] + (1.0f - os) * (mom[{idx}] - eq[{idx}] - bulk_neq) + 0.5f * bulk_src + (1.0f - 0.5f * os) * (src_mom[{idx}] - bulk_src);"
        );
    }
    let cols = n + 1;
    let _ = writeln!(s, "    var a: array<array<f32, {cols}>, {n}>;");
    for m in 0..n {
        for q in 0..n {
            let phi = phi_expr(L::C[q], basis[m]);
            let _ = writeln!(s, "    a[{m}][{q}] = {phi};");
        }
        let _ = writeln!(s, "    a[{m}][{n}] = post[{m}];");
    }
    let _ = writeln!(s, "    for (var col: u32 = 0u; col < {n}u; col = col + 1u) {{");
    *s += "        var pivot = col;\n";
    let _ = writeln!(s, "        for (var row: u32 = col + 1u; row < {n}u; row = row + 1u) {{");
    *s += "            if (abs(a[row][col]) > abs(a[pivot][col])) { pivot = row; }\n";
    *s += "        }\n";
    *s += "        if (pivot != col) {\n";
    let _ = writeln!(s, "            for (var j: u32 = col; j <= {n}u; j = j + 1u) {{");
    *s += "                let tmp = a[col][j]; a[col][j] = a[pivot][j]; a[pivot][j] = tmp;\n";
    *s += "            }\n";
    *s += "        }\n";
    *s += "        let inv = 1.0f / a[col][col];\n";
    let _ = writeln!(s, "        for (var j: u32 = col; j <= {n}u; j = j + 1u) {{ a[col][j] = a[col][j] * inv; }}");
    let _ = writeln!(s, "        for (var row: u32 = 0u; row < {n}u; row = row + 1u) {{");
    *s += "            if (row == col) { continue; }\n";
    *s += "            let factor = a[row][col];\n";
    let _ = writeln!(s, "            for (var j: u32 = col; j <= {n}u; j = j + 1u) {{ a[row][j] = a[row][j] - factor * a[col][j]; }}");
    *s += "        }\n";
    *s += "    }\n";
    for q in 0..n {
        let w = lit(L::W[q] as f32);
        let _ = writeln!(s, "    let fc{q} = a[{q}][{n}] - {w};");
    }
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
    s += "    flags: u32,\n    omega_shear: f32,\n    collision: u32,\n";
    s += "}\n\n";
    s += "struct BcParams {\n";
    for (name, ty) in BC_PARAMS_FIELDS {
        let _ = writeln!(s, "    {name}: {ty},");
    }
    s += "}\n\n";
    s += "@group(0) @binding(0) var<uniform> P: Params;\n";
    let f_ty = if storage == Storage::F16 {
        "f16"
    } else {
        "f32"
    };
    let _ = writeln!(
        s,
        "@group(0) @binding(1) var<storage, read> f_in: array<{f_ty}>;"
    );
    let _ = writeln!(
        s,
        "@group(0) @binding(2) var<storage, read_write> f_out: array<{f_ty}>;"
    );
    s += "@group(0) @binding(3) var<storage, read> mask: array<u32>;\n";
    s += "@group(0) @binding(4) var<storage, read> wall_u: array<vec3<f32>>;\n";
    s += "@group(0) @binding(5) var<storage, read> force_field: array<vec3<f32>>;\n";
    let _ = writeln!(
        s,
        "@group(0) @binding(6) var<storage, read> stash_in: array<{f_ty}>;"
    );
    let _ = writeln!(
        s,
        "@group(0) @binding(7) var<storage, read_write> stash_out: array<{f_ty}>;"
    );
    s += "@group(0) @binding(8) var<storage, read_write> probe_acc: array<atomic<u32>, 3>;\n";
    s += "@group(0) @binding(9) var<storage, read_write> rho_out: array<f32>;\n";
    s += "@group(0) @binding(10) var<storage, read_write> ux_out: array<f32>;\n";
    s += "@group(0) @binding(11) var<storage, read_write> uy_out: array<f32>;\n";
    s += "@group(0) @binding(14) var<storage, read_write> uz_out: array<f32>;\n";
    s += "@group(0) @binding(15) var<storage, read_write> omega_out: array<f32>;\n";
    s += "@group(0) @binding(12) var<uniform> B: BcParams;\n";
    s += "@group(0) @binding(13) var<storage, read> profile: array<vec3<f32>>;\n\n";
    let _ = writeln!(s, "const FLAG_FF: u32 = {FLAG_FORCE_FIELD}u;");
    let _ = writeln!(s, "const FLAG_WALE: u32 = {FLAG_WALE}u;");
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
    emit_step_entry::<L>(&mut s, "step", false, false, storage, false);
    emit_step_entry::<L>(&mut s, "step_cached", true, false, storage, false);
    emit_step_entry::<L>(&mut s, "step_wale", false, true, storage, false);
    emit_step_entry::<L>(&mut s, "step_cached_wale", true, true, storage, false);
    emit_step_entry::<L>(&mut s, "step_cumulant", false, false, storage, true);
    emit_step_entry::<L>(&mut s, "step_cached_cumulant", true, false, storage, true);
    emit_step_entry::<L>(&mut s, "step_wale_cumulant", false, true, storage, true);
    emit_step_entry::<L>(
        &mut s,
        "step_cached_wale_cumulant",
        true,
        true,
        storage,
        true,
    );

    // ---------------------------------------------------------- moments
    let _ = writeln!(s, "@compute @workgroup_size({wgx}, {wgy}, 1)");
    s += "fn moments(@builtin(global_invocation_id) gid: vec3<u32>) {\n";
    emit_cell_prologue::<L>(&mut s, false, storage);
    s += "    rho_out[i] = rho;\n";
    s += "    ux_out[i] = ux;\n";
    s += "    uy_out[i] = uy;\n";
    s += "    uz_out[i] = uz;\n";
    s += "}\n\n";

    // ---------------------------------------------------------- wale_omega
    s += "struct CellVel { v: vec3<f32>, solid: bool, valid: bool }\n\n";
    s += "fn cell_vel(j: u32) -> CellVel {\n";
    s += "    if ((mask[j] & 1u) != 0u) { return CellVel(wall_u[j].xyz, true, true); }\n";
    s += "    return CellVel(vec3<f32>(ux_out[j], uy_out[j], uz_out[j]), false, true);\n";
    s += "}\n\n";
    s += "fn neighbor_vel(x: u32, y: u32, z: u32, axis: u32, positive: bool) -> CellVel {\n";
    s += "    var p = vec3<i32>(i32(x), i32(y), i32(z));\n";
    s += "    let dims = vec3<i32>(i32(P.nx), i32(P.ny), i32(P.nz));\n";
    s += "    let da = select(-1, 1, positive);\n";
    s += "    p[axis] = p[axis] + da;\n";
    s += "    if (p[axis] < 0 || p[axis] >= dims[axis]) {\n";
    s += "        var periodic = false;\n";
    s += "        if (axis == 0u) { periodic = ((P.flags & 1u) != 0u) || ((P.flags & 2u) != 0u); }\n";
    s += "        if (axis == 1u) { periodic = ((P.flags & 4u) != 0u) || ((P.flags & 8u) != 0u); }\n";
    s += "        if (axis == 2u) { periodic = ((P.flags & 16u) != 0u) || ((P.flags & 32u) != 0u); }\n";
    s += "        if (!periodic) { return CellVel(vec3<f32>(0.0f), false, false); }\n";
    s += "        p[axis] = (p[axis] + dims[axis]) % dims[axis];\n";
    s += "    }\n";
    s += "    let j = u32(p.z) * P.nx * P.ny + u32(p.y) * P.nx + u32(p.x);\n";
    s += "    return cell_vel(j);\n";
    s += "}\n\n";
    s += "fn fd_comp(own: f32, plus: CellVel, minus: CellVel, comp: u32) -> f32 {\n";
    s += "    if (plus.valid && !plus.solid && minus.valid && !minus.solid) { return 0.5f * (plus.v[comp] - minus.v[comp]); }\n";
    s += "    if (plus.valid && !plus.solid && minus.valid && minus.solid) { return -1.3333333333333333f * minus.v[comp] + own + 0.3333333333333333f * plus.v[comp]; }\n";
    s += "    if (plus.valid && plus.solid && minus.valid && !minus.solid) { return 1.3333333333333333f * plus.v[comp] - own - 0.3333333333333333f * minus.v[comp]; }\n";
    s += "    if (plus.valid && !plus.solid && !minus.valid) { return plus.v[comp] - own; }\n";
    s += "    if (!plus.valid && minus.valid && !minus.solid) { return own - minus.v[comp]; }\n";
    s += "    return 0.0f;\n";
    s += "}\n\n";
    let _ = writeln!(s, "@compute @workgroup_size({wgx}, {wgy}, 1)");
    s += "fn wale_omega(@builtin(global_invocation_id) gid: vec3<u32>) {\n";
    s += "    let nx = P.nx;\n";
    s += "    let ny = P.ny;\n";
    s += "    let nz = P.nz;\n";
    s += "    let x = gid.x;\n";
    s += "    let y = gid.y;\n";
    s += "    let z = gid.z;\n";
    s += "    if (x >= nx || y >= ny || z >= nz) { return; }\n";
    s += "    let xy = nx * ny;\n";
    s += "    let n = xy * nz;\n";
    s += "    let i = z * xy + y * nx + x;\n";
    s += "    let base_omega = P.omega_p;\n";
    s += "    if ((mask[i] & 1u) != 0u) { omega_out[i] = base_omega; return; }\n";
    s += "    let rho = rho_out[i];\n";
    s += "    let ux = ux_out[i];\n";
    s += "    let uy = uy_out[i];\n";
    s += "    let uz = uz_out[i];\n";
    s += "    var fvx = P.fx;\n";
    s += "    var fvy = P.fy;\n";
    s += "    var fvz = P.fz;\n";
    s += "    if ((P.flags & FLAG_FF) != 0u) {\n";
    s += "        let ffv = force_field[i];\n";
    s += "        fvx = fvx + ffv.x;\n";
    s += "        fvy = fvy + ffv.y;\n";
    s += "        fvz = fvz + ffv.z;\n";
    s += "    }\n";
    s += "    let usq = ux * ux + uy * uy + uz * uz;\n";
    s += "    let drho = rho - 1.0f;\n";
    s += "    var pxx = 0.0f;\n";
    s += "    var pyy = 0.0f;\n";
    s += "    var pzz = 0.0f;\n";
    s += "    var pxy = 0.0f;\n";
    s += "    var pxz = 0.0f;\n";
    s += "    var pyz = 0.0f;\n";
    for q in 0..L::Q {
        let c = L::C[q];
        let w = lit(L::W[q] as f32);
        let cu = dot_expr(c, ["ux", "uy", "uz"]);
        let _ = writeln!(s, "    let wq{q} = {w};");
        let _ = writeln!(s, "    let cuw{q} = {cu};");
        let _ = writeln!(
            s,
            "    let feq_w{q} = wq{q} * (drho + rho * (3.0f * cuw{q} + 4.5f * cuw{q} * cuw{q} - 1.5f * usq));"
        );
        let _ = writeln!(
            s,
            "    let fnq{q} = {} - feq_w{q};",
            storage_load(storage, &format!("f_in[{q}u * n + i]"))
        );
        if c[0] != 0 {
            let _ = writeln!(s, "    pxx = pxx + fnq{q};");
        }
        if c[1] != 0 {
            let _ = writeln!(s, "    pyy = pyy + fnq{q};");
        }
        if c[2] != 0 {
            let _ = writeln!(s, "    pzz = pzz + fnq{q};");
        }
        match c[0] * c[1] {
            1 => {
                let _ = writeln!(s, "    pxy = pxy + fnq{q};");
            }
            -1 => {
                let _ = writeln!(s, "    pxy = pxy - fnq{q};");
            }
            _ => {}
        }
        match c[0] * c[2] {
            1 => {
                let _ = writeln!(s, "    pxz = pxz + fnq{q};");
            }
            -1 => {
                let _ = writeln!(s, "    pxz = pxz - fnq{q};");
            }
            _ => {}
        }
        match c[1] * c[2] {
            1 => {
                let _ = writeln!(s, "    pyz = pyz + fnq{q};");
            }
            -1 => {
                let _ = writeln!(s, "    pyz = pyz - fnq{q};");
            }
            _ => {}
        }
    }
    s += "    pxx = pxx + ux * fvx;\n";
    s += "    pyy = pyy + uy * fvy;\n";
    s += "    pzz = pzz + uz * fvz;\n";
    s += "    pxy = pxy + 0.5f * (ux * fvy + uy * fvx);\n";
    s += "    pxz = pxz + 0.5f * (ux * fvz + uz * fvx);\n";
    s += "    pyz = pyz + 0.5f * (uy * fvz + uz * fvy);\n";
    s += "    let stress_scale = -1.5f * P.omega_p / rho;\n";
    s += "    var sxx = pxx * stress_scale;\n";
    s += "    var syy = pyy * stress_scale;\n";
    s += "    var szz = pzz * stress_scale;\n";
    s += "    var sxy = pxy * stress_scale;\n";
    s += "    var sxz = pxz * stress_scale;\n";
    s += "    var syz = pyz * stress_scale;\n";
    if L::D == 2 {
        s += "    szz = 0.0f; sxz = 0.0f; syz = 0.0f;\n";
    }
    s += "    let own = vec3<f32>(ux, uy, uz);\n";
    s += "    let xp = neighbor_vel(x, y, z, 0u, true);\n";
    s += "    let xm = neighbor_vel(x, y, z, 0u, false);\n";
    s += "    let yp = neighbor_vel(x, y, z, 1u, true);\n";
    s += "    let ym = neighbor_vel(x, y, z, 1u, false);\n";
    s += "    let zp = neighbor_vel(x, y, z, 2u, true);\n";
    s += "    let zm = neighbor_vel(x, y, z, 2u, false);\n";
    s += "    let f00 = fd_comp(own.x, xp, xm, 0u);\n";
    s += "    let f01 = fd_comp(own.x, yp, ym, 0u);\n";
    s += "    let f02 = fd_comp(own.x, zp, zm, 0u);\n";
    s += "    let f10 = fd_comp(own.y, xp, xm, 1u);\n";
    s += "    let f11 = fd_comp(own.y, yp, ym, 1u);\n";
    s += "    let f12 = fd_comp(own.y, zp, zm, 1u);\n";
    s += "    let f20 = fd_comp(own.z, xp, xm, 2u);\n";
    s += "    let f21 = fd_comp(own.z, yp, ym, 2u);\n";
    s += "    let f22 = fd_comp(own.z, zp, zm, 2u);\n";
    s += "    let g00 = f00;\n";
    s += "    let g11 = f11;\n";
    s += "    let g22 = f22;\n";
    s += "    let g01 = sxy + 0.5f * (f01 - f10);\n";
    s += "    let g10 = sxy + 0.5f * (f10 - f01);\n";
    s += "    let g02 = sxz + 0.5f * (f02 - f20);\n";
    s += "    let g20 = sxz + 0.5f * (f20 - f02);\n";
    s += "    let g12 = syz + 0.5f * (f12 - f21);\n";
    s += "    let g21 = syz + 0.5f * (f21 - f12);\n";
    s += "    let ss = g00 * g00 + g11 * g11 + g22 * g22 + 2.0f * (sxy * sxy + sxz * sxz + syz * syz);\n";
    s += "    let g200 = g00 * g00 + g01 * g10 + g02 * g20;\n";
    s += "    let g201 = g00 * g01 + g01 * g11 + g02 * g21;\n";
    s += "    let g202 = g00 * g02 + g01 * g12 + g02 * g22;\n";
    s += "    let g210 = g10 * g00 + g11 * g10 + g12 * g20;\n";
    s += "    let g211 = g10 * g01 + g11 * g11 + g12 * g21;\n";
    s += "    let g212 = g10 * g02 + g11 * g12 + g12 * g22;\n";
    s += "    let g220 = g20 * g00 + g21 * g10 + g22 * g20;\n";
    s += "    let g221 = g20 * g01 + g21 * g11 + g22 * g21;\n";
    s += "    let g222 = g20 * g02 + g21 * g12 + g22 * g22;\n";
    s += "    let tr = g200 + g211 + g222;\n";
    s += "    let sd00 = g200 - tr / 3.0f;\n";
    s += "    let sd11 = g211 - tr / 3.0f;\n";
    s += "    let sd22 = g222 - tr / 3.0f;\n";
    s += "    let sd01 = 0.5f * (g201 + g210);\n";
    s += "    let sd02 = 0.5f * (g202 + g220);\n";
    s += "    let sd12 = 0.5f * (g212 + g221);\n";
    s += "    let sdsd = sd00 * sd00 + sd11 * sd11 + sd22 * sd22 + 2.0f * (sd01 * sd01 + sd02 * sd02 + sd12 * sd12);\n";
    s += "    let denom = pow(ss, 2.5f) + pow(sdsd, 1.25f);\n";
    s += "    var nut = 0.0f;\n";
    s += "    if (denom <= 1.0e-30f) { omega_out[i] = base_omega; return; }\n";
    s += "    nut = 0.105625f * pow(sdsd, 1.5f) / denom;\n";
    s += "    let nu0 = (1.0f / base_omega - 0.5f) / 3.0f;\n";
    s += "    omega_out[i] = 1.0f / (3.0f * (nu0 + nut) + 0.5f);\n";
    s += "}\n\n";

    // --------------------------------------------- open-face moment fixup
    s += "fn fix_bc_moments(i: u32, n: u32) {\n";
    let terms: Vec<String> = (0..L::Q)
        .map(|q| storage_load(storage, &format!("f_out[{q}u * n + i]")))
        .collect();
    let mx_terms: Vec<String> = (0..L::Q)
        .filter_map(|q| match L::C[q][0] {
            1 => Some(storage_load(storage, &format!("f_out[{q}u * n + i]"))),
            -1 => Some(format!(
                "-{}",
                storage_load(storage, &format!("f_out[{q}u * n + i]"))
            )),
            0 => None,
            _ => unreachable!(),
        })
        .collect();
    let my_terms: Vec<String> = (0..L::Q)
        .filter_map(|q| match L::C[q][1] {
            1 => Some(storage_load(storage, &format!("f_out[{q}u * n + i]"))),
            -1 => Some(format!(
                "-{}",
                storage_load(storage, &format!("f_out[{q}u * n + i]"))
            )),
            0 => None,
            _ => unreachable!(),
        })
        .collect();
    let mz_terms: Vec<String> = (0..L::Q)
        .filter_map(|q| match L::C[q][2] {
            1 => Some(storage_load(storage, &format!("f_out[{q}u * n + i]"))),
            -1 => Some(format!(
                "-{}",
                storage_load(storage, &format!("f_out[{q}u * n + i]"))
            )),
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
    let _ = writeln!(
        s,
        "        let ft1 = {};",
        storage_load(storage, "f_out[B.q_t1 * n + i]")
    );
    let _ = writeln!(
        s,
        "        let fmt1 = {};",
        storage_load(storage, "f_out[B.q_mt1 * n + i]")
    );
    let _ = writeln!(
        s,
        "        let ft2 = {};",
        storage_load(storage, "f_out[B.q_t2 * n + i]")
    );
    let _ = writeln!(
        s,
        "        let fmt2 = {};",
        storage_load(storage, "f_out[B.q_mt2 * n + i]")
    );
    let _ = writeln!(
        s,
        "        let fpp = {};",
        storage_load(storage, "f_out[B.q_pp * n + i]")
    );
    let _ = writeln!(
        s,
        "        let fpm = {};",
        storage_load(storage, "f_out[B.q_pm * n + i]")
    );
    let _ = writeln!(
        s,
        "        let fmp = {};",
        storage_load(storage, "f_out[B.q_mp * n + i]")
    );
    let _ = writeln!(
        s,
        "        let fmm = {};",
        storage_load(storage, "f_out[B.q_mm * n + i]")
    );
    let _ = writeln!(
        s,
        "        var s0 = {} + ft1 + fmt1;",
        storage_load(storage, &format!("f_out[{rest}u * n + i]"))
    );
    let _ = writeln!(
        s,
        "        var sneg = {} + {} + {};",
        storage_load(storage, "f_out[B.o_n * n + i]"),
        storage_load(storage, "f_out[B.o_p1 * n + i]"),
        storage_load(storage, "f_out[B.o_m1 * n + i]")
    );
    s += "        if (B.unk_count == 5u) {\n";
    s += "            s0 = s0 + ft2 + fmt2 + fpp + fpm + fmp + fmm;\n";
    let _ = writeln!(
        s,
        "            sneg = sneg + {} + {};",
        storage_load(storage, "f_out[B.o_p2 * n + i]"),
        storage_load(storage, "f_out[B.o_m2 * n + i]")
    );
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
        "        f_out[B.q_n * n + i] = {};",
        storage_store(
            storage,
            &format!(
                "{} + {c23} * r * un",
                storage_load(storage, "f_out[B.o_n * n + i]")
            )
        )
    );
    let _ = writeln!(
        s,
        "            f_out[B.q_p1 * n + i] = {};",
        storage_store(
            storage,
            &format!(
                "{} + {c16} * r * un + tcorr",
                storage_load(storage, "f_out[B.o_p1 * n + i]")
            )
        )
    );
    let _ = writeln!(
        s,
        "            f_out[B.q_m1 * n + i] = {};",
        storage_store(
            storage,
            &format!(
                "{} + {c16} * r * un - tcorr",
                storage_load(storage, "f_out[B.o_m1 * n + i]")
            )
        )
    );
    s += "        } else {\n";
    s += "            let qt1 = ft1 - fmt1 + fpp + fpm - fmp - fmm;\n";
    s += "            let qt2 = ft2 - fmt2 + fpp - fpm + fmp - fmm;\n";
    s += "            let n1 = (1.0f / 3.0f) * r * ut1 - 0.5f * qt1;\n";
    s += "            let n2 = (1.0f / 3.0f) * r * ut2 - 0.5f * qt2;\n";
    let _ = writeln!(
        s,
        "            f_out[B.q_n * n + i] = {};",
        storage_store(
            storage,
            &format!(
                "{} + (1.0f / 3.0f) * r * un",
                storage_load(storage, "f_out[B.o_n * n + i]")
            )
        )
    );
    let _ = writeln!(
        s,
        "            f_out[B.q_p1 * n + i] = {};",
        storage_store(
            storage,
            &format!(
                "{} + (1.0f / 6.0f) * r * (un + ut1) + n1",
                storage_load(storage, "f_out[B.o_p1 * n + i]")
            )
        )
    );
    let _ = writeln!(
        s,
        "            f_out[B.q_m1 * n + i] = {};",
        storage_store(
            storage,
            &format!(
                "{} + (1.0f / 6.0f) * r * (un - ut1) - n1",
                storage_load(storage, "f_out[B.o_m1 * n + i]")
            )
        )
    );
    let _ = writeln!(
        s,
        "            f_out[B.q_p2 * n + i] = {};",
        storage_store(
            storage,
            &format!(
                "{} + (1.0f / 6.0f) * r * (un + ut2) + n2",
                storage_load(storage, "f_out[B.o_p2 * n + i]")
            )
        )
    );
    let _ = writeln!(
        s,
        "            f_out[B.q_m2 * n + i] = {};",
        storage_store(
            storage,
            &format!(
                "{} + (1.0f / 6.0f) * r * (un - ut2) - n2",
                storage_load(storage, "f_out[B.o_m2 * n + i]")
            )
        )
    );
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
    let _ = writeln!(
        s,
        "        f_out[B.unk0 * n + i] = {};",
        storage_store(storage, &storage_load(storage, "f_out[B.unk0 * n + j]"))
    );
    let _ = writeln!(
        s,
        "        f_out[B.unk1 * n + i] = {};",
        storage_store(storage, &storage_load(storage, "f_out[B.unk1 * n + j]"))
    );
    let _ = writeln!(
        s,
        "        f_out[B.unk2 * n + i] = {};",
        storage_store(storage, &storage_load(storage, "f_out[B.unk2 * n + j]"))
    );
    s += "        if (B.unk_count == 5u) {\n";
    let _ = writeln!(
        s,
        "            f_out[B.unk3 * n + i] = {};",
        storage_store(storage, &storage_load(storage, "f_out[B.unk3 * n + j]"))
    );
    let _ = writeln!(
        s,
        "            f_out[B.unk4 * n + i] = {};",
        storage_store(storage, &storage_load(storage, "f_out[B.unk4 * n + j]"))
    );
    s += "        }\n";
    s += "        fix_bc_moments(i, n);\n";
    s += "        return;\n";
    s += "    }\n";
    let _ = writeln!(s, "    if (B.kind == {BC_CONVECTIVE}u) {{");
    s += "        let lam = B.p0;\n";
    for k in 0..5 {
        let _ = writeln!(
            s,
            "        if (B.unk_count > {k}u) {{ f_out[B.unk{k} * n + i] = {}; }}",
            storage_store(
                storage,
                &format!(
                    "({} + lam * {}) * B.cinv",
                    storage_load(storage, &format!("f_out[B.unk{k} * n + i]")),
                    storage_load(storage, &format!("f_out[B.unk{k} * n + j]"))
                )
            )
        );
    }
    // Mass pinning: rho(edge) := rho(neighbour), deficit spread over the
    // unknowns by weight (convective_face, q-ascending sums).
    let di: Vec<String> = (0..L::Q)
        .map(|q| storage_load(storage, &format!("f_out[{q}u * n + i]")))
        .collect();
    let dj: Vec<String> = (0..L::Q)
        .map(|q| storage_load(storage, &format!("f_out[{q}u * n + j]")))
        .collect();
    let _ = writeln!(s, "        let di = {};", di.join(" + "));
    let _ = writeln!(s, "        let dj = {};", dj.join(" + "));
    s += "        let corr = dj - di;\n";
    for k in 0..5 {
        let _ = writeln!(
            s,
            "        if (B.unk_count > {k}u) {{ f_out[B.unk{k} * n + i] = {}; }}",
            storage_store(
                storage,
                &format!(
                    "{} + corr * B.cw{k} / B.wsum",
                    storage_load(storage, &format!("f_out[B.unk{k} * n + i]"))
                )
            )
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
    use crate::lattice::{D2Q9, D3Q19, D3Q27};

    fn validate_wgsl(source: &str, capabilities: wgpu::naga::valid::Capabilities) {
        let module = wgpu::naga::front::wgsl::parse_str(source).expect("WGSL parse failed");
        let mut validator = wgpu::naga::valid::Validator::new(
            wgpu::naga::valid::ValidationFlags::all(),
            capabilities,
        );
        validator.validate(&module).expect("WGSL validation failed");
    }

    #[test]
    fn f32_storage_generation_matches_default_byte_for_byte() {
        assert_eq!(
            generate_with_storage::<D2Q9>(Storage::F32),
            generate::<D2Q9>()
        );
        assert_eq!(
            generate_with_storage::<D3Q19>(Storage::F32),
            generate::<D3Q19>()
        );
    }

    #[test]
    fn generated_wgsl_parses_and_validates_with_naga() {
        let source = generate::<D2Q9>();
        validate_wgsl(&source, wgpu::naga::valid::Capabilities::empty());
    }

    #[test]
    fn generated_d3q19_wgsl_parses_and_validates_with_naga() {
        let source = generate::<D3Q19>();
        validate_wgsl(&source, wgpu::naga::valid::Capabilities::empty());
    }

    #[test]
    fn generated_f16_wgsl_parses_and_validates_with_naga() {
        let capabilities = wgpu::naga::valid::Capabilities::SHADER_FLOAT16;
        validate_wgsl(&generate_with_storage::<D2Q9>(Storage::F16), capabilities);
        validate_wgsl(&generate_with_storage::<D3Q19>(Storage::F16), capabilities);
        validate_wgsl(&generate_with_storage::<D3Q27>(Storage::F16), capabilities);
    }

    #[test]
    fn generated_cumulant_entries_parse_and_validate_with_naga() {
        for source in [
            generate_with_storage::<D2Q9>(Storage::F32),
            generate_with_storage::<D3Q19>(Storage::F32),
            generate_with_storage::<D3Q27>(Storage::F32),
        ] {
            assert!(source.contains("fn step_cumulant("));
            assert!(source.contains("fn step_cached_cumulant("));
            validate_wgsl(&source, wgpu::naga::valid::Capabilities::empty());
        }
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
