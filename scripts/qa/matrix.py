"""Physics-QA sweep matrix, pass 1.

Each entry is one runnable (scenario x config) with the anomaly checks that
apply to it. Reference values come ONLY from docs/VALIDATION.md frozen bands
or analytic solutions; the `source` strings cite the exact section.

Scenario JSONs are generated verbatim from the `scenario` dicts (camelCase,
matching crates/lbm-scenario). Check names map to functions in qa_checks.py;
`args` are passed through.

Geometry notes replicated from the validation tests (do not "fix" these):
- Cavity Ghia sampling: N=129, L=N-2=127, centerline interpolation
  pos = 0.5 + frac*L over fluid rows 1..=N-2 (validation_cavity.rs).
- T8 cylinder D=40 cases use an EXCLUSIVE staircase (d^2 < r^2). The scenario
  circle is inclusive (d^2 <= r^2), so r=19.995 reproduces the exclusive
  staircase of r=20 exactly (integer lattice: no point has 399.6 < d^2 < 400).
- Body force in scenarios is a constant force DENSITY (VALIDATION T6:
  momentum grows N_fluid*F per step), not rho-proportional gravity. A buoyant
  bubble is therefore NOT expressible (uniform force density is exactly
  balanced by a hydrostatic pressure gradient) -- filed as a coverage gap.
"""

# Ghia, Ghia & Shin (1982) centerline tables, copied from
# crates/lbm-core/tests/validation_cavity.rs (single source of the frozen
# values used by T7). Re=400 v(x=0.9063) is a known typo -> excluded there.
GHIA_Y = [1.0000, 0.9766, 0.9688, 0.9609, 0.9531, 0.8516, 0.7344, 0.6172,
          0.5000, 0.4531, 0.2813, 0.1719, 0.1016, 0.0703, 0.0625, 0.0547, 0.0000]
GHIA_X = [1.0000, 0.9688, 0.9609, 0.9531, 0.9453, 0.9063, 0.8594, 0.8047,
          0.5000, 0.2344, 0.2266, 0.1563, 0.0938, 0.0781, 0.0703, 0.0625, 0.0000]
U_RE100 = [1.00000, 0.84123, 0.78871, 0.73722, 0.68717, 0.23151, 0.00332,
           -0.13641, -0.20581, -0.21090, -0.15662, -0.10150, -0.06434,
           -0.04775, -0.04192, -0.03717, 0.00000]
V_RE100 = [0.00000, -0.05906, -0.07391, -0.08864, -0.10313, -0.16914,
           -0.22445, -0.24533, 0.05454, 0.17527, 0.17507, 0.16077, 0.12317,
           0.10890, 0.10091, 0.09233, 0.00000]
U_RE400 = [1.00000, 0.75837, 0.68439, 0.61756, 0.55892, 0.29093, 0.16256,
           0.02135, -0.11477, -0.17119, -0.32726, -0.24299, -0.14612,
           -0.10338, -0.09266, -0.08186, 0.00000]
V_RE400 = [0.00000, -0.12146, -0.15663, -0.19254, -0.22847, -0.23827,
           -0.44993, -0.38598, 0.05186, 0.30174, 0.30203, 0.28124, 0.22965,
           0.20920, 0.19713, 0.18360, 0.00000]
U_RE1000 = [1.00000, 0.65928, 0.57492, 0.51117, 0.46604, 0.33304, 0.18719,
            0.05702, -0.06080, -0.10648, -0.27805, -0.38289, -0.29730,
            -0.22220, -0.20196, -0.18109, 0.00000]
V_RE1000 = [0.00000, -0.21388, -0.27669, -0.33714, -0.39188, -0.51550,
            -0.42665, -0.31966, 0.02526, 0.32235, 0.33075, 0.37095, 0.32627,
            0.30353, 0.29012, 0.27485, 0.00000]

GHIA = {
    100: (U_RE100, V_RE100, 0.02, False),
    400: (U_RE400, V_RE400, 0.02, True),   # True = exclude v(0.9063) typo
    1000: (U_RE1000, V_RE1000, 0.03, False),
}

BB = {"type": "bounceBack"}
PER = {"type": "periodic"}


def mw(ux, uy=0.0):
    return {"type": "movingWall", "u": [ux, uy]}


def _sc(name, nx, ny, nu, edges, steps, *, nz=1, collision="trt",
        force=(0.0, 0.0), steady=None, init=None, obstacles=None,
        inlet_profile=None, multiphase=None, probes=None, outputs=None):
    sc = {
        "version": 0,
        "name": name,
        "grid": {"nx": nx, "ny": ny},
        "physics": {
            "nu": nu,
            "collision": {"type": collision},
            "force": list(force),
            "precision": "f64",
        },
        "edges": edges,
        "init": init or {"kind": "rest"},
        "run": {"steps": steps},
    }
    if nz > 1:
        sc["grid"]["nz"] = nz
    if steady is not None:
        sc["run"]["stopWhenSteady"] = {"epsilon": steady[0], "checkEvery": steady[1]}
    if obstacles:
        sc["obstacles"] = obstacles
    if inlet_profile:
        sc["inletProfile"] = inlet_profile
    if multiphase:
        sc["multiphase"] = multiphase
    if probes:
        sc["probes"] = probes
    if outputs:
        sc["outputs"] = outputs
    return sc


def _cavity_ghia(re):
    n, u = 129, 0.1
    nu = u * (n - 2) / re
    return _sc(
        f"qa-cavity-re{re}", n, n, nu,
        {"left": BB, "right": BB, "bottom": BB, "top": mw(u)},
        300_000, steady=(1e-8, 1000),
        outputs=[
            {"field": "ux", "format": "csv", "every": 0},
            {"field": "uy", "format": "csv", "every": 0},
            {"field": "rho", "format": "csv", "every": 10_000},
            {"field": "speed", "format": "png", "every": 0},
        ],
    )


def _cylinder_2d2(name, right_edge, steps):
    return _sc(
        name, 880, 164, 0.04,
        {"left": {"type": "velocityInlet", "u": [0.15, 0.0]},
         "right": right_edge, "bottom": BB, "top": BB},
        steps,
        inlet_profile={"edge": "left", "kind": "parabolic", "umax": 0.15},
        obstacles=[{"shape": "circle", "cx": 80.0, "cy": 80.0, "r": 19.995}],
        probes=[{"type": "force", "every": 10}],
        outputs=[
            {"field": "vorticity", "format": "png", "every": 50_000},
            {"field": "speed", "format": "png", "every": 0},
            {"field": "ux", "format": "csv", "every": 20_000},
            {"field": "rho", "format": "csv", "every": 20_000},
        ],
    )


CONFIGS = [
    # ---------------------------------------------------------- conservation
    {
        "id": "t6-momentum-periodic",
        "track": "conservation",
        "scenario": _sc(
            "qa-t6-momentum", 64, 64, 0.02,
            {"left": PER, "right": PER, "bottom": PER, "top": PER},
            2000, force=(1e-6, 0.0),
            outputs=[{"field": "ux", "format": "csv", "every": 500},
                     {"field": "rho", "format": "csv", "every": 500}],
        ),
        "checks": [
            {"name": "momentum_growth", "args": {"fx": 1e-6, "band": 1e-10},
             "source": "VALIDATION T6 (momentum grows N_fluid*F per step, rel 1e-10)"},
            {"name": "field_uniform", "args": {"field": "ux", "band": 1e-12},
             "source": "analytic (uniform force on periodic box keeps u uniform)"},
            {"name": "mass_drift", "args": {"region": "all", "band_per_1e4": 1e-11},
             "source": "VALIDATION T6 (rel 1e-11 / 1e4 steps)"},
        ],
    },
    # ------------------------------------------------------------- Poiseuille
    {
        "id": "poiseuille-trt-h32",
        "track": "channel-analytic",
        "scenario": _sc(
            "qa-poiseuille-trt", 16, 34, 0.1,
            {"left": PER, "right": PER, "bottom": BB, "top": BB},
            200_000, force=(1e-6, 0.0), steady=(1e-11, 500),
            outputs=[{"field": "ux", "format": "csv", "every": 0},
                     {"field": "rho", "format": "csv", "every": 5000},
                     {"field": "rho", "format": "csv", "every": 0}],
        ),
        "checks": [
            {"name": "poiseuille_exact",
             "args": {"g": 1e-6, "nu": 0.1, "band": 1e-10},
             "source": "VALIDATION T2 (TRT Lambda=3/16 exact, LInf_rel <= 1e-10)"},
            {"name": "profile_symmetry", "args": {"band": 1e-13},
             "source": "VALIDATION T2 (top/bottom symmetry <= 1e-13)"},
            {"name": "mass_drift", "args": {"region": "interior", "band_per_1e4": 1e-11},
             "source": "VALIDATION T6"},
        ],
    },
    {
        "id": "poiseuille-bgk-h8",
        "track": "channel-analytic",
        "scenario": _sc(
            "qa-poiseuille-bgk8", 12, 10, 0.1,
            {"left": PER, "right": PER, "bottom": BB, "top": BB},
            100_000, force=(1e-6, 0.0), collision="bgk", steady=(1e-11, 500),
            outputs=[{"field": "ux", "format": "csv", "every": 0}],
        ),
        "checks": [
            {"name": "poiseuille_l2", "args": {"g": 1e-6, "nu": 0.1},
             "source": "VALIDATION T2 (BGK: order input, paired with h16)"},
        ],
    },
    {
        "id": "poiseuille-bgk-h16",
        "track": "channel-analytic",
        "scenario": _sc(
            "qa-poiseuille-bgk16", 12, 18, 0.1,
            {"left": PER, "right": PER, "bottom": BB, "top": BB},
            150_000, force=(1e-6, 0.0), collision="bgk", steady=(1e-11, 500),
            outputs=[{"field": "ux", "format": "csv", "every": 0}],
        ),
        "checks": [
            {"name": "poiseuille_l2", "args": {"g": 1e-6, "nu": 0.1},
             "source": "VALIDATION T2 (BGK: order input, paired with h8)"},
        ],
        # cross-config: order = log2(err_h8 / err_h16) >= 1.7 (run_sweep)
        "pair_order_with": "poiseuille-bgk-h8",
    },
    # ---------------------------------------------------------------- Couette
    {
        "id": "couette-trt-tau08",
        "track": "channel-analytic",
        "scenario": _sc(
            "qa-couette-trt", 16, 34, 0.1,
            {"left": PER, "right": PER, "bottom": BB, "top": mw(0.1)},
            100_000, steady=(1e-11, 500),
            outputs=[{"field": "ux", "format": "csv", "every": 0},
                     {"field": "rho", "format": "csv", "every": 10_000}],
        ),
        "checks": [
            {"name": "couette_exact", "args": {"u_wall": 0.1, "band": 1e-10},
             "source": "VALIDATION T3 (LInf_rel <= 1e-10)"},
            {"name": "mass_drift", "args": {"region": "interior", "band_per_1e4": 1e-12},
             "source": "VALIDATION T3 (moving-wall mass drift <= 1e-12 rel / 1e4 steps)"},
        ],
    },
    {
        "id": "couette-bgk-tau06",
        "track": "channel-analytic",
        "scenario": _sc(
            "qa-couette-bgk", 16, 34, (0.6 - 0.5) / 3.0,
            {"left": PER, "right": PER, "bottom": BB, "top": mw(0.1)},
            150_000, collision="bgk", steady=(1e-11, 500),
            outputs=[{"field": "ux", "format": "csv", "every": 0}],
        ),
        "checks": [
            {"name": "couette_exact", "args": {"u_wall": 0.1, "band": 1e-10},
             "source": "VALIDATION T3 (holds for BGK too, tau in {0.6,1.0,1.4})"},
        ],
    },
    # -------------------------------------------------------- Zou-He channel
    {
        "id": "channel-zouhe-t4",
        "track": "open-boundary",
        "scenario": _sc(
            # frozen T4 params (validation_open_bc.rs): nu=0.02, eps=1e-11.
            # The 2e-3 profile band is calibrated THERE; nu=0.05 measured
            # 2.75e-3 (O(Ma^2) axial-gradient growth) -> S3 note in the log.
            "qa-channel-t4", 96, 34, 0.02,
            {"left": {"type": "velocityInlet", "u": [0.05, 0.0]},
             "right": {"type": "pressureOutlet", "rho": 1.0},
             "bottom": BB, "top": BB},
            160_000, steady=(1e-11, 500),
            inlet_profile={"edge": "left", "kind": "parabolic", "umax": 0.05},
            outputs=[{"field": "ux", "format": "csv", "every": 10_000},
                     {"field": "rho", "format": "csv", "every": 10_000},
                     {"field": "ux", "format": "csv", "every": 0},
                     {"field": "rho", "format": "csv", "every": 0}],
        ),
        "checks": [
            {"name": "t4_flow_rate", "args": {"band": 1e-4, "outlet_margin": 24},
             "source": "VALIDATION T4 (bulk mass-flux constancy <= 1e-4; measured 2.4e-5)"},
            {"name": "t4_profile", "args": {"umax": 0.05, "band": 2e-3},
             "source": "VALIDATION T4 (central profile L2rel <= 2e-3 vs parabola)"},
        ],
    },
    # ------------------------------------------------------------ cavity/Ghia
    {"id": "cavity-re100", "track": "cavity",
     "scenario": _cavity_ghia(100),
     "checks": [
         {"name": "ghia_rms", "args": {"re": 100},
          "source": "VALIDATION T7 (RMS <= 0.02*U, Ghia 1982 17 pts)"},
         {"name": "mass_drift", "args": {"region": "interior", "band_per_1e4": 1e-11},
          "source": "VALIDATION T6"},
         {"name": "speed_scale", "args": {"u_ref": 0.1, "factor": 1.2},
          "source": "sanity (|u| should not exceed lid speed by >20%)"},
     ]},
    {"id": "cavity-re400", "track": "cavity",
     "scenario": _cavity_ghia(400),
     "checks": [
         {"name": "ghia_rms", "args": {"re": 400},
          "source": "VALIDATION T7 (RMS <= 0.02*U, v(0.9063) typo excluded)"},
         {"name": "mass_drift", "args": {"region": "interior", "band_per_1e4": 1e-11},
          "source": "VALIDATION T6"},
     ]},
    {"id": "cavity-re1000", "track": "cavity",
     "scenario": _cavity_ghia(1000),
     "expect_warnings": ["physics.nu"],  # tau = 0.538 < 0.55 advisory
     "checks": [
         {"name": "ghia_rms", "args": {"re": 1000},
          "source": "VALIDATION T7 (RMS <= 0.03*U at Re=1000)"},
         {"name": "mass_drift", "args": {"region": "interior", "band_per_1e4": 1e-11},
          "source": "VALIDATION T6"},
     ]},
    # ----------------------------------------------------------- cylinder T8
    {
        "id": "cylinder-re20-t8",
        "track": "cylinder",
        "scenario": _sc(
            "qa-cylinder-re20", 440, 82, 0.05,
            {"left": {"type": "velocityInlet", "u": [0.075, 0.0]},
             "right": {"type": "pressureOutlet", "rho": 1.0},
             "bottom": BB, "top": BB},
            30_000,
            inlet_profile={"edge": "left", "kind": "parabolic", "umax": 0.075},
            obstacles=[{"shape": "circle", "cx": 40.0, "cy": 40.0, "r": 10.0}],
            probes=[{"type": "force", "every": 10}],
            outputs=[{"field": "speed", "format": "png", "every": 0}],
        ),
        "checks": [
            {"name": "cd_cl_steady",
             "args": {"d": 20.0, "u_mean": 0.05, "sample_start": 20_000,
                      "cd_band": [5.2, 6.0], "cl_band": [-0.05, 0.08]},
             "source": "VALIDATION T8 2D-1 D=20 (staircase coarse-grid band)"},
        ],
    },
    {
        "id": "cylinder-re100-karman-t8",
        "track": "cylinder",
        "scenario": _cylinder_2d2("qa-karman-re100",
                                  {"type": "pressureOutlet", "rho": 1.0},
                                  150_000),
        "checks": [
            {"name": "karman",
             "args": {"d": 40.0, "u_mean": 0.1, "window_start": 110_000,
                      "st_band": [0.28, 0.32], "cdmax_band": [3.0, 3.5],
                      "clmax_band": [0.8, 1.2], "period_var": 0.02},
             "source": "VALIDATION T8 2D-2 D=40 (St/Cd_max/Cl_max/periodicity)"},
        ],
        "note": "spec test runs 120k with a seeded perturbation; from rest we run "
                "150k and measure on the last 40k with a saturation guard",
    },
    {
        "id": "outflow-karman-t9",
        "track": "open-boundary",
        "scenario": _cylinder_2d2("qa-outflow-t9", {"type": "outflow"}, 60_000),
        "checks": [
            {"name": "reverse_flow", "args": {"band": 0.05},
             "source": "VALIDATION T9 (reverse-flow mass flux <= 5% of inflow)"},
        ],
        "note": "T9 spec horizon is 1e5 steps; pass-1 runs 60k (time cap), "
                "full horizon deferred to pass 2",
    },
    # ------------------------------------------------------------------- 3D
    {
        "id": "duct3d-t15",
        "track": "3d",
        "scenario": _sc(
            "qa-duct3d", 12, 34, 34, {  # placeholder replaced below
            }, 0),
        "checks": [],  # filled below (needs nz plumbing)
    },
    {
        "id": "cavity3d-re1000-n64",
        "track": "3d",
        "scenario": _sc(
            "qa-cavity3d", 64, 64, 0.1 * 62 / 1000,
            {"left": BB, "right": BB, "bottom": BB, "top": mw(0.1),
             "front": BB, "back": BB},
            40_000, nz=64, steady=(1e-8, 1000),
            outputs=[{"field": "ux", "format": "vtk", "every": 0},
                     {"field": "uy", "format": "vtk", "every": 0},
                     {"field": "speed", "format": "vtk", "every": 0},
                     {"field": "rho", "format": "vtk", "every": 10_000}],
        ),
        "expect_warnings": ["physics.nu", "physics"],  # tau=0.5186, U/nu=16.1
        "checks": [
            {"name": "mass_drift", "args": {"region": "interior3d", "band_per_1e4": 1e-11},
             "source": "VALIDATION T15.5 (N=64 sentinel: mass drift ~1e-16)"},
            {"name": "z_mirror_symmetry", "args": {"field": "ux", "band": 1e-8, "u_ref": 0.1},
             "source": "VALIDATION T15.5 (midplane symmetry ~2e-15 at N=64)"},
        ],
        "note": "qualitative sentinel per T15.5 default suite; spec-grade profile "
                "RMS needs N=72 + full steady (pass 2 / heavy)",
    },
    # ------------------------------------------------------------ multiphase
    {
        "id": "droplet-t11",
        "track": "multiphase",
        "scenario": _sc(
            "qa-droplet-t11", 128, 128, 1.0 / 6.0,
            {"left": PER, "right": PER, "bottom": PER, "top": PER},
            30_000,
            init={"kind": "droplet", "cx": 64.0, "cy": 64.0, "r": 20.0,
                  "rhoLiquid": 2.0, "rhoVapor": 0.15},
            multiphase={"g": -5.0, "gWall": 0.0},
            outputs=[{"field": "rho", "format": "csv", "every": 5000},
                     {"field": "rho", "format": "csv", "every": 0},
                     {"field": "speed", "format": "csv", "every": 0},
                     {"field": "rho", "format": "png", "every": 0}],
        ),
        "checks": [
            {"name": "spurious_current", "args": {"band": 5e-3},
             "source": "VALIDATION T11 (max|u| <= 5e-3; measured 1.26e-3 flat)"},
            {"name": "laplace_sigma", "args": {"g": -5.0, "band_rel": 0.15,
                                               "sigma_ref": 3.32e-2},
             "source": "VALIDATION T11 Laplace (sigma = 3.32e-2 +-10% slope, "
                       "+-5% per droplet; +-15% single-droplet combined)"},
            {"name": "mass_drift", "args": {"region": "all", "band_per_1e4": 1e-10},
             "source": "VALIDATION T11 (total-mass drift <= 1e-10 rel)"},
            {"name": "xy_mirror_symmetry", "args": {"field": "rho", "band": 1e-9},
             "source": "symmetric IC + symmetric stencil (observation band)"},
        ],
    },
    {
        "id": "droplet-on-wall-t11c",
        "track": "multiphase",
        "scenario": _sc(
            "qa-droplet-wall", 160, 100, 1.0 / 6.0,
            {"left": PER, "right": PER, "bottom": BB, "top": BB},
            30_000,
            init={"kind": "droplet", "cx": 80.0, "cy": 1.0, "r": 22.0,
                  "rhoLiquid": 2.0, "rhoVapor": 0.15},
            multiphase={"g": -5.0, "gWall": 0.0, "wallRho": 1.0},
            outputs=[{"field": "rho", "format": "csv", "every": 0},
                     {"field": "rho", "format": "png", "every": 0}],
        ),
        "checks": [
            {"name": "contact_angle", "args": {"theta_ref": 63.0, "tol": 8.0},
             "source": "VALIDATION T11c (wallRho=1.0 -> 63 deg +- 8)"},
        ],
    },
    # ------------------------------------------------------------ robustness
    {
        "id": "tau-margin-cavity-t10",
        "track": "robustness",
        "scenario": _sc(
            "qa-tau-margin", 128, 128, 1.0 / 300.0,
            {"left": BB, "right": BB, "bottom": BB, "top": mw(0.05)},
            10_000,
            outputs=[{"field": "speed", "format": "png", "every": 0},
                     {"field": "rho", "format": "csv", "every": 2000}],
        ),
        "expect_warnings": ["physics.nu"],  # tau = 0.51
        "checks": [
            {"name": "speed_scale", "args": {"u_ref": 0.05, "factor": 2.0},
             "source": "VALIDATION T10 (tau=0.51 N=128 U=0.05 stable, max|u|=0.046)"},
        ],
    },
]

# -- duct3d needs nz > 1 through _sc: patch the entry properly ---------------
_DUCT = _sc(
    "qa-duct3d", 12, 34, 0.1,
    {"left": PER, "right": PER, "bottom": BB, "top": BB, "front": BB, "back": BB},
    100_000, nz=34, force=(1e-6, 0.0), steady=(1e-11, 500),
    outputs=[{"field": "ux", "format": "vtk", "every": 0},
             {"field": "rho", "format": "vtk", "every": 4000},
             {"field": "rho", "format": "vtk", "every": 0}],
)
for _c in CONFIGS:
    if _c["id"] == "duct3d-t15":
        _c["scenario"] = _DUCT
        _c["checks"] = [
            {"name": "duct_exact", "args": {"g": 1e-6, "nu": 0.1, "band": 1e-3},
             "source": "VALIDATION T15.2 (TRT LInf_rel <= 1e-3; measured 2.3e-4 at 32^2)"},
            {"name": "duct_flow_rate", "args": {"g": 1e-6, "nu": 0.1, "band": 5e-3},
             "source": "VALIDATION T15.2 (Q within +-0.5%; measured 0.094%)"},
            {"name": "duct_yz_symmetry", "args": {"band": 1e-11},
             "source": "square-duct y<->z exchange symmetry (discrete-exact)"},
            {"name": "mass_drift", "args": {"region": "interior3d", "band_per_1e4": 1e-11},
             "source": "VALIDATION T6"},
        ]

# Every config implicitly gets: finite-fields / manifest-status check,
# |u| <= 0.3 hard ceiling, and warning-expectation audit (run_sweep.py).
