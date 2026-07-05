#!/usr/bin/env python3
"""Assemble the multi-dataset BUNDLE (geometry + N stir-speed volumes) and inline
three.js/OrbitControls into dash2_template.html. Emits a full local HTML.
Usage: build_dash2.py <template> <three> <orbit> <out.html> <label:dir> [label:dir ...]"""
import sys, json, base64, pathlib
tmpl, three, orbit, out = map(pathlib.Path, sys.argv[1:5])
pairs = sys.argv[5:]

datasets, geom = [], None
for p in pairs:
    label, d = p.split(":", 1)
    d = pathlib.Path(d)
    m = json.loads((d / "volume.json").read_text())
    if geom is None:
        geom = {k: m[k] for k in ("vn","n","cx","cy","r_tank","zc","tip_r","disk_r","hub_r",
                "shaft_r","blade_hh","disk_hh","n_blades","blade_hw","baffle_len","baffle_hw")}
    b64 = base64.b64encode((d / "volume.bin").read_bytes()).decode()
    datasets.append({"label": label, "u_tip": m["u_tip"], "omega": m["omega"],
                     "speed_max": m["speed_max"], "shear_max": m["shear_max"], "vol": b64})

bundle = dict(geom); bundle["datasets"] = datasets
bundle_js = json.dumps(bundle)

html = tmpl.read_text()
html = html.replace("/*__THREE__*/", three.read_text())
html = html.replace("/*__ORBIT__*/", orbit.read_text())
html = html.replace("/*__BUNDLE__*/", bundle_js)
out.write_text(html)
mb = len(html)/1e6
print(f"wrote {out}  ({mb:.2f} MB, {len(datasets)} datasets, vn={geom['vn']})")
