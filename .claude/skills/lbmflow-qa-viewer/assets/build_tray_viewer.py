#!/usr/bin/env python3
"""Inline three.js/OrbitControls/datasets into tray_template.html.
Usage: build_tray_viewer.py <template> <three> <orbit> <out.html> <label:dir> ..."""
import sys, json, base64, pathlib
tmpl, three, orbit, out = map(pathlib.Path, sys.argv[1:5])
datasets = []
for p in sys.argv[5:]:
    label, d = p.split(":", 1)
    d = pathlib.Path(d)
    datasets.append({
        "label": label,
        "meta": json.loads((d / "meta.json").read_text()),
        "vol": base64.b64encode((d / "vel.bin").read_bytes()).decode(),
    })
html = tmpl.read_text()
html = html.replace("/*__THREE__*/", three.read_text())
html = html.replace("/*__ORBIT__*/", orbit.read_text())
html = html.replace("/*__BUNDLE__*/", json.dumps({"datasets": datasets}))
out.write_text(html)
print(f"wrote {out} ({len(html)/1e6:.2f} MB, {len(datasets)} datasets)")
