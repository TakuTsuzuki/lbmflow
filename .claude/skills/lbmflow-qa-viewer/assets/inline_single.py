#!/usr/bin/env python3
"""Inline three.js + OrbitControls + the LBM volume into a self-contained
viewer.html.  Usage: build_viewer.py <data_dir> <template> <three.min.js>
<OrbitControls.js> <out.html>"""
import sys, base64, pathlib

data_dir, tmpl, three, orbit, out = map(pathlib.Path, sys.argv[1:6])

meta = (data_dir / "volume.json").read_text().strip()
vol_b64 = base64.b64encode((data_dir / "volume.bin").read_bytes()).decode()

html = tmpl.read_text()
html = html.replace("/*__THREE__*/", three.read_text())
html = html.replace("/*__ORBIT__*/", orbit.read_text())
html = html.replace("/*__META__*/", meta)
html = html.replace("/*__DATA__*/", vol_b64)

out.write_text(html)
print(f"wrote {out}  ({len(html)/1e6:.2f} MB, volume {len(vol_b64)/1e6:.2f} MB b64)")
