#!/usr/bin/env python3
"""Stream-parse tray_velocity.vtk (ASCII STRUCTURED_POINTS, VECTORS), block-average
downsample by DS, write vel.bin (f32 LE, (k*ny+j)*nx+i interleaved vx,vy,vz) +
meta.json (dims, spacing, speed_max, density 12x12, metrics)."""
import sys, json, struct, csv, pathlib
import numpy as np

src = pathlib.Path(sys.argv[1])          # dataset dir
dst = pathlib.Path(sys.argv[2])          # out dir
DS = int(sys.argv[3]) if len(sys.argv) > 3 else 2
dst.mkdir(parents=True, exist_ok=True)

# --- parse VTK header + stream vectors ---
f = open(src / "tray_velocity.vtk")
dims = spacing = None
for line in f:
    t = line.split()
    if not t: continue
    if t[0] == "DIMENSIONS": dims = tuple(map(int, t[1:4]))
    elif t[0] == "SPACING": spacing = float(t[1])
    elif t[0] == "VECTORS": break
nx, ny, nz = dims
ox, oy, oz = nx // DS, ny // DS, nz // DS
acc = np.zeros((oz, oy, ox, 3), dtype=np.float64)
cnt = np.zeros((oz, oy, ox, 1), dtype=np.float64)
plane = nx * ny
idx = 0
buf = np.loadtxt(f, dtype=np.float64)    # (N,3) — loads the numeric block
assert buf.shape == (nx * ny * nz, 3), buf.shape
buf = buf.reshape(nz, ny, nx, 3)
# block average
buf = buf[:oz*DS, :oy*DS, :ox*DS, :]
v = buf.reshape(oz, DS, oy, DS, ox, DS, 3).mean(axis=(1, 3, 5)).astype(np.float32)
speed = np.sqrt((v.astype(np.float64) ** 2).sum(axis=-1))
speed_max = float(speed.max())

(dst / "vel.bin").write_bytes(v.astype("<f4").tobytes())

# --- density 12x12 ---
dens = {}
with open(src / "density.csv") as fh:
    for row in csv.DictReader(fh):
        dens[(int(row["bin_i"]), int(row["bin_j"]))] = (float(row["normalized_density"]), int(row["count"]))
pi = max(k[0] for k in dens) + 1
pj = max(k[1] for k in dens) + 1
dgrid = [[dens[(i, j)][0] for i in range(pi)] for j in range(pj)]
cgrid = [[dens[(i, j)][1] for i in range(pi)] for j in range(pj)]

metrics = json.loads((src / "metrics.json").read_text())
meta = dict(nx=ox, ny=oy, nz=oz, ds=DS, spacing=spacing * DS,
            speed_max=speed_max, density=dgrid, counts=cgrid,
            pi=pi, pj=pj, metrics=metrics)
(dst / "meta.json").write_text(json.dumps(meta))
print(f"{src.name}: {ox}x{oy}x{oz} speed_max={speed_max:.4e} m/s -> {dst}")
