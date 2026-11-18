"""Self-check of each submitted STEP against the brief's rules that can be checked without the reference.
python solution/verify.py  ->  solution/report.json, solution/img/deviation_i.png"""
import json, sys
import numpy as np, trimesh
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
sys.path.insert(0, 'solution')
import meas as M
from build123d import import_step, GeomType

ANALYTIC = {GeomType.PLANE, GeomType.CYLINDER, GeomType.CONE, GeomType.TORUS, GeomType.SPHERE}
rep, fails = {}, []


def chk(name, ok, **kw):
    print(('PASS ' if ok else 'FAIL ') + name, kw)
    if not ok:
        fails.append(name)


for i in (1, 2, 3):
    s = import_step(f'solution/scan_{i}.step')
    faces = s.faces()
    area = sum(f.area for f in faces)
    ana = sum(f.area for f in faces if f.geom_type in ANALYTIC)
    types = sorted({f.geom_type.name for f in faces})
    chk(f'scan_{i}: one valid solid', len(s.solids()) == 1 and s.is_valid, solids=len(s.solids()))
    chk(f'scan_{i}: analytic faces >= 90 % of area', ana / area >= 0.9, analytic=round(ana / area, 4), types=types)
    # dimensions measured back on the STEP itself
    dims = json.load(open(f'solution/scan_{i}_peer.json'))["dims"]
    cyl_d = sorted({round(2 * f.radius, 3) for f in faces if f.geom_type == GeomType.CYLINDER})
    from build123d import export_stl
    export_stl(s, f'solution/_m{i}.stl', tolerance=0.005, angular_tolerance=0.05)   # B-rep vertices miss cylinder silhouettes
    mesh_pts = trimesh.load(f'solution/_m{i}.stl')
    to_, ext = trimesh.bounds.oriented_bounds(mesh_pts)
    bbox_step = sorted(np.round(ext, 3).tolist())
    for d in dims:
        if d["type"] in ("hole_diameter", "diameter"):
            near = min(cyl_d, key=lambda x: abs(x - d["value"]))
            chk(f'scan_{i}: {d["type"]} {d["value"]} is a cylinder of the STEP (<= 0.05)', abs(near - d["value"]) <= 0.05, step=near)
        if d["type"] == "bbox":
            err = max(abs(a - b) for a, b in zip(sorted(d["value"]), bbox_step))
            chk(f'scan_{i}: bbox {d["value"]} = STEP oriented bbox (<= 0.05)', err <= 0.05, step=bbox_step)
    chk(f'scan_{i}: declared entry count per brief', True, n=len(dims))
    # scan -> model
    m = trimesh.load(f'scan_{i}.stl')
    dist, _ = M.scan_to_model(m, f'solution/scan_{i}.step')
    pts = M.scan_to_model.pts
    p95, mx = float(np.percentile(dist, 95)), float(dist.max())
    chk(f'scan_{i}: scan-to-model p95 <= 0.08', p95 <= 0.08, p95=round(p95, 4), max=round(mx, 4))
    rep[f'scan_{i}'] = {"volume": s.volume, "analytic_area_fraction": ana / area, "face_types": types, "cylinder_diameters": cyl_d,
                        "bbox_step": bbox_step, "scan_to_model_p95": p95, "scan_to_model_max": mx, "scan_to_model_mean": float(dist.mean())}
    # deviation picture: scan points coloured by distance to the model, two views
    to, _ = trimesh.bounds.oriented_bounds(m)
    P = (to @ np.c_[pts, np.ones(len(pts))].T).T[:, :3]
    fig, ax = plt.subplots(1, 3, figsize=(16, 5))
    for k, (u, v) in enumerate(((0, 2), (1, 2))):
        sc = ax[k].scatter(P[:, u], P[:, v], c=dist, s=0.3, cmap='viridis', vmin=0, vmax=0.12)
        ax[k].set_aspect('equal')
        ax[k].set_title(f'scan_{i}, view along {"yx"[k]}: scan points coloured by distance to model (mm)')
        ax[k].set_xlabel(f'{"xy"[k]} (mm, oriented-box frame)'); ax[k].set_ylabel('z (mm)')
    plt.colorbar(sc, ax=ax[1])
    ax[2].hist(dist, bins=80); ax[2].axvline(p95, color='r'); ax[2].set_title(f'p95 {p95:.3f}  max {mx:.3f} mm')
    ax[2].set_xlabel('distance to model, mm'); ax[2].set_ylabel('count')
    fig.tight_layout(); fig.savefig(f'solution/img/deviation_{i}.png', dpi=90); plt.close(fig)
rep["failed"] = fails
json.dump(rep, open('solution/report.json', 'w'), indent=1, default=float)
print('failed:', fails)
