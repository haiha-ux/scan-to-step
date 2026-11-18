"""scan_3: stepped turned shaft with a groove. Measure, rebuild as a parametric solid, export STEP + JSON."""
import json, sys
import numpy as np
sys.path.insert(0, 'solution')
import meas as M
from build123d import Cylinder, Pos, Align, export_step, Location, Matrix

m, mm, to, good = M.load(3)
v = mm.vertices
# 1. axis from the 20 mm section, well away from its ends
band = v[(v[:, 2] > -15) & (v[:, 2] < 0)]
band = band[np.abs(np.hypot(band[:, 0], band[:, 1]) - 10) < 0.5]
c, d, _, _ = M.fit_cyl(band, np.array([0, 0, 1.0]), np.zeros(3), 10)
w = v - c; zz = w @ d; rr = np.linalg.norm(w - np.outer(zz, d), axis=1)
n = mm.vertex_normals; na = n @ d                     # axial component of the normal


def radius(z0, z1, guess):
    s = (zz > z0 + 0.4) & (zz < z1 - 0.4) & (np.abs(na) < 0.3) & (np.abs(rr - guess) < 0.5)
    return float(np.mean(rr[s])), int(s.sum()), float(np.std(rr[s]))


fc = mm.triangles_center[good]; fn = mm.face_normals[good]; fa = mm.area_faces[good]
fz = (fc - c) @ d; fna = fn @ d


def plane(z_guess, sign):
    """End/step faces: area-weighted mean axial position of the triangles facing that way."""
    s = (fna * sign > 0.98) & (np.abs(fz - z_guess) < 0.4)
    return float(np.average(fz[s], weights=fa[s])), int(s.sum())


zb, _ = plane(-37.5, -1); zs1, _ = plane(-17.5, -1); zg0, _ = plane(2.5, 1); zg1, _ = plane(4.5, -1)
zs2, _ = plane(12.5, -1); zt, _ = plane(37.5, 1)
r1 = radius(zb, zs1, 7.0); r2a = radius(zs1, zg0, 10.0); r2b = radius(zg1, zs2, 10.0)
rg = radius(zg0, zg1, 8.5); r3 = radius(zs2, zt, 14.0)
r2 = (r2a[0] * r2a[1] + r2b[0] * r2b[1]) / (r2a[1] + r2b[1])
meas = {"axis_point": c.tolist(), "axis_dir": d.tolist(), "z": [zb, zs1, zg0, zg1, zs2, zt],
        "r": {"d14": r1, "d20_lower": r2a, "d20_upper": r2b, "groove": rg, "d28": r3}}
print(json.dumps(meas, indent=1))

# 2. model, in the axis frame: z from the bottom face
L = Location
seg = lambda z0, z1, r: Pos(0, 0, z0) * Cylinder(r, z1 - z0, align=(Align.CENTER, Align.CENTER, Align.MIN))
part = seg(zb, zs1, r1[0]) + seg(zs1, zg0, r2) + seg(zg0, zg1, rg[0]) + seg(zg1, zs2, r2) + seg(zs2, zt, r3[0])
T = M.frame_from_axis(c, d)
part = M.to_build123d_location(np.linalg.inv(to) @ T) * part        # back into the scan's own frame
assert len(part.solids()) == 1 and part.is_valid
export_step(part, 'solution/scan_3.step')

# 3. declared dimensions, measured on the scan; the same values drive the model
dims = [{"type": "diameter", "value": round(2 * r1[0], 3)}, {"type": "diameter", "value": round(2 * r2, 3)},
        {"type": "diameter", "value": round(2 * rg[0], 3)}, {"type": "diameter", "value": round(2 * r3[0], 3)},
        {"type": "width", "value": round(zg1 - zg0, 3)},
        {"type": "bbox", "value": sorted([round(2 * r3[0], 3), round(2 * r3[0], 3), round(zt - zb, 3)])}]
json.dump({"dims": dims}, open('solution/scan_3_peer.json', 'w'), indent=1)
print(json.dumps(dims))
dist, s = M.scan_to_model(m, 'solution/scan_3.step')
print('scan->model p95 %.4f max %.4f mean %.4f  volume %.2f' % (np.percentile(dist, 95), dist.max(), dist.mean(), s.volume))
