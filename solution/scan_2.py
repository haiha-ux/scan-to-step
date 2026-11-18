"""scan_2: block with a rectangular pocket and a counterbored through hole on the +x face; +x edges chamfered."""
import json, sys
import numpy as np
sys.path.insert(0, 'solution')
import meas as M
from build123d import Box, Cylinder, Pos, Rot, Align, export_step, Axis, chamfer, fillet

m, mm, to, good = M.load(2)
V = mm.vertices[np.unique(mm.faces[good].ravel())]
fp = lambda ax, sg, g: M.face_plane(mm, good, ax, sg, g)[0]
X_p, X_m, Y_p, Y_m, Z_p, Z_m = fp(0, 1, 15), fp(0, -1, -15), fp(1, 1, 20), fp(1, -1, -20), fp(2, 1, 30), fp(2, -1, -30)
X_pf, X_cb = fp(0, 1, 7.0), fp(0, 1, 8.5)                       # pocket floor, counterbore floor
PY0, PY1, PZ0, PZ1 = fp(1, 1, -10.1), fp(1, -1, 9.9), fp(2, 1, -22.0), fp(2, -1, 2.0)

# counterbore and through hole: circles in (y, z)
s = (V[:, 0] > X_cb + 0.4) & (V[:, 0] < X_p - 0.4) & (np.abs(np.hypot(V[:, 1] - 0.2, V[:, 2] - 16.8) - 5.4) < 0.6)
(cb_y, cb_z, cb_r), cb_sd, cb_n = M.fit_circle2d(V[s][:, [1, 2]], [0.2, 16.8, 5.4])
s = (V[:, 0] > X_m + 0.4) & (V[:, 0] < X_cb - 0.4) & (np.abs(np.hypot(V[:, 1] - 0.2, V[:, 2] - 16.8) - 3.2) < 0.6)
(th_y, th_z, th_r), th_sd, th_n = M.fit_circle2d(V[s][:, [1, 2]], [0.2, 16.8, 3.2])

# chamfers on the four +x edges: leg along x and along the side face
legs = []
for ax, sg, edge in ((1, 1, Y_p), (1, -1, Y_m), (2, 1, Z_p), (2, -1, Z_m)):
    other = 2 if ax == 1 else 1
    s = (V[:, 0] > X_p - 1.6) & (sg * V[:, ax] > sg * edge - 1.6) & (np.abs(V[:, other]) < (12 if ax == 1 else 16))
    u, w = X_p - V[s, 0], sg * (edge - V[s, ax])                       # depth below the face, inset from the side
    band = (u > 0.2) & (w > 0.2) & (u < 1.4) & (w < 1.4)
    # chamfer plane: u + w = c (45 deg) -> c = leg; least squares on the band
    c = float(np.median((u + w)[band]))
    legs.append(c)
C = float(np.mean(legs))

# pocket corner radius (edges along x): tangent rounds in (y, z) at the 4 corners
pr = []
for (yc, dy), (zc, dz) in [((PY0, 1), (PZ0, 1)), ((PY1, -1), (PZ0, 1)), ((PY0, 1), (PZ1, -1)), ((PY1, -1), (PZ1, -1))]:
    s = (V[:, 0] > X_pf + 0.5) & (V[:, 0] < X_p - 1.5) & (np.abs(V[:, 1] - yc) < 4) & (np.abs(V[:, 2] - zc) < 4) \
        & ((V[:, 1] - yc) * dy > -0.1) & ((V[:, 2] - zc) * dz > -0.1)
    pr.append(M.fit_tangent_round(V[s][:, [1, 2]], (yc, zc), dy, dz, 2.0))
meas = dict(planes=dict(X_p=X_p, X_m=X_m, Y_p=Y_p, Y_m=Y_m, Z_p=Z_p, Z_m=Z_m, X_pocket_floor=X_pf, X_cbore_floor=X_cb,
                        pocket=[PY0, PY1, PZ0, PZ1]), cbore=[cb_y, cb_z, cb_r, cb_sd, cb_n], thru=[th_y, th_z, th_r, th_sd, th_n],
            chamfer_legs=legs, pocket_corner_rounds=pr)
print(json.dumps(meas, indent=1))
R_pk = float(np.mean([p[0] for p in pr]))
use_pk_round = all(p[3] > 100 for p in pr) and np.mean([p[2] for p in pr]) < 0.05   # too few points = sharp corner

# model
A0 = (Align.MIN, Align.MIN, Align.MIN)
part = Pos(X_m, Y_m, Z_m) * Box(X_p - X_m, Y_p - Y_m, Z_p - Z_m, align=A0)
edges = part.edges().filter_by(lambda e: abs(e.center().X - X_p) < 1e-6)
part = chamfer(edges, C)
pocket = Pos(X_pf, PY0, PZ0) * Box(X_p - X_pf + 1, PY1 - PY0, PZ1 - PZ0, align=A0)
if use_pk_round:
    pocket = fillet(pocket.edges().filter_by(Axis.X), R_pk)
part -= pocket
cyl = lambda y, z, r, x0, x1: Pos(x0, y, z) * Rot(0, 90, 0) * Cylinder(r, x1 - x0, align=(Align.CENTER, Align.CENTER, Align.MIN))
cy, cz = (cb_y * cb_n + th_y * th_n) / (cb_n + th_n), (cb_z * cb_n + th_z * th_n) / (cb_n + th_n)   # one axis for both bores
part -= cyl(cy, cz, cb_r, X_cb, X_p + 1) + cyl(cy, cz, th_r, X_m - 1, X_cb + 0.01)
part = M.to_build123d_location(np.linalg.inv(to)) * part
assert len(part.solids()) == 1 and part.is_valid
export_step(part, 'solution/scan_2.step')

pocket_c = ((PY0 + PY1) / 2, (PZ0 + PZ1) / 2)
dist = abs(cz - pocket_c[1])   # projection on the long axis (z, 60 mm), as defined by the brief owner
dims = [{"type": "hole_diameter", "value": round(2 * cb_r, 3), "name": "counterbore"},
        {"type": "hole_diameter", "value": round(2 * th_r, 3), "name": "through hole"},
        {"type": "depth", "value": round(X_p - X_pf, 3), "name": "pocket depth"},
        {"type": "depth", "value": round(X_p - X_cb, 3), "name": "counterbore depth"},
        {"type": "distance", "value": round(dist, 3), "name": "pocket centre to hole axis, projected on the long (60 mm) axis"},
        {"type": "bbox", "value": sorted([round(X_p - X_m, 3), round(Y_p - Y_m, 3), round(Z_p - Z_m, 3)])}]
json.dump({"dims": dims}, open('solution/scan_2_peer.json', 'w'), indent=1)
print(json.dumps(dims), 'chamfer', round(C, 3), 'pocket round', round(R_pk, 3), use_pk_round)
d_, s_ = M.scan_to_model(m, 'solution/scan_2.step')

print('scan->model p95 %.4f max %.4f mean %.4f volume %.1f' % (np.percentile(d_, 95), d_.max(), d_.mean(), s_.volume))
