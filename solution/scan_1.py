"""scan_1: L-bracket, inner bend fillet, two rounded top corners, two through holes in the upright."""
import json, sys
import numpy as np
sys.path.insert(0, 'solution')
import meas as M
from build123d import (Box, Cylinder, Pos, Rot, Align, export_step, Axis, fillet, Edge, GeomType)

m, mm, to, good = M.load(1, refine=False)   # frame choice: OBB fits this scan better than plane re-alignment (p95 0.059 vs 0.067)
V = mm.vertices[np.unique(mm.faces[good].ravel())]
fp = lambda ax, sg, g: M.face_plane(mm, good, ax, sg, g)[0]
X_out, X_in, X_end = fp(0, 1, 14.99), fp(0, -1, 10.99), fp(0, -1, -15.0)
Y_p, Y_m = fp(1, 1, 20.0), fp(1, -1, -20.0)
Z_bot, Z_base, Z_top = fp(2, -1, -30.0), fp(2, 1, -26.0), fp(2, 1, 30.0)

# inner fillet and top corner rounds: radius fitted with the arc tangent to both neighbouring planes
s = (V[:, 0] > X_in - 4.5) & (V[:, 0] < X_in - 0.05) & (V[:, 2] > Z_base + 0.05) & (V[:, 2] < Z_base + 4.5) & (np.abs(V[:, 1]) < 19)
R_in, mean_in, sd_in, n_in = M.fit_tangent_round(V[s][:, [0, 2]], (X_in, Z_base), -1, 1, 3.0)
corners = {}
for sg in (1, -1):
    yedge = Y_p if sg > 0 else Y_m
    s = (V[:, 0] > X_in + 0.3) & (V[:, 0] < X_out - 0.3) & (sg * V[:, 1] > sg * yedge - 7) & (V[:, 2] > Z_top - 7)
    corners[sg] = M.fit_tangent_round(V[s][:, [1, 2]], (yedge, Z_top), -sg, -1, 5.0)
# holes: circles in (y, z) from points inside the upright slab, away from the outer boundary
holes = []
for g in ((-10.0, 15.0, 3.2), (10.0, 15.0, 4.2)):
    s = (V[:, 0] > X_in + 0.4) & (V[:, 0] < X_out - 0.4) & (np.abs(np.hypot(V[:, 1] - g[0], V[:, 2] - g[1]) - g[2]) < 0.6)
    holes.append(M.fit_circle2d(V[s][:, [1, 2]], list(g)))
meas = {"planes": dict(X_out=X_out, X_in=X_in, X_end=X_end, Y_p=Y_p, Y_m=Y_m, Z_bot=Z_bot, Z_base=Z_base, Z_top=Z_top),
        "inner_fillet": [R_in, mean_in, sd_in, n_in], "top_corners": {str(k): list(v) for k, v in corners.items()},
        "holes": [[*h[0], h[1], h[2]] for h in holes]}
print(json.dumps(meas, indent=1))
R_top = float(np.mean([corners[1][0], corners[-1][0]]))

# model in the OBB frame, then back to the scan frame
base = Pos(X_end, Y_m, Z_bot) * Box(X_out - X_end, Y_p - Y_m, Z_base - Z_bot, align=(Align.MIN, Align.MIN, Align.MIN))
up = Pos(X_in, Y_m, Z_bot) * Box(X_out - X_in, Y_p - Y_m, Z_top - Z_bot, align=(Align.MIN, Align.MIN, Align.MIN))
part = base + up
inner = part.edges().filter_by(Axis.Y).filter_by(lambda e: abs(e.center().X - X_in) < 1e-3 and abs(e.center().Z - Z_base) < 1e-3)
part = fillet(inner, R_in)
tops = part.edges().filter_by(Axis.X).filter_by(lambda e: abs(e.center().Z - Z_top) < 1e-3)
part = fillet(tops, R_top)
for (hy, hz, hr), _, _ in holes:
    part -= Pos(X_in - 1, hy, hz) * Rot(0, 90, 0) * Cylinder(hr, X_out - X_in + 2, align=(Align.CENTER, Align.CENTER, Align.MIN))
part = M.to_build123d_location(np.linalg.inv(to)) * part
assert len(part.solids()) == 1 and part.is_valid
export_step(part, 'solution/scan_1.step')

(h1y, h1z, h1r), (h2y, h2z, h2r) = holes[0][0], holes[1][0]
dims = [{"type": "hole_diameter", "value": round(2 * h1r, 3)}, {"type": "hole_diameter", "value": round(2 * h2r, 3)},
        {"type": "distance", "value": round(float(np.hypot(h2y - h1y, h2z - h1z)), 3)},
        {"type": "thickness", "value": round(X_out - X_in, 3)},
        {"type": "bbox", "value": sorted([round(X_out - X_end, 3), round(Y_p - Y_m, 3), round(Z_top - Z_bot, 3)])}]
json.dump({"dims": dims, "note": "plate thickness = upright (x); base thickness measured %.3f" % (Z_base - Z_bot)},
          open('solution/scan_1_peer.json', 'w'), indent=1)
print(json.dumps(dims), 'base t', round(Z_base - Z_bot, 3), 'R_in', round(R_in, 3), 'R_top', round(R_top, 3))
dist, s_ = M.scan_to_model(m, 'solution/scan_1.step')
print('scan->model p95 %.4f max %.4f mean %.4f volume %.1f' % (np.percentile(dist, 95), dist.max(), dist.mean(), s_.volume))
