"""Measurement helpers: OBB frame, plane peaks per axis, least-squares plane and cylinder fits."""
import numpy as np, trimesh
from scipy.optimize import least_squares


def load(i, refine=True):
    m = trimesh.load(f'scan_{i}.stl')
    to, ext = trimesh.bounds.oriented_bounds(m)
    mm = m.copy(); mm.apply_transform(to)
    # drop fill triangles: anything much larger than the median triangle is a patched scan hole
    a = mm.area_faces
    good = a < 20 * np.median(a)
    if refine:
        # the minimum-volume OBB need not line up with the faces (chamfers, pockets): re-align it on the big planes
        for _ in range(3):
            n = mm.face_normals; aa = mm.area_faces
            axes = []
            for k in range(3):
                s_ = good & (np.abs(n[:, k]) > 0.995)          # flat faces only: partly scanned hole walls bias a looser cut
                nk = (n[s_] * np.sign(n[s_, k])[:, None] * aa[s_, None]).sum(0)
                axes.append(nk / np.linalg.norm(nk))
            A = np.array(axes)
            U, _, Wt = np.linalg.svd(A)
            Rm = U @ Wt                                    # nearest rotation to the measured axes
            T = np.eye(4); T[:3, :3] = Rm
            mm.apply_transform(T); to = T @ to
    return m, mm, to, good


def fit_cyl(pts, axis0, c0, r0):
    """Cylinder: axis direction (2 angles), point on axis (2 coords), radius."""
    axis0 = axis0 / np.linalg.norm(axis0)
    # local frame
    t = np.array([1, 0, 0]) if abs(axis0[0]) < 0.9 else np.array([0, 1, 0])
    e1 = np.cross(axis0, t); e1 /= np.linalg.norm(e1); e2 = np.cross(axis0, e1)

    def res(p):
        a, b, x, y, r = p
        d = axis0 + a * e1 + b * e2; d /= np.linalg.norm(d)
        q = c0 + x * e1 + y * e2
        v = pts - q
        return np.linalg.norm(v - np.outer(v @ d, d), axis=1) - r
    s = least_squares(res, [0, 0, 0, 0, r0], loss='soft_l1', f_scale=0.05)
    a, b, x, y, r = s.x
    d = axis0 + a * e1 + b * e2; d /= np.linalg.norm(d)
    return c0 + x * e1 + y * e2, d, abs(r), float(np.std(res(s.x)))


def frame_from_axis(c, d):
    """4x4 transform: local (z along d, origin c) -> OBB frame."""
    d = d / np.linalg.norm(d)
    t = np.array([1.0, 0, 0]) if abs(d[0]) < 0.9 else np.array([0, 1.0, 0])
    e1 = np.cross(t, d); e1 /= np.linalg.norm(e1); e2 = np.cross(d, e1)
    T = np.eye(4); T[:3, 0], T[:3, 1], T[:3, 2], T[:3, 3] = e1, e2, d, c
    return T


def to_build123d_location(T):
    from build123d import Location, Plane, Vector
    pl = Plane(origin=Vector(*T[:3, 3]), x_dir=Vector(*T[:3, 0]), z_dir=Vector(*T[:3, 2]))
    return Location(pl)


def scan_to_model(scan_mesh, step_path, samples=60000):
    """Distances from scan vertices to the model surface (model tessellated finely)."""
    from build123d import import_step, export_stl
    import tempfile, os
    s = import_step(step_path)
    f = os.path.join(tempfile.mkdtemp(), 'm.stl')
    export_stl(s, f, tolerance=0.003, angular_tolerance=0.05)
    mod = trimesh.load(f)
    pts = scan_mesh.vertices
    if len(pts) > samples:
        pts = pts[np.random.default_rng(0).choice(len(pts), samples, replace=False)]
    d = np.concatenate([trimesh.proximity.closest_point(mod, pts[k:k + 5000])[1] for k in range(0, len(pts), 5000)])
    scan_to_model.pts = pts
    return d, s


def fit_circle2d(Q, guess):
    """Geometric least-squares circle in 2D, robust loss. Q: (n, 2)."""
    def res(p):
        return np.hypot(Q[:, 0] - p[0], Q[:, 1] - p[1]) - p[2]
    s = least_squares(res, guess, loss='soft_l1', f_scale=0.05)
    r = res(s.x)
    return s.x, float(np.std(r)), int(len(Q))


def face_plane(mm, good, axis, sign, guess, band=0.3):
    """Area-weighted mean position of triangles facing (axis, sign) within band of guess."""
    c = mm.triangles_center; n = mm.face_normals; a = mm.area_faces
    s = good & (n[:, axis] * sign > 0.98) & (np.abs(c[:, axis] - guess) < band)
    return float(np.average(c[s, axis], weights=a[s])), int(s.sum())


def fit_tangent_round(Q, corner, du, dv, r0=3.0):
    """Round tangent to two perpendicular planes meeting at `corner` (2D), material on the (du, dv) side.
    Centre = corner + R*(du, dv); only R is free. Points Q (n, 2) near the arc."""
    corner = np.asarray(corner, float)

    def stats(R):
        cc = corner + R * np.array([du, dv])
        d = np.hypot(Q[:, 0] - cc[0], Q[:, 1] - cc[1]) - R
        inside = ((Q[:, 0] - corner[0]) * du < R) & ((Q[:, 1] - corner[1]) * dv < R)   # points in the arc sector
        d = d[inside & (np.abs(d) < 0.3)]
        return np.sqrt(np.mean(d ** 2)) if len(d) > 20 else 9.0
    grid = np.arange(0.6 * r0, 1.4 * r0, 0.002)        # grid search on RMS: a bounded optimiser found a false minimum
    R = float(grid[np.argmin([stats(g) for g in grid])])

    class _S: pass
    s = _S(); s.x = R
    R = float(s.x); cc = corner + R * np.array([du, dv])
    d = np.hypot(Q[:, 0] - cc[0], Q[:, 1] - cc[1]) - R
    sel = ((Q[:, 0] - corner[0]) * du < R) & ((Q[:, 1] - corner[1]) * dv < R) & (np.abs(d) < 0.3)
    return R, float(np.mean(d[sel])), float(np.std(d[sel])), int(sel.sum())
