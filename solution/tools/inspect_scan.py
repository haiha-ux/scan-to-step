import sys, numpy as np, trimesh
for i in (1,2,3):
    m=trimesh.load(f'scan_{i}.stl')
    to, ext = trimesh.bounds.oriented_bounds(m)
    mm = m.copy(); mm.apply_transform(to)
    print(f'scan_{i}: faces {len(m.faces)} watertight {m.is_watertight} area {m.area:.1f} obb {np.round(np.sort(ext),3)}')
    # cluster face normals (area weighted) in the OBB frame
    n = mm.face_normals; a = mm.area_faces
    q = np.round(n, 1)
    keys, inv = np.unique(q, axis=0, return_inverse=True)
    s = np.bincount(inv.ravel(), weights=a)
    order = np.argsort(-s)[:12]
    for k in order: print('   n', keys[k], 'area', round(s[k],1))
