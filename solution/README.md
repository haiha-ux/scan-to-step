# Solution: scan to parametric STEP

Each part is rebuilt in build123d from planes and cylinders only. No mesh or faceted face is used. Each model is
placed back in the scan's own position and rotation, so the scan-to-model deviation can be measured directly.

| Part | Rebuilt as | Scan → model (mm) |
|---|---|---|
| scan_1 | L-bracket 29.99 × 40.00 × 60.00. Base and upright 4.00 thick. Inner bend round R 2.99, top corners R 4.99. Through holes Ø5.998 and Ø8.000 in the upright, 20.001 apart | p95 0.059, max 0.133 |
| scan_2 | Block 30.00 × 40.00 × 60.00. 45° chamfer 1.50 on the four edges of the machined face. Pocket 20.00 × 24.00 × 8.00 deep. Counterbore Ø10.999 × 6.50 deep, coaxial with a Ø6.600 through hole | p95 0.059, max 0.126 |
| scan_3 | Turned shaft 75.00 long: Ø13.999 × 20.0, Ø19.996 (split by a 2.003 wide groove, groove bottom Ø16.990), Ø27.991 × 25.0 | p95 0.059, max 0.134 |

The noise in the brief is 0.03 mm. A p95 of 0.059 mm is about 2σ, which is what a correct surface gives on its own.

Every declared value is the measured value, not rounded to a nominal size. The model is built from the same numbers,
so each value in `scan_i_peer.json` is also a dimension of `scan_i.step`. `verify.py` reads them back from the STEP.

## Method

1. Frame. The scan is moved into its oriented bounding box (trimesh). For scan_2, the box frame was then
   re-aligned on the measured normals of the large flat faces: the minimum-volume box was 0.8° off because of the
   chamfers and the pocket, which put scan points 0.4 mm outside the model. For scan_1 the plain box frame fits the
   scan better (p95 0.059 against 0.067), so it was kept. In both cases the frame was chosen by fit to the scan, not
   by any knowledge of the reference.
2. Scanner holes. The patched holes in the scan show up as fill triangles more than 20 times the median area.
   They are left out of every fit.
3. Planes. Take the triangles whose normal is within 11° of an axis. Their area-weighted histogram gives the
   plane candidates. Each plane position is the area-weighted mean of the triangle centres within ±0.3 mm of a peak.
4. Holes and turned diameters. Geometric least-squares circle (soft-L1 loss) on the vertices of each bore,
   taken away from its ends. For scan_3 the axis is a 3D least-squares cylinder fit on the Ø20 section. Each
   diameter is then the mean radial distance of the vertices in its section, keeping only vertices whose normal is
   radial.
5. Rounds. Fit with the arc tangent to both neighbouring planes, so only R is free. A free circle fit gave
   R 2.87 for the inner bend with an off-tangent centre. The tangent fit gives 2.99 with an RMS residual of 0.030 mm,
   equal to the scan noise. Grid search on the RMS, because a bounded optimiser stopped at a false minimum on one corner.
6. Chamfers (scan_2). The 45° plane u + w = c is fitted on the points between depth and inset 0.2–1.4 mm, one
   per edge. The four legs are 1.51 / 1.49 / 1.50 / 1.57 mm; their mean is 1.50.
7. Pocket corners (scan_2). A tangent-round fit found fewer than 10 points per corner, so the corners are taken
   as sharp. A cutter radius smaller than about 0.5 mm cannot be resolved at this noise.

## Declared dimensions

- scan_1:
  - Plate thickness is the upright, along x. The base is 4.001 and is not declared.
  - The distance is centre to centre in the plane of the upright.
- scan_2:
  - Distance = pocket centre to the hole axis, projected on the long (60 mm) axis, as the brief owner defined it.
  - Depths are measured from the machined face to the pocket floor and to the counterbore floor.
- All parts: bbox is measured in the part's own axes. `verify.py` checks it against the oriented bounding box of
  the STEP itself. Use the STL, not B-rep vertices, which miss the cylinder silhouettes.

## Files

| File | What |
|---|---|
| `scan_i.step` | the rebuilt solid (one solid, 100 % planes and cylinders) |
| `scan_i_peer.json` | the declared dimensions |
| `scan_i.py` | measure and rebuild one part; prints the fits |
| `meas.py` | frame, plane, circle, cylinder and tangent-round fits |
| `verify.py` → `report.json` | one solid, analytic area ≥ 90 %, declared diameters and bbox found on the STEP, scan-to-model p95 ≤ 0.08 |
| `img/deviation_i.png` | scan points coloured by distance to the model, with a histogram |
| `img/scan_i_iso*.png`, `slices.py`, `view_scan.py`, `inspect_scan.py` | how the parts were read |

## What could not be checked here

- Deviation to the hidden reference, hole-centre and depth errors against it, and volume error. Only the grader has
  the reference.
- Small edge breaks or cutter radii below about 0.5 mm. They are within the noise and are modelled as sharp.
