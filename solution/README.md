# Solution: scan to parametric STEP

The three scans are simulated. The brief's author generated them from CAD parts, added about 0.03 mm of Gaussian
noise, punched one or two holes, and rotated them at random. They are not real scanner data.

Each part is rebuilt in build123d from planes and cylinders only, with no mesh or faceted face. Each model is put back
in the scan's own position and rotation, so the scan-to-model deviation can be measured directly.

| Part | Rebuilt as | Scan → model (mm) |
|---|---|---|
| scan_1 | L-bracket 29.99 × 40.00 × 60.00. Base and upright 4.00 thick. Inner bend round R 2.99, top corners R 4.99. Through holes Ø5.998 and Ø8.000 in the upright, 20.001 apart | p95 0.059, max 0.133 |
| scan_2 | Block 30.00 × 40.00 × 60.00. 45° chamfer 1.50 on the four edges of the machined face. Pocket 20.00 × 24.00 × 8.00 deep. Counterbore Ø10.999 × 6.50 deep, coaxial with a Ø6.600 through hole | p95 0.059, max 0.126 |
| scan_3 | Turned shaft 75.00 long: Ø13.999 × 20.0, Ø19.996 (split by a 2.003 wide groove, groove bottom Ø16.990), Ø27.991 × 25.0 | p95 0.059, max 0.134 |

With 0.03 mm noise, a p95 of 0.059 mm is about 2σ. That shows the fit is consistent with the noise. The error against
the true part is known only from the scoring below.

Every declared value is the measured value, not rounded to a nominal size, and the model is built from the same
numbers. `verify.py` checks the declared diameters and bbox against the STEP. It is a consistency check, not an
accuracy check.

## How it was scored

The solution was built blind: I never saw the reference models. The brief's author scored it once against the hidden
references, and all three parts passed.

- Hausdorff distance to the reference was 0.0069 / 0.0023 / 0.0051 mm (max).
- Declared dimensions were within 0.01 mm.
- Volume was within 0.05 %, and 100 % of the faces were analytic.

These figures come from simulated scans with known noise. They show that the method works on clean synthetic data, not
that a real scanner would reach 0.01 mm.

## Method

1. Frame.
   - Each scan is moved into its oriented bounding box (trimesh).
   - For scan_2 and scan_3, the box frame is then re-aligned on the measured normals of the flat faces. On scan_2
     the minimum-volume box was 0.8° off because of the chamfers and the pocket, which left scan points 0.4 mm
     outside the model.
   - For scan_1 the plain box frame fits the scan better (p95 0.059 against 0.067), so it was kept.
   - In every case the frame was chosen by its fit to the scan, never by any knowledge of the reference.
2. Scanner holes. The patched holes show up as fill triangles more than 20 times the median area. They are left
   out of every fit.
3. Planes.
   - I read the rough position of each plane off section plots (`tools/slices.py`) and the per-axis plane list
     from `tools/inspect_scan.py`.
   - The measured position is the area-weighted mean of the centres of the triangles that face along that axis
     (normal within 11°), within ±0.3 mm of the rough value.
4. Holes and turned diameters.
   - Holes: a geometric least-squares circle (soft-L1 loss) on the vertices of each bore, taken away from its ends.
   - scan_3: the axis comes from a 3D least-squares cylinder fit on the Ø20 section. Each diameter is then the mean
     radial distance of the vertices in its section, keeping only vertices whose normal is radial.
5. Rounds.
   - Each round is fitted with the arc tangent to both neighbouring planes, so only R is free.
   - A free circle fit gave R 2.87 for the inner bend, with a centre that was not tangent. The tangent fit gives 2.99,
     with an RMS residual of 0.030 mm, which equals the scan noise.
   - R is found by grid search on the RMS, because a bounded optimiser stopped at a false minimum on one corner.
6. Chamfers (scan_2). The 45° plane u + w = c is fitted on the points between 0.2 and 1.4 mm in depth and inset,
   one fit per edge. The four legs are 1.51, 1.49, 1.50 and 1.57 mm, with a mean of 1.50.
7. Pocket corners (scan_2).
   - The tangent-round fit finds about R 1.50 at all four corners, but on only 21 points per corner. The code needs
     100 points before it accepts a round, so the corners are modelled sharp.
   - That choice passed the scoring. Whether the reference corners are sharp or have a small cutter radius is not
     shown.

## Declared dimensions

For scan_1, the plate thickness is the upright's, measured along x. The base is 4.001 thick and is not declared. The
hole distance is measured centre to centre in the plane of the upright.

For scan_2, the brief's author clarified the definition of the distance before I solved the part: it runs from the
pocket centre to the hole axis, projected on the long (60 mm) axis. The depths are measured from the machined face to
the pocket floor and to the counterbore floor.

For every part, the bbox is measured in the part's own axes. `verify.py` compares it with the oriented bounding box of
the STEP's own STL; B-rep vertices are not used, because they miss the cylinder silhouettes.

## Files

| File | What |
|---|---|
| `scan_i.step` | the rebuilt solid (one solid, planes and cylinders only), exactly as scored |
| `scan_i_peer.json` | the declared dimensions |
| `scan_i.py` | measures and rebuilds one part, and prints the fits |
| `meas.py` | frame, plane, circle, cylinder and tangent-round fits |
| `verify.py` → `report.json` | one solid, analytic area ≥ 90 %, declared diameters and bbox found on the STEP, scan-to-model p95 ≤ 0.08 |
| `img/deviation_i.png` | scan points coloured by distance to the model, two views and a histogram |
| `img/scan_i_iso*.png` | the raw scans; the scanner holes appear filled by long triangles from the mesher |
| `tools/` | `inspect_scan.py`, `slices.py`, `view_scan.py`: how the parts were read |
| `requirements.txt` | Python packages |

`scan_1.py` keeps its original plane seeds (14.99 / 10.99). Rounding them to 15.0 / 11.0 gives the same declared values
but moves one STEP coordinate by 6e-6 mm, and the scored file should stay reproducible.
