# Scan-to-STEP brief

Three simulated scans of machined parts, `scan_1.stl` to `scan_3.stl` (units mm). Each scan has about 0.03 mm of measurement noise, an uneven mesh, one or two holes in the surface where the scanner saw nothing, and an arbitrary position and rotation. There is no reference coordinate system.

## Task, for each scan
1. Rebuild the part as a parametric solid (planes, cylinders, cones, fillets, chamfers, holes). Export `scan_i.step`. A meshed or faceted solid is rejected.
2. Declare the key dimensions you measured in `scan_i_peer.json`: `{"dims": [{"type": "hole_diameter", "value": 6.0}, ...]}`. Units mm. `bbox` takes a list of the three overall sizes in ascending order.

Declare exactly these dimension types (count in brackets). Extra or missing entries fail:
- scan_1: hole_diameter (2), distance between the two hole centres (1), plate thickness (1), bbox (1)
- scan_2: hole_diameter (2), depth (2: pocket depth, counterbore depth), distance (1), bbox (1)
- scan_3: diameter (4: each turned diameter including the groove bottom), width (1: groove width), bbox (1)

## Scoring
The rebuilt solid is aligned to the hidden reference model and scored on:
- surface deviation to the reference: p95 at most 0.06 mm, max at most 0.2 mm; scan to model p95 at most 0.08 mm
- at least 90% of the surface area made of analytic faces
- hole centres within 0.1 mm, depths within 0.1 mm, volume within 1%
- every declared dimension within 0.05 mm of the reference and of the value measured on your own STEP
