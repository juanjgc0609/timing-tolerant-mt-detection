"""diag_projection.py - equirectangular local frame (r0_lib.to_xy) vs WGS84 geodesic: station-event and event-event
epicentral distances over the R0 stations and EVAL M>=2 events."""
import csv
import numpy as np
from pyproj import Geod
import verif_analysis as V
cat, (lat0, lon0) = V.catalogue()
rows = {r["station"]: r for r in csv.DictReader(open(V.DATA + "/stations.csv")) if r["channel"] == "HHZ"}
names = sorted({n for s in ("CLC SRT CCC TOW2 LRL WRC2 WBM SLA JRC2 WCS2 MPM".split()) for n in [s]})
g = Geod(ellps="WGS84")
raw = list(csv.DictReader(open(V.DATA + "/catalog.csv")))
k = 111.19
err = []
for i, e in enumerate(raw):
    if not (900 <= cat[i, 2] < 3560 and cat[i, 3] >= 2):
        continue
    for n in names:
        la, lo = float(rows[n]["latitude"]), float(rows[n]["longitude"])
        x, y = (lo - lon0) * k * np.cos(np.radians(lat0)), (la - lat0) * k
        dxy = np.hypot(x - cat[i, 0], y - cat[i, 1])
        _, _, dg = g.inv(float(e["longitude"]), float(e["latitude"]), lo, la)
        err.append(dxy - dg / 1000.0)
err = np.array(err)
print(f"station-event epicentral distance, frame minus geodesic: mean {err.mean():+.3f} km, max |.| {np.abs(err).max():.3f} km "
      f"(-> max travel-time error {np.abs(err).max()/ (6/1.73):.3f} s for S)")
