"""QC of the NN StationXML (metadata only) with ttmt.qc.station_qc; adds the distance to each disc centre.

Run on the laptop 2026-11-15 (output SHA-256 1e776ecb972e00cadaeb8b59be94ec191cbed1ae8d39b6039f8d5e213855ed3a):
Antelope Valley 10 of 11 NN stations usable (HH at WAK, BFC, YER, PNT; the rest HN), Bodie 6 of 14 (8 are EHZ only)."""
import json, math, os, sys
from obspy import read_inventory, UTCDateTime
from ttmt import qc

META = r"C:\Users\juanj\ttmt\real_meta"
DISCS = {
    "conf_antelope_valley_2021": (38.5, -119.5, "2021-07-08T22:00:00Z", "2021-08-07T22:00:00Z"),
    "conf_bodie_2020": (38.3, -118.9, "2020-04-11T14:00:00Z", "2020-05-11T14:00:00Z"),
}


def km(a, b, c, d):
    p1, p2 = math.radians(a), math.radians(c)
    h = math.sin((p2 - p1) / 2) ** 2 + math.cos(p1) * math.cos(p2) * math.sin(math.radians(d - b) / 2) ** 2
    return 2 * 6371.0 * math.asin(math.sqrt(h))


out_all = {"ttmt_qc_file": qc.__file__}
for sid, (lat0, lon0, t0, t1) in DISCS.items():
    path = os.path.join(META, "stationxml", f"NN_{sid}.xml")
    inv = read_inventory(path)
    res = qc.report_dict(qc.station_qc(inv, UTCDateTime(t0), UTCDateTime(t1)))
    for net in inv:
        for sta in net:
            key = f"{net.code}.{sta.code}"
            if key in res:
                res[key]["dist_km"] = round(km(lat0, lon0, sta.latitude, sta.longitude), 2)
                res[key]["lat"], res[key]["lon"] = sta.latitude, sta.longitude
    out_all[sid] = {"file": path, "t0": t0, "t1": t1, "centre": [lat0, lon0], "stations": res}
    print(f"== {sid} ({len(res)} stations, usable {sum(v['usable'] for v in res.values())})")
    for k in sorted(res, key=lambda k: res[k]["dist_km"]):
        v = res[k]
        ch_iss = {c: i for c, i in v["channels"].items() if i}
        print(f"{k} {v['dist_km']} km usable={v['usable']} 3C={v['three_component']} issues={v['issues']} "
              f"channels={len(v['channels'])} channel_issues={ch_iss}")
os.makedirs(os.path.join(META, "qc"), exist_ok=True)
dst = os.path.join(META, "qc", "nn_station_qc.json")
with open(dst, "w", encoding="utf-8", newline="\n") as f:
    json.dump(out_all, f, indent=1, default=str)
print("written", dst)
