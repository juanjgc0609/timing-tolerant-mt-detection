"""merge run/parts/det_eval_<sub>_<methods>.json into run/det_eval_<sub>.json (atoms and timings per method; the
calibration products are identical across parts and are checked equal).  N8 was written whole before the restart."""
import glob, json, os
for sub in ("N6", "N4", "U8"):
    fs = sorted(glob.glob(f"run/parts/det_eval_{sub}_*.json"))
    if not fs:
        continue
    out = None
    for f in fs:
        r = json.load(open(f))
        if out is None:
            out = r
        else:
            for k in ("corr_sc", "corr_s2", "scale", "STA"):
                assert r[k] == out[k], (sub, k)
            out["atoms"].update(r["atoms"])
            out["sec"].update(r["sec"])
    json.dump(out, open(f"run/det_eval_{sub}.json", "w"))
    print(sub, sorted(out["atoms"]))
