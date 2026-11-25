"""Stop condition of EQNet2 stage 1 (laptop runbook, step 1; rule 40), declared before the run: for each set, the median
over records of the largest event score (over stations and samples) within +-3 s of a true origin must exceed the
median over records of the largest score in 10 s windows more than 30 s from every true origin. If it does not in
some set, stage 2 does not start and the result is reported.

    python eqnet2_signal_check.py C:/Users/juanj/ttmt/tune/*/tune_eqnet2.json [--out check.json]

Reads the score files written by ``eqnet2_stage1.py scores`` and the truth of the 'tune' partition (SetReader, as
stage 1). It is a check of the inputs of the tuning, not a tuning trial: it reads no detection and no objective.
"""
import argparse
import json
import os
import sys

import numpy as np

NEAR_S = 3.0        # half width around a true origin
FAR_S = 30.0        # minimum distance of a noise window from every true origin
WIN_S = 10.0        # noise window length


def record_maxima(score, dt, t0, origins, near=NEAR_S, far=FAR_S, win=WIN_S):
    """(max score within +-near of any origin, max score over the 10 s windows at least far from every origin); nan
    when the record has no such sample. score: (n_sta, m) on times t0 + i * dt."""
    s = np.asarray(score, float).max(axis=0)
    t = t0 + dt * np.arange(s.size)
    o = np.asarray([x for x in origins if np.isfinite(x)], float)
    if o.size:
        d = np.abs(t[:, None] - o[None, :]).min(axis=1)
    else:
        d = np.full(t.size, np.inf)
    hit = s[d <= near]
    near_max = float(hit.max()) if hit.size else float("nan")
    w = max(1, int(round(win / dt)))
    noise = []
    for a in range(0, s.size - w + 1, w):                     # non-overlapping windows, all samples far
        if (d[a:a + w] > far).all():
            noise.append(s[a:a + w].max())
    noise_max = float(max(noise)) if noise else float("nan")
    return near_max, noise_max


def decide(pairs):
    """median of each column over records (nan ignored) and whether the near median exceeds the noise median."""
    a = np.asarray(pairs, float).reshape(-1, 2)
    near = float(np.nanmedian(a[:, 0])) if np.isfinite(a[:, 0]).any() else float("nan")
    noise = float(np.nanmedian(a[:, 1])) if np.isfinite(a[:, 1]).any() else float("nan")
    return dict(n_records=int(len(a)), n_near=int(np.isfinite(a[:, 0]).sum()), n_noise=int(np.isfinite(a[:, 1]).sum()),
                median_near=near, median_noise=noise, passes=bool(np.isfinite(near) and np.isfinite(noise) and near > noise))


def check_manifest(path, reader_of=None):
    man = json.load(open(path))
    if reader_of is None:
        from ttmt.synth.sets import SetReader

        def reader_of(set_id):
            return SetReader(set_id, "tune", verify="fingerprint")
    reader = reader_of(man["set_id"])
    if reader.partition != "tune":
        raise SystemExit("the check reads the tuning partition only")
    root = os.path.dirname(os.path.abspath(path))
    pairs = []
    for lab in sorted({r["label"] for r in man["records"]}):
        name = {r["label"]: r["eqnet2"] for r in man["records"]}[lab]
        z = np.load(os.path.join(root, name))
        tr = reader.truth(lab)
        pairs.append(record_maxima(z["score"], float(z["dt"]), float(z["t0"]), np.asarray(tr.t, float)))
    return dict(set_id=man["set_id"], **decide(pairs))


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("manifests", nargs="+", help="tune_eqnet2.json of each set")
    ap.add_argument("--out", default=None)
    a = ap.parse_args(argv)
    res = [check_manifest(m) for m in a.manifests]
    print(json.dumps(res, indent=1))
    if a.out:
        with open(a.out, "w") as fh:
            json.dump(res, fh, indent=1)
    return 0 if all(r["passes"] for r in res) else 1


if __name__ == "__main__":
    sys.exit(main())
