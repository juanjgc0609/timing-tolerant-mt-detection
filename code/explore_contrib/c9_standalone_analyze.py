"""C9 P1-P3 analysis: epicentral error of each travel-time-free estimator on single events; bootstrap over networks."""
import json
import sys

import numpy as np

rows = json.load(open("out/c9_standalone.json"))
CELLS = ("ladder_4", "ladder_6", "unilateral_c40")
EST = ["P1-azim", "P1-grid", "P1-grid-gradmodel", "P3-Strans", "P13-pol", "P3-gainfree", "P2-mt", "P2-mt-cap6",
       "P2b-amp", "TIME-raw", "TIME-truestatics"]
rng = np.random.default_rng(5)


def err(r, mode, e):
    xy = r[mode][e]
    v = float(np.hypot(xy[0] - r["src"][0], xy[1] - r["src"][1]))
    return v if np.isfinite(v) else 99.0


def stat(sub, mode, e, B=500):
    nets = sorted({r["net"] for r in sub})
    by = {n: np.array([err(r, mode, e) for r in sub if r["net"] == n]) for n in nets}
    allv = np.concatenate([by[n] for n in nets])
    med, p5 = np.median(allv), np.mean(allv <= 5.0)
    bm, bp = [], []
    for _ in range(B):
        v = np.concatenate([by[nets[k]] for k in rng.integers(0, len(nets), len(nets))])
        bm.append(np.median(v)); bp.append(np.mean(v <= 5.0))
    return med, np.percentile(bm, [2.5, 97.5]), p5, np.percentile(bp, [2.5, 97.5])


def table(title, sel, modes=("true", "blind"), ests=EST):
    print("\n" + title)
    print("  median epicentral error km [95% network bootstrap] | fraction within 5 km [95%]")
    for mode in modes:
        for e in ests:
            line = f"  {mode:5s} {e:18s}"
            for c in CELLS:
                sub = [r for r in rows if r["cell"] == c and sel(r)]
                m, cm, p, cp = stat(sub, mode, e)
                line += f" | {m:5.1f} [{cm[0]:4.1f},{cm[1]:5.1f}] {p:.2f} [{cp[0]:.2f},{cp[1]:.2f}]"
            print(line)


n = len([r for r in rows if r["dir_deg"] == 10 and r["site_sd"] == 0.5])
print(f"C9 standalone: {n} single events (20 networks x 10 events x 3 cells), thin clutter, each under 7 generator settings")
print("columns: ladder_4 (4 stations) | ladder_6 (6) | unilateral_c40 (8, one-sided)")
src = np.array([r["src"] for r in rows if r["dir_deg"] == 10 and r["site_sd"] == 0.5])
print(f"reference: guessing the centre of the source region gives median {np.median(np.hypot(src[:, 0], src[:, 1])):.1f} km, "
      f"within 5 km {np.mean(np.hypot(src[:, 0], src[:, 1]) <= 5):.2f}")
for c in CELLS:
    sub = [r for r in rows if r["cell"] == c and r["dir_deg"] == 10 and r["site_sd"] == 0.5]
    aP = np.concatenate([r["ampP"] for r in sub]); aS = np.concatenate([r["ampS"] for r in sub])
    bo = np.mean([r["blind_ok"] for r in sub], axis=0)
    print(f"  {c:15s} nominal: |v_P| median {np.median(aP):.1f} noise sd (P25 {np.percentile(aP, 25):.1f}), |v_S| median "
          f"{np.median(aS):.1f}; blind reader within 3 samples of the true arrival: P {bo[0]:.2f}, S {bo[1]:.2f}")
table("NOMINAL generator (polarization deviation 10 deg, amplitude-factor sd 0.5)",
      lambda r: r["dir_deg"] == 10 and r["site_sd"] == 0.5)
for dd in (0.0, 5.0, 20.0):
    table(f"P1/P3 versus polarization deviation: {dd:g} deg (amplitude-factor sd 0.5)",
          lambda r, dd=dd: r["dir_deg"] == dd and r["site_sd"] == 0.5, modes=("true",),
          ests=["P1-azim", "P1-grid", "P3-Strans", "P13-pol", "P3-gainfree", "P2-mt"])
for ss in (0.0, 0.25):
    table(f"P2 versus amplitude-factor sd: {ss:g} (polarization deviation 10 deg)",
          lambda r, ss=ss: r["dir_deg"] == 10 and r["site_sd"] == ss, modes=("true",),
          ests=["P2-mt", "P2-mt-cap6", "P2b-amp", "P3-gainfree", "P13-pol"])
table("ORACLE generator: polarization deviation 0, amplitude-factor sd 0 (waveform noise, coda, clutter, pulse and "
      "medium mismatch remain)", lambda r: r["dir_deg"] == 0 and r["site_sd"] == 0)
