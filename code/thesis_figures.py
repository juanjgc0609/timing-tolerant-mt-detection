"""
thesis_figures.py - figures of docs/thesis/Thesis.md from the stored outputs of Task R1 (no new analysis):
  fig1_r1_geometry.png   stations and truth events of the Hi-net swarm (frame of R1)
  fig2_r1_froc.png       recall of the 54 EVAL targets (LOC+T) against false detections per hour, from hyb_eval.json
  fig3_r1_ghost.png      matched-filter energy at three stations around a CAL event (M 2.6), with the predicted arrivals
                         of the event and of the highest-scoring pick-free atom that followed it (a coda ghost)
    python code/thesis_figures.py RUN_DIR OUT_DIR       (RUN_DIR = results/taskR1/run)
"""
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

sys.path.insert(0, "/home/claude/work/r1")
import r1_lib as R  # noqa: E402

INK, INK2, SURF = "#0b0b0b", "#52514e", "#fcfcfb"
C = {"HYB": "#2a78d6", "BASE": "#eb6834", "BESTPICK": "#1baf7a", "P1": "#eda100"}
LS = {"HYB": "-", "BASE": "--", "BESTPICK": "-.", "P1": ":"}
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK,
                     "xtick.color": INK2, "ytick.color": INK2, "axes.spines.top": False, "axes.spines.right": False,
                     "figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF})


def froc_curve(S, tg, M):
    """recall against the number of false detections allowed (exact, ties never split)."""
    lab = S.lab["LOC+T"]
    order = np.argsort(-S.s, kind="stable")
    fal, got, xs, ys = 0, set(), [0], [0.0]
    tgs = set(int(i) for i in tg)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and S.s[order[j + 1]] == S.s[order[i]]:
            j += 1
        for k in order[i:j + 1]:
            if lab[k] < 0:
                fal += 1
            elif M[lab[k]] >= R.M_TARGET:
                got.add(int(lab[k]))
        xs.append(fal)
        ys.append(len(got & tgs) / len(tgs))
        i = j + 1
    xs, ys = np.array(xs), np.array(ys)
    best = np.array([ys[xs <= f].max() for f in range(int(xs.max()) + 1)])
    return best


def fig1(out):
    c = R.frame()
    cat = R.catalogue_local(c)
    STA, _ = R.stations(R.SUBNETS["H7"], c)
    fig, ax = plt.subplots(figsize=(4.6, 4.2))
    ax.scatter(STA[:, 0], STA[:, 1], marker="^", s=70, color=INK, zorder=3)
    for (x, y, _z), n in zip(STA, R.SUBNETS["H7"]):
        ax.annotate(n, (x, y), xytext=(5, 4), textcoords="offset points", color=INK, fontsize=8)
    x = np.array([e["x"] for e in cat]); y = np.array([e["y"] for e in cat]); m = np.array([e["M"] for e in cat])
    j = np.array([e["jma"] for e in cat])
    ax.scatter(x[~j], y[~j], s=6 + 10 * np.clip(m[~j], 0, None), color=C["HYB"], alpha=0.6, lw=0, label="Match and Locate detections")
    ax.scatter(x[j], y[j], s=6 + 10 * m[j], facecolor="none", edgecolor=C["BASE"], lw=1.2, label="JMA events")
    ax.set_aspect("equal"); ax.set_xlabel("east (km)"); ax.set_ylabel("north (km)")
    ax.set_title("Task R1: seven Hi-net stations and the swarm", color=INK, fontsize=10, loc="left")
    ax.legend(frameon=False, fontsize=8, loc="lower right")
    fig.tight_layout(); fig.savefig(out, dpi=200); plt.close(fig)


def fig2(run, out):
    R.install()
    import r0_analyze as A
    d = json.load(open(os.path.join(run, "hyb_eval.json")))
    lo, hi = R.T_EVAL
    hours = (hi - lo) / 3600.0
    cat = R.catalogue_local(R.frame())
    M = np.array([e["M"] for e in cat])
    tg = np.array([i for i, e in enumerate(cat) if lo <= e["t"] < hi and e["M"] >= R.M_TARGET])
    lists = {"HYB": d["HYB"], "BASE": d["BASE"], "BESTPICK": d["BESTPICK"]['pyocto_1d|{"n_picks": 3, "tol": 0.5}'],
             "P1": d["P1"]}
    names = {"HYB": "HYB (hybrid)", "BASE": "BASE (GaMMA-1D)", "BESTPICK": "BESTPICK (PyOcto-1D)", "P1": "P1 (pick-free)"}
    fig, ax = plt.subplots(figsize=(5.6, 3.6))
    for k, v in lists.items():
        b = froc_curve(A.Set(v, cat, lo, hi), tg, M)
        f = np.arange(len(b)) / hours
        ax.step(f, b, where="post", color=C[k], ls=LS[k], lw=2, label=names[k])
    ax.axvline(30, color=INK2, lw=0.8, ls=(0, (2, 2)))
    ax.text(31, 0.02, "30 per hour (primary)", color=INK2, fontsize=8)
    ax.set_xlim(0, 130); ax.set_ylim(0, 0.6)
    ax.set_xlabel("false detections per hour (EVAL, 1120 s)"); ax.set_ylabel("recall of the 54 targets (LOC+T)")
    ax.set_title("Task R1: recall against false detections", color=INK, fontsize=10, loc="left")
    ax.grid(axis="y", color="#e4e3df", lw=0.6)
    ax.legend(frameon=False, fontsize=8, loc="center right")
    fig.tight_layout(); fig.savefig(out, dpi=200); plt.close(fig)


def fig3(run, out):
    R.setup()
    import c18_tt as TT
    import c18j_layered as LJ
    from xc import matched
    net = R.make_net("cal")
    p1 = json.load(open(os.path.join(run, "p1_cal.json")))["atoms"]["NOCAP"]
    ev = [e for e in net.cat if abs(e["t"] - 145.5) < 1.0][0]
    ghosts = sorted([a for a in p1 if 0.5 < a[3] - ev["t"] < 6.0], key=lambda a: -a[4])
    g = ghosts[0]
    corr = np.array(json.load(open(os.path.join(run, "p1_cal_chunks", "meta.json")))["corr_layered"])
    t0 = ev["t"] - 6.5
    W = R.window(net.d16, t0)
    MF, sd = matched(W)
    E = np.sum((MF / sd) ** 2, axis=1)
    tt = np.arange(E.shape[1]) * R.DT + t0 - ev["t"]
    te = TT.times(LJ.MODEL, LJ.R_VPVS, [ev["x"], ev["y"], ev["z"]], net.STA) + corr.T
    tg_ = TT.times(LJ.MODEL, LJ.R_VPVS, g[:3], net.STA) + corr.T + (g[3] - ev["t"])
    sel = [6, 4, 5]
    fig, axs = plt.subplots(len(sel), 1, figsize=(6.0, 5.8), sharex=True)
    for ax, i in zip(axs, sel):
        ax.semilogy(tt, np.maximum(E[i], 1e-2), color=INK2, lw=1)
        for b, lab in ((0, "P"), (1, "S")):
            ax.axvline(te[b, i], color=C["HYB"], lw=1.8)
            ax.text(te[b, i], 3e4, lab, color=INK, fontsize=8, ha="center")
            ax.axvline(tg_[b, i], color=C["BASE"], lw=1.8, ls="--")
        ax.set_ylim(1e-1, 1e5)
        ax.text(0.01, 0.82, f"{net.names[i]}", transform=ax.transAxes, color=INK, fontsize=9)
        ax.set_ylabel("MF energy")
    axs[-1].set_xlim(-1, 14)
    axs[-1].set_xlabel(f"time after the origin of the M {ev['M']:.1f} event (s)")
    axs[0].set_title("Task R1 (CAL): where a coda ghost sits", color=INK, fontsize=10, loc="left")
    axs[0].plot([], [], color=C["HYB"], lw=1.8, label="predicted P and S of the event")
    axs[0].plot([], [], color=C["BASE"], lw=1.8, ls="--",
                label=f"predicted P and S of the top pick-free atom (+{g[3] - ev['t']:.1f} s, "
                      f"{np.hypot(g[0] - ev['x'], g[1] - ev['y']):.0f} km away)")
    fig.legend(*axs[0].get_legend_handles_labels(), frameon=False, fontsize=8, loc="lower center", ncol=1)
    fig.tight_layout(rect=(0, 0.07, 1, 1)); fig.savefig(out, dpi=200); plt.close(fig)
    return ev, g


if __name__ == "__main__":
    run, outd = sys.argv[1], sys.argv[2]
    os.makedirs(outd, exist_ok=True)
    fig1(os.path.join(outd, "fig1_r1_geometry.png"))
    fig2(run, os.path.join(outd, "fig2_r1_froc.png"))
    print(fig3(run, os.path.join(outd, "fig3_r1_ghost.png")))
