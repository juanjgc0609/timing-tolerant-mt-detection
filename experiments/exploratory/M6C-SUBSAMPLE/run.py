"""M6C-SUBSAMPLE: gate G-M6c of the v2 specification (step 1 of its build order, plan T4.1).

Does the sub-sample deflation of E62 (ttmt.detector.deflate, deflation="fractional") leave at most 5 % of an arrival
after the greedy loop deflates it, in noise-free synthetics, where the grid deflation leaves a median of 16 %
(Appendix AB, Proposition B.6; diag_greedy_out.txt)?

Noise-free records of one or two double-couple sources, continuous positions and origin times (never on a node or on
the sample grid), block statics uniform in +-0.15 s, the 8-station ring, the registered 2 Hz kernel. The loop is the
package MT-FW loop (``greedy_loop``), first step only; the fitted atom is the least-squares tensor of the read vectors
at the selected node. Three deflations of that atom are compared:
  grid        the registered convention (the response centred on the read sample);
  fractional  E62 correction (the response delayed by the fitted sub-sample delay, amplitude of the fitted tensor);
  refit       M6.3 of the v2 specification: the delayed response with its amplitude re-fitted per block.
Two v2 loops (``ttmt.v2.loop.v2_loop``) are measured on the same scenes, with a synthetic picker (a Gaussian bump of
width 0.04 s at onset = peak - 0.25 s, onset error N(0, 0.02 s), no false maxima) and the lag known:
  v2-onset    M1 onset reading, amplitude at onset + lag, deflation at the read time (``deflate_at``);
  v2-joint    the same with M4 joint reading of the 20 best hypotheses (tolerance sigma = 0.15 s).
Three read-outs per block with ||g_ib|| >= 0.1 max_b ||g_ib|| of the found source:
  peak residual   max over +-0.5 s of the true arrival of ||residual - other source|| / max of ||block alone||
                  (the G-M6c quantity, the one of Proposition B.6 and diag_greedy);
  energy left     energy of (residual - other source) in +-0.5 s / energy of the block alone there;
  other removed   (two sources) fraction of the energy of the OTHER source, in +-0.5 s of its arrivals, that the
                  deflation took away: 1 - <R, other> / ||other||^2 there (0 = untouched).
A "found" source is the first atom within 5 km and 1 s of a true source (LOC+T would use 0.5 s; here the loop step,
not the hit rule, is measured).

    python experiments/exploratory/M6C-SUBSAMPLE/run.py > experiments/exploratory/M6C-SUBSAMPLE/run_out.txt
"""
import sys

import numpy as np

from ttmt.detector import Engine, ProjectedEnergy, NOCAP, ACAP, greedy_loop
from ttmt.kernel import LEGACY as K, ricker
from ttmt.radiation import node_operators, radiation_operator, _rand_dc
from ttmt.station_terms import OnsetDelay, PeakDelay
from ttmt.traveltime import homogeneous
from ttmt.v2.loop import v2_loop
from ttmt.v2.m1 import OnsetEngine, OnsetReader

SEED = 20261108
VP, VS = 6.0, 6.0 / 1.73
NT, DT = 640, K.dt
NREP = int(sys.argv[1]) if len(sys.argv) > 1 else 60
DEFL = ("grid", "fractional", "refit")
V2 = ("v2-onset", "v2-joint")
LAG = 0.25                 # s, onset-to-peak lag of the synthetic arrivals (R1 measured 0.24-0.27 s); known to the loop
PICK_SD = 0.02             # s, Gaussian error of the synthetic picker's onset maxima
Q_SD = 0.04                # s, width of the synthetic onset-probability bump
DTAUS = (None, 3.0, 1.0, 0.5)          # origin-time separation of a second source (s); None: one source


def network():
    a = np.arange(8) * 2 * np.pi / 8 + 0.2
    return np.stack([20 * np.cos(a), 20 * np.sin(a), np.zeros(8)], 1)


def nodes(z=8.0):
    ax = np.arange(-15.0, 15.0 + 1e-9, 1.5) + 0.37
    return np.array([[x, y, z] for x in ax for y in ax])


def field(STA, src, tau, m, stat):
    """raw 3C traces of one source (P and S Ricker pulses at continuous times) and its matched-filter output."""
    g = (radiation_operator(src, STA) @ m).reshape(2, len(STA), 3)
    d = np.linalg.norm(STA - src, axis=1)
    tarr = np.stack([d / VP, d / VS]) + tau + stat                  # (2, ns)
    t = np.arange(NT) * DT
    D = np.zeros((len(STA), 3, NT))
    for b in range(2):
        for i in range(len(STA)):
            D[i] += g[b, i][:, None] * ricker(t - tarr[b, i], K.pulse_fc)[None]
    MF = np.array([[np.convolve(D[i, c], K.taps[::-1], mode="same") for c in range(3)] for i in range(len(STA))])
    return MF, g, tarr


def measure(R, other, MFk, g, tarr):
    pk, en = [], []
    gmax = np.linalg.norm(g, axis=-1).max()
    t = np.arange(NT) * DT
    for b in range(2):
        for i in range(g.shape[1]):
            if np.linalg.norm(g[b, i]) < 0.1 * gmax:
                continue
            w = np.abs(t - tarr[b, i]) <= 0.5
            # the block alone: the field of this source with only this block's pulse
            alone = MFk[b][i][:, w]
            res = (R[i] - other[i])[:, w]
            ref = np.linalg.norm(alone, axis=0).max()
            pk.append(np.linalg.norm(res, axis=0).max() / ref)
            en.append((res ** 2).sum() / (alone ** 2).sum())
    return pk, en


def other_removed(R, other, tarr_o, g_o):
    """fraction of the other source's energy, in +-0.5 s of its strong arrivals, removed by the deflation."""
    t = np.arange(NT) * DT
    gmax = np.linalg.norm(g_o, axis=-1).max()
    num = den = 0.0
    for b in range(2):
        for i in range(g_o.shape[1]):
            if np.linalg.norm(g_o[b, i]) < 0.1 * gmax:
                continue
            w = np.abs(t - tarr_o[b, i]) <= 0.5
            o = other[i][:, w]
            kept = R[i][:, w]                          # record left after the deflation
            num += float((kept * o).sum())
            den += float((o * o).sum())
    return 1.0 - num / den


def per_block_fields(STA, src, tau, m, stat):
    """matched-filter output of each phase of one source separately (P only, S only)."""
    out = []
    for b in range(2):
        g = (radiation_operator(src, STA) @ m).reshape(2, len(STA), 3)
        d = np.linalg.norm(STA - src, axis=1)
        tarr = np.stack([d / VP, d / VS]) + tau + stat
        t = np.arange(NT) * DT
        D = np.zeros((len(STA), 3, NT))
        for i in range(len(STA)):
            D[i] += g[b, i][:, None] * ricker(t - tarr[b, i], K.pulse_fc)[None]
        out.append(np.array([[np.convolve(D[i, c], K.taps[::-1], mode="same") for c in range(3)]
                             for i in range(len(STA))]))
    return out


def synth_q(srcs_f, ns, r):
    """onset evidence of the synthetic picker for the arrivals (peak times tarr (2, ns)) of every source."""
    t = np.arange(NT) * DT
    q = np.zeros((ns, 2, NT))
    for tarr in srcs_f:
        on = tarr - LAG + r.normal(0, PICK_SD, tarr.shape)
        for b in range(2):
            for i in range(ns):
                q[i, b] = np.maximum(q[i, b], 0.9 * np.exp(-0.5 * ((t - on[b, i]) / Q_SD) ** 2))
    return q


def main():
    rng = np.random.default_rng(SEED)
    STA = network()
    XN = nodes()
    taus = np.round(np.arange(4.0, 16.0 + 1e-9, 0.1), 10)
    tt, _ = homogeneous(XN, STA, VP, VS)
    Gs, Qst = node_operators(XN, STA)
    E = Engine(XN, STA, tt, taus, PeakDelay.zeros(len(STA)), Gs, Qst, rounding="nearest")
    E2 = OnsetEngine(XN, STA, tt, taus, OnsetDelay(np.full((len(STA), 2), -LAG)), Gs, Qst, rounding="nearest")
    lag = np.full((len(STA), 2), LAG / DT)
    stats = {"MT-FW": ProjectedEnergy(NOCAP), "ACAP": ACAP(3.0)}
    print(f"M6C-SUBSAMPLE  seed {SEED}  {NREP} replicates per configuration  numpy {np.__version__}")
    print("8-station ring r = 20 km, nodes 21 x 21 at 1.5 km (offset 0.37 km), z = 8 km, origins 4-16 s step 0.1 s,")
    print("rounding 'nearest', 2 Hz kernel, block statics U(-0.15, 0.15) s, noise-free. First loop step only.")
    print("Peak residual = max ||residual|| in +-0.5 s / max ||block alone||; G-M6c passes if its median <= 0.05.\n")
    hdr = (f"{'statistic':9s} {'dtau':>5s} {'deflation':>10s} {'found':>6s} {'blocks':>6s} "
           f"{'pk med':>7s} {'pk p90':>7s} {'pk max':>7s} {'E med':>7s} {'E p90':>7s} {'oth med':>8s} {'oth p90':>8s}")
    print(hdr)
    for sname, stat_ in stats.items():
        for dtau in DTAUS:
            res = {k: ([], [], []) for k in DEFL + V2}
            found = {k: 0 for k in DEFL + V2}
            r2 = np.random.default_rng([SEED, int(10 * (dtau or 0)), len(sname)])
            for rep in range(NREP):
                srcs = []
                n_src = 1 if dtau is None else 2
                for k in range(n_src):
                    while True:
                        x = np.r_[r2.uniform(-10, 10, 2), 8.0]
                        if all(np.linalg.norm(x[:2] - s[0][:2]) >= 6.0 for s in srcs):
                            break
                    tau = 6.0 + r2.uniform(0, DT) + (0.0 if k == 0 else dtau)
                    m = _rand_dc(r2).ravel()[[0, 4, 8, 1, 2, 5]] * 2000.0
                    stat = r2.uniform(-0.15, 0.15, (2, len(STA)))
                    srcs.append((x, tau, m, stat))
                MF = [field(STA, *s) for s in srcs]
                Y = sum(f[0] for f in MF)
                q = synth_q([f[2] for f in MF], len(STA), r2)
                for dfl in DEFL + V2:
                    if dfl in DEFL:
                        atoms, R = greedy_loop(E, Y, 1, stat_, deflation=dfl, return_residual=True)
                    else:
                        jt = {} if dfl == "v2-joint" else None
                        atoms, R = v2_loop(E2, Y, q, 1, lag, statistic=stat_, reader=OnsetReader(4), joint=jt,
                                           n_joint=20, sig=0.15 / DT, return_residual=True)
                    a = atoms[0]
                    k = next((k for k, s in enumerate(srcs) if np.linalg.norm(a["x"][:2] - s[0][:2]) <= 5.0
                              and abs(a["tau"] - s[1]) <= 1.0), None)
                    if k is None:
                        continue
                    found[dfl] += 1
                    other = sum((MF[j][0] for j in range(len(srcs)) if j != k), np.zeros_like(Y))
                    MFk = per_block_fields(STA, *srcs[k])
                    pk, en = measure(R, other, MFk, MF[k][1], MF[k][2])
                    res[dfl][0].extend(pk)
                    res[dfl][1].extend(en)
                    if len(srcs) == 2:
                        j = 1 - k
                        res[dfl][2].append(other_removed(R, MF[j][0], MF[j][2], MF[j][1]))
            for dfl, (pk, en, orm) in res.items():
                pk, en = np.array(pk), np.array(en)
                d = "single" if dtau is None else f"{dtau:.1f}"
                o = f"{np.median(orm):8.3f} {np.quantile(orm, 0.9):8.3f}" if orm else f"{'-':>8s} {'-':>8s}"
                print(f"{sname:9s} {d:>5s} {dfl:>10s} {found[dfl]:6d} {len(pk):6d} {np.median(pk):7.3f} "
                      f"{np.quantile(pk, 0.9):7.3f} {pk.max():7.3f} {np.median(en):7.3f} {np.quantile(en, 0.9):7.3f} {o}")
        print()


if __name__ == "__main__":
    sys.exit(main())
