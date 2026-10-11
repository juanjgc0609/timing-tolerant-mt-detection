"""
c13_lib.py - C13: target-decoy / empirical-null control of false discoveries for the pick-free MT-FW detector (cap 6,
registered +SC terms, hard +-0.25 s peak rule).  SEED 999 ONLY.

The detector (c10_lib.mt_loop, unchanged) is run once on the real hypotheses (TARGET) and once per decoy construction
on hypotheses that cannot be a real source.  Decoys (a fresh random draw per window, from a stream of its own):

  D1p  station-permuted moveout: the travel-time table of station i is that of station pi(i), pi a random derangement;
       the station term stays with the physical station; amplitude operators (G, Q) unchanged.
  D1o  scrambled moveout by offsets: every (station, phase) column of the sample table is displaced by its own offset,
       circularly over the 40 s trace, offsets pairwise at least 1 s apart; G, Q unchanged.  The P-S pairing of a
       station is destroyed as well.
  D2s  sign-flipped templates: the 3-component block of half of the (station, phase) arrivals changes sign in G and Q;
       timing unchanged.
  D2r  polarisation-scrambled templates: an independent random rotation of the 3-component block of every station
       (both phases) in G and Q; timing unchanged.
       (Time reversal of the template is not a decoy here: the detector reads matched-filtered traces whose wavelet is
       an autocorrelation, symmetric in time.)
  D3   data surrogate: every station trace is shifted circularly by its own random lag, lags pairwise at least 3 s
       apart; hypotheses unchanged.  Equivalent to D1o with one offset per station (P-S pairing kept).

Calibration of the acceptance threshold from a batch of windows (targets and decoys pooled over the batch):
  tdc    separate search: accept targets with record score >= t*, t* the smallest t with
         (1 + #decoys >= t) / max(1, #targets >= t) <= q.
  bh     conformal p-value p = (1 + #decoys >= s) / (n_decoys + 1) of every target atom, Benjamini-Hochberg at level q
         over the m target atoms of the batch.  With n_decoys = m this is the tdc rule with the factor m / (m + 1).
  comp   rank-paired competition (knockoff+): atom r of the target list against atom r of the decoy list of the same
         window; W = max of the two record scores; accept the target-won pairs with W >= t*, t* the smallest t with
         (1 + #decoy-won pairs >= t) / max(1, #target-won pairs >= t) <= q.
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import copy
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.dirname(HERE)
for p in (HERE, WORK):
    if p not in sys.path:
        sys.path.insert(0, p)

import c10_lib as C  # noqa: E402
from task0_common import DT, NT, VP, VS  # noqa: E402

SEED = 999
NETS = tuple(range(1400, 1412))
TAG_STD, TAG_DENSE, TAG_CLUT, TAG_NOISE, TAG_DECOY, TAG_SMOKE = 1401, 1402, 1403, 1404, 1405, 1406
DECOYS = ("D1p", "D1o", "D2s", "D2r", "D3")
QS = (0.05, 0.1, 0.2, 0.3, 0.5)
GAP_D1O = int(round(1.0 / DT))
GAP_D3 = int(round(3.0 / DT))


# ============================================================================================ decoy draws
def derangement(rng, n):
    while True:
        p = rng.permutation(n)
        if not np.any(p == np.arange(n)):
            return p


def spaced_lags(rng, n, gap, lo=0):
    """n integer lags in [lo, NT - lo) whose pairwise circular distances are all >= gap."""
    while True:
        v = rng.integers(lo, NT - lo, n)
        d = np.abs(v[:, None] - v[None, :])
        d = np.minimum(d, NT - d)[np.triu_indices(n, 1)]
        if n == 1 or d.min() >= gap:
            return v


def rand_rot(rng):
    q, r = np.linalg.qr(rng.normal(size=(3, 3)))
    q = q * np.sign(np.diag(r))[None, :]
    if np.linalg.det(q) < 0:
        q[:, 0] = -q[:, 0]
    return q


class Decoys:
    """decoy engines of a network built from its target engine E (c10_lib.Eng, registered terms corr (ns, 2))."""

    def __init__(self, E, corr):
        self.E = E
        self.ns = ns = E.ns
        d = np.linalg.norm(E.XN[:, None, :] - E.STA[None, :, :], axis=2)
        tt = np.stack([d / VP, d / VS], axis=1)                                   # (N, 2, ns)
        self.x0 = (tt[:, None] + E.taus[None, :, None, None]) / DT                # no station term
        self.c = np.asarray(corr, float).T[None, None] / DT                       # (1, 1, 2, ns)
        assert np.array_equal(np.clip((self.x0 + self.c).astype(int), 0, NT - 1), E.IDX)
        self.G5 = E.G.reshape(E.N, 2, ns, 3, -1)
        self.Q5 = E.Q.reshape(E.N, 2, ns, 3, -1)

    def _with(self, IDX=None, G=None, Q=None):
        E2 = copy.copy(self.E)
        E2.cache = {}
        if IDX is not None:
            E2.IDX = IDX
        if G is not None:
            E2.G, E2.Q = G.reshape(self.E.G.shape), Q.reshape(self.E.Q.shape)
        return E2

    def d1p(self, rng):
        pi = derangement(rng, self.ns)
        return self._with(IDX=np.clip((self.x0[..., pi] + self.c).astype(int), 0, NT - 1))

    def d1o(self, rng):
        off = spaced_lags(rng, 2 * self.ns, GAP_D1O).reshape(1, 1, 2, self.ns)
        return self._with(IDX=(self.E.IDX + off) % NT)

    def d2s(self, rng):
        sg = np.ones(2 * self.ns)
        sg[rng.permutation(2 * self.ns)[:self.ns]] = -1.0
        sg = sg.reshape(1, 2, self.ns, 1, 1)
        return self._with(G=self.G5 * sg, Q=self.Q5 * sg)

    def d2r(self, rng):
        Rm = np.stack([rand_rot(rng) for _ in range(self.ns)])                    # (ns, 3, 3)
        return self._with(G=np.einsum("icd,nbidk->nbick", Rm, self.G5), Q=np.einsum("icd,nbidk->nbick", Rm, self.Q5))

    def d3_data(self, rng, MFn):
        lag = spaced_lags(rng, self.ns, GAP_D3)
        return np.stack([np.roll(MFn[i], int(lag[i]), axis=-1) for i in range(self.ns)])

    def run_all(self, rng, MFn, kmax):
        """decoy atom lists of one window (one draw of every construction)."""
        out = {}
        out["D1p"] = C.mt_loop(self.d1p(rng), MFn, kmax)
        out["D1o"] = C.mt_loop(self.d1o(rng), MFn, kmax)
        out["D2s"] = C.mt_loop(self.d2s(rng), MFn, kmax)
        out["D2r"] = C.mt_loop(self.d2r(rng), MFn, kmax)
        out["D3"] = C.mt_loop(self.E, self.d3_data(rng, MFn), kmax)
        return out


# ============================================================================================ calibration
def _thr_from_counts(sT, sD, q, add=1.0, scale=1.0):
    """smallest threshold t among the target scores with scale * (add + #D >= t) / max(1, #T >= t) <= q for every q in
    q (array).  Returns thresholds (len(q),), +inf when nothing is accepted."""
    q = np.atleast_1d(q)
    if not len(sT):
        return np.full(len(q), np.inf)
    t = np.unique(sT)                                                              # ascending
    sTs, sDs = np.sort(sT), np.sort(sD)
    nT = len(sTs) - np.searchsorted(sTs, t, side="left")
    nD = len(sDs) - np.searchsorted(sDs, t, side="left")
    est = scale * (add + nD) / np.maximum(1, nT)
    out = np.full(len(q), np.inf)
    for k, qq in enumerate(q):
        ok = np.where(est <= qq)[0]
        if len(ok):
            out[k] = t[ok[0]]
    return out


def thr_tdc(sT, sD, q):
    return _thr_from_counts(sT, sD, q)


def thr_bh(sT, sD, q):
    """BH on conformal p-values (1 + #D >= s) / (nD + 1): largest k with p_(k) <= q k / m, i.e. the smallest t with
    m (1 + #D >= t) / ((nD + 1) #T >= t) <= q."""
    return _thr_from_counts(sT, sD, q, scale=len(sT) / (len(sD) + 1.0))


def thr_comp(sT, sD, q):
    """rank-paired competition.  sT, sD aligned (same window, same rank).  Returns thresholds on W; the accepted set is
    {target-won pairs with W >= t}."""
    W = np.maximum(sT, sD)
    win = sT > sD
    return _thr_from_counts(W[win], W[~win], q)


def fdp_hat_curve(sT, sD):
    """estimated FDP (1 + #D >= t) / max(1, #T >= t) at every target score (descending)."""
    t = np.sort(np.unique(sT))[::-1]
    sTs, sDs = np.sort(sT), np.sort(sD)
    nT = len(sTs) - np.searchsorted(sTs, t, side="left")
    nD = len(sDs) - np.searchsorted(sDs, t, side="left")
    return t, (1.0 + nD) / np.maximum(1, nT)


# ============================================================================================ records
def pack(atoms):
    return [[float(a["x"][0]), float(a["x"][1]), float(a["x"][2]), float(a["tau"]), float(a["s"])] for a in atoms]


def unpack(lst):
    return [{"x": np.array(a[:3]), "tau": a[3], "s": a[4]} for a in lst]


def runmin(s):
    return np.minimum.accumulate(np.asarray(s, float)) if len(s) else np.zeros(0)
