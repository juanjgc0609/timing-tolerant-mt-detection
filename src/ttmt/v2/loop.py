"""The v2 greedy loop: onset reading (M1) or joint reading (M4), ACAP ranking (M5), deflation at the read time.

Per step:
  1. read every block of every hypothesis at the onset (``OnsetReader`` on the remaining onset evidence), and the
     amplitude vector at onset + lag with sub-sample interpolation;
  2. score every (node, origin) with the ranking statistic (ACAP, c = 3, by default: M5, Appendix AB C.4);
  3. optional M4: re-read the n_joint best hypotheses jointly and take the best joint objective;
  4. deflate the fitted atom at the sub-sample time (``deflate_at``): onset + lag, refined on the waveform in a
     narrow window (``refine_times``, +-0.25 s; E62's "fractional" fits over the whole +-3 s response support, where
     a neighbouring arrival can pull the delay). Without the refinement the picker's error sets the residual
     (M6C-SUBSAMPLE: 0.02 s of onset error leaves a median 22 % peak residual);
  5. optionally suppress, in q, the basin of every onset maximum the atom used (``peak_basin``); off by default
     (an onset shared by two overlapping sources must stay available to the second).
The 5 km exclusion of the registered loop is kept (its replacement, the resolvability test of M6, is step 7 of the
build order).
"""
import numpy as np

from ..detector import ACAP, subsample_delay
from ..kernel import LEGACY as LEGACY_KERNEL
from .m1 import OnsetReader, peak_basin, read_vectors
from .m4 import joint_read

RANK_DEFAULT = ACAP(3.0)
REFINE_HALF = 4            # samples (0.25 s): about the main lobe of the 2 Hz matched-filter response


def refine_times(R, t, kernel=LEGACY_KERNEL, bound=1.0, half=REFINE_HALF):
    """waveform refinement of read times t (2, ns): the onset reading decides WHICH arrival (no lobe jump); the
    matched filter decides WHERE, to sub-sample: the delay of the analytic response that best fits the residual in
    +-half samples of the nearest sample, within +-bound samples (``subsample_delay``). A narrow window keeps a
    neighbouring arrival out of the fit. Returns (times (2, ns), amplitudes (2, ns, 3)): the least-squares 3-vector
    of the delayed response over the same window, in the units of a vector read at the peak (NaN times: zeros)."""
    ns, nt, H = R.shape[0], R.shape[-1], kernel.ac_half
    out = np.array(t, float, copy=True)
    amp = np.zeros((2, ns, 3))
    for b in range(2):
        for i in range(ns):
            if np.isfinite(t[b, i]):
                k = int(np.clip(np.floor(t[b, i] + 0.5), 0, nt - 1))
                d = subsample_delay(R[i], k, kernel, bound=bound, half=half)
                out[b, i] = k + d
                lo, hi = max(0, k - half), min(nt, k + half + 1)
                r = kernel.response(d)[H - (k - lo):H + (hi - k)]
                amp[b, i] = R[i, :, lo:hi] @ r / float(r @ r)
    return out, amp


def deflate_at(R, pred, t, kernel=LEGACY_KERNEL):
    """subtract pred (2, ns, 3) times the deflation response centred at the float sample times t (2, ns); NaN times
    or zero vectors are skipped. In place. The response is the analytic one (Kernel.response), so an arrival read
    at its true time and amplitude is removed exactly."""
    ns, nt = R.shape[0], R.shape[-1]
    H = kernel.ac_half
    cache = {}
    for b in range(2):
        for i in range(ns):
            if not np.isfinite(t[b, i]) or not np.any(pred[b, i]):
                continue
            k = int(np.floor(t[b, i] + 0.5))
            d = round(float(t[b, i] - k), 6)
            if d not in cache:
                cache[d] = kernel.response(d)
            r = cache[d]
            lo, hi = max(0, k - H), min(nt, k + H + 1)
            if hi <= lo:
                continue
            R[i, :, lo:hi] -= np.outer(pred[b, i], r[H - (k - lo):H + (hi - k)])
    return R


def v2_loop(E, MFn, q, kmax, lag, statistic=RANK_DEFAULT, reader=None, joint=None, n_joint=20, sig=None,
            min_sep=5.0, kernel=LEGACY_KERNEL, refine=True, suppress=False, return_residual=False):
    """E: OnsetEngine; MFn (ns, 3, NT) normalised matched filter (whitened or not; E.G and E.Q must match);
    q (ns, 2, NT) onset evidence; lag (ns, 2) onset-to-peak lag (samples); statistic: ranking BlockStatistic;
    reader: OnsetReader (default half width 4 samples); joint: None or a dict of joint_read keyword arguments
    (enables M4 on the n_joint best hypotheses; sig (2, ns) samples is then required, the tolerance of M2, or a
    constant); refine: refine the deflation times and amplitudes on the waveform (``refine_times``) and refit the
    tensor on them; False: deflate the ranking fit at onset + lag, whose error is the picker's. Atoms: dicts x, tau, s (ranking score), obj (joint objective or None), onsets
    (2, ns) and deflation times t_defl (2, ns), in samples. suppress: remove from q the basin of every onset the
    atom used (off by default: two sources whose onsets merge into one maximum at a station would lose that block
    for the second source; the 5 km exclusion already prevents a repeat, and a deflated onset reads a near-zero
    vector)."""
    reader = reader or OnsetReader()
    ns = E.ns
    R = np.array(MFn, float, copy=True)
    qr = np.array(q, float, copy=True)
    lagT = np.asarray(lag, float).T
    excluded = np.zeros(E.N, bool)
    atoms = []
    for _ in range(kmax):
        if excluded.all():
            break
        s_on = reader.read(qr, E.tpred)                                   # (N, T, 2, ns)
        t_amp = s_on + lagT[None, None]
        vf = read_vectors(R, t_amp).reshape(E.N, E.T, 6 * ns)
        sc, aux = statistic.score(vf, E.Q)
        sc[excluded] = -np.inf
        if not np.isfinite(sc).any():
            break
        obj = None
        if joint is not None:
            flat = np.argsort(sc, axis=None)[::-1][:n_joint]
            sg = np.broadcast_to(np.asarray(sig, float), (2, ns))
            best = (-np.inf, None)
            for f in flat:
                gi, tj = divmod(int(f), E.T)
                if not np.isfinite(sc[gi, tj]):
                    continue
                s, v, o, _d = joint_read(qr, R, E.tpred[gi, tj], sg, lag, E.Q[gi], statistic, **joint)
                if o > best[0]:
                    best = (o, (gi, tj, s, v))
            obj, (gi, tj, s_sel, v_sel) = best[0], best[1]
            scj, auxj = statistic.score(v_sel.reshape(1, 1, -1), E.Q[gi][None])
            s_val = float(scj[0, 0])
            m, *_ = np.linalg.lstsq(E.G[gi], auxj[0, 0], rcond=None)
            pred = (E.G[gi] @ m).reshape(2, ns, 3)
        else:
            gi = int(np.argmax(sc.max(axis=1)))
            tj = int(np.argmax(sc[gi]))
            s_val = float(sc[gi, tj])
            s_sel = s_on[gi, tj]
            pred = statistic.deflation(E, gi, tj, vf, aux)
        t_sel = s_sel + lagT
        if refine:
            t_sel, amp = refine_times(R, t_sel, kernel)
            # amplitude re-read at the refined times, then the moment tensor refitted on it: the deflated atom stays
            # a moment-tensor atom (it cannot absorb another arrival block by block, unlike the "refit" convention)
            m, *_ = np.linalg.lstsq(E.G[gi], amp.reshape(-1), rcond=None)
            pred = (E.G[gi] @ m).reshape(2, ns, 3)
        src = np.asarray(E.XN[gi], float)
        atoms.append({"s": s_val, "obj": obj, "x": src, "tau": float(E.taus[tj]), "onsets": s_sel.copy(),
                      "t_defl": t_sel.copy()})
        excluded = excluded | (np.linalg.norm(E.XN[:, :2] - src[None, :2], axis=1) <= min_sep)
        deflate_at(R, pred, t_sel, kernel)
        for b in range(2):
            for i in range(ns):
                if suppress and np.isfinite(s_sel[b, i]):
                    k = int(np.clip(np.round(s_sel[b, i]), 0, qr.shape[-1] - 1))
                    # the read maximum is at k or a neighbour (sub-sample refinement); take the local max nearby
                    w = np.arange(max(0, k - 1), min(qr.shape[-1], k + 2))
                    k = int(w[np.argmax(qr[i, b, w])])
                    lo, hi = peak_basin(qr[i, b], k)
                    qr[i, b, lo:hi] = 0.0
    return (atoms, R) if return_residual else atoms
