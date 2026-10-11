"""
certified_fw.py - Certified linear-minimisation oracle for the Frank-Wolfe detector MT-FW
(task1x_mtglrt.fw_record) by branch-and-bound over an octree of (x, y, z, tau) cells.

THE STATISTIC.  For a hypothesis theta = (grid node g, origin time tau_j), station i and phase
b in {P, S}, node_vectors() takes the predicted time t_bi(theta) = |x_g - s_i| / V_b + tau_j,
the sample idx = clip(floor(t / DT), 0, NT-1), the window W = {clip(idx + o, 0, NT-1) : |o| <= SH}
(SH = 4 samples = 0.25 s), picks k* = argmax_{k in W} E_i(k) with E_i(k) = ||R_i[:, k]||^2, and
keeps the signed vector v_bi = R_i[:, k*] capped in norm at CAPV:  vt_bi = v_bi min(1, CAPV/|v_bi|).
The certificate is T(theta) = ||Q_g^T vt(theta)||^2 / 2, Q_g an orthonormal basis of the
moment-tensor subspace of node g.

THE BOUND (proved; verified numerically by check_bound()).
 (a) Projection.  Q_g has orthonormal columns, so ||Q_g^T w|| <= ||w|| for every w:
         T(theta) <= ||vt(theta)||^2 / 2 = (1/2) sum_bi ||vt_bi||^2.
 (b) Cap and argmax.  ||vt_bi||^2 = min(CAPV^2, ||v_bi||^2) = min(CAPV^2, max_{k in W(theta)} E_i(k)),
     because the selected sample is the energy argmax of the window and the cap only rescales.
     The map e -> min(CAPV^2, e) is non-decreasing, so for ANY index set W' containing W(theta)
         ||vt_bi||^2 <= min(CAPV^2, max_{k in W'} E_i(k)).
     Neither the argmax selection nor the cap can break the bound: the selected vector's energy is
     at most the largest energy of any superset window, and capping never increases a norm.
 (c) Eikonal Lipschitz.  For a cell C (a box of lattice nodes x tau-range) with centre (c, tau_c),
     spatial radius h = max_{x in C} |x - c| and half-width dtau, the triangle inequality gives
     | |x - s_i| - |c - s_i| | <= |x - c| <= h, hence |t_bi(theta) - t_bi(c, tau_c)| <= h / V_b + dtau
     =: Delta_b for every theta in C (in a heterogeneous medium V_b is replaced by V_min).  floor and
     clip are non-decreasing, so every window W(theta), theta in C, lies inside
         W_C,bi = [clip(floor((t_c - Delta_b)/DT) - SH), clip(floor((t_c + Delta_b)/DT) + SH)].
     Therefore  UB(C) = (1/2) sum_bi min(CAPV^2, max_{k in W_C,bi} E_i(k))  >=  max_{theta in C} T(theta).
     The optional 'box' mode replaces [|c-s|-h, |c-s|+h] by the exact range [dist(s, box), farthest
     corner] of distances from s_i to the cell's bounding box, which is contained in the former
     (every point of the box is within h of c), hence also valid and tighter.
 (d) Deflation.  Nothing above uses any property of R other than that T is computed from R by
     the node_vectors rule: E_i(k) is the energy of the CURRENT residual.  After each Frank-Wolfe
     step the fitted contribution is subtracted from R, E and its range-maximum table are rebuilt
     from the new R, and the bound holds for the deflated residual verbatim.
 (e) Floating point.  Cell windows use Delta_b + 1e-6 s, and every bound is inflated by a relative
     1e-9 (plus 1e-12), which dominates the rounding of ||Q^T w||^2 <= ||w||^2 and of the cap.

EXACTNESS OF THE ORACLE.  Best-first search: pop the cell with the largest bound; stop when that
bound is < the incumbent (strict).  A cell containing a hypothesis with T = T_max has bound
>= T_max >= incumbent, so it is never pruned: every node carrying a hypothesis whose value
equals the global maximum is evaluated exactly.  Exact evaluation of a node uses THE SAME CALL as
fw_record (node_vectors over the full tau array and vf @ Qs[g]), because BLAS matmul rounds rows
differently for different row counts; values are therefore bit-identical to the dense search, and
the winner is chosen with fw_record's tie rule (lowest node index, then first tau).  Nodes inside
the minimum-separation disc of accepted events are excluded with fw_record's own test.  The
deflation and the record are fw_record's code.  Hence certified_fw_record returns exactly the
records of fw_record (verified in main() on 72 windows).

    python certified_fw.py [out.txt]  # verification (72 windows); tables of certified_fw_report.txt
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import heapq
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import task1x_mtglrt as X  # noqa: E402
from task0_common import DT, NT, VP, VS  # noqa: E402
from task1o_consistent import TOL_KM  # noqa: E402

SH = X.SH
REL_INFL, ABS_INFL, T_MARGIN = 1e-9, 1e-12, 1e-6


# ----------------------------------------------------------------------------------------------
# range-maximum queries on the energy of the current residual
# ----------------------------------------------------------------------------------------------
class EnergyRMQ:
    """Sparse table: max_{lo <= k <= hi} E_i(k) in O(1) per query, vectorised over (b, i)."""

    def __init__(self, R):
        E = np.sum(R ** 2, axis=1)                        # (ns, NT), same formula as node_vectors
        self.E = E
        L = int(np.floor(np.log2(E.shape[1]))) + 1
        tab = [E]
        for lv in range(1, L):
            p = tab[-1]
            w = 1 << (lv - 1)
            q = p.copy()
            q[:, :-w] = np.maximum(p[:, :-w], p[:, w:])
            tab.append(q)
        self.tab = np.stack(tab)                          # (L, ns, NT)
        self.log2 = np.zeros(E.shape[1] + 1, int)
        for n in range(2, E.shape[1] + 1):
            self.log2[n] = self.log2[n // 2] + 1
        # exact 9-sample window maximum of node_vectors (window clip(t0 + o), |o| <= SH)
        t0 = np.arange(E.shape[1])
        self.WM = self.query(np.broadcast_to(np.arange(E.shape[0])[:, None], E.shape),
                             np.clip(t0 - SH, 0, NT - 1)[None].repeat(E.shape[0], 0),
                             np.clip(t0 + SH, 0, NT - 1)[None].repeat(E.shape[0], 0))

    def query(self, i, lo, hi):
        lv = self.log2[hi - lo + 1]
        return np.maximum(self.tab[lv, i, lo], self.tab[lv, i, hi - (1 << lv) + 1])


# ----------------------------------------------------------------------------------------------
# lattice bookkeeping
# ----------------------------------------------------------------------------------------------
class Lattice:
    def __init__(self, grid, taus):
        g = np.asarray(grid, float)
        self.garr = g
        self.ax = [np.unique(g[:, k]) for k in range(3)]
        nx, ny, nz = (len(a) for a in self.ax)
        if nx * ny * nz != len(g):
            raise ValueError("grid is not a full lattice")
        ix = [np.searchsorted(self.ax[k], g[:, k]) for k in range(3)]
        self.gidx = np.full((nx, ny, nz), -1, int)
        self.gidx[ix[0], ix[1], ix[2]] = np.arange(len(g))
        if (self.gidx < 0).any():
            raise ValueError("grid is not a full lattice")
        self.shape = (nx, ny, nz, len(taus))
        self.taus = np.asarray(taus, float)


def _cell_windows(STA, lat, cells, mode, tcorr):
    """cells (n, 8) half-open index ranges -> window limits lo, hi (n, 2, ns)."""
    ax, taus = lat.ax, lat.taus
    x0, x1, y0, y1, z0, z1, t0, t1 = cells.T
    lo_xyz = np.stack([ax[0][x0], ax[1][y0], ax[2][z0]], 1)
    hi_xyz = np.stack([ax[0][x1 - 1], ax[1][y1 - 1], ax[2][z1 - 1]], 1)
    tlo, thi = taus[t0], taus[t1 - 1]
    V = np.array([VP, VS])[None, :, None]
    if mode == "lipschitz":
        c = (lo_xyz + hi_xyz) / 2.0
        h = np.linalg.norm((hi_xyz - lo_xyz) / 2.0, axis=1)
        d = np.linalg.norm(STA[None] - c[:, None], axis=2)                 # (n, ns)
        tc = d[:, None, :] / V + ((tlo + thi) / 2.0)[:, None, None]
        dl = h[:, None, None] / V + ((thi - tlo) / 2.0)[:, None, None] + T_MARGIN
        tmin, tmax = tc - dl, tc + dl
    else:  # exact distance range to the cell's bounding box
        cl = np.clip(STA[None], lo_xyz[:, None], hi_xyz[:, None])          # (n, ns, 3)
        dmin = np.linalg.norm(STA[None] - cl, axis=2)
        far = np.where(np.abs(STA[None] - lo_xyz[:, None]) > np.abs(STA[None] - hi_xyz[:, None]),
                       lo_xyz[:, None], hi_xyz[:, None])
        dmax = np.linalg.norm(STA[None] - far, axis=2)
        tmin = dmin[:, None, :] / V + tlo[:, None, None] - T_MARGIN
        tmax = dmax[:, None, :] / V + thi[:, None, None] + T_MARGIN
    if tcorr is not None:
        tmin = tmin + tcorr.T[None]
        tmax = tmax + tcorr.T[None]
    lo = np.clip(np.floor(tmin / DT).astype(int) - SH, 0, NT - 1)
    hi = np.clip(np.floor(tmax / DT).astype(int) + SH, 0, NT - 1)
    return lo, hi


def cell_bounds(STA, lat, rmq, cells, mode="lipschitz", tcorr=None):
    """Certified upper bound of T over each cell (vectorised)."""
    cells = np.atleast_2d(cells)
    lo, hi = _cell_windows(STA, lat, cells, mode, tcorr)
    ns = len(STA)
    i = np.broadcast_to(np.arange(ns)[None, None, :], lo.shape)
    e = np.minimum(X.CAPV ** 2, rmq.query(i, lo, hi))
    return e.sum(axis=(1, 2)) / 2.0 * (1 + REL_INFL) + ABS_INFL


def _idx_node(STA, src, taus, tcorr):
    d = np.linalg.norm(STA - src, axis=1)
    tt = np.stack([d / VP, d / VS])[None] + taus[:, None, None]
    if tcorr is not None:
        tt = tt + tcorr.T[None]
    return np.clip((tt / DT).astype(int), 0, NT - 1)


def node_tau_bounds(STA, src, taus, rmq, tcorr=None):
    """Bound for every single hypothesis (node, tau): the same windows as node_vectors, so it is
    ||vt(theta)||^2 / 2 exactly (inflated)."""
    idx = _idx_node(STA, src, taus, tcorr)                                 # (nt, 2, ns)
    e = np.minimum(X.CAPV ** 2, rmq.WM[np.arange(len(STA))[None, None, :], idx])
    return e.sum(axis=(1, 2)) / 2.0 * (1 + REL_INFL) + ABS_INFL


def node_vectors_c(MFn, STA, src, taus, tcorr=None):
    """X.node_vectors, with optional station/phase time corrections tcorr (ns, 2) added to the
    predicted times.  With tcorr None it IS X.node_vectors (bit-identical)."""
    if tcorr is None:
        return X.node_vectors(MFn, STA, src, taus)
    ns = len(STA)
    st = np.arange(ns)
    idx = _idx_node(STA, src, taus, tcorr)
    ii = np.clip(idx[..., None] + X.OFFS, 0, NT - 1)
    seg = MFn[st[None, None, :, None], :, ii]
    ener = np.sum(seg ** 2, axis=-1)
    k = ener.argmax(axis=-1)
    v = np.take_along_axis(seg, k[..., None, None], axis=3)[:, :, :, 0, :]
    nv = np.linalg.norm(v, axis=-1, keepdims=True)
    v = v * np.minimum(1.0, X.CAPV / np.maximum(nv, 1e-12))
    tpk = np.take_along_axis(ii, k[..., None], axis=3)[..., 0]
    return v, tpk, seg, ener


def stack_Q(Qs):
    """(n_nodes, 6ns, 6) stack of the bases, zero-padded if some node has rank < 6."""
    r = max(Q.shape[1] for Q in Qs)
    out = np.zeros((len(Qs), Qs[0].shape[0], r))
    for g, Q in enumerate(Qs):
        out[g, :, :Q.shape[1]] = Q
    return out


class LeafTables:
    """Per-residual tables for the subspace-aware leaf enclosure: for every station and window
    centre t0 the sample chosen by node_vectors' window argmax, the capped vector there, and a
    flag when the top two energies of the window are within 1e-12 (argmax not certifiable)."""

    def __init__(self, R, rmq):
        ns = R.shape[0]
        Et = np.sum(np.ascontiguousarray(R.transpose(0, 2, 1)) ** 2, axis=-1)   # (ns, NT)
        II = np.clip(np.arange(NT)[:, None] + X.OFFS[None, :], 0, NT - 1)        # (NT, W)
        W = Et[:, II]                                                             # (ns, NT, W)
        a = W.argmax(axis=-1)
        self.K = II[np.arange(NT)[None, :], a]                                    # (ns, NT)
        # near-tie between the chosen sample and a DIFFERENT sample of the window (clipping at the
        # record ends duplicates the edge sample, which is not a tie)
        top = W.max(axis=-1)
        other = np.where(II[None] == self.K[..., None], -np.inf, W).max(axis=-1)
        self.amb = top - other <= 1e-12 * np.maximum(top, 1e-300)
        V = np.ascontiguousarray(R.transpose(0, 2, 1))                            # (ns, NT, 3)
        nv = np.linalg.norm(V, axis=-1, keepdims=True)
        self.VC = V * np.minimum(1.0, X.CAPV / np.maximum(nv, 1e-12))
        self.rmq = rmq


def leaf_enclosure(STA, garr, G, taus, Qst, lt, tcorr=None):
    """Rigorous upper bound, for every hypothesis (g in G, tau in taus), of the value fw_record
    computes in floating point.  The statistic is re-evaluated in batch (einsum); the rounding of
    both evaluations of ||Q^T v||^2 / 2 is < 60 u ||v||^2 (u = 1.1e-16, 6ns <= 60 terms), so
    adding 1e-11 (||v||^2 + 1) is a certified enclosure.  Hypotheses whose sample index or window
    argmax cannot be certified (t/DT within 1e-9 of an integer, or a top-2 energy tie) fall back
    to the energy bound over the window widened by one sample."""
    ns = len(STA)
    d = np.linalg.norm(STA[None] - garr[G][:, None], axis=2)                     # (nG, ns)
    tt = np.stack([d / VP, d / VS], axis=1)[:, None] + taus[None, :, None, None]  # (nG, nt, 2, ns)
    if tcorr is not None:
        tt = tt + tcorr.T[None, None]
    q = tt / DT
    idx = np.clip(q.astype(int), 0, NT - 1)
    fr = q - np.floor(q)
    st = np.arange(ns)[None, None, None, :]
    amb = (fr < 1e-9) | (fr > 1 - 1e-9) | lt.amb[st, idx]
    k = lt.K[st, idx]
    vf = lt.VC[st, k].reshape(len(G), len(taus), 6 * ns)                          # (nG, nt, 6ns)
    y = np.einsum("gtr,grk->gtk", vf, Qst[G])
    en2 = np.einsum("gtr,gtr->gt", vf, vf)
    enc = (y * y).sum(-1) / 2.0 + 1e-11 * (en2 + 1.0)
    if amb.any():
        bad = amb.any(axis=(2, 3))
        lo = np.clip(idx - SH - 1, 0, NT - 1)
        hi = np.clip(idx + SH + 1, 0, NT - 1)
        stb = np.broadcast_to(st, idx.shape)
        eb = np.minimum(X.CAPV ** 2, lt.rmq.query(stb, lo, hi)).sum(axis=(2, 3)) / 2.0 * (1 + REL_INFL) + ABS_INFL
        enc = np.where(bad, np.maximum(enc, eb), enc)
    return enc                                                                    # (nG, nt)


def _exclusion_update(mask, garr, a, min_sep):
    """fw_record's test any(norm(src[:2]-a[:2]) <= min_sep), vectorised with an exact recheck of
    the nodes whose distance is within 1e-9 km of the radius."""
    dd = np.linalg.norm(garr[:, :2] - a[:2], axis=1)
    new = dd <= min_sep
    amb = np.where(np.abs(dd - min_sep) < 1e-9)[0]
    for gi in amb:
        new[gi] = bool(np.linalg.norm(garr[gi][:2] - a[:2]) <= min_sep)
    return mask | new


def _split_dim(lat, cell):
    ax = lat.ax
    ext = []
    for k in range(3):
        a, b = cell[2 * k], cell[2 * k + 1]
        ext.append((ax[k][b - 1] - ax[k][a]) / 2.0 / VS if b - a > 1 else -1.0)
    a, b = cell[6], cell[7]
    ext.append((lat.taus[b - 1] - lat.taus[a]) / 2.0 if b - a > 1 else -1.0)
    return int(np.argmax(ext))


# ----------------------------------------------------------------------------------------------
# the oracle and the certified Frank-Wolfe record
# ----------------------------------------------------------------------------------------------
def certified_oracle(STA, R, lat, taus, Qs, excluded, mode="lipschitz", tcorr=None, stats=None,
                     Qst=None, leaf_nodes=16):
    """Global maximiser of T over non-excluded hypotheses, with fw_record's tie rule.
    Returns (s, gi, vf_row, tpk_row) or (-inf, None, None, None).
    mode 'lipschitz' / 'box': energy bound of the module docstring down to single nodes.
    mode 'enclosure': the 'box' energy bound on coarse cells; a cell with <= leaf_nodes nodes is
    resolved by leaf_enclosure() (subspace-aware, certified), and only the nodes whose enclosure
    reaches the incumbent are evaluated with fw_record's exact call."""
    rmq = EnergyRMQ(R)
    enclosure = mode == "enclosure"
    bmode = "box" if enclosure else mode
    if enclosure:
        lt = LeafTables(R, rmq)
        if Qst is None:
            Qst = stack_Q(Qs)
    n_hyp_enc = 0
    nx, ny, nz, nt = lat.shape
    free = (~excluded[lat.gidx]).astype(int)                                # (nx, ny, nz)
    P = np.zeros((nx + 1, ny + 1, nz + 1), int)
    P[1:, 1:, 1:] = free.cumsum(0).cumsum(1).cumsum(2)

    def n_free(c):
        x0, x1, y0, y1, z0, z1 = c[:6]
        return (P[x1, y1, z1] - P[x0, y1, z1] - P[x1, y0, z1] - P[x1, y1, z0]
                + P[x0, y0, z1] + P[x0, y1, z0] + P[x1, y0, z0] - P[x0, y0, z0])

    root = np.array([0, nx, 0, ny, 0, nz, 0, nt])
    if n_free(root) == 0:
        return (-np.inf, None, None, None)
    ub0 = float(cell_bounds(STA, lat, rmq, root[None], bmode, tcorr)[0])
    heap = [(-ub0, 0, tuple(root))]                   # entries: (-key, counter, cell) or (-key, counter, node)
    cnt, n_bound, n_nodebound = 1, 1, 0
    inc = -np.inf
    cache = {}

    def exact(gi):
        v, tpk, _s, _e = node_vectors_c(R, STA, lat.garr[gi], taus, tcorr)       # fw_record's call
        vf = v.reshape(len(taus), -1)
        sc = np.sum((vf @ Qs[gi]) ** 2, axis=1) / 2.0
        cache[gi] = (sc, vf, tpk)
        return float(sc.max())

    while heap:
        nub, _c, item = heapq.heappop(heap)
        if -nub < inc:
            break
        if isinstance(item, int):                     # a node certified by its enclosure
            if item not in cache:
                inc = max(inc, exact(item))
            continue
        cell = item
        x0, x1, y0, y1, z0, z1, t0, t1 = cell
        nn = (x1 - x0) * (y1 - y0) * (z1 - z0)
        if enclosure and nn <= leaf_nodes:
            G = [int(g) for g in lat.gidx[x0:x1, y0:y1, z0:z1].ravel() if not excluded[g] and g not in cache]
            if not G:
                continue
            enc = leaf_enclosure(STA, lat.garr, G, taus[t0:t1], Qst, lt, tcorr)
            n_hyp_enc += enc.size
            keys = enc.max(axis=1)
            for g, kk in zip(G, keys):
                if kk >= inc:
                    heapq.heappush(heap, (-float(kk), cnt, g))
                    cnt += 1
            continue
        if nn == 1:
            gi = int(lat.gidx[x0, y0, z0])
            if excluded[gi] or gi in cache:
                continue
            ubt = node_tau_bounds(STA, lat.garr[gi], taus[t0:t1], rmq, tcorr)
            n_nodebound += 1
            if ubt.max() < inc:
                continue
            inc = max(inc, exact(gi))
            continue
        k = _split_dim(lat, cell)
        a, b = cell[2 * k], cell[2 * k + 1]
        m = (a + b) // 2
        kids = []
        for part in ((a, m), (m, b)):
            c2 = list(cell)
            c2[2 * k], c2[2 * k + 1] = part
            if n_free(c2) > 0:
                kids.append(c2)
        if not kids:
            continue
        ubs = cell_bounds(STA, lat, rmq, np.array(kids), bmode, tcorr)
        n_bound += len(kids)
        for c2, u in zip(kids, ubs):
            if u >= inc:
                heapq.heappush(heap, (-float(u), cnt, tuple(c2)))
                cnt += 1
    # fw_record's tie rule: ascending node index, first tau, strict improvement
    best = (-np.inf, None, None, None)
    for gi in sorted(cache):
        sc, vf, tpk = cache[gi]
        j = int(sc.argmax())
        if sc[j] > best[0]:
            best = (float(sc[j]), gi, vf[j], tpk[j])
    if stats is not None:
        stats["node_evals"].append(len(cache))
        stats["cell_bounds"].append(n_bound)
        stats["node_bounds"].append(n_nodebound)
        stats["dense_nodes"].append(int((~excluded).sum()))
        stats["value"].append(best[0])
        stats.setdefault("hyp_enclosed", []).append(n_hyp_enc)
    return best


def certified_fw_record(STA, MFn, grid, taus, Gs, Qs, truth, kmax=5, min_sep=5.0,
                        mode="enclosure", tcorr=None, stats=None):
    """Exactly the records of task1x_mtglrt.fw_record, with the oracle by certified
    branch-and-bound.  tcorr (ns, 2): optional station/phase time corrections (s), added to the
    predicted P and S times (None = fw_record)."""
    ns = len(STA)
    R = MFn.copy()
    lat = Lattice(grid, taus)
    garr = lat.garr
    taus = np.asarray(taus)
    excluded = np.zeros(len(garr), bool)
    Qst = stack_Q(Qs) if mode == "enclosure" else None
    acc, rec, run = [], [], np.inf
    for _ in range(kmax):
        s, gi, vf, tpk = certified_oracle(STA, R, lat, taus, Qs, excluded, mode, tcorr, stats, Qst)
        if gi is None:
            break
        src = garr[gi]
        run = min(run, s)
        hit = -1
        if truth:
            dd = [float(np.linalg.norm(src[:2] - np.asarray(t_)[:2])) for t_ in truth]
            j = int(np.argmin(dd))
            if dd[j] <= TOL_KM:
                hit = j
        rec.append((run, hit))
        acc.append(src)
        excluded = _exclusion_update(excluded, garr, src, min_sep)
        G = Gs[gi]
        m, *_ = np.linalg.lstsq(G, vf, rcond=None)
        pred = (G @ m).reshape(2, ns, 3)
        for b in range(2):
            for i in range(ns):
                t0 = int(tpk[b, i])
                lo, hi = max(0, t0 - X.AC_HALF), min(NT, t0 + X.AC_HALF + 1)
                R[i, :, lo:hi] -= np.outer(pred[b, i], X.AC[X.AC_HALF - (t0 - lo):X.AC_HALF + (hi - t0)])
    return rec


def dense_fw_record_c(STA, MFn, grid, taus, Gs, Qs, truth, kmax=5, min_sep=5.0, tcorr=None):
    """fw_record with optional time corrections (reference implementation for tcorr != None:
    the loop of task1x_mtglrt.fw_record with node_vectors_c)."""
    if tcorr is None:
        return X.fw_record(STA, MFn, grid, taus, Gs, Qs, truth, kmax, min_sep)
    ns = len(STA)
    R = MFn.copy()
    garr = np.asarray(grid)
    acc, rec, run = [], [], np.inf
    for _ in range(kmax):
        best = (-np.inf, None, None, None)
        for gi, src in enumerate(garr):
            if any(np.linalg.norm(src[:2] - a[:2]) <= min_sep for a in acc):
                continue
            v, tpk, _seg, _e = node_vectors_c(R, STA, src, taus, tcorr)
            vf = v.reshape(len(taus), -1)
            sc = np.sum((vf @ Qs[gi]) ** 2, axis=1) / 2.0
            j = int(sc.argmax())
            if sc[j] > best[0]:
                best = (float(sc[j]), gi, vf[j], tpk[j])
        s, gi, vf, tpk = best
        if gi is None:
            break
        src = garr[gi]
        run = min(run, s)
        hit = -1
        if truth:
            dd = [float(np.linalg.norm(src[:2] - np.asarray(t_)[:2])) for t_ in truth]
            j = int(np.argmin(dd))
            if dd[j] <= TOL_KM:
                hit = j
        rec.append((run, hit))
        acc.append(src)
        G = Gs[gi]
        m, *_ = np.linalg.lstsq(G, vf, rcond=None)
        pred = (G @ m).reshape(2, ns, 3)
        for b in range(2):
            for i in range(ns):
                t0 = int(tpk[b, i])
                lo, hi = max(0, t0 - X.AC_HALF), min(NT, t0 + X.AC_HALF + 1)
                R[i, :, lo:hi] -= np.outer(pred[b, i], X.AC[X.AC_HALF - (t0 - lo):X.AC_HALF + (hi - t0)])
    return rec


# ----------------------------------------------------------------------------------------------
# direct test of the bound on random cells (current and deflated residuals)
# ----------------------------------------------------------------------------------------------
def dense_field(STA, R, lat, taus, Qs, tcorr=None):
    F = np.empty((len(lat.garr), len(taus)))
    for gi, src in enumerate(lat.garr):
        v, _t, _s, _e = node_vectors_c(R, STA, src, taus, tcorr)
        F[gi] = np.sum((v.reshape(len(taus), -1) @ Qs[gi]) ** 2, axis=1) / 2.0
    return F


def check_bound(STA, R, lat, taus, Qs, rng, n_cells=300, tcorr=None):
    """Returns (violations lipschitz, violations box, violations node-tau, median UB/max ratio)."""
    rmq = EnergyRMQ(R)
    F = dense_field(STA, R, lat, taus, Qs, tcorr)
    F4 = F[lat.gidx]                                                        # (nx, ny, nz, nt)
    nx, ny, nz, nt = lat.shape
    cells = []
    for _ in range(n_cells):
        c = []
        for n in (nx, ny, nz, nt):
            a = int(rng.integers(0, n))
            b = min(n, a + int(rng.choice([1, 2, 4, 8, n])))
            c += [a, b]
        cells.append(c)
    cells = np.array(cells)
    truemax = np.array([F4[c[0]:c[1], c[2]:c[3], c[4]:c[5], c[6]:c[7]].max() for c in cells])
    ub_l = cell_bounds(STA, lat, rmq, cells, "lipschitz", tcorr)
    ub_b = cell_bounds(STA, lat, rmq, cells, "box", tcorr)
    vl = int((ub_l < truemax).sum())
    vb = int((ub_b < truemax).sum())
    vn = 0
    for gi, src in enumerate(lat.garr):
        vn += int((node_tau_bounds(STA, src, taus, rmq, tcorr) < F[gi]).sum())
    lt = LeafTables(R, rmq)
    enc = leaf_enclosure(STA, lat.garr, list(range(len(lat.garr))), np.asarray(taus), stack_Q(Qs), lt, tcorr)
    ve = int((enc < F).sum())
    gap = float(np.max(enc - F))                     # enclosure slack (should be ~1e-9)
    return vl, vb, vn, float(np.median(ub_l / truemax)), float(np.median(ub_b / truemax)), F, ve, gap


# ----------------------------------------------------------------------------------------------
# verification protocol (task1y_clutter windows, seed 999 streams)
# ----------------------------------------------------------------------------------------------
def make_windows(scen, clut, n, seed=999, evbase=25000, moff=233):
    import task1y_clutter as Y
    from task0_common import make_network
    from task1o_consistent import draw_sources
    from task1r_fair import SCEN
    kind, nst, nclut = SCEN[scen]
    rng = np.random.default_rng(seed + 41 + evbase)
    rng_m = np.random.default_rng(seed + moff)
    rng_h = np.random.default_rng(seed + 61 + evbase)
    out = []
    for k in range(n):
        S = make_network(kind, nst, seed=seed + evbase + k)
        srcs, taus = draw_sources(rng)
        D, _p, tr, _a = Y.generate_y(S, srcs, taus, nclut, rng, rng_m, rng_h, Y.CLUT_AMP[clut])
        out.append((S, D, tr))
    return out


MODES = ("lipschitz", "box", "enclosure")


def main():
    from task0_common import DEFAULT_TAUS, default_grid
    from task1o_consistent import matched
    grid = default_grid()
    garr = np.asarray(grid)
    lat = Lattice(grid, DEFAULT_TAUS)
    NWIN = int(os.environ.get("CFW_NWIN", "12"))
    log = []
    P = lambda s="": (print(s, flush=True), log.append(s))  # noqa: E731
    rows, per_iter = [], {m: {} for m in MODES}
    brng = np.random.default_rng(12345)
    bc = {"n": 0, "vl": 0, "vb": 0, "vn": 0, "ve": 0, "rl": [], "rb": [], "gap": 0.0}
    mism = []
    for scen in ("ladder_4", "ladder_6", "unilateral_c40"):
        for clut in ("thin", "strong"):
            wins = make_windows(scen, clut, NWIN)
            a = {"t_dense": 0.0, "dense": 0, "hyp_dense": 0}
            for m in MODES:
                a.update({f"same_{m}": 0, f"t_{m}": 0.0, f"ev_{m}": 0, f"cb_{m}": 0, f"enc_{m}": 0})
            for w, (S, D, tr) in enumerate(wins):
                Gs = [X.radiation_operator(g, S) for g in garr]
                Qs = [X.basis(G) for G in Gs]
                MF, sd = matched(D)
                MFn = MF / sd
                t = time.perf_counter()
                r_d = X.fw_record(S, MFn, grid, DEFAULT_TAUS, Gs, Qs, tr)
                a["t_dense"] += time.perf_counter() - t
                for m in MODES:
                    st = {k: [] for k in ("node_evals", "cell_bounds", "node_bounds", "dense_nodes", "value")}
                    t = time.perf_counter()
                    r_c = certified_fw_record(S, MFn, grid, DEFAULT_TAUS, Gs, Qs, tr, mode=m, stats=st)
                    a[f"t_{m}"] += time.perf_counter() - t
                    a[f"same_{m}"] += int(r_c == r_d)
                    if r_c != r_d:
                        mism.append((scen, clut, w, m, r_d, r_c))
                    a[f"ev_{m}"] += sum(st["node_evals"])
                    a[f"cb_{m}"] += sum(st["cell_bounds"]) + sum(st["node_bounds"])
                    a[f"enc_{m}"] += sum(st["hyp_enclosed"])
                    if m == MODES[0]:
                        a["dense"] += sum(st["dense_nodes"])
                        a["hyp_dense"] += sum(st["dense_nodes"]) * len(DEFAULT_TAUS)
                    for it, (ne, dn, val) in enumerate(zip(st["node_evals"], st["dense_nodes"], st["value"])):
                        per_iter[m].setdefault((len(S), it), []).append((ne / dn, val))
                if w % 3 == 0:          # direct test of every bound on the raw and a deflated residual
                    R = MFn.copy()
                    for _it in range(2):
                        vl, vb, vn, rl, rb, F, ve, gap = check_bound(S, R, lat, DEFAULT_TAUS, Qs, brng)
                        bc["n"] += 1
                        bc["vl"] += vl
                        bc["vb"] += vb
                        bc["vn"] += vn
                        bc["ve"] += ve
                        bc["rl"].append(rl)
                        bc["rb"].append(rb)
                        bc["gap"] = max(bc["gap"], gap)
                        gi, j = np.unravel_index(int(np.argmax(F)), F.shape)
                        v, tpk, _s, _e = X.node_vectors(R, S, garr[gi], DEFAULT_TAUS)
                        vf = v.reshape(len(DEFAULT_TAUS), -1)[j]
                        mm, *_ = np.linalg.lstsq(Gs[gi], vf, rcond=None)
                        pred = (Gs[gi] @ mm).reshape(2, len(S), 3)
                        for b in range(2):
                            for i in range(len(S)):
                                t0 = int(tpk[j][b, i])
                                lo, hi = max(0, t0 - X.AC_HALF), min(NT, t0 + X.AC_HALF + 1)
                                R[i, :, lo:hi] -= np.outer(pred[b, i], X.AC[X.AC_HALF - (t0 - lo):X.AC_HALF + (hi - t0)])
            rows.append((scen, clut, NWIN, a))
            P(f"{scen:15s} {clut:6s}: identical " + " ".join(f"{m} {a['same_' + m]}/{NWIN}" for m in MODES)
              + f" | dense {a['t_dense']:.1f}s " + " ".join(f"{m} {a['t_' + m]:.1f}s" for m in MODES))
    tot = {k: sum(r[3][k] for r in rows) for k in rows[0][3]}
    nwin = sum(r[2] for r in rows)
    L = []
    L.append("CERTIFIED FRANK-WOLFE ORACLE (Component B) - certified_fw.py")
    L.append("=" * 110)
    L.append("")
    L.append("1. THE BOUND (full proof in the module docstring of certified_fw.py)")
    L.append("   T(theta) = ||Q_g^T vt(theta)||^2/2, vt = per-(station, phase) vector at the energy argmax of the")
    L.append("   +-0.25 s window around t_bi(theta), capped at |v| <= CAPV = 6.")
    L.append("   (a) ||Q^T w|| <= ||w|| (orthonormal columns)  =>  T <= (1/2) sum_bi ||vt_bi||^2.")
    L.append("   (b) ||vt_bi||^2 = min(CAPV^2, max_{k in W(theta)} E_i(k)): the argmax sample has the largest energy")
    L.append("       of its window and the cap only rescales; e -> min(36, e) is non-decreasing, so for any window")
    L.append("       W' containing W(theta): ||vt_bi||^2 <= min(36, max_{k in W'} E_i(k)).  Cap and argmax cannot break it.")
    L.append("   (c) |t_bi(theta) - t_bi(centre)| <= h/V_b + dtau (triangle inequality = eikonal Lipschitz; V_min in a")
    L.append("       heterogeneous medium); floor and clip are monotone, so all windows of the cell lie in")
    L.append("       [floor((t_c - D)/DT) - 4, floor((t_c + D)/DT) + 4]; range maxima by a sparse table (O(1) per query).")
    L.append("       'box' variant: the exact distance range from station to the cell's box (subset of centre +- h).")
    L.append("   (d) The bound only uses the energies of the CURRENT residual R; E and the sparse table are rebuilt")
    L.append("       after each deflation, so it holds verbatim for the deflated residual.")
    L.append("   (e) +1e-6 s on window limits and x(1 + 1e-9) + 1e-12 on bounds absorb floating-point rounding.")
    L.append("   Subspace-aware leaf enclosure ('enclosure' mode): for cells of <= 16 nodes the statistic itself is")
    L.append("   re-evaluated in batch (einsum) from the same argmax/cap tables; the two floating-point evaluations of")
    L.append("   ||Q^T v||^2/2 differ by < 60u||v||^2 (u = 1.1e-16), so value + 1e-11(||v||^2 + 1) is a certified upper")
    L.append("   bound of fw_record's computed value.  Hypotheses whose sample index (t/DT within 1e-9 of an integer)")
    L.append("   or window argmax (top-2 energy tie within 1e-12) is not certifiable fall back to the energy bound.")
    L.append("")
    L.append("2. EXACTNESS.  Best-first search stops when the largest remaining bound is < incumbent (strict), so every")
    L.append("   node carrying a hypothesis equal to the global maximum is evaluated exactly.  Exact evaluation of a node")
    L.append("   is fw_record's own call (node_vectors over all 31 taus, then vf @ Qs[g]): BLAS matmul rounds a row")
    L.append("   differently for different row counts (measured: mismatch in 5-90% of random trials for subsets of 1-16")
    L.append("   rows), so evaluating tau subsets would NOT be bit-identical.  The winner is chosen with fw_record's tie")
    L.append("   rule (ascending node, first tau, strict >); exclusion, deflation and records are fw_record's code.")
    L.append("")
    L.append("3. VERIFICATION: task1y_clutter.generate_y (pert), seed 999 streams (rng 999+41+25000, rng_m 999+233,")
    L.append(f"   rng_h 999+61+25000, networks make_network(kind, n, 999+25000+k)), {NWIN} windows per cell, 0-2 sources.")
    L.append("")
    L.append(f"{'cell':15s} {'clutter':7s} {'win':>4s} | {'identical records (lip / box / encl)':>37s}")
    for scen, clut, n, a in rows:
        L.append(f"{scen:15s} {clut:7s} {n:4d} | {a['same_lipschitz']:>13d} / {a['same_box']:>3d} / {a['same_enclosure']:>3d}  of {n}")
    L.append(f"{'TOTAL':15s} {'':7s} {nwin:4d} | {tot['same_lipschitz']:>13d} / {tot['same_box']:>3d} / {tot['same_enclosure']:>3d}  of {nwin}")
    if mism:
        L.append("MISMATCHES:")
        for x in mism:
            L.append("  " + str(x))
    L.append("")
    L.append("4. COST.  Exact evaluations = nodes evaluated with fw_record's call (x31 taus each).  'Batched hyp.' =")
    L.append("   hypotheses whose statistic was recomputed in batch for the enclosure (not bit-exact, certified by the")
    L.append("   error bound).  Bounds = cell + node-level bound evaluations.  Wall time: one core, same windows.")
    L.append(f"{'cell':15s} {'clutter':7s} | {'dense':>6s} {'lip':>6s} {'box':>6s} {'encl':>5s} | {'batched hyp.':>12s} {'(of dense)':>10s} | "
             f"{'t_dense':>7s} {'t_lip':>6s} {'t_box':>6s} {'t_encl':>6s} | {'speed-up lip/box/encl':>21s}")
    for scen, clut, n, a in rows + [("TOTAL", "", nwin, tot)]:
        L.append(f"{scen:15s} {clut:7s} | {a['dense']:6d} {a['ev_lipschitz']:6d} {a['ev_box']:6d} {a['ev_enclosure']:5d} | "
                 f"{a['enc_enclosure']:12d} {a['enc_enclosure'] / a['hyp_dense']:10.2f} | {a['t_dense']:6.1f}s {a['t_lipschitz']:5.1f}s "
                 f"{a['t_box']:5.1f}s {a['t_enclosure']:5.1f}s | {a['t_dense'] / a['t_lipschitz']:5.2f}x {a['t_dense'] / a['t_box']:5.2f}x "
                 f"{a['t_dense'] / a['t_enclosure']:5.2f}x")
    L.append(f"   Reduction of exact node evaluations: Lipschitz {tot['dense'] / tot['ev_lipschitz']:.2f}x, box "
             f"{tot['dense'] / tot['ev_box']:.2f}x, enclosure {tot['dense'] / tot['ev_enclosure']:.0f}x.")
    L.append("")
    L.append("5. WHERE THE ENERGY BOUND PRUNES.  Fraction of non-excluded nodes evaluated exactly (box mode), per")
    L.append("   Frank-Wolfe iteration, with the median oracle value (global maximum of T at that iteration):")
    L.append(f"{'stations':>8s} " + " ".join(f"{'iteration ' + str(k + 1):>18s}" for k in range(5)))
    for nsx in sorted({k[0] for k in per_iter["box"]}):
        cells = []
        for it in range(5):
            v = per_iter["box"].get((nsx, it), [])
            cells.append(f"{np.mean([x for x, _ in v]):5.2f} (T={np.median([y for _, y in v]):5.1f})" if v else "")
        L.append(f"{nsx:8d} " + " ".join(f"{c:>18s}" for c in cells))
    L.append("")
    L.append("6. DIRECT TEST OF THE BOUNDS on random cells (sizes 1..full per axis), raw residual and after one")
    L.append(f"   deflation: {bc['n']} residual fields x 300 cells; violations: Lipschitz {bc['vl']}, box {bc['vb']};")
    L.append(f"   single-hypothesis energy bound {bc['vn']} and enclosure {bc['ve']} violations over "
             f"{bc['n'] * len(garr) * len(DEFAULT_TAUS)} hypotheses; max enclosure slack {bc['gap']:.1e}.")
    L.append(f"   Median (cell bound / true cell maximum): Lipschitz {np.median(bc['rl']):.2f}, box {np.median(bc['rb']):.2f}.")
    L.append("")
    return L


if __name__ == "__main__":
    out = main()
    if len(sys.argv) > 1:                       # optional: write the tables to the given path
        with open(sys.argv[1], "w") as fh:
            fh.write("\n".join(out) + "\n")
