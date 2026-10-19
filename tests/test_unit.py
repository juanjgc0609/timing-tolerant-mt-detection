"""Unit tests of the core functions (ROADMAP phase 1): sample index, subspace projection, radial maps and block
statistics, readers, deflation, hit rules, exact FROC, partition boundary and shifts, kernel tie, station-term typing."""
from pathlib import Path

import numpy as np
import pytest

from ttmt import scoring as S
from ttmt.conventions import sample_index, sample_index_clipped
from ttmt.detector import (CAP6, LOG6, NOCAP, ENERGY_ARGMAX, EnergyArgmax, Engine, Huber, ProjectedEnergy, RadialMap,
                           Reader, deflate, greedy_loop, mt_loop, peak_table, vectors)
from ttmt.kernel import LEGACY, Kernel
from ttmt.partition import Partition, PartitionError, SegmentReader, legacy_t_end
from ttmt.radiation import basis, bases, node_operators, operators, project, radiation_operator, stack_Q
from ttmt.station_terms import OnsetDelay, PeakDelay, anchored_peak_terms, require
from ttmt.traveltime import homogeneous

RNG = np.random.default_rng(20261018)


# ------------------------------------------------------------------------------------------- sample index (E52)
def test_sample_index_floor_is_legacy_truncation():
    x = np.array([0.0, 0.49, 0.5, 0.99, 1.0, 7.5, 12.999])
    assert sample_index(x, "floor").tolist() == [0, 0, 0, 0, 1, 7, 12]
    assert sample_index(x, "floor").tolist() == x.astype(int).tolist()


def test_sample_index_nearest():
    x = np.array([0.0, 0.49, 0.5, 0.99, 1.0, 7.5, 12.999])
    assert sample_index(x, "nearest").tolist() == [0, 0, 1, 1, 1, 8, 13]


def test_sample_index_fraction_and_bias():
    x = RNG.uniform(0, 600, 20000)
    i_f, f_f = sample_index_clipped(x, 640, "floor")
    i_n, f_n = sample_index_clipped(x, 640, "nearest")
    assert np.all((f_f >= 0) & (f_f < 1)) and np.all((f_n >= -0.5) & (f_n < 0.5))
    assert abs(f_f.mean() - 0.5) < 0.01 and abs(f_n.mean()) < 0.01        # floor reads half a sample early
    with pytest.raises(ValueError):
        sample_index(x, "ceil")


def test_sample_index_clipped():
    idx, fr = sample_index_clipped(np.array([-3.2, 700.7]), 640)
    assert idx.tolist() == [0, 639]
    assert np.allclose(fr, [-3.2, 61.7])


# ------------------------------------------------------------------------------------------- subspace
def _net(ns=6):
    a = np.linspace(0, 2 * np.pi, ns, endpoint=False)
    return np.stack([25 * np.cos(a), 25 * np.sin(a), np.zeros(ns)], axis=1)


def test_subspace_projection_idempotent_and_orthogonal():
    STA = _net()
    G = radiation_operator(np.array([1.0, -2.0, 7.0]), STA)
    Q = basis(G)
    assert Q.shape == (6 * len(STA), 6)
    assert np.allclose(Q.T @ Q, np.eye(6), atol=1e-12)
    v = RNG.normal(size=(5, G.shape[0]))
    p = project(Q, v)
    assert np.allclose(project(Q, p), p, atol=1e-12)                       # idempotent
    assert np.allclose((v - p) @ Q, 0.0, atol=1e-12)                       # residual orthogonal to the span
    assert np.allclose((v - p) @ G, 0.0, atol=1e-10)
    assert np.allclose(project(Q, G.T).T, G, atol=1e-12)                   # span(G) is fixed


def test_operators_vectorised_match_per_node():
    STA = _net(4)
    XN = np.array([[0.37, 3.37, 5.0], [-9.63, 12.37, 9.0]])
    Gv = operators(XN, STA)
    Gs, Qst = node_operators(XN, STA)
    assert np.allclose(Gv, Gs, rtol=1e-12, atol=1e-15)
    for g in range(2):
        assert np.allclose(np.abs(bases(Gv)[g].T @ Qst[g]), np.eye(6), atol=1e-8) or \
            np.allclose(project(Qst[g], Gv[g].T).T, Gv[g], atol=1e-12)
    assert stack_Q([np.ones((4, 2)), np.ones((4, 3))]).shape == (2, 4, 3)


# ------------------------------------------------------------------------------------------- radial maps and statistics
def test_radial_maps():
    v = np.array([[3.0, 4.0, 0.0], [30.0, 40.0, 0.0], [0.0, 0.0, 0.0]])
    n = np.linalg.norm
    assert np.allclose(n(CAP6(v), axis=1), [5.0, 6.0, 0.0])
    assert np.allclose(NOCAP(v), v)
    lg = LOG6(v)
    assert np.allclose(n(lg, axis=1) ** 2, 36 * np.log1p(n(v, axis=1) ** 2 / 36))
    assert np.allclose(lg[0] / n(lg[0]), v[0] / n(v[0]))                  # direction kept
    assert RadialMap.from_capv(None) == NOCAP and RadialMap.from_capv(6.0) == CAP6 and RadialMap.from_capv(-6.0) == LOG6
    with pytest.raises(ValueError):
        RadialMap("sqrt")


def test_block_statistic_influence_functions():
    r = np.array([0.5, 2.0, 2.79, 10.0, 100.0])
    cap, log, hub, ls = ProjectedEnergy(CAP6), ProjectedEnergy(LOG6), Huber(), ProjectedEnergy(NOCAP)
    assert np.allclose(ls.influence(r), r)
    assert np.allclose(cap.influence(r), np.where(r < 6, r, 0.0))
    assert np.all(hub.influence(r) <= hub.c if hub.c else hub.influence(r) <= np.sqrt(7.814727903251178) + 1e-12)
    assert np.allclose(log.influence(r), r / (1 + r * r / 36))
    for st in (cap, log, hub, ls):                                         # influence = d rho / d r
        h = 1e-6
        num = (st.rho(r + h) - st.rho(r - h)) / (2 * h)
        ok = np.abs(r - 6.0) > 1e-3
        assert np.allclose(num[ok], st.influence(r)[ok], rtol=1e-5, atol=1e-6), st.name
    assert np.all(np.diff(cap.rho(r)) >= 0) and cap.rho(1e9) == 18.0       # bounded influence of one block


def test_huber_equals_projection_without_outliers():
    """with every block below the Huber threshold, HUBER-N reduces to the least-squares projected energy."""
    STA = _net(5)
    G = radiation_operator(np.array([0.0, 0.0, 8.0]), STA)
    Q = basis(G)
    v = (Q @ RNG.normal(size=6)) * 0.05
    sh, _ = Huber().score(v[None, None], Q[None])
    sl, _ = ProjectedEnergy(NOCAP).score(v[None, None], Q[None])
    assert np.allclose(sh, sl, rtol=1e-6)


# ------------------------------------------------------------------------------------------- readers and deflation
def test_energy_argmax_reader():
    R = np.zeros((2, 3, 64))
    R[0, 1, 20] = 5.0
    R[1, 0, 40] = 3.0
    idx = np.array([[[17, 44], [25, 36]]])                                  # (1, 2, ns)
    out = ENERGY_ARGMAX.read(R, idx)
    assert out.tolist() == [[[20, 40], [21, 40]]]                          # no energy within +-4 of 25: first sample
    assert EnergyArgmax(1).read(R, idx).tolist() == [[[16, 43], [24, 35]]]  # outside +-1: first sample of the window
    assert isinstance(ENERGY_ARGMAX, Reader)
    assert np.array_equal(peak_table(R)[np.arange(2)[None, None], idx], out)


def _toy_engine(STA, taus, rounding="floor"):
    XN = np.array([[x + 0.37, y + 0.37, z] for x in (-9.0, 0.0, 9.0) for y in (-9.0, 0.0, 9.0) for z in (5.0,)])
    tt, _d = homogeneous(XN, STA, 6.0, 6.0 / 1.73)
    Gs, Qst = node_operators(XN, STA)
    return Engine(XN, STA, tt, taus, PeakDelay.zeros(len(STA)), Gs, Qst, rounding=rounding)


def _synth_atom(E, gi, tj, m):
    """noise-free matched-filter output of one atom: G m at the predicted samples, shaped by the deflation response."""
    R = np.zeros((E.ns, 3, 640))
    pred = (E.G[gi] @ m).reshape(2, E.ns, 3)
    return deflate(R, -pred, E.IDX[gi, tj], LEGACY), pred


def test_deflation_removes_an_exact_atom():
    STA = _net(6)
    E = _toy_engine(STA, np.round(np.arange(4.0, 8.0, 0.1), 10))
    MFn, _ = _synth_atom(E, 4, 15, RNG.normal(size=6) * 2000)
    atoms = greedy_loop(E, MFn, 1, ProjectedEnergy(NOCAP))
    assert len(atoms) == 1
    # origins within the +-4-sample reading read the same peaks: the first of the tied origins is taken
    assert np.allclose(atoms[0]["x"], E.XN[4]) and abs(atoms[0]["tau"] - E.taus[15]) <= 0.25 + 1e-9
    tj = int(np.argmin(np.abs(E.taus - atoms[0]["tau"])))
    R = MFn.copy()
    tpk = ENERGY_ARGMAX.read(R, E.IDX)
    assert np.array_equal(tpk[4, tj], E.IDX[4, 15])
    vf = vectors(R, tpk).reshape(E.N, E.T, -1)
    _s, aux = ProjectedEnergy(NOCAP).score(vf, E.Q)
    deflate(R, ProjectedEnergy(NOCAP).deflation(E, 4, tj, vf, aux), tpk[4, tj])
    assert np.abs(R).max() < 1e-9 * np.abs(MFn).max()


def test_mt_loop_legacy_signature_and_exclusion():
    STA = _net(6)
    E = _toy_engine(STA, np.round(np.arange(4.0, 8.0, 0.1), 10))
    MFn, _ = _synth_atom(E, 0, 10, RNG.normal(size=6) * 50)
    MFn = MFn + RNG.normal(size=MFn.shape)
    a1 = mt_loop(E, MFn, 3, capv=None)
    a2 = greedy_loop(E, MFn, 3, ProjectedEnergy(NOCAP))
    assert [a["s"] for a in a1] == [a["s"] for a in a2]
    xs = np.array([a["x"][:2] for a in a1])
    for i in range(len(xs)):
        for j in range(i):
            assert np.linalg.norm(xs[i] - xs[j]) > 5.0                        # epicentral exclusion
    assert all(a1[k]["run"] == min(a["s"] for a in a1[:k + 1]) for k in range(len(a1)))


def test_engine_rounding_convention():
    STA = _net(4)
    taus = np.array([5.0, 5.1])
    Ef, En = _toy_engine(STA, taus, "floor"), _toy_engine(STA, taus, "nearest")
    assert np.all((En.IDX - Ef.IDX >= 0) & (En.IDX - Ef.IDX <= 1))
    assert np.all(np.abs(En.frac) <= 0.5 + 1e-12)


# ------------------------------------------------------------------------------------------- kernel (E51)
def test_kernel_ties_filter_and_deflation():
    k2, k3 = Kernel(2.0), Kernel(3.0)
    assert k2.ac.argmax() == k2.ac_half and abs(k2.ac.max() - 1.0) < 1e-15
    ref = np.convolve(k2.taps, k2.taps[::-1], mode="full")                   # autocorrelation of the 2 Hz kernel
    assert np.allclose(np.convolve(np.r_[np.zeros(16), k2.taps, np.zeros(16)][16:-16], k2.taps[::-1], "same"),
                       np.convolve(k2.taps, k2.taps[::-1], "same"))
    assert ref.argmax() == len(ref) // 2
    assert not np.allclose(k2.ac, k3.ac)                                     # changing fc changes the deflation too
    assert Kernel(3.0) == k3 and Kernel(3.0, pulse_fc=2.0) != k3              # a mismatch must be explicit
    with pytest.raises(Exception):
        k3.fc = 2.0                                                          # frozen: cannot be split at run time
    D = np.zeros((1, 3, 640))
    D[0, 0, 300] = 1.0
    MF, _sd = k2.matched(D)
    assert np.allclose(MF[0, 0, 300 - 16:300 + 17], k2.taps)                 # symmetric kernel


def test_kernel_ac_at():
    k = LEGACY
    assert k.ac_at(0) == 1.0 and k.ac_at(k.ac_half + 1) == 0.0
    assert np.allclose(k.ac_at(np.array([-3, 3])), k.ac[[k.ac_half - 3, k.ac_half + 3]])


# ------------------------------------------------------------------------------------------- station terms (E54)
def test_station_term_typing():
    STA = _net(4)
    taus = np.array([5.0])
    XN = np.array([[0.37, 0.37, 5.0]])
    tt, _ = homogeneous(XN, STA, 6.0, 6.0 / 1.73)
    Gs, Qst = node_operators(XN, STA)
    with pytest.raises(TypeError):
        Engine(XN, STA, tt, taus, OnsetDelay.zeros(4), Gs, Qst)
    with pytest.raises(TypeError):
        Engine(XN, STA, tt, taus, np.zeros((4, 2)), Gs, Qst)
    E = Engine(XN, STA, tt, taus, PeakDelay.zeros(4), Gs, Qst)
    with pytest.raises(TypeError):
        E.set_times(tt, OnsetDelay.zeros(4))
    with pytest.raises(TypeError):
        anchored_peak_terms([], STA, PeakDelay.zeros(4), None, 6.0, 3.5)
    from ttmt.r1 import score_candidates
    with pytest.raises(TypeError):
        score_candidates([], STA, [], PeakDelay.zeros(4), 0.5, None)
    assert require(OnsetDelay.zeros(2), OnsetDelay).shape == (2, 2)
    with pytest.raises(ValueError):
        PeakDelay(np.zeros(3))


def test_peak_terms_shift_the_engine():
    STA = _net(4)
    taus = np.array([5.0])
    XN = np.array([[0.37, 0.37, 5.0]])
    tt, _ = homogeneous(XN, STA, 6.0, 6.0 / 1.73)
    Gs, Qst = node_operators(XN, STA)
    v = np.zeros((4, 2))
    v[2, 1] = 0.25
    E0 = Engine(XN, STA, tt, taus, PeakDelay.zeros(4), Gs, Qst)
    E1 = Engine(XN, STA, tt, taus, PeakDelay(v), Gs, Qst)
    d = E1.IDX - E0.IDX
    assert d[0, 0, 1, 2] == 4 and np.count_nonzero(d) == 1


# ------------------------------------------------------------------------------------------- hit rules
def _atoms(rows):
    return [{"x": np.array([x, y, 5.0]), "tau": t, "s": s} for x, y, t, s in rows]


def test_hit_rule_legacy_loc_has_no_time_test_and_duplicates_neutral():
    truth, taus = [np.array([0.0, 0.0, 5.0]), np.array([20.0, 0.0, 5.0])], [6.0, 6.5]
    at = _atoms([(1.0, 0.0, 9.0, 10.0), (0.5, 0.0, 6.0, 9.0), (30.0, 0.0, 6.0, 8.0)])
    assert S.label(at, truth, taus, None) == [0, 0, -1]                     # 3 s origin error still a hit under LOC
    assert S.label(at, truth, taus, 0.5) == [-1, 0, -1]                     # LOC+T
    rec = S.record(at, S.label(at, truth, taus, None))
    th, H, F = S.froc_curve([rec])
    assert H.tolist() == [1, 1, 1] and F.tolist() == [0, 0, 1]              # the duplicate is neutral


def test_hit_rule_one_to_one():
    truth, taus = [np.array([0.0, 0.0, 5.0]), np.array([3.0, 0.0, 5.0])], [6.0, 6.2]
    at = _atoms([(0.2, 0.0, 6.1, 5.0), (0.1, 0.0, 6.0, 9.0), (0.0, 0.1, 6.05, 7.0)])
    # record order (loop order): atom 0 takes event 0, atom 1 the nearest unmatched (event 1), atom 2 is false
    assert S.label(at, truth, taus, 0.5, match="one_to_one") == [0, 1, -1]
    # raw-score order: atom 1 (9) takes event 0, atom 2 (7) event 1, atom 0 (5) is false
    assert S.label(at, truth, taus, 0.5, match="one_to_one", order="score") == [-1, 0, 1]
    assert S.label(at, truth, taus, 0.5, match="legacy") == [0, 0, 0]
    at2 = _atoms([(0.0, 0.0, 6.0, 9.0), (0.1, 0.0, 6.0, 8.0)])
    one = [S.record(at2, S.label(at2, [truth[0]], [6.0], 0.5, match="one_to_one"))]
    leg = [S.record(at2, S.label(at2, [truth[0]], [6.0], 0.5, match="legacy"))]
    assert S.froc_curve(one)[2].tolist() == [0, 1] and S.froc_curve(leg)[2].tolist() == [0, 0]   # duplicate false
    with pytest.raises(ValueError):
        S.label(at2, truth, taus, 0.5, match="nearest")
    with pytest.raises(ValueError):
        S.label(at2, truth, taus, 0.5, match="one_to_one", order="time")


def test_one_to_one_record_order_keeps_the_prefix_property():
    """loop order a1 s=5 (eligible for A and B, nearest A), a2 s=3, a3 s=6 (eligible for A only). At threshold 4 only
    a1 is admitted (the record is the running minimum 5, 3, 3) and it must hit A."""
    A, B = np.array([0.0, 0.0, 5.0]), np.array([4.0, 0.0, 5.0])
    at = _atoms([(1.0, 0.0, 6.0, 5.0), (20.0, 0.0, 6.0, 3.0), (-2.0, 0.0, 6.0, 6.0)])
    lab = S.label(at, [A, B], [6.0, 6.0], None, match="one_to_one")
    assert lab == [0, -1, -1]
    rec = S.record(at, lab)
    assert [r[0] for r in rec] == [5.0, 3.0, 3.0]
    th, H, F = S.froc_curve([rec])
    assert H.tolist() == [1, 1] and F.tolist() == [0, 2]                     # at 5: a1 hits A
    assert S.label(at, [A, B], [6.0, 6.0], None, match="one_to_one", order="score") == [1, -1, 0]   # breaks it


def _brute_froc(windows, tol_t):
    """thresholding each window's loop-order prefix and matching that prefix one to one from scratch."""
    levels = sorted({min(a["s"] for a in at[:k + 1]) for at, _t, _ta in windows for k in range(len(at))}, reverse=True)
    out = []
    for tau in levels:
        h = f = 0
        for at, truth, taus in windows:
            run = np.minimum.accumulate([a["s"] for a in at]) if at else []
            pre = [a for a, r in zip(at, run) if r >= tau]
            lab = S.label(pre, truth, taus, tol_t, match="one_to_one")
            h += len({x for x in lab if x >= 0})
            f += sum(1 for x in lab if x < 0)
        out.append((tau, h, f))
    return out


@pytest.mark.parametrize("tol_t", [None, 0.5])
def test_one_to_one_froc_equals_brute_force_thresholding(tol_t):
    rng = np.random.default_rng(99)
    windows, recs = [], []
    for _ in range(60):
        k = int(rng.integers(0, 7))
        truth = [np.array([rng.uniform(-8, 8), rng.uniform(-8, 8), 5.0]) for _ in range(int(rng.integers(0, 4)))]
        taus = [float(rng.uniform(5.5, 7.5)) for _ in truth]
        at = [{"x": np.array([rng.uniform(-10, 10), rng.uniform(-10, 10), 5.0]), "tau": float(rng.uniform(5, 8)),
               "s": float(rng.integers(1, 30))} for _ in range(k)]                  # non-monotone, with ties
        windows.append((at, truth, taus))
        recs.append(S.record(at, S.label(at, truth, taus, tol_t, match="one_to_one")))
    th, H, F = S.froc_curve(recs)
    brute = _brute_froc(windows, tol_t)
    assert [round(b[0], 9) for b in brute] == [round(x, 9) for x in th]
    assert [b[1] for b in brute] == H.tolist() and [b[2] for b in brute] == F.tolist()


def test_one_to_one_matches_froc_strict_of_the_other_thread():
    """cross-check with froc_strict.relabel (PR #4, code/froc_strict.py; before the merge, branch
    claude/project-thread-6f81ru), which matches in record order."""
    import importlib.util
    import subprocess
    import tempfile
    local = Path(__file__).parents[1] / "code" / "froc_strict.py"
    if local.is_file():
        src = local.read_text()
    else:
        try:
            src = subprocess.run(["git", "show", "origin/claude/project-thread-6f81ru:code/froc_strict.py"],
                                 capture_output=True, text=True, check=True, cwd=str(Path(__file__).parents[1])).stdout
        except (subprocess.CalledProcessError, FileNotFoundError):
            pytest.skip("froc_strict.py not available")
    with tempfile.TemporaryDirectory() as d:
        f = Path(d) / "froc_strict.py"
        f.write_text(src)
        spec = importlib.util.spec_from_file_location("froc_strict", f)
        FS = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(FS)
    rng = np.random.default_rng(5)
    n = 0
    for _ in range(300):
        truth = [np.array([rng.uniform(-8, 8), rng.uniform(-8, 8), 5.0]) for _ in range(int(rng.integers(0, 5)))]
        taus = [float(rng.uniform(5.5, 7.5)) for _ in truth]
        at = [{"x": np.array([rng.uniform(-10, 10), rng.uniform(-10, 10), 5.0]), "tau": float(rng.uniform(5, 8)),
               "s": float(rng.uniform(0, 30))} for _ in range(int(rng.integers(0, 8)))]
        locs = [list(a["x"]) + [a["tau"]] for a in at]
        locs = [[l[0], l[1], l[2], l[3]] for l in locs]
        for tol_t in (None, 0.5):
            assert S.label(at, truth, taus, tol_t, match="one_to_one") == FS.relabel(locs, truth, taus, tol_t), (at, truth)
            assert S.label(at, truth, taus, tol_t, match="legacy") == FS.relabel(locs, truth, taus, tol_t, one_to_one=False)
            n += 1
    assert n == 600


def test_catalogue_rule_prefers_targets_and_one_to_one():
    cat = [{"x": 0.0, "y": 0.0, "t": 10.0, "M": 1.0}, {"x": 2.0, "y": 0.0, "t": 10.2, "M": 2.5}]
    det = [[0.1, 0.0, 5.0, 10.0, 3.0], [0.0, 0.0, 5.0, 10.1, 2.0]]
    assert S.label_catalogue(det, cat, 0.5, 2.0).tolist() == [1, 1]
    assert S.label_catalogue(det, cat, 0.5, 2.0, match="one_to_one").tolist() == [1, 0]


# ------------------------------------------------------------------------------------------- exact FROC
def test_exact_froc_monotone_and_matches_threshold_sweep():
    recs = []
    for _ in range(40):
        k = RNG.integers(0, 6)
        s = np.sort(RNG.uniform(0, 10, k))[::-1]
        recs.append([(float(x), int(RNG.integers(-1, 3))) for x in s])
    th, H, F = S.froc_curve(recs)
    assert np.all(np.diff(th) < 0) and np.all(np.diff(H) >= 0) and np.all(np.diff(F) >= 0)
    thr, HH, FF = S.exact_HF(recs)
    for j, t in enumerate(thr):                                             # brute force at every threshold
        hits = sum(len({h for s_, h in r if s_ >= t and h >= 0}) for r in recs)
        fals = sum(sum(1 for s_, h in r if s_ >= t and h < 0) for r in recs)
        assert HH[:, j].sum() == hits and FF[:, j].sum() == fals
        k = np.where(np.isclose(th, t))[0]
        if len(k):
            assert H[k[0]] == hits and F[k[0]] == fals
    nte = 60
    r_lo, _ = S.recall_caps(recs, nte, (0.05,))
    r_hi, _ = S.recall_caps(recs, nte, (0.5,))
    assert r_lo[0] <= r_hi[0]


def test_greedy_record_prefix_property():
    grid = [np.array([x, y, 5.0]) for x in range(0, 40, 3) for y in range(0, 40, 3)]
    field = RNG.uniform(0, 1, len(grid))
    rec = S.greedy_record(grid, field, [np.array([10.0, 10.0, 5.0])])
    assert [r[0] for r in rec] == sorted([r[0] for r in rec], reverse=True)


# ------------------------------------------------------------------------------------------- partition (E56) and shifts
def _reader(shifts=None, strict=True):
    dt = 1.0 / 16
    X = np.arange(4 * 3 * 1600, dtype=float).reshape(4, 3, 1600)            # 100 s at 16 Hz
    part = Partition({"cal": (0.0, 50.0), "eval": (50.0, 100.0)})
    return SegmentReader(X, dt, part, "cal", shifts=shifts, strict=strict), X, part


def test_partition_reader_raises_across_the_boundary():
    rd, X, part = _reader()
    W = rd.window(10.0, 640)
    assert np.array_equal(W, X[:, :, 160:800])
    rd.window(10.0, 640)                                                     # ends at 50 s exactly: allowed
    with pytest.raises(PartitionError):
        rd.window(10.0 + 1.0 / 16, 640)                                     # one sample into eval
    with pytest.raises(PartitionError):
        rd.window(-1.0, 32)
    assert legacy_t_end(part, "cal") == 90.0                                # the legacy margin crossed it (E56)
    ev = SegmentReader(X, 1.0 / 16, part, "eval")
    with pytest.raises(PartitionError):
        ev.window(49.0, 32)
    with pytest.raises(ValueError):
        Partition({"a": (0.0, 10.0), "b": (5.0, 20.0)})


def test_zero_shift_is_identity():
    rd0, X, _ = _reader()
    rdz, _, _ = _reader(shifts=np.zeros(4, int))
    for t in (0.0, 7.25, 10.0):
        assert np.array_equal(rd0.window(t, 640), rdz.window(t, 640))


def test_shifts_stay_inside_the_segment():
    shifts = np.array([0, 5, -37, 799])                                     # 799 samples ~ one full segment
    rd, X, _ = _reader(shifts=shifts)
    seg = X[:, :, 0:800]
    W = rd.window(0.0, 800)
    for i, s in enumerate(shifts):
        assert np.array_equal(W[i], np.roll(seg[i], s, axis=-1))
    assert W.min() >= seg.min() and W.max() <= seg.max()                    # no eval sample is ever read
    assert set(np.unique(W)) <= set(np.unique(seg))
    rd_ns, _, _ = _reader(shifts=shifts, strict=False)
    with pytest.raises(PartitionError):
        rd_ns.window(20.0, 640)                                             # shifted reads are always guarded
    with pytest.raises(ValueError):
        _reader(shifts=np.zeros(3, int))


# ------------------------------------------------------------------------------------------- legacy harness
def test_legacy_install_keeps_explicit_patches():
    """run_script() calls install() again; the defaults must not drop a constant registered explicitly (this dropped
    r1_lib.HERE once and let PyOcto write its model file into code/r1/)."""
    from ttmt import legacy
    f = legacy.install(patches={"r1_lib": {"HERE": "/elsewhere"}})
    legacy.install()
    assert f.patches["r1_lib"]["HERE"] == "/elsewhere" and "DATA" in f.patches["r1_lib"]
    legacy.install(patches={"r1_lib": {"DATA": "/data"}})
    assert f.patches["r1_lib"] == {"HERE": "/elsewhere", "DATA": "/data"}
    del f.patches["r1_lib"]["HERE"]
