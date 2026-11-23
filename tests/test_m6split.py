"""Tests of ttmt.v2.m6split: the M6-SPLIT2b fix (review of PR #75). Both hypotheses refit their isolated events on
copies, so the master isolated fit is unchanged by the pair test, and converge_iso returns the RSS of its final state."""
import numpy as np

from test_m6 import _engine
from test_v2 import LAG_S, K, _q, _scene
from ttmt.v2 import m6, m6split as ms


def _setup():
    STA, _srcs, F = _scene(1, 0.0, 11)
    E = _engine(STA)
    lag = np.full((8, 2), LAG_S / K.dt)
    q = _q([F[0][2] - LAG_S], 8)
    Y = F[0][0]
    Yw = F[0][0] * 0.8
    aw = m6.adcg_loop(E, Yw, q, 1, lag, thr=None)
    for a in aw:
        a.update(anchor=a["T"].copy(), bound=np.full(a["T"].shape, 2.0), phi=np.zeros(E.ns), dl=np.zeros(E.ns))
    atoms = m6.adcg_loop(E, Y, q, 1, lag, thr=None)
    for a in atoms:
        a.update(anchor=a["T"].copy(), bound=np.full(a["T"].shape, 2.0), phi=np.zeros(E.ns), dl=np.zeros(E.ns))
    resp = ms._Resp(K)
    frames = ms.Frames(E.XN, E.STA)
    return E, Y, atoms, [(np.asarray(Yw, float), aw)], resp, frames


def _snap(iso):
    return [[(a["T"].copy(), a["m"].copy(), a["phi"].copy(), a["dl"].copy()) for a in aw] for _Yw, aw in iso]


def test_both_hypotheses_leave_the_master_isolated_fit_unchanged():
    E, Y, atoms, iso, resp, frames = _setup()
    dmax = ms.DMAX_S / K.dt
    before = _snap(iso)
    ms.free_one_event(Y, [], atoms[0], E.G, resp, frames, "station", dmax, iso=iso, converge=True)
    iso_t = ms.copy_iso(iso)
    trial = ms._copy_atoms(atoms)
    ms.joint_refit(Y, trial, E.G, resp, frames, "station", dmax, iso=iso_t)
    ms.converge_iso(iso_t, E.G, resp, frames)
    after = _snap(iso)
    for wb, wa in zip(before, after):
        for ab, aa in zip(wb, wa):
            for xb, xa in zip(ab, aa):
                assert np.array_equal(xb, xa)


def test_converge_iso_returns_the_rss_of_its_final_state():
    E, _Y, _atoms, iso, resp, frames = _setup()
    out = ms.converge_iso(iso, E.G, resp, frames)
    assert np.isclose(out, ms.iso_rss(iso, E.G, resp, frames), rtol=1e-12, atol=0)


def test_free_one_event_with_atom_start_nests_the_split_atom():
    """M6-SPLIT audit (23 Nov 2026): on noise-free data made by the split model itself, the free one-event alternative
    (a larger class) must fit at least as well as the atom; without atom_start its time search could stay at the
    unsplit optimum and fit much worse, which biased the pair gain upwards."""
    E, _Y, atoms, _iso, resp, frames = _setup()
    a = ms._copy_atoms(atoms)[0]
    a["anchor"], a["bound"] = a["T"].copy(), np.full(a["T"].shape, 2.0)
    r = np.random.default_rng(3)
    a["phi"], a["dl"] = r.uniform(0, np.pi, E.ns), r.uniform(0.05, 0.2, E.ns) / K.dt
    Y = ms.model([a], E.G, (E.ns, 3, _Y.shape[-1]), resp, frames)
    a0 = ms._copy_atoms([a])[0]
    a0["anchor"], a0["bound"] = a["anchor"], a["bound"]
    dmax = ms.DMAX_S / K.dt
    rss = ms.free_one_event(Y, [], a0, E.G, resp, frames, "fixed", dmax, atom_start=True)
    assert rss <= 1e-6 * float((Y * Y).sum())
