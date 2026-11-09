"""M7 time-slide null (ttmt.v2.m7): exact on independent stationary stations, broken by a shared transient."""
import numpy as np

from ttmt.v2.m7 import pvalue, slide, time_slide_pvalue


def _beam(w):
    """a zero-moveout beam: max over time of the sum over stations of the energy (sensitive to coherence)."""
    return float((w * w).sum(axis=1).sum(axis=0).max())


def test_slide_keeps_each_station():
    rng = np.random.default_rng(0)
    w = rng.normal(size=(4, 3, 64))
    s, off = slide(w, rng)
    for i in range(4):
        assert np.allclose(np.roll(w[i], int(off[i]), axis=-1), s[i])
        assert np.allclose(np.sort(w[i].ravel()), np.sort(s[i].ravel()))


def test_pvalue_nonfinite_counts_as_rejection():
    assert pvalue(np.nan, np.zeros(9)) == 0.1
    assert pvalue(np.inf, np.zeros(9)) == 0.1
    assert pvalue(1.0, np.array([0.0, 2.0, 1.0])) == 0.75


def test_level_on_iid_null():
    rng = np.random.default_rng(1)
    p = np.array([time_slide_pvalue(rng.standard_t(4, size=(5, 3, 128)), _beam, 49, rng) for _ in range(800)])
    for a in (0.04, 0.1, 0.2):
        assert np.mean(p <= a) <= a + 3 * np.sqrt(a * (1 - a) / len(p)), (a, np.mean(p <= a))


def test_shared_transient_breaks_the_null():
    """a burst at the same time on every station (a coherent transient, not an earthquake) is rejected: the slides
    test 'no coherent source', so a coherent non-seismic transient is a detection for D.3."""
    rng = np.random.default_rng(2)
    p = []
    for _ in range(300):
        w = rng.normal(size=(5, 3, 128))
        k = rng.integers(128)
        w[:, :, k] += 3.0 * rng.normal(size=(5, 3))
        p.append(time_slide_pvalue(w, _beam, 49, rng))
    assert np.mean(np.array(p) <= 0.05) > 0.5
