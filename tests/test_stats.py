"""Tests of ttmt.stats (ADR-0008, phase 6 analysis plan) and of the sample-size tool of experiments/exploratory/P6-SAP."""
import csv
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

from ttmt import stats as st

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "experiments" / "exploratory" / "P6-SAP"))
import sample_size as ss  # noqa: E402


class HitRule(unittest.TestCase):
    def setUp(self):
        # a close pair: same epicentre, origins 0.30 s apart (separation below the 0.5 s tolerance)
        self.ev_xy, self.ev_t = np.array([[0.0, 0.0], [0.0, 0.0]]), np.array([10.0, 10.30])

    def test_duplicate_of_one_event_is_not_two_hits(self):
        # two detections of the first event, 0.05 s apart; plain LOC+T one-to-one would hit both events
        lab = st.label_voronoi([2.0, 1.0], [[0, 0], [0, 0]], [10.0, 10.05], self.ev_xy, self.ev_t)
        self.assertEqual(lab.tolist(), [0, st.FALSE])

    def test_resolved_pair_gets_two_hits(self):
        lab = st.label_voronoi([2.0, 1.0], [[0, 0], [0, 0]], [10.02, 10.28], self.ev_xy, self.ev_t)
        self.assertEqual(lab.tolist(), [0, 1])

    def test_isolated_event_reduces_to_loc_t(self):
        lab = st.label_voronoi([1.0, 0.5], [[3.0, 0], [0, 0]], [5.4, 7.0], [[0, 0]], [5.0])
        self.assertEqual(lab.tolist(), [0, st.FALSE])

    def test_duplicator_control_gains_nothing(self):
        # a detector that hits only the first event of each pair; its duplicator (+-sep/2) must not resolve the pair
        rng = np.random.default_rng(4)
        seps = rng.uniform(0.1, 1.0, 200)
        ev_t = np.c_[np.arange(200) * 100.0, np.arange(200) * 100.0 + seps].ravel()
        ev_xy = np.zeros((400, 2))
        det_t = ev_t[0::2] + rng.normal(0, 0.02, 200)
        scores = rng.random(200)
        base = st.label_voronoi(scores, np.zeros((200, 2)), det_t, ev_xy, ev_t)
        s2, xy2, t2 = st.duplicator(scores, np.zeros((200, 2)), det_t, seps / 2)
        dup = st.label_voronoi(s2, xy2, t2, ev_xy, ev_t)
        hits = lambda lab: np.isin(np.arange(400), lab[lab >= 0])
        pairs = lambda h: int((h[0::2] & h[1::2]).sum())
        self.assertEqual(pairs(hits(base)), 0)
        # at most the pairs near 0.1 s whose 0.02 s jitter is a fifth of the separation (1 of 200 with this seed)
        self.assertLessEqual(pairs(hits(dup)), pairs(hits(base)) + 2)
        # with plain nearest-event eligibility (dominance 1) the duplicator resolves about half of the pairs
        dup1 = st.label_voronoi(s2, xy2, t2, ev_xy, ev_t, dominance=1.0)
        self.assertGreater(pairs(hits(dup1)), 60)


class OperatingPoint(unittest.TestCase):
    def test_budget_and_ties(self):
        scores = [9, 8, 8, 7, 6, 5]
        labels = [0, st.FALSE, 1, st.FALSE, 2, st.FALSE]
        # cap 1 false: the tie group at 8 holds one false, the group at 7 adds a second -> threshold 8
        self.assertEqual(st.operating_point(scores, labels, hours=1.0, budget=1.0), 8.0)
        self.assertEqual(st.operating_point(scores, labels, hours=1.0, budget=2.0), 6.0)
        self.assertEqual(st.operating_point(scores, labels, hours=1.0, budget=3.0), 5.0)
        self.assertEqual(st.operating_point([5], [st.FALSE], hours=1.0, budget=0.5), np.inf)

    def test_hits_and_false(self):
        scores, labels = [9, 8, 7, 6], [0, st.FALSE, st.NEUTRAL, 2]
        self.assertEqual(st.hits_at(scores, labels, 7, 3).tolist(), [True, False, False])
        self.assertEqual(st.hits_at(scores, labels, 6, 3).tolist(), [True, False, True])
        self.assertEqual(st.false_count(scores, labels, 6), 1)


class PairedTests(unittest.TestCase):
    def test_signflip_equals_exact_mcnemar_with_one_target_per_block(self):
        for b, c in [(5, 1), (3, 3), (7, 0), (0, 4)]:
            d = np.r_[np.ones(b), -np.ones(c)]
            self.assertAlmostEqual(st.signflip_pvalue(d), st.exact_mcnemar(b, c), places=12)

    def test_monte_carlo_close_to_exact(self):
        rng = np.random.default_rng(1)
        d = rng.integers(-3, 4, 14).astype(float)
        exact = st.signflip_pvalue(d, exact_max=16)
        mc = st.signflip_pvalue(d, exact_max=0, n_flip=200_000, rng=rng)
        self.assertLess(abs(exact - mc), 0.005)

    def test_clustered_mcnemar_reduces_to_mcnemar(self):
        d = np.r_[np.ones(9), -np.ones(4)]
        self.assertAlmostEqual(st.clustered_mcnemar_z(d), (9 - 4) / np.sqrt(13))

    def test_level_under_clustered_null(self):
        r = ss.run_design(200, delta_r=0.0, delta_c=0.0, icc=0.3, rho=0.6, m=4, k=0.5, sims=400, n_flip=999, seed=7)
        self.assertLessEqual(r["power_H1"], 0.025 + 2.6 * np.sqrt(0.025 * 0.975 / 400))

    def test_h2_zero_discordance_is_not_non_inferior(self):
        """30 blocks of one target, none discordant, cannot exclude a difference below -0.05 (review of PR #13)."""
        z = np.zeros(30)
        self.assertGreater(st.ni_score_pvalue(z, z, np.ones(30), 0.05), 0.05)
        self.assertLess(st.ni_score_pvalue(np.zeros(1000), np.zeros(1000), np.ones(1000), 0.05), 1e-6)
        c = z.copy(); c[0] = 1
        self.assertGreater(st.ni_score_pvalue(z, c, np.ones(30), 0.05), 0.2)

    def test_h2_size_at_least_favourable_configuration(self):
        """H2 (score test) at 0.025 with one comparator at the margin, 30 blocks, the worst simulated configuration of
        part 3 (recall 0.85, 10 targets per block, rho 0.9, icc 0.2)."""
        r = ss.run_design(30, m=10, k=0.5, p0=0.85, delta_r=0.05, delta_c=-0.05, rho=0.9, icc=0.2, n_comp=1,
                          sims=2000, n_flip=99, seed=41)
        self.assertLessEqual(r["power_H2"], 0.025 + 1.5 * np.sqrt(0.025 * 0.975 / 2000))

    def test_tau_keeps_marginal_recall(self):
        rng = np.random.default_rng(8)
        _, reg, ha, _ = ss.simulate_hits(rng, 20000, 3.0, 0, 0.5, 0.40, 0.10, -0.05, 0.7, 0.1, 0.6, 1)
        self.assertAlmostEqual(ha[reg].mean(), 0.50, delta=0.01)
        self.assertAlmostEqual(ha[~reg].mean(), 0.35, delta=0.01)

    def test_power_increases_with_blocks(self):
        a = ss.run_design(60, sims=200, seed=3)["power_H1"]
        b = ss.run_design(400, sims=200, seed=3)["power_H1"]
        self.assertLess(a, b)


class Multiplicity(unittest.TestCase):
    def test_rules(self):
        self.assertEqual(st.iut([0.01, 0.03]), 0.03)
        self.assertEqual(st.fixed_sequence([0.01, 0.2, 0.001], 0.025), [True, False, False])
        self.assertEqual(st.holm([0.01, 0.04, 0.02], 0.05), [True, True, True])
        self.assertEqual(st.holm([0.03, 0.04], 0.05), [False, False])


class Interval(unittest.TestCase):
    def test_bootstrap_flags_few_blocks(self):
        d, n = np.array([1.0, 0, 0, 0, 0]), np.array([3.0, 3, 3, 3, 3])
        r = st.block_bootstrap(st.recall_diff_stat(d, n), [np.arange(5)], n_boot=2000, rng=np.random.default_rng(0))
        self.assertTrue(r["indicative"])
        self.assertEqual(r["lo"], 0.0)
        self.assertGreater(r["atom_lo"], 0.05)

    def test_budget_check_level(self):
        """At a true rate equal to the null mean (1 + tol) * budget the check passes at most about 2.5-5 % of the
        time, with 40 blocks of unequal length and overdispersed counts."""
        rng = np.random.default_rng(5)
        passed = 0
        for _ in range(2000):
            h = rng.uniform(0.3, 1.0, 40)
            f = rng.negative_binomial(2, 2 / (2 + 36 * h))
            passed += st.budget_check(f, h, 30.0, tol=0.2)["pass"]
        self.assertLess(passed / 2000, 0.06)

    def test_budget_check_passes_at_budget_with_enough_hours(self):
        rng = np.random.default_rng(6)
        hrs = ss.hours_needed(30.0, phi=1.0)
        nb = 40
        passed = 0
        for _ in range(2000):
            h = np.full(nb, hrs / nb)
            passed += st.budget_check(rng.poisson(30 * h), h, 30.0)["pass"]
        # the floor phi >= 1 on an estimated dispersion makes the closed form slightly optimistic (about 0.75)
        self.assertTrue(0.70 <= passed / 2000 <= 0.85, passed / 2000)


class SampleSize(unittest.TestCase):
    def test_phi2(self):
        rng = np.random.default_rng(2)
        z = rng.multivariate_normal([0, 0], [[1, 0.6], [0.6, 1]], 400_000)
        mc = np.mean((z[:, 0] < -0.25) & (z[:, 1] < 0.4))
        self.assertAlmostEqual(ss.phi2(-0.25, 0.4, 0.6), mc, delta=0.003)

    def test_connor_matches_p6_power_grid(self):
        # P6-power (alpha 0.05, simulated exact McNemar): 80 % power at +0.10 needs 200 targets with rho 0.7
        n = ss.connor_n(0.40, 0.10, 0.7, alpha=0.05)
        self.assertTrue(150 <= n <= 200, n)

    def test_estimate_recovers_parameters(self):
        rng = np.random.default_rng(11)
        block, reg, ha, hc = ss.simulate_hits(rng, 3000, 4.0, 0, 1.0, 0.4, 0.0, 0.0, 0.7, 0.15, 0.0, 1)
        with tempfile.NamedTemporaryFile("w", suffix=".csv", delete=False, newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["block", "regime", "A", "B"])
            for row in zip(block, reg.astype(int), ha.astype(int), hc[0].astype(int)):
                w.writerow(row)
        e = ss.estimate(fh.name)
        self.assertAlmostEqual(e["rho"]["A|B"], 0.7, delta=0.04)
        self.assertAlmostEqual(e["icc"]["A"], 0.15, delta=0.05)


if __name__ == "__main__":
    unittest.main()
