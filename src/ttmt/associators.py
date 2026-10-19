"""Configuration grids of the pick-based comparators (task1r_baselines, frozen in Task 2A). The associators themselves
(GaMMA, PyOcto) and the picker (PhaseNet) live in the picks environment (environments/DEVIATIONS.md)."""
import itertools

GAMMA_GRID = [dict(min_picks=m, sigma=s) for m, s in itertools.product([3, 4, 5, 6, 8], [1.0, 2.0])]
PYOCTO_GRID = [dict(n_picks=n, tol=t) for n, t in itertools.product([3, 4, 6, 8, 10], [0.5, 1.0, 2.0])]
