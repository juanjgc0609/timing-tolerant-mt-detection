"""Search spaces and the shared random search of the ATLAS-0 tuning (card, section "Tuning"; ADR-0008 item 9).

Every method gets the same budget: N_TRIALS configurations, trial i drawn with SEEDS[i] by :func:`draw`, that is with
``numpy.random.default_rng(ttmt.atlas.seeds.tuning_seeds()[i])``: the anchored seed list of ATLAS-0 (ADR-0004 item 4,
stream "tuning"; also in experiments/exploratory/ATLAS-0/seeds.json), the same for every method, the prototype and the
model included. No other seed is used for tuning. A draw
visits the parameters in sorted name order, one generator per trial, so adding a method never changes the draws of
another and a trial can be redrawn from its seed alone.

Parameter names carry the stage they belong to (``pk_`` picker, ``as_`` associator, ``wf_`` waveform method), so a
pick-based pipeline is tuned jointly: its space is the union of its picker's and its associator's (card: "picker and
associator parameters in one draw").

Where the card's wording leaves the mapping to a program parameter open, the mapping is stated next to the parameter
(see the comparator README, ``experiments/exploratory/COMPARATORS/README.md``, "Card mappings").
"""
import math
from dataclasses import dataclass

import numpy as np

from ..atlas.seeds import N_TRIALS, tuning_seeds

SEEDS = tuple(tuning_seeds())


@dataclass(frozen=True)
class LogUniform:
    lo: float
    hi: float

    def draw(self, rng):
        return float(math.exp(rng.uniform(math.log(self.lo), math.log(self.hi))))

    def contains(self, v):
        return self.lo <= v <= self.hi


@dataclass(frozen=True)
class Uniform:
    lo: float
    hi: float

    def draw(self, rng):
        return float(rng.uniform(self.lo, self.hi))

    def contains(self, v):
        return self.lo <= v <= self.hi


@dataclass(frozen=True)
class IntRange:
    """integers lo, ..., hi (both included), uniform."""
    lo: int
    hi: int

    def draw(self, rng):
        return int(rng.integers(self.lo, self.hi + 1))

    def contains(self, v):
        return int(v) == v and self.lo <= v <= self.hi


@dataclass(frozen=True)
class Choice:
    options: tuple

    def draw(self, rng):
        return self.options[int(rng.integers(len(self.options)))]

    def contains(self, v):
        return v in self.options or (isinstance(v, list) and tuple(v) in self.options)


@dataclass(frozen=True)
class Pinned:
    """a parameter held at one value that still consumes the draw of the registered Choice it replaces
    (rng.integers(len(registered))), so every other parameter of a trial keeps the value of the registered draw."""
    value: object
    registered: tuple

    def draw(self, rng):
        rng.integers(len(self.registered))
        return self.value

    def contains(self, v):
        return v == self.value


# pickers (probability traces are computed once per weights; these act on the traces)
PHASENET = {"pk_peak_distance": LogUniform(0.05, 0.5), "pk_threshold": LogUniform(0.03, 0.5),
            "pk_weights": Choice(("original", "instance"))}
PHASENO = {"pk_peak_distance": LogUniform(0.05, 0.5), "pk_threshold": LogUniform(0.03, 0.5)}
EQTRANSFORMER = {"pk_peak_distance": LogUniform(0.05, 0.5), "pk_threshold": LogUniform(0.03, 0.5),
                 "pk_weights": Choice(("original", "instance"))}

# associators
# GaMMA as registered (registration 1); its ledger rows (tuning_ledger.csv) are redrawn with this space
GAMMA_REG1 = {"as_min_picks": IntRange(3, 8), "as_max_sigma11": Uniform(0.3, 2.0),
              "as_method": Choice(("BGMM", "GMM"))}
# ATLAS-0 amendment 3 (E63): GMM fails in GaMMA 80394dd before fitting (gamma/_base.py:140 reads covariance_prior,
# which only BayesianGaussianMixture has), so the method is pinned to BGMM. Pinned consumes the registered draw, so
# trial i keeps the registered as_min_picks, as_max_sigma11 and picker values; only the 37 GMM trials change method
GAMMA = {"as_min_picks": IntRange(3, 8), "as_max_sigma11": Uniform(0.3, 2.0),
         "as_method": Pinned("BGMM", ("BGMM", "GMM"))}
# ATLAS-0 amendment 1 (ADR-0009 item 1): REAL tint and PyOcto min_interevent_time log-uniform in [0.1, 5] s, so neither
# is kept from reporting two events inside regime T by construction
PYOCTO = {"as_min_picks": IntRange(3, 8), "as_tolerance": Uniform(0.3, 2.0),
          "as_pick_match_tolerance": Uniform(0.5, 2.0),        # card: "association cutoff in [0.5, 2] s"
          "as_min_interevent_time": LogUniform(0.1, 5.0)}
REAL = {"as_np0": IntRange(3, 6), "as_nps0": IntRange(5, 12), "as_tint": LogUniform(0.1, 5.0)}

# waveform methods
QUAKEMIGRATE = {"wf_threshold": Uniform(1.2, 3.0), "wf_stalta": Choice(((0.2, 1.0), (0.5, 5.0))),
                "wf_marginal_window": Uniform(0.3, 1.0)}
MALMI = {"wf_threshold": Uniform(0.1, 0.5), "wf_picker": Choice(("PhaseNet", "EQTransformer")),
         "wf_stack_window": Uniform(0.3, 1.0)}
EQNET = {"wf_threshold": Uniform(0.05, 0.7)}
EQCORRSCAN = {"wf_mad": Uniform(8.0, 12.0), "wf_min_stations": IntRange(3, 6)}

PICKERS = {"phasenet": PHASENET, "phaseno": PHASENO, "eqtransformer": EQTRANSFORMER}
ASSOCIATORS = {"gamma": GAMMA, "pyocto": PYOCTO, "real": REAL}
ASSOCIATORS_REG1 = {**ASSOCIATORS, "gamma": GAMMA_REG1}
WAVEFORM = {"quakemigrate": QUAKEMIGRATE, "malmi": MALMI, "eqnet": EQNET, "eqcorrscan": EQCORRSCAN}


def space_of(method, reg1=False):
    """the search space of a method: "picker+associator" for a pick-based pipeline (the union of both spaces), or
    the name of a waveform method. reg1: the associator spaces of registration 1 (GaMMA with method {BGMM, GMM}),
    which drew the rows of the registered ledger."""
    if "+" in method:
        pk, asc = method.split("+", 1)
        sp = {**PICKERS[pk], **(ASSOCIATORS_REG1 if reg1 else ASSOCIATORS)[asc]}
    else:
        sp = dict(WAVEFORM[method])
    return sp


def draw(space, seed):
    """one configuration: every parameter of space drawn, in sorted name order, from one generator seeded with seed."""
    rng = np.random.default_rng(int(seed))
    return {k: space[k].draw(rng) for k in sorted(space)}


def trials(method, n=N_TRIALS, reg1=False):
    """the n configurations of a method's budget: [(trial, seed, config)]."""
    sp = space_of(method, reg1)
    return [(i, SEEDS[i], draw(sp, SEEDS[i])) for i in range(n)]


def check(space, config):
    """raise unless config holds exactly the parameters of space, each inside its range."""
    if set(config) != set(space):
        raise ValueError(f"parameters {sorted(config)} differ from the space {sorted(space)}")
    bad = [k for k, v in config.items() if not space[k].contains(v)]
    if bad:
        raise ValueError(f"outside the space: {bad}")
