# T4.20 detectability ceiling in C, T4.22 decision rule, and Hi-net station terms (T4.24 / T4.25)

Status: card written and committed **before** the run (8-9 Nov 2026, commit de99a63). Run on 14 Nov 2026; the
sections above "Results" are unchanged since then except for this paragraph and the motivation note below.
Exploratory, tuning partition only (b60d99b3d0bbbbfa, the STAGE-LOSS inputs). No read-partition data, no ledger, no
anchored card is touched.

Motivation added on 14 Nov (supervision check X7): with oracle picks placed at the associator's own model times
PyOcto separates pairs (T 0.60-0.67), so its merging of pairs is travel-time misfit and not structural. A fair
comparison therefore gives the same station corrections to every comparator (PyOcto, GaMMA, and HARPA once
amendment 5 is anchored; HARPA is not run here) and to our model.

## Questions

1. **T4.20.** STAGE-LOSS (PR #38) found that C is lost mostly in the picker: about half of the C targets have fewer
   than 6 picked arrivals and station SNR 1-4. How much of that loss could any network detector recover at the
   objective's budget (30 false detections per hour)? This is the regime-C analogue of Corollary A.7.
2. **T4.24 / T4.25.** PyOcto with oracle picks reaches 0.96 isolated recall on Ridgecrest and 0.50 on Hi-net; the
   injected Hi-net arrivals deviate from JMA2001 (median residual P -0.15 s, S -0.35 s). Do static station terms,
   given to the comparator, close that gap? Any claim in C must be measured against the corrected comparator.
3. **T4.22.** Should C become a secondary superiority hypothesis (by amendment, before reading atlas data)?

## Method

**Ceiling** (`src/ttmt/detectability.py`). Ideal network subspace detector, Gaussian noise independent between
stations, propagation known (clairvoyant Green's functions), source unknown in a p-dimensional subspace. Statistic
chi^2_p under H0, noncentral chi^2_p(lambda) under H1. Three detectors: p = 1 (coherent oracle), p = 6 (moment-tensor
GLRT, our model's class), p = K (one free amplitude per station). Signal energy lambda = kappa^2 * sum of squared card
SNR over the target's subset. Bandpass-then-sample is one linear filter, so card SNR^2 <= matched-filter SNR^2 at
each station: kappa = 1 is a lower bound of the energy (P on the vertical only; S and horizontals ignored); kappa = 2
is a sensitivity value only. Threshold by the union bound: alpha = 30 / (3600 * 100 Hz * n_loc) per trial, n_loc in
{1e2, 1e4, 1e6}. The ceiling ignores the hit rule's location test (LOC+T), so it is optimistic there and conservative
in the energy; it is a ceiling for the detection stage, not a forecast of our recall.

**Station terms** (`src/ttmt/associate/station_terms.py`). Per-station, per-phase median residual of the injected
arrivals (arr_t_obs) relative to catalogue origin and the frame's 1-D model, one median per source first (sources
identified by catalogue position). Cross-fitted by realisation: the terms of a realisation come only from sources that
do not appear in it, so no target informs its own correction. Applied to the picks (t - delta) before association,
because PyOcto 0.2.0 and GaMMA take no station terms; the same correction goes to every comparator and to our model.
Variant `-origin` removes each source's median over stations (a constant trades off with origin time).

Conditions: PN (reproduces the frozen objective), OR, PN+ST, OR+ST, OR+ST-origin; recall per network, isolated recall
included.

## Decision rule for T4.22 (fixed before the run)

Primary cell: p = 6, kappa = 1, n_loc = 1e4. R_loss = mean detection probability of the ideal detector over the C
targets the comparator misses (the recoverable fraction of the loss).

Recommend C as a secondary superiority hypothesis only if all three hold:

- **(a)** R_loss >= 0.30 against PN, and >= 0.20 in each network separately;
- **(b)** at most 10 % of the comparator's own C hits have ideal detection probability below 0.5 (otherwise the
  kappa = 1 energy bound is too loose to read, and the decision is **deferred** to a waveform-based ceiling with
  empirical noise thresholds, which needs the noise blocks on the owner's machine);
- **(c)** (a) also holds against PN+ST (the comparator with the same station terms; equity, T4.25).

If (a) or (c) fails, C stays a non-inferiority hypothesis (H2) and the result is recorded as negative. If all pass,
the statistical design thread computes the power of the clustered test for half of the recoverable gain before any
amendment, and the owner anchors the amendment; nothing is read from the atlas before that.

Station terms are accepted for use in C (by us and by the comparators) if OR+ST raises Hi-net isolated recall to
>= 0.85 without lowering Ridgecrest isolated recall by more than 0.02. If not, the misfit is not static (distance or
depth dependent) and the next step is a model-error term (tolerance as a function of distance), not more terms.

## Run (owner's machine, CPU only, about 2-3 min for PyOcto)

```
python experiments/exploratory/T420-CEILING/ceiling_c.py experiments/exploratory/COMPARATORS/frozen/phasenet+pyocto.json <tune>/ss_ridgecrest_2019_egf_v1/tune_phasenet.json <tune>/ss_hinet_2012_egf_v1/tune_phasenet.json --out results/t420_ceiling/phasenet+pyocto.json
```

`<tune>` is `C:\Users\juanj\ttmt\tune`. GaMMA (frozen `phasenet+gamma_a4.json`) is the same command, slower (30-400 s
per condition). Smoke test on a synthetic two-network scene: the code runs end to end, the terms recover injected
delays (in-sample residual 0.01 s from 0.19 / 0.42 s), unit tests in `tests/test_detectability.py`.

## Results (14 Nov 2026)

Run in the cloud, CPU, on a read-only copy of the tuning partition (62 files, 77 MB). PyOcto 0.2.0 with pyrocko
2026.11.05; GaMMA at 80394dd with numba 0.68.1; numpy 2.5.3, scipy 1.18.1, scikit-learn 1.8.0, Python 3.13. Both
baselines reproduce the frozen objective exactly (PyOcto T 0.2063, C 0.3115; GaMMA a4 T 0.1794, C 0.2918), so the
version differences from the frozen environment do not move the numbers. Wall time: PyOcto 63 s for the whole script,
GaMMA about 8 min. Raw output: `results/t420_ceiling/phasenet+pyocto.json`, `phasenet+gamma_a4.json`,
`terms_variants_pyocto.txt`, `terms_variants_gamma_a4.txt`.

### 1. Ceiling in C (T4.20)

305 C targets (116 Hi-net, 189 Ridgecrest). Network SNR (root of the summed squared card SNR over the target's
subset): quartiles 10.6, 19.0, 38.6. Stations per target: 2 (28), 3 (20), 4 (115), 5 (24), 6 (73), 7 (8), 8 (37).
Primary cell: p = 6, kappa = 1, n_loc = 1e4 (the detector needs network SNR 6.61 for a detection probability of 0.5).

| Mean detection probability of the ideal detector | p = 1 | **p = 6** | p = K |
|:--|--:|--:|--:|
| all C targets (ceiling recall) | 0.857 | **0.828** | 0.837 |
| C targets PyOcto misses (R_loss, 210) | 0.792 | **0.750** | 0.763 |
| of those, lost at the picker (148) | 0.729 | **0.680** | 0.698 |
| Hi-net misses (108) | 0.811 | **0.769** | 0.792 |
| Ridgecrest misses (102) | 0.772 | **0.729** | 0.733 |
| C targets GaMMA a4 misses (216) | 0.798 | **0.757** | 0.770 |

Sensitivity of R_loss at p = 6: n_loc 1e2 0.790, 1e6 0.717; kappa = 2 0.897. No comparator hit has an ideal
detection probability below 0.5 (0 % in every cell, both comparators), so the kappa = 1 energy bound is readable.

### 2. Decision rule of T4.22

| Condition | Value | Threshold | Result |
|:--|--:|--:|:--|
| (a) R_loss against PN, PyOcto | 0.750 | >= 0.30 | passes |
| (a) per network (Hi-net, Ridgecrest) | 0.769, 0.729 | >= 0.20 | passes |
| (b) comparator hits with Pd < 0.5 | 0 % | <= 10 % | passes |
| (c) R_loss against PN+ST (registered terms) | 0.774 (0.785, 0.762) | as (a) | passes |

The same holds against GaMMA a4 (0.757; against PN+ST 0.762). Condition (c) holds for every estimate of the terms
tried in section 4: the best C recall of a comparator with PhaseNet picks is 0.354, against a ceiling of 0.83.

**By the rule fixed before the run, C is recommended as a secondary superiority hypothesis.** What this does and does
not say: about three quarters of the C targets the comparators miss carry enough energy for an ideal network detector
at 30 false detections per hour. It is a ceiling of the detection stage with known Green's functions and Gaussian
noise independent between stations; it ignores the location test of the hit rule and real noise is not Gaussian
(G-M5c), so it is not a forecast of our recall. The next steps of the rule stand: the power of the clustered test
for half of the recoverable gain is computed first, and the owner anchors the amendment; nothing is read from the
atlas before that.

### 3. Station terms: the acceptance criterion fails

| Recall, frozen settings | T | C | iso | Ridgecrest iso | Hi-net iso |
|:--|--:|--:|--:|--:|--:|
| PyOcto PN (frozen objective) | 0.206 | 0.311 | 0.450 | 0.680 | 0.130 |
| PyOcto PN+ST | 0.126 | 0.239 | 0.302 | 0.520 | 0.000 |
| PyOcto OR | 0.170 | 0.538 | 0.767 | 0.960 | 0.500 |
| PyOcto OR+ST | 0.296 | 0.620 | 0.783 | 0.960 | **0.537** |
| PyOcto OR+ST-origin | 0.314 | 0.587 | 0.806 | 1.000 | 0.537 |
| GaMMA a4 PN (frozen objective) | 0.179 | 0.292 | 0.395 | 0.680 | 0.000 |
| GaMMA a4 PN+ST | 0.161 | 0.275 | 0.372 | 0.640 | 0.000 |
| GaMMA a4 OR | 0.143 | 0.610 | 0.822 | 1.000 | 0.574 |
| GaMMA a4 OR+ST | 0.184 | 0.580 | 0.822 | 1.000 | 0.574 |

Criterion: OR+ST raises Hi-net isolated recall to >= 0.85. Measured 0.537 (PyOcto) and 0.574 (GaMMA). **Not
accepted.** Two causes, both measured (`terms_variants.py`, first line of its output):

1. **The registered cross-fitting gives Hi-net no terms at all.** Hi-net has 7 distinct sources and every realisation
   holds 5 to 7 of them, so after leaving out the sources of the realisation fewer than 3 remain and all 14 terms are
   0 in all 10 realisations (median |residual| 0.338 s before and after). The Hi-net rows above change only through
   the pooled operating point. This is a design error of this card: the number of distinct sources per network was
   not counted before fixing the rule. Ridgecrest (12 sources) gets 17 to 19 of 22 terms and its out-of-sample
   median |residual| falls from 0.172 to 0.093 s.
2. **Even with terms estimated in-sample the misfit on Hi-net is not static.** With every source used (optimistic: a
   target informs its own correction) PyOcto OR reaches Hi-net iso 0.648, still below 0.85, and the in-sample RMS
   residual stays at 0.20 s (P) and 0.14 s (S) from 0.26 and 0.44 s.

By the card, the next step for Hi-net is a model-error term (source or distance dependent), not more static terms.

### 4. Estimates of the terms and tolerances (added after the run, not registered)

Added because of result 3. Estimates: cf (registered), cf-fill (missing terms take the median of the phase, so that
no relative shift appears between stations), all (in-sample) and loro (the other realisations, the scheme of X7);
the last two are optimistic. Recall pooled over both networks.

| PyOcto | PN: T | PN: C | OR: T | OR: C |
|:--|--:|--:|--:|--:|
| no terms | **0.206** | 0.311 | 0.170 | 0.538 |
| cf (registered) | 0.126 | 0.239 | 0.296 | 0.620 |
| cf-fill | 0.152 | 0.279 | 0.296 | 0.593 |
| all (optimistic) | 0.143 | 0.315 | 0.327 | 0.646 |
| loro (optimistic) | 0.143 | 0.354 | 0.327 | 0.636 |

| PyOcto, PhaseNet picks, T (C) | frozen tolerances (0.65 s, 0.92 s) | 0.3 s, 0.1 s | 0.1 s, 0.1 s |
|:--|--:|--:|--:|
| no terms | **0.206** (0.311) | 0.175 (0.226) | 0.126 (0.230) |
| cf-fill | 0.152 (0.279) | 0.090 (0.230) | 0.170 (0.226) |
| all (optimistic) | 0.143 (0.315) | 0.126 (0.315) | 0.202 (0.331) |

| GaMMA a4 | PN: T | PN: C | OR: T | OR: C |
|:--|--:|--:|--:|--:|
| no terms | **0.179** | 0.292 | 0.143 | 0.610 |
| cf-fill | 0.130 | 0.298 | 0.139 | 0.525 |
| all (optimistic) | 0.103 | 0.289 | 0.135 | 0.561 |

Reading, with its basis:

- **Measured:** with oracle picks the terms raise PyOcto in T (0.170 to 0.30-0.33) and in C (0.54 to 0.59-0.65),
  with every estimate. This agrees in sign with X7; the size is smaller than X7's per-network figures (0.35-0.38)
  because here the operating point is pooled over both networks, as the objective does.
- **Measured:** GaMMA a4 does not gain from the terms even with oracle picks (T 0.143 to 0.135-0.184, C falls),
  which agrees with X8: its merging comes from its own clustering.
- **Measured:** with PhaseNet picks no combination of estimate and tolerance tried reaches the frozen T of 0.206 (best
  0.202, optimistic terms, tolerances 0.1 s), and GaMMA's T falls from 0.179 to 0.10-0.13. In C the optimistic terms
  give PyOcto at most 0.354. So the prediction of the supervision document (Part 7: "T rises from 0.21 towards
  0.25-0.35 with PhaseNet picks") is **not supported** by the terms alone at these settings.
- **Measured:** the picker is the first loss in T as well: 120 of 223 T targets have fewer than the minimum number
  of picked arrivals (PyOcto PN), and terms cannot return a pick that was never made.
- **Hypothesis, not tested:** with the frozen tolerances the misfit splits some pairs by accident (the picks of the
  second event do not fit the first), and removing the misfit removes those accidental splits. Test: count, per T
  pair, the detections with and without terms at fixed picks.
- **Hypothesis, not tested:** a full retuning with the terms (picker threshold and peak distance included, 60 trials
  per method under an amendment of ATLAS-0) moves the comparators above the frozen values. The three tolerance
  points above are a sensitivity, not a tuning. This retuning is the test that settles the fair bar.

Limits: one tuning partition (223 T and 305 C targets, 20 realisations); the terms come from the injected sources
themselves, 7 and 12 per network, and not from an independent set of catalogue events; `fill_missing` was added after
seeing the registered result; E65 (`first_arrival` loses the head wave for a source at the top of a layer) is not
fixed on this branch, and the residuals use the same travel-time function, so a term can absorb part of that error
for the sources it affects; HARPA and EQNet2 are not run (amendment 5 not anchored).

### 5. What follows (proposals; the owner decides)

1. Terms from an independent calibration set: catalogue events of each region with catalogue picks, outside every
   realisation of the tuning and read partitions, at least 30 sources per network. Hypothesis to test: they reach the
   in-sample residual (0.05-0.10 s on Ridgecrest).
2. For Hi-net, a source-region or distance-dependent model-error term, since static terms do not reach 0.85.
3. Amendment of ATLAS-0 (number 6 if 5 is anchored first): the terms as a fixed preprocessing step for every method,
   and a retuning of PyOcto, GaMMA and HARPA with them. The frozen configurations stay as they are until then.
4. Our model reads the same terms (typed as `OnsetDelay`, `ttmt.station_terms`), never its own.
