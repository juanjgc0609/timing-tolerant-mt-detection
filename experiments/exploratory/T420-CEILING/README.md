# T4.20 detectability ceiling in C, T4.22 decision rule, and Hi-net station terms (T4.24 / T4.25)

Status: card written and committed **before** the run (8-9 Nov 2026). Exploratory, tuning partition only
(b60d99b3d0bbbbfa, the STAGE-LOSS inputs). No read-partition data, no ledger, no anchored card is touched.

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

## Results

Pending the run.
