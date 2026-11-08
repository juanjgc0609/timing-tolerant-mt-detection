# Stage-loss diagnostic of the frozen comparators (tuning partition)

Question: where do the frozen pick-based comparators lose their targets, in the picker or in the associator, and how
far are they from the information limit of Corollary A.7? The answer decides which lever of the model pays most.

Data: the tuning partition b60d99b3d0bbbbfa only (Ridgecrest and Hi-net EGF sets, 50 records, n_T = 223,
n_C = 305), the same truth and PhaseNet `original` probability files the registered tuning read (copied read-only
from the owner's machine; nothing was recomputed). No read-partition data, no ledger was touched. Objective,
budget (30 false detections per hour) and hit rule (LOC+T-V) are those of `ttmt.associate.tune.objective`.

```
python experiments/exploratory/STAGE-LOSS/stage_loss.py experiments/exploratory/COMPARATORS/frozen/phasenet+pyocto.json <tune>/ss_ridgecrest_2019_egf_v1/tune_phasenet.json <tune>/ss_hinet_2012_egf_v1/tune_phasenet.json --out results/stage_loss/phasenet+pyocto.json
```

Cost: one associator run per realisation (the k = 4, 6, 8 records of a realisation share picks and detections), about
15 s per condition for PyOcto and 30-400 s for GaMMA on one cloud core. Each condition is one run of the frozen
configuration, not a search.

## Conditions

| Code | Picks | Associator |
|:--|:--|:--|
| PN | PhaseNet, frozen threshold and peak distance | frozen |
| PN-clean | PN minus every pick that matches no true arrival (0.5 s, nearest arrival, one pick per arrival) | frozen |
| OR | every finite injected arrival of every catalogued event (perfect picker) | frozen |
| OR-snr1, OR-snr2 | OR restricted to arrivals with per-station SNR >= 1, >= 2 | frozen |
| OR-mit | OR | PyOcto with min_interevent_time 0.1 s |
| OR-tight | OR | PyOcto with min_interevent_time, tolerance and pick match tolerance 0.1 s |

A missed target is attributed to the first stage that fails: picker (fewer than min_picks = 6 of its arrivals picked),
association (no detection within LOC+T), resolution (a detection within LOC+T, but not eligible under the dominance
rule or taken by the partner), budget (its eligible detection is below the operating point).

## Results

PN reproduces both frozen objectives exactly (PyOcto 0.2063 / 0.3115, GaMMA registered 0.1659 / 0.3016).

| Condition | PyOcto T | PyOcto C | GaMMA (registered, provisional) T | GaMMA C |
|:--|--:|--:|--:|--:|
| PN | 0.206 | 0.311 | 0.166 | 0.302 |
| PN-clean | 0.166 | 0.423 | 0.211 | 0.426 |
| OR (perfect picks) | 0.170 | 0.538 | 0.143 | 0.675 |
| OR-snr1 | 0.233 | 0.502 | 0.152 | 0.528 |
| OR-snr2 | 0.229 | 0.502 | 0.202 | 0.466 |
| OR-mit | 0.161 | 0.492 | | |
| OR-tight | 0.220 | 0.498 | | |
| A.7 limit (S2, conservative) | 0.93 | | 0.93 | |

Stages of the T targets (of 223):

| | hit | picker | association | resolution | budget |
|:--|--:|--:|--:|--:|--:|
| PyOcto PN | 46 | 120 | 32 | 25 | 0 |
| PyOcto OR | 38 | 0 | 39 | 139 | 7 |
| GaMMA PN | 37 | 112 | 27 | 47 | 0 |
| GaMMA OR | 32 | 0 | 40 | 151 | 0 |

Of the 146 T targets PyOcto leaves unresolved with perfect picks, 139 have a single detection for the whole pair
(merged), 7 have two detections that are too imprecise for the dominance rule. With PhaseNet picks 59 are merged.

Picker: PhaseNet picks 4 % of the arrivals at station SNR < 0.5, 19 % at 1-2, 53 % at 4-8 and 65 % above 8
(`pick_recall_by_snr` in the result files). A.7 (scenario S2, same mechanism, total error 0.1,
dtau_min = 0.345 / sqrt(SNR_tot), SNR_tot the root sum of squares of the subset's station SNR; K_mec is nan in the
EGF sets, so the conservative scenario): median SNR_tot of the T targets 16.4, median dtau_min 0.085 s, T separation
quartiles 0.30 / 0.46 / 0.67 s; 93 % of the T targets are above the limit.

## Reading

1. **T is lost in the associator, not in the picker.** A perfect picker does not raise T recall for either program
   (PyOcto 0.21 -> 0.17, GaMMA 0.17 -> 0.14): both emit one event for the two members of a pair. Tightening PyOcto's
   time tolerances to 0.1 s with perfect picks gives 0.22, so the merge is structural (one event per pick cluster), not
   a tuning artefact. The gap to the information limit in T is about 70 points (0.2 -> 0.93), and it is all
   association and resolution. This is the lever: an association that can explain one pick cluster with two
   overlapping events (the model's joint waveform likelihood) is where a method can beat both programs by a wide
   margin.
2. **C is lost mostly in the picker.** Perfect picks raise C from 0.31 to 0.54 (PyOcto) and from 0.30 to 0.68 (GaMMA);
   removing only the false picks already gives 0.42-0.43. Low-SNR arrivals are the bulk: about half of the C targets
   fail with fewer than 6 picked arrivals. A better detector at station SNR 1-4 (waveform-domain stacking rather than
   per-station picks) is the second lever, and the one that protects H2.
3. **PyOcto's min_interevent_time is not the cause** (OR-mit and PN-mit change nothing in T).
4. GaMMA here is the registered configuration (provisional, E63); the a4 configuration replaces it when it is frozen.

## Isolated events (oracle check)

Isolated targets (separation > 3 s, n = 129, part of C) test whether the oracle reaches 1 when no pair can merge.

| Condition | PyOcto iso | GaMMA (registered) iso |
|:--|--:|--:|
| PN | 0.450 | |
| OR (perfect picks) | 0.767 (Ridgecrest 0.96, Hi-net 0.50) | 0.930 |

With perfect picks GaMMA is near 1 (9 of 129 missed, all resolution), and so is PyOcto on Ridgecrest. PyOcto's 30
misses are Hi-net events: none is outside the search box, every nearest detection is within 2 km in epicentre (median
1.5 km), and half are more than 0.5 s off in origin time; on Hi-net its origin times are early by 0.31 s (median,
quartiles -0.40 / -0.27 s). The injected arrivals do not follow the associators' 1-D model there (residuals of
arr_t_obs against JMA2001, median P -0.15 s and S -0.35 s, interquartile ranges 0.39 and 0.56 s; Ridgecrest against
Hadley-Kanamori +0.11 and +0.15 s). So the isolated-event limit is travel-time mismatch with the 1-D model as PyOcto's
locator absorbs it, not the picks: GaMMA, with the same model and picks, reaches 0.93. For our model this argues for
station terms (or a model-error term in the likelihood) before any claim on C. Script: `stage_loss.py` reports
`recall_iso` and `stages_iso`; GaMMA's OR row is `results/stage_loss/phasenet+gamma_registered_provisional_OR_iso.json`.

## Limits

Tuning partition only (50 records, 20 realisations, one development objective, no intervals). The oracle uses the
injected arrival times, so it is an upper bound for a pick-based method, not for a waveform method. Attribution is
first-failure, so a target counted as "picker" may also have failed later.
