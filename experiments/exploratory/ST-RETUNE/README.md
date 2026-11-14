# ST-RETUNE: retuning the comparators with independent station terms (T4.25, design before the run)

Status: design and pilot only (14 Nov 2026). **No retuning has been run and no ledger is touched.** The retuning needs
an amendment of ATLAS-0 anchored by the owner first (draft below). Owner's approval of the proposal: 14 Nov 2026,
as a hypothesis to test, with a pivot criterion. Depends on PR #47 (station terms, `T420-CEILING`).

## Why

T420-CEILING (PR #47) measured, on the tuning partition, that terms estimated from the injected sources do not raise
the comparators with PhaseNet picks at their frozen settings (PyOcto T 0.206 at best 0.202; GaMMA falls) while they
do with oracle picks (PyOcto T 0.17 to 0.30-0.33). Two things were left open: the terms came from the targets' own
sources (7 and 12 per network), and the frozen settings were tuned with the misfit present. A fair bar needs terms
from an independent set and a retuning with them.

## Question and hypothesis

H: with static station terms from independent catalogue events applied to the picks, a retuning of each comparator
(60 trials, the same spaces and seeds) reaches an objective above its frozen one. Stated as a hypothesis: the pilot
below does not test it. What it predicts: PyOcto objective above 0.2589, GaMMA a4 above 0.2356. Either outcome is
recorded; the comparison with our model uses the better configuration of each comparator (with or without terms), so
a negative result cannot be used in our favour.

## Pilot (run, cloud, 4 minutes): are independent terms available and sound?

Calibration set: SCSN analyst picks of 9 July 2019 (day 190) from the SCEDC public bucket
(`tools/fetch_scedc_phases.py OUT 190 190`; 1893 phase files, 720 events of M >= 1.5 in the Ridgecrest box, 26 729
picks; SHA-256 of `events.csv` fa384bee..., `picks.csv` 14db43df...). The EGF sources of the Ridgecrest sets are
recorded events of 4 July 2019 (the R0 hour, `data/cards/synth_v2.md`), so no calibration event is a target or the
source of a target. Selection fixed before looking at any association result: pick weight >= 0.5, epicentral distance
<= 60 km, |residual| < 3 s, at least 30 picks per term, median; missing terms filled (`station_terms.fill_missing`).
Residuals against the frame's model (Hadley-Kanamori). Script `catalogue_terms.py`, output
`results/st_retune/pilot_ridgecrest_day190.json`.

| | value |
|:--|--:|
| stations of the frame with a P and an S term | 11 of 11 (183 to 600 picks per term) |
| correlation with the in-sample terms of the injected sources (19 terms) | 0.913 |
| median difference, catalogue minus injected | +0.025 s |
| median abs. residual of the injected arrivals: no terms / catalogue / in-sample | 0.172 / 0.093 / 0.062 s |
| same after removing each event's median (free origin time) | 0.147 / 0.068 / 0.053 s |

| Frozen PyOcto, pooled (Ridgecrest corrected, Hi-net unchanged) | T | C | Ridgecrest T | Ridgecrest C |
|:--|--:|--:|--:|--:|
| PhaseNet picks, no terms (frozen objective) | 0.206 | 0.311 | 0.255 | 0.460 |
| PhaseNet picks, catalogue terms | 0.202 | 0.328 | 0.248 | 0.466 |
| oracle picks, no terms | 0.170 | 0.538 | 0.139 | 0.672 |
| oracle picks, catalogue terms | 0.309 | 0.590 | 0.303 | 0.704 |

Reading (measured): independent terms exist for Ridgecrest, agree with the injected ones although the pick
definitions differ (analyst picks against AIC onsets), and halve the residual of arrivals they never saw. Unlike the
terms of PR #47 that left some stations uncorrected, they do not lower the comparator with PhaseNet picks (T 0.202
against 0.206, C 0.328 against 0.311: no change beyond a few targets) and they give the same gain with oracle picks.
The pilot used one day; the registered estimate uses days 190-196 (the M2-GP set, 3095 events).

**Hi-net is not covered.** The cloud does not reach the JMA unified catalogue or NIED; no independent picks were
found. PR #47 also showed that static terms do not explain the Hi-net misfit (in-sample they reach isolated recall
0.65, needed 0.85). Options, for the owner: (i) terms only where an independent set exists (Ridgecrest now, Alto
Tiberina to be checked), Hi-net tuned without terms and flagged; (ii) the owner downloads JMA arrival data with a
NIED account. The draft below takes (i).

## Design of the retuning (to run only after the amendment is anchored)

1. Terms: `station_terms.from_catalogue` on days 190-196 with the selection above, saved as a JSON file whose SHA-256
   enters the amendment. Typed as onset terms (`ttmt.station_terms.OnsetDelay`) when read by our model (E54).
2. Runner: `ttmt.associate.tune` applies `station_terms.apply` to the picks of every record of a network that has a
   term file, before association, for every method. One code path for all comparators and for our model.
3. Ledger: a new file `tuning_ledger_a6.csv`; the ledger module refuses any other file when terms are on, as it
   does for GaMMA after amendment 4. The registered ledgers and frozen files are not changed.
4. Trials: 60 per method from `seeds.tuning_seeds(60)`, the spaces of amendments 1 and 4 unchanged, the tuning
   partition b60d99b3 and the objective of amendment 1. A smoke run of 3 trials per method first (rule of E63).
5. Cost: PyOcto about 10 s per trial (10 min); GaMMA 3 to 6 min per trial (3 to 6 h, one process per core, BLAS
   pinned). CPU only, cloud or the owner's machine; no GPU (PhaseNet probabilities are cached).
6. HARPA joins when amendment 5 is anchored (PR #51), with the same terms.

## Reading rule (fixed now)

- For each method the configuration used against our model is the better of {frozen without terms, retuned with
  terms} by the tuning objective. The read partition is not touched by this experiment.
- H holds for a method if the retuned objective exceeds the frozen one. With one search per method there is no
  interval; a difference below 0.01 (about 3 targets) is reported as "no change".

## Risks and pivot criterion

| Risk | Sign | Action |
|:--|:--|:--|
| The gain needs a better picker, not better times | retuned objective within 0.01 of the frozen one for PyOcto | Record as negative; the fair bar stays the frozen configurations; stop spending on terms for the comparators and keep them only for our model (M6 needs them, X7) |
| Terms tuned on the same days as something else | overlap of calibration days with any partition | Checked: days 190-196 are outside the R0 hour; recheck against the atlas manifests before anchoring |
| Analyst picks and picker onsets differ by a phase-dependent lag (E54) | median difference above 0.1 s at registration | Estimate the lag on the calibration events with PhaseNet picks (needs waveforms: owner's machine) |
| Hi-net stays without terms | no independent catalogue | Option (i) above, flagged in every table; our model gets no terms there either |

Pivot: if neither PyOcto nor GaMMA gains at least 0.01 in the objective, the line "station terms for the
comparators" closes as a negative result after this one run; no second search space and no further estimates.

## Draft of amendment 6 of ATLAS-0 (not anchored; number 6 assumes amendment 5 is anchored first)

Task ATLAS-0, amendment 6 (date), written after the exploratory runs of T420-CEILING and the pilot above, and before
any tuning trial with station terms. The registered text, amendments 1 to 5 and every ledger are kept unchanged.

1. Station terms as a fixed preprocessing step. For a network with an independent calibration set, static
   per-station, per-phase terms are subtracted from the picks before association, for every method of the comparator
   set and for the model. Ridgecrest: SCSN analyst picks of 9 to 15 July 2019 (days 190-196), selection of
   `tools/fetch_scedc_phases.py`, estimator `ttmt.associate.station_terms.from_catalogue` (weight >= 0.5, distance
   <= 60 km, |residual| < 3 s, at least 30 picks, median, missing terms filled); file SHA-256 (to fill). Hi-net: no
   independent set, no terms, flagged. Alto Tiberina: (to fill after the check).
2. Ledger and budget. 60 trials per method from `seeds.tuning_seeds(60)` over the registered spaces, on the tuning
   partition and objective of amendment 1, in `tuning_ledger_a6.csv`; the best trial is frozen next to the existing
   frozen files, which stay.
3. Use. Against the model, each comparator runs with the better of its two frozen configurations by the tuning
   objective. This can only raise the bar of ADR-0008; it is not chosen on any read-partition result.

Disclosure: the frozen objectives (PyOcto 0.2589, GaMMA 0.2356), the exploratory results of T420-CEILING (terms from
the injected sources, three tolerance points) and the pilot above were known when this was written. Those runs used
the tuning partition only.
