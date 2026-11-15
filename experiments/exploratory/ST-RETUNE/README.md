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
5. Cost (estimate of 14 Nov; corrected with the ledger times in the last section): PyOcto about 10 s per trial (10 min); GaMMA 3 to 6 min per trial (3 to 6 h, one process per core, BLAS
   pinned). CPU only, cloud or the owner's machine; no GPU (PhaseNet probabilities are cached).
6. HARPA joins when amendment 5 is anchored (PR #51), with the same terms.

## Reading rule (fixed now)

- For each method the configuration used against our model is the one of {frozen without terms, retuned with
  terms} with the higher tuning objective (plain maximum; the 0.01 below is only for reading H and for the pivot). The read partition is not touched by this experiment.
- H holds for a method if the retuned objective exceeds the frozen one. With one search per method there is no
  interval; a difference below 0.01 (about 3 targets) is reported as "no change".

## Risks and pivot criterion

| Risk | Sign | Action |
|:--|:--|:--|
| The gain needs a better picker, not better times | retuned objective within 0.01 of the frozen one for PyOcto | Record as negative; the fair bar stays the frozen configurations; stop spending on terms for the comparators and keep them only for our model (M6 needs them, X7) |
| Terms tuned on the same days as something else | overlap of calibration days with any partition | Checked: days 190-196 are outside the R0 hour; recheck against the atlas manifests before anchoring |
| Analyst picks and picker onsets differ by a phase-dependent lag (E54) | median difference above 0.1 s at registration | Estimate the lag on the calibration events with PhaseNet picks (needs waveforms: owner's machine) |
| Hi-net stays without terms | no independent catalogue | Option (i) above, flagged in every table; our model gets no terms there either; leaving it uncorrected favours our model, so results go per network and a Hi-net-only advantage is not claimed without that caveat (amendment, item 4) |

Pivot: if neither PyOcto nor GaMMA gains at least 0.01 in the objective, the line "station terms for the
comparators" closes as a negative result after this one run; no second search space and no further estimates.

## Registered estimate and review (added 14 Nov 2026, after the independent review)

Terms from days 190-196 (3133 events, 135 102 picks): `terms_ridgecrest_2019.json`, 22 of 22 terms with 677 to 2627
picks. They differ from the one-day terms by at most 0.028 s; correlation with the injected in-sample terms 0.915,
median difference +0.017 s; median residual of the injected arrivals 0.172 to 0.086 s. Frozen PyOcto, pooled
(`results/st_retune/check_ridgecrest_days190_196.json`): PhaseNet picks T 0.215 and C 0.348 (0.206 and 0.311 without
terms); oracle picks T 0.296 and C 0.600 (0.170 and 0.538). These are checks at the frozen configuration, not tuning
trials.

TOW2, flagged by the review (P +0.39 s, S -0.33 s): the S residuals fall in two groups, about -0.54 s (439 picks,
median distance 11.6 km) and +0.77 s (238 picks, 16.6 km), so the median is not a static delay. The interquartile
range of the S residuals is 0.99 s at TOW2 and 0.67 s at CCC, against 0.27 to 0.45 s at the other nine stations. The
cause was not found (hypotheses, not tested: two source clusters with different paths, or S picks on a converted
phase). The term is kept as the rule gives it, because changing the rule after seeing the terms would add a choice;
this is declared in the amendment.

Selection values: weight >= 0.5, distance <= 60 km and at least 30 picks were set when the one-day check was written
and no other values were tried; |residual| < 3 s comes from the M2-GP card.

Runner: `run_pipeline(..., terms=)` subtracts the terms for every method; `tune_cli.py run --terms ridgecrest=<file>
--ledger .../tuning_ledger_a6.csv`; `ledger.require_ledger` refuses any other combination (tests in
`tests/test_associate.py`). No trial has been run.

## Amendment 6 of ATLAS-0

Draft in `amendment-6-draft.md` (not anchored). Steps to anchor, after amendment 5 is anchored and this PR is
reviewed: the text is appended to `experiments/exploratory/ATLAS-0/card.md` with the date, the PR is merged, and the
owner tags the merge commit:

    git checkout main && git pull
    git tag -a reg/ATLAS-0/amendment-6 -m "ATLAS-0 amendment 6: station terms as fixed preprocessing, retuning in tuning_ledger_a6.csv"
    git push origin reg/ATLAS-0/amendment-6

Before the tag, the text enters the card with `python append_amendment.py "<date>" --write` (it refuses while
amendment 5 is not in the card), and that diff is reviewed.

Then the retuning, with `run_retune.py` (it refuses any trial unless the tag exists and is an ancestor of HEAD):

    python experiments/exploratory/ST-RETUNE/run_retune.py check <tune>     # no trial, writes nothing
    python experiments/exploratory/ST-RETUNE/run_retune.py smoke <tune>     # first 3 trials of each method
    python experiments/exploratory/ST-RETUNE/run_retune.py full <tune>      # resumes to 60 per method
    python experiments/exploratory/ST-RETUNE/run_retune.py read             # freeze, reading rule, pivot

## What the run can and cannot show (added 15 Nov 2026, before the run)

Dry run, cloud, 15 Nov 2026 (`run_retune.py check`, no association run, nothing written): both methods load 50
records, 30 of them (Ridgecrest) with terms, partition b60d99b3d0bbbbfa, terms file hash as item 1, a6 ledger empty;
`smoke`, `full` and `read` refused because the tag does not exist. The 60 trials of each method are those of its
ledger without terms (same seeds and parameters; two parameters differ in the last digit, 1 part in 1e16, between the
Windows ledger and Linux).

Resolution of the objective, from the ledgers (facts):

- The tuning partition has 223 targets in T and 305 in C. At the frozen recalls the binomial standard error is 0.027
  for each recall and 0.019 for the objective, taking targets as independent; they are not (7 and 12 sources per
  network), so the true error is larger. The 0.01 of the reading rule and of the pivot is about 3 targets, below
  this error.
- The search itself is noisy: the best PyOcto trial (0.2589) is 0.027 above the second (0.2323); in random halves of
  its 60 trials the best of 30 falls more than 0.01 below 0.2589 half of the time. For GaMMA a4 the top five lie
  within 0.012 and that happens 11 % of the time.

Consequence: a difference of 0.01 between the two best objectives is not evidence that terms help or hurt. The
decisions of amendment 6 do not need that evidence: item 3 takes the maximum (it can only raise the bar) and item 5
closes the line when the gain is small. H is read as registered, and the statement "terms help the comparators"
is made only with the support of the paired comparison below.

Paired comparison (secondary, descriptive, fixed here before the run; it changes no rule of the amendment): the
same 60 parameter sets are scored with and without terms, so `run_retune.py read` reports, per method, the number
of trials whose objective rises, falls or stays, and the median difference. Two-sided sign test at 0.05 on the
trials that change: with 60 changing trials it needs 39 rises; power 0.56, 0.84 and 0.97 if terms raise 65, 70 and
75 % of the configurations. PyOcto has 60 trials with a non-zero objective without terms and GaMMA a4 has 40, so
GaMMA has fewer informative pairs and less power. The trials share the targets, so the test speaks about this
partition and not about other data.

Cost, measured in the ledgers on the owner's machine: PyOcto 2.6 h for 60 trials (median 39 s), GaMMA a4 3.4 h
(median 44 s). The earlier figure of 10 min for PyOcto was wrong. Cloud times are not measured; the smoke run gives
them.
