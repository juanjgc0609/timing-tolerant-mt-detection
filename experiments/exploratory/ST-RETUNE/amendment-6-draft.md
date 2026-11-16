## Amendment 6 (station terms; draft, to be anchored as `reg/ATLAS-0/amendment-6`)

Draft revised on 14 November 2026 after the independent review of PR #53. It is appended to
`experiments/exploratory/ATLAS-0/card.md` and anchored by the owner only after amendment 5 is anchored; the date
below is the date of anchoring.

Owner approval. On 14 November 2026 the owner approved retuning the comparators with station terms under an
amendment ("apruebo todo lo demás"). On 15 November 2026 (09:17 Bogotá) the owner approved the pending steps,
which the project coordinator reads as approval of this amendment: "verifica y apruebo todo lo que haya que
hacer". The anchoring and the tag remain the owner's own actions.

Task ATLAS-0, amendment 6 (⟨date of anchoring⟩), written after the exploratory runs of T420-CEILING (PR #47) and of
ST-RETUNE (PR #53), and before any tuning trial with station terms. The registered text, amendments 1 to 5, every
ledger and every frozen configuration are kept unchanged; these items add to them.

1. Station terms as a fixed preprocessing step. For a network with an independent calibration set, static
   per-station, per-phase terms are subtracted from the picks before association
   (`ttmt.associate.pipeline.run_pipeline`, argument `terms`), the same file for every method of the comparator set
   and for the model.
   - Ridgecrest: SCSN analyst picks of 9 to 15 July 2019 (days 190-196) from the SCEDC public bucket, selection of
     `tools/fetch_scedc_phases.py` (type eq, M >= 1.5, box 35.45-36.15 N, 117.95-117.25 W): 3133 events, 135 102
     picks, SHA-256 of `events.csv` aa87e042b84b47ad8a998b09816401ae92dba57966a051cf895d97bc20a061ed and of
     `picks.csv` 90a2146a986958ba4e341436d398411220d9037272b8fbfe48ee323e5ad860ff. Estimator
     `ttmt.associate.station_terms.from_catalogue` against the frame's model (Hadley-Kanamori): pick weight >= 0.5,
     epicentral distance <= 60 km, |residual| < 3 s, at least 30 picks per term, median; a missing term takes the
     median of its phase. Terms file `experiments/exploratory/ST-RETUNE/terms_ridgecrest_2019.json`, SHA-256
     1ba1719e07a6193ca490009a7d488efba8469a881d7007bdd91a50360ab898bc (22 terms, 677 to 2627 picks each).
   - Hi-net: no independent calibration set is available; no terms for any method, the model included.
   - Alto Tiberina: no terms for any method, the model included (no calibration set was sought; it is not part of
     the tuning partition).
   A network gains terms later only by a further amendment, before any run on it.
2. Ledger and budget. 60 trials per method (PyOcto, GaMMA; HARPA under amendment 5) from `seeds.tuning_seeds(60)`
   over the registered spaces (amendments 1, 4 and 5), on the tuning partition of amendment 1, item 4, with the
   objective of amendment 1, item 2. Every trial goes to `tuning_ledger_a6.csv`; the ledger module refuses any other
   ledger for a trial with terms and refuses a trial without terms in that one (`ledger.require_ledger`). A smoke run
   of the first 3 trials of each method comes first; its rows stay in the ledger and count in the 60. The best trial
   of each method is frozen as `frozen/<method>_a6.json`, next to the existing frozen files, which stay.
   HARPA's trials with terms also go to `tuning_ledger_a6.csv`, and none may run before amendment 5 is anchored and
   its 60 trials without terms are in `tuning_ledger_a5.csv` (`ledger.require_ledger` refuses them until then).
3. Use against the model. Each comparator runs with whichever of its two frozen configurations (without terms, and
   with terms under this amendment) has the higher tuning objective; an exact tie keeps the configuration without
   terms. The choice for each comparator is fixed and recorded in the card before any read-partition run. This can
   only raise the bar of ADR-0008. The model has a single configuration: terms on wherever an independent set exists
   (item 1), off elsewhere. It is not chosen among options.
4. Reporting. Results are reported per network. On Hi-net no method is corrected, the misfit of the comparators'
   1-D model is large there (PyOcto with oracle picks reaches 0.50 isolated recall), and the model is built to
   tolerate timing misfit, so leaving Hi-net uncorrected favours the model structurally: an advantage on Hi-net alone
   is not claimed as a gain without this caveat.
5. One search. This retuning is run once. If neither PyOcto nor GaMMA gains at least 0.01 in the tuning objective
   over its frozen value, the result is recorded as negative and the line closes: no second search space, no other
   estimate of the terms and no further retuning with terms for the comparators.

Disclosure. Known when this was written: the tuned objectives of PyOcto (0.2589) and GaMMA amendment 4 (0.2356).
The exploratory results of T420-CEILING on the tuning partition: terms estimated from the injected sources failed
their acceptance criterion (Hi-net isolated recall with oracle picks 0.54, needed 0.85), and with PhaseNet picks no
estimate or tolerance tried raised PyOcto above its frozen recall in T. The ST-RETUNE checks at the frozen PyOcto
configuration: with the terms of 9 July alone, T 0.202 and C 0.328 against 0.206 and 0.311; with the terms of this
amendment, T 0.215 and C 0.348; so a result of "no change" for PyOcto is a likely outcome. Days 190-196 were the
development days of M2-GP (PR #46) and days 197-203 its confirmation days; that use was for the model's own station
terms and involved no comparator. The M2-GP fetch of 9 November held 3095 events and 133 952 picks; the bucket has
since gained 38 events, and the hashes of item 1 fix the set used here. The three selection values (weight >= 0.5,
distance <= 60 km, at least 30 picks) were set when the one-day check was written, before any association result
with catalogue terms existed; no other values were tried, and the repository cannot prove the order because they
entered in the same commit as that check. The S term of TOW2 (-0.33 s, 677 picks) is the median of two groups of
residuals (about -0.54 s, 439 picks, and +0.77 s, 238 picks), so it is not a static delay; it is kept because the
rule was fixed before it was seen, and it is 1 of 22 terms.

Freeze: `ttmt.associate.station_terms`, `ttmt.associate.pipeline`, `ttmt.associate.ledger`, `ttmt.associate.tune`,
`experiments/exploratory/COMPARATORS/tune_cli.py`, the terms file of item 1 and their tests
(`tests/test_associate.py`, `tests/test_detectability.py`); environment as amendment 4.
