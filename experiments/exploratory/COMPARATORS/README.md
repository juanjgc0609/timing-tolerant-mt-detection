# Comparators of the failure atlas: wrappers, tuning and reproduction

Owner: thread "Comparadores y ajuste" (plan role F2). The atlas thread keeps the card
(`experiments/exploratory/ATLAS-0/card.md`), the hit rule and every analysis; this thread installs and wraps the
comparators, tunes them on the development tuning partition under the card's budget, and hands the atlas frozen
configurations with the ledger that produced them. It does not read atlas results.

Code: `src/ttmt/associate/` (tests in `tests/test_associate.py`). Scripts here: `tune_cli.py` (the search on the
device) and `repro_puente.py` (gate G3a, target 1).

## Method status

| Method (card) | Status |
|:--|:--|
| PhaseNet + GaMMA, PhaseNet + PyOcto, PhaseNet + REAL | ready: wrappers tested with the real programs |
| EQTransformer + GaMMA / PyOcto (extended set) | ready (same picker path as PhaseNet) |
| PhaseNO + GaMMA / PyOcto | pending: PhaseNO inference written to the probability-file layout (notebook code, GPU machine) |
| EQcorrscan | pending: wrapper (separate `ttmt-eqc` env) |
| QuakeMigrate, MALMI | pending: install and wrapper (not in the environments of PR #8) |
| EQNet | not runnable: no released weights for the 2022 model (`environments/DEVIATIONS.md`) |

`ttmt.associate.pipeline.STATUS` holds the same table; `tune_cli.py run` refuses a method that is not ready.

## The tuning protocol as implemented (card "Tuning", ADR-0008 item 9, CR-07)

- Budget: 60 configurations per method, trial i drawn with `numpy.random.default_rng(tuning_seeds()[i])`, the
  anchored seed list of ATLAS-0 (`ttmt.atlas.seeds.tuning_seeds`, ADR-0004 item 4; `seeds.json`), the same for every
  method. A draw visits the parameters in sorted name order with
  one generator per trial, so a trial is redrawn from its seed alone.
- Pick-based pipelines are tuned jointly: one draw covers the picker (`pk_*`) and the associator (`as_*`).
- Objective: development recall in T and in C (C includes the isolated events) at 30 false detections per hour under
  LOC+T-V, with the operating point set on the tuning records (`ttmt.atlas.froc`), events weighted as in the
  endpoint. The search maximises `(recall_T + recall_C) / 2`, ties to the lower trial (ATLAS-0 amendment 1,
  ADR-0009 item 2); both recalls go to the ledger. A tuning partition without targets in T or in C is refused before
  the first trial.
- Only records whose partition is `"tune"` are accepted. The partition comes from the generator's manifest (ADR-0009
  item 4: 300 s noise blocks, 30/70 by `ttmt.synth.realdata.tune_read`); the record's
  `partition` field is read as the generator wrote it (per-dataset fallback included), never recomputed;
  `tune.tuning_split` only calls `tune_read`.
- Ledger: `experiments/exploratory/ATLAS-0/tuning_ledger.csv`, append-only (a repeated (method, trial) is refused),
  one row per trial including failed ones (status `error`, objective nan; a failed trial still spends its place).
  Columns: method, trial, seed, params, objective, recall_T, recall_C, n_T, n_C, hours, status, error, wall_s,
  time_utc, code (git sha), env, partition (hash of the tuning records). `ledger.freeze` writes the frozen
  configuration with the ledger's SHA-256; `ledger.sensitivity` gives the sensitivity table.
- Methods whose best development recall in T or C lies outside [0.2, 0.85] are listed by `tune_cli.py flags`, for
  the design check of H2.
- Picker probability traces do not depend on the tuned parameters, so they are computed once per (weights, record) on
  the GPU and stored (`pickers.save_probabilities`); each trial only extracts picks with `find_peaks` (height,
  minimum distance), as R0 and R1 did.
- One 1-D model for every method that takes one (`velocity.Velocity1D`): GaMMA gets it as its eikonal profile, PyOcto
  as a `VelocityModel1D` table, REAL as a `ttdb.txt` computed from `ttmt.traveltime.first_arrival` (true derivatives
  for REAL's linear correction, checked in the tests).

### GaMMA worker count (handoff of the R1 pick stage)

Without `ncpu` in its configuration GaMMA 80394dd sets `ncpu = max(1, min(clusters // 4, min(32, cpu_count - 1)))`
and, with more than one worker, shuffles the DBSCAN clusters with the global numpy generator and runs them in a
process pool reseeded with 42 per worker; results then depend on the core count (R1: GaMMA-1D counts changed on 16
cores and a worker hung). The wrapper always passes `ncpu = 1` (cannot be overridden), seeds numpy with 42 before each
call, and passes a fresh configuration each time, because GaMMA overwrites `config["eikonal"]` with the initialised
solver in place. The test checks that two calls give identical detections and pick assignments.

### Card mappings

| Card wording | Program parameter |
|:--|:--|
| GaMMA "time sigma or tolerance in [0.3, 2] s" | `max_sigma11` |
| PyOcto "minimum picks {3..8}" | `n_picks`; `n_p_picks = n_s_picks = n_p_and_s_picks = max(1, n_picks // 3)` as in R0/R1 |
| PyOcto "tolerance in [0.3, 2] s" | velocity-model `tolerance` |
| PyOcto "association cutoff in [0.5, 2] s" | `pick_match_tolerance` (s). PyOcto's `association_cutoff_distance` is in km, so it cannot be meant |
| REAL "minimum P picks", "minimum total picks" | `np0`, `nps0` |
| REAL "window in [1, 5] s" | `tint`, REAL's minimum interval between two events, log-uniform in [0.1, 5] s (amendment 1) |
| PyOcto (amendment 1) | `min_interevent_time` log-uniform in [0.1, 5] s (default 3 s kept two events < 3 s apart from being reported) |
| ranges without a stated distribution | uniform (log-uniform only where the card says so) |

### Fixed minimum time between output events (ADR-0009 item 1)

After amendment 1 neither REAL nor PyOcto has a fixed one. GaMMA has no explicit minimum interval (its mixture
components are separated by the data). The waveform methods are listed here when their wrappers are written.

## Pretrained weights against the atlas data (CR-22)

`ttmt.associate.provenance` (shared with "Datos reales y control de calidad", which checks the periods marked
"verify" against the dataset metadata). "period": same region and overlapping period, so the evaluation data may be in
the training set; "region": same region, other period.

| Weights | Training data | Region | Period | ridgecrest_2019 | hinet_2012 | altotiberina_c16 | altotiberina_days |
|:--|:--|:--|:--|:--|:--|:--|:--|
| phasenet/original | NCEDC catalogue recordings (Zhu and Beroza, 2019, GJI 216, 261-273, Section 2: 779,514 recordings) | northern California | 1984-2019 | - | - | - | - |
| phasenet/190703-214543 | same model as phasenet/original, TensorFlow checkpoint of AI4EPS/PhaseNet at 62005c6 (R0, R1) | northern California | 1984-2019 | - | - | - | - |
| phasenet/instance | INSTANCE (Michelini et al., 2021, ESSD 13, 5509) | Italy | 2005-2020 | - | - | period | period |
| phasenet/stead | STEAD (Mousavi et al., 2019, IEEE Access) | global | 1984-2018 | region | period | period | period |
| phasenet/scedc | SCEDC (SeisBench benchmark dataset) | southern California | 2000-2020 | period | - | - | - |
| phasenet/ethz | ETHZ (SeisBench benchmark dataset) | Switzerland | not stated | - | - | - | - |
| eqtransformer/original | STEAD (Mousavi et al., 2020, Nat. Commun.) | global | 1984-2018 | region | period | period | period |
| eqtransformer/instance | INSTANCE | Italy | 2005-2020 | - | - | period | period |
| phaseno/epoch=19-step=1140000 | NCEDC events 1984-2019 with STEAD noise (Sun et al., 2023, GRL; arXiv 2305.03269) | northern California | 1984-2019 | - | - | - | - |
| eqnet/none | no released weights for the 2022 EQNet model (environments/DEVIATIONS.md) | - | not stated | - | - | - | - |

Declared overlaps that matter for the tuned space: PhaseNet `instance` (one of the two weights the card tunes) was
trained on Italian data from 2005-2020, which covers the Alto Tiberina 2010 noise days (ADR-0009 item 3: on those
records `instance` results are reported apart, flagged "training-period overlap"); the STEAD-trained weights are
global and may contain any of the three regions. PhaseNet `original` and PhaseNO (northern California) do not cover
Ridgecrest, Hi-net or Alto Tiberina.

## Gate G3a (reproduction before any atlas data are read)

Target 1, implemented (`repro_puente.py`): GaMMA, PyOcto and REAL on the crustal synthetic scenario of Puente et al.
(2025), parameters of their Supplementary Tables S3-S5, against the event-level recalls of Table S10 (the card says
S8-S9; in the supplement S8-S9 are the subduction results). Deviations are listed in the script's docstring: travel
times from `ttmt.traveltime` instead of NonLinLoc, the SCSN station list fetched from SCEDC (the paper gives no list),
events with fewer than 10 picks not generated, GaMMA with ncpu = 1.

Cloud dry run (24 October 2026, not the gate): the SCEDC station service is not reachable from the cloud container,
so the run used a stand-in layout, the 23 Ridgecrest stations of GaMMA's tests (`tests/comparison/stations_ridgecrest.csv`
at 80394dd), with about 36 picks per event. Recall against Table S10 (5 % relative tolerance):

| Events | False picks | PyOcto (S10) | REAL (S10) | GaMMA (S10) |
|:--|:--|:--|:--|:--|
| 100 | 30 % | 0.990 (1.00) | 0.980 (0.98) | 0.530 (1.00) |
| 100 | 100 % | 1.000 (0.98) | 0.980 (0.98) | 0.430 (0.99) |
| 100 | 300 % | 1.000 (0.99) | 0.990 (0.99) | 0.520 (0.99) |
| 500 | 30 % | 0.986 (0.97) | 0.902 (0.93) | 0.566 (0.99) |
| 500 | 100 % | 0.982 (0.98) | 0.894 (0.93) | 0.630 (1.00) |
| 500 | 300 % | 0.974 (0.98) | | |

PyOcto and REAL are within 5 % in every condition run. GaMMA is not, and the cause is the stand-in layout: Table S3
sets `dbscan_min_samples = 20`, sized for the full SCSN; with 23 stations GaMMA's DBSCAN discards many events as
noise. On the same 100-event, 30 % record, `dbscan_min_samples = 5` gives recall 1.00 (precision 0.98) and 20 gives
0.53, with everything else unchanged. The gate run uses the SCEDC station list on the GPU machine
(`repro_puente.py fetch-stations`).

`gate.py write [--puente run.json] [--ethz ethz.json] [--quakemigrate qm.json]` writes
`experiments/exploratory/ATLAS-0/gate.json`, which the atlas reads before any data: one entry per target
(`puente_pyocto`, `puente_gamma`, `ethz_phasenet`, `quakemigrate_example`), passed when every published condition
run is within 5 %. A Puente run counts only on a list written by `fetch-stations` from SCEDC (CI, HH/BH, the
scenario box plus 50 km, 2020-01-01; provenance in `<list>.source.json`, checked by SHA-256). The paper gives no
station list, so the network is fixed by its source and filter, not by a count; any other list, such as the cloud
stand-in, never passes. Order of events: `gate.py` first required 100 stations (PR #22, a number with no source); the first
device fetch (31 October 2026) returned 54 stations, and the floor was replaced by the query check (PR #23) while the
first 100-event run on that list was still running, with no result yet. The reference recalls are Table S10 of the
supplement: it is the crustal scenario's event-level table, the quantity the card names; its Tables S8-S9 are the
subduction scenario. The committed file has
no target passed: the Puente targets wait on the device run, and targets 2 and 3 are not implemented yet.

Targets 2 (PhaseNet on the SeisBench ETHZ benchmark with `pick-benchmark`) and 3 (QuakeMigrate on its published
example) need the GPU machine and, for 3, an install; they are not done yet.

## The fine-tuned picker (CR-03)

When the phase-4 picker fine-tuned on pairs exists, its probability traces are written in the same file layout under
the weights name `finetuned`; the best pick-based pipeline is re-tuned with `pk_weights` fixed to `finetuned`, with
the same 60 seeds, and joins the comparator set. The 2 x 2 factorial (original or fine-tuned picker, x best associator
or the model) is an analysis of the atlas thread on those frozen configurations.

## Running on the device

```
conda activate ttmt-gpu
pip install --no-deps -e .
set TTMT_REAL_BIN=C:\Users\juanj\ttmt\bin\REAL.exe      (built from Dal-mzhang/REAL src/REAL/REAL.c)
python experiments/exploratory/COMPARATORS/repro_puente.py fetch-stations C:\Users\juanj\ttmt\g3a\scsn.csv
python experiments/exploratory/COMPARATORS/repro_puente.py run --stations C:\Users\juanj\ttmt\g3a\scsn.csv --out C:\Users\juanj\ttmt\g3a\puente.json
python experiments/exploratory/COMPARATORS/repro_puente.py report C:\Users\juanj\ttmt\g3a\puente.json
python tools/build_semisynth.py --dataset ridgecrest_2019 --source egf --per-block 10      (generator v2, PR #20)
python tools/build_semisynth.py --dataset hinet_2012 --source egf --per-block 10
python experiments/exploratory/COMPARATORS/tune_cli.py records ss_ridgecrest_2019_egf_v1 ss_hinet_2012_egf_v1 --out C:\Users\juanj\ttmt\tune
python experiments/exploratory/COMPARATORS/tune_cli.py run phasenet+pyocto C:\Users\juanj\ttmt\tune\ss_ridgecrest_2019_egf_v1\tune_phasenet.json C:\Users\juanj\ttmt\tune\ss_hinet_2012_egf_v1\tune_phasenet.json --env ttmt-gpu
```
`records` (`ttmt.associate.records.write_manifest`) opens each set through `SetReader(set_id, "tune")` only, checks
the rebuilt waveforms by fingerprint, runs PhaseNet with the `original` and `instance` weights on every realisation
and writes one manifest per picker model in the `ttmt.associate.tune.load_manifest` format. Each realisation gives
one record per station count of the card that the network has (Ridgecrest, 11 stations: k = 4, 6, 8; Hi-net, 7
stations: k = 4, 6), named by the atlas run name `label__k<k>`, and is its own bootstrap block. Tuning uses the base
runs only (the contamination variants are reading only). The partition label is the generator's. The distance axis
uses the set manifest's `catalogue_centroid` (card item 4); a manifest without it is refused.
