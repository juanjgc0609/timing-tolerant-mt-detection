# CR-04 prevalence, Ridgecrest 2019, QTM catalogue (registered measurement, `docs/protocol/CR04_PREVALENCE.md`)

Run 2026-11-15 on the author's laptop (CPU, about 185 s). Pre-run decisions written before the catalogue was used:
catalogue columns, inventory, availability, duration and reading (below). No waveform and no reserved window were read.

## Inputs
- QTM catalogue (Ross et al. 2019): `https://service.scedc.caltech.edu/ftp/QTMcatalog-ridgecrest/ridgecrest_qtm.tar.gz`,
  SHA-256 dcb0035f69c6ea9a960e45544f31bdacb38cb7f502c3d29af827da4ff6850849. It holds one file, `ridgecrest_qtm.cat`
  (111,918 rows, 25 whitespace columns, 2019-07-04 03:14 to 2019-07-25 14:59 UTC, so 21 days and not 30). The file is
  not sorted by time. SCEDC terms: open data; the raw file is not committed, only its checksum (`data/manifests/laptop_meta_2026-11-15.tsv`).
- `cr04_prep.py` checks the ranges, takes columns 1-6 (UTC time), 8, 9 and 10 (the relocated lat, lon and depth),
  sorts by time and writes `qtm_cat.csv` (SHA-256 8b006ada4bd139557f4e37bcfbc18628ae0142a393f3f303eca718065e257ad8,
  not committed). It also writes `sta_atlas.csv`. No magnitude cut is applied.
- `sta_atlas.csv`: the 11 CI stations of the atlas (union of subnets N8 and U8 of `ttmt.datasets.RIDGECREST`), with
  coordinates from the R0 dataset's `stations.csv`. With 11 stations, k = 12 is undefined.
- Availability: epochs only (protocol 3.2; the table is flagged). Model `hk`, seed 20261101, 2,000 block-bootstrap
  replicates (1 h blocks).

Declared deviation: the pre-run check expected column 11 (magnitude) to be "roughly -1 to 7.1". The observed range is
-2.82 to 7.1 (160 rows below -1). The maxima, 6.4 and 7.1, are the Ridgecrest foreshock and mainshock. The check only
identified the column; magnitude is not used.

## Command
`python experiments/exploratory/CR04-PREVALENCE/run.py --catalogue qtm_cat.csv --stations sta_atlas.csv --model hk
--name ridgecrest_2019_qtm --out prevalence_ridgecrest_2019_qtm.json`, from a `git archive` of b4513390 with the
`traveltime` fix of this PR applied, and `PYTHONPATH=<export>/src`. Without `-I` because under `-I` the package
resolved to the laptop's editable install.

## Result (`prevalence_ridgecrest_2019_qtm.json`, SHA-256 99429c3a6d94bc9492e03dd2b5fbb81027258771687f5d1e0434c7146bab031a)
Y = share of the catalogued events whose separation is in T = [0.1, 1) s; 95 % block-bootstrap interval:

| k | n | Y | 95 % interval | Poisson reference |
|:--|--:|--:|:--|--:|
| 4 | 111,918 | 1.50 % | 1.40-1.59 % | 18.7 % |
| 6 | 111,918 | 1.67 % | 1.56-1.78 % | 16.9 % |
| 8 | 111,918 | 1.42 % | 1.32-1.52 % | 17.2 % |
| 12 | 0 | undefined | | |

Catalogue dead time (minimum, 1st percentile of inter-origin times): 0 s and 5.17 s.

Reading. Y is a lower bound (protocol 5). The Poisson reference exceeds Y by more than a factor of 2 (10 to 12), so
the claim states both values. Hypothesis, not tested: the template-matching catalogue suppresses events within a few
seconds of each other, so the true share lies between about 1.5 % and 17 %. The Poisson reference ignores spatial and
temporal clustering, so its upper value is not firm either. CR-04's usefulness bar of about 1 % is met by the lower bound.

## First run (`run1_before_fix/`, kept per rule 28)
The first run on unmodified b4513390 counted 111,893 events (Y = 1.494 %, 1.667 % and 1.417 %). The 217 QTM sources
shallower than 0.05 km (depth rounds to 0.0 km in `SourceTimes`; 196 are at or above 0 km) had no direct wave in
`traveltime.first_arrival`, so they got infinite times within about 20 km. 25 of them became NaN and were dropped;
the rest had inflated separations. The fix changes Y by at most 0.004 points. The first run's JSON hash covers its LF
text; on Windows `run.py` wrote CRLF, which this PR also fixes (`newline="\n"`).
