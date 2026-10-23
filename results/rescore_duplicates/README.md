# Rescoring with a one-to-one hit rule (2A, 2B, 2F, 2I, 2M, R0, R1)

Post hoc, descriptive. The registered results and verdicts stand as registered; this folder reports how they move
under a stricter hit rule (error register E58). Coverage by task is in the last table (critique CR-21).

## What was wrong with the registered rule

1. **Repeat detections are neutral.** Every atom is labelled with the nearest true event within 5 km
   (`froc_exact.greedy_record`, `task2c_methods.hit_of`, `c10_lib.label`), and the FROC counts each event once
   (`newhit` in `analyze_2f.Pool` and `c10_analyze.Boot`). A second atom on an already-hit event is therefore neither
   a hit nor a false event, so a method that emits several atoms per earthquake pays nothing for them. In 2M (REAL,
   DENSE) 5 to 7 % of all atoms are such repeats (`2m_dup_stats_out.txt`).
2. **The headline rule (LOC) has no origin-time test.** Thesis Section 5.2 says LOC means "within 5 km and within 3 s".
   The code applies no time test at all. On these synthetic windows the two coincide, because the true origins and
   the search grid (`DEFAULT_TAUS`, 5.0 to 8.0 s) both lie in a 3 s span. In dense windows (3 to 5 events per 40 s)
   an atom at the right epicentre can carry another event's origin time. Of the LOC hits in 2M DENSE, only 58 to 62 %
   have |origin error| <= 0.5 s, and the median error is 0.4 s.

## Corrected rule (`code/froc_strict.py`)

Atoms are matched in record order (the loop order, along which the running-minimum score of the record never
increases) to the nearest eligible true event that is **still unmatched**. An atom
with no unmatched eligible event is **false**. This is the JAFROC convention (marks after the first on a lesion are
non-lesion localisations) and the COCO-style detection convention. The matching only looks at earlier atoms, so every
threshold is still a prefix of the record and the exact FROC and the bootstrap code are unchanged; only the labels
change. Recommended headline: one-to-one with LOC+T (epicentre <= 5 km and |origin error| <= 0.5 s, the tolerance
already registered as secondary in 2I and 2M).

Rules reported: `legacy` (registered), `legacy-loct` (registered nearest-event rule plus 0.5 s, repeats neutral),
`strict` (one-to-one, LOC), `strict-loct` (one-to-one, LOC+T).

## Reproduction checks

- 2M: `rescore/rescore_2m.py legacy` on the stored parts reproduces `results/task2M/analyze_2m_out.txt` byte for
  byte (the parts-directory path in two header lines aside).
- 2F: the raw parts are not in the repository. They were regenerated with the frozen code (`sha256sum -c
  task2f_code_hashes.txt` passes; `OPENBLAS_CORETYPE=Prescott`, the same OpenBLAS builds as the run) with
  `task2f_run.py run --seeds 4242,8888`. `rescore/rescore_2f.py legacy` then reproduces every number of
  `results/task2f/analyze_2f_out.txt` except the run times (1209 of 1249 numbers identical; the other 40 are the timing section).
- In every window and method, `froc_strict.relabel(..., one_to_one=False)` on the stored locations reproduces the
  stored labels (asserted in the rescoring scripts).
- 2I: the primary arm (`prim`, 1440 parts) was regenerated the same way (`task2i_code_hashes.txt` passes) with
  `task2i_run.py run --seeds 5151,6262 --arms prim`; `rescore/rescore_2i.py legacy` reproduces every number of
  `results/task2i/analyze_2i_primary_out.txt` except the run times (2288 of 2310 numbers identical; the other 22 are
  the timing section). The secondary arms thr, s25 and s80 were not regenerated.
- The 2F and 2I parts archived in `results/task2f/records/` and `results/task2i/records/` (`*.tar.xz` with
  `SHA256SUMS`; 5.8 and 5.7 MB) are a **regeneration** with the frozen code, not the original registered records,
  which were never stored (E49). Their equivalence to the registered runs rests on reproducing every registered
  number above, not on the bytes of the original files. `tests/test_regression.py` (d) and (e) extract them and check the
  reproduction above on every CI run.
- 2A, 2B: `rescore/rescore_2a2b.py legacy` reproduces `analyze_2a_out.txt` and `analyze_2b_out.txt` (every number).
- R0, R1: `python -m ttmt.analysis.r0 eval results/taskR0/verif_r0/analysis 2000` and `python -m ttmt.analysis.r1 eval
  results/taskR1/run/hyb_eval.json results/taskR1/run/hyb_cal.json OUT` reproduce the stored analyses byte for byte
  (the R0 folder with every method; `results/taskR0/analysis_eval.txt` itself is reproduced by `run_close` for M1 and
  M3, `tests/test_regression.py` (b)); `--match one_to_one` gives the rescoring.

## 2M (REAL arm, cells pooled, @0.25; seed 2718 primary)

| Contrast | Registered (LOC) | One-to-one LOC | Registered LOC+T | One-to-one LOC+T |
|---|---|---|---|---|
| 2M-a MT@S0 - BOMP@S0, 2718 | +0.228 [+0.148, +0.298] | +0.205 [+0.133, +0.279] | +0.085 [+0.053, +0.120] | +0.083 [+0.051, +0.118] |
| 2M-a, 3141 | +0.057 [+0.004, +0.121] | +0.065 [+0.013, +0.114] | +0.037 [-0.003, +0.076] | +0.036 [-0.003, +0.073] |
| 2M-b HUBER-N@S0 - BOMP@S0, 2718 | +0.241 [+0.179, +0.293] | +0.232 [+0.178, +0.286] | +0.104 [+0.080, +0.134] | +0.102 [+0.078, +0.133] |
| 2M-b, 3141 | +0.120 [+0.070, +0.168] | +0.128 [+0.074, +0.170] | +0.068 [+0.038, +0.108] | +0.067 [+0.037, +0.105] |
| 2M-d HUBER-N@S0 - MT@S0, 2718 (margin -0.03) | +0.014 [-0.038, +0.056], not shown | +0.027 [-0.021, +0.070], **non-inferior** | | |

2M-a and 2M-b stay met under the one-to-one LOC rule. MT@S0 recall in 2718 DENSE: 0.325 (LOC) -> 0.302 (one-to-one
LOC) -> 0.124 (one-to-one LOC+T). Full output: `2m_strict_out.txt` (DENSE and STD, LOC and LOC+T, per cell).

## 2F (MT-FW+SC - BOMP-PK+SC, dense arm, @0.25; seeds 4242 / 8888)

| Reading | legacy | strict | legacy-loct | strict-loct |
|---|---|---|---|---|
| 2F-a thin | +0.072 / +0.074, met | +0.066 / +0.073, met | +0.045 / +0.042, met | +0.045 / +0.042, met |
| 2F-a strong | +0.051 / +0.043, met | +0.043 / +0.044, met | +0.024 / +0.027, **not met** | +0.023 / +0.027, **not met** |
| 2F-b MT-FW+SC - NOCAP | +0.060 / +0.059, met | +0.057 / +0.056, met | +0.035 / +0.034, met | +0.035 / +0.034, met |
| 2F-c TOST (std arm), 90% CI | [-0.039, +0.034] / [-0.048, +0.027], met | [-0.045, +0.026] / [-0.047, +0.029], met | [-0.021, +0.057] / [-0.040, +0.053], **not met** | [-0.021, +0.056] / [-0.040, +0.053], **not met** |

Under LOC+T, 2F-a strong fails only on the point criterion (>= 0.03): the lower bound stays above 0. 2F-c fails on the
upper side (MT-FW+SC better than BOMP-PK+SC by up to +0.057), which is a failure of equivalence, not evidence of harm.
The repeat-hit correction alone moves every 2F contrast by at most 0.008; the origin-time test is what changes the
verdicts. In 2F, BOMP-PK-CAP+SC (descriptive) stores no locations: under `strict` it gets
`froc_strict.drop_duplicates` of its stored labels, and it is left out of the LOC+T runs.

## 2I (HUBER-N+SC versus BOMP-PK+SC and MT-FW+SC, primary arm, @0.25; seeds 5151 / 6262)

| Reading | legacy | strict | legacy-loct | strict-loct |
|---|---|---|---|---|
| 2I-a HUBER-N - BOMP-PK, dense | +0.058 / +0.053, met | +0.059 / +0.056, met | +0.031 / +0.019, **not met** | +0.032 / +0.019, **not met** |
| 2I-b HUBER-N - MT-FW (margin -0.03), lower bound | -0.023 / -0.018, met | -0.022 / -0.018, met | -0.015 / -0.019, met | -0.015 / -0.018, met |
| 2I-c TOST (standard windows), 90% CI | [-0.010, +0.027] / [-0.013, +0.024], met | [-0.010, +0.028] / [-0.013, +0.024], met | [+0.011, +0.049] / [-0.008, +0.030], met | [+0.012, +0.049] / [-0.007, +0.031], met |

Under LOC+T, 2I-a fails in seed 6262 on the point criterion only (+0.019 against 0.03; lower bound +0.008 > 0), so the
registered interpretation flips to "a fails". The one-to-one rule alone moves the three registered 2I contrasts by at most 0.003.
Every 2I method stores locations, so the full one-to-one rule applies.

## 2A and 2B (stored records without locations)

The parts of 2A and 2B store (score, label) per atom but no locations or origin times, so only the weaker variant
applies: `froc_strict.drop_duplicates` turns a repeated event index into a false event (no reassignment, no LOC+T).
GaMMA and PyOcto store per-window counts per configuration and keep their registered scoring. The 2A-e operating
points are re-frozen on the seed-999 calibration under the same rule (`rescore/rescore_2a2b.py dedup`).
Every registered 2A verdict (2A-a to 2A-e, both regimes) and 2B (both regimes) is unchanged. The cell differences of
the rule lines move by at most 0.06 (2A-b, 6 st., seed 2024, thin: +0.24/+0.29 -> +0.21/+0.28; one-sided 8 rises +0.31/+0.24 -> +0.37/+0.30).
The strict reading of 2A-e in the thin regime goes from not met to met (seed 7331 rate within 2x the budget in every
cell). Outputs: `2a_dedup_out.txt`, `2b_dedup_out.txt`.

## R0 and R1 (real data, stored detection lists; `ttmt.analysis.r0` / `r1 --match one_to_one`, ADR-0007)

- R0: the registered primary contrast (thesis detector minus BESTPICK, LOC+T, 30/h) stays -0.383 [-0.481, -0.310],
  and no LOC+T reading changes. Under LOC, M3 on U8 at 10/h drops from 0.067 to 0.044 (M1 - M3 at 10/h: -0.022 ->
  -0.017). Output: `r0_one_to_one_out.txt`.
- R1: the analysis is identical under both rules (R1-a 23 vs 22 of 54, p = 0.50; R1-b 23 vs 26, p = 0.97: both not
  met). Output: `r1_one_to_one_out.txt`.

## Coverage by task (critique CR-21)

| Task | One-to-one applied | Origin time (LOC+T) | Verdict change |
|---|---|---|---|
| 2A, 2B | repeats as false only (no locations stored) | not possible | none (2A-e strict reading, thin, becomes met) |
| 2C, 2D, 2E, 2L | **no**: their raw records are not stored and were not regenerated | no | open |
| 2F | yes (regenerated parts) | yes | 2F-a strong and 2F-c fail under LOC+T |
| 2I | yes (regenerated primary arm) | yes | 2I-a fails under LOC+T (seed 6262) |
| 2M | yes (stored parts) | yes (registered secondary) | 2M-d (2718) becomes non-inferior; primaries hold |
| R0 | yes | yes (registered primary) | none |
| R1 | yes | yes (registered primary) | none |

2C, 2D, 2E and 2L would need their runs regenerated from the frozen code, as was done for 2F and 2I.

## Reproduce

```bash
cd code      # run from code/: the archives extract to code/task2f_parts and code/task2i_parts, the paths the
            # frozen scripts read (rescore_2f.py and rescore_2i.py take no parts argument; they are append-only)
python rescore/rescore_2m.py strict OUT.txt 2000 2718=../results/task2M/parts_2718 3141=../results/task2M/parts_3141
python rescore/dup_stats_2m.py
tar -xJf ../results/task2f/records/task2f_parts.tar.xz       # or regenerate (about 6.4 CPU-hours):
OPENBLAS_CORETYPE=Prescott python task2f_run.py run --seeds 4242,8888 --workers 4
python rescore/rescore_2f.py strict-loct OUT.txt                                     # legacy | legacy-loct | strict
tar -xJf ../results/task2i/records/task2i_parts.tar.xz       # or regenerate (about 7.7 CPU-hours):
OPENBLAS_CORETYPE=Prescott python task2i_run.py run --seeds 5151,6262 --arms prim --workers 4
python rescore/rescore_2i.py strict-loct OUT.txt
python rescore/rescore_2a2b.py dedup WORKDIR
TTMT_DATA_ROOT=... python -m ttmt.analysis.r0 eval ../results/taskR0/verif_r0/analysis 2000 --out OUT.txt --match one_to_one
TTMT_DATA_ROOT=... python -m ttmt.analysis.r1 eval ../results/taskR1/run/hyb_eval.json ../results/taskR1/run/hyb_cal.json OUT.txt --match one_to_one
```

Environment: Python 3.11, numpy 1.26.4, scipy 1.17.1, pandas 3.0.6, scikit-learn 1.8.0, pyproj 3.7.2,
quakemigrate 1.2.2, obspy 1.5.1, threadpoolctl. The extracted parts folders are ignored by git; the archives are tracked.
