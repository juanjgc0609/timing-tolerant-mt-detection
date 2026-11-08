# Experiment card: M5-ACAP. Gates G-M5a (robustness) and G-M5b (no saturation) of the ACAP ranking

Stage: S0 (scripts written, 8 November 2026; not run yet). Plan: T4.13, step 3 of the build order of
`especificacion-modelo-v2.md` §8. Exploratory, real data; no registered number changes (the registered runs keep the
cap 6).

## Question
Can the relative cap ACAP (each block norm limited to c = 3 x the median block norm of the hypothesis,
`ttmt.detector.ACAP`) replace the registered absolute cap 6 for ranking, without losing what the cap was for?
- **G-M5a, robustness:** on the 2M transient day, ACAP's recall is at most 0.03 below the cap's, with LOC+T.
- **G-M5b, no saturation:** on the continuous Ridgecrest development segment, ACAP's AUC is >= 0.65. The AUC
  compares true hypocentres with the maxima of windows that have no event. The cap gave 0.39-0.49 (C18c).

## Design (fixed before the run)
- **G-M5a (`run_a.py`).**
  - Pipeline: the registered Task 2M one (`code/t2m_lib.py`, `code/t2m_run.py`, `results/task2M/cfg_2718.json`).
    REAL arm, noise of 2010-12-31 (the day with transients), 24 networks x 3 cells.
  - Methods: three that differ only in the block statistic of the MT loop, all with station terms S0, hard reading
    and W = 4: `MT@S0` (cap 6, identical to the registered method), `ACAP@S0` and `NOCAP@S0` (the control).
  - Endpoint: DENSE, LOC+T, recall at <= 0.25 false events per window. Hierarchical bootstrap of `t2m_analyze`,
    B = 2000.
  - Gate: the 2.5th percentile of ACAP@S0 - MT@S0 is > -0.03 (the non-inferiority reading of 2M-d).
- **G-M5b (`run_b.py`).**
  - Pipeline: the C18c oracle (`code/explore_contrib/c18c_oracle.py`) unchanged except for the statistics, on the
    Ridgecrest CAL segment.
  - Statistics: CAP6 and NOCAP must reproduce `results/explorations/c18/c18c_*.json` exactly, which checks the
    pipeline. Also ACAP3 (the gate statistic) and ACAP2 / ACAP5 (sensitivity to c, declared, outside the gate).
  - Gate: ACAP3 AUC >= 0.65 on every one of the four subnets N4, N6, N8, U8. Bootstrap 95 % interval, B = 2000:
    there are only 16 events, so the interval is wide and is reported.
- The ACAP map in `run_a.py` patches `c10_lib._vectors` with a marker value of `capv`. It is identical to
  `ttmt.detector.ACAP`'s map; this was checked on random vectors, maximum difference 0.

## Cost (Ockham)
- G-M5b: about 4 min per subnet on one CPU core (C18c: 151-261 s).
- G-M5a: three MT-loop methods (0.25 s per window each, from the 2M runtime) instead of the five registered ones,
  about 40 s per part, 72 parts.
- Small first run: G-M5b on N4; G-M5a on network 2800 (3 parts, 78 windows).
- G-M5a needs the noise of 101231. If the laptop does not have it, `code/t2m_prep.sh` step "new" builds it from
  `pympa37_days`; its sha256 is checked against `results/task2M/noise_arrays.sha256`.

## Commands
    python experiments/exploratory/M5-ACAP/run_b.py N4 > experiments/exploratory/M5-ACAP/run_b_out.txt
    python experiments/exploratory/M5-ACAP/run_a.py run NOISE_DIR work/m5a_parts 2800
    python experiments/exploratory/M5-ACAP/run_a.py analyze work/m5a_parts > experiments/exploratory/M5-ACAP/run_a_out.txt

## Result
Pending.
