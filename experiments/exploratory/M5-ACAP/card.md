# Experiment card: M5-ACAP. Gates G-M5a (robustness) and G-M5b (no saturation) of the ACAP ranking

Stage: S2 (full run done, 8 November 2026: G-M5b passes, G-M5a fails; the variant ACAPF is being tested). Plan: T4.13, step 3 of the build order of
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

## Commands (current scripts)
    python experiments/exploratory/M5-ACAP/run_b.py N4 N6 N8 U8 [--eval] > experiments/exploratory/M5-ACAP/run_b_out.txt
    python experiments/exploratory/M5-ACAP/run_a.py run 2718 NOISE_101231 PARTS_2718 all
    python experiments/exploratory/M5-ACAP/run_a.py run 3141 NOISE_100210,NOISE_100211 PARTS_3141 all
    python experiments/exploratory/M5-ACAP/run_a.py analyze PARTS_2718 2000 > experiments/exploratory/M5-ACAP/run_a_out.txt

## Small run (8 November 2026, laptop, CPU, numpy 1.26.4)
- **G-M5b on N4** (87 s): CAP6 0.426 and NOCAP 0.740, both identical to c18c.
  - ACAP3 0.663 [0.500, 0.807], so it passes on N4, but the lower end of the interval reaches 0.5.
  - ACAP2 0.575, ACAP5 0.726.
- **G-M5a on network 2800** (3 parts, 48 DENSE windows, 43 s).
  - Pipeline check: MT@S0 0.042 / 0.182 (LOC+T @0.1 / @0.25) is identical to the registered parts of network 2800.
  - The noise sha256 differs from the registered one (a77cd786... against 6d90338b...). The preparation log is the
    same line for line (1826 eligible slots, 5 events, the same level quantiles), and the detector output is
    identical, so the difference is in the last bits.
  - ACAP@S0 - MT@S0 = -0.062 [-0.159, +0.038]; NOCAP@S0 - MT@S0 = -0.109 [-0.196, +0.000].
    ACAP recovers half of what removing the cap loses, but not all of it, at the point estimate.
  - With 3 networks the interval cannot meet the gate either way, so the small run decides nothing. It does point
    to a risk: ladder_6 gives -0.22, against +0.09 for unilateral_c40.
- Decision (Ockham): the full run costs about 20 min of CPU, which is cheaper than a diagnosis on 3 networks.
  ACAP5@S0 is added as a declared sensitivity. It was chosen because c = 5 was the best ACAP on Ridgecrest, and it
  stays outside the gate, so the gate keeps c = 3 as registered in the specification.
- `run_a.py analyze` now checks MT@S0 against the registered parts (`results/task2M/parts_2718`) part by part.

## Result of the full run (commit ab04b1c; `run_b_out.txt`, `run_a_out.txt`)
Laptop, CPU, numpy 1.26.4. G-M5a: 72 parts in 380 s on 4 processes. G-M5b: 541 s for the four subnets.
- **Pipeline:** MT@S0 is identical to the registered part in 72 of 72 parts, and CAP6 and NOCAP equal c18c on all four
  subnets.
- **G-M5b passes.** ACAP3 AUC: N4 0.663, N6 0.666, N8 0.700, U8 0.687, against 0.39-0.49 for the cap 6.
  - The lower 95 % ends are 0.500-0.546. With 16 events, the pass rests on the point estimate.
- **G-M5a fails.**

  | DENSE LOC+T, @0.25 (MT@S0 0.126) | difference to MT@S0 | 95 % interval |
  |:--|--:|:--|
  | ACAP@S0 (c = 3) | -0.020 | [-0.041, +0.009] |
  | ACAP5@S0 (c = 5) | -0.049 | [-0.070, -0.021] |
  | NOCAP@S0 | -0.084 | [-0.120, -0.052] |

  - The lower bound -0.041 is below the -0.03 margin.
  - The loss is in ladder_4: -0.046 [-0.079, -0.007]. unilateral_c40 gives +0.020; STD LOC+T gives +0.013.
- **Diagnosis.**
  - In a quiet record, after the normalisation by the noise sd, the median block norm is about 1.5: chi with 3
    degrees of freedom, since few blocks hold signal. ACAP3's cap is then about 4.6, which is tighter than the
    cap 6 and removes real evidence. This is worst with 8 blocks (ladder_4).
  - A looser cap lets the transients in: ACAP5 loses more, and NOCAP the most.
  - So on this day the absolute cap near 6 is close to the best bound. The saturation that ACAP corrects happens in
    loud records only (Ridgecrest), where the median grows.
- **Variant ACAPF**, cap = max(6, 3 x median) (`ttmt.detector.ACAP(3, floor=6)`).
  - In a quiet record it is the cap 6, exactly when 3 x median < 6. In a loud one it is ACAP3, so in principle it
    keeps both properties.
  - It was declared after seeing G-M5a, so it is post hoc. It is tested as follows, fixed before its run:
    - G-M5a on seed 2718: development. It is close to MT@S0 by construction, so this tests that it does not break.
    - G-M5a on seed 3141 (February days): replicate. These days have little transient noise; this checks that
      nothing is lost.
    - G-M5b on Ridgecrest CAL: development.
    - G-M5b on Ridgecrest EVAL (origins 900-3560 s), which was not used to choose the variant: held-out
      confirmation, same >= 0.65 threshold on every subnet.
  - ACAPF replaces ACAP3 as the M5 ranking statistic only if it passes all four. Otherwise the cap 6 stays, and M5
    keeps only the e-value for deciding.
- Ablation of note: NOCAP has the best Ridgecrest AUC (0.69-0.74), but loses 0.084 on the transient day. The two
  gates pull in opposite directions; only a bound that adapts to the record level satisfies both.
