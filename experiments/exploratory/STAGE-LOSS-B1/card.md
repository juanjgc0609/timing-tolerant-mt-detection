# Experiment card: STAGE-LOSS-B1. Is the pair merge of the comparators structural? And the duplicator control

Stage: exploratory (supervision review of 14 November 2026, items X6, X7, X8). No registered number changes; no
ledger touched; no configuration search (each condition is one run of a frozen configuration with declared changes).

Data: the tuning partition b60d99b3d0bbbbfa only (Ridgecrest and Hi-net EGF sets, 50 records, n_T = 223, n_C = 305),
the truth and PhaseNet `original` probability files the registered tuning read, staged read-only from the owner's
machine. Objective, budget (30 false detections per hour) and hit rule (LOC+T-V) are those of
`ttmt.associate.tune.objective`, through `STAGE-LOSS/stage_loss.py`. Environment: cloud, PyOcto 0.2.0, pyrocko;
GaMMA at commit 80394dd on PYTHONPATH. HARPA and EQNet2 were not run (amendment 5 of ATLAS-0 is not anchored).

    python experiments/exploratory/STAGE-LOSS-B1/duplicator.py <tune>      > duplicator_out.txt
    python experiments/exploratory/STAGE-LOSS-B1/oracle_pyocto.py <tune>   > oracle_pyocto_out.txt
    python experiments/exploratory/STAGE-LOSS-B1/static_pyocto.py <tune>   > static_pyocto_out.txt
    python experiments/exploratory/STAGE-LOSS-B1/oracle_gamma.py <tune>    > oracle_gamma_out.txt

`oracle_gamma_out.txt` is the output of the same script before its paths were made relative (run of 14 November 2026
in `/mnt/project-files/supervision/`); the other three were regenerated from the scripts as committed and are
identical to the first run.

## X6. Duplicator control of ADR-0008 item 3, at the registered operating point (measured)

Every detection of frozen PhaseNet+PyOcto is replaced by two copies at t -+ h with the same score.

| | T | C |
|:--|--:|--:|
| frozen comparator | 0.206 (46 of 223) | 0.311 |
| duplicator, h from 0.05 to 0.5 s and U[0.05, 0.5] | 0.07 to 0.12 (16 to 26 of 223) | 0.11 |

Reading: with the false-detection budget of the registered endpoint the duplicator loses about half the recall, so
the endpoint is not gamed by splitting detections. ADR-0008 item 3 stays as written. The duplicator "passing" in
the pair-only cells of M6-CAL came from scenes with no budget and is descriptive only.
Limits: one comparator (PyOcto), one partition, the copies share the location of the original.

## X7. PyOcto with oracle picks: the merge is travel-time misfit, not structure (measured)

Conditions: OR, picks at the injected arrival times (as STAGE-LOSS); ORS, picks at injected times minus a static
per-station term (median residual against the associator's 1-D model); ORM, picks at the times of the associator's
own model (zero residual). Suffixes: mit, min_interevent_time 0.1 s; tight, all time tolerances 0.1 s.

| T recall, best condition of each row | Ridgecrest (n = 165) | Hi-net (n = 58) |
|:--|--:|--:|
| PhaseNet picks, frozen (PN) | 0.236 | 0.069 |
| OR, frozen | 0.139 | 0.328 |
| ORS-tight | 0.461 | 0.517 |
| ORM-tight | 0.667 | 0.603 |

Median |model - injected| is 0.17 s (Ridgecrest) and 0.34 s (Hi-net); with static terms 0.06 and 0.05 s.
Reading: when the picks are consistent with its travel-time model, PyOcto separates about two of three pairs.
Reading 1 of `STAGE-LOSS/README.md` ("the merge is structural") is retracted for PyOcto. What remains is the picker:
with PhaseNet picks on overlapping waveforms T is 0.24 / 0.07.
Limits: ORM is an upper bound no real picker reaches (zero residual and every arrival picked); the "tight" settings
were chosen on this same partition, so 0.46-0.67 are optimistic for PyOcto and are not a comparator number; n = 58
on Hi-net gives a 95 % interval of about +-0.13.

## X8. GaMMA a4 with the same oracle picks (measured)

T recall stays at or below 0.29 in every condition (Ridgecrest 0.02 to 0.23, Hi-net 0.05 to 0.29), with about 83
detections for 165 Ridgecrest targets: GaMMA a4 keeps merging pairs whatever the picks. For GaMMA the merge may be
structural (hypothesis: the mixture's time covariance floor; not tested here).

## What follows (design, not result)

- A fair comparison gives the pick-based comparators station terms before any superiority claim (plan T4.25, PR #47).
  PR #47 measured that with PhaseNet picks the terms do not raise PyOcto or GaMMA, which agrees with the picker being
  the binding stage.
- Our lever in T is the waveform stage (two arrivals closer than 1 s that a picker reports as one), not the
  association.
