# Data-use ledger and attempt log (rule 28; calendar task T-4)

Started 20 November 2026 from the experiment cards, the data cards and the error log; each line cites its source.
Every new read of a held-out or evaluation set, and every new configuration tried, adds a line here in the same PR
as the run. Reads of the reserved confirmatory windows are logged by code (`ttmt.registry`, `unseal_log.jsonl`) and
are not repeated here; none has happened.

## 1. Held-out and evaluation sets: every read so far

| Dataset, set | Read by | Purpose | Status after | Source |
|:--|:--|:--|:--|:--|
| Ridgecrest 2019 hour, R0 EVAL (origins [900, 3560] s, 45 targets) | R0 (registered run) | registered evaluation, read once | spent | `data/cards/ridgecrest_2019.md`; Thesis §6.4 |
| same | C18j (post hoc, after R1) | out-of-sample check of the C18j gate candidate: 2 and 5 of 45 | spent; declared post hoc | error log E59; Thesis §6.5 |
| same | M5-ACAP, gate G-M5b (8 Nov 2026) | held-out test of ACAP: fails, AUC 0.616-0.721 | spent ("any new c would have no held-out data left") | `M5-ACAP/card.md` |
| Ridgecrest 2019 hour, CAL (origins [0, 900) s) | R0 development; C18a-C18j; WHITEN-REVIEW (20 Nov 2026) | development | development data (any use exploratory) | data card; `WHITEN-REVIEW/card.md` |
| Hi-net 2012 swarm, evaluation (origins >= 840 s) | R1 (registered run) | registered evaluation, read once | spent | `data/cards/hinet_swarm_2012.md` |
| same | M1A-ONSET step 2 (8 Nov 2026) | confirmation of G-M1a on held-out events, read blocks only: passes (0.92, n 61) | spent for timing readers | `M1A-ONSET/card.md` |
| same, 40 s past the development end | an R1 development run | accidental read | recorded as error E56 | error log E56 |
| Hi-net 2012 swarm, development | R1 development; R1 post hoc timing; M1A-ONSET step 1 | development | development data | data card |
| Alto Tiberina 2010, day 100723 and the 2M days | C16, Task 2M, G-M5a (ACAP, both 2M seeds) | noise for injection; ACAP choices | used for M5 choices (not for S_mix) | `M5-EVALUE/card.md` (lines 68-69) |
| Alto Tiberina 2010, M5 confirmation day (gate G-M5c-2) | nobody | held-out day for G-M5c | **unused** (M5 closed before it) | `M5-EVALUE/card.md` |
| SCSN analyst picks of the Ridgecrest 2019 sequence (card: 2019-07-09 to 07-15; it calls them "days 190-203", to be reconciled by the M2 owner) | M2-GP | G-M2 (fails) | spent for M2 choices | `M2-GP/card.md` |
| Ridgecrest 2019, 19:00-20:00 as truth | R0 | — | spent as truth (CR-23): may serve as noise, not as truth for a gate | `ATLAS-0/card.md` |
| M6-CAL set A (seed 20261115) | M6-CAL step 4 | calibration of route A | used for calibration | `M6-CAL/card.md` |
| M6-CAL set B (seed 20261116) | M6-CAL step 5 (and a 4-scene smoke run) | one confirmation of route B: AUC 0.980 | **spent** | `M6-CAL/card.md` |
| tune partition b60d99b3d0bbbbfa (semisynthetic EGF sets `ss_ridgecrest_2019_egf_v1`, `ss_hinet_2012_egf_v1`; 50 records, n_T 223, n_C 305) | comparator tuning (ledgers a1-a6), STAGE-LOSS, STAGE-LOSS-B1, T420-CEILING, N1-STACK, N2-REPICK, FE-FIX | tuning and diagnosis only | tuning data: any number read on it is optimistic for every method alike; never a confirmatory number | `COMPARATORS/README.md`; each card |
| Reserved confirmatory windows (PR #50; Kīlauea 2024 replaces Antelope Valley and Bodie; Magna 2020 spare) | nobody | — | sealed, never read | `data/reserved/reservation.json` |

Open items: (1) the pilot P1 reads the tune partition again; its card must say that P1 is a tuning-partition number.
(2) No independent Hi-net set exists for station terms (risk register).

## 2. Attempts per line of work (rule 28)

Counts as the cards state them. "Not stated" means the card does not give a count; the owning thread fills it in.

| Line | Variants or configurations tried | Outcome | Source |
|:--|:--|:--|:--|
| Comparators, PyOcto | 60 trials (a1), 60 with station terms (a6) | frozen `_a6`, objective 0.2814 | `COMPARATORS/README.md`, `ST-RETUNE/README.md` |
| Comparators, GaMMA | 60 registered (37 failed by E63, 23 valid), 3 (a3, stopped), 60 (a4), 60 with terms (a6) | frozen `_a6`, objective 0.2838 | same |
| Comparators, HARPA | 60 trials (a5) | ledger `harpa-a5/tuning_ledger_a5.csv` | ATLAS-0 amendment 5 |
| Comparators, EQNet2 | stage 1 scores; stage 2 not run | pending | CHANGELOG 2026-11-17 |
| C18 (Ridgecrest CAL) | about twenty configurations on the same 16 events | gate met in sample, fails out of sample | error log E53, E59 |
| M5 e-value | 3 (product, heavier tails, station average) + `--weighted` fit + `--local` | negative | `M5-EVALUE/card.md` |
| M5 ACAP | not stated (c values on 2M seeds; one held-out read) | negative | `M5-ACAP/card.md` |
| M2-GP | not stated | G-M2 fails (gate at the noise floor) | `M2-GP/card.md` |
| M6-CAL | steps 1-5; threshold grid of 5 values on set A; routes A and B | route B adopted and confirmed | `M6-CAL/card.md` |
| N1-STACK | 2 runs (uncapped, then the v2 statistic after a decision taken after seeing data) | negative (instrument not valid for N-1) | learning log 16 Nov; PR #63 |
| N2-REPICK | 1 | low zone, parked | PR #60 |
| FE-FIX | step 1, step 2, step 3 (bands: not stated here, see `fe-fix/card.md`) | weights fix only; band unchanged | PR #68 |
| WHITEN-REVIEW | 1 (2 arms x 2 statistics, declared) | uncapped: defect does not explain the R0 failure (+0.011); registered capped: explains part (+0.079, to about chance); neither separates after the fix | this PR |
