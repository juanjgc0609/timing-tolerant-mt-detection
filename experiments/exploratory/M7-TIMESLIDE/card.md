# Experiment card: M7-TIMESLIDE. Step 0 of M7: does the time-slide null hold on real noise?

Stage: S1 (design fixed, 9 November 2026). Specification `especificacion-modelo-v2.md` M7.1; theory `teoria-fase5.md`
D.3. Comes from the M5-EVALUE pivot (PR #45): the decision rests on M7. Exploratory, real noise; no registered
number changes.

## Question
Theorem D.3 makes the time-slide p-value exact under "no coherent source": each station's noise is circularly
shift-invariant within the window and the stations are independent. Day 101231 (seed 2718) has transients whose
tails coincide across stations (M5-EVALUE step 2: 1.7x the independent rate). Independent per-station slides scramble
that coherence. **On the noise windows of 101231, is P(p <= a) <= a for the time-slide p-value of the registered
detector's maximum?** If not, how large is the excess, and is it the transient day only (3141 quiet days as contrast)?

## Design (fixed before the run)
- Networks: the Task 2M REAL-arm networks of `results/task2M/cfg_SEED.json` (`code/t2m_lib.py` unchanged).
  Seed 2718: day 101231, 24 networks x 3 cells. Seed 3141 (days 100210 + 100211) is the quiet contrast.
- Windows: every 40 s slot of the network's noise block except its catalogue slots (slides need no fit, so no
  FIT/TEST split). About 56 per network.
- Detector: the registered MT@S0 first atom, `c10_lib.mt_loop(E, MFn, 1, mode="hard", pr, W=4)` with the coarse
  lattice, the catalogue station corrections and pr = E.prior(sigma_catalogue). T = its score s.
- Slides: B = 99 (small run) or 199 (full run) per window. Each station's MFn (three channels together) is rolled by
  an independent uniform offset in [0, NT), seeded by (seed, 7701, cell, net, slot), independent of the data. The MFn
  scale is one global sd per window, so a roll leaves it unchanged. p = (1 + #{T_b >= T}) / (B + 1).
- Positive control: the network's catalogue windows (synthetic events in real noise): p must be small.
- Reported: P(p <= a), a in {0.01, 0.05, 0.1, 0.2}, with 95 % intervals by a bootstrap of noise blocks (the three
  cells share noise); the pooled tail (rate of real T above the pooled slid quantile q(1 - a), a in {1e-2, 1e-3},
  which is what a global threshold would use); per cell and per seed.
- Control before real data: a fake i.i.d. Student t4 day in the container (`run_ts.py fake`), where D.3 holds exactly.

## Decision rule (fixed now)
- **Null holds on 101231** if the upper 95 % bound of P(p <= a) is <= 1.5 a for a = 0.05 and 0.1, and the point
  estimate at 0.01 is <= 0.02. Then M7.1 is built on independent slides as specified.
- **Null broken** otherwise. Next step (declared, development data only): split the excess windows by the detector
  atom's coherence type, (i) consistent P and S moveout from a lattice node (an uncatalogued earthquake: a true
  coherent source, not a failure of the null) or (ii) near-simultaneous energy at all stations / moveout without an
  S-P pattern (a coherent transient). Only (ii) requires a change to the null (for example slides that keep a common
  transient track, or a transient veto before the test). The confirmation noise day stays unused at this step.

## Cost (Ockham)
CPU only: one detector call is about 0.02-0.03 s (container, 4-8 stations). Small run: networks 2800-2801, 3 cells,
B = 99, about 5 min on the laptop. Full run on 2718: 72 parts x 56 windows x 200 calls, about 45 min on 8 processes.

## Commands
    python experiments/exploratory/M7-TIMESLIDE/run_ts.py fake FAKE_DIR 10 20.3
    python experiments/exploratory/M7-TIMESLIDE/run_ts.py run 2718 FAKE_DIR TS_FAKE 99 2800,2801,2802,2803 --procs 4
    python experiments/exploratory/M7-TIMESLIDE/run_ts.py run 2718 NOISE_101231 TS_2718_SMALL 99 2800,2801 --procs 6
    python experiments/exploratory/M7-TIMESLIDE/run_ts.py run 3141 NOISE_100210,NOISE_100211 TS_3141_SMALL 99 2900,2910 --procs 6
    python experiments/exploratory/M7-TIMESLIDE/run_ts.py analyze TS_DIR [nboot]

## Results
**Container smoke (fake i.i.d. t4 day, not a result: a check of the code).** Networks 2800-2803, 3 cells, B = 99,
681 noise windows, 12 min on 4 processes (39 CPU-min; about 3.3 CPU-min per part, more than the first estimate).
- P(p <= a): 0.0073 (0.01), 0.043 [0.034, 0.050] (0.05), 0.103 [0.095, 0.114] (0.1), 0.203 (0.2). Pooled tail:
  0.0073 at 1e-2, 0.0015 at 1e-3. The null holds where D.3's premise holds, as it must.
- Positive control: 0.975 of the 240 catalogue windows have p <= 0.01 (the floor with B = 99).
- Real-data cost: small run (2 networks per seed) about 40 CPU-min; full 2718 run at B = 199 about 8 CPU-hours
  (about 1 h on 8 processes).

**Day 101231 and quiet contrast:** pending (laptop).
