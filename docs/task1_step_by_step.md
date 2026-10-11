# TASK 1 — Polarization coherence of real clutter

Complete step-by-step guide. When you finish you will have the number that decides whether the absolute margins of
Appendices H, J and M hold outside the simulator.

## Result — task closed

Executed on the CM network, five stations. Full detail in Appendix T.1 of the proposal.

**The number.** At ZAR, the only station that passes the control, the chance-corrected fraction of coherent clutter
is $\eta=0.054\ [-0.049,0.164]$ over 159 false picks (90 days from
2024-10-01, 1-3 Hz band, 2 s centered window, P control 0.789 with n=109). It stays below 0.20:
Appendices H, J and M retain their validity where the polarization channel is usable, in this
precise sense: the measured false picks are not coherent with the preceding catalogued event. Coherence
with **uncatalogued** sources (microearthquakes below the ISC completeness) is not measured by
$\eta$ as defined, and it is registered (Task 1m).

**What the procedure below did not anticipate.** The control passes at one station out of five: CHI
0.256 (equal to chance, n=78), RUS 0.244 (n=45), OCA 0.303 (n=33, not distinguishable from chance) and SJC
0.474 (n=19, inconclusive owing to sample size). At CHI, metadata, fixed rotation, fixed reflection,
misassociation and an instrument change in mid-period were ruled out. One candidate explanation remains that
has not been tested: with near-vertical incidence the P has little horizontal component and the SV arrives
almost horizontal and radial, which would give P at chance and S aligned (CHI: S hits 0.677). The test is
the control restricted to incidences of 30° or more (Task 1l).

**What changed with respect to the procedure below**, all of it decided before the number it affects:

- The criterion is applied to $\eta=(f-f_0)/(f_{\rm ctrl}-f_0)$ and not to the raw fraction $f$. A
  clutter with no coherence at all yields in expectation the chance hit rate $f_0$, which is
  20/90 = 0.22 with uniform azimuths and 0.24-0.33 in the measured sets: the 0.20 threshold on
  $f$ lay below it and measured the geometry, not the coherence.
- The control uses only P picks (the motion of an S is not required to point at the event: SH is
  transverse).
- The 0.6 s centered window was replaced, per station, by the configuration that passes the control
  in two disjoint halves of events. The 0.6 s window is half pre-arrival noise.
- The reference event of a false pick is the most recent one within 200 km whose predicted P
  precedes it, not the closest in origin time (with one event every 28 minutes, that one was often a different one).
- An interval that crosses a single threshold is reported as "inconclusive", not as "intermediate".
- Two looks at ZAR (10 and 20 coda windows per event); the second was declared final before
  running it.

The scripts of this guide already include those changes. What follows in this document is the
procedure as it was registered, and it is kept as such.

## Why this task is decisive

The sixteen pre-studies (Appendices G to V) use **randomly polarized** clutter, which is the configuration
most favorable to the channel that turned out to be dominant. Appendix T measured that with **coherent** clutter
—a real neighboring event— the spurious rate of polarization rises from ~0.00 to 0.40. If real
clutter is mostly coherent, the channel hierarchy of Appendix M has to be re-measured and A2
probably does not survive.

The quantity to be measured is a single one: **what fraction of the real false picks is polarized
coherently with a nearby event.**

---

## Files

| File | What it does |
|:--|:--|
| `task1_verify.py` | verifies the polarization mathematics on synthetic data. **Run first** |
| `task1_fetch.py` | downloads continuous data, inventory and catalog |
| `task1_analyze.py` | picker, polarization and associated / unassociated labeling |
| `task1_decide.py` | coherence distribution and preregistered verdict |

Chained through `t1_data/`. Each stage leaves an inspectable artifact.

## Environment

```bash
conda activate sismo
pip install obspy seisbench pandas scipy numpy torch
python -c "import obspy, seisbench, torch; print('ok')"
```

---

## Step 0 — Verify before touching real data

```bash
python task1_verify.py
```

It must print **ALL PASS**. It checks that the estimator recovers a known back-azimuth, that the
rectilinearity responds, and that the error degrades gradually with noise. Verified here: error
0.00° without noise at eight azimuths, 3.0° with noise 0.20.

**Why it is mandatory.** On real data there is no known truth, so an error of the estimator
would be undetectable. If this step fails, the number of Task 1 means nothing.

---

## Step 1 — Fetch the data

```bash
python task1_fetch.py --start 2024-03-01 --days 2 \
    --stations OCA,RUS,BRR,PTB,ZAR,HEL,CHI,SJC
```

Start with **two days and eight stations**, not with ten days: you will want to iterate the pipeline
several times before scaling it up.

Download **continuous** data, not event windows. That is deliberate: the clutter is the object of
study, and event windows erase it before it is measured.

**The first thing to look at is the orientation table.** If channels `1`/`2` or azimuths
other than 0 and 90 appear, the script flags `needs_rotation` and stage 2 rotates to ZNE with the inventory. If
the horizontals are not oriented, **all back-azimuths come out rotated and the measurement is not
valid**.

Criterion to proceed: at least six stations with three components, coverage above 90%, and a
catalog with twenty events or more.

**If nothing downloads:** check the date —it cannot be in the future—, the `location` code, and the
channel. The script uses `location=None` and the priority `HH/BH/EH` precisely because `location=""`
fails on networks that use `00` or `10`.

**If the catalog comes back empty:** the script tries USGS, ISC, EMSC and INGV in that order. If none
responds, widen the period or lower `--minmag`. For Colombia the **SGC** catalog is more
complete; download it separately from its portal if performance is poor.

---

## Step 2 — Picks, polarization and labeling

```bash
python task1_analyze.py --band 1.0 3.0 --win 0.6
```

The band and the window are fixed to **match the synthetic generator**: 2 Hz center frequency
and a 0.6 s covariance window. Do not change them without declaring it, because the comparison with
Appendix T would no longer be valid.

The picker thresholds are deliberately low, 0.30. **We want clutter**, not clean detections.

The labeling marks a pick as associated if it falls in the window
$[t_0 + d/V_P - \tau,\ t_0 + d/V_S + \tau]$ of some catalog event at that station. The
rest are the clutter whose coherence is of interest.

Output: `t1_data/picks_polarization.csv`, one row per pick with rectilinearity, measured back-azimuth,
true azimuth to the nearest event and angular deviation.

---

## Step 3 — The number

```bash
python task1_decide.py
```

It computes the fraction of unassociated picks with an angular deviation below 20°, with a bootstrap
confidence interval, and prints the full histogram and the breakdown by distance.

### The control to look at before the result

The script computes the coherent fraction among the **associated** picks. Those come from known
events, so **they should point to their event**. If that fraction does not exceed 0.60, something is wrong
in the component orientation, the band or the window, and **the clutter number is not
interpretable**. The verdict automatically declares it INVALID in that case.

### Preregistered criterion

| Coherent fraction | Consequence |
|:--|:--|
| **< 0.20** | Appendices H, J and M retain their quantitative validity; the random clutter model is adequate |
| **0.20 to 0.50** | rescale the absolute margins with the Appendix T curve and re-measure the ladder with mixed clutter |
| **> 0.50** | re-measure the channel hierarchy of Appendix M; A2 probably does not survive in its current form |

The verdict uses the **confidence interval**, not the point estimate: it requires the whole interval to fall on
one side of the threshold. If it straddles the threshold, the result is indeterminate and more data are needed. That is
deliberate and must not be relaxed after seeing the number.

---

## The four pitfalls of this task

**Component orientation.** This is the critical one. Everything depends on the back-azimuth, and if the horizontals
are not in N and E the whole measurement is rotated. The script detects this and rotates, but check the
step 1 table with your own eyes.

**Instrument response.** For *relative* polarization between components it does not need to be removed,
provided that all three have the same one. Compare the response stages of the three components of
a station; if they differ, a correction is needed and the current script does not make it.

**180° ambiguity.** The sign of a covariance eigenvector is arbitrary, so the
polarization back-azimuth is defined modulo 180°. All comparisons are made folded to
$[0,90]$ and that is declared in the code. It is not a defect; it is the property of the estimator.

**The coherent fraction is not a single number.** It almost certainly depends on the regime: in
background seismicity the clutter is mostly incoherent noise, and in an aftershock sequence it is coda from neighbors,
that is, coherent. **Measure two periods**, one of background and one with an active sequence. That difference
is itself a result, and it is the one that decides in which regime Appendix T bites.

For Colombia the natural target for the second case is the **Bucaramanga Nest**, very active
intermediate-depth seismicity at ~150 km depth, surrounded by OCA, RUS, BRR, PTB and ZAR.

---

## What to report

The coherent fraction with its interval, the histogram, the breakdown by distance, the value of the
control with associated picks, and the exact versions:

```bash
python -c "import obspy, seisbench, torch, sys; print(sys.version); \
print('obspy', obspy.__version__, '| seisbench', seisbench.__version__, '| torch', torch.__version__)"
```

The picker model as well: `PhaseNet.from_pretrained("instance")`. If someone cannot
reproduce your number because of a difference in weights, the number is not valid.

---

## After Task 1

The fraction turned out low where it could be measured, and the control turned out to be station-dependent. That
opens Tasks 1c to 1n of `tasks_and_real_data_plan.md`. The synthetic ones (1e, 1j, 1j-b, 1k) have already been
executed and changed the question: in synthetic data, polarization has not been shown to add rejection
over the best continuous kinematic score (Appendix Y). The first ones with real data, in order, are
1n (real PhaseNet residuals, with the data already downloaded), 1l (orientation and incidence), 1f, 1d
and 1m.

If the fraction is low, the next step is **Task 1b**: characterize the validity regime of the
polarization channel with `task0_curves.py`, which takes one week and closes the only thing that remained
open from Appendix X.

If the fraction is high, the next step changes: **re-measure Appendix M with a clutter model
derived from real coda**, and that reopens the channel hierarchy. It is not bad news; it would be the
most important finding of the project to date, because it would mean that the
dominant discriminant is dominant only against noise and not against physics.
