# Exploration C19: the double-couple bias of the synthetic evidence

Seed 999 only, exploratory, not registered. Code: `code/explore_contrib/c19_lib.py`, `c19_run.py` and `c19_analyze.py`.
Full output: `c19_out.txt` in this folder.

## Question

Every synthetic task draws its sources as pure double couples (`task1w_mech.rand_dc`). Every detector (BOMP-PK,
MT-FW with cap 6, MT-FW without cap, HUBER-N) scores the full 6-D moment-tensor subspace of the homogeneous
straight-ray operator. C19 asks two things:

1. Does recall change when the sources are not double couples?
2. Does the 2F/2I result survive? That result is the bounded-influence advantage of MT-FW+SC and HUBER-N+SC over
   BOMP-PK+SC in dense windows.

## Design

Only the mechanism changes. All other parts of the window are the same across families: sources, origin times,
statics, polarisation jitter, site factors, wavelets, coda, clutter and noise. The self-test checks that the DC
family is bitwise equal to `c8_lib.generate_site`. Every family is scaled to the double couple's mean radiated
energy over the sphere, so the comparison is between radiation patterns at equal energy.

| Family | Source | Inside the point-MT subspace? |
|:--|:--|:--|
| DC | random double couple, as the generator draws it | yes |
| CLVD | randomly oriented CLVD, random sign | yes |
| ISO | explosion or implosion: P only, no S | yes |
| MIX | random symmetric matrix (generic full MT) | yes |
| FORCE | single force | no |
| DC2 | two independent DCs at the same hypocentre, 0.3 s apart, half the energy each | no (two origin times) |
| SCRAMBLE | the DC's amplitudes with a random direction per station and phase | no (diagnostic bound) |

The harness is the 2I one: `c8_lib.WNet` with iid site factors (sd 0.5), a DC catalogue of 20 events and +SC
corrections, and the `task2i_run.PRIMARY` engines.

- Networks: seed 999, j = 400..429. Cells are ladder_4, ladder_6 and unilateral_c40; regimes are thin and strong.
  That gives 180 parts.
- Windows per network and family: 5 dense windows (`draw_dense`, TAG 7101, kmax 7) and 5 standard windows
  (`draw_sources`, TAG 7102, kmax 5). In total: 900 dense windows with 3592 events and 900 standard windows with
  741 events.
- Statistics: recall at 0.1 and 0.25 false events per window. The CIs come from a hierarchical bootstrap with
  B = 1000. Every family and method uses the same resamples, so differences are paired.
- Hit rule: epicentre within 5 km.

numpy 2.5.3 and scipy 1.18.1 were used, not the frozen 1.26.4 and 1.17.1. The DC arm reproduces the 2I level:
dense recall at 0.25 is 0.186 for BOMP and 0.259 for MT-FW here, against 0.177 to 0.198 and 0.243 to 0.253 in 2I.

## Results

### Fraction of the noise-free arrival energy inside the detector's 6-D subspace at the true source

Dense windows; mean, with the 10th percentile in brackets.

| | DC | CLVD | ISO | MIX | FORCE | DC2 | SCRAMBLE |
|:--|:--|:--|:--|:--|:--|:--|:--|
| clean, gradient medium | 0.95 | 0.95 | 0.97 | 0.95 | 0.77 (0.49) | 0.95 | 0.26 |
| with jitter and site factors | 0.81 (0.70) | 0.81 | 0.81 | 0.81 | 0.67 (0.42) | 0.75* | 0.25 |

\*For DC2 this fraction ignores the 0.3 s split, which is where most of its mismatch lies.

Even the DC sources were never fully inside the subspace. About 19% of their energy falls outside it because of the
gradient medium, polarisation jitter and site factors. The "inverse crime" was therefore partial already.

### Dense windows, recall change against DC (same method, same windows), cap 0.25

| Family | BOMP | MT-FW | HUBER-N | MT-FW − BOMP | HUBER − BOMP |
|:--|:--|:--|:--|:--|:--|
| DC | — | — | — | +0.073 [+0.035, +0.126] | +0.054 [+0.022, +0.096] |
| CLVD | +0.006 | +0.006 | +0.026 | +0.073 | +0.073 |
| MIX | +0.016 | −0.005 | +0.015 | +0.053 | +0.054 |
| ISO | −0.064 [−0.100, −0.015] | −0.062 [−0.110, −0.027] | −0.062 [−0.101, −0.026] | +0.075 | +0.056 |
| FORCE | −0.012 | −0.017 | −0.014 | +0.069 | +0.052 |
| DC2 | −0.031 [−0.064, +0.013] | −0.061 [−0.096, −0.021] | −0.050 [−0.083, −0.014] | +0.043 [+0.013, +0.079] | +0.035 [+0.007, +0.064] |
| SCRAMBLE | −0.149 | −0.214 | −0.195 | +0.008 [−0.004, +0.023] | +0.008 [−0.003, +0.017] |

At cap 0.1 the pattern is the same. The exception is that the contrast also shrinks for ISO (MT-FW − BOMP goes
from +0.102 to +0.050, a change of −0.052 [−0.097, +0.008]) and for DC2 (+0.044, a change of −0.058 [−0.101, +0.007]).

### Standard windows

The bounded-influence statistics never won here, so there is no contrast to lose. DC2 lowers every method by about
0.07 to 0.11. CLVD, MIX, ISO and FORCE show no clear change (all CIs cover 0 at 0.25).

## Reading

1. **Non-DC moment tensors do not change the result.** For CLVD and generic full moment tensors, every change is
   within ±0.03 and every CI covers 0. The advantage of MT-FW and HUBER-N over BOMP-PK stays at +0.05 to +0.07. This
   is expected, because the detectors fit all 6 MT components, not a DC.
2. **Isotropic sources cost every method the same, about −0.06.** At equal energy, ISO has no S arrival. The
   detectors lose the S time and polarisation constraint, but the bounded-influence contrast at 0.25 is unchanged.
3. **A single force (about one third outside the subspace) costs little**: −0.012 to −0.017, with the contrast intact.
4. **A two-subevent rupture is the real risk.** DC2 costs MT-FW and HUBER-N about twice what it costs BOMP-PK, and
   the dense advantage shrinks from +0.073 to +0.043 (at 0.1, from +0.102 to +0.044). The change of the contrast is
   about −0.03 to −0.06, with CIs that just cover 0. The likely mechanism: the second subevent's energy is not
   deflated by the single-block fit, and the per-block bound limits how much of it either detection can use.
5. **The advantage needs moment-tensor coherence.** With SCRAMBLE all methods fall to the floor and the advantage
   disappears (+0.008). The bounded-influence gain is a property of MT-consistent radiation, not a generic robustness
   gain.

**For the thesis:** the "double couples inside the subspace" threat to validity (§7.4) can be narrowed. The
G1-R/2I advantage is robust to the source mechanism (CLVD, isotropic, generic MT, single force) at equal energy. It
is fragile to source complexity: two subevents shrink it by about 40% in dense windows.

## Limits

- This is one exploration on seed 999, not registered. The contrast changes for DC2 and ISO at 0.1 have CIs that
  just touch 0.
- Point sources, the same wavelet for force and MT sources, and equal-energy scaling. A moment-matched comparison
  would also change SNR.
- The DC2 lag (0.3 s) and energy split were fixed, not swept.
- The catalogue stays DC, so the +SC corrections are identical across families.
