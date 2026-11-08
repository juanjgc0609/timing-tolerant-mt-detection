# Experiment card: M6C-SUBSAMPLE. Gate G-M6c, sub-sample deflation in the greedy loop

Stage: S1 (done, 8 November 2026). Plan: T4.1 (step 1 of the build order of `especificacion-modelo-v2.md` §8).
Exploratory, simulation only; no registered number changes (the registered runs keep `deflation="grid"`, E62).

## Question
After the greedy loop deflates the atom it selected, is the peak residual within +-0.5 s of each strong arrival of
that source at most 5 % of the arrival (gate G-M6c), in noise-free synthetics? The grid deflation leaves a median of
16-19 % in the phase-5 simulation (Appendix AB, B.6; `THEORY-5/diag_greedy_out.txt`).

## Mechanism tested
Four ways of subtracting the selected atom, on the same scenes:
- `grid`: registered (response centred on the read sample);
- `fractional`: E62's correction (response delayed by the sub-sample delay fitted over the +-3 s response support);
- `refit`: M6.3 with the amplitude re-fitted per block (free 3-vector per block);
- `v2-onset` / `v2-joint`: the v2 loop (`ttmt.v2.loop.v2_loop`): onset reading (M1) or joint reading (M4), ACAP or
  MT-FW ranking, and deflation at the time refined on the waveform in +-0.25 s with the amplitude re-read there and
  the moment tensor refitted on it (the atom stays a moment-tensor atom).

## Data
Simulation: 8-station ring (r = 20 km), node lattice 21 x 21 at 1.5 km offset 0.37 km, z = 8 km, homogeneous
VP 6, VS 6/1.73, continuous source positions and origin times, block statics U(-0.15, 0.15) s, 2 Hz kernel, rounding
"nearest". One source, or two with origin separation 3, 1 or 0.5 s (epicentres >= 6 km apart). Synthetic picker for
the v2 loops: Gaussian bump (0.04 s) at onset = peak - 0.25 s, onset error N(0, 0.02 s), lag known. 60 replicates per
configuration, seed 20261108. 10.7 min on 4 cloud cores (a 4-replicate smoke run came first).

## Endpoint
Median, 90th percentile and maximum of the peak residual over blocks with ||g_ib|| >= 0.1 max; energy left; for two
sources, the fraction of the OTHER source's energy that the deflation removed.

## Result (`run_out.txt`)
| scene | grid | fractional | refit | v2-onset | v2-joint |
|:--|--:|--:|--:|--:|--:|
| one source, MT-FW | 0.197 | 0.043 | 0.000 | **0.035** | 0.035 |
| one source, ACAP | 0.199 | 0.043 | 0.000 | **0.036** | 0.037 |
| two, 3 s, ACAP | 0.239 | 0.084 | 0.000 | 0.062 | 0.077 |
| two, 1 s, ACAP | 0.269 | 0.123 | 0.004 | 0.110 | 0.110 |
| two, 0.5 s, ACAP | 0.274 | 0.124 | 0.006 | 0.111 | 0.120 |

Other source removed (median, ACAP, 0.5 s): grid 0.018, fractional 0.024, refit 0.176, v2-onset 0.016.

- **G-M6c passes for isolated events** with `fractional` (0.043) and with the v2 loop (0.035), against 0.197 for the
  registered grid deflation.
- **G-M6c fails with two sources** for every greedy deflation that keeps the other event intact (0.06 to 0.12). The
  only one that passes, `refit`, does it by eating 15-20 % of the other event's energy: it is not admissible.
- Diagnosis: with two sources the residual of a greedy step is set by the arrivals of the second event that fall
  under the first one's response (a shared onset at one station leaves that block unresolved, see
  `tests/test_v2.py` and the seed-3 scene). It cannot be fixed one atom at a time; it needs the joint refit of all
  accepted atoms (M6, ADCG). This agrees with the stage-loss diagnosis of the atlas (PR #38: in the T regime 139 of
  146 PyOcto failures are merged pairs).
- The joint reading (M4) does not lower the deflation residual, as expected: it chooses which maximum to read, it does
  not change how an atom is subtracted.
- Two v2 corrections came out of this run and are in the code: (i) deflating at onset + lag inherits the picker's
  error (0.02 s left a median 0.22 in the 4-replicate smoke run); the time is now refined on the waveform; (ii) re-reading the amplitude at the
  refined time and refitting the tensor brought the isolated median from 0.10 (smoke run) to 0.035. Suppressing a used onset in q
  was turned off by default: when two events share an onset maximum at a station it took that block from the second.

## Environment
Python 3.13.18, NumPy 2.5.3, SciPy 1.18.1 (cloud container; not the frozen core environment).

    python experiments/exploratory/M6C-SUBSAMPLE/run.py [NREP] > experiments/exploratory/M6C-SUBSAMPLE/run_out.txt
